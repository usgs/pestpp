/*
 * PEST++ C ABI implementation - walking skeleton.
 *
 * Everything a running tool needs lives in one PestppSession, owned by the handle. That is
 * what makes repeated create/destroy in a single process safe: no stack-scoped FileManager,
 * no process-global scenario, and the working directory is captured per session rather than
 * chdir'd into (which is MODFLOW 6's most-complained-about wart).
 *
 * Every entry point is wrapped by CAPI_TRY/CAPI_CATCH. An exception crossing a C frame is
 * undefined behavior, so this is mandatory rather than stylistic - and it is only possible
 * because the library's exit() calls are now throws.
 */
#include "pestpp_capi.h"

#include <exception>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <iostream>
#include <filesystem>

#include "Pest.h"
#include "FileManager.h"
#include "OutputFileWriter.h"
#include "PerformanceLog.h"
#include "RunManagerSerial.h"
#include "RunManagerPanther.h"
#include "EnsembleSmoother.h"
#include "DataAssimilator.h"
#include "MOEA.h"
#include "SQP.h"
#include "utilities.h"
#include "system_variables.h"

using namespace std;

const int PESTPP_NAME_LEN = 200;

namespace {

/** What the C ABI needs from a tool, independent of which tool it is.
 *
 * The four tools do not share a base class - ies and da derive from EnsembleMethod, mou and
 * sqp are their own - and their iteration shapes genuinely differ: ies solves lambdas, mou
 * runs a generation, sqp may issue several run batches per iteration while it line-searches.
 * Rather than force a common base onto the algorithms, the integration layer adapts them.
 * Anything a tool cannot do throws with a message saying so.
 */
struct ToolAdapter
{
    virtual ~ToolAdapter() {}

    virtual void initialize() = 0;
    /// First half of initialize(), returning how many runs the caller must service before
    /// initialize_finish(). 0 means none - either the tool supplies results instead of
    /// computing them, or this tool initializes atomically and has no batch to hand over.
    virtual int  initialize_prepare() = 0;
    virtual void initialize_finish() = 0;
    /// One iteration/generation. Returns PESTPP_RETRY when the step was rejected and the
    /// algorithm wants another attempt - an outcome, not a fault.
    virtual pestpp_status advance() = 0;
    virtual void finalize() = 0;
    virtual int  iteration() = 0;   // not const: EnsembleMethod::get_iter() is not const
    virtual bool should_terminate() = 0;

    /// The live ensemble for a pestpp_ensemble_id, or null if this tool has no such thing.
    virtual Ensemble* ensemble(int id) = 0;
    /// The parameter-side ensemble, used for queue/harvest and the CTL snapshot.
    virtual ParameterEnsemble* par_ensemble() = 0;
    virtual ObservationEnsemble* obs_ensemble() = 0;

    virtual const char* name() const = 0;
};

/** One PEST++ session: scenario, io, run manager and tool, with a stable lifetime. */
struct PestppSession
{
    pestpp_tool tool;
    string working_dir;
    string last_error;

    // Order matters: destruction runs bottom-up, and the tool refers to everything above it.
    unique_ptr<FileManager>       file_manager;
    unique_ptr<ofstream>          perf_stream;
    unique_ptr<ofstream>          rmr_stream;   // panther master log; unused when serial
    unique_ptr<PerformanceLog>    performance_log;
    unique_ptr<Pest>              pest_scenario;
    // da only: the per-cycle scenario. da's parameter and observation sets are cycle
    // dependent, so the tool must be built against the child, not the parent.
    unique_ptr<Pest>              child_scenario;
    unique_ptr<OutputFileWriter>  output_file_writer;
    unique_ptr<RunManagerAbstract> run_manager;
    unique_ptr<ToolAdapter>       adapter;

    bool initialized = false;
    // set between initialize_prepare() and initialize_finish(): the tool is half-initialized,
    // its ensembles drawn but the prior results not yet processed
    bool init_pending = false;
    // runs queued by pestpp_queue_ensemble(), awaiting pestpp_harvest_ensemble(). keyed by
    // realization name, which is what lets membership change while they are in flight.
    map<string,int> pending_runs;
    bool pending_runs_valid = false;
};

/** ies: begin_iteration -> solve (glm or mda) -> end_iteration. */
struct IesAdapter : public ToolAdapter
{
    IterEnsembleSmoother tool;
    Pest& scen;
    IesAdapter(Pest& p, FileManager& fm, OutputFileWriter& ofw, PerformanceLog* pl,
               RunManagerAbstract* rm)
        : tool(p, fm, ofw, pl, rm), scen(p) {}

    void initialize() override { tool.initialize(); }
    int  initialize_prepare() override { return tool.initialize_prepare(); }
    void initialize_finish() override { tool.initialize_finish(); }
    pestpp_status advance() override
    {
        bool use_mda = scen.get_pestpp_options().get_ies_use_mda();
        tool.begin_iteration();
        UpgradeStatus st = use_mda ? tool.solve_mda(false) : tool.solve_glm();
        tool.end_iteration();
        return (st == UpgradeStatus::REJECTED_RETRY) ? PESTPP_RETRY : PESTPP_OK;
    }
    void finalize() override { tool.finalize(); }
    int  iteration() override { return tool.get_iter(); }
    bool should_terminate() override { return tool.should_terminate(); }
    Ensemble* ensemble(int id) override
    {
        switch (id)
        {
        case PESTPP_PAR_EN:     return tool.get_pe_ptr();
        case PESTPP_OBS_EN:     return tool.get_oe_ptr();
        case PESTPP_NOISE_EN:   return tool.get_noise_oe_ptr();
        case PESTPP_WEIGHTS_EN: return tool.get_weights_ptr();
        default: return nullptr;
        }
    }
    ParameterEnsemble* par_ensemble() override { return tool.get_pe_ptr(); }
    ObservationEnsemble* obs_ensemble() override { return tool.get_oe_ptr(); }
    const char* name() const override { return "ies"; }
};

/** da, single cycle.
 *
 * pestpp-da.cpp drives a sequence of cycles, building a fresh child scenario and a fresh
 * DataAssimilator for each one. That cycle machinery is not exposed here: this adapter runs
 * one cycle, which is what a caller experimenting with the update itself wants. Driving the
 * whole cycle sequence through the API needs child-scenario construction and per-cycle run
 * manager re-initialization, and is deliberately left out rather than half-done.
 */
struct DaAdapter : public ToolAdapter
{
    DataAssimilator tool;
    int cycle;
    DaAdapter(Pest& p, FileManager& fm, OutputFileWriter& ofw, PerformanceLog* pl,
              RunManagerAbstract* rm, int _cycle)
        : tool(p, fm, ofw, pl, rm), cycle(_cycle) {}

    void initialize() override { tool.initialize(cycle, true, false); }
    int  initialize_prepare() override { return tool.initialize_prepare(cycle, true, false); }
    void initialize_finish() override { tool.initialize_finish(cycle); }
    pestpp_status advance() override { tool.da_update(cycle); return PESTPP_OK; }
    void finalize() override { tool.finalize(); }
    int  iteration() override { return tool.get_iter(); }
    bool should_terminate() override { return tool.should_terminate(); }
    Ensemble* ensemble(int id) override
    {
        switch (id)
        {
        case PESTPP_PAR_EN:     return tool.get_pe_ptr();
        case PESTPP_OBS_EN:     return tool.get_oe_ptr();
        case PESTPP_NOISE_EN:   return tool.get_noise_oe_ptr();
        case PESTPP_WEIGHTS_EN: return tool.get_weights_ptr();
        default: return nullptr;
        }
    }
    ParameterEnsemble* par_ensemble() override { return tool.get_pe_ptr(); }
    ObservationEnsemble* obs_ensemble() override { return tool.get_oe_ptr(); }
    const char* name() const override { return "da"; }
};

/** mou: one generation is generate -> run -> evaluate -> report over a GenerationContext. */
struct MouAdapter : public ToolAdapter
{
    MOEA tool;
    Pest& scen;
    int iter = 0;
    MouAdapter(Pest& p, FileManager& fm, OutputFileWriter& ofw, PerformanceLog* pl,
               RunManagerAbstract* rm)
        : tool(p, fm, ofw, pl, rm), scen(p) {}

    void initialize() override { tool.initialize(); }
    // mou's initialize() issues several population evaluations rather than one, so there is
    // no single batch to hand to a caller; it initializes atomically.
    int  initialize_prepare() override { tool.initialize(); return 0; }
    void initialize_finish() override {}
    pestpp_status advance() override
    {
        // one generation, phase by phase - the same four calls, on the same generator, that
        // iterate_to_solution() makes. the context is per-generation and non-copyable.
        GenerationContext ctx(&scen, tool.get_rand_gen_ptr());
        tool.generate_generation(ctx);
        tool.run_generation(ctx);
        tool.evaluate_generation(ctx);
        tool.report_generation(ctx);
        iter++;
        return PESTPP_OK;
    }
    void finalize() override { tool.finalize(); }
    int  iteration() override { return iter; }
    // mou has no convergence test - it runs the generations it was asked for
    bool should_terminate() override { return iter >= scen.get_control_info().noptmax; }
    Ensemble* ensemble(int id) override
    {
        switch (id)
        {
        case PESTPP_PAR_EN: return tool.get_dp_ptr();
        case PESTPP_OBS_EN: return tool.get_op_ptr();
        default: return nullptr;   // mou has no noise/weights ensembles
        }
    }
    ParameterEnsemble* par_ensemble() override { return tool.get_dp_ptr(); }
    ObservationEnsemble* obs_ensemble() override { return tool.get_op_ptr(); }
    const char* name() const override { return "mou"; }
};

/** sqp: one iteration is solve_new_ensemble(), which may issue several run batches. */
struct SqpAdapter : public ToolAdapter
{
    SeqQuadProgram tool;
    int iter = 0;
    SqpAdapter(Pest& p, FileManager& fm, OutputFileWriter& ofw, PerformanceLog* pl,
               RunManagerAbstract* rm)
        : tool(p, fm, ofw, pl, rm) {}

    void initialize() override { tool.initialize(); }
    // sqp's only batch inside initialize() is a single control-file-values run, not an
    // ensemble draw, so there is nothing worth handing over; it initializes atomically.
    int  initialize_prepare() override { tool.initialize(); return 0; }
    void initialize_finish() override {}
    pestpp_status advance() override
    {
        iter++;
        // false means the step was not accepted - an outcome, like ies REJECTED_RETRY
        return tool.solve_new_ensemble() ? PESTPP_OK : PESTPP_RETRY;
    }
    void finalize() override { tool.finalize(); }
    int  iteration() override { return iter; }
    bool should_terminate() override { return tool.should_terminate(); }
    Ensemble* ensemble(int id) override
    {
        switch (id)
        {
        case PESTPP_PAR_EN: return tool.get_dv_ptr();
        case PESTPP_OBS_EN: return tool.get_oe_ptr();
        default: return nullptr;   // sqp has no noise/weights ensembles
        }
    }
    ParameterEnsemble* par_ensemble() override { return tool.get_dv_ptr(); }
    ObservationEnsemble* obs_ensemble() override { return tool.get_oe_ptr(); }
    const char* name() const override { return "sqp"; }
};

string g_create_error;

PestppSession* as_session(pestpp_handle h)
{
    return static_cast<PestppSession*>(h);
}

/** Run with the session's working directory current, then restore it.
 *
 * Deliberately NOT OperSys::chdir(): that is a no-op on posix and always has been, and its
 * narrow-string handling on windows is load-bearing for the non-ascii cwd guard in CmdLine.
 * std::filesystem is the right tool here and is safe in this file specifically because the
 * previous directory is held as a path, never as a narrow string - so no encoding conversion
 * happens in either direction. */
struct ScopedWorkingDir
{
    std::filesystem::path prev;
    bool changed = false;
    explicit ScopedWorkingDir(const string& dir)
    {
        if (dir.empty())
            return;
        std::error_code ec;
        std::filesystem::path target(dir);
        prev = std::filesystem::current_path(ec);
        if (ec)
            throw runtime_error("cannot read the current working directory: " + ec.message());
        if (std::filesystem::equivalent(prev, target, ec))
            return;
        std::filesystem::current_path(target, ec);
        if (ec)
            throw runtime_error("cannot enter working directory '" + dir + "': " + ec.message());
        changed = true;
    }
    ~ScopedWorkingDir()
    {
        if (changed)
        {
            std::error_code ec;
            std::filesystem::current_path(prev, ec);   // never throw from a dtor
        }
    }
};

} // namespace

/* Wrap every entry point. Records the message on the handle so the caller can retrieve it,
   and guarantees nothing propagates past the C boundary. */
#define CAPI_BEGIN(h)                                                          \
    PestppSession* s = as_session(h);                                          \
    if (s == nullptr) return PESTPP_INVALID_HANDLE;                            \
    try {                                                                      \
        s->last_error.clear();                                                 \
        ScopedWorkingDir _swd(s->working_dir);

#define CAPI_END()                                                             \
    } catch (const std::exception& e) {                                        \
        s->last_error = e.what();                                              \
        return PESTPP_ERROR;                                                   \
    } catch (...) {                                                            \
        s->last_error = "unknown error";                                       \
        return PESTPP_ERROR;                                                   \
    }

extern "C" {

/* Shared by pestpp_create() and pestpp_create_panther(); `port` empty means serial. */
static pestpp_status create_session(int tool, const char* ctl_file, const char* working_dir,
                                    const char* port, pestpp_handle* out)
{
    g_create_error.clear();
    if ((out == nullptr) || (ctl_file == nullptr))
    {
        g_create_error = "null argument to pestpp_create";
        return PESTPP_ERROR;
    }
    *out = nullptr;
    string panther_port = (port == nullptr) ? string() : string(port);
    unique_ptr<PestppSession> s(new PestppSession());
    try
    {
        s->tool = static_cast<pestpp_tool>(tool);
        s->working_dir = (working_dir == nullptr) ? string() : string(working_dir);
        ScopedWorkingDir swd(s->working_dir);

        PestppOptions::ToolType tool_type;
        switch (s->tool)
        {
        case PESTPP_IES: tool_type = PestppOptions::ToolType::IES; break;
        case PESTPP_DA:  tool_type = PestppOptions::ToolType::DA;  break;
        case PESTPP_MOU: tool_type = PestppOptions::ToolType::MOU; break;
        case PESTPP_SQP: tool_type = PestppOptions::ToolType::SQP; break;
        default: throw runtime_error("unknown tool id");
        }

        string ctl(ctl_file);
        string base = pest_utils::get_filename_without_ext(ctl);

        s->file_manager.reset(new FileManager());
        s->file_manager->initialize_path(base, ".");
        s->file_manager->open_default_files();

        s->perf_stream.reset(new ofstream(s->file_manager->build_filename("log")));
        s->performance_log.reset(new PerformanceLog(*s->perf_stream));

        ofstream& fout_rec = s->file_manager->rec_ofstream();

        s->pest_scenario.reset(new Pest());
        s->pest_scenario->process_ctl_file(s->file_manager->open_ifile_ext("pst"),
                                           s->file_manager->build_filename("pst"), fout_rec);
        // follow pestpp-ies.cpp's order exactly: closing the pst stream and running the
        // scenario report are part of how the scenario gets fully realized
        s->file_manager->close_file("pst");
        // da assigns per-cycle metadata to parameters/observations/interface files straight
        // after parsing (pestpp-da.cpp does this before anything else). Without it every
        // cycle filter is empty and get_child_pest() yields a scenario with no tpl files.
        if (s->tool == PESTPP_DA)
            s->pest_scenario->assign_da_cycles(fout_rec);
        s->pest_scenario->check_inputs(fout_rec);
        s->pest_scenario->get_pestpp_options_ptr()->set_iter_summary_flag(false);

        s->output_file_writer.reset(new OutputFileWriter(*s->file_manager, *s->pest_scenario, false));
        s->output_file_writer->scenario_report(fout_rec, false);

        s->pest_scenario->get_pestpp_options_ptr()->apply_tool_defaults(tool_type, fout_rec);

        s->pest_scenario->check_io(fout_rec);
        const ModelExecInfo& exi = s->pest_scenario->get_model_exec_info();
        const PestppOptions& opt = s->pest_scenario->get_pestpp_options();
        if (panther_port.empty())
        {
            s->run_manager.reset(new RunManagerSerial(
                exi.comline_vec, exi.tplfile_vec, exi.inpfile_vec, exi.insfile_vec, exi.outfile_vec,
                s->file_manager->build_filename("rns"), ".",
                opt.get_max_run_fail(), opt.get_fill_tpl_zeros(),
                opt.get_additional_ins_delimiters(), opt.get_num_tpl_ins_threads(),
                opt.get_tpl_force_decimal(), opt.get_panther_echo()));
        }
        else
        {
            // the master writes its own .rmr log, same as the executables do
            s->rmr_stream.reset(new ofstream(s->file_manager->build_filename("rmr")));
            s->run_manager.reset(new RunManagerPanther(
                s->file_manager->build_filename("rns"), panther_port, *s->rmr_stream,
                opt.get_max_run_fail(), opt.get_overdue_reched_fac(), opt.get_overdue_giveup_fac(),
                opt.get_overdue_giveup_minutes(), opt.get_panther_echo(),
                s->pest_scenario->get_ctl_ordered_par_names(),
                s->pest_scenario->get_ctl_ordered_obs_names()));
        }
        s->run_manager->set_save_all_runs(
            s->pest_scenario->get_pestpp_options().get_save_all_runs());

        // Allocate run storage before the tool exists. This writes the .rns header; without
        // it the first add_run() finds an empty file and RunStorage::increment_nruns()
        // fails on a short read. The parameter/observation names are fixed for the run.
        {
            const ParamTransformSeq& base_trans_seq = s->pest_scenario->get_base_par_tran_seq();
            Parameters cur_ctl_parameters = s->pest_scenario->get_ctl_parameters();
            s->run_manager->initialize(base_trans_seq.ctl2model_cp(cur_ctl_parameters),
                                       s->pest_scenario->get_ctl_observations());
        }

        Pest& scen = *s->pest_scenario;
        FileManager& fm = *s->file_manager;
        OutputFileWriter& ofw = *s->output_file_writer;
        PerformanceLog* pl = s->performance_log.get();
        RunManagerAbstract* rm = s->run_manager.get();
        switch (s->tool)
        {
        case PESTPP_IES: s->adapter.reset(new IesAdapter(scen, fm, ofw, pl, rm)); break;
        case PESTPP_DA:
        {
            // da's parameter/observation sets are cycle dependent, so the tool has to be
            // built against the child scenario for the cycle, exactly as pestpp-da.cpp does.
            // Cycle 0 only - see DaAdapter for why the cycle sequence is not exposed.
            vector<int> cycles_in_tables;
            vector<int> cycles = scen.get_assim_dci_cycles(fout_rec, cycles_in_tables);
            int cycle = cycles.empty() ? 0 : *min_element(cycles.begin(), cycles.end());
            s->child_scenario.reset(new Pest(scen.get_child_pest(cycle)));
            s->child_scenario->check_inputs(fout_rec, false, true, cycle);
            s->child_scenario->check_io(fout_rec);
            s->adapter.reset(new DaAdapter(*s->child_scenario, fm, ofw, pl, rm, cycle));
            break;
        }
        case PESTPP_MOU: s->adapter.reset(new MouAdapter(scen, fm, ofw, pl, rm)); break;
        case PESTPP_SQP: s->adapter.reset(new SqpAdapter(scen, fm, ofw, pl, rm)); break;
        default: throw runtime_error("unknown tool id");
        }
    }
    catch (const std::exception& e) { g_create_error = e.what(); return PESTPP_ERROR; }
    catch (...) { g_create_error = "unknown error in pestpp_create"; return PESTPP_ERROR; }

    *out = s.release();
    return PESTPP_OK;
}

pestpp_status pestpp_create(int tool, const char* ctl_file,
                            const char* working_dir, pestpp_handle* out)
{
    return create_session(tool, ctl_file, working_dir, nullptr, out);
}

pestpp_status pestpp_create_panther(int tool, const char* ctl_file, const char* working_dir,
                                    const char* port, pestpp_handle* out)
{
    if ((port == nullptr) || (*port == '\0'))
    {
        g_create_error = "pestpp_create_panther requires a port";
        return PESTPP_ERROR;
    }
    return create_session(tool, ctl_file, working_dir, port, out);
}

pestpp_status pestpp_destroy(pestpp_handle h)
{
    PestppSession* s = as_session(h);
    if (s == nullptr)
        return PESTPP_INVALID_HANDLE;
    try { ScopedWorkingDir swd(s->working_dir); delete s; }
    catch (...) { return PESTPP_ERROR; }
    return PESTPP_OK;
}

const char* pestpp_last_error(pestpp_handle h)
{
    PestppSession* s = as_session(h);
    return (s == nullptr) ? "invalid handle" : s->last_error.c_str();
}

const char* pestpp_last_create_error(void) { return g_create_error.c_str(); }

pestpp_status pestpp_initialize(pestpp_handle h)
{
    CAPI_BEGIN(h)
        s->adapter->initialize();
        s->initialized = true;
        s->init_pending = false;
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_initialize_prepare(pestpp_handle h, int* n_runs)
{
    CAPI_BEGIN(h)
        if (s->init_pending)
            throw runtime_error("initialization is already in progress; call "
                                "pestpp_initialize_finish() before preparing again");
        int n = s->adapter->initialize_prepare();
        s->init_pending = true;
        if (n_runs != nullptr)
            *n_runs = n;
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_initialize_finish(pestpp_handle h)
{
    CAPI_BEGIN(h)
        if (!s->init_pending)
            throw runtime_error("nothing to finish; call pestpp_initialize_prepare() first");
        s->adapter->initialize_finish();
        s->init_pending = false;
        s->initialized = true;
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_solve_iteration(pestpp_handle h)
{
    CAPI_BEGIN(h)
        if (s->init_pending)
            throw runtime_error("initialization is incomplete: pestpp_initialize_finish() has "
                                "not been called since pestpp_initialize_prepare()");
        if (!s->initialized)
            throw runtime_error("pestpp_initialize() must be called before pestpp_solve_iteration()");
        // each tool's own notion of one step; RETRY is an algorithmic outcome, not a fault
        return s->adapter->advance();
    CAPI_END()
}

pestpp_status pestpp_finalize(pestpp_handle h)
{
    CAPI_BEGIN(h)
        s->adapter->finalize();
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_get_iteration(pestpp_handle h, int* iter)
{
    CAPI_BEGIN(h)
        if (iter == nullptr) throw runtime_error("null out-param");
        *iter = s->adapter->iteration();
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_should_terminate(pestpp_handle h, int* out)
{
    CAPI_BEGIN(h)
        if (out == nullptr) throw runtime_error("null out-param");
        *out = s->adapter->should_terminate() ? 1 : 0;
        return PESTPP_OK;
    CAPI_END()
}

namespace {

/** Resolve an ensemble id to the live object on the tool. */
Ensemble* pick_ensemble(PestppSession* s, int id)
{
    if ((id < PESTPP_PAR_EN) || (id > PESTPP_WEIGHTS_EN))
        throw runtime_error("unknown ensemble id");
    Ensemble* e = s->adapter->ensemble(id);
    if (e == nullptr)
        throw runtime_error(string("tool '") + s->adapter->name() +
                            "' has no ensemble with that id");
    return e;
}

/** Pack names as fixed-width space-padded blocks, MODFLOW-6 style. */
pestpp_status pack_names(const vector<string>& names, char* buf, int buf_len, int* count)
{
    if (count != nullptr)
        *count = (int)names.size();
    if (buf == nullptr)
        return PESTPP_OK;               /* query-only call */
    if (buf_len < (int)names.size() * PESTPP_NAME_LEN)
        throw runtime_error("name buffer too small; call with buf=NULL to size it first");
    for (size_t i = 0; i < names.size(); i++)
    {
        char* dst = buf + (i * PESTPP_NAME_LEN);
        for (int j = 0; j < PESTPP_NAME_LEN; j++)
            dst[j] = ' ';
        const string& n = names[i];
        size_t len = n.size() < (size_t)PESTPP_NAME_LEN ? n.size() : (size_t)PESTPP_NAME_LEN;
        for (size_t j = 0; j < len; j++)
            dst[j] = n[j];
    }
    return PESTPP_OK;
}

} // namespace

pestpp_status pestpp_get_ensemble_view(pestpp_handle h, int ensemble_id,
                                       double** data, int* nrow, int* ncol)
{
    CAPI_BEGIN(h)
        if ((data == nullptr) || (nrow == nullptr) || (ncol == nullptr))
            throw runtime_error("null out-param");
        Ensemble* ens = pick_ensemble(s, ensemble_id);
        Eigen::MatrixXd* m = ens->get_eigen_ptr_4_mod();
        *data = m->data();
        *nrow = (int)m->rows();
        *ncol = (int)m->cols();
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_get_ensemble_row_names(pestpp_handle h, int ensemble_id,
                                            char* buf, int buf_len, int* count)
{
    CAPI_BEGIN(h)
        return pack_names(pick_ensemble(s, ensemble_id)->get_real_names(), buf, buf_len, count);
    CAPI_END()
}

pestpp_status pestpp_get_ensemble_col_names(pestpp_handle h, int ensemble_id,
                                            char* buf, int buf_len, int* count)
{
    CAPI_BEGIN(h)
        return pack_names(pick_ensemble(s, ensemble_id)->get_var_names(), buf, buf_len, count);
    CAPI_END()
}

pestpp_status pestpp_get_par_transform_status(pestpp_handle h, int* tstat)
{
    CAPI_BEGIN(h)
        if (tstat == nullptr) throw runtime_error("null out-param");
        *tstat = (int)s->adapter->par_ensemble()->get_trans_status();
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_set_option(pestpp_handle h, const char* key, const char* value)
{
    CAPI_BEGIN(h)
        if ((key == nullptr) || (value == nullptr)) throw runtime_error("null argument");
        PestppOptions::ARG_STATUS st =
            s->pest_scenario->get_pestpp_options_ptr()->set_option(key, value);
        if (st == PestppOptions::ARG_STATUS::ARG_NOTFOUND)
            throw runtime_error(string("unknown option '") + key + "'");
        if (st == PestppOptions::ARG_STATUS::ARG_INVALID)
            throw runtime_error(string("invalid value for option '") + key + "': " + value);
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_get_option(pestpp_handle h, const char* key,
                                char* buf, int buf_len, int* needed)
{
    CAPI_BEGIN(h)
        if (key == nullptr) throw runtime_error("null argument");
        string v = s->pest_scenario->get_pestpp_options().get_option(key);
        if (needed != nullptr)
            *needed = (int)v.size() + 1;
        if (buf == nullptr)
            return PESTPP_OK;
        if (buf_len < (int)v.size() + 1)
            throw runtime_error("option buffer too small; call with buf=NULL to size it first");
        for (size_t i = 0; i < v.size(); i++)
            buf[i] = v[i];
        buf[v.size()] = '\0';
        return PESTPP_OK;
    CAPI_END()
}

/* ---- run management ------------------------------------------------------------------ */

namespace {

/** The panther master behind this handle, or a clear error if there isn't one. */
RunManagerPanther* panther(PestppSession* s)
{
    RunManagerPanther* p = dynamic_cast<RunManagerPanther*>(s->run_manager.get());
    if (p == nullptr)
        throw runtime_error("this call needs a PANTHER run manager - the serial manager has no "
                            "workers and cannot yield mid-batch. Create the handle with "
                            "pestpp_create_panther(), or check pestpp_supports_live_control()");
    return p;
}

int status_to_c(PantherRunStatus st)
{
    switch (st)
    {
    case PantherRunStatus::QUEUED:    return PESTPP_RUN_QUEUED;
    case PantherRunStatus::RUNNING:   return PESTPP_RUN_RUNNING;
    case PantherRunStatus::COMPLETED: return PESTPP_RUN_COMPLETED;
    case PantherRunStatus::FAILED:    return PESTPP_RUN_FAILED;
    case PantherRunStatus::TIMED_OUT: return PESTPP_RUN_TIMED_OUT;
    case PantherRunStatus::CANCELLED: return PESTPP_RUN_CANCELLED;
    }
    throw runtime_error("unknown run status");
}

/** Copy one name into a fixed-width space-padded slot. */
void pack_one_name(const string& n, char* dst)
{
    for (int j = 0; j < PESTPP_NAME_LEN; j++)
        dst[j] = ' ';
    size_t len = n.size() < (size_t)PESTPP_NAME_LEN ? n.size() : (size_t)PESTPP_NAME_LEN;
    for (size_t j = 0; j < len; j++)
        dst[j] = n[j];
}

/** NUL-terminated copy into a caller buffer, truncating rather than overrunning. */
void copy_to_buf(const string& v, char* buf, int buf_len)
{
    if ((buf == nullptr) || (buf_len <= 0))
        return;
    size_t n = v.size() < (size_t)(buf_len - 1) ? v.size() : (size_t)(buf_len - 1);
    for (size_t i = 0; i < n; i++)
        buf[i] = v[i];
    buf[n] = '\0';
}

/** Fill the int-vector out-params shared by the run-history calls. */
pestpp_status emit_ids(const vector<int>& ids, int* run_ids, int max_n, int* n_out)
{
    if (n_out != nullptr)
        *n_out = (int)ids.size();
    if (run_ids == nullptr)
        return PESTPP_OK;                 /* count-only call */
    if (max_n < (int)ids.size())
        throw runtime_error("run id buffer too small; call with run_ids=NULL to size it first");
    for (size_t i = 0; i < ids.size(); i++)
        run_ids[i] = ids[i];
    return PESTPP_OK;
}

} // namespace

pestpp_status pestpp_supports_live_control(pestpp_handle h, int* out)
{
    CAPI_BEGIN(h)
        if (out == nullptr) throw runtime_error("null out-param");
        *out = (dynamic_cast<RunManagerPanther*>(s->run_manager.get()) != nullptr) ? 1 : 0;
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_begin_batch(pestpp_handle h)
{
    CAPI_BEGIN(h)
        s->run_manager->begin_batch();
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_run_slice(pestpp_handle h, double max_seconds, int* all_done)
{
    CAPI_BEGIN(h)
        if (all_done == nullptr) throw runtime_error("null out-param");
        RunManagerAbstract::RUN_SLICE_STATUS st = s->run_manager->run_slice(max_seconds);
        *all_done = (st == RunManagerAbstract::RUN_SLICE_STATUS::ALL_DONE) ? 1 : 0;
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_end_batch(pestpp_handle h)
{
    CAPI_BEGIN(h)
        s->run_manager->end_batch();
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_get_run_time_stats(pestpp_handle h, double* avg_run_sec,
                                        int* n_completed, int* n_failed, int* n_timed_out,
                                        int* n_queued, int* n_running)
{
    CAPI_BEGIN(h)
        PantherRunTimeStats st = panther(s)->get_run_time_stats();
        if (avg_run_sec  != nullptr) *avg_run_sec  = st.global_avg_run_sec;
        if (n_completed  != nullptr) *n_completed  = st.n_completed;
        if (n_failed     != nullptr) *n_failed     = st.n_failed;
        if (n_timed_out  != nullptr) *n_timed_out  = st.n_timed_out;
        if (n_queued     != nullptr) *n_queued     = st.n_queued;
        if (n_running    != nullptr) *n_running    = st.n_running;
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_get_run_states(pestpp_handle h, const int* want_ids, int n_want,
                                    int* run_ids, int* statuses, double* elapsed_sec,
                                    int* n_failures, char* hosts, int max_n, int* n_out)
{
    CAPI_BEGIN(h)
        RunManagerPanther* p = panther(s);
        vector<PantherRunState> states;
        if ((want_ids != nullptr) && (n_want > 0))
            states = p->get_run_states(vector<int>(want_ids, want_ids + n_want));
        else
            states = p->get_run_states();

        if (n_out != nullptr)
            *n_out = (int)states.size();
        // every column optional, so a caller can size first with all-NULL then ask for only
        // the columns it cares about
        bool want_any = (run_ids != nullptr) || (statuses != nullptr) ||
                        (elapsed_sec != nullptr) || (n_failures != nullptr) || (hosts != nullptr);
        if (!want_any)
            return PESTPP_OK;
        if (max_n < (int)states.size())
            throw runtime_error("run state buffers too small; call with NULL arrays to size first");
        for (size_t i = 0; i < states.size(); i++)
        {
            if (run_ids     != nullptr) run_ids[i]     = states[i].run_id;
            if (statuses    != nullptr) statuses[i]    = status_to_c(states[i].status);
            if (elapsed_sec != nullptr) elapsed_sec[i] = states[i].elapsed_sec;
            if (n_failures  != nullptr) n_failures[i]  = states[i].n_failures;
            if (hosts       != nullptr) pack_one_name(states[i].host, hosts + (i * PESTPP_NAME_LEN));
        }
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_cancel_runs(pestpp_handle h, const int* run_ids, int n, int* n_cancelled)
{
    CAPI_BEGIN(h)
        if ((run_ids == nullptr) || (n <= 0))
            throw runtime_error("pestpp_cancel_runs needs at least one run id");
        int c = panther(s)->cancel_runs(vector<int>(run_ids, run_ids + n));
        if (n_cancelled != nullptr)
            *n_cancelled = c;
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_get_worker_count(pestpp_handle h, int* n)
{
    CAPI_BEGIN(h)
        if (n == nullptr) throw runtime_error("null out-param");
        *n = (int)panther(s)->get_worker_states().size();
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_get_worker_state(pestpp_handle h, int idx,
                                      char* host_buf, int host_buf_len,
                                      char* state_buf, int state_buf_len,
                                      int* current_run_id, double* current_elapsed_sec,
                                      double* avg_runtime_sec, int* n_failed_pings)
{
    CAPI_BEGIN(h)
        vector<PantherWorkerState> ws = panther(s)->get_worker_states();
        if ((idx < 0) || (idx >= (int)ws.size()))
            throw runtime_error("worker index out of range");
        const PantherWorkerState& w = ws[idx];
        copy_to_buf(w.hostname, host_buf, host_buf_len);
        copy_to_buf(w.state, state_buf, state_buf_len);
        if (current_run_id      != nullptr) *current_run_id      = w.current_run_id;
        if (current_elapsed_sec != nullptr) *current_elapsed_sec = w.current_elapsed_sec;
        if (avg_runtime_sec     != nullptr) *avg_runtime_sec     = w.avg_runtime_sec;
        if (n_failed_pings      != nullptr) *n_failed_pings      = w.n_failed_pings;
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_get_worker_run_history(pestpp_handle h, int idx, int which,
                                            int* run_ids, int max_n, int* n_out)
{
    CAPI_BEGIN(h)
        vector<PantherWorkerState> ws = panther(s)->get_worker_states();
        if ((idx < 0) || (idx >= (int)ws.size()))
            throw runtime_error("worker index out of range");
        const PantherWorkerState& w = ws[idx];
        switch (which)
        {
        case PESTPP_WORKER_COMPLETED: return emit_ids(w.completed_runs, run_ids, max_n, n_out);
        case PESTPP_WORKER_FAILED:    return emit_ids(w.failed_runs,    run_ids, max_n, n_out);
        case PESTPP_WORKER_TIMED_OUT: return emit_ids(w.timed_out_runs, run_ids, max_n, n_out);
        default: throw runtime_error("unknown worker history selector");
        }
    CAPI_END()
}

/* ---- queue / harvest ------------------------------------------------------------------ */

namespace {

/** Unpack fixed-width space-padded names into trimmed strings. */
vector<string> unpack_names(const char* buf, int n)
{
    vector<string> names;
    if ((buf == nullptr) || (n <= 0))
        return names;
    for (int i = 0; i < n; i++)
    {
        string v(buf + (i * PESTPP_NAME_LEN), PESTPP_NAME_LEN);
        size_t end = v.find_last_not_of(' ');
        names.push_back(end == string::npos ? string() : v.substr(0, end + 1));
    }
    return names;
}

} // namespace

pestpp_status pestpp_queue_ensemble(pestpp_handle h, int* n_queued)
{
    CAPI_BEGIN(h)
        if (s->pending_runs_valid)
            throw runtime_error("runs are already queued; harvest them before queueing again");
        s->pending_runs = queue_ensemble_util(s->performance_log.get(),
                                              s->file_manager->rec_ofstream(),
                                              *s->adapter->par_ensemble(), s->run_manager.get());
        s->pending_runs_valid = true;
        if (n_queued != nullptr)
            *n_queued = (int)s->pending_runs.size();
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_harvest_ensemble(pestpp_handle h, int* n_failed)
{
    CAPI_BEGIN(h)
        if (!s->pending_runs_valid)
            throw runtime_error("no queued runs to harvest; call pestpp_queue_ensemble first");
        vector<int> failed = harvest_ensemble_util(s->performance_log.get(),
                                                   s->file_manager->rec_ofstream(),
                                                   *s->adapter->par_ensemble(), *s->adapter->obs_ensemble(),
                                                   s->run_manager.get(), false, vector<int>(),
                                                   s->pending_runs);
        // clear even on the way out, so a failed harvest cannot leave a stale map behind
        s->pending_runs.clear();
        s->pending_runs_valid = false;
        if (n_failed != nullptr)
            *n_failed = (int)failed.size();
        return PESTPP_OK;
    CAPI_END()
}

/* ---- membership ----------------------------------------------------------------------- */

namespace {

/** Apply a keep-list to every coupled ensemble at once, pairing obs rows by position. */
void keep_across_coupled(PestppSession* s, const vector<string>& par_keep)
{
    ParameterEnsemble* pe = s->adapter->par_ensemble();
    ObservationEnsemble* oe = s->adapter->obs_ensemble();
    vector<string> pnames = pe->get_real_names(), onames = oe->get_real_names();

    map<string,int> pos;
    for (int i = 0; i < (int)pnames.size(); i++)
        pos[pnames[i]] = i;

    vector<string> keep_p, keep_o;
    for (auto& n : par_keep)
    {
        map<string,int>::iterator it = pos.find(n);
        if (it == pos.end())
            throw runtime_error("no such realization: '" + n + "'");
        keep_p.push_back(pnames[it->second]);
        if (it->second < (int)onames.size())
            keep_o.push_back(onames[it->second]);
    }
    if (keep_p.empty())
        throw runtime_error("refusing to leave the ensemble empty");

    // true -> the fixed-parameter store is culled along with the rows
    pe->keep_rows(keep_p, true);
    oe->keep_rows(keep_o);
    // noise and weights are obs-side and may legitimately be absent (ies_no_noise)
    // mou and sqp have neither; ies/da may legitimately have them empty (ies_no_noise)
    Ensemble* noise = s->adapter->ensemble(PESTPP_NOISE_EN);
    Ensemble* wts = s->adapter->ensemble(PESTPP_WEIGHTS_EN);
    if ((noise != nullptr) && (noise->shape().first > 0)) noise->keep_rows(keep_o);
    if ((wts != nullptr) && (wts->shape().first > 0))     wts->keep_rows(keep_o);
}

} // namespace

pestpp_status pestpp_keep_realizations(pestpp_handle h, const char* names, int n)
{
    CAPI_BEGIN(h)
        if ((names == nullptr) || (n <= 0))
            throw runtime_error("pestpp_keep_realizations needs at least one name");
        keep_across_coupled(s, unpack_names(names, n));
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_drop_realizations(pestpp_handle h, const char* names, int n)
{
    CAPI_BEGIN(h)
        if ((names == nullptr) || (n <= 0))
            throw runtime_error("pestpp_drop_realizations needs at least one name");
        vector<string> drop = unpack_names(names, n);
        vector<string> current = s->adapter->par_ensemble()->get_real_names();
        set<string> present(current.begin(), current.end());
        // name a bad request rather than silently dropping fewer than asked
        for (auto& nm : drop)
            if (present.find(nm) == present.end())
                throw runtime_error("no such realization: '" + nm + "'");
        // express drop as keep, so both go through one code path
        set<string> dset(drop.begin(), drop.end());
        vector<string> keep;
        for (auto& nm : current)
            if (dset.find(nm) == dset.end())
                keep.push_back(nm);
        keep_across_coupled(s, keep);
        return PESTPP_OK;
    CAPI_END()
}

/* ---- parameter snapshot ---------------------------------------------------------------- */

pestpp_status pestpp_get_par_snapshot(pestpp_handle h, double* data, int max_nrow, int max_ncol,
                                      int* nrow, int* ncol, char* row_names, char* col_names)
{
    CAPI_BEGIN(h)
        ParameterSnapshot snap = s->adapter->par_ensemble()->get_ctl_snapshot();
        int nr = (int)snap.values.rows(), nc = (int)snap.values.cols();
        if (nrow != nullptr) *nrow = nr;
        if (ncol != nullptr) *ncol = nc;
        if (row_names != nullptr)
            for (int i = 0; i < nr; i++)
                pack_one_name(snap.row_names[i], row_names + (i * PESTPP_NAME_LEN));
        if (col_names != nullptr)
            for (int j = 0; j < nc; j++)
                pack_one_name(snap.col_names[j], col_names + (j * PESTPP_NAME_LEN));
        if (data == nullptr)
            return PESTPP_OK;                  /* size-only call */
        if ((max_nrow < nr) || (max_ncol < nc))
            throw runtime_error("snapshot buffer too small; call with data=NULL to size it first");
        // column-major, matching the zero-copy view and numpy order="F"
        for (int j = 0; j < nc; j++)
            for (int i = 0; i < nr; i++)
                data[i + (j * nr)] = snap.values(i, j);
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_set_par_snapshot(pestpp_handle h, const double* data, int nrow, int ncol,
                                      const char* row_names, const char* col_names)
{
    CAPI_BEGIN(h)
        if ((data == nullptr) || (row_names == nullptr) || (col_names == nullptr))
            throw runtime_error("pestpp_set_par_snapshot needs data plus row and column names");
        if ((nrow <= 0) || (ncol <= 0))
            throw runtime_error("snapshot must have at least one row and one column");
        ParameterSnapshot snap;
        snap.row_names = unpack_names(row_names, nrow);
        snap.col_names = unpack_names(col_names, ncol);
        snap.values.resize(nrow, ncol);
        for (int j = 0; j < ncol; j++)
            for (int i = 0; i < nrow; i++)
                snap.values(i, j) = data[i + (j * nrow)];
        s->adapter->par_ensemble()->set_from_ctl_snapshot(snap);
        return PESTPP_OK;
    CAPI_END()
}

} // extern "C"

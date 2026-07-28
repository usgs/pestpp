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
#include <iostream>

#include "Pest.h"
#include "FileManager.h"
#include "OutputFileWriter.h"
#include "PerformanceLog.h"
#include "RunManagerSerial.h"
#include "EnsembleSmoother.h"
#include "utilities.h"
#include "system_variables.h"

using namespace std;

const int PESTPP_NAME_LEN = 200;

namespace {

/** One PEST++ session: scenario, io, run manager and tool, with a stable lifetime. */
struct PestppSession
{
    pestpp_tool tool;
    string working_dir;
    string last_error;

    // Order matters: destruction runs bottom-up, and the tool refers to everything above it.
    unique_ptr<FileManager>       file_manager;
    unique_ptr<ofstream>          perf_stream;
    unique_ptr<PerformanceLog>    performance_log;
    unique_ptr<Pest>              pest_scenario;
    unique_ptr<OutputFileWriter>  output_file_writer;
    unique_ptr<RunManagerAbstract> run_manager;
    unique_ptr<IterEnsembleSmoother> ies;

    bool initialized = false;
};

string g_create_error;

PestppSession* as_session(pestpp_handle h)
{
    return static_cast<PestppSession*>(h);
}

/** Run `body` with the session's working directory current, then restore. */
struct ScopedWorkingDir
{
    string prev;
    bool changed = false;
    explicit ScopedWorkingDir(const string& dir)
    {
        if (dir.empty())
            return;
        prev = OperSys::getcwd();
        if (prev != dir)
        {
            OperSys::chdir(dir.c_str());
            changed = true;
        }
    }
    ~ScopedWorkingDir()
    {
        if (changed)
        {
            try { OperSys::chdir(prev.c_str()); } catch (...) { /* never throw from a dtor */ }
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

pestpp_status pestpp_create(int tool, const char* ctl_file,
                            const char* working_dir, pestpp_handle* out)
{
    g_create_error.clear();
    if ((out == nullptr) || (ctl_file == nullptr))
    {
        g_create_error = "null argument to pestpp_create";
        return PESTPP_ERROR;
    }
    *out = nullptr;
    unique_ptr<PestppSession> s(new PestppSession());
    try
    {
        s->tool = static_cast<pestpp_tool>(tool);
        s->working_dir = (working_dir == nullptr) ? string() : string(working_dir);
        ScopedWorkingDir swd(s->working_dir);

        if (s->tool != PESTPP_IES)
            throw runtime_error("only PESTPP_IES is wired up in this skeleton");

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
        s->pest_scenario->check_inputs(fout_rec);
        s->pest_scenario->get_pestpp_options_ptr()->set_iter_summary_flag(false);

        s->output_file_writer.reset(new OutputFileWriter(*s->file_manager, *s->pest_scenario, false));
        s->output_file_writer->scenario_report(fout_rec, false);

        s->pest_scenario->get_pestpp_options_ptr()->apply_tool_defaults(
            PestppOptions::ToolType::IES, fout_rec);

        s->pest_scenario->check_io(fout_rec);
        const ModelExecInfo& exi = s->pest_scenario->get_model_exec_info();
        s->run_manager.reset(new RunManagerSerial(
            exi.comline_vec, exi.tplfile_vec, exi.inpfile_vec, exi.insfile_vec, exi.outfile_vec,
            s->file_manager->build_filename("rns"), ".",
            s->pest_scenario->get_pestpp_options().get_max_run_fail(),
            s->pest_scenario->get_pestpp_options().get_fill_tpl_zeros(),
            s->pest_scenario->get_pestpp_options().get_additional_ins_delimiters(),
            s->pest_scenario->get_pestpp_options().get_num_tpl_ins_threads(),
            s->pest_scenario->get_pestpp_options().get_tpl_force_decimal(),
            s->pest_scenario->get_pestpp_options().get_panther_echo()));
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

        s->ies.reset(new IterEnsembleSmoother(*s->pest_scenario, *s->file_manager,
                                              *s->output_file_writer, s->performance_log.get(),
                                              s->run_manager.get()));
    }
    catch (const std::exception& e) { g_create_error = e.what(); return PESTPP_ERROR; }
    catch (...) { g_create_error = "unknown error in pestpp_create"; return PESTPP_ERROR; }

    *out = s.release();
    return PESTPP_OK;
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
        s->ies->initialize();
        s->initialized = true;
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_solve_iteration(pestpp_handle h)
{
    CAPI_BEGIN(h)
        if (!s->initialized)
            throw runtime_error("pestpp_initialize() must be called before pestpp_solve_iteration()");
        bool use_mda = s->pest_scenario->get_pestpp_options().get_ies_use_mda();
        s->ies->begin_iteration();
        UpgradeStatus st = use_mda ? s->ies->solve_mda(false) : s->ies->solve_glm();
        s->ies->end_iteration();
        // REJECTED_RETRY is an algorithmic outcome, not a fault - surface it as its own code
        return (st == UpgradeStatus::REJECTED_RETRY) ? PESTPP_RETRY : PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_finalize(pestpp_handle h)
{
    CAPI_BEGIN(h)
        s->ies->finalize();
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_get_iteration(pestpp_handle h, int* iter)
{
    CAPI_BEGIN(h)
        if (iter == nullptr) throw runtime_error("null out-param");
        *iter = s->ies->get_iter();
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_should_terminate(pestpp_handle h, int* out)
{
    CAPI_BEGIN(h)
        if (out == nullptr) throw runtime_error("null out-param");
        *out = s->ies->should_terminate() ? 1 : 0;
        return PESTPP_OK;
    CAPI_END()
}

namespace {

/** Resolve an ensemble id to the live object on the tool. */
Ensemble* pick_ensemble(PestppSession* s, int id)
{
    switch (id)
    {
    case PESTPP_PAR_EN:     return s->ies->get_pe_ptr();
    case PESTPP_OBS_EN:     return s->ies->get_oe_ptr();
    case PESTPP_NOISE_EN:   return s->ies->get_noise_oe_ptr();
    case PESTPP_WEIGHTS_EN: return s->ies->get_weights_ptr();
    default: throw runtime_error("unknown ensemble id");
    }
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
        *tstat = (int)s->ies->get_pe_ptr()->get_trans_status();
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

} // extern "C"

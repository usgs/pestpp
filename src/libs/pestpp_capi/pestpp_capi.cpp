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
#include "pestpp-api.h"

#include <exception>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <iostream>
#include <cstdio>
#include <cstddef>
#if defined(_WIN32)
#  include <io.h>
#  include <fcntl.h>
#  define PESTPP_DUP    _dup
#  define PESTPP_DUP2   _dup2
#  define PESTPP_CLOSE  _close
#  define PESTPP_OPEN   _open
#  define PESTPP_O_FLAGS (_O_WRONLY | _O_CREAT | _O_APPEND | _O_TEXT)
#  define PESTPP_O_MODE  (_S_IREAD | _S_IWRITE)
#  include <sys/stat.h>
#else
#  include <unistd.h>
#  include <fcntl.h>
#  define PESTPP_DUP    dup
#  define PESTPP_DUP2   dup2
#  define PESTPP_CLOSE  close
#  define PESTPP_OPEN   open
#  define PESTPP_O_FLAGS (O_WRONLY | O_CREAT | O_APPEND)
#  define PESTPP_O_MODE  0644
#endif
#include <filesystem>

#include "config_os.h"
#include "Pest.h"
#include "FileManager.h"
#include "OutputFileWriter.h"
#include "PerformanceLog.h"
#include "RunManagerSerial.h"
#include "RunManagerPanther.h"
#include "RunManagerExternal.h"
#include "EnsembleView.h"
#include "EnsembleSmoother.h"
#include "DataAssimilator.h"
#include "MOEA.h"
#include "SQP.h"
#include "utilities.h"
#include "system_variables.h"

using namespace std;

const int PESTPP_NAME_LEN = 200;
const int PESTPP_MESSAGE_LEN = 4096;

namespace {

/** An exception that carries the status code it should become at the boundary.
 *
 * Everything in here reports failure by throwing, which is right - the alternative is
 * threading a status through every helper - but a bare runtime_error can only ever come out
 * as PESTPP_ERROR. This lets a throw site say "this one is a too-small buffer" and have that
 * survive the trip to the caller, which is what makes PESTPP_BUFFER_TOO_SMALL usable as a
 * retry signal rather than as prose.
 */
struct capi_error : public std::runtime_error
{
    pestpp_status status;
    capi_error(pestpp_status st, const string& what)
        : std::runtime_error(what), status(st) {}
};

/// Shorthand for the four that come up constantly.
[[noreturn]] void bad_arg(const string& what)   { throw capi_error(PESTPP_INVALID_ARGUMENT, what); }
[[noreturn]] void bad_state(const string& what) { throw capi_error(PESTPP_INVALID_STATE, what); }
[[noreturn]] void unsupported(const string& what) { throw capi_error(PESTPP_NOT_SUPPORTED, what); }
[[noreturn]] void too_small(const string& what)
{
    throw capi_error(PESTPP_BUFFER_TOO_SMALL, what);
}

/// Keep the promise PESTPP_MESSAGE_LEN makes. Truncation marked, so nobody hunts a missing tail.
string clamp_message(const string& msg)
{
    const size_t cap = (size_t)PESTPP_MESSAGE_LEN - 1;
    if (msg.size() <= cap)
        return msg;
    return msg.substr(0, cap - 3) + "...";
}

/** A name too long for a packed slot is an error, not something to shorten.
 *
 * Truncating would hand back a name that is not the name - and the caller would then send it
 * in again through set_par_snapshot() or set_obs_weights(), where it matches nothing, or
 * worse matches a DIFFERENT parameter that happens to share the first PESTPP_NAME_LEN
 * characters. Nothing in pest++ enforces the limit, so this is where it gets enforced. */
void check_name_fits(const string& n)
{
    if (n.size() > (size_t)PESTPP_NAME_LEN)
        bad_arg("the name '" + n.substr(0, 40) + "...' is " + to_string(n.size()) +
                " characters, longer than PESTPP_NAME_LEN (" + to_string(PESTPP_NAME_LEN) +
                "); it cannot be packed without becoming a different name");
}

/** Copy one name into a fixed-width space-padded slot. Refuses to shorten. */
void pack_one_name(const string& n, char* dst)
{
    check_name_fits(n);
    for (int j = 0; j < PESTPP_NAME_LEN; j++)
        dst[j] = ' ';
    for (size_t j = 0; j < n.size(); j++)
        dst[j] = n[j];
}

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

    /// mean/std/min/max of phi. Tools without a phi throw.
    virtual void phi_summary(int phi_type, double& mean, double& sd, double& mn, double& mx) = 0;
    /// phi per realization, paired with the realization names it is defined for.
    virtual void phi_vector(int phi_type, vector<string>& names, vector<double>& vals) = 0;
    /// squared weighted residual per (realization, observation), column-major.
    virtual void phi_residuals(int phi_type, vector<string>& rnames, vector<string>& cnames,
                               vector<double>& data) = 0;
    /// recompute the cached phi from the current ensembles and weights.
    virtual void update_phi() = 0;

    /// The da cycle runs should be tagged with. Only da has one; everything else runs
    /// untagged, which is what NULL_DA_CYCLE means.
    virtual int da_cycle() const { return NetPackage::NULL_DA_CYCLE; }

    virtual const char* name() const = 0;
};

/// Map a pestpp_phi_type onto the handler's enum.
static L2PhiHandler::phiType to_phi_type(int phi_type)
{
    switch (phi_type)
    {
    case PESTPP_PHI_MEAS:      return L2PhiHandler::phiType::MEAS;
    case PESTPP_PHI_COMPOSITE: return L2PhiHandler::phiType::COMPOSITE;
    case PESTPP_PHI_REGUL:     return L2PhiHandler::phiType::REGUL;
    case PESTPP_PHI_ACTUAL:    return L2PhiHandler::phiType::ACTUAL;
    case PESTPP_PHI_NOISE:     return L2PhiHandler::phiType::NOISE;
    default: bad_arg("unknown phi type");
    }
}

/// Shared by the two EnsembleMethod tools.
static void ensemble_phi_vector(EnsembleMethod& tool, int phi_type,
                                vector<string>& names, vector<double>& vals)
{
    names.clear();
    vals.clear();
    map<string, double> m = tool.get_phi_handler().get_phi_map(to_phi_type(phi_type));
    for (auto& kv : m)
    {
        names.push_back(kv.first);
        vals.push_back(kv.second);
    }
}

/// Shared by the two EnsembleMethod tools.
static void ensemble_phi_residuals(EnsembleMethod& tool, int phi_type,
                                   vector<string>& rnames, vector<string>& cnames,
                                   vector<double>& data)
{
    if ((phi_type != PESTPP_PHI_MEAS) && (phi_type != PESTPP_PHI_ACTUAL))
        unsupported("phi residuals are only defined for PESTPP_PHI_MEAS and "
                            "PESTPP_PHI_ACTUAL");
    map<string, map<string, double>> swr = tool.get_phi_handler().get_swr_real_map(
        *tool.get_oe_ptr(), *tool.get_weights_ptr(), to_phi_type(phi_type));

    rnames.clear();
    cnames.clear();
    data.clear();
    if (swr.empty())
        return;
    for (auto& kv : swr)
        rnames.push_back(kv.first);
    for (auto& kv : swr.begin()->second)
        cnames.push_back(kv.first);

    // column-major, matching the ensemble views
    data.assign(rnames.size() * cnames.size(), 0.0);
    for (size_t i = 0; i < rnames.size(); i++)
    {
        const map<string, double>& row = swr[rnames[i]];
        for (size_t j = 0; j < cnames.size(); j++)
        {
            map<string, double>::const_iterator it = row.find(cnames[j]);
            if (it != row.end())
                data[i + (j * rnames.size())] = it->second;
        }
    }
}

/// Shared by the two EnsembleMethod tools.
static void ensemble_phi_summary(EnsembleMethod& tool, int phi_type,
                                 double& mean, double& sd, double& mn, double& mx)
{
    L2PhiHandler::phiType pt = to_phi_type(phi_type);
    L2PhiHandler& ph = tool.get_phi_handler();
    mean = ph.get_mean(pt);
    sd   = ph.get_std(pt);
    mn   = ph.get_min(pt);
    mx   = ph.get_max(pt);
}

/** One PEST++ session: scenario, io, run manager and tool, with a stable lifetime. */
/** Stamped at the head of every live session and zeroed on destroy.
 *
 * The point is that pestpp_handle is a void*, so a caller can hand us anything at all - a
 * stale pointer, a destroyed handle, an int that was never a pointer. Without a marker the
 * only bad handle we can detect is NULL and everything else is undefined behavior, with
 * double-destroy being a double-free. Reading four bytes through a wild pointer can still
 * fault, so this is not a proof of validity; it catches the two mistakes that actually
 * happen - use-after-destroy and destroy-twice - and turns them into an error code. */
const unsigned int PESTPP_SESSION_MAGIC = 0x50455354u;   // 'PEST'

/** A handed-out zero-copy view, kept so pestpp_view_is_valid() can answer honestly.
 *
 * EnsembleView already does the hard part - it holds a weak guard token from the ensemble and
 * re-checks the buffer address and dimensions on every access, which is how it notices a
 * reallocation that no bookkeeping saw. All this adds is an integer name for one, so the same
 * question can be asked from C. */
struct ViewRecord
{
    unique_ptr<EnsembleView> view;
    int ensemble_id = -1;
};

struct PestppSession
{
    unsigned int magic = PESTPP_SESSION_MAGIC;   // must stay the first member
    pestpp_tool tool;
    int run_manager_type = 0;      // pestpp_run_manager
    string working_dir;
    // RunManagerExternal takes its storage filename by non-const reference, so it needs to
    // outlive the constructor call
    string rns_filename;
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
    // runs queued by pestpp_queue_runs(), awaiting pestpp_process_runs(). keyed by
    // realization name, which is what lets membership change while they are in flight.
    map<string,int> pending_runs;
    bool pending_runs_valid = false;

    // outstanding zero-copy views, by token. tokens never repeat within a session, so a
    // released token answers "not valid" rather than aliasing a later view.
    map<int, ViewRecord> views;
    int next_view_token = 1;
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
    void phi_summary(int t, double& mean, double& sd, double& mn, double& mx) override
    { ensemble_phi_summary(tool, t, mean, sd, mn, mx); }
    void phi_vector(int t, vector<string>& names, vector<double>& vals) override
    { ensemble_phi_vector(tool, t, names, vals); }
    void phi_residuals(int t, vector<string>& rn, vector<string>& cn, vector<double>& d) override
    { ensemble_phi_residuals(tool, t, rn, cn, d); }
    void update_phi() override
    { tool.get_phi_handler().update(*tool.get_oe_ptr(), *tool.get_pe_ptr(), *tool.get_weights_ptr()); }
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
    int  da_cycle() const override { return cycle; }
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
    void phi_summary(int t, double& mean, double& sd, double& mn, double& mx) override
    { ensemble_phi_summary(tool, t, mean, sd, mn, mx); }
    void phi_vector(int t, vector<string>& names, vector<double>& vals) override
    { ensemble_phi_vector(tool, t, names, vals); }
    void phi_residuals(int t, vector<string>& rn, vector<string>& cn, vector<double>& d) override
    { ensemble_phi_residuals(tool, t, rn, cn, d); }
    void update_phi() override
    { tool.get_phi_handler().update(*tool.get_oe_ptr(), *tool.get_pe_ptr(), *tool.get_weights_ptr()); }
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
        // ONE call into the tool, not a re-implementation of its loop body. Doing the four
        // phase calls here instead left MOEA's own generation counter untouched, so every
        // generation reported itself as the same number and wrote over the same .N. files.
        tool.solve_generation();
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
    void phi_summary(int, double&, double&, double&, double&) override
    { unsupported("mou optimizes objectives rather than minimizing a phi"); }
    void phi_vector(int, vector<string>&, vector<double>&) override
    { unsupported("mou optimizes objectives rather than minimizing a phi"); }
    void phi_residuals(int, vector<string>&, vector<string>&, vector<double>&) override
    { unsupported("mou optimizes objectives rather than minimizing a phi"); }
    void update_phi() override { unsupported("mou optimizes objectives rather than minimizing a phi"); }
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
        // ONE call into the tool. solve_new_ensemble() on its own is not an sqp iteration -
        // it is the first statement of one, and skips advancing the tool's own `iter` (which
        // names every .N. output file), the gradient runs, the CMA and hessian updates, the
        // constraint report and the pcs summary.
        bool accept = false;
        tool.solve_iteration(accept);
        // not accepted means the step was rejected - an outcome, like ies REJECTED_RETRY
        return accept ? PESTPP_OK : PESTPP_RETRY;
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
    void phi_summary(int, double&, double&, double&, double&) override
    { unsupported("sqp has an objective function, not a phi over realizations"); }
    void phi_vector(int, vector<string>&, vector<double>&) override
    { unsupported("sqp has an objective function, not a phi over realizations"); }
    void phi_residuals(int, vector<string>&, vector<string>&, vector<double>&) override
    { unsupported("sqp has an objective function, not a phi over realizations"); }
    void update_phi() override { unsupported("sqp has an objective function, not a phi over realizations"); }
    const char* name() const override { return "sqp"; }
};

/// The message for the handle-less entry points, which have nowhere else to put one.
string g_create_error;

/** Outstanding output redirects, innermost last.
 *
 * There is exactly one file descriptor 1 in the process, so redirecting it is not a per-handle
 * operation however much it looks like one. Nesting works - each redirect saves what stdout
 * was, so unwinding in the reverse order puts every layer back - but UNWINDING OUT OF ORDER
 * does not, and it fails silently and permanently:
 *
 *     A redirects -> saved_a = the console
 *     B redirects -> saved_b = A's log file
 *     A restores  -> stdout = the console        (B is still "redirected", to nothing)
 *     B restores  -> stdout = A's log file       <- and stays there forever
 *
 * The old signature handed the caller the raw saved descriptor, so there was nothing to check
 * against: any int was as plausible as any other. Now the caller gets an opaque TOKEN, the
 * descriptor stays in here, and restoring anything but the innermost redirect is refused.
 * Enforcing the order beats guessing at a repair - silently unwinding somebody else's redirect
 * is the corruption, not the fix. */
struct RedirectRecord
{
    int token;
    int saved_fd;
};
vector<RedirectRecord> g_redirect_stack;
int g_next_redirect_token = 1;

/** A session pointer, or null if this handle is not one we handed out and still own. */
PestppSession* as_session(pestpp_handle h)
{
    if (h == nullptr)
        return nullptr;
    PestppSession* s = static_cast<PestppSession*>(h);
    if (s->magic != PESTPP_SESSION_MAGIC)
        return nullptr;
    return s;
}

/// Set when a working-directory restore failed; makes the next entry point refuse.
string g_cwd_restore_error;

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
        if (!changed)
            return;
        std::error_code ec;
        std::filesystem::current_path(prev, ec);
        if (ec)
        {
            // one retry: on windows a transient holder (indexer, antivirus, a virus scanner
            // walking files the model just wrote) is the usual reason this fails
            ec.clear();
            std::filesystem::current_path(prev, ec);
        }
        if (ec)
        {
            // A destructor cannot throw, but neither can this be swallowed: the process is now
            // sitting in a directory nobody expects, and every relative path afterwards - in
            // this library and in the host program - silently resolves somewhere else. Record
            // it so the next entry point refuses rather than continuing in an unknown state.
            try
            {
                g_cwd_restore_error =
                    "failed to restore the working directory to '" + prev.string() + "': " +
                    ec.message() + ". The process is in an unknown directory and further "
                    "calls are unsafe.";
            }
            catch (...)
            {
                // prev.string() can itself throw on windows for a path the active code page
                // cannot represent; a message without the path still beats silence
                try { g_cwd_restore_error = "failed to restore the working directory; the "
                                            "process is in an unknown directory and further "
                                            "calls are unsafe."; }
                catch (...) { /* nothing left to do */ }
            }
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
        if (!g_cwd_restore_error.empty())                                      \
            bad_state(g_cwd_restore_error);                                    \
        ScopedWorkingDir _swd(s->working_dir);

/* capi_error is caught first so a throw site's chosen status survives; everything else is a
   plain failure. Both bands clamp the message to PESTPP_MESSAGE_LEN, because that constant is
   exported as "big enough for anything this writes" and a run-on what() must not make a liar
   of it. */
#define CAPI_END()                                                             \
    } catch (const capi_error& e) {                                            \
        s->last_error = clamp_message(e.what());                               \
        return e.status;                                                       \
    } catch (const std::exception& e) {                                        \
        s->last_error = clamp_message(e.what());                               \
        return PESTPP_ERROR;                                                   \
    } catch (...) {                                                            \
        s->last_error = "unknown error";                                       \
        return PESTPP_ERROR;                                                   \
    }

extern "C" {

pestpp_status pestpp_create(const pestpp_create_options* opts, pestpp_handle* out)
{
    g_create_error.clear();
    if ((out == nullptr) || (opts == nullptr))
    {
        g_create_error = "null argument to pestpp_create";
        return PESTPP_ERROR;
    }
    *out = nullptr;
    // struct_size is how a caller built against an older header stays workable: anything past
    // the size it declares is a field it never knew about, so it is NOT READ - reading it
    // would run off the end of the caller's allocation - and takes its documented default
    // instead. HAVE() is that test, one field at a time.
    const size_t declared = (opts->struct_size < 0) ? 0u : (size_t)opts->struct_size;
#define HAVE(field) \
    (declared >= offsetof(pestpp_create_options, field) + sizeof(opts->field))
    if (!HAVE(tool))
    {
        g_create_error = "pestpp_create_options.struct_size not set; it must be "
                         "sizeof(pestpp_create_options)";
        return PESTPP_ERROR;
    }
    if (declared > sizeof(pestpp_create_options))
    {
        g_create_error = "pestpp_create_options.struct_size is larger than this library "
                         "understands; the caller was built against a newer header";
        return PESTPP_ERROR;
    }
    if ((!HAVE(ctl_file)) || (opts->ctl_file == nullptr))
    {
        g_create_error = "pestpp_create_options.ctl_file is required";
        return PESTPP_ERROR;
    }
    const char* ctl_file    = opts->ctl_file;
    const char* working_dir = HAVE(working_dir) ? opts->working_dir : nullptr;
    int rm_type             = HAVE(run_manager) ? opts->run_manager : PESTPP_RM_SERIAL;
    string panther_port     = (HAVE(panther_port) && (opts->panther_port != nullptr))
                              ? string(opts->panther_port) : string();
#undef HAVE
    if ((rm_type == PESTPP_RM_PANTHER) && panther_port.empty())
    {
        g_create_error = "PESTPP_RM_PANTHER requires panther_port";
        return PESTPP_ERROR;
    }
    if ((rm_type < PESTPP_RM_SERIAL) || (rm_type > PESTPP_RM_EXTERNAL))
    {
        g_create_error = "unknown run manager id";
        return PESTPP_ERROR;
    }
    if (!g_cwd_restore_error.empty())
    {
        g_create_error = g_cwd_restore_error;
        return PESTPP_ERROR;
    }
    unique_ptr<PestppSession> s(new PestppSession());
    try
    {
        s->tool = static_cast<pestpp_tool>(opts->tool);
        s->run_manager_type = rm_type;
        s->working_dir = (working_dir == nullptr) ? string() : string(working_dir);
        ScopedWorkingDir swd(s->working_dir);

        PestppOptions::ToolType tool_type;
        switch (s->tool)
        {
        case PESTPP_IES: tool_type = PestppOptions::ToolType::IES; break;
        case PESTPP_DA:  tool_type = PestppOptions::ToolType::DA;  break;
        case PESTPP_MOU: tool_type = PestppOptions::ToolType::MOU; break;
        case PESTPP_SQP: tool_type = PestppOptions::ToolType::SQP; break;
        default: bad_arg("unknown tool id");
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

        const ModelExecInfo& exi = s->pest_scenario->get_model_exec_info();
        const PestppOptions& opt = s->pest_scenario->get_pestpp_options();
        switch (rm_type)
        {
        case PESTPP_RM_SERIAL:
            // check_io() belongs to the serial branch only, exactly as in the executables: it
            // requires the template and instruction files to be present locally, which is true
            // where the model runs and need not be true for a master or an external driver
            s->pest_scenario->check_io(fout_rec);
            s->run_manager.reset(new RunManagerSerial(
                exi.comline_vec, exi.tplfile_vec, exi.inpfile_vec, exi.insfile_vec, exi.outfile_vec,
                s->file_manager->build_filename("rns"), ".",
                opt.get_max_run_fail(), opt.get_fill_tpl_zeros(),
                opt.get_additional_ins_delimiters(), opt.get_num_tpl_ins_threads(),
                opt.get_tpl_force_decimal(), opt.get_panther_echo()));
            break;

        case PESTPP_RM_PANTHER:
            // the master writes its own .rmr log, same as the executables do
            s->rmr_stream.reset(new ofstream(s->file_manager->build_filename("rmr")));
            s->run_manager.reset(new RunManagerPanther(
                s->file_manager->build_filename("rns"), panther_port, *s->rmr_stream,
                opt.get_max_run_fail(), opt.get_overdue_reched_fac(), opt.get_overdue_giveup_fac(),
                opt.get_overdue_giveup_minutes(), opt.get_panther_echo(),
                // Empty name vectors, matching every tool except pestpp-da: these become the
                // par/obs names the master asks each worker to validate, and switching that on
                // for a tool whose executable leaves it off would be a behavior change.
                (s->tool == PESTPP_DA) ? s->pest_scenario->get_ctl_ordered_par_names() : vector<string>{},
                (s->tool == PESTPP_DA) ? s->pest_scenario->get_ctl_ordered_obs_names() : vector<string>{},
                // and the four the executables pass that were previously left at their defaults
                opt.get_panther_timeout_milliseconds(),
                opt.get_panther_echo_interval_milliseconds(),
                opt.get_panther_persistent_workers(),
                opt.get_panther_ping_interval_secs()));
            break;

        case PESTPP_RM_EXTERNAL:
            s->run_manager.reset(new RunManagerExternal(
                exi.comline_vec, exi.tplfile_vec, exi.inpfile_vec, exi.insfile_vec, exi.outfile_vec,
                s->rns_filename = s->file_manager->build_filename("rns"),
                opt.get_max_run_fail()));
            break;

        default:
            bad_arg("unknown run manager id");
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
        default: bad_arg("unknown tool id");
        }
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
    // Stamp it dead FIRST. A second pestpp_destroy() on the same pointer then answers
    // PESTPP_INVALID_HANDLE instead of double-freeing, and so does any other call that races
    // in behind it.
    s->magic = 0;
    pestpp_status rc = PESTPP_OK;
    // The chdir is a convenience - some destructors write files relative to the working
    // directory - but it must never decide WHETHER we free. Failing to enter the directory
    // used to leak the whole session, including the panther socket, which then held its port.
    try { ScopedWorkingDir swd(s->working_dir); delete s; }
    catch (...)
    {
        try { delete s; } catch (...) { }
        rc = PESTPP_ERROR;
    }
    return rc;
}

const char* pestpp_last_error(pestpp_handle h)
{
    PestppSession* s = as_session(h);
    return (s == nullptr) ? "invalid handle" : s->last_error.c_str();
}

pestpp_status pestpp_get_last_error(pestpp_handle h, char* buf, int buf_len, int* needed)
{
    PestppSession* s = as_session(h);
    // deliberately NOT CAPI_BEGIN: this is the call you make after a failure, and clearing
    // the message on the way in would empty the very thing being asked for
    const string& msg = (s == nullptr) ? string() : s->last_error;
    if (needed != nullptr)
        *needed = (int)msg.size() + 1;
    if (s == nullptr)
        return PESTPP_INVALID_HANDLE;
    if (buf == nullptr)
        return PESTPP_OK;
    if (buf_len < (int)msg.size() + 1)
        return PESTPP_BUFFER_TOO_SMALL;
    for (size_t i = 0; i < msg.size(); i++)
        buf[i] = msg[i];
    buf[msg.size()] = '\0';
    return PESTPP_OK;
}

const char* pestpp_last_global_error(void) { return g_create_error.c_str(); }
const char* pestpp_last_create_error(void) { return g_create_error.c_str(); }

const char* pestpp_get_fatal_error(void) { return g_cwd_restore_error.c_str(); }

pestpp_status pestpp_clear_fatal_error(void)
{
    // The caller is asserting they have put the working directory back. We cannot verify that
    // - we no longer know where "back" was - so this is an acknowledgement, not a repair.
    try { g_cwd_restore_error.clear(); }
    catch (...) { return PESTPP_ERROR; }
    return PESTPP_OK;
}

pestpp_status pestpp_redirect_output(const char* path, int* redirect_token)
{
    g_create_error.clear();
    if ((path == nullptr) || (redirect_token == nullptr))
    {
        g_create_error = "pestpp_redirect_output needs a path and a place to put the token";
        return PESTPP_INVALID_ARGUMENT;
    }
    *redirect_token = 0;
    try
    {
        cout.flush();
        fflush(stdout);
        int sink = PESTPP_OPEN(path, PESTPP_O_FLAGS, PESTPP_O_MODE);
        if (sink < 0)
        {
            g_create_error = string("could not open '") + path + "' for output capture";
            return PESTPP_ERROR;
        }
        int saved = PESTPP_DUP(1);
        if ((saved < 0) || (PESTPP_DUP2(sink, 1) < 0))
        {
            PESTPP_CLOSE(sink);
            if (saved >= 0) PESTPP_CLOSE(saved);
            g_create_error = "could not redirect the library's output";
            return PESTPP_ERROR;
        }
        PESTPP_CLOSE(sink);
        RedirectRecord rec;
        rec.token = g_next_redirect_token++;
        rec.saved_fd = saved;
        g_redirect_stack.push_back(rec);
        *redirect_token = rec.token;
        return PESTPP_OK;
    }
    catch (const std::exception& e) { g_create_error = clamp_message(e.what()); return PESTPP_ERROR; }
    catch (...) { g_create_error = "unknown error in pestpp_redirect_output"; return PESTPP_ERROR; }
}

pestpp_status pestpp_restore_output(int redirect_token)
{
    g_create_error.clear();
    if (redirect_token <= 0)
        return PESTPP_OK;                 /* nothing was redirected */
    try
    {
        // Not on the stack at all: never issued, or already restored. Either way there is no
        // descriptor to put back and nothing safe to do with the number.
        bool known = false;
        for (size_t i = 0; i < g_redirect_stack.size(); i++)
            if (g_redirect_stack[i].token == redirect_token) { known = true; break; }
        if (!known)
        {
            g_create_error = "pestpp_restore_output was given a token this library never "
                             "handed out, or one that has already been restored";
            return PESTPP_INVALID_ARGUMENT;
        }
        if (g_redirect_stack.back().token != redirect_token)
        {
            stringstream ss;
            ss << "output redirects must be undone innermost first: token " << redirect_token
               << " is still covered by " << (g_redirect_stack.size() - 1)
               << " later redirect(s). Restoring it now would leave stdout pointing at "
                  "another session's log file permanently.";
            g_create_error = ss.str();
            return PESTPP_INVALID_STATE;
        }
        int saved = g_redirect_stack.back().saved_fd;
        g_redirect_stack.pop_back();
        cout.flush();
        fflush(stdout);
        if (PESTPP_DUP2(saved, 1) < 0)
        {
            // stdout is now whatever the redirect pointed at, and the caller has no way to
            // recover it. Say so loudly rather than returning a bare status.
            PESTPP_CLOSE(saved);
            g_create_error = "could not restore the library's output; stdout is still "
                             "redirected and the saved descriptor has been closed";
            return PESTPP_ERROR;
        }
        PESTPP_CLOSE(saved);
        return PESTPP_OK;
    }
    catch (const std::exception& e) { g_create_error = clamp_message(e.what()); return PESTPP_ERROR; }
    catch (...) { g_create_error = "unknown error in pestpp_restore_output"; return PESTPP_ERROR; }
}

pestpp_status pestpp_get_redirect_depth(int* depth)
{
    if (depth == nullptr)
        return PESTPP_INVALID_ARGUMENT;
    *depth = (int)g_redirect_stack.size();
    return PESTPP_OK;
}

pestpp_status pestpp_flush_output(void)
{
    g_create_error.clear();
    try
    {
        cout.flush();     // C++ stream -> the CRT's stdout
        cerr.flush();
        fflush(stdout);   // the CRT's buffer -> the OS handle
        return PESTPP_OK;
    }
    catch (const std::exception& e) { g_create_error = clamp_message(e.what()); return PESTPP_ERROR; }
    catch (...) { g_create_error = "unknown error in pestpp_flush_output"; return PESTPP_ERROR; }
}

pestpp_status pestpp_get_version(char* buf, int buf_len, int* needed)
{
    g_create_error.clear();
    try
    {
        // NB: the PESTPP_VERSION macro carries a trailing semicolon, so it only works as a
        // bare initializer - string(PESTPP_VERSION) does not compile
        string v = PESTPP_VERSION;
        if (needed != nullptr)
            *needed = (int)v.size() + 1;
        if (buf == nullptr)
            return PESTPP_OK;
        if (buf_len < (int)v.size() + 1)
        {
            g_create_error = "version buffer too small; call with buf=NULL to size it first";
            return PESTPP_BUFFER_TOO_SMALL;
        }
        for (size_t i = 0; i < v.size(); i++)
            buf[i] = v[i];
        buf[v.size()] = '\0';
        return PESTPP_OK;
    }
    catch (...) { g_create_error = "unknown error in pestpp_get_version"; return PESTPP_ERROR; }
}

pestpp_status pestpp_get_api_version(int* major, int* minor, int* patch)
{
    if (major != nullptr) *major = PESTPP_API_VERSION_MAJOR;
    if (minor != nullptr) *minor = PESTPP_API_VERSION_MINOR;
    if (patch != nullptr) *patch = PESTPP_API_VERSION_PATCH;
    return PESTPP_OK;
}

pestpp_status pestpp_get_run_manager(pestpp_handle h, int* run_manager)
{
    CAPI_BEGIN(h)
        if (run_manager == nullptr) bad_arg("null out-param");
        *run_manager = s->run_manager_type;
        return PESTPP_OK;
    CAPI_END()
}

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
            bad_state("initialization is already in progress; call "
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
            bad_state("nothing to finish; call pestpp_initialize_prepare() first");
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
            bad_state("initialization is incomplete: pestpp_initialize_finish() has "
                                "not been called since pestpp_initialize_prepare()");
        if (!s->initialized)
            bad_state("pestpp_initialize() must be called before pestpp_solve_iteration()");
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

pestpp_status pestpp_get_phi_summary(pestpp_handle h, int phi_type,
                                     double* mean, double* std, double* min, double* max)
{
    CAPI_BEGIN(h)
        double m = 0.0, sd = 0.0, mn = 0.0, mx = 0.0;
        s->adapter->phi_summary(phi_type, m, sd, mn, mx);
        if (mean != nullptr) *mean = m;
        if (std  != nullptr) *std  = sd;
        if (min  != nullptr) *min  = mn;
        if (max  != nullptr) *max  = mx;
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_get_iteration(pestpp_handle h, int* iter)
{
    CAPI_BEGIN(h)
        if (iter == nullptr) bad_arg("null out-param");
        *iter = s->adapter->iteration();
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_should_terminate(pestpp_handle h, int* out)
{
    CAPI_BEGIN(h)
        if (out == nullptr) bad_arg("null out-param");
        *out = s->adapter->should_terminate() ? 1 : 0;
        return PESTPP_OK;
    CAPI_END()
}

namespace {

/** Resolve an ensemble id to the live object on the tool. */
Ensemble* pick_ensemble(PestppSession* s, int id)
{
    if ((id < PESTPP_PAR_EN) || (id > PESTPP_WEIGHTS_EN))
        bad_arg("unknown ensemble id");
    Ensemble* e = s->adapter->ensemble(id);
    if (e == nullptr)
        unsupported(string("tool '") + s->adapter->name() +
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
        too_small("name buffer too small; call with buf=NULL to size it first");
    for (size_t i = 0; i < names.size(); i++)
        pack_one_name(names[i], buf + (i * PESTPP_NAME_LEN));
    return PESTPP_OK;
}

} // namespace

pestpp_status pestpp_get_ensemble_view(pestpp_handle h, int ensemble_id,
                                       double** data, int* nrow, int* ncol,
                                       int* view_token)
{
    CAPI_BEGIN(h)
        if ((data == nullptr) || (nrow == nullptr) || (ncol == nullptr))
            bad_arg("null out-param");
        Ensemble* ens = pick_ensemble(s, ensemble_id);
        Eigen::MatrixXd* m = ens->get_eigen_ptr_4_mod();
        *data = m->data();
        *nrow = (int)m->rows();
        *ncol = (int)m->cols();
        if (view_token != nullptr)
        {
            // Record the buffer identity so pestpp_view_is_valid() can answer from evidence -
            // EnsembleView re-checks the address and the dimensions, which is what catches a
            // reallocation nobody told us about.
            int tok = s->next_view_token++;
            ViewRecord rec;
            rec.view.reset(new EnsembleView(*ens));
            rec.ensemble_id = ensemble_id;
            s->views[tok] = std::move(rec);
            *view_token = tok;
        }
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_view_is_valid(pestpp_handle h, int view_token, int* out)
{
    CAPI_BEGIN(h)
        if (out == nullptr)
            bad_arg("null out-param");
        map<int, ViewRecord>::const_iterator it = s->views.find(view_token);
        // An unknown token is not an error: released and never-issued both mean "that pointer
        // is not something you may use", which is exactly what 0 says.
        *out = ((it == s->views.end()) || (!it->second.view->valid())) ? 0 : 1;
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_release_view(pestpp_handle h, int view_token)
{
    CAPI_BEGIN(h)
        s->views.erase(view_token);       // idempotent by construction
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
        if (tstat == nullptr) bad_arg("null out-param");
        *tstat = (int)s->adapter->par_ensemble()->get_trans_status();
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_set_option(pestpp_handle h, const char* key, const char* value)
{
    CAPI_BEGIN(h)
        if ((key == nullptr) || (value == nullptr)) bad_arg("null argument");
        PestppOptions::ARG_STATUS st =
            s->pest_scenario->get_pestpp_options_ptr()->set_option(key, value);
        // fall through to the control data section, which exposes the same interface. noptmax
        // is the most-set quantity in pest and it lives there, not in PestppOptions, so a
        // caller that could not reach it would be missing the obvious thing.
        if (st == PestppOptions::ARG_STATUS::ARG_NOTFOUND)
            st = s->pest_scenario->get_control_info_4_mod().set_option(key, value);
        if (st == PestppOptions::ARG_STATUS::ARG_NOTFOUND)
            bad_arg(string("unknown option '") + key + "'");
        if (st == PestppOptions::ARG_STATUS::ARG_INVALID)
            bad_arg(string("invalid value for option '") + key + "': " + value);
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_get_option(pestpp_handle h, const char* key,
                                char* buf, int buf_len, int* needed, int* found)
{
    CAPI_BEGIN(h)
        if (key == nullptr) bad_arg("null argument");
        // "" is a legitimate value for a ++ option, so an empty string cannot be the signal
        // for "no such option" - ask the registry whether the key exists at all, separately
        // from what it holds.
        const PestppOptions& ppo = s->pest_scenario->get_pestpp_options();
        bool known = ppo.is_valid_arg(key);
        string v = known ? ppo.get_option(key) : string();
        if (!known)
        {
            // fall through to the * control data section, which exposes the same interface.
            // Its get_option() is an if-chain returning "" for anything it does not recognise,
            // and every key it DOES recognise formats to at least one character (a number, or
            // a pestmode word), so non-empty is a sound test for "known" here.
            v = s->pest_scenario->get_control_info().get_option(key);    // e.g. NOPTMAX
            known = !v.empty();
        }
        if (found != nullptr)
            *found = known ? 1 : 0;
        if (!known)
        {
            if (needed != nullptr)
                *needed = 1;                     // just the NUL
            if (buf != nullptr && buf_len > 0)
                buf[0] = '\0';
            // A caller that passed `found` is probing on purpose, so absence is an answer.
            // One that did not gets an error rather than an empty string it would read as a
            // value - which is also what pestpp_set_option does with an unknown key.
            if (found != nullptr)
                return PESTPP_OK;
            bad_arg(string("unknown option '") + key + "'; pass a non-NULL `found` if you "
                    "meant to test whether it exists");
        }
        if (needed != nullptr)
            *needed = (int)v.size() + 1;
        if (buf == nullptr)
            return PESTPP_OK;
        if (buf_len < (int)v.size() + 1)
            too_small("option buffer too small; call with buf=NULL to size it first");
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
        unsupported("this call needs a PANTHER run manager - the serial manager has no "
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
    bad_arg("unknown run status");
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
        too_small("run id buffer too small; call with run_ids=NULL to size it first");
    for (size_t i = 0; i < ids.size(); i++)
        run_ids[i] = ids[i];
    return PESTPP_OK;
}

} // namespace

pestpp_status pestpp_get_phi_vector(pestpp_handle h, int phi_type,
                                    double* phi, char* names, int max_n, int* n_out)
{
    CAPI_BEGIN(h)
        vector<string> rnames;
        vector<double> vals;
        s->adapter->phi_vector(phi_type, rnames, vals);
        if (n_out != nullptr)
            *n_out = (int)vals.size();
        if ((phi == nullptr) && (names == nullptr))
            return PESTPP_OK;                 /* count-only call */
        if (max_n < (int)vals.size())
            too_small("phi buffers too small; call with NULL arrays to size first");
        for (size_t i = 0; i < vals.size(); i++)
        {
            if (phi != nullptr)   phi[i] = vals[i];
            if (names != nullptr) pack_one_name(rnames[i], names + (i * PESTPP_NAME_LEN));
        }
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_get_phi_residuals(pestpp_handle h, int phi_type,
                                       double* data, int max_nrow, int max_ncol,
                                       int* nrow, int* ncol,
                                       char* row_names, char* col_names)
{
    CAPI_BEGIN(h)
        vector<string> rnames, cnames;
        vector<double> vals;
        s->adapter->phi_residuals(phi_type, rnames, cnames, vals);
        int nr = (int)rnames.size(), nc = (int)cnames.size();
        if (nrow != nullptr) *nrow = nr;
        if (ncol != nullptr) *ncol = nc;
        // Size-only call: report the shape and touch nothing. The name buffers are sized from
        // nrow/ncol by the caller, so writing them here - before the caller knows how big to
        // make them - would overflow.
        if ((data == nullptr) && (row_names == nullptr) && (col_names == nullptr))
            return PESTPP_OK;
        // Everything is bounds-checked BEFORE the first write. max_nrow/max_ncol describe the
        // caller's allocation for the names as well as the data.
        if ((max_nrow < nr) || (max_ncol < nc))
            too_small("phi residual buffers too small; call with all pointers NULL "
                                "to size them first");
        if (row_names != nullptr)
            for (int i = 0; i < nr; i++)
                pack_one_name(rnames[i], row_names + (i * PESTPP_NAME_LEN));
        if (col_names != nullptr)
            for (int j = 0; j < nc; j++)
                pack_one_name(cnames[j], col_names + (j * PESTPP_NAME_LEN));
        if (data != nullptr)
            // stride is the CALLER's row count, not ours: an over-allocated buffer is read
            // back at max_nrow, so writing at nr would silently transpose the block
            for (int j = 0; j < nc; j++)
                for (int i = 0; i < nr; i++)
                    data[i + (j * max_nrow)] = vals[i + ((size_t)j * nr)];
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_get_obs_groups(pestpp_handle h, char* buf, int buf_len, int* count)
{
    CAPI_BEGIN(h)
        // aligned with the observation ensemble's columns, so a caller can zip the two
        vector<string> onames = s->adapter->obs_ensemble()->get_var_names();
        const Observations& ctl_obs = s->pest_scenario->get_ctl_observations();
        const ObservationInfo* oi = s->pest_scenario->get_ctl_observation_info_ptr();
        vector<string> groups;
        for (auto& n : onames)
        {
            // see pestpp_get_obs_weights: an unchecked find() here is a crash, not an error
            if (ctl_obs.find(n) == ctl_obs.end())
                throw runtime_error("observation '" + n + "' is in the ensemble but not in "
                                    "the control file");
            groups.push_back(oi->get_group(n));
        }
        return pack_names(groups, buf, buf_len, count);
    CAPI_END()
}

pestpp_status pestpp_supports_live_control(pestpp_handle h, int* out)
{
    CAPI_BEGIN(h)
        if (out == nullptr) bad_arg("null out-param");
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
        if (all_done == nullptr) bad_arg("null out-param");
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
            too_small("run state buffers too small; call with NULL arrays to size first");
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
            bad_arg("pestpp_cancel_runs needs at least one run id");
        int c = panther(s)->cancel_runs(vector<int>(run_ids, run_ids + n));
        if (n_cancelled != nullptr)
            *n_cancelled = c;
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_get_worker_count(pestpp_handle h, int* n)
{
    CAPI_BEGIN(h)
        if (n == nullptr) bad_arg("null out-param");
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
            bad_arg("worker index out of range");
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
            bad_arg("worker index out of range");
        const PantherWorkerState& w = ws[idx];
        switch (which)
        {
        case PESTPP_WORKER_COMPLETED: return emit_ids(w.completed_runs, run_ids, max_n, n_out);
        case PESTPP_WORKER_FAILED:    return emit_ids(w.failed_runs,    run_ids, max_n, n_out);
        case PESTPP_WORKER_TIMED_OUT: return emit_ids(w.timed_out_runs, run_ids, max_n, n_out);
        default: bad_arg("unknown worker history selector");
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

pestpp_status pestpp_get_obs_weights(pestpp_handle h, double* weights, int max_n, int* n_out)
{
    CAPI_BEGIN(h)
        vector<string> onames = s->adapter->obs_ensemble()->get_var_names();
        if (n_out != nullptr)
            *n_out = (int)onames.size();
        if (weights == nullptr)
            return PESTPP_OK;
        if (max_n < (int)onames.size())
            too_small("weight buffer too small; call with weights=NULL to size it");
        // ObservationInfo::get_weight/get_group deref find() without checking end(), which
        // is a segfault rather than an exception - and these names come from an ENSEMBLE, so
        // a weights csv whose columns do not match the scenario would crash the host process
        const Observations& ctl_obs = s->pest_scenario->get_ctl_observations();
        const ObservationInfo* oi = s->pest_scenario->get_ctl_observation_info_ptr();
        for (size_t i = 0; i < onames.size(); i++)
        {
            if (ctl_obs.find(onames[i]) == ctl_obs.end())
                throw runtime_error("observation '" + onames[i] + "' is in the ensemble but "
                                    "not in the control file");
            weights[i] = oi->get_weight(onames[i]);
        }
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_set_obs_weights(pestpp_handle h, const char* names,
                                     const double* weights, int n)
{
    CAPI_BEGIN(h)
        if ((names == nullptr) || (weights == nullptr) || (n <= 0))
            bad_arg("pestpp_set_obs_weights needs names and values");
        ObservationInfo& oi = s->pest_scenario->get_ctl_observation_info_4_mod();
        vector<string> want = unpack_names(names, n);
        // validate everything before changing anything, so a bad name cannot leave the
        // weights half-updated
        const Observations& obs = s->pest_scenario->get_ctl_observations();
        for (auto& nm : want)
            if (obs.find(nm) == obs.end())
                bad_arg("no such observation: '" + nm + "'");
        for (int i = 0; i < n; i++)
        {
            if (weights[i] < 0.0)
                bad_arg("negative weight for '" + want[i] + "'");
            oi.set_weight(want[i], weights[i]);
        }
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_broadcast_weights(pestpp_handle h)
{
    CAPI_BEGIN(h)
        Ensemble* w = s->adapter->ensemble(PESTPP_WEIGHTS_EN);
        if (w == nullptr)
            unsupported(string("tool '") + s->adapter->name() +
                                "' has no weights ensemble");
        if (w->shape().first == 0)
            return PESTPP_OK;                 /* not built yet; initialize will read the vector */
        const ObservationInfo* oi = s->pest_scenario->get_ctl_observation_info_ptr();
        vector<string> wnames = w->get_var_names();
        const Observations& ctl_obs = s->pest_scenario->get_ctl_observations();
        Eigen::VectorXd wvec(wnames.size());
        for (size_t j = 0; j < wnames.size(); j++)
        {
            if (ctl_obs.find(wnames[j]) == ctl_obs.end())
                throw runtime_error("weights ensemble column '" + wnames[j] + "' is not an "
                                    "observation in the control file");
            wvec[j] = oi->get_weight(wnames[j]);
        }
        Eigen::MatrixXd* m = w->get_eigen_ptr_4_mod();
        for (int i = 0; i < m->rows(); i++)
            m->row(i) = wvec;
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_update_phi(pestpp_handle h)
{
    CAPI_BEGIN(h)
        s->adapter->update_phi();
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_queue_runs(pestpp_handle h, int* n_queued)
{
    CAPI_BEGIN(h)
        if (s->pending_runs_valid)
            bad_state("runs are already queued; harvest them before queueing again");
        // The cycle has to be passed explicitly: it defaults to NULL_DA_CYCLE, so a da handle
        // driving its own queue/harvest would tag runs with no cycle at all while the in-tree
        // da path tagged them correctly. Invisible on a control file whose first cycle is 0.
        s->pending_runs = queue_ensemble_util(s->performance_log.get(),
                                              s->file_manager->rec_ofstream(),
                                              *s->adapter->par_ensemble(), s->run_manager.get(),
                                              false, vector<int>(), s->adapter->da_cycle());
        s->pending_runs_valid = true;
        if (n_queued != nullptr)
            *n_queued = (int)s->pending_runs.size();
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_process_runs(pestpp_handle h, int* n_failed)
{
    CAPI_BEGIN(h)
        if (!s->pending_runs_valid)
            bad_state("no queued runs to harvest; call pestpp_queue_runs first");
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
            bad_arg("no such realization: '" + n + "'");
        keep_p.push_back(pnames[it->second]);
        if (it->second < (int)onames.size())
            keep_o.push_back(onames[it->second]);
    }
    if (keep_p.empty())
        bad_arg("refusing to leave the ensemble empty");

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
            bad_arg("pestpp_keep_realizations needs at least one name");
        keep_across_coupled(s, unpack_names(names, n));
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_drop_realizations(pestpp_handle h, const char* names, int n)
{
    CAPI_BEGIN(h)
        if ((names == nullptr) || (n <= 0))
            bad_arg("pestpp_drop_realizations needs at least one name");
        vector<string> drop = unpack_names(names, n);
        vector<string> current = s->adapter->par_ensemble()->get_real_names();
        set<string> present(current.begin(), current.end());
        // name a bad request rather than silently dropping fewer than asked
        for (auto& nm : drop)
            if (present.find(nm) == present.end())
                bad_arg("no such realization: '" + nm + "'");
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
        // size-only: report the shape without writing into buffers the caller has not sized
        if ((data == nullptr) && (row_names == nullptr) && (col_names == nullptr))
            return PESTPP_OK;
        if ((max_nrow < nr) || (max_ncol < nc))
            too_small("snapshot buffers too small; call with all pointers NULL to "
                                "size them first");
        if (row_names != nullptr)
            for (int i = 0; i < nr; i++)
                pack_one_name(snap.row_names[i], row_names + (i * PESTPP_NAME_LEN));
        if (col_names != nullptr)
            for (int j = 0; j < nc; j++)
                pack_one_name(snap.col_names[j], col_names + (j * PESTPP_NAME_LEN));
        // column-major at the CALLER's stride, matching the zero-copy view and numpy order="F"
        if (data != nullptr)
            for (int j = 0; j < nc; j++)
                for (int i = 0; i < nr; i++)
                    data[i + (j * max_nrow)] = snap.values(i, j);
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_set_par_snapshot(pestpp_handle h, const double* data, int nrow, int ncol,
                                      const char* row_names, const char* col_names)
{
    CAPI_BEGIN(h)
        if ((data == nullptr) || (row_names == nullptr) || (col_names == nullptr))
            bad_arg("pestpp_set_par_snapshot needs data plus row and column names");
        if ((nrow <= 0) || (ncol <= 0))
            bad_arg("snapshot must have at least one row and one column");
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

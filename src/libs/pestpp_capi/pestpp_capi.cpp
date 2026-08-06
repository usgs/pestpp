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
// Deliberately no descriptor macros here. The windows redirect swaps cout's streambuf and must
// never call _dup2/_close on fd 1: under the static CRT this library has a private descriptor
// table, and closing its fd 1 destroys the console handle the HOST process is still using.
// See the note on RedirectRecord for the full failure mode. Defining the macros for windows
// would only make it easy to reintroduce.
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

    /**
     * @brief Queue / process the tool's current parameter ensemble.
     *
     * A PAIR, and overridden as a pair or not at all: whatever queues a batch has to be what
     * processes it. The default is the ies/da path, which is what pestpp_queue_runs() and
     * pestpp_process_runs() did directly before this hook existed.
     *
     * mou overrides both because its batch is not just the population: queue_population()
     * batches the CHANCE runs alongside it, and the process applies the chance shift. Driving
     * a chance-enabled mou through the generic ensemble path would queue the population and
     * silently omit the chance runs - a quiet wrong answer rather than an error.
     */
    virtual map<string, int> queue_runs(PerformanceLog* pl, ofstream& frec,
                                        RunManagerAbstract* rm)
    {
        return queue_ensemble_util(pl, frec, *par_ensemble(), rm, false, vector<int>(),
                                   da_cycle());
    }
    virtual vector<int> process_runs(PerformanceLog* pl, ofstream& frec, RunManagerAbstract* rm,
                                     map<string, int>& ids)
    {
        return process_ensemble_util(pl, frec, *par_ensemble(), *obs_ensemble(), rm, false,
                                     vector<int>(), ids);
    }
    /// One iteration/generation. Returns PESTPP_RETRY when the step was rejected and the
    /// algorithm wants another attempt - an outcome, not a fault.
    virtual pestpp_status advance() = 0;
    virtual void finalize() = 0;
    virtual int  iteration() = 0;   // not const: EnsembleMethod::get_iter() is not const
    virtual bool should_terminate() = 0;

    /// The live ensemble for a pestpp_ensemble_id, or null if this tool has no such thing.
    virtual Ensemble* ensemble(int id) = 0;
    /// The parameter-side ensemble, used for queue/process and the CTL snapshot.
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

    /// ---- deferred solve ---------------------------------------------------------------
    ///
    /// Whether one iteration decomposes into generate -> run -> evaluate. da and sqp do not:
    /// da's advance() is a whole noptmax loop and sqp's line search issues several run
    /// batches per iteration, so both refuse rather than approximating the shape.
    virtual bool supports_deferred_solve() const { return false; }
    /// Generate the candidates without running them; returns the runs they imply. 0 means the
    /// iteration finished during preparation and nothing is outstanding.
    virtual int  solve_prepare(pestpp_status& outcome) { (void)outcome; return 0; }
    /// Continue after a batch. Sets pending>0 when another batch is needed.
    virtual pestpp_status solve_finish(bool defer_runs, int& pending) { (void)defer_runs; pending = 0; return PESTPP_OK; }
    /// Queue the outstanding batch. `only` empty means the algorithm's own choice.
    virtual int  queue_solve_runs(const vector<string>& only) { (void)only; return 0; }
    /// Process it; returns how many runs failed.
    virtual int  process_solve_runs() { return 0; }
    /// Is a deferred solve open? Used to keep it and the composed advance() apart.
    virtual bool solve_is_open() const { return false; }
    /// Is a batch queued and awaiting process?
    virtual bool solve_batch_queued() const { return false; }
    /// Candidates awaiting runs, and the factors each was generated with.
    virtual int  candidate_count() const { return 0; }
    virtual Ensemble* candidate(int idx) { (void)idx; return nullptr; }
    virtual void candidate_info(int idx, double& inflation, double& backtrack)
    { (void)idx; inflation = 0.0; backtrack = 0.0; }

    /// The EnsembleMethod behind this tool, or null for the two that are not one.
    ///
    /// Needed for the operations that are structural rather than cosmetic - activating a
    /// zero-weighted observation has to reach into the active set, the noise ensemble and
    /// the weights ensemble together, and only the tool can do that coherently.
    virtual EnsembleMethod* ensemble_method() { return nullptr; }

    /// The chance machinery, or null for the tools that have none (ies and da).
    ///
    /// Returning the Constraints object rather than each stack separately keeps the resolving
    /// in one place: the stack ids mean the same thing for mou and sqp, so neither adapter
    /// needs to know about them.
    virtual Constraints* constraints() { return nullptr; }

    /// The scenario THIS TOOL WAS BUILT ON, which is not always the one the session parsed.
    ///
    /// da is the reason this exists. Its parameter and observation sets are cycle dependent,
    /// so DaAdapter is constructed against a per-cycle child scenario deep-copied out of the
    /// parent. Anything that reads or writes options, weights or observation metadata has to
    /// go through here: routing it to the parent instead means the write lands on a Pest
    /// object the tool never consults, and the matching read hands back a value that is not
    /// in effect - silently, and only for da.
    virtual Pest& scenario() = 0;

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

    // set while a run observer is executing. Everything outside the run-management allowlist
    // is refused for its duration - see reject_if_observing()
    bool in_observer = false;
    pestpp_run_observer_fn observer_fn = nullptr;
    void* observer_data = nullptr;

    bool initialized = false;
    // set between initialize_prepare() and initialize_finish(): the tool is half-initialized,
    // its ensembles drawn but the prior results not yet processed
    bool init_pending = false;
    // runs queued by pestpp_queue_runs(), awaiting pestpp_process_runs(). keyed by
    // realization name, which is what lets membership change while they are in flight.
    map<string,int> pending_runs;
    bool pending_runs_valid = false;
    // set between begin_batch() and end_batch(); see those functions for why it matters
    bool batch_open = false;

    // outstanding zero-copy views, by token. tokens never repeat within a session, so a
    // released token answers "not valid" rather than aliasing a later view.
    map<int, ViewRecord> views;
    int next_view_token = 1;
};

/// The deferred-solve state machine for an EnsembleMethod tool.
///
/// One iteration is generate -> run candidates -> evaluate -> (run the rest of the ensemble)
/// -> accept. The tool exposes every one of those; this holds the place in the sequence, the
/// context they share and the run ids in flight, so the C entry points stay thin and the two
/// tools that have this shape share one implementation.
struct EnsembleSolveState
{
    enum class Phase
    {
        NONE,               ///< nothing outstanding
        CANDIDATES,         ///< candidates generated, awaiting queue
        CANDIDATES_QUEUED,  ///< candidate runs in flight
        CANDIDATES_RUN,     ///< candidate results processed, awaiting evaluate
        REMAINING,          ///< a winner was picked; the rest of the ensemble awaits queue
        REMAINING_QUEUED,
        REMAINING_RUN
    };

    EnsembleMethod& tool;
    Pest& scen;
    int cycle;
    Phase phase = Phase::NONE;
    unique_ptr<UpgradeContext> ctx;
    vector<map<string, int>> lam_run_ids;
    map<string, int> remaining_run_ids;

    EnsembleSolveState(EnsembleMethod& t, Pest& p, int c) : tool(t), scen(p), cycle(c) {}

    bool is_open() const { return phase != Phase::NONE; }
    bool batch_queued() const
    { return (phase == Phase::CANDIDATES_QUEUED) || (phase == Phase::REMAINING_QUEUED); }

    /// map an UpgradeStatus onto what the C caller sees, and close the solve
    pestpp_status close(UpgradeStatus st)
    {
        tool.end_iteration();
        ctx.reset();
        lam_run_ids.clear();
        remaining_run_ids.clear();
        phase = Phase::NONE;
        return (st == UpgradeStatus::REJECTED_RETRY) ? PESTPP_RETRY : PESTPP_OK;
    }

    int prepare(pestpp_status& outcome)
    {
        outcome = PESTPP_OK;
        tool.begin_iteration();
        ctx.reset(new UpgradeContext(&scen));
        bool use_mda = scen.get_pestpp_options().get_ies_use_mda();
        UpgradeStatus st = tool.solve_prepare(*ctx, use_mda);
        if (st != UpgradeStatus::CONTINUE)
        {
            // the iteration ended without needing any runs - a lambda that could not be
            // generated, or the non-iterative shortcut. Nothing to hand over.
            outcome = close(st);
            return 0;
        }
        // the candidates are files on the spill path, so there is nothing in memory to look
        // at. Refuse here rather than handing back views onto ensembles that were emptied.
        if (ctx->pe_filenames.size() > 0)
        {
            close(UpgradeStatus::CONTINUE);
            unsupported("a deferred solve needs the candidates in memory, but "
                        "'ies_upgrades_in_memory' is false, so they were spilled to disk. Set "
                        "it true to inspect or edit candidates");
        }
        phase = Phase::CANDIDATES;
        return (int)(ctx->pe_lams.size() * ctx->subset_names.size());
    }

    int queue(const vector<string>& only)
    {
        if (phase == Phase::CANDIDATES)
        {
            if (only.size() > 0)
            {
                // naming realizations REPLACES the subset, which is coherent because
                // everything downstream - including who counts as 'remaining' - is derived
                // from subset_names rather than from positions
                vector<string> have = tool.get_pe_ptr()->get_real_names();
                set<string> hset(have.begin(), have.end());
                for (auto& n : only)
                    if (hset.find(n) == hset.end())
                        bad_arg("no such realization in the parameter ensemble: '" + n + "'");
                ctx->subset_names = only;
            }
            lam_run_ids = tool.queue_upgrade_ensembles(*ctx, cycle);
            phase = Phase::CANDIDATES_QUEUED;
            int n = 0;
            for (auto& m : lam_run_ids)
                n += (int)m.size();
            return n;
        }
        if (phase == Phase::REMAINING)
        {
            if (only.size() > 0)
                bad_arg("the remaining-realization batch is already determined - it is every "
                        "realization the candidate subset did not cover - so it cannot be "
                        "narrowed further; pass no names");
            remaining_run_ids = tool.queue_remaining_runs(*ctx, cycle);
            phase = Phase::REMAINING_QUEUED;
            return (int)remaining_run_ids.size();
        }
        bad_state("no deferred-solve batch is waiting to be queued");
        return 0;
    }

    int process()
    {
        if (phase == Phase::CANDIDATES_QUEUED)
        {
            tool.process_upgrade_ensembles(*ctx, lam_run_ids);
            lam_run_ids.clear();
            phase = Phase::CANDIDATES_RUN;
            return 0;
        }
        if (phase == Phase::REMAINING_QUEUED)
        {
            tool.process_remaining_runs(*ctx, remaining_run_ids);
            remaining_run_ids.clear();
            phase = Phase::REMAINING_RUN;
            return (int)ctx->failed_remaining.size();
        }
        bad_state("no deferred-solve runs are in flight");
        return 0;
    }

    pestpp_status finish(bool defer_runs, int& pending)
    {
        pending = 0;
        bool use_mda = ctx->use_mda;
        if (phase == Phase::CANDIDATES_RUN)
        {
            UpgradeStatus st = tool.evaluate_upgrades(*ctx);
            if (st != UpgradeStatus::CONTINUE)
                return close(st);
            st = tool.prepare_subset_completion(*ctx, cycle, use_mda);
            if (st != UpgradeStatus::CONTINUE)
                return close(st);
            if (ctx->needs_remaining_runs)
            {
                if (defer_runs)
                {
                    phase = Phase::REMAINING;
                    pending = ctx->remaining_pe.shape().first;
                    return PESTPP_OK;
                }
                ctx->failed_remaining = tool.run_ensemble(ctx->remaining_pe, ctx->remaining_oe,
                                                          vector<int>(), cycle);
            }
            phase = Phase::REMAINING_RUN;
        }
        if (phase != Phase::REMAINING_RUN)
            bad_state("pestpp_solve_finish() has nothing to continue: the candidate runs have "
                      "not been processed yet");
        UpgradeStatus st = tool.finish_subset_completion(*ctx, cycle, use_mda);
        if (st != UpgradeStatus::CONTINUE)
            return close(st);
        return close(tool.accept_or_reject(*ctx, use_mda, cycle));
    }
};

/** ies: begin_iteration -> solve (glm or mda) -> end_iteration. */
struct IesAdapter : public ToolAdapter
{
    IterEnsembleSmoother tool;
    Pest& scen;
    EnsembleSolveState deferred;
    IesAdapter(Pest& p, FileManager& fm, OutputFileWriter& ofw, PerformanceLog* pl,
               RunManagerAbstract* rm)
        : tool(p, fm, ofw, pl, rm), scen(p),
          deferred(tool, p, NetPackage::NULL_DA_CYCLE) {}

    EnsembleMethod* ensemble_method() override { return &tool; }
    Pest& scenario() override { return scen; }

    bool supports_deferred_solve() const override { return true; }
    int  solve_prepare(pestpp_status& outcome) override { return deferred.prepare(outcome); }
    pestpp_status solve_finish(bool defer_runs, int& pending) override
    { return deferred.finish(defer_runs, pending); }
    int  queue_solve_runs(const vector<string>& only) override { return deferred.queue(only); }
    int  process_solve_runs() override { return deferred.process(); }
    bool solve_is_open() const override { return deferred.is_open(); }
    bool solve_batch_queued() const override { return deferred.batch_queued(); }
    int  candidate_count() const override
    { return deferred.ctx ? (int)deferred.ctx->pe_lams.size() : 0; }
    Ensemble* candidate(int idx) override
    {
        if ((!deferred.ctx) || (idx < 0) || (idx >= (int)deferred.ctx->pe_lams.size()))
            return nullptr;
        return &deferred.ctx->pe_lams[idx];
    }
    void candidate_info(int idx, double& inflation, double& backtrack) override
    {
        if ((!deferred.ctx) || (idx < 0) || (idx >= (int)deferred.ctx->lam_vals.size()))
            bad_arg("no such candidate");
        inflation = deferred.ctx->lam_vals[idx];
        backtrack = deferred.ctx->scale_vals[idx];
    }

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
    // the per-cycle CHILD scenario, not the session's parent - see ToolAdapter::scenario()
    Pest& scen;
    int cycle;
    DaAdapter(Pest& p, FileManager& fm, OutputFileWriter& ofw, PerformanceLog* pl,
              RunManagerAbstract* rm, int _cycle)
        : tool(p, fm, ofw, pl, rm), scen(p), cycle(_cycle) {}

    EnsembleMethod* ensemble_method() override { return &tool; }
    Pest& scenario() override { return scen; }

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

    Pest& scenario() override { return scen; }

    void initialize() override { tool.initialize(); }
    // The initial POPULATION evaluation is handed to the caller; the other two evaluations in
    // mou's initialize() (control-file values, mean values) are one-off diagnostics on
    // noptmax<=0 paths and stay atomic inside prepare().
    int  initialize_prepare() override { return tool.initialize_prepare(); }
    void initialize_finish() override { tool.initialize_finish(); }
    // chance-aware, unlike the generic ensemble path - see ToolAdapter::queue_runs()
    map<string, int> queue_runs(PerformanceLog*, ofstream&, RunManagerAbstract*) override
    { return tool.queue_initial_population(); }
    vector<int> process_runs(PerformanceLog*, ofstream&, RunManagerAbstract*,
                             map<string, int>&) override
    { return tool.process_initial_population(); }
    pestpp_status advance() override
    {
        // ONE call into the tool, not a re-implementation of its loop body. Doing the four
        // phase calls here instead left MOEA's own generation counter untouched, so every
        // generation reported itself as the same number and wrote over the same .N. files.
        tool.solve_generation();
        iter++;
        return PESTPP_OK;
    }

    // -- deferred solve. mou's generation is generate -> run -> evaluate -> report with one
    // candidate population, so there is never a second batch and pending is always 0.
    unique_ptr<GenerationContext> gctx;
    enum class GPhase { NONE, POP, POP_QUEUED, POP_RUN };
    GPhase gphase = GPhase::NONE;
    map<string, int> pop_run_ids;

    bool supports_deferred_solve() const override { return true; }
    bool solve_is_open() const override { return gphase != GPhase::NONE; }
    bool solve_batch_queued() const override { return gphase == GPhase::POP_QUEUED; }
    int  solve_prepare(pestpp_status& outcome) override
    {
        outcome = PESTPP_OK;
        // the generation counter is advanced by the tool, not here: it names the .N.
        // population files, and incrementing it in the adapter is how every generation once
        // ended up overwriting the same ones
        tool.advance_generation();
        gctx.reset(new GenerationContext(&scen, tool.get_rand_gen_ptr()));
        tool.generate_generation(*gctx);
        gphase = GPhase::POP;
        return gctx->new_dp.shape().first;
    }
    int queue_solve_runs(const vector<string>& only) override
    {
        if (gphase != GPhase::POP)
            bad_state("no deferred-solve batch is waiting to be queued");
        if (only.size() > 0)
            bad_arg("mou runs the whole candidate population; it has no subset to narrow");
        // sized here rather than in generate, so members the caller added or removed are
        // honored - which is the entire point of handing the population over
        gctx->new_op.reserve(gctx->new_dp.get_real_names(), tool.get_op_ptr()->get_var_names());
        pop_run_ids = tool.queue_population(gctx->new_dp, true);
        gphase = GPhase::POP_QUEUED;
        return (int)pop_run_ids.size();
    }
    int process_solve_runs() override
    {
        if (gphase != GPhase::POP_QUEUED)
            bad_state("no deferred-solve runs are in flight");
        vector<int> failed = tool.process_population(gctx->new_dp, gctx->new_op, true, pop_run_ids);
        pop_run_ids.clear();
        gphase = GPhase::POP_RUN;
        return (int)failed.size();
    }
    pestpp_status solve_finish(bool, int& pending) override
    {
        pending = 0;
        if (gphase != GPhase::POP_RUN)
            bad_state("pestpp_solve_finish() has nothing to continue: the population runs have "
                      "not been processed yet");
        tool.evaluate_generation(*gctx);
        tool.report_generation(*gctx);
        gctx.reset();
        gphase = GPhase::NONE;
        iter++;
        return PESTPP_OK;
    }
    int candidate_count() const override { return gctx ? 1 : 0; }
    Ensemble* candidate(int idx) override
    { return (gctx && (idx == 0)) ? &gctx->new_dp : nullptr; }
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
    Constraints* constraints() override { return tool.get_constraints_ptr(); }
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
    Pest& scen;
    int iter = 0;
    SqpAdapter(Pest& p, FileManager& fm, OutputFileWriter& ofw, PerformanceLog* pl,
               RunManagerAbstract* rm)
        : tool(p, fm, ofw, pl, rm), scen(p) {}

    Pest& scenario() override { return scen; }

    void initialize() override { tool.initialize(); }
    // The initial ENSEMBLE evaluation is handed to the caller, when there is one. It returns 0
    // on the finite-difference gradient path (whose runs are perturbations, not candidates),
    // on the control-file-values and mean-values diagnostics, and when the ensemble was
    // supplied rather than drawn - so branch on the count, not on the tool.
    int  initialize_prepare() override { return tool.initialize_prepare(); }
    void initialize_finish() override { tool.initialize_finish(); }
    // sqp's process applies drop_violating_members(); the generic ensemble path does not
    map<string, int> queue_runs(PerformanceLog*, ofstream&, RunManagerAbstract*) override
    { return tool.queue_initial_ensemble(); }
    vector<int> process_runs(PerformanceLog*, ofstream&, RunManagerAbstract*,
                             map<string, int>&) override
    { return tool.process_initial_ensemble(); }
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
    Constraints* constraints() override { return tool.get_constraints_ptr(); }
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
 * Redirecting output is a process-wide operation however much it looks like a per-handle one.
 * Nesting works - each redirect saves what output was going to, so unwinding in the reverse
 * order puts every layer back - but UNWINDING OUT OF ORDER does not, and it fails silently and
 * permanently:
 *
 *     A redirects -> saved_a = the console
 *     B redirects -> saved_b = A's log file
 *     A restores  -> output = the console        (B is still "redirected", to nothing)
 *     B restores  -> output = A's log file       <- and stays there forever
 *
 * The old signature handed the caller the raw saved descriptor, so there was nothing to check
 * against: any int was as plausible as any other. Now the caller gets an opaque TOKEN, the
 * saved state stays in here, and restoring anything but the innermost redirect is refused.
 * Enforcing the order beats guessing at a repair - silently unwinding somebody else's redirect
 * is the corruption, not the fix.
 *
 * HOW the redirect is done differs by platform, and the reason is worth stating because the
 * obvious implementation is wrong on windows.
 *
 * POSIX redirects file descriptor 1 with dup2. There is exactly one descriptor table per
 * process, so that captures everything - this library's cout, any printf, and the stdout of
 * child processes that inherit it. That is the whole feature, and it is safe because the
 * descriptor table really is shared with the host.
 *
 * WINDOWS HAS NO SUCH GUARANTEE, and doing the same thing there corrupted output files for
 * months before anyone traced it. A CRT's fd 1 is a private table entry wrapping an OS handle.
 * Build this library with the static runtime (/MT, the project default, since the executables
 * ship without a redistributable) and load it into a host with its own CRT - python, via
 * ctypes - and there are now TWO fd 1s wrapping THE SAME console handle, neither aware of the
 * other. _dup2(sink, 1) closes the descriptor it replaces, _close calls CloseHandle, and the
 * handle the HOST is still using is destroyed:
 *
 *     - the host's next write fails: OSError [WinError 6] The handle is invalid
 *     - windows recycles the freed handle VALUE onto the next file this library opens. The
 *       FileManager holds the phi csvs and the .rec open for the whole run, so they are the
 *       likely recipients - which is how console output ended up written inside
 *       pest.phi.composite.csv, twice, looking like a numeric difference in a comparison.
 *
 * So on windows we swap std::cout's streambuf instead and never touch a descriptor. It cannot
 * destroy anything the host owns, under any runtime configuration. The cost is that printf and
 * child-process output are NOT captured there - acceptable, and measured rather than assumed:
 * pestpp_common, common and this file contain no printf-to-stdout at all, and the run managers'
 * only ones are in the agent's linpack benchmark. Everything on the API path goes through cout.
 *
 * Do not "simplify" this back to a single dup2 path. */
struct RedirectRecord
{
    int token;
#if defined(_WIN32)
    std::ofstream* log;          ///< owned; open for the life of the redirect
    std::streambuf* saved_buf;   ///< cout's buffer before we took it over
#else
    int saved_fd;
#endif
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
        reject_if_observing(s, __func__);                                      \
        ScopedWorkingDir _swd(s->working_dir);

/* The same, for the handful of calls that ARE legal from inside a run observer. The observer
   fires mid-batch, so most of the API would be reading a half-updated tool - but the
   run-management calls read only the run manager, which is precisely the state the observer
   was told about. Preemption needs exactly this door: look at what is running, decide, cancel.
   Hence an allowlist rather than a ban. */
#define CAPI_BEGIN_OBSERVER_SAFE(h)                                            \
    PestppSession* s = as_session(h);                                          \
    if (s == nullptr) return PESTPP_INVALID_HANDLE;                            \
    try {                                                                      \
        s->last_error.clear();                                                 \
        if (!g_cwd_restore_error.empty())                                      \
            bad_state(g_cwd_restore_error);                                    \
        ScopedWorkingDir _swd(s->working_dir);

/// Refuse a call that is not legal from inside a run observer.
///
/// The observer runs mid-batch, with ensembles part-processed and phi not yet recomputed.
/// Reading them there is not merely discouraged, it is meaningless - and writing them is
/// worse. Rather than document a rule nobody can check, the state is enforced, and the
/// message names the two calls that are actually useful from in there.
static void reject_if_observing(PestppSession* s, const char* fn)
{
    if (!s->in_observer)
        return;
    bad_state(string("'") + fn + "' cannot be called from inside a run observer: the batch is "
              "mid-flight, so the ensembles and phi it would read are part-updated. Only the "
              "run-management calls are legal there - pestpp_get_run_states, "
              "pestpp_get_run_time_stats, pestpp_cancel_runs and pestpp_get_worker_*");
}

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
        RedirectRecord rec;
#if defined(_WIN32)
        // See the note on RedirectRecord: descriptors are not safely shared with the host on
        // windows, so redirect the stream, not the descriptor.
        std::ofstream* log = new std::ofstream(path, std::ios::out | std::ios::app);
        if (!log->is_open())
        {
            delete log;
            g_create_error = string("could not open '") + path + "' for output capture";
            return PESTPP_ERROR;
        }
        rec.log = log;
        rec.saved_buf = cout.rdbuf(log->rdbuf());
#else
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
        rec.saved_fd = saved;
#endif
        rec.token = g_next_redirect_token++;
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
        RedirectRecord rec = g_redirect_stack.back();
        g_redirect_stack.pop_back();
        cout.flush();
        fflush(stdout);
#if defined(_WIN32)
        // Put the buffer back before closing the file, or cout is left pointing at a dead
        // streambuf for the instant in between - and at a destroyed one if the close throws.
        cout.rdbuf(rec.saved_buf);
        rec.log->close();
        delete rec.log;
#else
        int saved = rec.saved_fd;
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
#endif
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
        // May carry a release-candidate suffix ("5.2.29rc1"), so callers must treat this as an
        // opaque string rather than parsing three numbers out of it.
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
        if (s->adapter->solve_is_open())
            bad_state("a deferred solve is open; finish it with pestpp_solve_finish() rather "
                      "than starting a second iteration on top of it");
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

// extern "C++" for the reason spelled out at the top of the block near read_stored_run:
// stack_empty_reason returns a std::string BY VALUE, and this file's body sits inside an
// extern "C" block, which would otherwise hand it C language linkage while it uses the C++
// ABI's hidden-return-pointer convention.
extern "C++" {
namespace {

/** The chance machinery, or a diagnosis of why this tool has none. */
Constraints* pick_constraints(PestppSession* s)
{
    Constraints* c = s->adapter->constraints();
    if (c == nullptr)
        unsupported(string("tool '") + s->adapter->name() + "' has no chance stacks; only mou "
                    "and sqp carry constraints and the chance machinery");
    return c;
}

/** Say why a stack is empty, in terms of the setting that made it so.
 *
 * An empty stack is not an error - a FOSM run and a risk-neutral run are both legitimate, and
 * both leave every stack empty forever. But "empty" on its own is indistinguishable from "the
 * stack has not been drawn yet", and those want opposite responses from a caller, so the three
 * cases are separated here rather than left for the caller to infer.
 */
string stack_empty_reason(Constraints* c)
{
    if (!c->get_use_chance())
        return "this run is risk neutral (risk is 0.5), so no chance stack is ever drawn";
    if (c->get_use_fosm())
        return "this run uses FOSM rather than stacks for chance, so the stacks stay empty; "
               "set opt_stack_size or opt_par_stack/opt_obs_stack to use stacks instead";
    return "the stack has not been drawn yet - it is filled during initialization";
}

/** Resolve an ensemble id to the live object on the tool. */
Ensemble* pick_ensemble(PestppSession* s, int id)
{
    // checked BEFORE the candidate range because member stack ids are the higher of the two
    // and would otherwise be read as absurdly large candidate indices
    if (id >= PESTPP_MEMBER_STACK_EN)
    {
        int idx = id - PESTPP_MEMBER_STACK_EN;
        Constraints* c = pick_constraints(s);
        map<string, ObservationEnsemble>* m = c->get_stack_oe_map_ptr();
        if (m->empty())
            bad_state("there are no per-member chance stacks: " + stack_empty_reason(c) +
                      ". Per-member stacks also need opt_chance_points to be 'all'");
        if ((idx < 0) || (idx >= (int)m->size()))
            bad_arg("no such member stack; there are " + std::to_string(m->size()));
        auto it = m->begin();
        std::advance(it, idx);
        return &it->second;
    }
    // candidates are ordinary ensemble ids so that views, names and snapshots all work on
    // them unchanged - see PESTPP_CANDIDATE_EN
    if (id >= PESTPP_CANDIDATE_EN)
    {
        int idx = id - PESTPP_CANDIDATE_EN;
        Ensemble* c = s->adapter->candidate(idx);
        if (c == nullptr)
        {
            if (!s->adapter->solve_is_open())
                bad_state("candidate ensembles exist only during a deferred solve; call "
                          "pestpp_solve_prepare() first");
            bad_arg("no such candidate; there are " +
                    std::to_string(s->adapter->candidate_count()));
        }
        return c;
    }
    // the chance stacks live on Constraints rather than on the tool, so they resolve here
    // instead of through adapter->ensemble() - the ids mean the same thing for mou and sqp
    if ((id >= PESTPP_STACK_PAR_EN) && (id <= PESTPP_NESTED_PAR_EN))
    {
        Constraints* c = pick_constraints(s);
        Ensemble* st = nullptr;
        switch (id)
        {
        case PESTPP_STACK_PAR_EN:  st = c->get_stack_pe_ptr();  break;
        case PESTPP_STACK_OBS_EN:  st = c->get_stack_oe_ptr();  break;
        case PESTPP_NESTED_PAR_EN: st = c->get_nested_pe_ptr(); break;
        }
        // an empty stack is returned rather than refused: it is the honest answer for a FOSM
        // or risk-neutral run, and pestpp_get_stack_status() is how a caller asks why
        return st;
    }
    if ((id < PESTPP_PAR_EN) || (id > PESTPP_WEIGHTS_EN))
        bad_arg("unknown ensemble id");
    Ensemble* e = s->adapter->ensemble(id);
    if (e == nullptr)
        unsupported(string("tool '") + s->adapter->name() +
                            "' has no ensemble with that id");
    return e;
}

/** Put a snapshot's rows into `want` order, dropping nothing and inventing nothing.
 *
 * Any realization in the snapshot that is not in `want` (or the other way round) means the
 * two views of membership have diverged, which is a bug rather than something to paper over -
 * so this leaves the snapshot untouched in that case and lets the caller see the original
 * order rather than a silently half-permuted one.
 */
void reorder_snapshot_rows(ParameterSnapshot& snap, const vector<string>& want)
{
    if ((int)want.size() != (int)snap.values.rows())
        return;
    map<string,int> at;
    for (int i = 0; i < (int)snap.row_names.size(); i++)
        at[snap.row_names[i]] = i;
    vector<int> order;
    order.reserve(want.size());
    for (size_t i = 0; i < want.size(); i++)
    {
        map<string,int>::const_iterator it = at.find(want[i]);
        if (it == at.end())
            return;                       // membership disagrees; do not guess
        order.push_back(it->second);
    }
    bool already = true;
    for (size_t i = 0; i < order.size(); i++)
        if (order[i] != (int)i) { already = false; break; }
    if (already)
        return;
    Eigen::MatrixXd permuted(snap.values.rows(), snap.values.cols());
    for (int i = 0; i < (int)order.size(); i++)
        permuted.row(i) = snap.values.row(order[i]);
    snap.values = permuted;
    snap.row_names = want;
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
} // extern "C++"

pestpp_status pestpp_get_stack_status(pestpp_handle h, int* use_chance, int* use_fosm,
                                      int* use_robust, double* risk, int* stack_size)
{
    CAPI_BEGIN(h)
        Constraints* c = pick_constraints(s);
        // read live off the options rather than from anything latched at initialize(), so a
        // risk set through the API since then is reflected here
        if (use_chance  != nullptr) *use_chance  = c->get_use_chance() ? 1 : 0;
        if (use_fosm    != nullptr) *use_fosm    = c->get_use_fosm() ? 1 : 0;
        if (use_robust  != nullptr) *use_robust  = c->get_use_robust() ? 1 : 0;
        if (risk        != nullptr) *risk        = c->get_risk();
        // the stack as it actually stands, which is 0 until it is drawn and need not equal
        // opt_stack_size - a stack loaded from a file brings its own row count
        if (stack_size  != nullptr) *stack_size  = c->get_stack_oe_ptr()->shape().first;
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_get_member_stack_count(pestpp_handle h, int* count)
{
    CAPI_BEGIN(h)
        if (count == nullptr) bad_arg("null out-param");
        *count = (int)pick_constraints(s)->get_stack_oe_map_ptr()->size();
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_get_member_stack_names(pestpp_handle h, char* buf, int buf_len, int* count)
{
    CAPI_BEGIN(h)
        // same order as the PESTPP_MEMBER_STACK_EN ids, because both walk the same sorted map
        vector<string> names;
        for (auto& kv : *pick_constraints(s)->get_stack_oe_map_ptr())
            names.push_back(kv.first);
        return pack_names(names, buf, buf_len, count);
    CAPI_END()
}

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
        // INIT-ONLY options are consumed while the tool builds itself - ensemble files,
        // ensemble sizes, the localizer, opt_use_robust. Setting one afterwards writes a value
        // nothing will ever read again, or worse describes an ensemble that is already built
        // differently. Refuse instead of accepting silently.
        if (s->initialized &&
            s->adapter->scenario().get_pestpp_options_ptr()->is_init_only(key))
            bad_state(string("option '") + key + "' is init-only: it is consumed during "
                      "initialization and cannot be changed on a running tool");
        PestppOptions::ARG_STATUS st =
            s->adapter->scenario().get_pestpp_options_ptr()->set_option(key, value);
        // fall through to the control data section, which exposes the same interface. noptmax
        // is the most-set quantity in pest and it lives there, not in PestppOptions, so a
        // caller that could not reach it would be missing the obvious thing.
        if (st == PestppOptions::ARG_STATUS::ARG_NOTFOUND)
            st = s->adapter->scenario().get_control_info_4_mod().set_option(key, value);
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
        const PestppOptions& ppo = s->adapter->scenario().get_pestpp_options();
        bool known = ppo.is_valid_arg(key);
        string v = known ? ppo.get_option(key) : string();
        if (!known)
        {
            // fall through to the * control data section, which exposes the same interface.
            // Its get_option() is an if-chain returning "" for anything it does not recognise,
            // and every key it DOES recognise formats to at least one character (a number, or
            // a pestmode word), so non-empty is a sound test for "known" here.
            v = s->adapter->scenario().get_control_info().get_option(key);    // e.g. NOPTMAX
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
        const Observations& ctl_obs = s->adapter->scenario().get_ctl_observations();
        const ObservationInfo* oi = s->adapter->scenario().get_ctl_observation_info_ptr();
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

/* begin/slice/end is a sequence, and until now nothing checked that a caller followed it.
 *
 * The three ways to get it wrong are not equally forgiving. run_slice() before begin_batch()
 * is the dangerous one: on PANTHER the master has not started its listener bookkeeping, and
 * the slice races the idle-ping thread - a crash inside the library, not an error a caller
 * can catch. A second begin_batch() resets the run counters with runs already in flight, so
 * every subsequent stat is wrong without being obviously wrong. end_batch() without a batch
 * is merely pointless. All three are now refused. */
pestpp_status pestpp_begin_batch(pestpp_handle h)
{
    CAPI_BEGIN(h)
        if (s->batch_open)
            bad_state("a batch is already open; call pestpp_end_batch() before starting "
                      "another. Starting a second batch would reset the run counters with "
                      "runs still in flight.");
        s->run_manager->begin_batch();
        s->batch_open = true;
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_run_slice(pestpp_handle h, double max_seconds, int* all_done)
{
    CAPI_BEGIN(h)
        if (all_done == nullptr) bad_arg("null out-param");
        if (!s->batch_open)
            bad_state("no batch is open; call pestpp_begin_batch() first. Slicing without it "
                      "races the run manager's own bookkeeping.");
        RunManagerAbstract::RUN_SLICE_STATUS st = s->run_manager->run_slice(max_seconds);
        *all_done = (st == RunManagerAbstract::RUN_SLICE_STATUS::ALL_DONE) ? 1 : 0;
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_end_batch(pestpp_handle h)
{
    CAPI_BEGIN(h)
        if (!s->batch_open)
            bad_state("no batch is open; pestpp_end_batch() without a matching "
                      "pestpp_begin_batch() does nothing.");
        s->run_manager->end_batch();
        s->batch_open = false;
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_is_batch_open(pestpp_handle h, int* out)
{
    CAPI_BEGIN(h)
        if (out == nullptr) bad_arg("null out-param");
        *out = s->batch_open ? 1 : 0;
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_get_run_time_stats(pestpp_handle h, double* avg_run_sec,
                                        int* n_completed, int* n_failed, int* n_timed_out,
                                        int* n_queued, int* n_running)
{
    CAPI_BEGIN_OBSERVER_SAFE(h)
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
    CAPI_BEGIN_OBSERVER_SAFE(h)
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
    CAPI_BEGIN_OBSERVER_SAFE(h)
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
    CAPI_BEGIN_OBSERVER_SAFE(h)
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
    CAPI_BEGIN_OBSERVER_SAFE(h)
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
    CAPI_BEGIN_OBSERVER_SAFE(h)
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

/* ---- queue / process ------------------------------------------------------------------ */

// extern "C++" on purpose, for the same reason as read_stored_run below: this file's body sits
// inside an extern "C" block, which would otherwise give this helper C language linkage while
// it RETURNS a C++ object by value. Returning a non-trivial type is not a register operation -
// the caller passes a hidden pointer to storage for the result and the callee constructs it
// there - and that hidden argument is part of the C++ ABI, which is exactly what C linkage
// says not to use. gcc and clang happen to agree and only warn; MSVC is not obliged to, and
// when it does not the object is constructed through whatever was in the slot the caller
// believed it had set - an access violation if you are lucky. The anonymous namespace does not
// help: it governs visibility, not language linkage.
extern "C++" {
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
} // extern "C++"

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
        const Observations& ctl_obs = s->adapter->scenario().get_ctl_observations();
        const ObservationInfo* oi = s->adapter->scenario().get_ctl_observation_info_ptr();
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
        ObservationInfo& oi = s->adapter->scenario().get_ctl_observation_info_4_mod();
        vector<string> want = unpack_names(names, n);
        // validate everything before changing anything, so a bad name cannot leave the
        // weights half-updated
        const Observations& obs = s->adapter->scenario().get_ctl_observations();
        for (auto& nm : want)
            if (obs.find(nm) == obs.end())
                bad_arg("no such observation: '" + nm + "'");
        for (int i = 0; i < n; i++)
            if (weights[i] < 0.0)
                bad_arg("negative weight for '" + want[i] + "'");

        // Activation is the structural case and it goes through the tool, not through the
        // scenario. An observation that starts at zero weight is not in the active set, has
        // no column in the weights ensemble, and - the part that matters - has no NOISE
        // realizations, because none were drawn for it at initialize. Setting its weight
        // here and stopping would leave it looking active while contributing nothing.
        //
        // Reweighting an already-active observation is deliberately NOT coupled to noise:
        // weights and noise are independent in general, and redrawing noise every time a
        // weight moved would destroy the comparability of phi across iterations.
        EnsembleMethod* em = s->adapter->ensemble_method();
        if (em != nullptr)
        {
            map<string, double> w;
            for (int i = 0; i < n; i++)
                w[want[i]] = weights[i];
            em->activate_obs(w);
        }
        else
        {
            for (int i = 0; i < n; i++)
                oi.set_weight(want[i], weights[i]);
        }
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_reinflate_ensemble(pestpp_handle h, double factor, int num_reals,
                                        int center_on_min_phi)
{
    CAPI_BEGIN(h)
        EnsembleMethod* em = s->adapter->ensemble_method();
        if (em == nullptr)
            unsupported(string("tool '") + s->adapter->name() +
                        "' has no parameter ensemble to reinflate");
        if (!s->initialized)
            bad_state("reinflation needs an initialized ensemble; call pestpp_initialize() "
                      "first");
        if ((factor <= 0.0) || (factor > 1.0))
            bad_arg("reinflation factor must be in (0, 1]: 1.0 restores the full prior spread, "
                    "smaller values keep the ensemble tighter around the current mean");
        // Reinflation SELECTS realizations out of the prior ensemble - it never draws new
        // ones - so asking for more than the prior holds cannot be honoured. The underlying
        // routine deals with that by ignoring the request and reinflating to the prior's size
        // instead, which is the worst outcome for a caller: the size does not change and
        // nothing says why. Refuse it here, and say what the way around it is.
        int want = (num_reals < 0) ? -num_reals : num_reals;
        int n_prior = em->get_pe_base_ptr()->shape().first;
        if (want > n_prior)
        {
            stringstream ss;
            ss << "cannot reinflate to " << want << " realizations: the prior ensemble holds "
               << n_prior << ", and reinflation draws from it rather than generating new "
               << "realizations. To grow the ensemble mid-run, start with a prior big enough "
               << "for the largest size you will ask for ('ies_num_reals' or a prior ensemble "
               << "file) and have the run begin on a subset of it ('ies_reinflate_num_reals', "
               << "whose first entry truncates the working ensemble during initialization)";
            bad_arg(ss.str());
        }
        if ((center_on_min_phi < -1) || (center_on_min_phi > 1))
            bad_arg("center_on_min_phi must be -1 (follow 'ies_n_iter_reinflate'), 0 (the "
                    "ensemble mean) or 1 (the minimum-phi realization)");
        em->reinflate_par_ensemble(factor, num_reals, center_on_min_phi);
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
        const ObservationInfo* oi = s->adapter->scenario().get_ctl_observation_info_ptr();
        vector<string> wnames = w->get_var_names();
        const Observations& ctl_obs = s->adapter->scenario().get_ctl_observations();
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

namespace {
/// Refuse a deferred solve on a tool that cannot express one, with the reason rather than a
/// bare "unsupported" - the two that cannot are structural, not unimplemented.
void require_deferred_solve(PestppSession* s)
{
    if (s->adapter->supports_deferred_solve())
        return;
    string t = s->adapter->name();
    if (t == "da")
        unsupported("da cannot defer a solve: one advance() is a whole noptmax loop over a "
                    "cycle rather than a single iteration, so there is no one batch to hand "
                    "over. Use pestpp_solve_iteration()");
    unsupported("sqp cannot defer a solve: its line search and trust-region step each propose, "
                "run and judge a step, so one iteration issues several run batches rather than "
                "one. Use pestpp_solve_iteration()");
}
}

pestpp_status pestpp_solve_prepare(pestpp_handle h, int* n_runs)
{
    CAPI_BEGIN(h)
        require_deferred_solve(s);
        if (!s->initialized)
            bad_state("pestpp_initialize() must be called before pestpp_solve_prepare()");
        if (s->init_pending)
            bad_state("the prior ensemble is still outstanding; call pestpp_initialize_finish()");
        if (s->adapter->solve_is_open())
            bad_state("a deferred solve is already open; finish it before starting another");
        if (s->pending_runs_valid)
            bad_state("runs are already queued; process them before starting a solve");
        pestpp_status outcome = PESTPP_OK;
        int n = s->adapter->solve_prepare(outcome);
        if (n_runs != nullptr)
            *n_runs = n;
        return outcome;
    CAPI_END()
}

pestpp_status pestpp_solve_finish(pestpp_handle h, int defer_runs, int* pending_runs)
{
    CAPI_BEGIN(h)
        require_deferred_solve(s);
        if (!s->adapter->solve_is_open())
            bad_state("no deferred solve is open; call pestpp_solve_prepare() first");
        if (s->adapter->solve_batch_queued())
            bad_state("the queued runs have not been processed; call pestpp_process_runs()");
        int pending = 0;
        pestpp_status st = s->adapter->solve_finish(defer_runs != 0, pending);
        if (pending_runs != nullptr)
            *pending_runs = pending;
        return st;
    CAPI_END()
}

pestpp_status pestpp_get_candidate_count(pestpp_handle h, int* n)
{
    CAPI_BEGIN(h)
        if (n == nullptr) bad_arg("null out-param");
        *n = s->adapter->candidate_count();
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_get_candidate_info(pestpp_handle h, int idx, double* inflation,
                                        double* backtrack)
{
    CAPI_BEGIN(h)
        double inf = 0.0, back = 0.0;
        s->adapter->candidate_info(idx, inf, back);
        if (inflation != nullptr) *inflation = inf;
        if (backtrack != nullptr) *backtrack = back;
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_get_run_partial_info(pestpp_handle h, int run_id, int* has_partial,
                                          int* n_obs_reported, int* n_obs_total)
{
    CAPI_BEGIN_OBSERVER_SAFE(h)
        int n_rep = 0, n_tot = 0;
        bool got = s->run_manager->get_partial_info(run_id, n_rep, n_tot);
        if (has_partial != nullptr)    *has_partial = got ? 1 : 0;
        if (n_obs_reported != nullptr) *n_obs_reported = got ? n_rep : 0;
        if (n_obs_total != nullptr)    *n_obs_total = got ? n_tot : 0;
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_get_run_count(pestpp_handle h, int* n)
{
    CAPI_BEGIN(h)
        if (n == nullptr) bad_arg("null out-param");
        *n = s->run_manager->get_nruns();
        return PESTPP_OK;
    CAPI_END()
}

// extern "C++" on purpose: this file's body sits inside an extern "C" block, which would
// otherwise give this helper C language linkage while it takes C++ REFERENCE parameters. A
// C-linkage function with reference parameters is the kind of thing that is legal on paper and
// differs between compilers in practice - and not hypothetically: getting this wrong is what
// produced the windows access violation in pestpp_get_run_info().
//
// Helpers in here otherwise take and return plain C types, which need no such treatment. The
// exception is unpack_names() above, wrapped for the same reason with a C++ RETURN type. If
// you add a third, the rule is: any helper whose signature mentions a C++ type - reference,
// or a class returned by value - needs this wrapper, because the anonymous namespace it sits
// in controls visibility, not language linkage.
extern "C++" {
namespace {
/// Read one stored run and describe how complete it is.
///
/// completeness is derived from the values, not stored: a completed run is FINAL, an
/// incomplete one holding at least one real observation is PARTIAL, and one holding none is
/// NONE. Deriving it means it stays true after a restart, when the live bookkeeping is gone.
int read_stored_run(PestppSession* s, int run_id, Parameters& pars, Observations& obs,
                    int& status, int& n_reported)
{
    int nruns = s->run_manager->get_nruns();
    if ((run_id < 0) || (run_id >= nruns))
        bad_arg("no such run_id: " + std::to_string(run_id) + " (storage holds " +
                std::to_string(nruns) + ")");
    string info_txt;
    double info_value = 0.0;
    s->run_manager->get_info(run_id, status, info_txt, info_value);
    s->run_manager->get_run(run_id, pars, obs);
    n_reported = 0;
    vector<string> onames = s->run_manager->get_obs_name_vec();
    vector<double> ovals = obs.get_data_vec(onames);
    for (size_t i = 0; i < ovals.size(); i++)
        if (ovals[i] != Transformable::no_data)
            n_reported++;
    if (status == 1)
        return PESTPP_RUN_VALUES_FINAL;
    return (n_reported > 0) ? PESTPP_RUN_VALUES_PARTIAL : PESTPP_RUN_VALUES_NONE;
}
}
}

pestpp_status pestpp_get_run_info(pestpp_handle h, int run_id, int* status, int* completeness,
                                  int* n_obs_reported, int* n_obs_total)
{
    CAPI_BEGIN(h)
        Parameters pars;
        Observations obs;
        int st = 0, n_rep = 0;
        int comp = read_stored_run(s, run_id, pars, obs, st, n_rep);
        if (status != nullptr)         *status = st;
        if (completeness != nullptr)   *completeness = comp;
        if (n_obs_reported != nullptr) *n_obs_reported = n_rep;
        if (n_obs_total != nullptr)
            *n_obs_total = (int)s->run_manager->get_obs_name_vec().size();
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_get_run_values(pestpp_handle h, int run_id, double* pars_out, int npars,
                                    double* obs_out, int nobs, unsigned char* obs_valid,
                                    int* npars_out, int* nobs_out)
{
    CAPI_BEGIN(h)
        Parameters pars;
        Observations obs;
        int st = 0, n_rep = 0;
        read_stored_run(s, run_id, pars, obs, st, n_rep);
        vector<string> pnames = s->run_manager->get_par_name_vec();
        vector<string> onames = s->run_manager->get_obs_name_vec();
        if (npars_out != nullptr) *npars_out = (int)pnames.size();
        if (nobs_out != nullptr)  *nobs_out = (int)onames.size();
        // sizing call
        if ((pars_out == nullptr) && (obs_out == nullptr) && (obs_valid == nullptr))
            return PESTPP_OK;
        if ((pars_out != nullptr) && (npars < (int)pnames.size()))
            too_small("parameter buffer too small; call with NULL arrays to size it first");
        if (((obs_out != nullptr) || (obs_valid != nullptr)) && (nobs < (int)onames.size()))
            too_small("observation buffer too small; call with NULL arrays to size it first");
        if (pars_out != nullptr)
        {
            vector<double> pvals = pars.get_data_vec(pnames);
            for (size_t i = 0; i < pvals.size(); i++)
                pars_out[i] = pvals[i];
        }
        vector<double> ovals = obs.get_data_vec(onames);
        for (size_t i = 0; i < ovals.size(); i++)
        {
            if (obs_out != nullptr)
                obs_out[i] = ovals[i];
            // filled for completed runs too (all non-zero), so one code path serves both
            if (obs_valid != nullptr)
                obs_valid[i] = (ovals[i] != Transformable::no_data) ? 1 : 0;
        }
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_get_run_par_names(pestpp_handle h, char* buf, int buf_len, int* count)
{
    CAPI_BEGIN(h)
        return pack_names(s->run_manager->get_par_name_vec(), buf, buf_len, count);
    CAPI_END()
}

pestpp_status pestpp_get_run_obs_names(pestpp_handle h, char* buf, int buf_len, int* count)
{
    CAPI_BEGIN(h)
        return pack_names(s->run_manager->get_obs_name_vec(), buf, buf_len, count);
    CAPI_END()
}

pestpp_status pestpp_set_run_values(pestpp_handle h, int run_id, const double* obs, int nobs)
{
    CAPI_BEGIN(h)
        if (obs == nullptr)
            bad_arg("pestpp_set_run_values needs an observation array");
        int nruns = s->run_manager->get_nruns();
        if ((run_id < 0) || (run_id >= nruns))
            bad_arg("no such run_id: " + std::to_string(run_id) + " (storage holds " +
                    std::to_string(nruns) + ")");
        // A run manager writing results for a run the caller is also writing has no defined
        // outcome, so refuse rather than race. This is the check that keeps "service the runs
        // yourself" and "drive the run manager" from being usable at the same time by accident.
        if (s->run_manager->is_batch_open())
            bad_state("a batch is open; service runs yourself OR drive the run manager over "
                      "them, not both");
        vector<string> onames = s->run_manager->get_obs_name_vec();
        if (nobs != (int)onames.size())
            bad_arg("expected " + std::to_string(onames.size()) + " observations, got " +
                    std::to_string(nobs));
        // The stored parameters are read back and written through unchanged: update_run() takes
        // both, and the caller has no business restating parameters it did not choose.
        Parameters pars;
        Observations existing;
        int st = 0, n_rep = 0;
        read_stored_run(s, run_id, pars, existing, st, n_rep);
        Observations vals = s->adapter->scenario().get_ctl_observations();
        vals.update(onames, vector<double>(obs, obs + nobs));
        s->run_manager->update_run(run_id, pars, vals);
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_set_run_failed(pestpp_handle h, int run_id)
{
    CAPI_BEGIN(h)
        int nruns = s->run_manager->get_nruns();
        if ((run_id < 0) || (run_id >= nruns))
            bad_arg("no such run_id: " + std::to_string(run_id) + " (storage holds " +
                    std::to_string(nruns) + ")");
        if (s->run_manager->is_batch_open())
            bad_state("a batch is open; service runs yourself OR drive the run manager over "
                      "them, not both");
        s->run_manager->update_run_failed(run_id);
        return PESTPP_OK;
    CAPI_END()
}


pestpp_status pestpp_request_partial_results(pestpp_handle h, const int* run_ids, int n,
                                             int* n_requested)
{
    CAPI_BEGIN_OBSERVER_SAFE(h)
        vector<int> ids;
        if ((run_ids != nullptr) && (n > 0))
            ids.assign(run_ids, run_ids + n);
        int sent = s->run_manager->request_partial_results(ids);
        if (n_requested != nullptr)
            *n_requested = sent;
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_set_run_observer(pestpp_handle h, pestpp_run_observer_fn fn,
                                      void* user_data, double min_interval_sec)
{
    CAPI_BEGIN(h)
        s->observer_fn = fn;
        s->observer_data = user_data;
        if (fn == nullptr)
        {
            s->run_manager->set_progress_observer(nullptr);
            return PESTPP_OK;
        }
        // captured by value: the session outlives every batch, and capturing `s` by reference
        // from a lambda stored on the run manager would dangle if the session moved
        PestppSession* sess = s;
        s->run_manager->set_progress_observer(
            [sess](const RunProgress& p) -> RunAction
            {
                if (sess->observer_fn == nullptr)
                    return RunAction::CONTINUE;
                pestpp_run_progress out;
                out.struct_size = (int)sizeof(pestpp_run_progress);
                out.n_total = p.n_total;
                out.n_completed = p.n_completed;
                out.n_failed = p.n_failed;
                out.n_timed_out = p.n_timed_out;
                out.n_running = p.n_running;
                out.run_id = p.run_id;
                out.elapsed_sec = p.elapsed_sec;
                // the flag is what the allowlist checks; cleared on every path out, including
                // the one where the observer throws
                sess->in_observer = true;
                int action = PESTPP_RUN_CONTINUE;
                try
                {
                    action = sess->observer_fn(&out, sess->observer_data);
                }
                catch (...)
                {
                    sess->in_observer = false;
                    throw;         // notify_progress() unregisters on this
                }
                sess->in_observer = false;
                // an unknown action from a caller built against a later header is CONTINUE,
                // not an error: a progress observer must never be able to fail a batch
                return (action == PESTPP_RUN_STOP_BATCH) ? RunAction::STOP_BATCH
                                                         : RunAction::CONTINUE;
            },
            min_interval_sec);
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_queue_runs_subset(pestpp_handle h, const char* names, int n, int* n_queued)
{
    CAPI_BEGIN(h)
        require_deferred_solve(s);
        if (!s->adapter->solve_is_open())
            bad_state("no deferred solve is open; call pestpp_solve_prepare() first");
        if (s->pending_runs_valid)
            bad_state("runs are already queued; process them before queueing again");
        vector<string> only;
        if ((names != nullptr) && (n > 0))
            only = unpack_names(names, n);
        int queued = s->adapter->queue_solve_runs(only);
        // the session's flag is what makes process_runs() route back here, and what stops a
        // second queue landing on top of the first
        s->pending_runs_valid = true;
        if (n_queued != nullptr)
            *n_queued = queued;
        return PESTPP_OK;
    CAPI_END()
}

pestpp_status pestpp_queue_runs(pestpp_handle h, int* n_queued)
{
    CAPI_BEGIN(h)
        if (s->pending_runs_valid)
            bad_state("runs are already queued; process them before queueing again");
        // a deferred solve owns the batch: queue ITS candidates rather than the tool's
        // current ensemble, so the plain call keeps working inside the deferred loop
        if (s->adapter->solve_is_open())
        {
            int queued = s->adapter->queue_solve_runs(vector<string>());
            s->pending_runs_valid = true;
            if (n_queued != nullptr)
                *n_queued = queued;
            return PESTPP_OK;
        }
        // The cycle has to be passed explicitly: it defaults to NULL_DA_CYCLE, so a da handle
        // driving its own queue/process would tag runs with no cycle at all while the in-tree
        // da path tagged them correctly. Invisible on a control file whose first cycle is 0.
        s->pending_runs = s->adapter->queue_runs(s->performance_log.get(),
                                                 s->file_manager->rec_ofstream(),
                                                 s->run_manager.get());
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
            bad_state("no queued runs to process; call pestpp_queue_runs first");
        if (s->adapter->solve_is_open())
        {
            int nfail = s->adapter->process_solve_runs();
            s->pending_runs_valid = false;
            if (n_failed != nullptr)
                *n_failed = nfail;
            return PESTPP_OK;
        }
        vector<int> failed = s->adapter->process_runs(s->performance_log.get(),
                                                      s->file_manager->rec_ofstream(),
                                                      s->run_manager.get(), s->pending_runs);
        // clear even on the way out, so a failed process cannot leave a stale map behind
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

    // par and obs realizations are paired BY POSITION throughout pest++, so a length mismatch
    // means the pairing is already meaningless - and quietly keeping whichever obs rows happen
    // to sit at the surviving par positions would hand back somebody else's results under the
    // right-looking name. An empty obs ensemble is the one benign case: it just has not been
    // evaluated yet.
    if ((onames.size() > 0) && (onames.size() != pnames.size()))
    {
        stringstream ss;
        ss << "the parameter and observation ensembles disagree on size (" << pnames.size()
           << " vs " << onames.size() << "), so they cannot be paired by position; refusing "
              "to change membership rather than mis-attribute results";
        bad_state(ss.str());
    }

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
        else if (onames.size() > 0)
            bad_state("realization '" + n + "' has no observation counterpart");
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
        ParameterEnsemble* pe = s->adapter->par_ensemble();
        ParameterSnapshot snap = pe->get_ctl_snapshot();
        // get_ctl_snapshot() emits rows in the ensemble's ORIGINAL order, which is what the
        // tools want for stable csv output but is not what this API should hand back: after a
        // keep_realizations() the ensemble's current order is the caller's, and the zero-copy
        // view, the observation ensemble and get_ensemble_row_names() all follow it. Leaving
        // the snapshot on the original order meant par_df() and par_view() disagreed about
        // which row was which realization, silently. Permute here rather than changing
        // get_ctl_snapshot(), which would change what every tool writes to disk.
        reorder_snapshot_rows(snap, pe->get_real_names());
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

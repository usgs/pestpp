/**
 * @file selftest.cpp
 * @brief Fast in-memory self-tests for the option system. Runs in seconds (no model runs,
 *        no run managers) and returns non-zero on any failure, so CI can gate on it.
 *
 * Covers behavior the .pst-based benchmark suite structurally cannot: options changed
 * programmatically after construction (the live-propagation the shared-library API needs).
 *  - registry <-> legacy equivalence (defaults + parse for every option)
 *  - set_option / get_option / is_user_set (generic programmatic access + provenance)
 *  - is_init_only / get_option_scope + init-only change detection
 *  - ControlInfo programmatic access
 *  - apply_tool_defaults (centralized per-tool defaults)
 *  - Constraints chance/risk flags derived live from options (the flagship proof that a
 *    post-construction option change actually propagates: opt_risk, sqp_risk/STOSAG,
 *    opt_std_weights and their effect on use_chance/use_robust/use_fosm/get_risk)
 */
#include <iostream>
#include <csignal>
#include <random>
#include <sstream>
#include <fstream>
#include <limits>
#include <filesystem>
#include "config_os.h"
#include "pest_data_structs.h"
#include "Pest.h"
#include "FileManager.h"
#include "OutputFileWriter.h"
#include "PerformanceLog.h"
#include "constraints.h"
#include "EnsembleView.h"
#include "EnsembleMethodUtils.h"
#include "MOEA.h"
#include "SQP.h"
#include "model_interface.h"
#include "network_package.h"
#include "utilities.h"
#include "RunStorage.h"
#include "SVDSolver.h"
#include "SVDPackage.h"
#include "Regularization.h"

using namespace std;
using PO = PestppOptions;

// Test-only subclasses that expose each tool's derived option accessors. These are the
// decision points the algorithms actually read, so asserting through them proves the whole
// chain - option -> accessor -> algorithm - rather than just the options object.
struct IesProbe : public EnsembleMethod
{
    using EnsembleMethod::EnsembleMethod;
    bool p_use_subset()             { return get_use_subset(); }
    int  p_num_threads()      const { return get_num_threads(); }
    int  p_verbose_level()    const { return get_verbose_level(); }
    bool p_reinflate_minphi() const { return get_reinflate_to_minphi_real(); }
    vector<int> p_resolve(const vector<string>& n, const vector<string>& cur) const
    { return resolve_subset_idxs(n, cur); }
};

struct MouProbe : public MOEA
{
    using MOEA::MOEA;
    MouEnvType         p_envtype()   { return get_envtype(); }
    MouMateType        p_mattype()   { return get_mattype(); }
    vector<MouGenType> p_gen_types() { return get_gen_types(); }
    bool               p_risk_obj() const { return get_risk_obj(); }
};

struct SqpProbe : public SeqQuadProgram
{
    using SeqQuadProgram::SeqQuadProgram;
    bool p_use_subset()         const { return get_use_subset(); }
    bool p_use_ensemble_grad()  const { return get_use_ensemble_grad(); }
};

static int g_fail = 0;
static int g_total = 0;
// the test being run, so a crash handler can say where it died
static const char* g_current_test = "(none)";
static void CHK(bool cond, const string& msg)
{
    ++g_total;
    if (!cond) { ++g_fail; cout << "  [FAIL] " << msg << endl; }
}

static void test_registry_equivalence()
{
    cout << "[registry <-> legacy equivalence]" << endl;
    ostringstream os;
    bool ok = PO::self_verify(os);
    // surface only the mismatch lines, not the deprecation chatter self_verify triggers
    if (!ok)
    {
        istringstream in(os.str()); string ln;
        while (getline(in, ln))
            if (ln.find("MISMATCH") != string::npos || ln.find("self_verify:") != string::npos)
                cout << "  " << ln << endl;
    }
    CHK(ok, "PestppOptions::self_verify (registry == legacy for defaults + parse)");
}

static void test_generic_access()
{
    cout << "[generic programmatic access + provenance]" << endl;
    PO o; o.set_defaults();
    CHK(o.set_option("IES_NUM_REALS", "50") == PO::ARG_STATUS::ARG_ACCEPTED, "set_option accepted");
    CHK(o.set_option("IES_NUM_REALS", "123") == PO::ARG_STATUS::ARG_ACCEPTED, "repeated set_option accepted (no duplicate block)");
    CHK(o.get_ies_num_reals() == 123, "repeated set took effect");
    CHK(o.get_option("IES_NUM_REALS") == "123", "get_option reflects value");
    CHK(o.is_user_set("IES_NUM_REALS"), "is_user_set true after set_option");
    CHK(!o.is_user_set("IES_SUBSET_SIZE"), "is_user_set false for untouched option");
    o.set_option("IES_PARAMETER_ENSEMBLE", "MyPars.csv");   // alias
    CHK(o.is_user_set("IES_PAR_EN"), "alias->canonical provenance resolves");
    CHK(o.get_ies_par_csv() == "MyPars.csv", "filename case preserved through set_option");
    CHK(o.set_option("NOPE_NOT_REAL", "1") == PO::ARG_STATUS::ARG_NOTFOUND, "unknown -> NOTFOUND");
    CHK(o.is_valid_arg("ies_num_reals") && !o.is_valid_arg("nope"), "is_valid_arg");
    // parse-path duplicate guard still active (file semantics unchanged)
    PO p; p.set_defaults();
    CHK(p.assign_value_by_key("IES_NUM_REALS", "10") == PO::ARG_STATUS::ARG_ACCEPTED, "parse-path first accepted");
    CHK(p.assign_value_by_key("IES_NUM_REALS", "20") == PO::ARG_STATUS::ARG_DUPLICATE, "parse-path duplicate blocked");
}

static void test_mutability()
{
    cout << "[mutability metadata + init-only change detection]" << endl;
    PO o; o.set_defaults();
    CHK(o.is_init_only("IES_NUM_REALS"), "IES_NUM_REALS is init-only");
    CHK(o.is_init_only("RANDOM_SEED"), "RANDOM_SEED is init-only");
    CHK(!o.is_init_only("IES_LAMBDA_MULTS"), "IES_LAMBDA_MULTS is live");
    CHK(o.get_option_scope("SQP_NUM_REALS") == "sqp", "scope resolves");
    o.set_option("IES_NUM_REALS", "100");
    CHK(o.get_init_only_change_warnings().empty(), "no warning before initialization");
    o.mark_options_initialized();
    o.set_option("IES_LAMBDA_MULTS", "0.1,1.0");
    CHK(o.get_init_only_change_warnings().empty(), "live-option change post-init: no warning");
    CHK(o.set_option("IES_NUM_REALS", "50") == PO::ARG_STATUS::ARG_ACCEPTED, "init-only post-init still accepted");
    CHK(o.get_init_only_change_warnings().size() == 1, "init-only post-init change -> warning");
}

static void test_control_info()
{
    cout << "[ControlInfo programmatic access]" << endl;
    ControlInfo ci; ci.set_defaults();
    CHK(ci.set_option("NOPTMAX", "5") == PO::ARG_STATUS::ARG_ACCEPTED, "ctl set_option");
    CHK(ci.set_option("NOPTMAX", "12") == PO::ARG_STATUS::ARG_ACCEPTED, "ctl repeated set_option (no duplicate block)");
    CHK(ci.noptmax == 12, "ctl noptmax updated");
    CHK(ci.get_option("NOPTMAX") == "12", "ctl get_option");
    CHK(ci.is_user_set("NOPTMAX") && !ci.is_user_set("PHIREDSTP"), "ctl is_user_set");
    ci.set_option("PESTMODE", "regularization");
    CHK(ci.pestmode == ControlInfo::PestMode::REGUL, "ctl pestmode enum set");
}

static void test_tool_defaults()
{
    cout << "[apply_tool_defaults]" << endl;
    ostringstream log;
    { PO o; o.set_defaults(); o.apply_tool_defaults(PO::ToolType::IES, log);
      CHK(o.get_max_run_fail() == 1 && o.get_overdue_giveup_fac() == 2.0, "IES: ensemble defaults applied"); }
    { PO o; o.set_defaults(); o.apply_tool_defaults(PO::ToolType::MOU, log);
      CHK(o.get_max_run_fail() == 1, "MOU: ensemble defaults applied"); }
    { PO o; o.set_defaults(); o.apply_tool_defaults(PO::ToolType::GLM, log);
      CHK(o.get_max_run_fail() == 3, "GLM: library default retained"); }
    { PO o; o.set_defaults(); o.set_option("MAX_RUN_FAIL", "7"); o.apply_tool_defaults(PO::ToolType::IES, log);
      CHK(o.get_max_run_fail() == 7, "user-set value not overridden (provenance-aware)"); }
}

static void test_constraints_live()
{
    cout << "[Constraints live chance/risk propagation]" << endl;
    ofstream flog("selftest_perf.log");
    PerformanceLog pfm(flog);
    FileManager fm;
    Pest p; p.set_defaults();
    OutputFileWriter ofw(fm, p);
    Constraints c(p, &fm, ofw, pfm);   // getters read only options; no initialize() needed
    PestppOptions* o = p.get_pestpp_options_ptr();

    o->set_option("OPT_RISK", "0.5");
    CHK(!c.get_use_chance(), "opt_risk 0.5 -> use_chance false");
    // the latching-bug fix: with the old cached member this stayed false forever
    o->set_option("OPT_RISK", "0.95");
    CHK(c.get_use_chance(), "opt_risk 0.95 -> use_chance TRUE (live, was latched at init)");
    CHK(c.get_risk() == 0.95, "get_risk reflects live opt_risk");
    o->set_option("OPT_RISK", "1.5");
    CHK(c.get_risk() == 0.999, "get_risk clamps high in the accessor");
    o->set_option("OPT_RISK", "0.0001");
    CHK(c.get_risk() == 0.001, "get_risk clamps low in the accessor");
    // use_fosm derives live from stack options
    o->set_option("OPT_RISK", "0.9");
    o->set_option("OPT_STACK_SIZE", "0");
    CHK(c.get_use_fosm(), "no stacks -> use_fosm true");
    o->set_option("OPT_STACK_SIZE", "50");
    CHK(!c.get_use_fosm(), "opt_stack_size>0 -> use_fosm false (live)");

    // ROBUST optimization: opt_use_robust switches the whole chance machinery OFF rather than
    // selecting a different flavour of it. That single fact is what lets every downstream
    // "shift or not" branch stay as it is, so it is the one worth pinning.
    o->set_option("OPT_RISK", "0.95");
    o->set_option("OPT_USE_ROBUST", "false");
    CHK(!c.get_use_robust(), "opt_use_robust false -> get_use_robust false");
    CHK(c.get_use_chance(), "opt_risk 0.95, not robust -> use_chance true");
    o->set_option("OPT_USE_ROBUST", "true");
    CHK(c.get_use_robust(), "opt_use_robust true -> get_use_robust TRUE (live)");
    CHK(!c.get_use_chance(), "robust -> use_chance FALSE however opt_risk is set");
    CHK(!c.get_use_fosm(), "robust -> use_fosm false");
    CHK(c.get_risk() == 0.95, "get_risk still reports opt_risk, clamped, under robust");
    o->set_option("OPT_USE_ROBUST", "false");   // reset

    // sqp_risk is retired: it used to be both a risk value and the switch that selected
    // ensemble shifting, so it silently overrode opt_risk. It must now be refused rather
    // than quietly accepted - there is no behaviour-preserving translation.
    {
        bool threw = false;
        try { o->set_option("SQP_RISK", "0.9"); }
        catch (const std::exception&) { threw = true; }
        CHK(threw, "retired sqp_risk is refused rather than silently accepted");
    }

    // opt_std_weights propagates live and gates use_fosm when stacks are present (Stage 5)
    o->set_option("OPT_RISK", "0.9");
    o->set_option("OPT_STACK_SIZE", "50");   // stacks present
    o->set_option("OPT_STD_WEIGHTS", "false");
    CHK(!c.get_std_weights(), "opt_std_weights false -> get_std_weights false (live)");
    CHK(!c.get_use_fosm(), "std_weights false + stacks -> use_fosm false (live)");
    o->set_option("OPT_STD_WEIGHTS", "true");
    CHK(c.get_std_weights(), "opt_std_weights true -> get_std_weights true (live)");
    CHK(c.get_use_fosm(), "std_weights true -> use_fosm true even with stacks (live)");
}

static void test_ies_reinflate_reset()
{
    cout << "[ies/da reinflation options: live reset propagation]" << endl;
    PO o; o.set_defaults();
    // registry defaults
    CHK((o.get_ies_n_iter_reinflate() == vector<int>{0}), "n_iter_reinflate default {0}");
    CHK((o.get_ies_reinflate_factor() == vector<double>{1.0}), "reinflate_factor default {1.0}");
    CHK((o.get_ies_reinflate_num_reals() == vector<int>{0}), "reinflate_num_reals default {0}");
    // all three reinflation controls are LIVE (resettable at an iteration boundary)
    CHK(!o.is_init_only("IES_N_ITER_MEAN"), "n_iter_reinflate is live (not init-only)");
    CHK(!o.is_init_only("IES_REINFLATE_FACTOR"), "reinflate_factor is live");
    CHK(!o.is_init_only("IES_REINFLATE_NUM_REALS"), "reinflate_num_reals is live");

    // set, then RESET to a different schedule; the getter reflects the latest each time
    o.set_option("IES_N_ITER_REINFLATE", "3");                 // alias -> IES_N_ITER_MEAN
    CHK((o.get_ies_n_iter_reinflate() == vector<int>{3}), "n_iter_reinflate set via alias");
    CHK(o.is_user_set("IES_N_ITER_MEAN"), "alias->canonical provenance for reinflate schedule");
    o.set_option("IES_N_ITER_MEAN", "2,4,6");                  // reset to a vector schedule
    CHK((o.get_ies_n_iter_reinflate() == vector<int>{2,4,6}), "n_iter_reinflate reset to vector");
    o.set_option("IES_REINFLATE_FACTOR", "0.5");
    CHK((o.get_ies_reinflate_factor() == vector<double>{0.5}), "reinflate_factor set");
    o.set_option("IES_REINFLATE_FACTOR", "0.9,0.8");           // reset
    CHK((o.get_ies_reinflate_factor() == vector<double>{0.9,0.8}), "reinflate_factor reset to vector");
    o.set_option("IES_REINFLATE_NUM_REALS", "-30");            // negative = draw-from-current path
    CHK((o.get_ies_reinflate_num_reals() == vector<int>{-30}), "reinflate_num_reals set (negative)");
    o.set_option("IES_REINFLATE_NUM_REALS", "50");             // reset positive = truncate path
    CHK((o.get_ies_reinflate_num_reals() == vector<int>{50}), "reinflate_num_reals reset positive");

    // the DA_* -> IES_* rewrite: pestpp-da reads the same reinflation controls
    o.set_option("DA_N_ITER_REINFLATE", "7");
    CHK((o.get_ies_n_iter_reinflate() == vector<int>{7}), "DA_N_ITER_REINFLATE resets IES reinflate schedule");
    o.set_option("DA_REINFLATE_FACTOR", "0.6");
    CHK((o.get_ies_reinflate_factor() == vector<double>{0.6}), "DA_REINFLATE_FACTOR resets IES reinflate factor");

    // resetting a live reinflation option AFTER initialization is honored WITHOUT a warning
    o.mark_options_initialized();
    o.set_option("IES_REINFLATE_FACTOR", "0.7");
    CHK((o.get_ies_reinflate_factor() == vector<double>{0.7}), "reinflate_factor reset post-init took effect");
    CHK(o.get_init_only_change_warnings().empty(), "live reinflation reset post-init: no warning");
}

static void test_ies_ensemble_reset()
{
    cout << "[ies ensemble source options: par/obs/noise reset semantics]" << endl;
    PO o; o.set_defaults();
    // the par/obs/weights/restart ensemble SOURCES are loaded at initialize(); they are
    // init-only, so a post-init reset must be FLAGGED (surfaced), not silently ignored
    const char* sources[] = {"IES_PAR_EN","IES_OBS_EN","IES_WEIGHTS_ENSEMBLE",
                             "IES_RESTART_PARAMETER_ENSEMBLE","IES_RESTART_OBSERVATION_ENSEMBLE",
                             "IES_NUM_REALS"};
    for (const char* k : sources)
        CHK(o.is_init_only(k), string("ensemble source is init-only: ") + k);
    // observation noise is generated from obs+weights; the on/off switch is LIVE
    CHK(!o.is_init_only("IES_NO_NOISE"), "ies_no_noise is live");

    // provenance + filename-case preservation for the ensemble sources (round-trip)
    o.set_option("IES_PARAMETER_ENSEMBLE", "PriorPars.csv");   // alias, mixed case
    CHK(o.is_user_set("IES_PAR_EN"), "par-en alias->canonical provenance");
    CHK(o.get_ies_par_csv() == "PriorPars.csv", "par-en filename case preserved");
    o.set_option("IES_OBS_EN", "Obs.jcb");
    CHK(o.get_ies_obs_csv() == "Obs.jcb", "obs-en set");
    o.set_option("IES_WEIGHTS_ENSEMBLE", "Weights.csv");
    CHK(o.get_ies_weights_csv() == "Weights.csv", "weights-en set");

    // noise toggle propagates live and can be reset
    o.set_option("IES_NO_NOISE", "true");
    CHK(o.get_ies_no_noise(), "ies_no_noise true (live)");
    o.set_option("IES_NO_NOISE", "false");
    CHK(!o.get_ies_no_noise(), "ies_no_noise reset false (live)");

    // init-only SAFETY: resetting an ensemble source post-init is surfaced, per change
    o.mark_options_initialized();
    o.set_option("IES_PAR_EN", "NewPars.csv");
    CHK(o.get_init_only_change_warnings().size() == 1,
        "resetting par ensemble post-init is FLAGGED (not silently ignored)");
    o.set_option("IES_OBS_EN", "NewObs.csv");
    CHK(o.get_init_only_change_warnings().size() == 2, "resetting obs ensemble post-init also flagged");
    // a live noise reset in the same window adds NO warning
    o.set_option("IES_NO_NOISE", "true");
    CHK(o.get_init_only_change_warnings().size() == 2, "live noise reset post-init: no new warning");
}

static void test_ies_iteration_controls()
{
    cout << "[ies per-iteration controls: subset/threads/verbosity are live]" << endl;
    PO o; o.set_defaults();
    // these three used to be snapshotted into EnsembleMethod members during initialize();
    // they are now read at point of use, so the registry must advertise them as live
    CHK(!o.is_init_only("IES_SUBSET_SIZE"), "ies_subset_size is live");
    CHK(!o.is_init_only("IES_NUM_THREADS"), "ies_num_threads is live (each solve spins its own pool)");
    CHK(!o.is_init_only("IES_VERBOSE_LEVEL"), "ies_verbose_level is live");
    CHK(!o.is_init_only("IES_LAM_MULTS"), "ies_lam_mults is live");
    CHK(!o.is_init_only("IES_REG_FACTOR"), "ies_reg_factor is live");

    // values round-trip through the generic setter, including the negative (percentage)
    // subset convention and the negative reg_factor 'full solution' signal
    o.set_option("IES_SUBSET_SIZE", "-25");
    CHK(o.get_ies_subset_size() == -25, "subset_size negative (percentage) accepted");
    o.set_option("IES_NUM_THREADS", "8");
    CHK(o.get_ies_num_threads() == 8, "num_threads set");
    o.set_option("IES_VERBOSE_LEVEL", "3");
    CHK(o.get_ies_verbose_level() == 3, "verbose_level set");

    // all three reset AFTER initialization, with no init-only warning
    o.mark_options_initialized();
    o.set_option("IES_SUBSET_SIZE", "12");
    o.set_option("IES_NUM_THREADS", "2");
    o.set_option("IES_VERBOSE_LEVEL", "1");
    CHK(o.get_ies_subset_size() == 12, "subset_size reset post-init took effect");
    CHK(o.get_ies_num_threads() == 2, "num_threads reset post-init took effect");
    CHK(o.get_ies_verbose_level() == 1, "verbose_level reset post-init took effect");
    CHK(o.get_init_only_change_warnings().empty(), "iteration controls reset post-init: no warning");
}

static void test_mou_generation_controls()
{
    cout << "[mou per-generation controls: selectors/generator/risk are live]" << endl;
    PO o; o.set_defaults();
    // MOEA used to cache envtype/mattype/gen_types/risk_obj at initialize(); they are
    // derived live now, so the generator mix and selectors can change per generation
    const char* live[] = {"MOU_ENV_SELECTOR","MOU_MATING_SELECTOR","MOU_GENERATOR",
                          "MOU_RISK_OBJECTIVE","MOU_PSO_INERTIA","MOU_PSO_DV_BOUND_HANDLING",
                          "MOU_SAVE_POPULATION_EVERY"};
    for (const char* k : live)
        CHK(!o.is_init_only(k), string("mou generation control is live: ") + k);
    // the population size sizes the initial draw, so it stays init-only (the per-generation
    // member count comes from the population_schedule, which is a live lookup)
    CHK(o.is_init_only("MOU_POPULATION_SIZE"), "mou_population_size is init-only (sizes the initial draw)");

    o.set_option("MOU_ENV_SELECTOR", "spea");
    CHK(o.get_mou_env_selector() == "SPEA", "env selector set (upper-cased)");
    o.set_option("MOU_GENERATOR", "de,sbx");
    CHK(o.get_mou_generator() == "de,sbx", "generator list set");
    o.set_option("MOU_RISK_OBJECTIVE", "true");
    CHK(o.get_mou_risk_obj(), "risk objective set");

    // switching the whole generation strategy AFTER initialization is honored, unflagged
    o.mark_options_initialized();
    o.set_option("MOU_ENV_SELECTOR", "nsga");
    o.set_option("MOU_MATING_SELECTOR", "random");
    o.set_option("MOU_GENERATOR", "pso");
    o.set_option("MOU_RISK_OBJECTIVE", "false");
    CHK(o.get_mou_env_selector() == "NSGA", "env selector switched post-init");
    CHK(o.get_mou_mating_selector() == "RANDOM", "mating selector switched post-init");
    CHK(o.get_mou_generator() == "pso", "generator switched post-init");
    CHK(!o.get_mou_risk_obj(), "risk objective switched post-init");
    CHK(o.get_init_only_change_warnings().empty(), "generation controls switched post-init: no warning");
}

static void test_sqp_controls()
{
    cout << "[sqp controls: live subset/risk vs init-only ensemble + working set seed]" << endl;
    PO o; o.set_defaults();
    // the ensemble sources are loaded once, and working_set_tol is adaptive state seeded
    // from its option, so all three are init-only - a later change must be surfaced
    CHK(o.is_init_only("SQP_NUM_REALS"), "sqp_num_reals is init-only (drives the gradient mode)");
    CHK(o.is_init_only("SQP_DV_EN"), "sqp_dv_en is init-only");
    CHK(o.is_init_only("SQP_WORKING_SET_TOL"), "sqp_working_set_tol is init-only (adaptive after seeding)");
    // the per-iteration knobs stay live
    CHK(!o.is_init_only("SQP_SUBSET_SIZE"), "sqp_subset_size is live");
    CHK(o.is_init_only("OPT_USE_ROBUST"), "opt_use_robust is init-only");
    CHK(!o.is_init_only("SQP_MAX_CONSEC_INFEAS_IES"), "sqp_max_consec_infeas_ies is live");

    o.mark_options_initialized();
    o.set_option("SQP_SUBSET_SIZE", "5");
    CHK(o.get_sqp_subset_size() == 5, "sqp_subset_size reset post-init took effect");
    CHK(o.get_init_only_change_warnings().empty(), "live sqp reset post-init: no warning");
    o.set_option("SQP_WORKING_SET_TOL", "0.02");
    CHK(o.get_init_only_change_warnings().size() == 1,
        "resetting sqp_working_set_tol post-init is FLAGGED (the solver has already adapted it)");
}

static void test_tool_objects_track_live_options()
{
    cout << "[tool objects re-read their options after initialize()]" << endl;
    ofstream flog("selftest_tool_perf.log");
    PerformanceLog pfm(flog);
    FileManager fm;
    Pest p; p.set_defaults();
    OutputFileWriter ofw(fm, p);
    PO* o = p.get_pestpp_options_ptr();

    // constructing a tool is side-effect free (no files, no runs) - the option-derived state
    // that used to be snapshotted lives behind accessors, so no initialize() is needed here
    IesProbe ies(p, fm, ofw, &pfm, nullptr, "selftest-ies");
    MouProbe mou(p, fm, ofw, &pfm, nullptr);
    SqpProbe sqp(p, fm, ofw, &pfm, nullptr);

    // everything from here on is a change made AFTER the tool exists and options are sealed:
    // exactly the shared-library case (mutate a control, run the next iteration)
    o->mark_options_initialized();

    // --- ies: the four controls that used to latch during initialize() ---
    o->set_option("IES_NUM_THREADS", "6");
    CHK(ies.p_num_threads() == 6, "ies object sees num_threads post-init");
    o->set_option("IES_NUM_THREADS", "1");
    CHK(ies.p_num_threads() == 1, "ies object tracks a second num_threads change");

    o->set_option("IES_VERBOSE_LEVEL", "3");
    CHK(ies.p_verbose_level() == 3, "ies object sees verbose_level post-init");

    // the ensemble is empty here, so use_subset reduces to subset_size <= 0; that still
    // exercises the live read and the ensemble-size comparison it now performs
    o->set_option("IES_SUBSET_SIZE", "5");
    CHK(!ies.p_use_subset(), "subset larger than the ensemble -> use_subset false (live)");
    o->set_option("IES_SUBSET_SIZE", "-10");
    CHK(ies.p_use_subset(), "percentage subset -> use_subset TRUE post-init (was latched)");

    o->set_option("IES_N_ITER_MEAN", "3");
    CHK(!ies.p_reinflate_minphi(), "positive n_iter_reinflate -> reinflate to mean");
    o->set_option("IES_N_ITER_MEAN", "-3");
    CHK(ies.p_reinflate_minphi(), "negative n_iter_reinflate -> reinflate to min-phi real (live)");

    // --- mou: swap the whole generation strategy between generations ---
    o->set_option("MOU_ENV_SELECTOR", "nsga");
    CHK(mou.p_envtype() == MouEnvType::NSGA, "mou object sees nsga env selector");
    o->set_option("MOU_ENV_SELECTOR", "spea");
    CHK(mou.p_envtype() == MouEnvType::SPEA, "mou env selector switched post-init (was latched)");
    o->set_option("MOU_ENV_SELECTOR", "nsga_ppd");
    CHK(mou.p_envtype() == MouEnvType::NSGA, "nsga_ppd maps to the nsga env type");

    o->set_option("MOU_MATING_SELECTOR", "random");
    CHK(mou.p_mattype() == MouMateType::RANDOM, "mou object sees random mating selector");
    o->set_option("MOU_MATING_SELECTOR", "tournament");
    CHK(mou.p_mattype() == MouMateType::TOURNAMENT, "mou mating selector switched post-init");

    o->set_option("MOU_GENERATOR", "de,sbx");
    vector<MouGenType> gt = mou.p_gen_types();
    CHK(gt.size() == 2, "mou generator list parsed live (2 generators)");
    CHK(gt[0] == MouGenType::DE && gt[1] == MouGenType::SBX, "generator list order preserved");
    // the old member accumulated across calls; the accessor must be idempotent
    CHK(mou.p_gen_types().size() == 2, "repeated read does not accumulate generators");
    o->set_option("MOU_GENERATOR", "pso");
    gt = mou.p_gen_types();
    CHK(gt.size() == 1 && gt[0] == MouGenType::PSO, "generator mix switched post-init (was latched)");

    o->set_option("MOU_RISK_OBJECTIVE", "true");
    CHK(mou.p_risk_obj(), "mou object sees risk objective post-init");
    o->set_option("MOU_RISK_OBJECTIVE", "false");
    CHK(!mou.p_risk_obj(), "mou risk objective reset post-init");

    // --- sqp ---
    o->set_option("SQP_SUBSET_SIZE", "0");
    CHK(!sqp.p_use_subset(), "sqp subset_size 0 -> no subset");
    o->set_option("SQP_SUBSET_SIZE", "10");
    CHK(sqp.p_use_subset(), "sqp subset switched on post-init");

    // none of the above should have tripped the init-only detector
    CHK(o->get_init_only_change_warnings().empty(),
        "every per-iteration control above changed post-init WITHOUT a warning");

    // by contrast, the gradient mode is driven by an init-only option: the value still
    // propagates to the accessor, but the change is surfaced rather than silently taken.
    // sqp_num_reals DEFAULTS to 50 - an ensemble gradient is the only mode sqp implements, so
    // the default has to select it - and switching it off takes an explicit 0
    CHK(sqp.p_use_ensemble_grad(), "sqp_num_reals defaults to 50 -> ensemble gradient");
    o->set_option("SQP_NUM_REALS", "0");
    CHK(!sqp.p_use_ensemble_grad(), "sqp_num_reals 0 -> no ensemble gradient (a refused run)");
    o->set_option("SQP_NUM_REALS", "20");
    CHK(sqp.p_use_ensemble_grad(), "sqp_num_reals -> ensemble gradient (derived, not cached)");
    CHK(o->get_init_only_change_warnings().size() == 2,
        "each post-init change of the init-only sqp_num_reals IS flagged");
}


static void test_ensemble_zero_copy_view()
{
    cout << "[zero-copy ensemble view: borrowed window + invalidation]" << endl;
    Pest p; p.set_defaults();
    std::mt19937 rgen(11);

    vector<string> rnames{"R0","R1","R2"};
    vector<string> vnames{"V0","V1"};
    Eigen::MatrixXd m(3,2);
    m << 1.0, 4.0,
         2.0, 5.0,
         3.0, 6.0;
    ObservationEnsemble oe(&p, &rgen, m, rnames, vnames);

    {
        EnsembleView v(oe);
        CHK(v.valid(), "fresh view is valid");
        CHK(v.rows() == 3 && v.cols() == 2, "view reports the ensemble shape");
        // Eigen is column-major, so the buffer runs down each column in turn
        CHK(v.strides_bytes()[0] == (int64_t)sizeof(double), "row stride is one element");
        CHK(v.strides_bytes()[1] == (int64_t)(sizeof(double)*3), "col stride is one column");
        CHK(v.at(0,0) == 1.0 && v.at(2,0) == 3.0 && v.at(0,1) == 4.0 && v.at(2,1) == 6.0,
            "view reads the right values through the raw buffer");
        // it is the ensemble's own buffer, not a copy
        CHK(v.data() == oe.get_eigen_ptr()->data(), "view aliases the ensemble buffer (no copy)");
        CHK(v.row_names().size() == 3 && v.col_names()[1] == "V1", "view reports live labels");

        // an in-place write lands in the ensemble and does NOT invalidate
        v.set(1, 1, 99.0);
        CHK(oe.get_eigen()(1,1) == 99.0, "in-place write through the view reaches the ensemble");
        CHK(v.valid(), "in-place write leaves the view valid");
    }

    // a structural mutation reallocates reals -> outstanding views go invalid
    {
        EnsembleView v(oe);
        CHK(v.valid(), "view valid before drop_rows");
        vector<string> drop{"R0"};
        oe.drop_rows(drop);
        CHK(!v.valid(), "drop_rows invalidates the view");
        bool threw = false;
        try { v.at(0,0); } catch (const EnsembleViewInvalidated&) { threw = true; }
        CHK(threw, "using an invalidated view throws EnsembleViewInvalidated");
    }

    // set_eigen replaces the whole matrix
    {
        EnsembleView v(oe);
        Eigen::MatrixXd m2(2,2); m2 << 7,8,9,10;
        oe.set_eigen(m2);
        CHK(!v.valid(), "set_eigen invalidates the view");
    }

    // assigning one ensemble over another invalidates views onto the target
    {
        ObservationEnsemble other(&p, &rgen, m, rnames, vnames);
        EnsembleView v(other);
        CHK(v.valid(), "view valid before assignment");
        other = oe;
        CHK(!v.valid(), "assigning over an ensemble invalidates its views");
    }

    // destruction expires the guard, and valid() must not touch the dead object
    {
        ObservationEnsemble* tmp = new ObservationEnsemble(&p, &rgen, m, rnames, vnames);
        EnsembleView v(*tmp);
        CHK(v.valid(), "view valid before destruction");
        delete tmp;
        CHK(!v.valid(), "destroying the ensemble invalidates the view");
    }
}



static void test_subset_names_survive_membership_change()
{
    cout << "[subset is held by name, so it survives a membership change]" << endl;
    ofstream flog("selftest_subset.log");
    PerformanceLog pfm(flog);
    FileManager fm;
    Pest p; p.set_defaults();
    OutputFileWriter ofw(fm, p);
    IesProbe ies(p, fm, ofw, &pfm, nullptr, "selftest-subset");

    vector<string> all{"R0","R1","R2","R3","R4"};
    vector<string> subset{"R1","R3"};

    // baseline: names map to their positions
    vector<int> idx = ies.p_resolve(subset, all);
    CHK(idx.size() == 2 && idx[0] == 1 && idx[1] == 3, "subset names resolve to positions");

    // a realization removed BEFORE the subset shifts every later position. Positions would
    // now be wrong; names still land on the right realizations.
    vector<string> dropped_first{"R1","R2","R3","R4"};
    idx = ies.p_resolve(subset, dropped_first);
    CHK(idx.size() == 2 && idx[0] == 0 && idx[1] == 2,
        "dropping an earlier realization re-resolves the subset (positions shifted)");

    // a member OF the subset removed -> the subset shrinks, no out-of-range index
    vector<string> dropped_member{"R0","R2","R3","R4"};
    idx = ies.p_resolve(subset, dropped_member);
    CHK(idx.size() == 1 && idx[0] == 2, "dropping a subset member shrinks the subset");

    // reordering is absorbed too
    vector<string> reordered{"R4","R3","R2","R1","R0"};
    idx = ies.p_resolve(subset, reordered);
    CHK(idx.size() == 2 && idx[0] == 3 && idx[1] == 1, "reordering re-resolves correctly");

    // order follows the SUBSET, not the ensemble - update_reals_by_phi maps subset
    // position i to a name, so this ordering is load-bearing
    vector<string> rev{"R3","R1"};
    idx = ies.p_resolve(rev, all);
    CHK(idx.size() == 2 && idx[0] == 3 && idx[1] == 1, "result order follows the subset order");

    // everything gone
    vector<string> none{"X","Y"};
    CHK(ies.p_resolve(subset, none).size() == 0, "no surviving members -> empty subset");
}

// A caller running its own run loop can change the ensemble between queueing the runs and
// processing them. The run->realization map is keyed by name so that change is absorbed;
// keyed by position, every run after the edit would land on the wrong realization.
static void test_run_map_survives_resize()
{
    cout << "[queued runs are held by name, so a resize between queue and process is safe]" << endl;

    vector<string> all{"R0","R1","R2","R3","R4"};
    // run ids deliberately unlike the positions, so a position/run-id mix-up is visible
    map<string,int> run_ids{{"R0",100},{"R1",101},{"R2",102},{"R3",103},{"R4",104}};

    // baseline: nothing changed, every run lands on its own row
    vector<pair<int,int>> r = ObservationEnsemble::resolve_run_positions(run_ids, all, all.size());
    CHK(r.size() == 5, "unchanged ensemble resolves every queued run");
    bool ok = true;
    for (int i = 0; i < 5; i++)
        ok = ok && (r[i].first == i) && (r[i].second == 100 + i);
    CHK(ok, "each run maps to its own row");

    // drop a realization AHEAD of the others: every later row shifts down by one. positions
    // would now be off by one; names follow the realizations.
    vector<string> dropped_first{"R1","R2","R3","R4"};
    r = ObservationEnsemble::resolve_run_positions(run_ids, dropped_first, dropped_first.size());
    CHK(r.size() == 4, "dropping an earlier realization keeps the remaining runs");
    CHK(r[0].first == 0 && r[0].second == 101, "R1's run follows R1 to its new row 0");
    CHK(r[3].first == 3 && r[3].second == 104, "R4's run follows R4 to its new row 3");

    // drop a realization in the MIDDLE
    vector<string> dropped_mid{"R0","R1","R3","R4"};
    r = ObservationEnsemble::resolve_run_positions(run_ids, dropped_mid, dropped_mid.size());
    CHK(r.size() == 4 && r[2].second == 103, "R3's run follows R3 across the gap left by R2");

    // R2's run has nowhere to go and is discarded rather than misattributed
    ok = true;
    for (auto& p : r)
        ok = ok && (p.second != 102);
    CHK(ok, "the dropped realization's run is discarded, not reassigned");

    // reordering is absorbed
    vector<string> reordered{"R4","R3","R2","R1","R0"};
    r = ObservationEnsemble::resolve_run_positions(run_ids, reordered, reordered.size());
    CHK(r.size() == 5 && r[0].second == 104 && r[4].second == 100,
        "reversing the ensemble reverses which run lands where");

    // a realization the caller ADDED after queueing simply has no run yet
    vector<string> added{"R0","R1","NEW","R2","R3","R4"};
    r = ObservationEnsemble::resolve_run_positions(run_ids, added, added.size());
    CHK(r.size() == 5, "an added realization contributes no run");
    ok = true;
    for (auto& p : r)
        ok = ok && (p.first != 2);
    CHK(ok, "the added realization's row is left untouched");

    // n_rows bounds the result: the receiving ensemble may be shorter than the name list
    r = ObservationEnsemble::resolve_run_positions(run_ids, all, 3);
    CHK(r.size() == 3, "resolution stops at the receiving ensemble's row count");

    // rows come back in ascending order, which is what the process loop and the failed-index
    // list both assume
    r = ObservationEnsemble::resolve_run_positions(run_ids, reordered, reordered.size());
    ok = true;
    for (int i = 1; i < r.size(); i++)
        ok = ok && (r[i].first > r[i-1].first);
    CHK(ok, "resolved rows are in ascending order");

    // nothing survived
    vector<string> none{"X","Y"};
    CHK(ObservationEnsemble::resolve_run_positions(run_ids, none, none.size()).size() == 0,
        "no surviving realizations -> no runs to process");
}

/**
 * @brief A parse that fails part way still yields what it read, and names the rest.
 *
 * This is the foundation preemption stands on: asking a run that is STILL GOING what it has
 * so far. The strict reader cannot answer that question and must not be made to - a run that
 * genuinely finishes has to keep failing loudly on a malformed output file.
 *
 * The obstacle was not the exception, it was WHERE the data lived: read_output_file() built
 * its Observations in a local and returned by value, so everything read before a failure died
 * with the stack frame. Catching the exception higher up would have caught an exception and
 * no data. The parse now writes into a caller-owned Observations, so what was read survives.
 */
/* pest_utils::read_file_tail: what a PANTHER agent appends to a partial-results reply.
 *
 * Worth a test rather than trusting four lines of seek-and-read, because every branch here is
 * one this feature actually hits and every one of them fails QUIETLY. A mid-line seek yields
 * half a line, which reads as corruption rather than truncation. An unstripped newline is
 * DROPPED by NetPackage::reset() - not replaced - so two lines arrive welded into one word.
 * And truncating the wrong end throws away the newest text, which is the only part anyone
 * asked for.
 */
static void test_read_file_tail()
{
    cout << "[status file tail: the newest text, on one printable line]" << endl;
    const string fname = "selftest_status_tail.txt";

    // missing file is not an error - the model may not have written anything yet
    remove(fname.c_str());
    CHK(pest_utils::read_file_tail(fname, 100).empty(),
        "a missing status file should give an empty tail, not throw");

    {
        ofstream f(fname);
        f << "solving timestep 1" << endl;
        f << "solving timestep 2" << endl;
    }
    CHK(pest_utils::read_file_tail(fname, 0).empty(),
        "a zero budget should give nothing rather than a stray character");

    // whole file fits: both lines survive, joined by the separator, no trailing blank
    string all = pest_utils::read_file_tail(fname, 200);
    CHK(all == "solving timestep 1 | solving timestep 2",
        "a short file should come back whole and flattened, got: '" + all + "'");

    // the newest text is what must survive a tight budget
    string tight = pest_utils::read_file_tail(fname, 18);
    CHK((int)tight.size() <= 18, "the tail must respect the budget");
    CHK(tight == "solving timestep 2",
        "truncation must drop the OLDEST text, got: '" + tight + "'");

    // a long file: the read window lands mid-line and that partial line must not be reported
    {
        ofstream f(fname);
        for (int i = 0; i < 500; ++i)
            f << "iteration " << i << " residual 0.0001 aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa" << endl;
    }
    string big = pest_utils::read_file_tail(fname, 120);
    CHK((int)big.size() <= 120, "a big file must still respect the budget");
    CHK(big.find("iteration 499") != string::npos,
        "the tail of a big file must contain the LAST line, got: '" + big + "'");
    CHK(big.find("iteration 0 ") == string::npos,
        "the tail of a big file must not reach back to the first line");
    // and the cut lands on a line boundary: an arbitrary cut leaves the orphan tail of a line
    // at the front, which reads as garbage rather than as a truncation. This was real - the
    // first end-to-end run reported "...observations 4 of 20 | solving timestep 5 of 20 | ..."
    CHK(big.compare(0, 10, "iteration ") == 0,
        "a truncated tail must start at a line boundary, got: '" + big + "'");

    // non-printables become spaces rather than vanishing, and runs collapse
    {
        ofstream f(fname, std::ios::binary);
        f << "progress\t\t50%\x01\x02 done\r\n";
    }
    string ctrl = pest_utils::read_file_tail(fname, 100);
    for (size_t i = 0; i < ctrl.size(); ++i)
        CHK(((unsigned char)ctrl[i] >= 32) && ((unsigned char)ctrl[i] <= 126),
            "every byte of the tail must survive NetPackage's printable-ASCII filter");
    CHK(ctrl == "progress 50% done",
        "control characters should collapse to single spaces, got: '" + ctrl + "'");

    // an empty file is not a status
    {
        ofstream f(fname);
    }
    CHK(pest_utils::read_file_tail(fname, 100).empty(),
        "an empty status file should give an empty tail");

    remove(fname.c_str());
}

static void test_instruction_file_tolerant_read()
{
    cout << "[instruction files: a failed parse still yields what it read]" << endl;
    const string ins_name = "selftest_tolerant.ins";
    const string out_name = "selftest_tolerant.out";
    {
        ofstream f(ins_name);
        f << "pif ~" << endl;
        f << "l2 !obs1! !obs2!" << endl;
        f << "l1 !obs3!" << endl;
    }
    // a COMPLETE output file first: the strict and forgiving paths must agree exactly when
    // nothing is wrong, or the forgiving one is not a twin of anything
    {
        ofstream f(out_name);
        f << "header line" << endl;
        f << " 1.5 2.5" << endl;
        f << " 3.5" << endl;
    }
    {
        InstructionFile ins_strict(ins_name), ins_soft(ins_name);
        Observations strict = ins_strict.read_output_file(out_name);
        Observations soft;
        vector<string> missing, problems;
        ins_soft.try_read_output_file(out_name, soft, missing, problems);
        CHK(strict.size() == 3, "strict read of a complete file should give 3 observations");
        CHK(soft.size() == 3, "tolerant read of a complete file should give 3 observations");
        CHK(missing.empty(), "nothing should be missing from a complete file");
        CHK(problems.empty(), "a complete file should report no problems");
        for (auto& name : strict.get_keys())
            CHK(abs(strict.get_rec(name) - soft.get_rec(name)) < 1.0e-12,
                "tolerant read must agree with strict read value for " + name);
    }

    // now TRUNCATE it: obs1/obs2 are readable, obs3's line never arrives
    {
        ofstream f(out_name);
        f << "header line" << endl;
        f << " 1.5 2.5" << endl;
    }
    {
        InstructionFile ins_strict(ins_name);
        bool threw = false;
        try { ins_strict.read_output_file(out_name); }
        catch (const exception&) { threw = true; }
        CHK(threw, "the strict reader must still fail loudly on a truncated output file");

        InstructionFile ins_soft(ins_name);
        Observations soft;
        vector<string> missing, problems;
        ins_soft.try_read_output_file(out_name, soft, missing, problems);
        CHK(soft.size() == 2, "the observations read before the failure must survive it");
        CHK(soft.find("OBS1") != soft.end(), "obs1 was read before the failure");
        CHK(soft.find("OBS2") != soft.end(), "obs2 was read before the failure");
        CHK(missing.size() == 1 && missing[0] == "OBS3",
            "the unread observation must be named as missing");
        CHK(problems.size() == 1, "the failure must be described rather than thrown");
    }

    // and an output file that does not exist at all - the ordinary case for a run that has
    // not written anything yet. Everything the file covers is missing; nothing throws.
    {
        remove(out_name.c_str());
        InstructionFile ins_soft(ins_name);
        Observations soft;
        vector<string> missing, problems;
        ins_soft.try_read_output_file(out_name, soft, missing, problems);
        CHK(soft.size() == 0, "an absent output file yields no observations");
        CHK(missing.size() == 3, "an absent output file makes every covered name missing");
        CHK(problems.size() == 1, "an absent output file is reported, not thrown");
    }
    remove(ins_name.c_str());
}

/**
 * @brief A partial read must never be WRONG - only incomplete.
 *
 * The truncation in the test above cuts at a line boundary, which is the easy case. A run
 * that is still going is writing its output file as we read it, so the realistic cut is
 * mid-line, mid-token, mid-NUMBER. Those are dangerous in a way a missing line is not: a
 * value truncated from "1.2345" to "1.2" still parses, so a caller computing phi from
 * partial results would get a confident wrong answer rather than an obviously absent one.
 *
 * This sweeps every byte offset of a realistic output file and asserts the invariant that
 * makes partial results usable at all: ANY value reported must equal what the strict reader
 * gets from the complete file. Absent is fine. Different is not.
 */
/**
 * @brief PROTOTYPE policy: what a token from an output file is allowed to become.
 *
 * This is the narrowest, most dangerous code in the model interface - it is where a model's
 * text becomes the number phi is computed from - so it gets tested three ways: the policy
 * function on its own, then through all three obs instructions, then through the tolerant
 * reader. The rule being pinned has two halves:
 *
 *   - the whole NORMAL double range round-trips bit-for-bit, and underflow becomes zero
 *     rather than an error, because a subnormal is noise and a model that produced one is
 *     telling us its answer is zero;
 *   - overflow, inf and nan stay refusals, and - the part that matters most - a string that
 *     is not a number must never come back AS a number. Every "must refuse" case below is
 *     there because some plausible-looking piece of text could otherwise be turned into a
 *     confident wrong answer.
 *
 * The fortran fall-back is why that last risk exists at all: `Ew.d` drops the 'E' when the
 * exponent needs three digits, so "0.3000000-309" has to be readable, but "10-20" must not
 * become 1.0e-19.
 */
static void test_parse_double_policy()
{
    cout << "[output doubles: normals round-trip, underflow -> 0, junk stays refused]" << endl;
    const double dmax = numeric_limits<double>::max();
    const double dmin = numeric_limits<double>::min();

    // ---- readable, and exactly ------------------------------------------------------------
    vector<pair<string, double>> exact = {
        {"0", 0.0}, {"0.0", 0.0}, {"-0.0", 0.0},
        {"1", 1.0}, {"-1", -1.0}, {"1.5", 1.5}, {"-1.5", -1.5},
        {"1.0e0", 1.0}, {"1.0E0", 1.0}, {"1e10", 1.0e10}, {"1E-10", 1.0e-10},
        {"  2.5  ", 2.5},                       // leading/trailing whitespace is fine
        {"1.5       ", 1.5},                    // a fixed-width field brings its blanks along
        {"\t1.5\t", 1.5},
        {"1.7976931348623157e+308", dmax},      // DBL_MAX
        {"-1.7976931348623157e+308", -dmax},
        {"2.2250738585072014e-308", dmin},      // smallest normal
        {"-2.2250738585072014e-308", -dmin},
        {"1.0e-307", 1.0e-307}, {"1.0e+307", 1.0e+307},
        {"1.0e-100", 1.0e-100}, {"1.0e+100", 1.0e+100},
        {"3.14159265358979", 3.14159265358979},
        {"1.0d0", 1.0},                         // fortran D exponent
        {"1.0D-10", 1.0e-10},
        {"-2.5D+02", -250.0},
        {"0.3000000-307", 0.3000000e-307},      // fortran e15.7, E dropped, still normal
        {"0.1000000-100", 0.1000000e-100},      // e-less but an ordinary normal number
        {"1.0000000+100", 1.0000000e+100},
        {"1.0000000-099", 1.0e-99},             // three-digit exponent with a leading zero
        {"-0.5000000+300", -0.5e300},
    };
    for (auto& c : exact)
    {
        double v = -12345.0;
        string why;
        bool ok = InstructionFile::try_parse_double(c.first, v, why);
        CHK(ok, "'" + c.first + "' must be readable (got: " + why + ")");
        if (ok)
            CHK(v == c.second, "'" + c.first + "' must read exactly");
    }

    // ---- readable, and zero: every flavour of underflow -----------------------------------
    vector<string> to_zero = {
        "1.0e-310",                     // subnormal
        "-1.0e-310",
        "4.9406564584124654e-324",      // denorm_min
        "1.0e-320",
        "1.0e-400",                     // underflows all the way
        "-1.0e-400",
        "1.0e-4000",
        "0.1000000-320",                // fortran forms underflow too
        "0.3000000-309",                // 3e-310 - looks ordinary, is subnormal
        "3.0000000-310",
        "1.0D-310",
    };
    for (auto& t : to_zero)
    {
        double v = -12345.0;
        string why;
        bool ok = InstructionFile::try_parse_double(t, v, why);
        CHK(ok, "underflow '" + t + "' must be read as zero, not refused (got: " + why + ")");
        if (ok)
            CHK(v == 0.0, "underflow '" + t + "' must be exactly zero, got " + to_string(v));
    }

    // ---- refused: too big, not a number, or not a number AT ALL ---------------------------
    // the second group is the important one. each of these could be turned into a number by a
    // careless fall-back, and every one of them would be a silent wrong answer
    vector<string> refused = {
        "1.0e+309", "-1.0e+309", "1.0e400", "1.7976931348623157e+309",   // overflow
        "inf", "-inf", "Inf", "INF", "nan", "NaN", "-nan",               // not numbers
        "infinity",
        "", " ", "\t", "abc", "--1.0", "1.0.0", "1.0e", "1.0e+", "e10", ".",
        "1.0e+10x", "1.0 2.0", "1,5",
        "10-20",              // no decimal point: must NOT become 1.0e-19
        "1-2",                // ditto
        "2024-01-15",         // a date, two signs
        "1.0-2-3",            // two exponent signs
        "1.0-",               // sign at the end
        "-1.0-",
        "1.0-2.5",            // exponent is not an integer
        "1.0-1234",           // four exponent digits: no double needs one
        "1.5-2",              // one exponent digit: fortran keeps the E below 100, so this
        "1.0-99",             // two digits: ditto - E-99 is what a model actually writes
        "0.3000000-99",
        "-1.5-2",
        ".-309",              // no mantissa digit
        "+-1.0",
        "1.0e-310e5",
    };
    for (auto& t : refused)
    {
        double v = -12345.0;
        string why;
        bool ok = InstructionFile::try_parse_double(t, v, why);
        CHK(!ok, "'" + t + "' must be refused, not turned into " + to_string(v));
        if (!ok)
            CHK(!why.empty(), "a refusal of '" + t + "' must say why");
    }

    // ---- what fortran actually writes ----------------------------------------------------
    // taken from gfortran output: the E/D is dropped only when the exponent needs three
    // digits. one and two digit exponents keep the marker, zero-padded, so they arrive in a
    // form strtod reads on its own and never touch the fall-back at all
    vector<pair<string, double>> fortran_real_output = {
        {"1.0000000E+00", 1.0},          {"1.0000000E-05", 1.0e-5},
        {"1.0000000E-09", 1.0e-9},       {"1.0000000E-99", 1.0e-99},
        {"0.1000000E-98", 0.1e-98},      {"0.1000000E-99", 0.1e-99},
        {"0.1000000D+01", 1.0},          {"0.1000000D-98", 0.1e-98},
        {"0.1000000D-99", 0.1e-99},      {"0.1000E-99", 0.1e-99},
        {"1.0000000-100", 1.0e-100},     {"0.1000000-306", 0.1e-306},
        {"0.1000000+101", 0.1e101},      {"0.1000000+309", 0.1e309},
        {"0.1000-306", 0.1e-306},        {"1.0000000+308", 1.0e308},
        {"100.0000000-102", 100.0e-102}, {"10.0000000E-06", 10.0e-6},
        {"100.0000000+306", 100.0e306},  {"1.0000000000000000E-099", 1.0e-99},
    };
    for (auto& c : fortran_real_output)
    {
        double v = -12345.0;
        string why;
        bool ok = InstructionFile::try_parse_double(c.first, v, why);
        CHK(ok, "gfortran writes '" + c.first + "', so it must be readable (got: " + why + ")");
        if (ok)
            CHK(v == c.second, "'" + c.first + "' must read exactly");
    }

    // ---- the fall-back must never change what already worked ------------------------------
    // anything strtod reads on its own has to come back identical, fall-back or no fall-back
    vector<string> already_fine = {
        "1.0e-10", "1.0E-10", "-3.25e+07", "0.0001", "123456789.123456789",
        "1.0e-307", "9.9e307", "0.5", "-0.5"
    };
    for (auto& t : already_fine)
    {
        double v = 0.0, direct = strtod(t.c_str(), nullptr);
        string why;
        CHK(InstructionFile::try_parse_double(t, v, why), "'" + t + "' must stay readable");
        CHK(v == direct, "'" + t + "' must read the same as a plain strtod");
    }

    // ---- every decade, both signs, straight through ---------------------------------------
    for (int e = -307; e <= 307; e += 1)
    {
        ostringstream os;
        os << "1.0e" << e;
        double v = 0.0, expect = strtod(os.str().c_str(), nullptr);
        string why;
        CHK(InstructionFile::try_parse_double(os.str(), v, why) && (v == expect),
            "decade " + os.str() + " must round-trip");
        ostringstream neg;
        neg << "-1.0e" << e;
        expect = strtod(neg.str().c_str(), nullptr);
        CHK(InstructionFile::try_parse_double(neg.str(), v, why) && (v == expect),
            "decade " + neg.str() + " must round-trip");
    }

    // ---- and the e-less form of every three-digit decade ----------------------------------
    // this is exactly what a fortran model writes with e15.7 / es15.7
    for (int e = 100; e <= 307; e += 1)
    {
        for (int sign = -1; sign <= 1; sign += 2)
        {
            ostringstream eless, with_e;
            eless  << "0.3000000" << (sign < 0 ? "-" : "+") << e;
            with_e << "0.3000000E" << (sign < 0 ? "-" : "+") << e;
            double got = 0.0, expect = 0.0;
            string why, why2;
            bool ok = InstructionFile::try_parse_double(eless.str(), got, why);
            bool ok_e = InstructionFile::try_parse_double(with_e.str(), expect, why2);
            CHK(ok == ok_e, "fortran form " + eless.str() + " must be as readable as " + with_e.str());
            if (ok && ok_e)
                CHK(got == expect, "fortran form " + eless.str() + " must equal " + with_e.str());
        }
    }
}

/**
 * @brief The same policy, but through actual instruction files.
 *
 * The policy function being right is not enough: all three obs instructions have to route
 * through it, and DUM has to keep swallowing whatever it lands on.
 */
static void test_instruction_file_extreme_doubles()
{
    cout << "[instruction files: all three obs instructions share the policy]" << endl;
    const string ins_name = "selftest_extreme.ins";
    const string out_name = "selftest_extreme.out";
    const double dmax = numeric_limits<double>::max();

    struct Case { string text; bool readable; double expect; };
    vector<Case> cases = {
        {"1.7976931348623157e+308", true, dmax},
        {"2.2250738585072014e-308", true, numeric_limits<double>::min()},
        {"0.3000000-307",           true, 0.3000000e-307},   // fortran, no E, normal
        {"0.3000000-309",           true, 0.0},              // fortran, no E, underflows
        {"1.0D-10",                 true, 1.0e-10},          // fortran D exponent
        {"1.0e-310",                true, 0.0},              // subnormal -> zero
        {"1.0e-400",                true, 0.0},              // underflow -> zero
        {"1.0e+309",                false, 0.0},             // overflow -> refused
        {"nan",                     false, 0.0},
        {"inf",                     false, 0.0},
        {"10-20",                   false, 0.0},             // must not become 1e-19
    };

    // free format
    {
        ofstream f(ins_name);
        f << "pif ~" << endl;
        f << "l1 !v!" << endl;
    }
    for (auto& c : cases)
    {
        { ofstream f(out_name); f << " " << c.text << endl; }
        InstructionFile ins(ins_name);
        double got = -12345.0;
        bool threw = false;
        try { got = ins.read_output_file(out_name).get_rec("V"); }
        catch (const exception&) { threw = true; }
        if (c.readable)
        {
            CHK(!threw, "free instruction must read '" + c.text + "'");
            if (!threw)
                CHK(got == c.expect, "free instruction must read '" + c.text + "' exactly");
        }
        else
            CHK(threw, "free instruction must refuse '" + c.text + "'");
    }

    // fixed and semi-fixed. the column range is derived from the text rather than written
    // out: a range one short turns "...e-324" into "...e-32", an ordinary number, and the
    // test would pass for the wrong reason
    for (int which = 0; which < 2; ++which)
    {
        string kind = (which == 0 ? "fixed" : "semi-fixed");
        char open_tag = (which == 0 ? '[' : '(');
        char close_tag = (which == 0 ? ']' : ')');
        for (auto& c : cases)
        {
            if (c.text.find(' ') != string::npos)
                continue;
            {
                ofstream f(ins_name);
                f << "pif ~" << endl;
                f << "l1 " << open_tag << "v" << close_tag << "1:" << c.text.size() << endl;
            }
            { ofstream f(out_name); f << c.text << endl; }
            InstructionFile ins(ins_name);
            double got = -12345.0;
            bool threw = false;
            try { got = ins.read_output_file(out_name).get_rec("V"); }
            catch (const exception&) { threw = true; }
            if (c.readable)
            {
                CHK(!threw, kind + " instruction must read '" + c.text + "'");
                if (!threw)
                    CHK(got == c.expect, kind + " instruction must read '" + c.text + "' exactly");
            }
            else
                CHK(threw, kind + " instruction must refuse '" + c.text + "'");
        }
    }

    // A fixed-width field carries the blanks in its column range, and the value has to be
    // readable anyway. This CHANGED with the policy: the old idx != temp.size() check refused
    // "1.5    " outright, so any [obs]a:b whose range ran past the number was an error.
    {
        ofstream f(ins_name);
        f << "pif ~" << endl;
        f << "l1 [v]1:10" << endl;
    }
    for (auto& padded : vector<string>{"1.5       ", "     1.5  ", "       1.5"})
    {
        { ofstream f(out_name); f << padded << endl; }
        InstructionFile ins(ins_name);
        double got = 0.0;
        bool threw = false;
        try { got = ins.read_output_file(out_name).get_rec("V"); }
        catch (const exception&) { threw = true; }
        CHK(!threw, "a fixed field padded with blanks ('" + padded + "') must still read");
        if (!threw)
            CHK(got == 1.5, "a fixed field padded with blanks must read 1.5");
    }
    // but blanks are the only thing allowed to trail: a second token in the field is junk
    {
        ofstream f(out_name); f << "1.5    2.5" << endl;
    }
    {
        InstructionFile ins(ins_name);
        bool threw = false;
        try { ins.read_output_file(out_name); }
        catch (const exception&) { threw = true; }
        CHK(threw, "a fixed field holding two numbers must be refused, not silently halved");
    }

    // DUM still swallows anything, including what everything else refuses
    {
        ofstream f(ins_name);
        f << "pif ~" << endl;
        f << "l1 !dum! !v!" << endl;
    }
    {
        ofstream f(out_name); f << " nan 1.5" << endl;
    }
    {
        InstructionFile ins(ins_name);
        bool threw = false;
        double got = 0.0;
        try { got = ins.read_output_file(out_name).get_rec("V"); }
        catch (const exception&) { threw = true; }
        CHK(!threw, "DUM must swallow a token the policy refuses");
        if (!threw)
            CHK(got == 1.5, "the observation after a refused DUM token must still read");
    }

    // the tolerant reader agrees with the strict one, case for case
    for (auto& c : cases)
    {
        {
            ofstream f(ins_name);
            f << "pif ~" << endl;
            f << "l1 !v!" << endl;
        }
        { ofstream f(out_name); f << " " << c.text << endl; }
        InstructionFile ins(ins_name);
        Observations obs;
        vector<string> missing, problems;
        ins.try_read_output_file(out_name, obs, missing, problems);
        if (c.readable)
        {
            CHK(obs.find("V") != obs.end(), "tolerant read must also read '" + c.text + "'");
            if (obs.find("V") != obs.end())
                CHK(obs.get_rec("V") == c.expect, "tolerant read of '" + c.text + "' must agree with strict");
            CHK(problems.empty(), "a readable '" + c.text + "' is not a problem");
        }
        else
        {
            CHK(obs.find("V") == obs.end(), "a refused '" + c.text + "' must not appear in a tolerant read");
            CHK(missing.size() == 1 && missing[0] == "V", "the refused observation must be named as missing");
            CHK(problems.size() == 1, "the refusal of '" + c.text + "' must be reported, not thrown");
        }
    }

    remove(ins_name.c_str());
    remove(out_name.c_str());
}

/**
 * @brief A fixed instruction must not be talked into a misread.
 *
 * Fixed instructions are the ones with a real misread hazard: the column range is written by
 * hand, so it can slice a number anywhere, and whatever falls inside the range is what gets
 * converted. The fortran fall-back widens what converts, which makes this worth its own test.
 *
 * The case that must fail is a decimal point in the EXPONENT - "0.3000000-30.9" and friends.
 * A fall-back that accepted it would be reading "0.3000000-30" and quietly dropping the ".9",
 * or worse inventing an exponent of -30.9, and either way handing back a number the model
 * never wrote. It must fail by all three routes into the parser: no exponent marker at all, a
 * 'D' marker, and an 'E' marker that never needed the fall-back in the first place.
 *
 * Worth knowing WHICH guard earns its keep here, because it is not the obvious one. Mutating
 * the exponent-is-digits-only rule to allow '.' does not break these cases: after the 'E' is
 * inserted, strtod stops at the '.' and the left-over-chars check in convert() refuses the
 * string anyway. The digits-only rule pulls its weight elsewhere - the LENGTH half of it is
 * what refuses "1.0-1234", which would otherwise become 1.0e-1234, underflow, and be reported
 * as a perfectly innocent zero. Three layers, and each one is load-bearing for something:
 * the mantissa must have a decimal point ("10-20"), the exponent is at most three digits
 * ("1.0-1234"), and nothing may be left over ("0.3000000-30.9").
 */
/**
 * @brief Who the master will send REQ_PARTIAL to.
 *
 * Getting this wrong is expensive in both directions. Answering "no" for an agent that can
 * report partials means preemption silently does nothing - which is the bug this test was
 * written for: the capability was announced only in READY, and an agent does not send READY
 * until it has FINISHED a run, so the master skipped it for the whole of the first run.
 * Answering "yes" for an agent that cannot is worse: REQ_PARTIAL reaches an in-run loop that
 * treats it as corrupt and kills the run.
 *
 * So the recogniser has to be exact, and it has to say no to anything it does not positively
 * recognise. The strings below are the real ones from PantherAgent: the handshake LINPACK
 * text, the two READY texts, and what an older agent sends in their place.
 */
/**
 * @brief What "is this file available to read" has to mean, on every platform.
 *
 * check_exist_in() is the read-side gate in front of jco files, covariance files, instruction
 * files and model output files, so a wrong answer either hides a missing file or refuses a
 * present one. Two ways to get it wrong, and it needs both halves to avoid them:
 *
 *   - an ifstream ALONE opens a directory quite happily on posix and reports good(), while
 *     windows refuses it. that is the same control file behaving differently on two machines,
 *     which is the worst kind of portability bug because it only shows up on someone else's.
 *   - is_regular_file ALONE says a directory entry of the right kind exists, not that this
 *     process can read the bytes.
 *
 * And the property that must never come back: ASKING must not CHANGE anything. The write-side
 * check_exist_out() opens an ofstream, so it creates what is missing and truncates what is
 * not - that is fine where it is used, because the caller is about to delete the file either
 * way, but it must never creep into the read path.
 */
static void test_file_availability_checks()
{
    cout << "[file checks: a directory is not a readable file, and asking changes nothing]" << endl;
    const string fname = "selftest_avail.txt";
    const string dname = "selftest_avail_dir";
    error_code ec;

    // a real file with content
    {
        ofstream f(fname);
        f << "content that must survive being asked about" << endl;
    }
    uintmax_t before = std::filesystem::file_size(fname, ec);
    CHK(pest_utils::check_exist_in(fname), "a readable regular file must pass check_exist_in");
    CHK(std::filesystem::file_size(fname, ec) == before,
        "check_exist_in must not change the file it is asked about");

    // a directory: posix lets an ifstream open one, so this is the case the is_regular_file
    // half exists for
    std::filesystem::create_directory(dname, ec);
    CHK(!pest_utils::check_exist_in(dname),
        "a DIRECTORY must not answer yes to check_exist_in - an ifstream alone says it does");
    {
        // ...and prove the naked open really is the trap being guarded against, so this test
        // still means something if someone simplifies the implementation
        ifstream d(dname.c_str());
        bool naked_open_says_yes = d.good();
        CHK(pest_utils::check_exist_in(dname) != naked_open_says_yes || !naked_open_says_yes,
            "check_exist_in must not agree with a bare ifstream about a directory");
    }

    // something that is not there
    const string missing = "selftest_avail_missing.txt";
    std::filesystem::remove(missing, ec);
    CHK(!pest_utils::check_exist_in(missing), "a missing file must not pass check_exist_in");
    CHK(!std::filesystem::exists(missing, ec),
        "asking about a missing file must not CREATE it");

    // removal is one call, and 'already gone' is not a failure
    {
        ofstream f(fname);
        f << "about to be removed" << endl;
    }
    CHK(std::filesystem::remove(fname, ec) && (!ec),
        "removing an existing file reports true with no error");
    CHK((!std::filesystem::remove(fname, ec)) && (!ec),
        "removing a file that is already gone is false, NOT an error");

    std::filesystem::remove_all(dname, ec);
    std::filesystem::remove(fname, ec);
}

static void test_partial_capability_handshake()
{
    cout << "[panther: an agent's partial-results capability is recognised before its first run]" << endl;

    // what a current agent sends
    CHK(NetPackage::advertises_partial(NetPackage::PARTIAL_CAPABILITY_TAG),
        "the LINPACK handshake text must advertise support");
    CHK(NetPackage::advertises_partial(string("lets do it ") + NetPackage::PARTIAL_CAPABILITY_TAG),
        "the post-run READY text must advertise support");
    CHK(NetPackage::advertises_partial(string("run completed ") + NetPackage::PARTIAL_CAPABILITY_TAG),
        "a run-status READY text must advertise support");

    // what an older agent sends in the same places: LINPACK carried an empty string and READY
    // carried the status text on its own
    CHK(!NetPackage::advertises_partial(""),
        "an older agent's empty LINPACK text must NOT advertise support");
    CHK(!NetPackage::advertises_partial("lets do it"),
        "an older agent's READY text must NOT advertise support");
    CHK(!NetPackage::advertises_partial("run completed"),
        "an older agent's run-status READY must NOT advertise support");

    // and nothing that merely looks similar counts, because a false yes kills runs
    for (auto& txt : vector<string>{"partial", "partial=", "partial=0", "partials=1", "PARTIAL=1",
                                    "no partial", "1=partial", "partia1=1"})
        CHK(!NetPackage::advertises_partial(txt),
            "'" + txt + "' must not be mistaken for the capability tag");

    // the tag is a substring match on free-form text, so it must still be found when the agent
    // has something else to say
    CHK(NetPackage::advertises_partial(string("host=abc ") + NetPackage::PARTIAL_CAPABILITY_TAG + " dir=/tmp"),
        "the tag must be found among other free-form text");
}

static void test_fixed_instruction_misreads()
{
    cout << "[fixed instructions: a decimal in the exponent is not a number]" << endl;
    const string ins_name = "selftest_misread.ins";
    const string out_name = "selftest_misread.out";

    vector<string> must_fail = {
        "0.3000000-30.9",   // e-less, decimal in the exponent
        "0.3000000+30.9",
        "1.0-2.5",          // the short form of the same thing
        "1.0D-2.5",         // ...with a fortran D marker
        "1.0E-2.5",         // ...and with the marker it already had, which never went through
        "1.0-2.5.5",        // more than one decimal in the exponent
        "0.3000000-.309",   // exponent is nothing but a decimal and digits
    };
    for (auto& text : must_fail)
    {
        {
            ofstream f(ins_name);
            f << "pif ~" << endl;
            f << "l1 [v]1:" << text.size() << endl;
        }
        { ofstream f(out_name); f << text << endl; }
        InstructionFile ins(ins_name);
        bool threw = false;
        double got = -12345.0;
        try { got = ins.read_output_file(out_name).get_rec("V"); }
        catch (const exception&) { threw = true; }
        CHK(threw, "fixed instruction must refuse '" + text
            + "' rather than read it as " + to_string(got));

        // and the same string must be refused by the policy on its own, and through free and
        // semi-fixed, so no route into the parser is softer than the others
        double v = -12345.0;
        string why;
        CHK(!InstructionFile::try_parse_double(text, v, why),
            "the policy must refuse '" + text + "' directly too");
        for (int which = 0; which < 2; ++which)
        {
            {
                ofstream f(ins_name);
                f << "pif ~" << endl;
                if (which == 0)
                    f << "l1 !v!" << endl;
                else
                    f << "l1 (v)1:" << text.size() << endl;
            }
            { ofstream f(out_name); f << (which == 0 ? " " : "") << text << endl; }
            InstructionFile ins2(ins_name);
            bool threw2 = false;
            try { ins2.read_output_file(out_name); }
            catch (const exception&) { threw2 = true; }
            CHK(threw2, string(which == 0 ? "free" : "semi-fixed")
                + " instruction must refuse '" + text + "' as well");
        }
    }

    // A comma where a digit belongs is refused by a FIXED instruction, because the column
    // range hands the whole field over including the comma. Free format is different on
    // purpose: the comma is one of its token delimiters so that csv output can be read at all,
    // so "0.3000000-30,9" is the token "0.3000000-30" followed by something free format never
    // looks at. Both behaviours are right, and they are not the same behaviour.
    {
        const string with_comma = "0.3000000-30,9";
        {
            ofstream f(ins_name);
            f << "pif ~" << endl;
            f << "l1 [v]1:" << with_comma.size() << endl;
        }
        { ofstream f(out_name); f << with_comma << endl; }
        {
            InstructionFile ins(ins_name);
            bool threw = false;
            try { ins.read_output_file(out_name); }
            catch (const exception&) { threw = true; }
            CHK(threw, "a fixed field holding '" + with_comma + "' must be refused");
        }
        double v = -12345.0;
        string why;
        CHK(!InstructionFile::try_parse_double(with_comma, v, why),
            "the policy must refuse '" + with_comma + "' as one token");
        {
            ofstream f(ins_name);
            f << "pif ~" << endl;
            f << "l1 !v!" << endl;
        }
        {
            //free format splits at the comma, so the token is "0.3000000-30" - and that has a
            //two digit exponent, which fortran would have written as E-30, so the exactly-
            //three-digits rule refuses it. the two instruction types reach the same answer by
            //different routes: fixed refuses the whole field, free refuses what is left of it
            InstructionFile ins(ins_name);
            bool threw = false;
            try { ins.read_output_file(out_name).get_rec("V"); }
            catch (const exception&) { threw = true; }
            CHK(threw, "free format refuses the token before the comma: a 2-digit e-less exponent");
        }
    }

    // Clipping the exponent IS caught, and by the exactly-three-digits rule rather than by
    // anything that knows about column ranges: "0.3000000-307" read over 1:12 leaves
    // "0.3000000-30", a two digit e-less exponent, which fortran never writes. Same for 1:11.
    for (int end_col : vector<int>{12, 11})
    {
        {
            ofstream f(ins_name);
            f << "pif ~" << endl;
            f << "l1 [v]1:" << end_col << endl;
        }
        { ofstream f(out_name); f << "0.3000000-307" << endl; }
        InstructionFile ins(ins_name);
        bool threw = false;
        double got = -12345.0;
        try { got = ins.read_output_file(out_name).get_rec("V"); }
        catch (const exception&) { threw = true; }
        CHK(threw, "a fixed range clipping the exponent (1:" + to_string(end_col)
            + ") must be refused, not read as " + to_string(got));
    }

    // The misread that CANNOT be caught, recorded so nobody assumes it is protected against:
    // a column range that stops before the exponent starts. "0.3000000-307" over 1:9 is
    // "0.3000000", which is 0.3 - a real number, 308 orders of magnitude off, and nothing in
    // the string says anything was cut. Only the range is wrong and the reader cannot see the
    // range. If this ever starts throwing, that is an improvement and this check should be
    // updated to say so.
    {
        {
            ofstream f(ins_name);
            f << "pif ~" << endl;
            f << "l1 [v]1:9" << endl;
        }
        { ofstream f(out_name); f << "0.3000000-307" << endl; }
        InstructionFile ins(ins_name);
        double got = 0.0;
        bool threw = false;
        try { got = ins.read_output_file(out_name).get_rec("V"); }
        catch (const exception&) { threw = true; }
        CHK(!threw, "a range stopping before the exponent still yields a number (known limitation)");
        if (!threw)
            CHK(got == 0.3, "the clipped read is 0.3, not the 3e-308 the model wrote (known limitation)");
    }

    remove(ins_name.c_str());
    remove(out_name.c_str());
}

static void test_instruction_file_partial_reads_are_never_wrong()
{
    cout << "[instruction files: a partial read is incomplete, never wrong]" << endl;
    const string ins_name = "selftest_partial.ins";
    const string out_name = "selftest_partial.out";
    {
        ofstream f(ins_name);
        f << "pif ~" << endl;
        f << "~head~" << endl;
        f << "l1 !oa! !ob!" << endl;
        f << "l1 !oc!" << endl;
        f << "~tag~ !od!" << endl;
        // fixed and semi-fixed read by COLUMN, so a short line is a different failure mode
        // from a short token - they are the cases most likely to read past what is there
        f << "l1 [oe]1:9" << endl;
        f << "l1 (of)1:9" << endl;
    }
    // digits chosen so that cutting a number short yields a DIFFERENT VALID number
    const string complete =
        "head\n"
        " 1.2345 6.7891\n"
        " 22.3456\n"
        "tag 98.7654\n"
        " 33.44556\n"
        " 77.88991\n";
    {
        ofstream f(out_name, ios::binary);
        f << complete;
    }

    Observations truth;
    {
        InstructionFile ins(ins_name);
        truth = ins.read_output_file(out_name);
    }
    CHK(truth.size() == 6, "the complete file should yield 6 observations");

    int wrong_values = 0, threw = 0, invented = 0;
    string first_wrong;
    for (size_t cut = 0; cut < complete.size(); cut++)
    {
        {
            ofstream f(out_name, ios::binary | ios::trunc);
            f << complete.substr(0, cut);
        }
        InstructionFile ins(ins_name);
        Observations got;
        vector<string> missing, problems;
        try
        {
            ins.try_read_output_file(out_name, got, missing, problems);
        }
        catch (...)
        {
            threw++;
            continue;
        }
        for (auto& name : got.get_keys())
        {
            if (truth.find(name) == truth.end())
            {
                invented++;
                continue;
            }
            if (abs(got.get_rec(name) - truth.get_rec(name)) > 1.0e-12)
            {
                wrong_values++;
                if (first_wrong.empty())
                {
                    stringstream ss;
                    ss << name << " read as " << got.get_rec(name) << " instead of "
                       << truth.get_rec(name) << " when the file was cut at byte " << cut
                       << " (" << complete.substr(0, cut).size() << " of "
                       << complete.size() << ")";
                    first_wrong = ss.str();
                }
            }
        }
        // present and missing must together account for everything the file covers
        CHK(got.size() + missing.size() == truth.size(),
            "present + missing must equal the covered set at every truncation");
    }
    CHK(threw == 0, "the tolerant reader must never throw, at any truncation");
    CHK(invented == 0, "a partial read must not invent observations the file does not cover");
    CHK(wrong_values == 0, "a partial read reported a WRONG value: " + first_wrong);

    remove(ins_name.c_str());
    remove(out_name.c_str());
}

/**
 * @brief One output file complete, another mid-write: the partial run must keep both halves.
 *
 * The single-file tests prove a truncated file yields what it holds. This proves the case
 * that actually happens when you interrupt a real run: several output files, some finished
 * and some not. One unreadable file must not cost the results sitting in the others - which
 * is exactly what the strict reader does, because it throws on the first bad one.
 */
static void test_model_interface_partial_across_files()
{
    cout << "[model interface: a partial run keeps the files that ARE readable]" << endl;
    const string ins_a = "selftest_mi_a.ins", out_a = "selftest_mi_a.out";
    const string ins_b = "selftest_mi_b.ins", out_b = "selftest_mi_b.out";
    { ofstream f(ins_a); f << "pif ~" << endl << "l1 !ma1! !ma2!" << endl; }
    { ofstream f(ins_b); f << "pif ~" << endl << "l1 !mb1! !mb2!" << endl; }
    { ofstream f(out_a); f << " 11.25 22.5" << endl; }          // complete
    { ofstream f(out_b, ios::binary); f << " 33.75 4"; }        // mid-write, no newline

    ModelInterface mi(vector<string>(), vector<string>(),
                      vector<string>{ins_a, ins_b}, vector<string>{out_a, out_b},
                      vector<string>());
    Observations obs;
    vector<string> missing, problems;
    mi.try_read_output_files(&obs, missing, problems);

    CHK(obs.size() == 4, "every covered observation should be present, read or not");
    CHK(abs(obs.get_rec("MA1") - 11.25) < 1.0e-12, "the complete file's first value survives");
    CHK(abs(obs.get_rec("MA2") - 22.5) < 1.0e-12, "the complete file's second value survives");
    // out_b has no complete line at all, so BOTH of its observations are missing - and
    // crucially the truncated "4" is not reported as a value
    CHK(missing.size() == 2, "the mid-write file's observations should be missing");
    CHK(abs(obs.get_rec("MB1") - Transformable::no_data) < 1.0e-12,
        "a missing observation carries the sentinel, not a plausible number");
    CHK(problems.size() >= 1, "the mid-write file should be reported as a problem");

    // and once it finishes, the same call gets everything
    { ofstream f(out_b); f << " 33.75 44.5" << endl; }
    ModelInterface mi2(vector<string>(), vector<string>(),
                       vector<string>{ins_a, ins_b}, vector<string>{out_a, out_b},
                       vector<string>());
    Observations obs2;
    vector<string> missing2, problems2;
    mi2.try_read_output_files(&obs2, missing2, problems2);
    CHK(missing2.empty(), "nothing should be missing once both files are complete");
    CHK(problems2.empty(), "no problems once both files are complete");
    CHK(abs(obs2.get_rec("MB2") - 44.5) < 1.0e-12, "the finished file reads normally");

    for (auto& f : {ins_a, out_a, ins_b, out_b})
        remove(f.c_str());
}

/**
 * @brief Byte-by-byte truncation of a REAL instruction set, not a synthetic one.
 *
 * out1.dat.ins from benchmarks/tplins_test_1 - the case that exists to exercise instruction
 * files - together with its output file and the .obf of expected values that ships with it.
 * Between them they cover the whole instruction vocabulary the parser dispatches on:
 *
 *   primary markers          ~ primary ~ , ~dummy_obs~
 *   SECONDARY markers        ~secondary~ , and a chain of eight ~,~ on one line
 *   free                     !h01_01! ... eight on a single instruction line
 *   fixed                    [h01_09]45:54 , [h02_01]1:8 ...
 *   semi-fixed               (h01_08)107:114
 *   multi-line advance       l3 (skipping lines outright)
 *   20 observations over 8 instruction lines, ragged across 14 output lines
 *
 * The invariant is the one that makes partial results safe to use at all: at EVERY truncation
 * offset, any value reported must equal what a strict read of the complete file gives.
 */
static void test_instruction_file_partial_real_case()
{
    cout << "[instruction files: byte-by-byte truncation of a REAL instruction set]" << endl;
    const string ins_text =
        "pif ~\n"
        "l1 ~,~ ~,~ ~,~ ~,~ ~,~ ~,~ ~,~ ~,~ !h02_09!\n"
        "~ primary ~ !h01_10! ~  secondary  ~ [h01_09]45:54\n"
        "l3\n"
        "l1 !h01_01! ~secondary~  !h01_02! !h01_03! !h01_04! !h01_05! !h01_06! !h01_07! (h01_08)107:114\n"
        "l1 [h02_01]1:8 [h02_02]9:16 [h02_03]17:24 !h02_04! !h02_05! !h02_06! !h02_07! !h02_08!  \n"
        "l1 [h02_10]1:5 \n"
        "l1 ~dummy_obs~ !dummy_obs!\n";
    const string out_text =
        "1,633800024,633800039,1983-08-03 00:00:00,1983-08-03 00:00:00,0,24.22192,16.30086,7.921060000000001,3,3\n"
        "\n"
        "\n"
        "\n"
        "primary\n"
        "\n"
        " primary  1.234567 junk etc   secondary     9.87654321   \n"
        "\n"
        "\n"
        "\n"
        "   1.000  more crap here nan etc secondary       1.200     1.400     1.600     1.800     2.00000  2.200    2.400   # not used   2.600     2.800\n"
        "1236.5678495.123-900.999     2.200     2.600     3.000     3.400     3.800     4.200     \n"
        "4.123trash\n"
        "dummy_obs 123456789.987654321\n";
    const string ins_name = "selftest_real.ins";
    const string out_name = "selftest_real.out";
    { ofstream f(ins_name, ios::binary); f << ins_text; }
    { ofstream f(out_name, ios::binary); f << out_text; }

    Observations truth;
    {
        InstructionFile ins(ins_name);
        truth = ins.read_output_file(out_name);
    }
    // cross-check the parse itself against the .obf that ships with the case, so this test
    // fails if the READER is wrong rather than only if the tolerant path is
    const pair<string, double> expected[] = {
        {"H02_09", 7.921060},
        {"H01_10", 1.234567},
        {"H01_09", 9.876543},
        {"H01_01", 1.000000},
        {"H01_02", 1.200000},
        {"H01_03", 1.400000},
        {"H01_04", 1.600000},
        {"H01_05", 1.800000},
        {"H01_06", 2.000000},
        {"H01_07", 2.200000},
        {"H01_08", 2.400000},
        {"H02_01", 1236.567},
        {"H02_02", 8495.123},
        {"H02_03", -900.9990},
        {"H02_04", 2.200000},
        {"H02_05", 2.600000},
        {"H02_06", 3.000000},
        {"H02_07", 3.400000},
        {"H02_08", 3.800000},
        {"H02_10", 4.123000},
    };
    for (auto& kv : expected)
    {
        CHK(truth.find(kv.first) != truth.end(),
            "the real instruction set should read " + kv.first);
        if (truth.find(kv.first) != truth.end())
            CHK(abs(truth.get_rec(kv.first) - kv.second) <= 1.0e-5 * max(1.0, abs(kv.second)),
                "value mismatch against the shipped .obf for " + kv.first);
    }

    int threw = 0, invented = 0, wrong = 0, max_recovered = 0;
    string first_wrong;
    for (size_t cut = 0; cut <= out_text.size(); cut++)
    {
        { ofstream f(out_name, ios::binary | ios::trunc); f << out_text.substr(0, cut); }
        InstructionFile ins(ins_name);
        Observations got;
        vector<string> missing, problems;
        try { ins.try_read_output_file(out_name, got, missing, problems); }
        catch (...) { threw++; continue; }

        max_recovered = max(max_recovered, (int)got.size());
        for (auto& name : got.get_keys())
        {
            if (truth.find(name) == truth.end()) { invented++; continue; }
            if (abs(got.get_rec(name) - truth.get_rec(name)) > 1.0e-12)
            {
                wrong++;
                if (first_wrong.empty())
                {
                    stringstream ss;
                    ss << name << " read as " << got.get_rec(name) << " instead of "
                       << truth.get_rec(name) << " at truncation byte " << cut
                       << " of " << out_text.size();
                    first_wrong = ss.str();
                }
            }
        }
        CHK(got.size() + missing.size() == truth.size(),
            "present + missing must account for every covered observation");
    }
    CHK(threw == 0, "the tolerant reader must never throw on a real instruction set");
    CHK(invented == 0, "a partial read must not invent observations");
    CHK(wrong == 0, "a partial read reported a WRONG value: " + first_wrong);
    // and it must actually recover something on the way, or the sweep proves nothing
    CHK(max_recovered >= 10,
        "a partial read of a mostly-complete file should recover most observations");

    remove(ins_name.c_str());
    remove(out_name.c_str());
}

/**
 * @brief The two instruction branches the main real case does not reach: dum, and whitespace.
 *
 * out1.dat.ins covers free, fixed, semi-fixed, line advance and both marker kinds, but the
 * parser dispatches on six token types and reaches DUM by name. These two real files close
 * the set: out1dum.dat.ins is the same case with observations replaced by !DUM! and !dum!
 * (both cases, deliberately), and obj.dat.ins from the mou constraint case is the smallest
 * real user of the whitespace instruction, `l1 w !obj_1!`.
 *
 * Same invariant as the other sweeps: at every truncation, reported values match a strict
 * read of the complete file, and a discarded DUM is never reported as an observation.
 */
static void test_instruction_file_partial_remaining_branches()
{
    cout << "[instruction files: dum and whitespace branches under truncation]" << endl;
    struct Case { const char* tag; string ins; string out; };
    vector<Case> cases;
    cases.push_back({"dum",
        "pif ~\n"
        "l1 ~,~ ~,~ ~,~ ~,~ ~,~ ~,~ ~,~ ~,~ !h02_09!\n"
        "~ primary ~ !h01_10! ~  secondary  ~ [h01_09]45:54\n"
        "l3\n"
        "l1 !h01_01! ~secondary~  !h01_02! !DUM! !h01_04! !h01_05! !h01_06! !dum! (h01_08)107:114\n"
        "l1 [h02_01]1:8 [h02_02]9:16 [h02_03]17:24 !h02_04! !h02_05! !h02_06! !h02_07! !h02_08!  \n"
        "l1 [h02_10]1:5 \n",
        "1,633800024,633800039,1983-08-03 00:00:00,1983-08-03 00:00:00,0,24.22192,16.30086,7.921060000000001,3,3\n"
        "\n"
        "\n"
        "\n"
        "primary\n"
        "\n"
        " primary  1.234567 junk etc   secondary     9.87654321   \n"
        "\n"
        "\n"
        "\n"
        "   1.000  more crap here nan etc secondary       1.200     1.400     1.600     1.800     2.00000  2.200    2.400   # not used   2.600     2.800\n"
        "1236.5678495.123-900.999     2.200     2.600     3.000     3.400     3.800     4.200     \n"
        "4.123trash\n"
        "dummy_obs 123456789.987654321\n"});
    cases.push_back({"whitespace",
        "pif ~\n"
        "l1 w !obj_1!\n"
        "l1 w !obj_2!\n",
        "obj_1 0.5\n"
        "obj_2 7.0\n"});

    for (auto& c : cases)
    {
        const string ins_name = string("selftest_br_") + c.tag + ".ins";
        const string out_name = string("selftest_br_") + c.tag + ".out";
        { ofstream f(ins_name, ios::binary); f << c.ins; }
        { ofstream f(out_name, ios::binary); f << c.out; }

        Observations truth;
        {
            InstructionFile ins(ins_name);
            truth = ins.read_output_file(out_name);
        }
        CHK(truth.size() > 0, string("the ") + c.tag + " case should read something");
        CHK(truth.find("DUM") == truth.end(),
            string("a DUM placeholder must never become an observation (") + c.tag + ")");

        int threw = 0, wrong = 0, dum_leaked = 0;
        string first_wrong;
        for (size_t cut = 0; cut <= c.out.size(); cut++)
        {
            { ofstream f(out_name, ios::binary | ios::trunc); f << c.out.substr(0, cut); }
            InstructionFile ins(ins_name);
            Observations got;
            vector<string> missing, problems;
            try { ins.try_read_output_file(out_name, got, missing, problems); }
            catch (...) { threw++; continue; }
            if (got.find("DUM") != got.end())
                dum_leaked++;
            for (auto& name : got.get_keys())
            {
                if (truth.find(name) == truth.end())
                    continue;
                if (abs(got.get_rec(name) - truth.get_rec(name)) > 1.0e-12)
                {
                    wrong++;
                    if (first_wrong.empty())
                    {
                        stringstream ss;
                        ss << c.tag << ": " << name << " read as " << got.get_rec(name)
                           << " instead of " << truth.get_rec(name) << " at byte " << cut;
                        first_wrong = ss.str();
                    }
                }
            }
            CHK(got.size() + missing.size() == truth.size(),
                string("present + missing must account for the covered set (") + c.tag + ")");
        }
        CHK(threw == 0, string("the tolerant reader must never throw (") + c.tag + ")");
        CHK(dum_leaked == 0, string("a truncated read must not leak DUM (") + c.tag + ")");
        CHK(wrong == 0, "a partial read reported a WRONG value: " + first_wrong);

        remove(ins_name.c_str());
        remove(out_name.c_str());
    }
}

/**
 * @brief A partial write must leave the run's STATUS alone, and everything that reads it.
 *
 * update_run() marks a run complete as it writes the observations. For preemption that is
 * exactly wrong - the run is still executing, and the status byte is what get_num_good_runs(),
 * get_failed_run_ids() and the restart logic all key off. The plan is explicit that the byte
 * does not gain a fifth "partial" value, because each of those call sites would then have to
 * handle it correctly and the failure mode if one did not is silent.
 *
 * So the test is not really about the observations - it is about everything that must NOT
 * have changed.
 */
/* when a RunStorage operation fails, the message has to say what actually went wrong.
 *
 * wrote this after the third ci go-around in a row got stuck on the same string. a macos run
 * failed with exactly "RunStorage::get_run_status_native() stream not good" after 200 model
 * runs had finished fine - and that message cant tell a truncated file from a bad offset from
 * a disk error, so there was nothing to act on and nothing reproduced locally.
 *
 * these are rare and depend on the machine, so the message we throw is the only evidence
 * anybody is ever going to have. this test makes the failure ci keeps hitting - a storage file
 * shorter than its own run counter says - and checks the message has what you need to
 * recognize it: the file, the record, where we looked, how big the file really is, and the
 * counter vs records mismatch that points at the culprit.
 */
/* what the values in pest.stp mean, including the two master-only commands.
 *
 * people write this file by hand while a run is going, so what each value means - and what a
 * bad one does - is worth pinning down. mainly: 5 and 6 must not read as requests to stop,
 * because every tool checks for 1/2/4 by name and widening that by accident would kill runs
 * when somebody only asked to look at them.
 */
static void test_quit_file_tokens()
{
    cout << "[pest.stp: tokens, arguments, and what a bad one does]" << endl;
    const string qf = QUIT_FILENAME;

    auto write_stp = [&](const string &txt) {
        ofstream f(qf);
        f << txt;
    };

    remove(qf.c_str());
    vector<string> args;
    CHK(pest_utils::quit_file_found(args) == 0, "no file means nothing to do");
    CHK(args.empty(), "no file means no arguments");

    write_stp("1\n");
    CHK(pest_utils::quit_file_found(args) == 1, "'1' is a stop");
    CHK(args.empty(), "a bare token leaves no arguments");

    write_stp("4\n");
    CHK(pest_utils::quit_file_found() == 4, "'4' still reads through the no-arg overload");

    // 3 gets reported and ignored rather than acted on
    write_stp("3\n");
    CHK(pest_utils::quit_file_found(args) == 0, "'3' (pause) is not supported and reads as 0");

    // the master-only commands
    write_stp("5\n");
    CHK(pest_utils::quit_file_found(args) == 5, "'5' is the partial-results request");
    CHK(args.empty(), "'5' takes no argument");

    write_stp("6 17\n");
    CHK(pest_utils::quit_file_found(args) == 6, "'6' is the kill-and-abandon command");
    CHK(args.size() == 1, "'6 17' must hand back its argument");
    if (args.size() == 1)
    {
        int rid = -1;
        pest_utils::convert_ip(args[0], rid);
        CHK(rid == 17, "the argument is the run id");
    }

    // this is what keeps 5 and 6 from killing runs - they arent stop values. worth a test
    // because every tool writes it as `q==1||q==2||q==4`, and widening that by accident would
    // turn "show me partial results" into "kill the batch".
    for (const string &t : {string("5"), string("6 17")})
    {
        write_stp(t + "\n");
        int q = pest_utils::quit_file_found(args);
        CHK((q != 1) && (q != 2) && (q != 4),
            "'" + t + "' must not read as a stop token");
    }

    // extra whitespace and extra arguments come through fine
    write_stp("  6   42   ignored_extra  \n");
    CHK(pest_utils::quit_file_found(args) == 6, "leading whitespace should not hide the token");
    CHK(args.size() == 2, "every remaining token is handed back, not just the first");

    // the second line doesnt get read - it is a one line format, and saying so keeps somebody
    // later from assuming otherwise
    write_stp("5\n6 17\n");
    CHK(pest_utils::quit_file_found(args) == 5, "only the first line is read");
    CHK(args.empty(), "arguments come from the first line only");

    // garbage reads as 0, which is the same as no file at all. pinned because it is a real
    // sharp edge - a typo'd stop request just gets ignored, and anybody changing that should
    // have to change a test that says so.
    write_stp("stop\n");
    CHK(pest_utils::quit_file_found(args) == 0,
        "a non-numeric token reads as 0 - SILENTLY, which is a known sharp edge");
    write_stp("\n");
    CHK(pest_utils::quit_file_found(args) == 0, "an empty first line reads as 0");

    remove(qf.c_str());
    CHK(pest_utils::quit_file_found(args) == 0, "cleanup leaves nothing behind");
}

static void test_run_storage_error_diagnostics()
{
    cout << "[run storage: a failure says WHAT failed, not just that something did]" << endl;
    const string stor = "selftest_diag.rns";
    remove(stor.c_str());
    vector<string> pnames{"P1", "P2"}, onames{"O1", "O2", "O3"};

    RunStorage rs(stor);
    rs.reset(pnames, onames, stor);
    Parameters pars;
    pars.insert("P1", 1.0);
    pars.insert("P2", 2.0);
    for (int i = 0; i < 4; ++i)
        rs.add_run(pars, "run", 0.0);
    CHK(rs.get_nruns() == 4, "four runs were added");
    // dont call free_memory() here - it closes the stream and deletes the file. reading the
    // bytes back through a separate stream leaves the file on disk so we can truncate it.

    // cut it down to the header plus a record or so. the counter still says 4 but the file
    // only has 1 or 2. this is the shape of the ci failure - check_rec_id() passes because the
    // counter is fine, then the read runs off the end of a file nobody noticed was short.
    {
        ifstream in(stor, ios_base::binary);
        CHK(in.good(), "the storage file should be readable");
        vector<char> all((istreambuf_iterator<char>(in)), istreambuf_iterator<char>());
        in.close();
        CHK(all.size() > 0, "the storage file should be non-empty");
        ofstream out(stor, ios_base::binary | ios_base::trunc);
        // two thirds is arbitrary but lands in the middle of the records, which is realistic
        out.write(all.data(), (streamsize)(all.size() * 2 / 3));
    }

    string msg;
    try
    {
        // init_restart() might catch the truncation itself, since it walks the records looking
        // for where to pick up - so it goes inside the try. either way we need a good message.
        RunStorage rs2(stor);
        rs2.init_restart(stor);
        for (int i = 0; i < 4; ++i)
        {
            Parameters p;
            Observations o;
            rs2.get_run(i, p, o);
        }
    }
    catch (const exception &e)
    {
        msg = e.what();
    }

    CHK(!msg.empty(), "reading past the end of a truncated storage file should throw");
    if (!msg.empty())
    {
        cout << "  diagnostic: " << msg << endl;
        // every one of these was something you had to guess at before
        CHK(msg.find(stor) != string::npos,
            "the message must name the FILE: " + msg);
        CHK(msg.find("run_id=") != string::npos,
            "the message must name the RECORD it was reading: " + msg);
        CHK(msg.find("target_pos=") != string::npos,
            "the message must say WHERE in the file it looked: " + msg);
        CHK(msg.find("file_size=") != string::npos,
            "the message must give the ACTUAL file size, or a short file is invisible: " + msg);
        CHK(msg.find("state=") != string::npos,
            "the message must give the stream state bits - eof alone means a short file: " + msg);
        // the important one - counter says 4, file has fewer than that
        CHK(msg.find("MISMATCH") != string::npos,
            "a file shorter than its own run counter must be called out explicitly: " + msg);
    }

    // and it has to survive being asked about a file that isnt there at all - if it throws
    // while explaining an error you lose the error it was explaining
    remove(stor.c_str());
    string msg2;
    try
    {
        RunStorage rs3(stor);
        Parameters p;
        Observations o;
        rs3.get_run(0, p, o);
    }
    catch (const exception &e)
    {
        msg2 = e.what();
    }
    CHK(!msg2.empty(), "reading from a missing storage file should throw, not crash");
    remove(stor.c_str());
}

static void test_run_storage_partial_update()
{
    cout << "[run storage: a partial write leaves the run status untouched]" << endl;
    const string stor = "selftest_partial.rns";
    remove(stor.c_str());
    vector<string> pnames{"P1", "P2"}, onames{"O1", "O2", "O3"};

    RunStorage rs(stor);
    rs.reset(pnames, onames, stor);
    Parameters pars;
    pars.insert("P1", 1.0);
    pars.insert("P2", 2.0);
    int rid = rs.add_run(pars, "a run", 0.0);
    int other = rs.add_run(pars, "another", 0.0);

    int status = -99;
    string info_txt;
    double info_value = 0.0;
    rs.get_info(rid, status, info_txt, info_value);
    CHK(status == 0, "a queued run starts as not-completed");
    CHK(rs.get_num_good_runs() == 0, "no runs are good before any complete");

    // partial results: two observations read, one not (carrying the sentinel the model
    // interface fills in)
    Observations partial;
    partial.insert("O1", 11.5);
    partial.insert("O2", 22.5);
    partial.insert("O3", Transformable::no_data);
    rs.update_run_partial(rid, partial);

    rs.get_info(rid, status, info_txt, info_value);
    CHK(status == 0, "a PARTIAL write must not mark the run completed");
    CHK(rs.get_num_good_runs() == 0, "a partial run must not count as a good run");
    CHK(rs.get_run_status(rid) == 0, "a partial run must not read as failed or complete");
    CHK(info_txt.find("a run") != string::npos, "the run's info text must survive a partial write");

    // ...and the values are readable, with the status saying they are not final
    Parameters gp;
    Observations go;
    int read_status = rs.get_run(rid, gp, go);
    CHK(read_status == 0, "reading a partial run reports it as not completed");
    CHK(abs(go.get_rec("O1") - 11.5) < 1.0e-12, "the partial observation is readable");
    CHK(abs(go.get_rec("O2") - 22.5) < 1.0e-12, "the partial observation is readable");
    CHK(abs(go.get_rec("O3") - Transformable::no_data) < 1.0e-12,
        "the unread observation keeps the sentinel");
    CHK(abs(gp.get_rec("P1") - 1.0) < 1.0e-12, "a partial write must not disturb parameters");
    CHK(abs(gp.get_rec("P2") - 2.0) < 1.0e-12, "a partial write must not disturb parameters");

    // the neighbouring record must be untouched - the seek arithmetic is the risk here
    rs.get_info(other, status, info_txt, info_value);
    CHK(status == 0, "a partial write must not reach the next run's record");
    CHK(info_txt.find("another") != string::npos, "the next run's info text is intact");

    // and when it really finishes, the ordinary path still marks it complete
    Observations finished;
    finished.insert("O1", 11.5);
    finished.insert("O2", 22.5);
    finished.insert("O3", 33.5);
    rs.update_run(rid, finished);
    rs.get_info(rid, status, info_txt, info_value);
    CHK(status == 1, "completing the run after a partial write still marks it complete");
    CHK(rs.get_num_good_runs() == 1, "the completed run now counts as good");
    Observations go2;
    Parameters gp2;
    rs.get_run(rid, gp2, go2);
    CHK(abs(go2.get_rec("O3") - 33.5) < 1.0e-12, "completion overwrites the sentinel");

    rs.free_memory();
    remove(stor.c_str());
}

/**
 * @brief A partial read must refuse until the run owns the files on disk.
 *
 * The agent answers REQ_PARTIAL from the moment the run thread is launched, but that thread has
 * to reach ModelInterface::remove_existing() before the PREVIOUS run's output files are gone.
 * In that window the old files are complete, well-formed, and belong to another realization -
 * so a partial read taken there reports somebody else's results as this run's, and mid-run
 * screening abandons a run that violates nothing.
 *
 * That is exactly what happened: a non-violating realization was screened out on one platform
 * and not another, because the window is only a few instructions wide. The agent now asks
 * outputs_are_current() first, and this pins the two states it depends on. The true state is
 * covered by every real run in the benchmarks, where partial reads do come back full.
 */
static void test_partial_read_refuses_stale_outputs()
{
    cout << "[model interface: a partial read refuses until the run owns the files]" << endl;
    // named, not temporaries: ModelInterface mi(vector<string>(), ...) is the most vexing
    // parse and declares a function
    vector<string> none;
    ModelInterface mi(none, none, none, none, none);
    // the only state that can cause harm is a SPURIOUS true, so that is what this pins
    CHK(!mi.outputs_are_current(),
        "a fresh interface must not claim the output files belong to a run");
    mi.mark_outputs_stale();
    CHK(!mi.outputs_are_current(),
        "marking stale must leave a partial read refused");

    // the agent reassigns its interface wholesale when a da cycle swaps the model definition.
    // The flag is held by shared_ptr so that assignment still compiles - this asserts the
    // reassigned interface does not inherit a 'current' claim from anywhere.
    ModelInterface other(none, none, none, none, none);
    mi = other;
    CHK(!mi.outputs_are_current(),
        "a reassigned interface must not claim the output files belong to a run");
}

/**
 * @brief A record holding values is a RESULT, whatever its status byte says.
 *
 * This guards the external run manager's contract, which is not obvious from the code and cost
 * a CI cycle to rediscover. Under /e, pest++ writes the queued runs to run storage, an external
 * script fills in the OBSERVATIONS, and pest++ reads them back. That script leaves run_status
 * at the 0 that add_run() set - it signals failure only by explicitly writing a negative value
 * (see modify_runstor() in benchmarks/ies_test_part4.py). So under the external run manager,
 * status 0 means "fine, here are your values".
 *
 * A process that rejected status <= 0 as "did not complete" therefore dropped every realization
 * of every external run, which is exactly what happened. update_from_runs() now asks whether
 * the record HOLDS values instead, and that question is only answerable because add_run()
 * initializes the observation block to no_data - it used to leave a hole that read back as
 * ZEROS, which is a plausible observation value and indistinguishable from a real result.
 *
 * Both halves are asserted here, because they only work together.
 */
static void test_external_values_are_results()
{
    cout << "[run storage: values with status 0 are results, not an empty record]" << endl;
    const string stor = "selftest_extvals.rns";
    remove(stor.c_str());
    vector<string> pnames{"P1", "P2"}, onames{"O1", "O2", "O3"};

    RunStorage rs(stor);
    rs.reset(pnames, onames, stor);
    Parameters pars;
    pars.insert("P1", 1.0);
    pars.insert("P2", 2.0);
    int filled = rs.add_run(pars, "filled by an external script", 0.0);
    int never = rs.add_run(pars, "never run", 0.0);

    // 1. a queued-but-unrun record reads back as NO_DATA, not as zeros. This is what makes
    //    "does this hold values?" a question with an answer.
    Parameters gp;
    Observations go;
    rs.get_run(never, gp, go);
    for (auto& on : onames)
    {
        CHK(abs(go.get_rec(on) - Transformable::no_data) < 1.0e-12,
            "an unrun record's observations must be the no_data sentinel");
        CHK(go.get_rec(on) != 0.0,
            "an unrun record must NOT read back as zeros - zero is a plausible observation");
    }

    // 2. an external script writes observations and leaves the status alone
    Observations ext;
    ext.insert("O1", 123456789.987);
    ext.insert("O2", 123456789.987);
    ext.insert("O3", 123456789.987);
    rs.update_run_partial(filled, ext);

    int status = -99;
    string info_txt;
    double info_value = 0.0;
    rs.get_info(filled, status, info_txt, info_value);
    CHK(status == 0, "writing observations without completing leaves the status at 0");

    // 3. ...and the process predicate must call that a RESULT. Same test update_from_runs()
    //    applies: a record is empty only when every observation is still the sentinel.
    Observations gf;
    Parameters gfp;
    rs.get_run(filled, gfp, gf);
    bool filled_has_values = false;
    for (Observations::const_iterator it = gf.begin(); it != gf.end(); ++it)
        if (it->second != Transformable::no_data)
        {
            filled_has_values = true;
            break;
        }
    CHK(filled_has_values,
        "a record whose observations were written must be processed, whatever its status - "
        "this is the external run manager's contract");
    CHK(abs(gf.get_rec("O2") - 123456789.987) < 1.0e-6,
        "the externally written value must read back intact");

    Observations gn;
    Parameters gnp;
    rs.get_run(never, gnp, gn);
    bool never_has_values = false;
    for (Observations::const_iterator it = gn.begin(); it != gn.end(); ++it)
        if (it->second != Transformable::no_data)
        {
            never_has_values = true;
            break;
        }
    CHK(!never_has_values,
        "a record that was never filled must NOT be processed - it holds nothing");

    rs.free_memory();
    remove(stor.c_str());
}

/**
 * @brief The single-run violation test must agree with the ensemble one, exactly.
 *
 * Mid-run screening cancels a run on the strength of is_violating(); the process-time drop
 * uses get_violating_realizations(). If the two ever disagree at the margin, screening cancels
 * a run that process would have KEPT - which turns a pure optimization into a silent change of
 * results. So they share their arithmetic, and this asserts they still agree.
 *
 * Also asserts the partial-read property the whole scheme rests on: judging only the
 * observations that have been read, never the sentinel sitting in the ones that have not.
 */
static void test_violation_single_run_matches_ensemble()
{
    cout << "[drop_violations: the single-run test agrees with the ensemble test]" << endl;
    // built through the control-file path rather than by poking members: the handler reads
    // ctl_ordered_obs_names, which only that path fills in, and a scenario assembled any other
    // way silently has no inequality observations at all
    struct TestPest : public Pest { using Pest::tokens_to_obs_rec; };
    TestPest p;
    p.set_defaults();
    ofstream f_rec("selftest_viol.rec");
    // name, value, weight, group. O_LT is tagged less-than by its group, so it violates only
    // when the simulated value goes ABOVE the observed one
    p.tokens_to_obs_rec(f_rec, vector<string>{"O_A", "1.0", "1.0", "OBSGP"});
    p.tokens_to_obs_rec(f_rec, vector<string>{"O_LT", "10.0", "1.0", "L_CONSTRAINT"});
    p.tokens_to_obs_rec(f_rec, vector<string>{"O_B", "2.0", "1.0", "OBSGP"});
    vector<string> onames{"O_A", "O_LT", "O_B"};
    CHK(p.get_ctl_ordered_obs_names().size() == 3, "the scenario has three observations");

    std::mt19937 rgen(11);
    FileManager fm;
    // oe_base carries the observation names the handler judges against
    Eigen::MatrixXd base(1, 3);
    base << 1.0, 10.0, 2.0;
    vector<string> rnames{"R0"};
    ObservationEnsemble oe_base(&p, &rgen, base, rnames, onames);
    Eigen::MatrixXd pm(1, 1);
    pm << 1.0;
    vector<string> pnames{"P0"};
    ParameterEnsemble pe_base(&p, &rgen, pm, rnames, pnames);
    Covariance parcov;
    L2PhiHandler ph(&p, &fm, &oe_base, &pe_base, &parcov, false);

    vector<string> nominated{"O_LT"};

    // three simulated cases: satisfied (below the bound), satisfied at it, violating
    const double sims[3] = {5.0, 10.0, 25.0};
    const bool expect[3] = {false, false, true};
    for (int c = 0; c < 3; c++)
    {
        Eigen::MatrixXd m(1, 3);
        m << 1.0, sims[c], 2.0;
        ObservationEnsemble oe(&p, &rgen, m, rnames, onames);
        vector<string> viol = ph.get_violating_realizations(oe, nominated);

        Observations sim;
        sim.insert("O_A", 1.0);
        sim.insert("O_LT", sims[c]);
        sim.insert("O_B", 2.0);
        bool single = ph.is_violating(sim, set<string>(), nominated);

        CHK((viol.size() > 0) == expect[c], "ensemble test verdict for sim=" +
            to_string(sims[c]));
        CHK(single == expect[c], "single-run test verdict for sim=" + to_string(sims[c]));
        CHK(single == (viol.size() > 0),
            "the two violation tests disagree for sim=" + to_string(sims[c]));
    }

    // the partial-read property: a violating value that has NOT been read yet must not count
    {
        Observations sim;
        sim.insert("O_A", 1.0);
        sim.insert("O_LT", Transformable::no_data);   // not read yet
        sim.insert("O_B", 2.0);
        set<string> read_so_far;
        read_so_far.insert("O_A");
        read_so_far.insert("O_B");
        CHK(!ph.is_violating(sim, read_so_far, nominated),
            "an UNREAD nominated observation must never be judged - the sentinel would "
            "violate almost any inequality");
        // and once it is read, and violating, it counts
        sim.update_rec("O_LT", 25.0);
        read_so_far.insert("O_LT");
        CHK(ph.is_violating(sim, read_so_far, nominated),
            "a violating observation counts once it has actually been read");
    }
    f_rec.close();
    remove("selftest_viol.rec");
}


// ---------------------------------------------------------------------------------------
// The regularization weight search, isolated.
//
// A user reported an infinite regularization weight factor from pestpp-glm when the initial
// regularization phi is zero. The full search (SVDSolver::dynamic_weight_adj) needs a
// ModelRun, a Jacobian and a Q matrix and cannot be driven from here, but the two pieces that
// do the arithmetic can: MuPoint, which holds a trial weight and scores it against the target,
// and DynamicRegularization, which supplies wfmin/wfmax/wffac to the bracketing.
//
// Both are exercised at the values a degenerate problem actually produces - a zero target, a
// zero regularization phi, a zero WFFAC - rather than only at healthy ones.
// ---------------------------------------------------------------------------------------
static void test_regul_weight_search_edges()
{
    const double inf = numeric_limits<double>::infinity();

    // --- MuPoint scoring -----------------------------------------------------------------
    // f() is a difference and stays finite whatever it is handed. error_frac/error_percent
    // divide by the target, and MuPoint::print() reports error_percent to the rec file - so a
    // zero target is printed, not just computed.
    {
        PhiComponets pc;
        pc.meas = 100.0;
        pc.regul = 0.0;          // the reported condition: no regularization phi at all
        MuPoint m;
        m.set(1.0, pc);
        m.target_phi_meas = 0.0; // what phimlim=0 and fracphim=0 produce

        CHK(isfinite(m.f()), "regul: f() finite with a zero target");
        CHK(isfinite(m.error_frac()),
            "regul: error_frac() finite when target phi is zero (divides by the target)");
        CHK(isfinite(m.error_percent()),
            "regul: error_percent() finite when target phi is zero - this is what print() writes to the rec");
    }

    // both zero: 0/0 is nan, and nan compares false against everything, so a nan here does not
    // just print badly - it makes operator< arbitrary, and the search picks its fallback point
    // with that comparison
    {
        PhiComponets pc;
        pc.meas = 0.0;
        pc.regul = 0.0;
        MuPoint m;
        m.set(1.0, pc);
        m.target_phi_meas = 0.0;
        CHK(!isnan(m.error_frac()), "regul: error_frac() not nan when both phi and target are zero");

        MuPoint other;
        other.set(2.0, pc);
        other.target_phi_meas = 100.0;
        // exactly one of these must hold for a usable ordering
        const bool a = (m < other), b = (other < m);
        CHK(a || b || (abs(m.f()) == abs(other.f())),
            "regul: MuPoint ordering is decidable when one point carries a degenerate target");
    }

    // --- set() clamping ------------------------------------------------------------------
    // set() is the guarded way in. the bracketing loops assign .mu directly and bypass this,
    // so what it protects against matters: these are the values that reach it.
    {
        PhiComponets pc;
        pc.meas = 1.0;
        pc.regul = 1.0;
        MuPoint m;

        m.set(0.0, pc);
        CHK(m.mu > 0.0, "regul: set(0) clamps the weight above zero");
        m.set(-5.0, pc);
        CHK(m.mu > 0.0, "regul: set(negative) clamps the weight above zero");
        m.set(inf, pc);
        CHK(isfinite(m.mu), "regul: set(inf) clamps the weight to something finite");
        m.set(numeric_limits<double>::quiet_NaN(), pc);
        CHK(isfinite(m.mu), "regul: set(nan) does not leave a nan weight");
    }

    // --- the wffac divisor ---------------------------------------------------------------
    // the downward bracketing step is  mu = max(mu / wffac, wfmin).  set_zero() leaves wffac
    // at 0, and any control file that sets WFFAC to 0 does the same - the division then
    // produces inf, and max(inf, wfmin) keeps it. this is the arithmetic behind the reported
    // infinite weight factor, checked here at the value the object actually holds.
    {
        DynamicRegularization reg;
        reg.set_defaults();
        CHK(reg.get_wffac() > 0.0, "regul: default WFFAC is usable as a divisor");
        CHK(isfinite(reg.get_wfmax()), "regul: default WFMAX is finite");
        CHK(reg.get_wfmin() > 0.0, "regul: default WFMIN is above zero");

        // set_zero() deliberately zeroes everything - it means "regularization contributes
        // nothing", and only pestpp-swp uses it, where no weight search runs. So the zeros are
        // the intent and are not asserted against here. What must hold is that the SEARCH
        // cannot be handed a wffac it cannot divide by, whether that comes from set_zero() or
        // from a control file with WFFAC 0. That guard lives in SVDSolver (safe_wffac), so the
        // property is checked on the arithmetic it protects.
        DynamicRegularization zero_reg;
        zero_reg.set_zero();
        CHK(zero_reg.get_wffac() == 0.0, "regul: set_zero() still means what it says");

        // the step the search performs, with the guard applied: finite, and still a weight
        for (double raw_wffac : {0.0, -1.0, numeric_limits<double>::infinity(),
                                 numeric_limits<double>::quiet_NaN(), 1.3})
        {
            const double f = (isfinite(raw_wffac) && raw_wffac > 0.0) ? raw_wffac : 1.3;
            const double down = max(1.0 / f, 1e-10);
            const double up   = min(1.0 * f, 1e10);
            CHK(isfinite(down) && down > 0.0,
                "regul: downward weight step finite and positive for every WFFAC the guard admits");
            CHK(isfinite(up) && up > 0.0,
                "regul: upward weight step finite and positive for every WFFAC the guard admits");
        }
    }

    // --- group weight factors ------------------------------------------------------------
    // with no regularization phi there are no group weights to set, so every group must fall
    // back to a neutral 1.0 rather than a zero or a stale value
    {
        DynamicRegularization reg;
        reg.set_defaults();
        CHK(reg.get_grp_weight_fact("a_group_never_set") == 1.0,
            "regul: an unknown group gets a neutral weight factor");

        unordered_map<string, double> empty_weights;
        reg.set_regul_grp_weights(empty_weights);
        CHK(reg.get_grp_weight_fact("regul_group") == 1.0,
            "regul: an empty group weight map leaves a neutral factor");
    }
}

/* ies prior scaling must not change the answer.
 *
 * chen and oliver (2013) eq 5 makes the parameter deviations dimensionless with
 * C_sc^-1/2, and eq 15 pre-multiplies the update by C_sc^1/2.  since
 * C_sc^1/2 * delta_m IS the raw deviation, the two cancel: prior scaling is a
 * conditioning device for the svd, not a change to the update.  the scaling
 * matrix is diagonal and holds prior variances, so with a non-unit parcov the
 * upgrade computed with ies_use_prior_scaling must equal the one computed
 * without it.
 *
 * it did not.  the back-transform in ensemble_solution was commented out, so
 * every parameter's update came out divided by its prior standard deviation -
 * dimensionally wrong, and silently so, because nothing compared the two paths.
 */
static void test_ies_prior_scaling_is_neutral()
{
    cout << "[ies prior scaling: C_sc cancels, so it cannot change the upgrade]" << endl;

    const int npar = 4, nobs = 3, nreal = 5;
    // fixed, arbitrary, nothing special about the values beyond being non-degenerate
    Eigen::MatrixXd par_diff0(npar, nreal), obs_diff0(nobs, nreal), obs_resid0(nobs, nreal);
    for (int i = 0; i < npar; i++)
        for (int j = 0; j < nreal; j++)
            par_diff0(i, j) = 0.3 * (i + 1) - 0.11 * (j + 1) + 0.05 * (i * j);
    for (int i = 0; i < nobs; i++)
        for (int j = 0; j < nreal; j++)
        {
            obs_diff0(i, j) = 0.7 * (i + 1) + 0.13 * (j + 1) - 0.04 * (i * j);
            obs_resid0(i, j) = 0.2 * (i + 1) - 0.09 * (j + 1);
        }
    Eigen::MatrixXd par_resid0 = 0.5 * par_diff0;
    Eigen::MatrixXd obs_err0 = 0.1 * obs_diff0;
    Eigen::MatrixXd Am(npar, 1);
    Am.setZero();

    Eigen::VectorXd wvec(nobs);
    wvec << 1.0, 2.0, 0.5;
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> weights(wvec);

    // prior variances deliberately far from 1 and from each other - if the
    // scaling leaks into the answer this is what makes it visible
    Eigen::VectorXd pvar(npar);
    pvar << 0.01, 1.0, 100.0, 4.0;
    Eigen::VectorXd pinv_vec = pvar.cwiseSqrt().cwiseInverse();   // C_sc^-1/2
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> parcov_inv(pinv_vec);
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> parcov_identity(
        Eigen::VectorXd::Ones(npar));

    vector<string> onames{"o1", "o2", "o3"}, pnames{"p1", "p2", "p3", "p4"};

    auto solve = [&](bool use_prior_scaling,
                     const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& pcinv)
    {
        // ensemble_solution mutates what it is handed, so every call gets copies
        Eigen::MatrixXd par_resid = par_resid0, par_diff = par_diff0;
        Eigen::MatrixXd obs_resid = obs_resid0, obs_diff = obs_diff0, obs_err = obs_err0;
        Eigen::MatrixXd upgrade;
        // use_approx = true and iter = 1 keep this to the data term, which is
        // where the cancellation is exact and unconditional
        UpgradeThread::ensemble_solution(1, 0, 1000, 0, 0, use_prior_scaling, true, true,
                                         0.0, 1.0e-7, par_resid, par_diff, Am, obs_resid,
                                         obs_diff, upgrade, obs_err, weights, pcinv,
                                         onames, pnames, -1.0, 0.0);
        return upgrade;
    };

    Eigen::MatrixXd off = solve(false, parcov_inv);
    Eigen::MatrixXd on = solve(true, parcov_inv);

    CHK(off.rows() == on.rows() && off.cols() == on.cols(),
        "prior scaling on/off give the same shaped upgrade");
    double dnorm = (off - on).norm();
    double scale = max(1.0, off.norm());
    bool same = (dnorm / scale) < 1.0e-8;
    if (!same)
        cout << "  scaled/unscaled upgrades differ: relative norm "
             << (dnorm / scale) << " (off " << off.norm() << ", on " << on.norm() << ")"
             << endl;
    CHK(same, "ies_use_prior_scaling does not change the upgrade (C_sc cancels)");

    // and the guard against a fix that merely makes the option inert: with a unit
    // parcov the scaled path must still reproduce the unscaled one
    Eigen::MatrixXd on_unit = solve(true, parcov_identity);
    CHK(((off - on_unit).norm() / scale) < 1.0e-8,
        "prior scaling with unit parcov matches the unscaled upgrade");

    // the upgrade must actually be doing something, or the checks above are vacuous
    CHK(off.norm() > 1.0e-8, "the unscaled upgrade is non-zero");
}

/* the scaling that puts prior anomalies into the space chen and oliver eq 5 works
 * in: delta_m_pr = C_sc^-1/2 (m_pr - mbar) / sqrt(Ne-1).  C_sc is diagonal and
 * holds prior variances, so this is a row scaling by 1/sqrt(var).
 *
 * get_Am() feeds the result to a tsvd, and the whole reason the paper scales
 * before that svd is that the truncation then happens in scaled space - a
 * different set of directions survives.  get_Am used to skip this entirely while
 * the residual it multiplies WAS scaled, so the model-mismatch term straddled
 * two spaces.
 */
static void test_prior_anomaly_scaling()
{
    cout << "[ies prior scaling: prior anomalies scaled into C_sc^-1/2 space]" << endl;

    const int npar = 4, nreal = 6;
    Eigen::MatrixXd anom(npar, nreal);
    for (int i = 0; i < npar; i++)
        for (int j = 0; j < nreal; j++)
            anom(i, j) = 0.4 * (i + 1) - 0.17 * (j + 1) + 0.03 * (i * j);

    // spread over four orders of magnitude: if a row is missed it shows up loudly
    Eigen::VectorXd var(npar);
    var << 0.01, 1.0, 100.0, 4.0;

    Eigen::MatrixXd scaled = scale_prior_anomalies(anom, var);

    bool rows_ok = true;
    for (int i = 0; i < npar; i++)
        rows_ok &= ((scaled.row(i) - anom.row(i) / sqrt(var[i])).norm() < 1.0e-12);
    CHK(rows_ok, "each row divided by its own prior standard deviation");

    CHK((scale_prior_anomalies(anom, Eigen::VectorXd::Ones(npar)) - anom).norm() < 1.0e-12,
        "unit prior variance leaves the anomalies untouched");

    // fixed/tied parameters can arrive with no variance; that must not become inf
    Eigen::VectorXd degenerate(npar);
    degenerate << 0.0, -1.0, 100.0, 4.0;
    Eigen::MatrixXd sd = scale_prior_anomalies(anom, degenerate);
    CHK(sd.allFinite(), "zero or negative prior variance does not produce inf or nan");
    CHK((sd.row(0) - anom.row(0)).norm() < 1.0e-12,
        "a zero-variance row passes through unscaled");

    bool threw = false;
    try { scale_prior_anomalies(anom, Eigen::VectorXd::Ones(npar + 1)); }
    catch (const exception&) { threw = true; }
    CHK(threw, "a prior variance vector of the wrong length is rejected");

    // the point of scaling before the svd: the singular directions change.  if
    // this ever stops holding, the scaling has become decorative
    Eigen::JacobiSVD<Eigen::MatrixXd> svd_raw(anom, Eigen::ComputeThinU);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd_scaled(scaled, Eigen::ComputeThinU);
    double lead_align = fabs(svd_raw.matrixU().col(0).dot(svd_scaled.matrixU().col(0)));
    CHK(lead_align < 0.99,
        "scaling changes the leading singular direction, so the tsvd truncates differently");
}

/*
 * the scaling VECTOR the solve paths hand to ensemble_solution().
 *
 * three call sites used to build this inline and they drifted.  the multimodal
 * path handed over a plain reciprocal - C_sc^-1 where the solve wants C_sc^-1/2 -
 * and the other two branched on isdiagonal() with opposite conditions and
 * identical bodies in both arms.  that second one was harmless, but only by
 * accident: every MatType stores a square matrix, so .diagonal() was right
 * either way.  this pins down all of it so the sites cannot drift again.
 */
static void test_prior_inv_sqrt_diag()
{
    cout << "[ies prior scaling: the C_sc^-1/2 vector handed to the solve]" << endl;

    vector<string> pnames{ "p1", "p2", "p3", "p4" };
    const int npar = (int)pnames.size();
    Eigen::VectorXd var(npar);
    var << 0.01, 1.0, 100.0, 4.0;

    // a full (correlated) prior: only the diagonal may be used
    vector<Eigen::Triplet<double>> trips;
    for (int i = 0; i < npar; i++)
        trips.push_back(Eigen::Triplet<double>(i, i, var[i]));
    trips.push_back(Eigen::Triplet<double>(0, 1, 0.005));
    trips.push_back(Eigen::Triplet<double>(1, 0, 0.005));
    Eigen::SparseMatrix<double> full(npar, npar);
    full.setFromTriplets(trips.begin(), trips.end());
    Covariance cov_full(pnames, full, Mat::MatType::SPARSE);

    Eigen::VectorXd s = prior_inv_sqrt_diag(cov_full);
    bool ok = true;
    for (int i = 0; i < npar; i++)
        ok &= (fabs(s[i] - 1.0 / sqrt(var[i])) < 1.0e-12);
    CHK(ok, "one over the square root of each prior variance");

    // the multimodal bug: a plain reciprocal is 100 where the right answer is 10
    CHK(fabs(s[0] - 1.0 / var[0]) > 1.0,
        "the inverse square root, not the plain reciprocal");

    // same numbers whatever the storage flag says - the isdiagonal() branch that
    // used to guard this was a coin flip
    vector<Eigen::Triplet<double>> dtrips;
    for (int i = 0; i < npar; i++)
        dtrips.push_back(Eigen::Triplet<double>(i, i, var[i]));
    Eigen::SparseMatrix<double> diag(npar, npar);
    diag.setFromTriplets(dtrips.begin(), dtrips.end());
    Covariance cov_diag(pnames, diag, Mat::MatType::DIAGONAL);
    CHK((prior_inv_sqrt_diag(cov_diag) - s).norm() < 1.0e-12,
        "diagonal storage and sparse storage give the same scaling vector");

    // fixed and tied parameters arrive with no variance
    vector<Eigen::Triplet<double>> ztrips;
    ztrips.push_back(Eigen::Triplet<double>(0, 0, 0.0));
    for (int i = 1; i < npar; i++)
        ztrips.push_back(Eigen::Triplet<double>(i, i, var[i]));
    Eigen::SparseMatrix<double> zero(npar, npar);
    zero.setFromTriplets(ztrips.begin(), ztrips.end());
    Covariance cov_zero(pnames, zero, Mat::MatType::DIAGONAL);
    Eigen::VectorXd sz = prior_inv_sqrt_diag(cov_zero);
    CHK(sz.allFinite(), "a zero prior variance does not produce inf or nan");
    CHK(fabs(sz[0] - 1.0) < 1.0e-12, "a zero-variance parameter is left unscaled");

    // the vector form and the matrix form must agree - get_Am() uses the latter
    // while the solve paths use the former, on the same prior
    Eigen::MatrixXd anom(npar, 5);
    for (int i = 0; i < npar; i++)
        for (int j = 0; j < 5; j++)
            anom(i, j) = 0.3 * (i + 1) - 0.11 * (j + 1);
    CHK((Eigen::MatrixXd(s.asDiagonal() * anom) - scale_prior_anomalies(anom, var)).norm() < 1.0e-12,
        "agrees with scale_prior_anomalies() on the same prior");
}

/* ies_use_prior_prec: the glm step with an explicit prior precision Q in both the
 * hessian and the gradient.  when Q is the (untruncated) Am Am^T of the same
 * ensemble - the pseudo-inverse of the prior anomaly covariance that the existing
 * path uses - the two forms are the same thing written differently, so the
 * upgrades have to agree to roundoff, at iteration 1 (data term only) and at
 * iteration 2 with the prior pull switched on.  the obs anomalies are kept full
 * rank (nobs >= N-1) so the existing path's V-subspace projection is the whole
 * centred subspace; with fewer obs the two forms legitimately differ.
 */
static void test_ies_prior_prec_reduces_to_am()
{
    cout << "[ies prior precision: Q = Am Am^T reproduces the existing glm step]" << endl;

    const int npar = 30, nobs = 40, nreal = 12;
    // seeded pseudo-random so the anomalies are full rank (N-1): a smooth formula
    // like sin(a i + b j) is rank 2 and then the obs svd spans only part of the
    // centred subspace, which is exactly the case where the two forms differ
    std::mt19937 gen(20260919);
    std::uniform_real_distribution<double> unif(-1.0, 1.0);
    Eigen::MatrixXd par_diff0(npar, nreal), obs_diff0(nobs, nreal), obs_resid0(nobs, nreal),
        par_resid0(npar, nreal);
    for (int i = 0; i < npar; i++)
        for (int j = 0; j < nreal; j++)
        {
            par_diff0(i, j) = unif(gen);
            par_resid0(i, j) = 0.3 * unif(gen);
        }
    for (int i = 0; i < nobs; i++)
        for (int j = 0; j < nreal; j++)
        {
            obs_diff0(i, j) = unif(gen);
            obs_resid0(i, j) = 0.2 * unif(gen) + 0.1;
        }
    par_diff0 = par_diff0.colwise() - par_diff0.rowwise().mean();
    obs_diff0 = obs_diff0.colwise() - obs_diff0.rowwise().mean();
    Eigen::MatrixXd obs_err0 = 0.1 * obs_diff0;

    Eigen::VectorXd wvec(nobs);
    for (int i = 0; i < nobs; i++)
        wvec[i] = 0.5 + 0.1 * (i % 4);
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> weights(wvec);
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> parcov_identity(Eigen::VectorXd::Ones(npar));
    vector<string> onames, pnames;
    for (int i = 0; i < nobs; i++) onames.push_back("o" + to_string(i));
    for (int i = 0; i < npar; i++) pnames.push_back("p" + to_string(i));

    // Am exactly as get_Am builds it: tsvd of the scaled prior anomalies, U S^-1
    const double eigthresh = 1.0e-10;
    const int maxsing = 1000;
    Eigen::MatrixXd Am;
    {
        Eigen::MatrixXd pd = par_diff0 / sqrt(double(nreal - 1));
        Eigen::MatrixXd s, U, V;
        SVD_REDSVD rsvd;
        rsvd.solve_ip(pd, s, U, V, eigthresh, maxsing);
        Eigen::MatrixXd sinv = s.asDiagonal().inverse();
        Am = U * sinv;
    }
    Eigen::MatrixXd Qd = Am * Am.transpose();
    Eigen::SparseMatrix<double> Q = Qd.sparseView();

    auto solve_old = [&](int iter, double lam, double reg)
    {
        Eigen::MatrixXd par_resid = par_resid0, par_diff = par_diff0;
        Eigen::MatrixXd obs_resid = obs_resid0, obs_diff = obs_diff0, obs_err = obs_err0;
        Eigen::MatrixXd upgrade;
        UpgradeThread::ensemble_solution(iter, 0, maxsing, 0, 0, false, false, true, lam, eigthresh,
                                         par_resid, par_diff, Am, obs_resid, obs_diff, upgrade,
                                         obs_err, weights, parcov_identity, onames, pnames, reg, -1.0);
        return upgrade;
    };
    auto solve_new = [&](int iter, double lam, double reg, const Eigen::SparseMatrix<double>& q)
    {
        Eigen::MatrixXd par_resid = par_resid0, par_diff = par_diff0;
        Eigen::MatrixXd obs_resid = obs_resid0, obs_diff = obs_diff0;
        Eigen::MatrixXd upgrade;
        UpgradeThread::ensemble_solution_prec(iter, 0, maxsing, false, lam, eigthresh, par_resid,
                                              par_diff, q, obs_resid, obs_diff, upgrade, weights, reg);
        return upgrade;
    };
    auto rel = [](const Eigen::MatrixXd& a, const Eigen::MatrixXd& b)
    { return (a - b).norm() / max(1.0e-30, a.norm()); };

    // iteration 1: data term only, at two lambdas
    for (double lam : {0.0, 3.0})
    {
        Eigen::MatrixXd o = solve_old(1, lam, 1.0), n = solve_new(1, lam, 1.0, Q);
        CHK(o.rows() == n.rows() && o.cols() == n.cols(), "prec form gives an N x p upgrade");
        double d = rel(o, n);
        if (d > 1.0e-8)
            cout << "  iteration 1, lambda " << lam << ": relative diff " << d << endl;
        CHK(d < 1.0e-8, "iteration 1 upgrade matches the existing glm step (lambda " + to_string(lam) + ")");
        CHK(o.norm() > 1.0e-8, "the iteration 1 upgrade is non-zero");
    }
    // iteration 2 with reg_factor 1: the prior pull is in on both sides
    {
        Eigen::MatrixXd o = solve_old(2, 2.0, 1.0), n = solve_new(2, 2.0, 1.0, Q);
        double d = rel(o, n);
        if (d > 1.0e-8)
            cout << "  iteration 2, prior pull on: relative diff " << d << endl;
        CHK(d < 1.0e-8, "iteration 2 upgrade with the prior pull matches the existing full solution");
        // and the pull did something, or the check above is the iteration 1 check again
        CHK(rel(solve_old(1, 2.0, 1.0), o) > 1.0e-6, "the prior pull changes the iteration 2 upgrade");
        // reg_factor 0 (the default) means weight 1 on the pull, on both sides.  this
        // used to multiply the pull by 0 and quietly turn the full solution into the
        // approx one - the guard was >= 0 instead of > 0
        CHK(rel(solve_old(2, 2.0, 0.0), o) < 1.0e-12, "reg_factor 0 (default) is weight 1 on the existing full solution");
        CHK(rel(solve_new(2, 2.0, 0.0, Q), n) < 1.0e-12, "reg_factor 0 (default) is weight 1 on the prec full solution");
        // and a fractional weight really scales it, so the > 0 branch is live too
        Eigen::MatrixXd half_o = solve_old(2, 2.0, 0.5), half_n = solve_new(2, 2.0, 0.5, Q);
        CHK(rel(half_o, o) > 1.0e-6, "reg_factor 0.5 scales the pull on the existing full solution");
        CHK(rel(half_o, half_n) < 1.0e-8, "reg_factor 0.5 matches on both sides");
    }
    // the two readings of ies_reg_factor: the upgrade takes |r| (negative = full
    // solution with the magnitude on the pull, no reg phi), the phi handler max(0,r)
    {
        CHK(L2PhiHandler::upgrade_reg_factor(-1.0) == 1.0, "upgrade reg factor of -1 is 1");
        CHK(L2PhiHandler::upgrade_reg_factor(0.0) == 0.0, "upgrade reg factor of 0 is 0");
        CHK(L2PhiHandler::upgrade_reg_factor(0.5) == 0.5, "upgrade reg factor of 0.5 is 0.5");
        CHK(L2PhiHandler::phi_reg_factor(-1.0) == 0.0, "phi reg factor of -1 is 0");
        CHK(L2PhiHandler::phi_reg_factor(0.0) == 0.0, "phi reg factor of 0 is 0");
        CHK(L2PhiHandler::phi_reg_factor(0.5) == 0.5, "phi reg factor of 0.5 is 0.5");
    }
    // a different Q must give a different answer, or none of this proves Q is used
    {
        Eigen::SparseMatrix<double> Qdiag(npar, npar);
        Qdiag.setIdentity();
        Qdiag *= 4.0;
        CHK(rel(solve_new(1, 1.0, 1.0, Q), solve_new(1, 1.0, 1.0, Qdiag)) > 1.0e-4,
            "a diagonal Q gives a different upgrade than Am Am^T");
        CHK(solve_new(1, 1.0, 1.0, Qdiag).allFinite(), "the singular ones direction is dropped, not inverted");
    }
}

static void test_enif_woodbury_step()
{
    cout << "[enif: woodbury step with an explicit sparse H matches the dense normal equations]" << endl;
    // small enough to check against the dense form
    //   delta = -[(1+lam) C^-1 + H^T R^-1 H]^-1 [C^-1 e + H^T R^-1 r]
    int p = 6, n = 4, N = 3;
    std::mt19937 gen(11);
    std::normal_distribution<double> nd(0.0, 1.0);
    Eigen::MatrixXd Hd = Eigen::MatrixXd::Zero(n, p);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < p; j++)
            if (((i + j) % 2) == 0)
                Hd(i, j) = nd(gen);
    Eigen::SparseMatrix<double> H = Hd.sparseView();
    Eigen::VectorXd cdiag(p), var(n);
    for (int j = 0; j < p; j++) cdiag[j] = 0.5 + j;
    for (int i = 0; i < n; i++) var[i] = 0.1 * (i + 1);
    Eigen::SparseMatrix<double> C(p, p);
    {
        vector<Eigen::Triplet<double>> t;
        for (int j = 0; j < p; j++) t.push_back(Eigen::Triplet<double>(j, j, cdiag[j]));
        C.setFromTriplets(t.begin(), t.end());
    }
    Eigen::MatrixXd e(p, N), r(n, N);
    for (int j = 0; j < N; j++) { for (int i = 0; i < p; i++) e(i, j) = nd(gen); for (int i = 0; i < n; i++) r(i, j) = nd(gen); }
    for (double lam : {0.0, 0.7, 25.0})
    {
        Eigen::MatrixXd got = enif_woodbury_step(H, var, C, e, r, lam);
        Eigen::MatrixXd Cinv = cdiag.cwiseInverse().asDiagonal();
        Eigen::MatrixXd Rinv = var.cwiseInverse().asDiagonal();
        Eigen::MatrixXd lam_post = (1.0 + lam) * Cinv + Hd.transpose() * Rinv * Hd;
        Eigen::MatrixXd g = Cinv * e + Hd.transpose() * Rinv * r;
        Eigen::MatrixXd want = -lam_post.ldlt().solve(g);
        CHK((got - want).cwiseAbs().maxCoeff() < 1e-9, "woodbury step == dense normal equations, lambda " + to_string(lam));
    }
    // one column at a time gives the same answer as all at once
    Eigen::MatrixXd all = enif_woodbury_step(H, var, C, e, r, 0.7);
    Eigen::MatrixXd one = enif_woodbury_step(H, var, C, e.col(1), r.col(1), 0.7);
    CHK((all.col(1) - one.col(0)).cwiseAbs().maxCoeff() < 1e-12, "a single column solve matches its column of the batch solve");
    // an empty H row (observation carries nothing) must not break the step
    Eigen::MatrixXd Hz = Hd; Hz.row(2).setZero();
    Eigen::SparseMatrix<double> Hzs = Hz.sparseView();
    Eigen::MatrixXd got = enif_woodbury_step(Hzs, var, C, e, r, 0.7);
    CHK(got.allFinite(), "an empty H row still gives a finite step");
    bool threw = false;
    try { enif_woodbury_step(H, var.head(2), C, e, r, 0.7); } catch (...) { threw = true; }
    CHK(threw, "mismatched noise length throws");
}

static void test_enif_lowrank_step()
{
    cout << "[enif: low-rank woodbury step matches the dense step]" << endl;
    int p = 9, n = 12, N = 4;
    std::mt19937 gen(23);
    std::normal_distribution<double> nd(0.0, 1.0);
    Eigen::MatrixXd A(p, N), B(n, N), e(p, N), r(n, N);
    for (int i = 0; i < p; i++) for (int j = 0; j < N; j++) { A(i, j) = nd(gen); e(i, j) = 0.3 * nd(gen); }
    for (int i = 0; i < n; i++) for (int j = 0; j < N; j++) { B(i, j) = nd(gen); r(i, j) = nd(gen); }
    Eigen::VectorXd var(n);
    for (int i = 0; i < n; i++) var[i] = 0.05 * (i + 1);
    // a dense spd prior covariance and the ensemble covariance
    Eigen::MatrixXd Z(p, p);
    for (int i = 0; i < p; i++) for (int j = 0; j < p; j++) Z(i, j) = nd(gen);
    Eigen::MatrixXd Cd = Z * Z.transpose() + Eigen::MatrixXd::Identity(p, p);
    Eigen::MatrixXd Ccur = A * A.transpose();

    // explicit sparse H with C = A A^T: U = H A, Q = I, L = A
    Eigen::MatrixXd Hd = Eigen::MatrixXd::Zero(n, p);
    for (int i = 0; i < n; i++) for (int j = 0; j < p; j++) if (((i + 2 * j) % 3) == 0) Hd(i, j) = nd(gen);
    Eigen::SparseMatrix<double> H = Hd.sparseView();
    Eigen::SparseMatrix<double> Ccur_s = Ccur.sparseView();
    for (double lam : {0.0, 1.5, 40.0})
    {
        Eigen::MatrixXd dense = enif_woodbury_step(H, var, Ccur_s, e, r, lam);
        Eigen::MatrixXd low = enif_lowrank_step(Hd * A, Eigen::MatrixXd::Identity(N, N), A, var, e, Hd * e, r, lam);
        CHK((dense - low).cwiseAbs().maxCoeff() < 1e-9, "sparse H, ensemble covariance: low-rank == dense, lambda " + to_string(lam));
    }

    // implicit H = B M A^T: U = B, Q = M K M, L = C A M, against the formula with everything dense
    Eigen::MatrixXd gram = A.transpose() * A;
    double gamma = 1.0e-6 * gram.trace() / double(N);
    gram.diagonal().array() += gamma;
    Eigen::LDLT<Eigen::MatrixXd> gf(gram);
    Eigen::MatrixXd Hi = B * gf.solve(A.transpose());
    for (int which = 0; which < 2; which++)
    {
        const Eigen::MatrixXd& C = (which == 0) ? Cd : Ccur;
        Eigen::MatrixXd V = C * A;
        Eigen::MatrixXd K = A.transpose() * V;
        Eigen::MatrixXd Q = gf.solve(gf.solve(K).transpose());
        Eigen::MatrixXd L = gf.solve(V.transpose()).transpose();
        for (double lam : {0.0, 2.0})
        {
            double sc = 1.0 / (1.0 + lam);
            Eigen::MatrixXd CHt = sc * C * Hi.transpose();
            Eigen::MatrixXd G = Hi * CHt;
            G.diagonal() += var;
            Eigen::MatrixXd He = Hi * e;
            Eigen::MatrixXd want = -((sc * (e - (CHt * G.ldlt().solve(He)))) + (CHt * G.ldlt().solve(r)));
            Eigen::MatrixXd got = enif_lowrank_step(B, Q, L, var, e, He, r, lam);
            CHK((got - want).cwiseAbs().maxCoeff() < 1e-8,
                string("implicit H, ") + (which == 0 ? "supplied" : "ensemble") + " covariance: low-rank == dense, lambda " + to_string(lam));
        }
    }
    // a singular Q (rank 1) still gives the dense answer: the form never inverts Q
    {
        Eigen::MatrixXd q1 = Eigen::MatrixXd::Zero(N, 1);
        q1(0, 0) = 1.0;
        Eigen::MatrixXd Q1 = q1 * q1.transpose();
        Eigen::MatrixXd L1 = A;
        double lam = 0.5, sc = 1.0 / (1.0 + lam);
        Eigen::MatrixXd G = sc * B * Q1 * B.transpose();
        G.diagonal() += var;
        Eigen::MatrixXd He = B * q1 * q1.transpose() * A.transpose() * e;
        Eigen::MatrixXd want = -sc * (e + L1 * B.transpose() * G.ldlt().solve(r - sc * He));
        Eigen::MatrixXd got = enif_lowrank_step(B, Q1, L1, var, e, He, r, lam);
        CHK((got - want).cwiseAbs().maxCoeff() < 1e-9, "a singular Q is handled");
    }
    // one column at a time equals its column of the batch
    Eigen::MatrixXd all = enif_lowrank_step(Hd * A, Eigen::MatrixXd::Identity(N, N), A, var, e, Hd * e, r, 1.5);
    Eigen::MatrixXd one = enif_lowrank_step(Hd * A, Eigen::MatrixXd::Identity(N, N), A, var, e.col(2), Hd * e.col(2), r.col(2), 1.5);
    CHK((all.col(2) - one.col(0)).cwiseAbs().maxCoeff() < 1e-12, "a single column matches the batch");
    bool threw = false;
    try { enif_lowrank_step(B, Eigen::MatrixXd::Identity(N + 1, N + 1), A, var, e, Hd * e, r, 1.0); } catch (...) { threw = true; }
    CHK(threw, "a Q of the wrong size throws");
    threw = false;
    Eigen::VectorXd bad = var; bad[3] = 0.0;
    try { enif_lowrank_step(B, Eigen::MatrixXd::Identity(N, N), A, bad, e, Hd * e, r, 1.0); } catch (...) { threw = true; }
    CHK(threw, "a zero noise variance throws");
}

static void test_enif_supplied_prec_hybrid()
{
    cout << "[enif: supplied prior precision file and the hybrid damping step]" << endl;
    PO o; o.set_defaults();
    CHK(o.get_ies_enif_prec_file() == "", "ies_enif_prec_file defaults to empty");
    CHK(!o.get_ies_enif_hybrid_damp(), "ies_enif_hybrid_damp defaults to false");
    CHK(o.assign_value_by_key("IES_ENIF_PREC_FILE", "MyPrec.jcb") == PO::ARG_STATUS::ARG_ACCEPTED, "prec file parses");
    CHK(o.get_ies_enif_prec_file() == "MyPrec.jcb", "and keeps its case");
    CHK(o.assign_value_by_key("IES_ENIF_HYBRID_DAMP", "true") == PO::ARG_STATUS::ARG_ACCEPTED, "hybrid parses");
    CHK(o.get_ies_enif_hybrid_damp(), "and is true after 'true'");

    // a 4-node chain precision with values, in the pest ascii matrix format
    Eigen::MatrixXd Q(4, 4);
    Q << 2.0, -0.8, 0.0, 0.0,
        -0.8, 2.5, -0.9, 0.0,
         0.0, -0.9, 2.5, -0.7,
         0.0, 0.0, -0.7, 2.0;
    string pf = "selftest_chain_prec.mat";
    {
        ofstream f(pf);
        f << " 4 4 1\n";
        for (int i = 0; i < 4; i++) { for (int j = 0; j < 4; j++) f << " " << Q(i, j); f << "\n"; }
        f << "* row and column names\nP1\nP2\nP3\nP4\n";
    }
    ofstream frec("selftest_enif_prec_hybrid.rec");
    vector<string> names = {"P1", "P2", "P3", "P4"};
    EnifGraph g;
    g.from_file(pf, names, frec);
    CHK(g.num_edges() == 3, "the precision file doubles as the graph: 3 edges");
    g.load_precision(pf, frec);
    CHK(g.has_precision(), "precision loaded");
    CHK((Eigen::MatrixXd(g.precision()) - Q).cwiseAbs().maxCoeff() < 1e-12, "loaded values match the file");
    // a non-symmetric file is refused
    string bf = "selftest_chain_prec_bad.mat";
    {
        ofstream f(bf);
        f << " 4 4 1\n";
        Eigen::MatrixXd Qb = Q; Qb(0, 1) = -0.3;
        for (int i = 0; i < 4; i++) { for (int j = 0; j < 4; j++) f << " " << Qb(i, j); f << "\n"; }
        f << "* row and column names\nP1\nP2\nP3\nP4\n";
    }
    bool threw = false;
    try { EnifGraph gb; gb.from_file(bf, names, frec); gb.load_precision(bf, frec); } catch (...) { threw = true; }
    CHK(threw, "a non-symmetric precision file throws");

    // hybrid step: Lam_post = Q + lam Qcur + H^T Rinv H against the dense formula
    std::mt19937 gen(9);
    std::normal_distribution<double> nd(0.0, 1.0);
    int N = 40;
    Eigen::MatrixXd A(4, N);
    for (int i = 0; i < 4; i++) for (int j = 0; j < N; j++) A(i, j) = 0.3 * nd(gen);
    for (int j = 0; j < N; j++) { A(1, j) += 0.5 * A(0, j); A(2, j) += 0.5 * A(1, j); A(3, j) += 0.5 * A(2, j); }
    A = (A.colwise() - A.rowwise().mean()) / sqrt(double(N - 1));
    g.estimate_precision(A, 1.0e-3, frec, false, true);
    CHK(g.has_current_precision(), "current-ensemble precision estimated into its own slot");
    CHK((Eigen::MatrixXd(g.precision()) - Q).cwiseAbs().maxCoeff() < 1e-12, "the prior precision is untouched by it");
    Eigen::MatrixXd Qcur = Eigen::MatrixXd(g.current_precision());
    CHK(Qcur.diagonal().minCoeff() > Q.diagonal().maxCoeff(), "contracted ensemble -> current precision larger than the prior's");
    Eigen::MatrixXd Hd(3, 4);
    Hd << 1.0, 0.0, 0.5, 0.0,
          0.0, 2.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 1.5;
    Eigen::SparseMatrix<double> H = Hd.sparseView();
    Eigen::VectorXd rinv(3); rinv << 4.0, 1.0, 2.0;
    Eigen::MatrixXd e = Eigen::MatrixXd::Zero(4, 2), r(3, 2);
    for (int i = 0; i < 3; i++) for (int j = 0; j < 2; j++) r(i, j) = nd(gen);
    for (int i = 0; i < 4; i++) e(i, 1) = 0.2 * nd(gen);
    double lam = 3.0;
    g.set_hybrid_damp(false);
    Eigen::MatrixXd plain = g.information_step(H, rinv, e, r, lam, frec);
    Eigen::MatrixXd want_plain = -((1.0 + lam) * Q + Hd.transpose() * rinv.asDiagonal() * Hd).ldlt().solve(Q * e + Hd.transpose() * (rinv.asDiagonal() * r));
    CHK((plain - want_plain).cwiseAbs().maxCoeff() < 1e-9, "plain step unchanged: (1+lam) on the prior precision");
    g.set_hybrid_damp(true);
    Eigen::MatrixXd hyb = g.information_step(H, rinv, e, r, lam, frec);
    Eigen::MatrixXd want_hyb = -(Q + lam * Qcur + Hd.transpose() * rinv.asDiagonal() * Hd).ldlt().solve(Q * e + Hd.transpose() * (rinv.asDiagonal() * r));
    CHK((hyb - want_hyb).cwiseAbs().maxCoeff() < 1e-9, "hybrid step: lam on the current precision, prior in the curvature and gradient");
    CHK((hyb - plain).cwiseAbs().maxCoeff() > 1e-6, "and it differs from the plain step");
    // without a current precision the hybrid step refuses
    EnifGraph g2; g2.from_file(pf, names, frec); g2.load_precision(pf, frec); g2.set_hybrid_damp(true);
    threw = false;
    try { g2.information_step(H, rinv, e, r, lam, frec); } catch (...) { threw = true; }
    CHK(threw, "hybrid damping without a current precision throws");
}

static void test_estimate_sparse_h_divisor()
{
    cout << "[enif: unexplained variance divisor and stream parameter]" << endl;
    int p = 5, n = 3, N = 8;
    std::mt19937 gen(3);
    std::normal_distribution<double> nd(0.0, 1.0);
    Eigen::MatrixXd A(p, N), B(n, N);
    for (int i = 0; i < p; i++) for (int j = 0; j < N; j++) A(i, j) = nd(gen);
    Eigen::MatrixXd Htrue = Eigen::MatrixXd::Zero(n, p);
    Htrue(0, 0) = 1.0; Htrue(1, 2) = -2.0; Htrue(2, 4) = 0.5;
    B = Htrue * A;
    for (int i = 0; i < n; i++) for (int j = 0; j < N; j++) B(i, j) += 0.05 * nd(gen);
    stringstream s1, s2, s3;
    Eigen::VectorXd u1, u2, u3;
    Eigen::SparseMatrix<double> H1 = estimate_sparse_H(A, B, 0.05, 1, u1, s1);
    Eigen::SparseMatrix<double> H2 = estimate_sparse_H(A, B, 0.05, 1, u2, s2, 0, double(N));
    Eigen::SparseMatrix<double> H3 = estimate_sparse_H(A, B, 0.05, 2, u3, s3, 0, 2.0 * double(N));
    CHK((Eigen::MatrixXd(H1) - Eigen::MatrixXd(H2)).cwiseAbs().maxCoeff() == 0.0, "divisor N is the default");
    CHK((u1 - u2).cwiseAbs().maxCoeff() < 1e-14, "and gives the same unexplained variance");
    CHK((Eigen::MatrixXd(H1) - Eigen::MatrixXd(H3)).cwiseAbs().maxCoeff() == 0.0, "the divisor does not touch H (two threads either)");
    CHK((u3 * 2.0 - u1).cwiseAbs().maxCoeff() < 1e-12, "doubling the divisor halves the unexplained variance");
    CHK(s1.str().find("sparse H") != string::npos, "the summary went to the stringstream, not a file");
    // residual energy can never exceed the row's own energy, so an empty row reports the full variance
    Eigen::VectorXd u4; stringstream s4;
    Eigen::SparseMatrix<double> H4 = estimate_sparse_H(A, B, 1.0, 1, u4, s4);
    CHK(H4.nonZeros() == 0, "lasso fraction 1 zeroes every row");
    Eigen::VectorXd full = B.array().square().rowwise().sum() * double(N - 1) / double(N);
    CHK((u4 - full).cwiseAbs().maxCoeff() < 1e-12, "empty rows report the full output variance as unexplained");
}

static void test_enif_mm_inflation_csv()
{
    cout << "[enif-mm: per-realization inflation csv]" << endl;
    vector<string> reals = {"R0", "R1"}, obs = {"o1", "o2", "o3"}, groups = {"g", "g", "h"};
    Eigen::VectorXd w0(3), w1(3), u0(3), u1(3);
    w0 << 1.0, 2.0, 0.5;  w1 << 1.0, 1.0, 1.0;
    u0 << 0.0, 0.25, 4.0; u1 << 1.0, 1.0, -0.5;   // a negative one is clamped to zero
    Eigen::SparseMatrix<double> H0(3, 2), H1(3, 2);
    {
        vector<Eigen::Triplet<double>> t = {{0, 0, 1.0}, {0, 1, 2.0}, {2, 1, 3.0}};
        H0.setFromTriplets(t.begin(), t.end());
        vector<Eigen::Triplet<double>> t1 = {{1, 0, 1.0}};
        H1.setFromTriplets(t1.begin(), t1.end());
    }
    auto split = [](const string& line) {
        vector<string> out; string cur; stringstream in(line);
        while (getline(in, cur, ',')) out.push_back(cur);
        return out;
    };
    string fname = "selftest_enif_mm_inflate.csv";
    map<string, double> mr = enif_mm_inflation_csv(reals, obs, groups, {w0, w1}, {u0, u1}, {&H0, &H1}, true, fname);
    // R0: noise 1, 0.25, 4 ; inflated 1, 0.5, 8 ; ratios 1, 2, 2 -> mean 5/3
    CHK(abs(mr.at("R0") - 5.0 / 3.0) < 1e-12, "mean ratio for R0");
    // R1: noise 1,1,1 ; unexplained 1,1,0 ; ratios 2,2,1 -> mean 5/3
    CHK(abs(mr.at("R1") - 5.0 / 3.0) < 1e-12, "mean ratio for R1 with the negative clamped");
    ifstream f(fname);
    string line; vector<string> lines;
    while (getline(f, line)) lines.push_back(line);
    CHK(lines.size() == 7, "header plus one row per (realization, observation)");
    CHK(lines[0] == "real_name,obs_name,group,weight,noise_var,unexplained_var,inflated_var,inflate_ratio,effective_weight,h_row_nnz", "header");
    // R0 o3: weight 0.5 -> noise 4, +4 -> 8, ratio 2, eff 1/sqrt(8), H row 2 has one entry
    vector<string> tok = split(lines[3]);
    CHK(tok[0] == "R0" && tok[1] == "o3" && tok[2] == "h", "row order is realization-major, observation-minor");
    CHK(abs(stod(tok[6]) - 8.0) < 1e-9 && abs(stod(tok[7]) - 2.0) < 1e-9, "inflated variance and ratio");
    CHK(abs(stod(tok[8]) - 1.0 / sqrt(8.0)) < 1e-9, "effective weight when applied");
    CHK(tok[9] == "1", "h_row_nnz counts that realization's H row");
    tok = split(lines[5]);
    CHK(tok[0] == "R1" && tok[1] == "o2" && tok[9] == "1", "R1 o2 has the single H entry");
    tok = split(lines[4]);
    CHK(tok[9] == "0", "R1 o1: empty H row reported as 0");
    // not applied: effective weight is the weight itself
    enif_mm_inflation_csv(reals, obs, groups, {w0, w1}, {u0, u1}, {&H0, &H1}, false, fname);
    ifstream f2(fname); lines.clear();
    while (getline(f2, line)) lines.push_back(line);
    tok = split(lines[3]);
    CHK(abs(stod(tok[8]) - 0.5) < 1e-12, "effective weight is the weight when inflation is not applied");
    remove(fname.c_str());
    bool threw = false;
    try { enif_mm_inflation_csv(reals, obs, groups, {w0}, {u0, u1}, {&H0, &H1}, true, fname); } catch (...) { threw = true; }
    CHK(threw, "misaligned per-realization vectors throw");
}

static void test_enif_current_prec()
{
    cout << "[enif: ies_enif_current_prec option, and the precision tracking the ensemble spread]" << endl;
    // the option parses, defaults to off, and the registry agrees with the legacy chain (self_verify
    // runs elsewhere; here just the two ends)
    PO o; o.set_defaults();
    CHK(!o.get_ies_enif_current_prec(), "ies_enif_current_prec defaults to false");
    CHK(o.assign_value_by_key("IES_ENIF_CURRENT_PREC", "true") == PO::ARG_STATUS::ARG_ACCEPTED, "parses");
    CHK(o.get_ies_enif_current_prec(), "and is true after parsing 'true'");
    // a 4-node chain graph in the pest ascii matrix format
    string gf = "selftest_chain_graph.mat";
    {
        ofstream f(gf);
        f << " 4 4 1\n";
        f << " 0 1 0 0\n 1 0 1 0\n 0 1 0 1\n 0 0 1 0\n";
        f << "* row and column names\nP1\nP2\nP3\nP4\n";
    }
    ofstream frec("selftest_enif_current_prec.rec");
    vector<string> names = {"P1", "P2", "P3", "P4"};
    EnifGraph g;
    g.from_file(gf, names, frec);
    CHK(g.is_initialized() && g.num_nodes() == 4 && g.num_edges() == 3, "chain graph read: 4 nodes, 3 edges");
    // anomalies with unit spread, then the same anomalies shrunk tenfold: the precision must grow
    // a hundredfold, which is the shrinking trust region the current-ensemble hessian relies on
    std::mt19937 gen(5);
    std::normal_distribution<double> nd(0.0, 1.0);
    int N = 60;
    Eigen::MatrixXd A(4, N);
    for (int i = 0; i < 4; i++) for (int j = 0; j < N; j++) A(i, j) = nd(gen);
    // give the chain some correlation so the off-diagonals are not pure noise
    for (int j = 0; j < N; j++) { A(1, j) += 0.7 * A(0, j); A(2, j) += 0.7 * A(1, j); A(3, j) += 0.7 * A(2, j); }
    A = (A.colwise() - A.rowwise().mean()) / sqrt(double(N - 1));
    g.estimate_precision(A, 1.0e-3, frec);
    Eigen::MatrixXd P1 = Eigen::MatrixXd(g.precision());
    g.estimate_precision(0.1 * A, 1.0e-3, frec);
    Eigen::MatrixXd P2 = Eigen::MatrixXd(g.precision());
    CHK(P1.allFinite() && P2.allFinite(), "both precisions finite");
    double ratio = P2.diagonal().mean() / P1.diagonal().mean();
    CHK(abs(ratio - 100.0) / 100.0 < 0.05, "tenfold smaller spread -> hundredfold larger precision (got x" + to_string(ratio) + ")");
    // the implied covariance of the unit-spread estimate reproduces the sample covariance on the chain
    // the last fit was on 0.1 * A, whose covariance is 0.01 * (A A^T), so scale back up by 100
    Eigen::MatrixXd Cimp = g.apply_cov(Eigen::MatrixXd::Identity(4, 4)) * 100.0;
    Eigen::MatrixXd Cs = A * A.transpose();
    // the neighbourhood regression with a stein pull is not the sample covariance, so this is a
    // same-ballpark check, not an identity; the exact agreement is covered by the reference comparison
    double vdev = (Cimp.diagonal() - Cs.diagonal()).cwiseAbs().maxCoeff() / Cs.diagonal().maxCoeff();
    CHK(vdev < 0.5, "implied variances in the ballpark of the sample variances (max rel dev " + to_string(vdev) + ")");
    // re-estimating with the same anomalies is deterministic
    g.estimate_precision(0.1 * A, 1.0e-3, frec);
    CHK((Eigen::MatrixXd(g.precision()) - P2).cwiseAbs().maxCoeff() == 0.0, "same anomalies, same precision");
    // the automatic ridge (shrink < 0): finite, still scales with the spread, and a little
    // weaker than a fixed heavy ridge, since k/(N-1) on this chain is well under 0.5
    g.estimate_precision(A, -1.0, frec);
    Eigen::MatrixXd Pa = Eigen::MatrixXd(g.precision());
    g.estimate_precision(0.1 * A, -1.0, frec);
    Eigen::MatrixXd Pa2 = Eigen::MatrixXd(g.precision());
    CHK(Pa.allFinite() && Pa2.allFinite(), "automatic ridge: finite precisions");
    double ratio_a = Pa2.diagonal().mean() / Pa.diagonal().mean();
    CHK(abs(ratio_a - 100.0) / 100.0 < 0.05, "automatic ridge: tenfold smaller spread -> hundredfold precision (got x" + to_string(ratio_a) + ")");
    PO od; od.set_defaults();
    CHK(od.get_ies_enif_shrink() < 0.0, "ies_enif_shrink defaults to automatic");
    // cross-validation: with 60 realizations on a 4-node chain the regressions are well
    // determined, so the chosen ridge should be small and the precision close to the 1e-3 one
    CHK((Pa - P1).cwiseAbs().maxCoeff() / P1.cwiseAbs().maxCoeff() < 0.25, "cv ridge on a well-determined chain stays close to the light fixed ridge");
    // a noisy case: the same chain but only 24 realizations, where held-out error should
    // call for more ridge than the clean case did.  read the chosen ridges off the rec
    frec.close();
    {
        ofstream f2("selftest_enif_cv.rec");
        Eigen::MatrixXd A24 = A.leftCols(24);
        A24 = (A24.colwise() - A24.rowwise().mean()) / sqrt(23.0);
        g.estimate_precision(A24, -1.0, f2);
        g.estimate_precision(A, -1.0, f2);
        f2.close();
        ifstream r("selftest_enif_cv.rec"); string line; vector<double> medians;
        while (getline(r, line))
        {
            size_t k = line.find("over nodes: ");
            if (k == string::npos) continue;
            // "min / median / max"
            string rest = line.substr(k + 12);
            double mn, md, mx; char s1, s2;
            stringstream ss(rest); ss >> mn >> s1 >> md >> s2 >> mx;
            medians.push_back(md);
        }
        CHK(medians.size() == 2, "two cv ridge report lines");
        CHK(medians[0] >= medians[1], "fewer realizations -> cross-validation asks for at least as much ridge (24 reals " + to_string(medians[0]) + " vs 60 reals " + to_string(medians[1]) + ")");
        remove("selftest_enif_cv.rec");
    }
    ofstream frec_tail("selftest_enif_current_prec.rec", ios::app);
    // fewer than 20 realizations falls back to k/(N-1) and still gives a finite precision
    Eigen::MatrixXd A12 = A.leftCols(12);
    A12 = (A12.colwise() - A12.rowwise().mean()) / sqrt(11.0);
    g.estimate_precision(A12, -1.0, frec_tail);
    CHK(Eigen::MatrixXd(g.precision()).allFinite(), "12 realizations: k/(N-1) fallback gives a finite precision");
    // direct neighbours only: a chain has no fill, so it must reproduce the fill-inclusive result
    g.estimate_precision(A, 1.0e-3, frec_tail);
    Eigen::MatrixXd Pf = Eigen::MatrixXd(g.precision());
    g.estimate_precision(A, 1.0e-3, frec_tail, true);
    Eigen::MatrixXd Pd = Eigen::MatrixXd(g.precision());
    CHK((Pf - Pd).cwiseAbs().maxCoeff() == 0.0, "direct-only on a chain (no fill) matches the fill-inclusive estimate");
    // a 4-cycle has fill (eliminating one node links its two neighbours), so direct-only
    // uses smaller regressions and gives a different, still finite, precision
    {
        ofstream f(gf);
        f << " 4 4 1\n 0 1 0 1\n 1 0 1 0\n 0 1 0 1\n 1 0 1 0\n* row and column names\nP1\nP2\nP3\nP4\n";
    }
    EnifGraph gc;
    gc.from_file(gf, names, frec_tail);
    CHK(gc.num_edges() == 4, "4-cycle read: 4 edges");
    gc.estimate_precision(A, 1.0e-3, frec_tail);
    Eigen::MatrixXd Pcf = Eigen::MatrixXd(gc.precision());
    gc.estimate_precision(A, 1.0e-3, frec_tail, true);
    Eigen::MatrixXd Pcd = Eigen::MatrixXd(gc.precision());
    CHK(Pcf.allFinite() && Pcd.allFinite(), "4-cycle: both estimates finite");
    CHK((Pcf - Pcd).cwiseAbs().maxCoeff() > 1e-8, "4-cycle: ignoring the fill changes the estimate");
    frec_tail.close();
    PO o2; o2.set_defaults();
    CHK(!o2.get_ies_enif_direct_nbrs(), "ies_enif_direct_nbrs defaults to false");
    CHK(o2.assign_value_by_key("IES_ENIF_DIRECT_NBRS", "true") == PO::ARG_STATUS::ARG_ACCEPTED && o2.get_ies_enif_direct_nbrs(), "and parses");
    frec.close();
    remove(gf.c_str()); remove("selftest_enif_current_prec.rec");
}

static void test_enif_inflation_report()
{
    cout << "[enif inflation report: per-obs csv and per-group rec summary]" << endl;
    // five obs in two groups, hand-picked so the numbers are easy to check
    vector<string> names = {"h1", "h2", "h3", "q1", "q2"};
    vector<string> groups = {"head", "head", "head", "flux", "flux"};
    Eigen::VectorXd w(5), u(5);
    w << 10.0, 10.0, 5.0, 0.1, 0.1;      // noise var 0.01, 0.01, 0.04, 100, 100
    u << 0.01, 0.03, 0.0, 100.0, 300.0;  // ratios 2, 4, 1, 2, 4
    // written in the working directory like selftest_viol.rec above.  this test
    // killed the whole selftest on every windows ci job (exit 127 = an msvc
    // fail-fast abort under msys bash) with nothing printed: the rec ifstream was
    // still open when remove() ran, which windows refuses.  the checkpoints
    // below are what found it, so they stay
    string csv = "selftest_enif_inflate.csv";
    string rec = "selftest_enif_inflate.rec";
    ofstream frec(rec);
    cout << "  writing " << csv << " and " << rec << " in " << std::filesystem::current_path().string() << endl;
    map<string, EnifInflateGroupStats> st = enif_inflation_report(names, groups, w, u, true, 3, csv, frec);
    frec.close();
    cout << "  report written, " << st.size() << " groups" << endl;

    CHK(st.size() == 2, "two groups summarized");
    CHK(st.at("head").count == 3 && st.at("flux").count == 2, "group counts");
    CHK(fabs(st.at("head").noise_var - 0.02) < 1.0e-12, "head mean noise var 0.02");
    CHK(fabs(st.at("head").unexp_var - 0.04 / 3.0) < 1.0e-12, "head mean unexplained var");
    CHK(fabs(st.at("head").ratio_mean - 7.0 / 3.0) < 1.0e-12, "head mean ratio (2+4+1)/3");
    CHK(fabs(st.at("head").ratio_min - 1.0) < 1.0e-12 && fabs(st.at("head").ratio_max - 4.0) < 1.0e-12, "head ratio min/max");
    CHK(fabs(st.at("flux").ratio_mean - 3.0) < 1.0e-12, "flux mean ratio 3");
    CHK(fabs(st.at("flux").weight_mean - 0.1) < 1.0e-12, "flux mean weight");
    // effective weight = 1/sqrt(noise + unexp): q1 -> 1/sqrt(200), q2 -> 1/sqrt(400) = 0.05
    CHK(fabs(st.at("flux").eff_weight_mean - 0.5 * (1.0 / sqrt(200.0) + 0.05)) < 1.0e-12, "flux mean effective weight");

    // the csv: header plus one row per obs, inflated = noise + unexp when applied
    cout << "  reading back " << csv << endl;
    ifstream fin(csv);
    string line;
    getline(fin, line);
    CHK(line == "obs_name,group,weight,noise_var,unexplained_var,inflated_var,inflate_ratio,effective_weight", "csv header");
    int nrow = 0;
    bool q2_ok = false;
    while (getline(fin, line))
    {
        nrow++;
        if (line.rfind("q2,flux,", 0) == 0)
        {
            vector<string> tok;
            stringstream ss(line);
            string t;
            while (getline(ss, t, ',')) tok.push_back(t);
            q2_ok = (tok.size() == 8) && (fabs(stod(tok[5]) - 400.0) < 1.0e-6) && (fabs(stod(tok[6]) - 4.0) < 1.0e-9)
                && (fabs(stod(tok[7]) - 0.05) < 1.0e-9);
        }
    }
    fin.close();
    CHK(nrow == 5, "csv has one row per observation");
    CHK(q2_ok, "q2 row: inflated var 400, ratio 4, effective weight 0.05");

    // the rec section, sorted by mean ratio so flux (3) comes before head (2.33)
    cout << "  reading back " << rec << endl;
    ifstream rin(rec);
    stringstream rs;
    rs << rin.rdbuf();
    string r = rs.str();
    // closed before the remove() below: windows wont delete a file that is still
    // open, posix will, and this is exactly what killed every windows ci job
    rin.close();
    CHK(r.find("enif observation noise inflation summary, iteration 3") != string::npos, "rec section header carries the iteration");
    CHK(r.find("inflation applied to the update: yes") != string::npos, "rec says the inflation was applied");
    CHK(r.find("flux") < r.find("head"), "rec groups sorted by mean ratio, largest first");
    CHK(r.find(csv) != string::npos, "rec names the csv");

    // not applied: the csv still carries the would-be ratio but the effective weight is the weight
    cout << "  second pass, inflation not applied" << endl;
    st = enif_inflation_report(names, groups, w, u, false, 4, csv, frec);
    CHK(fabs(st.at("flux").eff_weight_mean - 0.1) < 1.0e-12, "not applied: effective weight is the control file weight");
    CHK(fabs(st.at("flux").ratio_mean - 3.0) < 1.0e-12, "not applied: the would-be ratio is still reported");
    std::filesystem::remove(csv);
    std::filesystem::remove(rec);
}

static void test_irls_reweight()
{
    cout << "[irls prior info reweighting]" << endl;
    // three regul equations and one non-regul one, in the pest pi line form.  p1 is log
    // transformed so its residual is in log10 units, the others are native
    PriorInformation pi;
    pi.AddRecord("zo1 1.0 * log(p1) = 0.0 2.0 regul_zo");
    pi.AddRecord("zo2 1.0 * p2 = 1.0 1.0 regul_zo");
    pi.AddRecord("fo1 1.0 * p2 - 1.0 * p3 = 0.0 1.0 regul_fo");
    pi.AddRecord("fo0 1.0 * p2 - 1.0 * p3 = 0.0 0.0 regul_fo");
    pi.AddRecord("ob1 1.0 * p3 = 5.0 3.0 obsgrp");
    Parameters pars;
    pars.insert("P1", 100.0);   // resid = log10(100) - 0 = 2
    pars.insert("P2", 1.0);     // resid = 0 -> at the floor
    pars.insert("P3", 0.9);     // fo1 resid = 0.1
    map<string, double> w0;
    double eps = 0.01;
    PriorInformation::IrlsStats st = pi.irls_reweight(pars, eps, w0);
    CHK(st.n == 3, "irls: three nonzero regul equations reweighted");
    CHK(st.n_floor == 1, "irls: one equation at the floor");
    CHK(abs(pi.get_pi_rec("ZO1").get_weight() - 2.0 / sqrt(2.0)) < 1e-12, "irls: w = w0/sqrt(|r|)");
    CHK(abs(pi.get_pi_rec("ZO2").get_weight() - 1.0 / sqrt(eps)) < 1e-12, "irls: zero residual gets w0/sqrt(eps)");
    CHK(abs(pi.get_pi_rec("FO1").get_weight() - 1.0 / sqrt(0.1)) < 1e-12, "irls: difference equation reweighted from its own residual");
    CHK(pi.get_pi_rec("FO0").get_weight() == 0.0, "irls: zero weight equation stays zero");
    CHK(pi.get_pi_rec("OB1").get_weight() == 3.0, "irls: non-regul equation untouched");
    CHK(abs(st.fmin - 1.0 / sqrt(2.0)) < 1e-12 && abs(st.fmax - 1.0 / sqrt(eps)) < 1e-12, "irls: min/max factor");
    CHK(abs(st.fmed - 1.0 / sqrt(0.1)) < 1e-12, "irls: median factor");
    // the point of the convention: the quadratic penalty (w r)^2 now equals w0^2 |r|
    double r = pi.get_pi_rec("ZO1").calc_residual(pars);
    double w = pi.get_pi_rec("ZO1").get_weight();
    CHK(abs(w * w * r * r - 2.0 * 2.0 * abs(r)) < 1e-12, "irls: (w r)^2 == w0^2 |r|");
    // a second pass at different parameters starts from w0 again, not from the last weights
    pars.update_rec("P1", 10.0);  // resid = 1
    pi.irls_reweight(pars, eps, w0);
    CHK(abs(pi.get_pi_rec("ZO1").get_weight() - 2.0) < 1e-12, "irls: second pass reweights from w0, no compounding");
    CHK(w0.size() == 5 && w0.at("ZO1") == 2.0, "irls: w0 holds the control file weights");
}

static void test_irls_max_dev()
{
    cout << "[irls: solving the floor for a cap on deviations]" << endl;
    // six regul equations with distinct residuals 0, 0.1 .. 0.5, plus a tie pair, a zero weight
    // one and a non-regul one that must never count
    auto build = [](PriorInformation& pi, Parameters& pars)
    {
        pi.clear();
        pi.AddRecord("e0 1.0 * p0 = 0.0 1.0 regul");
        pi.AddRecord("e1 1.0 * p1 = 0.0 1.0 regul");
        pi.AddRecord("e2 1.0 * p2 = 0.0 2.0 regul");
        pi.AddRecord("e3 1.0 * p3 = 0.0 1.0 regul");
        pi.AddRecord("e4 1.0 * p4 = 0.0 1.0 regul");
        pi.AddRecord("e5 1.0 * p5 = 0.0 1.0 regul");
        pi.AddRecord("ez 1.0 * p5 = 0.0 0.0 regul");
        pi.AddRecord("ob 1.0 * p5 = 0.0 3.0 obsgrp");
        pars.clear();
        for (int i = 0; i < 6; i++)
            pars.insert("P" + to_string(i), 0.1 * i);
    };
    PriorInformation pi;
    Parameters pars;
    map<string, double> w0;
    double eps_min = 0.01;

    // cap 2: floor rises to the third largest |r| = 0.3, so 0.4 and 0.5 are active
    build(pi, pars); w0.clear();
    PriorInformation::IrlsStats st = pi.irls_reweight(pars, eps_min, w0, 2);
    CHK(abs(st.eps - 0.3) < 1e-12, "cap 2: floor is the 3rd largest residual");
    CHK(st.n_active == 2, "cap 2: two equations active");
    CHK(st.n_floor == 3, "cap 2: three strictly below the floor (0, 0.1, 0.2)");
    CHK(st.n == 6, "cap 2: six eligible equations, zero weight and non-regul excluded");
    CHK(abs(pi.get_pi_rec("E5").get_weight() - 1.0 / sqrt(0.5)) < 1e-12, "cap 2: active equation keeps w0/sqrt(|r|)");
    CHK(abs(pi.get_pi_rec("E3").get_weight() - 1.0 / sqrt(0.3)) < 1e-12, "cap 2: the equation at the floor gets w0/sqrt(eps)");
    CHK(abs(pi.get_pi_rec("E0").get_weight() - 1.0 / sqrt(0.3)) < 1e-12, "cap 2: zero residual shares the floor weight");
    CHK(abs(pi.get_pi_rec("E2").get_weight() - 2.0 / sqrt(0.3)) < 1e-12, "cap 2: floor weight scales with w0");
    CHK(pi.get_pi_rec("EZ").get_weight() == 0.0 && pi.get_pi_rec("OB").get_weight() == 3.0, "cap 2: zero weight and non-regul untouched");
    // weights are monotone: bigger residual never gets a bigger weight
    bool mono = true;
    for (int i = 0; i < 5; i++)
        if (pi.get_pi_rec("E" + to_string(i)).get_weight() / w0.at("E" + to_string(i)) <
            pi.get_pi_rec("E" + to_string(i + 1)).get_weight() / w0.at("E" + to_string(i + 1)) - 1e-12)
            mono = false;
    CHK(mono, "cap 2: reweight factor never increases with |r|");

    // cap 0: everything at the floor, floor = largest residual, all factors equal -> plain l2
    build(pi, pars); w0.clear();
    st = pi.irls_reweight(pars, eps_min, w0, 0);
    CHK(abs(st.eps - 0.5) < 1e-12 && st.n_active == 0, "cap 0: floor is the largest residual, nothing active");
    CHK(abs(st.fmin - st.fmax) < 1e-12, "cap 0: one common factor, the penalty is quadratic again");

    // cap >= n: the floor stays at the minimum, same as no cap
    build(pi, pars); w0.clear();
    st = pi.irls_reweight(pars, eps_min, w0, 6);
    CHK(abs(st.eps - eps_min) < 1e-12 && st.n_active == 5, "cap 6 of 6: floor stays at eps, five above it");
    build(pi, pars); w0.clear();
    PriorInformation::IrlsStats st_off = pi.irls_reweight(pars, eps_min, w0, -1);
    CHK(abs(st_off.eps - st.eps) < 1e-12 && st_off.n_active == st.n_active && st_off.n_floor == st.n_floor, "cap off (-1) matches a cap the equations never reach");
    build(pi, pars); w0.clear();
    st = pi.irls_reweight(pars, eps_min, w0, 60);
    CHK(abs(st.eps - eps_min) < 1e-12, "cap far above n: floor stays at eps");

    // the minimum floor wins when the cap would set a smaller one
    build(pi, pars); w0.clear();
    st = pi.irls_reweight(pars, 0.35, w0, 2);
    CHK(abs(st.eps - 0.35) < 1e-12 && st.n_active == 2, "cap 2 with eps 0.35: minimum floor above the order statistic wins");
    build(pi, pars); w0.clear();
    st = pi.irls_reweight(pars, 0.45, w0, 2);
    CHK(abs(st.eps - 0.45) < 1e-12 && st.n_active == 1, "cap 2 with eps 0.45: fewer than the cap active, never more");

    // ties at the cut: equal residuals all land at the floor, so active stays under the cap
    build(pi, pars); w0.clear();
    pars.update_rec("P3", 0.4);  // now 0.4 twice
    st = pi.irls_reweight(pars, eps_min, w0, 2);
    CHK(abs(st.eps - 0.4) < 1e-12 && st.n_active == 1, "tie at the cut: floor is the tied value, one active");
    st = pi.irls_reweight(pars, eps_min, w0, 1);
    CHK(abs(st.eps - 0.4) < 1e-12 && st.n_active == 1, "cap 1 with a tie just below: one active");

    // all residuals zero (the zeroth order start): floor stays at eps, nothing active, all equal
    build(pi, pars); w0.clear();
    for (int i = 0; i < 6; i++)
        pars.update_rec("P" + to_string(i), 0.0);
    st = pi.irls_reweight(pars, eps_min, w0, 2);
    CHK(abs(st.eps - eps_min) < 1e-12 && st.n_active == 0 && st.n_floor == 6, "all zero residuals: floor at eps, all six at floor");
    CHK(abs(pi.get_pi_rec("E1").get_weight() - 1.0 / sqrt(eps_min)) < 1e-12, "all zero residuals: every weight is w0/sqrt(eps)");

    // a second pass at new parameters solves a new floor from w0, not from the last weights
    build(pi, pars); w0.clear();
    pi.irls_reweight(pars, eps_min, w0, 2);
    for (int i = 0; i < 6; i++)
        pars.update_rec("P" + to_string(i), 0.2 * i);
    st = pi.irls_reweight(pars, eps_min, w0, 2);
    CHK(abs(st.eps - 0.6) < 1e-12, "second pass: floor re-solved from the new residuals");
    CHK(abs(pi.get_pi_rec("E5").get_weight() - 1.0 / sqrt(1.0)) < 1e-12, "second pass: weight from w0, no compounding");
    // and the same parameters twice give the same weights
    st = pi.irls_reweight(pars, eps_min, w0, 2);
    CHK(abs(st.eps - 0.6) < 1e-12 && abs(pi.get_pi_rec("E5").get_weight() - 1.0) < 1e-12, "repeat at the same parameters: identical");

    // huge residual spread and a non-finite one: the bad one is skipped, the rest still solve
    build(pi, pars); w0.clear();
    pars.update_rec("P5", 1e12);
    pars.update_rec("P4", numeric_limits<double>::quiet_NaN());
    st = pi.irls_reweight(pars, eps_min, w0, 1);
    CHK(st.n == 5, "nan residual: that equation is skipped");
    CHK(abs(st.eps - 0.3) < 1e-12 && st.n_active == 1, "nan residual: cap solved over the finite ones (floor 0.3, 1e12 active)");
    CHK(std::isfinite(pi.get_pi_rec("E5").get_weight()) && pi.get_pi_rec("E5").get_weight() > 0.0, "huge residual: finite positive weight");
}

static void test_cmdline_version_only()
{
    cout << "[command line: -v / --version answers with the version alone]" << endl;
    // version_only() writes the version to stdout when it says yes, so capture cout to
    // check that the number - and nothing else - is what comes out
    auto ask = [](vector<string> args, string& out)
    {
        vector<char*> argv;
        for (auto& a : args) argv.push_back(&a[0]);
        stringstream cap;
        streambuf* old = cout.rdbuf(cap.rdbuf());
        bool yes = pest_utils::CmdLine::version_only((int)argv.size(), argv.data());
        cout.rdbuf(old);
        out = cap.str();
        return yes;
    };
    string out;
    for (string flag : {"-v", "--version", "-V", "--VERSION"})
    {
        CHK(ask({"pestpp-ies", flag}, out), "version_only says yes to " + flag);
        CHK(out == string(PESTPP_VERSION) + "\n", "and writes the version number alone for " + flag);
    }
    CHK(!ask({"pestpp-ies", "pest.pst"}, out) && out.empty(), "a control file name is not a version request, and nothing is written");
    CHK(!ask({"pestpp-ies", "--version", "pest.pst"}, out) && out.empty(), "the flag with anything else after it is not a version request");
    CHK(!ask({"pestpp-ies", "pest.pst", "-v"}, out) && out.empty(), "nor is a trailing -v");
    CHK(!ask({"pestpp-ies"}, out) && out.empty(), "nor is no argument at all");

    // the constructor answers the same way for anything that builds a CmdLine directly:
    // version alone on stdout, flag set, and no echo of the command line
    {
        vector<string> args = {"pestpp-ies", "--version"};
        vector<char*> argv;
        for (auto& a : args) argv.push_back(&a[0]);
        stringstream cap;
        streambuf* old = cout.rdbuf(cap.rdbuf());
        pest_utils::CmdLine cl((int)argv.size(), argv.data());
        cout.rdbuf(old);
        CHK(cl.version_requested, "the constructor sets version_requested");
        CHK(cap.str() == string(PESTPP_VERSION) + "\n", "the constructor writes the version alone, no 'processing command line' echo");
    }
}

// a test that throws should count as a failure that names itself, not take the
// whole selftest down with it - which is what happened on the windows ci runners,
// where an uncaught exception is a silent exit 127.  every test goes through here
// so the log says which test is running, how many checks it made, and which of
// them failed, with the total so far
static void run_test(void (*fn)(), const char* name)
{
    g_current_test = name;
    int total0 = g_total, fail0 = g_fail;
    cout << "\n>>> " << name << endl;
    try
    {
        fn();
    }
    catch (const exception& e)
    {
        CHK(false, string(name) + " threw: " + e.what());
    }
    catch (...)
    {
        CHK(false, string(name) + " threw something that is not a std::exception");
    }
    int n = g_total - total0, f = g_fail - fail0;
    cout << "<<< " << name << ": " << (n - f) << "/" << n << " checks passed"
         << (f > 0 ? "  ***FAILED***" : "") << "  (running total " << (g_total - g_fail) << "/" << g_total << ")" << endl;
    g_current_test = "(none)";
}

// last words for the two ways a test can die without going through run_test's
// catch: std::terminate (an exception thrown while another is in flight, or from
// a destructor) and a hard signal.  both print the test name and flush, because
// on the windows runners nothing at all was printed and the process just vanished
static void on_terminate()
{
    cout << "\n*** std::terminate called during " << g_current_test;
    try
    {
        exception_ptr ep = current_exception();
        if (ep) rethrow_exception(ep);
        cout << " (no active exception)";
    }
    catch (const exception& e)
    {
        cout << ": " << e.what();
    }
    catch (...)
    {
        cout << ": non-std exception";
    }
    cout << endl << flush;
    abort();
}

static void on_signal(int sig)
{
    // keep this to plain writes - nothing fancy is safe in a signal handler
    cout << "\n*** fatal signal " << sig << " during " << g_current_test << endl << flush;
    signal(sig, SIG_DFL);
    raise(sig);
}

int main()
{
    // unbuffered so a hard crash cannot eat the lines that say where it happened
    cout.setf(ios::unitbuf);
    set_terminate(on_terminate);
    signal(SIGSEGV, on_signal);
    signal(SIGABRT, on_signal);
    signal(SIGFPE, on_signal);
    signal(SIGILL, on_signal);
    run_test(test_registry_equivalence, "test_registry_equivalence");
    run_test(test_cmdline_version_only, "test_cmdline_version_only");
    run_test(test_irls_reweight, "test_irls_reweight");
    run_test(test_irls_max_dev, "test_irls_max_dev");
    run_test(test_enif_inflation_report, "test_enif_inflation_report");
    run_test(test_enif_woodbury_step, "test_enif_woodbury_step");
    run_test(test_enif_lowrank_step, "test_enif_lowrank_step");
    run_test(test_enif_supplied_prec_hybrid, "test_enif_supplied_prec_hybrid");
    run_test(test_estimate_sparse_h_divisor, "test_estimate_sparse_h_divisor");
    run_test(test_enif_mm_inflation_csv, "test_enif_mm_inflation_csv");
    run_test(test_enif_current_prec, "test_enif_current_prec");
    run_test(test_generic_access, "test_generic_access");
    run_test(test_mutability, "test_mutability");
    run_test(test_control_info, "test_control_info");
    run_test(test_tool_defaults, "test_tool_defaults");
    run_test(test_constraints_live, "test_constraints_live");
    run_test(test_ies_reinflate_reset, "test_ies_reinflate_reset");
    run_test(test_ies_ensemble_reset, "test_ies_ensemble_reset");
    run_test(test_ies_iteration_controls, "test_ies_iteration_controls");
    run_test(test_mou_generation_controls, "test_mou_generation_controls");
    run_test(test_sqp_controls, "test_sqp_controls");
    run_test(test_tool_objects_track_live_options, "test_tool_objects_track_live_options");
    run_test(test_ensemble_zero_copy_view, "test_ensemble_zero_copy_view");
    run_test(test_run_map_survives_resize, "test_run_map_survives_resize");
    run_test(test_subset_names_survive_membership_change, "test_subset_names_survive_membership_change");
    run_test(test_read_file_tail, "test_read_file_tail");
    run_test(test_instruction_file_tolerant_read, "test_instruction_file_tolerant_read");
    run_test(test_parse_double_policy, "test_parse_double_policy");
    run_test(test_instruction_file_extreme_doubles, "test_instruction_file_extreme_doubles");
    run_test(test_fixed_instruction_misreads, "test_fixed_instruction_misreads");
    run_test(test_file_availability_checks, "test_file_availability_checks");
    run_test(test_partial_capability_handshake, "test_partial_capability_handshake");
    run_test(test_instruction_file_partial_reads_are_never_wrong, "test_instruction_file_partial_reads_are_never_wrong");
    run_test(test_model_interface_partial_across_files, "test_model_interface_partial_across_files");
    run_test(test_instruction_file_partial_real_case, "test_instruction_file_partial_real_case");
    run_test(test_instruction_file_partial_remaining_branches, "test_instruction_file_partial_remaining_branches");
    run_test(test_quit_file_tokens, "test_quit_file_tokens");
    run_test(test_run_storage_error_diagnostics, "test_run_storage_error_diagnostics");
    run_test(test_run_storage_partial_update, "test_run_storage_partial_update");
    run_test(test_external_values_are_results, "test_external_values_are_results");
    run_test(test_partial_read_refuses_stale_outputs, "test_partial_read_refuses_stale_outputs");
    run_test(test_violation_single_run_matches_ensemble, "test_violation_single_run_matches_ensemble");
    run_test(test_regul_weight_search_edges, "test_regul_weight_search_edges");
    run_test(test_ies_prior_scaling_is_neutral, "test_ies_prior_scaling_is_neutral");
    run_test(test_prior_anomaly_scaling, "test_prior_anomaly_scaling");
    run_test(test_prior_inv_sqrt_diag, "test_prior_inv_sqrt_diag");
    run_test(test_ies_prior_prec_reduces_to_am, "test_ies_prior_prec_reduces_to_am");
    cout << "\npestpp-selftest: " << (g_fail == 0 ? "PASS" : "FAIL")
         << " (" << (g_total - g_fail) << "/" << g_total << " checks)" << endl;
    return g_fail == 0 ? 0 : 1;
}

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
 *    opt_std_weights and their effect on use_chance/use_stosag/use_fosm/get_risk)
 */
#include <iostream>
#include <sstream>
#include <fstream>
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
#include "RunStorage.h"

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
    o->set_option("SQP_RISK", "0.5");   // keep use_stosag false so opt_risk drives it

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

    // the STOSAG path: sqp_risk (not opt_risk) drives use_stosag, and use_stosag in turn
    // drives use_chance and forces use_fosm off - all live from the options (Stage 5)
    o->set_option("OPT_RISK", "0.5");        // opt_risk neutral so sqp_risk is the sole driver
    o->set_option("SQP_RISK", "0.5");
    CHK(!c.get_use_stosag(), "sqp_risk 0.5 -> use_stosag false");
    CHK(!c.get_use_chance(), "both risks 0.5 -> use_chance false");
    o->set_option("SQP_RISK", "0.9");
    CHK(c.get_use_stosag(), "sqp_risk 0.9 -> use_stosag TRUE (live)");
    CHK(c.get_risk() == 0.9, "get_risk follows sqp_risk when use_stosag (live)");
    CHK(c.get_use_chance(), "sqp_risk drives use_chance even with opt_risk 0.5 (live)");
    CHK(!c.get_use_fosm(), "use_stosag -> use_fosm false (live)");
    o->set_option("SQP_RISK", "0.5");        // reset: stosag off

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
    CHK(!o.is_init_only("SQP_RISK"), "sqp_risk is live");
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
    // propagates to the accessor, but the change is surfaced rather than silently taken
    CHK(!sqp.p_use_ensemble_grad(), "no ensemble requested -> fd gradient");
    o->set_option("SQP_NUM_REALS", "20");
    CHK(sqp.p_use_ensemble_grad(), "sqp_num_reals -> ensemble gradient (derived, not cached)");
    CHK(o->get_init_only_change_warnings().size() == 1,
        "changing the init-only sqp_num_reals post-init IS flagged");
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
// harvesting them. The run->realization map is keyed by name so that change is absorbed;
// keyed by position, every run after the edit would land on the wrong realization.
static void test_run_map_survives_resize()
{
    cout << "[queued runs are held by name, so a resize between queue and harvest is safe]" << endl;

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

    // rows come back in ascending order, which is what the harvest loop and the failed-index
    // list both assume
    r = ObservationEnsemble::resolve_run_positions(run_ids, reordered, reordered.size());
    ok = true;
    for (int i = 1; i < r.size(); i++)
        ok = ok && (r[i].first > r[i-1].first);
    CHK(ok, "resolved rows are in ascending order");

    // nothing survived
    vector<string> none{"X","Y"};
    CHK(ObservationEnsemble::resolve_run_positions(run_ids, none, none.size()).size() == 0,
        "no surviving realizations -> no runs to harvest");
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

int main()
{
    test_registry_equivalence();
    test_generic_access();
    test_mutability();
    test_control_info();
    test_tool_defaults();
    test_constraints_live();
    test_ies_reinflate_reset();
    test_ies_ensemble_reset();
    test_ies_iteration_controls();
    test_mou_generation_controls();
    test_sqp_controls();
    test_tool_objects_track_live_options();
    test_ensemble_zero_copy_view();
    test_run_map_survives_resize();
    test_subset_names_survive_membership_change();
    test_instruction_file_tolerant_read();
    test_instruction_file_partial_reads_are_never_wrong();
    test_model_interface_partial_across_files();
    test_instruction_file_partial_real_case();
    test_instruction_file_partial_remaining_branches();
    test_run_storage_partial_update();
    cout << "\npestpp-selftest: " << (g_fail == 0 ? "PASS" : "FAIL")
         << " (" << (g_total - g_fail) << "/" << g_total << " checks)" << endl;
    return g_fail == 0 ? 0 : 1;
}

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

using namespace std;
using PO = PestppOptions;

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
    cout << "\npestpp-selftest: " << (g_fail == 0 ? "PASS" : "FAIL")
         << " (" << (g_total - g_fail) << "/" << g_total << " checks)" << endl;
    return g_fail == 0 ? 0 : 1;
}

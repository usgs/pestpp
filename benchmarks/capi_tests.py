"""Tests that drive PEST++ through the C ABI, the way a python caller would.

Everything else in CI drives the *executables*, which only ever run the built-in loop in the
built-in order. These tests are the only ones that exercise the shared library: loading it,
creating a handle, and running an iteration the caller controls -- queue the runs, drive the
run manager, change the ensemble mid-flight, harvest.

That middle step is the point. The built-in loop never changes ensemble membership between
queueing runs and harvesting them, so no benchmark can catch a regression in the name-keyed
run map; only a caller doing something the built-in loop never does can.

The library is located in the build tree rather than installed, matching how CI builds it.
"""
import os
import platform
import shutil
import subprocess
import sys
import time

import numpy as np
import pyemu

# anchored to this file, not the cwd, so the thin ctypes layer that ships with this checkout
# is the one under test even if the runner is invoked from elsewhere
_BENCH = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_BENCH)
sys.path.insert(0, os.path.join(_REPO, "python"))
from pestpp_lib import (  # noqa: E402
    PestppLib, PestppError, PESTPP_OK, PAR_EN, OBS_EN, NOISE_EN, WEIGHTS_EN,
    RM_SERIAL, RM_PANTHER, RM_EXTERNAL,
    TOOL_IES, TOOL_DA, TOOL_MOU, TOOL_SQP, WORKER_COMPLETED,
)

plat = "unknown"
if "linux" in platform.platform().lower():
    plat = "linux"
elif "darwin" in platform.platform().lower() or "macos" in platform.platform().lower():
    plat = "apple"
else:
    plat = "windows"

exe = ".exe" if plat == "windows" else ""
_sub = {"windows": "win", "apple": "mac", "linux": "linux"}[plat]

# The model binary (mfnwt) has to be findable by the forward run, and this has to be an
# ABSOLUTE path to the platform subdirectory. CI puts "../../test_bin/<plat>" on PATH, which
# is relative and therefore resolves against each process's cwd - it only works for case
# directories exactly two levels below benchmarks/, which is what every other benchmark script
# happens to use. The serial run manager runs the model from the session working directory, so
# a case directory at any other depth silently gets "execv() failed for command: mfnwt".
os.environ["PATH"] += os.pathsep + os.path.join(_BENCH, "test_bin", _sub)


def _find_agent_exe():
    """pestpp-ies, used to start PANTHER workers.

    The other benchmark scripts hardcode one relative path that only resolves under the CI
    directory layout. Try the candidates instead, so this file also runs from a plain clone.
    """
    roots = [os.path.join("pestpp", "bin"),
             os.path.join("..", "..", "pestpp", "bin"),
             os.path.join("..", "..", "..", "..", "pestpp", "bin"),
             os.path.join(_REPO, "bin")]
    for root in roots:
        for cand in (os.path.join(root, _sub, "pestpp-ies" + exe),
                     os.path.join(root, "pestpp-ies" + exe)):
            if os.path.exists(cand):
                return os.path.abspath(cand)
    raise RuntimeError("could not find pestpp-ies under any of {0}".format(roots))


port = 4062


def _find_library():
    """Locate the built C ABI shared library.

    Searched rather than hardcoded because the path depends on the generator: single-config
    builds drop it beside the sources, MSVC adds a per-config subdirectory.
    """
    names = {"windows": "pestpp-api.dll",
             "apple": "pestpp-api.dylib",
             "linux": "pestpp-api.so"}
    want = names[plat]
    roots = [os.path.join(_REPO, "build"), os.path.join("..", "..", "pestpp", "build"),
             os.path.join("..", "..", "..", "..", "pestpp", "build")]
    found = []
    for root in roots:
        if not os.path.isdir(root):
            continue
        for dirpath, dirnames, filenames in os.walk(root):
            # skip CPack staging: it holds a full copy of the install tree, so a stale one
            # from an earlier build looks exactly like the freshly compiled library
            dirnames[:] = [d for d in dirnames if d != "_CPack_Packages"]
            if want in filenames:
                found.append(os.path.abspath(os.path.join(dirpath, want)))
    if found:
        return max(found, key=os.path.getmtime)      # newest wins
    raise RuntimeError(
        "could not find {0} under any of {1}. The C ABI is built by the normal cmake build "
        "(always SHARED); if this fires, the build did not produce it.".format(want, roots))


# ---- case setup ---------------------------------------------------------------------------

def _copy_base(test_d, base_d=os.path.join("ies_10par_xsec", "template")):
    """Copy a checked-in base case into a fresh working directory.

    Both paths are anchored to this file's directory rather than the process cwd. The C ABI
    legitimately changes the working directory while a session is open (see ScopedWorkingDir),
    so anything cwd-relative here is one restore failure away from resolving somewhere else -
    which surfaces as a base case that has mysteriously gone missing.

    Only ever copy from a directory that is CHECKED IN. Several case directories under
    benchmarks/ are created dynamically by other benchmark scripts, so they exist in a working
    tree that has run those scripts and not in a fresh CI checkout.
    """
    base_abs = base_d if os.path.isabs(base_d) else os.path.join(_BENCH, base_d)
    test_abs = test_d if os.path.isabs(test_d) else os.path.join(_BENCH, test_d)
    if not os.path.isdir(base_abs):
        raise RuntimeError(
            "base case '{0}' not found (cwd is '{1}'). It must be a checked-in directory - "
            "directories created dynamically by other benchmark scripts are not available "
            "here.".format(base_abs, os.getcwd()))
    if os.path.exists(test_abs):
        shutil.rmtree(test_abs)
    shutil.copytree(base_abs, test_abs)
    return test_abs


def _setup(test_d, noptmax=1, num_reals=6):
    """A small, fast ies case in its own directory."""
    test_d = _copy_base(test_d)
    pst = pyemu.Pst(os.path.join(test_d, "pest.pst"))
    pst.control_data.noptmax = noptmax
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["random_seed"] = 11
    pst.write(os.path.join(test_d, "pest.pst"), version=2)
    return test_d


def _setup_da(test_d, noptmax=1, num_reals=6):
    """A single-cycle da case.

    Every quantity is marked cycle -1 ("all cycles"), which makes this the batch form. The
    C ABI exposes one cycle only - see the DaAdapter comment in pestpp_capi.cpp - so a
    multi-cycle case would not be driveable through the API yet.
    """
    test_d = _copy_base(test_d)
    pst = pyemu.Pst(os.path.join(test_d, "pest.pst"))
    for df in (pst.parameter_data, pst.observation_data,
               pst.model_input_data, pst.model_output_data):
        df.loc[:, "cycle"] = -1
    pst.control_data.noptmax = noptmax
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["da_num_reals"] = num_reals
    pst.pestpp_options["random_seed"] = 11
    pst.write(os.path.join(test_d, "pest.pst"), version=2)
    return test_d


def _setup_sqp(test_d, noptmax=1, num_reals=8):
    """An sqp case: one observation becomes the objective, another a constraint.

    sqp refuses to initialize without a recognized constraint/objective group, so the group
    names matter here - 'l_' marks a less-than constraint.
    """
    test_d = _copy_base(test_d)
    pst = pyemu.Pst(os.path.join(test_d, "pest.pst"))
    obs = pst.observation_data
    obj, con = pst.nnz_obs_names[0], pst.nnz_obs_names[1]
    obs.loc[obj, "obgnme"] = "obj_fn"
    obs.loc[con, "obgnme"] = "l_head"
    obs.loc[con, "obsval"] = float(obs.loc[con, "obsval"]) * 1.5
    pst.pestpp_options["opt_obj_func"] = obj
    # a real dv ensemble rather than a single point
    pst.pestpp_options["sqp_num_reals"] = num_reals
    pst.pestpp_options["random_seed"] = 11
    pst.control_data.noptmax = noptmax
    pst.write(os.path.join(test_d, "pest.pst"), version=2)
    return test_d


def _setup_mou(test_d, noptmax=1, pop_size=10):
    """A small mou case built from g07.

    g07 rather than the mou benchmarks' constr case, because constr_template is created
    dynamically by mou_tests.py and so does not exist in this job's checkout. g07 is checked
    in, already carries obj_fn and l_constraint groups, and its forward run is pure python -
    which also makes this the one tool test that needs no model binary at all.
    """
    test_d = _copy_base(test_d, os.path.join("g07", "template"))
    pst = pyemu.Pst(os.path.join(test_d, "g07.pst"))
    # g07 ships configured for sqp; drop those so the options in play are mou's
    for k in [k for k in pst.pestpp_options if k.startswith("sqp_")]:
        pst.pestpp_options.pop(k)
    pst.pestpp_options["mou_population_size"] = pop_size
    pst.pestpp_options["mou_generator"] = "de"
    pst.pestpp_options["random_seed"] = 11
    pst.control_data.noptmax = noptmax
    pst.write(os.path.join(test_d, "g07.pst"), version=2)
    return test_d


def _drive_batch(ies):
    """Queue, run and harvest, with the caller owning the run loop."""
    n_queued = ies.queue_runs()
    ies.begin_batch()
    while not ies.run_slice(0.05):
        pass
    ies.end_batch()
    return n_queued


# ---- ies ----------------------------------------------------------------------------------

def capi_smoke_test():
    """The library loads, a handle runs iterations, and views are live windows on the data."""
    wd = _setup("capi_smoke", noptmax=2)
    with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd) as ies:
        assert ies.get_option("IES_NUM_REALS") == "6", ies.get_option("IES_NUM_REALS")

        ies.initialize()
        assert ies.get_iteration() == 0, ies.get_iteration()

        pe, pe_tok = ies.get_ensemble_view(PAR_EN)
        rnames = ies.get_ensemble_row_names(PAR_EN)
        cnames = ies.get_ensemble_col_names(PAR_EN)
        assert pe.shape == (len(rnames), len(cnames)), (pe.shape, len(rnames), len(cnames))
        # Eigen is column-major; a view that is not F-contiguous means the numpy wrapper is
        # transposing silently, which would corrupt every value a caller writes back
        assert pe.flags["F_CONTIGUOUS"], "ensemble view lost its column-major layout"

        # it is a window, not a copy: a write through numpy must be visible to C++
        before = pe[0, 0]
        pe[0, 0] = before + 1234.5
        assert ies.get_ensemble_view(PAR_EN)[0][0, 0] == before + 1234.5, "view is a copy, not a view"
        pe[0, 0] = before

        for _ in range(2):
            ies.solve_iteration()
        assert ies.get_iteration() == 2, ies.get_iteration()
        oe, _ = ies.get_ensemble_view(OBS_EN)
        assert oe.shape[0] > 0 and np.isfinite(oe).all(), oe.shape
        ies.finalize()


def capi_snapshot_roundtrip_test():
    """The CTL-space snapshot round-trips exactly, and is matched by name rather than order."""
    wd = _setup("capi_snapshot")
    with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd) as ies:
        ies.initialize()
        vals, rows, cols = ies.get_par_snapshot()
        assert vals.shape == (len(rows), len(cols)), (vals.shape, len(rows), len(cols))

        ies.set_par_snapshot(vals, rows, cols)
        again, _, _ = ies.get_par_snapshot()
        assert np.array_equal(vals, again), "snapshot round trip changed values"

        # feed the rows back reversed: matching is by name, so the values must land in the
        # same realizations they came from
        perm = list(reversed(range(len(rows))))
        ies.set_par_snapshot(vals[perm, :], [rows[i] for i in perm], cols)
        after, _, _ = ies.get_par_snapshot()
        assert np.array_equal(vals, after), "snapshot was matched by position, not by name"


def capi_resize_between_queue_and_harvest_test():
    """The case no benchmark can reach: change membership while runs are in flight.

    Evaluate once with no membership change to learn what each realization's results are.
    Then evaluate the *same* parameters again, but drop a realization after the runs are made
    and before harvesting, and require every survivor to receive the identical results.

    That equality is the real assertion. Runs are tracked by realization name, so results
    follow their realizations to new row positions. Were they tracked by position, every
    realization after the dropped one would silently receive its neighbour's results - and
    the values would still be finite and nonzero, which is why a structural check is not
    enough to catch it.
    """
    wd = _setup("capi_resize")
    with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd) as ies:
        ies.initialize()
        names = ies.get_ensemble_row_names(PAR_EN)
        assert len(names) >= 4, "need a few realizations for this to mean anything"

        # -- reference pass: no membership change --
        assert _drive_batch(ies) == len(names)
        assert ies.process_runs() == 0, "reference pass had failed runs"
        oe_names = ies.get_ensemble_row_names(OBS_EN)
        ref_oe, _ = ies.get_ensemble_view(OBS_EN)
        reference = {n: ref_oe[i].copy() for i, n in enumerate(oe_names)}

        # the comparison below is only meaningful if realizations are actually distinguishable
        # - if they all produced the same results, misattribution would pass unnoticed
        closest = min(np.abs(ref_oe[i] - ref_oe[j]).max()
                      for i in range(len(oe_names)) for j in range(i + 1, len(oe_names)))
        assert closest > 1.0e-6, (
            "realizations are indistinguishable (closest pair differs by {0}), so this test "
            "cannot detect misattributed runs".format(closest))

        # -- second pass: same parameters, but resize between queue and harvest --
        # harvest does not touch pe when nothing fails, so these are the same runs
        assert _drive_batch(ies) == len(names)

        victim_par = names[1]
        victim_obs = oe_names[1]          # par and obs realizations pair by position
        ies.drop_realizations([victim_par])

        left = ies.get_ensemble_row_names(PAR_EN)
        assert len(left) == len(names) - 1, (len(left), len(names))
        assert victim_par not in left, "dropped realization survived in the par ensemble"
        assert len(ies.get_ensemble_row_names(OBS_EN)) == len(left), \
            "obs ensemble was not dropped in step with the par ensemble"

        n_failed = ies.process_runs()
        assert n_failed == 0, "{0} runs failed unexpectedly".format(n_failed)

        oe, _ = ies.get_ensemble_view(OBS_EN)
        surviving = ies.get_ensemble_row_names(OBS_EN)
        assert oe.shape[0] == len(left), (oe.shape, len(left))
        assert victim_obs not in surviving, "dropped realization survived in the obs ensemble"

        for i, name in enumerate(surviving):
            assert np.allclose(oe[i], reference[name]), (
                "realization '{0}' (now row {1}) did not get its own results back after the "
                "resize - runs are being mapped by position, not by name".format(name, i))


def capi_error_reporting_test():
    """Misuse produces a clear message rather than a crash across the C boundary."""
    wd = _setup("capi_errors")
    with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd) as ies:
        ies.initialize()

        try:
            ies.process_runs()
            raise AssertionError("harvesting with nothing queued should raise")
        except PestppError as e:
            assert "no queued runs" in str(e), str(e)

        try:
            ies.drop_realizations(["NOT_A_REALIZATION"])
            raise AssertionError("dropping an unknown realization should raise")
        except PestppError as e:
            assert "no such realization" in str(e), str(e)

        # the run-manager queries are panther-only, and say so
        assert not ies.supports_live_control(), "serial manager claimed live control"
        try:
            ies.get_run_time_stats()
            raise AssertionError("run stats on the serial manager should raise")
        except PestppError as e:
            assert "PANTHER" in str(e), str(e)

        try:
            ies.set_option("NOT_AN_OPTION", "1")
            raise AssertionError("an unknown option should raise")
        except PestppError as e:
            assert "unknown option" in str(e), str(e)

        # queueing twice without harvesting is a caller bug, not a silent overwrite
        ies.queue_runs()
        try:
            ies.queue_runs()
            raise AssertionError("double queue should raise")
        except PestppError as e:
            assert "already queued" in str(e), str(e)


def capi_caller_owned_initial_batch_test():
    """The caller owns the prior-ensemble evaluation - and can replace it before it runs.

    Before initialize() was split, it drew the prior ensemble and evaluated it in one
    uninterruptible call, so the realizations that got run were always the ones pest++ drew.
    Here the caller substitutes its own parameter values in the window between prepare and
    queue, and the assertion is that *those* are what the model saw.

    Proving it: every realization is set to identical parameter values, so every harvested
    observation row must be identical too. With the drawn ensemble they are emphatically not
    - the closest pair differs by ~1 in this case - so this cannot pass by accident.
    """
    wd = _setup("capi_init_split")
    with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd) as ies:
        n = ies.initialize_prepare()
        names = ies.get_ensemble_row_names(PAR_EN)
        assert n == len(names), (n, len(names))
        assert n > 1, "need several realizations for this to mean anything"

        # the window that did not exist before: replace the drawn prior with our own
        vals, rows, cols = ies.get_par_snapshot()
        mine = vals.copy()
        mine[:, :] = vals[0, :]          # every realization identical to the first
        ies.set_par_snapshot(mine, rows, cols)

        assert ies.queue_runs() == n
        ies.begin_batch()
        while not ies.run_slice(0.05):
            pass
        ies.end_batch()
        assert ies.process_runs() == 0, "runs failed on the substituted ensemble"
        ies.initialize_finish()

        # our values are what survived
        after, _, _ = ies.get_par_snapshot()
        assert np.allclose(after, mine), "the substituted parameter values were not kept"

        # and they are what the model saw
        oe, _ = ies.get_ensemble_view(OBS_EN)
        assert oe.shape[0] == n, (oe.shape, n)
        spread = np.abs(oe - oe[0]).max()
        assert spread < 1.0e-8, (
            "identical parameters produced different observations (max spread {0}) - the "
            "substituted ensemble was not the one evaluated".format(spread))

        # and the tool is properly initialized: it can still take a step
        ies.solve_iteration()
        assert ies.get_iteration() == 1, ies.get_iteration()
        ies.finalize()


def capi_initialize_split_guardrails_test():
    """A half-initialized tool refuses to be stepped, and the two halves must be paired."""
    wd = _setup("capi_init_guard")
    with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd) as ies:
        try:
            ies.initialize_finish()
            raise AssertionError("finish without prepare should raise")
        except PestppError as e:
            assert "nothing to finish" in str(e), str(e)

        ies.initialize_prepare()
        try:
            ies.initialize_prepare()
            raise AssertionError("a second prepare should raise")
        except PestppError as e:
            assert "already in progress" in str(e), str(e)

        # stepping a half-initialized tool would use an unevaluated ensemble
        try:
            ies.solve_iteration()
            raise AssertionError("solve_iteration mid-initialization should raise")
        except PestppError as e:
            assert "initialization is incomplete" in str(e), str(e)


def capi_initialize_prepare_reports_zero_test():
    """Tools that hand over no initial batch say so, so callers need no special-casing.

    mou evaluates several populations inside initialize() and sqp's only batch there is a
    single control-file-values run, so neither exposes a caller-owned batch. Both must report
    0 rather than leaving a caller to discover it.
    """
    for tool, wd, pst_name, tag in (
            (TOOL_MOU, _setup_mou("capi_init_mou"), "g07.pst", "mou"),
            (TOOL_SQP, _setup_sqp("capi_init_sqp"), "pest.pst", "sqp")):
        with PestppLib(_find_library(), tool, pst_name, wd) as t:
            assert t.initialize_prepare() == 0, "{0} should report no caller-owned batch".format(tag)
            t.initialize_finish()
            # still fully initialized: it can take a step
            t.solve_iteration()
            assert t.get_iteration() > 0, tag


# ---- da, mou, sqp -------------------------------------------------------------------------

def _drive_tool(tool, wd, pst_name, tag):
    """The caller-owned loop, for any tool: initialize, queue, run, harvest, advance.

    Every tool is driven through the same entry points. That is the point of the API - the
    differences between ies, da, mou and sqp live behind pestpp_solve_iteration(), not in
    the shape of the calls a caller has to make.
    """
    with PestppLib(_find_library(), tool, pst_name, wd) as t:
        t.initialize()

        pe, _ = t.get_ensemble_view(PAR_EN)
        oe, _ = t.get_ensemble_view(OBS_EN)
        pnames = t.get_ensemble_row_names(PAR_EN)
        assert pe.shape[0] > 0, "{0}: empty parameter ensemble".format(tag)
        assert pe.shape == (len(pnames), len(t.get_ensemble_col_names(PAR_EN))), \
            "{0}: par view {1} disagrees with its name lists".format(tag, pe.shape)
        assert pe.flags["F_CONTIGUOUS"], "{0}: par view is not column-major".format(tag)
        assert oe.shape[0] == pe.shape[0], \
            "{0}: par/obs ensembles disagree on size ({1} vs {2})".format(tag, pe.shape, oe.shape)

        # the caller owns the run loop
        n_queued = t.queue_runs()
        assert n_queued == pe.shape[0], (tag, n_queued, pe.shape[0])
        t.begin_batch()
        while not t.run_slice(0.05):
            pass
        t.end_batch()
        n_failed = t.process_runs()
        assert n_failed == 0, "{0}: {1} runs failed".format(tag, n_failed)

        oe, _ = t.get_ensemble_view(OBS_EN)
        assert np.isfinite(oe).all(), "{0}: harvested obs ensemble has non-finite values".format(tag)

        # and one algorithm step through the same entry point every tool uses
        before = t.get_iteration()
        t.solve_iteration()
        assert t.get_iteration() > before, \
            "{0}: iteration counter did not advance ({1} -> {2})".format(
                tag, before, t.get_iteration())
        t.finalize()


def capi_da_test():
    """da, single cycle, driven through the caller-owned loop."""
    _drive_tool(TOOL_DA, _setup_da("capi_da"), "pest.pst", "da")


def capi_sqp_test():
    """sqp, driven through the caller-owned loop."""
    _drive_tool(TOOL_SQP, _setup_sqp("capi_sqp"), "pest.pst", "sqp")


def capi_mou_test():
    """mou, driven through the caller-owned loop - a generation rather than an iteration."""
    _drive_tool(TOOL_MOU, _setup_mou("capi_mou"), "g07.pst", "mou")


def capi_da_resize_between_queue_and_harvest_test():
    """The resize guard, for da.

    da shares EnsembleMethod's queue/harvest with ies, so the same misattribution is possible
    and the same reference comparison catches it.
    """
    wd = _setup_da("capi_da_resize")
    with PestppLib(_find_library(), TOOL_DA, "pest.pst", wd) as t:
        t.initialize()
        names = t.get_ensemble_row_names(PAR_EN)
        assert len(names) >= 4, "need a few realizations for this to mean anything"

        assert _drive_batch(t) == len(names)
        assert t.process_runs() == 0, "da reference pass had failed runs"
        oe_names = t.get_ensemble_row_names(OBS_EN)
        ref_oe, _ = t.get_ensemble_view(OBS_EN)
        reference = {n: ref_oe[i].copy() for i, n in enumerate(oe_names)}
        closest = min(np.abs(ref_oe[i] - ref_oe[j]).max()
                      for i in range(len(oe_names)) for j in range(i + 1, len(oe_names)))
        assert closest > 1.0e-6, \
            "da realizations are indistinguishable, so this test proves nothing"

        assert _drive_batch(t) == len(names)
        victim = names[1]
        t.drop_realizations([victim])
        assert victim not in t.get_ensemble_row_names(PAR_EN)
        assert t.process_runs() == 0

        oe, _ = t.get_ensemble_view(OBS_EN)
        surviving = t.get_ensemble_row_names(OBS_EN)
        for i, name in enumerate(surviving):
            assert np.allclose(oe[i], reference[name]), (
                "da realization '{0}' did not get its own results back after the resize - "
                "runs are being mapped by position, not by name".format(name))


def capi_tool_ensemble_availability_test():
    """Each tool reports only the ensembles it actually has, rather than inventing them.

    ies and da carry a noise ensemble and a weights ensemble; mou and sqp do not. Asking for
    one that does not exist must say so, not hand back an empty array that reads as real.
    """
    cases = [
        (TOOL_MOU, _setup_mou("capi_avail_mou"), "g07.pst", "mou"),
        (TOOL_SQP, _setup_sqp("capi_avail_sqp"), "pest.pst", "sqp"),
    ]
    for tool, wd, pst_name, tag in cases:
        with PestppLib(_find_library(), tool, pst_name, wd) as t:
            t.initialize()
            # the two it does have
            assert t.get_ensemble_view(PAR_EN)[0].shape[0] > 0, tag
            assert t.get_ensemble_view(OBS_EN)[0].shape[0] > 0, tag
            for missing in (NOISE_EN, WEIGHTS_EN):
                try:
                    t.get_ensemble_view(missing)
                    raise AssertionError(
                        "{0} should not report an ensemble id {1}".format(tag, missing))
                except PestppError as e:
                    assert "has no ensemble" in str(e), str(e)


def capi_run_manager_selection_test():
    """All three run managers can be asked for, and each reports what it can do.

    Serial and external finish a batch in one slice and cannot be introspected; only PANTHER
    can. A caller has to be able to find that out rather than discover it from an error, which
    is what supports_live_control() is for.
    """
    for rm, name in ((RM_SERIAL, "serial"), (RM_EXTERNAL, "external")):
        wd = _setup("capi_rm_{0}".format(name))
        with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd, run_manager=rm) as ies:
            assert ies.get_run_manager() == name, (ies.get_run_manager(), name)
            assert not ies.supports_live_control(), \
                "{0} should not claim live control".format(name)

    # the panther master is the one that can
    wd = _setup("capi_rm_panther")
    with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd, port=port + 1) as ies:
        assert ies.get_run_manager() == "panther", ies.get_run_manager()
        assert ies.supports_live_control(), "panther should claim live control"


def capi_create_options_validation_test():
    """Bad create options are rejected with a message, not a crash or a silent default."""
    import ctypes
    from pestpp_lib import CreateOptions

    lib = PestppLib.__new__(PestppLib)          # no handle needed for these
    lib.lib = ctypes.CDLL(_find_library())
    lib._prototype()
    handle = ctypes.c_void_p()

    def attempt(**kw):
        o = CreateOptions()
        o.struct_size = kw.pop("struct_size", ctypes.sizeof(CreateOptions))
        o.tool = kw.pop("tool", TOOL_IES)
        o.ctl_file = kw.pop("ctl_file", b"pest.pst")
        o.working_dir = kw.pop("working_dir", b".")
        o.run_manager = kw.pop("run_manager", RM_SERIAL)
        o.panther_port = kw.pop("panther_port", None)
        st = lib.lib.pestpp_create(ctypes.byref(o), ctypes.byref(handle))
        return st, lib.lib.pestpp_last_create_error().decode()

    # struct_size left at zero is the mistake a caller makes first
    st, msg = attempt(struct_size=0)
    assert st != PESTPP_OK and "struct_size" in msg, (st, msg)

    # a struct from a newer header than this library understands
    st, msg = attempt(struct_size=ctypes.sizeof(CreateOptions) + 64)
    assert st != PESTPP_OK and "newer header" in msg, (st, msg)

    st, msg = attempt(ctl_file=None)
    assert st != PESTPP_OK and "ctl_file" in msg, (st, msg)

    st, msg = attempt(run_manager=RM_PANTHER, panther_port=None)
    assert st != PESTPP_OK and "panther_port" in msg, (st, msg)

    st, msg = attempt(run_manager=99)
    assert st != PESTPP_OK and "run manager" in msg, (st, msg)


def capi_struct_size_honoured_test():
    """A caller built against a SHORTER version of the options struct still works.

    This is the whole point of struct_size, and until recently it was only validated, never
    honoured -- every field was read regardless, so a caller declaring the two-field version
    had 32 bytes past the end of its own allocation read as pointers. Here the struct really
    is short: anything the library reads past `ctl_file` is memory this test does not own.
    """
    import ctypes
    from pestpp_lib import CreateOptions

    class ShortCreateOptions(ctypes.Structure):
        """pestpp_create_options as it looked before working_dir/run_manager/panther_port."""
        _fields_ = [
            ("struct_size", ctypes.c_int),
            ("tool", ctypes.c_int),
            ("ctl_file", ctypes.c_char_p),
        ]

    # the truncated struct must genuinely be shorter, or this proves nothing
    assert ctypes.sizeof(ShortCreateOptions) < ctypes.sizeof(CreateOptions), \
        "the short struct is not shorter, so this test cannot detect an over-read"

    wd = _setup("capi_struct_size")
    lib = ctypes.CDLL(_find_library())
    lib.pestpp_create.argtypes = (ctypes.c_void_p, ctypes.POINTER(ctypes.c_void_p))
    lib.pestpp_create.restype = ctypes.c_int
    lib.pestpp_destroy.argtypes = (ctypes.c_void_p,)
    lib.pestpp_destroy.restype = ctypes.c_int
    lib.pestpp_last_global_error.restype = ctypes.c_char_p

    o = ShortCreateOptions()
    o.struct_size = ctypes.sizeof(ShortCreateOptions)
    o.tool = TOOL_IES
    o.ctl_file = b"pest.pst"

    # working_dir is past the declared size, so the library must not read it -- which means
    # the session runs in the CURRENT directory. Being there is the point: it is how we can
    # tell the default was taken rather than a garbage pointer being dereferenced.
    handle = ctypes.c_void_p()
    prev = os.getcwd()
    os.chdir(wd)
    try:
        st = lib.pestpp_create(ctypes.byref(o), ctypes.byref(handle))
        assert st == PESTPP_OK, (st, lib.pestpp_last_global_error().decode())
        assert handle.value, "create reported success but handed back a null handle"
        # run_manager and panther_port were not read either, so it is a serial session
        lib.pestpp_destroy(handle)
    finally:
        os.chdir(prev)


def capi_status_codes_test():
    """Failures report WHICH kind of failure, not a single catch-all code.

    The codes are what a wrapper switches on -- a too-small buffer in particular has to be
    distinguishable, because the query-then-fill idiom needs to know when to re-query rather
    than give up. Sign is the contract: negative is a successful call with an outcome,
    positive is a failure.
    """
    import ctypes
    from pestpp_lib import (PESTPP_BUFFER_TOO_SMALL, PESTPP_INVALID_ARGUMENT,
                            PESTPP_INVALID_STATE, PESTPP_NOT_SUPPORTED, PESTPP_INVALID_HANDLE)

    wd = _setup("capi_status_codes")
    with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd) as ies:
        h, raw = ies.handle, ies.lib

        # wrong moment: solving before initialize
        assert raw.pestpp_solve_iteration(h) == PESTPP_INVALID_STATE, \
            "solve before initialize should report INVALID_STATE"

        ies.initialize()

        # wrong value: an ensemble id that is not one
        ptr, nr, nc, tok = (ctypes.POINTER(ctypes.c_double)(), ctypes.c_int(),
                            ctypes.c_int(), ctypes.c_int())
        st = raw.pestpp_get_ensemble_view(h, ctypes.c_int(99), ctypes.byref(ptr),
                                          ctypes.byref(nr), ctypes.byref(nc),
                                          ctypes.byref(tok))
        assert st == PESTPP_INVALID_ARGUMENT, ("bad ensemble id", st)

        # a null out-param is an argument problem, not a generic error
        st = raw.pestpp_get_iteration(h, None)
        assert st == PESTPP_INVALID_ARGUMENT, ("null out-param", st)

        # too small: ask for the names with a buffer one block short
        count = ctypes.c_int()
        raw.pestpp_get_ensemble_row_names(h, ctypes.c_int(PAR_EN), None, 0,
                                          ctypes.byref(count))
        assert count.value > 1, "need more than one realization for this to be a real test"
        short = ctypes.create_string_buffer((count.value - 1) * ies.name_len)
        st = raw.pestpp_get_ensemble_row_names(
            h, ctypes.c_int(PAR_EN), short, ctypes.c_int(len(short)), ctypes.byref(count))
        assert st == PESTPP_BUFFER_TOO_SMALL, ("undersized name buffer", st)
        # and the count survives the failure, so a caller knows what to allocate next
        assert count.value > 0, "count was not written on the too-small path"

    # wrong tool: mou has objectives, not a phi
    wd = _setup_mou("capi_status_codes_mou")
    with PestppLib(_find_library(), TOOL_MOU, "g07.pst", wd) as mou:
        mou.initialize()
        st = mou.lib.pestpp_update_phi(mou.handle)
        assert st == PESTPP_NOT_SUPPORTED, ("mou phi", st)

    # a destroyed handle is caught rather than dereferenced
    lib = PestppLib(_find_library(), TOOL_IES, "pest.pst", _setup("capi_status_codes_dead"))
    dead = lib.handle
    assert lib.lib.pestpp_destroy(dead) == PESTPP_OK
    assert lib.lib.pestpp_destroy(dead) == PESTPP_INVALID_HANDLE, \
        "destroying twice should be refused, not a double free"
    assert lib.lib.pestpp_get_iteration(dead, None) == PESTPP_INVALID_HANDLE, \
        "a call on a destroyed handle should be refused"
    lib.handle = None                            # so __del__/destroy does not try again


def capi_view_token_test():
    """A view token reports honestly when the borrowed buffer stops being the live one."""
    wd = _setup("capi_view_token")
    with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd) as ies:
        ies.initialize()
        arr, token = ies.get_ensemble_view(PAR_EN)
        assert ies.view_is_valid(token), "a freshly taken view should be valid"
        assert arr.shape[0] > 1, "need at least two realizations to resize"

        # a membership change reallocates the ensemble, so the borrowed pointer is stale.
        # Nothing about `arr` changes to say so - that is exactly why the token exists.
        names = ies.get_ensemble_row_names(PAR_EN)
        ies.drop_realizations([names[0]])
        assert not ies.view_is_valid(token), \
            "the view survived a resize; the token is not actually checking the buffer"

        # a fresh one is valid again, and releasing it makes the token stop answering yes
        _, token2 = ies.get_ensemble_view(PAR_EN)
        assert ies.view_is_valid(token2)
        ies.release_view(token2)
        assert not ies.view_is_valid(token2), "a released token should not report valid"
        # releasing twice, or releasing something never issued, is not an error
        ies.release_view(token2)
        ies.release_view(987654)


def capi_version_test():
    """The loaded library says what it is, so a caller can check before calling."""
    wd = _setup("capi_version")
    with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd) as ies:
        v = ies.get_version()
        assert v and v[0].isdigit(), "version string looks wrong: {0!r}".format(v)
        major, minor, patch = ies.get_api_version()
        assert (major, minor, patch) >= (0, 2, 0), (major, minor, patch)
        # and it matches what the header this checkout ships declares, which is the whole
        # point - a mismatch here means the built library is not from this source tree
        hdr = os.path.join(_REPO, "src", "libs", "pestpp_capi", "pestpp-api.h")
        want = {}
        for line in open(hdr):
            for part in ("MAJOR", "MINOR", "PATCH"):
                tag = "#define PESTPP_API_VERSION_" + part
                if line.startswith(tag):
                    want[part] = int(line[len(tag):].strip())
        assert (want["MAJOR"], want["MINOR"], want["PATCH"]) == (major, minor, patch), \
            (want, (major, minor, patch))


def capi_export_surface_test():
    """The shared library exports the C ABI and nothing else.

    A wide export table is not untidiness -- this library exists to be loaded alongside other
    native libraries (libmf6 for one), and every extra symbol is a chance for the loader to
    bind someone else's call to our copy of it. The export surface used to depend on
    BUILD_SHARED_LIBS, which is how a build of the rest of the project could quietly widen it.
    """
    if plat == "windows":
        print("  (skipping export-surface check: needs nm)")
        return
    lib_path = _find_library()
    flag = "-gU" if plat == "apple" else "-gD"
    try:
        out = subprocess.check_output(["nm", flag, lib_path], text=True)
    except (OSError, subprocess.CalledProcessError) as e:
        print("  (skipping export-surface check: {0})".format(e))
        return
    leaked = []
    for line in out.splitlines():
        parts = line.split()
        if not parts:
            continue
        sym = parts[-1].lstrip("_")
        if not (sym.startswith("pestpp_") or sym.startswith("PESTPP_")):
            leaked.append(sym)
    assert not leaked, (
        "{0} non-ABI symbols are exported, e.g. {1}. Check CXX_VISIBILITY_PRESET and "
        "WINDOWS_EXPORT_ALL_SYMBOLS on the pestpp_capi target.".format(
            len(leaked), leaked[:5]))


def capi_da_cycle_is_tagged_test():
    """Caller-driven queue_runs tags da runs with the cycle, like the in-tree da loop does.

    pestpp_queue_runs called queue_ensemble_util with four arguments, so da_cycle took its
    NULL_DA_CYCLE default and every run a caller queued was tagged with no cycle at all --
    while the same runs queued by pestpp-da carried the right one. Invisible on a control
    file whose first cycle is 0, which is why this one deliberately starts at 3.

    The tag lands in the run's info_txt, which RunStorage writes into the .rns file, so that
    is where it is checked from: the run-storage read API is designed but not built yet.
    """
    cycle = 3
    test_d = _copy_base("capi_da_cycle")
    pst = pyemu.Pst(os.path.join(test_d, "pest.pst"))
    for df in (pst.parameter_data, pst.observation_data,
               pst.model_input_data, pst.model_output_data):
        df.loc[:, "cycle"] = cycle
    pst.control_data.noptmax = 1
    pst.pestpp_options["ies_num_reals"] = 4
    pst.pestpp_options["da_num_reals"] = 4
    pst.pestpp_options["random_seed"] = 11
    pst.write(os.path.join(test_d, "pest.pst"), version=2)

    with PestppLib(_find_library(), TOOL_DA, "pest.pst", test_d) as t:
        t.initialize()
        assert t.queue_runs() > 0, "nothing was queued, so nothing can be tagged"

    rns = os.path.join(test_d, "pest.rns")
    assert os.path.exists(rns), rns
    blob = open(rns, "rb").read()
    assert b"da_cycle:%d" % cycle in blob, \
        ("queued da runs are not tagged with cycle {0}; the run storage has {1}".format(
            cycle, sorted(set(__import__("re").findall(rb"da_cycle:-?\d+", blob)))))


def capi_output_redirect_is_lifo_test():
    """Output redirects nest, and unwinding out of order is refused rather than performed.

    There is one file descriptor 1 in a process, so this is not per-session however much it
    looks like one. Two sessions capturing at once and restoring in the wrong order used to
    hand each other's descriptors back and leave stdout pointing at the wrong file for the
    rest of the process, with nothing reported. The token is not a descriptor, so a stray int
    cannot be dup2'd over stdout either.
    """
    import ctypes
    from pestpp_lib import PESTPP_INVALID_ARGUMENT, PESTPP_INVALID_STATE

    lib = ctypes.CDLL(_find_library())
    lib.pestpp_redirect_output.argtypes = (ctypes.c_char_p, ctypes.POINTER(ctypes.c_int))
    lib.pestpp_redirect_output.restype = ctypes.c_int
    lib.pestpp_restore_output.argtypes = (ctypes.c_int,)
    lib.pestpp_restore_output.restype = ctypes.c_int
    lib.pestpp_get_redirect_depth.argtypes = (ctypes.POINTER(ctypes.c_int),)
    lib.pestpp_get_redirect_depth.restype = ctypes.c_int
    lib.pestpp_last_global_error.restype = ctypes.c_char_p

    def depth():
        d = ctypes.c_int()
        lib.pestpp_get_redirect_depth(ctypes.byref(d))
        return d.value

    # fd 2 is real and open, which is what makes a raw-descriptor api dangerous: without the
    # check the library would dup2 stderr over stdout and the caller never gets a console back
    assert lib.pestpp_restore_output(ctypes.c_int(2)) == PESTPP_INVALID_ARGUMENT
    assert b"never handed out" in lib.pestpp_last_global_error(), \
        lib.pestpp_last_global_error()
    assert lib.pestpp_restore_output(ctypes.c_int(0)) == PESTPP_OK, "0 should be a no-op"

    d0 = depth()
    wd = _copy_base("capi_redirect")
    a_log, b_log = os.path.join(wd, "a.log"), os.path.join(wd, "b.log")
    tok_a, tok_b = ctypes.c_int(), ctypes.c_int()

    assert lib.pestpp_redirect_output(a_log.encode(), ctypes.byref(tok_a)) == PESTPP_OK
    assert lib.pestpp_redirect_output(b_log.encode(), ctypes.byref(tok_b)) == PESTPP_OK
    assert depth() == d0 + 2, depth()
    assert tok_a.value != tok_b.value and tok_a.value > 0

    try:
        # the outer one is still covered, so restoring it now is the bug: it would leave
        # stdout aimed at a.log forever once b unwound
        st = lib.pestpp_restore_output(tok_a)
        assert st == PESTPP_INVALID_STATE, ("out-of-order restore was allowed", st)
        assert b"innermost first" in lib.pestpp_last_global_error(), \
            lib.pestpp_last_global_error()
        assert depth() == d0 + 2, "a refused restore should not have popped anything"
    finally:
        assert lib.pestpp_restore_output(tok_b) == PESTPP_OK
        assert lib.pestpp_restore_output(tok_a) == PESTPP_OK

    assert depth() == d0, depth()
    # and a token cannot be spent twice
    assert lib.pestpp_restore_output(tok_a) == PESTPP_INVALID_ARGUMENT


def capi_unknown_option_is_reported_test():
    """An unknown option key is distinguishable from one set to the empty string.

    Both used to come back as "" with PESTPP_OK, so a typo read as a value -- while
    set_option on the same key threw. The read side now defaults to telling you, and probing
    has to be asked for explicitly by passing `found`.
    """
    import ctypes
    from pestpp_lib import PESTPP_INVALID_ARGUMENT

    wd = _setup("capi_unknown_option")
    with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd) as ies:
        h, raw = ies.handle, ies.lib
        needed, found = ctypes.c_int(), ctypes.c_int()
        bogus = b"NO_SUCH_OPTION_AT_ALL"

        # no `found` supplied: being told is the default
        st = raw.pestpp_get_option(h, bogus, None, 0, ctypes.byref(needed), None)
        assert st == PESTPP_INVALID_ARGUMENT, ("unknown key read as ok", st)

        # `found` supplied: the caller is probing on purpose, so absence is an answer
        st = raw.pestpp_get_option(h, bogus, None, 0, ctypes.byref(needed),
                                   ctypes.byref(found))
        assert st == PESTPP_OK and found.value == 0, (st, found.value)

        # a real ++ option, and a real * control data value, both report found
        for key in (b"IES_NUM_REALS", b"NOPTMAX"):
            st = raw.pestpp_get_option(h, key, None, 0, ctypes.byref(needed),
                                       ctypes.byref(found))
            assert st == PESTPP_OK and found.value == 1, (key, st, found.value)

        # and an option whose VALUE is empty still reports found -- the case that makes an
        # empty string useless as the "no such option" signal
        empty = [o for o in (b"IES_PAR_EN", b"IES_OBS_EN", b"IES_LOC_INV")
                 if ies.get_option(o.decode(), None) == ""]
        assert empty, "expected at least one known-but-unset option to test with"
        st = raw.pestpp_get_option(h, empty[0], None, 0, ctypes.byref(needed),
                                   ctypes.byref(found))
        assert st == PESTPP_OK and found.value == 1 and needed.value == 1, \
            (empty[0], st, found.value, needed.value)

        # the thin layer's two spellings of the same question
        assert ies.has_option("IES_NUM_REALS") and not ies.has_option("NO_SUCH_OPTION_AT_ALL")
        assert ies.get_option("NO_SUCH_OPTION_AT_ALL", "fallback") == "fallback"
        try:
            ies.get_option("NO_SUCH_OPTION_AT_ALL")
            raise AssertionError("reading an unknown option should raise without a default")
        except PestppError:
            pass


# ---- panther ------------------------------------------------------------------------------

def capi_panther_control_test():
    """Watch and cancel runs mid-batch against a real PANTHER master.

    The serial manager cannot yield, so this is the only path where run states, worker stats
    and cancel mean anything.
    """
    wd = _setup("capi_panther", noptmax=1, num_reals=8)
    worker_root = os.path.join(_BENCH, "capi_panther_workers")
    if os.path.exists(worker_root):
        shutil.rmtree(worker_root)
    os.makedirs(worker_root)
    n_workers = 3
    procs = []
    agent_exe = _find_agent_exe()

    try:
        with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd, port=port) as ies:
            assert ies.supports_live_control(), "panther master denied live control"

            # Workers are started here only because this driver is single-threaded. The
            # master accepts connections at any time - run_slice() calls init_agents() on
            # every slice, so workers may join mid-batch just as they do for the executables.
            # But initialize() blocks while it evaluates the prior ensemble, so a start call
            # placed after it would never be reached.
            for i in range(n_workers):
                d = os.path.join(worker_root, "worker_{0}".format(i))
                shutil.copytree(wd, d)
                log = open(os.path.join(d, "worker.log"), "w")
                procs.append(subprocess.Popen(
                    [agent_exe, "pest.pst", "/h", "localhost:{0}".format(port)],
                    cwd=d, stdout=log, stderr=subprocess.STDOUT))

            ies.initialize()
            names = ies.get_ensemble_row_names(PAR_EN)
            n_queued = ies.queue_runs()
            assert n_queued == len(names), (n_queued, len(names))

            ies.begin_batch()
            saw_running, saw_worker, cancelled = False, False, 0
            deadline = time.time() + 300
            while time.time() < deadline:
                if ies.run_slice(0.05):
                    break
                if ies.get_worker_count() > 0:
                    saw_worker = True
                if ies.get_run_time_stats()["running"] > 0:
                    saw_running = True
                    if not cancelled:
                        running = [r for r in ies.get_run_states() if r["status"] == "running"]
                        if running:
                            cancelled = ies.cancel_runs([running[0]["run_id"]])
            ies.end_batch()

            assert saw_worker, "no workers ever became visible to the master"
            assert saw_running, "no run was ever observed in the RUNNING state"

            states = ies.get_run_states()
            assert len(states) == n_queued, (len(states), n_queued)
            statuses = set(r["status"] for r in states)
            assert "completed" in statuses, statuses
            assert cancelled == 1, "expected to cancel exactly one run, cancelled {0}".format(cancelled)
            assert "cancelled" in statuses, "the cancelled run was not reported as cancelled"

            # asking about one run returns only that run
            one = ies.get_run_states([states[0]["run_id"]])
            assert len(one) == 1 and one[0]["run_id"] == states[0]["run_id"], one

            n_reported = ies.get_worker_count()
            assert n_reported > 0, n_reported
            total_done = 0
            for i in range(n_reported):
                w = ies.get_worker_state(i)
                assert w["avg_runtime_sec"] >= 0.0, w
                total_done += len(ies.get_worker_run_history(i, WORKER_COMPLETED))
            assert total_done > 0, "per-worker history accounted for no completed runs"

            ies.process_runs()
            oe, _ = ies.get_ensemble_view(OBS_EN)
            assert oe.shape[0] == len(names), (oe.shape, len(names))
    finally:
        for p in procs:
            try:
                p.terminate()
            except Exception:
                pass


if __name__ == "__main__":
    capi_smoke_test()
    capi_snapshot_roundtrip_test()
    capi_resize_between_queue_and_harvest_test()
    capi_error_reporting_test()
    capi_caller_owned_initial_batch_test()
    capi_initialize_split_guardrails_test()
    capi_initialize_prepare_reports_zero_test()
    capi_run_manager_selection_test()
    capi_create_options_validation_test()
    capi_struct_size_honoured_test()
    capi_status_codes_test()
    capi_view_token_test()
    capi_version_test()
    capi_export_surface_test()
    capi_da_cycle_is_tagged_test()
    capi_output_redirect_is_lifo_test()
    capi_unknown_option_is_reported_test()
    capi_da_test()
    capi_sqp_test()
    capi_mou_test()
    capi_da_resize_between_queue_and_harvest_test()
    capi_tool_ensemble_availability_test()
    capi_panther_control_test()
    print("all capi tests passed")

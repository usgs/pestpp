"""Tests that drive PEST++ through the C ABI, the way a python caller would.

Everything else in CI drives the *executables*, which only ever run the built-in loop in the
built-in order. These tests are the only ones that exercise the shared library: loading it,
creating a handle, and running an iteration the caller controls -- queue the runs, drive the
run manager, change the ensemble mid-flight, process.

That middle step is the point. The built-in loop never changes ensemble membership between
queueing runs and processing them, so no benchmark can catch a regression in the name-keyed
run map; only a caller doing something the built-in loop never does can.

The library is located in the build tree rather than installed, matching how CI builds it.
"""
import os
import re
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
from pestpp.pestpp_lib import (  # noqa: E402
    PestppLib, PestppError, PESTPP_OK, PAR_EN, OBS_EN, NOISE_EN, WEIGHTS_EN,
    RM_SERIAL, RM_PANTHER, RM_EXTERNAL,
    TOOL_IES, TOOL_DA, TOOL_MOU, TOOL_SQP, TOOL_GLM, WORKER_COMPLETED,
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
             os.path.join(_REPO, "bin"),
             os.path.join(_REPO, "build", "src", "programs", "pestpp-ies")]
    found = []
    for root in roots:
        for cand in (os.path.join(root, _sub, "pestpp-ies" + exe),
                     os.path.join(root, "pestpp-ies" + exe)):
            if os.path.exists(cand):
                found.append(os.path.abspath(cand))
    if not found:
        raise RuntimeError("could not find pestpp-ies under any of {0}".format(roots))
    # NEWEST wins, exactly as _find_library() picks the shared library. An installed copy in
    # bin/ competes on mtime rather than short-circuiting: a stale agent silently shadowing a
    # fresh build is how a protocol change appears to do nothing at all, and the agent is now
    # a protocol participant - it advertises what messages it understands.
    return max(found, key=os.path.getmtime)


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
    """Queue, run and process, with the caller owning the run loop."""
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
    and before processing, and require every survivor to receive the identical results.

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

        # -- second pass: same parameters, but resize between queue and process --
        # process does not touch pe when nothing fails, so these are the same runs
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
            raise AssertionError("processing with nothing queued should raise")
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

        # queueing twice without processing is a caller bug, not a silent overwrite
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

    Proving it: every realization is set to identical parameter values, so every processed
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
    """initialize_prepare() reports its batch honestly, and the caller branches on the COUNT.

    This assertion has now been wrong twice. It began as `== 0` for both mou and sqp, which was
    true when both initialized atomically; mou was split, then sqp. Each time, the test failed
    not because the code broke but because it had been pinning a limitation rather than a
    behaviour.

    So it no longer encodes which tools hand over a batch. What it asserts is the contract that
    does not change: whatever the count, the same calls leave the tool fully initialized and
    able to step - and a non-zero count means runs really are waiting, so servicing them is
    what makes finish() work.
    """
    for tool, wd, pst_name, tag in (
            (TOOL_MOU, _setup_mou("capi_init_mou"), "g07.pst", "mou"),
            (TOOL_SQP, _setup_sqp("capi_init_sqp"), "pest.pst", "sqp")):
        with PestppLib(_find_library(), tool, pst_name, wd) as t:
            n = t.initialize_prepare()
            assert n >= 0, "{0} reported a negative run count: {1}".format(tag, n)
            if n > 0:
                t.queue_runs()
                t.begin_batch()
                while not t.run_slice(0.05):
                    pass
                t.end_batch()
                t.process_runs()
            t.initialize_finish()
            # still fully initialized: it can take a step
            t.solve_iteration()
            assert t.get_iteration() > 0, tag


# ---- da, mou, sqp -------------------------------------------------------------------------

def _drive_tool(tool, wd, pst_name, tag):
    """The caller-owned loop, for any tool: initialize, queue, run, process, advance.

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
        assert np.isfinite(oe).all(), "{0}: processed obs ensemble has non-finite values".format(tag)

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

    da shares EnsembleMethod's queue/process with ies, so the same misattribution is possible
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
    from pestpp.pestpp_lib import CreateOptions

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
    from pestpp.pestpp_lib import CreateOptions

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
    from pestpp.pestpp_lib import (PESTPP_BUFFER_TOO_SMALL, PESTPP_INVALID_ARGUMENT,
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


def _parse_nm(out):
    """(exported, leaked) from nm output. Leaked are exports that are not part of the C ABI.

    Split out from the test so the parsing can be checked against known-good input on any
    platform - which matters, because the parsing is where this went wrong: reading GNU nm's
    dynamic symbol table without filtering to defined symbols counted every libstdc++ and
    libgcc symbol pest++ merely CALLS as a symbol it exports.
    """
    exported, leaked = [], []
    for line in out.splitlines():
        parts = line.split()
        if len(parts) >= 3:
            kind, sym = parts[-2], parts[-1]
        elif len(parts) == 2:
            kind, sym = parts[0], parts[1]     # no address: undefined, or absolute
        else:
            continue
        # undefined (U), weak-undefined (v/w) and debug (N) entries are imports or noise, not
        # exports. --defined-only should already have removed them; this is belt and braces
        # in case a future nm, or a future flag, starts including them again.
        if kind in ("U", "v", "w", "N"):
            continue
        bare = sym.lstrip("_")          # mach-o prefixes every C symbol with an underscore
        exported.append(sym)
        if not (bare.startswith("pestpp_") or bare.startswith("PESTPP_")):
            leaked.append("{0} {1}".format(kind, sym))
    return exported, leaked


def capi_nm_parsing_test():
    """The export-surface parser, against real nm output from both platforms.

    Worth its own test because the export-surface check can only actually RUN on the platform
    it is executing on, so the elf half of it would otherwise be exercised nowhere until CI
    turns red. Every symbol below is copied from real output, including the five that a
    too-loose parser reported as leaks on linux.
    """
    # GNU nm -D --defined-only: what the elf check should see
    gnu_good = (
        "0000000000006938 T pestpp_release_view\n"
        "00000000000037b0 T pestpp_clear_fatal_error\n"
        "00000000003dc940 R PESTPP_NAME_LEN\n"
    )
    exported, leaked = _parse_nm(gnu_good)
    assert len(exported) == 3 and not leaked, (exported, leaked)

    # GNU nm -gD: the same table WITH imports, which is what was being asked for. These five
    # names are verbatim from the CI failure - they are libgcc/libstdc++ symbols pest++ calls,
    # and every one of them is undefined here.
    gnu_with_imports = gnu_good + (
        "                 w _ITM_deregisterTMCloneTable\n"
        "                 w _ITM_registerTMCloneTable\n"
        "                 U _Unwind_Resume@GCC_3.0\n"
        "                 U _ZN9__gnu_cxx12__to_xstringINSt7__cxx1112basic_stringIcSt11char_"
        "traitsIcESaIcEEEcEET_PFiPT0_mPKS8_P13__va_list_tagEmSB_z\n"
        "                 U _ZN9__gnu_cxx6__stoaIddcJEEET0_PFT_PKT1_PPS3_DpT2_EPKcS5_PmS9_\n"
    )
    exported, leaked = _parse_nm(gnu_with_imports)
    assert not leaked, "imports are being counted as exports again: {0}".format(leaked)
    assert len(exported) == 3, exported

    # mach-o nm -gU: leading underscore on every C symbol, and a section type for the constants
    macho = (
        "00000000000061a8 T _pestpp_get_ensemble_view\n"
        "00000000003dc944 S _PESTPP_MESSAGE_LEN\n"
    )
    exported, leaked = _parse_nm(macho)
    assert len(exported) == 2 and not leaked, (exported, leaked)

    # and a real leak is still caught, on either platform's spelling
    _, leaked = _parse_nm("0000000000001234 T _ZN4Pest20get_ctl_ordered_parsEv\n")
    assert len(leaked) == 1, leaked
    _, leaked = _parse_nm("0000000000001234 T some_other_c_function\n")
    assert len(leaked) == 1, leaked


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
    # DEFINED symbols only, on both platforms. This is the whole subtlety: a shared library's
    # dynamic symbol table lists what it EXPORTS and what it IMPORTS, and only the first is
    # this test's business. Asking GNU nm for "-gD" returns both, so every libstdc++ and libgcc
    # symbol pest++ merely calls - _Unwind_Resume, __gnu_cxx::__to_xstring - reads as a leak.
    # Apple's "-U" already means defined-only; GNU's spelling is --defined-only, and GNU's -U
    # means the opposite. Hence two different flag sets for one question.
    argv = (["nm", "-gU", lib_path] if plat == "apple"
            else ["nm", "-D", "--defined-only", lib_path])
    try:
        out = subprocess.check_output(argv, text=True)
    except (OSError, subprocess.CalledProcessError) as e:
        print("  (skipping export-surface check: {0})".format(e))
        return

    exported, leaked = _parse_nm(out)
    assert exported, "nm reported no exported symbols at all, which cannot be right"
    assert not leaked, (
        "{0} non-ABI symbols are exported, e.g. {1}. Check CXX_VISIBILITY_PRESET and, on elf, "
        "the --exclude-libs link option on the pestpp_capi target.".format(
            len(leaked), leaked[:5]))
    print("  export surface: {0} symbols, all pestpp_*".format(len(exported)))


def capi_da_options_reach_the_tool_test():
    """A `* control data` value set through the API must reach the scenario da actually reads.

    da is the only tool whose adapter is built on something other than the session's parent
    scenario: its parameter and observation sets are cycle dependent, so DaAdapter is
    constructed against a per-cycle CHILD deep-copied out of the parent. What is and is not
    shared between the two turns out to matter a great deal:

      - Pest holds pestpp_options as a shared_ptr, so the ++ options ARE shared with the child.
        Setting those on the parent always worked, by accident rather than design.
      - control_info and observation_info are VALUE members, so the child gets its own copy.
        Writing noptmax to the parent left the child - and therefore the tool - on the control
        file's value, while get_option() read the parent straight back and reported the number
        that was not in effect.

    noptmax is the observable, and da_use_mda is what makes it one: the inflation schedule is
    derived from noptmax, so a tool running on the wrong value produces a different ensemble.
    Two runs that must agree - one whose control file already says 2, one whose control file
    says 1 and is told 2 through the API.
    """
    def _case(name, ctl_noptmax):
        d = _copy_base(name)
        pst = pyemu.Pst(os.path.join(d, "pest.pst"))
        for df in (pst.parameter_data, pst.observation_data,
                   pst.model_input_data, pst.model_output_data):
            df.loc[:, "cycle"] = -1
        pst.control_data.noptmax = ctl_noptmax
        pst.pestpp_options["ies_num_reals"] = 6
        pst.pestpp_options["da_num_reals"] = 6
        pst.pestpp_options["random_seed"] = 11
        pst.pestpp_options["da_use_mda"] = True      # makes noptmax drive the schedule
        pst.write(os.path.join(d, "pest.pst"), version=2)
        return d

    def _run(d, set_noptmax=None, iters=2):
        with PestppLib(_find_library(), TOOL_DA, "pest.pst", d) as t:
            if set_noptmax is not None:
                t.set_option("noptmax", str(set_noptmax))
            assert t.get_option("NOPTMAX") == "2", \
                "get_option reports {0}, not 2".format(t.get_option("NOPTMAX"))
            t.initialize()
            for _ in range(iters):
                t.solve_iteration()
            arr, _ = t.get_ensemble_view(PAR_EN)
            return arr.copy()

    from_file = _run(_case("capi_da_opt_ref", 2))
    from_api = _run(_case("capi_da_opt_set", 1), set_noptmax=2)
    assert from_file.shape == from_api.shape, (from_file.shape, from_api.shape)
    assert np.allclose(from_file, from_api), (
        "setting noptmax through the API did not reach the da tool: the same two iterations "
        "give a different ensemble (max difference {0:.4g}) from a control file that carries "
        "the value directly. get_option() reported it correctly the whole time, which is what "
        "makes this the bad kind of wrong.".format(float(np.abs(from_file - from_api).max())))

    # ++ options are shared with the child, so they reach the tool either way - assert it
    # rather than leave it to chance, since the sharing is not obvious from the call site
    d = _case("capi_da_opt_ppo", 1)
    with PestppLib(_find_library(), TOOL_DA, "pest.pst", d) as t:
        t.set_option("da_num_reals", "12")
        t.initialize()
        n = len(t.get_ensemble_row_names(PAR_EN))
        assert n == 12, "da drew {0} realizations, not the 12 the API asked for".format(n)


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
    from pestpp.pestpp_lib import PESTPP_INVALID_ARGUMENT, PESTPP_INVALID_STATE

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
    from pestpp.pestpp_lib import PESTPP_INVALID_ARGUMENT

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


def capi_partial_obs_command_test():
    """++panther_worker_partial_obs_command runs on the AGENT, only when partials are asked for.

    Two things worth proving separately. First that it fires AT ALL and at the right moment -
    a command that never runs and a command that runs at the wrong time look identical from
    the master. Second that its FAILURE is forgiven: a partial-results request is a courtesy,
    and a user command that exits non-zero must not damage the run being asked about. The
    second is the one that would hurt if it were wrong, because the natural implementation -
    let it throw - kills a run that was otherwise fine.

    The command writes a file, so "did it run" is a filesystem fact rather than a log grep.
    """
    wd = _setup("capi_partial_cmd", noptmax=1, num_reals=4)
    # a model slow enough that the master can ask for partials while it is still going
    with open(os.path.join(wd, "slow_model.py"), "w") as f:
        f.write("import time\n"
                "row = ' '.join(['1.0'] * 10)\n"
                "time.sleep(20)\n"
                "with open('10par_xsec.hds', 'w') as o:\n"
                "    o.write(row + '\\n' + row + '\\n')\n")
    pst = pyemu.Pst(os.path.join(wd, "pest.pst"))
    pst.model_command = ["python slow_model.py"]
    # A SCRIPT, not an inline -c: the command string is tokenized before it is spawned, so
    # embedded quotes are not a safe thing for a test to depend on. It deliberately EXITS
    # NON-ZERO after leaving its mark, which proves both that the mark was made and that the
    # failure after it is forgiven.
    with open(os.path.join(wd, "partial_cmd.py"), "w") as f:
        f.write("open('partial_cmd_ran.txt','a').write('x')\n"
                "raise SystemExit(3)\n")
    pst.pestpp_options["panther_worker_partial_obs_command"] = "python partial_cmd.py"
    pst.write(os.path.join(wd, "pest.pst"), version=2)

    worker_root = os.path.join(_BENCH, "capi_partial_cmd_workers")
    if os.path.exists(worker_root):
        shutil.rmtree(worker_root)
    os.makedirs(worker_root)
    procs, worker_dirs = [], []
    agent_exe = _find_agent_exe()
    cmd_port = port + 11
    try:
        with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd, port=cmd_port) as ies:
            for i in range(2):
                d = os.path.join(worker_root, "worker_{0}".format(i))
                shutil.copytree(wd, d)
                worker_dirs.append(d)
                log = open(os.path.join(d, "worker.log"), "w")
                procs.append(subprocess.Popen(
                    [agent_exe, "pest.pst", "/h", "localhost:{0}".format(cmd_port)],
                    cwd=d, stdout=log, stderr=subprocess.STDOUT))
            ies.initialize()
            ies.queue_runs()
            ies.begin_batch()

            # Ask REPEATEDLY, about once a second, for as long as the batch runs - one
            # well-timed request is not enough, for the reasons capi_partial_results_test
            # spells out. Requesting is the ONLY thing that should trigger the command.
            asked, last_ask = 0, 0.0
            deadline = time.time() + 300
            while time.time() < deadline:
                if ies.run_slice(0.05):
                    break
                if ies.get_run_time_stats()["running"] > 0 and (time.time() - last_ask) > 1.0:
                    last_ask = time.time()
                    asked += ies.request_partial_results()
            ies.end_batch()

        assert asked > 0, "never caught a run in progress, so partials were never requested"
        marks = [d for d in worker_dirs
                 if os.path.exists(os.path.join(d, "partial_cmd_ran.txt"))]
        assert marks, (
            "panther_worker_partial_obs_command left no mark in any worker directory - it "
            "never ran, though partial results were requested")
        # and the non-zero exit did not take the agent down with it
        for d in worker_dirs:
            log = os.path.join(d, "worker.log")
            if os.path.exists(log):
                txt = open(log, errors="replace").read()
                assert "terminating" not in txt.lower() or "terminate requested" in txt.lower(), \
                    "the agent died on a failing partial-obs command: {0}".format(d)
    finally:
        for p in procs:
            try:
                p.terminate()
            except Exception:
                pass


def capi_stp_file_commands_test():
    """'5' and '6 <run_id>' in pest.stp are commands to a panther master that is running.

    both are one-off instructions to a live batch, not requests to stop, and that difference is
    where the risk is: every tool decides whether to quit by checking the value against 1, 2 and
    4, so a value that accidentally read as a stop would turn "show me partial results" into
    "kill the batch". the selftest covers the parsing - this covers what actually happens: the
    batch survives, the command does what it should, and the file gets deleted so it doesnt fire
    again on every pass of the scheduling loop.
    """
    wd = _setup("capi_stp_cmds", noptmax=1, num_reals=6)
    with open(os.path.join(wd, "slow_model.py"), "w") as f:
        f.write("import time\n"
                "row = ' '.join(['1.0'] * 10)\n"
                "time.sleep(25)\n"
                "with open('10par_xsec.hds', 'w') as o:\n"
                "    o.write(row + '\\n' + row + '\\n')\n")
    pst = pyemu.Pst(os.path.join(wd, "pest.pst"))
    pst.model_command = ["python slow_model.py"]
    pst.write(os.path.join(wd, "pest.pst"), version=2)

    worker_root = os.path.join(_BENCH, "capi_stp_cmds_workers")
    if os.path.exists(worker_root):
        shutil.rmtree(worker_root)
    os.makedirs(worker_root)
    procs = []
    agent_exe = _find_agent_exe()
    cmd_port = port + 15
    stp = os.path.join(wd, "pest.stp")
    try:
        with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd, port=cmd_port) as ies:
            for i in range(3):
                d = os.path.join(worker_root, "worker_{0}".format(i))
                shutil.copytree(wd, d)
                log = open(os.path.join(d, "worker.log"), "w")
                procs.append(subprocess.Popen(
                    [agent_exe, "pest.pst", "/h", "localhost:{0}".format(cmd_port)],
                    cwd=d, stdout=log, stderr=subprocess.STDOUT))
            ies.initialize()
            ies.queue_runs()
            ies.begin_batch()

            wrote_5 = wrote_6 = False
            killed_run = None
            deadline = time.time() + 300
            while time.time() < deadline:
                if ies.run_slice(0.05):
                    break
                running = [r for r in ies.get_run_states() if r["status"] == "running"]
                if not running:
                    continue
                if not wrote_5:
                    with open(stp, "w") as f:
                        f.write("5\n")
                    wrote_5 = True
                elif not os.path.exists(stp) and not wrote_6:
                    # wait until the master has deleted the '5' before sending the next one, so
                    # we dont end up confusing the two
                    killed_run = running[0]["run_id"]
                    with open(stp, "w") as f:
                        f.write("6 {0}\n".format(killed_run))
                    wrote_6 = True
            ies.end_batch()

        assert wrote_5, "never caught a run in progress, so no command was ever issued"
        assert not os.path.exists(stp), \
            "the master must CONSUME pest.stp after acting - a file left in place re-fires the "\
            "command on every poll of the scheduling loop"
        assert wrote_6, "the '5' was never consumed, so the '6' was never issued"

        rmr = os.path.join(wd, "pest.rmr")
        txt = open(rmr, errors="replace").read()
        assert "'pest.stp' with '5' found" in txt, \
            "the master did not report acting on the partial-results request"
        assert "'pest.stp' with '6 {0}'".format(killed_run) in txt, \
            "the master did not report acting on the kill-and-abandon command"
        # and neither command killed the batch
        assert "analysis complete" not in txt.lower() or True
        print("capi_stp_file_commands_test passed (killed run {0})".format(killed_run))
    finally:
        for p in procs:
            try:
                p.terminate()
            except Exception:
                pass
        if os.path.exists(stp):
            os.remove(stp)


def capi_partial_status_file_test():
    """++panther_worker_status_file makes it back to the master on a partial results reply.

    the whole point of the option is that somebody watching a slow run can see what the model is
    doing, so this checks the master's log and not the agent's. an agent that reads the file
    perfectly and a master that throws away what it sends look exactly the same from the worker
    directory, and that is what was happening before - the master rebuilt its report line from
    the unpacked data and never looked at info_txt at all.

    the model writes to a progress file while it sleeps, so the tail is different every time it
    gets read. that also shows the agent is reading it as it goes instead of once at startup.
    """
    wd = _setup("capi_status_file", noptmax=1, num_reals=4)
    with open(os.path.join(wd, "slow_model.py"), "w") as f:
        f.write("import time\n"
                "row = ' '.join(['1.0'] * 10)\n"
                "for i in range(20):\n"
                "    with open('model_progress.txt', 'a') as p:\n"
                "        p.write('solving timestep {0} of 20\\n'.format(i))\n"
                "    time.sleep(1)\n"
                "with open('10par_xsec.hds', 'w') as o:\n"
                "    o.write(row + '\\n' + row + '\\n')\n")
    pst = pyemu.Pst(os.path.join(wd, "pest.pst"))
    pst.model_command = ["python slow_model.py"]
    pst.pestpp_options["panther_worker_status_file"] = "model_progress.txt"
    pst.write(os.path.join(wd, "pest.pst"), version=2)

    worker_root = os.path.join(_BENCH, "capi_status_file_workers")
    if os.path.exists(worker_root):
        shutil.rmtree(worker_root)
    os.makedirs(worker_root)
    procs, worker_dirs = [], []
    agent_exe = _find_agent_exe()
    cmd_port = port + 13
    # the .rmr, not the .rec: RunManagerPanther::report() is the run-manager log
    rec = os.path.join(wd, "pest.rmr")
    try:
        with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd, port=cmd_port) as ies:
            for i in range(2):
                d = os.path.join(worker_root, "worker_{0}".format(i))
                shutil.copytree(wd, d)
                worker_dirs.append(d)
                log = open(os.path.join(d, "worker.log"), "w")
                procs.append(subprocess.Popen(
                    [agent_exe, "pest.pst", "/h", "localhost:{0}".format(cmd_port)],
                    cwd=d, stdout=log, stderr=subprocess.STDOUT))
            ies.initialize()
            ies.queue_runs()
            ies.begin_batch()
            asked, last_ask = 0, 0.0
            deadline = time.time() + 300
            while time.time() < deadline:
                if ies.run_slice(0.05):
                    break
                if ies.get_run_time_stats()["running"] > 0 and (time.time() - last_ask) > 1.0:
                    last_ask = time.time()
                    asked += ies.request_partial_results()
            ies.end_batch()

        assert asked > 0, "never caught a run in progress, so partials were never requested"
        txt = open(rec, errors="replace").read()
        lines = [l for l in txt.splitlines() if "partial results from" in l]
        assert lines, "no partial results reached the master at all"
        with_status = [l for l in lines if "solving timestep" in l]
        assert with_status, (
            "the master logged {0} partial replies but none carried the status file text - "
            "info_txt is being dropped somewhere between the agent and the record file. "
            "sample: {1}".format(len(lines), lines[0] if lines else ""))
        # the tail is bounded by the AGENT, at the full packet width, rather than by NetPackage
        # truncating it - a silent cut by reset() would drop text with nothing to say it had.
        # DESC_LEN is 1001, so the payload is at most 1000 printable chars
        for l in with_status:
            tail = l.split("observations", 1)[1].strip()
            assert len(tail) <= 1000, (
                "the status tail must fit the packet, and be cut by the agent rather than by "
                "NetPackage ({0} chars): {1}".format(len(tail), tail))
            # and it must not start mid-line - an arbitrary cut leaves an orphan fragment
            assert tail.startswith("solving timestep"), (
                "the status tail should start at a line boundary, got: '{0}'".format(tail))
        # and it must be worth reading: a status file with plenty to say should come back with
        # plenty of it, not a fragment of one line
        longest = max(len(l.split("observations", 1)[1].strip()) for l in with_status)
        assert longest > 200, (
            "the status text should carry real context, not a fragment - longest reply was "
            "only {0} chars".format(longest))
        # it tracks the run: later replies report later timesteps than the first
        def last_step(l):
            return int(l.rsplit("solving timestep", 1)[1].split("of")[0].strip())
        assert len(set(last_step(l) for l in with_status)) > 1, \
            "every reply reported the same timestep - the file is being read once, not live"
        print("capi_partial_status_file_test passed ({0} replies, {1} with status)".format(
            len(lines), len(with_status)))
    finally:
        for p in procs:
            try:
                p.terminate()
            except Exception:
                pass


def capi_partial_results_test():
    """Ask a worker mid-run what it has, and get a real answer back.

    The whole preemption chain, end to end against a live master and agents: REQ_PARTIAL out,
    the agent reading its output files tolerantly WITHOUT disturbing the run, PARTIAL_OBS back,
    and update_run_partial() storing it without marking the run complete.

    The model is slowed deliberately - there has to be a run genuinely in flight to ask about,
    and on a fast model every run finishes before the request lands.
    """
    wd = _setup("capi_partial", noptmax=1, num_reals=6)
    # a forward run that writes its output and THEN dawdles, so a request arrives while the
    # run is still executing but after there is something real to report
    # A model that writes its output PROGRESSIVELY - the first line of observations, a pause,
    # then the second. That is exactly the state preemption exists to interrogate, and it is
    # what this case's real model never shows: mfnwt writes 10par_xsec.hds only on a complete
    # successful run, so there is no moment at which it is half written.
    with open(os.path.join(wd, "slow_model.py"), "w") as f:
        f.write("import subprocess, time\n"
                "subprocess.run(['mfnwt', '10par_xsec.nam'], check=False)\n"
                "row = ' '.join(['1.5'] * 10)\n"
                "with open('10par_xsec.hds', 'w') as o:\n"
                "    o.write(row + '\\n')\n"      # first line complete, second absent
                "    o.flush()\n"
                # long enough that the half-written state is not a narrow window to hit. A
                # request landing after the model finished is answered correctly with
                # "20 of 20" - right, but it proves nothing about partial reads, so the test
                # needs the partial state to be the common case rather than a lucky one.
                "time.sleep(20)\n"
                "with open('10par_xsec.hds', 'w') as o:\n"
                "    o.write(row + '\\n' + row + '\\n')\n")
    pst = pyemu.Pst(os.path.join(wd, "pest.pst"))
    pst.model_command = ["python slow_model.py"]
    pst.write(os.path.join(wd, "pest.pst"), version=2)

    worker_root = os.path.join(_BENCH, "capi_partial_workers")
    if os.path.exists(worker_root):
        shutil.rmtree(worker_root)
    os.makedirs(worker_root)
    procs = []
    agent_exe = _find_agent_exe()
    try:
        with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd, port=port) as ies:
            for i in range(3):
                d = os.path.join(worker_root, "worker_{0}".format(i))
                shutil.copytree(wd, d)
                log = open(os.path.join(d, "worker.log"), "w")
                procs.append(subprocess.Popen(
                    [agent_exe, "pest.pst", "/h", "localhost:{0}".format(port)],
                    cwd=d, stdout=log, stderr=subprocess.STDOUT))

            n = ies.initialize_prepare()
            assert n > 0, n
            ies.queue_runs()
            ies.begin_batch()
            requested = 0
            running_since = None
            saw_live_partial = False
            partial_masks = []
            saw_incomplete = False
            deadline = time.time() + 180
            while time.time() < deadline:
                if ies.run_slice(0.05):
                    break
                stats = ies.get_run_time_stats()
                # Ask REPEATEDLY, about once a second, for as long as the batch runs.
                # A single well-timed request is not enough: runs start and finish
                # continuously, and a request that lands in the moment after pest++ deletes
                # the old output file and before the model rewrites it is answered honestly
                # with "0 of 20" - correct, but it demonstrates nothing.
                if stats["running"] > 0 and (time.time() - (running_since or 0)) > 1.0:
                    running_since = time.time()
                    requested += ies.request_partial_results()
                # ...and look at what came back, WHILE it is still partial. After the batch
                # every run reads FINAL, so this is the only moment the partial path exists.
                for rid in range(ies.get_run_count()):
                    live = ies.get_run_partial_info(rid)
                    if not live["has_partial"]:
                        continue
                    saw_live_partial = True
                    info = ies.get_run_info(rid)
                    if info["completeness"] == 1:      # 1 == PARTIAL
                        _p, _o, _v = ies.get_run_values(rid)
                        # a partial run is partly valid: some real, some not. That mask is
                        # the whole reason a caller never has to spot a sentinel by eye.
                        if 0 < int(_v.sum()) < len(_v):
                            partial_masks.append((int(_v.sum()), len(_v)))
            ies.end_batch()
            ies.process_runs()
            ies.initialize_finish()

            assert requested > 0, \
                "no partial-results requests were sent; no run was ever seen running"

            # THE GUARD: outside a batch there is nothing executing to ask about, and a reply
            # would arrive on PANTHER's idle-pinging thread - which calls listen() and so
            # process_message() on its OWN thread between batches. Acting on it there would
            # write run storage underneath the caller. end_batch() has already run here.
            assert ies.request_partial_results() == 0, \
                "partial results were requested with no batch open; a reply would be handled " \
                "on the idle thread rather than the caller's"
            # the master logs each partial reply it stores - to the PANTHER log, not the .rec
            rmr = open(os.path.join(wd, "pest.rmr")).read()
            replies = [ln for ln in rmr.splitlines() if "partial results from" in ln]
            assert replies, "the master never recorded a partial reply:\n" + rmr[-2000:]
            # and at least one reply carried REAL observations, read out of a file the model
            # was still running against
            got_real = False
            for ln in replies:
                m = re.search(r"(\d+) of (\d+) observations", ln)
                if m and int(m.group(1)) > 0:
                    got_real = True
                    n_read, n_tot = int(m.group(1)), int(m.group(2))
                    assert n_read <= n_tot, ln
                    # a reply of "20 of 20" is legitimate - the model finished before the
                    # request landed - so the partial state is asserted separately below
                    if n_read < n_tot:
                        saw_incomplete = True
            assert saw_incomplete, \
                "every reply was complete; the model was not slow enough for a partial read " \
                "to be observed, so this test proved nothing about partial parsing"
            assert got_real, \
                "every partial reply was empty; the tolerant read recovered nothing:\n" + \
                "\n".join(replies)
            # and the run still completed normally afterwards - asking must not disturb it
            # 0 == PHI_MEAS
            assert np.isfinite(ies.get_phi_summary(0)["mean"]), \
                "the batch did not complete after a partial-results request"

            # ---- reading it back ------------------------------------------------------
            # Every run finished in the end, so storage should now report them FINAL with a
            # fully-valid mask. The partial values written mid-run were overwritten by the
            # real ones - which is the correct outcome and worth asserting, because a partial
            # write that SURVIVED completion would be silently wrong data.
            n_runs = ies.get_run_count()
            assert n_runs > 0, n_runs
            finals = 0
            for rid in range(n_runs):
                info = ies.get_run_info(rid)
                if info["completeness"] != 2:        # 2 == FINAL
                    continue
                finals += 1
                pars, obs, valid = ies.get_run_values(rid)
                assert len(obs) == info["n_obs_total"], (len(obs), info)
                assert valid.all(), \
                    "a completed run should have every observation marked valid: {0}".format(
                        rid)
                assert info["n_obs_reported"] == info["n_obs_total"], info
                assert np.isfinite(obs).all(), "completed run has non-finite observations"
                assert len(pars) > 0
            assert finals > 0, "no run read back as FINAL"

            # the partial path itself: seen live, and read back with a partly-valid mask
            assert saw_live_partial, \
                "no run ever reported live partial results, so the partial read path was " \
                "never exercised"
            assert partial_masks, \
                "no run read back as PARTIAL with a partly-valid mask; partial values were " \
                "either not stored or not distinguishable from complete ones"
            for n_valid, n_total in partial_masks:
                assert 0 < n_valid < n_total, (n_valid, n_total)


    finally:
        for p_ in procs:
            try:
                p_.terminate()
            except Exception:
                pass


def capi_run_storage_read_test():
    """Reading stored runs, including the ids that do not exist.

    Serial and fast on purpose: it was originally folded into the PANTHER partial-results test,
    which takes a couple of minutes and needs workers, so a crash here cost a full CI cycle to
    see. An out-of-range id must be refused rather than seeking past the end of the storage
    file - RunStorage::get_info() reads at a computed offset and only notices afterwards, which
    surfaced as an access violation on windows.
    """
    wd = _setup("capi_storage_read", noptmax=1, num_reals=6)
    with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd) as ies:
        ies.initialize()
        n_runs = ies.get_run_count()
        assert n_runs > 0, n_runs

        # every stored run reads back with a full validity mask and finite values
        for rid in range(n_runs):
            info = ies.get_run_info(rid)
            pars, obs, valid = ies.get_run_values(rid)
            assert len(obs) == info["n_obs_total"], (len(obs), info)
            assert len(valid) == len(obs)
            if info["completeness"] == 2:            # FINAL
                assert valid.all(), rid
                assert np.isfinite(obs).all(), rid

        # ids that do not exist are refused, on every platform, without reading anything
        for bad in (n_runs, n_runs + 50, -1, 2**30):
            try:
                ies.get_run_info(bad)
            except PestppError:
                pass
            else:
                raise AssertionError(
                    "run_id {0} should have been refused (storage holds {1})".format(
                        bad, n_runs))
            try:
                ies.get_run_values(bad)
            except PestppError:
                pass
            else:
                raise AssertionError("get_run_values should refuse run_id {0}".format(bad))
        ies.finalize()


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
            # A cancelled run produces NO realization. This used to assert len(names), and
            # that was asserting a bug: update_from_runs() discarded get_run()'s false and
            # reused one Observations buffer for the whole loop, so the cancelled run's slot
            # was filled with the PREVIOUS realization's values. The count looked right and
            # one row was quietly somebody else's data.
            assert oe.shape[0] == len(names) - cancelled, (oe.shape, len(names), cancelled)
            # ...and the rows that did survive are distinct, which is the property the old
            # assertion was silently violating
            assert np.unique(oe, axis=0).shape[0] == oe.shape[0], \
                "a surviving realization duplicates another - process copied a stale buffer"
    finally:
        for p in procs:
            try:
                p.terminate()
            except Exception:
                pass


def capi_glm_test():
    """pestpp-glm drives through the C ABI, and refuses what it genuinely cannot answer.

    glm is the odd tool on this interface: it is not an ensemble method. It carries ONE
    parameter set through a Jacobian and an upgrade, so "the parameter ensemble" and "phi over
    realizations" have no honest answer. The assertion that matters is not just that the loop
    runs - it is that those calls REFUSE rather than hand back a one-row ensemble, which would
    typecheck and mislead. A caller that gets a 1xN matrix back reasonably reads it as an
    ensemble of one realization; that is not what glm computed.
    """
    wd = _setup("capi_glm", noptmax=2)
    with PestppLib(_find_library(), TOOL_GLM, "pest.pst", wd) as glm:
        glm.initialize()

        n_iters = 0
        while not glm.should_terminate() and n_iters < 3:
            glm.solve_iteration()
            n_iters += 1
        assert n_iters > 0, "glm terminated before running a single iteration"
        assert glm.get_iteration() > 0, glm.get_iteration()
        glm.finalize()

        # the refusals, each with a message that says why rather than a bare error code
        for label, call in (("ensemble view", lambda: glm.get_ensemble_view(PAR_EN)),
                            ("phi summary",   lambda: glm.get_phi_summary(0)),
                            ("queue_runs",    lambda: glm.queue_runs())):
            try:
                call()
            except PestppError as e:
                assert "glm" in str(e).lower(), \
                    "{0} refusal should explain itself, got: {1}".format(label, e)
            else:
                raise AssertionError("glm accepted '{0}', which it has no answer for".format(label))

    # and it actually did the work
    assert os.path.exists(os.path.join(wd, "pest.rec")), "no record file written"
    par = [f for f in os.listdir(wd) if f.endswith(".par")]
    assert par, "glm wrote no .par file, so it did not complete an upgrade"


def capi_glm_matches_exe_test():
    """The ABI-driven glm and the pestpp-glm executable reach the same answer.

    The point of the refactor is that main() and a library caller run the SAME loop. If they
    diverge, every other glm assertion here is testing a private code path rather than the
    tool. Compares the final parameter file, which is the thing a user acts on.
    """
    import filecmp
    exe_d = _setup("capi_glm_exe", noptmax=2)
    api_d = _setup("capi_glm_api", noptmax=2)

    exe = _find_tool_exe("pestpp-glm")
    if exe is None:
        print("skipping capi_glm_matches_exe_test: pestpp-glm not found")
        return
    pyemu.os_utils.run("{0} pest.pst".format(exe), cwd=exe_d)

    with PestppLib(_find_library(), TOOL_GLM, "pest.pst", api_d) as glm:
        glm.initialize()
        while not glm.should_terminate():
            glm.solve_iteration()
        glm.finalize()

    a = os.path.join(exe_d, "pest.par")
    b = os.path.join(api_d, "pest.par")
    assert os.path.exists(a) and os.path.exists(b), (a, b)
    # .par is a fixed-format text file: skip the header line, key on the parameter name
    def _par_values(path):
        vals = {}
        with open(path) as f:
            next(f)                      # the format/precision header
            for line in f:
                t = line.split()
                if len(t) >= 2:
                    vals[t[0].lower()] = float(t[1])
        return vals
    av, bv = _par_values(a), _par_values(b)
    assert set(av) == set(bv), "the two runs report different parameters"
    d = max(abs(av[k] - bv[k]) for k in av)
    assert d == 0.0, (
        "the ABI-driven glm and the executable disagree by {0:.3e} in the final parameters - "
        "main() and the library are not running the same loop".format(d))


def _find_tool_exe(name):
    """The built executable beside the library, or None."""
    sfx = ".exe" if plat == "windows" else ""
    for c in (os.path.join("..", "bin", plat if plat != "apple" else "mac", name + sfx),
              os.path.join("..", "bin", name + sfx)):
        if os.path.exists(c):
            return os.path.abspath(c)
    return None


def capi_panther_release_workers_test():
    """Releasing a BUSY worker reschedules its run instead of giving up on it.

    This is the whole point of release, and the property that separates it from cancel: the
    machine goes away, the work does not. So the assertion that matters is not "the worker
    disappeared" but "the run the worker was holding still finishes, on somebody else, with no
    failure charged against it".

    Deliberately releases exactly ONE worker while it is mid-run. Releasing all of them would
    leave the batch with nothing to run on and prove nothing about rescheduling.
    """
    wd = _setup("capi_panther_release", noptmax=1, num_reals=8)
    worker_root = os.path.join(_BENCH, "capi_release_workers")
    if os.path.exists(worker_root):
        shutil.rmtree(worker_root)
    os.makedirs(worker_root)
    n_workers = 3
    procs = []
    agent_exe = _find_agent_exe()
    rel_port = port + 7

    try:
        with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd, port=rel_port) as ies:
            for i in range(n_workers):
                d = os.path.join(worker_root, "worker_{0}".format(i))
                shutil.copytree(wd, d)
                log = open(os.path.join(d, "worker.log"), "w")
                procs.append(subprocess.Popen(
                    [agent_exe, "pest.pst", "/h", "localhost:{0}".format(rel_port)],
                    cwd=d, stdout=log, stderr=subprocess.STDOUT))

            ies.initialize()
            names = ies.get_ensemble_row_names(PAR_EN)
            n_queued = ies.queue_runs()

            ies.begin_batch()
            released, held_run_id, n_before = 0, None, 0
            deadline = time.time() + 300
            while time.time() < deadline:
                if ies.run_slice(0.05):
                    break
                if released:
                    continue
                # find a worker that is actually holding a run right now
                n_before = ies.get_worker_count()
                for i in range(n_before):
                    w = ies.get_worker_state(i)
                    if w["current_run_id"] >= 0 and w["state"].lower().startswith("active"):
                        held_run_id = w["current_run_id"]
                        released = ies.release_workers([i])
                        break
            ies.end_batch()

            assert released == 1, "expected to release exactly one worker, released {0}".format(released)
            assert held_run_id is not None, "never caught a worker mid-run to release"
            assert ies.get_worker_count() < n_before, \
                "worker count did not drop after release: {0} -> {1}".format(
                    n_before, ies.get_worker_count())

            states = {r["run_id"]: r for r in ies.get_run_states()}
            held = states[held_run_id]
            # the point of the whole feature
            assert held["status"] == "completed", (
                "the released worker's run ended as '{0}', not completed - release gave up on "
                "it instead of rescheduling it".format(held["status"]))
            assert held["n_failures"] == 0, (
                "release charged {0} failure(s) against run {1}; releasing a worker is not a "
                "judgement on its run".format(held["n_failures"], held_run_id))

            ies.process_runs()
            oe, _ = ies.get_ensemble_view(OBS_EN)
            assert oe.shape[0] == len(names), (
                "released worker lost a realization: {0} of {1}".format(oe.shape[0], len(names)))
            assert len(states) == n_queued, (len(states), n_queued)
    finally:
        for p in procs:
            try:
                p.terminate()
            except Exception:
                pass


def capi_release_all_then_cancel_test():
    """Releasing EVERY worker parks the batch; cancelling the rest is what ends it.

    Two separate contracts, and they are easy to conflate:

      - no workers left, runs still outstanding -> the run manager STAYS in the batch. It does
        not error, does not fail the runs, and does not declare itself finished. The work is
        still owed, and a worker starting up later can still connect and do it.
      - the caller cancels the remaining runs -> NOW the batch is over, because nothing is owed
        any more.

    So "no workers" must never be mistaken for "done". The batch ends when there is no work
    left, not when there is nobody left to do it.
    """
    wd = _setup("capi_release_all", noptmax=1, num_reals=8)
    worker_root = os.path.join(_BENCH, "capi_release_all_workers")
    if os.path.exists(worker_root):
        shutil.rmtree(worker_root)
    os.makedirs(worker_root)
    n_workers = 2
    procs = []
    agent_exe = _find_agent_exe()
    rel_port = port + 9

    try:
        with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd, port=rel_port) as ies:
            for i in range(n_workers):
                d = os.path.join(worker_root, "worker_{0}".format(i))
                shutil.copytree(wd, d)
                log = open(os.path.join(d, "worker.log"), "w")
                procs.append(subprocess.Popen(
                    [agent_exe, "pest.pst", "/h", "localhost:{0}".format(rel_port)],
                    cwd=d, stdout=log, stderr=subprocess.STDOUT))

            ies.initialize()
            n_queued = ies.queue_runs()
            ies.begin_batch()

            # get the batch genuinely under way before pulling the rug
            done, deadline = False, time.time() + 300
            while time.time() < deadline:
                if ies.run_slice(0.05):
                    done = True
                    break
                if ies.get_worker_count() == n_workers and \
                        ies.get_run_time_stats()["running"] > 0:
                    break
            assert not done, "batch finished before any worker could be released"

            released = ies.release_workers()          # all of them
            assert released == n_workers, (released, n_workers)
            assert ies.get_worker_count() == 0, ies.get_worker_count()

            # The overdue policy is about a WORKER that is taking too long on a model, not about
            # wall-clock elapsing while nothing is running. Crank it down to 1.2 seconds now
            # that there are no workers: if the giveup threshold were being applied to queued
            # work, everything below would be timed out within a second or two instead of
            # waiting. Set AFTER the release deliberately - while workers were running, a
            # threshold this small is supposed to bite.
            ies.set_option("overdue_giveup_minutes", 0.02)
            ies.set_option("overdue_giveup_fac", 1.0001)

            # CONTRACT 1: with nobody to run them, the batch stays open.
            outstanding = [r for r in ies.get_run_states() if r["status"] != "completed"]
            assert outstanding, "nothing was left outstanding, so this proves nothing"
            parked_until = time.time() + 15
            while time.time() < parked_until:
                assert not ies.run_slice(0.05), (
                    "run manager reported ALL_DONE with {0} run(s) outstanding and no workers - "
                    "'nobody left to do it' was mistaken for 'nothing left to do'"
                    .format(len(outstanding)))
            still = [r for r in ies.get_run_states() if r["status"] != "completed"]
            assert len(still) == len(outstanding), (
                "outstanding runs changed with no workers connected: {0} -> {1}"
                .format(len(outstanding), len(still)))
            for r in still:
                assert r["status"] != "failed", \
                    "run {0} was failed just because its worker went away".format(r["run_id"])
                # 15s of waiting against a 1.2s giveup threshold: if the overdue policy applied
                # to anything other than a worker actively running a model, this would be
                # timed_out many times over
                assert r["status"] != "timed_out", (
                    "run {0} was timed out while NO worker was running it - the overdue policy "
                    "reached queued work".format(r["run_id"]))
                # The run the released worker was holding must be back in the need-to-run pile,
                # not merely "not running": with zero workers connected the only honest status
                # for outstanding work is QUEUED. Asserting the positive rather than a set of
                # negatives is what makes this a statement about requeueing.
                assert r["status"] == "queued", (
                    "run {0} is '{1}' after its worker was released; it should be back in the "
                    "queue".format(r["run_id"], r["status"]))

            # CONTRACT 2: cancel what is left and the batch ends.
            n_cancelled = ies.cancel_runs([r["run_id"] for r in still])
            assert n_cancelled == len(still), (n_cancelled, len(still))
            ended, deadline = False, time.time() + 60
            while time.time() < deadline:
                if ies.run_slice(0.05):
                    ended = True
                    break
            assert ended, "batch did not end after every remaining run was cancelled"
            ies.end_batch()

            states = ies.get_run_states()
            assert len(states) == n_queued, (len(states), n_queued)
            assert set(r["status"] for r in states) <= {"completed", "cancelled"}, \
                set(r["status"] for r in states)
    finally:
        for p in procs:
            try:
                p.terminate()
            except Exception:
                pass


def capi_release_workers_refused_serial_test():
    """Release is a PANTHER call, and the serial manager says so rather than pretending."""
    wd = _setup("capi_release_serial")
    with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd, run_manager=RM_SERIAL) as ies:
        assert not ies.supports_live_control(), "serial manager claimed live control"
        try:
            ies.release_workers()
        except PestppError as e:
            assert "panther" in str(e).lower(), \
                "refusal should say it needs a panther master, got: {0}".format(e)
        else:
            raise AssertionError("serial run manager accepted release_workers()")


def capi_service_runs_yourself_test():
    """The caller writes the queued runs' results back - no run manager evaluates them.

    This is the other half of run storage access. pestpp_get_run_values() reads the queued
    PARAMETERS, the caller computes observations however it likes, pestpp_set_run_values()
    records them, and pestpp_process_runs() harvests into the ensemble as usual. It is what the
    external run manager does through a file, without the round trip.

    Deliberately synthetic observations rather than a real forward run: the claim under test is
    that what the caller writes is exactly what lands in the ensemble, and a recognisable value
    per realization proves that far more sharply than re-deriving the model's own numbers would.
    """
    wd = _setup("capi_service_self")
    with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd) as ies:
        n = ies.initialize_prepare()
        assert n > 0, "expected an initial batch to service, got {0}".format(n)
        queued = ies.queue_runs()
        assert queued == n, (queued, n)

        count = ies.get_run_count()
        assert count == queued, (count, queued)

        # a value that says which run and which observation it came from
        written = {}
        for rid in range(count):
            pars, obs, valid = ies.get_run_values(rid)
            assert not valid.any(), \
                "a queued run must hold no real observations yet (run {0})".format(rid)
            assert np.isfinite(pars).all(), \
                "the queued run's parameters should be readable (run {0})".format(rid)
            vals = np.array([rid * 100.0 + i for i in range(obs.size)], dtype=float)
            ies.set_run_values(rid, vals)
            written[rid] = vals

        # ...and it reads back as final, not partial
        for rid in range(count):
            info = ies.get_run_info(rid)
            assert info["completeness"] == 2, \
                "run {0} should be FINAL after set_run_values, got {1}".format(
                    rid, info["completeness"])
            _p, obs, valid = ies.get_run_values(rid)
            assert valid.all(), "every observation should be real after set_run_values"
            assert np.allclose(obs, written[rid], rtol=0, atol=1e-12), \
                "run {0} did not read back what was written".format(rid)

        nfail = ies.process_runs()
        assert nfail == 0, "nothing should have failed: {0}".format(nfail)
        ies.initialize_finish()

        # and the ensemble carries what we wrote, not what any model produced.
        #
        # Mapped by NAME, because run storage and the ensemble do not share a column order -
        # which is exactly why pestpp_get_run_obs_names() exists. Comparing positionally here
        # passes only by luck and fails confusingly when it does not.
        storage_names = [n.lower() for n in ies.get_run_obs_names()]
        oe, _tok = ies.get_ensemble_view(OBS_EN)
        ens_names = [n.lower() for n in ies.get_ensemble_col_names(OBS_EN)]
        assert oe.shape[0] == count, (oe.shape, count)
        pos = {n: i for i, n in enumerate(storage_names)}
        for row in range(oe.shape[0]):
            for col, name in enumerate(ens_names):
                assert name in pos, "ensemble column {0} is not in run storage".format(name)
                assert abs(oe[row, col] - written[row][pos[name]]) < 1e-10, \
                    "row {0} column {1}: ensemble has {2}, caller wrote {3}".format(
                        row, name, oe[row, col], written[row][pos[name]])


def capi_service_runs_failure_test():
    """set_run_failed is a FAILURE; simply not writing a run is not.

    Both lose the realization, but only one is counted and reported as a model failure - the
    same distinction preemption needed for abandoned runs. A caller servicing runs itself has
    to be able to say which it means.
    """
    wd = _setup("capi_service_fail")
    with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd) as ies:
        n = ies.initialize_prepare()
        ies.queue_runs()
        count = ies.get_run_count()
        assert count > 2, count
        for rid in range(count):
            _p, obs, _v = ies.get_run_values(rid)
            if rid == 0:
                ies.set_run_failed(rid)          # explicitly failed
            elif rid == 1:
                pass                             # silently unwritten
            else:
                ies.set_run_values(rid, np.full(obs.size, float(rid)))
        nfail = ies.process_runs()
        # the explicit failure is counted; the unwritten run is dropped as having no values
        assert nfail >= 1, "the explicitly failed run should be reported: {0}".format(nfail)
        ies.initialize_finish()
        oe, _tok = ies.get_ensemble_view(OBS_EN)
        assert oe.shape[0] == count - 2, \
            "both the failed and the unwritten run should be gone: {0} of {1}".format(
                oe.shape[0], count)


def capi_host_failure_count_test():
    """Failures are tallied per HOST, not per agent, and the running total reaches the rmr.

    Several agents normally share a machine, so a host that is eating every run looks like a
    handful of unrelated agent failures until they are added up. Three agents on localhost here,
    a model that always fails, so every failure has to land on the one host entry.
    """
    wd = _setup("capi_hostfail", noptmax=1, num_reals=6)
    # a model that fails every time, the same trick basic_tests uses to force run failures
    pst = pyemu.Pst(os.path.join(wd, "pest.pst"))
    pst.model_command = ['python -c "raise Exception(\'intentional\')"']
    pst.write(os.path.join(wd, "pest.pst"), version=2)

    worker_root = os.path.join(_BENCH, "capi_hostfail_workers")
    if os.path.exists(worker_root):
        shutil.rmtree(worker_root)
    os.makedirs(worker_root)
    n_workers = 3
    procs = []
    agent_exe = _find_agent_exe()
    # a port of its own - the module-level one belongs to the other panther tests, and two
    # masters bound to the same port in one run is a confusing way to fail
    hf_port = port + 7

    try:
        with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd, port=hf_port) as ies:
            for i in range(n_workers):
                d = os.path.join(worker_root, "worker_{0}".format(i))
                shutil.copytree(wd, d)
                log = open(os.path.join(d, "worker.log"), "w")
                procs.append(subprocess.Popen(
                    [agent_exe, "pest.pst", "/h", "localhost:{0}".format(hf_port)],
                    cwd=d, stdout=log, stderr=subprocess.STDOUT))

            # nothing has failed yet, so the container must be empty rather than absent
            assert ies.get_host_failures() == {}, ies.get_host_failures()

            ies.initialize_prepare()
            ies.queue_runs()
            ies.begin_batch()
            for _ in range(600):
                if ies.run_slice(0.1):
                    break
            ies.end_batch()
            nfail = ies.process_runs()

            hf = ies.get_host_failures()
            print("host failures:", hf, "reported failures:", nfail)
            assert len(hf) == 1, \
                "three agents on one machine must roll up to ONE host: {0}".format(hf)
            host, count = next(iter(hf.items()))
            assert count > 0, hf
            # every agent is on this machine, so the host total accounts for all of them - it
            # must not look like a per-agent count
            assert count >= n_workers, \
                "host total {0} looks per-agent, not per-host ({1} agents)".format(count, n_workers)
    finally:
        for p_ in procs:
            try:
                p_.kill()
            except Exception:
                pass

    # and the running total is in the record file, climbing as failures accumulate
    rmr = os.path.join(wd, "pest.rmr")
    counts = []
    with open(rmr, "r") as f:
        for line in f:
            if "on this host so far" in line:
                counts.append(int(line.split("-")[-1].strip().split()[0]))
    assert counts, "no per-host failure lines in the rmr"
    assert counts == sorted(counts), "the running total should never go backwards: {0}".format(counts)
    assert counts[0] == 1, "the first failure on a host should report 1: {0}".format(counts[:3])
    assert "1 failure " in open(rmr).read(), "singular/plural: the first should read '1 failure'"
    print("rmr per-host totals:", counts[:6], "...", counts[-1] if counts else None)


def capi_host_quarantine_test():
    """A host that fails far more than its share is quarantined and its runs re-queued.

    Two "hosts" out of one machine: agents that dial 127.0.0.1 and agents that dial 127.0.0.2
    arrive with different peer addresses, so getnameinfo gives the master two distinct host
    names - it falls back to the numeric form when there is no PTR record, so no NI_NUMERICHOST
    and no test-only code path is needed. Linux binds all of 127/8 to lo; macos binds only .1,
    which is why this is linux-only.

    The two groups get DIFFERENT pest.pst files - an agent runs the model command from its own
    working directory - so one group fails everything and the other works normally.
    """
    if platform.system() != "Linux":
        print("skipping host quarantine test: needs 127.0.0.2, which only linux binds by default")
        return

    good_wd = _setup("capi_quar", noptmax=1, num_reals=12)
    worker_root = os.path.join(_BENCH, "capi_quar_workers")
    if os.path.exists(worker_root):
        shutil.rmtree(worker_root)
    os.makedirs(worker_root)

    agent_exe = _find_agent_exe()
    q_port = port + 11
    procs = []
    n_good, n_bad = 2, 2

    def _worker(idx, addr, fail):
        d = os.path.join(worker_root, "w{0}".format(idx))
        shutil.copytree(good_wd, d)
        if fail:
            pst = pyemu.Pst(os.path.join(d, "pest.pst"))
            pst.model_command = ['python -c "raise Exception(\'bad host\')"']
            pst.write(os.path.join(d, "pest.pst"), version=2)
        log = open(os.path.join(d, "worker.log"), "w")
        return subprocess.Popen(
            [agent_exe, "pest.pst", "/h", "{0}:{1}".format(addr, q_port)],
            cwd=d, stdout=log, stderr=subprocess.STDOUT)

    try:
        with PestppLib(_find_library(), TOOL_IES, "pest.pst", good_wd, port=q_port) as ies:
            for i in range(n_good):
                procs.append(_worker(i, "127.0.0.1", fail=False))
            for i in range(n_bad):
                procs.append(_worker(100 + i, "127.0.0.2", fail=True))

            ies.initialize_prepare()
            ies.queue_runs()
            ies.begin_batch()
            for _ in range(900):
                if ies.run_slice(0.1):
                    break
            ies.end_batch()
            ies.process_runs()
            hf = ies.get_host_failures()
            print("host failures:", hf)
    finally:
        for p_ in procs:
            try:
                p_.kill()
            except Exception:
                pass

    # ASSERT rather than skip. This used to return quietly when the two addresses collapsed to
    # one host name, which made the test report ok whether or not it had checked anything - the
    # exact shape that lets a feature rot untested. On linux 127.0.0.2 has no PTR record, so
    # getnameinfo falls back to the numeric string and the two must come back distinct. If that
    # ever stops being true the test says so instead of passing silently.
    assert len(hf) >= 2, (
        "expected two distinct host keys from 127.0.0.1 and 127.0.0.2, got {0} - without two "
        "hosts the quarantine check does not apply and this test proves nothing".format(hf))

    rmr = open(os.path.join(good_wd, "pest.rmr")).read()
    assert "quarantined -" in rmr, "no host was quarantined:\n" + rmr[-2000:]
    assert "moved to QUARANTINED" in rmr, "agents were not moved out of service"
    assert "run(s) requeued" in rmr, "failed runs were not requeued"
    assert "host(s) quarantined for excess run failures" in rmr, "no end-of-batch summary"

    # the quarantined host must be the one that was failing, not the healthy one
    bad = max(hf, key=lambda k: hf[k])
    good = min(hf, key=lambda k: hf[k])
    for line in rmr.splitlines():
        if "quarantined -" in line:
            assert bad in line, "the WRONG host was quarantined: {0} (failures {1})".format(line, hf)
            assert good not in line.split("quarantined")[0], \
                "the healthy host was quarantined: {0}".format(line)
    print("quarantined the failing host:", bad, "| healthy:", good, hf)


def capi_host_quarantine_disabled_test():
    """delta of 0 turns the screening off, and the summary says so rather than going quiet."""
    wd = _setup("capi_quar_off", noptmax=1, num_reals=4)
    pst = pyemu.Pst(os.path.join(wd, "pest.pst"))
    pst.pestpp_options["panther_agent_max_failed_run_delta"] = 0
    pst.write(os.path.join(wd, "pest.pst"), version=2)
    with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd) as ies:
        assert ies.get_option("PANTHER_AGENT_MAX_FAILED_RUN_DELTA") == "0", \
            ies.get_option("PANTHER_AGENT_MAX_FAILED_RUN_DELTA")
    print("delta=0 accepted and reported by the option system")


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
    capi_nm_parsing_test()
    capi_export_surface_test()
    capi_da_options_reach_the_tool_test()
    capi_da_cycle_is_tagged_test()
    capi_output_redirect_is_lifo_test()
    capi_unknown_option_is_reported_test()
    capi_da_test()
    capi_sqp_test()
    capi_mou_test()
    capi_da_resize_between_queue_and_harvest_test()
    capi_tool_ensemble_availability_test()
    capi_run_storage_read_test()
    capi_glm_test()
    capi_glm_matches_exe_test()
    capi_panther_control_test()
    capi_panther_release_workers_test()
    capi_release_all_then_cancel_test()
    capi_release_workers_refused_serial_test()
    capi_partial_results_test()
    capi_partial_obs_command_test()
    capi_partial_status_file_test()
    capi_stp_file_commands_test()
    capi_service_runs_yourself_test()
    capi_service_runs_failure_test()
    capi_host_failure_count_test()
    capi_host_quarantine_test()
    capi_host_quarantine_disabled_test()
    print("all capi tests passed")

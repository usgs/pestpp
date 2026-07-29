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
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_REPO, "python"))
from pestpp_lib import (  # noqa: E402
    PestppLib, PestppError, PAR_EN, OBS_EN, TOOL_IES, WORKER_COMPLETED,
)

plat = "unknown"
if "linux" in platform.platform().lower():
    plat = "linux"
elif "darwin" in platform.platform().lower() or "macos" in platform.platform().lower():
    plat = "apple"
else:
    plat = "windows"

# the model binary (mfnwt) has to be findable, same as the other benchmark scripts
os.environ["PATH"] += os.pathsep + os.path.abspath("test_bin")

exe = ".exe" if plat == "windows" else ""
_sub = {"windows": "win", "apple": "mac", "linux": "linux"}[plat]


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
    names = {"windows": "pestpp_capi.dll",
             "apple": "libpestpp_capi.dylib",
             "linux": "libpestpp_capi.so"}
    want = names[plat]
    roots = [os.path.join(_REPO, "build"), os.path.join("..", "..", "pestpp", "build"),
             os.path.join("..", "..", "..", "..", "pestpp", "build")]
    for root in roots:
        if not os.path.isdir(root):
            continue
        for dirpath, _, filenames in os.walk(root):
            if want in filenames:
                return os.path.abspath(os.path.join(dirpath, want))
    raise RuntimeError(
        "could not find {0} under any of {1}. The C ABI is built by the normal cmake build "
        "(always SHARED); if this fires, the build did not produce it.".format(want, roots))


def _setup(test_d, noptmax=1, num_reals=6):
    """A small, fast ies case in its own directory."""
    base_d = os.path.join("ies_10par_xsec", "template")
    if os.path.exists(test_d):
        shutil.rmtree(test_d)
    shutil.copytree(base_d, test_d)
    pst = pyemu.Pst(os.path.join(test_d, "pest.pst"))
    pst.control_data.noptmax = noptmax
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["random_seed"] = 11
    pst.write(os.path.join(test_d, "pest.pst"), version=2)
    return os.path.abspath(test_d)


def capi_smoke_test():
    """The library loads, a handle runs iterations, and views are live windows on the data."""
    wd = _setup("capi_smoke", noptmax=2)
    with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd) as ies:
        assert ies.get_option("IES_NUM_REALS") == "6", ies.get_option("IES_NUM_REALS")

        ies.initialize()
        assert ies.get_iteration() == 0, ies.get_iteration()

        pe = ies.get_ensemble_view(PAR_EN)
        rnames = ies.get_ensemble_row_names(PAR_EN)
        cnames = ies.get_ensemble_col_names(PAR_EN)
        assert pe.shape == (len(rnames), len(cnames)), (pe.shape, len(rnames), len(cnames))
        # Eigen is column-major; a view that is not F-contiguous means the numpy wrapper is
        # transposing silently, which would corrupt every value a caller writes back
        assert pe.flags["F_CONTIGUOUS"], "ensemble view lost its column-major layout"

        # it is a window, not a copy: a write through numpy must be visible to C++
        before = pe[0, 0]
        pe[0, 0] = before + 1234.5
        assert ies.get_ensemble_view(PAR_EN)[0, 0] == before + 1234.5, "view is a copy, not a view"
        pe[0, 0] = before

        for _ in range(2):
            ies.solve_iteration()
        assert ies.get_iteration() == 2, ies.get_iteration()
        oe = ies.get_ensemble_view(OBS_EN)
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


def _drive_batch(ies):
    """Queue, run and harvest, with the caller owning the run loop."""
    n_queued = ies.queue_ensemble()
    ies.begin_batch()
    while not ies.run_slice(0.05):
        pass
    ies.end_batch()
    return n_queued


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
        assert ies.harvest_ensemble() == 0, "reference pass had failed runs"
        oe_names = ies.get_ensemble_row_names(OBS_EN)
        ref_oe = ies.get_ensemble_view(OBS_EN)
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

        n_failed = ies.harvest_ensemble()
        assert n_failed == 0, "{0} runs failed unexpectedly".format(n_failed)

        oe = ies.get_ensemble_view(OBS_EN)
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
            ies.harvest_ensemble()
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
        ies.queue_ensemble()
        try:
            ies.queue_ensemble()
            raise AssertionError("double queue should raise")
        except PestppError as e:
            assert "already queued" in str(e), str(e)


def capi_panther_control_test():
    """Watch and cancel runs mid-batch against a real PANTHER master.

    The serial manager cannot yield, so this is the only path where run states, worker stats
    and cancel mean anything.
    """
    wd = _setup("capi_panther", noptmax=1, num_reals=8)
    worker_root = os.path.abspath("capi_panther_workers")
    if os.path.exists(worker_root):
        shutil.rmtree(worker_root)
    os.makedirs(worker_root)
    n_workers = 3
    procs = []
    agent_exe = _find_agent_exe()

    try:
        with PestppLib(_find_library(), TOOL_IES, "pest.pst", wd, port=port) as ies:
            assert ies.supports_live_control(), "panther master denied live control"

            # workers must be up before initialize(): it evaluates the prior ensemble, and
            # with no workers the master waits for them forever
            for i in range(n_workers):
                d = os.path.join(worker_root, "worker_{0}".format(i))
                shutil.copytree(wd, d)
                log = open(os.path.join(d, "worker.log"), "w")
                procs.append(subprocess.Popen(
                    [agent_exe, "pest.pst", "/h", "localhost:{0}".format(port)],
                    cwd=d, stdout=log, stderr=subprocess.STDOUT))

            ies.initialize()
            names = ies.get_ensemble_row_names(PAR_EN)
            n_queued = ies.queue_ensemble()
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

            ies.harvest_ensemble()
            oe = ies.get_ensemble_view(OBS_EN)
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
    capi_panther_control_test()
    print("all capi tests passed")

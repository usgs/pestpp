"""Tests for the ergonomic layer, python/pestpp.py.

capi_tests.py covers the C ABI and the thin binding. This file covers the layer users
actually touch: DataFrames, properties, the iteration generator, view scoping, and the
one-call run_* functions.

Worth testing separately because the helper layer can drift from the binding underneath it
without anything failing to compile - a renamed C symbol shows up here as an AttributeError,
and a semantic change (an iteration count, a transform space) shows up nowhere else at all.
"""
import gc
import os
import platform
import shutil
import sys

import numpy as np
import pandas as pd
import pyemu

_BENCH = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_BENCH)
sys.path.insert(0, os.path.join(_REPO, "python"))
from pestpp_lib import NOISE_EN  # noqa: E402
from pestpp import (  # noqa: E402
    Ies, Mou, Sqp, IterationStep, PestppError, ExpiredViewError, find_library,
    run_ies, PHI_ACTUAL,
)

plat = platform.platform().lower()
_sub = "win" if os.name == "nt" else ("mac" if ("darwin" in plat or "macos" in plat) else "linux")
# absolute, and anchored to this file: the library chdirs while a session is open
os.environ["PATH"] += os.pathsep + os.path.join(_BENCH, "test_bin", _sub)


def _case(name, noptmax=2, num_reals=6, **options):
    """A small ies case in its own directory."""
    base = os.path.join(_BENCH, "ies_10par_xsec", "template")
    d = os.path.join(_BENCH, name)
    if os.path.exists(d):
        shutil.rmtree(d)
    shutil.copytree(base, d)
    pst = pyemu.Pst(os.path.join(d, "pest.pst"))
    pst.control_data.noptmax = noptmax
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["random_seed"] = 11
    for k, v in options.items():
        pst.pestpp_options[k] = v
    pst.write(os.path.join(d, "pest.pst"), version=2)
    return d


def api_smoke_test():
    """A session opens, reports itself sensibly, and produces labelled DataFrames."""
    wd = _case("api_smoke", noptmax=1)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        assert ies.run_manager == "serial", ies.run_manager
        assert not ies.supports_live_control
        assert ies.noptmax == 1, ies.noptmax
        assert not ies.is_initialized

        ies.initialize()
        assert ies.is_initialized
        assert ies.n_reals == 6, ies.n_reals
        assert ies.iteration == 0, ies.iteration
        assert np.isfinite(ies.phi) and ies.phi > 0, ies.phi

        par, obs = ies.par_df(), ies.obs_df()
        assert isinstance(par, pd.DataFrame) and isinstance(obs, pd.DataFrame)
        assert par.shape[0] == ies.n_reals and obs.shape[0] == ies.n_reals
        # the labels have to be real, not positional: this is what makes the frames usable
        assert list(par.index) == list(ies.real_names), "par_df index is not the realizations"
        assert list(obs.columns) == ies.obs_names, "obs_df columns are not the observations"
        ies.finalize()


def api_iterations_respect_noptmax_test():
    """The generator stops at noptmax, not just when phi flattens.

    should_terminate() reports only the phi-based criteria, so a loop driven by it alone runs
    until phi stops improving - which is emphatically not what noptmax=2 asked for. This is a
    regression guard for exactly that.
    """
    wd = _case("api_noptmax", noptmax=2)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        assert ies.noptmax == 2, ies.noptmax
        ies.initialize()
        steps = list(ies.iterations())
        assert len(steps) == 2, "ran {0} iterations for noptmax=2".format(len(steps))
        assert all(isinstance(s, IterationStep) for s in steps)
        assert [s.iter for s in steps] == [1, 2], [s.iter for s in steps]
        # phi should have moved; these are the numbers a user reads off
        assert steps[-1].phi_mean < steps[0].phi_mean, [s.phi_mean for s in steps]
        assert np.isfinite(steps[0].phi_std)
        ies.finalize()

    # max_iter overrides it, for a caller who wants a couple of steps and a look around
    wd = _case("api_noptmax_override", noptmax=5)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        assert len(list(ies.iterations(max_iter=1))) == 1
        ies.finalize()


def api_par_df_roundtrip_test():
    """set_par_df matches by name, so a reordered frame still lands correctly."""
    wd = _case("api_roundtrip", noptmax=1)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        df = ies.par_df()

        ies.set_par_df(df)
        assert np.allclose(ies.par_df().values, df.values), "round trip changed values"

        # hand it back reversed: matched by name, so the values must land in the same rows
        ies.set_par_df(df.iloc[::-1])
        back = ies.par_df()
        assert np.allclose(back.values, df.values), "set_par_df matched by position, not name"
        assert list(back.index) == list(df.index)
        ies.finalize()


def api_par_view_is_live_test():
    """par_view yields the tool's own buffer, and only inside the block."""
    wd = _case("api_view", noptmax=1)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        assert ies.par_transform in ("ctl", "num", "model"), ies.par_transform

        with ies.par_view() as arr:
            assert arr.flags["F_CONTIGUOUS"], "view lost Eigen's column-major layout"
            assert arr.shape[0] == ies.n_reals
            before = arr[0, 0]
            arr[0, 0] = before + 321.0

        # the title of this test claims "only inside the block", so prove it. The array is the
        # tool's live Eigen buffer; using it afterwards is a use-after-free that reads as data.
        assert not arr.valid, "the view still reports valid after its block ended"
        for attempt in (lambda: arr[0, 0],
                        lambda: arr.__setitem__((0, 0), 1.0),
                        lambda: arr.shape,
                        lambda: np.asarray(arr)):
            try:
                attempt()
                raise AssertionError(
                    "an expired view was still usable - the context manager is decorative")
            except ExpiredViewError:
                pass

        with ies.par_view() as arr:
            assert np.isclose(arr[0, 0], before + 321.0), "write through the view did not stick"
            arr[0, 0] = before
        ies.finalize()


def api_view_invalidated_by_resize_test():
    """A view goes invalid when the ensemble reallocates, not just when the block ends.

    This is the harder half: the block is still open, the array still has its shape, and the
    numbers it returns come from memory the tool has moved on from. Nothing about the array
    itself can tell the caller that, which is why the library is asked instead.
    """
    wd = _case("api_view_resize", noptmax=1)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        names = ies.real_names
        assert len(names) > 2, "need enough realizations to drop one"

        with ies.par_view() as arr:
            assert arr.valid
            ies.drop_realizations([names[0]])       # reallocates the ensemble underneath
            assert not arr.valid, \
                "the view claims to be valid after a resize; it is pointing at freed memory"
            try:
                arr[0, 0]
                raise AssertionError("a stale view was readable after a resize")
            except ExpiredViewError:
                pass
        ies.finalize()


def api_own_the_initial_batch_test():
    """The prior ensemble can be replaced before it is ever evaluated.

    Every realization is set identical, so every observation row must come back identical.
    With the drawn ensemble they differ substantially, so this cannot pass by accident.
    """
    wd = _case("api_own_batch", noptmax=1)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        n = ies.initialize(defer_runs=True)
        assert n == ies.n_reals and n > 1, (n, ies.n_reals)
        assert not ies.is_initialized, "deferred initialize should not report initialized"

        drawn = ies.par_df()
        mine = drawn.copy()
        mine.loc[:, :] = drawn.iloc[0].values
        ies.set_par_df(mine)

        assert ies.run_ensemble() == 0, "runs failed on the substituted ensemble"
        ies.finish_initialize()
        assert ies.is_initialized

        obs = ies.obs_df().values
        spread = np.abs(obs - obs[0]).max()
        assert spread < 1.0e-8, (
            "identical parameters gave different observations (spread {0}); the substituted "
            "ensemble was not the one evaluated".format(spread))
        ies.finalize()


def api_phi_across_realizations_test():
    """phi is available per realization and per type, not just as one number."""
    wd = _case("api_phi", noptmax=1, ies_no_noise=False)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()

        df = ies.phi_df()
        assert list(df.columns) == ["meas", "actual", "composite", "regul", "noise"], list(df.columns)
        assert len(df) == ies.n_reals, (len(df), ies.n_reals)
        # the index is realizations, and meas/actual are real numbers that differ per row
        assert df["meas"].notna().all() and df["actual"].notna().all()
        assert df["meas"].std() > 0, "every realization has the same phi?"
        # a type not in play is NaN, not missing and not a sentinel
        assert df["regul"].isna().all(), "regul should be absent for this case"

        # the scalar agrees with the vector
        assert np.isclose(ies.phi, df["meas"].mean()), (ies.phi, df["meas"].mean())

        summ = ies.phi_summary_df()
        assert list(summ.columns) == ["mean", "std", "min", "max"], list(summ.columns)
        assert np.isclose(summ.loc["meas", "min"], df["meas"].min())
        # absent types must be NaN rather than the handler's DBL_MAX/-1e30 sentinels
        assert summ.loc["regul"].isna().all(), summ.loc["regul"].to_dict()

        one = ies.phi_df(phi_type=PHI_ACTUAL)
        assert list(one.columns) == ["actual"], list(one.columns)
        ies.finalize()


def api_phi_by_obs_group_test():
    """Phi decomposes by observation group, and the pieces add up.

    The case is rewritten to put its two non-zero-weighted observations in different groups,
    because a single-group case would not actually exercise the grouping.
    """
    base = os.path.join(_BENCH, "ies_10par_xsec", "template")
    d = os.path.join(_BENCH, "api_phigrp")
    if os.path.exists(d):
        shutil.rmtree(d)
    shutil.copytree(base, d)
    pst = pyemu.Pst(os.path.join(d, "pest.pst"))
    nz = pst.nnz_obs_names
    assert len(nz) >= 2, nz
    # pest++ upper-cases names internally; pyemu keeps them lower. Anything joining the two
    # has to reconcile that, and this test is the first place it bites.
    NZ = [n.upper() for n in nz]
    pst.observation_data.loc[nz[0], "obgnme"] = "grp_a"
    pst.observation_data.loc[nz[1], "obgnme"] = "grp_b"
    pst.control_data.noptmax = 1
    pst.pestpp_options["ies_num_reals"] = 6
    pst.pestpp_options["random_seed"] = 11
    pst.write(os.path.join(d, "pest.pst"), version=2)

    with Ies.from_pst("pest.pst", workdir=d) as ies:
        ies.initialize()

        groups = ies.obs_groups
        assert set(groups.loc[NZ]) == {"GRP_A", "GRP_B"}, groups.loc[NZ].to_dict()

        gdf = ies.phi_group_df()
        assert sorted(gdf.columns) == ["GRP_A", "GRP_B"], list(gdf.columns)
        assert len(gdf) == ies.n_reals, (len(gdf), ies.n_reals)

        odf = ies.phi_obs_df()
        assert len(odf) == ies.n_reals
        assert set(odf.columns) <= set(ies.obs_names), "residual columns are not observations"

        # the decomposition has to add up, in both directions
        assert np.allclose(gdf.sum(axis=1), ies.phi_df()["meas"]), \
            "group phi does not sum to the realization phi"
        assert np.allclose(odf.sum(axis=1), gdf.sum(axis=1)), \
            "per-observation terms do not sum to the group totals"

        # and the terms are non-negative, being squared weighted residuals
        assert (odf.values >= 0).all(), "a squared weighted residual came back negative"

        # ACTUAL is a different decomposition of the same shape
        gact = ies.phi_group_df(phi_type=PHI_ACTUAL)
        assert list(gact.columns) == list(gdf.columns)
        assert np.allclose(gact.sum(axis=1), ies.phi_df()["actual"])
        ies.finalize()


def api_weights_test():
    """Weights are changeable both ways, and the change actually reaches phi.

    Two distinct things share the name: the control-file VECTOR (one weight per observation)
    and the weights ENSEMBLE (one per observation per realization). initialize() seeds the
    ensemble from the vector once and phi is computed from the ensemble thereafter, so a
    vector change only takes effect once broadcast - and phi is cached on top of that. Both
    steps are covered here because either one missing looks like "setting weights did
    nothing".
    """
    wd = _case("api_weights", noptmax=1)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()

        w = ies.obs_weights
        assert len(w) == len(ies.obs_names), (len(w), len(ies.obs_names))
        nz = w[w > 0]
        assert len(nz) >= 2, nz.to_dict()
        name, w0 = nz.index[0], float(nz.iloc[0])
        phi0 = ies.phi

        # -- the vector path --
        ies.set_obs_weights({name: w0 * 4.0})
        assert np.isclose(ies.obs_weights[name], w0 * 4.0), ies.obs_weights[name]
        assert ies.phi > phi0, "quadrupling a weight did not raise phi ({0} -> {1})".format(
            phi0, ies.phi)
        # the broadcast reached the ensemble
        assert np.isclose(ies.weights_df().iloc[0][name], w0 * 4.0)

        # zero weights mean zero phi, which is the cleanest possible check that phi is
        # actually being recomputed rather than served from cache
        ies.set_obs_weights(0.0)
        assert np.isclose(ies.phi, 0.0), ies.phi

        # and it round trips exactly
        ies.set_obs_weights(nz)
        assert np.isclose(ies.phi, phi0), (ies.phi, phi0)

        # -- the ensemble path: weights differing between realizations --
        wdf = ies.weights_df()
        assert wdf.shape[0] == ies.n_reals, (wdf.shape, ies.n_reals)
        with ies.weights_view() as arr:
            arr[0, :] = 0.0
        ies.update_phi()
        per_real = ies.phi_df()["meas"]
        assert np.isclose(per_real.iloc[0], 0.0), per_real.iloc[0]
        assert per_real.iloc[1] > 0, "zeroing one realization should not affect the others"
        ies.finalize()


def api_weights_reject_bad_input_test():
    """A bad observation name or a negative weight is refused, not silently applied."""
    wd = _case("api_weights_bad", noptmax=1)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        before = ies.obs_weights.copy()

        try:
            ies.set_obs_weights({"NOT_AN_OBSERVATION": 1.0})
            raise AssertionError("an unknown observation should raise")
        except PestppError as e:
            assert "no such observation" in str(e), str(e)

        name = before[before > 0].index[0]
        try:
            ies.set_obs_weights({name: -1.0})
            raise AssertionError("a negative weight should raise")
        except PestppError as e:
            assert "negative weight" in str(e), str(e)

        # nothing was applied by either failed call
        assert np.allclose(ies.obs_weights.values, before.values), "a rejected call still wrote"
        ies.finalize()


def api_name_case_is_absorbed_test():
    """Names are accepted in any case, and frames can come back pyemu-cased.

    pest++ upper-cases names, pyemu keeps them lower. That mismatch is a constant papercut
    when the two are used together, so the helper layer absorbs it in both directions.
    """
    wd = _case("api_case", noptmax=1)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()

        # the library's own case is upper - if that ever changes this test is the canary
        w = ies.obs_weights
        name = w[w > 0].index[0]
        assert name == name.upper(), name

        # -- lower-case input is accepted everywhere names are taken --
        ies.set_obs_weights({name.lower(): 25.0})
        assert np.isclose(ies.obs_weights[name], 25.0), ies.obs_weights[name]
        assert ies.get_obs_group(name.lower()) == ies.get_obs_group(name)

        reals = list(ies.real_names)
        ies.drop_realizations([reals[1].lower()])
        assert len(ies.real_names) == len(reals) - 1
        assert reals[1] not in ies.real_names

        # -- output can be lower-cased for joining against pyemu --
        up, low = ies.par_df(), ies.par_df(lower=True)
        assert list(low.columns) == [c.lower() for c in up.columns]
        assert list(low.index) == [i.lower() for i in up.index]
        assert np.allclose(low.values, up.values), "lower-casing changed the values"

        # and a lower-cased frame still round trips, because input is normalized
        ies.set_par_df(low)
        assert np.allclose(ies.par_df().values, up.values), \
            "a lower-cased frame did not round trip"

        for f in (ies.obs_df, ies.phi_obs_df, ies.weights_df):
            assert list(f(lower=True).index) == [i.lower() for i in f().index], f.__name__
        ies.finalize()


def api_drop_realizations_test():
    """Dropping applies across the coupled ensembles, and the run continues."""
    wd = _case("api_drop", noptmax=1)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        before = list(ies.real_names)
        victims = before[1:3]

        ies.drop_realizations(victims)
        after = list(ies.real_names)
        assert len(after) == len(before) - 2, (len(after), len(before))
        assert not set(victims) & set(after), "dropped realizations survived"
        assert ies.obs_df().shape[0] == len(after), "obs ensemble was not dropped in step"
        assert ies.par_df().shape[0] == len(after)

        step = ies.solve()          # still usable afterwards
        assert step.n_reals == len(after), (step.n_reals, len(after))
        ies.finalize()


def api_quiet_captures_output_test():
    """The library's console output is captured to a file rather than flooding the session.

    The content check matters more than it looks. On windows the library links the static CRT
    and so buffers its output privately: child processes inherit the redirected descriptor and
    land in the file, but the library's own text escapes to the console unless it is flushed
    before the redirect unwinds. A size check alone passes in that case - the model output is
    there - so this asserts on text the library itself prints.
    """
    wd = _case("api_quiet", noptmax=1)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        assert os.path.exists(ies.log_file), ies.log_file
        assert os.path.getsize(ies.log_file) > 0, "quiet mode captured nothing"
        with open(ies.log_file, errors="replace") as f:
            captured = f.read()
        low = captured.lower()
        # phrases only the library prints, not the model
        assert any(tok in low for tok in ("initializing", "pest++", "run manager")), (
            "the library's own output was not captured - only {0} bytes, starting: {1!r}. "
            "On windows this means pestpp_flush_output() is not reaching the private CRT "
            "buffer before the descriptor is restored.".format(len(captured), captured[:300]))
        ies.finalize()


def api_tools_without_phi_say_so_test():
    """mou and sqp have no phi, and asking says why rather than returning a wrong number."""
    base = os.path.join(_BENCH, "g07", "template")
    d = os.path.join(_BENCH, "api_mou")
    if os.path.exists(d):
        shutil.rmtree(d)
    shutil.copytree(base, d)
    pst = pyemu.Pst(os.path.join(d, "g07.pst"))
    for k in [k for k in pst.pestpp_options if k.startswith("sqp_")]:
        pst.pestpp_options.pop(k)
    pst.pestpp_options["mou_population_size"] = 10
    pst.pestpp_options["mou_generator"] = "de"
    pst.pestpp_options["random_seed"] = 11
    pst.control_data.noptmax = 1
    pst.write(os.path.join(d, "g07.pst"), version=2)

    with Mou.from_pst("g07.pst", workdir=d) as mou:
        mou.initialize()
        assert mou.n_reals == 10, mou.n_reals
        try:
            mou.phi
            raise AssertionError("mou should not report a phi")
        except PestppError as e:
            assert "no phi" in str(e), str(e)
        assert np.isnan(mou.solve().phi_mean), "mou steps should carry no phi"
        mou.finalize()


def api_run_ies_one_call_test():
    """The whole-job function returns the three frames it promises."""
    wd = _case("api_oneshot", noptmax=1)
    out = run_ies("pest.pst", workdir=wd)
    assert set(out) == {"par", "obs", "steps"}, sorted(out)
    assert out["par"].shape[0] == 6 and out["obs"].shape[0] == 6
    assert len(out["steps"]) == 1, len(out["steps"])
    assert out["steps"].iloc[0]["phi_mean"] > 0


def api_find_library_test():
    """The library is discoverable, and PESTPP_LIB wins when set."""
    path = find_library()
    assert os.path.exists(path), path

    saved = os.environ.get("PESTPP_LIB")
    try:
        os.environ["PESTPP_LIB"] = path
        assert find_library() == path
        os.environ["PESTPP_LIB"] = os.path.join(_BENCH, "definitely_not_here.so")
        try:
            find_library()
            raise AssertionError("a bad PESTPP_LIB should raise")
        except FileNotFoundError as e:
            assert "PESTPP_LIB" in str(e), str(e)
    finally:
        if saved is None:
            os.environ.pop("PESTPP_LIB", None)
        else:
            os.environ["PESTPP_LIB"] = saved


def api_notebook_runs_test():
    """The example notebook executes cleanly, start to finish.

    An example that is not run is a claim, not a demonstration - and this one is the first
    thing a new user reads, so a stale cell is worse than no notebook. Executed in place with
    the notebook's own directory as the working directory, because its paths are relative to
    itself.
    """
    import nbformat
    from nbclient import NotebookClient

    # The _flat variants are the same demos written without `with`, and they are covered here
    # for the same reason as the original: they are what a notebook user is pointed at, so a
    # stale cell in one is as bad as a stale cell in the other. They deliberately use their own
    # scratch directories and ports, so running both in one session cannot have them collide.
    for nb_name in ("pestpp_api_demo.ipynb", "pestpp_api_demo_flat.ipynb"):
        nb_path = os.path.join(_REPO, "examples", nb_name)
        assert os.path.exists(nb_path), nb_path
        nb_dir = os.path.dirname(nb_path)

        nb = nbformat.read(nb_path, as_version=4)
        client = NotebookClient(nb, timeout=1800, kernel_name="python3",
                                allow_errors=True,
                                resources={"metadata": {"path": nb_dir}})
        client.execute()

        failures = []
        for i, cell in enumerate(nb.cells):
            if cell.cell_type != "code":
                continue
            for out in cell.get("outputs", []):
                if out.output_type == "error":
                    failures.append("cell {0}: {1}: {2}".format(
                        i, out.ename, str(out.evalue)[:200]))
        assert not failures, "{0} cells raised:\n  ".format(nb_name) + "\n  ".join(failures)

        # scratch directories the notebook creates alongside itself
        for stem in ("ies", "prior", "cull", "parallel", "parallel_workers",
                     "defer", "progress"):
            for prefix in ("nb_", "nbflat_"):
                shutil.rmtree(os.path.join(nb_dir, prefix + stem), ignore_errors=True)


def api_pyemu_objects_test():
    """The pyemu surface hands back real pyemu objects, not lookalikes.

    The test is not that the types are right - it is that the objects WORK: the Pst carries
    live weights, the ensembles compute phi, and the phi they compute agrees with the phi the
    tool reports. A DataFrame dressed up as a ParameterEnsemble would pass a type check and
    fail every one of these.
    """
    wd = _case("api_pyemu", noptmax=1)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()

        pst = ies.pst
        assert isinstance(pst, pyemu.Pst)
        assert pst is ies.pst, "the Pst should be parsed once and cached"

        pe, oe = ies.pe, ies.oe
        assert isinstance(pe, pyemu.ParameterEnsemble), type(pe)
        assert isinstance(oe, pyemu.ObservationEnsemble), type(oe)
        assert pe.shape[0] == ies.n_reals and oe.shape[0] == ies.n_reals

        # lowercase and pyemu-shaped, so these line up with the Pst without any renaming
        assert set(oe.columns) >= set(pst.nnz_obs_names), \
            "obs ensemble columns do not match the Pst's observation names"
        assert set(pe.columns) >= set(pst.adj_par_names), \
            "par ensemble columns do not match the Pst's parameter names"

        # the cross-check that makes the integration worth anything: pyemu computes phi from
        # the ensemble and the Pst's weights, the library computes it internally, and they
        # are the same quantity
        theirs = oe.phi_vector
        mine = ies.phi_df(PHI_ACTUAL, lower=True)["actual"]
        common = [n for n in theirs.index if n in mine.index]
        assert len(common) == ies.n_reals, (len(common), ies.n_reals)
        assert np.allclose(theirs.loc[common].values, mine.loc[common].values, rtol=1e-6), \
            "pyemu's phi_vector disagrees with the library's PHI_ACTUAL:\n{0}".format(
                pd.DataFrame({"pyemu": theirs.loc[common], "api": mine.loc[common]}))

        # axis names, so a merge or reset_index produces something a pyemu user recognises
        assert ies.par_df().index.name == "realization"
        assert ies.obs_df().columns.name == "obsnme"
        ies.finalize()


def api_pyemu_weights_stay_live_test():
    """pst.observation_data.weight tracks weights changed through the API.

    A cached Pst read off disk would go stale the moment a weight changed, and everything
    derived from it - phi_vector most of all - would quietly disagree with the tool.
    """
    wd = _case("api_pyemu_w", noptmax=1)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        target = ies.pst.nnz_obs_names[0]
        before = float(ies.pst.observation_data.loc[target, "weight"])

        ies.set_obs_weights({target: before * 3.0}, broadcast=True)
        after = float(ies.pst.observation_data.loc[target, "weight"])
        assert abs(after - before * 3.0) < 1.0e-9, \
            "the Pst still reports the old weight ({0} vs {1}) - it is a stale copy".format(
                after, before * 3.0)

        # and the ensembles built from it agree with the tool's own recomputed phi
        ies.update_phi()
        theirs = ies.oe.phi_vector
        mine = ies.phi_df(PHI_ACTUAL, lower=True)["actual"]
        common = [n for n in theirs.index if n in mine.index]
        assert np.allclose(theirs.loc[common].values, mine.loc[common].values, rtol=1e-6), \
            "after reweighting, pyemu and the library disagree on phi"
        ies.finalize()


def api_pyemu_bounds_enforced_test():
    """An out-of-bounds parameter is refused rather than silently turned into NaN.

    This is the sharpest consequence of the helper layer not knowing the Pst: pest++ maps a
    control value into the tool's transform space, and out-of-bounds on a LOG parameter
    becomes NaN, which then gets run through the model. pyemu draws go out of bounds
    routinely, so this is a live path.
    """
    wd = _case("api_pyemu_bounds", noptmax=1)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        pst = ies.pst
        adj = pst.adj_par_names[0]
        ub = float(pst.parameter_data.loc[adj, "parubnd"])

        bad = ies.par_df(lower=True)
        bad.loc[:, adj] = ub * 1000.0

        try:
            ies.set_par_df(bad)
            raise AssertionError("an out-of-bounds parameter was accepted silently")
        except PestppError as e:
            assert adj in str(e), str(e)

        # reset clips instead, which is what pyemu's ParameterEnsemble.enforce() does
        ies.set_par_df(bad, enforce="reset")
        back = ies.par_df(lower=True)
        assert np.allclose(back.loc[:, adj].values, ub), \
            "enforce='reset' did not clip to the upper bound"

        # and opting out is still possible, for a caller who means it
        ies.set_par_df(back, enforce=False)
        ies.finalize()


def api_pyemu_from_pst_object_test():
    """from_pst accepts a pyemu.Pst, matching pyemu's own from_* convention."""
    wd = _case("api_pyemu_obj", noptmax=1)
    pst = pyemu.Pst(os.path.join(wd, "pest.pst"))
    pst.pestpp_options["ies_num_reals"] = 5
    # the object is renamed so the file it lands in is unmistakably the one this wrote, not
    # the pest.pst already sitting in the directory
    pst.filename = os.path.join(wd, "handed_over.pst")

    with Ies.from_pst(pst, workdir=wd) as ies:
        ies.initialize()
        assert ies.n_reals == 5, ies.n_reals
        assert os.path.exists(os.path.join(wd, "handed_over.pst")), \
            "the Pst object was not written to the working directory"
        ies.finalize()


def api_pyemu_results_test():
    """`.results` is a pyemu.Results over the working directory, with pyemu's vocabulary."""
    wd = _case("api_pyemu_res", noptmax=1)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        for _ in ies.iterations():
            pass
        ies.finalize()
        r = ies.results
        assert isinstance(r, pyemu.Results), type(r)
        paren0 = r.ies.paren0
        assert paren0 is not None and paren0.shape[0] == ies.n_reals, paren0.shape
        # the on-disk iteration-0 ensemble and the live one agree on their realizations
        assert set(str(i).lower() for i in paren0.index) == \
            set(str(i).lower() for i in ies.real_names)


def api_activate_zero_weighted_obs_test():
    """An observation that starts at zero weight can be switched ON, and actually contributes.

    This used to return cleanly and do nothing. The active observation set is fixed at
    initialize from the non-zero-weighted names, and everything downstream is sized from it:
    the weights ensemble has no column for a zero-weighted observation, and - the part that
    is easy to miss - neither does the NOISE ensemble, because no noise was ever drawn for it.
    Setting a weight and stopping there left the observation looking active while contributing
    nothing to phi. Staged history matching (heads first, then fluxes) needs this to work.

    Note what is deliberately NOT coupled: reweighting an already-active observation does not
    touch its noise. Weights and noise are independent in general, and redrawing noise on every
    weight change would make phi incomparable across iterations. Only the zero-to-nonzero
    transition is structural, because only then is there no noise to preserve.
    """
    for no_noise in (True, False):
        tag = "nonoise" if no_noise else "noise"
        wd = _case("api_activate_" + tag, noptmax=1, ies_no_noise=no_noise)
        pst = pyemu.Pst(os.path.join(wd, "pest.pst"))
        nz = pst.nnz_obs_names
        assert len(nz) > 1, "need at least two weighted observations to switch one off"
        victim = nz[-1]
        pst.observation_data.loc[victim, "weight"] = 0.0
        pst.write(os.path.join(wd, "pest.pst"), version=2)

        with Ies.from_pst("pest.pst", workdir=wd) as ies:
            ies.initialize()
            before_cols = list(ies.weights_df(lower=True).columns)
            assert victim not in before_cols, \
                "{0}: a zero-weighted observation should not be in the weights ensemble".format(tag)
            phi_before = ies.phi

            ies.set_obs_weights({victim: 5.0})

            after_cols = list(ies.weights_df(lower=True).columns)
            assert victim in after_cols, \
                "{0}: activating did not add the observation to the weights ensemble".format(tag)

            # the noise ensemble has to gain a column too, or phi is computed against nothing
            noise_cols = [c.lower() for c in ies._lib.get_ensemble_col_names(NOISE_EN)]
            assert victim in noise_cols, \
                "{0}: no noise realizations were generated for the activated observation - "\
                "it would contribute a residual against a value that does not exist".format(tag)

            # and it actually contributes
            ies.update_phi()
            contrib = ies.phi_obs_df(lower=True)
            assert victim in contrib.columns, \
                "{0}: the activated observation contributes no phi term".format(tag)
            assert ies.phi != phi_before, \
                "{0}: phi did not change after activating an observation".format(tag)
            ies.finalize()


def api_reweight_leaves_noise_alone_test():
    """Changing the weight of an ALREADY-active observation must not disturb its noise.

    The counterpart to the test above, and the reason activation is special-cased rather than
    every weight change triggering a redraw: noise realizations are what phi is measured
    against, so regenerating them mid-run would make this iteration's phi incomparable with
    the last one's. GMDSI workflows reweight repeatedly and rely on that.
    """
    wd = _case("api_reweight_noise", noptmax=1, ies_no_noise=False)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        target = ies.pst.nnz_obs_names[0]
        noise_before = ies.noise_df(lower=True) if hasattr(ies, "noise_df") else None
        if noise_before is None:
            import pandas as _pd
            arr, tok = ies._lib.get_ensemble_view(NOISE_EN)
            noise_before = _pd.DataFrame(
                arr.copy(),
                index=ies._lib.get_ensemble_row_names(NOISE_EN),
                columns=[c.lower() for c in ies._lib.get_ensemble_col_names(NOISE_EN)])
            ies._lib.release_view(tok)

        ies.set_obs_weights({target: 99.0})

        arr, tok = ies._lib.get_ensemble_view(NOISE_EN)
        after = arr.copy()
        ies._lib.release_view(tok)
        assert np.allclose(noise_before.values, after), \
            "reweighting an active observation redrew its noise; phi is no longer comparable "\
            "with the previous iteration"
        ies.finalize()


def api_reinflate_size_is_capped_by_prior_test():
    """Asking to reinflate to more realizations than the prior holds is refused, not ignored.

    Reinflation SELECTS realizations out of pe_base rather than generating new ones, so the
    prior's row count is a hard ceiling. The underlying routine deals with an over-large
    request by silently reinflating to the prior's size instead - which from the outside looks
    exactly like the call having done nothing, with no way to find out why. The error message
    is the feature here: it has to name the ceiling and the way around it.
    """
    wd = _case("api_reinflate_range", noptmax=1, num_reals=6)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        # both signs, because the sign selects the anomaly source, not the size
        for n in (20, -20):
            try:
                ies.reinflate(factor=1.0, num_reals=n)
            except PestppError as e:
                assert "prior ensemble holds 6" in str(e), str(e)
                assert "ies_reinflate_num_reals" in str(e), \
                    "the error should say how to grow the ensemble: {0}".format(e)
            else:
                raise AssertionError("num_reals={0} exceeds the prior and should raise".format(n))

        # within the ceiling it still works, and the size really does change
        ies.reinflate(factor=1.0, num_reals=4)
        assert ies.n_reals == 4, ies.n_reals
        ies.finalize()


def api_reinflate_grows_from_bigger_prior_test():
    """The documented route to a LARGER ensemble mid-run: big prior, start on a subset.

    This is what the executable does with ies_reinflate_num_reals, and the only way to grow
    without inventing realizations from nowhere: draw the prior at the largest size the run
    will ever need, truncate the working ensemble during initialize, then reinflate back up
    into the spare rows. Worth a test because it is the answer to "why did my ensemble not
    get bigger", and it depends on two options cooperating.
    """
    wd = _case("api_reinflate_grow", noptmax=1, num_reals=20, ies_reinflate_num_reals="8")
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        assert ies.n_reals == 8, \
            "ies_reinflate_num_reals[0] should have truncated the working ensemble: {0}".format(
                ies.n_reals)
        ies.reinflate(factor=1.0, num_reals=16)
        assert ies.n_reals == 16, \
            "reinflation did not grow the ensemble into the prior's spare rows: {0}".format(
                ies.n_reals)
        ies.finalize()


def api_reinflate_centring_is_per_call_test():
    """center_on_min_phi is reachable per call, instead of only through an option's SIGN.

    In the shipped loop, centring on the min-phi realization is derived from a negative entry
    in ies_n_iter_reinflate - a schedule-wide flag on an unrelated setting, which a caller
    driving one reinflation at a time cannot express. The argument is checked here against the
    .rec file, because the choice of offset leaves no other visible trace.
    """
    wd = _case("api_reinflate_center", noptmax=1, num_reals=6)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        ies.solve()          # so the phi map has something to pick a minimum from
        ies.reinflate(factor=1.0, num_reals=6, center_on_min_phi=True)
        rec = open(os.path.join(wd, "pest.rec")).read()
        assert "using min-phi realization for offset" in rec, \
            "center_on_min_phi=True did not reach the min-phi branch"
        ies.finalize()

    wd = _case("api_reinflate_center_bad", noptmax=1, num_reals=6)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        try:
            ies._lib.reinflate_ensemble(1.0, 0, 7)      # only -1, 0, 1 are meaningful
        except PestppError as e:
            assert "center_on_min_phi must be" in str(e), str(e)
        else:
            raise AssertionError("center_on_min_phi=7 should have raised")
        ies.finalize()


def api_reinflate_keeps_caller_weights_test():
    """Reinflation must not revert a weight the caller set.

    Reinflation restores the scenario's weights from org_obs_info when it is done, and that
    was captured once at initialize. So the staged workflow - activate observations, then
    reinflate so the ensemble has the spread to respond to them - silently put the newly
    activated observations back to zero weight, undoing the switch it was called to support.

    Note what to assert on: the weights ENSEMBLE survived this even before the fix, because
    it is restored from weights_base, which activation does update. Only the scenario weights
    reverted, so that is what this checks.
    """
    wd = _case("api_reinflate_weights", noptmax=1, num_reals=6)
    pst = pyemu.Pst(os.path.join(wd, "pest.pst"))
    victim = pst.nnz_obs_names[0]
    pst.observation_data.loc[victim, "weight"] = 0.0
    pst.write(os.path.join(wd, "pest.pst"), version=2)

    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        ies.set_obs_weights({victim: 5.0})

        def scenario_weight():
            w = ies.obs_weights
            w.index = [n.lower() for n in w.index]
            return w[victim]

        assert scenario_weight() == 5.0, "activation did not take"
        ies.reinflate(factor=1.0, num_reals=6)
        assert scenario_weight() == 5.0, \
            "reinflation reverted the caller-set weight to {0}".format(scenario_weight())
        assert ies.weights_df(lower=True)[victim].mean() == 5.0, \
            "the weights ensemble lost the caller-set weight"
        ies.finalize()


def api_reinflate_grows_after_activation_test():
    """The staged sequence end to end: solve, activate observations, THEN grow the ensemble.

    Each piece worked alone; together they did not. Reinflation restores the weight ensemble
    from weights_base, and three things have to line up for that to be sound:

      * initialize() builds the weight ensemble at the PRIOR's size, before the working
        ensemble is truncated to ies_reinflate_num_reals[0], so it starts with spare rows
      * the first solve resizes the LIVE weight ensemble down to the working ensemble
      * activation splices new columns into both

    Activation used to do that last step by assigning the trimmed live ensemble over
    weights_base, which discarded the spare rows. Nothing failed at that point - it failed at
    the next reinflation that grew the ensemble, as Ensemble::get_eigen() reporting missing
    REALIZATION names, which points nowhere near weights.
    """
    wd = _case("api_reinflate_after_activate", noptmax=3, num_reals=20,
               ies_reinflate_num_reals="10")
    pst = pyemu.Pst(os.path.join(wd, "pest.pst"))
    victim = pst.nnz_obs_names[0]
    pst.observation_data.loc[victim, "weight"] = 0.0
    pst.write(os.path.join(wd, "pest.pst"), version=2)

    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        assert ies.n_reals == 10, ies.n_reals
        ies.solve()                              # trims the live weight ensemble to 10
        ies.set_obs_weights({victim: 5.0})       # activation splices a column into both
        ies.reinflate(factor=0.5, num_reals=20)  # grows into weights_base's spare rows

        assert ies.n_reals == 20, \
            "reinflation after activation did not grow the ensemble: {0}".format(ies.n_reals)
        w = ies.weights_df(lower=True)
        assert w.shape[0] == 20, \
            "the weight ensemble did not follow the parameter ensemble: {0}".format(w.shape)
        assert victim in w.columns, "the activated observation fell out of the weight ensemble"
        assert (w[victim] == 5.0).all(), \
            "the activated weight did not reach every realization: {0}".format(w[victim].unique())
        # and the run is still driveable afterwards
        ies.solve()
        ies.finalize()


def api_option_values_are_native_types_test():
    """Options take python values, not control-file spellings - lists especially.

    A bare str() is right for scalars by luck and wrong for sequences by construction: the
    vector options are tokenized on commas, so a list arrived as the literal "[0.1, 1.0]" and
    was rejected under the OPTION's name, which reads like the option was wrong rather than
    the brackets. Sequences are the whole point of this - the lambda multipliers and the
    reinflation schedule are the options a caller most wants to change mid-run.
    """
    wd = _case("api_option_types", noptmax=1, num_reals=6)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        for key, value, expect in (
                ("ies_lambda_mults", [0.1, 1.0, 10.0], "0.1,1,10"),
                ("ies_reinflate_factor", (1.0, 0.5), "1,0.5"),
                ("ies_n_iter_reinflate", np.array([3, 999]), "3,999"),
                ("ies_n_iter_reinflate", range(2, 5), "2,3,4"),
                ("ies_num_reals", 20, "20"),
                ("ies_num_reals", np.int64(30), "30"),
                ("ies_subset_size", -10, "-10"),
                ("ies_par_en", "prior.jcb", "prior.jcb"),
        ):
            ies.set_option(key, value)
            assert ies.get_option(key) == expect, \
                "{0}={1!r} read back as {2!r}, expected {3!r}".format(
                    key, value, ies.get_option(key), expect)

        # bools go in as bools rather than "True", which the parser only accepts by accident
        ies.set_option("ies_no_noise", True)
        assert ies.get_option("ies_no_noise") == "1", ies.get_option("ies_no_noise")
        ies.set_option("ies_no_noise", False)
        assert ies.get_option("ies_no_noise") == "0", ies.get_option("ies_no_noise")

        # and a value with no sensible spelling is refused HERE, naming the type, rather than
        # being stringified into something the parser rejects later under the option's name
        for bad, exc in (({0.1, 1.0}, TypeError), ([], ValueError), (object(), TypeError)):
            try:
                ies.set_option("ies_lambda_mults", bad)
            except exc:
                pass
            else:
                raise AssertionError("{0!r} should have raised {1}".format(bad, exc.__name__))
        ies.finalize()


def api_option_sequence_reaches_the_algorithm_test():
    """A list option is not just parsed, it takes effect - checked through behaviour.

    get_option() round-tripping proves the string was formed correctly and nothing more. This
    drives the one sequence option with a visible consequence: the first entry of
    ies_reinflate_num_reals truncates the working ensemble during initialize.
    """
    wd = _case("api_option_seq_effect", noptmax=1, num_reals=20,
               ies_reinflate_num_reals=[8, 16])
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        assert ies.n_reals == 8, \
            "the list-valued option did not reach the algorithm: {0} reals".format(ies.n_reals)
        ies.finalize()


def _mou_case(name, noptmax=1, pop=10):
    """A small mou case in its own directory, from the g07 template."""
    base = os.path.join(_BENCH, "g07", "template")
    d = os.path.join(_BENCH, name)
    if os.path.exists(d):
        shutil.rmtree(d)
    shutil.copytree(base, d)
    pst = pyemu.Pst(os.path.join(d, "g07.pst"))
    for k in [k for k in pst.pestpp_options if k.startswith("sqp_")]:
        pst.pestpp_options.pop(k)
    pst.pestpp_options["mou_population_size"] = pop
    pst.pestpp_options["mou_generator"] = "de"
    pst.pestpp_options["random_seed"] = 11
    pst.control_data.noptmax = noptmax
    pst.write(os.path.join(d, "g07.pst"), version=2)
    return d


def api_deferred_solve_drives_the_runs_test():
    """solve -> inspect candidates -> run -> finish, with the caller owning every batch.

    The whole point is that no model run happens inside a library call. ies generates one
    candidate per lambda x scale-factor combination, runs a SUBSET of the ensemble against all
    of them to pick a winner, then runs the REST of the ensemble at that winner - two batches,
    and finish_solve(defer_runs=True) hands back the second rather than running it.
    """
    wd = _case("api_defer_ies", noptmax=2, num_reals=8, ies_subset_size=4)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        n = ies.solve(defer_runs=True)
        cands = ies.candidates()
        assert len(cands) > 1, "expected a candidate per lambda x scale combination"
        assert n == len(cands) * 4, \
            "{0} runs implied by {1} candidates over a subset of 4".format(n, len(cands))
        # the candidates are ordinary ensembles: same frames, same views
        assert cands[0].par_df().shape[0] == ies.n_reals
        assert list(cands[0].par_df().index) == list(ies.par_df().index)
        # every candidate carries the factors it was generated with
        assert all(c.inflation > 0 for c in cands), [c.inflation for c in cands]

        assert ies.queue_runs() == n
        ies.run()
        assert ies.process_runs() == 0

        step = ies.finish_solve(defer_runs=True)
        assert step.pending_runs == 4, \
            "the 4 realizations outside the subset should be handed back: {0}".format(
                step.pending_runs)
        ies.queue_runs()
        ies.run()
        ies.process_runs()
        step = ies.finish_solve(defer_runs=True)
        assert step.pending_runs == 0, step.pending_runs
        assert np.isfinite(step.phi_mean)

        # and the composed path still works on the same session afterwards
        after = ies.solve()
        assert after.pending_runs == 0
        assert after.iter == step.iter + 1, (after.iter, step.iter)
        ies.finalize()


def api_deferred_solve_matches_composed_test():
    """The deferred path and solve() produce the same iteration, to the last digit.

    This is the test that matters: the decomposition is only safe if driving the stages by
    hand is the same computation as letting solve() do it. Same seed, same case, one iteration
    each way, compared on the resulting ensembles rather than on phi alone.
    """
    def composed():
        wd = _case("api_defer_cmp_a", noptmax=1, num_reals=8, ies_subset_size=4)
        with Ies.from_pst("pest.pst", workdir=wd) as ies:
            ies.initialize()
            ies.solve()
            out = (ies.par_df().copy(), ies.obs_df().copy(), ies.phi)
            ies.finalize()
        return out

    def deferred():
        wd = _case("api_defer_cmp_b", noptmax=1, num_reals=8, ies_subset_size=4)
        with Ies.from_pst("pest.pst", workdir=wd) as ies:
            ies.initialize()
            ies.solve(defer_runs=True)
            ies.queue_runs(); ies.run(); ies.process_runs()
            step = ies.finish_solve(defer_runs=True)
            while step.pending_runs:
                ies.queue_runs(); ies.run(); ies.process_runs()
                step = ies.finish_solve(defer_runs=True)
            out = (ies.par_df().copy(), ies.obs_df().copy(), ies.phi)
            ies.finalize()
        return out

    par_a, obs_a, phi_a = composed()
    par_b, obs_b, phi_b = deferred()
    assert par_a.shape == par_b.shape, (par_a.shape, par_b.shape)
    assert np.allclose(par_a.values, par_b.values), \
        "deferred solve produced a different parameter ensemble than solve()"
    assert np.allclose(obs_a.values, obs_b.values), \
        "deferred solve produced a different observation ensemble than solve()"
    assert phi_a == phi_b, (phi_a, phi_b)


def api_deferred_solve_edit_reaches_the_runs_test():
    """A write through a candidate view is what gets run - otherwise inspection is all it is."""
    wd = _case("api_defer_edit", noptmax=1, num_reals=6, ies_subset_size=6)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        ies.solve(defer_runs=True)
        cand = ies.candidates()[0]
        before = cand.par_df()
        with cand.par_view() as v:
            v[:, 0] = v[:, 0] * 0.5
        after = cand.par_df()
        assert not np.allclose(before.values, after.values), \
            "the write did not land in the candidate ensemble"
        # and it survives to the queue: the run manager takes what is in the ensemble now
        ies.queue_runs()
        ies.run()
        ies.process_runs()
        step = ies.finish_solve()
        assert step.pending_runs == 0
        ies.finalize()


def api_deferred_solve_subset_override_test():
    """Naming realizations replaces the algorithm's subset, and the rest become the remainder."""
    wd = _case("api_defer_subset", noptmax=1, num_reals=8, ies_subset_size=4)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        n_default = ies.solve(defer_runs=True)
        cands = len(ies.candidates())
        reals = list(ies.par_df().index)[:2]
        n = ies.queue_runs(reals=reals)
        assert n == cands * 2, \
            "naming 2 realizations should queue 2 per candidate, got {0}".format(n)
        ies.run(); ies.process_runs()
        step = ies.finish_solve(defer_runs=True)
        # 8 realizations, 2 named -> 6 left over
        assert step.pending_runs == 6, \
            "the unnamed realizations should become the remainder: {0}".format(
                step.pending_runs)
        assert n_default > n, "the default subset was larger, so it should have implied more runs"
        ies.queue_runs(); ies.run(); ies.process_runs()
        assert ies.finish_solve(defer_runs=True).pending_runs == 0
        ies.finalize()

        # a name that is not in the ensemble is refused rather than silently skipped
        ies2 = Ies.from_pst("pest.pst", workdir=wd)
        ies2.initialize()
        ies2.solve(defer_runs=True)
        try:
            ies2.queue_runs(reals=["not_a_realization"])
        except PestppError as e:
            assert "no such realization" in str(e), str(e)
        else:
            raise AssertionError("an unknown realization name should raise")
        ies2.close()


def api_deferred_solve_mou_test():
    """mou has the same shape with one candidate population and never a second batch."""
    # noptmax=2 because the test runs a second generation below: mou past its noptmax throws
    # from deep in the archive bookkeeping, on the composed path just the same
    d = _mou_case("api_defer_mou", noptmax=2, pop=10)
    with Mou.from_pst("g07.pst", workdir=d) as mou:
        mou.initialize()
        n = mou.solve(defer_runs=True)
        cands = mou.candidates()
        assert len(cands) == 1, "mou generates one candidate population: {0}".format(len(cands))
        assert n == cands[0].par_df().shape[0], (n, cands[0].par_df().shape)
        mou.queue_runs()
        mou.run()
        mou.process_runs()
        step = mou.finish_solve(defer_runs=True)
        assert step.pending_runs == 0, "mou never needs a second batch"
        # the generation counter has to advance, or every generation overwrites the same files
        assert step.iter == 1, step.iter
        assert mou.solve().iter == 2, "the composed path should carry on from the deferred one"
        mou.finalize()


def api_deferred_solve_refused_where_it_does_not_fit_test():
    """da and sqp say why they cannot defer, rather than approximating the shape.

    Both are structural: da's one advance() is a whole noptmax loop over a cycle, and sqp's
    line search proposes, runs and judges a step and may try again, so one iteration issues
    several run batches. Neither is one generate -> run -> evaluate.
    """
    d = _mou_case("api_defer_sqp", noptmax=1)
    pst = pyemu.Pst(os.path.join(d, "g07.pst"))
    for k in [k for k in pst.pestpp_options if k.startswith("mou_")]:
        pst.pestpp_options.pop(k)
    pst.pestpp_options["sqp_num_reals"] = 8
    pst.write(os.path.join(d, "g07.pst"), version=2)
    with Sqp.from_pst("g07.pst", workdir=d) as sqp:
        sqp.initialize()
        try:
            sqp.solve(defer_runs=True)
        except PestppError as e:
            assert "sqp cannot defer" in str(e), str(e)
            assert "several run batches" in str(e), str(e)
        else:
            raise AssertionError("sqp should refuse a deferred solve")
        sqp.finalize()


def api_deferred_solve_state_guards_test():
    """The sequence is enforced, because every way of getting it wrong is silent otherwise."""
    wd = _case("api_defer_guards", noptmax=2, num_reals=6, ies_subset_size=3)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()

        # finishing before starting
        try:
            ies.finish_solve()
        except PestppError as e:
            assert "no deferred solve is open" in str(e), str(e)
        else:
            raise AssertionError("finish_solve before solve_prepare should raise")

        ies.solve(defer_runs=True)
        # a second solve on top of an open one
        try:
            ies.solve(defer_runs=True)
        except PestppError as e:
            assert "already open" in str(e), str(e)
        else:
            raise AssertionError("a second deferred solve should raise")
        # the composed path must not run underneath an open deferred solve either
        try:
            ies.solve()
        except PestppError as e:
            assert "deferred solve is open" in str(e), str(e)
        else:
            raise AssertionError("solve() during a deferred solve should raise")
        # finishing before the runs have been processed
        ies.queue_runs()
        try:
            ies.finish_solve()
        except PestppError as e:
            assert "not been processed" in str(e), str(e)
        else:
            raise AssertionError("finish_solve with runs in flight should raise")
        ies.run(); ies.process_runs()
        step = ies.finish_solve()
        assert step.pending_runs == 0, "defer_runs=False should complete in one call"
        ies.finalize()


def api_mou_runs_past_noptmax_test():
    """A caller driving its own loop can run more generations than noptmax.

    The population schedule is built once at initialize() for generations 0..noptmax, and
    both readers of it were wrong past the end: generate_population() used at() and threw a
    bare "map::at: key not found" - nothing in the message connects that to a population
    schedule - while update_archive_spea() used operator[], which inserted a size of ZERO and
    carried on with an empty archive. noptmax is a bound on the SHIPPED loop, not on what the
    API can do, so the schedule holds its last value instead.
    """
    d = _mou_case("api_mou_past_noptmax", noptmax=1, pop=10)
    with Mou.from_pst("g07.pst", workdir=d) as mou:
        mou.initialize()
        for expect in (1, 2, 3):
            step = mou.solve()
            assert step.iter == expect, (step.iter, expect)
            assert mou.n_reals > 0, \
                "generation {0} produced an empty population".format(step.iter)
        mou.finalize()


def api_deferred_solve_needs_upgrades_in_memory_test():
    """Deferring is refused when the candidates were spilled to disk, and says why.

    With 'ies_upgrades_in_memory' false the candidate ensembles are files, so there is nothing
    in memory to inspect or edit - the whole reason to defer. Refusing beats handing back
    views onto ensembles that were emptied on the way to disk.
    """
    wd = _case("api_defer_spill", noptmax=1, num_reals=6, ies_upgrades_in_memory=False)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        try:
            ies.solve(defer_runs=True)
        except PestppError as e:
            assert "ies_upgrades_in_memory" in str(e), str(e)
        else:
            raise AssertionError("a deferred solve should be refused on the spill path")
        # and the composed path still works, on the same session, right after the refusal
        step = ies.solve()
        assert np.isfinite(step.phi_mean), "the refusal left the session unusable"
        ies.finalize()


def api_progress_renders_in_a_terminal_test():
    """A tty gets one rewriting line; a redirected stream gets whole readable lines.

    The distinction is the point. Carriage returns in a log file collapse the whole run into
    one unreadable line, and `quiet=True` sends output to a log by default - so the same
    renderer has to behave differently depending on where it is pointed.
    """
    import io
    sys.path.insert(0, os.path.join(_REPO, "python"))
    import pestpp_progress as pp

    class _Tty(io.StringIO):
        def isatty(self):
            return True

    tty = _Tty()
    bar = pp.TextProgress(stream=tty)
    bar.min_interval = 0
    bar.start("running", total=8)
    for i in range(1, 9):
        bar.update(done=i, failed=1 if i > 5 else 0)
    bar.close()
    out = tty.getvalue()
    assert "\r" in out, "a tty should be redrawn in place"
    assert out.endswith("\n"), "the final frame should end the line rather than leaving it open"
    final = out.split("\r")[-1]
    assert "8/8" in final and "100%" in final, final
    assert "failed=1" in final, final

    plain = io.StringIO()            # isatty() is False
    bar = pp.TextProgress(stream=plain)
    bar.min_interval = 0
    bar.start("running", total=8)
    for i in range(1, 9):
        bar.update(done=i)
    bar.note("lambda rejected, retrying")
    bar.close()
    text = plain.getvalue()
    assert "\r" not in text, "a redirected stream must not be given carriage returns"
    assert "lambda rejected, retrying" in text
    lines = text.splitlines()
    assert 3 <= len(lines) <= 12, "expected sparse whole lines, got {0}".format(len(lines))


def api_progress_is_inert_when_disabled_test():
    """progress=False costs nothing and every hook is still safe to call."""
    sys.path.insert(0, os.path.join(_REPO, "python"))
    import pestpp_progress as pp

    bar = pp.auto(False)
    assert type(bar) is pp.Progress, type(bar)
    # the no-op has to accept the same calls in any order, or callers end up branching
    bar.update(done=1)
    bar.start("x", total=2)
    bar.note("y")
    bar.close()
    bar.close()

    os.environ["PESTPP_NO_PROGRESS"] = "1"
    try:
        assert type(pp.auto(True)) is pp.Progress, "PESTPP_NO_PROGRESS should switch it off"
    finally:
        os.environ.pop("PESTPP_NO_PROGRESS")


def api_progress_renders_in_a_notebook_test():
    """In a real kernel it picks the notebook renderer and updates ONE output in place.

    Executed in a kernel rather than mocked, because the thing being tested is a property of
    the display protocol: many updates must collapse to a single output area. A mock would
    assert that update_display was called and prove nothing about what the notebook shows.
    """
    import nbformat
    from nbclient import NotebookClient

    src = (
        "import sys\n"
        "sys.path.insert(0, {0!r})\n"
        "import pestpp_progress as pp\n"
        "print('in_notebook:', pp.in_notebook())\n"
        "bar = pp.auto(True)\n"
        "print('renderer:', type(bar).__name__)\n"
        "bar.min_interval = 0\n"
        "bar.start('running', total=6)\n"
        "for i in range(1, 7):\n"
        "    bar.update(done=i, phi=100.0 / i)\n"
        "bar.close()\n"
    ).format(os.path.join(_REPO, "python"))

    nb = nbformat.v4.new_notebook(cells=[nbformat.v4.new_code_cell(src)])
    NotebookClient(nb, timeout=300, kernel_name="python3", allow_errors=True).execute()

    outs = nb.cells[0].get("outputs", [])
    errors = [o for o in outs if o.output_type == "error"]
    assert not errors, "{0}: {1}".format(errors[0].ename, errors[0].evalue)
    stream = "".join(o.get("text", "") for o in outs if o.output_type == "stream")
    assert "in_notebook: True" in stream, stream
    assert "renderer: NotebookProgress" in stream, stream

    displays = [o for o in outs if o.output_type in ("display_data", "update_display_data")]
    assert len(displays) == 1, \
        "six updates should collapse into one output area, got {0}".format(len(displays))
    html = displays[0].get("data", {}).get("text/html", "")
    assert "6/6" in html and "100%" in html, html[:400]


def api_progress_during_a_real_run_test():
    """progress=True on a live run reports the run manager's counters and gets out of the way.

    Also the check that it cannot break a run: the counters come from the run manager, and a
    manager that does not keep them must cost a bar, not the batch.
    """
    import io
    sys.path.insert(0, os.path.join(_REPO, "python"))
    import pestpp_progress as pp

    wd = _case("api_progress_run", noptmax=1, num_reals=6)
    sink = io.StringIO()
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        n = ies.queue_runs()
        ies.run(progress=pp.TextProgress(stream=sink))
        failed = ies.process_runs()
        assert failed == 0, failed
        text = sink.getvalue()
        assert "running" in text, text
        assert str(n) in text, "the batch size should appear in the bar: {0}".format(text)

        # and across iterations, annotated with phi
        sink2 = io.StringIO()
        steps = list(ies.iterations(max_iter=1, progress=pp.TextProgress(stream=sink2)))
        assert len(steps) == 1
        assert "phi=" in sink2.getvalue(), sink2.getvalue()
        ies.finalize()


def api_run_observer_fires_during_a_composed_solve_test():
    """The observer sees runs inside solve() - the call that is otherwise completely silent.

    Every other way of watching runs requires the caller to own them (queue/run/process, or a
    deferred solve). This is the one that works for someone driving `iterations()`, which is
    most people.
    """
    wd = _case("api_observer_solve", noptmax=1, num_reals=8)
    seen = []
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        ies._lib.set_run_observer(lambda p: seen.append(dict(p)) or True, 0.0)
        ies.solve()
        ies._lib.set_run_observer(None)
        assert seen, "the observer was never called during solve()"
        assert max(p["n_completed"] for p in seen) > 0, "no completions were reported"
        assert max(p["n_total"] for p in seen) > 0, "the batch size was never reported"
        # the struct the library filled in must be the one we declared
        assert seen[0]["run_id"] is not None
        ies.finalize()


def api_run_observer_is_throttled_in_the_library_test():
    """Throttling is at the source, because an observer cannot decline a call it was handed.

    The run managers poll in a hot loop, so an unthrottled observer pays a cross-ABI call per
    poll rather than per run - hundreds of thousands of them on a case this size.
    """
    wd = _case("api_observer_throttle", noptmax=1, num_reals=8)
    counts = {}
    for interval in (0.0, 0.25):
        n = [0]
        with Ies.from_pst("pest.pst", workdir=wd) as ies:
            ies.initialize()
            ies._lib.set_run_observer(lambda p: n.__setitem__(0, n[0] + 1) or True, interval)
            ies.solve()
            ies._lib.set_run_observer(None)
            ies.finalize()
        counts[interval] = n[0]
    assert counts[0.25] < counts[0.0], counts
    assert counts[0.25] < 500, \
        "a throttled observer should be called a handful of times, got {0}".format(counts[0.25])


def api_run_observer_reentrancy_is_an_allowlist_test():
    """Mid-batch, run management is legal and everything else is refused BY NAME.

    Not a ban: preemption needs exactly this door - look at what is running, decide, cancel -
    so the rule has to admit the run-management calls while keeping out the ones that would
    read a part-updated ensemble. Enforced rather than documented, because a rule nobody can
    check is not a rule.
    """
    wd = _case("api_observer_reentry", noptmax=1, num_reals=6)
    result = {}
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()

        def probe(p):
            try:
                ies._lib.get_phi_summary(0)
                result["phi"] = "ALLOWED"
            except PestppError as e:
                result["phi"] = str(e)
            try:
                ies._lib.cancel_runs([])
                result["cancel"] = "ALLOWED"
            except PestppError as e:
                result["cancel"] = str(e)
            return True

        ies._lib.set_run_observer(probe, 1e9)
        ies.solve()
        ies._lib.set_run_observer(None)
        ies.finalize()

    assert "cannot be called from inside a run observer" in result["phi"], result["phi"]
    assert "pestpp_get_phi_summary" in result["phi"], "the message should name the call"
    # cancel_runs is on the allowlist, so it reaches its OWN argument check rather than the
    # observer guard - which is what proves it got through
    assert "cannot be called from inside a run observer" not in result["cancel"], result["cancel"]
    assert "run id" in result["cancel"], result["cancel"]


def api_run_observer_can_stop_a_batch_test():
    """Returning False stops the batch, keeping the runs that already finished.

    The observer returns an action rather than nothing precisely so this - and later,
    preemption - is a new return VALUE rather than a new callback.
    """
    wd = _case("api_observer_stop", noptmax=1, num_reals=12)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize(defer_runs=True)
        queued = ies.queue_runs()
        assert queued == 12, queued
        seen = []

        def stop_after_3(p):
            seen.append(p["n_completed"])
            return p["n_completed"] < 3

        ies._lib.set_run_observer(stop_after_3, 0.0)
        ies.run()
        ies._lib.set_run_observer(None)
        assert max(seen) == 3, \
            "the batch should have stopped at 3 of {0}, saw {1}".format(queued, max(seen))
        ies.close()


def api_run_observer_survives_a_raising_callback_test():
    """An observer that raises is dropped, not allowed to unwind through runs in flight.

    A python exception crossing back into C++ through a function pointer is undefined
    behaviour on some ABIs, and a progress bar has no business failing a batch either way.
    """
    wd = _case("api_observer_raises", noptmax=1, num_reals=6)
    calls = [0]
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()

        def bad(p):
            calls[0] += 1
            raise ValueError("observer blew up")

        ies._lib.set_run_observer(bad, 0.0)
        step = ies.solve()          # must complete regardless
        ies._lib.set_run_observer(None)
        assert calls[0] > 0, "the observer was never called"
        assert np.isfinite(step.phi_mean), "the batch did not survive a raising observer"
        ies.finalize()


def api_run_observer_thunk_is_retained_test():
    """The ctypes thunk is held by the session, or C ends up pointing at collected memory.

    The classic ctypes footgun: a CFUNCTYPE kept only by a local is garbage-collected while
    the library still holds the pointer, and the crash lands at some unrelated later run.
    """
    wd = _case("api_observer_thunk", noptmax=1, num_reals=6)
    with Ies.from_pst("pest.pst", workdir=wd) as ies:
        ies.initialize()
        assert ies._lib._observer_thunk is None
        ies._lib.set_run_observer(lambda p: True, 0.1)
        assert ies._lib._observer_thunk is not None, \
            "the thunk must be referenced by the session for as long as C holds it"
        ies._lib.set_run_observer(None)
        assert ies._lib._observer_thunk is None, "unregistering should release it"
        ies.finalize()


def api_session_cleanup_does_not_need_with_test():
    """A dropped session cleans itself up: `with` is a convenience, not the contract.

    This is the notebook case, and `with` is no help for it. Re-running a setup cell rebinds
    the name and drops the previous session - so if cleanup only happened in __exit__ or in an
    explicit close(), every re-run would silently leak the handle AND leave that session's
    panther agents running. Dropping the reference is what has to be safe.
    """
    wd = _case("api_no_with", noptmax=1, num_reals=6)

    # 1. flat style, no context manager anywhere, runs to completion
    ies = Ies.from_pst("pest.pst", workdir=wd)
    ies.initialize()
    assert ies.phi_df() is not None
    fin = ies._finalizer
    assert fin.alive, "a live session should have a live finalizer"

    # 2. close() runs it, and is idempotent - the repeat calls a flat style invites are free
    ies.close()
    assert not fin.alive, "close() should have run the finalizer"
    ies.close()

    # 3. and dropping the last reference is enough on its own
    ies2 = Ies.from_pst("pest.pst", workdir=wd)
    ies2.initialize()
    fin2 = ies2._finalizer
    lib_fin = ies2._lib._finalizer
    del ies2
    gc.collect()
    assert not fin2.alive, "dropping the session should have finalized it"
    assert not lib_fin.alive, "...and released the underlying handle with it"


def api_deferred_initialize_sqp_test():
    """sqp's initialize() splits too, and the split must equal the atomic one exactly.

    sqp was the awkward one: its initial ensemble is not evaluated at the top of initialize()
    but deeper, inside prep_4_ensemble_grad(), which also does the gradient setup that depends
    on the results. So that function had to split as well, and initialize() composes around it.

    Same claim as mou: splitting changes WHO runs the batch and nothing else. Note the count is
    legitimately 0 on several paths - finite-difference gradients (whose runs are perturbations,
    not candidates), the diagnostics, and a supplied rather than drawn ensemble - so this drives
    whichever path the case takes and asserts equivalence either way.
    """
    import filecmp

    def _case(name):
        d = _mou_case(name, noptmax=1)
        pst = pyemu.Pst(os.path.join(d, "g07.pst"))
        for k in [k for k in pst.pestpp_options if k.startswith("mou_")]:
            pst.pestpp_options.pop(k)
        pst.pestpp_options["sqp_num_reals"] = 8
        pst.pestpp_options["random_seed"] = 11
        pst.write(os.path.join(d, "g07.pst"), version=2)
        return d

    d_atomic = _case("api_init_sqp_atomic")
    with Sqp.from_pst("g07.pst", workdir=d_atomic) as sqp:
        sqp.initialize()
        atomic_dv = sqp.par_df().copy()
        atomic_obs = sqp.obs_df().copy()

    d_split = _case("api_init_sqp_split")
    with Sqp.from_pst("g07.pst", workdir=d_split) as sqp:
        n = sqp.initialize(defer_runs=True)
        assert n > 0, \
            "sqp should hand back its initial ensemble runs on the ensemble-gradient " \
            "path, got {0}".format(n)
        drawn = sqp.par_df()
        assert drawn.shape[0] >= n, (drawn.shape, n)
        sqp.queue_runs()
        sqp.run()
        sqp.process_runs()
        sqp.finish_initialize()
        split_dv = sqp.par_df().copy()
        split_obs = sqp.obs_df().copy()

    assert list(split_dv.index) == list(atomic_dv.index), \
        "the split path produced a different ensemble:\n atomic: {0}\n split : {1}".format(
            list(atomic_dv.index), list(split_dv.index))
    assert np.allclose(split_dv.values, atomic_dv.values, rtol=0, atol=1e-10), \
        "the split path produced different decision variable values"
    assert np.allclose(split_obs.values, atomic_obs.values, rtol=0, atol=1e-10), \
        "the split path produced different observation values"

    for f in sorted(os.listdir(d_atomic)):
        if not f.endswith(".csv") or not f.startswith("g07."):
            continue
        b = os.path.join(d_split, f)
        if not os.path.exists(b):
            continue
        assert filecmp.cmp(os.path.join(d_atomic, f), b, shallow=False), \
            "{0} differs between the atomic and split initialize".format(f)


def api_deferred_initialize_mou_test():
    """mou's initialize() splits, and the split must equal the atomic one exactly.

    This is the feature that was impossible for mou until now: initialize() drew the initial
    population AND evaluated it in one uninterruptible call, so the members that got run were
    always the ones mou drew. defer_runs=True stops in between, and that window is the only
    point at which the population can be replaced.

    The claim is that splitting changes WHO runs the batch and nothing else, so this drives
    both paths on the same case and seed and compares the written populations. A restart case
    hands back 0 runs and must still finish - that path is asserted too, because "returns 0"
    and "does nothing" are easy to conflate and only one of them is correct.
    """
    import filecmp

    # 1. the composed path, for reference
    d_atomic = _mou_case("api_init_mou_atomic", noptmax=1, pop=10)
    with Mou.from_pst("g07.pst", workdir=d_atomic) as mou:
        mou.initialize()
        atomic_dv = mou.par_df().copy()
        atomic_obs = mou.obs_df().copy()

    # 2. the split path: prepare hands back the runs, we service them, then finish
    d_split = _mou_case("api_init_mou_split", noptmax=1, pop=10)
    with Mou.from_pst("g07.pst", workdir=d_split) as mou:
        n = mou.initialize(defer_runs=True)
        assert n > 0, \
            "mou's initialize should now hand back the initial population's runs, got {0}".format(n)
        # the window that makes this worth having: the population is in memory and unevaluated
        drawn = mou.par_df()
        assert drawn.shape[0] == n, (drawn.shape, n)
        mou.queue_runs()
        mou.run()
        mou.process_runs()
        mou.finish_initialize()
        split_dv = mou.par_df().copy()
        split_obs = mou.obs_df().copy()

    # 3. same members, same values - splitting changes who runs the batch, nothing else
    assert list(split_dv.index) == list(atomic_dv.index), \
        "the split path produced a different population:\n atomic: {0}\n split : {1}".format(
            list(atomic_dv.index), list(split_dv.index))
    assert np.allclose(split_dv.values, atomic_dv.values, rtol=0, atol=1e-10), \
        "the split path produced different decision variable values"
    assert np.allclose(split_obs.values, atomic_obs.values, rtol=0, atol=1e-10), \
        "the split path produced different observation values"

    # 4. and the files agree too - the reporting in finish() is half the work it does
    for f in sorted(os.listdir(d_atomic)):
        if not f.endswith(".csv") or not f.startswith("g07."):
            continue
        b = os.path.join(d_split, f)
        if not os.path.exists(b):
            continue
        assert filecmp.cmp(os.path.join(d_atomic, f), b, shallow=False), \
            "{0} differs between the atomic and split initialize".format(f)


def _stack_case(name, risk=0.95, stack_size=10, chance_points="all", tool="mou"):
    """A case that uses STACKS for chance, built from the tracked g07 template.

    g07's parameters are ALL decision variables, and a stack needs something else to vary -
    "adjustable but not a decision variable" is what gets drawn. So one parameter is moved to
    its own group, which is the same shape the rosenbrock chance case uses.

    sqp needs ``sqp_num_reals`` for any of this: without an ensemble it takes the
    finite-difference gradient path, where chance dies with "no stack runs have been
    processed". That is a real gap, not a configuration rule.
    """
    base = os.path.join(_BENCH, "g07", "template")
    d = os.path.join(_BENCH, name)
    if os.path.exists(d):
        shutil.rmtree(d)
    shutil.copytree(base, d)
    pst = pyemu.Pst(os.path.join(d, "g07.pst"))
    drop = "sqp_" if tool == "mou" else "mou_"
    for k in [k for k in pst.pestpp_options if k.startswith(drop)]:
        pst.pestpp_options.pop(k)
    # the uncertain parameter the stack is drawn over
    pst.parameter_data.loc["x10", "pargp"] = "uncertain"
    o = pst.pestpp_options
    o["opt_dec_var_groups"] = "decvar"
    o["opt_stack_size"] = stack_size
    o["opt_chance_points"] = chance_points
    o["random_seed"] = 11
    o["opt_risk"] = risk
    if tool == "mou":
        o["mou_population_size"] = 10
        o["mou_generator"] = "de"
    else:
        o["sqp_num_reals"] = 8
    pst.control_data.noptmax = 1
    pst.write(os.path.join(d, "g07.pst"), version=2)
    return d


def api_chance_stacks_are_reachable_test():
    """The chance stacks, live, for mou.

    mou and sqp carry uncertainty through a PARAMETER STACK that is actually run, and until now
    the stack was unreachable: stack_pe/stack_oe are private members of Constraints. That made
    the one thing a surrogate-assisted or bayesian-optimization workflow needs - the per-member
    uncertainty behind each objective value - invisible from python.

    The shapes are the assertion that matters. stack_oe must have one row per stack realization
    and one column per observation, and each per-member stack must match it, because they are
    the same stack evaluated at a different point in decision variable space.
    """
    d = _stack_case("api_chance_stacks", risk=0.95, stack_size=10)
    with Mou.from_pst("g07.pst", workdir=d) as mou:
        mou.initialize()

        st = mou.stack_status
        assert st["use_chance"], "risk is 0.95, so chance must be on: {0}".format(st)
        assert not st["use_fosm"], \
            "opt_stack_size is set, so stacks must beat fosm: {0}".format(st)
        assert abs(st["risk"] - 0.95) < 1.0e-8, st
        assert st["stack_size"] == 10, \
            "the drawn stack should have 10 realizations, got {0}".format(st["stack_size"])

        # the derived flags must agree with the status dict they are drawn from
        assert mou.risk == st["risk"] and mou.use_chance == st["use_chance"] \
            and mou.use_fosm == st["use_fosm"] and mou.use_robust == st["use_robust"]

        pe, oe = mou.stack_pe(), mou.stack_oe()
        assert pe.shape[0] == 10, "stack_pe should have 10 realizations, got {0}".format(pe.shape)
        assert oe.shape[0] == 10, "stack_oe should have 10 realizations, got {0}".format(oe.shape)
        assert list(oe.columns) == [n.upper() for n in mou.pst.obs_names], \
            "stack_oe should span the observations: {0}".format(list(oe.columns))
        assert not oe.isnull().values.any(), "the stack was run, so it should hold values"

        # opt_chance_points="all" means every member carries its own stack
        names = mou.member_stacks()
        assert len(names) > 0, "opt_chance_points is 'all', so there should be member stacks"
        one = mou.member_stack_oe(names[0])
        assert one.shape[0] == oe.shape[0], \
            "a member stack has the same realizations as the stack it came from: {0} vs {1}"\
            .format(one.shape, oe.shape)
        # but NOT the same columns: the per-member stacks carry the CONSTRAINTS only, while
        # the single stack spans the objective too. Worth pinning - it is invisible on a case
        # where every observation happens to be a constraint.
        assert set(one.columns) < set(oe.columns), \
            "a member stack should be a strict column subset of the single stack: {0} vs {1}"\
            .format(list(one.columns), list(oe.columns))
        # name and index address the same stack - the ids walk the same sorted map
        assert (mou.member_stack_oe(0).values == one.values).all(), \
            "member_stack_oe by index and by name disagree"

        try:
            mou.member_stack_oe("not-a-member")
            raise AssertionError("a bad member name should raise")
        except KeyError:
            pass
        try:
            mou.member_stack_oe(len(names))
            raise AssertionError("an out of range member index should raise")
        except IndexError:
            pass

        # lower= must reach the columns, as everywhere else in the friendly layer
        assert list(mou.stack_oe(lower=True).columns) == [c.lower() for c in oe.columns]

        mou.finalize()
    print("api_chance_stacks_are_reachable_test passed")


def api_empty_stack_says_why_test():
    """An empty stack is an ANSWER, not a failure - and the reason is retrievable.

    A fosm run and a risk-neutral run both leave every stack empty forever, and neither is an
    error. But "empty" alone is indistinguishable from "not drawn yet", and those want opposite
    responses, so the status flags have to separate them. This asserts the fosm case: chance is
    on, the stacks are empty anyway, and use_fosm says why.
    """
    # risk on, but no stack configured -> fosm, and the stacks stay empty
    d = _stack_case("api_chance_fosm", risk=0.95, stack_size=0, chance_points="single")
    with Mou.from_pst("g07.pst", workdir=d) as mou:
        mou.initialize()
        st = mou.stack_status
        assert st["use_chance"], "risk is still 0.95, so chance is on: {0}".format(st)
        assert st["use_fosm"], \
            "with no stack configured this must fall back to fosm: {0}".format(st)
        assert mou.stack_pe().empty and mou.stack_oe().empty, \
            "a fosm run draws no stack, so both stacks must be empty"
        assert mou.member_stacks() == [], "fosm has no per-member stacks either"
        # and asking for one names the reason rather than raising something opaque
        try:
            mou.member_stack_oe(0)
            raise AssertionError("there are no member stacks, so this should raise")
        except RuntimeError as e:
            assert "fosm" in str(e).lower(), \
                "the error should say WHY there are no stacks, got: {0}".format(e)
        mou.finalize()

    # risk neutral: chance off entirely, whatever else is configured
    d = _stack_case("api_chance_neutral", risk=0.5, stack_size=10)
    with Mou.from_pst("g07.pst", workdir=d) as mou:
        mou.initialize()
        st = mou.stack_status
        assert not st["use_chance"], \
            "risk 0.5 is risk neutral and switches chance off: {0}".format(st)
        assert mou.stack_oe().empty, "a risk neutral run draws no stack"
        assert mou.risk == 0.5
        mou.finalize()
    print("api_empty_stack_says_why_test passed")


def api_risk_is_clamped_and_live_test():
    """set_risk() lands where the tool reads it, and the tool clamps what it uses.

    The clamp is the point: opt_risk stays at whatever was set, while the value the tool
    actually uses is bounded to [0.001, 0.999]. Reading the option instead of `risk` gives a
    number the tool never uses, which is why `risk` is derived on the C++ side rather than
    recomputed here.
    """
    d = _stack_case("api_chance_clamp", risk=0.95, stack_size=10)
    with Mou.from_pst("g07.pst", workdir=d) as mou:
        mou.initialize()
        mou.set_risk(1.5)
        # get_option is string-valued by contract, so this is the raw text of what was set
        assert float(mou.get_option("opt_risk")) == 1.5, "the option should hold what was set"
        assert mou.risk == 0.999, \
            "the tool clamps to 0.999, got {0}".format(mou.risk)
        mou.set_risk(0.0)
        assert mou.risk == 0.001, "clamped at the bottom too, got {0}".format(mou.risk)
        mou.set_risk(0.5)
        assert not mou.use_chance, "0.5 switches chance back off"
        mou.finalize()
    print("api_risk_is_clamped_and_live_test passed")


def api_stacks_are_mou_and_sqp_only_test():
    """ies and da have no chance machinery, and say so rather than returning empties."""
    d = _case("api_chance_not_ies", noptmax=1)
    with Ies.from_pst("pest.pst", workdir=d) as ies:
        ies.initialize()
        # the friendly layer does not put the mixin on ies at all, so this is an attribute
        # error rather than a runtime one - the surface simply is not there
        assert not hasattr(ies, "stack_oe"), \
            "ies should not carry the chance surface at all"
        # and the ABI underneath refuses rather than inventing an empty stack
        try:
            ies._lib.get_stack_status()
            raise AssertionError("ies has no constraints, so this should refuse")
        except PestppError as e:
            assert "chance" in str(e).lower() or "constraint" in str(e).lower(), \
                "the refusal should name what is missing, got: {0}".format(e)
        ies.finalize()
    print("api_stacks_are_mou_and_sqp_only_test passed")


def api_sqp_chance_stacks_are_reachable_test():
    """The same stack surface on sqp, driven by the SAME option as mou and pestpp-opt.

    sqp used to reach chance through ``sqp_risk``, which was both a risk value and a switch
    selecting a different mechanism, and which silently overrode ``opt_risk``. It is retired -
    there is one risk option now, and ``use_robust`` names the other mechanism explicitly.
    """
    d = _stack_case("api_chance_sqp", risk=0.95, stack_size=8, chance_points="SINGLE",
                    tool="sqp")
    with Sqp.from_pst("g07.pst", workdir=d) as sqp:
        sqp.initialize()

        st = sqp.stack_status
        assert not st["use_robust"], "this is a chance run, not a robust one: {0}".format(st)
        assert st["use_chance"] and not st["use_fosm"], \
            "opt_risk with a stack means stack-based chance: {0}".format(st)
        assert abs(st["risk"] - 0.95) < 1.0e-8, st

        pe, oe = sqp.stack_pe(), sqp.stack_oe()
        assert pe.shape[0] == 8, "stack_pe should hold 8 realizations, got {0}".format(pe.shape)
        assert oe.shape[0] == 8, "stack_oe should hold 8 realizations, got {0}".format(oe.shape)
        assert not oe.isnull().values.any(), \
            "the stack runs were made, so the obs stack should hold values"
        assert sqp.member_stacks() == [], \
            "SINGLE chance points means no per-member stacks"
        sqp.finalize()
    print("api_sqp_chance_stacks_are_reachable_test passed")


def api_robust_does_no_risk_shifting_test():
    """opt_use_robust turns the chance machinery OFF rather than choosing a flavour of it.

    This is the assertion the whole refactor rests on. Robust optimization pairs each
    decision-variable realization with its own uncertain parameter realization and optimizes
    the ensemble as it stands - so there is nothing to shift, no stack to draw, and no runs to
    spend on one. If any of those were still happening, every downstream "shift or not" branch
    would need to know about robust mode; because they are not, none of them do.
    """
    d = _stack_case("api_robust_sqp", risk=0.5, stack_size=8, chance_points="SINGLE",
                    tool="sqp")
    pst = pyemu.Pst(os.path.join(d, "g07.pst"))
    pst.pestpp_options["opt_use_robust"] = True
    pst.write(os.path.join(d, "g07.pst"), version=2)

    with Sqp.from_pst("g07.pst", workdir=d) as sqp:
        sqp.initialize()
        st = sqp.stack_status
        assert st["use_robust"], "opt_use_robust was set: {0}".format(st)
        assert not st["use_chance"], \
            "robust does no risk shifting, so chance must be off: {0}".format(st)
        assert not st["use_fosm"], "robust does not fall back to fosm either: {0}".format(st)
        assert st["stack_size"] == 0, \
            "a robust run should not draw a stack at all, got {0}".format(st["stack_size"])
        assert sqp.stack_pe().empty and sqp.stack_oe().empty, \
            "no stack means both stack ensembles are empty"

        # the pairing: the dv ensemble carries the uncertain parameter alongside the dec vars
        cols = set(sqp.par_df().columns)
        assert "X10" in cols, \
            "robust pairs each dv realization with an uncertain parameter draw, so the "\
            "uncertain parameter should be in the ensemble: {0}".format(sorted(cols))

        # and set_risk is refused rather than silently doing nothing
        try:
            sqp.set_risk(0.95)
            raise AssertionError("set_risk should be refused on a robust run")
        except RuntimeError as e:
            assert "robust" in str(e).lower(), e
        sqp.finalize()
    print("api_robust_does_no_risk_shifting_test passed")


def api_robust_and_risk_are_mutually_exclusive_test():
    """Setting both is refused, rather than one silently winning.

    Silent precedence between the two mechanisms is exactly what retiring sqp_risk removed, so
    reintroducing it here would defeat the point: a robust run ignores opt_risk completely, and
    a user who sets both has asked for two different things.
    """
    d = _stack_case("api_robust_conflict", risk=0.95, stack_size=8, chance_points="SINGLE",
                    tool="sqp")
    pst = pyemu.Pst(os.path.join(d, "g07.pst"))
    pst.pestpp_options["opt_use_robust"] = True     # ...alongside opt_risk 0.95
    pst.write(os.path.join(d, "g07.pst"), version=2)
    try:
        with Sqp.from_pst("g07.pst", workdir=d) as sqp:
            sqp.initialize()
        raise AssertionError("opt_use_robust with a non-neutral opt_risk should be refused")
    except PestppError as e:
        assert "mutually exclusive" in str(e).lower(), \
            "the refusal should say why, got: {0}".format(e)
    print("api_robust_and_risk_are_mutually_exclusive_test passed")


def api_mou_robust_pairs_and_records_test():
    """mou supports opt_use_robust: each member gets one stack realization, and it is recorded.

    This test replaced one asserting mou REJECTED the flag - which it did until robust support
    landed. CI caught the staleness on all three platforms, which is the useful kind of failure.

    Three things have to hold together, and only the third is hard to get right:

      1. `dp` stays decision-variable only. The pairing rides the per-realization FixedParInfo
         side channel instead, so the generators cannot cross and mutate parameter realizations
         into arithmetic averages that are draws from nothing.
      2. Every generation is recorded, INCLUDING generation 0 - the initial population is the
         one every later generation descends from, and leaving it unpaired would mean the run
         started from members carrying no uncertainty at all.
      3. The values actually reach the model. A record can look perfect while every member runs
         at its control value, so the stack here is tagged with unmistakable marker values and
         the model's own input file is read back.
    """
    d = _stack_case("api_mou_robust", risk=0.5, stack_size=0, chance_points="SINGLE",
                    tool="mou")
    pst = pyemu.Pst(os.path.join(d, "g07.pst"))
    names = [p_.upper() for p_ in pst.par_names]
    ctl_x10 = float(pst.parameter_data.loc["x10", "parval1"])
    # markers nothing else in g07 comes near, so a hit cannot be coincidence
    markers = [1000.0 * (i + 1) for i in range(5)]
    rng = np.random.default_rng(3)
    cols = {c: rng.uniform(-1.0, 1.0, 5) for c in names}
    cols["X10"] = np.array(markers)
    pd.DataFrame(cols, index=["s{0}".format(i) for i in range(5)]).to_csv(
        os.path.join(d, "unc_stack.csv"))
    pst.parameter_data.loc["x10", "parubnd"] = 1.0e30     # the markers must be legal values
    pst.pestpp_options["opt_use_robust"] = True
    pst.pestpp_options["opt_risk"] = 0.5
    pst.pestpp_options["opt_par_stack"] = "unc_stack.csv"
    pst.pestpp_options["mou_population_size"] = 6
    pst.control_data.noptmax = 2
    pst.write(os.path.join(d, "g07.pst"), version=2)

    with Mou.from_pst("g07.pst", workdir=d) as mou:
        mou.initialize()
        for _ in mou.iterations():
            pass
        mou.finalize()

    rec = os.path.join(d, "g07.par_real.csv")
    assert os.path.exists(rec), "robust mou should write the par-real record"
    df = pd.read_csv(rec)
    assert list(df.columns) == ["member", "generation", "par_real"], list(df.columns)
    assert df.member.is_unique, "member names are globally unique, so the record joins cleanly"
    gens = sorted(df.generation.unique())
    assert 0 in gens, \
        "generation 0 must be paired too - it is what every later generation descends from: "\
        "{0}".format(gens)
    assert len(gens) >= 2, "each generation's new members should be recorded: {0}".format(gens)

    # 3. the values reached the model: par.dat is what it was handed on the final run
    par_dat = os.path.join(d, "par.dat")
    assert os.path.exists(par_dat), "expected the model input file to be left behind"
    vals = [float(x) for x in open(par_dat).read().split("\n")[1:] if x.strip()]
    x10 = vals[9]                                  # par.tpl writes x1..x10 in order
    assert any(abs(x10 - m) < 1.0e-6 for m in markers), \
        "x10 reached the model as {0}, which is not one of the stack realizations {1} - the "\
        "pairing is being recorded but not applied".format(x10, markers)
    assert abs(x10 - ctl_x10) > 1.0e-6, \
        "x10 is still at its control value, so no uncertainty was propagated"
    print("api_mou_robust_pairs_and_records_test passed")


def api_mou_robust_needs_a_stack_test():
    """Robust mou without a stack is refused, naming both ways to supply one.

    The pairing IS the mechanism, so there is nothing sensible to do without realizations to
    pair with - and silently running every member at its control value would look like a
    working robust run while propagating no uncertainty whatsoever.
    """
    d = _stack_case("api_mou_robust_nostack", risk=0.5, stack_size=0, chance_points="SINGLE",
                    tool="mou")
    pst = pyemu.Pst(os.path.join(d, "g07.pst"))
    pst.pestpp_options["opt_use_robust"] = True
    pst.pestpp_options["opt_risk"] = 0.5
    pst.write(os.path.join(d, "g07.pst"), version=2)
    try:
        with Mou.from_pst("g07.pst", workdir=d) as mou:
            mou.initialize()
        raise AssertionError("robust mou with no stack should be refused")
    except PestppError as e:
        msg = str(e).lower()
        assert "opt_stack_size" in msg and "opt_par_stack" in msg, \
            "the refusal should name both ways to supply a stack, got: {0}".format(e)
    print("api_mou_robust_needs_a_stack_test passed")


def api_retired_sqp_risk_is_refused_test():
    """A control file still saying sqp_risk gets told what to use instead.

    There is no behaviour-preserving translation - the mode sqp_risk selected did shift, and
    robust mode does not - so guessing on the user's behalf would silently change results.
    """
    d = _stack_case("api_retired_sqp_risk", risk=0.5, stack_size=8, chance_points="SINGLE",
                    tool="sqp")
    pst = pyemu.Pst(os.path.join(d, "g07.pst"))
    pst.pestpp_options["sqp_risk"] = 0.95
    pst.write(os.path.join(d, "g07.pst"), version=2)
    try:
        with Sqp.from_pst("g07.pst", workdir=d) as sqp:
            sqp.initialize()
        raise AssertionError("a retired option should be refused")
    except PestppError as e:
        msg = str(e).lower()
        assert "retired" in msg, "the error should say it is retired, got: {0}".format(e)
        assert ("opt_risk" in msg) and ("opt_use_robust" in msg), \
            "the error should name BOTH replacements, got: {0}".format(e)
    print("api_retired_sqp_risk_is_refused_test passed")


def api_init_only_options_are_refused_live_test():
    """init-only options are refused once the tool is running, not silently ignored.

    opt_use_robust decides how the dv ensemble is BUILT. Accepting it afterwards would leave a
    robust ensemble being evaluated as a chance run, or the reverse - a disagreement between
    the options and the object that no later call could detect.
    """
    d = _stack_case("api_init_only", risk=0.95, stack_size=8, chance_points="SINGLE",
                    tool="sqp")
    with Sqp.from_pst("g07.pst", workdir=d) as sqp:
        # before initialize the option is settable - it has not been consumed yet
        sqp.set_option("opt_use_robust", False)
        sqp.initialize()
        for key in ("opt_use_robust", "sqp_num_reals"):
            try:
                sqp.set_option(key, 4 if key == "sqp_num_reals" else True)
                raise AssertionError("'{0}' is init-only and should be refused".format(key))
            except PestppError as e:
                assert "init-only" in str(e).lower(), \
                    "the refusal should say why, got: {0}".format(e)
        # a live option still works, so this is not a blanket freeze
        sqp.set_option("sqp_subset_size", 3)
        sqp.finalize()
    print("api_init_only_options_are_refused_live_test passed")


def _sqp_dv_file(d, pst, path, include_uncertain, n=8, seed=11):
    """Write a dv ensemble file, with or without the uncertain parameter columns."""
    rng = np.random.default_rng(seed)
    names = [p.upper() for p in pst.par_names]
    keep = names if include_uncertain else [n_ for n_ in names if n_ != "X10"]
    lo = dict(zip(names, pst.parameter_data.parlbnd.values))
    hi = dict(zip(names, pst.parameter_data.parubnd.values))
    df = pd.DataFrame({c: rng.uniform(lo[c], hi[c], n) for c in keep},
                      index=["real{0}".format(i) for i in range(n)])
    full = os.path.join(d, path)
    if path.lower().endswith(".csv"):
        df.to_csv(full)
    else:
        pyemu.ParameterEnsemble(pst=pst, df=df).to_binary(full)
    return path


def api_sqp_ensemble_sources_test():
    """sqp reaches the same state from every way of supplying its ensembles.

    Four independent choices, and they have to compose: the dv ensemble may be DRAWN or read
    from CSV or from BINARY, and the stack may be DRAWN or read from a par-stack file. A tool
    that only works when everything is drawn is a tool that cannot be restarted or driven from
    a previous run's output, so the combinations are the contract, not the individual paths.
    """
    for dv_src in ("drawn", "csv", "jcb"):
        for stack_src in ("drawn", "file"):
            name = "api_src_{0}_{1}".format(dv_src, stack_src)
            d = _stack_case(name, risk=0.95, stack_size=8, chance_points="SINGLE", tool="sqp")
            pst = pyemu.Pst(os.path.join(d, "g07.pst"))
            if dv_src != "drawn":
                pst.pestpp_options["sqp_dv_en"] = _sqp_dv_file(
                    d, pst, "dv.csv" if dv_src == "csv" else "dv.jcb", include_uncertain=False)
            if stack_src == "file":
                rng = np.random.default_rng(7)
                names = [p.upper() for p in pst.par_names]
                pd.DataFrame({c: rng.uniform(-1.0, 1.0, 6) for c in names},
                             index=["s{0}".format(i) for i in range(6)]).to_csv(
                    os.path.join(d, "par_stack.csv"))
                pst.pestpp_options["opt_par_stack"] = "par_stack.csv"
            pst.write(os.path.join(d, "g07.pst"), version=2)

            with Sqp.from_pst("g07.pst", workdir=d) as sqp:
                sqp.initialize()
                st = sqp.stack_status
                assert st["use_chance"] and not st["use_fosm"], \
                    "{0}: opt_risk with a stack means stack chance: {1}".format(name, st)
                # a file-supplied stack brings its OWN row count - it is not opt_stack_size
                want = 6 if stack_src == "file" else 8
                assert st["stack_size"] == want, \
                    "{0}: expected a {1}-realization stack, got {2}".format(
                        name, want, st["stack_size"])
                assert sqp.stack_oe().shape[0] == want, \
                    "{0}: the obs stack should match".format(name)
                assert not sqp.stack_oe().isnull().values.any(), \
                    "{0}: the stack was run, so it should hold values".format(name)
                sqp.finalize()
    print("api_sqp_ensemble_sources_test passed")


def api_robust_requires_paired_parameters_test():
    """A supplied dv ensemble must carry the uncertain parameters when robust is on.

    Robust optimisation IS the pairing, so an ensemble of decision variables alone has nothing
    to pair with. Unchecked this fails silently rather than loudly: the uncertain columns get
    filled from one replicated row, every realization sees identical parameter values, and the
    run reports itself robust while optimising against no uncertainty whatsoever. Measured
    before the guard existed, the uncertain parameter came back with a single distinct value
    across every non-base realization.

    The drawn path cannot hit this - it concatenates the two draws - so this is specifically
    about user-supplied ensembles.
    """
    # drawn: sqp pairs for you, and the uncertain parameter is in the ensemble
    d = _stack_case("api_robust_drawn", risk=0.5, stack_size=0, chance_points="SINGLE",
                    tool="sqp")
    pst = pyemu.Pst(os.path.join(d, "g07.pst"))
    pst.pestpp_options["opt_use_robust"] = True
    pst.write(os.path.join(d, "g07.pst"), version=2)
    with Sqp.from_pst("g07.pst", workdir=d) as sqp:
        sqp.initialize()
        df = sqp.par_df()
        assert "X10" in df.columns, "the drawn path should pair the uncertain parameter in"
        assert df["X10"].nunique() > 1, \
            "a paired draw must actually VARY across realizations, got {0} distinct".format(
                df["X10"].nunique())
        sqp.finalize()

    # supplied WITH the uncertain columns: accepted, and it still varies
    for ext in ("csv", "jcb"):
        d = _stack_case("api_robust_paired_" + ext, risk=0.5, stack_size=0,
                        chance_points="SINGLE", tool="sqp")
        pst = pyemu.Pst(os.path.join(d, "g07.pst"))
        pst.pestpp_options["opt_use_robust"] = True
        pst.pestpp_options["sqp_dv_en"] = _sqp_dv_file(d, pst, "dv." + ext,
                                                       include_uncertain=True)
        pst.write(os.path.join(d, "g07.pst"), version=2)
        with Sqp.from_pst("g07.pst", workdir=d) as sqp:
            sqp.initialize()
            assert sqp.par_df()["X10"].nunique() > 1, \
                "{0}: the supplied pairing should be used as given".format(ext)
            sqp.finalize()

    # supplied WITHOUT them: refused, naming the parameter and both ways out
    d = _stack_case("api_robust_unpaired", risk=0.5, stack_size=0, chance_points="SINGLE",
                    tool="sqp")
    pst = pyemu.Pst(os.path.join(d, "g07.pst"))
    pst.pestpp_options["opt_use_robust"] = True
    pst.pestpp_options["sqp_dv_en"] = _sqp_dv_file(d, pst, "dv_only.csv",
                                                   include_uncertain=False)
    pst.write(os.path.join(d, "g07.pst"), version=2)
    try:
        with Sqp.from_pst("g07.pst", workdir=d) as sqp:
            sqp.initialize()
        raise AssertionError("a dv-only ensemble under robust should be refused")
    except PestppError as e:
        msg = str(e)
        assert "X10" in msg.upper(), \
            "the error should name the missing parameter, got: {0}".format(e)
        assert "sqp_dv_en" in msg, \
            "the error should say how to fix it, got: {0}".format(e)
    print("api_robust_requires_paired_parameters_test passed")


def api_run_manager_settings_are_live_test():
    """Run-manager tuning set through the API reaches the RUN MANAGER, not just the options.

    RunManagerPanther is deliberately decoupled from PestppOptions - it never reads the
    scenario, every value arrives through its constructor - and the API builds it during
    create(), before a caller can touch anything. So a set_option afterwards used to reach the
    options object and stop there: the call returned OK, get_option read the new value back,
    and the master went on using the value it was built with. Nothing anywhere said otherwise.

    These four are consulted inside the scheduling loop on every pass, so changing one mid-run
    changes what happens to the runs already in flight - which is the point, since "this batch
    is dragging, give up sooner" is a decision you can only make once you can see it dragging.
    """
    d = _case("api_rm_live", noptmax=1, num_reals=5)
    pst = pyemu.Pst(os.path.join(d, "pest.pst"))
    pst.pestpp_options["overdue_giveup_fac"] = 2.0
    pst.pestpp_options["overdue_resched_fac"] = 1.15
    pst.pestpp_options["overdue_giveup_minutes"] = 3000.0
    pst.pestpp_options["max_run_fail"] = 3
    pst.write(os.path.join(d, "pest.pst"), version=2)

    with Ies.from_pst("pest.pst", workdir=d, workers=2, port=4601) as ies:
        ies.initialize()
        built = ies.run_manager_settings()
        assert built["overdue_giveup_fac"] == 2.0, \
            "the manager should be built with what the control file asked for: {0}".format(built)
        assert built["max_run_fail"] == 3, built

        # the assertion that matters: the RUN MANAGER changes, not merely the option
        ies.set_option("overdue_giveup_fac", 99.0)
        ies.set_option("overdue_resched_fac", 7.5)
        ies.set_option("overdue_giveup_minutes", 11.0)
        ies.set_option("max_run_fail", 9)
        live = ies.run_manager_settings()
        assert live["overdue_giveup_fac"] == 99.0, \
            "set_option must reach the run manager, got {0}".format(live)
        assert live["overdue_resched_fac"] == 7.5, live
        assert live["overdue_giveup_minutes"] == 11.0, live
        assert live["max_run_fail"] == 9, live
        # ...and the options agree, so the two views cannot silently diverge
        assert float(ies.get_option("overdue_giveup_fac")) == 99.0

        # the construction-bound panther options are refused rather than silently ignored,
        # which is the honest half of the same problem
        for key, val in (("panther_persistent_workers", False),
                         ("panther_master_timeout_milliseconds", 50),
                         ("panther_ping_interval_secs", 15)):
            try:
                ies.set_option(key, val)
                raise AssertionError("'{0}' is consumed at construction and should be "
                                     "refused on a running session".format(key))
            except PestppError as e:
                assert "init-only" in str(e).lower(), \
                    "the refusal should say why, got: {0}".format(e)
        ies.finalize()
    print("api_run_manager_settings_are_live_test passed")


def api_overdue_policy_is_panther_only_test():
    """A serial session has no overdue policy, and says so rather than inventing numbers."""
    d = _case("api_rm_serial", noptmax=1, num_reals=5)
    with Ies.from_pst("pest.pst", workdir=d) as ies:      # serial: no workers
        ies.initialize()
        # max_run_fail is shared by every run manager, so it still answers
        assert ies.run_manager_settings(overdue=False)["max_run_fail"] >= 0
        try:
            ies.run_manager_settings()
            raise AssertionError("the overdue policy should be refused on a serial session")
        except PestppError as e:
            assert "panther" in str(e).lower(), \
                "the refusal should name the manager that has one, got: {0}".format(e)
        ies.finalize()
    print("api_overdue_policy_is_panther_only_test passed")


def api_missing_par_columns_get_control_values_test():
    """A parameter absent from a supplied ensemble gets its parval1, not zero.

    ParameterEnsemble::from_csv sets var_names to EVERY control-file parameter and then
    read_csv_by_reals does resize()+setZero(), so a parameter the csv did not mention is not
    absent from the ensemble - it is present, and it is ZERO. fill_fixed() had long since fixed
    this for FIXED parameters; adjustable ones got nothing.

    Zero is a perfectly plausible parameter value, so nothing downstream could notice: the
    model simply ran with it. On a log-transformed parameter it is not even plausible - it is
    out of bounds and becomes NaN once transformed.

    Only two callers forgive a missing column, and the second is the one that matters:
    sqp_dv_en, and opt_par_stack - which mou, sqp AND pestpp-opt all use to supply a chance
    stack. A stack realization silently carrying 0 feeds straight into the risk shift.

    Measured on this case before the fix: every stack realization had X10 = 0.0, against a
    parval1 of 5.07.
    """
    d = _stack_case("api_fill_stack", risk=0.95, stack_size=6, chance_points="SINGLE",
                    tool="mou")
    pst = pyemu.Pst(os.path.join(d, "g07.pst"))
    names = [p_.upper() for p_ in pst.par_names]
    drop = "X10"                                # the uncertain parameter the stack is over
    parval = float(pst.parameter_data.loc[drop.lower(), "parval1"])
    assert parval != 0.0, "this test is meaningless if parval1 is itself zero"

    rng = np.random.default_rng(5)
    keep = [n for n in names if n != drop]
    pd.DataFrame({c: rng.uniform(-1.0, 1.0, 6) for c in keep},
                 index=["s{0}".format(i) for i in range(6)]).to_csv(
        os.path.join(d, "par_stack.csv"))
    pst.pestpp_options["opt_par_stack"] = "par_stack.csv"
    pst.write(os.path.join(d, "g07.pst"), version=2)

    with Mou.from_pst("g07.pst", workdir=d) as mou:
        mou.initialize()
        pe = mou.stack_pe()
        assert drop in pe.columns, \
            "the column exists either way - that is what made this invisible"
        vals = pe[drop].values
        assert not (vals == 0.0).any(), \
            "{0} was left at the allocation zero rather than its control value".format(drop)
        assert np.allclose(vals, parval), \
            "{0} should take parval1 ({1}) in every realization, got {2}".format(
                drop, parval, sorted(set(vals))[:4])
        # the fill must not blanket-overwrite; the columns the file supplied still hold real
        # values. Note they are CONSTANT here and legitimately so - a stack replaces the
        # decision-variable entries of every realization with the current dv values, so only
        # the uncertain parameter actually varies down a stack column.
        supplied = pe[keep[0]].values
        assert not (supplied == 0.0).all(), \
            "a column the file supplied should not have been zeroed"
        mou.finalize()
    print("api_missing_par_columns_get_control_values_test passed")


def api_sqp_needs_an_ensemble_gradient_test():
    """sqp runs out of the box, and refuses clearly when it cannot form a gradient.

    PESTPP-SQP implements the ensemble gradient ONLY - `iterate_2_solution()` throws "only
    ensemble gradient currently implemented" for anything else, and the finite-difference entry
    point it names exists solely as a commented-out call, declared and defined nowhere. With
    sqp_num_reals defaulting to -1, that made the DEFAULT configuration the unimplemented one:
    running sqp without setting sqp_num_reals or sqp_dv_en failed outright.

    The default is now 50, so a plain run works. Switching the ensemble off explicitly is still
    a legitimate thing to ask for and still cannot be served, so it is refused during prepare -
    with a message naming what to set, rather than failing later in whatever machinery happens
    to run first. With chance enabled that used to surface as a complaint about stack runs
    never having been processed, which points nowhere near the actual problem.
    """
    # no sqp ensemble options at all: the default has to be a workable one
    d = _stack_case("api_sqp_default_grad", risk=0.5, stack_size=0, chance_points="SINGLE",
                    tool="sqp")
    pst = pyemu.Pst(os.path.join(d, "g07.pst"))
    pst.pestpp_options.pop("sqp_num_reals", None)
    pst.write(os.path.join(d, "g07.pst"), version=2)
    with Sqp.from_pst("g07.pst", workdir=d) as sqp:
        assert int(sqp.get_option("sqp_num_reals")) == 50, \
            "sqp_num_reals should default to 50, got {0}".format(
                sqp.get_option("sqp_num_reals"))
        sqp.initialize()
        assert sqp.n_reals > 1, \
            "the default should give sqp a real ensemble to work from, got {0}".format(
                sqp.n_reals)
        sqp.finalize()

    # explicitly switched off: refused, and the message says what to set
    d = _stack_case("api_sqp_no_grad", risk=0.5, stack_size=0, chance_points="SINGLE",
                    tool="sqp")
    pst = pyemu.Pst(os.path.join(d, "g07.pst"))
    pst.pestpp_options["sqp_num_reals"] = 0
    pst.write(os.path.join(d, "g07.pst"), version=2)
    try:
        with Sqp.from_pst("g07.pst", workdir=d) as sqp:
            sqp.initialize()
        raise AssertionError("sqp with no ensemble gradient should be refused")
    except PestppError as e:
        msg = str(e).lower()
        assert "sqp_num_reals" in msg and "sqp_dv_en" in msg, \
            "the refusal should name both ways to fix it, got: {0}".format(e)
        assert "finite-difference" in msg, \
            "the refusal should say finite differences are not implemented, got: {0}".format(e)
    print("api_sqp_needs_an_ensemble_gradient_test passed")


if __name__ == "__main__":
    api_smoke_test()
    api_iterations_respect_noptmax_test()
    api_par_df_roundtrip_test()
    api_par_view_is_live_test()
    api_view_invalidated_by_resize_test()
    api_own_the_initial_batch_test()
    api_phi_across_realizations_test()
    api_phi_by_obs_group_test()
    api_weights_test()
    api_weights_reject_bad_input_test()
    api_name_case_is_absorbed_test()
    api_drop_realizations_test()
    api_quiet_captures_output_test()
    api_tools_without_phi_say_so_test()
    api_run_ies_one_call_test()
    api_find_library_test()
    api_notebook_runs_test()
    api_pyemu_objects_test()
    api_pyemu_weights_stay_live_test()
    api_pyemu_bounds_enforced_test()
    api_pyemu_from_pst_object_test()
    api_pyemu_results_test()
    api_activate_zero_weighted_obs_test()
    api_reweight_leaves_noise_alone_test()
    api_reinflate_size_is_capped_by_prior_test()
    api_reinflate_grows_from_bigger_prior_test()
    api_reinflate_centring_is_per_call_test()
    api_reinflate_keeps_caller_weights_test()
    api_reinflate_grows_after_activation_test()
    api_option_values_are_native_types_test()
    api_option_sequence_reaches_the_algorithm_test()
    api_deferred_solve_drives_the_runs_test()
    api_deferred_solve_matches_composed_test()
    api_deferred_solve_edit_reaches_the_runs_test()
    api_deferred_solve_subset_override_test()
    api_deferred_solve_mou_test()
    api_deferred_solve_refused_where_it_does_not_fit_test()
    api_deferred_solve_state_guards_test()
    api_mou_runs_past_noptmax_test()
    api_deferred_solve_needs_upgrades_in_memory_test()
    api_progress_renders_in_a_terminal_test()
    api_progress_is_inert_when_disabled_test()
    api_progress_renders_in_a_notebook_test()
    api_progress_during_a_real_run_test()
    api_run_observer_fires_during_a_composed_solve_test()
    api_run_observer_is_throttled_in_the_library_test()
    api_run_observer_reentrancy_is_an_allowlist_test()
    api_run_observer_can_stop_a_batch_test()
    api_run_observer_survives_a_raising_callback_test()
    api_run_observer_thunk_is_retained_test()
    api_session_cleanup_does_not_need_with_test()
    api_deferred_initialize_mou_test()
    api_deferred_initialize_sqp_test()
    api_chance_stacks_are_reachable_test()
    api_empty_stack_says_why_test()
    api_risk_is_clamped_and_live_test()
    api_sqp_chance_stacks_are_reachable_test()
    api_robust_does_no_risk_shifting_test()
    api_robust_and_risk_are_mutually_exclusive_test()
    api_mou_robust_pairs_and_records_test()
    api_mou_robust_needs_a_stack_test()
    api_retired_sqp_risk_is_refused_test()
    api_init_only_options_are_refused_live_test()
    api_sqp_needs_an_ensemble_gradient_test()
    api_missing_par_columns_get_control_values_test()
    api_run_manager_settings_are_live_test()
    api_overdue_policy_is_panther_only_test()
    api_sqp_ensemble_sources_test()
    api_robust_requires_paired_parameters_test()
    api_stacks_are_mou_and_sqp_only_test()
    print("all helper tests passed")

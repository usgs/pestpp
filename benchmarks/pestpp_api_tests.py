"""Tests for the ergonomic layer, python/pestpp.py.

capi_tests.py covers the C ABI and the thin binding. This file covers the layer users
actually touch: DataFrames, properties, the iteration generator, view scoping, and the
one-call run_* functions.

Worth testing separately because the helper layer can drift from the binding underneath it
without anything failing to compile - a renamed C symbol shows up here as an AttributeError,
and a semantic change (an iteration count, a transform space) shows up nowhere else at all.
"""
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

    nb_path = os.path.join(_REPO, "examples", "pestpp_api_demo.ipynb")
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
    assert not failures, "notebook cells raised:\n  " + "\n  ".join(failures)

    # scratch directories the notebook creates alongside itself
    for name in ("nb_ies", "nb_prior", "nb_cull", "nb_parallel", "nb_parallel_workers",
                 "nb_defer", "nb_progress"):
        shutil.rmtree(os.path.join(nb_dir, name), ignore_errors=True)


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
        # finishing before the runs have been harvested
        ies.queue_runs()
        try:
            ies.finish_solve()
        except PestppError as e:
            assert "not been harvested" in str(e), str(e)
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
    print("all helper tests passed")

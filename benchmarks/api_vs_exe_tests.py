"""Does driving a tool through the API give the same answer as running the executable?

This is the question the rest of the API test suite cannot ask. capi_tests.py checks that the
C ABI behaves sensibly and pestpp_api_tests.py checks the ergonomic layer, but both judge the
API against itself: if an adapter quietly skipped half of a tool's iteration, every one of
those tests would still pass, because the numbers they assert on would be self-consistently
wrong.

So each test here runs the SAME control file with the SAME options twice - once with the
executable, once through the API - into two directories, and compares every output file byte
for byte. Options are deliberately non-default, because an adapter that only handles the easy
path is exactly what this is looking for.

It is not hypothetical. Written against the adapters as they were, these tests found that
mou never advanced MOEA's generation counter (so every generation reported itself as the same
number and wrote over the same files) and that sqp ran only the first statement of its
iteration - skipping the gradient runs, the CMA and hessian updates, the constraint report and
the pcs summary. Both are fixed; these tests are what keeps them fixed.

Byte comparison rather than a tolerance on purpose. The tools are deterministic given a seed,
so "close enough" would hide exactly the small divergence that means an adapter has started
taking a different path.

da is the exception, and its test says so at length: pestpp-da is a multi-cycle driver that
draws a global prior and injects it into each cycle, while the API exposes one cycle that
draws its own. That test pins the gap rather than pretending it is not there.
"""
import os
import platform
import shutil
import subprocess
import sys

import pyemu

_BENCH = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_BENCH)
sys.path.insert(0, os.path.join(_REPO, "python"))
from pestpp import Ies, Da, Mou, Sqp  # noqa: E402

plat = platform.platform().lower()
_sub = "win" if os.name == "nt" else ("mac" if ("darwin" in plat or "macos" in plat) else "linux")
_exe = ".exe" if os.name == "nt" else ""
# absolute, and anchored to this file: the library chdirs while a session is open, and the
# executables are launched with cwd set to the case directory
os.environ["PATH"] += os.pathsep + os.path.join(_BENCH, "test_bin", _sub)


def _find_exe(name):
    """The built executable, wherever this checkout put it."""
    roots = [os.path.join(_REPO, "bin"),
             os.path.join("..", "..", "pestpp", "bin"),
             os.path.join("..", "..", "..", "..", "pestpp", "bin")]
    for root in roots:
        for cand in (os.path.join(root, _sub, name + _exe), os.path.join(root, name + _exe)):
            if os.path.exists(cand):
                return os.path.abspath(cand)
    raise RuntimeError("could not find {0} under any of {1}".format(name, roots))


def _pair(tag, base_d, pst_name, edit):
    """Two identical case directories, one for the exe and one for the API."""
    out = []
    for which in ("exe", "api"):
        d = os.path.join(_BENCH, "{0}_{1}".format(tag, which))
        if os.path.exists(d):
            shutil.rmtree(d)
        shutil.copytree(os.path.join(_BENCH, base_d), d)
        pst = pyemu.Pst(os.path.join(d, pst_name))
        edit(pst)
        pst.write(os.path.join(d, pst_name), version=2)
        out.append(d)
    return out


def _run_exe(name, d, pst_name):
    r = subprocess.run([_find_exe(name), pst_name], cwd=d, capture_output=True, text=True)
    assert r.returncode == 0, "{0} failed:\n{1}".format(name, r.stdout[-3000:])


def _compare_csvs(d_exe, d_api, tag):
    """Every .csv the two runs have in common, plus which files only one of them wrote.

    csv rather than every file because the .rec carries timings and a tool banner, and the
    binary run-storage files carry paths - none of which say anything about whether the
    algorithm did the same thing.
    """
    names_e = {f for f in os.listdir(d_exe) if f.endswith(".csv")}
    names_a = {f for f in os.listdir(d_api) if f.endswith(".csv")}
    differing = []
    for f in sorted(names_e & names_a):
        with open(os.path.join(d_exe, f), "rb") as fp:
            x = fp.read()
        with open(os.path.join(d_api, f), "rb") as fp:
            y = fp.read()
        if x != y:
            differing.append(f)
    return sorted(names_e & names_a), differing, sorted(names_e - names_a), sorted(names_a - names_e)


#: text that can only have come from a model's console output, never from a pest++ csv
_MODEL_OUTPUT_MARKERS = ("MODFLOW", "Solving:  Stress period", "Normal termination",
                         "Run end date and time")


def looks_corrupted(path):
    """Whether a .csv contains something that cannot have been written by pest++.

    Descriptor corruption on windows has twice put a model's console output INSIDE one of the
    library's own output files - once in pest.phi.meas.csv, once in pest.phi.composite.csv.
    Both times the only symptom was "1 of 25 files differ", which reads like a numeric wobble
    and is why it was written off as a flake the first time.

    The exposure is long-lived streams: the phi csvs and the .rec are held open by the
    FileManager for the whole run, and the run storage fstream likewise, so they are open
    across every model spawn. The ensemble csvs are opened, written and closed between runs
    and have never been seen corrupted.

    Returns a reason, or None if the file looks like a csv.
    """
    try:
        with open(path, "r", errors="replace") as fp:
            text = fp.read()
    except OSError as e:
        return "could not be read: {0}".format(e)
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        return None
    for marker in _MODEL_OUTPUT_MARKERS:
        if marker in text:
            return ("contains model console output ({0!r}) - a descriptor the library holds "
                    "open was written to by something else".format(marker))
    # a csv whose rows disagree wildly on column count is not a csv any more
    widths = [ln.count(",") for ln in lines]
    if widths and max(widths) > 0:
        odd = [w for w in widths if w == 0]
        if len(odd) > max(1, len(widths) // 10):
            return "{0} of {1} non-empty lines have no commas at all".format(
                len(odd), len(widths))
    return None


def assert_outputs_well_formed(d, tag):
    """Fail loudly, and by name, if any csv in this directory has been corrupted."""
    bad = []
    for f in sorted(os.listdir(d)):
        if not f.endswith(".csv"):
            continue
        why = looks_corrupted(os.path.join(d, f))
        if why:
            bad.append("{0} {1}".format(f, why))
    assert not bad, (
        "{0}: output file(s) corrupted, which is NOT a numerical difference:\n  {1}".format(
            tag, "\n  ".join(bad)))


def _first_difference(d_exe, d_api, fname, max_lines=3):
    """The first few differing lines of one file, exe against api.

    "these files differ" is not a diagnosis, and on a platform you cannot reproduce locally it
    is the difference between fixing something and guessing at it. A last-digit wobble in one
    column and a structurally different run look identical in the file list and nothing alike
    here.
    """
    out = []
    with open(os.path.join(d_exe, fname)) as fp:
        a = fp.read().splitlines()
    with open(os.path.join(d_api, fname)) as fp:
        b = fp.read().splitlines()
    if len(a) != len(b):
        out.append("      line count differs: exe {0}, api {1}".format(len(a), len(b)))
    for i, (x, y) in enumerate(zip(a, b)):
        if x == y:
            continue
        out.append("      line {0}:".format(i))
        out.append("        exe: {0}".format(x[:300]))
        out.append("        api: {0}".format(y[:300]))
        # name the columns that actually moved, so a one-column wobble is obvious
        xs, ys = x.split(","), y.split(",")
        if len(xs) == len(ys):
            moved = [j for j in range(len(xs)) if xs[j] != ys[j]]
            head = a[0].split(",") if a else []
            named = [(head[j] if j < len(head) else str(j)) for j in moved]
            out.append("        columns differing: {0}".format(named[:12]))
        if len(out) >= max_lines * 4:
            break
    return out


def _assert_identical(d_exe, d_api, tag):
    # checked FIRST and separately: a corrupted file is not the same finding as a different
    # one, and reporting it as "the API took a different path" sends you looking at the
    # algorithm when the problem is a file descriptor
    assert_outputs_well_formed(d_exe, tag + " (exe)")
    assert_outputs_well_formed(d_api, tag + " (api)")
    shared, differing, only_e, only_a = _compare_csvs(d_exe, d_api, tag)
    detail = []
    for f in differing[:3]:
        detail.append("    {0}:".format(f))
        detail.extend(_first_difference(d_exe, d_api, f))
    assert shared, "{0}: the two runs share no csv output at all".format(tag)
    assert not only_e, (
        "{0}: the executable wrote {1} file(s) the API did not: {2}. The API is skipping part "
        "of the tool's work.".format(tag, len(only_e), only_e[:6]))
    assert not only_a, (
        "{0}: the API wrote {1} file(s) the executable did not: {2}".format(
            tag, len(only_a), only_a[:6]))
    assert not differing, (
        "{0}: {1} of {2} output files differ between the executable and the API: {3}. Same "
        "control file, same options, same seed - so the API is taking a different path "
        "through the algorithm.\n{4}".format(
            tag, len(differing), len(shared), differing[:6], "\n".join(detail)))
    return len(shared)


# ---- ies ------------------------------------------------------------------------------------

def api_vs_exe_ies_test():
    """ies: caller-driven loop must reproduce pestpp-ies exactly.

    Stressed with a lambda ladder and a subset - the subset path in particular evaluates only
    part of the ensemble per lambda, so an adapter that processed runs in the wrong order or
    lost the lambda sequence would show up here immediately.
    """
    def edit(pst):
        pst.control_data.noptmax = 2
        pst.pestpp_options["ies_num_reals"] = 8
        pst.pestpp_options["random_seed"] = 11
        pst.pestpp_options["ies_lambda_mults"] = "0.5,1.0,2.0"
        pst.pestpp_options["ies_subset_size"] = 4
        pst.pestpp_options["ies_initial_lambda"] = 10.0
        pst.pestpp_options["ies_bad_phi_sigma"] = 2.5

    d_exe, d_api = _pair("apiexe_ies", os.path.join("ies_10par_xsec", "template"),
                         "pest.pst", edit)
    _run_exe("pestpp-ies", d_exe, "pest.pst")
    with Ies.from_pst("pest.pst", workdir=d_api) as t:
        t.initialize()
        for _ in t.iterations():
            pass
        t.finalize()
    n = _assert_identical(d_exe, d_api, "ies")
    print("  ies: {0} output files identical".format(n))


# ---- mou ------------------------------------------------------------------------------------

def api_vs_exe_mou_test():
    """mou: one API generation must equal one executable generation.

    This is the test that caught MouAdapter re-implementing the loop body instead of calling
    it, which left MOEA's own generation counter at its initial value - so every generation
    wrote over the same .N. population files and the archive diverged.
    """
    def edit(pst):
        # g07 ships configured for sqp; drop those so the options in play are mou's
        for k in [k for k in pst.pestpp_options if k.startswith("sqp_")]:
            pst.pestpp_options.pop(k)
        pst.control_data.noptmax = 2
        pst.pestpp_options["mou_population_size"] = 10
        pst.pestpp_options["mou_generator"] = "de"
        pst.pestpp_options["mou_env_selector"] = "nsga"
        pst.pestpp_options["mou_de_f"] = 0.7
        pst.pestpp_options["mou_crossover_probability"] = 0.8
        pst.pestpp_options["random_seed"] = 11

    # g07 is checked in (constr_template is created dynamically by mou_tests.py and so is not
    # present in a fresh checkout), already carries obj_fn and l_constraint groups, and its
    # forward run is pure python - so this case needs no model binary at all
    d_exe, d_api = _pair("apiexe_mou", os.path.join("g07", "template"), "g07.pst", edit)
    _run_exe("pestpp-mou", d_exe, "g07.pst")
    with Mou.from_pst("g07.pst", workdir=d_api) as t:
        t.initialize()
        for _ in t.iterations():
            pass
        t.finalize()
    n = _assert_identical(d_exe, d_api, "mou")
    print("  mou: {0} output files identical".format(n))


# ---- sqp ------------------------------------------------------------------------------------

def api_vs_exe_sqp_test():
    """sqp: one API iteration must equal one executable iteration.

    This is the test that caught SqpAdapter calling solve_new_ensemble() and nothing else -
    the first statement of an sqp iteration, missing the gradient runs, the CMA and hessian
    updates, the constraint report, the pcs summary, and the advance of the tool's own
    iteration counter.
    """
    def edit(pst):
        obs = pst.observation_data
        obj, con = pst.nnz_obs_names[0], pst.nnz_obs_names[1]
        obs.loc[obj, "obgnme"] = "obj_fn"
        obs.loc[con, "obgnme"] = "l_head"          # 'l_' marks a less-than constraint
        obs.loc[con, "obsval"] = float(obs.loc[con, "obsval"]) * 1.5
        pst.pestpp_options["opt_obj_func"] = obj
        pst.pestpp_options["sqp_num_reals"] = 8
        pst.pestpp_options["random_seed"] = 11
        pst.control_data.noptmax = 2

    d_exe, d_api = _pair("apiexe_sqp", os.path.join("ies_10par_xsec", "template"),
                         "pest.pst", edit)
    _run_exe("pestpp-sqp", d_exe, "pest.pst")
    with Sqp.from_pst("pest.pst", workdir=d_api) as t:
        t.initialize()
        for _ in t.iterations():
            pass
        t.finalize()
    n = _assert_identical(d_exe, d_api, "sqp")
    print("  sqp: {0} output files identical".format(n))


# ---- da -------------------------------------------------------------------------------------

def api_vs_exe_da_test():
    """da: the API is NOT equivalent to pestpp-da, and this pins exactly how.

    Stated plainly because it would otherwise be assumed: driving da through the API today
    gives DIFFERENT NUMBERS from running pestpp-da on the same control file. Two reasons, both
    structural rather than a bug in the adapter:

      1. pestpp-da is a multi-cycle driver. It walks every assimilation cycle, rebuilding the
         scenario and the tool for each and carrying state forward. The API exposes a single
         cycle - the header says so - so none of the .global.* output exists.
      2. Even for one cycle the prior differs. The driver draws a GLOBAL ensemble spanning all
         cycles (generate_global_ensembles) and injects it with set_pe/set_oe/set_noise_oe,
         after forcing ies_include_base off. The API lets the cycle draw its own. Different
         draw, different ensemble, different everything downstream.

    So this test asserts the two things that are true and useful: the API produces a
    structurally valid single-cycle result with the same parameters and observations, and the
    executable additionally writes the .global.* family that marks the multi-cycle layer.

    It also asserts the values still DIFFER. That reads backwards, and it is deliberate: if
    someone closes this gap, this assertion fires and points at the docs and tests that then
    need updating - rather than the improvement landing silently and this test continuing to
    pass while describing a world that no longer exists.
    """
    def edit(pst):
        # every quantity marked cycle -1 ("all cycles"), which is the batch form - the only
        # shape the API can drive at all
        for df in (pst.parameter_data, pst.observation_data,
                   pst.model_input_data, pst.model_output_data):
            df.loc[:, "cycle"] = -1
        pst.control_data.noptmax = 1
        pst.pestpp_options["ies_num_reals"] = 8
        pst.pestpp_options["da_num_reals"] = 8
        pst.pestpp_options["random_seed"] = 11
        pst.pestpp_options["da_use_mda"] = True
        pst.pestpp_options["da_mda_init_fac"] = 5.0

    d_exe, d_api = _pair("apiexe_da", os.path.join("ies_10par_xsec", "template"),
                         "pest.pst", edit)
    _run_exe("pestpp-da", d_exe, "pest.pst")
    with Da.from_pst("pest.pst", workdir=d_api) as t:
        t.initialize()
        for _ in t.iterations():
            pass
        t.finalize()

    shared, differing, only_e, only_a = _compare_csvs(d_exe, d_api, "da")

    # (1) the single-cycle result is structurally what the executable produced for that cycle
    import pandas as pd
    for f in ("pest.0.1.par.csv", "pest.0.1.obs.csv"):
        assert f in shared, "da: expected both runs to write {0}, shared={1}".format(f, shared)
        a = pd.read_csv(os.path.join(d_exe, f), index_col=0)
        b = pd.read_csv(os.path.join(d_api, f), index_col=0)
        assert list(a.columns) == list(b.columns), \
            "da: {0} columns differ between exe and api - this is a real defect, not the " \
            "known prior-draw gap".format(f)
        assert a.shape == b.shape, \
            "da: {0} shape differs, {1} vs {2}".format(f, a.shape, b.shape)

    # (2) the exe writes the multi-cycle layer, the API does not
    global_only = [f for f in only_e if ".global." in f]
    assert global_only, (
        "da: the executable no longer writes any .global.* files. Either pestpp-da changed or "
        "this case stopped being multi-cycle; either way this test needs revisiting.")
    assert not only_a, "da: the API wrote files the executable did not: {0}".format(only_a)

    # (3) the values differ, for the reasons in the docstring
    assert differing, (
        "da: the API and pestpp-da now produce IDENTICAL output. That is good news and this "
        "test is now wrong. Replace it with _assert_identical(), and update the DA notes in "
        "pestpp-api.h and docs/api_part1/review_findings.md, which both say they diverge.")
    print("  da: single-cycle result structurally matches; {0} exe-only .global.* file(s), "
          "{1} file(s) differ in value (expected - see docstring)"
          .format(len(global_only), len(differing)))


if __name__ == "__main__":
    api_vs_exe_ies_test()
    api_vs_exe_mou_test()
    api_vs_exe_sqp_test()
    api_vs_exe_da_test()
    print("all api-vs-exe tests passed")

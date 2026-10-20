"""finite difference jacobian checks for pestpp-glm on an analytical model.

the model is three parameters and a handful of outputs with known
derivatives.  a central difference is exact for a quadratic and a one sided
difference is not, so the quadratic outputs tell the two apart with no
tolerance games: with a relative increment h on a parameter at 1.0,

    output       true d/dp   central      forward      backward
    lin  = p        1          1            1            1
    quad = p^2      2          2            2 + h        2 - h
    const = 1       0          0            0            0

so a run with forcen always_3 (two perturbation runs per parameter) must give
exactly 2 on the quadratic rows, and a run with always_2 must give 2 + h.
the parabolic method fits a parabola through the three points, which is also
exact for a quadratic, and must give zero on the constant row.

this was written for the report that pestpp-glm reads the first perturbation
run twice when building a central difference jacobian
(Jacobian_1to1::process_runs uses par_run.second[0] inside the loop over rid),
so every always_3 jacobian was one sided.  no model binaries needed."""
import os
import shutil
import platform
import numpy as np
import pandas as pd
import pyemu

bin_path = os.path.join("test_bin")
if "linux" in platform.platform().lower():
    bin_path = os.path.join(bin_path, "linux")
elif "darwin" in platform.platform().lower() or "macos" in platform.platform().lower():
    bin_path = os.path.join(bin_path, "mac")
else:
    bin_path = os.path.join(bin_path, "win")
bin_path = os.path.abspath("test_bin")
os.environ["PATH"] += os.pathsep + bin_path

# same convention as the other benchmark scripts: the exe path is relative to a
# run dir two levels below this one, and the repo dir is called pestpp
if os.path.exists(os.path.join("pestpp", "bin")):
    bin_path = os.path.join("..", "..", "pestpp", "bin")
else:
    bin_path = os.path.join("..", "..", "..", "..", "pestpp", "bin")
if "windows" in platform.platform().lower():
    exe_path = os.path.join(bin_path, "win", "pestpp-glm.exe")
elif "darwin" in platform.platform().lower() or "macos" in platform.platform().lower():
    exe_path = os.path.join(bin_path, "mac", "pestpp-glm")
else:
    exe_path = os.path.join(bin_path, "linux", "pestpp-glm")


def forward_run(ws="."):
    """the model: reads par.dat, writes obs.dat.  if a file named fail_below.txt
    is present it holds a parameter index and a threshold, and the model fails
    (no output, non zero exit) when that parameter is below the threshold - a
    stand in for a model that crashes on one side of a perturbation"""
    p = np.loadtxt(os.path.join(ws, "par.dat"))
    fb = os.path.join(ws, "fail_below.txt")
    if os.path.exists(fb):
        k, thresh = open(fb).read().split()
        if p[int(k)] < float(thresh):
            if os.path.exists(os.path.join(ws, "obs.dat")):
                os.remove(os.path.join(ws, "obs.dat"))
            raise SystemExit("deliberate failure for p{0} = {1}".format(k, p[int(k)]))
    out = {}
    for k in range(3):
        out["lin{0}".format(k)] = p[k]
        out["quad{0}".format(k)] = p[k] ** 2
    out["cross01"] = p[0] * p[1]
    out["const"] = 1.0
    with open(os.path.join(ws, "obs.dat"), "w") as f:
        for n, v in out.items():
            f.write("{0} {1:25.16E}\n".format(n, v))


def setup(template_d, derinc=0.02):
    """write the control file and model into template_d"""
    if os.path.exists(template_d):
        shutil.rmtree(template_d)
    os.makedirs(template_d)
    with open(os.path.join(template_d, "par.dat.tpl"), "w") as f:
        f.write("ptf ~\n")
        for k in range(3):
            f.write("~  p{0}  ~\n".format(k))
    np.savetxt(os.path.join(template_d, "par.dat"), np.ones(3))
    forward_run(template_d)
    onames = [l.split()[0] for l in open(os.path.join(template_d, "obs.dat"))]
    with open(os.path.join(template_d, "obs.dat.ins"), "w") as f:
        f.write("pif ~\n")
        for n in onames:
            f.write("l1 w !{0}!\n".format(n))
    shutil.copy2(os.path.abspath(__file__), os.path.join(template_d, "glm_jacobian_tests.py"))
    with open(os.path.join(template_d, "forward_run.py"), "w") as f:
        f.write("import glm_jacobian_tests\nif __name__ == '__main__':\n    glm_jacobian_tests.forward_run()\n")
    pst = pyemu.Pst.from_io_files([os.path.join(template_d, "par.dat.tpl")], [os.path.join(template_d, "par.dat")],
                                  [os.path.join(template_d, "obs.dat.ins")], [os.path.join(template_d, "obs.dat")],
                                  pst_path=".")
    par = pst.parameter_data
    par.loc[:, "partrans"] = "none"
    par.loc[:, "parchglim"] = "relative"
    par.loc[:, "parval1"] = 1.0
    par.loc[:, "parlbnd"] = 0.1
    par.loc[:, "parubnd"] = 10.0
    par.loc[:, "pargp"] = "p"
    pst.rectify_pgroups()
    pg = pst.parameter_groups
    pg.loc["p", "inctyp"] = "relative"
    pg.loc["p", "derinc"] = derinc
    pg.loc["p", "derinclb"] = 0.0
    pg.loc["p", "derincmul"] = 1.0
    pst.observation_data.loc[:, "weight"] = 1.0
    pst.observation_data.loc[:, "obgnme"] = "obs"
    pst.model_command = "python forward_run.py"
    pst.control_data.noptmax = -1
    pst.write(os.path.join(template_d, "jac.pst"), version=2)
    return pst


def run_jacobian(template_d, m_d, forcen, dermthd, derinc=0.02, fail_below=None):
    """one pestpp-glm jacobian run (noptmax -1), returns the jco as a dataframe.
    fail_below = (par index, threshold) makes the model fail on that side"""
    if os.path.exists(m_d):
        shutil.rmtree(m_d)
    shutil.copytree(template_d, m_d)
    if fail_below is not None:
        with open(os.path.join(m_d, "fail_below.txt"), "w") as f:
            f.write("{0} {1}\n".format(*fail_below))
    pst = pyemu.Pst(os.path.join(m_d, "jac.pst"))
    pst.parameter_groups.loc["p", "forcen"] = forcen
    pst.parameter_groups.loc["p", "dermthd"] = dermthd
    pst.parameter_groups.loc["p", "derinc"] = derinc
    pst.write(os.path.join(m_d, "jac.pst"), version=2)
    pyemu.os_utils.run("{0} jac.pst".format(exe_path), cwd=m_d)
    jco = pyemu.Jco.from_binary(os.path.join(m_d, "jac.jcb")).to_dataframe()
    runs = None
    for line in open(os.path.join(m_d, "jac.rec")):
        if "total model runs" in line.lower() or "model runs completed" in line.lower():
            runs = line.strip()
    return jco, runs


def expected(forcen, derinc=0.02, p0=1.0):
    """the exact derivatives each scheme should produce at p0"""
    h = derinc * p0
    quad = {"always_2": 2.0 * p0 + h, "always_3": 2.0 * p0}[forcen]
    return {"lin": 1.0, "quad": quad, "cross": 1.0, "const": 0.0}


def glm_central_difference_test(derinc=0.02, tol=1.0e-8):
    """always_3 must be a central difference: exact on the quadratic rows.
    always_2 is the control, and shows the one sided value the bug produces"""
    model_d = "glm_jacobian"
    template_d = os.path.join(model_d, "template")
    os.makedirs(model_d, exist_ok=True)
    setup(template_d, derinc)
    results = {}
    for forcen, dermthd in (("always_2", "outside_pts"), ("always_3", "outside_pts"),
                            ("always_3", "parabolic"), ("always_3", "best_fit")):
        tag = forcen + "_" + dermthd
        jco, runs = run_jacobian(template_d, os.path.join(model_d, "master_" + tag), forcen, dermthd, derinc)
        results[tag] = jco
        print("\n{0}  ({1})".format(tag, runs))
        print(jco.loc[["lin0", "quad0", "cross01", "const"], ["p0", "p1"]].to_string(float_format=lambda v: "{0:.6f}".format(v)))

    def check(tag, forcen):
        jco = results[tag]
        want = expected(forcen, derinc)
        errs = []
        for k in range(3):
            pk = "p{0}".format(k)
            for oname, key in (("lin{0}".format(k), "lin"), ("quad{0}".format(k), "quad"), ("const", "const")):
                got = jco.loc[oname, pk]
                if abs(got - want[key]) > tol:
                    errs.append("{0}: d {1}/d {2} = {3:.6f}, want {4:.6f}".format(tag, oname, pk, got, want[key]))
            # a parameter the output does not depend on must have a zero column entry
            for other in range(3):
                if other == k:
                    continue
                got = jco.loc["quad{0}".format(k), "p{0}".format(other)]
                if abs(got) > tol:
                    errs.append("{0}: d quad{1}/d p{2} = {3:.6f}, want 0".format(tag, k, other, got))
        got = jco.loc["cross01", "p0"]
        if abs(got - want["cross"]) > tol:
            errs.append("{0}: d cross01/d p0 = {1:.6f}, want {2:.6f}".format(tag, got, want["cross"]))
        return errs

    errs = check("always_2_outside_pts", "always_2")
    assert not errs, "the one sided control is wrong, the test itself is broken:\n" + "\n".join(errs)
    failures = []
    for tag in ("always_3_outside_pts", "always_3_parabolic", "always_3_best_fit"):
        failures += check(tag, "always_3")
    assert not failures, "always_3 is not a central difference:\n" + "\n".join(failures)
    return results


def glm_failed_perturbation_test(derinc=0.02, tol=1.0e-8):
    """always_3 where the model fails on the negative perturbation of p1 only.
    the surviving run must give a one sided (forward) quotient for p1, 2 + h on
    its quadratic row, while p0 and p2 stay central, and the run must not fail
    outright or drop p1.  before the fix the failed run was never read, so
    this case looked like a success with a backward quotient"""
    model_d = "glm_jacobian"
    template_d = os.path.join(model_d, "template")
    os.makedirs(model_d, exist_ok=True)
    setup(template_d, derinc)
    m_d = os.path.join(model_d, "master_always_3_fail_p1")
    jco, runs = run_jacobian(template_d, m_d, "always_3", "outside_pts", derinc, fail_below=(1, 1.0))
    print("\nalways_3 with the negative perturbation of p1 failing")
    print(jco.loc[["quad0", "quad1", "quad2", "lin1", "const"], ["p0", "p1", "p2"]].to_string(
        float_format=lambda v: "{0:.6f}".format(v)))
    h = derinc
    errs = []
    if abs(jco.loc["quad0", "p0"] - 2.0) > tol or abs(jco.loc["quad2", "p2"] - 2.0) > tol:
        errs.append("p0 or p2 lost its central difference")
    if abs(jco.loc["quad1", "p1"] - (2.0 + h)) > tol:
        errs.append("p1 should be the forward quotient {0:.4f} from the surviving run, got {1:.6f}".format(
            2.0 + h, jco.loc["quad1", "p1"]))
    if abs(jco.loc["lin1", "p1"] - 1.0) > tol:
        errs.append("lin1/p1 should be 1, got {0:.6f}".format(jco.loc["lin1", "p1"]))
    rec = open(os.path.join(m_d, "jac.rec")).read().lower()
    note = ("one sided" in rec) or ("one-sided" in rec) or ("degrad" in rec)
    print("rec file notes the degraded column:", note)
    assert not errs, "\n".join(errs)
    return jco, note


if __name__ == "__main__":
    glm_central_difference_test()
    glm_failed_perturbation_test()

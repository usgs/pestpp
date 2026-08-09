"""Keep the tests off the tracked fixtures.

The benchmarks used to read AND write the same `<model>/template` directory that is tracked in
git: a test would load `template/pest.pst`, change options, and write it straight back before
starting workers. That works once. The next run then reads a control file the previous run
edited, so it is no longer testing the fixture - it is testing whatever the last run left
behind.

Three separate false failures in one day came from exactly this, and every one of them cost
real diagnosis time because the symptom appears somewhere else entirely:

  * mf6_freyberg/template          - pst.write() reformatting showing up as a dirty tree
  * secondary_marker_test/template - same
  * glm_10par_xsec/template        - tenpar_superpar_restart_test failing on a FOSM
                                     comparison, because the prior it built came from a
                                     mutated pst. From a clean fixture the diff is exactly 0.

CI never sees any of it, because every job clones fresh. So the failure mode is specifically
"passes on CI, fails locally, for reasons unrelated to the change under test" - the most
expensive kind.

scratch_template() hands back a private copy instead.
"""
import os
import shutil


def scratch_template(template_d, suffix="_work"):
    """A private, writable copy of a tracked fixture directory. Returns the copy's path.

    The copy is a SIBLING of the original, at the same depth, and that is not cosmetic: the
    benchmarks reach the built binaries through relative paths like
    ../../../../pestpp/bin/<plat>/pestpp-ies, resolved from wherever the case directory sits.
    Put the copy anywhere else - a temp dir, a subdirectory - and every one of those breaks.

    Recreated from scratch on each call, so a test never inherits the previous test's leavings
    either.
    """
    if not os.path.isdir(template_d):
        raise Exception("fixture template not found: {0}".format(template_d))
    work_d = template_d + suffix
    if os.path.exists(work_d):
        shutil.rmtree(work_d)
    shutil.copytree(template_d, work_d)
    return work_d

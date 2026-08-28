"""Launchers for the bundled pest++ command line tools.

The wheel carries the real executables in _binaries/ next to the shared library. They cannot be
console_scripts directly - an entry point has to name a python callable, and these are native
binaries - so pip installs a small stub for each one and the stub hands over to the binary here.

Two details that are easy to get wrong:

  - the executable bit does not survive the wheel. Files installed as package data come out
    mode 644, because a wheel's RECORD does not carry permissions. So it is set here, once, on
    first use. Only the scripts/ directory of a wheel gets +x from pip, and these are not there.

  - posix replaces the process with execv, windows does not have that, so it spawns and waits.
    Replacing matters more than it looks: a PANTHER agent that gets ctrl+c should receive the
    signal itself rather than have a python parent intercept it, and a wrapper process would
    also confuse anything watching the pid.
"""
import os
import stat
import subprocess
import sys
from pathlib import Path

BIN_DIR = Path(__file__).resolve().parent / "_binaries"


def _binary(name):
    exe = BIN_DIR / (name + ".exe" if os.name == "nt" else name)
    if not exe.exists():
        raise SystemExit(
            "{0} is not in this install of pestpp.\n"
            "The wheel should carry it at {1}. If you installed from source rather than a\n"
            "wheel, the command line tools are not part of that - build them with cmake, or\n"
            "install the wheel from pypi.".format(name, exe))
    if os.name != "nt" and not os.access(str(exe), os.X_OK):
        # +x for owner/group/other, mirroring what pip does for real scripts
        mode = exe.stat().st_mode
        try:
            exe.chmod(mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        except OSError as e:
            raise SystemExit("cannot make {0} executable: {1}".format(exe, e))
    return str(exe)


def _run(name):
    exe = _binary(name)
    argv = [exe] + sys.argv[1:]
    if os.name == "nt":
        # no execv worth having on windows - it returns instead of replacing, which would exit
        # the stub while the tool is still running
        sys.exit(subprocess.call(argv))
    os.execv(exe, argv)


# one entry point per tool. pestpp-selftest is deliberately absent: it is a build gate, not
# something a user of the package has any reason to run.
def glm():  _run("pestpp-glm")
def ies():  _run("pestpp-ies")
def opt():  _run("pestpp-opt")
def sen():  _run("pestpp-sen")
def swp():  _run("pestpp-swp")
def mou():  _run("pestpp-mou")
def sqp():  _run("pestpp-sqp")
def da():   _run("pestpp-da")

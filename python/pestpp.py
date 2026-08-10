"""Ergonomic layer over the PEST++ shared library.

The thin binding is ``pestpp_lib``: one method per C symbol, same names, no cleverness. This
module is the opinionated layer on top of it -- the same split pypestutils uses
(``pestutilslib.py`` / ``helpers.py``) and MODFLOW 6 uses (``xmipy`` / ``modflowapi``).

Naming follows ``docs/api_part1/python_api_conventions.md``: properties with no prefix for
cheap derived state, ``get_*`` only when it takes arguments, ``from_*`` to construct, plain
verbs to act, ``*_df`` for DataFrames, ``*_view`` context managers for borrowed memory, and
module-level ``run_*`` functions for the whole-job case.

Quick start::

    from pestpp import Ies

    with Ies.from_pst("pest.pst", workdir="master") as ies:
        ies.initialize()
        for step in ies.iterations():
            print(step.iter, step.phi_mean)
        df = ies.par_df()
"""

from __future__ import annotations

import os
import platform
import shutil
import subprocess
import sys
import time
import weakref
from contextlib import contextmanager, nullcontext
from dataclasses import dataclass
from typing import TypeVar

import numpy as np
import pandas as pd

# pyemu is a REQUIRED dependency of this module, deliberately.
#
# The thin binding (pestpp_lib) stays free of it - ctypes and numpy only - so anyone writing a
# binding in another language has a reference implementation with no python ecosystem in it.
# This layer is the opposite: its whole job is to be comfortable for someone who already
# thinks in pyemu, which means handing back real Pst and Ensemble objects rather than frames
# that merely resemble them. Half-integrating would be worse than not integrating - a
# DataFrame that looks like a ParameterEnsemble but cannot enforce() is a trap.
try:
    import pyemu
except ImportError as e:                                     # pragma: no cover
    raise ImportError(
        "pestpp.py requires pyemu (the thin layer, pestpp_lib.py, does not). "
        "Install it with `pip install pyemu`.") from e

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pestpp_lib import (  # noqa: E402
    PestppLib, PestppError, _UNSET,
    PAR_EN, OBS_EN, NOISE_EN, WEIGHTS_EN,
    STACK_PAR_EN, STACK_OBS_EN, NESTED_PAR_EN, MEMBER_STACK_EN,
    TOOL_IES, TOOL_DA, TOOL_MOU, TOOL_SQP, TOOL_GLM,
    RM_SERIAL, RM_PANTHER, RM_EXTERNAL,
    PHI_MEAS, PHI_COMPOSITE, PHI_REGUL, PHI_ACTUAL, PHI_NOISE,
    WORKER_COMPLETED, WORKER_FAILED, WORKER_TIMED_OUT,
)

from pestpp_progress import Progress, auto as _progress_auto  # noqa: E402

# Bound to _Tool so from_pst() reports the SUBCLASS it was called on: Ies.from_pst() infers
# as Ies, Mou.from_pst() as Mou. Without it a classmethod returning cls() is opaque to jedi
# and pylance, which is what makes `Ies.from_pst("pest.pst").<TAB>` come back empty in a
# notebook - the single most useful completion there is.
_ToolT = TypeVar("_ToolT", bound="_Tool")

__all__ = [
    "Ies", "Da", "Mou", "Sqp", "Glm", "IterationStep", "Candidate", "PestppError", "ExpiredViewError",
    "Progress", "run_ies", "run_da", "run_mou", "run_sqp", "find_library",
    "PHI_MEAS", "PHI_COMPOSITE", "PHI_REGUL", "PHI_ACTUAL", "PHI_NOISE",
]


# ---- locating the shared library ----------------------------------------------------------

def find_library(lib_path: str | None = None) -> str:
    """Locate the built C ABI shared library.

    Searched rather than hardcoded because the path depends on the build generator. Set
    ``PESTPP_LIB`` to skip the search entirely.
    """
    if lib_path is not None:
        if not os.path.exists(lib_path):
            raise FileNotFoundError(lib_path)
        return lib_path
    env = os.environ.get("PESTPP_LIB")
    if env:
        if not os.path.exists(env):
            raise FileNotFoundError("PESTPP_LIB points at {0}, which does not exist".format(env))
        return env

    plat = platform.platform().lower()
    # pestpp-api.<so|dll|dylib> - no "lib" prefix, matching the pestpp-* executables
    name = ("pestpp-api.dll" if ("windows" in plat or os.name == "nt")
            else "pestpp-api.dylib" if "darwin" in plat or "macos" in plat
            else "pestpp-api.so")
    here = os.path.dirname(os.path.abspath(__file__))
    roots = [os.path.join(os.path.dirname(here), "build"), os.path.join(here, "..", "..", "build")]

    # An installed copy next to the executables COMPETES on mtime; it does not win outright.
    # Short-circuiting to it meant a stale install silently shadowed a fresh build - the
    # symptom is an AttributeError from ctypes about a symbol that plainly exists in the
    # source, which sends you looking in exactly the wrong place.
    found = []
    for cand in (os.path.join(os.path.dirname(here), "bin", name),):
        if os.path.exists(cand):
            found.append(os.path.abspath(cand))
    for root in roots:
        if not os.path.isdir(root):
            continue
        for dirpath, dirnames, filenames in os.walk(root):
            # Never a packaging staging area. CPack leaves a full copy of the install tree
            # under _CPack_Packages, so a stale one from an earlier build sits there looking
            # exactly like the real thing - and being picked would mean silently testing a
            # library that is not the one just compiled.
            dirnames[:] = [d for d in dirnames if d != "_CPack_Packages"]
            if name in filenames:
                found.append(os.path.abspath(os.path.join(dirpath, name)))
    if found:
        # newest wins, so a fresh build beats anything left over
        return max(found, key=os.path.getmtime)
    raise FileNotFoundError(
        "could not find {0} under {1}. Build pest++ first, or set PESTPP_LIB.".format(name, roots))


# ---- name case ----------------------------------------------------------------------------
#
# pest++ upper-cases every name internally; pyemu keeps them lower. That mismatch is a
# constant papercut when the two are used together - a name lifted from a Pst does not match
# one read back from the library. The helper layer absorbs it: anything you pass in is
# accepted in any case, and anything handed back can be lower-cased on request.

def _up(name) -> str:
    """A single name in the case the library uses."""
    return str(name).upper()


def _up_all(names):
    return [_up(n) for n in names]


def _named(df, index_name=None, columns_name=None):
    """Give a frame the axis names pyemu uses, so merges and reset_index just work.

    pyemu names its axes (`parnme`, `obsnme`) and pest++'s own ensemble csv header is
    `real_name`; pyemu's result handler calls the ensemble index `realization`. Matching that
    is free and it is the difference between `df.reset_index()` producing a column a pyemu
    user recognises and producing one called "index".
    """
    if index_name is not None:
        df.index.name = index_name
    if columns_name is not None:
        df.columns.name = columns_name
    return df


def _maybe_lower(df, lower: bool, axis="both"):
    """Lower-case a frame's labels when asked, for joining against pyemu objects."""
    if not lower:
        return df
    if isinstance(df, pd.Series):
        return df.rename(index=str.lower)
    if axis in ("both", "index"):
        df = df.rename(index=str.lower)
    if axis in ("both", "columns"):
        df = df.rename(columns=str.lower)
    return df


# ---- keeping the library's console output out of the way ----------------------------------

@contextmanager
def _capture_output(lib, path):
    """Send the library's console output to a file while this block runs.

    The redirect is performed by the LIBRARY, not here. Python's own
    ``contextlib.redirect_stdout`` rebinds ``sys.stdout`` and does nothing to C++ ``cout``;
    even redirecting file descriptor 1 from python is not enough, because on windows the
    library links the static CRT and so has its own descriptor table -- moving python's
    descriptor 1 relocates the process std handle (so child processes like the model do land
    in the file) while the library's own descriptor 1 still points at the console. Asking the
    library to redirect its own descriptor is the only thing that captures both.

    Captured to a file rather than discarded: the .rec file does not carry everything the
    console does, and swallowing a run manager's complaints outright would be worse than the
    noise.
    """
    sys.stdout.flush()
    token = lib.redirect_output(path)
    try:
        yield
    finally:
        # try/finally makes this strictly LIFO, which is what the library now requires: fd 1
        # is process-global, so restoring out of order would leave stdout pointing at another
        # session's log permanently. Nested `with` blocks unwind correctly by construction.
        lib.restore_output(token)


# ---- what an iteration reports ------------------------------------------------------------

@dataclass
class IterationStep:
    """One step of the algorithm, as yielded by :meth:`_Tool.iterations`."""
    iter: int
    n_reals: int
    retried: bool          # the step was rejected and the algorithm wants another attempt
    #: runs still needed before this iteration can complete. Always 0 from :meth:`_Tool.solve`;
    #: non-zero only from :meth:`_Tool.finish_solve` with ``defer_runs=True``, where it means
    #: the remaining realizations are waiting to be queued and run.
    pending_runs: int = 0
    phi_mean: float = float("nan")
    phi_std: float = float("nan")
    phi_min: float = float("nan")
    phi_max: float = float("nan")

    def __repr__(self):
        return ("IterationStep(iter={0}, n_reals={1}, phi_mean={2:.6g}{3})"
                .format(self.iter, self.n_reals, self.phi_mean, ", RETRY" if self.retried else ""))


# ---- borrowed views -------------------------------------------------------------------------

class ExpiredViewError(PestppError):
    """Raised when a zero-copy view is used after its ``with`` block, or after a resize."""


class EnsembleViewProxy:
    """An ndarray you may only use while the view is live.

    The array itself is the tool's Eigen buffer. Two things can make it stop being yours:
    leaving the ``with`` block, and the ensemble reallocating underneath you (a resize, a
    membership change, an algorithm step). Neither is visible in the array - it keeps its
    shape and returns numbers either way, and those numbers are read from freed memory.

    So this stands in front of it. Every operation asks two questions first - has the block
    ended, and does the library still recognise the buffer - and raises
    :class:`ExpiredViewError` rather than reading. That is a real guard, not a convention:
    the array is never handed out, so a caller cannot keep a reference past the block by
    assigning it elsewhere.

    Work in BULK. Each operation costs a call into the library, so ``a[:, 3] = x`` is one
    check and a python loop over elements is one per element. That is the natural numpy
    style anyway.

    The one deliberate escape is ``np.asarray(proxy)`` (and ``.copy()``, which is the safe
    version): it hands back the underlying array while the view is live, because that is what
    makes the proxy work with the rest of numpy. Copy it if you want to keep it.
    """

    __slots__ = ("_arr", "_lib", "_token", "_live")

    def __init__(self, arr, lib, token):
        object.__setattr__(self, "_arr", arr)
        object.__setattr__(self, "_lib", lib)
        object.__setattr__(self, "_token", token)
        object.__setattr__(self, "_live", True)

    # -- the guard --
    def _live_array(self):
        if not self._live:
            raise ExpiredViewError(
                "this ensemble view has expired: it is only valid inside the 'with' block "
                "that produced it. Use .copy() inside the block to keep the values, or "
                "par_df()/obs_df() for a labelled copy.")
        if not self._lib.view_is_valid(self._token):
            raise ExpiredViewError(
                "this ensemble view is no longer valid: the ensemble's storage was replaced "
                "or resized (a membership change or an algorithm step will do it). Take a "
                "fresh view.")
        return self._arr

    def _expire(self):
        object.__setattr__(self, "_live", False)
        try:
            self._lib.release_view(self._token)
        except Exception:
            # the handle may already be gone; expiry still stands
            pass

    @property
    def valid(self) -> bool:
        """True while this view may still be used. Never raises."""
        if not self._live:
            return False
        try:
            return self._lib.view_is_valid(self._token)
        except Exception:
            return False

    # -- ndarray surface --
    def __array__(self, dtype=None, copy=None):
        arr = self._live_array()
        if copy is False:
            return arr if dtype is None else arr.astype(dtype, copy=False)
        if dtype is None and not copy:
            return arr
        return arr.astype(dtype) if dtype is not None else arr.copy()

    def __getitem__(self, key):
        return self._live_array()[key]

    def __setitem__(self, key, value):
        self._live_array()[key] = value

    def __len__(self):
        return len(self._live_array())

    def __iter__(self):
        return iter(self._live_array())

    def __repr__(self):
        if not self._live:
            return "<EnsembleViewProxy (expired)>"
        return "<EnsembleViewProxy shape={0}>".format(self._arr.shape)

    def copy(self):
        """A normal ndarray you own. The way to keep values past the block."""
        return self._live_array().copy()

    def __getattr__(self, name):
        # shape, dtype, ndim, size, T, mean, sum, ... all forward, guarded
        if name.startswith("_"):
            raise AttributeError(name)
        return getattr(self._live_array(), name)


# ---- the tools ----------------------------------------------------------------------------

class _Tool:
    """Shared behaviour. Use :class:`Ies`, :class:`Da`, :class:`Mou` or :class:`Sqp`.

    **Name case.** pest++ upper-cases names internally, pyemu keeps them lower. Every method
    here that takes names accepts either, so a list lifted straight from a ``Pst`` works.
    Frames come back in the library's case by default -- which is what the tool's own .rec and
    .csv output use -- and every ``*_df`` takes ``lower=True`` to hand back pyemu-cased labels
    for joining the other way.
    """

    _tool_id = None
    _has_phi = True

    def __init__(self, lib: PestppLib, workdir: str, workers=None, quiet: bool = True,
                 pst_file: str = "pest.pst"):
        self._lib = lib
        self._workdir = workdir
        self._workers = workers or []
        self._initialized = False
        self._queued = 0
        self._quiet = quiet
        self._log = os.path.join(workdir, "pestpp.stdout.log")
        self._pst_file = pst_file
        self._pst = None            # lazily parsed by the .pst property

        # Cleanup backstop - see the note in PestppLib. It matters more here, because this
        # layer owns SUBPROCESSES as well as the handle: the panther agents started by
        # create(workers=N). `with` cannot cover the case this is really for, which is a
        # notebook re-running its setup cell -
        #
        #     ies = Ies.create(wd, workers=8)     # run this three times
        #
        # rebinding the name drops the old session with no cleanup, and there are now 24 agents
        # alive and three sessions holding ports, silently. Dropping the reference is what
        # triggers the finalizer, so that case is exactly the one it catches.
        self._finalizer = weakref.finalize(self, _release_session, self._workers, self._lib)

    def _q(self):
        """Capture the library's console output for the duration of a call.

        """
        return _capture_output(self._lib, self._log) if self._quiet else nullcontext()

    @property
    def log_file(self) -> str:
        """Where the library's console output went, when quiet."""
        return self._log

    # -- construction ----------------------------------------------------------------------

    @classmethod
    def from_pst(cls: type[_ToolT], pst_file, workdir: str = ".", workers: int = 0,
                 port: int | str = 4004, lib_path: str | None = None,
                 worker_root: str | None = None, exe_path: str | None = None,
                 run_manager: int | None = None, quiet: bool = True, **options) -> _ToolT:
        """Open a session on a control file.

        ``workers`` greater than zero starts a PANTHER master and that many worker processes,
        each in its own copy of ``workdir`` under ``worker_root``. Leave it at zero for the
        serial run manager, which is simpler and right for small problems.

        Extra keyword arguments are set as pest++ options before initialization, in their
        native python types, so ``Ies.from_pst("pest.pst", ies_num_reals=50,
        ies_lambda_mults=[0.1, 1.0, 10.0], ies_no_noise=True)`` does what it looks like.
        See :meth:`set_option` for how each type is spelled.

        ``quiet`` (the default) captures the library's console output to
        ``<workdir>/pestpp.stdout.log`` instead of letting it flood the session. Pass
        ``quiet=False`` to watch it live.

        ``pst_file`` may be a path OR a :class:`pyemu.Pst`. Passing the object matches
        pyemu's own ``from_*`` convention and is the natural end of a PstFrom workflow -
        build it, hand it over, never write it out yourself. It IS written to ``workdir``,
        because pest++ reads a file and the workers need one too; ``pst.filename`` names it
        if it has one, otherwise it becomes ``pest.pst``.
        """
        workdir = os.path.abspath(workdir)
        if hasattr(pst_file, "write") and hasattr(pst_file, "parameter_data"):
            pst_obj = pst_file
            name = os.path.basename(getattr(pst_obj, "filename", None) or "pest.pst")
            if not name.lower().endswith(".pst"):
                name += ".pst"
            os.makedirs(workdir, exist_ok=True)
            pst_obj.write(os.path.join(workdir, name), version=2)
            pst_file = name
        parallel = workers > 0
        if run_manager is None:
            run_manager = RM_PANTHER if parallel else RM_SERIAL
        # the library has to exist before it can redirect itself, so create is not captured;
        # everything after it is
        lib = PestppLib(find_library(lib_path), cls._tool_id, pst_file, workdir,
                        port=str(port) if run_manager == RM_PANTHER else None,
                        run_manager=run_manager)
        for key, value in options.items():
            lib.set_option(key, value)

        procs = []
        if parallel:
            # started before any blocking call: the master accepts workers whenever, but this
            # process is single-threaded, so anything started after initialize() would never
            # be reached
            procs = _start_workers(workdir, workers, port, worker_root, exe_path,
                                   cls._agent_exe(), pst_file)
        return cls(lib, workdir, procs, quiet=quiet, pst_file=pst_file)

    @staticmethod
    def _agent_exe():
        return "pestpp-ies"

    # -- lifecycle -------------------------------------------------------------------------

    def initialize(self, defer_runs: bool = False) -> int:
        """Set the tool up, evaluating the prior ensemble.

        ``defer_runs=True`` stops short of running it and returns how many runs are waiting,
        so the caller can drive them (and, more usefully, replace the drawn ensemble first).
        Follow with :meth:`queue_runs`, :meth:`run`, :meth:`process_runs`, then
        :meth:`finish_initialize`.
        A return of 0 means there was nothing to hand over.
        """
        with self._q():
            if defer_runs:
                return self._lib.initialize_prepare()
            self._lib.initialize()
        self._initialized = True
        return 0

    def finish_initialize(self) -> None:
        """Process the prior-ensemble results after a deferred :meth:`initialize`."""
        with self._q():
            self._lib.initialize_finish()
        self._initialized = True

    def solve(self, defer_runs: bool = False) -> IterationStep | int:
        """One iteration (a generation, for mou).

        ``defer_runs=True`` generates the upgrade candidates and stops WITHOUT running them,
        returning how many runs they imply, so the candidates can be inspected or changed
        before they are evaluated::

            n = ies.solve(defer_runs=True)
            for c in ies.candidates():
                c.par_df()                       # or edit through c.par_view()
            ies.queue_runs()                     # or queue_runs(reals=[...])
            ies.run(); ies.process_runs()
            step = ies.finish_solve(defer_runs=True)
            while step.pending_runs:
                ies.queue_runs(); ies.run(); ies.process_runs()
                step = ies.finish_solve(defer_runs=True)

        A return of 0 means the iteration finished during preparation - a lambda that could
        not be generated, or the non-iterative shortcut - and there is nothing to run or
        finish.

        ies and mou only. da's one iteration is a whole noptmax loop over a cycle and sqp's
        line search issues several run batches per iteration, so neither can be split this way
        and both raise.
        """
        from pestpp_lib import PESTPP_RETRY
        if defer_runs:
            with self._q():
                return self._lib.solve_prepare()
        with self._q():
            status = self._lib.solve_iteration()
        return self._step(retried=(status == PESTPP_RETRY))

    def finish_solve(self, defer_runs: bool = False) -> IterationStep:
        """Continue a deferred :meth:`solve` after its runs have been processed.

        ``defer_runs=True`` stops rather than running the remaining realizations itself: after
        a subset of the ensemble has picked the winning lambda, the REST of the ensemble still
        has to be run at it. The count comes back as ``step.pending_runs``; queue, run and
        process them, then call this again. ``pending_runs`` of 0 means the iteration is done.

        ``defer_runs=False`` (the default) runs them internally and completes the iteration in
        one call.
        """
        from pestpp_lib import PESTPP_RETRY
        with self._q():
            status, pending = self._lib.solve_finish(defer_runs)
        return self._step(retried=(status == PESTPP_RETRY), pending_runs=pending)

    def candidates(self) -> list:
        """The upgrade candidates of an open deferred :meth:`solve`.

        ies generates one per lambda x scale-factor combination; mou generates one. Each is a
        :class:`Candidate` exposing the same frame and view calls the tool's own ensembles do,
        so editing one is the same code as editing :meth:`par_df` / :meth:`par_view`.
        """
        return [Candidate(self, i) for i in range(self._lib.get_candidate_count())]

    def run_manager_settings(self, overdue: bool = True) -> dict:
        """What the LIVE run manager is using - not what :meth:`get_option` reports.

        The run manager takes its tuning values when the session is created, so these two can
        disagree. :meth:`set_option` now pushes the four below onto the running manager, and
        this is how you confirm it: ``get_option`` answers "what did I ask for", this answers
        "what will the master actually do with a run that is running long".

        ``max_run_fail``, ``overdue_resched_fac``, ``overdue_giveup_fac`` and
        ``overdue_giveup_minutes`` can all be changed mid-run - the scheduling loop re-reads
        them every pass, so a change reaches the runs already in flight. The remaining panther
        options are consumed when the manager is built and are refused once it is running.

        ``overdue`` must be False on a serial or external session, which has no overdue policy.
        """
        return self._lib.get_run_manager_settings(overdue=overdue)

    def set_option(self, key: str, value) -> None:
        """Set a ++ option or a * control data value. An unknown key raises.

        Values are given in their NATIVE python form and formatted here - bools as
        ``true``/``false``, ints and floats as themselves, and lists, tuples or arrays as the
        comma-separated form the vector options are parsed from::

            ies.set_option("ies_num_reals", 50)
            ies.set_option("ies_no_noise", True)
            ies.set_option("ies_lambda_mults", [0.1, 1.0, 10.0])
            ies.set_option("ies_n_iter_reinflate", [3, 999])

        A str is passed through untouched, so anything already in control-file spelling still
        works. A type with no sensible spelling raises rather than being stringified into
        something the parser rejects later under the option's name.

        Some options are consumed once during initialize() and cannot change the current run
        afterwards; setting one of those late is accepted but has no effect.
        """
        self._lib.set_option(key, value)

    def get_option(self, key: str, default=_UNSET):
        """One option's value, as a string. An unknown key raises unless a default is given.

        The default is what separates the two questions a caller might be asking. An option
        set to the empty string returns ``""``; an option that does not exist returns the
        default, or raises if none was supplied.
        """
        return (self._lib.get_option(key) if default is _UNSET
                else self._lib.get_option(key, default))

    def has_option(self, key: str) -> bool:
        """Whether this build knows the option at all, whatever its value.

        The way to feature-detect against a library you did not build -- ``set_option`` on an
        unknown key raises, so probing with it is not an option.
        """
        return self._lib.has_option(key)

    @property
    def noptmax(self) -> int:
        """Maximum iterations from the control data.

        Negative values are the usual pest special cases (-1 prior evaluation only, -2 a
        single base-realization run) and mean ZERO solution iterations -- the tools loop
        ``for (i = 0; i < noptmax; i++)``, which does not execute for a negative bound.
        """
        # deliberately not swallowed: a value we cannot read must not silently become a
        # loop bound, which is how "no iterations" turns into "iterate forever"
        return int(self._lib.get_option("NOPTMAX"))

    def iterations(self, max_iter: int | None = None, progress=False):
        """Yield an :class:`IterationStep` per iteration.

        Stops at ``noptmax``, or when the tool's own convergence test fires, whichever comes
        first. Both matter: should_terminate() only reports the phi-based criteria, so a loop
        driven by it alone would ignore noptmax entirely and run until phi stopped improving -
        which is not what a control file saying noptmax=3 asked for.

        A NEGATIVE noptmax yields nothing at all, matching the tools: they loop
        ``for (i = 0; i < noptmax; i++)``, which does not execute. noptmax=-1 (evaluate the
        prior and stop) is the standard way to run prior monte carlo with pestpp-ies, so
        treating negative as "unlimited" would turn the most common invocation into an
        endless loop.

        ``max_iter`` overrides noptmax, for a caller that wants a few steps and a look around.

        ``progress=True`` draws a bar across the iterations, annotated with phi as it moves.
        Each iteration is one blocking call into the library, so this advances per ITERATION,
        not per run - use :meth:`run` with ``progress=True``, or a deferred :meth:`solve`, to
        watch the model runs inside one.
        """
        limit = max_iter if max_iter is not None else self.noptmax
        if limit < 0:
            limit = 0
        n = 0
        with _progress_auto(progress) as bar:
            bar.start("iterating", total=limit)
            state = {"iter": 0}

            def _observe(p):
                # inside a solve the runs are the library's, so this is the only place they
                # can be seen. Reported as runs, with the iteration as context.
                bar.update(done=p["n_completed"], total=p["n_total"] or None,
                           iter=state["iter"],
                           **({"failed": p["n_failed"]} if p["n_failed"] else {}))
                return True

            observing = not isinstance(bar, Progress) or type(bar) is not Progress
            if observing:
                self._lib.set_run_observer(_observe, min_interval_sec=0.1)
            try:
                while n < limit and not self.should_terminate:
                    state["iter"] = self.iteration + 1
                    step = self.solve()
                    n += 1
                    fields = {"iter": step.iter}
                    if step.phi_mean == step.phi_mean:      # not NaN
                        fields["phi"] = step.phi_mean
                    if step.retried:
                        fields["retry"] = "yes"
                    # the observer has been driving `done` in run units; hand the line back to
                    # iteration units now that one is complete
                    bar.start("iterating", total=limit)
                    bar.update(done=n, **fields)
                    yield step
            finally:
                if observing:
                    self._lib.set_run_observer(None)

    def finalize(self) -> None:
        with self._q():
            self._lib.finalize()

    def close(self) -> None:
        """Release the handle and stop any workers this session started.

        Idempotent. Calling it is good practice and `with` calls it for you, but neither is
        required: a dropped session is cleaned up by its finalizer.
        """
        # One code path for explicit close(), __exit__ and collection; finalize() runs the
        # callback at most once, so the repeat calls this invites are free.
        self._finalizer()
        self._workers = []
        self._lib = None

    def __enter__(self: _ToolT) -> _ToolT:
        return self

    def __exit__(self, *exc):
        self.close()

    # -- state (properties: cheap, no arguments) -------------------------------------------

    @property
    def lib(self) -> PestppLib:
        """The thin binding underneath, for anything this layer does not wrap."""
        return self._lib

    @property
    def workdir(self) -> str:
        return self._workdir

    @property
    def iteration(self) -> int:
        return self._lib.get_iteration()

    @property
    def is_initialized(self) -> bool:
        return self._initialized

    @property
    def should_terminate(self) -> bool:
        # wrapped because the tools print their phi-termination report from this call
        with self._q():
            return self._lib.should_terminate()

    @property
    def real_names(self) -> np.ndarray:
        """Realization names, as an array so it can index boolean masks."""
        return np.array(self._lib.get_ensemble_row_names(PAR_EN))

    @property
    def adj_par_names(self) -> list:
        """The ADJUSTABLE parameters - the raw ensemble's columns.

        Named to match pyemu, where ``Pst.adj_par_names`` is the adjustable subset and
        ``Pst.par_names`` is everything. This used to be called ``par_names``, which meant
        the opposite of what a pyemu user would read it as - and once ``.pst`` appeared on
        this object, ``ies.par_names`` and ``ies.pst.par_names`` sat one dot apart giving
        different answers. Same word, two lists, no warning.
        """
        return self._lib.get_ensemble_col_names(PAR_EN)

    @property
    def par_names(self) -> list:
        """EVERY control-file parameter, including fixed and tied - like ``Pst.par_names``.

        This is the snapshot's column set, so it pairs with :meth:`par_df`.
        See :attr:`adj_par_names` for the adjustable subset the raw ensemble holds.
        """
        return list(self.par_df().columns)

    @property
    def obs_names(self) -> list:
        return self._lib.get_ensemble_col_names(OBS_EN)

    @property
    def n_reals(self) -> int:
        return len(self._lib.get_ensemble_row_names(PAR_EN))

    @property
    def par_transform(self) -> str:
        """Which space the raw parameter view is in: 'ctl', 'num' or 'model'."""
        return self._lib.get_par_transform_status()

    @property
    def run_manager(self) -> str:
        return self._lib.get_run_manager()

    @property
    def supports_live_control(self) -> bool:
        """Whether runs can be watched and cancelled mid-batch. PANTHER only."""
        return self._lib.supports_live_control()

    @property
    def phi(self) -> float:
        """Mean measurement phi. ies and da only."""
        return self.get_phi()["mean"]

    # -- things that take arguments --------------------------------------------------------

    #: phi types, in the order phi_df() lays them out
    PHI_TYPES = (("meas", PHI_MEAS), ("actual", PHI_ACTUAL), ("composite", PHI_COMPOSITE),
                 ("regul", PHI_REGUL), ("noise", PHI_NOISE))

    def get_phi(self, phi_type: int = PHI_MEAS) -> dict:
        """mean/std/min/max of phi. ies and da only."""
        if not self._has_phi:
            raise PestppError("{0} has no phi".format(type(self).__name__.lower()))
        return self._lib.get_phi_summary(phi_type)

    def phi_df(self, phi_type: int | None = None, lower: bool = False) -> pd.DataFrame:
        """Phi for every realization, one column per phi type.

        Index is the realization, columns are meas / actual / composite / regul / noise. A
        type that is not in play - regul outside a regularized run, noise under
        ``ies_no_noise`` - comes back as an all-NaN column rather than being missing, so the
        frame's shape does not depend on how the case was set up.

        Pass ``phi_type`` for a single named column instead of all five.
        """
        if not self._has_phi:
            raise PestppError("{0} has no phi".format(type(self).__name__.lower()))
        wanted = self.PHI_TYPES
        if phi_type is not None:
            wanted = [(n, t) for n, t in self.PHI_TYPES if t == phi_type]
            if not wanted:
                raise PestppError("unknown phi type {0}".format(phi_type))
        cols = {}
        for name, t in wanted:
            vals, names = self._lib.get_phi_vector(t)
            cols[name] = pd.Series(vals, index=names, dtype=float)
        return _maybe_lower(pd.DataFrame(cols), lower, axis="index")

    @property
    def obs_groups(self) -> pd.Series:
        """Observation name -> group, for every observation."""
        return pd.Series(self._lib.get_obs_groups(), index=self.obs_names, name="obgnme")

    def get_obs_group(self, name) -> str:
        """The group one observation belongs to. Case-insensitive."""
        return self.obs_groups.loc[_up(name)]

    def phi_obs_df(self, phi_type: int = PHI_MEAS, lower: bool = False) -> pd.DataFrame:
        """The individual terms phi is the sum of, per realization per observation.

        Index is the realization, columns the observations, values the squared weighted
        residual. Summing a row gives that realization's phi. MEAS and ACTUAL only.

        This is the finer-grained view behind :meth:`phi_group_df` - use it when the question
        is which *observations* are fitting badly rather than which group.
        """
        if not self._has_phi:
            raise PestppError("{0} has no phi".format(type(self).__name__.lower()))
        arr, rows, cols = self._lib.get_phi_residuals(phi_type)
        return _maybe_lower(pd.DataFrame(arr, index=rows, columns=cols), lower)

    def phi_group_df(self, phi_type: int = PHI_MEAS, lower: bool = False) -> pd.DataFrame:
        """Phi by observation group, per realization.

        Index is the realization, columns the observation groups - the same decomposition the
        tools print as their "observation group phi summary", but as a frame you can sort and
        plot. Grouping happens here rather than in C++ because the residual terms are the
        primitive and pandas is better at the aggregation.

        Only groups with a non-zero-weighted observation contribute, since a zero weight makes
        the term zero by construction.
        """
        obs = self.phi_obs_df(phi_type)
        if obs.empty:
            return obs
        groups = self.obs_groups.reindex(obs.columns)
        # a residual column with no group mapping would silently vanish in the groupby
        missing = groups.isna().sum()
        if missing:
            raise PestppError("{0} residual columns have no observation group".format(missing))
        return _maybe_lower(obs.T.groupby(groups.values).sum().T, lower)

    # -- weights ---------------------------------------------------------------------------

    @property
    def obs_weights(self) -> pd.Series:
        """The control-file weight vector, indexed by observation name."""
        return pd.Series(self._lib.get_obs_weights(), index=self.obs_names, name="weight")

    def set_obs_weights(self, weights, broadcast: bool = True) -> None:
        """Set the control-file weight vector. Accepts a Series, a dict, or a scalar.

        A scalar sets every observation. A Series or dict sets only the names it carries, so
        you can adjust a few without restating the rest.

        ``broadcast`` (the default) also pushes the new vector into every row of the weights
        ensemble, which is what makes the change take effect: initialize() seeds that ensemble
        from the vector once, and phi is computed from the ENSEMBLE thereafter. Pass
        ``broadcast=False`` if you are managing per-realization weights yourself through
        :meth:`weights_view`, since broadcasting would overwrite them.
        """
        if np.isscalar(weights):
            names = list(self.obs_names)
            vals = [float(weights)] * len(names)
        else:
            series = pd.Series(weights, dtype=float)
            names, vals = _up_all(series.index), list(series.values)
        self._lib.set_obs_weights(names, vals)
        if broadcast:
            self._lib.broadcast_weights()
            # phi is cached from the last update, so without this the caller reads the old
            # number back and concludes the weight change did nothing
            self.update_phi()

    def reinflate(self, factor: float = 1.0, num_reals: int = 0,
                  center_on_min_phi: bool | None = None) -> None:
        """Rebuild the parameter ensemble from the PRIOR's spread, re-centred on the current one.

        A smoother narrows the ensemble every iteration. After a few, the realizations can
        agree with each other far more than the data justifies - and an over-tight ensemble
        cannot respond to anything new, because there is no variance left to move. Reinflation
        keeps the location the assimilation has found and puts back some of the variance it
        started with.

        The tool does this on a schedule (``ies_n_iter_reinflate``). Calling it explicitly is
        for the case the schedule cannot express: reinflating at the moment new observations
        are brought in, so the ensemble has room to react to them.

        ``factor`` is in (0, 1] - 1.0 restores the full prior spread.

        ``num_reals`` of 0 keeps the current realization count; otherwise it is how many
        realizations to end up with. Realizations are SELECTED FROM THE PRIOR rather than
        generated, so the prior's size is a hard ceiling and asking for more raises. To grow
        the ensemble mid-run, start with a prior big enough for the largest size you will ask
        for and begin the run on a subset of it::

            ies = Ies.from_pst("pest.pst", ies_num_reals=50,
                               ies_reinflate_num_reals="10")   # prior 50, working ensemble 10
            ies.initialize()
            ...
            ies.reinflate(num_reals=20)                        # 10 -> 20, out of the 50

        The SIGN of ``num_reals`` chooses where the spread comes from: positive uses the
        prior's own anomalies scaled by ``factor``, negative resamples the CURRENT ensemble's
        anomalies instead, adding prior anomalies scaled by ``factor`` when it is below 1.

        ``center_on_min_phi`` is what the new spread is centred on: ``None`` follows
        ``ies_n_iter_reinflate`` the way the built-in loop does, ``False`` forces the ensemble
        mean, ``True`` forces the minimum-phi realization - the aggressive form.

        This RUNS the reinflated ensemble, so it costs ``num_reals`` model runs.
        """
        center = -1 if center_on_min_phi is None else int(bool(center_on_min_phi))
        with self._q():
            self._lib.reinflate_ensemble(factor, num_reals, center)

    def update_phi(self) -> None:
        """Recompute phi from the current ensembles and weights.

        Phi is cached - the algorithm refreshes it at its own points. Call this after writing
        an ensemble or weights view if you want the phi accessors to reflect the change.
        """
        with self._q():
            self._lib.update_phi()

    def weights_df(self, lower: bool = False) -> pd.DataFrame:
        """The weights ENSEMBLE - one weight per observation per realization, as a copy."""
        arr, token = self._lib.get_ensemble_view(WEIGHTS_EN)
        self._lib.release_view(token)          # a copy is taken below; nothing stays borrowed
        return _named(_maybe_lower(pd.DataFrame(
            arr.copy(),
            index=self._lib.get_ensemble_row_names(WEIGHTS_EN),
            columns=self._lib.get_ensemble_col_names(WEIGHTS_EN)), lower),
            "realization", "obsnme")

    def noise_df(self, lower: bool = False) -> pd.DataFrame:
        """The measurement-noise ensemble - one noise realization per observation, as a copy.

        These are what phi is measured AGAINST: `PHI_MEAS` is the residual between the
        simulated ensemble and this one, rather than against the single observed value. Its
        columns are the ACTIVE observations, so an observation at zero weight is absent - and
        activating one generates its noise here, which is the thing worth checking after a
        staged reweight.

        Empty when the run was configured with ``ies_no_noise``.
        """
        arr, token = self._lib.get_ensemble_view(NOISE_EN)
        self._lib.release_view(token)          # a copy is taken below; nothing stays borrowed
        return _named(_maybe_lower(pd.DataFrame(
            arr.copy(),
            index=self._lib.get_ensemble_row_names(NOISE_EN),
            columns=self._lib.get_ensemble_col_names(NOISE_EN)), lower),
            "realization", "obsnme")

    @contextmanager
    def noise_view(self):
        """Zero-copy view of the noise ensemble, valid only in this block."""
        yield from self._view(NOISE_EN)

    @contextmanager
    def weights_view(self):
        """Zero-copy view of the weights ENSEMBLE, valid only in this block.

        This is the per-realization path: write here when weights should differ between
        realizations. :meth:`set_obs_weights` is the one-weight-per-observation path.
        """
        yield from self._view(WEIGHTS_EN)

    def phi_summary_df(self) -> pd.DataFrame:
        """mean/std/min/max for every phi type, one row per type."""
        if not self._has_phi:
            raise PestppError("{0} has no phi".format(type(self).__name__.lower()))
        nan = {k: float("nan") for k in ("mean", "std", "min", "max")}
        rows = {}
        for name, t in self.PHI_TYPES:
            try:
                # a phi type with no realizations reports the handler's sentinels from
                # get_min/get_max (DBL_MAX and -1e30) rather than anything meaningful, so
                # decide emptiness from the vector and report NaN across the row
                vals, _ = self._lib.get_phi_vector(t)
                rows[name] = nan if len(vals) == 0 else self._lib.get_phi_summary(t)
            except PestppError:
                rows[name] = nan
        return pd.DataFrame(rows).T[["mean", "std", "min", "max"]]

    def get_run_states(self, run_ids=None) -> pd.DataFrame:
        """One row per run: id, status, elapsed, failures, host. PANTHER only."""
        return pd.DataFrame(self._lib.get_run_states(run_ids))

    def get_workers(self) -> pd.DataFrame:
        """One row per connected worker, with its completed/failed/timed-out run counts."""
        rows = []
        for i in range(self._lib.get_worker_count()):
            w = dict(self._lib.get_worker_state(i))
            w["index"] = i
            for which, label in ((WORKER_COMPLETED, "n_completed"),
                                 (WORKER_FAILED, "n_failed"),
                                 (WORKER_TIMED_OUT, "n_timed_out")):
                w[label] = len(self._lib.get_worker_run_history(i, which))
            rows.append(w)
        return pd.DataFrame(rows)

    def cancel_runs(self, run_ids) -> int:
        """Give up on these runs. Returns how many were actually cancelled."""
        return self._lib.cancel_runs(list(run_ids))

    def release_workers(self, worker_idxs=None, reap_timeout: float = 5.0) -> int:
        """Hand workers back so their compute can go elsewhere. Returns how many actually went.

        `worker_idxs` are rows of get_workers(); None releases every worker. Releasing is a
        request, not a command, so the returned count is the authority - asking for eight and
        getting three is a normal outcome.

        A worker that is MID-RUN is released too, and its run is rescheduled: put back at the
        front of the queue for another worker to pick up, not failed and not counted against
        max_n_failure. Use cancel_runs() when the judgement is about the run; use this when it
        is about the machine.

        Legal from inside a run observer, which is the case that motivates it - watching a
        batch drain, seeing workers go idle with a couple of runs left, and handing the idle
        ones back without waiting for the batch to end.
        """
        n = self._lib.release_workers(list(worker_idxs) if worker_idxs else None)
        if n:
            self._reap_released_workers(n, reap_timeout)
        return n

    def _reap_released_workers(self, n_released: int, timeout: float) -> None:
        """Forget worker processes that have exited, so teardown does not chase ghosts.

        There is no mapping from a run-manager worker index to one of OUR subprocesses: agents
        connect in whatever order they happen to come up, and some of them may not be ours at
        all. What IS reliable is that a released agent exits, so wait briefly for up to
        n_released of them to do so and drop the ones that have. Anything still alive stays
        tracked and is still terminated at teardown - the safe direction to be wrong in.

        Mutates the list IN PLACE: the finalizer registered in __init__ holds this exact list
        object, so rebinding it would leave teardown working from a stale copy.
        """
        deadline = time.monotonic() + max(0.0, timeout)
        while True:
            alive = [p for p in self._workers if p.poll() is None]
            if (len(self._workers) - len(alive) >= n_released) or (time.monotonic() >= deadline):
                self._workers[:] = alive
                return
            time.sleep(0.05)

    # -- data ------------------------------------------------------------------------------

    def par_df(self, lower: bool = False) -> pd.DataFrame:
        """Every control-file parameter in CTL space, as a copy. Ready for pyemu."""
        vals, rows, cols = self._lib.get_par_snapshot()
        return _named(_maybe_lower(pd.DataFrame(vals, index=rows, columns=cols), lower),
                      "realization", "parnme")

    def set_par_df(self, df, enforce: str = "raise") -> None:
        """Push parameter values back. Matched by NAME, so order does not matter.

        Accepts a DataFrame or a :class:`pyemu.ParameterEnsemble`.

        ``enforce`` decides what happens to values outside a parameter's bounds, and the
        default is to refuse rather than to proceed. That is not fussiness: pest++ maps a
        control-file value into the tool's transform space, and an out-of-bounds value on a
        LOG parameter becomes **NaN**, which is then run through the model. pyemu draws go
        out of bounds routinely - that is what ``ParameterEnsemble.enforce()`` exists for -
        so this is a live path, not a hypothetical one.

            "raise"  refuse, naming the offenders          (default)
            "reset"  clip to the bound, pyemu's how="reset"
            False    send exactly what you passed
        """
        if hasattr(df, "to_dataframe"):        # a pyemu Ensemble
            df = df.to_dataframe()
        df = df.copy()
        if enforce:
            df = self._enforce_bounds(df, how=enforce)
        self._lib.set_par_snapshot(df.values, _up_all(df.index), _up_all(df.columns))

    def _enforce_bounds(self, df, how="raise"):
        """Bounds check/clip against the control file, in the frame's own case."""
        par = self.pst.parameter_data
        lb = par.parlbnd.copy()
        ub = par.parubnd.copy()
        # the frame may be in either case; the Pst is lowercase, so meet it there
        cols = [str(c) for c in df.columns]
        lower_cols = [c.lower() for c in cols]
        known = [i for i, c in enumerate(lower_cols) if c in lb.index]
        if not known:
            return df
        idx = [lower_cols[i] for i in known]
        sub = df.iloc[:, known]
        lo = lb.loc[idx].values.astype(float)
        hi = ub.loc[idx].values.astype(float)
        below = sub.values < lo
        above = sub.values > hi
        n_bad = int(below.sum() + above.sum())
        if n_bad == 0:
            return df
        if str(how).lower() == "reset":
            vals = np.clip(sub.values, lo, hi)
            df.iloc[:, known] = vals
            return df
        bad_cols = sorted({idx[j] for j in np.where(below.any(axis=0) | above.any(axis=0))[0]})
        raise PestppError(
            "{0} value(s) across {1} parameter(s) are outside their control-file bounds, e.g. "
            "{2}. pest++ would transform an out-of-bounds LOG parameter to NaN and run it. "
            "Pass enforce='reset' to clip (pyemu's ParameterEnsemble.enforce() does the same), "
            "or enforce=False to send them anyway.".format(n_bad, len(bad_cols), bad_cols[:5]))

    # -- pyemu ------------------------------------------------------------------------------
    #
    # The point of this block is that someone who already knows pyemu should not have to learn
    # a second vocabulary. `pst`, `pe`, `oe` and `results` mean here exactly what they mean
    # there, the frames carry the axis names pyemu's own result handler gives them, and
    # everything is lowercase because that is what pest++ writes to its ensemble csvs and
    # therefore what `pyemu.Results` hands back.

    @property
    def pst(self) -> "pyemu.Pst":
        """The control file as a :class:`pyemu.Pst` - bounds, transforms, obsvals, groups.

        Parsed once and cached, because the static metadata cannot change mid-run. The one
        thing that CAN change is the weights, and those are re-synced from the library on
        every access, so ``pst.observation_data.weight`` is never stale and
        ``oe.phi_vector`` agrees with the phi the tool reports.

        Note this is the WEIGHT VECTOR, not the weights ensemble - see :meth:`weights_df` for
        the per-realization form, which a Pst has no way to represent.
        """
        if self._pst is None:
            self._pst = pyemu.Pst(os.path.join(self._workdir, self._pst_file))
        try:
            live = dict(zip([n.lower() for n in self.obs_names], self.obs_weights))
            obs = self._pst.observation_data
            hit = [n for n in obs.index if n in live]
            if hit:
                obs.loc[hit, "weight"] = [live[n] for n in hit]
        except PestppError:
            # before initialize() there may be no observation ensemble to read weights from;
            # the control-file values already in the Pst are the right answer then
            pass
        return self._pst

    @property
    def pe(self) -> "pyemu.ParameterEnsemble":
        """The parameter ensemble as a :class:`pyemu.ParameterEnsemble`, in CTL space.

        Untransformed (``istransformed=False``), which is what CTL space is - so
        ``.enforce()``, ``.to_binary()``, ``.covariance_matrix()`` and the plotting helpers
        all behave the way they do on an ensemble read off disk.
        """
        return pyemu.ParameterEnsemble(pst=self.pst, df=self.par_df(lower=True),
                                       istransformed=False)

    @property
    def oe(self) -> "pyemu.ObservationEnsemble":
        """The observation ensemble as a :class:`pyemu.ObservationEnsemble`.

        ``oe.phi_vector`` is then computable in pyemu and comparable with :attr:`phi_actual`.
        They are the same quantity by different routes, which makes it a useful cross-check.
        """
        return pyemu.ObservationEnsemble(pst=self.pst, df=self.obs_df(lower=True))

    def set_pe(self, pe, enforce: str = "reset") -> None:
        """Push a :class:`pyemu.ParameterEnsemble` (or DataFrame) back into the tool.

        Defaults to ``enforce="reset"`` rather than "raise", because the thing a caller most
        often hands this is a fresh pyemu draw - and a draw being out of bounds is expected,
        not exceptional. This is the API-side equivalent of ``pe.enforce()``.
        """
        self.set_par_df(pe, enforce=enforce)

    @property
    def results(self):
        """A :class:`pyemu.Results` over this session's working directory.

        The on-disk history, with the vocabulary a pyemu user already has: ``.ies.paren0``,
        ``.ies.phiactual``, ``.mou.dvpop``, and so on. The live in-memory state is what
        :attr:`pe`, :attr:`oe` and :meth:`par_df` give you; this is everything the tool has
        written so far, including the iterations already behind you.
        """
        case = self._pst_file[:-4] if self._pst_file.lower().endswith(".pst") else self._pst_file
        return pyemu.Results(self._workdir, case=case)

    def obs_df(self, lower: bool = False) -> pd.DataFrame:
        """The observation ensemble, as a copy."""
        arr, token = self._lib.get_ensemble_view(OBS_EN)
        self._lib.release_view(token)          # a copy is taken below; nothing stays borrowed
        return _named(_maybe_lower(pd.DataFrame(
            arr.copy(),
            index=self._lib.get_ensemble_row_names(OBS_EN),
            columns=self._lib.get_ensemble_col_names(OBS_EN)), lower),
            "realization", "obsnme")

    @contextmanager
    def par_view(self):
        """A zero-copy numpy view of the raw parameter ensemble, valid only in this block.

        Writes through it reach the tool directly, which is the point. It is a context manager
        rather than a property on purpose: the array points into the tool's live buffer, so
        holding one past a resize is a use-after-free. Keep anything you need with ``.copy()``,
        or use :meth:`par_df`.

        Note this is the RAW ensemble - usually NUM space, adjustable parameters only. See
        :attr:`par_transform`. :meth:`par_df` is the CTL-space, everything-included form.
        """
        yield from self._view(PAR_EN)

    @contextmanager
    def obs_view(self):
        """A zero-copy numpy view of the observation ensemble, valid only in this block."""
        yield from self._view(OBS_EN)

    def _view(self, which):
        arr, token = self._lib.get_ensemble_view(which)
        proxy = EnsembleViewProxy(arr, self._lib, token)
        try:
            yield proxy
        finally:
            proxy._expire()

    # -- membership ------------------------------------------------------------------------

    def drop_realizations(self, names) -> None:
        """Drop these from every coupled ensemble at once. Names are case-insensitive."""
        self._lib.drop_realizations(_up_all(names))

    def keep_realizations(self, names) -> None:
        """Keep only these, in this order. Also a reorder. Names are case-insensitive."""
        self._lib.keep_realizations(_up_all(names))

    # -- running things yourself -----------------------------------------------------------

    def queue_runs(self, reals=None) -> int:
        """Queue runs with the run manager. Returns how many were queued.

        Normally the current ensemble. During a deferred :meth:`solve` it is that solve's
        outstanding batch instead - the candidates, or the remaining realizations once
        :meth:`finish_solve` has asked for them - so the same call works inside the loop.

        ``reals`` names the realizations to run, and is meaningful only for the candidate
        batch: it REPLACES the subset the algorithm picked, so whatever is left unnamed
        becomes the remainder that :meth:`finish_solve` asks for later. ``None`` keeps the
        algorithm's own choice.
        """
        with self._q():
            if reals is not None:
                self._queued = self._lib.queue_runs_subset(_up_all(reals))
            else:
                self._queued = self._lib.queue_runs()
        return self._queued

    def run(self, slice_seconds: float = 0.05, callback=None, progress=False) -> None:
        """Drive the queued runs to completion.

        ``callback`` is called with this tool after each slice, which is where a cancel
        decision goes. PANTHER only - serial and external finish the whole batch in the first
        slice, so a callback fires at most once there.

        ``progress=True`` draws a live bar of completed/failed runs, in a terminal or a
        notebook. Pass a :class:`~pestpp_progress.Progress` instead to render it yourself.
        The counts come from the run manager, so they are only as live as it is: PANTHER
        reports them as workers land, the serial manager reports the whole batch at once.
        """
        bar = _progress_auto(progress)
        # Asked ONCE, and only when there is a bar to feed. The counters are PANTHER-only, so
        # on the serial manager every one of these calls throws - and a call that throws is
        # not free: it still enters and leaves the library's working directory on the way. A
        # per-run guaranteed exception on the most common path is not something to leave in
        # for the sake of a progress bar nobody asked for.
        live = (not isinstance(bar, Progress) or type(bar) is not Progress) \
            and self.supports_live_control
        with self._q():
            self._lib.begin_batch()
            try:
                bar.start("running", total=self._queued or None)
                while not self._lib.run_slice(slice_seconds):
                    if live:
                        self._update_run_progress(bar)
                    if callback is not None:
                        callback(self)
            finally:
                self._lib.end_batch()
                if live:
                    self._update_run_progress(bar)
                bar.close()

    def _update_run_progress(self, bar) -> None:
        """Feed the run manager's counters to a progress renderer, if it can supply them."""
        try:
            stats = self._lib.get_run_time_stats()
        except PestppError:
            return                      # not every run manager keeps these; a bar is not
        fields = {}                     # worth failing a run over
        if stats.get("failed"):
            fields["failed"] = stats["failed"]
        if stats.get("running"):
            fields["running"] = stats["running"]
        bar.update(done=stats.get("completed", 0), **fields)

    def process_runs(self) -> int:
        """Process the completed runs into the observation ensemble. Returns how many failed."""
        with self._q():
            n = self._lib.process_runs()
        self._queued = 0
        return n

    def run_ensemble(self, callback=None) -> int:
        """queue_runs, run and process_runs in one call. Returns the number of failed runs."""
        self.queue_runs()
        self.run(callback=callback)
        return self.process_runs()

    # -- internals -------------------------------------------------------------------------

    def _step(self, retried: bool, pending_runs: int = 0) -> IterationStep:
        step = IterationStep(iter=self.iteration, n_reals=self.n_reals, retried=retried,
                             pending_runs=pending_runs)
        if self._has_phi:
            try:
                p = self.get_phi()
                step.phi_mean, step.phi_std = p["mean"], p["std"]
                step.phi_min, step.phi_max = p["min"], p["max"]
            except PestppError:
                pass          # phi is not meaningful yet; leave the NaNs
        return step


class Candidate:
    """One upgrade candidate of an open deferred solve.

    A thin handle rather than a copy: the frame and view calls go straight at the tool's
    candidate ensemble, so a write through :meth:`par_view` is what gets run. Valid only until
    the solve moves on - :meth:`_Tool.finish_solve` releases the candidates, and using a handle
    afterwards raises rather than reading freed memory.
    """

    def __init__(self, tool, index: int):
        self._tool = tool
        self._lib = tool._lib
        self.index = index
        #: the lambda (or mda factor) and backtrack scale this candidate was generated with.
        #: mou has neither and reports 0.
        self.inflation, self.backtrack = self._lib.get_candidate_info(index)

    @property
    def _en_id(self) -> int:
        from pestpp_lib import CANDIDATE_EN
        return CANDIDATE_EN + self.index

    def par_df(self, lower: bool = False) -> pd.DataFrame:
        """This candidate's parameters, as a copy."""
        arr, token = self._lib.get_ensemble_view(self._en_id)
        try:
            out = pd.DataFrame(arr.copy(),
                               index=self._lib.get_ensemble_row_names(self._en_id),
                               columns=self._lib.get_ensemble_col_names(self._en_id))
        finally:
            self._lib.release_view(token)
        return _named(_maybe_lower(out, lower), "realization", "parnme")

    @contextmanager
    def par_view(self):
        """Zero-copy view of this candidate, valid only in this block. Writes reach the run."""
        yield from self._tool._view(self._en_id)

    @property
    def shape(self):
        return self.par_df().shape

    def __repr__(self):
        return ("Candidate(index={0}, inflation={1:.6g}, backtrack={2:.6g})"
                .format(self.index, self.inflation, self.backtrack))


class Ies(_Tool):
    """Iterative ensemble smoother."""
    _tool_id = TOOL_IES


class Da(_Tool):
    """Data assimilation. A SINGLE cycle - multi-cycle is not exposed yet."""
    _tool_id = TOOL_DA

    @staticmethod
    def _agent_exe():
        return "pestpp-da"


class Glm(_Tool):
    """pestpp-glm: Gauss-Levenberg-Marquardt, through the same loop the executable runs.

    The odd tool in this package, and worth being explicit about why. The other four carry a
    POPULATION - an ensemble of realizations - so ``par_df()``, the ensemble views and the
    phi-over-realizations calls all mean something. glm carries ONE parameter vector through a
    Jacobian and an upgrade. Those inherited calls therefore refuse rather than hand back a
    one-row frame, which would read as "an ensemble of one realization" and is not what glm
    computed.

    What glm offers instead:

    ``par_vector()`` / ``set_par_vector()``   the single parameter vector, as a pandas Series
    ``jco()``                                 the Jacobian, as a DataFrame, ndarray or pyemu Jco
    """
    _tool_id = TOOL_GLM
    #: one phi, not a phi over realizations, so the inherited phi helpers stay out of the way
    _has_phi = False

    @staticmethod
    def _agent_exe():
        return "pestpp-glm"

    @property
    def n_reals(self) -> int:
        """Always 0. glm carries a parameter VECTOR, so there are no realizations to count.

        0 rather than a refusal because this feeds IterationStep, which solve() returns on
        every iteration - raising here would make the ordinary loop unusable to report the
        absence of a thing glm never had.
        """
        return 0

    def jacobian_prepare(self, calc_init_obs: bool = False) -> int:
        """Queue the Jacobian runs WITHOUT running them. Returns how many are outstanding.

        The three stages - prepare, run, process - let you own the jco batch::

            n = glm.jacobian_prepare()
            glm.jacobian_run()          # or drive the run manager yourself
            glm.jacobian_process()
            jco = glm.jco()

        ``solve()`` composes exactly these internally, so both paths run the same code.
        """
        return self._lib.jacobian_prepare(calc_init_obs)

    def jacobian_run(self) -> None:
        """Run the queued Jacobian batch."""
        with self._q():
            self._lib.jacobian_run()

    def jacobian_process(self) -> bool:
        """Harvest the completed Jacobian runs. The jco is available from :meth:`jco` after."""
        with self._q():
            return self._lib.jacobian_process()

    def par_vector(self, which: str = "current", lower: bool = False) -> pd.Series:
        """The current (or ``"optimum"``) parameter vector in CTL space, as a labelled Series."""
        which_id = {"current": 0, "optimum": 1}.get(str(which).lower())
        if which_id is None:
            raise ValueError("which must be 'current' or 'optimum', not {0!r}".format(which))
        names, vals = self._lib.get_par_vector(which_id)
        idx = [n.lower() for n in names] if lower else names
        return pd.Series(vals, index=idx, name=str(which).lower())

    def set_par_vector(self, values) -> None:
        """Push a parameter vector back. Accepts a Series/dict, matched BY NAME.

        Partial is fine: only what you name is changed. A name glm does not hold is an error
        rather than a silent skip, because it means you and the tool disagree about the
        parameter set - and that is worth hearing about at the call, not three iterations later.
        """
        if isinstance(values, pd.Series):
            names, vals = list(values.index), values.values
        elif isinstance(values, dict):
            names, vals = list(values.keys()), list(values.values())
        else:
            raise TypeError("set_par_vector wants a pandas Series or a dict of name -> value, "
                            "not {0}; a bare array has no names to match on"
                            .format(type(values).__name__))
        self._lib.set_par_vector([str(n) for n in names], vals)

    def jco(self, as_: str = "pandas", lower: bool = False):
        """The Jacobian. ``as_`` is ``"pandas"``, ``"numpy"`` or ``"pyemu"``.

        Rows are simulated observations, columns the adjustable parameters it was built over.
        A dense COPY in every case - the matrix is sparse internally, so there is no view to
        hand out - and it exists only once ``initialize()`` has built it.

        ``"numpy"`` returns ``(values, row_names, col_names)`` rather than a bare array, because
        an unlabelled Jacobian is the easiest thing in this package to line up wrongly.
        """
        data, rnames, cnames = self._lib.get_jacobian()
        if lower:
            rnames = [r.lower() for r in rnames]
            cnames = [c.lower() for c in cnames]
        kind = str(as_).lower()
        if kind == "numpy":
            return data, rnames, cnames
        if kind == "pandas":
            return pd.DataFrame(data, index=rnames, columns=cnames)
        if kind == "pyemu":
            import pyemu
            # pyemu keeps names lower-cased; hand it that regardless of `lower`, which governs
            # what YOU get back, not what pyemu needs internally
            return pyemu.Jco(x=data, row_names=[r.lower() for r in rnames],
                             col_names=[c.lower() for c in cnames])
        raise ValueError("as_ must be 'pandas', 'numpy' or 'pyemu', not {0!r}".format(as_))


class _ChanceMixin:
    """Live chance-constraint / risk control, shared by the tools that support it.

    mou and sqp both carry uncertainty through PARAMETER STACKS: an ensemble of the uncertain
    parameters is run, the resulting spread shifts the constraints and objectives, and `risk`
    says how conservatively to shift them. So mou does have parameters and realizations - they
    live here rather than in the population.

    Everything below is asked of the tool, live, so it reflects what `Constraints` will do -
    including the clamp. That matters: the tool does not use the number you set, it
    uses the number these accessors report, and the two differ at the edges.

    The four derived flags below are ASKED OF THE TOOL rather than worked out here. The rule
    that turns options into behaviour is genuinely intricate - robust switches chance off
    entirely, stacks beat fosm but only when ``opt_std_weights`` is off, risk is clamped - and
    a second copy of it in python would be a copy that drifts. ``Constraints`` derives it live on every call, so
    what these report is what the tool will do.
    """

    @property
    def stack_status(self) -> dict:
        """The whole derived chance state in one call - the flags below share this."""
        return self._lib.get_stack_status()

    @property
    def use_robust(self) -> bool:
        """True for ROBUST optimization - ``opt_use_robust``, pestpp-sqp only.

        Each decision-variable realization is paired with its own uncertain parameter
        realization and the ensemble is optimized as it stands. NOTHING is risk shifted, so
        :attr:`use_chance` and :attr:`use_fosm` are both False whenever this is True, and the
        tool refuses ``opt_use_robust`` alongside a non-neutral ``opt_risk``.
        """
        return self.stack_status["use_robust"]

    @property
    def risk(self) -> float:
        """The risk value the tool will USE - ``opt_risk``, clamped to [0.001, 0.999].

        Read this rather than the option: setting ``opt_risk`` to 1.5 leaves the option at 1.5
        and the tool using 0.999.
        """
        return self.stack_status["risk"]

    @property
    def use_chance(self) -> bool:
        """True when chance constraints are active at all - i.e. :attr:`risk` is not 0.5.

        0.5 is risk-neutral and switches the whole machinery off, which is why a run with
        stacks configured but risk left at 0.5 does no chance work at all.
        """
        return self.stack_status["use_chance"]

    @property
    def use_fosm(self) -> bool:
        """True when uncertainty comes from FOSM rather than from stacks.

        Stacks win: any of ``opt_stack_size > 0``, ``opt_par_stack`` or ``opt_obs_stack`` set
        turns FOSM off, unless ``opt_std_weights`` is on. On a robust run it is always off.
        """
        return self.stack_status["use_fosm"]

    # -- the stacks themselves -------------------------------------------------------------

    def _stack_df(self, ensemble_id: int, lower: bool, what: str) -> pd.DataFrame:
        """One stack as a copy. Empty frame when the run does not use stacks."""
        rows = self._lib.get_ensemble_row_names(ensemble_id)
        cols = self._lib.get_ensemble_col_names(ensemble_id)
        if (len(rows) == 0) or (len(cols) == 0):
            # legitimately empty on a fosm or risk-neutral run - stack_status says which,
            # and an empty frame with no columns is the honest shape for "there is no stack"
            return _named(pd.DataFrame(index=rows, columns=cols, dtype=float), "realization",
                          what)
        arr, token = self._lib.get_ensemble_view(ensemble_id)
        try:
            df = pd.DataFrame(np.array(arr, copy=True), index=rows, columns=cols)
        finally:
            self._lib.release_view(token)
        return _named(_maybe_lower(df, lower), "realization", what)

    def stack_pe(self, lower: bool = False) -> pd.DataFrame:
        """The PARAMETER stack - the realizations that get run to measure uncertainty.

        Empty unless this is a stack-based chance run; :attr:`stack_status` says why.
        """
        return self._stack_df(STACK_PAR_EN, lower, "parnme")

    def stack_oe(self, lower: bool = False) -> pd.DataFrame:
        """The OBSERVATION stack - the results of running :meth:`stack_pe`.

        This is the ensemble the risk shift is computed from, so it is the one to look at when
        asking why a constraint was shifted as far as it was.
        """
        return self._stack_df(STACK_OBS_EN, lower, "obsnme")

    def nested_pe(self, lower: bool = False) -> pd.DataFrame:
        """The nested parameter stack. Only filled when chance is evaluated at several points
        in decision-variable space (``opt_chance_points`` = "all")."""
        return self._stack_df(NESTED_PAR_EN, lower, "parnme")

    def member_stacks(self) -> list:
        """Which members have their own observation stack.

        Empty unless ``opt_chance_points`` is "all" - otherwise one stack serves every point
        and :meth:`stack_oe` is the only stack there is.
        """
        if self._lib.get_member_stack_count() == 0:
            return []
        return self._lib.get_member_stack_names()

    def member_stack_oe(self, member, lower: bool = False) -> pd.DataFrame:
        """The observation stack belonging to one member, by name or by index.

        The per-point view of uncertainty: with ``opt_chance_points`` = "all" every population
        member carries its own stack, so the risk shift differs across the population rather
        than being one shift applied everywhere.

        Columns are the CONSTRAINTS, not every observation - narrower than :meth:`stack_oe`,
        which also spans the objective. Rows match :meth:`stack_oe`, since it is the same
        stack evaluated at a different point in decision-variable space.
        """
        names = self.member_stacks()
        if len(names) == 0:
            raise RuntimeError(
                "there are no per-member stacks - " +
                ("this run is risk neutral" if not self.use_chance else
                 "this run uses fosm" if self.use_fosm else
                 'opt_chance_points is not "all", so one stack serves every point'))
        if isinstance(member, str):
            if member not in names:
                raise KeyError(f"no stack for member '{member}'; have {names}")
            idx = names.index(member)
        else:
            idx = int(member)
            if not (0 <= idx < len(names)):
                raise IndexError(f"member stack {idx} out of range, have {len(names)}")
        return self._stack_df(MEMBER_STACK_EN + idx, lower, "obsnme")

    @property
    def chance_config(self) -> dict:
        """Every chance/stack option with its live value, plus the DERIVED state.

        Gathered because they are meaningless apart - stacks with risk at 0.5 do nothing, and
        a risk without stacks silently means FOSM. The derived keys are what the tool will
        actually do; the rest is what was asked for.
        """
        opts = ["opt_risk", "opt_use_robust", "opt_stack_size", "opt_par_stack", "opt_obs_stack",
                "opt_chance_points", "opt_chance_schedule", "opt_recalc_chance_every",
                "std_weights"]
        out = {}
        for k in opts:
            try:
                out[k] = self.get_option(k)
            except Exception:
                out[k] = None
        out["risk (effective)"] = self.risk
        out["use_chance"] = self.use_chance
        out["use_fosm"] = self.use_fosm
        out["use_robust"] = self.use_robust
        return out

    def set_risk(self, risk: float) -> None:
        """Set the risk level, live.

        0.5 is risk-neutral and disables chance entirely; above 0.5 is risk-averse for a
        less-than constraint. Values outside [0.001, 0.999] are accepted by the option and
        clamped by the tool - :attr:`risk` reports what will actually be used.

        Refused on a robust run: robust optimization does no risk shifting, so there is no
        risk to set.
        """
        if self.use_robust:
            raise RuntimeError(
                "this is a robust run (opt_use_robust): it evaluates each realization against "
                "its own paired parameter draw and does no risk shifting, so opt_risk has no "
                "effect and the tool refuses the combination")
        self.set_option("opt_risk", float(risk))

    def set_stack_size(self, n: int) -> None:
        """Set the parameter stack size, live. Non-zero turns FOSM off in favour of stacks."""
        self.set_option("opt_stack_size", int(n))

    def recalc_chance_every(self, n: int) -> None:
        """How often to re-run the stacks. Every iteration is expensive; never is stale."""
        self.set_option("opt_recalc_chance_every", int(n))


class Mou(_ChanceMixin, _Tool):
    """Multi-objective optimization under uncertainty. One 'iteration' is a generation.

    mou speaks of MEMBERS, DECISION VARIABLES and OBJECTIVES, not realizations, parameters and
    phi. The shared surface uses the latter because it is shared; the aliases here let mou code
    read like mou. They are aliases, not copies - ``dv_pop()`` is ``par_df()``.
    """
    _tool_id = TOOL_MOU
    _has_phi = False

    @staticmethod
    def _agent_exe():
        return "pestpp-mou"

    # -- vocabulary ------------------------------------------------------------------------

    def dv_pop(self, lower: bool = False) -> pd.DataFrame:
        """The decision-variable population. Alias of :meth:`par_df`."""
        return self.par_df(lower=lower)

    def obs_pop(self, lower: bool = False) -> pd.DataFrame:
        """The observation population. Alias of :meth:`obs_df`."""
        return self.obs_df(lower=lower)

    @property
    def population_size(self) -> int:
        """Alias of :attr:`n_reals`."""
        return self.n_reals

    @property
    def members(self) -> np.ndarray:
        """Alias of :attr:`real_names`."""
        return self.real_names

    @property
    def ppd_config(self) -> dict:
        """The options that define a probabilistic-dominance run, with their live values.

        Gathered in one place because they are meaningless apart: a beta without NSGA_PPD does
        nothing, and an infill size without an outer repository has nowhere to send its points.
        """
        self._require_live("ppd_config")
        keys = ["mou_env_selector", "mou_ppd_beta", "mou_fit_gamma", "mou_fit_epsilon",
                "mou_infill_size", "mou_max_archive_size", "mou_generator",
                "mou_outer_repo_obs_file", "mou_dv_population_file",
                "mou_obs_population_restart_file"]
        out = {}
        for k in keys:
            try:
                out[k] = self.get_option(k)
            except Exception:
                out[k] = None
        return out


class Sqp(_ChanceMixin, _Tool):
    """Sequential quadratic programming.

    sqp is a search along a trajectory, not an ensemble method - the ensemble is machinery for
    estimating a gradient. It has an OBJECTIVE FUNCTION and a CONSTRAINT VIOLATION, which is
    why ``phi`` is unavailable: phi has no room for the second number, and the two together are
    what say whether a run worked.
    """
    _tool_id = TOOL_SQP
    _has_phi = False

    @staticmethod
    def _agent_exe():
        return "pestpp-sqp"

    def dv_df(self, lower: bool = False) -> pd.DataFrame:
        """The decision-variable ensemble. Alias of :meth:`par_df`."""
        return self.par_df(lower=lower)



# ---- worker management --------------------------------------------------------------------

def _find_exe(name: str, explicit: str | None = None) -> str:
    if explicit is not None:
        return explicit
    exe = name + (".exe" if os.name == "nt" else "")
    found = shutil.which(exe)
    if found:
        return found
    here = os.path.dirname(os.path.abspath(__file__))
    plat = platform.platform().lower()
    sub = "win" if os.name == "nt" else ("mac" if ("darwin" in plat or "macos" in plat) else "linux")
    for cand in (os.path.join(os.path.dirname(here), "bin", sub, exe),
                 os.path.join(os.path.dirname(here), "bin", exe)):
        if os.path.exists(cand):
            return os.path.abspath(cand)
    raise FileNotFoundError(
        "could not find {0}; pass exe_path= or put it on PATH".format(exe))


def _release_session(workers, lib) -> None:
    """Stop this session's agents and hand the handle back. The counterpart to _start_workers.

    Module level so a finalizer can hold it without holding the session (see the note in
    _Tool.__init__), and it never raises: this runs from garbage collection and from the
    interpreter's atexit unwind, where a traceback would be attached to no particular line of
    anyone's code.
    """
    # Signal everything first, THEN collect. Agents shut down concurrently, so two passes cost
    # about one timeout in total rather than one per worker - which matters when a notebook
    # cell that asked for workers=8 is being cleaned up on rebind.
    for p in workers:
        try:
            p.terminate()
        except Exception:
            pass
    for p in workers:
        try:
            p.wait(timeout=2)
        except Exception:
            pass        # already gone, or refusing to; either way do not block a GC pass
    try:
        if lib is not None:
            lib.destroy()
    except Exception:
        pass


def _start_workers(workdir, n, port, worker_root, exe_path, exe_name, pst_file):
    worker_root = os.path.abspath(worker_root or (workdir + "_workers"))
    if os.path.exists(worker_root):
        shutil.rmtree(worker_root)
    os.makedirs(worker_root)
    exe = _find_exe(exe_name, exe_path)
    # the control file this session was opened with, NOT whichever .pst happens to sort
    # first - a working directory with more than one is normal after a restart, and picking
    # the wrong one starts every worker on a different problem than the master
    pst = os.path.basename(pst_file)
    procs = []
    for i in range(n):
        d = os.path.join(worker_root, "worker_{0}".format(i))
        shutil.copytree(workdir, d)
        log = open(os.path.join(d, "worker.log"), "w")
        procs.append(subprocess.Popen([exe, pst, "/h", "localhost:{0}".format(port)],
                                      cwd=d, stdout=log, stderr=subprocess.STDOUT))
    return procs


# ---- the whole job in one call ------------------------------------------------------------

def _run_tool(cls, pst_file, workdir=".", workers=0, noptmax=None, **kwargs):
    opts = dict(kwargs)
    if noptmax is not None:
        opts["noptmax"] = noptmax
    with cls.from_pst(pst_file, workdir=workdir, workers=workers, **opts) as tool:
        tool.initialize()
        steps = [s for s in tool.iterations()]
        out = {"par": tool.par_df(), "obs": tool.obs_df(),
               "steps": pd.DataFrame([s.__dict__ for s in steps])}
        tool.finalize()
    return out


def run_ies(pst_file, workdir=".", workers=0, noptmax=None, **kwargs) -> dict:
    """Run pestpp-ies start to finish. Returns {'par', 'obs', 'steps'} as DataFrames."""
    return _run_tool(Ies, pst_file, workdir, workers, noptmax, **kwargs)


def run_da(pst_file, workdir=".", workers=0, noptmax=None, **kwargs) -> dict:
    """Run pestpp-da (single cycle) start to finish."""
    return _run_tool(Da, pst_file, workdir, workers, noptmax, **kwargs)


def run_mou(pst_file, workdir=".", workers=0, noptmax=None, **kwargs) -> dict:
    """Run pestpp-mou start to finish."""
    return _run_tool(Mou, pst_file, workdir, workers, noptmax, **kwargs)


def run_sqp(pst_file, workdir=".", workers=0, noptmax=None, **kwargs) -> dict:
    """Run pestpp-sqp start to finish."""
    return _run_tool(Sqp, pst_file, workdir, workers, noptmax, **kwargs)

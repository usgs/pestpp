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
from contextlib import contextmanager, nullcontext
from dataclasses import dataclass

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pestpp_lib import (  # noqa: E402
    PestppLib, PestppError,
    PAR_EN, OBS_EN, NOISE_EN, WEIGHTS_EN,
    TOOL_IES, TOOL_DA, TOOL_MOU, TOOL_SQP,
    RM_SERIAL, RM_PANTHER, RM_EXTERNAL,
    PHI_MEAS, PHI_COMPOSITE, PHI_REGUL, PHI_ACTUAL, PHI_NOISE,
    WORKER_COMPLETED, WORKER_FAILED, WORKER_TIMED_OUT,
)

__all__ = [
    "Ies", "Da", "Mou", "Sqp", "IterationStep", "PestppError",
    "run_ies", "run_da", "run_mou", "run_sqp", "find_library",
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
    name = ("pestpp_capi.dll" if ("windows" in plat or os.name == "nt")
            else "libpestpp_capi.dylib" if "darwin" in plat or "macos" in plat
            else "libpestpp_capi.so")
    here = os.path.dirname(os.path.abspath(__file__))
    roots = [os.path.join(os.path.dirname(here), "build"), os.path.join(here, "..", "..", "build")]
    for root in roots:
        if not os.path.isdir(root):
            continue
        for dirpath, _, filenames in os.walk(root):
            if name in filenames:
                return os.path.abspath(os.path.join(dirpath, name))
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
def _redirect_c_stdout(path, flush=None):
    """Send the library's stdout to a file while this block runs.

    PEST++ writes to C++ ``cout``, which goes to file descriptor 1 directly -- Python's
    ``contextlib.redirect_stdout`` only rebinds ``sys.stdout`` and does nothing about it. So
    this redirects at the fd level, which is the only thing that works.

    Captured to a file rather than discarded: the .rec file does not carry everything the
    console does, and swallowing a run manager's complaints outright would be worse than the
    noise. Only fd 1 is touched -- fd 2 is left alone so python tracebacks still surface.

    ``flush`` is called before the descriptor is restored. On windows the library links the
    static CRT and so buffers its output privately; without flushing it, the file gets the
    model's output (child processes inherit the descriptor) but not the library's own, which
    escapes to the console afterwards instead.
    """
    sys.stdout.flush()
    saved = os.dup(1)
    sink = open(path, "a")
    try:
        os.dup2(sink.fileno(), 1)
        yield
    finally:
        if flush is not None:
            try:
                flush()
            except Exception:
                pass          # a failed flush must not mask whatever the block raised
        sys.stdout.flush()
        os.dup2(saved, 1)
        os.close(saved)
        sink.close()


# ---- what an iteration reports ------------------------------------------------------------

@dataclass
class IterationStep:
    """One step of the algorithm, as yielded by :meth:`_Tool.iterations`."""
    iter: int
    n_reals: int
    retried: bool          # the step was rejected and the algorithm wants another attempt
    phi_mean: float = float("nan")
    phi_std: float = float("nan")
    phi_min: float = float("nan")
    phi_max: float = float("nan")

    def __repr__(self):
        return ("IterationStep(iter={0}, n_reals={1}, phi_mean={2:.6g}{3})"
                .format(self.iter, self.n_reals, self.phi_mean, ", RETRY" if self.retried else ""))


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

    def __init__(self, lib: PestppLib, workdir: str, workers=None, quiet: bool = True):
        self._lib = lib
        self._workdir = workdir
        self._workers = workers or []
        self._initialized = False
        self._queued = 0
        self._quiet = quiet
        self._log = os.path.join(workdir, "pestpp.stdout.log")

    def _q(self):
        """Capture the library's console output for the duration of a call."""
        return (_redirect_c_stdout(self._log, flush=self._lib.flush_output)
                if self._quiet else nullcontext())

    @property
    def log_file(self) -> str:
        """Where the library's console output went, when quiet."""
        return self._log

    # -- construction ----------------------------------------------------------------------

    @classmethod
    def from_pst(cls, pst_file: str, workdir: str = ".", workers: int = 0,
                 port: int | str = 4004, lib_path: str | None = None,
                 worker_root: str | None = None, exe_path: str | None = None,
                 run_manager: int | None = None, quiet: bool = True, **options):
        """Open a session on a control file.

        ``workers`` greater than zero starts a PANTHER master and that many worker processes,
        each in its own copy of ``workdir`` under ``worker_root``. Leave it at zero for the
        serial run manager, which is simpler and right for small problems.

        Extra keyword arguments are set as pest++ options before initialization, so
        ``Ies.from_pst("pest.pst", ies_num_reals=50)`` does what it looks like.

        ``quiet`` (the default) captures the library's console output to
        ``<workdir>/pestpp.stdout.log`` instead of letting it flood the session. Pass
        ``quiet=False`` to watch it live.
        """
        workdir = os.path.abspath(workdir)
        parallel = workers > 0
        if run_manager is None:
            run_manager = RM_PANTHER if parallel else RM_SERIAL
        log = os.path.join(workdir, "pestpp.stdout.log")
        with (_redirect_c_stdout(log) if quiet else nullcontext()):
            lib = PestppLib(find_library(lib_path), cls._tool_id, pst_file, workdir,
                            port=str(port) if run_manager == RM_PANTHER else None,
                            run_manager=run_manager)
            if quiet:
                lib.flush_output()          # before the redirect unwinds
        for key, value in options.items():
            lib.set_option(key, value)

        procs = []
        if parallel:
            # started before any blocking call: the master accepts workers whenever, but this
            # process is single-threaded, so anything started after initialize() would never
            # be reached
            procs = _start_workers(workdir, workers, port, worker_root, exe_path, cls._agent_exe())
        return cls(lib, workdir, procs, quiet=quiet)

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

    def solve(self) -> IterationStep:
        """One iteration (a generation, for mou)."""
        from pestpp_lib import PESTPP_RETRY
        with self._q():
            status = self._lib.solve_iteration()
        return self._step(retried=(status == PESTPP_RETRY))

    @property
    def noptmax(self) -> int:
        """Maximum iterations from the control data. -1 and 0 have their usual pest meanings."""
        try:
            return int(self._lib.get_option("NOPTMAX"))
        except (ValueError, PestppError):
            return -1

    def iterations(self, max_iter: int | None = None):
        """Yield an :class:`IterationStep` per iteration.

        Stops at ``noptmax``, or when the tool's own convergence test fires, whichever comes
        first. Both matter: should_terminate() only reports the phi-based criteria, so a loop
        driven by it alone would ignore noptmax entirely and run until phi stopped improving -
        which is not what a control file saying noptmax=3 asked for.

        ``max_iter`` overrides noptmax, for a caller that wants a few steps and a look around.
        """
        limit = max_iter if max_iter is not None else self.noptmax
        n = 0
        while not self.should_terminate:
            if limit is not None and limit >= 0 and n >= limit:
                break
            yield self.solve()
            n += 1

    def finalize(self) -> None:
        with self._q():
            self._lib.finalize()

    def close(self) -> None:
        """Release the handle and stop any workers this session started."""
        for p in self._workers:
            try:
                p.terminate()
            except Exception:
                pass
        self._workers = []
        if self._lib is not None:
            self._lib.destroy()
            self._lib = None

    def __enter__(self):
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
    def par_names(self) -> list:
        return self._lib.get_ensemble_col_names(PAR_EN)

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

    def update_phi(self) -> None:
        """Recompute phi from the current ensembles and weights.

        Phi is cached - the algorithm refreshes it at its own points. Call this after writing
        an ensemble or weights view if you want the phi accessors to reflect the change.
        """
        with self._q():
            self._lib.update_phi()

    def weights_df(self, lower: bool = False) -> pd.DataFrame:
        """The weights ENSEMBLE - one weight per observation per realization, as a copy."""
        arr = self._lib.get_ensemble_view(WEIGHTS_EN)
        return _maybe_lower(pd.DataFrame(
            arr.copy(),
            index=self._lib.get_ensemble_row_names(WEIGHTS_EN),
            columns=self._lib.get_ensemble_col_names(WEIGHTS_EN)), lower)

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

    # -- data ------------------------------------------------------------------------------

    def par_df(self, lower: bool = False) -> pd.DataFrame:
        """Every control-file parameter in CTL space, as a copy. Ready for pyemu."""
        vals, rows, cols = self._lib.get_par_snapshot()
        return _maybe_lower(pd.DataFrame(vals, index=rows, columns=cols), lower)

    def set_par_df(self, df: pd.DataFrame) -> None:
        """Push a parameter DataFrame back. Matched by name, so order does not matter."""
        self._lib.set_par_snapshot(df.values, _up_all(df.index), _up_all(df.columns))

    def obs_df(self, lower: bool = False) -> pd.DataFrame:
        """The observation ensemble, as a copy."""
        arr = self._lib.get_ensemble_view(OBS_EN)
        return _maybe_lower(pd.DataFrame(
            arr.copy(),
            index=self._lib.get_ensemble_row_names(OBS_EN),
            columns=self._lib.get_ensemble_col_names(OBS_EN)), lower)

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
        arr = self._lib.get_ensemble_view(which)
        try:
            yield arr
        finally:
            # drop our reference so the caller's name is the only one left; a caller who
            # stashed it elsewhere is on their own, which is what the docstring says
            del arr

    # -- membership ------------------------------------------------------------------------

    def drop_realizations(self, names) -> None:
        """Drop these from every coupled ensemble at once. Names are case-insensitive."""
        self._lib.drop_realizations(_up_all(names))

    def keep_realizations(self, names) -> None:
        """Keep only these, in this order. Also a reorder. Names are case-insensitive."""
        self._lib.keep_realizations(_up_all(names))

    # -- running things yourself -----------------------------------------------------------

    def queue_runs(self) -> int:
        """Queue the current ensemble with the run manager. Returns how many runs were queued."""
        with self._q():
            self._queued = self._lib.queue_runs()
        return self._queued

    def run(self, slice_seconds: float = 0.05, callback=None) -> None:
        """Drive the queued runs to completion.

        ``callback`` is called with this tool after each slice, which is where progress
        reporting or a cancel decision goes. PANTHER only - serial and external finish the
        whole batch in the first slice.
        """
        with self._q():
            self._lib.begin_batch()
            try:
                while not self._lib.run_slice(slice_seconds):
                    if callback is not None:
                        callback(self)
            finally:
                self._lib.end_batch()

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

    def _step(self, retried: bool) -> IterationStep:
        step = IterationStep(iter=self.iteration, n_reals=self.n_reals, retried=retried)
        if self._has_phi:
            try:
                p = self.get_phi()
                step.phi_mean, step.phi_std = p["mean"], p["std"]
                step.phi_min, step.phi_max = p["min"], p["max"]
            except PestppError:
                pass          # phi is not meaningful yet; leave the NaNs
        return step


class Ies(_Tool):
    """Iterative ensemble smoother."""
    _tool_id = TOOL_IES


class Da(_Tool):
    """Data assimilation. A SINGLE cycle - multi-cycle is not exposed yet."""
    _tool_id = TOOL_DA

    @staticmethod
    def _agent_exe():
        return "pestpp-da"


class Mou(_Tool):
    """Multi-objective optimization under uncertainty. One 'iteration' is a generation."""
    _tool_id = TOOL_MOU
    _has_phi = False

    @staticmethod
    def _agent_exe():
        return "pestpp-mou"


class Sqp(_Tool):
    """Sequential quadratic programming."""
    _tool_id = TOOL_SQP
    _has_phi = False

    @staticmethod
    def _agent_exe():
        return "pestpp-sqp"


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


def _start_workers(workdir, n, port, worker_root, exe_path, exe_name):
    worker_root = os.path.abspath(worker_root or (workdir + "_workers"))
    if os.path.exists(worker_root):
        shutil.rmtree(worker_root)
    os.makedirs(worker_root)
    exe = _find_exe(exe_name, exe_path)
    pst = [f for f in os.listdir(workdir) if f.lower().endswith(".pst")][0]
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

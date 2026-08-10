"""Thin ctypes layer over the PEST++ C ABI.

This mirrors ``pestpp-api.h`` one-to-one on purpose. Method names match the C symbols so
the two can be read side by side, and nothing here is made Pythonic -- no properties, no
DataFrames, no cleverness. That is the same split MODFLOW 6 uses (``xmipy`` thin,
``modflowapi`` opinionated) and pypestutils uses (``pestutilslib.py`` thin, ``helpers.py``
opinionated), and it exists so this file stays auditable against the header.

The ergonomic layer belongs on top of this, not inside it.
"""

from __future__ import annotations

import os
import weakref
from ctypes import (
    CDLL,
    c_ubyte,
    CFUNCTYPE,
    POINTER,
    Structure,
    byref,
    c_char,
    c_char_p,
    c_double,
    c_int,
    c_void_p,
    create_string_buffer,
    sizeof,
)

import numpy as np

# ---- enums, mirroring the header ----------------------------------------------------

# Banded by sign, exactly as the header defines them: negative is a successful call that has
# an algorithmic outcome worth reporting, zero is plain success, positive is failure. Codes
# get added inside the bands, so _check tests the sign rather than listing values.
PESTPP_RETRY = -1

PESTPP_OK = 0

PESTPP_ERROR = 1
PESTPP_INVALID_HANDLE = 2
PESTPP_INVALID_ARGUMENT = 4
PESTPP_BUFFER_TOO_SMALL = 5
PESTPP_NOT_SUPPORTED = 6
PESTPP_INVALID_STATE = 7

TOOL_IES, TOOL_DA, TOOL_MOU, TOOL_SQP = 0, 1, 2, 3
#: pestpp-glm. Not an ensemble method: it carries one parameter set through a Jacobian and an
#: upgrade, so the ensemble and phi-over-realizations calls refuse rather than hand back a
#: one-row ensemble. initialize / advance / finalize work as they do for the others.
TOOL_GLM = 4
PAR_EN, OBS_EN, NOISE_EN, WEIGHTS_EN = 0, 1, 2, 3
#: the chance stacks, mou and sqp only. Empty - not an error - on a fosm or risk-neutral run;
#: get_stack_status() is how you tell those apart from "not drawn yet"
STACK_PAR_EN, STACK_OBS_EN, NESTED_PAR_EN = 4, 5, 6
#: candidate ensembles of a deferred solve are ordinary ensemble ids, CANDIDATE_EN + i, so
#: views, names and snapshots all work on them unchanged
CANDIDATE_EN = 1000
#: per-member chance stacks, MEMBER_STACK_EN + i, likewise ordinary ensemble ids. Only when
#: opt_chance_points is "all" - otherwise one stack serves every point and the count is 0
MEMBER_STACK_EN = 2000
TSTAT = {0: "ctl", 1: "num", 2: "model"}

RUN_QUEUED, RUN_RUNNING, RUN_COMPLETED, RUN_FAILED, RUN_TIMED_OUT, RUN_CANCELLED = range(6)
RUN_STATUS = {
    RUN_QUEUED: "queued", RUN_RUNNING: "running", RUN_COMPLETED: "completed",
    RUN_FAILED: "failed", RUN_TIMED_OUT: "timed_out", RUN_CANCELLED: "cancelled",
}
WORKER_COMPLETED, WORKER_FAILED, WORKER_TIMED_OUT = 0, 1, 2
PHI_MEAS, PHI_COMPOSITE, PHI_REGUL, PHI_ACTUAL, PHI_NOISE = range(5)

RM_SERIAL, RM_PANTHER, RM_EXTERNAL = 0, 1, 2
RUN_MANAGER = {RM_SERIAL: "serial", RM_PANTHER: "panther", RM_EXTERNAL: "external"}


#: what a run observer asks for next; mirrors pestpp_run_action
PESTPP_RUN_CONTINUE, PESTPP_RUN_STOP_BATCH = 0, 1


class RunProgress(Structure):
    """Mirrors pestpp_run_progress. struct_size is the library's, not ours to set."""
    _fields_ = [("struct_size", c_int),
                ("n_total", c_int),
                ("n_completed", c_int),
                ("n_failed", c_int),
                ("n_timed_out", c_int),
                ("n_running", c_int),
                ("run_id", c_int),
                ("elapsed_sec", c_double)]


_RUN_OBSERVER_FN = CFUNCTYPE(c_int, POINTER(RunProgress), c_void_p)


class CreateOptions(Structure):
    """Mirrors pestpp_create_options. Field order and types must match the header exactly.

    struct_size is set for you in PestppLib.__init__; it is how the library tells which
    version of the struct it was handed, so fields added later stay source-compatible.
    """
    from ctypes import c_void_p as _vp  # noqa: F401  (kept local; see _fields_ below)
    _fields_ = [
        ("struct_size", c_int),
        ("tool", c_int),
        ("ctl_file", c_char_p),
        ("working_dir", c_char_p),
        ("run_manager", c_int),
        ("panther_port", c_char_p),
        # reserved must match the header exactly: it is what keeps a future field from hiding
        # in tail padding, where sizeof() would not change and struct_size would stop
        # distinguishing versions. ctypes zeroes it for us.
        ("reserved", _vp * 4),
    ]


class _Unset:
    """Sentinel for "no default given", since None is a legitimate default."""
    def __repr__(self):
        return "<unset>"


_UNSET = _Unset()


class PestppError(Exception):
    """Raised when a C ABI call returns a non-zero status.

    Its own type rather than RuntimeError, so callers can catch precisely.
    """


def format_option_value(value) -> str:
    """Render a python value the way the option parser reads it back.

    Options cross the C ABI as strings, but a caller should not have to know the spelling the
    parser wants - and the failure when they get it wrong is not obvious. This used to be a
    bare ``str()``, which is right for scalars by luck and wrong for sequences by construction:
    the VECTOR options (``ies_lambda_mults``, ``ies_n_iter_reinflate``, the reinflation
    schedule) are tokenized on commas, so ``[0.1, 1.0]`` arrived as the literal text
    ``"[0.1, 1.0]"`` and came back as "invalid value for option", naming the option rather
    than the brackets.

    Handled: str (verbatim), bool, int, float, numpy scalars, os.PathLike, and any list, tuple,
    range or array-like of those, joined with commas. Anything else raises TypeError here
    rather than being stringified into something the parser will reject later - including sets,
    because these are ORDERED schedules and a set does not have an order to offer.
    """
    if isinstance(value, str):
        return value
    # before the int branch: bool is a subclass of int, and "1"/"0" would be accepted but
    # reads back as a number, which is not what the caller wrote
    if isinstance(value, (bool, np.bool_)):
        return "true" if value else "false"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        return str(float(value))
    if isinstance(value, os.PathLike):
        return os.fspath(value)
    if isinstance(value, (set, frozenset)):
        raise TypeError(
            "option values cannot be given as a {0}: the vector options are ordered "
            "schedules and a set has no order. Pass a list or tuple instead".format(
                type(value).__name__))
    # ndarray and pandas Series both answer tolist(); range and tuple are already sequences.
    # tolist() on a 0-d array hands back a plain scalar, so route that back through here
    # rather than falling to the TypeError below
    seq = value.tolist() if hasattr(value, "tolist") else value
    if hasattr(value, "tolist") and not isinstance(seq, (list, tuple, range)):
        return format_option_value(seq)
    if isinstance(seq, (list, tuple, range)):
        if len(seq) == 0:
            raise ValueError("option values cannot be an empty sequence")
        return ",".join(format_option_value(v) for v in seq)
    raise TypeError(
        "cannot use {0} as an option value; pass a str, bool, int, float, path or a "
        "sequence of those".format(type(value).__name__))


def _destroy_handle(lib, handle) -> None:
    """Hand one session back to the library. Module level so a finalizer can hold it.

    Never raises. This runs from garbage collection and from the interpreter's atexit unwind,
    where an exception would surface as noise on stderr attached to no particular line of the
    caller's code - and where there is nothing useful to do about it anyway, the process being
    on its way out. A double free is already impossible: finalize() calls this at most once,
    and pestpp_destroy() answers PESTPP_INVALID_HANDLE rather than crashing besides.
    """
    try:
        if handle:
            lib.pestpp_destroy(handle)
    except Exception:
        pass


class PestppLib:
    """One loaded shared library plus one session handle."""

    def __init__(self, lib_path: str, tool: int, ctl_file: str, working_dir: str = ".",
                 port: str | None = None, run_manager: int | None = None):
        """run_manager is RM_SERIAL / RM_PANTHER / RM_EXTERNAL.

        Left as None it is inferred: a port means PANTHER, otherwise serial. Only PANTHER can
        be observed or interrupted mid-batch -- see supports_live_control().
        """
        if not os.path.exists(lib_path):
            raise FileNotFoundError(lib_path)
        self.lib = CDLL(lib_path)
        self._prototype()
        # Read the name width from the library rather than hardcoding it -- the same trick
        # xmipy and pypestutils use, and it removes a whole class of drift.
        self.name_len = c_int.in_dll(self.lib, "PESTPP_NAME_LEN").value

        from ctypes import c_void_p

        # holds the ctypes thunk for a registered run observer; C keeps a raw pointer to it,
        # so it must stay referenced for as long as the session lives
        self._observer_thunk = None
        self.handle = c_void_p()
        if run_manager is None:
            run_manager = RM_PANTHER if port is not None else RM_SERIAL
        opts = CreateOptions()
        opts.struct_size = sizeof(CreateOptions)
        opts.tool = tool
        opts.ctl_file = ctl_file.encode()
        opts.working_dir = working_dir.encode()
        opts.run_manager = run_manager
        opts.panther_port = None if port is None else str(port).encode()
        # keep the encoded strings alive: ctypes does not own what c_char_p points at
        self._opts = opts
        status = self.lib.pestpp_create(byref(opts), byref(self.handle))
        if status != PESTPP_OK:
            raise PestppError(self.last_global_error())

        # Cleanup backstop. `with` is a convenience here, not the contract - nothing gates the
        # methods on having entered a context - so a caller working flat (a notebook, where the
        # object has to outlive a cell) must not be punished for it. Without this, dropping the
        # last reference leaks the session: the run manager, RunStorage and the FileManager's
        # open output files all stay alive in the library for the rest of the process.
        #
        # weakref.finalize rather than __del__, deliberately. __del__ runs during interpreter
        # shutdown when module globals may already be torn down, so `self.lib.pestpp_destroy`
        # can fail confusingly at the worst possible moment; finalize keeps its own strong
        # references to exactly what the callback needs and is invoked in a defined order.
        #
        # The first argument is watched WEAKLY, so naming self there is right and costs
        # nothing. What must never appear is self among the callback ARGUMENTS (or a bound
        # method of self, which smuggles one in): finalize holds those strongly, so that would
        # keep the object alive forever - precisely the leak this exists to close. Hence a
        # module-level function taking the handle, not a method.
        self._finalizer = weakref.finalize(self, _destroy_handle, self.lib, self.handle)

    # -- plumbing ---------------------------------------------------------------------

    def _prototype(self) -> None:
        """Declare argtypes/restype for every function.

        Without this ctypes silently truncates pointers on 64-bit, which fails in ways that
        look like data corruption rather than a type error.
        """
        from ctypes import c_void_p

        lib = self.lib
        lib.pestpp_create.argtypes = (POINTER(CreateOptions), POINTER(c_void_p))
        lib.pestpp_create.restype = c_int
        lib.pestpp_get_run_manager.argtypes = (c_void_p, POINTER(c_int))
        lib.pestpp_get_run_manager.restype = c_int
        lib.pestpp_destroy.argtypes = (c_void_p,)
        lib.pestpp_destroy.restype = c_int
        lib.pestpp_last_error.argtypes = (c_void_p,)
        lib.pestpp_last_error.restype = c_char_p
        lib.pestpp_get_last_error.argtypes = (c_void_p, c_char_p, c_int, POINTER(c_int))
        lib.pestpp_get_last_error.restype = c_int
        lib.pestpp_last_global_error.argtypes = ()
        lib.pestpp_last_global_error.restype = c_char_p
        lib.pestpp_last_create_error.argtypes = ()
        lib.pestpp_last_create_error.restype = c_char_p
        lib.pestpp_get_fatal_error.argtypes = ()
        lib.pestpp_get_fatal_error.restype = c_char_p
        lib.pestpp_clear_fatal_error.argtypes = ()
        lib.pestpp_clear_fatal_error.restype = c_int
        lib.pestpp_get_version.argtypes = (c_char_p, c_int, POINTER(c_int))
        lib.pestpp_get_version.restype = c_int
        lib.pestpp_get_api_version.argtypes = (POINTER(c_int), POINTER(c_int), POINTER(c_int))
        lib.pestpp_get_api_version.restype = c_int
        lib.pestpp_flush_output.argtypes = ()
        lib.pestpp_flush_output.restype = c_int
        lib.pestpp_redirect_output.argtypes = (c_char_p, POINTER(c_int))
        lib.pestpp_redirect_output.restype = c_int
        lib.pestpp_restore_output.argtypes = (c_int,)
        lib.pestpp_restore_output.restype = c_int
        lib.pestpp_get_redirect_depth.argtypes = (POINTER(c_int),)
        lib.pestpp_get_redirect_depth.restype = c_int

        lib.pestpp_initialize_prepare.argtypes = (c_void_p, POINTER(c_int))
        lib.pestpp_initialize_prepare.restype = c_int
        lib.pestpp_initialize_finish.argtypes = (c_void_p,)
        lib.pestpp_initialize_finish.restype = c_int

        for name in ("pestpp_initialize", "pestpp_solve_iteration", "pestpp_finalize"):
            getattr(lib, name).argtypes = (c_void_p,)
            getattr(lib, name).restype = c_int

        lib.pestpp_solve_prepare.argtypes = (c_void_p, POINTER(c_int))
        lib.pestpp_solve_prepare.restype = c_int
        lib.pestpp_solve_finish.argtypes = (c_void_p, c_int, POINTER(c_int))
        lib.pestpp_solve_finish.restype = c_int
        lib.pestpp_get_candidate_count.argtypes = (c_void_p, POINTER(c_int))
        lib.pestpp_get_candidate_count.restype = c_int
        lib.pestpp_get_run_manager_settings.argtypes = (c_void_p, POINTER(c_int),
                                                        POINTER(c_double), POINTER(c_double),
                                                        POINTER(c_double))
        lib.pestpp_get_run_manager_settings.restype = c_int
        lib.pestpp_get_stack_status.argtypes = (c_void_p, POINTER(c_int), POINTER(c_int),
                                                POINTER(c_int), POINTER(c_double), POINTER(c_int))
        lib.pestpp_get_stack_status.restype = c_int
        lib.pestpp_get_member_stack_count.argtypes = (c_void_p, POINTER(c_int))
        lib.pestpp_get_member_stack_count.restype = c_int
        lib.pestpp_get_member_stack_names.argtypes = (c_void_p, c_char_p, c_int, POINTER(c_int))
        lib.pestpp_get_member_stack_names.restype = c_int
        lib.pestpp_get_candidate_info.argtypes = (
            c_void_p, c_int, POINTER(c_double), POINTER(c_double))
        lib.pestpp_get_candidate_info.restype = c_int
        lib.pestpp_queue_runs_subset.argtypes = (c_void_p, c_char_p, c_int, POINTER(c_int))
        lib.pestpp_queue_runs_subset.restype = c_int

        for name in ("pestpp_get_iteration", "pestpp_should_terminate",
                     "pestpp_get_par_transform_status"):
            getattr(lib, name).argtypes = (c_void_p, POINTER(c_int))
            getattr(lib, name).restype = c_int

        lib.pestpp_get_ensemble_view.argtypes = (
            c_void_p, c_int, POINTER(POINTER(c_double)), POINTER(c_int), POINTER(c_int),
            POINTER(c_int))
        lib.pestpp_get_ensemble_view.restype = c_int
        lib.pestpp_view_is_valid.argtypes = (c_void_p, c_int, POINTER(c_int))
        lib.pestpp_view_is_valid.restype = c_int
        lib.pestpp_release_view.argtypes = (c_void_p, c_int)
        lib.pestpp_release_view.restype = c_int

        for name in ("pestpp_get_ensemble_row_names", "pestpp_get_ensemble_col_names"):
            getattr(lib, name).argtypes = (c_void_p, c_int, c_char_p, c_int, POINTER(c_int))
            getattr(lib, name).restype = c_int

        lib.pestpp_set_option.argtypes = (c_void_p, c_char_p, c_char_p)
        lib.pestpp_set_option.restype = c_int
        lib.pestpp_get_option.argtypes = (
            c_void_p, c_char_p, c_char_p, c_int, POINTER(c_int), POINTER(c_int))
        lib.pestpp_get_option.restype = c_int

        # -- run management --
        lib.pestpp_supports_live_control.argtypes = (c_void_p, POINTER(c_int))
        lib.pestpp_supports_live_control.restype = c_int
        lib.pestpp_begin_batch.argtypes = (c_void_p,)
        lib.pestpp_begin_batch.restype = c_int
        lib.pestpp_run_slice.argtypes = (c_void_p, c_double, POINTER(c_int))
        lib.pestpp_run_slice.restype = c_int
        lib.pestpp_end_batch.argtypes = (c_void_p,)
        lib.pestpp_end_batch.restype = c_int
        lib.pestpp_get_run_partial_info.argtypes = (
            c_void_p, c_int, POINTER(c_int), POINTER(c_int), POINTER(c_int))
        lib.pestpp_get_run_partial_info.restype = c_int
        lib.pestpp_get_run_count.argtypes = (c_void_p, POINTER(c_int))
        lib.pestpp_get_run_count.restype = c_int
        lib.pestpp_get_run_info.argtypes = (
            c_void_p, c_int, POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int))
        lib.pestpp_get_run_info.restype = c_int
        lib.pestpp_get_run_values.argtypes = (
            c_void_p, c_int, POINTER(c_double), c_int, POINTER(c_double), c_int,
            POINTER(c_ubyte), POINTER(c_int), POINTER(c_int))
        lib.pestpp_get_run_par_names.argtypes = (c_void_p, c_char_p, c_int, POINTER(c_int))
        lib.pestpp_get_run_par_names.restype = c_int
        lib.pestpp_get_run_obs_names.argtypes = (c_void_p, c_char_p, c_int, POINTER(c_int))
        lib.pestpp_get_run_obs_names.restype = c_int
        lib.pestpp_set_run_values.argtypes = (c_void_p, c_int, POINTER(c_double), c_int)
        lib.pestpp_set_run_values.restype = c_int
        lib.pestpp_set_run_failed.argtypes = (c_void_p, c_int)
        lib.pestpp_set_run_failed.restype = c_int
        lib.pestpp_get_run_values.restype = c_int
        lib.pestpp_request_partial_results.argtypes = (
            c_void_p, POINTER(c_int), c_int, POINTER(c_int))
        lib.pestpp_request_partial_results.restype = c_int
        lib.pestpp_set_run_observer.argtypes = (
            c_void_p, _RUN_OBSERVER_FN, c_void_p, c_double)
        lib.pestpp_set_run_observer.restype = c_int
        lib.pestpp_get_run_time_stats.argtypes = (
            c_void_p, POINTER(c_double), POINTER(c_int), POINTER(c_int), POINTER(c_int),
            POINTER(c_int), POINTER(c_int))
        lib.pestpp_get_run_time_stats.restype = c_int
        lib.pestpp_get_run_states.argtypes = (
            c_void_p, POINTER(c_int), c_int, POINTER(c_int), POINTER(c_int), POINTER(c_double),
            POINTER(c_int), c_char_p, c_int, POINTER(c_int))
        lib.pestpp_get_run_states.restype = c_int
        lib.pestpp_cancel_runs.argtypes = (c_void_p, POINTER(c_int), c_int, POINTER(c_int))
        lib.pestpp_cancel_runs.restype = c_int
        lib.pestpp_release_workers.argtypes = (c_void_p, POINTER(c_int), c_int, POINTER(c_int))
        lib.pestpp_release_workers.restype = c_int
        lib.pestpp_jacobian_prepare.argtypes = (c_void_p, c_int, POINTER(c_int))
        lib.pestpp_jacobian_prepare.restype = c_int
        lib.pestpp_jacobian_run.argtypes = (c_void_p,)
        lib.pestpp_jacobian_run.restype = c_int
        lib.pestpp_jacobian_process.argtypes = (c_void_p, POINTER(c_int))
        lib.pestpp_jacobian_process.restype = c_int
        lib.pestpp_get_jacobian.argtypes = (
            c_void_p, POINTER(c_double), c_int, c_int, POINTER(c_int), POINTER(c_int),
            c_char_p, c_char_p)
        lib.pestpp_get_jacobian.restype = c_int
        lib.pestpp_get_par_vector.argtypes = (
            c_void_p, c_int, POINTER(c_double), c_int, POINTER(c_int), c_char_p)
        lib.pestpp_get_par_vector.restype = c_int
        lib.pestpp_set_par_vector.argtypes = (c_void_p, POINTER(c_double), c_int, c_char_p)
        lib.pestpp_set_par_vector.restype = c_int
        lib.pestpp_get_worker_count.argtypes = (c_void_p, POINTER(c_int))
        lib.pestpp_get_worker_count.restype = c_int
        lib.pestpp_get_worker_state.argtypes = (
            c_void_p, c_int, c_char_p, c_int, c_char_p, c_int, POINTER(c_int),
            POINTER(c_double), POINTER(c_double), POINTER(c_int))
        lib.pestpp_get_worker_state.restype = c_int
        lib.pestpp_get_worker_run_history.argtypes = (
            c_void_p, c_int, c_int, POINTER(c_int), c_int, POINTER(c_int))
        lib.pestpp_get_worker_run_history.restype = c_int

        # -- queue/process, membership, snapshot --
        lib.pestpp_queue_runs.argtypes = (c_void_p, POINTER(c_int))
        lib.pestpp_queue_runs.restype = c_int
        lib.pestpp_process_runs.argtypes = (c_void_p, POINTER(c_int))
        lib.pestpp_process_runs.restype = c_int
        for name in ("pestpp_drop_realizations", "pestpp_keep_realizations"):
            getattr(lib, name).argtypes = (c_void_p, c_char_p, c_int)
            getattr(lib, name).restype = c_int
        lib.pestpp_get_phi_summary.argtypes = (
            c_void_p, c_int, POINTER(c_double), POINTER(c_double), POINTER(c_double),
            POINTER(c_double))
        lib.pestpp_get_phi_summary.restype = c_int
        lib.pestpp_get_phi_vector.argtypes = (
            c_void_p, c_int, POINTER(c_double), c_char_p, c_int, POINTER(c_int))
        lib.pestpp_get_phi_vector.restype = c_int
        lib.pestpp_get_phi_residuals.argtypes = (
            c_void_p, c_int, POINTER(c_double), c_int, c_int, POINTER(c_int), POINTER(c_int),
            c_char_p, c_char_p)
        lib.pestpp_get_phi_residuals.restype = c_int
        lib.pestpp_get_obs_groups.argtypes = (c_void_p, c_char_p, c_int, POINTER(c_int))
        lib.pestpp_get_obs_groups.restype = c_int
        lib.pestpp_get_obs_weights.argtypes = (c_void_p, POINTER(c_double), c_int,
                                               POINTER(c_int))
        lib.pestpp_get_obs_weights.restype = c_int
        lib.pestpp_set_obs_weights.argtypes = (c_void_p, c_char_p, POINTER(c_double), c_int)
        lib.pestpp_set_obs_weights.restype = c_int
        lib.pestpp_broadcast_weights.argtypes = (c_void_p,)
        lib.pestpp_broadcast_weights.restype = c_int
        lib.pestpp_reinflate_ensemble.argtypes = (c_void_p, c_double, c_int, c_int)
        lib.pestpp_reinflate_ensemble.restype = c_int
        lib.pestpp_update_phi.argtypes = (c_void_p,)
        lib.pestpp_update_phi.restype = c_int
        lib.pestpp_get_par_snapshot.argtypes = (
            c_void_p, POINTER(c_double), c_int, c_int, POINTER(c_int), POINTER(c_int),
            c_char_p, c_char_p)
        lib.pestpp_get_par_snapshot.restype = c_int
        lib.pestpp_set_par_snapshot.argtypes = (
            c_void_p, POINTER(c_double), c_int, c_int, c_char_p, c_char_p)
        lib.pestpp_set_par_snapshot.restype = c_int

    def _check(self, status: int, what: str) -> int:
        """Raise on failure; pass the negative band through, since it is an outcome.

        Tests the sign rather than enumerating codes: the header reserves the bands and will
        add to them, and a wrapper that listed values would mis-handle the first new one.
        """
        if status <= PESTPP_OK:
            return status
        # read the message before anything else touches the handle -- every entry point
        # clears it on the way in, so even a getter here would erase what we are reporting
        msg = self.get_last_error()
        raise PestppError("{0}: {1} (status {2})".format(what, msg or "unknown error", status))

    def get_last_error(self) -> str:
        """Why the most recent call on this handle failed; "" if it did not.

        Uses the copying form so the result is a python string the caller owns, rather than a
        pointer into storage the next call overwrites.
        """
        needed = c_int()
        if self.lib.pestpp_get_last_error(self.handle, None, 0, byref(needed)) != PESTPP_OK:
            return ""
        if needed.value <= 1:
            return ""
        buf = create_string_buffer(needed.value)
        if self.lib.pestpp_get_last_error(
                self.handle, buf, c_int(needed.value), byref(needed)) != PESTPP_OK:
            return ""
        return buf.value.decode(errors="replace")

    def _unpack_names(self, raw: bytes, count: int) -> list[str]:
        return [
            raw[i * self.name_len:(i + 1) * self.name_len].decode().strip()
            for i in range(count)
        ]

    # -- lifecycle --------------------------------------------------------------------

    def flush_output(self) -> None:
        """Flush the library's console output.

        Matters on windows, where the library has its own static CRT and so its own stdout
        buffer: without this, output captured by redirecting fd 1 misses everything the
        library itself printed.
        """
        self.lib.pestpp_flush_output()

    def redirect_output(self, path: str) -> int:
        """Send the library's console output to a file. Returns a token to restore with.

        The token is not a file descriptor -- the saved descriptor stays inside the library,
        so there is nothing here to close or to pass to the wrong call.
        """
        token = c_int()
        if self.lib.pestpp_redirect_output(str(path).encode(), byref(token)) != PESTPP_OK:
            raise PestppError("could not redirect output to {0}: {1}".format(
                path, self.last_global_error()))
        return token.value

    def restore_output(self, redirect_token: int) -> None:
        """Undo redirect_output(). Strictly innermost-first.

        Raises on failure rather than returning quietly: if this does not work the process's
        stdout is gone, and swallowing the status means nobody ever finds out why. Restoring
        out of order is a failure too -- see the header for why doing it anyway would leave
        stdout pointing somewhere permanently wrong.
        """
        if self.lib.pestpp_restore_output(c_int(redirect_token)) != PESTPP_OK:
            raise PestppError("could not restore output: {0}".format(self.last_global_error()))

    def get_redirect_depth(self) -> int:
        """How many output redirects are outstanding process-wide. 0 when not capturing."""
        depth = c_int()
        self.lib.pestpp_get_redirect_depth(byref(depth))
        return depth.value

    def last_global_error(self) -> str:
        """Why the most recent handle-less call failed -- create, redirect, restore, flush."""
        msg = self.lib.pestpp_last_global_error()
        return msg.decode(errors="replace") if msg else ""

    def get_fatal_error(self) -> str:
        """The process-wide latched failure, or "" when healthy. See clear_fatal_error()."""
        msg = self.lib.pestpp_get_fatal_error()
        return msg.decode(errors="replace") if msg else ""

    def clear_fatal_error(self) -> None:
        """Acknowledge the latched failure and resume.

        Only meaningful once the working directory has actually been put back -- the library
        cannot check that, so this is an assertion by the caller, not a repair.
        """
        self.lib.pestpp_clear_fatal_error()

    def get_version(self) -> str:
        """The pest++ release this library was built from, e.g. "5.2.28"."""
        needed = c_int()
        if self.lib.pestpp_get_version(None, 0, byref(needed)) != PESTPP_OK:
            return ""
        buf = create_string_buffer(needed.value)
        if self.lib.pestpp_get_version(buf, c_int(needed.value), byref(needed)) != PESTPP_OK:
            return ""
        return buf.value.decode(errors="replace")

    def get_api_version(self) -> tuple[int, int, int]:
        """The C ABI's own (major, minor, patch), which moves apart from the release above."""
        major, minor, patch = c_int(), c_int(), c_int()
        self.lib.pestpp_get_api_version(byref(major), byref(minor), byref(patch))
        return major.value, minor.value, patch.value

    def destroy(self) -> None:
        """Release the session. Idempotent, and safe to leave to the finalizer."""
        # Routed through the finalizer rather than duplicating the call, so explicit destroy(),
        # __exit__ and garbage collection are one code path that runs at most once.
        fin = getattr(self, "_finalizer", None)
        if fin is not None:
            fin()
        self.handle = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.destroy()

    # -- driving ----------------------------------------------------------------------

    def initialize(self) -> None:
        self._check(self.lib.pestpp_initialize(self.handle), "pestpp_initialize")

    def initialize_prepare(self) -> int:
        """Draw the ensembles without running them.

        Returns how many runs the caller must service before initialize_finish(). 0 means
        none -- a restart supplies results instead, or the tool initializes atomically.

        Between this call and queue_runs() is the only point at which the prior ensemble
        can be replaced with your own; see set_par_snapshot().
        """
        n = c_int()
        self._check(self.lib.pestpp_initialize_prepare(self.handle, byref(n)),
                    "pestpp_initialize_prepare")
        return n.value

    def initialize_finish(self) -> None:
        """Process the prior-ensemble results: phi, drop failures, write the .0. files."""
        self._check(self.lib.pestpp_initialize_finish(self.handle), "pestpp_initialize_finish")

    def solve_iteration(self) -> int:
        """Run one iteration. Returns PESTPP_OK or PESTPP_RETRY."""
        return self._check(self.lib.pestpp_solve_iteration(self.handle), "pestpp_solve_iteration")

    def solve_prepare(self) -> int:
        """Generate the candidates without running them. Returns the runs they imply."""
        n = c_int()
        self._check(self.lib.pestpp_solve_prepare(self.handle, byref(n)), "pestpp_solve_prepare")
        return n.value

    def solve_finish(self, defer_runs: bool = False) -> tuple[int, int]:
        """Continue after a batch. Returns (status, pending_runs)."""
        pending = c_int()
        st = self._check(
            self.lib.pestpp_solve_finish(self.handle, c_int(1 if defer_runs else 0),
                                         byref(pending)), "pestpp_solve_finish")
        return st, pending.value

    def get_candidate_count(self) -> int:
        n = c_int()
        self._check(self.lib.pestpp_get_candidate_count(self.handle, byref(n)),
                    "pestpp_get_candidate_count")
        return n.value

    def get_run_manager_settings(self, overdue: bool = True) -> dict:
        """What the LIVE run manager is using, not what the options say.

        The overdue policy is panther-only; pass overdue=False on a serial or external session.
        """
        mrf = c_int()
        rs, gf, gm = c_double(), c_double(), c_double()
        self._check(self.lib.pestpp_get_run_manager_settings(
            self.handle, byref(mrf),
            byref(rs) if overdue else None,
            byref(gf) if overdue else None,
            byref(gm) if overdue else None), "pestpp_get_run_manager_settings")
        out = {"max_run_fail": mrf.value}
        if overdue:
            out.update({"overdue_resched_fac": rs.value, "overdue_giveup_fac": gf.value,
                        "overdue_giveup_minutes": gm.value})
        return out

    def get_stack_status(self) -> dict:
        """How chance is being accounted for right now - mou and sqp only.

        Derived from the options at the moment it is asked, not latched at initialize(), so a
        risk set through set_option() since then is reflected. Keys: use_chance, use_fosm,
        use_robust, risk (clamped to [0.001, 0.999]), stack_size (the stack AS IT STANDS, 0
        before it is drawn).
        """
        chance, fosm, robust, size = c_int(), c_int(), c_int(), c_int()
        risk = c_double()
        self._check(self.lib.pestpp_get_stack_status(self.handle, byref(chance), byref(fosm),
                                                     byref(robust), byref(risk), byref(size)),
                    "pestpp_get_stack_status")
        return {"use_chance": bool(chance.value), "use_fosm": bool(fosm.value),
                "use_robust": bool(robust.value), "risk": risk.value,
                "stack_size": size.value}

    def get_member_stack_count(self) -> int:
        """How many per-member stacks exist - 0 unless opt_chance_points is "all"."""
        n = c_int()
        self._check(self.lib.pestpp_get_member_stack_count(self.handle, byref(n)),
                    "pestpp_get_member_stack_count")
        return n.value

    def get_member_stack_names(self) -> list[str]:
        """The members owning each per-member stack, in MEMBER_STACK_EN + i order."""
        return self._get_storage_names(self.lib.pestpp_get_member_stack_names)

    def get_candidate_info(self, idx: int) -> tuple[float, float]:
        """(inflation, backtrack) factors candidate `idx` was generated with."""
        inf, back = c_double(), c_double()
        self._check(self.lib.pestpp_get_candidate_info(self.handle, c_int(idx), byref(inf),
                                                       byref(back)), "pestpp_get_candidate_info")
        return inf.value, back.value

    def queue_runs_subset(self, names=None) -> int:
        """Queue the deferred solve's outstanding batch. None means the algorithm's choice."""
        n = c_int()
        if names:
            names = list(names)
            self._check(self.lib.pestpp_queue_runs_subset(
                self.handle, self._pack_names(names), len(names), byref(n)),
                "pestpp_queue_runs_subset")
        else:
            self._check(self.lib.pestpp_queue_runs_subset(self.handle, None, 0, byref(n)),
                        "pestpp_queue_runs_subset")
        return n.value

    def finalize(self) -> None:
        self._check(self.lib.pestpp_finalize(self.handle), "pestpp_finalize")

    def get_iteration(self) -> int:
        v = c_int()
        self._check(self.lib.pestpp_get_iteration(self.handle, byref(v)), "pestpp_get_iteration")
        return v.value

    def should_terminate(self) -> bool:
        v = c_int()
        self._check(self.lib.pestpp_should_terminate(self.handle, byref(v)),
                    "pestpp_should_terminate")
        return bool(v.value)

    # -- ensembles --------------------------------------------------------------------

    def get_ensemble_view(self, ensemble_id: int) -> tuple[np.ndarray, int]:
        """A numpy view onto the live buffer -- no copy. Returns (array, view_token).

        Eigen is column-major, hence order="F". The array is only valid while that
        ensemble's storage is unchanged; ask view_is_valid(token) rather than guessing, and
        release_view(token) when done.
        """
        ptr = POINTER(c_double)()
        nrow, ncol, token = c_int(), c_int(), c_int()
        self._check(
            self.lib.pestpp_get_ensemble_view(
                self.handle, c_int(ensemble_id), byref(ptr), byref(nrow), byref(ncol),
                byref(token)),
            "pestpp_get_ensemble_view")
        n = nrow.value * ncol.value
        if n == 0:
            # An empty ensemble has no buffer to borrow. Return a READ-ONLY array: a writable
            # empty one is a trap, because writes to it go nowhere and report no error.
            empty = np.empty((nrow.value, ncol.value), order="F")
            empty.flags.writeable = False
            return empty, token.value
        buf = np.ctypeslib.as_array(ptr, shape=(n,))
        return buf.reshape((nrow.value, ncol.value), order="F"), token.value

    def view_is_valid(self, view_token: int) -> bool:
        """Is the buffer that token was issued for still the ensemble's live storage?"""
        out = c_int()
        self._check(
            self.lib.pestpp_view_is_valid(self.handle, c_int(view_token), byref(out)),
            "pestpp_view_is_valid")
        return bool(out.value)

    def release_view(self, view_token: int) -> None:
        self._check(self.lib.pestpp_release_view(self.handle, c_int(view_token)),
                    "pestpp_release_view")

    def get_ensemble_row_names(self, ensemble_id: int) -> list[str]:
        return self._get_names(self.lib.pestpp_get_ensemble_row_names, ensemble_id)

    def get_ensemble_col_names(self, ensemble_id: int) -> list[str]:
        return self._get_names(self.lib.pestpp_get_ensemble_col_names, ensemble_id)

    def _get_names(self, fn, ensemble_id: int) -> list[str]:
        count = c_int()
        # query the count first, then size the buffer -- no library-allocated memory
        self._check(fn(self.handle, c_int(ensemble_id), None, 0, byref(count)), fn.__name__)
        if count.value == 0:
            return []
        buf = create_string_buffer(count.value * self.name_len)
        self._check(fn(self.handle, c_int(ensemble_id), buf, count.value * self.name_len,
                       byref(count)), fn.__name__)
        return self._unpack_names(buf.raw, count.value)

    def get_par_transform_status(self) -> str:
        v = c_int()
        self._check(self.lib.pestpp_get_par_transform_status(self.handle, byref(v)),
                    "pestpp_get_par_transform_status")
        return TSTAT.get(v.value, str(v.value))

    # -- options ----------------------------------------------------------------------

    def set_option(self, key: str, value) -> None:
        self._check(
            self.lib.pestpp_set_option(self.handle, key.encode(),
                                       format_option_value(value).encode()),
            f"pestpp_set_option({key})")

    def get_option(self, key: str, default=_UNSET):
        """One option's value.

        An unknown key raises, matching set_option -- unless a `default` is supplied, which is
        how you ask "does this library know about X?" without an exception. Note the two are
        genuinely different questions: an option set to the empty string returns "", while an
        option that does not exist returns the default.
        """
        needed, found = c_int(), c_int()
        probing = default is not _UNSET
        fp = byref(found) if probing else None
        self._check(
            self.lib.pestpp_get_option(self.handle, key.encode(), None, 0, byref(needed), fp),
            f"pestpp_get_option({key})")
        if probing and not found.value:
            return default
        buf = create_string_buffer(needed.value)
        self._check(self.lib.pestpp_get_option(self.handle, key.encode(), buf, needed.value,
                                               byref(needed), fp), f"pestpp_get_option({key})")
        return buf.value.decode()

    def has_option(self, key: str) -> bool:
        """Whether this library knows the option at all, regardless of its value."""
        found = c_int()
        self._check(
            self.lib.pestpp_get_option(self.handle, key.encode(), None, 0, None, byref(found)),
            f"pestpp_get_option({key})")
        return bool(found.value)

    # -- run management ----------------------------------------------------------------

    def get_run_manager(self) -> str:
        """Which run manager this handle got: 'serial', 'panther' or 'external'."""
        v = c_int()
        self._check(self.lib.pestpp_get_run_manager(self.handle, byref(v)),
                    "pestpp_get_run_manager")
        return RUN_MANAGER.get(v.value, str(v.value))

    def supports_live_control(self) -> bool:
        """False for the serial run manager, whose runs cannot be watched or cancelled."""
        v = c_int()
        self._check(self.lib.pestpp_supports_live_control(self.handle, byref(v)),
                    "pestpp_supports_live_control")
        return bool(v.value)

    def begin_batch(self) -> None:
        self._check(self.lib.pestpp_begin_batch(self.handle), "pestpp_begin_batch")

    def run_slice(self, max_seconds: float = 0.05) -> bool:
        """Advance the run manager for about max_seconds. Returns True when all runs are done."""
        done = c_int()
        self._check(self.lib.pestpp_run_slice(self.handle, c_double(max_seconds), byref(done)),
                    "pestpp_run_slice")
        return bool(done.value)

    def end_batch(self) -> None:
        self._check(self.lib.pestpp_end_batch(self.handle), "pestpp_end_batch")

    def set_run_observer(self, fn, min_interval_sec: float = 0.0) -> None:
        """Watch every batch this session runs, from inside the library.

        ``fn`` is called with a dict of counters and returns True to carry on or False to stop
        the batch. Pass None to stop observing.

        The thunk is stored on THIS OBJECT deliberately. A ctypes callback held only by a
        local is garbage-collected while C still has the pointer, and the crash lands at some
        unrelated later run - the single most common way to get this wrong.
        """
        if fn is None:
            # a NULL of the right function-pointer type: ctypes will not accept a bare None
            # where it has been told to expect a CFUNCTYPE
            self._check(self.lib.pestpp_set_run_observer(
                self.handle, _RUN_OBSERVER_FN(), None, c_double(0.0)),
                "pestpp_set_run_observer")
            self._observer_thunk = None      # only after C has stopped pointing at it
            return

        def _trampoline(progress_ptr, _user_data):
            # nothing may escape into C: an exception here would unwind through the run
            # manager, and on some ABIs through a C function pointer at all
            try:
                p = progress_ptr.contents
                out = fn({"n_total": p.n_total, "n_completed": p.n_completed,
                          "n_failed": p.n_failed, "n_timed_out": p.n_timed_out,
                          "n_running": p.n_running, "run_id": p.run_id,
                          "elapsed_sec": p.elapsed_sec})
                return PESTPP_RUN_CONTINUE if (out is None or out) else PESTPP_RUN_STOP_BATCH
            except Exception:
                return PESTPP_RUN_CONTINUE

        self._observer_thunk = _RUN_OBSERVER_FN(_trampoline)
        self._check(self.lib.pestpp_set_run_observer(
            self.handle, self._observer_thunk, None, c_double(float(min_interval_sec))),
            "pestpp_set_run_observer")

    def get_run_partial_info(self, run_id: int) -> dict:
        """Live: what a worker has reported for this run in the CURRENT batch."""
        has, rep, tot = c_int(), c_int(), c_int()
        self._check(self.lib.pestpp_get_run_partial_info(
            self.handle, c_int(run_id), byref(has), byref(rep), byref(tot)),
            "pestpp_get_run_partial_info")
        return {"has_partial": bool(has.value), "n_obs_reported": rep.value,
                "n_obs_total": tot.value}

    def get_run_count(self) -> int:
        n = c_int()
        self._check(self.lib.pestpp_get_run_count(self.handle, byref(n)),
                    "pestpp_get_run_count")
        return n.value

    def get_run_info(self, run_id: int) -> dict:
        """Stored: status, and how complete the VALUES are (0 none, 1 partial, 2 final)."""
        st, comp, rep, tot = c_int(), c_int(), c_int(), c_int()
        self._check(self.lib.pestpp_get_run_info(
            self.handle, c_int(run_id), byref(st), byref(comp), byref(rep), byref(tot)),
            "pestpp_get_run_info")
        return {"status": st.value, "completeness": comp.value,
                "n_obs_reported": rep.value, "n_obs_total": tot.value}

    def get_run_values(self, run_id: int):
        """Stored (pars, obs, obs_valid) for one run. obs_valid marks the REAL values."""
        npar, nobs = c_int(), c_int()
        self._check(self.lib.pestpp_get_run_values(
            self.handle, c_int(run_id), None, 0, None, 0, None, byref(npar), byref(nobs)),
            "pestpp_get_run_values")
        p = (c_double * npar.value)()
        o = (c_double * nobs.value)()
        v = (c_ubyte * nobs.value)()
        self._check(self.lib.pestpp_get_run_values(
            self.handle, c_int(run_id), p, npar.value, o, nobs.value, v,
            byref(npar), byref(nobs)), "pestpp_get_run_values")
        return (np.ctypeslib.as_array(p).copy(), np.ctypeslib.as_array(o).copy(),
                np.ctypeslib.as_array(v).copy().astype(bool))

    def get_run_par_names(self) -> list[str]:
        """Parameter names in RUN STORAGE order - what get_run_values() returns."""
        return self._get_storage_names(self.lib.pestpp_get_run_par_names)

    def get_run_obs_names(self) -> list[str]:
        """Observation names in RUN STORAGE order - what set_run_values() expects.

        Not the same as the ensemble's column order, which may differ and may be a subset.
        """
        return self._get_storage_names(self.lib.pestpp_get_run_obs_names)

    def _get_storage_names(self, fn) -> list[str]:
        n = c_int()
        self._check(fn(self.handle, None, 0, byref(n)), "get run storage names")
        buf = create_string_buffer(n.value * self.name_len)
        self._check(fn(self.handle, buf, n.value * self.name_len, byref(n)),
                    "get run storage names")
        return self._unpack_names(buf.raw, n.value)

    def set_run_values(self, run_id: int, obs) -> None:
        """Record this run's observations and mark it complete - service the run yourself.

        The mirror of get_run_values(): read the parameters, evaluate them however you like,
        write the results back. `obs` is in get_obs_names() order and must be the full set.
        """
        arr = np.ascontiguousarray(np.asarray(obs, dtype=np.float64))
        buf = (c_double * arr.size)(*arr.tolist())
        self._check(self.lib.pestpp_set_run_values(
            self.handle, c_int(run_id), buf, c_int(arr.size)), "pestpp_set_run_values")

    def set_run_failed(self, run_id: int) -> None:
        """Mark this run failed. Distinct from not writing it: this COUNTS as a failure."""
        self._check(self.lib.pestpp_set_run_failed(self.handle, c_int(run_id)),
                    "pestpp_set_run_failed")

    def request_partial_results(self, run_ids=None) -> int:
        """Ask the workers running these runs for whatever they have. Returns requests sent."""
        n_req = c_int()
        if run_ids:
            ids = list(run_ids)
            arr = (c_int * len(ids))(*ids)
            self._check(self.lib.pestpp_request_partial_results(
                self.handle, arr, len(ids), byref(n_req)), "pestpp_request_partial_results")
        else:
            self._check(self.lib.pestpp_request_partial_results(
                self.handle, None, 0, byref(n_req)), "pestpp_request_partial_results")
        return n_req.value

    def get_run_time_stats(self) -> dict:
        avg = c_double()
        cs = [c_int() for _ in range(5)]
        self._check(self.lib.pestpp_get_run_time_stats(
            self.handle, byref(avg), *[byref(c) for c in cs]), "pestpp_get_run_time_stats")
        keys = ("completed", "failed", "timed_out", "queued", "running")
        out = {"avg_run_sec": avg.value}
        out.update({k: c.value for k, c in zip(keys, cs)})
        return out

    def get_run_states(self, run_ids: list[int] | None = None) -> list[dict]:
        """One dict per run. Pass run_ids to ask about specific runs only."""
        want, n_want = None, 0
        if run_ids:
            want = (c_int * len(run_ids))(*run_ids)
            n_want = len(run_ids)
        n = c_int()
        # size first, with every column NULL
        self._check(self.lib.pestpp_get_run_states(
            self.handle, want, n_want, None, None, None, None, None, 0, byref(n)),
            "pestpp_get_run_states")
        if n.value == 0:
            return []
        ids = (c_int * n.value)()
        sts = (c_int * n.value)()
        elapsed = (c_double * n.value)()
        nfail = (c_int * n.value)()
        hosts = create_string_buffer(n.value * self.name_len)
        self._check(self.lib.pestpp_get_run_states(
            self.handle, want, n_want, ids, sts, elapsed, nfail, hosts, n.value, byref(n)),
            "pestpp_get_run_states")
        host_names = self._unpack_names(hosts.raw, n.value)
        return [
            {"run_id": ids[i], "status": RUN_STATUS.get(sts[i], sts[i]),
             "elapsed_sec": elapsed[i], "n_failures": nfail[i], "host": host_names[i]}
            for i in range(n.value)
        ]

    def cancel_runs(self, run_ids: list[int]) -> int:
        """Give up on these runs. Returns how many were actually cancelled."""
        arr = (c_int * len(run_ids))(*run_ids)
        n = c_int()
        self._check(self.lib.pestpp_cancel_runs(self.handle, arr, len(run_ids), byref(n)),
                    "pestpp_cancel_runs")
        return n.value

    def release_workers(self, worker_idxs: list[int] | None = None) -> int:
        """Hand workers back. None (or an empty list) releases every one.

        A worker holding a run is released too and its run is rescheduled - put back at the
        front of the queue for someone else - not failed and not cancelled. Returns how many
        were actually released, which is the authority: asking for eight and getting three is
        a normal outcome, not an error.
        """
        if worker_idxs:
            arr = (c_int * len(worker_idxs))(*worker_idxs)
            count = len(worker_idxs)
        else:
            arr, count = None, 0
        n = c_int()
        self._check(self.lib.pestpp_release_workers(self.handle, arr, count, byref(n)),
                    "pestpp_release_workers")
        return n.value

    def jacobian_prepare(self, calc_init_obs: bool = False) -> int:
        """Queue the Jacobian runs without running them. Returns how many."""
        n = c_int()
        self._check(self.lib.pestpp_jacobian_prepare(
            self.handle, 1 if calc_init_obs else 0, byref(n)), "pestpp_jacobian_prepare")
        return n.value

    def jacobian_run(self) -> None:
        """Run the queued Jacobian batch."""
        self._check(self.lib.pestpp_jacobian_run(self.handle), "pestpp_jacobian_run")

    def jacobian_process(self) -> bool:
        """Harvest the completed runs into the Jacobian."""
        ok = c_int()
        self._check(self.lib.pestpp_jacobian_process(self.handle, byref(ok)),
                    "pestpp_jacobian_process")
        return bool(ok.value)

    def get_jacobian(self):
        """The Jacobian as (values, row_names, col_names). Rows are observations.

        A dense COPY in numpy F order, matching how the library writes it. glm only; the
        ensemble tools refuse, because they never form the matrix.
        """
        nr, nc = c_int(), c_int()
        self._check(self.lib.pestpp_get_jacobian(self.handle, None, 0, 0, byref(nr), byref(nc),
                                                 None, None), "pestpp_get_jacobian")
        n_r, n_c = nr.value, nc.value
        data = np.empty((n_r, n_c), dtype=np.float64, order="F")
        rbuf = create_string_buffer(n_r * self.name_len)
        cbuf = create_string_buffer(n_c * self.name_len)
        self._check(self.lib.pestpp_get_jacobian(
            self.handle, data.ctypes.data_as(POINTER(c_double)), n_r, n_c,
            byref(nr), byref(nc), rbuf, cbuf), "pestpp_get_jacobian")
        return data, self._unpack_names(rbuf.raw, n_r), self._unpack_names(cbuf.raw, n_c)

    def get_par_vector(self, which: int = 0):
        """The tool's single parameter vector as (names, values). which: 0 current, 1 optimum."""
        n = c_int()
        self._check(self.lib.pestpp_get_par_vector(self.handle, which, None, 0, byref(n), None),
                    "pestpp_get_par_vector")
        cnt = n.value
        vals = np.empty(cnt, dtype=np.float64)
        nbuf = create_string_buffer(cnt * self.name_len)
        self._check(self.lib.pestpp_get_par_vector(
            self.handle, which, vals.ctypes.data_as(POINTER(c_double)), cnt, byref(n), nbuf),
            "pestpp_get_par_vector")
        return self._unpack_names(nbuf.raw, cnt), vals

    def set_par_vector(self, names, values) -> None:
        """Push values onto the current parameter vector, matched by name. Partial is fine."""
        names = list(names)
        vals = np.asarray(values, dtype=np.float64)
        if len(names) != vals.size:
            raise ValueError("names and values differ in length: {0} vs {1}"
                             .format(len(names), vals.size))
        buf = self._pack_names(names)
        self._check(self.lib.pestpp_set_par_vector(
            self.handle, vals.ctypes.data_as(POINTER(c_double)), len(names), buf),
            "pestpp_set_par_vector")

    def get_worker_count(self) -> int:
        v = c_int()
        self._check(self.lib.pestpp_get_worker_count(self.handle, byref(v)),
                    "pestpp_get_worker_count")
        return v.value

    def get_worker_state(self, idx: int) -> dict:
        host = create_string_buffer(self.name_len)
        state = create_string_buffer(self.name_len)
        run_id, npings = c_int(), c_int()
        elapsed, avg = c_double(), c_double()
        self._check(self.lib.pestpp_get_worker_state(
            self.handle, c_int(idx), host, self.name_len, state, self.name_len,
            byref(run_id), byref(elapsed), byref(avg), byref(npings)),
            "pestpp_get_worker_state")
        return {"host": host.value.decode(), "state": state.value.decode(),
                "current_run_id": run_id.value, "current_elapsed_sec": elapsed.value,
                "avg_runtime_sec": avg.value, "n_failed_pings": npings.value}

    def get_worker_run_history(self, idx: int, which: int) -> list[int]:
        """which is WORKER_COMPLETED / WORKER_FAILED / WORKER_TIMED_OUT."""
        n = c_int()
        self._check(self.lib.pestpp_get_worker_run_history(
            self.handle, c_int(idx), c_int(which), None, 0, byref(n)),
            "pestpp_get_worker_run_history")
        if n.value == 0:
            return []
        arr = (c_int * n.value)()
        self._check(self.lib.pestpp_get_worker_run_history(
            self.handle, c_int(idx), c_int(which), arr, n.value, byref(n)),
            "pestpp_get_worker_run_history")
        return list(arr)

    # -- queue / process ---------------------------------------------------------------

    def queue_runs(self) -> int:
        """Queue the current parameter ensemble. Returns the number of runs queued."""
        n = c_int()
        self._check(self.lib.pestpp_queue_runs(self.handle, byref(n)), "pestpp_queue_runs")
        return n.value

    def process_runs(self) -> int:
        """Collect the queued runs into the obs ensemble. Returns the number that failed."""
        n = c_int()
        self._check(self.lib.pestpp_process_runs(self.handle, byref(n)),
                    "pestpp_process_runs")
        return n.value

    # -- membership --------------------------------------------------------------------

    def _pack_names(self, names: list[str]):
        buf = create_string_buffer(len(names) * self.name_len)
        for i, n in enumerate(names):
            enc = n.encode().ljust(self.name_len)[:self.name_len]
            buf[i * self.name_len:(i + 1) * self.name_len] = enc
        return buf

    def drop_realizations(self, names: list[str]) -> None:
        """Drop these realizations from every coupled ensemble at once."""
        self._check(self.lib.pestpp_drop_realizations(
            self.handle, self._pack_names(names), len(names)), "pestpp_drop_realizations")

    def keep_realizations(self, names: list[str]) -> None:
        """Keep only these, in this order -- the inverse of drop, and also a reorder."""
        self._check(self.lib.pestpp_keep_realizations(
            self.handle, self._pack_names(names), len(names)), "pestpp_keep_realizations")

    # -- parameter snapshot ------------------------------------------------------------

    def get_phi_summary(self, phi_type: int = 0) -> dict:
        """mean/std/min/max of phi across realizations. ies and da only."""
        vals = [c_double() for _ in range(4)]
        self._check(self.lib.pestpp_get_phi_summary(
            self.handle, c_int(phi_type), *[byref(v) for v in vals]), "pestpp_get_phi_summary")
        return dict(zip(("mean", "std", "min", "max"), [v.value for v in vals]))

    def get_phi_vector(self, phi_type: int = 0):
        """Phi per realization: (values, names). Empty when that phi type is not in play."""
        n = c_int()
        self._check(self.lib.pestpp_get_phi_vector(
            self.handle, c_int(phi_type), None, None, 0, byref(n)), "pestpp_get_phi_vector")
        if n.value == 0:
            return np.empty(0), []
        vals = (c_double * n.value)()
        names = create_string_buffer(n.value * self.name_len)
        self._check(self.lib.pestpp_get_phi_vector(
            self.handle, c_int(phi_type), vals, names, n.value, byref(n)),
            "pestpp_get_phi_vector")
        return np.array(vals), self._unpack_names(names.raw, n.value)

    def get_phi_residuals(self, phi_type: int = 0):
        """Squared weighted residuals: (array, row_names, col_names). MEAS and ACTUAL only."""
        nrow, ncol = c_int(), c_int()
        self._check(self.lib.pestpp_get_phi_residuals(
            self.handle, c_int(phi_type), None, 0, 0, byref(nrow), byref(ncol), None, None),
            "pestpp_get_phi_residuals")
        nr, nc = nrow.value, ncol.value
        if nr == 0 or nc == 0:
            return np.empty((0, 0)), [], []
        data = (c_double * (nr * nc))()
        rows = create_string_buffer(nr * self.name_len)
        cols = create_string_buffer(nc * self.name_len)
        self._check(self.lib.pestpp_get_phi_residuals(
            self.handle, c_int(phi_type), data, nr, nc, byref(nrow), byref(ncol), rows, cols),
            "pestpp_get_phi_residuals")
        arr = np.ctypeslib.as_array(data, shape=(nr * nc,)).reshape((nr, nc), order="F").copy()
        return arr, self._unpack_names(rows.raw, nr), self._unpack_names(cols.raw, nc)

    def get_obs_groups(self) -> list:
        """The group each observation belongs to, aligned with the obs ensemble columns."""
        n = c_int()
        self._check(self.lib.pestpp_get_obs_groups(self.handle, None, 0, byref(n)),
                    "pestpp_get_obs_groups")
        if n.value == 0:
            return []
        buf = create_string_buffer(n.value * self.name_len)
        self._check(self.lib.pestpp_get_obs_groups(
            self.handle, buf, n.value * self.name_len, byref(n)), "pestpp_get_obs_groups")
        return self._unpack_names(buf.raw, n.value)

    def get_obs_weights(self) -> np.ndarray:
        """The control-file weight for each observation, in obs-ensemble column order."""
        n = c_int()
        self._check(self.lib.pestpp_get_obs_weights(self.handle, None, 0, byref(n)),
                    "pestpp_get_obs_weights")
        if n.value == 0:
            return np.empty(0)
        w = (c_double * n.value)()
        self._check(self.lib.pestpp_get_obs_weights(self.handle, w, n.value, byref(n)),
                    "pestpp_get_obs_weights")
        return np.array(w)

    def set_obs_weights(self, names, weights) -> None:
        """Set control-file weights by name. See broadcast_weights()."""
        names = list(names)
        arr = (c_double * len(names))(*[float(x) for x in weights])
        buf = create_string_buffer(len(names) * self.name_len)
        for i, nm in enumerate(names):
            buf[i * self.name_len:(i + 1) * self.name_len] = \
                nm.encode().ljust(self.name_len)[:self.name_len]
        self._check(self.lib.pestpp_set_obs_weights(self.handle, buf, arr, len(names)),
                    "pestpp_set_obs_weights")

    def reinflate_ensemble(self, factor: float = 1.0, num_reals: int = 0,
                           center_on_min_phi: int = -1) -> None:
        """Rebuild the parameter ensemble from the prior's spread, re-centred on the current."""
        self._check(self.lib.pestpp_reinflate_ensemble(
            self.handle, c_double(float(factor)), c_int(int(num_reals)),
            c_int(int(center_on_min_phi))),
            "pestpp_reinflate_ensemble")

    def broadcast_weights(self) -> None:
        """Push the weight vector into every row of the weights ensemble."""
        self._check(self.lib.pestpp_broadcast_weights(self.handle), "pestpp_broadcast_weights")

    def update_phi(self) -> None:
        """Recompute the cached phi from the current ensembles and weights."""
        self._check(self.lib.pestpp_update_phi(self.handle), "pestpp_update_phi")

    def get_par_snapshot(self):
        """CTL-space copy of every control-file parameter: (values, row_names, col_names).

        A copy, unlike get_ensemble_view -- safe to hold, and the form set_par_snapshot takes.
        """
        nrow, ncol = c_int(), c_int()
        self._check(self.lib.pestpp_get_par_snapshot(
            self.handle, None, 0, 0, byref(nrow), byref(ncol), None, None),
            "pestpp_get_par_snapshot")
        nr, nc = nrow.value, ncol.value
        data = (c_double * (nr * nc))()
        rows = create_string_buffer(nr * self.name_len)
        cols = create_string_buffer(nc * self.name_len)
        self._check(self.lib.pestpp_get_par_snapshot(
            self.handle, data, nr, nc, byref(nrow), byref(ncol), rows, cols),
            "pestpp_get_par_snapshot")
        arr = np.ctypeslib.as_array(data, shape=(nr * nc,)).reshape((nr, nc), order="F").copy()
        return arr, self._unpack_names(rows.raw, nr), self._unpack_names(cols.raw, nc)

    def set_par_snapshot(self, values, row_names: list[str], col_names: list[str]) -> None:
        """Push CTL-space values back. Matched by name, so order need not match."""
        arr = np.asfortranarray(np.asarray(values, dtype=np.float64))
        nr, nc = arr.shape
        if nr != len(row_names) or nc != len(col_names):
            raise PestppError(
                f"snapshot shape {arr.shape} does not match "
                f"{len(row_names)} row names x {len(col_names)} col names")
        flat = (c_double * (nr * nc))(*arr.flatten(order="F"))
        self._check(self.lib.pestpp_set_par_snapshot(
            self.handle, flat, nr, nc, self._pack_names(row_names), self._pack_names(col_names)),
            "pestpp_set_par_snapshot")

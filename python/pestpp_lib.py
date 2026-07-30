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
from ctypes import (
    CDLL,
    POINTER,
    Structure,
    byref,
    c_char,
    c_char_p,
    c_double,
    c_int,
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
PAR_EN, OBS_EN, NOISE_EN, WEIGHTS_EN = 0, 1, 2, 3
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
        lib.pestpp_get_worker_count.argtypes = (c_void_p, POINTER(c_int))
        lib.pestpp_get_worker_count.restype = c_int
        lib.pestpp_get_worker_state.argtypes = (
            c_void_p, c_int, c_char_p, c_int, c_char_p, c_int, POINTER(c_int),
            POINTER(c_double), POINTER(c_double), POINTER(c_int))
        lib.pestpp_get_worker_state.restype = c_int
        lib.pestpp_get_worker_run_history.argtypes = (
            c_void_p, c_int, c_int, POINTER(c_int), c_int, POINTER(c_int))
        lib.pestpp_get_worker_run_history.restype = c_int

        # -- queue/harvest, membership, snapshot --
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
        if getattr(self, "handle", None) is not None and self.handle:
            self.lib.pestpp_destroy(self.handle)
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

    # -- queue / harvest ---------------------------------------------------------------

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

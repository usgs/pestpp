"""Thin ctypes layer over the PEST++ C ABI.

This mirrors ``pestpp_capi.h`` one-to-one on purpose. Method names match the C symbols so
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
    byref,
    c_char,
    c_char_p,
    c_double,
    c_int,
    create_string_buffer,
)

import numpy as np

# ---- enums, mirroring the header ----------------------------------------------------

PESTPP_OK = 0
PESTPP_ERROR = 1
PESTPP_INVALID_HANDLE = 2
PESTPP_RETRY = 3

TOOL_IES, TOOL_DA, TOOL_MOU, TOOL_SQP = 0, 1, 2, 3
PAR_EN, OBS_EN, NOISE_EN, WEIGHTS_EN = 0, 1, 2, 3
TSTAT = {0: "ctl", 1: "num", 2: "model"}

RUN_QUEUED, RUN_RUNNING, RUN_COMPLETED, RUN_FAILED, RUN_TIMED_OUT, RUN_CANCELLED = range(6)
RUN_STATUS = {
    RUN_QUEUED: "queued", RUN_RUNNING: "running", RUN_COMPLETED: "completed",
    RUN_FAILED: "failed", RUN_TIMED_OUT: "timed_out", RUN_CANCELLED: "cancelled",
}
WORKER_COMPLETED, WORKER_FAILED, WORKER_TIMED_OUT = 0, 1, 2


class PestppError(Exception):
    """Raised when a C ABI call returns a non-zero status.

    Its own type rather than RuntimeError, so callers can catch precisely.
    """


class PestppLib:
    """One loaded shared library plus one session handle."""

    def __init__(self, lib_path: str, tool: int, ctl_file: str, working_dir: str = ".",
                 port: str | None = None):
        """port=None uses the serial run manager; a port starts a PANTHER master.

        Only the PANTHER master can be observed or interrupted mid-batch -- see
        supports_live_control().
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
        if port is None:
            status = self.lib.pestpp_create(
                c_int(tool), ctl_file.encode(), working_dir.encode(), byref(self.handle))
        else:
            status = self.lib.pestpp_create_panther(
                c_int(tool), ctl_file.encode(), working_dir.encode(), str(port).encode(),
                byref(self.handle))
        if status != PESTPP_OK:
            raise PestppError(self.lib.pestpp_last_create_error().decode())

    # -- plumbing ---------------------------------------------------------------------

    def _prototype(self) -> None:
        """Declare argtypes/restype for every function.

        Without this ctypes silently truncates pointers on 64-bit, which fails in ways that
        look like data corruption rather than a type error.
        """
        from ctypes import c_void_p

        lib = self.lib
        lib.pestpp_create.argtypes = (c_int, c_char_p, c_char_p, POINTER(c_void_p))
        lib.pestpp_create.restype = c_int
        lib.pestpp_create_panther.argtypes = (c_int, c_char_p, c_char_p, c_char_p,
                                              POINTER(c_void_p))
        lib.pestpp_create_panther.restype = c_int
        lib.pestpp_destroy.argtypes = (c_void_p,)
        lib.pestpp_destroy.restype = c_int
        lib.pestpp_last_error.argtypes = (c_void_p,)
        lib.pestpp_last_error.restype = c_char_p
        lib.pestpp_last_create_error.argtypes = ()
        lib.pestpp_last_create_error.restype = c_char_p

        for name in ("pestpp_initialize", "pestpp_solve_iteration", "pestpp_finalize"):
            getattr(lib, name).argtypes = (c_void_p,)
            getattr(lib, name).restype = c_int

        for name in ("pestpp_get_iteration", "pestpp_should_terminate",
                     "pestpp_get_par_transform_status"):
            getattr(lib, name).argtypes = (c_void_p, POINTER(c_int))
            getattr(lib, name).restype = c_int

        lib.pestpp_get_ensemble_view.argtypes = (
            c_void_p, c_int, POINTER(POINTER(c_double)), POINTER(c_int), POINTER(c_int))
        lib.pestpp_get_ensemble_view.restype = c_int

        for name in ("pestpp_get_ensemble_row_names", "pestpp_get_ensemble_col_names"):
            getattr(lib, name).argtypes = (c_void_p, c_int, c_char_p, c_int, POINTER(c_int))
            getattr(lib, name).restype = c_int

        lib.pestpp_set_option.argtypes = (c_void_p, c_char_p, c_char_p)
        lib.pestpp_set_option.restype = c_int
        lib.pestpp_get_option.argtypes = (c_void_p, c_char_p, c_char_p, c_int, POINTER(c_int))
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
        lib.pestpp_queue_ensemble.argtypes = (c_void_p, POINTER(c_int))
        lib.pestpp_queue_ensemble.restype = c_int
        lib.pestpp_harvest_ensemble.argtypes = (c_void_p, POINTER(c_int))
        lib.pestpp_harvest_ensemble.restype = c_int
        for name in ("pestpp_drop_realizations", "pestpp_keep_realizations"):
            getattr(lib, name).argtypes = (c_void_p, c_char_p, c_int)
            getattr(lib, name).restype = c_int
        lib.pestpp_get_par_snapshot.argtypes = (
            c_void_p, POINTER(c_double), c_int, c_int, POINTER(c_int), POINTER(c_int),
            c_char_p, c_char_p)
        lib.pestpp_get_par_snapshot.restype = c_int
        lib.pestpp_set_par_snapshot.argtypes = (
            c_void_p, POINTER(c_double), c_int, c_int, c_char_p, c_char_p)
        lib.pestpp_set_par_snapshot.restype = c_int

    def _check(self, status: int, what: str) -> int:
        """Raise on error; pass RETRY through, since it is an outcome and not a fault."""
        if status in (PESTPP_OK, PESTPP_RETRY):
            return status
        msg = self.lib.pestpp_last_error(self.handle)
        raise PestppError(f"{what}: {msg.decode() if msg else 'unknown error'}")

    def _unpack_names(self, raw: bytes, count: int) -> list[str]:
        return [
            raw[i * self.name_len:(i + 1) * self.name_len].decode().strip()
            for i in range(count)
        ]

    # -- lifecycle --------------------------------------------------------------------

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

    def get_ensemble_view(self, ensemble_id: int) -> np.ndarray:
        """A numpy view onto the live buffer -- no copy.

        Eigen is column-major, hence order="F". The array is only valid while that
        ensemble's storage is unchanged; re-fetch after anything that could mutate it.
        """
        ptr = POINTER(c_double)()
        nrow, ncol = c_int(), c_int()
        self._check(
            self.lib.pestpp_get_ensemble_view(
                self.handle, c_int(ensemble_id), byref(ptr), byref(nrow), byref(ncol)),
            "pestpp_get_ensemble_view")
        n = nrow.value * ncol.value
        if n == 0:
            return np.empty((nrow.value, ncol.value), order="F")
        buf = np.ctypeslib.as_array(ptr, shape=(n,))
        return buf.reshape((nrow.value, ncol.value), order="F")

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
        self._check(self.lib.pestpp_set_option(self.handle, key.encode(), str(value).encode()),
                    f"pestpp_set_option({key})")

    def get_option(self, key: str) -> str:
        needed = c_int()
        self._check(self.lib.pestpp_get_option(self.handle, key.encode(), None, 0, byref(needed)),
                    f"pestpp_get_option({key})")
        buf = create_string_buffer(needed.value)
        self._check(self.lib.pestpp_get_option(self.handle, key.encode(), buf, needed.value,
                                               byref(needed)), f"pestpp_get_option({key})")
        return buf.value.decode()

    # -- run management ----------------------------------------------------------------

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

    def queue_ensemble(self) -> int:
        """Queue the current parameter ensemble. Returns the number of runs queued."""
        n = c_int()
        self._check(self.lib.pestpp_queue_ensemble(self.handle, byref(n)), "pestpp_queue_ensemble")
        return n.value

    def harvest_ensemble(self) -> int:
        """Collect the queued runs into the obs ensemble. Returns the number that failed."""
        n = c_int()
        self._check(self.lib.pestpp_harvest_ensemble(self.handle, byref(n)),
                    "pestpp_harvest_ensemble")
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

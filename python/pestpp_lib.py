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


class PestppError(Exception):
    """Raised when a C ABI call returns a non-zero status.

    Its own type rather than RuntimeError, so callers can catch precisely.
    """


class PestppLib:
    """One loaded shared library plus one session handle."""

    def __init__(self, lib_path: str, tool: int, ctl_file: str, working_dir: str = "."):
        if not os.path.exists(lib_path):
            raise FileNotFoundError(lib_path)
        self.lib = CDLL(lib_path)
        self._prototype()
        # Read the name width from the library rather than hardcoding it -- the same trick
        # xmipy and pypestutils use, and it removes a whole class of drift.
        self.name_len = c_int.in_dll(self.lib, "PESTPP_NAME_LEN").value

        self.handle = c_int()  # placeholder; replaced by the void* below
        from ctypes import c_void_p

        self.handle = c_void_p()
        status = self.lib.pestpp_create(
            c_int(tool),
            ctl_file.encode(),
            working_dir.encode(),
            byref(self.handle),
        )
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

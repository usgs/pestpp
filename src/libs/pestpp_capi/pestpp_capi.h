/*
 * PEST++ C ABI - walking skeleton.
 *
 * A flat, C-callable surface over the PEST++ tools, intended to be driven from Python via
 * ctypes (see python/pestpp_lib.py) or from any language that can call C. Design rationale
 * and the full intended surface are in docs/api_part1/c_abi_design.md; this header is the
 * first vertical slice of that design - enough to create a handle, drive iterations, and
 * read an ensemble without copying.
 *
 * Conventions, all deliberate:
 *   - Every function returns pestpp_status. Nothing throws across the boundary; exceptions
 *     are caught at every entry point and turned into a status plus a per-handle message.
 *   - Strings out are written into caller-allocated buffers, with a `needed` out-param so a
 *     caller can size then fill. The library never allocates something the caller must free.
 *   - The working directory is an argument, never process-global state.
 */
#ifndef PESTPP_CAPI_H_
#define PESTPP_CAPI_H_

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32)
#  define PESTPP_API __declspec(dllexport)
#else
#  define PESTPP_API __attribute__((visibility("default")))
#endif

/* Opaque session handle. Owns the scenario, file manager, run manager and tool object. */
typedef void* pestpp_handle;

typedef enum {
    PESTPP_OK = 0,
    PESTPP_ERROR = 1,          /* something failed; see pestpp_last_error()          */
    PESTPP_INVALID_HANDLE = 2,
    PESTPP_RETRY = 3           /* algorithmic "no upgrade taken, try again" - not a fault */
} pestpp_status;

typedef enum {
    PESTPP_IES = 0,
    PESTPP_DA  = 1,
    PESTPP_MOU = 2,
    PESTPP_SQP = 3
} pestpp_tool;

/* Which ensemble to look at. The raw parameter ensemble is in whatever transform space the
   tool is currently using (usually NUM during a solve) and holds only adjustable
   parameters - see pestpp_get_par_transform_status(). */
typedef enum {
    PESTPP_PAR_EN     = 0,
    PESTPP_OBS_EN     = 1,
    PESTPP_NOISE_EN   = 2,
    PESTPP_WEIGHTS_EN = 3
} pestpp_ensemble_id;

/* Length of a single name in the packed name buffers used below. Read this from the library
   rather than hard-coding it, the way xmipy and pypestutils do. */
PESTPP_API extern const int PESTPP_NAME_LEN;

/* ---- lifecycle -------------------------------------------------------------------- */

/* Create a session. `ctl_file` is resolved relative to `working_dir`. */
PESTPP_API pestpp_status pestpp_create(int tool, const char* ctl_file,
                                       const char* working_dir, pestpp_handle* out);
PESTPP_API pestpp_status pestpp_destroy(pestpp_handle h);

/* Message from the most recent failed call on this handle; "" if none. Valid until the
   next call on the same handle. Never returns NULL. */
PESTPP_API const char* pestpp_last_error(pestpp_handle h);
/* For pestpp_create() failures, where there is no handle to ask. */
PESTPP_API const char* pestpp_last_create_error(void);

/* ---- driving the algorithm --------------------------------------------------------- */

PESTPP_API pestpp_status pestpp_initialize(pestpp_handle h);

/* One iteration. Returns PESTPP_RETRY where the algorithm rejected the upgrade and wants
   another attempt with a new lambda - a distinct outcome from failure. */
PESTPP_API pestpp_status pestpp_solve_iteration(pestpp_handle h);

PESTPP_API pestpp_status pestpp_finalize(pestpp_handle h);
PESTPP_API pestpp_status pestpp_get_iteration(pestpp_handle h, int* iter);
PESTPP_API pestpp_status pestpp_should_terminate(pestpp_handle h, int* out);

/* ---- reading ensembles ------------------------------------------------------------- */

/* Zero-copy. `data` points at the tool's live Eigen buffer, which is COLUMN-major:
   element (i,j) is data[i + j*(*nrow)]. It stays valid only until something changes that
   ensemble's storage, so treat it as borrowed for the current step and re-fetch after any
   call that could mutate. */
PESTPP_API pestpp_status pestpp_get_ensemble_view(pestpp_handle h, int ensemble_id,
                                                  double** data, int* nrow, int* ncol);

/* Names packed as fixed-width PESTPP_NAME_LEN blocks, space padded, not NUL terminated -
   the same shape MODFLOW 6 uses for its variable-name blocks. Pass buf=NULL to query
   `count` only. */
PESTPP_API pestpp_status pestpp_get_ensemble_row_names(pestpp_handle h, int ensemble_id,
                                                       char* buf, int buf_len, int* count);
PESTPP_API pestpp_status pestpp_get_ensemble_col_names(pestpp_handle h, int ensemble_id,
                                                       char* buf, int buf_len, int* count);

/* 0=CTL, 1=NUM, 2=MODEL - tells a caller which space a PESTPP_PAR_EN view is in. */
PESTPP_API pestpp_status pestpp_get_par_transform_status(pestpp_handle h, int* tstat);

/* ---- options ----------------------------------------------------------------------- */

PESTPP_API pestpp_status pestpp_set_option(pestpp_handle h, const char* key, const char* value);
PESTPP_API pestpp_status pestpp_get_option(pestpp_handle h, const char* key,
                                           char* buf, int buf_len, int* needed);

#ifdef __cplusplus
}
#endif

#endif /* PESTPP_CAPI_H_ */

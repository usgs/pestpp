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

/* All four are driveable. What one "iteration" means differs by tool - ies solves lambdas,
   mou runs a generation, sqp may issue several run batches while it line-searches - but the
   calls a caller makes do not: initialize, then solve_iteration, then finalize.

   PESTPP_DA runs a SINGLE cycle. pestpp-da.cpp drives a sequence of cycles, rebuilding the
   scenario and the tool for each; that machinery is not exposed yet, so a multi-cycle control
   file will run only its first cycle here. */
typedef enum {
    PESTPP_IES = 0,
    PESTPP_DA  = 1,
    PESTPP_MOU = 2,
    PESTPP_SQP = 3
} pestpp_tool;

/* Which ensemble to look at. The raw parameter ensemble is in whatever transform space the
   tool is currently using (usually NUM during a solve) and holds only adjustable
   parameters - see pestpp_get_par_transform_status().

   PAR/OBS mean the tool's parameter and result ensembles: pe/oe for ies and da, the decision
   variable population dp and its results op for mou, dv/oe for sqp. NOISE and WEIGHTS exist
   only for ies and da; asking mou or sqp for them is an error rather than an empty array. */
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

/* Which run manager backs the session.
 *
 * Only PANTHER can be watched or interrupted mid-batch. Serial and external finish the whole
 * batch in the first run_slice() and answer the run/worker queries with an error - ask
 * pestpp_supports_live_control() rather than discovering it. */
typedef enum {
    PESTPP_RM_SERIAL   = 0,
    PESTPP_RM_PANTHER  = 1,   /* needs panther_port */
    PESTPP_RM_EXTERNAL = 2
} pestpp_run_manager;

/* Create-time settings.
 *
 * A struct rather than a parameter list because this will grow - da multi-cycle wants a start
 * cycle, a restart/hotstart file is an obvious next ask - and every addition to a flat
 * signature is either a new symbol or a break. `struct_size` is how the library knows which
 * version of the struct it was handed: set it to sizeof(pestpp_create_options) and zero the
 * rest, and fields added later default sensibly.
 *
 * Everything else a run manager needs is already a regular option, so set it with
 * pestpp_set_option() rather than looking for it here. */
typedef struct {
    int         struct_size;     /* = sizeof(pestpp_create_options); required             */
    int         tool;            /* pestpp_tool                                           */
    const char* ctl_file;        /* resolved relative to working_dir                      */
    const char* working_dir;     /* NULL or "" means the current directory                */
    int         run_manager;     /* pestpp_run_manager; 0 (serial) if left zeroed         */
    const char* panther_port;    /* PANTHER only; NULL otherwise                          */
} pestpp_create_options;

PESTPP_API pestpp_status pestpp_create(const pestpp_create_options* opts, pestpp_handle* out);
PESTPP_API pestpp_status pestpp_destroy(pestpp_handle h);

/* Which run manager this handle actually got. */
PESTPP_API pestpp_status pestpp_get_run_manager(pestpp_handle h, int* run_manager);

/* Message from the most recent failed call on this handle; "" if none. Valid until the
   next call on the same handle. Never returns NULL. */
PESTPP_API const char* pestpp_last_error(pestpp_handle h);
/* For pestpp_create() failures, where there is no handle to ask. */
PESTPP_API const char* pestpp_last_create_error(void);

/* Flush the library's console output.
 *
 * Needed because on windows the library links the static CRT (/MT), so it owns a COPY of the
 * C runtime with its own stdout buffer, independent of the host program's. A caller that
 * redirects file descriptor 1 to capture output gets the model's output - child processes
 * inherit the descriptor - but the library's own text stays in that private buffer and
 * surfaces later, after the redirect has been undone. Call this before restoring, and the
 * buffer lands where it was meant to. Harmless everywhere else. */
PESTPP_API pestpp_status pestpp_flush_output(void);

/* Send the library's console output to a file, and put it back.
 *
 * This has to happen INSIDE the library. On windows it links the static CRT (/MT), so it owns
 * a private copy of the C runtime with its own file descriptor table: a host program that
 * redirects its own descriptor 1 moves the process std handle - so child processes such as
 * the model land in the file - while the library's descriptor 1 still refers to the original
 * console. Flushing does not help; the descriptor itself is the wrong one. Doing the dup2 in
 * here operates on the table that actually matters, and catches the child processes too.
 *
 * redirect() writes the descriptor to restore into `saved_fd`; hand that back to restore().
 * Appends, so repeated redirects to one path accumulate rather than truncating. */
PESTPP_API pestpp_status pestpp_redirect_output(const char* path, int* saved_fd);
PESTPP_API pestpp_status pestpp_restore_output(int saved_fd);

/* ---- driving the algorithm --------------------------------------------------------- */

PESTPP_API pestpp_status pestpp_initialize(pestpp_handle h);

/* initialize(), split so the caller owns the prior-ensemble evaluation.
 *
 *     int n;
 *     pestpp_initialize_prepare(h, &n);   // draws the ensembles; runs nothing
 *     if (n) { queue -> begin_batch/run_slice/end_batch -> harvest }
 *     pestpp_initialize_finish(h);        // phi, drop failures, write the .0. files
 *
 * `n_runs` is how many runs the caller must service. 0 is normal and needs no branch beyond
 * skipping the batch: a restart or `ies_obs_restart_csv` supplies results instead of
 * computing them, and mou and sqp initialize atomically and hand over nothing.
 *
 * The window between prepare and queue is the useful part: it is the only point at which a
 * caller can REPLACE the drawn prior ensemble with its own - via pestpp_set_par_snapshot()
 * or a write through the zero-copy view - and have those realizations be what gets run.
 *
 * Calling pestpp_solve_iteration() while a prepare is outstanding is an error: the tool
 * would be half-initialized, its ensembles drawn but failures and phi unprocessed. */
PESTPP_API pestpp_status pestpp_initialize_prepare(pestpp_handle h, int* n_runs);
PESTPP_API pestpp_status pestpp_initialize_finish(pestpp_handle h);

/* One iteration. Returns PESTPP_RETRY where the algorithm rejected the upgrade and wants
   another attempt with a new lambda - a distinct outcome from failure. */
PESTPP_API pestpp_status pestpp_solve_iteration(pestpp_handle h);

PESTPP_API pestpp_status pestpp_finalize(pestpp_handle h);
PESTPP_API pestpp_status pestpp_get_iteration(pestpp_handle h, int* iter);

/* Which phi. MEAS is the measurement objective function most people mean by "phi"; ACTUAL
   drops the noise realizations. ies and da only - mou has objectives rather than a phi, and
   sqp has an objective function, so both report an error. */
typedef enum {
    PESTPP_PHI_MEAS      = 0,
    PESTPP_PHI_COMPOSITE = 1,
    PESTPP_PHI_REGUL     = 2,
    PESTPP_PHI_ACTUAL    = 3,
    PESTPP_PHI_NOISE     = 4
} pestpp_phi_type;

/* Summary of phi across the realizations. Any out-param may be NULL. Only meaningful once
   the ensemble has been evaluated. */
PESTPP_API pestpp_status pestpp_get_phi_summary(pestpp_handle h, int phi_type,
                                                double* mean, double* std,
                                                double* min, double* max);

/* Phi for every realization, not just the summary. `names` are packed PESTPP_NAME_LEN wide,
   the same shape as the ensemble name buffers, and are the realizations phi is defined for -
   which need not be every row of the parameter ensemble, so pair the two rather than assuming
   the order. Pass phi=NULL, names=NULL, max_n=0 to learn the count first.

   A phi type that is not in play (REGUL outside a regularized run, NOISE with ies_no_noise)
   reports a count of 0 rather than failing: it is absent, not an error. */
PESTPP_API pestpp_status pestpp_get_phi_vector(pestpp_handle h, int phi_type,
                                               double* phi, char* names,
                                               int max_n, int* n_out);

/* The individual terms phi is the sum of: the squared weighted residual for every
   (realization, observation) pair. Summing a row gives that realization's phi; summing the
   columns belonging to an observation group gives that group's contribution, which is what
   the tools' "observation group phi summary" reports.

   Column-major like the ensemble views: element (i,j) is data[i + j*(*nrow)]. Rows are
   realizations, columns observations, and both name lists are returned so a caller never has
   to assume an order. Pass data=NULL to learn the shape first.

   PESTPP_PHI_MEAS and PESTPP_PHI_ACTUAL only - the residual is not defined for the others. */
PESTPP_API pestpp_status pestpp_get_phi_residuals(pestpp_handle h, int phi_type,
                                                  double* data, int max_nrow, int max_ncol,
                                                  int* nrow, int* ncol,
                                                  char* row_names, char* col_names);

/* The observation group each observation belongs to, in the same order as
   pestpp_get_ensemble_col_names(PESTPP_OBS_EN). Packed PESTPP_NAME_LEN wide. */
PESTPP_API pestpp_status pestpp_get_obs_groups(pestpp_handle h, char* buf, int buf_len,
                                               int* count);

/* ---- weights --------------------------------------------------------------------------
 *
 * There are two of them and they are not the same thing:
 *
 *   the WEIGHT VECTOR   - one weight per observation, from the control file
 *   the WEIGHTS ENSEMBLE - one weight per observation PER REALIZATION, reachable as
 *                          PESTPP_WEIGHTS_EN through the zero-copy ensemble view
 *
 * initialize() reads the vector once and broadcasts it into every row of the ensemble, and
 * from then on phi is computed from the ENSEMBLE. So setting the vector after initialization
 * changes nothing on its own - call pestpp_broadcast_weights() to push it through, or write
 * the ensemble directly when you want weights to differ between realizations. */

/* Weights in the same order as pestpp_get_ensemble_col_names(PESTPP_OBS_EN). */
PESTPP_API pestpp_status pestpp_get_obs_weights(pestpp_handle h, double* weights,
                                                int max_n, int* n_out);

/* Set weights by name, so a caller can change a few without restating the rest. */
PESTPP_API pestpp_status pestpp_set_obs_weights(pestpp_handle h, const char* names,
                                                const double* weights, int n);

/* Push the current weight vector into every row of the weights ensemble, so a vector change
   takes effect. Overwrites per-realization weights - do not use it if you have been writing
   PESTPP_WEIGHTS_EN yourself. No-op before the ensemble exists. */
PESTPP_API pestpp_status pestpp_broadcast_weights(pestpp_handle h);

/* Recompute phi from the current ensembles and weights.
 *
 * Phi is CACHED - the values pestpp_get_phi_* report are whatever the last update computed,
 * which the algorithm does at its own points. Change weights or write an ensemble and the
 * cached phi is stale until this is called. */
PESTPP_API pestpp_status pestpp_update_phi(pestpp_handle h);
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

/* ---- run management ------------------------------------------------------------------
 *
 * A caller driving its own loop queues runs, drives the run manager in slices, and harvests
 * when it is ready. Slicing is what makes the run manager observable and interruptible: each
 * slice returns after roughly max_seconds so the caller can inspect or cancel between them.
 *
 * The introspection and cancel calls below are PANTHER-only. The serial run manager has no
 * workers, cannot yield mid-batch, and finishes everything in the first slice. Rather than
 * silently pretending otherwise, those calls fail with a message - ask pestpp_supports_live_
 * control() first.
 */

typedef enum {
    PESTPP_RUN_QUEUED    = 0,
    PESTPP_RUN_RUNNING   = 1,
    PESTPP_RUN_COMPLETED = 2,
    PESTPP_RUN_FAILED    = 3,
    PESTPP_RUN_TIMED_OUT = 4,
    PESTPP_RUN_CANCELLED = 5
} pestpp_run_status;

/* Which of a worker's three run histories to read. */
typedef enum {
    PESTPP_WORKER_COMPLETED = 0,
    PESTPP_WORKER_FAILED    = 1,
    PESTPP_WORKER_TIMED_OUT = 2
} pestpp_worker_history;

/* Nonzero when this handle's run manager can yield mid-batch and answer the queries below. */
PESTPP_API pestpp_status pestpp_supports_live_control(pestpp_handle h, int* out);

/* Drive the queued runs. begin_batch() once, then run_slice() until all_done, then
   end_batch(). max_seconds bounds one slice; it is ignored by the serial manager. */
PESTPP_API pestpp_status pestpp_begin_batch(pestpp_handle h);
PESTPP_API pestpp_status pestpp_run_slice(pestpp_handle h, double max_seconds, int* all_done);
PESTPP_API pestpp_status pestpp_end_batch(pestpp_handle h);

/* Aggregate counts for the batch in flight. Any out-param may be NULL. */
PESTPP_API pestpp_status pestpp_get_run_time_stats(pestpp_handle h, double* avg_run_sec,
                                                   int* n_completed, int* n_failed,
                                                   int* n_timed_out, int* n_queued,
                                                   int* n_running);

/* Per-run state. Pass want_ids=NULL to get every run, or an array of n_want ids to ask about
   specific ones. Pass max_n=0 with NULL arrays to learn the count first, then size and refill.
   Each out array may be NULL if that column is not wanted. `statuses` holds pestpp_run_status.
   `host` names are packed PESTPP_NAME_LEN-wide, like the ensemble name buffers. */
PESTPP_API pestpp_status pestpp_get_run_states(pestpp_handle h,
                                               const int* want_ids, int n_want,
                                               int* run_ids, int* statuses, double* elapsed_sec,
                                               int* n_failures, char* hosts,
                                               int max_n, int* n_out);

/* Give up on these runs. They come back as PESTPP_RUN_CANCELLED, and the count actually
   cancelled is written to n_cancelled (a run that already finished cannot be cancelled). */
PESTPP_API pestpp_status pestpp_cancel_runs(pestpp_handle h, const int* run_ids, int n,
                                            int* n_cancelled);

/* Workers. get_worker_count() first, then index 0..n-1. Any out-param may be NULL. */
PESTPP_API pestpp_status pestpp_get_worker_count(pestpp_handle h, int* n);
PESTPP_API pestpp_status pestpp_get_worker_state(pestpp_handle h, int idx,
                                                 char* host_buf, int host_buf_len,
                                                 char* state_buf, int state_buf_len,
                                                 int* current_run_id, double* current_elapsed_sec,
                                                 double* avg_runtime_sec, int* n_failed_pings);

/* Which runs this worker completed / failed / timed out. Same query-then-fill shape:
   pass run_ids=NULL, max_n=0 to get the count. */
PESTPP_API pestpp_status pestpp_get_worker_run_history(pestpp_handle h, int idx, int which,
                                                       int* run_ids, int max_n, int* n_out);

/* ---- queue / harvest -----------------------------------------------------------------
 *
 * The two halves of an evaluation, so a caller can own what happens in between: queue the
 * current parameter ensemble, drive and watch the run manager, then harvest into the
 * observation ensemble.
 *
 * Runs are tracked by realization NAME, so membership may change between the two calls.
 * A realization dropped after queueing has its run discarded rather than misattributed; one
 * added after queueing simply has no result yet. */
PESTPP_API pestpp_status pestpp_queue_runs(pestpp_handle h, int* n_queued);
PESTPP_API pestpp_status pestpp_process_runs(pestpp_handle h, int* n_failed);

/* ---- membership ----------------------------------------------------------------------
 *
 * Realizations do not travel alone: an ies handle holds pe, oe, the noise ensemble and the
 * weights, and they must agree on their realization set. These apply the change to all of
 * them at once, because making a caller keep four containers in step is a trap.
 *
 * Names are PAR-side realization names, packed PESTPP_NAME_LEN wide. The obs-side rows are
 * taken by position, which is how pest++ pairs par and obs realizations throughout. */
PESTPP_API pestpp_status pestpp_drop_realizations(pestpp_handle h, const char* names, int n);
PESTPP_API pestpp_status pestpp_keep_realizations(pestpp_handle h, const char* names, int n);

/* ---- parameter snapshot ---------------------------------------------------------------
 *
 * A CTL-space copy of every control-file parameter, including fixed and tied - the same
 * values to_csv() writes. This is the round-trippable form; the zero-copy view above is the
 * raw (usually NUM-space, adjustable-only) buffer. Data is column-major, like the view.
 *
 * Call with data=NULL to learn nrow/ncol, then size and refill. */
PESTPP_API pestpp_status pestpp_get_par_snapshot(pestpp_handle h, double* data,
                                                 int max_nrow, int max_ncol,
                                                 int* nrow, int* ncol,
                                                 char* row_names, char* col_names);

/* Push values back. Realizations and parameters are matched by NAME, so row and column order
   need not match what get_par_snapshot() returned - but every realization and every
   adjustable and fixed parameter currently held must be covered. Tied parameters are ignored
   because they are derived from their parents. This changes values, never membership. */
PESTPP_API pestpp_status pestpp_set_par_snapshot(pestpp_handle h, const double* data,
                                                 int nrow, int ncol,
                                                 const char* row_names, const char* col_names);

#ifdef __cplusplus
}
#endif

#endif /* PESTPP_CAPI_H_ */

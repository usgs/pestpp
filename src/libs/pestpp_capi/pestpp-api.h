/*
 * pest++ c api
 *
 * a flat set of c functions for driving the pest++ tools, meant to be called from python with
 * ctypes (see python/pestpp_lib.py and python/README.md) or from anything else that can call c.
 * why it is put together this way is in docs/api_part1/c_abi_design.md.
 *
 * stability: this is 0.x, so expect breaking changes. it is out for feedback and that feedback
 * is going to change it. until 1.0, check the PESTPP_API_VERSION_* values below - what you
 * built against - against what pestpp_get_api_version() reports - what you actually loaded -
 * and dont assume they match. the python files move with the library, so use them from the
 * same commit. if you need something that wont move on you, run the exes.
 *
 * how it works, on purpose:
 *   - every function returns pestpp_status. nothing throws across the boundary - exceptions
 *     get caught at every entry point and turned into a status plus a message you can ask for.
 *   - strings come back in buffers you allocate, with a `needed` out-param so you can size it
 *     first and then fill it. the library never hands you something you have to free.
 *   - everything is __cdecl on every platform. there is no __stdcall entry point, so on
 *     windows you want ctypes.CDLL, same as posix, not WinDLL.
 *
 * what this library does to your process
 *
 * it gets loaded into your process and it does things. before you hook it up to anything else,
 * know all five:
 *
 *   1. it changes the working directory. every function chdir's into the session's working_dir
 *      on the way in and puts it back on the way out. that is process wide while the call is
 *      running, so if another thread does relative-path io at the same time, that io ends up in
 *      the wrong place. see threading below.
 *   2. it can redirect file descriptor 1 if you ask it to (pestpp_redirect_output).
 *   3. it writes .rec / .log / .rns / .rmr files next to the control file.
 *   4. it binds a tcp port for a panther session.
 *   5. it starts child processes any time the run manager runs the forward model.
 *
 * threading
 *
 * you can have two handles in one process, but dont call them at the same time - the chdir
 * above is process wide, so two handles with different working directories cant run at once.
 * other than that, separate handles dont share anything that changes. the two exceptions are
 * the fatal-error flag (pestpp_get_fatal_error) and the fd-1 redirect, and both say so where
 * they are declared. one handle is not safe to use from two threads at once.
 *
 * error handling
 *
 * every call reports through pestpp_status. switch on the sign, not the number: negative
 * means the algorithm did something normal that isnt success, zero means success, positive
 * means it failed. new codes get added inside those bands, so check the sign.
 */
#ifndef PESTPP_API_H_
#define PESTPP_API_H_

#ifdef __cplusplus
extern "C" {
#endif

/* exported when building the library, imported when including the installed header.
 * PESTPP_API_BUILD is only defined on the pestpp_capi target. this matters on windows - if you
 * get it wrong and compile with dllexport, your binary re-exports these symbols, and for the
 * data symbol PESTPP_NAME_LEN below that is a link error. */
#if defined(_WIN32)
#  if defined(PESTPP_API_BUILD)
#    define PESTPP_API __declspec(dllexport)
#  else
#    define PESTPP_API __declspec(dllimport)
#  endif
#else
#  define PESTPP_API __attribute__((visibility("default")))
#endif

/* the version of this header, so you can compare what you built against with what
 * pestpp_get_api_version() says you loaded. bumped by hand: minor for stuff that gets added
 * without changing what is already there, major for anything else.
 *
 * while major is still 0, a breaking change might show up as a minor bump too - that is what
 * 0.x means. pin to a commit if that matters to you. */
#define PESTPP_API_VERSION_MAJOR 0
#define PESTPP_API_VERSION_MINOR 5
#define PESTPP_API_VERSION_PATCH 0

/* Opaque session handle. Owns the scenario, file manager, run manager and tool object. */
typedef void* pestpp_handle;

/* Banded by sign, so a wrapper can classify a code it has never seen:
 *
 *      < 0   the call worked; this is an ALGORITHMIC OUTCOME worth reporting
 *      = 0   the call worked
 *      > 0   the call FAILED; pestpp_last_error() says why
 *
 * Test the sign. New codes will be added inside the bands, and an `if (rc) goto fail` that
 * treats every nonzero as a failure will mis-handle the negative ones. */
typedef enum {
    PESTPP_RETRY = -1,         /* algorithmic "no upgrade taken, try again" - not a fault */

    PESTPP_OK = 0,

    PESTPP_ERROR = 1,          /* something failed; see pestpp_last_error()          */
    PESTPP_INVALID_HANDLE = 2, /* the handle is NULL, destroyed, or not a pestpp handle */
    PESTPP_INVALID_ARGUMENT = 4,   /* a NULL out-param, an unknown enum value, a bad shape */
    PESTPP_BUFFER_TOO_SMALL = 5,   /* re-query the size and try again - see below           */
    PESTPP_NOT_SUPPORTED    = 6,   /* right call, wrong tool or wrong run manager           */
    PESTPP_INVALID_STATE    = 7    /* right call, wrong moment (not initialized, and so on) */
    /* 3 is retired: it was PESTPP_RETRY before the bands existed. Left unused deliberately,
       so a stale caller comparing against the old value cannot silently match something. */
} pestpp_status;

/* PESTPP_BUFFER_TOO_SMALL is the one worth handling rather than merely reporting. The
   query-then-fill calls below tell you a count and then expect a buffer that big - but the
   count can CHANGE between the two calls, because dropping failed realizations changes the
   ensemble. A robust caller loops:

       for (;;) {
           rc = pestpp_get_phi_vector(h, t, NULL, NULL, 0, &n);   if (rc) break;
           buf = realloc(buf, n * sizeof(double));
           rc = pestpp_get_phi_vector(h, t, buf, NULL, n, &n);
           if (rc != PESTPP_BUFFER_TOO_SMALL) break;
       }

   n_out / count / needed is always written before the size check, so it is valid to read
   after a PESTPP_BUFFER_TOO_SMALL return and tells you what to allocate next. */

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
    PESTPP_SQP = 3,
    /* pestpp-glm. The odd one out: it is not an ensemble method and has no population, so it
       carries ONE parameter set through a Jacobian and an upgrade rather than an ensemble
       through a solve. The ensemble and phi-over-realizations calls therefore refuse rather
       than inventing a single-row ensemble - see pestpp_get_ensemble_view. What it does share
       with the others is the loop: initialize, iterate, finalize. */
    PESTPP_GLM = 4,
    /* pestpp-opt: sequential linear programming under chance constraints. Like glm it carries
       one decision-variable vector rather than a population, so the ensemble and
       phi-over-realizations calls refuse; unlike glm it has the chance machinery, so the stack
       and risk calls work as they do for mou and sqp. */
    PESTPP_OPT = 5
} pestpp_tool;

/* Which ensemble to look at. The raw parameter ensemble is in whatever transform space the
   tool is currently using (usually NUM during a solve) and holds only adjustable
   parameters - see pestpp_get_par_transform_status().

   PAR/OBS mean the tool's parameter and result ensembles: pe/oe for ies and da, the decision
   variable population dp and its results op for mou, dv/oe for sqp. NOISE and WEIGHTS exist
   only for ies and da; asking mou or sqp for them is an error rather than an empty array.

   STACK_PAR/STACK_OBS/NESTED_PAR are the CHANCE STACKS, and exist only for mou and sqp. They
   are empty rather than an error when the run is not stack-based - see the stack section
   below for how to tell "no stack because FOSM" from "no stack because no chance" from
   "stack not drawn yet". */
typedef enum {
    PESTPP_PAR_EN        = 0,
    PESTPP_OBS_EN        = 1,
    PESTPP_NOISE_EN      = 2,
    PESTPP_WEIGHTS_EN    = 3,
    PESTPP_STACK_PAR_EN  = 4,
    PESTPP_STACK_OBS_EN  = 5,
    PESTPP_NESTED_PAR_EN = 6
} pestpp_ensemble_id;

/* The PER-MEMBER chance stacks, addressed as PESTPP_MEMBER_STACK_EN + i for i in
   [0, pestpp_get_member_stack_count()). mou evaluates chance at several points in decision
   variable space when opt_chance_points is "all", and then each population member has its own
   observation stack; pestpp_get_member_stack_name() says which member index i belongs to.

   Like the candidate ids these are ordinary ensemble ids on purpose - views, row and column
   names and the snapshot calls all work on them unchanged.

   When chance is evaluated at a single point the count is 0 and PESTPP_STACK_OBS_EN is the
   only stack there is. */
#define PESTPP_MEMBER_STACK_EN 2000

/* The candidate ensembles of a DEFERRED SOLVE, addressed as PESTPP_CANDIDATE_EN + i for i in
   [0, pestpp_get_candidate_count()). They are ordinary ensemble ids on purpose: views, row and
   column names and the snapshot calls all work on them unchanged, so inspecting a candidate is
   the same code as inspecting the parameter ensemble.

   Valid only between pestpp_solve_prepare() and the call that consumes them. A candidate id
   outside an open deferred solve is an error, not an empty ensemble. */
#define PESTPP_CANDIDATE_EN 1000

/* Buffer sizes the library owns. Read these from the library rather than hard-coding them,
   the way xmipy and pypestutils do - and note they are DATA exports, which is why the
   dllimport half of PESTPP_API above matters.

   NAME_LEN is the width of one block in the packed name buffers used throughout. A name
   longer than this is an ERROR, not a truncation: a silently shortened name is a different
   name, and would fail to match on the way back in.
   MESSAGE_LEN is enough for anything pestpp_get_last_error() will write. */
PESTPP_API extern const int PESTPP_NAME_LEN;
PESTPP_API extern const int PESTPP_MESSAGE_LEN;

/* ---- version ------------------------------------------------------------------------
 *
 * Ask the LOADED library what it is, and compare against the PESTPP_API_VERSION_* macros the
 * caller compiled against. A caller that skips this and calls a symbol the library does not
 * export gets a loader error with no explanation, which is a bad first experience. */

/* Human readable, e.g. "5.2.28". Pass buf=NULL to learn `needed` (including the NUL). */
PESTPP_API pestpp_status pestpp_get_version(char* buf, int buf_len, int* needed);

/* The C ABI's own version, which moves independently of the pest++ release above. Any
   out-param may be NULL. */
PESTPP_API pestpp_status pestpp_get_api_version(int* major, int* minor, int* patch);

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
 * signature is either a new symbol or a break.
 *
 * HOW TO FILL IT IN, exactly:
 *
 *     pestpp_create_options opts;
 *     memset(&opts, 0, sizeof(opts));            // zero FIRST - every field defaults to 0
 *     opts.struct_size = (int)sizeof(opts);      // then declare which version you compiled
 *     opts.ctl_file = "pest.pst";
 *
 * struct_size is how the library knows which version of the struct it was handed, and it is
 * honoured rather than merely checked: a field that lies past the size you declare is one
 * your build did not have, so the library never reads it and uses the documented default
 * instead. That is what lets a binary built against this header keep working against a later
 * library that added fields.
 *
 * `reserved` exists so those later fields cannot hide in tail padding. Without it, adding an
 * int to the end of this struct would not change sizeof() on LP64 - the padding would absorb
 * it - and struct_size would silently stop distinguishing the two versions. Leave it zeroed.
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
    void*       reserved[4];     /* must be zero; room for growth outside the padding     */
} pestpp_create_options;

PESTPP_API pestpp_status pestpp_create(const pestpp_create_options* opts, pestpp_handle* out);

/* Tears the session down: the tool, the run manager, the open files, and a PANTHER master's
   listening socket. The handle is dead afterwards and passing it anywhere - including to
   pestpp_destroy() again - returns PESTPP_INVALID_HANDLE rather than crashing. Destroying
   with runs in flight abandons them. */
PESTPP_API pestpp_status pestpp_destroy(pestpp_handle h);

/* Which run manager this handle actually got. */
PESTPP_API pestpp_status pestpp_get_run_manager(pestpp_handle h, int* run_manager);

/* Why the MOST RECENT call on this handle failed; "" if it succeeded. Never returns NULL.
 *
 * READ IT IMMEDIATELY. Every entry point clears this on the way in, so the message does not
 * survive one intervening call - not even a getter:
 *
 *     if (pestpp_solve_iteration(h) > 0) {
 *         pestpp_get_iteration(h, &i);       // clears it
 *         puts(pestpp_last_error(h));        // ""
 *     }
 *
 * The returned pointer aims into storage the next call on this handle overwrites and may
 * reallocate, so copy it if you need to keep it. pestpp_get_last_error() is the copying form
 * and is the better choice for a binding; this one stays for convenience in C.
 *
 * pestpp_last_error(NULL) returns a static string and is always safe to call. */
PESTPP_API const char* pestpp_last_error(pestpp_handle h);

/* Same message, copied into a caller buffer that the caller controls the lifetime of.
   PESTPP_MESSAGE_LEN is always large enough. Pass buf=NULL to learn `needed` (with the NUL). */
PESTPP_API pestpp_status pestpp_get_last_error(pestpp_handle h, char* buf, int buf_len,
                                               int* needed);

/* Why the most recent HANDLE-LESS call failed - pestpp_create(), pestpp_redirect_output(),
   pestpp_restore_output(), pestpp_flush_output(). There is no handle to hang the message on,
   so they share one process-global slot, which the next handle-less call overwrites. */
PESTPP_API const char* pestpp_last_global_error(void);
/* Deprecated spelling of pestpp_last_global_error(); it was never create-only. */
PESTPP_API const char* pestpp_last_create_error(void);

/* ---- the fatal flag -------------------------------------------------------------------
 *
 * One failure is not confined to a handle: if the library cannot restore the working
 * directory after a call, the PROCESS is left somewhere nobody expects and every relative
 * path afterwards - in this library and in the host program - resolves to the wrong place.
 * Continuing would corrupt data quietly, so the library latches a fatal flag and every
 * subsequent call on every handle refuses with PESTPP_INVALID_STATE.
 *
 * It is process-global and it is sticky. Once you have put the working directory back
 * yourself, clear it to resume; there is no other way out. */

/* "" when healthy. Never returns NULL. */
PESTPP_API const char* pestpp_get_fatal_error(void);
/* Acknowledge and resume. Returns PESTPP_OK whether or not anything was set. */
PESTPP_API pestpp_status pestpp_clear_fatal_error(void);

/* Flush the library's console output.
 *
 * Needed because on windows the library links the static CRT (/MT), so it owns a COPY of the
 * C runtime with its own output buffer, independent of the host program's. Text the library
 * has written but not flushed sits in that private buffer and surfaces later, after a redirect
 * has been undone - landing wherever output points by then rather than in the log. Call this
 * before restoring. Harmless everywhere else. */
PESTPP_API pestpp_status pestpp_flush_output(void);

/* Send the library's console output to a file, and put it back.
 *
 * This has to happen INSIDE the library, and WHAT GETS CAPTURED DIFFERS BY PLATFORM.
 *
 * On posix it redirects file descriptor 1. There is one descriptor table per process, so that
 * catches everything the library emits and the stdout of child processes too, since the model
 * inherits the descriptor when it is spawned.
 *
 * On windows it redirects the library's C++ output stream and DOES NOT TOUCH ANY DESCRIPTOR.
 * Two consequences worth knowing before you rely on this:
 *
 *   - the model's console output is NOT captured. It never was on windows: the model is
 *     spawned with bInheritHandles=false, so it inherits nothing and writes to the console.
 *     Only the library's own text lands in the file.
 *   - it cannot damage the host. Under the static CRT (/MT, the default) this library owns a
 *     private descriptor table, so calling _dup2 on ITS descriptor 1 closed a console handle
 *     the HOST was still using - a python caller's next print died with WinError 6, and the
 *     freed handle value was recycled onto the next file opened here, which is how model
 *     output was found written inside pest.phi.composite.csv. Swapping the stream instead has
 *     no such failure mode under any runtime configuration.
 *
 * If you need the model's output on windows, capture it in the host around the model command,
 * not here.
 *
 * redirect() writes an opaque TOKEN into `redirect_token`; hand that back to restore().
 * Appends, so repeated redirects to one path accumulate rather than truncating.
 *
 * PROCESS-GLOBAL, and strictly LIFO. Output is a process-wide resource, so this is not really
 * a per-session operation however much it looks like one. Nesting is fine - each redirect
 * remembers where output was going - but they must be undone INNERMOST FIRST, and restoring
 * out of order is refused with PESTPP_INVALID_STATE rather than performed. That refusal is the
 * feature: unwinding out of order used to leave output pointing at another session's log file
 * permanently, with nothing reported.
 *
 * The token is NOT a file descriptor - the saved descriptor stays inside the library, so a
 * caller cannot close it, dup2 it, or pass a stray int that happens to be a live descriptor.
 * A token that was never issued, or has already been restored, is rejected.
 *
 * restore(0) is a no-op, so a caller can use 0 as "nothing redirected" without branching.
 * Failures report through pestpp_last_global_error(). */
PESTPP_API pestpp_status pestpp_redirect_output(const char* path, int* redirect_token);
PESTPP_API pestpp_status pestpp_restore_output(int redirect_token);

/* How many redirects are outstanding. 0 means the library is not capturing anything. Mostly
   for tests and for a host program checking it unwound everything it opened. */
PESTPP_API pestpp_status pestpp_get_redirect_depth(int* depth);

/* ---- driving the algorithm --------------------------------------------------------- */

PESTPP_API pestpp_status pestpp_initialize(pestpp_handle h);

/* initialize(), split so the caller owns the prior-ensemble evaluation.
 *
 *     int n;
 *     pestpp_initialize_prepare(h, &n);   // draws the ensembles; runs nothing
 *     if (n) { queue -> begin_batch/run_slice/end_batch -> process }
 *     pestpp_initialize_finish(h);        // phi, drop failures, write the .0. files
 *
 * `n_runs` is how many runs the caller must service. 0 is normal and needs no branch beyond
 * skipping the batch: a restart or `ies_obs_restart_csv` supplies results instead of computing
 * them, and sqp hands over nothing (its initial ensemble is evaluated deeper, inside its
 * gradient preparation, so there is no single batch to give away at this point).
 *
 * ies, da and mou all hand over their initial batch. Branch on the COUNT, never on the tool:
 * whether a given tool hands one over has already changed once - mou initialized atomically
 * until it was split - and 0 has always meant "nothing to service", never "unsupported".
 *
 * The window between prepare and queue is the useful part: it is the only point at which a
 * caller can REPLACE the drawn prior ensemble or population with its own - via
 * pestpp_set_par_snapshot() or a write through the zero-copy view - and have those
 * realizations be what gets run.
 *
 * Calling pestpp_solve_iteration() while a prepare is outstanding is an error: the tool
 * would be half-initialized, its ensembles drawn but failures and phi unprocessed. */
PESTPP_API pestpp_status pestpp_initialize_prepare(pestpp_handle h, int* n_runs);
PESTPP_API pestpp_status pestpp_initialize_finish(pestpp_handle h);

/* One iteration. Returns PESTPP_RETRY (negative - the call SUCCEEDED) where the algorithm
   rejected the upgrade and wants another attempt with a new lambda. That is an outcome, not a
   fault, which is why it is in the negative band: `if (rc) goto fail` would get it wrong. */
PESTPP_API pestpp_status pestpp_solve_iteration(pestpp_handle h);

/* ---- deferred solve --------------------------------------------------------------------
 *
 * The same iteration, split so the caller owns the model runs and can see - or change - the
 * candidates before they are evaluated:
 *
 *     pestpp_solve_prepare(h, &n);        // generate candidates, run nothing
 *     ... inspect/edit PESTPP_CANDIDATE_EN + i ...
 *     pestpp_queue_runs_subset(h, NULL, 0, &nq);   // NULL = the algorithm's own subset
 *     begin_batch -> run_slice -> end_batch
 *     pestpp_process_runs(h, &nfail);
 *     pestpp_solve_finish(h, 1, &pending);
 *     while (pending) { queue -> run -> process; pestpp_solve_finish(h, 1, &pending); }
 *
 * `n_runs` from prepare is how many runs the candidates imply. ZERO means the iteration
 * finished during preparation - a lambda that could not be generated, or the non-iterative
 * shortcut - and there is nothing to run or finish.
 *
 * ies and mou only. da's advance() is a whole noptmax loop rather than one iteration, and
 * sqp's line search issues several run batches per iteration, so neither can be expressed as
 * one generate -> run -> evaluate; both refuse rather than approximating it.
 *
 * Not available when 'ies_upgrades_in_memory' is false: the candidates are files on that
 * path, so there is nothing in memory to inspect. */
PESTPP_API pestpp_status pestpp_solve_prepare(pestpp_handle h, int* n_runs);

/* Continue after the candidate runs. `defer_runs` nonzero stops rather than running the
 * remaining realizations itself and reports how many are waiting in `pending_runs`; queue,
 * run and process them, then call again. Zero runs them internally and completes the
 * iteration in one call.
 *
 * `pending_runs` is 0 when the iteration is complete. The return value carries the same
 * PESTPP_RETRY as pestpp_solve_iteration() when the upgrade was rejected. */
PESTPP_API pestpp_status pestpp_solve_finish(pestpp_handle h, int defer_runs, int* pending_runs);

/* How many candidate ensembles are waiting. ies generates one per lambda x scale-factor
   combination; mou generates one. 0 when no deferred solve is open. */
PESTPP_API pestpp_status pestpp_get_candidate_count(pestpp_handle h, int* n);

/* The factors candidate `idx` was generated with - lambda (or mda factor) and the backtrack
   scale. Both out-params may be NULL. mou has neither and reports 0. */
PESTPP_API pestpp_status pestpp_get_candidate_info(pestpp_handle h, int idx,
                                                   double* inflation, double* backtrack);

/* Queue the outstanding batch of a deferred solve - the candidates, or the remaining
   realizations after pestpp_solve_finish() asked for them.

   `names`/`n` name the realizations to run, packed PESTPP_NAME_LEN wide like every other name
   buffer. Pass names=NULL, n=0 for the algorithm's own choice: the subset it picked for
   candidate testing, or every remaining realization. Naming realizations explicitly REPLACES
   the subset, so the realizations you do not name become the remainder. */
PESTPP_API pestpp_status pestpp_queue_runs_subset(pestpp_handle h, const char* names, int n,
                                                  int* n_queued);

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
   to assume an order.

   Pass data=NULL, row_names=NULL, col_names=NULL to learn the shape and touch nothing.

   max_nrow and max_ncol DESCRIBE THE NAME BUFFERS TOO - there is deliberately no separate
   capacity for them. row_names must hold max_nrow * PESTPP_NAME_LEN bytes and col_names
   max_ncol * PESTPP_NAME_LEN, whether or not you also pass `data`. Everything is checked
   before the first byte is written, and an undersized buffer returns
   PESTPP_BUFFER_TOO_SMALL with *nrow and *ncol telling you what to allocate.

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

/* ---- reinflation ------------------------------------------------------------------------
 *
 * Rebuild the current parameter ensemble from the PRIOR's spread, re-centred on where the
 * ensemble has got to. An ensemble smoother narrows the ensemble every iteration, and after
 * enough of them the spread can collapse to the point where there is nothing left to learn
 * from - the realizations agree with each other far more than the data justifies. Reinflation
 * is the deliberate remedy: keep the location you have fought for, put back some of the
 * variance you started with.
 *
 * The tools do this on a schedule (ies_n_iter_reinflate and friends). This is the same
 * operation as an explicit call, which is what a caller needs when the decision depends on
 * something the schedule cannot see - most obviously, reinflating AT THE MOMENT new
 * observations are brought in, so the ensemble has the spread to respond to them.
 *
 * Reinflation SELECTS realizations out of the PRIOR ensemble; it does not generate new ones.
 * The prior's size is therefore a hard ceiling, and asking for more is an error rather than
 * something quietly rounded down. To grow the ensemble mid-run, do what the executable does:
 * start with a prior big enough for the largest size you will ask for ('ies_num_reals', or a
 * prior ensemble file) and have the run begin on a subset of it - the first entry of
 * 'ies_reinflate_num_reals' truncates the working ensemble during initialization while the
 * prior keeps every realization.
 *
 *   factor      in (0, 1]. 1.0 restores the full prior spread; smaller keeps it tighter.
 *   num_reals   0 to keep the current realization count; otherwise how many realizations the
 *               reinflated ensemble should have, at most the prior's count. The SIGN selects
 *               where the spread comes from: POSITIVE uses the prior's own anomalies scaled
 *               by factor, NEGATIVE resamples the CURRENT ensemble's anomalies instead (and,
 *               when factor < 1, adds prior anomalies scaled by it on top).
 *   center_on_min_phi
 *               what the new spread is centred on. -1 follows 'ies_n_iter_reinflate' the way
 *               the built-in loop does (a negative entry there means min-phi), 0 forces the
 *               ensemble mean, 1 forces the minimum-phi realization - the aggressive form.
 *
 * ies and da only. */
PESTPP_API pestpp_status pestpp_reinflate_ensemble(pestpp_handle h, double factor,
                                                   int num_reals, int center_on_min_phi);

/* Recompute phi from the current ensembles and weights.
 *
 * Phi is CACHED - the values pestpp_get_phi_* report are whatever the last update computed,
 * which the algorithm does at its own points. Change weights or write an ensemble and the
 * cached phi is stale until this is called. */
PESTPP_API pestpp_status pestpp_update_phi(pestpp_handle h);
PESTPP_API pestpp_status pestpp_should_terminate(pestpp_handle h, int* out);

/* ---- reading ensembles ------------------------------------------------------------- */

/* Zero-copy. `data` points at the tool's live Eigen buffer, which is COLUMN-major:
   element (i,j) is data[i + j*(*nrow)].
 *
 * IT IS BORROWED. The pointer stops being yours the moment anything changes that ensemble's
 * storage - dropping or keeping realizations, an algorithm step that reallocates, destroying
 * the handle. Reading or writing through a stale pointer is a use-after-free, and the numbers
 * that come back from one look exactly like data.
 *
 * `view_token` is how you find out rather than guess. Pass a non-NULL int and the library
 * records what the view was made from; pestpp_view_is_valid() then answers whether the buffer
 * is still the ensemble's current storage - checking the address and the dimensions, not
 * merely whether the handle is alive. Pass NULL if you do not want a token, and the view is
 * untracked and entirely your problem.
 *
 * The recommended shape, and what the python layer does:
 *
 *     pestpp_get_ensemble_view(h, PESTPP_PAR_EN, &p, &nr, &nc, &tok);
 *     ... read and write p ...
 *     pestpp_release_view(h, tok);       // and stop using p
 *
 * Tokens are small and cheap, but they are held until released or until the handle is
 * destroyed, so release them. */

/* ---------------------------------------------------------------------------------------- *
 *  CHANCE STACKS - mou and sqp only
 *
 *  Both tools can account for uncertainty in the constraints and objectives, by one of two
 *  mechanisms that are configured through the same options and are NOT interchangeable:
 *
 *    FOSM    - a linear propagation through a jacobian. Cheap, no extra model runs, and the
 *              stacks below stay EMPTY.
 *    STACKS  - an ensemble of parameter realizations that is actually run. The stacks below
 *              hold it.
 *
 *  Which one is in force is derived from the options every time it is asked, not latched at
 *  initialization, so changing opt_risk or opt_stack_size through pestpp_set_option() takes
 *  effect on the next chance calculation. pestpp_get_stack_status() reports the derivation
 *  rather than the raw options, which is the difference that matters: risk is CLAMPED to
 *  [0.001, 0.999], and a stack-based run is one with stacks configured AND opt_std_weights
 *  off, neither of which is readable from any single option.
 *
 *  An empty stack is never an error - it is the correct answer for FOSM and for a risk
 *  neutral run. Ask pestpp_get_stack_status() to tell those apart from "not drawn yet". */

/* All out-params are optional; pass NULL for the ones you do not want.
 *
 * use_robust reports ROBUST optimization (opt_use_robust, pestpp-sqp only), where each
 * decision-variable realization is paired with its own parameter realization and NOTHING is
 * risk shifted - so use_chance and use_fosm are both 0 whenever use_robust is 1, and
 * opt_use_robust is refused alongside a non-neutral opt_risk. risk is the clamped, in-force
 * value, and is 0.5 exactly when chance is off. stack_size is the stack AS IT STANDS - 0 before it is drawn,
 * and not necessarily opt_stack_size, since a stack loaded from a file brings its own rows. */
PESTPP_API pestpp_status pestpp_get_stack_status(pestpp_handle h, int* use_chance, int* use_fosm,
                                                 int* use_robust, double* risk, int* stack_size);

/* How many per-member stacks exist - 0 unless opt_chance_points is "all". */
PESTPP_API pestpp_status pestpp_get_member_stack_count(pestpp_handle h, int* count);

/* The member names owning each per-member stack, packed PESTPP_NAME_LEN wide, in the same
   order as the PESTPP_MEMBER_STACK_EN + i ids. */
PESTPP_API pestpp_status pestpp_get_member_stack_names(pestpp_handle h, char* buf, int buf_len,
                                                       int* count);

PESTPP_API pestpp_status pestpp_get_ensemble_view(pestpp_handle h, int ensemble_id,
                                                  double** data, int* nrow, int* ncol,
                                                  int* view_token);

/* Nonzero when the buffer that token was issued for is still the ensemble's live storage.
   An unknown or released token is not an error - it answers 0, because it is certainly not
   valid. */
PESTPP_API pestpp_status pestpp_view_is_valid(pestpp_handle h, int view_token, int* out);

/* Give the token back. Idempotent, and harmless on a token that was never issued. */
PESTPP_API pestpp_status pestpp_release_view(pestpp_handle h, int view_token);

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

/* Set a ++ option or a * control data value by name. An unknown key is an ERROR: a set that
   quietly does nothing is the worst outcome here, because the caller goes on believing the
   option took. */
/* What the LIVE run manager is using, which is not always what the options say.
 *
 * The run manager takes its tuning values at construction. pestpp_set_option() now pushes the
 * four below onto the running manager as well, so the two agree - but reading them from the
 * manager rather than the options is how you CONFIRM that, and it is the honest answer to
 * "what will this master actually do with a run that is running long".
 *
 * All out-params are optional. The overdue policy belongs to the panther manager alone;
 * asking for it on a serial or external session is an error rather than a made-up number,
 * because nothing runs concurrently there for a run to be late against.
 *
 * These four ARE settable mid-run, including from inside a run observer - the observer
 * executes on the run loop's own thread, and the scheduling loop re-reads all four every
 * pass, so a change lands on the runs already in flight. The other panther options
 * (panther_persistent_workers, the timeout/echo/ping cadences) are consumed when the manager
 * is built and are marked init-only; setting one on a running session is refused. */
PESTPP_API pestpp_status pestpp_get_run_manager_settings(pestpp_handle h, int* max_run_fail,
                                                         double* overdue_resched_fac,
                                                         double* overdue_giveup_fac,
                                                         double* overdue_giveup_minutes);

PESTPP_API pestpp_status pestpp_set_option(pestpp_handle h, const char* key, const char* value);

/* Read one back. Pass buf=NULL to learn `needed` (including the NUL).
 *
 * `found` is how you tell "this option is set to the empty string" from "there is no such
 * option" - both used to come back as "" with PESTPP_OK, so a typo'd key read as a value.
 * It also doubles as the probe: pass a non-NULL `found` and an unknown key is answered with
 * *found = 0 and PESTPP_OK, which is what you want when asking whether this library supports
 * something. Pass NULL and an unknown key is PESTPP_INVALID_ARGUMENT instead, matching
 * pestpp_set_option - so the DEFAULT is to be told, and silence has to be asked for. */
PESTPP_API pestpp_status pestpp_get_option(pestpp_handle h, const char* key,
                                           char* buf, int buf_len, int* needed, int* found);

/* ---- run management ------------------------------------------------------------------
 *
 * A caller driving its own loop queues runs, drives the run manager in slices, and processes
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

/* Nonzero between begin_batch() and end_batch().
 *
 * The three calls above are a sequence and it is now enforced: slicing without a batch open
 * is refused rather than performed, because on PANTHER it races the master's own bookkeeping
 * and takes the process down rather than returning an error. Starting a second batch while
 * one is open is refused too - it would reset the run counters with runs in flight. */
PESTPP_API pestpp_status pestpp_is_batch_open(pestpp_handle h, int* out);

/* ---- watching a batch from inside -------------------------------------------------------
 *
 * The counters below are a POLL: they answer when you ask. This is the push side, and it is
 * the only way to see anything during pestpp_solve_iteration(), which owns its runs and does
 * not return until they are done.
 *
 * The observer is also the only place a caller sees runs WHILE THEY ARE IN FLIGHT, which is
 * why it returns an action rather than void. Stopping a batch is the action that exists
 * today; preemption - ask the workers what they have, keep the runs worth finishing - is a
 * new RETURN VALUE here rather than a new callback, and `struct_size` is what lets the struct
 * gain the fields that will need (see docs/api_part1/panther_preemption.md). Both would have
 * been ABI breaks against a void observer taking a bare struct. */

typedef struct {
    int    struct_size;      /* = sizeof(pestpp_run_progress); required, and checked */
    int    n_total;          /* runs in this batch */
    int    n_completed;
    int    n_failed;
    int    n_timed_out;
    int    n_running;
    int    run_id;           /* the run this is about, or -1 for a periodic tick */
    double elapsed_sec;
} pestpp_run_progress;

/* What the observer asks for next. Unknown values are treated as CONTINUE, so a caller built
   against a later header cannot wedge an older library. */
typedef enum {
    PESTPP_RUN_CONTINUE   = 0,
    PESTPP_RUN_STOP_BATCH = 1
    /* 2 reserved: request partial results from the runs still going */
} pestpp_run_action;

typedef int (*pestpp_run_observer_fn)(const pestpp_run_progress* progress, void* user_data);

/* Observe every batch this session runs. Pass fn=NULL to stop observing.
 *
 * `min_interval_sec` throttles inside the library, because an observer cannot decline a call
 * it has already been handed and a serial batch of thousands of quick runs would otherwise
 * pay a cross-ABI call for each one. A run that finishes, fails or times out is reported
 * whatever the interval - throttling drops periodic ticks, never events.
 *
 * Three rules, and the first is the one that bites:
 *
 * 1. The observer is called with the tool's state MID-BATCH. Only the run-management calls
 *    are legal from inside it - pestpp_get_run_states, pestpp_get_run_time_stats,
 *    pestpp_cancel_runs, pestpp_release_workers, pestpp_get_worker_*. Anything touching
 *    ensembles, phi or options is refused with PESTPP_INVALID_STATE rather than allowed to
 *    read a half-updated ensemble.
 *    That allowlist is exactly what preemption needs, which is why it is a list and not a ban.
 * 2. It is called on the thread that called into the library, never a worker thread.
 * 3. It must not throw across the boundary. One that does is unregistered rather than being
 *    allowed to unwind through a batch. */
PESTPP_API pestpp_status pestpp_set_run_observer(pestpp_handle h, pestpp_run_observer_fn fn,
                                                 void* user_data, double min_interval_sec);

/* Ask the workers running these runs to report whatever results they have RIGHT NOW.
 *
 * Asynchronous: this returns as soon as the requests are sent. A worker answers when it can,
 * or not at all - nothing blocks and no run is interrupted. `n_requested` receives how many
 * requests actually went out, which can be fewer than asked for: a request is only sent to an
 * agent that advertised support for it, because an older agent treats an unrecognised message
 * during a run as corrupt and TERMINATES the run rather than ignoring it.
 *
 * Pass run_ids=NULL, n=0 for every run currently executing.
 *
 * PANTHER only - there is nobody to ask on the serial run manager, which reports 0 rather
 * than failing. Legal from inside a run observer. */
PESTPP_API pestpp_status pestpp_request_partial_results(pestpp_handle h, const int* run_ids,
                                                        int n, int* n_requested);


/* ---- partial results ---------------------------------------------------------------------
 *
 * Two views of the same idea, and the difference matters. The first is LIVE: what a worker has
 * reported for a run in the batch that is executing now. The second reads it back out of run
 * storage, and works after the fact.
 *
 * Neither reports partial-ness as a run STATUS. The status byte has four meanings and three of
 * them are load-bearing in the failed-run, good-run and restart logic; a fifth would have to be
 * handled correctly at every one of those sites, silently wrong if it were not. Partial-ness is
 * a property of the VALUES, so it is reported alongside them. */

/* Live: has a worker reported partial results for this run in the current batch?
 *
 * `has_partial` is 0/1; the counts say how much of the run is real. PANTHER only - the serial
 * manager has no worker to ask - and it reports has_partial=0 rather than failing. Cleared when
 * a new batch begins, because a partial result describes a process that is still running. */
PESTPP_API pestpp_status pestpp_get_run_partial_info(pestpp_handle h, int run_id,
                                                     int* has_partial, int* n_obs_reported,
                                                     int* n_obs_total);

/* How complete a run's stored VALUES are - see pestpp_get_run_values(). */
typedef enum {
    PESTPP_RUN_VALUES_NONE    = 0,  /* nothing recorded yet                            */
    PESTPP_RUN_VALUES_PARTIAL = 1,  /* preemption wrote what a worker could parse      */
    PESTPP_RUN_VALUES_FINAL   = 2   /* the run finished and these are its real results */
} pestpp_run_completeness;

/* How many runs are in storage. */
PESTPP_API pestpp_status pestpp_get_run_count(pestpp_handle h, int* n);

/* Status and completeness for one stored run; any out-param may be NULL.
 *
 * `completeness` is DERIVED from the values rather than stored: a completed run is FINAL, a
 * run that is not complete but has at least one real observation is PARTIAL, and one with
 * none is NONE. That keeps it truthful after a restart, when the live bookkeeping is gone and
 * ought to be - a partial result from a process that is no longer running is not meaningful. */
PESTPP_API pestpp_status pestpp_get_run_info(pestpp_handle h, int run_id, int* status,
                                             int* completeness, int* n_obs_reported,
                                             int* n_obs_total);

/* Values for one stored run. Pass npars/nobs of 0 with NULL arrays to learn the sizes.
 *
 * `obs_valid` is the point of this signature: one byte per observation, non-zero where the
 * value is real. A partial result fills the rest with the no-data sentinel, and a caller must
 * never have to recognise that sentinel by eye - that is how a magic number ends up in
 * somebody's phi. It may be NULL, and it is filled for completed runs too (all non-zero), so
 * one code path serves both. */
PESTPP_API pestpp_status pestpp_get_run_values(pestpp_handle h, int run_id,
                                               double* pars, int npars,
                                               double* obs, int nobs,
                                               unsigned char* obs_valid,
                                               int* npars_out, int* nobs_out);

/* ---- servicing runs yourself ----------------------------------------------------------
 *
 * The mirror of pestpp_get_run_values(): evaluate the queued runs however you like and write
 * the results back, instead of handing the batch to a run manager.
 *
 *     pestpp_queue_runs(h, &n);                       // parameters land in run storage
 *     for (id = 0; id < n; id++) {
 *         pestpp_get_run_values(h, id, pars, np, NULL, 0, NULL, NULL, NULL);
 *         ... evaluate, in this process, on a cluster, wherever ...
 *         pestpp_set_run_values(h, id, obs, nobs);    // or pestpp_set_run_failed(h, id)
 *     }
 *     pestpp_process_runs(h, &n_failed);              // into the ensembles, as usual
 *
 * This is what the EXTERNAL run manager does, without the round trip: /e writes the .rns, a
 * script fills in observations, pest++ reads them back. From the API that detour is pointless
 * when the caller is already in the process.
 *
 * The parameters you read back are in MODEL space - transformed, exactly as a forward run
 * would receive them - because that is what was written when the runs were queued.
 *
 * Use INSTEAD of driving a run manager over the same batch, not alongside one. Both calls are
 * refused while a batch is open, because a run manager writing results for a run you are also
 * writing has no defined outcome.
 *
 * A run you never write stays incomplete, and pestpp_process_runs() drops its realization
 * rather than inventing one - the observation block of an unwritten record is the no-data
 * sentinel throughout, not zeros. Silence is therefore safe, if wasteful. */

/* The run-storage name order, which is what pestpp_get_run_values() returns and what
 * pestpp_set_run_values() expects. Packed as fixed-width PESTPP_NAME_LEN blocks, like the
 * ensemble name calls; pass buf=NULL to query the count.
 *
 * Needed because run storage does NOT share the ensemble's column order - an ensemble may be
 * ordered differently and can hold a subset - so the positional arrays either side of this
 * interface are meaningless without it. Ask here, not from pestpp_get_ensemble_col_names(). */
PESTPP_API pestpp_status pestpp_get_run_par_names(pestpp_handle h, char* buf, int buf_len,
                                                  int* count);
PESTPP_API pestpp_status pestpp_get_run_obs_names(pestpp_handle h, char* buf, int buf_len,
                                                  int* count);

/* Record this run's observations and mark it complete. `nobs` must be the full observation
   count, in pestpp_get_run_obs_names() order. The run's stored parameters are left
   untouched. */
PESTPP_API pestpp_status pestpp_set_run_values(pestpp_handle h, int run_id,
                                               const double* obs, int nobs);

/* Mark this run failed - the model did not produce usable output. Distinct from simply not
   writing it: this is recorded as a failure, and counts as one wherever failures are counted. */
PESTPP_API pestpp_status pestpp_set_run_failed(pestpp_handle h, int run_id);

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

/* Hand workers back so their compute can go elsewhere. Indices address the same ordering
   pestpp_get_worker_state() uses, so release(idx) frees the worker get_worker_state(idx) just
   described. Pass worker_idxs=NULL with n=0 to release every worker. The count actually
   released is written to n_released; asking for 8 and getting 3 is visible rather than silent.

   A BUSY worker is released too, and its run is RESCHEDULED - put back at the front of the
   queue for another worker to pick up. It is not failed, not cancelled, and nothing is charged
   against max_n_failure. Use pestpp_cancel_runs() when the judgement is about the run; use
   this when it is about the machine.

   Legal from inside a run observer, like the other worker calls. PANTHER only. */
PESTPP_API pestpp_status pestpp_release_workers(pestpp_handle h, const int* worker_idxs, int n,
                                                int* n_released);

/* ---- Jacobian and the single parameter vector (glm) ------------------------------------
   The ensemble tools refuse these: they carry a population, not one parameter set, and never
   form the Jacobian. Use pestpp_get_par_snapshot() and the ensemble views for those. */

/* Which parameter vector pestpp_get_par_vector() should return. */
typedef enum {
    PESTPP_PAR_CURRENT = 0,   /* what the next iteration will work from */
    PESTPP_PAR_OPTIMUM = 1    /* the best found so far */
} pestpp_par_vector;

/* The Jacobian as a DENSE COLUMN-MAJOR COPY, rows = simulated observations, columns = the
   adjustable parameters it was built over. Call with all pointers NULL to learn the shape,
   size the buffers, then call again. Not a view: the matrix is sparse internally, so there is
   no contiguous buffer to hand out. Available only after initialize() has built it. */
PESTPP_API pestpp_status pestpp_get_jacobian(pestpp_handle h, double* data,
                                             int max_nrow, int max_ncol,
                                             int* nrow, int* ncol,
                                             char* row_names, char* col_names);

/* ---- pestpp-da: running the assimilation cycles -------------------------------------------
   a DataAssimilator only handles one cycle. what makes it data assimilation is the sequence:
   every cycle gets its own child pest object and its own tool, and the posterior from cycle N
   becomes the prior for cycle N+1. any realizations a cycle loses get dropped from the global
   ensembles too, so the two stay lined up.

   without these calls a session only runs the one cycle it was made with, so a control file
   with several cycles would quietly give you a one cycle answer.

   two ways to run it. pestpp_da_run_all_cycles() does the whole thing, same as the exe. or you
   can step through it:

       pestpp_da_cycles_initialize(h);
       int begun;
       while (pestpp_da_cycle_begin(h, &begun) == PESTPP_OK && begun) {
           pestpp_da_cycle_drive(h);      // or drive the cycle yourself, see below
           pestpp_da_cycle_end(h);        // NOT optional - this is what carries the cycle forward
       }

   in between begin and end, every other call on the session is about the cycle that is open -
   its ensembles, its parameter and observation names, its phi. so instead of calling
   pestpp_da_cycle_drive() you can drive the cycle yourself with the normal calls and get
   iteration by iteration control inside each cycle. either way you still have to call
   pestpp_da_cycle_end().

   all of these give an error for any tool other than da. */

/* build the global ensembles and work out the cycle list. calling it twice doesnt rebuild
   anything - the global ensembles are the state of the run, so rebuilding them partway through
   would throw away everything assimilated so far. */
PESTPP_API pestpp_status pestpp_da_cycles_initialize(pestpp_handle h);

/* how many cycles are in the control file. */
PESTPP_API pestpp_status pestpp_da_n_cycles(pestpp_handle h, int* n);

/* the cycle numbers in the order they run. `cycles` needs room for pestpp_da_n_cycles() ints. */
PESTPP_API pestpp_status pestpp_da_get_cycles(pestpp_handle h, int* cycles, int n);

/* start the next cycle. *begun is 1 if a cycle is now open and 0 if you are done - being done
   isnt an error, so the status is PESTPP_OK either way and you loop on *begun. cycles before
   ++da_hotstart_cycle get skipped here, same as the exe skips them. */
PESTPP_API pestpp_status pestpp_da_cycle_begin(pestpp_handle h, int* begun);

/* the cycle that is open, or -1 if none is. */
PESTPP_API pestpp_status pestpp_da_current_cycle(pestpp_handle h, int* cycle);

/* run the open cycle - initialize the tool, assimilate, report phi. */
PESTPP_API pestpp_status pestpp_da_cycle_drive(pestpp_handle h);

/* pull the open cycle's posterior back into the global ensembles and close the cycle. */
PESTPP_API pestpp_status pestpp_da_cycle_end(pestpp_handle h);

/* run every cycle, start to finish. */
PESTPP_API pestpp_status pestpp_da_run_all_cycles(pestpp_handle h);

/* ---- constraints --------------------------------------------------------------------------
   mou, sqp and opt all carry constraints, so these are shared instead of being opt-only.
   "what is the optimum" and "which constraints are binding there" are two halves of the same
   answer, and until now you could only get the first half.

   the simulated constraint VALUES come back through pestpp_get_obs_vector() - they are just the
   simulated observations at the current point. */

/* the constraint names, obs constraints first then prior information ones.
   call with buf=NULL to get the count. */
PESTPP_API pestpp_status pestpp_get_constraint_names(pestpp_handle h, char* buf, int buf_len,
                                                     int* count);

/* the sense of each one (less_than / greater_than / equal_to), same order as the names. */
PESTPP_API pestpp_status pestpp_get_constraint_senses(pestpp_handle h, char* buf, int buf_len,
                                                      int* count);

/* how far the current point violates its constraints, added up. 0 means feasible. */
PESTPP_API pestpp_status pestpp_get_sum_of_violations(pestpp_handle h, double* total);

/* the simulated observations at the current point, for tools that carry one set rather than an
   ensemble of them (glm and opt). for opt these are the constraint values. `vals` and `names`
   may each be NULL - pass both NULL with a non-NULL count to size the buffers first. */
PESTPP_API pestpp_status pestpp_get_obs_vector(pestpp_handle h, double* vals, int max_n,
                                               char* names, int names_len, int* count);

/* ---- pestpp-opt: objective values and decision variables ----------------------------------
   opt has one objective function value that changes from iteration to iteration. no other tool
   works that way - the ensemble methods report phi over realizations and mou reports a pareto
   front - so opt gets its own calls instead of trying to squeeze this into the phi calls.

   the decision variable VALUES come back through pestpp_get_par_vector(), since they are opt's
   parameter vector and a second call would just be another name for the same thing. what does
   need its own call is WHICH parameters are decision variables, because when you turn on
   chance constraints the adjustable set is bigger than the decision variable set.

   all of these give an error for any tool other than opt. */

/* how many objective values there are - also how big to make the array below. */
PESTPP_API pestpp_status pestpp_opt_n_objectives(pestpp_handle h, int* n);

/* the objective value at each iteration, oldest first. */
PESTPP_API pestpp_status pestpp_opt_get_objectives(pestpp_handle h, double* vals, int max_n);

/* the best objective value seen, and the one it started from. either can be NULL. */
PESTPP_API pestpp_status pestpp_opt_objective_bounds(pestpp_handle h, double* best,
                                                     double* initial);

/* the decision variable names. call with buf=NULL first to get the count. */
PESTPP_API pestpp_status pestpp_opt_dec_var_names(pestpp_handle h, char* buf, int buf_len,
                                                  int* count);

/* ---- Tikhonov regularization (glm) -------------------------------------------------------
   The dynamically adjusted regularization weight and the dials that govern it. These arrive
   from the control file's `* regularization` section rather than as ++ options, which is why
   pestpp_set_option() does not reach them.

   The ensemble tools refuse both calls: they regularize through the prior ensemble and have no
   weight to read or set. */

/* Read the whole state. Any out-param may be NULL. `use_dynamic` is 0/1. */
PESTPP_API pestpp_status pestpp_get_regularization(pestpp_handle h, int* use_dynamic,
                                                   double* weight, double* phimlim,
                                                   double* phimaccept, double* fracphim,
                                                   double* wfmin, double* wfmax, double* wffac,
                                                   double* wftol, double* wfinit,
                                                   int* max_reg_iter, int* adj_grp_weights);

/* Set ONE dial by name: use_dynamic_reg, weight, phimlim, phimaccept, fracphim, wfmin, wfmax,
   wffac, wftol, wfinit, max_reg_iter, adj_grp_weights (accepted as `iregadj` too, which is
   the control-file spelling). Keyed rather than wide so a caller changes one without
   restating the rest, and so an unknown name is an error rather than a silently dropped
   argument. Live - the solver reads these at point of use, so a change takes effect on the
   next iteration. */
PESTPP_API pestpp_status pestpp_set_regularization(pestpp_handle h, const char* key,
                                                   double value);

/* The upgraded parameters a given LAMBDA would produce, WITHOUT running anything. glm only.

   The query the lambda loop makes internally, exposed so a caller can make it directly and
   compare lambdas before spending runs on them. Needs a Jacobian, so call it after an
   iteration or after the jacobian stages below.

   Values come back in CTL space, bound- and change-limited exactly as the solver's own
   candidates are; `limit_buf` receives "none", "bound", "factor" or "relative", because a
   truncated upgrade is a different answer from an untruncated one. Names are sorted and packed
   PESTPP_NAME_LEN-wide, matching pestpp_get_par_vector(). NULL pointers size it.

   Does NOT adjust the regularization weight - see the note on SVDSolver::compute_upgrade. */
PESTPP_API pestpp_status pestpp_compute_upgrade(pestpp_handle h, double lambda,
                                                double* vals, int max_n, int* n, char* names,
                                                char* limit_buf, int limit_buf_len);

/* The Jacobian batch as three stages: queue, run, harvest. glm only.

   prepare() builds the finite-difference runs and QUEUES them with the run manager without
   running them, returning how many. Run them however you like - pestpp_jacobian_run() is the
   in-tree way, but driving the run manager yourself is the point of the split - then
   process() harvests them into the Jacobian and writes the jco/sen output.

   The tool's own iteration composes exactly these three calls, so there is one code path. */
PESTPP_API pestpp_status pestpp_jacobian_prepare(pestpp_handle h, int calc_init_obs, int* n_runs);
PESTPP_API pestpp_status pestpp_jacobian_run(pestpp_handle h);
PESTPP_API pestpp_status pestpp_jacobian_process(pestpp_handle h, int* ok);

/* The tool's parameter vector in CTL space, names packed PESTPP_NAME_LEN-wide and sorted, so
   names and values always correspond. NULL pointers size it, as above. */
PESTPP_API pestpp_status pestpp_get_par_vector(pestpp_handle h, int which, double* vals,
                                               int max_n, int* n, char* names);

/* Push values back onto the CURRENT vector, matched BY NAME - so a partial vector updates only
   what it names and leaves the rest alone. A name the tool does not hold is an error, not a
   silent skip: it means the caller and the tool disagree about the parameter set. */
PESTPP_API pestpp_status pestpp_set_par_vector(pestpp_handle h, const double* vals, int n,
                                               const char* names);

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

/* ---- queue / process -----------------------------------------------------------------
 *
 * The two halves of an evaluation, so a caller can own what happens in between: queue the
 * current parameter ensemble, drive and watch the run manager, then process into the
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
 * Call with data=NULL, row_names=NULL, col_names=NULL to learn nrow/ncol, then size and
 * refill. As with pestpp_get_phi_residuals(), max_nrow and max_ncol bound the NAME buffers as
 * well as `data` - row_names needs max_nrow * PESTPP_NAME_LEN bytes and col_names
 * max_ncol * PESTPP_NAME_LEN. Undersized returns PESTPP_BUFFER_TOO_SMALL before writing. */
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

#endif /* PESTPP_API_H_ */

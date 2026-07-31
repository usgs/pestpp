/*


	This file is part of PEST++.

	PEST++ is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	PEST++ is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with PEST++.  If not, see<http://www.gnu.org/licenses/>.
*/
#ifndef RUNMANAGERABSTRACT_H
#define RUNMANAGERABSTRACT_H

#include <string>
#include <vector>
#include <set>
#include <map>
#include "RunStorage.h"
#include <Eigen/Dense>
#include <chrono>
#include <functional>



class ModelExecInfo;
class Parameters;
class Observations;


/**
 * @brief Abstract base class for all PEST++ run managers.
 *
 * Manages the lifecycle of model runs including parameter distribution,
 * run execution, result collection, and failure tracking.  Three RunStorage
 * instances are maintained: \c file_stor for the primary run storage file
 * (*case.rns*), \c failed_file_stor for the failed-run storage file
 * (*case.rnf*), and \c all_file_stor for the optional "all runs" storage file
 * (*case.allruns.rns*).  The failed-run file persists after a successful run so
 * that users can inspect which parameter combinations caused model failure.
 * The all-runs file, enabled via the SAVE_ALL_RUNS control-file option, records
 * every completed model run (both successes and failures) across all iterations
 * and also persists after a graceful exit.
 */
/**
 * @brief What a progress observer is told, and what it can ask for in return.
 *
 * Deliberately a plain value with no references into run-manager internals: it is copied
 * across the C ABI to a caller that may be in another language, and it is emitted from inside
 * a batch, when very little else is in a consistent state.
 */
struct RunProgress
{
	int n_total = 0;        ///< runs in this batch
	int n_completed = 0;
	int n_failed = 0;
	int n_timed_out = 0;
	int n_running = 0;
	/// the run this notification is about, or -1 for a periodic tick with no particular run
	int run_id = -1;
	double elapsed_sec = 0.0;
};

/**
 * @brief What an observer asks the run manager to do next.
 *
 * The observer is the only place a caller sees runs WHILE THEY ARE IN FLIGHT, so it is also
 * the only place a decision about them can be made in time to matter. That is why this is an
 * action rather than void: preemption - ask the workers what they have, kill the runs that
 * are not worth finishing (see docs/api_part1/panther_preemption.md) - is a new value here,
 * not a new callback. A void observer would have made that an ABI break.
 */
enum class RunAction
{
	CONTINUE = 0,    ///< carry on
	STOP_BATCH = 1   ///< stop scheduling new runs and bring the batch to an orderly end
	// REQUEST_PARTIAL = 2 - reserved for preemption
};

class RunManagerAbstract
{
public:
	enum class RUN_UNTIL_COND { NORMAL, NO_OPS, TIME, NO_OPS_OR_TIME };
	enum class RUN_MGR_TYPE {NOTDEFINED, PANTHER, SERIAL};

	/// Observe a batch as it runs. Pass an empty function to stop observing.
	///
	/// `min_interval_sec` throttles at the SOURCE rather than in the observer: a serial batch
	/// of thousands of sub-second runs would otherwise pay for a cross-ABI call per run, and
	/// an observer that only draws a bar cannot decline the call it has already been handed.
	/// A run that finishes, fails or times out is always reported, whatever the interval -
	/// throttling may drop periodic ticks, never events.
	void set_progress_observer(std::function<RunAction(const RunProgress&)> fn,
		double min_interval_sec = 0.0)
	{
		progress_fn = fn;
		progress_min_interval = min_interval_sec;
		progress_last = std::chrono::system_clock::time_point::min();
		progress_stop_requested = false;
	}
	/// Set once an observer has asked for STOP_BATCH; the run managers check it in their
	/// scheduling loops. Cleared by the next set_progress_observer() or begin_batch().
	bool progress_stop_requested = false;
	RunManagerAbstract(const std::vector<std::string> _comline_vec,
		const std::vector<std::string> _tplfile_vec, const std::vector<std::string> _inpfile_vec,
		const std::vector<std::string> _insfile_vec, const std::vector<std::string> _outfile_vec,
		const std::string &stor_filename, int _max_n_failure=1);
	virtual void initialize(const std::vector<std::string> &model_par_names, std::vector<std::string> &obs_names, const std::string &_filename = std::string(""));
	virtual void initialize(const Parameters &model_pars, const Observations &obs, const std::string &_filename = std::string(""));
	virtual void initialize_restart(const std::string &_filename);
	virtual void reinitialize(const std::string &_filename = std::string(""));
	virtual void free_memory();
	virtual int add_run(const Parameters &model_pars, const std::string &info_txt="", double info_value=RunStorage::no_data);
	virtual int add_run(const std::vector<double> &model_pars, const std::string &info_txt="", double info_valuee=RunStorage::no_data);
	virtual int add_run(const Eigen::VectorXd &model_pars, const std::string &info_txt="", double info_valuee=RunStorage::no_data);
	virtual void update_run(int run_id, const Parameters &pars, const Observations &obs);
	virtual void run() = 0;
	virtual RunManagerAbstract::RUN_UNTIL_COND run_until(RUN_UNTIL_COND condition, int n_nops = 0, double sec = 0.0);

	/// Whether a run_slice() slice left work outstanding.
	enum class RUN_SLICE_STATUS { RUNNING, ALL_DONE };

	// Sliced batch driving: begin_batch(), run_slice() until ALL_DONE, end_batch(). Lets a caller
	// - the tools here, or an API user - get control back between slices to watch progress,
	// query run states or cancel runs. Only PANTHER can actually yield mid-batch; the
	// defaults below run the whole batch in the first run_slice(), so every run manager can be
	// driven the same way and callers do not have to special-case.
	virtual void begin_batch() {}
	virtual RUN_SLICE_STATUS run_slice(double max_seconds = 0.05) { run(); return RUN_SLICE_STATUS::ALL_DONE; }
	virtual void end_batch(RUN_UNTIL_COND terminate_reason = RUN_UNTIL_COND::NORMAL) {}
	virtual const std::vector<std::string> &get_par_name_vec() const;
	virtual const std::vector<std::string> &get_obs_name_vec() const;
	virtual void get_info(int run_id, int &run_status, std::string &info_txt, double &info_value);
	virtual bool run_finished(int run_id);
	virtual bool get_run(int run_id, Parameters &pars, Observations &obs, bool clear_old=true);
	virtual bool get_run(int run_id, Parameters &pars, Observations &obs, std::string &info_txt, double &info_value, bool clear_old=true);
	virtual bool get_run(int run_id, double *pars, size_t npars, double *obs, size_t nobs, std::string &info_txt, double &info_value);
	virtual bool get_run(int run_id, double *pars, size_t npars, double *obs, size_t nobs);
	virtual bool get_run(int run_id, std::vector<double> &pars_vec, std::vector<double> &obs_vec, std::string &info_txt, double &info_value);
	virtual bool get_run(int run_id, std::vector<double> &pars_vec, std::vector<double> &obs_vec);
	/** @brief Return the set of run IDs whose failure count has exceeded max_run_fail. */
	virtual const std::set<int> get_failed_run_ids();
    virtual const std::map<std::string,std::vector<int>> get_run_info_map();
	virtual bool get_model_parameters(int run_num, Parameters &pars);
	virtual bool get_observations_vec(int run_id, std::vector<double> &data_vec);
	virtual Observations get_obs_template(double value = -9999.0) const;
	virtual int get_total_runs(void) const {return total_runs;}
	virtual int get_num_good_runs(void);
	/** @brief Return the number of runs whose failure count has exceeded max_run_fail. */
	virtual int get_num_failed_runs(void);
	/** @brief Return true if the run identified by @p id has failed more than max_n_failure times. */
	virtual bool n_run_failures_exceeded(int id);
	virtual int get_nruns(void) {return file_stor.get_nruns();}
	virtual int get_cur_groupid(void);
	virtual std::vector<int> get_outstanding_run_ids();
	virtual ~RunManagerAbstract(void) {}
	virtual std::string get_run_filename() { return file_stor.get_filename(); }
	/** @brief Return a const reference to the primary run storage (*case.rns*). */
	virtual const RunStorage& get_runstorage_ref() const;
	/** @brief Return a const reference to the failed-run storage (*case.rnf*). */
	virtual const RunStorage& get_failed_runstorage_ref() const;
	/** @brief Return the number of runs recorded in the failed-run storage file. */
	virtual int get_num_failed_stored() { return failed_file_stor.get_nruns(); }
	/** @brief Return the filename of the failed-run storage file (*case.rnf*). */
	virtual std::string get_failed_run_filename() { return failed_file_stor.get_filename(); }
	/** @brief Enable/disable recording every completed run to the all-runs storage file (*case.allruns.rns*). */
	virtual void set_save_all_runs(bool _save_all_runs) { save_all_runs = _save_all_runs; }
	/** @brief Return true if all completed runs are being recorded to the all-runs storage file. */
	virtual bool get_save_all_runs() const { return save_all_runs; }
	/** @brief Return a const reference to the all-runs storage (*case.allruns.rns*). */
	virtual const RunStorage& get_all_runstorage_ref() const;
	/** @brief Return the number of runs recorded in the all-runs storage file. */
	virtual int get_num_all_stored() { return all_file_stor.get_nruns(); }
	/** @brief Return the filename of the all-runs storage file (*case.allruns.rns*). */
	virtual std::string get_all_run_filename() { return all_file_stor.get_filename(); }
	virtual void print_run_summary(std::ostream &fout) { file_stor.print_run_summary(fout); }
	//virtual Observations get_init_run_obs() { return init_run_obs; }
	virtual std::vector<double> get_init_sim() { return init_sim;  }
	virtual void set_init_sim(std::vector<double> _init_sim) { init_sim = _init_sim; }
	virtual RUN_MGR_TYPE get_mgr_type() { return mgr_type; }

protected:
	std::function<RunAction(const RunProgress&)> progress_fn;
	double progress_min_interval = 0.0;
	std::chrono::system_clock::time_point progress_last =
		std::chrono::system_clock::time_point::min();

	/**
	 * @brief Tell the observer where the batch has got to.
	 *
	 * Called from the run managers' own scheduling loops. `is_event` marks a run that
	 * finished, failed or timed out - those are never dropped by the throttle, because a
	 * caller deciding whether to keep waiting must not miss the thing it is waiting for.
	 *
	 * MUST be called only from the thread that entered the run manager. Both scheduling
	 * loops run there, so this costs nothing to honour - and it is what keeps a cross-ABI
	 * observer out of a foreign thread, where taking python's GIL while the calling thread
	 * sits blocked inside the library is a deadlock rather than an error.
	 */
	void notify_progress(const RunProgress& p, bool is_event = false)
	{
		if (!progress_fn)
			return;
		std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
		if ((!is_event) && (progress_min_interval > 0.0) &&
			(progress_last != std::chrono::system_clock::time_point::min()))
		{
			std::chrono::duration<double> since = now - progress_last;
			if (since.count() < progress_min_interval)
				return;
		}
		progress_last = now;
		// an observer that throws must not unwind through a batch: a C ABI observer cannot
		// let an exception cross the boundary at all, and a C++ one has no business killing
		// runs that are already in flight
		try
		{
			if (progress_fn(p) == RunAction::STOP_BATCH)
				progress_stop_requested = true;
		}
		catch (...)
		{
			progress_fn = nullptr;   // it failed once; do not keep calling it
		}
	}

	int total_runs;
	int max_n_failure; // maximum number of times to retry a failed model run
	int cur_group_id;  // used in some of the derived classes (ie PANTHER)
	RunStorage file_stor;          ///< Primary run storage (*case.rns*).
	RunStorage failed_file_stor;   ///< Failed-run storage (*case.rnf*). Persists after graceful exit.
	RunStorage all_file_stor;      ///< Optional all-runs storage (*case.allruns.rns*). Persists after graceful exit.
	bool save_all_runs = false;    ///< When true, every completed run is copied into all_file_stor.
	RUN_MGR_TYPE mgr_type;
	std::vector<std::string> comline_vec;
	std::vector<std::string> tplfile_vec;
	std::vector<std::string> inpfile_vec;
	std::vector<std::string> insfile_vec;
	std::vector<std::string> outfile_vec;
	bool run_requried(int run_id);
	//Observations init_run_obs;
	std::vector<double> init_sim;
	/**
	 * @brief Mark a run as failed and, if the failure threshold is exceeded,
	 *        record its parameters in the failed-run storage file (*case.rnf*).
	 *
	 * @param run_id  The run manager ID of the failed run.
	 */
	virtual void update_run_failed(int run_id);
	/**
	 * @brief If SAVE_ALL_RUNS is enabled, copy a completed run (its parameters,
	 *        observations, metadata and status) from the primary storage into the
	 *        all-runs storage file (*case.allruns.rns*).  No-op when disabled.
	 *
	 * @param run_id  The run manager ID of the completed run.
	 */
	void record_all_run(int run_id);
};

#endif /*  RUNMANAGERABSTRACT_H */

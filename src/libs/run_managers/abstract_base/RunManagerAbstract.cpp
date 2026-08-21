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
/**
 * @file RunManagerAbstract.cpp
 * @brief Implementation of RunManagerAbstract.
 */

#include "RunManagerAbstract.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <iterator>
#include <cassert>
#include <cstring>
#include "Transformable.h"
#include "utilities.h"

/**
 * @brief Derive the failed-run storage filename from the primary storage filename.
 *
 * Replaces the *.rns* extension with *.rnf*, or appends *.rnf* if no
 * *.rns* extension is found.
 *
 * @param stor_filename  The primary run storage filename (*case.rns*).
 * @return The corresponding failed-run storage filename (*case.rnf*).
 */
static string make_failed_filename(const string& stor_filename) {
	string fn = stor_filename;
	size_t pos = fn.rfind(".rns");
	if (pos != string::npos) fn.replace(pos, 4, ".rnf");
	else fn += ".rnf";
	return fn;
}

/**
 * @brief Derive the all-runs storage filename from the primary storage filename.
 *
 * Replaces the *.rns* extension with *.allruns.rns*, or appends *.allruns.rns*
 * if no *.rns* extension is found.  The resulting name is distinct from the
 * primary *case.rns* file, so it is not removed by the run-storage cleanup.
 *
 * @param stor_filename  The primary run storage filename (*case.rns*).
 * @return The corresponding all-runs storage filename (*case.allruns.rns*).
 */
static string make_all_filename(const string& stor_filename) {
	string fn = stor_filename;
	size_t pos = fn.rfind(".rns");
	if (pos != string::npos) fn.replace(pos, 4, ".allruns.rns");
	else fn += ".allruns.rns";
	return fn;
}

RunManagerAbstract::RunManagerAbstract(const vector<string> _comline_vec,
	const vector<string> _tplfile_vec, const vector<string> _inpfile_vec,
	const vector<string> _insfile_vec, const vector<string> _outfile_vec,
	const string &stor_filename, int _max_n_failure)
  : total_runs(0), max_n_failure(_max_n_failure), file_stor(stor_filename),
    failed_file_stor(make_failed_filename(stor_filename)),
    all_file_stor(make_all_filename(stor_filename)),
    comline_vec(_comline_vec), tplfile_vec(_tplfile_vec),
    inpfile_vec(_inpfile_vec), insfile_vec(_insfile_vec), outfile_vec(_outfile_vec)
{
	/*cout << endl;
	cout << "             Generalized Run Manager Interface" << endl;
	cout << "                             by:" << endl << endl;
	cout << "               The PEST++ Development Team" << endl;

	cout << endl << endl;*/
	cur_group_id = -1;
	mgr_type = RUN_MGR_TYPE::NOTDEFINED;
}

/**
 * @brief Initialize.
 *
 * @param model_pars Description.
 * @param obs Description.
 * @param _filename Description.
 */
void RunManagerAbstract::initialize(const Parameters &model_pars, const Observations &obs, const string &_filename)
{
	file_stor.reset(model_pars.get_keys(), obs.get_keys(), _filename);
	failed_file_stor.reset(model_pars.get_keys(), obs.get_keys());
	if (save_all_runs)
		all_file_stor.reset(model_pars.get_keys(), obs.get_keys());
}

/**
 * @brief Initialize.
 *
 * @param par_names Description.
 * @param obs_names Description.
 * @param _filename Description.
 */
void RunManagerAbstract::initialize(const std::vector<std::string> &par_names, std::vector<std::string> &obs_names, const string &_filename)
{
	file_stor.reset(par_names, obs_names, _filename);
	failed_file_stor.reset(par_names, obs_names);
	if (save_all_runs)
		all_file_stor.reset(par_names, obs_names);
}

/**
 * @brief Reinitialize.
 *
 * @param _filename Description.
 */
void RunManagerAbstract::reinitialize(const string &_filename)
{
	vector<string> par_names = get_par_name_vec();
	vector<string> obs_names = get_obs_name_vec();
	file_stor.reset(par_names, obs_names, _filename);
	// Note: failed_file_stor and all_file_stor are intentionally NOT reset here
	// so that failed runs and the full run history accumulate across all
	// iterations/batches.
}

/**
 * @brief Initialize restart.
 *
 * @param _filename Description.
 */
void RunManagerAbstract::initialize_restart(const std::string &_filename)
{
	file_stor.init_restart(_filename);
	string failed_filename = make_failed_filename(_filename);
	if (pest_utils::check_exist_in(failed_filename))
		failed_file_stor.init_restart(failed_filename);
	else
		failed_file_stor.reset(file_stor.get_par_name_vec(), file_stor.get_obs_name_vec());
	if (save_all_runs)
	{
		string all_filename = make_all_filename(_filename);
		if (pest_utils::check_exist_in(all_filename))
			all_file_stor.init_restart(all_filename);
		else
			all_file_stor.reset(file_stor.get_par_name_vec(), file_stor.get_obs_name_vec());
	}
}

/**
 * @brief Add run.
 *
 * @param model_pars Description.
 * @param info_txt Description.
 * @param info_value Description.
 *
 * @return Description.
 */
int RunManagerAbstract::add_run(const vector<double> &model_pars, const string &info_txt, double info_value)
{
	int run_id = file_stor.add_run(model_pars, info_txt, info_value);
	return run_id;
}

/**
 * @brief Add run.
 *
 * @param model_pars Description.
 * @param info_txt Description.
 * @param info_value Description.
 *
 * @return Description.
 */
int RunManagerAbstract::add_run(const Parameters &model_pars, const string &info_txt, double info_value)
{
	int run_id = file_stor.add_run(model_pars, info_txt, info_value);
	return run_id;
}

/**
 * @brief Add run.
 *
 * @param model_pars Description.
 * @param info_txt Description.
 * @param info_value Description.
 *
 * @return Description.
 */
int RunManagerAbstract::add_run(const Eigen::VectorXd &model_pars, const string &info_txt, double info_value)
{
	int run_id = file_stor.add_run(model_pars, info_txt, info_value);
	return run_id;
}

/**
 * @brief Update run.
 *
 * @param run_id Description.
 * @param pars Description.
 * @param obs Description.
 */
void RunManagerAbstract::update_run(int run_id, const Parameters &pars, const Observations &obs)
{

	file_stor.update_run(run_id, pars, obs);
}

 const vector<string>& RunManagerAbstract::get_par_name_vec() const
 {
	return file_stor.get_par_name_vec();
 }

 const vector<string>& RunManagerAbstract::get_obs_name_vec() const
{
	return file_stor.get_obs_name_vec();
}

/**
 * @brief Run finished.
 *
 * @param run_id Description.
 *
 * @return Description.
 */
 bool RunManagerAbstract::run_finished(int run_id)
 {
	 int run_status;
	 string info_txt;
	 double info_value;
	 get_info(run_id, run_status, info_txt, info_value);
	 bool run_finished = (run_status > 0) ? true : false;
	 return run_finished;
 }

/**
 * @brief Get info.
 *
 * @param run_id Description.
 * @param run_status Description.
 * @param info_txt Description.
 * @param info_value Description.
 */
 void RunManagerAbstract::get_info(int run_id, int &run_status, std::string &info_txt, double &info_value)
 {
	  file_stor.get_info(run_id, run_status, info_txt, info_value);
 }

/**
 * @brief Get run.
 *
 * @param run_id Description.
 * @param pars Description.
 * @param obs Description.
 * @param info_txt Description.
 * @param info_value Description.
 * @param clear_old Description.
 *
 * @return Description.
 */
bool RunManagerAbstract::get_run(int run_id, Parameters &pars, Observations &obs, string &info_txt, double &info_value, bool clear_old)
{
	bool success = false;
	int status = file_stor.get_run(run_id, pars, obs, info_txt, info_value, clear_old);
	if (status > 0) success = true;
	return success;
}

/**
 * @brief Get run.
 *
 * @param run_id Description.
 * @param pars Description.
 * @param obs Description.
 * @param clear_old Description.
 *
 * @return Description.
 */
bool RunManagerAbstract::get_run(int run_id, Parameters &pars, Observations &obs,  bool clear_old)
{
	string info_txt;
	double info_value;

	return get_run(run_id, pars, obs, info_txt, info_value, clear_old);
}

/**
 * @brief Get run.
 *
 * @param run_id Description.
 * @param pars_vec Description.
 * @param obs_vec Description.
 * @param info_txt Description.
 * @param info_value Description.
 *
 * @return Description.
 */
bool RunManagerAbstract::get_run(int run_id, vector<double> &pars_vec, vector<double> &obs_vec, string &info_txt, double &info_value)
{
	bool success = false;
	int status = file_stor.get_run(run_id, pars_vec, obs_vec, info_txt, info_value);
	if (status > 0) success = true;
	return success;
}

/**
 * @brief Get run.
 *
 * @param run_id Description.
 * @param pars_vec Description.
 * @param obs_vec Description.
 *
 * @return Description.
 */
bool RunManagerAbstract::get_run(int run_id, vector<double> &pars_vec, vector<double> &obs_vec)
{
	string info_txt;
	double info_value;

	return get_run(run_id, pars_vec, obs_vec, info_txt, info_value);
}

/**
 * @brief Get run.
 *
 * @param run_id Description.
 * @param pars Description.
 * @param npars Description.
 * @param obs Description.
 * @param nobs Description.
 * @param info_txt Description.
 * @param info_value Description.
 *
 * @return Description.
 */
bool  RunManagerAbstract::get_run(int run_id, double *pars, size_t npars, double *obs, size_t nobs, string &info_txt, double &info_value)
{
	bool success = false;
	int status = file_stor.get_run(run_id, pars, npars, obs, nobs, info_txt, info_value);
	if (status > 0) success = true;
	return success;
}

/**
 * @brief Get run.
 *
 * @param run_id Description.
 * @param pars Description.
 * @param npars Description.
 * @param obs Description.
 * @param nobs Description.
 *
 * @return Description.
 */
bool  RunManagerAbstract::get_run(int run_id, double *pars, size_t npars, double *obs, size_t nobs)
{
	string info_txt;
	double info_value;

	return get_run(run_id, pars, npars, obs, nobs, info_txt, info_value);
}


/**
 * @brief Free memory.
 */
void  RunManagerAbstract::free_memory()
{
}

/**
 * @brief N run failures exceeded.
 *
 * @param id Description.
 *
 * @return Description.
 */
bool RunManagerAbstract::n_run_failures_exceeded(int id)
{
	bool ret_val;
	int istatus = file_stor.get_run_status(id);
	if (istatus<= -max_n_failure  &&  istatus > -100)
	 {
		 ret_val = true;
	 }
	else
	{
		ret_val = false;
	}
	return ret_val;
}

const std::map<std::string,std::vector<int>> RunManagerAbstract::get_run_info_map()
{
    map<string,std::vector<int>> run_info_map;
    int n_runs = file_stor.get_nruns();
    string info;
    double info_val;
    int status;
    for (int id=0; id<n_runs; ++id)
    {
        file_stor.get_info(id,status,info,info_val);
        if (run_info_map.find(info)==run_info_map.end())
        {
            run_info_map[info] = vector<int>();
        }
        run_info_map.at(info).push_back(id);
    }
    return run_info_map;
}

const std::set<int> RunManagerAbstract::get_failed_run_ids()
{
	std::set<int> failed_runs;
	int n_runs = file_stor.get_nruns();
	for (int id=0; id<n_runs; ++id)
	{
		if(n_run_failures_exceeded(id))
		 {
			 failed_runs.insert(id);
		 }
	}
	return failed_runs;
}

/**
 * @brief Get num good runs.
 *
 * @return Description.
 */
int RunManagerAbstract::get_num_good_runs(void)
{
	int n_runs_ok = file_stor.get_num_good_runs();
	return n_runs_ok;
}


/**
 * @brief Get num failed runs.
 *
 * @return Description.
 */
int RunManagerAbstract::get_num_failed_runs(void)
{
	int n_failed = 0;
	int n_runs = file_stor.get_nruns();
	for (int id=0; id<n_runs; ++id)
	 {
		if(n_run_failures_exceeded(id))
		{
			++n_failed;
		 }
	 }
	 return n_failed;
}

/**
 * @brief Get model parameters.
 *
 * @param run_id Description.
 * @param pars Description.
 *
 * @return Description.
 */
bool RunManagerAbstract::get_model_parameters(int run_id, Parameters &pars)
 {
	bool success = false;
	int status = file_stor.get_parameters(run_id, pars);
	if (status > 0) success = true;
        return success;
 }

/**
 * @brief Get observations vec.
 *
 * @param run_id Description.
 * @param data_vec Description.
 *
 * @return Description.
 */
bool RunManagerAbstract::get_observations_vec(int run_id, vector<double> &data_vec)
{
	bool success = false;
	int status = file_stor.get_observations_vec(run_id, data_vec);
	if (status > 0) success = true;
	return success;
}

/**
 * @brief Get obs template.
 *
 * @param value Description.
 *
 * @return Description.
 */
 Observations RunManagerAbstract::get_obs_template(double value) const
 {
	Observations ret_obs;
	const vector<string> &obs_name_vec = file_stor.get_obs_name_vec();
	int nobs = obs_name_vec.size();
	for(int i=0; i<nobs; ++i)
	{
		ret_obs[obs_name_vec[i]] = value;
	}
	return ret_obs;
 }
/**
 * @brief Get cur groupid.
 *
 * @return Description.
 */
 int RunManagerAbstract::get_cur_groupid()
 {
	 return cur_group_id;
 }

/**
 * @brief Run requried.
 *
 * @param run_id Description.
 *
 * @return Description.
 */
 bool RunManagerAbstract::run_requried(int run_id)
 {
	 bool ret_val;
	 int istatus = file_stor.get_run_status(run_id);
	 if (istatus <=0 && istatus > -max_n_failure)
	 {
		 ret_val = true;
	 }
	 else
	 {
		 ret_val = false;
	 }
	 return ret_val;
 }

/**
 * @brief Get outstanding run ids.
 *
 * @return Description.
 */
 vector<int> RunManagerAbstract::get_outstanding_run_ids()
 {
	 vector<int> run_ids;
	 int n_runs = file_stor.get_nruns();
	 for (int id=0; id<n_runs; ++id)
	 {
		 if(run_requried(id))
		 {
			 run_ids.push_back(id);
		 }
	 }
	 return run_ids;
 }

/**
 * @brief Mark a run as failed in the primary storage and, if the failure
 *        count has reached max_run_fail, copy its parameter values into
 *        the failed-run storage file (*case.rnf*).
 *
 * The failure count is stored in the run-status byte of the *.rnf* record
 * as the negative of the attempt count (e.g. -3 means the run was attempted
 * three times).  The parameter values and metadata (info_txt, info_value) are
 * copied from the primary storage so that users can later inspect which
 * parameter combinations caused persistent model failure.
 *
 * @param run_id  The run manager ID of the failed run.
 */
 void  RunManagerAbstract::update_run_failed(int run_id)
 {
	 file_stor.update_run_failed(run_id);

	 // If this run just crossed the failure threshold, record it in the failed runs file
	 if (n_run_failures_exceeded(run_id))
	 {
		 vector<double> pars_vec, obs_vec;
		 string info_txt;
		 double info_value;
		 file_stor.get_run(run_id, pars_vec, obs_vec, info_txt, info_value);
		 int failed_id = failed_file_stor.add_run(pars_vec, info_txt, info_value);
		 int nfail = -file_stor.get_run_status(run_id);
		 failed_file_stor.set_run_nfailed(failed_id, nfail);
		 // also record the permanently-failed run in the all-runs file
		 record_all_run(run_id);
	 }
 }

/**
 * @brief Copy a completed run from the primary storage into the all-runs
 *        storage file (*case.allruns.rns*) when SAVE_ALL_RUNS is enabled.
 *
 * The run's parameter values, metadata (info_txt, info_value) and completion
 * status are read back from the primary storage and appended to the all-runs
 * storage.  Successful runs (status > 0) additionally carry their observation
 * values; failed runs (status < 0) preserve the failure count in the run-status
 * byte.  When SAVE_ALL_RUNS is disabled this is a no-op.
 *
 * @param run_id  The run manager ID of the completed run.
 */
 void RunManagerAbstract::record_all_run(int run_id)
 {
	 if (!save_all_runs)
		 return;
	 vector<double> pars_vec, obs_vec;
	 string info_txt;
	 double info_value;
	 file_stor.get_run(run_id, pars_vec, obs_vec, info_txt, info_value);
	 int all_id = all_file_stor.add_run(pars_vec, info_txt, info_value);
	 int status = file_stor.get_run_status(run_id);
	 if (status > 0)
	 {
		 Observations obs;
		 obs.insert(all_file_stor.get_obs_name_vec(), obs_vec);
		 all_file_stor.update_run(all_id, obs);
	 }
	 else if (status < 0)
	 {
		 all_file_stor.set_run_nfailed(all_id, -status);
	 }
	 // status == 0 (never completed) is left in the added (status 0) state
 }

 const RunStorage& RunManagerAbstract::get_runstorage_ref() const
 {
	 return file_stor;
 }

/**
 * @brief Return a const reference to the failed-run storage (*case.rnf*).
 *
 * The returned RunStorage contains one entry for each run that exceeded
 * the max_run_fail threshold.  Each entry holds the parameter values that
 * were supplied to the model, associated metadata, and the number of
 * failed attempts encoded in the run-status byte.
 */
 const RunStorage& RunManagerAbstract::get_failed_runstorage_ref() const
 {
	 return failed_file_stor;
 }

/**
 * @brief Return a const reference to the all-runs storage (*case.allruns.rns*).
 *
 * When SAVE_ALL_RUNS is enabled the returned RunStorage contains one entry for
 * every completed model run (both successes and permanent failures) across all
 * iterations/batches.
 */
 const RunStorage& RunManagerAbstract::get_all_runstorage_ref() const
 {
	 return all_file_stor;
 }

/**
 * @brief Run until.
 *
 * @param condition Description.
 * @param n_nops Description.
 * @param sec Description.
 *
 * @return Description.
 */
 RunManagerAbstract::RUN_UNTIL_COND RunManagerAbstract::run_until(RUN_UNTIL_COND condition, int n_nops, double sec)
 {
	 run();
	 return RUN_UNTIL_COND::NORMAL;
 }

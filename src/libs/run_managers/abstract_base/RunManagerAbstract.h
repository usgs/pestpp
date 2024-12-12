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



class ModelExecInfo;
class Parameters;
class Observations;


class RunManagerAbstract
{
public:
	enum class RUN_UNTIL_COND { NORMAL, NO_OPS, TIME, NO_OPS_OR_TIME };
	enum class RUN_MGR_TYPE {NOTDEFINED, PANTHER, SERIAL};
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
	virtual const std::set<int> get_failed_run_ids();
    virtual const std::map<std::string,std::vector<int>> get_run_info_map();
	virtual bool get_model_parameters(int run_num, Parameters &pars);
	virtual bool get_observations_vec(int run_id, std::vector<double> &data_vec);
	virtual Observations get_obs_template(double value = -9999.0) const;
	virtual int get_total_runs(void) const {return total_runs;}
	virtual int get_num_good_runs(void);
	virtual int get_num_failed_runs(void);
	virtual bool n_run_failures_exceeded(int id);
	virtual int get_nruns(void) {return file_stor.get_nruns();}
	virtual int get_cur_groupid(void);
	virtual std::vector<int> get_outstanding_run_ids();
	virtual ~RunManagerAbstract(void) {}
	virtual std::string get_run_filename() { return file_stor.get_filename(); }
	virtual const RunStorage& get_runstorage_ref() const;
	virtual void print_run_summary(std::ostream &fout) { file_stor.print_run_summary(fout); }
	//virtual Observations get_init_run_obs() { return init_run_obs; }
	virtual std::vector<double> get_init_sim() { return init_sim;  }
	virtual void set_init_sim(std::vector<double> _init_sim) { init_sim = _init_sim; }
	virtual RUN_MGR_TYPE get_mgr_type() { return mgr_type; }

protected:
	int total_runs;
	int max_n_failure; // maximum number of times to retry a failed model run
	int cur_group_id;  // used in some of the derived classes (ie PANTHER)
	RunStorage file_stor;
	RUN_MGR_TYPE mgr_type;
	std::vector<std::string> comline_vec;
	std::vector<std::string> tplfile_vec;
	std::vector<std::string> inpfile_vec;
	std::vector<std::string> insfile_vec;
	std::vector<std::string> outfile_vec;
	bool run_requried(int run_id);
	//Observations init_run_obs;
	std::vector<double> init_sim;
	virtual void update_run_failed(int run_id);
};

#endif /*  RUNMANAGERABSTRACT_H */

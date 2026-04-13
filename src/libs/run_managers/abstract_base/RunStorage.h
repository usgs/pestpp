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

#ifndef RUN_STORAGE_H_
#define RUN_STORAGE_H_

#include <string>
#include <fstream>
#include <ostream>
#include <vector>
#include <cstdint>
#include <Eigen/Dense>
#include "network_package.h"
#include "Pest.h"

class Parameters;
class Observations;

class RunStorage {
	// This class stores a sequence of model runs in a single binary file using the following format:
	//     nruns (number of model runs stored in file)                       int_64_t
	//     run_size (number of bytes required to store each model run)       int_64_t
	//     par_name_vec_size (number of bytes required to store parameter names)  int_64_t
	//     obes_name_vec_size (number of bytes required to store observation names)  int_64_t
	//     parameter names (serialized parameter names)                               char*par_name_vec_size
	//     observation names (serialized observation names)                               char*obes_name_vec_size
	//   The following structure is repeated nruns times (ie once for each model run)
	//        run_status (flag indicating status of mode run)                         int_8_t
	//              run_status=0   this run has not yet been completed
	//			    run_status=-100   run was canceled
	//				run_status<0 and >-100  run completed and failed.  This is the number of times it failed
	//				run_status=1   run and been successfully completed
	//       info_txt  (description of model run)                                     char*
	//       info_value (variable used to store an important value.  The variable     double
	//                   depends on the type of model run being stored  )
	//       parameter_values  (parameters values for model runs)                     double*number of parameters
	//       observationn_values( observations results produced by the model run)     double*number of observations

public:
	static const double no_data;
	RunStorage(const std::string &_filename);
	void reset(const std::vector<std::string> &par_names, const std::vector<std::string> &obs_names, const std::string &_filename = std::string(""));
	void init_restart(const std::string &_filename);
	virtual int add_run(const std::vector<double> &model_pars, const std::string &info_txt="", double info_value=no_data);
	virtual int add_run(const Parameters &pars, const std::string &info_txt="", double info_value=no_data);
	virtual int add_run(const Eigen::VectorXd &model_pars, const std::string &info_txt="", double info_value=no_data);
	void copy(const RunStorage &rhs_rs);
	void update_run(int run_id, const Parameters &pars, const Observations &obs);
	void update_run(int run_id, const Observations &obs);
	void update_run(int run_id, const std::vector<char> serial_data);
	void update_run_failed(int run_id);
	void set_run_nfailed(int run_id, int nfail);
	int get_nruns();
	int get_num_good_runs();
	int increment_nruns();
	const std::vector<std::string>& get_par_name_vec()const;
	const std::vector<std::string>& get_obs_name_vec()const;
	int get_run_status(int run_id);
	void get_info(int run_id, int &run_status, std::string &info_txt, double &info_value);
	int get_run(int run_id, Parameters &pars, Observations &obs, bool clear_old=true);
	int get_run(int run_id, Parameters &pars, Observations &obs, std::string &info_txt, double &info_value, bool clear_old=true);
	int get_run(int run_id, double *pars, size_t npars, double *obs, size_t nobs);
	int get_run(int run_id, double *pars, size_t npars, double *obs, size_t nobs, std::string &info_txt, double &info_value);
	int get_run(int run_id, std::vector<double> &pars_vec, std::vector<double> &obs_vec,
		    std::string &info_txt, double &info_value);
	int get_run(int run_id, std::vector<double> &pars_vec, std::vector<double> &obs_vec);
	int get_parameters(int run_id, Parameters &pars);
	std::vector<char> get_serial_pars(int run_id);
	int get_observations_vec(int run_id, std::vector<double> &data_vec);
	int get_observations(int run_id, Observations &obs);
	static void export_diff_to_text_file(const std::string &in1_filename, const std::string &in2_filename, const std::string &out_filename);
	void free_memory();
	std::string get_filename() { return filename; }
	void print_run_summary(std::ostream &fout);

	~RunStorage();
private:
	static const int info_txt_length = NetPackage::DESC_LEN;
	std::string filename;
	mutable std::fstream buf_stream;
	std::streamoff beg_run0;
	std::streamoff run_byte_size;
	std::streamoff run_par_byte_size;
	std::streamoff run_data_byte_size;
	std::vector<std::string> par_names;
	std::vector<std::string> obs_names;
	void check_rec_size(const std::vector<char> &serial_data) const;
	void check_rec_id(int run_id);
	std::int8_t get_run_status_native(int run_id);
	std::streamoff get_stream_pos(int run_id);
};

#endif //RUN_STORAGE_H_

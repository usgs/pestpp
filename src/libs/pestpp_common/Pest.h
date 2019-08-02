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

#ifndef PEST_H_
#define PEST_H_

#include <string>
#include <map>
#include <set>
#include <vector>
#include "ParamTransformSeq.h"
#include "pest_data_structs.h"
#include "PriorInformation.h"
#include "Regularization.h"

using namespace std;
class FileManager;
class PriorInformation;

class Pest {
public:
	enum class LimitType { NONE, LBND, UBND, REL, FACT };
	friend ostream& operator<< (ostream &os, const Pest& val);
	Pest();
	void set_defaults();
	void check_inputs(ostream &f_rec);
	void check_io();
	int process_ctl_file(ifstream &fin, string pst_filename, ofstream &f_rec);
	int process_ctl_file(ifstream &fin, string pst_filename);
	int get_n_adj_par(){ return n_adj_par; }
	const Parameters& get_ctl_parameters() const {return ctl_parameters;}
	const Observations& get_ctl_observations() const {return observation_values;}
	const ParameterInfo& get_ctl_parameter_info()const {return ctl_parameter_info;}
	const ParameterGroupInfo& get_base_group_info() const {return  base_group_info;}
	ParameterGroupInfo* get_base_group_info_ptr() { return  &base_group_info; }
	const ObservationInfo &get_ctl_observation_info() const {return observation_info;}
	const std::map<std::string, std::string> get_observation_groups() const;
	const PriorInformation &get_prior_info() {return prior_info;}
	PriorInformation *get_prior_info_ptr() {return &prior_info;}
	const SVDInfo& get_svd_info() const {return svd_info;}
	const ControlInfo&  get_control_info() const {return control_info;}
	const ParamTransformSeq& get_base_par_tran_seq() const {return base_par_transform;}
	const vector<string> &get_ctl_ordered_par_names() const {return ctl_ordered_par_names;}
	const vector<string> &get_ctl_ordered_obs_names() const {return ctl_ordered_obs_names;}
	const vector<string> &get_ctl_ordered_par_group_names() { return ctl_ordered_par_group_names; }
	const vector<string> &get_ctl_ordered_obs_group_names() { return ctl_ordered_obs_group_names; }
	vector<string> *get_ctl_ordered_obs_group_names_ptr() { return &ctl_ordered_obs_group_names; }
	const vector<string> &get_ctl_ordered_pi_names() const { return ctl_ordered_pi_names; }
	vector<string> *get_ctl_ordered_pi_names_ptr() { return &ctl_ordered_pi_names; }
	const vector<string> get_ctl_ordered_nz_obs_names();
	const vector<string> get_ctl_ordered_adj_par_names();
	const ModelExecInfo &get_model_exec_info() {return model_exec_info;}
	const  vector<string> &get_comline_vec();
	const  vector<string> &get_tplfile_vec();
	const  vector<string> &get_inpfile_vec();
	const  vector<string> &get_insfile_vec();
	const  vector<string> &get_outfile_vec();
	const PestppOptions &get_pestpp_options() const {return pestpp_options;}
    PestppOptions *get_pestpp_options_ptr() { return &pestpp_options; }
	ObservationInfo *get_observation_info_ptr() { return &observation_info; }
	DynamicRegularization* get_regul_scheme_ptr() { return regul_scheme_ptr; }
	const ParetoInfo &get_pareto_info() const { return pareto_info; }
	vector<string> get_nonregul_obs() const;
	string get_pst_filename() { return pst_filename; }
	void enforce_par_limits(Parameters &update_active_ctl_pars, const Parameters &last_active_ctl_pars, bool enforce_chglim=true, bool enforce_bounds=false);
	map<string,double> get_pars_at_near_bounds(const Parameters &pars, double tol=0.0);
	pair<Parameters,Parameters> get_effective_ctl_lower_upper_bnd(Parameters &pars);
	
	virtual ~Pest();
	
private:
	int n_adj_par = 0;
	ControlInfo control_info;
	ParetoInfo pareto_info;
	SVDInfo svd_info;
	Parameters ctl_parameters;
	ParameterInfo ctl_parameter_info;
	ParameterGroupInfo base_group_info;
	Observations observation_values;
	ObservationInfo observation_info;
	PriorInformation prior_info;
	ModelExecInfo model_exec_info;
	PestppOptions pestpp_options;
	ParamTransformSeq base_par_transform;
	vector<string> ctl_ordered_par_names;
	vector<string> ctl_ordered_obs_names;
	vector<string> ctl_ordered_par_group_names;
	vector<string> ctl_ordered_obs_group_names;
	vector<string> ctl_ordered_pi_names;
	DynamicRegularization *regul_scheme_ptr;
	map<int,string> other_lines;
	string pst_filename;
};
ostream& operator<< (ostream &os, const Pest& val);
#endif /* PEST_H_ */

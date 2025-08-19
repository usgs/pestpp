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
#include "utilities.h"

using namespace std;
class FileManager;
class PriorInformation;


class Pest {
public:
	enum class LimitType { NONE, LBND, UBND, REL, FACT };
	friend ostream& operator<< (ostream &os, const Pest& val);
	Pest();
	void set_defaults();
	void set_default_dynreg();
	void check_inputs(ostream &f_rec, bool forgive=false, bool forgive_parchglim=true,int cycle=NetPackage::NULL_DA_CYCLE);
	void check_io(ofstream& f_rec, bool echo_errors=true);
	//int process_ctl_file_old(ifstream &fin, string pst_filename, ofstream &f_rec);
	//int process_ctl_file_old(ifstream &fin, string pst_filename);
	int process_ctl_file(ifstream& fin, string pst_filename, ofstream& f_rec);
	int process_ctl_file(ifstream& fin, string pst_filename);

	int get_n_adj_par(){ return n_adj_par; }
	const Parameters& get_ctl_parameters() const {return ctl_parameters;}
	Parameters& get_ctl_parameters_4_mod() { return ctl_parameters; }
	const Observations& get_ctl_observations() const {return observation_values;}
	Observations& get_ctl_observations_4_mod() { return observation_values;  }
	const ParameterInfo& get_ctl_parameter_info()const {return ctl_parameter_info;}
	void set_ctl_parameter_info(ParameterInfo pi) { ctl_parameter_info = pi; }
	void set_ctl_observation_info(ObservationInfo oi) { observation_info = oi; }

	ParameterInfo get_ctl_parameter_info_copy() { return ctl_parameter_info; }
	ParameterInfo* get_ctl_parameter_info_ptr_4_mod() { return &ctl_parameter_info; }
	const ParameterGroupInfo& get_base_group_info() const {return  base_group_info;}
	ParameterGroupInfo* get_base_group_info_ptr_4_mod() { return  &base_group_info; }
	ParameterGroupInfo* get_base_group_info_ptr() { return  &base_group_info; }
	const ObservationInfo &get_ctl_observation_info() const {return observation_info;}
	ObservationInfo get_ctl_observation_info_copy() { return observation_info; }
	ObservationInfo& get_ctl_observation_info_4_mod() { return observation_info; }
	const ObservationInfo* get_ctl_observation_info_ptr() const { return &observation_info; }
	void set_observation_info(ObservationInfo oi) { observation_info = oi; }

	const std::map<std::string, std::string> get_observation_groups() const;
	const PriorInformation &get_prior_info() {return prior_info;}
	PriorInformation *get_prior_info_ptr() {return &prior_info;}
	const SVDInfo& get_svd_info() const {return svd_info;}
	const ControlInfo&  get_control_info() const {return control_info;}
	ControlInfo& get_control_info_4_mod() {return control_info;}
	const ParamTransformSeq& get_base_par_tran_seq() const {return base_par_transform;}
	ParamTransformSeq& get_base_par_tran_seq_4_mod() { return base_par_transform; }
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
	pair<string,double> enforce_par_limits(PerformanceLog* perfomance_log, Parameters &update_active_ctl_pars, const Parameters &last_active_ctl_pars, bool enforce_chglim=true, bool enforce_bounds=false);
	map<string,double> get_pars_at_near_bounds(const Parameters &pars, double tol=0.0);
	pair<Parameters,Parameters> get_effective_ctl_lower_upper_bnd(Parameters &pars);
	map<string, double> calc_par_dss(const Jacobian& jac, ParamTransformSeq& par_transform);

	Pest get_child_pest(int icycle);
	void child_pest_update(int icycle, bool keep_order=true);
	//vector<int> get_assim_cycles(ofstream& f_rec, vector<int> unique_cycles);
	vector<int> get_assim_dci_cycles(ofstream& f_rec, vector<int> unique_cycles);
	//vector<int> get_assim_cycles(ofstream& f_rec);
	void assign_da_cycles(ofstream& f_rec);
	vector<pair<string, int>> extract_cycle_numbers2(ofstream& f_rec, string section_name, vector<string> possible_name_cols);
    vector<pair<string, DaCycleInfo>> extract_cycle_info(ofstream& f_rec, string section_name, vector<string> possible_name_cols);
	map<string, double> get_ext_file_double_map(const string& section_name, const string& col_name);
	map<string, string> get_ext_file_string_map(const string& section_name, const string& col_name);
    int get_num_ext_file_maps() const {return efiles_map.size();}
    void clear_ext_files() { efiles_map.clear(); }
	map<string, vector<pest_utils::ExternalCtlFile>>& get_efiles_map() { return efiles_map; }
	void release_unused_for_agent();
	virtual ~Pest();
	
protected:
	//this is the list of external file cols that have meaning...
	set<string> efile_keep_cols{ "standard_deviation", "obsnme","parnme","name", "upper_bound","lower_bound", "cycle", "state_par_link","drop_violations","greater_than","less_than","link_to" };
	int n_adj_par = 0;
	string prior_info_string;
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
	vector<int> assimilation_cycles;
	DynamicRegularization *regul_scheme_ptr;
	map<int,string> other_lines;
	string pst_filename;

	set<string> s_pargp, s_obgnme; 

	pair<string, string> parse_keyword_line(ofstream &f_rec, const string &line);
	void throw_control_file_error(ofstream& f_rec, const string &message, bool should_throw=true, bool echo=true);
	void check_report_assignment(ofstream& f_rec, PestppOptions::ARG_STATUS stat, const string &key, const string &org_value);

	void tokens_to_par_group_rec(ofstream &f_rec, const vector<string>& tokens);
	void tokens_to_par_rec(ofstream &f_rec, const vector<string>& tokens,TranFixed *t_fixed, TranLog10 *t_log, TranScale *t_scale, TranOffset *t_offset);
	void tokens_to_obs_group_rec(ofstream& f_rec, const vector<string>& tokens);
	void tokens_to_obs_rec(ofstream& f_rec, const vector<string> &tokens);
	void tokens_to_pi_rec(ofstream& f_rec, const vector<string>& tokens);
	void tokens_to_pi_rec(ofstream& f_rec, const string& line_upper);
	void rectify_par_groups();
	map<string, vector<pest_utils::ExternalCtlFile>> efiles_map;
	DaCycleInfo parse_cycle_str(string& raw_cycle_val, string& efilename, int row, ofstream& f_rec);
};
ostream& operator<< (ostream &os, const Pest& val);
bool every_cycle(const DaCycleInfo& dci);
bool cycle_in_range(int cycle,const DaCycleInfo& dci);

#endif /* PEST_H_ */

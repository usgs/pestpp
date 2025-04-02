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
#ifndef PEST_DATAS_STRUCTS_H_
#define PEST_DATAS_STRUCTS_H_

#include <string>
#include <map>
#include <unordered_map>
#include <set>
#include <vector>
#include <random>
#include "Transformable.h"


//class TransformSeq<;
class LaTridiagMatDouble;



struct DaCycleInfo{
    int start=0;
    int stop = -999;
    int stride = 1;
};

class ParameterGroupRec {
public:
	string name;
	string inctyp;
	double derinc;
	double derinclb;
	string forcen;
	double derincmul;
	string dermthd;
	double splitthresh;
	double splitreldiff;
	ParameterGroupRec(const string &_name="", const string &_inctyp="", double _derinc=0.0, double _derinclb=0.0,
		const string &_forcen="", double _derincmul=0.0,
		const string &_dermthd = "", double _splitthresh = 0.0, double _splitreldiff=.50)
		: name(_name), inctyp(_inctyp), derinc(_derinc), derinclb(_derinclb), forcen(_forcen),
		derincmul(_derincmul), dermthd(_dermthd), splitthresh(_splitthresh),
		splitreldiff(_splitreldiff){}
	ParameterGroupRec(const ParameterGroupRec &rhs) {*this=rhs;}
	ParameterGroupRec& operator=(const ParameterGroupRec &rhs);
	void set_defaults();
	void set_name(string& _name) { name = _name; }
};

ostream& operator<< (ostream &os, const ParameterGroupRec& val);
ostream& operator<< (ostream &os, const map<string, ParameterGroupRec> &val);

class ParameterGroupInfo {
	friend ostream& operator<< (ostream &os, const ParameterGroupInfo &val);

public:
	ParameterGroupInfo() {}
	ParameterGroupInfo(const ParameterGroupInfo&rhs) {*this=rhs;}
	void insert_group(const string &group_name, ParameterGroupRec &rec);
    void clear() {groups.clear();parameter2group.clear();}
	/** @brief Creates a link from a parameter to its parameters group.

	This method add a record to a hash table to link the specified parameter
	to the specified group
	*/
	void insert_parameter_link(const string &parameter_name, const string & group_name);
	const ParameterGroupRec* get_group_rec_ptr(const string &par_name) const;
	const ParameterGroupRec* get_group_by_groupname(const string &group_name) const { return groups.at(group_name); }
	ParameterGroupRec* get_group_by_groupname_4_mod(const string &group_name) { return groups.at(group_name); }
	ParameterGroupRec* get_group_rec_ptr_4_mod(const string &par_name);
	string get_group_name(const string &par_name) const;
	const ParameterGroupInfo& operator=(const ParameterGroupInfo &rhs);
	bool have_switch_derivative() const;
	vector<string> get_group_names() const;
	void par_erase(const string& par_name) { parameter2group.erase(par_name); }
	void grp_erase(const string& grp_name) { groups.erase(grp_name);
	}
	~ParameterGroupInfo();
private:
	unordered_map<string, ParameterGroupRec*> groups;
	unordered_map<string, ParameterGroupRec*> parameter2group;

};


class ParameterRec {
public:
	enum class TRAN_TYPE {NONE, FIXED, TIED, LOG};
	string chglim;
	double lbnd;
	double ubnd;
	double init_value;
	double scale;
	double offset;
	string group;
	int dercom;
	//int cycle;
	DaCycleInfo dci;
	TRAN_TYPE tranform_type;
	ParameterRec() : chglim(""), lbnd(0.0), ubnd(0.0), init_value(0.0), group(""),
		dercom(1), tranform_type(TRAN_TYPE::NONE), scale(1.0), offset(0.0),
		dci(DaCycleInfo()){}
	bool is_active() const { return !(tranform_type == TRAN_TYPE::FIXED || tranform_type == TRAN_TYPE::TIED); }
};
ostream& operator<< (ostream &os, const ParameterRec& val);

class ParameterInfo {
	friend ostream& operator<< (ostream &os, const ParameterInfo& val);
public:
	Parameters get_low_bnd(const vector<string> &keys) const;
	Parameters get_up_bnd(const vector<string> &keys) const;
	Parameters get_init_value(const vector<string> &keys) const;
	const ParameterRec* get_parameter_rec_ptr(const string &name) const;
	ParameterRec* get_parameter_rec_ptr_4_mod(const string &name);
	void insert(const string &name, const ParameterRec &rec) {parameter_info[name] = rec;}
	void erase(const string& name) { parameter_info.erase(name);}
	ParameterInfo() {}
	~ParameterInfo() {}
private:
	unordered_map<string, ParameterRec> parameter_info;
};

ostream& operator<< (ostream &os, const ParameterInfo& val);

class ObservationGroupRec {
public:
	double gtarg;  // optional
	string covfile; // optional
	ObservationGroupRec() : gtarg(0.0), covfile(""){};
	static bool is_regularization(const std::string &grp_name);

};
ostream& operator<< (ostream &os, const ObservationGroupRec& val);

class ObservationRec {
public:
	double weight;
	string group;
	//int cycle;
	DaCycleInfo dci;
	ObservationRec(double _weight=0.0,const string &_group="", int _cycle = 0)
		: weight(_weight), group(_group), dci(DaCycleInfo())  {}
	bool is_regularization() const;
};
ostream& operator<< (ostream &os, const ObservationRec& val);

class ObservationInfo {
public:
	ObservationInfo() {}
	ObservationInfo(const ObservationInfo & rhs) : groups(rhs.groups), observations(rhs.observations){}
	bool is_regularization(const string &obs_name) const;
	unordered_map<string, ObservationGroupRec> groups;
	unordered_map<string, ObservationRec> observations;
	double get_weight(const string &obs_name) const;
	void set_weight(const string &obs_name, double value);
	string get_group(const string &obs_name) const;
	const ObservationRec* get_observation_rec_ptr(const string &name) const;
	ObservationRec* get_observation_rec_ptr_4_mod(const string& name);
	const ObservationGroupRec* get_group_rec_ptr(const string &name) const;
	Observations get_regulatization_obs(const Observations &obs_in);
	int get_nnz_obs() const;
	int get_nnz_obs_and_reg() const;
	vector<string> get_groups() const;
	void reset_group_weights(string &group, double val);
	void scale_group_weights(string &group, double scale_val);
	void erase_ob(const string& name) { observations.erase(name);}
	void erase_gp(const string& name) { groups.erase(name); }

};

class ModelExecInfo {
public:
	std::vector<std::string> comline_vec;
	std::vector<std::string> tplfile_vec;
	std::vector<std::string> inpfile_vec;
	std::vector<std::string> insfile_vec;
	std::vector<std::string> outfile_vec;
	//std::vector<int> incycle_vec;
	//std::vector<int> outcycle_vec;
	std::vector<DaCycleInfo> incycle_dci_vec;
    std::vector<DaCycleInfo> outcycle_dci_vec;

};

class ParetoInfo {
public:
	double wf_start, wf_fin, wf_inc;
	string obsgroup;
	int niter_start, niter_gen, niter_fin;
	ParetoInfo() : wf_start(1.0), wf_fin(1.0), wf_inc(0.0),
		obsgroup(""), niter_start(1), niter_gen(1), niter_fin(1) {};
};


class PestppOptions {
public:
	enum SVD_PACK { EIGEN, PROPACK, REDSVD };
	enum MAT_INV { Q12J, JTQJ };
	enum GLOBAL_OPT { NONE, OPT_DE, OPT_MOEA};
	enum GLMNormalForm { IDENT,DIAG, PRIOR };
	enum ARG_STATUS {ARG_ACCEPTED, ARG_DUPLICATE, ARG_NOTFOUND, ARG_INVALID};
	PestppOptions() { use_da_args=false; }



	//void parce_line(const string &line);
	map<string,ARG_STATUS> parse_plusplus_line(const string &line);
	vector<string> notfound_args;
	ARG_STATUS assign_value_by_key(string key, const string org_value);
	void rectify_ies_da_args();

	bool assign_value_by_key_sqp(const string& key, const string& value, const string& org_value);
	bool assign_mou_value_by_key(const string& key, const string& value, const string& org_value);
	bool assign_da_value_by_key(const string& key, const string& value, const string& org_value);
	bool assign_ies_value_by_key(const string& key, const string& value, const string& org_value);
	bool assign_value_by_key_continued(const string& key, const string& value, const string& org_value);

	int get_max_n_super() const { return max_n_super; }
	double get_super_eigthres() const { return super_eigthres; }
	int get_n_iter_base() const { return n_iter_base; }
	int get_n_iter_super() const { return n_iter_super; }
	SVD_PACK get_svd_pack() const { return svd_pack; }
	double get_super_relparmax() const { return super_relparmax; }
	int get_max_run_fail() const { return max_run_fail; }


	int get_max_super_frz_iter()const { return max_super_frz_iter; }
	int get_max_reg_iter()const { return max_reg_iter; }
	const vector<double>& get_base_lambda_vec() const { return base_lambda_vec; }
	void set_base_lambda_vec(vector<double> _vals) { base_lambda_vec = _vals; }
	const vector<double>& get_lambda_scale_vec() const { return lambda_scale_vec; }
	void set_lambda_scale_vec(vector<double> sv) {lambda_scale_vec = sv; }
	bool get_iter_summary_flag() const { return iter_summary_flag; }
	bool get_der_forgive() const { return der_forgive; }
	void set_der_forgive(bool _flag) { der_forgive = _flag; }
	GLOBAL_OPT get_global_opt() const { return global_opt; }
	double get_de_f() const { return de_f; }
	void set_de_f(double _val) { de_f = _val; }
	double get_de_cr() const { return de_cr; }
	void set_de_cr(double _val) { de_cr = _val; }
	int get_de_npopulation() const { return de_npopulation; }
	void set_de_npopulation(int _val) { de_npopulation = _val; }
	int get_de_max_gen() const { return de_max_gen; }
	void set_de_max_gen(int _val ) { de_max_gen = _val; }
	bool get_de_dither_f() const { return de_dither_f; }
	void set_de_dither_f(bool _val) {de_dither_f = _val; }
	void set_global_opt(const GLOBAL_OPT _global_opt) { global_opt = _global_opt; }
	void set_max_n_super(int _max_n_super) { max_n_super = _max_n_super; }
	void set_super_eigthres(double _super_eigthres) { super_eigthres = _super_eigthres; }
	void set_n_iter_base(int _n_iter_base) { n_iter_base = _n_iter_base; }
	void set_n_iter_super(int _n_iter_super) { n_iter_super = _n_iter_super; }
	void set_svd_pack(const SVD_PACK _svd_pack) { svd_pack = _svd_pack; }
	void set_super_relparmax(double _super_relparmax) { super_relparmax = _super_relparmax; };
	void set_max_run_fail(int _max_run_fail) { max_run_fail = _max_run_fail; }
	void set_max_super_frz_iter(int n) { max_super_frz_iter = n; }
	void set_max_reg_iter(int n) { max_reg_iter = n; }
	void set_iter_summary_flag(const bool _iter_summary_flag) { iter_summary_flag = _iter_summary_flag; }
	void set_uncert_flag(bool _flag) { uncert = _flag; }
	bool get_uncert_flag()const { return uncert; }
	void set_prediction_names(vector<string> _names) { prediction_names = _names; }
	vector<string> get_prediction_names()const { return prediction_names; }
	void set_parcov_filename(string _filename) { parcov_filename = _filename; }
	string get_parcov_filename()const { return parcov_filename; }
	void set_obscov_filename(string _filename) { obscov_filename = _filename; }
	string get_obscov_filename()const { return obscov_filename; }
	void set_basejac_filename(string _filename) { basejac_filename = _filename; }
	string get_basejac_filename()const { return basejac_filename; }
	int get_glm_num_reals() const { return glm_num_reals; }
	void set_glm_num_reals(int _glm_num_reals) { glm_num_reals = _glm_num_reals; }
	GLMNormalForm get_glm_normal_form() const { return glm_normal_form;}
	void set_glm_normal_form(GLMNormalForm form) { glm_normal_form = form; }
	bool get_glm_debug_der_fail() const { return glm_debug_der_fail; }
	void set_glm_debug_der_fail(bool _flag) { glm_debug_der_fail = _flag;}
	bool get_glm_debug_lamb_fail() const { return glm_debug_lamb_fail; }
	void set_glm_debug_lamb_fail(bool _flag) { glm_debug_lamb_fail = _flag; }
	bool get_glm_debug_real_fail() const { return glm_debug_real_fail; }
	void set_glm_debug_real_fail(bool _flag) { glm_debug_real_fail = _flag; }
	bool get_glm_accept_mc_phi() const { return glm_accept_mc_phi; }
	void set_glm_accept_mc_phi(bool _flag) { glm_accept_mc_phi = _flag; }
	bool get_glm_rebase_super() const { return glm_rebase_super; }
	void set_glm_rebase_super(bool _flag) { glm_rebase_super = _flag; }
    bool get_glm_iter_mc() const { return glm_iter_mc; }
    void set_glm_iter_mc(bool _flag) { glm_iter_mc = _flag; }
    bool get_glm_debug_high_2nd_iter_phi() const {return glm_debug_high_2nd_iter_phi;}
    void set_glm_debug_high_2nd_iter_phi(bool _flag) {glm_debug_high_2nd_iter_phi = _flag;}





    string get_sweep_parameter_csv_file()const { return sweep_parameter_csv_file; }
	void set_sweep_parameter_csv_file(string _file) { sweep_parameter_csv_file = _file; }
	string get_sweep_output_csv_file()const { return sweep_output_csv_file; }
	void set_sweep_output_csv_file(string _file) { sweep_output_csv_file = _file; }
	int get_sweep_chunk()const { return sweep_chunk; }
	void set_sweep_chunk(int _chunk) { sweep_chunk = _chunk; }
	bool get_sweep_forgive()const { return sweep_forgive; }
	void set_sweep_forgive(bool _forgive) { sweep_forgive = _forgive; }
	bool get_sweep_base_run()const { return sweep_base_run; }
	void set_sweep_base_run(bool _base) { sweep_base_run = _base; }
	bool get_sweep_include_regul_phi() const {return sweep_include_regul_phi;}
	void set_sweep_include_regul_phi(bool _flag) {sweep_include_regul_phi = _flag;}

	bool get_jac_scale()const { return jac_scale; }
	void set_jac_scale(bool _jac_scale) { jac_scale = _jac_scale; }

	void set_hotstart_resfile(string _res_file) { hotstart_resfile = _res_file; }
	string get_hotstart_resfile() const { return hotstart_resfile; }


    bool get_tie_by_group() const { return tie_by_group; }
	void set_tie_by_group(bool _flag) { tie_by_group = _flag; }

	string get_opt_obj_func()const { return opt_obj_func; }
	void set_opt_obj_func(string _opt_obj_func) { opt_obj_func = _opt_obj_func; }
	bool get_opt_coin_log()const { return opt_coin_log; }
	void set_opt_coin_log(bool _log) { opt_coin_log = _log; }
	bool get_opt_skip_final()const { return opt_skip_final; }
	void set_opt_skip_final(bool _skip_final) { opt_skip_final = _skip_final; }


	vector<string> get_opt_dec_var_groups()const { return opt_dec_var_groups; }
	void set_opt_dec_var_groups(vector<string> _grps) { opt_dec_var_groups = _grps; }
	vector<string> get_opt_ext_var_groups()const { return opt_external_var_groups; }
	void set_opt_ext_var_groups(vector<string> _grps) { opt_external_var_groups = _grps; }
	vector<string> get_opt_constraint_groups()const { return opt_constraint_groups; }
	void set_opt_constraint_groups(vector<string> _grps) { opt_constraint_groups = _grps; }
	double get_opt_risk()const { return opt_risk; }
	void set_opt_risk(double _risk) { opt_risk = _risk; }
	double get_opt_direction()const { return opt_direction; }
	void set_opt_direction(double _direction) { opt_direction = _direction; }
	double get_opt_iter_tol()const { return opt_iter_tol; }
	void set_opt_iter_tol(double _tol) { opt_iter_tol = _tol; }
	int get_opt_recalc_fosm_every()const { return opt_recalc_fosm_every; }
	void set_opt_recalc_fosm_every(int _every) { opt_recalc_fosm_every = _every; }
	double get_opt_iter_derinc_fac() const { return opt_iter_derinc_fac; }
	void set_opt_iter_derinc_fac(double _opt_iter_derinc_fac) { opt_iter_derinc_fac = _opt_iter_derinc_fac; }
	bool get_opt_include_bnd_pi()const { return opt_include_bnd_pi; }
	void set_opt_include_bnd_pi(bool _include_bnd_pi) { opt_include_bnd_pi = _include_bnd_pi; }
	bool get_opt_std_weights()const { return opt_std_weights; }
	void set_opt_std_weights(bool _opt_std_weights) { opt_std_weights = _opt_std_weights; }
	int get_opt_stack_size()const { return opt_stack_size; }
	void set_opt_stack_size(int _size) { opt_stack_size = _size; }
	string get_opt_par_stack()const { return opt_par_stack; }
	void set_opt_par_stack(string _stack) { opt_par_stack = _stack; }
	string get_opt_obs_stack()const { return opt_obs_stack; }
	void set_opt_obs_stack(string _stack) { opt_obs_stack = _stack; }
	string get_opt_chance_points() const { return opt_chance_points; }
	void set_opt_chance_points(string chance_points) { opt_chance_points = chance_points; }


	string get_sqp_dv_en()const { return sqp_dv_en; }
	void set_sqp_dv_en(string _file) { sqp_dv_en = _file; }
	string get_sqp_obs_restart_en()const { return sqp_obs_restart_en; }
	void set_sqp_obs_restart_en(string _file) { sqp_obs_restart_en = _file; }
	int get_sqp_num_reals()const { return sqp_num_reals; }
	void set_sqp_num_reals(int _num_reals) { sqp_num_reals = _num_reals; }
	bool get_sqp_update_hessian()const { return sqp_update_hessian; }
	void set_sqp_update_hessian(bool _flag) { sqp_update_hessian = _flag; }
	vector<double> get_sqp_scale_facs() const { return sqp_scale_facs; }  // perhaps change arg name to sqp_alpha_mults
	void set_sqp_scale_facs(vector<double> _mults) { sqp_scale_facs = _mults; }


	string get_mou_generator() const { return mou_generator; }
	void set_mou_generator(string name) { mou_generator = name; }
	int get_mou_population_size() const { return mou_population_size; }
	void set_mou_population_size(int size) { mou_population_size = size; }
	string get_mou_dv_population_file() const { return mou_dv_population_file; }
	void set_mou_dv_population_file(string name) { mou_dv_population_file = name; }
	string get_mou_obs_population_restart_file() const { return mou_obs_population_restart_file; }
	void set_mou_obs_population_restart_file(string name) { mou_obs_population_restart_file = name; }
	vector<string> get_mou_objectives() const { return mou_objectives; }
	void set_mou_objectives(const vector<string>& objs) { mou_objectives = objs; }
	int get_mou_max_archive_size() const { return mou_max_archive_size; }
	void set_mou_max_archive_size(int size) { mou_max_archive_size = size; }
	bool get_mou_risk_obj()const { return mou_risk_obj; }
	void set_mou_risk_obj(bool _flag) { mou_risk_obj = _flag; }
	int get_mou_verbose_level() const { return mou_verbose_level; }
	void set_mou_verbose_level(int size) { mou_verbose_level = size; }
	string get_mou_env_selector() const { return mou_env_selector; }
	void set_mou_env_selector(string _env) { mou_env_selector = _env; }
	string get_mou_mating_selector() const { return mou_mating_selector; }
	void set_mou_mating_selector(string _mating) { mou_mating_selector = _mating; }
	double get_mou_crossover_probability() const { return mou_crossover_prob; }
	void set_mou_crossover_probability(double _val) { mou_crossover_prob = _val; }
	double get_mou_mutation_probability() const { return mou_mutation_prob; }
	void set_mou_mutation_probability(double _val) { mou_mutation_prob = _val; }
	double get_mou_de_f() const { return mou_de_f; }
	void set_mou_de_f(double _val) { mou_de_f = _val; }
	int get_mou_save_population_every() const { return mou_save_population_every; }
	void set_mou_save_population_every(int every) { mou_save_population_every = every; }
	double get_mou_pso_omega() const { return mou_pso_omega; }
	void set_mou_pso_omega(double val) { mou_pso_omega = val; }
	vector<double> get_mou_pso_social_const() const { return mou_pso_social_const; }
	void set_mou_pso_social_const(vector<double> _vals) { mou_pso_social_const = _vals; }
	vector<double> get_mou_pso_cognitive_const() const { return mou_pso_cognitive_const; }
	void set_mou_pso_cognitive_const(vector<double> _vals) { mou_pso_cognitive_const = _vals; }
	double get_mou_pso_alpha() const { return mou_pso_alpha; }
	void set_mou_pso_alpha(double val) { mou_pso_alpha = val; }
	double get_mou_pso_rramp() const { return mou_pso_rramp; }
	void set_mou_pso_rramp(double val) { mou_pso_rramp = val; }
	double get_mou_pso_rfit() const { return mou_pso_rfit; }
	void set_mou_pso_rfit(double val) { mou_pso_rfit = val; }
	vector<double> get_mou_pso_inertia() const { return mou_pso_inertia; }
	void set_mou_pso_inertia(vector<double> _vals) { mou_pso_inertia = _vals; }
	double get_mou_pso_vmax_factor() const { return mou_pso_vmax_factor; }
	void set_mou_pso_vmax_factor(double _val) { mou_pso_vmax_factor = _val; }
	string get_mou_pso_dv_bound_restoration() const { return mou_pso_dv_bound_restoration; }
	void set_mou_pso_dv_bound_restoration(string name) { mou_pso_dv_bound_restoration = name; }
	int get_mou_max_nn_search() const { return mou_max_nn_search;}
	void set_mou_max_nn_search(int val) { mou_max_nn_search = val; }
	string get_mou_outer_repo_obs_file() const { return mou_outer_repo_obs_file; }
	void set_mou_outer_repo_obs_file(string name) { mou_outer_repo_obs_file = name; }
	double get_mou_hypervolume_extreme() const { return mou_hypervolume_extreme; }
	void set_mou_hypervolume_extreme(double val) { mou_hypervolume_extreme = val; }
	int get_mou_infill_size() const { return mou_infill_size; }
	void set_mou_infill_size(int size) { mou_infill_size = size; }
	double get_mou_ppd_beta() const { return mou_ppd_beta; }
	void set_mou_ppd_beta(double val) { mou_ppd_beta = val; }
	double get_mou_fit_epsilon() const { return mou_fit_epsilon; }
	void set_mou_fit_epsilon(double val) { mou_fit_epsilon = val; }
	double get_mou_fit_gamma() const { return mou_fit_gamma; }
	void set_mou_fit_gamma(double val) { mou_fit_gamma = val; }
	int get_mou_resample_every()const { return mou_resample_every; }
	void set_mou_resample_every(int _every) { mou_resample_every = _every; }
	string get_mou_resample_command()const { return mou_resample_command; }
	void set_mou_resample_command(string _rescmd) { mou_resample_command = _rescmd; }
	string get_mou_population_schedule() const {return mou_population_schedule;}
    void set_mou_population_schedule(string fname) {mou_population_schedule = fname;}
	int get_mou_simplex_reflections() const { return mou_simplex_reflections; }
	void set_mou_simplex_reflections(int val) { mou_simplex_reflections = val; }
	vector<double> get_mou_simplex_factors() const { return mou_simplex_factors; }
	void set_mou_simplex_factors(vector<double> _factors) { mou_simplex_factors = _factors; }
    bool get_mou_simplex_mutation() const {return mou_simplex_mutation;}
    void set_mou_simplex_mutation(bool _flag) {mou_simplex_mutation = _flag;}
    bool get_mou_use_multigen() const {return mou_use_multigen;}
    void set_mou_use_multigen(bool _flag) {mou_use_multigen = _flag;}
    bool get_mou_shuffle_fixed_pars() const {return mou_shuffle_fixed_pars;}
    void set_mou_shuffle_fixed_pars(bool _flag) {mou_shuffle_fixed_pars = _flag;}

	string get_ies_par_csv()const { return ies_par_csv; }
	void set_ies_par_csv(string _ies_par_csv) { ies_par_csv = _ies_par_csv; }
	string get_ies_obs_csv()const { return ies_obs_csv; }
	void set_ies_obs_csv(string _ies_obs_csv) { ies_obs_csv = _ies_obs_csv; }
	string get_ies_obs_restart_csv() const { return ies_obs_restart_csv; }
	void set_ies_obs_restart_csv(string _ies_obs_restart_csv) { ies_obs_restart_csv = _ies_obs_restart_csv; }
	string get_ies_weights_csv() const {return ies_weights_csv;}
	void set_ies_weights_csv(string _ies_weights_csv) {ies_weights_csv = _ies_weights_csv;}
	string get_ies_par_restart_csv() const { return ies_par_restart_csv; }
	void set_ies_par_restart_csv(string _ies_par_restart_csv) { ies_par_restart_csv = _ies_par_restart_csv; }
	vector<double> get_ies_lam_mults() const { return ies_lam_mults; }
	void set_ies_lam_mults(vector<double> _ies_lam_mults) { ies_lam_mults = _ies_lam_mults; }
	double get_ies_init_lam() const {return ies_init_lam;}
	void set_ies_init_lam(double _ies_init_lam) { ies_init_lam = _ies_init_lam; }
	bool get_ies_use_approx() const { return ies_use_approx; }
	void set_ies_use_approx(bool _ies_use_approx) { ies_use_approx = _ies_use_approx; }
	int get_ies_subset_size() const { return ies_subset_size; }
	void set_ies_subset_size(int _ies_subset_size) { ies_subset_size = _ies_subset_size; }
	double get_ies_reg_factor() const { return ies_reg_factor; }
	void set_ies_reg_factor(double _ies_reg_factor) { ies_reg_factor = _ies_reg_factor; }
	int get_ies_verbose_level() const { return ies_verbose_level; }
	void set_ies_verbose_level(int _ies_verbose_level) { ies_verbose_level = _ies_verbose_level; }
	bool get_ies_use_prior_scaling() const { return ies_use_prior_scaling; }
	void set_ies_use_prior_scaling(bool _ies_use_prior_scaling) { ies_use_prior_scaling = _ies_use_prior_scaling; }
	int get_ies_num_reals() const { return ies_num_reals; }
	void set_ies_num_reals(int _ies_num_reals) { ies_num_reals = _ies_num_reals; }
	double get_ies_bad_phi() const { return ies_bad_phi; }
	void set_ies_bad_phi(double _ies_bad_phi) { ies_bad_phi = _ies_bad_phi; }
	double get_ies_bad_phi_sigma() const { return ies_bad_phi_sigma; }
	void set_ies_bad_phi_sigma(double _ies_bad_phi_sigma) { ies_bad_phi_sigma = _ies_bad_phi_sigma; }
	bool get_ies_include_base() const { return ies_include_base; }
	void set_ies_include_base(bool _ies_include_base) { ies_include_base = _ies_include_base; }
	bool get_ies_use_empirical_prior() const { return ies_use_empirical_prior; }
	void set_ies_use_empirical_prior(bool _ies_use_empirical_prior) { ies_use_empirical_prior = _ies_use_empirical_prior; }
	bool get_ies_group_draws() const { return ies_group_draws; }
	void set_ies_group_draws(bool _ies_group_draws) { ies_group_draws = _ies_group_draws; }
	bool get_ies_enforce_bounds() const { return ies_enforce_bounds; }
	void set_ies_enforce_bounds(bool _ies_enforce_bounds) { ies_enforce_bounds = _ies_enforce_bounds; }

	double get_par_sigma_range() const { return par_sigma_range; }
	void set_par_sigma_range(double _par_sigma_range) { par_sigma_range = _par_sigma_range; }
	bool get_save_binary() const { return save_binary; }
	void set_save_binary(bool _ies_save_binary) { save_binary = _ies_save_binary; }
	string get_ies_localizer() const { return ies_localizer; }
	void set_ies_localizer(string _ies_localizer) { ies_localizer = _ies_localizer; }
	double get_ies_accept_phi_fac() const { return ies_accept_phi_fac; }
	void set_ies_accept_phi_fac(double _acc_phi_fac) { ies_accept_phi_fac = _acc_phi_fac; }
	double get_ies_lambda_inc_fac() const { return ies_lambda_inc_fac; }
	void set_ies_lambda_inc_fac(double _inc_fac) { ies_lambda_inc_fac = _inc_fac; }
	double get_ies_lambda_dec_fac() const { return ies_lambda_dec_fac; }
	void set_ies_lambda_dec_fac(double _dec_fac) { ies_lambda_dec_fac = _dec_fac; }
	bool get_ies_save_lambda_en() const { return ies_save_lambda_en; }
	void set_ies_save_lambda_en(bool _ies_save_lambda_en) { ies_save_lambda_en = _ies_save_lambda_en; }
	string get_ies_subset_how() const { return ies_subset_how; }
	void set_ies_subset_how(string _ies_subset_how) { ies_subset_how = _ies_subset_how; }
	void set_ies_localize_how(string _how) { ies_localize_how = _how; }
	string get_ies_localize_how() const { return ies_localize_how; }
	bool get_ies_localizer_forgive_missing() const {return ies_localizer_forgive_missing;}
	void set_ies_localizer_forgive_missing(bool _flag) {ies_localizer_forgive_missing = _flag;}
	string get_ies_phi_fractions_file() const {return ies_phi_fractions_file;}
	void set_ies_phi_fractions_files(string _file) {ies_phi_fractions_file = _file;}
    bool get_ies_phi_factors_by_real() const {return ies_phi_factors_by_real;}
    void set_ies_phi_factors_by_real(bool _flag) {ies_phi_factors_by_real = _flag;}



	int get_ies_num_threads() const { return ies_num_threads; }
	void set_ies_num_threads(int _threads) { ies_num_threads = _threads; }

	bool get_ies_debug_fail_subset() const { return ies_debug_fail_subset; }
	void set_ies_debug_fail_subset(bool _fail) { ies_debug_fail_subset = _fail; }
	bool get_ies_debug_fail_remainder() const { return ies_debug_fail_remainder; }
	void set_ies_debug_fail_remainder(bool _fail) { ies_debug_fail_remainder = _fail; }
	bool get_ies_debug_bad_phi() const { return ies_debug_bad_phi; }
	void set_ies_debug_bad_phi(bool _fail) { ies_debug_bad_phi = _fail; }
	bool get_ies_debug_upgrade_only() const { return ies_debug_upgrade_only; }
	void set_ies_debug_upgrade_only(bool _flag) { ies_debug_upgrade_only = _flag; }
	bool get_ies_debug_high_subset_phi() const { return ies_debug_high_subset_phi; }
	void set_ies_debug_high_subset_phi(bool _flag) { ies_debug_high_subset_phi = _flag; }
	bool get_ies_debug_high_upgrade_phi() const { return ies_debug_high_upgrade_phi; }
	void set_ies_debug_high_upgrade_phi(bool _flag) { ies_debug_high_upgrade_phi = _flag; }

	bool get_ies_csv_by_reals() const { return ies_csv_by_reals; }
	void set_ies_csv_by_reals(bool _flag) { ies_csv_by_reals = _flag; }
	bool get_ies_autoadaloc() const { return ies_autoadaloc; }
	void set_ies_autoadaloc(bool _flag) { ies_autoadaloc = _flag; }
	double get_ies_autoadaloc_sigma_dist() const { return ies_autoadaloc_sigma_dist; }
	void set_ies_autoadaloc_sigma_dist(double _dist) { ies_autoadaloc_sigma_dist = _dist; }
	bool get_ies_enforce_chglim() const { return ies_enforce_chglim; }
	void set_ies_enforce_chglim(bool _flag) { ies_enforce_chglim = _flag; }
	string get_ies_center_on()const { return ies_center_on; }
	void set_ies_center_on(string _value) { ies_center_on = _value; }
	bool get_ies_no_noise() const { return ies_no_noise; }
	void set_ies_no_noise(bool _flag) { ies_no_noise = _flag; }
	bool get_ies_drop_conflicts() const { return ies_drop_conflicts; }
	void set_ies_drop_conflicts(bool _flag) { ies_drop_conflicts = _flag; }
	bool get_ies_save_rescov() const { return ies_save_rescov; }
	void set_ies_save_rescov(bool _flag) { ies_save_rescov = _flag; }
	double get_ies_pdc_sigma_distance() const { return ies_pdc_sigma_distance; }
	void set_ies_pdc_sigma_distance(double distance) { ies_pdc_sigma_distance = distance; }
	bool get_ies_use_mda() const { return ies_use_mda; }
	void set_ies_use_mda(bool _flag) { ies_use_mda = _flag; }
	double get_ies_mda_init_fac() const { return ies_mda_init_fac; }
	void set_ies_mda_init_fac(double fac) { ies_mda_init_fac = fac; }
	double get_ies_mda_dec_fac() const { return ies_mda_dec_fac; }
	void set_ies_mda_dec_fac(double fac) { ies_mda_dec_fac = fac; }
	string get_ies_loc_type() const { return ies_loc_type; }
	void set_ies_loc_type(string typ) { ies_loc_type = typ; }
	bool get_ies_upgrades_in_memory() const { return ies_upgrades_in_memory; }
	void set_ies_upgrades_in_memory(bool _flag) { ies_upgrades_in_memory = _flag; }
	bool get_ies_ordered_binary() const { return ies_ordered_binary; }
	void set_ies_ordered_binary(bool _flag) { ies_ordered_binary = _flag; }
    double get_ies_multimodal_alpha() const { return ies_multimodal_alpha; }
    void set_ies_multimodal_alpha(double _flag) { ies_multimodal_alpha = _flag; }
    void set_ensemble_output_precision(int prec) { ensemble_output_precision = prec;}
    int get_ensemble_output_precision() const {return ensemble_output_precision;}

	string get_gsa_method() const { return gsa_method; }
	void set_gsa_method(string _m) { gsa_method = _m; }
	bool get_gsa_morris_pooled_obs() const { return gsa_morris_pooled_obs; }
	void set_gsa_morris_pooled_obs(bool _flag) {gsa_morris_pooled_obs = _flag; }
	bool get_gsa_morris_obs_sen() const { return gsa_morris_obs_sen; }
	void set_gsa_morris_obs_sen(bool _flag) { gsa_morris_obs_sen = _flag; }
	double get_gsa_morris_p() const { return gsa_morris_p; }
	void set_gsa_morris_p(int _p) { gsa_morris_p = _p; }
	double get_gsa_morris_r() const { return gsa_morris_r; }
	void set_gsa_morris_r(int _r) { gsa_morris_r = _r; }
	double get_gsa_morris_delta() const { return gsa_morris_delta; }
	void set_gsa_morris_delta(double _d) { gsa_morris_delta = _d; }
	int get_gsa_sobol_samples() const { return gsa_sobol_samples; }
	void set_gsa_sobol_samples(int _s) { gsa_sobol_samples = _s; }
	string get_gsa_sobol_par_dist() const { return gsa_sobol_par_dist; }
	void set_gsa_sobol_par_dist(string _d) { gsa_sobol_par_dist = _d; }

	set<string> get_passed_args() const { return passed_args; }
	map<string, string> get_arg_map()const { return arg_map; }

	void set_enforce_tied_bounds(bool _flag) { enforce_tied_bounds = _flag; }
	bool get_enforce_tied_bounds() const { return enforce_tied_bounds; }

	void set_debug_parse_only(bool _flag) { debug_parse_only = _flag; }
	bool get_debug_parse_only() const { return debug_parse_only; }

	void set_check_tplins(bool _flag) { check_tplins = _flag; }
	bool get_check_tplins() const { return check_tplins; }
	void set_fill_tpl_zeros(bool _flag) { fill_tpl_zeros = _flag; }
	bool get_fill_tpl_zeros() const { return fill_tpl_zeros; }
    void set_tpl_force_decimal(bool _flag) { tpl_force_decimal = _flag; }
    bool get_tpl_force_decimal() const { return tpl_force_decimal; }
	void set_additional_ins_delimiters(string _delims) { additional_ins_delimiters = _delims; }
	string get_additional_ins_delimiters() const { return additional_ins_delimiters; }
	void set_random_seed(int seed) { random_seed = seed; }
	int get_random_seed()const { return random_seed; }

	string get_da_par_cycle_table() const { return da_par_cycle_table; }
	void set_da_par_cycle_table(string _filename) { da_par_cycle_table = _filename; }
	string get_da_obs_cycle_table() const { return da_obs_cycle_table; }
	void set_da_obs_cycle_table(string _filename) { da_obs_cycle_table = _filename; }
	string get_da_weight_cycle_table() const { return da_weight_cycle_table; }
	void set_da_weight_cycle_table(string _filename) { da_weight_cycle_table = _filename; }
	int get_da_hotstart_cycle() const { return da_hotstart_cycle; }
	void set_da_hotstart_cycle(int val) { da_hotstart_cycle = val; }
    int get_da_stop_cycle() const { return da_stop_cycle; }
    void set_da_stop_cycle(int val) { da_stop_cycle = val; }
    bool get_da_use_simulated_states() const {return da_use_simulated_states; }
    void set_da_use_simulated_states(bool _flag) {da_use_simulated_states = _flag;}
    string get_da_noptmax_schedule() const {return da_noptmax_schedule;}
	void set_da_noptmax_schedule(string _file) {da_noptmax_schedule = _file;}


    void set_forgive_unknown_args(bool _flag) { forgive_unknown_args = _flag; }
    bool get_forgive_unknown_args() const { return forgive_unknown_args; }

    bool get_debug_check_par_en_consistency() const { return debug_check_paren_consistency; }
	void set_debug_check_par_en_consistency(bool _flag) { debug_check_paren_consistency = _flag; }
	void set_num_tpl_ins_threads(int num) { num_tpl_ins_threads = num; }
	int get_num_tpl_ins_threads()const { return num_tpl_ins_threads; }

    void set_worker_poll_interval(double _val) { worker_poll_interval = _val; }
    double get_worker_poll_interval() const { return worker_poll_interval; }
    double get_overdue_reched_fac()const { return overdue_reched_fac; }
    void set_overdue_reched_fac(double _val) { overdue_reched_fac = _val; }
    double get_overdue_giveup_fac()const { return overdue_giveup_fac; }
    void set_overdue_giveup_fac(double _val) { overdue_giveup_fac = _val; }
    double get_overdue_giveup_minutes() const { return overdue_giveup_minutes; }
    void set_overdue_giveup_minutes(double overdue_minutes) { overdue_giveup_minutes = overdue_minutes; }
    string get_condor_submit_file() const { return condor_submit_file; }
    void set_condor_submit_file(string _condor_submit_file) { condor_submit_file = _condor_submit_file; }
    void set_panther_agent_restart_on_error(bool _flag) { panther_agent_restart_on_error = _flag; }
    bool get_panther_agent_restart_on_error() const { return panther_agent_restart_on_error; }
    void set_panther_agent_no_ping_timeout_secs(int _timeout_secs) { panther_agent_no_ping_timeout_secs = _timeout_secs; }
    int get_panther_agent_no_ping_timeout_secs() const { return panther_agent_no_ping_timeout_secs; }
    void set_panther_debug_loop(bool _flag) { panther_debug_loop = _flag; }
    bool get_panther_debug_loop() const { return panther_debug_loop; }
    void set_panther_debug_fail_freeze(bool _flag) { panther_debug_fail_freeze = _flag; }
    bool get_panther_debug_fail_freeze() const { return panther_debug_fail_freeze; }
    bool get_panther_echo() const { return panther_echo; }
    void set_panther_echo(bool _flag) { panther_echo = _flag; }
    const vector<string>& get_panther_transfer_on_finish() const {return panther_transfer_on_finish;}
    const vector<string>& get_panther_transfer_on_fail() const {return panther_transfer_on_fail;}
    void set_panther_transfer_on_finish(vector<string> _files) {panther_transfer_on_finish = _files;}
    void set_panther_transfer_on_fail(vector<string> _files) {panther_transfer_on_fail = _files;}



    //bool get_use_da_args() const { return use_da_args; }
	//void set_use_dat_args(bool _flag) { use_da_args = _flag; }

	void set_defaults();
	void summary(ostream& os) const;

private:

	set<string> passed_args;
	set<string> passed_da_ies_args;
	bool use_da_args;
	bool forgive_unknown_args;
	int n_iter_base;
	int n_iter_super;
	int max_n_super;
	double super_eigthres;
	SVD_PACK svd_pack;
	double super_relparmax;
	int max_run_fail;
	int max_super_frz_iter;
	int max_reg_iter;
	int glm_num_reals;
	GLMNormalForm glm_normal_form;
	bool glm_debug_der_fail;
	bool glm_debug_lamb_fail;
	bool glm_debug_real_fail;
	bool glm_accept_mc_phi;
	bool glm_rebase_super;
	bool glm_iter_mc;
	bool glm_debug_high_2nd_iter_phi;

	vector<double> base_lambda_vec;
	vector<double> lambda_scale_vec;
	bool iter_summary_flag;
	bool der_forgive;
	bool uncert;
	vector<string> prediction_names;
	string basejac_filename;
	bool jac_scale;
	string hotstart_resfile;
	string parcov_filename;
	string obscov_filename;
	
	bool tie_by_group;
	bool enforce_tied_bounds;
	bool debug_parse_only;
	bool check_tplins;
	bool fill_tpl_zeros;
	bool tpl_force_decimal;
	string additional_ins_delimiters;

	int random_seed;

	double overdue_reched_fac;
	double overdue_giveup_fac;
	double overdue_giveup_minutes;
	double worker_poll_interval;
	string condor_submit_file;
	int num_tpl_ins_threads;
    int ensemble_output_precision;
	

	string sweep_parameter_csv_file;
	string sweep_output_csv_file;
	bool sweep_forgive;
	int sweep_chunk;
	bool sweep_base_run;
	bool sweep_include_regul_phi;

	GLOBAL_OPT global_opt;
	string moea_name;
	double de_f;
	double de_cr;
	int de_npopulation;
	int de_max_gen;
	bool de_dither_f;

	string opt_obj_func;
	bool opt_coin_log;
	bool opt_skip_final;
	vector<string> opt_dec_var_groups;
	vector<string> opt_external_var_groups;
	vector<string> opt_constraint_groups;
	double opt_risk;
	double opt_direction;
	double opt_iter_tol;
	int opt_recalc_fosm_every;
	double opt_iter_derinc_fac;
	bool opt_include_bnd_pi;
	bool opt_std_weights;
	int opt_stack_size;
	string opt_par_stack;
	string opt_obs_stack;
	string opt_chance_points;

	string sqp_dv_en;
	string sqp_obs_restart_en;
	int sqp_num_reals;
	bool sqp_update_hessian;
	vector<double> sqp_scale_facs;

	int mou_population_size;
	string mou_generator; 
	string mou_dv_population_file;
	string mou_obs_population_restart_file;
	vector<string> mou_objectives;
	int mou_max_archive_size;
	bool mou_risk_obj;
	int mou_verbose_level;
	string mou_env_selector;
	double mou_crossover_prob;
	double mou_mutation_prob;
	string mou_mating_selector;
	double mou_de_f;
	int mou_save_population_every;
	double mou_pso_omega;
	vector<double> mou_pso_social_const;
	vector<double> mou_pso_cognitive_const;
	double mou_pso_alpha;
	double mou_pso_rramp;
	double mou_pso_rfit;
	string mou_pso_dv_bound_restoration;
	vector<double> mou_pso_inertia;
	double mou_pso_vmax_factor;
	double mou_ppd_beta;
	double mou_fit_gamma;
	double mou_fit_epsilon;
	string mou_outer_repo_obs_file;
	int mou_max_nn_search;
	int mou_infill_size;
	double mou_hypervolume_extreme;
	int mou_resample_every;
	string mou_resample_command;
	string mou_population_schedule;
	int mou_simplex_reflections;
	vector<double> mou_simplex_factors;
	bool mou_simplex_mutation;
	bool mou_use_multigen;
	bool mou_shuffle_fixed_pars;

	int ies_subset_size;
	string ies_par_csv;
	string ies_obs_csv;
	string ies_obs_restart_csv;
	string ies_par_restart_csv;
	string ies_weights_csv;
	double ies_init_lam;
	bool ies_use_approx;
	double ies_reg_factor;
	vector<double> ies_lam_mults;
	int ies_verbose_level;
	bool ies_use_prior_scaling;
	int ies_num_reals;	
	double ies_bad_phi;
	double ies_bad_phi_sigma;
	bool ies_include_base;
	bool ies_use_empirical_prior;
	bool ies_group_draws;
	//bool ies_num_reals_passed;
	bool ies_enforce_bounds;
	double par_sigma_range;
	bool save_binary;
	string ies_localizer;
	double ies_accept_phi_fac;
	double ies_lambda_inc_fac;
	double ies_lambda_dec_fac;
	bool ies_save_lambda_en;
	

	map<string, string> arg_map;
	string ies_subset_how;
	string ies_localize_how;
	int ies_num_threads;
	bool ies_debug_fail_subset;
	bool ies_debug_fail_remainder;
	bool ies_debug_bad_phi;
	bool ies_debug_upgrade_only;
	bool ies_debug_high_subset_phi;
	bool ies_debug_high_upgrade_phi;
	bool ies_csv_by_reals;
	bool ies_autoadaloc;
	double ies_autoadaloc_sigma_dist;
	bool ies_enforce_chglim;
	string ies_center_on;
	bool ies_no_noise;
	bool ies_drop_conflicts;
	bool ies_save_rescov;
	double ies_pdc_sigma_distance;
	bool ies_use_mda;
	double ies_mda_init_fac;
	double ies_mda_dec_fac;
	string ies_loc_type;
	bool ies_upgrades_in_memory;
	bool ies_ordered_binary;
	double ies_multimodal_alpha;
	bool ies_localizer_forgive_missing;
	string ies_phi_fractions_file;
	bool ies_phi_factors_by_real;


	// Data Assimilation parameters
	/*string da_mode;
	bool da_use_ies;
	bool use_da = false;
	string da_par_csv;
	int da_num_reals;
	string da_obs_csv;*/
	string da_par_cycle_table;
	string da_obs_cycle_table;
	string da_weight_cycle_table;
	int da_hotstart_cycle;
	int da_stop_cycle;
	bool da_use_simulated_states;
	string da_noptmax_schedule;

	// End of Data Assimilation Parameters

	string gsa_method;
	int gsa_morris_p;
	int gsa_morris_r;
	int gsa_sobol_samples;
	bool gsa_morris_pooled_obs;
	bool gsa_morris_obs_sen;
	double gsa_morris_delta;
	string gsa_sobol_par_dist;


	bool panther_agent_restart_on_error;
	int panther_agent_no_ping_timeout_secs;
	bool panther_debug_loop;
	bool debug_check_paren_consistency;		
	bool panther_debug_fail_freeze;
	bool panther_echo;
	vector<string> panther_transfer_on_finish, panther_transfer_on_fail;

};
//ostream& operator<< (ostream &os, const PestppOptions& val);
ostream& operator<< (ostream &os, const ObservationInfo& val);

class ControlInfo {

public:
	enum PestMode { ESTIMATION, REGUL, PARETO, UNKNOWN };
	double relparmax;
	double facparmax;
	double facorig;
	double phiredswh;
	int noptmax;
	int jacfile;
	int numcom;
	double phiredstp;
	int nphistp;
	int nphinored;
	double relparstp;
	int nrelpar;
	int noptswitch;
	double splitswh;
	PestMode pestmode;
	PestppOptions::ARG_STATUS assign_value_by_key(const string key, const string org_value);
	ControlInfo() { ; }
	void set_defaults();

private:
	set<string> passed_args;
};
ostream& operator<< (ostream& os, const ControlInfo& val);

class SVDInfo {
public:
	int maxsing;
	int eigwrite;
	double eigthresh;
	SVDInfo();
	void set_defaults();
	PestppOptions::ARG_STATUS assign_value_by_key(const string key, const string org_value);
private:
	set<string> passed_args;
};
ostream& operator<< (ostream& os, const SVDInfo& val);

double draw_standard_normal(std::mt19937& rand_gen);
vector<double> uniform_draws(int num_reals, double lower_bound, double upper_bound, std::mt19937& rand_gen);





#endif  /* PEST_DATAS_STRUCTS_H_ */

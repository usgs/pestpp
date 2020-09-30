#ifndef CONSTRAINTS_H_
#define CONSTRAINTS_H_

#include <string>
#include <sstream>
#include <vector>
#include <random>
#include<Eigen/Sparse>
#include <limits>
#include "Pest.h"
#include "logger.h"
#include "FileManager.h"
#include "OutputFileWriter.h"
#include "RunManagerAbstract.h"
#include "Transformable.h"
#include "covariance.h"
#include "Jacobian_1to1.h"
#include "Ensemble.h"

using namespace std;


class OptObjFunc
{
public:
	OptObjFunc(Pest& _pest_scenario, FileManager* _file_mgr_ptr, PerformanceLog& _pfm);
	void initialize(vector<string> _constraint_names, vector<string> _dv_names);
	double get_obj_func_value(Parameters& pars, Observations& obs);
	void report();
	void update_coef_map_from_jacobian(Jacobian& jco);

	const string get_obj_sense() { return obj_sense;  }
	const bool get_use_obs_obj() { return use_obj_obs; }

	const map<string, double> get_obj_func_coef_map() { return obj_func_coef_map;  }


private:
	Pest& pest_scenario;
	FileManager* file_mgr_ptr;
	PerformanceLog& pfm;

	string obj_func_str;
	string obj_obs;
	string obj_sense;
	bool use_obj_obs;
	map<string, double> obj_func_coef_map;
	vector<string> constraint_names;
	vector<string> dv_names;

	void throw_optobjfunc_error(string message);
};

class Constraints
{
	

public:
	enum ConstraintSense { less_than, greater_than, equal_to, undefined };
	enum ConstraintType { deter, pi, fosm, stack };
	Constraints(Pest& _pest_scenario, FileManager* _file_mgr_ptr, OutputFileWriter& _of_wr, PerformanceLog& _pfm);
	void initialize(vector<string>& ctl_ord_dec_var_names, double _dbl_max);

	void presolve_report(int iter, Parameters& current_pars, Observations& current_obs);
	void presolve_chance_report(int iter, Observations& current_obs);
	void postsolve_obs_constraints_report(Observations& old_obs, Observations& new_obs, string tag, int iter,
		map<string, string> status_map=map<string,string>(), map<string, double> price_map=map<string,double>());
	void postsolve_pi_constraints_report(Parameters& old_pars, Parameters& new_pars, int iter,
		map<string, string> status_map = map<string, string>(), map<string, double> price_map = map<string, double>());
	void initial_report();
	
	//queue up chance related runs
	void add_runs(int iter, Parameters& current_pars, Observations& current_obs, RunManagerAbstract* run_mgr_ptr);

	void add_runs(int iter, ParameterEnsemble& current_pe, Observations& current_obs, RunManagerAbstract* run_mgr_ptr);


	//queue up chance related runs at several points in dev var space
	//void add_runs_at_multiple_points(RunManagerAbstract* run_mgr_ptr, ParameterEnsemble& pe, vector<string> only_reals = vector<string>());

	//shift the constraint columns of an obs ensemble in place, potentially accounting for multiple-point dec var evals
	//void risk_shift_obs_en_ip(ObservationEnsemble& oe);
	
	//get the (risk-shifted) distance to constraints (upper and lower bounds)
	pair<vector<double>,vector<double>> get_constraint_bound_vectors(Parameters& current_pars, Observations &current_obs);
	
	//setters
	void set_jco(Jacobian_1to1& _jco) { jco = _jco; }
	void set_initial_constraints_sim(Observations _initial_constraints_sim) { initial_constraints_sim = _initial_constraints_sim;  }
	
	//get fosm-based parameter names
	vector<string> get_fosm_par_names();

	//process chance-related runs
	void process_runs(RunManagerAbstract* run_mgr_ptr,int iter);

	void write_res_files(Observations& constraints_sim, Parameters& pars_and_dec_vars, string tag, int iter);
	
	PriorInformation get_pi_constraints() { return constraints_pi; }

	//get maps of obs and pi constraints that are not satified - the value is the distance to cosntraint RHS
	map<string, double> get_unsatified_pi_constraints(Parameters& par_and_dec_vars, double tol=0.0);
	map<string, double> get_unsatified_obs_constraints(Observations& constraints_sim, double tol=0.0, bool do_shift = true);
	
	//get the number of non-zero Prior info constraint elements
	int get_num_nz_pi_constraint_elements();
	
	//update the chance offset calcs - doesnt make runs, mostly for FOSM calcs
	void update_chance_offsets();
	
	//get the max scale constraint change against upgrade_obs - used for convergence testing
	double get_max_constraint_change(Observations& current_obs, Observations& upgrade_obs);
	
	//get the flags related to chance constraints
	bool get_std_weights() { return std_weights; }
	bool get_use_chance() { return use_chance; }
	bool get_use_fosm() { return use_fosm; }

	//get the dimensions
	int num_obs_constraints() { return ctl_ord_obs_constraint_names.size(); }
	int num_pi_constraints() { return ctl_ord_pi_constraint_names.size(); }
	int num_constraints() { return num_obs_constraints() + num_pi_constraints(); }
	int num_adj_pars() { return adj_par_names.size(); }
	int num_nz_obs() { return nz_obs_names.size(); }

	//get some names
	vector<string> get_obs_constraint_names() { return ctl_ord_obs_constraint_names; }
	vector<string> get_nz_obs_names() { return nz_obs_names; }
	vector<string> get_pi_constraint_names() { return ctl_ord_pi_constraint_names; }
	vector<string> get_adj_par_names() { return adj_par_names; }

	//decide whether it is time to update the chance constraints
	bool should_update_chance(int iter);

	//workout a constraints sense
	static pair<ConstraintSense, string> get_sense_from_group_name(const string& name);

	//get risk-shifted simulated constraint values using current_constraints_sim_ptr
	//Observations get_chance_shifted_constraints();
	//get risk-shifted simulated constraint values using _constraints_sim arg
	Observations get_chance_shifted_constraints(Observations& _constraints_sim);

	ObservationEnsemble get_chance_shifted_constraints(ObservationEnsemble& oe);

private:
	Pest& pest_scenario;
	PerformanceLog& pfm;
	std::mt19937 rand_gen;
	FileManager* file_mgr_ptr;
	OutputFileWriter& of_wr;
	int stack_size;
	bool use_chance;
	bool use_fosm;
	bool std_weights;
	double risk;
	double probit_val;
	double dbl_max;

	Covariance obscov;
	Covariance parcov;
	Jacobian_1to1 jco;
	
	ParameterEnsemble stack_pe;
	ObservationEnsemble stack_oe;
	
	map<string, ObservationEnsemble> stack_oe_map;


	PriorInformation* null_prior = new PriorInformation();
	PriorInformation constraints_pi;

	map<string, ConstraintSense> constraint_sense_map;
	map<string, string> constraint_sense_name;
	map<string, double> prior_const_var;
	map<string, double> post_const_var;
	map<string, double> prior_constraint_offset;
	map<string, double> prior_constraint_stdev;
	map<string, double> post_constraint_offset;
	map<string, double> post_constraint_stdev;

	map<int, int> stack_pe_run_map;
	map<string, map<int, int>> population_stack_pe_run_map;

	vector<string> dec_var_names;
	vector<string> nz_obs_names;
	vector<string> adj_par_names;
	vector<string> ctl_ord_obs_constraint_names;
	vector<string> ctl_ord_pi_constraint_names;

	Observations constraints_obs;
	//Observations* current_constraints_sim_ptr;
	Observations initial_constraints_sim;
	//Parameters* current_pars_and_dec_vars_ptr;
	pair<vector<double>, vector<double>> current_bounds;

	//get the (risk-shifted) residual (distance) vector between constraints RHS and sim arg
	vector<double> get_constraint_residual_vec(Observations& sim);

	//error handlers
	void throw_constraints_error(string message);
	void throw_constraints_error(string message, const vector<string>& messages);
	void throw_constraints_error(string message, const set<string>& messages);

	//get the risk-shift value
	double get_probit();
	double ErfInv2(double x);

	//write residual file
	void write_res_file(Observations& constraints, Parameters& pars_and_dec_vars, string tag, int iter, bool include_chance);

	void process_stack_runs(RunManagerAbstract* run_mgr_ptr, int iter);
	pair<vector<int>,ObservationEnsemble> process_stack_runs(string real_name, int iter, map<int, int> _stack_pe_run_map, 
		RunManagerAbstract* run_mgr_ptr, bool drop_fails=true);

	void save_oe_stack(int iter, string real_name, ObservationEnsemble& _stack_oe);
	void save_pe_stack(int iter, string real_name, ParameterEnsemble& _stack_pe);


};
#endif

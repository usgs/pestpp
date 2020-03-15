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

class Constraints
{
	enum ConstraintSense { less_than, greater_than, equal_to, undefined };
	enum ConstraintType { deter, pi, fosm, stack };

public:
	Constraints(Pest& _pest_scenario, FileManager* _file_mgr_ptr, OutputFileWriter& _of_wr, PerformanceLog& _pfm);
	void initialize(vector<string>& ctl_ord_dec_var_names, Parameters* _current_pars_and_dec_vars_ptr,
		Observations* _current_constraints_sim_ptr, double _dbl_max);

	//map<string, double> get_residual_map(Observations& sim);
	//map<string, double> get_chance_map(Observations& sim);
	//void report(Observations& sim);
	void presolve_report(int iter);
	void presolve_chance_report(int iter);

	void postsolve_obs_constraints_report(Observations& constraints_sim, string tag, int iter,
		map<string, string> status_map=map<string,string>(), map<string, double> price_map=map<string,double>());
	void postsolve_pi_constraints_report(Parameters& pars_and_dec_vars, int iter,
		map<string, string> status_map = map<string, string>(), map<string, double> price_map = map<string, double>());

	
	void initial_report();
	//void chance_report(Observations& sim);
	
	void add_runs(RunManagerAbstract* run_mgr_ptr);
	//void update_obs_and_pi_constraints(Observations& _constraints_sim, Parameters& _par_and_dec_vars);
	pair<vector<double>,vector<double>> get_constraint_bound_vectors();
	void set_jco(Jacobian_1to1& _jco) { jco = _jco; }
	void set_initial_constraints_sim(Observations _initial_constraints_sim) { initial_constraints_sim = _initial_constraints_sim;  }
	vector<string> get_fosm_par_names();

	void process_runs(RunManagerAbstract* run_mgr_ptr,int iter);

	void write_res_files(Observations& constraints_sim, Parameters& pars_and_dec_vars, string tag, int iter);
	
	PriorInformation get_pi_constraints() { return constraints_pi; }

	map<string, double> get_unsatified_pi_constraints(Parameters& par_and_dec_vars, double tol=0.0);
	map<string, double> get_unsatified_obs_constraints(Observations& constraints_sim, double tol=0.0);
	
	int get_num_nz_pi_constraint_elements();
	
	void update_chance_offsets();
	
	double get_max_constraint_change(Observations& upgrade_obs);
	bool get_std_weights() { return std_weights; }
	bool get_use_chance() { return use_chance; }
	bool get_use_fosm() { return use_fosm; }

	int num_obs_constraints() { return ctl_ord_obs_constraint_names.size(); }
	int num_pi_constraints() { return ctl_ord_pi_constraint_names.size(); }
	int num_constraints() { return num_obs_constraints() + num_pi_constraints(); }
	int num_adj_pars() { return adj_par_names.size(); }
	int num_nz_obs() { return nz_obs_names.size(); }

	vector<string> get_obs_constraint_names() { return ctl_ord_obs_constraint_names; }
	vector<string> get_nz_obs_names() { return nz_obs_names; }
	vector<string> get_pi_constraint_names() { return ctl_ord_pi_constraint_names; }
	vector<string> get_adj_par_names() { return adj_par_names; }
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

	pair<ConstraintSense, string> get_sense_from_group_name(const string& name);
	Observations get_chance_shifted_constraints();
	vector<double> get_constraint_residual_vec(Observations& sim);
	
	vector<string> dec_var_names;
	vector<string> nz_obs_names;
	vector<string> adj_par_names;
	vector<string> ctl_ord_obs_constraint_names;
	vector<string> ctl_ord_pi_constraint_names;

	Observations constraints_obs;
	Observations* current_constraints_sim_ptr;
	//Observations current_constraints_chance;
	Observations initial_constraints_sim;
	Parameters* current_pars_and_dec_vars_ptr;
	pair<vector<double>, vector<double>> current_bounds;
	
	

	void throw_constraints_error(string message);
	void throw_constraints_error(string message, const vector<string>& messages);
	void throw_constraints_error(string message, const set<string>& messages);

	double get_probit();
	double ErfInv2(double x);


	void write_res_file(Observations& constraints, Parameters& pars_and_dec_vars, string tag, int iter, bool include_chance);


	

};
#endif

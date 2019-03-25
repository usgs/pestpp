#pragma once
#ifndef SEQUENTIAL_LP_H
#define SEQUENTIAL_LP_H

#include "Pest.h"
#include "Jacobian_1to1.h"
#include "ClpSimplex.hpp"
#include "RunManagerAbstract.h"
#include "TerminationController.h"
#include "covariance.h"
#include "FileManager.h"
#include "OutputFileWriter.h"
#include "Transformable.h"
#include "ModelRunPP.h"

class sequentialLP
{
	enum ConstraintSense {less_than,greater_than,equal_to,undefined};
public:
	sequentialLP(Pest &_pest_scenario, RunManagerAbstract* _run_mgr_ptr,
		         Covariance &_parcov, FileManager* _file_mgr_ptr, OutputFileWriter of_wr);
	void initialize_and_check();
	void solve();

	~sequentialLP();
	//ModelRun get_optimum_run() { return optimum_run; }

private:
	string obj_func_str;
	bool use_chance;
	bool terminate;
	bool super_secret_option;
	bool use_obj_obs;
	bool std_weights;
	string obj_obs;
	int slp_iter;

	double* dec_var_lb;
	double* dec_var_ub;
	double* constraint_lb;
	double* constraint_ub;
	double* ctl_ord_obj_func_coefs;
	const double* row_price;
	double risk;
	double iter_derinc_fac;
	double obj_best;
	double obj_init;
	double probit_val;

	string obj_sense;
	ClpSimplex model;
	CoinMessageHandler coin_hr;
	FILE* coin_log_ptr;
	Jacobian_1to1 jco;
	OutputFileWriter of_wr;

	vector<double> probit_inputs = { 0.0010,0.0100,0.0200,0.0300,0.0400,0.0500,0.0600,0.0700,0.0800,0.0900,
		0.1000,0.1100,0.1200,0.1300,0.1400,0.1500,0.1600,0.1700,0.1800,0.1900,
		0.2000,0.2100,0.2200,0.2300,0.2400,0.2500,0.2600,0.2700,0.2800,0.2900,
		0.3000,0.3100,0.3200,0.3300,0.3400,0.3500,0.3600,0.3700,0.3800,0.3900,
		0.4000,0.4100,0.4200,0.4300,0.4400,0.4500,0.4600,0.4700,0.4800,0.4900,
		0.5000,0.5100,0.5200,0.5300,0.5400,0.5500,0.5600,0.5700,0.5800,0.5900,
		0.6000,0.6100,0.6200,0.6300,0.6400,0.6500,0.6600,0.6700,0.6800,0.6900,
		0.7000,0.7100,0.7200,0.7300,0.7400,0.7500,0.7600,0.7700,0.7800,0.7900,
		0.8000,0.8100,0.8200,0.8300,0.8400,0.8500,0.8600,0.8700,0.8800,0.8900,
		0.9000,0.9100,0.9200,0.9300,0.9400,0.9500,0.9600,0.9700,0.9800,0.9990 };
	vector<double> probit_outputs = { -3.0902,-2.3263,-2.0537,-1.8808,-1.7507,-1.6449,-1.5548,-1.4758,-1.4051,-1.3408,
		-1.2816,-1.2265,-1.1750,-1.1264,-1.0803,-1.0364,-0.9945,-0.9542,-0.9154,-0.8779,
		-0.8416,-0.8064,-0.7722,-0.7388,-0.7063,-0.6745,-0.6433,-0.6128,-0.5828,-0.5534,
		-0.5244,-0.4959,-0.4677,-0.4399,-0.4125,-0.3853,-0.3585,-0.3319,-0.3055,-0.2793,
		-0.2533,-0.2275,-0.2019,-0.1764,-0.1510,-0.1257,-0.1004,-0.0753,-0.0502,-0.0251,
		0.0000,0.0251,0.0502,0.0753,0.1004,0.1257,0.1510,0.1764,0.2019,0.2275,
		0.2533,0.2793,0.3055,0.3319,0.3585,0.3853,0.4125,0.4399,0.4677,0.4959,
		0.5244,0.5534,0.5828,0.6128,0.6433,0.6745,0.7063,0.7388,0.7722,0.8064,
		0.8416,0.8779,0.9154,0.9542,0.9945,1.0364,1.0803,1.1264,1.1750,1.2265,
		1.2816,1.3408,1.4051,1.4758,1.5548,1.6449,1.7507,1.8808,2.0537,3.0902 };


	map<ClpSimplex::Status, string> status_name_map = { {ClpSimplex::Status::atLowerBound,"at lower bound"},
	{ ClpSimplex::Status::atUpperBound,"at upper bound"},{ClpSimplex::Status::basic,"basic"},
	{ ClpSimplex::Status::isFree,"free"},{ ClpSimplex::Status::isFixed,"fixed"}};

	map<string, double> obj_func_coef_map;
	map<string, ConstraintSense> constraint_sense_map;
	map <string, string> constraint_sense_name;
	map<string, double> prior_constraint_stdev;
	map<string, double> post_constraint_stdev;
	map<string, double> prior_constraint_offset;
	map<string, double> post_constraint_offset;
	map<string, double> prior_const_var;
	map<string, double> post_const_var;
	//map<string, map<string, double>> pi_constraint_factors;
	//map<string, double> pi_constraint_rhs;

	vector<double> iter_obj_values;
	vector<string> ctl_ord_dec_var_names;
	vector<string> ctl_ord_obs_constraint_names;
	vector<string> ctl_ord_pi_constraint_names;
	vector<string> ctl_ord_ext_var_names;
	vector<string> nz_obs_names;
	vector<string> adj_par_names;

	PriorInformation* null_prior = new PriorInformation();
	Parameters all_pars_and_dec_vars;
	Parameters all_pars_and_dec_vars_initial;
	Parameters all_pars_and_dec_vars_best;
	ParamTransformSeq par_trans;
	Observations constraints_obs;
	Observations constraints_sim;
	Observations constraints_fosm;
	Observations constraints_sim_initial;
	PriorInformation constraints_pi;
	Observations obj_func_obs;
	ObservationInfo obj_func_info;
	Pest pest_scenario;
	RunManagerAbstract* run_mgr_ptr;
	Covariance parcov;
	Covariance obscov;
	FileManager* file_mgr_ptr;
	//OutputFileWriter* out_wtr_ptr;

	int num_dec_vars() { return ctl_ord_dec_var_names.size(); }
	int num_obs_constraints() { return ctl_ord_obs_constraint_names.size(); }
	int num_pi_constraints() { return ctl_ord_pi_constraint_names.size(); }
	int num_constraints() { return num_obs_constraints() + num_pi_constraints(); }
	int num_adj_pars() { return adj_par_names.size(); }
	int num_nz_obs() { return nz_obs_names.size(); }

	void build_dec_var_bounds();

	//get the interpolated probit value
	double get_probit();

	//report the infeasible info
	void iter_infeasible_report();

	//get the number of non zero elements in the prior information constraints
	int num_nz_pi_constraint_elements();

	//parse the obs or pi group name to get the constraint sense
	pair<ConstraintSense,string> get_sense_from_group_name(const string &name);

	//solve the current LP problem
	void iter_solve();

	//report initial conditions to rec file
	void initial_report();

	//report the constraint info before the solving the current LP problem
	void presolve_constraint_report();

	//report the fosm chance constraint info before solving the current LP problem
	void presolve_fosm_report();

	//report dec var info the newly solved LP solution.  returns the current and new obj func
	pair<double,double> postsolve_decision_var_report(Parameters &upgrade_pars);

	//report the current and newly solved LP constraint info
	void postsolve_constraint_report(Observations &upgrade_obs, Parameters &upgrade_pars);

	//check that all constraints and dec vars are satified
	pair < map < string, double > , map<string,double >> postsolve_check(Observations &upgrade_obs, Parameters &upgrade_pars);

	//prepare for LP solution, including filling response matrix
	void iter_presolve();

	//run the model with dec var values from the newly solved LP problem
	bool make_upgrade_run(Parameters &upgrade_pars, Observations &upgrade_obs);

	//process the LP solve, including check for convergence
	void iter_postsolve();

	//convert the jacobian to a coin packed matrix instance
	CoinPackedMatrix jacobian_to_coinpackedmatrix();

	//convert the constraint info from Transformable to double*
	void build_constraint_bound_arrays();

	//error handlers
	void throw_sequentialLP_error(string message);
	void throw_sequentialLP_error(string message,const vector<string> &messages);
	void throw_sequentialLP_error(string message, const set<string> &messages);

	//get the current constraint residual vector
	vector<double> get_constraint_residual_vec();

	//get a residual vector comparing constraints_obs and sim_vals
	vector<double> get_constraint_residual_vec(Observations &sim_vals);

	//set the double* obj_func array
	void build_obj_func_coef_array();

	//calc FOSM-based chance constraint offsets
	void calc_chance_constraint_offsets();

	double obj_func_report();
};




#endif

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
#include "constraints.h"

class sequentialLP
{
	enum ConstraintSense {less_than,greater_than,equal_to,undefined};
public:
	sequentialLP(Pest &_pest_scenario, RunManagerAbstract* _run_mgr_ptr,
		         Covariance &_parcov, FileManager* _file_mgr_ptr, OutputFileWriter of_wr,
				PerformanceLog& _pfm);
	void initialize_and_check();
	void solve();

	~sequentialLP();

private:
	std::mt19937 rand_gen;
	PerformanceLog& pfm;
	//string obj_func_str;
	bool terminate;
	bool super_secret_option;
	//bool use_obj_obs;
	//string obj_obs;
	int slp_iter;

	double* dec_var_lb;
	double* dec_var_ub;
	double* ctl_ord_obj_func_coefs;
	const double* row_price;
	double iter_derinc_fac;
	double obj_best;
	double obj_init;

	//string obj_sense;
	ClpSimplex model;
	CoinMessageHandler coin_hr;
	FILE* coin_log_ptr;
	Jacobian_1to1 jco;
	OutputFileWriter of_wr;

	Constraints constraints;
	OptObjFunc optobjfunc;

	map<ClpSimplex::Status, string> status_name_map = { {ClpSimplex::Status::atLowerBound,"at lower bound"},
	{ ClpSimplex::Status::atUpperBound,"at upper bound"},{ClpSimplex::Status::basic,"basic"},
	{ ClpSimplex::Status::isFree,"free"},{ ClpSimplex::Status::isFixed,"fixed"}};

	//map<string, double> obj_func_coef_map;

	vector<double> iter_obj_values;
	vector<string> dv_names;
	vector<string> ext_dv_names;

	PriorInformation* null_prior = new PriorInformation();
	Parameters current_pars;
	Parameters initial_pars;
	Parameters best_pars;
	ParamTransformSeq par_trans;
	Observations current_constraints_sim;
	
	Observations obj_func_obs;
	ObservationInfo obj_func_info;
	Pest pest_scenario;
	RunManagerAbstract* run_mgr_ptr;
	Covariance parcov;
	Covariance obscov;
	FileManager* file_mgr_ptr;
	//OutputFileWriter* out_wtr_ptr;

	int num_dec_vars() { return dv_names.size(); }

	void build_dec_var_bounds();

	//report the infeasible info
	void iter_infeasible_report();

	//parse the obs or pi group name to get the constraint sense
	pair<ConstraintSense,string> get_sense_from_group_name(const string &name);

	//solve the current LP problem
	void iter_solve();

	//report initial conditions to rec file
	void initial_report();

	//report dec var info the newly solved LP solution.  returns the current and new obj func
	pair<double,double> postsolve_decision_var_report(Parameters &upgrade_pars);

	//check that all constraints and dec vars are satified
	map<string, double> get_out_of_bounds_dec_vars(Parameters &upgrade_pars);

	//prepare for LP solution, including filling response matrix
	void iter_presolve();

	//run the model with dec var values from the newly solved LP problem
	bool make_upgrade_run(Parameters &upgrade_pars, Observations &upgrade_obs);

	//process the LP solve, including check for convergence
	void iter_postsolve();

	//convert the jacobian to a coin packed matrix instance
	CoinPackedMatrix jacobian_to_coinpackedmatrix();

	//error handlers
	void throw_sequentialLP_error(string message);
	void throw_sequentialLP_error(string message,const vector<string> &messages);
	void throw_sequentialLP_error(string message, const set<string> &messages);

	//set the double* obj_func array
	void build_obj_func_coef_array();

	//double obj_func_report();

};




#endif

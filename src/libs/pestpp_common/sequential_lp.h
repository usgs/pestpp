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

	// -- the SLP loop, as separately callable steps ----------------------------------------
	//
	// solve() is the in-tree composition of these, so there is one code path rather than an
	// in-tree one and an API one that drift apart. A caller drives it as
	//     while (!slp.should_terminate()) slp.solve_iteration();
	//     slp.finalize();

	/// ONE LP iteration: presolve, LP solve, postsolve, and the objective-sequence bookkeeping.
	/// Advances slp_iter by exactly one UNLESS this iteration converged, in which case the
	/// counter stays put - which is what the loop did before, and what the .par/.rei file
	/// names are keyed on. Returns whether the caller should keep going.
	bool solve_iteration();

	/// Whether the loop is over: convergence, noptmax reached, or 'pest.stp' seen. False
	/// before the first solve_iteration(), so the while-loop above enters.
	bool should_terminate() const { return solve_finished; }

	/// End-of-run reporting: the objective function sequence, the best value, and (unless
	/// ++opt_skip_final) the final model run at the best decision variables.
	void finalize();

	void solve();

	/// Which SLP iteration the tool is on. 0 before the first solve_iteration().
	int get_slp_iter() const { return slp_iter; }

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

	/// Whether solve_begin() has run. The slp.iobj.csv stream and slp_iter used to be set up
	/// at the top of solve(); a caller driving iterations itself never goes through there, so
	/// the first solve_iteration() does it instead - lazily, to keep the file's creation time
	/// where it was rather than moving it into the constructor.
	bool solve_begun;
	/// Sticky: the loop has decided to stop. Separate from `terminate`, which means only
	/// "converged" - noptmax and 'pest.stp' also end the loop and neither sets that flag.
	bool solve_finished;

	/// Open slp.iobj.csv, write its header, and start the iteration counter.
	void solve_begin();

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
	// reference, not a copy: the Constraints and OptObjFunc members below are built from
	// the caller's Pest, so holding a copy here gave sequentialLP a second, divergent
	// PestppOptions. Binding the reference keeps one source of truth.
	Pest& pest_scenario;
	RunManagerAbstract* run_mgr_ptr;
	Covariance parcov;
	Covariance obscov;
	FileManager* file_mgr_ptr;
	//OutputFileWriter* out_wtr_ptr;

	int num_dec_vars() { return dv_names.size(); }

	void build_dec_var_bounds();

	//report the infeasible info
	void iter_infeasible_report();
	
	//solve the current LP problem
	void iter_solve();

	//report initial conditions to rec file
	void initial_report();

	//report dec var info the newly solved LP solution.  returns the current and new obj func
	pair<double,double> postsolve_decision_var_report(Parameters &upgrade_pars);

	//check that all constraints and dec vars are satisfied
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

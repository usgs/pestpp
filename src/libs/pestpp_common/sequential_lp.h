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

	// -- the slp loop, broken into steps you can call yourself ------------------------------
	//
	// solve() just calls these in order, so there is one code path instead of an in-tree one and
	// an api one that end up drifting apart. you drive it like this:
	//     while (!slp.should_terminate()) slp.solve_iteration();
	//     slp.finalize();

	/// one lp iteration - presolve, lp solve, postsolve, and keeping track of the objective
	/// values. bumps slp_iter by one unless this iteration converged, in which case the counter
	/// stays where it is. that is what the loop did before and the .par/.rei file names are
	/// keyed on it. returns whether you should keep going.
	bool solve_iteration();

	/// whether the loop is done - converged, hit noptmax, or saw 'pest.stp'. false before the
	/// first solve_iteration() so the while loop above gets in.
	bool should_terminate() const { return solve_finished; }

	/// end of run reporting - the objective function sequence, the best value, and unless you
	/// set ++opt_skip_final, one last model run at the best decision variables.
	void finalize();

	void solve();

	/// which slp iteration we are on. 0 before the first solve_iteration().
	int get_slp_iter() const { return slp_iter; }

	// -- current state, for somebody driving the loop themselves -----------------------------
	//
	// all const refs. these are the tool's own objects and you read them between iterations.
	// handing back copies would be safe enough but then somebody modifies the copy and wonders
	// why nothing happened.

	/// the decision variables as they are right now. this is opt's parameter vector - the thing
	/// being optimized. until you could get at this you could drive opt all the way to an answer
	/// and then not be able to read the answer.
	const Parameters& get_current_pars() const { return current_pars; }

	/// the decision variables at the best objective we have seen, which usually isnt the current
	/// one - the last slp iteration might have come out worse and still be what we are sitting
	/// on.
	const Parameters& get_best_pars() const { return best_pars; }

	/// the simulated constraint values that go with get_current_pars().
	const Observations& get_current_constraints_sim() const { return current_constraints_sim; }

	/// the objective values, oldest first. element 0 is the objective at the starting decision
	/// variables, so the size is (iterations + 1), not the iteration count. this is how you see
	/// opt converging - same idea as the phi sequence in the ensemble tools - and it is the
	/// sequence the rec file prints at the end.
	const std::vector<double>& get_obj_values() const { return iter_obj_values; }

	/// best and starting objective function values.
	double get_obj_best() const { return obj_best; }
	double get_obj_init() const { return obj_init; }

	/// which parameters are decision variables. not every adjustable parameter is - the rest are
	/// the uncertain ones that feed the chance constraints. so "the parameters" and "the
	/// decision variables" are two different questions with two different answers.
	const std::vector<std::string>& get_dv_names() const { return dv_names; }
	/// decision variables that came from an external file instead of the control file.
	const std::vector<std::string>& get_ext_dv_names() const { return ext_dv_names; }

	/// the chance stuff, so opt gets the same stack/risk calls mou and sqp have through the c
	/// api. public because the adapter needs it - this class still owns the object.
	Constraints* get_constraints_ptr() { return &constraints; }

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

	/// whether solve_begin() has run yet. the slp.iobj.csv stream and slp_iter used to get set
	/// up at the top of solve(), but somebody driving the iterations themselves never goes
	/// through there. so the first solve_iteration() does it instead - done lazily so the file
	/// still gets created at the same point instead of in the constructor.
	bool solve_begun;
	/// sticky - the loop has decided to quit. separate from `terminate`, which only means
	/// converged. noptmax and 'pest.stp' also end the loop and neither one sets that flag.
	bool solve_finished;

	/// open slp.iobj.csv, write the header, start the iteration counter.
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
	/// declared before jco and constraints, and the order matters. Jacobian_1to1 keeps an
	/// OutputFileWriter* and Constraints keeps an OutputFileWriter&, and both point at this
	/// member. members get initialized in declaration order so of_wr has to exist first. before,
	/// they were bound to the constructor's argument instead, which left two references to a
	/// stack object that got destroyed as soon as the constructor returned.
	OutputFileWriter of_wr;
	Jacobian_1to1 jco;

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

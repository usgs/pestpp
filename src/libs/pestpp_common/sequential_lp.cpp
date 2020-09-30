#include "sequential_lp.h"
#include "Pest.h"
#include "Jacobian_1to1.h"
#include "ClpSimplex.hpp"
#include "RunManagerAbstract.h"
#include "TerminationController.h"
#include "covariance.h"
#include "linear_analysis.h"
#include "FileManager.h"
#include "OutputFileWriter.h"
#include <Eigen/Sparse>
#include "CoinFinite.hpp"
#include "ClpPresolve.hpp"
#include <iomanip>
#include "utilities.h"

sequentialLP::sequentialLP(Pest &_pest_scenario, RunManagerAbstract* _run_mgr_ptr,
	Covariance &_parcov, FileManager* _file_mgr_ptr, OutputFileWriter _of_wr, PerformanceLog& _pfm) 
	: pest_scenario(_pest_scenario), run_mgr_ptr(_run_mgr_ptr),
	parcov(_parcov), file_mgr_ptr(_file_mgr_ptr),jco(*_file_mgr_ptr,_of_wr), of_wr(_of_wr), pfm(_pfm),
	constraints(_pest_scenario,_file_mgr_ptr,_of_wr,_pfm), optobjfunc(_pest_scenario, _file_mgr_ptr,_pfm)
{
	rand_gen = std::mt19937(pest_scenario.get_pestpp_options().get_random_seed());

	try
	{
		initialize_and_check();
	}
	catch (const runtime_error& error)
	{
		cout << "error initializing sequentialLP process: " << error.what() << endl;
		exit(1);
	}

}

sequentialLP::~sequentialLP()
{
	delete[] dec_var_lb;
	delete[] dec_var_ub;
	//delete[] constraint_lb;
	//delete[] constraint_ub;
	delete[] ctl_ord_obj_func_coefs;
	//delete[] row_price;
}

void sequentialLP::throw_sequentialLP_error(string message,const vector<string> &messages)
{
	stringstream ss;
	for (auto &mess : messages)
		ss << mess + ',';
	throw_sequentialLP_error(message + ss.str());
}

void sequentialLP::throw_sequentialLP_error(string message, const set<string> &messages)
{
	stringstream ss;
	for (auto &mess : messages)
		ss << mess + ',';
	throw_sequentialLP_error(message + ss.str());
}

void sequentialLP::throw_sequentialLP_error(string message)
{
	string error_message = "error in sequentialLP process: " + message;
	file_mgr_ptr->rec_ofstream() << error_message << endl;
	file_mgr_ptr->close_file("rec");
	cout << endl << endl << error_message << endl << endl;
	throw runtime_error(error_message);
}

void sequentialLP::initial_report()
{
	ofstream &f_rec = file_mgr_ptr->rec_ofstream();
	f_rec << endl << "  -------------------------------------------------------------" << endl;
	f_rec << "  ---  sequential linear programming problem information  ---  " << endl;
	f_rec << "  -------------------------------------------------------------" << endl << endl;

	f_rec << "-->number of iterations of sequential linear programming (noptmax): " << pest_scenario.get_control_info().noptmax << endl;

	f_rec << "-->objective function sense (direction): " << optobjfunc.get_obj_sense() << endl;


	f_rec << "-->number of decision variable: " << num_dec_vars() << endl;
	f_rec << "-->number of observation constraints: " << constraints.num_obs_constraints() << endl;
	f_rec << "-->number of prior information constraints: " << constraints.num_pi_constraints() << endl;
	if (iter_derinc_fac != 1.0)
		f_rec << "-->iteration DERINC reduction factor (++opt_iter_derinc_fac): " << iter_derinc_fac << endl;
	
	string bj = pest_scenario.get_pestpp_options().get_basejac_filename();
	if (!bj.empty())
		f_rec << "-->start with existing jacobian: " << bj << endl;
	string hotstart = pest_scenario.get_pestpp_options().get_hotstart_resfile();
	if (!hotstart.empty())
		f_rec << "-->hot start with residual file: " << hotstart << endl;
	bool sf = pest_scenario.get_pestpp_options().get_opt_skip_final();
	super_secret_option = false;
	if (sf)
	{
		f_rec << "-->skipping final, optimal model run" << endl;
		if ((!bj.empty()) && (pest_scenario.get_control_info().noptmax == 1))
		{
			f_rec << "-->super secrect option to skip final run and upgrade run activated..." << endl;
			super_secret_option = true;
		}

	}
	f_rec << endl << endl << "  ---  decision variables active in SLP  ---  " << endl;
	
	optobjfunc.report();

	if (!optobjfunc.get_use_obs_obj())
	{
		obj_init = optobjfunc.get_obj_func_value(initial_pars,current_constraints_sim);
		f_rec << endl << "  ---  objective function value (using initial dec var values): " << obj_init << endl << endl;
		cout << endl << "  ---  objective function value (using initial dec var values): " << obj_init << endl << endl;

	}
	f_rec << " note: bound and initial value info reported in 'parameter data' section" << endl << endl;

	constraints.initial_report();

	return;
}


//double sequentialLP::obj_func_report()
//{
//	vector<string> missing;
//	ofstream &f_rec = file_mgr_ptr->rec_ofstream();
//	map<string, double>::iterator end = obj_func_coef_map.end();
//	f_rec << setw(20) << left << "name" << setw(25) << "obj func coefficient" << endl;
//	for (auto &name : ctl_ord_dec_var_names)
//	{
//		f_rec << setw(20) << left << name;
//		if (obj_func_coef_map.find(name) != end)
//		{
//			f_rec << setw(25) << obj_func_coef_map.at(name) << endl;
//		}
//		else
//		{
//			f_rec << setw(25) << "not listed" << endl;
//			missing.push_back(name);
//		}
//	}
//	double obj = 0.0;
//	for (auto &p : ctl_ord_dec_var_names)
//	{
//		obj += all_pars_and_dec_vars_initial[p] * obj_func_coef_map[p];
//	}
//	if (missing.size() > 0)
//	{
//		f_rec << endl << endl << "WARNING: the following decision variables have '0.0' objective function coef:" << endl;
//		cout << endl << endl << "WARNING: the following decision variables have '0.0' objective function coef:" << endl;
//
//		for (auto &name : missing)
//		{
//			f_rec << "    " << name << endl;
//			f_rec << "    " << name << endl;
//		}
//	}
//	return obj;
//
//}

map<string,double> sequentialLP::get_out_of_bounds_dec_vars(Parameters &upgrade_pars)
{
	double opt_tol = pest_scenario.get_pestpp_options().get_opt_iter_tol();
	Parameters ubnd = pest_scenario.get_ctl_parameter_info().get_up_bnd(dv_names);
	Parameters lbnd = pest_scenario.get_ctl_parameter_info().get_low_bnd(dv_names);
	double sim_val;
	map<string, double> invalid_dec_vars;
	for (auto &name : dv_names)
	{
		sim_val = upgrade_pars[name];
		if (sim_val > ubnd[name])
			invalid_dec_vars[name] = sim_val - ubnd[name];
		else if (sim_val < lbnd[name])
			invalid_dec_vars[name] = lbnd[name] - sim_val;
	}

	return invalid_dec_vars;
}

pair<double,double> sequentialLP::postsolve_decision_var_report(Parameters &upgrade_pars)
{
	ofstream &f_rec = file_mgr_ptr->rec_ofstream();
	const double *reduced_cost = model.getReducedCost();
	f_rec << endl << endl << "     decision variable information at end of SLP iteration " << slp_iter << endl << endl;
	f_rec << setw(20) << left << "name" << right << setw(15) << "current" << setw(15)  << "new";
	f_rec << setw(15) << "objfunc coef" << setw(15) << "cur contrib" << setw(15) << "new contrib" << setw(15) << "reduced cost";
	f_rec << setw(25) << "simplex status" << endl;
	string name;
	ClpSimplex::Status status;
	double obj_coef, cur_val, new_val, upgrade, init_val;
	double cur_obj=0.0, new_obj=0.0;
	Parameters actual_pars = upgrade_pars;
	for (int i = 0; i < num_dec_vars(); ++i)
	{
		name = dv_names[i];
		status = model.getColumnStatus(i);
		obj_coef = ctl_ord_obj_func_coefs[i];
		init_val = initial_pars[name];
		cur_val = current_pars[name];

		new_val = upgrade_pars[name];
		actual_pars.update_rec(name, new_val);
		f_rec << setw(20) << left << name;
		f_rec << setw(15) << right << cur_val;
		f_rec << setw(15) << new_val;
		f_rec << setw(15) << obj_coef;
		f_rec << setw(15) << cur_val * obj_coef;
		f_rec << setw(15) << new_val * obj_coef;
		f_rec << setw(15) << reduced_cost[i];
		f_rec << setw(25) << status_name_map[status] << endl;
		cur_obj += cur_val * obj_coef;
		new_obj += new_val * obj_coef;

	}
	stringstream ss;
	ss << slp_iter << ".par";

	of_wr.write_par(file_mgr_ptr->open_ofile_ext(ss.str()),actual_pars,*par_trans.get_offset_ptr(),*par_trans.get_scale_ptr());
	file_mgr_ptr->close_file(ss.str());
	of_wr.write_par(file_mgr_ptr->open_ofile_ext("par"), actual_pars, *par_trans.get_offset_ptr(), *par_trans.get_scale_ptr());
	file_mgr_ptr->close_file("par");
	return pair<double,double>(cur_obj,new_obj);
}

pair<sequentialLP::ConstraintSense,string> sequentialLP::get_sense_from_group_name(const string &name)
{

	if ((name.compare(0, 2, "L_") == 0) || (name.compare(0, 4, "LESS")==0))
		return pair<ConstraintSense,string>(ConstraintSense::less_than,"less_than");
	else if ((name.compare(0, 2, "G_") == 0) || (name.compare(0, 7, "GREATER")==0))
		return pair<ConstraintSense, string>(ConstraintSense::greater_than,"greater_than");
	else if ((name.compare(0, 2, "N_") == 0) || (name.compare(0, 2, "E_") == 0) || (name.compare(0, 5, "EQUAL")==0))
		return pair<ConstraintSense,string>(ConstraintSense::equal_to,"equal_to");
	else
		return pair<ConstraintSense,string>(ConstraintSense::undefined,"undefined");
}


void sequentialLP::initialize_and_check()
{
	ofstream &f_rec = file_mgr_ptr->rec_ofstream();
	//TODO: handle restart condition

	if (pest_scenario.get_control_info().pestmode != ControlInfo::PestMode::ESTIMATION)
	{
		string mess = "'pestmode' != 'estimation'.  pestpp-opt really only operates in kind-of 'estimation' mode.  ignoring";
		cout << endl << mess << endl;
		f_rec << endl << mess << endl;
	}


	if (pest_scenario.get_control_info().noptmax < 1)
		throw_sequentialLP_error("noptmax must be greater than 0");

	coin_log_ptr = fopen(file_mgr_ptr->build_filename("coin_log").c_str(), "w");
	coin_hr = CoinMessageHandler(coin_log_ptr);
	model.passInMessageHandler(&coin_hr);

	terminate = false;

	iter_derinc_fac = pest_scenario.get_pestpp_options().get_opt_iter_derinc_fac();
	if ((iter_derinc_fac > 1.0) || (iter_derinc_fac <= 0.0))
		throw_sequentialLP_error("++opt_iter_derinc_fac must be greater than 0.0 and less than or equal to 1.0");

	//-----------------------------
	//  ---  decision vars  ---
	//-----------------------------

	current_pars = pest_scenario.get_ctl_parameters();
	initial_pars = Parameters(current_pars);
	best_pars = Parameters(current_pars);
	par_trans = pest_scenario.get_base_par_tran_seq();
	//set ordered dec var name vec
	//and check for illegal parameter transformations
	vector<string> dec_var_groups = pest_scenario.get_pestpp_options().get_opt_dec_var_groups();
	vector<string> ext_var_groups = pest_scenario.get_pestpp_options().get_opt_ext_var_groups();
	dec_var_groups.insert(dec_var_groups.begin(), ext_var_groups.begin(), ext_var_groups.end());
	dv_names.clear();
	//if the ++opt_dec_var_groups arg was passed
	if (dec_var_groups.size() != 0)
	{
		//first make sure all the groups are actually listed in the control file
		vector<string> missing;
		vector<string> pst_groups = pest_scenario.get_ctl_ordered_par_group_names();
		vector<string>::iterator end = pst_groups.end();
		vector<string>::iterator start = pst_groups.begin();
		for (auto grp : dec_var_groups)
			if (find(start, end, grp) == end)
				missing.push_back(grp);
		if (missing.size() > 0)
			throw_sequentialLP_error("the following ++opt_dec_var_groups were not found: ", missing);

		//find the parameter in the dec var groups
		ParameterGroupInfo pinfo = pest_scenario.get_base_group_info();
		string group;
		end = dec_var_groups.end();
		start = dec_var_groups.begin();
		for (auto &par_name : pest_scenario.get_ctl_ordered_par_names())
		{
			group = pinfo.get_group_name(par_name);
			if (find(start, end, group) != end)
			{
				dv_names.push_back(par_name);
				//check if this is an ext var
				if (find(ext_var_groups.begin(), ext_var_groups.end(), group) != ext_var_groups.end())
					ext_dv_names.push_back(par_name);
			}
		}

		if (num_dec_vars() == 0)
			throw_sequentialLP_error("no decision variables found in groups: ", dec_var_groups);
	}
	//if not ++opt_dec_var_names was passed, use all parameter as decision variables
	else
		dv_names = pest_scenario.get_ctl_ordered_par_names();

	//if any decision vars have a transformation that is not allowed
	vector<string> problem_trans;
	for (auto &name : dv_names)
		if (pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(name)->tranform_type != ParameterRec::TRAN_TYPE::NONE)
			problem_trans.push_back(name);
	if (problem_trans.size() > 0)
		throw_sequentialLP_error("the following decision variables don't have 'none' type parameter transformation: ", problem_trans);

	current_constraints_sim = pest_scenario.get_ctl_observations();

 	constraints.initialize(dv_names, COIN_DBL_MAX);

	
	//--------------------------------
	//  ---  objective function  ---
	//--------------------------------

	//initialize the objective function
	//obj_func_str = pest_scenario.get_pestpp_options().get_opt_obj_func();
	//obj_sense = (pest_scenario.get_pestpp_options().get_opt_direction() == 1) ? "minimize" : "maximize";


	//check if the obj_str is an observation
	//use_obj_obs = false;
	//if (pest_scenario.get_ctl_observations().find(obj_func_str) != pest_scenario.get_ctl_observations().end())
	//{
	//	use_obj_obs = true;
	//	obj_obs = obj_func_str;
	//	//check
	//	vector<string> cnames = constraints.get_obs_constraint_names();
	//	set<string> names(cnames.begin(), cnames.end());
	//	if (names.find(obj_obs) != names.end())
	//	{
	//		throw runtime_error("objective function obs is a constraint, #sad");
	//	}
	//	names.clear();
	//	cnames = constraints.get_nz_obs_names();
	//	names.insert(cnames.begin(),cnames.end());
	//	if (names.find(obj_obs) != names.end())
	//	{
	//		throw runtime_error("objective function obs has non-zero weight and chance constraints are active");
	//	}

	//}

	//else
	//{
	//	if (obj_func_str.size() == 0)
	//	{
	//		f_rec << " warning: no ++opt_objective_function-->forming a generic objective function (1.0 coef for each decision var)" << endl;
	//		for (auto &name : ctl_ord_dec_var_names)
	//			obj_func_coef_map[name] = 1.0;
	//	}

	//	//or if it is a prior info equation
	//	else if (pest_scenario.get_prior_info().find(obj_func_str) != pest_scenario.get_prior_info().end())
	//	{
	//		obj_func_coef_map = pest_scenario.get_prior_info().get_pi_rec_ptr(obj_func_str).get_atom_factors();
	//		//throw_sequentialLP_error("prior-information-based objective function not implemented");
	//	}
	//	else
	//	{
	//		//check if this obj_str is a filename
	//		ifstream if_obj(obj_func_str);
	//		if (!if_obj.good())
	//			throw_sequentialLP_error("unrecognized ++opt_objective_function arg: " + obj_func_str);
	//		else
	//			obj_func_coef_map = pest_utils::read_twocol_ascii_to_map(obj_func_str);
	//	}


	//	//check that all obj_coefs are decsision vars
	//	vector<string> missing_vars;
	//	for (auto &coef : obj_func_coef_map)
	//		if (find(ctl_ord_dec_var_names.begin(), ctl_ord_dec_var_names.end(), coef.first) == ctl_ord_dec_var_names.end())
	//			missing_vars.push_back(coef.first);
	//	if (missing_vars.size() > 0)
	//		throw_sequentialLP_error("the following objective function components are not decision variables: ", missing_vars);

	//}
	optobjfunc.initialize(constraints.get_obs_constraint_names(), dv_names);

	jco.set_base_numeric_pars(current_pars);
	jco.set_base_sim_obs(pest_scenario.get_ctl_observations());
	if (pest_scenario.get_pestpp_options().get_opt_coin_log())
		model.setLogLevel(4 + 8 + 16 + 32);
	initial_report();
	return;
}

void sequentialLP::build_obj_func_coef_array()
{
	ctl_ord_obj_func_coefs = new double[num_dec_vars()];
	double coef;
	map<string, double> obj_func_coef_map = optobjfunc.get_obj_func_coef_map();
	map<string, double>::iterator end = obj_func_coef_map.end();
	int i = 0;
	for (auto &name : dv_names)
	{
		if (obj_func_coef_map.find(name) != end)
			ctl_ord_obj_func_coefs[i] = obj_func_coef_map.at(name);
		else
			ctl_ord_obj_func_coefs[i] = 0.0;
		i++;
	}
	return;
}

void sequentialLP::iter_infeasible_report()
{
	ofstream &f_rec = file_mgr_ptr->rec_ofstream();
	int num_inf = model.numberPrimalInfeasibilities();

	double * inf_ray = new double[num_dec_vars()];
	inf_ray = model.infeasibilityRay(true);
	//inf_ray = model.unboundedRay();

	double sum_inf_primal = model.sumPrimalInfeasibilities();
	double sum_inf_dual = model.sumDualInfeasibilities();
	double inf_cost = model.infeasibilityCost();
	f_rec << "----------------------------------" << endl;
	f_rec << "  ---  infeasibility report  ---  " << endl;
	f_rec << "----------------------------------" << endl;

	f_rec << "-->number primary infeasible constraints: " << num_inf << endl;
	f_rec << "-->sum of primal infeasibilites: " << sum_inf_primal << endl;
	//ClpSimplex::Status stat = model.getRowStatus(0);
	//f_rec << "-->infeasibility cost: " << inf_cost << endl;
	if (model.rayExists())
		for (int i=0;i<num_dec_vars();++i)
			f_rec << inf_ray[i] << endl;
	f_rec << endl << endl;
	return;
}

void sequentialLP::build_dec_var_bounds()
{
	//set the decision var lower and upper bound arrays
	dec_var_lb = new double[num_dec_vars()];
	dec_var_ub = new double[num_dec_vars()];
	Parameters parlbnd = pest_scenario.get_ctl_parameter_info().get_low_bnd(dv_names);
	Parameters parubnd = pest_scenario.get_ctl_parameter_info().get_up_bnd(dv_names);
	for (int i = 0; i < num_dec_vars(); ++i)
	{
		dec_var_lb[i] = parlbnd.get_rec(dv_names[i]) - current_pars.get_rec(dv_names[i]);
		dec_var_ub[i] = parubnd.get_rec(dv_names[i]) - current_pars.get_rec(dv_names[i]);
	}
}

void sequentialLP::iter_solve()
{

	ofstream &f_rec = file_mgr_ptr->rec_ofstream();

	//convert Jacobian_1to1 to CoinPackedMatrix
	cout << "  ---  forming LP model  --- " << endl;
	CoinPackedMatrix matrix = jacobian_to_coinpackedmatrix();

	build_dec_var_bounds();

	pair<vector<double>, vector<double>> bounds = constraints.get_constraint_bound_vectors(current_pars, current_constraints_sim);
	constraints.presolve_report(slp_iter,current_pars, current_constraints_sim);
	//load the linear simplex model
	//model.loadProblem(matrix, dec_var_lb, dec_var_ub, ctl_ord_obj_func_coefs, constraint_lb, constraint_ub);
	model.loadProblem(matrix, dec_var_lb, dec_var_ub, ctl_ord_obj_func_coefs, bounds.first.data(), bounds.second.data());
	vector<string> obs_constraint_names = constraints.get_obs_constraint_names();
	vector<string> pi_constraint_names = constraints.get_pi_constraint_names();
	for (int i = 0; i < constraints.num_obs_constraints(); ++i)
		model.setRowName(i, obs_constraint_names[i]);
	for (int i = 0; i < pi_constraint_names.size(); ++i)
		model.setRowName(i+constraints.num_obs_constraints(), pi_constraint_names[i]);
	for (int i = 0; i < num_dec_vars(); ++i)
		model.setColumnName(i, dv_names[i]);

	model.setOptimizationDirection(pest_scenario.get_pestpp_options().get_opt_direction());
	//if maximum ++opt_coin_loglev, then also write iteration specific mps files
	if (pest_scenario.get_pestpp_options().get_opt_coin_log())
	{
		stringstream ss;
		ss << slp_iter << ".mps";
		string mps_name = file_mgr_ptr->build_filename(ss.str());
		model.writeMps(mps_name.c_str(),0,1);
	}
	f_rec << "  ---  solving linear program for iteration " << slp_iter << "  ---  " << endl;
	cout << "  ---  solving linear program for iteration " << slp_iter << "  ---  " << endl;

	//solve the linear program
	ClpPresolve presolve_info;
	ClpSimplex* presolved_model = presolve_info.presolvedModel(model);

	//if presolvedModel is Null, then it is primal infeasible, so
	//try the dual
	if (!presolved_model)
	{
		f_rec << "  ---  primal presolve model infeasible, crashing solution with additional dual and primal solves..." << endl;
		cout << "  ---  primal presolve model infeasible, crashing solution..." << endl;

		model.moveTowardsPrimalFeasible();
		model.dual(1);
		model.checkSolution();
		model.primal(1);
	}

	//update the status arrays of both the presolve and original models
	presolve_info.postsolve(true);

	//this seems to help with some test problems
	model.dual(1);
	//int all_slack_basis = model.crash(0.001,0);

	//model.primal();
	model.checkSolution();
	if (!model.primalFeasible())
	{
		model.dual(1);
		model.checkSolution();
		model.primal(1);
	}

	//check the solution
	model.checkSolution();
	if (model.isProvenOptimal())
	{
		f_rec << " iteration " << slp_iter << " linear solution is proven optimal" << endl << endl;
		cout << " iteration " << slp_iter << " linear solution is proven optimal" << endl << endl;
	}

	else if (!model.primalFeasible())
	{
		f_rec << "  ---  warning: primal solution infeasible, terminating iterations ---  " << endl;
		cout << "  ---  warning: primal solution infeasible, terminating iterations  ---  " << endl;
		iter_infeasible_report();
		terminate = true;
	}

	else
	{
		f_rec << endl << "iteration " << slp_iter << " linear solution is not proven optimal...continuing" << endl << endl;
		cout << endl << "iteration " << slp_iter << " linear solution is not proven optimal...continuing" << endl << endl;
	}

	f_rec << endl << "  ---  linear program solution complete for iteration " << slp_iter << "  ---  " << endl;
	cout << endl << "  ---  linear program solution complete for iteration " << slp_iter << "  ---  " << endl;

	return;
}

CoinPackedMatrix sequentialLP::jacobian_to_coinpackedmatrix()
{

	Eigen::SparseMatrix<double> eig_ord_jco = jco.get_matrix(constraints.get_obs_constraint_names(), dv_names);

	file_mgr_ptr->rec_ofstream() << "number of nonzero elements in response matrix: " << eig_ord_jco.nonZeros() << " of " << eig_ord_jco.size() << endl;
	cout << "number of nonzero elements in response matrix: " << eig_ord_jco.nonZeros() << " of " << eig_ord_jco.size() << endl;

	if (eig_ord_jco.nonZeros() == 0)
		throw_sequentialLP_error("response matrix all zeros...cannot continue");

	//cout << eig_ord_jco << endl;

	//the triplet elements to pass to the coinpackedmatrix constructor
	int num_elems = eig_ord_jco.nonZeros() + constraints.get_num_nz_pi_constraint_elements();
	int * row_idx = new int[num_elems];
	int * col_idx = new int[num_elems];
	double * elems = new double[num_elems];
	int elem_count = 0;

	//iterate through the eigen sparse matrix
	int elems_par;
	int npar_zelems = 0;
	for (int i = 0; i < eig_ord_jco.outerSize(); ++i)
	{
		elems_par = 0;
		for (Eigen::SparseMatrix<double>::InnerIterator it(eig_ord_jco, i); it; ++it)
		{
			row_idx[elem_count] = it.row();
			col_idx[elem_count] = it.col();
			elems[elem_count] = it.value();
			elems_par++;
			elem_count++;
		}
		if (elems_par == 0)
		{
			//cout << "all zero elements for decision variable: " << ctl_ord_dec_var_names[i] << endl;
			npar_zelems++;
		}
	}
	cout << "number of decision variables with all zero elements: " << npar_zelems << endl;
	if (elem_count != eig_ord_jco.nonZeros())
	{
		throw_sequentialLP_error("sequentialLP::jacobian_to_coinpackedmatrix() error: wrong number of triplet components");
	}

	if (elem_count == 0)
	{
		throw_sequentialLP_error("sequentialLP::jacobian_to_coinpackedmatrix() error: zero triplets found");
	}

	//prior information constraints
	int irow = constraints.num_obs_constraints();
	int jcol;
	vector<string>::iterator start = dv_names.begin();
	vector<string>::iterator end = dv_names.end();
	PriorInformation constraints_pi = constraints.get_pi_constraints();
	for (auto &pi_name : constraints.get_pi_constraint_names())
	{
		for (auto &pi_factor : constraints_pi.get_pi_rec_ptr(pi_name).get_atom_factors())
		{
			jcol = find(start, end, pi_factor.first) - start;
			row_idx[elem_count] = irow;
			col_idx[elem_count] = jcol;
			elems[elem_count] = pi_factor.second;
			elem_count++;
		}
		irow++;
	}
	if (elem_count != num_elems)
		throw_sequentialLP_error("problem packing prior information constraints into CoinPackedMatrix...");

	CoinPackedMatrix matrix(true,row_idx,col_idx,elems,elem_count);

	//this is useful for debugging
	//cout << eig_ord_jco << endl << endl;
	//matrix.dumpMatrix();
	return matrix;
}

void sequentialLP::solve()
{
	ofstream &f_rec = file_mgr_ptr->rec_ofstream();

	slp_iter = 1;
	while (true)
	{
		f_rec << endl << endl << "  ---------------------------------" << endl;
		f_rec <<         "  --- starting LP iteration " << slp_iter << "  ---  " << endl;
		f_rec << "  ---------------------------------" << endl << endl << endl;
		cout << endl << endl << "  ----------------------------------" << endl;
		cout << "  --- starting LP iteration " << slp_iter << "  ---  " << endl;
		cout << "  ---------------------------------" << endl << endl << endl;

		iter_presolve();
		iter_solve();
		iter_postsolve();
		if (terminate) break;
		slp_iter++;
		if (slp_iter > pest_scenario.get_control_info().noptmax)
			break;
	}
	f_rec << endl << "  ---  objective function sequence  ---   " << endl << setw(10) << "iteration" << setw(15) << "obj func" << endl;
	int i = 0;
	for (auto &obj : iter_obj_values)
	{
		f_rec << setw(10) << i << setw(15) << obj << endl;
		i++;
	}
	f_rec << "  ---  best objective function value: " << obj_best << endl;
	cout << "  ---  best objective function value: " << obj_best << endl;


	if (!pest_scenario.get_pestpp_options().get_opt_skip_final())
	{
		f_rec << "  ---  running model one last time with best decision variables  ---  " << endl;
		cout << "  ---  running model one last time with best decision variables  ---  " << endl;
	
		bool success = make_upgrade_run(best_pars, current_constraints_sim);
		if (!success)
		{
			throw_sequentialLP_error("error running model with best decision variable values");
		}
		//write 'best' rei
		of_wr.write_opt_constraint_rei(file_mgr_ptr->open_ofile_ext("res"), slp_iter, current_pars,
			pest_scenario.get_ctl_observations(), current_constraints_sim);
		file_mgr_ptr->close_file("res");
		
		
	}
}


void sequentialLP::iter_postsolve()
{
	ofstream &f_rec = file_mgr_ptr->rec_ofstream();
	f_rec << "  ---  processing results for iteration " << slp_iter << " LP solution  ---  " << endl << endl;
	
	//extract (optimal) decision vars
	//and track some info for convergence checking
	const double *dec_var_vals = model.getColSolution();
	double max_abs_dec_var_change = -1.0E+10;
	double max_abs_dec_var_val = -1.0E+10;

	double diff, val;
	Parameters upgrade_pars(current_pars);
	if (!model.primalFeasible())
	{
		for (auto &name : dv_names)
		{
			upgrade_pars.update_rec(name, numeric_limits<double>::lowest());

		}
		stringstream ss;
		ss << slp_iter << ".par";
		of_wr.write_par(file_mgr_ptr->open_ofile_ext(ss.str()), upgrade_pars, *par_trans.get_offset_ptr(), *par_trans.get_scale_ptr());
		file_mgr_ptr->close_file(ss.str());
		of_wr.write_par(file_mgr_ptr->open_ofile_ext("par"), upgrade_pars, *par_trans.get_offset_ptr(), *par_trans.get_scale_ptr());
		file_mgr_ptr->close_file("par");
		f_rec << " --- warning: parameter file " << ss.str() << "contains double min (extreme) values" << endl;
		f_rec << "     and no res / rei files are being written b/c solution is infeasible" << endl;

		return;
	}

	//check for denormal solution values and values (slightly) out of bounds - reset accordingly...
	string name;
	vector<double> optimal_vals;

	for (int i = 0; i < num_dec_vars(); ++i)
	{
		name = dv_names[i];
		val = dec_var_vals[i];
		
		
		if (val < dec_var_lb[i])
		{
			f_rec << "resetting optimal value for " << name << "from " << val << "to lower bound of" << dec_var_lb[i] << endl;
			optimal_vals.push_back(dec_var_lb[i]);
		}
		else if (val > dec_var_ub[i])
		{
			f_rec << "resetting optimal value for " << name << "from " << val << "to upper bound of" << dec_var_ub[i] << endl;
			optimal_vals.push_back(dec_var_ub[i]);
		}
		else
		{
			if (abs(val) < 1.0e-30) // in case forward model is single precision
			{
				f_rec << "flushing optimal value for " << name << "from " << val << "to 0.0" << endl;
				optimal_vals.push_back(0.0);
			}
			else
				optimal_vals.push_back(val);
		}
	}
	//delete dec_var_vals;

	Parameters dv_changes = upgrade_pars;
	for (int i = 0; i < num_dec_vars(); ++i)
	{
		name = dv_names[i];
		val = current_pars[name];
		diff = abs(optimal_vals[i] - current_pars[name]);
		upgrade_pars.update_rec(name,optimal_vals[i] + val);
		dv_changes.update_rec(name,upgrade_pars[name] - val);
		max_abs_dec_var_change = (diff > max_abs_dec_var_change) ? diff : max_abs_dec_var_change;
		max_abs_dec_var_val = (abs(val) > max_abs_dec_var_val) ? val : max_abs_dec_var_val;
	}
	max_abs_dec_var_change /= max(max_abs_dec_var_val,1.0);

	pair<double, double> cur_new_obj = postsolve_decision_var_report(upgrade_pars);

	

	double obj_val = model.getObjValue();

	f_rec << endl << endl << "  ---  iteration " << slp_iter << " objective function value: " << setw(15) << cur_new_obj.second << "  ---  " << endl << endl;
	cout << endl << endl << "  ---  iteration " << slp_iter << " objective function value: " << setw(15) << cur_new_obj.second << "  ---  " << endl << endl;

	row_price = model.getRowPrice();

	constraints.postsolve_pi_constraints_report(current_pars, upgrade_pars,slp_iter);

	Observations upgrade_obs = current_constraints_sim;
	
	Eigen::VectorXd est_obs_vec = current_constraints_sim.get_data_eigen_vec(constraints.get_obs_constraint_names()) +  
		jco.get_matrix(constraints.get_obs_constraint_names(), dv_names) *
		dv_changes.get_partial_data_eigen_vec(dv_names);
	upgrade_obs.update_without_clear(constraints.get_obs_constraint_names(), est_obs_vec);
	constraints.postsolve_obs_constraints_report(current_constraints_sim, upgrade_obs, "estimated",slp_iter);
	constraints.write_res_files(upgrade_obs, upgrade_pars, "est",slp_iter);
	
	if (!super_secret_option)
	{
		f_rec << "  ---  running the model once with optimal decision var values" << endl;
		bool success = make_upgrade_run(upgrade_pars, upgrade_obs);
		if (!success)
		{
			f_rec << " --- optimal decision var value run failed, cannot continue... " << endl;
			return;
		}
		
		//postsolve_constraint_report(upgrade_obs, upgrade_pars, "simulated");
		constraints.postsolve_obs_constraints_report(current_constraints_sim, upgrade_obs, "simulated", slp_iter);
		constraints.write_res_files(upgrade_obs, upgrade_pars, "sim", slp_iter);
		
	}

	//track the objective function values
	if (slp_iter == 1)
	{
		iter_obj_values.push_back(cur_new_obj.first);
		obj_best = cur_new_obj.second;
	}
	iter_obj_values.push_back(cur_new_obj.second);
	double obj_func_change = abs(cur_new_obj.first - cur_new_obj.second) / abs(max(max(cur_new_obj.first,cur_new_obj.second),1.0));
	//if this is a max problem
	if (pest_scenario.get_pestpp_options().get_opt_direction() == -1)
	{
		if (cur_new_obj.second > obj_best)
		{
			best_pars.update_without_clear(dv_names, upgrade_pars.get_data_vec(dv_names));
			obj_best = cur_new_obj.second;
		}
	}
	else
		if (cur_new_obj.second < obj_best)
		{
			best_pars.update_without_clear(dv_names, upgrade_pars.get_data_vec(dv_names));
			obj_best = cur_new_obj.second;
		}

	//check for changes for in constraints
	double max_abs_constraint_change = constraints.get_max_constraint_change(current_constraints_sim, upgrade_obs);
	double max_abs_constraint_val = -1.0E+10;
	Observations constraints_sim = current_constraints_sim;
	for (auto &name : constraints.get_obs_constraint_names())
	{

		diff = abs(constraints_sim[name] - upgrade_obs[name]);
		max_abs_constraint_change = (diff > max_abs_constraint_change) ? diff : max_abs_constraint_change;
		max_abs_constraint_val = (val > max_abs_constraint_val) ? val : max_abs_constraint_val;
	}
	max_abs_constraint_change /= max(max_abs_constraint_val,1.0);

	map<string,double> invalid_dec_vars = get_out_of_bounds_dec_vars(upgrade_pars);

	if (super_secret_option)
	{
		f_rec << "super secret option active...done" << endl;
		return;
	}

	//convergence check
	double opt_iter_tol = pest_scenario.get_pestpp_options().get_opt_iter_tol();

	f_rec << endl << "  ---  convergence check iteration " << slp_iter << "  ---  " << endl << endl;
	f_rec << "-->                     ++opt_iter_tol:" << setw(15) << opt_iter_tol << endl;
	f_rec << "-->scaled max decision variable change:" << setw(15) << max_abs_dec_var_change << endl;
	f_rec << "-->      scaled  max constraint change:" << setw(15) << max_abs_constraint_change << endl;
	f_rec << "-->   scaled objective function change:" << setw(15) << obj_func_change << endl;

	bool valid = true;
	if (invalid_dec_vars.size() > 0)
	{
		valid = false;
		f_rec << "-->the following decision variables are out of bounds (distance shown):" << endl;
		for (auto &name : invalid_dec_vars)
			f_rec << "-->   " << name.first << setw(15) << name.second << endl;
	}
	map<string, double> invalid_constraints = constraints.get_unsatified_obs_constraints(upgrade_obs, opt_iter_tol);
	map<string, double> invalid_pi_constraints = constraints.get_unsatified_pi_constraints(upgrade_pars, opt_iter_tol);
	invalid_constraints.insert(invalid_pi_constraints.begin(), invalid_pi_constraints.end());
	if (invalid_constraints.size() > 0)
	{
		valid = false;
		f_rec << "-->the following constraints are not satisfied (distance shown):" << endl;
		for (auto &name : invalid_constraints)
			f_rec << "-->   " << name.first << setw(15) << name.second << endl;
	}

	//if three or more iters have passed, start testing the last three
	//obj func vals to see if we have stagnated
	bool obj_func_stag = false;

	int num_obj_vals = iter_obj_values.size();
	if (num_obj_vals > 3)
	{
		obj_func_stag = true;
		double last = iter_obj_values[num_obj_vals-1];
		for (int i = num_obj_vals-2; i >= num_obj_vals-4; --i)
		{
			diff = abs((iter_obj_values[i] - last) / max(1.0, last));
			if (diff > opt_iter_tol)
			{
				obj_func_stag = false;
				break;
			}
			last = iter_obj_values[i];
		}
	}

	//if the obj func has stopped changing or all the other criteria are met
	if ((obj_func_stag) || ((valid) && (max_abs_dec_var_change <= opt_iter_tol) && (max_abs_constraint_change <= opt_iter_tol) && (obj_func_change <= opt_iter_tol)))
	{
		f_rec << endl << "  ---  SLP convergence  ---  " << endl << endl;
		cout << endl << "  ---  SLP convergence  ---  " << endl << endl;
		terminate = true;
	}
	else
	{

	}


	//if continuing, update the master decision var instance
	current_pars.update_without_clear(dv_names, upgrade_pars.get_data_vec(dv_names));
	current_constraints_sim = upgrade_obs;
	return;
}

bool sequentialLP::make_upgrade_run(Parameters &upgrade_pars, Observations &upgrade_obs)
{

	cout << "  ---  running the model once with optimal decision variables  ---  " << endl;
	int run_id = run_mgr_ptr->add_run(par_trans.ctl2model_cp(upgrade_pars));
	run_mgr_ptr->run();
	bool success = run_mgr_ptr->get_run(run_id, upgrade_pars, upgrade_obs);
	if (success)
		par_trans.model2ctl_ip(upgrade_pars);
	return success;
}


void sequentialLP::iter_presolve()
{
	ofstream &f_rec = file_mgr_ptr->rec_ofstream();
	f_rec << "  ---  calculating response matrix for iteration " << slp_iter << "  ---  " << endl;
	cout << "  ---  calculating response matrix for iteration " << slp_iter << "  ---  " << endl;

	//read an existing jacobain
	string basejac_filename = pest_scenario.get_pestpp_options().get_basejac_filename();
	Parameters pars = best_pars;
	if ((slp_iter == 1) && (basejac_filename.size() > 0))
	{
		jco.read(basejac_filename);
		//check to make sure decision vars and constraints are found
		vector<string> temp = jco.get_base_numeric_par_names();
		set<string> names(temp.begin(),temp.end());
		set<string>::iterator start = names.begin();
		set<string>::iterator end = names.end();
		//vector<string>::iterator ext_start = ctl_ord_ext_var_names.begin();
		//vector<string>::iterator ext_end = ctl_ord_ext_var_names.end();
		set<string> ext_set(ext_dv_names.begin(), ext_dv_names.end());
		vector<string> missing;
		for (auto &name : dv_names)
			//if this dec var is not in the jco and is not an external var
			//if ((find(start, end, name) == end) && (find(ext_start,ext_end,name) == ext_end))
			if ((names.find(name) == end) && (ext_set.find(name) == ext_set.end()))
				missing.push_back(name);
		if (missing.size() > 0)
			throw_sequentialLP_error("the following decision vars were not found in the jacobian " + basejac_filename + " : ", missing);
		if ((constraints.get_use_chance()) && (constraints.get_use_fosm()))
		{
			for (auto& name : constraints.get_adj_par_names())
				//if (find(start, end, name) == end)
				if ((names.find(name) == end))
					missing.push_back(name);
			if (missing.size() > 0)
				throw_sequentialLP_error("the following adjustable parameters were not found in the jacobian " + basejac_filename + " : ", missing);
		}


		names.clear();
		temp = jco.get_sim_obs_names();
		names.insert(temp.begin(),temp.end());
		start = names.begin();
		end = names.end();
		for (auto &name : constraints.get_obs_constraint_names())
			if (find(start, end, name) == end)
				missing.push_back(name);
		if (missing.size() > 0)
			throw_sequentialLP_error("the following constraints were not found in the jacobian " + basejac_filename + " : ", missing);

		for (auto &name : constraints.get_nz_obs_names())
			if (find(start, end, name) == end)
				missing.push_back(name);
		if (missing.size() > 0)
			throw_sequentialLP_error("the following non-zero weight observations were not found in the jacobian " + basejac_filename + " : ", missing);
		if ((constraints.get_use_chance()) && (constraints.get_use_fosm()))
			constraints.set_jco(jco);
		string res_filename = pest_scenario.get_pestpp_options().get_hotstart_resfile();
		//Observations constraints_sim = constraints.get_current_constaints_sim();
		if (!res_filename.empty())
		{
			stringstream message;
			message << "  reading  residual file " << res_filename << " for hot-start...";
			cout << message.str();
			f_rec << message.str();
			for (auto &oname : pest_scenario.get_ctl_ordered_obs_names())
				current_constraints_sim[oname] = -1.0e+30;
			pest_utils::read_res(res_filename, current_constraints_sim);
			f_rec << "done" << endl;
			cout << "done" << endl;
			if ((constraints.should_update_chance(slp_iter)) && (!constraints.get_use_fosm()))
				constraints.add_runs(slp_iter, current_pars,current_constraints_sim, run_mgr_ptr);
			

		}
		else
		{
			//make the intial base run
			cout << "  ---  running the model once with initial decision variables  ---  " << endl;
			int run_id = run_mgr_ptr->add_run(par_trans.ctl2model_cp(current_pars));
			//this would be only for stack runs since the fosm runs should have been in the jco
			if ((constraints.should_update_chance(slp_iter)) && (!constraints.get_use_fosm()))
				constraints.add_runs(slp_iter, current_pars, current_constraints_sim, run_mgr_ptr);
			/*else
			{
				cout << "  ---  running the model once with initial decision variables  ---  " << endl;
			}*/
			run_mgr_ptr->run();
			Parameters pars;
			bool success = run_mgr_ptr->get_run(run_id, pars, current_constraints_sim);
			if (!success)
				throw_sequentialLP_error("initial (base) run with initial decision vars failed...cannot continue");

		}
		constraints.process_runs(run_mgr_ptr, slp_iter);
	}

	//otherwise, fill the jacobian
	else
	{
		set<string> out_of_bounds;
		vector<string> names_to_run;
		for (auto& name : dv_names)
		{
			if (find(ext_dv_names.begin(), ext_dv_names.end(), name) == ext_dv_names.end())
				names_to_run.push_back(name);
		}
		if ((constraints.should_update_chance(slp_iter)) && (constraints.get_use_fosm()))
		{
				//names_to_run.insert(names_to_run.end(), adj_par_names.begin(), adj_par_names.end());
				vector<string> fosm_par_names = constraints.get_fosm_par_names();
				names_to_run.insert(names_to_run.end(), fosm_par_names.begin(), fosm_par_names.end());
		}

		//turn down the purb value each iteration
		if ((slp_iter > 1) && (iter_derinc_fac != 1.0))
		{
			cout << "  ---  decreasing derinc for active decision variables " << endl;
			file_mgr_ptr->rec_ofstream() << "  ---  decreasing derinc for active decision variables " << endl;

			//find group names for active dec vars
			vector<string> act_dv_grps;
			ParameterGroupInfo* pinfo = pest_scenario.get_base_group_info_ptr();

			const ParameterGroupRec *gr_ptr;
			for (auto &name : dv_names)
			{
				gr_ptr = pinfo->get_group_rec_ptr(name);
				if (find(act_dv_grps.begin(), act_dv_grps.end(), gr_ptr->name) == act_dv_grps.end())
					act_dv_grps.push_back(gr_ptr->name);
			}
			double org_derinc, new_derinc;
			ofstream &frec = file_mgr_ptr->rec_ofstream();
			frec << setw(20) << left << "group_name" << setw(20) << left << "old_derinc" << setw(20) << left << "new_derinc" << endl;
			for (auto &name : act_dv_grps)
			{
				org_derinc = pinfo->get_group_by_groupname_4_mod(name)->derinc;
				new_derinc = org_derinc * iter_derinc_fac;
				pinfo->get_group_by_groupname_4_mod(name)->derinc = new_derinc;
				frec << setw(20) << left << name << setw(20) << left << org_derinc << setw(20) << left << new_derinc << endl;
			}
		}
		
		bool init_obs = false;
		if (slp_iter == 1) init_obs = true;
		
		bool success = jco.build_runs(current_pars, current_constraints_sim, names_to_run, par_trans,
			pest_scenario.get_base_group_info(), pest_scenario.get_ctl_parameter_info(),
			*run_mgr_ptr, out_of_bounds,false,init_obs);
		if (!success)
		{
			const set<string> failed = jco.get_failed_parameter_names();
			throw_sequentialLP_error("failed to calc derviatives for the following decision vars: ", failed);
		}

		if ((constraints.should_update_chance(slp_iter)) && (!constraints.get_use_fosm()))
		{
			constraints.add_runs(slp_iter, current_pars, current_constraints_sim, run_mgr_ptr);
		}


		//jco.make_runs(*run_mgr_ptr);
		run_mgr_ptr->run();
		set<int> failed = run_mgr_ptr->get_failed_run_ids();

		//process the remaining responses
		success = jco.process_runs(par_trans, pest_scenario.get_base_group_info(), *run_mgr_ptr, 
			*null_prior, false,false);
		if (!success)
			throw_sequentialLP_error("error processing response matrix runs ", jco.get_failed_parameter_names());
		constraints.process_runs(run_mgr_ptr, slp_iter);
		stringstream ss;
		ss << slp_iter << ".jcb";
		string rspmat_file = file_mgr_ptr->build_filename(ss.str());
		f_rec << endl << "saving iteration " << slp_iter << " reponse matrix to file: " << rspmat_file << endl;
		jco.save(ss.str());

		//check for failed runs
		//TODO: something better than just dying
		if (failed.size() > 0)
			throw_sequentialLP_error("failed runs when filling decision var response matrix...cannot continue ");


		//get the base run and update simulated constraint values
		
		if (init_obs)
		{
			//Observations temp_obs;
			run_mgr_ptr->get_run(0, current_pars, current_constraints_sim, false);
		}
	}

	if (optobjfunc.get_use_obs_obj())
	{
		optobjfunc.update_coef_map_from_jacobian(jco);
		f_rec << "  ---  objective function coefficients for iteration " << slp_iter << "  ---  " << endl;
		optobjfunc.report();
		if (slp_iter == 1)
		{
			obj_init = optobjfunc.get_obj_func_value(current_pars,current_constraints_sim);
			f_rec << endl << "  ---  objective function value (using initial dec var values): " << obj_init << endl << endl;
			cout << endl << "  ---  objective function value (using initial dec var values): " << obj_init << endl << endl;
		}
	}
	
	//build the objective function
	build_obj_func_coef_array();

	stringstream ss;
	if ((constraints.get_use_chance()) && (constraints.get_use_fosm()))
		constraints.set_jco(jco);
	constraints.update_chance_offsets();
	constraints.presolve_chance_report(slp_iter,current_constraints_sim);
	constraints.write_res_files(current_constraints_sim, pars,"jcb",slp_iter);
	return;
}



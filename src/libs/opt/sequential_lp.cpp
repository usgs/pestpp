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
	Covariance &_parcov, FileManager* _file_mgr, OutputFileWriter _of_wr) : pest_scenario(_pest_scenario), run_mgr_ptr(_run_mgr_ptr),
	parcov(_parcov), file_mgr_ptr(_file_mgr),jco(*_file_mgr,_of_wr), of_wr(_of_wr)
{
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
	delete[] constraint_lb;
	delete[] constraint_ub;
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

vector<double> sequentialLP::get_constraint_residual_vec()
{
	if (use_chance)
		return get_constraint_residual_vec(constraints_fosm);
	else
		return get_constraint_residual_vec(constraints_sim);
}

vector<double> sequentialLP::get_constraint_residual_vec(Observations &sim_vals)
{
	vector<double> residuals_vec;
	residuals_vec.resize(num_constraints(), 0.0);

	Observations::const_iterator found_obs;
	Observations::const_iterator not_found_obs = sim_vals.end();
	PriorInformation::const_iterator found_prior_info;

	int i = 0;
	for (vector<string>::iterator b = ctl_ord_obs_constraint_names.begin(), e = ctl_ord_obs_constraint_names.end(); b != e; ++b, ++i)
	{
		found_obs = sim_vals.find(*b);
		//double fosm_offset = post_constraint_offset[*b];
		if (found_obs != not_found_obs)
		{
			residuals_vec[i] = constraints_obs.get_rec(*b) - ((*found_obs).second);
		}

	}
	return residuals_vec;
}

void sequentialLP::initial_report()
{
	ofstream &f_rec = file_mgr_ptr->rec_ofstream();
	f_rec << endl << "  -------------------------------------------------------------" << endl;
	f_rec << "  ---  sequential linear programming problem information  ---  " << endl;
	f_rec << "  -------------------------------------------------------------" << endl << endl;

	f_rec << "-->number of iterations of sequential linear programming (noptmax): " << pest_scenario.get_control_info().noptmax << endl;

	f_rec << "-->objective function sense (direction): " << obj_sense << endl;


	f_rec << "-->number of decision variable: " << num_dec_vars() << endl;
	f_rec << "-->number of observation constraints: " << num_obs_constraints() << endl;
	f_rec << "-->number of prior information constraints: " << num_pi_constraints() << endl;
	if (iter_derinc_fac != 1.0)
		f_rec << "-->iteration DERINC reduction factor (++opt_iter_derinc_fac): " << iter_derinc_fac << endl;
	if (use_chance)
	{
		f_rec << "-->using FOSM-based chance constraints - good choice!" << endl;
		f_rec << "-->++opt_risk and corresponding probit function value: " << setw(10) << risk << setw(20) << probit_val << endl;
		f_rec << "-->number of adjustable parameters for FOSM calcs: " << num_adj_pars() << endl;
		f_rec << "-->number of non-zero weight observations for FOSM calcs: " << num_nz_obs() << endl;
		f_rec << "-->repeat FOSM calcs every: " << pest_scenario.get_pestpp_options().get_opt_recalc_fosm_every() << " iterations" << endl << endl;
	}
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
	map<string, double>::iterator end = obj_func_coef_map.end();
	vector<string> missing;
	if (obj_func_coef_map.size() == 0)
		f_rec << "objective function coefficients defined by observation: " << pest_scenario.get_pestpp_options().get_opt_obj_func() << endl;
	else
	{
		f_rec << "  ---  objective function coefficients  ---  " << endl;
		obj_init = obj_func_report();
		f_rec << endl << "  ---  objective function value (using initial dec var values): " << obj_init << endl << endl;
		cout << endl << "  ---  objective function value (using initial dec var values): " << obj_init << endl << endl;

	}
	f_rec << " note: bound and initial value info reported in 'parameter data' section" << endl << endl;


	f_rec << endl << "  ---  observation constraints in SLP  ---  " << endl;
	f_rec << setw(20) << "name" << setw(20) << "sense" << setw(20) << "value" << endl;
	for (auto &name : ctl_ord_obs_constraint_names)
	{
		f_rec << setw(20) << left << name;
		f_rec << setw(20) << constraint_sense_name[name];
		f_rec << setw(20) << constraints_obs.get_rec(name) << endl;
	}

	if (num_pi_constraints() > 0)
	{
		f_rec << endl << "  ---  prior information constraints in SLP  ---  " << endl;
		f_rec << setw(20) << "name" << setw(20) << "sense" << setw(20) << "value" << endl;
		for (auto &name : ctl_ord_pi_constraint_names)
		{
			f_rec << setw(20) << left << name;
			f_rec << setw(20) << constraint_sense_name[name];
			f_rec << setw(20) << constraints_pi.get_pi_rec_ptr(name).get_obs_value() << endl;
		}
	}

	if (use_chance)
	{
		f_rec << endl << endl << "  ---  chance constraint FOSM information  ---  " << endl;
		f_rec << "   adjustable parameters used in FOSM calculations:" << endl;
		int i = 1;
		for (auto &name : adj_par_names)
		{
			f_rec << setw(15) << name;
			if (i % 6 == 0)
				f_rec << endl;
			i++;

		}
		f_rec << endl;
		if (num_nz_obs() == 0)
		{
			f_rec << endl << endl << "  ---  Note: No nonzero weight observations found." << endl;
			f_rec << "           Prior constraint uncertainty will be used in chance constraint calculations" << endl;
		}
		else
		{
			f_rec << "  non-zero weight observations used for conditioning in FOSM calculations: " << endl;
			int i = 0;
			for (auto &name : nz_obs_names)
			{
				f_rec << setw(15) << name;
				if (i % 6 == 0)
					f_rec << endl;
				i++;
			}
			f_rec << endl;

		}
	}
	return;
}


double sequentialLP::obj_func_report()
{
	vector<string> missing;
	ofstream &f_rec = file_mgr_ptr->rec_ofstream();
	map<string, double>::iterator end = obj_func_coef_map.end();
	f_rec << setw(20) << left << "name" << setw(25) << "obj func coefficient" << endl;
	for (auto &name : ctl_ord_dec_var_names)
	{
		f_rec << setw(20) << left << name;
		if (obj_func_coef_map.find(name) != end)
		{
			f_rec << setw(25) << obj_func_coef_map.at(name) << endl;
		}
		else
		{
			f_rec << setw(25) << "not listed" << endl;
			missing.push_back(name);
		}
	}
	double obj = 0.0;
	for (auto &p : ctl_ord_dec_var_names)
	{
		obj += all_pars_and_dec_vars_initial[p] * obj_func_coef_map[p];
	}
	if (missing.size() > 0)
	{
		f_rec << endl << endl << "WARNING: the following decision variables have '0.0' objective function coef:" << endl;
		cout << endl << endl << "WARNING: the following decision variables have '0.0' objective function coef:" << endl;

		for (auto &name : missing)
		{
			f_rec << "    " << name << endl;
			f_rec << "    " << name << endl;
		}
	}
	return obj;

}

void sequentialLP::presolve_fosm_report()
{
	ofstream &f_rec = file_mgr_ptr->rec_ofstream();
	vector<double> residuals = get_constraint_residual_vec();
	f_rec << endl << "  FOSM-based chance constraint information at start of iteration " << slp_iter << endl;
	f_rec << setw(20) << left << "name" << right << setw(10) << "sense" << setw(15) << "required" << setw(15) << "sim value";
	f_rec << setw(15) << "prior stdev" << setw(15) << "post stdev" << setw(15) << "offset";
	f_rec << setw(15) << "new sim value" << endl;
	vector<string> out_of_bounds;
	for (int i = 0; i<num_obs_constraints(); ++i)
	{
		string name = ctl_ord_obs_constraint_names[i];
		f_rec << setw(20) << left << name;
		f_rec << setw(10) << right << constraint_sense_name[name];
		f_rec << setw(15) << constraints_obs[name];
		f_rec << setw(15) << constraints_sim_initial.get_rec(name);
		f_rec << setw(15) << prior_constraint_stdev[name];
		f_rec << setw(15) << post_constraint_stdev[name];
		f_rec << setw(15) << post_constraint_offset[name];
		f_rec << setw(15) << constraints_fosm[name] << endl;
	}
	f_rec << "  note: 'offset' is the value added to the simulated constraint value to account" << endl;
	f_rec << "        for the uncertainty in the constraint value arsing from uncertainty in the " << endl;
	f_rec << "        adjustable parameters identified in the control file." << endl << endl;

	return;
}

void sequentialLP::presolve_constraint_report()
{
	ofstream &f_rec = file_mgr_ptr->rec_ofstream();
	vector<double> residuals = get_constraint_residual_vec();
	f_rec << endl << "  observation constraint information at start of iteration " << slp_iter << endl;
	f_rec << setw(20) << left << "name" << right << setw(15) << "sense" << setw(15) << "required" << setw(15) << "sim value";
	f_rec << setw(15) << "residual" << setw(15) << "lower bound" << setw(15) << "upper bound" << endl;

	for (int i=0;i<num_obs_constraints();++i)
	{
		string name = ctl_ord_obs_constraint_names[i];
		f_rec << setw(20) << left << name;
		f_rec << setw(15) << right << constraint_sense_name[name];
		f_rec << setw(15) << constraints_obs.get_rec(name);
		if (use_chance)
			f_rec << setw(15) << constraints_fosm.get_rec(name);
		else
			f_rec << setw(15) << constraints_sim.get_rec(name);
		f_rec << setw(15) << residuals[i];
		f_rec << setw(15) << constraint_lb[i];
		f_rec << setw(15) << constraint_ub[i] << endl;

	}

	if (num_pi_constraints() == 0) return;

	//report prior information constraints
	f_rec << endl << "  prior information constraint information at start of iteration " << slp_iter << endl;
	f_rec << setw(20) << left << "name" << right << setw(15) << "sense" << setw(15) << "required" << setw(15) << "sim value";
	f_rec << setw(15) << "residual" << setw(15) << "lower bound" << setw(15) << "upper bound" << endl;
	for (int i = 0; i<num_pi_constraints(); ++i)
	{
		string name = ctl_ord_pi_constraint_names[i];
		PriorInformationRec pi_rec = constraints_pi.get_pi_rec_ptr(name);
		f_rec << setw(20) << left << name;
		f_rec << setw(15) << right << constraint_sense_name[name];
		f_rec << setw(15) << pi_rec.get_obs_value();
		f_rec << setw(15) << pi_rec.calc_residual_and_sim_val(all_pars_and_dec_vars).first;
		f_rec << setw(15) << pi_rec.calc_residual(all_pars_and_dec_vars);
		f_rec << setw(15) << constraint_lb[num_obs_constraints() + i];
		f_rec << setw(15) << constraint_ub[num_obs_constraints() + i] << endl;

	}
	return;
}

void sequentialLP::postsolve_constraint_report(Observations &upgrade_obs,Parameters &upgrade_pars)
{
	ofstream &f_rec = file_mgr_ptr->rec_ofstream();
	f_rec << endl << endl << "     constraint information at end of SLP iteration " << slp_iter << endl << endl;
	f_rec << setw(20) << left << "name" << right << setw(15) << "sense" << setw(15) << "required" << setw(25) << "simplex status";
	f_rec << setw(15) << "price";
	if (use_chance)
		f_rec << setw(15) << "fosm offset";
	f_rec << setw(15) << "current" << setw(15) << "residual";
	f_rec << setw(15) << "new" << setw(15) << "residual" << endl;
	vector<double> cur_residuals = get_constraint_residual_vec();
	vector<double> new_residuals = get_constraint_residual_vec(upgrade_obs);
	double sim_val;
	for (int i = 0; i<num_obs_constraints(); ++i)
	{
		string name = ctl_ord_obs_constraint_names[i];
		sim_val = upgrade_obs[name];
		f_rec << setw(20) << left << name;
		f_rec << setw(15) << right << constraint_sense_name[name];
		f_rec << setw(15) << constraints_obs.get_rec(name);
		f_rec << setw(25) << status_name_map[model.getRowStatus(i)];
		f_rec << setw(15) << row_price[i];
		if (use_chance)
		{
			double offset = post_constraint_offset[name];
			sim_val += offset;
			f_rec << setw(15) << offset;
			f_rec << setw(15) << constraints_sim.get_rec(name) + offset;
			f_rec << setw(15) << cur_residuals[i];
			f_rec << setw(15) << sim_val;
			f_rec << setw(15) << new_residuals[i] - offset << endl;
		}
		else
		{
			f_rec << setw(15) << constraints_sim.get_rec(name);
			f_rec << setw(15) << cur_residuals[i];
			f_rec << setw(15) << sim_val;
			f_rec << setw(15) << new_residuals[i] << endl;

		}
	}

	//report prior information constraints
	if (num_pi_constraints() > 0)
	{
		f_rec << endl << endl << "     prior information constraint information at end of SLP iteration " << slp_iter << endl << endl;
		f_rec << setw(20) << left << "name" << right << setw(15) << "sense" << setw(15) << "required" << setw(25) << "simplex status";
		f_rec << setw(15) << "price" << setw(15) << "current" << setw(15) << "residual";
		f_rec << setw(15) << "new" << setw(15) << "residual" << endl;
		int i = 0;
		for (auto &name : ctl_ord_pi_constraint_names)
		{
			PriorInformationRec pi_rec = constraints_pi.get_pi_rec_ptr(name);
			pair<double,double> cur_sim_resid = pi_rec.calc_residual_and_sim_val(all_pars_and_dec_vars);
			pair<double,double> new_sim_resid = pi_rec.calc_residual_and_sim_val(upgrade_pars);
			f_rec << setw(20) << left << name;
			f_rec << setw(15) << right << constraint_sense_name[name];
			f_rec << setw(15) << pi_rec.get_obs_value();
			f_rec << setw(25) << status_name_map[model.getRowStatus(i+num_obs_constraints())];
			f_rec << setw(15) << row_price[i + num_obs_constraints()];
			f_rec << setw(15) << cur_sim_resid.first;
			f_rec << setw(15) << cur_sim_resid.second;
			f_rec << setw(15) << new_sim_resid.first;
			f_rec << setw(15) << new_sim_resid.second << endl;
			i++;
		}
	}
	stringstream ss;
	ss << slp_iter << ".rei";
	of_wr.write_opt_constraint_rei(file_mgr_ptr->open_ofile_ext(ss.str()), slp_iter, upgrade_pars,
		pest_scenario.get_ctl_observations(), upgrade_obs);
	file_mgr_ptr->close_file(ss.str());
	if (use_chance)
	{
		f_rec << "  ---  note: residual file " << ss.str() << " reports the simulated" << endl;
		f_rec << "       constraint values from the model outputs without FOSM offsets" << endl << endl;
	}
	return;
}

pair<map<string,double>, map<string,double>> sequentialLP::postsolve_check(Observations &upgrade_obs, Parameters &upgrade_pars)
{
	map<string,double> invalid_constraints;
	map<string,double> invalid_dec_vars;

	double sim_val, obs_val, scaled_diff;
	double opt_tol = pest_scenario.get_pestpp_options().get_opt_iter_tol();
	for (int i = 0; i<num_obs_constraints(); ++i)
	{
		string name = ctl_ord_obs_constraint_names[i];
		sim_val = upgrade_obs[name];
		if (use_chance)
		{
			double offset = post_constraint_offset[name];
			sim_val += offset;
		}
		obs_val = constraints_obs[name];
		scaled_diff = abs((obs_val - sim_val) / obs_val);
		//check for invalid obs constraints (e.g. satified)
		if ((constraint_sense_map[name] == ConstraintSense::less_than) && (sim_val > obs_val) && (scaled_diff > opt_tol))
			invalid_constraints[name] = sim_val - obs_val;
		else if ((constraint_sense_map[name] == ConstraintSense::greater_than) && (sim_val < obs_val) && (scaled_diff > opt_tol))
			invalid_constraints[name] = obs_val - sim_val;
		else if ((constraint_sense_map[name] == ConstraintSense::equal_to) && (sim_val != obs_val) && (scaled_diff > opt_tol))
			invalid_constraints[name] = abs(sim_val - obs_val);
	}

	//report prior information constraints
	if (num_pi_constraints() > 0)
	{
		for (auto &name : ctl_ord_pi_constraint_names)
		{
			PriorInformationRec pi_rec = constraints_pi.get_pi_rec_ptr(name);
			pair<double, double> cur_sim_resid = pi_rec.calc_residual_and_sim_val(all_pars_and_dec_vars);
			pair<double, double> new_sim_resid = pi_rec.calc_residual_and_sim_val(upgrade_pars);
			//check for invalid pi constraints
			sim_val = new_sim_resid.first;
			obs_val = pi_rec.get_obs_value();
			scaled_diff = abs((obs_val - sim_val) / obs_val);
			if ((constraint_sense_map[name] == ConstraintSense::less_than) && (sim_val > obs_val) && (scaled_diff > opt_tol))
				invalid_constraints[name] = sim_val - obs_val;
			else if ((constraint_sense_map[name] == ConstraintSense::greater_than) && (sim_val < obs_val) && (scaled_diff > opt_tol))
				invalid_constraints[name] = obs_val - sim_val;
			else if ((constraint_sense_map[name] == ConstraintSense::equal_to) && (sim_val != obs_val) && (scaled_diff > opt_tol))
				invalid_constraints[name] = abs(sim_val - obs_val);
		}
	}

	Parameters ubnd = pest_scenario.get_ctl_parameter_info().get_up_bnd(ctl_ord_dec_var_names);
	Parameters lbnd = pest_scenario.get_ctl_parameter_info().get_low_bnd(ctl_ord_dec_var_names);
	for (auto &name : ctl_ord_dec_var_names)
	{
		sim_val = upgrade_pars[name];
		if (sim_val > ubnd[name])
			invalid_dec_vars[name] = sim_val - ubnd[name];
		else if (sim_val < lbnd[name])
			invalid_dec_vars[name] = lbnd[name] - sim_val;
	}

	return pair<map<string,double>,map<string,double>>(invalid_dec_vars,invalid_constraints);
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
		name = ctl_ord_dec_var_names[i];
		status = model.getColumnStatus(i);
		obj_coef = ctl_ord_obj_func_coefs[i];
		init_val = all_pars_and_dec_vars_initial[name];
		cur_val = all_pars_and_dec_vars[name];

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


double  ErfInv2(double x) 
{
	float tt1, tt2, lnx, sgn;
	sgn = (x < 0) ? -1.0f : 1.0f;

	x = (1 - x)*(1 + x);        // x = 1 - x*x;
	lnx = logf(x);
	//double PI = 3.14159265358979323846;
	tt1 = 2 / (M_PI*0.147) + 0.5f * lnx;
	tt2 = 1 / (0.147) * lnx;

	return(sgn*sqrtf(-tt1 + sqrtf(tt1*tt1 - tt2)));
}

double sequentialLP::get_probit()
{

	//double probit = 0.0;
	//vector<double>::iterator start = probit_inputs.begin();
	//vector<double>::iterator end = probit_inputs.end();

	double output = sqrt(2.0) * ErfInv2((2.0 * risk) - 1.0);

	//find the nearest values in the probit_input vector
	/*auto const it = std::lower_bound(start, end, risk);
	int idx = find(start, end, *it) - start;
	if (idx >= probit_inputs.size())
		throw_sequentialLP_error("error looking up probit function value");
	double output = probit_outputs[idx];
	*/
	return output;


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

	all_pars_and_dec_vars = pest_scenario.get_ctl_parameters();
	all_pars_and_dec_vars_initial = Parameters(all_pars_and_dec_vars);
	all_pars_and_dec_vars_best = Parameters(all_pars_and_dec_vars);
	par_trans = pest_scenario.get_base_par_tran_seq();
	//set ordered dec var name vec
	//and check for illegal parameter transformations
	vector<string> dec_var_groups = pest_scenario.get_pestpp_options().get_opt_dec_var_groups();
	vector<string> ext_var_groups = pest_scenario.get_pestpp_options().get_opt_ext_var_groups();
	dec_var_groups.insert(dec_var_groups.begin(), ext_var_groups.begin(), ext_var_groups.end());
	ctl_ord_dec_var_names.clear();
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
				ctl_ord_dec_var_names.push_back(par_name);
				//check if this is an ext var
				if (find(ext_var_groups.begin(), ext_var_groups.end(), group) != ext_var_groups.end())
					ctl_ord_ext_var_names.push_back(par_name);
			}
		}

		if (num_dec_vars() == 0)
			throw_sequentialLP_error("no decision variables found in groups: ", dec_var_groups);
	}
	//if not ++opt_dec_var_names was passed, use all parameter as decision variables
	else
		ctl_ord_dec_var_names = pest_scenario.get_ctl_ordered_par_names();

	//if any decision vars have a transformation that is not allowed
	vector<string> problem_trans;
	for (auto &name : ctl_ord_dec_var_names)
		if (pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(name)->tranform_type != ParameterRec::TRAN_TYPE::NONE)
			problem_trans.push_back(name);
	if (problem_trans.size() > 0)
		throw_sequentialLP_error("the following decision variables don't have 'none' type parameter transformation: ", problem_trans);

	if (pest_scenario.get_pestpp_options().get_opt_include_bnd_pi())
	{
		PriorInformation *pi_ptr = pest_scenario.get_prior_info_ptr();
		ParameterInfo par_info = pest_scenario.get_ctl_parameter_info();
		stringstream ss;
		string s;
		for (auto &dname : ctl_ord_dec_var_names)
		{
			ss.str("");
			ss << "_lb_" << dname << " 1.0 * " << dname << " = " << par_info.get_parameter_rec_ptr(dname)->lbnd << " 1.0 greater_pi";
			pi_ptr->AddRecord(pest_utils::upper_cp(ss.str()));
			ss.str("");
			ss << "_LB_" << pest_utils::upper_cp(dname);
			pest_scenario.get_ctl_ordered_pi_names_ptr()->push_back(ss.str());
			ss.str("");
			ss << "_ub_" << dname << " 1.0 * " << dname << " = " << par_info.get_parameter_rec_ptr(dname)->ubnd << " 1.0 less_pi";
			pi_ptr->AddRecord(pest_utils::upper_cp(ss.str()));
			ss.str("");
			ss << "_UB_" << pest_utils::upper_cp(dname);
			pest_scenario.get_ctl_ordered_pi_names_ptr()->push_back(ss.str());
			ss.str("");

		}
		pest_scenario.get_ctl_ordered_obs_group_names_ptr()->push_back("GREATER_PI");
		pest_scenario.get_ctl_ordered_obs_group_names_ptr()->push_back("LESS_PI");
	}

	//make sure all decision variables have an initial value of zero
	//jwhite 10june2017 - don't need this.  GWM accounts for base pumping in streamflow constraints
	//but here we leave it up to the user to make sure stream flow depletion constraints are
	//in sync with the initial pumping rates
	//vector<string> problem_vars;
	//for (auto &name : ctl_ord_dec_var_names)
	//	if (pest_scenario.get_ctl_parameters().get_rec(name) != 0.0)
	//		problem_vars.push_back(name);
	//if (problem_vars.size() > 0)
	//	throw_sequentialLP_error("the following decision variables have non-zero initial values", problem_vars);

	//set the decision var lower and upper bound arrays
	/*dec_var_lb = new double[num_dec_vars()];
	dec_var_ub = new double[num_dec_vars()];
	Parameters parlbnd = pest_scenario.get_ctl_parameter_info().get_low_bnd(ctl_ord_dec_var_names);
	Parameters parubnd = pest_scenario.get_ctl_parameter_info().get_up_bnd(ctl_ord_dec_var_names);
	for (int i = 0; i < num_dec_vars(); ++i)
	{
		dec_var_lb[i] = parlbnd.get_rec(ctl_ord_dec_var_names[i]);
		dec_var_ub[i] = parubnd.get_rec(ctl_ord_dec_var_names[i]);
	}*/



	//---------------------------
	//  ---  obs constraints  ---
	//---------------------------
	//set the two constraints attribs and ordered constraint name vec
	vector<string> constraint_groups = pest_scenario.get_pestpp_options().get_opt_constraint_groups();
	//if not constraint groups set as ++ arg, then look for all compatible obs groups
	if (constraint_groups.size() == 0)
	{
		for (auto &group : pest_scenario.get_ctl_ordered_obs_group_names())
		{
			pair<ConstraintSense, string> sense = get_sense_from_group_name(group);
			if (sense.first != ConstraintSense::undefined)
			{
				constraint_groups.push_back(group);
			}
		}
	}

	//if we still don't have any constraint groups, something is wrong
	if (constraint_groups.size() == 0)
		throw_sequentialLP_error("no viable observation or prior information constraint groups found.  Constraint group names must start with the following: {'l_','less','g_','greater'}");

	ctl_ord_obs_constraint_names.clear();
	//if the ++opt_constraint_groups arg was passed
	if (constraint_groups.size() != 0)
	{
		//first make sure all the groups are actually listed in the control file
		vector<string> missing;
		vector<string> pst_groups = pest_scenario.get_ctl_ordered_obs_group_names();
		vector<string>::iterator end = pst_groups.end();
		vector<string>::iterator start = pst_groups.begin();
		for (auto grp : constraint_groups)
			if (find(start, end, grp) == end)
				missing.push_back(grp);
		if (missing.size() > 0)
			throw_sequentialLP_error("the following ++opt_constraint_groups were not found: ", missing);

		//find the observations in constraints groups
		ObservationInfo oinfo = pest_scenario.get_ctl_observation_info();
		string group;
		double weight;
		end = constraint_groups.end();
		start = constraint_groups.begin();
		int nobszw = 0;
		for (auto &obs_name : pest_scenario.get_ctl_ordered_obs_names())
		{
			group = oinfo.get_observation_rec_ptr(obs_name)->group;
			weight = oinfo.get_weight(obs_name);
			if (find(start, end, group) != end)
			{
				if (weight == 0.0)
				{
					//cout << "Warning: observation constraint " << obs_name << " has 0.0 weight, skipping" << endl;
					nobszw++;
					f_rec << "Warning: observation constraint " << obs_name << " has 0.0 weight, skipping" << endl;
				}
				else
				{
					ctl_ord_obs_constraint_names.push_back(obs_name);
				}

			}
		}
		cout << "Warning: " << nobszw << " of the observation constraints (see rec file for list) have 0.0 weight, skipping" << endl;


		//look for prior information constraints
		const PriorInformation* pinfo = pest_scenario.get_prior_info_ptr();
		PriorInformationRec pi_rec;
		for (auto &pi_name : pest_scenario.get_ctl_ordered_pi_names())
		{
			group = pinfo->get_pi_rec_ptr(pi_name).get_group();
			if (find(start, end, group) != end)
			{
				ctl_ord_pi_constraint_names.push_back(pi_name);
				pi_rec = pinfo->get_pi_rec_ptr(pi_name);
				//pi_constraint_factors[pi_name] = pi_rec.get_atom_factors();
				//pi_constraint_rhs[pi_name] = pi_rec.get_obs_value();
				constraints_pi.AddRecord(pi_name, &pi_rec);

			}
		}

		//check the pi constraint factors for compatibility with available
		start = ctl_ord_dec_var_names.begin();
		end = ctl_ord_dec_var_names.end();
		map<string, vector<string>> missing_map;
		if (num_pi_constraints() > 0)
		{
			for (auto &pi_name : ctl_ord_pi_constraint_names)
			{
				pi_rec = constraints_pi.get_pi_rec_ptr(pi_name);
				missing.clear();
				for (auto &pi_factor : pi_rec.get_atom_factors())
					if (find(start, end, pi_factor.first) == end)
						missing.push_back(pi_factor.first);
				if (missing.size() > 0)
					missing_map[pi_name] = missing;
			}
		}
		if (missing_map.size() > 0)
		{
			stringstream ss;
			ss << " the following prior information constraints reference parameters that are not treated as decision variables:" << endl;
			for (auto &missing_pi : missing_map)
			{
				ss << missing_pi.first << ": ";
				for (auto &par_name : missing_pi.second)
					ss << par_name << ",";
			}
			throw_sequentialLP_error("errors in prior information constraints:" + ss.str());
		}

		//TODO: investigate a pi constraint only formulation
		if (num_obs_constraints() == 0)
			throw_sequentialLP_error("no constraints found in groups: ", constraint_groups);
	}


	//if not ++opt_constraint_names was passed, use all observations and prior information as constraints
	/*else
	{
		ctl_ord_obs_constraint_names = pest_scenario.get_ctl_ordered_obs_names();
		ctl_ord_pi_constraint_names = pest_scenario.get_ctl_ordered_pi_names();
	}*/


	constraints_obs = pest_scenario.get_ctl_observations().get_subset(ctl_ord_obs_constraint_names.begin(), ctl_ord_obs_constraint_names.end());
	constraints_sim = Observations(constraints_obs);

	//build map of obs constraint sense
	vector<string> problem_constraints;
	for (auto &name : ctl_ord_obs_constraint_names)
	{
		string group = pest_scenario.get_ctl_observation_info().get_observation_rec_ptr(name)->group;
		pair<ConstraintSense,string> sense = get_sense_from_group_name(group);
		if (sense.first == ConstraintSense::undefined)
		{
			problem_constraints.push_back(name + ',' + group);
		}
		constraint_sense_map[name] = sense.first;
		constraint_sense_name[name] = sense.second;

	}
	if (problem_constraints.size() > 0)
	{
		throw_sequentialLP_error("the following obs constraints do not have a correct group name prefix {'l_','less','g_','greater','e_','equal'}: ", problem_constraints);
	}

	//build map of pi constraint sense
	problem_constraints.clear();
	for (auto &name : ctl_ord_pi_constraint_names)
	{
		string group = pest_scenario.get_prior_info_ptr()->get_pi_rec_ptr(name).get_group();
		pair<ConstraintSense, string> sense = get_sense_from_group_name(group);
		if (sense.first == ConstraintSense::undefined)
		{
			problem_constraints.push_back(name + ',' + group);
		}
		constraint_sense_map[name] = sense.first;
		constraint_sense_name[name] = sense.second;
	}
	if (problem_constraints.size() > 0)
	{
		throw_sequentialLP_error("the following prior info constraints do not have a correct group name prefix {'l_','less','g_','greater','e_','equal'}: ", problem_constraints);
	}

	//allocate the constraint bound arrays
	constraint_lb = new double[num_constraints()];
	constraint_ub = new double[num_constraints()];










	//------------------------------------------
	//  ---  chance constratints and fosm  ---
	//------------------------------------------
	risk = pest_scenario.get_pestpp_options().get_opt_risk();
	if (risk != 0.5)
	{
		use_chance = true;
		std_weights = pest_scenario.get_pestpp_options().get_opt_std_weights();

		//make sure risk value is valid
		if ((risk > 1.0) || (risk < 0.0))
			throw_sequentialLP_error("++opt_risk parameter must between 0.0 and 1.0");

		//reset risk extreme risk values
		if (risk > 0.999)
		{
			f_rec << endl << "  ---  note: resetting risk value of " << risk << " to a practical value of " << 0.999 << endl << endl;
			risk = 0.999;
		}
		if (risk < 0.001)
		{
			f_rec << endl << "  ---  note: resetting risk value of " << risk << " to a practical value of " << 0.001 << endl << endl;
			risk = 0.001;
		}

	    probit_val = get_probit();

		if (std_weights)
		{
			double std, var;
			cout << "++opt_std_weights = True, using weights as chance constraint uncertainty" << endl;
			f_rec << "++opt_std_weights = True, using the following weights as prior and posterior chance constraint uncertainty: " << endl;
			f_rec << setw(25) << "model-based constraint" << setw(15) << "std (weight)" << setw(15) << "variance" << endl;
			for (auto &cname : ctl_ord_obs_constraint_names)
			{
				std = pest_scenario.get_observation_info_ptr()->get_weight(cname);
				var = std * std;
				f_rec << setw(25) << cname << setw(15) << std << setw(15) << var << endl;
				prior_const_var[cname] = var;
				post_const_var[cname] = var;


			}
		}
		else
		{
			//make sure there is at least one non-decision var adjustable parameter
			vector<string>::iterator start = ctl_ord_dec_var_names.begin();
			vector<string>::iterator end = ctl_ord_dec_var_names.end();
			set<string> dec_set(ctl_ord_dec_var_names.begin(), ctl_ord_dec_var_names.end());
			for (auto &name : pest_scenario.get_ctl_ordered_par_names())
			{
				//if this parameter is not a decision var
				//if (find(start, end, name) == end)
				if (dec_set.find(name) == dec_set.end())
				{
					ParameterRec::TRAN_TYPE tt = pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(name)->tranform_type;
					if ((tt == ParameterRec::TRAN_TYPE::LOG) || (tt == ParameterRec::TRAN_TYPE::NONE))
						adj_par_names.push_back(name);
				}
			}
			if (num_adj_pars() == 0)
				throw_sequentialLP_error("++opt_risk != 0.5, but no adjustable parameters found in control file");

			//look for non-zero weighted obs
			start = ctl_ord_obs_constraint_names.begin();
			end = ctl_ord_obs_constraint_names.end();
			for (auto &name : pest_scenario.get_ctl_ordered_obs_names())
			{
				if (find(start, end, name) == end)
				{
					if (pest_scenario.get_ctl_observation_info().get_observation_rec_ptr(name)->weight > 0.0)
						nz_obs_names.push_back(name);
				}
			}


			//string parcov_filename = pest_scenario.get_pestpp_options().get_parcov_filename();
			////build the adjustable parameter parcov
			////from filename

			//if (parcov_filename.size() > 0)
			//{
			//	//throw_sequentialLP_error("parcov from filename not implemented");
			//	Covariance temp_parcov ;
			//	temp_parcov.from_ascii(parcov_filename);
			//	//check that all adj par names in temp_parcov
			//	vector<string> temp_names = temp_parcov.get_col_names();
			//	start = temp_names.begin();
			//	end = temp_names.end();
			//	vector<string> missing;
			//	for (auto &name : adj_par_names)
			//		if (find(start, end, name) == end)
			//			missing.push_back(name);
			//	if (missing.size() > 0)
			//		throw_sequentialLP_error("the following adjustable parameters were not found in the ++parcov_filename covaraince matrix: ", missing);

			//	parcov = temp_parcov.get(adj_par_names);

			//}
			////from parameter bounds
			//else
			//{
			//	parcov.from_parameter_bounds(adj_par_names, pest_scenario.get_ctl_parameter_info());
			//}
			vector<string> drop;
			set<string> sadj(adj_par_names.begin(), adj_par_names.end());
			for (auto &n : parcov.get_row_names())
			{
				if (sadj.find(n) == sadj.end())
					drop.push_back(n);
			}
			parcov.drop(drop);


			//build the nz_obs obs_cov
			if (num_nz_obs() != 0)
				obscov.from_observation_weights(nz_obs_names, pest_scenario.get_ctl_observation_info(), vector<string>(), null_prior);
		}
	}
	else use_chance = false;


	//--------------------------------
	//  ---  objective function  ---
	//--------------------------------

	//initialize the objective function
	obj_func_str = pest_scenario.get_pestpp_options().get_opt_obj_func();
	obj_sense = (pest_scenario.get_pestpp_options().get_opt_direction() == 1) ? "minimize" : "maximize";


	//check if the obj_str is an observation
	use_obj_obs = false;
	if (pest_scenario.get_ctl_observations().find(obj_func_str) != pest_scenario.get_ctl_observations().end())
	{
		use_obj_obs = true;
		obj_obs = obj_func_str;
		//check
		set<string> names(ctl_ord_obs_constraint_names.begin(), ctl_ord_obs_constraint_names.end());
		if (names.find(obj_obs) != names.end())
		{
			throw runtime_error("objective function obs is a constraint, #sad");
		}
		names.clear();
		names.insert(nz_obs_names.begin(),nz_obs_names.end());
		if (names.find(obj_obs) != names.end())
		{
			throw runtime_error("objective function obs has non-zero weight and chance constraints are active");
		}

	}

	else
	{
		if (obj_func_str.size() == 0)
		{
			f_rec << " warning: no ++opt_objective_function-->forming a generic objective function (1.0 coef for each decision var)" << endl;
			for (auto &name : ctl_ord_dec_var_names)
				obj_func_coef_map[name] = 1.0;
		}

		//or if it is a prior info equation
		else if (pest_scenario.get_prior_info().find(obj_func_str) != pest_scenario.get_prior_info().end())
		{
			obj_func_coef_map = pest_scenario.get_prior_info().get_pi_rec_ptr(obj_func_str).get_atom_factors();
			//throw_sequentialLP_error("prior-information-based objective function not implemented");
		}
		else
		{
			//check if this obj_str is a filename
			ifstream if_obj(obj_func_str);
			if (!if_obj.good())
				throw_sequentialLP_error("unrecognized ++opt_objective_function arg: " + obj_func_str);
			else
				obj_func_coef_map = pest_utils::read_twocol_ascii_to_map(obj_func_str);
		}


		//check that all obj_coefs are decsision vars
		vector<string> missing_vars;
		for (auto &coef : obj_func_coef_map)
			if (find(ctl_ord_dec_var_names.begin(), ctl_ord_dec_var_names.end(), coef.first) == ctl_ord_dec_var_names.end())
				missing_vars.push_back(coef.first);
		if (missing_vars.size() > 0)
			throw_sequentialLP_error("the following objective function components are not decision variables: ", missing_vars);

	}


	jco.set_base_numeric_pars(all_pars_and_dec_vars);
	jco.set_base_sim_obs(pest_scenario.get_ctl_observations());
	if (pest_scenario.get_pestpp_options().get_opt_coin_log())
		model.setLogLevel(4 + 8 + 16 + 32);
	initial_report();
	return;
}

void sequentialLP::calc_chance_constraint_offsets()
{
	ofstream & f_rec = file_mgr_ptr->rec_ofstream();
	/*if ((slp_iter != 1) && ((slp_iter+1) % pest_scenario.get_pestpp_options().get_opt_recalc_fosm_every() != 0))
	{

		f_rec << endl << "  ---  reusing fosm offsets from previous iteration  ---  " << endl;
		cout << endl << "  ---  reusing fosm offsets from previous iteration  ---  " << endl;
		return;
	}*/
	cout << "  ---  calculating FOSM-based chance constraint components  ---  " << endl;
	f_rec << "  ---  calculating FOSM-based chance constraint components  ---  " << endl;
	prior_constraint_offset.clear();
	prior_constraint_stdev.clear();
	post_constraint_offset.clear();
	post_constraint_stdev.clear();
	if ((!std_weights) && ((slp_iter == 1) || ((slp_iter + 1) % pest_scenario.get_pestpp_options().get_opt_recalc_fosm_every() == 0)))
	{
		//the rows of the fosm jacobian include nonzero weight obs (for schur comp)
		//plus the names of the names of constraints, which get treated as forecasts
		vector<string> fosm_row_names(nz_obs_names);
		fosm_row_names.insert(fosm_row_names.end(), ctl_ord_obs_constraint_names.begin(), ctl_ord_obs_constraint_names.end());

		//extract the part of the full jco we need for fosm
		Eigen::SparseMatrix<double> fosm_mat = jco.get_matrix(fosm_row_names, adj_par_names);

		Mat fosm_jco(fosm_row_names, adj_par_names, fosm_mat);

		//create a linear object
		Logger logger(file_mgr_ptr->get_ofstream("pfm"), false);
		linear_analysis la(&fosm_jco, &pest_scenario, &obscov, &logger);

		//set the prior parameter covariance matrix
		la.set_parcov(&parcov);

		//set the predictions (the constraints)
		la.set_predictions(ctl_ord_obs_constraint_names);

		//get the prior and posterior variance of the constraints
		prior_const_var = la.prior_prediction_variance();

		//if at least one nz obs was found, then use schur complment, otherwise,
		//just use the prior constraint uncertainty
		if (num_nz_obs() > 0)
			post_const_var = la.posterior_prediction_variance();
		else
			post_const_var = prior_const_var;
	}
	//work out the offset for each constraint
	//and set the values in the constraints_fosm Obseravtions
	double new_constraint_val, old_constraint_val, required_val;
	double pr_offset, pt_offset;
	//constraints_fosm.clear();
	Observations iter_fosm;

	//map<string,double> out_of_bounds;
	for (auto &name : ctl_ord_obs_constraint_names)
	{
		prior_constraint_stdev[name] = sqrt(prior_const_var[name]);
		post_constraint_stdev[name] = sqrt(post_const_var[name]);
		pr_offset = probit_val * prior_constraint_stdev[name];
		pt_offset = probit_val * post_constraint_stdev[name];
		//important: using the initial simulated constraint values
		old_constraint_val = constraints_sim[name];
		//old_constraint_val = constraints_sim_initial[name];
		required_val = constraints_obs[name];

		//if less_than constraint, then add to the sim value, to move positive
		// WRT the required constraint value
		if (constraint_sense_map[name] == ConstraintSense::less_than)
		{

			//if the old val was in bounds and fosm is out
			/*if ((old_constraint_val <= required_val) && (new_constraint_val > required_val))
				out_of_bounds[name] = abs(new_constraint_val - required_val);*/
			/*if (old_constraint_val < 0.0)
			{
				new_constraint_val = old_constraint_val - pt_offset;
				post_constraint_offset[name] = -pt_offset;
				prior_constraint_offset[name] = -pr_offset;
			}
			else
			{*/
				new_constraint_val = old_constraint_val + pt_offset;
				post_constraint_offset[name] = pt_offset;
				prior_constraint_offset[name] = pr_offset;
			//}

		}
		//if greater_than constraint, the substract from the sim value to move
		//negative WRT the required constraint value
		else if (constraint_sense_map[name] == ConstraintSense::greater_than)
		{

			//if the old val was inbounds and fosm is out
			/*if ((old_constraint_val >= required_val) && (new_constraint_val < required_val))
				out_of_bounds[name] = abs(new_constraint_val - required_val);*/
			/*if (old_constraint_val < 0.0)
			{
				new_constraint_val = old_constraint_val + pt_offset;
				post_constraint_offset[name] = pt_offset;
				prior_constraint_offset[name] = pr_offset;
			}
			else
			{*/
				new_constraint_val = old_constraint_val - pt_offset;
				post_constraint_offset[name] = -pt_offset;
				prior_constraint_offset[name] = -pr_offset;
			//}
		}

		else
			new_constraint_val = constraints_sim[name];
		iter_fosm.insert(name, new_constraint_val);
	}

	vector<string> names = iter_fosm.get_keys();
	constraints_fosm.update_without_clear(names, iter_fosm.get_data_vec(names));
	cout << "  ---   done with FOSM-based chance constraint calculations  ---  " << endl << endl;
	f_rec << "  ---   done with FOSM-based chance constraint calculations  ---  " << endl << endl;
	return;
}

void sequentialLP::build_constraint_bound_arrays()
{

	//if needed update the fosm related attributes
	//before calculating the residual vector!
	if (use_chance)
		calc_chance_constraint_offsets();

	//these residuals will include FOSM-based posterior offsets
	vector<double> residuals = get_constraint_residual_vec();

	for (int i = 0; i < num_obs_constraints(); ++i)
	{
		string name = ctl_ord_obs_constraint_names[i];
		if (constraint_sense_map[name] == ConstraintSense::less_than)
			constraint_ub[i] = residuals[i];
		else
			constraint_ub[i] = COIN_DBL_MAX;

		if (constraint_sense_map[name] == ConstraintSense::greater_than)
			constraint_lb[i] = residuals[i];
		else
			constraint_lb[i] = -COIN_DBL_MAX;
		if (constraint_sense_map[name] == ConstraintSense::equal_to)
		{
			constraint_ub[i] = residuals[i];
			constraint_lb[i] = residuals[i];
		}
	}

	int noc = num_obs_constraints();
	for (int i = 0; i < num_pi_constraints(); ++i)
	{
		string name = ctl_ord_pi_constraint_names[i];
		double residual = -constraints_pi.get_pi_rec_ptr(name).calc_residual(all_pars_and_dec_vars_initial);
		if (constraint_sense_map[name] == ConstraintSense::less_than)
			constraint_ub[i+noc] = residual;
		else
			constraint_ub[i+noc] = COIN_DBL_MAX;

		if (constraint_sense_map[name] == ConstraintSense::greater_than)
			constraint_lb[i+noc] = residual;
		else
			constraint_lb[i+noc] = -COIN_DBL_MAX;
		if (constraint_sense_map[name] == ConstraintSense::equal_to)
		{
			constraint_ub[i+noc] = residual;
			constraint_lb[i+noc] = residual;
		}
	}
	return;
}


void sequentialLP::build_obj_func_coef_array()
{
	ctl_ord_obj_func_coefs = new double[num_dec_vars()];
	double coef;
	map<string, double>::iterator end = obj_func_coef_map.end();
	int i = 0;
	for (auto &name : ctl_ord_dec_var_names)
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
	Parameters parlbnd = pest_scenario.get_ctl_parameter_info().get_low_bnd(ctl_ord_dec_var_names);
	Parameters parubnd = pest_scenario.get_ctl_parameter_info().get_up_bnd(ctl_ord_dec_var_names);
	for (int i = 0; i < num_dec_vars(); ++i)
	{
		dec_var_lb[i] = parlbnd.get_rec(ctl_ord_dec_var_names[i]) - all_pars_and_dec_vars.get_rec(ctl_ord_dec_var_names[i]);
		dec_var_ub[i] = parubnd.get_rec(ctl_ord_dec_var_names[i]) - all_pars_and_dec_vars.get_rec(ctl_ord_dec_var_names[i]);
	}
}

void sequentialLP::iter_solve()
{

	ofstream &f_rec = file_mgr_ptr->rec_ofstream();

	//convert Jacobian_1to1 to CoinPackedMatrix
	cout << "  ---  forming LP model  --- " << endl;
	CoinPackedMatrix matrix = jacobian_to_coinpackedmatrix();

	build_dec_var_bounds();

	//load the linear simplex model
	model.loadProblem(matrix, dec_var_lb, dec_var_ub, ctl_ord_obj_func_coefs, constraint_lb, constraint_ub);
	for (int i = 0; i < num_obs_constraints(); ++i)
		model.setRowName(i, ctl_ord_obs_constraint_names[i]);
	for (int i = 0; i < ctl_ord_pi_constraint_names.size(); ++i)
		model.setRowName(i+num_obs_constraints(), ctl_ord_pi_constraint_names[i]);
	for (int i = 0; i < num_dec_vars(); ++i)
		model.setColumnName(i, ctl_ord_dec_var_names[i]);

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

int sequentialLP::num_nz_pi_constraint_elements()
{
	int num = 0;

	for (auto &pi_name : ctl_ord_pi_constraint_names)
	{
		num += constraints_pi.get_pi_rec_ptr(pi_name).get_atom_factors().size();
	}
	return num;
}

CoinPackedMatrix sequentialLP::jacobian_to_coinpackedmatrix()
{

	Eigen::SparseMatrix<double> eig_ord_jco = jco.get_matrix(ctl_ord_obs_constraint_names, ctl_ord_dec_var_names);

	file_mgr_ptr->rec_ofstream() << "number of nonzero elements in response matrix: " << eig_ord_jco.nonZeros() << " of " << eig_ord_jco.size() << endl;
	cout << "number of nonzero elements in response matrix: " << eig_ord_jco.nonZeros() << " of " << eig_ord_jco.size() << endl;

	if (eig_ord_jco.nonZeros() == 0)
		throw_sequentialLP_error("response matrix all zeros...cannot continue");

	//cout << eig_ord_jco << endl;

	//the triplet elements to pass to the coinpackedmatrix constructor
	int num_elems = eig_ord_jco.nonZeros() + num_nz_pi_constraint_elements();
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
	int irow = num_obs_constraints();
	int jcol;
	vector<string>::iterator start = ctl_ord_dec_var_names.begin();
	vector<string>::iterator end = ctl_ord_dec_var_names.end();

	for (auto &pi_name : ctl_ord_pi_constraint_names)
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
//#ifdef _DEBUG
//	cout << eig_ord_jco << endl << endl;
//	matrix.dumpMatrix();
//#endif
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
		bool success = make_upgrade_run(all_pars_and_dec_vars_best, constraints_sim);
		if (!success)
		{
			throw_sequentialLP_error("error running model with best decision variable values");
		}
		//write 'best' rei
		of_wr.write_opt_constraint_rei(file_mgr_ptr->open_ofile_ext("res"), slp_iter, all_pars_and_dec_vars_best,
			pest_scenario.get_ctl_observations(), constraints_sim);
		file_mgr_ptr->close_file("res");
		//write 'best' parameters file
		/*of_wr.write_par(file_mgr_ptr->open_ofile_ext("par"), all_pars_and_dec_vars_best,
			*par_trans.get_offset_ptr(), *par_trans.get_scale_ptr());
		file_mgr_ptr->close_file("par");*/
	}
}

void sequentialLP::iter_postsolve()
{

	ofstream &f_rec = file_mgr_ptr->rec_ofstream();

	row_price = model.getRowPrice();

	//extract (optimal) decision vars
	//and track some info for convergence checking
	const double *dec_var_vals = model.getColSolution();
	double max_abs_dec_var_change = -1.0E+10;
	double max_abs_dec_var_val = -1.0E+10;

	double diff, val;
	Parameters upgrade_pars(all_pars_and_dec_vars);
	string name;
	for (int i = 0; i < num_dec_vars(); ++i)
	{
		name = ctl_ord_dec_var_names[i];
		val = all_pars_and_dec_vars[name];
		diff = abs(dec_var_vals[i] - all_pars_and_dec_vars[name]);
		upgrade_pars.update_rec(name,dec_var_vals[i] + val);

		max_abs_dec_var_change = (diff > max_abs_dec_var_change) ? diff : max_abs_dec_var_change;
		max_abs_dec_var_val = (abs(val) > max_abs_dec_var_val) ? val : max_abs_dec_var_val;
	}
	max_abs_dec_var_change /= max(max_abs_dec_var_val,1.0);

	//run the model with optimal decision var values

	Observations upgrade_obs = constraints_sim;
	if (!super_secret_option)
	{
		bool success = make_upgrade_run(upgrade_pars, upgrade_obs);

		f_rec << "  ---  processing results for iteration " << slp_iter << " LP solution  ---  " << endl << endl;
		postsolve_constraint_report(upgrade_obs, upgrade_pars);
	}
	double obj_val = model.getObjValue();

	pair<double,double> cur_new_obj = postsolve_decision_var_report(upgrade_pars);


	f_rec << endl << endl <<  "  ---  iteration " << slp_iter << " objective function value: " << setw(15) << cur_new_obj.second << "  ---  " << endl << endl;
	cout << endl << endl << "  ---  iteration " << slp_iter << " objective function value: " << setw(15) << cur_new_obj.second << "  ---  " << endl << endl;

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
			all_pars_and_dec_vars_best.update_without_clear(ctl_ord_dec_var_names, upgrade_pars.get_data_vec(ctl_ord_dec_var_names));
			obj_best = cur_new_obj.second;
		}
	}
	else
		if (cur_new_obj.second < obj_best)
		{
			all_pars_and_dec_vars_best.update_without_clear(ctl_ord_dec_var_names, upgrade_pars.get_data_vec(ctl_ord_dec_var_names));
			obj_best = cur_new_obj.second;
		}

	//check for changes for in constraints
	double max_abs_constraint_change = -1.0E+10;
	double max_abs_constraint_val = -1.0E+10;
	for (auto &name : ctl_ord_obs_constraint_names)
	{

		diff = abs(constraints_sim[name] - upgrade_obs[name]);
		max_abs_constraint_change = (diff > max_abs_constraint_change) ? diff : max_abs_constraint_change;
		max_abs_constraint_val = (val > max_abs_constraint_val) ? val : max_abs_constraint_val;
	}
	max_abs_constraint_change /= max(max_abs_constraint_val,1.0);

	pair<map<string,double>, map<string,double>> invalid_vars_const = postsolve_check(upgrade_obs, upgrade_pars);

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
	if (invalid_vars_const.first.size() > 0)
	{
		valid = false;
		f_rec << "-->the following decision variables are out of bounds (distance shown):" << endl;
		for (auto &name : invalid_vars_const.first)
			f_rec << "-->   " << name.first << setw(15) << name.second << endl;
	}
	if (invalid_vars_const.second.size() > 0)
	{
		valid = false;
		f_rec << "-->the following constraints are not satisfied (distance shown):" << endl;
		for (auto &name : invalid_vars_const.second)
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
	all_pars_and_dec_vars.update_without_clear(ctl_ord_dec_var_names, upgrade_pars.get_data_vec(ctl_ord_dec_var_names));

	//constraints_sim.update(ctl_ord_obs_constraint_names,upgrade_obs.get_data_vec(ctl_ord_obs_constraint_names));
	constraints_sim.update_without_clear(ctl_ord_obs_constraint_names, upgrade_obs.get_data_vec(ctl_ord_obs_constraint_names));
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
		set<string> ext_set(ctl_ord_ext_var_names.begin(), ctl_ord_ext_var_names.end());
		vector<string> missing;
		for (auto &name : ctl_ord_dec_var_names)
			//if this dec var is not in the jco and is not an external var
			//if ((find(start, end, name) == end) && (find(ext_start,ext_end,name) == ext_end))
			if ((names.find(name) == end) && (ext_set.find(name) == ext_set.end()))
				missing.push_back(name);
		if (missing.size() > 0)
			throw_sequentialLP_error("the following decision vars were not found in the jacobian " + basejac_filename + " : ", missing);

		for (auto &name : adj_par_names)
			//if (find(start, end, name) == end)
			if ((names.find(name) == end))
				missing.push_back(name);
		if (missing.size() > 0)
			throw_sequentialLP_error("the following adjustable parameters were not found in the jacobian " + basejac_filename + " : ", missing);


		names.clear();
		temp = jco.get_sim_obs_names();
		names.insert(temp.begin(),temp.end());
		start = names.begin();
		end = names.end();
		for (auto &name : ctl_ord_obs_constraint_names)
			if (find(start, end, name) == end)
				missing.push_back(name);
		if (missing.size() > 0)
			throw_sequentialLP_error("the following constraints were not found in the jacobian " + basejac_filename + " : ", missing);

		for (auto &name : nz_obs_names)
			if (find(start, end, name) == end)
				missing.push_back(name);
		if (missing.size() > 0)
			throw_sequentialLP_error("the following non-zero weight observations were not found in the jacobian " + basejac_filename + " : ", missing);

		string res_filename = pest_scenario.get_pestpp_options().get_hotstart_resfile();
		if (!res_filename.empty())
		{
			stringstream message;
			message << "  reading  residual file " << res_filename << " for hot-start...";
			cout << message.str();
			f_rec << message.str();
			for (auto &oname : pest_scenario.get_ctl_ordered_obs_names())
				constraints_sim[oname] = -1.0e+30;
			pest_utils::read_res(res_filename, constraints_sim);
			f_rec << "done" << endl;
			cout << "done" << endl;

		}
		else
		{
			//make the intial base run
			cout << "  ---  running the model once with initial decision variables  ---  " << endl;
			int run_id = run_mgr_ptr->add_run(par_trans.ctl2model_cp(all_pars_and_dec_vars));
			run_mgr_ptr->run();
			Parameters pars;
			bool success = run_mgr_ptr->get_run(run_id, pars, constraints_sim);
			if (!success)
				throw_sequentialLP_error("initial (base) run with initial decision vars failed...cannot continue");

		}
	}

	//otherwise, fill the jacobian
	else
	{
		set<string> out_of_bounds;
		vector<string> names_to_run;
		for (auto &name : ctl_ord_dec_var_names)
			if (find(ctl_ord_ext_var_names.begin(), ctl_ord_ext_var_names.end(), name) == ctl_ord_ext_var_names.end())
				names_to_run.push_back(name);

		if ((!std_weights) && ((slp_iter == 1) || ((slp_iter+1) % pest_scenario.get_pestpp_options().get_opt_recalc_fosm_every() == 0)))
		{
			names_to_run.insert(names_to_run.end(), adj_par_names.begin(), adj_par_names.end());
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
			for (auto &name : ctl_ord_dec_var_names)
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
		//ParameterGroupInfo pinfo = pest_scenario.get_base_group_info();
		//for (auto &name : pinfo.get_group_names())
		//{
		//	//pinfo.get_group_by_groupname_4_mod(name)->derinc *= 0.75;
		//	double inc = pinfo.get_group_by_groupname_4_mod(name)->derinc;
		//	cout << name << ',' << inc << endl;
		//}

		bool init_obs = false;
		if (slp_iter == 1) init_obs = true;
		bool success = jco.build_runs(all_pars_and_dec_vars, constraints_sim, names_to_run, par_trans,
			pest_scenario.get_base_group_info(), pest_scenario.get_ctl_parameter_info(),
			*run_mgr_ptr, out_of_bounds,false,init_obs);
		if (!success)
		{
			const set<string> failed = jco.get_failed_parameter_names();
			throw_sequentialLP_error("failed to calc derviatives for the following decision vars: ", failed);
		}

		jco.make_runs(*run_mgr_ptr);
		set<int> failed = run_mgr_ptr->get_failed_run_ids();

		//process the remaining responses
		success = jco.process_runs(par_trans, pest_scenario.get_base_group_info(), *run_mgr_ptr, *null_prior, false);
		if (!success)
			throw_sequentialLP_error("error processing response matrix runs ", jco.get_failed_parameter_names());

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
			Parameters temp_pars;
			Observations temp_obs;
			run_mgr_ptr->get_run(0, temp_pars, constraints_sim, false);

		}
		//par_trans.model2ctl_ip(temp_pars);


	}
	//if this is the first time through, set the initial constraint simulated values
	if (slp_iter == 1)
		constraints_sim_initial = Observations(constraints_sim);

	if (use_obj_obs)
	{
		obj_func_coef_map.clear();
		vector<string> onames = jco.get_sim_obs_names();
		vector<string> pnames = jco.get_base_numeric_par_names();
		set<string> sdecvar(ctl_ord_dec_var_names.begin(), ctl_ord_dec_var_names.end());

		int idx = find(onames.begin(), onames.end(), obj_obs) - onames.begin();
		if (idx >= onames.size())
			throw runtime_error("obj function obs not found in jco row names, #sad");
		Eigen::VectorXd vec = Eigen::VectorXd(jco.get_matrix_ptr()->row(idx));
		for (int i = 0; i < vec.size(); i++)
		{
			if (sdecvar.find(pnames[i]) == sdecvar.end())
				continue;
			obj_func_coef_map[pnames[i]] = vec[i];
		}
		f_rec << "  ---  objective function coefficients for iteration " << slp_iter << "  ---  " << endl;
		double obj = obj_func_report();
		if (slp_iter == 1)
		{
			obj_init = obj;
			f_rec << endl << "  ---  objective function value (using initial dec var values): " << obj_init << endl << endl;
			cout << endl << "  ---  objective function value (using initial dec var values): " << obj_init << endl << endl;
		}

	}

	//set/update the constraint bound arrays
	build_constraint_bound_arrays();

	if (use_chance)
		presolve_fosm_report();

	//report to rec file
	presolve_constraint_report();

	//build the objective function
	build_obj_func_coef_array();

	stringstream ss;
	ss << "jcb." << slp_iter << ".rei";
	of_wr.write_opt_constraint_rei(file_mgr_ptr->open_ofile_ext(ss.str()), slp_iter, pest_scenario.get_ctl_parameters(),
		pest_scenario.get_ctl_observations(), constraints_sim);

	return;
}



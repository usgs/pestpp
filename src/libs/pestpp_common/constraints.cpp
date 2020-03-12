
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <random>
#include <iterator>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>

#include "SVDPackage.h"
#include "Pest.h"
#include "utilities.h"
#include "covariance.h"
#include "FileManager.h"
#include "constraints.h"

using namespace std;

Constraints::Constraints(Pest& _pest_scenario, FileManager* _file_mgr_ptr, OutputFileWriter& _of_wr, vector<string>& ctl_ord_dec_var_names)
	:pest_scenario(_pest_scenario), file_mgr_ptr(_file_mgr_ptr), of_wr(_of_wr), jco(*_file_mgr_ptr,_of_wr)
{
	ofstream& f_rec = file_mgr_ptr->rec_ofstream();
	if (pest_scenario.get_pestpp_options().get_opt_include_bnd_pi())
	{
		PriorInformation* pi_ptr = pest_scenario.get_prior_info_ptr();
		ParameterInfo par_info = pest_scenario.get_ctl_parameter_info();
		stringstream ss;
		string s;
		for (auto& dname : ctl_ord_dec_var_names)
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


	//---------------------------
	//  ---  obs constraints  ---
	//---------------------------
	//set the two constraints attribs and ordered constraint name vec
	vector<string> constraint_groups = pest_scenario.get_pestpp_options().get_opt_constraint_groups();
	//if not constraint groups set as ++ arg, then look for all compatible obs groups
	if (constraint_groups.size() == 0)
	{
		for (auto& group : pest_scenario.get_ctl_ordered_obs_group_names())
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
		throw_constraints_error("no viable observation or prior information constraint groups found.  Constraint group names must start with the following: {'l_','less','g_','greater'}");

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
			throw_constraints_error("the following ++opt_constraint_groups were not found: ", missing);

		//find the observations in constraints groups
		ObservationInfo oinfo = pest_scenario.get_ctl_observation_info();
		string group;
		double weight;
		end = constraint_groups.end();
		start = constraint_groups.begin();
		int nobszw = 0;
		for (auto& obs_name : pest_scenario.get_ctl_ordered_obs_names())
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
		for (auto& pi_name : pest_scenario.get_ctl_ordered_pi_names())
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
			for (auto& pi_name : ctl_ord_pi_constraint_names)
			{
				pi_rec = constraints_pi.get_pi_rec_ptr(pi_name);
				missing.clear();
				for (auto& pi_factor : pi_rec.get_atom_factors())
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
			for (auto& missing_pi : missing_map)
			{
				ss << missing_pi.first << ": ";
				for (auto& par_name : missing_pi.second)
					ss << par_name << ",";
			}
			throw_constraints_error("errors in prior information constraints:" + ss.str());
		}

		//TODO: investigate a pi constraint only formulation
		if (num_obs_constraints() == 0)
			throw_constraints_error("no constraints found in groups: ", constraint_groups);
	}


	//if not ++opt_constraint_names was passed, use all observations and prior information as constraints
	/*else
	{
		ctl_ord_obs_constraint_names = pest_scenario.get_ctl_ordered_obs_names();
		ctl_ord_pi_constraint_names = pest_scenario.get_ctl_ordered_pi_names();
	}*/


	constraints_obs = pest_scenario.get_ctl_observations().get_subset(ctl_ord_obs_constraint_names.begin(), ctl_ord_obs_constraint_names.end());
	//constraints_sim = Observations(constraints_obs);

	//build map of obs constraint sense
	vector<string> problem_constraints;
	for (auto& name : ctl_ord_obs_constraint_names)
	{
		string group = pest_scenario.get_ctl_observation_info().get_observation_rec_ptr(name)->group;
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
		throw_constraints_error("the following obs constraints do not have a correct group name prefix {'l_','less','g_','greater','e_','equal'}: ", problem_constraints);
	}

	//build map of pi constraint sense
	problem_constraints.clear();
	for (auto& name : ctl_ord_pi_constraint_names)
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
		throw_constraints_error("the following prior info constraints do not have a correct group name prefix {'l_','less','g_','greater','e_','equal'}: ", problem_constraints);
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
			throw_constraints_error("++opt_risk parameter must between 0.0 and 1.0");

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
			for (auto& cname : ctl_ord_obs_constraint_names)
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
			for (auto& name : pest_scenario.get_ctl_ordered_par_names())
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
			if (adj_par_names.size() == 0)
				throw_constraints_error("++opt_risk != 0.5, but no adjustable parameters found in control file");

			//look for non-zero weighted obs
			start = ctl_ord_obs_constraint_names.begin();
			end = ctl_ord_obs_constraint_names.end();
			for (auto& name : pest_scenario.get_ctl_ordered_obs_names())
			{
				if (find(start, end, name) == end)
				{
					if (pest_scenario.get_ctl_observation_info().get_observation_rec_ptr(name)->weight > 0.0)
						nz_obs_names.push_back(name);
				}
			}

			vector<string> drop;
			set<string> sadj(adj_par_names.begin(), adj_par_names.end());
			for (auto& n : parcov.get_row_names())
			{
				if (sadj.find(n) == sadj.end())
					drop.push_back(n);
			}
			parcov.drop(drop);

			//build the nz_obs obs_cov
			if (nz_obs_names.size() != 0)
			{
				ObservationInfo oi = pest_scenario.get_ctl_observation_info();
				obscov.from_observation_weights(nz_obs_names, oi, vector<string>(), null_prior);
			}
		}
	}
	else use_chance = false;


}

pair<Constraints::ConstraintSense, string> Constraints::get_sense_from_group_name(const string& name)
{

	if ((name.compare(0, 2, "L_") == 0) || (name.compare(0, 4, "LESS") == 0))
		return pair<ConstraintSense, string>(ConstraintSense::less_than, "less_than");
	else if ((name.compare(0, 2, "G_") == 0) || (name.compare(0, 7, "GREATER") == 0))
		return pair<ConstraintSense, string>(ConstraintSense::greater_than, "greater_than");
	else if ((name.compare(0, 2, "N_") == 0) || (name.compare(0, 2, "E_") == 0) || (name.compare(0, 5, "EQUAL") == 0))
		return pair<ConstraintSense, string>(ConstraintSense::equal_to, "equal_to");
	else
		return pair<ConstraintSense, string>(ConstraintSense::undefined, "undefined");
}


void Constraints::throw_constraints_error(string message, const vector<string>& messages)
{
	stringstream ss;
	for (auto& mess : messages)
		ss << mess + ',';
	throw_constraints_error(message + ss.str());
}

void Constraints::throw_constraints_error(string message, const set<string>& messages)
{
	stringstream ss;
	for (auto& mess : messages)
		ss << mess + ',';
	throw_constraints_error(message + ss.str());
}

void Constraints::throw_constraints_error(string message)
{
	string error_message = "error in Constraints: " + message;
	file_mgr_ptr->rec_ofstream() << error_message << endl;
	file_mgr_ptr->close_file("rec");
	cout << endl << endl << error_message << endl << endl;
	throw runtime_error(error_message);
}

double  ErfInv2(double x)
{
	float tt1, tt2, lnx, sgn;
	sgn = (x < 0) ? -1.0f : 1.0f;

	x = (1 - x) * (1 + x);        // x = 1 - x*x;
	lnx = logf(x);
	//double PI = 3.14159265358979323846;
	tt1 = 2 / (M_PI * 0.147) + 0.5f * lnx;
	tt2 = 1 / (0.147) * lnx;

	return(sgn * sqrtf(-tt1 + sqrtf(tt1 * tt1 - tt2)));
}


double Constraints::get_probit()
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
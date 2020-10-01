
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
#include "linear_analysis.h"

using namespace std;

OptObjFunc::OptObjFunc(Pest& _pest_scenario, FileManager* _file_mgr_ptr, PerformanceLog& _pfm):
	pest_scenario(_pest_scenario), file_mgr_ptr(_file_mgr_ptr), pfm(_pfm)
{

}

void OptObjFunc::update_coef_map_from_jacobian(Jacobian& jco)
{
	if (!use_obj_obs)
		return;

	obj_func_coef_map.clear();
	vector<string> onames = jco.get_sim_obs_names();
	vector<string> pnames = jco.get_base_numeric_par_names();
	set<string> sdecvar(dv_names.begin(), dv_names.end());

	int idx = find(onames.begin(), onames.end(), obj_obs) - onames.begin();
	if (idx >= onames.size())
		throw_optobjfunc_error("obj function obs name '" + obj_obs +"' not found in jco row names, #sad");
	Eigen::VectorXd vec = Eigen::VectorXd(jco.get_matrix_ptr()->row(idx));
	for (int i = 0; i < vec.size(); i++)
	{
		if (sdecvar.find(pnames[i]) == sdecvar.end())
			continue;
		obj_func_coef_map[pnames[i]] = vec[i];
	}
}

double OptObjFunc::get_obj_func_value(Parameters& pars, Observations& obs)
{
	double obj_val = 0.0;
	if (use_obj_obs)

	{
		if (obj_func_coef_map.size() == 0)
			obj_val = obs.get_rec(obj_obs);
		else
		{
			for (auto dv_name : dv_names)
			{
				obj_val += pars.get_rec(dv_name) * obj_func_coef_map[dv_name];
			}
		}
	}
	else if (obj_func_coef_map.size() == 0)
		throw_optobjfunc_error("get_obj_func_value: not using observation-based objective and obj coef map is empty");
	else
	{
		for (auto dv_name : dv_names)
		{
			obj_val += pars.get_rec(dv_name) * obj_func_coef_map[dv_name];
		}
	}
	return obj_val;
}

void OptObjFunc::report()
{
	ofstream& f_rec = file_mgr_ptr->rec_ofstream();
	map<string, double>::iterator end = obj_func_coef_map.end();
	vector<string> missing;
	if (use_obj_obs)
		f_rec << "objective function coefficients defined by observation: " << pest_scenario.get_pestpp_options().get_opt_obj_func() << endl;
	else
	{
		f_rec << "  ---  objective function coefficients  ---  " << endl;
		vector<string> missing;
		ofstream& f_rec = file_mgr_ptr->rec_ofstream();
		map<string, double>::iterator end = obj_func_coef_map.end();
		f_rec << setw(20) << left << "name" << setw(25) << "obj func coefficient" << endl;
		for (auto& name : dv_names)
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

		if (missing.size() > 0)
		{
			f_rec << endl << endl << "WARNING: the following decision variables have '0.0' objective function coef:" << endl;
			cout << endl << endl << "WARNING: the following decision variables have '0.0' objective function coef:" << endl;

			for (auto& name : missing)
			{
				f_rec << "    " << name << endl;
				f_rec << "    " << name << endl;
			}
		}
	}
}

void OptObjFunc::throw_optobjfunc_error(string message)
{
	string error_message = "error in sequentialLP process: " + message;
	file_mgr_ptr->rec_ofstream() << error_message << endl;
	file_mgr_ptr->close_file("rec");
	cout << endl << endl << error_message << endl << endl;
	throw runtime_error(error_message);
}

void OptObjFunc::initialize(vector<string> _constraint_names, vector<string> _dv_names)
{
	//initialize the objective function
	obj_func_str = pest_scenario.get_pestpp_options().get_opt_obj_func();
	obj_sense = (pest_scenario.get_pestpp_options().get_opt_direction() == 1) ? "minimize" : "maximize";

	ofstream& f_rec = file_mgr_ptr->rec_ofstream();

	dv_names = _dv_names;
	constraint_names = _constraint_names;

	//check if the obj_str is an observation
	use_obj_obs = false;
	if (pest_scenario.get_ctl_observations().find(obj_func_str) != pest_scenario.get_ctl_observations().end())
	{
		use_obj_obs = true;
		obj_obs = obj_func_str;
		//check
		set<string> names(constraint_names.begin(), constraint_names.end());
		if (names.find(obj_obs) != names.end())
		{
			throw runtime_error("objective function obs is a constraint, #sad");
		}
		names.clear();
		vector<string> cnames = pest_scenario.get_ctl_ordered_nz_obs_names();
		names.insert(cnames.begin(), cnames.end());
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
			for (auto& name : dv_names)
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
				throw_optobjfunc_error("unrecognized ++opt_objective_function arg: " + obj_func_str);
			else
				obj_func_coef_map = pest_utils::read_twocol_ascii_to_map(obj_func_str);
		}


		//check that all obj_coefs are decsision vars
		vector<string> missing_vars;
		set<string> s_dv_names(dv_names.begin(), dv_names.end());
		for (auto& coef : obj_func_coef_map)
			if (s_dv_names.find(coef.first) == s_dv_names.end())
				missing_vars.push_back(coef.first);
		if (missing_vars.size() > 0)
		{
			stringstream ss;
			ss << "the following objective function components are not decision variables: ";
			for (auto m : missing_vars)
				ss << m << ",";
			throw_optobjfunc_error(ss.str());
		}
			

	}
}


Constraints::Constraints(Pest& _pest_scenario, FileManager* _file_mgr_ptr, OutputFileWriter& _of_wr, PerformanceLog& _pfm)
	:pest_scenario(_pest_scenario), file_mgr_ptr(_file_mgr_ptr), of_wr(_of_wr), pfm(_pfm), jco(*_file_mgr_ptr, _of_wr)
{
	;
}


void Constraints::initialize(vector<string>& ctl_ord_dec_var_names, double _dbl_max)
{
	/* initialize the constraint class, most of this is to deal with chances*/
	dec_var_names = ctl_ord_dec_var_names;
	//store a pointer back to the original simulated constrain values 
	//this is so that as these are updated, the constraints instance
	//has access to the recent
	//current_constraints_sim_ptr = _current_constraints_sim_ptr;
	//same for dec var values
	//current_pars_and_dec_vars_ptr = _current_pars_and_dec_vars_ptr;
	//check for a stack size arg - if postive, then use stacks for chances
	stack_size = pest_scenario.get_pestpp_options().get_opt_stack_size();
	//an existing parameter stack for chances
	string par_stack_name = pest_scenario.get_pestpp_options().get_opt_par_stack();
	//maybe even an existing observations (e.g. constraints) stack!
	string obs_stack_name = pest_scenario.get_pestpp_options().get_opt_obs_stack();
	//by default, we want to use fosm for chances, but it any
	//of those stack options were passed, then use stacks instead
	use_fosm = true;
	if ((stack_size > 0) || (par_stack_name.size() > 0) || (obs_stack_name.size() > 0))
		use_fosm = false;
	//initialize the stack constainers (ensemble class instances)
	stack_pe.set_pest_scenario(&pest_scenario);
	stack_pe.set_rand_gen(&rand_gen);
	stack_oe.set_pest_scenario(&pest_scenario);
	stack_oe.set_rand_gen(&rand_gen);
	//initialize some more things
	dbl_max = _dbl_max;
	rand_gen = std::mt19937(pest_scenario.get_pestpp_options().get_random_seed());
	ctl_ord_obs_constraint_names.clear();
	ctl_ord_pi_constraint_names.clear();
	constraint_sense_map.clear();
	constraint_sense_name.clear();
	nz_obs_names.clear();
	adj_par_names.clear();
	prior_const_var.clear();
	post_const_var.clear();

	ofstream& f_rec = file_mgr_ptr->rec_ofstream();
	pfm.log_event("initializing constraints");
	//for the LP solve, we need to augment the PI based constraints with dec var bounds
	//this is to make sure we dont end up with an unbounded problem
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
		throw_constraints_error("no viable observation or prior information constraint/objective groups found.  Constraint/objective group names must start with the following: {'l_','less','g_','greater'}");

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
					f_rec << "Warning: observation constraint/objective " << obs_name << " has 0.0 weight, skipping" << endl;
				}
				else
				{
					ctl_ord_obs_constraint_names.push_back(obs_name);
				}

			}
		}
		if (nobszw > 0)
			cout << "Warning: " << nobszw << " of the observation constraints/objectives (see rec file for list) have 0.0 weight, skipping" << endl;

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
			throw_constraints_error("no constraints/objectives found in groups: ", constraint_groups);
	}


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
		throw_constraints_error("the following obs constraints/objectives do not have a correct group name prefix {'l_','less','g_','greater','e_','equal'}: ", problem_constraints);
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


	//------------------------------------------
	//  ---  chance constratints  ---
	//------------------------------------------
	risk = pest_scenario.get_pestpp_options().get_opt_risk();
	if (risk != 0.5)
	{
		pfm.log_event("initializing chance constraints/objectives");
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
		//if the std weight options was selected, use it - it overrides all other options
		if (std_weights)
		{
			double std, var;
			cout << "++opt_std_weights = True, using weights as chance constraint/objective uncertainty" << endl;
			f_rec << "++opt_std_weights = True, using the following weights as prior and posterior chance constraint/objective uncertainty: " << endl;
			f_rec << setw(25) << "observation constraint/objective" << setw(15) << "std (weight)" << setw(15) << "variance" << endl;
			for (auto& cname : ctl_ord_obs_constraint_names)
			{
				std = pest_scenario.get_observation_info_ptr()->get_weight(cname);
				var = std * std;
				f_rec << setw(25) << cname << setw(15) << std << setw(15) << var << endl;
				prior_const_var[cname] = var;
				post_const_var[cname] = var;
			}
		}
		//otherwise, if no stack options, use fosm
		else if (use_fosm)
		{
			//make sure there is at least one non-decision var adjustable parameter
			if (adj_par_names.size() == 0)
				throw_constraints_error("++opt_risk != 0.5, but no adjustable parameters found in control file");

			//look for non-zero weighted obs
			vector<string>::iterator start = ctl_ord_obs_constraint_names.begin();
			vector<string>::iterator end = ctl_ord_obs_constraint_names.end();
			for (auto& name : pest_scenario.get_ctl_ordered_obs_names())
			{
				if (find(start, end, name) == end)
				{
					if (pest_scenario.get_ctl_observation_info().get_observation_rec_ptr(name)->weight > 0.0)
						nz_obs_names.push_back(name);
				}
			}

			f_rec << "  using FOSM-based chance constraints/objectives with risk = " << risk << endl;

			parcov.try_from(pest_scenario, *file_mgr_ptr);
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
				map<string, double> obs_std = pest_scenario.get_ext_file_double_map("observation data external", "standard_deviation");
				obscov.from_observation_weights(file_mgr_ptr->rec_ofstream(), nz_obs_names, oi, vector<string>(), null_prior, obs_std);
			}
		}
		//otherwise, stack time baby!
		else
		{
			//make sure there is at least one non-decision var adjustable parameter
			if (adj_par_names.size() == 0)
				throw_constraints_error("++opt_risk != 0.5, but no adjustable parameters found in control file");
			bool size_passed = true;
			set<string> passed = pest_scenario.get_pestpp_options().get_passed_args();
			if (passed.find("OPT_STACK_SIZE") == passed.end())
				size_passed = false;
			string par_csv = pest_scenario.get_pestpp_options().get_opt_par_stack();
			
			//if a par stack wasnt passed, draw one
			if (par_csv.size() == 0)
			{
				if ((stack_size == 0) && (obs_stack_name.size() == 0))
					throw_constraints_error("++opt_stack_size is zero");
				else if (stack_size > 0)
				{
					f_rec << "drawing " << stack_size << "stack realizations" << endl;
					pfm.log_event("loading parcov");
					parcov.try_from(pest_scenario, *file_mgr_ptr);
					pfm.log_event("drawing stack realizations");
					stack_pe.draw(stack_size, pest_scenario.get_ctl_parameters(), parcov, &pfm,
						pest_scenario.get_pestpp_options().get_ies_verbose_level(), file_mgr_ptr->rec_ofstream());
				}
			}
			//otherwise, load the par stack
			else
			{
				string par_ext = pest_utils::lower_cp(par_csv).substr(par_csv.size() - 3, par_csv.size());
				pfm.log_event("processing par stack file " + par_csv);
				if (par_ext.compare("csv") == 0)
				{
					pfm.log_event("loading par stack from csv file: " + par_csv);
					try
					{
						stack_pe.from_csv(par_csv);
					}
					catch (const exception& e)
					{

						throw_constraints_error("error processing par stack: " + string(e.what()));
					}
					catch (...)
					{
						throw_constraints_error(string("error processing par stack"));
					}
				}
				else if ((par_ext.compare("jcb") == 0) || (par_ext.compare("jco") == 0))
				{
					pfm.log_event("loading par stack from binary file+ " + par_csv);
					try
					{
						stack_pe.from_binary(par_csv);
					}
					catch (const exception& e)
					{

						throw_constraints_error("error processing par stack binary file: " + string(e.what()));
					}
					catch (...)
					{
						throw_constraints_error("error processing par stack binary file");
					}
				}
				else
				{
					throw_constraints_error("unrecognized par stack extension, looking for csv, jcb, jco.  found: " + par_ext);
				}
				vector<string> snames = stack_pe.get_var_names();
				set<string> pset(snames.begin(), snames.end());
				vector<string> missing;
				for (auto name : adj_par_names)
					if (pset.find(name) == pset.end())
						missing.push_back(name);
				if (missing.size() > 0)
					throw_constraints_error("par stack missing the following adjustable parameters: ", missing);

				missing.clear();
				for (auto name : dec_var_names)
					if (pset.find(name) == pset.end())
						missing.push_back(name);
				if (missing.size() > 0)
				{
					Eigen::MatrixXd dmat(stack_pe.shape().first, missing.size());
					dmat.setZero();
					stack_pe.extend_cols(dmat,missing);
				}
					

				if ((size_passed) && (stack_pe.shape().first > stack_size))
				{
					vector<int> drop_rows;
					for (int i = stack_size - 1; i < stack_pe.shape().first; i++)
						drop_rows.push_back(i);
					f_rec << "droppping " << drop_rows.size() << "realizations from par stack b/c of ++opt_stack_size req" << endl;
					stack_pe.drop_rows(drop_rows);
				}
			}
			string filename = file_mgr_ptr->get_base_filename() + ".0.par_stack";
			if (pest_scenario.get_pestpp_options().get_ies_save_binary())
			{
				filename = filename + ".jcb";
				stack_pe.to_binary(filename);
			}
			else
			{
				filename = filename + ".csv";
				stack_pe.to_csv(filename);
			}
			f_rec << "saved initial parameter stack to " << filename << endl;
			
			//maybe there was an obs stack passed?
			string obs_csv = pest_scenario.get_pestpp_options().get_opt_obs_stack();
			if (obs_csv.size() == 0)
			{
				stack_oe.reserve(stack_pe.get_real_names(), pest_scenario.get_ctl_ordered_obs_names());
			}
			//if so, load the obs (e.g. constraint) stack
			else
			{
				string obs_ext = pest_utils::lower_cp(obs_csv).substr(obs_csv.size() - 3, obs_csv.size());
				pfm.log_event("processing obs stack file " + obs_csv);
				if (obs_ext.compare("csv") == 0)
				{
					pfm.log_event("loading obs stack from csv file: " + obs_csv);
					try
					{
						stack_oe.from_csv(obs_csv);
					}
					catch (const exception& e)
					{

						throw_constraints_error("error processing obs stack: " + string(e.what()));
					}
					catch (...)
					{
						throw_constraints_error(string("error processing obs stack"));
					}
				}
				else if ((obs_ext.compare("jcb") == 0) || (obs_ext.compare("jco") == 0))
				{
					pfm.log_event("loading obs stack from binary file+ " + obs_csv);
					try
					{
						stack_oe.from_binary(obs_csv);
					}
					catch (const exception& e)
					{

						throw_constraints_error("error processing obs stack binary file: " + string(e.what()));
					}
					catch (...)
					{
						throw_constraints_error("error processing obs stack binary file");
					}
				}
				else
				{
					throw_constraints_error("unrecognized obs stack extension, looking for csv, jcb, jco.  found: " + obs_ext);
				}
				vector<string> snames = stack_oe.get_var_names();
				set<string> pset(snames.begin(), snames.end());
				vector<string> missing;
				for (auto name : ctl_ord_obs_constraint_names)
					if (pset.find(name) == pset.end())
						missing.push_back(name);
				if (missing.size() > 0)
					throw_constraints_error("obs stack missing the following constraints: ", missing);

				if (((size_passed) || (stack_pe.shape().first > 0)) && (stack_oe.shape().first > stack_size))
				{
					vector<int> drop_rows;
					for (int i = stack_size - 1; i < stack_oe.shape().first; i++)
						drop_rows.push_back(i);
					f_rec << "droppping " << drop_rows.size() << "realizations from obs stack b/c of ++opt_stack_size req" << endl;
					stack_oe.drop_rows(drop_rows);
				}
			}
		}
	}
	//otherwise, dont use chance constraints
	else use_chance = false;
}

void Constraints::initial_report()
{
	/*make a rec file report before any iterations have happened*/
	ofstream& f_rec = file_mgr_ptr->rec_ofstream();

	f_rec << endl << "  ---  observation constraints and/or objectives ---  " << endl;
	f_rec << setw(20) << "name" << setw(20) << "sense" << setw(20) << "value" << endl;
	for (auto& name : ctl_ord_obs_constraint_names)
	{
		f_rec << setw(20) << left << name;
		f_rec << setw(20) << constraint_sense_name[name];
		f_rec << setw(20) << constraints_obs.get_rec(name) << endl;
	}

	if (num_pi_constraints() > 0)
	{
		f_rec << endl << "  ---  prior information constraints   ---  " << endl;
		f_rec << setw(20) << "name" << setw(20) << "sense" << setw(20) << "value" << endl;
		for (auto& name : ctl_ord_pi_constraint_names)
		{
			f_rec << setw(20) << left << name;
			f_rec << setw(20) << constraint_sense_name[name];
			f_rec << setw(20) << constraints_pi.get_pi_rec_ptr(name).get_obs_value() << endl;
		}
	}

	if (use_chance)
	{
		f_rec << endl << endl << "  ---  chance constraint/objective FOSM information  ---  " << endl;
		if (use_fosm)
		{
			f_rec << "-->using FOSM-based chance constraints/objectives" << endl;
			f_rec << "-->++opt_risk and corresponding probit function value: " << setw(10) << risk << setw(20) << probit_val << endl;
			f_rec << "-->number of non-zero weight observations for FOSM calcs: " << num_nz_obs() << endl;
		}
		else
		{ 
			f_rec << "-->using stack-based chance constraints/objectives" << endl;
		}
		f_rec << "-->number of adjustable parameters for chance constraint/objective calcs: " << num_adj_pars() << endl;
		
		f_rec << "-->repeat chance constraint/objective calculations every: " << pest_scenario.get_pestpp_options().get_opt_recalc_fosm_every() << " iterations" << endl << endl;
		if (use_fosm)
		{
			f_rec << "   adjustable parameters used in FOSM calculations:" << endl;
			int i = 1;
			for (auto& name : adj_par_names)
			{
				f_rec << setw(15) << name;
				if (i % 6 == 0)
					f_rec << endl;
				i++;
			}
			f_rec << endl;
			if ((num_nz_obs() == 0) && (use_fosm))
			{
				f_rec << endl << endl << "  ---  Note: No nonzero weight observations found." << endl;
				f_rec << "           Prior constraint/objective uncertainty will be used in FOSM-based calculations" << endl;
			}
			else
			{
				f_rec << "  non-zero weight observations used for conditioning in FOSM-based chance constraint/objective calculations: " << endl;
				int i = 0;
				for (auto& name : nz_obs_names)
				{
					f_rec << setw(15) << name;
					if (i % 6 == 0)
						f_rec << endl;
					i++;
				}
				f_rec << endl;
			}
		}
		cout << "...opt_risk = " << risk;
		if (std_weights)
			cout << ", using FOSM-based chance constraints with the weights of non-zero-weighted constraints as standard deviation" << endl;
		else if (use_fosm)
			cout << ", using FOSM-based chance constraints with " << adj_par_names.size() << " adjustable parameters" << endl;
		else
			cout << ", using stack-based chance constraints with " << stack_pe.shape().first << " realizations" << endl;

	}
}


void Constraints::update_chance_offsets()
{
	/* recalculate the chance bits*/
	if ((!use_chance) || (std_weights))
		return;
	prior_const_var.clear();
	post_const_var.clear();
	//if we are using FOSM basd chance constraints, we need to do some new fosm calcs
	if (use_fosm)
	{
		ofstream& f_rec = file_mgr_ptr->rec_ofstream();
		cout << "  ---  calculating FOSM-based chance constraint/objective components  ---  " << endl;
		f_rec << "  ---  calculating FOSM-based chance constraint/objective components  ---  " << endl;

		//the rows of the fosm jacobian include nonzero weight obs (for schur comp)
		//plus the names of the names of constraints, which get treated as forecasts
		vector<string> fosm_row_names(nz_obs_names);
		fosm_row_names.insert(fosm_row_names.end(), ctl_ord_obs_constraint_names.begin(), ctl_ord_obs_constraint_names.end());

		//extract the part of the full jco we need for fosm
		Eigen::SparseMatrix<double> fosm_mat = jco.get_matrix(fosm_row_names, adj_par_names);

		if (fosm_mat.size() == 0)
			throw_constraints_error("FOSM-based chance constraint/objective-to-parameter vectors are all zeros");

		Mat fosm_jco(fosm_row_names, adj_par_names, fosm_mat);

		//create a linear object
		//Logger logger(file_mgr_ptr->get_ofstream("log"), false);
		LinearAnalysis la(fosm_jco, pest_scenario, *file_mgr_ptr, pfm, parcov, &rand_gen);
		la.set_obscov(obscov);

		//set the prior parameter covariance matrix
		la.set_parcov(parcov);

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
		cout << "  ---   done with FOSM-based chance constraint/objective calculations  ---  " << endl << endl;
		f_rec << "  ---   done with FOSM-based chance constraint/objective calculations  ---  " << endl << endl;
	}
	// or if we are using stacks, we just calculate a new token variance
	//this isnt actually used in the shifting since we do a non-parametric shift along the PDF, but 
	//we can still report this variance.
	else
	{
		pair<map<string, double>, map<string, double>> stack_oe_mean_stdev = stack_oe.get_moment_maps();
		for (auto cname : ctl_ord_obs_constraint_names)
		{
			prior_const_var[cname] = stack_oe_mean_stdev.second[cname] * stack_oe_mean_stdev.second[cname];
			post_const_var[cname] = prior_const_var[cname];
		}
	}
}

double Constraints::get_max_constraint_change(Observations& current_obs, Observations& upgrade_obs)
{
	/* work out which constraint has changed the most between the pointer to 
	current sim constraint values and upgrade_obs constraint values*/
	double max_abs_constraint_change = -1.0E+10;
	double max_abs_constraint_val = -1.0E+10;
	double val, diff;
	for (auto& name : ctl_ord_obs_constraint_names)
	{
		val = current_obs.get_rec(name);
		diff = abs(current_obs.get_rec(name) - upgrade_obs[name]);
		max_abs_constraint_change = (diff > max_abs_constraint_change) ? diff : max_abs_constraint_change;
		max_abs_constraint_val = (val > max_abs_constraint_val) ? val : max_abs_constraint_val;
	}
	max_abs_constraint_change /= max(max_abs_constraint_val, 1.0);
	return max_abs_constraint_change;
}

//Observations Constraints::get_chance_shifted_constraints()
//{
//	/* one version of this method that doesnt take any args just use the pointer to the
//	current sim constraint values*/
//	return get_chance_shifted_constraints(*current_constraints_sim_ptr);
//}

ObservationEnsemble Constraints::get_chance_shifted_constraints(ObservationEnsemble& oe)
{
	ObservationEnsemble shifted_oe(oe);//copy
	ofstream& frec = file_mgr_ptr->rec_ofstream();
	if (stack_oe_map.size() == 0)
	{
		
		pfm.log_event("risk-shifting observation population using a single, 'optimal' chance point");
		frec << "risk - shifting observation population using a single, 'optimal' chance point" << endl;

		
		//constraints.presolve_chance_report(iter,);
		Observations sim, sim_shifted;
		Eigen::VectorXd real_vec;
		vector<string> onames = shifted_oe.get_var_names();
		for (int i = 0; i < shifted_oe.shape().first; i++)
		{
			real_vec = shifted_oe.get_real_vector(i);
			sim.update_without_clear(onames, real_vec);
			sim_shifted = get_chance_shifted_constraints(sim);
			shifted_oe.replace(i, sim_shifted);
		}

	}
	else
	{
		pfm.log_event("risk-shifting observation population using 'ALL' decision variable solutions as chance points");
		frec << "risk-shifting observation population using 'ALL' decision variable solutions as chance points" << endl;
		//vector<string> real_names = shifted_oe.get_real_names();
		vector<string> real_names = oe.get_real_names();
		set<string> snames(real_names.begin(), real_names.end());
		map<string, int> real_map = oe.get_real_map();
		Observations sim, sim_shifted;
		Eigen::VectorXd real_vec;
		vector<string> onames = shifted_oe.get_var_names();
		for (auto& real_info : stack_oe_map)
		{
			if (snames.find(real_info.first) == snames.end())
			{
				throw_constraints_error("dec var population realization name '" + real_info.first + "' not found in observation population realization names");
			}
			//overwrite the class stack_oe attribute with the realization stack
			stack_oe = real_info.second;
			real_vec = shifted_oe.get_real_vector(real_info.first);
			sim.update_without_clear(onames, real_vec);
			//this call uses the class stack_oe attribute;
			sim_shifted = get_chance_shifted_constraints(sim);
			shifted_oe.replace(real_map[real_info.first], sim_shifted);
		}
	}
	return shifted_oe;
}

Observations Constraints::get_chance_shifted_constraints(Observations& current_obs)
{
	/* get the simulated constraint values with the chance shift applied*/
	ofstream& f_rec = file_mgr_ptr->rec_ofstream();
	prior_constraint_offset.clear();
	prior_constraint_stdev.clear();
	post_constraint_offset.clear();
	post_constraint_stdev.clear();
	double new_constraint_val, old_constraint_val, required_val;
	double pr_offset, pt_offset;
	Observations iter_fosm;

	if (use_fosm)
	{
		for (auto& name : ctl_ord_obs_constraint_names)
		{
			prior_constraint_stdev[name] = sqrt(prior_const_var[name]);
			post_constraint_stdev[name] = sqrt(post_const_var[name]);
			//the offset (shift) is just the stdev * the risk-based probit value
			pr_offset = probit_val * prior_constraint_stdev[name];
			pt_offset = probit_val * post_constraint_stdev[name];
			old_constraint_val = current_obs.get_rec(name);
			required_val = constraints_obs[name];

			//if less_than constraint, then add to the sim value, to move positive
			// WRT the required constraint value
			if (constraint_sense_map[name] == ConstraintSense::less_than)
			{

				new_constraint_val = old_constraint_val + pt_offset;
				post_constraint_offset[name] = pt_offset;
				prior_constraint_offset[name] = pr_offset;
			}
			//if greater_than constraint, the substract from the sim value to move
			//negative WRT the required constraint value
			else if (constraint_sense_map[name] == ConstraintSense::greater_than)
			{
				new_constraint_val = old_constraint_val - pt_offset;
				post_constraint_offset[name] = -pt_offset;
				prior_constraint_offset[name] = -pr_offset;
			}
			//if its not an ineq constraint, then there is no direction to shift
			else
				new_constraint_val = current_obs.get_rec(name);
			iter_fosm.insert(name, new_constraint_val);
		}
	}
	else
	{
		//work out which realization index corresponds to the risk value
		if (stack_oe.shape().first < 3)
			throw_constraints_error("too few (<3) stack members, cannot continue with stack-based chance constraints/objectives");
		int cur_num_reals = stack_oe.shape().first;
		//get the inde (which realization number) represents the risk value according to constraint sense
		//
		//the "less than" realization index is just the risk value times the number of reals
		int lt_idx = int(risk * cur_num_reals);
		//the greater than index is the opposite direction
		int gt_idx = int((1.0-risk) * cur_num_reals);
		//the equality contraint risk index
		int eq_idx = int(0.5 * cur_num_reals);
		//get the mean-centered anomalies - we want to subtract off the mean 
		//in case these stack values are being re-used from a previous iteration
		Eigen::MatrixXd anom = stack_oe.get_eigen_anomalies();
		//get the map of realization name to index location in the stack
		stack_oe.update_var_map();
		map<string, int> var_map = stack_oe.get_var_map();
		//get the mean and stdev summary containters, summarized by observation (e.g. constraint) name
		pair<map<string, double>, map<string, double>> mm = stack_oe.get_moment_maps();
		for (auto& name : ctl_ord_obs_constraint_names)
		{
			old_constraint_val = current_obs.get_rec(name);
			//the value that must be statified (from the control file)
			required_val = constraints_obs[name];
			// the realized values of this stack are the anomalies added to the
			//current constraint value - this assumes the current value 
			//is the mean of the stack distribution
			Eigen::VectorXd cvec = anom.col(var_map[name]).array() + old_constraint_val;
			//now sort the anomolies + current (mean) value vector
			sort(cvec.data(), cvec.data() + cvec.size());
			//set the stdev container info - this isnt used in 
			//calculations but gets reported
			prior_constraint_stdev[name] = mm.second[name];
			post_constraint_stdev[name] = mm.second[name];
			
			if (prior_constraint_stdev[name] == 0.0)
			{
				throw_constraints_error("model-based constraint '" + name + "' has empirical (stack) standard deviation of 0.0, something is wrong");
			}

			// if this is a "less than" constraint
			if (constraint_sense_map[name] == ConstraintSense::less_than)
			{
				//the posterior shifted constraint value is the stack value at the less than index location minus the 
				//current constraint value
				pt_offset = cvec[lt_idx] - old_constraint_val;
				new_constraint_val = cvec[lt_idx];
				post_constraint_offset[name] = pt_offset;
				prior_constraint_offset[name] = pt_offset;
				
			}
			
			else if (constraint_sense_map[name] == ConstraintSense::greater_than)
			{
				//the posterior shifted constraint value is the stack value at the greater than index location minus the 
				//current constraint value
				pt_offset = cvec[gt_idx] - old_constraint_val;
				new_constraint_val = cvec[gt_idx];
				post_constraint_offset[name] = -pt_offset;
				prior_constraint_offset[name] = -pt_offset;

			}
			else
				new_constraint_val = old_constraint_val;
			iter_fosm.insert(name, new_constraint_val);
		}
	}
	vector<string> names = iter_fosm.get_keys();
	Observations constraints_chance(current_obs);
	constraints_chance.update_without_clear(names, iter_fosm.get_data_vec(names));
	
	return constraints_chance;
}


vector<double> Constraints::get_constraint_residual_vec(Observations& sim)
{
	/* get the current distance to the the constraint bound - the edge of feasible*/
	vector<double> residuals_vec;
	residuals_vec.resize(num_constraints(), 0.0);

	Observations::const_iterator found_obs;
	Observations::const_iterator not_found_obs = sim.end();
	PriorInformation::const_iterator found_prior_info;

	int i = 0;
	for (vector<string>::iterator b = ctl_ord_obs_constraint_names.begin(), e = ctl_ord_obs_constraint_names.end(); b != e; ++b, ++i)
	{
		found_obs = sim.find(*b);
		if (found_obs != not_found_obs)
		{
			residuals_vec[i] = constraints_obs.get_rec(*b) - ((*found_obs).second);
		}
	}
	return residuals_vec;
}

pair<vector<double>,vector<double>> Constraints::get_constraint_bound_vectors(Parameters& current_pars, Observations& current_obs)
{
	/* get the upper and lower bound constraint vectors. For less than constraints, the lower bound is 
	set to double max, for greater than constraints, the upper bound is set to double max.
	These are needed for the simplex solve*/
	vector<double> residuals;
	if (use_chance)
	{
		Observations current_constraints_chance = get_chance_shifted_constraints(current_obs);
		residuals = get_constraint_residual_vec(current_constraints_chance);
	}
	else
	{
		residuals = get_constraint_residual_vec(current_obs);
	}

	vector<double> constraint_ub, constraint_lb;
	for (int i = 0; i < num_obs_constraints(); ++i)
	{
		string name = ctl_ord_obs_constraint_names[i];
		if (constraint_sense_map[name] == ConstraintSense::less_than)
			constraint_ub.push_back(residuals[i]);
		else
			constraint_ub.push_back(dbl_max);

		if (constraint_sense_map[name] == ConstraintSense::greater_than)
			constraint_lb.push_back(residuals[i]);
		else
			constraint_lb.push_back(-dbl_max);
		if (constraint_sense_map[name] == ConstraintSense::equal_to)
		{
			constraint_ub.push_back(residuals[i]);
			constraint_lb.push_back(residuals[i]);
		}
	}


	int noc = num_obs_constraints();
	for (int i = 0; i < num_pi_constraints(); ++i)
	{
		string name = ctl_ord_pi_constraint_names[i];
		double residual = -constraints_pi.get_pi_rec_ptr(name).calc_residual(current_pars);
		if (constraint_sense_map[name] == ConstraintSense::less_than)
			constraint_ub.push_back(residual);
		else
			constraint_ub.push_back(dbl_max);

		if (constraint_sense_map[name] == ConstraintSense::greater_than)
			constraint_lb.push_back(residual);
		else
			constraint_lb.push_back(-dbl_max);
		if (constraint_sense_map[name] == ConstraintSense::equal_to)
		{
			constraint_ub.push_back(residual);
			constraint_lb.push_back(residual);
		}
	}
	current_bounds = pair<vector<double>, vector<double>>(constraint_lb, constraint_ub);
	return current_bounds;
}


pair<Constraints::ConstraintSense, string> Constraints::get_sense_from_group_name(const string& name)
{
	/* work out the constraint sense from the obs group name*/
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

double  Constraints::ErfInv2(double x)
{
	/* the inverse error function, needed for the FOSM-based chance constraints*/
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
	/* the probit function estimate - needed for the fosm-basd chance constraints*/
	double output = sqrt(2.0) * ErfInv2((2.0 * risk) - 1.0);
	return output;
}

void Constraints::presolve_report(int iter, Parameters& current_pars, Observations& current_obs)
{
	/* this is a report to the rec file for the status of the constraints before solving the current iteration*/
	
	ofstream& f_rec = file_mgr_ptr->rec_ofstream();
	vector<double> residuals;
	f_rec << endl << "  observation constraint/objective information at start of iteration " << iter << endl;
	f_rec << setw(20) << left << "name" << right << setw(15) << "sense" << setw(15) << "required" << setw(15) << "sim value";
	f_rec << setw(15) << "residual" << setw(15) << "lower bound" << setw(15) << "upper bound" << endl;
	
	//if we are using chances, then we need to account for that here
	Observations current;
	if (use_chance)
	{
		current = get_chance_shifted_constraints(current_obs);	
	}
	//otherwise, just use the pointer to the current constraint info
	else
	{
		current = current_obs;
	}

	residuals = get_constraint_residual_vec(current);
	for (int i = 0; i < num_obs_constraints(); ++i)
	{
		string name = ctl_ord_obs_constraint_names[i];
		f_rec << setw(20) << left << name;
		f_rec << setw(15) << right << constraint_sense_name[name];
		f_rec << setw(15) << constraints_obs.get_rec(name);	
		f_rec << setw(15) << current.get_rec(name);
		f_rec << setw(15) << residuals[i];
		f_rec << setw(15) << current_bounds.first[i];
		f_rec << setw(15) << current_bounds.second[i] << endl;
	}
	if (num_pi_constraints() > 0)
	{

		//report prior information constraints
		f_rec << endl << "  prior information constraint information at start of iteration " << iter << endl;
		f_rec << setw(20) << left << "name" << right << setw(15) << "sense" << setw(15) << "required" << setw(15) << "sim value";
		f_rec << setw(15) << "residual" << setw(15) << "lower bound" << setw(15) << "upper bound" << endl;
		for (int i = 0; i < num_pi_constraints(); ++i)
		{
			string name = ctl_ord_pi_constraint_names[i];
			PriorInformationRec pi_rec = constraints_pi.get_pi_rec_ptr(name);
			f_rec << setw(20) << left << name;
			f_rec << setw(15) << right << constraint_sense_name[name];
			f_rec << setw(15) << pi_rec.get_obs_value();
			f_rec << setw(15) << pi_rec.calc_sim_and_resid(current_pars).first;
			f_rec << setw(15) << pi_rec.calc_residual(current_pars);
			f_rec << setw(15) << current_bounds.first[num_obs_constraints() + i];
			f_rec << setw(15) << current_bounds.second[num_obs_constraints() + i] << endl;

		}
	}

	if (use_chance)
		presolve_chance_report(iter,current_obs);	 
	
	return;
}

void Constraints::write_res_files(Observations& constraints, Parameters& pars_and_dec_vars, string tag, int iter)
{
	/* one form of this method that doesnt require advanced knowledge of whether or not chances are being used.
	write a pest-style residuals file, and, if using chances, write a sep file that includes the chance offset info*/
	write_res_file(constraints, pars_and_dec_vars, tag, iter, false);
	if (use_chance)
		write_res_file(constraints, pars_and_dec_vars, tag, iter, true);
}

void Constraints::write_res_file(Observations& constraints, Parameters& pars_and_dec_vars, string tag, int iter, bool include_chance)
{
	/* write a pest-style residuals file.  If chances are in use, then add "chance" to the file name
	and write chance-shifted values for the "modelled" column*/
	stringstream ss;
	Observations temp = constraints;
	if (include_chance)
	{
		ss << iter << "." << tag << "+chance.rei";
		for (auto& o : ctl_ord_obs_constraint_names)
			temp[o] = temp[o] + post_constraint_offset[o];
	}
	else
	{

		ss << iter << "." << tag << ".rei";
	}
	of_wr.write_opt_constraint_rei(file_mgr_ptr->open_ofile_ext(ss.str()), iter, pars_and_dec_vars,
		pest_scenario.get_ctl_observations(), temp);
	file_mgr_ptr->close_file(ss.str());

}

void Constraints::presolve_chance_report(int iter, Observations& current_obs)
{
	/* write chance info to the rec file before undertaking the current iteration process*/
	if (!use_chance)
		return;
	
	ofstream& f_rec = file_mgr_ptr->rec_ofstream();
	//vector<double> residuals = get_constraint_residual_vec();
	f_rec << endl << "  FOSM-based chance constraint/objective information at start of iteration " << iter << endl;
	f_rec << setw(20) << left << "name" << right << setw(10) << "sense" << setw(15) << "required" << setw(15) << "sim value";
	f_rec << setw(15) << "prior stdev" << setw(15) << "post stdev" << setw(15) << "offset";
	f_rec << setw(15) << "new sim value" << endl;
	vector<string> out_of_bounds;
	Observations current_constraints_chance = get_chance_shifted_constraints(current_obs);
	for (int i = 0; i < num_obs_constraints(); ++i)
	{
		string name = ctl_ord_obs_constraint_names[i];
		f_rec << setw(20) << left << name;
		f_rec << setw(10) << right << constraint_sense_name[name];
		f_rec << setw(15) << constraints_obs[name];
		f_rec << setw(15) << current_obs.get_rec(name);
		f_rec << setw(15) << prior_constraint_stdev[name];
		f_rec << setw(15) << post_constraint_stdev[name];
		f_rec << setw(15) << post_constraint_offset[name];
		f_rec << setw(15) << current_constraints_chance[name] << endl;
	}
	f_rec << "  note: 'offset' is the value added to the simulated constraint/objective value to account" << endl;
	f_rec << "        for the uncertainty in the constraint/objective value arsing from uncertainty in the " << endl;
	f_rec << "        adjustable parameters identified in the control file." << endl << endl;
	if (!use_fosm)
	{
		f_rec << "  note: the above standard deviations are empirical estimates from the stack" << endl;
	}
	return;
	
}


bool Constraints::should_update_chance(int iter)
{
	/* this is a total hack - it tries to determine if it is time to update the chance
	information, like should we queue up some JCO pertubation runs for FOSM chances or 
	should we queue up some stack runs for stack-based chances*/
	if (!use_chance)
		return false;
	if (std_weights)
		return false;
	if (iter == 1)
	{
		if (use_fosm)
			return true;
		if (pest_scenario.get_pestpp_options().get_opt_obs_stack().size() > 0)
			return false;
		else
			return true;
	}
	else if ((iter + 1) % pest_scenario.get_pestpp_options().get_opt_recalc_fosm_every() == 0)
		return true;
	return false;
}


void Constraints::postsolve_obs_constraints_report(Observations& old_obs, Observations& new_obs, string tag, int iter, 
									map<string,string> status_map, map<string,double> price_map)
{
	/* write out constraint info after the current iteration solution process is over to the rec file.  
	
	This really only applies to the simplex solution because we can, with the linear assumption, 
	project what the resulting constraint values should be*/

	ofstream& f_rec = file_mgr_ptr->rec_ofstream();
	f_rec << endl << endl << "     " << tag << " constraint/objective information at end of iteration " << iter << endl << endl;
	f_rec << setw(20) << left << "name" << right << setw(15) << "sense" << setw(15) << "required" << setw(25);
	if (status_map.size() > 0)
		f_rec << "simplex status";
	if (price_map.size() > 0)
		f_rec << setw(15) << "price";
	if (use_chance)
		f_rec << setw(15) << "fosm offset";
	f_rec << setw(15) << "current" << setw(15) << "residual";
	f_rec << setw(15) << "new" << setw(15) << "residual" << endl;

	vector<double> cur_residuals = get_constraint_residual_vec(old_obs);
	vector<double> new_residuals;
	if (use_chance)
	{
		Observations constraints_shift = get_chance_shifted_constraints(new_obs);
		new_residuals = get_constraint_residual_vec(constraints_shift);
	}
	else
		new_residuals = get_constraint_residual_vec(new_obs);
	double sim_val;
	for (int i = 0; i < num_obs_constraints(); ++i)
	{
		string name = ctl_ord_obs_constraint_names[i];
		sim_val = new_obs[name];
		f_rec << setw(20) << left << name;
		f_rec << setw(15) << right << constraint_sense_name[name];
		f_rec << setw(15) << constraints_obs.get_rec(name);
		if (status_map.size() > 0)
			f_rec << setw(25) << status_map[name];
		if (price_map.size() > 0)
			f_rec << setw(15) << price_map[name];
		if (use_chance)
		{
			double offset = post_constraint_offset[name];
			sim_val += offset;
			f_rec << setw(15) << offset;
			f_rec << setw(15) << old_obs.get_rec(name);
			f_rec << setw(15) << cur_residuals[i];
			f_rec << setw(15) << sim_val;
			f_rec << setw(15) << new_residuals[i] << endl;
		}
		else
		{
			f_rec << setw(15) << old_obs.get_rec(name);
			f_rec << setw(15) << cur_residuals[i];
			f_rec << setw(15) << sim_val;
			f_rec << setw(15) << new_residuals[i] << endl;

		}

	}
}

void Constraints::postsolve_pi_constraints_report(Parameters& old_pars, Parameters& new_pars, int iter, map<string,string> status_map, map<string,double> price_map)
{
	
	/*report prior information constraints stats to the rec file after the current iteration solution process.  

	This really only applies to the simplex solution because we can, with the linear assumption,
	project what the resulting constraint values should be
	*/
	ofstream &f_rec = file_mgr_ptr->rec_ofstream();
	if (num_pi_constraints() > 0)
	{
		f_rec << endl << endl << "     prior information constraint information at end of iteration " << iter << endl << endl;
		f_rec << setw(20) << left << "name" << right << setw(15) << "sense" << setw(15) << "required" << setw(25);
		if (status_map.size() > 0)
			f_rec << "simplex status";
		if (price_map.size() > 0)
			f_rec << setw(15);
		f_rec << "price" << setw(15) << "current" << setw(15) << "residual";
		f_rec << setw(15) << "new" << setw(15) << "residual" << endl;
		int i = 0;
		for (auto &name : ctl_ord_pi_constraint_names)
		{
			PriorInformationRec pi_rec = constraints_pi.get_pi_rec_ptr(name);
			pair<double,double> cur_sim_resid = pi_rec.calc_sim_and_resid(old_pars);
			pair<double,double> new_sim_resid = pi_rec.calc_sim_and_resid(new_pars);
			f_rec << setw(20) << left << name;
			f_rec << setw(15) << right << constraint_sense_name[name];
			f_rec << setw(15) << pi_rec.get_obs_value();
			if (status_map.size() > 0)
				f_rec << setw(25) << status_map[name];
			if (price_map.size() > 0)
				f_rec << setw(15) << price_map[name];
			f_rec << setw(15) << cur_sim_resid.first;
			f_rec << setw(15) << cur_sim_resid.second;
			f_rec << setw(15) << new_sim_resid.first;
			f_rec << setw(15) << new_sim_resid.second << endl;
			i++;
		}
	}
	
	return;
}

pair<vector<int>,ObservationEnsemble> Constraints::process_stack_runs(string real_name, int iter, map<int, int> _stack_pe_run_map, 
	RunManagerAbstract* run_mgr_ptr, bool drop_fails)
{
	ObservationEnsemble _stack_oe(stack_oe);//copy
	vector<int> failed_runs = _stack_oe.update_from_runs(_stack_pe_run_map, run_mgr_ptr);
	stringstream ss;
	if (failed_runs.size() > 0)
	{
		ss.str("");
		ss << "WARNING: " << failed_runs.size() << " stack runs failed for realization " << real_name;
		pfm.log_event(ss.str());
		cout << ss.str() << endl;
	}
		
	return pair<vector<int>,ObservationEnsemble>(failed_runs,_stack_oe);
}

void Constraints::save_oe_stack(int iter, string real_name, ObservationEnsemble& _stack_oe)
{
	stringstream ss;
	ss.str("");
	if (real_name.size() > 0)
		ss << file_mgr_ptr->get_base_filename() << "." << iter << "." << real_name << ".obs_stack";
	else
		ss << file_mgr_ptr->get_base_filename() << "." << iter << ".obs_stack";

	if (pest_scenario.get_pestpp_options().get_ies_save_binary())
	{
		ss << ".jcb";
		_stack_oe.to_binary(ss.str());
	}
	else
	{
		ss << ".csv";
		_stack_oe.to_csv(ss.str());
	}
	ss.str("");
	if (real_name.size() > 0)
		pfm.log_event("saved realization '" + real_name + "' stack_oe to " + ss.str());
	else
		pfm.log_event("saved stack_oe to " + ss.str());
}


void Constraints::save_pe_stack(int iter, string real_name, ParameterEnsemble& _stack_pe)
{
	stringstream ss;
	ss.str("");
	if (real_name.size() > 0)
		ss << file_mgr_ptr->get_base_filename() << "." << iter << "." << real_name << ".par_stack";
	else
		ss << file_mgr_ptr->get_base_filename() << "." << iter << ".par_stack";

	if (pest_scenario.get_pestpp_options().get_ies_save_binary())
	{
		ss << ".jcb";
		_stack_pe.to_binary(ss.str());
	}
	else
	{
		ss << ".csv";
		_stack_pe.to_csv(ss.str());
	}
	ss.str("");
	if (real_name.size() > 0)
		pfm.log_event("saved realization '" + real_name + "' stack_pe to " + ss.str());
	else
		pfm.log_event("saved stack_pe to " + ss.str());
}


void Constraints::process_stack_runs(RunManagerAbstract* run_mgr_ptr, int iter)
{
	pair<vector<int>, ObservationEnsemble> stack_info;
	if (population_stack_pe_run_map.size() > 0)
	{
		stack_oe_map.clear();
		//work out what var names for the stack_oe we dont need so we can drop them!
		vector<string> drop_names;
		pair<vector<int>, ObservationEnsemble> stack_info;
		set<string> s_obs_constraint_names(ctl_ord_obs_constraint_names.begin(), ctl_ord_obs_constraint_names.end());
		for (auto name : stack_oe.get_var_names())
		{
			if (s_obs_constraint_names.find(name) == s_obs_constraint_names.end())
				drop_names.push_back(name);
		}
		for (auto& real_info : population_stack_pe_run_map)
		{
			pfm.log_event("processing stack runs for realization " + real_info.first);
			stack_info = process_stack_runs(real_info.first, iter, real_info.second, run_mgr_ptr, false);
			save_oe_stack(iter, real_info.first, stack_info.second);
			stack_info.second.drop_cols(drop_names);
			stack_oe_map[real_info.first] = stack_info.second;
			
		}
	}
	else
	{
		pfm.log_event("processing stack runs");
		stack_info = process_stack_runs("", iter, stack_pe_run_map, run_mgr_ptr, true);
		if (stack_info.first.size() > 0)
		{
			file_mgr_ptr->rec_ofstream() << "dropping " << stack_info.first.size() << " failed realizations for par stack";
			pfm.log_event("dropping failed realizations from par stack");
			stack_pe.drop_rows(stack_info.first);
		}
		save_oe_stack(iter, "", stack_info.second);
		save_pe_stack(iter, "", stack_pe);
		stack_oe = stack_info.second;
	}
}

void Constraints::process_runs(RunManagerAbstract* run_mgr_ptr,int iter)
{
	/* using the passed in run mgr pointer, process any runs that were queued up for chance-based 
	calculations. handles FOSM, single stack and population nested stacks
	
	*/
	stringstream ss;
	ofstream& f_rec = file_mgr_ptr->rec_ofstream();
	if (!use_chance)
		return;
	//if we are using fosm chances, then process some jco runs
	if (use_fosm)
	{
		if (jco.get_par_run_map().size() == 0)
			return;
		pfm.log_event("reading FOSM-based parameter pertubation runs into JCO");
		ParamTransformSeq par_trans = pest_scenario.get_base_par_tran_seq();
		bool success = jco.process_runs(par_trans, pest_scenario.get_base_group_info(), *run_mgr_ptr,
			*null_prior, false, false);
		if (!success)
			throw_constraints_error("error processing FOSM JCO matrix runs ", jco.get_failed_parameter_names());
		ss.str("");
		ss << file_mgr_ptr->get_base_filename() << "." << iter << ".fosm.jcb";
		jco.save(ss.str());
	}
	//otherwise, process some stack runs
	else
	{

		process_stack_runs(run_mgr_ptr,iter);

	}

}

void Constraints::add_runs(int iter, Parameters& current_pars, Observations& current_obs, RunManagerAbstract* run_mgr_ptr)
{
	/* using the passed run mgr pointer, queue up chance runs
	
	*/
	if (!use_chance)
		return;
	//for fosm, we need to queue up parameter pertubation runs, centered at the current dec var values
	if (use_fosm)
	{
		pfm.log_event("building FOSM-based parameter pertubation runs");
		ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();
		set<string> out_of_bounds;
		//the  second-to-last true arg is to run the parvals listed in _current_pars_and_dec_vars as the "base" run
		//this is to make sure we get good pertubations at the current point in par/dec var space.
		//but this means we might have to run the base run twice at least for the SLP process.
		//the way around this is to just add the runs all at once wihin the SLP add run
		//bit using Constraints::get_fosm_par_names()
		//the last false says not to reinitialize the run mgr since the calling process may have also
		//added runs
		bool success = jco.build_runs(current_pars, current_obs, adj_par_names, pts,
			pest_scenario.get_base_group_info(), pest_scenario.get_ctl_parameter_info(),
			*run_mgr_ptr, out_of_bounds, false, true, false);
		if (!success)
		{
			const set<string> failed = jco.get_failed_parameter_names();
			throw_constraints_error("failed to calc derviatives for the following FOSM parameters: ", failed);
		}
		cout << "...adding " << jco.get_par_run_map().size() << " model runs for FOSM-based chance constraints" << endl;

	}
	//for stacks, we need to queue up the stack realizations, but replace the dec var entries in each realization
	//with the current dec var values.
	else
	{
		
		//update stack_pe parameter values for decision variables using current_pars_and_dec_vars_ptr
		int num_reals = stack_pe.shape().first;
		map<string, int> var_map = stack_pe.get_var_map();
		for (auto dname : dec_var_names)
		{
			Eigen::VectorXd dvec(num_reals);
			dvec.setConstant(current_pars.get_rec(dname));
			stack_pe.replace_col(dname, dvec);
		}
		pfm.log_event("building stack-based parameter runs");
		cout << "...adding " << stack_pe.shape().first << " model runs for stack-based chance constraints" << endl;
		stack_pe_run_map.clear();
		stack_pe_run_map = stack_pe.add_runs(run_mgr_ptr);
	}
}

void Constraints::add_runs(int iter, ParameterEnsemble& current_pe, Observations& current_obs, RunManagerAbstract* run_mgr_ptr)
{
	if (!use_chance)
		return;
	if (use_fosm)
	{
		throw_constraints_error("add_runs() error: FOSM-based chance constraints not supported for ensemble/population-based algorithms");
	}
	population_stack_pe_run_map.clear();
	pfm.log_event("queuing up nested-sets of chance runs for decision-variable ensemble ");
	
	Eigen::VectorXd par_vec;
	vector<string> par_names = current_pe.get_var_names();
	Parameters real_pars = pest_scenario.get_ctl_parameters();
	pest_scenario.get_base_par_tran_seq().ctl2numeric_ip(real_pars);
	current_pe.transform_ip(ParameterEnsemble::transStatus::NUM);
	ofstream& f_rec = file_mgr_ptr->rec_ofstream();
	f_rec << "queuing up nested-sets of chance runs ('opt_chance_points' == 'all') for decision-variable ensemble " << endl;
	cout << "queuing up nested - sets of chance runs for decision-variable ensemble " << endl;
	for (auto real_info : current_pe.get_real_map())
	{

		par_vec = current_pe.get_real_vector(real_info.first);
		real_pars.update_without_clear(par_names, par_vec);
		add_runs(iter, real_pars, current_obs, run_mgr_ptr);
		//relying on class attribute being set in add_runs():
		population_stack_pe_run_map[real_info.first] = stack_pe_run_map;
		//relying on class attribute being set in add_runs():
		save_pe_stack(iter, real_info.first, stack_pe);
		f_rec << "...added " << stack_pe.shape().first << " runs for decision variable solution '" << real_info.first << "'" << endl;
	}
}

vector<string> Constraints::get_fosm_par_names()
{
	/* get a vector of the fosm-based "parameter" names - adjustable, non-dec-var parameters in the control file*/
	if (use_fosm)
		return adj_par_names;
	else
		return vector<string>();
}

map<string, double> Constraints::get_unsatified_pi_constraints(Parameters& par_and_dec_vars, double tol)
{
	/* get a map of name, distance for each of the prior info constraints that are not satisfied in the par_and_dec_vars container.
	tol is a percent-based tolerance to accont for constraints that are very near their required (rhs) value

	*/
	double sim_val, obs_val, scaled_diff;
	map<string, double> unsatisfied;
	
	for (auto& name : ctl_ord_pi_constraint_names)
	{
		PriorInformationRec pi_rec = constraints_pi.get_pi_rec_ptr(name);
		//pair<double, double> cur_sim_resid = pi_rec.calc_residual_and_sim_val(*current_pars_and_dec_vars_ptr);
		pair<double, double> new_sim_resid = pi_rec.calc_sim_and_resid(par_and_dec_vars);
		//check for invalid pi constraints
		sim_val = new_sim_resid.first;
		obs_val = pi_rec.get_obs_value();
		scaled_diff = abs((obs_val - sim_val) / obs_val);
		if ((constraint_sense_map[name] == ConstraintSense::less_than) && (sim_val > obs_val) && (scaled_diff > tol))
			unsatisfied[name] = sim_val - obs_val;
		else if ((constraint_sense_map[name] == ConstraintSense::greater_than) && (sim_val < obs_val) && (scaled_diff > tol))
			unsatisfied[name] = obs_val - sim_val;
		else if ((constraint_sense_map[name] == ConstraintSense::equal_to) && (sim_val != obs_val) && (scaled_diff > tol))
			unsatisfied[name] = abs(sim_val - obs_val);
	}
	return unsatisfied;
}

map<string, double> Constraints::get_unsatified_obs_constraints(Observations& constraints_sim, double tol, bool do_shift)
{
	/* get a map of name, distance for each of the obs-based (e.g. model-based) constraints that are not satisfied in the constraint_obs container.
	tol is a percent-based tolerance to accont for constraints that are very near their required (rhs) value

	*/
	double sim_val, obs_val, scaled_diff;
	map<string, double> unsatisfied;
	Observations _constraints_sim(constraints_sim);
	if ((do_shift) && (use_chance))
		_constraints_sim = get_chance_shifted_constraints(constraints_sim);
	for (int i = 0; i < num_obs_constraints(); ++i)
	{
		string name = ctl_ord_obs_constraint_names[i];
		sim_val = constraints_sim[name];
		if ((do_shift) && (use_chance))
		{
			double offset = post_constraint_offset[name];
			sim_val += offset;
		}
		obs_val = constraints_obs[name];
		scaled_diff = abs((obs_val - sim_val) / obs_val);
		//check for invalid obs constraints (e.g. satified)
		if ((constraint_sense_map[name] == ConstraintSense::less_than) && (sim_val > obs_val) && (scaled_diff > tol))
			unsatisfied[name] = sim_val - obs_val;
		else if ((constraint_sense_map[name] == ConstraintSense::greater_than) && (sim_val < obs_val) && (scaled_diff > tol))
			unsatisfied[name] = obs_val - sim_val;
		else if ((constraint_sense_map[name] == ConstraintSense::equal_to) && (sim_val != obs_val) && (scaled_diff > tol))
			unsatisfied[name] = abs(sim_val - obs_val);
	}
	return unsatisfied;

}


int Constraints::get_num_nz_pi_constraint_elements()
{
	/* work out how many elements are in the PI constraints.  This is needed for dimensioning the 
	simplex response matrix*/
	int num = 0;

	for (auto &pi_name : ctl_ord_pi_constraint_names)
	{
		num += constraints_pi.get_pi_rec_ptr(pi_name).get_atom_factors().size();
	}
	return num;
}



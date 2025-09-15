
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include <iterator>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>

#include "SVDPackage.h"
#include "Pest.h"
#include "utilities.h"
#include "covariance.h"
#include "FileManager.h"
#include "constraints.h"
#include "linear_analysis.h"
#include "eigen_tools.h"
#include "MOEA.h"

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
    stringstream ss;
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
        ss.str("");
        ss << "...objective function defined by observation '" << obj_func_str << "'" << endl;
        cout << ss.str();
        f_rec << ss.str();
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
			f_rec << " note: no ++opt_objective_function-->forming a generic objective function (1.0 coef for each decision var)" << endl;
			cout << " note: no ++opt_objective_function-->forming a generic objective function (1.0 coef for each decision var)" << endl;
            for (auto& name : dv_names)
				obj_func_coef_map[name] = 1.0;
		}

		//or if it is a prior info equation
		else if (pest_scenario.get_prior_info().find(obj_func_str) != pest_scenario.get_prior_info().end())
		{
            ss.str("");
            ss << "...objective function defined by prior information equation '" << obj_func_str << "'" << endl;
            cout << ss.str();
            f_rec << ss.str();
			obj_func_coef_map = pest_scenario.get_prior_info().get_pi_rec(obj_func_str).get_atom_factors();
			//throw_sequentialLP_error("prior-information-based objective function not implemented");
		}
		else
		{
			//check if this obj_str is a filename
            obj_func_str = pest_scenario.get_pestpp_options().get_org_opt_obj_func();
            ss.str("");
            ss << "...objective function defined by 2-column external file '" << obj_func_str << "'" << endl;
            cout << ss.str();
            f_rec << ss.str();
            if (!pest_utils::check_exist_in(obj_func_str))
            {
                throw_optobjfunc_error("unable to open objective function file '"+obj_func_str+"' for reading");
            }
			/*ifstream if_obj(obj_func_str);
			if (!if_obj.good())
				throw_optobjfunc_error("unrecognized ++opt_objective_function arg: " + obj_func_str);
			else*/
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
	:pest_scenario(_pest_scenario), file_mgr_ptr(_file_mgr_ptr), of_wr(_of_wr), pfm(_pfm), 
	jco(*_file_mgr_ptr, _of_wr)
{

	use_fosm = false;
	std_weights = false;
	probit_val = 0.0;
	stack_runs_processed = false;
	stack_pe.set_pest_scenario_ptr(&_pest_scenario);
	nested_pe.set_pest_scenario_ptr(&_pest_scenario);

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
	//check for a stack size arg - if positive, then use stacks for chances
	stack_size = pest_scenario.get_pestpp_options().get_opt_stack_size();
	//an existing parameter stack for chances
	string par_stack_name = pest_scenario.get_pestpp_options().get_opt_par_stack();
	//maybe even an existing observations (e.g. constraints) stack!
	string obs_stack_name = pest_scenario.get_pestpp_options().get_opt_obs_stack();
	//by default, we want to use fosm for chances, but it any
	//of those stack options were passed, then use stacks instead
	use_fosm = true;
	std_weights = pest_scenario.get_pestpp_options().get_opt_std_weights();
	if ((!std_weights) && ((stack_size > 0) || (par_stack_name.size() > 0) || (obs_stack_name.size() > 0)))
		use_fosm = false;
	//initialize the stack containers (ensemble class instances)
	stack_pe.set_pest_scenario(&pest_scenario);
	stack_pe.set_rand_gen(&rand_gen);
	stack_oe.set_pest_scenario_ptr(&pest_scenario);
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
		f_rec << "...augmenting prior information constraints with decision variable bounds" << endl;
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
			group = pinfo->get_pi_rec(pi_name).get_group();
			if (find(start, end, group) != end)
			{
				ctl_ord_pi_constraint_names.push_back(pi_name);
				pi_rec = pinfo->get_pi_rec(pi_name);
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
				pi_rec = constraints_pi.get_pi_rec(pi_name);
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
			ss << " the following prior information constraints/objectives reference parameters that are not treated as decision variables:" << endl;
			for (auto& missing_pi : missing_map)
			{
				ss << missing_pi.first << ": ";
				for (auto& par_name : missing_pi.second)
					ss << par_name << ",";
			}
			throw_constraints_error("errors in prior information constraints/objectives:" + ss.str());
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
		string group = pest_scenario.get_prior_info_ptr()->get_pi_rec(name).get_group();
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
		throw_constraints_error("the following prior info constraints/objectives do not have a correct group name prefix {'l_','less','g_','greater','e_','equal'}: ", problem_constraints);
	}

	set<string> dec_set(ctl_ord_dec_var_names.begin(), ctl_ord_dec_var_names.end());
	vector<string> prob_pars;
    //disabled the exception for adj pars tied to dec vars
    //string partied;
    //map<string,pair<string,double>> tied_info = pest_scenario.get_base_par_tran_seq().get_tied_ptr()->get_items();

    for (auto& name : pest_scenario.get_ctl_ordered_par_names())
	{
		//if this parameter is not a decision var
		//if (find(start, end, name) == end)
		if (dec_set.find(name) == dec_set.end())
		{
			ParameterRec::TRAN_TYPE tt = pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(name)->tranform_type;
			if ((tt == ParameterRec::TRAN_TYPE::LOG) || (tt == ParameterRec::TRAN_TYPE::NONE))
            {
                adj_par_names.push_back(name);
            }
//            else if (tt == ParameterRec::TRAN_TYPE::TIED)
//            {
//                partied = tied_info[name].first;
//                if (dec_set.find(partied) != dec_set.end())
//                    prob_pars.push_back(name);
//            }
		}

	}
//     if (!prob_pars.empty())
//     {
//         throw_constraints_error("the following adjustable pars are 'tied' to decision variables",prob_pars);
//     }


	//------------------------------------------
	//  ---  chance constraints  ---
	//------------------------------------------
	risk = pest_scenario.get_pestpp_options().get_opt_risk();
	if (risk != 0.5)
	{
		pfm.log_event("initializing chance constraints/objectives");
		use_chance = true;
		

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
        initialize_chance_schedule(f_rec);
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
			string obs_csv = pest_scenario.get_pestpp_options().get_opt_obs_stack();
			//if a par stack wasn't passed, draw one
			if (par_csv.size() == 0)
			{
				if ((stack_size == 0) && (obs_stack_name.size() == 0))
					throw_constraints_error("++opt_stack_size is zero");
				else if (stack_size > 0)
				{
					f_rec << "drawing " << stack_size << " stack realizations" << endl;
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
						stack_pe.from_csv(par_csv, true);
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
					pfm.log_event("loading par stack from binary file " + par_csv);
					try
					{
						stack_pe.from_binary(par_csv, true);
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
				stack_pe.clear_fixed_map();

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
				
				//a restart with "all"
				if ((obs_csv.size() > 0) && (pest_scenario.get_pestpp_options().get_opt_chance_points() == "ALL"))
				{

					stack_pe.transform_ip(ParameterEnsemble::transStatus::CTL);
					vector<string> tokens, keep_org_names, keep_new_names, var_names = stack_pe.get_var_names();
					string member_name, last_member_name = "";
					for (auto n : stack_pe.get_real_names())
					{
						if (n.find_first_of(NESTED_STACK_NAME_SEP) == string::npos)
							throw_constraints_error("no nested stack name sep chars ('" +
								NESTED_STACK_NAME_SEP + "') in par stack name: " + n);
						else
						{
							tokens.clear();
							pest_utils::tokenize(n,tokens, NESTED_STACK_NAME_SEP);
							if (tokens.size() > 2)
							{
								throw_constraints_error("too many nested stack name sep chars ('" +
									NESTED_STACK_NAME_SEP + "') in par stack name: " + n);
							}
							if (tokens.size() < 2)
							{
								throw_constraints_error("too few entries in split nested stack name: "+n+
									", expecting <stack_name>"+NESTED_STACK_NAME_SEP+"<member_name>");
							}
							member_name = tokens[1];
							if (member_name != last_member_name)
							{
								Parameters ctl_pars = pest_scenario.get_ctl_parameters();
								Eigen::VectorXd real = stack_pe.get_real_vector(n);
								ctl_pars.update_without_clear(var_names, real);
								stack_pe_map[member_name] = ctl_pars;
								last_member_name = member_name;
								keep_org_names.clear();
								keep_new_names.clear();
								keep_org_names.push_back(n);
								keep_new_names.push_back(tokens[0]);
							}
							else
							{
								keep_org_names.push_back(n);
								keep_new_names.push_back(tokens[0]);
							}

						}
					}
					//set the stack pe attribute in case the chances need to be recalc'd
					Eigen::MatrixXd real_mat = stack_pe.get_eigen(keep_org_names, vector<string>());
					stack_pe = ParameterEnsemble(&pest_scenario, &rand_gen, real_mat, keep_new_names, var_names);
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
			if (pest_scenario.get_pestpp_options().get_save_binary())
			{

                if (pest_scenario.get_pestpp_options().get_save_dense())
                {
                    filename = filename + ".bin";
                    stack_pe.to_dense(filename);
                }
                else {
                    filename = filename + ".jcb";
                    stack_pe.to_binary(filename);
                }
			}
			else
			{
				filename = filename + ".csv";
				stack_pe.to_csv(filename);
			}
			f_rec << "saved initial parameter stack to " << filename << endl;
			
			//maybe there was an obs stack passed?
			
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
                if (stack_size == 0)
                {
                    stack_size = stack_oe.shape().first;
                }
				//a restart with "all"
				if ((stack_pe_map.size() > 0) && (pest_scenario.get_pestpp_options().get_opt_chance_points() == "ALL"))
				{
					if (par_csv.size() == 0)
					{
						throw_constraints_error("restarting with opt_chance_points=='all' requires both the par and obs nested stack files");

					}
					vector<string> tokens, keep_org_names, keep_new_names, var_names = stack_pe.get_var_names();
					missing.clear();
					string member_name, last_member_name = "";
					for (auto n : stack_oe.get_real_names())
					{
						if (n.find_first_of(NESTED_STACK_NAME_SEP) == string::npos)
							throw_constraints_error("no nested stack name sep chars ('" +
								NESTED_STACK_NAME_SEP + "') in obs stack name: " + n);
						else
						{
							tokens.clear();
							pest_utils::tokenize(n, tokens, NESTED_STACK_NAME_SEP);
							if (tokens.size() > 2)
							{
								throw_constraints_error("too many nested stack name sep chars ('" +
									NESTED_STACK_NAME_SEP + "') in obs stack name: " + n);
							}
							if (tokens.size() < 2)
							{
								throw_constraints_error("too few entries in split nested stack name: " + n +
									", expecting <stack_name>" + NESTED_STACK_NAME_SEP + "<member_name>");
							}
							member_name = tokens[1];
							if (last_member_name == "")
								last_member_name = member_name;
							if (stack_pe_map.find(member_name) == stack_pe_map.end())
								missing.push_back(member_name);
							else if (member_name != last_member_name)
							{
								Eigen::MatrixXd stack_mat = stack_oe.get_eigen(keep_org_names,vector<string>());
								ObservationEnsemble _oe_stack(&pest_scenario, &rand_gen, stack_mat, keep_new_names, stack_oe.get_var_names());
								stack_oe_map[last_member_name] = _oe_stack;
								last_member_name = member_name;
								keep_org_names.clear();
								keep_org_names.push_back(n);
								keep_new_names.clear();
								keep_new_names.push_back(tokens[0]);
							}
							else
							{
								keep_new_names.push_back(tokens[0]);
								keep_org_names.push_back(n);
							}

						}
					}
					if (keep_org_names.size() > 0)
					{
						Eigen::MatrixXd stack_mat = stack_oe.get_eigen(keep_org_names, vector<string>());
						ObservationEnsemble _oe_stack(&pest_scenario, &rand_gen, stack_mat, keep_new_names, stack_oe.get_var_names());
						stack_oe_map[last_member_name] = _oe_stack;
					}
					
					if (missing.size() > 0)
					{
						if (missing.size() == stack_pe_map.size())
							throw_constraints_error("all nested par stack names are missing from nested obs stack");
						for (auto m : missing)
							stack_pe_map.erase(m);
					}
					if (stack_oe_map.size() == 0)
						throw_constraints_error("no nested obs stack names could be matched to nested par stack names...");
					

				}

				else if (((size_passed) || (stack_pe.shape().first > 0)) && (stack_oe.shape().first > stack_size))
				{
					vector<int> drop_rows;
					for (int i = stack_size - 1; i < stack_oe.shape().first; i++)
						drop_rows.push_back(i);
					f_rec << "droppping " << drop_rows.size() << "realizations from obs stack b/c of ++opt_stack_size req" << endl;
					stack_oe.drop_rows(drop_rows);
				}
				stack_runs_processed = true;
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

	int n = 21;
    for (auto& name : ctl_ord_obs_constraint_names)
        n = max(n,(int)name.size()+1);

	f_rec << endl << "  ---  observation constraints and/or objectives ---  " << endl;
	f_rec << left << setw(n) << "name" << setw(15) << "sense" << setw(15) << "value" << endl;
	for (auto& name : ctl_ord_obs_constraint_names)
	{
		f_rec << setw(n) << left << name;
		f_rec << setw(15) << constraint_sense_name[name];
		f_rec << setw(15) << constraints_obs.get_rec(name) << endl;
	}
	cout << "..." << ctl_ord_obs_constraint_names.size() << " obs-based constraints/objectives, see rec file for listing" << endl;

	if (num_pi_constraints() > 0)
	{
        int n = 21;
        for (auto& name : ctl_ord_pi_constraint_names)
            n = max(n,(int)name.size()+1);
		f_rec << endl << "  ---  prior information constraints   ---  " << endl;
		f_rec << setw(n) << "name" << setw(15) << "sense" << setw(15) << "value" << endl;
		for (auto& name : ctl_ord_pi_constraint_names)
		{
			f_rec << setw(n) << left << name;
			f_rec << setw(15) << constraint_sense_name[name];
			f_rec << setw(15) << constraints_pi.get_pi_rec(name).get_obs_value() << endl;
		}
		cout << "..." << ctl_ord_pi_constraint_names.size() << " pi-based constraints, see rec file for listing" << endl;
	}

	if (use_chance)
	{
		f_rec << endl << endl << "  ---  chance constraint/objective FOSM information  ---  " << endl;
		if (use_fosm)
		{
			f_rec << "...using FOSM-based chance constraints/objectives" << endl;
			f_rec << "...++opt_risk and corresponding probit function value: " << setw(10) << risk << setw(20) << probit_val << endl;
			f_rec << "...number of non-zero weight observations for FOSM calcs: " << num_nz_obs() << endl;
		}
		else
		{ 
			f_rec << "...using stack-based chance constraints/objectives" << endl;
		}
		f_rec << "...number of adjustable parameters for chance constraint/objective calcs: " << num_adj_pars() << endl;
		
		f_rec << "...repeat chance constraint/objective calculations every: " << pest_scenario.get_pestpp_options().get_opt_recalc_fosm_every() << " iterations" << endl << endl;
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
				f_rec << endl << endl << "...Note: No nonzero weight observations found." << endl;
				f_rec << "...Prior constraint/objective uncertainty will be used in FOSM-based calculations" << endl;
			}
			else
			{
				f_rec << "  non-zero weight observations used for conditioning in FOSM-based chance constraint/objective calculations: " << endl;
				i = 0;
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
	else
	{
		cout << "...not using chance constraints" << endl;
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

		if (fosm_mat.nonZeros() == 0)
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
	//this isn't actually used in the shifting since we do a non-parametric shift along the PDF, but 
	//we can still report this variance.  So we dont need to worry about accounting for 
	//every stack if there are some failed runs, etc.
	else
	{
		//pair<map<string, double>, map<string, double>> stack_oe_mean_stdev;
		map<string, double> mean_map, std_map;
		if (stack_oe_map.size() > 0)
		{
			
			ObservationEnsemble _mean_stack = get_stack_mean(stack_oe_map);
			_mean_stack.fill_moment_maps(mean_map, std_map);
		}
		else
			stack_oe.fill_moment_maps(mean_map, std_map);
		for (auto cname : ctl_ord_obs_constraint_names)
		{
			prior_const_var[cname] = std_map[cname] * std_map[cname];
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

//Observations Constraints::get_stack_shifted_chance_constraints()
//{
//	/* one version of this method that doesn't take any args just use the pointer to the
//	current sim constraint values*/
//	return get_stack_shifted_chance_constraints(*current_constraints_sim_ptr);
//}

double Constraints::get_sum_of_violations(Parameters& pars, Observations& obs)
{
	double sum = 0.0;
	for (auto& c : get_unsatified_obs_constraints(obs))
		sum += c.second;
	for (auto& c : get_unsatified_pi_constraints(pars))
		sum += c.second;
	if (sum < 1e-10)
	    sum = 0.0;
	return sum;
}


vector<double> Constraints::get_sum_of_violations(ParameterEnsemble& pe, ObservationEnsemble& oe)
{
	if (pe.shape().first != oe.shape().first)
		throw_constraints_error("get_sum_of_violations(): pe rows != oe rows");
	vector<double> sums;
	for (int i = 0; i < pe.shape().first; i++)
	{

	}
	return sums;
}


ObservationEnsemble Constraints::get_chance_shifted_constraints(ParameterEnsemble& pe, ObservationEnsemble& oe, int gen, string risk_obj, string opt_member)
{
    stringstream ss;
    ss.str("");
    ss << file_mgr_ptr->get_base_filename() << "." << gen << ".population_stack_summary.csv";
    ofstream csv;
    if (!use_fosm) {
        csv.open(ss.str());

        if (csv.bad()) {
            throw_constraints_error("error opening '" + ss.str() + "' for writing");
        }
        csv << "member,map_to_member,constraint,risk,stack_mean,stack_stdev,pre_shift,post_shift" << endl;
    }
	map<string, double> risk_map;
	for (auto& rname : pe.get_real_names())
	{
		risk_map[rname] = risk;
	}
	pe.update_var_map();
	map<string, int> vm = pe.get_var_map();
	if (risk_obj.size() > 0)
	{
		if (vm.find(risk_obj) == vm.end())
			throw_constraints_error("couldn't find 'risk_obj' " + risk_obj + " in pe var names");
		Eigen::VectorXd risk_vec = pe.get_var_vector(risk_obj);
		int i = 0;
		for (auto& rname : pe.get_real_names())
		{
			risk_map[rname] = risk_vec[i];
			i++;
		}
	}
	ObservationEnsemble shifted_oe(oe);//copy
	ofstream& frec = file_mgr_ptr->rec_ofstream();
	map<string,double> sim_mean,sim_std,stack_mean,stack_std;
	if (stack_oe_map.size() == 0)
	{
		
		pfm.log_event("risk-shifting observation population using a single, 'optimal' chance point");
		frec << "risk - shifting observation population using a single, 'optimal' chance point" << endl;
	    stack_oe.fill_moment_maps(stack_mean,stack_std);
		//constraints.presolve_chance_report(iter,);
		Observations sim, sim_shifted;
		Eigen::VectorXd real_vec;
		vector<string> onames = shifted_oe.get_var_names();
		vector<string> rnames = shifted_oe.get_real_names();
		for (int i = 0; i < shifted_oe.shape().first; i++)
		{
			real_vec = shifted_oe.get_real_vector(i);
			sim.update_without_clear(onames, real_vec);
			sim_shifted = get_chance_shifted_constraints(sim, risk_map[rnames[i]]);
			shifted_oe.replace(i, sim_shifted);
			//cout << *shifted_oe.get_eigen_ptr() << endl;
			//cout << endl;
			if (!use_fosm) {
                for (auto &constraint : ctl_ord_obs_constraint_names) {
                    csv << rnames[i] << "," << opt_member << "," << constraint << ",";
                    csv << risk_map[rnames[i]] << "," << stack_mean.at(constraint) << "," << stack_std.at(constraint)
                        << ",";
                    csv << sim.get_rec(constraint) << "," << sim_shifted.get_rec(constraint) << endl;
                }
            }
		}

	}
	else
	{
		pfm.log_event("risk-shifting observation population using 'ALL' decision variable solutions as chance points");
		frec << "risk-shifting observation population using 'ALL' decision variable solutions as chance points" << endl;
		//vector<string> real_names = shifted_oe.get_real_names();
		vector<string> real_names = oe.get_real_names();
		//set<string> snames(real_names.begin(), real_names.end());
		map<string, int> real_map = oe.get_real_map();
		Observations sim, sim_shifted;
		Eigen::VectorXd real_vec;
		vector<string> onames = shifted_oe.get_var_names();
		vector<string> missing;

		for (auto& real_name : real_names)
		{
			//if (snames.find(real_info.first) == snames.end())
			if (stack_oe_map.find(real_name) == stack_oe_map.end())
			{
				missing.push_back(real_name);
			}
			else
			{
				real_vec = shifted_oe.get_real_vector(real_name);
				sim.update_without_clear(onames, real_vec);
				sim_shifted = get_stack_shifted_chance_constraints(sim, stack_oe_map[real_name], risk_map[real_name],true, false);
                stack_oe_map[real_name].fill_moment_maps(stack_mean,stack_std);
				shifted_oe.replace(real_map[real_name], sim_shifted);
                for (auto& constraint : ctl_ord_obs_constraint_names) {
                    csv << real_name << "," << real_name << "," << constraint << ",";
                    csv << risk_map[real_name] << "," << stack_mean.at(constraint) << "," << stack_std.at(constraint)
                        << ",";
                    csv << sim.get_rec(constraint) << "," << sim_shifted.get_rec(constraint) << endl;
                }
			}
		}

		if (missing.size() > 0)
		{
			frec << missing.size() << " members not in current oe stack runs, mapping to nearest point in decision variable space" << endl;
			cout << missing.size() << " members not in current oe stack runs, mapping to nearest point in decision variable space" << endl;




			//so the idea here is to figure out which pe in the pe stack map is closest to the missing solutions in dec var space. 
			if (stack_pe_map.size() != stack_oe_map.size())
				throw_constraints_error("stack_pe_map size != stack_oe_map size");
			//remove the risk dv name if found
			vector<string> dvnames = dec_var_names;
			vector<string>::iterator rit = find(dvnames.begin(),dvnames.end(),RISK_NAME);
			if (rit != dvnames.end())
            {
                dvnames.erase(rit);
            }
			Eigen::MatrixXd missing_dv_mat = pe.get_eigen(missing, dvnames);
			Eigen::VectorXd missing_dv_vec, stack_dv_vec;
			double dist, min_dist, max_dist,cutoff;
			string min_real_name,missing_real_name;
			vector<double> distances;
			vector<double> factors,temp, temp2;
			vector<string> dreal_names;

            Eigen::MatrixXd shifts;
            double factor_sum, factor, factor_sum2;
			for (int i=0;i<missing.size();i++)
			{
                missing_real_name = missing[i];
				missing_dv_vec = missing_dv_mat.row(i);
				min_dist = numeric_limits<double>::max();
				min_real_name = "";
				distances.clear();
                double _risk = risk;

                if (risk_obj.size() > 0)
                {
                    //_risk = stack_pe_map[min_real_name][risk_obj];
                    _risk = risk_map[missing_real_name];
                }

				for (auto& p : stack_pe_map)
				{
					stack_dv_vec = p.second.get_data_eigen_vec(dvnames);
					dist = (stack_dv_vec - missing_dv_vec).squaredNorm();
					distances.push_back(dist);
					dreal_names.push_back(p.first);
					if (dist < min_dist)
					{
						min_dist = dist;
						min_real_name = p.first;
					}
				}
                temp.clear();
				max_dist = *max_element(distances.begin(),distances.end());
                for (auto& d : distances)
                    temp.push_back(d / max_dist);
                max_dist = *max_element(temp.begin(),temp.end());
                temp2 = temp;
                sort(temp2.begin(),temp2.end());
                cutoff = temp2[temp2.size()-1];
                if (temp2.size() > 5)
                {
                    cutoff = temp2[4];
                }
                factors.clear();
                factor_sum = 0.0;
                for (auto& t : temp)
                {

                    if (t == 0)
                        factor = 10000.0;
                    else if (t > cutoff)
                    {
                        factor = 0.0;
                    }
                    else
                        factor = 1.0/t;
                    factor_sum += factor;
                    factors.push_back(factor);

                }
                temp = factors;
                factors.clear();
                factor_sum2 = 0.0;
                for (auto& t : temp)
                {
                    factor = t / factor_sum;
                    factors.push_back(factor);
                    factor_sum2 += factor;

                }

                shifts.resize(factors.size(),shifted_oe.shape().second);
                shifts.setZero();
                real_vec = shifted_oe.get_real_vector(missing[i]);
                sim.update_without_clear(onames, real_vec);
                for (int ii=0;ii<factors.size();ii++)
                {
                    if (factors[ii] == 0.0)
                    {
                        continue;
                    }
                    //this call uses the class stack_oe attribute;
                    sim_shifted = get_stack_shifted_chance_constraints(sim, stack_oe_map[dreal_names[ii]], _risk, true, false);
                    shifts.row(ii) = sim_shifted.get_data_eigen_vec(shifted_oe.get_var_names()) * factors[ii];
                }
                real_vec = shifts.colwise().sum();
                shifted_oe.get_eigen_ptr_4_mod()->row(real_map.at(missing[i])) = real_vec;

				if (min_real_name.size() == 0)
					//do something here
					throw_constraints_error("couldn't find a nearest dv real for stack mapping");
				if (stack_oe_map.find(min_real_name) == stack_oe_map.end())
					throw_constraints_error("nearest dv real '" + min_real_name +"' not in stack oe map");
				if (min_dist > 0.0) min_dist = sqrt(min_dist);


                if (pest_scenario.get_pestpp_options().get_mou_verbose_level() > 3) {
                    for (int ii = 0; ii < factors.size(); ii++) {
                        frec << "member '" << missing[i] << "' mapped to '" << dreal_names[ii] << "' at a distance of "
                             << distances[ii] << " with a IDW factor of " << factors[ii] << " for stack-based chances"
                             << endl;
                    }

                }

				real_vec = shifted_oe.get_real_vector(missing[i]);
				sim.update_without_clear(onames, real_vec);
				//this call uses the class stack_oe attribute;
				sim_shifted = get_stack_shifted_chance_constraints(sim, stack_oe_map[min_real_name], _risk, true, false);
				//shifted_oe.replace(real_map.at(missing[i]), sim_shifted);
                stack_oe_map[min_real_name].fill_moment_maps(stack_mean,stack_std);
                for (auto& constraint : ctl_ord_obs_constraint_names) {
                    csv << missing[i] << "," << min_real_name << "," << constraint << ",";
                    csv << _risk << "," << stack_mean.at(constraint) << "," << stack_std.at(constraint)
                        << ",";
                    csv << sim.get_rec(constraint) << "," << sim_shifted.get_rec(constraint) << endl;
                }
			}
		}
	}
    if (!use_fosm)
    {
        csv.close();
    }
	return shifted_oe;
}

Observations Constraints::get_chance_shifted_constraints(Observations& current_obs)
{
	return get_chance_shifted_constraints(current_obs, risk);

}



Observations Constraints::get_chance_shifted_constraints(Observations& current_obs, double _risk, bool use_stack_anomalies)
{
	/* get the simulated constraint values with the chance shift applied*/
	ofstream& f_rec = file_mgr_ptr->rec_ofstream();
	prior_constraint_offset.clear();
	prior_constraint_stdev.clear();
	post_constraint_offset.clear();
	post_constraint_stdev.clear();
	double new_constraint_val, old_constraint_val, required_val;
	double pr_offset, pt_offset;
	Observations shifted_obs;
	double _probit_val = get_probit(_risk);
	if (use_fosm)
	{
		for (auto& name : ctl_ord_obs_constraint_names)
		{
			prior_constraint_stdev[name] = sqrt(prior_const_var.at(name));
			post_constraint_stdev[name] = sqrt(post_const_var.at(name));
			//the offset (shift) is just the stdev * the risk-based probit value
			pr_offset = _probit_val * prior_constraint_stdev.at(name);
			pt_offset = _probit_val * post_constraint_stdev.at(name);
			old_constraint_val = current_obs.get_rec(name);
			required_val = constraints_obs.get_rec(name);

			//if less_than constraint, then add to the sim value, to move positive
			// WRT the required constraint value
			if (constraint_sense_map[name] == ConstraintSense::less_than)
			{

				new_constraint_val = old_constraint_val + pt_offset;
				post_constraint_offset[name] = pt_offset;
				prior_constraint_offset[name] = pr_offset;
			}
			//if greater_than constraint, the subtract from the sim value to move
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
			shifted_obs.insert(name, new_constraint_val);
		}
	}
	else
	{
		if (!stack_runs_processed)
			throw_constraints_error("no stack runs have been processed, something is wrong");
		//if we are using nested stacks, then calculate the average stack
		if (stack_oe_map.size() > 0)
		{
			ObservationEnsemble _mean_stack = get_stack_mean(stack_oe_map);
			shifted_obs = get_stack_shifted_chance_constraints(current_obs, _mean_stack, _risk, false, false);
		}

		else
			shifted_obs = get_stack_shifted_chance_constraints(current_obs, stack_oe, _risk,false, use_stack_anomalies);
		
	}
	vector<string> names = shifted_obs.get_keys();
	Observations constraints_chance(current_obs);
	constraints_chance.update_without_clear(names, shifted_obs.get_data_vec(names));
	
	return constraints_chance;
}

pair<Eigen::VectorXd, Eigen::VectorXd> Constraints::get_obs_resid_constraint_vectors(Parameters& _current_ctl_dv_vals, Observations& _constraints_sim, vector<string>& cnames)
{
	Eigen::VectorXd rhs(cnames.size()), resid(cnames.size());
	Observations sim = _constraints_sim;
	if (use_chance)
		sim = get_chance_shifted_constraints(_constraints_sim);
	string cname;
	Observations obs = pest_scenario.get_ctl_observations();
	for (int i=0;i<cnames.size();i++)
	{
		cname = cnames[i];
		if (sim.find(cname) != sim.end())
		{
			rhs[i] = obs[cname];
			resid[i] = obs[cname] - sim[cname];
		}
		else if (constraints_pi.find(cname) != constraints_pi.end())
		{
			resid[i] = constraints_pi.get_pi_rec(cname).calc_sim_and_resid(_current_ctl_dv_vals).second;
			rhs[i] = constraints_pi.get_pi_rec(cname).get_obs_value();
		}
		else // todo Jdub else needed here right?
		{
			throw_constraints_error("get_obs_resid_cosntraint_vectors(): constraint '" + cname + "' not found");
		}
	}
	
	return pair<Eigen::VectorXd, Eigen::VectorXd>(rhs,resid);
}


vector<string> Constraints::get_working_set_ineq_names(vector<string>& cnames)
{
	vector<string> names;
	for (auto& name : cnames)
		if (constraint_sense_map[name] != ConstraintSense::equal_to)
			names.push_back(name);
	return names;
}

ObservationEnsemble Constraints::get_stack_mean(map<string, ObservationEnsemble>& _stack_oe_map)
{
	int max_nrow = 0;
	for (auto& o : _stack_oe_map)
	{
		max_nrow = max(max_nrow, o.second.shape().first);
	}
	int count = 0;
	Eigen::MatrixXd oe_stack_mean;
	bool first = true;
	vector<string> real_names;
	for (auto& o : _stack_oe_map)
	{
		//only use stacks that are full size (no failed runs)
		if (o.second.shape().first != max_nrow)
			continue;

		if (first)
		{
			oe_stack_mean = o.second.get_eigen(vector<string>(), ctl_ord_obs_constraint_names);
			real_names = o.second.get_real_names();
			first = false;
		}

		else
		{
			oe_stack_mean = oe_stack_mean + o.second.get_eigen(vector<string>(), ctl_ord_obs_constraint_names);
		}
		count += 1;
	}
	if (first)
		throw_constraints_error("unable to compute oe_stack_mean - all stacks are empty");
	oe_stack_mean = oe_stack_mean / double(count);
	ObservationEnsemble _mean_stack(&pest_scenario, &rand_gen, oe_stack_mean, real_names, ctl_ord_obs_constraint_names);
	return _mean_stack;
}

Observations Constraints::get_stack_shifted_chance_constraints(Observations& current_obs, ObservationEnsemble& _stack_oe,
                                                               double _risk, bool full_obs, bool use_stack_anomalies)
{
	double old_constraint_val, required_val, pt_offset, new_constraint_val;

	Observations shifted_obs;
	if (full_obs)
		shifted_obs = Observations(current_obs);
	//work out which realization index corresponds to the risk value
	if (_stack_oe.shape().first < 3)
		throw_constraints_error("too few (<3) stack members, cannot continue with stack-based chance constraints/objectives");
	int cur_num_reals = _stack_oe.shape().first;
	//get the index (which realization number) represents the risk value according to constraint sense
	//
	//the "less than" realization index is just the risk value times the number of reals
	int lt_idx = int(_risk * cur_num_reals);
	//the greater than index is the opposite direction
	int gt_idx = int((1.0 - _risk) * cur_num_reals);
	//the equality constraint risk index
	int eq_idx = int(0.5 * cur_num_reals);
	//get the mean-centered anomalies - we want to subtract off the mean 
	//in case these stack values are being reused from a previous iteration
	//Eigen::MatrixXd anom = _stack_oe.get_eigen_anomalies();
    Eigen::MatrixXd anom;
    if (use_stack_anomalies)
        anom = _stack_oe.get_eigen_anomalies();
    else
        anom = _stack_oe.get_eigen();

    //get the map of realization name to index location in the stack
	_stack_oe.update_var_map();
	map<string, int> var_map = _stack_oe.get_var_map();
	//get the mean and stdev summary containers, summarized by observation (e.g. constraint) name
	//pair<map<string, double>, map<string, double>> mm = _stack_oe.get_moment_maps();
	map<string, double> mean_map, std_map;
	_stack_oe.fill_moment_maps(mean_map, std_map);
    Eigen::VectorXd cvec;
	for (auto& name : ctl_ord_obs_constraint_names)
	{
		old_constraint_val = current_obs.get_rec(name);
		//the value that must be statified (from the control file)
		required_val = constraints_obs[name];
		// the realized values of this stack are the anomalies added to the
		//current constraint value - this assumes the current value 
		//is the mean of the stack distribution
        if (use_stack_anomalies)
            cvec = anom.col(var_map[name]).array() + old_constraint_val;
        else
		    cvec = anom.col(var_map[name]).array();// + old_constraint_val;
		//now sort the anomalies + current (mean) value vector
		sort(cvec.data(), cvec.data() + cvec.size());
		/*cout << cvec << endl << endl;
        Eigen::VectorXd temp = _stack_oe.get_var_vector(name);
        sort(temp.data(),temp.data() + temp.size());
        cout << temp << endl << endl;*/

        //set the stdev container info - this isn't used in
		//calculations but gets reported
		prior_constraint_stdev[name] = std_map[name];
		post_constraint_stdev[name] = std_map[name];

		if ((prior_constraint_stdev[name] == 0.0) && (pest_scenario.get_pestpp_options().get_mou_verbose_level() > 3))
		{
			throw_constraints_error("model-based constraint '" + name + "' has empirical (stack) standard deviation of 0.0", false);
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
		shifted_obs[name] = new_constraint_val;
	}
	return shifted_obs;
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

pair<vector<double>,vector<double>> Constraints::get_constraint_bound_vectors(Parameters& current_pars, Observations& current_obs, bool use_stack_anomalies)
{
	/* get the upper and lower bound constraint vectors. For less than constraints, the lower bound is 
	set to double max, for greater than constraints, the upper bound is set to double max.
	These are needed for the simplex solve*/
	vector<double> residuals;
	if (use_chance)
	{
		Observations current_constraints_chance = get_chance_shifted_constraints(current_obs, risk, use_stack_anomalies);
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
		double residual = -constraints_pi.get_pi_rec(name).calc_residual(current_pars);
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


void Constraints::throw_constraints_error(string message, const vector<string>& messages, bool should_throw)
{
	stringstream ss;
	for (auto& mess : messages)
		ss << mess + ',';
	throw_constraints_error(message + ss.str(), should_throw);
}

void Constraints::throw_constraints_error(string message, const set<string>& messages, bool should_throw)
{
	stringstream ss;
	for (auto& mess : messages)
		ss << mess + ',';
	throw_constraints_error(message + ss.str(), should_throw);
}

void Constraints::throw_constraints_error(string message, bool should_throw)
{
	string error_message;
	if (should_throw)
		error_message = "error in Constraints: " + message;
	else
		error_message = "Constraints Warning: " + message;
	file_mgr_ptr->rec_ofstream() << error_message << endl;
	
	
	if (should_throw)
	{
		cout << endl << endl << error_message << endl << endl;
		throw runtime_error(error_message);
		//file_mgr_ptr->close_file("rec");
	}
	else
	{
		cout << error_message << endl;
	}
}

double  Constraints::ErfInv2(double x)
{
	/* the inverse error function, needed for the FOSM-based chance constraints*/
	float tt1, tt2, lnx, sgn;
	sgn = (x < 0) ? -1.0f : 1.0f;

	x = (1 - x) * (1 + x);        // x = 1 - x*x;
	lnx = logf(x);
	double PI = 3.14159265358979323846;
	tt1 = 2 / (PI * 0.147) + 0.5f * lnx;
	tt2 = 1 / (0.147) * lnx;

	return(sgn * sqrtf(-tt1 + sqrtf(tt1 * tt1 - tt2)));
}

double Constraints::get_probit()
{
	return get_probit(risk);
}

double Constraints::get_probit(double _risk)
{
	/* the probit function estimate
	 * - needed for the fosm-basd chance constraints*/
	double output = sqrt(2.0) * ErfInv2((2.0 * _risk) - 1.0);
	return output;
}

void Constraints::sqp_report(int iter, Parameters& current_pars, Observations& current_obs, bool echo, string tag)
{
	
	ofstream& f_rec = file_mgr_ptr->rec_ofstream();
	vector<double> residuals;
	map<string, double> infeas_dist;
	stringstream ss;
	ss.str("");
	
	int nsize = 20;
	for (auto name : ctl_ord_obs_constraint_names)
	{
		nsize = max(nsize, int(name.size()));
	}
	ss << endl << "  observation constraint information at iteration " << iter;
	if (tag.size() > 0)
		ss << " for " << tag;
	ss << endl;
	ss << setw(nsize) << left << "name" << right << setw(14) << "sense" << setw(12) << "required" << setw(15) << "sim value";
	ss << setw(11) << "satisfied" << setw(15) << "distance" << endl;

	infeas_dist = get_unsatified_obs_constraints(current_obs, 0.0, false);
	for (int i = 0; i < num_obs_constraints(); ++i)
	{
		string name = ctl_ord_obs_constraint_names[i];
		ss << setw(nsize) << left << name;
		ss << setw(14) << right << constraint_sense_name[name];
		ss << setw(12) << constraints_obs.get_rec(name);
		ss << setw(15) << current_obs.get_rec(name);
		if (infeas_dist.find(name) != infeas_dist.end())
		{
			ss << setw(11) << "false" << setw(15) << infeas_dist[name];
		}
		else
		{
			ss << setw(11) << "true" << setw(15) << 0.0;
		}
		ss << endl;

	}
	
	if (ctl_ord_pi_constraint_names.size() > 0)
	{


		nsize = 20;
		for (auto name : ctl_ord_pi_constraint_names)
		{
			nsize = max(nsize, int(name.size()));
		}

		//report prior information constraints
		infeas_dist = get_unsatified_pi_constraints(current_pars);
		ss << endl << " --- prior information constraint summary at iteration " << iter;
		if (tag.size() > 0)
			ss << " for " << tag;
		ss << " --- " <<  endl;
		ss << setw(nsize) << left << "name" << right << setw(14) << "sense" << setw(12) << "required" << setw(15) << "sim value";
		ss << setw(15) << "satisfied" << setw(15) << "distance" << endl;
		for (int i = 0; i < num_pi_constraints(); ++i)
		{
			string name = ctl_ord_pi_constraint_names[i];
			PriorInformationRec pi_rec = constraints_pi.get_pi_rec(name);
			ss << setw(nsize) << left << name;
			ss << setw(14) << right << constraint_sense_name[name];
			ss << setw(12) << pi_rec.get_obs_value();
			ss << setw(15) << pi_rec.calc_sim_and_resid(current_pars).first;
			if (infeas_dist.find(name) != infeas_dist.end())
			{
				ss << setw(11) << "false" << setw(15) << infeas_dist[name];
			}
			else
			{
				ss << setw(11) << "true" << setw(15) << 0.0;
			}
			ss << endl;
		}
	}
	ss << endl;
	f_rec << ss.str();
	if (echo)
		cout << ss.str();

	return;
}

void Constraints::mou_report(int iter, Parameters& current_pars, Observations& current_obs, const vector<string>& obs_obj_names,
	const vector<string>& pi_obj_names, bool echo)
{

	set<string> skip_names(obs_obj_names.begin(), obs_obj_names.end());

	ofstream& f_rec = file_mgr_ptr->rec_ofstream();
	vector<double> residuals;
	map<string, double> infeas_dist;
	stringstream ss;
	ss.str("");
	int obs_infeas = 0, pi_infeas = 0;
	if (skip_names.size() < ctl_ord_obs_constraint_names.size())
	{
		int nsize = 20;
		for (auto name : ctl_ord_obs_constraint_names)
		{
			if (skip_names.find(name) != skip_names.end())
				continue;
			nsize = max(nsize, int(name.size()));
		}
		ss << endl << "  observation constraint information at iteration " << iter << endl;
		ss << setw(nsize) << left << "name" << right << setw(42) << "sense" << setw(12) << "required" << setw(15) << "sim value";
		ss << setw(11) << "satisfied" << setw(15) << "distance" << endl;

		

		infeas_dist = get_unsatified_obs_constraints(current_obs, 0.0, false);
		obs_infeas = infeas_dist.size();
		for (int i = 0; i < num_obs_constraints(); ++i)
		{
			string name = ctl_ord_obs_constraint_names[i];
			if (skip_names.find(name) != skip_names.end())
				continue;
			ss << setw(nsize) << left << name;
			ss << setw(14) << right << constraint_sense_name[name];
			ss << setw(12) << constraints_obs.get_rec(name);
			ss << setw(15) << current_obs.get_rec(name);
			if (infeas_dist.find(name) != infeas_dist.end())
			{
				ss << setw(11) << "false" << setw(15) << infeas_dist[name];
			}
			else
			{
				ss << setw(11) << "true" << setw(15) << 0.0;
			}
			ss << endl;

		}
	}
	skip_names.clear();
	skip_names.insert(pi_obj_names.begin(), pi_obj_names.end());
	if (num_pi_constraints() > skip_names.size())
	{

		int nsize = 20;
		for (auto name : ctl_ord_pi_constraint_names)
		{
			if (skip_names.find(name) != skip_names.end())
				continue;
			nsize = max(nsize, int(name.size()));
		}

		//report prior information constraints
		infeas_dist = get_unsatified_pi_constraints(current_pars);
		pi_infeas = infeas_dist.size();
		ss << endl << " --- prior information constraint information at iteration " << iter << " --- " << endl;
		ss << setw(nsize) << left << "name" << right << setw(14) << "sense" << setw(12) << "required" << setw(15) << "sim value";
		ss << setw(15) << "satisfied" << setw(15) << "distance" << endl;
		for (int i = 0; i < num_pi_constraints(); ++i)
		{
			string name = ctl_ord_pi_constraint_names[i];
			if (skip_names.find(name) != skip_names.end())
				continue;
			PriorInformationRec pi_rec = constraints_pi.get_pi_rec(name);
			ss << setw(nsize) << left << name;
			ss << setw(14) << right << constraint_sense_name[name];
			ss << setw(12) << pi_rec.get_obs_value();
			ss << setw(15) << pi_rec.calc_sim_and_resid(current_pars).first;
			if (infeas_dist.find(name) != infeas_dist.end())
			{
				ss << setw(11) << "false" << setw(15) << infeas_dist[name];
			}
			else
			{
				ss << setw(11) << "true" << setw(15) << 0.0;
			}
			ss << endl;
		}
	}
	f_rec << ss.str();
	if (echo)
		cout << ss.str();
	else
		cout << "iteration " << iter << ": " << obs_infeas << " obs constraints and " << pi_infeas << " prior info constraints not satisfied, see rec file for listing" << endl;
	return;
}

string Constraints::mou_population_observation_constraint_summary(int iter, ObservationEnsemble& oe, string tag, const vector<string>& obs_obj_names)
{
    set<string> skip_names(obs_obj_names.begin(), obs_obj_names.end());
    vector<double> residuals;
    map<string, double> infeas_dist;
    stringstream ss;
    ss.str("");
    vector<string> names = oe.get_var_names();
    Observations current_obs;
    map<string, int> infeas_count;
    if (skip_names.size() < ctl_ord_obs_constraint_names.size())
    {
        int nsize = 20;
        for (auto o : ctl_ord_obs_constraint_names)
        {
            infeas_count[o] = 0;
            nsize = max(nsize, int(o.size()));
        }
        for (auto real : oe.get_real_names())
        {
            current_obs.update(names, eigenvec_2_stlvec(oe.get_real_vector(real)));
            infeas_dist = get_unsatified_obs_constraints(current_obs, 0.0, false);
            for (auto in : infeas_dist)
                infeas_count[in.first]++;
        }

        //pair<map<string,double>,map<string,double>> mm = oe.get_moment_maps();
        map<string, double> mean_map, std_map;
        oe.fill_moment_maps(mean_map, std_map);

        ss << endl << "  ---  " << tag << " population observation constraint summary at iteration " << iter << "  ---  " << endl;
        ss << setw(nsize) << left << "name" << right << setw(14) << "sense" << setw(12) << "required" << setw(15) << "sim min";
        ss << setw(15) << "sim mean" << setw(15) << "sim max" << setw(15) << "% unsatisfied" << setw(15) << endl;

        infeas_dist = get_unsatified_obs_constraints(current_obs, 0.0, false);
        int num_reals = oe.shape().first;
        for (int i = 0; i < num_obs_constraints(); ++i)
        {
            string name = ctl_ord_obs_constraint_names[i];
            if (skip_names.find(name) != skip_names.end())
                continue;
            ss << setw(nsize) << left << name;
            ss << setw(14) << right << constraint_sense_name[name];
            ss << setw(12) << constraints_obs.get_rec(name);

            ss << setw(15) << oe.get_var_vector(name).minCoeff();
            ss << setw(15) << mean_map[name];
            ss << setw(15) << oe.get_var_vector(name).maxCoeff();
            ss << setw(15) << 100.0 * double(infeas_count[name]) / double(num_reals);
            ss << endl;

        }
    }
    return ss.str();

}


void Constraints::mou_report(int iter, ParameterEnsemble& pe, ObservationEnsemble& oe, const vector<string>& obs_obj_names,
	const vector<string>& pi_obj_names, bool echo)
{

	set<string> skip_names;

	ofstream& f_rec = file_mgr_ptr->rec_ofstream();
	vector<double> residuals;
	map<string, double> infeas_dist;
	map<string,int> infeas_count;
	stringstream ss;
	ss.str("");
	string tag = "";
	if (use_chance)
	    tag = "chance-shifted";
	string op_sum = mou_population_observation_constraint_summary(iter,oe,tag,obs_obj_names);
	ss << op_sum;
	skip_names.clear();
	skip_names.insert(pi_obj_names.begin(), pi_obj_names.end());
	if (num_pi_constraints() > skip_names.size())
	{
		
		Parameters current_pars;
		vector<string> names = pe.get_var_names();
		map<string, double> tots;
		map<string, double> mins,maxs;
		PriorInformationRec pi_rec;
		int nsize = 20;
		double v;
		for (auto p : ctl_ord_pi_constraint_names)
		{
			tots[p] = 0.0;
			mins[p] = 1.0e+300;
			maxs[p] = -1.0e+300;
			infeas_count[p] = 0;
			nsize = max(nsize, int(p.size()));
		}
		for (auto real : pe.get_real_names())
		{
			current_pars.update(names, eigenvec_2_stlvec(pe.get_real_vector(real)));
			infeas_dist = get_unsatified_pi_constraints(current_pars);
			for (int i = 0; i < num_pi_constraints(); ++i)
			{
				string name = ctl_ord_pi_constraint_names[i];
				pi_rec = constraints_pi.get_pi_rec(name);
				v = pi_rec.calc_sim_and_resid(current_pars).first;
				tots[name] = tots[name] + v;
				mins[name] = min(mins[name], v);
				maxs[name] = max(maxs[name], v);



			}
            for (auto in : infeas_dist)
                infeas_count[in.first]++;
		}
		//report prior information constraints
		ss << endl << " --- prior information constraint summary at iteration " << iter << " --- " << endl;
		ss << setw(nsize) << left << "name" << right << setw(14) << "sense" << setw(12) << "required" << setw(15) << "sim min";
		ss << setw(15) << "sim mean" << setw(15) << "sim max" << setw(15) << "% unsatisfied" << endl;
		int num_reals = pe.shape().first;
		for (int i = 0; i < num_pi_constraints(); ++i)
		{
			string name = ctl_ord_pi_constraint_names[i];
			if (skip_names.find(name) != skip_names.end())
				continue;
			pi_rec = constraints_pi.get_pi_rec(name);
			ss << setw(nsize) << left << name;
			ss << setw(14) << right << constraint_sense_name[name];
			ss << setw(12) << pi_rec.get_obs_value();
			ss << setw(15) << mins[name];
			ss << setw(15) << double(tots[name] / num_reals);
			ss << setw(15) << maxs[name];
			ss << setw(15) << 100.0 * double(infeas_count[name]) / double(num_reals);
			ss << endl;
		}
	}
	f_rec << ss.str();
	if (echo)
		cout << ss.str();

	return;
}


void Constraints::presolve_report(int iter, Parameters& current_pars, Observations& current_obs)
{
	/* this is a report to the rec file for the status of the constraints before solving the current iteration*/
	int nsize = 20;
	for (auto o : ctl_ord_obs_constraint_names)
		nsize = max(nsize, int(o.size()));

	ofstream& f_rec = file_mgr_ptr->rec_ofstream();
	vector<double> residuals;
	f_rec << endl << "  observation constraint/objective information at start of iteration " << iter << endl;
	f_rec << setw(nsize) << left << "name" << right << setw(14) << "sense" << setw(12) << "required" << setw(15) << "sim value";
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
		f_rec << setw(nsize) << left << name;
		f_rec << setw(14) << right << constraint_sense_name[name];
		f_rec << setw(12) << constraints_obs.get_rec(name);	
		f_rec << setw(15) << current.get_rec(name);
		f_rec << setw(15) << residuals[i];
		f_rec << setw(15) << current_bounds.first[i];
		f_rec << setw(15) << current_bounds.second[i] << endl;
	}
	if (num_pi_constraints() > 0)
	{
		nsize = 20;
		for (auto o : ctl_ord_pi_constraint_names)
			nsize = max(nsize, int(o.size()));
		//report prior information constraints
		f_rec << endl << " --- prior information constraint summary at start of iteration " << iter << " --- " << endl;
		f_rec << setw(nsize) << left << "name" << right << setw(14) << "sense" << setw(12) << "required" << setw(15) << "sim value";
		f_rec << setw(15) << "residual" << setw(15) << "lower bound" << setw(15) << "upper bound" << endl;
		for (int i = 0; i < num_pi_constraints(); ++i)
		{
			string name = ctl_ord_pi_constraint_names[i];
			PriorInformationRec pi_rec = constraints_pi.get_pi_rec(name);
			f_rec << setw(nsize) << left << name;
			f_rec << setw(14) << right << constraint_sense_name[name];
			f_rec << setw(12) << pi_rec.get_obs_value();
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
	/* one form of this method that doesn't require advanced knowledge of whether or not chances are being used.
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

void Constraints::stack_summary(int iter, Observations& shifted_obs, bool echo, string header)
{
	/* write chance info to the rec file before undertaking the current iteration process*/
	if (!use_chance)
		return;

	//Observations current_constraints_chance = get_stack_shifted_chance_constraints(current_obs);

	int nsize = 20;
	for (auto o : ctl_ord_obs_constraint_names)
		nsize = max(nsize, int(o.size()));

	stringstream ss;
	ss << endl << "   ";
	if (header.size() == 0)
		//dont change this text - its used in one of the opt tests
		ss << "Chance constraint/objective information at start of iteration " << iter;
	else
		ss << header;
	ss << endl;

	//vector<double> residuals = get_constraint_residual_vec();

	ss << setw(nsize) << left << "name" << right << setw(14) << "sense" << setw(12) << "required";
	ss << setw(12) << "stdev" << setw(12) << "offset";
	ss << setw(12) << "shifted val" << endl;
	vector<string> out_of_bounds;
	for (int i = 0; i < num_obs_constraints(); ++i)
	{
		string name = ctl_ord_obs_constraint_names[i];
		ss << setw(nsize) << left << name;
		ss << setw(14) << right << constraint_sense_name[name];
		ss << setw(12) << constraints_obs[name];
		ss << setw(12) << post_constraint_stdev[name];
		ss << setw(12) << post_constraint_offset[name];
		ss << setw(12) << shifted_obs[name] << endl;
	}

	ss << "  note: The above standard deviations and offset are empirical estimates from the stack." << endl;
	ss <<"         The 'required' value only applies to constraint quantities." << endl;
	
	file_mgr_ptr->rec_ofstream() << ss.str();
	if (echo)
		cout << ss.str();
	return;

}


void Constraints::presolve_chance_report(int iter, Observations& current_obs, bool echo, string header)
{
	/* write chance info to the rec file before undertaking the current iteration process*/
	if (!use_chance)
		return;
	
	if (!use_fosm)
	{
		stack_summary(iter, current_obs, echo, header);
		return;
	}

	int nsize = 20;
	for (auto o : ctl_ord_obs_constraint_names)
		nsize = max(nsize, int(o.size()));

	stringstream ss;
	ss << endl << "   ";
	if (header.size() == 0)
		//dont change this text - its used in one of the opt tests
		ss << "Chance constraint/objective information at start of iteration " << iter;
	else
		ss << header;
	ss << endl;

	//vector<double> residuals = get_constraint_residual_vec();
	
	ss << setw(nsize) << left << "name" << right << setw(14) << "sense" << setw(12) << "required" << setw(12) << "sim value";
	ss << setw(12) << "prior stdev" << setw(12) << "post stdev" << setw(12) << "offset";
	ss << setw(14) << "new sim value" << endl;
	vector<string> out_of_bounds;
	Observations current_constraints_chance = get_chance_shifted_constraints(current_obs);
	for (int i = 0; i < num_obs_constraints(); ++i)
	{
		string name = ctl_ord_obs_constraint_names[i];
		ss << setw(nsize) << left << name;
		ss << setw(14) << right << constraint_sense_name[name];
		ss << setw(12) << constraints_obs[name];
		ss << setw(12) << current_obs.get_rec(name);
		ss << setw(12) << prior_constraint_stdev[name];
		ss << setw(12) << post_constraint_stdev[name];
		ss << setw(12) << post_constraint_offset[name];
		ss << setw(14) << current_constraints_chance[name] << endl;
	}
	ss << "  note: 'offset' is the value added to the simulated constraint/objective value to account" << endl;
	ss << "        for the uncertainty in the constraint/objective value arsing from uncertainty in the " << endl;
	ss << "        adjustable parameters identified in the control file." << endl;
    ss <<"         The 'required' value only applies to constraint quantities." << endl << endl;
	if (!use_fosm)
	{
		ss << "  note: the above standard deviations are empirical estimates from the stack" << endl;
	}
	file_mgr_ptr->rec_ofstream() << ss.str();
	if (echo)
		cout << ss.str();
	return;
	
}


bool Constraints::should_update_chance(int iter)
{
	/* this is a total hack - it tries to determine if it is time to update the chance
	information, like should we queue up some JCO perturbation runs for FOSM chances or 
	should we queue up some stack runs for stack-based chances*/
	if (!use_chance)
		return false;
	if (std_weights)
		return false;
	if (iter == 0)
	{
		if (use_fosm)
			return true;
		if (pest_scenario.get_pestpp_options().get_opt_obs_stack().size() > 0)
			return false;
		else
			return true;
	}
	//else if ((iter) % pest_scenario.get_pestpp_options().get_opt_recalc_fosm_every() == 0)
    else if (chance_schedule[iter])
		return true;
	return false;
}

void Constraints::initialize_chance_schedule(ofstream& frec)
{
    stringstream ss;
    chance_schedule.clear();

    string fname = pest_scenario.get_pestpp_options().get_opt_chance_schedule();
    string line;
    vector<string> tokens;
    int lcount = 0, gen;
    bool should_eval = false;
    int recalc_every = pest_scenario.get_pestpp_options().get_opt_recalc_fosm_every();
    for (int i=0;i<max(1,pest_scenario.get_control_info().noptmax+1);i++)
        if (i%recalc_every == 0)
            chance_schedule[i] = true;
        else
        chance_schedule[i] = false;
    chance_schedule[0] = true;
    if (fname.size() > 0)
    {
        ss.str("");
        ss << "...reading population schedule from file '" << fname << "' " << endl;
        frec << ss.str();
        cout << ss.str();
        ifstream in(fname);
        if (in.bad())
        {
            throw runtime_error("error opening opt_chance_schedule file '"+fname+"'");
        }
        while (getline(in,line))
        {
            lcount++;
            tokens.clear();
            pest_utils::tokenize(line,tokens,"\t ,");
            if (tokens.size() < 2)
            {
                ss.str("");
                ss << "opt_chance_schedule file '" << fname << "' line " << lcount << " needs at least two entries";
                throw runtime_error(ss.str());
            }
            try
            {
                gen = stoi(tokens[0]);
            }
            catch (...)
            {
                ss.str("");
                ss << "error casting '" << tokens[0] << "' to generation integer on line " << lcount << " in opt_chance_schedule file";
                throw runtime_error(ss.str());
            }

            try
            {
                should_eval = pest_utils::parse_string_arg_to_bool(tokens[1]);
            }
            catch (...)
            {
                ss.str("");
                ss << "error casting '" << tokens[1] << "' to bool on line " << lcount << " in opt_chance_schedule file";
                throw runtime_error(ss.str());
            }
            chance_schedule[gen] = should_eval;
        }
        in.close();
    }

    if (!chance_schedule[0])
    {
        ss.str("");
        ss << "...chances must be evaluated during the first generation, resetting chance_schedule" << endl;
        frec << ss.str();
        cout << ss.str();
        chance_schedule[0] = true;
    }

    frec << "...chance schedule: generation,should_eval_chances:" << endl;
    for (int i=0;i<pest_scenario.get_control_info().noptmax;i++)
    {
        frec << "...   " << i << ", " << chance_schedule.at(i) << endl;
    }
}

void Constraints::postsolve_obs_constraints_report(Observations& old_obs, Observations& new_obs, string tag, int iter, 
									map<string,string> status_map, map<string,double> price_map)
{
	/* write out constraint info after the current iteration solution process is over to the rec file.  
	
	This really only applies to the simplex solution because we can, with the linear assumption, 
	project what the resulting constraint values should be*/

	int width = 20;
	for(auto& n : ctl_ord_obs_constraint_names)
	    width = max(width,(int)n.size());
	width = width + 2;
	ofstream& f_rec = file_mgr_ptr->rec_ofstream();
	f_rec << endl << endl << "     " << tag << " constraint/objective information at end of iteration " << iter << endl << endl;
	f_rec << setw(width) << left << "name" << right << setw(15) << "sense" << setw(15) << "required" << setw(25);
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
		f_rec << setw(width) << left << name;
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
	    int width = 20;
	    for (auto& n : ctl_ord_pi_constraint_names)
	        width = max(width,(int)n.size());
	    width = width + 2;
		f_rec << endl << endl << " --- prior information constraint summary at end of iteration " << iter << " --- " << endl << endl;
		f_rec << setw(width) << left << "name" << right << setw(15) << "sense" << setw(15) << "required" << setw(25);
		if (status_map.size() > 0)
			f_rec << "simplex status";
		if (price_map.size() > 0)
			f_rec << setw(15) << "price";
		f_rec  << setw(15) << "current" << setw(15) << "residual";
		f_rec << setw(15) << "new" << setw(15) << "residual" << endl;
		int i = 0;
		for (auto &name : ctl_ord_pi_constraint_names)
		{
			PriorInformationRec pi_rec = constraints_pi.get_pi_rec(name);
			pair<double,double> cur_sim_resid = pi_rec.calc_sim_and_resid(old_pars);
			pair<double,double> new_sim_resid = pi_rec.calc_sim_and_resid(new_pars);
			f_rec << setw(width) << left << name;
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
	RunManagerAbstract* run_mgr_ptr, bool drop_fails, bool debug_fail)
{
	if (_stack_pe_run_map.size() == 0)
		pfm.log_event("process_stack_runs() was passed an empty run map");
	else
		stack_runs_processed = true;
	ObservationEnsemble _stack_oe(stack_oe);//copy
	vector<int> failed_runs = _stack_oe.update_from_runs(_stack_pe_run_map, run_mgr_ptr);
	if (debug_fail)
	{

		failed_runs.push_back(0);
		cout << "ies_debug_fail_subset = true, failing first stack realization" << endl;
		file_mgr_ptr->rec_ofstream() << "debug = true, failing first stack realization for " << real_name <<  endl;
	}
	stringstream ss;
	if (failed_runs.size() > 0)
	{
		ss.str("");
		ss << "WARNING: " << failed_runs.size() << " stack runs failed for realization " << real_name;
		
		if ((drop_fails) && (failed_runs.size() < _stack_oe.shape().first))
		{
			_stack_oe.drop_rows(failed_runs);
			ss << ", dropped...";
		}
		pfm.log_event(ss.str());
		file_mgr_ptr->rec_ofstream() << ss.str() << endl;
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

	if (pest_scenario.get_pestpp_options().get_save_binary())
	{
        if (pest_scenario.get_pestpp_options().get_save_dense())
        {
            ss << ".bin";
            _stack_oe.to_dense(ss.str());
        }
        else {
            ss << ".jcb";
            _stack_oe.to_binary(ss.str());
        }
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

	if (pest_scenario.get_pestpp_options().get_save_binary())
	{
        if (pest_scenario.get_pestpp_options().get_save_dense())
        {
            ss << ".bin";
            _stack_pe.to_dense(ss.str());
        }
        else {
            ss << ".jcb";
            _stack_pe.to_binary(ss.str());
        }

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
		cout << "...processing nested stack runs" << endl;
		ObservationEnsemble nested_oe;
		vector<string> keep_nested_names;
		stack_oe_map.clear();
		int i = 0;
		stringstream ss;
		

		//work out what var names for the stack_oe we dont need so we can drop them!
		vector<string> drop_names;
		set<string> s_obs_constraint_names(ctl_ord_obs_constraint_names.begin(), ctl_ord_obs_constraint_names.end());
		vector<string> names1, names2;
		for (auto name : stack_oe.get_var_names())
		{
			if (s_obs_constraint_names.find(name) == s_obs_constraint_names.end())
				drop_names.push_back(name);
		}
		vector<string> drop_members;
		bool test_failed = false;
		bool should_fail = false;
	
		for (auto& real_info : population_stack_pe_run_map)
		{
			pfm.log_event("processing stack runs for realization " + real_info.first);
			should_fail = false;
			if ((!test_failed) && (pest_scenario.get_pestpp_options().get_ies_debug_fail_subset()))
			{
				should_fail = true;
				test_failed = true;
			}
			stack_info = process_stack_runs(real_info.first, iter, real_info.second, run_mgr_ptr, true, should_fail);

			//test all but one stack runs failed
			if ((!test_failed) && (pest_scenario.get_pestpp_options().get_ies_debug_fail_remainder()))
			{
				for (int ii = 1; ii < real_info.second.size(); ii++)
					stack_info.first.push_back(ii);
				test_failed = true;
				cout << "ies_debug_fail_remainder = true, failing full stack for member " << real_info.first << endl;
				file_mgr_ptr->rec_ofstream() << "ies_debug_fail_remainder = true, failing most of stack for member " << real_info.first << endl;
			}
			if (stack_info.first.size() > 0)
			{
				if (stack_info.first.size() == stack_pe.shape().first)
				{
					ss.str("");
					ss << "WARNING: all stack runs failed for population member '" << real_info.first << 
						", removing this member from chance calcs" << endl;
					file_mgr_ptr->rec_ofstream() << ss.str();
					cout << ss.str();
					drop_members.push_back(real_info.first);
					continue;
				}
				else
				{
					//nothing to do here since we have only stored the dec var values, not the full stack of par values.
				}
			}

			names1 = stack_info.second.get_real_names();
			names2.clear();
			for (auto n : names1)
				names2.push_back(n + NESTED_STACK_NAME_SEP + real_info.first);
			stack_info.second.set_real_names(names2,true);
			if (i == 0)
			{
				//stack_info.second.to_csv_by_reals(fstack, true);
				nested_oe = stack_info.second;
			}
			else
			{
				//stack_info.second.to_csv_by_reals(fstack, false);
				nested_oe.append_other_rows(stack_info.second, true);
			}
			stack_info.second.set_real_names(names1,true);
			//save_oe_stack(iter, real_info.first, stack_info.second);
			stack_info.second.drop_cols(drop_names);
			stack_oe_map[real_info.first] = stack_info.second;
			cout << real_info.first << "\r" << flush;
			i++;
			
		}
		for (auto d : drop_members)
		{
			population_stack_pe_run_map.erase(d);
			stack_pe_map.erase(d);
		}
		//update nested pe for any failed runs
		nested_pe.keep_rows(nested_oe.get_real_names());
		nested_pe.reset_org_real_names();
		if (pest_scenario.get_pestpp_options().get_save_binary())
		{
            if (pest_scenario.get_pestpp_options().get_save_dense())
            {
                ss.str("");
                ss << file_mgr_ptr->get_base_filename() << "." << iter << ".nested.obs_stack.bin";
                nested_oe.to_dense(ss.str());
                ss.str("");
                ss << file_mgr_ptr->get_base_filename() << "." << iter << ".nested.par_stack.bin";
                nested_pe.to_dense(ss.str());

            }
            else {
                ss.str("");
                ss << file_mgr_ptr->get_base_filename() << "." << iter << ".nested.obs_stack.jcb";
                nested_oe.to_binary(ss.str());
                ss.str("");
                ss << file_mgr_ptr->get_base_filename() << "." << iter << ".nested.par_stack.jcb";
                nested_pe.to_binary(ss.str());
            }
			
		}
		else
		{
			ss.str("");
			ss << file_mgr_ptr->get_base_filename() << "." << iter << ".nested.obs_stack.csv";
			nested_oe.to_csv(ss.str());
			ss.str("");
			ss << file_mgr_ptr->get_base_filename() << "." << iter << ".nested.par_stack.csv";
			nested_pe.to_csv(ss.str());

		}
		//clear these two out
		nested_oe = ObservationEnsemble(&pest_scenario);
		nested_pe = ParameterEnsemble(&pest_scenario);
		nested_stack_stdev_summary(stack_oe_map);
		

	}
	else
	{
		pfm.log_event("processing stack runs");
		stack_info = process_stack_runs("", iter, stack_pe_run_map, run_mgr_ptr, true);
		
		if (stack_info.first.size() > 0)
		{
			if (stack_info.first.size() == stack_pe.shape().first)
			{
				throw_constraints_error("all stack runs failed, cannot continue");
			}
			file_mgr_ptr->rec_ofstream() << "dropping " << stack_info.first.size() << " failed realizations for par stack";
			pfm.log_event("dropping failed realizations from par stack");
			stack_pe.drop_rows(stack_info.first);
		}
		save_oe_stack(iter, "", stack_info.second);
		save_pe_stack(iter, "", stack_pe);
		stack_oe = stack_info.second;
	}
}

void Constraints::nested_stack_stdev_summary(map<string, ObservationEnsemble>& _stack_oe_map)
{
	Eigen::VectorXd cvals(ctl_ord_obs_constraint_names.size());
	pair < map<string, double>, map<string, double>> mm;
	map <string, vector<double>> nested_std_map;
	int mxlen = 11;
	for (auto& cname : ctl_ord_obs_constraint_names)
	{
		nested_std_map[cname] = vector<double>();
		mxlen = max(mxlen, (int)cname.size());
	}

	//string cname;
	map<string, double> mean_map, std_map;
	for (auto& m : _stack_oe_map)
	{
		mean_map.clear();
		std_map.clear();
		m.second.fill_moment_maps(mean_map, std_map);
		for (auto& cname : ctl_ord_obs_constraint_names)
		{
			nested_std_map[cname].push_back(std_map[cname]);
		}
	}
	
	stringstream ss;
	ss.str("");
	ss << "   Change point = 'all' constraint standard deviation summary" << endl;
	ss << left << setw(mxlen+1) << "constraint" << setw(9) << "count" << setw(13) << "mean" << setw(13) << "stdev" << setw(13) << "min" << setw(13) << "max" << endl;
	double mean, st,mn,mx;
	for (auto& cname : ctl_ord_obs_constraint_names)
	{
		cvals = stlvec_2_eigenvec(nested_std_map[cname]);
		mean = cvals.mean();
		mn = cvals.minCoeff();
		mx = cvals.maxCoeff();
		cvals = cvals.array() - mn;
		st = cvals.array().pow(2).sum();
		if (st > 0)
			st = sqrt(st);
		ss << setw(mxlen) << cname << " ";
		ss << setw(8) << cvals.size() << " ";
		ss << setw(12) << mean << " ";
		ss << setw(12) << st << " ";
		ss << setw(12) << mn << " ";
		ss << setw(12) << mx << " ";
		ss << endl;
	}
	ss << "   note: the above summary is an indication of how the chance estimates interact with" << endl;
	ss << "         the decision variables. A large variation in the standard deviation across" << endl;
	ss << "         the chance point locations indicates that there is non-trivial interaction" << endl;
	ss << "         between the chance estimates and the decision variables." << endl;
	cout << ss.str() << endl;
	file_mgr_ptr->rec_ofstream() << ss.str() << endl;
 	
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
		pfm.log_event("reading FOSM-based parameter perturbation runs into JCO");
		ParamTransformSeq par_trans = pest_scenario.get_base_par_tran_seq();
		bool success = jco.process_runs(par_trans, pest_scenario.get_base_group_info(), *run_mgr_ptr,
			*null_prior, false, false);
		if (!success)
			throw_constraints_error("error processing FOSM JCO matrix runs ", jco.get_failed_parameter_names());
		ss.str("");
		ss << iter << ".fosm.jcb";
		jco.save(ss.str());
	}
	//otherwise, process some stack runs
	else
	{
		if (stack_runs_processed)
			return;
		process_stack_runs(run_mgr_ptr,iter);

	}

}

void Constraints::add_runs(int iter, Parameters& current_pars, Observations& current_obs, RunManagerAbstract* run_mgr_ptr)
{
	/* using the passed run mgr pointer, queue up chance runs
	
	*/
	if (!use_chance)
		return;
	//for fosm, we need to queue up parameter perturbation runs, centered at the current dec var values
	if (use_fosm)
	{
		pfm.log_event("building FOSM-based parameter perturbation runs");
		ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();
		set<string> out_of_bounds;
		//the  second-to-last true arg is to run the parvals listed in _current_pars_and_dec_vars as the "base" run
		//this is to make sure we get good perturbations at the current point in par/dec var space.
		//but this means we might have to run the base run twice at least for the SLP process.
		//the way around this is to just add the runs all at once within the SLP add run
		//bit using Constraints::get_fosm_par_names()
		//the last false says not to reinitialize the run mgr since the calling process may have also
		//added runs
		bool success = jco.build_runs(current_pars, current_obs, adj_par_names, pts,
			pest_scenario.get_base_group_info(), pest_scenario.get_ctl_parameter_info(),
			*run_mgr_ptr, out_of_bounds, false, true, false);
		if (!success)
		{
			const set<string> failed = jco.get_failed_parameter_names();
			throw_constraints_error("failed to calc derivatives for the following FOSM parameters: ", failed);
		}
		cout << "...adding " << jco.get_par_run_map().size() << " model runs for FOSM-based chance constraints" << endl;

	}
	//for stacks, we need to queue up the stack realizations, but replace the dec var entries in each realization
	//with the current dec var values.
	else
	{
		stack_pe_run_map = add_stack_runs(iter, stack_pe, current_pars, current_obs, run_mgr_ptr);
		cout << "...adding " << stack_pe_run_map.size() << " model runs for stack-based chance constraints" << endl;
		stack_runs_processed = false;
		
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
	stack_pe_map.clear();
	pfm.log_event("queuing up nested-sets of chance runs for decision-variable ensemble ");
	
	Eigen::VectorXd par_vec;
	vector<string> par_names = current_pe.get_var_names();
	Parameters real_pars = pest_scenario.get_ctl_parameters();
	pest_scenario.get_base_par_tran_seq().ctl2numeric_ip(real_pars);
	current_pe.transform_ip(ParameterEnsemble::transStatus::NUM);
	ofstream& f_rec = file_mgr_ptr->rec_ofstream();
	f_rec << "queuing up nested sets of chance runs ('opt_chance_points' == 'all') for decision-variable ensemble " << endl;
	cout << "queuing up nested sets of chance runs for decision-variable ensemble " << endl;
	map<int, int> _stack_pe_run_map;
	int count = 0;
	int i = 0;
	stringstream ss;
	
	vector<string> names1, names2;
	for (auto real_info : current_pe.get_real_map())
	{

		par_vec = current_pe.get_real_vector(real_info.first);
		real_pars.update_without_clear(par_names, par_vec);
		_stack_pe_run_map = add_stack_runs(iter, stack_pe, real_pars, current_obs, run_mgr_ptr);
		population_stack_pe_run_map[real_info.first] = _stack_pe_run_map;
		
		stack_pe_map[real_info.first] = real_pars;
		names1 = stack_pe.get_real_names();
		names2.clear();
		for (auto n : names1)
			names2.push_back(n + NESTED_STACK_NAME_SEP + real_info.first);
		stack_pe.set_real_names(names2,true);
		//TODO: add a flag to decide if saving
		if (i == 0)
		{
			//stack_pe.to_csv_by_reals(fstack, true);
			nested_pe = stack_pe;
		}
		else
		{
			//stack_pe.to_csv_by_reals(fstack, false);
			nested_pe.append_other_rows(stack_pe,true);
		}
		stack_pe.set_real_names(names1,true);
		//save_pe_stack(iter, real_info.first, stack_pe);
		f_rec << "...added " << stack_pe.shape().first << " runs for decision variable solution '" << real_info.first << "'" << endl;
		count = count + _stack_pe_run_map.size();
		cout << real_info.first << "\r" << flush;
		i++;
	}

	//save this as bin here - it will be resaved after runs are processed...
	if (pest_scenario.get_pestpp_options().get_save_binary())
    {
        if (pest_scenario.get_pestpp_options().get_save_dense())
        {
            ss.str("");
            ss << file_mgr_ptr->get_base_filename() << "." << iter << ".nested.par_stack.bin";
            nested_pe.to_dense(ss.str());
        }
        else
        {
            ss.str("");
            ss << file_mgr_ptr->get_base_filename() << "." << iter << ".nested.par_stack.jcb";
            nested_pe.to_binary(ss.str());
        }
    }
    else
    {
        ss.str("");
        ss << file_mgr_ptr->get_base_filename() << "." << iter << ".nested.par_stack.csv";
        nested_pe.to_csv(ss.str());

    }

	cout << "...adding " << count << " runs nested stack-based chance constraints" << endl;
	stack_runs_processed = false;
	//reset stack_oe to use the same real names as stack_pe
	stack_oe.reserve(names1, stack_oe.get_var_names());
	stack_oe.set_real_names(names1, true);
	
	
}

map<int, int> Constraints::add_stack_runs(int iter, ParameterEnsemble& _stack_pe, Parameters& current_pars, Observations& current_obs, RunManagerAbstract* run_mgr_ptr)
{
	//update _stack_pe parameter values for decision variables using current_pars
	int num_reals = _stack_pe.shape().first;
	_stack_pe.update_var_map();
	for (auto dname : dec_var_names)
	{
		Eigen::VectorXd dvec(num_reals);
		dvec.setConstant(current_pars.get_rec(dname));
		_stack_pe.replace_col(dname, dvec,false);
	}
	pfm.log_event("building stack-based parameter runs");
	map<int,int> _stack_pe_run_map = _stack_pe.add_runs(run_mgr_ptr);
	return _stack_pe_run_map;
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
	tol is a percent-based tolerance to account for constraints that are very near their required (rhs) value

	*/
	double sim_val, obs_val, scaled_diff;
	map<string, double> unsatisfied;
	
	for (auto& name : ctl_ord_pi_constraint_names)
	{
		PriorInformationRec pi_rec = constraints_pi.get_pi_rec(name);
		//pair<double, double> cur_sim_resid = pi_rec.calc_residual_and_sim_val(*current_pars_and_dec_vars_ptr);
		pair<double, double> new_sim_resid = pi_rec.calc_sim_and_resid(par_and_dec_vars);
		//check for invalid pi constraints
		sim_val = new_sim_resid.first;
		obs_val = pi_rec.get_obs_value();
		if (obs_val != 0)
			scaled_diff = abs((obs_val - sim_val) / obs_val);
		else
			scaled_diff = abs(obs_val - sim_val);
		if ((constraint_sense_map[name] == ConstraintSense::less_than) && (sim_val > obs_val) && (scaled_diff > tol))
			unsatisfied[name] = sim_val - obs_val;
		else if ((constraint_sense_map[name] == ConstraintSense::greater_than) && (sim_val < obs_val) && (scaled_diff > tol))
			unsatisfied[name] = obs_val - sim_val;
		else if ((constraint_sense_map[name] == ConstraintSense::equal_to) && (sim_val != obs_val) && (scaled_diff > tol))
			unsatisfied[name] = abs(sim_val - obs_val);
	}
	return unsatisfied;
}

map<string, map<string, double>> Constraints::get_ensemble_violations_map(ParameterEnsemble& pe, ObservationEnsemble& oe, double tol, bool include_weight)
{
	//make sure pe and oe share realizations
	vector<string> pe_names = pe.get_real_names();
	vector<string> oe_names = oe.get_real_names();
	if (pe_names.size() != oe_names.size())
		throw_constraints_error("get_ensemble_violations_map(): pe reals != oe reals");
	for (int i=0;i<pe_names.size();i++)
		if (pe_names[i] != oe_names[i])
			throw_constraints_error("get_ensemble_violations_map(): pe reals != oe reals");
	
	Observations obs = pest_scenario.get_ctl_observations();
	Eigen::VectorXd v;
	vector<string> vnames = oe.get_var_names();
	map<string, double> vmap;
	map<string, map<string, double>> violations;
	for (auto& name : oe_names)
	{
		v = oe.get_real_vector(name);
		obs.update_without_clear(vnames, v);
		vmap = get_unsatified_obs_constraints(obs,tol,true,include_weight);
		violations[name] = vmap;
	}

	Parameters pars = pest_scenario.get_ctl_parameters();
	ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();
	pts.ctl2numeric_ip(pars);
	pe.transform_ip(ParameterEnsemble::transStatus::NUM);
	vnames = pe.get_var_names();
	for (auto& name : pe_names)
	{
		v = pe.get_real_vector(name);
		pars.update_without_clear(vnames, v);
		vmap = get_unsatified_pi_constraints(pars);
		violations[name].insert(vmap.begin(), vmap.end());
	}
	return violations;
}


map<string, double> Constraints::get_unsatified_obs_constraints(Observations& constraints_sim, double tol, bool do_shift, bool include_weight)
{
	/* get a map of name, distance for each of the obs-based (e.g. model-based) constraints that are not satisfied in the constraint_obs container.
	tol is a percent-based tolerance to account for constraints that are very near their required (rhs) value

	*/
	double sim_val, obs_val, scaled_diff;
	map<string, double> unsatisfied;
	Observations _constraints_sim(constraints_sim);
	if ((do_shift) && (use_chance))
		_constraints_sim = get_chance_shifted_constraints(constraints_sim);
	ObservationInfo oi = pest_scenario.get_ctl_observation_info();
	double weight;
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
		if (obs_val != 0)
			scaled_diff = abs((obs_val - sim_val) / obs_val);
		else
			scaled_diff = abs(obs_val - sim_val);
        weight = 1.0;
        if (include_weight)
            weight = oi.get_weight(name);
		//check for invalid obs constraints (e.g. satisfied)
		if ((constraint_sense_map[name] == ConstraintSense::less_than) && (sim_val > obs_val) && (scaled_diff > tol))
			unsatisfied[name] = weight * (sim_val - obs_val);
		else if ((constraint_sense_map[name] == ConstraintSense::greater_than) && (sim_val < obs_val) && (scaled_diff > tol))
			unsatisfied[name] = weight * (obs_val - sim_val);
		else if ((constraint_sense_map[name] == ConstraintSense::equal_to) && (sim_val != obs_val) && (scaled_diff > tol))
			unsatisfied[name] = weight * abs(sim_val - obs_val);
	}
	return unsatisfied;

}

map<string, double> Constraints::get_constraint_map(Parameters& par_and_dec_vars, Observations& constraints_sim, bool do_shift)
{
	double sim_val, obs_val;
	map<string, double> constraint_map;
	Observations _constraints_sim(constraints_sim);
	if ((do_shift) && (use_chance))
		_constraints_sim = get_chance_shifted_constraints(constraints_sim);
	//for (int i = 0; i < num_obs_constraints(); ++i)
	for (auto& name : ctl_ord_obs_constraint_names)
	{
		//string name = ctl_ord_obs_constraint_names[i];
		sim_val = constraints_sim[name];
		if ((do_shift) && (use_chance))
		{
			double offset = post_constraint_offset[name];
			sim_val += offset;
		}
		obs_val = constraints_obs[name];
		constraint_map[name] = sim_val - obs_val;
	}

	for (auto& name : ctl_ord_pi_constraint_names)
	{
		PriorInformationRec pi_rec = constraints_pi.get_pi_rec(name);
		//pair<double, double> cur_sim_resid = pi_rec.calc_residual_and_sim_val(*current_pars_and_dec_vars_ptr);
		pair<double, double> new_sim_resid = pi_rec.calc_sim_and_resid(par_and_dec_vars);
		//check for invalid pi constraints
		sim_val = new_sim_resid.first;
		obs_val = pi_rec.get_obs_value();
		constraint_map[name] = sim_val - obs_val;
		
	}
	return constraint_map;
}

Mat Constraints::get_working_set_constraint_matrix(Parameters& par_and_dec_vars, Observations& constraints_sim, ParameterEnsemble& dv, ObservationEnsemble& oe, bool do_shift, double working_set_tol)
{
    pair<vector<string>,vector<string>> working_set = get_working_set(par_and_dec_vars,constraints_sim,do_shift,working_set_tol);
    Mat mat;
    if (working_set.first.size() > 0) {
        Covariance cov = dv.get_empirical_cov_matrices(file_mgr_ptr).second;
        Eigen::MatrixXd delta_dv = *cov.inv().e_ptr() * dv.get_eigen_anomalies().transpose();

        cov = oe.get_empirical_cov_matrices(file_mgr_ptr).second;
        Eigen::MatrixXd delta_oe = *cov.inv().e_ptr() * oe.get_eigen_anomalies().transpose();
        //todo: pseudo inv for delta_dv - will almost certainly be singular for large problems...
        Eigen::MatrixXd s, V, U;
        Eigen::BDCSVD<Eigen::MatrixXd> svd_fac(delta_dv, Eigen::DecompositionOptions::ComputeFullU |
                                                         Eigen::DecompositionOptions::ComputeFullV);
        s = svd_fac.singularValues();
        U = svd_fac.matrixU();
        V = svd_fac.matrixV();
        s.col(0).asDiagonal().inverse();
        Eigen::MatrixXd full_s_inv(V.rows(), U.cols());
        full_s_inv.setZero();
        for (int i = 0; i < s.size(); i++)
            full_s_inv(i, i) = s(i);
        delta_dv = V.transpose() * full_s_inv * U;
        //delta_oe.transposeInPlace();
        //delta_dv.transposeInPlace();
        Eigen::MatrixXd approx_jco = delta_oe * delta_dv;
        delta_dv.resize(0, 0);
        delta_oe.resize(0, 0);
        V.resize(0, 0);
        U.resize(0, 0);
        full_s_inv.resize(0, 0);
        Eigen::MatrixXd working_mat(working_set.first.size(), dv.shape().second);
        oe.update_var_map();
        map<string, int> vmap = oe.get_var_map();
        int i = 0;
        for (auto &n : working_set.first) {
            working_mat.row(i) = approx_jco.row(vmap[n]);
            i++;
        }
        mat = Mat(working_set.first, dv.get_var_names(), working_mat.sparseView());
        //return Mat(working_set, dv.get_var_names(), working_mat.sparseView());
    }
    if (working_set.second.size() > 0)
    {
        //deal with pi constraints in the working set
        augment_constraint_mat_with_pi(mat,working_set.second);
    }
    return Mat();
}

Mat Constraints::get_working_set_constraint_matrix(Parameters& par_and_dec_vars, Observations& constraints_sim, const Jacobian_1to1& _jco, bool do_shift, double working_set_tol)
{
	pair<vector<string>,vector<string>> working_set = get_working_set(par_and_dec_vars,constraints_sim,do_shift,working_set_tol);
	Mat mat;
    if (working_set.first.size() > 0) {
        Eigen::SparseMatrix<double> t = _jco.get_matrix(working_set.first, dec_var_names);
        Mat(working_set.first, dec_var_names, t);
    }
	if (working_set.second.size() > 0) {
        augment_constraint_mat_with_pi(mat,working_set.second);


    }
	return mat;
}

void Constraints::augment_constraint_mat_with_pi(Mat& mat, vector<string>& pi_names)
{
    

}

pair<vector<string>,vector<string>> Constraints::get_working_set(Parameters& par_and_dec_vars, Observations& constraints_sim, bool do_shift, double working_set_tol) {
    map<string, double> constraint_map = get_constraint_map(par_and_dec_vars, constraints_sim, do_shift);
    vector<string> working_set,working_set_pi;
    for (auto &name : ctl_ord_obs_constraint_names) {
        if (constraint_sense_map[name] == ConstraintSense::equal_to)
            working_set.push_back(name);

        else if (abs(constraint_map[name]) < working_set_tol) {
            working_set.push_back(name);
        }
    }
    for (auto &name : ctl_ord_pi_constraint_names) {
        if (constraint_sense_map[name] == ConstraintSense::equal_to)
            working_set_pi.push_back(name);

        else if (abs(constraint_map[name]) < working_set_tol) {
            working_set_pi.push_back(name);
        }
    }
    return pair<vector<string>,vector<string>>(working_set,working_set_pi);
}

int Constraints::get_num_nz_pi_constraint_elements()
{
	/* work out how many elements are in the PI constraints.  This is needed for dimensioning the 
	simplex response matrix*/
	int num = 0;

	for (auto &pi_name : ctl_ord_pi_constraint_names)
	{
		num += constraints_pi.get_pi_rec(pi_name).get_atom_factors().size();
	}
	return num;
}



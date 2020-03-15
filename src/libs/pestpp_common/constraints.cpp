
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

Constraints::Constraints(Pest& _pest_scenario, FileManager* _file_mgr_ptr, OutputFileWriter& _of_wr, PerformanceLog& _pfm)
	:pest_scenario(_pest_scenario), file_mgr_ptr(_file_mgr_ptr), of_wr(_of_wr), pfm(_pfm), jco(*_file_mgr_ptr, _of_wr)
{
	;
}


void Constraints::initialize(vector<string>& ctl_ord_dec_var_names, Parameters* _current_pars_and_dec_vars_ptr,
	Observations* _current_constraints_sim_ptr, double _dbl_max)
{
	current_constraints_sim_ptr = _current_constraints_sim_ptr;
	current_pars_and_dec_vars_ptr = _current_pars_and_dec_vars_ptr;
	stack_size = pest_scenario.get_pestpp_options().get_opt_stack_size();
	use_fosm = true;
	if (stack_size > 0)
		use_fosm = false;
	stack_pe.set_pest_scenario(&pest_scenario);
	stack_pe.set_rand_gen(&rand_gen);
	dbl_max = _dbl_max;
	rand_gen = std::mt19937(pest_scenario.get_pestpp_options().get_random_seed());
	ctl_ord_obs_constraint_names.clear();
	ctl_ord_pi_constraint_names;
	constraint_sense_map.clear();
	constraint_sense_name.clear();
	nz_obs_names.clear();
	adj_par_names.clear();
	prior_const_var.clear();
	post_const_var.clear();

	ofstream& f_rec = file_mgr_ptr->rec_ofstream();
	pfm.log_event("initializing constraints");
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
	//constraint_lb = new double[num_constraints()];
	//constraint_ub = new double[num_constraints()];


	//------------------------------------------
	//  ---  chance constratints  ---
	//------------------------------------------
	risk = pest_scenario.get_pestpp_options().get_opt_risk();
	if (risk != 0.5)
	{
		pfm.log_event("initializing chance constraints");
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
		else if (use_fosm)
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

			f_rec << "  using FOSM-based chance constraints with risk = " << risk << endl;

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
				obscov.from_observation_weights(nz_obs_names, oi, vector<string>(), null_prior);
			}
		}
		else
		{
			string par_csv = pest_scenario.get_pestpp_options().get_ies_par_csv();
			
			if (par_csv.size() == 0)
			{
			
				f_rec << "drawing " << stack_size << "stack realizations" << endl;
				pfm.log_event("loading parcov");
				parcov.try_from(pest_scenario, *file_mgr_ptr);
				pfm.log_event("drawing stack realizations");
				stack_pe.draw(stack_size, pest_scenario.get_ctl_parameters(), parcov, &pfm, 
					pest_scenario.get_pestpp_options().get_ies_verbose_level());
			}
			else
			{
				string par_ext = pest_utils::lower_cp(par_csv).substr(par_csv.size() - 3, par_csv.size());
				pfm.log_event("processing par csv " + par_csv);
				if (par_ext.compare("csv") == 0)
				{
					pfm.log_event("loading par ensemble from csv file: " + par_csv);
					try
					{
						stack_pe.from_csv(par_csv);
					}
					catch (const exception& e)
					{

						throw_constraints_error("error processing par csv: " + string(e.what()));
					}
					catch (...)
					{
						throw_constraints_error(string("error processing par csv"));
					}
				}
				else if ((par_ext.compare("jcb") == 0) || (par_ext.compare("jco") == 0))
				{
					pfm.log_event("loading par ensemble from binary file+ " + par_csv);
					try
					{
						stack_pe.from_binary(par_csv);
					}
					catch (const exception& e)
					{

						throw_constraints_error("error processing par ensemble binary file: " + string(e.what()));
					}
					catch (...)
					{
						throw_constraints_error("error processing par ensemble binary file");
					}
				}
				else
				{
					throw_constraints_error("unrecognized par ensemble extension, looking for csv, jcb, jco.  found: " + par_ext);
				}
				if (stack_pe.shape().first > stack_size)
				{
					vector<int> drop_rows;
					for (int i = stack_size - 1; i < stack_pe.shape().first; i++)
						drop_rows.push_back(i);
					f_rec << "droppping " << drop_rows.size() << "realizations from stack b/c of ++opt_stack_size req" << endl;
					stack_pe.drop_rows(drop_rows);
				}
			}
			string filename = file_mgr_ptr->get_base_filename() + ".0.stack";
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
			f_rec << "saved initial stack to " << filename << endl;
		}
	}
	else use_chance = false;
}

void Constraints::initial_report()
{
	ofstream& f_rec = file_mgr_ptr->rec_ofstream();

	f_rec << endl << "  ---  observation constraints in SLP  ---  " << endl;
	f_rec << setw(20) << "name" << setw(20) << "sense" << setw(20) << "value" << endl;
	for (auto& name : ctl_ord_obs_constraint_names)
	{
		f_rec << setw(20) << left << name;
		f_rec << setw(20) << constraint_sense_name[name];
		f_rec << setw(20) << constraints_obs.get_rec(name) << endl;
	}

	if (num_pi_constraints() > 0)
	{
		f_rec << endl << "  ---  prior information constraints in SLP  ---  " << endl;
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
		f_rec << endl << endl << "  ---  chance constraint FOSM information  ---  " << endl;
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
		if (num_nz_obs() == 0)
		{
			f_rec << endl << endl << "  ---  Note: No nonzero weight observations found." << endl;
			f_rec << "           Prior constraint uncertainty will be used in chance constraint calculations" << endl;
		}
		else
		{
			f_rec << "  non-zero weight observations used for conditioning in FOSM calculations: " << endl;
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
}


void Constraints::update_chance_offsets()
{
	if ((!use_chance) || (std_weights))
		return;

	if (use_fosm)
	{
		ofstream& f_rec = file_mgr_ptr->rec_ofstream();
		cout << "  ---  calculating FOSM-based chance constraint components  ---  " << endl;
		f_rec << "  ---  calculating FOSM-based chance constraint components  ---  " << endl;

		//the rows of the fosm jacobian include nonzero weight obs (for schur comp)
		//plus the names of the names of constraints, which get treated as forecasts
		vector<string> fosm_row_names(nz_obs_names);
		fosm_row_names.insert(fosm_row_names.end(), ctl_ord_obs_constraint_names.begin(), ctl_ord_obs_constraint_names.end());

		//extract the part of the full jco we need for fosm
		Eigen::SparseMatrix<double> fosm_mat = jco.get_matrix(fosm_row_names, adj_par_names);

		if (fosm_mat.size() == 0)
			throw_constraints_error("FOSM-based chance constraint-to-parameter vectors are all zeros");

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
		cout << "  ---   done with FOSM-based chance constraint calculations  ---  " << endl << endl;
		f_rec << "  ---   done with FOSM-based chance constraint calculations  ---  " << endl << endl;
	}
	else
	{
		throw_constraints_error("constraints::get_chance() not implemented for non-fosm chances");
	}
}

double Constraints::get_max_constraint_change(Observations& upgrade_obs)
{
	double max_abs_constraint_change = -1.0E+10;
	double max_abs_constraint_val = -1.0E+10;
	double val, diff;
	for (auto& name : ctl_ord_obs_constraint_names)
	{
		val = current_constraints_sim_ptr->get_rec(name);
		diff = abs(current_constraints_sim_ptr->get_rec(name) - upgrade_obs[name]);
		max_abs_constraint_change = (diff > max_abs_constraint_change) ? diff : max_abs_constraint_change;
		max_abs_constraint_val = (val > max_abs_constraint_val) ? val : max_abs_constraint_val;
	}
	max_abs_constraint_change /= max(max_abs_constraint_val, 1.0);
	return max_abs_constraint_change;
}

Observations Constraints::get_chance_shifted_constraints()
{
	ofstream& f_rec = file_mgr_ptr->rec_ofstream();
	prior_constraint_offset.clear();
	prior_constraint_stdev.clear();
	post_constraint_offset.clear();
	post_constraint_stdev.clear();
	//if ((!std_weights) && ((slp_iter == 1) || ((slp_iter + 1) % pest_scenario.get_pestpp_options().get_opt_recalc_fosm_every() == 0)))
	
	//work out the offset for each constraint
	//and set the values in the constraints_fosm Obseravtions
	double new_constraint_val, old_constraint_val, required_val;
	double pr_offset, pt_offset;
	//constraints_fosm.clear();
	Observations iter_fosm;

	//map<string,double> out_of_bounds;
	for (auto& name : ctl_ord_obs_constraint_names)
	{
		prior_constraint_stdev[name] = sqrt(prior_const_var[name]);
		post_constraint_stdev[name] = sqrt(post_const_var[name]);
		pr_offset = probit_val * prior_constraint_stdev[name];
		pt_offset = probit_val * post_constraint_stdev[name];
		//important: using the initial simulated constraint values
		old_constraint_val = current_constraints_sim_ptr->get_rec(name);
		//old_constraint_val = constraints_sim_initial[name];
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
			//}
		}
		else
			new_constraint_val = current_constraints_sim_ptr->get_rec(name);
		iter_fosm.insert(name, new_constraint_val);
	}

	vector<string> names = iter_fosm.get_keys();
	Observations constraints_chance(*current_constraints_sim_ptr);
	constraints_chance.update_without_clear(names, iter_fosm.get_data_vec(names));
	
	return constraints_chance;
}


vector<double> Constraints::get_constraint_residual_vec(Observations& sim)
{
	vector<double> residuals_vec;
	residuals_vec.resize(num_constraints(), 0.0);

	Observations::const_iterator found_obs;
	Observations::const_iterator not_found_obs = sim.end();
	PriorInformation::const_iterator found_prior_info;

	int i = 0;
	for (vector<string>::iterator b = ctl_ord_obs_constraint_names.begin(), e = ctl_ord_obs_constraint_names.end(); b != e; ++b, ++i)
	{
		found_obs = sim.find(*b);
		//double fosm_offset = post_constraint_offset[*b];
		if (found_obs != not_found_obs)
		{
			residuals_vec[i] = constraints_obs.get_rec(*b) - ((*found_obs).second);
		}
	}
	return residuals_vec;
}

//void Constraints::update_obs_and_pi_constraints(Observations& _constraints_sim, Parameters& _pars_and_dec_vars) 
//{
//	current_pars_and_dec_vars = _pars_and_dec_vars;
//	current_constraints_sim = _constraints_sim;
//}

pair<vector<double>,vector<double>> Constraints::get_constraint_bound_vectors()
{
	vector<double> residuals;
	if (use_chance)
	{
		/*bool update = false;
		if ((!std_weights) && ((iter == 1) || ((iter + 1) % pest_scenario.get_pestpp_options().get_opt_recalc_fosm_every() == 0)))
			update = true;*/
		Observations current_constraints_chance = get_chance_shifted_constraints();
		residuals = get_constraint_residual_vec(current_constraints_chance);
	}
	else
	{
		//current_constraints_sim = constraints_sim;
		residuals = get_constraint_residual_vec(*current_constraints_sim_ptr);
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
		double residual = -constraints_pi.get_pi_rec_ptr(name).calc_residual(*current_pars_and_dec_vars_ptr);
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

void Constraints::presolve_report(int iter)
{
	
	ofstream& f_rec = file_mgr_ptr->rec_ofstream();
	vector<double> residuals = get_constraint_residual_vec(*current_constraints_sim_ptr);
	f_rec << endl << "  observation constraint information at start of iteration " << iter << endl;
	f_rec << setw(20) << left << "name" << right << setw(15) << "sense" << setw(15) << "required" << setw(15) << "sim value";
	f_rec << setw(15) << "residual" << setw(15) << "lower bound" << setw(15) << "upper bound" << endl;
	Observations current;
	if (use_chance)
		current = get_chance_shifted_constraints();
	else
		current = *current_constraints_sim_ptr;
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
			f_rec << setw(15) << pi_rec.calc_residual_and_sim_val(*current_pars_and_dec_vars_ptr).first;
			f_rec << setw(15) << pi_rec.calc_residual(*current_pars_and_dec_vars_ptr);
			f_rec << setw(15) << current_bounds.first[num_obs_constraints() + i];
			f_rec << setw(15) << current_bounds.second[num_obs_constraints() + i] << endl;

		}
	}

	if (use_chance)
		presolve_chance_report(iter);	 
	
	return;
}

void Constraints::write_res_files(Observations& constraints, Parameters& pars_and_dec_vars, string tag, int iter)
{
	write_res_file(constraints, pars_and_dec_vars, tag, iter, false);
	if (use_chance)
		write_res_file(constraints, pars_and_dec_vars, tag, iter, true);
}

void Constraints::write_res_file(Observations& constraints, Parameters& pars_and_dec_vars, string tag, int iter, bool include_chance)
{
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

void Constraints::presolve_chance_report(int iter)
{
	if (!use_chance)
		return;
	
	ofstream& f_rec = file_mgr_ptr->rec_ofstream();
	//vector<double> residuals = get_constraint_residual_vec();
	f_rec << endl << "  FOSM-based chance constraint information at start of iteration " << iter << endl;
	f_rec << setw(20) << left << "name" << right << setw(10) << "sense" << setw(15) << "required" << setw(15) << "sim value";
	f_rec << setw(15) << "prior stdev" << setw(15) << "post stdev" << setw(15) << "offset";
	f_rec << setw(15) << "new sim value" << endl;
	vector<string> out_of_bounds;
	Observations current_constraints_chance = get_chance_shifted_constraints();
	for (int i = 0; i < num_obs_constraints(); ++i)
	{
		string name = ctl_ord_obs_constraint_names[i];
		f_rec << setw(20) << left << name;
		f_rec << setw(10) << right << constraint_sense_name[name];
		f_rec << setw(15) << constraints_obs[name];
		f_rec << setw(15) << current_constraints_sim_ptr->get_rec(name);
		f_rec << setw(15) << prior_constraint_stdev[name];
		f_rec << setw(15) << post_constraint_stdev[name];
		f_rec << setw(15) << post_constraint_offset[name];
		f_rec << setw(15) << current_constraints_chance[name] << endl;
	}
	f_rec << "  note: 'offset' is the value added to the simulated constraint value to account" << endl;
	f_rec << "        for the uncertainty in the constraint value arsing from uncertainty in the " << endl;
	f_rec << "        adjustable parameters identified in the control file." << endl << endl;

	return;
	
}


void Constraints::postsolve_obs_constraints_report(Observations& constraints_sim, string tag, int iter, 
									map<string,string> status_map, map<string,double> price_map)
{
	
	ofstream& f_rec = file_mgr_ptr->rec_ofstream();
	f_rec << endl << endl << "     " << tag << " constraint information at end of iteration " << iter << endl << endl;
	f_rec << setw(20) << left << "name" << right << setw(15) << "sense" << setw(15) << "required" << setw(25);
	if (status_map.size() > 0)
		f_rec << "simplex status";
	if (price_map.size() > 0)
		f_rec << setw(15) << "price";
	if (use_chance)
		f_rec << setw(15) << "fosm offset";
	f_rec << setw(15) << "current" << setw(15) << "residual";
	f_rec << setw(15) << "new" << setw(15) << "residual" << endl;
	vector<double> cur_residuals = get_constraint_residual_vec(*current_constraints_sim_ptr);
	vector<double> new_residuals = get_constraint_residual_vec(constraints_sim);
	double sim_val;
	for (int i = 0; i < num_obs_constraints(); ++i)
	{
		string name = ctl_ord_obs_constraint_names[i];
		sim_val = constraints_sim[name];
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
}

void Constraints::postsolve_pi_constraints_report(Parameters& pars_and_dec_vars, int iter, map<string,string> status_map, map<string,double> price_map)
{
	
	//report prior information constraints
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
			pair<double,double> cur_sim_resid = pi_rec.calc_residual_and_sim_val(*current_pars_and_dec_vars_ptr);
			pair<double,double> new_sim_resid = pi_rec.calc_residual_and_sim_val(pars_and_dec_vars);
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

void Constraints::process_runs(RunManagerAbstract* run_mgr_ptr,int iter)
{
	if (!use_chance)
		return;
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
	}
	else
	{
		//throw runtime_error("non-FOSM update_from_runs() not implemented");
		//string filename = file_mgr_ptr->get_base_filename() + "..stack";
		stringstream ss;
		ss << file_mgr_ptr->get_base_filename() << "." << iter << ".stack";
		if (pest_scenario.get_pestpp_options().get_ies_save_binary())
		{
			ss << ".jcb";
			stack_pe.to_binary(ss.str());
		}
		else
		{
			ss << ".csv";
			stack_pe.to_csv(ss.str());
		}


	}

}

void Constraints::add_runs(RunManagerAbstract* run_mgr_ptr)
{
	if (!use_chance)
		return;
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
		cout << "  ---  running " << num_adj_pars() << " model runs for FOSM-based chance constraints  ---  " << endl;
		bool success = jco.build_runs(*current_pars_and_dec_vars_ptr, *current_constraints_sim_ptr, adj_par_names, pts,
			pest_scenario.get_base_group_info(), pest_scenario.get_ctl_parameter_info(),
			*run_mgr_ptr, out_of_bounds, false, true, false);
		if (!success)
		{
			const set<string> failed = jco.get_failed_parameter_names();
			throw_constraints_error("failed to calc derviatives for the following FOSM parameters: ", failed);
		}
	}
	else
	{
		//throw_constraints_error("non-fosm add_runs() not implemented");
		//todo: update stack_pe parameter values for decision variables using current_pars_and_dec_vars_ptr
		pfm.log_event("building stack-based parameter runs");
		cout << "  ---  running " << stack_pe.shape().first << " model runs for stack-based chance constraints  ---  " << endl;
		stack_pe_run_map.clear();
		stack_pe_run_map = stack_pe.add_runs(run_mgr_ptr);
	}
}

vector<string> Constraints::get_fosm_par_names()
{
	if (use_fosm)
		return adj_par_names;
	else
		return vector<string>();
}

map<string, double> Constraints::get_unsatified_pi_constraints(Parameters& par_and_dec_vars, double tol)
{
	double sim_val, obs_val, scaled_diff;
	map<string, double> unsatisfied;
	
	for (auto& name : ctl_ord_pi_constraint_names)
	{
		PriorInformationRec pi_rec = constraints_pi.get_pi_rec_ptr(name);
		//pair<double, double> cur_sim_resid = pi_rec.calc_residual_and_sim_val(*current_pars_and_dec_vars_ptr);
		pair<double, double> new_sim_resid = pi_rec.calc_residual_and_sim_val(par_and_dec_vars);
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

map<string, double> Constraints::get_unsatified_obs_constraints(Observations& constraints_sim, double tol)
{
	double sim_val, obs_val, scaled_diff;
	map<string, double> unsatisfied;
	for (int i = 0; i < num_obs_constraints(); ++i)
	{
		string name = ctl_ord_obs_constraint_names[i];
		sim_val = constraints_sim[name];
		if (use_chance)
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
	int num = 0;

	for (auto &pi_name : ctl_ord_pi_constraint_names)
	{
		num += constraints_pi.get_pi_rec_ptr(pi_name).get_atom_factors().size();
	}
	return num;
}



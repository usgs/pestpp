/*


	This file is part of PEST++.

	PEST++ is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	PEST++ is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with PEST++.  If not, see<http://www.gnu.org/licenses/>.
*/

#include <cstdlib>
#include <vector>
#include <set>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "Jacobian_1to1.h"
#include "Transformable.h"
#include "ParamTransformSeq.h"
#include "pest_error.h"
#include "pest_data_structs.h"
#include "ModelRunPP.h"
#include "RunManagerAbstract.h"
#include "ObjectiveFunc.h"
#include "utilities.h"
#include "FileManager.h"
#include "PriorInformation.h"
#include "debug.h"
#include "OutputFileWriter.h"

using namespace std;
using namespace pest_utils;

Jacobian_1to1::Jacobian_1to1(FileManager &_file_manager, OutputFileWriter &_output_file_writer) : Jacobian(_file_manager)
{
	output_file_writer_ptr = &_output_file_writer;
}

Jacobian_1to1::~Jacobian_1to1() {
}

bool Jacobian_1to1::build_runs(Parameters &ctl_pars, Observations &ctl_obs, vector<string> numeric_par_names, ParamTransformSeq &par_transform,
	const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info,
	RunManagerAbstract &run_manager, set<string> &out_of_bound_par, bool phiredswh_flag, bool calc_init_obs, bool reinitiailize)
{
	Parameters model_parameters(par_transform.ctl2model_cp(ctl_pars));
	base_numeric_parameters = par_transform.ctl2numeric_cp(ctl_pars);
	if (reinitiailize)
		run_manager.reinitialize(file_manager.build_filename("rnj"));
	debug_msg("Jacobian_1to1::build_runs begin");

	failed_parameter_names.clear();
	failed_ctl_parameters.clear();
	par_run_map.clear();
	// add base run
	int run_id = run_manager.add_run(model_parameters, "par_name:__base__", 0);
	par_run_map["__base__"] = vector<int>{ run_id };
	//if base run is has already been complete, update it and mark it as complete
	// compute runs for to jacobain calculation as it is influenced by derivative type( forward or central)
	if (!calc_init_obs) {
		const Observations &init_obs = ctl_obs;
		run_manager.update_run(run_id, model_parameters, init_obs);
	}

	Parameters new_derivative_pars;
	bool success;
	Parameters base_derivative_parameters = par_transform.numeric2active_ctl_cp(base_numeric_parameters);
	Parameters base_model_parameters = par_transform.numeric2model_cp(base_numeric_parameters);
	//Loop through derivative parameters and build the parameter sets necessary for computing the jacobian
	
	for (auto &i_name : numeric_par_names)
	{
		assert(base_derivative_parameters.find(i_name) != base_derivative_parameters.end());
		vector<double> tmp_del_numeric_par_vec;
		double derivative_par_value = base_derivative_parameters.get_rec(i_name);
		success = get_derivative_parameters(i_name, derivative_par_value, par_transform, group_info, ctl_par_info,
			tmp_del_numeric_par_vec, phiredswh_flag);
		if (success && !tmp_del_numeric_par_vec.empty())
		{
			// update changed model parameters in model_parameters
			for (const auto &par : tmp_del_numeric_par_vec)
			{
				Parameters new_pars;
				new_pars.insert(make_pair(i_name, par));
				par_transform.active_ctl2model_ip(new_pars);
				for (auto &ipar : new_pars)
				{
					model_parameters[ipar.first] = ipar.second;
				}
				int id = run_manager.add_run(model_parameters, "par_name:"+i_name, par);
				if (par_run_map.count(i_name) == 0)
					par_run_map[i_name] = vector<int>{ id };
				else
					par_run_map[i_name].push_back(id);
				//reset the perturbed parameters back to the values associated with the base condition
				for (const auto &ipar : new_pars)
				{
					model_parameters[ipar.first] = base_model_parameters[ipar.first];
				}
			}
		}
		else
		{
			cout << endl << " warning: failed to compute parameter derivative for " << i_name << endl;
			file_manager.rec_ofstream() << " warning: failed to compute parameter derivative for " << i_name << endl;
			failed_parameter_names.insert(i_name);
			failed_to_increment_parmaeters.insert(i_name, derivative_par_value);
		}
	}
	output_file_writer_ptr->write_jco_run_id(run_manager.get_cur_groupid(), par_run_map);
	ofstream &fout_restart = file_manager.get_ofstream("rst");
	debug_print(failed_parameter_names);
	debug_msg("Jacobian_1to1::build_runs end");
	if (failed_parameter_names.size() > 0)
		return false;
	return true;
}


bool Jacobian_1to1::build_runs(ModelRun &init_model_run, vector<string> numeric_par_names, ParamTransformSeq &par_transform,
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info,
		RunManagerAbstract &run_manager, set<string> &out_of_bound_par, bool phiredswh_flag, bool calc_init_obs)
{
	//Parameters pars = init_model_run.get_ctl_pars();
	//Observations obs = init_model_run.get_obs_template();
	//return build_runs(pars, obs, numeric_par_names, par_transform, group_info, ctl_par_info, run_manager, out_of_bound_par, phiredswh_flag, calc_init_obs);
	Parameters model_parameters(par_transform.ctl2model_cp(init_model_run.get_ctl_pars()));
	base_numeric_parameters = par_transform.ctl2numeric_cp(init_model_run.get_ctl_pars());
	run_manager.reinitialize(file_manager.build_filename("rnj"));
	debug_msg("Jacobian_1to1::build_runs begin");

	failed_parameter_names.clear();
	failed_ctl_parameters.clear();
	par_run_map.clear();
	// add base run
	int run_id = run_manager.add_run(model_parameters, "par_name:__base__", 0);
	par_run_map["__base__"] = vector<int>{ run_id };
	//if base run is has already been complete, update it and mark it as complete
	// compute runs for to jacobain calculation as it is influenced by derivative type( forward or central)
	if (!calc_init_obs) {
		const Observations &init_obs = init_model_run.get_obs();
		run_manager.update_run(run_id, model_parameters, init_obs);
	}

	Parameters new_derivative_pars;
	bool success;
	Parameters base_derivative_parameters = par_transform.numeric2active_ctl_cp(base_numeric_parameters);
	Parameters base_model_parameters = par_transform.numeric2model_cp(base_numeric_parameters);
	//Loop through derivative parameters and build the parameter sets necessary for computing the jacobian
	
	for (auto &i_name : numeric_par_names)
	{
		assert(base_derivative_parameters.find(i_name) != base_derivative_parameters.end());
		vector<double> tmp_del_numeric_par_vec;
		double derivative_par_value = base_derivative_parameters.get_rec(i_name);
		success = get_derivative_parameters(i_name, derivative_par_value, par_transform, group_info, ctl_par_info,
			tmp_del_numeric_par_vec, phiredswh_flag);
		if (success && !tmp_del_numeric_par_vec.empty())
		{
			// update changed model parameters in model_parameters
			for (const auto &par : tmp_del_numeric_par_vec)
			{
				Parameters new_pars;
				new_pars.insert(make_pair(i_name, par));
				par_transform.active_ctl2model_ip(new_pars);
				for (auto &ipar : new_pars)
				{
					model_parameters[ipar.first] = ipar.second;
				}
				int id = run_manager.add_run(model_parameters, "par_name:"+i_name, par);

				if (par_run_map.count(i_name) == 0)
					par_run_map[i_name] = vector<int>{ id };
				else
					par_run_map[i_name].push_back(id);
				//reset the perturbed parameters back to the values associated with the base condition
				for (const auto &ipar : new_pars)
				{
					model_parameters[ipar.first] = base_model_parameters[ipar.first];
				}
			}
		}
		else
		{
			cout << endl << " warning: failed to compute parameter derivative for " << i_name << endl;
			file_manager.rec_ofstream() << " warning: failed to compute parameter derivative for " << i_name << endl;
			failed_parameter_names.insert(i_name);
			failed_to_increment_parmaeters.insert(i_name, derivative_par_value);
		}
	}
	output_file_writer_ptr->write_jco_run_id(run_manager.get_cur_groupid(), par_run_map);

	ofstream &fout_restart = file_manager.get_ofstream("rst");
	debug_print(failed_parameter_names);
	debug_msg("Jacobian_1to1::build_runs end");

	output_file_writer_ptr->write_jco_run_id(run_manager.get_cur_groupid(), par_run_map);
	if (failed_parameter_names.size() == numeric_par_names.size())
	{
		throw runtime_error("Jacobian_1to1::build_runs() error: parameter derivative calculations failed for all parameters");
	}
	if (failed_parameter_names.size() > 0)
		return false;
	return true;
}


void Jacobian_1to1::make_runs(RunManagerAbstract &run_manager)
{
	// make model runs
	run_manager.run();
}

bool Jacobian_1to1::process_runs(ParamTransformSeq &par_transform,
		const ParameterGroupInfo &group_info,
		RunManagerAbstract &run_manager, 
	const PriorInformation &prior_info, bool splitswh_flag,
	bool debug_fail)
{
	debug_msg("Jacobian_1to1::process_runs begin");
       base_sim_obs_names = run_manager.get_obs_name_vec();
	vector<string> prior_info_name = prior_info.get_keys();
	base_sim_obs_names.insert(base_sim_obs_names.end(), prior_info_name.begin(), prior_info_name.end());
	std::vector<Eigen::Triplet<double> > triplet_list;

	unordered_map<string, int> par2col_map = get_par2col_map();
	unordered_map<string, int>::iterator found;
	unordered_map<string, int>::iterator not_found = par2col_map.end();

	JacobianRun base_run;
	int i_run = 0;
	// get base run parameters and observation for initial model run from run manager storage
	/*run_manager.get_model_parameters(i_run,  base_run.ctl_pars);
	bool success = run_manager.get_observations_vec(i_run, base_run.obs_vec);
	if (!success)
	{
		throw(PestError("Error: Base parameter run failed.  Can not compute the Jacobian"));
	}
	par_transform.model2ctl_ip(base_run.ctl_pars);
	base_numeric_parameters = par_transform.ctl2numeric_cp(base_run.ctl_pars);
	++i_run;
	*/
	string base = "__base__";
	if (par_run_map.find(base) != par_run_map.end())
	{
		run_manager.get_model_parameters(par_run_map[base][0], base_run.ctl_pars);
		bool success = run_manager.get_observations_vec(par_run_map[base][0], base_run.obs_vec);
		if (!success)
		{
			throw(PestError("Error: Base parameter run failed.  Can not compute the Jacobian"));
		}
		par_transform.model2ctl_ip(base_run.ctl_pars);
		base_numeric_parameters = par_transform.ctl2numeric_cp(base_run.ctl_pars);
		par_run_map.erase(base);
		i_run++;
	}
	else
	{
		throw runtime_error("'__base__' run tag not in Jacobian.par_run_map, something is wrong");
	}
		

		

	// process the parameter perturbation runs
	int nruns = run_manager.get_nruns();
	base_numeric_par_names.clear();
	int icol = 0;
	int r_status;
	vector<string>par_name_vec;
	string cur_par_name;
	string par_name_next;
	int run_status_next;
	double par_value_next;
	double cur_numeric_par_value;
	list<JacobianRun> run_list;
    int nfailed = 0;
	
	for (auto par_run : par_run_map)
    {
		for (auto rid : par_run.second)
		{
			run_list.push_back(JacobianRun());
			run_manager.get_info(par_run.second[0], r_status, cur_par_name, cur_numeric_par_value);
			cur_par_name = cur_par_name.substr(9,cur_par_name.size());
			run_manager.get_model_parameters(par_run.second[0], run_list.back().ctl_pars);
			bool success = run_manager.get_observations_vec(par_run.second[0], run_list.back().obs_vec);
			run_list.back().numeric_derivative_par = cur_numeric_par_value;
			if ((debug_fail) && (nfailed < 3))
			{
				file_manager.rec_ofstream() << "NOTE: 'GLM_DEBUG_DER_FAIL' is true, failing jco run for parameter '" << cur_par_name << "'" << endl;
				success = false;
				nfailed++;
			}
			if (success)
			{
				par_transform.model2ctl_ip(run_list.back().ctl_pars);
				// get the updated parameter value which reflects roundoff errors

				par_name_vec.clear();
				par_name_vec.push_back(cur_par_name);
				Parameters numeric_pars(run_list.back().ctl_pars, par_name_vec);
				par_transform.ctl2numeric_ip(numeric_pars);
				run_list.back().numeric_derivative_par = numeric_pars.get_rec(cur_par_name);
			}
			else
			{
				run_list.pop_back();
			}
		}

		
		if (!run_list.empty())
		{
			base_numeric_par_names.push_back(cur_par_name);
			base_run.numeric_derivative_par = base_numeric_parameters.get_rec(cur_par_name);
			double cur_numeric_value = base_run.numeric_derivative_par;
			run_list.push_front(base_run);
			std::vector<Eigen::Triplet<double> > tmp_triplet_vec = calc_derivative(cur_par_name, cur_numeric_value, icol, run_list, group_info, prior_info, splitswh_flag);
			triplet_list.insert(triplet_list.end(), tmp_triplet_vec.begin(), tmp_triplet_vec.end());
			icol++;
		}
		else
		{
			failed_parameter_names.insert(cur_par_name);
			failed_ctl_parameters.insert(cur_par_name, cur_numeric_par_value);
		}
		run_list.clear();
		

	}
	par_run_map.clear();
	matrix.resize(base_sim_obs_names.size(), base_numeric_par_names.size());
	matrix.setZero();
	matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
	// clean up
	ofstream &fout_restart = file_manager.get_ofstream("rst");
	run_manager.free_memory();
	debug_print(failed_parameter_names);
	debug_msg("Jacobian_1to1::process_runs end");
	return true;
}

bool Jacobian_1to1::get_derivative_parameters(const string &par_name, double par_value, const ParamTransformSeq &par_trans, const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info,
		vector<double> &delta_numeric_par_vec, bool phiredswh_flag)
{
	bool success = false;
	const ParameterGroupRec g_rec = group_info.get_group_rec(par_name);

	if (g_rec.forcen != "ALWAYS_2"  && (g_rec.forcen == "ALWAYS_3" || phiredswh_flag == true) ) {
		// Central Difference
		vector<double> new_par_vec;
		vector<Parameters> dir_numeric_pars_vec;
		success = central_diff(par_name, par_value, group_info, ctl_par_info, par_trans, new_par_vec, dir_numeric_pars_vec);
		if (success)
		{
			for (auto & ipar : new_par_vec)
			{
				delta_numeric_par_vec.push_back(ipar);
			}
		}
	}
	if (!success) {
		// Forward Difference
		success = forward_diff(par_name, par_value, group_info, ctl_par_info, par_trans, par_value);
		if(success)
		{
			delta_numeric_par_vec.push_back(par_value);
		}
	}
	return success;
}

bool Jacobian_1to1::forward_diff(const string &par_name, double base_derivative_val,
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, const ParamTransformSeq &par_trans, double &new_par_val)
{
	const ParameterRec *par_info_ptr = ctl_par_info.get_parameter_rec_ptr(par_name);
	Parameters new_par;
	bool out_of_bound_forward;
	bool out_of_bound_backward;
	vector<string> out_of_bound__forward_par_vec;
	vector<string> out_of_bound__backard_par_vec;
	string tmp_name;

	// perturb derivative parameters
	double incr = derivative_inc(par_name, group_info, base_derivative_val, false);
	if (incr == 0.0)
        return false;
	new_par_val = new_par[par_name] = base_derivative_val + incr;
	// try forward derivative
	out_of_bound_forward = out_of_bounds(new_par, par_info_ptr);
	if (!out_of_bound_forward) {
		return true;
	}
	// try backward derivative if forward derivative didn't work
	new_par.clear();
	new_par_val = new_par[par_name] = base_derivative_val - incr;
	out_of_bound_backward = out_of_bounds(new_par, par_info_ptr);
	if (!out_of_bound_backward)
	{
		return true;
	}
	return false;
}

bool Jacobian_1to1::central_diff(const string &par_name, double base_derivative_val,
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, const ParamTransformSeq &par_trans, vector<double> &new_par_vec,
		vector<Parameters>  &perturb_derivative_par_vec)
{
	double new_par;
	bool out_of_bnds_forward, out_of_bnds_back, out_of_bnds;
	Parameters perturb_derivative_pars;
	string tmp_name;

	const ParameterRec *par_info_ptr = ctl_par_info.get_parameter_rec_ptr(par_name);
	double incr = derivative_inc(par_name, group_info, base_derivative_val, true);
	if (incr == 0.0) return false;
	// try backward difference
	new_par = perturb_derivative_pars[par_name] = base_derivative_val - incr;
	out_of_bnds_back = out_of_bounds(perturb_derivative_pars, par_info_ptr);

	if (!out_of_bnds_back) {
		new_par_vec.push_back(new_par);
		perturb_derivative_par_vec.push_back(perturb_derivative_pars);
	}
	// try forward derivative
	new_par = perturb_derivative_pars[par_name] = base_derivative_val + incr;
	out_of_bnds_forward = out_of_bounds(perturb_derivative_pars, par_info_ptr);
	if (!out_of_bnds_forward) {
		new_par_vec.push_back(new_par);
		perturb_derivative_par_vec.push_back(perturb_derivative_pars);
	}
	// if backward difference was out of bounds do a second forward derivative
	if (out_of_bnds_back) {
		new_par = perturb_derivative_pars[par_name] = base_derivative_val + 2.0 * incr;
		out_of_bnds = out_of_bounds(perturb_derivative_pars, par_info_ptr);
		if (!out_of_bnds) {
			new_par_vec.push_back(new_par);
			perturb_derivative_par_vec.push_back(perturb_derivative_pars);
		}
		else
		{
			return false;  // can't do central difference without going out of bounds
		}
	}
	// if forward difference was out of bounds do a second backward derivative
	if (out_of_bnds_forward) {
		new_par = perturb_derivative_pars[par_name] = base_derivative_val - 2.0 * incr;
		out_of_bnds = out_of_bounds(perturb_derivative_pars, par_info_ptr);
		if (!out_of_bnds) {
			new_par_vec.insert(new_par_vec.begin(), new_par);
			perturb_derivative_par_vec.push_back(perturb_derivative_pars);
		}
		else
		{
			return false;  // can't do central difference without going out of bounds
		}
	}
	return true;
}

bool Jacobian_1to1::out_of_bounds(const Parameters &ctl_parameters,
	const ParameterRec *par_info_ptr) const
{
	bool out_of_bounds=false;

        // This will always only contain one entry one 1 to 1 Jacobians
	for (auto &p : ctl_parameters)
	{
		double max = par_info_ptr->ubnd;
		double min = par_info_ptr->lbnd;
		if (p.second > max || p.second < min) {
			out_of_bounds = true;
		}
	}
	return out_of_bounds;
}

void Jacobian_1to1::report_errors(std::ostream &fout)
{
	if (failed_to_increment_parmaeters.size() > 0)
	{
		fout << "    Parameters that went out of bounds while computing jacobian" << endl;
		fout << "      Parameter" << endl;
		fout << "        Name" << endl;
		fout << "      ----------" << endl;
	}
	for (const auto & ipar : failed_to_increment_parmaeters)
	{
		fout << right;
		fout << "  " << setw(12) << ipar.first << endl;
	}

	if (failed_to_increment_parmaeters.size() > 0)
	{
		fout << endl;
	}

	if (failed_ctl_parameters.size() > 0)
	{
		fout << "    Parameters whose perturbation runs failed while computing jacobian" << endl;
		fout << "      Parameter     Failed" << endl;
		fout << "        Name        Value" << endl;
		fout << "      ----------  ------------" << endl;
	}

	for (const auto & ipar : failed_ctl_parameters)
	{
		fout << right;
		fout << "  " << setw(12) << ipar.first;
		fout << right;
		fout << "  " << setw(12) << ipar.second << endl;
	}


}

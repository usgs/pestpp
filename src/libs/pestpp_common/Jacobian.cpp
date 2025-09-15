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
#include <fstream>
#include <iomanip>
#include "Jacobian.h"
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
#include "eigen_tools.h"
#include "Pest.h"

using namespace std;
using namespace pest_utils;
using namespace Eigen;

Jacobian::Jacobian(FileManager &_file_manager) : file_manager(_file_manager)
{
}


void Jacobian::set_base_numeric_pars(Parameters _base_numeric_pars)
{
	base_numeric_parameters = _base_numeric_pars;
}

void Jacobian::set_base_sim_obs(Observations _base_sim_obs)
{
	base_sim_observations = _base_sim_obs;
}

Jacobian::~Jacobian() {
}


void Jacobian::remove_cols(std::set<string> &rm_parameter_names)
{
	vector<size_t> del_col_ids;

	// build list of columns that needs to be removed from the matrix
	auto iter_rm_end = rm_parameter_names.end();
	size_t npar = base_numeric_par_names.size();
	for (int i = 0; i<npar; ++i)
	{
		if (rm_parameter_names.find(base_numeric_par_names[i]) != iter_rm_end)
		{
			del_col_ids.push_back(i);
		}
	}

	//remove frozen parameters from base_parameter_names
	auto end_iter = std::remove_if(base_numeric_par_names.begin(), base_numeric_par_names.end(),
		[&rm_parameter_names](string &str)->bool{return rm_parameter_names.find(str) != rm_parameter_names.end(); });
	base_numeric_par_names.resize(std::distance(base_numeric_par_names.begin(), end_iter));
	matrix_del_rows_cols(matrix, del_col_ids,false,true);
}

void Jacobian::remove_rows(std::set<string>& rm_obs_names)
{
	vector<size_t> del_row_ids;

	// build list of columns that needs to be removed from the matrix
	auto iter_rm_end = rm_obs_names.end();
	vector<string> row_names = observation_list();
	size_t nobs = row_names.size();

	for (int i = 0; i < nobs; ++i)
	{
		if (rm_obs_names.find(row_names[i]) != iter_rm_end)
		{
			del_row_ids.push_back(i);
		}
	}

	auto end_iter = std::remove_if(base_sim_obs_names.begin(), base_sim_obs_names.end(),
		[&rm_obs_names](string& str)->bool {return rm_obs_names.find(str) != rm_obs_names.end(); });
	base_sim_obs_names.resize(std::distance(base_sim_obs_names.begin(), end_iter));
	matrix_del_rows_cols(matrix, del_row_ids, true, false);
}

void Jacobian::add_cols(set<string> &new_pars_names)
{
	//Note:  This method does not add the parameters in the base_numeric_parameters container.
	//       The values must already be in that container
	// check if any of new_par parameters are already in the jacobian
	set<string> par_name_set;
	set<string> repeated_pars;
	par_name_set.insert(base_numeric_par_names.begin(), base_numeric_par_names.end());
	for_each(new_pars_names.begin(), new_pars_names.end(), [&par_name_set, &repeated_pars](const string &p) {
		if ((par_name_set.find(p)) != par_name_set.end())
		{
			repeated_pars.insert(p);
		}
	});
	if (repeated_pars.size() > 0)
	{
		ostringstream str;
		str << " Jacobian::add_cols - parameters already present in jacobian: ";
		for (auto &ipar : repeated_pars)
		{
			str << " " << ipar;
		}
		throw PestError(str.str());
	}


	for (const auto &ipar : new_pars_names)
	{
		base_numeric_par_names.push_back(ipar);
	}
// add empty columns for new parameter.  sensitivities for new parameters will all = 0.
matrix.resize(matrix.rows(), matrix.cols() + new_pars_names.size());
}



const vector<string>& Jacobian::obs_and_reg_list() const
{
	return base_sim_obs_names;
}

unordered_map<string, int> Jacobian::get_par2col_map() const
{
	unordered_map<string, int> par2col_map;
	int icol_old = 0;
	for (vector<string>::const_iterator b = base_numeric_par_names.begin(), e = base_numeric_par_names.end();
		b != e; ++b, ++icol_old)
	{
		par2col_map[(*b)] = icol_old;
	}
	return par2col_map;
}

unordered_map<string, int> Jacobian::get_obs2row_map() const
{
	unordered_map<string, int> obs2row_map;
	int irow_old = 0;
	for (vector<string>::const_iterator b = base_sim_obs_names.begin(), e = base_sim_obs_names.end();
		b != e; ++b, ++irow_old)
	{
		obs2row_map[(*b)] = irow_old;
	}
	return obs2row_map;
}


Eigen::SparseMatrix<double> Jacobian::get_matrix(const vector<string> &obs_names, const vector<string> & par_names, bool forgive_missing, int n_cols) const
{
    /* the n_cols arg is so you can reserve a sparse matrix with columns that are all zeros - this is for the LP solver in pestpp-opt when you have
     * external dec vars - they won't be in the Jacobian instance but they need to be in the LP solution matrix.
     */
    stringstream ss;
	int n_rows = obs_names.size();
	if (n_cols == 0)
	    n_cols = par_names.size();
	int irow_new;
	int icol_new;

	if (!forgive_missing) {

        set<string> sobs_names(base_sim_obs_names.begin(), base_sim_obs_names.end());
        set<string> spar_names(base_numeric_par_names.begin(), base_numeric_par_names.end());

        set<string>::iterator send = sobs_names.end();
        vector<string> missing;
        for (auto &obs_name : obs_names)
            if (sobs_names.find(obs_name) == send)
                missing.push_back(obs_name);
        if (missing.size() > 0) {
            ss.str("");

            ss << "Jco::get_matrix(): the following obs names are not in the matrix:";
            for (auto &m : missing)
                ss << " " << m;
            throw runtime_error(ss.str());
        }

        send = spar_names.end();
        for (auto &par_name : par_names)
            if (spar_names.find(par_name) == send)
                missing.push_back(par_name);
        if (missing.size() > 0){
            ss.str("");

            ss << "Jco::get_matrix(): the following par names are not in the matrix:";
            for (auto &m : missing)
                ss << " " << m;
            throw runtime_error(ss.str());
        }
    }

	unordered_map<string, int> obs_name2newindex_map;
	unordered_map<string, int> par_name2new_index_map;

	// Build mapping of parameter names to column number in new matrix to be returned
	icol_new = 0;
	for (vector<string>::const_iterator b = par_names.begin(), e = par_names.end();
		b != e; ++b, ++icol_new) {
		par_name2new_index_map[(*b)] = icol_new;
	}

	// Build mapping of observation names to row  number in new matrix to be returned
	irow_new = 0;
	for (vector<string>::const_iterator b = obs_names.begin(), e = obs_names.end();
		b != e; ++b, ++irow_new) {
		obs_name2newindex_map[(*b)] = irow_new;
	}

	//build jacobian
	unordered_map<string, int>::const_iterator found_par;
	unordered_map<string, int>::const_iterator found_obs;
	unordered_map<string, int>::const_iterator not_found_par_map = par_name2new_index_map.end();
	unordered_map<string, int>::const_iterator not_found_obs_map = obs_name2newindex_map.end();

	//unordered_map<string, int>::const_iterator not_found_par2col_map = par2col_map.end();
	//map<string, map<string, double>>::const_iterator found_prior_info;
	//map<string, map<string, double>>::const_iterator not_found_prior_info = prior_info_sen.end();

	double data;
	const string *obs_name;
	const string *par_name;
	std::vector<Eigen::Triplet<double> > triplet_list;
	for (int icol = 0; icol < matrix.outerSize(); ++icol)
	{
		for (SparseMatrix<double>::InnerIterator it(matrix, icol); it; ++it)
		{
			data = it.value();
			par_name = &base_numeric_par_names[it.col()];
			obs_name = &base_sim_obs_names[it.row()];
			found_par = par_name2new_index_map.find(*par_name);
			found_obs = obs_name2newindex_map.find(*obs_name);

			if (found_par != not_found_par_map && found_obs != not_found_obs_map)
			{
				triplet_list.push_back(Eigen::Triplet<double>(found_obs->second, found_par->second, data));
			}
		}
	}
	Eigen::SparseMatrix<double> new_matrix(n_rows, n_cols);
	new_matrix.setZero();  // initialize all entries to 0
	new_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
	return new_matrix;
}


bool Jacobian::build_runs(Parameters &ctl_pars, Observations &ctl_obs, vector<string> numeric_par_names, ParamTransformSeq &par_transform,
	const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info,
	RunManagerAbstract &run_manager, set<string> &out_of_bound_par, bool phiredswh_flag, bool calc_init_obs)
{

	run_manager.reinitialize(file_manager.build_filename("rnj"));
	failed_parameter_names.clear();

	debug_msg("Jacobian::build_runs method: begin");
	// add base run
	Parameters model_pars = par_transform.ctl2model_cp(ctl_pars);
	int run_id = run_manager.add_run(model_pars, "par_name:__base__", 0);

	if (!calc_init_obs) {
		const Observations &observations = ctl_obs;
		run_manager.update_run(run_id, model_pars, observations);
	}

	// compute runs for to jacobain calculation as it is influenced by derivative type( forward or central)
	out_of_bound_par.clear();
	Parameters numeric_pars = par_transform.ctl2numeric_cp(ctl_pars);

	vector<double> del_numeric_par_vec;
	for (const auto &ipar_name : numeric_par_names)
	{
		debug_print(ipar_name);
		del_numeric_par_vec.clear();
		// need to optimize already computing model pars in get_derivative_parameters.  should not need to compute them again
		bool success = get_derivative_parameters(ipar_name, numeric_pars, par_transform, group_info, ctl_par_info, del_numeric_par_vec, phiredswh_flag, out_of_bound_par);
		debug_print(out_of_bound_par);
		if (success && !del_numeric_par_vec.empty())
		{
			debug_msg("success");
			Parameters numeric_parameters = par_transform.ctl2numeric_cp(ctl_pars);
			for (double ipar_val : del_numeric_par_vec)
			{
				numeric_parameters.update_rec(ipar_name, ipar_val);
				Parameters model_parameters = par_transform.numeric2model_cp(numeric_parameters);
				run_manager.add_run(model_parameters, "par_name:"+ipar_name, ipar_val);
			}
		}
		else
		{
			debug_msg("fail");
			//cout << endl << " warning: failed to compute parameter derivative for " << ipar_name << endl;
			file_manager.rec_ofstream() << " warning: failed to compute parameter derivative for " << ipar_name << endl;
			failed_parameter_names.insert(ipar_name);
		}
	}
	debug_print(failed_parameter_names);
	debug_msg("Jacobian::build_runs method: end");
	
	if (failed_parameter_names.size() > 0)
	{
		return false;
	}
	return true;
}



bool Jacobian::build_runs(ModelRun &init_model_run, vector<string> numeric_par_names, ParamTransformSeq &par_transform,
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info,
		RunManagerAbstract &run_manager, set<string> &out_of_bound_par, bool phiredswh_flag, bool calc_init_obs)
{
	Parameters pars = init_model_run.get_ctl_pars();
	Observations obs = init_model_run.get_obs();
	return build_runs(pars, obs, numeric_par_names, par_transform, group_info, ctl_par_info, run_manager,
		out_of_bound_par, phiredswh_flag, calc_init_obs);
}

void Jacobian::make_runs(RunManagerAbstract &run_manager)
{
	// make model runs
	run_manager.run();
}

bool Jacobian::process_runs(ParamTransformSeq &par_transform,
		const ParameterGroupInfo &group_info,
		RunManagerAbstract &run_manager, 
		const PriorInformation &prior_info, bool splitswh_flag,
		bool debug_fail)
{
	// calculate jacobian
    base_sim_obs_names = run_manager.get_obs_name_vec();
	vector<string> prior_info_name = prior_info.get_keys();
	base_sim_obs_names.insert(base_sim_obs_names.end(), prior_info_name.begin(), prior_info_name.end());
	std::vector<Eigen::Triplet<double> > triplet_list;

	JacobianRun base_run;
	int i_run = 0;
	// get base run parameters and observation for initial model run from run manager storage
	{
			run_manager.get_model_parameters(i_run,  base_run.ctl_pars);
			bool success = run_manager.get_observations_vec(i_run, base_run.obs_vec);
		if (!success)
		{
			throw(PestError("Error: Super-parameter base parameter run failed.  Can not compute the Jacobian"));
		}
		par_transform.model2ctl_ip(base_run.ctl_pars);
		base_numeric_parameters = par_transform.ctl2numeric_cp(base_run.ctl_pars);
		++i_run;
	}

	// process the parameter perturbation runs
	int nruns = run_manager.get_nruns();
	int icol = 0;
	int r_status;
	string cur_par_name;
	string par_name_next;
	int run_status_next;
	double par_value_next;
	double cur_numeric_par_value;
	list<JacobianRun> run_list;
	base_numeric_par_names.clear();
	int nfailed = 0;

	for(; i_run<nruns; ++i_run)
	{
		run_list.push_back(JacobianRun());
				run_manager. get_info(i_run, r_status, cur_par_name, cur_numeric_par_value);
            //this is to strip off the "par_name:" tag
            cur_par_name = cur_par_name.substr(9,cur_par_name.size());
            run_manager.get_model_parameters(i_run,  run_list.back().ctl_pars);
			bool success = run_manager.get_observations_vec(i_run, run_list.back().obs_vec);
			if ((debug_fail) && (nfailed < 2))
            {
				file_manager.rec_ofstream() << "NOTE: 'GLM_DEBUG_DER_FAIL' is true, failing jco run for parameter '" << cur_par_name << "'" << endl;
				success = false;
				nfailed++;
			}
		if (success)
		{
			par_transform.model2ctl_ip(run_list.back().ctl_pars);
			run_list.back().numeric_derivative_par = cur_numeric_par_value;
		}
		else
		{
			run_list.pop_back();
		}

		// read information associated with the next model run;
		if (i_run+1<nruns)
		{
			run_manager.get_info(i_run+1, run_status_next, par_name_next, par_value_next);
			//again the par_name: tag
			par_name_next = par_name_next.substr(9,par_name_next.size());
		}

		if (i_run + 1 >= nruns || (cur_par_name != par_name_next))
		{
			if (!run_list.empty())
			{
				base_numeric_par_names.push_back(cur_par_name);
				double base_numeric_par_value = base_numeric_parameters.get_rec(cur_par_name);
				base_run.numeric_derivative_par = base_numeric_par_value;
				run_list.push_front(base_run);
				std::vector<Eigen::Triplet<double> > tmp_triplet_vec = calc_derivative(cur_par_name, base_numeric_par_value, icol, run_list, group_info, prior_info, splitswh_flag);
				triplet_list.insert(triplet_list.end(), tmp_triplet_vec.begin(), tmp_triplet_vec.end());
				icol++;
				run_list.clear();
			}
			else
			{
				failed_parameter_names.insert(cur_par_name);
				throw(PestError("Error: All runs for parameter: " + cur_par_name
					+ " failed.  Cannot compute the Jacobian"));
			}
		}
	}
	matrix.resize(base_sim_obs_names.size(), base_numeric_par_names.size());
	matrix.setZero();
	matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
	// clean up
	run_manager.free_memory();
	return true;
}

bool Jacobian::get_derivative_parameters(const string &par_name, Parameters &numeric_pars, ParamTransformSeq &par_transform, 
	const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info,
		vector<double> &delta_numeric_par_vec, bool phiredswh_flag, set<string> &out_of_bound_par)
{
	//first check that all active pars are in bounds
	bool already_out = out_of_bounds(par_transform.numeric2ctl_cp(numeric_pars), ctl_par_info, out_of_bound_par);
	if (already_out)
	{
		stringstream ss;
		ss << "Jacobian::get_derivative_parameters() error: the following parameters are already out of bounds: " << endl;
		for (auto p : out_of_bound_par)
			ss << p << endl;
		//throw runtime_error(ss.str());

	}
	bool success = false;
	ParameterGroupRec g_rec;
	debug_msg("Jacobian::get_derivative_parameters begin");
	debug_print(par_name);
	g_rec = group_info.get_group_rec(par_name);

	if (g_rec.forcen == "ALWAYS_3" || phiredswh_flag == true) {
		debug_msg("trying central");
		// Central Difference
		vector<double> new_par_vec;
		vector<Parameters> dir_numeric_pars_vec;
		success = central_diff(par_name, numeric_pars, group_info, ctl_par_info, par_transform, new_par_vec, dir_numeric_pars_vec, out_of_bound_par);
		if (success)
		{
			debug_msg("central successful");
			for (auto & ipar : new_par_vec)
			{
				delta_numeric_par_vec.push_back(ipar);
			}
		}
	}
	if (!success) {
	  double par_1;
		// Forward Difference
	  debug_msg("trying forward");
		success = forward_diff(par_name, numeric_pars, group_info, ctl_par_info, par_transform, par_1, out_of_bound_par);
		if (success)
		{
			debug_msg("forward successful");
			delta_numeric_par_vec.push_back(par_1);
		}
	}
	debug_print(success);
	debug_msg("Jacobian::get_derivative_parameters end");
	return success;
}


std::vector<Eigen::Triplet<double> >  Jacobian::calc_derivative(const string &numeric_par_name, double base_numeric_par_value, int jcol, list<JacobianRun> &run_list,
	const ParameterGroupInfo &group_info, const PriorInformation &prior_info, bool splitswh_flag)
{
	ParameterGroupRec g_rec;
	double del_par;
	double del_obs;
	double der;
	int irow;
	Parameters::const_iterator par_iter;
	std::vector<Eigen::Triplet<double> > triplet_list;

	// sort run_list the parameter numeric_par_name;
	auto compare = [](const JacobianRun &a, const JacobianRun &b)
	{return a.numeric_derivative_par > b.numeric_derivative_par; };
	run_list.sort(compare);
	auto &run_first = run_list.front();
	auto  &run_last = run_list.back();

	//p_rec = group_info.get_parameter_rec_ptr(*par_name);
	g_rec = group_info.get_group_rec(numeric_par_name);
	double splitthresh = g_rec.splitthresh;
	double splitreldiff = g_rec.splitreldiff;

	irow = 0;
	vector<double> sen_vec;
	for (auto &iobs_name : base_sim_obs_names)
	{
		// Check if this is not prior information
		if (prior_info.find(iobs_name) == prior_info.end())
		{
			//Apply Split threshold on derivative if applicable
			bool success = false;
			if (run_list.size() == 3 && splitswh_flag == true && splitswh_flag)
			{
				sen_vec.clear();
				list<JacobianRun>::const_iterator iter_run = run_list.begin();
				list<JacobianRun>::const_iterator run_2 = run_list.begin();
				for (++iter_run; iter_run != run_list.end(); ++iter_run)
				{
					list<JacobianRun>::const_iterator run_1 = run_2;
					run_2 = iter_run;
					del_par = run_2->numeric_derivative_par - run_1->numeric_derivative_par;
					del_obs = run_2->obs_vec[irow] - run_1->obs_vec[irow];
					sen_vec.push_back(del_obs / del_par);
				}
				std::sort(sen_vec.begin(), sen_vec.end(), [](double a, double b) {
					return std::abs(a) < std::abs(b);});
				if (abs(sen_vec.back()) >= splitthresh &&
					abs(sen_vec.back() - sen_vec.front()) / sen_vec.front() > splitreldiff )
				{
					success = true;
					if (sen_vec.front() != 0)
					{
						triplet_list.push_back(Eigen::Triplet<double>(irow, jcol, sen_vec.front() ));
					}
				}
			}

			if (run_list.size() == 3 && g_rec.dermthd == "PARABOLIC" && !success)
			{
				// Central Difference Parabola
				// Solve Ac = o for c to get the equation for a parabola where:
				//        | p0**2  p0  1 |               | c0 |            | o0 |
				//   A =  | p1**2  p1  1 |          c =  | c1 |        y = | o1 |
				//        | p2**2  p2  1 |               | c2 |            | o2 |
				// then compute the derivative as:
				//   dy/dx = 2 * c0 * base_numeric_par_value + c1
				auto &run2 = (*(++run_list.begin()));
				MatrixXd a_mat(3,3);
				VectorXd c(3), y(3);
				// assemble A matrix
				int i=0;
				for (auto &irun_pair : run_list)
				{
					// assemble A matrix
					double par_value = irun_pair.numeric_derivative_par;
					a_mat(i,0) = par_value*par_value; a_mat(i,1) = par_value; a_mat(i,2) = 1;
					// assemble y vector
					y(i) =  irun_pair.obs_vec[irow];
					++i;
				}
				c = a_mat.colPivHouseholderQr().solve(y);
				//derivative is calculated around "base_numeric_par_value"
				der = 2.0 * c(0) *  base_numeric_par_value + c(1);
				if (der != 0)
				{
					triplet_list.push_back(Eigen::Triplet<double>(irow, jcol, der));
				}
			}
			else if (!success)
			{
				// Forward Difference and Central Difference Outer
				del_par = run_last.numeric_derivative_par - run_first.numeric_derivative_par;
				del_obs = run_last.obs_vec[irow] - run_first.obs_vec[irow];
				if (del_obs != 0)
				{
					triplet_list.push_back(Eigen::Triplet<double>(irow, jcol, del_obs / del_par));
				}
			}
		}
		else
		{
			// Prior Information always calculated using outer model runs even for central difference
			del_par = run_last.numeric_derivative_par - run_first.numeric_derivative_par;
			double del_prior_info;

			const PriorInformationRec *pi_rec;
			const Parameters &ctl_pars_1 = run_first.ctl_pars;
			const Parameters &ctl_pars_2 = run_last.ctl_pars;
			const auto prior_info_it = prior_info.find(iobs_name);
			if (prior_info_it != prior_info.end())
			{
				pi_rec = &(prior_info_it->second);
				del_prior_info = pi_rec->calc_residual(ctl_pars_2) - pi_rec->calc_residual(ctl_pars_1);
				if (del_prior_info != 0) {
					triplet_list.push_back(Eigen::Triplet<double>(irow, jcol, del_prior_info / del_par));
				}
			}
		}
		++irow;
	}
	return triplet_list;
}


bool Jacobian::forward_diff(const string &par_name, const Parameters &numeric_parameters,
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info,
		const ParamTransformSeq &par_trans, double &new_par, set<string> &out_of_bound_par)
{
	// if the transformation is one to one, call the simpler and more effiecient version of this routine
	// designed specifically for this case
	bool out_of_bound_forward;
	bool out_of_bound_backward;
	string tmp_name;

	double incr = derivative_inc(par_name, group_info, numeric_parameters.get_rec(par_name), false);
	// try forward derivative
	// perturb numeric paramateres
	Parameters numeric_derivative_pars(numeric_parameters);
	new_par = numeric_derivative_pars[par_name] += incr;
	out_of_bound_forward = out_of_bounds(par_trans.numeric2ctl_cp(numeric_derivative_pars), ctl_par_info, out_of_bound_par);
	if (!out_of_bound_forward) {
		return true;
	}
	// try backward derivative if forward derivative didn't work
	numeric_derivative_pars = numeric_parameters;
	new_par = numeric_derivative_pars[par_name] -= incr;
	out_of_bound_backward = out_of_bounds(par_trans.numeric2ctl_cp(numeric_derivative_pars), ctl_par_info, out_of_bound_par);
	if (!out_of_bound_backward)
	{
		return true;
	}
	return false;
}


bool Jacobian::central_diff(const string &par_name, const Parameters &pest_parameters,
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, const ParamTransformSeq &par_trans, vector<double> &new_par_vec,
		vector<Parameters>  &numeric_dir_par_vec, set<string> &out_of_bound_par)
{
	double new_par;
	bool out_of_bnds_forward, out_of_bnds_back, out_of_bnds;
	Parameters numeric_dir_pars;
	string tmp_name;

	double incr = derivative_inc(par_name, group_info, pest_parameters.get_rec(par_name), true);
	// try backward difference
	numeric_dir_pars = pest_parameters;
	new_par = numeric_dir_pars[par_name] -= incr;
	out_of_bnds_back = out_of_bounds(par_trans.numeric2ctl_cp(numeric_dir_pars), ctl_par_info, out_of_bound_par);

	if (!out_of_bnds_back) {
		new_par_vec.push_back(new_par);
		numeric_dir_par_vec.push_back(numeric_dir_pars);
	}
	// try forward derivative
	numeric_dir_pars = pest_parameters;
	new_par = numeric_dir_pars[par_name] += incr;
	out_of_bnds_forward = out_of_bounds(par_trans.numeric2ctl_cp(numeric_dir_pars), ctl_par_info, out_of_bound_par);
	if (!out_of_bnds_forward) {
		new_par_vec.push_back(new_par);
		numeric_dir_par_vec.push_back(numeric_dir_pars);
	}
	// if backward difference was out of bounds do a second forward derivative
	if (out_of_bnds_back) {
		set<string> tmp_out_of_bound_par;
		numeric_dir_pars = pest_parameters;
		new_par = numeric_dir_pars[par_name] += 2.0 * incr;
		out_of_bnds = out_of_bounds(par_trans.numeric2ctl_cp(numeric_dir_pars), ctl_par_info, tmp_out_of_bound_par);
		if (!out_of_bnds) {
			new_par_vec.push_back(new_par);
			numeric_dir_par_vec.push_back(numeric_dir_pars);
		}
		else
		{
			return false;  // can't do central difference without going out of bounds
		}
	}
	// if forward difference was out of bounds do a second backward derivative
	if (out_of_bnds_forward) {
		set<string> tmp_out_of_bound_par;
		numeric_dir_pars = pest_parameters;
		new_par = numeric_dir_pars[par_name] -= 2.0 * incr;
		out_of_bnds = out_of_bounds(par_trans.numeric2ctl_cp(numeric_dir_pars), ctl_par_info, tmp_out_of_bound_par);
		if (!out_of_bnds) {
			new_par_vec.insert(new_par_vec.begin(), new_par);
			numeric_dir_par_vec.insert(numeric_dir_par_vec.begin(), numeric_dir_pars);
		}
		else
		{
			return false;  // can't do central difference without going out of bounds
		}
	}
	return true;
}

bool Jacobian::out_of_bounds(const Parameters &ctl_parameters,
	const ParameterInfo &ctl_par_info, set<string> &out_of_bound_par) const
{
	const string *par_name;
	double min, max;
	const ParameterRec *par_info_ptr;
	bool out_of_bounds=false;
	double bnd_tol = 0.001;
	for(Parameters::const_iterator b=ctl_parameters.begin(), e=ctl_parameters.end();
		b!=e; ++b) {
			par_name = &(*b).first;
			par_info_ptr = ctl_par_info.get_parameter_rec_ptr(*par_name);
			max = par_info_ptr->ubnd + abs(par_info_ptr->ubnd * bnd_tol);
			min = par_info_ptr->lbnd - abs(par_info_ptr->lbnd * bnd_tol);
		if ((*b).second > max || (*b).second < min) {
			out_of_bounds = true;
			out_of_bound_par.insert(*par_name);
		}
	}
	return out_of_bounds;
}

double Jacobian::derivative_inc(const string &name, const ParameterGroupInfo &group_info, double cur_par_value, bool central)
{
	ParameterGroupRec g_rec;
	double incr = 0.0;

	//// to do add error checking
	g_rec = group_info.get_group_rec(name);
	if (g_rec.inctyp == "ABSOLUTE") {
		incr = g_rec.derinc;}
	else if (g_rec.inctyp == "RELATIVE") {
		incr =  g_rec.derinc * abs(cur_par_value);
	}
	// apply derincmul for central derivatives
	if (central) {
		incr *= g_rec.derincmul;
	}
	// apply lower bound
	if ((g_rec.inctyp != "ABSOLUTE") && (incr < g_rec.derinclb)) {
		incr = g_rec.derinclb;
	}
	return incr;
}

const set<string>& Jacobian::get_failed_parameter_names() const
{
	return failed_parameter_names;
}

Jacobian& Jacobian::operator=(const Jacobian &rhs)
{
	base_numeric_par_names = rhs.base_numeric_par_names;
	base_numeric_parameters = rhs.base_numeric_parameters;
	failed_parameter_names = rhs.failed_parameter_names;
	base_sim_obs_names = rhs.base_sim_obs_names;
	base_sim_observations = rhs.base_sim_observations;
	matrix = rhs.matrix;
	file_manager = rhs.file_manager;
	return *this;
}
void Jacobian::transform(const ParamTransformSeq &par_trans, void(ParamTransformSeq::*meth_prt)(Jacobian &jac) const)
{
	(par_trans.*meth_prt)(*this);
}

void Jacobian::print(std::ostream &fout) const
{
	fout << "Jacobian:" << endl;
	fout << "base_numeric_par_names: " << base_numeric_par_names << endl;
	fout << "base_numeric_parameters: " << base_numeric_parameters << endl;
	fout << "failed_parameter_names: " << failed_parameter_names << endl;
	fout << "base_sim_obs_names: " << base_sim_obs_names << endl;
	fout << "base_sim_observations: " << base_sim_observations << endl;
	fout << "matrix: " << matrix << endl;
}

//void Jacobian::save_old(const string &ext) const
//{
//	ofstream &jout = file_manager.open_ofile_ext(ext, ios::out |ios::binary);
//	int n_par = base_numeric_par_names.size();
//	int n_obs_and_pi = base_sim_obs_names.size();
//	int n;
//	int tmp;
//	double data;
//	char par_name[12];
//	char obs_name[20];
//
//	// write header
//	tmp  = -n_par;
//	jout.write((char*) &tmp, sizeof(tmp));
//	tmp = -n_obs_and_pi;
//	jout.write((char*) &tmp, sizeof(tmp));
//
//	//write number nonzero elements in jacobian (includes prior information)
//	n = matrix.nonZeros();
//	jout.write((char*)&n, sizeof(n));
//
//	//write matrix
//	n = 0;
//	map<string, double>::const_iterator found_pi_par;
//	map<string, double>::const_iterator not_found_pi_par;
//
//	Eigen::SparseMatrix<double> matrix_T(matrix);
//	matrix_T.transpose();
//	for (int icol=0; icol<matrix.outerSize(); ++icol)
//	{
//		for (SparseMatrix<double>::InnerIterator it(matrix_T, icol); it; ++it)
//		{
//			data = it.value();
//			n = it.row() + 1 + it.col() * matrix_T.rows();
//			jout.write((char*) &(n), sizeof(n));
//			jout.write((char*) &(data), sizeof(data));
//			}
//	}
//	//save parameter names
//	for(vector<string>::const_iterator b=base_numeric_par_names.begin(), e=base_numeric_par_names.end();
//		b!=e; ++b) {
//		string l = lower_cp(*b);
//		string_to_fortran_char(l, par_name, 12);
//		jout.write(par_name, 12);
//	}
//
//	//save observation and Prior information names
//	for(vector<string>::const_iterator b=base_sim_obs_names.begin(), e=base_sim_obs_names.end();
//		b!=e; ++b) {
//		string l = lower_cp(*b);
//		string_to_fortran_char(l, obs_name, 20);
//		jout.write(obs_name, 20);
//	}
//	//save observation names (part 2 prior information)
//	file_manager.close_file(ext);
//}

void Jacobian::save(const string &ext) const
{
	string filename = file_manager.build_filename(ext);
	pest_utils::save_binary(filename,  base_sim_obs_names, base_numeric_par_names, matrix);

	//ofstream &jout = file_manager.open_ofile_ext(ext, ios::out | ios::binary);
	//int n_par = base_numeric_par_names.size();
	//int n_obs_and_pi = base_sim_obs_names.size();

	//int n;
	//int tmp;
	//double data;
	//char par_name[200];
	//char obs_name[200];

	//// write header
	//tmp = n_par;
	//jout.write((char*)&tmp, sizeof(tmp));
	//tmp = n_obs_and_pi;
	//jout.write((char*)&tmp, sizeof(tmp));

	////write number nonzero elements in jacobian (includes prior information)
	//n = matrix.nonZeros();
	//jout.write((char*)&n, sizeof(n));

	////write matrix
	//n = 0;
	//map<string, double>::const_iterator found_pi_par;
	//map<string, double>::const_iterator not_found_pi_par;

	//Eigen::SparseMatrix<double> matrix_T(matrix);
	//matrix_T.transpose();
	//for (int icol = 0; icol<matrix.outerSize(); ++icol)
	//{
	//	for (Eigen::SparseMatrix<double>::InnerIterator it(matrix_T, icol); it; ++it)
	//	{
	//		data = it.value();
	//		n = it.row() - 1;
	//		jout.write((char*) &(n), sizeof(n));
	//		n = it.col() - 1;
	//		jout.write((char*) &(n), sizeof(n));

	//		jout.write((char*) &(data), sizeof(data));
	//	}
	//}
	////save parameter names
	//for (vector<string>::const_iterator b = base_numeric_par_names.begin(), e = base_numeric_par_names.end();
	//	b != e; ++b) {
	//	string l = lower_cp(*b);
	//	string_to_fortran_char(l, par_name, 200);
	//	jout.write(par_name, 200);
	//}

	////save observation and Prior information names
	//for (vector<string>::const_iterator b = base_sim_obs_names.begin(), e = base_sim_obs_names.end();
	//	b != e; ++b) {
	//	string l = lower_cp(*b);
	//	string_to_fortran_char(l, obs_name, 200);
	//	jout.write(obs_name, 200);
	//}
	////save observation names (part 2 prior information)
	//file_manager.close_file(ext);
}

void Jacobian::read(const string &filename)
{
	pest_utils::read_binary(filename,base_sim_obs_names, base_numeric_par_names, matrix);
	//ifstream fin;
	//fin.open(filename.c_str(), ifstream::binary|ios::in);

	//if (!fin)
	//{
	//	throw runtime_error("unable to open binary jacobian file: " + filename + " for reading");
	//}

	//int n_par;
	//int n_nonzero;
	//int n_obs_and_pi;
	//int i,j,n;
	//double data;
	//char par_name[12];
	//char obs_name[20];

	//// read header
	//fin.read((char*) &n_par, sizeof(n_par));
	//fin.read((char*) &n_obs_and_pi, sizeof(n_obs_and_pi));
	//n_par = -n_par;
	//n_obs_and_pi = -n_obs_and_pi;
	//////read number nonzero elements in jacobian (observations + prior information)
	//fin.read((char*)&n_nonzero, sizeof(n_nonzero));

	//file_manager.get_ofstream("rec") << "  -->reading " << n_nonzero <<
	//	" jacobian elements for " << n_par << " parameters and" << endl <<
	//	"  --> " << n_obs_and_pi << " observations and prior info" << endl;

	//cout << "  -->reading " << n_nonzero <<
	//	" jacobian elements for " << n_par << " parameters and" << endl <<
	//	"  --> " << n_obs_and_pi << " observations and prior info" << endl;

	//// record current position in file

	//streampos begin_sen_pos = fin.tellg();

	////advance to parameter names section
	//fin.seekg(n_nonzero*(sizeof(double)+sizeof(int)), ios_base::cur);

	////read parameter names
	//base_numeric_par_names.clear();
	//for (int i_rec=0; i_rec<n_par; ++i_rec)
	//{
	//	fin.read(par_name, 12);
	//	string temp_par = string(par_name, 12);
	//	strip_ip(temp_par);
	//	upper_ip(temp_par);
	//	base_numeric_par_names.push_back(temp_par);
	//}
	////read observation and Prior info names
	//base_sim_obs_names.clear();

	//for (int i_rec=0; i_rec<n_obs_and_pi; ++i_rec)
	//{
	//	fin.read(obs_name, 20);
	//	string tmp_obs_name = strip_cp(string(obs_name, 20));
	//	upper_ip(tmp_obs_name);
	//	base_sim_obs_names.push_back(tmp_obs_name);
	//}

	////return to sensitivity section of file
	//fin.seekg(begin_sen_pos, ios_base::beg);

	//// read matrix
	//std::vector<Eigen::Triplet<double> > triplet_list;
	//triplet_list.reserve(n_nonzero);
	//for (int i_rec=0; i_rec<n_nonzero; ++ i_rec)
	//{
	//	fin.read((char*) &(n), sizeof(n));
	//	--n;
	//	fin.read((char*) &(data), sizeof(data));
	//	j = int(n/(n_obs_and_pi)); // parameter index
	//	i = (n-n_obs_and_pi*j) % n_obs_and_pi;  //observation index
	//	triplet_list.push_back(Eigen::Triplet<double>(i, j, data));
	//}
	//matrix.resize(n_obs_and_pi, n_par);
	//matrix.setZero();
	//matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
	//fin.close();

	//cout << "   ....done" << endl;
	//file_manager.get_ofstream("rec") << "   ....done" << endl;

}

void Jacobian::report_errors(std::ostream &fout)
{
	if (failed_parameter_names.size() > 0)
	{
		fout << "    Parameters whose perturbation runs failed while computing jacobian: " << endl;

		for (const auto & ipar : failed_parameter_names)
		{
			fout << right << "  " << setw(12) << ipar << endl;
		}
		cout << "WARNING:  " << failed_parameter_names.size() << " parameter perturbation runs failed while computing jacobian, see rec file for listing" << endl;
	}

}

Eigen::SparseMatrix<double>* Jacobian::get_matrix_ptr()
{
	Eigen::SparseMatrix<double>* ptr = &matrix;
	return ptr;
}

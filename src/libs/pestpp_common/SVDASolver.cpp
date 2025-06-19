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
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include "SVDASolver.h"
#include "ModelRunPP.h"
#include "QSqrtMatrix.h"
#include "eigen_tools.h"
#include "ObjectiveFunc.h"
#include "utilities.h"
#include "FileManager.h"
#include "TerminationController.h"
#include "ParamTransformSeq.h"
#include "Transformation.h"
#include "PriorInformation.h"
#include "debug.h"
#include "covariance.h"
#include "Ensemble.h"
#include "linear_analysis.h"

using namespace std;
using namespace pest_utils;
using namespace Eigen;

const string SVDASolver::svda_solver_type_name = "svda_super_par";
/*
SVDSolver(Pest &_pest_scenario, FileManager &_file_manager, ObjectiveFunc *_obj_func,
		const ParamTransformSeq &_par_transform, Jacobian &_jacobian,
		OutputFileWriter &_output_file_writer, SVDSolver::MAT_INV _mat_inv,
		PerformanceLog *_performance_log, const string &description = string("base parameter solution"),
		bool _phiredswh_flag = false, bool _splitswh_flag = false, bool _save_next_jacobian = true);
*/
SVDASolver::SVDASolver(Pest &_pest_scenario, FileManager &_file_manager, ObjectiveFunc *_obj_func_ptr,
	const ParamTransformSeq &_par_transform, Jacobian &_jacobian,	OutputFileWriter &_output_file_writer,
	PerformanceLog *_performance_log, Covariance& _parcov, std::mt19937* _rand_gen_ptr, bool _phiredswh_flag, bool _splitswh_flag)
	: SVDSolver(_pest_scenario, _file_manager, _obj_func_ptr, _par_transform, _jacobian,
		_output_file_writer, _performance_log,_parcov,_rand_gen_ptr, "super parameter solution", _phiredswh_flag, _splitswh_flag, false),
		max_super_frz_iter(_pest_scenario.get_pestpp_options().get_max_super_frz_iter())
{
	//PestppOptions::GLMNormalForm nf = pest_scenario.get_pestpp_options().get_glm_normal_form();
	glm_normal_form = pest_scenario.get_pestpp_options().get_glm_normal_form();
	if (glm_normal_form == PestppOptions::GLMNormalForm::PRIOR)
	{
		file_manager.rec_ofstream() << "Note: using 'GLM_NORMAL_FORM' = 'DIAG' for super iterations" << endl;
		glm_normal_form = PestppOptions::GLMNormalForm::DIAG;
	}
}


Parameters SVDASolver::limit_parameters_freeze_all_ip(const Parameters &init_active_ctl_pars,
		Parameters &upgrade_active_ctl_pars, const Parameters &prev_frozen_active_ctl_pars)
{
	const string *name;
	double val_init;
	double val_upgrade;

	// superparameters are always limited by RELPARMAX
	// if any base parameters go outof bounds, they should be frozen,
	// a new SVD factorization should be computed, and then this routine should be called again

	//Assume Super parameters are independent and apply limits one at a time

	par_transform.active_ctl2numeric_ip(upgrade_active_ctl_pars);
	Parameters init_numeric_pars = par_transform.active_ctl2numeric_cp(init_active_ctl_pars);

	for (auto &ipar : upgrade_active_ctl_pars)
	{
		name = &(ipar.first);  // parameter name
		val_upgrade = ipar.second; // initial parameter value
		auto iter_pu = init_numeric_pars.find(*name);
		assert(iter_pu != init_numeric_pars.end());
		val_init = (*iter_pu).second;  // upgrade parameter value
		// only use relative change limits
		if (abs((val_upgrade - val_init) / val_init) > ctl_info->relparmax)
		{
			double delta = abs(ctl_info->relparmax * val_init) * sign(val_upgrade - val_init);
			ipar.second = val_init + delta;
		}
	}

	//convert parameters to their ctl form and check any any that have exceeded their bounds
	par_transform.numeric2active_ctl_ip(upgrade_active_ctl_pars);
	Parameters freeze_active_ctl_par;
	for (auto &ipar : upgrade_active_ctl_pars)
	{
		name = &(ipar.first);
		const ParameterRec *p_info = ctl_par_info_ptr->get_parameter_rec_ptr(*name);
		if (prev_frozen_active_ctl_pars.find(*name) == prev_frozen_active_ctl_pars.end())
		{
			const ParameterRec *p_info = ctl_par_info_ptr->get_parameter_rec_ptr(*name);
			if (ipar.second > p_info->ubnd)
			{
				ipar.second = freeze_active_ctl_par[*name] = p_info->ubnd;
			}
			// Check parameter lower bound
			else if (ipar.second < p_info->lbnd)
			{
				ipar.second = freeze_active_ctl_par[*name] = p_info->lbnd;
			}
		}
	}
	for (auto &ipar : prev_frozen_active_ctl_pars)
	{
		upgrade_active_ctl_pars[ipar.first] = ipar.second;
	}

	return freeze_active_ctl_par;
}

void SVDASolver::calc_upgrade_vec(double i_lambda, Parameters &prev_frozen_active_ctl_pars, QSqrtMatrix &Q_sqrt,
	const DynamicRegularization &regul, VectorXd &residuals_vec,
	vector<string> &obs_names_vec, const Parameters &base_run_active_ctl_pars, Parameters &upgrade_active_ctl_pars,
	Pest::LimitType &limit_type)
{
	Parameters upgrade_ctl_del_pars;
	Parameters grad_ctl_del_pars;
	int num_upgrade_out_grad_in;
	Parameters new_frozen_ctl_pars;

	// define a function type for upgrade methods
	typedef void(SVDSolver::*UPGRADE_FUNCTION) (const Jacobian &jacobian, const QSqrtMatrix &Q_sqrt, const DynamicRegularization &regul,
		const Eigen::VectorXd &Residuals, const vector<string> &obs_name_vec,
		const Parameters &base_ctl_pars, const Parameters &prev_frozen_ctl_pars,
		double lambda, Parameters &ctl_upgrade_pars, Parameters &upgrade_ctl_del_pars,
		Parameters &grad_ctl_del_pars);

	UPGRADE_FUNCTION calc_lambda_upgrade = &SVDASolver::calc_lambda_upgrade_vec_JtQJ;
	

		// need to remove parameters frozen due to failed jacobian runs when calling calc_lambda_upgrade_vec
		//Freeze Parameters at the boundary whose upgrade vector and gradient both head out of bounds
	(*this.*calc_lambda_upgrade)(jacobian, Q_sqrt, regul, residuals_vec, obs_names_vec,
			base_run_active_ctl_pars, prev_frozen_active_ctl_pars, i_lambda, upgrade_active_ctl_pars, upgrade_ctl_del_pars,
			grad_ctl_del_pars);
		num_upgrade_out_grad_in = check_bnd_par(new_frozen_ctl_pars, base_run_active_ctl_pars, upgrade_ctl_del_pars, grad_ctl_del_pars);
		prev_frozen_active_ctl_pars.insert(new_frozen_ctl_pars.begin(), new_frozen_ctl_pars.end());
		//Recompute the upgrade vector without the newly frozen parameters and freeze those at the boundary whose upgrade still goes heads out of bounds
		if (num_upgrade_out_grad_in > 0)
		{
			new_frozen_ctl_pars.clear();
			(*this.*calc_lambda_upgrade)(jacobian, Q_sqrt, regul, residuals_vec, obs_names_vec,
				base_run_active_ctl_pars, prev_frozen_active_ctl_pars, i_lambda, upgrade_active_ctl_pars, upgrade_ctl_del_pars,
				grad_ctl_del_pars);
			check_bnd_par(new_frozen_ctl_pars, base_run_active_ctl_pars, upgrade_active_ctl_pars);
			prev_frozen_active_ctl_pars.insert(new_frozen_ctl_pars.begin(), new_frozen_ctl_pars.end());
			new_frozen_ctl_pars.clear();
		}
		//If there are newly frozen parameters recompute the upgrade vector
		if (new_frozen_ctl_pars.size() > 0)
		{
			(*this.*calc_lambda_upgrade)(jacobian, Q_sqrt, regul, residuals_vec, obs_names_vec,
				base_run_active_ctl_pars, prev_frozen_active_ctl_pars, i_lambda, upgrade_active_ctl_pars, upgrade_ctl_del_pars,
				grad_ctl_del_pars);
		}
		//Freeze any new parameters that want to go out of bounds
		new_frozen_ctl_pars.clear();
		new_frozen_ctl_pars = limit_parameters_freeze_all_ip(base_run_active_ctl_pars, upgrade_active_ctl_pars, prev_frozen_active_ctl_pars);
		prev_frozen_active_ctl_pars.insert(new_frozen_ctl_pars.begin(), new_frozen_ctl_pars.end());
}


ModelRun SVDASolver::update_run(RunManagerAbstract &run_manager, ModelRun &base_run)
{
	ModelRun new_base_run = base_run;
	Parameters base_ctl_pars = new_base_run.get_ctl_pars();
	// make sure these are all in bounds
	for (auto &ipar : base_ctl_pars)
	{
		const string &name = ipar.first;
		const ParameterRec *p_info = ctl_par_info_ptr->get_parameter_rec_ptr(name);
		if (ipar.second > p_info->ubnd)
		{
			ipar.second = p_info->ubnd;
		}
		// Check parameter lower bound
		else if (ipar.second < p_info->lbnd)
		{
			ipar.second = p_info->lbnd;
		}
	}

	run_manager.reinitialize(file_manager.build_filename("rnr"));
	int run_id = run_manager.add_run(par_transform.ctl2model_cp(new_base_run.get_ctl_pars()));

	run_manager.run();
	Parameters tmp_pars;
	Observations tmp_obs;
	bool success = run_manager.get_run(run_id, tmp_pars, tmp_obs);
	if (!success)
	{
		throw(PestError("Error: Base super parameter run failed."));
	}
	par_transform.model2ctl_ip(tmp_pars);
	new_base_run.update_ctl(tmp_pars, tmp_obs);
	return new_base_run;
}


ModelRun SVDASolver::iteration_reuse_jac(RunManagerAbstract &run_manager, TerminationController &termination_ctl, ModelRun &base_run, bool rerun_base, const string &filename)
{
	ModelRun new_base_run = base_run;

	ostream &fout_restart = file_manager.get_ofstream("rst");
	ostream &os = file_manager.rec_ofstream();
	vector<string> numeric_par_names_vec;



	string jac_filename = filename;
	if (filename.empty()) jac_filename = file_manager.build_filename("jcs");
	cout << "  reading previously computed jacobian: " << jac_filename << endl;
	file_manager.get_ofstream("rec") << "  reading previously computed jacobian: " << jac_filename << endl;
	jacobian.read(jac_filename);

	if (rerun_base)
	{
		new_base_run = update_run(run_manager, base_run);
	}

	//jacobian.save("jcs");
	output_file_writer.write_jco(false, "jcs", jacobian);
	// sen file for this iteration
	output_file_writer.append_sen(file_manager.sen_ofstream(), termination_ctl.get_iteration_number() + 1,
		jacobian, *(new_base_run.get_obj_func_ptr()), get_parameter_group_info(), *regul_scheme_ptr, true, par_transform);
	cout << endl;
	return new_base_run;
}

bool SVDASolver::iteration_jac(RunManagerAbstract &run_manager, TerminationController &termination_ctl, ModelRun &base_run, bool calc_init_obs, bool restart_runs)
{
	ostream &fout_restart = file_manager.get_ofstream("rst");
	ostream &os = file_manager.rec_ofstream();
	vector<string> numeric_par_names_vec;

	// save state of termination controller
	//termination_ctl.save_state(fout_restart);

	Parameters base_ctl_pars = base_run.get_ctl_pars();
	// make sure these are all in bounds
	for (auto &ipar : base_ctl_pars)
	{
		const string &name = ipar.first;
		const ParameterRec *p_info = ctl_par_info_ptr->get_parameter_rec_ptr(name);
		if (ipar.second > p_info->ubnd)
		{
			ipar.second = p_info->ubnd;
		}
		// Check parameter lower bound
		else if (ipar.second < p_info->lbnd)
		{
			ipar.second = p_info->lbnd;
		}
	}

	if (restart_runs)
	{
		super_parameter_group_info = par_transform.get_svda_ptr()->build_par_group_info(*par_group_info_ptr);
		par_transform.get_svda_fixed_ptr()->reset(par_transform.get_svda_ptr()->get_frozen_derivative_pars());
	}
	else
	{
		set<string> out_of_bound_pars;
		// Calculate Jacobian
		// build model runs
		performance_log->log_event("commencing to build jacobian parameter sets");
		bool success_build_runs = false;
		cout << "  calculating jacobian... ";
		int n_freeze_iter = 0;
		while (true) //loop and feeeze any base parameters that go out of bounds when computing the jacobian
		{
			++n_freeze_iter;
			if (n_freeze_iter > max_super_frz_iter)
			{
				terminate_local_iteration = true;
				cout << "Terminating super parameter iterations." << endl;
				cout << "Max number of iterations to freeze parameters to compute jacobian exceeded" << endl;
				os << "Terminating super parameter iterations." << endl;
				os << "Max number of iterations to freeze parameters to compute jacobian exceeded" << endl;
				return false;
			}
			// fix frozen parameters in SVDA transformation
			//check one more time for base parameters at bounds...
			map<string, double> at_bounds = pest_scenario.get_pars_at_near_bounds(base_run.get_ctl_pars());
			if (at_bounds.size() > 0)
			{
				Parameters frz;
				for (auto item : at_bounds)
				{
					frz.insert(item.first, item.second);
				}
				base_run.add_frozen_ctl_parameters(frz);
			}
			debug_print(base_run.get_frozen_ctl_pars());
			
			par_transform.get_svda_ptr()->update_add_frozen_pars(base_run.get_frozen_ctl_pars());
			par_transform.get_svda_fixed_ptr()->reset(par_transform.get_svda_ptr()->get_frozen_derivative_pars());
			Parameters numeric_pars = par_transform.ctl2numeric_cp(base_run.get_ctl_pars());
			numeric_par_names_vec = numeric_pars.get_keys();

			calc_init_obs = true;

			super_parameter_group_info = par_transform.get_svda_ptr()->build_par_group_info(*par_group_info_ptr);
			performance_log->log_event("commencing to build jacobian parameter sets");
			out_of_bound_pars.clear();
			success_build_runs = jacobian.build_runs(base_run, numeric_par_names_vec, par_transform,
				super_parameter_group_info, *ctl_par_info_ptr, run_manager, out_of_bound_pars,
				phiredswh_flag, calc_init_obs);
			if (success_build_runs)
			{
				break;
			}
			else if (!success_build_runs && out_of_bound_pars.size() > 0)
			{

				cout << "  can not compute super-par jacobian without the following parameters going out of bounds... " << endl;
				os << "  can not compute super-par jacobian without the following parameters going out of bounds... " << endl;
				print(out_of_bound_pars, cout, 4);
				print(out_of_bound_pars, os, 4);

				cout << "  freezing parameters and recalculating the jacobian" << endl;
				os << "  freezing parameters and recalculating the jacobian" << endl;
				//add out of bound parameters to frozen parameter list
				Parameters new_frz_derivative_pars;
				for (auto &ipar : out_of_bound_pars)
				{
					const auto iter = base_ctl_pars.find(ipar);
					assert(iter != base_ctl_pars.end());
					new_frz_derivative_pars.insert(iter->first, iter->second);
				}
				if (new_frz_derivative_pars.size() > 0)
				{
					Parameters frz_pars = base_run.get_frozen_ctl_pars();
					for (auto p : new_frz_derivative_pars)
					{
						if (frz_pars.find(p.first) != frz_pars.end())
						{
							cout << "SDVA::iteration_jac() WARNING: parameter " << p.first << " already in frozen par names" << endl;
							os << "SDVA::iteration_jac() WARNING: parameter " << p.first << " already in frozen par names" << endl;
						}
					}

					base_run.add_frozen_ctl_parameters(new_frz_derivative_pars);
				}
				else
				{

				}
			}
			else
			{
				break;
			}
		}

		if (!success_build_runs)
		{
			throw PestError("Error in SVDASolver::iteration: Can not compute super parameter derivatives without base parameters going out of bounds");
		}
		performance_log->log_event("jacobian parameter sets built, commencing model runs");

	}
	// save super parameter transformation
	ofstream &fout_rst = file_manager.open_ofile_ext("rtj", ios_base::out | ios_base::binary);
	par_transform.get_svda_ptr()->save(fout_rst);
	file_manager.close_file("rtj");
	RestartController::write_jac_runs_built(fout_restart);
	//make model runs
	jacobian.make_runs(run_manager);
	performance_log->log_event("jacobian runs complete, processing runs");
	bool success_process_runs = jacobian.process_runs(par_transform,
		super_parameter_group_info, run_manager, *prior_info_ptr, splitswh_flag,
		pest_scenario.get_pestpp_options().get_glm_debug_der_fail());
	if (!success_process_runs)
	{
		throw PestError("Error in SVDASolver::iteration: Can not compute super parameter derivatives");
	}
	performance_log->log_event("saving jacobian and sen files");
	// save jacobian
	//jacobian.save("jcs");
	output_file_writer.write_jco(false,"jcs", jacobian);

	//Update parameters and observations for base run
	{
		Parameters tmp_pars;
		Observations tmp_obs;
		bool success = run_manager.get_run(0, tmp_pars, tmp_obs);
		par_transform.model2ctl_ip(tmp_pars);
		base_run.update_ctl(tmp_pars, tmp_obs);
	}

	// sen file for this iteration
	output_file_writer.append_sen(file_manager.sen_ofstream(), termination_ctl.get_iteration_number() + 1,
		jacobian, *(base_run.get_obj_func_ptr()), get_parameter_group_info(), *regul_scheme_ptr, true, par_transform);
	return true;
}

ModelRun SVDASolver::iteration_upgrd(RunManagerAbstract &run_manager, TerminationController &termination_ctl, ModelRun &base_run, bool restart_runs)
{
	ostream &os = file_manager.rec_ofstream();
	ostream &fout_restart = file_manager.get_ofstream("rst");
	int num_lam_runs = 0;
	if (restart_runs)
	{
		run_manager.initialize_restart(file_manager.build_filename("rnu"));
		Parameters tmp_pars;
		Observations tmp_obs;
		bool success = run_manager.get_run(0, tmp_pars, tmp_obs);
		if (!success)
		{
			throw(PestError("Error: Cannot retrieve the base run to compute upgrade vectors."));
		}
		par_transform.model2ctl_ip(tmp_pars);
		base_run.update_ctl(tmp_pars, tmp_obs);
	}
	else
	{
		vector<string> obs_names_vec = base_run.get_obs_template().get_keys();
		performance_log->log_event("computing upgrade vectors");

		// populate vectors with sorted observations (standard and prior info) and parameters
		{
			//vector<string> prior_info_names = prior_info_ptr->get_keys();
			//obs_names_vec.insert(obs_names_vec.end(), prior_info_names.begin(), prior_info_names.end());
		}

		//build residuals vector
		VectorXd residuals_vec = -1.0 * stlvec_2_eigenvec(base_run.get_residuals_vec(obs_names_vec));

		//Build model runs
		run_manager.reinitialize(file_manager.build_filename("rnu"));

		// Save base run as first model run so it is eassily accessible
		Parameters base_model_pars = par_transform.ctl2model_cp(base_run.get_ctl_pars());
		int run_id = run_manager.add_run(base_model_pars, "base_run");
		num_lam_runs++;
		run_manager.update_run(run_id, base_model_pars, base_run.get_obs());

		//Marquardt Lambda Update Vector
		Upgrade ml_upgrade;

		// build weights matrix sqrt(Q)
		QSqrtMatrix Q_sqrt(obs_info_ptr, prior_info_ptr);
		{
			const Parameters base_run_active_ctl_pars = par_transform.ctl2active_ctl_cp(base_run.get_ctl_pars());
			Parameters frozen_active_ctl_pars;
			//If running in regularization mode, adjust the regularization weights
			// define a function type for upgrade methods
			if (regul_scheme_ptr->get_use_dynamic_reg())
			{
				dynamic_weight_adj(base_run, jacobian, Q_sqrt, residuals_vec, obs_names_vec,
					base_run_active_ctl_pars, frozen_active_ctl_pars);
			}
		}

		vector<double> lambda_vec = pest_scenario.get_pestpp_options().get_base_lambda_vec();
		std::sort(lambda_vec.begin(), lambda_vec.end());
		auto iter = std::unique(lambda_vec.begin(), lambda_vec.end());
		lambda_vec.resize(std::distance(lambda_vec.begin(), iter));
		stringstream message;
		stringstream prf_message;

		ofstream &fout_frz = file_manager.open_ofile_ext("fpr");
		ofstream& fout_rec = file_manager.rec_ofstream();
		int i_update_vec = 0;
		
		for (double i_lambda : lambda_vec)
		{
			prf_message.str("");
			prf_message << "beginning upgrade vector calculations, lambda = " << i_lambda;
			performance_log->log_event(prf_message.str());
			fout_rec << "   ...calculating lambda upgrade vector: " << i_lambda << endl;

			//std::cout << string(message.str().size(), '\b');
			message.str("");
			message << "  computing upgrade vector (lambda = " << i_lambda << ")  " << ++i_update_vec << " / " << lambda_vec.size() << "             " << endl;
			std::cout << message.str();
			cout.flush();

			Parameters new_pars;
			const Parameters base_run_active_ctl_pars = par_transform.ctl2active_ctl_cp(base_run.get_ctl_pars());
			Parameters base_numeric_pars = par_transform.ctl2numeric_cp(base_run.get_ctl_pars());

			Parameters frzn_pars;
			Pest::LimitType limit_type;
			calc_upgrade_vec(i_lambda, frzn_pars, Q_sqrt, *regul_scheme_ptr, residuals_vec,
				obs_names_vec, base_run_active_ctl_pars,
				new_pars,limit_type);

			//transform new_pars to model parameters
			par_transform.active_ctl2model_ip(new_pars);
			int run_id = run_manager.add_run(new_pars, "upgrade_nrm", i_lambda);
			num_lam_runs++;
			output_file_writer.write_upgrade(termination_ctl.get_iteration_number(), 1, i_lambda, 1.0, new_pars);
			save_frozen_pars(fout_frz, frzn_pars, run_id);
			
			//Add Scaled Upgrade Vectors
			Parameters new_numeric_pars = par_transform.model2numeric_cp(new_pars);
			Parameters base_numeric_pars1 = par_transform.model2numeric_cp(base_model_pars);
			Parameters del_numeric_pars = new_numeric_pars - base_numeric_pars1;
			for (double i_scale : lambda_scale_vec)
			{
				if (i_scale >= 1.0)
					continue;
				Parameters scaled_pars = base_numeric_pars1 + del_numeric_pars * i_scale;
				
				Parameters scaled_ctl_pars = par_transform.numeric2ctl_cp(scaled_pars);
				output_file_writer.write_upgrade(termination_ctl.get_iteration_number(),
					0, i_lambda, i_scale, scaled_ctl_pars);

				stringstream ss;
				ss << "scale:" << std::fixed << std::setprecision(2) << i_scale << "";
				par_transform.numeric2model_ip(scaled_pars);
				int run_id = run_manager.add_run(scaled_pars, ss.str(), i_lambda);
				num_lam_runs++;
				fout_rec << "   ...calculating scaled lambda vector-scale factor: " << i_lambda << ", " << i_scale << endl;
				save_frozen_pars(fout_frz, frzn_pars, run_id);
			}


		}
		file_manager.close_file("fpr");
		RestartController::write_upgrade_runs_built(fout_restart);
	}

	//instance of a Mat for the jco
	pair<ParameterEnsemble, map<int, int>> fosm_real_info;
	Jacobian base_jco = jacobian; //copy
	par_transform.jac_numeric2active_ctl_ip(base_jco);
	pest_scenario.get_base_par_tran_seq().jac_active_ctl_ip2numeric_ip(base_jco);
	Mat j(base_jco.get_sim_obs_names(), base_jco.get_base_numeric_par_names(),
		base_jco.get_matrix_ptr());
	if (pest_scenario.get_prior_info_ptr()->get_nnz_pi() > 0)
	{
		vector<string> pi_names = pest_scenario.get_ctl_ordered_pi_names();
		j.drop_rows(pi_names);
	}
	LinearAnalysis la(j, pest_scenario, file_manager, *performance_log, parcov, rand_gen_ptr);
	if (pest_scenario.get_pestpp_options().get_uncert_flag())
	{
		cout << "-->starting iteration FOSM process..." << endl;

		performance_log->log_event("LinearAnalysis::glm_iter_fosm");

		try
		{
			la.glm_iter_fosm(base_run, output_file_writer, termination_ctl.get_iteration_number(), &run_manager);
			if (pest_scenario.get_pestpp_options().get_glm_iter_mc())
				fosm_real_info = la.draw_fosm_reals(&run_manager, termination_ctl.get_iteration_number(), base_run);
		}
		catch (exception& e)
		{
			os << "Error in GLM iteration FOSM process:" << e.what() << ", continuing..." << endl;

		}
		cout << "-->finished iteration FOSM process..." << endl;
	}



	cout << endl;

	// make upgrade model runs
	cout << "  performing upgrade vector model runs... ";
	cout.flush();
	run_manager.run();

	// process model runs
	cout << "  testing upgrade vectors... " << endl;
	ifstream &fin_frz = file_manager.open_ifile_ext("fpr");
	bool best_run_updated_flag = false;
	Parameters base_run_active_ctl_par_tmp = par_transform.ctl2active_ctl_cp(base_run.get_ctl_pars());
	ModelRun best_upgrade_run(base_run);

	os << "  Summary of upgrade runs:" << endl;
	Parameters new_frozen_pars;

	int n_runs = run_manager.get_nruns();
	bool one_success = false;
	//for (int i = 1; i < n_runs; ++i) {
	for (int i = 1; i < num_lam_runs; ++i) 
	{
		ModelRun upgrade_run(base_run);
		Parameters tmp_pars;
		Observations tmp_obs;
		string lambda_type;
		double i_lambda;
		bool success = run_manager.get_run(i, tmp_pars, tmp_obs, lambda_type, i_lambda);
		if ((pest_scenario.get_pestpp_options().get_glm_debug_lamb_fail()) && (i == 1))
		{
			file_manager.rec_ofstream() << "'GLM_DEBUG_LAMB_FAIL' is true, failing first lambda run" << endl;
			success = false;
		}
		if (success)
		{
			one_success = true;
			par_transform.model2ctl_ip(tmp_pars);
			upgrade_run.update_ctl(tmp_pars, tmp_obs);
			par_transform.ctl2active_ctl_ip(tmp_pars);
			double magnitude = Transformable::l2_norm(base_run_active_ctl_par_tmp, tmp_pars);

			streamsize n_prec = os.precision(2);
			os << "    Lambda = ";
			os << setiosflags(ios::fixed) << setw(8) << i_lambda;
			os << "; Type: " << setw(12) << lambda_type;
			os << "; length = " << std::scientific << magnitude;
			os << setiosflags(ios::fixed);
			os.precision(n_prec);
			os.unsetf(ios_base::floatfield); // reset all flags to default
			os << ";   phi = " << upgrade_run.get_phi(*regul_scheme_ptr);
			os.precision(2);
			os << setiosflags(ios::fixed);
			os << " (" << upgrade_run.get_phi(*regul_scheme_ptr) / base_run.get_phi(*regul_scheme_ptr) * 100 << "% of starting phi)" << endl;
			os.precision(n_prec);
			os.unsetf(ios_base::floatfield); // reset all flags to default

			if (upgrade_run.obs_valid() && (!best_run_updated_flag ||
				ModelRun::cmp_lt(upgrade_run, best_upgrade_run, *regul_scheme_ptr)))
			{
				best_run_updated_flag = true;
				best_upgrade_run = upgrade_run;
				new_frozen_pars = read_frozen_pars(fin_frz, i);
				best_lambda = i_lambda;
			}
		}
		else
		{
			streamsize n_prec = os.precision(2);
			os << "     Marquardt Lambda = ";
			os << setiosflags(ios::fixed) << setw(8) << i_lambda;
			os << ";   length = " << "NA";
			os.precision(n_prec);
			os.unsetf(ios_base::floatfield); // reset all flags to default
			os << ";    run failed" << endl;
		}
	}

	if (fosm_real_info.second.size() > 0)
	{
		pair<ObservationEnsemble, map<string, double>> fosm_obs_info = la.process_fosm_reals(&run_manager, fosm_real_info,
			termination_ctl.get_iteration_number(), base_run.get_phi(*regul_scheme_ptr));
	}

	// Print frozen parameter information for parameters frozen in SVD transformation
	const Parameters &frz_ctl_pars_svd = best_upgrade_run.get_frozen_ctl_pars();
	if (frz_ctl_pars_svd.size() > 0)
	{
		vector<string> keys = frz_ctl_pars_svd.get_keys();
		std::sort(keys.begin(), keys.end());
		os << endl;
		os << "    Parameters previously frozen in SVD transformation:" << endl;
		for (auto &ikey : keys)
		{
			auto iter = frz_ctl_pars_svd.find(ikey);
			if (iter != frz_ctl_pars_svd.end())
			{
				os << "      " << iter->first << " frozen at " << iter->second << endl;
			}
		}
	}

	if (new_frozen_pars.size() > 0)
	{
		vector<string> keys = new_frozen_pars.get_keys();
		std::sort(keys.begin(), keys.end());
		os << endl;
		os << "    Parameters frozen during upgrades:" << endl;
		for (auto &ikey : keys)
		{
			auto iter = new_frozen_pars.find(ikey);
			if (iter != new_frozen_pars.end())
			{
				os << "      " << iter->first << " frozen at " << iter->second << endl;
			}
		}
	}
	if (new_frozen_pars.size() > 0) best_upgrade_run.add_frozen_ctl_parameters(new_frozen_pars);
	// clean up run_manager memory
	run_manager.free_memory();

	if (!one_success)
	{
		throw runtime_error("all upgrade runs failed");
	}

	// Check if best_lambda is at the edge of lambda_vec

	// regrab lambda_vec
	vector<double> lambda_vec = pest_scenario.get_pestpp_options().get_base_lambda_vec();	
	std::sort(lambda_vec.begin(), lambda_vec.end());
	auto iter = std::unique(lambda_vec.begin(), lambda_vec.end());
	lambda_vec.resize(std::distance(lambda_vec.begin(), iter));

	auto last_lambda = lambda_vec.back();
	auto first_lambda = lambda_vec.front();
	file_manager.rec_ofstream() << "Checking to see if best lambda is at the edge" << std::endl;
	double lambda_spacing_factor = 10.0; // doing powers of 10 for now
	bool extended = false;

	if (best_lambda == last_lambda)
	{
		// Add a new larger lambda
		double new_lambda = last_lambda * lambda_spacing_factor;
		if (std::find(lambda_vec.begin(), lambda_vec.end(), new_lambda) == lambda_vec.end())
		{
			lambda_vec.push_back(new_lambda);
			extended = true;
			file_manager.rec_ofstream() << "*** Extending lambda_vec: added larger lambda " << new_lambda << std::endl;
		}
	}
	else if (best_lambda == first_lambda)
	{
		// Add a new smaller lambda
		double new_lambda = first_lambda / lambda_spacing_factor;
		if (std::find(lambda_vec.begin(), lambda_vec.end(), new_lambda) == lambda_vec.end())
		{
			lambda_vec.push_back(new_lambda);
			extended = true;
			file_manager.rec_ofstream() << "*** Extending lambda_vec: added smaller lambda " << new_lambda << std::endl;
		}
	}
	std::sort(lambda_vec.begin(), lambda_vec.end());
	pest_scenario.get_pestpp_options_ptr()->set_base_lambda_vec(lambda_vec);
	return best_upgrade_run;
}



SVDASolver::~SVDASolver(void)
{
}

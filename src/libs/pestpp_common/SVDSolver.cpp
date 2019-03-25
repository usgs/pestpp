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
#include <algorithm>
#include <sstream>
#include "SVDSolver.h"
#include "RunManagerAbstract.h"
#include "QSqrtMatrix.h"
#include "eigen_tools.h"
#include "ObjectiveFunc.h"
#include "utilities.h"
#include "FileManager.h"
#include "TerminationController.h"
#include "ParamTransformSeq.h"
#include "Transformation.h"
#include "PriorInformation.h"
#include "Regularization.h"
#include "SVD_PROPACK.h"
#include "OutputFileWriter.h"
#include "debug.h"
#include "covariance.h"

using namespace std;
using namespace pest_utils;
using namespace Eigen;

const string SVDSolver::svd_solver_type_name = "svd_base_par";


void MuPoint::set(double _mu, const PhiComponets &_phi_comp)
{
	mu = max(numeric_limits<double>::min(), _mu);
	mu = min(numeric_limits<double>::max(), mu);
	phi_comp = _phi_comp;

}

double MuPoint::f() const
{
	return phi_comp.meas - target_phi_meas;
}

double MuPoint::error_frac()
{
	return abs((phi_comp.meas - target_phi_meas) / target_phi_meas);
}

double MuPoint::error_percent()
{
	return (phi_comp.meas - target_phi_meas) / target_phi_meas;
}

void MuPoint::print(ostream &os)
{
		//streamsize n = os.precision(numeric_limits<double>::digits10 + 1);
		streamsize n = os.precision(6);
		//os << "    recalculating regularization weight factor         :" << endl;
		os << "      updated regularization  weight factor            : " << mu << endl;
		os << "      FRACPHIM adjusted measurement objective function : " << phi_comp.meas << endl;
		os.precision(3);
		os << "      discrepancy                                      : " << abs(error_percent() * 100) << "%" << endl;
		os.precision(n);

}
bool MuPoint::operator< (const MuPoint &rhs) const
{
	return abs(f()) < abs(rhs.f());
}



SVDSolver::SVDSolver(Pest &_pest_scenario, FileManager &_file_manager, ObjectiveFunc *_obj_func,
	const ParamTransformSeq &_par_transform, Jacobian &_jacobian,
	OutputFileWriter &_output_file_writer, SVDSolver::MAT_INV _mat_inv,
	PerformanceLog *_performance_log, const string &_description, Covariance _parcov, bool _phiredswh_flag, bool _splitswh_flag, bool _save_next_jacobian)
	: ctl_info(&_pest_scenario.get_control_info()), svd_info(_pest_scenario.get_svd_info()), par_group_info_ptr(&_pest_scenario.get_base_group_info()),
	ctl_par_info_ptr(&_pest_scenario.get_ctl_parameter_info()), obs_info_ptr(&_pest_scenario.get_ctl_observation_info()), obj_func(_obj_func),
	file_manager(_file_manager), observations_ptr(&_pest_scenario.get_ctl_observations()), par_transform(_par_transform), der_forgive(_pest_scenario.get_pestpp_options().get_der_forgive()), phiredswh_flag(_phiredswh_flag),
	splitswh_flag(_splitswh_flag), save_next_jacobian(_save_next_jacobian), prior_info_ptr(_pest_scenario.get_prior_info_ptr()), jacobian(_jacobian),
	regul_scheme_ptr(_pest_scenario.get_regul_scheme_ptr()), output_file_writer(_output_file_writer), mat_inv(_mat_inv), description(_description), best_lambda(20.0),
	performance_log(_performance_log), base_lambda_vec(_pest_scenario.get_pestpp_options().get_base_lambda_vec()), lambda_scale_vec(_pest_scenario.get_pestpp_options().get_lambda_scale_vec()),
	terminate_local_iteration(false), reg_frac(_pest_scenario.get_pestpp_options().get_reg_frac()),
		parcov(_parcov),parcov_scale_fac(_pest_scenario.get_pestpp_options().get_parcov_scale_fac()),upgrade_augment(_pest_scenario.get_pestpp_options().get_upgrade_augment())
{
	if (_pest_scenario.get_pestpp_options().get_jac_scale())
	{
		mar_mat = MarquardtMatrix::IDENT;
	}
	else
		mar_mat = MarquardtMatrix::JTQJ;
	if (_pest_scenario.get_pestpp_options().get_upgrade_bounds() == "ROBUST")
		upgrade_bounds = UpgradeBounds::ROBUST;
	else
		upgrade_bounds = UpgradeBounds::CHEAP;
	svd_package = new SVD_EIGEN();
}

void SVDSolver::set_svd_package(PestppOptions::SVD_PACK _svd_pack)
{
	if (_svd_pack == PestppOptions::PROPACK){
		delete svd_package;
		svd_package = new SVD_PROPACK;
	}
	else if (_svd_pack == PestppOptions::EIGEN){
		delete svd_package;
		svd_package = new SVD_EIGEN;
	}
	else if (_svd_pack == PestppOptions::REDSVD)
	{
		delete svd_package;
		svd_package = new SVD_REDSVD;

	}

	svd_package->set_max_sing(svd_info.maxsing);
	svd_package->set_eign_thres(svd_info.eigthresh);
	svd_package->set_performance_log(performance_log);
}

SVDSolver::~SVDSolver(void)
{
	delete svd_package;
}


ModelRun SVDSolver::solve(RunManagerAbstract &run_manager, TerminationController &termination_ctl, int max_iter,
	ModelRun &cur_run, ModelRun &optimum_run, RestartController &restart_controller, bool calc_first_jacobian)
{
	ostream &os = file_manager.rec_ofstream();
	ostream &fout_restart = file_manager.get_ofstream("rst");
	ModelRun best_upgrade_run(cur_run);
	// Start Solution iterations
	bool save_nextjac = false;
	string matrix_inv = (mat_inv == MAT_INV::Q12J) ? "\"Q 1/2 J\"" : "\"Jt Q J\"";
	terminate_local_iteration = false;

	bool calc_jacobian = calc_first_jacobian;

	if (restart_controller.get_restart_option() == RestartController::RestartOption::RESUME_NEW_ITERATION)
	{
		restart_controller.get_restart_option() = RestartController::RestartOption::NONE;
	}

	for (int iter_num = 1; iter_num <= max_iter && !terminate_local_iteration; ++iter_num)
	{
		//only processed when debugging is turned on
		int global_iter_num = termination_ctl.get_iteration_number() + 1;
		int nruns_start_iter = run_manager.get_total_runs();
		ostringstream tmp_str;

		debug_msg("========================================================================");
		debug_msg("Iteration Number");
		debug_print(global_iter_num);
		cout << endl << endl;
		os << endl;
		output_file_writer.iteration_report(cout, global_iter_num, run_manager.get_total_runs(), get_description(), svd_package->description, matrix_inv);

		if (restart_controller.get_restart_option() == RestartController::RestartOption::NONE)
		{
			output_file_writer.iteration_report(os, global_iter_num, run_manager.get_total_runs(), get_description(), svd_package->description, matrix_inv);
		}

		if (restart_controller.get_restart_option() != RestartController::RestartOption::RESUME_UPGRADE_RUNS)
		{
			//
			//Save information necessary for restart
			RestartController::write_start_iteration(fout_restart, this->get_solver_type(), iter_num, global_iter_num);
			// save state of termination controller
			termination_ctl.save_state(fout_restart);
			//write current parameters so we have a backup for restarting
			RestartController::write_start_parameters_updated(fout_restart, file_manager.build_filename("parb", false));
			output_file_writer.write_par(file_manager.open_ofile_ext("parb"), best_upgrade_run.get_ctl_pars(), *(par_transform.get_offset_ptr()),
				*(par_transform.get_scale_ptr()));
			file_manager.close_file("parb");
			RestartController::write_finish_parameters_updated(fout_restart, file_manager.build_filename("parb", false));

			// write header for SVD file
			output_file_writer.write_svd_iteration(global_iter_num);

			performance_log->log_blank_lines();
			performance_log->add_indent(-10);

			tmp_str << "beginning iteration " << global_iter_num;
			performance_log->log_event(tmp_str.str(), 0, "start_iter");
			performance_log->add_indent();
			if (!calc_jacobian)
			{
				calc_jacobian = true;
			}
			else
			{
				bool restart_runs = (restart_controller.get_restart_option() == RestartController::RestartOption::RESUME_JACOBIAN_RUNS);
				iteration_jac(run_manager, termination_ctl, best_upgrade_run, false, restart_runs);
				if (restart_runs) restart_controller.get_restart_option() = RestartController::RestartOption::NONE;
			}

			// Update Regularization weights if REG_FRAC is used
			ModelRun prev_run(best_upgrade_run);
			if (regul_scheme_ptr->get_use_dynamic_reg() && reg_frac > 0.0)
			{
				//dynamic_weight_adj_percent(best_upgrade_run, 0.10);
				dynamic_weight_adj_percent(best_upgrade_run, reg_frac);
			}
			os << endl;
			// write out report for starting phi
			PhiData phi_data = obj_func->phi_report(best_upgrade_run.get_obs(), best_upgrade_run.get_ctl_pars(), *regul_scheme_ptr);
			output_file_writer.phi_report(os, termination_ctl.get_iteration_number() + 1, run_manager.get_total_runs(), phi_data, regul_scheme_ptr->get_weight());
			//print out number of zero terms in jacobian
			long jac_num_nonzero = jacobian.get_nonzero();
			long jac_num_total = jacobian.get_size();
			long jac_num_zero = jac_num_total - jac_num_nonzero;
			streamsize n_prec = os.precision(2);
			os << "  Number of terms in the jacobian equal to zero: " << jac_num_zero << " / " << jac_num_total
				<< " (" << double(jac_num_zero) / double(jac_num_total) * 100 << "%)" << endl << endl;
			os.precision(n_prec);
		}

		{
			const set<string> &failed_parameter_runs = jacobian.get_failed_parameter_names();
			if (!failed_parameter_runs.empty())
			{
				os << endl;
				jacobian.report_errors(os);
				os << endl;
				cout << endl;
				jacobian.report_errors(cout);
				cout << endl;
				os.flush();
				cout.flush();
				if (!der_forgive) exit(0);
			}
		}

		cout << "  computing upgrade vectors... " << endl;
		cout.flush();

		// This must be located before the call to iteration_upgrd.  As iteration_upgrd will update the base run
		//observations when performing a restart
		ModelRun prev_run(best_upgrade_run);
		bool upgrade_start = (restart_controller.get_restart_option() == RestartController::RestartOption::RESUME_UPGRADE_RUNS);
		best_upgrade_run = iteration_upgrd(run_manager, termination_ctl, prev_run, upgrade_start);
		// reload best parameters and set flag to switch to central derivatives next iteration
		double prev_phi = prev_run.get_phi(*regul_scheme_ptr);
		double best_new_phi = best_upgrade_run.get_phi(*regul_scheme_ptr);
		double phi_ratio = best_new_phi / prev_phi;

		cout << endl << "  ...Lambda testing complete for iteration " << termination_ctl.get_iteration_number() + 1 << endl;
		cout << "    starting phi = " << prev_phi << ";  ending phi = " << best_new_phi <<
			"  (" << phi_ratio * 100 << "% of starting phi)" << endl;

		if (prev_phi != 0 && par_group_info_ptr->have_switch_derivative() && !phiredswh_flag &&
			termination_ctl.get_iteration_number() + 1 > ctl_info->noptswitch &&
			(prev_phi - best_new_phi) / prev_phi < ctl_info->phiredswh)
		{
			phiredswh_flag = true;
			os << endl << "      Switching to central derivatives:" << endl;
			cout << endl << "    Switching to central derivatives:" << endl;
		}

		else if (!splitswh_flag && phiredswh_flag && par_group_info_ptr->have_switch_derivative() && prev_phi != 0 &&
			ctl_info->splitswh > 0 && phi_ratio >= ctl_info->splitswh)
		{
			splitswh_flag = true;
			os << endl << "    Switching to split threshold derivatives" << endl << endl;
			cout << endl << "  Switching to split threshold derivatives" << endl << endl;
		}

		restart_controller.get_restart_option() = RestartController::RestartOption::NONE;

		int nruns_end_iter = run_manager.get_total_runs();
		os << endl;
		os << "  Model calls in iteration " << global_iter_num << ": " << nruns_end_iter - nruns_start_iter << endl;
		os << "  Total model calls at end of iteration " << global_iter_num << ": " << nruns_end_iter << endl;
		if (terminate_local_iteration)
		{
			break;
		}
		tmp_str.str("");
		tmp_str.clear();
		tmp_str << "completed iteration " << global_iter_num;
		performance_log->log_event(tmp_str.str(), 0, "end_iter");
		tmp_str.str("");
		tmp_str.clear();
		tmp_str << "time to complete iteration " << global_iter_num;
		performance_log->log_summary(tmp_str.str(), "end_iter", "start_iter");
		// write files that get wrtten at the end of each iteration
		stringstream filename;
		string complete_filename;

		// rei file for this iteration
		filename << "rei" << global_iter_num;
		output_file_writer.write_rei(file_manager.open_ofile_ext(filename.str()), global_iter_num,
			*(best_upgrade_run.get_obj_func_ptr()->get_obs_ptr()),
			best_upgrade_run.get_obs(), *(best_upgrade_run.get_obj_func_ptr()),
			best_upgrade_run.get_ctl_pars());
		file_manager.close_file(filename.str());
		// par file for this iteration
		output_file_writer.write_par(file_manager.open_ofile_ext("par"), best_upgrade_run.get_ctl_pars(), *(par_transform.get_offset_ptr()),
			*(par_transform.get_scale_ptr()));
		file_manager.close_file("par");

		filename.str(""); // reset the stringstream
		filename << "par" << global_iter_num;
		output_file_writer.write_par(file_manager.open_ofile_ext(filename.str()), best_upgrade_run.get_ctl_pars(), *(par_transform.get_offset_ptr()),
			*(par_transform.get_scale_ptr()));
		file_manager.close_file(filename.str());
		if (save_nextjac) {
			if (description.find("base") != string::npos)
				output_file_writer.write_jco(true, "jco",jacobian);
			else
				output_file_writer.write_jco(false, "jco",jacobian);

			//jacobian.save();
		}
		if (!optimum_run.obs_valid() || ModelRun::cmp_lt(best_upgrade_run, optimum_run, *regul_scheme_ptr))
		{
			optimum_run.set_ctl_parameters(best_upgrade_run.get_ctl_pars());
			optimum_run.set_observations(best_upgrade_run.get_obs());
			// save new optimum parameters to .par file
			output_file_writer.write_par(file_manager.open_ofile_ext("bpa"), optimum_run.get_ctl_pars(), *(par_transform.get_offset_ptr()),
				*(par_transform.get_scale_ptr()));
			file_manager.close_file("bpa");
			// save new optimum residuals to .rei file
			output_file_writer.write_rei(file_manager.open_ofile_ext("rei"), global_iter_num,
				*(optimum_run.get_obj_func_ptr()->get_obs_ptr()),
				optimum_run.get_obs(), *(optimum_run.get_obj_func_ptr()),
				optimum_run.get_ctl_pars());
			file_manager.close_file("rei");
			if (description.find("base") != string::npos)
				output_file_writer.write_jco(true, "jco", jacobian);
			else
				output_file_writer.write_jco(false, "jco", jacobian);
			//jacobian.save();
			// jacobian calculated next iteration will be at the current parameters and
			// will be more accurate than the one caluculated at the begining of this iteration
			save_nextjac = true;
		}
		os << endl;
		iteration_update_and_report(os, prev_run, best_upgrade_run, termination_ctl, run_manager);

		if (termination_ctl.check_last_iteration()){
			break;
		}
	}
	return best_upgrade_run;
}

VectorXd SVDSolver::calc_residual_corrections(const Jacobian &jacobian, const Parameters &del_numeric_pars,
	const vector<string> obs_name_vec)
{
	VectorXd del_residuals;
	if (del_numeric_pars.size() > 0)
	{
		vector<string>frz_par_name_vec = del_numeric_pars.get_keys();
		//remove the parameters for which the jaocbian could not be computed
		const set<string> &failed_jac_par_names = jacobian.get_failed_parameter_names();
		auto end_iter = remove_if(frz_par_name_vec.begin(), frz_par_name_vec.end(),
			[&failed_jac_par_names](string &str)->bool{return failed_jac_par_names.find(str) != failed_jac_par_names.end(); });
		frz_par_name_vec.resize(std::distance(frz_par_name_vec.begin(), end_iter));

		VectorXd frz_del_par_vec = del_numeric_pars.get_data_eigen_vec(frz_par_name_vec);

		MatrixXd jac_frz = jacobian.get_matrix(obs_name_vec, frz_par_name_vec);
		del_residuals = (jac_frz)*  frz_del_par_vec;
	}
	else
	{
		del_residuals = VectorXd::Zero(obs_name_vec.size());
	}
	return del_residuals;
}

void SVDSolver::calc_lambda_upgrade_vec_JtQJ(const Jacobian &jacobian, const QSqrtMatrix &Q_sqrt, const DynamicRegularization &regul,
	const Eigen::VectorXd &Residuals, const vector<string> &obs_name_vec,
	const Parameters &base_active_ctl_pars, const Parameters &prev_frozen_active_ctl_pars,
	double lambda, Parameters &active_ctl_upgrade_pars, Parameters &upgrade_active_ctl_del_pars,
	Parameters &grad_active_ctl_del_pars, MarquardtMatrix marquardt_type, bool scale_upgrade)
{
	Parameters base_numeric_pars = par_transform.active_ctl2numeric_cp(base_active_ctl_pars);
	//Create a set of Derivative Parameters which does not include the frozen Parameters
	Parameters pars_nf = base_active_ctl_pars;
	pars_nf.erase(prev_frozen_active_ctl_pars);
	//Transform these parameters to numeric parameters
	par_transform.active_ctl2numeric_ip(pars_nf);
	vector<string> numeric_par_names = pars_nf.get_keys();

	//Compute effect of frozen parameters on the residuals vector
	Parameters delta_freeze_pars = prev_frozen_active_ctl_pars;
	Parameters base_freeze_pars(base_active_ctl_pars, delta_freeze_pars.get_keys());
	par_transform.ctl2numeric_ip(delta_freeze_pars);
	par_transform.ctl2numeric_ip(base_freeze_pars);
	delta_freeze_pars -= base_freeze_pars;
	VectorXd del_residuals = calc_residual_corrections(jacobian, delta_freeze_pars, obs_name_vec);
	VectorXd corrected_residuals = Residuals + del_residuals;
	VectorXd Sigma;
	VectorXd Sigma_trunc;
	Eigen::SparseMatrix<double> U;
	Eigen::SparseMatrix<double> Vt;
	// the last boolean arguement is an instruction to compute the square weights
	Eigen::SparseMatrix<double> q_mat = Q_sqrt.get_sparse_matrix(obs_name_vec, regul, true);
	// removed this line when true added to end of the previous call to get_sparce_matrix
	//q_mat = (q_mat * q_mat).eval();
	Eigen::SparseMatrix<double> jac = jacobian.get_matrix(obs_name_vec, numeric_par_names);
	Eigen::SparseMatrix<double> ident;
	ident.resize(jac.cols(), jac.cols());
	ident.setIdentity();
	performance_log->log_event("forming JtQJ matrix");
	Eigen::SparseMatrix<double> JtQJ = jac.transpose() * q_mat * jac;

	if (parcov.nrow() > 0)
	//if (parcov_scale_fac > 0.0)
	{
		cout << parcov_scale_fac << endl;
		performance_log->log_event("JtQJ plus parcov.inv");
		JtQJ = JtQJ + (parcov_scale_fac * *parcov.get(numeric_par_names).inv().e_ptr());
	}

	Eigen::VectorXd upgrade_vec;
	if (marquardt_type == MarquardtMatrix::IDENT)
	{
		stringstream info_str;
		bool recalc_js = false;
		Eigen::SparseMatrix<double> S;

		if (true)
		{
			//Compute Scaling Matrix Sii
			performance_log->log_event("commencing to scale JtQJ matrix");
			svd_package->solve_ip(JtQJ, Sigma, U, Vt, Sigma_trunc, 0.0);
			VectorXd Sigma_inv_sqrt = Sigma.array().inverse().sqrt();
			S = Vt.transpose() * Sigma_inv_sqrt.asDiagonal() * U.transpose();
			VectorXd S_diag = S.diagonal();
			MatrixXd S_tmp = S_diag.asDiagonal();
			S = S_tmp.sparseView();
			stringstream info_str1;
			info_str1 << "S info: " << "rows = " << S.rows() << ": cols = " << S.cols() << ": size = " << S.size() << ": nonzeros = " << S.nonZeros();
			performance_log->log_event(info_str1.str());
			performance_log->log_event("JS");
			JS = jac * S;
			performance_log->log_event("JS.transpose() * q_mat * JS + lambda * S.transpose() * S");
			JtQJ = JS.transpose() * q_mat * JS + lambda * S.transpose() * S;
			info_str.str("");
			info_str << "S info: " << "rows = " << S.rows() << ": cols = " << S.cols() << ": size = " << S.size() << ": nonzeros = " << S.nonZeros();
			performance_log->log_event(info_str.str());
			recalc_js = true;
		}



		// Returns truncated Sigma, U and Vt arrays with small singular parameters trimed off
		performance_log->log_event("commencing SVD factorization");
		svd_package->solve_ip(JtQJ, Sigma, U, Vt, Sigma_trunc);
		performance_log->log_event("SVD factorization complete");

		if (!recalc_js)
		{
			VectorXd Sigma_inv_sqrt = Sigma.array().inverse().sqrt();
			S = Vt.transpose() * Sigma_inv_sqrt.asDiagonal() * U.transpose();
			VectorXd S_diag = S.diagonal();
			MatrixXd S_tmp = S_diag.asDiagonal();
			S = S_tmp.sparseView();
		}

		output_file_writer.write_svd(Sigma, Vt, lambda, prev_frozen_active_ctl_pars, Sigma_trunc);

		VectorXd Sigma_inv = Sigma.array().inverse();
		performance_log->log_event("commencing linear algebra multiplication to compute ugrade");

		info_str << "Vt info: " << "rows = " << Vt.rows() << ": cols = " << Vt.cols() << ": size = " << Vt.size() << ": nonzeros = " << Vt.nonZeros();
		performance_log->log_event(info_str.str());
		info_str.str("");
		info_str << "U info: " << "rows = " << U.rows() << ": cols = " << U.cols() << ": size = " << U.size() << ": nonzeros = " << U.nonZeros();
		performance_log->log_event(info_str.str());
		info_str.str("");
		info_str << "jac info: " << "rows = " << jac.rows() << ": cols = " << jac.cols() << ": size = " << jac.size() << ": nonzeros = " << jac.nonZeros();
		performance_log->log_event(info_str.str());
		upgrade_vec = S * (Vt.transpose() * (Sigma_inv.asDiagonal() * (U.transpose() * ((jac * S).transpose()* (q_mat  * (corrected_residuals))))));
	}
	else
	{
		performance_log->log_event("commencing SVD factorization");
		svd_package->solve_ip(JtQJ, Sigma, U, Vt, Sigma_trunc);
		performance_log->log_event("SVD factorization complete");
		//Only add lambda to singular values above the threshhold
		Sigma = Sigma.array() + (Sigma.cwiseProduct(Sigma).array() * lambda).sqrt();
		output_file_writer.write_svd(Sigma, Vt, lambda, prev_frozen_active_ctl_pars, Sigma_trunc);
		VectorXd Sigma_inv = Sigma.array().inverse();

		performance_log->log_event("commencing linear algebra multiplication to compute ugrade");
		stringstream info_str;
		info_str << "Vt info: " << "rows = " << Vt.rows() << ": cols = " << Vt.cols() << ": size = " << Vt.size() << ": nonzeros = " << Vt.nonZeros();
		performance_log->log_event(info_str.str());
		info_str.str("");
		info_str << "U info: " << "rows = " << U.rows() << ": cols = " << U.cols() << ": size = " << U.size() << ": nonzeros = " << U.nonZeros();
		performance_log->log_event(info_str.str());
		info_str.str("");
		info_str << "jac info: " << "rows = " << jac.rows() << ": cols = " << jac.cols() << ": size = " << jac.size() << ": nonzeros = " << jac.nonZeros();
		performance_log->log_event(info_str.str());
		upgrade_vec = Vt.transpose() * (Sigma_inv.asDiagonal() * (U.transpose() * (jac.transpose() * (q_mat  * corrected_residuals))));
	}

	// scale the upgrade vector using the technique described in the PEST manual
	if (scale_upgrade)
	{
		double beta = 1.0;
		Eigen::VectorXd gama = jac * upgrade_vec;
		Eigen::SparseMatrix<double> Q_mat_tmp = Q_sqrt.get_sparse_matrix(obs_name_vec, regul);
		Eigen::SparseMatrix<double> Q_diag = get_diag_matrix(Q_mat_tmp);
		Q_diag = (Q_diag * Q_diag).eval();
		double top = corrected_residuals.transpose() * Q_diag * gama;
		double bot = gama.transpose() * Q_diag * gama;
		if (bot != 0)
		{
			beta = top / bot;
		}
		upgrade_vec *= beta;
	}


	Eigen::VectorXd grad_vec;
	grad_vec = -2.0 * (jac.transpose() * (q_mat * Residuals));
	performance_log->log_event("linear algebra multiplication to compute ugrade complete");

	//tranfere newly computed componets of the ugrade vector to upgrade.svd_uvec
	upgrade_active_ctl_del_pars.clear();
	grad_active_ctl_del_pars.clear();

	string *name_ptr;
	auto it_nf_end = pars_nf.end();
	for (size_t i = 0; i < numeric_par_names.size(); ++i)
	{
		name_ptr = &(numeric_par_names[i]);
		upgrade_active_ctl_del_pars[*name_ptr] = upgrade_vec(i);
		grad_active_ctl_del_pars[*name_ptr] = grad_vec(i);
		auto it_nf = pars_nf.find(*name_ptr);
		if (it_nf != it_nf_end)
		{
			it_nf->second += upgrade_vec(i);
		}
	}
	// Transform upgrade_pars back to derivative parameters
	active_ctl_upgrade_pars = par_transform.numeric2active_ctl_cp(pars_nf);
	Parameters tmp_pars(base_numeric_pars);
	par_transform.del_numeric_2_del_active_ctl_ip(upgrade_active_ctl_del_pars, tmp_pars);
	tmp_pars = base_numeric_pars;
	par_transform.del_numeric_2_del_active_ctl_ip(grad_active_ctl_del_pars, tmp_pars);

	//tranfere previously frozen componets of the ugrade vector to upgrade.svd_uvec
	for (auto &ipar : prev_frozen_active_ctl_pars)
	{
		active_ctl_upgrade_pars[ipar.first] = ipar.second;
	}
}

void SVDSolver::calc_lambda_upgrade_vecQ12J(const Jacobian &jacobian, const QSqrtMatrix &Q_sqrt, const DynamicRegularization &regul,
	const Eigen::VectorXd &Residuals, const vector<string> &obs_name_vec,
	const Parameters &base_active_ctl_pars, const Parameters &prev_frozen_active_ctl_pars,
	double lambda, Parameters &active_ctl_upgrade_pars, Parameters &upgrade_active_ctl_del_pars,
	Parameters &grad_active_ctl_del_pars, MarquardtMatrix marquardt_type, bool scale_upgrade)
{
	Parameters base_numeric_pars = par_transform.active_ctl2numeric_cp(base_active_ctl_pars);
	//Create a set of Ctl Parameters which does not include the frozen Parameters
	Parameters pars_nf = base_active_ctl_pars;
	pars_nf.erase(prev_frozen_active_ctl_pars);
	//Transform these parameters to numeric parameters
	par_transform.active_ctl2numeric_ip(pars_nf);
	vector<string> numeric_par_names = pars_nf.get_keys();

	//Compute effect of frozen parameters on the residuals vector
	Parameters delta_freeze_pars = prev_frozen_active_ctl_pars;
	Parameters base_freeze_pars(base_active_ctl_pars, delta_freeze_pars.get_keys());
	par_transform.active_ctl2numeric_ip(delta_freeze_pars);
	par_transform.active_ctl2numeric_ip(base_freeze_pars);
	delta_freeze_pars -= base_freeze_pars;
	VectorXd del_residuals = calc_residual_corrections(jacobian, delta_freeze_pars, obs_name_vec);
	VectorXd corrected_residuals = Residuals + del_residuals;

	VectorXd Sigma;
	VectorXd Sigma_trunc;
	Eigen::SparseMatrix<double> U;
	Eigen::SparseMatrix<double> Vt;
	Eigen::SparseMatrix<double> q_sqrt = Q_sqrt.get_sparse_matrix(obs_name_vec, regul);
	Eigen::SparseMatrix<double> jac = jacobian.get_matrix(obs_name_vec, numeric_par_names);
	Eigen::SparseMatrix<double> SqrtQ_J = q_sqrt * jac;
	// Returns truncated Sigma, U and Vt arrays with small singular parameters trimed off
	performance_log->log_event("commencing SVD factorization");
	svd_package->solve_ip(SqrtQ_J, Sigma, U, Vt, Sigma_trunc);
	performance_log->log_event("SVD factorization complete");
	//Only add lambda to singular values above the threshhold
	if (marquardt_type == MarquardtMatrix::IDENT)
	{
		Sigma = Sigma.array() + lambda;
	}
	else
	{
		//this needs checking
		Sigma = Sigma.array() + (Sigma.cwiseProduct(Sigma).array() * lambda).sqrt();
	}
	output_file_writer.write_svd(Sigma, Vt, lambda, prev_frozen_active_ctl_pars, Sigma_trunc);
	VectorXd Sigma_inv = Sigma.array().inverse();

	performance_log->log_event("commencing linear algebra multiplication to compute ugrade");
	stringstream info_str;
	info_str << "Vt info: " << "rows = " << Vt.rows() << ": cols = " << Vt.cols() << ": size = " << Vt.size() << ": nonzeros = " << Vt.nonZeros();
	performance_log->log_event(info_str.str());
	info_str.str("");
	info_str << "U info: " << "rows = " << U.rows() << ": cols = " << U.cols() << ": size = " << U.size() << ": nonzeros = " << U.nonZeros();
	performance_log->log_event(info_str.str());
	Eigen::VectorXd upgrade_vec;
	upgrade_vec = Vt.transpose() * (Sigma_inv.asDiagonal() * (U.transpose() * (q_sqrt  * corrected_residuals)));



	// scale the upgrade vector using the technique described in the PEST manual
	if (scale_upgrade)
	{
		double beta = 1.0;
		Eigen::VectorXd gama = jac * upgrade_vec;
		Eigen::SparseMatrix<double> Q_diag = get_diag_matrix(q_sqrt);
		Q_diag = (Q_diag * Q_diag).eval();
		double top = corrected_residuals.transpose() * Q_diag * gama;
		double bot = gama.transpose() * Q_diag * gama;
		if (bot != 0)
		{
			beta = top / bot;
		}
		upgrade_vec *= beta;
	}


	Eigen::VectorXd grad_vec;
	grad_vec = -2.0 * (jac.transpose() * (q_sqrt * (q_sqrt * Residuals)));
	performance_log->log_event("linear algebra multiplication to compute ugrade complete");

	//tranfere newly computed componets of the ugrade vector to upgrade.svd_uvec
	upgrade_active_ctl_del_pars.clear();
	grad_active_ctl_del_pars.clear();

	string *name_ptr;
	auto it_nf_end = pars_nf.end();
	for (size_t i = 0; i < numeric_par_names.size(); ++i)
	{
		name_ptr = &(numeric_par_names[i]);
		upgrade_active_ctl_del_pars[*name_ptr] = upgrade_vec(i);
		grad_active_ctl_del_pars[*name_ptr] = grad_vec(i);
		auto it_nf = pars_nf.find(*name_ptr);
		if (it_nf != it_nf_end)
		{
			it_nf->second += upgrade_vec(i);
		}
	}
	// Transform upgrade_pars back to ctl parameters
	active_ctl_upgrade_pars = par_transform.numeric2active_ctl_cp(pars_nf);
	Parameters tmp_pars(base_numeric_pars);
	par_transform.del_numeric_2_del_active_ctl_ip(upgrade_active_ctl_del_pars, tmp_pars);
	tmp_pars = base_numeric_pars;
	par_transform.del_numeric_2_del_active_ctl_ip(grad_active_ctl_del_pars, tmp_pars);

	//tranfere previously frozen componets of the ugrade vector to upgrade.svd_uvec
	for (auto &ipar : prev_frozen_active_ctl_pars)
	{
		active_ctl_upgrade_pars[ipar.first] = ipar.second;
	}
}



void SVDSolver::calc_upgrade_vec(double i_lambda, Parameters &prev_frozen_active_ctl_pars, QSqrtMatrix &Q_sqrt,
	const DynamicRegularization &regul, VectorXd &residuals_vec, vector<string> &obs_names_vec,
	const Parameters &base_run_active_ctl_pars, Parameters &upgrade_active_ctl_pars,
	MarquardtMatrix marquardt_type, LimitType &limit_type, bool scale_upgrade)
{
	Parameters upgrade_ctl_del_pars;
	Parameters grad_ctl_del_pars;
	int num_upgrade_out_grad_in;
	Parameters new_frozen_active_ctl_pars;

	upgrade_active_ctl_pars.clear();
	// define a function type for upgrade methods
	typedef void(SVDSolver::*UPGRADE_FUNCTION) (const Jacobian &jacobian, const QSqrtMatrix &Q_sqrt,
		const DynamicRegularization &regul,
		const Eigen::VectorXd &Residuals, const vector<string> &obs_name_vec,
		const Parameters &base_ctl_pars, const Parameters &prev_frozen_ctl_pars,
		double lambda, Parameters &ctl_upgrade_pars, Parameters &upgrade_ctl_del_pars,
		Parameters &grad_ctl_del_pars, MarquardtMatrix marquardt_type, bool scale_upgrade);

	UPGRADE_FUNCTION calc_lambda_upgrade = &SVDSolver::calc_lambda_upgrade_vec_JtQJ;

	if (mat_inv == MAT_INV::Q12J)
	{
		calc_lambda_upgrade = &SVDSolver::calc_lambda_upgrade_vecQ12J;
	}

	// need to remove parameters frozen due to failed jacobian runs when calling calc_lambda_upgrade_vec
	//Freeze Parameters at the boundary whose ugrade vector and gradient both head out of bounds
	performance_log->log_event("commencing calculation of upgrade vector");
	(*this.*calc_lambda_upgrade)(jacobian, Q_sqrt, regul, residuals_vec, obs_names_vec,
		base_run_active_ctl_pars, prev_frozen_active_ctl_pars, i_lambda, upgrade_active_ctl_pars, upgrade_ctl_del_pars,
		grad_ctl_del_pars, marquardt_type, scale_upgrade);
	if (upgrade_active_ctl_pars.get_notnormal_keys().size() > 0)
	{
		stringstream ss;
		for (auto &n : upgrade_active_ctl_pars.get_notnormal_keys())
			ss << n << ',';
		string message = "not normal floating point in upgrade pars: " + ss.str();
		file_manager.rec_ofstream() << message << endl;
		cout << message;
		throw runtime_error(message);

	}
	performance_log->log_event("commencing check of parameter bounds");
	num_upgrade_out_grad_in = check_bnd_par(new_frozen_active_ctl_pars, base_run_active_ctl_pars, upgrade_ctl_del_pars, grad_ctl_del_pars);
	prev_frozen_active_ctl_pars.insert(new_frozen_active_ctl_pars.begin(), new_frozen_active_ctl_pars.end());

	//Recompute the ugrade vector without the newly frozen parameters and freeze those at the boundary whose upgrade still goes heads out of bounds
	if ((upgrade_bounds == UpgradeBounds::ROBUST) && (num_upgrade_out_grad_in > 0))
	{
		new_frozen_active_ctl_pars.clear();
		performance_log->log_event("commencing recalculation of upgrade vector freezing parameters whose upgrade and gradient point out of bounds");
		(*this.*calc_lambda_upgrade)(jacobian, Q_sqrt, regul, residuals_vec, obs_names_vec,
			base_run_active_ctl_pars, prev_frozen_active_ctl_pars, i_lambda, upgrade_active_ctl_pars, upgrade_ctl_del_pars,
			grad_ctl_del_pars, marquardt_type, scale_upgrade);
		performance_log->log_event("commencing check of parameter bounds with new parameters");
		check_bnd_par(new_frozen_active_ctl_pars, prev_frozen_active_ctl_pars, upgrade_active_ctl_pars);
		prev_frozen_active_ctl_pars.insert(new_frozen_active_ctl_pars.begin(), new_frozen_active_ctl_pars.end());
		new_frozen_active_ctl_pars.clear();
	}
	//If there are newly frozen parameters recompute the upgrade vector
	if ((upgrade_bounds == UpgradeBounds::ROBUST) && (new_frozen_active_ctl_pars.size() > 0))
	{
		performance_log->log_event("commencing recalculation of upgrade vector freezing parameters whose upgrade heads out of bounds");
		(*this.*calc_lambda_upgrade)(jacobian, Q_sqrt, regul, residuals_vec, obs_names_vec,
			base_run_active_ctl_pars, prev_frozen_active_ctl_pars, i_lambda, upgrade_active_ctl_pars, upgrade_ctl_del_pars,
			grad_ctl_del_pars, marquardt_type, scale_upgrade);
	}
	//Freeze any new parameters that want to go out of bounds
	performance_log->log_event("limiting out of bounds pars");
	limit_parameters_ip(base_run_active_ctl_pars, upgrade_active_ctl_pars,
		limit_type, prev_frozen_active_ctl_pars);
	performance_log->log_event("checking for denormal floating point values");
	if (upgrade_active_ctl_pars.get_notnormal_keys().size() > 0)
	{
		stringstream ss;
		for (auto &n : upgrade_active_ctl_pars.get_notnormal_keys())
			ss << n << ',';
		string message = "not normal floating point in upgrade pars after applying limits: " + ss.str();
		file_manager.rec_ofstream() << message << endl;
		cout << message;
		throw runtime_error(message);

	}
		//	limit_type, freeze_active_ctl_pars, true);
	//new_frozen_active_ctl_pars.clear();
	//new_frozen_active_ctl_pars = limit_parameters_freeze_all_ip(base_run_active_ctl_pars, base_run_active_ctl_pars, prev_frozen_active_ctl_pars);
	//prev_frozen_active_ctl_pars.insert(new_frozen_active_ctl_pars.begin(), new_frozen_active_ctl_pars.end());
}

void SVDSolver::calc_upgrade_vec_freeze(double i_lambda, Parameters &prev_frozen_active_ctl_pars, QSqrtMatrix &Q_sqrt,
		const DynamicRegularization &regul, VectorXd &residuals_vec, vector<string> &obs_names_vec,
		const Parameters &base_run_active_ctl_pars, Parameters &upgrade_active_ctl_pars,
		MarquardtMatrix marquardt_type, bool scale_upgrade)
{
	Parameters upgrade_ctl_del_pars;
	Parameters grad_ctl_del_pars;
	int num_upgrade_out_grad_in;
	Parameters new_frozen_active_ctl_pars;

	upgrade_active_ctl_pars.clear();
	// define a function type for upgrade methods
	typedef void(SVDSolver::*UPGRADE_FUNCTION) (const Jacobian &jacobian, const QSqrtMatrix &Q_sqrt,
		const DynamicRegularization &regul,
		const Eigen::VectorXd &Residuals, const vector<string> &obs_name_vec,
		const Parameters &base_ctl_pars, const Parameters &prev_frozen_ctl_pars,
		double lambda, Parameters &ctl_upgrade_pars, Parameters &upgrade_ctl_del_pars,
		Parameters &grad_ctl_del_pars, MarquardtMatrix marquardt_type, bool scale_upgrade);

	UPGRADE_FUNCTION calc_lambda_upgrade = &SVDSolver::calc_lambda_upgrade_vec_JtQJ;

	if (mat_inv == MAT_INV::Q12J)
	{
		calc_lambda_upgrade = &SVDSolver::calc_lambda_upgrade_vecQ12J;
	}

	// need to remove parameters frozen due to failed jacobian runs when calling calc_lambda_upgrade_vec
	//Freeze Parameters at the boundary whose ugrade vector and gradient both head out of bounds
	performance_log->log_event("commencing calculation of upgrade vector");
	(*this.*calc_lambda_upgrade)(jacobian, Q_sqrt, regul, residuals_vec, obs_names_vec,
		base_run_active_ctl_pars, prev_frozen_active_ctl_pars, i_lambda, upgrade_active_ctl_pars, upgrade_ctl_del_pars,
		grad_ctl_del_pars, marquardt_type, scale_upgrade);
	performance_log->log_event("commencing check of parameter bounds");
	num_upgrade_out_grad_in = check_bnd_par(new_frozen_active_ctl_pars, base_run_active_ctl_pars, upgrade_ctl_del_pars, grad_ctl_del_pars);
	prev_frozen_active_ctl_pars.insert(new_frozen_active_ctl_pars.begin(), new_frozen_active_ctl_pars.end());
	//Recompute the ugrade vector without the newly frozen parameters and freeze those at the boundary whose upgrade still goes heads out of bounds
	if ((upgrade_bounds == UpgradeBounds::ROBUST) && (num_upgrade_out_grad_in > 0))
	{
		new_frozen_active_ctl_pars.clear();
		performance_log->log_event("commencing recalculation of upgrade vector freezing parameters whose upgrade and gradient point out of bounds");
		(*this.*calc_lambda_upgrade)(jacobian, Q_sqrt, regul, residuals_vec, obs_names_vec,
			base_run_active_ctl_pars, prev_frozen_active_ctl_pars, i_lambda, upgrade_active_ctl_pars, upgrade_ctl_del_pars,
			grad_ctl_del_pars, marquardt_type, scale_upgrade);
		performance_log->log_event("commencing check of parameter bounds with new parameters");
		check_bnd_par(new_frozen_active_ctl_pars, prev_frozen_active_ctl_pars, upgrade_active_ctl_pars);
		prev_frozen_active_ctl_pars.insert(new_frozen_active_ctl_pars.begin(), new_frozen_active_ctl_pars.end());
		new_frozen_active_ctl_pars.clear();
	}
	//If there are newly frozen parameters recompute the upgrade vector
	if ((upgrade_bounds == UpgradeBounds::ROBUST) && (new_frozen_active_ctl_pars.size() > 0))
	{
		performance_log->log_event("commencing recalculation of upgrade vector freezing parameters whose upgrade heads out of bounds");
		(*this.*calc_lambda_upgrade)(jacobian, Q_sqrt, regul, residuals_vec, obs_names_vec,
			base_run_active_ctl_pars, prev_frozen_active_ctl_pars, i_lambda, upgrade_active_ctl_pars, upgrade_ctl_del_pars,
			grad_ctl_del_pars, marquardt_type, scale_upgrade);
	}
	//Freeze any new parameters that want to go out of bounds
	new_frozen_active_ctl_pars.clear();
	new_frozen_active_ctl_pars = limit_parameters_freeze_all_ip(base_run_active_ctl_pars, upgrade_active_ctl_pars, prev_frozen_active_ctl_pars);
	prev_frozen_active_ctl_pars.insert(new_frozen_active_ctl_pars.begin(), new_frozen_active_ctl_pars.end());
}


	ModelRun SVDSolver::compute_jacobian(RunManagerAbstract &run_manager, TerminationController &termination_ctl, ModelRun &cur_run, bool restart_runs)
	{
		ostream &os = file_manager.rec_ofstream();
		ostream &fout_restart = file_manager.get_ofstream("rst");
		ModelRun best_upgrade_run(cur_run);
		// Start Solution iterations
		terminate_local_iteration = false;

		RestartController::write_start_iteration(fout_restart, this->get_solver_type(), -9999, -9999);

		//write current parameters so we have a backup for restarting
		RestartController::write_start_parameters_updated(fout_restart, file_manager.build_filename("parb", false));
		output_file_writer.write_par(file_manager.open_ofile_ext("parb"), best_upgrade_run.get_ctl_pars(), *(par_transform.get_offset_ptr()),
			*(par_transform.get_scale_ptr()));
		file_manager.close_file("parb");
		RestartController::write_finish_parameters_updated(fout_restart, file_manager.build_filename("parb", false));

		cout << "COMPUTING JACOBIAN:" << endl << endl;
		os << "COMPUTING JACOBIAN:" << endl << endl;
		cout << "  Iteration type: " << get_description() << endl;
		os << "    Iteration type: " << get_description() << endl;
		os << "    Model calls so far : " << run_manager.get_total_runs() << endl << endl << endl;
		iteration_jac(run_manager, termination_ctl, best_upgrade_run, false, restart_runs);

		//write final parameters for restarting
		output_file_writer.write_par(file_manager.open_ofile_ext("par"), best_upgrade_run.get_ctl_pars(), *(par_transform.get_offset_ptr()),
			*(par_transform.get_scale_ptr()));
		file_manager.close_file("par");
		return best_upgrade_run;
	}


ModelRun SVDSolver::iteration_reuse_jac(RunManagerAbstract &run_manager, TerminationController &termination_ctl, ModelRun &base_run,
	bool rerun_base, const string &jco_filename, const string &res_filename)
{
	ModelRun new_base_run(base_run);
	ostream &os = file_manager.rec_ofstream();

	vector<string> numeric_parname_vec = par_transform.ctl2numeric_cp(new_base_run.get_ctl_pars()).get_keys();


	string jac_filename = jco_filename;
	if (jco_filename.empty()) jac_filename = file_manager.build_filename("jcb");
	cout << "  reading previously computed jacobian:  " << jac_filename << endl;
	file_manager.get_ofstream("rec") << "  reading previously computed jacobian:  " << jac_filename << endl;
	jacobian.read(jac_filename);
	//todo: make sure the jco has the right pars and obs

	if (!res_filename.empty())
	{
		stringstream message;

		message << "  reading  residual file " << res_filename << " for hot-start...";
		cout << message.str();
		file_manager.rec_ofstream() << message.str();

		Observations temp_obs = new_base_run.get_obs_template();
		string rfile = res_filename;
		read_res(rfile, temp_obs);
		Parameters temp_pars = new_base_run.get_ctl_pars();
		new_base_run.update_ctl(temp_pars, temp_obs);
		rerun_base = false;
		//message.clear();
		//message << "done" << endl;
		cout << "done" << endl;
		file_manager.rec_ofstream() << "done" << endl;
	}

	if (rerun_base)
	{
		cout << endl << endl;
		cout << "  running the model once with the current parameters... " << endl;
		run_manager.reinitialize(file_manager.build_filename("rnr"));
		int run_id = run_manager.add_run(par_transform.ctl2model_cp(new_base_run.get_ctl_pars()));
		run_manager.run();
		Parameters tmp_pars;
		Observations tmp_obs;
		bool success = run_manager.get_run(run_id, tmp_pars, tmp_obs);
		if (!success)
		{
			throw(PestError("Error: Base parameter run failed.  Can not continue."));
		}

		//Update parameters and observations for base run
		par_transform.model2ctl_ip(tmp_pars);
		new_base_run.update_ctl(tmp_pars, tmp_obs);

		//if at least one parameter in jacobian par names is in the base run, then update the jacobian pars and obs
		bool found = false;
		for (auto &par : tmp_pars)
		{
			auto i = find(jacobian.parameter_list().begin(), jacobian.parameter_list().end(), par.first);

			if (i != jacobian.parameter_list().end())
			{
				found = true;
				break;
			}
		}

		if (found)
		{
			jacobian.set_base_numeric_pars(par_transform.ctl2numeric_cp(new_base_run.get_ctl_pars()));
			jacobian.set_base_sim_obs(new_base_run.get_obs());
		}


	}

	//jacobian.save("jcb");
	output_file_writer.write_jco(true, "jcb", jacobian);
	// sen file for this iteration
	output_file_writer.append_sen(file_manager.sen_ofstream(), termination_ctl.get_iteration_number() + 1,
		jacobian, *(new_base_run.get_obj_func_ptr()), get_parameter_group_info(), *regul_scheme_ptr, false, par_transform);
	cout << endl;
	return new_base_run;
}

void SVDSolver::iteration_jac(RunManagerAbstract &run_manager, TerminationController &termination_ctl, ModelRun &base_run, bool calc_init_obs, bool restart_runs)
{
	ostream &os = file_manager.rec_ofstream();
	ostream &fout_restart = file_manager.get_ofstream("rst");

	set<string> out_ofbound_pars;

	vector<string> numeric_parname_vec = par_transform.ctl2numeric_cp(base_run.get_ctl_pars()).get_keys();

	if (!restart_runs)
	{

		// Calculate Jacobian
		if (!base_run.obs_valid() || calc_init_obs == true) {
			calc_init_obs = true;
		}
		cout << "  calculating jacobian... ";
		performance_log->log_event("commencing to build jacobian parameter sets");
		jacobian.build_runs(base_run, numeric_parname_vec, par_transform,
			*par_group_info_ptr, *ctl_par_info_ptr, run_manager, out_ofbound_pars,
			phiredswh_flag, calc_init_obs);

		RestartController::write_jac_runs_built(fout_restart);
	}

	performance_log->log_event("jacobian parameter sets built, commencing model runs");
	jacobian.make_runs(run_manager);
	performance_log->log_event("jacobian runs complete, processing runs");
	jacobian.process_runs(par_transform,
		*par_group_info_ptr, run_manager, *prior_info_ptr, splitswh_flag);
	performance_log->log_event("processing jacobian runs complete");

	performance_log->log_event("saving jacobian and sen files");
	// save jacobian
	//jacobian.save("jcb");
	output_file_writer.write_jco(true,"jcb", jacobian);

	//Update parameters and observations for base run
	{
		Parameters tmp_pars;
		Observations tmp_obs;
		bool success = run_manager.get_run(0, tmp_pars, tmp_obs);
		par_transform.model2ctl_ip(tmp_pars);
		base_run.update_ctl(tmp_pars, tmp_obs);
	}
	// sen file for this iteration
	output_file_writer.append_sen(file_manager.sen_ofstream(), termination_ctl.get_iteration_number() + 1, jacobian,
		*(base_run.get_obj_func_ptr()), get_parameter_group_info(), *regul_scheme_ptr, false, par_transform);
}

ModelRun SVDSolver::iteration_upgrd(RunManagerAbstract &run_manager, TerminationController &termination_ctl, ModelRun &base_run, bool restart_runs)
{
	ostream &os = file_manager.rec_ofstream();
	ostream &fout_restart = file_manager.get_ofstream("rst");

	if (restart_runs)
	{
		run_manager.initialize_restart(file_manager.build_filename("rnu"));
		{
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
	}
	else
	{
		vector<string> obs_names_vec;
		Observations obs = base_run.get_obs();
		for (auto &o : base_run.get_obs_template().get_keys())
		{
			if (base_run.get_obj_func_ptr()->get_obs_info_ptr()->get_weight(o) > 0.0)
				obs_names_vec.push_back(o);
		}

		//Freeze Parameter for which the jacobian could not be calculated
		auto &failed_jac_pars_names = jacobian.get_failed_parameter_names();
		auto  failed_jac_pars = base_run.get_ctl_pars().get_subset(failed_jac_pars_names.begin(), failed_jac_pars_names.end());

		// populate vectors with sorted observations (standard and prior info) and parameters
		{
			vector<string> prior_info_names = prior_info_ptr->get_keys();
			obs_names_vec.insert(obs_names_vec.end(), prior_info_names.begin(), prior_info_names.end());
		}

		// build weights matrix sqrt(Q)
		QSqrtMatrix Q_sqrt(obs_info_ptr, prior_info_ptr);
		//build residuals vector
		VectorXd residuals_vec = -1.0 * stlvec_2_egienvec(base_run.get_residuals_vec(obs_names_vec));

		Parameters base_run_active_ctl_par = par_transform.ctl2active_ctl_cp(base_run.get_ctl_pars());

		//If running in regularization mode, adjust the regularization weights
		// define a function type for upgrade methods
		{
			Parameters tmp_new_par;
			Parameters frozen_active_ctl_pars = failed_jac_pars;
			LimitType limit_type;
			//use call to calc_upgrade_vec to compute frozen parameters
			calc_upgrade_vec_freeze(0, frozen_active_ctl_pars, Q_sqrt, *regul_scheme_ptr, residuals_vec,
				obs_names_vec, base_run_active_ctl_par,
				tmp_new_par, mar_mat, false);
			if (regul_scheme_ptr->get_use_dynamic_reg() && reg_frac < 0.0)
			{
				dynamic_weight_adj(base_run, jacobian, Q_sqrt, residuals_vec, obs_names_vec,
					base_run_active_ctl_par, frozen_active_ctl_pars);
			}
		}

		//Build model runs
		run_manager.reinitialize(file_manager.build_filename("rnu"));
		// Save base run as first model run so it is eassily accessible
		Parameters base_model_pars = par_transform.ctl2model_cp(base_run.get_ctl_pars());
		int run_id = run_manager.add_run(base_model_pars, "base_run");
		run_manager.update_run(run_id, base_model_pars, base_run.get_obs());
		//Marquardt Lambda Update Vector
		vector<double> lambda_vec = base_lambda_vec;
		if (upgrade_augment)
		{
			lambda_vec.push_back(best_lambda);
			lambda_vec.push_back(best_lambda / 2.0);
			lambda_vec.push_back(best_lambda * 2.0);
		}
		std::sort(lambda_vec.begin(), lambda_vec.end());
		auto iter = std::unique(lambda_vec.begin(), lambda_vec.end());
		lambda_vec.resize(std::distance(lambda_vec.begin(), iter));
		int i_update_vec = 0;
		stringstream message;
		stringstream prf_message;

		ofstream &fout_frz = file_manager.open_ofile_ext("fpr");
		ofstream &fout_rec = file_manager.rec_ofstream();
		for (double i_lambda : lambda_vec)
		{
			prf_message.str("");
			prf_message << "beginning upgrade vector calculations, lambda = " << i_lambda;
			performance_log->log_event(prf_message.str());
			performance_log->add_indent();
			std::cout << string(message.str().size(), '\b');
			message.str("");
			message << "  computing upgrade vector (lambda = " << i_lambda << ")  " << ++i_update_vec << " / " << lambda_vec.size() << "             ";
			std::cout << message.str();
			fout_rec << message.str();

			Parameters new_pars;
			// reset frozen_active_ctl_pars
			LimitType limit_type;
			Parameters frozen_active_ctl_pars = failed_jac_pars;
			try
			{
				calc_upgrade_vec(i_lambda, frozen_active_ctl_pars, Q_sqrt, *regul_scheme_ptr, residuals_vec,
					obs_names_vec, base_run_active_ctl_par, new_pars, mar_mat, limit_type, false);
			}
			catch (const runtime_error& error)
			{
				stringstream ss;
				ss << "error computing upgrade vector for lambda " << i_lambda << ", skipping : " << error.what();
				fout_rec << endl << endl << ss.str() << endl;
				std::cout << endl << ss.str() << endl;
				performance_log->log_event(ss.str());
				performance_log->add_indent(-1);
				continue;
			}
			performance_log->log_event("transforming upgrade pars to model pars");
			Parameters new_par_model = par_transform.active_ctl2model_cp(new_pars);
			performance_log->log_event("writing upgrade vector to csv file");
			output_file_writer.write_upgrade(termination_ctl.get_iteration_number(), 0, i_lambda, 1.0, new_pars);
			int run_id = run_manager.add_run(new_par_model, "normal", i_lambda);
			save_frozen_pars(fout_frz, frozen_active_ctl_pars, run_id);

			//Add Scaled Upgrade Vectors
			Parameters new_numeric_pars = par_transform.model2numeric_cp(new_par_model);
			Parameters base_numeric_pars = par_transform.model2numeric_cp(base_model_pars);
			Parameters del_numeric_pars = new_numeric_pars - base_numeric_pars;
			for (double i_scale : lambda_scale_vec)
			{
				if (i_scale == 1.0)
					continue;
				Parameters scaled_pars = base_numeric_pars + del_numeric_pars * i_scale;
				par_transform.numeric2model_ip(scaled_pars);
				Parameters scaled_ctl_pars = par_transform.numeric2ctl_cp(scaled_pars);
				output_file_writer.write_upgrade(termination_ctl.get_iteration_number(),
					0, i_lambda, i_scale, scaled_ctl_pars);

				stringstream ss;
				ss << "scale(" << std::fixed << std::setprecision(2) << i_scale << ")";
				int run_id = run_manager.add_run(scaled_pars, ss.str(), i_lambda);
				fout_rec << "   ...calculating scaled lambda vector-scale factor: " << i_scale << endl;
				save_frozen_pars(fout_frz, frozen_active_ctl_pars, run_id);
			}

			////Try to extend the previous upgrade vector
			if ((upgrade_augment) &&  (limit_type == LimitType::LBND || limit_type == LimitType::UBND))
			{
				Parameters frozen_active_ctl_pars = failed_jac_pars;
				calc_upgrade_vec_freeze(i_lambda, frozen_active_ctl_pars, Q_sqrt, *regul_scheme_ptr, residuals_vec,
					obs_names_vec, base_run_active_ctl_par,
					new_pars, mar_mat, false);

				Parameters new_par_model = par_transform.active_ctl2model_cp(new_pars);
				int run_id = run_manager.add_run(new_par_model, "extended", i_lambda);
				save_frozen_pars(fout_frz, frozen_active_ctl_pars, run_id);
			}
			performance_log->add_indent(-1);

		}
		file_manager.close_file("fpr");
		RestartController::write_upgrade_runs_built(fout_restart);
	}

	cout << endl;
	performance_log->add_indent(-1);
	cout << "  performing upgrade vector model runs... ";
	run_manager.run();

	// process model runs
	cout << "  testing upgrade vectors... ";

	ifstream &fin_frz = file_manager.open_ifile_ext("fpr");
	bool best_run_updated_flag = false;

	Parameters base_run_active_ctl_par_tmp = par_transform.ctl2active_ctl_cp(base_run.get_ctl_pars());
	ModelRun best_upgrade_run(base_run);

	os << "  Summary of upgrade runs:" << endl;

	int n_runs = run_manager.get_nruns();
	bool one_success = false;
	for (int i = 1; i < n_runs; ++i) {
		ModelRun upgrade_run(base_run);
		Parameters tmp_pars;
		Observations tmp_obs;
		string lambda_type;
		double i_lambda;
		//This must be outside the loop to insure all parrameter sets are read in order
		bool success = run_manager.get_run(i, tmp_pars, tmp_obs, lambda_type, i_lambda);
		if (success)
		{
			one_success = true;
			par_transform.model2ctl_ip(tmp_pars);
			upgrade_run.update_ctl(tmp_pars, tmp_obs);

			Parameters frozen_pars = read_frozen_pars(fin_frz, i);
			upgrade_run.set_frozen_ctl_parameters(frozen_pars);

			par_transform.ctl2active_ctl_ip(tmp_pars);
			double magnitude = Transformable::l2_norm(base_run_active_ctl_par_tmp, tmp_pars);
			streamsize n_prec = os.precision(2);
			os << "    Lambda = ";
			os << setiosflags(ios::fixed) << setw(8) << i_lambda;
			os << "; Type: " << setw(12) << lambda_type;
			os << "; length = " << std::scientific << magnitude;
			os << setiosflags(ios::fixed);
			os.unsetf(ios_base::floatfield); // reset all flags to default
			os << ";  phi = " <<  setw(9) << std::setprecision(4) << upgrade_run.get_phi(*regul_scheme_ptr);
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
				best_lambda = i_lambda;
			}
		}
		else
		{
			streamsize n_prec = os.precision(2);
			os << "      Marquardt Lambda = ";
			os << setiosflags(ios::fixed) << setw(4) << i_lambda;
			os << "; length = " << "NA";
			os.precision(n_prec);
			os.unsetf(ios_base::floatfield); // reset all flags to default
			os << ";  run failed" << endl;
		}
	}
	file_manager.close_file("fpr");

	// Print frozen parameter information
	const Parameters &frz_ctl_pars = best_upgrade_run.get_frozen_ctl_pars();

	if (frz_ctl_pars.size() > 0)
	{
		vector<string> keys = frz_ctl_pars.get_keys();
		std::sort(keys.begin(), keys.end());
		os << endl;
		os << "    Parameters frozen during best upgrade:" << endl;
		for (auto &ikey : keys)
		{
			auto iter = frz_ctl_pars.find(ikey);
			if (iter != frz_ctl_pars.end())
			{
				os << "      " << iter->first << " frozen at " << iter->second << endl;
			}
		}
	}

	// clean up run_manager memory
	run_manager.free_memory();

	if (!one_success)
	{
		throw runtime_error("all upgrade runs failed.");
	}

	return best_upgrade_run;
}

void SVDSolver::check_limits(const Parameters &init_active_ctl_pars, const Parameters &upgrade_active_ctl_pars,
	map<string, LimitType> &limit_type_map, Parameters &active_ctl_parameters_at_limit)
{
	const string *name;
	double p_init;
	double p_upgrade;
	double b_facorg_lim;
	pair<bool, double> par_limit;
	const ParameterRec *p_info;

	for (auto &ipar : upgrade_active_ctl_pars)
	{
		name = &(ipar.first);  // parameter name
		p_info = ctl_par_info_ptr->get_parameter_rec_ptr(*name);
		if (p_info->is_active())
		{
			par_limit = pair<bool, double>(false, 0.0);
			p_upgrade = ipar.second;  // upgrade parameter value
			p_init = init_active_ctl_pars.get_rec(*name);

			double init_value = ctl_par_info_ptr->get_parameter_rec_ptr(*name)->init_value;
			if (init_value == 0.0)
			{
				init_value = ctl_par_info_ptr->get_parameter_rec_ptr(*name)->ubnd / 4.0;
			}
			b_facorg_lim = ctl_info->facorig * init_value;
			if (abs(p_init) >= b_facorg_lim) {
				b_facorg_lim = p_init;
			}

			// Check Relative Chanage Limit
			if (p_info->chglim == "RELATIVE" && abs((p_upgrade - p_init) / b_facorg_lim) > ctl_info->relparmax)
			{
				par_limit.first = true;
				par_limit.second = p_init + sign(p_upgrade - p_init) * ctl_info->relparmax *  abs(b_facorg_lim);
				limit_type_map[*name] = LimitType::REL;
			}

			// Check Factor Change Limit
			else if (p_info->chglim == "FACTOR") {
				if (b_facorg_lim > 0 && p_upgrade < b_facorg_lim / ctl_info->facparmax)
				{
					par_limit.first = true;
					par_limit.second = b_facorg_lim / ctl_info->facparmax;
					limit_type_map[*name] = LimitType::FACT;
				}
				else if (b_facorg_lim > 0 && p_upgrade > b_facorg_lim*ctl_info->facparmax)
				{
					par_limit.first = true;
					par_limit.second = b_facorg_lim * ctl_info->facparmax;
					limit_type_map[*name] = LimitType::FACT;
				}
				else if (b_facorg_lim < 0 && p_upgrade < b_facorg_lim*ctl_info->facparmax)
				{
					par_limit.first = true;
					par_limit.second = b_facorg_lim * ctl_info->facparmax;
					limit_type_map[*name] = LimitType::FACT;
				}
				else if (b_facorg_lim < 0 && p_upgrade > b_facorg_lim / ctl_info->facparmax)
				{
					par_limit.first = true;
					par_limit.second = b_facorg_lim / ctl_info->facparmax;
					limit_type_map[*name] = LimitType::FACT;
				}
			}
			// Check parameter upper bound
			if ((!par_limit.first && p_upgrade > p_info->ubnd) ||
				(par_limit.first && par_limit.second > p_info->ubnd)) {
				par_limit.first = true;
				par_limit.second = p_info->ubnd;
				limit_type_map[*name] = LimitType::UBND;
			}
			// Check parameter lower bound
			else if ((!par_limit.first && p_upgrade < p_info->lbnd) ||
				(par_limit.first && par_limit.second < p_info->lbnd)) {
				par_limit.first = true;
				par_limit.second = p_info->lbnd;
				limit_type_map[*name] = LimitType::LBND;
			}
			// Add any limited parameters to model_parameters_at_limit
			if (par_limit.first) {
				active_ctl_parameters_at_limit.insert(*name, par_limit.second);
			}
		}
	}
}

Parameters SVDSolver::limit_parameters_freeze_all_ip(const Parameters &init_active_ctl_pars,
	Parameters &upgrade_active_ctl_pars, const Parameters &prev_frozen_active_ctl_pars)
{
	map<string, LimitType> limit_type_map;
	Parameters limited_ctl_parameters;
	const string *name;
	double p_init;
	double p_upgrade;
	double p_limit;
	Parameters new_frozen_active_ctl_parameters;
	pair<bool, double> par_limit;

	//remove frozen parameters
	upgrade_active_ctl_pars.erase(prev_frozen_active_ctl_pars);

	check_limits(init_active_ctl_pars, upgrade_active_ctl_pars, limit_type_map, limited_ctl_parameters);
	// Remove parameters at their upper and lower bound limits as these will be frozen
	vector<string> pars_at_bnds;
	for (auto ipar : limit_type_map)
	{
		if (ipar.second == LimitType::LBND || ipar.second == LimitType::UBND)
		{
			pars_at_bnds.push_back(ipar.first);
		}
	}
	limited_ctl_parameters.erase(pars_at_bnds);

	// Calculate most stringent limit factor on a PEST parameter
	double limit_factor = 1.0;
	double tmp_limit;
	string limit_parameter_name = "";
	Parameters init_numeric_pars = par_transform.active_ctl2numeric_cp(init_active_ctl_pars);
	Parameters upgrade_numeric_pars = par_transform.active_ctl2numeric_cp(upgrade_active_ctl_pars);
	Parameters numeric_parameters_at_limit = par_transform.active_ctl2numeric_cp(limited_ctl_parameters);
	for (auto &ipar : numeric_parameters_at_limit)
	{
		name = &(ipar.first);
		p_limit = ipar.second;
		p_init = init_numeric_pars.get_rec(*name);
		p_upgrade = upgrade_numeric_pars.get_rec(*name);
		tmp_limit = (p_limit - p_init) / (p_upgrade - p_init);
		if (tmp_limit < limit_factor)  {
			limit_factor = tmp_limit;
			limit_parameter_name = *name;
		}
	}
	// Apply limit factor to PEST upgrade parameters
	if (limit_factor != 1.0)
	{
		for (auto &ipar : upgrade_numeric_pars)
		{
			name = &(ipar.first);
			p_init = init_numeric_pars.get_rec(*name);
			ipar.second = p_init + (ipar.second - p_init) *  limit_factor;
		}
	}

	//Transform parameters back their ative control state and freeze any that violate their bounds
	upgrade_active_ctl_pars = par_transform.numeric2active_ctl_cp(upgrade_numeric_pars);

	check_limits(init_active_ctl_pars, upgrade_active_ctl_pars, limit_type_map, limited_ctl_parameters);
	for (auto &ipar : upgrade_active_ctl_pars)
	{
		name = &(ipar.first);
		if (limit_type_map[*name] == LimitType::UBND)
		{
			double limit_value = ctl_par_info_ptr->get_parameter_rec_ptr(*name)->ubnd;
			new_frozen_active_ctl_parameters[*name] = limit_value;
		}
		else if (limit_type_map[*name] == LimitType::LBND)
		{
			double limit_value = ctl_par_info_ptr->get_parameter_rec_ptr(*name)->lbnd;
			new_frozen_active_ctl_parameters[*name] = limit_value;
		}
	}

	// Impose frozen Parameters
	for (auto &ipar : prev_frozen_active_ctl_pars)
	{
		upgrade_active_ctl_pars[ipar.first] = ipar.second;
	}

	for (auto &ipar : new_frozen_active_ctl_parameters)
	{
		upgrade_active_ctl_pars[ipar.first] = ipar.second;
	}
	return new_frozen_active_ctl_parameters;
}

void SVDSolver::param_change_stats(double p_old, double p_new, bool &have_fac, double &fac_change, bool &have_rel, double &rel_change)
{
	have_rel = have_fac = true;
	double a = max(abs(p_new), abs(p_old));
	double b = min(abs(p_new), abs(p_old));
	// compute relative change
	if (p_old == 0) {
		have_rel = false;
		rel_change = -9999;
	}
	else
	{
		rel_change = (p_old - p_new) / p_old;
	}
	//compute factor change
	if (p_old == 0.0 || p_new == 0.0) {
		have_fac = false;
		fac_change = -9999;
	}
	else {
		fac_change = a / b;
	}
}

void SVDSolver::iteration_update_and_report(ostream &os, const ModelRun &base_run, ModelRun &upgrade, TerminationController &termination_ctl, RunManagerAbstract &run_manager)
{
	const string *p_name;
	double p_old, p_new;
	double rel_change = -9999;
	bool have_fac = false, have_rel = false;
	double max_rel_change = 0;
	const Parameters &old_ctl_pars = base_run.get_ctl_pars();
	const Parameters &new_ctl_pars = upgrade.get_ctl_pars();

	if (termination_ctl.get_iteration_number() == 0)
	{
		PhiData phi_data = obj_func->phi_report(base_run.get_obs(), base_run.get_ctl_pars(), *regul_scheme_ptr);
		output_file_writer.write_obj_iter(0, 0, phi_data);
	}
	PhiData phi_data = obj_func->phi_report(upgrade.get_obs(), upgrade.get_ctl_pars(), *regul_scheme_ptr);

	output_file_writer.phi_report(os, termination_ctl.get_iteration_number(), run_manager.get_total_runs(), phi_data, regul_scheme_ptr->get_weight(), false, "Final");
	output_file_writer.par_report(os, termination_ctl.get_iteration_number()+1, new_ctl_pars, old_ctl_pars, "Control File");
	output_file_writer.write_obj_iter(termination_ctl.get_iteration_number() + 1, run_manager.get_total_runs(), phi_data);
	for (const auto &ipar : new_ctl_pars)
	{
		p_name = &(ipar.first);
		p_new = ipar.second;
		p_old = old_ctl_pars.get_rec(*p_name);
		rel_change = (p_old - p_new) / p_old;
		if ((p_old != 0.0) && abs(rel_change) >= abs(max_rel_change))
		{
			max_rel_change = rel_change;
		}
	}
	const Parameters old_numeric_pars = par_transform.ctl2numeric_cp(base_run.get_ctl_pars());
	const Parameters new_numeric_pars = par_transform.ctl2numeric_cp(upgrade.get_ctl_pars());
	output_file_writer.par_report(os, termination_ctl.get_iteration_number()+1, new_numeric_pars, old_numeric_pars, "Transformed Numeric");

	//Need to pass defualt regularization weights so this comparision is consistent across iterations
	termination_ctl.process_iteration(upgrade.get_phi_comp(DynamicRegularization::get_unit_reg_instance()), max_rel_change);
}

bool SVDSolver::par_heading_out_bnd(double p_org, double p_del, double lower_bnd, double upper_bnd)
{
	bool out_of_bnd = false;
	double tolerance = 1.0e-5;
	if (((1.0 + tolerance) * p_org > upper_bnd && p_del > 0) || ((1.0 - tolerance) * p_org < lower_bnd && p_del < 0))
	{
		out_of_bnd = true;
	}
	return out_of_bnd;
}

int SVDSolver::check_bnd_par(Parameters &new_freeze_active_ctl_pars, const Parameters &current_active_ctl_pars,
	const Parameters &del_upgrade_active_ctl_pars, const Parameters &del_grad_active_ctl_pars)
{
	int num_upgrade_out_grad_in = 0;
	double p_org;
	double p_del;
	double upper_bnd;
	double lower_bnd;
	const string *name_ptr;
	const auto it_end = del_upgrade_active_ctl_pars.end();
	for (const auto &ipar : current_active_ctl_pars)
	{
		name_ptr = &(ipar.first);
		const auto it = del_upgrade_active_ctl_pars.find(*name_ptr);

		if (it != it_end)
		{
			//first check upgrade parameters
			p_del = it->second;
			p_org = current_active_ctl_pars.get_rec(*name_ptr);
			upper_bnd = ctl_par_info_ptr->get_parameter_rec_ptr(*name_ptr)->ubnd;
			lower_bnd = ctl_par_info_ptr->get_parameter_rec_ptr(*name_ptr)->lbnd;
			//these are active parameters so this is not really necessary - just being extra safe
			bool par_active = ctl_par_info_ptr->get_parameter_rec_ptr(*name_ptr)->is_active();
			bool par_going_out = par_heading_out_bnd(p_org, p_del, lower_bnd, upper_bnd);
			//if gradient parameters are provided, also check these
			if (par_active &&  par_going_out && del_grad_active_ctl_pars.size() > 0)
			{
				const auto it_grad = del_grad_active_ctl_pars.find(*name_ptr);
				if (it_grad != del_grad_active_ctl_pars.end())
				{
					p_del = it_grad->second;
					par_going_out = par_heading_out_bnd(p_org, p_del, lower_bnd, upper_bnd);
				}
				else
				{
					++num_upgrade_out_grad_in;
				}
			}
			if (par_going_out)
				new_freeze_active_ctl_pars.insert(*name_ptr, p_org);
		}
	}
	return num_upgrade_out_grad_in;
}

void SVDSolver::limit_parameters_ip(const Parameters &init_active_ctl_pars, Parameters &upgrade_active_ctl_pars,
	LimitType &limit_type, const Parameters &frozen_active_ctl_pars)
{
	map<string, LimitType> limit_type_map;
	limit_type = LimitType::NONE;
	const string *name;
	double p_init;
	double p_upgrade;
	double p_limit;

	Parameters limited_active_ctl_parameters;
	//remove forozen parameters from upgrade pars
	upgrade_active_ctl_pars.erase(frozen_active_ctl_pars);

	check_limits(init_active_ctl_pars, upgrade_active_ctl_pars, limit_type_map, limited_active_ctl_parameters);

	////delete any limits cooresponding to ignored types
	//for (auto it = limited_active_ctl_parameters.begin(); it != limited_active_ctl_parameters.end();)
	//{
	//	const string &name = (*it).first;
	//	const LimitType l_type = limit_type_map[name];

	//	auto temp_it = it;
	//	++temp_it;

	//	if (l_type == LimitType::LBND || l_type == LimitType::UBND)
	//	{
	//		limited_active_ctl_parameters.erase(it);
	//	}

	//	it = temp_it;
	//}
	// Calculate most stringent limit factor on a numeric PEST parameters
	double limit_factor = 1.0;
	double tmp_limit;
	string limit_parameter_name = "";
	Parameters limited_numeric_parameters = par_transform.active_ctl2numeric_cp(limited_active_ctl_parameters);
	//this can be optimized to just compute init_numeric_parameters for those parameters at their limits
	Parameters init_numeric_pars = par_transform.active_ctl2numeric_cp(init_active_ctl_pars);
	Parameters upgrade_numeric_pars = par_transform.active_ctl2numeric_cp(upgrade_active_ctl_pars);
	for (auto &ipar : limited_numeric_parameters)
	{
		name = &(ipar.first);
		p_limit = ipar.second;
		p_init = init_numeric_pars.get_rec(*name);
		p_upgrade = upgrade_numeric_pars.get_rec(*name);
		tmp_limit = (p_limit - p_init) / (p_upgrade - p_init);
		if (tmp_limit < limit_factor)
		{
			limit_factor = tmp_limit;
			limit_parameter_name = *name;
			limit_type = limit_type_map[*name];
		}
	}
	if (limit_factor == 0.0)
	{
		cout << endl <<  "WARNING: zero-length upgrade vector resulting from parameter bounds enforcement" << endl;
		file_manager.rec_ofstream() << "WARNING: zero-length upgrade vector resulting from parameter bounds enforcement" << endl;
	}
	// Apply limit factor to numeric PEST upgrade parameters
	if (limit_factor != 1.0)
	{
		for (auto &ipar : upgrade_numeric_pars)
		{
			name = &(ipar.first);
			p_init = init_numeric_pars.get_rec(*name);
			ipar.second = p_init + (ipar.second - p_init) *  limit_factor;
		}
	}
	//Convert newly limited parameters to their derivative state
	upgrade_active_ctl_pars = par_transform.numeric2active_ctl_cp(upgrade_numeric_pars);
	// Impose frozen Parameters as they were removed in the beginning
	for (auto &ipar : frozen_active_ctl_pars)
	{
		upgrade_active_ctl_pars[ipar.first] = ipar.second;
	}
	return;
}

PhiComponets SVDSolver::phi_estimate(const ModelRun &base_run, const Jacobian &jacobian, QSqrtMatrix &Q_sqrt, const DynamicRegularization &regul,
	const Eigen::VectorXd &residuals_vec, const vector<string> &obs_names_vec,
	const Parameters &base_run_active_ctl_par, const Parameters &freeze_active_ctl_pars,
	DynamicRegularization &regul_scheme, bool scale_upgrade)
{
	Parameters upgrade_ctl_del_pars;
	Parameters grad_ctl_del_pars;
	Parameters new_pars;

	typedef void(SVDSolver::*UPGRADE_FUNCTION) (const Jacobian &jacobian, const QSqrtMatrix &Q_sqrt, const DynamicRegularization &regul,
		const Eigen::VectorXd &Residuals, const vector<string> &obs_name_vec,
		const Parameters &base_ctl_pars, const Parameters &prev_frozen_ctl_pars,
		double lambda, Parameters &ctl_upgrade_pars, Parameters &upgrade_ctl_del_pars,
		Parameters &grad_ctl_del_pars, MarquardtMatrix marquardt_type, bool scale_upgrade);

	UPGRADE_FUNCTION calc_lambda_upgrade = &SVDSolver::calc_lambda_upgrade_vec_JtQJ;

	if (mat_inv == MAT_INV::Q12J)
	{
		calc_lambda_upgrade = &SVDSolver::calc_lambda_upgrade_vecQ12J;
	}

	(*this.*calc_lambda_upgrade)(jacobian, Q_sqrt, regul, residuals_vec, obs_names_vec,
		base_run_active_ctl_par, freeze_active_ctl_pars, 0, new_pars, upgrade_ctl_del_pars,
		grad_ctl_del_pars, MarquardtMatrix::IDENT, scale_upgrade);

	//Don't limit parameters as this is just an estimate
	//limit_parameters_ip(base_run_active_ctl_par, new_pars,
	//	limit_type, freeze_active_ctl_pars, true);

	Parameters new_ctl_pars = par_transform.active_ctl2ctl_cp(new_pars);
	Parameters delta_par = par_transform.active_ctl2numeric_cp(new_pars)
		- par_transform.active_ctl2numeric_cp(base_run_active_ctl_par);
	vector<string> numeric_par_names = delta_par.get_keys();
	VectorXd delta_par_vec = transformable_2_egien_vec(delta_par, numeric_par_names);
	Eigen::SparseMatrix<double> jac = jacobian.get_matrix(obs_names_vec, numeric_par_names);
	VectorXd delta_obs_vec = jac * delta_par_vec;
	//Transformable delta_obs(obs_names_vec, delta_obs_vec);
	Observations projected_obs = base_run.get_obs();
	Observations full_delta_obs(projected_obs);
	full_delta_obs.update_without_clear(obs_names_vec,delta_obs_vec);
	projected_obs += full_delta_obs;

	PhiComponets proj_phi_comp = base_run.get_obj_func_ptr()->get_phi_comp(projected_obs, new_ctl_pars, regul_scheme);
	return proj_phi_comp;
}


void SVDSolver::dynamic_weight_adj_percent(const ModelRun &base_run, double reg_frac)
{
	ostream &os = file_manager.rec_ofstream();

	double phimlim = regul_scheme_ptr->get_phimlim();
	double fracphim = regul_scheme_ptr->get_fracphim();
	double wfmin = regul_scheme_ptr->get_wfmin();
	double wfmax = regul_scheme_ptr->get_wfmax();
	double wffac = regul_scheme_ptr->get_wffac();
	double wf_cur = regul_scheme_ptr->get_weight();
	double wftol = regul_scheme_ptr->get_wftol();

	Parameters new_pars;
	vector<MuPoint> mu_vec;
	mu_vec.resize(4);

	// Equalize Reqularization Groups if IREGADJ = 1
	if (regul_scheme_ptr->get_adj_grp_weights())
	{
		std::unordered_map<std::string, double> regul_grp_weights;
		auto reg_grp_phi = base_run.get_obj_func_ptr()->get_group_phi(base_run.get_obs(), base_run.get_ctl_pars(),
			DynamicRegularization::get_unit_reg_instance(), PhiComponets::OBS_TYPE::REGUL);
		double avg_reg_grp_phi = 0;
		for (const auto &igrp : reg_grp_phi)
		{
			avg_reg_grp_phi += igrp.second;
		}
		if (reg_grp_phi.size() > 0)
		{
			avg_reg_grp_phi /= reg_grp_phi.size();
		}
		if (avg_reg_grp_phi > 0)
		{
			for (const auto &igrp : reg_grp_phi)
			{
				if (igrp.second > 0)
				{
					regul_grp_weights[igrp.first] = sqrt(avg_reg_grp_phi / igrp.second);
				}
			}
		}
		regul_scheme_ptr->set_regul_grp_weights(regul_grp_weights);
	}

	DynamicRegularization tmp_regul_scheme = *regul_scheme_ptr;
	PhiComponets phi_comp_cur = base_run.get_obj_func_ptr()->get_phi_comp(base_run.get_obs(), base_run.get_ctl_pars(), *regul_scheme_ptr);

	double wf_new;
	//if (phi_comp_cur.regul <= wfmin || reg_frac == 1.0)
	if (reg_frac >= 1.0)
	{
		wf_new = wfmin;
	}
	else if (phi_comp_cur.regul < wfmin)
	{
		wf_new = reg_frac * phi_comp_cur.meas / (1.0 - reg_frac);
		wf_new *= wf_cur;
	}
	else
	{
		wf_new = reg_frac * phi_comp_cur.meas / (phi_comp_cur.regul * (1.0 - reg_frac));
		wf_new *= wf_cur;
	}

	wf_new = max(wf_new, wfmin);
	wf_new = min(wf_new, wfmax);
	regul_scheme_ptr->set_weight(wf_new);

	PhiComponets phi_comp_new = base_run.get_obj_func_ptr()->get_phi_comp(base_run.get_obs(), base_run.get_ctl_pars(), *regul_scheme_ptr);

	os << "  Regularization weight factor update (REG_FRAC = " << reg_frac << ")" << endl;
	os << "    Initial values:" << endl;
	os << "      regularization weight factor =  " << wf_cur << endl;
	os << "      PHI (total=" << phi_comp_cur.meas + phi_comp_cur.regul <<
		"; meas=" << phi_comp_cur.meas << "; reg=" << phi_comp_cur.regul << ")" << endl;
	os << "    Updated values:" << endl;
	os << "      regularization weight factor = " << wf_new << endl;
	os << "      PHI (total=" << phi_comp_new.meas + phi_comp_new.regul <<
		"; meas=" << phi_comp_new.meas << "; reg=" << phi_comp_new.regul << ")" << endl;
}


void SVDSolver::dynamic_weight_adj(const ModelRun &base_run, const Jacobian &jacobian, QSqrtMatrix &Q_sqrt,
	const Eigen::VectorXd &residuals_vec, const vector<string> &obs_names_vec,
	const Parameters &base_run_active_ctl_par, const Parameters &freeze_active_ctl_pars)
{
	ostream &os = file_manager.rec_ofstream();

	double phimlim = regul_scheme_ptr->get_phimlim();
	double fracphim = regul_scheme_ptr->get_fracphim();
	double wfmin = regul_scheme_ptr->get_wfmin();
	double wfmax = regul_scheme_ptr->get_wfmax();
	double wffac = regul_scheme_ptr->get_wffac();
	double mu_cur = regul_scheme_ptr->get_weight();
	double wftol = regul_scheme_ptr->get_wftol();
	int max_iter = regul_scheme_ptr->get_max_reg_iter();
	Parameters new_pars;
	vector<MuPoint> mu_vec;
	mu_vec.resize(4);

	// Equalize Reqularization Groups if IREGADJ = 1
	std::unordered_map<std::string, double> regul_grp_weights;
	auto reg_grp_phi = base_run.get_obj_func_ptr()->get_group_phi(base_run.get_obs(), base_run.get_ctl_pars(),
		DynamicRegularization::get_unit_reg_instance(), PhiComponets::OBS_TYPE::REGUL);
	double avg_reg_grp_phi = 0;
	for (const auto &igrp : reg_grp_phi)
	{
		avg_reg_grp_phi += igrp.second;
	}
	if (reg_grp_phi.size() > 0)
	{
		avg_reg_grp_phi /= reg_grp_phi.size();
	}
	if (avg_reg_grp_phi>0)
	{
		for (const auto &igrp : reg_grp_phi)
		{
			if (igrp.second > 0)
			{
				regul_grp_weights[igrp.first] = sqrt(avg_reg_grp_phi / igrp.second);
			}
		}
	}
	regul_scheme_ptr->set_regul_grp_weights(regul_grp_weights);

	DynamicRegularization tmp_regul_scheme = *regul_scheme_ptr;
	PhiComponets phi_comp_cur = base_run.get_obj_func_ptr()->get_phi_comp(base_run.get_obs(), base_run.get_ctl_pars(), *regul_scheme_ptr);
	double target_phi_meas_frac = phi_comp_cur.meas * fracphim;
	double target_phi_meas = max(phimlim, target_phi_meas_frac);

	os << endl << "    ---  Solving for regularization weight factor   ---   " << endl;
	os << "    Starting regularization weight factor      : " << mu_cur << endl;
	os << "    Starting measurement objective function    : " << phi_comp_cur.meas << endl;
	os << "    Starting regularization objective function : " << phi_comp_cur.regul << endl;
	os << endl << "    Target measurement objective function      : " << target_phi_meas << endl << endl;
	cout << "    ---  Solving for regularization weight factor   ---   " << endl;
	cout << "    Starting regularization weight factor      : " << mu_cur << endl;
	cout << "    Starting measurement objective function    : " << phi_comp_cur.meas << endl;
	cout << "    Starting regularization objective function : " << phi_comp_cur.regul << endl;
	cout << "    Target measurement objective function      : " << target_phi_meas << endl << endl;

	for (auto &i_mu : mu_vec)
	{
		i_mu.target_phi_meas = target_phi_meas;
	}
	PhiComponets proj_phi_cur = phi_estimate(base_run, jacobian, Q_sqrt, *regul_scheme_ptr, residuals_vec, obs_names_vec,
		base_run_active_ctl_par, freeze_active_ctl_pars, *regul_scheme_ptr);
	double f_cur = proj_phi_cur.meas - target_phi_meas;

	if (f_cur < 0)
	{
		mu_vec[0].set(mu_cur, proj_phi_cur);
		double mu_new = mu_vec[0].mu * wffac;
		mu_new = max(wfmin, mu_new);
		mu_new = min(mu_new, wfmax);
		tmp_regul_scheme.set_weight(mu_new);
		PhiComponets phi_proj_new = phi_estimate(base_run, jacobian, Q_sqrt, tmp_regul_scheme, residuals_vec, obs_names_vec,
			base_run_active_ctl_par, freeze_active_ctl_pars, tmp_regul_scheme);
		mu_vec[3].set(mu_new, phi_proj_new);
	}
	else
	{
		mu_vec[3].set(mu_cur, proj_phi_cur);
		double mu_new = mu_vec[3].mu * wffac;
		mu_new = max(wfmin, mu_new);
		mu_new = min(mu_new, wfmax);
		tmp_regul_scheme.set_weight(mu_new);
		PhiComponets phi_proj_new = phi_estimate(base_run, jacobian, Q_sqrt, tmp_regul_scheme, residuals_vec, obs_names_vec,
			base_run_active_ctl_par, freeze_active_ctl_pars, tmp_regul_scheme);
		mu_vec[0].set(mu_new, phi_proj_new);
	}

	// make sure f[0] and f[3] bracket the solution

	int i;
	for (i = 0; i < max_iter; ++i)
	{
		if (mu_vec[0].f() < 0)
		{
			break;
		}
		else
		{
			mu_vec[3] = mu_vec[0];
			mu_vec[0].mu = max(mu_vec[0].mu / wffac, wfmin);
			tmp_regul_scheme.set_weight(mu_vec[0].mu);
			mu_vec[0].phi_comp = phi_estimate(base_run, jacobian, Q_sqrt, tmp_regul_scheme, residuals_vec, obs_names_vec,
				base_run_active_ctl_par, freeze_active_ctl_pars, tmp_regul_scheme);
			mu_vec[0].print(os);
			os << endl;
			cout << "    ...solving for optimal weight factor : " << setw(6) << mu_vec[0].mu << endl << flush;
			//mu_vec[0].print(cout);
			//cout << endl;
		}
		if (mu_vec[0].mu <= wfmin) break;
	}

	for (; i < max_iter; ++i)
	{
		if (mu_vec[3].f() > 0)
		{
			break;
		}
		else
		{
			mu_vec[0] = mu_vec[3];
			mu_vec[3].mu = min(mu_vec[3].mu * wffac, wfmax);
			tmp_regul_scheme.set_weight(mu_vec[3].mu);
			mu_vec[3].phi_comp = phi_estimate(base_run, jacobian, Q_sqrt, tmp_regul_scheme, residuals_vec, obs_names_vec,
				base_run_active_ctl_par, freeze_active_ctl_pars, tmp_regul_scheme);
			mu_vec[3].print(os);
			os << endl;
			//mu_vec[0].print(cout);
			//cout << endl;
			cout << "    ...solving for optimal weight factor : " << setw(6) << mu_vec[3].mu << endl << flush;
		}
		if (mu_vec[3].mu >= wfmax) break;
	}

	double tau = (sqrt(5.0) - 1.0) / 2.0;
	double lw = mu_vec[3].mu - mu_vec[0].mu;
	mu_vec[1].mu = mu_vec[0].mu + (1.0 - tau) * lw;
	tmp_regul_scheme.set_weight(mu_vec[1].mu);
	mu_vec[1].phi_comp = phi_estimate(base_run, jacobian, Q_sqrt, tmp_regul_scheme, residuals_vec, obs_names_vec,
		base_run_active_ctl_par, freeze_active_ctl_pars, tmp_regul_scheme);

	mu_vec[2].mu = mu_vec[3].mu - (1.0 - tau) * lw;
	tmp_regul_scheme.set_weight(mu_vec[2].mu);
	mu_vec[2].phi_comp = phi_estimate(base_run, jacobian, Q_sqrt, tmp_regul_scheme, residuals_vec, obs_names_vec,
		base_run_active_ctl_par, freeze_active_ctl_pars, tmp_regul_scheme);

	if (mu_vec[0].f() > 0)
	{
		auto min_mu = (mu_vec[0] < mu_vec[3]) ? mu_vec[0] : mu_vec[3];
		os << "    ---  optimal regularization weight factor found  ---   " << endl;
		min_mu.print(os);
		os << endl;
		cout << "    ---  optimal regularization weight factor found  ---   " << endl;
		min_mu.print(cout);
		cout << endl;
		double m = (mu_vec[0] < mu_vec[3]) ? mu_vec[0].mu : mu_vec[3].mu;
		regul_scheme_ptr->set_weight(m);
	}
	else if (mu_vec[3].f() < 0)
	{
		auto min_mu = (mu_vec[0] < mu_vec[3]) ? mu_vec[0] : mu_vec[3];
		os << "    ---  optimal regularization weight factor found  ---   " << endl;
		min_mu.print(os);
		os << endl;
		cout << "    ---  optimal regularization weight factor found  ---   " << endl;
		min_mu.print(cout);
		cout << endl;
		double m = (mu_vec[0] < mu_vec[3]) ? mu_vec[0].mu : mu_vec[3].mu;
		regul_scheme_ptr->set_weight(m);
	}

	else
	{
		for (; i < max_iter; ++i)
		{
			if (abs(mu_vec[0].f()) > abs(mu_vec[3].f()) && mu_vec[0].f() < 0)
			{
				mu_vec[0] = mu_vec[1];
				mu_vec[1] = mu_vec[2];
				lw = mu_vec[3].mu - mu_vec[0].mu;
				mu_vec[2].mu = mu_vec[3].mu - (1.0 - tau) * lw;
				tmp_regul_scheme.set_weight(mu_vec[2].mu);
				mu_vec[2].phi_comp = phi_estimate(base_run, jacobian, Q_sqrt, tmp_regul_scheme, residuals_vec, obs_names_vec,
					base_run_active_ctl_par, freeze_active_ctl_pars, tmp_regul_scheme);
			}
			else
			{
				mu_vec[3] = mu_vec[2];
				mu_vec[2] = mu_vec[1];
				lw = mu_vec[3].mu - mu_vec[0].mu;
				mu_vec[1].mu = mu_vec[0].mu + (1.0 - tau) * lw;
				tmp_regul_scheme.set_weight(mu_vec[1].mu);
				mu_vec[1].phi_comp = phi_estimate(base_run, jacobian, Q_sqrt, tmp_regul_scheme, residuals_vec, obs_names_vec,
					base_run_active_ctl_par, freeze_active_ctl_pars, tmp_regul_scheme);
			}

			auto min_mu = std::min_element(mu_vec.begin(), mu_vec.end());
			min_mu->print(os);
			os << endl;
			//mu_vec[0].print(cout);
			//cout << endl;
			cout << "    ...solving for optimal weight factor : " << setw(6) << min_mu->mu << endl << flush;
			if (min_mu->error_frac() <= wftol) break;
		}
		auto min_mu = std::min_element(mu_vec.begin(), mu_vec.end());

		os << "    ---  optimal regularization weight factor found  ---   " << endl;
		min_mu->print(os);
		os << endl;
		cout << "    ---  optimal regularization weight factor found  ---   " << endl;
		min_mu->print(cout);
		cout << endl;
		double new_mu = std::min_element(mu_vec.begin(), mu_vec.end())->mu;
		regul_scheme_ptr->set_weight(new_mu);
	}
	os << endl;
}

 void SVDSolver::save_frozen_pars(std::ostream &fout, const Parameters &frozen_pars, int id)
 {
		 fout << "frozen_parameter_set_begin  " << id << endl;
		 fout << frozen_pars;
		 fout << "frozen_parameter_set_end  " << id << endl;
		 fout.flush();
 }

 Parameters SVDSolver::read_frozen_pars(std::istream &fin, int id)
 {
	 Parameters fz_pars;
	 string line;
	 vector<string> tokens;
	 double value;
	 int cur_read_id = 0;
	 bool failed = true;

	 while (getline(fin, line))
	 {
		strip_ip(line);
		tokens.clear();
		tokenize(line, tokens);

		if (tokens[0] == "frozen_parameter_set_begin")
		{
			fz_pars.clear();
		}
		else if (tokens[0] == "frozen_parameter_set_end")
		{
			failed = false;
			convert_ip(tokens[1], cur_read_id);
			if (cur_read_id == id)
			{
				break;
			}
			else
			{
				failed = true;
			}
		}
		else
		{
			convert_ip(tokens[2], value);
			fz_pars.insert(tokens[0], value);
		}
	 }
	 if (failed)
	 {
		 throw PestError("Error: reading past end of frozen parameter file (SVDSolver::read_frozen_pars)");
	 }
	 return fz_pars;
 }

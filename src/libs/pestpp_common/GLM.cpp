#include "GLM.h"

#include <iostream>
#include <sstream>

#include "eigen_tools.h"
#include "linear_analysis.h"
#include "utilities.h"

using namespace std;

/// The regularization scheme, creating the default one first if the control file had no
/// regularization section.
///
/// Pest leaves regul_scheme_ptr NULL until a '* regularization' section sets it, and the
/// termination controller below reads through that pointer - so on a control file without one
/// this is a null dereference. pest++.cpp used to paper over it by calling set_default_dynreg()
/// in main() before building anything, which meant the precondition was the CALLER's to know:
/// the executable satisfied it and a library caller segfaulted. Owning it here is what makes
/// "main and the API run the same loop" true rather than aspirational.
static DynamicRegularization& ensure_dynreg(Pest& scenario)
{
	if (scenario.get_regul_scheme_ptr() == nullptr)
		scenario.set_default_dynreg();
	return *scenario.get_regul_scheme_ptr();
}

GLM::GLM(Pest& _pest_scenario, FileManager& _file_manager,
	OutputFileWriter& _output_file_writer, PerformanceLog* _performance_log,
	RunManagerAbstract* _run_mgr_ptr,
	RestartController* _restart_ctl, bool _restart_flag, bool _save_restart_rec_header)
	: pest_scenario(_pest_scenario), file_manager(_file_manager),
	output_file_writer(_output_file_writer), performance_log(_performance_log),
	run_mgr_ptr(_run_mgr_ptr), restart_ctl(_restart_ctl), restart_flag(_restart_flag),
	save_restart_rec_header(_save_restart_rec_header),
	obj_func(&(_pest_scenario.get_ctl_observations()),
		&(_pest_scenario.get_ctl_observation_info()),
		&(_pest_scenario.get_prior_info())),
	rand_gen(_pest_scenario.get_pestpp_options().get_random_seed()),
	termination_ctl(_pest_scenario.get_control_info().noptmax,
		_pest_scenario.get_control_info().phiredstp,
		_pest_scenario.get_control_info().nphistp,
		_pest_scenario.get_control_info().nphinored,
		_pest_scenario.get_control_info().relparstp,
		_pest_scenario.get_control_info().nrelpar,
		ensure_dynreg(_pest_scenario).get_use_dynamic_reg(),
		ensure_dynreg(_pest_scenario).get_phimaccept()),
	cur_run(&obj_func, _pest_scenario.get_ctl_observations()),
	optimum_run(&obj_func, _pest_scenario.get_ctl_observations()),
	noptmax(_pest_scenario.get_control_info().noptmax),
	initialized(false)
{
}

void GLM::throw_glm_error(const string& message)
{
	// flush, never close - see the note on the declaration
	file_manager.rec_ofstream() << "GLM error: " << message << endl;
	file_manager.rec_ofstream().flush();
	cout << "GLM error: " << message << endl;
	throw runtime_error("GLM error: " + message);
}

void GLM::initialize()
{
	int n_runs = initialize_prepare();
	if (n_runs > 0)
		run_mgr_ptr->run();
	initialize_finish();
}

int GLM::initialize_prepare()
{
	const ParamTransformSeq& base_trans_seq = pest_scenario.get_base_par_tran_seq();

	// The glm output files (.sen, .svd, the iteration summaries) and the svd output option.
	// main() used to do this before constructing anything, which made it another precondition
	// the executable happened to satisfy and a library caller did not - the first upgrade then
	// died on 'Error accessing file: "pest.svd"'. Same reasoning as ensure_dynreg().
	output_file_writer.prep_glm_files(restart_flag);
	output_file_writer.set_svd_output_opt(pest_scenario.get_svd_info().eigwrite);
	// iteration 0's parameter record, which has to follow the prep above: .ipar does not
	// exist until prepare_iteration_summary_files() opens it. main() used to do both, in that
	// order, and separating them is what broke pestpp-glm with 'Error accessing file:
	// "pest.ipar"' - so they travel together.
	// write_par_iter() checks ++iteration_summary() itself, so no gate here
	output_file_writer.write_par_iter(0, pest_scenario.get_ctl_parameters());

	base_jacobian_ptr.reset(new Jacobian_1to1(file_manager, output_file_writer));

	if (restart_flag && (restart_ctl != nullptr))
		restart_ctl->update_termination_ctl(termination_ctl);

	file_manager.rec_ofstream() << "...loading prior parameter covariance matrix" << endl << endl;
	performance_log->log_event("loading parcov");
	parcov.try_from(pest_scenario, file_manager);

	base_svd_ptr.reset(new SVDSolver(pest_scenario, file_manager, &obj_func, base_trans_seq,
		*base_jacobian_ptr, output_file_writer, performance_log, parcov, &rand_gen,
		"base parameter solution"));
	base_svd_ptr->set_svd_package(pest_scenario.get_pestpp_options().get_svd_pack());
	cout << endl << endl;

	cur_ctl_parameters = pest_scenario.get_ctl_parameters();

	// Allocates space for the run manager. This initializes the model parameter and
	// observation names, neither of which changes over the course of the simulation.
	if ((restart_ctl != nullptr) &&
		(restart_ctl->get_restart_option() == RestartController::RestartOption::RESUME_JACOBIAN_RUNS))
	{
		run_mgr_ptr->initialize_restart(file_manager.build_filename("rns"));
	}
	else
	{
		run_mgr_ptr->initialize(base_trans_seq.ctl2model_cp(cur_ctl_parameters),
			pest_scenario.get_ctl_observations());
	}

	// noptmax==0 is one forward run with the initial parameters, queued here so a caller can
	// service it. Every other path computes its first batch inside the first iteration, so
	// there is nothing to hand over.
	if (noptmax == 0)
	{
		Parameters init_model_pars = base_trans_seq.ctl2model_cp(cur_ctl_parameters);
		optimum_run.set_ctl_parameters(init_model_pars);
		run_mgr_ptr->reinitialize();
		run_mgr_ptr->add_run(init_model_pars);
		return 1;
	}
	return 0;
}

void GLM::initialize_finish()
{
	if (noptmax == 0)
	{
		run_initial_parameters_only();
		initialized = true;
		return;
	}

	// Define the model run for the base parameters, using the base parameter transformations
	cur_run.set_ctl_parameters(cur_ctl_parameters);
	if ((restart_ctl != nullptr) &&
		(restart_ctl->get_restart_option() != RestartController::RestartOption::NONE))
	{
		Parameters restart_pars = restart_ctl->get_restart_parameters(
			file_manager.build_filename("parb"), file_manager.build_filename("par"));
		if (restart_pars.size() > 0)
			cur_run.set_ctl_parameters(restart_pars);
	}
	if (!restart_flag || save_restart_rec_header)
		file_manager.rec_ofstream() << "   -----    Starting pestpp-glm Iterations    ----    "
			<< endl << endl;
	initialized = true;
}

void GLM::run_initial_parameters_only()
{
	const ParamTransformSeq& base_trans_seq = pest_scenario.get_base_par_tran_seq();
	ofstream& fout_rec = file_manager.rec_ofstream();

	Parameters tmp_pars;
	Observations tmp_obs;
	bool success = run_mgr_ptr->get_run(0, tmp_pars, tmp_obs);
	base_trans_seq.model2ctl_ip(tmp_pars);

	if (!success)
	{
		// used to be exit(1), which a library caller cannot survive
		throw_glm_error("model run failed - no results were recorded");
	}

	termination_ctl.check_last_iteration();
	optimum_run.update_ctl(tmp_pars, tmp_obs);
	// save parameters to .par file
	output_file_writer.write_par(file_manager.open_ofile_ext("par"), optimum_run.get_ctl_pars(),
		*(base_trans_seq.get_offset_ptr()), *(base_trans_seq.get_scale_ptr()));
	file_manager.close_file("par");
	// save new residuals to .rei file
	output_file_writer.write_rei(file_manager.open_ofile_ext("rei"), 0,
		*(optimum_run.get_obj_func_ptr()->get_obs_ptr()),
		optimum_run.get_obs(), *(optimum_run.get_obj_func_ptr()),
		optimum_run.get_ctl_pars());
	PhiData pd = optimum_run.get_obj_func_ptr()->phi_report(optimum_run.get_obs(),
		optimum_run.get_ctl_pars(), *(pest_scenario.get_regul_scheme_ptr()));
	output_file_writer.write_obj_iter(0, run_mgr_ptr->get_nruns(), pd);
	file_manager.close_file("rei");
	run_mgr_ptr->free_memory();

	termination_ctl.set_terminate(true);
	(void)fout_rec;
}

int GLM::jacobian_prepare(bool calc_init_obs)
{
	if (!initialized)
		throw_glm_error("jacobian_prepare() called before initialize()");
	return base_svd_ptr->jacobian_prepare(*run_mgr_ptr, cur_run, calc_init_obs, false);
}

void GLM::jacobian_run()
{
	if (!initialized)
		throw_glm_error("jacobian_run() called before initialize()");
	base_svd_ptr->jacobian_run(*run_mgr_ptr);
}

bool GLM::jacobian_process()
{
	if (!initialized)
		throw_glm_error("jacobian_process() called before initialize()");
	return base_svd_ptr->jacobian_process(*run_mgr_ptr, termination_ctl, cur_run);
}

void GLM::iterate_2_solution()
{
	if (!initialized)
		throw_glm_error("iterate_2_solution() called before initialize()");
	while (!termination_ctl.terminate())
	{
		if (!solve_iteration())
			break;
	}
}

bool GLM::solve_iteration()
{
	int q = pest_utils::quit_file_found();
	if ((q == 1) || (q == 2))
	{
		termination_ctl.set_terminate(true);
		termination_ctl.set_reason("'pest.stp' found");
	}
	if (termination_ctl.terminate())
		return false;

	SVDSolver& base_svd = *base_svd_ptr;
	RestartController::RestartOption restart_option = (restart_ctl != nullptr)
		? restart_ctl->get_restart_option() : RestartController::RestartOption::NONE;

	if (noptmax < 0)
	{
		// a single Jacobian and then stop - no upgrade is attempted
		if (restart_option == RestartController::RestartOption::REUSE_JACOBIAN)
		{
			try
			{
				string jco_filename = pest_scenario.get_pestpp_options().get_basejac_filename();
				jco_filename = ((jco_filename.empty()) ? file_manager.build_filename("jco") : jco_filename);
				string res_filename = pest_scenario.get_pestpp_options().get_hotstart_resfile();

				cur_run = base_svd.iteration_reuse_jac(*run_mgr_ptr, termination_ctl, cur_run,
					true, jco_filename, res_filename);
				if (!cur_run.obs_valid())
					cur_run = base_svd.solve(*run_mgr_ptr, termination_ctl, noptmax, cur_run,
						optimum_run, *restart_ctl, false);
			}
			catch (const exception& e)
			{
				throw_glm_error(string("error restarting with existing jco: ") + e.what());
			}
		}
		else
		{
			try
			{
				bool restart_runs =
					(restart_option == RestartController::RestartOption::RESUME_JACOBIAN_RUNS);
				cur_run = base_svd.compute_jacobian(*run_mgr_ptr, termination_ctl, cur_run,
					restart_runs);
				if (restart_runs && (restart_ctl != nullptr))
					restart_ctl->get_restart_option() = RestartController::RestartOption::NONE;
			}
			catch (const exception& e)
			{
				throw_glm_error(string("error filling base Jacobian: ") + e.what());
			}
		}
		optimum_run = cur_run;
		output_file_writer.write_rei(file_manager.open_ofile_ext("rei"), -1,
			pest_scenario.get_ctl_observations(), cur_run.get_obs(), *cur_run.get_obj_func_ptr(),
			pest_scenario.get_ctl_parameters());
		termination_ctl.set_terminate(true);
		termination_ctl.set_reason("NOPTMAX criterion met");
	}
	else if (restart_option == RestartController::RestartOption::REUSE_JACOBIAN)
	{
		try
		{
			bool calc_first_jacobian = false;
			string jco_filename = pest_scenario.get_pestpp_options().get_basejac_filename();
			jco_filename = ((jco_filename.empty()) ? file_manager.build_filename("jco") : jco_filename);
			string res_filename = pest_scenario.get_pestpp_options().get_hotstart_resfile();

			cur_run = base_svd.iteration_reuse_jac(*run_mgr_ptr, termination_ctl, cur_run, true,
				jco_filename, res_filename);
			// run the model once with the current parameters to compute the observations
			cur_run = base_svd.solve(*run_mgr_ptr, termination_ctl, noptmax, cur_run,
				optimum_run, *restart_ctl, calc_first_jacobian);
			termination_ctl.check_last_iteration();
		}
		catch (const exception& e)
		{
			throw_glm_error(string("error in iteration process with existing Jacobian: ") + e.what());
		}
	}
	else
	{
		try
		{
			if (restart_option == RestartController::RestartOption::RESUME_UPGRADE_RUNS)
				base_svd.iteration_reuse_jac(*run_mgr_ptr, termination_ctl, cur_run, false,
					file_manager.build_filename("jcb"));
			bool calc_first_jacobian = true;
			// solve() needs a RestartController even when there is no restart; a default-
			// constructed one is RestartOption::NONE, which is exactly "no restart"
			RestartController none_ctl;
			cur_run = base_svd.solve(*run_mgr_ptr, termination_ctl, noptmax, cur_run,
				optimum_run, (restart_ctl != nullptr) ? *restart_ctl : none_ctl,
				calc_first_jacobian);
			termination_ctl.check_last_iteration();
		}
		catch (const exception& e)
		{
			throw_glm_error(string("error in iteration process: ") + e.what());
		}
	}

	cur_ctl_parameters = cur_run.get_ctl_pars();
	return !termination_ctl.terminate();
}

void GLM::finalize()
{
	ofstream& fout_rec = file_manager.rec_ofstream();

	cout << endl;
	termination_ctl.termination_summary(cout);
	cout << endl;
	termination_ctl.termination_summary(fout_rec);
	fout_rec << endl;

	if ((pest_scenario.get_ctl_ordered_par_names().size() < 1000) &&
		(pest_scenario.get_ctl_ordered_obs_names().size() < 1000))
	{
		cout << "FINAL OPTIMISATION RESULTS" << endl << endl;
		fout_rec << "FINAL OPTIMISATION RESULTS" << endl << endl;

		fout_rec << "  Optimal parameter values  " << endl;
		output_file_writer.par_report(fout_rec, optimum_run.get_ctl_pars());

		fout_rec << endl << "  Observations with optimal model-simulated equivalents and residuals" << endl;
		ObservationInfo oi = pest_scenario.get_ctl_observation_info();
		output_file_writer.obs_report(fout_rec, *obj_func.get_obs_ptr(), optimum_run.get_obs(), oi);
	}

	fout_rec << endl << "Final composite objective function " << endl;
	PhiData phi_data = obj_func.phi_report(optimum_run.get_obs(), optimum_run.get_ctl_pars(),
		*(pest_scenario.get_regul_scheme_ptr()));
	output_file_writer.phi_report(fout_rec, termination_ctl.get_iteration_number() + 1,
		run_mgr_ptr->get_total_runs(), phi_data, 0.0, true);
	output_file_writer.phi_report(cout, termination_ctl.get_iteration_number() + 1,
		run_mgr_ptr->get_total_runs(), phi_data, 0.0, true);
	fout_rec << endl << endl;
	fout_rec << "Number of forward model runs performed during optimization: "
		<< run_mgr_ptr->get_total_runs() << endl;

	if ((noptmax != 0) && (pest_scenario.get_pestpp_options().get_uncert_flag()))
		write_fosm();
}

void GLM::write_fosm()
{
	ofstream& fout_rec = file_manager.rec_ofstream();

	cout << endl << endl << "...starting posterior FOSM calculations..." << endl;

	fout_rec << endl << endl << endl << endl;
	fout_rec << "-----------------------------------------------------------------------" << endl;
	fout_rec << "Note: The following uncertainty estimates were calculated using " << endl;
	fout_rec << "      Schur's complement for linear-based conditional uncertainty " << endl;
	fout_rec << "      propagation.  For a derivation from Bayes equation, see " << endl;
	fout_rec << "      M. N. Fienen, J. E. Doherty, R. J. Hunt, and H. W. Reeves. " << endl;
	fout_rec << "      2010. 'Using Prediction Uncertainty Analysis to Design " << endl;
	fout_rec << "      Hydrologic Monitoring Networks: Example Applications " << endl;
	fout_rec << "      from the Great Lakes Water Availability Pilot Project'. " << endl;
	fout_rec << "      See PEST++ V3 documentation for implementation details." << endl;
	fout_rec << "-----------------------------------------------------------------------" << endl;
	fout_rec << endl;

	fout_rec << "Note: Any observations or prior information equations with a group name" << endl;
	fout_rec << "      starting with 'regul' are dropped from the Jacobian and observation" << endl;
	fout_rec << "      covariance matrices before uncertainty calculations.  Please" << endl;
	fout_rec << "      make sure that all expert knowledge is expressed in the prior " << endl;
	fout_rec << "      parameter bounds or through a covariance matrix, which can be " << endl;
	fout_rec << "      supplied as a ++ option as '++parcov(<matrix_file_name>)'," << endl;
	fout_rec << "      where <matrix_file_name> can be an ASCII PEST-compatible matrix file (.mat) or" << endl;
	fout_rec << "      a PEST-compatible uncertainty file (.unc)." << endl << endl;

	performance_log->log_event("FOSM-based posterior unc calcs");

	if (base_jacobian_ptr->get_base_numeric_par_names().size() == 0)
	{
		// used to `return 0` straight out of main(), which also skipped the cleanup below it
		cout << "WARNING: no parameters in base Jacobian, can't calculate uncertainty with FOSM" << endl;
		fout_rec << "WARNING: no parameters in base Jacobian, can't calculate uncertainty with FOSM" << endl;
		return;
	}

	//instance of a Mat for the jco
	Mat j(base_jacobian_ptr->get_sim_obs_names(), base_jacobian_ptr->get_base_numeric_par_names(),
		base_jacobian_ptr->get_matrix_ptr());
	if (pest_scenario.get_prior_info_ptr()->get_nnz_pi() > 0)
	{
		vector<string> pi_names = pest_scenario.get_ctl_ordered_pi_names();
		j.drop_rows(pi_names);
	}
	LinearAnalysis la(j, pest_scenario, file_manager, *performance_log, parcov, &rand_gen);
	la.glm_iter_fosm(optimum_run, output_file_writer, -999, run_mgr_ptr);
	if (pest_scenario.get_pestpp_options().get_glm_num_reals() > 0)
	{
		cout << endl << "...drawing and running "
			<< pest_scenario.get_pestpp_options().get_glm_num_reals()
			<< " FOSM-based posterior realizations" << endl;

		pair<ParameterEnsemble, map<string, int>> fosm_real_info =
			la.draw_fosm_reals(run_mgr_ptr, -999, optimum_run);
		run_mgr_ptr->run();
		DynamicRegularization ptr;
		ptr.set_defaults();
		ptr.set_weight(0.0);
		double phi = optimum_run.get_phi(ptr);
		pair<ObservationEnsemble, map<string, double>> fosm_obs_info =
			la.process_fosm_reals(run_mgr_ptr, fosm_real_info, -999, phi);
	}
}

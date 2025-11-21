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
#include "RunManagerPanther.h" //needs to be first because it includes winsock2.h
//#include <vld.h> // Memory Leak Detection using "Visual Leak Detector"
#include <iostream>
#include <fstream>
#include <algorithm>
#include "config_os.h"
#include "Pest.h"
#include "Jacobian_1to1.h"
#include "Transformable.h"
#include "Transformation.h"
#include "ParamTransformSeq.h"
#include "utilities.h"
#include "pest_error.h"
#include "ModelRunPP.h"
#include "SVDASolver.h"
#include  "QSqrtMatrix.h"
#include "FileManager.h"
#include "TerminationController.h"
#include "RunManagerSerial.h"
#include "RunManagerExternal.h"
#include "SVD_PROPACK.h"
#include "OutputFileWriter.h"
#include "PantherAgent.h"
#include "Serialization.h"
#include "system_variables.h"
#include "pest_error.h"
#include "RestartController.h"
#include "PerformanceLog.h"
#include "debug.h"
#include "DifferentialEvolution.h"

#include "linear_analysis.h"
#include "logger.h"
#include "covariance.h"
#include "Ensemble.h"
#include "eigen_tools.h"


using namespace std;
using namespace pest_utils;

//using namespace pest_utils;

int main(int argc, char* argv[])
{
#ifndef _DEBUG
	try
	{
#endif
		string version = PESTPP_VERSION;
		cout << endl << endl;
		cout << "             pestpp-glm: a tool for GLM parameter estimation and FOSM uncertainty analysis" << endl << endl;
		cout << "                                   by The PEST++ Development Team" << endl;
		cout << endl;
        auto start = chrono::steady_clock::now();

        string start_string = get_time_string();
		CmdLine cmdline(argc, argv);

        if (quit_file_found())
        {
            cerr << "'pest.stp' found, please remove this file " << endl;
            return 1;
        }
		
		FileManager file_manager;
		string filename = cmdline.ctl_file_name;
		string pathname = ".";
		file_manager.initialize_path(get_filename_without_ext(filename), pathname);

	
		if (cmdline.runmanagertype == CmdLine::RunManagerType::PANTHER_WORKER)
		{
			try
			{
				ofstream frec("panther_worker.rec");
				if (frec.bad())
					throw runtime_error("error opening 'panther_worker.rec'");
				cmdline.startup_report(frec,start_string);
				cmdline.startup_report(cout,start_string);

				PANTHERAgent yam_agent(frec);
				string ctl_file = "";
				try 
				{	
					// process traditional PEST control file
					ctl_file = file_manager.build_filename("pst");
					yam_agent.process_ctl_file(ctl_file);
				}
				catch (exception &e)
				{
                    frec << "Error processing control file: " << ctl_file << endl << endl;
                    frec << e.what() << endl << endl;
					cerr << "Error prococessing control file: " << ctl_file << endl << endl;
					cerr << e.what() << endl << endl;
					throw(e);
				}

				yam_agent.start(cmdline.panther_host_name,cmdline.panther_port);
			}
			catch (PestError &perr)
			{
				cerr << perr.what();
				throw(perr);
			}
			cout << endl << "Work Done..." << endl;
			exit(0);
		}
		//Check for PANTHER Master
		if (cmdline.runmanagertype == CmdLine::RunManagerType::GENIE)
		{
			cerr << "Genie run manager ('/g') no longer supported, please use PANTHER instead" << endl;
			exit(1);
		}
		
		RestartController restart_ctl;
		bool restart_flag = false;
		bool save_restart_rec_header = true;

		debug_initialize(file_manager.build_filename("dbg"));
		if (cmdline.jac_restart)
		{
			cout << endl << "ERROR: '/j' restart option is deprecated.  Please use ++base_jacobian() instead." << endl << endl;
			//restart_ctl.get_restart_option() = RestartController::RestartOption::REUSE_JACOBIAN;
			//file_manager.open_default_files();
		}
		else if (cmdline.restart)
		{
			ifstream &fin_rst = file_manager.open_ifile_ext("rst");
			if (fin_rst.bad())
            {
			    throw runtime_error("restart error: error opening rst file '"+file_manager.get_base_filename()+".rst'");
            }
			restart_ctl.process_rst_file(fin_rst);
			file_manager.close_file("rst");
			restart_flag = true;
			file_manager.open_default_files(true);
			ofstream &fout_rec_tmp = file_manager.rec_ofstream();
			fout_rec_tmp << endl << endl;
			if (cmdline.runmanagertype == CmdLine::RunManagerType::EXTERNAL)
			{
				save_restart_rec_header = false;
			}
			else
			{
				fout_rec_tmp << "Restarting pestpp-glm ....." << endl << endl;
				cout << "    Restarting pestpp-glm ....." << endl << endl;
			}
		}
		else
		{
			restart_ctl.get_restart_option() = RestartController::RestartOption::NONE;
			file_manager.open_default_files();
		}

		ofstream &fout_rec = file_manager.rec_ofstream();
		PerformanceLog performance_log(file_manager.open_ofile_ext("log"));


		fout_rec << "             pestpp-glm " << version << endl << endl;
		fout_rec << "    by The PEST++ Development Team" << endl;

		fout_rec << endl;
		cmdline.startup_report(fout_rec,start_string);
		cmdline.startup_report(cout,start_string);




		// create pest run and process control file to initialize it
		Pest pest_scenario;
        pest_scenario.set_default_dynreg();
#ifndef _NDEBUG
		try {
#endif
			performance_log.log_event("starting to process control file");
			pest_scenario.process_ctl_file(file_manager.open_ifile_ext("pst"), file_manager.build_filename("pst"),fout_rec);
			file_manager.close_file("pst");
			performance_log.log_event("finished processing control file");
#ifndef _NDEBUG
		}
		catch (exception &e)
		{
            fout_rec << "Error processing control file: " << filename << endl << endl;
            fout_rec << e.what() << endl << endl;
			cerr << "Error processing control file: " << filename << endl << endl;
			cerr << e.what() << endl << endl;
			throw(e);
	 	}
#endif

		pest_scenario.check_inputs(fout_rec, false, false);
		// reset this here because we want to draw from the FOSM posterior as a whole matrix
		pest_scenario.get_pestpp_options_ptr()->set_ies_group_draws(false);
		
		if (pest_scenario.get_pestpp_options().get_glm_normal_form() == PestppOptions::GLMNormalForm::PRIOR)
		{
			if (pest_scenario.get_control_info().pestmode == ControlInfo::PestMode::REGUL)
			{
				throw runtime_error("'GLM_NORMAL_FORM' = 'PRIOR' is incompatible with 'PESTMODE' = 'REGULARIZATION'");
			}
		}

		if (pest_scenario.get_control_info().pestmode == ControlInfo::PestMode::REGUL)
        {
		    if (pest_scenario.get_prior_info().get_nnz_pi() == 0)
            {
		        throw runtime_error("regularization mode requires at least one non-zero weighted prior info equation");
            }
		    if ((pest_scenario.get_pestpp_options().get_glm_iter_mc()) &&
		    (pest_scenario.get_pestpp_options().get_glm_accept_mc_phi()))
            {
		        stringstream ss;
		        ss << endl << "WARNING 'regularization' mode is not conceptually compatible with 'glm_accept_mc_phi'" << endl;
		        cout << ss.str();
		        fout_rec << ss.str();
            }

        }

		//Initialize OutputFileWriter to handle IO of supplementary files (.par, .par, .svd)
		//bool save_eign = pest_scenario.get_svd_info().eigwrite > 0;
		OutputFileWriter output_file_writer(file_manager, pest_scenario, restart_flag);
		
		if (!restart_flag)
		{
			output_file_writer.scenario_report(fout_rec);
		}
		if (pest_scenario.get_pestpp_options().get_debug_parse_only())
		{
			cout << endl << endl << "DEBUG_PARSE_ONLY is true, exiting..." << endl << endl;
			exit(0);
		}
		output_file_writer.prep_glm_files(restart_flag);
		output_file_writer.set_svd_output_opt(pest_scenario.get_svd_info().eigwrite);
		//if base jco arg read from control file, reset restart controller
		if (!pest_scenario.get_pestpp_options().get_basejac_filename().empty())
		{
			restart_ctl.get_restart_option() = RestartController::RestartOption::REUSE_JACOBIAN;
		}

		if (pest_scenario.get_pestpp_options().get_iter_summary_flag())
		{
			output_file_writer.write_par_iter(0, pest_scenario.get_ctl_parameters());
		}
		RunManagerAbstract *run_manager_ptr;
		if (cmdline.runmanagertype == CmdLine::RunManagerType::PANTHER_MASTER)
		{
			// using PANTHER run manager
			if (pest_scenario.get_control_info().noptmax == 0)
			{
				cout << endl << endl << "WARNING: 'noptmax' = 0 but using parallel run mgr.  This prob isn't what you want to happen..." << endl << endl;
			}
			const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
			//check for condor wrapper
			string csf = pest_scenario.get_pestpp_options().get_condor_submit_file();
			if (csf.size() > 0)
			{
				if (!pest_utils::check_exist_in(csf))
					throw runtime_error("++condor_submit_file '" + csf + "' not found");
				run_manager_ptr = new RunManagerYAMRCondor(
					file_manager.build_filename("rns"), cmdline.panther_port,
					file_manager.open_ofile_ext("rmr"),
					pest_scenario.get_pestpp_options().get_max_run_fail(),
					pest_scenario.get_pestpp_options().get_overdue_reched_fac(),
					pest_scenario.get_pestpp_options().get_overdue_giveup_fac(),
					pest_scenario.get_pestpp_options().get_overdue_giveup_minutes(),
					csf);
			}
			else
			{
				run_manager_ptr = new RunManagerPanther(
					file_manager.build_filename("rns"), cmdline.panther_port,
					file_manager.open_ofile_ext("rmr"),
					pest_scenario.get_pestpp_options().get_max_run_fail(),
					pest_scenario.get_pestpp_options().get_overdue_reched_fac(),
					pest_scenario.get_pestpp_options().get_overdue_giveup_fac(),
					pest_scenario.get_pestpp_options().get_overdue_giveup_minutes(),
					pest_scenario.get_pestpp_options().get_panther_echo(),
                    vector<string>{}, vector<string>{},
                    pest_scenario.get_pestpp_options().get_panther_timeout_milliseconds(),
                    pest_scenario.get_pestpp_options().get_panther_echo_interval_milliseconds(),
                    pest_scenario.get_pestpp_options().get_panther_persistent_workers());
			}
		}
		
		else if (cmdline.runmanagertype == CmdLine::RunManagerType::EXTERNAL)
		{
			const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
			string rns_file = file_manager.build_filename("rns");
			run_manager_ptr = new RunManagerExternal(exi.comline_vec,
				exi.tplfile_vec, exi.inpfile_vec, exi.insfile_vec, exi.outfile_vec,
				rns_file,
				pest_scenario.get_pestpp_options().get_max_run_fail());
		}
		else
		{
			performance_log.log_event("starting basic model IO error checking");
			cout << "checking model IO files...";
			pest_scenario.check_io(fout_rec);
			performance_log.log_event("finished basic model IO error checking");
			cout << "done" << endl;
			const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
			run_manager_ptr = new RunManagerSerial(exi.comline_vec,
				exi.tplfile_vec, exi.inpfile_vec, exi.insfile_vec, exi.outfile_vec,
				file_manager.build_filename("rns"), pathname,
				pest_scenario.get_pestpp_options().get_max_run_fail(),
				pest_scenario.get_pestpp_options().get_fill_tpl_zeros(),
				pest_scenario.get_pestpp_options().get_additional_ins_delimiters(),
				pest_scenario.get_pestpp_options().get_num_tpl_ins_threads(),
				pest_scenario.get_pestpp_options().get_tpl_force_decimal());
		}

		const ParamTransformSeq &base_trans_seq = pest_scenario.get_base_par_tran_seq();
		ObjectiveFunc obj_func(&(pest_scenario.get_ctl_observations()), &(pest_scenario.get_ctl_observation_info()), &(pest_scenario.get_prior_info()));
		Jacobian *base_jacobian_ptr = new Jacobian_1to1(file_manager,output_file_writer);
		std::mt19937 rand_gen(pest_scenario.get_pestpp_options().get_random_seed());
		TerminationController termination_ctl(pest_scenario.get_control_info().noptmax, pest_scenario.get_control_info().phiredstp,
			pest_scenario.get_control_info().nphistp, pest_scenario.get_control_info().nphinored, pest_scenario.get_control_info().relparstp,
			pest_scenario.get_control_info().nrelpar, pest_scenario.get_regul_scheme_ptr()->get_use_dynamic_reg(),
			pest_scenario.get_regul_scheme_ptr()->get_phimaccept());

		//if we are doing a restart, update the termination_ctl
		if (restart_flag)
		{
			restart_ctl.update_termination_ctl(termination_ctl);
		}

		file_manager.rec_ofstream() << "...loading prior parameter covariance matrix" << endl << endl;
		performance_log.log_event("loading parcov");
		Covariance parcov;
		parcov.try_from(pest_scenario, file_manager);
		
		SVDSolver base_svd(pest_scenario, file_manager, &obj_func, base_trans_seq,
			*base_jacobian_ptr, output_file_writer, &performance_log, parcov, &rand_gen, "base parameter solution");

		base_svd.set_svd_package(pest_scenario.get_pestpp_options().get_svd_pack());
		//Build Super-Parameter problem
		Jacobian *super_jacobian_ptr = new Jacobian(file_manager);
		ParamTransformSeq trans_svda;
		// method must be involked as pointer as the transformation sequence it is added to will
		// take responsibility for destroying it
		TranSVD *tran_svd = new TranSVD(pest_scenario.get_pestpp_options().get_max_n_super(),
			pest_scenario.get_pestpp_options().get_super_eigthres(), "SVD Super Parameter Transformation");

		if (pest_scenario.get_pestpp_options().get_svd_pack() != PestppOptions::SVD_PACK::REDSVD)
		{
			tran_svd->set_SVD_pack();
		}
		tran_svd->set_performance_log(&performance_log);

		TranFixed *tr_svda_fixed = new TranFixed("SVDA Fixed Parameter Transformation");
		trans_svda = base_trans_seq;
		trans_svda.push_back_ctl2active_ctl(tr_svda_fixed);
		trans_svda.add_custom_tran_seq(string("svda_derv2basenumeric"), trans_svda.get_ctl2active_ctl_tranformations());
		trans_svda.push_back_active_ctl2numeric(tran_svd);

		ControlInfo svd_control_info = pest_scenario.get_control_info();
		svd_control_info.relparmax = pest_scenario.get_pestpp_options().get_super_relparmax();
		// Start Solution iterations
		cout << endl << endl;
		int n_base_iter = pest_scenario.get_pestpp_options().get_n_iter_base();
		int n_super_iter = pest_scenario.get_pestpp_options().get_n_iter_super();
		int max_n_super = pest_scenario.get_pestpp_options().get_max_n_super();
		double super_eigthres = pest_scenario.get_pestpp_options().get_super_eigthres();

		Parameters cur_ctl_parameters = pest_scenario.get_ctl_parameters();
		//Allocates Space for Run Manager.  This initializes the model parameter names and observations names.
		//Neither of these will change over the course of the simulation

		if (restart_ctl.get_restart_option() == RestartController::RestartOption::RESUME_JACOBIAN_RUNS)
		{
			run_manager_ptr->initialize_restart(file_manager.build_filename("rnj"));
		}
		else
		{
			run_manager_ptr->initialize(base_trans_seq.ctl2model_cp(cur_ctl_parameters), pest_scenario.get_ctl_observations());
		}


		ModelRun optimum_run(&obj_func, pest_scenario.get_ctl_observations());
		//for tracking the initial model simulated equivalents
		//vector<double> init_sim;
		// if noptmax=0 make one run with the initial parameters
		if (pest_scenario.get_control_info().noptmax == 0) {
			Parameters init_model_pars = base_trans_seq.ctl2model_cp(cur_ctl_parameters);
			optimum_run.set_ctl_parameters(init_model_pars);
			run_manager_ptr->reinitialize();
			run_manager_ptr->add_run(init_model_pars);
			try
			{
				run_manager_ptr->run();
			}
			catch (exception &e)
			{
				cout << "Model run failed.  No results were recorded." << endl << e.what() << endl;
				fout_rec << "Model run failed.  No results were recorded." << endl << e.what() << endl;
				exit(1);
			}
			Parameters tmp_pars;
			Observations tmp_obs;
			bool success = run_manager_ptr->get_run(0, tmp_pars, tmp_obs);
			base_trans_seq.model2ctl_ip(tmp_pars);
			//termination_ctl.set_terminate(true);

			if (success)
			{
				termination_ctl.check_last_iteration();
				optimum_run.update_ctl(tmp_pars, tmp_obs);
				// save parameters to .par file
				output_file_writer.write_par(file_manager.open_ofile_ext("par"), optimum_run.get_ctl_pars(), *(base_trans_seq.get_offset_ptr()),
					*(base_trans_seq.get_scale_ptr()));
				file_manager.close_file("par");
				// save new residuals to .rei file
				output_file_writer.write_rei(file_manager.open_ofile_ext("rei"), 0,
					*(optimum_run.get_obj_func_ptr()->get_obs_ptr()),
					optimum_run.get_obs(), *(optimum_run.get_obj_func_ptr()),
					optimum_run.get_ctl_pars());
				PhiData pd = optimum_run.get_obj_func_ptr()->phi_report(optimum_run.get_obs(), optimum_run.get_ctl_pars(),
					*(pest_scenario.get_regul_scheme_ptr()));
				output_file_writer.write_obj_iter(0, run_manager_ptr->get_nruns(), pd);
				file_manager.close_file("rei");
				run_manager_ptr->free_memory();
			}
			else
			{
				cout << "Model run failed.  No results were recorded." << endl << endl;
				fout_rec << "Model run failed.  No results were recorded." << endl << endl;
				exit(1);
			}
			termination_ctl.set_terminate(true);
		}


		// Differential Evolution
		if (pest_scenario.get_pestpp_options().get_global_opt() ==  PestppOptions::OPT_DE)
		{
		    throw runtime_error("DE-based global optimization is deprecated in pestpp-glm. please use pestpp-mou");
			/*int rand_seed = 1;
			int np = pest_scenario.get_pestpp_options().get_de_npopulation();
			int max_gen = pest_scenario.get_pestpp_options().get_de_max_gen();
			double f = pest_scenario.get_pestpp_options().get_de_f();
			double cr = pest_scenario.get_pestpp_options().get_de_cr();
			bool dither_f = pest_scenario.get_pestpp_options().get_de_dither_f();
			ModelRun init_run(&obj_func, pest_scenario.get_ctl_observations());
			Parameters cur_ctl_parameters = pest_scenario.get_ctl_parameters();
			run_manager_ptr->reinitialize();
			DifferentialEvolution de_solver(pest_scenario, file_manager, &obj_func,
				base_trans_seq, output_file_writer, &performance_log, rand_seed);
			de_solver.initialize_population(*run_manager_ptr, np);
			de_solver.solve(*run_manager_ptr, restart_ctl, max_gen, f, cr, dither_f, init_run);
			run_manager_ptr->free_memory();
			exit(0);*/
		}


		//Define model Run for Base Parameters (uses base parameter transformations)
		ModelRun cur_run(&obj_func, pest_scenario.get_ctl_observations());

		cur_run.set_ctl_parameters(cur_ctl_parameters);
		//If this is a restart we need to get the latest ctl parameters
		if (restart_ctl.get_restart_option() != RestartController::RestartOption::NONE)
		{
			Parameters restart_pars = restart_ctl.get_restart_parameters(file_manager.build_filename("parb"), file_manager.build_filename("par"));
			if (restart_pars.size() > 0)
			{
				cur_run.set_ctl_parameters(restart_pars);
			}
		}
		if (!restart_flag || save_restart_rec_header)
		{
			fout_rec << "   -----    Starting pestpp-glm Iterations    ----    " << endl << endl;
		}

		while (!termination_ctl.terminate())
		{
            int q = pest_utils::quit_file_found();
            if ((q == 1) || (q == 2))
            {
		        termination_ctl.set_terminate(true);
		        termination_ctl.set_reason("'pest.stp' found");
            }
			//base parameter iterations
			try
			{
				if (restart_ctl.get_restart_option() != RestartController::RestartOption::NONE  &&
					restart_ctl.get_iteration_type() == RestartController::IterationType::SUPER)
				{
					try
					{
						string filename = pest_scenario.get_pestpp_options().get_basejac_filename();
						filename = ((filename.empty()) ? file_manager.build_filename("jcb") : filename);
						base_svd.iteration_reuse_jac(*run_manager_ptr, termination_ctl, cur_run, false, filename);
					}
					catch (exception &e)
					{
						cout << "error restarting super parameter process: " << e.what() << endl;
						throw runtime_error(e.what());
					}
				}
				else if (restart_ctl.get_restart_option() == RestartController::RestartOption::REUSE_JACOBIAN && n_base_iter < 0)
				{
					try
					{
						string filename = pest_scenario.get_pestpp_options().get_basejac_filename();
						string res_filename = pest_scenario.get_pestpp_options().get_hotstart_resfile();
						filename = ((filename.empty()) ? file_manager.build_filename("jco") : filename);
						cur_run = base_svd.iteration_reuse_jac(*run_manager_ptr, termination_ctl, cur_run, true, filename, res_filename);
					}
					catch (exception &e)
					{
						cout << "error restarting with existing jco and n_iter_base < 0: " << e.what() << endl;
						throw runtime_error(e.what());
					}
				}
				else if (n_base_iter < 0)
				{
					bool restart_runs = (restart_ctl.get_restart_option() == RestartController::RestartOption::REUSE_JACOBIAN);
					try
					{
						cur_run = base_svd.compute_jacobian(*run_manager_ptr, termination_ctl, cur_run, restart_runs);
						if (pest_scenario.get_control_info().noptmax < 0)
						{
							optimum_run = cur_run;
							output_file_writer.write_rei(file_manager.open_ofile_ext("rei"), -1, pest_scenario.get_ctl_observations(),
								cur_run.get_obs(), *cur_run.get_obj_func_ptr(), pest_scenario.get_ctl_parameters());
							termination_ctl.set_terminate(true);
							termination_ctl.set_reason("NOPTMAX criterion met");
							//bool success = run_manager_ptr->get_observations_vec(0, init_sim);
						}
					}
					catch (exception &e)
					{
						cout << "error initializing with n_iter_base < 0: " << e.what() << endl;
						throw runtime_error(e.what());
					}
				}
				else if (pest_scenario.get_control_info().noptmax < 0)
				{
					if (restart_ctl.get_restart_option() == RestartController::RestartOption::REUSE_JACOBIAN)
					{
						try
						{
							string jco_filename = pest_scenario.get_pestpp_options().get_basejac_filename();
							jco_filename = ((jco_filename.empty()) ? file_manager.build_filename("jco") : jco_filename);
							string res_filename = pest_scenario.get_pestpp_options().get_hotstart_resfile();

							cur_run = base_svd.iteration_reuse_jac(*run_manager_ptr, termination_ctl, cur_run, true, jco_filename, res_filename);
							if (!cur_run.obs_valid())
								cur_run = base_svd.solve(*run_manager_ptr, termination_ctl, n_base_iter, cur_run, optimum_run, restart_ctl, false);
						}
						catch (exception &e)
						{
							cout << "error restarting with existing jco: " << e.what() << endl;
							throw runtime_error(e.what());
						}
					}
					else
					{
						try
						{


							bool restart_runs = (restart_ctl.get_restart_option() == RestartController::RestartOption::RESUME_JACOBIAN_RUNS);
							cur_run = base_svd.compute_jacobian(*run_manager_ptr, termination_ctl, cur_run, restart_runs);
							if (restart_runs) restart_ctl.get_restart_option() = RestartController::RestartOption::NONE;

							//bool success = run_manager_ptr->get_observations_vec(0, init_sim);
						}
						catch (exception &e)
						{
							cout << "error filling base Jacobian: " << e.what() << endl;
							throw runtime_error(e.what());
						}
					}
					optimum_run = cur_run;
					output_file_writer.write_rei(file_manager.open_ofile_ext("rei"), -1, pest_scenario.get_ctl_observations(),
						cur_run.get_obs(), *cur_run.get_obj_func_ptr(), pest_scenario.get_ctl_parameters());
					termination_ctl.set_terminate(true);
					termination_ctl.set_reason("NOPTMAX criterion met");
				}
				else if (restart_ctl.get_restart_option() == RestartController::RestartOption::REUSE_JACOBIAN)
				{
					try
					{
						bool calc_first_jacobian = false;
						string jco_filename = pest_scenario.get_pestpp_options().get_basejac_filename();
						jco_filename = ((jco_filename.empty()) ? file_manager.build_filename("jco") : jco_filename);
						string res_filename = pest_scenario.get_pestpp_options().get_hotstart_resfile();

						cur_run = base_svd.iteration_reuse_jac(*run_manager_ptr, termination_ctl, cur_run, true, jco_filename, res_filename);
						// Run the model once with the current parameters to compute the observations
						cur_run = base_svd.solve(*run_manager_ptr, termination_ctl, n_base_iter, cur_run, optimum_run, restart_ctl, calc_first_jacobian);
						termination_ctl.check_last_iteration();
					}
					catch (exception &e)
					{
						cout << "error in base parameter iteration process with existing Jacobian: " << e.what() << endl;
						throw runtime_error(e.what());
					}
				}
				else
				{
					try
					{
						if (restart_ctl.get_restart_option() == RestartController::RestartOption::RESUME_UPGRADE_RUNS)
						{
							base_svd.iteration_reuse_jac(*run_manager_ptr, termination_ctl, cur_run, false, file_manager.build_filename("jcb"));
						}
						bool calc_first_jacobian = true;
						cur_run = base_svd.solve(*run_manager_ptr, termination_ctl, n_base_iter, cur_run, optimum_run, restart_ctl, calc_first_jacobian);
						termination_ctl.check_last_iteration();
					}
					catch (exception &e)
					{
						cout << "error in base parameter iteration process: " << e.what() << endl;
						throw runtime_error(e.what());
					}

				}
				//if (termination_ctl.get_iteration_number() == 1)
					//bool success = run_manager_ptr->get_observations_vec(0, init_sim);
				cur_ctl_parameters = cur_run.get_ctl_pars();
				if (termination_ctl.terminate())  break;
			}
			catch (exception &e)
			{
				cout << endl << endl;;
				cout << e.what() << endl;
				fout_rec << endl << endl;
				fout_rec << e.what() << endl;
				cout << "Error encountered, cannot continue" << endl;
				fout_rec << "Error encountered, cannot continue" << endl;
				exit(1);
			}
			// Build Super Parameter or SVDA problem
			try
			{
				const vector<string> &nonregul_obs = pest_scenario.get_nonregul_obs();
				Parameters base_numeric_pars = base_trans_seq.ctl2numeric_cp(cur_ctl_parameters);
				const vector<string> &pars = base_numeric_pars.get_keys();
				QSqrtMatrix Q_sqrt(&(pest_scenario.get_ctl_observation_info()), &pest_scenario.get_prior_info());
				
				if (restart_ctl.get_restart_option() == RestartController::RestartOption::RESUME_JACOBIAN_RUNS)
				{
					//read previously computed super parameter transformation
					ifstream &fin_rst = file_manager.open_ifile_ext("rtj", ios_base::in | ios_base::binary);
					(*tran_svd).read(fin_rst);
					file_manager.close_file("rtj");
				}
				else
				{
					cout << "...forming super parameter transformation, requires forming and factoring JTQJ..." << endl;
					fout_rec << "...forming super parameter transformation, requires forming and factoring JTQJ..." << endl;
					Eigen::SparseMatrix<double> parcov_inv;
					ParamTransformSeq par_transform = pest_scenario.get_base_par_tran_seq();
					map<string, double> dss = pest_scenario.calc_par_dss(*base_jacobian_ptr,par_transform);
					if (pest_scenario.get_pestpp_options().get_glm_normal_form() == PestppOptions::GLMNormalForm::PRIOR)
					{
						
						parcov_inv = *parcov.get(base_jacobian_ptr->get_base_numeric_par_names()).inv().e_ptr();
					}
					performance_log.log_event("updating super parameter transformation, requires formation and SVD of JtQJ");
					
					(*tran_svd).update_reset_frozen_pars(*base_jacobian_ptr, Q_sqrt, base_numeric_pars, max_n_super, super_eigthres, 
						pars, nonregul_obs, parcov_inv, dss, cur_run.get_frozen_ctl_pars());
					(*tr_svda_fixed).reset((*tran_svd).get_frozen_derivative_pars());
				}
				SVDASolver super_svd(pest_scenario, file_manager, &obj_func,
					trans_svda, *super_jacobian_ptr,
					output_file_writer, &performance_log,parcov,&rand_gen,
					base_svd.get_phiredswh_flag(), base_svd.get_splitswh_flag());
				super_svd.set_svd_package(pest_scenario.get_pestpp_options().get_svd_pack());
				//use base jacobian to compute first super jacobian if there was not a super upgrade
				bool calc_first_jacobian = true;
				if (n_base_iter == -1)
				{
					//transform base jacobian to super jacobian
					super_svd.get_jacobian() = base_svd.get_jacobian();
					super_svd.get_jacobian().transform(base_trans_seq, &ParamTransformSeq::jac_numeric2active_ctl_ip);
					super_svd.get_jacobian().transform(trans_svda, &ParamTransformSeq::jac_active_ctl_ip2numeric_ip);
					//rerun base run to account for round off error in super parameters
					if ((cur_run.obs_valid()) && (!pest_scenario.get_pestpp_options().get_glm_rebase_super()))
					{
						fout_rec << "...glm_rebase_super is false, using existing residuals as super-par-truncated base residuals..." << endl;
						cout << "...glm_rebase_super is false, using existing residuals as super-par-truncated base residuals..." << endl;
					}
					else
					{
						cout << "...running super-par-truncated base parameter values once to account for roundoff in super par transformation" << endl;
						fout_rec << "...running super-par-truncated base parameter values once to account for roundoff in super par transformation" << endl;
						cur_run = super_svd.update_run(*run_manager_ptr, cur_run);
					}
					calc_first_jacobian = false;
					//bool success = run_manager_ptr->get_observations_vec(0, init_sim);
				}
				cur_run = super_svd.solve(*run_manager_ptr, termination_ctl, n_super_iter, cur_run, optimum_run, restart_ctl, calc_first_jacobian);
				cur_ctl_parameters = cur_run.get_ctl_pars();
				base_svd.set_phiredswh_flag(super_svd.get_phiredswh_flag());
				base_svd.set_splitswh_flag(super_svd.get_splitswh_flag());
				if (super_svd.local_iteration_terminatated() && n_base_iter == -1)
				{
					n_base_iter = 1;
				}
			}
			catch (exception &e)
			{
				cout << e.what() << endl;
				fout_rec << e.what() << endl;
				cout << "WARNING: super parameter process failed.  Switching to base parameters" << endl << endl;
				fout_rec << "WARNING: super parameter process failed.  Switching to base parameters" << endl << endl;
				ofstream &fout_restart = file_manager.get_ofstream("rst");
				RestartController::write_start_failed_super(fout_restart);
				restart_ctl.get_restart_option() = RestartController::RestartOption::NONE;
				if (pest_scenario.get_pestpp_options().get_n_iter_base() == -1)
				{
					cout << "resetting n_iter_base to 1 since super parameter process failed" << endl;
					fout_rec << "resetting n_iter_base to 1 since super parameter process failed" << endl;
					pest_scenario.get_pestpp_options_ptr()->set_n_iter_base(1);
					n_base_iter = -1;
				}	
			}
		}
		cout << endl;
		termination_ctl.termination_summary(cout);
		cout << endl;
		termination_ctl.termination_summary(fout_rec);
		fout_rec << endl;
        if ((pest_scenario.get_ctl_ordered_par_names().size() < 1000) && (pest_scenario.get_ctl_ordered_obs_names().size() < 1000)) {
            cout << "FINAL OPTIMISATION RESULTS" << endl << endl;
            fout_rec << "FINAL OPTIMISATION RESULTS" << endl << endl;

            fout_rec << "  Optimal parameter values  " << endl;
            output_file_writer.par_report(fout_rec, optimum_run.get_ctl_pars());

            fout_rec << endl << "  Observations with optimal model-simulated equivalents and residuals" << endl;
            ObservationInfo oi = pest_scenario.get_ctl_observation_info();
            output_file_writer.obs_report(fout_rec, *obj_func.get_obs_ptr(), optimum_run.get_obs(), oi);
        }

		fout_rec << endl << "Final composite objective function " << endl;
		PhiData phi_data = obj_func.phi_report(optimum_run.get_obs(), optimum_run.get_ctl_pars(), *(pest_scenario.get_regul_scheme_ptr()));
		output_file_writer.phi_report(fout_rec, termination_ctl.get_iteration_number() + 1, run_manager_ptr->get_total_runs(), phi_data, 0.0, true);
		output_file_writer.phi_report(cout, termination_ctl.get_iteration_number() + 1, run_manager_ptr->get_total_runs(), phi_data, 0.0, true);
		fout_rec << endl << endl;
		fout_rec << "Number of forward model runs performed during optimization: " << run_manager_ptr->get_total_runs() << endl;

		//linear analysis stuff
		if ((pest_scenario.get_control_info().noptmax != 0) &&
			(pest_scenario.get_pestpp_options().get_uncert_flag()))

		{
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


			performance_log.log_event("FOSM-based posterior unc calcs");

			if (base_jacobian_ptr->get_base_numeric_par_names().size() == 0)
			{
				cout << "WARNING: no parameters in base Jacobian, can't calculate uncertainty with FOSM" << endl;
				fout_rec << "WARNING: no parameters in base Jacobian, can't calculate uncertainty with FOSM" << endl;
				return 0;
			}

			//instance of a Mat for the jco
			Mat j(base_jacobian_ptr->get_sim_obs_names(), base_jacobian_ptr->get_base_numeric_par_names(),
				base_jacobian_ptr->get_matrix_ptr());
			if (pest_scenario.get_prior_info_ptr()->get_nnz_pi() > 0)
			{
				vector<string> pi_names = pest_scenario.get_ctl_ordered_pi_names();
				j.drop_rows(pi_names);
			}
			LinearAnalysis la(j, pest_scenario, file_manager, performance_log,parcov,&rand_gen);
			ObservationInfo reweight = la.glm_iter_fosm(optimum_run, output_file_writer, -999, run_manager_ptr);
			if (pest_scenario.get_pestpp_options().get_glm_num_reals() > 0)
			{
				cout << endl << "...drawing and running " << pest_scenario.get_pestpp_options().get_glm_num_reals() << " FOSM-based posterior realizations" << endl;
				
				pair<ParameterEnsemble, map<int, int>> fosm_real_info = la.draw_fosm_reals(run_manager_ptr, -999, optimum_run);
				run_manager_ptr->run();
				DynamicRegularization ptr;
				ptr.set_defaults();
				ptr.set_weight(0.0);
				double phi = optimum_run.get_phi(ptr);
				pair<ObservationEnsemble,map<string,double>> fosm_obs_info = la.process_fosm_reals(run_manager_ptr, fosm_real_info, -999, phi);
				
				//here is the adjustment process for each realization - one lambda each for now
				//todo: make sure to handle failed realizations - use oe real names to retrieve pe rows
				//todo: what about lambda?
				//QSqrtMatrix q(&reweight, pest_scenario.get_prior_info_ptr());
				//DynamicRegularization reg;
				//reg.set_zero();
				//Observations ctl_obs = pest_scenario.get_ctl_observations();
				//vector<string> nz_obs_names = pest_scenario.get_ctl_ordered_nz_obs_names();
				//Eigen::MatrixXd sim_reals = oe.get_eigen(oe.get_real_names(), nz_obs_names);
				//vector<string> active_par_names = base_jacobian_ptr->get_base_numeric_par_names();
				////Eigen::MatrixXd par_reals = fosm_real_info.first.get_eigen(fosm_real_info.first.get_real_names(), active_par_names);
				//Parameters freeze_pars,upgrade_pars,grad_upgrade_pars,del_upgrade_pars;
				//Eigen::MatrixXd upgraded_reals(oe.shape().first, active_par_names.size());
				//if (true)
				//{
				//	cout << "...making upgrade calculations for each realization..." << endl;
				//	vector<string> oe_real_names = oe.get_real_names();
				//	vector<double> lamb = pest_scenario.get_pestpp_options().get_base_lambda_vec();
				//	double min_lamb = lamb[0];
				//	for (int i = 0; i < oe.shape().first; i++)
				//	{
				//		
				//		Eigen::VectorXd res = ctl_obs.get_data_eigen_vec(nz_obs_names) - sim_reals.row(i).transpose();
				//		Parameters real_pars = pest_scenario.get_ctl_parameters();
				//		Eigen::VectorXd pe_vals = fosm_real_info.first.get_real_vector(oe_real_names[i]);
				//		vector<double> vals(pe_vals.data(), pe_vals.data() + pe_vals.size());
				//		real_pars.update(active_par_names, vals);
				//		base_trans_seq.numeric2active_ctl_ip(real_pars);
				//		base_svd.calc_lambda_upgrade_vec_JtQJ(*base_jacobian_ptr, q, reg, res, nz_obs_names,
				//			real_pars, freeze_pars, 0.0, upgrade_pars,
				//			del_upgrade_pars, grad_upgrade_pars);
				//		upgraded_reals.row(i) = upgrade_pars.get_data_eigen_vec(active_par_names);
				//	}

				//	ParameterEnsemble upgrade_pe(&pest_scenario, upgraded_reals, oe_real_names, active_par_names);
				//	map<int, int> run_id_map = upgrade_pe.add_runs(run_manager_ptr);
				//	run_manager_ptr->run();
				//	/*ObservationEnsemble oe_upgrade(&pest_scenario);
				//	oe_upgrade.reserve(oe_real_names, pest_scenario.get_ctl_ordered_obs_names());
				//	oe_upgrade.update_from_runs(run_id_map, run_manager_ptr);*/
				//	pair<ParameterEnsemble, map<int, int>> upgrade_fosm_real_info(upgrade_pe, run_id_map);
				//	ObservationEnsemble upgrade_oe = la.process_fosm_reals(run_manager_ptr, upgrade_fosm_real_info,-999, optimum_run.get_phi(ptr));
				//}
			}
		}

		// clean up
		//fout_rec.close();
		delete base_jacobian_ptr;
		delete super_jacobian_ptr;
		delete run_manager_ptr;

		string case_name = file_manager.get_base_filename();
		file_manager.close_file("rst");
		pest_utils::try_clean_up_run_storage_files(case_name);

		cout << endl << endl << "pestpp-glm analysis complete..." << endl;
        fout_rec << endl << endl << "pestpp-glm analysis complete..." << endl;
        auto end = chrono::steady_clock::now();
        cout << "started at " << start_string << endl;
        cout << "finished at " << get_time_string() << endl;
        cout << "took " << setprecision(6) << (double)chrono::duration_cast<chrono::seconds>(end - start).count()/60.0 << " minutes" << endl;
        fout_rec << "started at " << start_string << endl;
        fout_rec << "finished at " << get_time_string() << endl;
        fout_rec << "took " << setprecision(6) << (double)chrono::duration_cast<chrono::seconds>(end - start).count()/60.0 << " minutes" << endl;
        fout_rec.close();
        cout << flush;
		return 0;
#ifndef _DEBUG
	}
	catch (exception &e)
	{
		cout << "Error condition prevents further execution: " << endl << e.what() << endl;
		return 1;
		//cout << "press enter to continue" << endl;
		//char buf[256];
		//OperSys::gets_s(buf, sizeof(buf));
	}
	catch (...)
	{
		cout << "Error condition prevents further execution" << endl;
		return 1;
	}
#endif
}

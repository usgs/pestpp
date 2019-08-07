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
		cout << endl << endl << "version: " << version << endl;
		cout << "binary compiled on " << __DATE__ << " at " << __TIME__ << endl << endl;
		// build commandline
		string commandline = "";
		for (int i = 0; i < argc; ++i)
		{
			commandline.append(" ");
			commandline.append(argv[i]);
		}

		vector<string> cmd_arg_vec(argc);
		copy(argv, argv + argc, cmd_arg_vec.begin());
		for (vector<string>::iterator it = cmd_arg_vec.begin(); it != cmd_arg_vec.end(); ++it)
		{
			transform(it->begin(), it->end(), it->begin(), ::tolower);
		}

		string complete_path;
		enum class RunManagerType { SERIAL, PANTHER, GENIE, EXTERNAL };

		if (argc >= 2) {
			complete_path = argv[1];
		}
		else {
			cerr << "--------------------------------------------------------" << endl;
			cerr << "usage:" << endl << endl;
			cerr << "    serial run manager:" << endl;
			cerr << "        pest++ pest_ctl_file.pst" << endl << endl;
			cerr << "    PANTHER master:" << endl;
			cerr << "        pest++ control_file.pst /H :port" << endl << endl;
			cerr << "    PANTHER runner:" << endl;
			cerr << "        pest++ control_file.pst /H hostname:port " << endl << endl;
			cerr << "    external run manager:" << endl;
			cerr << "        pest++ control_file.pst /E" << endl << endl;
			cerr << " additional options can be found in the PEST++ manual" << endl;
			cerr << "--------------------------------------------------------" << endl;
			exit(0);
		}


		FileManager file_manager;
		string filename = complete_path;
		string pathname = ".";
		file_manager.initialize_path(get_filename_without_ext(filename), pathname);

		//by default use the serial run manager.  This will be changed later if another
		//run manger is specified on the command line.
		RunManagerType run_manager_type = RunManagerType::SERIAL;

		vector<string>::const_iterator it_find, it_find_next;
		string next_item;
		string socket_str = "";
		//Check for external run manager
		it_find = find(cmd_arg_vec.begin(), cmd_arg_vec.end(), "/e");
		if (it_find != cmd_arg_vec.end() )
		{
			run_manager_type = RunManagerType::EXTERNAL;
		}
		//Check for PANTHER Slave
		it_find = find(cmd_arg_vec.begin(), cmd_arg_vec.end(), "/h");
		next_item.clear();
		if (it_find != cmd_arg_vec.end() && it_find + 1 != cmd_arg_vec.end())
		{
			next_item = *(it_find + 1);
			strip_ip(next_item);
		}
		if (it_find != cmd_arg_vec.end() && !next_item.empty() && next_item[0] != ':')
		{
			// This is a PANTHER Slave, start PEST++ as a PANTHER Slave
			vector<string> sock_parts;
			vector<string>::const_iterator it_find_PANTHER_ctl;
			string file_ext = get_filename_ext(filename);
			tokenize(next_item, sock_parts, ":");
			try
			{
				if (sock_parts.size() != 2)
				{
					cerr << "PANTHER agent requires the master be specified as /H hostname:port" << endl << endl;
					throw(PestCommandlineError(commandline));
				}
				PANTHERAgent yam_agent;
				string ctl_file = "";
				try {
					string ctl_file;
					if (upper_cp(file_ext) == "YMR")
					{
						ctl_file = file_manager.build_filename("ymr");
						yam_agent.process_panther_ctl_file(ctl_file);
					}
					else
					{
						// process traditional PEST control file
						ctl_file = file_manager.build_filename("pst");
						yam_agent.process_panther_ctl_file(ctl_file);
					}
				}
				catch (PestError e)
				{
					cerr << "Error prococessing control file: " << ctl_file << endl << endl;
					cerr << e.what() << endl << endl;
					throw(e);
				}

				yam_agent.start(sock_parts[0], sock_parts[1]);
			}
			catch (PestError &perr)
			{
				cerr << perr.what();
				throw(perr);
			}
			cout << endl << "Simulation Complete..." << endl;
			exit(0);
		}
		//Check for PANTHER Master
		else if (it_find != cmd_arg_vec.end())
		{
			// using PANTHER run manager
			run_manager_type = RunManagerType::PANTHER;
			socket_str = next_item;
		}

		it_find = find(cmd_arg_vec.begin(), cmd_arg_vec.end(), "/g");
		next_item.clear();
		if (it_find != cmd_arg_vec.end())
		{
			cerr << "Genie run manager ('/g') no longer supported, please use PANTHER instead" << endl;
			return 1;

		}
		

		RestartController restart_ctl;

		//process restart and reuse jacobian directives
		vector<string>::const_iterator it_find_j = find(cmd_arg_vec.begin(), cmd_arg_vec.end(), "/j");
		vector<string>::const_iterator it_find_r = find(cmd_arg_vec.begin(), cmd_arg_vec.end(), "/r");
		bool restart_flag = false;
		bool save_restart_rec_header = true;

		debug_initialize(file_manager.build_filename("dbg"));
		if (it_find_j != cmd_arg_vec.end())
		{
			restart_ctl.get_restart_option() = RestartController::RestartOption::REUSE_JACOBIAN;
			file_manager.open_default_files();
		}
		else if (it_find_r != cmd_arg_vec.end())
		{
			ifstream &fin_rst = file_manager.open_ifile_ext("rst");
			restart_ctl.process_rst_file(fin_rst);
			file_manager.close_file("rst");
			restart_flag = true;
			file_manager.open_default_files(true);
			ofstream &fout_rec_tmp = file_manager.rec_ofstream();
			fout_rec_tmp << endl << endl;
			if (run_manager_type == RunManagerType::EXTERNAL)
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
		PerformanceLog performance_log(file_manager.open_ofile_ext("pfm"));

		if (!restart_flag || save_restart_rec_header)
		{
			fout_rec << "             pestpp-glm " << version << endl << endl;
			fout_rec << "    by The PEST++ Development Team" << endl;

			fout_rec << endl;
			fout_rec << endl << endl << "version: " << version << endl;
			fout_rec << "binary compiled on " << __DATE__ << " at " << __TIME__ << endl << endl;
			fout_rec << "using control file: \"" << complete_path << "\"" << endl << endl;
			fout_rec << "in directory: \"" << OperSys::getcwd() << "\"" << endl << endl;

		}

		cout << endl;
		cout << "using control file: \"" << complete_path << "\"" << endl << endl;
		cout << "in directory: \"" << OperSys::getcwd() << "\"" << endl << endl;

		// create pest run and process control file to initialize it
		Pest pest_scenario;
		pest_scenario.set_defaults();

		try {
			performance_log.log_event("starting to process control file", 1);
			pest_scenario.process_ctl_file(file_manager.open_ifile_ext("pst"), file_manager.build_filename("pst"),fout_rec);
			file_manager.close_file("pst");
			performance_log.log_event("finished processing control file");
		}
		catch (PestError e)
		{
			cerr << "Error prococessing control file: " << filename << endl << endl;
			cerr << e.what() << endl << endl;
			throw(e);
	 	}
		pest_scenario.check_inputs(fout_rec);

		//if base jco arg read from control file, reset restart controller
		if (!pest_scenario.get_pestpp_options().get_basejac_filename().empty())
		{
			restart_ctl.get_restart_option() = RestartController::RestartOption::REUSE_JACOBIAN;
		}

		//Initialize OutputFileWriter to hadle IO of suplementary files (.par, .par, .svd)
		//bool save_eign = pest_scenario.get_svd_info().eigwrite > 0;
		OutputFileWriter output_file_writer(file_manager, pest_scenario, restart_flag);
		output_file_writer.prep_glm_files(restart_flag);
		output_file_writer.set_svd_output_opt(pest_scenario.get_svd_info().eigwrite);
		if (!restart_flag)
		{
			output_file_writer.scenario_report(fout_rec);
		}
		if (pest_scenario.get_pestpp_options().get_iter_summary_flag())
		{
			output_file_writer.write_par_iter(0, pest_scenario.get_ctl_parameters());
		}
		RunManagerAbstract *run_manager_ptr;
		if (run_manager_type == RunManagerType::PANTHER)
		{
			// using PANTHER run manager
			if (pest_scenario.get_control_info().noptmax == 0)
			{
				cout << endl << endl << "WARNING: 'noptmax' = 0 but using parallel run mgr.  This prob isn't what you want to happen..." << endl << endl;
			}
			string port = socket_str;
			strip_ip(port);
			strip_ip(port, "front", ":");
			const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
			//check for condor wrapper
			string csf = pest_scenario.get_pestpp_options().get_condor_submit_file();
			if (csf.size() > 0)
			{
				if (!pest_utils::check_exist_in(csf))
					throw runtime_error("++condor_submit_file '" + csf + "' not found");
				run_manager_ptr = new RunManagerYAMRCondor(
					file_manager.build_filename("rns"), port,
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
					file_manager.build_filename("rns"), port,
					file_manager.open_ofile_ext("rmr"),
					pest_scenario.get_pestpp_options().get_max_run_fail(),
					pest_scenario.get_pestpp_options().get_overdue_reched_fac(),
					pest_scenario.get_pestpp_options().get_overdue_giveup_fac(),
					pest_scenario.get_pestpp_options().get_overdue_giveup_minutes());
			}
		}
		
		else if (run_manager_type == RunManagerType::EXTERNAL)
		{
			const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
			run_manager_ptr = new RunManagerExternal(exi.comline_vec,
				exi.tplfile_vec, exi.inpfile_vec, exi.insfile_vec, exi.outfile_vec,
				file_manager.build_filename("rns"), file_manager.build_filename("ext"),
				file_manager.build_filename("exi"),
				pest_scenario.get_pestpp_options().get_max_run_fail());
		}
		else
		{
			performance_log.log_event("starting basic model IO error checking", 1);
			cout << "checking model IO files...";
			pest_scenario.check_io();
			performance_log.log_event("finished basic model IO error checking");
			cout << "done" << endl;
			const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
			run_manager_ptr = new RunManagerSerial(exi.comline_vec,
				exi.tplfile_vec, exi.inpfile_vec, exi.insfile_vec, exi.outfile_vec,
				file_manager.build_filename("rns"), pathname,
				pest_scenario.get_pestpp_options().get_max_run_fail());
		}

		const ParamTransformSeq &base_trans_seq = pest_scenario.get_base_par_tran_seq();

		ObjectiveFunc obj_func(&(pest_scenario.get_ctl_observations()), &(pest_scenario.get_ctl_observation_info()), &(pest_scenario.get_prior_info()));
		Jacobian *base_jacobian_ptr = new Jacobian_1to1(file_manager,output_file_writer);

		TerminationController termination_ctl(pest_scenario.get_control_info().noptmax, pest_scenario.get_control_info().phiredstp,
			pest_scenario.get_control_info().nphistp, pest_scenario.get_control_info().nphinored, pest_scenario.get_control_info().relparstp,
			pest_scenario.get_control_info().nrelpar, pest_scenario.get_regul_scheme_ptr()->get_use_dynamic_reg(),
			pest_scenario.get_regul_scheme_ptr()->get_phimaccept());

		//if we are doing a restart, update the termination_ctl
		if (restart_flag)
		{
			restart_ctl.update_termination_ctl(termination_ctl);
		}

		//SVDSolver::MAT_INV mat_inv = SVDSolver::MAT_INV::JTQJ;
		SVDSolver base_svd(pest_scenario, file_manager, &obj_func, base_trans_seq,
			*base_jacobian_ptr, output_file_writer, &performance_log, "base parameter solution");

		base_svd.set_svd_package(pest_scenario.get_pestpp_options().get_svd_pack());
		//Build Super-Parameter problem
		Jacobian *super_jacobian_ptr = new Jacobian(file_manager);
		ParamTransformSeq trans_svda;
		// method must be involked as pointer as the transformation sequence it is added to will
		// take responsibility for destroying it
		TranSVD *tran_svd = new TranSVD(pest_scenario.get_pestpp_options().get_max_n_super(),
			pest_scenario.get_pestpp_options().get_super_eigthres(), "SVD Super Parameter Tranformation");

		if (pest_scenario.get_pestpp_options().get_svd_pack() == PestppOptions::PROPACK)
		{
			tran_svd->set_SVD_pack_propack();
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
		// if noptmax=0 make one run with the intital parameters
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
			int rand_seed = 1;
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
			exit(1);
		}


		//Define model Run for Base Parameters (uses base parameter tranformations)
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
			//base parameter iterations
			try
			{
				if (restart_ctl.get_restart_option() != RestartController::RestartOption::NONE  && restart_ctl.get_iteration_type() == RestartController::IterationType::SUPER)
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
						base_svd.iteration_reuse_jac(*run_manager_ptr, termination_ctl, cur_run, true, filename, res_filename);
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
							cout << "error filling base jacobian: " << e.what() << endl;
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
						cout << "error in base parameter iteration process with existing jacobian: " << e.what() << endl;
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
					(*tran_svd).update_reset_frozen_pars(*base_jacobian_ptr, Q_sqrt, base_numeric_pars, max_n_super, super_eigthres, pars, nonregul_obs, cur_run.get_frozen_ctl_pars());
					(*tr_svda_fixed).reset((*tran_svd).get_frozen_derivative_pars());
				}
				SVDASolver super_svd(pest_scenario, file_manager, &obj_func,
					trans_svda, *super_jacobian_ptr,
					output_file_writer, &performance_log,
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
					cur_run = super_svd.update_run(*run_manager_ptr, cur_run);
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
			}
		}
		cout << endl;
		termination_ctl.termination_summary(cout);
		cout << endl;
		termination_ctl.termination_summary(fout_rec);
		fout_rec << endl;
		cout << "FINAL OPTIMISATION RESULTS" << endl << endl;
		fout_rec << "FINAL OPTIMISATION RESULTS" << endl << endl;

		fout_rec << "  Optimal parameter values  " << endl;
		output_file_writer.par_report(fout_rec, optimum_run.get_ctl_pars());

		fout_rec << endl << "  Observations with optimal model-simulated equivalents and residuals" << endl;
		ObservationInfo oi = pest_scenario.get_ctl_observation_info();
		output_file_writer.obs_report(fout_rec, *obj_func.get_obs_ptr(), optimum_run.get_obs(),oi);

		fout_rec << endl << "Final composite objective function " << endl;
		PhiData phi_data = obj_func.phi_report(optimum_run.get_obs(), optimum_run.get_ctl_pars(), *(pest_scenario.get_regul_scheme_ptr()));
		output_file_writer.phi_report(fout_rec, termination_ctl.get_iteration_number() + 1, run_manager_ptr->get_total_runs(), phi_data, 0.0, true);
		output_file_writer.phi_report(cout, termination_ctl.get_iteration_number() + 1, run_manager_ptr->get_total_runs(), phi_data, 0.0, true);
		fout_rec << endl << endl;
		fout_rec << "Number of forward model runs performed during optimiztion: " << run_manager_ptr->get_total_runs() << endl;

		//linear analysis stuff
		if ((pest_scenario.get_control_info().noptmax != 0) &&
			(pest_scenario.get_pestpp_options().get_uncert_flag()))

		{
			cout << endl << endl << endl;
			cout << endl << endl << endl;
			cout << "  ---  starting uncertainty analysis calculations  ---  " << endl << endl << endl;
			cout << "  uncertainty estimates calculated using Schur's " << endl;
			cout << "  complement for linear-based conditional uncertainty " << endl;
			cout << "  propogation.  For a derviation from Bayes equation, see " << endl;
			cout << "  M. N. Fienen, J. E. Doherty, R. J. Hunt, and H. W. Reeves. " << endl;
			cout << "  2010. 'Using Prediction Uncertainty Analysis to Design " << endl;
			cout << "  Hydrologic Monitoring Networks : Example Applications " << endl;
			cout << "  from the Great Lakes Water Availability Pilot Project'. " << endl;
			cout << "  See PEST++ V3 documentation for implementation details." << endl;
			cout << endl << endl << endl;

			fout_rec << endl << endl << endl << endl;
			fout_rec << "-----------------------------------------------------------------------" << endl;
			fout_rec << "Note: The following uncertainty estimates were calculated using " << endl;
			fout_rec << "      Schur's complement for linear-based conditional uncertainty " << endl;
			fout_rec << "      propogation.  For a derviation from Bayes equation, see " << endl;
			fout_rec << "      M. N. Fienen, J. E. Doherty, R. J. Hunt, and H. W. Reeves. " << endl;
			fout_rec << "      2010. 'Using Prediction Uncertainty Analysis to Design " << endl;
			fout_rec << "      Hydrologic Monitoring Networks : Example Applications " << endl;
			fout_rec << "      from the Great Lakes Water Availability Pilot Project'. " << endl;
			fout_rec << "      See PEST++ V3 documentation for implementation details." << endl;
			fout_rec << "-----------------------------------------------------------------------" << endl;
			fout_rec << endl;

			fout_rec << "Note: Any observations or prior information equations with a group name" << endl;
			fout_rec << "      starting with 'regul' are dropped from the jacobian and observation" << endl;
			fout_rec << "      covariance matrices before uncertainty calculations.  Please" << endl;
			fout_rec << "      make sure that all expert knowledge is expressed in the prior " << endl;
			fout_rec << "      parameter bounds or through a covariance matix, which can be " << endl;
			fout_rec << "      supplied as a ++ option as 'parameter_covariance(<matrix_file_name>)," << endl;
			fout_rec << "      where <matrix_file_name> can be an ASCII PEST-compatible matrix file (.mat) or" << endl;
			fout_rec << "      a PEST-compatible uncertainty file (.unc)." << endl << endl;


			ofstream &pfm = file_manager.get_ofstream("pfm");
			pfm << endl << endl << "-----------------------------------" << endl;
			pfm << "starting linear uncertainty analyses" << endl;
			pfm << "-----------------------------------" << endl << endl;
			Logger unc_log(pfm);

			//instance of a Mat for the jco
			Mat j(base_jacobian_ptr->get_sim_obs_names(), base_jacobian_ptr->get_base_numeric_par_names(),
				base_jacobian_ptr->get_matrix_ptr());

			//get a new obs info instance that accounts for residual phi (and expected objection value if passed)
			// and report new weights to the rec file
			fout_rec << endl;
			ObservationInfo reweight;
			Observations sim = optimum_run.get_obs();
			reweight = normalize_weights_by_residual(pest_scenario, sim);
			fout_rec << "Note: The observation covariance matrix has been constructed from " << endl;
			fout_rec << "      weights listed in the pest control file that have been scaled by " << endl;
			fout_rec << "      by the final residuals to account for " << endl;
			fout_rec << "      the level of measurement noise implied by the original weights so" << endl;
			fout_rec << "      the total objective function is equal to the number of  " << endl;
			fout_rec << "      non-zero weighted observations." << endl;


			fout_rec << endl;
			/*fout_rec << "Scaled observation weights used to form observation noise covariance matrix written to residual file " <<  endl;
			fout_rec << endl << setw(20) << "observation" << setw(20) << "group" << setw(20) << "scaled_weight" << endl;
			for (auto &oi : reweight.observations)
			if (oi.second.weight > 0.0)
				fout_rec << setw(20) << oi.first << setw(20) << oi.second.group << setw(20) << oi.second.weight << endl;
			fout_rec << endl << endl;*/
			string reres_filename = file_manager.get_base_filename() + ".fosm_reweight.rei";
			ofstream &reres_of = file_manager.open_ofile_absolute("fosm_reweight.rei",reres_filename);
			
			Observations obs = pest_scenario.get_ctl_observations();
			output_file_writer.obs_report(reres_of, obs, sim, reweight);
			fout_rec << "Scaled observation weights used to form observation noise covariance matrix written to residual file '" << reres_filename << "'" << endl << endl;
			
			//instance of linear analysis
			linear_analysis la(j, pest_scenario, file_manager, &unc_log);

			//if needed, set the predictive sensitivity vectors
			const vector<string> pred_names = pest_scenario.get_pestpp_options().get_prediction_names();
			//make sure prediction weights are zero
			for (auto &pname : pred_names)
			{
				if (pest_scenario.get_ctl_observation_info().get_weight(pname) != 0.0)
				{
					cout << endl << "WARNING: prediction: " << pname << " has a non-zero weight" << endl << endl;
					fout_rec << endl << "WARNING: prediction: " << pname << " has a non-zero weight" << endl << endl;
				}
			}
			if (pred_names.size() > 0)
				la.set_predictions(pred_names);

			//drop all 'regul' obs and equations
			la.drop_prior_information(pest_scenario);

			//write the posterior covariance matrix
			string postcov_filename = file_manager.get_base_filename() + ".post.cov";
			la.posterior_parameter_ptr()->to_ascii(postcov_filename);
			fout_rec << "Note : posterior parameter covariance matrix written to file '" + postcov_filename +
				"'" << endl << endl;

			//write a parameter prior and posterior summary to the rec file
			const ParamTransformSeq trans = pest_scenario.get_base_par_tran_seq();
			Parameters pars = pest_scenario.get_ctl_parameters();
			string parsum_filename = file_manager.get_base_filename() + ".par.usum.csv";
			la.write_par_credible_range(fout_rec, parsum_filename, pest_scenario.get_ctl_parameter_info(),
				trans.active_ctl2numeric_cp(pest_scenario.get_ctl_parameters()),
				trans.active_ctl2numeric_cp(optimum_run.get_ctl_pars()),
				pest_scenario.get_ctl_ordered_par_names());
			fout_rec << "Note : the above parameter uncertainty summary was written to file '" + parsum_filename +
				"'" << endl << endl;
			//if predictions were defined, write a prior and posterior summary to the rec file
			if (pred_names.size() > 0)
			{
				map<string, pair<double, double>> init_final_pred_values;
				double ival, fval;
				for (auto &pred_name : pred_names)
				{
					fval = optimum_run.get_obs().get_rec(pred_name);
					if (run_manager_ptr->get_init_sim().size() > 0)
					{
						int idx = distance(run_manager_ptr->get_obs_name_vec().begin(), find(run_manager_ptr->get_obs_name_vec().begin(),
							run_manager_ptr->get_obs_name_vec().end(), pred_name));
						ival = run_manager_ptr->get_init_sim()[idx];
					}
					else
					{
						cout << "WARNING: initial simulation results not available, falling back to optimum run outputs for prior forecast mean" << endl;
						fout_rec << "WARNING: initial simulation results not available, falling back to optimum run outputs for prior forecast mean" << endl;
						ival = fval;
					}
					
					init_final_pred_values[pred_name] = pair<double, double>(ival, fval);
				}
				string predsum_filename = file_manager.get_base_filename() + ".pred.usum.csv";
				la.write_pred_credible_range(fout_rec, predsum_filename, init_final_pred_values);
				fout_rec << "Note : the above prediction uncertainty summary was written to file '" + predsum_filename +
					"'" << endl << endl;
			}
			set<string> args = pest_scenario.get_pestpp_options().get_passed_args();

			if (args.find("NUM_REALS") != args.end() && pest_scenario.get_pestpp_options().get_ies_num_reals() > 0)
			{
				bool binary = pest_scenario.get_pestpp_options().get_ies_save_binary();
				int num_reals = pest_scenario.get_pestpp_options().get_ies_num_reals();
				fout_rec << "drawing " << num_reals << " posterior parameter realizations";
				ParameterEnsemble pe(&pest_scenario);
				Covariance cov = la.posterior_parameter_matrix();
				pe.draw(num_reals,optimum_run.get_ctl_pars(),cov, &performance_log,1);
				if (binary)
					pe.to_binary(file_manager.get_base_filename() + ".post.paren.jcb");
				else
					pe.to_csv(file_manager.get_base_filename() + ".post.paren.csv");
				map<int,int> run_map = pe.add_runs(run_manager_ptr);
				run_manager_ptr->run();
				ObservationEnsemble oe(&pest_scenario);
				Covariance obscov = la.get_obscov();
				oe.draw(num_reals, obscov, &performance_log, 1);
				oe.update_from_runs(run_map, run_manager_ptr);
				if (binary)
					oe.to_binary(file_manager.get_base_filename() + ".post.obsen.jcb");
				else
					oe.to_csv(file_manager.get_base_filename() + ".post.obsen.csv");
			}
			cout << "  ---  finished uncertainty analysis calculations  ---  " << endl << endl << endl;
		}

		// clean up
		//fout_rec.close();
		delete base_jacobian_ptr;
		delete super_jacobian_ptr;
		delete run_manager_ptr;
		cout << endl << endl << "Simulation Complete..." << endl;
		cout << flush;
#ifndef _DEBUG
	}
	catch (exception &e)
	{
		cout << "Error condition prevents further execution: " << endl << e.what() << endl;
		//cout << "press enter to continue" << endl;
		//char buf[256];
		//OperSys::gets_s(buf, sizeof(buf));
	}
#endif
}

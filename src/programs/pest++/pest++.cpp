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
/**
 * @file pest++.cpp
 * @brief Implementation of pest++.
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
#include "GLM.h"
#include "SVDSolver.h"
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
#include "linear_analysis.h"
#include "logger.h"
#include "covariance.h"
#include "Ensemble.h"
#include "eigen_tools.h"


using namespace std;
using namespace pest_utils;

//using namespace pest_utils;

/**
 * @brief Main.
 *
 * @param argc Description.
 * @param argv Description.
 *
 * @return Description.
 */
int main(int argc, char* argv[]) {
#ifndef _DEBUG
	try {
#endif
		string version = PESTPP_VERSION;
		cout << endl << endl;
		cout << "             pestpp-glm: a tool for GLM parameter estimation and FOSM uncertainty analysis" << endl << endl;
		cout << "                                   by The PEST++ Development Team" << endl;
		cout << endl;
		auto start = chrono::steady_clock::now();

		string start_string = get_time_string();
		CmdLine cmdline(argc, argv);
		// -v/--version has already printed the version; nothing else is set up
		if (cmdline.version_requested)
			return 0;

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
		//if base jco arg read from control file, reset restart controller
		if (!pest_scenario.get_pestpp_options().get_basejac_filename().empty())
		{
			restart_ctl.get_restart_option() = RestartController::RestartOption::REUSE_JACOBIAN;
		}

		//the svd-assist options that used to be caught here are now refused by name when the
		//control file is parsed - see get_retired_message() - so there is nothing left to check

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
                    pest_scenario.get_pestpp_options().get_panther_persistent_workers(),
					pest_scenario.get_pestpp_options().get_panther_ping_interval_secs());
				// host failure screening is a master-side policy, so it is set here rather than threaded
				// through the constructor's already long parameter list
				dynamic_cast<RunManagerPanther*>(run_manager_ptr)->set_max_failed_run_delta(
				    pest_scenario.get_pestpp_options().get_panther_agent_max_failed_run_delta());
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

		run_manager_ptr->set_save_all_runs(pest_scenario.get_pestpp_options().get_save_all_runs());

		// pestpp-glm's loop lives in GLM, the peer of MOEA / SeqQuadProgram / EnsembleMethod.
		// What is left here is what the other tools' mains keep: arguments, the scenario, the
		// run manager, and the restart policy the command line decided.
		GLM glm(pest_scenario, file_manager, output_file_writer, &performance_log,
			run_manager_ptr, &restart_ctl, restart_flag, save_restart_rec_header);
		glm.initialize();
		glm.iterate_2_solution();
		glm.finalize();

		// clean up
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
        if (pest_scenario.get_pestpp_options().get_save_all_runs())
        {
            RunStorage::print_persistent_summary(case_name + ".allruns.rns", cout);
            RunStorage::print_persistent_summary(case_name + ".allruns.rns", fout_rec);
        }
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

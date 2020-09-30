// pestpp-mou.cpp : Defines the entry point for the console application.
//

#include "RunManagerPanther.h" //needs to be first because it includes winsock2.h
//#include <vld.h> // Memory Leak Detection using "Visual Leak Detector"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include "config_os.h"
#include "Pest.h"
#include "Transformable.h"
#include "Transformation.h"
#include "ParamTransformSeq.h"
#include "utilities.h"
#include "pest_error.h"
#include "ModelRunPP.h"
#include "FileManager.h"
#include "RunManagerSerial.h"
#include "OutputFileWriter.h"
#include "PantherAgent.h"
#include "Serialization.h"
#include "system_variables.h"
#include "pest_error.h"
#include "RestartController.h"
#include "PerformanceLog.h"
#include "debug.h"
#include "logger.h"
#include "Ensemble.h"
#include "MOEA.h"

using namespace std;
using namespace pest_utils;


int main(int argc, char* argv[])
{
#ifndef _DEBUG
	try
	{
#endif
		string version = PESTPP_VERSION;
		cout << endl << endl;
		cout << "             pestpp-mou: multi-objective optimization with uncertainty for PEST++ datasets" << endl << endl;
		//cout << "                     for PEST(++) datasets " << endl << endl;
		cout << "                   by the PEST++ development team" << endl;
		cout << endl << endl << "version: " << version << endl;
		cout << "binary compiled on " << __DATE__ << " at " << __TIME__ << endl << endl;

		CmdLine cmdline(argc, argv);

		FileManager file_manager;
		string filename = cmdline.ctl_file_name;

		string pathname = ".";
		file_manager.initialize_path(get_filename_without_ext(filename), pathname);
		//jwhite - something weird is happening with the machine is busy and an existing
		//rns file is really large. so let's remove it explicitly and wait a few seconds before continuing...
		string rns_file = file_manager.build_filename("rns");
		int flag = remove(rns_file.c_str());

		
		if (cmdline.runmanagertype == CmdLine::RunManagerType::EXTERNAL)
		{
			cerr << "External run manager ('/e') not supported by pestpp-mou, please use panther instead" << endl;
			exit(1);
		}
		if (cmdline.runmanagertype == CmdLine::RunManagerType::GENIE)
		{
			cerr << "genie run manager ('/g') not supported by pestpp-mou, please use panther instead" << endl;
			exit(1);
		}
		
		if (cmdline.runmanagertype == CmdLine::RunManagerType::PANTHER_WORKER)
		{
			// This is a PANTHER worker, start PEST++ as a PANTHER worker
			
			try
			{

				ofstream frec("panther_worker.rec");
				if (frec.bad())
					throw runtime_error("error opening 'panther_worker.rec'");
				PANTHERAgent yam_agent(frec);
				string ctl_file = "";
				try {
					
					// process traditional PEST control file
					ctl_file = file_manager.build_filename("pst");
					yam_agent.process_ctl_file(ctl_file);
					
				}
				catch (PestError e)
				{
					cerr << "Error processing control file: " << ctl_file << endl << endl;
					cerr << e.what() << endl << endl;
					throw(e);
				}
				catch (...)
				{
					cerr << "Error processing control file" << endl;
					throw runtime_error("error processing control file");
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

		RestartController restart_ctl;

		//process restart and reuse jacobian directives

		bool restart_flag = false;
		bool save_restart_rec_header = true;

		debug_initialize(file_manager.build_filename("dbg"));

		if (cmdline.jac_restart)
		{
			throw runtime_error("/j option not supported by pestpp-mou");
		}
		else if (cmdline.restart)

		{
			throw runtime_error("/r option not supported by pestpp-mou");
		}
		else
		{
			restart_ctl.get_restart_option() = RestartController::RestartOption::NONE;
			file_manager.open_default_files();
		}

		ofstream &fout_rec = file_manager.rec_ofstream();
		PerformanceLog performance_log(file_manager.open_ofile_ext("log"));

		if (!restart_flag || save_restart_rec_header)
		{
			fout_rec << "              pestpp-mou: multi-objective optimization with uncertainty" << endl;
			fout_rec << "                         by the PEST++ developement team" << endl << endl << endl;
			fout_rec << endl;
			fout_rec << endl << endl << "version: " << version << endl;
			fout_rec << "binary compiled on " << __DATE__ << " at " << __TIME__ << endl << endl;

			fout_rec << "using control file: \"" << cmdline.ctl_file_name << "\"" << endl;

			fout_rec << "in directory: \"" << OperSys::getcwd() << "\"" << endl << endl;
		}

		cout << endl;

		cout << "using control file: \"" << cmdline.ctl_file_name << "\"" << endl;

		cout << "in directory: \"" << OperSys::getcwd() << "\"" << endl << endl;

		// create pest run and process control file to initialize it
		Pest pest_scenario;
		pest_scenario.set_defaults();
		try {
			performance_log.log_event("starting to process control file");
			pest_scenario.process_ctl_file(file_manager.open_ifile_ext("pst"), file_manager.build_filename("pst"),fout_rec);
			file_manager.close_file("pst");
			performance_log.log_event("finished processing control file");
		}
		catch (PestError e)
		{
			cerr << "Error prococessing control file: " << filename << endl << endl;
			cerr << e.what() << endl << endl;
			fout_rec << "Error prococessing control file: " << filename << endl << endl;
			fout_rec << e.what() << endl << endl;
			fout_rec.close();
			throw(e);
		}
		pest_scenario.check_inputs(fout_rec);
		
		//Initialize OutputFileWriter to handle IO of suplementary files (.par, .par, .svd)
		//bool save_eign = pest_scenario.get_svd_info().eigwrite > 0;
		pest_scenario.get_pestpp_options_ptr()->set_iter_summary_flag(false);
		OutputFileWriter output_file_writer(file_manager, pest_scenario, restart_flag);
		output_file_writer.scenario_report(fout_rec,false);
		//output_file_writer.scenario_io_report(fout_rec);
		//ZAK: define mou scenario 
		if (pest_scenario.get_pestpp_options().get_ies_verbose_level() > 1)
		{
			output_file_writer.scenario_pargroup_report(fout_rec);
			output_file_writer.scenario_par_report(fout_rec);
			output_file_writer.scenario_obs_report(fout_rec);
		}

		//reset some default args here:
		/*PestppOptions *ppo = pest_scenario.get_pestpp_options_ptr();
		set<string> pp_args = ppo->get_passed_args();
		if (pp_args.find("MAX_RUN_FAIL") == pp_args.end())
			ppo->set_max_run_fail(1);
		if (pp_args.find("OVERDUE_GIVEUP_FAC") == pp_args.end())
			ppo->set_overdue_giveup_fac(2.0);
		if (pp_args.find("OVERDUE_resched_FAC") == pp_args.end())
			ppo->set_overdue_reched_fac(1.15);
		*/
		if (pest_scenario.get_pestpp_options().get_debug_parse_only())
		{
			cout << endl << endl << "DEBUG_PARSE_ONLY is true, exiting..." << endl << endl;
			exit(0);
		}

		RunManagerAbstract *run_manager_ptr;



		if (cmdline.runmanagertype == CmdLine::RunManagerType::PANTHER_MASTER)

		{
			if (pest_scenario.get_control_info().noptmax == 0)
			{
				cout << endl << endl << "WARNING: 'noptmax' = 0 but using parallel run mgr.  This prob isn't what you want to happen..." << endl << endl;
			}

			const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
			run_manager_ptr = new RunManagerPanther(
				rns_file, cmdline.panther_port,

				file_manager.open_ofile_ext("rmr"),
				pest_scenario.get_pestpp_options().get_max_run_fail(),
				pest_scenario.get_pestpp_options().get_overdue_reched_fac(),
				pest_scenario.get_pestpp_options().get_overdue_giveup_fac(),
				pest_scenario.get_pestpp_options().get_overdue_giveup_minutes());
		}
		else
		{
			performance_log.log_event("starting basic model IO error checking");
			cout << "checking model IO files...";
			pest_scenario.check_io(fout_rec);
			//pest_scenario.check_par_obs();
			performance_log.log_event("finished basic model IO error checking");
			cout << "done" << endl;
			const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
			run_manager_ptr = new RunManagerSerial(exi.comline_vec,
				exi.tplfile_vec, exi.inpfile_vec, exi.insfile_vec, exi.outfile_vec,
				file_manager.build_filename("rns"), pathname,
				pest_scenario.get_pestpp_options().get_max_run_fail(),
				pest_scenario.get_pestpp_options().get_fill_tpl_zeros(),
				pest_scenario.get_pestpp_options().get_additional_ins_delimiters());
		}


		const ParamTransformSeq &base_trans_seq = pest_scenario.get_base_par_tran_seq();
		ObjectiveFunc obj_func(&(pest_scenario.get_ctl_observations()), &(pest_scenario.get_ctl_observation_info()), &(pest_scenario.get_prior_info()));

		Parameters cur_ctl_parameters = pest_scenario.get_ctl_parameters();
		//Allocates Space for Run Manager.  This initializes the model parameter names and observations names.
		//Neither of these will change over the course of the simulation


		run_manager_ptr->initialize(base_trans_seq.ctl2model_cp(cur_ctl_parameters), pest_scenario.get_ctl_observations());
		
		MOEA moea(pest_scenario, file_manager,output_file_writer, &performance_log, run_manager_ptr);
		
		//ZAK: Initialize random generator here
		moea.initialize();
		moea.iterate_to_solution();
		moea.finalize();


		// clean up
		fout_rec.close();
		delete run_manager_ptr;
		cout << endl << endl << "pestpp-mou analysis complete..." << endl;
		cout << flush;
		return 0;
#ifndef _DEBUG
	}
	catch (exception &e)
	{
		cout << "Error condition prevents further execution: " << endl << e.what() << endl;
		//cout << "press enter to continue" << endl;
		//char buf[256];
		//OperSys::gets_s(buf, sizeof(buf));
		return 1;
	}
	catch (...)
	{
		cout << "Error condition prevents further execution: " << endl;
		return 1;
	}
#endif
}

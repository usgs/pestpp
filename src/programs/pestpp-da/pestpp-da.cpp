// pestpp-da.cpp : Defines the entry point for the console application.
//

#include "RunManagerPanther.h" //needs to be first because it includes winsock2.h
//#include <vld.h> // Memory Leak Detection using "Visual Leak Detector"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include "config_os.h"
#include "Pest.h"
#include "MultiPest.h"
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
#include "EnsembleSmoother.h"
#include "DataAssimilator.h"


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
		cout << "             ==============================================" << endl;
		cout << "             pestpp-da: Model Independent Data Assimilation" << endl;
		cout << "             ==============================================" << endl << endl;

		//cout << "                     for PEST(++) datasets " << endl << endl;
		cout << "               Developed by the PEST++ development team" << endl;
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
			cerr << "        pestpp-da control_file.pst" << endl << endl;
			cerr << "    PANTHER master:" << endl;
			cerr << "        pestpp-da control_file.pst /H :port" << endl << endl;
			cerr << "    PANTHER worker:" << endl;
			cerr << "        pestpp-da control_file.pst /H hostname:port " << endl << endl;
			
			cerr << " additional options can be found in the PEST++ manual" << endl;
			cerr << "--------------------------------------------------------" << endl;
			exit(0);
		}


		FileManager file_manager;
		string filename = complete_path;
		string pathname = ".";
		file_manager.initialize_path(get_filename_without_ext(filename), pathname);
		//jwhite - something weird is happening with the machine is busy and an existing
		//rns file is really large. so let's remove it explicitly and wait a few seconds before continuing...
		string rns_file = file_manager.build_filename("rns");
		int flag = remove(rns_file.c_str());
		//w_sleep(2000);
		//by default use the serial run manager.  This will be changed later if another
		//run manger is specified on the command line.
		RunManagerType run_manager_type = RunManagerType::SERIAL;

		vector<string>::const_iterator it_find, it_find_next;
		string next_item;
		string socket_str = "";
		//Check for external run manager
		it_find = find(cmd_arg_vec.begin(), cmd_arg_vec.end(), "/e");
		if (it_find != cmd_arg_vec.end())
		{
			throw runtime_error("External run manager not supported by pestpp-ies");
		}
		//Check for PANTHER worker
		it_find = find(cmd_arg_vec.begin(), cmd_arg_vec.end(), "/h");
		next_item.clear();
		if (it_find != cmd_arg_vec.end() && it_find + 1 != cmd_arg_vec.end())
		{
			next_item = *(it_find + 1);
			strip_ip(next_item);
		}
		if (it_find != cmd_arg_vec.end() && !next_item.empty() && next_item[0] != ':')
		{
			// This is a PANTHER worker, start PEST++ as a PANTHER worker
			vector<string> sock_parts;
			vector<string>::const_iterator it_find_yamr_ctl;
			string file_ext = get_filename_ext(filename);
			tokenize(next_item, sock_parts, ":");
			try
			{
				if (sock_parts.size() != 2)
				{
					cerr << "PANTHER worker requires the master be specified as /H hostname:port" << endl << endl;
					throw(PestCommandlineError(commandline));
				}
				ofstream frec("panther_agent.rec");
				if (frec.bad())
					throw runtime_error("error opening 'panther_agent.rec'");
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
				yam_agent.start(sock_parts[0], sock_parts[1]);
			}
			catch (PestError &perr)
			{
				cerr << perr.what();
				throw(perr);
			}
			
			cout << endl << "Work Done..." << endl;
			exit(0);
		}
		//Check for PANTHER master
		else if (it_find != cmd_arg_vec.end())
		{
			// using PANTHER run manager
			run_manager_type = RunManagerType::PANTHER;
			socket_str = next_item;
		}

		it_find = find(cmd_arg_vec.begin(), cmd_arg_vec.end(), "/g");
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
			throw runtime_error("/j option not supported by pestpp-ies");
		}
		else if (it_find_r != cmd_arg_vec.end())
		{
			throw runtime_error("/r option not supported by pestpp-ies");
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
			fout_rec << "             pestpp-da.exe - a Data Assimilation Tool" << endl << "for PEST(++) datasets " << endl << endl;
			fout_rec << "                 by the PEST++ developement team" << endl << endl << endl;
			fout_rec << endl;
			fout_rec << endl << endl << "version: " << version << endl;
			fout_rec << "binary compiled on " << __DATE__ << " at " << __TIME__ << endl << endl;
			fout_rec << "using control file: \"" << complete_path << "\"" << endl;
			fout_rec << "in directory: \"" << OperSys::getcwd() << "\"" << endl << endl;
		}

		cout << endl;
		cout << "using control file: \"" << complete_path << "\"" << endl;
		cout << "in directory: \"" << OperSys::getcwd() << "\"" << endl << endl;

		// create pest run and process control file to initialize it
		Pest pest_scenario;
		pest_scenario.set_defaults();
		try {
			performance_log.log_event("starting to process control file");
			pest_scenario.process_ctl_file(file_manager.open_ifile_ext("pst"), file_manager.build_filename("pst"),fout_rec);
			file_manager.close_file("pst");
			pest_scenario.assign_da_cycles(fout_rec);
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
		//output_file_writer.scenario_report(fout_rec);
		output_file_writer.scenario_io_report(fout_rec);
		if (pest_scenario.get_pestpp_options().get_ies_verbose_level() > 1)
		{
			output_file_writer.scenario_pargroup_report(fout_rec);
			output_file_writer.scenario_par_report(fout_rec);
			output_file_writer.scenario_obs_report(fout_rec);
		}
		
		

		//reset some default args for da here:
		PestppOptions *ppo = pest_scenario.get_pestpp_options_ptr();
		set<string> pp_args = ppo->get_passed_args();
		if (pp_args.find("MAX_RUN_FAIL") == pp_args.end())
			ppo->set_max_run_fail(1);
		if (pp_args.find("OVERDUE_GIVEUP_FAC") == pp_args.end())
			ppo->set_overdue_giveup_fac(2.0);
		if (pp_args.find("OVERDUE_resched_FAC") == pp_args.end())
			ppo->set_overdue_reched_fac(1.15);
		
		if (pest_scenario.get_pestpp_options().get_debug_parse_only())
		{
			cout << endl << endl << "DEBUG_PARSE_ONLY is true, exiting..." << endl << endl;
			exit(0);
		}

		RunManagerAbstract* run_manager_ptr;

		if (run_manager_type == RunManagerType::PANTHER)
		{
			if (pest_scenario.get_control_info().noptmax == 0)
			{
				cout << endl << endl << "WARNING: 'noptmax' = 0 but using parallel run mgr.  This prob isn't what you want to happen..." << endl << endl;
			}
			string port = socket_str;
			strip_ip(port);
			strip_ip(port, "front", ":");
			run_manager_ptr = new RunManagerPanther(
				rns_file, port,
				file_manager.open_ofile_ext("rmr"),
				pest_scenario.get_pestpp_options().get_max_run_fail(),
				pest_scenario.get_pestpp_options().get_overdue_reched_fac(),
				pest_scenario.get_pestpp_options().get_overdue_giveup_fac(),
				pest_scenario.get_pestpp_options().get_overdue_giveup_minutes(),
				pest_scenario.get_ctl_ordered_par_names(),
				pest_scenario.get_ctl_ordered_obs_names());
			run_manager_ptr->initialize(pest_scenario.get_ctl_parameters(), pest_scenario.get_ctl_observations());
		}
		else
		{
			performance_log.log_event("starting basic model IO error checking");
			cout << "checking model IO files...";
			pest_scenario.check_io(fout_rec);
			//pest_scenario.check_par_obs();
			performance_log.log_event("finished basic model IO error checking");
			cout << "done" << endl;
			const ModelExecInfo& exi = pest_scenario.get_model_exec_info();
			run_manager_ptr = new RunManagerSerial(exi.comline_vec,
				exi.tplfile_vec, exi.inpfile_vec, exi.insfile_vec, exi.outfile_vec,
				file_manager.build_filename("rns"), pathname,
				pest_scenario.get_pestpp_options().get_max_run_fail(),
				pest_scenario.get_pestpp_options().get_fill_tpl_zeros(),
				pest_scenario.get_pestpp_options().get_additional_ins_delimiters());
		}



		//process da par cycle table
		map<int, map<string, double>> par_cycle_info = process_da_par_cycle_table(pest_scenario, fout_rec);
		// process da obs cycle table
		set<string> obs_in_tbl; //we need this so we can set weights to zero in childpest if a value isnt listed for a given cycle
		map<int, map<string, double>> obs_cycle_info = process_da_obs_cycle_table(pest_scenario, fout_rec, obs_in_tbl);
		// loop over assimilation cycles

		
		vector<int> assimilation_cycles;
		pest_scenario.assign_da_cycles(fout_rec); 
		assimilation_cycles = pest_scenario.get_assim_cycles(fout_rec);
		ParameterEnsemble curr_pe;

		// loop over assimilation cycles
		for (auto icycle = assimilation_cycles.begin(); icycle != assimilation_cycles.end(); icycle++)
		{
			cout << endl;
			cout << " =======================================" << endl;
			cout << " >>>> Assimilating data in cycle " << *icycle << endl;
			cout << " =======================================" << endl;

			fout_rec << endl;
			fout_rec << " =======================================" << endl;
			fout_rec << " >>>> Assimilating data in cycle " << *icycle << endl;
			fout_rec << " =======================================" << endl;

			Pest childPest;
			childPest = pest_scenario.get_child_pest(*icycle);
			vector <string> xxxx=childPest.get_ctl_ordered_par_names();
			//childPest.get_pestpp_options.set_check_tplins(false);

			// -----------------------------  
			OutputFileWriter output_file_writer(file_manager, childPest, restart_flag);			
			output_file_writer.scenario_io_report(fout_rec);
			if (pest_scenario.get_pestpp_options().get_ies_verbose_level() > 1)
			{
				output_file_writer.scenario_pargroup_report(fout_rec);
				output_file_writer.scenario_par_report(fout_rec);
				output_file_writer.scenario_obs_report(fout_rec);
			}
			//------------------------------

			if (run_manager_type == RunManagerType::PANTHER)
			{
				run_manager_ptr->initialize(childPest.get_ctl_parameters(), childPest.get_ctl_observations());
			}
			else
			{
				performance_log.log_event("starting basic model IO error checking");
				cout << "checking model IO files...";
				childPest.check_io(fout_rec);
				//pest_scenario.check_par_obs();
				performance_log.log_event("finished basic model IO error checking");
				cout << "done" << endl;
				const ModelExecInfo& exi = childPest.get_model_exec_info();
				run_manager_ptr = new RunManagerSerial(exi.comline_vec,
					exi.tplfile_vec, exi.inpfile_vec, exi.insfile_vec, exi.outfile_vec,
					file_manager.build_filename("rns"), pathname,
					childPest.get_pestpp_options().get_max_run_fail(),
					childPest.get_pestpp_options().get_fill_tpl_zeros(),
					childPest.get_pestpp_options().get_additional_ins_delimiters());
			}



			ParamTransformSeq& base_trans_seq = childPest.get_base_par_tran_seq_4_mod();
			
			//check for entries in the par cycle table
			if (par_cycle_info.find(*icycle) != par_cycle_info.end())
			{
				Parameters pars = childPest.get_ctl_parameters_4_mod();
				TranFixed* fixed_ptr = base_trans_seq.get_fixed_ptr_4_mod();
				map<string, double> cycle_map = par_cycle_info[*icycle];
				for (auto item : cycle_map)
				{
					fixed_ptr->insert(item.first, item.second);
					pars.update_rec(item.first, item.second);
				}
			}
			//check for entries in the obs cycle table
			if (obs_cycle_info.find(*icycle) != obs_cycle_info.end())
			{
				Observations obs = childPest.get_ctl_observations();
				ObservationInfo* oi_ptr = childPest.get_observation_info_ptr();
				map<string, double> cycle_map = obs_cycle_info[*icycle];
				for (auto tbl_obs_name : obs_in_tbl)
				{
					if (cycle_map.find(tbl_obs_name) == cycle_map.end())
					{
						oi_ptr->set_weight(tbl_obs_name,0.0);
					}
					else
					{
						obs.update_rec(tbl_obs_name, cycle_map[tbl_obs_name]);
					}
				}
			}

			Parameters par1 = childPest.get_ctl_parameters();
			base_trans_seq.ctl2numeric_ip(par1);
			base_trans_seq.numeric2model_ip(par1);

			ObjectiveFunc obj_func(&(childPest.get_ctl_observations()), &(childPest.get_ctl_observation_info()), &(childPest.get_prior_info()));

			Parameters cur_ctl_parameters = childPest.get_ctl_parameters();
			//Allocates Space for Run Manager.  This initializes the model parameter names and observations names.
			//Neither of these will change over the course of the simulation
			//make sure we use the vector-based initializer so that the pars and obs are in order on the 
			//workers - PantherAgent uses this same strategy (child pest with cycle number, then sorted par and 
			//obs names)
			vector<string> par_names = base_trans_seq.ctl2model_cp(cur_ctl_parameters).get_keys();
			sort(par_names.begin(), par_names.end());
			vector<string> obs_names = childPest.get_ctl_observations().get_keys();
			sort(obs_names.begin(), obs_names.end());
			run_manager_ptr->initialize(par_names,obs_names);

			DataAssimilator da(childPest, file_manager, output_file_writer, &performance_log, run_manager_ptr);
			// use ies or da?
			da.use_ies = pest_scenario.get_pestpp_options_ptr()->get_da_use_ies();
			if (*icycle > 0)
			{
				da.set_pe(curr_pe);
			}

			
			da.initialize(*icycle);

			if (da.use_ies) // use ies
			{
				da.iterate_2_solution();
				curr_pe = da.get_pe();
				curr_pe.to_csv("cncnc.csv");
			}
			else // use da
			{
				da.kf_upate();
				curr_pe = da.get_pe();
				curr_pe.to_csv("cncnc.csv");

			}


					
		} // end cycle loop
		fout_rec.close();
		cout << endl << endl << "pestpp-da analysis complete..." << endl;
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
	}
#endif
}

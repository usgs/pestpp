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
        auto start = chrono::steady_clock::now();
        string start_string = get_time_string();
        cout << "started at " << start_string << endl;
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
			cerr << "External run manager ('/e') not supported by pestpp-da, please use panther instead" << endl;
			exit(1);
		}
		if (cmdline.runmanagertype == CmdLine::RunManagerType::GENIE)
		{
			cerr << "genie run manager ('/g') not supported by pestpp-da, please use panther instead" << endl;
			exit(1);
		}

		if (cmdline.runmanagertype == CmdLine::RunManagerType::PANTHER_WORKER)
		{
			// This is a PANTHER worker, start PEST++ as a PANTHER worker

			try
			{

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
				yam_agent.start(cmdline.panther_host_name, cmdline.panther_port);
			}
			catch (PestError & perr)
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
			throw runtime_error("/j option not supported by pestpp-da");
		}
		else if (cmdline.restart)
		{
			throw runtime_error("/r option not supported by pestpp-da");
		}
		else
		{
			restart_ctl.get_restart_option() = RestartController::RestartOption::NONE;
			file_manager.open_default_files();
		}

		ofstream& fout_rec = file_manager.rec_ofstream();
		PerformanceLog performance_log(file_manager.open_ofile_ext("pfm"));

		if (!restart_flag || save_restart_rec_header)
		{
			fout_rec << "             pestpp-da.exe - a Data Assimilation Tool" << endl << "for PEST(++) datasets " << endl << endl;
			fout_rec << "                 by the PEST++ developement team" << endl << endl << endl;
			fout_rec << endl;
			fout_rec << endl << endl << "version: " << version << endl;
			fout_rec << "binary compiled on " << __DATE__ << " at " << __TIME__ << endl << endl;
			fout_rec << "using control file: \"" << cmdline.ctl_file_name << "\"" << endl;
			fout_rec << "in directory: \"" << OperSys::getcwd() << "\"" << endl;
			fout_rec << "on host: \"" << w_get_hostname() << "\"" << endl;
            fout_rec << "started at " << start_string << endl << endl;
		}

		cout << endl;
		cout << "using control file: \"" << cmdline.ctl_file_name << "\"" << endl;
		cout << "in directory: \"" << OperSys::getcwd() << "\"" << endl;
		cout << "on host: \"" << w_get_hostname() << "\"" << endl << endl;

		// create pest run and process control file to initialize it
		Pest pest_scenario;
		try {
			performance_log.log_event("starting to process control file");
			pest_scenario.process_ctl_file(file_manager.open_ifile_ext("pst"), file_manager.build_filename("pst"), fout_rec);
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
		//pest_scenario.get_pestpp_options_ptr()->set_use_da(true);
		OutputFileWriter output_file_writer(file_manager, pest_scenario, restart_flag);
		//output_file_writer.scenario_report(fout_rec);
		fout_rec << endl << endl << "...Global Interface Summary:" << endl;
		output_file_writer.scenario_io_report(fout_rec);
		int verbose_level;
		PestppOptions* ppo = pest_scenario.get_pestpp_options_ptr();
		verbose_level = ppo->get_ies_verbose_level();
		output_file_writer.scenario_pargroup_report(fout_rec);
        output_file_writer.scenario_par_report(fout_rec);
        output_file_writer.scenario_obs_report(fout_rec);

		//reset some default args for da here:		
		set<string> pp_args = ppo->get_passed_args();
		if (pp_args.find("MAX_RUN_FAIL") == pp_args.end())
			ppo->set_max_run_fail(1);
		if (pp_args.find("OVERDUE_GIVEUP_FAC") == pp_args.end())
			ppo->set_overdue_giveup_fac(2.0);
		if (pp_args.find("OVERDUE_resched_FAC") == pp_args.end())
			ppo->set_overdue_reched_fac(1.15);

		pest_scenario.get_pestpp_options().summary(fout_rec);

		//do all this up here so we can use the parse only option
		//to check the interface and then quit
		//process da (recurrent??) par cycle table
		vector <int> cycles_in_tables;
		map<int, map<string, double>> par_cycle_info = process_da_par_cycle_table(pest_scenario, cycles_in_tables, fout_rec);
		// process da obs cycle table
		set<string> obs_in_tbl; //we need this so we can set weights to zero in childpest of a value isnt listed for a given cycle
		map<int, map<string, double>> obs_cycle_info = process_da_obs_cycle_table(pest_scenario, cycles_in_tables, fout_rec, obs_in_tbl);
		//process weights table
		set<string> weights_in_tbl;
		map<int, map<string, double>> weight_cycle_info = process_da_weight_cycle_table(pest_scenario, 
			cycles_in_tables, fout_rec, weights_in_tbl);
		vector<int> assimilation_cycles;
		//pest_scenario.get_pestpp_options
		int n_da_cycles = -1;// ppo->da_ctl_params.get_ivalue("DA_CYCLES_NUMBER");
		if (n_da_cycles > 0)
		{
			for (int i = 0; i < n_da_cycles; i++)
			{
				assimilation_cycles.push_back(i);
			}
		}
		else
		{
			//assimilation_cycles = pest_scenario.get_assim_cycles(fout_rec, cycles_in_tables);
			assimilation_cycles = pest_scenario.get_assim_dci_cycles(fout_rec,cycles_in_tables);
		}

		std::sort(assimilation_cycles.begin(), assimilation_cycles.end());
        stringstream ss;

		ss.str("");
		ss << "...assimilating over " << assimilation_cycles.size() << " cycles from " << assimilation_cycles[0] << " to " << assimilation_cycles[assimilation_cycles.size()-1] << endl;
		cout << ss.str();
		fout_rec << ss.str();

		int start_cycle = pest_scenario.get_pestpp_options().get_da_hotstart_cycle();
		int max_cycle = assimilation_cycles[assimilation_cycles.size() - 1];
		if (start_cycle > max_cycle)
		{
			ss.str("");
			ss << "'da_hotstart_cycle' (" << start_cycle << ") greater than max cycle (" << max_cycle << ")";
			throw runtime_error(ss.str());
		}

		if (pest_scenario.get_control_info().noptmax != 0)
		{

			for (auto icycle = assimilation_cycles.begin(); icycle != assimilation_cycles.end(); icycle++)
			{
				cout << endl;


				cout << " >>>> Checking data in cycle " << *icycle << endl;
				fout_rec << endl;
				fout_rec << " >>>> Checking data in cycle " << *icycle << endl;


				performance_log.log_event("instantiating child pest object");

				//Pest childPest;
				Pest childPest = pest_scenario.get_child_pest(*icycle);
				childPest.check_io(fout_rec);

				if (obs_cycle_info.find(*icycle) != obs_cycle_info.end())
				{
					performance_log.log_event("updating obs using da obs cycle table info");
					map<string, double> cycle_map = obs_cycle_info[*icycle];
					map<string, double> weight_cycle_map;
					if (weight_cycle_info.find(*icycle) != weight_cycle_info.end())
						weight_cycle_map = weight_cycle_info[*icycle];
					ObservationInfo oi = pest_scenario.get_ctl_observation_info();
					childPest.set_observation_info(oi);
					for (auto tbl_obs_name : obs_in_tbl)
					{
						//if this obs is not in this cycle, give it a zero weight
						if (cycle_map.find(tbl_obs_name) == cycle_map.end())
						{
							oi.set_weight(tbl_obs_name, 0.0);
						}
						else
						{
							childPest.get_ctl_observations_4_mod().update_rec(tbl_obs_name, cycle_map[tbl_obs_name]);
							//check if this obs is in this cycle's weight info
							if (weight_cycle_map.find(tbl_obs_name) != weight_cycle_map.end())
							{
								oi.set_weight(tbl_obs_name, weight_cycle_map[tbl_obs_name]);
								//pest_scenario.get_observation_info_ptr()->set_weight(tbl_obs_name, weight_cycle_map[tbl_obs_name]);
							}
						}
					}
					childPest.set_observation_info(oi);

				}

				Parameters par1 = childPest.get_ctl_parameters();
				ParamTransformSeq& base_trans_seq = childPest.get_base_par_tran_seq_4_mod();
				base_trans_seq.ctl2numeric_ip(par1);
				base_trans_seq.numeric2model_ip(par1);
				ParameterInfo pi = pest_scenario.get_ctl_parameter_info();
				int nadj_par = 0;
				for (auto par : par1)

				{
//					if ((pi.get_parameter_rec_ptr(par.first)->cycle == *icycle) ||
//						(pi.get_parameter_rec_ptr(par.first)->cycle < 0))
					if (cycle_in_range(*icycle,pi.get_parameter_rec_ptr(par.first)->dci))
					{

						if ((pi.get_parameter_rec_ptr(par.first)->tranform_type != ParameterRec::TRAN_TYPE::FIXED) &&
							(pi.get_parameter_rec_ptr(par.first)->tranform_type != ParameterRec::TRAN_TYPE::TIED))
							nadj_par++;
					}
				}

				cout << "...number of adjustable parameters in cycle " << *icycle << ": " << nadj_par << endl;
				fout_rec << "...number of adjustable parameters in cycle " << *icycle << ": " << nadj_par << endl;

				ObservationInfo oi = childPest.get_ctl_observation_info();
				int nnz_obs = 0;
				for (auto o : childPest.get_ctl_observations())
				{
					if (oi.get_observation_rec_ptr(o.first)->weight != 0.0)
						nnz_obs++;
				}
				cout << "...number of non-zero weighted observations in cycle " << *icycle << ": " << nnz_obs << endl;
				fout_rec << "...number of non-zero weighted observations in cycle " << *icycle << ": " << nnz_obs << endl;



				stringstream ss;

				/*ss << "num adj par: " << nadj_par << ", ";
				ss << "num nz obs: " << nnz_obs << ", ";*/
				ss << "num tpl files: " << childPest.get_tplfile_vec().size() << ", ";
				ss << "num ins files: " << childPest.get_insfile_vec().size() << endl;
				cout << ss.str();
				fout_rec << ss.str();

			}
		}
		if (pest_scenario.get_pestpp_options().get_debug_parse_only())
		{
			cout << endl << endl << "DEBUG_PARSE_ONLY is true, exiting..." << endl << endl;
			exit(0);
		}

		RunManagerAbstract* run_manager_ptr;

		if (cmdline.runmanagertype == CmdLine::RunManagerType::PANTHER_MASTER)
		{
			if (pest_scenario.get_control_info().noptmax == 0)
			{
				cout << endl << endl << "WARNING: 'noptmax' = 0 but using parallel run mgr.  This prob isn't what you want to happen..." << endl << endl;
			}

			run_manager_ptr = new RunManagerPanther(
				rns_file, cmdline.panther_port,
				file_manager.open_ofile_ext("rmr"),
				pest_scenario.get_pestpp_options().get_max_run_fail(),
				pest_scenario.get_pestpp_options().get_overdue_reched_fac(),
				pest_scenario.get_pestpp_options().get_overdue_giveup_fac(),
				pest_scenario.get_pestpp_options().get_overdue_giveup_minutes(),
				pest_scenario.get_pestpp_options().get_panther_echo(),
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
				pest_scenario.get_pestpp_options().get_additional_ins_delimiters(),
				pest_scenario.get_pestpp_options().get_num_tpl_ins_threads());
		}
		
		//generate a parent ensemble which includes all parameters across all cycles

		DataAssimilator da(pest_scenario, file_manager, output_file_writer, &performance_log, run_manager_ptr);
		ParameterEnsemble curr_pe(&pest_scenario);
		ObservationEnsemble curr_oe(&pest_scenario);
		ObservationEnsemble curr_noise(&pest_scenario);
		generate_global_ensembles(da, fout_rec, curr_pe, curr_oe, curr_noise);
		
		//prepare a phi csv file for all cycles
		string phi_file = file_manager.get_base_filename() + ".global.phi.actual.csv";
		ofstream f_phi(phi_file);
		if (f_phi.bad())
			throw runtime_error("error opening " + phi_file + " for writing");
		f_phi << "cycle,iteration,mean,standard_deviation,min,max";
		vector<string> init_real_names = curr_oe.get_real_names();
		for (auto real : init_real_names)
			f_phi << "," << pest_utils::lower_cp(real);
		f_phi << endl;
	
		cout << "initializing 'global' localizer" << endl;
		fout_rec << "initializing 'global' localizer" << endl;
		Localizer global_loc(&pest_scenario);
		global_loc.initialize(&performance_log);

		//ParameterEnsemble *_base_pe_ptr, FileManager *_file_manager_ptr, OutputFileWriter* _output_file_writer_ptr
		ParChangeSummarizer pcs(&curr_pe, &file_manager, &output_file_writer);

		for (auto icycle = assimilation_cycles.begin(); icycle != assimilation_cycles.end(); icycle++)
		{
			if (*icycle < start_cycle)
			{
				cout << "fast-forwarding past cycle " << *icycle << endl;
				continue;
			}
			// da_start_cycle, da_end_cycle
			cout << endl;
			cout << " =======================================" << endl;
			cout << " >>>> Assimilating data in cycle " << *icycle << endl;
			cout << " =======================================" << endl;

			fout_rec << endl;
			fout_rec << " =======================================" << endl;
			fout_rec << " >>>> Assimilating data in cycle " << *icycle << endl;
			fout_rec << " =======================================" << endl;

			performance_log.log_event("instantiating child pest object");

			Pest childPest = pest_scenario.get_child_pest(*icycle);

			OutputFileWriter output_file_writer(file_manager, childPest, restart_flag);

			if (verbose_level > 1) // 
			{
                output_file_writer.scenario_io_report(fout_rec);
				output_file_writer.scenario_pargroup_report(fout_rec);
				output_file_writer.scenario_par_report(fout_rec);
				output_file_writer.scenario_obs_report(fout_rec);
			}

			//------------------------------

			if (cmdline.runmanagertype == CmdLine::RunManagerType::PANTHER_MASTER)
			{
				performance_log.log_event("reinitializing panther master");
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
				performance_log.log_event("updating pars using da par cycle table info");
				//Parameters pars = childPest.get_ctl_parameters_4_mod();
				//TranFixed* fixed_ptr = base_trans_seq.get_fixed_ptr_4_mod();
				map<string, double> cycle_map = par_cycle_info[*icycle];
				for (auto item : cycle_map)
				{
					base_trans_seq.get_fixed_ptr_4_mod()->insert(item.first, item.second);
					childPest.get_ctl_parameters_4_mod().update_rec(item.first, item.second);
				}
			}
			//check for entries in the obs cycle table
			//cout << childPest.get_ctl_observations().get_rec("GAGE_1") << ", " << pest_scenario.get_observation_info_ptr()->get_weight("GAGE_1") << endl;
			if (obs_cycle_info.find(*icycle) != obs_cycle_info.end())
			{
				performance_log.log_event("updating obs using da obs cycle table info");
				map<string, double> cycle_map = obs_cycle_info[*icycle];
				map<string, double> weight_cycle_map;
				if (weight_cycle_info.find(*icycle) != weight_cycle_info.end())
					weight_cycle_map = weight_cycle_info[*icycle];
				ObservationInfo oi = pest_scenario.get_ctl_observation_info();
				childPest.set_observation_info(oi);
				for (auto tbl_obs_name : obs_in_tbl)
				{
					//if this obs is not in this cycle, give it a zero weight
					if (cycle_map.find(tbl_obs_name) == cycle_map.end())
					{
						oi.set_weight(tbl_obs_name, 0.0);
					}
					else
					{
						childPest.get_ctl_observations_4_mod().update_rec(tbl_obs_name, cycle_map[tbl_obs_name]);
						//check if this obs is in this cycle's weight info
						if (weight_cycle_map.find(tbl_obs_name) != weight_cycle_map.end())
						{
							oi.set_weight(tbl_obs_name, weight_cycle_map[tbl_obs_name]);
							//pest_scenario.get_observation_info_ptr()->set_weight(tbl_obs_name, weight_cycle_map[tbl_obs_name]);
						}
					}
				}
				childPest.set_observation_info(oi);

			}
			//cout << childPest.get_ctl_observations().get_rec("GAGE_1") << ", " << childPest.get_ctl_observation_info().get_weight("GAGE_1") << endl;

			if (pest_scenario.get_pestpp_options().get_ies_verbose_level() > 1) // todo: add da verbose
			{
			    fout_rec << "...verbose level > 1, cycle " << *icycle << " Pest Interface Summary:" << endl;
                output_file_writer.scenario_io_report(fout_rec);
				output_file_writer.scenario_pargroup_report(fout_rec);
				output_file_writer.scenario_par_report(fout_rec);
				output_file_writer.scenario_obs_report(fout_rec);
			}
			

			Parameters par1 = childPest.get_ctl_parameters();
			base_trans_seq.ctl2numeric_ip(par1);
			base_trans_seq.numeric2model_ip(par1);
			ParameterInfo pi = pest_scenario.get_ctl_parameter_info();//I change this as well. TODO use pest_scenrio instead of child

			int nadj_par = 0;
			for (auto par : par1)

			{
				// ayman:  base_trans_seq above was copied from parent pest without any changes; the following statement temporarly fix
				// the issue; permenat solution should occur during the creation of childpest
//				if ((pi.get_parameter_rec_ptr(par.first)->cycle == *icycle) ||
//					(pi.get_parameter_rec_ptr(par.first)->cycle < 0))
				if (cycle_in_range(*icycle,pi.get_parameter_rec_ptr(par.first)->dci))
				{
					if ((pi.get_parameter_rec_ptr(par.first)->tranform_type != ParameterRec::TRAN_TYPE::FIXED) &&
						(pi.get_parameter_rec_ptr(par.first)->tranform_type != ParameterRec::TRAN_TYPE::TIED))
						nadj_par++;
				}
			}

			cout << "...number of adjustable parameters in cycle " << *icycle << ": " << nadj_par << endl;
			fout_rec << "...number of adjustable parameters in cycle " << *icycle << ": " << nadj_par << endl;

			ObservationInfo oi = childPest.get_ctl_observation_info();
			int nnz_obs = 0;
			for (auto o : childPest.get_ctl_observations())
			{
				if (oi.get_observation_rec_ptr(o.first)->weight != 0.0)
					nnz_obs++;
			}
			cout << "...number of non-zero weighted observations in cycle " << *icycle << ": " << nnz_obs << endl;
			fout_rec << "...number of non-zero weighted observations in cycle " << *icycle << ": " << nnz_obs << endl;

			//ObjectiveFunc obj_func(&(childPest.get_ctl_observations()), &(childPest.get_ctl_observation_info()), &(childPest.get_prior_info()));

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
			run_manager_ptr->initialize(par_names, obs_names);
			performance_log.log_event("instantiating DA instance");
			DataAssimilator da(childPest, file_manager, output_file_writer, &performance_log, run_manager_ptr);
			da.initialize_dynamic_states();
			//da.da_type = pest_scenario.get_pestpp_options_ptr()->da_ctl_params.get_svalue("DA_TYPE");
			std::mt19937 rand_gen = da.get_rand_gen();
			vector<string> act_par_names = childPest.get_ctl_ordered_adj_par_names();
			performance_log.log_event("instantiating cycle parameter ensemble instance");
			ParameterEnsemble cycle_curr_pe(&childPest, &rand_gen, curr_pe.get_eigen(vector<string>(), act_par_names), curr_pe.get_real_names(), act_par_names);
			cycle_curr_pe.set_trans_status(curr_pe.get_trans_status());
			cycle_curr_pe.set_fixed_info(curr_pe.get_fixed_map());
			if (par_cycle_info.find(*icycle) != par_cycle_info.end())
            {
                map<pair<string, string>, double> fm = cycle_curr_pe.get_fixed_map();

			    for (auto& info : par_cycle_info.at(*icycle))
                {
			        for (auto& fi: fm)
                    {
			            if (fi.first.second == info.first)
                        {
			                fi.second = info.second;
                        }
                    }

                }
			    cycle_curr_pe.set_fixed_info(fm);
            }
			da.set_pe(cycle_curr_pe);
			obs_names = childPest.get_ctl_ordered_obs_names();
			ObservationEnsemble cycle_curr_oe(&childPest, &rand_gen, curr_oe.get_eigen(vector<string>(),obs_names), curr_oe.get_real_names(), obs_names);
			da.set_oe(cycle_curr_oe);
			if (nnz_obs > 0)
			{
				obs_names = childPest.get_ctl_ordered_nz_obs_names();
				ObservationEnsemble cycle_curr_noise(&childPest, &rand_gen, curr_noise.get_eigen(cycle_curr_oe.get_real_names(), obs_names), cycle_curr_oe.get_real_names(), obs_names);
				//correct for obs cycle table
				if ((pest_scenario.get_pestpp_options().get_da_weight_cycle_table().size() > 0) && (
				        !pest_scenario.get_pestpp_options().get_ies_no_noise())
				        )
                {
				    cout << "...'da_weight_cycle_table' passed, redrawing noise realizations" << endl;
				    fout_rec << "...'da_weight_cycle_table' passed, redrawing noise realizations" << endl;
				    da.initialize_obscov();
				    cycle_curr_noise.draw(cycle_curr_oe.shape().first,*da.get_obscov_ptr(),&performance_log,
                              pest_scenario.get_pestpp_options().get_ies_verbose_level(),fout_rec);
				    vector<string> names = cycle_curr_oe.get_real_names();
				    cycle_curr_noise.set_real_names(names);
                }

				else if (pest_scenario.get_pestpp_options().get_da_obs_cycle_table().size() > 0)
                {
                    cout << "...'da_obs_cycle_table' passed, translating noise realizations" << endl;
                    fout_rec << "...'da_obs_cycle_table' passed, translating noise realizations" << endl;
				    cycle_curr_noise.update_var_map();
				    Observations org_obs = pest_scenario.get_ctl_observations();
				    map<string,int> vmap = cycle_curr_noise.get_var_map();
                    for (auto o : childPest.get_ctl_observations())
                    {
                        if (oi.get_observation_rec_ptr(o.first)->weight != 0.0)
                        {
                            cycle_curr_noise.get_eigen_ptr_4_mod()->col(vmap[o.first]).array() +=o.second -  org_obs[o.first];
                        }
                    }
                }
				da.set_noise_oe(cycle_curr_noise);
			}
			da.set_localizer(global_loc);
			da.initialize(*icycle);
			
			write_global_phi_info(*icycle, f_phi, da, init_real_names);

			if (childPest.get_ctl_ordered_nz_obs_names().size() > 0)
			{

				if (pest_scenario.get_control_info().noptmax > 0) // 
				{
					da.da_update(*icycle);
					ss.str("");
					ss << file_manager.get_base_filename() << ".global." << *icycle << "." << da.get_iter() << ".pcs.csv";
					pcs.summarize(cycle_curr_pe, ss.str());
				}
				write_global_phi_info(*icycle, f_phi, da, init_real_names);
			}
			else
			{
				cout << "...Note: no non-zero-weighted observations in cycle " << *icycle << ", continuing..." << endl;
				fout_rec << "...Note: no non-zero-weighted observations in cycle " << *icycle << ", continuing..." << endl;
				performance_log.log_event("no non-zero-weighted observations in current cycle");

			}

			//replace all the pars used in this cycle in the parent parameter ensemble
			performance_log.log_event("updating global ensemble with cycle ensemble columns");
			cycle_curr_pe = da.get_pe();
			//if we lost some realizations...
			if (curr_pe.shape().first > cycle_curr_pe.shape().first)
			{
				vector<string> missing;
				vector<string> t = cycle_curr_pe.get_real_names();
				set<string> ccp_reals(t.begin(), t.end());
				for (auto r : curr_pe.get_real_names())
					if (ccp_reals.find(r) == ccp_reals.end())
						missing.push_back(r);
				fout_rec << "...dropping the following " << missing.size() << " realizations from global parameter ensemble:" << endl;
				cout << "...dropping " << missing.size() << " realizations from global parameter ensemble, see rec file for listing" << endl;
				int i = 0;
				for (auto m : missing)
				{
					fout_rec << m << ",";
					if (i > 10)
					{
						fout_rec << endl;
						i = 0;
					}
				}
				curr_pe.drop_rows(missing);
				curr_pe.reorder(t, vector<string>());
			}

			cycle_curr_pe.transform_ip(curr_pe.get_trans_status());
			curr_pe.replace_col_vals(cycle_curr_pe.get_var_names(), *cycle_curr_pe.get_eigen_ptr());
			
			cycle_curr_oe = da.get_oe();
			//if we lost some realizations...
			if (curr_oe.shape().first > cycle_curr_oe.shape().first)
			{
				vector<string> missing;
				vector<string> t = cycle_curr_oe.get_real_names();
				set<string> ccp_reals(t.begin(), t.end());
				for (auto r : curr_oe.get_real_names())
					if (ccp_reals.find(r) == ccp_reals.end())
						missing.push_back(r);
				fout_rec << "...dropping the following " << missing.size() << " realizations from global observation ensemble:" << endl;
				cout << "...dropping " << missing.size() << " realizations from global observation ensemble, see rec file for listing" << endl;
				int i = 0;
				for (auto m : missing)
				{
					fout_rec << m << ",";
					if (i > 10)
					{
						fout_rec << endl;
						i = 0;
					}
				}
				curr_oe.drop_rows(missing);
				curr_oe.reorder(t, vector<string>());
			}

			curr_oe.replace_col_vals(cycle_curr_oe.get_var_names(), *cycle_curr_oe.get_eigen_ptr());
			ss.str("");
			ss << ".global." << *icycle << ".oe.csv";
			curr_oe.to_csv(file_manager.get_base_filename() + ss.str());
			ss.str("");
			ss << ".global." << *icycle << ".pe.csv";
			curr_pe.to_csv(file_manager.get_base_filename() + ss.str());

			file_manager.close_all_files_ending_with("phi");

			//transfer the best (current) simulated final states to the inital states pars in the pe for the cycle
			//is the place to do this?
			da.transfer_dynamic_state_from_oe_to_pe(curr_pe, curr_oe);
			/*ss.str("");
			ss << "test_" << *icycle << ".csv";
			curr_pe.to_csv(ss.str());
			cout << curr_oe.get_eigen() << endl;
			cout << endl;*/

		} // end cycle loop
		fout_rec.close();
		cout << endl << endl << "pestpp-da analysis complete..." << endl;

        auto end = chrono::steady_clock::now();
        cout << "started at " << start_string << endl;
        cout << "finished at " << get_time_string() << endl;
        cout << "took " << chrono::duration_cast<chrono::seconds>(end - start).count() << " seconds" << endl;
        cout << flush;
        fout_rec << "started at " << start_string << endl;
        fout_rec << "finished at " << get_time_string() << endl;
        fout_rec << "took " << chrono::duration_cast<chrono::seconds>(end - start).count() << " seconds" << endl;

        return 0;
#ifndef _DEBUG
	}
	catch (exception & e)
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

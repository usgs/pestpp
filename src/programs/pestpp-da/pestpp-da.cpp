/**
 * @file pestpp-da.cpp
 * @brief Implementation of pestpp-da.
 */
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
#include "DaCycleDriver.h"
#include "RunManagerExternal.h"


using namespace std;
using namespace pest_utils;


/**
 * @brief Main.
 *
 * @param argc Description.
 * @param argv Description.
 *
 * @return Description.
 */
int main(int argc, char* argv[])
{
#ifndef _DEBUG
	try
	{
#endif

        string version = PESTPP_VERSION;
        cout << endl << endl;
        cout << "             pestpp-da: generalized iterative sequential/batch data assimilation" << endl << endl;
        //cout << "                     for PEST(++) datasets " << endl << endl;
        cout << "                            by the PEST++ development team" << endl;
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
		//jwhite - something weird is happening with the machine is busy and an existing
		//rns file is really large. so let's remove it explicitly and wait a few seconds before continuing...
		string rns_file = file_manager.build_filename("rns");
		int flag = remove(rns_file.c_str());


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
				cmdline.startup_report(frec,start_string);
				cmdline.startup_report(cout,start_string);
				PANTHERAgent yam_agent(frec);
				string ctl_file = "";
				try {

					// process traditional PEST control file
					ctl_file = file_manager.build_filename("pst");
					yam_agent.process_ctl_file(ctl_file);

				}
				catch (exception &e)
				{
                    frec << "Error processing control file: " << ctl_file << endl << endl;
                    frec << e.what() << endl << endl;
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
		PerformanceLog performance_log(file_manager.open_ofile_ext("log"));


		fout_rec << "             pestpp-da: generalized iterative sequential/batch data assimilation" << endl << endl;
		fout_rec << "                            by the PEST++ development team" << endl;
		fout_rec << endl;
		cmdline.startup_report(fout_rec,start_string);
		cmdline.startup_report(cout,start_string);

		// create pest run and process control file to initialize it
		Pest pest_scenario;
		try {
			performance_log.log_event("starting to process control file");
			pest_scenario.process_ctl_file(file_manager.open_ifile_ext("pst"), file_manager.build_filename("pst"), fout_rec);
			file_manager.close_file("pst");
			pest_scenario.assign_da_cycles(fout_rec);
			performance_log.log_event("finished processing control file");
		}
		catch (exception &e)
		{
			cerr << "Error processing control file: " << filename << endl << endl;
			cerr << e.what() << endl << endl;
			fout_rec << "Error processing control file: " << filename << endl << endl;
			fout_rec << e.what() << endl << endl;

			throw(e);
		}
		pest_scenario.check_inputs(fout_rec);
        stringstream ss;


		//Initialize OutputFileWriter to handle IO of supplementary files (.par, .par, .svd)
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

		ppo->apply_tool_defaults(PestppOptions::ToolType::DA, fout_rec);


        pest_scenario.get_pestpp_options().summary(fout_rec);

		//do all this up here so we can use the parse only option
		//to check the interface and then quit
		//process da (recurrent??) par cycle table
		vector <int> cycles_in_tables;
		map<int, map<string, double>> par_cycle_info = process_da_par_cycle_table(pest_scenario, cycles_in_tables, fout_rec);
		// process da obs cycle table
		set<string> obs_in_tbl; //we need this so we can set weights to zero in childpest of a value isn't listed for a given cycle
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


		ss.str("");
		ss << "...assimilating over " << assimilation_cycles.size() << " cycles from " << assimilation_cycles[0] << " to " << assimilation_cycles[assimilation_cycles.size()-1] << endl;
		cout << ss.str();
		fout_rec << ss.str();

        if ((pest_scenario.get_num_ext_file_maps() == 0) && (n_da_cycles > 1))
        {
            ss.str("");
            ss << endl;
            ss << "WARNING: no external control file sections found and multiple cycles detected.  ";
            ss << "         PESTPP-DA usually requires these external sections for sequential assimilation..." << endl;
            cout << ss.str();
            fout_rec << ss.str();
        }

		int start_cycle = pest_scenario.get_pestpp_options().get_da_hotstart_cycle();
		int max_cycle = assimilation_cycles[assimilation_cycles.size() - 1];
		if (start_cycle > max_cycle)
		{
			ss.str("");
			ss << "'da_hotstart_cycle' (" << start_cycle << ") greater than max cycle (" << max_cycle << ")";
			throw runtime_error(ss.str());
		}

		vector<string> errors;
		if ((pest_scenario.get_pestpp_options().get_check_tplins()) &&
		((pest_scenario.get_control_info().noptmax != 0) ||
		(pest_scenario.get_pestpp_options().get_debug_parse_only())))
		{
		    //create a DA instance here to check the dynamic state entries
		    RunManagerAbstract* run_manager_ptr = 0;
            DataAssimilator da(pest_scenario, file_manager, output_file_writer, &performance_log, run_manager_ptr);
            da.initialize_dynamic_states();
            da.sanity_checks();
            delete(run_manager_ptr);

            //now check the cycles
			for (auto icycle = assimilation_cycles.begin(); icycle != assimilation_cycles.end(); icycle++)
			{
				cout << endl;
				cout << " >>>> Checking data in cycle " << *icycle << endl;
				fout_rec << endl;
				fout_rec << " >>>> Checking data in cycle " << *icycle << endl;


				performance_log.log_event("instantiating child pest object");

				//Pest childPest;
				Pest childPest = pest_scenario.get_child_pest(*icycle);
                OutputFileWriter child_ofw(file_manager,childPest);

				childPest.check_inputs(fout_rec,false,true,*icycle);

				if (childPest.get_ctl_observations().size() == 0)
                {
				    ss.str("");
				    ss << "Error: no observations found for cycle " << *icycle << endl;
                    errors.push_back(ss.str());
                    continue;
                }
                if (childPest.get_ctl_parameters().size() == 0)
                {
                    ss.str("");
                    ss << "Error: no parameters found for cycle " << *icycle << endl;
                    errors.push_back(ss.str());
                    continue;
                }
                if (childPest.get_model_exec_info().tplfile_vec.size() == 0)
                {
                    ss.str("");
                    ss << "Error: no template files found for cycle " << *icycle << endl;
                    errors.push_back(ss.str());
                    continue;
                }
                if (childPest.get_model_exec_info().insfile_vec.size() == 0)
                {
                    ss.str("");
                    ss << "Error: no instruction files found for cycle " << *icycle << endl;
                    errors.push_back(ss.str());
                    continue;
                }
                try {
                    childPest.check_io(fout_rec, false);
                }
                catch (exception& e)
                {
                    ss.str("");
                    ss << "interface error for cycle " << *icycle << ": " << e.what() << endl;
                    errors.push_back(ss.str());
                    continue;
                }


				if ((obs_cycle_info.find(*icycle) != obs_cycle_info.end()) || (weight_cycle_info.find(*icycle) != weight_cycle_info.end()))
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
						}
					}
					for (auto tbl_obs_name : weights_in_tbl)
					{
						if (weight_cycle_map.find(tbl_obs_name) == weight_cycle_map.end())
						{
							oi.set_weight(tbl_obs_name, 0.0);
						}
						else
						{	
							oi.set_weight(tbl_obs_name, weight_cycle_map[tbl_obs_name]);
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



				ss.str("");

				ss << "...number of  template files in cycle " << *icycle << ": " << childPest.get_tplfile_vec().size() << endl;
				ss << "...number of instruction files in cycle " << *icycle << ": " << childPest.get_insfile_vec().size() << endl;
				cout << ss.str();
				fout_rec << ss.str();
				child_ofw.scenario_io_report(fout_rec);


			}
		}
		if (errors.size() > 0)
        {
		    ss.str("");
		    ss << errors.size() << " errors detected in cycle interface checking:" << endl;
		    cout << ss.str();
		    fout_rec << ss.str();
		    for (auto& e: errors)
            {
                cout << e;
                fout_rec << e;
            }
		    throw runtime_error("Errors detected in cycle interface checking");
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
				pest_scenario.get_ctl_ordered_obs_names(),
                pest_scenario.get_pestpp_options().get_panther_timeout_milliseconds(),
                pest_scenario.get_pestpp_options().get_panther_echo_interval_milliseconds(),
                pest_scenario.get_pestpp_options().get_panther_persistent_workers(),
				pest_scenario.get_pestpp_options().get_panther_ping_interval_secs());
			// host failure screening is a master-side policy, so it is set here rather than threaded
			// through the constructor's already long parameter list
			dynamic_cast<RunManagerPanther*>(run_manager_ptr)->set_max_failed_run_delta(
			    pest_scenario.get_pestpp_options().get_panther_agent_max_failed_run_delta());
			run_manager_ptr->initialize(pest_scenario.get_ctl_parameters(), pest_scenario.get_ctl_observations());
		}
        else if (cmdline.runmanagertype == CmdLine::RunManagerType::EXTERNAL)
        {

        	const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
        	run_manager_ptr = new RunManagerExternal(exi.comline_vec,
			exi.tplfile_vec, exi.inpfile_vec,
		   exi.insfile_vec, exi.outfile_vec,
			rns_file);
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
				pest_scenario.get_pestpp_options().get_num_tpl_ins_threads(),
				pest_scenario.get_pestpp_options().get_tpl_force_decimal());
			run_manager_ptr->initialize(pest_scenario.get_ctl_parameters(), pest_scenario.get_ctl_observations());
		}
		
		// DaCycleDriver::initialize() builds the global ensembles, the noptmax schedule, the phi
		// file, the localizer and the par change summarizer now. all of that is per-run state
		// that every caller needs, not setup for this one program, and having a copy here is
		// exactly how the two entry points would end up drifting apart.

		// the cycle sequence lives in DaCycleDriver now, not here. it used to be this loop, which
		// is why the c api could only ever run cycle one - somebody with a handle to a
		// DataAssimilator has one cycle, and everything that made it assimilation (the child
		// pest object per cycle, the posterior becoming the next prior, the dropped
		// realizations) was stuck in main() where nobody else could get at it.
		DaCycleDriver::RunManagerKind rm_kind = DaCycleDriver::RunManagerKind::OTHER;
		if (cmdline.runmanagertype == CmdLine::RunManagerType::PANTHER_MASTER)
			rm_kind = DaCycleDriver::RunManagerKind::PANTHER_MASTER;
		else if (cmdline.runmanagertype == CmdLine::RunManagerType::SERIAL)
			rm_kind = DaCycleDriver::RunManagerKind::SERIAL;
		else if (cmdline.runmanagertype == CmdLine::RunManagerType::EXTERNAL)
			rm_kind = DaCycleDriver::RunManagerKind::EXTERNAL;

		DaCycleDriver driver(pest_scenario, file_manager, &performance_log, run_manager_ptr,
			rm_kind, restart_flag);
		driver.set_pathname(pathname);
		driver.run_all_cycles();

		cout << endl << endl << "pestpp-da analysis complete..." << endl;
        fout_rec << endl << endl << "pestpp-da analysis complete..." << endl;

        auto end = chrono::steady_clock::now();
        cout << "started at " << start_string << endl;
        cout << "finished at " << get_time_string() << endl;


        fout_rec << "started at " << start_string << endl;
        fout_rec << "finished at " << get_time_string() << endl;
        fout_rec << "took " << setprecision(6) << (double)chrono::duration_cast<chrono::seconds>(end - start).count()/60.0 << " minutes" << endl;
        fout_rec << flush;
        cout << "took " << setprecision(6) << (double)chrono::duration_cast<chrono::seconds>(end - start).count()/60.0 << " minutes" << endl;
        cout << flush;
        if (pest_scenario.get_pestpp_options().get_save_all_runs())
        {
            string all_runs_file = file_manager.get_base_filename() + ".allruns.rns";
            RunStorage::print_persistent_summary(all_runs_file, cout);
            RunStorage::print_persistent_summary(all_runs_file, fout_rec);
        }
        fout_rec.close();
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
		return 1;
	}
#endif
}

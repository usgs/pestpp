
#include "RunManagerPanther.h" //needs to be first because it includes winsock2.h
//#include <vld.h> // Memory Leak Detection using "Visual Leak Detector"
#include <iostream>
#include <iomanip>
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
#include "FileManager.h"
#include "TerminationController.h"
#include "RunManagerSerial.h"
#include "RunManagerExternal.h"
#include "OutputFileWriter.h"
#include "PantherAgent.h"
#include "Serialization.h"
#include "system_variables.h"
#include "pest_error.h"
#include "RestartController.h"
#include "PerformanceLog.h"
#include "debug.h"
#include "DifferentialEvolution.h"
#include "RunManagerExternal.h"

#include "linear_analysis.h"
#include "logger.h"
#include "covariance.h"
#include "sequential_lp.h"

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
		cout << "             pestpp-opt - a tool for chance-constrained linear programming" << endl << endl;// , version " << version << endl << endl;
		cout << "                             by the PEST++ development team" << endl;
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

		vector<string>::const_iterator it_find, it_find_next;
		string next_item;
		string socket_str = "";
		
		if (cmdline.runmanagertype == CmdLine::RunManagerType::EXTERNAL)
		{
			cerr << "External run manager ('/e') not supported, please use panther instead" << endl;
			exit(1);
		}
		
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
				try {
					
					// process traditional PEST control file
					ctl_file = file_manager.build_filename("pst");
					yam_agent.process_ctl_file(ctl_file);
					
				}
				catch (PestError e)
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
		
		
		if (cmdline.runmanagertype == CmdLine::RunManagerType::GENIE)
		{
			cerr << "Genie run manager ('/g') no longer supported, please use PANTHER instead" << endl;
			exit(1);

		}


		RestartController restart_ctl;

		debug_initialize(file_manager.build_filename("dbg"));
		if (cmdline.restart)
		{
			cerr << "ERROR: pestpp-opt does not support restart options" << endl <<endl;

		}
		
		restart_ctl.get_restart_option() = RestartController::RestartOption::NONE;
		file_manager.open_default_files();
		
		ofstream &fout_rec = file_manager.rec_ofstream();
		PerformanceLog performance_log(file_manager.open_ofile_ext("log"));

		
		fout_rec << "             pestpp-opt version " << endl << endl;
		fout_rec << "       by the pestpp development team" << endl;
		fout_rec << endl;
		cmdline.startup_report(fout_rec,start_string);
		cmdline.startup_report(cout,start_string);
		// create pest run and process control file to initialize it
		Pest pest_scenario;
#ifndef _DEBUG
		try {
#endif
			performance_log.log_event("starting to process control file");
			pest_scenario.process_ctl_file(file_manager.open_ifile_ext("pst"), file_manager.build_filename("pst"),fout_rec);
			file_manager.close_file("pst");
			performance_log.log_event("finished processing control file");
#ifndef _DEBUG
		}
		catch (exception &e)
		{
			cerr << "Error processing control file: " << filename << endl << endl;
			cerr << e.what() << endl << endl;
			throw(e);
		}
#endif
		pest_scenario.check_inputs(fout_rec);
		if (pest_scenario.get_pestpp_options().get_debug_parse_only())
		{
			cout << endl << endl << "DEBUG_PARSE_ONLY is true, exiting..." << endl << endl;
			exit(0);
		}

		pest_scenario.get_pestpp_options_ptr()->set_iter_summary_flag(false);

		//if base jco arg read from control file, reset restart controller
		if (!pest_scenario.get_pestpp_options().get_basejac_filename().empty())
		{
			restart_ctl.get_restart_option() = RestartController::RestartOption::REUSE_JACOBIAN;
		}

		//Initialize OutputFileWriter to handle IO of supplementary files (.par, .par, .svd)
		//bool save_eign = pest_scenario.get_svd_info().eigwrite > 0;	=
		OutputFileWriter output_file_writer(file_manager, pest_scenario,false);
		
		output_file_writer.scenario_report(fout_rec, false);
		/*output_file_writer.scenario_io_report(fout_rec);
		output_file_writer.scenario_pargroup_report(fout_rec);
		output_file_writer.scenario_par_report(fout_rec);
		output_file_writer.scenario_obs_report(fout_rec);
		output_file_writer.scenario_pi_report(fout_rec);*/
		
		/*if (pest_scenario.get_pestpp_options().get_iter_summary_flag())
		{
			output_file_writer.write_par_iter(0, pest_scenario.get_ctl_parameters());
		}*/
		RunManagerAbstract *run_manager_ptr;
		if (cmdline.runmanagertype == CmdLine::RunManagerType::PANTHER_MASTER)
		{
			// using PANTHER run manager
			if (pest_scenario.get_control_info().noptmax == 0)
			{
				cout << endl << endl << "WARNING: 'noptmax' = 0 but using parallel run mgr.  This prob isn't what you want to happen..." << endl << endl;
			}
			const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
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

		else if (cmdline.runmanagertype == CmdLine::RunManagerType::EXTERNAL)
		{
			string rns_file = file_manager.build_filename("rns");
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
			if ((pest_scenario.get_pestpp_options().get_opt_skip_final()) &&
				(pest_scenario.get_pestpp_options().get_basejac_filename().size()) > 0 &&
				(pest_scenario.get_pestpp_options().get_hotstart_resfile().size()) > 0 &&
				(pest_scenario.get_control_info().noptmax == 1))
			{
				try
				{
					pest_scenario.check_io(fout_rec);
				}
				catch (...)
				{
					cout << "error checking I/O files...continuing" << endl;
				}

			}
			else
			{
				pest_scenario.check_io(fout_rec);
			}

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

		//setup the parcov, if needed
		//Covariance parcov;
		//if (pest_scenario.get_pestpp_options().get_use_parcov_scaling())
		/*double parcov_scale_fac = pest_scenario.get_pestpp_options().get_parcov_scale_fac();
		if (parcov_scale_fac > 0.0)
		{
			parcov.try_from(pest_scenario, file_manager);
		}*/
		const ParamTransformSeq &base_trans_seq = pest_scenario.get_base_par_tran_seq();

		ObjectiveFunc obj_func(&(pest_scenario.get_ctl_observations()), &(pest_scenario.get_ctl_observation_info()), &(pest_scenario.get_prior_info()));

//		TerminationController termination_ctl(pest_scenario.get_control_info().noptmax, pest_scenario.get_control_info().phiredstp,
//			pest_scenario.get_control_info().nphistp, pest_scenario.get_control_info().nphinored, pest_scenario.get_control_info().relparstp,
//			pest_scenario.get_control_info().nrelpar, pest_scenario.get_regul_scheme_ptr()->get_use_dynamic_reg(),
//			pest_scenario.get_regul_scheme_ptr()->get_phimaccept());

		

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
				//termination_ctl.check_last_iteration();
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
				file_manager.close_file("rei");
				run_manager_ptr->free_memory();
			}
			else
			{
				cout << "Model run failed.  No results were recorded." << endl << endl;
				fout_rec << "Model run failed.  No results were recorded." << endl << endl;
				exit(1);
			}
			//termination_ctl.set_terminate(true);
		}


		
		else
		{
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
			
			fout_rec << "   -----    Starting Optimization Iterations    ----    " << endl << endl;
			try {
                Covariance parcov;
                parcov.try_from(pest_scenario, file_manager);
                sequentialLP slp(pest_scenario, run_manager_ptr, parcov, &file_manager, output_file_writer,
                                 performance_log);
                slp.solve();
            }
            catch (exception &e)
            {
                fout_rec << "ERROR: " << e.what() << endl;
                throw runtime_error(e.what());

            }
            catch (...)
            {
                throw runtime_error("ERROR in sLP process");
            }
            fout_rec << "Number of forward model runs performed during optimization: " << run_manager_ptr->get_total_runs() << endl;
		}
		// clean up
		//fout_rec.close();
		delete run_manager_ptr;

		string case_name = file_manager.get_base_filename();
		file_manager.close_file("rst");
		pest_utils::try_clean_up_run_storage_files(case_name);

		cout << endl << endl << "pestpp-opt analysis complete..." << endl;
        auto end = chrono::steady_clock::now();
        cout << "started at " << start_string << endl;
        cout << "finished at " << get_time_string() << endl;
        cout << "took " << setprecision(6) << (double)chrono::duration_cast<chrono::seconds>(end - start).count()/60.0 << " minutes" << endl;
        cout << flush;
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
		//cout << "press enter to continue" << endl;
		//char buf[256];
		//OperSys::gets_s(buf, sizeof(buf));
		return 1;
	}
	catch (...)
	{
		cout << "Error condition prevents further execution" << endl;
		return 1;
	}
#endif
}

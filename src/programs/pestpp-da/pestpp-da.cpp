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

		ofstream &fout_rec = file_manager.rec_ofstream();
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
			fout_rec << "on host: \"" << w_get_hostname() << "\"" << endl << endl;
		}

		cout << endl;
		cout << "using control file: \"" << cmdline.ctl_file_name << "\"" << endl;
		cout << "in directory: \"" << OperSys::getcwd() << "\"" << endl;
		cout << "on host: \"" << w_get_hostname() << "\"" << endl << endl;

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



		//process da par cycle table
		map<int, map<string, double>> par_cycle_info = process_da_par_cycle_table(pest_scenario, fout_rec);
		// process da obs cycle table
		set<string> obs_in_tbl; //we need this so we can set weights to zero in childpest if a value isnt listed for a given cycle
		map<int, map<string, double>> obs_cycle_info = process_da_obs_cycle_table(pest_scenario, fout_rec, obs_in_tbl);
		//process weights table
		set<string> weights_in_tbl; 
		map<int, map<string, double>> weight_cycle_info = process_da_weight_cycle_table(pest_scenario, fout_rec, weights_in_tbl);
		
		vector<int> assimilation_cycles;
		pest_scenario.assign_da_cycles(fout_rec); 
		assimilation_cycles = pest_scenario.get_assim_cycles(fout_rec);

		//generate a parent ensemble which includes all parameters across all cycles
	
		DataAssimilator da(pest_scenario, file_manager, output_file_writer, &performance_log, run_manager_ptr);
		ParameterEnsemble curr_pe;
		ObservationEnsemble curr_oe;
		generate_parent_ensembles(da, fout_rec, curr_pe, curr_oe);
		
		//prepare a phi csv file for all cycles
		string phi_file = file_manager.get_base_filename() + ".parent.actual.phi.csv";
		ofstream f_phi(phi_file);
		if (f_phi.bad())
			throw runtime_error("error opening " + phi_file + " for writing");
		f_phi << "cycle,iteration,mean,standard_deviation,min,max";
		vector<string> init_real_names = curr_oe.get_real_names();
		for (auto real : init_real_names)
			f_phi << "," << pest_utils::lower_cp(real);
		f_phi << endl;


		// loop over assimilation cycles
		stringstream ss;
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

			performance_log.log_event("instantiating child pest object");
			
			Pest childPest;
			childPest = pest_scenario.get_child_pest(*icycle);
			//vector <string> xxxx=childPest.get_ctl_ordered_par_names();
			//childPest.get_pestpp_options.set_check_tplins(false);

			// -----------------------------  
			
			
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
			OutputFileWriter output_file_writer(file_manager, childPest, restart_flag);
			output_file_writer.scenario_io_report(fout_rec);
			if (pest_scenario.get_pestpp_options().get_ies_verbose_level() > 1) // todo: add da verbose
			{
				output_file_writer.scenario_pargroup_report(fout_rec);
				output_file_writer.scenario_par_report(fout_rec);
				output_file_writer.scenario_obs_report(fout_rec);
			}
			
			

			Parameters par1 = childPest.get_ctl_parameters();
			base_trans_seq.ctl2numeric_ip(par1);
			base_trans_seq.numeric2model_ip(par1);
			ParameterInfo pi = childPest.get_ctl_parameter_info();

			int nadj_par = 0;
			for (auto par : par1)
			{
				if ((pi.get_parameter_rec_ptr(par.first)->tranform_type != ParameterRec::TRAN_TYPE::FIXED) &&
					(pi.get_parameter_rec_ptr(par.first)->tranform_type != ParameterRec::TRAN_TYPE::TIED))
					nadj_par++;

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
			performance_log.log_event("instantiating DA instance");
			DataAssimilator da(childPest, file_manager, output_file_writer, &performance_log, run_manager_ptr);
			// use ies or da?
			da.use_ies = pest_scenario.get_pestpp_options_ptr()->get_da_use_ies();
			da.da_type = pest_scenario.get_pestpp_options_ptr()->da_ctl_params.get_svalue("DA_TYPE");
			std::mt19937 rand_gen = da.get_rand_gen();
			vector<string> act_par_names = childPest.get_ctl_ordered_adj_par_names();
			performance_log.log_event("instantiating cycle parameter ensemble instance");
			ParameterEnsemble cycle_curr_pe(&childPest, &rand_gen, curr_pe.get_eigen(vector<string>(),act_par_names), curr_pe.get_real_names(), act_par_names);
			cycle_curr_pe.set_trans_status(curr_pe.get_trans_status());

			//if (*icycle > 0)
			//{

			da.set_pe(cycle_curr_pe);
			//}
			
			if (da.use_ies)
			{
				da.initialize(*icycle);
			}
			else
			{
				da.da_initialize(*icycle);
			}

			write_parent_phi_info(*icycle, f_phi, da, init_real_names);

			if (childPest.get_ctl_ordered_nz_obs_names().size() > 0)
			{
				if (da.use_ies) // use ies
				{
					da.iterate_2_solution();
					/*curr_pe = da.get_pe();
					curr_pe.to_csv("cncnc.csv");*/
				}
				else // use da
				{
					if (pest_scenario.get_control_info().noptmax > 0) // 
					{
						da.da_upate();
					}
					/*curr_pe = da.get_pe();
					curr_pe.to_csv("cncnc.csv");*/
				}
				write_parent_phi_info(*icycle, f_phi, da, init_real_names);
			}
			else
			{
				cout << "...Note: no non-zero-weighted observations in cycle " << *icycle << ", continuing..." << endl;
				fout_rec << "...Note: no non-zero-weighted observations in cycle " << *icycle << ", continuing..." << endl;
				performance_log.log_event("no non-zero-weighted observations in current cycle");
				
			}
			//replace all the pars used in this cycle in the parent parameter ensemble
			performance_log.log_event("updating parent ensemble with cycle ensemble columns");
			
			//if we lost some realizations...
			if (curr_pe.shape().first > cycle_curr_pe.shape().first)
			{
				vector<string> missing;
				vector<string> t = cycle_curr_pe.get_real_names();
				set<string> ccp_reals(t.begin(), t.end());
				for (auto r : curr_pe.get_real_names())
					if (ccp_reals.find(r) == ccp_reals.end())
						missing.push_back(r);
				fout_rec << "...dropping the following " << missing.size() << " realizations from parent parameter ensemble:" << endl;
				cout << "...dropping " << missing.size() << " realizations from parent parameter ensemble, see rec file for listing" << endl;
				int i = 0;
				for (auto m : missing)
				{
					fout_rec << m << ",";
					if (i > 10)
					{
						fout_rec << endl;
						i == 0;
					}
				}
				curr_pe.drop_rows(missing);
				curr_pe.reorder(t, vector<string>());
			}

			cycle_curr_pe.transform_ip(curr_pe.get_trans_status());
			curr_pe.replace_col_vals(cycle_curr_pe.get_var_names(), *cycle_curr_pe.get_eigen_ptr());
			curr_pe.to_csv("cncnc.csv");

			ObservationEnsemble cycle_curr_oe = da.get_oe();
			//if we lost some realizations...
			if (curr_oe.shape().first > cycle_curr_oe.shape().first)
			{
				vector<string> missing;
				vector<string> t = cycle_curr_oe.get_real_names();
				set<string> ccp_reals(t.begin(), t.end());
				for (auto r : curr_oe.get_real_names())
					if (ccp_reals.find(r) == ccp_reals.end())
						missing.push_back(r);
				fout_rec << "...dropping the following " << missing.size() << " realizations from parent observation ensemble:" << endl;
				cout << "...dropping " << missing.size() << " realizations from parent observation ensemble, see rec file for listing" << endl;
				int i = 0;
				for (auto m : missing)
				{
					fout_rec << m << ",";
					if (i > 10)
					{
						fout_rec << endl;
						i == 0;
					}
				}
				curr_oe.drop_rows(missing);
				curr_oe.reorder(t, vector<string>());
			}

			curr_oe.replace_col_vals(cycle_curr_oe.get_var_names(), *cycle_curr_oe.get_eigen_ptr());
			ss.str("");
			ss << "." << *icycle << ".parent.oe.csv";
			curr_oe.to_csv(file_manager.get_base_filename()+ss.str());

			

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

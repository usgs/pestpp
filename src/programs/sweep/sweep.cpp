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
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <unordered_set>
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
#include "Jacobian.h"


using namespace std;
using namespace pest_utils;


map<string,int> prepare_parameter_csv(Parameters pars, ifstream &csv, bool forgive)
{
	if (!csv.good())
	{
		throw runtime_error("ifstream not good");
	}

	//process the header
	//any missing header labels will be marked to ignore those columns later
	string line;
	vector<string> header_tokens;
	if (!getline(csv, line))
		throw runtime_error("error reading header (first) line from csv file :");
	strip_ip(line);
	upper_ip(line);
	tokenize(line, header_tokens, ",", false);

	for (auto &t : header_tokens)
	{
		strip_ip(t);
	}
	//cout << tokens << endl;
	//vector<string> header_tokens = tokens;

	// check for parameter names that in the pest control file but that are missing from the csv file
	vector<string> missing_names;
	string name;
	for (auto &p : pars)
	if (find(header_tokens.begin(), header_tokens.end(), p.first) == header_tokens.end())
		missing_names.push_back(p.first);

	if (missing_names.size() > 0)
	{
		stringstream ss;
		ss << " the following pest control file parameter names were not found in the parameter csv file:" << endl;
		for (auto &n : missing_names) ss << n << endl;
		if (!forgive)
			throw runtime_error(ss.str());
		else
			cout << ss.str() << endl << "continuing anyway..." << endl;
	}

	if (header_tokens[header_tokens.size() - 1].size() == 0)
		header_tokens.pop_back();


	//build up a list of idxs to use
	vector<string> ctl_pnames = pars.get_keys();
	unordered_set<string> s_pnames(ctl_pnames.begin(), ctl_pnames.end());
	unordered_set<string>::iterator end = s_pnames.end();
	ctl_pnames.resize(0);
	vector<int> header_idxs;
	map<string, int> header_info;
	for (int i = 0; i < header_tokens.size(); i++)
	{
		//if (find(ctl_pnames.begin(), ctl_pnames.end(), header_tokens[i]) != ctl_pnames.end())
		if (s_pnames.find(header_tokens[i]) != end)
		{
			//header_idxs.push_back(i);
			header_info[header_tokens[i]] = i;
		}
	}
	return header_info;
}

//pair<vector<string>,vector<Parameters>> load_parameters_from_csv(map<string,int> &header_info, ifstream &csv, int chunk, const Parameters &ctl_pars, vector<string> &run_ids, vector<Parameters> &sweep_pars)
void load_parameters_from_csv(map<string, int>& header_info, ifstream& csv, int chunk, const Parameters& ctl_pars, vector<string>& run_ids, vector<Parameters>& sweep_pars)

{
	cout << endl;
	//process each parameter value line in the csv file
	int lcount = 1;
	//map<int,Parameters> sweep_pars;
	//vector<string> run_ids;
	//vector<Parameters> sweep_pars;
	run_ids.clear();
	sweep_pars.clear();
	sweep_pars.reserve(chunk);
	double val;
	string line;
	vector<string> tokens,names;
	vector<double> vals;
	Parameters pars = ctl_pars;
	string run_id;
	istringstream i;
	char c;
	
	while (getline(csv, line))
	{
		strip_ip(line);
		tokens.clear();
		vals.clear();
		names.clear();
		tokenize(line, tokens, ",", false);
		if (tokens[tokens.size() - 1].size() == 0)
			tokens.pop_back();

		if (tokens.size() != header_info.size()+1) // +1 for run id in first column
		{
			stringstream ss;
			ss << "error parsing csv file on line " << lcount << ": wrong number of tokens, ";
			ss << "expecting " << header_info.size() + 1 << ", found " << tokens.size();
			throw runtime_error(ss.str());
		}
		try
		{
			convert_ip(tokens[0], run_id);
		}
		catch (exception &e)
		{
			stringstream ss;
			ss << "error converting token '" << tokens[0] << "' to <int> run_id on line " << lcount << ": " << line << endl << e.what();
			throw runtime_error(ss.str());
		}

		//for (int i = 0; i < header_tokens.size(); i++)
		for (auto hi : header_info)
		{
			// if the par name of this column is in the passed-in par names, then replace the
			// value in pars
			//if (find(ctl_pnames.begin(), ctl_pnames.end(), header_tokens[i]) != ctl_pnames.end())
			//{
			try
			{
				//val = convert_cp<double>(tokens[hi.second]);
				val = stod(tokens[hi.second]);
			}
			catch (exception &e)
			{
				stringstream ss;
				ss << "error converting token '" << tokens[hi.second] << "' at location "<< hi.first << " to double on line " << lcount << ": " << line << endl << e.what();
				throw runtime_error(ss.str());
			}
			//vals.push_back(val);
			//names.push_back(hi.first);
			//pars[i] = val;
			pars.update_rec(hi.first, val);
			//}
		}
		//pars.update_without_clear(names, vals);
		sweep_pars.push_back(pars);
		//sweep_pars.emplace_back(pars);
		run_ids.push_back(run_id);
		lcount++;
		if (lcount > chunk)
			break;
		if (pars.size() > 10000)
			cout << lcount << "\r" << flush;
	}
	//csv.close();
	//return pair<vector<string>,vector<Parameters>> (run_ids,sweep_pars);
}


void prep_sweep_output_file(Pest &pest_scenario, ofstream &csv)
{
	//ofstream csv(pest_scenario.get_pestpp_options().get_sweep_output_csv_file());
	csv.open(pest_scenario.get_pestpp_options().get_sweep_output_csv_file());	
	if (!csv.good())
	{
		throw runtime_error("could not open sweep_output_csv_file for writing: " +
			pest_scenario.get_pestpp_options().get_sweep_output_csv_file());
	}
	csv.precision(numeric_limits<double>::digits10);
	csv << "run_id,input_run_id,failed_flag";
	csv << ",phi,meas_phi,regul_phi";
	for (auto &ogrp : pest_scenario.get_ctl_ordered_obs_group_names())
	{
		csv << ',' << pest_utils::lower_cp(ogrp);
	}
	for (auto &oname : pest_scenario.get_ctl_ordered_obs_names())
		csv << ',' << pest_utils::lower_cp(oname);
	csv << endl;
	csv.flush();
	//return csv;

}


void process_sweep_runs(ofstream &csv, Pest &pest_scenario, RunManagerAbstract* run_manager_ptr, vector<int> run_ids, vector<string> listed_run_ids, ObjectiveFunc obj_func,int total_runs_done)
{
	Parameters pars;
	Observations obs;
	double fail_val = -1.0E+10;
	int run_id;
	string listed_run_id;
	//for (auto &run_id : run_ids)
	for (int i = 0;i <run_ids.size();++i)
	{
		run_id = run_ids[i];
		listed_run_id = listed_run_ids[i];
		csv << run_id + total_runs_done;
		csv << ',' << listed_run_id;
		// if the run was successful
		if (run_manager_ptr->get_run(run_id, pars, obs))
		{
			PhiData phi_data = obj_func.phi_report(obs, pars, *(pest_scenario.get_regul_scheme_ptr()));
			csv << ",0";

			csv << ',' << phi_data.total();
			csv << ',' << phi_data.meas;
			csv << ',' << phi_data.regul;
			for (auto &obs_grp : pest_scenario.get_ctl_ordered_obs_group_names())
			{
				csv << ',' << phi_data.group_phi.at(obs_grp);
			}
			for (auto oname : pest_scenario.get_ctl_ordered_obs_names())
			{
				csv << ',' << obs[oname];
			}
			csv << endl;
		}
		//if the run bombed
		else
		{
			csv << ",1";
			csv << ",,,";
			for (auto &ogrp : pest_scenario.get_ctl_ordered_obs_group_names())
			{
				csv << ',';
			}
			for (int i = 0; i < pest_scenario.get_ctl_ordered_obs_names().size(); i++)
			{
				csv << ',' << fail_val;
			}
			csv << endl;
		}
	}
}



int main(int argc, char* argv[])
{
#ifndef _DEBUG
	try
	{
#endif
		string version = PESTPP_VERSION;
		cout << endl << endl;
		cout << "             pestpp-swp - a parameteric sweep utility, version " << version << endl;
		cout << "                     for PEST(++) datasets " << endl << endl;
		cout << "                 by the PEST++ development team" << endl << endl << endl;
		

		CmdLine cmdline(argc, argv);
		
		

		FileManager file_manager;
		string filename = cmdline.ctl_file_name;
		string pathname = ".";
		file_manager.initialize_path(get_filename_without_ext(filename), pathname);
		//jwhite - something weird is happening with the machine is busy and an existing
		//rns file is really large. so let's remove it explicitly and wait a few seconds before continuing...
		string rns_file = file_manager.build_filename("rns");
		int flag = remove(rns_file.c_str());
		//w_sleep(2000);
		//by default use the serial run manager.  This will be changed later if another
		//run manger is specified on the command line.
		
		if (cmdline.runmanagertype == CmdLine::RunManagerType::EXTERNAL)
		{
			cerr << "External run manager ('/e') not supported by sweep, please use panther instead" << endl;
			exit(1);
		}
		if (cmdline.runmanagertype == CmdLine::RunManagerType::GENIE)
		{
			cerr << "Genie run manager ('/e') deprecated, please use panther instead" << endl;
			exit(1);
		}
		
		if (cmdline.runmanagertype == CmdLine::RunManagerType::PANTHER_WORKER)
		{
			
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
		

		RestartController restart_ctl;

		//process restart and reuse jacobian directives
		
		bool restart_flag = false;
		bool save_restart_rec_header = true;

		debug_initialize(file_manager.build_filename("dbg"));
		if (cmdline.jac_restart)
		{
			throw runtime_error("/j option not supported by sweep");
		}
		if (cmdline.restart)
		{
			throw runtime_error("/r option not supported by sweep");
		}
		
		
		restart_ctl.get_restart_option() = RestartController::RestartOption::NONE;
		file_manager.open_default_files();
		

		ofstream &fout_rec = file_manager.rec_ofstream();
		PerformanceLog performance_log(file_manager.open_ofile_ext("log"));

		if (!restart_flag || save_restart_rec_header)
		{
			fout_rec << "             pestpp-swp.exe - a parameteric sweep utility" << endl << "for PEST(++) datasets " << endl << endl;
			fout_rec << "                 by the PEST++ developement team" << endl << endl << endl;
			fout_rec << endl;
			fout_rec << endl << endl << "version: " << version << endl;
			fout_rec << "binary compiled on " << __DATE__ << " at " << __TIME__ << endl << endl;
			fout_rec << "using control file: \"" << cmdline.ctl_file_name << "\"" << endl << endl;
			fout_rec << "in directory: \"" << OperSys::getcwd() << "\"" << endl;
			fout_rec << "on host: \"" << w_get_hostname() << "\"" << endl << endl;
		}

		cout << endl;
		cout << endl << endl << "version: " << version << endl;
		cout << "binary compiled on " << __DATE__ << " at " << __TIME__ << endl << endl;
		cout << "using control file: \"" << cmdline.ctl_file_name << "\"" << endl << endl;
		cout << "in directory: \"" << OperSys::getcwd() << "\"" << endl;
		cout << "on host: \"" << w_get_hostname() << "\"" << endl << endl;

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
		pest_scenario.check_inputs(fout_rec, true);
		
		OutputFileWriter ofw(file_manager, pest_scenario, false, false, 0);
		ofw.scenario_report(fout_rec, false);
		PestppOptions ppopt = pest_scenario.get_pestpp_options();

		fout_rec << "    sweep parameter csv file = " << left << setw(50) << ppopt.get_sweep_parameter_csv_file() << endl;
		fout_rec << "    sweep output csv file = " << left << setw(50) << ppopt.get_sweep_output_csv_file() << endl;
		fout_rec << "    sweep chunk size = " << left << setw(10) << ppopt.get_sweep_chunk() << endl;
		//fout_rec << "    sweep base run = " << left << setw(10) << ppopt.get_sweep_base_run() << endl;
		fout_rec << "    sweep forgive failed runs = " << left << setw(10) << ppopt.get_sweep_forgive() << endl;

		if (pest_scenario.get_pestpp_options().get_debug_parse_only())
		{
			cout << endl << endl << "DEBUG_PARSE_ONLY is true, exiting..." << endl << endl;
			exit(0);
		}

		// process the parameter csv file
		string par_csv_file;
		if (pest_scenario.get_pestpp_options().get_sweep_parameter_csv_file().empty())
		{
			//throw runtime_error("control file pest++ type argument 'SWEEP_PARAMETER_CSV_FILE' is required for sweep");
			cout << "pest++ arg SWEEP_PARAMETER_CSV_FILE not set, using 'SWEEP_IN.CSV'" << endl;
			par_csv_file = "sweep_in.csv";
		}
		else
			par_csv_file = pest_scenario.get_pestpp_options().get_sweep_parameter_csv_file();

		ifstream par_stream(par_csv_file);
		if (!par_stream.good())
		{
			throw runtime_error("could not open parameter csv file " + par_csv_file);
		}

		RunManagerAbstract *run_manager_ptr;
		if (cmdline.runmanagertype == CmdLine::RunManagerType::PANTHER_MASTER)
		{
			const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
			run_manager_ptr = new RunManagerPanther(
				file_manager.build_filename("rns"), cmdline.panther_port,
				file_manager.open_ofile_ext("rmr"),
				pest_scenario.get_pestpp_options().get_max_run_fail(),
				pest_scenario.get_pestpp_options().get_overdue_reched_fac(),
				pest_scenario.get_pestpp_options().get_overdue_giveup_fac(),
				pest_scenario.get_pestpp_options().get_overdue_giveup_minutes(),
				pest_scenario.get_pestpp_options().get_panther_echo());
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
				pest_scenario.get_pestpp_options().get_additional_ins_delimiters(),
				pest_scenario.get_pestpp_options().get_num_tpl_ins_threads());
		}


		const ParamTransformSeq &base_trans_seq = pest_scenario.get_base_par_tran_seq();
		ObjectiveFunc obj_func(&(pest_scenario.get_ctl_observations()), &(pest_scenario.get_ctl_observation_info()), &(pest_scenario.get_prior_info()));

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

		map<string, int> header_info;

		string par_ext = par_csv_file.substr(par_csv_file.size() - 3, par_csv_file.size());
		pest_utils::lower_ip(par_ext);

		Eigen::MatrixXd jco_mat;
		bool use_jco = false;
		vector<string> jco_col_names;
		if ((par_ext.compare("jcb") == 0) || (par_ext.compare("jco") == 0))
		{
			cout << "  ---  binary jco-type file detected for par_csv" << endl;
			use_jco = true;
			Jacobian jco(file_manager);
			jco.read(par_csv_file);
			cout << jco.get_matrix_ptr()->rows() << " runs found in binary jco-type file" << endl;
			//check that the jco is compatible with the control file
			vector<string> names = jco.get_base_numeric_par_names();
			jco_col_names = jco.get_sim_obs_names();
			set<string> jset(names.begin(), names.end());
			set<string> pset(pest_scenario.get_ctl_ordered_par_names().begin(), pest_scenario.get_ctl_ordered_par_names().end());
			vector<string> missing;
			set_symmetric_difference(jset.begin(), jset.end(), pset.begin(), pset.end(), back_inserter(missing));
			if (missing.size() > 0)
			{
				stringstream ss;
				ss << "binary jco file does not have the same parameters as pest control file.  The following parameters are not the same between the two: ";
				for (auto &m : missing)
					ss << m << " , ";
				throw runtime_error(ss.str());
			}
			/*for (int i = 0; i < names.size(); i++)
				header_info[names[i]] = i;*/
			cout << "  --- converting sparse JCO matrix to dense" << endl;
			jco_mat = jco.get_matrix(jco.get_sim_obs_names(), pest_scenario.get_ctl_ordered_par_names()).toDense();
		}

		else
		{
			header_info = prepare_parameter_csv(pest_scenario.get_ctl_parameters(),
				par_stream, pest_scenario.get_pestpp_options().get_sweep_forgive());
		}

		// prepare the output file
		ofstream obs_stream;
		prep_sweep_output_file(pest_scenario,obs_stream);

		int chunk = pest_scenario.get_pestpp_options().get_sweep_chunk();
		//pair<vector<string>,vector<Parameters>> sweep_par_info;

		//if desired, add the base run to the list of runs
		if (pest_scenario.get_pestpp_options().get_sweep_base_run())
		{
			throw runtime_error("base runs no longer supported by sweep");
			//sweep_pars[-999] = pest_scenario.get_ctl_parameters();
		}
		int total_runs_done = 0;
		vector<string> run_ids;
		vector<Parameters> sweep_pars;
		vector<int> irun_ids;
		while (true)
		{
			//read some realizations
			//sweep_pars.clear();
			cout << "reading par values...";
			if (use_jco)
			{
				//just use the total_runs_done counter as the run id
				vector<string> par_names = pest_scenario.get_ctl_ordered_par_names();
				Parameters par;
				//vector<Parameters> pars;
				sweep_pars.clear();
				sweep_pars.reserve(chunk);
				run_ids.clear();
				for (int i = 0; i < chunk; i++)
				{
					if (total_runs_done + i >= jco_mat.rows())
						break;
					par.update_without_clear(par_names, jco_mat.row(total_runs_done + i));
					sweep_pars.push_back(par);
					//run_ids.push_back(total_runs_done + i);
					run_ids.push_back(jco_col_names[total_runs_done + i]);

				}
				//sweep_par_info = pair<vector<string>, vector<Parameters>>(run_ids, pars);

			}
			else
			{
				try {
					performance_log.log_event("starting to read parameter csv file");
					load_parameters_from_csv(header_info, par_stream, chunk, pest_scenario.get_ctl_parameters(), run_ids,sweep_pars);
					performance_log.log_event("finished reading parameter csv file");
				}
				catch (exception &e)
				{
					stringstream ss;
					ss << "error processing parameter csv file: " << e.what();
					performance_log.log_event(ss.str());
					fout_rec << endl << ss.str() << endl;
					fout_rec.close();

					throw runtime_error(ss.str());
				}
			}
			cout << "done" << endl;
			// if there are no parameters to run, we are done
			//if (sweep_par_info.first.size() == 0)
			if (run_ids.size() == 0)
			{
				cout << "no more runs...done" << endl;
				break;
			}
			run_manager_ptr->reinitialize();

			cout << "starting runs " << total_runs_done << " --> " << total_runs_done + run_ids.size() << endl;

			// queue up some runs
			irun_ids.clear();
			for (auto &par : sweep_pars)
			{
				//Parameters temp = base_trans_seq.active_ctl2model_cp(par);
		        irun_ids.push_back(run_manager_ptr->add_run(base_trans_seq.active_ctl2model_cp(par)));

			}

			//make some runs

			run_manager_ptr->run();

			//process the runs
			cout << "processing runs...";
			//process_sweep_runs(obs_stream, pest_scenario, run_manager_ptr, run_ids, obj_func,total_runs_done);
			process_sweep_runs(obs_stream, pest_scenario, run_manager_ptr,irun_ids, run_ids, obj_func, total_runs_done);

			cout << "done" << endl;
			total_runs_done += run_ids.size();
		}

		// clean up
		fout_rec.close();
		obs_stream.close();
		delete run_manager_ptr;

		string case_name = file_manager.get_base_filename();
		file_manager.close_file("rst");
		pest_utils::try_clean_up_run_storage_files(case_name);

		cout << endl << endl << "PESTPP-SWP Analysis Complete..." << endl;
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

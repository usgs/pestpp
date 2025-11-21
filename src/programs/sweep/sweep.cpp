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
#include "RunManagerExternal.h"


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
	set<string> stokens(header_tokens.begin(),header_tokens.end());
	for (auto &p : pars)
	if (stokens.find(p.first) == stokens.end())
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

map<string,int> prepare_parameter_dense_binary(Parameters pars, ifstream &in, bool forgive, vector<string>& header_tokens)
{
    stringstream ss;
    if (!in.good())
    {
        throw runtime_error("ifstream not good in prepare_parameter_dense_binary()");
    }

    //process the header
    //any missing header labels will be marked to ignore those columns later
    int tmp1,n_col,tmp3;
    read_binary_matrix_header(in,tmp1,n_col,tmp3);
    if (!is_dense_binary_matrix(tmp1,n_col,tmp3))
    {
        throw runtime_error("prepare_parameter_dense_binary() file does not contain a dense binary matrix");
    }
    n_col *= -1;
    header_tokens.clear();
    header_tokens = read_dense_binary_col_names(in,n_col);
    ss.str("");
    ss << "..." << header_tokens.size() << " valid parameter entries found in dense binary parameter file ("
    << n_col << " total columns)" << endl;
    cout << ss.str();
    // check for parameter names that in the pest control file but that are missing from the csv file
    vector<string> missing_names;
    string name;
    set<string> stokens(header_tokens.begin(),header_tokens.end());
    for (auto &p : pars)
        if (stokens.find(p.first) == stokens.end())
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

void load_parameters_from_dense_binary(map<string, int>& header_info, ifstream& in, int chunk, const Parameters& ctl_pars, vector<string>& run_ids, vector<Parameters>& sweep_pars)

{
    cout << endl;
    //process each parameter value line in the csv file
    streampos current_pos = in.tellg();
    in.seekg(0,std::ios::beg);
    int tmp1,n_col,tmp3;
    read_binary_matrix_header(in,tmp1,n_col,tmp3);
    if (!is_dense_binary_matrix(tmp1,n_col,tmp3))
    {
        throw runtime_error("prepare_parameter_dense_binary() file does not contain a dense binary matrix");
    }
    n_col *= -1;
    in.seekg(current_pos,std::ios::beg);


    vector<vector<double>> rec_vecs;

    run_ids.clear();
    sweep_pars.clear();

    bool success = read_dense_binary_records(in,chunk,n_col,run_ids,rec_vecs);
    if (!success)
    {
        throw runtime_error("error processing dense binary file parameter records");
    }
    sweep_pars.reserve(rec_vecs.size());

    double val;
    string line;
    vector<string> tokens,names;
    vector<double> vals;
    Parameters pars = ctl_pars;
    string run_id;

    char c;

    for (int irec=0;irec < run_ids.size();irec++)
    {
        //for (int i = 0; i < header_tokens.size(); i++)
        for (auto hi : header_info)
        {
            pars.update_rec(hi.first, rec_vecs[irec][hi.second]);
            //}
        }
        //pars.update_without_clear(names, vals);
        sweep_pars.push_back(pars);
        //sweep_pars.emplace_back(pars);
        //run_ids.push_back(run_id);

        if (pars.size() > 10000)
            cout << irec << "\r" << flush;
    }
    //csv.close();
    //return pair<vector<string>,vector<Parameters>> (run_ids,sweep_pars);
}

void prep_sweep_output_file(Pest &pest_scenario, ofstream &out, bool& is_binary)
{
	//ofstream out(pest_scenario.get_pestpp_options().get_sweep_output_csv_file());
    vector<string> obs_names = pest_scenario.get_ctl_ordered_obs_names();
    string filename = pest_scenario.get_pestpp_options().get_sweep_output_csv_file();

    string obs_ext = filename.substr(filename.size() - 3, filename.size());
    pest_utils::lower_ip(obs_ext);
    if (obs_ext.compare("bin") == 0)
    {
        out.open(filename,ios::binary);
        if (!out.good()) {
            throw runtime_error("could not open sweep_output_file for writing: " +
                                filename);
        }
        prep_save_dense_binary(out,obs_names);
        is_binary = true;

    }
    else {
        out.open(filename);


        if (!out.good()) {
            throw runtime_error("could not open sweep_output_file for writing: " +
                                filename);
        }
        out << setprecision(numeric_limits<double>::digits10);
        out << "run_id,input_run_id,failed_flag";
        out << ",phi,meas_phi,regul_phi";
        for (auto &ogrp: pest_scenario.get_ctl_ordered_obs_group_names()) {
            out << ',' << pest_utils::lower_cp(ogrp);
        }
        for (auto &oname: pest_scenario.get_ctl_ordered_obs_names())
            out << ',' << pest_utils::lower_cp(oname);
        out << endl;
        out.flush();
        is_binary = false;
        //return out;
    }
}

void process_sweep_runs(ofstream &out, Pest &pest_scenario, RunManagerAbstract* run_manager_ptr, vector<int> run_ids, vector<string> listed_run_ids,
                        ObjectiveFunc obj_func,int total_runs_done,bool is_binary_output, PerformanceLog& performance_log)
{
	Parameters pars;
	Observations obs;
	double fail_val = -1.0E+100;
	int run_id;
	string listed_run_id;
    bool success;
	//for (auto &run_id : run_ids)

    vector<string> obs_names = pest_scenario.get_ctl_ordered_obs_names();
    vector<string> obs_group_names = pest_scenario.get_ctl_ordered_obs_group_names();
    Eigen::VectorXd vec(obs_names.size());
    vec.setConstant(fail_val);
    stringstream ss;
	for (int i = 0;i <run_ids.size();++i)
	{
		run_id = run_ids[i];
        ss.str("");
        ss << "processing run_id:" << run_id;
        performance_log.log_event(ss.str());
		listed_run_id = listed_run_ids[i];
        success = run_manager_ptr->get_run(run_id, pars, obs);
        if (is_binary_output)
        {
            vec.setConstant(fail_val);
            if (success)
            {
                vec = obs.get_data_eigen_vec(obs_names);
            }
            save_dense_binary(out,listed_run_id,vec);
        }
        else {
            out << run_id + total_runs_done;
            out << ',' << listed_run_id;
            // if the run was successful
            if (success) {
                performance_log.log_event("calculating phi info");
                PhiData phi_data = obj_func.phi_report(obs, pars, *(pest_scenario.get_regul_scheme_ptr()));

                out << ",0";

                out << ',' << phi_data.total();
                out << ',' << phi_data.meas;
                out << ',' << phi_data.regul;
                for (auto &obs_grp: obs_group_names) {
                    out << ',' << phi_data.group_phi.at(obs_grp);
                }
                for (auto &oname: obs_names) {
                    out << ',' << obs[oname];
                }
                out << endl;
            }
                //if the run bombed
            else {
                out << ",1";
                out << ",,,";
                for (auto &ogrp: obs_group_names) {
                    out << ',';
                }
                for (int i = 0; i < obs_names.size(); i++) {
                    out << ',' << fail_val;
                }
                out << endl;
            }
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
		cout << "             pestpp-swp - a parametric sweep utility, version " << version << endl;
		cout << "                     for PEST(++) datasets " << endl << endl;
		cout << "                 by the PEST++ development team" << endl << endl << endl;
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
		//jwhite - something weird is happening with the machine is busy and an existing
		//rns file is really large. so let's remove it explicitly and wait a few seconds before continuing...
		string rns_file = file_manager.build_filename("rns");
		int flag = remove(rns_file.c_str());
		//w_sleep(2000);
		//by default use the serial run manager.  This will be changed later if another
		//run manager is specified on the command line.

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

		fout_rec << "             pestpp-swp.exe - a parametric sweep utility" << endl << "for PEST(++) datasets " << endl << endl;
		fout_rec << "                 by the PEST++ development team" << endl << endl << endl;
		fout_rec << endl;
		cmdline.startup_report(fout_rec,start_string);
		cmdline.startup_report(cout,start_string);

		// create pest run and process control file to initialize it
		Pest pest_scenario;
		pest_scenario.set_default_dynreg();
		try {
			performance_log.log_event("starting to process control file");
			pest_scenario.process_ctl_file(file_manager.open_ifile_ext("pst"), file_manager.build_filename("pst"),fout_rec);
			file_manager.close_file("pst");
			performance_log.log_event("finished processing control file");
		}
		catch (exception &e)
		{
			cerr << "Error processing control file: " << filename << endl << endl;
			cerr << e.what() << endl << endl;
			fout_rec << "Error processing control file: " << filename << endl << endl;
			fout_rec << e.what() << endl << endl;
			fout_rec.close();
			throw(e);
		}
		pest_scenario.check_inputs(fout_rec, true);

		if (!pest_scenario.get_pestpp_options().get_sweep_include_regul_phi())
        {
		    pest_scenario.get_regul_scheme_ptr()->set_zero();
            pest_scenario.get_prior_info_ptr()->clear();
        }

		OutputFileWriter ofw(file_manager, pest_scenario, false, false, 0);
		ofw.scenario_report(fout_rec, false);
		PestppOptions ppopt = pest_scenario.get_pestpp_options();
        set<string> passed = ppopt.get_passed_args();
        if (passed.find("SWEEP_OUTPUT_FILE")==passed.end())
        {
            if (ppopt.get_save_dense())
            {
                fout_rec << "'save_dense' is true, resetting 'sweep_output_file' to 'sweep_output.bin'" << endl;
                pest_scenario.get_pestpp_options_ptr()->set_sweep_output_csv_file("sweep_output.bin");
            }
        }
        //get a fresh copy
        ppopt = pest_scenario.get_pestpp_options();
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


		fout_rec << "    sweep parameter file = " << left << setw(50) << par_csv_file << endl;
		fout_rec << "    sweep output file = " << left << setw(50) << ppopt.get_sweep_output_csv_file() << endl;
		fout_rec << "    sweep chunk size = " << left << setw(10) << ppopt.get_sweep_chunk() << endl;
		//fout_rec << "    sweep base run = " << left << setw(10) << ppopt.get_sweep_base_run() << endl;
		fout_rec << "    sweep forgive failed runs = " << left << setw(10) << ppopt.get_sweep_forgive() << endl;


		if (pest_scenario.get_pestpp_options().get_debug_parse_only())
		{
			cout << endl << endl << "DEBUG_PARSE_ONLY is true, exiting..." << endl << endl;
			exit(0);
		}



		ifstream par_stream(par_csv_file);
		if (!par_stream.good())
		{
			throw runtime_error("could not open parameter sweep file " + par_csv_file);
		}

		RunManagerAbstract *run_manager_ptr;
		if (cmdline.runmanagertype == CmdLine::RunManagerType::PANTHER_MASTER)
		{
			const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
            run_manager_ptr = new RunManagerPanther(
                    rns_file, cmdline.panther_port,
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
            fout_rec << "  ---  binary jco-type file detected for par_csv" << endl;
			use_jco = true;
			Jacobian jco(file_manager);
			jco.read(par_csv_file);
			cout << jco.get_matrix_ptr()->rows() << " runs found in binary jco-type file" << endl;
            fout_rec << jco.get_matrix_ptr()->rows() << " runs found in binary jco-type file" << endl;
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

		else if (par_ext.compare("csv") == 0)
		{
            cout << "  ---  csv file detected for par_csv" << endl;
            fout_rec << "  ---  csv file detected for par_csv" << endl;
			header_info = prepare_parameter_csv(pest_scenario.get_ctl_parameters(),
				par_stream, pest_scenario.get_pestpp_options().get_sweep_forgive());
		}
        else if (par_ext.compare("bin")==0)
        {
            cout << "  ---  dense binary file detected for par_csv" << endl;
            fout_rec << "  ---  dense binary file detected for par_csv" << endl;
			par_stream.close();
			par_stream.open(par_csv_file, ifstream::binary);
            vector<string> col_names;
            header_info = prepare_parameter_dense_binary(pest_scenario.get_ctl_parameters(),par_stream,
                                                         pest_scenario.get_pestpp_options().get_sweep_forgive(),col_names);
            vector<string> all_row_names = read_dense_binary_remaining_row_names(par_stream,col_names);
            cout << "..." << all_row_names.size() << " parameter sets found in dense binary file" << endl;
            fout_rec << "..." << all_row_names.size() << " parameter sets found in dense binary file" << endl;


        }
        else
        {
            throw runtime_error("unrecognized parameter sweep input file extension: '"+par_ext+"'");
        }

		// prepare the output file
		ofstream obs_stream;
        bool is_binary_output = false;

		prep_sweep_output_file(pest_scenario,obs_stream,is_binary_output);

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
                if (par_ext.compare("bin")==0)
                {
                    try {
                        performance_log.log_event("starting to read parameter dense binary file");
                        load_parameters_from_dense_binary(header_info, par_stream, chunk, pest_scenario.get_ctl_parameters(), run_ids,sweep_pars);
                        performance_log.log_event("finished reading parameter dense binary file");
                    }
                    catch (exception &e)
                    {
                        stringstream ss;
                        ss << "error processing parameter dense binary file: " << e.what();
                        performance_log.log_event(ss.str());
                        fout_rec << endl << ss.str() << endl;
                        fout_rec.close();

                        throw runtime_error(ss.str());
                    }

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

            performance_log.log_event("queuing runs");
			// queue up some runs
			irun_ids.clear();
			for (auto &par : sweep_pars)
			{
				//Parameters temp = base_trans_seq.active_ctl2model_cp(par);
		        irun_ids.push_back(run_manager_ptr->add_run(base_trans_seq.active_ctl2model_cp(par)));

			}

			//make some runs
            performance_log.log_event("calling run mgr");
			run_manager_ptr->run();

			//process the runs
			cout << "processing runs...";
            performance_log.log_event("processing runs");
			//process_sweep_runs(obs_stream, pest_scenario, run_manager_ptr, run_ids, obj_func,total_runs_done);
			process_sweep_runs(obs_stream, pest_scenario, run_manager_ptr,irun_ids, run_ids, obj_func, total_runs_done,is_binary_output, performance_log);

			cout << "done" << endl;
			total_runs_done += run_ids.size();
            int q = pest_utils::quit_file_found();
            if ((q == 1) || (q == 2))
            {
			    cout << "'pest.stp' found, quitting" << endl;
			    fout_rec << "'pest.stp' found, quitting" << endl;
			    break;
            }
            else if (q == 4) {
                cout << "...pest.stp found with '4'.  run mgr has returned control, removing file." << endl;
                fout_rec << "...pest.stp found with '4'.  run mgr has returned control, removing file." << endl;

                if (!pest_utils::try_remove_quit_file()) {
                    cout << "...error removing pest.stp file, bad times ahead..." << endl;
                    fout_rec << "...error removing pest.stp file, bad times ahead..." << endl;

                }
            }
		}

		// clean up

		obs_stream.close();
        par_stream.close();
		delete run_manager_ptr;

		string case_name = file_manager.get_base_filename();
		file_manager.close_file("rst");
		pest_utils::try_clean_up_run_storage_files(case_name);

		cout << endl << endl << "pestpp-swp analysis complete..." << endl;
        fout_rec << endl << endl << "pestpp-swp analysis complete..." << endl;
        auto end = chrono::steady_clock::now();
        cout << "started at " << start_string << endl;
        cout << "finished at " << get_time_string() << endl;
        cout << "took " << setprecision(6) << (double)chrono::duration_cast<chrono::seconds>(end - start).count()/60.0 << " minutes" << endl;
        cout << flush;
        fout_rec << "started at " << start_string << endl;
        fout_rec << "finished at " << get_time_string() << endl;
        fout_rec << "took " << setprecision(6) << (double)chrono::duration_cast<chrono::seconds>(end - start).count()/60.0 << " minutes" << endl;
        fout_rec.close();

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

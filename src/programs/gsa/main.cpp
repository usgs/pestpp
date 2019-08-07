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
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iomanip>
#include "config_os.h"
#include "MorrisMethod.h"
#include "sobol.h"
//#include "TornadoPlot.h"
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



using namespace std;
using namespace pest_utils;
using Eigen::MatrixXd;
using Eigen::VectorXd;


int main(int argc, char* argv[])
{
#ifndef _DEBUG
	try
	{
#endif
	string version = PESTPP_VERSION;
	cout << endl << endl;
	cout << "             pestpp-sen: a tool for global sensitivity analysis" << endl << endl;
	cout << "                       by The PEST++ Development Team" << endl;
	cout << endl << endl << "version: " << version << endl;
	cout << "binary compiled on " << __DATE__ << " at " << __TIME__ << endl << endl;
	
	// build commandline
	string commandline = "";
	for(int i=0; i<argc; ++i)
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
	enum class RunManagerType {SERIAL, PANTHER, GENIE};

	if (argc >=2) {
		complete_path = argv[1];
	}
	else {
		cerr << "--------------------------------------------------------" << endl;
		cerr << "usage:" << endl << endl;
		cerr << "    serial run manager:" << endl;
		cerr << "        gsa pest_ctl_file.pst" << endl << endl;
		cerr << "    PANTHER master:" << endl;
		cerr << "        gsa control_file.pst /H :port" << endl;
		cerr << "    PANTHER worker:" << endl;
		cerr << "        gsa control_file.pst /H hostname:port " << endl << endl;
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
	//Check for PANTHER worker
	vector<string>::const_iterator it_find, it_find_next;
	string next_item;
	string socket_str = "";
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
					yam_agent.process_ctl_file(ctl_file);
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
	//Check for PANTHER master
	else if (it_find != cmd_arg_vec.end())
	{
		// using PANTHER run manager
		run_manager_type = RunManagerType::PANTHER;
		socket_str = next_item;
	}

	it_find = find(cmd_arg_vec.begin(), cmd_arg_vec.end(), "/g");
	next_item.clear();
	it_find = find(cmd_arg_vec.begin(), cmd_arg_vec.end(), "/g");
	next_item.clear();
	if (it_find != cmd_arg_vec.end())
	{
		cerr << "Genie run manager ('/g') no longer supported, please use PANTHER instead" << endl;
		return 1;

	}

	ofstream &fout_rec = file_manager.open_ofile_ext("rec");
	fout_rec << "             pestpp-sen: a tool for global sensitivity analysis" << endl << endl;
	fout_rec << "                         by The PEST++ Development Team" << endl << endl;
	fout_rec << endl << endl << "version: " << version << endl;
	fout_rec << "binary compiled on " << __DATE__ << " at " << __TIME__ << endl << endl;
	fout_rec << "using control file: \"" << complete_path << "\"" << endl;
	fout_rec << "in directory: \"" << OperSys::getcwd() << "\"" << endl << endl;

	
	cout << endl;
	cout << "using control file: \"" << complete_path << "\"" << endl;
	cout << "in directory: \"" << OperSys::getcwd() << "\"" << endl << endl;

	// create pest run and process control file to initialize it
	Pest pest_scenario;
	pest_scenario.set_defaults();
	try {
		pest_scenario.process_ctl_file(file_manager.open_ifile_ext("pst"), file_manager.build_filename("pst"),fout_rec);
		file_manager.close_file("pst");
		pest_scenario.check_inputs(fout_rec);
	}
	catch(PestError e)
	{
		cerr << "Error prococessing control file: " << filename << endl << endl;
		cerr << e.what() << endl << endl;
		fout_rec << "Error prococessing control file: " << filename << endl << endl;
		fout_rec << e.what() << endl;
		fout_rec.close();
		//throw(e);
		return 1;
	}


	RunManagerAbstract *run_manager_ptr;
	if (run_manager_type == RunManagerType::PANTHER)
	{
		string port = argv[3];
		strip_ip(port);
		strip_ip(port, "front", ":");
		const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
		run_manager_ptr = new RunManagerPanther (
			file_manager.build_filename("rns"), port,
			file_manager.open_ofile_ext("rmr"),
			pest_scenario.get_pestpp_options().get_max_run_fail(),
			pest_scenario.get_pestpp_options().get_overdue_reched_fac(),
			pest_scenario.get_pestpp_options().get_overdue_giveup_fac(),
			pest_scenario.get_pestpp_options().get_overdue_giveup_minutes());
	}
	else
	{
		const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
		run_manager_ptr = new RunManagerSerial(exi.comline_vec,
		exi.tplfile_vec, exi.inpfile_vec, exi.insfile_vec, exi.outfile_vec,
		file_manager.build_filename("rns"), pathname);
	}

	cout << endl;
	fout_rec << endl;
	cout << "using control file: \"" <<  complete_path << "\"" << endl;
	fout_rec << "using control file: \"" <<  complete_path << "\"" << endl;


	enum class GSA_RESTART { NONE, RESTART };
	GSA_RESTART gsa_restart = GSA_RESTART::NONE;
	//process restart and  reuse jacibian directives
	vector<string>::const_iterator it_find_r = find(cmd_arg_vec.begin(), cmd_arg_vec.end(), "/r");
	if (it_find_r != cmd_arg_vec.end())
	{
		gsa_restart = GSA_RESTART::RESTART;
	}
	else
	{
		gsa_restart = GSA_RESTART::NONE;
	}
	//OutputFileWriter(FileManager &_file_manager, Pest &_pest_scenario, bool restart_flag = false, bool _save_rei = true, int _eigenwrite = 0);
	OutputFileWriter output_writer(file_manager, pest_scenario, false, false);
	ofstream &frec = file_manager.rec_ofstream();
	output_writer.scenario_report(frec);
	output_writer.scenario_io_report(frec);
	output_writer.scenario_par_report(frec);
	output_writer.scenario_obs_report(frec);

	pest_scenario.check_inputs(frec);
	pest_scenario.check_io();

	//map<string, string> gsa_opt_map;
	//process .gsa file
	string gsa_filename = file_manager.get_base_filename() + ".gsa";
	/*if (!check_exist_in(gsa_filename))
	{
		cout << "WARNING: " << gsa_filename << " not found, using standard settings and Method of Morris:" << endl;
		cout << "     MORRIS_P: 4" << endl;
		cout << "     MORRIS_R: 4" << endl;
		cout << "     MORRIS_DELTA: 0.666" << endl;
		gsa_opt_map["METHOD"] = "MORRIS";
		gsa_opt_map["MORRIS_DELTA"] = ".666666";
		gsa_opt_map["MORRIS_P"] = "4";
		gsa_opt_map["MORRIS_R"] = "4";
		gsa_opt_map["RAND_SEED"] = "2";
	}
	else
	{
		try
		{
			gsa_opt_map = GsaAbstractBase::process_gsa_file(file_manager.open_ifile_ext("gsa"), file_manager);
			file_manager.close_file("gsa");
		}
		catch (PestError e)
		{
			cerr << "Error prococessing .gsa file: " << file_manager.build_filename("gsa") << endl << endl;
			cerr << e.what() << endl << endl;
			throw(e);
		}
	}*/
	if (check_exist_in(gsa_filename))
	{
		cout << "WARNING: use of .gsa files is deprecated - .gsa file '" << gsa_filename << "' is being ignored, please use '++' args";
		return 1;
	}

	PestppOptions *pp_ptr = pest_scenario.get_pestpp_options_ptr();
	map<string, string> gsa_opt_map = pp_ptr->get_arg_map();
	/*gsa_opt_map["METHOD"] = pp_ptr->get_gsa_method();
	gsa_opt_map["MORRIS_DELTA"] = pp_ptr->get_gsa_morris_delta();
	gsa_opt_map["MORRIS_P"] = pp_ptr->get_gsa_morris_p();
	gsa_opt_map["MORRIS_R"] = pp_ptr->get_gsa_morris_r;
	gsa_opt_map["MORRIS_OBS_SEN"] = pp_ptr->get_gsa_morris_obs_sen();
*/
	//Build Transformation with ctl_2_numberic
	ParamTransformSeq base_partran_seq(pest_scenario.get_base_par_tran_seq());
	Parameters ctl_par = pest_scenario.get_ctl_parameters();


	//Build Transformation with ctl_2_numberic
	ObjectiveFunc obj_func(&(pest_scenario.get_ctl_observations()), &(pest_scenario.get_ctl_observation_info()), &(pest_scenario.get_prior_info()));
	ModelRun model_run(&obj_func, pest_scenario.get_ctl_observations());
	const set<string> &log_trans_pars = base_partran_seq.get_log10_ptr()->get_items();
	auto method = gsa_opt_map.find("GSA_METHOD");

	GsaAbstractBase* gsa_method = nullptr;
	if (method == gsa_opt_map.end() || method->second == "MORRIS")
	{
		int morris_r = 4;
		int morris_p = 4;

		double morris_delta = .666;
		bool default_delta = true;
		bool calc_pooled_obs = false;
		bool calc_morris_obs_sen = true;
		auto morris_r_it = gsa_opt_map.find("GSA_MORRIS_R");
		if (morris_r_it != gsa_opt_map.end())
		{
			convert_ip(morris_r_it->second, morris_r);
		}
		auto morris_p_it = gsa_opt_map.find("GSA_MORRIS_P");
		if (morris_p_it != gsa_opt_map.end())
		{
			convert_ip(morris_p_it->second, morris_p);
		}
		auto morris_d_it = gsa_opt_map.find("GSA_MORRIS_DELTA");
		if (morris_d_it != gsa_opt_map.end())
		{
			
			convert_ip(morris_d_it->second, morris_delta);
			default_delta = false;
		}
		auto morris_pool_it = gsa_opt_map.find("GSA_MORRIS_POOLED_OBS");
		if (morris_pool_it != gsa_opt_map.end())
		{
			string pooled_obs_flag = morris_pool_it->second;
			upper_ip(pooled_obs_flag);
			if (pooled_obs_flag == "TRUE") calc_pooled_obs = true;
		}

		auto morris_obs_sen_it = gsa_opt_map.find("GSA_MORRIS_OBS_SEN");
		if (morris_obs_sen_it != gsa_opt_map.end())
		{
			string obs_sen_flag = morris_obs_sen_it->second;
			upper_ip(obs_sen_flag);
			if (obs_sen_flag == "FALSE") calc_morris_obs_sen = false;
		}

		if (default_delta) morris_delta = morris_p / (2.0 * (morris_p - 1));

		MorrisMethod *m_ptr = new MorrisMethod(pest_scenario, file_manager, &obj_func,
			base_partran_seq, morris_p, morris_r, morris_delta, calc_pooled_obs,
			calc_morris_obs_sen, GsaAbstractBase::PARAM_DIST::uniform, 1.0);
		gsa_method = m_ptr;
		m_ptr->process_pooled_var_file();
		
		frec << endl << endl << endl << "Method of Morris settings:" << endl;

		frec << scientific << left << setw(30) << " morris_r " << morris_r << endl;
		frec << scientific << left << setw(30) << " morris_p " << morris_p << endl;
		frec << scientific << left << setw(30) << " morris_delta " << morris_delta << endl << endl;
		

	}
	
	else if (method != gsa_opt_map.end() && method->second == "SOBOL")
	{
		GsaAbstractBase::PARAM_DIST par_dist = GsaAbstractBase::PARAM_DIST::uniform;
		string par_dist_str = "UNIFORM";
		int n_sample = 100;
		auto sob_n_sam_it = gsa_opt_map.find("GSA_SOBOL_SAMPLES");
		if (sob_n_sam_it != gsa_opt_map.end())
		{
			convert_ip(sob_n_sam_it->second, n_sample);
		}
		auto sob_p_dist_it = gsa_opt_map.find("GSA_SOBOL_PAR_DIST");
		if (sob_p_dist_it != gsa_opt_map.end())
		{
			par_dist_str = sob_p_dist_it->second;
			upper_ip(par_dist_str);
			if (par_dist_str == "NORM") par_dist = GsaAbstractBase::PARAM_DIST::normal;
			else if (par_dist_str == "UNIF") par_dist = GsaAbstractBase::PARAM_DIST::uniform;
			else
			{
				ostringstream str;
				str << "SOBOL_PAR_DIST(" << par_dist_str << "):  \"" << par_dist_str << "\" is an invalid distribuation type";
				throw PestError(str.str());
			}
		}

		gsa_method = new Sobol(pest_scenario, file_manager, &obj_func,
			base_partran_seq, n_sample, par_dist, 1.0);

		frec << endl << endl << endl << "Method of Sobol settings:" << endl;

		frec << scientific << left << setw(30) << " n_sample " << n_sample << endl;
		frec << scientific << left << setw(30) << " sobol par dist " <<par_dist_str << endl;
		
	}
	else
	{
		throw PestError("A valid method for computing the sensitivity must be specified in the control file");
	}

	auto morris_r_it = gsa_opt_map.find("RAND_SEED");
	if (morris_r_it != gsa_opt_map.end())
	{
		unsigned int seed = convert_cp<unsigned int>(morris_r_it->second);
		gsa_method->set_seed(seed);
	}
	frec << scientific << left << setw(30) << " gsa random seed " << gsa_method->get_seed() << endl;
	// make model runs
	if (gsa_restart == GSA_RESTART::NONE)
	{
		//Allocates Space for Run Manager.  This initializes the model parameter names and observations names.
		//Neither of these will change over the course of the simulation
		cout << endl;
		cout << "Building model run parameter sets..." << endl;
		run_manager_ptr->initialize(base_partran_seq.ctl2model_cp(ctl_par), pest_scenario.get_ctl_observations());

		Parameters model_pars = base_partran_seq.ctl2model_cp(ctl_par);
		run_manager_ptr->reinitialize();
		gsa_method->assemble_runs(*run_manager_ptr);
	}
	else
	{
		run_manager_ptr->initialize_restart(file_manager.build_filename("rns"));
	}
	cout << endl;
	cout << "Performing model runs..." << endl;
	run_manager_ptr->run();

	cout << "Calculating sensitivities..." << endl;
	gsa_method->calc_sen(*run_manager_ptr, model_run);
	file_manager.close_file("srw");
	file_manager.close_file("msn");
	file_manager.close_file("orw");
	delete run_manager_ptr;
	cout << endl << endl << "Simulation Complete..." << endl;
	//cout << endl << "Simulation Complete - Press RETURN to close window" << endl;
	//char buf[256];
	//OperSys::gets_s(buf, sizeof(buf));
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

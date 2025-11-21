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
#include "RunManagerExternal.h"
#include "RunManagerPanther.h"
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
	string pathname = ".";
	file_manager.initialize_path(get_filename_without_ext(cmdline.ctl_file_name), pathname);

	if (cmdline.runmanagertype == CmdLine::RunManagerType::PANTHER_WORKER)
	{
		// This is a PANTHER worker, start PEST++ as a PANTHER worker
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
	/*if (cmdline.runmanagertype == CmdLine::RunManagerType::EXTERNAL)
	{
		cerr << "external run manager ('/e') no longer supported, please use PANTHER instead" << endl;
		exit(1);

	}*/

	ofstream &fout_rec = file_manager.open_ofile_ext("rec");
	fout_rec << "             pestpp-sen: a tool for global sensitivity analysis" << endl << endl;
	fout_rec << "                         by The PEST++ Development Team" << endl << endl;
	fout_rec << endl;
	cmdline.startup_report(fout_rec,start_string);
	cmdline.startup_report(cout, start_string);

	// create pest run and process control file to initialize it
	Pest pest_scenario;
	try 
	{
		pest_scenario.process_ctl_file(file_manager.open_ifile_ext("pst"), file_manager.build_filename("pst"),fout_rec);
		file_manager.close_file("pst");
		//pest_scenario.check_inputs(fout_rec);
	}
	catch(PestError e)
	{
		cerr << "Error prococessing control file: " << cmdline.ctl_file_name << endl << endl;
		cerr << e.what() << endl << endl;
		fout_rec << "Error prococessing control file: " << cmdline.ctl_file_name << endl << endl;
		fout_rec << e.what() << endl;
		fout_rec.close();
		//throw(e);
		return 1;
	}
	pest_scenario.check_inputs(fout_rec);
	//OutputFileWriter(FileManager &_file_manager, Pest &_pest_scenario, bool restart_flag = false, bool _save_rei = true, int _eigenwrite = 0);
	
	OutputFileWriter ofw(file_manager,pest_scenario,false,false,0);
	//ofw.scenario_report(fout_rec, false);

	if (pest_scenario.get_pestpp_options().get_debug_parse_only())
	{
		cout << endl << endl << "DEBUG_PARSE_ONLY is true, exiting..." << endl << endl;
		exit(0);
	}


	RunManagerAbstract *run_manager_ptr;
	if (cmdline.runmanagertype == CmdLine::RunManagerType::PANTHER_MASTER)
	{
		const ModelExecInfo &exi = pest_scenario.get_model_exec_info();
		run_manager_ptr = new RunManagerPanther (
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

	cout << endl;
	fout_rec << endl;
	cout << "using control file: \"" <<  cmdline.ctl_file_name << "\"" << endl;
	fout_rec << "using control file: \"" <<  cmdline.ctl_file_name << "\"" << endl;


	enum class GSA_RESTART { NONE, RESTART };
	GSA_RESTART gsa_restart = GSA_RESTART::NONE;
	//process restart and  reuse jacibian directives
	if (cmdline.jac_restart)
	{
		cerr << "jacobian restart ('/j') not supported in PESTPP-SEN";
		exit(1);
	}
	if (cmdline.restart)
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
	//output_writer.scenario_io_report(frec);
	//output_writer.scenario_par_report(frec);
	//output_writer.scenario_obs_report(frec);

	pest_scenario.check_inputs(frec);
	pest_scenario.check_io(frec);

	//check for zero-weighted obs groups and warn
	stringstream ss;


    if (pest_scenario.get_ctl_ordered_nz_obs_names().size() == 0)
    {
        ss.str("");
        ss << "WARNING: all observations are zero weighted - resetting all weights to 1.0";
        cout << ss.str() << endl;
        fout_rec << ss.str() << endl;
        ObservationInfo* oi_ptr = pest_scenario.get_observation_info_ptr();

        for (auto& o : pest_scenario.get_ctl_observations())
        {
            oi_ptr->set_weight(o.first,1.0);
        }
    }
    else
    {
        ss.str("");
        ss << endl << "Note: only non-zero weighted observations contribute to" << endl;
        ss << "      the phi and group phi sensitivity metrics.  Please" << endl;
        ss << "      make sure this is what you want..." << endl << endl;
        cout << ss.str();
        fout_rec << ss.str();
    }

	//map<string, string> gsa_opt_map;
	//process .gsa file
	string gsa_filename = file_manager.get_base_filename() + ".gsa";
	
	if (check_exist_in(gsa_filename))
	{
		cout << "ERROR: use of .gsa files is deprecated - .gsa file '" << gsa_filename << "' is being ignored, please use '++' args";
		fout_rec<< "ERROR: use of .gsa files is deprecated - .gsa file '" << gsa_filename << "' is being ignored, please use '++' args";
		return 1;
	}

	PestppOptions *pp_ptr = pest_scenario.get_pestpp_options_ptr();
	map<string, string> gsa_opt_map = pp_ptr->get_arg_map();
	
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
		bool calc_pooled_obs = true;
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
		int n_sample = pest_scenario.get_pestpp_options().get_gsa_sobol_samples();
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
				str << "SOBOL_PAR_DIST(" << par_dist_str << "):  \"" << par_dist_str << "\" is an invalid distribution type";
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
	//gsa_method->set_seed(2);
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

	string case_name = file_manager.get_base_filename();
	file_manager.close_file("rst");
	pest_utils::try_clean_up_run_storage_files(case_name);

	cout << endl << endl << "pestpp-sen analysis complete..." << endl;
	fout_rec << endl << endl << "pestpp-sen analysis complete..." << endl;
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
		return 1;
	}
	catch (...)
	{
		cout << "Error condition prevents further execution" << endl;
		return 1;
	}

#endif
}

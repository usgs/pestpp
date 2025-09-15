// pestpp-ies.cpp : Defines the entry point for the console application.
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
#include "EnsembleSmoother.h"

using namespace std;
using namespace pest_utils;


int main(int argc, char* argv[])
{
#ifndef _DEBUG
	try {
#endif
        string version = PESTPP_VERSION;
        cout << endl << endl;
        cout << "             pestpp-ies: a GLM iterative ensemble smoother" << endl << endl;
        //cout << "                     for PEST(++) datasets " << endl << endl;
        cout << "                   by the PEST++ development team" << endl;
        cout << endl << endl << "version: " << version << endl;
        cout << "binary compiled on " << __DATE__ << " at " << __TIME__ << endl << endl;
        auto start = chrono::steady_clock::now();
        string start_string = get_time_string();
        cout << "started at " << start_string << endl;
        CmdLine cmdline(argc, argv);

        if (quit_file_found()) {
            cerr << "'pest.stp' found, please remove this file " << endl;
            return 1;
        }


        FileManager file_manager;
        string pathname = ".";
        file_manager.initialize_path(get_filename_without_ext(cmdline.ctl_file_name), pathname);
        //jwhite - something weird is happening with the machine is busy and an existing
        //rns file is really large. so let's remove it explicitly and wait a few seconds before continuing...
        string rns_file = file_manager.build_filename("rns");
        int flag = remove(rns_file.c_str());
        //w_sleep(2000);
        //by default use the serial run manager.  This will be changed later if another
        //run manager is specified on the command line.

        if (cmdline.runmanagertype == CmdLine::RunManagerType::PANTHER_WORKER) {
            try {
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
                catch (exception &e) {
                    frec << "Error processing control file: " << ctl_file << endl << endl;
                    frec << e.what() << endl << endl;
                    cerr << "Error processing control file: " << ctl_file << endl << endl;
                    cerr << e.what() << endl << endl;
                    throw (e);
                }
                catch (...) {
                    cerr << "Error processing control file" << endl;
                    throw runtime_error("error processing control file");
                }
                yam_agent.start(cmdline.panther_host_name, cmdline.panther_port);
            }
            catch (PestError &perr) {
                cerr << perr.what();
                throw (perr);
            }

            cout << endl << "Work Done..." << endl;
            exit(0);
        }

        if (cmdline.runmanagertype == CmdLine::RunManagerType::GENIE) {
            cerr << "Genie run manager ('/g') no longer supported, please use PANTHER instead" << endl;
            return 1;

        }
        if (cmdline.runmanagertype == CmdLine::RunManagerType::EXTERNAL) {
            cerr << "external run manager ('/e') no supported in PESTPP-IES, please use PANTHER instead" << endl;
            return 1;

        }

        RestartController restart_ctl;

        //process restart and reuse jacobian directives
        bool restart_flag = false;
        bool save_restart_rec_header = true;

        debug_initialize(file_manager.build_filename("dbg"));
        if (cmdline.jac_restart) {
            throw runtime_error("/j option not supported by pestpp-ies");
        } else if (cmdline.restart) {
            throw runtime_error("/r option not supported by pestpp-ies");
        } else {
            restart_ctl.get_restart_option() = RestartController::RestartOption::NONE;
            file_manager.open_default_files();
        }

        ofstream &fout_rec = file_manager.rec_ofstream();
        PerformanceLog performance_log(file_manager.open_ofile_ext("log"));

        if (!restart_flag || save_restart_rec_header) {
            fout_rec << "             pestpp-ies - a GLM iterative Ensemble Smoother" << endl
                     << "                      for PEST(++) datasets " << endl << endl;
            fout_rec << "                 by the PEST++ development team" << endl << endl << endl;
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
        //try {
        performance_log.log_event("starting to process control file");
        try {
            pest_scenario.process_ctl_file(file_manager.open_ifile_ext("pst"), file_manager.build_filename("pst"),
                                           fout_rec);
        }
        catch (exception &e)
        {
            fout_rec << "Error processing control file: " << file_manager.build_filename("pst") << endl << endl;
            fout_rec << e.what() << endl << endl;
            cerr << "Error processing control file: " << file_manager.build_filename("pst") << endl << endl;
            cerr << e.what() << endl << endl;
            throw(e);
        }
        file_manager.close_file("pst");
        performance_log.log_event("finished processing control file");
        /*}
        catch (PestError e)
        {
            cerr << "Error prococessing control file: " << cmdline.ctl_file_name << endl << endl;
            cerr << e.what() << endl << endl;
            fout_rec << "Error prococessing control file: " << cmdline.ctl_file_name << endl << endl;
            fout_rec << e.what() << endl << endl;
            fout_rec.close();
            throw(e);
        }*/
        //pest_scenario.clear_ext_files();
        pest_scenario.check_inputs(fout_rec);

        //Initialize OutputFileWriter to handle IO of supplementary files (.par, .par, .svd)
        //bool save_eign = pest_scenario.get_svd_info().eigwrite > 0;
        pest_scenario.get_pestpp_options_ptr()->set_iter_summary_flag(false);
        OutputFileWriter output_file_writer(file_manager, pest_scenario, restart_flag);
        output_file_writer.scenario_report(fout_rec, false);
        //output_file_writer.scenario_io_report(fout_rec);
        if (pest_scenario.get_pestpp_options().get_ies_verbose_level() > 1) {
            output_file_writer.scenario_pargroup_report(fout_rec);
            output_file_writer.scenario_par_report(fout_rec);
            output_file_writer.scenario_obs_report(fout_rec);
        }

        //reset some default args for ies here:
        PestppOptions *ppo = pest_scenario.get_pestpp_options_ptr();
        set<string> pp_args = ppo->get_passed_args();
        if (pp_args.find("MAX_RUN_FAIL") == pp_args.end()) {
            ppo->set_max_run_fail(1);
            fout_rec << "...resetting max_run_fail to 1" << endl;
        }

        if (pp_args.find("OVERDUE_GIVEUP_FAC") == pp_args.end()) {
            ppo->set_overdue_giveup_fac(2.0);
            fout_rec << "...resetting overdue_giveup_fac to 2.0" << endl;
        }

        if (pp_args.find("OVERDUE_RESCHED_FAC") == pp_args.end()) {
            ppo->set_overdue_reched_fac(1.15);
            fout_rec << "...resetting overdue_resched_fac to 1.15" << endl;
        }


        RunManagerAbstract *run_manager_ptr;


        if (cmdline.runmanagertype == CmdLine::RunManagerType::PANTHER_MASTER) {
            if (pest_scenario.get_control_info().noptmax == 0) {
                cout << endl << endl
                     << "WARNING: 'noptmax' = 0 but using parallel run mgr.  This prob isn't what you want to happen..."
                     << endl << endl;
            }
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
        } else {
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
                                                   pest_scenario.get_pestpp_options().get_tpl_force_decimal(),
                                                   pest_scenario.get_pestpp_options().get_panther_echo());
        }

        const ParamTransformSeq &base_trans_seq = pest_scenario.get_base_par_tran_seq();
        ObjectiveFunc obj_func(&(pest_scenario.get_ctl_observations()), &(pest_scenario.get_ctl_observation_info()),
                               &(pest_scenario.get_prior_info()));

        Parameters cur_ctl_parameters = pest_scenario.get_ctl_parameters();
        //Allocates Space for Run Manager.  This initializes the model parameter names and observations names.
        //Neither of these will change over the course of the simulation


        run_manager_ptr->initialize(base_trans_seq.ctl2model_cp(cur_ctl_parameters),
                                    pest_scenario.get_ctl_observations());

        IterEnsembleSmoother ies(pest_scenario, file_manager, output_file_writer, &performance_log, run_manager_ptr);
        ies.initialize();
        if (pest_scenario.get_pestpp_options().get_debug_parse_only()) {
            cout << endl << endl << "DEBUG_PARSE_ONLY is true, exiting..." << endl << endl;
            exit(0);
        }

        int q = pest_utils::quit_file_found();
        if ((q == 1) || (q == 2)) {
            cout << "...'pest.stp' found, quitting" << endl;
            fout_rec << "...'pest.stp' found, quitting" << endl;

        } else {
            if (q == 4) {
                cout << "...pest.stp found with '4'.  run mgr has returned control, removing file." << endl;
                fout_rec << "...pest.stp found with '4'.  run mgr has returned control, removing file." << endl;

                if (!pest_utils::try_remove_quit_file()) {
                    cout << "...error removing pest.stp file, bad times ahead..." << endl;
                    fout_rec << "...error removing pest.stp file, bad times ahead..." << endl;

                }
            }
            ies.iterate_2_solution();
        }

        ies.finalize();



        // clean up

        delete run_manager_ptr;
        string case_name = file_manager.get_base_filename();
        file_manager.close_file("rst");
        pest_utils::try_clean_up_run_storage_files(case_name);

        cout << endl << endl << "pestpp-ies analysis complete..." << endl;
        fout_rec << endl << endl << "pestpp-ies analysis complete..." << endl;
        auto end = chrono::steady_clock::now();
        cout << "started at " << start_string << endl;
        cout << "finished at " << get_time_string() << endl;
        //cout << "took " << chrono::duration_cast<chrono::seconds>(end - start).count() << " seconds" << endl;
        cout << "took " << setprecision(6)
             << (double) chrono::duration_cast<chrono::seconds>(end - start).count() / 60.0 << " minutes" << endl;
        cout << flush;
        fout_rec << "started at " << start_string << endl;
        fout_rec << "finished at " << get_time_string() << endl;
        fout_rec << "took " << setprecision(6)
                 << (double) chrono::duration_cast<chrono::seconds>(end - start).count() / 60.0 << " minutes" << endl;
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

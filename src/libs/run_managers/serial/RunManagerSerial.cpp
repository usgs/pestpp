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
#include "RunManagerSerial.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <iterator>
#include <cassert>
#include <cstring>
#include <map>
#include <algorithm>
#include "system_variables.h"
#include "Transformable.h"
#include "utilities.h"
#include "model_interface.h"

using namespace std;
using namespace pest_utils;


RunManagerSerial::RunManagerSerial(const vector<string> _comline_vec,
	const vector<string> _tplfile_vec, const vector<string> _inpfile_vec,
	const vector<string> _insfile_vec, const vector<string> _outfile_vec,
	const string &stor_filename, const string &_run_dir, int _max_run_fail,
	bool fill_tpl_zeros, string additional_ins_delimiters, int _num_threads,
	bool tpl_force_decimal)
	: RunManagerAbstract(_comline_vec, _tplfile_vec, _inpfile_vec,
	_insfile_vec, _outfile_vec, stor_filename, _max_run_fail),
	run_dir(_run_dir), mi(_tplfile_vec,_inpfile_vec,_insfile_vec,_outfile_vec, _comline_vec)
{
	mi.set_additional_ins_delimiters(additional_ins_delimiters);
	mi.set_fill_tpl_zeros(fill_tpl_zeros);
	mi.set_num_threads(_num_threads);
	mi.set_tpl_force_decimal(tpl_force_decimal);
    mi.set_sleep_ms(5);
	cout << "              starting serial run manager ..." << endl << endl;
	mgr_type = RUN_MGR_TYPE::SERIAL;
}

void RunManagerSerial::run_async(pest_utils::thread_flag* terminate, pest_utils::thread_flag* finished, exception_ptr& run_exception,
                             Parameters* pars, Observations* obs)
{
    mi.run(terminate,finished,run_exception, pars, obs);
}

void RunManagerSerial::run(Parameters* pars, Observations* obs)
{
    thread_flag f_terminate(false);
    thread_flag f_finished(false);
    exception_ptr run_exception = nullptr;
    stringstream ss;
    thread run_thread(&RunManagerSerial::run_async, this, &f_terminate, &f_finished, std::ref(run_exception),
                      pars, obs);


    while (true)
    {

        if (run_exception) {
            try {
                rethrow_exception(run_exception);
            }
            catch (exception &e)
            {
                ss.str("");
                ss << "exception raised by run thread: " << std::endl;
                ss << e.what() << std::endl;
                cout << ss.str();
            }
            f_terminate.set(true);
            run_thread.join();
            break;
        }
        //check if the runner thread has finished
        if (f_finished.get())
        {
            break;
        }
        int q = pest_utils::quit_file_found();
        if ((q == 1) || (q == 2) || (q == 4))
        {
            f_terminate.set(true);
            run_thread.join();
            break;
        }

    }
    run_thread.join();
    if (run_exception)
        rethrow_exception(run_exception);
}

void RunManagerSerial::run()
{
	int success_runs = 0;
	int prev_sucess_runs = 0;
	int failed_runs = 0;
	const vector<string> &par_name_vec = file_stor.get_par_name_vec();
	const vector<string> &obs_name_vec = file_stor.get_obs_name_vec();

	stringstream message;
	int nruns = get_outstanding_run_ids().size();
	message.str("");
	message << endl << endl << endl << "    ---  starting serial run manager for " << nruns << " runs ---    " << endl << endl << endl;
	std::cout << message.str();

	std::vector<double> obs_vec;
	vector<int> run_id_vec;
	
	std::chrono::system_clock::time_point start_time_all = std::chrono::system_clock::now();
	int run_status;
	string info_txt;
	int group_id = -1;
	double info_value;
	while (!(run_id_vec = get_outstanding_run_ids()).empty())
	{
		for (int i_run : run_id_vec)
		{
            int q = pest_utils::quit_file_found();
            if ((q == 1) || (q == 2) || (q == 4))
            {
		        cout << "'pest.stp' found, treating run " << i_run << " as a fail" << endl;
                update_run_failed(i_run);
                failed_runs++;
                continue;
            }

			file_stor.get_info(i_run, run_status, info_txt, info_value);
			try
			{
				ofstream fout("run.info");
				if (fout.good())
				{
					fout << "run_id, " << i_run << endl;
					fout << "group_id, " << group_id << endl;
					fout << "info_txt," << info_txt << endl;
				}
				fout.close();
			}
			catch (...)
			{

			}
			std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();
			try
			{	
				Observations obs;
				vector<double> par_values;
				Parameters pars;
				file_stor.get_parameters(i_run, pars);
				obs_vec.resize(obs_name_vec.size(), RunStorage::no_data);
				obs.clear();
				obs.insert(obs_name_vec, obs_vec);
				//mi.run(&pars, &obs);

				run(&pars, &obs);
				
				//OperSys::chdir(run_dir.c_str());
				success_runs += 1;
				//std::cout << string(message.str().size(), '\b');
				//message.str("");
				//message << "(" << success_runs << "/" << nruns << " runs complete)";
				//std::cout << message.str();
				file_stor.update_run(i_run, pars, obs);

			}
			catch (const std::exception& ex)
			{
				update_run_failed(i_run);
				failed_runs++;
				cerr << "  Error running model: " << ex.what() << endl;
				cerr << "  Aborting model run" << endl << endl;

			}
			catch (...)
			{
				update_run_failed(i_run);
				failed_runs++;
				cerr << "  Error running model" << endl;
				cerr << "  Aborting model run" << endl << endl;
			}
			message.str("");
			message << endl << endl << "-->" << pest_utils::get_time_string() << " run complete, took: " << pest_utils::get_duration_sec(start_time) << " seconds";
			message << endl << "-->" << success_runs << " of " << nruns << " complete, "<<  failed_runs << " failed" << endl << endl << endl;
			std::cout << message.str();
		}
	}
	total_runs += success_runs;
	message.str("");
	message << endl << endl << endl << "    ---  serial run manager runs summary:  ---    " << endl;
	message << "    " << success_runs << " of " << nruns << " complete, " << failed_runs << " failed" << endl;
	message << "    " << "process took : " << pest_utils::get_duration_sec(start_time_all) << " seconds" << endl << endl << endl;
	std::cout << message.str();

	if (prev_sucess_runs > 0)
	{
		message.str("");
		message << " and " << prev_sucess_runs << " additional run completed previously";
		std::cout << message.str();
	}
	
	if (failed_runs > 0)
	{			cout << endl << endl;
		cout << "WARNING: " << failed_runs << " out of " << nruns << " runs failed" << endl << endl;
		cout << "    failed run ids:" << endl;
		int i = 1;
		for (auto id : get_failed_run_ids())
		{
			cout << id << ",";
			i++;
			if (i > 10)
			{
				cout << endl << "    ";
				i = 1;
			}		
		}
	}
	std::cout << endl << endl;
	if (init_sim.size() == 0)
	{
		vector<double> pars;
		int status = file_stor.get_run(0, pars, init_sim);
	}
}


RunManagerSerial::~RunManagerSerial(void)
{
}

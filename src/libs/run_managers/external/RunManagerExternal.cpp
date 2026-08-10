/**
 * @file RunManagerExternal.cpp
 * @brief Implementation of RunManagerExternal.
 */
#include "network_wrapper.h"
#include <iostream>
#include "RunManagerExternal.h"
#include "utilities.h"
#include <system_variables.h>


using namespace std;

RunManagerExternal::RunManagerExternal(vector<string> _comline_vec,
	vector<string> _tplfile_vec, vector<string> _inpfile_vec,
	vector<string> _insfile_vec, vector<string> _outfile_vec,
	string &stor_filename, int _max_n_failure, int _sleep_ms)
	: RunManagerAbstract(_comline_vec, _tplfile_vec, _inpfile_vec,
	_insfile_vec, _outfile_vec, stor_filename, _max_n_failure), sleep_ms(_sleep_ms)
{
	cout << "              starting external run manager ..." << endl << endl;
}

/**
 * @brief Run.
 */
void RunManagerExternal::run()
{

	std::chrono::system_clock::time_point start_time = chrono::system_clock::now();

    cout << pest_utils::get_time_string() << " external run manager calling forward run command(s)" << endl;

	// The same pest_utils::run_command() the model interface uses. This file used to carry its
	// own copy of the spawn/poll/kill logic, and the two had already drifted - this one always
	// narrated, and had no cancellation at all. Sharing means the windows job object and the
	// posix process-group kill are written down once, which is what makes grandchildren die
	// with their parent on both platforms.
	//
	// nullptr for the terminate flag: the external run manager has nothing to cancel from. It
	// hands the whole batch to the model and waits, which is the one case the parameter was
	// made nullable for.
	pest_utils::run_commands(comline_vec, nullptr, true, sleep_ms);

	cout << pest_utils::get_time_string() << " forward run command(s) finished, took " << pest_utils::get_duration_sec(start_time) << " seconds" << endl;
	//cout << pest_utils::get_time_string() << " re-initializing run storage file";
	//file_stor.init_restart(file_stor.get_filename());

}


/**
 * @brief Destructor for .
 */
RunManagerExternal::~RunManagerExternal()
{
}

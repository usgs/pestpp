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

void RunManagerExternal::run()
{

	std::chrono::system_clock::time_point start_time = chrono::system_clock::now();

    cout << pest_utils::get_time_string() << " external run manager calling forward run command(s)" << endl;

#ifdef OS_WIN
	//a flag to track if the run was terminated
	bool term_break = false;
	//create a job object to track child and grandchild process
	HANDLE job = CreateJobObject(NULL, NULL);
	if (job == NULL) throw PestError("could not create job object handle");
	JOBOBJECT_EXTENDED_LIMIT_INFORMATION jeli = { 0 };
	jeli.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE;
	if (0 == SetInformationJobObject(job, JobObjectExtendedLimitInformation, &jeli, sizeof(jeli)))
	{
		throw PestError("could not assign job limit flag to job object");
	}
	for (auto &cmd_string : comline_vec)
	{

	    cout << pest_utils::get_time_string() << " calling forward run command: '" << cmd_string << "' " << endl;
		//start the command
		PROCESS_INFORMATION pi;
		try
		{
			pi = start(cmd_string);
		}
		catch (...)
		{
			throw std::runtime_error("start_command() failed for command: " + cmd_string);
		}
		if (0 == AssignProcessToJobObject(job, pi.hProcess))
		{
			throw PestError("could not add process to job object: " + cmd_string);
		}
		DWORD pid = pi.dwProcessId;
        cout << "...pid: " << pid << endl;
		DWORD exitcode;
		while (true)
		{
            //sleep
			std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
			//check if process is still active
			GetExitCodeProcess(pi.hProcess, &exitcode);
			//if the process ended, break
			if (exitcode != STILL_ACTIVE)
			{
				if (exitcode != 0)
				{
					cout << "exit_code: " << exitcode << endl;
					throw std::runtime_error("GetExitCodeProcess() returned error status for command: " + cmd_string);
				}
                cout << "...exit_code: " << exitcode << endl;
                break;
			}

		}
	}


#endif

#ifdef OS_LINUX
	//a flag to track if the run was terminated

	for (auto &cmd_string : comline_vec)
	{

        cout << pest_utils::get_time_string() << " calling forward run command: '" << cmd_string << "' " << endl;
		//start the command
		int command_pid = start(cmd_string);
		cout << "...pid: " << command_pid << endl;
		while (true)
		{
            std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
            //check if process is still active
			int status = 0;
			pid_t exit_code = waitpid(command_pid, &status, WNOHANG);

            //if the process ended, break
			if ((exit_code == -1) || (status != 0))
			{
                cout << "exit_code: " << exit_code << endl;
                cout << "status: " << status << endl;
				throw std::runtime_error("waitpid() returned error status for command: " + cmd_string);
			}
			else if (exit_code != 0)
			{

		        cout << "...exit_code: " << exit_code << endl;
		        cout << "...status: " << status << endl;
			    break;
			}


		}

	}
#endif

	cout << pest_utils::get_time_string() << " forward run command(s) finished, took " << pest_utils::get_duration_sec(start_time) << " seconds" << endl;
	//cout << pest_utils::get_time_string() << " re-initializing run storage file";
	//file_stor.init_restart(file_stor.get_filename());

}


RunManagerExternal::~RunManagerExternal()
{
}

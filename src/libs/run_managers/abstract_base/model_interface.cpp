#include "network_wrapper.h"
#include "network_package.h"
#include "utilities.h"
#include "system_variables.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <sstream>
#include <thread>
#include "model_interface.h"

using namespace std;

extern "C"
{
	void mio_initialise_w_(int *, int *, int *, int *, int *);
	void mio_put_file_w_(int *, int *, int *, char *, long *);
	void mio_get_file_w_(int *, int *, int *, char *);
	void mio_store_instruction_set_w_(int *);
	void mio_process_template_files_w_(int *, int *, char *);
	void mio_delete_output_files_w_(int *, char *);
	void mio_write_model_input_files_w_(int *, int *, char *, double *);
	void mio_read_model_output_files_w_(int *, int *, char *, double *);
	void mio_finalise_w_(int *);
	void mio_get_status_w_(int *, int *);
	void mio_get_dimensions_w_(int *, int *);
	void mio_get_message_string_w_(int *, int *, char *);

}


void ModelInterface::throw_mio_error(string base_message)
{
	int mess_len = 500;
	char message[500];
	int nerr_len = 500;
	char err_instruct[500];
	for (int i = 0; i < 500; i++)
		err_instruct[i] = ' ';
	//cout << endl << endl << " MODEL INTERFACE ERROR:" << endl;
	mio_get_message_string_w_(&ifail, &mess_len, message);
	string err = string(message);
	auto s_end = err.find_last_not_of(" \t", 500);
	err = err.substr(0, s_end);
	throw runtime_error("model input/output error:" + base_message + "\n" + err);
}


void ModelInterface::set_files()
{
	//put template files
	int inum = 1;
	int itype = 1;
	for (auto &file : tplfile_vec)
	{
		long f_name_len = 180;
		vector<char> f_name = pest_utils::string_as_fortran_char_ptr(file, f_name_len);
		mio_put_file_w_(&ifail, &itype, &inum, f_name.data(), &f_name_len);
		if (ifail != 0) throw_mio_error("putting template file" + file);
		inum++;
	}

	//put model in files
	inum = 1;
	itype = 2;
	for (auto &file : inpfile_vec)
	{
		long f_name_len = 180;
		vector<char> f_name = pest_utils::string_as_fortran_char_ptr(file, f_name_len);
		mio_put_file_w_(&ifail, &itype, &inum, f_name.data(), &f_name_len);
		if (ifail != 0) throw_mio_error("putting model input file" + file);
		inum++;
	}

	//put instructions files
	inum = 1;
	itype = 3;
	for (auto &file : insfile_vec)
	{
		long f_name_len = 180;
		vector<char> f_name = pest_utils::string_as_fortran_char_ptr(file, f_name_len);
		mio_put_file_w_(&ifail, &itype, &inum, f_name.data(), &f_name_len);
		if (ifail != 0) throw_mio_error("putting instruction file" + file);
		inum++;
	}

	//put model out files
	inum = 1;
	itype = 4;
	for (auto &file : outfile_vec)
	{
		long f_name_len = 180;
		vector<char> f_name = pest_utils::string_as_fortran_char_ptr(file, f_name_len);
		mio_put_file_w_(&ifail, &itype, &inum, f_name.data(), &f_name_len);
		if (ifail != 0) throw_mio_error("putting model output file" + file);
		inum++;
	}
}

ModelInterface::ModelInterface()
{
	initialized = false;
}

ModelInterface::ModelInterface(vector<string> _tplfile_vec, vector<string> _inpfile_vec,
	vector<string> _insfile_vec, vector<string> _outfile_vec, vector<string> _comline_vec)
{
	tplfile_vec = _tplfile_vec;
	inpfile_vec = _inpfile_vec;
	insfile_vec = _insfile_vec;
	outfile_vec = _outfile_vec;
	comline_vec = _comline_vec;

	initialized = false;
}

void ModelInterface::initialize(vector<string> _tplfile_vec, vector<string> _inpfile_vec,
	vector<string> _insfile_vec, vector<string> _outfile_vec, vector<string> _comline_vec,
	vector<string> &_par_name_vec, vector<string> &_obs_name_vec)
{
	tplfile_vec = _tplfile_vec;
	inpfile_vec = _inpfile_vec;
	insfile_vec = _insfile_vec;
	outfile_vec = _outfile_vec;
	comline_vec = _comline_vec;

	initialize(_par_name_vec,_obs_name_vec);
}


void ModelInterface::initialize(vector<string> &_par_name_vec, vector<string> &_obs_name_vec)
{
	par_name_vec = _par_name_vec;
	obs_name_vec = _obs_name_vec;
	int npar = par_name_vec.size();
	int nobs = obs_name_vec.size();
	int ntpl = tplfile_vec.size();
	int nins = insfile_vec.size();

	if (ntpl <= 0)
		throw runtime_error("number of template files <= 0");
	if (nins <= 0)
		throw runtime_error("number of instructino files <=0");

	mio_initialise_w_(&ifail, &ntpl, &nins, &npar, &nobs);
	if (ifail != 0) throw_mio_error("initializing mio module");

	set_files();

	//check template files
	mio_process_template_files_w_(&ifail, &npar, pest_utils::StringvecFortranCharArray(par_name_vec, 200, pest_utils::TO_LOWER).get_prt());
	if (ifail != 0)throw_mio_error("error in template files");

	////build instruction set
	mio_store_instruction_set_w_(&ifail);
	if (ifail != 0) throw_mio_error("error building instruction set");

	initialized = true;

}

void ModelInterface::finalize()
{
	mio_finalise_w_(&ifail);
	if (ifail != 0) ModelInterface::throw_mio_error("error finalizing model interface");
	initialized = false;
}

ModelInterface::~ModelInterface()
{
	finalize();

}

void ModelInterface::run(Parameters* pars, Observations* obs)
{

	pest_utils::thread_flag terminate(false);
	pest_utils::thread_flag finished(false);
	pest_utils::thread_exceptions shared_exceptions;



	run(&terminate, &finished, &shared_exceptions, pars, obs);
	if (shared_exceptions.size() > 0)
	{
		finalize();
		shared_exceptions.rethrow();
	}

}


void ModelInterface::run(pest_utils::thread_flag* terminate, pest_utils::thread_flag* finished, pest_utils::thread_exceptions *shared_execptions,
						Parameters* pars, Observations* obs)
{



	if (!initialized)
	{
		vector<string> pnames = pars->get_keys();
		vector<string> onames = obs->get_keys();
		initialize(pnames, onames);
	}
	//get par vals that are aligned with this::par_name_vec since the mio module was initialized with this::par_name_vec order
	par_vals = pars->get_data_vec(par_name_vec);

	try
	{
		//first delete any existing input and output files
		// This outer loop is a work around for a bug in windows.  Window can fail to release a file
		// handle quick enough when the external run executes very quickly
		bool failed_file_op = true;
		int n_tries = 0;
		while (failed_file_op)
		{
			vector<string> failed_file_vec;
			failed_file_op = false;
			for (auto &out_file : outfile_vec)
			{
				if ((pest_utils::check_exist_out(out_file)) && (remove(out_file.c_str()) != 0))
				{
					failed_file_vec.push_back(out_file);
					failed_file_op = true;
				}
			}
			for (auto &in_file : inpfile_vec)
			{
				if ((pest_utils::check_exist_out(in_file)) && (remove(in_file.c_str()) != 0))
				{
					failed_file_vec.push_back(in_file);
					failed_file_op = true;
				}
			}
			if (failed_file_op)
			{
				++n_tries;
				w_sleep(1000);
				if (n_tries > 5)
				{
					ostringstream str;
					str << "model interface error: Cannot delete existing following model files:";
					for (const string &ifile : failed_file_vec)
					{
						str << " " << ifile;
					}
					throw PestError(str.str());
				}
			}

		}

		//check for nans in par vals before continuing
		// vector<string> invalid;

		// vector<double> ivals;
		// for (int i = 0; i != par_name_vec.size(); i++)
		// {
		// 	if (OperSys::double_is_invalid(par_vals.at(i)))
		// 	{
		// 		invalid.push_back(par_name_vec.at(i));
		// 		ivals.push_back(par_vals.at(i));
		// 	}


		// for (int i = 0; i != par_name_vec.size(); i++)
		// {
		// 	if (OperSys::double_is_invalid(par_vals.at(i)))
		// 		invalid.push_back(par_name_vec.at(i));
		// }
		// if (invalid.size() > 0)
		// {
		// 	stringstream ss;
		// 	ss << "internal PEST++ error: invalid parameter values passed to model_interface for the following parameters: ";
		// 	for (auto &i : invalid)
		// 		ss << i << '\n';
		// 	for (auto &iv : ivals)
		// 		ss << iv << "\n";
		// 	throw PestError(ss.str());
		// }

		int npar = par_vals.size();
		try
		{
			mio_write_model_input_files_w_(&ifail, &npar,
				pest_utils::StringvecFortranCharArray(par_name_vec, 200, pest_utils::TO_LOWER).get_prt(),
				&par_vals[0]);
		}
		catch (exception &e)
		{
			string emess = e.what();
			throw_mio_error("uncaught error writing model input files from template files:" + emess);
		}
		if (ifail != 0) throw_mio_error("error writing model input files from template files");


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
			//start the command
			PROCESS_INFORMATION pi;
			try
			{
				pi = start(cmd_string);
			}
			catch (...)
			{
				finished->set(true);
				throw std::runtime_error("start_command() failed for command: " + cmd_string);
			}
			if (0 == AssignProcessToJobObject(job, pi.hProcess))
			{
				throw PestError("could not add process to job object: " + cmd_string);
			}
			DWORD exitcode;
			while (true)
			{
				//sleep
				std::this_thread::sleep_for(std::chrono::milliseconds(OperSys::thread_sleep_milli_secs));
				//check if process is still active
				GetExitCodeProcess(pi.hProcess, &exitcode);
				//if the process ended, break
				if (exitcode != STILL_ACTIVE)
				{
					break;
				}
				//else cout << exitcode << "...still waiting for command " << cmd_string << endl;
				//check for termination flag
				if (terminate->get())
				{
					std::cout << "received terminate signal" << std::endl;
					//try to kill the process
					bool success = (CloseHandle(job) != 0);

					//bool success = TerminateProcess(pi.hProcess, 0);
					if (!success)
					{
						finished->set(true);
						throw std::runtime_error("unable to terminate process for command: " + cmd_string);
					}
					term_break = true;

					break;
				}
			}
			//jump out of the for loop if terminated
			if (term_break) break;
		}


#endif

#ifdef OS_LINUX
		//a flag to track if the run was terminated
		bool term_break = false;
		for (auto &cmd_string : comline_vec)
		{
			//start the command
			int command_pid = start(cmd_string);
			while (true)
			{
				//sleep
				std::this_thread::sleep_for(std::chrono::milliseconds(OperSys::thread_sleep_milli_secs));
				//check if process is still active
				int status;
				pid_t exit_code = waitpid(command_pid, &status, WNOHANG);
				//if the process ended, break
				if (exit_code == -1)
				{
					finished->set(true);
					throw std::runtime_error("waitpid() returned error status for command: " + cmd_string);
				}
				else if (exit_code != 0)
				{
					break;
				}
				//check for termination flag
				if (terminate->get())
				{
					std::cout << "received terminate signal" << std::endl;
					//try to kill the process
					errno = 0;
					int success = kill(-command_pid, SIGKILL);
					if (success == -1)
					{
						finished->set(true);
						throw std::runtime_error("unable to terminate process for command: " + cmd_string);
					}
					term_break = true;
					break;
				}
			}
			//jump out of the for loop if terminated
			if (term_break) break;
		}
#endif

		if (term_break) return;

		// process instruction files
		int nins = insfile_vec.size();
		int nobs = obs_name_vec.size();
		obs_vals.resize(nobs, -9999.00);
		/*int nerr_len = 500;
		char err_instruct[500];
		for (int i = 0; i < 500; i++)
			err_instruct[i] = '|';*/
		try {
			mio_read_model_output_files_w_(&ifail, &nobs,
				pest_utils::StringvecFortranCharArray(obs_name_vec, 200, pest_utils::TO_LOWER).get_prt(),
				&obs_vals[0]);
		}
		catch (exception &e)
		{
			string emess = e.what();
			throw_mio_error("uncaught error processing model output files:" + emess);
		}
		if (ifail != 0)
		{
			/*int jfail;
			mio_get_message_string_w_(&jfail, &nerr_len, err_instruct);
			string err = string(err_instruct);
			auto s_end = err.find_last_not_of(' ',500);
			err = err.substr(0, s_end);*/

			throw_mio_error("error processing model output files");
		}

		// invalid.clear();
		// for (int i = 0; i != par_name_vec.size(); i++)
		// {
		// 	if (OperSys::double_is_invalid(par_vals.at(i)))
		// 		invalid.push_back(par_name_vec.at(i));
		// }
		// if (invalid.size() > 0)
		// {
		// 	stringstream ss;
		// 	ss << "invalid parameter values read for the following parameters: ";
		// 	for (auto &i : invalid)
		// 		ss << i << '\n';
		// 	throw PestError(ss.str());
		// }

		// for (int i = 0; i != obs_name_vec.size(); i++)
		// {
		// 	if (OperSys::double_is_invalid(obs_vals.at(i)))
		// 		invalid.push_back(obs_name_vec.at(i));
		// }
		// if (invalid.size() > 0)
		// {
		// 	stringstream ss;
		// 	ss << "invalid observation values read for the following observations: ";
		// 	for (auto &i : invalid)
		// 		ss << i << '\n';
		// 	throw PestError(ss.str());
		// }

		pars->update(par_name_vec, par_vals);
		obs->update(obs_name_vec, obs_vals);


		//set the finished flag for the listener thread
		finished->set(true);

	}
	catch (...)
	{
		shared_execptions->add(current_exception());
	}
	return;

}




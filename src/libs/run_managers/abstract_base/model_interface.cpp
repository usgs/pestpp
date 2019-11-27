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

	// old fortran processing routines...
	//mio_initialise_w_(&ifail, &ntpl, &nins, &npar, &nobs);
	//if (ifail != 0) throw_mio_error("initializing mio module");
	//set_files();
	//check template files
	//mio_process_template_files_w_(&ifail, &npar, pest_utils::StringvecFortranCharArray(par_name_vec, 200, pest_utils::TO_LOWER).get_prt());
	//if (ifail != 0)throw_mio_error("error in template files");

	for (auto tplfile : tplfile_vec)
	{
		TemplateFile tpl(tplfile);
		//tpl.parse_and_check();
		templatefiles.push_back(tpl);
	}

	for (auto insfile : insfile_vec)
	{
		InstructionFile ins(insfile);
		instructionfiles.push_back(ins);
	}

	////build instruction set
	//mio_store_instruction_set_w_(&ifail);
	//if (ifail != 0) throw_mio_error("error building instruction set");

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
	//par_vals = pars->get_data_vec(par_name_vec);
	//map<string, double> obs_map, isf_obs_map,par_map,tf_par_map;
	Observations pro_obs;
	vector<Parameters> pro_par_vec;
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

		//this is the old fortran tpl processing
		/*int npar = par_vals.size();
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
        */
		
		for (int i = 0; i < templatefiles.size(); i++)
			pro_par_vec.push_back(templatefiles[i].write_input_file(inpfile_vec[i], *pars));
		//update pars to account for possibly truncated par values...important for jco calcs
		for (auto pro_pars : pro_par_vec)
			pars->update_without_clear(pro_pars.get_keys(), pro_pars.get_data_vec(pro_pars.get_keys()));


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
	
		//this is the old fortran ins processing
		/*try {
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

			throw_mio_error("error processing model output files");
		}
		*/

		for (int i = 0; i < insfile_vec.size(); i++)
		{
			pro_obs = instructionfiles[i].read_output_file(outfile_vec[i]);
			obs->update_without_clear(pro_obs.get_keys(), pro_obs.get_data_vec(pro_obs.get_keys()));
		}

		//todo: return updated (possibly truncated) par values
		/*pars->update(par_name_vec, par_vals);
		obs->update(obs_name_vec, obs_vals);
		for (auto o : obs_map)
			obs->update_rec(o.first, o.second);
		for (auto p : par_map)
			pars->update_rec(p.first, p.second);*/

		//set the finished flag for the listener thread
		finished->set(true);

	}
	catch (...)
	{
		shared_execptions->add(current_exception());
	}
	return;

}

set<string> TemplateFile::parse_and_check()
{
	ifstream f(tpl_filename);
	prep_tpl_file_for_reading(f);
	return get_names(f);

}

Parameters TemplateFile::write_input_file(const string& input_filename, Parameters& pars)
{
	ifstream f_tpl(tpl_filename);
	prep_tpl_file_for_reading(f_tpl);
	ofstream f_in(input_filename);
	if (f_in.bad())
		throw_tpl_error("couldn't open model input file '" + input_filename + "' for writing");
	string line, val_str, name;
	double val;
	map<string, pair<int, int>> tpl_line_map;
	Parameters pro_pars;
	while (true)
	{
		if (f_tpl.eof())
			break;
		line = read_line(f_tpl);
		tpl_line_map = parse_tpl_line(line);
		
		for (auto t : tpl_line_map)
		{
			name = t.first;
			
			/*val = 1.23456789123456789123456789E+100;
			val_str = cast_to_fixed_len_string(f_rec, 200, val, name);

			val = 1.23456789123456789123456789E+100;
			val_str = cast_to_fixed_len_string(f_rec, 8, val, name);

			val = 1.23456789123456789123456789E-100;
			val_str = cast_to_fixed_len_string(f_rec, 8, val, name);

			val = -1.23456789123456789123456789E+100;
			val_str = cast_to_fixed_len_string(f_rec, 9, val, name);

			val = -1.23456789123456789123456789E-100;
			val_str = cast_to_fixed_len_string(f_rec, 9, val, name);

			val = 1.23456789123456789123456789E+10;
			val_str = cast_to_fixed_len_string(f_rec, 7, val, name);
				
			val = 1.23456789123456789123456789E-10;
			val_str = cast_to_fixed_len_string(f_rec, 7, val, name);

			val = 1.23456789123456789123456789;
			val_str = cast_to_fixed_len_string(f_rec, 1, val, name);

			val = -1.23456789123456789123456789;
			val_str = cast_to_fixed_len_string(f_rec, 2, val, name);*/

			val = pars.get_rec(t.first);
			val_str = cast_to_fixed_len_string(t.second.second, val, name);
			line.replace(t.second.first, t.second.second, val_str);
			pest_utils::convert_ip(val_str, val);
			pro_pars.insert(name, val);
		}
		f_in << line << endl;
	}
	return pro_pars;
}

void TemplateFile::prep_tpl_file_for_reading(ifstream& f_tpl)
{
	if (f_tpl.bad())
	{
		throw_tpl_error("couldn't open tpl file for reading");
	}
	string tag, line;
	vector<string> tokens;
	line = read_line(f_tpl);
	pest_utils::tokenize(line, tokens);
	if (tokens.size() < 2)
		throw_tpl_error("incorrect first line - expecting 'ptf <marker>'", line_num);
	if (tokens.size() > 2)
		throw_tpl_error("extra unused items on first line");
	tag = pest_utils::upper_cp(tokens[0]);
	if (tag != "PTF")
		throw_tpl_error("first line should start with 'PTF', not: " + tag);
	marker = tokens[1];
	if (marker.size() != 1)
		throw_tpl_error("marker on first line should be one character, not: " + marker);
}

set<string> TemplateFile::get_names(ifstream& f)
{
	set<string> names;
	string line;
	map<string, pair<int, int>> tpl_line_map;
	
	while (true)
	{
		if (f.eof())
			break;
		line = read_line(f);
		tpl_line_map = parse_tpl_line(line);
		for (auto t : tpl_line_map)
			names.insert(t.first);
	}
	return names;
}

vector<int> TemplateFile::find_all_marker_indices(const string& line)
{
	vector<int> indices;
	int pos = line.find(marker);
	while (pos != string::npos)
	{
		indices.push_back(pos);
		pos = line.find(marker, pos + marker.size());
	}
	return indices;
}

void TemplateFile::throw_tpl_error(const string& message, int lnum , bool warn)
{
	stringstream ss;
	if (warn)
		ss << "TemplateFile warning in " << tpl_filename;
	else
		ss << "TemplateFile error in " << tpl_filename;
	if (lnum != 0)
		ss << "on line: " << lnum;
	ss <<" : " << message;
	if (warn)
		cout << endl << ss.str() << endl;
	else
		throw runtime_error(ss.str());
}

map<string, pair<int, int>> TemplateFile::parse_tpl_line(const string& line)
{
	vector<int> indices = find_all_marker_indices(line);
	if (indices.size() % 2 != 0)
		throw_tpl_error("unbalanced marker ('" + marker + "') ", line_num);
	int s, e, len;
	string name;
	pair<int, int> se_idx;
	map<string, pair<int, int>> tpl_line_map;
	for (int i = 0; i < indices.size(); i = i + 2)
	{
		s = indices[i];
		e = indices[i + 1];
		len = (e - s) + 1;
		name = line.substr(s+1, len-2);
		pest_utils::upper_ip(name);
		pest_utils::strip_ip(name);
		tpl_line_map[name] = pair<int, int>(s, len);
	}
	return tpl_line_map;
}

string TemplateFile::cast_to_fixed_len_string(int size, double value, string& name)
{
	string val_str;
	stringstream ss;
	int precision = size;
	if (value < 0)
		precision--; // for the minus sign
	if ((abs(value) >= 10) || (abs(value) <= 1.0))
	{
		ss << scientific;
		precision = precision - 2; //for the "e" and (at least) 1 exponent digit
	}
	ss.width(size);
	
	ss.fill('0');
	while (true)
	{
		
		ss.str("");
		ss.precision(precision);
		ss << value;
		val_str = ss.str();
		if (val_str.size() <= size)
			break;
		if (val_str.size() > size)
			precision--;
		if (precision <= 0)
		{
			ss.str("");
			ss << "TemplateFile casting error: cant represent value " << value;
			ss << " for " << name << " in space that is only " << size << " chars wide";
			throw_tpl_error(ss.str());
		}
		
	}
	//occasionally, when reducing precision, rounding will cause an 
	// extra char to be dropped, so this left pads it back
	if (val_str.size() < size) 
	{
		ss.str("");
		for (int i = val_str.size(); i < size; i++)
			ss << "0";
		ss << val_str;
		val_str = ss.str();
	}
	
	/*int width = size;
	if (value < 0.0)
		width--;
		
	ss << value;
	val_str = ss.str();*/
	if (val_str.size() != size)
		throw_tpl_error("val_str != size: " + val_str);
	return val_str;
}

string TemplateFile::read_line( ifstream& f_tpl)
{
	if (f_tpl.bad())
		throw_tpl_error("cant read next line", line_num);
	string line;
	if (f_tpl.eof())
		throw_tpl_error("unexpected eof", line_num);
	
	getline(f_tpl, line);
	line_num++;
	return line;
}

string InstructionFile::read_ins_line(ifstream& f_ins)
{
	if (f_ins.bad())
		throw_ins_error("cant read next instruction file line", ins_line_num);
	string line;
	if (f_ins.eof())
		throw_ins_error("unexpected instruction file eof ", ins_line_num);

	getline(f_ins, line);
	last_ins_line = line;
	ins_line_num++;
	return line;
}


string InstructionFile::read_out_line(ifstream& f_out)
{
	if (f_out.bad())
		throw_ins_error("cant read next output file line", ins_line_num, out_line_num);
	string line;
	if (f_out.eof())
		throw_ins_error("unexpected output file eof ", ins_line_num, out_line_num);
	getline(f_out, line);
	last_out_line = line;
	out_line_num++;
	return line;
}


InstructionFile::InstructionFile(string _ins_filename): ins_filename(_ins_filename), ins_line_num(0),
out_line_num(0),last_ins_line(""),last_out_line("")
{
	obs_tags.push_back(pair<char, char>('{', '}'));
	obs_tags.push_back(pair<char, char>('[', ']'));
	obs_tags.push_back(pair<char, char>('!', '!'));	
}


set<string> InstructionFile::parse_and_check()
{
	ifstream f_ins(ins_filename);
	prep_ins_file_for_reading(f_ins);
	string line, name;
	vector<string> tokens, s_names;
	int spos,epos;
	while (true)
	{
		if (f_ins.eof())
			break;
		line = read_ins_line(f_ins);
		pest_utils::upper_ip(line);
		pest_utils::tokenize(line, tokens);
		//for (auto token : tokens)
		for (int i=0;i<tokens.size();i++)
		{
			parse_obs_name_from_token(tokens[i], name);
			if (name.size() > 0)
			{
				s_names.push_back(name);
				
			}
		}
	}
	f_ins.close();
	names.insert(s_names.begin(), s_names.end());
	return names;
}

void InstructionFile::prep_ins_file_for_reading(ifstream& f_ins)
{
	if (f_ins.bad())
	{
		throw_ins_error("couldn't open ins file for reading");
	}
	string tag, line;
	vector<string> tokens;
	line = read_ins_line(f_ins);
	pest_utils::tokenize(line, tokens);
	if (tokens.size() < 2)
		throw_ins_error("incorrect first line - expecting 'pif <marker>'", ins_line_num);
	if (tokens.size() > 2)
		throw_ins_error("extra unused items on first line");
	tag = pest_utils::upper_cp(tokens[0]);
	if (tag != "PIF")
		throw_ins_error("first line should start with 'PIF', not: " + tag);
	string s_marker = tokens[1];
	if (s_marker.size() != 1)
		throw_ins_error("marker on first line should be one character, not: " + s_marker);
	marker = s_marker.c_str()[0];
}


Observations InstructionFile::read_output_file(const string& output_filename)
{
	ifstream f_ins(ins_filename);
	ifstream f_out(output_filename);
	prep_ins_file_for_reading(f_ins);
	if (f_out.bad())
	{
		throw_ins_error("can't open output file'" + output_filename + "' for reading");
	}
	string ins_line, out_line;
	vector<string> tokens;
	Observations obs;
	pair<string, double> lhs;
	while (true)
	{

		if (f_ins.eof())
			break;
		tokens.clear();
		ins_line = read_ins_line(f_ins);
		pest_utils::upper_ip(ins_line);
		pest_utils::tokenize(ins_line, tokens);
		for (auto token : tokens)
		{

			if (token[0] == 'L')
			{
				execute_line_advance(token, out_line, f_out);
			}
			else if (token[0] == 'W')
			{
				execute_whitespace(token, out_line, f_out);
			}
			else if (token[0] == '[')
			{
				lhs = execute_fixed(token, out_line, f_out);
				if (lhs.first != "DUM")
					obs.insert(lhs.first,lhs.second);
			}
			else if (token[0] == '!')
			{
				lhs = execute_free(token, out_line, f_out);
				if (lhs.first != "DUM")
					obs.insert(lhs.first, lhs.second);
			}
			else if (token[0] == '(')
			{
				lhs = execute_free(token, out_line, f_out);
				if (lhs.first != "DUM")
					obs.insert(lhs.first, lhs.second);
			}
			else if (token[0] == marker)
			{
				//if this is the first instruction, its a primary search
				if (token == tokens[0])
				{
					execute_primary(token, out_line, f_out);
				}
				else
				{
					execute_secondary(token, out_line, f_out);
				}
			}
		}
	}
	return obs;	
}


void InstructionFile::throw_ins_error(const string& message, int ins_lnum, int out_lnum, bool warn)
{
	stringstream ss;
	if (warn)
		ss << "InstructionFile warning in " << ins_filename;
	else
		ss << "InstructionFile error in " << ins_filename;
	if (ins_lnum != 0)
		ss << "on instruction file line: " << ins_lnum;
	if (out_lnum != 0)
		ss << "on output file line: " << out_lnum;
	ss << " : " << message;
	if (warn)
		cout << endl << ss.str() << endl;
	else
		throw runtime_error(ss.str());
}

void InstructionFile::parse_obs_name_from_token(const string& token, string& name)
{
	int spos, epos;
	//for (auto ot : obs_tags)
	for (int i=0;i<obs_tags.size();i++)
	{
		if (token[0] == obs_tags[i].first)
		{
			if (token.find(obs_tags[i].second) == string::npos)
				throw_ins_error("unbalanced obs tag'" + string(1, obs_tags[i].first) + "'", ins_line_num);
			spos = token.find(obs_tags[i].first);
			epos = token.find(obs_tags[i].second, spos + 1);
			name = pest_utils::strip_cp(token.substr(spos + 1, (epos - spos) - 1));
			
			//break;
		}
	}
	//return name;
}

pair<string, double> InstructionFile::execute_fixed(const string& token, string& line, ifstream& f_out)
{
	return pair<string, double>();
}

pair<string, double> InstructionFile::execute_semi(const string& token, string& line, ifstream& f_out)
{
	return pair<string, double>();
}

pair<string, double> InstructionFile::execute_free(const string& token, string& line, ifstream& f_out)
{
	vector<string> tokens;
	pest_utils::tokenize(line, tokens);
	if (tokens.size() == 0)
		throw_ins_error("error tokenizing output line ('"+last_out_line+"') for whitespace instruction execution on line: "+last_ins_line, ins_line_num, out_line_num);
	double value;
	try
	{
		pest_utils::convert_ip(tokens[0], value);
	}
	catch (...)
	{
		throw_ins_error("error converting '" + tokens[0] + "' to double on output line '" + last_out_line + "'", ins_line_num, out_line_num);
	}
	string name = token.substr(1, token.size() - 2);
	int pos = line.find(tokens[0]);
	if (pos == string::npos)
	{
		throw_ins_error("internal error: could not find whitespace token", ins_line_num, out_line_num);
	}
	line = line.substr(pos + tokens[0].size());

	return pair<string, double>(name,value);
}

void InstructionFile::execute_primary(const string& token, string& line, ifstream& f_out)
{
	//check that a closing marker is found
	int pos = token.substr(1).find(marker);
	if (pos == string::npos)
		throw_ins_error("primary marker token '" + token + "' doesn't have a closing marker char", ins_line_num);

	string primary_tag = token.substr(1, pos - 1);
	while (true)
	{
		if (f_out.eof())
			throw_ins_error("EOF encountered while executing marker search ('" + token + "')", ins_line_num, out_line_num);
		line = read_out_line(f_out);
		pos = line.find(primary_tag);
		if (pos != string::npos)
		{
			break;
		}
	}
	pos = pos + primary_tag.size();
	line = line.substr(pos);
	return;
}


void InstructionFile::execute_secondary(const string& token, string& line, ifstream& f_out)
{
	//check that a closing marker is found
	int pos = token.substr(1).find(marker);
	if (pos == string::npos)
		throw_ins_error("secondary marker token '" + token + "' doesnt have a closing marker char");
	string secondary_tag = token.substr(1, pos-1);
	pos = line.find(secondary_tag);
	if (pos == string::npos)
	{
		throw_ins_error("EOL encountered while executing secondary marker ('" + secondary_tag + "') search on output line", ins_line_num,out_line_num);
	}
	line = line.substr(pos + secondary_tag.size());
}


void InstructionFile::execute_whitespace(const string& token, string& line, ifstream& f_out)
{
	int pos = line.find_first_of(" \t");
	if (pos == string::npos)
	{
		throw_ins_error("EOL encountered while executing whitespace instruction on output line", ins_line_num, out_line_num);
	}
	line = line.substr(pos + 1);
}


void InstructionFile::execute_line_advance(const string& token, string& line, ifstream& f_out)
{
	stringstream ss;
	int num;
	pest_utils::convert_ip(token.substr(1), num);
	for (int i = 0; i < num; i++)
	{
		if (f_out.bad())
		{	
			throw_ins_error("'bad' stream when executing line advance instruction", ins_line_num, out_line_num);
		}
		if (f_out.eof())
		{
			throw_ins_error("EOF encountered when executing line advance instruction", ins_line_num, out_line_num);
		}
		line = read_out_line(f_out);
	}
}

#include "network_wrapper.h"
#include "network_package.h"
#include "utilities.h"
#include "system_variables.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstring>
#include <sstream>
#include <thread>
#include <unordered_set>
#include "model_interface.h"

using namespace std;

//extern "C"
//{
//	void mio_initialise_w_(int *, int *, int *, int *, int *);
//	void mio_put_file_w_(int *, int *, int *, char *, long *);
//	void mio_get_file_w_(int *, int *, int *, char *);
//	void mio_store_instruction_set_w_(int *);
//	void mio_process_template_files_w_(int *, int *, char *);
//	void mio_delete_output_files_w_(int *, char *);
//	void mio_write_model_input_files_w_(int *, int *, char *, double *);
//	void mio_read_model_output_files_w_(int *, int *, char *, double *);
//	void mio_finalise_w_(int *);
//	void mio_get_status_w_(int *, int *);
//	void mio_get_dimensions_w_(int *, int *);
//	void mio_get_message_string_w_(int *, int *, char *);
//
//}


void ModelInterface::throw_mio_error(string base_message)
{
	throw runtime_error("model input/output error:" + base_message);
}



void ModelInterface::check_io_access()
{
	
	if (tplfile_vec.size() == 0)
	{
		throw_mio_error("number of template files = 0");
	}
	if (insfile_vec.size() == 0)
	{
		throw_mio_error("number of instruction files = 0");
	}
	vector<string> inaccessible_files;
	for (auto& file : insfile_vec)
		if (!pest_utils::check_exist_in(file)) inaccessible_files.push_back(file);
	for (auto& file : outfile_vec)
		if (!pest_utils::check_exist_out(file)) inaccessible_files.push_back(file);
	for (auto& file : tplfile_vec)
		if (!pest_utils::check_exist_in(file)) inaccessible_files.push_back(file);
	for (auto& file : inpfile_vec)
		if (!pest_utils::check_exist_out(file)) inaccessible_files.push_back(file);
	
	if (inaccessible_files.size() != 0)
	{
		string missing;
		for (auto& file : inaccessible_files)
			missing += file + " , ";
		cout << "Could not access the following model interface files: " << missing;
		throw PestError("Could not access the following model interface files: " + missing);
		
	}
	if (inaccessible_files.size() != 0)
	{
		string missing;
		for (auto& file : inaccessible_files)
			missing += file + " , ";
		throw PestError("Could not access the following model interface files: " + missing);
	}
}

void ModelInterface::check_tplins(const vector<string> &par_names, const vector<string> &obs_names)
{	
	//rigorous checking of names in tpl and ins files vs control file
	unordered_set<string> ins_obs_names, file_obs_names;
	for (auto ins_file : insfile_vec)
	{
		InstructionFile isf(ins_file);
		file_obs_names = isf.parse_and_check();
		ins_obs_names.insert(file_obs_names.begin(), file_obs_names.end());
		//isf.read_output_file(model_exec_info.outfile_vec[0]);
	}
	unordered_set<string> pst_obs_names, diff;
	pst_obs_names.insert(obs_names.begin(), obs_names.end());
	unordered_set<string>::iterator end = ins_obs_names.end();
	for (auto p : pst_obs_names)
		if (ins_obs_names.find(p) == end)
			diff.insert(p);
	if (diff.size() > 0)
	{
		stringstream ss;
		ss << "Error: the following observations were found in the control file but not in the instruction files:" << endl;
		for (auto d : diff)
			ss << d << endl;
		throw_mio_error(ss.str());
	}
	end = pst_obs_names.end();
	for (auto p : ins_obs_names)
		if (pst_obs_names.find(p) == end)
			diff.insert(p);
	if (diff.size() > 0)
	{
		stringstream ss;
		ss << "Error: the following observations were found in the instruction files but not in the control file:" << endl;
		for (auto d : diff)
			ss << d << endl;
		throw_mio_error(ss.str());
	}


	unordered_set<string> tpl_par_names, file_par_names;
	for (auto tpl_file : tplfile_vec)
	{
		TemplateFile tf(tpl_file);
		file_par_names = tf.parse_and_check();
		tpl_par_names.insert(file_par_names.begin(), file_par_names.end());
	}
	unordered_set<string> pst_par_names;
	pst_par_names.insert(par_names.begin(), par_names.end());
	end = tpl_par_names.end();
	for (auto p : pst_par_names)
		if (tpl_par_names.find(p) == end)
			diff.insert(p);
	if (diff.size() > 0)
	{
		stringstream ss;
		ss << "Error: the following parameters were found in the control file but not in the template files:" << endl;
		for (auto d : diff)
			ss << d << endl;
		throw_mio_error(ss.str());
	}
	end = pst_par_names.end();
	for (auto p : tpl_par_names)
		if (pst_par_names.find(p) == end)
			diff.insert(p);
	if (diff.size() > 0)
	{
		stringstream ss;
		ss << "Error: the following parameters were found in the template files but not in the control file:" << endl;
		for (auto d : diff)
			ss << d << endl;
		throw_mio_error(ss.str());
	}
}



void ModelInterface::run(Parameters* pars, Observations* obs)
{

	pest_utils::thread_flag terminate(false);
	pest_utils::thread_flag finished(false);
	pest_utils::thread_exceptions shared_exceptions;

	run(&terminate, &finished, &shared_exceptions, pars, obs);
	if (shared_exceptions.size() > 0)
	{
		//finalize();
		shared_exceptions.rethrow();
	}

}


void ModelInterface::run(pest_utils::thread_flag* terminate, pest_utils::thread_flag* finished, pest_utils::thread_exceptions *shared_execptions,
						Parameters* pars, Observations* obs)
{

	if (templatefiles.size() == 0)
		for (auto t : tplfile_vec)
		{
			TemplateFile tt(t);
			tt.set_fill_zeros(fill_tpl_zeros);
			templatefiles.push_back(tt);
		}
			

	if (instructionfiles.size() == 0)
		for (auto i : insfile_vec)
		{
			InstructionFile ii(i);
			ii.set_additional_delimiters(additional_ins_delimiters);
			instructionfiles.push_back(ii);
		}

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
		cout << "processing tpl files...";
		vector<string> notnormal = pars->get_notnormal_keys();
		if (notnormal.size() > 0)
		{
			throw runtime_error("denormal floating point parameter values found");
		}
		for (int i = 0; i < templatefiles.size(); i++)
		{
			string name = templatefiles[i].get_tpl_filename();
			//cout << name << endl;
			pro_par_vec.push_back(templatefiles[i].write_input_file(inpfile_vec[i], *pars));
		}
		//update pars to account for possibly truncated par values...important for jco calcs
		for (auto pro_pars : pro_par_vec)
			pars->update_without_clear(pro_pars.get_keys(), pro_pars.get_data_vec(pro_pars.get_keys()));
		cout << "done" << endl;

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

		
		cout << "processing ins files...";
		Observations temp_obs;
		for (int i = 0; i < instructionfiles.size(); i++)
		{
			pro_obs = instructionfiles[i].read_output_file(outfile_vec[i]);
			temp_obs.update_without_clear(pro_obs.get_keys(), pro_obs.get_data_vec(pro_obs.get_keys()));
		}
		unordered_set<string> ins_names, pst_names;
		vector<string> t, diff;
		t = obs->get_keys();
		pst_names.insert(t.begin(), t.end());
		t = temp_obs.get_keys();
		ins_names.insert(t.begin(), t.end());
		unordered_set<string>::iterator end = ins_names.end();
		for (auto o : pst_names)
		{
			if (ins_names.find(o) == end)
				diff.push_back(o);
		}
		if (diff.size() > 0)
		{
			stringstream ss;
			ss << "ModelInterace error: the following instruction observations are not in the control file:";
			for (auto d : diff)
				ss << d << ",";
			throw_mio_error(ss.str());
		}
		end = pst_names.end();
		for (auto o : ins_names)
		{
			if (pst_names.find(o) == end)
				diff.push_back(o);
		}
		if (diff.size() > 0)
		{
			stringstream ss;
			ss << "ModelInterace error: the following control file observations are not in the instruction files:";
			for (auto d : diff)
				ss << d << ",";
			throw_mio_error(ss.str());
		}
		t = temp_obs.get_keys();
		obs->update(t, temp_obs.get_data_vec(t));
		cout << "done" << endl;

		

		//set the finished flag for the listener thread
		finished->set(true);

	}
	catch (...)
	{
		shared_execptions->add(current_exception());
	}
	return;

}

unordered_set<string> TemplateFile::parse_and_check()
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
	vector<pair<string, pair<int, int>>> tpl_line_map;
	Parameters pro_pars;
	vector<string> t = pars.get_keys();
	unordered_set<string> pnames(t.begin(), t.end());
	unordered_set<string>::iterator end = pnames.end();
	t.resize(0);
	while (true)
	{
		if (f_tpl.eof())
			break;
		line = read_line(f_tpl);
		
		if (line.size() == 0)
		{
			if (f_tpl.eof())
				break;
			f_in << endl;
			continue;
		}
		tpl_line_map = parse_tpl_line(line);
		for (auto t : tpl_line_map)
		{
			name = t.first;
			if (pnames.find(name) == end)
				throw_tpl_error("parameter '" + name + "' not listed in control file");

			/*val = 1.23456789123456789123456789E+100;
			val_str = cast_to_fixed_len_string(200, val, name);
			pest_utils::convert_ip(val_str, val);

			val = 1.23456789123456789123456789E+100;
			val_str = cast_to_fixed_len_string(8, val, name);
			pest_utils::convert_ip(val_str, val);

			val = 1.23456789123456789123456789E-100;
			val_str = cast_to_fixed_len_string(8, val, name);
			pest_utils::convert_ip(val_str, val);

			val = -1.23456789123456789123456789E+100;
			val_str = cast_to_fixed_len_string(9, val, name);
			pest_utils::convert_ip(val_str, val);

			val = -1.23456789123456789123456789E-100;
			val_str = cast_to_fixed_len_string(9, val, name);
			pest_utils::convert_ip(val_str, val);

			val = 1.23456789123456789123456789E+10;
			val_str = cast_to_fixed_len_string(7, val, name);
			pest_utils::convert_ip(val_str, val);

			val = 1.23456789123456789123456789E-10;
			val_str = cast_to_fixed_len_string(7, val, name);
			pest_utils::convert_ip(val_str, val);

			val = 1.23456789123456789123456789;
			val_str = cast_to_fixed_len_string(1, val, name);
			pest_utils::convert_ip(val_str, val);

			val = -1.23456789123456789123456789;
			val_str = cast_to_fixed_len_string(2, val, name);
			pest_utils::convert_ip(val_str, val);
			*/

			val = pars.get_rec(t.first);
			val_str = cast_to_fixed_len_string(t.second.second, val, name);
			line.replace(t.second.first, t.second.second, val_str);
			//pest_utils::convert_ip(val_str, val);
			val = stod(val_str);
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
	if ((tag != "PTF") && (tag != "JTF"))
		throw_tpl_error("first line should start with 'PTF' or 'JTF', not: " + tag);
	marker = tokens[1];
	if (marker.size() != 1)
		throw_tpl_error("marker on first line should be one character, not: " + marker);
}

unordered_set<string> TemplateFile::get_names(ifstream& f)
{
	unordered_set<string> names;
	string line;
	vector<pair<string, pair<int, int>>> tpl_line_info;
	while (true)
	{
		if (f.eof())
			break;
		line = read_line(f);
		tpl_line_info = parse_tpl_line(line);
		for (auto t : tpl_line_info)
			names.insert(t.first);
	}
	return names;
}

vector<int> TemplateFile::find_all_marker_indices(const string& line, const string& marker)
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

vector<pair<string,pair<int, int>>> TemplateFile::parse_tpl_line(const string& line)
{
	vector<int> indices = find_all_marker_indices(line, marker);
	if (indices.size() % 2 != 0)
		throw_tpl_error("unbalanced marker ('" + marker + "') ", line_num);
	int s, e, len;
	string name;
	pair<int, int> se_idx;
	pair<string, pair<int, int>> entry;
	vector<pair<string,pair<int,int>>> tpl_line_info;
	for (int i = 0; i < indices.size(); i = i + 2)
	{
		s = indices[i];
		e = indices[i + 1];
		len = (e - s) + 1;
		name = line.substr(s+1, len-2);
		pest_utils::upper_ip(name);
		pest_utils::strip_ip(name);
		//tpl_line_map[name] = pair<int, int>(s, len);
		se_idx = pair<int, int>(s, len);
		entry = pair<string, pair<int, int>>(name, se_idx);
		tpl_line_info.push_back(entry);
	}
	return tpl_line_info;
}

string TemplateFile::cast_to_fixed_len_string(int size, double value, string& name)
{
	string val_str, fill_val=" ";
	stringstream ss;
	int precision = size;
	bool sci = false;
	if (value < 0)
		precision--; // for the minus sign
	if ((abs(value) >= 100) || (abs(value) < 0.01))
	{
		ss << scientific;
		precision = precision - 2; //for the "e" and (at least) 1 exponent digit
		sci = true;
	}
	else
	{
		ss << fixed;
	}
	ss.width(size);
	
	int size_last = -1;
	if (fill_zeros)
	{
		ss.fill('0');
		ss << internal;
		fill_val = "0";
	}
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
			//time for desparate measures:
			//if the exponent has a leading zero, drop it
			if (val_str.substr(val_str.size() - 2, 1) == "0")
			{
				string t = val_str.substr(0, val_str.size() - 2);
				val_str = t + val_str.substr(val_str.size() - 1, 1);
				if (val_str.size() <= size)
					break;
			}
			//if there is an unnesscary zero(s) between the radix and the exponent
			int r_idx = val_str.find_first_of(".")+1; // to skip past the radix
			int e_idx = val_str.find_first_of("Ee");
			if (r_idx != e_idx)
			{
				string t = val_str.substr(r_idx, e_idx - r_idx);
				if (stod(t) == 0.0)
				{
					t = val_str.substr(0, r_idx-1); // to skip the radix in the new number
					val_str = t + val_str.substr(e_idx);
					if (val_str.size() <= size)
						break;
				}
			}
			ss.str("");
		 	ss << "TemplateFile casting error: cant represent value " << value;
			ss << " for " << name << " in space that is only " << size << " chars wide";
			throw_tpl_error(ss.str());
		}
		if (val_str.size() == size_last)
		{
			if (sci)
				throw_tpl_error("internal error: val_str size not decreasing over successive attempts:" + val_str);
			else
			{
				val_str = val_str.substr(0, size);
				break;
			}
			

			
		}
		size_last = val_str.size();
	}
	//occasionally, when reducing precision, rounding will cause an 
	// extra char to be dropped, so this left pads it back
	//this also pads for really large par spaces
	if (val_str.size() < size)
	{
		ss.str("");
		//if the fill value isnt a space and its a negative value
		// we need to push the dash into the stringstream, then 
		//remove the dash from the val string
		int s = size;
		if ((fill_val != " ") && (val_str.at(0) == '-'))
		{
			ss << "-";
			val_str = val_str.substr(1, val_str.size() - 1);
			s--;
		}
		for (int i = val_str.size(); i < s; i++)
			ss << fill_val;
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
	pest_utils::strip_ip(line, "\n\r");
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


InstructionFile::InstructionFile(string _ins_filename, string _addtitional_delimiters): ins_filename(_ins_filename), ins_line_num(0),
out_line_num(0),last_ins_line(""),last_out_line(""), additional_delimiters(_addtitional_delimiters)
{
	obs_tags.push_back(pair<char, char>('(', ')'));
	obs_tags.push_back(pair<char, char>('[', ']'));	
}


unordered_set<string> InstructionFile::parse_and_check()
{
	unordered_set<string> names;
	ifstream f_ins(ins_filename);
	prep_ins_file_for_reading(f_ins);
	string line, name;
	vector<string> tokens;
	int spos,epos;
	char first;
	while (true)
	{
		if (f_ins.eof())
			break;
		line = read_ins_line(f_ins);
		pest_utils::upper_ip(line);
		tokens.clear();
		pest_utils::tokenize(line, tokens);
		
		for (int i = 0; i < tokens.size(); i++)
		{
			first = tokens[i].at(0);
			if ((first == '!') || (first == '(') || (first == '['))
			{
				name = parse_obs_name_from_token(tokens[i]);
				if (!(name.find("DUM") != std::string::npos))
				{
					if (names.find(name) != names.end())
					{
						cout << line << endl;
						throw_ins_error("observation '" + name + "' listed multiple times in ins file '" + ins_filename + "'");
					}
					names.emplace(name);
				}
			}
		}
	}
	f_ins.close();
	
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
	if ((tag != "PIF") && (tag != "JIF"))
		throw_ins_error("first line should start with 'PIF' or 'JIF', not: " + tag);
	string s_marker = tokens[1];
	if (s_marker.size() != 1)
		throw_ins_error("marker on first line should be one character, not: " + s_marker);
	marker = s_marker.c_str()[0];
}


Observations InstructionFile::read_output_file(const string& output_filename)
{
	if (!pest_utils::check_exist_in(output_filename))
		throw_ins_error("output file'" + output_filename + "' not found");
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
		tokens = tokenize_ins_line(ins_line);
		//check that the first token is either a marker or a line advance
		if (tokens.size() > 0)
		{
			char first = tokens[0][0];
			if ((first != 'L') && (first != marker))
			{
				stringstream ss;
				ss << "first token on each instruction file line must be either a primary marker ";
				ss << " or a line advance instruction, not '" << tokens[0] << "'";
				throw_ins_error(ss.str());
			}
		}
		//int itoken = 0;
		bool all_markers_so_far = true;
		//for (auto token : tokens)
		for (int itoken=0;itoken<tokens.size();itoken++)
		{
			string token = tokens[itoken];

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
				all_markers_so_far = false;
			}
			else if (token[0] == '!')
			{
				lhs = execute_free(token, out_line, f_out);
				if (lhs.first != "DUM")
					obs.insert(lhs.first, lhs.second);
				all_markers_so_far = false;
			}
			else if (token[0] == '(')
			{
				lhs = execute_semi(token, out_line, f_out);
				if (lhs.first != "DUM")
					obs.insert(lhs.first, lhs.second);
				all_markers_so_far = false;
			}
			else if (token[0] == marker)
			{
				if (token.size() == 1)
				{
					throw_ins_error("markers with spaces not supported...", ins_line_num);
				}
				//if this is the first instruction, its a primary search
				//if (token == tokens[0])
				if (itoken == 0)
				{
					execute_primary(token, out_line, f_out);
				}
				else
				{
					bool rewind = execute_secondary(token, out_line, f_out, all_markers_so_far);
					if (rewind)
					{
						itoken = -1; //-1 so that when the for loop increments we are back to zero
						continue;
					}

				}
			}
			else
			{
				throw_ins_error("unrecognized instruction '" + token + "'", ins_line_num);
			}
			//itoken++;
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
		ss << "InstructionFile error in file " << ins_filename;
	if (ins_lnum != 0)
		ss << " on line: " << ins_lnum;
	if (out_lnum != 0)
		ss << " on output file line: " << out_lnum;
	ss << " : " << message;
	if (warn)
		cout << endl << ss.str() << endl;
	else
		throw runtime_error(ss.str());
}

string InstructionFile::parse_obs_name_from_token(const string& token)
{
	int spos, epos;
	string name;
	//whitespace obs
	if (token.at(0) == '!')
	{
		return token.substr(1, token.size() - 2);
	}
	

	pair<string, pair<int, int>> info;
	for (int i=0;i<obs_tags.size();i++)
	{
		if (token[0] == obs_tags[i].first)
		{
			info = parse_obs_instruction(token, string(1,obs_tags[i].second));
			return info.first;
		}
	}
	throw_ins_error("instruction type not recognized for observation instruction '" + token + "'");
}

vector<string> InstructionFile::tokenize_ins_line(const string& ins_line)
{
	int s, e;
	vector<string> tokens, temp_tokens, marker_tags;
	vector<int> marker_indices;
	//check for markers - might need to tokenized differently..
	if (ins_line.find(marker) != string::npos)
	{
		//get the indices of all markers on the line
		marker_indices = TemplateFile::find_all_marker_indices(ins_line, string(1, marker));
		if (marker_indices.size() % 2 != 0)
			throw_ins_error("unbalanced marker '" + string(1, marker) + "'", ins_line_num);
		
		//extract the un-altered marker tags (the strings between the markers)
		int s, e;
		for (int i = 0; i < marker_indices.size(); i = i + 2)
		{
			s = marker_indices[i];
			e = marker_indices[i + 1];
			//include the marker in the string b/c that is used later to decide
			//how to handle the token
			marker_tags.push_back(last_ins_line.substr(s, (e - s)+1));
		}
		
		//now tokenize in pieces
		// anything before the first marker
		if (marker_indices[0] != 0)
		{
			temp_tokens.clear();
			pest_utils::tokenize(pest_utils::upper_cp(ins_line.substr(0, marker_indices[0])), temp_tokens);
			tokens.insert(tokens.end(), temp_tokens.begin(), temp_tokens.end());
		}
	
		int im = 0;
		for (int i = 1; i < marker_indices.size()-1; i = i + 2)
		{
			tokens.push_back(marker_tags[im]);
			im++;
			e = marker_indices[i+1];
			s = marker_indices[i];
			temp_tokens.clear();
			pest_utils::tokenize(pest_utils::upper_cp(ins_line.substr(s + 1, (e - s) - 1)), temp_tokens);
			tokens.insert(tokens.end(),temp_tokens.begin(), temp_tokens.end());
		}
		tokens.push_back(marker_tags[marker_tags.size() - 1]);
		temp_tokens.clear();
		s = marker_indices[marker_indices.size() - 1];
		e = ins_line.size();
		pest_utils::tokenize(pest_utils::upper_cp(ins_line.substr(s+1, (e - s) - 1)), temp_tokens);
		tokens.insert(tokens.end(), temp_tokens.begin(), temp_tokens.end());
	}
	else
		pest_utils::tokenize(pest_utils::upper_cp(ins_line), tokens);
	return tokens;
}

pair<string, pair<int, int>> InstructionFile::parse_obs_instruction(const string& token, const string& close_tag)
{
	string name, temp;
	int s, e, pos = token.find(close_tag);
	if (pos == string::npos)
	{
		throw_ins_error("unbalanced (semi-)fixed observation instruction for token '" + token + "'", ins_line_num);
	}
	name = token.substr(1, pos-1);
	temp = token.substr(pos+1);
	pos = temp.find(":");
	if (pos == string::npos)
		throw_ins_error("couldn't find ':' in (semi-)fixed observation token '" + token + "'", ins_line_num);
	try
	{
		//pest_utils::convert_ip(temp.substr(0, pos), s);
		s = stoi(temp.substr(0, pos));
	}
	catch (...)
	{
		throw_ins_error("error casting first index '" + temp.substr(0, pos) + "' from (semi-)fixed observation instruction '" + token + "'");
	}
	try
	{
		//pest_utils::convert_ip(temp.substr(pos+1), e);
		e = stoi(temp.substr(pos + 1));
	}
	catch (...)
	{
		throw_ins_error("error casting second index '" + temp.substr(pos) + "' from observation instruction '" + token + "'");
	}
	pair<int, int> se(s-1, e-1);
	return pair<string, pair<int, int>>(name,se);
}

pair<string, double> InstructionFile::execute_fixed(const string& token, string& line, ifstream& f_out)
{
	string temp;
	double value;
	pair<string, pair<int, int>> info = parse_obs_instruction(token, "]");
	//use the raw last_out_line since "line" has been getting progressively truncated
	if (last_out_line.size() < info.second.second)
	{
		//throw_ins_error("output line not long enough for fixed obs instruction '" + token + "',");
		info.second.second = last_out_line.size();
	}
	int len = (info.second.second - info.second.first) + 1;
	temp = last_out_line.substr(info.second.first, len);
	try
	{
		//pest_utils::convert_ip(temp, value);
		value = stod(temp);
	}
	catch (...)
	{
		throw_ins_error("error casting fixed observation '" + token + "' from output string '" + temp + "'");
	}
	int pos = line.find(temp);
	if (pos == string::npos)
		throw_ins_error("internal error: string t: '"+temp+"' not found in line: '"+line+"'",ins_line_num,out_line_num);
	line = line.substr(pos + temp.size());
	return pair<string, double>(info.first,value);
}

pair<string, double> InstructionFile::execute_semi(const string& token, string& line, ifstream& f_out)
{
	string temp;
	vector<string> tokens;
	double value;
	pair<string, pair<int, int>> info = parse_obs_instruction(token, ")");
	//use the raw last_out_line since "line" has been getting progressively truncated
	if (last_out_line.size() < info.second.second)
	{
		//throw_ins_error("output line not long enough for semi-fixed obs instruction '" + token + "',");
		info.second.second = last_out_line.size();
	}

	int len = (info.second.second - info.second.first) + 1;
	int pos = last_out_line.find_first_not_of(", \t\n\r"+additional_delimiters, info.second.first); //include the comma here for csv files
	if (pos == string::npos)
		throw_ins_error("EOL encountered when looking for non-whitespace char in semi-fixed instruction '" + token + "'",ins_line_num,out_line_num);
	if (pos > info.second.second)
		throw_ins_error("no non-whitespace char found before end index in semi-fixed instruction '" + token + "'", ins_line_num,out_line_num);
	pest_utils::tokenize(last_out_line.substr(pos), tokens);
	temp = tokens[0];
	try
	{
		//pest_utils::convert_ip(temp, value);
		value = stod(temp);
	}
	catch (...)
	{
		throw_ins_error("error casting string '" + temp + "' to double for semi-fixed instruction", ins_line_num, out_line_num);
	}
	pos = line.find(temp);
	if (pos == string::npos)
		throw_ins_error("internal error: temp '" + temp + "' not found in line: '" + line + "'", ins_line_num, out_line_num);
	line = line.substr(pos + temp.size());
	return pair<string, double>(info.first,value);
}

pair<string, double> InstructionFile::execute_free(const string& token, string& line, ifstream& f_out)
{
	vector<string> tokens;
	pest_utils::tokenize(line, tokens,", \t\n\r" + additional_delimiters) ; //include the comma in the delimiters here
	if (tokens.size() == 0)
		throw_ins_error("error tokenizing output line ('"+last_out_line+"') for instruction '"+token+"' on line: " +last_ins_line, ins_line_num, out_line_num);
	double value;
	try
	{
		//pest_utils::convert_ip(tokens[0], value);
		value = stod(tokens[0]);
	}
	catch (...)
	{
		throw_ins_error("error converting '" + tokens[0] + "' to double on output line '" + last_out_line + "' for instruciton '"+token+"'", ins_line_num, out_line_num);
	}
	string name = token.substr(1, token.size() - 2);
	int pos = line.find(tokens[0]);
	if (pos == string::npos)
	{
		throw_ins_error("internal error: could not find free obs token '"+tokens[0]+"'", ins_line_num, out_line_num);
	}
	line = line.substr(pos + tokens[0].size());

	return pair<string, double>(name,value);
}

void InstructionFile::execute_primary(const string& token, string& line, ifstream& f_out)
{
	//check that a closing marker is found
	//this shouldnt be a prob,but good to check
	if (token.substr(token.size()-1,1) != string(1,marker))
		throw_ins_error("primary marker token '" + token + "' doesn't have a closing marker char", ins_line_num);
	int pos;
	string primary_tag = token.substr(1, token.size() - 2);
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


bool InstructionFile::execute_secondary(const string& token, string& line, ifstream& f_out, bool all_markers_so_far)
{
	//check that a closing marker is found
	int pos;
	if (token.substr(token.size()-1,1) != string(1,marker))
		throw_ins_error("secondary marker token '" + token + "' doesnt have a closing marker char");
	string secondary_tag = token.substr(1, token.size()-2);
	pos = line.find(secondary_tag);
	if (pos == string::npos)
	{
		if (all_markers_so_far)
			return true;
		else
			throw_ins_error("EOL encountered while executing secondary marker ('" + secondary_tag + "') search on output line", ins_line_num,out_line_num);
	}
	line = line.substr(pos + secondary_tag.size());
	return false;
}


void InstructionFile::execute_whitespace(const string& token, string& line, ifstream& f_out)
{
	string delims = " \t" + additional_delimiters;

	int pos = line.find_first_not_of(delims);
	if (pos == string::npos)
	{
		throw_ins_error("EOL encountered while executing whitespace instruction on output line", ins_line_num, out_line_num);
	}
	//if the cursor is already on a non-delim char, we need to read past that and then apply
	//the search
	if (pos == 0)
	{
		vector<string> tokens;
		pest_utils::tokenize(line, tokens, delims);
		pos = line.find(tokens[0]);
		if (pos == string::npos)
		{
			throw_ins_error("internal error in execute_whitespace: couldnt find first token");
		}
		line = line.substr(pos + tokens[0].size());
		pos = line.find_first_not_of(delims);
		if (pos == string::npos)
		{
			throw_ins_error("EOL encountered while executing whitespace instruction on output line", ins_line_num, out_line_num);
		}
	}
	//place the "cursor" on the first char not in delims
	line = line.substr(pos);
}


void InstructionFile::execute_line_advance(const string& token, string& line, ifstream& f_out)
{
	stringstream ss;
	int num;
	//pest_utils::convert_ip(token.substr(1), num);
	num = stoi(token.substr(1));
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

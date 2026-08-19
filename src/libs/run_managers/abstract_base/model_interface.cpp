/**
 * @file model_interface.cpp
 * @brief Implementation of model_interface.
 */
#include "network_wrapper.h"
#include "network_package.h"
#include "utilities.h"
#include "system_variables.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstring>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <thread>
#include <unordered_set>
#include "model_interface.h"
#include <limits>
#include <filesystem>

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


/**
 * @brief Model interface.
 *
 * @param _tplfile_vec Description.
 * @param _inpfile_vec Description.
 * @param _insfile_vec Description.
 * @param _outfile_vec Description.
 * @param _comline_vec Description.
 *
 * @return Description.
 */
ModelInterface::ModelInterface(vector<string> _tplfile_vec, vector<string> _inpfile_vec, vector<string>
                               _insfile_vec, vector<string> _outfile_vec, vector<string> _comline_vec):
                               insfile_vec(_insfile_vec), outfile_vec(_outfile_vec), tplfile_vec(_tplfile_vec),
                               inpfile_vec(_inpfile_vec), comline_vec(_comline_vec), fill_tpl_zeros(false),
                               additional_ins_delimiters(""),num_threads(1),tpl_force_decimal(false)
{
    //scrub any os seps from the file names

    for (auto& fname : insfile_vec)
    {
        scrub_filename_ip(fname);
    }
    for (auto& fname : outfile_vec)
    {
        scrub_filename_ip(fname);
    }
    for (auto& fname : tplfile_vec)
    {
        scrub_filename_ip(fname);
    }
    for (auto& fname : inpfile_vec)
    {
        scrub_filename_ip(fname);
    }
    for (auto& fname : comline_vec)
    {
        scrub_filename_ip(fname);
    }

}

/**
 * @brief Scrub filename ip.
 *
 * @param fname Description.
 */
void ModelInterface::scrub_filename_ip(string& fname)
{
    vector<string> tokens;
    string tname;
    pest_utils::tokenize(fname,tokens,"/\\");
    if (tokens.size()> 1)
    {
        tname = tokens[0];

        for (int i=1;i<tokens.size();i++)
        {

            tname = tname + OS_SEP + tokens[i];
        }
        //cout << fname << "," << tname << endl;
        fname = tname;
    }
}


/**
 * @brief Throw mio error.
 *
 * @param base_message Description.
 */
void ModelInterface::throw_mio_error(string base_message)
{
	throw runtime_error("model input/output error:" + base_message);
}



/**
 * @brief Check io access.
 */
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

/**
 * @brief Check tplins.
 *
 * @param par_names Description.
 * @param obs_names Description.
 */
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



/**
 * @brief Run.
 *
 * @param pars Description.
 * @param obs Description.
 */
void ModelInterface::run(Parameters* pars, Observations* obs)
{
    exception_ptr run_exception = nullptr;
	pest_utils::thread_flag terminate(false);
	pest_utils::thread_flag finished(false);

	run(&terminate, &finished, run_exception, pars, obs);
	if (run_exception)
	{
        rethrow_exception(run_exception);
	}

}

/**
 * @brief Work.
 *
 * @param tid Description.
 * @param ins_idx Description.
 * @param obs Description.
 * @param additional_ins_delims Description.
 */
void ThreadedInstructionProcess::work(int tid, vector<int>& ins_idx, Observations& obs, string additional_ins_delims)
{
	int count = 0, i;
	unique_lock<mutex> obs_guard(obs_lock, defer_lock);
	unique_lock<mutex> idx_guard(idx_lock, defer_lock);
	while (true)
	{
		i = -999;
		while (true)
		{
			if (idx_guard.try_lock())
			{
				if (ins_idx.size() == 0)
				{
					cout << "thread " << tid << " processed " << count << " instruction files" << endl;
					return;
				}
				i = ins_idx[ins_idx.size() - 1];
				ins_idx.pop_back();
				idx_guard.unlock();
				break;
			}
		}
		InstructionFile ins(insfile_vec[i]);
		ins.set_additional_delimiters(additional_ins_delims);
		Observations oobs = ins.read_output_file(outfile_vec[i]);
		while (true)
		{
			if (obs_guard.try_lock())
			{
				//pro_par_vec.push_back(pro_pars);
				obs.update_without_clear(oobs.get_keys(), oobs.get_data_vec(oobs.get_keys()));
				obs_guard.unlock();
				break;
			}
		}
		count++;
	}
}

/**
 * @brief Work.
 *
 * @param tid Description.
 * @param tpl_idx Description.
 * @param pars Description.
 * @param pro_pars Description.
 */
void ThreadedTemplateProcess::work(int tid, vector<int>& tpl_idx, Parameters pars, Parameters& pro_pars)
{
	int count = 0, i;
	unique_lock<mutex> par_guard(par_lock, defer_lock);
	unique_lock<mutex> idx_guard(idx_lock, defer_lock);
	while (true)
	{
		i = -999;
		while (true)
		{
			if (idx_guard.try_lock())
			{
				if (tpl_idx.size() == 0)
				{
					//cout << "thread " << tid << " processed " << count << " template files" << endl;
					return;
				}
				i = tpl_idx[tpl_idx.size() - 1];
				tpl_idx.pop_back();
				idx_guard.unlock();
				break;
			}
		}
		TemplateFile tpl(tplfile_vec[i]);
		tpl.set_fill_zeros(fill);
		tpl.set_force_decimal(force_decimal);
		Parameters ppars = tpl.write_input_file(inpfile_vec[i], pars);
		while (true)
		{
			if (par_guard.try_lock())
			{
				//pro_par_vec.push_back(pro_pars);
				pro_pars.update_without_clear(ppars.get_keys(), ppars.get_data_vec(ppars.get_keys()));
				par_guard.unlock();
				break;
			}
		}
		count++;
	}
}

/**
 * @brief Process template file thread.
 *
 * @param tid Description.
 * @param tpl_idx Description.
 * @param ttp Description.
 * @param pars Description.
 * @param pro_pars Description.
 * @param eptr Description.
 */
void process_template_file_thread(int tid, vector<int>& tpl_idx, ThreadedTemplateProcess& ttp, Parameters pars, Parameters& pro_pars, exception_ptr& eptr)
{
	
	try
	{
		ttp.work(tid, tpl_idx, pars, pro_pars);
	}
	catch (const std::exception& e)
	{
		eptr = current_exception();
	}
	catch (...)
	{
		eptr = current_exception();
	}

	return;
}

/**
 * @brief Process instruction file thread.
 *
 * @param tid Description.
 * @param ins_idx Description.
 * @param tip Description.
 * @param obs Description.
 * @param additional_ins_delims Description.
 * @param eptr Description.
 */
void process_instruction_file_thread(int tid, vector<int>& ins_idx, ThreadedInstructionProcess& tip, Observations& obs, string additional_ins_delims, exception_ptr& eptr)
{

	try
	{
		tip.work(tid, ins_idx, obs, additional_ins_delims);
	}
	catch (const std::exception& e)
	{
		eptr = current_exception();
	}
	catch (...)
	{
		eptr = current_exception();
	}

	return;
}


/**
 * @brief Write input files.
 *
 * @param pars_ptr Description.
 */
void ModelInterface::write_input_files(Parameters *pars_ptr)
{
	int nnum_threads = num_threads;
	if (nnum_threads > tplfile_vec.size())
		nnum_threads = tplfile_vec.size();
	std::chrono::system_clock::time_point start_time = chrono::system_clock::now();
    if (should_echo)
        cout << pest_utils::get_time_string() << " processing template files with " << nnum_threads << " threads..." << endl;
	vector<thread> threads;
	vector<exception_ptr> exception_ptrs;
	Parameters pro_pars = *pars_ptr; //copy
	ThreadedTemplateProcess ttp(tplfile_vec, inpfile_vec, fill_tpl_zeros, tpl_force_decimal);

	for (int i = 0; i < nnum_threads; i++)
	{
		exception_ptrs.push_back(exception_ptr());
	}

	vector<int> tpl_idx;
	for (int i = 0; i < tplfile_vec.size(); i++)
		tpl_idx.push_back(i);

	for (int i = 0; i < nnum_threads; i++)
	{
		threads.push_back(thread(process_template_file_thread, i, std::ref(tpl_idx), std::ref(ttp), *pars_ptr, std::ref(pro_pars), std::ref(exception_ptrs[i])));
	}
	stringstream ss;
	int num_exp = 0;
	for (int i = 0; i < nnum_threads; ++i)
	{
		bool found = false;
		if (exception_ptrs[i])
		{
			found = true;
			try
			{
				rethrow_exception(exception_ptrs[i]);
			}
			catch (const std::exception& e)
			{
				//stringstream ss;
				ss << " thread processing template file raised an exception: " << e.what() << endl;
				num_exp++;
				//cout << "Error: " << ss.str();
				//throw runtime_error(ss.str());
			}
			catch (...)
			{
				//stringstream ss;
				ss << " thread processing template file raised an exception" << endl;
				num_exp++;
				//cout << "Error: " << ss.str();
				//throw runtime_error(ss.str());
			}
		}
		threads[i].join();
		if ((exception_ptrs[i]) && (!found))
		{
			try
			{
				rethrow_exception(exception_ptrs[i]);
			}
			catch (const std::exception& e)
			{
				//stringstream ss;
				ss << " thread processing template file raised an exception: " << e.what() << endl;
				num_exp++;
				//cout << "Error: " << ss.str();
				//throw runtime_error(ss.str());
			}
			catch (...)
			{
				//stringstream ss;
				ss << " thread processing template file raised an exception" << endl;
				num_exp++;
				//cout << "Error: " << ss.str();
				//throw runtime_error(ss.str());
			}
		}
	}

	//if (templatefiles.size() == 0)
	//{
	//	for (auto t : tplfile_vec)
	//	{
	//		TemplateFile tt(t);
	//		tt.set_fill_zeros(fill_tpl_zeros);
	//		templatefiles.push_back(tt);
	//	}
	//}

	//
	//vector<string> notnormal = pars_ptr->get_notnormal_keys();
	//if (notnormal.size() > 0)
	//{
	//	throw runtime_error("denormal floating point parameter values found");
	//}

	//vector<Parameters> pro_par_vec;
	//for (int i = 0; i < templatefiles.size(); i++)
	//{
	//	string name = templatefiles[i].get_tpl_filename();
	//	//cout << name << endl;
	//	pro_par_vec.push_back(templatefiles[i].write_input_file(inpfile_vec[i], *pars_ptr));
	//}
	//update pars to account for possibly truncated par values...important for jco calcs
	//for (auto pro_pars : pro_par_vec)
	//	pars_ptr->update_without_clear(pro_pars.get_keys(), pro_pars.get_data_vec(pro_pars.get_keys()));
	if (num_exp > 0)
	{
		//cout << "errors processing template files: " << endl << ss.str();
		throw runtime_error(ss.str());
	}

	pars_ptr->update_without_clear(pro_pars.get_keys(), pro_pars.get_data_vec(pro_pars.get_keys()));
    if (should_echo)
        cout << pest_utils::get_time_string() << " done, took " << pest_utils::get_duration_sec(start_time) << " seconds" << endl;
}

/**
 * @brief Read output files.
 *
 * @param obs Description.
 */
void ModelInterface::try_read_output_files(Observations *obs, vector<string>& missing,
	vector<string>& problems)
{
	// built from insfile_vec, the way ThreadedInstructionProcess::work does. The
	// `instructionfiles` member looks like the right thing to iterate and is always empty -
	// its only push_back sits inside a commented-out block.
	for (int i = 0; i < insfile_vec.size(); i++)
	{
		Observations file_obs;
		vector<string> file_missing, file_problems;
		InstructionFile ins(insfile_vec[i]);
		ins.set_additional_delimiters(additional_ins_delimiters);
		// never throws, by contract - so one unreadable output file cannot cost us the
		// results sitting in the others, which is the entire point of asking
		ins.try_read_output_file(outfile_vec[i], file_obs, file_missing, file_problems);
		vector<string> keys = file_obs.get_keys();
		obs->update_without_clear(keys, file_obs.get_data_vec(keys));
		missing.insert(missing.end(), file_missing.begin(), file_missing.end());
		problems.insert(problems.end(), file_problems.begin(), file_problems.end());
	}
	// fill the gaps with the sentinel so the caller gets a complete vector, and keep the
	// names alongside so it never has to detect absence by comparing against a magic number
	if (missing.size() > 0)
	{
		vector<double> fill(missing.size(), Transformable::no_data);
		obs->update_without_clear(missing, fill);
	}
	sort(missing.begin(), missing.end());
}

void ModelInterface::read_output_files(Observations *obs)
{
	int nnum_threads = num_threads;
	if (nnum_threads > insfile_vec.size())
		nnum_threads = insfile_vec.size();
	std::chrono::system_clock::time_point start_time = chrono::system_clock::now();
    if (should_echo)
        cout << pest_utils::get_time_string() <<  " processing instruction files with " << nnum_threads << " threads..." << endl;
	vector<thread> threads;
	vector<exception_ptr> exception_ptrs;
	Observations temp_obs;

	ThreadedInstructionProcess tip(insfile_vec, outfile_vec);

	for (int i = 0; i < nnum_threads; i++)
	{
		exception_ptrs.push_back(exception_ptr());
	}

	vector<int> ins_idx;
	for (int i = 0; i < insfile_vec.size(); i++)
		ins_idx.push_back(i);

	for (int i = 0; i < nnum_threads; i++)
	{
		threads.push_back(thread(process_instruction_file_thread, i, std::ref(ins_idx), std::ref(tip), std::ref(temp_obs), additional_ins_delimiters,
			std::ref(exception_ptrs[i])));
	}
	stringstream ss;
	int num_exp = 0;
	for (int i = 0; i < nnum_threads; ++i)
	{
		bool found = false;
		if (exception_ptrs[i])
		{
			found = true;
			try
			{
				rethrow_exception(exception_ptrs[i]);
			}
			catch (const std::exception& e)
			{
				//stringstream ss;
				ss << " thread processing instruction file raised an exception: " << e.what() << endl;
				//cout << "Error: " << ss.str() << endl;
				//throw runtime_error(ss.str());
				num_exp++;
			}
			catch (...)
			{
				//stringstream ss;
				ss << " thread processing instruction file raised an exception" << endl;
				//cout << "Error: " << ss.str() << endl;
				//throw runtime_error(ss.str());
				num_exp++;
			}
		}
		threads[i].join();
		if ((exception_ptrs[i]) && (!found))
		{
			try
			{
				rethrow_exception(exception_ptrs[i]);
			}
			catch (const std::exception& e)
			{
				//stringstream ss;
				ss << " thread processing instruction file raised an exception: " << e.what() << endl;
				//cout << "Error: " << ss.str() << endl;
				//throw runtime_error(ss.str());
				num_exp++;
			}
			catch (...)
			{
				//stringstream ss;
				ss << " thread processing instruction file raised an exception" << endl;
				//cout << "Error: " << ss.str() << endl;
				//throw runtime_error(ss.str());
				num_exp++;
			}
		}
	}
	
	
	
	/*if (instructionfiles.size() == 0)
	{
		for (auto i : insfile_vec)
		{
			InstructionFile ii(i);
			ii.set_additional_delimiters(additional_ins_delimiters);
			instructionfiles.push_back(ii);
		}
	}
	
	Observations temp_obs, pro_obs;
	for (int i = 0; i < instructionfiles.size(); i++)
	{
		pro_obs = instructionfiles[i].read_output_file(outfile_vec[i]);
		temp_obs.update_without_clear(pro_obs.get_keys(), pro_obs.get_data_vec(pro_obs.get_keys()));
	}*/

	if (num_exp > 0)
	{
		//cout << "errors processing instruction files: " << endl << ss.str();
		throw runtime_error(ss.str());
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
    if (should_echo)
	    cout << pest_utils::get_time_string() << " done, took " << pest_utils::get_duration_sec(start_time) << " seconds" << endl;

}


/**
 * @brief Remove existing.
 */
void ModelInterface::remove_existing()
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
		// One call, not check-then-act. std::filesystem::remove() returns true when it
		// deleted something and false when there was nothing there, and only sets the
		// error_code on a real failure - so "already gone" is not an error and there is no
		// window between the test and the removal for anything to change.
		//
		// It also stops the test from being destructive. The check_exist_out() this replaced
		// opens an ofstream, which CREATES the file when it is missing and TRUNCATES it when
		// it is not, so the old pair created every absent model file purely in order to
		// delete it again.
		std::error_code ec;
		for (auto& out_file : outfile_vec)
		{
			std::filesystem::remove(out_file, ec);
			if (ec)
			{
				failed_file_vec.push_back(out_file);
				failed_file_op = true;
				ec.clear();
			}
		}
		for (auto& in_file : inpfile_vec)
		{
			std::filesystem::remove(in_file, ec);
			if (ec)
			{
				failed_file_vec.push_back(in_file);
				failed_file_op = true;
				ec.clear();
			}
		}
		if (failed_file_op)
		{
			++n_tries;
			w_sleep(500);
			if (n_tries > 5)
			{
				ostringstream str;
				str << "model interface error: Cannot delete existing following model files:";
				for (const string& ifile : failed_file_vec)
				{
					str << " " << ifile;
				}
				throw PestError(str.str());
			}
		}
	}
}


void ModelInterface::run(pest_utils::thread_flag* terminate, pest_utils::thread_flag* finished, exception_ptr& eptr,
						Parameters* pars_ptr, Observations* obs_ptr)
{
			
	try
	{
		remove_existing();
		// from here the output files on disk can only be THIS run's: the previous run's are
		// gone, and anything that appears was written by the model we are about to call. Until
		// this point a partial read would have returned the previous realization's results.
		if (outputs_cleared) outputs_cleared->store(true);
		write_input_files(pars_ptr);
		
		std::chrono::system_clock::time_point start_time = chrono::system_clock::now();
		if (should_echo)
            cout << pest_utils::get_time_string() << " calling forward run command(s)" << endl;

		// The OS-specific spawn/poll/kill is pest_utils::run_commands() now - windows job
		// objects, posix process groups, and the grandchild-killing each exists for. It lived
		// here and in a drifted near-copy in RunManagerExternal; one of them is enough.
		//
		// It does NOT touch `finished`: that flag is this class's contract with its listener
		// thread, not something a general command runner should know about. The catch blocks
		// below set it instead, which also covers write_input_files()/read_output_files()
		// throwing - paths that previously left the listener waiting forever.
		bool term_break = !pest_utils::run_commands(comline_vec, terminate, should_echo, sleep_ms);
        if (should_echo)
		    cout << pest_utils::get_time_string() << " forward run command(s) finished, took " << pest_utils::get_duration_sec(start_time) << " seconds" << endl;

		if (term_break) return;

		read_output_files(obs_ptr);

		//set the finished flag for the listener thread
		finished->set(true);

	}
	catch (const std::exception& e)
	{
		cout << "exception raised by run thread: " << e.what() << endl;
		// the listener thread waits on this; without it a failed run hangs rather than reports
		finished->set(true);
		eptr = current_exception();
	}
	catch (...)
	{
		cout << "exception raised by run thread" << endl;
		finished->set(true);
		eptr = current_exception();
	}
	return;

}

/**
 * @brief Parse and check.
 *
 * @return Description.
 */
unordered_set<string> TemplateFile::parse_and_check()
{
	ifstream f(tpl_filename);
	prep_tpl_file_for_reading(f);
	return get_names(f);

}

/**
 * @brief Write input file.
 *
 * @param input_filename Description.
 * @param pars Description.
 *
 * @return Description.
 */
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
	//vector<string> t = pars.get_keys();
	//unordered_set<string> pnames(t.begin(), t.end());
	//unordered_set<string>::iterator end = pnames.end();
	//t.resize(0);
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
			//if (pnames.find(name) == end)
			//	throw_tpl_error("parameter '" + name + "' not listed in control file");

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
			try
			{
				val = pars.get_rec(t.first);
			}
			catch (...)
			{
				throw_tpl_error("parameter '" + name + "' not in parameters instance");
			}
			val_str = cast_to_fixed_len_string(t.second.second, val, name);
			line.replace(t.second.first, t.second.second, val_str);
			//pest_utils::convert_ip(val_str, val);
			val = stod(val_str);
			pro_pars.insert(name, val);
		}
		f_in << line << endl;
		if (f_in.bad())
		{
			throw_tpl_error("ofstream is bad after writing line '" + line + "'", line_num);
		}
	}
	f_tpl.close();
	f_in.close();
	if (f_in.bad())
	{
		throw_tpl_error("ofstream is bad after closing file, something is probably corrupt");
	}
	return pro_pars;
}

/**
 * @brief Prep tpl file for reading.
 *
 * @param f_tpl Description.
 */
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

/**
 * @brief Get names.
 *
 * @param f Description.
 *
 * @return Description.
 */
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

/**
 * @brief Find all marker indices.
 *
 * @param line Description.
 * @param marker Description.
 *
 * @return Description.
 */
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

/**
 * @brief Throw tpl error.
 *
 * @param message Description.
 * @param lnum Description.
 * @param warn Description.
 */
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

/**
 * @brief Parse tpl line.
 *
 * @param line Description.
 *
 * @return Description.
 */
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

/**
 * @brief Cast to fixed len string.
 *
 * @param size Description.
 * @param value Description.
 * @param name Description.
 *
 * @return Description.
 */
string TemplateFile::cast_to_fixed_len_string(int size, double value, string& name)
{
	string val_str, fill_val=" ";
	stringstream ss;
	int precision = size;
	bool sci = false;
	if (value < 0)
		precision--; // for the minus sign
	if ((!force_decimal) && ((abs(value) >= 100) || (abs(value) < 0.01)))
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
			//time for desperate measures:
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
			if ((e_idx != -1) && (r_idx != -1) && (r_idx != e_idx))
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
			//ok now just drop everything right of the radix and the radix
			else if ((e_idx == -1) && (r_idx != -1))
            {
			    string t = val_str.substr(0,r_idx-1);
			    if (t.size() <= size)
			        val_str = t;
			        break;
            }
			ss.str("");
		 	ss << "TemplateFile casting error: can't represent value " << value;
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
		//if the fill value isn't a space and its a negative value
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

/**
 * @brief Read line.
 *
 * @param f_tpl Description.
 *
 * @return Description.
 */
string TemplateFile::read_line( ifstream& f_tpl)
{
	if (f_tpl.bad())
		throw_tpl_error("can't read next line", line_num);
	string line;
	if (f_tpl.eof())
		throw_tpl_error("unexpected eof", line_num);
	
	getline(f_tpl, line);
	pest_utils::strip_ip(line, "\n\r");
	line_num++;
	return line;
}

/**
 * @brief Read ins line.
 *
 * @param f_ins Description.
 *
 * @return Description.
 */
/**
 * @brief PROTOTYPE: turn an output-file token into a double, or say why not.
 *
 * See the header for what this is for. The order matters: a straight parse first, so anything
 * that reads correctly today reads identically, and only then the fortran repair.
 *
 * The repair is narrow ON PURPOSE. The whole risk of a fall-back like this is that it turns
 * something that is not a number into a number - "10-20" into 1.0e-19, a date into a
 * magnitude - and hands back a confident wrong answer, which is worse than the refusal it
 * replaced. So it insists on all of:
 *
 *   - a mantissa with a decimal point in it, which is what kills "10-20" and "2024-01-15";
 *     every fortran E/ES/D descriptor writes the point
 *   - exactly one exponent sign, after at least one mantissa digit, not at the very end
 *   - nothing but digits after that sign, at most three of them (a double cannot need more:
 *     the extreme exponents are 308 and 324)
 *   - nothing left over
 *
 * A 'd' or 'D' exponent marker is handled by the same path, since it fails for the same
 * reason and is the same fortran-ism.
 */
bool InstructionFile::try_parse_double(const string& token, double& value, string& why)
{
	value = 0.0;
	why.clear();
	if (token.empty())
	{
		why = "empty token";
		return false;
	}

	//classify what strtod made of a candidate string.  errno rather than the return value
	//alone, because a subnormal and a genuine small number are not distinguishable after the
	//fact, and inf could have been written literally
	auto convert = [&](const string& text, double& out, string& problem) -> bool
	{
		errno = 0;
		const char* start = text.c_str();
		char* end = nullptr;
		double v = strtod(start, &end);
		if (end == start)
		{
			problem = "no number found";
			return false;
		}
		//trailing whitespace is fine; anything else is left-over junk
		while ((*end != '\0') && (isspace((unsigned char)*end) != 0))
			end++;
		if (*end != '\0')
		{
			problem = "left-over chars: '" + string(end) + "'";
			return false;
		}
		if (isnan(v))
		{
			problem = "not a number";
			return false;
		}
		if (isinf(v))
		{
			//either "inf" was written literally or the exponent is past DBL_MAX.  both mean
			//the model has no answer for us, and both must stay refusals
			problem = "value overflows a double";
			return false;
		}
		if ((errno == ERANGE) || ((v != 0.0) && (!isnormal(v))))
		{
			//underflow: subnormal, or all the way to zero.  ERANGE alone is not enough to
			//call it underflow - overflow sets it too - but inf was already handled above
			out = 0.0;
			return true;
		}
		out = v;
		return true;
	};

	string problem;
	if (convert(token, value, problem))
		return true;

	//straight parse failed - try the fortran forms
	string repaired = token;
	bool has_marker = false;
	for (size_t i = 0; i < repaired.size(); i++)
	{
		if ((repaired[i] == 'd') || (repaired[i] == 'D'))
		{
			repaired[i] = 'E';
			has_marker = true;
			break;
		}
		if ((repaired[i] == 'e') || (repaired[i] == 'E'))
		{
			has_marker = true;
			break;
		}
	}
	if (!has_marker)
	{
		//no exponent marker at all: look for the exponent's sign, which must sit after the
		//mantissa rather than in front of it
		size_t sign_pos = string::npos;
		int n_signs = 0;
		for (size_t i = 1; i < repaired.size(); i++)
		{
			if ((repaired[i] == '+') || (repaired[i] == '-'))
			{
				n_signs++;
				sign_pos = i;
			}
		}
		bool ok = (n_signs == 1) && (sign_pos != string::npos) && (sign_pos + 1 < repaired.size());
		if (ok)
		{
			string mantissa = repaired.substr(0, sign_pos);
			string exponent = repaired.substr(sign_pos + 1);
			//the mantissa must contain a digit AND a decimal point, and the exponent must be
			//nothing but one to three digits
			bool m_digit = false, m_point = false;
			for (auto c : mantissa)
			{
				if (isdigit((unsigned char)c) != 0) m_digit = true;
				else if (c == '.') m_point = true;
			}
			//EXACTLY three digits, not "up to three".  fortran drops the E/D only when the
			//exponent needs three digits - a one or two digit exponent always keeps its
			//marker, zero-padded (1e-99 goes out as E-99, never as -99), and that holds for
			//es/e/en/d/g editing and at any field width.  so anything with a shorter exponent
			//was not written by a fortran format, and reading "1.5-2" as 0.015 would be
			//inventing a number out of what is much more likely a subtraction or a range.
			bool e_digits = (exponent.size() == 3);
			for (auto c : exponent)
				if (isdigit((unsigned char)c) == 0)
					e_digits = false;
			if (m_digit && m_point && e_digits)
				repaired.insert(sign_pos, "E");
			else
				ok = false;
		}
		if (!ok)
		{
			why = problem;
			return false;
		}
	}
	if (repaired == token)
	{
		why = problem;
		return false;
	}
	string second_problem;
	if (convert(repaired, value, second_problem))
		return true;
	why = second_problem;
	return false;
}

string InstructionFile::read_ins_line(ifstream& f_ins)
{
	if (f_ins.bad())
		throw_ins_error("can't read next instruction file line", ins_line_num);
	string line;
	if (f_ins.eof())
		throw_ins_error("unexpected instruction file eof ", ins_line_num);

	getline(f_ins, line);
	last_ins_line = line;
	ins_line_num++;
	return line;
}


/**
 * @brief Read out line.
 *
 * @param f_out Description.
 *
 * @return Description.
 */
string InstructionFile::read_out_line(ifstream& f_out)
{
	if (f_out.bad())
		throw_ins_error("can't read next output file line", ins_line_num, out_line_num);
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


/**
 * @brief Parse and check.
 *
 * @return Description.
 */
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
		//pest_utils::tokenize(line, tokens);
        tokens = tokenize_ins_line(line);
		for (int i = 0; i < tokens.size(); i++)
		{
			first = tokens[i].at(0);
			if ((first == '!') || (first == '(') || (first == '['))
			{
				name = parse_obs_name_from_token(tokens[i]);
				if (name != "DUM")
				{
					if (names.find(name) != names.end())
					{
						cout << line << endl;
						throw_ins_error("observation '" + name + "' listed multiple times in ins file '" + ins_filename + "'");
					}
					names.emplace(name);
				}
			}
			else if ((first == marker) || (first == 'L') || (first == 'W'))
            {

            }
			else
            {
			    stringstream ss;
			    ss << "unrecognized instruction: '" << tokens[i];
                throw_ins_error(ss.str(),ins_line_num);
            }

		}
	}
	f_ins.close();
	
	return names;
}

/**
 * @brief Prep ins file for reading.
 *
 * @param f_ins Description.
 */
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
    if (marker == '!')
    {
        throw_ins_error("the bang ('!') can't be used as the marker bc it is used to indicate free-format instructions");
    }
}


/**
 * @brief Read output file.
 *
 * @param output_filename Description.
 *
 * @return Description.
 */
Observations InstructionFile::read_output_file(const string& output_filename)
{
	Observations obs;
	read_output_file(output_filename, obs);
	return obs;
}

void InstructionFile::read_output_file(const string& output_filename, Observations& obs)
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
	f_ins.close();
	f_out.close();
}

/**
 * @brief Copy the complete (newline-terminated) prefix of a file to a temp file.
 *
 * Returns the temp file's name, or "" when the file has no complete line yet. Used only by
 * the tolerant reader - see try_read_output_file() for why a partial final line must not be
 * parsed. Written to the system temp directory rather than beside the output file, because
 * the model may still be writing in that directory.
 */
string InstructionFile::write_complete_lines_to_temp(const string& output_filename)
{
	ifstream f_in(output_filename, ios::binary);
	if (!f_in.good())
		throw runtime_error("could not open '" + output_filename + "' for reading");
	stringstream buf;
	buf << f_in.rdbuf();
	f_in.close();
	string text = buf.str();
	size_t last = text.find_last_of('\n');
	if (last == string::npos)
		return string();                 // nothing complete to read yet
	// BESIDE the output file, never in the shared system temp directory.
	//
	// temp_directory_path() is machine-wide and this name derives only from the output file's
	// BASENAME - which is identical for every agent running the same case. So every agent on
	// one host wrote to, and read back, one single file. Whichever wrote last won, and all of
	// them then reported ITS observations as their own: four agents on four realizations
	// returned byte-identical values, and mid-run screening abandoned runs on the strength of
	// another realization's results.
	//
	// The output file's own directory is the right home because PANTHER already guarantees one
	// working directory per agent, so it is per-agent by construction. A relative
	// output_filename gives an empty parent_path(), which resolves to the agent's cwd - the
	// same directory - so both forms land in the right place.
	std::filesystem::path out_path(output_filename);
	std::filesystem::path tmp = out_path.parent_path() /
		(out_path.filename().string() + ".pestpp_partial");
	ofstream f_out(tmp, ios::binary | ios::trunc);
	if (!f_out.good())
		throw runtime_error("could not open a temp file for the partial read");
	f_out.write(text.data(), (streamsize)(last + 1));
	f_out.close();
	return tmp.string();
}

void InstructionFile::try_read_output_file(const string& output_filename, Observations& obs,
	vector<string>& missing, vector<string>& problems)
{
	// every name this instruction file is responsible for. Asked up front because it is the
	// only way to report what is MISSING after a parse that stopped early - and because a
	// file that could not be opened at all still has to name everything it was going to
	// supply rather than reporting nothing at all.
	unordered_set<string> covered;
	try
	{
		covered = parse_and_check();
	}
	catch (const exception& e)
	{
		// the instruction file itself is unreadable, which is not a partial-results problem:
		// there is nothing to be partial about, and nothing can be named as missing
		problems.push_back(string("could not parse instruction file '") + ins_filename +
			"': " + e.what());
		return;
	}
	catch (...)
	{
		problems.push_back(string("could not parse instruction file '") + ins_filename + "'");
		return;
	}

	// A missing or unopenable output file is NOT an error here. A run that has not written
	// its output yet is the ordinary case when asking a run in flight what it has.
	if (!pest_utils::check_exist_in(output_filename))
	{
		problems.push_back("output file '" + output_filename + "' not found (yet)");
	}
	else
	{
		// Only COMPLETE lines are trusted. A run that is still going may be mid-write, and an
		// unterminated final line is indistinguishable from a finished one to getline() - so
		// "1.2345" caught at "1.2" parses cleanly and would be reported as a confident WRONG
		// value rather than an absent one. That is the one outcome partial results must never
		// produce: a caller computing phi from them would get a plausible wrong answer.
		//
		// The cost is a legitimately-complete last line with no trailing newline being
		// skipped. That is the right trade here - a run that has genuinely FINISHED is read
		// by the strict path, which is unchanged and reads it.
		string safe_name;
		try
		{
			safe_name = write_complete_lines_to_temp(output_filename);
		}
		catch (const exception& e)
		{
			problems.push_back(string("could not stage output file '") + output_filename +
				"' for a partial read: " + e.what());
		}
		if (safe_name.empty())
		{
			if (problems.empty())
				problems.push_back("output file '" + output_filename +
					"' has no complete lines yet");
		}
		else
		{
			try
			{
				read_output_file(safe_name, obs);
			}
			catch (const exception& e)
			{
				// obs keeps whatever was read before this - that is the point of the reference
				problems.push_back(e.what());
			}
			catch (...)
			{
				problems.push_back("unknown error reading output file '" + output_filename + "'");
			}
			std::error_code ec;
			std::filesystem::remove(safe_name, ec);
		}
	}

	vector<string> have = obs.get_keys();
	unordered_set<string> hset(have.begin(), have.end());
	for (const string& name : covered)
	{
		if (name == "DUM")            // the throwaway name; never a real observation
			continue;
		if (hset.find(name) == hset.end())
			missing.push_back(name);
	}
	sort(missing.begin(), missing.end());
}


/**
 * @brief Throw ins error.
 *
 * @param message Description.
 * @param ins_lnum Description.
 * @param out_lnum Description.
 * @param warn Description.
 */
void InstructionFile::throw_ins_error(const string& message, int ins_lnum, int out_lnum, bool warn)
{
	stringstream ss;
	if (warn)
		ss << "InstructionFile warning in '" << ins_filename << "'";
	else
		ss << "InstructionFile error in file '" << ins_filename << "'";
	if (ins_lnum != 0)
		ss << " on instruction file line: " << ins_lnum;
	if (out_lnum != 0)
		ss << ", on output file line: " << out_lnum;
	ss << " : " << message;
	if (warn)
		cout << endl << ss.str() << endl;
	else
		throw runtime_error(ss.str());
}

/**
 * @brief Parse obs name from token.
 *
 * @param token Description.
 *
 * @return Description.
 */
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
	return "";
}

/**
 * @brief Tokenize ins line.
 *
 * @param ins_line Description.
 *
 * @return Description.
 */
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

/**
 * @brief Parse obs instruction.
 *
 * @param token Description.
 * @param close_tag Description.
 *
 * @return Description.
 */
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
		throw_ins_error("error casting first index '" + temp.substr(0, pos) + "' from (semi-)fixed observation instruction '" + token + "'", ins_line_num);
	}
	try
	{
		//pest_utils::convert_ip(temp.substr(pos+1), e);
		e = stoi(temp.substr(pos + 1));
	}
	catch (...)
	{
		throw_ins_error("error casting second index '" + temp.substr(pos) + "' from (semi)-fixed observation instruction '" + token + "'");
	}
	pair<int, int> se(s-1, e-1);
	return pair<string, pair<int, int>>(name,se);
}

/**
 * @brief Execute fixed.
 *
 * @param token Description.
 * @param line Description.
 * @param f_out Description.
 *
 * @return Description.
 */
pair<string, double> InstructionFile::execute_fixed(const string& token, string& line, ifstream& f_out)
{
	string temp;
	double value = 1.0e+30;
	pair<string, pair<int, int>> info = parse_obs_instruction(token, "]");
	//use the raw last_out_line since "line" has been getting progressively truncated
	if (last_out_line.size() < info.second.second)
	{
		//throw_ins_error("output line not long enough for fixed obs instruction '" + token + "',");
		info.second.second = last_out_line.size();
	}
	int len = (info.second.second - info.second.first) + 1;
	temp = last_out_line.substr(info.second.first, len);
	string why;
	if ((!try_parse_double(temp, value, why)) && (info.first != "DUM"))
	{
		throw_ins_error("error converting '" + temp + "' to double on output line '" + last_out_line + "' for fixed instruction: '" + token + "': " + why, ins_line_num, out_line_num);
	}
	int pos = line.find(temp);
	if (pos == string::npos)
		throw_ins_error("internal error: string t: '"+temp+"' not found in line: '"+line+"'",ins_line_num,out_line_num);
	line = line.substr(pos + temp.size());
	
	return pair<string, double>(info.first,value);
}

/**
 * @brief Execute semi.
 *
 * @param token Description.
 * @param line Description.
 * @param f_out Description.
 *
 * @return Description.
 */
pair<string, double> InstructionFile::execute_semi(const string& token, string& line, ifstream& f_out)
{
	string temp;
	vector<string> tokens;
	double value = 1.0e+30;
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
		throw_ins_error("EOL encountered when looking for non-whitespace char in semi-fixed instruction '" + token + "' on line: '" + line + "'",ins_line_num,out_line_num);
	if (pos > info.second.second)
		throw_ins_error("no non-whitespace char found before end index in semi-fixed instruction '" + token + "' on line: '" + line + "'", ins_line_num,out_line_num);
	pest_utils::tokenize(last_out_line.substr(pos), tokens);
	temp = tokens[0];
	string why;
	if ((!try_parse_double(temp, value, why)) && (info.first != "DUM"))
	{
		throw_ins_error("error converting '" + temp + "' to double on output line '" + last_out_line + "' for semi-fixed instruction: '" + token + "': " + why, ins_line_num, out_line_num);
	}
	pos = line.find(temp);
	if (pos == string::npos)
		throw_ins_error("internal error: temp '" + temp + "' not found in line: '" + line + "'", ins_line_num, out_line_num);
	line = line.substr(pos + temp.size());
	return pair<string, double>(info.first,value);
}



/**
 * @brief Execute free.
 *
 * @param token Description.
 * @param line Description.
 * @param f_out Description.
 *
 * @return Description.
 */
pair<string, double> InstructionFile::execute_free(const string& token, string& line, ifstream& f_out)
{
	int tsize = line.size() / 20;
	if (tsize > 50)
		tsize = 50;
	string name = token.substr(1, token.size() - 2);
	vector<string> tokens;
	tokens.reserve(tsize);
	tokenize(line, tokens,", \t\n\r" + additional_delimiters,true,1) ; //include the comma in the delimiters here
	if (tokens.size() == 0)
		throw_ins_error("error tokenizing output line ('"+last_out_line+"') for free instruction '"+token+"' on line: " +last_ins_line, ins_line_num, out_line_num);
	double value = 1.0e+30;
	string why;
	if ((!try_parse_double(tokens[0], value, why)) && (name != "DUM"))
	{
		throw_ins_error("error converting '" + tokens[0] + "' to double on output line '" + last_out_line + "' for free instruction: '" + token + "': " + why, ins_line_num, out_line_num);
	}
	int pos = line.find(tokens[0]);
	if (pos == string::npos)
	{
		throw_ins_error("internal error: could not find free obs token '"+tokens[0]+"'", ins_line_num, out_line_num);
	}
	line = line.substr(pos + tokens[0].size());

	return pair<string, double>(name,value);
}

/**
 * @brief Tokenize.
 *
 * @param str Description.
 * @param tokens Description.
 * @param delimiters Description.
 * @param trimEmpty Description.
 * @param mx_tokens Description.
 */
void InstructionFile::tokenize(const std::string& str, vector<string>& tokens, const std::string& delimiters, const bool trimEmpty, int mx_tokens)
{
	std::string::size_type pos, lastPos = 0;
	while (true)
	{
		pos = str.find_first_of(delimiters, lastPos);
		if (pos == std::string::npos)
		{
			pos = str.length();
			if (pos != lastPos || !trimEmpty)
			{
				tokens.push_back(string(str.data() + lastPos, string::size_type(pos - lastPos)));
			}
			break;
		}
		else
		{
			if (pos != lastPos || !trimEmpty)
				tokens.push_back(string(str.data() + lastPos, string::size_type(pos - lastPos)));
		}
		if ((mx_tokens > 0) && (tokens.size() > mx_tokens))
			break;

		lastPos = pos + 1;
	}
}

/**
 * @brief Execute primary.
 *
 * @param token Description.
 * @param line Description.
 * @param f_out Description.
 */
void InstructionFile::execute_primary(const string& token, string& line, ifstream& f_out)
{
	//check that a closing marker is found
	//this shouldn't be a prob,but good to check
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


/**
 * @brief Execute secondary.
 *
 * @param token Description.
 * @param line Description.
 * @param f_out Description.
 * @param all_markers_so_far Description.
 *
 * @return Description.
 */
bool InstructionFile::execute_secondary(const string& token, string& line, ifstream& f_out, bool all_markers_so_far)
{
	//check that a closing marker is found
	int pos;
	if (token.substr(token.size()-1,1) != string(1,marker))
		throw_ins_error("secondary marker token '" + token + "' doesn't have a closing marker char");
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


/**
 * @brief Execute whitespace.
 *
 * @param token Description.
 * @param line Description.
 * @param f_out Description.
 */
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
			throw_ins_error("internal error in execute_whitespace: couldn't find first token");
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


/**
 * @brief Execute line advance.
 *
 * @param token Description.
 * @param line Description.
 * @param f_out Description.
 */
void InstructionFile::execute_line_advance(const string& token, string& line, ifstream& f_out)
{
	stringstream ss;
	int num;
	//pest_utils::convert_ip(token.substr(1), num);
	num = stoi(token.substr(1));
	if (num < 1)
		throw_ins_error("line advance instruction error: number of lines must be greater or equal to 1, not '" + token.substr(1) + "'",ins_line_num,out_line_num);
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

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
/**
 * @file RunStorage.cpp
 * @brief Implementation of RunStorage.
 */



#include <sstream>
#include <cstdio>
#include <cassert>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "RunStorage.h"
#include "Serialization.h"
#include "Transformable.h"
#include <limits>
#include "utilities.h"

using std::numeric_limits;

using namespace std;

const double RunStorage::no_data = -9999.0;

/**
 * @brief everything we know about a stream operation that failed, in one message.
 *
 * see the header for why. short version - these show up in ci and dont reproduce, so the
 * message we throw right when it happens is the only evidence anybody is ever going to get.
 * "stream not good" wasnt enough to tell a short file from a bad offset from a disk error, and
 * three separate go-arounds on this died right there.
 *
 * gets the file size through its own ifstream instead of buf_stream. by the time we get called
 * buf_stream is already in a bad state, so seeking it to measure the file would either fail too
 * or, worse, clear the very state we are trying to report.
 */
string RunStorage::stream_diagnostic(const string &where, int run_id,
	std::streamoff expected_pos) const
{
	stringstream ss;
	try
	{
		const int saved_errno = errno;
		ss << "RunStorage::" << where << "() stream not good";
		ss << " [file='" << filename << "'";

		// which bit is set tells you a lot by itself. just eof means we asked for bytes past the
		// end, so the file is short. bad means the os said no.
		ss << " state=";
		if (buf_stream.good()) ss << "good";
		else
		{
			if (buf_stream.eof())  ss << "eof|";
			if (buf_stream.fail()) ss << "fail|";
			if (buf_stream.bad())  ss << "bad|";
		}
		ss << " open=" << (buf_stream.is_open() ? "yes" : "no");

		if (run_id >= 0)
			ss << " run_id=" << run_id;
		if (expected_pos >= 0)
			ss << " target_pos=" << (long long)expected_pos;

		// the record sizes, so you can tell a bad offset from an actual io problem
		ss << " beg_run0=" << (long long)beg_run0
		   << " run_byte_size=" << (long long)run_byte_size;

		// the counter out of the header, read without going through get_nruns() - that one
		// throws on a bad stream, and something that throws while explaining an error is no
		// help at all.
		long long nruns = -1, fsize = -1;
		{
			ifstream probe(filename, ios_base::in | ios_base::binary);
			if (probe.good())
			{
				probe.seekg(0, ios_base::end);
				fsize = (long long)probe.tellg();
				probe.clear();
				probe.seekg(sizeof(std::int64_t) * 0, ios_base::beg);
				std::int64_t tmp = -1;
				probe.read(reinterpret_cast<char*>(&tmp), sizeof(tmp));
				if (probe)
					nruns = (long long)tmp;
			}
		}
		ss << " file_size=" << (fsize < 0 ? string("?") : to_string(fsize));
		ss << " header_nruns=" << (nruns < 0 ? string("?") : to_string(nruns));

		// do the comparison here instead of making whoever reads this do it. if the file is
		// smaller than the offset we asked for then the file is short, and the question is who
		// didnt write it - not why the read failed.
		if ((fsize >= 0) && (expected_pos >= 0))
		{
			ss << " short_by=";
			if (expected_pos < fsize) ss << "0(offset is within the file)";
			else ss << (long long)(expected_pos - fsize + 1);
		}
		if ((fsize >= 0) && (run_byte_size > 0) && (beg_run0 > 0))
		{
			long long records = (long long)((fsize - beg_run0) / run_byte_size);
			ss << " records_in_file=" << records;
			if ((nruns >= 0) && (records != nruns))
				ss << " MISMATCH(counter says " << nruns << ", file holds " << records << ")";
		}

		// only report errno when badbit is set, meaning the os actually refused. a short read
		// sets eof/fail without touching errno, so printing it there just shows whatever
		// unrelated call failed last - "errno=2(No such file or directory)" against a file that
		// is obviously open and readable, which is worse than saying nothing.
		if (buf_stream.bad() && (saved_errno != 0))
			ss << " errno=" << saved_errno << "(" << strerror(saved_errno) << ")";
		ss << "]";
	}
	catch (...)
	{
		// if this throws while explaining an error you lose the original error completely
		return "RunStorage::" + where + "() stream not good [diagnostic unavailable]";
	}
	return ss.str();
}

/**
 * @brief Run storage.
 *
 * @param _filename Description.
 *
 * @return Description.
 */
// zero beg_run0 and the byte sizes here instead of leaving it to reset()/init_restart(). they
// feed the get_stream_pos() math, so before they were set up we were making an offset out of
// uninitialized memory. the new error message caught this right away and printed
// "beg_run0=49807605984", which would only ever confuse whoever read it.
RunStorage::RunStorage(const string &_filename) :filename(_filename), beg_run0(0),
	run_byte_size(0), run_par_byte_size(0), run_data_byte_size(0)
{
}

/**
 * @brief Reset.
 *
 * @param _par_names Description.
 * @param _obs_names Description.
 * @param _filename Description.
 */
void RunStorage::reset(const vector<string> &_par_names, const vector<string> &_obs_names, const string &_filename)
{
	par_names = _par_names;
	obs_names = _obs_names;
	// a file needs to exist before it can be opened it with read and write
	// permission.   So open it with write permission to crteate it, close
	// and then reopen it with read and write permission.
	if (_filename.size() > 0)
	{
		filename = _filename;
	}
	if (buf_stream.is_open())
	{
		buf_stream.close();
	}

	/*if ((pest_utils::check_exist_in(filename)) || (pest_utils::check_exist_out(filename)))
    {
	    int flag = remove(filename.c_str());
	    if (flag != 0)
        {
	        throw runtime_error("RunStorage::reset(): error removing existing file '"+filename+"'");

        }
    }*/

	buf_stream.open(filename.c_str(), ios_base::out |  ios_base::binary);
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("reset"));
    }
	buf_stream.close();

	// install a large stream buffer before (re)opening so that the bulk header
	// write (serialized par/obs names) drains to disk in a handful of syscalls
	// instead of thousands of tiny default-buffer flushes.  pubsetbuf must be
	// called while the stream is closed to take effect on the next open.
	io_buffer.resize(IO_BUFFER_SIZE);
	buf_stream.rdbuf()->pubsetbuf(io_buffer.data(), io_buffer.size());

	buf_stream.open(filename.c_str(), ios_base::out | ios_base::in | ios_base::binary);
	//assert(buf_stream.good() == true);
	if (!buf_stream)
	{
		throw runtime_error(stream_diagnostic("reset"));
	}
	// calculate the number of bytes required to store parameter names
	vector<int8_t> serial_pnames(Serialization::serialize(par_names));
	std::int64_t p_name_size_64 = serial_pnames.size() * sizeof(char);
	// calculate the number of bytes required to store observation names
	vector<int8_t> serial_onames(Serialization::serialize(obs_names));
	std::int64_t o_name_size_64 = serial_onames.size() * sizeof(char);
	// calculate the number of bytes required to store a model run
	run_par_byte_size = par_names.size() * sizeof(double);
	run_data_byte_size = run_par_byte_size + obs_names.size() * sizeof(double);
	//compute the amount of memory required to store a single model run
	// run_byte_size = size of run_status + size of info_txt + size of info_value + size of parameter oand observation data
	run_byte_size =  sizeof(std::int8_t) + info_txt_length*sizeof(char) * sizeof(double) + run_data_byte_size;
	std::int64_t  run_size_64 = run_byte_size;
	beg_run0 = 4 * sizeof(std::int64_t) + serial_pnames.size() + serial_onames.size();
	std::int64_t n_runs_64=0;
	// write header to file
	buf_stream.seekp(0, ios_base::beg);
	buf_stream.write((char*) &n_runs_64, sizeof(n_runs_64));
	buf_stream.write((char*) &run_size_64, sizeof(run_size_64));
	buf_stream.write((char*) &p_name_size_64, sizeof(p_name_size_64));
	buf_stream.write((char*) &o_name_size_64, sizeof(o_name_size_64));
	buf_stream.write((char*)serial_pnames.data(), serial_pnames.size());
	buf_stream.write((char*)serial_onames.data(), serial_onames.size());
	//add flag for double buffering
	std::int8_t buf_status = 0;
	int end_of_runs = get_nruns();
	buf_stream.seekp(get_stream_pos(end_of_runs), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&buf_status), sizeof(buf_status));
	buf_stream.flush();
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("reset", end_of_runs, get_stream_pos(end_of_runs)));
    }
}


/**
 * @brief Init restart.
 *
 * @param _filename Description.
 */
void RunStorage::init_restart(const std::string &_filename)
{
	filename = _filename;
	par_names.clear();
	obs_names.clear();

	if (buf_stream.is_open())
	{
		buf_stream.close();
	}

	buf_stream.open(filename.c_str(), ios_base::out | ios_base::in | ios_base::binary | ios_base::ate);
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("init_restart"));
    }
	// read header
	buf_stream.seekg(0, ios_base::beg);

	std::int64_t n_runs_64;
	buf_stream.read((char*) &n_runs_64, sizeof(n_runs_64));

	std::int64_t  run_size_64;
	buf_stream.read((char*) &run_size_64, sizeof(run_size_64));
	run_byte_size = run_size_64;

	std::int64_t p_name_size_64;
	buf_stream.read((char*) &p_name_size_64, sizeof(p_name_size_64));

	std::int64_t o_name_size_64;
	buf_stream.read((char*) &o_name_size_64, sizeof(o_name_size_64));

	vector<int8_t> serial_pnames;
	serial_pnames.resize(p_name_size_64);
	buf_stream.read((char *)serial_pnames.data(), serial_pnames.size());
	Serialization::unserialize(serial_pnames, par_names);

	vector<int8_t> serial_onames;
	serial_onames.resize(o_name_size_64);
	buf_stream.read((char *)serial_onames.data(), serial_onames.size());
	Serialization::unserialize(serial_onames, obs_names);

	beg_run0 = 4 * sizeof(std::int64_t) + serial_pnames.size() + serial_onames.size();
	run_par_byte_size = par_names.size() * sizeof(double);
	run_data_byte_size = run_par_byte_size + obs_names.size() * sizeof(double);

	//check buffer to see if a write was improperly terminated
	std::int8_t r_status = 0;
	std::int8_t buf_status = 0;
	std::int32_t buf_run_id = 0;

	int end_of_runs = get_nruns();
	buf_stream.seekg(get_stream_pos(end_of_runs), ios_base::beg);
	buf_stream.read(reinterpret_cast<char*>(&buf_status), sizeof(buf_status));
	if (buf_status == 1 || buf_status == 2)
	{
		buf_stream.read(reinterpret_cast<char*>(&buf_run_id), sizeof(buf_run_id));
		buf_stream.read(reinterpret_cast<char*>(&r_status), sizeof(r_status));
		check_rec_id(buf_run_id);
		size_t n_par = par_names.size();
		size_t n_obs = obs_names.size();
		vector<double> pars_vec(n_par, Parameters::no_data);
		vector<double> obs_vec(n_obs, Observations::no_data);

		buf_stream.read(reinterpret_cast<char*>(pars_vec.data()), n_par * sizeof(double));
		buf_stream.read(reinterpret_cast<char*>(obs_vec.data()), n_obs * sizeof(double));

		//write data
		buf_stream.seekp(get_stream_pos(buf_run_id), ios_base::beg);
		buf_stream.write(reinterpret_cast<char*>(&r_status), sizeof(r_status));
		//skip over info_txt and info_value fields
		buf_stream.seekp(sizeof(char)*info_txt_length + sizeof(double), ios_base::cur);
		buf_stream.write(reinterpret_cast<char*>(pars_vec.data()), pars_vec.size() * sizeof(double));
		buf_stream.write(reinterpret_cast<char*>(obs_vec.data()), obs_vec.size() * sizeof(double));
		buf_stream.flush();
		//reset flag for buffer at end of file to 0 to signal it is no longer relevant
		buf_status = 0;
		buf_stream.seekp(get_stream_pos(end_of_runs), ios_base::beg);
		buf_stream.write(reinterpret_cast<char*>(&buf_status), sizeof(buf_status));
		buf_stream.flush();
	}
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("init_restart", end_of_runs, get_stream_pos(end_of_runs)));
    }
}

/**
 * @brief Get nruns.
 *
 * @return Description.
 */
int RunStorage::get_nruns()
{

	// just read the run count from the file header.  do NOT tellg()/restore the position: tellg() is
	// an input op, and calling it right after the stream was used for output (which is constant here)
	// trips libstdc++'s "you must reposition when switching output<->input" rule and fails - clear()
	// does not satisfy that boundary.  clear()+seekg(0) IS a clean reposition; every caller absolute-
	// seeks before its own next I/O, so none rely on the position being preserved across this call.
	buf_stream.clear();
	buf_stream.seekg(0, ios_base::beg);
	std::int64_t n_runs_64;
	buf_stream.read((char*) &n_runs_64, sizeof(n_runs_64));
	int n_runs = n_runs_64;
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("get_nruns", -1, 0));
    }
    if (n_runs < 0)
    {
        cout << "RunStorage::get_nruns(): warning: nruns < 0: " << n_runs << endl;
    }

	return n_runs;
}

/**
 * @brief Get num good runs.
 *
 * @return Description.
 */
int RunStorage::get_num_good_runs()
{
	int n_ok = 0;
	int n_runs = get_nruns();
	for (int id = 0; id<n_runs; ++id)
	{
		std::int8_t tmp_r_status = get_run_status_native(id);
		if (tmp_r_status > 0)
		{
			++n_ok;
		}
	}
	return n_ok;
}
/**
 * @brief Increment nruns.
 *
 * @return Description.
 */
int RunStorage::increment_nruns()
{
	// reset stale error state before the read-then-write of the header on the shared fstream; the
	// seekp(0) below provides the required boundary when switching from reading to writing.
	buf_stream.clear();
	buf_stream.seekg(0, ios_base::beg);
	std::int64_t n_runs_64;
	buf_stream.read((char*) &n_runs_64, sizeof(n_runs_64));
	++n_runs_64;
	buf_stream.seekp(0, ios_base::beg);
	buf_stream.write((char*) &n_runs_64, sizeof(n_runs_64));
	int n_runs = n_runs_64;
	buf_stream.flush();
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("increment_nruns", -1, 0));
    }
	return n_runs;
}
const std::vector<string>& RunStorage::get_par_name_vec()const
{
	return par_names;
}

const std::vector<string>& RunStorage::get_obs_name_vec()const
{
	return obs_names;
}

/**
 * @brief Get stream pos.
 *
 * @param run_id Description.
 *
 * @return Description.
 */
streamoff RunStorage::get_stream_pos(int run_id)
{
	streamoff pos = beg_run0 + run_byte_size*run_id;
	return pos;
}

/**
 * @brief Add run.
 *
 * @param model_pars Description.
 * @param info_txt Description.
 * @param info_value Description.
 *
 * @return Description.
 */
 int RunStorage::add_run(const vector<double> &model_pars, const string &info_txt, double info_value)
 {
	std::int8_t r_status = 0;
	int run_id = increment_nruns() - 1;
	vector<char> info_txt_buf;
	info_txt_buf.resize(info_txt_length, '\0');
	copy_n(info_txt.begin(), min(info_txt.size(), size_t(info_txt_length)-1) , info_txt_buf.begin());
	buf_stream.seekp(get_stream_pos(run_id), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	buf_stream.write(reinterpret_cast<char*>(info_txt_buf.data()), sizeof(char)*info_txt_buf.size());
	buf_stream.write(reinterpret_cast<char*>(&info_value), sizeof(double));
	buf_stream.write(reinterpret_cast<const char*>(&model_pars[0]), model_pars.size()*sizeof(double));
	//and the observation block, explicitly as no_data. Leaving it unwritten made it a hole in
	//the file, and a hole reads back as ZEROS - a perfectly plausible observation value, so a
	//record read before its run completed handed back numbers nothing could tell were absent.
	vector<double> no_data_obs(obs_names.size(), Observations::no_data);
	buf_stream.write(reinterpret_cast<const char*>(no_data_obs.data()), no_data_obs.size()*sizeof(double));
	//add flag for double buffering
	std::int8_t buf_status = 0;
	int end_of_runs = get_nruns();
	buf_stream.seekp(get_stream_pos(end_of_runs), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&buf_status), sizeof(buf_status));
	buf_stream.flush();
     if (!buf_stream)
     {
         throw runtime_error(stream_diagnostic("add_run", end_of_runs, get_stream_pos(end_of_runs)));
     }
	return run_id;
 }

/**
 * @brief Add run.
 *
 * @param model_pars Description.
 * @param info_txt Description.
 * @param info_value Description.
 *
 * @return Description.
 */
 int RunStorage::add_run(const Eigen::VectorXd &model_pars, const string &info_txt, double info_value)
 {
	std::int8_t r_status = 0;
	int run_id = increment_nruns() - 1;
	vector<char> info_txt_buf;
	info_txt_buf.resize(info_txt_length, '\0');
	copy_n(info_txt.begin(), min(info_txt.size(), size_t(info_txt_length)-1) , info_txt_buf.begin());
	buf_stream.seekp(get_stream_pos(run_id), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	buf_stream.write(reinterpret_cast<char*>(info_txt_buf.data()), sizeof(char)*info_txt_buf.size());
	buf_stream.write(reinterpret_cast<char*>(&info_value), sizeof(double));
	buf_stream.write(reinterpret_cast<const char*>(&model_pars(0)), model_pars.size()*sizeof(model_pars(0)));
	//and the observation block as no_data, for the reason given in the overload above
	vector<double> no_data_obs(obs_names.size(), Observations::no_data);
	buf_stream.write(reinterpret_cast<const char*>(no_data_obs.data()), no_data_obs.size()*sizeof(double));
	//add flag for double buffering
	std::int8_t buf_status = 0;
	int end_of_runs = get_nruns();
	buf_stream.seekp(get_stream_pos(end_of_runs), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&buf_status), sizeof(buf_status));
	buf_stream.flush();
     if (!buf_stream)
     {
         throw runtime_error(stream_diagnostic("add_run", end_of_runs, get_stream_pos(end_of_runs)));
     }
	return run_id;
 }


/**
 * @brief Add run.
 *
 * @param pars Description.
 * @param info_txt Description.
 * @param info_value Description.
 *
 * @return Description.
 */
int RunStorage::add_run(const Parameters &pars, const string &info_txt, double info_value)
{
	vector<double> data(pars.get_data_vec(par_names));
	int run_id = add_run(data, info_txt, info_value);
	return run_id;
}

/**
 * @brief Copy.
 *
 * @param rhs_rs Description.
 */
void RunStorage::copy(const RunStorage &rhs_rs)
{
	if (buf_stream.is_open())
	{
		buf_stream.close();
	}
	// std::ofstream::trunc will delete the file if it already exist
	// a file needs to exist before it can be opened it with read and write
	// permission.   So open it with write permission to create it, close
	// and then reopen it with read and write permission.
	buf_stream.open(filename.c_str(), ios_base::out | ios_base::binary | std::ofstream::trunc);
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("copy"));
    }
	buf_stream.close();
	buf_stream.open(filename.c_str(), ios_base::out | ios_base::in | ios_base::binary);
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("copy"));
    }

	// copy rhs runstorage information
	std::streampos rhs_initial_pos = rhs_rs.buf_stream.tellg();
	rhs_rs.buf_stream.seekg(0, ios_base::beg);
	buf_stream << rhs_rs.buf_stream.rdbuf();
	rhs_rs.buf_stream.seekg(rhs_initial_pos);
	beg_run0 = rhs_rs.beg_run0;
	run_byte_size = rhs_rs.run_byte_size;
	run_par_byte_size = rhs_rs.run_par_byte_size;
	run_data_byte_size = rhs_rs.run_par_byte_size;
	par_names = rhs_rs.par_names;
	obs_names = rhs_rs.obs_names;
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("copy"));
    }
}

/**
 * @brief Update run.
 *
 * @param run_id Description.
 * @param pars Description.
 * @param obs Description.
 */
void RunStorage::update_run(int run_id, const Parameters &pars, const Observations &obs)
{

	//set run status flag to complete
	std::int8_t r_status = 1;
	check_rec_id(run_id);
	vector<double> par_data(pars.get_data_vec(par_names));
	vector<double> obs_data(obs.get_data_vec(obs_names));
	//write data to buffer at end of file and set buffer flag to 1
	std::int8_t buf_status = 0;
	std::int32_t buf_run_id = run_id;
	int end_of_runs = get_nruns();
	buf_stream.seekp(get_stream_pos(end_of_runs), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&buf_status), sizeof(buf_status));
	buf_stream.write(reinterpret_cast<char*>(&buf_run_id), sizeof(buf_run_id));
	buf_stream.write(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	buf_stream.write(reinterpret_cast<char*>(par_data.data()), par_data.size() * sizeof(double));
	buf_stream.write(reinterpret_cast<char*>(obs_data.data()), obs_data.size() * sizeof(double));
	buf_status = 1;
	buf_stream.seekp(get_stream_pos(end_of_runs), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&buf_status), sizeof(buf_status));
	buf_stream.flush();
	//write data
	buf_stream.seekp(get_stream_pos(run_id), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	//skip over info_txt and info_value fields
	buf_stream.seekp(sizeof(char)*info_txt_length+sizeof(double), ios_base::cur);
	buf_stream.write(reinterpret_cast<char*>(par_data.data()), par_data.size() * sizeof(double));
	buf_stream.write(reinterpret_cast<char*>(obs_data.data()), obs_data.size() * sizeof(double));
	buf_stream.flush();
	//reset flag for buffer at end of file to 0 to signal it is no longer relevant
	buf_status = 0;
	buf_stream.seekp(get_stream_pos(end_of_runs), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&buf_status), sizeof(buf_status));
	buf_stream.flush();
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("update_run", run_id, get_stream_pos(run_id)));
    }
}


/**
 * @brief Update run.
 *
 * @param run_id Description.
 * @param obs Description.
 */
void RunStorage::update_run_partial(int run_id, const Observations &obs)
{
	check_rec_id(run_id);
	vector<double> obs_data(obs.get_data_vec(obs_names));
	size_t n_pars = par_names.size();

	buf_stream.seekp(get_stream_pos(run_id), ios_base::beg);
	// step OVER the status byte rather than writing it - that is the whole point
	buf_stream.seekp(sizeof(std::int8_t), ios_base::cur);
	//skip over info_txt and info_value fields
	buf_stream.seekp(sizeof(char)*info_txt_length + sizeof(double), ios_base::cur);
	//skip over parameter section
	buf_stream.seekp(n_pars * sizeof(double), ios_base::cur);
	buf_stream.write(reinterpret_cast<char*>(obs_data.data()), obs_data.size() * sizeof(double));
	buf_stream.flush();
	if (!buf_stream)
	{
		throw runtime_error(stream_diagnostic("update_run_partial", run_id, get_stream_pos(run_id)));
	}
}

void RunStorage::update_run(int run_id, const Observations &obs)
{

	//set run status flag to complete
	std::int8_t r_status = 1;
	check_rec_id(run_id);
	vector<double> obs_data(obs.get_data_vec(obs_names));
	size_t n_pars = par_names.size();

	//write data to buffer at end of file and set buffer flag to 1
	std::int8_t buf_status = 0;
	std::int32_t buf_run_id = run_id;
	int end_of_runs = get_nruns();
	buf_stream.seekp(get_stream_pos(end_of_runs), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&buf_status), sizeof(buf_status));
	buf_stream.write(reinterpret_cast<char*>(&buf_run_id), sizeof(buf_run_id));
	buf_stream.write(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	//skip over parameter section
	buf_stream.seekp(n_pars * sizeof(double), ios_base::cur);
	buf_stream.write(reinterpret_cast<char*>(obs_data.data()), obs_data.size() * sizeof(double));
	buf_status = 1;
	buf_stream.seekp(get_stream_pos(end_of_runs), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&buf_status), sizeof(buf_status));
	buf_stream.flush();

	//write data to main part of file
	buf_stream.seekp(get_stream_pos(run_id), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	//skip over info_txt and info_value fields
	buf_stream.seekp(sizeof(char)*info_txt_length + sizeof(double), ios_base::cur);
	//skip over parameter section
	buf_stream.seekp(n_pars * sizeof(double), ios_base::cur);
	buf_stream.write(reinterpret_cast<char*>(obs_data.data()), obs_data.size() * sizeof(double));
	buf_stream.flush();
	//reset flag for buffer at end of file to 0 to signal it is no longer relevant
	buf_status = 0;
	buf_stream.seekp(get_stream_pos(end_of_runs), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&buf_status), sizeof(buf_status));
	buf_stream.flush();
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("update_run", run_id, get_stream_pos(run_id)));
    }
}

/**
 * @brief Update run.
 *
 * @param run_id Description.
 * @param serial_data Description.
 */
void RunStorage::update_run(int run_id, const vector<char> serial_data)
{

	//set run status flag to complete
	std::int8_t r_status = 1;
	check_rec_size(serial_data);
	check_rec_id(run_id);
	//write data to buffer at end of file and set buffer flag to 2
	std::int8_t buf_status = 0;
	std::int32_t buf_run_id = run_id;
	int end_of_runs = get_nruns();
	buf_stream.seekp(get_stream_pos(end_of_runs), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&buf_status), sizeof(buf_status));
	buf_stream.write(reinterpret_cast<char*>(&buf_run_id), sizeof(buf_run_id));
	buf_stream.write(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	buf_stream.write(serial_data.data(), serial_data.size());
	buf_status = 2;
	buf_stream.seekp(get_stream_pos(end_of_runs), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&buf_status), sizeof(buf_status));
	buf_stream.flush();
	//write data
	buf_stream.seekp(get_stream_pos(run_id), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	//skip over info_txt and info_value fields
	buf_stream.seekp(sizeof(char)*info_txt_length+sizeof(double), ios_base::cur);
	buf_stream.write(serial_data.data(), serial_data.size());
	buf_stream.flush();
	//reset flag for buffer at end of file to 0 to signal it is no longer relevant
	buf_status = 0;
	buf_stream.seekp(get_stream_pos(end_of_runs), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&buf_status), sizeof(buf_status));
	buf_stream.flush();
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("update_run", run_id, get_stream_pos(run_id)));
    }
}


/**
 * @brief Update run failed.
 *
 * @param run_id Description.
 */
void RunStorage::update_run_failed(int run_id)
{

	std::int8_t r_status = get_run_status_native(run_id);
	if (r_status < 1)
	{
		--r_status;
		check_rec_id(run_id);
		//update run status flag
		buf_stream.seekp(get_stream_pos(run_id), ios_base::beg);
		buf_stream.write(reinterpret_cast<char*>(&r_status), sizeof(r_status));
		buf_stream.flush();
	}
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("update_run_failed", run_id, get_stream_pos(run_id)));
    }
}

/**
 * @brief Set the failure count for a run record.
 *
 * Writes the negative of @p nfail into the run-status byte of the record
 * identified by @p run_id.  For example, if a run was attempted 3 times
 * before being declared permanently failed, @p nfail should be 3 and the
 * stored status will be -3.
 *
 * @param run_id  The zero-based index of the run record.
 * @param nfail   The number of failed attempts for this run.
 */
void RunStorage::set_run_nfailed(int run_id, int nfail)
{

	std::int8_t r_status = -nfail;
	check_rec_id(run_id);
	//update run status flag
	buf_stream.seekp(get_stream_pos(run_id), ios_base::beg);
	buf_stream.write(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	buf_stream.flush();
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("set_run_nfailed", run_id, get_stream_pos(run_id)));
    }
}

/**
 * @brief Get run status native.
 *
 * @param run_id Description.
 *
 * @return Description.
 */
std::int8_t RunStorage::get_run_status_native(int run_id)
{

	std::int8_t  r_status;
	check_rec_id(run_id);
	buf_stream.clear();  // reset stale error state before seeking the shared read/write fstream: libstdc++ can leave failbit after a write, making this seekg a silent no-op on good data (Linux-only). a genuinely bad read still throws below.
	buf_stream.seekg(get_stream_pos(run_id), ios_base::beg);
	buf_stream.read(reinterpret_cast<char*>(&r_status), sizeof(r_status));
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("get_run_status_native", run_id, get_stream_pos(run_id)));
    }
	return r_status;
}

/**
 * @brief Get run status.
 *
 * @param run_id Description.
 *
 * @return Description.
 */
int RunStorage::get_run_status(int run_id)
{
	int status = get_run_status_native(run_id);
	return status;
}

/**
 * @brief Get info.
 *
 * @param run_id Description.
 * @param run_status Description.
 * @param info_txt Description.
 * @param info_value Description.
 */
void RunStorage::get_info(int run_id, int &run_status, string &info_txt, double &info_value)
{
	// Bounds-checked like the other record accessors. Without this, an out-of-range id seeks
	// past the end of the storage file and READS ANYWAY - the stream reports the failure only
	// afterwards, and the info_txt read lands in a buffer sized for a record that is not
	// there. That surfaced as an access violation on windows rather than an error.
	check_rec_id(run_id);

	std::int8_t  r_status;
	vector<char> info_txt_buf;
	info_txt_buf.resize(info_txt_length, '\0');

	buf_stream.clear();  // reset stale error state before seeking the shared read/write fstream: libstdc++ can leave failbit after a write, making this seekg a silent no-op on good data (Linux-only). a genuinely bad read still throws below.
	buf_stream.seekg(get_stream_pos(run_id), ios_base::beg);
	buf_stream.read(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	buf_stream.read(reinterpret_cast<char*>(&info_txt_buf[0]), sizeof(char)*info_txt_length);
	buf_stream.read(reinterpret_cast<char*>(&info_value), sizeof(double));

	run_status = r_status;
	info_txt = info_txt_buf.data();
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("get_run_info", run_id, get_stream_pos(run_id)));
    }
}

/**
 * @brief Get run.
 *
 * @param run_id Description.
 * @param pars Description.
 * @param obs Description.
 * @param info_txt Description.
 * @param info_value Description.
 * @param clear_old Description.
 *
 * @return Description.
 */
int RunStorage::get_run(int run_id, Parameters &pars, Observations &obs, string &info_txt, double &info_value, bool clear_old)
{

	vector<double> par_data;
	vector<double> obs_data;

	int status = get_run(run_id, par_data, obs_data, info_txt, info_value);
	if (clear_old)
	{
	  pars.update(par_names, par_data);
	  obs.update(obs_names, obs_data);
	}
	else
	{
	  pars.update_without_clear(par_names, par_data);
	  obs.update_without_clear(obs_names, obs_data);
	}
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("get_run", run_id, get_stream_pos(run_id)));
    }
	return status;
}

/**
 * @brief Get run.
 *
 * @param run_id Description.
 * @param pars Description.
 * @param obs Description.
 * @param clear_old Description.
 *
 * @return Description.
 */
int RunStorage::get_run(int run_id, Parameters &pars, Observations &obs, bool clear_old)
{
	string info_txt;
	double info_value;
	return get_run(run_id, pars, obs, info_txt, info_value, clear_old);
}

/**
 * @brief Get run.
 *
 * @param run_id Description.
 * @param pars Description.
 * @param npars Description.
 * @param obs Description.
 * @param nobs Description.
 * @param info_txt Description.
 * @param info_value Description.
 *
 * @return Description.
 */
int RunStorage::get_run(int run_id, double *pars, size_t npars, double *obs, size_t nobs, string &info_txt, double &info_value)
{
    if (!buf_stream.good())
    {
        throw runtime_error(stream_diagnostic("get_run", run_id, get_stream_pos(run_id)));
    }
	std::int8_t r_status;
	vector<char> info_txt_buf;
	info_txt_buf.resize(info_txt_length, '\0');

	check_rec_id(run_id);

	size_t p_size = par_names.size();
	size_t o_size = obs_names.size();

	//assert(npars == p_size);
	//assert(nobs == o_size);

	if (npars != p_size) {
		throw(PestIndexError("RunStorage::get_run: parameter dimension in incorrect"));
	}
	if (nobs != o_size) {
		throw(PestIndexError("RunStorage::get_run: observation dimension in incorrect"));
	}

	p_size = min(p_size, npars);
	o_size = min(o_size, nobs);
	buf_stream.clear();  // reset stale error state before seeking the shared read/write fstream: libstdc++ can leave failbit after a write, making this seekg a silent no-op on good data (Linux-only). a genuinely bad read still throws below.
	buf_stream.seekg(get_stream_pos(run_id), ios_base::beg);
	buf_stream.read(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	buf_stream.read(reinterpret_cast<char*>(&info_txt_buf[0]), sizeof(char)*info_txt_length);
	buf_stream.read(reinterpret_cast<char*>(&info_value), sizeof(double));
	buf_stream.read(reinterpret_cast<char*>(pars), p_size * sizeof(double));
	buf_stream.read(reinterpret_cast<char*>(obs), o_size * sizeof(double));
	int status = r_status;
	info_txt = info_txt_buf.data();
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("get_run", run_id, get_stream_pos(run_id)));
    }
	return status;
}

/**
 * @brief Get run.
 *
 * @param run_id Description.
 * @param pars_vec Description.
 * @param obs_vec Description.
 * @param info_txt Description.
 * @param info_value Description.
 *
 * @return Description.
 */
int RunStorage::get_run(int run_id, vector<double> &pars_vec, vector<double> &obs_vec, string &info_txt, double &info_value)
{
	std::int8_t  r_status;
	vector<char> info_txt_buf;
	info_txt_buf.resize(info_txt_length, '\0');

	size_t n_par = par_names.size();
	size_t n_obs = obs_names.size();

	pars_vec.resize(n_par);
	obs_vec.resize(n_obs);

	check_rec_id(run_id);

	buf_stream.clear();  // reset stale error state before seeking the shared read/write fstream: libstdc++ can leave failbit after a write, making this seekg a silent no-op on good data (Linux-only). a genuinely bad read still throws below.
	buf_stream.seekg(get_stream_pos(run_id), ios_base::beg);
	buf_stream.read(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	buf_stream.read(reinterpret_cast<char*>(&info_txt_buf[0]), sizeof(char)*info_txt_length);
	buf_stream.read(reinterpret_cast<char*>(&info_value), sizeof(double));
	buf_stream.read(reinterpret_cast<char*>(&pars_vec[0]), n_par * sizeof(double));
	buf_stream.read(reinterpret_cast<char*>(&obs_vec[0]), n_obs * sizeof(double));
	int status = r_status;
	info_txt = info_txt_buf.data();
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("get_run", run_id, get_stream_pos(run_id)));
    }
	return status;
}

/**
 * @brief Get run.
 *
 * @param run_id Description.
 * @param pars_vec Description.
 * @param obs_vec Description.
 *
 * @return Description.
 */
int RunStorage::get_run(int run_id, vector<double> &pars_vec, vector<double> &obs_vec)
{
	string info_txt;
	double info_value;
	return get_run(run_id, pars_vec, obs_vec, info_txt, info_value);
}

/**
 * @brief Get run.
 *
 * @param run_id Description.
 * @param pars Description.
 * @param npars Description.
 * @param obs Description.
 * @param nobs Description.
 *
 * @return Description.
 */
int RunStorage::get_run(int run_id, double *pars, size_t npars, double *obs, size_t nobs)
{
	string info_txt;
	double info_value;
	return get_run(run_id, pars, npars, obs, nobs, info_txt, info_value);
}

/**
 * @brief Get serial pars.
 *
 * @param run_id Description.
 *
 * @return Description.
 */
vector<char> RunStorage::get_serial_pars(int run_id)
{

	check_rec_id(run_id);
	std::int8_t r_status;

	vector<char> serial_data;
	serial_data.resize(run_par_byte_size);
	buf_stream.clear();  // reset stale error state before seeking the shared read/write fstream: libstdc++ can leave failbit after a write, making this seekg a silent no-op on good data (Linux-only). a genuinely bad read still throws below.
	buf_stream.seekg(get_stream_pos(run_id), ios_base::beg);
	buf_stream.seekg(sizeof(r_status)+sizeof(char)*info_txt_length+sizeof(double), ios_base::cur);
	buf_stream.read(serial_data.data(), serial_data.size());
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("get_serial_pars", run_id, get_stream_pos(run_id)));
    }
	return serial_data;
}

/**
 * @brief Get parameters.
 *
 * @param run_id Description.
 * @param pars Description.
 *
 * @return Description.
 */
int  RunStorage::get_parameters(int run_id, Parameters &pars)
{
	std::int8_t r_status;
	vector<char> info_txt_buf;
	info_txt_buf.resize(info_txt_length, '\0');
	double info_value;

	check_rec_id(run_id);

	size_t n_par = par_names.size();
	vector<double> par_data;
	par_data.resize(n_par);
	buf_stream.clear();  // reset stale error state before seeking the shared read/write fstream: libstdc++ can leave failbit after a write, making this seekg a silent no-op on good data (Linux-only). a genuinely bad read still throws below.
	buf_stream.seekg(get_stream_pos(run_id), ios_base::beg);
	buf_stream.read(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	buf_stream.read(reinterpret_cast<char*>(&info_txt_buf[0]), sizeof(char)*info_txt_length);
	buf_stream.read(reinterpret_cast<char*>(&info_value), sizeof(double));

	buf_stream.read(reinterpret_cast<char*>(par_data.data()), n_par*sizeof(double));
	pars.update(par_names, par_data);
	int status = r_status;
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("get_parameters", run_id, get_stream_pos(run_id)));
    }
	return status;
}


/**
 * @brief Get observations.
 *
 * @param run_id Description.
 * @param obs Description.
 *
 * @return Description.
 */
int  RunStorage::get_observations(int run_id, Observations &obs)
{

	std::int8_t r_status;
	vector<char> info_txt_buf;
	info_txt_buf.resize(info_txt_length, '\0');
	double info_value;

	check_rec_id(run_id);

	size_t n_par = par_names.size();
	size_t n_obs = obs_names.size();
	vector<double> obs_data;
	obs_data.resize(n_obs);
	buf_stream.clear();  // reset stale error state before seeking the shared read/write fstream: libstdc++ can leave failbit after a write, making this seekg a silent no-op on good data (Linux-only). a genuinely bad read still throws below.
	buf_stream.seekg(get_stream_pos(run_id), ios_base::beg);
	buf_stream.read(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	buf_stream.read(reinterpret_cast<char*>(&info_txt_buf[0]), sizeof(char)*info_txt_length);
	buf_stream.read(reinterpret_cast<char*>(&info_value), sizeof(double));
	buf_stream.seekg(n_par*sizeof(double), ios_base::cur);
	buf_stream.read(reinterpret_cast<char*>(obs_data.data()), n_obs*sizeof(double));
	int status = r_status;
	obs.update(obs_names, obs_data);
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("get_observations", run_id, get_stream_pos(run_id)));
    }
	return status;
}


/**
 * @brief Get observations vec.
 *
 * @param run_id Description.
 * @param obs_data Description.
 *
 * @return Description.
 */
int  RunStorage::get_observations_vec(int run_id, vector<double> &obs_data)
{

	std::int8_t r_status;
	vector<char> info_txt_buf;
	info_txt_buf.resize(info_txt_length, '\0');
	double info_value;

	check_rec_id(run_id);

	size_t n_par = par_names.size();
	size_t n_obs = obs_names.size();
	obs_data.resize(n_obs);
	buf_stream.clear();  // reset stale error state before seeking the shared read/write fstream: libstdc++ can leave failbit after a write, making this seekg a silent no-op on good data (Linux-only). a genuinely bad read still throws below.
	buf_stream.seekg(get_stream_pos(run_id), ios_base::beg);
	buf_stream.read(reinterpret_cast<char*>(&r_status), sizeof(r_status));
	buf_stream.read(reinterpret_cast<char*>(&info_txt_buf[0]), sizeof(char)*info_txt_length);
	buf_stream.read(reinterpret_cast<char*>(&info_value), sizeof(double));
	buf_stream.seekg(n_par*sizeof(double), ios_base::cur);
	buf_stream.read(reinterpret_cast<char*>(obs_data.data()), n_obs*sizeof(double));
	int status = r_status;
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("get_observations_vec", run_id, get_stream_pos(run_id)));
    }
	return status;
}

/**
 * @brief Free memory.
 */
void RunStorage::free_memory()
{

	if (buf_stream.is_open()) {
		buf_stream.close();
		remove(filename.c_str());
	}
    if (!buf_stream)
    {
        throw runtime_error(stream_diagnostic("free_memory"));
    }
}

/**
 * @brief Check rec size.
 *
 * @param serial_data Description.
 */
void RunStorage::check_rec_size(const vector<char> &serial_data) const
{
	if (serial_data.size() != run_data_byte_size)
	{
		throw PestError("Error in RunStorage routine.  Size of serial data is different from what is expected");
	}
}

/**
 * @brief Check rec id.
 *
 * @param run_id Description.
 */
void RunStorage::check_rec_id(int run_id)
{
	int n_runs = get_nruns();
	if ( run_id + 1 > n_runs)
	{
		ostringstream msg;
		msg << "Error in RunStorage routine: run id = " << run_id << " is not valid.  Valid values are 0 to " << n_runs - 1 << endl;
		throw PestError(msg.str());
	}
}


/**
 * @brief Print run summary.
 *
 * @param fout Description.
 */
void RunStorage::print_run_summary(std::ostream &fout)
{
	int nruns = get_nruns();
	fout << "nruns = " << nruns << endl;
	int status;
	string info_text;
	double info_value;
	for (int irun = 0; irun < nruns; ++irun)
	{
		get_info(irun, status, info_text, info_value);
		fout << "run_id=" << irun << "  :status=" << status << "  :info_text=" << info_text << "  :info_value=" << info_value << endl;
	}

}

/**
 * @brief Open a run-storage file by name and print a concise run-count summary.
 *
 * Reads the header and per-run status bytes of @p filename and writes the total
 * number of recorded runs along with the successful and failed counts to
 * @p fout.  If the file does not exist this is a no-op (so callers can invoke it
 * unconditionally); if the file exists but cannot be parsed a warning is printed.
 *
 * @param filename  The run-storage file to summarize (e.g. *case.allruns.rns*).
 * @param fout      The output stream to write the summary to.
 */
void RunStorage::print_persistent_summary(const std::string &filename, std::ostream &fout)
{
	{
		ifstream f(filename, ios::binary);
		if (!f.good())
			return;
	}
	try
	{
		RunStorage rs("");
		rs.init_restart(filename);
		int nruns = rs.get_nruns();
		int ngood = rs.get_num_good_runs();
		int nfailed = 0;
		for (int irun = 0; irun < nruns; ++irun)
		{
			if (rs.get_run_status(irun) < 0)
				++nfailed;
		}
		fout << endl << "  all-runs storage file '" << filename << "':" << endl;
		fout << "    total completed runs recorded: " << nruns << endl;
		fout << "    successful runs: " << ngood << endl;
		fout << "    failed runs: " << nfailed << endl;
	}
	catch (const exception &e)
	{
		fout << endl << "  warning: unable to read all-runs storage file '" << filename << "': " << e.what() << endl;
	}
}

/**
 * @brief Export diff to text file.
 *
 * @param in1_filename Description.
 * @param in2_filename Description.
 * @param out_filename Description.
 */
void RunStorage::export_diff_to_text_file(const std::string &in1_filename, const std::string &in2_filename, const std::string &out_filename)
{
	RunStorage rs1("");
	RunStorage rs2("");
	rs1.init_restart(in1_filename);
	rs2.init_restart(in2_filename);

	ofstream fout;
	fout.open(out_filename);
	fout.precision(numeric_limits<double>::digits10 + 1);

	vector<string> par_name_vec_1 = rs1.get_par_name_vec();
	vector<string> par_name_vec_2 = rs2.get_par_name_vec();
	vector<string> obs_name_vec_1 = rs1.get_obs_name_vec();
	vector<string> obs_name_vec_2 = rs2.get_obs_name_vec();

	auto npar = par_name_vec_1.size();
	auto nobs = obs_name_vec_1.size();
	fout << in1_filename << endl;
	fout << in2_filename << endl;
	fout << "npar = " << npar << endl;
	fout << "nobs = " << nobs << endl;
	//check parameter names are the same
	if (par_name_vec_1.size() != par_name_vec_2.size())
		fout << "Parameter name arrays sizes differ: " << par_name_vec_1.size() << ", " << par_name_vec_2.size() << endl;

	for (size_t ipar = 0; ipar < npar; ++ipar)
	{
		if (par_name_vec_1[ipar] != par_name_vec_2[ipar])
			fout << "Parameter name diff index =  " << ipar << ": " << par_name_vec_1[ipar] << ", " << par_name_vec_2[ipar] << endl;
	}

	//check observation names are the same
	if (obs_name_vec_1.size() != obs_name_vec_2.size())
		fout << "Observation name arrays sizes differ: " << obs_name_vec_1.size() << ", " << obs_name_vec_2.size() << endl;

	for (size_t iobs = 0; iobs < nobs; ++iobs)
	{
		if (obs_name_vec_1[iobs] != obs_name_vec_2[iobs])
			fout << "Observation  name diff index =  " << iobs << ": " << obs_name_vec_1[iobs] << ", " << obs_name_vec_2[iobs] << endl;
	}

	vector<double> pars_vec1;
	vector<double> pars_vec2;
	vector<double> obs_vec1;
	vector<double> obs_vec2;
	string info_txt1;
	string info_txt2;
	double info_value1;
	double info_value2;

	int nruns = rs1.get_nruns();
	fout << "nruns = " << nruns << endl;
	for (int irun = 0; irun < nruns; ++irun)
	{
		fout << "Run Id = " << irun << " --------------------" << endl;
		fout << "Run Description " << info_txt1 << endl;
		rs1.get_run(irun, pars_vec1, obs_vec1, info_txt1, info_value1);
		rs2.get_run(irun, pars_vec2, obs_vec2, info_txt2, info_value2);
		if (info_txt1 != info_txt2)
		{
			fout << "infotext diff at run id = " << irun <<  endl;
			fout << info_txt1 << endl;
			fout << info_txt2 << endl;
		}
		if (info_value1 != info_value2)
		{
			fout << "info_value diff at run id = " << irun << endl;
			fout << info_value1 << endl;
			fout << info_value2 << endl;
		}

		for (size_t ipar = 0; ipar < npar; ++ipar)
		{
			if (pars_vec1[ipar] != pars_vec2[ipar])
				fout << "Parameter diff at index =  " << ipar << ": " << "name = " << par_name_vec_1[ipar] << ":  " << pars_vec1[ipar] << ", " << pars_vec2[ipar] << endl;
		}

		for (size_t iobs = 0; iobs < nobs; ++iobs)
		{
			if (obs_vec1[iobs] != obs_vec2[iobs])
				fout << "Observation diff at index =  " << iobs << ": " << "name = " << obs_name_vec_1[iobs] << ":  " << obs_vec1[iobs] << ", " << obs_vec2[iobs] << endl;
		}
	}
	fout.close();
}


/**
 * @brief Destructor for .
 */
RunStorage::~RunStorage()
{
  //free_memory();
}

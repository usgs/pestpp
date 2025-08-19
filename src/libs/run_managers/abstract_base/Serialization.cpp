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

#include "network_wrapper.h"
#include <string>
#include <sstream>
#include <memory>
#include <cassert>
#include "Serialization.h"
#include "Transformable.h"
#include "utilities.h"

using namespace std;
using namespace pest_utils;


vector<int8_t> Serialization::serialize(int64_t data)
{
	vector<int8_t> buf;
	int64_t buf_sz = sizeof(data);
    buf.resize(buf_sz);
	size_t i_start = 0;
	w_memcpy_s(&buf[i_start], buf_sz, &data, buf_sz);
	return buf;
}


unsigned long Serialization::unserialize(const vector<int8_t> &buf, int64_t &data, unsigned long start_loc)
{
	assert(buf.size()-start_loc >= sizeof(data));
	w_memcpy_s(&data, sizeof(data), &buf[start_loc], sizeof(data));
	unsigned bytes_read = sizeof(data);
	return bytes_read;
}


vector<int8_t> Serialization::serialize(const Transformable &tr_data)
{
	vector<int8_t> buf;
	int64_t buf_sz = 0;
	int64_t names_buf_sz = 0;;
	int64_t data_buf_sz = 0;
	// calculate buffer size
	for (auto &b : tr_data)
	{
		names_buf_sz += b.first.size() + 1;
	}
	data_buf_sz = sizeof(double)*tr_data.size();
	buf_sz = sizeof(names_buf_sz)+names_buf_sz + sizeof(data_buf_sz)+data_buf_sz;
	// allocate space
	buf.resize(buf_sz, '\0');
	// build string with space deliminated names and array of numbers
	vector<int8_t> names;
	names.reserve(names_buf_sz);
	vector<double> values;
	for (auto &b : tr_data)
	{
		vector<int8_t> tmp_str = NetPackage::pack_string(b.first.begin(), b.first.end());
		names.insert(names.end(), tmp_str.begin(), tmp_str.end());
		names.push_back(' ');
		values.push_back(b.second);
	}

	int64_t n_rec = values.size();
	//write information to buffer
	size_t i_start = 0;
	w_memcpy_s(&buf[i_start], buf_sz - i_start, &names_buf_sz, sizeof(names_buf_sz));
	i_start += sizeof(names_buf_sz);
	w_memcpy_s(&buf[i_start], buf_sz - i_start, names.data(), names_buf_sz);
	i_start += names_buf_sz;
	w_memcpy_s(&buf[i_start], buf_sz - i_start, &n_rec, sizeof(n_rec));
	i_start += sizeof(n_rec);
	if (n_rec > 0)
	{
		w_memcpy_s(&buf[i_start], buf_sz - i_start, &values[0], data_buf_sz);
	}
	return buf;
}

vector<int8_t> Serialization::serialize(const vector< const Transformable*> tr_vec)
{
	vector<int8_t> buf;
	for (auto &i : tr_vec)
	{
		vector<int8_t> serial_data = serialize(*i);
		buf.insert(buf.end(), serial_data.begin(), serial_data.end());
	}
	return buf;
}

vector<int8_t> Serialization::serialize(const std::vector<Transformable*> &tr_vec)
{
	vector<const Transformable*> const_data_vec;
	for (auto &i : tr_vec)
	{
		const_data_vec.push_back(i);
	}
	return serialize(const_data_vec);
}

vector<int8_t> Serialization::serialize(const Parameters &pars, const Observations &obs)
{
	 vector<const Transformable*> tr_vec;
	 tr_vec.push_back(&pars);
	 tr_vec.push_back(&obs);
	 return serialize(tr_vec);
}

vector<int8_t> Serialization::serialize(const Parameters &pars, const vector<string> &par_names_vec, const Observations &obs, const vector<string> &obs_names_vec, double run_time)
{

	assert(pars.size() == par_names_vec.size());
	assert(obs.size() == obs_names_vec.size());
	vector<int8_t> serial_data;
	size_t npar = par_names_vec.size();
	size_t nobs = obs_names_vec.size();
	size_t par_buf_sz = npar * sizeof(double);
	size_t obs_buf_sz = nobs * sizeof(double);
	size_t run_time_sz = sizeof(double);
	serial_data.resize(par_buf_sz + obs_buf_sz + run_time_sz, Parameters::no_data);

	int8_t *buf = &serial_data[0];
	vector<double> par_data = pars.get_data_vec(par_names_vec);
	w_memcpy_s(buf, par_buf_sz, &par_data[0], par_data.size() * sizeof(double));

	vector<double> obs_data = obs.get_data_vec(obs_names_vec);

	w_memcpy_s(buf+par_buf_sz, obs_buf_sz, &obs_data[0], obs_data.size() * sizeof(double));
	w_memcpy_s(buf + par_buf_sz + obs_buf_sz, run_time_sz, &run_time, sizeof(double));

	return serial_data;
}

vector<int8_t> Serialization::serialize(const vector<string> &string_vec)
{
	vector<int8_t> serial_data;
	for (auto &i : string_vec)
	{
		vector<int8_t> tmp_str = NetPackage::pack_string(i.begin(),i.end());
		serial_data.insert(serial_data.end(), tmp_str.begin(), tmp_str.end());
		serial_data.push_back('\0');
	}
	return serial_data;
}

vector<int8_t> Serialization::serialize(const vector<vector<string>const*> &string_vec_vec)
{
	vector<int8_t> serial_data;
	unsigned long buf_sz = 0;
	for (auto &i : string_vec_vec)
	{
		vector<int8_t> buf_ser = serialize(*i);
		buf_sz = buf_ser.size();
		vector<int8_t> buf_size = serialize(buf_sz);
		serial_data.insert(serial_data.end(), buf_size.begin(), buf_size.end());
		serial_data.insert(serial_data.end(), buf_ser.begin(), buf_ser.end());
	}
	return serial_data;
}


unsigned long Serialization::unserialize(const std::vector<int8_t> &buf, Transformable &tr_data, unsigned long start_loc)
{
	// delete all existing items
	tr_data.clear();
	int64_t names_arg_sz;
	int64_t n_rec;
	int64_t bytes_read;
	size_t i_start = 0;
	// get size of names record
	i_start = start_loc;
	w_memcpy_s(&names_arg_sz, sizeof(names_arg_sz), buf.data()+i_start, sizeof(names_arg_sz));
	i_start += sizeof(names_arg_sz);
	unique_ptr<int8_t[]> names_buf(new int8_t[names_arg_sz]);
	w_memcpy_s(names_buf.get(), sizeof(int8_t)*names_arg_sz, buf.data() + i_start, sizeof(int8_t)*names_arg_sz);
	i_start += sizeof(int8_t)*names_arg_sz;
	w_memcpy_s(&n_rec, sizeof(n_rec),  buf.data()+i_start, sizeof(n_rec));
	i_start += sizeof(n_rec);
	// build transformable data set
	double *value_ptr = (double *) ( buf.data()+i_start);
	vector<string> names_vec;
	tokenize(strip_cp(NetPackage::extract_string(names_buf.get(), names_arg_sz)), names_vec);
	//tokenize(strip_cp(string(names_buf.get(), names_arg_sz)), names_vec);
	assert(names_vec.size() == n_rec);
	for (int64_t i = 0; i< n_rec; ++i)
	{
		tr_data.insert(names_vec[i], *(value_ptr+i));
	}
	bytes_read = i_start + n_rec * sizeof(double);
	return bytes_read;
}

unsigned long Serialization::unserialize(const std::vector<int8_t> &ser_data, std::vector<Transformable*> &tr_vec, unsigned long start_loc)
{
	unsigned i_tr=start_loc;
	unsigned i_char=0;
	unsigned long bytes_read = 0;
	unsigned long total_bytes_read = 0;


	while(i_tr < tr_vec.size() && i_char < ser_data.size())
	{
		bytes_read = Serialization::unserialize(ser_data, *tr_vec[i_tr], i_char);
		i_char +=  bytes_read;
		++i_tr;
		total_bytes_read += bytes_read;
	}
	return total_bytes_read;
}

unsigned long Serialization::unserialize(const std::vector<int8_t> &data, Parameters &pars, Observations &obs, unsigned long start_loc)
{
	 unsigned total_bytes_read = 0;
	 vector<Transformable*> tr_vec;
	 tr_vec.push_back(&pars);
	 tr_vec.push_back(&obs);
	 total_bytes_read = unserialize(data, tr_vec, start_loc);
	 return total_bytes_read;
}

unsigned long Serialization::unserialize(const vector<int8_t> &ser_data, vector<string> &string_vec, unsigned long start_loc, unsigned long max_read_bytes)
{
	unsigned total_bytes_read = 0;
	assert(start_loc < ser_data.size());
	total_bytes_read = min<unsigned long>((unsigned long)(ser_data.size()-start_loc), max_read_bytes);
	string tmp_str = NetPackage::extract_string(ser_data, start_loc, total_bytes_read);
	//string tmp_str(ser_data.begin()+start_loc, ser_data.begin()+start_loc+total_bytes_read);
	string delm;
	delm.push_back('\0');
	strip_ip(tmp_str, "both", delm);
	tokenize(tmp_str, string_vec, delm);
	return total_bytes_read;
}

unsigned long Serialization::unserialize(const vector<int8_t> &ser_data, Transformable &items, const vector<string> &names_vec, unsigned long start_loc)
{
	unsigned long total_bytes_read = 0;
	size_t vec_size = (ser_data.size() - start_loc) / sizeof(double);
	assert(vec_size >= names_vec.size());

	items.clear();
	double *value = 0;
	size_t iloc = start_loc;
	size_t n_read = 0;
	size_t n_val = names_vec.size();
	for (n_read=0; n_read<n_val; ++n_read)
	{
		value = (double*)&ser_data[iloc];
		items.insert(names_vec[n_read], *value);
		iloc += sizeof(double);
	}
	total_bytes_read = n_read * sizeof(double);
	return total_bytes_read;
}

unsigned long Serialization::unserialize(const vector<int8_t> &ser_data, Parameters &pars, const vector<string> &par_names, Observations &obs, const vector<string> &obs_names, double &run_time)
{
	unsigned long bytes_read = 0;

	bytes_read = unserialize(ser_data, pars, par_names, 0);
	bytes_read += unserialize(ser_data, obs, obs_names, bytes_read);
	w_memcpy_s(&run_time, sizeof(double), ser_data.data() + bytes_read, sizeof(double));
	return bytes_read;
}

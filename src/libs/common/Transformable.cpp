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
#include <ostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <fstream>
#include <iomanip>
#include <utility>
#include <cassert>
#include <memory>
#include <Eigen/Dense>
#include <utility>
#include <cmath>
#include <cfloat>
#include "Transformable.h"
#include "pest_error.h"
#include "utilities.h"


using std::string;
using std::map;
using std::ostream;
using std::endl;
using std::isnormal;
using namespace pest_utils;
using namespace Eigen;

const double Transformable::no_data = -9.99E99;

Transformable::Transformable(const Transformable &copyin) : items(copyin.items)
{
}

Transformable::Transformable(const Transformable &&copyin)
{
	items = std::move(copyin.items);
}

Transformable::Transformable(const Transformable &copyin, const vector<string> &copy_names)
{
	for (vector<string>::const_iterator b=copy_names.begin(), e=copy_names.end(); b!=e; ++b) {
		items[*b] = copyin.get_rec(*b);
	}
}

Transformable::Transformable(const vector<string> &names, const Eigen::VectorXd &values)
{
	assert(names.size() == values.size());
	if(names.size() != values.size())
	{
		throw PestIndexError("Transformable::Transformable(const vector<string> &names, Eigen::VectorXd &values)",
			"size of names vector does not match the size of the values vector");
	}
	size_t len = min(size_t(names.size()), size_t(values.size()));
	for (size_t i=0; i<len; ++i)
	{
		items[names[i]] = values(i);
	}
}


const Transformable& Transformable::operator=(const Transformable &rhs)
{
	items = rhs.items;
	return *this;
}

bool Transformable::operator==(const Transformable &rhs) const
{
	return items.size() == rhs.items.size()
		&& std::equal(items.begin(), items.end(),
					  rhs.items.begin());
}

bool Transformable::operator!=(const Transformable &rhs) const
{
	return !(*this == rhs);
}

Transformable& Transformable::operator+=(const Transformable &rhs)
{
	for(auto &i : items)
	{
		auto iter = rhs.items.find(i.first);
		assert(iter != rhs.items.end());
		i.second += iter->second;
	}
	return *this;
 }


Transformable& Transformable::operator-=(const Transformable &rhs)
{
	for(auto &i : items)
	{
		auto iter = rhs.items.find(i.first);
		assert(iter != rhs.items.end());
		i.second -= iter->second;
	}
	return *this;
 }

Transformable& Transformable::operator*=(double scale)
{
	for(auto &i : items)
	{
		i.second *= scale;
	}
	 return *this;
}

Transformable Transformable::operator-(const Transformable &rhs) const
{
	Transformable ret_val(*this);
	ret_val -= rhs;
	return ret_val;
}

Transformable Transformable::operator+(const Transformable &rhs) const
{
	Transformable ret_val(*this);
	ret_val += rhs;
	return ret_val;


}

Transformable Transformable::operator*(double scale) const
{
	Transformable ret_val(*this);
	ret_val *= scale;
	return ret_val;
}

double &Transformable::operator[](const string &name)
{
	return items[name];
}

vector<string> Transformable::get_notnormal_keys()
{
	vector<string> not_normal;
	for (auto &i : items)
	{
		if (isnormal(i.second))
			continue;
		not_normal.push_back(i.first);
	}
	return not_normal;
}

pair<Transformable::iterator,bool> Transformable::insert(const string &name, double value)
{
	pair<string, double> rec(name, value);
	return items.insert(rec);
}

pair<Transformable::iterator,bool>  Transformable::insert(const pair<string, double> &x)
{
	return items.insert(x);
}

void Transformable::insert(const vector<string> &name_vec, const vector<double> &value_vec)
{
	assert(name_vec.size() == value_vec.size());
	int vec_size = name_vec.size();
	items.reserve(vec_size);
	for(int i=0; i<vec_size; ++i)
	{
		insert(pair<string, double>(name_vec[i], value_vec[i]));
	}
}

void Transformable::insert(iterator first, iterator last )
{
	items.insert(first, last);
}

void Transformable::insert(const Transformable &insert_items)
{
	for(auto &ipar : insert_items)
	{
		items[ipar.first] = ipar.second;
	}
}

size_t Transformable::erase(const string &name)
{
	return items.erase(name);
}

void Transformable::erase(iterator position)
{
	items.erase(position);
}

void Transformable::erase(const Parameters &erase_pars)
{
	auto end_erase_pars = erase_pars.end();
	auto end_items = items.end();
	auto item_iter = items.end();;
	for(auto iter = erase_pars.begin(); iter != end_erase_pars;) {
		item_iter = items.find((*iter).first);
		if (item_iter != end_items)
		{
		   items.erase(item_iter++);
		}
		else
		{
			++iter;
	   }
	}
}

void Transformable::erase(const vector<string> &erase_par_names)
{
	auto end_erase_pars = erase_par_names.end();
	auto end_items = items.end();
	auto item_iter = items.end();;
	for(auto iter = erase_par_names.begin(); iter != end_erase_pars;) {
		item_iter = items.find((*iter));
		if (item_iter != end_items)
		{
		   items.erase(item_iter++);
		}
		else
		{
			++iter;
	   }
	}
}

Transformable::iterator Transformable::find(const string &name)
{
	return items.find(name);
}

Transformable::const_iterator Transformable::find(const string &name) const
{
	return items.find(name);
}

const double* Transformable::get_rec_ptr(const string &name) const
{
	const double *ret_val = 0;
	Transformable::const_iterator iter = find(name);
	if (iter != end()) {
		ret_val = &((*iter).second);
	}
	return ret_val;
}

double Transformable::get_rec(const string &name) const
{
	double ret_val = no_data;
	Transformable::const_iterator iter = find(name);
	if (iter != end()) {
		ret_val = (*iter).second;
	}
	else {
		throw(Transformable_value_error(name));
	}
	return ret_val;
}

void Transformable::update_rec(const string &name, double value)
{
	Transformable::iterator iter = find(name);
	assert(iter != end());
	if (iter != end()) {
		(*iter).second = value;
	}
	else {
		throw(Transformable_value_error(name));
	}
}

void Transformable::update(const vector<string> &names, const vector<double> &values)
{
	items.clear();
	assert(names.size() == values.size());
	size_t n_rec = names.size();
	items.reserve(n_rec);
	for (size_t i=0; i<n_rec; ++i)
	{
		items.insert(make_pair(names[i], values[i]));
	}
}
void Transformable::update_without_clear(const vector<string> &names, const vector<double> &values)
{
	assert(names.size() == values.size());
	size_t n_rec = names.size();
	items.reserve(n_rec);
	for (size_t i=0; i<n_rec; ++i)
	{
		items[names[i]] = values[i];
	}
}

void Transformable::update_without_clear(const vector<string> &names, const Eigen::VectorXd &values)
{
	assert(names.size() == values.size());
	size_t n_rec = names.size();
	items.reserve(n_rec);
	for (size_t i = 0; i<n_rec; ++i)
	{
		items[names[i]] = values[i];
	}
}

vector<double> Transformable::get_data_vec(const vector<string> &keys) const
{
	vector<double> v;
	v.resize(keys.size(), 0.0);

	double value;
	int i = 0;

	for(vector<string>::const_iterator b=keys.begin(), e=keys.end();
		b!=e; ++b, ++i)
	{
		value = get_rec((*b));
		v[i] = value;
	}
	return v;
}


Eigen::VectorXd Transformable::get_data_eigen_vec(const vector<string> &keys) const
{
	VectorXd vec;
	vec.resize(keys.size());
	int i = 0;
	for (auto &k : keys)
	{
		vec(i++) = get_rec(k);
	}
	return vec;
}

Eigen::VectorXd Transformable::get_partial_data_eigen_vec(const vector<string> &keys) const
{
	VectorXd vec;
	vec.resize(keys.size());
	int i = 0;
	auto end = items.end();
	for (auto &k : keys)
	{
		auto iter = items.find(k);
		if (iter != end)
		{
			vec(i) = get_rec(k);
		}
		else
		{
			vec(i) = 0;
		}
		++i;
	}
	return vec;
}





vector<string> Transformable::get_keys() const
{
	vector<string> ret_val;
	ret_val.reserve(size());
	for(Transformable::const_iterator b=begin(), e=end();
		b!=e; ++b)
	{
		ret_val.push_back((*b).first);
	}
	return ret_val;
}

double Transformable::l2_norm() const
{
   double norm=0;
   for (auto &i : items)
   {
	norm += pow(i.second, 2.0);
   }
   norm = sqrt(norm);
   return norm;
}


double  Transformable::l2_norm(const Transformable &d1, const Transformable &d2)
{
	double norm = 0.0;
	const auto it_end = d1.end();
	assert(d1.size() == d2.size());
	for (const auto &i_d2 : d2.items)
	{
		auto it_d1 = d1.find(i_d2.first);
		if (it_d1 != it_end)
		{
			norm += pow(it_d1->second - i_d2.second, 2.0);
		}
		else
		{
			assert(true);
		}
	}
	norm = sqrt(norm);
	return norm;
}

ostream& operator<<(ostream& out, const Transformable &rhs)
{
	for(Transformable::const_iterator b=rhs.begin(), e=rhs.end();
		b!=e; ++b)
	{
		out << b->first << " = " << b->second << endl;
	}
	return out;
}




void Parameters::read_par_file(ifstream &fin,  map<string, double> &offset, map<string, double> &scale)
{
	clear();
	offset.clear();
	scale.clear();

	vector<string> tokens;
	string line;
	getline(fin, line);
	while (getline(fin, line))
	{
		strip_ip(line);
		if (line.length() > 0)
		{
			tokens.clear();
			tokenize(line, tokens);
			string name = upper_cp(tokens[0]);
			double par_val = convert_cp<double>(tokens[1]);
			double offset_val = convert_cp<double>(tokens[2]);
			double scale_val = convert_cp<double>(tokens[3]);
			insert(name, par_val);
			offset[name] = offset_val;
			scale[name] = scale_val;
		}
	}
}

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
 * @file Transformable.cpp
 * @brief Implementation of Transformable.
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

/**
 * @brief Transformable.
 *
 * @param copyin Description.
 *
 * @return Description.
 */
Transformable::Transformable(const Transformable &copyin) : items(copyin.items)
{
}

/**
 * @brief Transformable.
 *
 * @param copyin Description.
 *
 * @return Description.
 */
Transformable::Transformable(const Transformable &&copyin)
{
	items = std::move(copyin.items);
}

/**
 * @brief Transformable.
 *
 * @param copyin Description.
 * @param copy_names Description.
 *
 * @return Description.
 */
Transformable::Transformable(const Transformable &copyin, const vector<string> &copy_names)
{
	for (vector<string>::const_iterator b=copy_names.begin(), e=copy_names.end(); b!=e; ++b) {
		items[*b] = copyin.get_rec(*b);
	}
}

/**
 * @brief Transformable.
 *
 * @param names Description.
 * @param values Description.
 *
 * @return Description.
 */
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

/**
 * @brief Get notnormal keys.
 *
 * @return Description.
 */
vector<string> Transformable::get_notnormal_keys()
{
	vector<string> not_normal;
	for (auto &i : items)
	{
		if ((isnormal(i.second)) || (i.second ==0.0))
			continue;
		not_normal.push_back(i.first);
	}
	return not_normal;
}

/**
 * @brief Insert.
 *
 * @param name Description.
 * @param value Description.
 *
 * @return Description.
 */
pair<Transformable::iterator,bool> Transformable::insert(const string &name, double value)
{
	pair<string, double> rec(name, value);
	return items.insert(rec);
}

/**
 * @brief Insert.
 *
 * @param x Description.
 *
 * @return Description.
 */
pair<Transformable::iterator,bool>  Transformable::insert(const pair<string, double> &x)
{
	return items.insert(x);
}

/**
 * @brief Insert.
 *
 * @param name_vec Description.
 * @param value_vec Description.
 */
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

/**
 * @brief Insert.
 *
 * @param first Description.
 * @param last Description.
 */
void Transformable::insert(iterator first, iterator last )
{
	items.insert(first, last);
}

/**
 * @brief Insert.
 *
 * @param insert_items Description.
 */
void Transformable::insert(const Transformable &insert_items)
{
	for(auto &ipar : insert_items)
	{
		items[ipar.first] = ipar.second;
	}
}

/**
 * @brief Erase.
 *
 * @param name Description.
 *
 * @return Description.
 */
size_t Transformable::erase(const string &name)
{
	return items.erase(name);
}

/**
 * @brief Erase.
 *
 * @param position Description.
 */
void Transformable::erase(iterator position)
{
	items.erase(position);
}

/**
 * @brief Erase.
 *
 * @param erase_pars Description.
 */
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

/**
 * @brief Erase.
 *
 * @param erase_par_names Description.
 */
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

/**
 * @brief Find.
 *
 * @param name Description.
 *
 * @return Description.
 */
Transformable::iterator Transformable::find(const string &name)
{
	return items.find(name);
}

/**
 * @brief Find.
 *
 * @param name Description.
 *
 * @return Description.
 */
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

/**
 * @brief Get rec.
 *
 * @param name Description.
 *
 * @return Description.
 */
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

/**
 * @brief Update rec.
 *
 * @param name Description.
 * @param value Description.
 */
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

/**
 * @brief Update.
 *
 * @param names Description.
 * @param values Description.
 */
void Transformable::update(const vector<string>& names, const Eigen::VectorXd& values)
{
	items.clear();
	assert(names.size() == values.size());
	size_t n_rec = names.size();
	items.reserve(n_rec);
	for (size_t i = 0; i < n_rec; ++i)
	{
		items[names[i]] = values[i];
	}
}

/**
 * @brief Update.
 *
 * @param names Description.
 * @param values Description.
 */
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
/**
 * @brief Update without clear.
 *
 * @param names Description.
 * @param values Description.
 */
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

/**
 * @brief Update without clear.
 *
 * @param names Description.
 * @param values Description.
 */
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

/**
 * @brief Get data vec.
 *
 * @param keys Description.
 *
 * @return Description.
 */
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


/**
 * @brief Get data eigen vec.
 *
 * @param keys Description.
 *
 * @return Description.
 */
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

/**
 * @brief Get partial data eigen vec.
 *
 * @param keys Description.
 *
 * @return Description.
 */
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





/**
 * @brief Get keys.
 *
 * @return Description.
 */
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

/**
 * @brief L2 norm.
 *
 * @return Description.
 */
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


/**
 * @brief L2 norm.
 *
 * @param d1 Description.
 * @param d2 Description.
 *
 * @return Description.
 */
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

/**
 * @brief Overloaded operator << operator.
 *
 * @param out Description.
 * @param rhs Description.
 *
 * @return Description.
 */
ostream& operator<<(ostream& out, const Transformable &rhs)
{
	for(Transformable::const_iterator b=rhs.begin(), e=rhs.end();
		b!=e; ++b)
	{
		out << b->first << " = " << b->second << endl;
	}
	return out;
}




/**
 * @brief Read par file.
 *
 * @param fin Description.
 * @param offset Description.
 * @param scale Description.
 */
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

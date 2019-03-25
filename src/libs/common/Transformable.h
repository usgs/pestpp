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

#ifndef TRANSFORMABLE_H_
#define TRANSFORMABLE_H_

#include <unordered_map>
#include <string>
#include <vector>
#include <ostream>
#include <utility>
#include <Eigen/Dense>
#include <map>
#include "pest_error.h"

using namespace std;

class ParameterInfo;
class Parameters;
class Observations;

class Transformable_value_error : public PestError {
public:
	Transformable_value_error(const string &_str, const string &_message="")
		: PestError(_message) , str(_str){
		message = string("Transformable_value_error:  Error can not access tag: \"") + _str + "\"" + message;
	}
	virtual ~Transformable_value_error() throw () {};
	virtual const char* what() const throw()
	{
		return message.c_str();
	}
private:
	string str;
};



class Transformable {
public:
	static const double no_data;
	typedef unordered_map<string, double>::iterator iterator;
	typedef unordered_map<string, double>::const_iterator const_iterator;
	Transformable(){};
	Transformable(const Transformable &copyin);
	Transformable(const Transformable &&copyin);
	Transformable(const Transformable &copyin, const vector<string> &copy_names);
	Transformable(const std::vector<std::string> &names, const Eigen::VectorXd &values);
	const Transformable& operator=(const Transformable &rhs);
	bool operator==(const Transformable &rhs) const;
	bool operator!=(const Transformable &rhs) const;
	Transformable& operator+=(const Transformable &rhs);
	Transformable& operator-=(const Transformable &rhs);
	Transformable& operator*=(double scale);
	Transformable operator*(double scale) const;
	Transformable operator-(const Transformable &rhs) const;
	Transformable operator+(const Transformable &rhs) const;
	double &operator[](const string &name);
	pair<iterator,bool> insert(const string &name, double value);
	pair<iterator, bool> insert(const pair<string, double> &x);
	void insert(const vector<string> &name_vec, const vector<double> &value_vec);
	void insert (iterator first, iterator last );
	void insert(const Transformable &insert_pars);
	size_t erase(const string &name);
	void erase(iterator position);
	void erase(const Parameters &erase_pars);
	void erase(const vector<string> &erase_par_names);
	iterator find(const string &name);
	const double* get_rec_ptr(const string &name) const;
	double get_rec(const string &name) const;
	template <class NameIterator>
	Transformable get_subset (NameIterator first, NameIterator last) const;
	void update_rec(const string &name, double value);
	void update(const vector<string> &names, const vector<double> &values);
	void update_without_clear(const vector<string> &names, const vector<double> &values);
	void update_without_clear(const vector<string> &names, const Eigen::VectorXd &values);
	const_iterator find(const string &name) const;
	size_t size() const {return items.size();}
	void clear() {items.clear();}
	vector<string> get_notnormal_keys();
	vector<string> get_keys() const;
	vector<double> get_data_vec(const vector<string> &keys) const;
	Eigen::VectorXd get_data_eigen_vec(const vector<string> &keys) const;
	Eigen::VectorXd get_partial_data_eigen_vec(const vector<string> &keys) const;
	Transformable::iterator begin(){return items.begin();}
	Transformable::const_iterator begin() const {return items.begin();}
	Transformable::iterator end() {return items.end();}
	Transformable::const_iterator end() const {return items.end();}
	double l2_norm() const;
	static double l2_norm(const Transformable &d1, const Transformable &d2);
	virtual ~Transformable(){}
protected:
	unordered_map<string, double> items;
};

ostream& operator<< (ostream& out, const Transformable &rhs);


class Parameters : public Transformable {
public:
	Parameters() : Transformable(){}
	Parameters(const Transformable &copyin) : Transformable(copyin) {}
	Parameters(const Parameters &copyin) : Transformable(copyin) {}
	Parameters(const Parameters &copyin, const vector<string> &copy_names) : Transformable(copyin, copy_names){}
	Parameters(const std::vector<std::string> &names, const Eigen::VectorXd &values) : Transformable(names, values) {}
	template <class NameIterator>
	Parameters get_subset (NameIterator first, NameIterator last)const;
	void read_par_file(std::ifstream &fin, std::map<std::string, double> &offset, std::map<std::string, double> &scale);
	virtual ~Parameters(){}
private:
};


class Observations : public Transformable {
public:
	Observations() : Transformable(){}
	Observations(const Observations &copyin) : Transformable(copyin) {}
	Observations(const Observations &copyin, const vector<string> &copy_names) : Transformable(copyin, copy_names){}
	Observations(const vector<string> &copyin, const Eigen::VectorXd &copy_vec) : Transformable(copyin, copy_vec) {}
	template <class NameIterator>
	Observations get_subset(NameIterator first, NameIterator last)const;
	virtual ~Observations(){}
private:
};


template <class NameIterator>
Transformable Transformable::get_subset (const NameIterator first, const NameIterator last) const
{
	Transformable subset;
	for(auto i = first; i!=last; ++i)
	{
		auto t_iter = find(*i);
		if (t_iter == this->end())
		{
			throw(Transformable_value_error(*i));
		}
		subset.insert(t_iter);
	}
}


template <class NameIterator>
Parameters Parameters::get_subset (NameIterator first, NameIterator last) const
{
	Parameters subset;
	for(auto i = first; i!=last; ++i)
	{
		auto t_iter = find(*i);
		if (t_iter == this->end())
		{
			throw(Transformable_value_error(*i));
		}
		subset.insert(t_iter->first, t_iter->second);
	}
	return subset;
}

template <class NameIterator>
Observations Observations::get_subset(NameIterator first, NameIterator last) const
{
	Observations subset;
	for (auto i = first; i != last; ++i)
	{
		auto t_iter = find(*i);
		if (t_iter == this->end())
		{
			throw(Transformable_value_error(*i));
		}
		subset.insert(t_iter->first, t_iter->second);
	}
	return subset;
}


#endif /* TRANSFORMABLE_H_ */

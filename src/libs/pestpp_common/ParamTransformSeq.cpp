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

#include <iostream>
#include <cmath>
#include <vector>
#include <iterator>
#include <algorithm>
#include "ParamTransformSeq.h"
#include "Transformation.h"
#include "Transformable.h"
#include "Jacobian.h"

using namespace std;

//These routines handle reference count on the Transformation used by ParamTransformSeq.
// Using reference counting makes it possible to not store new copies of all the
// transformations each time a ParamTransformSeq is copied.

map<const Transformation*, int> ParamTransformSeq::tran_ref_count = map<const Transformation*, int>();

int ParamTransformSeq::tran_add_ref_count(const Transformation *new_tran)
{
	int ref_count = 0;
	map<const Transformation *,int>::iterator it;
	it = tran_ref_count.find(new_tran);
	if (it != tran_ref_count.end()) {
		(*it).second += 1;
		ref_count = (*it).second;
	}
	else {
		tran_ref_count.insert(pair<const Transformation *, int>(new_tran, 1));
		ref_count = 1;
	}
	return ref_count;
}

int ParamTransformSeq::tran_sub_ref_count(const Transformation *new_tran)
{
	int ref_count = 0;
	const Transformation *t_ptr;
	map<const Transformation *,int>::iterator it;
	it = tran_ref_count.find(new_tran);
	if (it != tran_ref_count.end()) {
		(*it).second -= 1;
		ref_count = (*it).second;
		if (ref_count<=0) {
			t_ptr = (*it).first;
			tran_ref_count.erase(it);
			delete t_ptr;
		}
	}
	return ref_count;
}


ParamTransformSeq::ParamTransformSeq(const ParamTransformSeq &rhs)
{
	copy(rhs);
}

ParamTransformSeq::ParamTransformSeq(const ParamTransformSeq &rhs, const set<Transformation *> &deep_copy_tran_set)
{
	copy(rhs, deep_copy_tran_set);
}

ParamTransformSeq::~ParamTransformSeq()
{
	clear();
}

void ParamTransformSeq::clear()
{
	clear_tranSeq_ctl2model();
	clear_tranSeq_ctl2active_ctl();
	clear_tranSeq_active_ctl2numeric();
	default_deep_copy_tran_set.clear();
	name = "empty";
}

void ParamTransformSeq::clear_tranSeq_ctl2model()
{
	for(vector<Transformation*>::iterator i = tranSeq_ctl2model.begin(),
		e=tranSeq_ctl2model.end(); i != e; ++i)
	{
		tran_sub_ref_count(*i);
		default_deep_copy_tran_set.erase(*i);
	}
	tranSeq_ctl2model.clear();
}

void  ParamTransformSeq::clear_tranSeq_ctl2active_ctl()
{
	for(vector<Transformation*>::iterator i = tranSeq_ctl2active_ctl.begin(),
		e=tranSeq_ctl2active_ctl.end(); i != e; ++i)
	{
		tran_sub_ref_count(*i);
		default_deep_copy_tran_set.erase(*i);
	}
	tranSeq_ctl2active_ctl.clear();
}


void ParamTransformSeq::clear_tranSeq_active_ctl2numeric()
{
	for(vector<Transformation*>::iterator i = tranSeq_active_ctl2numeric.begin(),
		e=tranSeq_active_ctl2numeric.end(); i != e; ++i)
	{
		tran_sub_ref_count(*i);
		default_deep_copy_tran_set.erase(*i);
	}
	tranSeq_active_ctl2numeric.clear();
}

void ParamTransformSeq::deep_copy(const ParamTransformSeq &rhs)
{
	set<Transformation *> deep_copy_set(rhs.tranSeq_ctl2model.begin(), rhs.tranSeq_ctl2model.end());
	deep_copy_set.insert(rhs.tranSeq_ctl2active_ctl.begin(), rhs.tranSeq_ctl2active_ctl.end());
	deep_copy_set.insert(rhs.tranSeq_active_ctl2numeric.begin(), rhs.tranSeq_active_ctl2numeric.end());
	copy(rhs, deep_copy_set);
}

ParamTransformSeq& ParamTransformSeq::operator=(const ParamTransformSeq &rhs)
{
	copy(rhs);
	return *this;
}

void ParamTransformSeq::copy(const ParamTransformSeq &rhs)
{
	copy(rhs, rhs.default_deep_copy_tran_set);
}

void ParamTransformSeq::copy(const ParamTransformSeq &rhs, const set<Transformation *> &deep_copy_tran_set)
{
	clear();
	name = "copy of " + rhs.name;

	Transformation *t_ptr;
	for(vector<Transformation*>::const_iterator i = rhs.tranSeq_ctl2model.begin(),
		e=rhs.tranSeq_ctl2model.end(); i != e; ++i)
	{
		if (deep_copy_tran_set.find(*i) ==  deep_copy_tran_set.end()) {
			t_ptr = *i;
		}
		else {
			t_ptr = (*i)->clone();
		}
		push_back_ctl2model(t_ptr);
		if (rhs.default_deep_copy_tran_set.find(*i) !=  rhs.default_deep_copy_tran_set.end()) {
			default_deep_copy_tran_set.insert(t_ptr);
		}
	}
	for(vector<Transformation*>::const_iterator i = rhs.tranSeq_ctl2active_ctl.begin(),
		e=rhs.tranSeq_ctl2active_ctl.end(); i != e; ++i)
	{
		if (deep_copy_tran_set.find(*i) ==  deep_copy_tran_set.end()) {
			t_ptr = *i;
		}
		else {
			t_ptr = (*i)->clone();
		}
		push_back_ctl2active_ctl(t_ptr);
		if (rhs.default_deep_copy_tran_set.find(*i) !=  rhs.default_deep_copy_tran_set.end()) {
			default_deep_copy_tran_set.insert(t_ptr);
		}
	}
	for(vector<Transformation*>::const_iterator i = rhs.tranSeq_active_ctl2numeric.begin(),
		e=rhs.tranSeq_active_ctl2numeric.end(); i != e; ++i)
	{
		if (deep_copy_tran_set.find(*i) ==  deep_copy_tran_set.end()) {
			t_ptr = *i;
		}
		else {
			t_ptr = (*i)->clone();
		}
		push_back_active_ctl2numeric(t_ptr);
		if (rhs.default_deep_copy_tran_set.find(*i) !=  rhs.default_deep_copy_tran_set.end()) {
			default_deep_copy_tran_set.insert(t_ptr);
		}
	}
	for (auto &itran_seq : rhs.custom_tran_seq )
	{
	  for (auto &itran : itran_seq.second )
	  {
		if (deep_copy_tran_set.find(itran) ==  deep_copy_tran_set.end()) {
			t_ptr = itran;
		}
		else {
			t_ptr = itran->clone();
		}
		custom_tran_seq[itran_seq.first].push_back(t_ptr);
		if (rhs.default_deep_copy_tran_set.find(itran) !=  rhs.default_deep_copy_tran_set.end()) {
			default_deep_copy_tran_set.insert(t_ptr);
		}
	  }
	}
}

ParamTransformSeq &ParamTransformSeq::operator+=(const ParamTransformSeq &rhs)
{
	append(rhs, default_deep_copy_tran_set);
	return *this;
}

void ParamTransformSeq::append(const ParamTransformSeq &rhs)
{
	append(rhs, rhs.default_deep_copy_tran_set);
}

void ParamTransformSeq::append(const ParamTransformSeq &rhs, const set<Transformation *> &deep_copy_tran_set )
{
	Transformation *t_ptr;

	for (vector<Transformation*>::const_iterator i=rhs.tranSeq_ctl2model.begin(), e=rhs.tranSeq_ctl2model.end();
		i!=e; ++i)
	{
		if (deep_copy_tran_set.find(*i) ==  deep_copy_tran_set.end()) {
			t_ptr = *i;
		}
		else {
			t_ptr = (*i)->clone();
		}
		push_back_ctl2model(t_ptr);
		if (rhs.default_deep_copy_tran_set.find(*i) !=  rhs.default_deep_copy_tran_set.end()) {
			default_deep_copy_tran_set.insert(t_ptr);
		}
	}
	for (vector<Transformation*>::const_iterator i=rhs.tranSeq_ctl2active_ctl.begin(), e=rhs.tranSeq_ctl2active_ctl.end();
		i!=e; ++i)
	{
		if (deep_copy_tran_set.find(*i) ==  deep_copy_tran_set.end()) {
			t_ptr = *i;
		}
		else {
			t_ptr = (*i)->clone();
		}
		push_back_ctl2active_ctl(t_ptr);
		if (rhs.default_deep_copy_tran_set.find(*i) !=  rhs.default_deep_copy_tran_set.end()) {
			default_deep_copy_tran_set.insert(t_ptr);
		}
	}
	for (vector<Transformation*>::const_iterator i=rhs.tranSeq_active_ctl2numeric.begin(), e=rhs.tranSeq_active_ctl2numeric.end();
		i!=e; ++i)
	{
		if (deep_copy_tran_set.find(*i) ==  deep_copy_tran_set.end()) {
			t_ptr = *i;
		}
		else {
			t_ptr = (*i)->clone();
		}
		push_back_active_ctl2numeric(t_ptr);
		if (rhs.default_deep_copy_tran_set.find(*i) !=  rhs.default_deep_copy_tran_set.end()) {
			default_deep_copy_tran_set.insert(t_ptr);
		}
	}
	for (auto &itran_seq : rhs.custom_tran_seq )
	{
	  for (auto &itran : itran_seq.second )
	  {
		if (deep_copy_tran_set.find(itran) ==  deep_copy_tran_set.end()) {
			t_ptr = itran;
		}
		else {
			t_ptr = itran->clone();
		}
		custom_tran_seq[itran_seq.first].push_back(t_ptr);
		if (rhs.default_deep_copy_tran_set.find(itran) !=  rhs.default_deep_copy_tran_set.end()) {
			default_deep_copy_tran_set.insert(t_ptr);
		}
	  }
	}
}

void ParamTransformSeq::push_back_ctl2model(Transformation *tr)
{
	tranSeq_ctl2model.push_back(tr);
	tran_add_ref_count(tr);
}

void ParamTransformSeq::push_back_ctl2active_ctl(Transformation *tr)
{
	tranSeq_ctl2active_ctl.push_back(tr);
	tran_add_ref_count(tr);
}

void ParamTransformSeq::insert_ctl2active_ctl(const string &location_name, Transformation *tr)
{
	const auto iter = find_in_ctl2active_ctl(location_name);
	if (iter != tranSeq_ctl2active_ctl.end())
	{
		tranSeq_ctl2active_ctl.insert(iter, tr);
	}
	else
	{
		stringstream err_str;
		err_str << "ParamTransformSeq::insert_ctl2active_ctl: location name = \"" << location_name << "\" not found";
		throw(PestIndexError(err_str.str()));
	}
}

void ParamTransformSeq::push_back_active_ctl2numeric(Transformation *tr)
{
	tranSeq_active_ctl2numeric.push_back(tr);
	tran_add_ref_count(tr);
}

void ParamTransformSeq::push_front_ctl2active_ctl(Transformation *tr)
{
	tranSeq_ctl2active_ctl.insert(tranSeq_ctl2active_ctl.begin(), tr);
}

void ParamTransformSeq::ctl2model_ip(Parameters &data) const
{
	vector<Transformation*>::const_iterator iter, e;

	for(iter = tranSeq_ctl2model.begin(), e = tranSeq_ctl2model.end();
		iter != e; ++iter)
	{
		(*iter)->forward(data);
	}
}

Parameters ParamTransformSeq::ctl2model_cp(const Parameters &data) const
{
	Parameters ret_val(data);
	ctl2model_ip(ret_val);
	return ret_val;
}

void ParamTransformSeq::ctl2active_ctl_ip(Parameters &data) const
{
	vector<Transformation*>::const_iterator iter, e;
	for(iter = tranSeq_ctl2active_ctl.begin(), e = tranSeq_ctl2active_ctl.end();
		iter != e; ++iter)
	{
		(*iter)->forward(data);

	}
}

Parameters ParamTransformSeq::ctl2active_ctl_cp(const Parameters &data) const
{
	Parameters ret_val(data);
	ctl2active_ctl_ip(ret_val);
	return ret_val;
}

void ParamTransformSeq::ctl2numeric_ip(Parameters &data) const
{
	 ctl2active_ctl_ip(data);
	 active_ctl2numeric_ip(data);
}

Parameters ParamTransformSeq::ctl2numeric_cp(const Parameters &data) const
{
	Parameters ret_val(data);
	ctl2numeric_ip(ret_val);
	return ret_val;
}

void ParamTransformSeq::model2ctl_ip(Parameters &data) const
{
	vector<Transformation*>::const_reverse_iterator iter, e;

	for(iter = tranSeq_ctl2model.rbegin(), e = tranSeq_ctl2model.rend();
		iter != e; ++iter)
	{
		(*iter)->reverse(data);
	}
}

Parameters ParamTransformSeq::model2ctl_cp(const Parameters &data) const
{
	Parameters ret_val(data);
	model2ctl_ip(ret_val);
	return ret_val;
}

void ParamTransformSeq::model2active_ctl_ip(Parameters &data) const
{
	model2ctl_ip(data);
	ctl2active_ctl_ip(data);
}

Parameters ParamTransformSeq::model2active_ctl_cp(const Parameters &data) const
{
	Parameters ret_val(data);
	model2active_ctl_ip(ret_val);
	return(data);
}


void ParamTransformSeq::numeric2active_ctl_ip(Parameters &data) const
{
	vector<Transformation*>::const_reverse_iterator iter, e;

	for(iter = tranSeq_active_ctl2numeric.rbegin(), e = tranSeq_active_ctl2numeric.rend();
		iter != e; ++iter)
	{
		(*iter)->reverse(data);
	}
}


Parameters ParamTransformSeq::numeric2active_ctl_cp(const Parameters &data) const
{
	Parameters ret_val(data);
	numeric2active_ctl_ip(ret_val);
	return ret_val;
}

void ParamTransformSeq::numeric2ctl_ip(Parameters &data) const
{
	numeric2active_ctl_ip(data);
	active_ctl2ctl_ip(data);
}

Parameters ParamTransformSeq::numeric2ctl_cp(const Parameters &data) const
{
	Parameters ret_val(data);
	numeric2ctl_ip(ret_val);
	return ret_val;
}


void ParamTransformSeq::numeric2model_ip(Parameters &data) const
{
	numeric2ctl_ip(data);
	ctl2model_ip(data);
}

Parameters ParamTransformSeq::numeric2model_cp(const Parameters &data) const
{
	Parameters ret_val(data);
	numeric2model_ip(ret_val);
	return ret_val;
}

void ParamTransformSeq::model2numeric_ip(Parameters &data) const
{
	model2ctl_ip(data);
	ctl2numeric_ip(data);
}

Parameters ParamTransformSeq::model2numeric_cp(const Parameters &data) const
{
	Parameters ret_val(data);
	model2numeric_ip(ret_val);
	return ret_val;
}

void ParamTransformSeq::active_ctl2numeric_ip(Parameters &data) const
{
	vector<Transformation*>::const_iterator iter, e;

	for(iter = tranSeq_active_ctl2numeric.begin(), e = tranSeq_active_ctl2numeric.end();
		iter != e; ++iter)
	{
		(*iter)->forward(data);
	}
}

Parameters ParamTransformSeq::active_ctl2numeric_cp(const Parameters &data) const
{
	Parameters ret_val(data);
	active_ctl2numeric_ip(ret_val);
	return ret_val;
}


void ParamTransformSeq::active_ctl2ctl_ip(Parameters &data) const
{
	vector<Transformation*>::const_reverse_iterator iter, e;

	for(iter = tranSeq_ctl2active_ctl.rbegin(), e = tranSeq_ctl2active_ctl.rend();
		iter != e; ++iter)
	{
		(*iter)->reverse(data);
	}
}

Parameters ParamTransformSeq::active_ctl2ctl_cp(const Parameters &data) const
{
	Parameters ret_val(data);
	active_ctl2ctl_ip(ret_val);
	return ret_val;
}

void ParamTransformSeq::active_ctl2model_ip(Parameters &data) const
{
	active_ctl2ctl_ip(data);
	ctl2model_ip(data);
}

Parameters ParamTransformSeq::active_ctl2model_cp(const Parameters &data) const
{
	Parameters ret_val(data);
	active_ctl2model_ip(ret_val);
	return ret_val;
}


void ParamTransformSeq::del_numeric_2_del_active_ctl_ip(Parameters &del_data, Parameters &data) const
{
	map<string, double> factors;
	vector<Transformation*>::const_reverse_iterator iter, e;
	for (iter = tranSeq_active_ctl2numeric.rbegin(), e = tranSeq_active_ctl2numeric.rend();
		iter != e; ++iter)
	{
		(*iter)->d2_to_d1(del_data, data);
	}
}



bool ParamTransformSeq::is_one_to_one() const
{
	for (vector<Transformation*>::const_iterator b=tranSeq_ctl2model.begin(), e=tranSeq_ctl2model.end();
		b!=e; ++b) {
			if((*b)->is_one_to_one() == false) {
				return false;
			}
	}
	for (vector<Transformation*>::const_iterator b=tranSeq_ctl2active_ctl.begin(), e=tranSeq_ctl2active_ctl.end();
		b!=e; ++b) {
			if((*b)->is_one_to_one() == false) {
				return false;
			}
	}
	for (vector<Transformation*>::const_iterator b=tranSeq_active_ctl2numeric.begin(), e=tranSeq_active_ctl2numeric.end();
		b!=e; ++b) {
			if((*b)->is_one_to_one() == false) {
				return false;
			}
	}
	return true;
}

Transformation* ParamTransformSeq::get_transformation(const string &name)
{
	Transformation *t_ptr = 0;
	vector<Transformation*>::iterator iter;

	iter = find_in_ctl2model(name);
	if (iter != tranSeq_ctl2model.end())
	{
		t_ptr = *iter;
		return t_ptr;
	}

	iter = find_in_ctl2active_ctl(name);
	if (iter != tranSeq_ctl2active_ctl.end())
	{
		t_ptr = *iter;
		return t_ptr;
	}

	iter = find_in_active_ctl2numeric(name);
	if (iter != tranSeq_active_ctl2numeric.end())
	{
		t_ptr = *iter;
		return t_ptr;
	}
	return t_ptr;
}

void ParamTransformSeq::print(ostream &os) const
{
	os << "ParamTransformSeq name = " << name << endl;
	os << "Control file to model transformations" << endl;
	for (vector<Transformation*>::const_iterator b=tranSeq_ctl2model.begin(), e=tranSeq_ctl2model.end();
		b!=e; ++b) {
			(*b)->print(os);
	}
	os << "Control file to derivative transformations" << endl;
	for (vector<Transformation*>::const_iterator b=tranSeq_ctl2active_ctl.begin(), e=tranSeq_ctl2active_ctl.end();
		b!=e; ++b) {
			(*b)->print(os);
	}
	os << "Derivative to numeric transformations" << endl;
	for (vector<Transformation*>::const_iterator b=tranSeq_active_ctl2numeric.begin(), e=tranSeq_active_ctl2numeric.end();
		b!=e; ++b) {
			(*b)->print(os);
	}
}

void ParamTransformSeq::add_custom_tran_seq(const std::string &name,  const vector<Transformation*> &tran_seq)
{
	custom_tran_seq[name] = tran_seq;
}

const vector<Transformation*> ParamTransformSeq::get_custom_tran_seq(const string &name) const
{
	auto iter = custom_tran_seq.find(name);
	assert(iter != custom_tran_seq.end());
	return iter->second;
}

void ParamTransformSeq::custom_tran_seq_forward_ip(const std::string &name, Parameters &data) const
{
	auto iter= custom_tran_seq.find(name);
	assert(iter != custom_tran_seq.end());
	const vector<Transformation*> &tran_vec = iter->second;

	for(const auto &itran : tran_vec)
	{
		itran->forward(data);
	}

}


vector<Transformation*>::iterator ParamTransformSeq::find_in_ctl2model(const string &name)
{
	auto iter = find_if(tranSeq_ctl2model.begin(), tranSeq_ctl2model.end(),
		[&name](Transformation *tr_ptr)->bool{return tr_ptr->get_name() == name;});
	return iter;
}

vector<Transformation*>::const_iterator ParamTransformSeq::find_in_ctl2model(const string &name) const
{
	vector<Transformation*>::const_iterator iter = find_if(tranSeq_ctl2model.cbegin(), tranSeq_ctl2model.cend(),
		[&name](Transformation *tr_ptr)->bool{return tr_ptr->get_name() == name;});
	return iter;
}

vector<Transformation*>::iterator ParamTransformSeq::find_in_ctl2active_ctl(const string &name)
{
	auto iter = find_if(tranSeq_ctl2active_ctl.begin(), tranSeq_ctl2active_ctl.end(),
		[&name](Transformation *tr_ptr)->bool{return tr_ptr->get_name() == name;});
	return iter;
}

vector<Transformation*>::const_iterator ParamTransformSeq::find_in_ctl2active_ctl(const string &name) const
{
	vector<Transformation*>::const_iterator iter = find_if(tranSeq_ctl2active_ctl.cbegin(), tranSeq_ctl2active_ctl.cend(),
		[&name](Transformation *tr_ptr)->bool{return tr_ptr->get_name() == name;});
	return iter;
}


vector<Transformation*>::iterator ParamTransformSeq::find_in_active_ctl2numeric(const string &name)
{
	auto iter = find_if(tranSeq_active_ctl2numeric.begin(), tranSeq_active_ctl2numeric.end(),
		[&name](Transformation *tr_ptr)->bool{return tr_ptr->get_name() == name;});
	return iter;
}

vector<Transformation*>::const_iterator ParamTransformSeq::find_in_active_ctl2numeric(const string &name) const
{
	vector<Transformation*>::const_iterator iter = find_if(tranSeq_active_ctl2numeric.cbegin(), tranSeq_active_ctl2numeric.cend(),
		[&name](Transformation *tr_ptr)->bool{return tr_ptr->get_name() == name;});
	return iter;
}

const TranOffset *ParamTransformSeq::get_offset_ptr() const
{
	const Transformation* ptr=0;
	const auto iter = find_in_ctl2model(string("PEST to model offset transformation"));
	if (iter != tranSeq_ctl2model.end())
	{
		ptr = (*iter);
	}
	return dynamic_cast<const TranOffset*>(ptr);
}

const TranScale *ParamTransformSeq::get_scale_ptr() const
{
	const Transformation* ptr=0;
	const auto iter = find_in_ctl2model(string("PEST to model scale transformation"));
	if (iter != tranSeq_ctl2model.end())
	{
		ptr = (*iter);
	}
	return dynamic_cast<const TranScale*>(ptr);
}

const TranFixed* ParamTransformSeq::get_fixed_ptr()const
{
	const Transformation* ptr=0;
	const auto iter = find_in_ctl2active_ctl(string("PEST to model fixed transformation"));
	if (iter != tranSeq_ctl2active_ctl.end())
	{
		ptr = (*iter);
	}
	return dynamic_cast<const TranFixed*>(ptr);
}

TranFixed* ParamTransformSeq::get_fixed_ptr_4_mod()
{
	Transformation* ptr = 0;
	auto iter = find_in_ctl2active_ctl(string("PEST to model fixed transformation"));
	if (iter != tranSeq_ctl2active_ctl.end())
	{
		ptr = (*iter);
	}
	return dynamic_cast<TranFixed*>(ptr);
}


TranTied* ParamTransformSeq::get_tied_ptr()const
{
	Transformation* ptr = 0;
	const auto iter = find_in_ctl2active_ctl(string("PEST to model tied transformation"));
	if (iter != tranSeq_ctl2active_ctl.end())
	{
		ptr = (*iter);
	}
	return dynamic_cast<TranTied*>(ptr);
}

const TranLog10 *ParamTransformSeq::get_log10_ptr() const
{
	const Transformation* ptr=0;
	const auto iter = find_in_active_ctl2numeric(string("PEST to model log transformation"));
	if (iter != tranSeq_active_ctl2numeric.end())
	{
		ptr = (*iter);
	}
	return dynamic_cast<const TranLog10*>(ptr);
}


TranSVD *ParamTransformSeq::get_svda_ptr()const
{
	Transformation* ptr=0;
	auto iter = find_in_active_ctl2numeric(string("SVD Super Parameter Transformation"));
	if (iter != tranSeq_active_ctl2numeric.end())
	{
		ptr = (*iter);
	}
	return dynamic_cast<TranSVD*>(ptr);

}

TranFixed *ParamTransformSeq::get_svda_fixed_ptr()const
{
	Transformation* ptr=0;
	const auto iter = find_in_ctl2active_ctl(string("SVDA Fixed Parameter Transformation"));
	if (iter != tranSeq_ctl2active_ctl.end())
	{
		ptr = (*iter);
	}
	return dynamic_cast<TranFixed*>(ptr);
}



void ParamTransformSeq::jac_test_ip(Jacobian &jac) const
{
	vector<Transformation*>::const_reverse_iterator iter, e;
	vector<Transformation*> test_seq;;
	test_seq.push_back(tranSeq_active_ctl2numeric[1]);

	jac.print(cout);
	for (iter = test_seq.rbegin(), e = test_seq.rend();
		iter != e; ++iter)
	{
		(*iter)->jacobian_reverse(jac);
		cout << (*iter)->get_name() << endl;
		jac.print(cout);
	}

	vector<Transformation*>::const_iterator iter2, e2;
	for (iter2 = test_seq.begin(), e2 = test_seq.end();
		iter2 != e2; ++iter2)
	{
		(*iter2)->jacobian_forward(jac);
		cout << (*iter2)->get_name() << endl;
		jac.print(cout);
	}
}

void ParamTransformSeq::jac_numeric2active_ctl_ip(Jacobian &jac) const
{
	vector<Transformation*>::const_reverse_iterator iter, e;

	for (iter = tranSeq_active_ctl2numeric.rbegin(), e = tranSeq_active_ctl2numeric.rend();
		iter != e; ++iter)
	{
		(*iter)->jacobian_reverse(jac);
	}
}

void ParamTransformSeq::jac_active_ctl_ip2numeric_ip(Jacobian &jac) const
{
	vector<Transformation*>::const_iterator iter, e;

	for (iter = tranSeq_active_ctl2numeric.begin(), e = tranSeq_active_ctl2numeric.end();
		iter != e; ++iter)
	{
		(*iter)->jacobian_forward(jac);
	}
}

ostream& operator<< (ostream &os, const ParamTransformSeq& val)
{
	val.print(os);
	return os;
}


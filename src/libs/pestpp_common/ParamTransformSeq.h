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
#ifndef ParamTransformSeq_H_
#define ParamTransformSeq_H_

/** @file
 @brief ParamTransformSeq class

 This file defines the ParamTransformSeq.  This class is used to define the sequence of
 transformations used to switch between model, control and numeric parameter representations
*/

#include <vector>
#include <string>
#include "Transformation.h"

class Parameters;
class JacobianAnalytic;

using namespace std;

/**
 @brief ParamTransformSeq class

 This class is used to define the sequence of transformations used to switch
 between model, control and numeric parameter representations
*/
class ParamTransformSeq {
public:
	ParamTransformSeq(const string &_name="unnamed ParamTransformSeq") : name(_name) {}
	ParamTransformSeq(const ParamTransformSeq &rhs);
	ParamTransformSeq(const ParamTransformSeq &rhs, const set<Transformation *> &deep_copy_tran_set);
	void copy(const ParamTransformSeq &rhs);
	ParamTransformSeq &operator=(const ParamTransformSeq &rhs);
	void copy(const ParamTransformSeq &rhs, const set<Transformation *> &deep_copy_tran_set);
	void deep_copy(const ParamTransformSeq &rhs);
	ParamTransformSeq &operator+=(const ParamTransformSeq &rhs);
	void append(const ParamTransformSeq &rhs);
	void append(const ParamTransformSeq &rhs, const set<Transformation *> &deep_copy_tran_set);
	virtual ~ParamTransformSeq();
	void clear();
	void clear_tranSeq_ctl2model();
	void clear_tranSeq_ctl2active_ctl();
	void clear_tranSeq_active_ctl2numeric();
	void push_back_ctl2model(Transformation *tr);
	void push_back_ctl2active_ctl(Transformation *tr);
	void insert_ctl2active_ctl(const string &location_name, Transformation *tr);
	void push_back_active_ctl2numeric(Transformation *tr);
	void push_front_ctl2active_ctl(Transformation *tr);
	void numeric2active_ctl_ip(Parameters &data) const;
	void numeric2ctl_ip(Parameters &data) const;
	void numeric2model_ip(Parameters &data) const;
	void ctl2active_ctl_ip(Parameters &data) const;
	void ctl2numeric_ip(Parameters &data) const;
	void ctl2model_ip(Parameters &data) const;
	void model2ctl_ip(Parameters &data) const;
	void model2active_ctl_ip(Parameters &data) const;
	void model2numeric_ip(Parameters &data) const;
	void active_ctl2numeric_ip(Parameters &data) const;
	void active_ctl2ctl_ip(Parameters &data) const;
	void active_ctl2model_ip(Parameters &data) const;
	Parameters numeric2active_ctl_cp(const Parameters &data) const;
	Parameters numeric2ctl_cp(const Parameters &data) const;
	Parameters numeric2model_cp(const Parameters &data) const;
	Parameters ctl2active_ctl_cp(const Parameters &data) const;
	Parameters ctl2numeric_cp(const Parameters &data) const;
	Parameters ctl2model_cp(const Parameters &data) const;
	Parameters model2ctl_cp(const Parameters &data) const;
	Parameters model2active_ctl_cp(const Parameters &data) const;
	Parameters model2numeric_cp(const Parameters &data) const;
	Parameters active_ctl2numeric_cp(const Parameters &data) const;
	Parameters active_ctl2ctl_cp(const Parameters &data) const;
	Parameters active_ctl2model_cp(const Parameters &data) const;
	void del_numeric_2_del_active_ctl_ip(Parameters &del_data, Parameters &data) const;
	Transformation* get_transformation(const string &name);
	const TranOffset *get_offset_ptr() const;
	const TranScale *get_scale_ptr() const;
	const TranFixed *get_fixed_ptr()const;
	const TranLog10 *get_log10_ptr() const;
	TranSVD *get_svda_ptr()const;
	TranFixed *get_svda_fixed_ptr()const;
	const vector<Transformation*> get_ctl2model_tranformations() const {return tranSeq_ctl2model;}
	const vector<Transformation*> get_ctl2active_ctl_tranformations() const {return tranSeq_ctl2active_ctl;}
	const vector<Transformation*> get_active_ctl2numeric_tranformations() const {return tranSeq_active_ctl2numeric;}
	const vector<Transformation*> get_custom_tran_seq(const string &name) const;
	void add_custom_tran_seq(const std::string &name,  const vector<Transformation*> &tran_seq);
	void custom_tran_seq_forward_ip(const std::string &name, Parameters &data) const;
	//void add_default_deep_copy(Transformation *tr){default_deep_copy_tran_set.insert(tr);}
	//void clear_default_deep_copies(Transformation *tr){default_deep_copy_tran_set.clear();}
	//set <Transformation *> get_default_deep_copy_vec() const {return default_deep_copy_tran_set;}
	bool is_one_to_one() const;
	void print(ostream &os) const;
	void jac_numeric2active_ctl_ip(Jacobian &jac) const;
	void jac_active_ctl_ip2numeric_ip(Jacobian &jac) const;
	void jac_test_ip(Jacobian &jac) const;
private:
	vector<Transformation*> tranSeq_ctl2model;
	vector<Transformation*> tranSeq_ctl2active_ctl;
	vector<Transformation*> tranSeq_active_ctl2numeric;
	map<string,vector<Transformation*> > custom_tran_seq;
	set <Transformation *> default_deep_copy_tran_set;
	static map<const Transformation*, int> tran_ref_count;
	static int tran_add_ref_count(const Transformation *);
	static int tran_sub_ref_count(const Transformation *);
	vector<Transformation*>::iterator find_in_ctl2model(const string &name);
	vector<Transformation*>::const_iterator find_in_ctl2model(const string &name) const;
	vector<Transformation*>::iterator find_in_ctl2active_ctl(const string &name);
	vector<Transformation*>::const_iterator find_in_ctl2active_ctl(const string &name) const;
	vector<Transformation*>::iterator find_in_active_ctl2numeric(const string &name);
	vector<Transformation*>::const_iterator find_in_active_ctl2numeric(const string &name) const;
	string name;
};

template <typename Parameters>
ostream& operator<< (ostream &os, const ParamTransformSeq& val);

#endif /* ParamTransformSeq_H_ */

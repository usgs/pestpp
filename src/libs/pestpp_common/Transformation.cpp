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
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <math.h>
#include <sstream>
#include <Eigen/Dense>
#include <cassert>
#include <iostream>
#include "Transformation.h"
#include "Transformable.h"
#include "eigen_tools.h"
#include "Jacobian.h"
#include "QSqrtMatrix.h"
#include "pest_data_structs.h"
#include "debug.h"
#include "Regularization.h"
#include "Serialization.h"
#include "eigen_tools.h"
#include "covariance.h"

using namespace std;
using namespace Eigen;

///////////////// Transformation Methods /////////////////


///////////////// TranMapBase Methods /////////////////
void TranMapBase::insert(const string &item_name, double item_value)
{
	items[item_name] = item_value;
}


void TranMapBase::insert(const Parameters &pars)
{
	for (const auto &ipar : pars)
	{
		items[ipar.first] = ipar.second;
	}
}

void TranMapBase::reset(const Parameters &pars)
{
	items.clear();
	for (const auto &ipar : pars)
	{
		items[ipar.first] = ipar.second;
	}
}

void TranMapBase::print(ostream &os) const
{
	os << "Transformation name = " << name << "; (type=TranMapBase)" << endl;
	for (map<string,double>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
			os << "  item name = " << (*b).first << ";  value = " << (*b).second << endl;
	}
}

pair<bool, double> TranMapBase::get_value(const string &name) const
{
	pair<bool, double> ret_val(false, 0.0);
	map<string, double>::const_iterator it;

	it = items.find(name);
	if (it !=items.end()) {
		ret_val = pair<bool, double>(true, (*it).second);
	}
	return ret_val;
}

///////////////// TranSetBase Methods /////////////////
void TranSetBase::insert(const string &item_name)
{
	items.insert(item_name);
}

void TranSetBase::print(ostream &os) const
{
	os << "Transformation name = " << name << "; (type=TranSetBase)" << endl;
	for (set<string>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
			os << "  item name = " << *b << endl;
	}
}

bool TranSetBase::has_value(const string &name) const
{
	bool ret_val = false;
	set<string>::const_iterator it;

	it = items.find(name);
	if (it !=items.end()) {
		ret_val = true;
	}
	return ret_val;
}

///////////////// TranOffset Methods /////////////////
void TranOffset::forward(Transformable &data)
{
	Transformable::iterator data_iter, data_end = data.end();
	for (map<string,double>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
		data_iter = data.find(b->first);
		if (data_iter != data_end)
		{
			(*data_iter).second += (*b).second;
		}
	}
}

void TranOffset::reverse(Transformable &data)
{
	Transformable::iterator data_iter, data_end = data.end();
	for (map<string,double>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
		data_iter = data.find(b->first);
		if (data_iter != data_end)
		{
			(*data_iter).second -= b->second;
		}
	}
}

void TranOffset::jacobian_forward(Jacobian &jac)
{
	Transformable &data = jac.base_numeric_parameters;
	forward(data);
}

void TranOffset::jacobian_reverse(Jacobian &jac)
{
	Transformable &data = jac.base_numeric_parameters;
	reverse(data);
}

void TranOffset::d2_to_d1(Transformable &del_data, Transformable &data)
{
	// Offset transformation does not affect derivatives.
	// Nothing to do.
	reverse(data);
}

void TranOffset::d1_to_d2(Transformable &del_data, Transformable &data)
{
	// Offset transformation does not affect derivatives.
	// Nothing to do.
	forward(data);
}

void TranOffset::print(ostream &os) const
{
	os << "Transformation name = " << name << "; (type=TranOffset)" << endl;
	for (map<string, double>::const_iterator b = items.begin(), e = items.end();
		b != e; ++b) {
		os << "  item name = " << (*b).first << ";  offset value = " << (*b).second << endl;
	}
}

///////////////// TranScale Methods /////////////////
void TranScale::forward(Transformable &data)
{
	Transformable::iterator data_iter, data_end = data.end();
	for (map<string,double>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
		data_iter = data.find(b->first);
		if (data_iter != data_end)
		{
			(*data_iter).second *= b->second;
		}
	}
}


void TranScale::reverse(Transformable &data)
{
	Transformable::iterator data_iter, data_end = data.end();
	for (map<string,double>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
		data_iter = data.find(b->first);
		if (data_iter != data_end)
		{
			(*data_iter).second /= b->second;
		}
	}
}

void TranScale::jacobian_forward(Jacobian &jac)
{
	size_t icol = 0;
	double factor = 0;
	Transformable &data = jac.base_numeric_parameters;
	unordered_map<string, int> par_2_col_map = jac.get_par2col_map();
	auto iter_end = par_2_col_map.end();
	for (const auto &irec : items)
	{
		auto iter = par_2_col_map.find(irec.first);
		if (iter != iter_end)
		{
			icol = iter->second;
			factor = irec.second;
			jac.matrix.col(icol) /= factor;
		}
	}
	forward(data);
}

void TranScale::jacobian_reverse(Jacobian &jac)
{
	size_t icol = 0;
	double factor = 0;
	Transformable &data = jac.base_numeric_parameters;
	unordered_map<string, int> par_2_col_map = jac.get_par2col_map();
	auto iter_end = par_2_col_map.end();
	for (const auto &irec : items)
	{
		auto iter = par_2_col_map.find(irec.first);
		if (iter != iter_end)
		{
			icol = iter->second;
			factor = irec.second;
			jac.matrix.col(icol) *= factor;
		}
	}
	reverse(data);
}

void TranScale::d1_to_d2(Transformable &del_data, Transformable &data)
{
	Transformable::iterator del_data_iter, del_data_end = del_data.end();
	for (map<string, double>::const_iterator b = items.begin(), e = items.end();
		b != e; ++b) {
		del_data_iter = del_data.find(b->first);
		if (del_data_iter != del_data_end)
		{
			(*del_data_iter).second /= b->second;
		}
	}
	forward(data);
}

void TranScale::d2_to_d1(Transformable &del_data, Transformable &data)
{
	Transformable::iterator del_data_iter, del_data_end = del_data.end();
	for (map<string, double>::const_iterator b = items.begin(), e = items.end();
		b != e; ++b) {
		del_data_iter = del_data.find(b->first);
		if (del_data_iter != del_data_end)
		{
			(*del_data_iter).second *= b->second;
		}
	}
	reverse(data);
}

void TranScale::print(ostream &os) const
{
	os << "Transformation name = " << name << "; (type=TranScale)" << endl;
	for (map<string,double>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
			os << "  item name = " << (*b).first << ";  scale value = " << (*b).second << endl;
	}
}

///////////////// TranLog10 Methods /////////////////
void TranLog10::forward(Transformable &data)
{
	Transformable::iterator data_iter, data_end = data.end();
	for (set<string>::const_iterator b=items.begin(), e=items.end(); b!=e; ++b)
	{
		data_iter = data.find(*b);
		if (data_iter != data_end)
		{
			(*data_iter).second = log10((*data_iter).second);
		}
	}
}

void TranLog10::reverse(Transformable &data)
{
	Transformable::iterator data_iter, data_end = data.end();
	for (set<string>::const_iterator b=items.begin(), e=items.end(); b!=e; ++b)
	{
		data_iter = data.find(*b);
		if (data_iter != data_end)
		{
			(*data_iter).second = pow(10.0, (*data_iter).second);
		}
	}
}



void TranLog10::jacobian_forward(Jacobian &jac)
{
	size_t icol = 0;
	double factor = 0;
	double d = 0;
	Transformable &data = jac.base_numeric_parameters;
	forward(data);
	unordered_map<string, int> par_2_col_map = jac.get_par2col_map();
	auto iter_end = par_2_col_map.end();
	for (const auto &ipar : items)
	{
		auto iter = par_2_col_map.find(ipar);
		if (iter != iter_end)
		{
			d = data.get_rec(ipar);
			icol = iter->second;
			factor = pow(10, d) * log(10.0);
		}
		jac.matrix.col(icol) *= factor;
	}
}

void TranLog10::jacobian_reverse(Jacobian &jac)
{
	size_t icol = 0;
	double factor = 0;
	double d = 0;
	Transformable &data = jac.base_numeric_parameters;
	reverse(data);
	unordered_map<string, int> par_2_col_map = jac.get_par2col_map();
	auto iter_end = par_2_col_map.end();
	for (const auto &ipar : items)
	{
		auto iter = par_2_col_map.find(ipar);
		if (iter != iter_end)
		{
			d = data.get_rec(ipar);
			icol = iter->second;
			factor = 1.0 / (d * log(10.0));
			jac.matrix.col(icol) *= factor;
		}
	}
}


void TranLog10::d1_to_d2(Transformable &del_data, Transformable &data)
{
	forward(data);
	Transformable::iterator del_data_iter, del_data_end = del_data.end();
	for (set<string>::const_iterator b = items.begin(), e = items.end(); b != e; ++b)
	{
		del_data_iter = del_data.find(*b);
		if (del_data_iter != del_data_end)
		{
			double d1 = data.get_rec(*b);
			double factor = pow(10.0, d1) * log(10.0);
			(*del_data_iter).second *= factor;
		}
	}
}

void TranLog10::d2_to_d1(Transformable &del_data, Transformable &data)
{
	reverse(data);
	Transformable::iterator del_data_iter, del_data_end = del_data.end();
	for (set<string>::const_iterator b = items.begin(), e = items.end(); b != e; ++b)
	{
		del_data_iter = del_data.find(*b);
		if (del_data_iter != del_data_end)
		{
			double d2 = data.get_rec(*b);
			double factor = 1.0 / (d2 * log(10.0));
			(*del_data_iter).second *= factor;
		}
	}
	reverse(data);
}

void TranLog10::print(ostream &os) const
{
	os << "Transformation name = " << name << "; (type=TranLog10)" << endl;
	for (set<string>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
			os << "  item name = " << *b << endl;
	}
}

///////////////// TranFixed Methods /////////////////
void TranFixed::forward(Transformable &data)
{
	for (map<string,double>::iterator b=items.begin(), e=items.end(); b!=e; ++b)
	{
		data.erase(b->first);
	}
}

void TranFixed::reverse(Transformable &data)
{
	for (map<string,double>::iterator b=items.begin(), e=items.end(); b!=e; ++b)
	{
		data.insert(b->first, b->second);
	}
}


void TranFixed::jacobian_forward(Jacobian &jac)
{
	Transformable &data = jac.base_numeric_parameters;
	set<string> rm_par_set;
	for (const auto &irec : items)
	{
		rm_par_set.insert(irec.first);
	}
	//remove Fixed parameters from base_parameter_names
	jac.remove_cols(rm_par_set);
	forward(data);
}

void TranFixed::jacobian_reverse(Jacobian &jac)
{
	Transformable &data = jac.base_numeric_parameters;
	set<string> new_pars;
	for (auto &i : items)
	{
		new_pars.insert(i.first);
	}
	jac.add_cols(new_pars);
	reverse(data);
}



void TranFixed::d1_to_d2(Transformable &del_data, Transformable &data)
{
	for (map<string, double>::iterator b = items.begin(), e = items.end(); b != e; ++b)
	{
		del_data.erase(b->first);
	}
	forward(data);
}

void TranFixed::d2_to_d1(Transformable &del_data, Transformable &data)
{
	for (map<string, double>::iterator b = items.begin(), e = items.end(); b != e; ++b)
	{
		del_data.insert(b->first, 0.0);
	}
	reverse(data);
}


void TranFixed::print(ostream &os) const
{
	os << "Transformation name = " << name << "; (type=TranFixed)" << endl;
	for (map<string,double>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
			os << "  item name = " << (*b).first << ";  imposed value = " << (*b).second << endl;
	}
}

void TranFrozen::print(ostream &os) const
{
	os << "Transformation name = " << name << "; (type=TranFrozen)" << endl;
	for (map<string,double>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
			os << "  item name = " << (*b).first << ";  imposed value = " << (*b).second << endl;
	}
}

void TranTied::insert(const string &item_name, const pair<string, double> &item_value)
{
	items[item_name] = item_value;
}


void TranTied::forward(Transformable &data)
{
	for (map<string, pair_string_double>::iterator ii = items.begin(); ii != items.end(); ++ii)
	{
		data.erase(ii->first);
	}
}


void TranTied::reverse(Transformable &data)
{
	string const *base_name;
	double *factor;
	Transformable::iterator base_iter;
	for (map<string, pair_string_double>::iterator b = items.begin(), e = items.end();
		b != e; ++b)
	{
		base_name = &(b->second.first);
		factor = &(b->second.second);
		base_iter = data.find(*base_name);
		if (base_iter != data.end())
		{
			//cout << b->first << ',' << data[b->first] <<  ',' << (*base_iter).second << endl;
			data.erase(b->first);
			data.insert(b->first, (*base_iter).second * (*factor));
			//cout << b->first << ',' << data[b->first] << ',' << (*base_iter).second << endl;
		}

	}
}

void TranTied::jacobian_forward(Jacobian &jac)
{
	throw(PestError("Error: TranTied::jacobian_forward - TranTied does not support Jacobian transformations"));
}



void TranTied::jacobian_reverse(Jacobian &jac)
{
	throw(PestError("Error: TranTied::jacobian_forward - TranTied does not support Jacobian transformations"));
}

void TranTied::d1_to_d2(Transformable &del_data, Transformable &data)
{
	throw(PestError("Error: TranTied::d1_to_d2 - TranTied does not support d1_to_d2 transformations"));
}

void TranTied::d2_to_d1(Transformable &del_data, Transformable &data)
{
	throw(PestError("Error: TranTied::d2_to_d1 - TranTied does not support d2_to_d1 transformations"));
}

void TranTied::print(ostream &os) const
{
	os << "Transformation name = " << name << "; (type=TranTied)" << endl;
	for (map<string, pair_string_double>::const_iterator b=items.begin(), e=items.end();
		b!=e; ++b) {
			os << "  item name = " << (*b).first << "   tied to \"" << (*b).second.first <<
				"\" with factor " <<  (*b).second.second<<  endl;
	}
}

const  Eigen::SparseMatrix<double>& TranSVD::get_vt() const
{
	return Vt;
}


TranSVD::TranSVD(int _max_sing, double _eign_thresh, const string &_name) : Transformation(_name)
{
	tran_svd_pack = new SVD_REDSVD(_max_sing, _eign_thresh);
}


TranSVD::TranSVD(const TranSVD& rhs)
	: Transformation(rhs), base_parameter_names(rhs.base_parameter_names),
	super_parameter_names(rhs.super_parameter_names),
	obs_names(rhs.obs_names),
	//SqrtQ_J(rhs.SqrtQ_J),
	jtqj(rhs.jtqj),
	Sigma(rhs.Sigma),
	U(rhs.U),
	Vt(rhs.Vt),
	init_base_numeric_parameters(rhs.init_base_numeric_parameters),
	frozen_derivative_parameters(rhs.frozen_derivative_parameters)
{
	tran_svd_pack = rhs.tran_svd_pack->clone();
}


void TranSVD::set_SVD_pack()
{
	int max_sing = tran_svd_pack->get_max_sing();
	double eigthresh = tran_svd_pack->get_eign_thres();
	delete tran_svd_pack;
	tran_svd_pack = new SVD_EIGEN(max_sing, eigthresh);
}

void TranSVD::set_performance_log(PerformanceLog *performance_log)
{
	tran_svd_pack->set_performance_log(performance_log);
}

void TranSVD::calc_svd()
{
	debug_msg("TranSVD::calc_svd begin");
	stringstream sup_name;
	VectorXd Sigma_trunc;
	//tran_svd_pack->solve_ip(SqrtQ_J, Sigma, U, Vt, Sigma_trunc);
	tran_svd_pack->solve_ip(jtqj, Sigma, U, Vt, Sigma_trunc);
	// calculate the number of singluar values above the threshold

	debug_print(Sigma);
	debug_print(U);
	debug_print(Vt);
	debug_print(Sigma_trunc);


	int n_sing_val = Sigma.size();

	super_parameter_names.clear();
	for(int i=0; i<n_sing_val; ++i) {
		sup_name.str("");
		sup_name << "SUP_";
		sup_name << i+1;
		super_parameter_names.push_back(sup_name.str());
	}
	if (n_sing_val <= 0 )
	{
		cout << jtqj.rows() << "," << jtqj.cols() << endl;
		throw PestError("TranSVD::update() - super parameter transformation returned 0 super parameters.  Jacobian must equal 0.");
	}

	Eigen::MatrixXd temp = Vt.toDense();
	//for (auto name : base_parameter_names)
	for (int i=0;i<base_parameter_names.size();i++)
	{
		if ((dss.find(base_parameter_names[i]) != dss.end()) && (abs(dss[base_parameter_names[i]]) < 1.0e-6))
		{
			temp.col(i).setZero();
		}
	}
	Vt = temp.sparseView();
	debug_print(super_parameter_names);
	debug_msg("TranSVD::calc_svd end");
}

void TranSVD::update_reset_frozen_pars(const Jacobian &jacobian, const QSqrtMatrix &Q_sqrt, const Parameters &base_numeric_pars,
		int maxsing, double _eigthresh, const vector<string> &par_names, const vector<string> &_obs_names, 
		Eigen::SparseMatrix<double>& parcov_inv, map<string,double> _dss,
		const Parameters &_frozen_derivative_pars)
{
	dss = _dss;
	debug_msg("TranSVD::update_reset_frozen_pars begin");
	debug_print(_frozen_derivative_pars);
	stringstream sup_name;
	super_parameter_names.clear();

	tran_svd_pack->set_max_sing(maxsing);
	tran_svd_pack->set_eign_thres(_eigthresh);
	obs_names = _obs_names;

	//these are where the derivative was computed so they can be different than the frozen values;
	init_base_numeric_parameters = base_numeric_pars;
	base_parameter_names = par_names;
	//remove frozen parameters from base_parameter_names
	auto end_iter = std::remove_if(base_parameter_names.begin(), base_parameter_names.end(),
		[&_frozen_derivative_pars](string &str)->bool{return _frozen_derivative_pars.find(str)!=_frozen_derivative_pars.end();});
	base_parameter_names.resize(std::distance(base_parameter_names.begin(), end_iter));
	frozen_derivative_parameters = _frozen_derivative_pars;
	//remove frozen derivatives from matrix parameter list
//	std::remove_if(base_parameter_names.begin(), base_parameter_names.end(),
//		[this](string &str)->bool{return this->frozen_derivative_parameters.find(str)!=this->frozen_derivative_parameters.end();});

	//SqrtQ_J = Q_sqrt.get_sparse_matrix(obs_names, DynamicRegularization::get_unit_reg_instance()) * jacobian.get_matrix(obs_names, base_parameter_names);
	Eigen::SparseMatrix<double> j = jacobian.get_matrix(obs_names, base_parameter_names);
	jtqj = j.transpose() * Q_sqrt.get_sparse_matrix(obs_names, DynamicRegularization::get_unit_reg_instance(), true) * j;
	//if (parcov.ncol() > 0)
	if (parcov_inv.rows() > 0)
	{
		vector<string> numeric_par_names = jacobian.get_base_numeric_par_names();
		Eigen::MatrixXd lamb = Eigen::MatrixXd::Ones(jtqj.rows(), jtqj.cols());
		lamb = lamb + parcov_inv.toDense();
		for (int i = 0; i < numeric_par_names.size(); i++)
		{
			if (abs(dss[numeric_par_names[i]]) < 1.0e-6)
			{
				lamb.col(i).setZero();
				lamb.row(i).setZero();
			}
		}
		lamb = lamb + jtqj.toDense();
		jtqj = lamb.sparseView();
		lamb.resize(0, 0);
	}
	calc_svd();

	debug_print(this->base_parameter_names);
	debug_print(this->frozen_derivative_parameters);
	debug_msg("TranSVD::update_reset_frozen_pars end");
}

void TranSVD::update_add_frozen_pars(const Parameters &frozen_pars)
{
	debug_msg("TranSVD::update_reset_frozen_pars begin");
	debug_print(frozen_pars);
	vector<size_t> del_col_ids;
	Parameters new_frozen_pars;
	for (auto &ipar : frozen_pars)
	{
		auto iter = frozen_derivative_parameters.find(ipar.first);
		if (iter == frozen_derivative_parameters.end())
		{
			new_frozen_pars.insert(ipar);
		}
	}

	frozen_derivative_parameters.insert(new_frozen_pars);

	// build list of columns that needs to be removed from the matrix
	for (int i=0; i<base_parameter_names.size(); ++i)
	{
		if(new_frozen_pars.find(base_parameter_names[i]) != new_frozen_pars.end())
		{
			del_col_ids.push_back(i);
		}
	}

	if (del_col_ids.size() == base_parameter_names.size())
	{
		throw PestError("TranSVD::update_add_frozen_pars - All parameters are frozen in SVD transformation");
	}
	if ((del_col_ids.size() == 0) && (frozen_pars.size() > 0))
		cout << "WARNING: TravSVD::update_add_froze_pars(): all pars to freeze already frozen" << endl;
	//remove frozen parameters from base_parameter_names
	auto end_iter = std::remove_if(base_parameter_names.begin(), base_parameter_names.end(),
		[&new_frozen_pars](string &str)->bool{return new_frozen_pars.find(str)!=new_frozen_pars.end();});
	base_parameter_names.resize(std::distance(base_parameter_names.begin(), end_iter));
	//matrix_del_rows_cols(SqrtQ_J, del_col_ids,false,true);
	matrix_del_rows_cols(jtqj, del_col_ids, true, true);
	calc_svd();
	debug_print(this->base_parameter_names);
	debug_print(this->frozen_derivative_parameters);
	debug_msg("TranSVD::update_reset_frozen_pars end");
}
void TranSVD::reverse(Transformable &data)
{
	// Transform super-parameters to base parameters
	assert(Vt.cols() == base_parameter_names.size());
	int n_base = Vt.cols();
	vector<double> super_par_vec = data.get_data_vec(super_parameter_names);
	vector<double>::iterator it;
	for (it=super_par_vec.begin(); it!=super_par_vec.end(); ++it)
	{
		(*it) -= 10.0;
	}
	Transformable ret_base_pars;
	int n_sing_val = Sigma.size();
	VectorXd delta_base_mat = Vt.block(0,0,n_sing_val, Vt.cols()).transpose() *  stlvec_2_eigenvec(super_par_vec);
	for (int i=0; i<n_base; ++i) {
		ret_base_pars.insert(base_parameter_names[i], delta_base_mat(i) + init_base_numeric_parameters.get_rec(base_parameter_names[i]));
	}

	data = ret_base_pars;
}

void TranSVD::forward(Transformable &data)
{
	//Transform base parameters to super-parameters
	Transformable super_pars;
	VectorXd value;

	Transformable delta_data = init_base_numeric_parameters;
	delta_data *= 0.0;

	for (auto &it : data)
	{
		delta_data[it.first] = it.second - init_base_numeric_parameters.get_rec(it.first);
	}
	int n_sing_val = Sigma.size();
	value = Vt * delta_data.get_data_eigen_vec(base_parameter_names);
	for (int i=0; i<n_sing_val; ++i) {
		super_pars.insert(super_parameter_names[i], value(i)+10.0);
	}
	data = super_pars;
}

void TranSVD::jacobian_forward(Jacobian &jac)
{
	Transformable &data = jac.base_numeric_parameters;
	Eigen::SparseMatrix<double> old_matrix = jac.get_matrix(jac.observation_list(), base_parameter_names);
	Eigen::SparseMatrix<double> super_jacobian;
	super_jacobian = old_matrix * Vt.transpose();
	jac.matrix = super_jacobian;
	jac.base_numeric_par_names = super_parameter_names;

	forward(data);
}

void TranSVD::jacobian_reverse(Jacobian &jac)
{
	Transformable &data = jac.base_numeric_parameters;
	Eigen::SparseMatrix<double> old_matrix = jac.get_matrix(jac.observation_list(), super_parameter_names);
	Eigen::SparseMatrix<double> base_jacobian;
	base_jacobian = old_matrix * Vt;
	jac.matrix = base_jacobian;
	jac.base_numeric_par_names = base_parameter_names;
	reverse(data);
}

void TranSVD::d1_to_d2(Transformable &del_data, Transformable &data)
{
	Transformable new_data;
	VectorXd d1_vec = del_data.get_partial_data_eigen_vec(base_parameter_names);
	VectorXd d2_vec = Vt * d1_vec;
	for (int i = 0; i<Sigma.size(); ++i) {
		new_data[super_parameter_names[i]] = d2_vec[i];
	}
	del_data = new_data;
	forward(data);
}


void TranSVD::d2_to_d1(Transformable &del_data, Transformable &data)
{
	Transformable new_data;
	VectorXd d2_vec = del_data.get_partial_data_eigen_vec(super_parameter_names);
	VectorXd d1_vec = Vt.transpose() * d2_vec;
	for (size_t i = 0; i < base_parameter_names.size(); ++i)
	{
		new_data[base_parameter_names[i]] = d1_vec[i];
	}
	del_data = new_data;
	reverse(data);
}

void TranSVD::save(ostream &fout) const
{
	size_t size;
	vector<int8_t> serial_data;
	serial_data = Serialization::serialize(base_parameter_names);
	size = serial_data.size();
	fout.write((char*)&size, sizeof(size));
	fout.write((char*)serial_data.data(), size);

	serial_data = Serialization::serialize(super_parameter_names);
	size = serial_data.size();
	fout.write((char*)&size, sizeof(size));
	fout.write((char*)serial_data.data(), size);

	serial_data = Serialization::serialize(obs_names);
	size = serial_data.size();
	fout.write((char*)&size, sizeof(size));
	fout.write((char*)serial_data.data(), size);

	//save_triplets_bin(SqrtQ_J, fout);
	save_triplets_bin(jtqj, fout);
	save_vector_bin(Sigma, fout);
	save_triplets_bin(U, fout);
	save_triplets_bin(Vt, fout);

	serial_data = Serialization::serialize(init_base_numeric_parameters);
	size = serial_data.size();
	fout.write((char*)&size, sizeof(size));
	fout.write((char*)serial_data.data(), size);

	serial_data = Serialization::serialize(frozen_derivative_parameters);
	size = serial_data.size();
	fout.write((char*)&size, sizeof(size));
	fout.write((char*)serial_data.data(), size);
}

void TranSVD::read(istream &fin)
{
	size_t size;
	vector<int8_t> serial_data;

	base_parameter_names.clear();
	fin.read((char*)&size, sizeof(size));
	serial_data.clear();
	serial_data.resize(size);
	fin.read((char*)serial_data.data(), size);
	Serialization::unserialize(serial_data, base_parameter_names);

	super_parameter_names.clear();
	fin.read((char*)&size, sizeof(size));
	serial_data.clear();
	serial_data.resize(size);
	fin.read((char*)serial_data.data(), size);
	Serialization::unserialize(serial_data, super_parameter_names);

	obs_names.clear();
	fin.read((char*)&size, sizeof(size));
	serial_data.clear();
	serial_data.resize(size);
	fin.read((char*)serial_data.data(), size);
	Serialization::unserialize(serial_data, obs_names);

	//load_triplets_bin(SqrtQ_J, fin);
	load_triplets_bin(jtqj, fin);
	load_vector_bin(Sigma, fin);
	load_triplets_bin(U, fin);
	load_triplets_bin(Vt, fin);

	init_base_numeric_parameters.clear();
	fin.read((char*)&size, sizeof(size));
	serial_data.clear();
	serial_data.resize(size);
	fin.read((char*)serial_data.data(), size);
	Serialization::unserialize(serial_data, init_base_numeric_parameters);

	frozen_derivative_parameters.clear();
	fin.read((char*)&size, sizeof(size));
	serial_data.clear();
	serial_data.resize(size);
	fin.read((char*)serial_data.data(), size);
	Serialization::unserialize(serial_data, frozen_derivative_parameters);
}

ParameterGroupInfo TranSVD::build_par_group_info(const ParameterGroupInfo &base_pg_info)
{
	double derinc_sup;
	double derinc_par;
	int max_col;
	double max_val;
	ParameterGroupInfo pg_info;
	stringstream grp_name;
	for (int i_sup=0, n_sup=super_parameter_names.size(); i_sup < n_sup; ++i_sup)
	{
		get_MatrixXd_row_abs_max(Vt, i_sup, &max_col, &max_val);
		derinc_par = base_pg_info.get_group_rec_ptr(base_parameter_names[max_col])->derinc;
		derinc_sup = .01;
		grp_name.str("");
		grp_name << "g_" << super_parameter_names[i_sup];
		ParameterGroupRec sup_rec(grp_name.str(), "ABSOLUTE", derinc_sup, 0.0, "SWITCH", 2.0, "PARABOLIC");
		//add new group
		pg_info.insert_group(grp_name.str(), sup_rec);
		// connect super parameter to new group
		pg_info.insert_parameter_link(super_parameter_names[i_sup], grp_name.str());
	}
	return pg_info;
}


Parameters TranSVD::map_basepar_to_super(const Parameters &base_pars)
{
	Parameters super_pars;
	VectorXd base_par_vec = base_pars.get_partial_data_eigen_vec(base_parameter_names);
	VectorXd super_par_vec = Vt * base_par_vec;
	for (size_t i=0; i<super_parameter_names.size(); ++i)
	{
		super_pars[super_parameter_names[i]] = super_par_vec(i);
	}
	return super_pars;
}

TranSVD::~TranSVD()
{
	delete tran_svd_pack;
}

void TranSVD::print(ostream &os) const
{
	os << "Transformation name = " << name << "; (type=TranSVD)" << endl;
	os << "  Singular Values = " << Sigma << endl;
}

//void TranNormalize::forward(Transformable &data)
//{
//	Transformable::iterator data_iter, data_end = data.end();
//	for (map<string,NormData>::const_iterator b=items.begin(), e=items.end();
//		b!=e; ++b) {
//		data_iter = data.find(b->first);
//		if (data_iter != data_end)
//		{
//			(*data_iter).second += b->second.offset;
//			(*data_iter).second *= b->second.scale;
//		}
//	}
//}
//
//
//void TranNormalize::reverse(Transformable &data)
//{
//	Transformable::iterator data_iter, data_end = data.end();
//	for (map<string,NormData>::const_iterator b=items.begin(), e=items.end();
//		b!=e; ++b) {
//		data_iter = data.find(b->first);
//		if (data_iter != data_end)
//		{
//			(*data_iter).second /= b->second.scale;
//			(*data_iter).second -= b->second.offset;
//		}
//	}
//}
//
//
//void TranNormalize::jacobian_forward(Jacobian &jac)
//{
//	size_t icol = 0;
//	double factor = 0;
//	Transformable &data = jac.base_numeric_parameters;
//	unordered_map<string, int> par_2_col_map = jac.get_par2col_map();
//	auto iter_end = par_2_col_map.end();
//	for (const auto &irec : items)
//	{
//		auto iter = par_2_col_map.find(irec.first);
//		if (iter != iter_end)
//		{
//			icol = iter->second;
//			factor = irec.second.scale;
//			jac.matrix.col(icol) /= factor;
//		}
//	}
//	forward(data);
//}
//
//void TranNormalize::jacobian_reverse(Jacobian &jac)
//{
//	size_t icol = 0;
//	double factor = 0;
//	Transformable &data = jac.base_numeric_parameters;
//	unordered_map<string, int> par_2_col_map = jac.get_par2col_map();
//	auto iter_end = par_2_col_map.end();
//	for (const auto &irec : items)
//	{
//		auto iter = par_2_col_map.find(irec.first);
//		if (iter != iter_end)
//		{
//			icol = iter->second;
//			factor = irec.second.scale;
//			jac.matrix.col(icol) *= factor;
//		}
//	}
//	reverse(data);
//}
//
//void TranNormalize::d1_to_d2(Transformable &del_data, Transformable &data)
//{
//	Transformable::iterator del_data_iter, del_data_end = del_data.end();
//	for (map<string, NormData>::const_iterator b = items.begin(), e = items.end();
//		b != e; ++b) {
//		del_data_iter = data.find(b->first);
//		if (del_data_iter != del_data_end)
//		{
//			(*del_data_iter).second /= b->second.scale;
//		}
//	}
//	forward(data);
//}
//
//void TranNormalize::d2_to_d1(Transformable &del_data, Transformable &data)
//{
//	Transformable::iterator del_data_iter, del_data_end = del_data.end();
//	for (map<string, NormData>::const_iterator b = items.begin(), e = items.end();
//		b != e; ++b) {
//		del_data_iter = del_data.find(b->first);
//		if (del_data_iter != del_data_end)
//		{
//			(*del_data_iter).second *= b->second.scale;
//		}
//	}
//	reverse(data);
//}
//
//void TranNormalize::insert(const string &item_name, double _offset, double _scale)
//{
//	items[item_name] = NormData(_offset, _scale);
//}
//
//void TranNormalize::print(ostream &os) const
//{
//	os << "Transformation name = " << name << "; (type=TranNormalize)" << endl;
//	for (map<string,NormData>::const_iterator b=items.begin(), e=items.end();
//		b!=e; ++b) {
//			os << "  item name = " << (*b).first << ";  scale value = " << (*b).second.scale
//				<< ";  offset value = " << (*b).second.offset <<endl;
//	}
//}

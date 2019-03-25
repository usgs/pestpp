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
#ifndef REGULARIZATION_H_
#define REGULARIZATION_H_

#include <unordered_map>
#include <string>
class ModelRun;

class DynamicRegularization
{
public:
	DynamicRegularization(bool _use_dynamic_reg = false, bool _grp_weight_adj = false, double _phi_m_lim = 0,
		double _phi_m_accept = 0, double _frac_phi_m = 1, double _wf_min = 1e-10, double _wf_max = 1e10,
		double _wffac = 0, double _wftol = 1000, double _wf_init = 1.0, int _max_reg_iter=20);
	DynamicRegularization(const DynamicRegularization &rhs);
	virtual double get_weight() const;
	virtual int get_max_reg_iter() const { return max_reg_iter; }
	virtual double get_phimlim() const { return phi_m_lim; }
	virtual double get_phimaccept() const { return phi_m_accept; }
	virtual double get_fracphim() const { return frac_phi_m; }
	virtual double get_wfmin() const { return wf_min; }
	virtual double get_wfmax() const { return wf_max; }
	virtual double get_wffac() const { return wffac; }
	virtual double get_wftol() const { return wftol; }
	virtual double get_wfinit() const { return wf_init; }
	virtual bool get_use_dynamic_reg() const { return use_dynamic_reg; }
	virtual bool get_adj_grp_weights() const { return adj_grp_weights; }
	virtual double get_grp_weight_fact(const std::string &grp_name) const;
	virtual void set_weight(double _tikhonov_weight) {tikhonov_weight = _tikhonov_weight;}
	virtual void set_max_reg_iter(int _max_reg_iter) { max_reg_iter = _max_reg_iter; }
	virtual void set_regul_grp_weights(const std::unordered_map<std::string, double> &_regul_grp_weights) { regul_grp_weights = _regul_grp_weights; }
	static DynamicRegularization get_unit_reg_instance() { return DynamicRegularization(); }
	static DynamicRegularization get_zero_reg_instance() { return DynamicRegularization(true, false, 0,0,0,0,0,0,0,0); }
	virtual ~DynamicRegularization(void){}
protected:
	bool use_dynamic_reg;
	bool adj_grp_weights;
	int max_reg_iter;
	double phi_m_lim;
	double phi_m_accept;
	double frac_phi_m;
	double wf_min;
	double wf_max;
	double wffac;
	double wftol;
	double wf_init;
	double tikhonov_weight;
	std::unordered_map<std::string, double> regul_grp_weights;
};


#endif //REGULARIZATIONABSTRACT_H_

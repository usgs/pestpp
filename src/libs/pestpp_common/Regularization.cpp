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
#include "Regularization.h"
#include "ObjectiveFunc.h"
#include "ModelRunPP.h"
#include "utilities.h"
#include "pest_data_structs.h"


DynamicRegularization::DynamicRegularization(bool _use_dynamic_reg, bool _adj_grp_weights, double _phi_m_lim,
	double _phi_m_accept, double _frac_phi_m, double _wf_min, double _wf_max,
	double _wffac, double _wftol, double _wf_init, int _max_reg_iter)
	: use_dynamic_reg(_use_dynamic_reg), adj_grp_weights(_adj_grp_weights), phi_m_lim(_phi_m_lim), phi_m_accept(_phi_m_accept), frac_phi_m(_frac_phi_m),
	wf_min(_wf_min), wf_max(_wf_max), wffac(_wffac), wftol(_wftol), wf_init(_wf_init),
	tikhonov_weight(_wf_init), max_reg_iter(_max_reg_iter)
{
}

DynamicRegularization::DynamicRegularization(const DynamicRegularization &rhs)
	: use_dynamic_reg(rhs.use_dynamic_reg), adj_grp_weights(rhs.adj_grp_weights), phi_m_lim(rhs.phi_m_lim), phi_m_accept(rhs.phi_m_accept), frac_phi_m(rhs.frac_phi_m),
	wf_min(rhs.wf_min), wf_max(rhs.wf_max), wffac(rhs.wffac), wftol(rhs.wftol), wf_init(rhs.wf_init),
	tikhonov_weight(rhs.wf_init), max_reg_iter(rhs.max_reg_iter), regul_grp_weights(rhs.regul_grp_weights)
{
}
double DynamicRegularization::get_weight() const
{
	return tikhonov_weight;
}

double DynamicRegularization::get_grp_weight_fact(const std::string &grp_name) const
{
	double fact = 1.0;
	auto iter = regul_grp_weights.find(grp_name);
	if (use_dynamic_reg && adj_grp_weights && iter != regul_grp_weights.end())
	{
		fact = iter->second;
	}
	return fact;
}

PestppOptions::ARG_STATUS DynamicRegularization::assign_value_by_key(const std::string key, const std::string org_value)
{
	/*DynamicRegularization(bool _use_dynamic_reg = false, bool _grp_weight_adj = false, double _phi_m_lim = 0,
		double _phi_m_accept = 0, double _frac_phi_m = 1, double _wf_min = 1e-10, double _wf_max = 1e10,
		double _wffac = 0, double _wftol = 1000, double _wf_init = 1.0, int _max_reg_iter = 20);*/
	std::string value = pest_utils::upper_cp(org_value);
	if (passed_args.find(key) != passed_args.end())
		return PestppOptions::ARG_STATUS::ARG_DUPLICATE;
	passed_args.insert(key);
	if (key == "PHIMLIM")
		pest_utils::convert_ip(value, phi_m_lim);


	return PestppOptions::ARG_STATUS::ARG_ACCEPTED;
}

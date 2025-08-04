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


//DynamicRegularization::DynamicRegularization(bool _use_dynamic_reg, bool _adj_grp_weights, double _phi_m_lim,
//	double _phi_m_accept, double _frac_phi_m, double _wf_min, double _wf_max,
//	double _wffac, double _wftol, double _wf_init, int _max_reg_iter)
//	: use_dynamic_reg(_use_dynamic_reg), adj_grp_weights(_adj_grp_weights), phi_m_lim(_phi_m_lim), phi_m_accept(_phi_m_accept), frac_phi_m(_frac_phi_m),
//	wf_min(_wf_min), wf_max(_wf_max), wffac(_wffac), wftol(_wftol), wf_init(_wf_init),
//	tikhonov_weight(_wf_init), max_reg_iter(_max_reg_iter)
//{
//}

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

DynamicRegularization DynamicRegularization::get_unit_reg_instance()
{
	DynamicRegularization new_reg;
	new_reg.set_defaults();
	new_reg.set_weight(1.0);
	return new_reg;
}

DynamicRegularization DynamicRegularization::get_zero_reg_instance()
{
	DynamicRegularization new_reg;
	new_reg.set_zero();
	//{ return DynamicRegularization(true, false, 0,0,0,0,0,0,0,0); }
	/*new_reg.assign_value_by_key("PHIMLIM", 0);
	new_reg.assign_value_by_key("PHIMACCEPT", 0);
	new_reg.assign_value_by_key("FRACPHIM", 0);
	new_reg.assign_value_by_key("PHIMLIM", 0);
	new_reg.assign_value_by_key("WFINIT", 0);
	new_reg.assign_value_by_key("WFMAX", 0);
	new_reg.assign_value_by_key("WFMIN", 0);
	new_reg.assign_value_by_key("WFFAC", 0);
	new_reg.assign_value_by_key("WFTOL", 0);*/

	//new_reg.use_dynamic_reg = true;
	//new_reg.adj_grp_weights = false;

	return new_reg;
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
	else if (key == "PHIMACCEPT")
		pest_utils::convert_ip(value, phi_m_accept);
	else if (key == "FRACPHIM")
		pest_utils::convert_ip(value, frac_phi_m);
	else if (key == "WFMIN")
		pest_utils::convert_ip(value, wf_min);
	else if (key == "WFMAX")
		pest_utils::convert_ip(value, wf_max);
	else if (key == "WFINIT") {
        pest_utils::convert_ip(value, wf_init);
        tikhonov_weight = wf_init;
    }

	else if (key == "WFTOL")
		pest_utils::convert_ip(value, wftol);
	else if (key == "WFFAC")
		pest_utils::convert_ip(value, wffac);
	else if (key == "MAX_REG_ITER")
		pest_utils::convert_ip(value, max_reg_iter);
	else if (key == "IREGADJ")
	{
		int temp;
		pest_utils::convert_ip(value, temp);
		if (temp == 1)
			adj_grp_weights = true;
	}
		else
	{
		return PestppOptions::ARG_STATUS::ARG_NOTFOUND;
	}

	return PestppOptions::ARG_STATUS::ARG_ACCEPTED;
}

void DynamicRegularization::set_defaults()
{
	/*DynamicRegularization(bool _use_dynamic_reg = false, bool _grp_weight_adj = false, double _phi_m_lim = 0,
		double _phi_m_accept = 0, double _frac_phi_m = 1, double _wf_min = 1e-10, double _wf_max = 1e10,
		double _wffac = 0, double _wftol = 1000, double _wf_init = 1.0, int _max_reg_iter = 20);*/
	use_dynamic_reg = false;
	adj_grp_weights = false;
	phi_m_lim = 1.0e-10;
	phi_m_accept = 1.1e-10;
	frac_phi_m = 1;
	wf_min = 1e-10;
	wf_max = 1e+10;
	wffac = 1.3;
	wftol = 0.01;
	wf_init = 1.0;
	tikhonov_weight = 1.0;
	max_reg_iter = 5;
}

void DynamicRegularization::set_zero()
{
	use_dynamic_reg = true;
	adj_grp_weights = false;
	phi_m_lim = 0;
	phi_m_accept = 0;
	frac_phi_m = 0;
	wf_min = 0;
	wf_max = 0;
	wffac = 0;
	wftol = 0;
	wf_init = 0;
	max_reg_iter = 0;
}

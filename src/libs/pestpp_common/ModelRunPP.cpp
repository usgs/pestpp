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
 * @file ModelRunPP.cpp
 * @brief Implementation of ModelRunPP.
 */

#include <vector>
#include <map>
#include "utilities.h"
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <utility>
#include "ModelRunPP.h"
#include "Transformable.h"
#include "system_variables.h"
#include "pest_error.h"
#include "Transformation.h"

using namespace std;
using namespace pest_utils;

/**
 * @brief Model run.
 *
 * @param _obj_func_ptr Description.
 * @param _sim_obs Description.
 *
 * @return Description.
 */
ModelRun::ModelRun(const ObjectiveFunc *_obj_func_ptr, const Observations &_sim_obs)
	: obj_func_ptr(_obj_func_ptr), sim_obs(_sim_obs), obs_is_valid(false)
{
}

ModelRun& ModelRun::operator=(const ModelRun &rhs)
{
	frozen_ctl_par_names = rhs.frozen_ctl_par_names;
	obj_func_ptr = rhs.obj_func_ptr;
	ctl_pars = rhs.ctl_pars;
	sim_obs = rhs.sim_obs;
	obs_is_valid = rhs.obs_is_valid;
	return *this;
}


/**
 * @brief Get frozen ctl pars.
 *
 * @return Description.
 */
Parameters ModelRun::get_frozen_ctl_pars() const
{
	Parameters frz_pars(ctl_pars, frozen_ctl_par_names);
	return frz_pars;
}

/**
 * @brief Set frozen ctl parameters.
 *
 * @param frz_pars Description.
 */
void ModelRun::set_frozen_ctl_parameters(const Parameters &frz_pars)
{
	frozen_ctl_par_names.clear();
	add_frozen_ctl_parameters(frz_pars);
}

/**
 * @brief Add frozen ctl parameters.
 *
 * @param frz_pars Description.
 */
void ModelRun::add_frozen_ctl_parameters(const Parameters &frz_pars)
{
	for (auto &ipar : frz_pars)
	{
		frozen_ctl_par_names.push_back(ipar.first);
		ctl_pars[ipar.first] = ipar.second;
	}
}

/**
 * @brief Clear frozen ctl parameters.
 */
void ModelRun::clear_frozen_ctl_parameters()
{
	frozen_ctl_par_names.clear();
}

/**
 * @brief Set ctl parameters.
 *
 * @param pars Description.
 */
void ModelRun::set_ctl_parameters(const Parameters &pars)
{
	ctl_pars = pars;
}

/**
 * @brief Update ctl.
 *
 * @param ctl_pars Description.
 * @param obs Description.
 */
void ModelRun::update_ctl(Parameters &ctl_pars, Observations &obs)
{
	set_ctl_parameters(ctl_pars);
	// Process Observations
	set_observations(obs);
}

/**
 * @brief Set observations.
 *
 * @param observations Description.
 */
void ModelRun::set_observations(const Observations &observations)
{
	sim_obs = observations;
	obs_is_valid = true;
}

const Parameters &ModelRun::get_ctl_pars() const
{
	return ctl_pars;
}

const Observations &ModelRun::get_obs() const
{
	if( !obs_is_valid) {
		throw PestError("ModelRun::get_obs() - observations is invalid");
	}
	return sim_obs;
}

/**
 * @brief Get obs template.
 *
 * @return Description.
 */
Observations ModelRun::get_obs_template() const
{
	Observations ret_val(sim_obs);
	for (Observations::iterator b=ret_val.begin(), e=ret_val.end();
		b!=e; ++b) {
			b->second = Observations::no_data;
	}
	return ret_val;
}

/**
 * @brief Get phi.
 *
 * @param dynamic_reg Description.
 * @param norm Description.
 *
 * @return Description.
 */
double ModelRun::get_phi(const DynamicRegularization &dynamic_reg, double norm)
{
	PhiComponets phi_comp_temp = get_phi_comp(dynamic_reg, norm);
	return phi_comp_temp.meas + phi_comp_temp.regul;
}


/**
 * @brief Get phi comp.
 *
 * @param dynamic_reg Description.
 * @param norm Description.
 *
 * @return Description.
 */
PhiComponets ModelRun::get_phi_comp(const DynamicRegularization &dynamic_reg, double norm) const
{
	PhiComponets phi_comp = obj_func_ptr->get_phi_comp(sim_obs, get_ctl_pars(), dynamic_reg, norm);
	return phi_comp;
}


/**
 * @brief Get residuals vec.
 *
 * @param obs_names Description.
 *
 * @return Description.
 */
vector<double> ModelRun::get_residuals_vec(const vector<string> &obs_names)
{
	return obj_func_ptr->get_residuals_vec(get_obs(), get_ctl_pars(), obs_names);
}

//void ModelRun::full_report(ostream &os, const DynamicRegularization &dynamic_reg,bool limit_par)
//{
//	if( !obs_is_valid) {
//		throw PestError("ModelRun::full_report() - Simulated observations are invalid.  Can not produce full report.");
//	}
//	PhiComponets phi_comp = obj_func_ptr->full_report(os, sim_obs, get_ctl_pars(), dynamic_reg,limit_par);
//	output
//
//}


/**
 * @brief Obs valid.
 *
 * @return Description.
 */
bool ModelRun::obs_valid() const
{
	return obs_is_valid;
}

/**
 * @brief Cmp lt.
 *
 * @param r1 Description.
 * @param r2 Description.
 * @param reg Description.
 *
 * @return Description.
 */
bool ModelRun::cmp_lt(const ModelRun &r1, const ModelRun &r2, const DynamicRegularization &reg)
{
	bool cmp_flg = 0;
	bool use_dyamic_reg = reg.get_use_dynamic_reg();
	double phi_accept = reg.get_phimaccept();
	PhiComponets phi_1 = r1.get_phi_comp(reg);
	PhiComponets phi_2 = r2.get_phi_comp(reg);
	if (!use_dyamic_reg)
	{
		PhiComponets phi_1n = r1.get_phi_comp(DynamicRegularization::get_unit_reg_instance());
		PhiComponets phi_2n = r2.get_phi_comp(DynamicRegularization::get_unit_reg_instance());
		cmp_flg = ((phi_1n.meas + phi_1n.regul) < (phi_2n.meas + phi_2n.regul));
	}
	else if (phi_1.meas > phi_accept && phi_2.meas > phi_accept)
	{
		cmp_flg = ((phi_1.meas + phi_1.regul) < (phi_2.meas + phi_2.regul));
	}
	else if (phi_1.meas > phi_accept || phi_2.meas > phi_accept)
	{
		cmp_flg = (phi_1.meas < phi_2.meas) ? true : false;
	}
	else
	{
		cmp_flg = (phi_1.regul < phi_2.regul) ? true : false;
	}
	return cmp_flg;
}

/**
 * @brief Destructor for .
 */
ModelRun::~ModelRun()
{
}

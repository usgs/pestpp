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
#ifndef OBJECTIVEFUNC_H_
#define OBJECTIVEFUNC_H_

#include <list>
#include <utility>
#include "pest_data_structs.h"
#include "Pest.h"
#include "Transformable.h"

using namespace std;

class PhiComponets
{
public:
	enum class OBS_TYPE { ALL, REGUL, MEAS };
	PhiComponets() : meas(0.0), regul(0.0) {}
	PhiComponets(const PhiComponets &rhs) {*this = rhs;}
	const PhiComponets& operator=(const PhiComponets &rhs);
	double meas;
	double regul;
	pair<double, double> as_pair(){ return pair<double, double>(meas, regul); }
};

class PhiData
{
public:
	double meas;
	double regul;
	double total() const { return meas + regul;}
	map<string, double> group_phi;
};

class ObjectiveFunc
{
public:
	ObjectiveFunc() { ; }

	ObjectiveFunc(const Observations *_observations_ptr, const ObservationInfo *_obs_info_ptr, const PriorInformation *_prior_info_ptr)
		: observations_ptr(_observations_ptr), obs_info_ptr(_obs_info_ptr), prior_info_ptr(_prior_info_ptr) {}

	ObjectiveFunc(const Observations *_observations_ptr, const ObservationInfo *_obs_info_ptr, const PriorInformation *_prior_info_ptr,
			      const Pest *_ctl_file_ptr)
		: observations_ptr(_observations_ptr), obs_info_ptr(_obs_info_ptr), prior_info_ptr(_prior_info_ptr),ctl_file_ptr(_ctl_file_ptr) {}

	double get_phi(const Observations &sim_obs, const Parameters &pars, const DynamicRegularization &dynamic_reg, double norm = 2.0) const;
	PhiComponets get_phi_comp(const Observations &sim_obs, const Parameters &pars, const DynamicRegularization &dynamic_reg, double norm = 2.0) const;
	map<string, double> get_group_phi(const Observations &sim_obs, const Parameters &pars, const DynamicRegularization &dynamic_reg,
		PhiComponets::OBS_TYPE obs_type = PhiComponets::OBS_TYPE::ALL) const;
	PhiData phi_report(const Observations &sim_obs, const Parameters &pars, const DynamicRegularization &dynamic_reg) const;
	//PhiComponets full_report(ostream &os, const Observations &sim_obs, const Parameters &pars, const DynamicRegularization &dynamic_reg,bool limit_par=false) const;
	vector<double> get_residuals_vec(const Observations &sim_obs, const Parameters &pars, const vector<string> &obs_names) const;
	const Observations* get_obs_ptr() const;
	const ObservationInfo* get_obs_info_ptr() const;
	const PriorInformation* get_prior_info_ptr() const;
	~ObjectiveFunc(void) {}
private:
	const Observations *observations_ptr;
	const ObservationInfo *obs_info_ptr;
	const PriorInformation *prior_info_ptr;
	const Pest *ctl_file_ptr;
};

class MultiObjFuncs
{
public:

	MultiObjFuncs(const Observations* _observations_ptr, const ObservationInfo* _obs_info_ptr, const PriorInformation* _prior_info_ptr,
		const Pest* _ctl_file_ptr)
		: observations_ptr(_observations_ptr), obs_info_ptr(_obs_info_ptr), prior_info_ptr(_prior_info_ptr), ctl_file_ptr(_ctl_file_ptr) {}

	double get_phi(const Observations& sim_obs, const Parameters& pars, const DynamicRegularization& dynamic_reg, double norm = 2.0) const;
	PhiComponets get_phi_comp(const Observations& sim_obs, const Parameters& pars, const DynamicRegularization& dynamic_reg, double norm = 2.0) const;
	map<string, double> get_group_phi(const Observations& sim_obs, const Parameters& pars, const DynamicRegularization& dynamic_reg,
		PhiComponets::OBS_TYPE obs_type = PhiComponets::OBS_TYPE::ALL) const;
	PhiData phi_report(const Observations& sim_obs, const Parameters& pars, const DynamicRegularization& dynamic_reg) const;
	//PhiComponets full_report(ostream &os, const Observations &sim_obs, const Parameters &pars, const DynamicRegularization &dynamic_reg,bool limit_par=false) const;
	vector<double> get_residuals_vec(const Observations& sim_obs, const Parameters& pars, const vector<string>& obs_names) const;
	const Observations* get_obs_ptr() const;
	const ObservationInfo* get_obs_info_ptr() const;
	const PriorInformation* get_prior_info_ptr() const;
	~MultiObjFuncs(void) {}
private:
	const Observations* observations_ptr;
	const ObservationInfo* obs_info_ptr;
	const PriorInformation* prior_info_ptr;
	const Pest* ctl_file_ptr;
};

#endif /* OBJECTIVEFUNC_H_ */

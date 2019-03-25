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

#ifndef MODELRUNPP_H_
#define MODELRUNPP_H_

#include <list>
#include <utility>
#include <string>
#include "Transformable.h"
#include "ObjectiveFunc.h"
#include "ParamTransformSeq.h"

using namespace std;

class ModelRun
{
public:
	ModelRun(const ObjectiveFunc *_objectiveFunc, const Observations &_sim_obs= Observations());
	ModelRun() : obj_func_ptr(nullptr){}
	ModelRun& operator=(const ModelRun &rhs);
	virtual Parameters get_frozen_ctl_pars() const;
	virtual void set_frozen_ctl_parameters(const Parameters &frz_pars);
	virtual void add_frozen_ctl_parameters(const Parameters &frz_pars);
	virtual void clear_frozen_ctl_parameters();
	virtual void set_ctl_parameters(const Parameters &pars);
	virtual void set_observations(const Observations &obs);
	virtual void update_ctl(Parameters &ctl_pars, Observations &obs);
	virtual const Parameters &get_ctl_pars() const;
	virtual const Observations &get_obs() const;
	virtual Observations get_obs_template() const;
	const ObjectiveFunc *get_obj_func_ptr() const {return obj_func_ptr;}
	virtual double get_phi(const DynamicRegularization &dynamic_reg, double norm=2);
	virtual PhiComponets get_phi_comp(const DynamicRegularization &dynamic_reg, double norm=2) const;
	virtual vector<double> get_residuals_vec(const vector<string> &obs_names);
	//void full_report(ostream &os, const DynamicRegularization &dynamic_reg=DynamicRegularization(false),bool limit_par=false);
	virtual bool obs_valid() const;
	static  bool cmp_lt(const ModelRun &r1, const ModelRun &r2, const DynamicRegularization &reg);
	virtual ~ModelRun();
protected:
	mutable bool obs_is_valid;
	const ObjectiveFunc *obj_func_ptr;
	Parameters ctl_pars;
	Observations sim_obs;
	std::vector<std::string> frozen_ctl_par_names;
};

#endif /* MODELRUNPP_H_ */

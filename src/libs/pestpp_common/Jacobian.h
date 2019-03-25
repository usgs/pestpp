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

#ifndef JACOBIAN_H_
#define JACOBIAN_H_
#include<map>
#include<unordered_map>
#include<vector>
#include<set>
#include<list>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include "Transformable.h"
#include "Transformation.h"

class ParamTransformSeq;
class ParameterInfo;
class ParameterGroupInfo;
class RunManagerAbstract;
class ModelRun;
class FileManager;
class PriorInformation;

class JacobianRun{
public:
 JacobianRun(std::vector<double> _obs_vec = std::vector<double>(), Parameters _ctl_pars = Parameters(),
	 double _numeric_derivative_par = Parameters::no_data) : obs_vec(_obs_vec), ctl_pars(_ctl_pars),
	 numeric_derivative_par(_numeric_derivative_par){}
	std::vector<double> obs_vec;
	Parameters ctl_pars;
	double numeric_derivative_par;
};

class Jacobian {
public:
	friend void TranOffset::jacobian_forward(Jacobian &jac);
	friend void TranOffset::jacobian_reverse(Jacobian &jac);
	friend void TranScale::jacobian_forward(Jacobian &jac);
	friend void TranScale::jacobian_reverse(Jacobian &jac);
	friend void TranFixed::jacobian_forward(Jacobian &jac);
	friend void TranFixed::jacobian_reverse(Jacobian &jac);
	friend void TranLog10::jacobian_forward(Jacobian &jac);
	friend void TranLog10::jacobian_reverse(Jacobian &jac);
	friend void TranSVD::jacobian_forward(Jacobian &jac);
	friend void TranSVD::jacobian_reverse(Jacobian &jac);
	friend void TranNormalize::jacobian_forward(Jacobian &jac);
	friend void TranNormalize::jacobian_reverse(Jacobian &jac);
	Jacobian(FileManager &_file_manager);
	virtual const vector<string>& parameter_list() const{return base_numeric_par_names;}
	virtual const vector<string>& observation_list() const {return  base_sim_obs_names;}
	virtual const vector<string>& obs_and_reg_list() const;
	virtual const Parameters &get_base_numeric_parameters() const{return base_numeric_parameters;};
	Eigen::SparseMatrix<double> get_matrix(const vector<string> &obs_names, const vector<string> & par_name_vec) const;

	virtual bool build_runs(ModelRun &model_run, vector<string> numeric_par_names, ParamTransformSeq &par_transform,
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info,
		RunManagerAbstract &run_manager, set<string> &out_of_bound_par, bool phiredswh_flag=false, bool calc_init_obs=true);
	bool build_runs(Parameters &ctl_pars, Observations &ctl_obs, vector<string> numeric_par_names, ParamTransformSeq &par_transform,
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info,
		RunManagerAbstract &run_manager, set<string> &out_of_bound_par, bool phiredswh_flag=false, bool calc_init_obs=true);

	virtual void make_runs(RunManagerAbstract &run_manager);
	virtual bool process_runs(ParamTransformSeq &par_transform,
		const ParameterGroupInfo &group_info,
		RunManagerAbstract &run_manager, const PriorInformation &prior_info, bool splitswh_flag);

	virtual void save(const std::string &ext="jco") const;
	void read(const std::string &filename);
	virtual void print(std::ostream &fout) const;
	virtual const set<string>& get_failed_parameter_names() const;
	virtual long get_nonzero() const { return matrix.nonZeros();}
	virtual long get_size() const { return matrix.size(); }
	virtual int get_outersize() const { return matrix.outerSize(); }
	virtual void report_errors(std::ostream &fout);
	virtual void remove_cols(std::set<string> &rm_parameter_names);
	virtual void add_cols(set<string> &new_pars_names);
	virtual void transform(const ParamTransformSeq &par_trans, void(ParamTransformSeq::*meth_prt)(Jacobian &jac) const);
	Jacobian& operator=(const Jacobian &rhs);
	virtual const std::set<std::string>&  failed_runs_par_names(){ return  failed_parameter_names; }
	virtual ~Jacobian();


	Eigen::SparseMatrix<double> get_matrix() const{ return matrix; }
	Eigen::SparseMatrix<double>* get_matrix_ptr();
	vector<string> get_base_numeric_par_names() const{ return base_numeric_par_names;  }
	vector<string> get_sim_obs_names() const{ return base_sim_obs_names;  }

	void set_base_numeric_pars(Parameters _base_numeric_pars);
	void set_base_sim_obs(Observations _base_sim_obs);

protected:
	vector<string> base_numeric_par_names;  //ordered names of base parameters used to calculate the jacobian
	Parameters base_numeric_parameters;  //values of base parameters used to calculate the jacobian
	set<string> failed_parameter_names;
	vector< string>  base_sim_obs_names;  //names of base observations used to calculate the jacobian
	Observations  base_sim_observations;  //values of base observations used to calculate the jacobian
	//const vector<string> &ctl_file_ordered_par_names;
	//const vector<string> &ctl_file_ordered_obs_names;
	//const vector<string> &ctl_file_ordered_pi_names;
	Eigen::SparseMatrix<double> matrix;
	FileManager &file_manager;  // filemanger used to get name of jaobian file

	virtual std::vector<Eigen::Triplet<double> > calc_derivative(const string &numeric_par_name, double base_numeric_par_value, int jcol, list<JacobianRun> &run_list, const ParameterGroupInfo &group_info,
		const PriorInformation &prior_info, bool splitswh_flag);
	virtual bool forward_diff(const string &par_name, const Parameters &pest_parameters,
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, const ParamTransformSeq &par_trans,
		double &new_par, set<string> &out_of_bound_par);
	virtual bool central_diff(const string &par_name, const Parameters &pest_parameters,
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, const ParamTransformSeq &par_trans, vector<double> &new_par,
		vector<Parameters> &model_par_vec, set<string> &out_of_bound_par);
	virtual bool out_of_bounds(const Parameters &model_parameters, const ParameterInfo &ctl_par_info, set<string> &out_of_bound_par) const;
	virtual double derivative_inc(const string &name, const ParameterGroupInfo &group_info,   double cur_par_value,  bool central = false);
	virtual bool get_derivative_parameters(const string &par_name, Parameters &numeric_pars, ParamTransformSeq &par_transform, const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info,
		vector<double> &delta_numeric_par_vec, bool phiredswh_flag, set<string> &out_of_bound_par);
	virtual unordered_map<string, int> get_par2col_map() const;
	virtual unordered_map<string, int> get_obs2row_map() const;
};

#endif /* JACOBIAN_H_ */

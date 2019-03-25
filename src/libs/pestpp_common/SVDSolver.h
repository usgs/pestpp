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
#ifndef SVDSOLVER_H_
#define SVDSOLVER_H_

#include <map>
#include <set>
#include <iomanip>
#include <Eigen/Dense>
#include "Transformable.h"
#include "ParamTransformSeq.h"
#include "Jacobian.h"
#include "pest_data_structs.h"
#include "ModelRunPP.h"
#include "TerminationController.h"
#include "RunManagerAbstract.h"
#include "OutputFileWriter.h"
#include "RestartController.h"
#include "PerformanceLog.h"
#include "covariance.h"


class FileManager;
class ModelRun;
class QSqrtMatrix;
class PriorInformation;
class DynamicRegularization;
class SVDPackage;


class MuPoint
{
public:
	double mu;
	PhiComponets phi_comp;
	double target_phi_meas;
	void set(double _mu, const PhiComponets &_phi_comp);
	double f() const;
	double error_frac();
	double error_percent();
	void print(ostream &os);
	bool operator< (const MuPoint &rhs) const;
};


class SVDSolver
{
public:
	enum class MAT_INV{ Q12J, JTQJ };
protected:
	enum class LimitType {NONE, LBND, UBND, REL, FACT};
	enum class MarquardtMatrix {IDENT, JTQJ};
	enum class UpgradeBounds {ROBUST, CHEAP};
public:
	SVDSolver(Pest &_pest_scenario, FileManager &_file_manager, ObjectiveFunc *_obj_func,
		const ParamTransformSeq &_par_transform, Jacobian &_jacobian,
		OutputFileWriter &_output_file_writer, SVDSolver::MAT_INV _mat_inv,
		PerformanceLog *_performance_log, const string &description = string("base parameter solution"),Covariance parcov=Covariance(),
		bool _phiredswh_flag = false, bool _splitswh_flag = false, bool _save_next_jacobian = true);
	virtual ModelRun compute_jacobian(RunManagerAbstract &run_manager, TerminationController &termination_ctl, ModelRun &cur_run, bool restart_runs = false);
	virtual ModelRun solve(RunManagerAbstract &run_manager, TerminationController &termination_ctl, int max_iter, ModelRun &cur_run,
		ModelRun &optimum_run, RestartController &restart_controller, bool calc_first_jacobian = true);
	virtual ModelRun iteration_reuse_jac(RunManagerAbstract &run_manager, TerminationController &termination_ctl, ModelRun &base_run, bool rerun_base = true, const std::string &jco_filename = "",const std::string &res_filename="");
	virtual void iteration_jac(RunManagerAbstract &run_manager, TerminationController &termination_ctl, ModelRun &base_run, bool calc_init_obs = false, bool restart_runs = false);
	virtual ModelRun iteration_upgrd(RunManagerAbstract &run_manager, TerminationController &termination_ctl, ModelRun &base_run, bool restart_runs = false);
	virtual void set_svd_package(PestppOptions::SVD_PACK _svd_pack);
	bool get_phiredswh_flag() const { return phiredswh_flag;}
	void set_phiredswh_flag(bool _phiredswh_flag) { phiredswh_flag = _phiredswh_flag;}
	bool get_splitswh_flag() const { return splitswh_flag; }
	void set_splitswh_flag(bool _splitswh_flag) { splitswh_flag = _splitswh_flag; }
	virtual ParameterGroupInfo get_parameter_group_info() const { return *par_group_info_ptr; }
	Jacobian & get_jacobian() {return jacobian; }
	bool local_iteration_terminatated()const { return terminate_local_iteration; }
	virtual ~SVDSolver(void);
	virtual string get_solver_type() const { return svd_solver_type_name; }
protected:
	class Upgrade {
	public:
		Eigen::VectorXd uvec;
		double norm;
		vector<string> par_name_vec;
		Parameters frozen_numeric_pars;
	};

	const static string svd_solver_type_name;
	SVDPackage *svd_package;
	MAT_INV mat_inv;
	MarquardtMatrix mar_mat;
	UpgradeBounds upgrade_bounds;
	const string description;
	const ControlInfo *ctl_info;
	SVDInfo svd_info;
	ObjectiveFunc *obj_func;
	const ParameterInfo *ctl_par_info_ptr;
	const ParameterGroupInfo *par_group_info_ptr;
	ParamTransformSeq par_transform;
	const Observations *observations_ptr;
	const ObservationInfo *obs_info_ptr;
	const PriorInformation *prior_info_ptr;
	DynamicRegularization *regul_scheme_ptr;
	FileManager &file_manager;
	Jacobian &jacobian;
	bool phiredswh_flag;
	bool splitswh_flag;
	bool save_next_jacobian;
	double best_lambda;
	OutputFileWriter &output_file_writer;
	PerformanceLog *performance_log;
	std::vector<double> base_lambda_vec;
	std::vector<double> lambda_scale_vec;
	bool terminate_local_iteration;
	bool der_forgive;
	bool upgrade_augment;
	double reg_frac;
	Covariance parcov;
	double parcov_scale_fac;
	Eigen::SparseMatrix<double> JS;
	virtual void limit_parameters_ip(const Parameters &init_active_ctl_pars, Parameters &upgrade_active_ctl_pars,
		LimitType &limit_type, const Parameters &frozen_ative_ctl_pars);
	virtual Parameters limit_parameters_freeze_all_ip(const Parameters &init_active_ctl_pars,
		Parameters &upgrade_active_ctl_pars, const Parameters &frozen_active_ctl_pars = Parameters());
	virtual const string &get_description(){return description;}
	virtual void iteration_update_and_report(ostream &os, const ModelRun &base_run, ModelRun &upgrade, TerminationController &termination_ctl, RunManagerAbstract &run_manager);
	void param_change_stats(double p_old, double p_new, bool &have_fac, double &fac_change, bool &have_rel,
		double &rel_change);
	void calc_upgrade_vec(double i_lambda, Parameters &frozen_ctl_pars, QSqrtMatrix &Q_sqrt, const DynamicRegularization &regul,
		Eigen::VectorXd &residuals_vec, vector<string> &obs_names_vec, const Parameters &base_run_ctl_pars,
		Parameters &new_ctl_pars, MarquardtMatrix marquardt_type, LimitType &limit_type, bool scale_upgrade=false);
	void calc_upgrade_vec_freeze(double i_lambda, Parameters &frozen_ctl_pars, QSqrtMatrix &Q_sqrt, const DynamicRegularization &regul,
		Eigen::VectorXd &residuals_vec, vector<string> &obs_names_vec, const Parameters &base_run_ctl_pars,
		Parameters &new_ctl_pars, MarquardtMatrix marquardt_type, bool scale_upgrade = false);
	void calc_lambda_upgrade_vecQ12J(const Jacobian &jacobian, const QSqrtMatrix &Q_sqrt, const DynamicRegularization &regul,
		const Eigen::VectorXd &Residuals, const vector<string> &obs_name_vec,
		const Parameters &base_active_ctl_pars, const Parameters &freeze_active_ctl_pars,
		double lambda, Parameters &active_ctl_upgrade_pars, Parameters &upgrade_active_ctl_del_pars,
		Parameters &grad_active_ctl_del_pars, MarquardtMatrix marquardt_type, bool scale_upgrade);
	void calc_lambda_upgrade_vec_JtQJ(const Jacobian &jacobian, const QSqrtMatrix &Q_sqrt, const DynamicRegularization &regul,
		const Eigen::VectorXd &Residuals, const vector<string> &obs_name_vec,
		const Parameters &active_base_ctl_pars, const Parameters &freeze_active_ctl_pars,
		double lambda, Parameters &active_ctl_upgrade_pars, Parameters &upgrade_active_ctl_del_pars,
		Parameters &grad_active_ctl_del_pars, MarquardtMatrix marquardt_type, bool scale_upgrade=false);
	void check_limits(const Parameters &init_ctl_pars, const Parameters &upgrade_ctl_pars,
		map<string, LimitType> &limit_type_map, Parameters &active_ctl_parameters_at_limit);
	Eigen::VectorXd calc_residual_corrections(const Jacobian &jacobian, const Parameters &del_numeric_pars,
							   const vector<string> obs_name_vec);
	void dynamic_weight_adj(const ModelRun &base_run, const Jacobian &jacobian, QSqrtMatrix &Q_sqrt,
		const Eigen::VectorXd &Residuals, const vector<string> &obs_name_vec,
		const Parameters &base_active_ctl_pars, const Parameters &freeze_active_ctl_pars);
	void dynamic_weight_adj_percent(const ModelRun &base_run, double reg_frac);
	bool par_heading_out_bnd(double org_par, double new_par, double lower_bnd, double upper_bnd);
	double sidi_method(const vector<double> &x, const vector<double> &y);
	double secant_method(double x0, double y0, double x1, double y1);
	void save_frozen_pars(std::ostream &fout, const Parameters &frozen_pars, int id);
	Parameters read_frozen_pars(std::istream &fin, int id);
	PhiComponets phi_estimate(const ModelRun &base_run, const Jacobian &jacobian, QSqrtMatrix &Q_sqrt, const DynamicRegularization &regul,
		const Eigen::VectorXd &residuals_vec, const vector<string> &obs_names_vec,
		const Parameters &base_run_active_ctl_par, const Parameters &freeze_active_ctl_pars,
		DynamicRegularization &tmp_regul_scheme, bool scale_upgrade = false);
	int check_bnd_par(Parameters &new_freeze_active_ctl_pars, const Parameters &current_active_ctl_pars, const Parameters &new_upgrade_active_ctl_pars, const Parameters &new_grad_active_ctl_pars = Parameters());
};

#endif /* SVDSOLVER_H_ */

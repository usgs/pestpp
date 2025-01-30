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
#ifndef OUTPUTFILEWRITER_H
#define OUTPUTFILEWRITER_H

#include <string>
#include <iostream>
#include <fstream>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include "FileManager.h"
#include "Pest.h"
#include "ObjectiveFunc.h"

class Observations;
class ObjectiveFunc;
class Parameters;
class TranOffset;
class TranScale;
class Jacobian;
class QSqrtMatrix;
class ParameterGroupInfo;
class DynamicRegularization;

class OutputFileWriter
{
public:
	OutputFileWriter(FileManager &_file_manager, Pest &_pest_scenario, bool restart_flag = false, bool _save_rei = true, int _eigenwrite = 0);
	void prep_glm_files(bool restart_flag);
	void write_rei(std::ofstream &fout, int iter_no, const Observations &obs,
		const Observations &sim, const ObjectiveFunc &obj_func, const Parameters &pars,
        map<string,double> alt_weights=map<string,double>());
	void write_par(std::ofstream &fout, const Parameters &pars, const TranOffset &offset_tran, const TranScale &scale_tran);
	void write_restart_header(std::ostream &fout);
	void write_sen_header(std::ostream &fout, const std::string &case_name);
	void set_svd_output_opt(int _eigenwrite);
	void append_sen(std::ostream &fout, int iter_no, const Jacobian &jac, const ObjectiveFunc &obj_func,
		const ParameterGroupInfo &par_grp_info, const DynamicRegularization &regul,bool is_super,
		const ParamTransformSeq &par_transform);
	void write_svd(Eigen::VectorXd &Sigma, Eigen::SparseMatrix<double> &Vt, double lambda, const Parameters &freeze_numeric_pars, Eigen::VectorXd &Sigma_trunc);
	void write_svd_iteration(int iteration_no);

	void write_opt_constraint_rei(std::ofstream &fout, int iter_no, const Parameters pars, const Observations &obs, const Observations &sim);

	void scenario_pargroup_report(std::ostream &os);
	void scenario_io_report(std::ostream &os);
	void scenario_par_report(std::ostream &os);
	void scenario_obs_report(std::ostream &os);
	void scenario_pi_report(std::ostream &os);

	void scenario_obs_csv(ostream& os, map<string,double> alt_weights=map<string,double>());

	void phi_report(std::ostream &os,int const iter, int const nruns, PhiData const &phi_comps,
		double const dynamic_reg_weight,bool final=false, string tag="Starting");
	void par_report(std::ostream &os, Parameters const &new_ctl_pars);
	void par_report(std::ostream &os, int const iter, Parameters const &new_pars, Parameters const &old_pars, string par_type);
	void iteration_report(std::ostream &os, int iter, int nruns, string iteration_type, string svd_type=string(""));
	void scenario_report(std::ostream &os, bool report_mode=true);
	void obs_report(std::ostream &os, const Observations &obs, const Observations &sim, ObservationInfo &oi, map<string,double> alt_weights=map<string,double>());

	void param_change_stats(double p_old, double p_new, bool &have_fac, double &fac_change, bool &have_rel, double &rel_change);
	void write_par_iter(int iter, Parameters const &ctl_pars);
	void write_obj_iter(int iter, int nruns, PhiData const &pph_data);
	void write_sen_iter(int iter, map<string, double> &ctl_par_sens);

	void write_jco(bool isBaseIter, string ext, Jacobian &jco);

	void write_upgrade(int iteration, int is_super, double lambda, double scale_factor, Parameters &pars);
	void write_jco_run_id(int groupid, std::map<string, vector<int>> &par_run_map);

	void set_pest_scenario(Pest& _pest_scenario);

private:
	FileManager &file_manager;
	Pest &pest_scenario;
	std::string case_name;
	int eigenwrite;
	bool save_rei;

	void prepare_iteration_summary_files(bool restart_flag);
	void prepare_upgrade_summary_files();
	//void prepare_jco_run_id_file();

};
#endif /* OUTPUTFILEWRITER_H */

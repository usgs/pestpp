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

#include <fstream>
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include "OutputFileWriter.h"
#include "Transformable.h"
#include "Pest.h"
#include "ObjectiveFunc.h"
#include "PriorInformation.h"
#include "Jacobian.h"
#include "QSqrtMatrix.h"
#include "ModelRunPP.h"
#include "eigen_tools.h"
#include "utilities.h"

using namespace std;
using namespace Eigen;
using namespace::pest_utils;


OutputFileWriter::OutputFileWriter(FileManager &_file_manager, Pest &_pest_scenario, bool restart_flag, bool _save_rei, int _eigenwrite)
	: file_manager(_file_manager), pest_scenario(_pest_scenario),case_name(_file_manager.get_base_filename()), save_rei(_save_rei), eigenwrite(_eigenwrite)
{
	
}


void OutputFileWriter::prep_glm_files(bool restart_flag)
{
	if (restart_flag)
	{
		ofstream &fout_sen = file_manager.open_ofile_ext("sen", ofstream::app);
		write_restart_header(fout_sen);
		ofstream &fout_svd = file_manager.open_ofile_ext("svd", ofstream::app);
		write_restart_header(fout_svd);
	}
	else
	{
		ofstream &fout_sen = file_manager.open_ofile_ext("sen");
		write_sen_header(fout_sen, case_name);
		file_manager.open_ofile_ext("svd");
	}
	if (pest_scenario.get_pestpp_options().get_iter_summary_flag())
	{
		prepare_iteration_summary_files(restart_flag);
		prepare_upgrade_summary_files();
	}
}

void OutputFileWriter::iteration_report(std::ostream &os, int iter, int nruns, string iteration_type, string svd_type)
{
	os << "OPTIMISATION ITERATION NUMBER: " << iter << endl << endl;
	os << "  Iteration type: " << iteration_type << endl;
	if (!svd_type.empty())
	{
		os << "  SVD Package: " << svd_type << endl;
	}
	os << "  Model calls so far : " << nruns << endl;
	os << endl;
}

//void OutputFileWriter::prepare_jco_run_id_file()
//{
//	file_manager.open_ofile_ext("irid.csv");
//	ofstream &os = file_manager.get_ofstream("irid.csv");
//	os << "group_id";
//	for (auto &pname : pest_scenario.get_ctl_ordered_par_names())
//	{
//		os << ',' << pname;
//	}
//	os << endl;
//}

void OutputFileWriter::write_jco_run_id(int group_id, std::map<string, vector<int>> &par_run_map)
{
	if (!pest_scenario.get_pestpp_options().get_iter_summary_flag())
		return;

	file_manager.open_ofile_ext("rid");
	ofstream &os = file_manager.get_ofstream("rid");
	os << endl << "group id: " << group_id << endl;
	for (auto &pname : pest_scenario.get_ctl_ordered_par_names())
	{
		if (!par_run_map.count(pname))
			continue;
		os << pname;
		for (auto &id : par_run_map[pname])
			os << ',' << id;
		os << endl;
	}
	//file_manager.close_file("rid");
}

void OutputFileWriter::prepare_upgrade_summary_files()
{
	file_manager.open_ofile_ext("upg.csv");
	ofstream &os_up = file_manager.get_ofstream("upg.csv");
	os_up << "iteration,_is_super,_lambda,_scale_factor";
	for (auto &par_name : pest_scenario.get_ctl_ordered_par_names())
	{
		os_up << ',' << lower_cp(par_name);
	}
	os_up << endl;
}

void OutputFileWriter::write_upgrade(int iteration, int is_super, double lambda, double scale_factor, Parameters &pars)
{
	ofstream &os_up = file_manager.get_ofstream("upg.csv");
	os_up << iteration << ',' << is_super << ',' << lambda << ',' << scale_factor;
	//vector<double> vals = pars.get_data_vec(pest_scenario.get_ctl_ordered_par_names());
	//for (auto &pval : vals)
	double pval;
	vector<string> keys = pars.get_keys();
	for (auto &name : pest_scenario.get_ctl_ordered_par_names())
	{
		os_up << ',';
		if (find(keys.begin(), keys.end(), name) != keys.end())
		{
			pval = pars.get_rec(name);
			os_up << pval;
		}
	}
	os_up << endl;

}

void OutputFileWriter::prepare_iteration_summary_files(bool restart_flag)
{
	if (restart_flag)
	{
		file_manager.open_ofile_ext("ipar", ofstream::app);
		file_manager.open_ofile_ext("iobj", ofstream::app);
		file_manager.open_ofile_ext("isen", ofstream::app);
	}
	else
	{
		file_manager.open_ofile_ext("ipar");
		file_manager.open_ofile_ext("iobj");
		file_manager.open_ofile_ext("isen");
		ofstream &os_ipar = file_manager.get_ofstream("ipar");
		ofstream &os_isen = file_manager.get_ofstream("isen");
		os_ipar << "iteration";
		os_isen << "iteration";
		for (auto &par_name : pest_scenario.get_ctl_ordered_par_names())
		{
			os_ipar << ',' << lower_cp(par_name);
			os_isen << ',' << lower_cp(par_name);
		}
		os_ipar << endl;
		os_isen << endl;
		ofstream &os_iobj = file_manager.get_ofstream("iobj");
		os_iobj << "iteration,model_runs_completed,total_phi,measurement_phi,regularization_phi";
		for (auto &obs_grp : pest_scenario.get_ctl_ordered_obs_group_names())
		{
			os_iobj << ',' << lower_cp(obs_grp);
		}
		os_iobj << endl;
	}
}

void OutputFileWriter::write_sen_iter(int iter, map<string, double> &ctl_par_sens)
{
	if (!pest_scenario.get_pestpp_options().get_iter_summary_flag())
		return;
	ofstream &os = file_manager.get_ofstream("isen");
	os << iter;
	map<string, double>::iterator is;
	for (auto &par_name : pest_scenario.get_ctl_ordered_par_names())
	{
		is = ctl_par_sens.find(par_name);
		if (is == ctl_par_sens.end())
			os << ',' << -999;
		else
			os << ',' << ctl_par_sens.at(par_name);
	}
	os << endl;
}

void OutputFileWriter::write_par_iter(int iter, Parameters const &ctl_pars)
{
	if (!pest_scenario.get_pestpp_options().get_iter_summary_flag())
		return;
	ofstream &os = file_manager.get_ofstream("ipar");
	os << iter;
	for (auto &par_name : pest_scenario.get_ctl_ordered_par_names())
	{
		os << ',' << ctl_pars.get_rec(par_name);
	}
	os << endl;
}

void OutputFileWriter::write_obj_iter(int iter, int nruns, PhiData const &phi_data)
{
	if (!pest_scenario.get_pestpp_options().get_iter_summary_flag())
		return;
	ofstream &os = file_manager.get_ofstream("iobj");
	os << iter << ',' << nruns;
	os << ',' << phi_data.total();
	os << ',' << phi_data.meas;
	os << ',' << phi_data.regul;
	for (auto &obs_grp : pest_scenario.get_ctl_ordered_obs_group_names())
	{
		os << ',' << phi_data.group_phi.at(obs_grp);
	}
	os << endl;
}

void OutputFileWriter::scenario_report(std::ostream &os, bool report_mode)
{
	string mode;// = "estimation";
	/*if (pest_scenario.get_regul_scheme_ptr()->get_use_dynamic_reg())
	{
		mode = "regularization (with a \"z\")";
	}*/

	os << endl << "    This software has been approved for release by the" << endl;
	os << "    U.S.Geological Survey(USGS).Although the software has " << endl;
	os << "    been subjected to rigorous review, the USGS reserves the" << endl;
	os << "    right to update the software as needed pursuant to further" << endl;
	os << "    analysisand review.No warranty, expressed or implied, is " << endl;
	os << "    made by the USGS or the U.S.Government as to the" << endl;
	os << "    functionality of the softwareand related material nor shall" << endl;
	os << "    the fact of release constitute any such warranty." << endl;
	os << "    Furthermore, the software is released on condition that" << endl;
	os << "    neither the USGS nor the U.S.Government shall be held" << endl;
	os << "    liable for any damages resulting from its authorized " << endl;
	os << "    or unauthorized use." << endl << endl;

	switch (pest_scenario.get_control_info().pestmode)
	{
	case ControlInfo::PestMode::ESTIMATION: 
	{
		mode = "estimation";
		break;
	}

	case ControlInfo::PestMode::REGUL: 
	{
		mode = "regularization";
		break;
	}

	case ControlInfo::PestMode::PARETO:
	{
		mode = "pareto";
		break;
	}
	case ControlInfo::PestMode::UNKNOWN: 
	{
		mode = "unknown";
		break;
	}


	default: mode = "WTF";
	}
	//os << endl << endl << "binary compiled on " << __DATE__ << " at " << __TIME__ << endl << endl;
	if (report_mode)
		os << endl << "pestmode:- " << endl << "   " << mode << endl << endl;

	os << endl << "Case dimensions:- " << endl;
	os << setw(0) << "    Number of parameters = " << pest_scenario.get_ctl_ordered_par_names().size() << endl;
	os << setw(0) << "    Number of adjustable parameters = " << pest_scenario.get_n_adj_par() << endl;
	os << setw(0) << "    Number of observations = " << pest_scenario.get_ctl_ordered_obs_names().size() << endl;
    os << setw(0) << "    Number of non-zero weighted observations = " << pest_scenario.get_ctl_ordered_nz_obs_names().size() << endl;
    os << setw(0) << "    Number of prior estimates = " << pest_scenario.get_ctl_ordered_pi_names().size() << endl << endl;

	os << pest_scenario.get_control_info() << endl;
	pest_scenario.get_pestpp_options().summary(os);

	scenario_io_report(os);
	scenario_pargroup_report(os);
	scenario_par_report(os);
	scenario_obs_report(os);
	scenario_pi_report(os);

	os << endl << pest_scenario.get_svd_info() << endl;
	os << endl;

	if (((report_mode) && (pest_scenario.get_control_info().pestmode == ControlInfo::PestMode::REGUL)) &&
		(pest_scenario.get_regul_scheme_ptr()))
	{
		os << "Regularization information:" << endl;
		os << setw(0) << "    phimlim = " << pest_scenario.get_regul_scheme_ptr()->get_phimlim() << endl;
		os << setw(0) << "    fracphim = " << pest_scenario.get_regul_scheme_ptr()->get_fracphim() << endl;
		os << setw(0) << "    phimaccept = " << pest_scenario.get_regul_scheme_ptr()->get_phimaccept() << endl;
		os << setw(0) << "    wfinit = " << pest_scenario.get_regul_scheme_ptr()->get_wfinit() << endl;

	}
	/*os << "PEST++ arguments:" << endl;
	for (auto pp_arg : pest_scenario.get_pestpp_options().get_arg_map())
		os << setw(0) << pp_arg.first << " = " << pp_arg.second << endl;*/


	os << endl;
	os << endl;
	os << endl;

}
void OutputFileWriter::set_pest_scenario(Pest& _pest_scenario) { pest_scenario = _pest_scenario; }

void OutputFileWriter::scenario_io_report(std::ostream &os)
{

	os << "Model command line(s):- " << endl;
	for (auto &cmd : pest_scenario.get_comline_vec())
	{
		os << "    " << cmd << endl;
	}
	os << endl;

	os << "Model interface files:-" << endl;
	os << "    template files:" << endl;
	for (auto &f : pest_scenario.get_tplfile_vec())
	{
		os << "      " << f << endl;
	}
	os << "    model input files:" << endl;
	for (auto &f : pest_scenario.get_inpfile_vec())
	{
		os << "      " << f << endl;
	}
	os << endl << "    instruction files:" << endl;
	for (auto &f : pest_scenario.get_insfile_vec())
	{
		os << "      " << f << endl;
	}
	os << "    model output files:" << endl;
	for (auto &f : pest_scenario.get_outfile_vec())
	{
		os << "      " << f << endl;
	}
	os << endl << endl;
}

void OutputFileWriter::scenario_pargroup_report(std::ostream &os)
{
	const ParameterGroupRec *grp_rec;
	int grp_len = 12;
	for (auto& par_name : pest_scenario.get_ctl_ordered_par_group_names())
		grp_len = max((int)par_name.size(), grp_len);
	grp_len++;
	os << "Parameter group information" << endl;
	os << left << setw(grp_len) << "NAME" << right << setw(15) << "INCREMENT TYPE" << setw(25) << "DERIVATIVE INCREMENT";
	os << setw(25) << "INCREMENT LOWER BOUND" << setw(15) << "FORCE CENTRAL" << setw(25) << "INCREMENT MULTIPLIER" << endl;
	for (auto &grp_name : pest_scenario.get_ctl_ordered_par_group_names())
	{
		grp_rec = pest_scenario.get_base_group_info().get_group_by_groupname(grp_name);
		os << left << setw(grp_len) << lower_cp(grp_rec->name) << right << setw(15) << grp_rec->inctyp << setw(25) << grp_rec->derinc;
		os << setw(25) << grp_rec->derinclb << setw(15) << grp_rec->forcen << setw(25) << grp_rec->derincmul << endl;
	}
	os << endl << endl;
}

void OutputFileWriter::scenario_par_report(std::ostream &os)
{
	if (pest_scenario.get_ctl_ordered_par_names().size() > 100000)
	{
		os << endl << "...more than 100,000 pars, not writing par data" << endl;
		return;
	}
	map<int, string> trans_type;
	trans_type[0] = "none";
	trans_type[1] = "fixed";
	trans_type[2] = "tied";
	trans_type[3] = "log";
	int grp_len = 12;
	for (auto& par_name : pest_scenario.get_ctl_ordered_par_group_names())
		grp_len = max((int)par_name.size(), grp_len);
	grp_len++;
	int par_len = 12;
	for (auto& par_name : pest_scenario.get_ctl_ordered_par_names())
		par_len = max((int)par_name.size(), par_len);
	par_len++;
	os << endl << "Parameter information" << endl;
	os << left << setw(par_len) << "NAME" << setw(10) << "TRANSFORMATION" << right << setw(20) << "CHANGE LIMIT" << setw(15) << "INITIAL VALUE";
	os << setw(15) << "LOWER BOUND";
	os << setw(15) << "UPPER BOUND" << setw(grp_len) << "GROUP";

	os << setw(15) << "SCALE" << setw(15) << "OFFSET" << setw(20) << "DERIVATIVE COMMAND" << endl;
	const ParameterRec* par_rec;
	for (auto &par_name : pest_scenario.get_ctl_ordered_par_names())
	{
		par_rec = pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(par_name);
		os << left << setw(par_len) << lower_cp(par_name);
		os << setw(10) << trans_type[static_cast<int>(par_rec->tranform_type)];
		os << right << setw(20) << par_rec->chglim;
		os << setw(15) << par_rec->init_value;
		os << setw(15) << par_rec->lbnd;
		os << setw(15) << par_rec->ubnd;
		os << setw(grp_len) << lower_cp(par_rec->group);
		os << setw(15) << par_rec->scale;
		os << setw(15) << par_rec->offset;
		os << setw(20) << par_rec->dercom << endl;
	}
	os << endl << endl;
}

void OutputFileWriter::scenario_obs_csv(ostream& os, map<string,double> alt_weights)
{
	if (os.bad())
		throw runtime_error("OutputFileWriter::scenario_obs_csv(): os is bad");
	os << "name,value,group,weight" << endl;
	const ObservationRec* obs_rec;
	const Observations& obs = pest_scenario.get_ctl_observations();
    double weight;
	for (auto& obs_name : pest_scenario.get_ctl_ordered_obs_names())
	{
		obs_rec = pest_scenario.get_ctl_observation_info().get_observation_rec_ptr(obs_name);
        weight = obs_rec->weight;
        if (alt_weights.find(obs_name) != alt_weights.end())
        {
            weight = alt_weights[obs_name];
        }
		os << lower_cp(obs_name) << "," << obs.get_rec(obs_name) << "," << lower_cp(obs_rec->group) << "," << weight << endl;
	}
}

void OutputFileWriter::scenario_obs_report(std::ostream &os)
{
	if (pest_scenario.get_ctl_ordered_obs_names().size() > 100000)
	{
		os << endl << "...more than 100,000 obs, not writing obs data" << endl;
		return;
	}
	int obs_len = 20;
	for (auto& obs_name : pest_scenario.get_ctl_ordered_obs_names())
		obs_len = max((int)obs_name.size(), obs_len);
	obs_len++;
	int grp_len = 20;
	for (auto& obs_name : pest_scenario.get_ctl_ordered_obs_group_names())
		grp_len = max((int)obs_name.size(), grp_len);
	grp_len++;
	os << endl << "Observation information" << endl;
	os << left << setw(obs_len) << "NAME" << right << setw(20) << "VALUE" << setw(grp_len) << "GROUP" << setw(20) << "WEIGHT" << endl;
	const ObservationRec* obs_rec;
	const Observations &obs = pest_scenario.get_ctl_observations();
	for (auto &obs_name : pest_scenario.get_ctl_ordered_obs_names())
	{
		obs_rec = pest_scenario.get_ctl_observation_info().get_observation_rec_ptr(obs_name);
		os << left << setw(obs_len) << lower_cp(obs_name);
		os << right << setw(20) << obs.get_rec(obs_name);
		os << setw(grp_len) << lower_cp(obs_rec->group);
		os << setw(20) << obs_rec->weight << endl;
	}
	os << endl << endl;
}


void OutputFileWriter::scenario_pi_report(std::ostream &os)
{
	const PriorInformation &pi = pest_scenario.get_prior_info();
	os << endl << "Prior information" << endl;
	if (pi.size() == 0)
	{
		os << endl << "   no prior information provided" << endl;
	}
	for (auto &pi_name : pi)
	{
		os << lower_cp(pi_name.first) << "  " << pi_name.second;
	}
	os << endl << endl;
}


void OutputFileWriter::par_report(std::ostream &os, Parameters const &new_ctl_pars)
{
	double val;
	os << endl;
	os << "     Parameter            " << endl;
	os << "        Name         Value" << endl;
	os << "    ------------  ------------" << endl;
	for (auto &p_name : pest_scenario.get_ctl_ordered_par_names())
	{
		val = new_ctl_pars.get_rec(p_name);
		os << left;
		os << "    " << setw(12) << lower_cp(p_name);
		os << right;
		os << "  " << setw(12) << val << endl;
	}
	os << endl;
}


void OutputFileWriter::par_report(std::ostream &os, int const iter, Parameters const &new_pars, Parameters const &old_pars,
	string par_type)
{
	double p_old, p_new;
	double fac_change = -9999, rel_change = -9999;
	bool have_fac = false, have_rel = false;
	double max_fac_change = 0;
	double max_rel_change = 0;
	string max_fac_par = "N/A";
	string max_rel_par = "N/A";

	os << "    Iteration "<<iter<<" Parameter Upgrades (" << par_type << " Parameters) " << endl;
	os << "      Parameter     Current       Previous       Factor       Relative" << endl;
	os << "        Name         Value         Value         Change        Change" << endl;
	os << "      ----------  ------------  ------------  ------------  ------------" << endl;
	vector<string> par_names;
	if (lower_cp(par_type) == "control file")
		par_names = pest_scenario.get_ctl_ordered_par_names();
	else
	{
		par_names = new_pars.get_keys();
		sort(par_names.begin(), par_names.end());
	}
	//for (const auto &ipar : new_ctl_pars)
	for (auto &p_name : par_names)
	{
		Parameters::const_iterator pi = new_pars.find(p_name);
		if (pi == new_pars.end()) continue;
		p_new = new_pars.get_rec(p_name);
		p_old = old_pars.get_rec(p_name);
		param_change_stats(p_old, p_new, have_fac, fac_change, have_rel, rel_change);
		if (have_fac && fac_change >= max_fac_change)
		{
			max_fac_change = fac_change;
			max_fac_par = p_name;
		}
		if (have_rel && abs(rel_change) >= abs(max_rel_change))
		{
			max_rel_change = rel_change;
			max_rel_par = p_name;
		}
		os << right;
		os << "    " << setw(12) << lower_cp(p_name);
		os << right;
		os << "  " << setw(12) << p_new;
		os << "  " << setw(12) << p_old;
		if (have_fac)
			os << "  " << setw(12) << fac_change;
		else
			os << "  " << setw(12) << "N/A";
		if (have_rel)
			os << "  " << setw(12) << rel_change;
		else
			os << "  " << setw(12) << "N/A";
		os << endl;
	}
	os << "       Maximum changes in \"" << par_type << "\" parameters:" << endl;
	os << "         Maximum relative change = " << max_rel_change << "   [" << lower_cp(max_rel_par) << "]" << endl;
	os << "         Maximum factor change = " << max_fac_change << "   [" << lower_cp(max_fac_par) << "]" << endl;
	os << endl;
	if ((lower_cp(par_type) == "control file") && (pest_scenario.get_pestpp_options().get_iter_summary_flag()))
	{
		write_par_iter(iter, new_pars);
	}
}


void OutputFileWriter::param_change_stats(double p_old, double p_new, bool &have_fac,
	double &fac_change, bool &have_rel, double &rel_change)
{
	have_rel = have_fac = true;
	double a = max(abs(p_new), abs(p_old));
	double b = min(abs(p_new), abs(p_old));
	// compute relative change
	if (p_old == 0) {
		have_rel = false;
		rel_change = -9999;
	}
	else
	{
		rel_change = (p_old - p_new) / p_old;
	}
	//compute factor change
	if (p_old == 0.0 || p_new == 0.0) {
		have_fac = false;
		fac_change = -9999;
	}
	else {
		fac_change = a / b;
	}
}




void OutputFileWriter::phi_report(std::ostream &os, int const iter, int const nruns, PhiData const &phi_comps, double const dynamic_reg_weight,bool final,string tag)
{
	map<string, double>::const_iterator it = phi_comps.group_phi.find("REGUL");
	if ((!dynamic_reg_weight) || (it == phi_comps.group_phi.end()))
	{
		if (final)
		{
			os << "  Final phi                                           Total : " << phi_comps.total() << endl;
		}
		else
		{
			os << endl << "  " << tag << " phi for this iteration                     Total : " << phi_comps.total() << endl << endl << endl;
		}
	}
	else
	{
		if (final)
		{
			os << "  Final regularization weight factor                        : " << dynamic_reg_weight << endl;
			os << "  Final phi                                           Total : " << phi_comps.total() << endl;
			os << "  Final measurement phi for this iteration            Total : " << phi_comps.meas << endl;
			os << "  Final regularization phi for this iteration         Total : " << phi_comps.regul << endl;
		}
		else
		{
			os << "  Current regularization weight factor                      : " << dynamic_reg_weight << endl;
			os << "  " << tag << " phi for this iteration                     Total : " << phi_comps.total() << endl;
			os << "  " << tag << " measurement phi for this iteration         Total : " << phi_comps.meas << endl;
			os << "  " << tag << " regularization phi for this iteration      Total : " << phi_comps.regul << endl;
		}
	}

	for (auto &gname : pest_scenario.get_ctl_ordered_obs_group_names())
	{
		os << "  Contribution to phi from observation group ";
		os << setw(17) << setiosflags(ios::right) << "\"" + lower_cp(gname) + "\" : ";
		os << phi_comps.group_phi.at(gname) << endl;
	}
	/*if (pest_scenario.get_pestpp_options().get_iter_summary_flag())
	{
		write_obj_iter(iter, nruns,phi_comps);
	}*/
}


void OutputFileWriter::obs_report(ostream &os, const Observations &obs, const Observations &sim, ObservationInfo &oi, map<string,double> alt_weights)
{
	vector<string> obs_name_vec = pest_scenario.get_ctl_ordered_obs_names();
	int nsize = 20;
	for (auto o : obs_name_vec)
	{
		nsize = max(nsize, int(o.size()));
	}
	os << setw(nsize + 1) << " Name" << setw(13) << " Group" << setw(21) << " Measured" << setw(21) << " Modelled" << setw(21) << " Residual" << setw(21) << " Weight" << endl;
	//vector<string> obs_name_vec = obs.get_keys();
	
	double obs_val, sim_val;
	//ObservationInfo oi = pest_scenario.get_ctl_observation_info();
	//for(vector<string>::const_iterator b = obs_name_vec.begin(),
	//	e = obs_name_vec.end(); b!=e; ++b)
    double weight;
	if (obs_name_vec.size() < 100000)
	{
		for (auto& b : obs_name_vec)
		{
            weight = oi.get_observation_rec_ptr(b)->weight;
            if (alt_weights.find(b) != alt_weights.end())
            {
                weight = alt_weights.at(b);
            }
			obs_val = obs.get_rec(b);
			sim_val = sim.get_rec(b);
			os << " " << setw(nsize) << lower_cp(b)
				<< " " << setw(12) << lower_cp(oi.get_observation_rec_ptr(b)->group)
				<< " " << showpoint << setw(20) << obs_val
				<< " " << showpoint << setw(20) << sim_val
				<< " " << showpoint << setw(20) << obs_val - sim_val
				<< " " << showpoint << setw(20) << weight << endl;
		}
	}
	else
	{
		for (auto& b : obs_name_vec)
		{
            weight = oi.get_observation_rec_ptr(b)->weight;
            if (alt_weights.find(b) != alt_weights.end())
            {
                weight = alt_weights.at(b);
            }
			obs_val = obs.get_rec(b);
			sim_val = sim.get_rec(b);
			os << " " << lower_cp(b) 
				<< " " << lower_cp(oi.get_observation_rec_ptr(b)->group)
				<< " " << obs_val
				<< " " << sim_val
				<< " " << obs_val - sim_val
				<< " " << weight << endl;
		}
	}


}


void OutputFileWriter::write_opt_constraint_rei(std::ofstream &fout, int iter_no, const Parameters pars, const Observations &obs, const Observations &sim)
{
	fout << setiosflags(ios::left);
	fout.unsetf(ios::floatfield);
	fout.precision(12);
	fout << " MODEL OUTPUTS AT END OF OPTIMISATION ITERATION NO. " << iter_no << ":-" << endl;
	fout << endl << endl;
	ObservationInfo oi = pest_scenario.get_ctl_observation_info();
	obs_report(fout, obs, sim, oi);
	//process prior information
	//const PriorInformation *prior_info_ptr = obj_func.get_prior_info_ptr();
	const PriorInformation *prior_info_ptr = pest_scenario.get_prior_info_ptr();
	const PriorInformationRec *pi_rec_ptr;
	PriorInformation::const_iterator ipi;
	//for(PriorInformation::const_iterator b = prior_info_ptr->begin(),
	//	e = prior_info_ptr->end(); b!=e; ++b)
	vector<string> obs_name_vec = pest_scenario.get_ctl_ordered_pi_names();
	double obs_val, residual, sim_val;

	int nsize = 20;
	for (auto& b : obs_name_vec)
	{
		nsize = max(nsize, int(b.size()));
	}

	for (auto &b : obs_name_vec)
	{
		ipi = prior_info_ptr->find(b);
		pi_rec_ptr = &(*ipi).second;
		obs_val = pi_rec_ptr->get_obs_value();
		residual = pi_rec_ptr->calc_residual(pars);
		sim_val = obs_val + residual;
		fout << " " << setw(nsize) << lower_cp(b)
			<< " " << setw(12) << lower_cp(pi_rec_ptr->get_group())
			<< " " << showpoint << setw(20) << obs_val
			<< " " << showpoint << setw(20) << sim_val
			<< " " << showpoint << setw(20) << residual
			<< " " << showpoint << setw(20) << sqrt(pi_rec_ptr->get_weight()) << endl;
	}
}


void OutputFileWriter::write_rei(ofstream &fout, int iter_no, const Observations &obs, const Observations &sim,
	const ObjectiveFunc &obj_func, const Parameters &pars, map<string,double> alt_weights)
{
	fout << setiosflags(ios::left);
	fout.unsetf(ios::floatfield);
	fout.precision(12);
	fout << " MODEL OUTPUTS AT END OF OPTIMISATION ITERATION NO. " << iter_no << ":-" << endl;
	fout << endl << endl;
	ObservationInfo oi = pest_scenario.get_ctl_observation_info();
	obs_report(fout, obs, sim, oi,alt_weights);
	//process prior information
	const PriorInformation *prior_info_ptr = pest_scenario.get_prior_info_ptr();
	const PriorInformationRec *pi_rec_ptr;
	PriorInformation::const_iterator ipi;
	//for(PriorInformation::const_iterator b = prior_info_ptr->begin(),
	//	e = prior_info_ptr->end(); b!=e; ++b)
	vector<string> obs_name_vec = pest_scenario.get_ctl_ordered_pi_names();
	double obs_val, residual, sim_val;
	for (auto &b : obs_name_vec)
	{
		ipi = prior_info_ptr->find(b);
		pi_rec_ptr = &(*ipi).second;
		obs_val = pi_rec_ptr->get_obs_value();
		residual = pi_rec_ptr->calc_residual(pars);
		sim_val = obs_val + residual;
		fout << " " << setw(20) << lower_cp(b)
		<< " " << setw(12) << lower_cp(pi_rec_ptr->get_group())
		<< " " << showpoint <<   setw(20) << obs_val
		<< " " << showpoint <<  setw(20) << sim_val
		<<  " " << showpoint <<  setw(20) << residual
		<< " " << showpoint << setw(20)  << pi_rec_ptr->get_weight() << endl;
	}
}



void OutputFileWriter::write_par(ofstream &fout, const Parameters &pars, const TranOffset &offset_tran, const TranScale &scale_tran)
{
	Parameters::const_iterator it;
	pair<bool, double> val_pair;
	double scale, offset;

	fout.unsetf(ios::floatfield);
	fout.precision(numeric_limits<double>::digits10 + 1);
	fout << "single point" << endl;
	//for(Parameters::const_iterator b=pars.begin(), e=pars.end();
	//	b!=e; ++b)
	vector<string> par_name_vec = pest_scenario.get_ctl_ordered_par_names();
	for (auto &b : par_name_vec)
	{
		//name_ptr = &(*b).first;
		val_pair = offset_tran.get_value(b);
		if (val_pair.first == true)
			offset = val_pair.second;
		else
			offset = 0.0;
		val_pair = scale_tran.get_value(b);
		if (val_pair.first == true)
			scale = val_pair.second;
		else
			scale = 1.0;

		fout << setw(14) << lower_cp(b) << setw(22) << " "
		<<  showpoint<< pars.get_rec(b) << " " << setw(20) << showpoint << scale << " " << setw(20) << showpoint << offset << endl;
	}
	fout.flush();
}

void OutputFileWriter::write_sen_header(std::ostream &fout, const string &case_name)
{
	fout << "                    PARAMETER SENSITIVITIES: CASE " << case_name << endl;
	fout << endl << endl;
}


void OutputFileWriter::write_restart_header(std::ostream &fout)
{
	fout << endl << endl;
	fout << "Restarting PEST++ .... " << endl;
	fout << endl << endl;
}


void OutputFileWriter::append_sen(std::ostream &fout, int iter_no, const Jacobian &jac,
	const ObjectiveFunc &obj_func, const ParameterGroupInfo &par_grp_info, const DynamicRegularization &regul,
	bool is_super, const ParamTransformSeq &par_transform)
{
	fout << setiosflags(ios::left);
	fout.unsetf(ios::floatfield);
	fout.precision(12);
	fout << "NUMERIC PARAMETER SENSITIVITIES FOR OPTIMISATION ITERATION NO. " << setw(3) << iter_no << " ----->" << endl;
	//fout << " Parameter name   Group        Current Value           CSS w/reg           CSS w/o reg" << endl;
	  //fout << "  Parameter name   Group        Current Value        HILL_CSS_w_reg       PEST_CSS_w_reg       HILL_CSS_wo_reg      PEST_CSS_wo_reg" << endl;
	fout << "  Parameter name   Group        Current Value        PEST_CSS_wo_reg" << endl;
	const vector<string> &par_list = jac.parameter_list();
	//const vector<string> &par_list = pest_scenario.get_ctl_ordered_par_names();
	const vector<string> &obs_list = jac.obs_and_reg_list();
	Parameters pars = jac.get_base_numeric_parameters();
	Parameters ctl_pars = par_transform.numeric2ctl_cp(pars);
	if (pars.size() == 0)
	{
		fout << "parameter values are not avaialble to compute CSS" << endl;
		fout << endl << endl;
	}
	else
	{
		VectorXd par_vec = pars.get_data_eigen_vec(par_list);
		MatrixXd par_mat_tmp = par_vec.asDiagonal();
		Eigen::SparseMatrix<double> par_mat = par_mat_tmp.sparseView();
		QSqrtMatrix Q_sqrt(obj_func.get_obs_info_ptr(), obj_func.get_prior_info_ptr());
		//Eigen::SparseMatrix<double> q_sqrt_reg = Q_sqrt.get_sparse_matrix(obs_list, regul);
		//Eigen::SparseMatrix<double> dss_mat_reg = q_sqrt_reg * jac.get_matrix(obs_list, par_list) * par_mat;
		Eigen::SparseMatrix<double> q_sqrt_no_reg = Q_sqrt.get_sparse_matrix(obs_list, DynamicRegularization::get_zero_reg_instance());
		//Eigen::SparseMatrix<double> dss_mat_no_reg = q_sqrt_no_reg * jac.get_matrix(obs_list, par_list) * par_mat;
		//Eigen::SparseMatrix<double> dss_mat_reg_pest = q_sqrt_reg * jac.get_matrix(obs_list, par_list);
		Eigen::SparseMatrix<double> dss_mat_no_reg_pest = q_sqrt_no_reg * jac.get_matrix(obs_list, par_list);

		//cout << q_sqrt_reg << endl << endl;
		//cout << q_sqrt_no_reg << endl << endl;
	    //cout << q_sqrt_reg << endl << endl;
		//cout << dss_mat_reg << endl << endl;
		//cout << q_sqrt_no_reg << endl << endl;
		//cout << dss_mat_no_reg << endl << endl;
		//cout << jac.get_matrix(obs_list, par_list) << endl;

		int n_par = par_list.size();
		//int n_nonzero_weights_reg = q_sqrt_reg.nonZeros();
		//int n_nonzero_weights_no_reg = q_sqrt_no_reg.nonZeros();

		int n_nonzero_weights_no_reg = obj_func.get_obs_info_ptr()->get_nnz_obs();
		int n_nonzero_weights_reg = n_nonzero_weights_no_reg + obj_func.get_prior_info_ptr()->get_nnz_pi();
		vector<string> par_names;
		if (is_super)
		{
			par_names = par_list;
			sort(par_names.begin(), par_names.end());
		}
		else
		{
			par_names = pest_scenario.get_ctl_ordered_par_names();
		}
		//for isen file
		map<string, double> par_sens;
		double val;
		//drop any names that aren't in par_list
		vector<string>::const_iterator is;
		int i;
		//for (int i = 0; i < n_par; ++i)
		for (auto &pname : par_names)
		{
			is = find(par_list.begin(), par_list.end(), pname);
			if (is == par_list.end()) continue;
			i = is - par_list.begin();
			fout << "   " << setw(15) << lower_cp(pname)
				<< " " << setw(12) << lower_cp(par_grp_info.get_group_name(pname));
			if (is_super)
				fout << " " << showpoint << setw(20) << pars.get_rec(pname);
			else
				fout << " " << showpoint << setw(20) << ctl_pars.get_rec(pname);

			/*if (n_nonzero_weights_reg > 0)
			{
				fout << " " << showpoint << setw(20) << dss_mat_reg.col(i).norm() / pow(n_nonzero_weights_reg, 2.0);
				fout << " " << showpoint << setw(20) << dss_mat_reg_pest.col(i).norm() / n_nonzero_weights_reg;
			}
			else
			{
				fout << " " << showpoint << setw(20) << "NA";
				fout << " " << showpoint << setw(20) << "NA";
			}*/
			if (n_nonzero_weights_no_reg > 0)
			{
				//val = dss_mat_no_reg.col(i).norm() / pow(n_nonzero_weights_no_reg, 2.0);
				//fout << " " << showpoint << setw(20) << val;
				val = dss_mat_no_reg_pest.col(i).norm() / n_nonzero_weights_no_reg;
				par_sens[pname] = val;
				fout << " " << showpoint << setw(20) << val;
			}
			else
			{
				fout << " " << showpoint << setw(20) << "NA";
				fout << " " << showpoint << setw(20) << "NA";
				par_sens[pname] = -999;
			}
			fout << endl;
		}
		fout << endl << endl;
		if (!is_super)
		{
			write_sen_iter(iter_no, par_sens);
		}
	}
}


void OutputFileWriter::set_svd_output_opt(int _eigenwrite)
{
	eigenwrite = _eigenwrite;
}


void OutputFileWriter::write_svd(VectorXd &Sigma, Eigen::SparseMatrix<double> &Vt, double lambda, const Parameters &freeze_numeric_pars, VectorXd &Sigma_trunc)
{
		ofstream &fout_svd = file_manager.get_ofstream("svd");
		fout_svd<< "CURRENT VALUE OF MARQUARDT LAMBDA = " << lambda << " --------->" << endl << endl;
		fout_svd << "FROZEN PARAMETERS-" << endl;
		fout_svd << freeze_numeric_pars << endl << endl;
		fout_svd << "SINGULAR VALUES IN SOLUTION:-" << endl;
		print(Sigma, fout_svd, 7);
		fout_svd << "TRUNCATED SINGULAR VALUES:-" << endl;
		print(Sigma_trunc, fout_svd, 7);
		fout_svd << endl << endl;
		if (eigenwrite > 0)
		{
			fout_svd << "MATRIX OF EIGENVECTORS IN SOLUTION:-" << endl;
			print(Vt.toDense(), fout_svd, 7);
			fout_svd << endl;
		}
		fout_svd << "Number of singular values used in solution = " << Sigma.size() << endl << endl << endl;
}


void OutputFileWriter::write_svd_iteration(int iteration_no)
{
		ofstream &fout_svd = file_manager.get_ofstream("svd");
	// write head for SVD file
		fout_svd << "------------------------------------------------------------------------------" << endl;
		fout_svd << "OPTIMISATION ITERATION NO.        : " << iteration_no << endl << endl;
}

void OutputFileWriter::write_jco(bool isBaseIter, string ext, Jacobian &jco)
{

	//ofstream &jout = file_manager.open_ofile_ext(ext, ios::out | ios::binary);
	if (jco.get_matrix_ptr()->nonZeros() == 0)
	{
		stringstream ss;
		ss << "WARNING: jacobian matrix has no non-zeros - the parameter pertubations " << endl;
		ss << "         have no effect on the control file observations." << endl;
		ss << "         This usually means something is not setup correctly." << endl;
		cout << ss.str();
		file_manager.rec_ofstream() << ss.str();

	}
	vector<string> obs_names;
	vector<string> par_names;

	if (isBaseIter)
	{
		obs_names = pest_scenario.get_ctl_ordered_obs_names();
		par_names = pest_scenario.get_ctl_ordered_par_names();
		vector<string> pi_names = pest_scenario.get_ctl_ordered_pi_names();
		obs_names.insert(obs_names.end(), pi_names.begin(), pi_names.end());
		vector<string> jco_par_names = jco.get_base_numeric_par_names();
		if (par_names.size() != jco_par_names.size())
		{
			/*cout << "warning: base parameters missing from binary jco..." << endl;
			cout << "         no need to write control file order." << endl;
			isBaseIter = false;
			par_names = jco.get_base_numeric_par_names();
			obs_names = jco.get_sim_obs_names();*/
			auto new_end = std::remove_if(par_names.begin(), par_names.end(), [&](string &pname)
			{
				return find(jco_par_names.begin(), jco_par_names.end(), pname) == jco_par_names.end();
			});
			par_names.erase(new_end,par_names.end());
		}
	}
	else
	{
		par_names = jco.get_base_numeric_par_names();
		obs_names = jco.get_sim_obs_names();
	}
	string filename = file_manager.build_filename(ext);
	Eigen::SparseMatrix<double> matrix = jco.get_matrix(obs_names,par_names);
	pest_utils::save_binary(filename, obs_names, par_names, matrix);
	//int n_par = par_names.size();
	//int n_obs_and_pi = obs_names.size();
	//int n;
	//int tmp;
	//double data;
	//char par_name[12];
	//char obs_name[20];

	//// write header
	//tmp = -n_par;
	//jout.write((char*)&tmp, sizeof(tmp));
	//tmp = -n_obs_and_pi;
	//jout.write((char*)&tmp, sizeof(tmp));

	////write number nonzero elements in jacobian (includes prior information)
	////n = matrix.nonZeros();
	//n = jco.get_nonzero();
	//jout.write((char*)&n, sizeof(n));

	////write matrix
	//n = 0;
	//map<string, double>::const_iterator found_pi_par;
	//map<string, double>::const_iterator not_found_pi_par;

	////Eigen::SparseMatrix<double> matrix_T(matrix);
	//Eigen::SparseMatrix<double> matrix_T;
	//if (isBaseIter)
	//	matrix_T = jco.get_matrix(obs_names,par_names);
	//else
	//	matrix_T = jco.get_matrix();
	//matrix_T.transpose();
	//for (int icol = 0; icol<jco.get_outersize(); ++icol)
	//{
	//	for (SparseMatrix<double>::InnerIterator it(matrix_T, icol); it; ++it)
	//	{
	//		data = it.value();
	//		n = it.row() + 1 + it.col() * matrix_T.rows();
	//		jout.write((char*)&(n), sizeof(n));
	//		jout.write((char*)&(data), sizeof(data));
	//	}
	//}
	////save parameter names
	//for (vector<string>::const_iterator b = par_names.begin(), e = par_names.end();
	//	b != e; ++b) {
	//	string_to_fortran_char(*b, par_name, 12);
	//	jout.write(par_name, 12);
	//}

	////save observation and Prior information names
	//for (vector<string>::const_iterator b = obs_names.begin(), e = obs_names.end();
	//	b != e; ++b) {
	//	string_to_fortran_char(*b, obs_name, 20);
	//	jout.write(obs_name, 20);
	//}
	////save observation names (part 2 prior information)
	//file_manager.close_file(ext);
}

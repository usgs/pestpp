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

#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <numeric>
#include "Pest.h"
#include "utilities.h"
#include "pest_error.h"
#include <sstream>
#include "pest_data_structs.h"
#include "Transformation.h"
#include "FileManager.h"


using namespace::std;
using namespace::pest_utils;



Pest::Pest() : base_par_transform("PEST base_par_transform"), regul_scheme_ptr(0)
{
}

void Pest::set_defaults()
{
	svd_info.maxsing = 0;
	svd_info.eigthresh = 1.0e-6;
}

void Pest::check_inputs(ostream &f_rec)
{
	if (other_lines.size() > 0)
	{
		stringstream ss;
		ss << "Note: " << other_lines.size() << " unused lines in pest control file, see rec file..." << endl;
		cout << ss.str();
		f_rec << "Note: " << other_lines.size() << " unused lines in pest control file:" << endl;
		for (auto &line : other_lines)
		{
			f_rec <<"  -->  line number " << line.first << ": '" << line.second << "' " << endl;
		}
	}

	if ((control_info.noptmax == 0) && (pestpp_options.get_max_run_fail() > 1))
	{
		cout << "noptmax = 0, resetting max_run_fail = 1" << endl;
		f_rec << "noptmax = 0, resetting max_run_fail = 1" << endl;
		pestpp_options.set_max_run_fail(1);

	}

	Parameters numeric_pars = base_par_transform.ctl2model_cp(ctl_parameters);
	if (svd_info.maxsing == 0) {
		svd_info.maxsing = min(numeric_pars.size(), observation_values.size());
	}

	if (control_info.facparmax <= 1.0)
		throw PestError("'facparmax' must be greater than 1.0");

	if (control_info.pestmode == ControlInfo::PestMode::PARETO)
	{
		bool found_pareto = false, found_other = false;
		if (pareto_info.obsgroup.substr(0, 5) == "REGUL")
			found_pareto = true;
		else if (find(ctl_ordered_obs_group_names.begin(), ctl_ordered_obs_group_names.end(), pareto_info.obsgroup) == ctl_ordered_obs_group_names.end())
			throw PestError("pareto obsgroup not found: " + pareto_info.obsgroup);
		//make sure at least one other obs group has a nonzero weight obs in it


		for (auto &on : ctl_ordered_obs_names)
		{
				if (observation_info.get_group(on) == pareto_info.obsgroup)
					found_pareto = true;
				else if (observation_info.get_weight(on) > 0.0)
					found_other = true;
		}

		if (!found_pareto || !found_other)
		{
			for (auto &pi : ctl_ordered_pi_names)
			{
				if (prior_info.get_pi_rec_ptr(pi).get_group() == pareto_info.obsgroup)
					found_pareto = true;
				else if (prior_info.get_pi_rec_ptr(pi).get_weight() > 0.0)
					found_other = true;
			}
		}
		if (!found_pareto)
			throw PestError("no non-zero weighted obs found in pareto obsgroup: " + pareto_info.obsgroup);
		if (!found_other)
			throw PestError("no non-zero weighted obs found outside of pareto obsgroup");
	}

	vector<string> par_warnings;
	vector<string> par_problems;
	bool unfixed_par = false;
	int par_ub = 0;
	int par_lb = 0;
	bool forgive_bound = false;
	if (control_info.noptmax == 0)
		forgive_bound = true;
	for (auto &pname : ctl_ordered_par_names)
	{
		//double pval = ctl_parameters[pname];
		//double lb = ctl_parameter_info.get_low_bnd(pname);
		const ParameterRec *prec = ctl_parameter_info.get_parameter_rec_ptr(pname);
		if (prec->tranform_type != ParameterRec::TRAN_TYPE::FIXED)
			unfixed_par = true;
		if (prec->init_value < prec->lbnd)
			if (forgive_bound)
				par_warnings.push_back(pname + " is less than lower bound, but noptmax=0, continuing...");
			else
				par_problems.push_back(pname + " is less than lower bound");
		else if (prec->init_value == prec->lbnd)
		{
			//par_warnings.push_back(pname + " is at lower bound");
			par_lb++;
		}
		if (prec->init_value > prec->ubnd)
			if (forgive_bound)
				par_warnings.push_back(pname + " is greater than upper bound, but noptmax=0, continuing...");
			else
				par_problems.push_back(pname + " is greater than upper bound");
		else if (prec->init_value == prec->ubnd)
		{
			//par_warnings.push_back(pname + " is at upper bound");
			par_ub++;
		}
		if (prec->dercom > 1)
		{
			par_warnings.push_back(pname + " has 'dercom' > 1, pestpp suite doesn't support 'dercom' > 1, ignoring");
		}
		if (((prec->tranform_type != ParameterRec::TRAN_TYPE::FIXED) && (prec->tranform_type != ParameterRec::TRAN_TYPE::TIED)) && 
			(prec->chglim != "RELATIVE") && (prec->chglim != "FACTOR"))
			par_problems.push_back(pname + " 'parchglim not in ['factor','relative']: " + prec->chglim);
		
		if ((prec->ubnd > 0.0) && (prec->lbnd < 0.0))
		{
			if (prec->chglim == "FACTOR")
				par_problems.push_back(pname + " 'factor' parchglim not compatible with bounds that cross zero");
			else if ((prec->chglim == "RELATIVE") && (control_info.relparmax < 1.0))
				par_problems.push_back(pname + "bounds cross zero, requires 'relparmax' > 1.0");
		}
	}

	bool err = false;

	if (par_lb > 0)
	{
		cout << "parameter warning: " << par_lb << " parameters are at lower bound" << endl;
		f_rec << "parameter warning: " << par_lb << " parameters are at lower bound" << endl;
	}

	if (par_ub > 0)
	{
		cout << "parameter warning: " << par_ub << " parameters are at upper bound" << endl;
		f_rec << "parameter warning: " << par_ub << " parameters are at upper bound" << endl;
	}
	for (auto &str : par_warnings)
	{
		cout << "parameter warning: " << str << endl;
		f_rec << "parameter warning: " << str << endl;
	}
	for (auto &str : par_problems)
	{
		cout << "parameter error: " << str << endl;
		f_rec << "parameter error: " << str << endl;

		err = true;
	}

	if (!unfixed_par)
	{
		cout << "parameter error: no adjustable parameters" << endl;
		f_rec << "parameter error: no adjustable parameters" << endl;
		err = true;
	}

	if (err)
		throw runtime_error("error in parameter data");

	int n_base = get_pestpp_options().get_n_iter_base();
	if (n_base == -1 || n_base > 0)
	{
	}
	else
	{
		stringstream ss;
		ss << "pest++ option 'n_iter_base' must either be -1 or greater than 0, not " << n_base;
		f_rec << "pest++ option 'n_iter_base' must either be -1 or greater than 0, not " << n_base;
		throw PestError(ss.str());
	}

	int n_super = get_pestpp_options().get_n_iter_super();
	if (n_super < 0)
	{
		stringstream ss;
		ss << "pest++ option 'n_iter_super' must be >= 0, not " << n_super;
		f_rec << "pest++ option 'n_iter_super' must be >= 0, not " << n_super;
		throw PestError(ss.str());
	}

	if ((n_base == -1) && (n_super == 0))
	{
		stringstream ss;
		ss << "pest++ option 'n_iter_base' == -1 so 'n_iter_super' must be > 0, not " << n_super;
		f_rec << "pest++ option 'n_iter_base' == -1 so 'n_iter_super' must be > 0, not " << n_super;
		throw PestError(ss.str());
	}

	//check that prediction names are list in obs
	if ((pestpp_options.get_uncert_flag()) && (pestpp_options.get_prediction_names().size() > 0))
	{
		vector<string> missing;
		for (auto &pred_name : pestpp_options.get_prediction_names())
		{
			auto it_obs = find(ctl_ordered_obs_names.begin(), ctl_ordered_obs_names.end(), pred_name);
			if (it_obs == ctl_ordered_obs_names.end())
			{
				missing.push_back(pred_name);
			}
		}
		if (missing.size() > 0)
		{
			stringstream ss;
			ss << "Pest::check_inputs() the following predictions were not found in the observation names: ";
			for (auto &m : missing)
				ss << m << ',';
			f_rec << ss.str() << endl;
			cout << ss.str() << endl;
			throw PestError(ss.str());
		}
	}

	if (pestpp_options.get_hotstart_resfile().size() > 0)
		if (pestpp_options.get_basejac_filename().size() == 0)
		{
			cout << "'base_jacobian' is none, so 'hotstart_resfile' being ignored..." << endl;
			f_rec << "'base_jacobian' is none, so 'hotstart_resfile' being ignored..." << endl;
		}

	

}

void Pest::check_io()
{

	if (model_exec_info.tplfile_vec.size() == 0)
	{
		cout << "Error: number of template files = 0" << endl;
		throw runtime_error("number of template files = 0");
	}
	if (model_exec_info.insfile_vec.size() == 0)
	{
		cout << "Error: number of instruction files = 0" << endl;
		throw runtime_error("number of instruction files = 0");
	}
	//make sure we can atleast access the model IO files
	vector<string> inaccessible_files;
	for (auto &file : model_exec_info.insfile_vec)
		if (!check_exist_in(file)) inaccessible_files.push_back(file);
	for (auto &file : model_exec_info.outfile_vec)
		if (!check_exist_out(file)) inaccessible_files.push_back(file);
	for (auto &file : model_exec_info.tplfile_vec)
		if (!check_exist_in(file)) inaccessible_files.push_back(file);
	for (auto &file : model_exec_info.inpfile_vec)
		if (!check_exist_out(file)) inaccessible_files.push_back(file);

	if (inaccessible_files.size() != 0)
	{
		string missing;
		for (auto &file : inaccessible_files)
			missing += file + " , ";
		cout << "Could not access the following model interface files: " << missing;
		throw PestError("Could not access the following model interface files: "+missing);
		//cout << "WARNING: could not access the following model interface files: " << missing << endl;

	}
}

const vector<string> Pest::get_ctl_ordered_nz_obs_names()
{
	vector<string> nz_obs;
	for (auto &oname : ctl_ordered_obs_names)
		if (observation_info.get_observation_rec_ptr(oname)->weight > 0.0)
			nz_obs.push_back(oname);
	return nz_obs;
}

const vector<string> Pest::get_ctl_ordered_adj_par_names()
{
	vector<string> adj_pars;
	ParameterRec::TRAN_TYPE ttype;
	for (auto &pname : ctl_ordered_par_names)
	{
		ttype = ctl_parameter_info.get_parameter_rec_ptr(pname)->tranform_type;
		if ((ttype != ParameterRec::TRAN_TYPE::FIXED) && (ttype != ParameterRec::TRAN_TYPE::TIED))
			adj_pars.push_back(pname);
	}
	return adj_pars;
}

const map<string, string> Pest::get_observation_groups() const
{
	map<string, string> obs_grp_map;
	for(auto &iobs : ctl_ordered_obs_names)
	{
		obs_grp_map[iobs] = (observation_info.get_observation_rec_ptr(iobs)->group);
	}
	return obs_grp_map;
}

vector<string> Pest::get_nonregul_obs() const
{
	vector<string>ret_val;
	for(Observations::const_iterator b=observation_values.begin(), e=observation_values.end();
		b!=e; ++b) {
		if (!observation_info.is_regularization((*b).first)) {
			ret_val.push_back((*b).first);
		}
	}
	return ret_val;
}

int Pest::process_ctl_file(ifstream &fin, string pst_filename)
{
	ofstream f_out("ctl_process.out");
	return process_ctl_file(fin, pst_filename, f_out);
}

int Pest::process_ctl_file(ifstream &fin, string _pst_filename, ofstream &f_rec)
{
	string line;
	string line_upper;
	string section("");
	vector<string> tokens;
	int sec_begin_lnum, sec_lnum;
	double value;
	string name;
	string *trans_type;
	string prior_info_string;
	pair<string, string> pi_name_group;
	int lnum;
	int num_par;
	int num_tpl_file;
	int dercom;
	int i_tpl_ins = 0;
	double phimlim;
	double phimaccept;
	double fracphim;
	double wfinit = 1.0;
	double wfmin = numeric_limits<double>::min();
	double wfmax = numeric_limits<double>::max();
	double wffac;
	double wftol;
	bool use_dynamic_reg = false;
	bool reg_adj_grp_weights = false;
	vector<string> pestpp_input;
	regul_scheme_ptr = new DynamicRegularization(use_dynamic_reg);

	//dont change these text names - they are used in ParamTransformSeq
	TranTied *t_tied = new TranTied("PEST to model tied transformation");
	TranOffset *t_offset = new TranOffset("PEST to model offset transformation");
	TranScale *t_scale = new TranScale("PEST to model scale transformation");
	TranLog10 *t_log = new TranLog10("PEST to model log transformation");
	TranFixed *t_fixed = new TranFixed("PEST to model fixed transformation");
	//TranNormalize *t_auto_norm = new TranNormalize("PEST auto-normalization transformation");

	base_par_transform.push_back_ctl2model(t_scale);
	base_par_transform.push_back_ctl2model(t_offset);
	base_par_transform.push_back_ctl2active_ctl(t_tied);
	base_par_transform.push_back_ctl2active_ctl(t_fixed);
	base_par_transform.push_back_active_ctl2numeric(t_log);
	

	set<string> tied_names;
	pst_filename = _pst_filename;

#ifndef _DEBUG
	try {
#endif
	prior_info_string = "";
	for(lnum=1, sec_begin_lnum=1; getline(fin, line); ++ lnum)
	{
		strip_ip(line);
		line_upper = upper_cp(line);
		tokens.clear();
		tokenize(line_upper, tokens);
		sec_lnum = lnum - sec_begin_lnum;
		
		if (lnum == 1)
		{
			if (tokens[0] != "PCF")
			{
				cout << "WARNING: fist line of control file should be 'PCF' not " << tokens[0] << endl;
			}
		}

		else if (tokens.empty())
		{
			//skip blank line
		}
		else if (line[0] == '#')
		{

		}
		else if (line_upper.substr(0,2) == "++")
		{
			pestpp_input.push_back(line);
		}

		else if (line_upper[0] == '*')
		{
			section = upper_cp(strip_cp(line_upper, "both", " *\t\n"));
			sec_begin_lnum = lnum;
		}
		else if (section == "CONTROL DATA")
		{
			if (sec_lnum == 1)
			{
				if (tokens[1] == "REGULARIZATION" || tokens[1] == "REGULARISATION")
				{
					use_dynamic_reg = true;
					control_info.pestmode = ControlInfo::PestMode::REGUL;

				}
				else if (tokens[1] == "ESTIMATION")
					control_info.pestmode = ControlInfo::PestMode::ESTIMATION;
				else if (tokens[1] == "PARETO")
					control_info.pestmode = ControlInfo::PestMode::PARETO;
			}


			else if (sec_lnum == 2)
			{
				convert_ip(tokens[0], num_par);
			}
			else if (sec_lnum == 3)
			{
				convert_ip(tokens[0], num_tpl_file);
				if(tokens.size() >= 5) {
					convert_ip(tokens[4], control_info.numcom);
				}
				else {
					control_info.numcom = 0;
				}
				if(tokens.size() >= 6) {
					convert_ip(tokens[5], control_info.jacfile);
				}
				else {
					control_info.jacfile = 0;
				}
			}
			else if (sec_lnum == 5)
			{
				convert_ip(tokens[0], control_info.relparmax);
				convert_ip(tokens[1], control_info.facparmax);
				convert_ip(tokens[2], control_info.facorig);
			}
			else if (sec_lnum == 6)
			{
				// remove text arguements from the line as these can be specified out of order
				// and PEST++ does not use them
				set<string> remove_tags = { "aui", "auid", "noaui", "senreuse", "nsenreuse", "boundscale", "noboundscale" };
				auto end_iter = std::remove_if(tokens.begin(), tokens.end(),
					[&remove_tags](string &str)->bool{return (remove_tags.find(upper_cp(str)) != remove_tags.end()
					|| remove_tags.find(lower_cp(str)) != remove_tags.end()); });
				tokens.resize(std::distance(tokens.begin(), end_iter));

				convert_ip(tokens[0], control_info.phiredswh);
				if (tokens.size() >= 2) {
					convert_ip(tokens[1], control_info.noptswitch);
				}
				if (tokens.size() >= 3) {
					convert_ip(tokens[2], control_info.splitswh);
				}

			}
			else if (sec_lnum == 7)
			{
				convert_ip(tokens[0], control_info.noptmax);
				convert_ip(tokens[1], control_info.phiredstp);
				convert_ip(tokens[2], control_info.nphistp);
				convert_ip(tokens[3], control_info.nphinored);
				convert_ip(tokens[4], control_info.relparstp);
				convert_ip(tokens[5], control_info.nrelpar);
			}
			else
			{
				other_lines[lnum] = line;
			}
		}
		else if (section == "SINGULAR VALUE DECOMPOSITION")
		{
			if (sec_lnum == 2) {
				convert_ip(tokens[0], svd_info.maxsing);
				convert_ip(tokens[1], svd_info.eigthresh);
			}
			else if (sec_lnum == 3) {
				convert_ip(tokens[0], svd_info.eigwrite);
			}
			else
			{
				other_lines[lnum] = line;
			}
		}
		else if (section == "PARAMETER GROUPS")
		{
			ParameterGroupRec pgi;
			name = tokens[0];
			size_t n_tokens = tokens.size();
			pgi.name = name;
			ctl_ordered_par_group_names.push_back(name);
			convert_ip(tokens[1], pgi.inctyp);
			convert_ip(tokens[2], pgi.derinc);
			convert_ip(tokens[3], pgi.derinclb);
			convert_ip(tokens[4], pgi.forcen);
			convert_ip(tokens[5], pgi.derincmul);
			convert_ip(tokens[6], pgi.dermthd);
			if (n_tokens >= 8) convert_ip(tokens[7], pgi.splitthresh);
			if (n_tokens >= 9) convert_ip(tokens[8], pgi.splitreldiff);
			base_group_info.insert_group(name, pgi);
		}
		else if (section == "PARAMETER DATA")
		{
			if (sec_lnum <= num_par) {
				double scale;
				double offset;
				ParameterRec pi;
				name = tokens[0];
				trans_type = &tokens[1];
				convert_ip(tokens[2], pi.chglim);
				convert_ip(tokens[3], pi.init_value);
				convert_ip(tokens[4], pi.lbnd);
				convert_ip(tokens[5], pi.ubnd);
				convert_ip(tokens[6], pi.group);
				convert_ip(tokens[7], scale);
				convert_ip(tokens[8], offset);
				if (control_info.numcom > 1)
					convert_ip(tokens[9], pi.dercom);
				else
					pi.dercom = 1;
				pi.scale = scale;
				pi.offset = offset;
				// add parameters to model parameter and paramter_info datasets
				ctl_ordered_par_names.push_back(name);
				if (*trans_type == "FIXED")
				{
					pi.tranform_type = ParameterRec::TRAN_TYPE::FIXED;
				}
				else if (*trans_type == "LOG")
				{
					pi.tranform_type = ParameterRec::TRAN_TYPE::LOG;
					n_adj_par++;
				}
				else if (*trans_type == "TIED")
				{
					pi.tranform_type = ParameterRec::TRAN_TYPE::TIED;
				}
				else if (*trans_type == "NONE")
				{
					pi.tranform_type = ParameterRec::TRAN_TYPE::NONE;
					n_adj_par++;
				}
				else
				{
					//pi.tranform_type = ParameterRec::TRAN_TYPE::NONE;
					//assert(true);
					//n_adj_par++;
					throw PestError("unrecognized partrans for par " + name + ": " + *trans_type);
				}
				ctl_parameter_info.insert(name, pi);
				ctl_parameters.insert(name, pi.init_value);
				base_group_info.insert_parameter_link(name, pi.group);

				// build appropriate transformations
				if (*trans_type == "FIXED") {
					t_fixed->insert(name, pi.init_value);}
				else if (*trans_type == "LOG") {
					t_log->insert(name);
				}
				if (offset!=0) {
					t_offset->insert(name, offset);
				}
				if (scale !=1) {
					t_scale->insert(name, scale);
				}
			}
			// Get rest of information for tied paramters
			else {
				name = tokens[0];
				string name_tied =  tokens[1];
				double ratio =  ctl_parameters[name] / ctl_parameters[name_tied];
				t_tied->insert(name, pair<string, double>(name_tied, ratio));
				tied_names.insert(name_tied);
			}
		}
		else if (section == "OBSERVATION GROUPS")
		{
			string name = tokens[0];
			if (tokens.size() > 1)
			{
				stringstream ss;
				ss << "observation covariance matrix detected for group '" << tokens[0] << "' - these are not supported...yet!";
				string s = ss.str();
				throw PestError(s);
			}
			ObservationGroupRec group_rec;
			observation_info.groups[name] = group_rec;
			vector<string>::iterator is = find(ctl_ordered_obs_group_names.begin(), ctl_ordered_obs_group_names.end(), name);
			if (is == ctl_ordered_obs_group_names.end())
			{
				ctl_ordered_obs_group_names.push_back(name);
			}
		}
		else if (section == "OBSERVATION DATA")
		{
			ObservationRec obs_i;
			name = tokens[0];
			convert_ip(tokens[1], value);
			convert_ip(tokens[2], obs_i.weight);
			obs_i.group = tokens[3];
			ctl_ordered_obs_names.push_back(name);
			observation_info.observations[name] = obs_i;
			observation_values.insert(name, value);
		}

		else if (section == "PRIOR INFORMATION")
		{
			//This section processes the prior information.  It does not write out the
			//last prior infomration.  THis is because it must check for line continuations
			if (!prior_info_string.empty() && tokens[0] != "&"){
				pi_name_group = prior_info.AddRecord(prior_info_string);
				ctl_ordered_pi_names.push_back(pi_name_group.first);
				vector<string>::iterator is = find(ctl_ordered_obs_group_names.begin(), ctl_ordered_obs_group_names.end(), pi_name_group.second);
				if (is == ctl_ordered_obs_group_names.end())
				{
					ctl_ordered_obs_group_names.push_back(pi_name_group.second);
				}
				prior_info_string.clear();
			}
			else if (tokens[0] == "&") {
				prior_info_string.append(" ");
			}
			prior_info_string.append(line_upper);
		}

		else if (section == "PRIOR INFORMATION" )
		{
			//This section processes the prior information.  It does not write out the
			//last prior infomration
			if (!prior_info_string.empty() && tokens[0] != "&") {
				pi_name_group = prior_info.AddRecord(prior_info_string);
				ctl_ordered_pi_names.push_back(pi_name_group.first);
				vector<string>::iterator is = find(ctl_ordered_obs_group_names.begin(), ctl_ordered_obs_group_names.end(), pi_name_group.second);
				if (is == ctl_ordered_obs_group_names.end())
				{
					ctl_ordered_obs_group_names.push_back(pi_name_group.second);
				}
				prior_info_string.clear();
			}
			else if (tokens[0] != "&") {
				prior_info_string.append(" ");
			}
			prior_info_string.append(line);
		}
		else if (section == "MODEL COMMAND LINE" )
		{
			model_exec_info.comline_vec.push_back(line);
		}
		else if (section == "MODEL INPUT/OUTPUT" )
		{
			vector<string> tokens_case_sen;
			tokenize(line, tokens_case_sen);
			if(i_tpl_ins < num_tpl_file)
			{
				model_exec_info.tplfile_vec.push_back(tokens_case_sen[0]);
				model_exec_info.inpfile_vec.push_back(tokens_case_sen[1]);
			}
			else
			{
				model_exec_info.insfile_vec.push_back(tokens_case_sen[0]);
				model_exec_info.outfile_vec.push_back(tokens_case_sen[1]);
			}
			++i_tpl_ins;
		}
		else if (section == "REGULARISATION" || section=="REGULARIZATION" )
		{
			if (sec_lnum == 1) {
				convert_ip(tokens[0], phimlim);
				convert_ip(tokens[1], phimaccept);
				fracphim = 0.0;
				if(tokens.size() >=3) convert_ip(tokens[2], fracphim);
			}
			else if (sec_lnum == 2) {
				convert_ip(tokens[0], wfinit);
				convert_ip(tokens[1], wfmin);
				convert_ip(tokens[2], wfmax);
			}
			else if (sec_lnum == 3) {
				int iregadj;
				convert_ip(tokens[0], wffac);
				convert_ip(tokens[1], wftol);
				if (tokens.size() > 2)
				{
					convert_ip(tokens[2], iregadj);
					if (iregadj == 1) reg_adj_grp_weights = true;
				}
				delete regul_scheme_ptr;
				regul_scheme_ptr = new DynamicRegularization(use_dynamic_reg, reg_adj_grp_weights, phimlim,
					phimaccept, fracphim, wfmin, wfmax, wffac, wftol, wfinit);
			}
			else
			{
				other_lines[lnum] = line;
			}
		}
		else if (section == "PARETO")
		{
			if (sec_lnum == 1)
			{
				convert_ip(tokens[0], pareto_info.obsgroup);
			}
			else if (sec_lnum == 2)
			{
				convert_ip(tokens[0], pareto_info.wf_start);
				convert_ip(tokens[1], pareto_info.wf_fin);
				convert_ip(tokens[2], pareto_info.wf_inc);
			}
			else if (sec_lnum == 3)
			{
				convert_ip(tokens[0], pareto_info.niter_start);
				convert_ip(tokens[1], pareto_info.niter_gen);
				convert_ip(tokens[2], pareto_info.niter_fin);
			}

		}
		else
		{
			other_lines[lnum] = line;
		}
	}

	// write out last prior information record
	if (!prior_info_string.empty())
	{
		pi_name_group = prior_info.AddRecord(prior_info_string);
		ctl_ordered_pi_names.push_back(pi_name_group.first);
		vector<string>::iterator is = find(ctl_ordered_obs_group_names.begin(), ctl_ordered_obs_group_names.end(), pi_name_group.second);
		if (is == ctl_ordered_obs_group_names.end())
		{
			ctl_ordered_obs_group_names.push_back(pi_name_group.second);
		}
		prior_info_string.clear();
	}
#ifndef _DEBUG
	}
	catch (PestConversionError &e) {
		std::stringstream out;
		out << "Error parsing \"" << pst_filename << "\" on line number " << lnum << endl;
		out << e.what() << endl;
		e.add_front(out.str());
		e.raise();
	}
#endif
	fin.close();
	// process pest++ options last
	pestpp_options.set_n_iter_super(0);
	pestpp_options.set_n_iter_base(max(1, control_info.noptmax));
	pestpp_options.set_super_eigthres(svd_info.eigthresh);
	pestpp_options.set_max_n_super(ctl_parameters.size());
	pestpp_options.set_max_super_frz_iter(5);
	pestpp_options.set_max_n_super(n_adj_par);
	pestpp_options.set_max_reg_iter(20);
	pestpp_options.set_uncert_flag(true);
	pestpp_options.set_prediction_names(vector<string>());
	pestpp_options.set_parcov_filename(string());
	pestpp_options.set_obscov_filename(string());
	pestpp_options.set_basejac_filename(string());
	pestpp_options.set_sweep_parameter_csv_file(string());
	pestpp_options.set_sweep_output_csv_file("sweep_out.csv");
	pestpp_options.set_sweep_base_run(false);
	pestpp_options.set_sweep_forgive(false);
	pestpp_options.set_sweep_chunk(500);
	pestpp_options.set_tie_by_group(false);
	pestpp_options.set_enforce_tied_bounds(false);
	
	pestpp_options.set_jac_scale(true);
	pestpp_options.set_opt_obj_func("");
	pestpp_options.set_opt_coin_log(true);
	pestpp_options.set_opt_skip_final(false);
	pestpp_options.set_opt_std_weights(false);
	pestpp_options.set_opt_dec_var_groups(vector<string>());
	pestpp_options.set_opt_ext_var_groups(vector<string>());
	pestpp_options.set_opt_constraint_groups(vector<string>());
	pestpp_options.set_opt_risk(0.5);
	pestpp_options.set_opt_direction(1.0);
	pestpp_options.set_opt_iter_tol(0.001);
	pestpp_options.set_opt_recalc_fosm_every(1);
	pestpp_options.set_opt_iter_derinc_fac(1.0);
	pestpp_options.set_opt_include_bnd_pi(true);
	pestpp_options.set_hotstart_resfile(string());
	pestpp_options.set_ies_par_csv("");
	pestpp_options.set_ies_obs_csv("");
	pestpp_options.set_ies_obs_restart_csv("");
	pestpp_options.set_ies_par_restart_csv("");
	pestpp_options.set_ies_lam_mults(vector<double>());
	pestpp_options.set_ies_init_lam(-999);
	pestpp_options.set_ies_use_approx(true);
	pestpp_options.set_ies_subset_size(5);
	pestpp_options.set_ies_reg_factor(0.0);
	pestpp_options.set_ies_verbose_level(0);
	pestpp_options.set_ies_use_prior_scaling(false);
	pestpp_options.set_ies_num_reals(50);
	pestpp_options.set_ies_bad_phi(1.0e+300);
	pestpp_options.set_ies_bad_phi_sigma(1.0e+300);
	pestpp_options.set_ies_include_base(true);
	pestpp_options.set_ies_use_empirical_prior(false);
	pestpp_options.set_ies_group_draws(true);
	//pestpp_options.set_ies_num_reals_passed(false);
	pestpp_options.set_ies_enforce_bounds(true);
	pestpp_options.set_par_sigma_range(4.0);
	pestpp_options.set_ies_save_binary(false);
	pestpp_options.set_ies_localizer("");
	pestpp_options.set_ies_accept_phi_fac(1.05);
	pestpp_options.set_ies_lambda_inc_fac(10.0);
	pestpp_options.set_ies_lambda_dec_fac(0.75);
	pestpp_options.set_ies_save_lambda_en(false);
	pestpp_options.set_ies_weight_csv("");
	pestpp_options.set_ies_subset_how("RANDOM");
	pestpp_options.set_ies_localize_how("PARAMETERS");
	pestpp_options.set_ies_num_threads(-1);
	pestpp_options.set_ies_debug_fail_subset(false);
	pestpp_options.set_ies_debug_fail_remainder(false);
	pestpp_options.set_ies_debug_bad_phi(false);
	pestpp_options.set_ies_debug_upgrade_only(false);
	pestpp_options.set_ies_debug_high_subset_phi(false);
	pestpp_options.set_ies_debug_high_upgrade_phi(false);
	pestpp_options.set_ies_csv_by_reals(true);
	pestpp_options.set_ies_autoadaloc(false);
	pestpp_options.set_ies_autoadaloc_sigma_dist(1.0);
	pestpp_options.set_ies_enforce_chglim(false);
	pestpp_options.set_ies_center_on("");

	pestpp_options.set_gsa_method("MORRIS");
	//many of these defaults are also redefined in gsa main
	pestpp_options.set_gsa_morris_p(4);
	pestpp_options.set_gsa_morris_r(4);
	pestpp_options.set_gsa_morris_delta(0.6666);
	pestpp_options.set_gsa_morris_obs_sen(true);
	pestpp_options.set_gsa_morris_pooled_obs(false);
	pestpp_options.set_gsa_sobol_par_dist("norm");
	pestpp_options.set_gsa_sobol_samples(4);
	pestpp_options.set_gsa_rand_seed(2);

	pestpp_options.set_condor_submit_file(string());
	pestpp_options.set_overdue_giveup_minutes(1.0e+30);

	for(vector<string>::const_iterator b=pestpp_input.begin(),e=pestpp_input.end();
		b!=e; ++b) {

			pestpp_options.parce_line(*b);
	}
	
	regul_scheme_ptr->set_max_reg_iter(pestpp_options.get_max_reg_iter());
//	//Make sure we use Q1/2J is PROPACK is chosen
//	if (pestpp_options.get_svd_pack() == PestppOptions::SVD_PACK::PROPACK)
//	{
//		pestpp_options.set_mat_inv(PestppOptions::MAT_INV::Q12J);
//	}

	if (pestpp_options.get_tie_by_group())
	{
		cout << "Note: ++tie_by_group(true) - tying adjustable parameters by groups" << endl;
		f_rec << "Note: ++tie_by_group(true) - tying adjustable parameters by groups" << endl;
		map<string, vector<string>> group_map;
		string gname;
		ParameterRec::TRAN_TYPE tlog = ParameterRec::TRAN_TYPE::LOG, 
			tnone = ParameterRec::TRAN_TYPE::NONE, 
			ttied = ParameterRec::TRAN_TYPE::TIED;
		int new_n_adj_par = 0;
		for (auto pname : ctl_ordered_par_names)
		{
			if ((ctl_parameter_info.get_parameter_rec_ptr(pname)->tranform_type == tlog) ||
				(ctl_parameter_info.get_parameter_rec_ptr(pname)->tranform_type == tnone))
			{

				if (tied_names.find(pname) != tied_names.end())
				{
					new_n_adj_par++;
					continue;
				}
				gname = ctl_parameter_info.get_parameter_rec_ptr(pname)->group;
				if (group_map.find(gname) == group_map.end())
					group_map[gname] = vector<string>();
				group_map[gname].push_back(pname);
			}
			
		}
		string tie_to_name;
		vector<string> to_tie_names;
		vector<string>::const_iterator first, last;
		for (auto gm : group_map)
		{
			tie_to_name = gm.second[0];
			new_n_adj_par++;
			first = gm.second.begin() + 1;
			last = gm.second.end();
			to_tie_names = vector<string>(first, last);
			for (auto pname : to_tie_names)
			{
				double ratio = ctl_parameters[pname] / ctl_parameters[tie_to_name];
				t_tied->insert(pname, pair<string, double>(tie_to_name, ratio));
				ctl_parameter_info.get_parameter_rec_ptr_4_mod(pname)->tranform_type = ttied;
			}

		}
		f_rec << "-->number of adjustable parameters reduced from " << n_adj_par << " to " << new_n_adj_par << endl;
		cout << "-->number of adjustable parameters reduced from " << n_adj_par << " to " << new_n_adj_par << endl;
		n_adj_par = new_n_adj_par;
		if (tied_names.size() > 0)
		{
			f_rec << "-->existing adjustable parameters that others tie to have been maintained:" << endl;
			for (auto tname : tied_names)
				f_rec << tname << endl;
		}
	}

	return 0;
}


const vector<string> &Pest::get_comline_vec()
{
	return model_exec_info.comline_vec;
}
const vector<string> &Pest::get_tplfile_vec()
{
	return model_exec_info.tplfile_vec;
}
const  vector<string> &Pest::get_inpfile_vec()
{
	return model_exec_info.inpfile_vec;
}
const  vector<string> &Pest::get_insfile_vec()
{
	return model_exec_info.insfile_vec;
}
const vector<string> &Pest::get_outfile_vec()
{
	return model_exec_info.outfile_vec;
}

void Pest::enforce_par_limits(Parameters & upgrade_active_ctl_pars, const Parameters &last_active_ctl_pars, bool enforce_chglim, bool enforce_bounds)
{
	if ((!enforce_chglim) && (!enforce_bounds))
		return;
	stringstream ss;
	double fpm = control_info.facparmax;
	double facorig = control_info.facorig;
	double rpm = control_info.relparmax;
	double orig_val, last_val, fac_lb, fac_ub, rel_lb, rel_ub, eff_ub, eff_lb,chg_lb, chg_ub;
	double chg_fac, chg_rel;
	double scaling_factor = 1.0;
	string parchglim;
	double bnd_tol = 0.001;
	double scaled_bnd_val;
	string controlling_par = "";
	ParameterInfo &p_info = ctl_parameter_info;
	const ParameterRec *p_rec;
	Parameters upgrade_ctl_pars;
	Parameters last_ctl_pars;

	if (pestpp_options.get_enforce_tied_bounds())
	{
		upgrade_ctl_pars = base_par_transform.active_ctl2ctl_cp(upgrade_active_ctl_pars);
		last_ctl_pars = base_par_transform.active_ctl2ctl_cp(last_active_ctl_pars);
	}
	else
	{
		upgrade_ctl_pars = upgrade_active_ctl_pars;
		last_ctl_pars = last_active_ctl_pars;

	}
	for (auto p : upgrade_ctl_pars)
	{
		last_val = last_ctl_pars.get_rec(p.first);
		
		p_rec = p_info.get_parameter_rec_ptr(p.first);
		parchglim = p_rec->chglim;

		if (p.second == 0.0)
			p.second = p_rec->ubnd / 4.0;
		orig_val = ctl_parameters.get_rec(p.first);
		if (orig_val == 0.0)
			orig_val = p_rec->ubnd / 4.0;

		//apply facorig correction if needed
		if (ctl_parameter_info.get_parameter_rec_ptr(p.first)->tranform_type == ParameterRec::TRAN_TYPE::NONE)
		{
			if (abs(p.second) < abs(orig_val) * facorig)
				p.second = orig_val * facorig;
			if (abs(last_val) < abs(orig_val * facorig))
				last_val = orig_val * facorig;
		}
			


		//calc fac lims
		if (abs(last_val) > abs(p.second))
			chg_fac = last_val / p.second;
		else
			chg_fac = p.second / last_val;
		if (p.second > 0.0)
		{
			fac_lb = last_val / fpm;
			fac_ub = last_val * fpm;
		}
	
		else
		{
			fac_lb = last_val * fpm;
			fac_ub = last_val / fpm;
			
		}

		//calc rel lims
		rel_lb = last_ctl_pars.get_rec(p.first) - (abs(last_val) * rpm);
		rel_ub = last_ctl_pars.get_rec(p.first) + (abs(last_val) * rpm);
		chg_rel = (last_val - p.second) / last_val;

		if (parchglim == "FACTOR")
		{
			chg_lb = fac_lb;
			chg_ub = fac_ub;
		}
		else if (parchglim == "RELATIVE")
		{
			chg_lb = rel_lb;
			chg_ub = rel_ub;
		}
		else
		{
			throw runtime_error("Pest::enforce_par_limits() error: unrecognized 'parchglim': " + parchglim);
		}


		double temp = 1.0;
		if (enforce_chglim)
		{		
			if (p.second > chg_ub)
			{
				temp = abs((chg_ub - last_val) / (p.second - last_val));
				if ((temp > 1.0) || (temp < 0.0))
				{
					ss.str("");
					ss << "Pest::enforce_par_limts() error: invalid upper parchglim scaling factor " << temp << " for par " << p.first << endl;
					ss << " chglim:" << chg_ub << ", last_val:" << last_val << ", current_val:" << p.second << endl;
					throw runtime_error(ss.str());
				}

				if (temp < scaling_factor)
				{
					scaling_factor = temp;
					controlling_par = p.first;
				}
			}

			else if (p.second < chg_lb)
				temp = abs((last_val - chg_lb) / (last_val - p.second));
			if ((temp > 1.0) || (temp < 0.0))
			{
				ss.str("");
				ss << "Pest::enforce_par_limts() error: invalid lower parchglim scaling factor " << temp << " for par " << p.first << endl;
				ss << " chglim:" << chg_lb << ", last_val:" << last_val << ", current_val:" << p.second << endl;
				throw runtime_error(ss.str());
			}
			if (temp < scaling_factor)
			{
				scaling_factor = temp;
				controlling_par = p.first;
			}
		}

		if (enforce_bounds)
		{
			/*if (last_val >= p_rec->ubnd)
			{
				ss.str("");
				ss << "Pest::enforce_par_limits() error: last value for parameter " << p.first << " at upper bound";
				throw runtime_error(ss.str());
			}

			else if (last_val <= p_rec->lbnd)
			{
				ss.str("");
				ss << "Pest::enforce_par_limits() error: last value for parameter " << p.first << " at lower bound";
				throw runtime_error(ss.str());
			}*/
			scaled_bnd_val = p_rec->ubnd + (p_rec->ubnd * bnd_tol);
			if (p.second > scaled_bnd_val)
			{
				temp = abs((p_rec->ubnd - last_val) / (p.second - last_val));
				if ((temp > 1.0) || (temp < 0.0))
				{
					
					ss << "Pest::enforce_par_limts() error: invalid upper bound scaling factor " << temp << " for par " << p.first << endl;
					ss << " ubnd:" << p_rec->ubnd << ", last_val:" << last_val << ", current_val:" << p.second << endl;
					throw runtime_error(ss.str());
				}
				if (temp < scaling_factor)
				{
					scaling_factor = temp;
					controlling_par = p.first;
				}
			}
			scaled_bnd_val = p_rec->lbnd - (p_rec->lbnd * bnd_tol);
			if (p.second < p_rec->lbnd)
			{
				temp = abs((last_val - p_rec->lbnd) / (last_val - p.second));
				if ((temp > 1.0) || (temp < 0.0))
				{
					stringstream ss;
					ss << "Pest::enforce_par_limts() error: invalid lower bound scaling factor " << temp << " for par " << p.first << endl;
					ss << " lbnd:" << p_rec->lbnd << ", last_val:" << last_val << ", current_val:" << p.second << endl;
					throw runtime_error(ss.str());
				}
				if (temp < scaling_factor)
				{
					scaling_factor = temp;
					controlling_par = p.first;
				}
			}
		}	
	}

	if (scaling_factor == 0.0)
	{
		throw runtime_error("Pest::enforce_par_change_limits error : zero length parameter vector");
	}

	if (scaling_factor != 1.0)
	{
		for (auto &p : upgrade_active_ctl_pars)
		{
			
			last_val = last_ctl_pars.get_rec(p.first);
			p.second =last_val + (p.second - last_val) *  scaling_factor;
		}
	}
}

pair<Parameters,Parameters> Pest::get_effective_ctl_lower_upper_bnd(Parameters &pars)
{
	vector<string> keys = pars.get_keys();
	Parameters lbnd = ctl_parameter_info.get_low_bnd(keys);
	Parameters ubnd = ctl_parameter_info.get_up_bnd(keys);
	if (base_par_transform.get_tied_ptr()->get_items().size() == 0)
		return pair<Parameters,Parameters>(lbnd,ubnd);
	auto tt_items = base_par_transform.get_tied_ptr()->get_items();
	double ref_bnd, tie_bnd, tie_ratio,dist_ratio,new_bnd;
	double tie_val, ref_val;
	for (auto &tt_item : tt_items)
	{
		ref_val = pars.get_rec(tt_item.second.first);
		tie_val = pars.get_rec(tt_item.first);
		
		//lower bound
		ref_bnd = ref_val - ctl_parameter_info.get_parameter_rec_ptr(tt_item.second.first)->lbnd;
		tie_bnd = tie_val - ctl_parameter_info.get_parameter_rec_ptr(tt_item.first)->lbnd;
		tie_ratio = tt_item.second.second;
		if ((ref_bnd <= 0) || (tie_bnd <= 0))
			new_bnd = ref_val;
		else
		{
			dist_ratio = tie_bnd / ref_bnd;
			new_bnd = ref_val - (dist_ratio * ref_bnd);
		}
		if (new_bnd > lbnd.get_rec(tt_item.second.first))
		{
			lbnd.update_rec(tt_item.second.first, new_bnd);
		}
		
		//upper bound
		ref_bnd = ctl_parameter_info.get_parameter_rec_ptr(tt_item.second.first)->ubnd - ref_val;
		tie_bnd = ctl_parameter_info.get_parameter_rec_ptr(tt_item.first)->ubnd - tie_val;
		tie_ratio = tt_item.second.second;
		if ((ref_bnd <= 0) || (tie_bnd <= 0))
			new_bnd = ref_val;
		else
		{
			dist_ratio = tie_bnd / ref_bnd;
			new_bnd = ref_val + (dist_ratio * ref_bnd);
		}
		if (new_bnd < ubnd.get_rec(tt_item.second.first))
		{
			ubnd.update_rec(tt_item.second.first, new_bnd);
		}
	}
	return pair<Parameters,Parameters>(lbnd,ubnd);
}

map<string, double> Pest::get_pars_at_near_bounds(const Parameters & pars, double tol)
{
	 map<string, double> bnd_map;
	 const ParameterRec *p_rec;
	 ParameterInfo &pinfo = ctl_parameter_info;
	 double v;
	 for (auto p : pars)
	 {
		 p_rec = pinfo.get_parameter_rec_ptr(p.first);
		 v = (p_rec->ubnd - (tol * p_rec->ubnd));
		 if (p.second >= v)
			 bnd_map[p.first] = v;
		 v = (p_rec->lbnd + (tol * p_rec->lbnd));
		 if (p.second <= v)
			 bnd_map[p.first] = v;

	 }
	 return bnd_map;
}


Pest::~Pest() {
	/*if (regul_scheme_ptr != 0)
	{
		try
		{
			delete regul_scheme_ptr;
		}
		catch (...)
		{
		}
	}*/
}

ostream& operator<< (ostream &os, const Pest& val)
{
	os << val.control_info;
	os << val.svd_info;
	return os;
}

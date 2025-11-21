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
#include <iterator>
#include <unordered_set>
#include <sstream>
#include <algorithm>
#include "Pest.h"
#include "utilities.h"
#include "pest_error.h"
#include "pest_data_structs.h"
#include "Transformation.h"
#include "FileManager.h"
#include "model_interface.h"
#include "Jacobian.h"
#include "QSqrtMatrix.h"
#include <limits>
#include "network_package.h"
#include <cmath>


using namespace::std;
using namespace::pest_utils;



Pest::Pest() : base_par_transform("PEST base_par_transform"), regul_scheme_ptr(0)
{
	set_defaults();
}

void Pest::set_defaults()
{
	pestpp_options.set_defaults();
	svd_info.set_defaults();
	regul_scheme_ptr = 0;//new DynamicRegularization;
	//regul_scheme_ptr->set_defaults();
	control_info.set_defaults();

}

void Pest::set_default_dynreg()
{
    regul_scheme_ptr = new DynamicRegularization;
    regul_scheme_ptr->set_defaults();
}

void Pest::check_inputs(ostream &f_rec, bool forgive, bool forgive_parchglim, int cycle)
{
	if (control_info.noptmax == 0)
	{
		if (!forgive)
		{
			cout << endl << "Note: 'NOPTMAX' == 0, switching to forgiveness mode when checking inputs" << endl << endl;
			f_rec << endl << "Note 'NOPTMAX' == 0, switching to forgiveness mode when checking inputs" << endl << endl;
		}

		forgive = true;
		forgive_parchglim = true;
	}

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
	{
		if (forgive)
		{
			cout << "WARNING: 'facparmax` should be greater than 1.0" << endl;
			f_rec << "WARNING: 'facparmax` should be greater than 1.0" << endl;
		}
		else
		{
			throw PestError("'facparmax' must be greater than 1.0");
		}
	}

	vector<string> par_warnings;
	vector<string> par_problems;
	bool adj_par;
	int par_ub = 0;
	int par_lb = 0;
	vector<string> adj_pnames = get_ctl_ordered_adj_par_names();
	set<string> sadj(adj_pnames.begin(), adj_pnames.end());
    const ParameterRec *prec;
    const map<string,pair<string,double>> tied_map = base_par_transform.get_tied_ptr()->get_items();
    ParameterRec::TRAN_TYPE tranfixed = ParameterRec::TRAN_TYPE::FIXED;
    ParameterRec::TRAN_TYPE trantied = ParameterRec::TRAN_TYPE::TIED;
    ParameterRec::TRAN_TYPE tranlog = ParameterRec::TRAN_TYPE::LOG;

    ParameterRec::TRAN_TYPE tran;
    stringstream ss;
	for (auto &pname : ctl_ordered_par_names) {
        //double pval = ctl_parameters[pname];
        //double lb = ctl_parameter_info.get_low_bnd(pname);
        prec = ctl_parameter_info.get_parameter_rec_ptr(pname);
        tran = prec->tranform_type;
        adj_par = true;
        if ((tran == tranfixed) || (tran == tranfixed))
            adj_par = false;
        if (tran == trantied) {
            string partied = tied_map.at(pname).first;
            if (sadj.find(partied) == sadj.end()) {
                par_problems.push_back("'tied' parameter '" + pname + "' is tied to '" + partied +
                                       "' which is not an adjustable parameter");
            }
        } else if (tran == tranlog)
        {
            if (prec->init_value <= 0.0)
            {
                ss.str("");
                ss << pname << ": log transform and initial value <= 0.0: " << prec->init_value;
                par_problems.push_back(ss.str());
            }
            if (prec->lbnd <= 0.0)
            {
                ss.str("");
                ss << pname << ": log transform and lower bound <= 0.0: " << prec->lbnd;
                par_problems.push_back(ss.str());
            }
            if (prec->ubnd <= 0.0)
            {
                ss.str("");
                ss << pname << ": log transform and upper bound <= 0.0: " << prec->ubnd;
                par_problems.push_back(ss.str());
            }
        }

		if (prec->lbnd >= prec->ubnd)
		{
            ss.str("");
            ss << pname << ": bounds are busted:" << prec->lbnd << " >= " << prec->ubnd;
			if ((forgive) || (!adj_par))
            {
				par_warnings.push_back(ss.str());
			}
			else
			{
				par_problems.push_back(ss.str());
			}
		}
		if (prec->init_value < prec->lbnd)
		{
            ss.str("");
            ss << pname << ": initial value is less than lower bound:" << prec->init_value << " < " << prec->lbnd;
			if ((forgive) || (!adj_par))
			{
				par_warnings.push_back(ss.str());
			}
			else
			{
				par_problems.push_back(ss.str());
			}
		}
		else if (prec->init_value == prec->lbnd)
		{
			par_lb++;
		}
		if (prec->init_value > prec->ubnd) {
            ss.str("");
            ss << pname << ": initial value is greater than upper bound:" << prec->init_value << " > " << prec->ubnd;

            if ((forgive) || (!adj_par)) {
                par_warnings.push_back(ss.str());
            } else {
                par_problems.push_back(ss.str());
            }
        }
		else if (prec->init_value == prec->ubnd)
		{
			par_ub++;
		}
		if (prec->dercom > 1)
		{
			par_warnings.push_back(pname + " has 'dercom' > 1, pestpp suite doesn't support 'dercom' > 1, ignoring");
		}
		if ((!adj_par) && (prec->chglim != "RELATIVE") && (prec->chglim != "FACTOR"))
				par_problems.push_back(pname + " 'parchglim not in ['factor','relative']: " + prec->chglim);
		
		if ((prec->ubnd > 0.0) && (prec->lbnd < 0.0))
		{
			if ((!forgive_parchglim) && (prec->chglim == "FACTOR"))
			{
				if ((forgive) || (!adj_par))
				{
					par_warnings.push_back(pname + " 'factor' parchglim not compatible with bounds that cross zero");
				}
				else
				{
					par_problems.push_back(pname + " 'factor' parchglim not compatible with bounds that cross zero");
			
				}
			}
			else if ((prec->chglim == "RELATIVE") && (control_info.relparmax < 1.0))
			{
				if ((forgive) || (!adj_par))
				{
					par_warnings.push_back(pname + "bounds cross zero, requires 'relparmax' > 1.0");
				}
				else
				{
					par_problems.push_back(pname + "bounds cross zero, requires 'relparmax' > 1.0");
				}
			}
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
	
	for (auto& str : par_warnings)
	{
			
		f_rec << "parameter warning: " << str << endl;
	}
	if (par_warnings.size() > 0)
		cout << par_warnings.size() << " parameter warnings, see rec file " << endl;
	
	for (auto &str : par_problems)
	{
		cout << "parameter error: " << str << endl;
		f_rec << "parameter error: " << str << endl;

		err = true;
	}

	if (get_n_adj_par() == 0)
	{
		if ((forgive) || (cycle != NetPackage::NULL_DA_CYCLE))
		{
			cout << "parameter warning: no adjustable parameters" << endl;
			f_rec << "parameter warning: no adjustable parameters" << endl;
		}
		else
		{
			cout << "parameter error: no adjustable parameters" << endl;
			f_rec << "parameter error: no adjustable parameters" << endl;
			err = true;
		}

	}
	vector<string> obs_problems;
	for (string& oname : ctl_ordered_obs_names) {
		if (observation_info.get_weight(oname) < 0) {
			ss.str("");
			ss << "observation " << oname << " weight less than zero";
			obs_problems.push_back(ss.str());
		}
	}

	for (auto &str : obs_problems)
	{
		cout << "observation error: " << str << endl;
		f_rec << "observation error: " << str << endl;
		err = true;
	}



	if (get_ctl_ordered_nz_obs_names().size() == 0)
	{
		if ((forgive) || (NetPackage::NULL_DA_CYCLE))
		{
			cout << "observation warning: no non-zero weighted observations" << endl;
			f_rec << "observation warning: no non-zero weighted observations" << endl;
		}
		else
		{
			cout << "observation error: no non-zero weighted observations" << endl;
			f_rec << "observation error: no non-zero weighted observations" << endl;
			err = true;
		}
	}

	if (err)
		throw runtime_error("error in inputs...");

	int n_base = get_pestpp_options().get_n_iter_base();
	if (n_base == -1 || n_base > 0)
	{
	}
	else
	{
		stringstream ss;
		ss << "pest++ option 'n_iter_base' must either be -1 or greater than 0, not " << n_base;
		f_rec << "pest++ option 'n_iter_base' must either be -1 or greater than 0, not " << n_base;
		if (!forgive)
			throw PestError(ss.str());
	}

	int n_super = get_pestpp_options().get_n_iter_super();
	if (n_super < 0)
	{
		stringstream ss;
		ss << "pest++ option 'n_iter_super' must be >= 0, not " << n_super;
		f_rec << "pest++ option 'n_iter_super' must be >= 0, not " << n_super;
		if (!forgive)
			throw PestError(ss.str());
	}

	if ((n_base == -1) && (n_super == 0))
	{
		stringstream ss;
		ss << "pest++ option 'n_iter_base' == -1 so 'n_iter_super' must be > 0, not " << n_super;
		f_rec << "pest++ option 'n_iter_base' == -1 so 'n_iter_super' must be > 0, not " << n_super;
		if (!forgive)
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
			if (!forgive)
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

void Pest::check_io(ofstream& f_rec, bool echo_errors)
{
	ModelInterface mi(model_exec_info.tplfile_vec,model_exec_info.inpfile_vec,
		model_exec_info.insfile_vec,model_exec_info.outfile_vec,model_exec_info.comline_vec);

	try
	{
		mi.check_io_access();

	}
	catch (exception& e)
	{
		string mess = e.what();
		throw_control_file_error(f_rec, "error in model interface file access:" + mess, true,echo_errors);
	}
	catch (...)
	{
		throw_control_file_error(f_rec, "unspecified error in model interface file access checking", true, echo_errors);
	}
	if (pestpp_options.get_check_tplins())
	{
		try
		{
			mi.check_tplins(ctl_ordered_par_names,ctl_ordered_obs_names);

		}
		catch (exception& e)
		{
			string mess = e.what();
			throw_control_file_error(f_rec, "error in model interface file checking:" + mess, true, echo_errors);
		}
		catch (...)
		{
			throw_control_file_error(f_rec, "unspecified error in model interface file checking", true, echo_errors);
		}
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

int Pest::process_ctl_file(ifstream& fin, string _pst_filename)
{
	ofstream f_out("ctl_process.out");
	return process_ctl_file(fin, _pst_filename, f_out);
}

bool IsQuote(char c)
{
    switch(c)
    {
        case '\"':
        case '\'':
            return true;
        default:
            return false;
    }
}

int Pest::process_ctl_file(ifstream& fin, string _pst_filename, ofstream& f_rec)
{
	cout << "processing control file " << _pst_filename << endl;
	if (!fin)
    {
	    throw PestError("control file stream is not good");
    }
	string line;
	string line_upper;
	string section("");
	vector<string> tokens;
	int sec_begin_lnum, sec_lnum;
	double value;
	string name;
	string* trans_type;
	pair<string, string> pi_name_group;
	int lnum;
	int num_par = -1;
	int num_tpl_file;
	int dercom;
	int i_tpl_ins = 0;
	bool use_dynamic_reg = false;
	bool reg_adj_grp_weights = false;
	vector<string> pestpp_input;
	
	//dont change these text names - they are used in ParamTransformSeq
	TranTied* t_tied = new TranTied("PEST to model tied transformation");
	TranOffset* t_offset = new TranOffset("PEST to model offset transformation");
	TranScale* t_scale = new TranScale("PEST to model scale transformation");
	TranLog10* t_log = new TranLog10("PEST to model log transformation");
	TranFixed* t_fixed = new TranFixed("PEST to model fixed transformation");
	//TranNormalize *t_auto_norm = new TranNormalize("PEST auto-normalization transformation");

	base_par_transform.push_back_ctl2model(t_scale);
	base_par_transform.push_back_ctl2model(t_offset);
	base_par_transform.push_back_ctl2active_ctl(t_tied);
	base_par_transform.push_back_ctl2active_ctl(t_fixed);
	base_par_transform.push_back_active_ctl2numeric(t_log);

	pestpp_options.set_defaults();
	vector<string> notfound_args;
	set<string> tied_names;
	pst_filename = _pst_filename;
	set<string> sections_found;
	set<string> nonkeyword_sections = { "SINGULAR VALUE DECOMPOSITION","REGULARIZATION",
		"REGULARISATION","PLUSPLUS","CONTROL DATA" };
	stringstream ss;
	PestppOptions::ARG_STATUS stat;
	vector<string> par_group_formal_names{ "PARGPNME","INCTYP","DERINC","DERINCLB","FORCEN","DERINCMUL","DERMTHD" };
	vector<string> optional_par_group_formal_names{ "SPLITTHRESH","SPLITRELDIFF" };
	vector<string> par_formal_names{ "PARNME","PARTRANS","PARCHGLIM","PARVAL1","PARLBND","PARUBND","PARGP","SCALE","OFFSET","DERCOM" };
	vector<string> par_easy_names{ "NAME","TRANSFORM","CHANGE_LIMIT","VALUE","LOWER_BOUND","UPPER_BOUND","GROUP","SCALE","OFFSET","DERCOM" };
	map <string,string> row_map, temp_tied_map;
	vector<string> obs_formal_names{ "OBSNME","OBSVAL","WEIGHT","OBGNME" };
	vector<string> obs_easy_names{ "NAME","VALUE","WEIGHT","GROUP" };
	vector<string> obs_group_formal_names{ "OBGNME" };
	vector<string> pi_formal_names{ "PILBL","EQUATION","WEIGHT","OBGNME" };
	vector<string> tokens_case_sen;
	vector<string> model_input_formal_names{ "PEST_FILE","MODEL_FILE" };
	vector<string> model_output_formal_names{ "PEST_FILE","MODEL_FILE" };
//#ifndef _DEBUG
//	try {
//#endif
		prior_info_string = "";
		
		for (lnum = 1, sec_begin_lnum = 1; getline(fin, line); ++lnum)
		{
			strip_ip(line);
			line_upper = upper_cp(line);
			tokens.clear();
			tokens_case_sen.clear();
			tokenize(line_upper, tokens);
			tokenize(line, tokens_case_sen);
			sec_lnum = lnum - sec_begin_lnum;

			if (lnum == 1)
			{
				if (tokens[0] != "PCF")
				{
					throw_control_file_error(f_rec, "fist line of control file should be 'PCF' not " + tokens[0]);
				}
			}

			else if (tokens.empty())
			{
				//skip blank line
				lnum--;
			}
			else if (line[0] == '#')
			{
				
				lnum--;
			}
			else if (line_upper.substr(0, 2) == "++")
			{
				if (sections_found.find("CONTROL DATA KEYWORD") != sections_found.end())
					throw_control_file_error(f_rec, "'* control data keyword' can't be used with '++' args");
				sections_found.insert("PLUSPLUS");
				pestpp_input.push_back(line);
				section = "PLUSPLUS";
				if (!prior_info_string.empty())
                {
                    tokens_to_pi_rec(f_rec,prior_info_string);
                    prior_info_string.clear();
                }
			}

			else if (line_upper[0] == '*')
			{
				section = upper_cp(strip_cp(line_upper, "both", " *\t\n"));
				if (sections_found.find(section) != sections_found.end())
				{
					ss.str("");
					ss << "control file error: duplicate entries for section: '" << section << "'";
					throw_control_file_error(f_rec, ss.str());
				}
				sections_found.insert(section);
				if ((nonkeyword_sections.find(section) != nonkeyword_sections.end()) &&
					(sections_found.find("CONTROL DATA KEYWORD") != sections_found.end()))
					
					throw_control_file_error(f_rec,"non-keyword section '" + section + "' not allowed to be used with '* control data keyword'");
				
				if (section == "MODEL INPUT/OUTPUT")
				{
					if (sections_found.find("MODEL INPUT") != sections_found.end())
						throw_control_file_error(f_rec, "'MODEL INPUT/OUTPUT section can't be used with 'MODEL INPUT' section");
					if (sections_found.find("MODEL OUTPUT") != sections_found.end())
						throw_control_file_error(f_rec, "'MODEL INPUT/OUTPUT section can't be used with 'MODEL OUTPUT' section");
				}
				if ((section == "MODEL INPUT") && (sections_found.find("MODEL INPUT/OUTPUT") != sections_found.end()))
					throw_control_file_error(f_rec, "'MODEL INPUT' section can't be used with 'MODEL INPUT/OUTPUT' section");
				if ((section == "MODEL OUTPUT") && (sections_found.find("MODEL INPUT/OUTPUT") != sections_found.end()))
					throw_control_file_error(f_rec, "'MODEL OUTPUT' section can't be used with 'MODEL INPUT/OUTPUT' section");


				sec_begin_lnum = lnum;
			}
			else if (section == "CONTROL DATA KEYWORD")
			{
				pair<string, string> kv = parse_keyword_line(f_rec, line);
				//cout << kv.first << ", " << kv.second << endl;
				//cout << endl;
				//PestppOptions::ARG_STATUS stat = pestpp_options.assign_value_by_key(kv.first, kv.second);
				stat = pestpp_options.assign_value_by_key(kv.first,kv.second);
				check_report_assignment(f_rec, stat, kv.first, kv.second);

				//try to use this as a control data arg
				if (stat == PestppOptions::ARG_STATUS::ARG_NOTFOUND)
				{
					stat = control_info.assign_value_by_key(kv.first,kv.second,f_rec);
					check_report_assignment(f_rec, stat, kv.first, kv.second);
				}
				
				//try to use as an SVD arg
				if (stat == PestppOptions::ARG_STATUS::ARG_NOTFOUND)
				{
					stat = svd_info.assign_value_by_key(kv.first,kv.second);
					check_report_assignment(f_rec, stat, kv.first, kv.second);
				}

				//try to use as a regul arg
				if (stat == PestppOptions::ARG_STATUS::ARG_NOTFOUND) {
					if (regul_scheme_ptr) {
						stat = regul_scheme_ptr->assign_value_by_key(kv.first, kv.second);
						check_report_assignment(f_rec, stat, kv.first, kv.second);

					}
					else
					{
						DynamicRegularization dr;
						stat = dr.assign_value_by_key(kv.first, kv.second);
						check_report_assignment(f_rec, stat, kv.first, kv.second);


					}
				}
				
				//ok, found no home for this line
				if (stat == PestppOptions::ARG_STATUS::ARG_NOTFOUND)
				{
					ss.str("");
					ss << "unrecognized '* control data keyword' key-value pair on line: " << endl << "    '" << line << "'" <<  endl;
					throw_control_file_error(f_rec, ss.str(), false);
					notfound_args.push_back(line);
				}		
			}

			else if (section == "CONTROL DATA")
			{

				if (sec_lnum == 1)
				{
					if (tokens[1] == "REGULARIZATION" || tokens[1] == "REGULARISATION")
					{
						control_info.pestmode = ControlInfo::PestMode::REGUL;

					}
					else if (tokens[1] == "ESTIMATION")
						control_info.pestmode = ControlInfo::PestMode::ESTIMATION;
					else if (tokens[1] == "PARETO")
						control_info.pestmode = ControlInfo::PestMode::PARETO;
					else
					{
						throw_control_file_error(f_rec, "unrecognized 'pestmode': " + tokens[1]);
					}
				}

				else if (sec_lnum == 2)
				{
					convert_ip(tokens[0], num_par);
				}
				else if (sec_lnum == 3)
				{
					convert_ip(tokens[0], num_tpl_file);
					if (tokens.size() >= 5) {
                        try {
                            convert_ip(tokens[4], control_info.numcom);
                        }
                        catch (...)
                        {
                            cout << "WARNING: error parsing '" << tokens[4] <<"' to numcom option...continuing" << endl;
                            control_info.numcom = 0;
                        }
					}
					else {
						control_info.numcom = 0;
					}
					if (tokens.size() >= 6) {
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
					// remove text arguments from the line as these can be specified out of order
					// and PEST++ does not use them
					set<string> remove_tags = { "aui", "auid", "noaui", "senreuse", "nsenreuse", "boundscale", "noboundscale" };
					auto end_iter = std::remove_if(tokens.begin(), tokens.end(),
						[&remove_tags](string& str)->bool {return (remove_tags.find(upper_cp(str)) != remove_tags.end()
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
			else if (section == "PARAMETER GROUPS EXTERNAL")
			{
				pest_utils::ExternalCtlFile efile(line);
				efile.read_file(f_rec);
				
				set<string> cnames = efile.get_col_set();
				for (auto n : par_group_formal_names)
				{
					if (cnames.find(n) == cnames.end())
					{
						ss.str("");
						ss << "external '* parameter group' file '" << efile.get_filename() << "' missing required column '" << n << "'";
						throw_control_file_error(f_rec, ss.str());
					}
				}
				efile.set_index_col_name(par_group_formal_names[0]);
				for (auto oname : optional_par_group_formal_names)
				{
					if (cnames.find(oname) != cnames.end())
						par_group_formal_names.push_back(oname);
				}
				vector<string> par_group_tokens;
				for (auto ro : efile.get_row_order())
				{
					par_group_tokens = efile.get_row_vector(ro, par_group_formal_names);
					tokens_to_par_group_rec(f_rec, par_group_tokens);
				}
				efile.keep_cols(efile_keep_cols);
				if (efiles_map.find(section) == efiles_map.end())
					efiles_map[section] = vector<pest_utils::ExternalCtlFile>{ efile };
				else
					efiles_map[section].push_back(efile);
 
			}
			else if (section == "PARAMETER GROUPS")
			{	
				tokens_to_par_group_rec(f_rec, tokens);				
			}

			else if (section == "PARAMETER DATA")
			{
			if (num_par == -1)
				throw_control_file_error(f_rec,"'* parameter data' section found before '* control data'");
				if (sec_lnum <= num_par) 
				{
					tokens_to_par_rec(f_rec, tokens, t_fixed, t_log, t_scale, t_offset);
				}
				// Get rest of information for tied parameters
				else 
				{
					name = tokens[0];
					string name_tied = tokens[1];
					double numer = ctl_parameters[name];
					double demon = ctl_parameters[name_tied];
					double ratio;
					if (demon == 0.0)
					{
						if (numer == 0.0)
							ratio = 1.0;
						else
						{
							ss.str();
							ss << "tied parameter '" << name << "' is tied to a parameter that has an initial value of 0.0, using tied ratio of 1.0";
							throw_control_file_error(f_rec, ss.str(),false);
							ratio = 1.0;
						}
					}
					else
						ratio = numer / demon;
					t_tied->insert(name, pair<string, double>(name_tied, ratio));
					tied_names.insert(name_tied);
				}
			}

			else if (section == "PARAMETER DATA EXTERNAL")
			{
				
				pest_utils::ExternalCtlFile efile(line);	
				efile.read_file(f_rec);
				set<string> cnames = efile.get_col_set();
				vector<string> get_names;
				//for (auto n : par_formal_names)
				for (int i=0;i<par_formal_names.size();i++)
				{
					string n = par_formal_names[i];
					if (cnames.find(n) == cnames.end())
					{
						string nn = par_easy_names[i];
						if (cnames.find(nn) == cnames.end())
						{
							ss.str("");
							ss << "external '* parameter data' file '" << efile.get_filename() << "' missing required column '";
							ss << n << "' (alias '" + nn + "' also not found";
							throw_control_file_error(f_rec, ss.str());
						}
						else
						{
							get_names.push_back(nn);
							if (i == 0)
								efile.set_index_col_name(par_easy_names[i]);
						}
					}
					else
					{
						get_names.push_back(n);
						if (i == 0)
							efile.set_index_col_name(par_formal_names[i]);
					}
						
				}
				
				string tcol,pcol;
				if (cnames.find("PARTRANS") != cnames.end())
					tcol = "PARTRANS";
				else
					tcol = "TRANSFORM";
                if (cnames.find("PARNME") != cnames.end())
                    pcol = "PARNME";
                else
                    pcol = "NAME";

				vector<string> partrans = efile.get_col_string_vector(tcol);
				set<string> s_partrans(partrans.begin(), partrans.end());
				if (s_partrans.find("TIED") != s_partrans.end())
					if (cnames.find("PARTIED") == cnames.end())
					{
						ss.str("");
						ss << "external '* parameter data' file '" << efile.get_filename() << "' included 'tied' parameters";
						ss << "but doesn't have 'PARTIED' column";
						throw_control_file_error(f_rec, ss.str());
					}

				vector<string> par_tokens;
				for (auto ro : efile.get_row_order())
				{
					par_tokens = efile.get_row_vector(ro, get_names);
					tokens_to_par_rec(f_rec, par_tokens, t_fixed, t_log, t_scale, t_offset);
					//save any tied pars for processing later bc the par its tied to
					//might not have been processed yet.
					row_map = efile.get_row_map(ro);
					if (row_map.at(tcol) == "TIED")
						temp_tied_map[row_map.at(pcol)] = row_map.at("PARTIED");
				}
				efile.keep_cols(efile_keep_cols);
				if (efiles_map.find(section) == efiles_map.end())
					efiles_map[section] = vector<pest_utils::ExternalCtlFile>();
				efiles_map[section].push_back(efile);
	
			}
			else if (section == "OBSERVATION GROUPS")
			{
				tokens_to_obs_group_rec(f_rec, tokens);
			}

			else if (section == "OBSERVATION GROUPS EXTERNAL")
			{
				pest_utils::ExternalCtlFile efile(line);		
				efile.read_file(f_rec);
					
				set<string> cnames = efile.get_col_set();
				for (auto n : obs_group_formal_names)
				{
					if (cnames.find(n) == cnames.end())
					{
						ss.str("");
						ss << "external '* observation group' file '" << efile.get_filename() << "' missing required column '" << n << "'";
						throw_control_file_error(f_rec, ss.str());
					}
				}
				efile.set_index_col_name(obs_group_formal_names[0]);
				vector<string> obs_group_tokens;
				for (auto ro : efile.get_row_order())
				{
					obs_group_tokens = efile.get_row_vector(ro, obs_group_formal_names);
					tokens_to_obs_group_rec(f_rec, obs_group_tokens);
				}
				efile.keep_cols(efile_keep_cols);
				if (efiles_map.find(section) == efiles_map.end())
				{
					efiles_map[section] = vector<pest_utils::ExternalCtlFile>();
				}	
				efiles_map[section].push_back(efile);
			}

			else if (section == "OBSERVATION DATA")
			{
				tokens_to_obs_rec(f_rec, tokens);
			}

			else if (section == "OBSERVATION DATA EXTERNAL")
			{
				
				pest_utils::ExternalCtlFile efile(line);		
				efile.read_file(f_rec);
					
				vector<string> get_names;
				set<string> cnames = efile.get_col_set();
				//for (auto n : obs_formal_names)
				for (int i = 0; i < obs_formal_names.size(); i++)
				{
					string n = obs_formal_names[i];
					if (cnames.find(n) == cnames.end())
					{
						string nn = obs_easy_names[i];
						if (cnames.find(nn) == cnames.end())
						{
							ss.str("");
							ss << "external '* observation data' file '" << efile.get_filename() << "' missing required column '" << n << "'";
							throw_control_file_error(f_rec, ss.str());
						}
						else
						{
							get_names.push_back(nn);
							if (i == 0)
								efile.set_index_col_name(obs_easy_names[i]);
						}

					}
					else
					{
						get_names.push_back(n);
						if (i == 0)
							efile.set_index_col_name(obs_formal_names[i]);
					}

				}

				vector<string> obs_tokens;
				for (auto ro : efile.get_row_order())
				{
					obs_tokens = efile.get_row_vector(ro, get_names);
					tokens_to_obs_rec(f_rec, obs_tokens);
				}
				efile.keep_cols(efile_keep_cols);
				if (efiles_map.find(section) == efiles_map.end())
					efiles_map[section] = vector<pest_utils::ExternalCtlFile>();
				efiles_map[section].push_back(efile);
			}

			else if (section == "PRIOR INFORMATION")
			{
				tokens_to_pi_rec(f_rec, line_upper);
			}

			else if (section == "PRIOR INFORMATION EXTERNAL")

			{
				
				pest_utils::ExternalCtlFile efile(line);
				efile.read_file(f_rec);
				set<string> cnames = efile.get_col_set();
				for (auto n : pi_formal_names)
				{
					if (cnames.find(n) == cnames.end())
					{
						ss.str("");
						ss << "external '* prior information' file '" << efile.get_filename() << "' missing required column '" << n << "'";
						throw_control_file_error(f_rec, ss.str());
					}
				}
				efile.set_index_col_name(pi_formal_names[0]);
				vector<string> pi_tokens;
				for (auto ro : efile.get_row_order())
				{
					pi_tokens = efile.get_row_vector(ro, pi_formal_names);
					tokens_to_pi_rec(f_rec, pi_tokens);
				}
				efile.keep_cols(efile_keep_cols);
				if (efiles_map.find(section) == efiles_map.end())
					efiles_map[section] = vector<pest_utils::ExternalCtlFile>{ efile };
				else
					efiles_map[section].push_back(efile);	
			}

		
			else if (section == "MODEL COMMAND LINE")
			{
                if ((line.find('\"') != std::string::npos) || (line.find('\'') != std::string::npos))
                {
                    ss.str("");
                    ss << "WARNING: single and/or double quote char(s) found in model command line :" << line << endl;

                    string temp_line = line;
                    temp_line.erase(std::remove_if(temp_line.begin(), temp_line.end(), IsQuote), temp_line.end());
                    //pest_utils::strip_ip(temp_line);

                    ss << "         new model command line: " << temp_line << endl;
                    cout << ss.str();
                    f_rec << ss.str();
                    model_exec_info.comline_vec.push_back(string(temp_line));
                }
                else
                {
                    model_exec_info.comline_vec.push_back(line);
                }
			}

			else if (section == "MODEL INPUT")
			{
				if (tokens.size() != 2)
					throw_control_file_error(f_rec, "wrong number of tokens on '* model input' line '" + line + "' expecting 2");
				for (auto& token : tokens_case_sen)
                {
                    if ((token.find('\"') != std::string::npos) || (token.find('\'') != std::string::npos))
                    {
                        ss.str("");
                        ss << "WARNING: single and/or double quote char(s) found in model interface file: " << token << endl;
                        cout << ss.str();
                        f_rec << ss.str();
                        string temp_line = token;
                        temp_line.erase(std::remove_if(temp_line.begin(), temp_line.end(), IsQuote), temp_line.end());
                        token = string(temp_line);
                    }
                }


                model_exec_info.tplfile_vec.push_back(tokens_case_sen[0]);
				model_exec_info.inpfile_vec.push_back(tokens_case_sen[1]);
			}

			else if (section == "MODEL INPUT EXTERNAL")
			{
				
				pest_utils::ExternalCtlFile efile(line, false);	
				efile.read_file(f_rec);
					
				set<string> cnames = efile.get_col_set();
				for (auto n : model_input_formal_names)
				{
					if (cnames.find(n) == cnames.end())
					{
						ss.str("");
						ss << "external '* model input' file '" << efile.get_filename() << "' missing required column '" << n << "'";
						throw_control_file_error(f_rec, ss.str());
					}
				}
				efile.set_index_col_name(model_input_formal_names[0]);
				vector<string> mi_tokens;
				for (auto ro : efile.get_row_order())
				{
					mi_tokens = efile.get_row_vector(ro, model_input_formal_names);
                    for (auto& token: mi_tokens) {
                        if ((token.find('\"') != std::string::npos) || (token.find('\'') != std::string::npos)) {
                            ss.str("");
                            ss << "WARNING: single and/or double quote char(s) found in model interface file: " << token
                               << endl;
                            cout << ss.str();
                            f_rec << ss.str();
                            string temp_line = token;
                            temp_line.erase(std::remove_if(temp_line.begin(), temp_line.end(), IsQuote),
                                            temp_line.end());
                            token = string(temp_line);
                        }
                    }
					model_exec_info.tplfile_vec.push_back(mi_tokens[0]);
					model_exec_info.inpfile_vec.push_back(mi_tokens[1]);
					//model_exec_info.incycle_vec.push_back(std::stoi(mi_tokens[2]));
				}
				efile.keep_cols(efile_keep_cols);
				if (efiles_map.find(section) == efiles_map.end())
					efiles_map[section] = vector<pest_utils::ExternalCtlFile>{ efile };
				else
					efiles_map[section].push_back(efile);	
			
			}

			else if (section == "MODEL OUTPUT")
			{
				if (tokens.size() != 2)
					throw_control_file_error(f_rec, "wrong number of tokens on '* model output' line '" + line + "' expecting 2");
                for (auto& token : tokens_case_sen) {
                    if ((token.find('\"') != std::string::npos) || (token.find('\'') != std::string::npos)) {
                        ss.str("");
                        ss << "WARNING: single and/or double quote char(s) found in model interface file: " << token
                           << endl;
                        cout << ss.str();
                        f_rec << ss.str();
                        string temp_line = token;
                        temp_line.erase(std::remove_if(temp_line.begin(), temp_line.end(), IsQuote), temp_line.end());
                        token = string(temp_line);
                    }
                }
				model_exec_info.insfile_vec.push_back(tokens_case_sen[0]);
				model_exec_info.outfile_vec.push_back(tokens_case_sen[1]);
			}


			else if (section == "MODEL OUTPUT EXTERNAL")
			{
				pest_utils::ExternalCtlFile efile(line, false);	
				efile.read_file(f_rec);
				set<string> cnames = efile.get_col_set();
				for (auto n : model_output_formal_names)
				{
					if (cnames.find(n) == cnames.end())
					{
						ss.str("");
						ss << "external '* model output' file '" << efile.get_filename() << "' missing required column '" << n << "'";
						throw_control_file_error(f_rec, ss.str());
					}
				}
				efile.set_index_col_name(model_output_formal_names[0]);
				vector<string> mo_tokens;
				for (auto ro : efile.get_row_order())
				{
					mo_tokens = efile.get_row_vector(ro, model_output_formal_names);
                    for (auto& token : mo_tokens)
                    {
                        if ((token.find('\"') != std::string::npos) || (token.find('\'') != std::string::npos))
                        {
                            ss.str("");
                            ss << "WARNING: single and/or double quote char(s) found in model interface file: " << token << endl;
                            cout << ss.str();
                            f_rec << ss.str();
                            string temp_line = token;
                            temp_line.erase(std::remove_if(temp_line.begin(), temp_line.end(), IsQuote), temp_line.end());
                            token = string(temp_line);
                        }
                    }
					model_exec_info.insfile_vec.push_back(mo_tokens[0]);
					model_exec_info.outfile_vec.push_back(mo_tokens[1]);
					//model_exec_info.outcycle_vec.push_back(std::stoi(mo_tokens[2]));
					
				}
				efile.keep_cols(efile_keep_cols);
				if (efiles_map.find(section) == efiles_map.end())
					efiles_map[section] = vector<pest_utils::ExternalCtlFile>{ efile };
				else
					efiles_map[section].push_back(efile);
			
			}
			else if (section == "MODEL INPUT/OUTPUT")
			{
			if (tokens.size() != 2)
				throw_control_file_error(f_rec, "wrong number of tokens on '* model input/output' line '" + line + "' expecting 2");
                for (auto& token : tokens_case_sen)
                {
                    if ((token.find('\"') != std::string::npos) || (token.find('\'') != std::string::npos))
                    {
                        ss.str("");
                        ss << "WARNING: single and/or double quote char(s) found in model interface file: " << token << endl;
                        cout << ss.str();
                        f_rec << ss.str();
                        string temp_line = token;
                        temp_line.erase(std::remove_if(temp_line.begin(), temp_line.end(), IsQuote), temp_line.end());
                        token = string(temp_line);
                    }
                }
				if (i_tpl_ins < num_tpl_file)
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
			else if (section == "REGULARISATION" || section == "REGULARIZATION")
			{
				if (sec_lnum == 1) {
					if (regul_scheme_ptr)
					{
						regul_scheme_ptr->assign_value_by_key("PHIMLIM", tokens[0]);
						regul_scheme_ptr->assign_value_by_key("PHIMACCEPT", tokens[1]);
					}
					if (tokens.size() >= 3)
						if (regul_scheme_ptr)
							regul_scheme_ptr->assign_value_by_key("FRACPHIM", tokens[2]);
					/*convert_ip(tokens[0], phimlim);
					convert_ip(tokens[1], phimaccept);
					fracphim = 0.0;
					if (tokens.size() >= 3) convert_ip(tokens[2], fracphim);*/
				}
				else if (sec_lnum == 2) {
					/*convert_ip(tokens[0], wfinit);
					convert_ip(tokens[1], wfmin);
					convert_ip(tokens[2], wfmax);*/
					if (regul_scheme_ptr)
					{
						regul_scheme_ptr->assign_value_by_key("WFINIT", tokens[0]);
						regul_scheme_ptr->assign_value_by_key("WFMIN", tokens[1]);
						regul_scheme_ptr->assign_value_by_key("WFMAX", tokens[2]);
					}
				}
				else if (sec_lnum == 3) {
					int iregadj;
					/*convert_ip(tokens[0], wffac);
					convert_ip(tokens[1], wftol);*/
					if (regul_scheme_ptr)
					{
						regul_scheme_ptr->assign_value_by_key("WFFAC", tokens[0]);
						regul_scheme_ptr->assign_value_by_key("WFTOL", tokens[1]);
					}

					if (tokens.size() > 2)
					{
						if (regul_scheme_ptr)
							regul_scheme_ptr->assign_value_by_key("IREGADJ", tokens[2]);
					}
					/*delete regul_scheme_ptr;
					regul_scheme_ptr = new DynamicRegularization(use_dynamic_reg, reg_adj_grp_weights, phimlim,
						phimaccept, fracphim, wfmin, wfmax, wffac, wftol, wfinit);*/
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
			tokens_to_pi_rec(f_rec, line_upper);
		}
//#ifndef _DEBUG
//	}
//	catch (PestConversionError& e) {
//		std::stringstream out;
//		out << "Error processing \"" << pst_filename << "\" on line number " << lnum << endl;
//		out << e.what() << endl;
//		e.add_front(out.str());
//		e.raise();
//	}
//#endif
	fin.close();

	
	
	// handle any tied pars found in external files
	double numer, demon, ratio;
	vector<string> missing;
	for (auto p: temp_tied_map)
	{
		name = p.first;
		string name_tied = p.second;
		numer = ctl_parameters[name];
//		if (ctl_parameters.find(name_tied) == ctl_parameters.end())
//        {
//		    missing.push_back(name_tied);
//		    continue;
//        }
		demon = ctl_parameters[name_tied];
		if (demon == 0.0)
		{
			if (numer == 0.0)
				ratio = 1.0;
			else
			{
				ss.str("");
				ss << "tied parameter '" << name << "' has an initial value of 0.0.  Using a tied ratio of 1.0";
				throw_control_file_error(f_rec, ss.str(),false);
				ratio = 1.0;
			}
		}
		else
			ratio = numer / demon;
		t_tied->insert(name, pair<string, double>(name_tied, ratio));
		tied_names.insert(name_tied);
	}

	if (missing.size() > 0)
    {
	    ss.str("");
	    ss << "Error: the following `partied` parameters were not found in the control file:";
	    for (auto& m: missing)
	        ss << m << ",";
	    f_rec << ss.str() << endl;
	    throw runtime_error(ss.str());
    }

	//process pestpp options
	map<string, PestppOptions::ARG_STATUS> arg_map, line_arg_map;
	vector<string> dup;
	for (vector<string>::const_iterator b = pestpp_input.begin(), e = pestpp_input.end();
		b != e; ++b) {

		try
		{
			line_arg_map = pestpp_options.parse_plusplus_line(*b);
		}
		catch (exception &e)
		{
			throw runtime_error("error parsing '++' line :'" + *b + "': "+ e.what());
		}
		catch (...)
		{
			throw runtime_error("error parsing '++' line :'" + *b + "'");
		}
		for (auto arg : line_arg_map)
		{
			if (arg_map.find(arg.first) != arg_map.end())
				dup.push_back(arg.first);

		}
		arg_map.insert(line_arg_map.begin(), line_arg_map.end());
	}
	pestpp_options.rectify_ies_da_args();

	if (dup.size() > 0)
	{
		ss.str("");
		ss << " the following '++' args are duplicates (possibly thru an alias):" << endl;
		for (auto n : dup)
			ss << n << ",";
		throw_control_file_error(f_rec, ss.str());
	}

	vector<string> invalid;
	for (auto kv : arg_map)
	{
		if (kv.second == PestppOptions::ARG_STATUS::ARG_INVALID)
		{
			invalid.push_back(kv.first);
		}
	}
	if (invalid.size() > 0)
	{
		ss.str("");
		ss << " the following '++' args have invalid values:" << endl;
		for (auto n : invalid)
			ss << n << ",";
		throw_control_file_error(f_rec, ss.str());
	}

	vector<string> not_accepted;
	for (auto kv : arg_map)
	{
		if (kv.second != PestppOptions::ARG_STATUS::ARG_ACCEPTED)
		{
			not_accepted.push_back(kv.first);
		}
	}

	if (not_accepted.size() > 0)
	{
		ss.str("");
		ss << " the following '++' args were not accepted:" << endl;
		for (auto n : not_accepted)
			ss << n << ",";
		if (pestpp_options.get_forgive_unknown_args())
		{
			ss << endl << "forgive_unknown_args is 'true', continuing" << endl;
			throw_control_file_error(f_rec, ss.str(), false);
		}
		else
		{
			ss << endl << "forgive_unknown_args is 'false' so this is treated as an error" << endl;
			throw_control_file_error(f_rec, ss.str(), true);
		}
	}
	if (notfound_args.size() > 0)
	{
		ss.str("");
		ss << " the following control data keyword lines were not accepted:" << endl;
		for (auto n : notfound_args)
			ss << n << ",";
		if (pestpp_options.get_forgive_unknown_args())
		{
			ss << endl << "forgive_unknown_args is 'true', continuing" << endl;
			throw_control_file_error(f_rec, ss.str(), false);
		}
		else
		{
			ss << endl << "forgive_unknown_args is 'false' so this is treated as an error" << endl;
			throw_control_file_error(f_rec, ss.str(), true);
		}
	}
	if (control_info.pestmode == ControlInfo::PestMode::REGUL)
	{
		if (regul_scheme_ptr)
			regul_scheme_ptr->set_use_dynamic_reg(true);
	}
	
	if (ctl_ordered_obs_group_names.size() == 0)
	{
		set<string> found;
		string gname;
		for (auto& _name : ctl_ordered_obs_names)
		{
			gname = observation_info.get_group(_name);
			if (found.find(gname) == found.end())
			{
				found.insert(gname);
				ctl_ordered_obs_group_names.push_back(gname);
			}
		}
	}
	if (ctl_ordered_par_group_names.size() == 0)
	{
		ctl_ordered_par_group_names = base_group_info.get_group_names();
	}



	//check if the predictions ++ arg might be a file name?
	vector<string> pred_arg = pestpp_options.get_prediction_names();
	if ((pestpp_options.get_uncert_flag()) && (pred_arg.size() == 1))
	{
		string fname = pred_arg[0];
		if (pest_utils::check_exist_in(fname))
		{
			f_rec << "filename '" << fname << "' detected for prediction names, reading...";
			vector<string> pred_names = pest_utils::read_onecol_ascii_to_vector(fname);
			f_rec << pred_names.size() << " predictions found" << endl;
			pestpp_options.set_prediction_names(pred_names);
		}
	}

	//since par groups are optional, make sure every par has a group...
	rectify_par_groups();

	if (regul_scheme_ptr)
		regul_scheme_ptr->set_max_reg_iter(pestpp_options.get_max_reg_iter());

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
				ratio = ctl_parameters[pname] / ctl_parameters[tie_to_name];
				t_tied->insert(pname, pair<string, double>(tie_to_name, ratio));
				ctl_parameter_info.get_parameter_rec_ptr_4_mod(pname)->tranform_type = ttied;
			}

		}
		f_rec << "-->number of adjustable parameters reduced from " << n_adj_par << " to " << new_n_adj_par << endl;
		cout << "-->number of adjustable parameters reduced from " << n_adj_par << " to " << new_n_adj_par << endl;
		n_adj_par = new_n_adj_par;
		if (tied_names.size() > 0)
		{
			f_rec << "-->the following existing adjustable parameters that others tied to it have been maintained:" << endl;
			for (auto tname : tied_names)
				f_rec << tname << endl;
		}
	}

	//clear temp containers
	s_pargp.clear();
	s_obgnme.clear();
	

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

pair<string,double> Pest::enforce_par_limits(PerformanceLog* performance_log, Parameters & upgrade_active_ctl_pars, const Parameters &last_active_ctl_pars, bool enforce_chglim, bool enforce_bounds)
{
	if ((!enforce_chglim) && (!enforce_bounds)) {
        pair<string, double> _control_info("no enforcement", 1.0);
        return _control_info;
    }
	stringstream ss;
	double fpm = control_info.facparmax;
	double facorig = control_info.facorig;
	double rpm = control_info.relparmax;
	double orig_val, last_val, fac_lb, fac_ub, rel_lb, rel_ub, eff_ub, eff_lb,chg_lb, chg_ub;
	double chg_fac, chg_rel;
	string parchglim;
	double scaling_factor = 1.0;
	double temp = 1.0;
	double bnd_tol = 0.001;
	double scaled_bnd_val;
	string controlling_par = "";
	string control_type = "";
	ParameterInfo &p_info = ctl_parameter_info;
	const ParameterRec *p_rec;
	Parameters upgrade_ctl_pars;
	Parameters last_ctl_pars;

	// if tied parameters exist, use the old scaling_factor code
	// this ensures compliance with tied parameters.
	if (pestpp_options.get_enforce_tied_bounds())
	{
		upgrade_ctl_pars = base_par_transform.active_ctl2ctl_cp(upgrade_active_ctl_pars);
		last_ctl_pars = base_par_transform.active_ctl2ctl_cp(last_active_ctl_pars);
		for (auto& p : upgrade_ctl_pars)
		{
			/*if (pest_utils::lower_cp(p.first) == "s_xomehgwat")
				cout << p.first << endl;*/
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
			//if (p.second > 0.0)
			if (last_val > 0.0)
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
						control_type = "upper change limit";
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
					control_type = "lower change limit";
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
				scaled_bnd_val = p_rec->ubnd + abs(p_rec->ubnd * bnd_tol);
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
						control_type = "upper bound";
					}
				}
				scaled_bnd_val = p_rec->lbnd - abs(p_rec->lbnd * bnd_tol);
				if (p.second < p_rec->lbnd)
				{
					temp = abs((last_val - p_rec->lbnd) / (last_val - p.second));
					if ((temp > 1.0) || (temp < 0.0))
					{
						ss.str("");
						ss << "Pest::enforce_par_limts() error: invalid lower bound scaling factor " << temp << " for par " << p.first << endl;
						ss << " lbnd:" << p_rec->lbnd << ", last_val:" << last_val << ", current_val:" << p.second << endl;
						throw runtime_error(ss.str());
					}
					if (temp < scaling_factor)
					{
						scaling_factor = temp;
						controlling_par = p.first;
						control_type = "lower bound";
					}
				}
			}	
		}
		ss.str("");
		ss << "change enforcement controlling par:" << controlling_par << ", control_type: " << control_type << ", scaling_factor: " << scaling_factor << endl;

		if (scaling_factor == 0.0)
		{
			ss.str("");
			ss << "Pest::enforce_par_change_limits error : zero length parameter vector" << endl;
			ss << "parameter: " << controlling_par << ", control type: " << control_type;
			throw runtime_error(ss.str());
		}

		if (scaling_factor != 1.0)
		{
			for (auto &p : upgrade_active_ctl_pars)
			{
				
				last_val = last_ctl_pars.get_rec(p.first);
				p.second =last_val + (p.second - last_val) *  scaling_factor;
			}
		}
		
		//check for slightly out of bounds
		for (auto &p : upgrade_ctl_pars)
		{
			p_rec = p_info.get_parameter_rec_ptr(p.first);
			if (p.second < p_rec->lbnd)
				p.second = p_rec->lbnd;
			else if (p.second > p_rec->ubnd)
				p.second = p_rec->ubnd;

		}
	}
	// if we don't have tied parameters, use clamping logic instead.
	else
	{
		for (auto &p : upgrade_active_ctl_pars)
		{
			
			last_val = last_active_ctl_pars.get_rec(p.first);
			p_rec = p_info.get_parameter_rec_ptr(p.first);
			parchglim = p_rec->chglim;
			if (parchglim == "RELATIVE" && last_val == 0.0)
			{
				throw runtime_error("Relative parchglim not defined for zero-valued parameter " + p.first);
			}

			

			if (p.second == 0.0)
				p.second = p_rec->ubnd / 4.0;
			orig_val = ctl_parameters.get_rec(p.first);
			if (orig_val == 0.0)
				orig_val = p_rec->ubnd / 4.0;

			//calc fac lims
			if (last_val > 0.0)
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
			rel_lb = last_active_ctl_pars.get_rec(p.first) - (abs(last_val) * rpm);
			rel_ub = last_active_ctl_pars.get_rec(p.first) + (abs(last_val) * rpm);

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


			// double temp = 1.0;

			// New logic for enforcing chglim and bounds
			// First, we'll check the change limits.
			// Next, we'll check parameter bounds.
			// If anything violates, clamp to the offending bound.

			if (enforce_chglim)
			// similar to below, clamp rather than shrink every parameter if a parameter violates change limits
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
					p.second = chg_ub;
				}
				else if (p.second < chg_lb)
				{
					temp = abs((last_val - chg_lb) / (last_val - p.second));
					if ((temp > 1.0) || (temp < 0.0))
					{
						ss.str("");
						ss << "Pest::enforce_par_limts() error: invalid lower parchglim scaling factor " << temp << " for par " << p.first << endl;
						ss << " chglim:" << chg_lb << ", last_val:" << last_val << ", current_val:" << p.second << endl;
						throw runtime_error(ss.str());
					}
					p.second = chg_lb;
				}

				//apply facorig correction if needed
				if (ctl_parameter_info.get_parameter_rec_ptr(p.first)->tranform_type == ParameterRec::TRAN_TYPE::NONE)
				{
					if (abs(p.second) < abs(orig_val) * facorig)
						p.second = orig_val * facorig;
					if (abs(last_val) < abs(orig_val * facorig))
						last_val = orig_val * facorig;
				}
			}

			if (enforce_bounds)
			{
				if (p.second > p_rec->ubnd)
				{
					p.second = p_rec->ubnd;
				}
				else if (p.second < p_rec->lbnd)
				{
					p.second = p_rec->lbnd;
				}
			}	
		}
	}

	pair<string, double> _control_info(ss.str(), scaling_factor);
	return _control_info;
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

Pest Pest::get_child_pest(int icycle)
{
	// TODO: insert return statement here

	Pest child_pest = Pest(*this);
	child_pest.child_pest_update(icycle);
	return child_pest;
//	Pest* child_pest = new Pest(*this);
//	child_pest->child_pest_update(icycle);
//	return *child_pest;
	
}


void Pest::assign_da_cycles(ofstream &f_rec)
{
	vector<string> str_values, str_names;
	stringstream ss, ss_tpl;
	set<string> col_names;
	string name_col;
	int cycle;
	if (efiles_map.find("PARAMETER DATA EXTERNAL") == efiles_map.end())
	{
		throw_control_file_error(f_rec, "could not find 'parameter data external' section for cycle info, all parameter quantities being assigned 'cycle'=-1", false);
//		for (auto pname : get_ctl_ordered_adj_par_names())
//		{
//			ctl_parameter_info.get_parameter_rec_ptr_4_mod(pname)->cycle = -1;
//		}
	}

	else
	{
        vector<pair<string, DaCycleInfo>> par_cycle_dci_map = extract_cycle_info(f_rec, "PARAMETER DATA EXTERNAL", vector<string>{"PARNME", "NAME"});
//		vector<pair<string, int>> par_cycle_map = extract_cycle_numbers2(f_rec, "PARAMETER DATA EXTERNAL", vector<string>{"PARNME", "NAME"});
//		if (par_cycle_map.size() == 0)
//		{
//			throw_control_file_error(f_rec, "could not find cycle info in external file(s), all parameter quantities being assigned 'cycle'=-1", false);
//			for (auto pname : get_ctl_ordered_adj_par_names())
//			{
//				ctl_parameter_info.get_parameter_rec_ptr_4_mod(pname)->cycle = -1;
//			}
//
//		}
//		else
		{
//			for (auto pc : par_cycle_map)
//			{
//				ctl_parameter_info.get_parameter_rec_ptr_4_mod(pc.first)->cycle = pc.second;
//			}
			for (auto& pc : par_cycle_dci_map) {
                ctl_parameter_info.get_parameter_rec_ptr_4_mod(pc.first)->dci = pc.second;
            }

		}

	}
		
	

	if (efiles_map.find("OBSERVATION DATA EXTERNAL") == efiles_map.end())
	{
		throw_control_file_error(f_rec, "could not find 'observation data external' section, assigning all observations to cycle '0'", false);
//		for (auto name : get_ctl_ordered_nz_obs_names())
//		{
//			observation_info.get_observation_rec_ptr_4_mod(name)->cycle = 0;
//		}
	}

	else
	{
        vector<pair<string, DaCycleInfo>> obs_cycle_dci_map = extract_cycle_info(f_rec, "OBSERVATION DATA EXTERNAL", vector<string>{"OBSNME", "NAME"});

        //vector<pair<string, int>> obs_cycle_map = extract_cycle_numbers2(f_rec, "OBSERVATION DATA EXTERNAL", vector<string>{"OBSNME", "NAME"});
		if (obs_cycle_dci_map.size() == 0)
		{
			throw_control_file_error(f_rec, "no observation cycle information was found in external file(s), assigning all observations to cycle '0'", false);

		}
		vector<string> missing, obs_vec;
		for (auto pp: obs_cycle_dci_map)
		{
			obs_vec.push_back(pp.first);
		}
		for (auto name : get_ctl_ordered_nz_obs_names())
		{
			if (std::find(obs_vec.begin(), obs_vec.end(), name) == obs_vec.end())			
				missing.push_back(name);
		}
		if (missing.size() > 0)
		{
			ss.str("");
			ss << "the following non-zero weighted observations do not have cycle information:";
			for (auto m : missing)
				ss << m << ",";
			throw_control_file_error(f_rec, ss.str());
		}
//		for (auto oc : obs_cycle_map) {
//            observation_info.get_observation_rec_ptr_4_mod(oc.first)->cycle = oc.second;
//        }
		for (auto& oc : obs_cycle_dci_map)
        {
            observation_info.get_observation_rec_ptr_4_mod(oc.first)->dci = oc.second;
        }
	}

	model_exec_info.incycle_dci_vec.clear();
	if (efiles_map.find("MODEL INPUT EXTERNAL") == efiles_map.end())
	{
		throw_control_file_error(f_rec, "could not find 'model input external' section, all template/model input files will be used every cycle", false);
		for (auto tpl : model_exec_info.tplfile_vec)
		{
			//model_exec_info.incycle_vec.push_back(-1);
			model_exec_info.incycle_dci_vec.push_back(DaCycleInfo());
		}
	}
	else
	{
		//vector<pair<string, int>> mi_cycle_map = extract_cycle_numbers2(f_rec, "MODEL INPUT EXTERNAL", vector<string>{"PEST_FILE"});
        vector<pair<string, DaCycleInfo>> mi_cycle_dci_map = extract_cycle_info(f_rec, "MODEL INPUT EXTERNAL", vector<string>{"PEST_FILE"});
		if (mi_cycle_dci_map.size() == 0)
		{
			throw_control_file_error(f_rec, "could not find cycle info in external file(s), all template/model input files will be used every cycle", false);
			for (auto tpl : model_exec_info.tplfile_vec)
			{
				//model_exec_info.incycle_vec.push_back(-1);
				model_exec_info.incycle_dci_vec.push_back(DaCycleInfo());

			}
		}	
		else
		{
			
//			for (auto tpl : mi_cycle_map)
//			{
//				model_exec_info.incycle_vec.push_back(tpl.second); // we can do that because row order does not change.
//
//			}
			for (auto& tpl : mi_cycle_dci_map)
            {
                model_exec_info.incycle_dci_vec.push_back(tpl.second);
            }
		}
	}


	model_exec_info.outcycle_dci_vec.clear();
	if (efiles_map.find("MODEL OUTPUT EXTERNAL") == efiles_map.end())
	{
		throw_control_file_error(f_rec, "could not find 'model output external' section, assigning all instruction/out files to cycle '0'", false);
		for (auto ins : model_exec_info.insfile_vec)
		{
			//model_exec_info.outcycle_vec.push_back(0);
			model_exec_info.outcycle_dci_vec.push_back(DaCycleInfo());
		}
	}
	else
	{
		//vector<pair<string, int>> mi_cycle_map = extract_cycle_numbers2(f_rec, "MODEL OUTPUT EXTERNAL", vector<string>{"PEST_FILE"});
        vector<pair<string, DaCycleInfo>> mi_cycle_dci_map = extract_cycle_info(f_rec, "MODEL OUTPUT EXTERNAL", vector<string>{"PEST_FILE"});
		if (mi_cycle_dci_map.size() == 0)
		{
			throw_control_file_error(f_rec, "no model output cycle information was found in external file(s), assigning all instruction/out files to cycle '0'",false);
			for (auto ins : model_exec_info.insfile_vec)
			{
				//model_exec_info.outcycle_vec.push_back(0);
                model_exec_info.outcycle_dci_vec.push_back(DaCycleInfo());
			}
		}
		else
		{		
//			for (auto ins : mi_cycle_map)
//			{
//				model_exec_info.outcycle_vec.push_back(ins.second); // we can do that because row order does not change.
//			}
			for (auto& ins : mi_cycle_dci_map)
            {
			    model_exec_info.outcycle_dci_vec.push_back(ins.second);
            }
		}
	}
	//TODO: prior info...with pi, need to check that parameters in each eqs are in the same cycle!
}

DaCycleInfo Pest::parse_cycle_str(string& raw_cycle_val, string& efilename, int row, ofstream& f_rec)
{
    stringstream ss;
    DaCycleInfo dci;
    dci.start = 0;
    dci.stop = -999;
    dci.stride = 1;
    string sub_str;
    int cycle;
    int idx;
    if (raw_cycle_val.find(':') != string::npos)
    {
        //no start
        if (raw_cycle_val[0] == ':')
        {
            //no stop, so just stride
            if (raw_cycle_val[1] == ':')
            {
                idx = raw_cycle_val.find_last_of(':');
                sub_str = raw_cycle_val.substr(idx+1,raw_cycle_val.size());
                try {
                    dci.stride = stoi(sub_str);
                }
                catch (...) {
                    ss.str("");
                    ss << "error casting cycle stride '" << sub_str << "' to int for cycle info string '" << raw_cycle_val << "' on row " << row << "of external file "
                       << efilename << " , Stopped...";
                    throw_control_file_error(f_rec, ss.str());
                }
            }

            else
            {
                //parse stop
                sub_str = raw_cycle_val.substr(1,raw_cycle_val.size());
                idx = sub_str.find_first_of(':');
                try {
                    dci.stop = stoi(sub_str.substr(0,idx));

                }
                catch (...) {
                    ss.str("");
                    ss << "error casting cycle stop '" << sub_str.substr(0,idx) << "' to int for cycle info string '" << raw_cycle_val << "' on row " << row << "of external file "
                       << efilename << " , Stopped...";
                    throw_control_file_error(f_rec, ss.str());
                }

                sub_str = sub_str.substr(idx,sub_str.size());
                idx = sub_str.find_first_of(':');
                //a stride too
                if (idx != string::npos)
                {
                    sub_str = sub_str.substr(idx+1,sub_str.size());
                    try {
                        dci.stride = stoi(sub_str);

                    }
                    catch (...) {
                        ss.str("");
                        ss << "error casting cycle stride '" << sub_str << "' to int for cycle info string '" << raw_cycle_val << "' on row " << row << "of external file "
                           << efilename << " , Stopped...";
                        throw_control_file_error(f_rec, ss.str());
                    }
                }
            }
        }
        else
        {
            //parse the start
            idx = raw_cycle_val.find_first_of(':');
            sub_str = raw_cycle_val.substr(0,idx);

            try {
                dci.start = stoi(sub_str);

            }
            catch (...) {
                ss.str("");
                ss << "error casting cycle start '" << sub_str << "' to int for cycle info string '" << raw_cycle_val << "' on row " << row << "of external file "
                   << efilename << " , Stopped...";
                throw_control_file_error(f_rec, ss.str());
            }

            sub_str = raw_cycle_val.substr(idx+1,raw_cycle_val.size());

            //no stop but a stride
            if (sub_str[0] == ':')
            {
                try {
                    dci.stride = stoi(sub_str.substr(1,sub_str.size()));

                }
                catch (...) {
                    ss.str("");
                    ss << "error casting cycle stride '" << sub_str.substr(1,sub_str.size()) << "' to int for cycle info string '" << raw_cycle_val << "' on row " << row << "of external file "
                       << efilename << " , Stopped...";
                    throw_control_file_error(f_rec, ss.str());
                }
            }
            else {
                idx = sub_str.find_first_of(':');

                if (idx != string::npos) {

                    try {
                        dci.stop = stoi(sub_str.substr(0, idx));

                    }
                    catch (...) {
                        ss.str("");
                        ss << "error casting cycle stop '" << sub_str.substr(0, idx)
                           << "' to int for cycle info string '" << raw_cycle_val << "' on row " << row
                           << "of external file "
                           << efilename << " , Stopped...";
                        throw_control_file_error(f_rec, ss.str());
                    }
                    sub_str = sub_str.substr(idx + 1, sub_str.size());
                    try {
                        dci.stride = stoi(sub_str);

                    }
                    catch (...) {
                        ss.str("");
                        ss << "error casting cycle stride '" << sub_str << "' to int for cycle info string '"
                           << raw_cycle_val << "' on row " << row << "of external file "
                           << efilename << " , Stopped...";
                        throw_control_file_error(f_rec, ss.str());
                    }

                } else {
                    try {
                        dci.stop = stoi(sub_str);

                    }
                    catch (...) {
                        ss.str("");
                        ss << "error casting cycle stop '" << sub_str << "' to int for cycle info string '"
                           << raw_cycle_val << "' on row " << row << "of external file "
                           << efilename << " , Stopped...";
                        throw_control_file_error(f_rec, ss.str());
                    }
                }
            }

        }
    }
    else {
        try {
            cycle = stoi(raw_cycle_val);

        }
        catch (...) {
            ss.str("");
            ss << "error casting cycle '" << raw_cycle_val << "' to int on row " << row << "of external file "
               << efilename << " , Stopped...";
            throw_control_file_error(f_rec, ss.str());
        }
        if (cycle != -1) {
            dci.start = cycle;
            dci.stop = cycle;
            dci.stride = 1;
        }
    }
    return dci;
}

vector<pair<string, DaCycleInfo>> Pest::extract_cycle_info(ofstream& f_rec, string section_name, vector<string> possible_name_cols)
{
    vector<string> str_values, str_names;
    stringstream ss;
    set<string> col_names;
    string name_col = "";
    vector<pair<string, DaCycleInfo>> cycle_map;
    DaCycleInfo dci;
    string efilename;
    if (this->efiles_map.find(section_name) == this->efiles_map.end())
        return cycle_map;
    for (auto efile : this->efiles_map[section_name])
    {
        col_names = efile.get_col_set();
        for (auto possible_name : possible_name_cols)
        {
            if (col_names.find(possible_name) != col_names.end())
            {
                name_col = possible_name;
                break;
            }
        }
        if (name_col.size() == 0)
        {
            ss.str("");
            ss << "could not find any possible name cols: ";
            for (auto name : possible_name_cols)
                ss << name << ",";
            ss << " in efile '" << efile.get_filename() << "' columns";

            throw_control_file_error(f_rec, ss.str());
        }
        if (col_names.find("CYCLE") == col_names.end())
            continue;
        str_values = efile.get_col_string_vector("CYCLE");
        str_names = efile.get_col_string_vector(name_col);
        efilename = efile.get_filename();
        for (int i = 0; i < str_values.size(); i++)
        {
            dci = parse_cycle_str(str_values[i],efilename,i,f_rec);
            cycle_map.push_back(make_pair(str_names[i],dci));
        }
    }
    return cycle_map;
}

vector<pair<string, int>> Pest::extract_cycle_numbers2(ofstream& f_rec, string section_name, vector<string> possible_name_cols)
{
	vector<string> str_values, str_names;
	stringstream ss;
	set<string> col_names;
	string name_col = "";
	int cycle;
	vector<pair<string, int>> cycle_map;
	if (this->efiles_map.find(section_name) == this->efiles_map.end())
		return cycle_map;
	for (auto efile : this->efiles_map[section_name])
	{
		col_names = efile.get_col_set();
		for (auto possible_name : possible_name_cols)
		{
			if (col_names.find(possible_name) != col_names.end())
			{
				name_col = possible_name;
				break;
			}

		}

		if (name_col.size() == 0)
		{
			ss.str("");
			ss << "could not find any possible name cols: ";
			for (auto name : possible_name_cols)
				ss << name << ",";
			ss << " in efile '" << efile.get_filename() << "' columns";

			throw_control_file_error(f_rec, ss.str());
		}
		if (col_names.find("CYCLE") == col_names.end())
			continue;
		str_values = efile.get_col_string_vector("CYCLE");
		str_names = efile.get_col_string_vector(name_col);
		for (int i = 0; i < str_values.size(); i++)
		{
			try
			{
				cycle = stoi(str_values[i]);
				cycle_map.push_back(make_pair(str_names[i], cycle));
			}
			catch (...)
			{
				ss.str("");
				ss << "error casting cycle '" << str_values[i] << "' to int on row " << i << "of external file " << efile.get_filename() << " , Stopped...";
				throw_control_file_error(f_rec, ss.str());
			}

		}
	}
	return cycle_map;

}

void Pest::child_pest_update(int icycle, bool keep_order)
{
	/*
	Update Pest members to reflect current cycle data only
	*/
	vector<string> cycle_grps, unique_cycle_grps, grps;	
	vector<string> parnames, obsnames;
	
	//prior_info_string
	//control_info
	//pareto_info
	//svd_info

	if (keep_order) {
        parnames = get_ctl_ordered_par_names();
        obsnames = get_ctl_ordered_obs_names();
    }
	else
    {
        //parnames = get_ctl_parameters().get_keys();
        //obsnames = get_ctl_observations().get_keys();
    }

	ParameterInfo pi = get_ctl_parameter_info();
	for (auto& p : ctl_parameters.get_keys()) {
        if (!cycle_in_range(icycle, pi.get_parameter_rec_ptr(p)->dci)) {
            ctl_parameter_info.erase(p);
            ctl_parameters.erase(p);
            base_group_info.par_erase(p);// .parameter2group.erase(p);
            if (keep_order) {
                parnames.erase(remove(parnames.begin(),
                                      parnames.end(), p), parnames.end());
            }
        }
//        if (pi.get_parameter_rec_ptr(p)->cycle != icycle &&
//            pi.get_parameter_rec_ptr(p)->cycle >= 0) {
//            ctl_parameter_info.erase(p);
//            ctl_parameters.erase(p);
//            base_group_info.par_erase(p);// .parameter2group.erase(p);
//            parnames.erase(remove(parnames.begin(),
//                                  parnames.end(), p), parnames.end());
//
//        }
        else {
            cycle_grps.push_back(base_group_info.get_group_name(p));

        }
    }
    // get unique groups
	for (auto curr = cycle_grps.begin(); curr != cycle_grps.end(); curr++) {
		if (find(unique_cycle_grps.begin(), unique_cycle_grps.end(), *curr) == unique_cycle_grps.end())
		{
			unique_cycle_grps.push_back(*curr);
		}
	}
	// remove groups not in current cycle
	grps = base_group_info.get_group_names();
	for (auto grp = grps.begin(); grp != grps.end(); grp++)
	{
		if (find(unique_cycle_grps.begin(), unique_cycle_grps.end(), *grp) == unique_cycle_grps.end())
			base_group_info.grp_erase(*grp);
	}
	// update observations
	vector<string> cycle_grps_o, unique_cycle_grps_o, grps_o;
	ObservationInfo obs_info = get_ctl_observation_info();
	for (auto& ob : observation_values.get_keys())
	    if (!cycle_in_range(icycle,obs_info.get_observation_rec_ptr(ob)->dci))
//		if (obs_info.get_observation_rec_ptr(ob)->cycle != icycle &&
//			obs_info.get_observation_rec_ptr(ob)->cycle >= 0)
		{
			observation_info.erase_ob(ob);
			observation_values.erase(ob);
			if (keep_order) {
                obsnames.erase(remove(obsnames.begin(),
                                      obsnames.end(), ob), obsnames.end());
            }
		}
		else
		{
			cycle_grps_o.push_back(obs_info.get_group(ob));
		}

	// get unique groups
	for (auto curr = cycle_grps_o.begin(); curr != cycle_grps_o.end(); curr++) {
		if (find(unique_cycle_grps_o.begin(), unique_cycle_grps_o.end(), *curr) == unique_cycle_grps_o.end())
		{
			unique_cycle_grps_o.push_back(*curr);
		}
	}
	// remove groups not in current cycle
	grps_o = observation_info.get_groups();
	for (auto grp = grps_o.begin(); grp != grps_o.end(); grp++)
	{
		if (find(unique_cycle_grps_o.begin(), unique_cycle_grps_o.end(), *grp) == unique_cycle_grps_o.end())
			observation_info.erase_gp(*grp);
	}

	// prior info

	//ctl_parameters
	//model_exec_info
	ModelExecInfo curr_mod_exe = model_exec_info;
	std::vector<std::string> _tpl_vec, _infile_vec, _ins_vec, _out_vec;
	std::vector<DaCycleInfo> incy, outcy;
	int index = 0;
//	for (auto ic : model_exec_info.incycle_vec)
//	{
//		if ((ic == icycle) || (ic < 0))
//		{
//			_tpl_vec.push_back(model_exec_info.tplfile_vec[index]);
//			_infile_vec.push_back(model_exec_info.inpfile_vec[index]);
//			incy.push_back(model_exec_info.incycle_vec[index]);
//		}
//
//		index = index + 1;
//	}

    for (auto ic : model_exec_info.incycle_dci_vec)
    {
        if (cycle_in_range(icycle,ic))
        {
            _tpl_vec.push_back(model_exec_info.tplfile_vec[index]);
            _infile_vec.push_back(model_exec_info.inpfile_vec[index]);
            incy.push_back(model_exec_info.incycle_dci_vec[index]);
        }

        index = index + 1;
    }

	//output files
	index = 0;
//	for (auto ic : model_exec_info.outcycle_vec)
//	{
//		if ((ic == icycle) || (ic < 0))
//		{
//			_ins_vec.push_back(model_exec_info.insfile_vec[index]);
//			_out_vec.push_back(model_exec_info.outfile_vec[index]);
//			outcy.push_back(model_exec_info.outcycle_vec[index]);
//		}
//
//		index = index + 1;
//	}
    for (auto ic : model_exec_info.outcycle_dci_vec)
    {
        if (cycle_in_range(icycle,ic))
        {
            _ins_vec.push_back(model_exec_info.insfile_vec[index]);
            _out_vec.push_back(model_exec_info.outfile_vec[index]);
            outcy.push_back(model_exec_info.outcycle_dci_vec[index]);
        }

        index = index + 1;
    }
	model_exec_info.tplfile_vec = _tpl_vec;
	model_exec_info.inpfile_vec = _infile_vec;
	model_exec_info.incycle_dci_vec = incy;
	model_exec_info.insfile_vec = _ins_vec;
	model_exec_info.outfile_vec = _out_vec;
	model_exec_info.outcycle_dci_vec = outcy;


	//pestpp_options
	//base_par_transform
	//ctl_ordered_par_names
	if (keep_order) {
        ctl_ordered_par_names = parnames;

        //ctl_ordered_obs_names
        ctl_ordered_obs_names = obsnames;

        //ctl_ordered_par_group_names ----> TODO: Check if the groups order is preserved..
        ctl_ordered_par_group_names = unique_cycle_grps;
        ctl_ordered_obs_group_names = unique_cycle_grps_o;
    }
	
	// get number of adj par for current cycle
	ParameterRec::TRAN_TYPE tfixed = ParameterRec::TRAN_TYPE::FIXED;
	ParameterRec::TRAN_TYPE ttied = ParameterRec::TRAN_TYPE::TIED;
	int new_n_adj_par = 0;
	for (auto pname : ctl_ordered_par_names)
	{
		if ((ctl_parameter_info.get_parameter_rec_ptr(pname)->tranform_type != tfixed) &&
			(ctl_parameter_info.get_parameter_rec_ptr(pname)->tranform_type != ttied))
		{
			new_n_adj_par++;
		}
	}
	n_adj_par = new_n_adj_par;
	//this->get_pestpp_options_ptr()->set_check_tplins(false);
	//get_pestpp_options_ptr->set_check_tplins(false);
	//this.check_inputs();
}

vector<int> Pest::get_assim_dci_cycles(ofstream& f_rec, vector<int> unique_cycles)
{
    stringstream ss;
    set<int> scycles(unique_cycles.begin(),unique_cycles.end());

    int start_cycle = std::numeric_limits<int>::max();
    int stop_cycle = std::numeric_limits<int>::min();
    for (auto &s : scycles)
    {
        if (s < start_cycle)
            start_cycle = s;
        if (s > stop_cycle)
            stop_cycle = s;
    }
    DaCycleInfo dci;
    for (auto pname : ctl_ordered_par_names)
    {
        dci = ctl_parameter_info.get_parameter_rec_ptr(pname)->dci;
        if (dci.start < start_cycle)
            start_cycle = dci.start;
        if (dci.stop > stop_cycle)
            stop_cycle = dci.stop;
    }
    for (auto oname : ctl_ordered_obs_names)
    {
        dci = observation_info.get_observation_rec_ptr(oname)->dci;
        if (dci.start < start_cycle)
            start_cycle = dci.start;
        if (dci.stop > stop_cycle)
            stop_cycle = dci.stop;
    }

    if (start_cycle < 0)
    {
        ss.str("");
        ss << "minimum start '" << start_cycle << "' cycle less than zero";
        throw_control_file_error(f_rec,ss.str());
    }

    if (stop_cycle < 0) {
        ss.str("");
        ss << "WARNING: didn't find any explicit 'stop' cycle values in control file info, assuming smoother formulation" << endl;
        ss << "         assigning a generic 'stop' value of " << start_cycle + 1 << " which is 'start' cycle plus 1" << endl;
        cout << ss.str();
        f_rec << ss.str();
        stop_cycle = start_cycle + 1;
        return vector<int>{start_cycle};

    }
    if (start_cycle > stop_cycle) {
        ss.str("");
        ss << "minimum start cycle '" << start_cycle << "' greater than maximum stop cycle '" << stop_cycle << "'";
        throw_control_file_error(f_rec, ss.str());
    }



    for (int i=start_cycle;i<stop_cycle+1;i++)
    {
        for (auto pname : ctl_ordered_par_names)
        {
            if (cycle_in_range(i,ctl_parameter_info.get_parameter_rec_ptr(pname)->dci))
                scycles.emplace(i);
        }
        for (auto oname : ctl_ordered_obs_names)
        {
            if (cycle_in_range(i,observation_info.get_observation_rec_ptr(oname)->dci))
                scycles.emplace(i);
        }
    }
    return vector<int>(scycles.begin(),scycles.end());
}

bool cycle_in_range(int cycle,const DaCycleInfo& dci) {
    if ((dci.start == dci.stop) && (dci.start == cycle))
        return true;
    if (dci.start > cycle)
        return false;
    if ((dci.stop >= 0) && (dci.stop < cycle))
        return false;
    if (dci.stride == 1)
        return true;
    if (abs(cycle - dci.start) % dci.stride == 0)
        return true;
    return false;

}

bool every_cycle(const DaCycleInfo& dci)
{
    if (dci.stop != -999)
        return false;
    if (dci.start != 0)
        return false;
    if (dci.stride != 1)
        return false;
    return true;
}



//vector<int> Pest::get_assim_cycles(ofstream& f_rec, vector<int> unique_cycles)
//{
//
//	int curr_cycle;
//	vector<int> cycles_ordered_list;
//
//	for (auto pname : ctl_ordered_par_names)
//	{
//		curr_cycle = ctl_parameter_info.get_parameter_rec_ptr(pname)->cycle;
//		if (curr_cycle>= 0)
//			cycles_ordered_list.push_back(curr_cycle);
//	}
//	for (auto oname : ctl_ordered_obs_names)
//	{
//		curr_cycle = observation_info.get_observation_rec_ptr(oname)->cycle;
//		if (curr_cycle >= 0)
//			cycles_ordered_list.push_back(curr_cycle);
//	}
//
//
//	// get unique groups
//	for (auto curr = cycles_ordered_list.begin(); curr != cycles_ordered_list.end(); curr++) {
//		if (find(unique_cycles.begin(), unique_cycles.end(), *curr) == unique_cycles.end())
//		{
//			unique_cycles.push_back(*curr);
//		}
//	}
//	if (unique_cycles.size() == 0)
//		throw_control_file_error(f_rec, "no valid cycle info found in par and obs data (possibly all '-1' cycles?)");
//
//	sort(unique_cycles.begin(), unique_cycles.end());
//	stringstream ss;
//	//to catch common non-zero-indexing
//	if (unique_cycles.front() != 0)  //recall a cycle of -1 is just a flag
//	{
//		ss.str("");
//		ss << "a cycle with index of zero does not exist; cycling needs to start at zero for initialization...";
//		throw_control_file_error(f_rec, ss.str());
//	}
//	return unique_cycles;
//
//
//}



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
	if (regul_scheme_ptr != 0)
	{
		try
		{
			delete regul_scheme_ptr;
		}
		catch (...)
		{
		}
	}
    base_group_info.free_mem();
}

pair<string, string> Pest::parse_keyword_line(ofstream &f_rec, const string &line)
{
	string key;
	string value;
	
	// look for a comment char
	size_t found = line.find_first_of("#");
	if (found == string::npos) {
		found = line.length();
	}
	string tmp_line = line.substr(0, found);
	//strip any leading or trailing whitespace
	strip_ip(tmp_line, "both", "\t\n\r+ ");
	//remove any quote chars
	tmp_line.erase(remove(tmp_line.begin(), tmp_line.end(), '\"'), tmp_line.end());
	tmp_line.erase(remove(tmp_line.begin(), tmp_line.end(), '\''), tmp_line.end());
	
	//split on whitespaces
	vector<string> tokens;
	tokenize(tmp_line,tokens,"\t ");
	if (tokens.size() < 2)
	{
		throw_control_file_error(f_rec, "Pest::parse_keyword_line() error: too few tokens on line '" + line + "', need at least 2");
	}
	key = tokens[0];
	upper_ip(key);
	//the other chars to the right are the "value" (could have spaces, commas, etc)
	value = tmp_line.substr(key.size(), tmp_line.size());
	strip_ip(value, "both", "\t\n\r+ ");

	return pair<string, string>(key,value);
}

void Pest::throw_control_file_error(ofstream& f_rec,const string &message, bool should_throw, bool echo)
{
	stringstream ss;
	if (should_throw)
		ss << "control file parsing error: " << message << endl;
	else
		ss << "control file parsing warning: " << message << endl;

	if (echo)
    {
	    cout << ss.str();
    }
	f_rec << ss.str();
	
	if (should_throw)
	{
	    if (echo){
            cerr << ss.str();
	    }

		//f_rec.close();
		throw runtime_error(ss.str());
	}	
}

void Pest::check_report_assignment(ofstream &f_rec, PestppOptions::ARG_STATUS stat, const string &key, const string &org_value)
{
	if (stat == PestppOptions::ARG_STATUS::ARG_INVALID)
	{
		throw_control_file_error(f_rec, "invalid value for key,value pair '" + key + ", " + org_value);
	}
	if (stat == PestppOptions::ARG_STATUS::ARG_DUPLICATE)
	{
		throw_control_file_error(f_rec, "duplicate entry for key,value pair '" + key + ", " + org_value);
	}
}

void Pest::tokens_to_par_group_rec(ofstream &f_rec, const vector<string>& tokens)
{
	ParameterGroupRec pgi;
	string name = tokens[0];
	size_t n_tokens = tokens.size();
	
	if (s_pargp.find(name) == s_pargp.end())
	{
		ctl_ordered_par_group_names.push_back(name);
		s_pargp.emplace(name);
	}

	pgi.name = name;
	convert_ip(tokens[1], pgi.inctyp);
	convert_ip(tokens[2], pgi.derinc);
	convert_ip(tokens[3], pgi.derinclb);
	convert_ip(tokens[4], pgi.forcen);
	convert_ip(tokens[5], pgi.derincmul);
	convert_ip(tokens[6], pgi.dermthd);
	if (n_tokens >= 8) convert_ip(tokens[7], pgi.splitthresh);
	if (n_tokens >= 9) convert_ip(tokens[8], pgi.splitreldiff);

    if ((pgi.derinc != 0) && (!isnormal(pgi.derinc)))
        throw_control_file_error(f_rec,"denormal derinc '"+tokens[4]+"' for parameter group "+tokens[0]);
    base_group_info.insert_group(name, pgi);

}

void Pest::tokens_to_par_rec(ofstream &f_rec, const vector<string>& tokens, TranFixed* t_fixed, TranLog10* t_log, TranScale* t_scale, TranOffset* t_offset)
{
	ParameterRec pi;
	double scale;
	double offset;
	string name = tokens[0];
	string trans_type = tokens[1];
	convert_ip(tokens[2], pi.chglim);
	convert_ip(tokens[3], pi.init_value);
	convert_ip(tokens[4], pi.lbnd);
	convert_ip(tokens[5], pi.ubnd);
	convert_ip(tokens[6], pi.group);
	convert_ip(tokens[7], scale);
	convert_ip(tokens[8], offset);

    if ((pi.init_value != 0) && (!isnormal(pi.init_value)))
        throw_control_file_error(f_rec,"denormal parval1 '"+tokens[3]+"' for parameter "+tokens[0]);
    if ((pi.lbnd != 0) && (!isnormal(pi.lbnd)))
        throw_control_file_error(f_rec,"denormal parlbnd '"+tokens[4]+"' for parameter "+tokens[0]);
    if ((pi.ubnd != 0) && (!isnormal(pi.ubnd)))
        throw_control_file_error(f_rec,"denormal parubnd '"+tokens[5]+"' for parameter "+tokens[0]);
    if ((pi.scale != 0) && (!isnormal(pi.scale)))
        throw_control_file_error(f_rec,"denormal scale '"+tokens[7]+"' for parameter "+tokens[0]);
    if ((pi.offset != 0) && (!isnormal(pi.offset)))
        throw_control_file_error(f_rec,"denormal offset '"+tokens[8]+"' for parameter "+tokens[0]);

    if (control_info.numcom > 1)
        try {

            convert_ip(tokens[9], pi.dercom);
        }
	    catch (...)
        {
	        float f;
	        convert_ip(tokens[9],f);
	        pi.dercom = (int)std::floor(f + 0.5);
        }
	else
		pi.dercom = 1;
	
	pi.scale = scale;
	pi.offset = offset;
	// add parameters to model parameter and paramter_info datasets
	ctl_ordered_par_names.push_back(name);
	if (trans_type == "FIXED")
	{
		pi.tranform_type = ParameterRec::TRAN_TYPE::FIXED;
	}
	else if (trans_type == "LOG")
	{
		pi.tranform_type = ParameterRec::TRAN_TYPE::LOG;
		n_adj_par++;
	}
	else if (trans_type == "TIED")
	{
		pi.tranform_type = ParameterRec::TRAN_TYPE::TIED;
	}
	else if (trans_type == "NONE")
	{
		pi.tranform_type = ParameterRec::TRAN_TYPE::NONE;
		n_adj_par++;
	}
	else
	{
		throw_control_file_error(f_rec, "unrecognized partrans for par " + name + ": " + trans_type);
	}
	ctl_parameter_info.insert(name, pi);
	if (ctl_parameters.find(name) != ctl_parameters.end())
        throw_control_file_error(f_rec,"duplicate parameter names in control file for: '"+name+"'");
	ctl_parameters.insert(name, pi.init_value);
	
	//if (find(ctl_ordered_par_group_names.begin(), ctl_ordered_par_group_names.end(), pi.group) == ctl_ordered_par_group_names.end())
	if (s_pargp.find(pi.group) == s_pargp.end())
	{
		ParameterGroupRec pgr;
		pgr.set_defaults();
		pgr.set_name(pi.group);
		base_group_info.insert_group(pi.group,pgr);
		ctl_ordered_par_group_names.push_back(pi.group);
		s_pargp.emplace(pi.group);
	}
	base_group_info.insert_parameter_link(name, pi.group);

	// build appropriate transformations
	if (trans_type == "FIXED") {
		t_fixed->insert(name, pi.init_value);
	}
	else if (trans_type == "LOG") {
		t_log->insert(name);
	}
	if (offset != 0) {
		t_offset->insert(name, offset);
	}
	if (scale != 1) {
		t_scale->insert(name, scale);
	}

}

void Pest::tokens_to_obs_group_rec(ofstream& f_rec, const vector<string>& tokens)
{
	string name = tokens[0];
	if (tokens.size() > 1)
	{
		stringstream ss;
		ss << "observation covariance matrix detected for group '" << tokens[0] << "' - these are not supported...yet!";
		throw_control_file_error(f_rec, ss.str());
	}
	ObservationGroupRec group_rec;
	observation_info.groups[name] = group_rec;
	//vector<string>::iterator is = find(ctl_ordered_obs_group_names.begin(), ctl_ordered_obs_group_names.end(), name);
	//if (is == ctl_ordered_obs_group_names.end())
	if (s_obgnme.find(name) == s_obgnme.end())
	{
		ctl_ordered_obs_group_names.push_back(name);
		s_obgnme.emplace(name);
	}
}

void Pest::tokens_to_obs_rec(ofstream& f_rec, const vector<string> &tokens)
{
	ObservationRec obs_i;
	string name = tokens[0];
	double value;
	size_t idx;
	//convert_ip(tokens[1], value);
    try
    {
        value = stod(tokens[1],&idx);
    }
	catch (...)
    {
	    throw_control_file_error(f_rec,"error parsing obsval '"+tokens[1]+"' for observation "+tokens[0]);
    }
	if (idx != tokens[1].size())
    {
        throw_control_file_error(f_rec,"error parsing obsval '"+tokens[1]+"' for observation "+tokens[0]);
    }
	//convert_ip(tokens[2], obs_i.weight);
	try {
        obs_i.weight = stod(tokens[2],&idx);
    }
	catch (...)
    {
        throw_control_file_error(f_rec,"error parsing weight '"+tokens[2]+"' for observation "+tokens[0]);
    }
    if (idx != tokens[2].size())
    {
        throw_control_file_error(f_rec,"error parsing weight '"+tokens[2]+"' for observation "+tokens[0]);
    }
	if ((value != 0) && (!isnormal(value)))
        throw_control_file_error(f_rec,"denormal obsval '"+tokens[1]+"' for observation "+tokens[0]);
    if (((obs_i.weight!= 0) && !isnormal(obs_i.weight)))
        throw_control_file_error(f_rec,"denormal weight '"+tokens[2]+"' for observation "+tokens[0]);

    obs_i.group = tokens[3];
	if (observation_values.find(name) != observation_values.end())
        throw_control_file_error(f_rec,"duplicate observation names in control file for: '" + name + "'");
	ctl_ordered_obs_names.push_back(name);
	observation_info.observations[name] = obs_i;
	observation_values.insert(name, value);
	name = obs_i.group;
	//vector<string>::iterator is = find(ctl_ordered_obs_group_names.begin(), ctl_ordered_obs_group_names.end(), name);
	//if (is == ctl_ordered_obs_group_names.end())
	if (s_obgnme.find(name) == s_obgnme.end())
	{
		ctl_ordered_obs_group_names.push_back(name);
		s_obgnme.emplace(name);
        ObservationGroupRec ogr;
        observation_info.groups[obs_i.group] = ogr;
	}
}

void Pest::tokens_to_pi_rec(ofstream& f_rec, const string& line_upper)
{
	string first = line_upper.substr(0, 1);
	if (!prior_info_string.empty() && first != "&") {
		pair<string, string> pi_name_group = prior_info.AddRecord(prior_info_string);
		ctl_ordered_pi_names.push_back(pi_name_group.first);
		vector<string>::iterator is = find(ctl_ordered_obs_group_names.begin(), ctl_ordered_obs_group_names.end(), pi_name_group.second);
		if (is == ctl_ordered_obs_group_names.end())
		{
			ctl_ordered_obs_group_names.push_back(pi_name_group.second);
		}
		prior_info_string.clear();
	}
	else if (first == "&") {
		prior_info_string.append(" ");
		
	}
	prior_info_string.append(line_upper);

}


void Pest::tokens_to_pi_rec(ofstream& f_rec, const vector<string>& tokens)
{
	
	if (prior_info.find(tokens[0]) != prior_info.end())
        throw_control_file_error(f_rec,"duplicate prior info names in control file for: '"+tokens[0]+"'");
	pair<string,string> pi_name_group = prior_info.AddRecord(tokens);
	ctl_ordered_pi_names.push_back(pi_name_group.first);
	//vector<string>::iterator is = find(ctl_ordered_obs_group_names.begin(), ctl_ordered_obs_group_names.end(), pi_name_group.second);
	//if (is == ctl_ordered_obs_group_names.end())
	if (s_obgnme.find(pi_name_group.second) == s_obgnme.end())
	{
		ctl_ordered_obs_group_names.push_back(pi_name_group.second);
        s_obgnme.emplace(pi_name_group.second);
	}
	
}

void Pest::rectify_par_groups()
{
	const ParameterRec *pr;
	vector<string> temp = base_group_info.get_group_names();
	set<string> group_names(temp.begin(), temp.end());
	for (auto p : ctl_ordered_par_names)
	{
		pr = ctl_parameter_info.get_parameter_rec_ptr(p);
		if (group_names.find(pr->group) == group_names.end())
		{
			ParameterGroupRec pgi;
			pgi.set_defaults();
			pgi.name = pr->group;
			
			//update containers
			ctl_ordered_par_group_names.push_back(pr->group);
			base_group_info.insert_group(pr->group, pgi);
			temp = base_group_info.get_group_names();
			group_names.clear();
			group_names.insert(temp.begin(), temp.end());
			ctl_ordered_par_group_names.push_back(pr->group);

		}
	}

}

map<string, double> Pest::calc_par_dss(const Jacobian& jac, ParamTransformSeq& par_transform)
{
	Parameters pars = jac.get_base_numeric_parameters();
	Parameters ctl_pars = par_transform.numeric2ctl_cp(pars);
	const vector<string>& par_list = jac.parameter_list();
	const vector<string>& obs_list = jac.obs_and_reg_list();
	Eigen::VectorXd par_vec = pars.get_data_eigen_vec(par_list);
	Eigen::MatrixXd par_mat_tmp = par_vec.asDiagonal();
	Eigen::SparseMatrix<double> par_mat = par_mat_tmp.sparseView();
	QSqrtMatrix Q_sqrt(&observation_info, &prior_info);
	Eigen::SparseMatrix<double> q_sqrt_no_reg = Q_sqrt.get_sparse_matrix(obs_list, DynamicRegularization::get_zero_reg_instance());
	Eigen::SparseMatrix<double> dss_mat_no_reg_pest = q_sqrt_no_reg * jac.get_matrix(obs_list, par_list);
	int n_nonzero_weights_no_reg = observation_info.get_nnz_obs();
	map<string, double> par_sens;
	double val;
	for (int i = 0; i < par_list.size(); ++i)
	{
		val = dss_mat_no_reg_pest.col(i).norm() / n_nonzero_weights_no_reg;
		par_sens[par_list[i]] = val;
	}
	return par_sens;
	return par_sens;
}

map<string,string> Pest::get_ext_file_string_map(const string& section_name, const string& col_name)
{
	string sname_upper = pest_utils::upper_cp(section_name);
	map<string, string> val_map;
	if (efiles_map.find(sname_upper) == efiles_map.end())
		return val_map;
	set<string> col_names;
	string cname_upper = pest_utils::upper_cp(col_name);
	string idx_col;
	vector<string> svals, sidx;
	double val;
	for (auto efile : efiles_map[sname_upper])
	{
		idx_col = efile.get_index_col_name();
		if (idx_col.size() == 0)
			continue;
		col_names = efile.get_col_set();
		if (col_names.find(cname_upper) != col_names.end())
		{
			svals = efile.get_col_string_vector(cname_upper);
			sidx = efile.get_col_string_vector(idx_col);
			for (int i = 0; i < svals.size(); i++)
			{
				if (svals[i].size() > 0)
					val_map.emplace(make_pair(sidx[i],svals[i]));
			}
		}
	}
	return val_map;
}

void Pest::release_unused_for_agent()
{
    ctl_ordered_obs_names.clear();
    ctl_ordered_obs_group_names.clear();
    ctl_ordered_par_names.clear();
    ctl_ordered_pi_names.clear();
    ctl_ordered_par_group_names.clear();
    //base_group_info.clear();
    prior_info.clear();
    base_group_info.free_mem();
    regul_scheme_ptr = NULL;

}


map<string, double> Pest::get_ext_file_double_map(const string& section_name, const string& col_name)
{
	string sname_upper = pest_utils::upper_cp(section_name);
	map<string, string> str_map = get_ext_file_string_map(section_name, col_name);
	map<string, double> val_map;
	double val;
	for (auto sval : str_map)
	{
		try
		{
			val = stod(sval.second);
			val_map[sval.first] = val;
		}
		catch (...)
		{
			continue;
		}
	}
	return val_map;
}


ostream& operator<< (ostream &os, const Pest& val)
{
	os << val.control_info;
	os << val.svd_info;
	return os;
}

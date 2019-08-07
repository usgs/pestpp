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
#include <iomanip>
#include <vector>
#include <map>
#include <sstream>
#include <list>
#include <regex>
#include "pest_data_structs.h"
#include "utilities.h"
#include "pest_error.h"
#include "ParamTransformSeq.h"
#include "Transformation.h"
#include "pest_error.h"
#include "Transformable.h"

using namespace::std;
using namespace::pest_utils;


ostream& operator<< (ostream &os, const ControlInfo& val)
{
	os << "PEST Control Information" << endl;
	os << "    relparmax = " << val.relparmax << endl;
	os << "    facparmax = " << val.facparmax << endl;
	os << "    facorig = " << val.facorig << endl;
	os << "    phiredswh = " << val.phiredswh << endl;
	os << "    noptmax = " << val.noptmax << endl;
	os << "    phiredstp = " << val.phiredstp << endl;
	os << "    nphistp = " << val.nphistp << endl;
	os << "    nphinored = " << val.nphinored << endl;
	os << "    relparstp = " << val.relparstp << endl;
	os << "    nrelpar = " << val.nrelpar << endl;
	return os;
}

//////////////////  ParameterGroupRec Methods ////////////////////////////////

ParameterGroupRec& ParameterGroupRec::operator=(const ParameterGroupRec &rhs)
{
	name = rhs.name;
	inctyp = rhs.inctyp;
	derinc = rhs.derinc;
	derinclb = rhs.derinclb;
	forcen = rhs.forcen;
	derincmul = rhs.derincmul;
	dermthd = rhs.dermthd;
	splitthresh = rhs.splitthresh;
	splitreldiff = rhs.splitreldiff;
	return *this;
}

ostream& operator<< (ostream &os, const ParameterGroupRec& val)
{
	os << "PEST Parameter Group Information" << endl;
	os << "    name   = " << val.name << endl;
	os << "    inctyp = " << val.inctyp << endl;
	os << "    derinc = " << val.derinc << endl;
	os << "    derinclb = " << val.derinclb << endl;
	os << "    forcen = " << val.forcen << endl;
	os << "    derincmul = " << val.derincmul << endl;
	return os;
}

/////////////////////////////////////  ParameterGroupInfo Methods ///////////////////////
const ParameterGroupRec* ParameterGroupInfo::get_group_rec_ptr(const string &name) const
{
	const ParameterGroupRec *ret_val = 0;
	unordered_map<string, ParameterGroupRec*>::const_iterator g_iter;

	g_iter = parameter2group.find(name);
	if(g_iter != parameter2group.end()) {
		ret_val = (*g_iter).second;
	}
	return ret_val;
}

ParameterGroupRec* ParameterGroupInfo::get_group_rec_ptr_4_mod(const string &name)
{
	ParameterGroupRec *ret_val = 0;
	unordered_map<string, ParameterGroupRec*>::const_iterator g_iter;

	g_iter = parameter2group.find(name);
	if (g_iter != parameter2group.end()) {
		ret_val = (*g_iter).second;
	}
	return ret_val;
}

string ParameterGroupInfo::get_group_name(const string &par_name) const
{
	return get_group_rec_ptr(par_name)->name;
}

void ParameterGroupInfo::insert_group(const string &group_name, ParameterGroupRec &rec)
{
	groups[group_name] = new ParameterGroupRec(rec);
}

void ParameterGroupInfo::insert_parameter_link(const string &parameter_name, const string & group_name)
{
	unordered_map<string, ParameterGroupRec*>::const_iterator g_iter;

	g_iter = groups.find(group_name);
	if(g_iter == groups.end()) {
		throw PestIndexError(group_name, "Invalid parameter group name");
	}
	else {
		parameter2group[parameter_name] = (*g_iter).second;
	}

}

const ParameterGroupInfo& ParameterGroupInfo::operator=(const ParameterGroupInfo &rhs)
{
	unordered_map<ParameterGroupRec*, ParameterGroupRec*> old2new;
	unordered_map<string, ParameterGroupRec*>::const_iterator it(rhs.groups.begin());
	unordered_map<string, ParameterGroupRec*>::const_iterator end(rhs.groups.end());
	for (; it != end; ++it) {
		ParameterGroupRec* new_ptr = new ParameterGroupRec(*(*it).second);
		groups[(*it).first] = new_ptr;
		old2new[(*it).second] = new_ptr;
	}
	unordered_map<ParameterGroupRec*, ParameterGroupRec*>::iterator it_find;
	it = rhs.parameter2group.begin();
	end =  rhs.parameter2group.end();
	for (; it != end; ++it) {
		it_find = old2new.find((*it).second);
		parameter2group[(*it).first] = (*it_find).second;
	}
	return *this;
}

vector<string> ParameterGroupInfo::get_group_names() const
{
	vector<string> group_names;
	for (auto &g : groups)
		group_names.push_back(g.first);
	return group_names;
}

bool ParameterGroupInfo::have_switch_derivative() const
{
	bool switch_der = false;
	for (const auto irec : groups)
	{
		if (lower_cp(irec.second->forcen) == "switch")
		{
			switch_der = true;
			break;
		}
	}
	return switch_der;
}

ParameterGroupInfo::~ParameterGroupInfo()
{
	unordered_map<string, ParameterGroupRec*>::iterator it(groups.begin());
	unordered_map<string, ParameterGroupRec*>::iterator end(groups.end());
	for (; it != end; ++it) {
		delete (*it).second;
	}
}

ostream& operator<< (ostream &os, const ParameterGroupInfo &val)
{
	unordered_map<string, ParameterGroupRec*>::const_iterator it(val.parameter2group.begin());
	unordered_map<string, ParameterGroupRec*>::const_iterator end(val.parameter2group.end());
	for (; it != end; ++it) {
		os << *(it->second) << endl;
	}
	return os;
}


/////////////////////////////////////  ParameterRec Methods ///////////////////////

ostream& operator<< (ostream &os, const ParameterRec& val)
{
	os << "PEST Parameter Information" << endl;
	os << "    chglim = " << val.chglim << endl;
	os << "    lbnd = " << val.lbnd << endl;
	os << "    ubnd = " << val.ubnd << endl;
	os << "    group = " << val.group << endl;
	os << "    dercom = " << val.dercom << endl;
	return os;
}

/////////////////////////////////////  CtlParameterInfo Methods ///////////////////////



const ParameterRec* ParameterInfo::get_parameter_rec_ptr(const string &name) const
{
	const ParameterRec *ret_val = 0;
	unordered_map<string, ParameterRec>::const_iterator p_iter;

	p_iter = parameter_info.find(name);
	if(p_iter != parameter_info.end()) {
		ret_val = &((*p_iter).second);
	}
	return ret_val;
}

ParameterRec* ParameterInfo::get_parameter_rec_ptr_4_mod(const string &name)
{
	ParameterRec *ret_val = 0;
	unordered_map<string, ParameterRec>::iterator p_iter;

	p_iter = parameter_info.find(name);
	if (p_iter != parameter_info.end()) {
		ret_val = &((*p_iter).second);
	}
	return ret_val;
}

Parameters ParameterInfo::get_low_bnd(const vector<string> &keys) const
{
	Parameters l_bnd;
	vector<string>::const_iterator iend=keys.end();
	ParameterRec const *v_ptr;
	for (vector<string>::const_iterator i=keys.begin(); i!=iend; ++i)
	{
		v_ptr = get_parameter_rec_ptr(*i);
		if (v_ptr) {
			l_bnd.insert(*i, v_ptr->lbnd);
		}
		else {
			l_bnd.insert(*i, Parameters::no_data);
		}
	}
	return l_bnd;
}

Parameters ParameterInfo::get_up_bnd(const vector<string> &keys) const
{
	Parameters u_bnd;
	vector<string>::const_iterator iend=keys.end();
	ParameterRec const *v_ptr;
	for (vector<string>::const_iterator i=keys.begin(); i!=iend; ++i)
	{
		v_ptr = get_parameter_rec_ptr(*i);
		if (v_ptr) {
			u_bnd.insert(*i, v_ptr->ubnd);
		}
		else {
			u_bnd.insert(*i, Parameters::no_data);
		}
	}
	return u_bnd;
}

Parameters ParameterInfo::get_init_value(const vector<string> &keys) const
{
	Parameters init_value;
	vector<string>::const_iterator iend=keys.end();
	ParameterRec const *v_ptr;
	for (vector<string>::const_iterator i=keys.begin(); i!=iend; ++i)
	{
		v_ptr = get_parameter_rec_ptr(*i);
		if (v_ptr) {
			init_value.insert(*i, v_ptr->init_value);
		}
		else {
			init_value.insert(*i, Parameters::no_data);
		}
	}
	return init_value;
}


ostream& operator<< (ostream &os, const PestppOptions& val)
{
	os << "PEST++ Options" << endl;
	os << "    n_iter_base = " << left << setw(20) << val.get_n_iter_base() << endl;
	os << "    n_iter_super = " << left << setw(20) << val.get_n_iter_super() << endl;
	os << "    max_n_super = " << left << setw(20) << val.get_max_n_super() << endl;
	os << "    super eigthres = " << left << setw(20) << val.get_super_eigthres() << endl;
	os << "    svd pack = " << left << setw(20) << val.get_svd_pack() << endl;
	os << "    super relparmax = " << left << setw(20) << val.get_super_relparmax() << endl;
	os << "    max super frz iter = " << left << setw(20) << val.get_max_super_frz_iter() << endl;
	os << "    max run fail = " << left << setw(20) << val.get_max_run_fail() << endl;
	os << "    max reg iter = " << left << setw(20) << val.get_max_reg_iter() << endl;
	os << "    use jacobian scaling a la PEST? = ";
	if (val.get_jac_scale())
		os << " yes" << endl;
	else
		os << " no" << endl;
	os << "    lambdas = " << endl;
	for (auto &lam : val.get_base_lambda_vec())
	{
		os << right << setw(15) << lam << endl;
	}
	os << "    lambda scaling factors = " << endl;
	for (auto &ls : val.get_lambda_scale_vec())
		os << right << setw(15) << ls << endl;

	if (!val.get_basejac_filename().empty())
	{
		os << "   restarting with existing jacobian matrix file: " << val.get_basejac_filename() << endl;
		if (!val.get_hotstart_resfile().empty())
		{
			os << "   and using existing residual file " << val.get_hotstart_resfile() << " to forego initial model run" << endl;
		}
	}
	if (val.get_uncert_flag())
	{
		os << "    using FOSM-based uncertainty estimation for parameters" << endl;
		os << "    parameter covariance file = " << left << setw(20) << val.get_parcov_filename() << endl;
		if (val.get_prediction_names().size() > 0)
		{
			os << "    using FOSM-based uncertainty for forecasts" << endl;
			os << "    forecast names = " << endl;
		}

	}
	for (auto &pname : val.get_prediction_names())
		os << right << setw(15) << pname << endl;
	os << "    derivative run failure forgive = " << left << setw(15) << val.get_der_forgive() << endl;
	os << "    run overdue reschedule factor = " << left << setw(20) << val.get_overdue_reched_fac() << endl;
	os << "    run overdue giveup factor = " << left << setw(20) << val.get_overdue_giveup_fac() << endl;
	os << "    base parameter jacobian filename = " << left << setw(20) << val.get_basejac_filename() << endl;
	if (val.get_global_opt() == PestppOptions::GLOBAL_OPT::OPT_DE)
	{
		os << "    global optimizer = differential evolution (DE)" << endl;
		os << "    DE CR = " << left << setw(10) << val.get_de_cr() << endl;
		os << "    DE F = " << left << setw(10) << val.get_de_f() << endl;
		os << "    DE population size = " << setw(10) << val.get_de_npopulation() << endl;
		os << "    DE max generations = " << setw(10) << val.get_de_max_gen() << endl;
		os << "    DE F dither = " << left << setw(10) << val.get_de_dither_f() << endl;
	}
	os << endl;
	return os;
}

PestppOptions::PestppOptions(int _n_iter_base, int _n_iter_super, int _max_n_super, double _super_eigthres,
	SVD_PACK _svd_pack, double _super_relparmax, int _max_run_fail,
	bool _iter_summary_flag, bool _der_forgive, double _overdue_reched_fac, double _overdue_giveup_fac,
	GLOBAL_OPT _global_opt, double _de_f, double _de_cr, int _de_npopulation, int _de_max_gen, bool _de_dither_f)
	: n_iter_base(_n_iter_base), n_iter_super(_n_iter_super), max_n_super(_max_n_super), super_eigthres(_super_eigthres),
	svd_pack(_svd_pack), super_relparmax(_super_relparmax),
	max_run_fail(_max_run_fail), max_super_frz_iter(50), max_reg_iter(50), base_lambda_vec({ 0.1, 1.0, 10.0, 100.0, 1000.0 }),
	lambda_scale_vec({1.0}),
	iter_summary_flag(_iter_summary_flag), der_forgive(_der_forgive), overdue_reched_fac(_overdue_reched_fac),
	overdue_giveup_fac(_overdue_giveup_fac), global_opt(_global_opt),
	de_f(_de_f), de_cr(_de_cr), de_npopulation(_de_npopulation), de_max_gen(_de_max_gen), de_dither_f(_de_dither_f)
{
}

void PestppOptions::parce_line(const string &line)
{
	string key;
	string value;
	regex lambda_reg("(\\w+)(?:\\s*\\()([^\\)]+)(?:\\))");
	const std::sregex_iterator end_reg;
	//regex lambda_reg("((\\w+\\s*\\([^\\)]+\\))+?)");
	cmatch mr;

	size_t found = line.find_first_of("#");
	if (found == string::npos) {
		found = line.length();
	}
	string tmp_line = line.substr(0, found);
	strip_ip(tmp_line, "both", "\t\n\r+ ");
	tmp_line.erase(remove(tmp_line.begin(), tmp_line.end(), '\"'),tmp_line.end());
	tmp_line.erase(remove(tmp_line.begin(), tmp_line.end(), '\''), tmp_line.end());
	//upper_ip(tmp_line);


	for (std::sregex_iterator i(tmp_line.begin(), tmp_line.end(), lambda_reg); i != end_reg; ++i)
	{
		string key = (*i)[1];
		string org_value = strip_cp((*i)[2]);
		upper_ip(key);

		string value = upper_cp(org_value);
		
		if (value.size() > 0)
			if (passed_args.find(key) != passed_args.end())
				throw PestParsingError(line, "Duplicate key word \"" + key + "\", possibly through an alias");
			passed_args.insert(key);
			arg_map[key] = value;

		if (key=="MAX_N_SUPER"){
			convert_ip(value, max_n_super);

		}
		else if ((key=="SUPER_EIGTHRESH") || (key=="SUPER_EIGTHRES"))
		{
			passed_args.insert("SUPER_EIGTHRESH");
			passed_args.insert("SUPER_EIGTHRES");
			convert_ip(value, super_eigthres);
		}
		else if (key=="N_ITER_BASE"){
			convert_ip(value, n_iter_base);
		}
		else if (key=="N_ITER_SUPER"){
			convert_ip(value, n_iter_super);
		}
		else if (key=="SVD_PACK"){

			if (value == "PROPACK")
				svd_pack = PROPACK;
			else if (value == "REDSVD")
				svd_pack = REDSVD;
			else if ((value == "EIGEN") || (value == "JACOBI"))
				svd_pack = EIGEN;
			else
				throw PestParsingError(line, "Invalid ++svd_pack: \"" + value + "\"");


		}
		
		else if (key == "SUPER_RELPARMAX"){
			convert_ip(value, super_relparmax);
		}
		else if (key == "MAX_SUPER_FRZ_ITER"){
			convert_ip(value, max_super_frz_iter);
		}
		
		else if (key == "MAX_RUN_FAIL"){
			convert_ip(value, max_run_fail);
		}
		else if (key == "MAX_REG_ITER"){
			convert_ip(value, max_reg_iter);
		}
		else if (key == "LAMBDAS")
		{
			base_lambda_vec.clear();
			vector<string> lambda_tok;
			tokenize(value, lambda_tok, ",");
			for (const auto &ilambda : lambda_tok)
			{
				base_lambda_vec.push_back(convert_cp<double>(ilambda));
			}
		}
		else if (key == "LAMBDA_SCALE_FAC")
		{
			lambda_scale_vec.clear();
			vector<string> scale_tok;
			tokenize(value, scale_tok, ",");
			for (const auto &iscale : scale_tok)
			{
				lambda_scale_vec.push_back(convert_cp<double>(iscale));
			}
		}
		else if (key == "ITERATION_SUMMARY")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> iter_summary_flag;
		}
		else if (key == "DER_FORGIVE")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> der_forgive;
		}
		else if (key == "UNCERTAINTY")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> uncert;
		}
		else if (key == "PREDICTIONS" || key == "FORECASTS")
		{
			passed_args.insert("PREDICTIONS");
			passed_args.insert("FORECASTS");
			prediction_names.clear();
			vector<string> prediction_tok;
			tokenize(value, prediction_tok, ",");
			for (const auto &pname : prediction_tok)
			{ 
				
				prediction_names.push_back(strip_cp(pname));
			}
		}
		else if ((key == "PARCOV") || (key == "PARAMETER_COVARIANCE")
			|| (key == "PARCOV_FILENAME"))
		{
			passed_args.insert("PARCOV");
			passed_args.insert("PARAMETER_COVARIANCE");
			passed_args.insert("PARCOV_FILENAME");
			//convert_ip(org_value, parcov_filename);
			parcov_filename = org_value;
		}

		else if ((key == "OBSCOV") || (key == "OBSERVATION_COVARIANCE")
			|| (key == "OBSCOV_FILENAME"))
		{
			passed_args.insert("OBSCOV");
			passed_args.insert("OBSERVATION_COVARIANCE");
			passed_args.insert("OBSCOV_FILENAME");
			//convert_ip(org_value, parcov_filename);
			obscov_filename = org_value;
		}

		else if ((key == "BASE_JACOBIAN") || (key == "BASE_JACOBIAN_FILENAME"))
		{
			passed_args.insert("BASE_JACOBIAN");
			passed_args.insert("BASE_JACOBIAN_FILENAME");
			
			//convert_ip(org_value, basejac_filename);
			basejac_filename = org_value;
		}

		else if (key == "HOTSTART_RESFILE")
		{
			//convert_ip(org_value, hotstart_resfile);
			hotstart_resfile = org_value;
		}

		else if (key == "OVERDUE_RESCHED_FAC"){
			convert_ip(value, overdue_reched_fac);
		}
		else if (key == "OVERDUE_GIVEUP_FAC"){
			convert_ip(value, overdue_giveup_fac);
		}
		else if (key == "OVERDUE_GIVEUP_MINUTES")
		{
			convert_ip(value, overdue_giveup_minutes);
		}
		else if (key == "CONDOR_SUBMIT_FILE")
		{
			//convert_ip(value, condor_submit_file);
			condor_submit_file = org_value;
		}
		else if ((key == "SWEEP_PARAMETER_CSV_FILE") || (key == "SWEEP_PAR_CSV"))
		{
			passed_args.insert("SWEEP_PARAMETER_CSV_FILE");
			passed_args.insert("SWEEP_PAR_CSV");
			
			//convert_ip(org_value, sweep_parameter_csv_file);
			sweep_parameter_csv_file = org_value;
		}
		else if ((key == "SWEEP_OUTPUT_CSV_FILE") || (key == "SWEEP_OBS_CSV"))
		{
			passed_args.insert("SWEEP_OUTPUT_CSV_FILE");
			passed_args.insert("SWEEP_OBS_CSV");
			
			//convert_ip(org_value, sweep_output_csv_file);
			sweep_output_csv_file = org_value;
		}
		else if (key == "SWEEP_CHUNK")
			convert_ip(value, sweep_chunk);
		else if (key == "SWEEP_FORGIVE")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> sweep_forgive;
		}
		else if (key == "SWEEP_BASE_RUN")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> sweep_base_run;
		}

		else if (key == "TIE_BY_GROUP")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> tie_by_group;
		}

		else if (key == "JAC_SCALE")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> jac_scale;

		}

		else if (key == "UPGRADE_AUGMENT")
		{
		cout << "++UPGRADE_AUGMENT is deprecated and no longer supported...ignoring" << endl;

		}

		else if (key == "UPGRADE_BOUNDS")
		{
		cout << "++UPGRADE_BOUNDS is deprecated and no longer supported...ignoring" << endl;

		}

		else if (key == "AUTO_NORM")
		{
			cout << "++AUTO_NORM is deprecated and no longer supported...ignoring" << endl;

		}
		else if (key == "MAT_INV")
		{
			cout << "++MAT_INV is deprecated (JtQJ only) and no longer supported...ignoring" << endl;

		}

		else if (key == "GLOBAL_OPT")
		{
			if (value == "DE") global_opt = OPT_DE;
		}
		else if (key == "DE_F")
		{
			convert_ip(value, de_f);
		}
		else if (key == "DE_CR")
		{
			convert_ip(value, de_cr);
		}
		else if (key == "DE_POP_SIZE")
		{
			convert_ip(value, de_npopulation);
		}
		else if (key == "DE_MAX_GEN")
		{
			convert_ip(value, de_max_gen);
		}
		else if (key == "DE_DITHER_F")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> de_dither_f;
		}
		else if ((key == "OPT_OBJ_FUNC") || (key == "OPT_OBJECTIVE_FUNCTION"))
		{
			passed_args.insert("OPT_OBJ_FUNC");
			passed_args.insert("OPT_OBJECTIVE_FUNCTION");
			convert_ip(value,opt_obj_func);
		}
		else if (key == "OPT_COIN_LOG")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> opt_coin_log;
		}
		else if (key == "OPT_SKIP_FINAL")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> opt_skip_final;
		}

		else if (key == "OPT_STD_WEIGHTS")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> opt_std_weights;
		}

		else if ((key == "OPT_DEC_VAR_GROUPS") || (key == "OPT_DECISION_VARIABLE_GROUPS"))
		{
			passed_args.insert("OPT_DEC_VAR_GROUPS");
			passed_args.insert("OPT_DECISION_VARIABLE_GROUPS");
			opt_dec_var_groups.clear();
			vector<string> tok;
			tokenize(value, tok, ", ");
			for (const auto &name : tok)
			{
				opt_dec_var_groups.push_back(strip_cp(name));
			}
		}

		else if ((key == "OPT_EXT_VAR_GROUPS") || (key == "OPT_EXTERNAL_VARIABLE_GROUPS"))
		{
			passed_args.insert("OPT_EXT_VAR_GROUPS");
			passed_args.insert("OPT_EXTERNAL_VARIABLE_GROUPS");
			opt_external_var_groups.clear();
			vector<string> tok;
			tokenize(value, tok, ", ");
			for (const auto &name : tok)
			{
				opt_external_var_groups.push_back(strip_cp(name));
			}
		}

		else if ((key == "OPT_CONSTRAINT_GROUPS"))
		{
			opt_constraint_groups.clear();
			vector<string> tok;
			tokenize(value, tok, ", ");
			for (const auto &name : tok)
			{
				opt_constraint_groups.push_back(strip_cp(name));
			}
		}

		else if (key == "OPT_RISK")
		{
			convert_ip(value, opt_risk);
		}

		else if (key == "OPT_ITER_DERINC_FAC")
		{
			convert_ip(value, opt_iter_derinc_fac);
		}

		else if (key == "OPT_DIRECTION")
		{
			string v;
			convert_ip(value,v);
			if (v == "MAX")
				opt_direction = -1;
			else if (v == "MIN")
				opt_direction = 1;
			else
				throw runtime_error("++opt_direction arg must be in {MAX,MIN}, not " + v);
		}

		else if (key == "OPT_ITER_TOL")
		{
			convert_ip(value, opt_iter_tol);
		}

		else if (key == "OPT_RECALC_FOSM_EVERY")
		{
			convert_ip(value, opt_recalc_fosm_every);
		}
		else if (key == "OPT_INCLUDE_BND_PI")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> opt_include_bnd_pi;
		}

		else if ((key == "IES_PAR_EN") || (key == "IES_PARAMETER_ENSEMBLE"))
		{
			passed_args.insert("IES_PARAMETER_ENSEMBLE");
			passed_args.insert("IES_PAR_EN");
			//convert_ip(value, ies_par_csv);
			ies_par_csv = org_value;
		}
		else if ((key == "IES_OBS_EN") || (key == "IES_OBSERVATION_ENSEMBLE"))
		{
			passed_args.insert("IES_OBSERVATION_ENSEMBLE");
			passed_args.insert("IES_OBS_EN");
			//convert_ip(value, ies_obs_csv);
			ies_obs_csv = org_value;
		}
		else if ((key == "IES_RESTART_PARAMETER_ENSEMBLE") || (key == "IES_RESTART_PAR_EN"))
		{
			passed_args.insert("IES_RESTART_PARAMETER_ENSEMBLE");
			passed_args.insert("IES_RESTART_PAR_EN");
			//convert_ip(value, ies_obs_restart_csv);
			ies_par_restart_csv = org_value;
		}
		else if ((key == "IES_RESTART_OBSERVATION_ENSEMBLE") || (key == "IES_RESTART_OBS_EN"))
		{
			passed_args.insert("IES_RESTART_OBSERVATION_ENSEMBLE");
			passed_args.insert("IES_RESTART_OBS_EN");
			//convert_ip(value, ies_obs_restart_csv);
			ies_obs_restart_csv = org_value;
		}

		else if ((key == "IES_USE_APPROXIMATE_SOLUTION") || (key == "IES_USE_APPROX"))
		{
			passed_args.insert("IES_USE_APPROXIMATE_SOLUTION");
			passed_args.insert("IES_USE_APPROX");
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> ies_use_approx;
		}
		else if (key == "IES_LAMBDA_MULTS")
		{
			ies_lam_mults.clear();
			vector<string> tok;
			tokenize(value, tok, ",");
			for (const auto &iscale : tok)
			{
				ies_lam_mults.push_back(convert_cp<double>(iscale));
			}
		}
		else if ((key == "IES_INIT_LAM") || (key == "IES_INITIAL_LAMBDA"))
		{
			passed_args.insert("IES_INIT_LAM");
			passed_args.insert("IES_INITIAL_LAMBDA");
			convert_ip(value, ies_init_lam);
		}
		else if (key == "IES_USE_APPROX")
		{
			convert_ip(value, ies_obs_restart_csv);
		}
		else if (key == "IES_SUBSET_SIZE")
		{
			convert_ip(value, ies_subset_size);
		}
		else if  ((key == "IES_REG_FACTOR") || (key == "IES_REG_FAC"))
		{
			passed_args.insert("IES_REG_FACTOR");
			passed_args.insert("IES_REG_FAC");
			convert_ip(value, ies_reg_factor);
		}
		else if (key == "IES_VERBOSE_LEVEL")
		{
			convert_ip(value, ies_verbose_level);
		}
		else if (key == "IES_USE_PRIOR_SCALING")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> ies_use_prior_scaling;
		}
		else if ((key == "IES_NUM_REALS") || (key == "NUM_REALS"))
		{
			passed_args.insert("IES_NUM_REALS");
			passed_args.insert("NUM_REALS");
			convert_ip(value, ies_num_reals);
			//ies_num_reals_passed = true;
		}
		else if (key == "IES_BAD_PHI")
		{
			convert_ip(value, ies_bad_phi);
		}
		else if (key == "IES_BAD_PHI_SIGMA")
		{
			convert_ip(value, ies_bad_phi_sigma);
		}
		else if ((key == "IES_INCLUDE_BASE") || (key == "IES_ADD_BASE"))
		{
			passed_args.insert("IES_INCLUDE_BASE");
			passed_args.insert("IES_ADD_BASE");
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> ies_include_base;
		}
		else if (key == "IES_USE_EMPIRICAL_PRIOR")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> ies_use_empirical_prior;
		}
		else if (key == "IES_GROUP_DRAWS")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> ies_group_draws;
		}
		else if (key == "IES_ENFORCE_BOUNDS")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> ies_enforce_bounds;
		}
		else if ((key == "IES_SAVE_BINARY") || (key == "SAVE_BINARY"))
		{
			passed_args.insert("IES_SAVE_BINARY");
			passed_args.insert("SAVE_BINARY");
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> ies_save_binary;
		}
		else if (key == "PAR_SIGMA_RANGE")
		{
			convert_ip(value, par_sigma_range);
		}
		else if (key == "YAMR_POLL_INTERVAL") {
			//doesn't apply here
		}
		else if (key == "IES_LOCALIZER")
		{
			//convert_ip(value, ies_localizer);
			ies_localizer = org_value;
		}
		else if (key == "IES_ACCEPT_PHI_FAC")
		{
			convert_ip(value, ies_accept_phi_fac);
		}
		else if (key == "IES_LAMBDA_INC_FAC")
		{
			convert_ip(value, ies_lambda_inc_fac);
		}
		else if (key == "IES_LAMBDA_DEC_FAC")
		{
			convert_ip(value, ies_lambda_dec_fac);
		}
		else if ((key == "IES_SAVE_LAMBDA_EN") || (key == "IES_SAVE_LAMBDA_ENSEMBLES"))
		{
			passed_args.insert("IES_SAVE_LAMBDA_EN");
			passed_args.insert("IES_SAVE_LAMBDA_ENSEMBLES");
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> ies_save_lambda_en;
		}
		else if ((key == "IES_WEIGHTS_EN") || (key == "IES_WEIGHTS_ENSEMBLE"))
		{
			passed_args.insert("IES_WEIGHTS_EN");
			passed_args.insert("IES_WEIGHTS_ENSEMBLE");
			ies_weight_csv = org_value;
		}
		else if (key == "IES_SUBSET_HOW")
		{
			convert_ip(value,ies_subset_how);
		}
		else if (key == "IES_LOCALIZE_HOW")
		{
			convert_ip(value, ies_localize_how);
		}
		else if (key == "IES_NUM_THREADS")
		{
			convert_ip(value, ies_num_threads);
		}
		else if (key == "IES_DEBUG_FAIL_SUBSET")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> ies_debug_fail_subset;
		}
		else if (key == "IES_DEBUG_FAIL_REMAINDER")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> ies_debug_fail_remainder;
		}
		else if (key == "IES_DEBUG_BAD_PHI")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> ies_debug_bad_phi;
		}
		else if (key == "IES_DEBUG_UPGRADE_ONLY")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> ies_debug_upgrade_only;
		}
		else if (key == "IES_DEBUG_HIGH_SUBSET_PHI")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> ies_debug_high_subset_phi;
		}
		else if (key == "IES_DEBUG_HIGH_UPGRADE_PHI")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> ies_debug_high_upgrade_phi;
		}
		else if (key == "IES_CSV_BY_REALS")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> ies_csv_by_reals;
		}
		else if (key == "IES_AUTOADALOC")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> ies_autoadaloc;
		}
		else if (key == "IES_AUTOADALOC_SIGMA_DIST")
		{
			convert_ip(value, ies_autoadaloc_sigma_dist);
		}
		else if (key == "IES_ENFORCE_CHGLIM")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> ies_enforce_chglim;
		}
		else if (key == "IES_CENTER_ON")
		{
			convert_ip(value, ies_center_on);
		}


		else if (key == "GSA_METHOD")
		{
			convert_ip(value, gsa_method);
		}
		else if (key == "GSA_MORRIS_POOLED_OBS")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> gsa_morris_pooled_obs;
		}
		else if (key == "GSA_MORRIS_OBS_SEN")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> gsa_morris_obs_sen;
		}
		else if (key == "GSA_MORRIS_P")
		{
			convert_ip(value, gsa_morris_p);
		}
		else if (key == "GSA_MORRIS_R")
		{
			convert_ip(value, gsa_morris_r);
		}
		else if (key == "GSA_MORRIS_DELTA")
		{
			convert_ip(value, gsa_morris_delta);
		}

		else if (key == "GSA_SOBOL_SAMPLES")
		{
			convert_ip(value, gsa_sobol_samples);
		}
		else if (key == "GSA_SOBOL_PAR_DIST")
		{
			convert_ip(value, gsa_sobol_par_dist);
		}

		else if (key == "GSA_SOBOL_PAR_DIST")
		{
			convert_ip(value, gsa_sobol_par_dist);
		}
		else if (key == "ENFORCE_TIED_BOUNDS")
		{
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			istringstream is(value);
			is >> boolalpha >> enforce_tied_bounds;
		 }
		else {

			throw PestParsingError(line, "Invalid key word \"" + key +"\"");
		}
	}
}

ostream& operator<< (ostream &os, const ParameterInfo& val)
{
	for(unordered_map<string, ParameterRec>::const_iterator s=val.parameter_info.begin(),
				e=val.parameter_info.end(); s!=e; ++s) {
			os << (*s).second;
		}
	return os;
}

ostream& operator<< (ostream &os, const ObservationGroupRec& val)
{
	os << "    gtarg = " << val.gtarg << endl;
	os << "    covfile = " << val.covfile << endl;
	return os;
}

ostream& operator<< (ostream &os, const ObservationRec& val)
{
	os << "    weight = " << val.weight << endl;
	os << "    group = " << val.group << endl;
	return os;
}

bool ObservationGroupRec::is_regularization(const string &grp_name)
{
	bool is_reg = false;
	int found = upper_cp(grp_name).find("REGUL");
	if (found == 0) is_reg = true;
	return is_reg;

}

void ObservationInfo::reset_group_weights(string &group, double val)
{
	for (auto &o : observations)
	{
		if (o.second.group == group)
		{
			o.second.weight = val;
		}
	}

}


void ObservationInfo::scale_group_weights(string &group, double scale_val)
{
	for (auto &o : observations)
	{
		if (o.second.group == group)
		{
			o.second.weight *= scale_val;
		}
	}

}

const ObservationRec* ObservationInfo::get_observation_rec_ptr(const string &name) const
{
	const ObservationRec *ret_val = 0;
	unordered_map<string, ObservationRec>::const_iterator p_iter;

	p_iter = observations.find(name);
	if(p_iter != observations.end()) {
		ret_val = &((*p_iter).second);
	}
	return ret_val;
}

const ObservationGroupRec* ObservationInfo::get_group_rec_ptr(const string &name) const
{
	const ObservationGroupRec *ret_val = 0;
	unordered_map<string, ObservationGroupRec>::const_iterator p_iter;

	p_iter = groups.find(name);
	if(p_iter != groups.end()) {
		ret_val = &((*p_iter).second);
	}
	return ret_val;
}

vector<string> ObservationInfo::get_groups()
{
	vector<string> ogroups;
	for (auto &g : groups)
	{
		ogroups.push_back(g.first);
	}
	return ogroups;
}

bool ObservationRec::is_regularization() const
{
	bool is_reg = false;
	int found = upper_cp(group).find("REGUL");
	if (found == 0) is_reg=true;
	return is_reg;
}

int ObservationInfo::get_nnz_obs() const
{
	int nnz = 0;
	for (auto &obs : observations)
	{
		if ((!obs.second.is_regularization()) && (obs.second.weight > 0.0))
			nnz++;
	}
	return nnz;
}

int ObservationInfo::get_nnz_obs_and_reg() const
{
	int nnz = 0;
	for (auto &obs : observations)
	{
		if (obs.second.weight > 0.0)
			nnz++;
	}
	return nnz;
}

double ObservationInfo::get_weight(const string &obs_name) const
{
	return observations.find(obs_name)->second.weight;
}

void ObservationInfo::set_weight(const string &obs_name, double &value)
{
	if (observations.find(obs_name) == observations.end())
		throw PestError("ObservationInfo::set_weight() error: observation\
			    " + obs_name + " not found");
	observations[obs_name].weight = value;
}

string ObservationInfo::get_group(const string &obs_name) const
{
	return observations.find(obs_name)->second.group;
}

bool ObservationInfo::is_regularization(const string &obs_name) const
{
	return observations.find(obs_name)->second.is_regularization();
}

Observations ObservationInfo::get_regulatization_obs(const Observations &obs_in)
{
	Observations reg_obs;

	for (Observations::const_iterator b=obs_in.begin(), e=obs_in.end();
			b!=e; ++b)
	{
		if(is_regularization((*b).first)) {
			reg_obs.insert(*b);
		}
	}
	return reg_obs;
}


ostream& operator<< (ostream &os, const ObservationInfo& val)
{
	os << "PEST Observation Information" << endl;
	os << "PEST Observation Groups" << endl;
	for(unordered_map<string, ObservationGroupRec>::const_iterator b=val.groups.begin(),
		e=val.groups.end(); b!=e; ++b) {
			os << "    name = " << (*b).first << endl;
	}
	os << "PEST Observation Information" << endl;
	for(unordered_map<string, ObservationRec>::const_iterator b=val.observations.begin(),
		e=val.observations.end(); b!=e; ++b) {
			os << "  Observation Name = " << (*b).first << endl;
			os << (*b).second;
	}
	return os;
}

ostream& operator<< (ostream &os, const SVDInfo& val)
{
	os << "PEST SVD Information" << endl;
	os << "    maxsing = " << val.maxsing << endl;
	os << "    eigthresh = " << val.eigthresh << endl;
	return os;
}

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
#include <random>
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

void ParameterGroupRec::set_defaults()
{
	name = "PARGP";
	inctyp = "RELATIVE";
	derinc = 0.01;
	derinclb = 0.0;
	forcen = "SWITCH";
	derincmul = 2.0;
	dermthd = "PARABOLIC";
	splitthresh = 1.0e-5;
	splitreldiff = 0.5;
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
const ParameterGroupRec ParameterGroupInfo::get_group_rec(const string &name) const
{
	ParameterGroupRec ret_val;
	unordered_map<string, ParameterGroupRec>::const_iterator g_iter;

	g_iter = parameter2group.find(name);
	if(g_iter != parameter2group.end()) {
		ret_val = g_iter->second;
	}
	return ret_val;
}

ParameterGroupRec* ParameterGroupInfo::get_group_rec_ptr_4_mod(const string &name)
{
	ParameterGroupRec* ret_val = 0;
	unordered_map<string, ParameterGroupRec>::iterator g_iter;

	g_iter = parameter2group.find(name);
	if (g_iter != parameter2group.end()) {
		ret_val = &g_iter->second;
	}
	return ret_val;
}

string ParameterGroupInfo::get_group_name(const string &par_name) const
{
	return get_group_rec(par_name).name;
}

void ParameterGroupInfo::insert_group(const string &group_name, ParameterGroupRec &rec)
{
	groups[group_name] = rec;
}

void ParameterGroupInfo::insert_parameter_link(const string &parameter_name, const string & group_name)
{
	unordered_map<string, ParameterGroupRec>::const_iterator g_iter;

	g_iter = groups.find(group_name);
	if(g_iter == groups.end()) {
		//throw PestIndexError(group_name, "Invalid parameter group name");
		ParameterGroupRec pgr;
		pgr.set_defaults();
		insert_group(group_name, pgr);
		g_iter = groups.find(group_name);
	}
	/*else {
		parameter2group[parameter_name] = (*g_iter).second;
	}*/
	parameter2group[parameter_name] = (*g_iter).second;

}

//const ParameterGroupInfo& ParameterGroupInfo::operator=(const ParameterGroupInfo &rhs)
//{
//	unordered_map<ParameterGroupRec*, ParameterGroupRec*> old2new;
//	unordered_map<string, ParameterGroupRec*>::const_iterator it(rhs.groups.begin());
//	unordered_map<string, ParameterGroupRec*>::const_iterator end(rhs.groups.end());
//	for (; it != end; ++it) {
//		ParameterGroupRec* new_ptr = new ParameterGroupRec(*(*it).second);
//		groups[(*it).first] = new_ptr;
//		old2new[(*it).second] = new_ptr;
//
//	}
//	unordered_map<ParameterGroupRec*, ParameterGroupRec*>::iterator it_find;
//	it = rhs.parameter2group.begin();
//	end =  rhs.parameter2group.end();
//	for (; it != end; ++it) {
//		it_find = old2new.find((*it).second);
//		if (it_find != old2new.end())
//			parameter2group[(*it).first] = (*it_find).second;
//	}
////    unordered_map<ParameterGroupRec*, ParameterGroupRec*>::iterator iit(old2new.begin());
////    unordered_map<ParameterGroupRec*, ParameterGroupRec*>::iterator eend(old2new.end());
////
////    for (; iit != eend; ++iit) {
////        delete (*iit).second;
////        delete (*iit).first;
////    }
//
//    return *this;
//}

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
		if (lower_cp(irec.second.forcen) == "switch")
		{
			switch_der = true;
			break;
		}
	}
	return switch_der;
}

ParameterGroupInfo::~ParameterGroupInfo()
{
//	unordered_map<string, ParameterGroupRec*>::iterator it(groups.begin());
//	unordered_map<string, ParameterGroupRec*>::iterator end(groups.end());
//	for (; it != end; ++it) {
//		delete (*it).second;
//	}
//
//    it = parameter2group.begin();
//	end = parameter2group.end();
//
//    for (; it != end; ++it) {
//        delete (*it).second;
//    }
}

void ParameterGroupInfo::free_mem()
{

    groups.clear();

    parameter2group.clear();

}

ostream& operator<< (ostream &os, const ParameterGroupInfo &val)
{
	unordered_map<string, ParameterGroupRec>::const_iterator it(val.parameter2group.begin());
	unordered_map<string, ParameterGroupRec>::const_iterator end(val.parameter2group.end());
	for (; it != end; ++it) {
		os << (it->second) << endl;
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

void PestppOptions::rectify_ies_da_args()
{
	for (auto& a : passed_da_ies_args)
		passed_args.emplace(a);
}

map<string,PestppOptions::ARG_STATUS> PestppOptions::parse_plusplus_line(const string& line)
{
	map<string, ARG_STATUS> arg_map;
	ARG_STATUS stat;

	map<string, string> spairs = pest_utils::parse_plusplus_line(line);
	for (auto& spair : spairs)
	{
		if (spair.second.size() == 0)
		{
			arg_map[spair.first] = ARG_STATUS::ARG_INVALID;
			continue;
		}
			
		try
		{
			stat = assign_value_by_key(spair.first, spair.second);
			arg_map[spair.first] = stat;
		}
		catch (...)
		{
			arg_map[spair.first] = ARG_STATUS::ARG_INVALID;
		}

	}
	return arg_map;
}

PestppOptions::ARG_STATUS PestppOptions::assign_value_by_key(string key, const string org_value)
{
	upper_ip(key);

	string value = upper_cp(org_value);
		
	if (value.size() == 0)
		return ARG_STATUS::ARG_INVALID;

	if (passed_args.find(key) != passed_args.end())
	{
		//throw PestParsingError(line, "Duplicate key word \"" + key + "\", possibly through an alias");
		//cout << "parse_plusplus_line() Error: Duplicate key word " << key << ", possibly through an alias" << endl;
		return ARG_STATUS::ARG_DUPLICATE;
	}
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
		{
			cout << "++SVD_PACK(PROPACK) is deprecated, resorting to REDSVD" << endl;
			svd_pack = REDSVD;			
		}
		else if (value == "REDSVD")
			svd_pack = REDSVD;
		else if ((value == "EIGEN") || (value == "JACOBI"))
			svd_pack = EIGEN;
		else
		{
			//throw PestParsingError(line, "Invalid ++svd_pack: \"" + value + "\"");
			return ARG_STATUS::ARG_INVALID;
		}
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
		iter_summary_flag = pest_utils::parse_string_arg_to_bool(value);
	}
	else if (key == "DER_FORGIVE")
	{
		der_forgive = pest_utils::parse_string_arg_to_bool(value);
	}
	else if (key == "UNCERTAINTY")
	{
		uncert = pest_utils::parse_string_arg_to_bool(value);
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
			
		//convert_ip(org_value, basejac_filename);f
		basejac_filename = org_value;
	}

	else if (key == "HOTSTART_RESFILE")
	{
		//convert_ip(org_value, hotstart_resfile);
		hotstart_resfile = org_value;
	}

	else if (key == "GLM_NUM_REALS")
	{
		convert_ip(value, glm_num_reals);
	}
	else if (key == "GLM_ACCEPT_MC_PHI")
	{
		glm_accept_mc_phi = pest_utils::parse_string_arg_to_bool(value);
	}
	else if (key == "GLM_REBASE_SUPER")
	{	
		glm_rebase_super = pest_utils::parse_string_arg_to_bool(value);
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
	else if ((key == "YAMR_POLL_INTERVAL") || (key == "PANTHER_POLL_INTERVAL"))
	{
	    passed_args.insert("YAMR_POLL_INTERVAL");
	    passed_args.insert("PANTHER_POLL_INTERVAL");
	    convert_ip(value, worker_poll_interval);
	}
	else if (key == "CONDOR_SUBMIT_FILE")
	{
		//convert_ip(value, condor_submit_file);
		condor_submit_file = org_value;
	}
	else if ((key == "SWEEP_PARAMETER_CSV_FILE") || (key == "SWEEP_PAR_CSV") || (key == "SWEEP_PARAMETER_FILE"))
	{
		passed_args.insert("SWEEP_PARAMETER_CSV_FILE");
		passed_args.insert("SWEEP_PAR_CSV");
        passed_args.insert("SWEEP_PARAMETER_FILE");


        //convert_ip(org_value, sweep_parameter_csv_file);
		sweep_parameter_csv_file = org_value;
	}
	else if ((key == "SWEEP_OUTPUT_CSV_FILE") || (key == "SWEEP_OBS_CSV") || (key == "SWEEP_OUTPUT_FILE"))
	{
		passed_args.insert("SWEEP_OUTPUT_CSV_FILE");
		passed_args.insert("SWEEP_OBS_CSV");
        passed_args.insert("SWEEP_OUTPUT_FILE");
			
		//convert_ip(org_value, sweep_output_csv_file);
		sweep_output_csv_file = org_value;
	}
	else if (key == "SWEEP_CHUNK")
		convert_ip(value, sweep_chunk);
	else if (key == "SWEEP_FORGIVE")
	{
		sweep_forgive = pest_utils::parse_string_arg_to_bool(value);
	}
	else if (key == "SWEEP_BASE_RUN")
	{
		sweep_base_run = pest_utils::parse_string_arg_to_bool(value);
	}
    else if (key == "SWEEP_INCLUDE_REGUL_PHI")
    {
        sweep_include_regul_phi = pest_utils::parse_string_arg_to_bool(value);
    }


    else if (key == "TIE_BY_GROUP")
	{
		tie_by_group = pest_utils::parse_string_arg_to_bool(value);
	}

	else if (key == "JAC_SCALE")
	{
		jac_scale = pest_utils::parse_string_arg_to_bool(value);

	}
	else if (key == "GLM_NORMAL_FORM")
	{
		if (value == "DIAG")
			glm_normal_form = GLMNormalForm::DIAG;
		else if (value == "IDENT")
			glm_normal_form = GLMNormalForm::IDENT;
		else if (value == "PRIOR")
			glm_normal_form = GLMNormalForm::PRIOR;
	}

	else if (key == "GLM_DEBUG_DER_FAIL")
	{
		glm_debug_der_fail = pest_utils::parse_string_arg_to_bool(value);
	}
	else if (key == "GLM_DEBUG_LAMB_FAIL")
	{
		glm_debug_lamb_fail = pest_utils::parse_string_arg_to_bool(value);
	}
	else if (key == "GLM_DEBUG_REAL_FAIL")
	{
		glm_debug_real_fail = pest_utils::parse_string_arg_to_bool(value);
	}
    else if (key == "GLM_DEBUG_HIGH_2ND_ITER_PHI")
    {
        glm_debug_high_2nd_iter_phi = pest_utils::parse_string_arg_to_bool(value);
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
		cout << "++MAT_INV is deprecated (JtQJ is the only form now supported) and no longer supported...ignoring" << endl;
	}

	else if (key == "GLOBAL_OPT")
	{
	if (value == "DE") global_opt = OPT_DE;
	else if (value == "MOEA") global_opt = OPT_MOEA;
	else
		throw runtime_error(value + "is not a supported global optimization option");
	}
	else if (key == "MOEA_NAME")
	{
	convert_ip(value, moea_name);
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
		de_dither_f = pest_utils::parse_string_arg_to_bool(value);
	}
	else if ((key == "OPT_OBJ_FUNC") || (key == "OPT_OBJECTIVE_FUNCTION"))
	{
		passed_args.insert("OPT_OBJ_FUNC");
		passed_args.insert("OPT_OBJECTIVE_FUNCTION");
		convert_ip(value,opt_obj_func);
        org_opt_obj_func = org_value;
	}
	else if (key == "OPT_COIN_LOG")
	{
		opt_coin_log = pest_utils::parse_string_arg_to_bool(value);
	}
	else if (key == "OPT_SKIP_FINAL")
	{
		opt_skip_final = pest_utils::parse_string_arg_to_bool(value);
	}

	else if (key == "OPT_STD_WEIGHTS")
	{
		opt_std_weights = pest_utils::parse_string_arg_to_bool(value);
	}

	else if (key == "OPT_STACK_SIZE")
	{
		convert_ip(value,opt_stack_size);
	}
	else if (key == "OPT_PAR_STACK")
	{
		opt_par_stack = org_value;
	}
	else if (key == "OPT_OBS_STACK")
	{
		opt_obs_stack = org_value;
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

	else if ((key == "OPT_RECALC_FOSM_EVERY") || (key == "OPT_RECALC_CHANCE_EVERY"))
	{
		convert_ip(value, opt_recalc_fosm_every);
	}
	else if (key == "OPT_INCLUDE_BND_PI")
	{
		opt_include_bnd_pi = pest_utils::parse_string_arg_to_bool(value);
	}

	
	
	else if (key == "GSA_METHOD")
	{
		convert_ip(value, gsa_method);
	}
	else if (key == "GSA_MORRIS_POOLED_OBS")
	{
		gsa_morris_pooled_obs = pest_utils::parse_string_arg_to_bool(value);
	}
	else if (key == "GSA_MORRIS_OBS_SEN")
	{
		gsa_morris_obs_sen = pest_utils::parse_string_arg_to_bool(value);
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
	else if (key == "ENFORCE_TIED_BOUNDS")
	{
		enforce_tied_bounds = pest_utils::parse_string_arg_to_bool(value);
	}

	else if ((key == "DEBUG_PARSE_ONLY") || (key == "PARSE_ONLY"))
	{
        passed_args.insert("DEBUG_PARSE_ONLY");
        passed_args.insert("PARSE_ONLY");
		debug_parse_only = pest_utils::parse_string_arg_to_bool(value);
	
	}
	else if (key == "PAR_SIGMA_RANGE")
	{
	convert_ip(value, par_sigma_range);
	}
    else if (key == "ENSEMBLE_OUTPUT_PRECISION")
    {
        convert_ip(value, ensemble_output_precision);
    }
	
	else if ((!assign_value_by_key_continued(key, value, org_value)) && 
	(!assign_value_by_key_sqp(key, value, org_value)) &&
	(!assign_mou_value_by_key(key, value, org_value)) && 
	(!assign_ies_value_by_key(key, value, org_value)))
	{
		//special treatment of the da args...
		if (!assign_da_value_by_key(key, value, org_value))
		{
			passed_args.insert(key);
			return ARG_STATUS::ARG_NOTFOUND;

		}
			
	}
	passed_args.insert(key);

	return ARG_STATUS::ARG_ACCEPTED;
}


bool PestppOptions::assign_ies_value_by_key(const string& key, const string& value, const string& org_value)
{
	if ((use_da_args) && (key.find("IES") != string::npos))
	{
		string da_key = "DA" + key.substr(3, key.size());
		// if the same da arg has been passed already, then we have already set 
		//the value in the ies arg for the da arg and the current arg, 
		//while valid, should not be used.
		if (passed_args.find(da_key) != passed_args.end())
			return true;
	}
	if ((key == "IES_PAR_EN") || (key == "IES_PARAMETER_ENSEMBLE"))
	{
	passed_args.insert("IES_PARAMETER_ENSEMBLE");
	passed_args.insert("IES_PAR_EN");
	//convert_ip(value, ies_par_csv);
	ies_par_csv = org_value;
	return true;
	}
	else if ((key == "IES_OBS_EN") || (key == "IES_OBSERVATION_ENSEMBLE"))
	{
	passed_args.insert("IES_OBSERVATION_ENSEMBLE");
	passed_args.insert("IES_OBS_EN");
	//convert_ip(value, ies_obs_csv);
	ies_obs_csv = org_value;
	return true;
	}
	else if ((key == "IES_RESTART_PARAMETER_ENSEMBLE") || (key == "IES_RESTART_PAR_EN"))
	{
	passed_args.insert("IES_RESTART_PARAMETER_ENSEMBLE");
	passed_args.insert("IES_RESTART_PAR_EN");
	//convert_ip(value, ies_obs_restart_csv);
	ies_par_restart_csv = org_value;
	return true;
	}
	else if ((key == "IES_RESTART_OBSERVATION_ENSEMBLE") || (key == "IES_RESTART_OBS_EN"))
	{
	passed_args.insert("IES_RESTART_OBSERVATION_ENSEMBLE");
	passed_args.insert("IES_RESTART_OBS_EN");
	//convert_ip(value, ies_obs_restart_csv);
	ies_obs_restart_csv = org_value;
	return true;
	}

    else if ((key == "IES_WEIGHTS_ENSEMBLE") || (key == "IES_WEIGHTS_EN") ||
            (key == "IES_WEIGHT_ENSEMBLE") || (key == "IES_WEIGHT_EN"))
    {
        passed_args.insert("IES_WEIGHTS_ENSEMBLE");
        passed_args.insert("IES_WEIGHTS_EN");
        passed_args.insert("IES_WEIGHT_ENSEMBLE");
        passed_args.insert("IES_WEIGHT_EN");
        //convert_ip(value, ies_obs_restart_csv);
        ies_weights_csv = org_value;
        return true;
    }

	else if ((key == "IES_USE_APPROXIMATE_SOLUTION") || (key == "IES_USE_APPROX"))
	{
	passed_args.insert("IES_USE_APPROXIMATE_SOLUTION");
	passed_args.insert("IES_USE_APPROX");
	ies_use_approx = pest_utils::parse_string_arg_to_bool(value);
	return true;
	}
	else if (key == "IES_LAMBDA_MULTS")
	{
	ies_lam_mults.clear();
	vector<string> tok;
	tokenize(value, tok, ",");
	for (const auto& iscale : tok)
	{
		ies_lam_mults.push_back(convert_cp<double>(iscale));
	}
	return true;
	}

	else if ((key == "IES_INIT_LAM") || (key == "IES_INITIAL_LAMBDA"))
	{
	passed_args.insert("IES_INIT_LAM");
	passed_args.insert("IES_INITIAL_LAMBDA");
	convert_ip(value, ies_init_lam);
	return true;
	}
	else if (key == "IES_USE_APPROX")
	{
	convert_ip(value, ies_use_approx);
	return true;
	}
	else if (key == "IES_SUBSET_SIZE")
	{
	convert_ip(value, ies_subset_size);
	return true;
	}
	else if ((key == "IES_REG_FACTOR") || (key == "IES_REG_FAC"))
	{
	passed_args.insert("IES_REG_FACTOR");
	passed_args.insert("IES_REG_FAC");
	convert_ip(value, ies_reg_factor);
	return true;
	}
	else if (key == "IES_VERBOSE_LEVEL")
	{
	convert_ip(value, ies_verbose_level);
	return true;
	}
	else if (key == "IES_USE_PRIOR_SCALING")
	{
	ies_use_prior_scaling = pest_utils::parse_string_arg_to_bool(value);
	return true;
	}

	else if (key == "IES_NUM_REALS")
	{
	convert_ip(value, ies_num_reals);
	return true;
	}
	else if (key == "IES_BAD_PHI")
	{
	convert_ip(value, ies_bad_phi);
	return true;
	}
	else if (key == "IES_BAD_PHI_SIGMA")
	{
	convert_ip(value, ies_bad_phi_sigma);
	return true;
	}
	else if ((key == "IES_INCLUDE_BASE") || (key == "IES_ADD_BASE"))
	{
	passed_args.insert("IES_INCLUDE_BASE");
	passed_args.insert("IES_ADD_BASE");
	ies_include_base = pest_utils::parse_string_arg_to_bool(value);
	return true;
	}
	else if (key == "IES_USE_EMPIRICAL_PRIOR")
	{
	ies_use_empirical_prior = pest_utils::parse_string_arg_to_bool(value);
	return true;
	}
	else if (key == "IES_GROUP_DRAWS")
	{
	ies_group_draws = pest_utils::parse_string_arg_to_bool(value);
	return true;
	}
	else if (key == "IES_ENFORCE_BOUNDS")
	{
	ies_enforce_bounds = pest_utils::parse_string_arg_to_bool(value);
	return true;
	}
	else if ((key == "IES_SAVE_BINARY") || (key == "SAVE_BINARY"))
	{
	passed_args.insert("IES_SAVE_BINARY");
	passed_args.insert("SAVE_BINARY");
	save_binary = pest_utils::parse_string_arg_to_bool(value);
	return true;
	}	
	else if (key == "IES_LOCALIZER")
	{
	//convert_ip(value, ies_localizer);
	ies_localizer = org_value;
	return true;
	}
	else if (key == "IES_ACCEPT_PHI_FAC")
	{
	convert_ip(value, ies_accept_phi_fac);
	return true;
	}
	else if (key == "IES_LAMBDA_INC_FAC")
	{
	convert_ip(value, ies_lambda_inc_fac);
	return true;
	}
	else if (key == "IES_LAMBDA_DEC_FAC")
	{
	convert_ip(value, ies_lambda_dec_fac);
	return true;
	}
	else if ((key == "IES_SAVE_LAMBDA_EN") || (key == "IES_SAVE_LAMBDA_ENSEMBLES"))
	{
	passed_args.insert("IES_SAVE_LAMBDA_EN");
	passed_args.insert("IES_SAVE_LAMBDA_ENSEMBLES");
	ies_save_lambda_en = pest_utils::parse_string_arg_to_bool(value);
	return true;
	}
	else if (key == "IES_SUBSET_HOW")
	{
	convert_ip(value, ies_subset_how);
	return true;
	}
	else if (key == "IES_LOCALIZE_HOW")
	{
	convert_ip(value, ies_localize_how);
	return true;
	}
	else if (key == "IES_NUM_THREADS")
	{
	convert_ip(value, ies_num_threads);
	return true;
	}
	else if (key == "IES_DEBUG_FAIL_SUBSET")
	{
	ies_debug_fail_subset = pest_utils::parse_string_arg_to_bool(value);
	return true;
	}
	else if (key == "IES_DEBUG_FAIL_REMAINDER")
	{
	ies_debug_fail_remainder = pest_utils::parse_string_arg_to_bool(value);
	return true;
	}
	else if (key == "IES_DEBUG_BAD_PHI")
	{
	ies_debug_bad_phi = pest_utils::parse_string_arg_to_bool(value);
	return true;
	}
	else if (key == "IES_DEBUG_UPGRADE_ONLY")
	{
	ies_debug_upgrade_only = pest_utils::parse_string_arg_to_bool(value);
	return true;
	}
	else if (key == "IES_DEBUG_HIGH_SUBSET_PHI")
	{
	ies_debug_high_subset_phi = pest_utils::parse_string_arg_to_bool(value);
	return true;
	}
	else if (key == "IES_DEBUG_HIGH_UPGRADE_PHI")
	{
	ies_debug_high_upgrade_phi = pest_utils::parse_string_arg_to_bool(value);
	return true;
	}
	else if (key == "IES_CSV_BY_REALS")
	{
	ies_csv_by_reals = pest_utils::parse_string_arg_to_bool(value);
	return true;
	}
	else if (key == "IES_AUTOADALOC")
	{
	ies_autoadaloc = pest_utils::parse_string_arg_to_bool(value);
	return true;
	}
	else if (key == "IES_AUTOADALOC_SIGMA_DIST")
	{
	convert_ip(value, ies_autoadaloc_sigma_dist);
	return true;
	}
	else if (key == "IES_ENFORCE_CHGLIM")
	{
	ies_enforce_chglim = pest_utils::parse_string_arg_to_bool(value);
	return true;
	}
	else if (key == "IES_CENTER_ON")
	{
	convert_ip(value, ies_center_on);
	return true;
	}
	else if (key == "IES_NO_NOISE")
	{
	ies_no_noise = pest_utils::parse_string_arg_to_bool(value);
	return true;
	}
	else if (key == "IES_DROP_CONFLICTS")
	{
	ies_drop_conflicts = pest_utils::parse_string_arg_to_bool(value);
	return true;
	}
	else if (key == "IES_SAVE_RESCOV")
	{
	ies_save_rescov = pest_utils::parse_string_arg_to_bool(value);
	return true;
	}
	else if (key == "IES_PDC_SIGMA_DISTANCE")
	{
	convert_ip(value, ies_pdc_sigma_distance);
	return true;
	}
	else if (key == "IES_USE_MDA")
	{
	ies_use_mda = pest_utils::parse_string_arg_to_bool(value);
	return true;
	}
	else if (key == "IES_MDA_INIT_FAC")
	{
	convert_ip(value, ies_mda_init_fac);
	return true;
	}
	else if (key == "IES_MDA_DEC_FAC")
	{
	convert_ip(value, ies_mda_dec_fac);
	return true;
	}
	else if ((key == "IES_LOC_TYPE") || (key == "IES_LOCALIZATION_TYPE"))
	{
	passed_args.insert("IES_LOC_TYPE");
	passed_args.insert("IES_LOCALIZATION_TYPE");
	convert_ip(value, ies_loc_type);
	return true;
	}
	else if (key == "IES_UPGRADES_IN_MEMORY")
	{
	ies_upgrades_in_memory = pest_utils::parse_string_arg_to_bool(value);
	return true;
	}
	else if (key == "IES_ORDERED_BINARY")
	{
	ies_ordered_binary = pest_utils::parse_string_arg_to_bool(value);
	return true;
	}
    else if (key == "IES_MULTIMODAL_ALPHA")
    {
        convert_ip(value,ies_multimodal_alpha);
        return true;
    }
    else if ((key == "IES_LOCALIZER_FORGIVE_MISSING") || (key == "IES_LOCALIZER_FORGIVE_EXTRA"))
    {
        passed_args.insert("IES_LOCALIZER_FORGIVE_MISSING");
        passed_args.insert("IES_LOCALIZER_FORGIVE_EXTRA");
        ies_localizer_forgive_missing = pest_utils::parse_string_arg_to_bool(value);
        return true;
    }
    else if (key == "IES_PHI_FACTOR_FILE")
    {
        ies_phi_fractions_file = org_value;
        return true;
    }
    else if (key == "IES_PHI_FACTORS_BY_REAL")
    {
        ies_phi_factors_by_real = pest_utils::parse_string_arg_to_bool(value);
        return true;
    }
    else if ((key == "IES_N_ITER_MEAN") || (key == "IES_N_ITER_REINFLATE"))
    {
        passed_args.insert("IES_N_ITER_MEAN");
        passed_args.insert("IES_N_ITER_REINFLATE");
        ies_n_iter_reinflate.clear();
        vector<string> tok;
        tokenize(value, tok, ",");
        for (const auto& fac : tok)
        {
            ies_n_iter_reinflate.push_back(convert_cp<int>(fac));
        }
        return true;
    }
    else if (key == "IES_REINFLATE_FACTOR")
    {
        vector<string> tok;
        tokenize(value, tok, ",");
        ies_reinflate_factor.clear();
        for (const auto& fac : tok)
        {
            ies_reinflate_factor.push_back(convert_cp<double>(fac));
        }
        return true;
    }
    else if (key == "IES_UPDATE_BY_REALS")
    {
        ies_updatebyreals = pest_utils::parse_string_arg_to_bool(value);
        return true;
    }

    else if (key == "SAVE_DENSE")
    {
        save_dense = pest_utils::parse_string_arg_to_bool(value);
        return true;
    }

    else if (key == "IES_AUTOADALOC_INDICATOR_PARS")
    {
        vector<string> tok;
        tokenize(value, tok, ",");
        ies_aal_indicator_pars.clear();
        for (const auto& fac : tok)
        {
            ies_aal_indicator_pars.push_back(upper_cp(fac));
        }
        return true;
    }

    else if (key == "IES_RUN_REALNAME")
    {
        ies_run_realname = value;
        return true;

    }



    return false;
}

bool PestppOptions::assign_da_value_by_key(const string& key, const string& value, const string& org_value)
{
	//DA parameters
	
	if (key == "DA_PARAMETER_CYCLE_TABLE")
	{
		da_par_cycle_table = org_value;
		return true;
	}
	else if (key == "DA_OBSERVATION_CYCLE_TABLE")
	{
		da_obs_cycle_table = org_value;
		return true;
	}
	else if (key == "DA_WEIGHT_CYCLE_TABLE")
	{
		da_weight_cycle_table = org_value;
		return true;
	}
	else if (key == "DA_HOTSTART_CYCLE")
	{
		convert_ip(value, da_hotstart_cycle);
		return true;
	}
	else if (key == "DA_STOP_CYCLE")
    {
        convert_ip(value, da_stop_cycle);
        return true;
    }
    else if (key == "DA_USE_SIMULATED_STATES")
    {
        da_use_simulated_states = pest_utils::parse_string_arg_to_bool(value);
        return true;
    }
    else if (key == "DA_NOPTMAX_SCHEDULE")
    {
        da_noptmax_schedule = org_value;
        return true;
    }
	//any additional da specific args must be before here!
	//some hackery!
	else
	{
		//
		//if (index != string::npos)
		set<string> ies_da_args;
		if (key.substr(0,2) == "DA")
		{
			string ies_key = "IES" + key.substr(2, key.size());
			//temp remove the ies arg if it has been found..
			bool already_found = false;
			if (passed_args.find(ies_key) != passed_args.end())
			{
				passed_args.erase(ies_key);
				already_found = true;
			}

			bool ies_accepted = assign_ies_value_by_key(ies_key, value, org_value);
			if (ies_accepted)
			{
				//passed_args.erase(ies_key);
				//passed_args.emplace(key);
				passed_da_ies_args.emplace(ies_key);
			}
			if (already_found)
			{
				passed_args.emplace(ies_key);
			}
			return ies_accepted;

		}

		return false;
	}

	
	
	
	
}


bool PestppOptions::assign_value_by_key_continued(const string& key, const string& value, const string& org_value)
{
	// This method was added as a workaround for a compiler limit of at most 128 nesting levels (MSVC); no more else if blocks could be added to assign_value_by_key()
	if (key == "PANTHER_AGENT_RESTART_ON_ERROR")
	{
		panther_agent_restart_on_error = pest_utils::parse_string_arg_to_bool(value);
		return true;
	}

	if (key == "PANTHER_AGENT_NO_PING_TIMEOUT_SECS")
	{
		convert_ip(value, panther_agent_no_ping_timeout_secs);
		return true;
	}
	else if (key == "ADDITIONAL_INS_DELIMITERS")
	{
		convert_ip(value, additional_ins_delimiters);
		return true;
	}
	else if ((key == "RANDOM_SEED") || (key == "RAND_SEED"))
	{
		passed_args.insert("RANDOM_SEED");
		passed_args.insert("RAND_SEED");
		convert_ip(value, random_seed);
		return true;
	}
	else if (key == "GLM_ITER_MC")
	{
		glm_iter_mc = pest_utils::parse_string_arg_to_bool(value);
		return true;
	}
	
	else if (key == "PANTHER_DEBUG_LOOP")
	{
		panther_debug_loop = pest_utils::parse_string_arg_to_bool(value);
		return true;
	}
	else if (key == "DEBUG_CHECK_PAR_EN_CONSISTENCY")
	{
		debug_check_paren_consistency = pest_utils::parse_string_arg_to_bool(value);
		return true;
	}	
	else if (key == "PANTHER_AGENT_FREEZE_ON_FAIL")
	{
		panther_debug_fail_freeze = pest_utils::parse_string_arg_to_bool(value);
		return true;
	}
	else if (key == "CHECK_TPLINS")
	{
		check_tplins = pest_utils::parse_string_arg_to_bool(value);
		return true;
	}
	else if (key == "FILL_TPL_ZEROS")
	{
		fill_tpl_zeros = pest_utils::parse_string_arg_to_bool(value);
		return true;
	}
    else if (key == "TPL_FORCE_DECIMAL")
    {
        tpl_force_decimal = pest_utils::parse_string_arg_to_bool(value);
        return true;
    }
	else if (key == "FORGIVE_UNKNOWN_ARGS")
	{
		forgive_unknown_args = pest_utils::parse_string_arg_to_bool(value);
		return true;
	}
	else if (key == "PANTHER_ECHO")
	{
		panther_echo = pest_utils::parse_string_arg_to_bool(value);
		return true;
	}
	else if (key == "NUM_TPL_INS_THREADS")
	{
		convert_ip(value, num_tpl_ins_threads);
		return true;
	}
	else if (key == "PANTHER_TRANSFER_ON_FINISH")
    {
        panther_transfer_on_finish.clear();
        vector<string> tokes;
        tokenize(org_value, tokes, ",");
        for (const auto& toke: tokes)
        {
            panther_transfer_on_finish.push_back(strip_cp(toke));
        }
        return true;
    }
    else if (key == "PANTHER_TRANSFER_ON_FAIL")
    {
        panther_transfer_on_fail.clear();
        vector<string> tokes;
        tokenize(org_value, tokes, ",");
        for (const auto& toke: tokes)
        {
            panther_transfer_on_fail.push_back(strip_cp(toke));
        }
        return true;
    }
    else if (key == "PANTHER_MASTER_TIMEOUT_MILLISECONDS")
    {
        convert_ip(value, panther_timeout_milliseconds);
        return true;

    }
    else if (key == "PANTHER_MASTER_ECHO_INTERVAL_MILLISECONDS")
    {
        convert_ip(value, panther_echo_interval_milliseconds);
        return true;

    }
    else if (key == "PANTHER_PERSISTENT_WORKERS")
    {
        panther_persistent_workers = pest_utils::parse_string_arg_to_bool(value);
        return true;
    }

	
	return false;
}


bool PestppOptions::assign_mou_value_by_key(const string& key, const string& value, const string& org_value)
{
	if (key == "MOU_GENERATOR")
	{
		mou_generator = org_value;
		return true;
	}

	else if (key == "MOU_POPULATION_SIZE")
	{
		convert_ip(value, mou_population_size);
		return true;
	}

	else if (key == "MOU_DV_POPULATION_FILE")
	{
		mou_dv_population_file = org_value;
		return true;
	}
	else if (key == "MOU_OBS_POPULATION_RESTART_FILE")
	{
		mou_obs_population_restart_file = org_value;
		return true;
	}

	else if (key == "MOU_OBJECTIVES")
	{
	
		mou_objectives.clear();
		vector<string> tok;
		tokenize(value, tok, ",		");
		for (const auto& obj : tok)
		{
			mou_objectives.push_back(upper_cp(strip_cp(obj)));
		}
		return true;
	}
	
	else if (key == "MOU_MAX_ARCHIVE_SIZE")
	{
		convert_ip(value, mou_max_archive_size);
		return true;
	}

	else if (key == "OPT_CHANCE_POINTS")
	{
		opt_chance_points = value;
		return true;
	}

	else if (key == "MOU_RISK_OBJECTIVE")
	{
		mou_risk_obj = pest_utils::parse_string_arg_to_bool(value);
		return true;
	}
	else if (key == "MOU_VERBOSE_LEVEL")
	{
		convert_ip(value, mou_verbose_level);
		return true;
	}
	else if (key == "MOU_ENV_SELECTOR")
	{
		mou_env_selector = upper_cp(strip_cp(value));
		return true;
	}
	else if (key == "MOU_MATING_SELECTOR")
	{
		mou_mating_selector = upper_cp(strip_cp(value));
		return true;
	}
	else if (key == "MOU_CROSSOVER_PROBABILITY")
	{
		convert_ip(value, mou_crossover_prob);
		return true;
	}
	else if (key == "MOU_MUTATION_PROBABILITY")
	{
		convert_ip(value, mou_mutation_prob);
		return true;
	}
	else if (key == "MOU_DE_F")
	{
		convert_ip(value, mou_de_f);
		return true;
	}
	else if (key == "MOU_SAVE_POPULATION_EVERY")
	{
		convert_ip(value, mou_save_population_every);
		return true;
	}
	else if (key == "MOU_PSO_OMEGA")
	{
		convert_ip(value, mou_pso_omega);
		return true;
	}
	else if (key == "MOU_PSO_SOCIAL_CONST")
	{
		mou_pso_social_const.clear();
		vector<string> tok;
		tokenize(value, tok, ",		");
		double v;
		for (const auto& t : tok)
		{
			convert_ip(t, v);
			mou_pso_social_const.push_back(v);
		}
		return true;
	}
	else if (key == "MOU_PSO_COGNITIVE_CONST")
	{
		mou_pso_cognitive_const.clear();
		vector<string> tok;
		tokenize(value, tok, ",		");
		double v;
		for (const auto& t : tok)
		{
			convert_ip(t, v);
			mou_pso_cognitive_const.push_back(v);
		}
		return true;
		}
	else if (key == "MOU_PSO_ALPHA")
	{
		convert_ip(value, mou_pso_alpha);
		return true;
		}
	else if (key == "MOU_PSO_RRAMP")
	{
		convert_ip(value, mou_pso_rramp);
		return true;
		}
	else if (key == "MOU_PSO_RFIT")
	{
		convert_ip(value, mou_pso_rfit);
		return true;
		}
	else if (key == "MOU_PSO_INERTIA")
	{
		mou_pso_inertia.clear();
		vector<string> tok;
		tokenize(value, tok, ",		");
		double v;
		for (const auto& t : tok)
		{
			convert_ip(t, v);
			mou_pso_inertia.push_back(v);
		}
		return true;
	}
	else if (key == "MOU_PSO_VMAX_FACTOR")
	{
		convert_ip(value, mou_pso_vmax_factor);
		return true;
		}
	if (key == "MOU_PSO_DV_BOUND_HANDLING")
	{
		mou_pso_dv_bound_handling = upper_cp(strip_cp(value));
		return true;
	}
	else if (key == "MOU_OUTER_REPO_OBS_FILE")
	{
		mou_outer_repo_obs_file = org_value;
		return true;
		}
	else if (key == "MOU_MAX_NN_SEARCH")
	{
		convert_ip(value, mou_max_nn_search);
		return true;
		}
	else if (key == "MOU_HYPERVOLUME_EXTREME")
	{
		convert_ip(value, mou_hypervolume_extreme);
		return true;
		}
	else if (key == "MOU_INFILL_SIZE")
	{
		convert_ip(value, mou_infill_size);
		return true;
		}
	else if (key == "MOU_RESAMPLE_EVERY")
	{
		convert_ip(value, mou_resample_every);
		return true;
		}
	else if (key == "MOU_RESAMPLE_COMMAND")
	{
		mou_resample_command = value;
		return true;
		}
	else if (key == "MOU_PPD_BETA")
	{
		convert_ip(value, mou_ppd_beta);
		return true;
	}
	else if (key == "MOU_FIT_GAMMA")
	{
		convert_ip(value, mou_fit_gamma);
		return true;
		}
	else if (key == "MOU_FIT_EPSILON")
	{
		convert_ip(value, mou_fit_epsilon);
		return true;
		}
	else if (key == "MOU_POPULATION_SCHEDULE")
    {
	    mou_population_schedule = org_value;
	    return true;
    }
    else if (key == "OPT_CHANCE_SCHEDULE")
    {
        opt_chance_schedule = org_value;
        return true;
    }

	else if (key == "MOU_SIMPLEX_REFLECTIONS")
	{
		convert_ip(value, mou_simplex_reflections);
		return true;
	}
	else if (key == "MOU_SIMPLEX_FACTORS")
	{
		mou_simplex_factors.clear();
		vector<string> tok;
		tokenize(value, tok, ",		");
		double v;
		for (const auto& t : tok)
		{
			convert_ip(t, v);
			mou_simplex_factors.push_back(v);
		}
		return true;
	}
	else if (key == "MOU_SIMPLEX_MUTATION")
    {
        mou_simplex_mutation = pest_utils::parse_string_arg_to_bool(value);
        return true;
    }
    else if (key == "MOU_USE_MULTIGEN_POPULATION")
    {
        mou_use_multigen = pest_utils::parse_string_arg_to_bool(value);
        return true;
    }
    else if (key == "MOU_SHUFFLE_FIXED_PARS")
    {
        mou_shuffle_fixed_pars = pest_utils::parse_string_arg_to_bool(value);
        return true;
    }
	else if (key == "MOU_DEBUG_DV_HANDLING")
	{
		mou_debug_dv_handling = pest_utils::parse_string_arg_to_bool(value);
		return true;
	}


	return false;
}


bool PestppOptions::assign_value_by_key_sqp(const string& key, const string& value, const string& org_value)
{
	if (key == "SQP_DV_EN")
	{
		sqp_dv_en = org_value;
		return true;
	}

	else if (key == "SQP_RESTART_OBS_EN")
	{
		sqp_obs_restart_en = org_value;
		return true;
	}
	else if (key == "SQP_NUM_REALS")
	{
		convert_ip(value, sqp_num_reals);
		return true;
	}
	else if (key == "SQP_UPDATE_HESSIAN")
	{
		sqp_update_hessian = pest_utils::parse_string_arg_to_bool(value);
		return true;
	}
	else if (key == "SQP_SCALE_FACS")
	{
		sqp_scale_facs.clear();
		vector<string> tok;
		tokenize(value, tok, ",		");
		double v;
		for (const auto& t : tok)
		{
			convert_ip(t, v);
			sqp_scale_facs.push_back(v);
		}
		return true;
	}
	
	return false;
}

void PestppOptions::summary(ostream& os) const
{

	os << endl << "    PEST++ OPTIONS: " << endl << endl;
	os << "...general options (used in multiple tools): " << endl;
	os << "svd_pack: ";
	if (svd_pack == PROPACK)
		os << "propack" << endl;
	if (svd_pack == REDSVD)
		os << "redsvd" << endl;
	if (svd_pack == EIGEN)
		os << "eigen" << endl;
	os << "lambda_scale_fac: ";
	for (auto s : lambda_scale_vec)
		os << s << ",";
	os << endl;
	os << "forgive_unknown_args: " << forgive_unknown_args << endl;
	os << "max_run_fail: " << max_run_fail << endl;
	os << "panther_poll_interval: " << worker_poll_interval << endl;
	os << "parameter_covariance: " << parcov_filename << endl;
	os << "observation_covariance: " << obscov_filename << endl;
	os << "hotstart_resfile: " << hotstart_resfile << endl;
	os << "overdue_resched_fac: " << overdue_reched_fac << endl;
	os << "overdue_giveup_fac: " << overdue_giveup_fac << endl;
	os << "overdue_giveup_minutes: " << overdue_giveup_minutes << endl;
	os << "condor_submit_file: " << condor_submit_file << endl;
	os << "tie_by_group: " << tie_by_group << endl;
	os << "par_sigma_range: " << par_sigma_range << endl;
	os << "enforce_tied_bounds: " << enforce_tied_bounds << endl;
	os << "debug_parse_only: " << debug_parse_only << endl;
	os << "check_tplins: " << check_tplins << endl;
	os << "fill_tpl_zeros: " << fill_tpl_zeros << endl;
	os << "tpl_force_decimal: " << tpl_force_decimal << endl;
	os << "additional_ins_delimiters: " << additional_ins_delimiters << endl;
	os << "random_seed: " << random_seed << endl;
	os << "num_tpl_ins_threads: " << num_tpl_ins_threads << endl;
	os << "save_binary: " << save_binary << endl;
    os << "save_dense: " << save_dense << endl;
    os << "ensemble_output_precision: " << ensemble_output_precision << endl;
	
	os << "panther_echo: " << panther_echo << endl;
	os << "panther_agent_restart_on_error: " << panther_agent_restart_on_error << endl;
	os << "panther_agent_no_ping_timeout_secs: " << panther_agent_no_ping_timeout_secs << endl;
	os << "panther_debug_loop: " << panther_debug_loop << endl;
	os << "debug_check_par_en_consistency: " << debug_check_paren_consistency << endl;
	os << "panther_agent_freeze_on_fail: " << panther_debug_fail_freeze << endl;
	os << "panther_transfer_on_finish: " << endl;
	for (auto& file : panther_transfer_on_finish)
	    os << file << endl;
    os << "panther_transfer_on_fail: " << endl;
    for (auto& file : panther_transfer_on_fail)
        os << file << endl;
    os << "panther_timeout_milliseconds: " << panther_timeout_milliseconds << endl;
    os << "panther_echo_interval_milliseconds: " << panther_echo_interval_milliseconds << endl;
    os << "panther_persistent_workers: " << panther_persistent_workers << endl;

    os << endl;

	os << "...pestpp-glm specific options:" << endl;
	os << "max_n_super: " << max_n_super << endl;
	os << "super_eigthresh: " << super_eigthres << endl;
	os << "n_iter_base: " << n_iter_base << endl;
	os << "n_iter_super: " << n_iter_super << endl;;
	os << "super_relparmax: " << super_relparmax << endl;
	os << "max_super_frz_iter: " << max_super_frz_iter << endl;	
	os << "max_reg_iter: " << max_reg_iter << endl;
	os << "lambdas: ";
	for (auto l : base_lambda_vec)
		os << l << ",";
	os << endl;
	os << "iteration_summary: " << iter_summary_flag << endl;
	os << "der_forgive: " << der_forgive << endl;
	os << "uncertainty: " << uncert << endl;
	os << "forecasts: ";
	for (auto f : prediction_names)
		os << f << ",";
	os << endl;
	os << "base_jacobian: " << basejac_filename << endl;
	os << "glm_num_reals: " << glm_num_reals << endl;
	os << "jac_scale: " << jac_scale << endl;
	string norm_str;
	if (glm_normal_form == GLMNormalForm::DIAG)
		norm_str = "DIAG";
	else if (glm_normal_form == GLMNormalForm::IDENT)
		norm_str = "IDENT";
	else if (glm_normal_form == GLMNormalForm::PRIOR)
		norm_str = "PRIOR";
	os << "glm_normal_form: " << norm_str << endl;
	os << "glm_debug_der_fail: " << glm_debug_der_fail << endl;
	os << "glm_debug_lamb_fail: " << glm_debug_lamb_fail << endl;
	os << "glm_debug_real_fail: " << glm_debug_real_fail << endl;
	os << "glm_accept_mc_phi: " << glm_accept_mc_phi << endl;
	os << "glm_rebase_super: " << glm_rebase_super << endl;
	os << "glm_iter_mc: " << glm_iter_mc << endl;
	os << "glm_high_2nd_iter_phi: " << glm_debug_high_2nd_iter_phi << endl;

//	if (global_opt == OPT_DE)
//	{
//		os << "global_opt: de" << endl;
//		os << "de_f: " << de_f << endl;
//		os << "de_cr: " << de_cr << endl;
//		os << "de_pop_size: " << de_npopulation << endl;
//		os << "de_max_gen: " << de_max_gen << endl;
//		os << "de_dither_f: " << de_dither_f << endl;
//	}
//
//	if (global_opt == OPT_MOEA)
//	{
//		os << "global_opt: MOEA" << endl;
//		os << "de_f: " << de_f << endl;
//		os << "de_cr: " << de_cr << endl;
//		os << "de_pop_size: " << de_npopulation << endl;
//		os << "de_max_gen: " << de_max_gen << endl;
//		os << "de_dither_f: " << de_dither_f << endl;
//	}

	os << endl << "...pestpp-swp options:" << endl;
	os << "sweep_parameter_csv_file: " << sweep_parameter_csv_file << endl;
	os << "sweep_output_csv_file: " << sweep_output_csv_file << endl;
	os << "sweep_chunk: " << sweep_chunk << endl;
	os << "sweep_forgive: " << sweep_forgive << endl;
	os << "sweep_base_run: " << sweep_base_run << endl;
	os << "sweep_include_regul_phi: " << sweep_include_regul_phi << endl;
	
	os << endl << "...pestpp-opt options:" << endl;
	os << "opt_objective_function: " <<  opt_obj_func << endl;
	os << "opt_coin_log: " << opt_coin_log << endl;
	os << "opt_skip_final: " << opt_skip_final << endl;
	os << "opt_std_weights: " << opt_std_weights << endl;
	os << "opt_stack_size: " << opt_stack_size << endl;
	os << "opt_par_stack: " << opt_par_stack << endl;
	os << "opt_obs_stack: " << opt_obs_stack << endl;
	os << "opt_decision_variable_groups: ";
	for (auto v : opt_dec_var_groups)
		os << v << ",";
	os << endl;
	os << "opt_external_variable_groups: ";
	for (auto v : opt_external_var_groups)
		os << v << ",";
	os << endl;
	os << "opt_constraint_groups: ";
	for (auto v : opt_constraint_groups)
		os << v << ",";
	os << endl;
	os << "opt_risk: " << opt_risk << endl;
	os << "opt_iter_derinc_fac: " << opt_iter_derinc_fac << endl;
	os << "opt_direction: " << opt_direction << endl;
	os << "opt_iter_tol: " << opt_iter_tol << endl;
	os << "opt_recalc_fosm_every: " << opt_recalc_fosm_every << endl;
	os << "opt_chance_points: " << opt_chance_points << endl;
    os << "opt_chance_schedule: " << opt_chance_schedule << endl;

	
	

	os << endl << "...pestpp-sqp options:" << endl;
	os << "sqp_dv_en: " << sqp_dv_en << endl;
	os << "sqp_obs_restart_en: " << sqp_obs_restart_en << endl;
	os << "sqp_num_reals: " << sqp_num_reals << endl;
	os << "sqp_update_hessian: " << sqp_update_hessian << endl;
	os << "sqp_scale_facs:" << endl;
	for (auto m : sqp_scale_facs)
		os << "  " << m << endl;

	os << endl << "...pestpp-mou options:" << endl;
	os << "mou_generator: " << mou_generator << endl;
	os << "mou_population_size: " << mou_population_size << endl;
	os << "mou_dv_population_file: " << mou_dv_population_file << endl;
	os << "mou_obs_population_restart_file: " << mou_obs_population_restart_file << endl;
	os << "mou_objectives: " << endl;
	for (auto obj : mou_objectives)
		os << "  " << obj << endl;
	os << "mou_max_archive_size: " << mou_max_archive_size << endl;
	os << "mou_risk_objective: " << mou_risk_obj << endl;
	os << "mou_verbose_level: " << mou_verbose_level << endl;
	os << "mou_env_selector: " << mou_env_selector << endl;
	os << "mou_mating_selector: " << mou_mating_selector << endl;
	os << "mou_crossover_prob: " << mou_crossover_prob << endl;
	os << "mou_mutation_prob: " << mou_mutation_prob << endl;
	os << "mou_de_f: " << mou_de_f << endl;
	os << "mou_save_population_every: " << mou_save_population_every << endl;
	os << "mou_pso_omega: " << mou_pso_omega << endl;
	os << "mou_pso_inertia (IINER, FINERT, INITER): " << endl;
	for (auto& f : mou_pso_inertia)
	{
		os << " " << f << " ";
	}
	os << endl;
	os << "mou_pso_social_const: " << endl;
	for (auto& f : mou_pso_social_const)
	{
		os << " " << f << endl;
	}
	os << "mou_pso_cognitive: " << endl;
	for (auto& f : mou_pso_cognitive_const)
	{
		os << " " << f << endl;
	}
	os << "mou_pso_alpha: " << mou_pso_alpha << endl;
	os << "mou_pso_rramp: " << mou_pso_rramp << endl;
	os << "mou_pso_rfit: " << mou_pso_rfit << endl;
	os << "mou_pso_vmax_factor: " << mou_pso_vmax_factor << endl;
	os << "mou_pso_dv_bound_handling: " << mou_pso_dv_bound_handling << endl;
	os << "mou_ppd_beta: " << mou_ppd_beta << endl;
	os << "mou_fit_gamma: " << mou_fit_gamma << endl;
	os << "mou_fit_epsilon: " << mou_fit_epsilon << endl;
	os << "mou_population_schedule: " << mou_population_schedule << endl;
	os << "mou_simplex_reflections:" << mou_simplex_reflections << endl;
	os << "mou_simplex_factors: " << endl;
	for (auto& f : mou_simplex_factors)
	{
		os << " " << f << endl;
	}
	os << "mou_simplex_mutation: " << mou_simplex_mutation << endl;
    os << "mou_use_multigen:" << mou_use_multigen << endl;
    os << "mou_shuffle_fixed_pars: " << mou_shuffle_fixed_pars << endl;
	os << "mou_debug_dv_handling: " << mou_debug_dv_handling << endl;

	os << endl << "...shared pestpp-ies/pestpp-da options:" << endl;
	os << "(note: 'da' args override 'ies' args when using pestpp-da)" << endl;
	os << "ies_parameter_ensemble: " << ies_par_csv << endl;
	os << "ies_observation_ensemble: " << ies_obs_csv << endl;
	os << "ies_weights_ensemble: " << ies_weights_csv << endl;
	os << "ies_restart_parameter_ensemble: " << ies_par_restart_csv << endl;
	os << "ies_restart_observation_ensemble: " << ies_obs_restart_csv << endl;
	os << "ies_use_approximate_solution: " << ies_use_approx << endl;
	os << "ies_lambda_mults: ";
	for (auto v : ies_lam_mults)
		os << v << ",";
	os << endl;
	os << "ies_initial_lambda: " << ies_init_lam << endl;
	os << "ies_use_approx: " << ies_use_approx << endl;
	os << "ies_subset_size: "<< ies_subset_size << endl;
	os << "ies_reg_factor: " << ies_reg_factor << endl;
	os << "ies_verbose_level: " << ies_verbose_level << endl;
	os << "ies_use_prior_scaling: " << ies_use_prior_scaling << endl;
	os << "ies_num_reals: " << ies_num_reals << endl;
	os << "ies_bad_phi: " << ies_bad_phi << endl;
	os << "ies_bad_phi_sigma: " << ies_bad_phi_sigma << endl;
	os << "ies_include_base: " << ies_include_base << endl;
	os << "ies_use_empirical_prior: " << ies_use_empirical_prior << endl;
	os << "ies_group_draws: " << ies_group_draws << endl;
	os << "ies_enforce_bounds: " << ies_enforce_bounds << endl;
	os << "ies_localizer: " << ies_localizer << endl;
	os << "ies_accept_phi_fac: " << ies_accept_phi_fac << endl;
	os << "ies_lambda_inc_fac: " << ies_lambda_inc_fac << endl;
	os << "ies_lambda_dec_fac: " << ies_lambda_dec_fac << endl;
	os << "ies_save_lambda_ensembles: " << ies_save_lambda_en << endl;
	os << "ies_subset_how: " << ies_subset_how << endl;
	os << "ies_localize_how: " << ies_localize_how << endl;
	os << "ies_num_threads: " << ies_num_threads << endl;
	os << "ies_debug_fail_subset: " << ies_debug_fail_subset << endl;
	os << "ies_debug_fail_remainder: " << ies_debug_fail_remainder << endl;
	os << "ies_debug_bad_phi: " << ies_debug_bad_phi << endl;
	os << "ies_debug_upgrade_only: " << ies_debug_upgrade_only << endl;
	os << "ies_debug_high_subset_phi: " << ies_debug_high_subset_phi << endl;
	os << "ies_debug_high_upgrade_phi: " << ies_debug_high_upgrade_phi << endl;
	os << "ies_csv_by_reals: " << ies_csv_by_reals << endl;
	os << "ies_autoadaloc: " << ies_autoadaloc << endl;
	os << "ies_autoadaloc_sigma_dist: " << ies_autoadaloc_sigma_dist << endl;
	os << "ies_enforce_chglim: " << ies_enforce_chglim << endl;
	os << "ies_center_on: " << ies_center_on << endl;
	os << "ies_no_noise: " << ies_no_noise << endl;
	os << "ies_drop_conflicts: " << ies_drop_conflicts << endl;
	os << "ies_save_rescov:" << ies_save_rescov << endl;
	os << "ies_pdc_sigma_distance: " << ies_pdc_sigma_distance << endl;
	os << "ies_use_mda: " << ies_use_mda << endl;
	os << "ies_mda_init_fac: " << ies_mda_init_fac << endl;
	os << "ies_mda_dec_fac: " << ies_mda_dec_fac << endl;
	os << "ies_localization_type: " << ies_loc_type << endl;
	os << "ies_upgrades_in_memory: " << ies_upgrades_in_memory << endl;
	os << "ies_ordered_binary: " << ies_ordered_binary << endl;
	os << "ies_multimodal_alpha: " << ies_multimodal_alpha << endl;
	os << "ies_localizer_forgive_extra: " << ies_localizer_forgive_missing << endl;
	os << "ies_phi_factors_file: " << ies_phi_fractions_file << endl;
    os << "ies_phi_factors_by_real: " << ies_phi_factors_by_real << endl;
    os << "ies_n_iter_reinflate: " << endl;
    for (auto v : ies_n_iter_reinflate)
        os << v << ",";
    os << endl;
    os << "ies_reinflate_factor: " << endl;
    for (auto v : ies_reinflate_factor)
        os << v << ",";
    os << endl;
    os << "ies_updatebyreals: " << ies_updatebyreals << endl;

    os << "ies_autoadaloc_indicator_pars: " << endl;
    for (auto v : ies_aal_indicator_pars)
        os << v << ",";
    os << endl;
    os << "ies_run_realname: " << ies_run_realname;

    os << endl << "pestpp-sen options: " << endl;
	os << "gsa_method: " << gsa_method << endl;
	os << "gsa_morris_pooled_obs: " << gsa_morris_pooled_obs << endl;
	os << "gsa_morris_obs_sen: " << gsa_morris_obs_sen << endl;
	os << "gsa_morris_p: " << gsa_morris_p << endl;
	os << "gsa_morris_r: " << gsa_morris_r << endl;
	os << "gsa_morris_delta: " <<  gsa_morris_delta << endl;
	os << "gsa_sobol_samples: " << gsa_sobol_samples << endl;
	os << "gsa_sobol_par_dist: " << gsa_sobol_par_dist << endl << endl;
	
	os << "pestpp-da options (those not shared with pestpp-ies):" << endl;
	os << "da_parameter_cycle_table: " << da_par_cycle_table << endl;
	os << "da_observation_cycle_table: " << da_obs_cycle_table << endl;
	os << "da_hotstart_cycle: " << da_hotstart_cycle << endl;
	os << "da_stop_cycle: " << da_stop_cycle << endl;
	os << "da_use_simulated_states: " << da_use_simulated_states << endl;
	os << "da_noptmax_schedule: " << da_noptmax_schedule << endl;

	
	os << endl << endl << endl;
}


void PestppOptions::set_defaults()
{

	set_svd_pack(PestppOptions::SVD_PACK::REDSVD);
	set_super_relparmax(0.1);
	
	set_iter_summary_flag(true);
	set_der_forgive(true);
	
	set_forgive_unknown_args(false);
	set_random_seed(358183147);
	set_base_lambda_vec(vector<double>{ 0.1, 1.0, 10.0, 100.0, 1000.0 });
	set_lambda_scale_vec(vector<double>{0.75, 1.0, 1.1});
	set_global_opt(PestppOptions::GLOBAL_OPT::NONE);
	set_de_cr(0.6);
	set_de_f(0.7);
	set_de_dither_f(true);
	set_de_npopulation(40);
	set_de_max_gen(100);
	set_n_iter_super(0);
	set_n_iter_base(1000000);
	set_super_eigthres(1.0e-6);
	set_max_n_super(1000000);
	set_max_super_frz_iter(20);
	set_max_reg_iter(20);
	set_uncert_flag(true);
	set_glm_num_reals(0);
	set_glm_normal_form(GLMNormalForm::DIAG);
	set_glm_debug_der_fail(false);
	set_glm_debug_lamb_fail(false);
	set_glm_debug_real_fail(false);
	set_glm_accept_mc_phi(false);
	set_glm_rebase_super(false);
	set_glm_iter_mc(false);
    set_glm_debug_high_2nd_iter_phi(false);
	set_prediction_names(vector<string>());
	set_parcov_filename(string());
	set_obscov_filename(string());
	set_basejac_filename(string());
    set_jac_scale(true);

    set_sweep_parameter_csv_file(string());
	set_sweep_output_csv_file("sweep_out.csv");
	set_sweep_base_run(false);
	set_sweep_forgive(false);
	set_sweep_chunk(500);
	set_tie_by_group(false);
	set_enforce_tied_bounds(false);
    set_sweep_include_regul_phi(false);

	set_opt_obj_func("");
    set_org_opt_obj_func("");
	set_opt_coin_log(true);
	set_opt_skip_final(false);
	set_opt_std_weights(false);
	set_opt_dec_var_groups(vector<string>());
	set_opt_ext_var_groups(vector<string>());
	set_opt_constraint_groups(vector<string>());
	set_opt_risk(0.5);
	set_opt_direction(1.0);
	set_opt_iter_tol(0.001);
	set_opt_recalc_fosm_every(1);
	set_opt_iter_derinc_fac(1.0);
	set_opt_include_bnd_pi(true);
	set_hotstart_resfile(string());
	set_opt_stack_size(0);
	set_opt_par_stack("");
	set_opt_obs_stack("");
	set_opt_chance_points("SINGLE");
    set_opt_chance_schedule("");
	

	set_sqp_dv_en("");
	set_sqp_obs_restart_en("");
	set_sqp_num_reals(-1);
	set_sqp_update_hessian(false);
	set_sqp_scale_facs(vector<double>{0.00001, 0.0001,0.0005, 0.001, 0.0025, 0.005, 0.01, 0.05, 0.075, 0.1, 0.25,0.5, 1.0,2.,5.,10.,});

	set_mou_generator("PSO");
	set_mou_population_size(100);
	set_mou_dv_population_file("");
	set_mou_obs_population_restart_file("");
	set_mou_objectives(vector<string>());
	set_mou_max_archive_size(500);
	set_mou_risk_obj(false);
	set_mou_verbose_level(1);
	set_mou_env_selector("NSGA");
	set_mou_mating_selector("TOURNAMENT");
	set_mou_crossover_probability(0.75);
	set_mou_mutation_probability(-999);
	set_mou_de_f(0.8);
	set_mou_save_population_every(-1);
	set_mou_pso_omega(0.7);
	set_mou_pso_cognitive_const(vector<double>{2.0});
	set_mou_pso_social_const(vector<double>{2.0});
	set_mou_pso_alpha(1.0);
	set_mou_pso_rramp(-5e+02);
	set_mou_pso_rfit(2.0);
	set_mou_pso_inertia(vector<double>{0.7, 0.4, 0});
	set_mou_pso_vmax_factor(0.8);
	set_mou_pso_dv_bound_handling("HYBRID");
	set_mou_outer_repo_obs_file("");
	set_mou_max_nn_search(get_mou_population_size());
	set_mou_hypervolume_extreme(1e+10);
	set_mou_infill_size(100);
	set_mou_ppd_beta(0.5);
	set_mou_fit_epsilon(0.05);
	set_mou_fit_gamma(0.25);
	set_mou_resample_every(-1);
	set_mou_resample_command("");
	set_mou_population_schedule("");
	set_mou_simplex_reflections(10);
	set_mou_simplex_factors(vector<double>{0.5, 0.6, 0.7, 0.8});
    set_mou_simplex_mutation(false);
    set_mou_use_multigen(false);
    set_mou_shuffle_fixed_pars(false);
	set_mou_debug_dv_handling(false);
	
	set_ies_par_csv("");
	set_ies_obs_csv("");
	set_ies_obs_restart_csv("");
	set_ies_par_restart_csv("");
    set_ies_weights_csv("");
	set_ies_lam_mults(vector<double>());
	set_ies_init_lam(0.0);
	set_ies_use_approx(true);
	set_ies_subset_size(-10);
	set_ies_reg_factor(0.0);
	set_ies_verbose_level(1);
	set_ies_use_prior_scaling(false);
	set_ies_num_reals(50);
	set_ies_bad_phi(std::numeric_limits<double>::max());
	set_ies_bad_phi_sigma(std::numeric_limits<double>::max());
	set_ies_include_base(true);
	set_ies_use_empirical_prior(false);
	set_ies_group_draws(true);
	set_ies_enforce_bounds(true);
	set_par_sigma_range(4.0);
	set_save_binary(false);
	set_ies_localizer("");
	set_ies_accept_phi_fac(1.05);
	set_ies_lambda_inc_fac(10.0);
	set_ies_lambda_dec_fac(0.75);
	set_ies_save_lambda_en(false);
	set_ies_subset_how("RANDOM");
	set_ies_localize_how("PARAMETERS");
	set_ies_num_threads(-1);
	set_ies_debug_fail_subset(false);
	set_ies_debug_fail_remainder(false);
	set_ies_debug_bad_phi(false);
	set_ies_debug_upgrade_only(false);
	set_ies_debug_high_subset_phi(false);
	set_ies_debug_high_upgrade_phi(false);
	set_ies_csv_by_reals(true);
	set_ies_autoadaloc(false);
	set_ies_autoadaloc_sigma_dist(1.0);
	set_ies_enforce_chglim(false);
	set_ies_center_on("");
	set_ies_lam_mults(vector<double>{0.1, 1.0, 10.0});
	set_ies_no_noise(false);
	set_ies_drop_conflicts(false);
	set_ies_save_rescov(false);
	set_ies_pdc_sigma_distance(-1.0);
	set_ies_use_mda(false);
	set_ies_mda_init_fac(10);
	set_ies_mda_dec_fac(0.5);
	set_ies_loc_type("LOCAL");
	set_ies_upgrades_in_memory(true);
	set_ies_ordered_binary(true);
    set_ies_multimodal_alpha(1.0);
    set_ensemble_output_precision(20);
    set_ies_localizer_forgive_missing(false);
    set_ies_phi_fractions_files("");
    set_ies_phi_factors_by_real(false);
    set_ies_n_iter_reinflate(vector < int > {0});
    set_ies_reinflate_factor(vector < double > {1.0});
    set_ies_run_realname("");

    set_ies_updatebyreals(false);
    set_save_dense(false);
    set_ies_aal_indicator_pars(vector<string>());

	// DA parameters
	//set_da_use_ies(false);
	set_da_par_cycle_table("");
	set_da_obs_cycle_table("");
	//set_da_par_csv("");
	//set_da_num_reals(50);
	set_da_hotstart_cycle(0);
    set_da_stop_cycle(1000000000);
    set_da_use_simulated_states(true);
    set_da_noptmax_schedule("");

	// End of DA parameters
	set_gsa_method("MORRIS");
	//many of these defaults are also redefined in gsa main
	set_gsa_morris_p(4);
	set_gsa_morris_r(4);
	set_gsa_morris_delta(0.6666);
	set_gsa_morris_obs_sen(true);
	set_gsa_morris_pooled_obs(false);
	set_gsa_sobol_par_dist("norm");
	set_gsa_sobol_samples(4);

	set_condor_submit_file(string());
	set_overdue_giveup_minutes(1.0e+30);
	set_overdue_reched_fac(1.15);
	set_overdue_giveup_fac(100);
	set_worker_poll_interval(1.0);
	set_max_run_fail(3);

	set_debug_parse_only(false);
	set_check_tplins(true);
	set_fill_tpl_zeros(false);
    set_tpl_force_decimal(false);
	set_additional_ins_delimiters("");
	set_num_tpl_ins_threads(1);	

	set_panther_agent_restart_on_error(false);
	set_panther_agent_no_ping_timeout_secs(-1);
	set_panther_debug_loop(false);
	set_debug_check_par_en_consistency(false);
	set_panther_debug_fail_freeze(false);
	set_panther_echo(true);
    set_panther_transfer_on_finish(vector<string>{});
    set_panther_transfer_on_fail(vector<string>{});
    set_panther_timeout_milliseconds(-999);
    set_panther_echo_interval_milliseconds(1000);
    set_panther_persistent_workers(true);

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

ObservationRec* ObservationInfo::get_observation_rec_ptr_4_mod(const string& name)
{
	ObservationRec* ret_val = 0;
	unordered_map<string, ObservationRec>::iterator p_iter;

	p_iter = observations.find(name);
	if (p_iter != observations.end()) {
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

vector<string> ObservationInfo::get_groups() const
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

void ObservationInfo::set_weight(const string &obs_name, double value)
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

PestppOptions::ARG_STATUS ControlInfo::assign_value_by_key(const string key, const string org_value, ofstream& f_rec)
{
	/*enum PestMode { ESTIMATION, REGUL, PARETO, UNKNOWN };
	double phiredswh;
	int noptmax;
	int jacfile;
	int numcom;
	double phiredstp;
	int nphistp;
	int nphinored;
	double relparstp;
	int nrelpar;
	int noptswitch;
	double splitswh;
	PestMode pestmode;*/
    stringstream ss;
	set<string> valid_args{"RSTFLE","PESTMODE","NPAR","NOBS","NPARGP","NPRIOR","NOBSGP","MAXCOMPDIM","NTPLFLE","NINSFLE",
                         "PRECIS","DPOINT","NUMCOM","JACFILE","MESSFILE","OBSREREF","RLAMBDA1","RLAMFAC","PHIREDLAM",
                         "PHIRATSUF","NUMLAM","JACUPDATE","LAMFORGIVE","DERFORGIVE","RELPARMAX","FACPARMAX",
                         "FACORIG","IBOUNDSTICK","UPVECBEND","PHIREDSWH","NOPTSWITCH","SPLITSWH","DOAUI","DOSENRESUSE",
                         "BOUNDSCALE","NOPTMAX","PHIREDSTP","NPHINORED","RELPARSTP","NRELPAR","PHISTOPTHRESH",
                         "LASTRUN","PHIABANDON","ICOV","ICOR","IEIG","IRES","JCOSAVE","VERBOSEREC","JCOSAVEITN",
                         "REISAVEITN","PARSAVEITN","PARSAVERUN"};
	string value = upper_cp(org_value);
	if (passed_args.find(key) != passed_args.end())
		return PestppOptions::ARG_STATUS::ARG_DUPLICATE;
	passed_args.insert(key);

	if (key == "NOPTMAX") 
		convert_ip(value, noptmax);
	else if (key == "RELPARMAX")
		convert_ip(value, relparmax);
	else if (key == "FACPARMAX")
		convert_ip(value, facparmax);
	else if (key == "FACORIG")
		convert_ip(value, facorig);
	else if (key == "PHIREDSTP")
		convert_ip(value, phiredstp);
	else if (key == "NPHISTP")
		convert_ip(value, nphistp);
	else if (key == "NPHINORED")
		convert_ip(value, nphinored);
	else if (key == "RELPARSTP")
		convert_ip(value, relparstp);
	else if (key == "NOPTSWITCH")
		convert_ip(value, noptswitch);
	else if (key == "SPLITSWH")
		convert_ip(value, splitswh);
    else if (key == "PHIREDSWH")
    {
        convert_ip(value,phiredswh);
    }
    else if (key == "NRELPAR")
    {
        convert_ip(value, nrelpar);
    }
	else if (key == "PESTMODE")
	{
		if (value == "ESTIMATION")
			pestmode = PestMode::ESTIMATION;
		else if ((value == "REGULARIZATION") || (value == "REGULARISATION"))
			pestmode = PestMode::REGUL;
		else if (value == "PARETO")
			pestmode = PestMode::PARETO;
		else
			return PestppOptions::ARG_STATUS::ARG_INVALID;
	}
	else if (valid_args.find(key) != valid_args.end())
    {
        ss.str("");
        ss << "WARNING: unused 'control data keyword' option,value:" << key << "," << value;
        f_rec << ss.str();
        cout << ss.str();
    }
	else
		return PestppOptions::ARG_STATUS::ARG_NOTFOUND;
	return PestppOptions::ARG_STATUS::ARG_ACCEPTED;
}

void ControlInfo::set_defaults()
{
	/*ControlInfo() : relparmax(0.0), facparmax(0.0), facorig(0.0), phiredswh(0.0), noptmax(0),
		phiredstp(0.0), nphistp(0), nphinored(0), relparstp(0.0), nrelpar(0), noptswitch(0),
		splitswh(0.0), pestmode(PestMode::ESTIMATION) {}*/

	facparmax = 1.1;
	relparmax = 1.0;
	facorig = 0.001;
	phiredswh = 0.1;
	noptmax = 0;
	phiredstp = 0.01;
	nphistp = 3;
	nphinored = 3;
	relparstp = 0.01;
	nrelpar = 3;
	noptswitch = 1;
	splitswh = 1.1;
	pestmode = PestMode::ESTIMATION;
	numcom = 1;
	jacfile = 0;

}

SVDInfo::SVDInfo()
{
	set_defaults();
}

void SVDInfo::set_defaults()
{
	maxsing = 1000000;
	eigthresh = 1.0e-6;
	eigwrite = 0;
}

PestppOptions::ARG_STATUS SVDInfo::assign_value_by_key(const string key, const string org_value)
{
	string value = upper_cp(org_value);
	if (passed_args.find(key) != passed_args.end())
		return PestppOptions::ARG_STATUS::ARG_DUPLICATE;
	passed_args.insert(key);
	if (key == "EIGTHRESH")
		convert_ip(value, eigthresh);
	else if (key == "MAXSING")
		convert_ip(value, maxsing);
	else if (key == "EIGWRITE")
		convert_ip(value, eigwrite);
	else if (key == "SVDMODE")
	{
	}
	else
		return PestppOptions::ARG_STATUS::ARG_NOTFOUND;
	return PestppOptions::ARG_STATUS::ARG_ACCEPTED;
}

double draw_standard_normal(std::mt19937& rand_gen)
{	
	using std::sqrt;
	using std::log;
	using std::cos;
	using std::sin;

	const double pi = 3.14159265358979323846264338327950288419716939937511;
	double scale = 1.0 / (rand_gen.max() - rand_gen.min() + 1.0);
	double v1 = (rand_gen() - rand_gen.min()) * scale;
	double v2 = (rand_gen() - rand_gen.min()) * scale;
	
	double r = sqrt(-2.0 * log(v1));
	return r * sin(2.0 * pi * v2);
}

vector<double> uniform_draws(int num_reals, double lower_bound, double upper_bound, std::mt19937& rand_gen)
{
	double scale = 1.0 / (rand_gen.max() - rand_gen.min() + 1.0);
	vector<double> vals;
	double bscale = upper_bound - lower_bound;
	double v1;
	for (int i = 0; i < num_reals; i++)
	{
		v1 = (rand_gen() - rand_gen.min()) * scale;
		v1 = lower_bound + (v1 * bscale);
		vals.push_back(v1);
	}
	return vals;
	
}




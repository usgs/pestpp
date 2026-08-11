/**
 * @file pest_options_registry.cpp
 * @brief The option registry: the single source of truth for the ++ options.
 * The table drives set_defaults, summary, and assign_value_by_key. The *_legacy
 * implementations are retained as an equivalence reference that self_verify() checks against.
 */
#include <string>
#include <vector>
#include <set>
#include <map>
#include <functional>
#include <sstream>
#include <iostream>
#include <limits>
#include <algorithm>
#include "pest_data_structs.h"
#include "utilities.h"

using namespace std;
using namespace pest_utils;

enum class OptType { INT, DOUBLE, BOOL, STRING, VEC_DOUBLE, VEC_INT, VEC_STRING, CUSTOM };

struct OptionSpec {
    string name;
    vector<string> aliases;
    OptType type;
    string scope;
    bool init_only;
    function<PestppOptions::ARG_STATUS(PestppOptions&, const string&, const string&)> parse;
    function<void(PestppOptions&)> apply_default;
    function<string(const PestppOptions&)> to_str;
};

const std::vector<OptionSpec>& PestppOptions::get_option_registry()
{
    static const std::vector<OptionSpec> R = {
    OptionSpec{ "MAX_N_SUPER", {}, OptType::INT, "glm", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_max_n_super(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_max_n_super(1000000); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_max_n_super()); } },
    OptionSpec{ "SUPER_EIGTHRESH", {"SUPER_EIGTHRES"}, OptType::DOUBLE, "glm", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_super_eigthres(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_super_eigthres(1.0e-6); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_super_eigthres(); return ss.str(); } },
    OptionSpec{ "N_ITER_BASE", {}, OptType::INT, "glm", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_n_iter_base(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_n_iter_base(1000000); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_n_iter_base()); } },
    OptionSpec{ "N_ITER_SUPER", {}, OptType::INT, "glm", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_n_iter_super(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_n_iter_super(0); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_n_iter_super()); } },
    OptionSpec{ "SVD_PACK", {}, OptType::CUSTOM, "general", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ if(value=="PROPACK"){cout<<"++SVD_PACK(PROPACK) is deprecated, resorting to REDSVD"<<endl;o.set_svd_pack(PestppOptions::REDSVD);}else if(value=="REDSVD")o.set_svd_pack(PestppOptions::REDSVD);else if(value=="EIGEN"||value=="JACOBI")o.set_svd_pack(PestppOptions::EIGEN);else return PestppOptions::ARG_STATUS::ARG_INVALID; return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_svd_pack(PestppOptions::SVD_PACK::REDSVD); },
        [](const PestppOptions& o)->string{ switch(o.get_svd_pack()){case PestppOptions::PROPACK:return "propack";case PestppOptions::REDSVD:return "redsvd";default:return "eigen";} } },
    OptionSpec{ "SUPER_RELPARMAX", {}, OptType::DOUBLE, "glm", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_super_relparmax(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_super_relparmax(0.1); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_super_relparmax(); return ss.str(); } },
    OptionSpec{ "MAX_SUPER_FRZ_ITER", {}, OptType::INT, "general", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_max_super_frz_iter(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_max_super_frz_iter(20); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_max_super_frz_iter()); } },
    OptionSpec{ "MAX_RUN_FAIL", {}, OptType::INT, "general", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_max_run_fail(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_max_run_fail(3); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_max_run_fail()); } },
    OptionSpec{ "MAX_REG_ITER", {}, OptType::INT, "general", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_max_reg_iter(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_max_reg_iter(20); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_max_reg_iter()); } },
    OptionSpec{ "LAMBDAS", {}, OptType::VEC_DOUBLE, "general", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ vector<double> v; vector<string> tok; tokenize(value,tok,","); for(const auto& t:tok){ v.push_back(convert_cp<double>(t)); } o.set_base_lambda_vec(v); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_base_lambda_vec(vector<double>{ 0.1, 1.0, 10.0, 100.0, 1000.0 }); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; auto v=o.get_base_lambda_vec(); for(size_t i=0;i<v.size();++i){if(i)ss<<",";ss<<v[i];} return ss.str(); } },
    OptionSpec{ "LAMBDA_SCALE_FAC", {}, OptType::VEC_DOUBLE, "general", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ vector<double> v; vector<string> tok; tokenize(value,tok,","); for(const auto& t:tok){ v.push_back(convert_cp<double>(t)); } o.set_lambda_scale_vec(v); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_lambda_scale_vec(vector<double>{0.75, 1.0, 1.1}); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; auto v=o.get_lambda_scale_vec(); for(size_t i=0;i<v.size();++i){if(i)ss<<",";ss<<v[i];} return ss.str(); } },
    // DEPRECATED. It gated only the OPENING of the .ipar/.iobj/.isen/.upg.csv files, while the
    // routines that write them were never gated, so setting it false did not suppress the
    // summaries - it made the run fail on a missing file. The files are now always prepared
    // (OutputFileWriter::prep_glm_files) and the value is recorded but not obeyed.
    OptionSpec{ "ITERATION_SUMMARY", {}, OptType::BOOL, "general", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ bool b=pest_utils::parse_string_arg_to_bool(value); if(!b) cout<<"++ITERATION_SUMMARY(false) is deprecated and no longer suppresses the summary files; they are always written"<<endl; o.set_iter_summary_flag(b); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_iter_summary_flag(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_iter_summary_flag()?1:0); } },
    OptionSpec{ "DER_FORGIVE", {}, OptType::BOOL, "general", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_der_forgive(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_der_forgive(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_der_forgive()?1:0); } },
    OptionSpec{ "UNCERTAINTY", {}, OptType::BOOL, "general", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_uncert_flag(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_uncert_flag(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_uncert_flag()?1:0); } },
    OptionSpec{ "PREDICTIONS", {"FORECASTS"}, OptType::VEC_STRING, "general", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ vector<string> v; vector<string> tok; tokenize(value,tok,","); for(const auto& t:tok){ v.push_back(strip_cp(t)); } o.set_prediction_names(v); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_prediction_names({}); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; auto v=o.get_prediction_names(); for(size_t i=0;i<v.size();++i){if(i)ss<<",";ss<<v[i];} return ss.str(); } },
    OptionSpec{ "PARCOV", {"PARAMETER_COVARIANCE","PARCOV_FILENAME"}, OptType::STRING, "general", true,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_parcov_filename(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_parcov_filename(string()); },
        [](const PestppOptions& o)->string{ return o.get_parcov_filename(); } },
    OptionSpec{ "OBSCOV", {"OBSERVATION_COVARIANCE","OBSCOV_FILENAME"}, OptType::STRING, "general", true,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_obscov_filename(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_obscov_filename(string()); },
        [](const PestppOptions& o)->string{ return o.get_obscov_filename(); } },
    OptionSpec{ "BASE_JACOBIAN", {"BASE_JACOBIAN_FILENAME"}, OptType::STRING, "general", true,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_basejac_filename(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_basejac_filename(string()); },
        [](const PestppOptions& o)->string{ return o.get_basejac_filename(); } },
    OptionSpec{ "HOTSTART_RESFILE", {}, OptType::STRING, "general", true,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_hotstart_resfile(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_hotstart_resfile(string()); },
        [](const PestppOptions& o)->string{ return o.get_hotstart_resfile(); } },
    OptionSpec{ "GLM_NUM_REALS", {}, OptType::INT, "glm", true,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_glm_num_reals(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_glm_num_reals(0); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_glm_num_reals()); } },
    OptionSpec{ "GLM_ACCEPT_MC_PHI", {}, OptType::BOOL, "glm", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_glm_accept_mc_phi(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_glm_accept_mc_phi(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_glm_accept_mc_phi()?1:0); } },
    OptionSpec{ "GLM_REBASE_SUPER", {}, OptType::BOOL, "glm", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_glm_rebase_super(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_glm_rebase_super(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_glm_rebase_super()?1:0); } },
    OptionSpec{ "OVERDUE_RESCHED_FAC", {}, OptType::DOUBLE, "general", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_overdue_reched_fac(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_overdue_reched_fac(1.15); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_overdue_reched_fac(); return ss.str(); } },
    OptionSpec{ "OVERDUE_GIVEUP_FAC", {}, OptType::DOUBLE, "general", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_overdue_giveup_fac(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_overdue_giveup_fac(100); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_overdue_giveup_fac(); return ss.str(); } },
    OptionSpec{ "OVERDUE_GIVEUP_MINUTES", {}, OptType::DOUBLE, "general", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_overdue_giveup_minutes(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_overdue_giveup_minutes(1.0e+30); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_overdue_giveup_minutes(); return ss.str(); } },
    OptionSpec{ "YAMR_POLL_INTERVAL", {"PANTHER_POLL_INTERVAL"}, OptType::DOUBLE, "panther", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_worker_poll_interval(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_worker_poll_interval(1.0); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_worker_poll_interval(); return ss.str(); } },
    OptionSpec{ "CONDOR_SUBMIT_FILE", {}, OptType::STRING, "general", false,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_condor_submit_file(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_condor_submit_file(string()); },
        [](const PestppOptions& o)->string{ return o.get_condor_submit_file(); } },
    OptionSpec{ "SWEEP_PARAMETER_CSV_FILE", {"SWEEP_PAR_CSV","SWEEP_PARAMETER_FILE"}, OptType::STRING, "sweep", false,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_sweep_parameter_csv_file(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sweep_parameter_csv_file(string()); },
        [](const PestppOptions& o)->string{ return o.get_sweep_parameter_csv_file(); } },
    OptionSpec{ "SWEEP_OUTPUT_CSV_FILE", {"SWEEP_OBS_CSV","SWEEP_OUTPUT_FILE"}, OptType::STRING, "sweep", false,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_sweep_output_csv_file(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sweep_output_csv_file("sweep_out.csv"); },
        [](const PestppOptions& o)->string{ return o.get_sweep_output_csv_file(); } },
    OptionSpec{ "SWEEP_CHUNK", {}, OptType::INT, "sweep", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_sweep_chunk(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sweep_chunk(500); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_sweep_chunk()); } },
    OptionSpec{ "SWEEP_FORGIVE", {}, OptType::BOOL, "sweep", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_sweep_forgive(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sweep_forgive(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_sweep_forgive()?1:0); } },
    OptionSpec{ "SWEEP_BASE_RUN", {}, OptType::BOOL, "sweep", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_sweep_base_run(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sweep_base_run(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_sweep_base_run()?1:0); } },
    OptionSpec{ "SWEEP_INCLUDE_REGUL_PHI", {}, OptType::BOOL, "sweep", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_sweep_include_regul_phi(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sweep_include_regul_phi(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_sweep_include_regul_phi()?1:0); } },
    OptionSpec{ "TIE_BY_GROUP", {}, OptType::BOOL, "general", true,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_tie_by_group(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_tie_by_group(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_tie_by_group()?1:0); } },
    OptionSpec{ "JAC_SCALE", {}, OptType::BOOL, "general", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_jac_scale(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_jac_scale(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_jac_scale()?1:0); } },
    OptionSpec{ "GLM_NORMAL_FORM", {}, OptType::CUSTOM, "glm", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ if(value=="DIAG")o.set_glm_normal_form(PestppOptions::GLMNormalForm::DIAG);else if(value=="IDENT")o.set_glm_normal_form(PestppOptions::GLMNormalForm::IDENT);else if(value=="PRIOR")o.set_glm_normal_form(PestppOptions::GLMNormalForm::PRIOR);else if(value=="HP")o.set_glm_normal_form(PestppOptions::GLMNormalForm::HP); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_glm_normal_form(PestppOptions::GLMNormalForm::DIAG); },
        [](const PestppOptions& o)->string{ switch(o.get_glm_normal_form()){case PestppOptions::GLMNormalForm::IDENT:return "ident";case PestppOptions::GLMNormalForm::DIAG:return "diag";case PestppOptions::GLMNormalForm::PRIOR:return "prior";default:return "hp";} } },
    OptionSpec{ "GLM_DEBUG_DER_FAIL", {}, OptType::BOOL, "glm", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_glm_debug_der_fail(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_glm_debug_der_fail(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_glm_debug_der_fail()?1:0); } },
    OptionSpec{ "GLM_DEBUG_LAMB_FAIL", {}, OptType::BOOL, "glm", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_glm_debug_lamb_fail(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_glm_debug_lamb_fail(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_glm_debug_lamb_fail()?1:0); } },
    OptionSpec{ "GLM_DEBUG_REAL_FAIL", {}, OptType::BOOL, "glm", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_glm_debug_real_fail(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_glm_debug_real_fail(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_glm_debug_real_fail()?1:0); } },
    OptionSpec{ "GLM_DEBUG_HIGH_2ND_ITER_PHI", {}, OptType::BOOL, "glm", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_glm_debug_high_2nd_iter_phi(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_glm_debug_high_2nd_iter_phi(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_glm_debug_high_2nd_iter_phi()?1:0); } },
    OptionSpec{ "GLM_HP_LAMBDAS", {}, OptType::BOOL, "glm", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_glm_hp_lambdas(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_glm_hp_lambdas(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_glm_hp_lambdas()?1:0); } },
    OptionSpec{ "GLM_PANTHER_LAMBDAS", {}, OptType::BOOL, "glm", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_glm_panther_lambdas(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_glm_panther_lambdas(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_glm_panther_lambdas()?1:0); } },
    OptionSpec{ "UPGRADE_AUGMENT", {}, OptType::CUSTOM, "general", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ cout<<"++UPGRADE_AUGMENT is deprecated and no longer supported...ignoring"<<endl; return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions&){},
        [](const PestppOptions& o)->string{ return "(deprecated)"; } },
    OptionSpec{ "UPGRADE_BOUNDS", {}, OptType::CUSTOM, "general", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ cout<<"++UPGRADE_BOUNDS is deprecated and no longer supported...ignoring"<<endl; return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions&){},
        [](const PestppOptions& o)->string{ return "(deprecated)"; } },
    OptionSpec{ "AUTO_NORM", {}, OptType::CUSTOM, "general", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ cout<<"++AUTO_NORM is deprecated and no longer supported...ignoring"<<endl; return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions&){},
        [](const PestppOptions& o)->string{ return "(deprecated)"; } },
    OptionSpec{ "MAT_INV", {}, OptType::CUSTOM, "general", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ cout<<"++MAT_INV is deprecated (JtQJ is the only form now supported) and no longer supported...ignoring"<<endl; return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions&){},
        [](const PestppOptions& o)->string{ return "(deprecated)"; } },
    OptionSpec{ "GLOBAL_OPT", {}, OptType::CUSTOM, "general", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ if(value=="DE")o.set_global_opt(PestppOptions::OPT_DE);else if(value=="MOEA")o.set_global_opt(PestppOptions::OPT_MOEA);else throw runtime_error(value+"is not a supported global optimization option"); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_global_opt(PestppOptions::GLOBAL_OPT::NONE); },
        [](const PestppOptions& o)->string{ switch(o.get_global_opt()){case PestppOptions::OPT_DE:return "de";case PestppOptions::OPT_MOEA:return "moea";default:return "none";} } },
    OptionSpec{ "MOEA_NAME", {}, OptType::CUSTOM, "general", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ convert_ip(value,o.moea_name); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.moea_name=string(); },
        [](const PestppOptions& o)->string{ return o.moea_name; } },
    OptionSpec{ "DE_F", {}, OptType::DOUBLE, "glm", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_de_f(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_de_f(0.7); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_de_f(); return ss.str(); } },
    OptionSpec{ "DE_CR", {}, OptType::DOUBLE, "glm", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_de_cr(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_de_cr(0.6); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_de_cr(); return ss.str(); } },
    OptionSpec{ "DE_POP_SIZE", {}, OptType::INT, "glm", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_de_npopulation(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_de_npopulation(40); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_de_npopulation()); } },
    OptionSpec{ "DE_MAX_GEN", {}, OptType::INT, "glm", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_de_max_gen(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_de_max_gen(100); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_de_max_gen()); } },
    OptionSpec{ "DE_DITHER_F", {}, OptType::BOOL, "glm", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_de_dither_f(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_de_dither_f(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_de_dither_f()?1:0); } },
    OptionSpec{ "OPT_OBJ_FUNC", {"OPT_OBJECTIVE_FUNCTION"}, OptType::CUSTOM, "opt", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ convert_ip(value,o.opt_obj_func); o.org_opt_obj_func=org_value; return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_opt_obj_func(string()); o.set_org_opt_obj_func(string()); },
        [](const PestppOptions& o)->string{ return o.get_opt_obj_func(); } },
    OptionSpec{ "OPT_COIN_LOG", {}, OptType::BOOL, "opt", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_opt_coin_log(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_opt_coin_log(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_opt_coin_log()?1:0); } },
    OptionSpec{ "OPT_SKIP_FINAL", {}, OptType::BOOL, "opt", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_opt_skip_final(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_opt_skip_final(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_opt_skip_final()?1:0); } },
    OptionSpec{ "OPT_STD_WEIGHTS", {}, OptType::BOOL, "opt", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_opt_std_weights(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_opt_std_weights(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_opt_std_weights()?1:0); } },
    OptionSpec{ "OPT_STACK_SIZE", {}, OptType::INT, "opt", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_opt_stack_size(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_opt_stack_size(0); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_opt_stack_size()); } },
    OptionSpec{ "OPT_PAR_STACK", {}, OptType::STRING, "opt", true,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_opt_par_stack(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_opt_par_stack(""); },
        [](const PestppOptions& o)->string{ return o.get_opt_par_stack(); } },
    OptionSpec{ "OPT_OBS_STACK", {}, OptType::STRING, "opt", true,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_opt_obs_stack(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_opt_obs_stack(""); },
        [](const PestppOptions& o)->string{ return o.get_opt_obs_stack(); } },
    OptionSpec{ "OPT_DEC_VAR_GROUPS", {"OPT_DECISION_VARIABLE_GROUPS"}, OptType::VEC_STRING, "opt", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ vector<string> v; vector<string> tok; tokenize(value,tok,", "); for(const auto& t:tok){ v.push_back(strip_cp(t)); } o.set_opt_dec_var_groups(v); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_opt_dec_var_groups({}); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; auto v=o.get_opt_dec_var_groups(); for(size_t i=0;i<v.size();++i){if(i)ss<<",";ss<<v[i];} return ss.str(); } },
    OptionSpec{ "OPT_EXT_VAR_GROUPS", {"OPT_EXTERNAL_VARIABLE_GROUPS"}, OptType::VEC_STRING, "opt", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ vector<string> v; vector<string> tok; tokenize(value,tok,", "); for(const auto& t:tok){ v.push_back(strip_cp(t)); } o.set_opt_ext_var_groups(v); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_opt_ext_var_groups({}); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; auto v=o.get_opt_ext_var_groups(); for(size_t i=0;i<v.size();++i){if(i)ss<<",";ss<<v[i];} return ss.str(); } },
    OptionSpec{ "OPT_CONSTRAINT_GROUPS", {}, OptType::VEC_STRING, "opt", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ vector<string> v; vector<string> tok; tokenize(value,tok,", "); for(const auto& t:tok){ v.push_back(strip_cp(t)); } o.set_opt_constraint_groups(v); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_opt_constraint_groups({}); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; auto v=o.get_opt_constraint_groups(); for(size_t i=0;i<v.size();++i){if(i)ss<<",";ss<<v[i];} return ss.str(); } },
    OptionSpec{ "OPT_USE_ROBUST", {}, OptType::BOOL, "opt", true,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_opt_use_robust(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_opt_use_robust(false); },
        [](const PestppOptions& o)->string{ return o.get_opt_use_robust() ? "1" : "0"; } },
    OptionSpec{ "OPT_RISK", {}, OptType::DOUBLE, "opt", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_opt_risk(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_opt_risk(0.5); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_opt_risk(); return ss.str(); } },
    OptionSpec{ "OPT_ITER_DERINC_FAC", {}, OptType::DOUBLE, "opt", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_opt_iter_derinc_fac(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_opt_iter_derinc_fac(1.0); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_opt_iter_derinc_fac(); return ss.str(); } },
    OptionSpec{ "OPT_DIRECTION", {}, OptType::CUSTOM, "opt", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ string v; convert_ip(value,v); if(v=="MAX")o.set_opt_direction(-1); else if(v=="MIN")o.set_opt_direction(1); else throw runtime_error("++opt_direction arg must be in {MAX,MIN}, not "+v); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_opt_direction(1.0); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_opt_direction()); } },
    OptionSpec{ "OPT_ITER_TOL", {}, OptType::DOUBLE, "opt", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_opt_iter_tol(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_opt_iter_tol(0.001); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_opt_iter_tol(); return ss.str(); } },
    OptionSpec{ "OPT_RECALC_FOSM_EVERY", {"OPT_RECALC_CHANCE_EVERY"}, OptType::INT, "opt", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_opt_recalc_fosm_every(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_opt_recalc_fosm_every(1); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_opt_recalc_fosm_every()); } },
    OptionSpec{ "OPT_INCLUDE_BND_PI", {}, OptType::BOOL, "opt", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_opt_include_bnd_pi(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_opt_include_bnd_pi(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_opt_include_bnd_pi()?1:0); } },
    OptionSpec{ "GSA_METHOD", {}, OptType::STRING, "gsa", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_gsa_method(value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_gsa_method("MORRIS"); },
        [](const PestppOptions& o)->string{ return o.get_gsa_method(); } },
    OptionSpec{ "GSA_MORRIS_POOLED_OBS", {}, OptType::BOOL, "gsa", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_gsa_morris_pooled_obs(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_gsa_morris_pooled_obs(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_gsa_morris_pooled_obs()?1:0); } },
    OptionSpec{ "GSA_MORRIS_OBS_SEN", {}, OptType::BOOL, "gsa", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_gsa_morris_obs_sen(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_gsa_morris_obs_sen(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_gsa_morris_obs_sen()?1:0); } },
    OptionSpec{ "GSA_MORRIS_P", {}, OptType::INT, "gsa", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_gsa_morris_p(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_gsa_morris_p(4); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_gsa_morris_p()); } },
    OptionSpec{ "GSA_MORRIS_R", {}, OptType::INT, "gsa", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_gsa_morris_r(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_gsa_morris_r(4); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_gsa_morris_r()); } },
    OptionSpec{ "GSA_MORRIS_DELTA", {}, OptType::DOUBLE, "gsa", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_gsa_morris_delta(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_gsa_morris_delta(0.6666); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_gsa_morris_delta(); return ss.str(); } },
    OptionSpec{ "GSA_SOBOL_SAMPLES", {}, OptType::INT, "gsa", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_gsa_sobol_samples(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_gsa_sobol_samples(4); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_gsa_sobol_samples()); } },
    OptionSpec{ "GSA_SOBOL_PAR_DIST", {}, OptType::STRING, "gsa", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_gsa_sobol_par_dist(value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_gsa_sobol_par_dist("norm"); },
        [](const PestppOptions& o)->string{ return o.get_gsa_sobol_par_dist(); } },
    OptionSpec{ "ENFORCE_TIED_BOUNDS", {}, OptType::BOOL, "general", true,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_enforce_tied_bounds(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_enforce_tied_bounds(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_enforce_tied_bounds()?1:0); } },
    OptionSpec{ "DEBUG_PARSE_ONLY", {"PARSE_ONLY"}, OptType::BOOL, "general", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_debug_parse_only(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_debug_parse_only(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_debug_parse_only()?1:0); } },
    OptionSpec{ "PAR_SIGMA_RANGE", {}, OptType::DOUBLE, "general", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_par_sigma_range(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_par_sigma_range(4.0); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_par_sigma_range(); return ss.str(); } },
    OptionSpec{ "ENSEMBLE_OUTPUT_PRECISION", {}, OptType::INT, "general", true,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_ensemble_output_precision(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ensemble_output_precision(20); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ensemble_output_precision()); } },
    OptionSpec{ "IES_PAR_EN", {"IES_PARAMETER_ENSEMBLE"}, OptType::STRING, "ies", true,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_ies_par_csv(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_par_csv(""); },
        [](const PestppOptions& o)->string{ return o.get_ies_par_csv(); } },
    OptionSpec{ "IES_OBS_EN", {"IES_OBSERVATION_ENSEMBLE"}, OptType::STRING, "ies", true,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_ies_obs_csv(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_obs_csv(""); },
        [](const PestppOptions& o)->string{ return o.get_ies_obs_csv(); } },
    OptionSpec{ "IES_RESTART_PARAMETER_ENSEMBLE", {"IES_RESTART_PAR_EN"}, OptType::STRING, "ies", true,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_ies_par_restart_csv(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_par_restart_csv(""); },
        [](const PestppOptions& o)->string{ return o.get_ies_par_restart_csv(); } },
    OptionSpec{ "IES_RESTART_OBSERVATION_ENSEMBLE", {"IES_RESTART_OBS_EN"}, OptType::STRING, "ies", true,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_ies_obs_restart_csv(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_obs_restart_csv(""); },
        [](const PestppOptions& o)->string{ return o.get_ies_obs_restart_csv(); } },
    OptionSpec{ "IES_WEIGHTS_ENSEMBLE", {"IES_WEIGHTS_EN","IES_WEIGHT_ENSEMBLE","IES_WEIGHT_EN"}, OptType::STRING, "ies", true,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_ies_weights_csv(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_weights_csv(""); },
        [](const PestppOptions& o)->string{ return o.get_ies_weights_csv(); } },
    OptionSpec{ "IES_USE_APPROXIMATE_SOLUTION", {"IES_USE_APPROX"}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_use_approx(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_use_approx(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_use_approx()?1:0); } },
    OptionSpec{ "IES_LAMBDA_MULTS", {}, OptType::VEC_DOUBLE, "ies", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ vector<double> v; vector<string> tok; tokenize(value,tok,","); for(const auto& t:tok){ v.push_back(convert_cp<double>(t)); } o.set_ies_lam_mults(v); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_lam_mults(vector<double>{0.1, 1.0, 10.0}); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; auto v=o.get_ies_lam_mults(); for(size_t i=0;i<v.size();++i){if(i)ss<<",";ss<<v[i];} return ss.str(); } },
    OptionSpec{ "IES_INIT_LAM", {"IES_INITIAL_LAMBDA"}, OptType::DOUBLE, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_ies_init_lam(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_init_lam(0.0); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_ies_init_lam(); return ss.str(); } },
    OptionSpec{ "IES_SUBSET_SIZE", {}, OptType::INT, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_ies_subset_size(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_subset_size(-10); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_subset_size()); } },
    OptionSpec{ "IES_REG_FACTOR", {"IES_REG_FAC"}, OptType::DOUBLE, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_ies_reg_factor(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_reg_factor(0.0); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_ies_reg_factor(); return ss.str(); } },
    OptionSpec{ "IES_VERBOSE_LEVEL", {}, OptType::INT, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_ies_verbose_level(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_verbose_level(1); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_verbose_level()); } },
    OptionSpec{ "IES_USE_PRIOR_SCALING", {}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_use_prior_scaling(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_use_prior_scaling(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_use_prior_scaling()?1:0); } },
    OptionSpec{ "IES_NUM_REALS", {}, OptType::INT, "ies", true,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_ies_num_reals(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_num_reals(50); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_num_reals()); } },
    OptionSpec{ "IES_BAD_PHI", {}, OptType::DOUBLE, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_ies_bad_phi(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_bad_phi(std::numeric_limits<double>::max()); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_ies_bad_phi(); return ss.str(); } },
    OptionSpec{ "PREEMPTION_POLL_INTERVAL_MINUTES", {}, OptType::DOUBLE, "", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_preemption_poll_interval_minutes(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_preemption_poll_interval_minutes(0.0); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_preemption_poll_interval_minutes(); return ss.str(); } },
    OptionSpec{ "IES_BAD_PHI_SIGMA", {}, OptType::DOUBLE, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_ies_bad_phi_sigma(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_bad_phi_sigma(std::numeric_limits<double>::max()); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_ies_bad_phi_sigma(); return ss.str(); } },
    OptionSpec{ "IES_INCLUDE_BASE", {"IES_ADD_BASE"}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_include_base(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_include_base(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_include_base()?1:0); } },
    OptionSpec{ "IES_USE_EMPIRICAL_PRIOR", {}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_use_empirical_prior(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_use_empirical_prior(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_use_empirical_prior()?1:0); } },
    OptionSpec{ "IES_GROUP_DRAWS", {}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_group_draws(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_group_draws(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_group_draws()?1:0); } },
    OptionSpec{ "IES_ENFORCE_BOUNDS", {}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_enforce_bounds(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_enforce_bounds(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_enforce_bounds()?1:0); } },
    OptionSpec{ "IES_SAVE_BINARY", {"SAVE_BINARY"}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_save_binary(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_save_binary(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_save_binary()?1:0); } },
    OptionSpec{ "IES_LOCALIZER", {}, OptType::STRING, "ies", true,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_ies_localizer(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_localizer(""); },
        [](const PestppOptions& o)->string{ return o.get_ies_localizer(); } },
    OptionSpec{ "IES_ACCEPT_PHI_FAC", {}, OptType::DOUBLE, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_ies_accept_phi_fac(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_accept_phi_fac(1.05); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_ies_accept_phi_fac(); return ss.str(); } },
    OptionSpec{ "IES_LAMBDA_INC_FAC", {}, OptType::DOUBLE, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_ies_lambda_inc_fac(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_lambda_inc_fac(10.0); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_ies_lambda_inc_fac(); return ss.str(); } },
    OptionSpec{ "IES_LAMBDA_DEC_FAC", {}, OptType::DOUBLE, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_ies_lambda_dec_fac(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_lambda_dec_fac(0.75); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_ies_lambda_dec_fac(); return ss.str(); } },
    OptionSpec{ "IES_SAVE_LAMBDA_EN", {"IES_SAVE_LAMBDA_ENSEMBLES"}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_save_lambda_en(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_save_lambda_en(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_save_lambda_en()?1:0); } },
    OptionSpec{ "IES_SUBSET_HOW", {}, OptType::STRING, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_subset_how(value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_subset_how("RANDOM"); },
        [](const PestppOptions& o)->string{ return o.get_ies_subset_how(); } },
    OptionSpec{ "IES_LOCALIZE_HOW", {}, OptType::STRING, "ies", true,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_localize_how(value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_localize_how("PARAMETERS"); },
        [](const PestppOptions& o)->string{ return o.get_ies_localize_how(); } },
    OptionSpec{ "IES_NUM_THREADS", {}, OptType::INT, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_ies_num_threads(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_num_threads(-1); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_num_threads()); } },
    OptionSpec{ "IES_DEBUG_FAIL_SUBSET", {}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_debug_fail_subset(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_debug_fail_subset(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_debug_fail_subset()?1:0); } },
    OptionSpec{ "IES_DEBUG_FAIL_REMAINDER", {}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_debug_fail_remainder(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_debug_fail_remainder(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_debug_fail_remainder()?1:0); } },
    OptionSpec{ "IES_DEBUG_BAD_PHI", {}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_debug_bad_phi(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_debug_bad_phi(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_debug_bad_phi()?1:0); } },
    OptionSpec{ "IES_DEBUG_UPGRADE_ONLY", {}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_debug_upgrade_only(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_debug_upgrade_only(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_debug_upgrade_only()?1:0); } },
    OptionSpec{ "IES_DEBUG_HIGH_SUBSET_PHI", {}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_debug_high_subset_phi(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_debug_high_subset_phi(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_debug_high_subset_phi()?1:0); } },
    OptionSpec{ "IES_DEBUG_HIGH_UPGRADE_PHI", {}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_debug_high_upgrade_phi(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_debug_high_upgrade_phi(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_debug_high_upgrade_phi()?1:0); } },
    OptionSpec{ "IES_CSV_BY_REALS", {}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_csv_by_reals(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_csv_by_reals(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_csv_by_reals()?1:0); } },
    OptionSpec{ "IES_AUTOADALOC", {}, OptType::BOOL, "ies", true,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_autoadaloc(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_autoadaloc(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_autoadaloc()?1:0); } },
    OptionSpec{ "IES_AUTOADALOC_SIGMA_DIST", {}, OptType::DOUBLE, "ies", true,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_ies_autoadaloc_sigma_dist(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_autoadaloc_sigma_dist(1.0); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_ies_autoadaloc_sigma_dist(); return ss.str(); } },
    OptionSpec{ "IES_ENFORCE_CHGLIM", {}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_enforce_chglim(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_enforce_chglim(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_enforce_chglim()?1:0); } },
    OptionSpec{ "IES_CENTER_ON", {}, OptType::STRING, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_center_on(value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_center_on(""); },
        [](const PestppOptions& o)->string{ return o.get_ies_center_on(); } },
    OptionSpec{ "IES_NO_NOISE", {}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_no_noise(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_no_noise(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_no_noise()?1:0); } },
    OptionSpec{ "IES_DROP_CONFLICTS", {}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_drop_conflicts(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_drop_conflicts(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_drop_conflicts()?1:0); } },
    OptionSpec{ "IES_SAVE_RESCOV", {}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_save_rescov(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_save_rescov(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_save_rescov()?1:0); } },
    OptionSpec{ "IES_PDC_SIGMA_DISTANCE", {}, OptType::DOUBLE, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_ies_pdc_sigma_distance(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_pdc_sigma_distance(-1.0); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_ies_pdc_sigma_distance(); return ss.str(); } },
    OptionSpec{ "IES_USE_MDA", {}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_use_mda(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_use_mda(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_use_mda()?1:0); } },
    OptionSpec{ "IES_MDA_INIT_FAC", {}, OptType::DOUBLE, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_ies_mda_init_fac(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_mda_init_fac(10); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_ies_mda_init_fac(); return ss.str(); } },
    OptionSpec{ "IES_MDA_DEC_FAC", {}, OptType::DOUBLE, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_ies_mda_dec_fac(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_mda_dec_fac(0.5); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_ies_mda_dec_fac(); return ss.str(); } },
    OptionSpec{ "IES_LOC_TYPE", {"IES_LOCALIZATION_TYPE"}, OptType::STRING, "ies", true,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_loc_type(value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_loc_type("LOCAL"); },
        [](const PestppOptions& o)->string{ return o.get_ies_loc_type(); } },
    OptionSpec{ "IES_UPGRADES_IN_MEMORY", {}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_upgrades_in_memory(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_upgrades_in_memory(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_upgrades_in_memory()?1:0); } },
    OptionSpec{ "IES_ORDERED_BINARY", {}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_ordered_binary(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_ordered_binary(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_ordered_binary()?1:0); } },
    OptionSpec{ "IES_MULTIMODAL_ALPHA", {}, OptType::DOUBLE, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_ies_multimodal_alpha(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_multimodal_alpha(0.0); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_ies_multimodal_alpha(); return ss.str(); } },
    OptionSpec{ "IES_MULTIMODAL_WEIGHT_EXPONENT", {}, OptType::DOUBLE, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_ies_multimodal_weight_exponent(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_multimodal_weight_exponent(0.0); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_ies_multimodal_weight_exponent(); return ss.str(); } },
    OptionSpec{ "IES_MULTIMODAL_PHI_WEIGHT", {}, OptType::DOUBLE, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_ies_multimodal_phi_weight(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_multimodal_phi_weight(0.5); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_ies_multimodal_phi_weight(); return ss.str(); } },
    OptionSpec{ "IES_LOCALIZER_FORGIVE_MISSING", {"IES_LOCALIZER_FORGIVE_EXTRA"}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_localizer_forgive_missing(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_localizer_forgive_missing(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_localizer_forgive_missing()?1:0); } },
    OptionSpec{ "IES_PHI_FACTOR_FILE", {}, OptType::STRING, "ies", false,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_ies_phi_fractions_files(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_phi_fractions_files(""); },
        [](const PestppOptions& o)->string{ return o.get_ies_phi_fractions_file(); } },
    OptionSpec{ "IES_PHI_FACTORS_BY_REAL", {}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_phi_factors_by_real(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_phi_factors_by_real(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_phi_factors_by_real()?1:0); } },
    OptionSpec{ "IES_N_ITER_MEAN", {"IES_N_ITER_REINFLATE"}, OptType::VEC_INT, "ies", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ vector<int> v; vector<string> tok; tokenize(value,tok,","); for(const auto& t:tok){ v.push_back(convert_cp<int>(t)); } o.set_ies_n_iter_reinflate(v); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_n_iter_reinflate(vector<int>{0}); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; auto v=o.get_ies_n_iter_reinflate(); for(size_t i=0;i<v.size();++i){if(i)ss<<",";ss<<v[i];} return ss.str(); } },
    OptionSpec{ "IES_REINFLATE_FACTOR", {}, OptType::VEC_DOUBLE, "ies", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ vector<double> v; vector<string> tok; tokenize(value,tok,","); for(const auto& t:tok){ v.push_back(convert_cp<double>(t)); } o.set_ies_reinflate_factor(v); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_reinflate_factor(vector<double>{1.0}); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; auto v=o.get_ies_reinflate_factor(); for(size_t i=0;i<v.size();++i){if(i)ss<<",";ss<<v[i];} return ss.str(); } },
    OptionSpec{ "IES_UPDATE_BY_REALS", {}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_updatebyreals(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_updatebyreals(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_updatebyreals()?1:0); } },
    OptionSpec{ "SAVE_DENSE", {}, OptType::BOOL, "general", true,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_save_dense(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_save_dense(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_save_dense()?1:0); } },
    OptionSpec{ "IES_AUTOADALOC_INDICATOR_PARS", {}, OptType::VEC_STRING, "ies", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ vector<string> v; vector<string> tok; tokenize(value,tok,","); for(const auto& t:tok){ v.push_back(upper_cp(t)); } o.set_ies_aal_indicator_pars(v); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_aal_indicator_pars({}); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; auto v=o.get_ies_aal_indicator_pars(); for(size_t i=0;i<v.size();++i){if(i)ss<<",";ss<<v[i];} return ss.str(); } },
    OptionSpec{ "IES_RUN_REALNAME", {}, OptType::CUSTOM, "ies", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_ies_run_realname(value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_run_realname(string()); },
        [](const PestppOptions& o)->string{ return o.get_ies_run_realname(); } },
    OptionSpec{ "IES_REINFLATE_NUM_REALS", {}, OptType::VEC_INT, "ies", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ vector<int> v; vector<string> tok; tokenize(value,tok,","); for(const auto& t:tok){ v.push_back(convert_cp<int>(t)); } o.set_ies_reinflate_num_reals(v); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_reinflate_num_reals(vector<int>{0}); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; auto v=o.get_ies_reinflate_num_reals(); for(size_t i=0;i<v.size();++i){if(i)ss<<",";ss<<v[i];} return ss.str(); } },
    OptionSpec{ "IES_USE_PHI_LAMBDA_ITERS", {}, OptType::BOOL, "ies", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_ies_use_phi_lambda_iters(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_ies_use_phi_lambda_iters(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_ies_use_phi_lambda_iters()?1:0); } },
    OptionSpec{ "DA_PARAMETER_CYCLE_TABLE", {}, OptType::STRING, "da", false,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_da_par_cycle_table(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_da_par_cycle_table(""); },
        [](const PestppOptions& o)->string{ return o.get_da_par_cycle_table(); } },
    OptionSpec{ "DA_OBSERVATION_CYCLE_TABLE", {}, OptType::STRING, "da", false,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_da_obs_cycle_table(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_da_obs_cycle_table(""); },
        [](const PestppOptions& o)->string{ return o.get_da_obs_cycle_table(); } },
    OptionSpec{ "DA_WEIGHT_CYCLE_TABLE", {}, OptType::STRING, "da", false,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_da_weight_cycle_table(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_da_weight_cycle_table(string()); },
        [](const PestppOptions& o)->string{ return o.get_da_weight_cycle_table(); } },
    OptionSpec{ "DA_HOTSTART_CYCLE", {}, OptType::INT, "da", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_da_hotstart_cycle(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_da_hotstart_cycle(0); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_da_hotstart_cycle()); } },
    OptionSpec{ "DA_STOP_CYCLE", {}, OptType::INT, "da", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_da_stop_cycle(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_da_stop_cycle(1000000000); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_da_stop_cycle()); } },
    OptionSpec{ "DA_USE_SIMULATED_STATES", {}, OptType::BOOL, "da", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_da_use_simulated_states(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_da_use_simulated_states(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_da_use_simulated_states()?1:0); } },
    OptionSpec{ "DA_NOPTMAX_SCHEDULE", {}, OptType::STRING, "da", false,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_da_noptmax_schedule(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_da_noptmax_schedule(""); },
        [](const PestppOptions& o)->string{ return o.get_da_noptmax_schedule(); } },
    OptionSpec{ "PANTHER_AGENT_RESTART_ON_ERROR", {}, OptType::BOOL, "panther", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_panther_agent_restart_on_error(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_panther_agent_restart_on_error(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_panther_agent_restart_on_error()?1:0); } },
    // Runs on the AGENT, only when the master asks for partial results, and only just before
    // the output files are read. Failures are forgiven by design - see the call site: a
    // partial-results request is a courtesy and must never be able to damage the run it is
    // asking about.
    // Closes FRACPHIM's open loop - see the note in SVDSolver::dynamic_weight_adj(). ON by
    // default: it only ever RELAXES a target, and only where FRACPHIM is the binding term and
    // phi has genuinely stalled, so a run whose phimlim is attainable is provably unaffected
    // (api_glm_dynreg_reachable_phimlim_test). Set it false to get the pre-2026 trajectory.
    OptionSpec{ "REG_USE_ACHIEVABLE_TARGET", {}, OptType::BOOL, "glm", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_reg_use_achievable_target(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_reg_use_achievable_target(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_reg_use_achievable_target()?1:0); } },
    OptionSpec{ "PANTHER_WORKER_PARTIAL_OBS_COMMAND", {}, OptType::STRING, "panther", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_panther_worker_partial_obs_command(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_panther_worker_partial_obs_command(""); },
        [](const PestppOptions& o)->string{ return o.get_panther_worker_partial_obs_command(); } },
    // Agent-only, and NOT init_only: the agent reads it at the moment a partial request is
    // answered, so a live set takes effect on the next request rather than needing a restart.
    OptionSpec{ "PANTHER_WORKER_STATUS_FILE", {}, OptType::STRING, "panther", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_panther_worker_status_file(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_panther_worker_status_file(""); },
        [](const PestppOptions& o)->string{ return o.get_panther_worker_status_file(); } },
    OptionSpec{ "PANTHER_AGENT_NO_PING_TIMEOUT_SECS", {}, OptType::INT, "panther", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_panther_agent_no_ping_timeout_secs(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_panther_agent_no_ping_timeout_secs(-1); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_panther_agent_no_ping_timeout_secs()); } },
    OptionSpec{ "ADDITIONAL_INS_DELIMITERS", {}, OptType::STRING, "general", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_additional_ins_delimiters(value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_additional_ins_delimiters(""); },
        [](const PestppOptions& o)->string{ return o.get_additional_ins_delimiters(); } },
    OptionSpec{ "RANDOM_SEED", {"RAND_SEED"}, OptType::INT, "general", true,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_random_seed(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_random_seed(358183147); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_random_seed()); } },
    OptionSpec{ "GLM_ITER_MC", {}, OptType::BOOL, "glm", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_glm_iter_mc(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_glm_iter_mc(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_glm_iter_mc()?1:0); } },
    OptionSpec{ "PANTHER_DEBUG_LOOP", {}, OptType::BOOL, "panther", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_panther_debug_loop(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_panther_debug_loop(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_panther_debug_loop()?1:0); } },
    OptionSpec{ "DEBUG_CHECK_PAR_EN_CONSISTENCY", {}, OptType::BOOL, "general", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_debug_check_par_en_consistency(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_debug_check_par_en_consistency(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_debug_check_par_en_consistency()?1:0); } },
    OptionSpec{ "PANTHER_AGENT_FREEZE_ON_FAIL", {}, OptType::BOOL, "panther", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_panther_debug_fail_freeze(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_panther_debug_fail_freeze(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_panther_debug_fail_freeze()?1:0); } },
    OptionSpec{ "SAVE_ALL_RUNS", {}, OptType::BOOL, "general", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_save_all_runs(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_save_all_runs(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_save_all_runs()?1:0); } },
    OptionSpec{ "CHECK_TPLINS", {}, OptType::BOOL, "general", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_check_tplins(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_check_tplins(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_check_tplins()?1:0); } },
    OptionSpec{ "FILL_TPL_ZEROS", {}, OptType::BOOL, "general", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_fill_tpl_zeros(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_fill_tpl_zeros(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_fill_tpl_zeros()?1:0); } },
    OptionSpec{ "TPL_FORCE_DECIMAL", {}, OptType::BOOL, "general", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_tpl_force_decimal(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_tpl_force_decimal(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_tpl_force_decimal()?1:0); } },
    OptionSpec{ "FORGIVE_UNKNOWN_ARGS", {}, OptType::BOOL, "general", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_forgive_unknown_args(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_forgive_unknown_args(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_forgive_unknown_args()?1:0); } },
    OptionSpec{ "PANTHER_ECHO", {}, OptType::BOOL, "panther", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_panther_echo(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_panther_echo(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_panther_echo()?1:0); } },
    OptionSpec{ "NUM_TPL_INS_THREADS", {}, OptType::INT, "general", true,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_num_tpl_ins_threads(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_num_tpl_ins_threads(1); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_num_tpl_ins_threads()); } },
    OptionSpec{ "PANTHER_TRANSFER_ON_FINISH", {}, OptType::VEC_STRING, "panther", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ vector<string> v; vector<string> tok; tokenize(org_value,tok,","); for(const auto& t:tok){ v.push_back(strip_cp(t)); } o.set_panther_transfer_on_finish(v); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_panther_transfer_on_finish({}); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; auto v=o.get_panther_transfer_on_finish(); for(size_t i=0;i<v.size();++i){if(i)ss<<",";ss<<v[i];} return ss.str(); } },
    OptionSpec{ "PANTHER_TRANSFER_ON_FAIL", {}, OptType::VEC_STRING, "panther", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ vector<string> v; vector<string> tok; tokenize(org_value,tok,","); for(const auto& t:tok){ v.push_back(strip_cp(t)); } o.set_panther_transfer_on_fail(v); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_panther_transfer_on_fail({}); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; auto v=o.get_panther_transfer_on_fail(); for(size_t i=0;i<v.size();++i){if(i)ss<<",";ss<<v[i];} return ss.str(); } },
    OptionSpec{ "PANTHER_MASTER_TIMEOUT_MILLISECONDS", {}, OptType::INT, "panther", true,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_panther_timeout_milliseconds(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_panther_timeout_milliseconds(-999); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_panther_timeout_milliseconds()); } },
    OptionSpec{ "PANTHER_MASTER_ECHO_INTERVAL_MILLISECONDS", {}, OptType::INT, "panther", true,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_panther_echo_interval_milliseconds(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_panther_echo_interval_milliseconds(1000); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_panther_echo_interval_milliseconds()); } },
    OptionSpec{ "PANTHER_PERSISTENT_WORKERS", {}, OptType::BOOL, "panther", true,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_panther_persistent_workers(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_panther_persistent_workers(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_panther_persistent_workers()?1:0); } },
    OptionSpec{ "PANTHER_PING_INTERVAL_SECS", {}, OptType::INT, "panther", true,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_panther_ping_interval_secs(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_panther_ping_interval_secs(60); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_panther_ping_interval_secs()); } },
    OptionSpec{ "MOU_GENERATOR", {}, OptType::STRING, "mou", false,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_mou_generator(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_generator("PSO"); },
        [](const PestppOptions& o)->string{ return o.get_mou_generator(); } },
    OptionSpec{ "MOU_POPULATION_SIZE", {}, OptType::INT, "mou", true,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_mou_population_size(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_population_size(100); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_mou_population_size()); } },
    OptionSpec{ "MOU_DV_POPULATION_FILE", {}, OptType::STRING, "mou", true,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_mou_dv_population_file(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_dv_population_file(""); },
        [](const PestppOptions& o)->string{ return o.get_mou_dv_population_file(); } },
    OptionSpec{ "MOU_OBS_POPULATION_RESTART_FILE", {}, OptType::STRING, "mou", true,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_mou_obs_population_restart_file(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_obs_population_restart_file(""); },
        [](const PestppOptions& o)->string{ return o.get_mou_obs_population_restart_file(); } },
    OptionSpec{ "MOU_OBJECTIVES", {}, OptType::VEC_STRING, "mou", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ vector<string> v; vector<string> tok; tokenize(value,tok,",\t\t"); for(const auto& t:tok){ v.push_back(strip_cp(t)); } o.set_mou_objectives(v); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_objectives({}); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; auto v=o.get_mou_objectives(); for(size_t i=0;i<v.size();++i){if(i)ss<<",";ss<<v[i];} return ss.str(); } },
    OptionSpec{ "MOU_MAX_ARCHIVE_SIZE", {}, OptType::INT, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_mou_max_archive_size(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_max_archive_size(500); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_mou_max_archive_size()); } },
    OptionSpec{ "OPT_CHANCE_POINTS", {}, OptType::CUSTOM, "opt", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_opt_chance_points(value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_opt_chance_points("SINGLE"); },
        [](const PestppOptions& o)->string{ return o.get_opt_chance_points(); } },
    OptionSpec{ "MOU_RISK_OBJECTIVE", {}, OptType::BOOL, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_mou_risk_obj(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_risk_obj(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_mou_risk_obj()?1:0); } },
    OptionSpec{ "MOU_VERBOSE_LEVEL", {}, OptType::INT, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_mou_verbose_level(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_verbose_level(1); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_mou_verbose_level()); } },
    OptionSpec{ "MOU_ENV_SELECTOR", {}, OptType::CUSTOM, "mou", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_mou_env_selector(upper_cp(strip_cp(value))); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_env_selector("NSGA"); },
        [](const PestppOptions& o)->string{ return o.get_mou_env_selector(); } },
    OptionSpec{ "MOU_MATING_SELECTOR", {}, OptType::CUSTOM, "mou", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_mou_mating_selector(upper_cp(strip_cp(value))); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_mating_selector("TOURNAMENT"); },
        [](const PestppOptions& o)->string{ return o.get_mou_mating_selector(); } },
    OptionSpec{ "MOU_CROSSOVER_PROBABILITY", {}, OptType::DOUBLE, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_mou_crossover_probability(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_crossover_probability(0.75); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_mou_crossover_probability(); return ss.str(); } },
    OptionSpec{ "MOU_MUTATION_PROBABILITY", {}, OptType::DOUBLE, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_mou_mutation_probability(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_mutation_probability(-999); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_mou_mutation_probability(); return ss.str(); } },
    OptionSpec{ "MOU_DE_F", {}, OptType::DOUBLE, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_mou_de_f(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_de_f(0.8); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_mou_de_f(); return ss.str(); } },
    OptionSpec{ "MOU_SAVE_POPULATION_EVERY", {}, OptType::INT, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_mou_save_population_every(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_save_population_every(-1); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_mou_save_population_every()); } },
    OptionSpec{ "MOU_PSO_OMEGA", {}, OptType::DOUBLE, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_mou_pso_omega(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_pso_omega(0.7); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_mou_pso_omega(); return ss.str(); } },
    OptionSpec{ "MOU_PSO_SOCIAL_CONST", {}, OptType::VEC_DOUBLE, "mou", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ vector<double> v; vector<string> tok; tokenize(value,tok,",\t\t"); for(const auto& t:tok){ v.push_back(convert_cp<double>(t)); } o.set_mou_pso_social_const(v); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_pso_social_const(vector<double>{2.0}); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; auto v=o.get_mou_pso_social_const(); for(size_t i=0;i<v.size();++i){if(i)ss<<",";ss<<v[i];} return ss.str(); } },
    OptionSpec{ "MOU_PSO_COGNITIVE_CONST", {}, OptType::VEC_DOUBLE, "mou", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ vector<double> v; vector<string> tok; tokenize(value,tok,",\t\t"); for(const auto& t:tok){ v.push_back(convert_cp<double>(t)); } o.set_mou_pso_cognitive_const(v); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_pso_cognitive_const(vector<double>{2.0}); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; auto v=o.get_mou_pso_cognitive_const(); for(size_t i=0;i<v.size();++i){if(i)ss<<",";ss<<v[i];} return ss.str(); } },
    OptionSpec{ "MOU_PSO_ALPHA", {}, OptType::DOUBLE, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_mou_pso_alpha(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_pso_alpha(1.0); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_mou_pso_alpha(); return ss.str(); } },
    OptionSpec{ "MOU_PSO_RRAMP", {}, OptType::DOUBLE, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_mou_pso_rramp(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_pso_rramp(-5e+02); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_mou_pso_rramp(); return ss.str(); } },
    OptionSpec{ "MOU_PSO_RFIT", {}, OptType::DOUBLE, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_mou_pso_rfit(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_pso_rfit(2.0); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_mou_pso_rfit(); return ss.str(); } },
    OptionSpec{ "MOU_PSO_INERTIA", {}, OptType::VEC_DOUBLE, "mou", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ vector<double> v; vector<string> tok; tokenize(value,tok,",\t\t"); for(const auto& t:tok){ v.push_back(convert_cp<double>(t)); } o.set_mou_pso_inertia(v); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_pso_inertia(vector<double>{0.7, 0.4, 0}); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; auto v=o.get_mou_pso_inertia(); for(size_t i=0;i<v.size();++i){if(i)ss<<",";ss<<v[i];} return ss.str(); } },
    OptionSpec{ "MOU_PSO_VMAX_FACTOR", {}, OptType::DOUBLE, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_mou_pso_vmax_factor(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_pso_vmax_factor(0.8); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_mou_pso_vmax_factor(); return ss.str(); } },
    OptionSpec{ "MOU_PSO_DV_BOUND_HANDLING", {}, OptType::CUSTOM, "mou", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_mou_pso_dv_bound_handling(upper_cp(strip_cp(value))); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_pso_dv_bound_handling("CLAMP"); },
        [](const PestppOptions& o)->string{ return o.get_mou_pso_dv_bound_handling(); } },
    OptionSpec{ "MOU_OUTER_REPO_OBS_FILE", {}, OptType::STRING, "mou", false,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_mou_outer_repo_obs_file(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_outer_repo_obs_file(""); },
        [](const PestppOptions& o)->string{ return o.get_mou_outer_repo_obs_file(); } },
    OptionSpec{ "MOU_MAX_NN_SEARCH", {}, OptType::INT, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_mou_max_nn_search(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_max_nn_search(o.get_mou_population_size()); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_mou_max_nn_search()); } },
    OptionSpec{ "MOU_HYPERVOLUME_EXTREME", {}, OptType::DOUBLE, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_mou_hypervolume_extreme(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_hypervolume_extreme(1e+10); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_mou_hypervolume_extreme(); return ss.str(); } },
    OptionSpec{ "MOU_INFILL_SIZE", {}, OptType::INT, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_mou_infill_size(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_infill_size(100); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_mou_infill_size()); } },
    OptionSpec{ "MOU_RESAMPLE_EVERY", {}, OptType::INT, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_mou_resample_every(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_resample_every(-1); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_mou_resample_every()); } },
    OptionSpec{ "MOU_RESAMPLE_COMMAND", {}, OptType::CUSTOM, "mou", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_mou_resample_command(value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_resample_command(string()); },
        [](const PestppOptions& o)->string{ return o.get_mou_resample_command(); } },
    OptionSpec{ "MOU_PPD_BETA", {}, OptType::DOUBLE, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_mou_ppd_beta(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_ppd_beta(0.5); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_mou_ppd_beta(); return ss.str(); } },
    OptionSpec{ "MOU_FIT_GAMMA", {}, OptType::DOUBLE, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_mou_fit_gamma(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_fit_gamma(0.25); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_mou_fit_gamma(); return ss.str(); } },
    OptionSpec{ "MOU_FIT_EPSILON", {}, OptType::DOUBLE, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_mou_fit_epsilon(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_fit_epsilon(0.05); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_mou_fit_epsilon(); return ss.str(); } },
    OptionSpec{ "MOU_POPULATION_SCHEDULE", {}, OptType::STRING, "mou", false,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_mou_population_schedule(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_population_schedule(""); },
        [](const PestppOptions& o)->string{ return o.get_mou_population_schedule(); } },
    OptionSpec{ "OPT_CHANCE_SCHEDULE", {}, OptType::STRING, "opt", false,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_opt_chance_schedule(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_opt_chance_schedule(""); },
        [](const PestppOptions& o)->string{ return o.get_opt_chance_schedule(); } },
    OptionSpec{ "MOU_SIMPLEX_REFLECTIONS", {}, OptType::INT, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_mou_simplex_reflections(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_simplex_reflections(10); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_mou_simplex_reflections()); } },
    OptionSpec{ "MOU_SIMPLEX_FACTORS", {}, OptType::VEC_DOUBLE, "mou", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ vector<double> v; vector<string> tok; tokenize(value,tok,",\t\t"); for(const auto& t:tok){ v.push_back(convert_cp<double>(t)); } o.set_mou_simplex_factors(v); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_simplex_factors(vector<double>{0.5, 0.6, 0.7, 0.8}); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; auto v=o.get_mou_simplex_factors(); for(size_t i=0;i<v.size();++i){if(i)ss<<",";ss<<v[i];} return ss.str(); } },
    OptionSpec{ "MOU_SIMPLEX_MUTATION", {}, OptType::BOOL, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_mou_simplex_mutation(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_simplex_mutation(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_mou_simplex_mutation()?1:0); } },
    OptionSpec{ "MOU_USE_MULTIGEN_POPULATION", {}, OptType::BOOL, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_mou_use_multigen(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_use_multigen(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_mou_use_multigen()?1:0); } },
    OptionSpec{ "MOU_SHUFFLE_FIXED_PARS", {}, OptType::BOOL, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_mou_shuffle_fixed_pars(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_shuffle_fixed_pars(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_mou_shuffle_fixed_pars()?1:0); } },
    OptionSpec{ "MOU_DEBUG_DV_HANDLING", {}, OptType::BOOL, "mou", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_mou_debug_dv_handling(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_mou_debug_dv_handling(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_mou_debug_dv_handling()?1:0); } },
    OptionSpec{ "SQP_DV_EN", {}, OptType::STRING, "sqp", true,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_sqp_dv_en(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_dv_en(""); },
        [](const PestppOptions& o)->string{ return o.get_sqp_dv_en(); } },
    OptionSpec{ "SQP_RESTART_OBS_EN", {}, OptType::STRING, "sqp", true,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_sqp_obs_restart_en(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_obs_restart_en(""); },
        [](const PestppOptions& o)->string{ return o.get_sqp_obs_restart_en(); } },
    OptionSpec{ "SQP_SEARCH_METHOD", {}, OptType::STRING, "sqp", false,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_sqp_search_method(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_search_method("LINE"); },
        [](const PestppOptions& o)->string{ return o.get_sqp_search_method(); } },
    OptionSpec{ "SQP_SOLVE_METHOD", {}, OptType::STRING, "sqp", false,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_sqp_solve_method(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_solve_method("NULL"); },
        [](const PestppOptions& o)->string{ return o.get_sqp_solve_method(); } },
    OptionSpec{ "SQP_CMA_BOUND_HANDLING", {}, OptType::STRING, "sqp", false,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_sqp_cma_bound_handling(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_cma_bound_handling("reject"); },
        [](const PestppOptions& o)->string{ return o.get_sqp_cma_bound_handling(); } },
    OptionSpec{ "SQP_NUM_REALS", {}, OptType::INT, "sqp", true,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_sqp_num_reals(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_num_reals(50); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_sqp_num_reals()); } },
    OptionSpec{ "SQP_NUM_REFINED_SEARCH_PTS", {}, OptType::INT, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_sqp_num_refined_search_pts(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_num_refined_search_pts(1.0); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_sqp_num_refined_search_pts()); } },
    OptionSpec{ "SQP_SUBSET_SIZE", {}, OptType::INT, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_sqp_subset_size(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_subset_size(-10); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_sqp_subset_size()); } },
    OptionSpec{ "SQP_UPDATE_HESSIAN", {}, OptType::BOOL, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_sqp_update_hessian(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_update_hessian(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_sqp_update_hessian()?1:0); } },
    OptionSpec{ "SQP_HESSIAN_UPDATE_METHOD", {}, OptType::STRING, "sqp", false,
        [](PestppOptions& o,const string&,const string& org_value)->PestppOptions::ARG_STATUS{ o.set_sqp_hessian_update_method(org_value); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_hessian_update_method("BFGS"); },
        [](const PestppOptions& o)->string{ return o.get_sqp_hessian_update_method(); } },
    OptionSpec{ "SQP_FILTER_TOL", {}, OptType::DOUBLE, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_sqp_filter_tol(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_filter_tol(0.001); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_sqp_filter_tol(); return ss.str(); } },
    // init-only: SeqQuadProgram seeds its adaptive working_set_tol from this once and then
    // tightens it per accepted iteration, so a later change would be silently ignored
    OptionSpec{ "SQP_WORKING_SET_TOL", {}, OptType::DOUBLE, "sqp", true,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_sqp_working_set_tol(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_working_set_tol(0.10); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_sqp_working_set_tol(); return ss.str(); } },
    OptionSpec{ "SQP_CMA_REINFLATION_FACTOR", {}, OptType::DOUBLE, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_sqp_cma_reinflation_factor(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_cma_reinflation_factor(-1.0); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_sqp_cma_reinflation_factor(); return ss.str(); } },
    OptionSpec{ "SQP_ALPHA_MULTS", {}, OptType::VEC_DOUBLE, "sqp", false,
        [](PestppOptions& o,const string& value,const string& org_value)->PestppOptions::ARG_STATUS{ vector<double> v; vector<string> tok; tokenize(value,tok,",\t\t"); for(const auto& t:tok){ v.push_back(convert_cp<double>(t)); } o.set_sqp_alpha_mults(v); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_alpha_mults(vector<double>{0.001, 0.005, 0.01, 0.1, 0.5, 1.0}); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; auto v=o.get_sqp_alpha_mults(); for(size_t i=0;i<v.size();++i){if(i)ss<<",";ss<<v[i];} return ss.str(); } },
    OptionSpec{ "SQP_CMA_C1", {}, OptType::DOUBLE, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_sqp_cma_c1(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_cma_c1(-1); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_sqp_cma_c1(); return ss.str(); } },
    OptionSpec{ "SQP_CMA_CMU", {}, OptType::DOUBLE, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_sqp_cma_cmu(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_cma_cmu(-1); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_sqp_cma_cmu(); return ss.str(); } },
    OptionSpec{ "SQP_CMA_CC", {}, OptType::DOUBLE, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_sqp_cma_cc(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_cma_cc(-1); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_sqp_cma_cc(); return ss.str(); } },
    OptionSpec{ "SQP_CMA_PARENT_NUM", {}, OptType::INT, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_sqp_cma_parent_num(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_cma_parent_num(-1); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_sqp_cma_parent_num()); } },
    OptionSpec{ "SQP_CMA_STEPSIZE_CONTROL", {}, OptType::BOOL, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_sqp_cma_stepsize_control(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_cma_stepsize_control(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_sqp_cma_stepsize_control()?1:0); } },
    OptionSpec{ "SQP_MAX_CONSEC_INFEAS_IES", {}, OptType::INT, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_sqp_max_consec_infeas_ies(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_max_consec_infeas_ies(3); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_sqp_max_consec_infeas_ies()); } },
    OptionSpec{ "SQP_MAX_REINFLATION_COND_NUM", {}, OptType::DOUBLE, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_sqp_max_reinflation_cond_num(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_max_reinflation_cond_num(500.0); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_sqp_max_reinflation_cond_num(); return ss.str(); } },
    OptionSpec{ "SQP_SCALE_DOWN_FACTOR", {}, OptType::DOUBLE, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_sqp_scale_down_factor(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_scale_down_factor(1.0); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_sqp_scale_down_factor(); return ss.str(); } },
    OptionSpec{ "SQP_HESS_MAX_COND_NUM", {}, OptType::DOUBLE, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_sqp_hess_max_cond_num(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_hess_max_cond_num(1E+8); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_sqp_hess_max_cond_num(); return ss.str(); } },
    OptionSpec{ "SQP_SAVE_COV_EVERY", {}, OptType::INT, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_sqp_save_cov_every(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_save_cov_every(-1); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_sqp_save_cov_every()); } },
    OptionSpec{ "SQP_ENFORCE_BOUNDS", {}, OptType::BOOL, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_sqp_enforce_bounds(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_enforce_bounds(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_sqp_enforce_bounds()?1:0); } },
    OptionSpec{ "SQP_VIOL_PAD", {}, OptType::DOUBLE, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_sqp_viol_pad(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_viol_pad(1E-4); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_sqp_viol_pad(); return ss.str(); } },
    OptionSpec{ "SQP_RESET_HESSIAN_EVERY", {}, OptType::INT, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_sqp_reset_hessian_every(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_reset_hessian_every(-1); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_sqp_reset_hessian_every()); } },
    OptionSpec{ "SQP_USE_ENSEMBLE_APPROX_HESSIAN", {}, OptType::BOOL, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_sqp_use_ensemble_approx_hessian(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_use_ensemble_approx_hessian(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_sqp_use_ensemble_approx_hessian()?1:0); } },
    OptionSpec{ "SQP_RESCALE_SEARCH_DIR", {}, OptType::BOOL, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_sqp_rescale_search_dir(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_rescale_search_dir(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_sqp_rescale_search_dir()?1:0); } },
    OptionSpec{ "SQP_SEEK_FEAS_MAX_ITER", {}, OptType::INT, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ int x; convert_ip(value,x); o.set_sqp_seek_feas_max_iter(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_seek_feas_max_iter(3); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_sqp_seek_feas_max_iter()); } },
    OptionSpec{ "SQP_RISK", {}, OptType::DOUBLE, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ (void)o; (void)value; throw runtime_error("++sqp_risk has been retired. For chance-constrained sqp use ++opt_risk (the same option mou and pestpp-opt use). For robust optimization over paired decision-variable/parameter realizations use ++opt_use_robust(true), which does no risk shifting and cannot be combined with opt_risk."); },
        [](PestppOptions& o){ o.set_sqp_risk(0.50); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_sqp_risk(); return ss.str(); } },
    OptionSpec{ "SQP_POWELL_DAMPING_FACTOR", {}, OptType::DOUBLE, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ double x; convert_ip(value,x); o.set_sqp_powell_damping_factor(x); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_powell_damping_factor(0.2); },
        [](const PestppOptions& o)->string{ std::ostringstream ss; ss<<o.get_sqp_powell_damping_factor(); return ss.str(); } },
    OptionSpec{ "SQP_DEBUG_ENABLE_CONSTRAINT_WEIGHTED_JCO", {}, OptType::BOOL, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_sqp_debug_enable_constraint_weighted_jco(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_debug_enable_constraint_weighted_jco(true); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_sqp_debug_enable_constraint_weighted_jco()?1:0); } },
    OptionSpec{ "SQP_DEBUG_HESSIAN", {}, OptType::BOOL, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_sqp_debug_hessian(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_debug_hessian(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_sqp_debug_hessian()?1:0); } },
    OptionSpec{ "SQP_DEBUG_CMA", {}, OptType::BOOL, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_sqp_debug_cma(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_debug_cma(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_sqp_debug_cma()?1:0); } },
    OptionSpec{ "SQP_DEBUG_STOSAG_GRAD", {}, OptType::BOOL, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_sqp_debug_stosag_grad(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_debug_stosag_grad(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_sqp_debug_stosag_grad()?1:0); } },
    OptionSpec{ "SQP_USE_IES_INFEAS", {}, OptType::BOOL, "sqp", false,
        [](PestppOptions& o,const string& value,const string&)->PestppOptions::ARG_STATUS{ o.set_sqp_use_ies_infeas(pest_utils::parse_string_arg_to_bool(value)); return PestppOptions::ARG_STATUS::ARG_ACCEPTED; },
        [](PestppOptions& o){ o.set_sqp_use_ies_infeas(false); },
        [](const PestppOptions& o)->string{ return std::to_string(o.get_sqp_use_ies_infeas()?1:0); } }
    };
    return R;
}

static const map<string, const OptionSpec*>& registry_index()
{
    static map<string, const OptionSpec*> idx;
    if (idx.empty())
    {
        for (const auto& s : PestppOptions::get_option_registry())
        {
            idx[s.name] = &s;
            for (const auto& a : s.aliases) idx[a] = &s;
        }
    }
    return idx;
}

// Resolve a user-supplied option name to its registry entry, honouring the DA_* -> IES_*
// rewrite that assign_registry_impl() performs below.
//
// Without this the same name means different things depending on direction:
// set_option("da_num_reals") is accepted, because the assign path rewrites it, while
// get_option("da_num_reals") returns "" and is_valid_arg() says it does not exist, because
// they looked the name up raw. A caller that sets an option and reads it back - which is
// exactly what a programmatic api does - got silence.
static const OptionSpec* lookup_spec(const string& name)
{
    const auto& idx = registry_index();
    string k = upper_cp(name);
    auto it = idx.find(k);
    if (it != idx.end())
        return it->second;
    // same test the assign path uses, deliberately, so the two cannot drift apart
    if (k.substr(0, 2) == "DA")
    {
        auto it2 = idx.find("IES" + k.substr(2, k.size()));
        if (it2 != idx.end())
            return it2->second;
    }
    return nullptr;
}

void PestppOptions::set_defaults_registry()
{
    for (const auto& s : get_option_registry()) s.apply_default(*this);
}

void PestppOptions::summary_registry(ostream& os) const
{
    os << endl << "    PEST++ OPTIONS (registry): " << endl << endl;
    for (const auto& s : get_option_registry())
        os << s.name << ": " << s.to_str(*this) << endl;
}

bool PestppOptions::is_valid_arg(const string& key) const
{
    return lookup_spec(key) != nullptr;
}

set<string> PestppOptions::get_registered_args() const
{
    set<string> out;
    for (const auto& s : get_option_registry()) out.insert(s.name);
    return out;
}

// Shared core for both the control-file parse path and the programmatic set_option() path.
// check_duplicate=true is the file-parse behavior (a repeated ++arg is an error); the
// programmatic path passes false so a caller can change an option repeatedly (e.g. between
// algorithm iterations). Both record provenance in passed_args so tool-default overrides
// treat the option as user-supplied.
// Options that have been REMOVED, and what to do instead. Checked before the registry lookup
// and deliberately OUTSIDE the parse try/catch below, which turns every exception into a bare
// ARG_INVALID ("invalid value for option 'X'") - useless for a retired option, where the value
// was never the problem and the caller needs to be told what replaced it.
string PestppOptions::get_retired_message(const string& name)
{
    //case-normalized here rather than at the call sites: the '++' parser leaves keys exactly as
    //written (its upper_ip is commented out), while the registry uppercases internally - so the
    //two callers hand this different spellings of the same option.
    string key = upper_cp(name);
    //the pestpp-glm differential evolution implementation was removed - pestpp-mou has carried a
    //maintained DE generator for a long time - and these options went with it.  They are refused
    //by NAME here rather than left to be accepted and ignored, which is what ++global_opt(MOEA),
    //++moea_name and the whole ++de_* family were doing: parsed, stored, read by nothing.
    if (key == "GLOBAL_OPT")
        return "++global_opt has been retired. Its only values were DE - whose implementation has been removed from pestpp-glm - and MOEA, which never did anything. Use pestpp-mou, which is the maintained global/multi-objective optimizer.";
    if (key == "MOEA_NAME")
        return "++moea_name has been retired; it was read by nothing. Use pestpp-mou and select the generator with ++mou_generator (de, sbx, pm, pso or smp).";
    if (key == "DE_F")
        return "++de_f has been retired along with the pestpp-glm differential evolution implementation. Use pestpp-mou: ++mou_generator(de) and ++mou_de_f.";
    if (key == "DE_CR")
        return "++de_cr has been retired along with the pestpp-glm differential evolution implementation. Use pestpp-mou: ++mou_generator(de) and ++mou_crossover_probability.";
    if (key == "DE_DITHER_F")
        return "++de_dither_f has been retired along with the pestpp-glm differential evolution implementation. Use pestpp-mou: ++mou_generator(de) and ++mou_de_f.";
    if (key == "DE_POP_SIZE")
        return "++de_pop_size has been retired along with the pestpp-glm differential evolution implementation. Use pestpp-mou: ++mou_population_size.";
    if (key == "DE_MAX_GEN")
        return "++de_max_gen has been retired along with the pestpp-glm differential evolution implementation. Use pestpp-mou and set noptmax to the number of generations.";
    if (key == "SQP_RISK")
        return "++sqp_risk has been retired. It was both a risk value AND the switch that "
               "selected ensemble-based shifting, so it silently overrode opt_risk. Use "
               "++opt_risk for chance-constrained sqp (the same option mou and pestpp-opt "
               "use), or ++opt_use_robust(true) for robust optimization over paired "
               "decision-variable/parameter realizations, which does no risk shifting and "
               "cannot be combined with opt_risk.";
    return string();
}

PestppOptions::ARG_STATUS PestppOptions::assign_registry_impl(string key, const string org_value, bool check_duplicate)
{
    upper_ip(key);
    {
        string retired = get_retired_message(key);
        if (retired.size() > 0)
            throw runtime_error(retired);
    }
    string value = upper_cp(org_value);
    if (value.size() == 0) return ARG_STATUS::ARG_INVALID;
    if (check_duplicate && passed_args.find(key) != passed_args.end()) return ARG_STATUS::ARG_DUPLICATE;
    arg_map[key] = value;
    const auto& idx = registry_index();
    auto it = idx.find(key);
    if (it != idx.end())
    {
        ARG_STATUS st;
        try { st = it->second->parse(*this, value, org_value); }
        catch (...) { return ARG_STATUS::ARG_INVALID; }
        if (st == ARG_STATUS::ARG_ACCEPTED)
        {
            passed_args.insert(it->second->name);
            for (const auto& a : it->second->aliases) passed_args.insert(a);
            if (it->second->init_only && options_initialized)
                init_only_change_warnings.push_back(
                    "option '" + it->second->name + "' is init-only and was changed after "
                    "initialization; the change will not affect the current run");
        }
        return st;
    }
    // DA_* -> IES_* rewrite (mirror of legacy assign_da_value_by_key "hackery")
    if (key.substr(0, 2) == "DA")
    {
        string ies_key = "IES" + key.substr(2, key.size());
        auto it2 = idx.find(ies_key);
        if (it2 != idx.end())
        {
            bool already_found = false;
            if (passed_args.find(ies_key) != passed_args.end())
            {
                passed_args.erase(ies_key);
                already_found = true;
            }
            ARG_STATUS st;
            try { st = it2->second->parse(*this, value, org_value); }
            catch (...) { return ARG_STATUS::ARG_INVALID; }
            bool ok = (st == ARG_STATUS::ARG_ACCEPTED);
            if (ok) passed_da_ies_args.emplace(ies_key);
            if (already_found) passed_args.emplace(ies_key);
            passed_args.insert(key);
            return ok ? ARG_STATUS::ARG_ACCEPTED : ARG_STATUS::ARG_NOTFOUND;
        }
    }
    passed_args.insert(key);
    return ARG_STATUS::ARG_NOTFOUND;
}

PestppOptions::ARG_STATUS PestppOptions::assign_value_by_key_registry(string key, const string org_value)
{
    return assign_registry_impl(key, org_value, true);   // file parse: duplicate is an error
}

// ---- generic programmatic option access ----

// String-keyed programmatic setter: the parse path without the duplicate guard, so a caller
// may change an option repeatedly at runtime. org_value is preserved for case-sensitive
// (filename) options exactly as the file-parse path does.
PestppOptions::ARG_STATUS PestppOptions::set_option(const string& name, const string& org_value)
{
    return assign_registry_impl(name, org_value, false);
}

// Current value of an option rendered as a string (via the registry to_str). Empty string
// if the name is not a registered option.
string PestppOptions::get_option(const string& name) const
{
    const OptionSpec* spec = lookup_spec(name);
    return (spec == nullptr) ? string() : spec->to_str(*this);
}

// Was this option explicitly supplied (by the user in the .pst OR via set_option)? Resolves
// aliases so any spelling of an aliased option answers consistently.
bool PestppOptions::is_user_set(const string& name) const
{
    string key = upper_cp(name);
    const OptionSpec* spec = lookup_spec(name);
    if (spec == nullptr) return passed_args.count(key) > 0;
    // the da spelling is recorded under its own name by the assign path, so check it too
    if (passed_args.count(key)) return true;
    if (passed_args.count(spec->name)) return true;
    for (const auto& a : spec->aliases) if (passed_args.count(a)) return true;
    return false;
}

// Is this option consumed once at setup (irreversible after initialization)?
bool PestppOptions::is_init_only(const string& name) const
{
    const OptionSpec* spec = lookup_spec(name);
    return (spec != nullptr) && spec->init_only;
}

// Which tool the option belongs to (ies, mou, sqp, opt, glm, da, gsa, panther, sweep,
// general); useful for grouping/introspection.
string PestppOptions::get_option_scope(const string& name) const
{
    const OptionSpec* spec = lookup_spec(name);
    return (spec == nullptr) ? string() : spec->scope;
}

// ---- public entry points: parsing + defaults route through the registry ----
void PestppOptions::set_defaults() { set_defaults_registry(); }
PestppOptions::ARG_STATUS PestppOptions::assign_value_by_key(string key, const string org_value)
{
    return assign_value_by_key_registry(key, org_value);
}
void PestppOptions::summary(ostream& os) const { summary_registry(os); }

// ---- self-verification: registry vs the *_legacy reference, byte-for-byte via to_str dumps ----
static string registry_dump(const PestppOptions& o)
{
    ostringstream ss;
    for (const auto& s : PestppOptions::get_option_registry())
        ss << s.name << "=" << s.to_str(o) << "\n";
    return ss.str();
}

bool PestppOptions::self_verify(ostream& os)
{
    int fails = 0;
    // 1) defaults: registry vs legacy
    {
        PestppOptions a; a.set_defaults_legacy();       // legacy impl
        PestppOptions b; b.set_defaults_registry();      // registry
        // compare per option so mismatches are named
        for (const auto& s : get_option_registry())
        {
            string va = s.to_str(a), vb = s.to_str(b);
            if (va != vb)
            {
                os << "DEFAULT MISMATCH " << s.name << ": legacy='" << va << "' registry='" << vb << "'" << endl;
                fails++;
            }
        }
    }
    // 2) parse: registry vs legacy for a per-type probe on every option
    map<OptType,string> probe = {
        {OptType::INT,"5"}, {OptType::DOUBLE,"2.5"}, {OptType::BOOL,"true"},
        {OptType::STRING,"TestFile.csv"}, {OptType::VEC_DOUBLE,"1.0,2.0"},
        {OptType::VEC_INT,"1,2"}, {OptType::VEC_STRING,"a,b"} };
    map<string,string> custom_probe = {
        {"SVD_PACK","REDSVD"},{"GLM_NORMAL_FORM","PRIOR"},{"GLOBAL_OPT","DE"},
        {"OPT_DIRECTION","MAX"},{"OPT_OBJ_FUNC","obj"},{"OPT_CHANCE_POINTS","ALL"},
        {"IES_RUN_REALNAME","base"},{"MOU_ENV_SELECTOR","nsga"},{"MOU_MATING_SELECTOR","tournament"},
        {"MOU_PSO_DV_BOUND_HANDLING","reflect"},{"MOU_RESAMPLE_COMMAND","cmd"},{"MOEA_NAME","de"},
        {"UPGRADE_AUGMENT","1"},{"UPGRADE_BOUNDS","1"},{"AUTO_NORM","1"},{"MAT_INV","1"} };
    for (const auto& s : get_option_registry())
    {
        string p = (s.type==OptType::CUSTOM) ? custom_probe[s.name] : probe[s.type];
        if (p.empty()) continue;
        PestppOptions a; a.set_defaults();
        PestppOptions b; b.set_defaults_registry();
        ARG_STATUS sa, sb;
        try { sa = a.assign_value_by_key_legacy(s.name, p); } catch (...) { sa = ARG_STATUS::ARG_INVALID; }
        try { sb = b.assign_value_by_key_registry(s.name, p); } catch (...) { sb = ARG_STATUS::ARG_INVALID; }
        if (sa != sb)
        {
            os << "PARSE STATUS MISMATCH " << s.name << ": legacy=" << (int)sa << " registry=" << (int)sb << endl;
            fails++; continue;
        }
        string da = s.to_str(a), db = s.to_str(b);
        if (da != db)
        {
            os << "PARSE VALUE MISMATCH " << s.name << " probe='" << p << "': legacy='" << da << "' registry='" << db << "'" << endl;
            fails++;
        }
    }
    os << "self_verify: " << (fails==0 ? "PASS" : "FAIL") << " (" << get_option_registry().size()
       << " options, " << fails << " mismatches)" << endl;
    return fails == 0;
}

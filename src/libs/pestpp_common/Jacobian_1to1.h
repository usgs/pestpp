#ifndef JACOBIAN_1TO1H_
#define JACOBIAN_1TO1H_
#include<map>
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include "Transformable.h"
#include "Jacobian.h"
#include "OutputFileWriter.h"

class ParamTransformSeq;
class ParameterInfo;
class ParameterGroupInfo;
class RunManagerAbstract;
class ModelRun;
class FileManager;
class PriorInformation;
class ParameterRec;

class Jacobian_1to1 : public Jacobian{

public:
	Jacobian_1to1(FileManager &_file_manager, OutputFileWriter &_output_file_writer);
	virtual bool build_runs(ModelRun &model_run, vector<string> numeric_par_names, ParamTransformSeq &par_transform,
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info,
		RunManagerAbstract &run_manager, set<string> &out_of_bound_par, bool phiredswh_flag=false, bool calc_init_obs=true);
	virtual bool build_runs(Parameters &ctl_pars, Observations &ctl_obs, vector<string> numeric_par_names, ParamTransformSeq &par_transform,
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info,
		RunManagerAbstract &run_manager, set<string> &out_of_bound_par, bool phiredswh_flag = false, bool calc_init_obs = true);
	virtual void make_runs(RunManagerAbstract &run_manager);
	virtual bool process_runs(ParamTransformSeq &par_transform,
		const ParameterGroupInfo &group_info,
		RunManagerAbstract &run_manager, const PriorInformation &prior_info, bool splitswh_flag);
	virtual void report_errors(std::ostream &fout);
	virtual ~Jacobian_1to1();
protected:
	Parameters failed_ctl_parameters;
	Parameters failed_to_increment_parmaeters;
	OutputFileWriter* output_file_writer_ptr;
	bool forward_diff(const string &par_name, double derivative_par_value,
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, const ParamTransformSeq &par_trans, double &new_par_val);
	bool central_diff(const string &par_name, double derivative_par_value,
		const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info, const ParamTransformSeq &par_trans, vector<double> &new_par_vec,
		vector<Parameters>  &numeric_dir_par_vec);
	bool out_of_bounds(const Parameters &model_parameters, const ParameterRec *par_info_ptr) const;
	bool get_derivative_parameters(const string &par_name, double derivative_par_value, const ParamTransformSeq &par_trans, const ParameterGroupInfo &group_info, const ParameterInfo &ctl_par_info,
		vector<double> &delta_numeric_par_vec, bool phiredswh_flag);
};

#endif /* JACOBIAN_1TO1H_ */

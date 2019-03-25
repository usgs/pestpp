#ifndef TORNADO_PLOT_H_
#define TORNADO_PLOT_H_

#include "GsaAbstractBase.h"

#include <vector>
#include <string>
#include <map>
#include <set>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include "GsaAbstractBase.h"
#include "Transformable.h"
#include "GsaAbstractBase.h"
#include "pest_data_structs.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class ParamTransformSeq;
class RunManagerAbstract;
class ModelRun;
class FileManager;
class RunningStats;
class TornadoPlot : public GsaAbstractBase
{
public:
	TornadoPlot(const std::vector<std::string> &_adj_par_name_vec, const Parameters &_fixed_ctl_pars, const Parameters &_init_pars,
		const Parameters &lower_bnd,
		const Parameters &upper_bnd, const std::set<std::string> &_log_trans_pars,
		ParamTransformSeq *base_partran_seq,
		const std::vector<std::string> &_obs_name_vec, FileManager *_file_manager_ptr,
		const ObservationInfo *_obs_info_ptr, bool _calc_obs_sen);
	~TornadoPlot();

	void assemble_runs(RunManagerAbstract &run_manager);
	void tornado_calc(RunManagerAbstract &run_manager, ModelRun model_run, std::ofstream &fout, const string obs_name = "");
	void calc_sen(RunManagerAbstract &run_manager, ModelRun model_run);
	bool calc_obs_sen;
	Parameters init_pars;
	const ObservationInfo *obs_info_ptr;
};
#endif //TORNADO_PLOT_H_


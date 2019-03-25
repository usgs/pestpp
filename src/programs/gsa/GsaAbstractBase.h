#ifndef GSA_ABSTRACT_BASE_H_
#define GSA_ABSTRACT_BASE_H_

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <set>
#include <random>
#include "Transformable.h"
#include "FileManager.h"
#include "ObjectiveFunc.h"
#include "OutputFileWriter.h"
#include "PerformanceLog.h"
#include "RunStorage.h"
using Eigen::MatrixXd;
using Eigen::VectorXd;

class ParamTransformSeq;
class RunManagerAbstract;
class ModelRun;
class FileManager;


class GsaAbstractBase
{
public:
	enum class PARAM_DIST{ normal, uniform };
	static const double MISSING_DATA;
	static mt19937_64 rand_engine;
	GsaAbstractBase(Pest &_pest_scenario,
		FileManager &_file_manager, ObjectiveFunc *_obj_func_ptr,
		const ParamTransformSeq &_par_transform, PARAM_DIST _par_dist, unsigned int seed);
	virtual void assemble_runs(RunManagerAbstract &run_manager) = 0;
	virtual void calc_sen(RunManagerAbstract &run_manager, ModelRun model_run) = 0;
	static std::map<std::string, std::string>  process_gsa_file(std::ifstream &fin, FileManager &file_manager);
	std::vector<double> calc_interval_midpoints(int n_interval, double min, double max);
	double ltqnorm(double p);
	static void set_seed(unsigned int _seed) { rand_engine.seed(_seed); }
	virtual ~GsaAbstractBase(void);

protected:
	PARAM_DIST par_dist;
	std::vector<std::string> adj_par_name_vec;
	std::vector<std::string> obs_name_vec;
	ParamTransformSeq gsa_parm_tran_seq;
	const ParamTransformSeq *base_partran_seq_ptr;
	FileManager *file_manager_ptr;
	ObjectiveFunc *obj_func_ptr;
	static void parce_line(const std::string &line, std::map<std::string, std::string> &arg_map);
	map<string, double> calc_parameter_norm_std_dev();
	map<string, double> calc_parameter_unif_std_dev();
	Parameters max_numeric_pars;
	Parameters min_numeric_pars;
	unsigned int seed;
};

std::ostream& operator<< (std::ostream& out, const std::vector<double> &rhs);
std::ostream& operator<< (std::ostream& out, const std::pair<double, double> &pr);

#endif /* GSA_ABSTRACT_BASE_H_ */

#ifndef SOBOL_H_
#define SOBOL_H_

#include <vector>
#include <string>
#include <Eigen/Dense>
#include "GsaAbstractBase.h"
#include "Transformable.h"
using Eigen::MatrixXd;
using Eigen::VectorXd;

class ParamTransformSeq;
class RunManagerAbstract;
class ModelRun;

class Sobol : public GsaAbstractBase
{
public:
	Sobol(Pest &_pest_scenario,
		FileManager &_file_manager, ObjectiveFunc *_obj_func_ptr,
		const ParamTransformSeq &_par_transform,
		int _n_sample, PARAM_DIST _par_dist, unsigned int _seed);
	void assemble_runs(RunManagerAbstract &run_manager);
	void calc_sen(RunManagerAbstract &run_manager, ModelRun model_run);
	pair<vector<double>,vector<double>> calc_sen_single(RunManagerAbstract &run_manager, ModelRun model_run, std::ofstream &fout_sbl, const std::string &obs_name);
	void calc_sen_single_old(RunManagerAbstract& run_manager, ModelRun model_run, std::ofstream& fout_sbl, const std::string& obs_name);

private:
	VectorXd gen_rand_vec(long nsample, double min, double max);
	void gen_m1_m2();
	MatrixXd gen_N_matrix(const MatrixXd &m1, const MatrixXd &m2, const vector<int> &idx_vec);
	void add_model_runs(RunManagerAbstract &run_manager, const MatrixXd &n, ofstream &f_out, string tag);
	vector<double> get_obs_vec(RunManagerAbstract &run_manager, int run_set, ModelRun &model_run, const string &obs_name );
	vector<double> get_phi_vec(RunManagerAbstract &run_manager, int run_set, ModelRun &model_run);
	int n_sample;
	void process_runs(RunManagerAbstract& run_manager, ModelRun &model_run);
	map<int, Observations> run_map;
	Eigen::MatrixXd m1;
	Eigen::MatrixXd m2;
};
#endif /* SOBOL_H_ */

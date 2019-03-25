#ifndef MORRISMETHOD_H_
#define MORRISMETHOD_H_

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

class MorrisObsSenFile
{
	friend class MorrisMethod;
public:
	void initialize(const std::vector<std::string> &par_names_vec, const std::vector<std::string> &obs_names_vec, double _no_data, const GsaAbstractBase *_gsa_abstract_base);
	void add_sen_run_pair(const std::string &par_name, double p1, Observations &obs1, double p2, Observations &obs2);
	void calc_pooled_obs_sen(std::ofstream &fout_obs_sen, map<string, double> &obs_2_sen_weight, map<string, double> &par_2_sen_weight);
private:
	double no_data;
	vector<string> par_names_vec;
	vector<string> obs_names_vec;
	const GsaAbstractBase *gsa_abstract_base;
	map< pair<string,string>, RunningStats> map_obs_stats;
	map<string, int> parname_to_indexmap;
};

class MorrisMethod : public GsaAbstractBase
{
public:
	MorrisMethod(Pest &_pest_scenario,
		FileManager &_file_manager, ObjectiveFunc *_obj_func_ptr,
		const ParamTransformSeq &_par_transform,
		int _p, int _r, double _delta, bool _calc_pooled_obs,
		bool _calc_morris_obs_sen, PARAM_DIST _par_dist, unsigned int _seed);
	void process_pooled_var_file();
	void initialize(int _p, int _r, double _delta);
	void assemble_runs(RunManagerAbstract &run_manager);
	void calc_sen(RunManagerAbstract &run_manager, ModelRun model_run);
	void calc_morris_obs(std::ostream &, MorrisObsSenFile &morris_sen_file);
	~MorrisMethod(void);
private:
	int p; // number of levels for each parameters
	int r;
	bool calc_morris_obs_sen;
	bool calc_obs_sen;
	double delta;
	MatrixXd b_star_mat;
	static MatrixXd create_B_mat(int k);
	static MatrixXd create_J_mat(int k);
	static MatrixXd create_D_mat(int k);
	static MatrixXd create_P_mat(int k);
		MatrixXd create_P_star_mat(int k);
	VectorXd create_x_vec(int k);
	Parameters get_numeric_parameters(int row);
	static int rand_plus_minus_1(void);
	MorrisObsSenFile obs_sen_file;
	const ObservationInfo *obs_info_ptr;
	map<std::string, std::string> group_2_pool_group_map;
};




#endif /* MORRISMETHOD_H_ */

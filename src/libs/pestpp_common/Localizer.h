#ifndef LOCALIZER_H_
#define LOCALIZER_H_

#include <map>
#include <random>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "FileManager.h"
#include "ObjectiveFunc.h"
#include "OutputFileWriter.h"
#include "PerformanceLog.h"
#include "RunStorage.h"
#include "covariance.h"
#include "RunManagerAbstract.h"
#include "PerformanceLog.h"
#include "Ensemble.h"

class AutoAdaLocThread
{
public:

	AutoAdaLocThread(PerformanceLog *_performance_log, ofstream *_f_out, int _iter, int _ies_verbose, int _npar, int _nobs, vector<int> &_par_indices,
		Eigen::MatrixXd &_pe_diff, Eigen::MatrixXd &_oe_diff, Eigen::ArrayXd &_par_std, Eigen::ArrayXd &_obs_std, vector<string> &_par_names, vector<string> &_obs_names,
		vector<Eigen::Triplet<double>> &_triplets, double _sigma_dist, map<string,set<string>> &_list_obs);

	//Eigen::DiagonalMatrix<double, Eigen::Dynamic> get_matrix_from_map(vector<string> &names, map<string, double> &dmap);	
	//Eigen::MatrixXd get_matrix_from_map(int num_reals, vector<string> &names, map<string, Eigen::VectorXd> &emap);


	void work(int thread_id);


private:
	int npar, nobs, ies_verbose, iter;
	int nzero_par, nzero_obs;
	double sigma_dist;
	Eigen::MatrixXd &pe_diff, &oe_diff;
	vector<int> &par_indices;
	Eigen::ArrayXd &par_std, &obs_std;
	vector<string> &par_names, &obs_names;
	vector<Eigen::Triplet<double>> &triplets;
	PerformanceLog *performance_log;
	ofstream *f_out;
	map<string, set<string>> list_obs;
	map<int, string> idx2obs;
	mutex pe_diff_lock, oe_diff_lock, par_indices_lock, pfm_lock;
	mutex f_out_lock, par_names_lock, obs_names_lock;
	mutex triplets_lock;
};



class Localizer
{
public:
	enum How { PARAMETERS, OBSERVATIONS};
	Localizer() { ; }
	Localizer(Pest *_pest_scenario_ptr) { pest_scenario_ptr = _pest_scenario_ptr; }
	bool initialize(PerformanceLog *performance_log);
	unordered_map<string, pair<vector<string>, vector<string>>> get_localizer_map(int iter, ObservationEnsemble &oe, ParameterEnsemble &pe, PerformanceLog *performance_log);// { return localizer_map; }
	void set_pest_scenario(Pest *_pest_scenario_ptr) { pest_scenario_ptr = _pest_scenario_ptr; }
	Eigen::MatrixXd get_localizing_obs_hadamard_matrix(int num_reals,string col_name,vector<string> &obs_names);
	Eigen::MatrixXd get_localizing_par_hadamard_matrix(int num_reals, string row_name, vector<string> &par_names);
	How get_how() { return how; }
	bool get_use() { return use; }
	bool get_autoadaloc() { return autoadaloc; }
	string get_filename() { return filename;  }
	int get_num_upgrade_steps() { return localizer_map.size(); }
	void report(ofstream &f_rec);
private:
	bool use;
	bool autoadaloc;
	double sigma_dist;
	How how;
	Pest * pest_scenario_ptr;
	Mat mat;
	string filename;
	unordered_map<string,pair<vector<string>, vector<string>>> localizer_map;
	map<string, set<string>> listed_obs;
	map<string, int> obs2row_map, par2col_map;

	void process_mat(PerformanceLog *performance_log);	
};

void aal_upgrade_thread_function(int id, AutoAdaLocThread &worker, exception_ptr &eptr);

#endif

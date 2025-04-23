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
	int npar, nobs, ies_verbose, iter, par_count;
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
	enum LocTyp {COVARIANCE, LOCALANALYSIS };
	Localizer() { initialized=false; }
	Localizer(Pest* _pest_scenario_ptr) { pest_scenario_ptr = _pest_scenario_ptr; initialized = false; }
	bool initialize(PerformanceLog *performance_log, ofstream& frec, bool forgive_missing=false);
	unordered_map<string, pair<vector<string>, vector<string>>> get_localanalysis_case_map(int iter, vector<string>& act_obs_names, vector<string>& act_par_names, 
		ObservationEnsemble &oe, ParameterEnsemble &pe, PerformanceLog *performance_log, ofstream& frec);// { return localizer_map; }
	
	void set_pest_scenario(Pest *_pest_scenario_ptr) { pest_scenario_ptr = _pest_scenario_ptr; }
	Eigen::MatrixXd get_obsdiff_hadamard_matrix(int num_reals,string col_name,vector<string> &obs_names);
	Eigen::MatrixXd get_pardiff_hadamard_matrix(int num_reals, string row_name, vector<string> &par_names);
	//Eigen::MatrixXd get_kalmangain_hadamard_matrix(vector<string>& obs_names, vector<string>& par_names);
	Eigen::VectorXd get_obs_hadamard_vector(string par_name, vector<string>& obs_names);

	How get_how() { return how; }
	bool get_use() { return use; }
	bool get_autoadaloc() { return autoadaloc; }
	string get_filename() { return filename;  }
	int get_num_upgrade_steps() { return _localizer_map.size(); }
	LocTyp get_loctyp() { return loctyp; }
	void report(ofstream &f_rec);
	bool is_initialized() { return initialized; }
	Mat* get_orgmat_ptr() {return &org_mat;}
private:
	bool use;
	bool autoadaloc;
	double sigma_dist;
	bool initialized;
	How how;
	LocTyp loctyp;
	Pest * pest_scenario_ptr;
	Mat org_mat, cur_mat;
	string filename;
	unordered_map<string,pair<vector<string>, vector<string>>> _localizer_map;
	//map<string, set<string>> listed_obs;
	map<string, int> obs2row_map, par2col_map,colname2col_map, rowname2row_map;

	unordered_map<string, pair<vector<string>, vector<string>>> process_mat(PerformanceLog *performance_log, Mat& mat, bool forgive_missing=false);
	void update_obs_info_from_mat(Mat& mat, vector<vector<string>>& obs_map, vector<string>& missing, vector<string>& dups, 
		set<string>& obs_names, map<string, vector<string>>& obgnme_map, vector<string>& not_allowed, bool forgive_missing=false);
	void update_par_info_from_mat(Mat& mat, vector<vector<string>>& par_map, vector<string>& missing, vector<string>& dups,
		set<string>& par_names, map<string, vector<string>>& pargp_map, vector<string>& not_allowed);

};

void aal_upgrade_thread_function(int id, AutoAdaLocThread &worker, exception_ptr &eptr);

#endif

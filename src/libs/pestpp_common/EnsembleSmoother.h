#ifndef ENSEMBLESMOOTHER_H_
#define ENSEMBLESMOOTHER_H_

#include <map>
#include <random>
#include <mutex>
#include <thread>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "FileManager.h"
#include "ObjectiveFunc.h"
#include "OutputFileWriter.h"
#include "PerformanceLog.h"
#include "RunStorage.h"
#include "covariance.h"
#include "RunManagerAbstract.h"
#include "ObjectiveFunc.h"
#include "Localizer.h"


class PhiHandler
{
public:

	enum phiType { MEAS, COMPOSITE, REGUL, ACTUAL };
	PhiHandler() { ; }
	PhiHandler(Pest *_pest_scenario, FileManager *_file_manager,
		       ObservationEnsemble *_oe_base, ParameterEnsemble *_pe_base,
		       Covariance *_parcov, double *_reg_factor, ObservationEnsemble *_weights);
	void update(ObservationEnsemble &oe, ParameterEnsemble &pe, bool include_regul=true);
	double get_mean(phiType pt);
	double get_std(phiType pt);
	double get_max(phiType pt);
	double get_min(phiType pt);

	double calc_mean(map<string, double> *phi_map);
	double calc_std(map<string, double> *phi_map);

	map<string, double>* get_phi_map(PhiHandler::phiType &pt);
	void report(bool echo=true, bool include_regul=true);
	void write(int iter_num, int total_runs, bool write_group = true);
	void write_group(int iter_num, int total_runs, vector<double> extra);
	vector<int> get_idxs_greater_than(double bad_phi, double bad_phi_sigma, ObservationEnsemble &oe);

	Eigen::MatrixXd get_obs_resid(ObservationEnsemble &oe);
	Eigen::MatrixXd get_obs_resid_subset(ObservationEnsemble &oe);

	Eigen::MatrixXd get_par_resid(ParameterEnsemble &pe);
	Eigen::MatrixXd get_par_resid_subset(ParameterEnsemble &pe);
	Eigen::MatrixXd get_actual_obs_resid(ObservationEnsemble &oe);
	Eigen::VectorXd get_q_vector();
	vector<string> get_lt_obs_names() { return lt_obs_names; }
	vector<string> get_gt_obs_names() { return gt_obs_names; }

	void apply_ineq_constraints(Eigen::MatrixXd &resid, vector<string> &names);

private:
	map<string, double> get_summary_stats(phiType pt);
	string get_summary_string(phiType pt);
	string get_summary_header();
	void prepare_csv(ofstream &csv,vector<string> &names);
	void prepare_group_csv(ofstream &csv, vector<string> extra = vector<string>());

	map<string, Eigen::VectorXd> calc_meas(ObservationEnsemble &oe, Eigen::VectorXd &_q_vec);
	map<string, Eigen::VectorXd> calc_regul(ParameterEnsemble &pe);// , double _reg_fac);
	map<string, Eigen::VectorXd> calc_actual(ObservationEnsemble &oe, Eigen::VectorXd &_q_vec);
	map<string, double> calc_composite(map<string,double> &_meas, map<string,double> &_regul);
	//map<string, double>* get_phi_map(PhiHandler::phiType &pt);
	void write_csv(int iter_num, int total_runs,ofstream &csv, phiType pt,
		           vector<string> &names);
	void write_group_csv(int iter_num, int total_runs, ofstream &csv,
		vector<double> extra = vector<double>());

	double *reg_factor;
	double org_reg_factor;
	Eigen::VectorXd org_q_vec;
	vector<string> oreal_names,preal_names;
	Pest* pest_scenario;
	FileManager* file_manager;
	ObservationEnsemble* oe_base;
	ParameterEnsemble* pe_base;
	ObservationEnsemble *weights;
	//Covariance parcov_inv;
	Eigen::VectorXd parcov_inv_diag;
	map<string, double> meas;
	map<string, double> regul;
	map<string, double> composite;
	map<string, double> actual;

	vector<string> lt_obs_names;
	vector<string> gt_obs_names;

	map<string, vector<int>> obs_group_idx_map;
	map<string, vector<int>> par_group_idx_map;
	map<string, map<string, double>> obs_group_phi_map, par_group_phi_map;

	map<string, double> get_obs_group_contrib(Eigen::VectorXd &phi_vec);
	map<string, double> get_par_group_contrib(Eigen::VectorXd &phi_vec);

};

class ParChangeSummarizer
{
public:
	ParChangeSummarizer() { ; }
	ParChangeSummarizer(ParameterEnsemble *_base_pe_ptr, FileManager *_file_manager_ptr);
	void summarize(ParameterEnsemble &pe);

private:
	ParameterEnsemble * base_pe_ptr;
	FileManager *file_manager_ptr;
	map<string, set<string>> pargp2par_map;
	pair<map<string,double>, map<string, double>> init_moments;


};


class LocalUpgradeThread
{
public:

	LocalUpgradeThread(PerformanceLog *_performance_log, unordered_map<string, Eigen::VectorXd> &_par_resid_map, unordered_map<string, Eigen::VectorXd> &_par_diff_map,
		unordered_map<string, Eigen::VectorXd> &_obs_resid_map, unordered_map<string, Eigen::VectorXd> &_obs_diff_map,
		Localizer &_localizer, unordered_map<string, double> &_parcov_inv_map,
		unordered_map<string, double> &_weight_map, ParameterEnsemble &_pe_upgrade,
		unordered_map<string, pair<vector<string>, vector<string>>> &_cases,
		unordered_map<string, Eigen::VectorXd> &_Am_map, Localizer::How &_how);

	//Eigen::DiagonalMatrix<double, Eigen::Dynamic> get_matrix_from_map(vector<string> &names, map<string, double> &dmap);	
	//Eigen::MatrixXd get_matrix_from_map(int num_reals, vector<string> &names, map<string, Eigen::VectorXd> &emap);


	void work(int thread_id, int iter, double cur_lam);


private:
	PerformanceLog * performance_log;
	Localizer::How how;
	vector<string> keys;
	int count, total;
	//double eigthresh, cur_lam;
	//int maxsing, num_reals,iter, thread_id;
	//bool use_approx, use_prior_scaling;

	unordered_map<string, pair<vector<string>, vector<string>>> &cases;

	ParameterEnsemble &pe_upgrade;
	//PhiHandler &ph;
	Localizer &localizer;
	unordered_map<string, double> &parcov_inv_map;
	unordered_map<string, double> &weight_map;

	unordered_map<string, Eigen::VectorXd> &par_resid_map, &par_diff_map, &Am_map;
	unordered_map<string, Eigen::VectorXd> &obs_resid_map, &obs_diff_map;

	mutex ctrl_lock, weight_lock, loc_lock, parcov_lock;
	mutex obs_resid_lock, obs_diff_lock, par_resid_lock;
	mutex par_diff_lock, am_lock, put_lock;
	mutex next_lock;
	
};


class IterEnsembleSmoother
{
public:
	IterEnsembleSmoother(Pest &_pest_scenario, FileManager &_file_manager,
		OutputFileWriter &_output_file_writer, PerformanceLog *_performance_log,
		RunManagerAbstract* _run_mgr_ptr);
	void initialize();
	void iterate_2_solution();
	void pareto_iterate_2_solution();
	void finalize();
	void throw_ies_error(string message);
	bool should_terminate();

private:
	int  verbose_level;
	Pest &pest_scenario;
	FileManager &file_manager;
	OutputFileWriter &output_file_writer;
	PerformanceLog *performance_log;
	RunManagerAbstract* run_mgr_ptr;
	PhiHandler ph;
	ParChangeSummarizer pcs;
	Covariance parcov, obscov;
	double reg_factor;

	string base_name = "BASE"; //this is also defined in Ensemble

	bool use_localizer;
	Localizer localizer;

	int num_threads;

	set<string> pp_args;

	int iter,subset_size;
	bool use_subset;

	double last_best_lam, last_best_mean,last_best_std;
	vector<double> best_mean_phis;
	double best_phi_yet;

	int consec_bad_lambda_cycles;

	double lambda_max, lambda_min;
	int warn_min_reals, error_min_reals;
	vector<double> lam_mults;
	map<string, double> pareto_obs;
	map<string, double> pareto_weights;
	//string fphi_name;
	//ofstream fphi;
	vector<string> oe_org_real_names, pe_org_real_names;
	vector<string> act_obs_names, act_par_names;
	vector<int> subset_idxs;

	ParameterEnsemble pe, pe_base;
	ObservationEnsemble oe, oe_base, weights;
	//Eigen::MatrixXd prior_pe_diff;
	//Eigen::MatrixXd Am;
	Eigen::DiagonalMatrix<double,Eigen::Dynamic> obscov_inv_sqrt, parcov_inv_sqrt;

	//bool solve_old();
	bool solve_new();
	void adjust_pareto_weight(string &obsgroup, double wfac);

	//ParameterEnsemble calc_upgrade(vector<string> &obs_names, vector<string> &par_names,double lamb, int num_reals);

	//ParameterEnsemble calc_localized_upgrade(double cur_lam);
	ParameterEnsemble calc_localized_upgrade_threaded(double cur_lam, unordered_map<string, pair<vector<string>, vector<string>>> &loc_map);

	//EnsemblePair run_ensemble(ParameterEnsemble &_pe, ObservationEnsemble &_oe);
	vector<int> run_ensemble(ParameterEnsemble &_pe, ObservationEnsemble &_oe, const vector<int> &real_idxs=vector<int>());
	vector<ObservationEnsemble> run_lambda_ensembles(vector<ParameterEnsemble> &pe_lams, vector<double> &lam_vals, vector<double> &scale_vals);
	//map<string, double> get_phi_vec_stats(map<string,PhiComponets> &phi_info);
	//map<string,PhiComponets> get_phi_info(ObservationEnsemble &_oe);
	void report_and_save();
	void save_mat(string prefix, Eigen::MatrixXd &mat);
	bool initialize_pe(Covariance &cov);
	bool initialize_oe(Covariance &cov);
	void initialize_restart();
	void initialize_weights();
	void initialize_parcov();
	void initialize_obscov();
	void drop_bad_phi(ParameterEnsemble &_pe, ObservationEnsemble &_oe, bool is_subset=false);
	//void check_ensembles(ObservationEnsemble &oe, ParameterEnsemble &pe);
	template<typename T, typename A>
	void message(int level, const string &_message, vector<T, A> _extras, bool echo=true);
	void message(int level, const string &_message);

	//template<typename T, typename A>
	//void message(int level, char* _message, vector<T, A> _extras);// { message(level, string(_message), _extras); }
	//void message(int level, char* _message);// { message(level, string(_message)); }

	template<typename T>
	void message(int level, const string &_message, T extra);

	//template<typename T>
	//void message(int level, char* _message, T extra);

	void sanity_checks();

	void add_bases();

	void update_reals_by_phi(ParameterEnsemble &_pe, ObservationEnsemble &_oe);


	//map<int,int> get_subset_idx_map();
	void set_subset_idx(int size);
	Eigen::MatrixXd get_Am(const vector<string> &real_names, const vector<string> &par_names);

};

#endif

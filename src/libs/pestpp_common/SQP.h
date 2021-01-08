#ifndef SQP_H_
#define SQP_H_

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
#include "EnsembleMethodUtils.h"
#include "constraints.h"



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

class SqpFilter
{
public:
	SqpFilter(bool _minimize=true,double _obj_tol = 0.01, double _viol_tol = 0.01) {
		minimize = _minimize; obj_tol = _obj_tol; viol_tol = _viol_tol;
	}
	bool accept(double obj_val, double violation_val);
	bool update(double obj_val, double violation_val);

private:
	bool minimize;
	double obj_tol;
	double viol_tol;

	set<pair<double, double>> obj_viol_pairs;

	bool first_dominates_second(const pair<double, double>& first, const pair<double, double>& second);
	
};

class SeqQuadProgram
{
public:
	SeqQuadProgram(Pest& _pest_scenario, FileManager& _file_manager,
		OutputFileWriter& _output_file_writer, PerformanceLog* _performance_log,
		RunManagerAbstract* _run_mgr_ptr);
	
	void initialize();
	void iterate_2_solution();
	void finalize();
	void throw_sqp_error(string message);
	bool should_terminate();

private:
	int  verbose_level;
	Pest &pest_scenario;
	FileManager &file_manager;
	std::mt19937 rand_gen;
	std::mt19937 subset_rand_gen;
	OutputFileWriter &output_file_writer;
	PerformanceLog *performance_log;
	RunManagerAbstract* run_mgr_ptr;
	L2PhiHandler ph;
	ParChangeSummarizer pcs;
	Covariance parcov, obscov;
	double reg_factor;
	chancePoints chancepoints;
	string obj_func_str;
	string obj_obs;
	string obj_sense;
	bool use_obj_obs;
	bool use_obj_pi;
	map<string, double> obj_func_coef_map;

	int num_threads;

	double eigthresh;
	vector<double> scale_vals;
	set<string> pp_args;

	int iter;

	double last_best;
	vector<double> best_phis;
	double best_phi_yet;

	int warn_min_reals, error_min_reals;
	
	vector<string> oe_org_real_names, pe_org_real_names;
	vector<string> act_obs_names, act_par_names;
	vector<string> dv_names;
	//vector<int> subset_idxs;
	
	Parameters current_ctl_dv_values;
	Observations current_obs;

	Parameters current_grad_vector;
	map<int, Parameters> grad_vector_map;

	ParameterEnsemble dv, dv_base;
	ObservationEnsemble oe, oe_base;

	//these are used so that we can update the constraints based on the current best values
	//Parameters best_mean_dv_values;
	//Observations best_mean_obs_values;

	Constraints constraints;

	bool oe_drawn, dv_drawn;

	bool use_ensemble_grad;

	Jacobian_1to1 jco;

	//store the hessian as a cov since it is symmetric...
	Covariance hessian;

	SqpFilter filter;

	void prep_4_ensemble_grad();
	void prep_4_fd_grad();

	bool update_hessian_and_grad_vector();

	bool solve_new();

	bool pick_candidate_and_update_current(ParameterEnsemble& dv_candidates, ObservationEnsemble& _oe);

	Parameters calc_gradient_vector(const Parameters& _current_dv_);

	Eigen::VectorXd get_obj_vector(ParameterEnsemble& _dv, ObservationEnsemble& _oe);
	
	double get_obj_value(Parameters& _current_ctl_dv_vals, Observations& _current_obs);
	map<string, double> get_obj_map(ParameterEnsemble& _dv, ObservationEnsemble& _oe);

	Eigen::VectorXd calc_search_direction_vector(const Parameters& _current_dv_, Eigen::VectorXd &grad_vector);  // to check

	Parameters fancy_solve_routine(double scale_val, const Parameters& _current_dv_);

	vector<int> run_ensemble(ParameterEnsemble &_pe, ObservationEnsemble &_oe, const vector<int> &real_idxs=vector<int>());
	ObservationEnsemble run_candidate_ensemble(ParameterEnsemble&dv_candidates, vector<double> &scale_vals);

	void run_jacobian(Parameters& _current_dv_vals,Observations& _current_obs, bool init_obs);

	void make_gradient_runs(Parameters& _current_dv_vals, Observations& _current_obs);

	void report_and_save_ensemble();
	void report_and_save_ensemble(ParameterEnsemble& _dv, ObservationEnsemble& _oe);
	void save(ParameterEnsemble& _dv, ObservationEnsemble& _oe, bool save_base=true);
	void save_mat(string prefix, Eigen::MatrixXd &mat);
	bool initialize_dv(Covariance &cov);
	//bool initialize_oe(Covariance &cov);
	bool initialize_restart();
	void initialize_parcov();
	void initialize_obscov();
	void initialize_objfunc();
	void drop_bad_phi(ParameterEnsemble &_pe, ObservationEnsemble &_oe, bool is_subset=false);
	
	void queue_chance_runs();
	
	template<typename T, typename A>
	void message(int level, const string &_message, vector<T, A> _extras, bool echo=true);
	void message(int level, const string &_message);

	template<typename T>
	void message(int level, const string &_message, T extra);

	void sanity_checks();

	void add_bases();

	void update_reals_by_phi(ParameterEnsemble &_pe, ObservationEnsemble &_oe);

	void set_subset_idx(int size);
	
};

#endif

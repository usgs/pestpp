#ifndef ENSEMBLEMETHODUTILS_H_
#define ENSEMBLEMETHODUTILS_H_

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
#include "network_package.h"

enum chancePoints { ALL, SINGLE };


class MmNeighborThread
{
public:
    MmNeighborThread(unordered_map<string,Eigen::VectorXd>& _real_vec_map,
                     unordered_map<string,vector<int>>& _mm_real_idx_map,
                     unordered_map<string,pair<vector<string>,vector<string>>>& _mm_real_name_map,
                     unordered_map<string,unordered_map<string,double>>& _neighbor_phi_map,
                     unordered_map<string,unordered_map<string,double>>& _neighbor_pardist_map);

    void work(int tid, int verbose_level, double mm_alpha, map<string,map<string,double>> weight_phi_map,
              vector<string> preal_names, vector<string> oreal_names,map<string,int> real_map,
              Eigen::SparseMatrix<double> parcov_inv);

protected:
    unordered_map<string, Eigen::VectorXd>& real_vec_map;
    unordered_map<string,unordered_map<string,double>>& neighbor_phi_map;
    unordered_map<string,unordered_map<string,double>>& neighbor_pardist_map;
    vector<int> indexes;
    int count, total;

    unordered_map<string,vector<int>>& mm_real_idx_map;
    unordered_map<string,pair<vector<string>,vector<string>>>& mm_real_name_map;

    mutex next_lock, results_lock;

};

class PhiThread
{
public:
    explicit PhiThread(vector<string> _oe_real_names);

    void work(int thread_id, Eigen::MatrixXd& weights, Eigen::MatrixXd& resid, vector<string>& oe_real_names, map<string,map<string,double>>& phi_map);

protected:
    vector<string> keys;
    int count, total;

    mutex next_lock, phi_map_lock;
};

class L2PhiHandler
{
public:

	enum phiType { MEAS, COMPOSITE, REGUL, ACTUAL, NOISE };
	L2PhiHandler() { ; }
	L2PhiHandler(Pest *_pest_scenario, FileManager *_file_manager,
		       ObservationEnsemble *_oe_base, ParameterEnsemble *_pe_base,
		       Covariance *_parcov, bool should_prep_csv = true, string _tag=string());
	void update(ObservationEnsemble &oe, ParameterEnsemble &pe);
    void update(ObservationEnsemble &oe, ParameterEnsemble &pe, ObservationEnsemble& weights);
	double get_mean(phiType pt);
	double get_std(phiType pt);
	double get_max(phiType pt);
	double get_min(phiType pt);
    double get_representative_phi(phiType pt);

	double calc_mean(map<string, double> *phi_map);
	double calc_std(map<string, double> *phi_map);

	static double calc_median(const std::vector<double> &values);

	double calc_iqr_thresh(map<string, double> *phi_map, double bad_phi_sigma);

	map<string, double>* get_phi_map_ptr(L2PhiHandler::phiType pt);
	map<string, double> get_phi_map(L2PhiHandler::phiType pt);
	void report(bool echo=true, bool group_report=true);
	void report_group(bool echo=true);

	void write(int iter_num, int total_runs, bool write_group = true);
	void write_group(int iter_num, int total_runs, vector<double> extra);
    //csv << "iteration,num_reals,current_lambda,accept_phi,lambda_mult,lambda_scale_fac,lambda,meas_phi" << endl;
    void write_lambda(int iteration,int num_reals,double current_lambda,double current_comp_mean_phi,double current_comp_std_phi,
                      double lambda_mult,double lambda, double comp_mean_phi, double comp_std_phi);
	vector<int> get_idxs_greater_than(double bad_phi, double bad_phi_sigma, ObservationEnsemble &oe, ObservationEnsemble& weights);

	Eigen::MatrixXd get_obs_resid(ObservationEnsemble &oe, bool apply_ineq=true);
	Eigen::MatrixXd get_obs_resid_subset(ObservationEnsemble &oe, bool apply_ineq=true,vector<string> real_names=vector<string>());

	Eigen::MatrixXd get_par_resid(ParameterEnsemble &pe);
	Eigen::MatrixXd get_par_resid_subset(ParameterEnsemble &pe,vector<string> real_names=vector<string>());
	Eigen::MatrixXd get_actual_obs_resid(ObservationEnsemble &oe);
    //Eigen::MatrixXd get_noise_resid(ObservationEnsemble &oe, bool apply_ineq=true);
	Eigen::VectorXd get_q_vector();
	vector<string> get_lt_obs_names() { return lt_obs_names; }
	vector<string> get_gt_obs_names() { return gt_obs_names; }
    map<string,double> get_lt_obs_bounds() {return lt_obs_bounds;}
    map<string,double> get_gt_obs_bounds() {return gt_obs_bounds;}
    map<string,pair<double,double>> get_double_obs_bounds() {return double_obs_bounds;}

	void apply_ineq_constraints(Eigen::MatrixXd &resid, Eigen::MatrixXd &sim_vals, vector<string> &names);

	void save_residual_cov(ObservationEnsemble& oe, int iter);

	//map<string,double> get_meas_phi(ObservationEnsemble& oe, Eigen::VectorXd& q_vec);

	map<string,map<string,double>> get_swr_real_map(ObservationEnsemble& oe, ObservationEnsemble& weights,phiType ptype=phiType::MEAS);

    map<string,double> get_swr_map(ObservationEnsemble& oe, string real_name= "",phiType ptype=phiType::MEAS);
	map<string,map<string,double>> get_meas_phi_weight_ensemble(ObservationEnsemble& oe, ObservationEnsemble& weights);

    vector<string> get_violating_realizations(ObservationEnsemble& oe, const vector<string>& viol_obs_names);
    vector<string> detect_simulation_data_conflict(ObservationEnsemble& _oe, string csv_tag);

private:
	string tag;
	map<string, double> get_summary_stats(phiType pt);
	string get_summary_string(phiType pt);
	string get_summary_header();
	void prepare_csv(ofstream &csv,vector<string> &names);
	void prepare_group_csv(ofstream &csv, vector<string> extra = vector<string>());
    void prepare_lambda_csv(ofstream &csv);

	map<string, Eigen::VectorXd> calc_meas(ObservationEnsemble &oe, Eigen::VectorXd& q_vec);
    map<string, Eigen::VectorXd> calc_meas(ObservationEnsemble &oe, ObservationEnsemble& weights);
    //map<string, Eigen::VectorXd> calc_noise(ObservationEnsemble &oe, ObservationEnsemble& weights);
    //map<string, Eigen::VectorXd> calc_noise(ObservationEnsemble &oe, Eigen::VectorXd& q_vec);
	map<string, Eigen::VectorXd> calc_regul(ParameterEnsemble &pe);// , double _reg_fac);
	map<string, Eigen::VectorXd> calc_actual(ObservationEnsemble &oe, Eigen::VectorXd& q_vec);
    map<string, Eigen::VectorXd> calc_actual(ObservationEnsemble & oe, ObservationEnsemble& weights);

	map<string, double> calc_composite(map<string,double> &_meas, map<string,double> &_regul);
	//map<string, double>* get_phi_map(PhiHandler::phiType &pt);
	void write_csv(int iter_num, int total_runs,ofstream &csv, phiType pt,
		           vector<string> &names);
	void write_group_csv(int iter_num, int total_runs, ofstream &csv,
		vector<double> extra = vector<double>());

	double org_reg_factor;
	vector<string> oreal_names,preal_names;
	Pest* pest_scenario;
	FileManager* file_manager;
	ObservationEnsemble* oe_base;
	ParameterEnsemble* pe_base;
	//Covariance parcov_inv;
	Eigen::VectorXd parcov_inv_diag;
	map<string, double> meas;
	map<string, double> regul;
	map<string, double> composite;
	map<string, double> actual;
    map<string, double> noise;
    map<string,int> num_conflict_group;

	vector<string> lt_obs_names;
	vector<string> gt_obs_names;
    map<string,double> lt_obs_bounds;
    map<string,double> gt_obs_bounds;
    map<string,pair<double,double>> double_obs_bounds;

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
	ParChangeSummarizer(ParameterEnsemble *_base_pe_ptr, FileManager *_file_manager_ptr, OutputFileWriter* _output_file_writer_ptr);
	void summarize(ParameterEnsemble &pe, string filename = string());
private:
	double cv_dec_threshold = 0.3;
	ParameterEnsemble * base_pe_ptr;
	FileManager *file_manager_ptr;
	OutputFileWriter* output_file_writer_ptr;
	map<string, set<string>> pargp2par_map;
	pair<map<string,double>, map<string, double>> init_moments;
	map<string, double> mean_change;
	map<string, double> std_change;
	map<string, double> init_cv;
	map<string, double> curr_cv;
	map<string, int> num_at_ubound;
	map<string, double> percent_at_ubound;
    map<string, int> num_at_lbound;
    map<string, double> percent_at_lbound;

	void update(ParameterEnsemble& pe);
	void write_to_csv(string& filename);
    map<string,int> get_npar_per_group_with_excess_std_reduction(ParameterEnsemble& _pe, double thresh=0.95);


};

pair<Parameters,Observations> save_real_par_rei(Pest& pest_scenario, ParameterEnsemble& pe, ObservationEnsemble& oe,
    OutputFileWriter& output_file_writer, FileManager& file_manager, int iter, string tag = BASE_REAL_NAME,
    int cycle = NetPackage::NULL_DA_CYCLE,map<string,double> base_weights=map<string,double>());

vector<int> run_ensemble_util(PerformanceLog* performance_log, ofstream& frec, ParameterEnsemble& _pe,
	ObservationEnsemble& _oe, RunManagerAbstract* run_mgr_ptr,
	bool check_pe_consistency = false, const vector<int>& real_idxs = vector<int>(),int da_cycle=NetPackage::NULL_DA_CYCLE,
	string additional_tag="");

class EnsembleSolver
{
public:
	EnsembleSolver(PerformanceLog* _performance_log, FileManager& _file_manager, Pest& _pest_scenario, ParameterEnsemble& _pe,
		ObservationEnsemble& _oe, ObservationEnsemble& _base_oe, ObservationEnsemble& _weights, Localizer& _localizer,
		Covariance& _parcov,Eigen::MatrixXd& _Am, L2PhiHandler& _ph,
		bool _use_localizer, int _iter, vector<string>& _act_par_names, vector<string> &_act_obs_names,
		double _reg_factor);

	void solve(int num_threads, double cur_lam, bool use_glm_form, ParameterEnsemble& pe_upgrade, unordered_map<string, pair<vector<string>, vector<string>>>& loc_map);
    void solve_multimodal(int num_threads, double cur_lam, bool use_glm_form, ParameterEnsemble& pe_upgrade, unordered_map<string,pair<vector<string>, vector<string>>>& loc_map, double mm_alpha);
    void update_multimodal_components(const double mm_alpha);


private:
	PerformanceLog* performance_log;
	FileManager& file_manager;
	int iter, verbose_level;
	bool use_localizer;
	double reg_factor;
	Pest& pest_scenario;
	ParameterEnsemble& pe;
	ObservationEnsemble& oe, base_oe, weights;
	Localizer& localizer;
	Covariance& parcov;
	Eigen::MatrixXd& Am;
	L2PhiHandler& ph;
	unordered_map<string, Eigen::VectorXd> par_resid_map, obs_resid_map, Am_map;
	unordered_map<string, Eigen::VectorXd> par_diff_map, obs_diff_map, obs_err_map;
	unordered_map<string, double> weight_map;
	unordered_map<string, double> parcov_inv_map;
	unordered_map<string,vector<int>> mm_real_idx_map;
	unordered_map<string,pair<vector<string>,vector<string>>> mm_real_name_map;
    unordered_map<string,Eigen::VectorXd> mm_q_vec_map;
	//unordered_map<string, pair<vector<string>, vector<string>>> loc_map;
	vector<string>& act_par_names, act_obs_names;
	template<typename T, typename A>
	void message(int level, const string& _message, vector<T, A> _extras, bool echo = true);
	void message(int level, const string& _message);

	template<typename T>
	void message(int level, const string& _message, T extra);

	void initialize_for_localized_solve(string center_on = string(), vector<int> real_idxs=vector<int>());
    void initialize_for_mm_solve();
	void nonlocalized_solve(double cur_lam,bool use_glm_form, ParameterEnsemble& pe_upgrade,
                         string center_on=string(), vector<int> real_idxs=vector<int>(),Eigen::VectorXd q_vec=Eigen::VectorXd());

};

class MmUpgradeThread
{
public:
    MmUpgradeThread(PerformanceLog* _performance_log, unordered_map<string, Eigen::VectorXd>& _par_resid_map,
                    unordered_map<string, Eigen::VectorXd>& _par_diff_map,
                  unordered_map<string, Eigen::VectorXd>& _obs_resid_map, unordered_map<string, Eigen::VectorXd>& _obs_diff_map,
                  unordered_map<string, Eigen::VectorXd>& _obs_err_map,
                  unordered_map<string, Eigen::VectorXd>& _weight_map, ParameterEnsemble& _pe_upgrade,
                  unordered_map<string, pair<vector<string>, vector<string>>>& _cases, double _reg_factor);

    void work(int thread_id, int iter, double cur_lam, bool use_glm_form, Eigen::VectorXd parcov_inv_vec, Eigen::MatrixXd Am);

protected:
    PerformanceLog* performance_log;
    vector<string> keys;
    int count, total;
	double reg_factor;
    unordered_map<string, pair<vector<string>, vector<string>>>& cases;

    ParameterEnsemble& pe_upgrade;
    unordered_map<string, Eigen::VectorXd>& weight_map;

    unordered_map<string, Eigen::VectorXd>& par_resid_map, & par_diff_map;
    unordered_map<string, Eigen::VectorXd>& obs_resid_map, & obs_diff_map, obs_err_map;

    mutex ctrl_lock, weight_lock, loc_lock, parcov_lock;
    mutex obs_resid_lock, obs_diff_lock, par_resid_lock;
    mutex par_diff_lock, am_lock, put_lock, obs_err_lock;
    mutex next_lock;

};

class UpgradeThread
{
public: 
	UpgradeThread(PerformanceLog* _performance_log, unordered_map<string, Eigen::VectorXd>& _par_resid_map, unordered_map<string, Eigen::VectorXd>& _par_diff_map,
		unordered_map<string, Eigen::VectorXd>& _obs_resid_map, unordered_map<string, Eigen::VectorXd>& _obs_diff_map, unordered_map<string, Eigen::VectorXd>& _obs_err_map,
		Localizer& _localizer, unordered_map<string, double>& _parcov_inv_map,
		unordered_map<string, double>& _weight_map, ParameterEnsemble& _pe_upgrade,
		unordered_map<string, pair<vector<string>, vector<string>>>& _cases,
		unordered_map<string, Eigen::VectorXd>& _Am_map, Localizer::How& _how, double _reg_factor);

	virtual void work(int thread_id, int iter, double cur_lam, bool use_glm_form, vector<string> par_names, vector<string> obs_names) { ; }


    static void ensemble_solution(const int iter, const int verbose_level,const int maxsing,  const int thread_id,
                           const int t_count, const bool use_prior_scaling,const bool use_approx, const bool use_glm,
                           const double cur_lam,const double eigthresh, Eigen::MatrixXd& par_resid, Eigen::MatrixXd& par_diff,
                           const Eigen::MatrixXd& Am, Eigen::MatrixXd& obs_resid,Eigen::MatrixXd& obs_diff, Eigen::MatrixXd& upgrade_1,
                           Eigen::MatrixXd& obs_err, const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& weights,
                           const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& parcov_inv,
                           const vector<string>& act_obs_names,const vector<string>& act_par_names, double _reg_factor);
protected:
	PerformanceLog* performance_log;
	Localizer::How how;
	vector<string> keys;
	int count, total;
	double reg_factor;
	unordered_map<string, pair<vector<string>, vector<string>>>& cases;

	ParameterEnsemble& pe_upgrade;
	Localizer& localizer;
	unordered_map<string, double>& parcov_inv_map;
	unordered_map<string, double>& weight_map;

	unordered_map<string, Eigen::VectorXd>& par_resid_map, & par_diff_map, & Am_map;
	unordered_map<string, Eigen::VectorXd>& obs_resid_map, & obs_diff_map, obs_err_map;

	mutex ctrl_lock, weight_lock, loc_lock, parcov_lock;
	mutex obs_resid_lock, obs_diff_lock, par_resid_lock;
	mutex par_diff_lock, am_lock, put_lock, obs_err_lock;
	mutex next_lock;



};

class LocalAnalysisUpgradeThread: public UpgradeThread
{
public:
	using UpgradeThread::UpgradeThread;

	void work(int thread_id, int iter, double cur_lam,bool use_glm_form, vector<string> par_names, vector<string> obs_names);
};

class EnsembleMethod
{

public:
	EnsembleMethod(Pest& _pest_scenario, FileManager& _file_manager,
		OutputFileWriter& _output_file_writer, PerformanceLog* _performance_log,
		RunManagerAbstract* _run_mgr_ptr, string _alg_tag="EnsembleMethod");

	virtual void throw_em_error(string message);
	bool should_terminate(int current_n_iter_mean=0);
	void sanity_checks();
	//template<typename T, typename A>
	//void message(int level, const string& _message, vector<T, A> _extras, bool echo = true);
	void message(int level, const string& _message, vector<string> _extras, bool echo = true);
	void message(int level, const string& _message, vector<int> _extras, bool echo = true);
	void message(int level, const string& _message, vector<double> _extras, bool echo = true);
	void message(int level, const string& _message);
	//template<typename T>
	//void message(int level, const string& _message, T extra);
	void message(int level, const string& _message, string extra);
	void message(int level, const string& _message, int extra);
	void message(int level, const string& _message, double extra);
	void message(int level, const string& _message, size_t extra);

	ParameterEnsemble get_pe() { return pe; }
	ParameterEnsemble* get_pe_ptr() { return &pe; }
	ObservationEnsemble* get_oe_ptr() {return &oe; }
	void set_pe(ParameterEnsemble& new_pe) { pe = new_pe; }
	void set_oe(ObservationEnsemble& new_oe) { oe = new_oe; }
	void set_noise_oe(ObservationEnsemble& new_noise) { oe_base = new_noise; }
	void set_localizer(Localizer& new_loc) { localizer = new_loc; }
	Localizer get_localizer() { return localizer;  }
	bool initialize_pe(Covariance& cov);
	void initialize_parcov();
	bool initialize_oe(Covariance& cov);
	void initialize_obscov();
	bool initialize_weights();
	Covariance* get_parcov_ptr() { return &parcov; }
	Covariance* get_obscov_ptr() { return &obscov; }
	std::mt19937& get_rand_gen() { return rand_gen; }
	vector<string> get_act_par_names() { return act_par_names; }
	ObservationEnsemble get_oe() { return oe; }
	ObservationEnsemble get_noise_oe() { return oe_base; }
	L2PhiHandler& get_phi_handler() { return ph; }
	int get_iter() { return iter; }
	FileManager& get_file_manager() { return file_manager; }
	Pest* get_pest_scenario_ptr() { return &pest_scenario; }
	
	void initialize(int cycle = NetPackage::NULL_DA_CYCLE, bool run = true, bool use_existing=false);

	//this is not called in the initialization - must be called before initialize() to trigger dynamic state handling...
	void initialize_dynamic_states(bool rec_report=true);

	void transfer_dynamic_state_from_oe_to_initial_pe(ParameterEnsemble& _pe, ObservationEnsemble& _oe);
    void transfer_dynamic_state_from_oe_to_final_pe(ParameterEnsemble& _pe, ObservationEnsemble& _oe);
	void transfer_par_dynamic_state_final_to_initial_ip(ParameterEnsemble& _pe);

	pair<string, string> save_ensembles(string tag, int cycle, ParameterEnsemble& _pe, ObservationEnsemble& _oe);
	vector<string>& get_par_dyn_state_names() { return par_dyn_state_names; }


protected:
	string alg_tag;
	int  verbose_level;
	Pest& pest_scenario;
	FileManager& file_manager;
	std::mt19937 rand_gen;
	std::mt19937 subset_rand_gen;
	OutputFileWriter& output_file_writer;
	PerformanceLog* performance_log;
	RunManagerAbstract* run_mgr_ptr;
	L2PhiHandler ph;
	ParChangeSummarizer pcs;
	Covariance parcov, obscov;
	double reg_factor;
	bool use_localizer;
	Localizer localizer;
	int num_threads;
	set<string> pp_args;
	int iter;
	bool use_subset;
	double last_best_lam, last_best_mean, last_best_std;
	vector<double> best_mean_phis;
	double best_phi_yet;
	vector<double> mda_lambdas;
	vector<string> obs_dyn_state_names, par_dyn_state_names;
	map<string,string> final2init_par_state_names;
	int consec_bad_lambda_cycles;
	double lambda_max, lambda_min;
	int warn_min_reals, error_min_reals;
	vector<double> lam_mults;
	vector<string> oe_org_real_names, pe_org_real_names;
	vector<string> act_obs_names, act_par_names;
	vector<string> violation_obs;
	ParameterEnsemble pe, pe_base;
	ObservationEnsemble oe, oe_base, weights, weights_base;
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> obscov_inv_sqrt, parcov_inv_sqrt;
	bool oe_drawn, pe_drawn;
    bool reinflate_to_minphi_real;
    ObservationInfo org_obs_info;
    string dense_file_ext = ".bin";
	bool solve_glm(int cycle = NetPackage::NULL_DA_CYCLE);

	bool solve_mda(bool last_iter, int cycle = NetPackage::NULL_DA_CYCLE);

	bool solve(bool use_mda, vector<double> inflation_factors, vector<double> backtrack_factors, int cycle=NetPackage::NULL_DA_CYCLE);

	vector<int> run_ensemble(ParameterEnsemble& _pe, ObservationEnsemble& _oe, const vector<int>& real_idxs = vector<int>(), int cycle=NetPackage::NULL_DA_CYCLE);

	vector<ObservationEnsemble> run_lambda_ensembles(vector<ParameterEnsemble>& pe_lams, vector<double>& lam_vals, vector<double>& scale_vals, int cycle, vector<int>& pe_subset_idxs, vector<int>& oe_subset_idxs);

	void report_and_save(int cycle);

	void save_mat(string prefix, Eigen::MatrixXd& mat);

	void initialize_restart();

	void drop_bad_reals(ParameterEnsemble& _pe, ObservationEnsemble& _oe, vector<int> subset_idxs = vector<int>());

	void add_bases();

	void update_reals_by_phi(ParameterEnsemble& _pe, ObservationEnsemble& _oe, vector<int> subset_idxs=vector<int>());

	vector<int> get_subset_idxs(int size, int _subset_size);

	Eigen::MatrixXd get_Am(const vector<string>& real_names, const vector<string>& par_names);


	void zero_weight_obs(vector<string>& obs_to_zero_weight, bool update_obscov = true, bool update_oe_base = true);

	void norm_map_report(map<string, double>& norm_map, string tag, double thres = 0.1);

	void adjust_weights(bool save=false);

    void adjust_weights_single(map<string,vector<string>>& group_to_obs_map, map<string,vector<string>>& group_map,
            map<string,double>& phi_fracs);

    void adjust_weights_by_real(map<string,vector<string>>& group_to_obs_map, map<string,vector<string>>& group_map,
                               map<string,map<string,double>>& phi_fracs_by_real,vector<string> index);

    void check_and_fill_phi_factors(map<string,vector<string>>& group_to_obs_map,map<string,vector<string>>& group_map,
                                    map<string,map<string,double>>& phi_fracs_by_real,
                                    vector<string>& index, bool check_reals);

    void prep_drop_violations();

    void remove_external_pe_filenames(vector<string>& pe_filenames);

    double get_lambda();

    void reinflate_par_ensemble(double reinflate_factor,int reinflate_num_reals);
};
#endif

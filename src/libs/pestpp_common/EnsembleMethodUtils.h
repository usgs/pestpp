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

    void work(int tid, int verbose_level, double mm_alpha, double phi_weight, map<string,map<string,double>> weight_phi_map,
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

/**
 * @brief Inequality observations, and the 'drop_violations' nomination, derived from the scenario.
 *
 * Extracted from L2PhiHandler so that mou and sqp can use it. Neither has a phi handler, but
 * the violation test needs only observed values, simulated values and the inequality treatment
 * - all scenario-level. Keeping ONE definition matters more than the sharing is convenient:
 * the harvest-time drop and the mid-run preemption screen must never disagree about whether a
 * run is violating, or screening cancels a run that harvest would have kept.
 *
 * L2PhiHandler holds one of these and forwards to it, so phi's inequality treatment and the
 * violation test are the same code.
 */
class ViolationDetector
{
public:
	ViolationDetector() : pest_scenario(nullptr) {}
	explicit ViolationDetector(Pest* _pest_scenario);

	/// Zero the residuals of satisfied inequalities, in place. Moved here verbatim from
	/// L2PhiHandler, which now forwards to it.
	void apply_ineq_constraints(Eigen::MatrixXd &resid, Eigen::MatrixXd &sim_vals, vector<string> &names);

	/**
	 * @brief Does ONE run violate a nominated observation?
	 *
	 * Sound on a PARTIAL read because the test is a sum of non-negative terms - absolute
	 * residuals, zero wherever the inequality is satisfied - so the sum over what has been
	 * read is a lower bound on the final sum. Once over the threshold, reading the rest
	 * cannot bring it back under. Anything that makes this a mean, a ratio or a variance
	 * breaks that argument and silently breaks mid-run screening with it.
	 *
	 * @param valid_names which observations in `sim` are REAL; empty means all of them. An
	 *        unread observation carries the no_data sentinel and must never be compared
	 *        against a bound - it would violate almost any inequality.
	 */
	/// `why`, when given, names the observation that condemned the run and its values.
	bool is_violating(Observations& sim, const set<string>& valid_names,
	                  const vector<string>& viol_obs_names, string* why = nullptr);

	vector<string> get_lt_obs_names() { return lt_obs_names; }
	vector<string> get_gt_obs_names() { return gt_obs_names; }
	map<string,double> get_lt_obs_bounds() { return lt_obs_bounds; }
	map<string,double> get_gt_obs_bounds() { return gt_obs_bounds; }
	map<string,pair<double,double>> get_double_obs_bounds() { return double_obs_bounds; }

	/// Observations nominated by the 'drop_violations' column of the external observation
	/// file. Empty means the feature is off, and every violation call is a no-op.
	const vector<string>& get_nominated() const { return nominated; }
	bool has_nominated() const { return !nominated.empty(); }

	/// The threshold every violation test uses.
	static constexpr double VIOLATION_TOL = 1.0e-7;

	/// Read the 'drop_violations' nomination straight from the scenario.
	///
	/// Static because a caller must NOT have to own a detector - or a phi handler - to ask.
	/// Routing this through L2PhiHandler made it depend on when that handler was constructed,
	/// and on the normal path (noptmax != 0) it is built AFTER the nomination is read, so the
	/// answer came back empty and drop_violations silently did nothing.
	static vector<string> read_nominated(Pest& scenario);

	/**
	 * @brief Enforce the preemption precondition. Returns "" if fine, else why not.
	 *
	 * 'preemption_poll_interval_minutes' says how often to ask workers what they have so far
	 * and abandon runs already violating. That is only meaningful if something has been
	 * nominated to judge them against, AND that observation carries a non-zero weight - a
	 * zero-weighted observation is excluded from every violation test, so nominating one and
	 * nothing else would poll the workers forever and never be able to abandon anything.
	 *
	 * Shared so all four tools refuse the same configuration for the same reason, rather than
	 * three of them refusing it and one polling pointlessly.
	 */
	static string check_preemption_config(Pest& scenario);

private:
	Pest* pest_scenario;
	vector<string> lt_obs_names, gt_obs_names, nominated;
	map<string,double> lt_obs_bounds, gt_obs_bounds;
	map<string,pair<double,double>> double_obs_bounds;
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
	vector<string> get_lt_obs_names() { return ineq.get_lt_obs_names(); }
	vector<string> get_gt_obs_names() { return ineq.get_gt_obs_names(); }
    map<string,double> get_lt_obs_bounds() {return ineq.get_lt_obs_bounds();}
    map<string,double> get_gt_obs_bounds() {return ineq.get_gt_obs_bounds();}
    map<string,pair<double,double>> get_double_obs_bounds() {return ineq.get_double_obs_bounds();}
    /// The shared evaluator; mou and sqp use one directly.
    ViolationDetector& get_ineq() { return ineq; }

	void apply_ineq_constraints(Eigen::MatrixXd &resid, Eigen::MatrixXd &sim_vals, vector<string> &names);

	void save_residual_cov(ObservationEnsemble& oe, int iter);

	//map<string,double> get_meas_phi(ObservationEnsemble& oe, Eigen::VectorXd& q_vec);

	map<string,map<string,double>> get_swr_real_map(ObservationEnsemble& oe, ObservationEnsemble& weights,phiType ptype=phiType::MEAS);

    map<string,double> get_swr_map(ObservationEnsemble& oe, string real_name= "",phiType ptype=phiType::MEAS);
	map<string,map<string,double>> get_meas_phi_weight_ensemble(ObservationEnsemble& oe, ObservationEnsemble& weights);

    vector<string> get_violating_realizations(ObservationEnsemble& oe, const vector<string>& viol_obs_names);
    /**
     * @brief Does ONE run already violate a nominated 'drop_violations' observation?
     *
     * The single-run twin of get_violating_realizations(), for judging a run from a PARTIAL
     * read while it is still executing (docs/api_part1/partial_screening_plan.md). Both go
     * through the same residual computation and the same threshold, because a screen that
     * disagreed with the harvest-time drop at the margin would cancel a run harvest would
     * have kept.
     *
     * Sound on partial data because the test is a SUM OF NON-NEGATIVE TERMS - absolute
     * residuals, zero wherever the inequality is satisfied - so the sum over what has been
     * read is a lower bound on the final sum. Once it is over the threshold, reading the rest
     * cannot bring it back under. Anything that later makes this a mean, a ratio or a
     * variance breaks that argument and silently breaks mid-run screening with it.
     *
     * @param valid_names which observations in `sim` are REAL. Empty means all of them, for a
     *        complete run. An unread observation carries the no_data sentinel and must never
     *        be compared against a bound - it would violate almost anything.
     */
    /// `why`, when given, names the observation that condemned the run and its values.
    bool is_violating(Observations& sim, const set<string>& valid_names,
                      const vector<string>& viol_obs_names, string* why = nullptr);
    /// The threshold both violation tests use.
    static constexpr double VIOLATION_TOL = 1.0e-7;
private:
    /// Nominated observations that actually count: present in the ensemble's variables and
    /// non-zero weighted, as name -> column index. Shared so the two tests cannot disagree
    /// about WHICH observations are being judged.
    map<string,int> get_violation_idx_map(const vector<string>& viol_obs_names,
                                          const vector<string>& act_obs_names);
public:
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

	// live reg factor for phi: negative option value means 'full solution', for which the
	// phi handler must ignore regularization (0.0); otherwise the option value. Was cached in
	// org_reg_factor after the option got mutated to 0.0 at init.
	double get_reg_factor() const { double r = pest_scenario->get_pestpp_options().get_ies_reg_factor(); return r < 0.0 ? 0.0 : r; }
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

	/// the shared inequality/violation evaluator - replaces the sets this class used to own
	ViolationDetector ineq;

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

// queue -> (drive the run manager) -> harvest, for a single ensemble. run_ensemble_util() is
// the in-tree composition; the halves are available for a caller running its own run loop.
// the returned map is keyed by realization NAME, so the caller may change ensemble
// membership between these two calls without the runs being misattributed
map<string, int> queue_ensemble_util(PerformanceLog* performance_log, ofstream& frec, ParameterEnsemble& _pe,
	RunManagerAbstract* run_mgr_ptr,
	bool check_pe_consistency = false, const vector<int>& real_idxs = vector<int>(), int da_cycle=NetPackage::NULL_DA_CYCLE,
	string additional_tag="");

vector<int> harvest_ensemble_util(PerformanceLog* performance_log, ofstream& frec, ParameterEnsemble& _pe,
	ObservationEnsemble& _oe, RunManagerAbstract* run_mgr_ptr,
	bool check_pe_consistency, const vector<int>& real_idxs, map<string, int>& real_run_ids);

vector<int> run_ensemble_util(PerformanceLog* performance_log, ofstream& frec, ParameterEnsemble& _pe,
	ObservationEnsemble& _oe, RunManagerAbstract* run_mgr_ptr,
	bool check_pe_consistency = false, const vector<int>& real_idxs = vector<int>(),int da_cycle=NetPackage::NULL_DA_CYCLE,
	string additional_tag="");

class EnsembleSolver
{
public:
	EnsembleSolver(PerformanceLog* _performance_log, FileManager& _file_manager, Pest& _pest_scenario, ParameterEnsemble& _pe,
		ObservationEnsemble& _oe, ObservationEnsemble& _base_oe, ObservationEnsemble& _weights, Localizer& _localizer,
		Covariance& _parcov,Eigen::MatrixXd& _am, L2PhiHandler& _ph,
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
    //per-center-realization weights (aligned to mm_real_name_map[key].first) used to scale
    //realization contributions when mm_alpha >= 1 and ies_multimodal_weight_exponent > 0; empty otherwise
    unordered_map<string,Eigen::VectorXd> mm_real_weight_map;
    //the weight vector for the center realization currently being solved in the localized
    //multimodal path; passed to LocalAnalysisUpgradeThread via solve(); empty when not weighting
    Eigen::VectorXd mm_current_weights;
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
                  unordered_map<string, pair<vector<string>, vector<string>>>& _cases, double _reg_factor,
                  unordered_map<string, Eigen::VectorXd>& _real_weight_map);

    void work(int thread_id, int iter, double cur_lam, bool use_glm_form, Eigen::VectorXd parcov_inv_vec, Eigen::MatrixXd Am);

protected:
    PerformanceLog* performance_log;
    vector<string> keys;
    int count, total;
	double reg_factor;
    unordered_map<string, pair<vector<string>, vector<string>>>& cases;

    ParameterEnsemble& pe_upgrade;
    unordered_map<string, Eigen::VectorXd>& weight_map;
    unordered_map<string, Eigen::VectorXd>& real_weight_map;

    unordered_map<string, Eigen::VectorXd>& par_resid_map, & par_diff_map;
    unordered_map<string, Eigen::VectorXd>& obs_resid_map, & obs_diff_map, & obs_err_map;

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

	//optional per-realization weights (aligned to the realization columns) for the localized
	//multimodal solve; empty means no realization weighting
	void set_real_weights(const Eigen::VectorXd& w) { real_weights = w; }


    static void ensemble_solution(const int iter, const int verbose_level,const int maxsing,  const int thread_id,
                           const int t_count, const bool use_prior_scaling,const bool use_approx, const bool use_glm,
                           const double cur_lam,const double eigthresh, Eigen::MatrixXd& par_resid, Eigen::MatrixXd& par_diff,
                           const Eigen::MatrixXd& Am, Eigen::MatrixXd& obs_resid,Eigen::MatrixXd& obs_diff, Eigen::MatrixXd& upgrade_1,
                           Eigen::MatrixXd& obs_err, const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& weights,
                           const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& parcov_inv,
                           const vector<string>& act_obs_names,const vector<string>& act_par_names, double _reg_factor,
                           double mm_weight_sum = -1.0);
protected:
	PerformanceLog* performance_log;
	Localizer::How how;
	vector<string> keys;
	int count, total;
	double reg_factor;
	unordered_map<string, pair<vector<string>, vector<string>>>& cases;

	ParameterEnsemble& pe_upgrade;
	Localizer& localizer;
	Eigen::VectorXd real_weights;
	unordered_map<string, double>& parcov_inv_map;
	unordered_map<string, double>& weight_map;

	unordered_map<string, Eigen::VectorXd>& par_resid_map, & par_diff_map, & Am_map;
	unordered_map<string, Eigen::VectorXd>& obs_resid_map, & obs_diff_map, &obs_err_map;

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

/**
 * @brief The par-ensemble reinflation schedule: when to reinflate, by how much, over how many
 *        realizations.
 *
 * ies_n_iter_reinflate / ies_reinflate_factor / ies_reinflate_num_reals are parallel vectors
 * walked one entry at a time as reinflations are consumed. That bookkeeping used to live as
 * six loose locals in IterEnsembleSmoother::iterate_2_solution(), which meant any caller
 * writing its own loop had to re-derive it - and re-derive the index quirks with it. Holding
 * it here keeps the built-in loop and an API caller on identical rules.
 *
 * Usage per solution iteration: tick(), then if due() reinflate with get_factor() /
 * get_num_reals() and call advance().
 */
class ReinflationSchedule
{
public:
	ReinflationSchedule(Pest& pest_scenario);

	/// count one completed solution iteration
	void tick() { iters_since++; }
	/// is a reinflation due now? (never when the current entry is 0 == reinflation off)
	bool due() const { return (current_n_iter != 0) && (iters_since >= current_n_iter); }
	/// consume this reinflation and move to the next schedule entry
	void advance();

	/// iterations between reinflations for the current entry; 0 means reinflation is off
	int get_n_iter() const { return current_n_iter; }
	double get_factor() const { return current_factor; }
	int get_num_reals() const { return current_num_reals; }
	/// is reinflation in use at this point in the schedule?
	bool is_active() const { return current_n_iter != 0; }

private:
	vector<int> n_iter_reinflate;
	vector<double> reinflate_factor;
	vector<int> reinflate_num_reals;
	int iters_since = 0;
	int idx = 0;
	int current_n_iter = 0;
	double current_factor = 1.0;
	int current_num_reals = 0;
};

/**
 * @brief Outcome of one stage of an upgrade solve.
 *
 * The old solve() signalled all of this through a bool, where `false` did not mean "error"
 * but "retry this iteration with a different lambda". Naming the outcomes lets the stages be
 * driven one at a time - by solve(), or by an API caller running its own loop.
 */
enum class UpgradeStatus
{
	CONTINUE,        ///< stage finished, go on to the next one
	ACCEPTED,        ///< the upgrade was taken; the iteration is done
	REJECTED_RETRY,  ///< no usable upgrade; lambda was adjusted, run the iteration again
	HALTED           ///< a 'pest.stp' quit file was seen; stop without further work
};

/**
 * @brief State carried across the stages of a single upgrade solve.
 *
 * These were locals inside the 688-line solve(). Holding them in one place is what lets the
 * stages be called separately: generate the upgrade ensembles, hand them to the run manager,
 * then evaluate what came back.
 *
 * Deliberately non-copyable: `pe_lams`/`oe_lams` are large, and stages take it by reference.
 */
struct UpgradeContext
{
	UpgradeContext(Pest* _pest_scenario) : oe_lam_best(_pest_scenario) {}
	UpgradeContext(const UpgradeContext&) = delete;
	UpgradeContext& operator=(const UpgradeContext&) = delete;

	// -- the factors this solve is running with, filled by solve_prepare(). Held here so the
	//    run and evaluate halves cannot drift onto a different set than the candidates were
	//    generated with, which is a real hazard once a caller can change options mid-solve.
	bool use_mda = false;
	vector<double> inflation_factors, backtrack_factors;

	// -- filled by prepare_upgrades()
	int local_subset_size = 0;
	/// Realizations chosen for subset testing, held by NAME rather than by position.
	/// Positions go stale the moment membership changes - which is exactly what an API
	/// caller does when it culls realizations mid-iteration - so the names are the state
	/// and positions are resolved fresh at each point of use.
	vector<string> subset_names;
	unordered_map<string, pair<vector<string>, vector<string>>> loc_map;
	Eigen::MatrixXd Am;

	// -- filled by generate_upgrades(): one entry per lambda x scale-factor combination
	vector<ParameterEnsemble> pe_lams;
	vector<double> lam_vals, scale_vals;
	vector<string> pe_filenames;   ///< non-empty only when ies_upgrades_in_memory is false

	// -- filled by queue_upgrade_ensembles(), reused by harvest_upgrade_ensembles() so both
	//    halves address the same rows even if membership moved in between
	vector<int> pe_subset_idxs, oe_subset_idxs;

	// -- filled by run_upgrade_ensembles()
	vector<ObservationEnsemble> oe_lams;

	// -- filled by evaluate_upgrades()
	ObservationEnsemble oe_lam_best;
	int best_idx = -1;
	double best_mean = 1.0e+300, best_std = 1.0e+300;

	// -- filled by prepare_subset_completion(), consumed by finish_subset_completion()
	//
	// The remaining-realization runs. Subset testing evaluates a few realizations against
	// every lambda; once a winner is picked the REST of the ensemble still has to be run at
	// that lambda. That run sat in the middle of complete_subset_runs() with a dozen locals
	// either side of it, which is precisely what stopped a caller from owning it.
	/// false when subset testing was not in play and there is nothing left to run
	bool needs_remaining_runs = false;
	ParameterEnsemble remaining_pe;
	ObservationEnsemble remaining_oe;
	/// the realization names remaining_pe/remaining_oe went in with, so failures can be
	/// identified positionally after the fact the way the in-tree path always has
	vector<string> org_pe_names, org_oe_names;
	/// the non-subset members, i.e. who the remaining runs are for
	vector<string> pe_keep_names, oe_keep_names;
	/// positions within remaining_pe that failed; filled by whoever performed the runs
	vector<int> failed_remaining;

	/**
	 * Who decides when the losing candidates die.
	 *
	 * Left false, they are reclaimed inline at the point solve() has always done it -
	 * immediately before the remaining-realization runs, which is where memory peaks. Set
	 * true and nothing is released automatically: the caller owns the lifetime and calls
	 * release_unused_candidates() / release_spilled_candidate_files() itself. That is what
	 * lets an API caller keep every candidate and response around to compare, and it also
	 * means a borrowed view into pe_lams stays valid until the caller says otherwise.
	 *
	 * Beware the asymmetry: unreleased memory comes back when the context dies, but
	 * unreleased spill files do not - those filenames are iteration-stamped and accumulate.
	 */
	bool defer_candidate_release = false;
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
	void message(int level, const string& _message, string extra, bool echo = true);
	void message(int level, const string& _message, int extra);
	void message(int level, const string& _message, double extra);
	void message(int level, const string& _message, size_t extra);

	ParameterEnsemble get_pe() { return pe; }
	ParameterEnsemble* get_pe_ptr() { return &pe; }
	ObservationEnsemble* get_oe_ptr() {return &oe; }
	/// Pointer accessors for the other two coupled ensembles, matching get_pe_ptr()/
	/// get_oe_ptr(). get_noise_oe() returns by value, which is a copy an API caller does
	/// not want for a large ensemble.
	ObservationEnsemble* get_noise_oe_ptr() { return &oe_base; }
	ObservationEnsemble* get_weights_ptr() { return &weights; }
	/// The PRIOR parameter ensemble. Reinflation draws its realizations from here, so this is
	/// also the hard ceiling on how many realizations a reinflation can produce - a caller
	/// needs to be able to see that number before asking for one.
	ParameterEnsemble* get_pe_base_ptr() { return &pe_base; }
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
	
	// initialize() splits at the prior-ensemble evaluation so a caller can own that batch -
	// and, more usefully, substitute its own prior ensemble before it is ever run.
	// initialize() is the in-tree composition of the two halves.
	int  initialize_prepare(int cycle = NetPackage::NULL_DA_CYCLE, bool run = true, bool use_existing=false);
	void initialize_finish(int cycle = NetPackage::NULL_DA_CYCLE);
	void initialize(int cycle = NetPackage::NULL_DA_CYCLE, bool run = true, bool use_existing=false);

	//this is not called in the initialization - must be called before initialize() to trigger dynamic state handling...
	void initialize_dynamic_states(bool rec_report=true);

	void transfer_dynamic_state_from_oe_to_initial_pe(ParameterEnsemble& _pe, ObservationEnsemble& _oe);
    void transfer_dynamic_state_from_oe_to_final_pe(ParameterEnsemble& _pe, ObservationEnsemble& _oe);
	void transfer_par_dynamic_state_final_to_initial_ip(ParameterEnsemble& _pe);

	pair<string, string> save_ensembles(string tag, int cycle, ParameterEnsemble& _pe, ObservationEnsemble& _oe);
	vector<string>& get_par_dyn_state_names() { return par_dyn_state_names; }

	// ---- per-iteration building blocks --------------------------------------------------
	// The loop in IterEnsembleSmoother::iterate_2_solution() and DataAssimilator::da_update()
	// is written against nothing but these, so if the built-in loop compiles, an API caller
	// can write any loop. solve()/solve_glm()/solve_mda() are the in-tree compositions of the
	// stages below and double as the worked example of how to sequence them.
	/// The inflation (lambda / mda) and backtrack (scale) factors for the current iteration -
	/// the prologues of solve_glm()/solve_mda(), so a caller driving the stages gets the same
	/// numbers. get_mda_factors() advances the mda schedule and must be called once per
	/// iteration, which is why solve_prepare() is the only thing that should call it.
	void get_glm_factors(vector<double>& inflation_factors, vector<double>& backtrack_factors);
	void get_mda_factors(bool last_iter, vector<double>& inflation_factors, vector<double>& backtrack_factors);
	/// solve() up to but NOT including the upgrade runs: factors, prepare, generate, shortcut
	/// check. Leaves the candidates in ctx.pe_lams for a caller to inspect or replace before
	/// they are run. CONTINUE means candidates are waiting; anything else ended the iteration.
	UpgradeStatus solve_prepare(UpgradeContext& ctx, bool use_mda, bool last_iter = false);
	UpgradeStatus solve_glm(int cycle = NetPackage::NULL_DA_CYCLE);

	UpgradeStatus solve_mda(bool last_iter, int cycle = NetPackage::NULL_DA_CYCLE);

	UpgradeStatus solve(bool use_mda, vector<double> inflation_factors, vector<double> backtrack_factors, int cycle=NetPackage::NULL_DA_CYCLE);
	// The stages solve() is built from, in call order. Each takes the shared UpgradeContext
	// by reference; a status other than CONTINUE means the iteration is over. Split out so an
	// API caller can drive generate -> run -> evaluate itself and own run management in
	// between; solve() is now just the in-tree composition of these.
	UpgradeStatus prepare_upgrades(UpgradeContext& ctx, bool use_mda,
		const vector<double>& inflation_factors, const vector<double>& backtrack_factors);
	void generate_upgrades(UpgradeContext& ctx, bool use_mda,
		const vector<double>& inflation_factors, const vector<double>& backtrack_factors);
	UpgradeStatus check_noniterative_shortcut(UpgradeContext& ctx);
	/// The two halves of run_upgrade_ensembles(), for a caller that wants to own the candidate
	/// batch. They live here rather than in the caller because resolving the subset against
	/// current membership - and the spill path's different indexing - are the tool's business.
	vector<map<string, int>> queue_upgrade_ensembles(UpgradeContext& ctx, int cycle);
	void harvest_upgrade_ensembles(UpgradeContext& ctx, vector<map<string, int>>& run_ids);
	/// Likewise for the remaining-realization batch prepare_subset_completion() sets up. Here
	/// rather than in the caller because the run manager and performance log are the tool's.
	map<string, int> queue_remaining_runs(UpgradeContext& ctx, int cycle);
	void harvest_remaining_runs(UpgradeContext& ctx, map<string, int>& run_ids);
	void run_upgrade_ensembles(UpgradeContext& ctx, int cycle,
		const vector<double>& inflation_factors, const vector<double>& backtrack_factors);
	UpgradeStatus evaluate_upgrades(UpgradeContext& ctx);
	UpgradeStatus complete_subset_runs(UpgradeContext& ctx, int cycle, bool use_mda);
	/// The two halves of complete_subset_runs(), either side of the remaining-realization
	/// runs. prepare_ leaves ctx.remaining_pe/remaining_oe ready and ctx.needs_remaining_runs
	/// set; whoever runs them records failures in ctx.failed_remaining; finish_ assembles the
	/// result. complete_subset_runs() is the in-tree composition and runs them itself.
	UpgradeStatus prepare_subset_completion(UpgradeContext& ctx, int cycle, bool use_mda);
	UpgradeStatus finish_subset_completion(UpgradeContext& ctx, int cycle, bool use_mda);
	UpgradeStatus accept_or_reject(UpgradeContext& ctx, bool use_mda, int cycle);

	// Candidate lifetime. Called inline by the stages above unless the context asks to defer
	// (UpgradeContext::defer_candidate_release), in which case they are the caller's to call.
	void release_unused_candidates(UpgradeContext& ctx, bool include_responses = false);
	void release_spilled_candidate_files(UpgradeContext& ctx);

	vector<int> run_ensemble(ParameterEnsemble& _pe, ObservationEnsemble& _oe, const vector<int>& real_idxs = vector<int>(), int cycle=NetPackage::NULL_DA_CYCLE);

	// queue -> (drive the run manager) -> harvest. run_lambda_ensembles() is the in-tree
	// composition of the two halves; a caller that wants to watch or cancel runs while they
	// are in flight calls the halves itself with its own run_slice() loop in between.
	vector<map<string, int>> queue_lambda_ensembles(vector<ParameterEnsemble>& pe_lams, vector<double>& lam_vals, vector<double>& scale_vals, int cycle, vector<int>& pe_subset_idxs, vector<int>& oe_subset_idxs);
	vector<ObservationEnsemble> harvest_lambda_ensembles(vector<ParameterEnsemble>& pe_lams, vector<double>& lam_vals, vector<double>& scale_vals, vector<map<string, int>>& real_run_ids_vec, vector<int>& pe_subset_idxs, vector<int>& oe_subset_idxs);
	vector<ObservationEnsemble> run_lambda_ensembles(vector<ParameterEnsemble>& pe_lams, vector<double>& lam_vals, vector<double>& scale_vals, int cycle, vector<int>& pe_subset_idxs, vector<int>& oe_subset_idxs);

	void report_and_save(int cycle);
	/// One iteration's bookkeeping, either side of the solve. The shipped loop and an API
	/// caller both go through these, so iteration numbering and per-iteration output match.
	void begin_iteration();
	void end_iteration(int cycle = NetPackage::NULL_DA_CYCLE);
	void adjust_weights(bool save=false);
	/// Rebuild the parameter ensemble from the prior's spread, re-centred on where the current
	/// ensemble has got to.
	///
	/// reinflate_num_reals SELECTS realizations from pe_base; it never draws new ones, so
	/// abs() of it cannot exceed pe_base's row count. Its SIGN picks where the anomalies come
	/// from: positive uses pe_base's own rows scaled by the factor, negative resamples the
	/// CURRENT ensemble's anomalies (adding scaled prior anomalies when the factor is < 1).
	///
	/// center_on_min_phi is a TRI-STATE: negative means "whatever ies_n_iter_reinflate says",
	/// which is what the shipped loop passes and leaves its behaviour unchanged; 0 forces the
	/// ensemble mean; positive forces the min-phi realization. The option version of this is a
	/// schedule-wide flag derived from the SIGN of an unrelated setting, which a caller driving
	/// one reinflation at a time has no way to express - hence the argument.
    void reinflate_par_ensemble(double reinflate_factor,int reinflate_num_reals,
                                int center_on_min_phi = -1);

protected:
	string alg_tag;
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
	// live reg factor magnitude for the upgrade calc (abs of the option; a negative option
	// value signals 'full solution' but still uses the magnitude). Was a cached member.
	double get_reg_factor() const { double r = pest_scenario.get_pestpp_options().get_ies_reg_factor(); return r < 0.0 ? -r : r; }
	// live verbosity - was cached at initialize(), so bumping it mid-run did nothing
	int get_verbose_level() const { return pest_scenario.get_pestpp_options().get_ies_verbose_level(); }
	// live thread count - each solve() spins its own pool, so this is safe to change per iteration
	int get_num_threads() const { return pest_scenario.get_pestpp_options().get_ies_num_threads(); }
	// derived: subset testing is off when the requested subset is at least the whole ensemble.
	// negative (percentage) subset sizes always compare below the ensemble size, so they enable it
	bool get_use_subset() { return pest_scenario.get_pestpp_options().get_ies_subset_size() <= pe.shape().first; }
	// derived: any negative entry in n_iter_reinflate selects the min-phi realization to reinflate around
	bool get_reinflate_to_minphi_real() const {
		for (auto n : pest_scenario.get_pestpp_options().get_ies_n_iter_reinflate())
			if (n < 0) return true;
		return false;
	}
	bool use_localizer;
	Localizer localizer;
	set<string> pp_args;
	int iter;
	double last_best_lam, last_best_mean, last_best_std;
	vector<double> best_mean_phis;
	double best_phi_yet;
	vector<double> mda_lambdas;
	vector<string> obs_dyn_state_names, par_dyn_state_names;
	map<string,string> final2init_par_state_names;
	int consec_bad_lambda_cycles;
	double lambda_max, lambda_min;
	int warn_min_reals, error_min_reals;
	vector<string> oe_org_real_names, pe_org_real_names;
	vector<string> act_obs_names, act_par_names;
	// carried from initialize_prepare() to initialize_finish(); meaningless outside that window
	vector<string> init_pnames, init_onames;
	bool init_needs_finish = false;
	vector<string> violation_obs;
	/// This class's OWN evaluator, not ph's. The phi handler is constructed at several points
	/// and is still default-constructed during the initial ensemble run - reaching through it
	/// from the screener closure dereferenced a null scenario and segfaulted the master. mou
	/// and sqp each hold their own for the same reason.
	ViolationDetector viol_detector;
	ParameterEnsemble pe, pe_base;
	ObservationEnsemble oe, oe_base, weights, weights_base;
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> obscov_inv_sqrt, parcov_inv_sqrt;
	bool oe_drawn, pe_drawn;
    ObservationInfo org_obs_info;
    // not an option cache: this records a one-time resolved decision (the save_binary +
    // oversized-problem fallback flips save_dense on, erasing the state it was derived from)
    string dense_file_ext = ".bin";


	void save_mat(string prefix, Eigen::MatrixXd& mat);

	void initialize_restart();

	void drop_bad_reals(ParameterEnsemble& _pe, ObservationEnsemble& _oe, vector<int> subset_idxs = vector<int>());

	void add_bases();

	void update_reals_by_phi(ParameterEnsemble& _pe, ObservationEnsemble& _oe, vector<int> subset_idxs=vector<int>());

	vector<int> get_subset_idxs(int size, int _subset_size);
	/// Map subset realization names onto their current positions in `current_names`.
	/// Names no longer present are silently skipped, so a caller that dropped realizations
	/// gets a smaller subset rather than an out-of-range index.
	vector<int> resolve_subset_idxs(const vector<string>& names, const vector<string>& current_names) const;

	Eigen::MatrixXd get_Am(const vector<string>& real_names, const vector<string>& par_names);


	void zero_weight_obs(vector<string>& obs_to_zero_weight, bool update_obscov = true, bool update_oe_base = true);
public:
	/// The mirror of zero_weight_obs(): bring observations that are currently zero-weighted
	/// back into the active set, generating the noise realizations they never had. Returns
	/// the names actually activated (those whose weight really was zero and is now not).
	/// Observations already carrying a weight are reweighted and nothing structural happens,
	/// which is the common case and deliberately cheap.
	vector<string> activate_obs(const map<string, double>& obs_to_activate);
protected:

	void norm_map_report(map<string, double>& norm_map, string tag, double thres = 0.1);


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

};
#endif

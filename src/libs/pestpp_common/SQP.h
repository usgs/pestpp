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
#include "Localizer.h"
#include "EnsembleMethodUtils.h"
#include "constraints.h"




struct FilterRec
{
	double obj_val;
	double viol_val;
	int iter;
	Parameters dp_val;
	Observations oe_val;
	string real_name;
	double viol_padded;
    // Strict weak ordering with a deterministic total tie-break. The previous
    // comparator (logical AND of two '<') encoded Pareto dominance, which is NOT
    // a strict weak ordering: incomparable records produced a non-transitive
    // "equivalence", which is undefined behavior for the std::multiset key and
    // ordered them differently on libstdc++ vs MSVC. Dominance/pruning is handled
    // separately (first_strictly_dominates_second); this only defines storage order.
    friend bool operator<(const FilterRec &k1, const FilterRec &k2) {
        if (k1.obj_val != k2.obj_val)
            return k1.obj_val < k2.obj_val;
        if (k1.viol_val != k2.viol_val)
            return k1.viol_val < k2.viol_val;
        return k1.real_name < k2.real_name;
    }
};

class SqpFilter
{
public:
	SqpFilter(bool _minimize=true,double _obj_tol = 0.001, double _viol_tol = 0.001) {
		minimize = _minimize; obj_tol = _obj_tol; viol_tol = _viol_tol;
	}
	bool accept(double obj_val, double violation_val, double violation_padded, Parameters p, Observations o, string rname, int iter=0, bool keep=false);
	bool update(double obj_val, double violation_val, double violation_padded, Parameters p, Observations o, string rname, int iter=0);
    void report(ofstream& frec,int iter);
    double get_viol_tol() {return viol_tol;}
	void set_tol(double tol) { 
		obj_tol = tol; 
		viol_tol = tol;}
	vector<FilterRec> get_feasible_solutions(bool padded = false) const;
	vector<FilterRec> get_filter_members() const;

private:
	bool minimize;
	double obj_tol;
	double viol_tol;
	multiset<FilterRec> obj_viol_pairs;
	bool first_partially_dominates_second(const FilterRec& first, const FilterRec& second);
    bool first_strictly_dominates_second(const FilterRec& first, const FilterRec& second);

};

class CovMatAdap
{

public:
	CovMatAdap(Pest* pest_ptr, std::mt19937* rand_gen, FileManager* file_mgr) : pest_scenario_ptr(pest_ptr), rand_gen_ptr(rand_gen), file_manager(file_mgr) {}
	CovMatAdap() : pest_scenario_ptr(nullptr), rand_gen_ptr(nullptr), file_manager(nullptr) {}

	void initialize(int n_params, int _num_reals);
	void update(Parameters prev_m, Parameters curr_m, int iter);
	void reinflate_C(double reinflation_factor = 1.0, bool reset_corr = false, double max_cond_num = -1.0);
	void set_covariance(const Eigen::MatrixXd& _C) { C = _C; }
	Eigen::MatrixXd get_covariance_matrix() { return C; }

	void set_mean(const Eigen::VectorXd& _m) { m = _m; }
	Eigen::VectorXd get_mean() const { return m; }
	double get_sigma() const { return sigma; }
	int get_parent_num() { return mu; }
	void set_parent_num(int mu_new) { mu = mu_new; }

	ParameterEnsemble generate_population(Parameters& _curr_m, ParameterEnsemble _dv);
	// Diagnostics for the last generate_population() call: total sample draws (a large
	// value relative to the population size signals bounds-rejection thrashing, which is
	// O(n_dv^2) per draw and degrades badly in high dimensions) and how many realizations
	// exhausted max_draws and had to be clipped to bounds.
	long get_last_generate_draws() const { return last_generate_draws; }
	int get_last_generate_clipped() const { return last_generate_clipped; }
	void update_archives(const ParameterEnsemble& pe, map<string, double> obj_map, map<string, double> viol_map, string tag, bool clear = true);
	void clear_archives();
	string get_cma_update_summary() const { return cma_update_summary; }
	bool should_terminate();

private:
	struct CovMetrics {
		double trace;
		double determinant;
		double frobenius_norm;
		double max_eigenvalue;
		double min_eigenvalue;
		double condition_number;
	};

	int lambda;                    // Population size
	int mu;                       // Number of parents
	double sigma = 1.0;           // Step size
	Eigen::VectorXd m;            // Mean vector
	Eigen::MatrixXd C;            // Covariance matrix
	Eigen::MatrixXd B;            // Eigenvectors of C
	Eigen::VectorXd D;            // Eigenvalues of C
	Eigen::MatrixXd pc;           // Evolution path
	Eigen::MatrixXd ps;           // Evolution path for sigma
	double c_sigma, c_c, c_1, c_mu, d_sigma, chi_n, mu_eff;
	vector<double> weights;

	double trace_ratio, det_ratio, frobenius_ratio, max_eigenval_ratio;
	double trace_ratio_0, det_ratio_0, frobenius_ratio_0, max_eigenval_ratio_0;
	bool possibly_converged = false;
	long last_generate_draws = 0;   // total sample draws in the last generate_population() call
	int last_generate_clipped = 0;  // realizations that hit max_draws and were clipped to bounds
	
	ParameterEnsemble sorted_dp_archive;
	map<string, double> sorted_obj_map;
	string cma_update_summary;
	CovMetrics metrics_init;
	CovMetrics compute_cov_metrics() const;
	string report_cma_metrics(const CovMetrics& prior, const CovMetrics& post, int iter);
	
	

protected:
	std::mt19937* rand_gen_ptr;
	Pest* pest_scenario_ptr;
	FileManager* file_manager;
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

	// ---- per-iteration building blocks -------------------------------------------------
	// SQP does not decompose into one generate -> run -> evaluate per iteration the way ies
	// and mou do: the line search and trust-region step each propose a step, run it, judge
	// it, and may try again, so a single iteration issues several run batches. These are the
	// units that shape actually has. solve_new_ensemble() is the in-tree composition of them.

	/// One iteration: search direction, then line search / trust region until a step is taken.
	bool solve_new_ensemble();

	/// The 'generate' half - the search direction for the current point.
	pair<Eigen::VectorXd, Eigen::VectorXd> calc_search_direction_vector(Parameters& _current_dv_, Observations& _current_obs_values, Eigen::VectorXd& grad_vector, Eigen::MatrixXd* _constraint_jco ,vector<string>* _cnames = nullptr);
	bool recalc_search_direction_vector(const string& realization, Parameters& dv_vals, Observations& obs_vals, Eigen::VectorXd& grad_vector);

	/// Propose-run-judge cycles. Each runs candidates internally via
	/// queue_candidate_ensemble()/harvest_candidate_ensemble().
	FilterRec line_search(map<string, Eigen::VectorXd>& search_d_map, Eigen::VectorXd& grad, map<string, double> current_obj_ens, ParameterEnsemble* dvs_subset = nullptr, bool recalc = false);
	bool trust_region_step(Parameters& current_dv_values, Eigen::VectorXd grad);
	FilterRec trust_region_step(Eigen::VectorXd& grad, map<string, double> current_obj_ens, map<string, vector<string>>& cnames_en,
		map<string, Eigen::MatrixXd>& constraint_jco_en, ParameterEnsemble* dvs_subset, bool recalc);

	// queue -> (drive the run manager) -> harvest. The run_* calls are the in-tree
	// compositions; the halves let a caller run its own run_slice() loop in between.
	map<int, int> queue_ensemble(ParameterEnsemble &_pe, const vector<int> &real_idxs=vector<int>());
	vector<int> harvest_ensemble(ParameterEnsemble &_pe, ObservationEnsemble &_oe, const vector<int> &real_idxs, map<int, int>& real_run_ids);
	vector<int> run_ensemble(ParameterEnsemble &_pe, ObservationEnsemble &_oe, const vector<int> &real_idxs=vector<int>());

	map<int, int> queue_candidate_ensemble(ParameterEnsemble& dv_candidates);
	ObservationEnsemble harvest_candidate_ensemble(ParameterEnsemble& dv_candidates, map<int, int>& real_run_ids);
	ObservationEnsemble run_candidate_ensemble(ParameterEnsemble&dv_candidates);

protected:
	// Derived live from the options.  Protected rather than private only so the selftest
	// harness can subclass and assert they track a post-initialize() option change; this is
	// not part of the intended public surface.

	// sqp_subset_size follows the ies convention: 0 => no subset (use the full ensemble);
	// any non-zero value (negative = percentage of the ensemble, positive = absolute count)
	// requests a subset.  Read live so a mid-run change to sqp_subset_size propagates.
	bool get_use_subset() const { return pest_scenario.get_pestpp_options().get_sqp_subset_size() != 0; }

	// an ensemble gradient is used whenever an ensemble was requested, either by size or by
	// file.  both drivers are init-only, so this is stable across a run
	bool get_use_ensemble_grad() const {
		return (pest_scenario.get_pestpp_options().get_sqp_num_reals() > 0) ||
		       (pest_scenario.get_pestpp_options().get_sqp_dv_en().size() > 0);
	}

private:
	int  verbose_level = 1;   // default until initialize() syncs it from ies_verbose_level
	Pest &pest_scenario;
	FileManager &file_manager;
	std::mt19937 rand_gen;
	std::mt19937 subset_rand_gen;
	OutputFileWriter &output_file_writer;
	PerformanceLog *performance_log;
	RunManagerAbstract* run_mgr_ptr;

	ParChangeSummarizer pcs;
	Covariance parcov, obscov, uncertain_parcov;
	CovMatAdap cma;
	chancePoints chancepoints;
	string obj_func_str;
	string obj_obs;
	string obj_sense;
	bool use_obj_obs;
	bool use_obj_pi;
	bool converged = false;
	bool found_feasible = false;
	map<string, double> obj_func_coef_map;
	bool reset = false;
	int recalc_attempt = 0;
	int n_consec_infeas;
    double BASE_SCALE_FACTOR = 1.0;
    bool SOLVE_EACH_REAL = false;
	bool reset_corr = false;

	//trust region parameters
	double trust_radius = 5.0;
	double trust_radius_max = 5.0;
	double trust_radius_min = 1e-4;
	
	vector<double> previous_obj_values;
	const int memory_length = 5;  // Number of previous objectives to remember
	const double c1 = 0.0001;  // Armijo condition parameter
	const double c2 = 0.9;     // curvature condition parameter

	Eigen::VectorXd diagonal_scaling;
	double adaptation_rate = 0.3;  // How quickly to adapt scaling (0-1)
	int scaling_update_frequency = 1;  // Update scaling every N iterations

	set<string> pp_args;

	int iter;

	double last_best;
	double last_viol;
	vector<double> best_phis, best_feas_phis;
	vector<double> best_violations;
	double best_phi_yet;
	double best_violation_yet;
	// adaptive state, not an option cache: seeded from sqp_working_set_tol at initialize()
	// and then tightened on each accepted fd iteration, so the option is an initial value
	// only (registered init-only, so a later change is surfaced rather than ignored)
	double working_set_tol;

	map<string, double> obj_map, total_viol_map;
	map<string, Eigen::VectorXd>  step_length_map;
	map<string, string> ls_parent_map;
	bool is_good_search = false;

	int warn_min_reals, error_min_reals;

	vector<string> oe_org_real_names, pe_org_real_names;
	vector<string> act_obs_names, act_par_names;
	vector<string> dv_names;
	vector<string> adj_par_names;
	// resolved init state: latched off during initialize() when cma was requested but could
	// not be set up (sqp_num_reals is init-only, so this cannot go stale mid-run)
	bool use_cma = true;

	Parameters current_ctl_dv_values, prev_ctl_dv_values, trial_ctl_dv_values, infeas_cand_dv_values;
	Observations current_obs, trial_obs, infeas_cand_obs;
	
	Parameters current_grad_vector;
	map<int, Parameters> grad_vector_map;

	Mat current_constraint_mat, prev_constraint_mat;
	Eigen::MatrixXd constraint_jco, base_constraint_jco;
	vector<string> cnames;

	ParameterEnsemble dv, dv_base, dv_to_save;
	ObservationEnsemble oe, oe_base, oe_to_save;
	map<string, string> constraint_sense;
	Eigen::VectorXd lambda;

	map <string, pair<Mat, bool>> constraint_mat_en;
	map<string, vector<string>> cnames_en;
	map<string, Eigen::VectorXd> search_d_en, lm_en;
	map<string, double> current_obj_en;
	map<string, Eigen::MatrixXd> constraint_jco_en;
	map<string, Covariance> hessian_en;
	vector<string> cnames_base;
	Eigen::VectorXd lm_base;

	void save_current_dv_obs();

	Constraints constraints;

	bool oe_drawn, dv_drawn;
	set<int> selected_dv_indices; 
	set<int> unselected_dv_indices;  
	bool sampling_tracking_initialized, cma_reset_archive = true;

	bool is_blocking_constraint = false;
	bool is_base_infeas = false;
	bool seek_ies = false;

	Jacobian_1to1 jco;
	Covariance hessian, used_hessian;
	// per-iteration cache of the regularized Hessian.  Within solve_new_ensemble's
	// per-realization loop every realization regularizes the identical shared Hessian,
	// so it is computed once per iteration and reused (avoids an O(n_dv^3) eigen-
	// decomposition per realization -- the dominant cost at high n_dv).
	Eigen::MatrixXd cached_reg_hessian;
	bool cached_reg_hessian_valid = false;
	// per-iteration cache of the UNCONSTRAINED Newton step (-H^-1 grad).  With an empty
	// working set and a shared objective gradient this is identical for every realization,
	// so the O(n_dv^3) LDLT factorization+solve is done once per iteration and reused.  The
	// cached gradient is compared so the cache is skipped if a realization's gradient differs.
	Eigen::VectorXd cached_unconstrained_search_d;
	Eigen::VectorXd cached_unconstrained_grad;
	Covariance cached_unconstrained_used_hessian;  // preserve regularize_hessian's used_hessian side effect on cache hit
	bool cached_unconstrained_valid = false;
	// cache of the dv-covariance pseudoinverse SVD (s,U,V).  It depends only on the dv
	// ensemble anomalies (which change only at make_gradient_runs), yet is needed identically
	// by calc_gradient_vector and calc_objective_hessian, so it is computed once per dv
	// ensemble and reused.  Keyed on the anomalies themselves, so it auto-invalidates when
	// the ensemble changes.  This SVD is ~350s per call at high n_dv.
	Eigen::MatrixXd cached_dvcov_s, cached_dvcov_U, cached_dvcov_V, cached_dv_anoms;
	bool cached_dvcov_valid = false;
	Eigen::VectorXd step_k;
	string selected_ls_parent, selected_ls_child;
	map<double, map<string, string>> sv_lineage_map;
	map<string, double> cname_sf_map;

	SqpFilter filter;

	void prep_4_ensemble_grad();
	void prep_4_fd_grad();

	bool update_hessian(string how);
	void update_scaling(const Eigen::VectorXd& step, const Eigen::VectorXd& grad);
	Eigen::MatrixXd regularize_hessian(const Eigen::MatrixXd& H, const string& context);
	bool try_modify_hessian();
	bool hessian_update_bfgs(Eigen::VectorXd s_k, Eigen::VectorXd y_k, Covariance old_hessian);
	bool hessian_update_sr1(Eigen::VectorXd s_k, Eigen::VectorXd y_k, Covariance old_hessian);
	Covariance calc_objective_hessian();
	// Cached SVD of the dv-covariance (V*s^-1*U^T is the StoSAG pseudoinverse).  Forms
	// dv_cov = dv_anoms^T*dv_anoms/(m-1) internally and caches the decomposition keyed on
	// the anomalies so repeated calls on the same dv ensemble reuse the O(n_dv^3) SVD.
	void dvcov_svd(const Eigen::MatrixXd& dv_anoms, Eigen::MatrixXd& s, Eigen::MatrixXd& U, Eigen::MatrixXd& V);

	bool seek_feasible();
	void generate_intermediate_candidates(const string& parent_name, double start_scale, double end_scale,	int num_points,	ParameterEnsemble* dvs_subset,	const map<string, Eigen::VectorXd>& search_d_map, ParameterEnsemble& dv_intermediate, vector<string>& intermediate_cand_names);
	FilterRec pick_upgrade_and_update_current(ParameterEnsemble& dv_candidates, ObservationEnsemble& _oe, bool cma_reset_arc = true, bool report = false, ParameterEnsemble* dvs_subset = nullptr, bool recalc = false);
	tuple<FilterRec, SqpFilter> pick_from_filter(ParameterEnsemble& dv_candidates, ObservationEnsemble& _oe, bool recalc = true);
	FilterRec pick_from_filter_by_merit(SqpFilter _filtered);

	double compute_actual_reduction(Parameters& trial_dv_values, Observations& trial_obs);
	double compute_predicted_reduction(const Eigen::VectorXd& step, const Eigen::VectorXd& grad);

	Eigen::VectorXd solve_trust_region_subproblem_dogleg(const Eigen::MatrixXd& B, const Eigen::VectorXd& g, double radius);
	Eigen::VectorXd solve_constrained_trust_region_step(const Eigen::MatrixXd& B, const Eigen::VectorXd& g, const Eigen::MatrixXd& A, double radius);


	Parameters calc_gradient_vector(const Parameters& _current_dv_values, string _center_on=string());
	Eigen::VectorXd calc_gradient_vector_from_coeffs(const Parameters & _current_dv_values);

	Eigen::VectorXd get_obj_vector(ParameterEnsemble& _dv, ObservationEnsemble& _oe);
	
	double get_obj_value(Parameters& _current_ctl_dv_vals, Observations& _current_obs);
	map<string, double> get_obj_map(ParameterEnsemble& _dv, ObservationEnsemble& _oe);
	pair<Mat, bool> get_constraint_mat(Parameters& _dv_vals, Observations&_obs_vals, double working_set_tol = 0.005, const Eigen::VectorXd* lagrange_mults = nullptr, vector<string> curr_ws = vector<string>());


	pair<Eigen::VectorXd, Eigen::VectorXd> _kkt_direct(const Eigen::MatrixXd& inv_hessian, Eigen::MatrixXd& _constraint_jco, Eigen::VectorXd& constraint_diff, Eigen::VectorXd& curved_grad, vector<string>* _cnames = nullptr);
	pair<Eigen::VectorXd, Eigen::VectorXd> _kkt_null_space(const Eigen::MatrixXd& inv_hessian, Eigen::MatrixXd& _constraint_jco, Eigen::VectorXd& constraint_diff, Eigen::VectorXd& curved_grad, vector<string>* _cnames = nullptr);

	FilterRec run_search_routine(Eigen::VectorXd& grad, ParameterEnsemble* dvs_subset = nullptr, bool recalc = false);

	void run_jacobian(Parameters& _current_dv_vals,Observations& _current_obs, bool init_obs);

	void make_gradient_runs(Parameters& _current_dv_vals, Observations& _current_obs);

	ObservationEnsemble combine_obs_and_pi(ObservationEnsemble& _oe, ParameterEnsemble& _pe);
	Ensemble get_pi_ensemble(ParameterEnsemble& _dv, vector<string>& pinames);

	void report_and_save_ensemble(ParameterEnsemble& _dv, ObservationEnsemble& _oe);
	void save(ParameterEnsemble& _dv, ObservationEnsemble& _oe, bool save_base=true);
	void save_mat(string prefix, Eigen::MatrixXd &mat);
	bool initialize_dv(Covariance &cov);
	bool initialize_restart();
	void initialize_parcov();
	void initialize_objfunc();
	void queue_chance_runs();

	template<typename T, typename A>
	void message(int level, const string &_message, vector<T, A> _extras, bool echo=true);
	void message(int level, const string &_message);

	template<typename T>
	void message(int level, const string &_message, T extra);

	void sanity_checks();
	bool isfullrank(const Eigen::MatrixXd& mat);

	void add_current_as_bases(ParameterEnsemble& _dv, ObservationEnsemble& _oe);

	vector<int> get_subset_idxs(int size, int nreal_subset);
};




#endif

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




struct FilterRec
{
	double obj_val;
	double viol_val;
	int iter;
	double alpha;
    friend bool operator<(const FilterRec &k1, const FilterRec &k2) {
        if ((k1.obj_val < k2.obj_val) && (k1.viol_val < k2.viol_val))
            return true;
        return false;
    }
};

//template<>
//struct std::less<FilterRec> {
//    bool operator()(const FilterRec &k1, const FilterRec &k2) const {
//        if ((k1.obj_val < k2.obj_val) && (k1.viol_val < k2.viol_val))
//            return true;
//        return false;
//    }
//};

class SqpFilter
{
public:
	SqpFilter(bool _minimize=true,double _obj_tol = 0.001, double _viol_tol = 0.001) {
		minimize = _minimize; obj_tol = _obj_tol; viol_tol = _viol_tol;
	}
	bool accept(double obj_val, double violation_val,int iter=0,double alpha=-1.0, bool keep=false);
	bool update(double obj_val, double violation_val, int iter=0,double alpha=-1.0);
    void report(ofstream& frec,int iter);
    double get_viol_tol() {return viol_tol;}
	void set_tol(double tol) { 
		obj_tol = tol; 
		viol_tol = tol;}

private:
	bool minimize;
	double obj_tol;
	double viol_tol;

	multiset<FilterRec> obj_viol_pairs;

	bool first_partially_dominates_second(const FilterRec& first, const FilterRec& second);
    bool first_strictly_dominates_second(const FilterRec& first, const FilterRec& second);

};

class CovMatAdapES 
{

public:
	CovMatAdapES(Pest* pest_ptr, std::mt19937* rand_gen) : pest_scenario_ptr(pest_ptr), rand_gen_ptr(rand_gen) {}
	CovMatAdapES() {}

	void initialize(int n_params, int _num_reals);
	void update(ParameterEnsemble curr_pe, ObservationEnsemble curr_oe, Parameters prev_m, Parameters curr_m);
	void set_covariance(const Eigen::MatrixXd& _C) { C = _C; }
	Covariance get_covariance_matrix(const vector<string>& par_names);

	void set_mean(const Eigen::VectorXd& _m) { m = _m; }
	Eigen::VectorXd get_mean() const { return m; }

	double get_sigma() const { return sigma; }

	ParameterEnsemble generate_population(Parameters& _curr_m, vector<string> _parname, ParameterEnsemble _dv);

private:
	// CMA-ES parameters
	int lambda;                    // Population size
	int mu;                       // Number of parents
	double sigma = 0.1;                 // Step size
	Eigen::VectorXd m;            // Mean vector
	Eigen::MatrixXd C;            // Covariance matrix
	Eigen::MatrixXd B;            // Eigenvectors of C
	Eigen::VectorXd D;            // Eigenvalues of C
	Eigen::MatrixXd pc;           // Evolution path
	Eigen::MatrixXd ps;           // Evolution path for sigma
	double c_sigma, c_c, c_1, c_mu, d_sigma, chi_n, c_m;
	vector<double> weights;

protected:
	std::mt19937* rand_gen_ptr;
	Pest* pest_scenario_ptr;
	Eigen::MatrixXd reals;

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
	//L2PhiHandler ph;
	ParChangeSummarizer pcs;
	Covariance parcov, obscov;
	CovMatAdapES cmaes;
	double reg_factor;
	chancePoints chancepoints;
	string obj_func_str;
	string obj_obs;
	string obj_sense;
	bool use_obj_obs;
	bool use_obj_pi;
	bool converged = false;
	map<string, double> obj_func_coef_map;

	int num_threads;
	int n_consec_infeas;
	//todo: make these ++ args
	int MAX_CONSEC_INFEAS = 100000;
    int MAX_CONSEC_INFEAS_IES = 3;
    int MAX_CONSEC_PHIINC = 100000;
    double SF_DEC_FAC = 0.5;
    double SF_INC_FAC = 1.1;
    double BASE_SCALE_FACTOR = 1.0;
    double PAR_SIGMA_DEC_FAC = 0.9;
    double PAR_SIGMA_INC_FAC = 2.0;
    bool SOLVE_EACH_REAL = false;
    double PHI_ACCEPT_FAC = 0.05;
    double par_sigma_max = 100;
	double penalty_param = 10.0;  // Base penalty parameter
    //todo add warning for par_sigma_range too low
    double par_sigma_min = 10;
	double eigthresh;
	bool reset_hessian = false;

	int n_consec_failures = 0;
	int max_consec_failures = 2; //put this somewhere else later
	int max_line_search_attempts = 3;

	//trust region parameters
	//TODO: Maybe put these as ++args later
	double trust_radius = 2.0;
	double trust_radius_max = 2.0;
	double trust_radius_min = 1e-4;
	double eta1 = 0.25;  // ratio threshold for radius reduction
	double eta2 = 0.75;  // ratio threshold for radius increase
	double gamma1 = 0.5; // radius reduction factor
	double gamma2 = 2.0; // radius increase factor
	const int batch_size = 10;
	
	vector<double> previous_obj_values;
	const int memory_length = 5;  // Number of previous objectives to remember
	double prev_successful_scale = 1.0;
	const double c1 = 0.0001;  // Armijo condition parameter
	const double c2 = 0.9;     // curvature condition parameter
	const double min_scale = 1e-8;

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
	double working_set_tol;

	int warn_min_reals, error_min_reals;

	vector<string> oe_org_real_names, pe_org_real_names;
	vector<string> act_obs_names, act_par_names;
	vector<string> dv_names;
	string best_name;
	bool use_subset, use_localization, use_cmaes = true;

	Parameters current_ctl_dv_values, prev_ctl_dv_values, trial_ctl_dv_values, infeas_cand_dv_values;
	Observations current_obs, trial_obs, infeas_cand_obs;

	Parameters current_grad_vector, prev_grad_vector;
	map<int, Parameters> grad_vector_map;

	Mat current_constraint_mat, prev_constraint_mat;
	Eigen::MatrixXd constraint_jco, base_constraint_jco;
	vector<string> cnames;

	ParameterEnsemble dv, dv_base;
	ObservationEnsemble oe, oe_base;
	map<string, string> constraint_sense;
	Eigen::VectorXd lambda, base_lambda;

	//these are used so that we can update the constraints based on the current best values
	//Parameters best_mean_dv_values;
	//Observations best_mean_obs_values;

	void save_current_dv_obs();

	Constraints constraints;

	bool oe_drawn, dv_drawn;
	set<int> selected_dv_indices; 
	set<int> unselected_dv_indices;  
	bool sampling_tracking_initialized;

	bool use_ensemble_grad;
	bool is_blocking_constraint = false;

	Jacobian_1to1 jco;

	//store the hessian as a cov since it is symmetric...
	Covariance hessian;

	SqpFilter filter;

	void prep_4_ensemble_grad();
	void prep_4_fd_grad();

	bool update_hessian();
	void update_scaling(const Eigen::VectorXd& step, const Eigen::VectorXd& grad);
	bool try_modify_hessian();
	bool hessian_update_bfgs(Eigen::VectorXd s_k, Eigen::VectorXd y_k, Covariance old_hessian);
	bool hessian_update_sr1(Eigen::VectorXd s_k, Eigen::VectorXd y_k, Covariance old_hessian);
	bool solve_new();
	bool solve_new_ensemble();

	bool seek_feasible();
	bool line_search(Eigen::VectorXd& search_d, const Parameters& _current_dv_values, Eigen::VectorXd& grad);

	bool line_search(map<string, Eigen::VectorXd>& search_d, Eigen::VectorXd& grad, ParameterEnsemble* dvs_subset = nullptr, bool first_pick = false);
	bool iterative_partial_step(const string& _blocking_constraint);
	bool pick_candidate_and_update_current(ParameterEnsemble& dv_candidates, ObservationEnsemble& _oe, map<string,double>& sf_map, bool first_pick = false);
	pair<bool, pair<Parameters, Observations>> pick_upgrade_and_update_current(ParameterEnsemble& dv_candidates, ObservationEnsemble& _oe);
	bool pick_partial_step(ParameterEnsemble& dv_candidates, ObservationEnsemble& _oe, map<string, double>& sf_map);
	bool check_wolfe_conditions(Parameters& trial_dv_values, Observations& trial_obs, const Eigen::VectorXd& search_d, 
		const Eigen::VectorXd& grad, double scale, double initial_obj, double initial_slope);
	double get_reference_obj();

	double compute_actual_reduction(Parameters& trial_dv_values, Observations& trial_obs);
	double compute_predicted_reduction(const Eigen::VectorXd& step, const Eigen::VectorXd& grad);

	Eigen::VectorXd compute_adaptive_localization_weights(const vector<string>& dv_names, const Eigen::MatrixXd& dv_anoms);
	
	bool trust_region_step(Parameters& current_dv_values, Eigen::VectorXd grad);
	Eigen::VectorXd solve_trust_region_subproblem_dogleg(const Eigen::MatrixXd& B, const Eigen::VectorXd& g, double radius);
	Eigen::VectorXd solve_trust_region_subproblem(const Eigen::MatrixXd& B, const Eigen::VectorXd& g, double radius);
	Eigen::VectorXd compute_boundary_solution(const Eigen::VectorXd& p,	const Eigen::VectorXd& d, double radius);
	Parameters calc_gradient_vector(const Parameters& _current_dv_values, string _center_on=string());
	
	Eigen::VectorXd calc_gradient_vector_from_coeffs(const Parameters & _current_dv_values);

	Eigen::VectorXd get_obj_vector(ParameterEnsemble& _dv, ObservationEnsemble& _oe);

	

	double get_obj_value(Parameters& _current_ctl_dv_vals, Observations& _current_obs);
	map<string, double> get_obj_map(ParameterEnsemble& _dv, ObservationEnsemble& _oe);
	pair<Mat, bool> get_constraint_mat(Parameters& _dv_vals, Observations&_obs_vals, double working_set_tol = 0.005, const Eigen::VectorXd* lagrange_mults = nullptr);

	//Eigen::VectorXd calc_search_direction_vector(const Parameters& _current_dv_, Eigen::VectorXd &
	// );
	pair<Eigen::VectorXd, Eigen::VectorXd> calc_search_direction_vector(Parameters& _current_dv_, Observations& _current_obs_values, Eigen::VectorXd& grad_vector, Eigen::MatrixXd* _constraint_jco ,vector<string>* _cnames = nullptr);

	pair<Eigen::VectorXd, Eigen::VectorXd> _kkt_direct(Eigen::MatrixXd& inv_hessian, Eigen::MatrixXd& _constraint_jco, Eigen::VectorXd& constraint_diff, Eigen::VectorXd& curved_grad, vector<string>& cnames);
	pair<Eigen::VectorXd, Eigen::VectorXd> _kkt_null_space(Eigen::MatrixXd& inv_hessian, Eigen::MatrixXd& _constraint_jco, Eigen::VectorXd& constraint_diff, Eigen::VectorXd& curved_grad);

	//Parameters fancy_solve_routine(double scale_val, const Parameters& _current_dv_);
	Eigen::VectorXd fancy_solve_routine(const Parameters& _current_dv_, const Parameters& _grad_vector);

	vector<int> run_ensemble(ParameterEnsemble &_pe, ObservationEnsemble &_oe, const vector<int> &real_idxs=vector<int>());
	ObservationEnsemble run_candidate_ensemble(ParameterEnsemble&dv_candidates);

	void run_jacobian(Parameters& _current_dv_vals,Observations& _current_obs, bool init_obs);

	void make_gradient_runs(Parameters& _current_dv_vals, Observations& _current_obs);

	void report_and_save_ensemble();
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

#ifndef DATAASSIMILATOR_H_
#define DATAASSIMILATOR_H_

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




class LocalUpgradeThread_da
{
public:

	LocalUpgradeThread_da(PerformanceLog* _performance_log, unordered_map<string, Eigen::VectorXd>& _par_resid_map, unordered_map<string, Eigen::VectorXd>& _par_diff_map,
		unordered_map<string, Eigen::VectorXd>& _obs_resid_map, unordered_map<string, Eigen::VectorXd>& _obs_diff_map,
		Localizer& _localizer, unordered_map<string, double>& _parcov_inv_map,
		unordered_map<string, double>& _weight_map, ParameterEnsemble& _pe_upgrade,
		unordered_map<string, pair<vector<string>, vector<string>>>& _cases,
		unordered_map<string, Eigen::VectorXd>& _Am_map, Localizer::How& _how,
		unordered_map<string, Eigen::VectorXd>& _obs_err_map);

	//Eigen::DiagonalMatrix<double, Eigen::Dynamic> get_matrix_from_map(vector<string> &names, map<string, double> &dmap);	
	//Eigen::MatrixXd get_matrix_from_map(int num_reals, vector<string> &names, map<string, Eigen::VectorXd> &emap);


	void work(int thread_id, int iter, double cur_lam);

	void work2da(int thread_id, int iter, double cur_lam);


private:
	PerformanceLog* performance_log;
	Localizer::How how;
	vector<string> keys;
	int count, total;
	//double eigthresh, cur_lam;
	//int maxsing, num_reals,iter, thread_id;
	//bool use_approx, use_prior_scaling;

	unordered_map<string, pair<vector<string>, vector<string>>>& cases;

	ParameterEnsemble& pe_upgrade;
	//PhiHandler_da &ph;
	Localizer& localizer;
	unordered_map<string, double>& parcov_inv_map;
	unordered_map<string, double>& weight_map;

	unordered_map<string, Eigen::VectorXd>& par_resid_map, & par_diff_map, & Am_map;
	unordered_map<string, Eigen::VectorXd>& obs_resid_map, & obs_diff_map, & obs_err_map;

	mutex ctrl_lock, weight_lock, loc_lock, parcov_lock;
	mutex obs_resid_lock, obs_diff_lock, par_resid_lock, obs_err_lock;
	mutex par_diff_lock, am_lock, put_lock;
	mutex next_lock;

};


class DataAssimilator: public EnsembleMethod
{
public:
	/*DataAssimilator(Pest& _pest_scenario, FileManager& _file_manager,
		OutputFileWriter& _output_file_writer, PerformanceLog* _performance_log,
		RunManagerAbstract* _run_mgr_ptr);*/

	using EnsembleMethod::EnsembleMethod;

	//void forward_run_noptmax_0(int icycle);
	//void initialize(int _icycle);
	//void da_save_ensemble_pe(string fprefix, string dtyp);
	//void da_save_ensemble_oe(string fprefix, string dtyp);
	//void add_dynamic_state_to_pe();
	//void add_dynamic_state_to_pe();
	
	
	void da_update();

	//void kf_upate();
	void finalize();
	//void throw_da_error(string message);
	//bool should_terminate();
	ParameterEnsemble get_pe() { return pe;}
	void set_pe(ParameterEnsemble new_pe) { pe = new_pe;}
	////bool use_ies; 
	//string da_type;

	//bool initialize_pe(Covariance& cov);
	//void initialize_parcov();
	Covariance* get_parcov_ptr() { return &parcov; }
	std::mt19937& get_rand_gen() { return rand_gen; }
	vector<string> get_act_par_names() { return act_par_names; }
	ObservationEnsemble& get_oe() { return oe; }
	L2PhiHandler& get_phi_handler() { return ph; }
	int get_iter() { return iter; }
	FileManager& get_file_manager() { return file_manager; }
	Pest& get_pest_scenario() { return pest_scenario; }
	string da_type;
private:
	int icycle;
	
	CtlPar_container da_ctl_params;

	//vector<string> obs_dyn_state_names;
	//vector<string> par_dyn_state_names;

	vector<double> infl_facs;
	int solution_iterations;
	ParameterEnsemble pe_post;
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> obscov_inv_sqrt, parcov_inv_sqrt;

	bool oe_drawn, pe_drawn;
	//void initialize_dynamic_states();
	bool solve_new_da();
	void update_starting_state();
	void return_post_dyn_state(vector<ParameterEnsemble>& pe_lams, vector<ParameterEnsemble> poterior_dyn_states);
	vector<ParameterEnsemble> temp_remove_dyn_state(vector<ParameterEnsemble>& pe_lams);
	
	//ParameterEnsemble calc_upgrade(vector<string> &obs_names, vector<string> &par_names,double lamb, int num_reals);

	//ParameterEnsemble calc_localized_upgrade(double cur_lam);
	//ParameterEnsemble calc_localized_upgrade_threaded(double cur_lam, unordered_map<string, pair<vector<string>, vector<string>>>& loc_map);

	//ParameterEnsemble calc_kf_upgrade(double cur_lam, unordered_map<string, pair<vector<string>, vector<string>>>& loc_map);

	void eig2csv(string name, Eigen::MatrixXd matrix);

	/*ParameterEnsemble kf_work(PerformanceLog* performance_log, unordered_map<string, 
		Eigen::VectorXd>& par_resid_map, unordered_map<string, Eigen::VectorXd>& par_diff_map,
		unordered_map<string, Eigen::VectorXd>& obs_resid_map, unordered_map<string, Eigen::VectorXd>& obs_diff_map,
		unordered_map<string, Eigen::VectorXd>& obs_err_map, Localizer& localizer,
		unordered_map<string, double>& parcov_inv_map, unordered_map<string, double>& weight_map,
		ParameterEnsemble& pe_upgrade, double cur_lam, unordered_map<string, pair<vector<string>, vector<string>>>& loc_map,
		unordered_map<string, Eigen::VectorXd>& Am_map, Localizer::How& how);*/
	/*
	ParameterEnsemble kf_work(PerformanceLog* _performance_log, unordered_map<string,
		Eigen::VectorXd>& _par_resid_map, unordered_map<string, Eigen::VectorXd>& _par_diff_map, 
		unordered_map<string, Eigen::VectorXd>& _obs_resid_map, unordered_map<string, 
		Eigen::VectorXd>& _obs_diff_map, unordered_map<string, Eigen::VectorXd>& obs_err_map,  Localizer& _localizer, unordered_map<string,
		double>& _parcov_inv_map, unordered_map<string,double>& _weight_map, ParameterEnsemble& _pe_upgrade, unordered_map<string, 
		pair<vector<string>, vector<string>>>& _cases, unordered_map<string, 
		Eigen::VectorXd>& _Am_map, Localizer::How& _how);*/

	
	//EnsemblePair run_ensemble(ParameterEnsemble &_pe, ObservationEnsemble &_oe);
	//vector<int> run_ensemble(ParameterEnsemble& _pe, ObservationEnsemble& _oe, const vector<int>& real_idxs = vector<int>());
	//vector<ObservationEnsemble> run_lambda_ensembles(vector<ParameterEnsemble>& pe_lams, vector<double>& lam_vals, vector<double>& scale_vals);
	//map<string, double> get_phi_vec_stats(map<string,PhiComponets> &phi_info);
	//map<string,PhiComponets> get_phi_info(ObservationEnsemble &_oe);
	//void report_and_save();
	//void save_mat(string prefix, Eigen::MatrixXd& mat);
	
	//bool initialize_oe(Covariance& cov);
	//void initialize_restart();
	
	//void initialize_obscov();
	//void drop_bad_phi(ParameterEnsemble& _pe, ObservationEnsemble& _oe, bool is_subset = false);
	//void check_ensembles(ObservationEnsemble &oe, ParameterEnsemble &pe);
	/*template<typename T, typename A>
	void message(int level, const string& _message, vector<T, A> _extras, bool echo = true);
	void message(int level, const string& _message);*/

	//template<typename T, typename A>
	//void message(int level, char* _message, vector<T, A> _extras);// { message(level, string(_message), _extras); }
	//void message(int level, char* _message);// { message(level, string(_message)); }

	//template<typename T>
	//void message(int level, const string& _message, T extra);

	//template<typename T>
	//void message(int level, char* _message, T extra);

	void sanity_checks();

	//void add_bases();

	//void update_reals_by_phi(ParameterEnsemble& _pe, ObservationEnsemble& _oe);

	//void initialize();

	//vector<string> detect_prior_data_conflict();

	//map<int,int> get_subset_idx_map();
	//void set_subset_idx(int size);
	//Eigen::MatrixXd get_Am(const vector<string>& real_names, const vector<string>& par_names);

};


map<int, map<string, double>> process_da_par_cycle_table(Pest& pest_scenario, vector <int>& ncycles_in_tables, ofstream& fout_rec);
map<int, map<string, double>> process_da_obs_cycle_table(Pest& pest_scenario, vector <int>& ncycles_in_tables, ofstream& fout_rec, set<string>& obs_in_tbl);
map<int, map<string, double>> process_da_weight_cycle_table(Pest& pest_scenario, vector <int>& ncycles_in_tables, ofstream& fout_rec, set<string>& obs_in_tbl);

void write_global_phi_info(int cycle, ofstream& f_phi, DataAssimilator& da, vector<string>& init_real_names);

void generate_global_ensembles(DataAssimilator& da, ofstream& fout_rec, ParameterEnsemble& curr_pe, ObservationEnsemble& curr_oe);



#endif

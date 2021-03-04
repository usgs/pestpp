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


enum chancePoints { ALL, SINGLE };




class L2PhiHandler
{
public:

	enum phiType { MEAS, COMPOSITE, REGUL, ACTUAL };
	L2PhiHandler() { ; }
	L2PhiHandler(Pest *_pest_scenario, FileManager *_file_manager,
		       ObservationEnsemble *_oe_base, ParameterEnsemble *_pe_base,
		       Covariance *_parcov, bool should_prep_csv = true, string _tag=string());
	void update(ObservationEnsemble &oe, ParameterEnsemble &pe);
	double get_mean(phiType pt);
	double get_std(phiType pt);
	double get_max(phiType pt);
	double get_min(phiType pt);

	double calc_mean(map<string, double> *phi_map);
	double calc_std(map<string, double> *phi_map);

	map<string, double>* get_phi_map_ptr(L2PhiHandler::phiType pt);
	map<string, double> get_phi_map(L2PhiHandler::phiType pt);
	void report(bool echo=true);
	void write(int iter_num, int total_runs, bool write_group = true);
	void write_group(int iter_num, int total_runs, vector<double> extra);
	vector<int> get_idxs_greater_than(double bad_phi, double bad_phi_sigma, ObservationEnsemble &oe);

	Eigen::MatrixXd get_obs_resid(ObservationEnsemble &oe, bool apply_ineq=true);
	Eigen::MatrixXd get_obs_resid_subset(ObservationEnsemble &oe, bool apply_ineq=true);

	Eigen::MatrixXd get_par_resid(ParameterEnsemble &pe);
	Eigen::MatrixXd get_par_resid_subset(ParameterEnsemble &pe);
	Eigen::MatrixXd get_actual_obs_resid(ObservationEnsemble &oe);
	Eigen::VectorXd get_q_vector();
	vector<string> get_lt_obs_names() { return lt_obs_names; }
	vector<string> get_gt_obs_names() { return gt_obs_names; }

	void apply_ineq_constraints(Eigen::MatrixXd &resid, vector<string> &names);

	void save_residual_cov(ObservationEnsemble& oe, int iter);


private:
	string tag;
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

	double org_reg_factor;
	Eigen::VectorXd org_q_vec;
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
	ParChangeSummarizer(ParameterEnsemble *_base_pe_ptr, FileManager *_file_manager_ptr, OutputFileWriter* _output_file_writer_ptr);
	void summarize(ParameterEnsemble &pe, int iiter, string filename = string());
	

private:
	ParameterEnsemble * base_pe_ptr;
	FileManager *file_manager_ptr;
	OutputFileWriter* output_file_writer_ptr;
	map<string, set<string>> pargp2par_map;
	pair<map<string,double>, map<string, double>> init_moments;
	map<string, double> mean_change;
	map<string, double> std_change;
	map<string, int> num_at_bounds;
	map<string, int> percent_at_bounds;

	void update(ParameterEnsemble& pe);
	void write_to_csv(string& filename);

};

pair<Parameters,Observations> save_real_par_rei(Pest& pest_scenario, ParameterEnsemble& pe, ObservationEnsemble& oe,
	OutputFileWriter& output_file_writer, FileManager& file_manager, int iter, string tag = BASE_REAL_NAME);

vector<int> run_ensemble_util(PerformanceLog* performance_log, ofstream& frec, ParameterEnsemble& _pe,
	ObservationEnsemble& _oe, RunManagerAbstract* run_mgr_ptr,
	bool check_pe_consistency = false, const vector<int>& real_idxs = vector<int>(),int da_cycle=NetPackage::NULL_DA_CYCLE);

class EnsembleSolver
{
public:
	EnsembleSolver(PerformanceLog* _performance_log, FileManager& _file_manager, Pest& _pest_scenario, ParameterEnsemble& _pe,
		ObservationEnsemble& _oe, ObservationEnsemble& _base_oe, Localizer& _localizer, Covariance& _parcov,Eigen::MatrixXd& _Am, L2PhiHandler& _ph, 
		bool _use_localizer, int _iter, vector<string>& _act_par_names, vector<string> &_act_obs_names);

	void solve(int num_threads, double cur_lam, bool use_glm_form, ParameterEnsemble& pe_upgrade, unordered_map<string, pair<vector<string>, vector<string>>>& loc_map);
	
private:
	PerformanceLog* performance_log;
	FileManager& file_manager;
	int iter, verbose_level;
	bool use_localizer;
	Pest& pest_scenario;
	ParameterEnsemble& pe;
	ObservationEnsemble& oe, base_oe;
	Localizer& localizer;
	Covariance& parcov;
	Eigen::MatrixXd& Am;
	L2PhiHandler& ph;
	unordered_map<string, Eigen::VectorXd> par_resid_map, obs_resid_map, Am_map;
	unordered_map<string, Eigen::VectorXd> par_diff_map, obs_diff_map, obs_err_map;
	unordered_map<string, double> weight_map;
	unordered_map<string, double> parcov_inv_map;
	//unordered_map<string, pair<vector<string>, vector<string>>> loc_map;

	template<typename T, typename A>
	void message(int level, const string& _message, vector<T, A> _extras, bool echo = true);
	void message(int level, const string& _message);

	template<typename T>
	void message(int level, const string& _message, T extra);

	
};


class LocalAnalysisUpgradeThread
{
public:

	LocalAnalysisUpgradeThread(PerformanceLog* _performance_log, unordered_map<string, Eigen::VectorXd>& _par_resid_map, unordered_map<string, Eigen::VectorXd>& _par_diff_map,
		unordered_map<string, Eigen::VectorXd>& _obs_resid_map, unordered_map<string, Eigen::VectorXd>& _obs_diff_map, unordered_map<string, Eigen::VectorXd>& _obs_err_map,
		Localizer& _localizer, unordered_map<string, double>& _parcov_inv_map,
		unordered_map<string, double>& _weight_map, ParameterEnsemble& _pe_upgrade,
		unordered_map<string, pair<vector<string>, vector<string>>>& _cases,
		unordered_map<string, Eigen::VectorXd>& _Am_map, Localizer::How& _how);

	void work(int thread_id, int iter, double cur_lam,bool use_glm_form);


private:
	PerformanceLog* performance_log;
	Localizer::How how;
	vector<string> keys;
	int count, total;

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

	void fill_maps();

};

#endif

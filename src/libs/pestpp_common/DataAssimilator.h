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
	
	using EnsembleMethod::EnsembleMethod;

	void da_update(int cycle);
	void finalize();
	ParameterEnsemble get_pe() { return pe;}
	void set_pe(ParameterEnsemble new_pe) { pe = new_pe;}
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
	CtlPar_container da_ctl_params;
	void eig2csv(string name, Eigen::MatrixXd matrix);	
	void sanity_checks();

};


map<int, map<string, double>> process_da_par_cycle_table(Pest& pest_scenario, vector <int>& ncycles_in_tables, ofstream& fout_rec);
map<int, map<string, double>> process_da_obs_cycle_table(Pest& pest_scenario, vector <int>& ncycles_in_tables, ofstream& fout_rec, set<string>& obs_in_tbl);
map<int, map<string, double>> process_da_weight_cycle_table(Pest& pest_scenario, vector <int>& ncycles_in_tables, ofstream& fout_rec, set<string>& obs_in_tbl);

void write_global_phi_info(int cycle, ofstream& f_phi, DataAssimilator& da, vector<string>& init_real_names);

void generate_global_ensembles(DataAssimilator& da, ofstream& fout_rec, ParameterEnsemble& curr_pe, ObservationEnsemble& curr_oe);



#endif

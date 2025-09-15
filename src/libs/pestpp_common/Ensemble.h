#ifndef ENSEMBLE_H_
#define ENSEMBLE_H_

#include <map>
#include <random>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <mutex>
#include "FileManager.h"
#include "ObjectiveFunc.h"
#include "OutputFileWriter.h"
#include "PerformanceLog.h"
#include "RunStorage.h"
#include "covariance.h"
#include "RunManagerAbstract.h"
#include "PerformanceLog.h"


const string BASE_REAL_NAME = "BASE";
const string MEDIAN_CENTER_ON_NAME = "_MEDIAN_";

class Ensemble
{
public:
	Ensemble(Pest* _pest_scenario, std::mt19937* _rand_gen_ptr);
	Ensemble(Pest* _pest_scenario);
	Ensemble(){ ; }

	void to_csv(string file_name);
	void to_csv_by_reals(ofstream& csv, bool write_header = true);
	void to_csv_by_vars(ofstream& csv, bool write_header = true);
	void to_binary(string file_name, bool transposed=false);
    void to_dense(string file_name);
	void from_eigen_mat(Eigen::MatrixXd mat, const vector<string> &_real_names, const vector<string> &_var_names);
	pair<int, int> shape() { return pair<int, int>(reals.rows(), reals.cols()); }
	void throw_ensemble_error(string message);
	void throw_ensemble_error(string message,vector<string> vec);
	const vector<string> get_var_names() const { return var_names; }
	const vector<string> get_real_names() const { return real_names; }

	const vector<string> get_real_names(vector<int> &indices);

	void extend_cols(Eigen::MatrixXd &_reals, const vector<string> &_var_names);
	void add_2_cols_ip(Ensemble &other);
	void add_2_cols_ip(const vector<string> &_var_names, const Eigen::MatrixXd &mat);
	void add_2_row_ip(const string& real_name,const Eigen::VectorXd& row_vec);

	Ensemble zeros_like(int nrow=-1);

	void broadcast_vec2mat(const vector<string>& other_var_names, const Eigen::MatrixXd& mat);

	void replace_col_vals(const vector<string>& other_var_names, const Eigen::MatrixXd& mat);

	vector<string> get_generic_real_names(int num_reals);

	void reserve(vector<string> _real_names, vector<string> _var_names);
	void replace_col(string var_name, Eigen::VectorXd & vec, bool update_map=true);
	Eigen::VectorXd get_real_vector(int ireal);
	Eigen::VectorXd get_real_vector(const string &real_name);
	Eigen::VectorXd get_var_vector(const string& var_name);
	void update_real_ip(const string &rname, Eigen::VectorXd &real);
	Eigen::MatrixXd get_eigen(vector<string> row_names, vector<string> col_names, bool update_vmap=true);
	const Eigen::MatrixXd get_eigen() const { return reals; }
	const Eigen::MatrixXd* get_eigen_ptr() const { return &reals; }
	Eigen::MatrixXd* get_eigen_ptr_4_mod() { return &reals; }
    map<string,double> get_real_map(string real_name,bool forgive=false);

	void set_eigen(Eigen::MatrixXd _reals);
	void reset_org_real_names() { org_real_names = real_names; }

	Eigen::MatrixXd get_eigen_anomalies(string on_real="");
	Eigen::MatrixXd get_eigen_anomalies(const vector<string> &_real_names, const vector<string> &_var_names, string on_real="");


	vector<double> get_mean_stl_var_vector();
	void fill_moment_maps(map<string, double>& mean_map, map<string, double>& std_map);
	void append_other_rows(Ensemble &other,bool reset_org_real_names=false);
	void append_other_rows(const vector<string>& _real_names, Eigen::MatrixXd& _reals);
	void append(string real_name, const Transformable &trans);
    void append(string real_name, const Eigen::VectorXd& vec);
	void replace(int idx, const Transformable &trans, string real_name="");

	Covariance get_diagonal_cov_matrix();
	pair<Covariance,Covariance> get_empirical_cov_matrices(FileManager* file_manager_ptr);
	void reorder(const vector<string> &_real_names, const vector<string> &_var_names, bool update_org_real_names=false);
	void drop_rows(const vector<int> &row_idxs, bool update_org_real_names=false);
	void drop_rows(const vector<string> &drop_names, bool update_org_real_names=false);
	void drop_cols(const vector<string>& drop_names);
	void keep_rows(const vector<int> &row_idxs);
	void keep_rows(const vector<string> &keep_names);
	

	Pest* get_pest_scenario_ptr() { return pest_scenario_ptr; }
	Pest get_pest_scenario() { return *pest_scenario_ptr; }
	void set_pest_scenario_ptr(Pest *_pest_scenario) { pest_scenario_ptr = _pest_scenario;}
	void set_real_names(vector<string> &_real_names, bool update_org_names=false);
	void check_for_dups();
	void check_for_normal(string context="");

	void draw(int num_reals, Covariance cov, Transformable &tran, const vector<string> &draw_names, const map<string,vector<string>> &grouper, PerformanceLog *plog, int level);
	void update_var_map();
	~Ensemble();
	void set_rand_gen(std::mt19937* _rand_gen_ptr) { rand_gen_ptr = _rand_gen_ptr; }
	std::mt19937* get_rand_gen_ptr() { return rand_gen_ptr; }
	map<string, int> get_var_map() { return var_map; }
	map<string, int> get_real_map();

	bool try_align_other_rows(PerformanceLog* performance_log, Ensemble& other);

protected:
	std::mt19937* rand_gen_ptr;
	Pest* pest_scenario_ptr;
	Eigen::MatrixXd reals;
	vector<string> var_names;
	vector<string> real_names;	
	vector<string> org_real_names;
	map<string, int> var_map;
	void read_csv_by_reals(int num_reals,ifstream &csv, map<string,int> &header_info, map<string,int> &index_info);
	void read_csv_by_vars(int num_reals, ifstream &csv, map<string, int> &header_info, map<string, int> &index_info);
	map<string,int> from_binary_old(string file_name, vector<string> &names,  bool transposed);
	map<string, int> from_binary(string file_name, vector<string> &names, bool transposed);
	pair<map<string, int>, map<string, int>> prepare_csv(const vector<string> &names, ifstream &csv, bool forgive);
	
};

class FixedParInfo
{
public:
	FixedParInfo(vector<string> _fixed_names);
	FixedParInfo() { initialized=false; }
	void set_fixed_names(vector<string>& _fixed_names) { fixed_names = _fixed_names; initialize(); }
	bool get_fixed_value(const string& pname, const string& rname, double& value);
	map<string, double> get_par_fixed_values(const string& pname);
	vector<double> get_real_fixed_values(const string& rname, vector<string>& pnames);
	map<string, double> get_real_fixed_values(const string& rname);
	void add_realization(string rname, Eigen::VectorXd& rvals, vector<string>& pnames);
    void add_realization(string rname, map<string, double>& rvals);
	void keep_realizations(const vector<string>& keep);
	void update_realizations(const vector<string>& other_var_names, const vector<string>& other_real_names, const Eigen::MatrixXd& other_mat);
	void add_realizations(map<string,map<string,double>>& other_fixed_info);
	void update_par_values(const map<string, double>& pval_map);
	void clear() { fixed_info.clear(); fixed_names.clear(); }
	void fill_fixed(map<string, double>& fixed_map, vector<string>& rnames);
	int get_map_size() {return fixed_info.size();}
	vector<string> get_real_names();
    vector<string> get_fixed_names() {return fixed_names;}
    map<string,map<string,double>> get_fixed_info_map() {return fixed_info;}


private:
	bool initialized;
	vector<string> fixed_names;
	map<string, map<string, double>> fixed_info;
	void initialize();
};

class ParameterEnsemble : public Ensemble
{

public:
	enum transStatus { CTL, NUM, MODEL };
	ParameterEnsemble(Pest *_pest_scenario_ptr, std::mt19937* rand_gen_ptr);
	ParameterEnsemble(Pest *_pest_scenario_ptr, std::mt19937* rand_gen_ptr, Eigen::MatrixXd _reals, vector<string> _real_names, vector<string> _var_names);
	ParameterEnsemble(Pest* _pest_scenario);
	ParameterEnsemble() { ; }
	ParameterEnsemble zeros_like(int nrows = -1);


	void set_zeros();
	void from_csv(string file_name, bool forgive = false);	
	void from_binary(string file_name, bool forgive = false);

	void from_eigen_mat(Eigen::MatrixXd mat, const vector<string> &_real_names, const vector<string> &_var_names,
		transStatus _tstat = transStatus::NUM);
	map<string,double>  enforce_bounds(PerformanceLog* plog, bool shrink);
	map<string,double> enforce_change_limits_and_bounds(PerformanceLog* plog, ParameterEnsemble& other);
	void to_csv(string file_name);
	void to_csv_by_reals(ofstream &csv, bool write_header=true);
	void to_csv_by_vars(ofstream &csv, bool write_header=true);
	//Pest* get_pest_scenario_ptr() { return &pest_scenario; }
	transStatus get_trans_status() const { return tstat; }
	void set_trans_status(transStatus _tstat) { tstat = _tstat; }
	void set_pest_scenario_ptr(Pest* _pest_scenario) { pest_scenario_ptr = _pest_scenario; par_transform = _pest_scenario->get_base_par_tran_seq(); }
	ParamTransformSeq get_par_transform() const { return par_transform; }
	void transform_ip(transStatus to_tstat);
	void set_pest_scenario(Pest *_pest_scenario);
	map<int,int> add_runs(RunManagerAbstract *run_mgr_ptr,const vector<int> &real_idxs=vector<int>(),
	        int da_cycle=NetPackage::NULL_DA_CYCLE, string additional_tag="");
	void set_fixed_names();
	void draw_uniform(int num_reals, vector<string> par_names, PerformanceLog* plog, int level, ofstream& frec);
	map<string,double> draw(int num_reals, Parameters par, Covariance &cov, PerformanceLog *plog, int level, ofstream& frec);
	Covariance get_diagonal_cov_matrix();
	void to_binary_ordered(string filename);
	void to_binary_unordered(string filename);
	void to_binary(string filename);
	void to_dense(string filename);
	void to_dense_unordered(string filename);
	void to_dense_ordered(string filename);
	void clear_fixed_map() { pfinfo.clear(); }
	void replace_col_vals_and_fixed(const vector<string>& other_var_names, const Eigen::MatrixXd& mat);
	FixedParInfo& get_fixed_info() { return pfinfo; }
	void set_fixed_info(FixedParInfo _pfinfo) { pfinfo = _pfinfo; }
	void keep_rows(const vector<int>& keep, bool update_fixed_map = false);
	void keep_rows(const vector<string>& keep, bool update_fixed_map = false);

private:
	ParamTransformSeq par_transform;
	transStatus tstat;
	void save_fixed(vector<string>& fixed_names);
	void fill_fixed(const map<string, int> &header_info, vector<string>& fixed_names);
	void replace_fixed(string real_name,Parameters &pars);
	void prep_par_ensemble_after_read(map<string,int>& header_info);
	FixedParInfo pfinfo;
};

class ObservationEnsemble : public Ensemble
{
public:
	ObservationEnsemble(Pest *_pest_scenario_ptr, std::mt19937* rand_gen_ptr);
	ObservationEnsemble(Pest *_pest_scenario_ptr, std::mt19937* rand_gen_ptr, Eigen::MatrixXd _reals, vector<string> _real_names, vector<string> _var_names);
	ObservationEnsemble(Pest* _pest_scenario);
	ObservationEnsemble() { ; }
	void to_binary(string filename) { Ensemble::to_binary(filename, true); }
    void to_dense(string filename) { Ensemble::to_dense(filename); }

	void update_from_obs(int row_idx, Observations &obs);
	void update_from_obs(string real_name, Observations &obs);
	void from_csv(string file_name);
	void from_eigen_mat(Eigen::MatrixXd mat, const vector<string> &_real_names, const vector<string> &_var_names);
	void from_binary(string file_name);
	vector<int> update_from_runs(map<int,int> &real_run_ids, RunManagerAbstract *run_mgr_ptr);
	vector<int> update_from_runs(map<int, int>& real_run_ids, RunManagerAbstract* run_mgr_ptr, ParameterEnsemble& run_mgr_pe);

	void draw(int num_reals, Covariance &cov, PerformanceLog *plog, int level, ofstream& frec);
	void initialize_without_noise(int num_reals);
};



class DrawThread
{
public:

	DrawThread(PerformanceLog *_performance_log, Covariance &_cov,Eigen::MatrixXd *_draws_ptr,
		vector<string> &_group_keys, const map<string,vector<string>> &_grouper);
	
	void work(int thread_id, int num_reals, int ies_verbose, map<string, int> idx_map, map<string, double> std_map);


private:
    //int num_real, ies_verbose;
	Eigen::MatrixXd *draws_ptr;
	PerformanceLog* performance_log;
	Covariance cov;
	//Eigen::MatrixXd draws;
	vector<string> group_keys;
	//map<string, int> idx_map;
	map<string, vector<string>> grouper;
	//map<string, double> std_map;
	mutex cov_lock, draw_lock, key_lock, pfm_lock;
};


#endif

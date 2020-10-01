#ifndef MOEA_H_
#define MOEA_H_

#include <unordered_map>
#include <random>
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <string>
#include "FileManager.h"
#include "ObjectiveFunc.h"
#include "OutputFileWriter.h"
#include "PerformanceLog.h"
#include "Ensemble.h"
#include "constraints.h"
#include "EnsembleMethodUtils.h"


class ParetoObjectives
{
public:
	ParetoObjectives(Pest& _pest_scenario, FileManager& _file_manager, 
		PerformanceLog* _performance_log,Constraints* _constraints_ptr = nullptr);

	pair<vector<string>, vector<string>> pareto_dominance_sort(ObservationEnsemble& op, ParameterEnsemble& dp, bool report=true);
	
	//this must be called at least once before the diversity metrixs can be called...
	void set_pointers(vector<string>& _obs_obj_names, vector<string>& _pi_obj_names, map<string, double>& _obj_dir_mult)
	{
		obs_obj_names_ptr = &_obs_obj_names; 
		pi_obj_names_ptr = &_pi_obj_names; 
		obj_dir_mult_ptr = &_obj_dir_mult;
	}
	void update_member_struct(ObservationEnsemble& oe, ParameterEnsemble& dp);
	
	map<string,double> get_crowding_distance(ObservationEnsemble& oe, ParameterEnsemble& dp);
	
	map<string, double> get_hypervolume(ObservationEnsemble& oe, ParameterEnsemble& dp);
	
	 
private:
	Pest& pest_scenario;
	FileManager& file_manager;
	PerformanceLog* performance_log;
	//vector<string> obj_names;
	Constraints* constraints_ptr;
	vector<string> sort_members_by_crowding_distance(vector<string>& members);
	bool first_dominates_second(map<string, double>& first, map<string, double>& second);

	void drop_duplicates(ObservationEnsemble& op, ParameterEnsemble& dp);
	bool first_equals_second(map<string, double>& first, map<string, double>& second);

	map<int, vector<string>> sort_members_by_dominance_into_fronts(map<string, map<string, double>>& member_struct);

	//sort specific members
	map<string, double> get_crowding_distance(vector<string>& members);
	//sort all members in member struct
	map<string, double> get_crowding_distance();

	//constraint dominance principal Fan (2017) Tian (2020)
	//void apply_constraints_cdp();
	//self-adaptive constraint penalty Fan (2017), Woldesenbet (2009)
	//void apply_constraints_sp();
	//in-feasibility driven evolution-ary algorithms Fan (2017)
	//void apply_constraints_idea();

	map<string, map<string, double>> member_struct;
	vector<string>* obs_obj_names_ptr;
	vector<string>* pi_obj_names_ptr;
	map<string, double>* obj_dir_mult_ptr;

	typedef std::function<bool(std::pair<std::string, double>, std::pair<std::string, double>)> Comparator;
	// Defining a lambda function to compare two pairs. It will compare two pairs using second field
	Comparator compFunctor = [](std::pair<std::string, double> elem1, std::pair<std::string, double> elem2)
	{
		return elem1.second < elem2.second;
	};

};


class MOEA
{
public:
	static mt19937_64 rand_engine;
	MOEA(Pest &_pest_scenario, FileManager &_file_manager, OutputFileWriter &_output_file_writer,
		PerformanceLog *_performance_log, RunManagerAbstract* _run_mgr_ptr);
	void initialize();
    void iterate_to_solution();
	void finalize();
	typedef pair<vector<string>, vector<string>> DomPair;
private:
	Pest& pest_scenario;
	set<string> pp_args;
	vector<string> act_obs_names, act_par_names;
	int iter, warn_min_members, error_min_members;
	int member_count;
	int archive_size;
	string population_dv_file, population_obs_restart_file;
	string dv_pop_file_tag = "dv_pop";
	string obs_pop_file_tag = "obs_pop";
	chancePoints chancepoints;
	FileManager &file_manager; 
	std::mt19937 rand_gen;
	vector<string> obs_obj_names, pi_obj_names;
	vector<string> dv_names;
	map<string, double> obj_dir_mult;

	map<string, map<string, double>> previous_obj_summary;


	//these two instances are passed as pointers to the constraints
	//Parameters effective_constraint_pars;
	//Observations effective_constraint_obs;

	ParetoObjectives objectives;
	Constraints constraints;
	const ParameterInfo *ctl_par_info_ptr;
	const ParameterGroupInfo *par_group_info_ptr;
	ParamTransformSeq par_transform;
	OutputFileWriter &output_file_writer;
	PerformanceLog *performance_log;
	RunManagerAbstract* run_mgr_ptr;
	const ObservationInfo *obs_info_ptr;


	ParameterEnsemble dp, dp_archive;
	ObservationEnsemble op, op_archive;

	void update_archive(ObservationEnsemble& _op, ParameterEnsemble& _dp);

	void throw_moea_error(const string& message);

	template<typename T, typename A>
	void message(int level, const string& _message, vector<T, A> _extras, bool echo = true);
	void message(int level, const string& _message);
	template<typename T>
	void message(int level, const string& _message, T extra);

	void sanity_checks();
	vector<int> run_population(ParameterEnsemble& _dp, ObservationEnsemble& _op);

	void queue_chance_runs(ParameterEnsemble& _dp);
	ObservationEnsemble get_chance_shifted_op(ObservationEnsemble& _op);

	bool initialize_dv_population();
	void initialize_obs_restart_population();

	ParameterEnsemble generate_population();

	ParameterEnsemble generate_diffevol_population(int num_members, ParameterEnsemble& _dp);
	ParameterEnsemble generate_nsga2_population(int num_members, ParameterEnsemble& _dp);

	string get_new_member_name(string tag = string());

	void save_populations(ParameterEnsemble& dp, ObservationEnsemble& op, string tag = string());
	void mutate(double probability, double eta_m, ParameterEnsemble& temp_dp);
	pair<Eigen::VectorXd, Eigen::VectorXd> crossover(double probability, double eta_m, int idx1, int idx2);

	pair<Parameters, Observations> get_optimal_solution(ParameterEnsemble& _dp, ObservationEnsemble& _oe, bool use_mean=false);

	map<string, map<string, double>> obj_func_report(ParameterEnsemble& _dp, ObservationEnsemble& _oe);
	map<string, map<string, double>> get_obj_func_summary_stats(ParameterEnsemble& _dp, ObservationEnsemble& _op);
	map<string, map<string, double>> obj_func_change_report(map<string, map<string, double>>& current_obj_summary);

	int get_max_len_obj_name();
};

#endif //MOEA_H_

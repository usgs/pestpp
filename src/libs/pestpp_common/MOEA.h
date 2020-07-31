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


class ParetoObjectives
{
public:
	ParetoObjectives(Pest& _pest_scenario, FileManager& _file_manager, PerformanceLog* _performance_log);
	vector<string> pareto_dominance_sort(const vector<string>& obj_names, ObservationEnsemble& op, ParameterEnsemble& dp);
private:
	Pest& pest_scenario;
	FileManager& file_manager;
	PerformanceLog* performance_log;
	vector<string> obj_names;

	vector<string> sort_members_by_crowding_distance(vector<string>& members, map<string, map<string, double>>& memeber_struct);
	bool first_dominates_second(map<string, double>& first, map<string, double>& second);
	map<int, vector<string>> sort_members_by_dominance_into_fronts(map<string, map<string, double>>& member_struct);

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

private:
	Pest& pest_scenario;
	set<string> pp_args;
	vector<string> act_obs_names, act_par_names;
	int iter, warn_min_members, error_min_members;
	int member_count;
	string population_dv_file, population_obs_restart_file;
	string dv_pop_file_tag = "dv_pop";
	string obs_pop_file_tag = "obs_pop";
	FileManager &file_manager; 
	std::mt19937 rand_gen;
	vector<string> obj_names;
	vector<string> dv_names;

	ParetoObjectives objectives;
	Constraints constraints;
	const ParameterInfo *ctl_par_info_ptr;
	const ParameterGroupInfo *par_group_info_ptr;
	ParamTransformSeq par_transform;
	OutputFileWriter &output_file_writer;
	PerformanceLog *performance_log;
	RunManagerAbstract* run_mgr_ptr;
	const ObservationInfo *obs_info_ptr;
	const PriorInformation *prior_info_ptr;

	ParameterEnsemble dp;
	ObservationEnsemble op;

	void throw_moea_error(const string& message);

	template<typename T, typename A>
	void message(int level, const string& _message, vector<T, A> _extras, bool echo = true);
	void message(int level, const string& _message);
	template<typename T>
	void message(int level, const string& _message, T extra);

	void sanity_checks();
	vector<int> run_population(ParameterEnsemble& _pe, ObservationEnsemble& _oe, const vector<int>& real_idxs = vector<int>());

	bool initialize_dv_population();
	void initialize_obs_restart_population();

	void generate_nsga2_member();
	void generate_de_member();
	void generate_spea2_member();

};

#endif //MOEA_H_

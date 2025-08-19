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
#include "RedSVD-h.h"


const string POP_SUM_TAG = "pareto.summary.csv";
const string ARC_SUM_TAG = "pareto.archive.summary.csv";
const string ARC_TRIM_SUM_TAG = "pareto.trimmed.archive.summary.csv";
const string RISK_NAME = "_RISK_";
const string DE_F_NAME = "_DE_F_";
const string CR_NAME = "_CR_";
const string MR_NAME = "_MR_";
const double CROWDING_EXTREME = 1.0e+30;
const double FLOAT_EPSILON = 1.0e-10;
enum MouGenType { DE, SBX, PM, PSO, SMP };
enum MouEnvType { NSGA, SPEA, NSGA_PPD }; //added NSGA_PPD for probabilistic Pareto dominance

enum MouMateType { RANDOM, TOURNAMENT };

class ParetoObjectives
{
public:
	ParetoObjectives(Pest& _pest_scenario, FileManager& _file_manager, 
		PerformanceLog* _performance_log);

	pair<vector<string>, vector<string>> get_nsga2_pareto_dominance(int generation, ObservationEnsemble& op, 
		ParameterEnsemble& dp, Constraints* constraints_ptr=nullptr, bool ppd=false, bool report=true, string sum_tag=string());
	
	map<string, map<string, double>> get_members(ObservationEnsemble& op, ParameterEnsemble& dp) { return get_member_struct(op, dp); };
	void set_ppd_beta() { ppd_beta = pest_scenario.get_pestpp_options().get_mou_ppd_beta(); }
	void set_prob_pareto(bool ppd) { prob_pareto = ppd; }
	void set_hypervolume_partitions(map<string, map<string, double>> _hv_parts);
	void get_ehvi(ObservationEnsemble& op, ParameterEnsemble& dp);

	//this must be called at least once before the diversity metrixs can be called...
	void set_pointers(vector<string>& _obj_names, vector<string>& _obs_obj_names, vector<string>& _obs_obj_sd_names, vector<string>& _pi_obj_names,
                      vector<string>& _pi_obj_sd_names, map<string, double>& _obj_dir_mult,map<string,string>& _ppd_obj_to_sd)
	{
		obj_names_ptr = &_obj_names;
		obs_obj_names_ptr = &_obs_obj_names; 
		obs_obj_sd_names_ptr = &_obs_obj_sd_names;
		pi_obj_names_ptr = &_pi_obj_names; 
		pi_obj_sd_names_ptr = &_pi_obj_sd_names;
		obj_dir_mult_ptr = &_obj_dir_mult;
		prep_pareto_summary_file(POP_SUM_TAG);
		prep_pareto_summary_file(ARC_SUM_TAG);
		prep_pareto_summary_file(ARC_TRIM_SUM_TAG);
        ppd_obj_to_sd_ptr = &_ppd_obj_to_sd;
	}
	
	void update(ObservationEnsemble& oe, ParameterEnsemble& dp, Constraints* constraints_ptr = nullptr);

	bool compare_two(string& first,string& second, MouEnvType envtyp);

	map<string, double> get_spea2_fitness(int generation, ObservationEnsemble& op, ParameterEnsemble& dp, 
		Constraints* constraints_ptr = nullptr, bool report = true, string sum_tag = string());
	map<string, double> get_spea2_kth_nn_crowding_distance(ObservationEnsemble& oe, ParameterEnsemble& dp);
	 
	void get_spea2_archive_names_to_keep(int num_members, vector<string>& keep, const ObservationEnsemble& op, const ParameterEnsemble& dp);

	void prep_pareto_summary_file(string summary_tag);
	void write_pareto_summary(string& sum_tag, int generation, ObservationEnsemble& op, ParameterEnsemble& dp, 
		Constraints* constr_ptr=nullptr);

	//sort specific members
	map<string, double> get_cuboid_crowding_distance(vector<string>& members);
	map<string, double> get_cluster_crowding_fitness(vector<string>& members);
	map<string, double> get_ehvi(vector<string>& members);
	map<string, double> get_mopso_fitness(vector<string> members, ObservationEnsemble& op, ParameterEnsemble& dp);

	double get_ei(map<string, double> phi, string obj, double curr_opt);

	set<string> get_duplicates() { return duplicates;  }

	int get_num_feasible(){ return feas_member_struct.size();}

private:
	
	Pest& pest_scenario;
	FileManager& file_manager;
	PerformanceLog* performance_log;
	//vector<string> obj_names;
	vector<string> sort_members_by_crowding_distance(int front, vector<string>& members, map<string, double>& crowd_map, map<string, map<string, double>>& _member_struct);
	bool first_dominates_second(map<string, double>& first, map<string, double>& second);
	map<string, map<string, double>> get_member_struct(ObservationEnsemble& oe, ParameterEnsemble& dp);
	void drop_duplicates(map<string, map<string, double>>& _member_struct);
	
	bool first_equals_second(map<string, double>& first, map<string, double>& second);

	map<int, vector<string>> sort_members_by_dominance_into_fronts(map<string, map<string, double>>& _member_struct);
	map<string, double> get_mopso_fitness(vector<string> members, map<string, map<string, double>>& _member_struct);
	pair<map<string, double>, map<string, double>> get_spea2_fitness(map<string, map<string, double>>& _member_struct);

	void fill_domination_containers(map<string, map<string, double>>& _member_struct, map<string,
		vector<string>>&solutions_dominated_map, map<string, int>& num_dominating_map, bool dup_as_dom=false);

	bool compare_two_nsga(string& first, string& second);
	bool compare_two_spea(string& first, string& second);

	//sort all members in member struct
	//map<string, double> get_cuboid_crowding_distance();
	map<string, double> get_cuboid_crowding_distance(map<string, map<string, double>>& _member_struct);
	map<string, double> get_cuboid_crowding_distance(vector<string>& members, map<string, map<string, double>>& _member_struct);

	map<string, double> get_spea2_kth_nn_crowding_distance(map<string, map<string, double>>& _member_struct);
	map<string, double> get_spea2_kth_nn_crowding_distance(vector<string>& members, map<string, map<string, double>>& _member_struct);	
	map<string, double> get_cuboid_crowding_distance(ObservationEnsemble& oe, ParameterEnsemble& dp);

	vector<double> get_euclidean_distance(map<string, double> first, map<string, double> second);
	double get_euclidean_fitness(double E, double V);
	map<string, double> get_cluster_crowding_fitness(vector<string>& members, map<string, map<string, double>>& _member_struct);

	map<string, map<string, double>> member_struct;
	vector<string>* obj_names_ptr;
	vector<string>* obs_obj_names_ptr;
	vector<string>* obs_obj_sd_names_ptr;
	vector<string>* pi_obj_names_ptr;
	vector<string>* pi_obj_sd_names_ptr;
	map<string, double>* obj_dir_mult_ptr;
    map<string,string>* ppd_obj_to_sd_ptr;
	set<string> duplicates;

	map<string, map<string, double>> feas_member_struct;
	map<int, vector<string>> front_map;
	map<string, double> crowd_map, fitness_map,  min_sd, nn_map;
	map<string, int> member_front_map;
	map<string, double> member_cvar;
	map<string, double> infeas;
	vector<string> infeas_ordered;
	map<string, double> spea2_constrained_fitness_map;
	map<string, double> spea2_unconstrained_fitness_map;
	

	//PPD-related stuff
	double dominance_probability(map<string, double>& first, map<string, double>& second);
	double dominance_prob_adhoc(map<string, double>& first, map<string, double>& second);
	bool prob_pareto, ppd_sort;
	double ppd_beta;
	vector<double> ppd_range;

	//EHVI-related stuff
	double std_norm_df(double x, double mu, double sd, bool cumdf);
	double psi_function(double aa, double bb, double mu, double sd);
	map<string, double> ehvi_member_map;
	map<string, map<string, double>> incumbent_front_extreme, fit_lookup;
	map<int, vector<double>> hypervolume_partitions;
	double EHVI;
	int iter;
	double get_ehvi(string& member, map<string, map<string, double>>& _member_struct);
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
	MouEnvType envtype;
	MouMateType mattype;
	double epsilon = 1.0e-15;
	Pest& pest_scenario;
	set<string> pp_args;
	vector<MouGenType> gen_types;
	vector<string> act_obs_names, act_par_names;
	int iter, warn_min_members, error_min_members;
	int member_count;
	int archive_size, infill_size;
	string population_dv_file, population_obs_restart_file;
	string dv_pop_file_tag = "dv_pop";
	string pi_pop_file_tag = "pi_pop";
	string obs_pop_file_tag = "obs_pop";
	string lineage_tag = "lineage.csv";
	chancePoints chancepoints;
	FileManager &file_manager; 
	std::mt19937 rand_gen;
	vector<string> obj_names, obs_obj_names, pi_obj_names, obs_obj_sd_names, pi_obj_sd_names;
	vector<string> dv_names;
    map<string,string> ppd_obj_to_sd;
	map<string, double> obj_dir_mult;
	int n_adaptive_dvs;
	map<string, map<string, double>> previous_obj_summary, previous_dv_summary;
	bool risk_obj;
	bool prob_pareto = false; //probabilistic pareto dominance
	bool ppd_sort;
	int restart_iter_offset;
	int save_every;
	map<int,int> population_schedule;
	vector<double> inertia_info, cog_const_range, social_const_range;
	double curr_omega;
	map<string, double> pso_vmax;
	string pso_dv_bound_handling;

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

	map<string,Eigen::VectorXd> par_sim_map, obs_sim_map, pso_velocity_map;

	ParameterEnsemble pso_velocity, pso_pbest_dp;
	ObservationEnsemble pso_pbest_op;

	void update_sim_maps(ParameterEnsemble& _dp, ObservationEnsemble& _op);
	void fill_populations_from_maps(ParameterEnsemble& new_dp, ObservationEnsemble& new_op );

	void update_archive_nsga(ObservationEnsemble& _op, ParameterEnsemble& _dp);
	void update_archive_spea(ObservationEnsemble& _op, ParameterEnsemble& _dp);

	void throw_moea_error(const string& message);

	template<typename T, typename A>
	void message(int level, const string& _message, vector<T, A> _extras, bool echo = true);
	void message(int level, const string& _message);
	template<typename T>
	void message(int level, const string& _message, T extra);

	void sanity_checks();
	vector<int> run_population(ParameterEnsemble& _dp, ObservationEnsemble& _op, bool allow_chance);

	void queue_chance_runs(ParameterEnsemble& _dp);
	ObservationEnsemble get_chance_shifted_op(ParameterEnsemble& _dp, ObservationEnsemble& _op, string& opt_member);

	void initialize_pso();
    ParameterEnsemble get_initial_pso_velocities(int num_members);
    void update_pso_velocity_map(ParameterEnsemble& _pso_velocity);
    void initialize_population_schedule();
	bool initialize_dv_population();
	void initialize_obs_restart_population();

	ParameterEnsemble generate_population();

	ParameterEnsemble generate_diffevol_population(int num_members, ParameterEnsemble& _dp);
	ParameterEnsemble generate_sbx_population(int num_members, ParameterEnsemble& _dp);
	ParameterEnsemble generate_pm_population(int num_members, ParameterEnsemble& _dp);
	ParameterEnsemble generate_pso_population(int num_members, ParameterEnsemble& _dp);
	ParameterEnsemble simplex_cceua_kn(ParameterEnsemble s, int k, int optbounds);																																		
    ParameterEnsemble generate_simplex_population(int num_members, ParameterEnsemble& _dp, ObservationEnsemble& _op);

    ParameterEnsemble generate_empcov_population(int num_members, ParameterEnsemble& _dp, ObservationEnsemble& _op);
	pair<ParameterEnsemble, ParameterEnsemble> get_updated_pso_velocity(ParameterEnsemble& _dp, vector<string>& gbest_solutions);

	vector<string> get_pso_gbest_solutions(int num_reals, ParameterEnsemble& _dp, ObservationEnsemble& _op);
	void update_pso_pbest(ParameterEnsemble& _dp, ObservationEnsemble& _op);

	map<string, string> current_pso_lineage_map, current_empcov_lineage_map;

	vector<int> selection(int num_to_select, ParameterEnsemble& _dp, MouMateType& matetype);

	string get_new_member_name(string tag = string());

	Ensemble save_pi_constraints(ParameterEnsemble &_dp, vector<string> &pinames);
	void save_populations(ParameterEnsemble& _dp, ObservationEnsemble& _op, string tag = string(), bool force_save=false);
	void gauss_mutation_ip(ParameterEnsemble& _dp);
	pair<Eigen::VectorXd, Eigen::VectorXd> sbx(double probability, double eta_m, int idx1, int idx2);
	pair<Eigen::VectorXd, Eigen::VectorXd> sbx_new(double crossover_probability, double di, Eigen::VectorXd& parent1,
		Eigen::VectorXd parent2, vector<string>& _dv_names, Parameters& lbnd, Parameters& ubnd);
	Eigen::VectorXd hybrid_pm(Eigen::VectorXd& parent1, double mutation_probability, double disrupt_probabilty, 
		vector<string>& _dv_names, Parameters& lbnd, Parameters& ubnd);
	pair<double, double> get_betas(double v1, double v2, double distribution_index);
	void get_sbx_child_values(const double& p1, const double& p2, const double& lbnd, 
		const double& ubnd, const double& eta, double& rnd, double& c1, double& c2);


	pair<Parameters, Observations> get_optimal_solution(ParameterEnsemble& _dp, ObservationEnsemble& _oe, string& opt_member_name);

	map<string, map<string, double>> obj_func_report(ParameterEnsemble& _dp, ObservationEnsemble& _oe);
	map<string, map<string, double>> get_obj_func_summary_stats(ParameterEnsemble& _dp, ObservationEnsemble& _op);
	map<string, map<string, double>> obj_func_change_report(map<string, map<string, double>>& current_obj_summary);
    map<string, map<string, double>> decvar_report(ParameterEnsemble& _dp);
    map<string, map<string, double>> decvar_change_report(map<string, map<string, double>>& current_dv_summary);


	int get_max_len_obj_name();
	bool should_use_multigen();

	void queue_resample_runs(ParameterEnsemble& _dp); //outer iters

};

#endif //MOEA_H_

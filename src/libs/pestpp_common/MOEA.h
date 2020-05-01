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
#include "RunStorage.h"
#include "MOUglobal.h"

class Pest;
class RunManagerAbstract;
class RestartController;
class Objectives;
class Constraints;

class MOEA
{
public:
	static mt19937_64 rand_engine;
	MOEA(Pest &_pest_scenario, FileManager &_file_manager, Objectives *_objs, Constraints *_cons,
		const ParamTransformSeq &_par_transform, OutputFileWriter &_output_file_writer,
		PerformanceLog *_performance_log, RunManagerAbstract* _run_mgr_ptr);
	void initialize(RunManagerAbstract &run_manager, int d) throw (MOEAexception);
	void solve(RunManagerAbstract &run_manager, RestartController &restart_controller,
		int max_gen, double f, double cr, bool _dither_f, ModelRun &cur_run);
    // the following is from NSGA2.h
	void advance();
	void evolve();
    // 
    //void set_seed(int s) {
    //    //seed = s;
    //    rgen.set_seed(s);
    //};
    void set_function(individual_config::funcType f) {
        this->function = f;
    };
    void set_popfunction(individual_config::popFuncType f) {
        this->popFunction = f;
    };
    void set_crowdobj(bool crowd) {
        this->crowd_obj = crowd;
    };
    void set_nreal(int nreal) {
        this->nreal = nreal;
    };
    void set_nbin(int nbin) {
        this->nbin = nbin;
    };
    void set_nobj(int nobj) {
        this->nobj = nobj;
    };
    void set_ncon(int ncon) {
        this->ncon = ncon;
    };
    void set_popsize(int popsize) {
        this->popsize = popsize;
    };
    void set_ngen(int ngen) {
        this->ngen = ngen;
    };
    void set_nreport(int nrep) {
        this->nreport = nrep;
    };
    void set_pcross_real(double pcross_real) {
        this->pcross_real = pcross_real;
    };
    void set_pcross_bin(double pcross_bin) {
        this->pcross_bin = pcross_bin;
    };
    void set_pmut_real(double pmut_real) {
        this->pmut_real = pmut_real;
    };
    void set_pmut_bin(double pmut_bin) {
        this->pmut_bin = pmut_bin;
    };
    void set_eta_c(double eta_c) {
        this->eta_c = eta_c;
    };
    void set_eta_m(double eta_m) {
        this->eta_m = eta_m;
    };
    void set_epsilon_c(double epsi_c) {
        this->epsilon_c = epsi_c;
    };
    void set_nbits(const std::vector<int>& nbits) {
        this->nbits = nbits;
    };
    void set_limits_realvar(const std::vector< std::pair<double, double> >& limits_realvar) {
        this->limits_realvar = limits_realvar;
    };
    void set_limits_binvar(const std::vector< std::pair<double, double> >& limits_binvar) {
        this->limits_binvar = limits_binvar;
    };
    void set_backup_filename(const std::string& filename) {
        this->backupFilename = filename;
    }
    void set_custom_report_function(individual_config::popFuncType f) {
        this->reportFunction = f;
    }
    //end copy from NSGA2.h
	~MOEA();
private:
	Pest& pest_scenario;
	const static string solver_type_name;
	FileManager &file_manager;
	Objectives *objs_ptr;
	Constraints *cons_ptr;
	const ParameterInfo *ctl_par_info_ptr;
	const ParameterGroupInfo *par_group_info_ptr;
	ParamTransformSeq par_transform;
	OutputFileWriter &output_file_writer;
	PerformanceLog *performance_log;
	RunManagerAbstract* run_mgr_ptr;
	const ObservationInfo *obs_info_ptr;
	const PriorInformation *prior_info_ptr;
	std::vector<std::string> par_list;
	Parameters max_numeric_pars;
	Parameters min_numeric_pars;
	RunStorage gen_1;
	int best_run_idx;
	int failed_runs_old;
	int failed_runs_new;
	double best_phi;
	double phi_avg_old;
	double phi_avg_new;

	void initialize_vector(Parameters &ctl_pars);
	void mutation(RunManagerAbstract &run_manager, double f, bool dither_f, double cr);
	int recombination(RunManagerAbstract &run_manager);
	void write_run_summary(std::ostream &os,
		int nrun_par, double avg_par, double min_par, double max_par,
		int nrun_can, double avg_can, double min_can, double max_can,
		int nrun_child, double avg_child, double min_child, double max_child);
    // copied from NSGA2.h
    // Parameters to be defined by the user
    int nreal;
    int nbin;
    int nobj;
    int ncon;
    int popsize;
    int ngen;
    int nreport;
    double pcross_real;
    double pcross_bin;
    double pmut_real;
    double pmut_bin;
    double eta_c;
    double eta_m;
    double epsilon_c;
    std::vector<int> nbits;
    std::vector< std::pair<double, double> > limits_realvar;
    // std::vector<double> min_realvar;
    // std::vector<double> max_realvar;
    // double *min_binvar;
    // double *max_binvar;
    std::vector< std::pair<double, double> > limits_binvar;
    // int choice; // to be added later, maybe.
    // int obj1;
    // int obj2;
    // int obj3;
    // int angle1;
    // int angle2;
    individual_config::funcType function;
    individual_config::popFuncType popFunction;
    individual_config::popFuncType reportFunction;
    int t;

    std::string backupFilename;

  public:
    //private:
    void init_streams();
    void report_parameters(std::ostream& os) const;
    void report_pop(const population& pop, std::ostream& os) const;
    void save_backup() const;
    bool load_backup();
    
    void selection(population& oldpop, population& newpop)
        throw (mou::MOEAexception);
    individual& tournament(individual& ind1, individual& ind2) const;
    void crossover(const individual& parent1, const individual& parent2,
        individual& child1, individual& child2);
    void realcross(const individual& parent1, const individual& parent2,
        individual& child1, individual& child2);
    void bincross(const individual& parent1, const individual& parent2,
        individual& child1, individual& child2);
    
    void custom_report(population& pop);
    
    int nbinmut;
    int nrealmut;
    int nbincross;
    int nrealcross;
    int bitlength;
    // random generator?
    std::ofstream fpt1;
    std::ofstream fpt2;
    std::ofstream fpt3;
    std::ofstream fpt4;
    std::ofstream fpt5;
    // FILE *gp;
    population* parent_pop;
    population* child_pop;
    population* mixed_pop;
    bool crowd_obj;
    // end copy
};

#endif //MOEA_H_

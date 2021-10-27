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


class DataAssimilator: public EnsembleMethod
{
public:
	
	using EnsembleMethod::EnsembleMethod;

	void da_update(int cycle);
	void finalize();
    void sanity_checks();
private:
	void eig2csv(string name, Eigen::MatrixXd matrix);	


};


map<int, map<string, double>> process_da_par_cycle_table(Pest& pest_scenario, vector <int>& ncycles_in_tables, ofstream& fout_rec);
map<int, map<string, double>> process_da_obs_cycle_table(Pest& pest_scenario, vector <int>& ncycles_in_tables, ofstream& fout_rec, set<string>& obs_in_tbl);
map<int, map<string, double>> process_da_weight_cycle_table(Pest& pest_scenario, vector <int>& ncycles_in_tables, ofstream& fout_rec, set<string>& obs_in_tbl);

void write_global_phi_info(int cycle, ofstream& f_phi, DataAssimilator& da, vector<string>& init_real_names);

void generate_global_ensembles(DataAssimilator& da, ofstream& fout_rec, ParameterEnsemble& curr_pe, 
	ObservationEnsemble& curr_oe, ObservationEnsemble& curr_noise);



#endif

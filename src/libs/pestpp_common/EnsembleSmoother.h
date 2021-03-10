#ifndef ENSEMBLESMOOTHER_H_
#define ENSEMBLESMOOTHER_H_

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




class IterEnsembleSmoother: public EnsembleMethod
{
public:

	using EnsembleMethod::EnsembleMethod;
	
	
	//void initialize();
	void iterate_2_solution();
	void finalize();
	void throw_ies_error(string message);
	//bool should_terminate();
	

private:
	
	bool use_mda;
	vector<double> mda_facs;

	
	void sanity_checks();

	//void add_bases();

	//void update_reals_by_phi(ParameterEnsemble &_pe, ObservationEnsemble &_oe);

	//vector<string> detect_prior_data_conflict();

	//void set_subset_idx(int size);
	//Eigen::MatrixXd get_Am(const vector<string> &real_names, const vector<string> &par_names);

	//void zero_weight_obs(vector<string>& obs_to_zero_weight, bool update_obscov=true,bool update_oe_base=true);

	//void norm_map_report(map<string, double>& norm_map, string tag, double thres = 0.1);

};

#endif

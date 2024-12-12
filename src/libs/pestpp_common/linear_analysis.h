#ifndef LINEAR_ANALYSIS_H_
#define LINEAR_ANALYSIS_H_

//#include <vector>
//#include <string>

#include "Pest.h"
#include "covariance.h"
#include "logger.h"
#include "ObjectiveFunc.h"
#include "ModelRunPP.h"
#include "OutputFileWriter.h"
#include "PerformanceLog.h"
#include "RunManagerAbstract.h"
#include "Ensemble.h"

using namespace std;
class LinearAnalysis
{
public:
	//empty constructor
	//linear_analysis():pest_scenario(Pest()),file_manager(FileManager()),jacobian(Mat()){;}

	//constructor for pest++ integration
	LinearAnalysis(Mat &_jacobian, Pest &_pest_scenario, FileManager& _file_manager, PerformanceLog &pfm, Covariance& _parcov, std::mt19937* _rand_gen_ptr);
	//linear_analysis(Mat* _jacobian, Pest* pest_scenario, Mat* _obscov, Logger* _log = new Logger());


	//directly from Mat objects
	//linear_analysis(Mat& _jacobian,Pest& _pest_scenario, FileManager& _file_manager, Mat& _parcov, Mat& _obscov, map<string, Mat> _predictions,Logger* _log = new Logger());

    ObservationInfo glm_iter_fosm(ModelRun& optimum_run, OutputFileWriter& output_file_writer, int iter, 
		RunManagerAbstract* run_mgr_ptr);
	pair<ParameterEnsemble,map<int,int>> draw_fosm_reals(RunManagerAbstract* run_mgr_ptr, int iter,ModelRun& optimum_run);
	pair<ObservationEnsemble,map<string,double>> process_fosm_reals(RunManagerAbstract* run_mgr_ptr, pair<ParameterEnsemble, map<int, int>>& fosm_real_info, int iter,
		double last_best_phi);

	void set_predictions(vector<string> preds,bool forgive=false);
	

	void  set_parcov(Mat &_parcov);
	void set_obscov(Mat &_obscov) { obscov = _obscov; }

	//get a new linear analysis object consisting of a subset of par and obs names
	//LinearAnalysis get(vector<string> &new_par_names, vector<string> &new_obs_names);

	//exposed schur functionality
	//from the diagonal of parcov
	double prior_parameter_variance(string &par_name);
	//a map of <par_name,variance> from diagonal of parcov
	map<string, double> prior_parameter_variance();
	Mat prior_parameter_matrix() { return parcov; }

	//from the diagonal of schur's complement
	double posterior_parameter_variance(string &par_name);
	//a map of <par_name,variance> from the diagonal of schur's complement
	map<string, double> posterior_parameter_variance();
	//the full matrix
	Mat posterior_parameter_matrix();
	Mat* posterior_parameter_ptr();
	Covariance posterior_parameter_covariance_matrix();
	//prior predictive variance from parcov
	double prior_prediction_variance(string &pred_name);
	//map <pred_name,variance> from parcov
	map<string, double> prior_prediction_variance();

	//posterior predictive variance from schur's complement
	double posterior_prediction_variance(string &pred_name);
	//map <pred_name,variance> from schur's complement
	map<string, double> posterior_prediction_variance();

	Mat get_jacobian(){ return jacobian; }
	Mat get_omitted_jacobian(){ return omitted_jacobian; }
	Covariance get_parcov(){ return parcov; }
	Covariance get_omitted_parcov(){ return omitted_parcov; }
	Covariance get_obscov(){ return obscov; }

	Mat get_S(){ return S; }
	Mat get_V(){ return V; }

	Mat* get_jacobian_ptr(){ return &jacobian; }
	Covariance* get_parcov_ptr(){ return &parcov; }
	Covariance* get_obscov_ptr(){ return &obscov; }
	Mat* get_omitted_jacobian_ptr(){ return &omitted_jacobian; }
	Covariance* get_omitted_parcov_ptr(){ return &omitted_parcov; }

	map<string,Mat> get_predictions(){ return predictions; }
	map<string,Mat> get_omitted_predictions(){ return omitted_predictions; }

	//aligns the different linear components
	void align();

	map<string, double> like_preds(double val);

	~LinearAnalysis();

	//some convenience methods for PEST++ integration
	void write_par_credible_range(ofstream &fout, string sum_filename, ParameterInfo parinfo,
		Parameters init_pars, Parameters opt_pars,vector<string> ordered_names);
	void write_pred_credible_range(ofstream &fout, string sum_filename, map<string,pair<double,double>> init_final_pred_values);
	void drop_prior_information(const Pest &pest_scenario);

private:

	//Logger* log;
	std::mt19937* rand_gen_ptr;
	PerformanceLog& pfm;
	FileManager& file_manager;
	Pest& pest_scenario;
	Mat& jacobian;
	Mat S, V;
	Mat normal;
	Mat R, G, ImR, V1;
	int R_sv, G_sv, ImR_sv,V1_sv;

	Mat omitted_jacobian;
	Covariance parcov;
	Covariance obscov;
	Covariance posterior;
	map<string, Mat> predictions;
	map<string,Mat> omitted_predictions;
	Covariance omitted_parcov;

	void calc_posterior();
	void svd();
	void build_normal();
	void build_R(int sv);
	void build_G(int sv);
	void build_ImR(int sv);
	void build_V1(int sv);


	Covariance condition_on(vector<string> &keep_par_names,vector<string> &cond_par_names);

	void throw_error(const string &message);
	void load_jco(Mat &jco, const string &jco_filename);
	void load_pst(Pest &pest_scenario, const string &pst_filename);
	void load_parcov(const string &parcov_filename);
	void load_obscov(const string &obscov_filename);

	//scale the jacobian by parcov
	void kl_scale();

	pair<double, double> get_range(double value, double variance, const ParameterRec::TRAN_TYPE &tt);









};

map<string, double> get_obj_comps(string &filename);
map<string, int> get_nnz_group(Pest &pest_scenario);
ObservationInfo normalize_weights_by_residual(Pest &pest_scenario, string &resid_filename);
ObservationInfo normalize_weights_by_residual(Pest &pest_scenario, PhiData obj);
ObservationInfo normalize_weights_by_residual(Pest &pest_scenario, Observations &sim);
#endif

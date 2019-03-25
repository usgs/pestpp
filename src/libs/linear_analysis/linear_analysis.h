#ifndef LINEAR_ANALYSIS_H_
#define LINEAR_ANALYSIS_H_

//#include <vector>
//#include <string>

#include "Pest.h"
#include "covariance.h"
#include "logger.h"
#include "ObjectiveFunc.h"

using namespace std;
class linear_analysis
{
public:
	//empty constructor
	linear_analysis(){;}

	//the easiest constructor, builds parcov and obscov from the pst associated with jco_filename
	linear_analysis(string &jco_filename,Logger* _log = new Logger());

	//loads parcov and obscov from files, can .pst, .mat or .unc files
	linear_analysis(string &jco_filename, string &parcov_filename, string &obscov_filename, Logger* _log = new Logger());

	//load parcov and obscov from parameter bounds and observation weights
	linear_analysis(Mat _jacobian, Pest pest_scenario, Logger* _log = new Logger());

	//pointer constructors for pest++ integration
	linear_analysis(Mat* _jacobian, Pest* pest_scenario, Logger* _log = new Logger());
	linear_analysis(Mat* _jacobian, Pest* pest_scenario, Mat* _obscov, Logger* _log = new Logger());


	//directly from Mat objects
	linear_analysis(Mat _jacobian, Mat _parcov, Mat _obscov, map<string, Mat> _predictions,Logger* _log = new Logger());

	void set_predictions(vector<string> preds);
	void set_predictions(vector<Mat> preds);
	void set_predictions(Mat* preds);

	void  set_parcov(Mat* _parcov);

	//get a new linear analysis object consisting of a subset of par and obs names
	linear_analysis get(vector<string> &new_par_names, vector<string> &new_obs_names);

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

	//the reduction in predictive variance from some obs
	double posterior_predictive_worth(string &pred_name, vector<string> &obs_names);
	//<pred_name,variance_reduction> from some obs
	map<string, double> worth(vector<string> &obs_names);

	//reduction in prior and posterior predictive variance from perfect knowledge of some pars
	//<pred_name,prior and posterior variance_reduction> from perfect knowledge of some pars
	map<string, pair<double, double>> contribution(vector<string> &par_names);



	//exposed error variance functionality

	//<err_var_component("null","solution","omitted"),error_variance> for a set of singular values and a parameter
	map<string, vector<double>> parameter_error_variance_components(vector<int> sing_vals, string &par_name);
	//<err_var_component("null","solution","omitted"),error_variance> for a set of singular values and a prediction
	map<string, vector<double>> prediction_error_variance_components(vector<int> sing_vals, string &pred_name);

	//extract elements from the jacobian, parcov, and predictions and set them as omitted
	void extract_omitted(vector<string> &omitted_par_names);
	void extract_omitted(string &omitted_par_name);

	Covariance first_parameter(int sv);
	Covariance second_parameter(int sv);
	Covariance third_parameter(int sv);

	map<string, double> first_prediction(int sv);
	map<string, double> second_prediction(int sv);
	map<string, double> third_prediction(int sv);

	//other stuff
	//<par_names,ident>
	map<string, double> parameter_ident(int sv);

	//<singular_value,vector<sup_obs>>
	map<int, vector<double>> super_obs(vector<int> sing_vals);
	map<int, vector<double>> super_obs(vector<int> sing_vals, vector<string> &sup_obs_names);

	//<singular_value,vector<sup_par>>
	map<int, vector<double>> super_par(vector<int> sing_vals);
	map<int, vector<double>> super_par(vector<int> sing_vals, vector<string> &sup_par_names);

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

	Mat* get_normal_ptr();
	Mat* get_S_ptr(int sv);
	Mat* get_V_ptr(int sv);

	Mat* get_R_ptr(int sv);
	Mat* get_G_ptr(int sv);
	Mat* get_ImR_ptr(int sv);
	Mat* get_V1_ptr(int sv);

	map<string,Mat> get_predictions(){ return predictions; }
	map<string,Mat> get_omitted_predictions(){ return omitted_predictions; }

	//aligns the different linear components
	void align();

	map<string, double> like_preds(double val);

	~linear_analysis();

	//some convience methods for PEST++ integration
	void write_par_credible_range(ofstream &fout, string sum_filename, ParameterInfo parinfo,
		Parameters init_pars, Parameters opt_pars,vector<string> ordered_names);
	void write_pred_credible_range(ofstream &fout, string sum_filename, map<string,pair<double,double>> init_final_pred_values);
	void drop_prior_information(const Pest &pest_scenario);

private:

	Logger* log;
	Mat jacobian;
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

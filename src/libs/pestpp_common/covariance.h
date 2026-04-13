#ifndef COVARIANCE_H_
#define COVARIANCE_H_
#include <string>
#include <sstream>
#include <vector>
#include <random>
#include<Eigen/Sparse>

#include "Pest.h"
#include "logger.h"
#include "FileManager.h"

using namespace std;

class Mat
{
public:
	enum class MatType{ DIAGONAL, SPARSE, DENSE };
	Mat(){ autoalign = true; };
	Mat(string filename);

	Mat(vector<string> _row_names, vector<string> _col_names,
		Eigen::SparseMatrix<double> _matrix);

	Mat(vector<string> _row_names, vector<string> _col_names,
		Eigen::SparseMatrix<double>* _matrix);

	Mat(vector<string> _row_names, vector<string> _col_names,
		Eigen::SparseMatrix<double> _matrix, MatType _matttype);

	vector<string> get_row_names(){ return row_names; }
	vector<string> get_col_names(){ return col_names; }
	void clear_names(){row_names.clear();col_names.clear();}
	const vector<string>* rn_ptr();
	const vector<string>* cn_ptr();
	Eigen::SparseMatrix<double> get_matrix(){ return matrix; }

	const Eigen::SparseMatrix<double>* e_ptr();
	const Eigen::SparseMatrix<double>* U_ptr();
	const Eigen::SparseMatrix<double>* V_ptr();
	const Eigen::VectorXd* s_ptr();

	Mat get_U();
	Mat get_V();
	Mat get_s();

	MatType get_mattype(){ return mattype; }

	void to_ascii(const string &filename);
	void from_ascii(const string &filename);
	void to_binary(const string &filename);
	void to_binary_new(const string &filename);
	void from_binary(const string &filename);

	void from_csv(const string &filename);
	//void to_csv(const string filename);

	void from_file(const string &filename);

	void from_triplets(const vector<string> &_row_names, const vector<string> &_col_names, const vector<Eigen::Triplet<double>> &triplets);

	void transpose_ip();
	Mat transpose();
	Mat T();
	Mat inv(Logger* log);
	Mat inv(bool echo=false);
	void inv_ip(PerformanceLog& pfm);
	void inv_ip(bool echo=false);
	void SVD();

	Mat identity();
	Mat zero();


	Mat get(const vector<string> &other_row_names, const vector<string> &other_col_names,bool update=true);
	Mat leftCols(const int idx);
	Mat rightCols(const int idx);
	Mat extract(const vector<string> &extract_row_names, const vector<string> &extract_col_names);
	Mat extract(const string &extract_row_name, const vector<string> &extract_col_names);
	Mat extract(const vector<string> &extract_row_names, const string &extract_col_name);
	void drop_rows(const vector<string> &drop_row_names);
	void drop_cols(const vector<string> &drop_col_names);

	int nrow(){ return row_names.size(); }
	int ncol(){ return col_names.size(); }

	bool isdiagonal();

	void update_sets();

	//void set_row_names(vector<string> _row_names) { row_names = _row_names; }
	//void set_col_names(vector<string> _col_names) { col_names = _col_names; }
	

protected:
	bool autoalign;
	Eigen::SparseMatrix<double> matrix;
	Eigen::SparseMatrix<double> U;
	Eigen::SparseMatrix<double> V;
	Eigen::VectorXd s;
	vector<string> row_names;
	vector<string> col_names;
	set<string> row_set;
	set<string> col_set;
	int icode = 2;
	MatType mattype;

	vector<string> read_namelist(ifstream &in, int &nitems);
};


class Covariance : public Mat
{
public:
	Covariance(vector<string> &names);
	Covariance();
	Covariance(string filename);
	Covariance(Mat _mat);
	Covariance(vector<string> _row_names, Eigen::SparseMatrix<double> _matrix, Mat::MatType _mattype = Mat::MatType::SPARSE);
	
	Covariance get(const vector<string> &other_names, bool update=true);
	double get(const string& row_name, const string& col_name);
	Mat get(vector<string> &other_row_names, vector<string> &other_col_names){ return Mat::get(other_row_names, other_col_names); }
	void drop(vector<string> &drop_names);
	Covariance extract(vector<string> &extract_names);

	Covariance diagonal(double val);
	void from_diagonal(Covariance &other);

	string try_from(Pest &pest_scenario, FileManager &file_manager,bool is_parcov=true, bool forgive_missing=false);
	void from_uncertainty_file(const string &filename, vector<string> &ordered_names);
	void from_parameter_bounds(Pest &pest_scenario, ofstream& frec);
	void from_parameter_bounds(ofstream& frec, const vector<string> &par_names, const ParameterInfo &par_info, map<string,
		double>& par_std,double sigma_range=4.0);

	void from_observation_weights(Pest &pest_scenario, ofstream& frec);
	void from_observation_weights(ofstream& frec, const vector<string>& obs_names, const ObservationInfo& obs_info,
		const vector<string>& pi_names, const PriorInformation* pi,map<string,double>& obs_std);
	void from_parameter_bounds(const string &pst_filename, ofstream& frec);
	void from_observation_weights(const string &pst_filename, ofstream& frec);

	void to_uncertainty_file(const string &filename);

	vector<Eigen::VectorXd> draw(int ndraws);
	vector<double> standard_normal(default_random_engine gen);
	void cholesky();


private:
	Eigen::SparseMatrix<double> lower_cholesky;
};

ostream& operator<< (std::ostream &os, Mat mat);
#endif

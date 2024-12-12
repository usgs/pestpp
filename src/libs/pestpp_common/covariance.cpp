
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <random>
#include <iterator>
#include <Eigen/Dense>

#include <Eigen/SparseCholesky>

#include "SVDPackage.h"
#include "Pest.h"
#include "utilities.h"
#include "covariance.h"
#include "FileManager.h"

using namespace std;

//---------------------------------------
//Mat constructors
//---------------------------------------

Mat::Mat(string filename)
{
	string ext = filename.substr(filename.find_last_of(".") + 1);
	pest_utils::upper_ip(ext);
	if (ext == "PST")
		throw runtime_error("Mat::Mat() error: cannot instantiate a Mat with PST");
	else if ((ext == "JCO") || (ext == "JCB"))
		from_binary(filename);
	else if (ext == "MAT")
		from_ascii(filename);
	else
		throw runtime_error("Mat::Mat() error: only .jco/.jcb or .mat\
							 files can be used to instantiate a Mat");

}

Mat::Mat(vector<string> _row_names, vector<string> _col_names,
	Eigen::SparseMatrix<double> _matrix)
{
	row_names = _row_names;
	col_names = _col_names;
	assert(row_names.size() == _matrix.rows());
	assert(col_names.size() == _matrix.cols());
	matrix = _matrix;
	mattype = MatType::SPARSE;
}

Mat::Mat(vector<string> _row_names, vector<string> _col_names,
	Eigen::SparseMatrix<double>* _matrix)
{
	row_names = _row_names;
	col_names = _col_names;
	assert(row_names.size() == _matrix->rows());
	assert(col_names.size() == _matrix->cols());
	matrix = *_matrix;
	mattype = MatType::SPARSE;
}

Mat::Mat(vector<string> _row_names, vector<string> _col_names,
	Eigen::SparseMatrix<double> _matrix,MatType _mattype)
{
	row_names = _row_names;
	col_names = _col_names;
	assert(row_names.size() == _matrix.rows());
	assert(col_names.size() == _matrix.cols());
	matrix = _matrix;
	mattype = _mattype;
}

void Mat::update_sets()
{
	row_set.clear();
	//row_set.emplace(row_names.begin(), row_names.end());
	row_set = set<string>(row_names.begin(), row_names.end());
	col_set.clear();
	//col_set.emplace(col_names.begin(), col_names.end());
	col_set = set<string>(col_names.begin(), col_names.end());
}


const Eigen::SparseMatrix<double>* Mat::e_ptr()
{
	const Eigen::SparseMatrix<double>* ptr = &matrix;
	return ptr;
}

const vector<string>* Mat::rn_ptr()
{
	const vector<string>* ptr = &row_names;
	return ptr;
}

const vector<string>* Mat::cn_ptr()
{
	const vector<string>* ptr = &col_names;
	return ptr;
}

Mat Mat::identity()
{
	Eigen::SparseMatrix<double> i(nrow(), ncol());
	i.setZero();
	i.setIdentity();
	return Mat(*rn_ptr(),*cn_ptr(),i);
}

Mat Mat::zero()
{
	Eigen::SparseMatrix<double> i(nrow(), ncol());
	i.setZero();
	return Mat(*rn_ptr(), *cn_ptr(), i);
}


const Eigen::SparseMatrix<double>* Mat::U_ptr()
{
	if (U.rows() == 0)
	{
		SVD();
	}
	const Eigen::SparseMatrix<double>* ptr = &U;
	return ptr;
}

const Eigen::SparseMatrix<double>* Mat::V_ptr()
{
	if (V.rows() == 0)
	{
		SVD();
	}
	const Eigen::SparseMatrix<double>* ptr = &V;
	return ptr;
}

const Eigen::VectorXd* Mat::s_ptr()
{
	if (s.size() == 0)
	{
		SVD();
	}
	const Eigen::VectorXd* ptr = &s;
	return ptr;
}

Mat Mat::get_U()
{
	if (U.rows() == 0) SVD();
	vector<string> u_col_names;
	stringstream ss;
	for (int i = 0; i < nrow(); i++)
	{
		ss.clear();
		ss.str(string());
		ss << "left_sing_vec_";
		ss << i + 1;
		u_col_names.push_back(ss.str());
	}
	return Mat(row_names, u_col_names, U);
}

Mat Mat::get_V()
{
	if (V.rows() == 0) SVD();
	vector<string> v_col_names;
	stringstream ss;
	for (int i = 0; i < ncol(); i++)
	{
		ss.clear();
		ss.str(string());
		ss << "right_sing_vec_";
		ss << i + 1;
		v_col_names.push_back(ss.str());
	}
	return Mat(col_names, v_col_names, V);
}


Mat Mat::get_s()
{
	if (V.rows() == 0) SVD();
	vector<string> s_names;
	vector<Eigen::Triplet<double>> triplet_list;
	stringstream ss;
	for (int i = 0; i < s.size(); i++)
	{
		ss.clear();
		ss.str(string());
		ss << "sing_val_";
		ss << i + 1;
		s_names.push_back(ss.str());
		triplet_list.push_back(Eigen::Triplet<double>(i, i, s[i]));
	}
	Eigen::SparseMatrix<double> s_mat(s.size(),s.size());
	s_mat.setZero();
	s_mat.setFromTriplets(triplet_list.begin(), triplet_list.end());
	return Mat(s_names, s_names, s_mat);
}

Mat Mat::transpose()
{
	return Mat(col_names, row_names, matrix.transpose());
}

Mat Mat::T()
{
	return Mat(col_names, row_names, matrix.transpose());
}

void Mat::transpose_ip()
{
	if (mattype != MatType::DIAGONAL)
	{
		matrix = matrix.transpose();
		vector<string> temp = row_names;
		row_names = col_names;
		col_names = temp;
	}
}

Mat Mat::inv(bool echo)
{
	Logger* log = new Logger();
	log->set_echo(echo);
	//inv_ip(log);
	Mat new_mat = inv(log);
	delete log;
	return new_mat;
}


Mat Mat::inv(Logger* log)
{
	if (nrow() != ncol()) throw runtime_error("Mat::inv() error: only symmetric positive definite matrices can be inverted with Mat::inv()");
	if (mattype == MatType::DIAGONAL)
	{
		log->log("inverting diagonal matrix in place");
		log->log("extracting diagonal");
		Eigen::VectorXd diag = matrix.diagonal().eval();
		log->log("inverting diagonal");
		log->log("building triplets");
		vector<Eigen::Triplet<double>> triplet_list;
		for (int i = 0; i != diag.size(); ++i)
		{
			triplet_list.push_back(Eigen::Triplet<double>(i, i, 1.0 / diag[i]));
		}
		Eigen::SparseMatrix<double> inv_mat;
		inv_mat.conservativeResize(triplet_list.size(), triplet_list.size());
		inv_mat.setZero();
		log->log("setting matrix from triplets");
		inv_mat.setFromTriplets(triplet_list.begin(), triplet_list.end());
		return Mat(row_names, col_names, inv_mat);
	}
	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	solver.compute(matrix);
	Eigen::SparseMatrix<double> I(nrow(), nrow());
	I.setIdentity();
	Eigen::SparseMatrix<double> inv_mat = solver.solve(I);
	return Mat(row_names, col_names, inv_mat);
}

void Mat::inv_ip(bool echo)
{
	ofstream flog("Mat.log");
	PerformanceLog pfm(flog);
	inv_ip(pfm);
	return;
}

void Mat::inv_ip(PerformanceLog& pfm)
{
	if (nrow() != ncol()) throw runtime_error("Mat::inv() error: only symmetric positive definite matrices can be inverted with Mat::inv()");
	if (mattype == MatType::DIAGONAL)
	{
		pfm.log_event("inverting diagonal matrix in place");
		Eigen::VectorXd diag = matrix.diagonal().eval();
		vector<Eigen::Triplet<double>> triplet_list;
		for (int i = 0; i != diag.size(); ++i)
		{
			triplet_list.push_back(Eigen::Triplet<double>(i, i, 1.0/diag[i]));
		}
		//log->log("resizeing matrix to size " + triplet_list.size());
		//matrix.conservativeResize(triplet_list.size(),triplet_list.size());
		matrix.setZero();
		matrix.setFromTriplets(triplet_list.begin(),triplet_list.end());
		//Eigen::SparseMatrix<double> inv_mat(triplet_list.size(), triplet_list.size());
		//inv_mat.setZero();
		//inv_mat.setFromTriplets(triplet_list.begin(), triplet_list.end());
		//matrix = inv_mat;
		//cout << "diagonal inv_ip()" << endl;
		/*matrix.resize(triplet_list.size(), triplet_list.size());
		matrix.setZero();
		matrix.setFromTriplets(triplet_list.begin(),triplet_list.end());*/
		return;
	}

	//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
	pfm.log_event("inverting non-diagonal matrix in place");
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	solver.compute(matrix);
	Eigen::SparseMatrix<double> I(nrow(), nrow());
	I.setIdentity();
	//Eigen::SparseMatrix<double> inv_mat = solver.solve(I);
	//matrix = inv_mat;
	matrix = solver.solve(I);
	//cout << "full inv_ip()" << endl;
	/*matrix.setZero();
	matrix = solver.solve(I);*/

}


void Mat::SVD()
{
	Eigen::JacobiSVD<Eigen::MatrixXd> svd_fac(matrix, Eigen::DecompositionOptions::ComputeFullU |
		Eigen::DecompositionOptions::ComputeFullV);
	s = svd_fac.singularValues();
	U = svd_fac.matrixU().sparseView();
	V = svd_fac.matrixV().sparseView();
}



//---------------------------------------
//Mat operator
//--------------------------------------

ostream& operator<< (ostream &os, Mat mat)
{
	cout << "row names : ";
	for (auto &name : mat.get_row_names())
		cout << name << ',';
	cout << endl << "col names : ";
	for (auto &name : mat.get_col_names())
		cout << name << ',';
	cout << endl;
	cout << *mat.e_ptr();
	return os;
}




//-----------------------------------------
//Mat IO
//-----------------------------------------




void Mat::to_ascii(const string &filename)
{
	ofstream out(filename);
	if (!out.good())
	{
		throw runtime_error("Mat::to_ascii() error: cannot open " + filename + "\
													 to write ASCII matrix");
	}
	out << setw(6) << nrow() << setw(6) << ncol() << setw(6) << icode << endl;
	out << matrix.toDense() << endl;
	if (icode == 1)
	{
		out<< "* row and column names" << endl;
		for (auto &name : row_names)
			out << pest_utils::lower_cp(name) << endl;

	}
	else
	{
		out << "* row names" << endl;
		for (auto &name : row_names)
			out << pest_utils::lower_cp(name) << endl;
		out << "* column names" << endl;
		for (auto &name : col_names)
			out << pest_utils::lower_cp(name) << endl;
	}
	out.close();
}

void Mat::from_file(const string &filename)
{
	stringstream ss;
	string ext = filename.substr(filename.find_last_of(".") + 1);
	pest_utils::upper_ip(ext);
	if ((ext == "JCB") || (ext == "JCO"))
	{
		from_binary(filename);
	}
	else if ((ext == "MAT") || (ext == "COV"))
	{
		from_ascii(filename);
	}
	else if (ext == "CSV")
	{
		from_csv(filename);
	}
	else
	{
		ss << "Mat::from_file() error: unrecognized extension'" << ext << "', should be JCB, JCO, MAT or CSV";
		throw runtime_error(ss.str());
	}
	mattype = MatType::SPARSE;

}

void Mat::from_csv(const string &filename)
{
	ifstream csv(filename);
	if (!csv.good())
		throw runtime_error("Mat::from_csv() error: cannot open " + filename + " \
												to read csv matrix");


	//process the header
	//any missing header labels will be marked to ignore those columns later
	string line;
	if (!getline(csv, line))
		throw runtime_error("error reading header (first) line from csv file :");
	pest_utils::strip_ip(line);
	pest_utils::upper_ip(line);


	pest_utils::tokenize(line, col_names, ",", false);
	col_names.erase(col_names.begin()); //drop the index label
	vector<Eigen::Triplet<double>> triplet_list;
	double val;

	//read a csv file to an Ensmeble
	int lcount = 0, irow = 0;
	//vector<vector<double>> vectors;
	vector<string> tokens;
	string row_name;

	while (getline(csv, line))
	{
		pest_utils::strip_ip(line);
		tokens.clear();
		pest_utils::tokenize(line, tokens, ",", false);
		if (tokens[tokens.size() - 1].size() == 0)
			tokens.pop_back();

		try
		{
			pest_utils::convert_ip(tokens[0], row_name);
		}
		catch (exception &e)
		{
			stringstream ss;
			ss << "error converting token '" << tokens[0] << "' to <int> run_id on line " << lcount << ": " << line << endl << e.what();
			throw runtime_error(ss.str());
		}
		tokens.erase(tokens.begin()); //drop the row name
		if (tokens.size() != col_names.size())
		{
			stringstream ss;
			ss << "Matrix.from_csv() error: wrong number of entries on line " << lcount << " , expecting " << col_names.size() << ", found " << tokens.size() << endl;
			throw runtime_error(ss.str());
		}
		row_names.push_back(pest_utils::upper_cp(row_name));
		for (int j = 0; j < col_names.size(); j++)
		{
			try
			{
				val = pest_utils::convert_cp<double>(tokens[j]);
			}
			catch (exception &e)
			{
				stringstream ss;
				ss << "error converting token '" << tokens[j] << "' to double for " << col_names[j] << " on line " << lcount << " : " << e.what();
				throw runtime_error(ss.str());
			}
			if (val != 0.0)
				triplet_list.push_back(Eigen::Triplet<double>(irow, j, val));
		}

	lcount++;
	irow++;
	}
	matrix.resize(row_names.size(), col_names.size());
	matrix.setZero();
	matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());

}


void Mat::from_triplets(const vector<string> &_row_names, const vector<string> &_col_names, const vector<Eigen::Triplet<double>> &triplets)
{
	matrix.resize(_row_names.size(), _col_names.size());
	//matrix.setZero();
	matrix.setFromTriplets(triplets.begin(), triplets.end());
	//cout << matrix.cols() << endl;
	row_names = _row_names;
	col_names = _col_names;

}

void Mat::from_ascii(const string &filename)
{
	ifstream in(filename);
	if (!in.good())
		throw runtime_error("Mat::from_ascii() error: cannot open " + filename + " \
												to read ASCII matrix");
	int nrow = -999, ncol = -999;
	if (in >> nrow >> ncol >> icode){}
	else
		throw runtime_error("Mat::from_ascii() error reading nrow ncol icode from first line\
							 of ASCII matrix file: " + filename);

	vector<Eigen::Triplet<double>> triplet_list;
	double val;
	int irow = 0, jcol = 0;
	for (int inode = 0; inode < nrow*ncol;inode++)
	{
		if (in >> val)
		{
			if (val != 0.0)
				triplet_list.push_back(Eigen::Triplet<double>(irow,jcol,val));
			jcol++;
			if (jcol >= ncol)
			{
				irow++;
				jcol = 0;
			}
		}
		else
		{
			string i_str = to_string(inode);
			throw runtime_error("Mat::from_ascii() error reading entry number "+i_str+" from\
								 ASCII matrix file: "+filename);
		}
	}

	string header;
	//read the newline char
	getline(in, header);
	if (!getline(in,header))
		throw runtime_error("Mat::from_ascii() error reading row/col description\
							 line from ASCII matrix file: " + filename);
	pest_utils::upper_ip(header);
	string name;
	if (icode == 1)
	{
		if (nrow != ncol)
			throw runtime_error("Mat::from_ascii() error: nrow != ncol for icode type 1 ASCII matrix file:" + filename);
		if((header.find("ROW") == string::npos) || (header.find("COLUMN") == string::npos))
			throw runtime_error("Mat::from_ascii() error: expecting row and column names header instead\
								 of:" + header + " in ASCII matrix file: " + filename);
		try
		{
			row_names = read_namelist(in, nrow);
		}
		catch (exception &e)
		{
			throw runtime_error("Mat::from_ascii() error reading row/column names from ASCII matrix file: " + filename + "\n" + e.what());
		}
		if ((nrow != row_names.size()) || (ncol != row_names.size()))
			throw runtime_error("Mat::from_ascii() error: number of row/col names does not match matrix dimensions");
		col_names = row_names;
	}
	else
	{
		if(header.find("ROW") == string::npos)
			throw runtime_error("Mat::from_ascii() error: expecting row names header instead of:" + header + " in ASCII matrix file: " + filename);
		try
		{
			row_names = read_namelist(in, nrow);
		}
		catch (exception &e)
		{
			throw runtime_error("Mat::from_ascii() error reading row names from ASCII matrix file: " + filename + "\n" + e.what());
		}
		if (!getline(in, header))
		{
			throw runtime_error("Mat::from_ascii() error reading column name descriptor from ASCII matrix file: " + filename);
		}
		pest_utils::upper_ip(header);
		if (header.find("COLUMN") == string::npos)
			throw runtime_error("Mat::from_ascii() error: expecting column names header instead of:" + header + " in ASCII matrix file: " + filename);
		try
		{
			col_names = read_namelist(in, ncol);
		}
		catch (exception &e)
		{
			throw runtime_error("Mat::from_ascii() error reading column names from ASCII matrix file: " + filename + "\n" + e.what());
		}
		if (nrow != row_names.size())
			throw runtime_error("Mat::from_ascii() error: nrow != row_names.size() in ASCII matrix file: " + filename);

		if(ncol != col_names.size())
			throw runtime_error("Mat::from_ascii() error: ncol != col_names.size() in ASCII matrix file: " + filename);
	}
	in.close();

	Eigen::SparseMatrix<double> new_matrix(nrow, ncol);
	new_matrix.setZero();  // initialize all entries to 0
	new_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
	
	matrix = new_matrix;
}

vector<string> Mat::read_namelist(ifstream &in, int &nitems)
{
	vector<string> names;
	string name;
	for (int i = 0; i < nitems; i++)
	{
		if (!getline(in, name))
		{
			string i_str = to_string(i);
			throw runtime_error("Mat::read_namelist() error reading name for entry " + i_str);
		}
		if (name.find("*") != string::npos)
		{
			string i_str = to_string(i);
			throw runtime_error("Mat::read_namelist() error: '*' found in item name: " + name+", item number: "+i_str);
		}
		pest_utils::strip_ip(name);
		pest_utils::upper_ip(name);
		if (find(names.begin(), names.end(), name) != names.end())
			throw runtime_error("Mat::read_namelist() error: duplicate name: " + name + " found in name list");
		names.push_back(name);
	}
	return names;
}

void Mat::to_binary(const string &filename)
{
	pest_utils::save_binary_orgfmt(filename, row_names, col_names, matrix);
//	ofstream jout(filename, ios::out | ios::binary);
//	int n_par = col_names.size();
//	int n_obs_and_pi = row_names.size();
//	int n;
//	int tmp;
//	double data;
//	char par_name[12];
//	char obs_name[20];
//
//	// write header
//	tmp = -n_par;
//	jout.write((char*)&tmp, sizeof(tmp));
//	tmp = -n_obs_and_pi;
//	jout.write((char*)&tmp, sizeof(tmp));
//
//	//write number nonzero elements in jacobian (includes prior information)
//	n = matrix.nonZeros();
//	jout.write((char*)&n, sizeof(n));
//
//	//write matrix
//	n = 0;
//	map<string, double>::const_iterator found_pi_par;
//	map<string, double>::const_iterator not_found_pi_par;
//
//	Eigen::SparseMatrix<double> matrix_T(matrix);
//	matrix_T.transpose();
//	for (int icol = 0; icol<matrix.outerSize(); ++icol)
//	{
//		for (Eigen::SparseMatrix<double>::InnerIterator it(matrix_T, icol); it; ++it)
//		{
//			data = it.value();
//			n = it.row() + 1 + it.col() * matrix_T.rows();
//			jout.write((char*) &(n), sizeof(n));
//			jout.write((char*) &(data), sizeof(data));
//		}
//	}
//	//save parameter names
//	for (vector<string>::const_iterator b = col_names.begin(), e = col_names.end();
//		b != e; ++b) {
//		string l = pest_utils::lower_cp(*b);
//		pest_utils::string_to_fortran_char(l, par_name, 12);
//		jout.write(par_name, 12);
//	}
//
//	//save observation and Prior information names
//	for (vector<string>::const_iterator b = row_names.begin(), e = row_names.end();
//		b != e; ++b) {
//		string l = pest_utils::lower_cp(*b);
//		pest_utils::string_to_fortran_char(l, obs_name, 20);
//		jout.write(obs_name, 20);
//	}
//	//save observation names (part 2 prior information)
//	jout.close();
}

void Mat::to_binary_new(const string &filename)
{
	pest_utils::save_binary_extfmt(filename, row_names, col_names, matrix);
//	ofstream jout(filename, ios::out | ios::binary);
//	int n_par = col_names.size();
//	int n_obs_and_pi = row_names.size();
//	int n;
//	int tmp;
//	double data;
//	char par_name[200];
//	char obs_name[200];
//
//	// write header
//	tmp = n_par;
//	jout.write((char*)&tmp, sizeof(tmp));
//	tmp = n_obs_and_pi;
//	jout.write((char*)&tmp, sizeof(tmp));
//
//	//write number nonzero elements in jacobian (includes prior information)
//	n = matrix.nonZeros();
//	jout.write((char*)&n, sizeof(n));
//
//	//write matrix
//	n = 0;
//	map<string, double>::const_iterator found_pi_par;
//	map<string, double>::const_iterator not_found_pi_par;
//
//	Eigen::SparseMatrix<double> matrix_T(matrix);
//	matrix_T.transpose();
//	for (int icol = 0; icol<matrix.outerSize(); ++icol)
//	{
//		for (Eigen::SparseMatrix<double>::InnerIterator it(matrix_T, icol); it; ++it)
//		{
//			data = it.value();
//			n = it.row() - 1;
//			jout.write((char*) &(n), sizeof(n));
//			n = it.col() - 1;
//			jout.write((char*) &(n), sizeof(n));
//
//			jout.write((char*) &(data), sizeof(data));
//		}
//	}
//	//save parameter names
//	for (vector<string>::const_iterator b = col_names.begin(), e = col_names.end();
//		b != e; ++b) {
//		string l = pest_utils::lower_cp(*b);
//		pest_utils::string_to_fortran_char(l, par_name, 200);
//		jout.write(par_name, 200);
//	}
//
//	//save observation and Prior information names
//	for (vector<string>::const_iterator b = row_names.begin(), e = row_names.end();
//		b != e; ++b) {
//		string l = pest_utils::lower_cp(*b);
//		pest_utils::string_to_fortran_char(l, obs_name, 200);
//		jout.write(obs_name, 200);
//	}
//	//save observation names (part 2 prior information)
//	jout.close();
}

void Mat::from_binary(const string &filename)
{
	pest_utils::read_binary(filename, row_names, col_names, matrix);
}



//-----------------------------------------
//Maninpulate the shape and ordering of Mats
//-----------------------------------------

Mat Mat::leftCols(const int idx)
{
	vector<string> cnames;
	vector<string> base_cnames = *cn_ptr();
	for (int i = 0; i < idx; i++)
		cnames.push_back(base_cnames[i]);
	return Mat(row_names,cnames,matrix.leftCols(idx));
}

Mat Mat::rightCols(const int idx)
{
	vector<string> cnames;
	vector<string> base_cnames = *cn_ptr();
	for (int i = ncol() - idx; i < ncol(); i++)
		cnames.push_back(base_cnames[i]);
	return Mat(row_names,cnames,matrix.rightCols(idx));
}

Mat Mat::get(const vector<string> &new_row_names, const vector<string> &new_col_names, bool update)
{
	//check that every row and col name is listed
	if (new_row_names.size() == 0) throw runtime_error("Mat::get() error: new_row_names is empty");
	if (new_col_names.size() == 0) throw runtime_error("Mat::get() error: new_col_names is empty");
	vector<string> row_not_found;
	
	if (update)
		update_sets();
	set<string>::iterator end = row_set.end();
	for (auto &n : new_row_names)
		if (row_set.find(n) == end)
			row_not_found.push_back(n);

	vector<string> col_not_found;
	//set<string> col_set(col_names.begin(), col_names.end());
	end = col_set.end();
	for (auto &n : new_col_names)
		if (col_set.find(n) == end)
			col_not_found.push_back(n);
	if (row_not_found.size() != 0)
	{
		cout << "Mat::get() error: the following row names were not found:" << endl;
		for (auto &name : row_not_found)
			cout << name << ",";
		cout << endl;
	}

	if (col_not_found.size() != 0)
	{
		cout << "Mat::get() error: the following col names were not found:" << endl;
		for (auto &name : col_not_found)
			cout << name << ",";
		cout << endl;
	}

	if ((row_not_found.size() != 0) || (col_not_found.size() != 0))
	{
		throw runtime_error("Mat::get() error: at least one row or col name not found in Mat::get()");
	}


	int nrow = new_row_names.size();
	int ncol = new_col_names.size();
	int irow_new;
	int icol_new;

	unordered_map<string, int> row_name2new_index_map;
	unordered_map<string, int> col_name2new_index_map;

	// Build mapping of parameter names to column number in new matrix to be returned
	icol_new = 0;
	for (vector<string>::const_iterator b = new_col_names.begin(), e = new_col_names.end();
		b != e; ++b, ++icol_new) {
		col_name2new_index_map[(*b)] = icol_new;
	}

	// Build mapping of observation names to row  number in new matrix to be returned
	irow_new = 0;
	for (vector<string>::const_iterator b = new_row_names.begin(), e = new_row_names.end();
		b != e; ++b, ++irow_new) {
		row_name2new_index_map[(*b)] = irow_new;
	}

	unordered_map<string, int>::const_iterator found_col;
	unordered_map<string, int>::const_iterator found_row;
	unordered_map<string, int>::const_iterator not_found_col_map = col_name2new_index_map.end();
	unordered_map<string, int>::const_iterator not_found_row_map = row_name2new_index_map.end();

	const string *row_name;
	const string *col_name;
	std::vector<Eigen::Triplet<double> > triplet_list;
	for (int icol = 0; icol<matrix.outerSize(); ++icol)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(matrix, icol); it; ++it)
		{
			col_name = &col_names[it.col()];
			row_name = &row_names[it.row()];
			found_col = col_name2new_index_map.find(*col_name);
			found_row = row_name2new_index_map.find(*row_name);
			if (found_col != not_found_col_map && found_row != not_found_row_map)
			{
				triplet_list.push_back(Eigen::Triplet<double>(found_row->second, found_col->second, it.value()));
			}
		}
	}
	//if (triplet_list.size() == 0)
		//throw runtime_error("Mat::get() error: triplet list is empty");

	Eigen::SparseMatrix<double> new_matrix(nrow, ncol);
	new_matrix.setZero();
	if (triplet_list.size() > 0)
		new_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
	return Mat(new_row_names,new_col_names,new_matrix,this->mattype);
}

Mat Mat::extract(const vector<string> &extract_row_names, const vector<string> &extract_col_names)
{
	Mat new_mat;
	if ((extract_row_names.size() == 0) && (extract_col_names.size() == 0))
		throw runtime_error("Mat::extract() error: extract_rows and extract_cols both empty");
	else if (extract_row_names.size() == 0)
	{
		new_mat = get(row_names, extract_col_names);
		drop_cols(extract_col_names);
	}
	else if (extract_col_names.size() == 0)
	{
		new_mat = get(extract_row_names, col_names);
		drop_rows(extract_row_names);
	}
	else
	{
		new_mat = get(extract_row_names, extract_col_names);
		drop_rows(extract_row_names);
		drop_cols(extract_col_names);
	}
	return new_mat;
}

Mat Mat::extract(const string &extract_row_name, const vector<string> &extract_col_names)
{
	vector<string> extract_row_names;
	extract_row_names.push_back(extract_row_name);
	return extract(extract_row_names, extract_col_names);
}
Mat Mat::extract(const vector<string> &extract_row_names, const string &extract_col_name)
{
	vector<string> extract_col_names;
	extract_col_names.push_back(extract_col_name);
	return extract(extract_row_names, extract_col_names);
}

bool Mat::isdiagonal()
{
	if (mattype == MatType::DIAGONAL)
		return true;
	return false;
}

void Mat::drop_cols(const vector<string> &drop_col_names)
{
	vector<string> missing_col_names;
	const set<string> snames(col_names.begin(), col_names.end());
	set<string>::const_iterator send = snames.end();
	for (auto &name : drop_col_names)
	{
		//if (find(col_names.begin(), col_names.end(), name) == col_names.end())
		if (snames.find(name) == send)
			missing_col_names.push_back(name);
	}

	if (missing_col_names.size() != 0)
	{
		cout << "Mat::drop_cols() error: the following drop_col_names were not found:" << endl;
		for (auto &name : drop_col_names)
			cout << name << ',';
		cout << endl;
		throw runtime_error("Mat::drop_cols() error: at least one drop col name not found");
	}
	vector<string> new_col_names;
	if (drop_col_names.size() == 0)
		new_col_names = col_names;
	else
	{

		const set<string> snames(drop_col_names.begin(), drop_col_names.end());
		send = snames.end();
		for (auto &name : col_names)
		{
			//if (find(drop_col_names.begin(), drop_col_names.end(), name) == drop_col_names.end())
			if (snames.find(name) == send)
				new_col_names.push_back(name);
		}
	}
	Mat new_mat = get(row_names, new_col_names);
	matrix = new_mat.get_matrix();
	col_names = new_col_names;
	mattype = new_mat.get_mattype();
}

void Mat::drop_rows(const vector<string> &drop_row_names)
{

	vector<string> missing_row_names;
	const set<string> snames(row_names.begin(), row_names.end());
	set<string>::const_iterator send = snames.end();

	for (auto &name : drop_row_names)
	{
		//if (find(row_names.begin(), row_names.end(), name) == row_names.end())
		if (snames.find(name) == send)
			missing_row_names.push_back(name);
	}

	if (missing_row_names.size() != 0)
	{
		cout << "Mat::drop_rows() error: the following drop_row_names were not found:" << endl;
		for (auto &name : drop_row_names)
			cout << name << ',';
		cout << endl;
		throw runtime_error("Mat::drop_rows() error: at least one drop row name not found");
	}

	vector<string> new_row_names;
	if (drop_row_names.size() == 0)
		new_row_names = row_names;
	else
	{
		const set<string> snames(drop_row_names.begin(), drop_row_names.end());
		send = snames.end();
		for (auto &name : row_names)
		{
			//if (find(drop_row_names.begin(), drop_row_names.end(), name) == drop_row_names.end())
			if (snames.find(name) == send)
				new_row_names.push_back(name);
		}
	}
	if (new_row_names.size() == 0)
		matrix = Eigen::SparseMatrix<double>();
	else
	{
		Mat new_mat = get(new_row_names, col_names);
		matrix = new_mat.get_matrix();
		mattype = new_mat.get_mattype();
	}
	row_names = new_row_names;

}



//-----------------------------------------
//covariance matrices
//-----------------------------------------
Covariance::Covariance(string filename)
{
	mattype = MatType::SPARSE;
	string ext = filename.substr(filename.find_last_of(".") + 1);
	vector<string> empty;
	pest_utils::upper_ip(ext);
	if (ext == "PST")
		throw runtime_error("Cov::Cov() error: cannot instantiate a cov with PST");
	else if (ext == "MAT")
		from_ascii(filename);
	else if (ext == "UNC")
		from_uncertainty_file(filename, empty);
	else
		throw runtime_error("Cov::Cov() error: only .unc or .mat\
														 files can be used to instantiate a Cov");
}

Covariance::Covariance(vector<string> &names)
{
	row_names = names;
	col_names = names;
	icode = 1;
	mattype = MatType::SPARSE;
}

Covariance::Covariance()
{
	icode = 1;
	mattype = MatType::SPARSE;
}

Covariance::Covariance(vector<string> _names, Eigen::SparseMatrix<double> _matrix, Mat::MatType _mattype)
{
	if ((_names.size() != _matrix.rows()) || (_names.size() != _matrix.cols()))
		throw runtime_error("Covariance::Covariance() error: names.size() does not match matrix dimensions");
	Eigen::SparseMatrix<double> test = _matrix;
	matrix = _matrix;
	row_names = _names;
	col_names = _names;
	icode = 1;
	mattype =_mattype;
}

Covariance::Covariance(Mat _mat)
{
	if (_mat.get_row_names() != _mat.get_col_names())
		throw runtime_error("Cov::Cov() error instantiating Covariance from Mat: row_names != col_names");
	row_names = _mat.get_row_names();
	col_names = _mat.get_col_names();
	matrix = _mat.get_matrix();
	icode = 1;
	mattype = _mat.get_mattype();
}


void Covariance::from_diagonal(Covariance &other)
{
	row_names = other.get_row_names();
	col_names = other.get_col_names();
	if (other.get_mattype() == Mat::MatType::DIAGONAL)
	{
		matrix = other.get_matrix();
	}
	else
	{
		Eigen::MatrixXd temp = other.e_ptr()->diagonal().asDiagonal();
		Eigen::SparseMatrix<double> temp2 = temp.sparseView();
		matrix = temp2;
	}

}
Covariance Covariance::diagonal(double val)
{
	vector<Eigen::Triplet<double>> triplet_list;
	for (int i = 0; i != nrow(); i++)
		triplet_list.push_back(Eigen::Triplet<double>(i, i, val));
	Eigen::SparseMatrix<double> i(nrow(), ncol());
	i.setZero();
	i.setFromTriplets(triplet_list.begin(), triplet_list.end());
	return Covariance(*rn_ptr(), i);
}

string Covariance::try_from(Pest &pest_scenario, FileManager &file_manager, bool is_parcov, bool forgive_missing)
{
	stringstream how;
	stringstream ss;
	string cov_fname;
	vector<string> ordered_names;
	if (is_parcov)
	{
		cov_fname = pest_scenario.get_pestpp_options().get_parcov_filename();
		ordered_names = pest_scenario.get_ctl_ordered_adj_par_names();
	}
	else
	{
		cov_fname = pest_scenario.get_pestpp_options().get_obscov_filename();
		ordered_names = pest_scenario.get_ctl_ordered_nz_obs_names();
	}
	if (!cov_fname.empty())
	{
		string ext = cov_fname.substr(cov_fname.size() - 3, 3);
		pest_utils::upper_ip(ext);
		if (ext == "UNC")
		{
			try
			{
				from_uncertainty_file(cov_fname, ordered_names);
				how << "from unc file " << cov_fname;
			}
			catch (exception &e)
			{
				ss << "Cov::try_from() error reading uncertainty file " << cov_fname << " :" << e.what();
				throw runtime_error(ss.str());
			}
		}
		else
		{
			try
			{
				Mat::from_file(cov_fname);
				how << " from file " << cov_fname;
			}
			catch (exception &e)
			{
				ss << "Cov:try_from() error reading from file " << cov_fname << " :" << e.what();
				throw runtime_error(ss.str());
			}
		}
	}
	else
	{
		if (is_parcov)
		{
			from_parameter_bounds(pest_scenario, file_manager.rec_ofstream());
			how << "from parameter bounds, using par_sigma_range " << pest_scenario.get_pestpp_options().get_par_sigma_range();
		}
		else
		{
			from_observation_weights(pest_scenario, file_manager.rec_ofstream());
			how << "from observation weights";
		}
		}
	//check that the parcov matrix has the right parameter names
	vector<string> missing;
	if ((is_parcov) && (cov_fname.size() > 0))
	{
		set<string> parcov_names(row_names.begin(), row_names.end());
		const ParameterRec *prec;
		for (auto &pname : pest_scenario.get_ctl_ordered_par_names())
		{
			prec = pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(pname);
			if ((prec->tranform_type == ParameterRec::TRAN_TYPE::LOG) ||
				(prec->tranform_type == ParameterRec::TRAN_TYPE::NONE))
			{
				if (parcov_names.find(pname) == parcov_names.end())
				{
					missing.push_back(pname);
				}
			}
		}
		if (missing.size() > 0)
		{
			ofstream& frec = file_manager.rec_ofstream();
			frec << "...Note: parcov missing the following " << missing.size() << " adjustable parameters:" << endl;
			int i = 0;
			for (auto& pname : missing)
			{
				frec << ',' << pname;
				i++;
				if (i > 10)
				{
					frec << endl;
					i = 0;
				}
			}
			frec << endl;
			if (forgive_missing)
			{
				cout << "WARNING: " << missing.size() << " adjustable parameters missing from parcov, continuing..." << endl;
			}
			else
			{
				ss.str("");
				ss << "parcov missing " << missing.size() << "adjustable parameters, see rec file for listing";
				throw PestError(ss.str());
			}
		}
		
	}
	if ((!is_parcov) && (cov_fname.size() > 0))
	{
		double weight;
		set<string> cov_names(row_names.begin(), row_names.end());
		for (auto &oname : pest_scenario.get_ctl_ordered_obs_names())
		{
			weight = pest_scenario.get_ctl_observation_info().get_weight(oname);
			if (weight == 0.0)
				continue;
			if (cov_names.find(oname) == cov_names.end())
			{
				missing.push_back(oname);
			}

		}
		if (missing.size() > 0)
		{
			ofstream& frec = file_manager.rec_ofstream();
			frec << "...Note: obscov missing the following " << missing.size() << " non-zero weighted obs:" << endl;
			int i = 0;
			for (auto& name : missing)
			{
				frec << ',' << name;
				i++;
				if (i > 10)
				{
					frec << endl;
					i = 0;
				}
			}
			frec << endl;
			if (forgive_missing)
			{
				cout << "WARNING: " << missing.size() << " non-zero weighted observations missing from obscov, continuing..." << endl;
			}
			else
			{
				ss.str("");
				ss << "obscov missing " << missing.size() << "non-zero weighted observations, see rec file for listing";
				throw PestError(ss.str());
			}
		}
	}

	vector<string> extra;
	set<string> allowed_names(ordered_names.begin(), ordered_names.end());
	ordered_names.clear();

	for (auto name : row_names)
	{
		if (allowed_names.find(name) == allowed_names.end())
			extra.push_back(name);
	}
	if (extra.size() > 0)
	{
		ss.str("");
		ss << "WARNING: Cov::try_from(): " << extra.size() << " extra elements in covariance matrix being drop - these are probably for fixed parameters and/or zero-weight observations";
		//for (auto name : extra)
		//	ss << " " << name;
		//throw PestError(ss.str());
		file_manager.rec_ofstream() << ss.str() << endl << endl;
		//cout << "WARNING: " << extra.size() << " unrecognized elements in covariance matrix being dropped, see .rec file for listing" << endl;
		drop(extra);
	}
	return how.str();
}


Covariance Covariance::get(const vector<string> &other_names, bool update)
{
	Covariance new_cov(Mat::get(other_names, other_names, update));
	return new_cov;
}

Covariance Covariance::extract(vector<string> &extract_names)
{
	Covariance new_cov(Mat::extract(extract_names, extract_names));
	return new_cov;
}

void Covariance::drop(vector<string> &drop_names)
{
	drop_rows(drop_names);
	drop_cols(drop_names);
}

void Covariance::from_uncertainty_file(const string &filename, vector<string> &ordered_names)
{
	ifstream in(filename);
	if (!in.good())
	{
		throw runtime_error("Cov::from_uncertainty_file() error: cannot open " + filename + " to read uncertainty file: "+filename);
	}
	mattype = MatType::DIAGONAL;
	vector<Eigen::Triplet<double>> triplet_list;
	vector<string> names;
	string line,name,word;
	double val;
	vector<string> tokens;
	int irow=0, jcol=0;

	while (getline(in, line))
	{
		if (line.substr(0, 1).find("#") != string::npos)
		{
			continue;
		}
		pest_utils::upper_ip(line);
		//if this is the start of some block
		if (line.find("START") != string::npos)
		{
			if (line.find("STANDARD_DEVIATION") != string::npos)
			{
				while (true)
				{
					if (!getline(in, line))
						throw runtime_error("Cov::from_uncertainty_file() error:EOF encountered while reading standard_deviation block\
							from uncertainty file:" + filename);
					pest_utils::upper_ip(line);
					if (line.find("END STANDARD_DEVIATION") != string::npos)
					{
						break;
					}
					
					tokens.clear();
					pest_utils::tokenize(line, tokens);
					pest_utils::convert_ip(tokens[1], val);
					name = tokens[0];
					double std_mlt = 1.0;
					if (name == "STD_MULTIPLIER")
					{
						std_mlt = val;
						continue;
					}
					if (find(names.begin(), names.end(), name) != names.end())
						throw runtime_error(name + " listed more than once in uncertainty file:" + filename);
					names.push_back(tokens[0]);
					triplet_list.push_back(Eigen::Triplet<double>(irow, jcol, (val * std_mlt)*(val * std_mlt)));
					irow++, jcol++;
				}

			}
			else if (line.find("COVARIANCE_MATRIX") != string::npos)
			{
				string cov_filename = "none";
				double var_mult = 1.0;
				string start_par = "", end_par = "";
				string par_list_file = "";
				while (true)
				{
					if (!getline(in, line))
						throw runtime_error("Cov::from_uncertainty_file() error:EOF encountered while reading covariance_matrix block\
							from uncertainty file:" + filename);
					// keep line in original case to preserve filename
					line.erase(remove(line.begin(), line.end(), '\"'), line.end());
					line.erase(remove(line.begin(), line.end(), '\''), line.end());
					string upper_line = line;
					pest_utils::upper_ip(upper_line);
					if (upper_line.find("END") != string::npos) break;

					tokens.clear();
					pest_utils::tokenize(line, tokens);
					word = pest_utils::upper_cp(tokens[0]);
					if (word.find("FILE") != string::npos)
						cov_filename = tokens[1];
					else if (word.find("VARIANCE_MULTIPLIER") != string::npos)
						pest_utils::convert_ip(tokens[1], var_mult);
					else if (word.find("FIRST_PARAMETER") != string::npos)
						start_par = pest_utils::upper_cp(tokens[1]);
					else if (word.find("LAST_PARAMETER") != string::npos)
						end_par = pest_utils::upper_cp(tokens[1]);
					else if (word.find("PARAMETER_LIST_FILE"))
						par_list_file = tokens[1];
					else
						throw runtime_error("Cov::from_uncertainty_file() error:unrecognized token:" + tokens[0] + " in covariance matrix block in uncertainty file:" + filename);
				}

				if ((start_par.size() > 0) && (end_par.size() == 0))
					throw runtime_error("Cov::from_uncertainty_file() error: 'FIRST PARAMETER' passed but 'LAST_PARAMETER' was not");
				if ((start_par.size() == 0) && (end_par.size() > 0))
					throw runtime_error("Cov::from_uncertainty_file() error: 'LAST PARAMETER' passed but 'FIRST_PARAMETER' was not");

				//read the covariance matrix
				Covariance cov;
				cov.from_ascii(cov_filename);

				if ((start_par.size() > 0) && (end_par.size() > 0))
				{
					if (par_list_file.size() > 0)
						throw runtime_error("Cov::from_uncertainty_file() error: both 'PARAMETER_LIST_FILE' AND 'FIRST_PARAMETER'/'LAST_PARAMETER' supplied");
					if (ordered_names.size() == 0)
						throw runtime_error("Cov::from_uncertainty_file() error: ordered_names arg req for first_par/last_par option");
					vector<string>::iterator first, last;
					first = find(ordered_names.begin(), ordered_names.end(), start_par);
					if (first == ordered_names.end())
						throw runtime_error("Cov::from_uncertainty_file() error: couldn't find 'FIRST_PARAMETER' " + start_par);
					last = find(ordered_names.begin(), ordered_names.end(), end_par);
					if (last == ordered_names.end())
						throw runtime_error("Cov::from_uncertainty_file() error: couldn't find 'LAST_PARAMETER' " + end_par);
					vector<string> ordered_names_matrix(first, last+1);
					if (ordered_names_matrix.size() != cov.col_names.size())
					{
						stringstream ss;
						ss << "Cov::from_uncertainty_file() error: number of elements in covariance matrix " << cov_filename << " (" << cov.col_names.size();
						ss << ") different from number of elements between 'FIRST_PARAMETER' and 'LAST_PARAMETER' (" << ordered_names_matrix.size() << ")";
						throw runtime_error(ss.str());
					}
					cov.row_names = ordered_names_matrix;
					cov.col_names = ordered_names_matrix;
				}
				else if (par_list_file.size() > 0)
				{
					vector<string> par_names;
					try
					{
						par_names = pest_utils::read_onecol_ascii_to_vector(par_list_file);
					}
					catch (...)
					{
						throw runtime_error("Cov::from_uncertainty_file() error reading 'PARAMETER_LIST_FILE' " + par_list_file);
					}
					if (par_names.size() != cov.col_names.size())
					{
						stringstream ss;
						ss << "Cov::from_uncertainty_file() error: number of elements in covariance matrix " << cov_filename;
						ss << " (" << cov.col_names.size() << ") different from number of elements in 'PARAMETER_LIST_FILE' (" << par_names.size() << ")";
						throw runtime_error(ss.str());
					}
					cov.row_names = par_names;
					cov.col_names = par_names;
				}



				//check that the names in the covariance matrix are not already listed
				vector<string> dup_names;
				for (auto &_name : cov.get_row_names())
				{
					if (find(names.begin(), names.end(), _name) != names.end())
						dup_names.push_back(_name);
					else
						names.push_back(_name);
				}
				if (dup_names.size() != 0)
				{
					cout << "the following names from covariance matrix file " << cov_filename << " have already be found in uncertainty file " << filename << endl;
					for (auto &_name : dup_names)
						cout << _name << ',';
					cout << endl;
					throw runtime_error("Cov::from_uncertainty_file() error:at least one name in covariance matrix " + cov_filename + " is already listed in uncertainty file: " + filename);
				}

				//build triplets from the covariance matrix
				int start_irow = irow;
				Eigen::SparseMatrix<double> cov_matrix = cov.get_matrix();
				for (int icol = 0; icol < cov_matrix.outerSize(); ++icol)
				{
					for (Eigen::SparseMatrix<double>::InnerIterator it(cov_matrix, icol); it; ++it)
					{
						triplet_list.push_back(Eigen::Triplet<double>(start_irow + it.row(), jcol, var_mult * it.value()));
						irow++;
					}
					jcol++;
					irow = start_irow;
				}
				mattype = MatType::SPARSE;
				irow = jcol;
			}
			else
				throw runtime_error("Cov::from_uncertainty_file() error:unrecognized block:" + line + " in uncertainty file:" + filename);
		}
	}

	Eigen::SparseMatrix<double> new_matrix(names.size(), names.size());
	new_matrix.setZero();  // initialize all entries to 0
	new_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
	//cout << new_matrix.diagonal() << endl;
	matrix = new_matrix;
	row_names = names;
	col_names = names;
}

void Covariance::from_parameter_bounds(ofstream& frec, const vector<string> &par_names,const ParameterInfo &par_info, 
	map<string, double>& par_std, double sigma_range)
{
	matrix.resize(0, 0);
	row_names.clear();
	col_names.clear();
	vector<Eigen::Triplet<double>> triplet_list;
	const ParameterRec* par_rec;
	int i = 0;
	double upper, lower;
	for (auto par_name : par_names)
	{
		pest_utils::upper_ip(par_name);
		par_rec = par_info.get_parameter_rec_ptr(par_name);
		if ((par_rec->tranform_type == ParameterRec::TRAN_TYPE::FIXED) || (par_rec->tranform_type == ParameterRec::TRAN_TYPE::TIED))
		{
			continue;
		}
		if (par_std.find(par_name) != par_std.end())
		{
			triplet_list.push_back(Eigen::Triplet<double>(i, i, pow(par_std[par_name], 2.0)));
			row_names.push_back(par_name);
			col_names.push_back(par_name);
			i++;
		}
		else
		{
			upper = par_rec->ubnd;
			lower = par_rec->lbnd;
			if (par_rec->tranform_type == ParameterRec::TRAN_TYPE::LOG)
			{
				upper = log10(upper);
				lower = log10(lower);
			}
			row_names.push_back(par_name);
			col_names.push_back(par_name);
			//double temp = pow((upper - lower) / 4.0,2.0);
			triplet_list.push_back(Eigen::Triplet<double>(i, i, pow((upper - lower) / sigma_range, 2.0)));
			i++;
		}
	}
	if (triplet_list.size() > 0)
	{
		matrix.resize(row_names.size(), row_names.size());
		matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
	}
	else
	{
		throw runtime_error("Cov::from_parameter_bounds() error:Error loading covariance from parameter bounds: no non-fixed/non-tied parameters found");
	}
	mattype = Mat::MatType::DIAGONAL;
}



void Covariance::from_parameter_bounds(Pest &pest_scenario, ofstream& frec)
{
	map<string, double> par_std = pest_scenario.get_ext_file_double_map("parameter data external", "standard_deviation");
	const ParameterInfo pi = pest_scenario.get_ctl_parameter_info();
    if (par_std.size() > 0)
	{
		frec << "Note: the following parameters have 'standard_deviation' defined - this will be used" << endl;
		frec << "      instead of bounds for the prior parameter covariance matrix : " << endl;
		vector<string> remove;
		for (auto pname : pest_scenario.get_ctl_ordered_par_names())
		{
			if (par_std.find(pname) != par_std.end())
			{
				if (par_std[pname] <= 0.0)
				{
                    ParameterRec::TRAN_TYPE t = pi.get_parameter_rec_ptr(pname)->tranform_type;
                    if ((t != ParameterRec::TRAN_TYPE::FIXED) && (t != ParameterRec::TRAN_TYPE::TIED)) {
                        frec << "Warning: parameter " << pname
                             << " 'standard_deviation' less than or equal to zero, using bounds instead" << endl;
                        remove.push_back(pname);
                    }
				}
				else
					frec << pname << ' ' << par_std[pname] << endl;
			}
				
		}
		for (auto r : remove)
			par_std.erase(r);
	}

	from_parameter_bounds(frec, pest_scenario.get_ctl_ordered_par_names(), pest_scenario.get_ctl_parameter_info(),
		par_std,pest_scenario.get_pestpp_options().get_par_sigma_range());
}

void Covariance::from_parameter_bounds(const string &pst_filename, ofstream& frec)
{
	ifstream ipst(pst_filename);
	if (!ipst.good()) throw runtime_error("Cov::from_parameter_bounds() error opening pst file: " + pst_filename);
	Pest pest_scenario;
	pest_scenario.process_ctl_file(ipst, pst_filename);
	from_parameter_bounds(pest_scenario, frec);
}

void Covariance::from_observation_weights(const string &pst_filename, ofstream& frec)
{
	ifstream ipst(pst_filename);
	if (!ipst.good()) throw runtime_error("Cov::from_observation_weights() error opening pst file: " + pst_filename);
	Pest pest_scenario;
	pest_scenario.process_ctl_file(ipst, pst_filename);
	from_observation_weights(pest_scenario,frec);

}


void Covariance::from_observation_weights(ofstream& frec, const vector<string>& obs_names, const ObservationInfo& obs_info, 
	const vector<string>& pi_names, const PriorInformation* pi, map<string,double>& obs_std)
{
	matrix.resize(0, 0);
	row_names.clear();
	col_names.clear();
	vector<Eigen::Triplet<double>> triplet_list;
	const ObservationRec* obs_rec;
	int i = 0;
	double weight = 0;
	for (auto obs_name : obs_names)
	{
		pest_utils::upper_ip(obs_name);
		obs_rec = obs_info.get_observation_rec_ptr(obs_name);
		if (obs_std.find(obs_name) != obs_std.end())
		{
			triplet_list.push_back(Eigen::Triplet<double>(i, i, pow(obs_std[obs_name], 2)));
			row_names.push_back(obs_name);
			col_names.push_back(obs_name);
			i++;
		}
		else
		{
			obs_rec = obs_info.get_observation_rec_ptr(obs_name);
			weight = obs_rec->weight;
			if (weight <= 0.0)
				weight = 1.0e+60;
			else
				weight = pow(1.0 / obs_rec->weight, 2.0);
			triplet_list.push_back(Eigen::Triplet<double>(i, i, weight));
			row_names.push_back(obs_name);
			col_names.push_back(obs_name);
			i++;
		}
	}

	/*PriorInformation::const_iterator pi_iter;
	PriorInformation::const_iterator not_pi_iter = pi->end();

	for (auto pi_name : pi_names)
	{
		pi_iter = pi->find(pi_name);
		if (pi_iter != not_pi_iter)
		{
			weight = pi_iter->second.get_weight();
			if (weight <= 0.0) weight = 1.0e-30;
			triplet_list.push_back(Eigen::Triplet<double>(i, i, pow(1.0 / weight, 2.0)));
			row_names.push_back(pi_name);
			col_names.push_back(pi_name);
			i++;
		}
	}*/
	if (row_names.size() > 0)
	{
		matrix.resize(row_names.size(), row_names.size());
		matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
	}
	else
	{
		throw runtime_error("Cov::from_observation_weights() error:Error loading covariance from obs weights: no non-zero weighted obs found");
	}
	mattype = Mat::MatType::DIAGONAL;
	/*if (mattype == Mat::MatType::DIAGONAL)
		cout << "diagonal" << endl;*/
}


void Covariance::from_observation_weights(Pest &pest_scenario, ofstream& frec)
{
	map<string, double> obs_std = pest_scenario.get_ext_file_double_map("observation data external", "standard_deviation");
	vector<string> remove;
    const ObservationInfo oi = pest_scenario.get_ctl_observation_info();
	if (obs_std.size() > 0)
	{
		frec << "Note: the following observations have 'standard_deviation' defined - this will be used" << endl;
		frec << "      instead of weight for the observation noise covariance matrix : " << endl;
		for (auto oname : pest_scenario.get_ctl_ordered_obs_names())
		{
			if (obs_std.find(oname) != obs_std.end())
			{
				if ((obs_std[oname] <= 0.0) && (oi.get_weight(oname) > 0.0))				{
					frec << "Warning: observation " << oname << " 'standard_deviation' less than or equal to zero, using weight" << endl;
					remove.push_back(oname);
				}
				else
					frec << oname << ' ' << obs_std[oname] << endl;
			}
		}
		for (auto r : remove)
			obs_std.erase(r);
	}

	from_observation_weights(frec, pest_scenario.get_ctl_ordered_obs_names(), pest_scenario.get_ctl_observation_info(),
		pest_scenario.get_ctl_ordered_pi_names(), pest_scenario.get_prior_info_ptr(), obs_std);

}

void Covariance::to_uncertainty_file(const string &filename)
{
	ofstream out(filename);
	if (!out.good())
	{
		throw runtime_error("Cov::to_uncertainty_file() error opening file: " + filename + " to write an uncertainty file");
	}

	//check if diagonal, write stdevs
	if (mattype == Mat::MatType::DIAGONAL)
	{
		Eigen::VectorXd vec(matrix.diagonal());
		out << "START STANDARD_DEVIATION" << endl;
		int i=0;
		for (vector<string>::iterator name = row_names.begin(); name != row_names.end(); ++name, i++)
		{
			out << "  " << setw(20) << left << *name << "  " << setw(20) << left << vec(i) << endl;
		}
		out << "END STANDARD_DEVIATION" << endl;
		out.close();
	}
	else
	{
		out << "START COVARIANCE_MATRIX" << endl;
		out << "  file emu_cov.mat" << endl;
		out << " variance multiplier 1.0" << endl;
		out << "END COVARIANCE_MATRIX" << endl;
		out.close();
		to_ascii("emu_cov.mat");
	}
}

void Covariance::cholesky()
{
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt;
	lower_cholesky = llt.matrixL();
}

vector<Eigen::VectorXd> Covariance::draw(int ndraws)
{
	throw runtime_error("Covariance::draw() not implemented");
}

vector<double> Covariance::standard_normal(default_random_engine gen)
{
	normal_distribution<double> stanard_normal(0.0, 1.0);
	vector<double> sn_vec;
	for (auto &name : row_names)
	{
		sn_vec.push_back(stanard_normal(gen));
	}
	return sn_vec;
}

#include <string>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <list>
#include <cstring>
#include <cmath>
#include <cassert>
#include <mutex>
#include "config_os.h"
#include "Transformable.h"
#include "network_package.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>


//template <class T> const T& min ( const T& a, const T& b );



#include "utilities.h"
#include "system_variables.h"


using namespace std;

std::ostream& operator<< (std::ostream &os, const std::set<std::string> val)
{
	for (const auto &i : val)
	{
		os << i << endl;
	}
	return os;
}

std::ostream& operator<< (std::ostream &os, const std::vector<std::string> val)
{
	for (const auto &i : val)
	{
		os << i << endl;
	}
	return os;
}

void print(std::set<std::string> val, std::ostream &os, int indent)
{
	string space(indent, ' ');
	for (const auto &i : val)
	{
		os << space << i << endl;
	}
}




namespace pest_utils
{


double get_duration_sec(std::chrono::system_clock::time_point start_time)
{
	chrono::system_clock::duration dt = chrono::system_clock::now() - start_time;
	return (double)std::chrono::duration_cast<std::chrono::milliseconds>(dt).count() / 1000.0;
}

template < class ContainerT >
void tokenize(const std::string& str, ContainerT& tokens, const std::string& delimiters, const bool trimEmpty)
{
/*
 I usually choose to use std::vector<std::string> types as my second parameter
  (ContainerT)... but list<> is way faster than vector<> for when direct access
   is not needed, and you can even create your own string class and use something
    like std::list<SubString> where SubString does not do any copies for incredible speed increases.
It's more than double as fast as the fastest tokenize on this page and almost 5 times
faster than some others. Also with the perfect parameter types you can eliminate all string and list copies.
Additionally it does not do the (extremely inefficient) return of result, but rather
it passes the tokens as a reference, thus also allowing you to build up tokens using
multiple calls if you so wished.Lastly it allows you to specify whether to trim empty
tokens from the results via a last optional parameter.All it needs is std::string...
the rest are optional. It does not use streams or the boost library, but is flexible
enough to be able to accept some of these foreign types naturally.
 */
	std::string::size_type pos, lastPos = 0;
	while(true)
	{
		pos = str.find_first_of(delimiters, lastPos);
		if(pos == std::string::npos)
		{
			pos = str.length();
			if(pos != lastPos || !trimEmpty)
			{
				tokens.push_back(typename ContainerT::value_type(str.data()+lastPos, typename ContainerT::value_type::size_type(pos-lastPos)));
			}
			break;
		}
		else
		{
			if(pos != lastPos || !trimEmpty)
				tokens.push_back(typename ContainerT::value_type(str.data()+lastPos, typename ContainerT::value_type::size_type(pos-lastPos )));
		}

		lastPos = pos + 1;
	}
}

// Instantiate tokenize for the explicitly specified vector<string> container
template void tokenize(const std::string& str, vector<string>& tokens, const std::string& delimiters, const bool trimEmpty);
template void tokenize(const std::string& str, list<string>& tokens, const std::string& delimiters, const bool trimEmpty);

std::string& strip_ip(string &s, const string &op, const string &delimiters)
{
	size_t i;

	if (op == "both" || op == "front")
	{
		i = s.find_first_not_of(delimiters);
		if (i>0) s.erase(0, i);
	}

	if (op == "both" || op == "back")
	{
		i = s.find_last_not_of(delimiters);
		if (i == string::npos) i=0;
		else i+=1;
		if (i<s.size()) s.erase(i);
	}
	return s;
}

string strip_cp(const string &s, const string &op, const string &delimiters)
{
	string s2(s);
	strip_ip(s2, op, delimiters);
	return s2;
}

void upper_ip(string &s)
{
	for(unsigned int i=0; i<s.length(); i++)
	{
      s[i] = toupper(s[i]);
   }
}

string upper_cp(const string &s)
{
	string s2(s);
	upper_ip(s2);
	return s2;
}

string upper(char *txt)
{
	string tmp = txt;
	upper_ip(tmp);
	return tmp;
}
void lower_ip(string &s)
{
	string new_s;
	for(unsigned int i=0; i<s.length(); i++)
	{
		s[i] = tolower(s[i]);
	}
}

string lower_cp(const string &s)
{
	string s2(s);
	lower_ip(s2);
	return s2;
}

string get_base_filename(const string &s)
{
	string base_name;
	size_t end;
	end = s.find_last_of(".");
	if (end == string::npos) {end = s.size();}
	base_name = s.substr(0, end);
	return base_name;
}

void string_to_fortran_char(string in, char out[], int length, CASE_CONV conv_type)
{
	int str_len = in.size();
	if (conv_type == TO_LOWER)
	{
		lower_ip(in);
	}
	else if (conv_type == TO_UPPER)
	{
		upper_ip(in);
	}
	str_len = min(str_len, length);
	memset(out, ' ', length);
	memcpy(out, in.data(), str_len);
}

vector<char> string_as_fortran_char_ptr(string in, int _size)
{
	vector<char> out;
	out.assign(_size, ' ');
	string_to_fortran_char(in, out.data(), _size);
	return out;
}

StringvecFortranCharArray::StringvecFortranCharArray(const vector<string> in, int length, CASE_CONV conv_type)
{
	fort_array = new char[in.size() * length];
	for(int i=0, n_str=in.size(); i<n_str; ++i)
	{
		string_to_fortran_char(in[i], &fort_array[i*length], length, conv_type);
	}
}

char *StringvecFortranCharArray::get_prt()
{
	return fort_array;
}

StringvecFortranCharArray::~StringvecFortranCharArray()
{
	delete [] fort_array;
}


string get_filename_without_ext(const string &filename)
{
	// remove .pst or .PST from the end of the filename
	string new_str = filename;
	string filename_lower = lower_cp(filename);
	//size_t found = filename_lower.find_last_of('.');
	size_t found = filename_lower.find(".pst");
	if (found != string::npos)
	{
		new_str = new_str.substr(0, found);
	}
	return new_str;
}

string get_filename_ext(const string &filename)
{
	// remove .pst or .PST from the end of the filename
	string new_str = "";
	size_t found = filename.find_last_of(".");
	if (found != string::npos)
	{
		new_str = filename.substr(found + 1);
	}
	return new_str;
}

string get_filename(const string &complete_path)
{
	vector<string> tokens;
	tokenize(complete_path, tokens, OperSys::DIR_SEP);
	string filename = tokens.back();
	strip_ip(filename);
	return filename;
}
string get_pathname(const string &complete_path)
{
	vector<string> tokens;
	stringstream ret_val;
	int ntokens;

	tokenize(complete_path, tokens, OperSys::DIR_SEP);
	ntokens = tokens.size();

	if (complete_path.find(OperSys::DIR_SEP)==0){
		ret_val << OperSys::DIR_SEP;
	}
	if (ntokens >1) {
		int len = ntokens - 1;
		for (int i=0; i<len; ++i) {
			ret_val << tokens[i] << OperSys::DIR_SEP;
		}
	}
	return ret_val.str();
}

String2CharPtr::String2CharPtr(const std::string &str)
{
	my_str.resize(str.size()+1);
	std::copy(str.begin(), str.end(), my_str.begin());
	my_str.back() = '\0';
}

char* String2CharPtr::get_char_ptr()
{
	return &(my_str[0]);
}

void copyfile(const string &from_file, const string &to_file)
{
	fstream source(from_file, ios::binary);
    ofstream dest(to_file, ios::binary);

    dest << source.rdbuf();

    source.close();
    dest.close();
}

template <class keyType, class dataType>
vector<keyType> get_map_keys(const map<keyType,dataType> &my_map)
{
	vector<keyType> keys;
	for(auto &i : my_map)
	{
		keys.push_back(i.first);
	}
	return keys;
}

template vector<string> get_map_keys(const map<string, map<string, double>> &my_map);



string fortran_str_2_string(char *fstr, int str_len)
{
    string new_str = string(fstr, str_len);
    strip_ip(new_str);
    return new_str;
}


vector<string> fortran_str_array_2_vec(char *fstr, int str_len, int array_len)
{
	vector<string> str_vec;
	str_vec.reserve(array_len);

	for (int ia=0; ia < array_len; ++ia)
	{
	    string new_str(fstr+ia*str_len, str_len);
            strip_ip(new_str);
	    str_vec.push_back(new_str);
	}
	return str_vec;
}

//template <class dataType>
//void read_twocol_ascii_to_map(map<string,dataType> &result, string filename, int header_lines, int data_col)
//{
//	map<string, dataType> result;
//	ifstream fin(filename);
//	if (!fin.good())
//		throw runtime_error("could not open file " + filename + " for reading");
//	string line;
//	dataType value;
//	vector<string> tokens;
//	for (int i = 0; i < header_lines; i++)
//		getline(fin, line);
//	while (getline(fin, line))
//	{
//		strip_ip(line);
//		if (line.at(0) == '#')
//			continue;
//		tokens.clear();
//		tokenize(line, tokens, "\t\r, ");
//		//only use the first two columns of file
//		if (tokens.size() < data_col + 1)
//			throw runtime_error("not enough entries on line :" + line);
//		convert_ip(tokens[data_col], value);
//		result[tokens[0]] = value;
//	}
//	return;
//}

map<string, double> read_twocol_ascii_to_map(string filename, int header_lines, int data_col)
{
	map<string, double> result;
	ifstream fin(filename);
	if (!fin.good())
		throw runtime_error("could not open file " + filename + " for reading");
	string line;
	double value;
	vector<string> tokens;
	for (int i = 0; i < header_lines; i++)
		getline(fin, line);
	while (getline(fin, line))
	{
		strip_ip(line);
		upper_ip(line);
		if ((line.size() == 0) || (line.at(0) == '#'))
			continue;
		tokens.clear();
		tokenize(line, tokens,"\t\r, ");
		//only use the first two columns of file
		if (tokens.size() < data_col + 1)
			throw runtime_error("not enough entries on line :" + line);
		convert_ip(tokens[data_col], value);
		result[tokens[0]] = value;
	}
	fin.close();
	return result;
}

vector<string> read_onecol_ascii_to_vector(std::string filename)
{
	vector<string> result;
	ifstream fin(filename);
	if (!fin.good())
		throw runtime_error("could not open file " + filename + " for reading");
	string line;
	vector<string> tokens;

	while (getline(fin, line))
	{
		strip_ip(line);
		upper_ip(line);
		if ((line.size() == 0) || (line.at(0) == '#'))
			continue;
		tokens.clear();
		tokenize(line, tokens, "\t\r, ");
		result.push_back(tokens[0]);
	}
	fin.close();
	return result;
}

void read_res(string &res_filename, Observations &obs)
{
	map<string, double> result;
	ifstream fin(res_filename);
	if (!fin.good())
		throw runtime_error("could not open residuals file " + res_filename + " for reading");
	vector<string> tokens;
	string line,name;
	double value;
	bool found = false;
	int mod_idx;
	while (getline(fin, line))
	{
		strip_ip(line);
		upper_ip(line);
		if (line.find("MODELLED") != std::string::npos)
		{
			tokens.clear();
			tokenize(line, tokens);
			mod_idx = find(tokens.begin(), tokens.end(), "MODELLED") - tokens.begin();
			found = true;
			break;
		}
	}

	if (!found)
		throw runtime_error("didn't find header line with 'modelled' in res file");
	vector<string> extra;
	vector<string> visited;
	vector<string> obs_names = obs.get_keys();
	while (getline(fin,line))
	{
		strip_ip(line);
		tokens.clear();
		tokenize(line, tokens);
		name = upper_cp(tokens[0]);
		convert_ip(tokens[mod_idx], value);
		if (find(obs_names.begin(), obs_names.end(), name) == obs_names.end())
			extra.push_back(name);
		else
		{
			obs[name] = value;
			visited.push_back(name);
		}

	}
	stringstream missing;
	missing << "the following obs were not found in the residual file: ";
	int i = 0;
	for (auto &oname : obs)
	{
		if (find(visited.begin(), visited.end(), oname.first) == visited.end())
		{
			missing << oname.first << ' ';
			i++;
			if (i % 5 == 0) missing << endl;
		}
	}
	if (i > 0)
		throw runtime_error(missing.str());

	if (extra.size() > 0)
	{
		stringstream ss;
		ss << "extra obs found res file...ignoring: ";
		int i = 0;
		for (auto &n : extra)
		{
			ss << n << ' ';
			i++;
			if (i % 5 == 0) ss << endl;
		}
		cout << ss.str();
	}
}


void read_par(ifstream &fin, Parameters &pars)
{
	string line;
	string name;
	double value;
	vector<string> tokens;
	getline(fin, line);
	while (getline(fin, line))
	{
		strip_ip(line);
		tokens.clear();
		tokenize(line, tokens);
		name = tokens[0];
		convert_ip(tokens[1], value);
		pars[name] = value;
	}
}

bool check_exist_in(std::string filename)
{
	std::ifstream f(filename.c_str());
	if (f.good())
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool check_exist_out(std::string filename)
{
	std::ofstream f(filename.c_str());
	if (f.good())
	{
		return true;
	}
	else
	{
		return false;
	}
}

thread_flag::thread_flag(bool _flag)
{
	flag = _flag;
}

bool thread_flag::set(bool f)
{
	std::lock_guard<std::mutex> lock(m);
	flag = f;
	return flag;
}
bool thread_flag::get()
{
	std::lock_guard<std::mutex> lock(m);
	if (flag) return true;
	else return false;
}

void thread_exceptions::add(std::exception_ptr ex_ptr)
{
	std::lock_guard<std::mutex> lock(m);
	shared_exception_vec.push_back(ex_ptr);
}

void  thread_exceptions::rethrow()
{
	std::lock_guard<std::mutex> lock(m);
	for (auto &iex : shared_exception_vec)
	{
		std::rethrow_exception(iex);
	}
}

bool read_binary(const string &filename, vector<string> &row_names, vector<string> &col_names, Eigen::SparseMatrix<double> &matrix)
{

	ifstream in;
	in.open(filename.c_str(), ifstream::binary);
	if (!in.good())
	{
		stringstream ss;
		ss << "Mat::from_binary() error opening binary file " << filename << " for reading";
		throw runtime_error(ss.str());
	}

	row_names.clear();
	col_names.clear();
	matrix.resize(0, 0);

	int n_par;
	int n_nonzero;
	int n_obs_and_pi;
	int i, j;
	unsigned int n;
	double data;
	char col_name[12];
	char row_name[20];

	// read header
	in.read((char*)&n_par, sizeof(n_par));
	in.read((char*)&n_obs_and_pi, sizeof(n_obs_and_pi));
	
	bool is_new_format = false;
	if (n_par > 0)
	{
		is_new_format = true;
		char col_name[200];
		char row_name[200];


		if (n_par > 100000000)
			throw runtime_error("pest_utils::read_binary() failed sanity check: npar > 100 mil");

		////read number nonzero elements in jacobian (observations + prior information)
		in.read((char*)&n_nonzero, sizeof(n_nonzero));

		if ((n_par == 0) || (n_obs_and_pi == 0) || (n_nonzero == 0))
		{
			throw runtime_error("pest_utils::read_binary() npar, nobs and/or nnz is zero");
		}

		cout << "reading " << n_nonzero << " elements, " << n_obs_and_pi << " rows, " << n_par << " columns" << endl;

		// record current position in file
		streampos begin_sen_pos = in.tellg();

		//advance to parameter names section
		in.seekg(n_nonzero*(sizeof(double) + sizeof(int) + sizeof(int)), ios_base::cur);

		//read parameter names
		for (int i_rec = 0; i_rec<n_par; ++i_rec)
		{
			in.read(col_name, 200);
			string temp_col = string(col_name, 200);
			pest_utils::strip_ip(temp_col);
			pest_utils::upper_ip(temp_col);
			col_names.push_back(temp_col);
		}
		//read observation and Prior info names
		for (int i_rec = 0; i_rec<n_obs_and_pi; ++i_rec)
		{
			in.read(row_name, 200);
			string temp_row = pest_utils::strip_cp(string(row_name, 200));
			pest_utils::upper_ip(temp_row);
			row_names.push_back(temp_row);
		}

		//return to sensitivity section of file
		in.seekg(begin_sen_pos, ios_base::beg);

		// read matrix
		std::vector<Eigen::Triplet<double> > triplet_list;
		triplet_list.reserve(n_nonzero);
		for (int i_rec = 0; i_rec<n_nonzero; ++i_rec)
		{
			in.read((char*)&(i), sizeof(i));
			in.read((char*)&(j), sizeof(j));
			in.read((char*)&(data), sizeof(data));

			if ((i >= n_obs_and_pi) || (i < 0))
				cout << "invalid 'i':" << i << " at " << n << " data:" << data << " j: " << j << endl;
			if ((j >= n_par) || (j < 0))
				cout << "invalid 'j':" << j << " at " << n << " data:" << data << " i: " << i << endl;
			triplet_list.push_back(Eigen::Triplet<double>(i, j, data));
		}
		matrix.resize(n_obs_and_pi, n_par);
		matrix.setZero();
		matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
		in.close();

	}
	else
	{
		char col_name[12];
		char row_name[20];
		n_par = -n_par;
		n_obs_and_pi = -n_obs_and_pi;

		if (n_par > 100000000)
			throw runtime_error("pest_utils::read_binary() failed sanity check: npar > 100 mil");

		////read number nonzero elements in jacobian (observations + prior information)
		in.read((char*)&n_nonzero, sizeof(n_nonzero));

		if ((n_par == 0) || (n_obs_and_pi == 0) || (n_nonzero == 0))
		{
			throw runtime_error("pest_utils::read_binary() npar, nobs and/or nnz is zero");
		}

		cout << "reading " << n_nonzero << " elements, " << n_obs_and_pi << " rows, " << n_par << " columns" << endl;

		// record current position in file
		streampos begin_sen_pos = in.tellg();

		//advance to parameter names section
		in.seekg(n_nonzero*(sizeof(double) + sizeof(int)), ios_base::cur);

		//read parameter names
		for (int i_rec = 0; i_rec<n_par; ++i_rec)
		{
			in.read(col_name, 12);
			string temp_col = string(col_name, 12);
			pest_utils::strip_ip(temp_col);
			pest_utils::upper_ip(temp_col);
			col_names.push_back(temp_col);
		}
		//read observation and Prior info names
		for (int i_rec = 0; i_rec<n_obs_and_pi; ++i_rec)
		{
			in.read(row_name, 20);
			string temp_row = pest_utils::strip_cp(string(row_name, 20));
			pest_utils::upper_ip(temp_row);
			row_names.push_back(temp_row);
		}

		//return to sensitivity section of file
		in.seekg(begin_sen_pos, ios_base::beg);

		// read matrix
		std::vector<Eigen::Triplet<double> > triplet_list;
		triplet_list.reserve(n_nonzero);
		for (int i_rec = 0; i_rec<n_nonzero; ++i_rec)
		{
			in.read((char*)&(n), sizeof(n));
			n = n - 1;
			in.read((char*)&(data), sizeof(data));
			j = int(n / (n_obs_and_pi)); // column index
			i = (n - n_obs_and_pi * j) % n_obs_and_pi;  //row index
			if ((i >= n_obs_and_pi) || (i < 0))
				cout << "invalid 'i':" << i << " at " << n << " data:" << data << " j: " << j << endl;
			if ((j >= n_par) || (j < 0))
				cout << "invalid 'j':" << j << " at " << n << " data:" << data << " i: " << i << endl;
			triplet_list.push_back(Eigen::Triplet<double>(i, j, data));
		}
		matrix.resize(n_obs_and_pi, n_par);
		matrix.setZero();
		matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
		in.close();
	}
	return is_new_format;
}

bool read_binary(const string &filename, vector<string> &row_names, vector<string> &col_names, Eigen::MatrixXd &matrix)
{
	Eigen::SparseMatrix<double> smatrix;
	bool is_new_format = read_binary(filename, row_names, col_names, smatrix);
	matrix = Eigen::MatrixXd(smatrix);
	return is_new_format;
}

void save_binary(const string &filename, const vector<string> &row_names, const vector<string> &col_names, const Eigen::SparseMatrix<double> &matrix)
{
	//check row name and col name lengths
	int mx_rlen = 0, mx_clen = 0, v;
	for (auto& n : row_names)
	{
		v = n.length();
		mx_rlen = max(mx_rlen, v);
	}
		
	for (auto& n : col_names)
	{
		v = n.length();
		mx_clen = max(mx_clen, v);
	}
		
	//mx_rlen = max_element(row_names.begin(), row_names.end()) - row_names.begin();
	//mx_clen = max_element(col_names.begin(), col_names.end()) - col_names.begin();
	if ((mx_rlen > 20) || (mx_clen > 12))
		save_binary_extfmt(filename, row_names, col_names, matrix);
	else
		save_binary_orgfmt(filename, row_names, col_names, matrix);
}

void save_binary_extfmt(const string &filename, const vector<string> &row_names, const vector<string> &col_names, const Eigen::SparseMatrix<double> &matrix)
{
	ofstream jout(filename, ios::out | ios::binary);
	int n_par = col_names.size();
	int n_obs_and_pi = row_names.size();
	int n;
	int tmp;
	double data;
	char par_name[200];
	char obs_name[200];

	// write header
	tmp = n_par;
	jout.write((char*)&tmp, sizeof(tmp));
	tmp = n_obs_and_pi;
	jout.write((char*)&tmp, sizeof(tmp));

	//write number nonzero elements in jacobian (includes prior information)
	n = matrix.nonZeros();
	jout.write((char*)&n, sizeof(n));

	//write matrix
	n = 0;
	map<string, double>::const_iterator found_pi_par;
	map<string, double>::const_iterator not_found_pi_par;

	Eigen::SparseMatrix<double> matrix_T(matrix);
	matrix_T.transpose();
	for (int icol = 0; icol<matrix.outerSize(); ++icol)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(matrix_T, icol); it; ++it)
		{
			data = it.value();
			n = it.row() - 1;
			jout.write((char*) &(n), sizeof(n));
			n = it.col() - 1;
			jout.write((char*) &(n), sizeof(n));

			jout.write((char*) &(data), sizeof(data));
		}
	}
	//save parameter names
	for (vector<string>::const_iterator b = col_names.begin(), e = col_names.end();
		b != e; ++b) {
		string l = pest_utils::lower_cp(*b);
		pest_utils::string_to_fortran_char(l, par_name, 200);
		jout.write(par_name, 200);
	}

	//save observation and Prior information names
	for (vector<string>::const_iterator b = row_names.begin(), e = row_names.end();
		b != e; ++b) {
		string l = pest_utils::lower_cp(*b);
		pest_utils::string_to_fortran_char(l, obs_name, 200);
		jout.write(obs_name, 200);
	}
	//save observation names (part 2 prior information)
	jout.close();
}


void save_binary_orgfmt(const string &filename, const vector<string> &row_names, const vector<string> &col_names, const Eigen::SparseMatrix<double> &matrix)
{
	ofstream jout(filename, ios::out | ios::binary);
	int n_par = col_names.size();
	int n_obs_and_pi = row_names.size();
	int n;
	int tmp;
	double data;
	char par_name[12];
	char obs_name[20];

	// write header
	tmp = -n_par;
	jout.write((char*)&tmp, sizeof(tmp));
	tmp = -n_obs_and_pi;
	jout.write((char*)&tmp, sizeof(tmp));

	//write number nonzero elements in jacobian (includes prior information)
	n = matrix.nonZeros();
	jout.write((char*)&n, sizeof(n));

	//write matrix
	n = 0;
	map<string, double>::const_iterator found_pi_par;
	map<string, double>::const_iterator not_found_pi_par;

	//cout << matrix.toDense() << endl;

	Eigen::SparseMatrix<double> matrix_T(matrix);
	matrix_T.transpose();
	for (int icol = 0; icol<matrix.outerSize(); ++icol)
	{
		for (Eigen::SparseMatrix<double>::InnerIterator it(matrix_T, icol); it; ++it)
		{
			data = it.value();
			//cout << icol << "," << col_names[icol] << " - " << it.row() << "," << row_names[it.row()] << ":" << data << endl;
			n = it.row() + 1 + it.col() * matrix_T.rows();
			jout.write((char*) &(n), sizeof(n));
			jout.write((char*) &(data), sizeof(data));
		}
	}
	//save parameter names
	for (vector<string>::const_iterator b = col_names.begin(), e = col_names.end();
		b != e; ++b) {
		string l = pest_utils::lower_cp(*b);
		pest_utils::string_to_fortran_char(l, par_name, 12);
		jout.write(par_name, 12);
	}

	//save observation and Prior information names
	for (vector<string>::const_iterator b = row_names.begin(), e = row_names.end();
		b != e; ++b) {
		string l = pest_utils::lower_cp(*b);
		pest_utils::string_to_fortran_char(l, obs_name, 20);
		jout.write(obs_name, 20);
	}
	//save observation names (part 2 prior information)
	jout.close();
}


} // end of namespace pest_utils





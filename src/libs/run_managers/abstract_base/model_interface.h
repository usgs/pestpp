#ifndef MODEL_INTERFACE_H_
#define MODEL_INTERFACE_H_

#include <vector>
#include <string>
#include "Transformable.h"
#include "utilities.h"

using namespace std;

class TemplateFile {
public:
	TemplateFile(string _tpl_filename): tpl_filename(_tpl_filename),line_num(0) { ; }
	set<string> parse_and_check(ofstream& f_rec);
	void write_input_file(ofstream &f_rec, const string& input_filename, Parameters& pars);
	void throw_tpl_error(ofstream& f_rec, const string& message, int lnum=0, bool warn=false);


private:
	int line_num;
	string marker;
	string tpl_filename;
	set<string> names;
	map<string, pair<int, int>> parse_tpl_line(ofstream& f_rec, const string& line);
	string cast_to_fixed_len_string(ofstream& f_rec, int size, double value, string& name);
	string read_line(ofstream& f_rec, ifstream& f);
	void prep_tpl_file_for_reading(ofstream& f_rec, ifstream& f);
	set<string> get_names(ofstream& f_rec, ifstream& f);
	vector<int> find_all_marker_indices(const string& line);

};


class InstructionFile {
public:
	InstructionFile(string _ins_filename): ins_filename(_ins_filename), line_num(0) { ; }
	set<string> parse_and_check(ofstream& f_rec);
	map<string, double> read_output_file(ofstream& f_rec, const string& output_filename);

private:
	int line_num;
	string marker;
	string ins_filename;
	set<string> names;
	double execute_fixed(ofstream& f_rec, ifstream& f);
	double execute_semi(ofstream& f_rec, ifstream& f);
	double execute_free(ofstream& f_rec, ifstream& f);
	void execute_primary(ofstream& f_rec, ifstream& f);
	void execute_secondary(ofstream& f_rec, ifstream& f);
	void execute_w(ofstream& f_rec, ifstream& f);
	void execute_dum(ofstream& f_rec, ifstream& f);

};


class ModelInterface{
public:
	ModelInterface();
	ModelInterface(vector<string> _tplfile_vec,vector<string> _inpfile_vec, vector<string> _insfile_vec, vector<string> _outfile_vec,vector<string> _comline_vec);
	void throw_mio_error(string base_message);
	void run(Parameters* pars, Observations* obs);
	void run(pest_utils::thread_flag* terminate, pest_utils::thread_flag* finished,
		pest_utils::thread_exceptions *shared_execptions,
		vector<string> &par_name_vec, vector<double> &par_values,
		vector<string> &obs_name_vec, vector<double> &obs_vec);
	void run(pest_utils::thread_flag* terminate, pest_utils::thread_flag* finished,
		pest_utils::thread_exceptions *shared_execptions,
		Parameters* par, Observations* obs);


	void initialize(vector<string> &_par_name_vec, vector<string> &_obs_name_vec);
	void initialize(vector<string> _tplfile_vec, vector<string> _inpfile_vec,
		vector<string> _insfile_vec, vector<string> _outfile_vec, vector<string> _comline_vec,
		vector<string> &_par_name_vec, vector<string> &_obs_name_vec);
	void finalize();
	~ModelInterface();
	bool get_initialized(){ return initialized; }
private:

	void set_files();
	void check();

	bool initialized;
	int ifail;
	vector<string> par_name_vec;
	vector<string> obs_name_vec;
	vector<string> tplfile_vec;
	vector<string> inpfile_vec;
	vector<string> outfile_vec;
	vector<string> insfile_vec;
	vector<string> comline_vec;

	vector<double> par_vals;
	vector<double> obs_vals;

};

#endif /* MODEL_INTERFACE_H_ */


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
	set<string> parse_and_check();
	void write_input_file(const string& input_filename, Parameters& pars);
	void throw_tpl_error(const string& message, int lnum=0, bool warn=false);


private:
	int line_num;
	string marker;
	string tpl_filename;
	set<string> names;
	map<string, pair<int, int>> parse_tpl_line(const string& line);
	string cast_to_fixed_len_string(int size, double value, string& name);
	string read_line(ifstream& f);
	void prep_tpl_file_for_reading(ifstream& f);
	set<string> get_names(ifstream& f);
	vector<int> find_all_marker_indices(const string& line);

};


class InstructionFile {
public:
	InstructionFile(string _ins_filename): ins_filename(_ins_filename), line_num(0) { ; }
	set<string> parse_and_check();
	map<string, double> read_output_file(const string& output_filename);

private:
	int line_num;
	string marker;
	string ins_filename;
	set<string> names;
	double execute_fixed(ifstream& f);
	double execute_semi(ifstream& f);
	double execute_free(ifstream& f);
	void execute_primary(ifstream& f);
	void execute_secondary(ifstream& f);
	void execute_w(ifstream& f);
	void execute_dum(ifstream& f);

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
	
	bool initialized;
	int ifail;
	vector<string> par_name_vec;
	vector<string> obs_name_vec;
	vector<string> tplfile_vec;
	vector<string> inpfile_vec;
	vector<string> outfile_vec;
	vector<string> insfile_vec;
	vector<string> comline_vec;
	vector<TemplateFile> templatefiles;

	vector<double> par_vals;
	vector<double> obs_vals;

};

#endif /* MODEL_INTERFACE_H_ */


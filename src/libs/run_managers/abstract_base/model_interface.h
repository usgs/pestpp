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
	Parameters write_input_file(const string& input_filename, Parameters& pars);
	void throw_tpl_error(const string& message, int lnum=0, bool warn=false);


private:
	int line_num;
	string marker;
	string tpl_filename;
	set<string> names;
	map<string, pair<int, int>> parse_tpl_line(const string& line);
	string cast_to_fixed_len_string(int size, double value, string& name);
	string read_line(ifstream& f_tpl);
	void prep_tpl_file_for_reading(ifstream& f_tpl);
	set<string> get_names(ifstream& f);
	vector<int> find_all_marker_indices(const string& line);
	


};


class InstructionFile {
	
public:
	InstructionFile(string _ins_filename);
	set<string> parse_and_check();
	Observations read_output_file(const string& output_filename);

private:
	int ins_line_num, out_line_num;
	char marker;
	string ins_filename, last_out_line, last_ins_line;
	set<string> names;
	vector<pair<char, char>> obs_tags;
	pair<string, double> execute_fixed(const string& token, string& line, ifstream& f_out);
	pair<string, double> execute_semi(const string& token, string& line, ifstream& f_out);
	pair<string, double> execute_free(const string& token, string& line, ifstream& f_out);
	void execute_primary(const string& token, string& line, ifstream& f_out);
	void execute_secondary(const string& token, string& line, ifstream& f_out);
	void execute_whitespace(const string& token, string& line, ifstream& f_out);
	void execute_line_advance(const string& token, string& line, ifstream& f_out);
	void prep_ins_file_for_reading(ifstream& f_ins);
	string read_ins_line(ifstream& f_ins);
	string read_out_line(ifstream& f_out);
	void throw_ins_error(const string& message, int ins_lnum = 0, int out_lnum=0, bool warn = false);
	void parse_obs_name_from_token(const string& token, string& obs_name);
};


class ModelInterface{
public:
	ModelInterface();
	ModelInterface(vector<string> _tplfile_vec,vector<string> _inpfile_vec, vector<string> _insfile_vec, vector<string> _outfile_vec,vector<string> _comline_vec);
	void throw_mio_error(string base_message);
	void run(Parameters* pars, Observations* obs);
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
	vector<InstructionFile> instructionfiles;
	vector<double> par_vals;
	vector<double> obs_vals;

};

#endif /* MODEL_INTERFACE_H_ */


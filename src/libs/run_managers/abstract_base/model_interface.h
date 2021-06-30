#ifndef MODEL_INTERFACE_H_
#define MODEL_INTERFACE_H_

#include <vector>
#include <string>
#include <unordered_set>
#include <mutex>
#include "Transformable.h"
#include "utilities.h"
#include "Pest.h"

using namespace std;

class TemplateFile {
public:
	static vector<int> find_all_marker_indices(const string& line, const string& marker);
	TemplateFile(string _tpl_filename, bool _fill_zeros=false): tpl_filename(_tpl_filename),line_num(0),
	fill_zeros(_fill_zeros){ ; }
	unordered_set<string> parse_and_check();
	Parameters write_input_file(const string& input_filename, Parameters& pars);
	void throw_tpl_error(const string& message, int lnum=0, bool warn=false);
	void set_fill_zeros(bool _flag) { fill_zeros = _flag; }
	string get_tpl_filename() { return tpl_filename; }
private:
	int line_num;
	string marker;
	string tpl_filename;
	vector<pair<string, pair<int, int>>> parse_tpl_line(const string& line);
	string cast_to_fixed_len_string(int size, double value, string& name);
	string read_line(ifstream& f_tpl);
	void prep_tpl_file_for_reading(ifstream& f_tpl);
	unordered_set<string> get_names(ifstream& f);
	bool fill_zeros;
	
};

class ThreadedTemplateProcess {
public:
	ThreadedTemplateProcess(vector<string> _tplfile_vec, vector<string> _inpfile_vec, bool _fill) : 
		tplfile_vec(_tplfile_vec), inpfile_vec(_inpfile_vec), fill(_fill) {;};
	void work(int tid, vector<int>& tpl_idx, Parameters pars, Parameters& pro_pars);
private:
	vector<string> tplfile_vec;
	vector<string> inpfile_vec;
	bool fill;
	mutex par_lock, idx_lock;
};

class ThreadedInstructionProcess {
public:
	ThreadedInstructionProcess(vector<string> _insfile_vec, vector<string> _outfile_vec) :
		insfile_vec(_insfile_vec), outfile_vec(_outfile_vec){;};
	void work(int tid, vector<int>& ins_idx, Observations& obs, string additional_ins_delims);
private:
	vector<string> insfile_vec;
	vector<string> outfile_vec;
	mutex obs_lock, idx_lock;
};


class InstructionFile {
	
public:
	InstructionFile(string _ins_filename, string _additional_delimiters="");
	unordered_set<string> parse_and_check();
	Observations read_output_file(const string& output_filename);
	void set_additional_delimiters(string delims) { additional_delimiters = delims; }
private:
	int ins_line_num, out_line_num;
	char marker;
	string ins_filename, last_out_line, last_ins_line;
	vector<pair<char, char>> obs_tags;
	pair<string, double> execute_fixed(const string& token, string& line, ifstream& f_out);
	pair<string, double> execute_semi(const string& token, string& line, ifstream& f_out);
	pair<string, double> execute_free(const string& token, string& line, ifstream& f_out);
	void execute_primary(const string& token, string& line, ifstream& f_out);
	bool execute_secondary(const string& token, string& line, ifstream& f_out,bool all_markers_so_far);
	void execute_whitespace(const string& token, string& line, ifstream& f_out);
	void execute_line_advance(const string& token, string& line, ifstream& f_out);
	void prep_ins_file_for_reading(ifstream& f_ins);
	string read_ins_line(ifstream& f_ins);
	string read_out_line(ifstream& f_out);
	void throw_ins_error(const string& message, int ins_lnum = 0, int out_lnum=0, bool warn = false);
	string parse_obs_name_from_token(const string& token);
	vector<string> tokenize_ins_line(const string& line);
	pair<string, pair<int, int>> parse_obs_instruction(const string& token, const string& close_tag);
	string additional_delimiters;
	
	void tokenize(const std::string& str, vector<string>& tokens, const std::string& delimiters, const bool trimEmpty=true, int mx_tokens=-1);
	

};


class ModelInterface{
public:
	ModelInterface() { ; }
	//ModelInterface(Pest* _pest_scenario_ptr) { pest_scenario_ptr = _pest_scenario_ptr; }
	ModelInterface(vector<string> _tplfile_vec, vector<string> _inpfile_vec, vector<string>
		_insfile_vec, vector<string> _outfile_vec, vector<string> _comline_vec) :
		insfile_vec(_insfile_vec), outfile_vec(_outfile_vec), tplfile_vec(_tplfile_vec),
		inpfile_vec(_inpfile_vec), comline_vec(_comline_vec), fill_tpl_zeros(false), additional_ins_delimiters(""),
		num_threads(1)
		{;}
	void throw_mio_error(string base_message);
	void run(Parameters* pars, Observations* obs);
	void run(pest_utils::thread_flag* terminate, pest_utils::thread_flag* finished,
		exception_ptr& eptr,
		Parameters* par, Observations* obs);
	void check_io_access();
	void check_tplins(const vector<string> &par_names, const vector<string> &obs_names);
	void set_additional_ins_delimiters(string delims) { additional_ins_delimiters = delims; }
	void set_fill_tpl_zeros(bool _flag) { fill_tpl_zeros = _flag; }
	void set_num_threads(int _num_threads) { num_threads = _num_threads; }

private:
	int num_threads;
	//Pest* pest_scenario_ptr;
	vector<TemplateFile> templatefiles;
	vector<InstructionFile> instructionfiles;
	vector<string> insfile_vec; 
	vector<string> inpfile_vec; 
	vector<string> outfile_vec; 
	vector<string> tplfile_vec; 
	vector<string> comline_vec; 
	bool fill_tpl_zeros;
	string additional_ins_delimiters;

	void write_input_files(Parameters *pars_ptr);
	void read_output_files(Observations *obs_ptr);
	void remove_existing();

};

#endif /* MODEL_INTERFACE_H_ */


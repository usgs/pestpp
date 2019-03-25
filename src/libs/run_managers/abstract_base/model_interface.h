#ifndef MODEL_INTERFACE_H_
#define MODEL_INTERFACE_H_

#include <vector>
#include <string>
#include "Transformable.h"
#include "utilities.h"

using namespace std;

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


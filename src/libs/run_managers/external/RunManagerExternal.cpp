#include <iostream>
#include "RunManagerExternal.h"
#include "utilities.h"

using namespace std;

RunManagerExternal::RunManagerExternal(const vector<string> _comline_vec,
	const vector<string> _tplfile_vec, const vector<string> _inpfile_vec,
	const vector<string> _insfile_vec, const vector<string> _outfile_vec,
	const string &stor_filename, const string &_ext_filename,
	const string &_exi_filename,int _max_n_failure)
	: RunManagerAbstract(_comline_vec, _tplfile_vec, _inpfile_vec,
	_insfile_vec, _outfile_vec, stor_filename, _max_n_failure),
	ext_filename(_ext_filename), exi_filename(_exi_filename)
{
	cout << "              starting external run manager ..." << endl << endl;
}

void RunManagerExternal::run()
{
	vector<int> waiting_run_ids = get_outstanding_run_ids();

	//if there are outstanding runs that need to be made
	//exit so the external run manager can be involked
	if (!waiting_run_ids.empty())
	{
		ofstream fout_ext(ext_filename);
		fout_ext << get_run_filename() << endl;
		fout_ext << max_n_failure << endl;
		fout_ext.close();

		cout << endl << endl;
		try{
			ifstream fin_exi(exi_filename);
			string line;
			getline(fin_exi, line);
			if (!line.empty())
			{
				pest_utils::strip_ip(line);
				cout << "  External run mananager running script \"" << line << "\"" << endl << endl;
				system(line.c_str());
			}
		}
		catch(...){

		}

		cout << endl;
		cout << "  External run mananager involked.  Leaving PEST++ ..." << endl;
		cout.flush();
		exit(0);
	}


}

RunManagerExternal::~RunManagerExternal()
{
}

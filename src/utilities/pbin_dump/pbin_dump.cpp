/*


This file is part of PEST++.

PEST++ is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PEST++ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PEST++.  If not, see<http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "RunStorage.h"

using namespace std;


void usage(ostream &fout)
{
	fout << "--------------------------------------------------------" << endl;
	fout << "usage:" << endl << endl;
	fout << "  pbin_dump.exe bin_file  ascii_file" << endl << endl;
	fout << " where:" << endl;
	fout << "  bin_file:    pest++ binary model run file to be read" << endl;
	fout << "  ascii_file:  name of ascii file where results will be" << endl;
	fout << "               saved" << endl;
	fout << "--------------------------------------------------------" << endl;
}


int main(int argc, char* argv[])
{
	string in_filename = argv[1];
	string out_filename = argv[2];
	ofstream fout(out_filename);
	RunStorage rs("");
	rs.init_restart(in_filename);

	vector<string> par_name_vec =  rs.get_par_name_vec();
	vector<string> obs_name_vec = rs.get_obs_name_vec();

	int n_runs = rs.get_nruns();
	cout << "processing " << n_runs << " runs"<< endl;

	fout << "#####################################################################################################" << endl;
	fout << "Header information" << endl;
	fout << "number of runs = " << n_runs << endl;
	fout << "parameter names:" << endl;
	const vector<string> &par_names_vec = rs.get_par_name_vec();
	for (const auto &name : par_names_vec)
	{
		fout << "  " << name << endl;
	}
	fout << "observation names:" << endl;
	const vector<string> &obs_names_vec = rs.get_obs_name_vec();
	for (const auto &name : obs_names_vec)
	{
		fout << "  " << name << endl;
	}
	int status = 0;
	string info_text;
	double info_value;
	vector<double> pars_vec;
	vector<double> obs_vec;
	fout << "#####################################################################################################" << endl;

	for (int i = 0; i < n_runs; ++i)
	{
		fout << "run_id = " << i << endl;
		rs.get_info(i, status, info_text, info_value);
		fout << "status = " << status << endl;
		fout << "info text = " << info_text << endl;
		fout << "info value = " << info_value << endl;
		obs_vec.clear();
		pars_vec.clear();
		rs.get_run(i, pars_vec, obs_vec);
		fout << "parameter values:" << endl;
		for (const auto &ipar : pars_vec)
		{
			fout << "  " << ipar << endl;
		}
		fout << "observation values:" << endl;
		for (const auto &iobs : obs_vec)
		{
			fout << "  " <<  iobs << endl;
		}
		fout << "#####################################################################################################" << endl;
	}
	fout.close();
}

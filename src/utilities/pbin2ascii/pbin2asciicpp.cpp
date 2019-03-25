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


/*
NOTE:
	This project can be compiled under different configurations, and produces a different binary
	based on the configuration chosen.  To change configurations, you will need to go into the
	configuration manager and choose the configuration that matches your need.

	Default Configuration:
	pbin2ascii.exe - Choose the Debug/Release configuration


	pbin2par.exe - Choose the pbin2par-Debug/pbin2par-Release configuration

	pbin2ascmat.exe - Choose the pbin2ascmat-Debug/pbin2ascmat-Release configuration

*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "RunStorage.h"
#include "utilities.h"

#ifdef ASCIIMATRIX
	// show headers
	#define ASCIIMATRIX_PAR_HEADERS
	#define ASCIIMATRIX_RUN_HEADERS

	#define EXPECTED_ARG_COUNT 3
#else
	#ifndef PARFILES
		#define ASCII
	#endif
	#define EXPECTED_ARG_COUNT 4
	#define GENERATE_INDIVIDUAL_FILES
#endif

using namespace std;

void usage(string exeName, ostream &fout)
{

	fout << "--------------------------------------------------------" << endl;
	fout << "usage:" << endl << endl;
	fout << " " << exeName << " ext_file runs_file";
	if(EXPECTED_ARG_COUNT > 3)
	{
		fout << " par_prefix";
	}
	fout << endl << endl;

	fout << " where:" << endl;
	fout << "  ext_file:   name of the pest++ .ext file containing the" << endl;
	fout << "              binary run file name" << endl;
	fout << "  runs_file:  name of file to be created with list of" << endl;
	fout << "              required runs" << endl;
	if(EXPECTED_ARG_COUNT > 3)
	{
		fout << "  par_prefix: prefix used to create parameter files" << endl;
	}
	fout << "--------------------------------------------------------" << endl;
}


int main(int argc, char* argv[])
{
	string exeName = argv[0];
	exeName = exeName.substr(exeName.find_last_of("/\\") + 1);

	if(argc != EXPECTED_ARG_COUNT)
	{
		cerr << "Error: incorrect number of command line arguments" << endl << endl;
		usage(exeName, cerr);
		return 1;
	}

	string ext_filename = argv[1];
	string out_filename = argv[2];
	string base_par_filename;
	if(EXPECTED_ARG_COUNT > 3)
	{
		base_par_filename = argv[3];
	}


	string pbin_filename;
	string errorMessage;
	try
	{
		ifstream fin_ext(ext_filename);
		if(!fin_ext.good())
		{
			fin_ext.close();
			throw PestFileError(ext_filename);
		}
		getline(fin_ext, pbin_filename);
		pest_utils::strip_ip(pbin_filename);
		fin_ext.close();

		if(pbin_filename == "")
		{
			throw runtime_error("PEST binary file not specified in *.ext \"" + ext_filename + "\": line 1");
		}
	}
	catch (exception &e)
	{
		cerr << "Error processing external run manager *.ext \"" << ext_filename << "\"" << endl;
		cerr << e.what() << endl << endl;
		//usage(exeName, cerr);
		return 1;
	//	throw(e);
	}


	RunStorage rs("");
	try
	{
		rs.init_restart(pbin_filename);
	}
	catch (exception &e)
	{
		cerr << "Error processing PEST++ binary file \"" << pbin_filename << "\"" << endl;
		cerr << e.what() << endl;
		//usage(exeName, cerr);
		return 1;
		//throw(e);
	}

	vector<string> par_name_vec = rs.get_par_name_vec();
	vector<string> obs_name_vec = rs.get_obs_name_vec();

	ofstream fout(out_filename);

	int n_runs = rs.get_nruns();
	std::cout << "processing " << n_runs << " runs" << endl;

	int status = 0;
	string info_text;
	double info_value;
	vector<double> pars_vec;
	vector<double> obs_vec;
	bool headersPrinted = false;
	int maxParNameLength = -1;

	for (int i = 0; i < n_runs; ++i)
	{
		obs_vec.clear();
		pars_vec.clear();
		rs.get_info(i, status, info_text, info_value);
		if (status == 0)
		{
			rs.get_run(i, pars_vec, obs_vec);

#ifdef GENERATE_INDIVIDUAL_FILES
			stringstream par_file_name;
			par_file_name << base_par_filename << "_" << i << ".par";
			ofstream fout_par;
			fout_par.open(par_file_name.str());
#endif
			assert(par_name_vec.size() == pars_vec.size());
			size_t n_par = par_name_vec.size();

			if(maxParNameLength < 0)
			{
				maxParNameLength = 19;  //set the baseline
				for(size_t ipar = 0; ipar < n_par; ++ipar)
				{
					maxParNameLength = (maxParNameLength < (int)par_name_vec[ipar].length()) ? (int)par_name_vec[ipar].length() : maxParNameLength;
				}
				maxParNameLength++;
			}

#ifdef PARFILES
			fout_par << "single point" << endl;
#endif
#ifdef ASCIIMATRIX
	#ifdef ASCIIMATRIX_PAR_HEADERS
			if(!headersPrinted)
			{
				headersPrinted = true;
		#ifdef ASCIIMATRIX_RUN_HEADERS
				fout << setw(10) << left << "RUN";
		#endif
				for(size_t ipar = 0; ipar < n_par; ++ipar)
				{
					fout << setw(maxParNameLength) << par_name_vec[ipar] << " ";
				}
				fout << endl;
			}
	#endif
	#ifdef ASCIIMATRIX_RUN_HEADERS
			fout << setw(9) << left << i << " ";
	#endif
#endif

			for (size_t ipar = 0; ipar < n_par; ++ipar)
			{
#ifdef ASCII
				fout_par << setw(maxParNameLength) << par_name_vec[ipar] << " ";
				int prec = fout_par.precision(numeric_limits<double>::digits10 + 1);
				fout_par.precision(prec);
				fout_par << pars_vec[ipar] << endl;
#else
	#ifdef PARFILES
				fout_par << setw(maxParNameLength) << par_name_vec[ipar] << " ";
				fout_par << setprecision(16) << fixed << setw(20) << pars_vec[ipar];
				fout_par << setprecision(2) << fixed << setw(7) << 1.0;				//scalar
				fout_par << setprecision(2) << fixed << setw(7) << 0.0 << endl;		//offset
	#else
		#ifdef ASCIIMATRIX
				fout << setprecision(16) << fixed << setw(20) << pars_vec[ipar] << " ";
		#endif
	#endif
#endif
			}
#ifdef GENERATE_INDIVIDUAL_FILES
			fout_par.close();
			fout << i << " " << par_file_name.str() << endl;
#else
	#ifdef ASCIIMATRIX
			fout << endl;
	#endif
#endif
		}
	}
	fout.close();
}

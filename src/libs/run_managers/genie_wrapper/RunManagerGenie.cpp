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
#include "RunManagerGenie.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <iterator>
#include <cassert>
#include <cstring>
#include "Transformable.h"
#include "utilities.h"
#include "config_os.h"

using namespace std;
using namespace pest_utils;

const int RunManagerGenie::LEN_PARAMETER_NAME = 12+1;
const int RunManagerGenie::LEN_OBSERVATION_NAME = 20+1;

extern "C" {
	int GENIE_INTERFACE(int*, int*, char*,
		int*, int*, char*, char*, double*, double*,
		int*, int*, char*, char*, char*,
		char*, char*, char*, int*);
	int GENIE_KILL_GMAN(char*,char*);
}



RunManagerGenie::RunManagerGenie(const vector<string> _comline_vec,
	const vector<string> _tplfile_vec, const vector<string> _inpfile_vec,
	const vector<string> _insfile_vec, const vector<string> _outfile_vec,
	const string &stor_filename, const std::string &_host,
	const std::string &_id)
	: RunManagerAbstract(_comline_vec, _tplfile_vec, _inpfile_vec,
		_insfile_vec, _outfile_vec, stor_filename),
		host(_host), id(_id)
{
	cout << "                starting GENIE run manager ..." << endl;
	cout << "   	    developed by Chris Muffels of SSPA" << endl << endl;
}

void RunManagerGenie::run()
{
	int nexec = comline_vec.size();
	const vector<string> &par_name_vec = file_stor.get_par_name_vec();
	const vector<string> &obs_name_vec = file_stor.get_obs_name_vec();
	int npar = par_name_vec.size();
	int nobs = obs_name_vec.size();
	int ntpl = tplfile_vec.size();
	int nins = insfile_vec.size();

	stringstream apar;
	stringstream aobs;
	stringstream execnames;
	stringstream tplfle;
	stringstream infle;
	stringstream insfle;
	stringstream outfle;

	 vector<double> par_val;
	 vector<double> obs_val;
	 vector<int> ifail;

	// This is necessary to support restart as some run many already be complete
	vector<int> run_id_vec = get_outstanding_run_ids();
	int nruns = run_id_vec.size();

	 par_val.clear();
	 obs_val.resize(nruns*nobs);
	 ifail.resize(nruns);
	for (int i_run : run_id_vec)
	{
	     Parameters tmp_pars;
	     file_stor.get_parameters(i_run, tmp_pars);
		 vector<double> tmp_vec = tmp_pars.get_data_vec(par_name_vec);
		 par_val.insert(par_val.end(), tmp_vec.begin(), tmp_vec.end());
	}

	vector<string> par_name_vec_lwr = par_name_vec;
	for (auto &ipar : par_name_vec_lwr)
	{
		lower_ip(ipar);
	}
	vector<string> obs_name_vec_lwr = obs_name_vec;
	for (auto &iobs : obs_name_vec_lwr)
	{
		lower_ip(iobs);
	}
	std::copy(par_name_vec_lwr.begin(), par_name_vec_lwr.end(),std::ostream_iterator<std::string>(apar,"\n"));
	std::copy(obs_name_vec_lwr.begin(), obs_name_vec_lwr.end(),std::ostream_iterator<std::string>(aobs,"\n"));
	std::copy(comline_vec.begin(), comline_vec.end(),std::ostream_iterator<std::string>(execnames,"\n"));
	std::copy(tplfile_vec.begin(), tplfile_vec.end(),std::ostream_iterator<std::string>(tplfle,"\n"));
	std::copy(inpfile_vec.begin(), inpfile_vec.end(),std::ostream_iterator<std::string>(infle,"\n"));
	std::copy(insfile_vec.begin(), insfile_vec.end(),std::ostream_iterator<std::string>(insfle,"\n"));
	std::copy(outfile_vec.begin(), outfile_vec.end(),std::ostream_iterator<std::string>(outfle,"\n"));

	#ifdef OS_WIN
	GENIE_INTERFACE(&nruns, &nexec, String2CharPtr(execnames.str()).get_char_ptr(), &npar, &nobs,
	String2CharPtr(apar.str()).get_char_ptr(), String2CharPtr(aobs.str()).get_char_ptr(),
					&par_val[0], &obs_val[0], &ntpl, &nins, String2CharPtr(tplfle.str()).get_char_ptr(),
					String2CharPtr(infle.str()).get_char_ptr(), String2CharPtr(insfle.str()).get_char_ptr(),
					String2CharPtr(outfle.str()).get_char_ptr(),
					String2CharPtr(host).get_char_ptr(), String2CharPtr(id).get_char_ptr(), &ifail[0]);
	#endif
	#ifdef OS_LINUX
	throw(PestError("Error: Genie run manager is not supported under linux"));
	#endif

	total_runs += nruns;
	Parameters pars;
	Observations obs;
	for (int i=0; i<nruns; ++i)
	{
		if (ifail[i] == 0)
		{
			int run_id = run_id_vec[i];
			pars.clear();
			vector<double> i_par_vec(par_val.begin() + i*npar, par_val.begin() + (i + 1)*npar);
			pars.insert(par_name_vec, i_par_vec);
			obs.clear();
			vector<double> i_obs_vec(obs_val.begin() + i*nobs, obs_val.begin() + (i + 1)*nobs);
			obs.insert(obs_name_vec, i_obs_vec);
			file_stor.update_run(run_id, pars, obs);
		}
		else
		{
			file_stor.set_run_nfailed(i, max_n_failure);
		}
	}
	/*if (init_run_obs.size() == 0)
		int status = file_stor.get_observations(0, init_run_obs);*/
	if (init_sim.size() == 0)
	{
		vector<double> pars;
		int status = file_stor.get_run(0, pars, init_sim);
	}
}

RunManagerGenie::~RunManagerGenie(void)
{
	//close GMAN run manager
	cout << "Closing GMAN run manager..." << endl;
	//GENIE_KILL_GMAN(String2CharPtr(id).get_char_ptr(), String2CharPtr(host).get_char_ptr());

	free_memory();
}


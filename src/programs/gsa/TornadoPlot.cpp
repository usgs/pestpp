#include <vector>
#include <cassert>
#include <numeric>
#include <math.h>
#include <errno.h>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <regex>
#include "TornadoPlot.h"
#include "Transformable.h"
#include "RunManagerAbstract.h"
#include "ParamTransformSeq.h"
#include "ModelRunPP.h"
#include "utilities.h"
#include "FileManager.h"

using namespace std;
using namespace pest_utils;
TornadoPlot::TornadoPlot(const vector<string> &_adj_par_name_vec, const Parameters &_fixed_pars, const Parameters &_init_pars,
	const Parameters &_lower_bnd, const Parameters &_upper_bnd, const set<string> &_log_trans_pars,
	ParamTransformSeq *_base_partran_seq_ptr, const std::vector<std::string> &_obs_name_vec,
	FileManager *file_manager_ptr, const ObservationInfo *_obs_info_ptr,
	bool _calc_obs_sen)
	//: GsaAbstractBase(_base_partran_seq_ptr, _adj_par_name_vec, _fixed_pars, _lower_bnd, _upper_bnd,
	//_obs_name_vec, file_manager_ptr), init_pars(_init_pars), obs_info_ptr(_obs_info_ptr), calc_obs_sen(_calc_obs_sen)
{
	log_trans_pars = _log_trans_pars;
}

void TornadoPlot::assemble_runs(RunManagerAbstract &run_manager)
{
	Parameters pars = fixed_ctl_pars;
	pars.insert(init_pars);

	// Assemble run Based on the initial Parameter values
	Parameters tmp_pars;
	tmp_pars = pars;
	base_partran_seq_ptr->ctl2model_ip(tmp_pars);
	int run_id = run_manager.add_run(tmp_pars, "base_run", Parameters::no_data);

	// Assemble runs which perturb each parameter to its max and min value
	for (const auto &ipar : adj_par_name_vec)
	{
		tmp_pars = pars;
		double low_value = lower_bnd.get_rec(ipar);
		tmp_pars[ipar] = low_value;
		base_partran_seq_ptr->ctl2model_ip(tmp_pars);
		run_id = run_manager.add_run(tmp_pars, ipar + " L", Parameters::no_data);

		tmp_pars = pars;
		double hi_value = upper_bnd.get_rec(ipar);
		tmp_pars[ipar] = hi_value;
		base_partran_seq_ptr->ctl2model_ip(tmp_pars);
		run_id = run_manager.add_run(tmp_pars, ipar + " U", Parameters::no_data);
	}
}
void  TornadoPlot::calc_sen(RunManagerAbstract &run_manager, ModelRun model_run)
{
	int n_runs = run_manager.get_nruns();
	bool run_ok = false;
	int run_status;
	string run_name;
	double par_value_not_used;
	cout << endl;
	if (n_runs <= 1)
	{
		cerr << "Can not perform tornado calculations: insufficient number of runs (" << n_runs << ")" << endl;
		return;
	}
	run_manager.get_info(0, run_status, run_name, par_value_not_used);
	if (run_status <= 0)
	{
		cerr << "Cannot perform tornado calculations: base run failed " << endl;
		return;
	}

	ofstream &fout_tor = file_manager_ptr->open_ofile_ext("tor");
	tornado_calc(run_manager, model_run, fout_tor, "");
	file_manager_ptr->close_file("tor");

	ofstream &fout_toi = file_manager_ptr->open_ofile_ext("toi");
	if (calc_obs_sen)
	{
		for (const string &iobs : obs_name_vec)
		{
			tornado_calc(run_manager, model_run, fout_toi, iobs);
		}
	}
	file_manager_ptr->close_file("toi");
}

void  TornadoPlot::tornado_calc(RunManagerAbstract &run_manager, ModelRun model_run, ofstream &fout, const string obs_name)
{
	// if obs_name == "" compute Tornado Polt for the global phi.
	// Otherwise, compute Tornada Plot for the named observation.
	ModelRun base_run = model_run;
	ModelRun run1 = model_run;
	Parameters pars_init;
	Parameters pars;
	Observations obs;

	const vector<string> &run_mngr_obs_name_vec = run_manager.get_obs_name_vec();

	int n_runs = run_manager.get_nruns();
	bool run_ok = false;

	stringstream message;
	cout << endl;
	if (n_runs <= 1)
	{
		cerr << "Can not perform tornado calculations: insufficient number of runs (" << n_runs << ")" << endl;
		return;
	}
	run_ok = run_manager.get_run(0, pars_init, obs);
	if (!run_ok)
	{
		cerr << "Cannot perform tornado calculations: base run failed " << endl;
		return;
	}
	base_partran_seq_ptr->model2ctl_ip(pars_init);
	model_run.update_ctl(pars_init, obs);
	double phi_base = 0;
	if (obs_name.empty())
	{
		 phi_base = model_run.get_phi(DynamicRegularization::get_zero_reg_instance());
	}
	else
	{
		phi_base = obs[obs_name];
	}

	int run_status;
	string par_name;
	string run_name;
	char run_type;
	double par_value_not_used;
	map<string, map<string, double> > phi_tornado_map;
	for (const auto &ipar : adj_par_name_vec)
	{
		phi_tornado_map[ipar] = map<string, double>();
	}

	for (int i_run = 1; i_run < n_runs; ++i_run)
	{
		std::cout << string(message.str().size(), '\b');
		message.str("");
		message << "processing run " << i_run + 1 << " / " << n_runs;
		std::cout << message.str();

		run_manager.get_info(i_run, run_status, run_name, par_value_not_used);
		strip_ip(run_name);
		run_type = run_name.back();
		par_name = run_name.substr(0, run_name.size() - 2);
		run_ok = run_manager.get_run(i_run, pars, obs);
		if (run_ok)
		{
			base_partran_seq_ptr->model2ctl_ip(pars);
			double phi;
			if (obs_name.empty())
			{
				model_run.update_ctl(pars, obs);
				phi= model_run.get_phi(DynamicRegularization::get_zero_reg_instance());
			}
			else
			{
				phi = obs[obs_name];
			}

			if (run_type == 'L')
			{
				phi_tornado_map[par_name]["lo"] = phi;
			}
			else if (run_type == 'U')
			{
				phi_tornado_map[par_name]["hi"] = phi;
			}
		}
	}
	//sort tornado plots
	vector< pair<string, double> > phi_range_vec;
	for (const auto &irec : phi_tornado_map)
	{
		double range = 0;
		map<string, double>::const_iterator map_end = irec.second.end();
		map<string, double>::const_iterator iter_lo = irec.second.find("lo");
		map<string, double>::const_iterator iter_hi = irec.second.find("hi");
		if (iter_lo != map_end && iter_hi != map_end)
		{
			range = max(iter_hi->second, max(phi_base, iter_lo->second)) - min(iter_hi->second, min(phi_base, iter_lo->second));
			range = abs(range);
		}
		else if (iter_hi != map_end)
		{
			range = iter_hi->second - phi_base;
			range = abs(range);
		}
		else if (iter_lo != map_end)
		{
			range = phi_base - iter_lo->second;
			range = abs(range);
		}
		phi_range_vec.push_back(make_pair(irec.first, range));
	}
	std::sort(phi_range_vec.begin(), phi_range_vec.end(), [](pair<string, double> a, pair<string, double> b) {
		return std::abs(a.second) > std::abs(b.second); });

	if (obs_name.empty())
	{
		fout << "Observation: Global Phi" << endl;
	}
	else
	{
		fout << "Observation: " << obs_name << endl;
	}
	fout << "parameter_name, phi_lo, phi_init, phi_hi, par_lo, par_init, par_hi" << endl;
	for (const auto &i : phi_range_vec)
	{
		const string &par_name = i.first;
		const auto irec = phi_tornado_map.find(par_name);
		if (irec != phi_tornado_map.end())
		{
			map<string, double>::const_iterator map_end = irec->second.end();
			map<string, double>::const_iterator iter_lo = irec->second.find("lo");
			map<string, double>::const_iterator iter_hi = irec->second.find("hi");

			stringstream phi_hi;
			if (iter_hi != map_end)
			{
				phi_hi << iter_hi->second;
			}
			else
			{
				phi_hi << "na";
			}
			stringstream phi_lo;
			if (iter_lo != map_end)
			{
				phi_lo << iter_lo->second;
			}
			else
			{
				phi_lo << "na";
			}
			fout << par_name << ", " << phi_lo.str() << ", " << phi_base << ", " << phi_hi.str()
				<< ", " << lower_bnd[par_name] << ", " << init_pars[par_name] << ", " << upper_bnd[par_name] << endl;
		}
	}
}



TornadoPlot::~TornadoPlot()
{
}

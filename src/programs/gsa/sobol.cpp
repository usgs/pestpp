#include <fstream>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <iterator>
#include "sobol.h"
#include "Transformable.h"
#include "RunManagerAbstract.h"
#include "ParamTransformSeq.h"
#include "ModelRunPP.h"
#include "Stats.h"
#include "FileManager.h"
#include "utilities.h"
#include "eigen_tools.h"

using namespace std;
using namespace Eigen;


Sobol::Sobol(Pest &_pest_scenario,
	FileManager &_file_manager, ObjectiveFunc *_obj_func_ptr,
	const ParamTransformSeq &_par_transform,
	int _n_sample, PARAM_DIST _par_dist, unsigned int _seed)
	: GsaAbstractBase(_pest_scenario, _file_manager, _obj_func_ptr, _par_transform,
		_par_dist, _seed), n_sample(_n_sample)
	{
	}

VectorXd Sobol::gen_rand_vec(long nsample, double min, double max)
{
	VectorXd v(nsample);
	long v_len = v.size();


	if (par_dist == PARAM_DIST::normal)
	{
		std::normal_distribution<> distribution((max + min) / 2.0, (max - min) / 4.0);
		for (long i = 0; i < v_len; ++i)
		{
			v[i] = distribution(rand_engine);
			while (v[i] < min || v[i] > max) v[i] = distribution(rand_engine);
		}
	}
	else
	{
		std::uniform_real_distribution<double> distribution(min, max);
		for (long i = 0; i < v_len; ++i)
		{
			v[i] = distribution(rand_engine);
		}
	}
	return v;
}

void Sobol::gen_m1_m2()
{
	long npar = adj_par_name_vec.size();
	//generate random matrices
	double par_min;
	double par_max;
	VectorXd v1;
	VectorXd v2;
	m1 = MatrixXd::Zero(n_sample, npar);
	m2 = MatrixXd::Zero(n_sample, npar);
	for (int i=0; i<npar; ++i)
	{
		string &p_name = adj_par_name_vec[i];
		par_min = min_numeric_pars[p_name];
		par_max = max_numeric_pars[p_name];
		v1 = gen_rand_vec(n_sample, par_min, par_max);
		v2 = gen_rand_vec(n_sample, par_min, par_max);
		m1.col(i) = v1;
		m2.col(i) = v2;
	}
}

MatrixXd Sobol::gen_N_matrix(const MatrixXd &m1, const MatrixXd &m2, const vector<int> &idx_vec)
{
  MatrixXd n = m1;
  for (int i : idx_vec)
  {
	n.col(i) = m2.col(i);
  }
  return n;
}

void Sobol::add_model_runs(RunManagerAbstract &run_manager, const MatrixXd &n, ofstream &f_out, string tag)
{
	int run_id;
    stringstream ss;
	for (int i=0; i<n_sample; ++i)
	{
		VectorXd tmp_vec =  n.row(i);
		Parameters tmp_pars(adj_par_name_vec, tmp_vec);
		base_partran_seq_ptr->numeric2model_ip(tmp_pars);
        ss.str("");
        ss << "tag:" << tag + "_isample:" << i;;
		run_id = run_manager.add_run(tmp_pars,ss.str());
		base_partran_seq_ptr->model2ctl_ip(tmp_pars);
		f_out << run_id << "," << ss.str();
		for (auto pname : adj_par_name_vec)
			f_out << "," << tmp_pars.get_rec(pname);
		f_out << endl;
	}
}

void Sobol::assemble_runs(RunManagerAbstract &run_manager)
{
	MatrixXd c;
	run_manager.reinitialize();
	gen_m1_m2();

	//open csv par file
	ofstream &f_out = file_manager_ptr->open_ofile_ext("sobol.par.csv");
	f_out << "run_id,info_txt";
	for (auto pname : adj_par_name_vec)
		f_out << "," << pest_utils::lower_cp(pname);
	f_out << endl;


	//calculate a0
	int n_adj_par = adj_par_name_vec.size();
	
	
	add_model_runs(run_manager, m1, f_out, "m1");
	add_model_runs(run_manager, m2, f_out,"m2");

	//calculate first order runs a1,....an
	vector<int> idx_vec;
    stringstream ss;
	for (int ai=0; ai<n_adj_par; ++ai)
	{
		idx_vec.clear();

		idx_vec.push_back(ai);
		c = gen_N_matrix(m1, m2, idx_vec);
		//cout << c << endl << endl;
        ss.str("");
        ss << "n-" << adj_par_name_vec[ai];
		add_model_runs(run_manager, c, f_out,ss.str());
	}
	f_out.close();
}


vector<double> Sobol::get_obs_vec(RunManagerAbstract &run_manager, int run_set, ModelRun &model_run, const string &obs_name)
{
	ModelRun run0 = model_run;
	int run_b = run_set * n_sample;
	int run_e = run_b + n_sample;
	Parameters pars0;
	Observations obs0;
	int nrun = 0;
	vector<double> obs_vec = vector<double>(n_sample, MISSING_DATA);
    int run_stat;
    string info_txt;
    double info_val;
	for (int run_id = run_b; run_id<run_e; ++run_id)
	{
		double obs = MISSING_DATA;
		
		//bool success = run_manager.get_run(run_id, pars0, obs0);
		//if (success)
		if (run_map.find(run_id) != run_map.end())
		{
			//obs0 = run_map[run_id];
			//run0.update_ctl(pars0, obs0);
            run_manager.get_info(run_id,run_stat,info_txt,info_val);
			obs = run_map[run_id].get_rec(obs_name);
			if (obs == Observations::no_data) obs = MISSING_DATA;
		}
		obs_vec[nrun] = obs;
		nrun++;
	}
	return obs_vec;
}


vector<double> Sobol::get_phi_vec(RunManagerAbstract &run_manager, int run_set, ModelRun &model_run)
{
	ModelRun run0 = model_run;

	int run_b = run_set * n_sample;
	int run_e = run_b + n_sample;

	Parameters pars0;
	Observations obs0;
	int nrun = 0;
	vector<double> phi_vec = vector<double>(n_sample, MISSING_DATA);
	DynamicRegularization zero_reg = DynamicRegularization::get_zero_reg_instance();
	for(int run_id=run_b; run_id<run_e; ++run_id)
	{
		double phi = MISSING_DATA;
		bool success = run_manager.get_run(run_id, pars0, obs0);
		if (success)
		{
			run0.update_ctl(pars0, obs0);
			phi = run0.get_phi(zero_reg);
		}
		phi_vec[nrun] = phi;
		nrun++;
	}
	return phi_vec;
}

void Sobol::process_runs(RunManagerAbstract& run_manager, ModelRun &model_run)
{
	run_map.clear();
	Parameters pars;
	Observations obs;
	ModelRun run = model_run;
	for (int run_id = 0; run_id < run_manager.get_nruns(); run_id++)
	{
		bool success = run_manager.get_run(run_id, pars, obs);
		if (success)
		{
			//run.update_ctl(pars, obs);
			run_map[run_id] = obs;
		}
	}
}

void Sobol::calc_sen(RunManagerAbstract &run_manager, ModelRun model_run)
{
	process_runs(run_manager, model_run);

	ofstream &fout_sbl = file_manager_ptr->open_ofile_ext("sbl");
	ofstream &f_out = file_manager_ptr->open_ofile_ext("sobol.obs.csv");
	ofstream& f_si = file_manager_ptr->open_ofile_ext("sobol.si.csv");
	ofstream& f_sti = file_manager_ptr->open_ofile_ext("sobol.sti.csv");

	f_si << "output";
	f_sti << "output";

	for (auto pname : pest_scenario_ptr->get_ctl_ordered_adj_par_names())
	{
		f_si << "," << pest_utils::lower_cp(pname);
		f_sti << "," << pest_utils::lower_cp(pname);
	}
	f_si << endl;
	f_sti << endl;

	//phi
	fout_sbl << "Sobol Sensitivity for PHI" << endl;
	pair<vector<double>, vector<double>> vals;
	vals = calc_sen_single(run_manager, model_run, fout_sbl, string());
	f_si << "phi";
	f_sti << "phi";
	for (auto v : vals.first)
		f_si << "," << v;
	f_si << endl;
	for (auto v : vals.second)
		f_sti << "," << v;
	f_sti << endl;


	vector<string> obs_names = run_manager.get_obs_name_vec();
	f_out << "run_id,failed_flag";
	for (auto oname : obs_names)
		f_out << "," << pest_utils::lower_cp(oname);
	f_out << endl;
	Observations obs0;
	Parameters pars0;
	for (int i=0; i < run_manager.get_nruns(); i++)
	{
		f_out << i;
		bool success = run_manager.get_run(i, pars0, obs0);
		if (success)
		{
			f_out << "," << 0;
			for (auto oname : obs_names)
				f_out << "," << obs0.get_rec(oname);
		}
		else
		{
			f_out << "," << 1;
			for (auto oname : obs_names)
				f_out << ",1e10";
		}
		f_out << endl;
	}
	file_manager_ptr->close_file("sobol.obs.csv");
	
	for (const string &iobs : obs_names)
	{
		fout_sbl << endl << endl;
		fout_sbl << "Sobol Sensitivity for observation \"" << iobs <<"\"" << endl;
		vals = calc_sen_single(run_manager, model_run, fout_sbl, iobs);
		string low_obs = pest_utils::lower_cp(iobs);
		f_si << low_obs;
		f_sti << low_obs;
		for (auto v : vals.first)
			f_si << "," << v;
		f_si << endl;
		for (auto v : vals.second)
			f_sti << "," << v;
		f_sti << endl;
	}

	file_manager_ptr->close_file("sbl");
	file_manager_ptr->close_file("sobol.si.csv");
	file_manager_ptr->close_file("sobol.sti.csv");

}


void Sobol::calc_sen_single_old(RunManagerAbstract& run_manager, ModelRun model_run, ofstream& fout_sbl, const string& obs_name)
{
	vector<double> ya;
	vector<double> yb;
	if (obs_name.empty())
	{
		ya = get_phi_vec(run_manager, 0, model_run);
		yb = get_phi_vec(run_manager, 1, model_run);
	}
	else
	{
		ya = get_obs_vec(run_manager, 0, model_run, obs_name);
		yb = get_obs_vec(run_manager, 1, model_run, obs_name);
	}

	vector<double> ya_yb_prod = vec_array_prod(ya, yb, MISSING_DATA);
	vector<double> y_ab;
	y_ab.reserve(ya.size() + yb.size()); // preallocate memory
	y_ab.insert(y_ab.end(), ya.begin(), ya.end());
	y_ab.insert(y_ab.end(), yb.begin(), yb.end());


	//Compute Mean for the S_i's
	double mean_sq_si = vec_mean_missing_data(ya_yb_prod, MISSING_DATA);
	// Compute Var for S_i's
	pair<double, size_t> data = sum_of_prod_missing_data(y_ab, y_ab, MISSING_DATA);
	//double var_si = data.first / (data.second - 2.0) - mean_sq_si;
	double var_si = (data.first / data.second) - mean_sq_si;

	//Compute Mean for the S_ti's
	double mean_sq_sti = pow(vec_mean_missing_data(y_ab, MISSING_DATA), 2.0);
	// Compute Var for S_ti's
	double sti_u = sobol_u_missing_data(yb, yb, MISSING_DATA);
	double var_sti = sti_u - mean_sq_sti;
	fout_sbl << "E(Y) = " << sqrt(mean_sq_si) << ";  Var(Y) = " << var_si << " (for S_i calculations)" << endl;
	fout_sbl << "E(Y) = " << sqrt(mean_sq_sti) << ";  Var(Y) = " << var_sti << " (for S_ti calculations)" << endl;
	size_t npar = adj_par_name_vec.size();

	fout_sbl << "parameter_name, s_i, st_i, n_runs" << endl;
	for (size_t i = 0; i < npar; ++i)
	{
		vector<double> yci;
		if (obs_name.empty())
		{
			yci = get_phi_vec(run_manager, i + 2, model_run);
		}
		else
		{
			yci = get_obs_vec(run_manager, i + 2, model_run, obs_name);
		}

		pair<double, int> sumprod_num = sum_of_prod_missing_data(ya, yci, MISSING_DATA);
		long int n_runs = sumprod_num.second;
		double sobol_uj = sumprod_num.first / (sumprod_num.second - 1.0);
		double si = (sobol_uj - mean_sq_si) / var_si;

		double sobol_umj = sobol_u_missing_data(yb, yci, MISSING_DATA);
		double sti = 1 - ((sobol_umj - mean_sq_sti) / var_sti);

		fout_sbl << adj_par_name_vec[i] << ", " << si << ", " << sti << ", " << n_runs << endl;
	}
}


pair<vector<double>, vector<double>> Sobol::calc_sen_single(RunManagerAbstract &run_manager, ModelRun model_run, ofstream &fout_sbl, const string &obs_name)
{
	vector<double> ya;
	vector<double> yb;
	if (obs_name.empty())
	{
		ya = get_phi_vec(run_manager, 0, model_run);
		yb = get_phi_vec(run_manager, 1, model_run);
	}
	else
	{
		ya = get_obs_vec(run_manager, 0, model_run, obs_name);
		yb = get_obs_vec(run_manager, 1, model_run, obs_name);
	}

	vector<double> y_ab;
	y_ab.reserve(ya.size() + yb.size()); // preallocate memory
	y_ab.insert(y_ab.end(), ya.begin(), ya.end());
	y_ab.insert(y_ab.end(), yb.begin(), yb.end());

	size_t npar = adj_par_name_vec.size();
	

	//Compute Mean for the S_i's
	double mean_sq_si = pow(vec_mean_missing_data(y_ab, MISSING_DATA),2.0);
	// Compute Var for S_i's
	pair<double, size_t> data = sum_of_prod_missing_data(y_ab, y_ab, MISSING_DATA);
	//double var_si = data.first / (data.second - 2.0) - mean_sq_si;
	double var_si = (data.first / data.second) - mean_sq_si;

	fout_sbl << "parameter_name, s_i, st_i, n_runs" << endl;
	vector<double> si_vals, sti_vals,r_yci, r_ya, r_yb;
	si_vals.reserve(npar);
	sti_vals.reserve(npar);
	double si, sti, numer;
	pair<double, int> sumprod;
	long int n_runs = 0;
	int r;
	for (size_t i=0; i<npar; ++i)
	{
		vector<double> yci;
		if (obs_name.empty())
		{
			yci = get_phi_vec(run_manager, i + 2, model_run);
		}
		else
		{
			yci = get_obs_vec(run_manager, i + 2, model_run, obs_name);
		}
		
		if (var_si == 0)
		{
			si = 0.0;
		}

		else
		{
			numer = si_saltelli_numer(ya, yb, yci, MISSING_DATA);
			si = numer / var_si;
		}

		
		if (var_si == 0.0)
			sti = 0.0;
		else
		{
			numer = sti_saltelli_numer(ya, yci, MISSING_DATA);
			sti = numer / var_si;
		}
		fout_sbl << adj_par_name_vec[i] << ", " << si << ", " << sti << ", " << n_runs << endl;
		si_vals.push_back(si);
		sti_vals.push_back(sti);
	}

	pair<vector<double>, vector<double>> vals(si_vals, sti_vals);
	return vals;

}

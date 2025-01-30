#include <vector>
#include <Eigen/Dense>
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
#include "MorrisMethod.h"
#include "Transformable.h"
#include "RunManagerAbstract.h"
#include "ParamTransformSeq.h"
#include "ModelRunPP.h"
#include "utilities.h"
#include "FileManager.h"
#include "Stats.h"
//#include "Ensemble.h"

using namespace std;
using namespace pest_utils;
using Eigen::MatrixXd;
using Eigen::VectorXd;

void MorrisObsSenFile::initialize(const vector<string> &_par_names_vec, const vector<string> &_obs_names_vec, double _no_data, const GsaAbstractBase *_gsa_abstract_base)
{
	par_names_vec =_par_names_vec;
	obs_names_vec = _obs_names_vec;
	gsa_abstract_base = _gsa_abstract_base;
	no_data = _no_data;
}

void MorrisObsSenFile::add_sen_run_pair(const std::string &par_name, double p1, Observations &obs1, double p2, Observations &obs2)
{
	assert (obs1.size() == obs2.size());

	// compute sensitivities of individual observations
	double isen;
	double del_par = p2 - p1;

	for (const auto &iobs : obs_names_vec)
	{
		isen = (obs2[iobs] - obs1[iobs]) / del_par;
		const auto id = make_pair(par_name, iobs);
		if (map_obs_stats.find(id) == map_obs_stats.end())
		{
			map_obs_stats[id] = RunningStats();
		}
		map_obs_stats[make_pair(par_name, iobs)].add(isen);
	}
}

void MorrisObsSenFile::calc_pooled_obs_sen(ofstream &fout_obs_sen, map<string, double> &obs_2_sen_weight,
	map<string, double> &par_2_sen_weight)
{
	fout_obs_sen << "par_name, n_samples, obs_name, mean, abs_mean, sigma, scaled_sen" << endl;
	for (const auto &ip : par_names_vec)
	{
		string ipar = ip;
		for (const auto &iobs : obs_names_vec)
		{
			double mean = no_data;
			double abs_mean = no_data;
			double sigma = no_data;
			int n_samples = 0;
			const auto id = make_pair(ipar, iobs);
			auto sen_itr = map_obs_stats.find(id);
			if (sen_itr != map_obs_stats.end() && sen_itr->second.comp_nsamples() > 0)
			{
				mean = sen_itr->second.comp_mean();
				abs_mean = sen_itr->second.comp_abs_mean();
				sigma = sen_itr->second.comp_sigma();
				n_samples = sen_itr->second.comp_nsamples();
			}
			string weighted_sen = "N/A";
			auto it_obs = obs_2_sen_weight.find(iobs);
			auto it_par = par_2_sen_weight.find(ipar);
			if (it_obs != obs_2_sen_weight.end() && it_par != par_2_sen_weight.end() && abs_mean != no_data)
			{
				stringstream sstr;
				double value = abs_mean * it_par->second / it_obs->second;
				sstr << value;
				weighted_sen = sstr.str();
			}
			fout_obs_sen << ipar << ", " << n_samples << ", " << iobs << ", " << mean << ", " << abs_mean << ", " << sigma << ", " << weighted_sen << endl;
		}
	}
}

void MorrisMethod::process_pooled_var_file()
{
    group_2_pool_group_map.clear();
	if (calc_obs_sen) {
        std::set<string> obs_group_names;
        for (const auto &imap : obs_info_ptr->groups)
            obs_group_names.insert(imap.first);
        string fname = file_manager_ptr->get_base_filename() + ".pgp";
        if (check_exist_in(fname)) {
            ifstream &fin = file_manager_ptr->open_ifile_ext("pgp");


            string line;
            string cur_pool_grp;
            regex reg_reg("regex\\s*\\(\"(.+)\"\\)", regex_constants::icase);
            regex reg_grp("pool_group\\s*\\((.+)\\)", regex_constants::icase);
            cmatch mr;
            while (getline(fin, line)) {
                if (regex_match(line.c_str(), mr, reg_reg)) {
                    regex inp_reg = regex(mr[1].str(), regex_constants::icase);

                    for (auto itr = obs_group_names.begin(); itr != obs_group_names.end();) {
                        auto here = itr++;
                        if (regex_match(*here, inp_reg)) {
                            group_2_pool_group_map[*here] = cur_pool_grp;
                            obs_group_names.erase(here);
                        }
                    }

                }
                if (regex_match(line.c_str(), mr, reg_grp)) {
                    cur_pool_grp = mr[1];
                }
            }
            file_manager_ptr->close_file("pgp");
        }
        else
        {
            for (auto& g : obs_group_names)
            {
                group_2_pool_group_map[g] = g;
            }


        }
    }
}

MatrixXd MorrisMethod::create_P_star_mat(int k)
{
	MatrixXd b_mat = create_B_mat(k);
	MatrixXd j_mat = create_J_mat(k);
	MatrixXd d_mat = create_D_mat(k);
	MatrixXd x_vec = create_x_vec(k);
	MatrixXd p_mat = create_P_mat(k);
	b_star_mat = j_mat.col(0) * x_vec.transpose() + (delta/2.0)*((2.0 * b_mat - j_mat) * d_mat + j_mat) ;
	return b_star_mat;
}

MorrisMethod::MorrisMethod(Pest &_pest_scenario,
	FileManager &_file_manager, ObjectiveFunc *_obj_func_ptr,
	const ParamTransformSeq &_par_transform,
	int _p, int _r, double _delta,
	bool _calc_pooled_obs, bool _calc_morris_obs_sen, PARAM_DIST _par_dist, unsigned int _seed)
	: GsaAbstractBase(_pest_scenario, _file_manager, _obj_func_ptr, _par_transform,
		_par_dist, _seed),
	calc_obs_sen(_calc_pooled_obs), calc_morris_obs_sen(_calc_morris_obs_sen)
{
    obs_info_ptr = pest_scenario_ptr->get_observation_info_ptr();
	rand_gen = mt19937(_pest_scenario.get_pestpp_options().get_random_seed());
	initialize(_p, _r, _delta);
}


void MorrisMethod::initialize(int _p, int _r, double _delta)
{
	p = _p;
	r = _r;
	delta = _delta;
	//delta = p / (2.0 * (p - 1));
}

MatrixXd MorrisMethod::create_B_mat(int k)
{
	MatrixXd b_mat = MatrixXd::Constant(k+1, k, 0.0);
	int nrow = k+1;
	int ncol = k;
	for (int icol=0; icol<ncol; ++icol)
	{
		for (int irow=icol+1; irow<nrow; ++irow)
		{
			b_mat(irow, icol) = 1.0;
		}
	}
	return b_mat;
}

MatrixXd MorrisMethod::create_J_mat(int k)
{
	MatrixXd j_mat = MatrixXd::Constant(k+1, k, 1.0);
	return j_mat;
}

MatrixXd MorrisMethod::create_D_mat(int k)
{
	MatrixXd d_mat=MatrixXd::Constant(k, k, 0.0);
	for (int i=0; i<k; ++i) {
		d_mat(i,i) = rand_plus_minus_1();
	}
	return d_mat;
}

MatrixXd MorrisMethod::create_P_mat(int k)
{
	// generate and return a random permutation matrix
	MatrixXd p_mat=MatrixXd::Constant(k, k, 0.0);
	vector<int> rand_idx;
	rand_idx.reserve(k);
	for (int i=0; i<k; ++i) {
		rand_idx.push_back(i);
	}
	// Shuffle random index vector
	shuffle(rand_idx.begin(), rand_idx.end(), rand_engine);

	for(int irow=0; irow<k; ++irow)
	{
	  p_mat(irow, rand_idx[irow]) = 1;
	}

	return p_mat;
}

VectorXd MorrisMethod::create_x_vec(int k)
{
	// Warning Satelli's book is wrong.  Need to use Morris's original paper
	// to compute x.  Maximim value of any xi shoud be 1.0 - delta.
	VectorXd x(k);
	double rnum;
	// Used with the randon number generator below, max_num will produce a
	// random interger which varies between 0 and (p-2)/2.  This will in turn
	// produce the correct xi's which vary from {0, 1/(p-1), 2/(p-1), .... 1 - delta)
	// when divided by (p-1).  See Morris's paper for the derivaition.
	// Satelli's book has this equation wrong!
	int max_num = 1+(p-2)/2;
	for (int i=0; i<k; ++i)
	{
		rnum = rand_gen() % max_num;
		x(i) = rnum / (p - 1);
	}
	return x;
}


int MorrisMethod::rand_plus_minus_1(void)
{
	//return (rand_engine() % 2 * 2) - 1;
	return (rand_gen() % 2 * 2) - 1;
}


Parameters MorrisMethod::get_numeric_parameters(int row)
{
	Parameters numeric_pars;
	auto e=adj_par_name_vec.end();
	size_t n_cols = b_star_mat.cols();
	for (int j=0; j<n_cols; ++j)
	{
		const string &p = adj_par_name_vec[j];
		auto it_lbnd = min_numeric_pars.find(p);
		assert(it_lbnd != min_numeric_pars.end());
		auto it_ubnd = max_numeric_pars.find(p);
		assert(it_ubnd != max_numeric_pars.end());
		numeric_pars[p] =  it_lbnd->second + (it_ubnd->second - it_lbnd->second) * b_star_mat(row, j);
	}
	return numeric_pars;
}

void MorrisMethod::assemble_runs(RunManagerAbstract &run_manager)
{
	map<int, Parameters> par_map;
	for (int tmp_r=0; tmp_r<r; ++tmp_r)
	{
		b_star_mat = create_P_star_mat(adj_par_name_vec.size());
		auto n_rows = b_star_mat.rows();
		int run_id;
		string par_name = "";
		for (int i=0; i<n_rows; ++i)
		{
			par_name.clear();
			//get control parameters
			Parameters pars = get_numeric_parameters(i);
			// convert control parameters to model parameters
			base_partran_seq_ptr->numeric2model_ip(pars);
			// convert control parameters to model parameters
			//base_partran_seq_ptr->ctl2model_ip(pars);
			if (i>0)
			{
				par_name = adj_par_name_vec[i-1];
			}
			//note: the "par_name:" prefix is hard coded later when processing runs...
			run_id = run_manager.add_run(pars, "par_name:"+par_name, Parameters::no_data);
			base_partran_seq_ptr->model2ctl_ip(pars);
			par_map[run_id] =  pars;

		}
	}
	//write the parameter sequence to a csv file
	string filename = file_manager_ptr->get_base_filename() + ".sen.par.csv";
	ofstream fout(filename);
	vector<string> var_names = pest_scenario_ptr->get_ctl_ordered_par_names();
	fout << "run_id";
	for (auto pname : var_names)
		fout << "," << pest_utils::lower_cp(pname);
	fout << endl;
	for (auto &p : par_map)
	{
		fout << p.first;
		for (auto &pname : var_names)
			fout << "," << p.second.get_rec(pname);
		fout << endl;
	}
	/*ParameterEnsemble pe(pest_scenario_ptr);
	vector<string> real_names,var_names=pest_scenario_ptr->get_ctl_ordered_par_names();
	stringstream ss;
	map<int, string> real_name_map;
	for (auto &p : par_map)
	{
		ss.str("");
		ss << p.first;
		real_names.push_back(ss.str());
		real_name_map[p.first] = ss.str();
	}
	pe.set_trans_status(ParameterEnsemble::transStatus::CTL);
	pe.reserve(real_names, var_names);
	for (auto &p : par_map)
	{
		pe.update_real_ip(real_name_map[p.first], p.second.get_data_eigen_vec(var_names));
	}
	string filename = file_manager_ptr->get_base_filename() + ".sen.par.csv";
	pe.to_csv(filename);*/


}

void  MorrisMethod::calc_sen(RunManagerAbstract &run_manager, ModelRun model_run)
{
	ofstream &fout_morris = file_manager_ptr->open_ofile_ext("msn");
	ofstream &fout_raw = file_manager_ptr->open_ofile_ext("raw.csv");
    ofstream &fout_grp = file_manager_ptr->open_ofile_ext("group.raw.csv");
    ofstream &fout_grpmorris = file_manager_ptr->open_ofile_ext("group.msn");

    ModelRun run0 = model_run;
	ModelRun run1 = model_run;
	Parameters pars0;
	Observations obs0;
	Parameters pars1;
	Observations obs1;

	map<string, RunningStats > sen_map;
	map<string, RunningStats> obs_stats_map;
    map<pair<string,string>, RunningStats> ogroup_sens_map;

	for (auto &it_p : adj_par_name_vec)
	{
		sen_map[it_p] = RunningStats();
        for (auto& gname : pest_scenario_ptr->get_ctl_ordered_obs_group_names())
        {
            ogroup_sens_map[make_pair(it_p,gname)] = RunningStats();
        }
	}

	for (auto &it_obs : obs_name_vec)
	{
		obs_stats_map[it_obs] = RunningStats();
	}

	const vector<string> &run_mngr_obs_name_vec = run_manager.get_obs_name_vec();
	obs_sen_file.initialize(adj_par_name_vec, run_mngr_obs_name_vec, Observations::no_data, this);

	fout_raw << "parameter_name,phi_0,phi_1,par_0,par_1,elem_effect" << endl;
    fout_grp << "parameter_name,obs_group,phi_0,phi_1,par_0,par_1,elem_effect" << endl;

    int n_runs = run_manager.get_nruns();
	bool run0_ok = false;
	bool run1_ok = false;
	string par_name_1;
	double null_value;
	stringstream message;
	cout << endl;
	run1_ok = run_manager.get_run(0, pars1, obs1);
	base_partran_seq_ptr->model2numeric_ip(pars1);
	for (int i_run=1; i_run<n_runs; ++i_run)
	{


		run0_ok = run1_ok;
		pars0 = pars1;
		obs0 = obs1;
		run1_ok = run_manager.get_run(i_run, pars1, obs1, par_name_1, null_value);
		if (par_name_1.empty())
        {
		    throw runtime_error("Morris: empty par name");
        }
		//remove the "par_name:" prefix
        //cout << par_name_1 << endl;
		par_name_1 = par_name_1.substr(9);
		//cout << par_name_1 << endl << endl;
		if ((!par_name_1.empty()) && (pars1.find(par_name_1) == pars1.end()))
        {
		    throw runtime_error("Morris:parameter name not found: "+par_name_1);
        }
		base_partran_seq_ptr->model2numeric_ip(pars1);
		// Add run0 to obs_stats
		if (run0_ok)
		{
			for (const auto &i_obs : run_mngr_obs_name_vec)
			{
				auto it = obs0.find(i_obs);
				if (it != obs0.end() && it->second != Observations::no_data)
				{
					obs_stats_map[i_obs].add(obs0[i_obs]);
				}
			}
		}

		if (run0_ok && run1_ok && !par_name_1.empty())
		{
			Parameters tmp_ctl_par = base_partran_seq_ptr->numeric2ctl_cp(pars0);
			run0.update_ctl(tmp_ctl_par, obs0);
			double phi0 = run0.get_phi(DynamicRegularization::get_zero_reg_instance());
			map<string,double> pg0 = run0.get_obj_func_ptr()->get_group_phi(obs0,tmp_ctl_par,DynamicRegularization::get_zero_reg_instance());
			tmp_ctl_par = base_partran_seq_ptr->numeric2ctl_cp(pars1);
			run1.update_ctl(tmp_ctl_par, obs1);
			double phi1 = run1.get_phi(DynamicRegularization::get_zero_reg_instance());
            map<string,double> pg1 = run1.get_obj_func_ptr()->get_group_phi(obs1,tmp_ctl_par,DynamicRegularization::get_zero_reg_instance());

            double p0 = pars0[par_name_1];
			double p1 = pars1[par_name_1];
			// compute standard Morris Sensitivity on the global objective function
			double sen = (phi1 - phi0) / delta;
			fout_raw << pest_utils::lower_cp(par_name_1) << "," << phi1 << "," << phi0 << "," << p1 << "," << p0 << "," << sen << endl;

			const auto &it_senmap = sen_map.find(par_name_1);
			if (it_senmap != sen_map.end())
			{
				it_senmap->second.add(sen);
			}

			//Compute sensitvities of indiviual observations
			obs_sen_file.add_sen_run_pair(par_name_1, p0, obs0, p1, obs1);
            for (auto& gphi : pg0)

            {
                phi0 = gphi.second;
                phi1 = pg1.at(gphi.first);
                sen = (phi1 - phi0) / delta;
                fout_grp << pest_utils::lower_cp(par_name_1) << "," << gphi.first << "," << phi1 << "," << phi0 << "," << p1 << "," << p0 << "," << sen << endl;
                ogroup_sens_map[make_pair(par_name_1,gphi.first)].add(sen);
            }
        }
	}
	// Add final run to obs_stats
	if (run1_ok)
	{
		for (const auto &i_obs : run_mngr_obs_name_vec)
		{
			auto it = obs1.find(i_obs);
			if (it != obs1.end() && it->second != Observations::no_data)
			{
				obs_stats_map[i_obs].add(obs0[i_obs]);
			}
		}
	}
	cout << endl;
	cout << "writing output files" << endl;
	// write standard Morris Sensitivity for the global objective function
	fout_morris << "parameter_name,n_samples,sen_mean,sen_mean_abs,sen_std_dev" << endl;
    fout_grpmorris << "parameter_name,obs_group_name,n_samples,sen_mean,sen_mean_abs,sen_std_dev" << endl;

    for (const auto &it_par : adj_par_name_vec)
	{
		const auto &it_senmap = sen_map.find(it_par);
		if (it_senmap != sen_map.end())
		{
			fout_morris << pest_utils::lower_cp(it_par) << "," << it_senmap->second.comp_nsamples() << "," << it_senmap->second.comp_mean() << "," << it_senmap->second.comp_abs_mean() << "," << sqrt(it_senmap->second.comp_var()) << endl;
		}
	}

    for (auto& pname : adj_par_name_vec)
    {
        for (auto& ogname : pest_scenario_ptr->get_ctl_ordered_obs_group_names()) {

            auto p = make_pair(pname,ogname);
            if (ogroup_sens_map.find(p) == ogroup_sens_map.end())
            {
                cout << "missing " << pname << "," << ogname << endl;
                continue;
            }
            //fout_morris << pest_utils::lower_cp(it_par) << "," << it_senmap->second.comp_nsamples() << ","
            // << it_senmap->second.comp_mean() << "," << it_senmap->second.comp_abs_mean() << ","
            // << sqrt(it_senmap->second.comp_var()) << endl;
            fout_grpmorris << pest_utils::lower_cp(pname) << "," << pest_utils::lower_cp(ogname);
            fout_grpmorris << "," << ogroup_sens_map.at(p).comp_nsamples() << "," << ogroup_sens_map.at(p).comp_mean();
            fout_grpmorris << "," << ogroup_sens_map.at(p).comp_abs_mean() << "," << sqrt(ogroup_sens_map.at(p).comp_var()) << endl;
         }
    }

	if (calc_morris_obs_sen)
	{
		// write standard Morris Sensitivity for individual observations
		ofstream &fout_mis = file_manager_ptr->open_ofile_ext("mio");
		calc_morris_obs(fout_mis, obs_sen_file);
		file_manager_ptr->close_file("mio");
	}

	if (calc_obs_sen)
	{
		ofstream &fout_mos = file_manager_ptr->open_ofile_ext("mos");
		////compute pooled standard deviations
		map<string, double> obs_2_sen_weight;
		{
			//Compute Pooled Standard Deviations
			map<string, vector<RunningStats> > tmp_pool_grps;
			for (const auto &it : obs_stats_map)
			{
				const string &obs_name = it.first;
				const string &obs_group = obs_info_ptr->get_group(obs_name);
				auto it_pg = group_2_pool_group_map.find(obs_group);
				if (it_pg != group_2_pool_group_map.end())
				{
					const string &pool_group = it_pg->second;
					if (tmp_pool_grps.find(pool_group) == tmp_pool_grps.end())
					{
						tmp_pool_grps[pool_group] = vector<RunningStats>();
					}
					tmp_pool_grps[pool_group].push_back(it.second);
				}
			}
			////compute pooled standard deviations
			map<string, double> obs_pooled_grp_std_dev;
			for (const auto &i_pgrp : tmp_pool_grps)
			{
				const string &pool_group = i_pgrp.first;
				double var_sum = 0;
				long int weight_sum = 0;
				for (const auto &istat : i_pgrp.second)
				{
					long int weight = istat.comp_nsamples() - 1;
					var_sum += weight * istat.comp_var();
					weight_sum += weight;
				}
				if (weight_sum > 0)
				{
					obs_pooled_grp_std_dev[pool_group] = sqrt(var_sum / weight_sum);
				}
			}
			for (const auto &iobs : obs_name_vec)
			{
				const string &obs_name = iobs;
				const string &obs_group = obs_info_ptr->get_group(obs_name);
				auto it_pg = group_2_pool_group_map.find(obs_group);
				if (it_pg != group_2_pool_group_map.end())
				{
					const string &pool_group = it_pg->second;
					if (obs_pooled_grp_std_dev.find(pool_group) != obs_pooled_grp_std_dev.end())
					{
						obs_2_sen_weight[iobs] = obs_pooled_grp_std_dev[pool_group];
					}
				}
			}
		}

		//compute parameter standard deviations
		map<string, double> par_std_dev;
		par_std_dev = calc_parameter_unif_std_dev();
		obs_sen_file.calc_pooled_obs_sen(fout_mos, obs_2_sen_weight, par_std_dev);
		file_manager_ptr->close_file("mos");
	}
	file_manager_ptr->close_file("msn");
    file_manager_ptr->close_file("group.msn");

    file_manager_ptr->close_file("group.raw.csv");
    file_manager_ptr->close_file("raw.csv");
}

void MorrisMethod::calc_morris_obs(ostream &fout, MorrisObsSenFile &morris_sen_file)
{
	// write standard Morris Sensitivity
	fout << "observation_name,parameter_name,n_samples,sen_mean,sen_mean_abs,sen_std_dev" << endl;

	for (const auto &i_obs : morris_sen_file.obs_names_vec)
	{
		//fout << "Method of Morris for observation: " << i_obs << endl;
		//fout << "parameter_name, n_samples, sen_mean, sen_mean_abs, sen_std_dev" << endl;
		//fout << i_obs;
		for (const auto &i_par : morris_sen_file.par_names_vec)
		{
			const auto &it_senmap = morris_sen_file.map_obs_stats.find(make_pair(i_par, i_obs));
			if (it_senmap != morris_sen_file.map_obs_stats.end())
			{
				fout << pest_utils::lower_cp(i_obs) <<"," << pest_utils::lower_cp(i_par) << "," << it_senmap->second.comp_nsamples() << "," << it_senmap->second.comp_mean() << "," << it_senmap->second.comp_abs_mean() << "," << sqrt(it_senmap->second.comp_var()) << endl;
			}
		}
		//fout << endl;
	}
}


MorrisMethod::~MorrisMethod(void)
{
}

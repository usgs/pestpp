#include <random>
#include <map>
#include <iomanip>
#include <mutex>
#include <thread>
#include "Ensemble.h"
#include "RestartController.h"
#include "utilities.h"
#include "Ensemble.h"
#include "EnsembleSmoother.h"
#include "ObjectiveFunc.h"
#include "covariance.h"
#include "RedSVD-h.h"
#include "SVDPackage.h"
#include "eigen_tools.h"


L2PhiHandler::L2PhiHandler(Pest *_pest_scenario, FileManager *_file_manager,
	ObservationEnsemble *_oe_base, ParameterEnsemble *_pe_base,
	Covariance *_parcov, bool should_prep_csv)
{
	pest_scenario = _pest_scenario;
	file_manager = _file_manager;
	oe_base = _oe_base;
	pe_base = _pe_base;
	
	//check for inequality constraints
	//for (auto &og : pest_scenario.get_ctl_ordered_obs_group_names())
	string og;
	double weight;
	const ObservationInfo* oi = pest_scenario->get_ctl_observation_info_ptr();
	for (auto &oname : pest_scenario->get_ctl_ordered_obs_names())
	{
		og = oi->get_group(oname);
		weight = oi->get_weight(oname);
		if (weight == 0)
			continue;
		if ((og.compare(0, 2, "L_") == 0) || (og.compare(0, 4, "LESS")==0))
		{
			lt_obs_names.push_back(oname);
		}
		else if ((og.compare(0, 2, "G_")==0) || (og.compare(0, 7, "GREATER")==0))
		{
			gt_obs_names.push_back(oname);
		}
	}

	//save the org reg factor and org q vector
	org_reg_factor = pest_scenario->get_pestpp_options().get_ies_reg_factor();
	org_q_vec = get_q_vector();
	//Eigen::VectorXd parcov_inv_diag = parcov_inv.e_ptr()->diagonal();
	parcov_inv_diag = _parcov->e_ptr()->diagonal();
	for (int i = 0; i < parcov_inv_diag.size(); i++)
		parcov_inv_diag(i) = 1.0 / parcov_inv_diag(i);

	//parcov_inv = _parcov->inv();
	//parcov.inv_ip();
	oreal_names = oe_base->get_real_names();
	preal_names = pe_base->get_real_names();
	if (should_prep_csv)
	{
		prepare_csv(file_manager->open_ofile_ext("phi.actual.csv"), oreal_names);
		prepare_csv(file_manager->open_ofile_ext("phi.meas.csv"), oreal_names);
		prepare_csv(file_manager->open_ofile_ext("phi.composite.csv"), oreal_names);
		prepare_csv(file_manager->open_ofile_ext("phi.regul.csv"), preal_names);
		prepare_group_csv(file_manager->open_ofile_ext("phi.group.csv"));
	}

}

Eigen::MatrixXd L2PhiHandler::get_obs_resid(ObservationEnsemble &oe, bool apply_ineq)
{
	vector<string> names = oe_base->get_var_names();
	Eigen::MatrixXd resid = oe.get_eigen(vector<string>(),names) -
		oe_base->get_eigen(oe.get_real_names(), vector<string>());
	
	if (apply_ineq)
		apply_ineq_constraints(resid,names);
	return resid;
}


Eigen::MatrixXd L2PhiHandler::get_obs_resid_subset(ObservationEnsemble &oe, bool apply_ineq)
{
	vector<string> names = oe.get_var_names();
	Eigen::MatrixXd resid = oe.get_eigen() - oe_base->get_eigen(oe.get_real_names(), names);
	if (apply_ineq)
		apply_ineq_constraints(resid, names);
	return resid;
}

Eigen::MatrixXd L2PhiHandler::get_par_resid(ParameterEnsemble &pe)
{
	Eigen::MatrixXd resid = pe.get_eigen(vector<string>(), pe_base->get_var_names()) -
		pe_base->get_eigen(pe.get_real_names(), vector<string>());
	return resid;
}

Eigen::MatrixXd L2PhiHandler::get_par_resid_subset(ParameterEnsemble &pe)
{
	Eigen::MatrixXd resid = pe.get_eigen() - pe_base->get_eigen(pe.get_real_names(),pe.get_var_names());
	return resid;
}

Eigen::MatrixXd L2PhiHandler::get_actual_obs_resid(ObservationEnsemble &oe)
{
	//vector<string> act_obs_names = pest_scenario->get_ctl_ordered_nz_obs_names();
	vector<string> act_obs_names = oe_base->get_var_names();
	Eigen::MatrixXd resid(oe.shape().first, act_obs_names.size());
	resid.setZero();
	Observations obs = pest_scenario->get_ctl_observations();
	Eigen::MatrixXd oe_vals = oe.get_eigen(vector<string>(), act_obs_names);
	Eigen::MatrixXd ovals = obs.get_data_eigen_vec(act_obs_names);
	ovals.transposeInPlace();
	for (int i = 0; i < resid.rows(); i++)
		resid.row(i) = oe_vals.row(i) - ovals;
	apply_ineq_constraints(resid, act_obs_names);
	return resid;
}

Eigen::VectorXd L2PhiHandler::get_q_vector()
{
	const ObservationInfo* oinfo = pest_scenario->get_ctl_observation_info_ptr();
	Eigen::VectorXd q;
	//vector<string> act_obs_names = pest_scenario->get_ctl_ordered_nz_obs_names();
	vector<string> act_obs_names = oe_base->get_var_names();

	/*if (act_obs_names.size() == 0)
		act_obs_names = oe_base->get_var_names();*/
	q.resize(act_obs_names.size());
	double w;
	for (int i = 0; i < act_obs_names.size(); i++)
	{
		q(i) = oinfo->get_weight(act_obs_names[i]);
	}
	return q;
}

map<string, double> L2PhiHandler::get_obs_group_contrib(Eigen::VectorXd &phi_vec)
{
	map<string, double> group_phi_map;
	double sum;
	for (auto &og : obs_group_idx_map)
	{
		sum = 0.0;
		for (auto i : og.second)
			sum += phi_vec[i];
		group_phi_map[og.first] = sum;
	}

	return group_phi_map;
}

map<string, double> L2PhiHandler::get_par_group_contrib(Eigen::VectorXd &phi_vec)
{
	map<string, double> group_phi_map;
	double sum;
	for (auto &pg : par_group_idx_map)
	{
		sum = 0.0;
		for (auto i : pg.second)
			sum += phi_vec[i];
		group_phi_map[pg.first] = sum;
	}
	return group_phi_map;
}

void L2PhiHandler::update(ObservationEnsemble & oe, ParameterEnsemble & pe)
{
	//build up obs group and par group idx maps for group reporting
	obs_group_idx_map.clear();
	vector<string> nnz_obs = oe_base->get_var_names();
	ObservationInfo oinfo = pest_scenario->get_ctl_observation_info();
	vector<int> idx;

	/*for (auto& og : pest_scenario->get_ctl_ordered_obs_group_names())
	{
		idx.clear();
		for (int i = 0; i < nnz_obs.size(); i++)
		{
			if (oinfo.get_group(nnz_obs[i]) == og)
			{
				idx.push_back(i);
			}
		}
		if (idx.size() > 0)
			obs_group_idx_map[og] = idx;
	}*/
	for (auto& og : pest_scenario->get_ctl_ordered_obs_group_names())
		obs_group_idx_map[og] = vector<int>();

	for (int i = 0; i < nnz_obs.size(); i++)
	{
		obs_group_idx_map[oinfo.get_group(nnz_obs[i])].push_back(i);
	}
	
	
	//update the various phi component vectors
	meas.clear();
	obs_group_phi_map.clear();
	Eigen::VectorXd q = get_q_vector();
	map<string, Eigen::VectorXd> meas_map = calc_meas(oe, q);
	for (auto &pv : meas_map)
	{
		meas[pv.first] = pv.second.sum();

	}
	if (org_reg_factor != 0.0)
	{
		par_group_idx_map.clear();
		vector<string> pars = pe_base->get_var_names();
		ParameterInfo pi = pest_scenario->get_ctl_parameter_info();
		for (auto& pg : pest_scenario->get_ctl_ordered_par_group_names())
		{
			idx.clear();
			for (int i = 0; i < pars.size(); i++)
			{
				if (pi.get_parameter_rec_ptr(pars[i])->group == pg)
					idx.push_back(i);
			}
			if (idx.size() > 0)
				par_group_idx_map[pg] = idx;
		}
		regul.clear();
		map<string, Eigen::VectorXd> reg_map = calc_regul(pe);//, *reg_factor);
		//for (auto &pv : calc_regul(pe))
		string name;
		//big assumption - if oe is a diff shape, then this
		//must be a subset, so just use the first X rows of pe
		for (int i = 0; i < oe.shape().first; i++)
		{
			//name = preal_names[i];
			name = pe.get_real_names()[i];
			//cout << name << endl;
			regul[name] = reg_map[name].sum();
			par_group_phi_map[name] = get_par_group_contrib(reg_map[name]);
		}
	}
	
	actual.clear();
	for (auto &pv : calc_actual(oe, q))
	{
		actual[pv.first] = pv.second.sum();
		obs_group_phi_map[pv.first] = get_obs_group_contrib(pv.second);
	}
 	composite.clear();
	composite = calc_composite(meas, regul);
}

void L2PhiHandler::save_residual_cov(ObservationEnsemble& oe, int iter)
{
	Eigen::MatrixXd rmat = get_obs_resid(oe, false); //dont apply ineq constraints
	ObservationEnsemble res(pest_scenario, oe.get_rand_gen_ptr());
	res.reserve(oe.get_real_names(), oe_base->get_var_names());
	res.set_eigen(rmat);
	pair<Covariance,Covariance> rcovs = res.get_empirical_cov_matrices(file_manager);
	stringstream ss;
	ss << file_manager->get_base_filename() << "." << iter << ".res.";
	if (pest_scenario->get_pestpp_options().get_ies_save_binary())
	{
		ss << "jcb";
		rcovs.first.to_binary_new(ss.str());
	}
	else
	{
		ss << "cov";
		rcovs.first.to_ascii(ss.str());
	}

	
	ss.str("");
	ss << file_manager->get_base_filename() << "." << iter << ".shrunk_res.";
	if (pest_scenario->get_pestpp_options().get_ies_save_binary())
	{
		ss << "jcb";
		rcovs.second.to_binary_new(ss.str());
	}
	else
	{
		ss << "cov";
		rcovs.second.to_ascii(ss.str());
	}

}

map<string, double>* L2PhiHandler::get_phi_map(L2PhiHandler::phiType &pt)
{
	switch (pt)
	{
	case L2PhiHandler::phiType::ACTUAL:
		return &actual;
	case L2PhiHandler::phiType::COMPOSITE:
		return &composite;
	case L2PhiHandler::phiType::MEAS:
		return &meas;
	case L2PhiHandler::phiType::REGUL:
		return &regul;
	}
	throw runtime_error("PhiHandler::get_phi_map() didn't find a phi map...");
}

double L2PhiHandler::calc_mean(map<string, double> *phi_map)
{
	double mean = 0.0;
	map<string, double>::iterator pi = phi_map->begin(), end = phi_map->end();
	for (; pi != end; ++pi)
		mean = mean + pi->second;
	return mean / phi_map->size();
}

double L2PhiHandler::calc_std(map<string, double> *phi_map)
{
	double mean = calc_mean(phi_map);
	double var = 0.0;
	map<string, double>::iterator pi = phi_map->begin(), end = phi_map->end();
	for (; pi != end; ++pi)
		var = var + (pow(pi->second - mean, 2));
	if (var == 0.0)
		return 0.0;
	return sqrt(var / (phi_map->size() - 1));
}

double L2PhiHandler::get_mean(phiType pt)
{
	//double mean = 0.0;
	map<string, double>* phi_map = get_phi_map(pt);
	/*map<string, double>::iterator pi = phi_map->begin(), end = phi_map->end();
	for (;pi != end; ++pi)
		mean = mean + pi->second;
	return mean / phi_map->size();*/
	return calc_mean(phi_map);
}

double L2PhiHandler::get_std(phiType pt)
{
	//double mean = get_mean(pt);
	//double var = 0.0;
	map<string, double>* phi_map = get_phi_map(pt);
	/*map<string, double>::iterator pi = phi_map->begin(), end = phi_map->end();
	for (; pi != end; ++pi)
		var = var + (pow(pi->second - mean,2));
	if (var == 0.0)
		return 0.0;
	return sqrt(var/(phi_map->size()-1));*/
	return calc_std(phi_map);
}

double L2PhiHandler::get_max(phiType pt)
{
	double mx = -1.0e+30;
	map<string, double>* phi_map = get_phi_map(pt);
	map<string, double>::iterator pi = phi_map->begin(), end = phi_map->end();
	for (; pi != end; ++pi)
		mx = (pi->second > mx) ? pi->second : mx;
	return mx;
}

double L2PhiHandler::get_min(phiType pt)
{
	double mn = 1.0e+30;
	map<string, double>* phi_map = get_phi_map(pt);
	map<string, double>::iterator pi = phi_map->begin(), end = phi_map->end();
	for (; pi != end; ++pi)
		mn = (pi->second < mn) ? pi->second : mn;
	return mn;
}

map<string, double> L2PhiHandler::get_summary_stats(L2PhiHandler::phiType pt)
{
	map<string, double> stats;
	stats["mean"] = get_mean(pt);
	stats["std"] = get_std(pt);
	stats["max"] = get_max(pt);
	stats["min"] = get_min(pt);
	return stats;
}

string L2PhiHandler::get_summary_string(L2PhiHandler::phiType pt)
{
	map<string, double> stats = get_summary_stats(pt);
	stringstream ss;
	string typ;
	switch (pt)
	{
	case L2PhiHandler::phiType::ACTUAL:
		typ = "actual";
		break;
	case L2PhiHandler::phiType::MEAS:
		typ = "measured";
		break;
	case L2PhiHandler::phiType::REGUL:
		typ = "regularization";
		break;
	case L2PhiHandler::phiType::COMPOSITE:
		typ = "composite";
		break;
	}
	ss << setw(15) << typ << setw(15) << stats["mean"] << setw(15) << stats["std"] << setw(15) << stats["min"] << setw(15) << stats["max"] << endl;
	return ss.str();
}

string L2PhiHandler::get_summary_header()
{
	stringstream ss;
	ss << setw(15) << "phi type" << setw(15) << "mean" << setw(15) << "std" << setw(15) << "min" << setw(15) << "max" << endl;
	return ss.str();
}


void L2PhiHandler::report(bool echo)
{
	ofstream& f = file_manager->rec_ofstream();
	string s;
	f << get_summary_header();
	if (echo)
		cout << get_summary_header();
	if (pest_scenario->get_pestpp_options().get_ies_no_noise())
	{
		if (org_reg_factor == 0)
		{
			s = get_summary_string(L2PhiHandler::phiType::ACTUAL);
			f << s;
			if (echo)
				cout << s;
		}
		else
		{
			
			string s = get_summary_string(L2PhiHandler::phiType::COMPOSITE);
			f << s;
			if (echo)
				cout << s;
			s = get_summary_string(L2PhiHandler::phiType::REGUL);
			f << s;
			if (echo)
				cout << s;
			s = get_summary_string(L2PhiHandler::phiType::ACTUAL);
			f << s;
			if (echo)
				cout << s;
		}
		
	}
	else
	{
		s = get_summary_string(L2PhiHandler::phiType::MEAS);
		f << s;
		if (echo)
			cout << s;
		if (org_reg_factor == 0.0)
		{
			s = get_summary_string(L2PhiHandler::phiType::ACTUAL);
			f << s;
			if (echo)
				cout << s;	
		}
		else
		{
			string s = get_summary_string(L2PhiHandler::phiType::COMPOSITE);
			f << s;
			if (echo)
				cout << s;
			s = get_summary_string(L2PhiHandler::phiType::REGUL);
			f << s;
			if (echo)
				cout << s;
			s = get_summary_string(L2PhiHandler::phiType::ACTUAL);
			f << s;
			if (echo)
				cout << s;
		}

	}
		
	if (org_reg_factor != 0.0)
	{

		f << "     note: 'regularization' phi reported above does not " << endl;
		f << "           include the effects of reg_factor, " << endl;
		f << "           but 'composite' phi does." << endl;
		if (echo)
		{
			cout << "     note: 'regularization' phi reported above does not " << endl;
			cout << "           include the effects of reg_factor, " << endl;
			cout << "           but 'composite' phi does." << endl;
		}
	}
	if (!pest_scenario->get_pestpp_options().get_ies_no_noise())
	{
		f << "     note: 'measured' phi reported above includes " << endl;
		f << "           realizations of measurement noise, " << endl;
		f << "           'actual' phi does not." << endl;
		if (echo)
		{
			cout << "     note: 'measured' phi reported above includes " << endl;
			cout << "           realizations of measurement noise, " << endl;
			cout << "           'actual' phi does not." << endl;
		}
	}
	

	f << endl << endl;
	f.flush();
}



//void PhiHandler::report(bool echo)
//{
//	ofstream &f = file_manager->rec_ofstream();
//	f << get_summary_header();
//	if (echo)
//		cout << get_summary_header();
//	string s = get_summary_string(PhiHandler::phiType::COMPOSITE);
//	f << s;
//	if (echo)
//		cout << s;
//	s = get_summary_string(PhiHandler::phiType::MEAS);
//	f << s;
//	if (echo)
//		cout << s;
//	if (org_reg_factor != 0.0)
//	{
//		s = get_summary_string(PhiHandler::phiType::REGUL);
//		f << s;
//		if (echo)
//			cout << s;
//	}
//	s = get_summary_string(PhiHandler::phiType::ACTUAL);
//	f << s;
//	if (echo)
//		cout << s;
//	if (org_reg_factor != 0.0)
//	{
//		if (*reg_factor == 0.0)
//		{
//			f << "    (note: reg_factor is zero; regularization phi reported but not used)" << endl;
//			if (echo)
//				cout << "    (note: reg_factor is zero; regularization phi reported but not used)" << endl;
//		}
//		else
//		{
//			f << "     current reg_factor: " << *reg_factor << endl;
//			if (echo)
//				cout << "     current reg_factor: " << *reg_factor << endl;
//		}
//		if (*reg_factor != 0.0)
//		{
//
//			f << "     note: regularization phi reported above does not " << endl;
//			f << "           include the effects of reg_factor, " << endl;
//			f << "           but composite phi does." << endl;
//			if (echo)
//			{
//				cout << "     note: regularization phi reported above does not " << endl;
//				cout << "           include the effects of reg_factor, " << endl;
//				cout << "           but composite phi does." << endl;
//			}
//		}
//	}
//	f << endl << endl;
//	f.flush();
//}

void L2PhiHandler::write(int iter_num, int total_runs, bool write_group)
{
	write_csv(iter_num, total_runs, file_manager->get_ofstream("phi.actual.csv"), phiType::ACTUAL,oreal_names);
	write_csv(iter_num, total_runs, file_manager->get_ofstream("phi.meas.csv"), phiType::MEAS, oreal_names);
	if (pest_scenario->get_pestpp_options().get_ies_reg_factor() != 0.0)
	{
		write_csv(iter_num, total_runs, file_manager->get_ofstream("phi.regul.csv"), phiType::REGUL, preal_names);	
	}
	write_csv(iter_num, total_runs, file_manager->get_ofstream("phi.composite.csv"), phiType::COMPOSITE, oreal_names);
	if (write_group)
		write_group_csv(iter_num, total_runs, file_manager->get_ofstream("phi.group.csv"));
}

void L2PhiHandler::write_group(int iter_num, int total_runs, vector<double> extra)
{
	write_group_csv(iter_num, total_runs, file_manager->get_ofstream("phi.group.csv"),extra);
}

void L2PhiHandler::write_csv(int iter_num, int total_runs, ofstream &csv, phiType pt, vector<string> &names)
{
	map<string, double>* phi_map = get_phi_map(pt);
	map<string, double>::iterator pmi = phi_map->end();
	csv << iter_num << ',' << total_runs;
	map<string, double> stats = get_summary_stats(pt);
	csv << ',' << stats["mean"] << ',' << stats["std"] << ',' << stats["min"] << ',' << stats["max"];
	for (auto &name : names)
	{
		csv << ',';
		if (phi_map->find(name) != pmi)
			csv << phi_map->at(name);
	}
	csv << endl;
	csv.flush();
}

void L2PhiHandler::prepare_csv(ofstream & csv,vector<string> &names)
{
	csv << "iteration,total_runs,mean,standard_deviation,min,max";
	for (auto &name : names)
		csv << ',' << pest_utils::lower_cp(name);
	csv << endl;
}


void L2PhiHandler::write_group_csv(int iter_num, int total_runs, ofstream &csv, vector<double> extra)
{
	//csv << "iteration,total_runs,realiation";
	string oreal, preal;
	for (int ireal = 0; ireal < oreal_names.size(); ireal++)
	{
		oreal = oreal_names[ireal];
		preal = preal_names[ireal];
		if (obs_group_phi_map.find(oreal) == obs_group_phi_map.end())
			continue;

		csv << iter_num << ',' << total_runs << ',' << pest_utils::lower_cp(oreal) << ',' << pest_utils::lower_cp(preal);
		for (auto &e : extra)
			csv << ',' << e;

		for (auto &name : pest_scenario->get_ctl_ordered_obs_group_names())
			if (obs_group_phi_map[oreal].find(name) == obs_group_phi_map[oreal].end())
				csv << ',' << 0.0;
			else
				csv  << ',' << obs_group_phi_map[oreal][name];
		if (org_reg_factor != 0.0)
		{
			for (auto& name : pest_scenario->get_ctl_ordered_par_group_names())
				if (par_group_phi_map[preal].find(name) == par_group_phi_map[preal].end())
					csv << ',' << 0.0;
				else
					csv << ',' << par_group_phi_map[preal][name];
		}
		csv << endl;;
		csv.flush();
	}
}

void L2PhiHandler::prepare_group_csv(ofstream &csv, vector<string> extra)
{
	csv << "iteration,total_runs,obs_realization,par_realization";
	for (auto &name : extra)
		csv << ',' << pest_utils::lower_cp(name);
	for (auto &name : pest_scenario->get_ctl_ordered_obs_group_names())
		csv << ',' << pest_utils::lower_cp(name);
	if (org_reg_factor != 0.0)
	{
		for (auto& name : pest_scenario->get_ctl_ordered_par_group_names())
			csv << ',' << pest_utils::lower_cp(name);
	}
	csv << endl;
}

vector<int> L2PhiHandler::get_idxs_greater_than(double bad_phi, double bad_phi_sigma, ObservationEnsemble &oe)
{
	map<string, double> _meas;
	Eigen::VectorXd q = get_q_vector();
	for (auto &pv : calc_meas(oe, q))
		_meas[pv.first] = pv.second.sum();
	double mean = calc_mean(&_meas);
	double std = calc_std(&_meas);
	vector<int> idxs;
	vector<string> names = oe.get_real_names();
	for (int i=0;i<names.size();i++)
		if ((_meas[names[i]] > bad_phi) || (_meas[names[i]] > mean + (std * bad_phi_sigma)))
			idxs.push_back(i);
	return idxs;
}

map<string, Eigen::VectorXd> L2PhiHandler::calc_meas(ObservationEnsemble & oe, Eigen::VectorXd &q_vec)
{
	map<string, Eigen::VectorXd> phi_map;
	Eigen::VectorXd oe_base_vec, oe_vec, diff, w_vec;
	//Eigen::VectorXd q = get_q_vector();
	vector<string> act_obs_names = pest_scenario->get_ctl_ordered_nz_obs_names();
	vector<string> base_real_names = oe_base->get_real_names(), oe_real_names = oe.get_real_names();
	vector<string>::iterator start = base_real_names.begin(), end = base_real_names.end();

	
	double phi;
	string rname;

	if (act_obs_names.size() == 0)
	{
		for (auto name : oe.get_real_names())
			phi_map[name] = Eigen::VectorXd();
		return phi_map;
	}

	Eigen::MatrixXd resid = get_obs_resid(oe);
	ObservationInfo oi = pest_scenario->get_ctl_observation_info();
	vector<string> names = oe_base->get_var_names();
	w_vec.resize(names.size());
	for (int i=0;i<names.size();i++)
	{
		w_vec(i) = oi.get_weight(names[i]);
	}
	
	assert(oe_real_names.size() == resid.rows());
	for (int i = 0; i<resid.rows(); i++)
	{
		rname = oe_real_names[i];
		if (find(start, end, rname) == end)
			continue;
		diff = resid.row(i);
		
		diff = diff.cwiseProduct(w_vec);
		
		phi = (diff.cwiseProduct(diff)).sum();
		phi_map[rname] = diff.cwiseProduct(diff);
	}
	return phi_map;
}

map<string, Eigen::VectorXd> L2PhiHandler::calc_regul(ParameterEnsemble & pe)
{
	map<string, Eigen::VectorXd> phi_map;
	vector<string> real_names = pe.get_real_names();
	pe_base->transform_ip(ParameterEnsemble::transStatus::NUM);
	pe.transform_ip(ParameterEnsemble::transStatus::NUM);
	Eigen::MatrixXd diff_mat = get_par_resid(pe);


	Eigen::VectorXd diff;
	for (int i = 0; i < real_names.size(); i++)
	{
		diff = diff_mat.row(i);
		diff = diff.cwiseProduct(diff);
		//cout << diff << endl;
		diff = diff.cwiseProduct(parcov_inv_diag);
		//cout << diff << endl;
		//cout << parcov_inv_diag << endl;
		//phi_map[real_names[i]] = _reg_fac * diff;
		phi_map[real_names[i]] = diff;
		
	}
	return phi_map;
}


void L2PhiHandler::apply_ineq_constraints(Eigen::MatrixXd &resid, vector<string> &names)
{
	
	//vector<string> names = oe_base->get_var_names();
	//vector<string> lt_names = get_lt_obs_names(), gt_names = get_gt_obs_names();
	//vector<string> act_obs_names = pest_scenario->get_ctl_ordered_nz_obs_names();
	
	assert(names.size() == resid.cols());

	map<string, double> lt_vals,gt_vals;
	Observations obs = pest_scenario->get_ctl_observations();
	for (auto &n : lt_obs_names)
		lt_vals[n] = obs.get_rec(n);
	for (auto &n : gt_obs_names)
		gt_vals[n] = obs.get_rec(n);
	if ((lt_vals.size() == 0) && (gt_vals.size() == 0))
		return;
	map<string, int> idxs;
	//for (int i = 0; i < act_obs_names.size(); i++)
	//	idxs[act_obs_names[i]] = i;
	for (int i = 0; i < names.size(); i++)
		idxs[names[i]] = i;
	int idx;
	double val;
	Eigen::VectorXd col;

	for (auto iv : lt_vals)
	{
		idx = idxs[iv.first];
		col = resid.col(idx);
		val = iv.second;
		//cout << resid.col(idx) << endl;
		for (int i = 0; i < resid.rows(); i++)
			col(i) = (col(i) < 0.0) ? 0.0 : col(i);
		//cout << resid.col(idx) << endl;
		resid.col(idx) = col;
		//cout << resid.col(idx) << endl;
	}

	for (auto iv : gt_vals)
	{
		idx = idxs[iv.first];
		col = resid.col(idx);
		val = iv.second;
		for (int i = 0; i < resid.rows(); i++)
			col(i) = (col(i) > 0.0) ? 0.0 : col(i);
		resid.col(idx) = col;
	}
}


map<string, Eigen::VectorXd> L2PhiHandler::calc_actual(ObservationEnsemble & oe, Eigen::VectorXd &q_vec)
{
	map<string, Eigen::VectorXd> phi_map;
	Eigen::MatrixXd resid = get_actual_obs_resid(oe);
	vector<string> base_real_names = oe_base->get_real_names(), oe_real_names = oe.get_real_names();
	vector<string>::iterator start = base_real_names.begin(), end = base_real_names.end();
	double phi;
	string rname;

	Eigen::MatrixXd oe_reals = oe.get_eigen(vector<string>(), oe_base->get_var_names());
	Eigen::VectorXd diff;
	for (int i = 0; i<oe.shape().first; i++)
	{
		rname = oe_real_names[i];
		if (find(start, end, rname) == end)
			continue;
		//diff = (oe_vec - obs_val_vec).cwiseProduct(q);
		diff = resid.row(i);
		diff = diff.cwiseProduct(q_vec);
		//phi = (diff.cwiseProduct(diff)).sum();
		phi_map[rname] = diff.cwiseProduct(diff);
	}
	return phi_map;
}


map<string, double> L2PhiHandler::calc_composite(map<string, double> &_meas, map<string, double> &_regul)
{
	map<string, double> phi_map;
	string prn, orn;
	double reg, mea;
	map<string, double>::iterator meas_end = _meas.end(), regul_end = _regul.end();
	for (int i = 0; i < oe_base->shape().first; i++)
	{
		prn = preal_names[i];
		orn = oreal_names[i];
		if (meas.find(orn) != meas_end)
		{
			mea = _meas[orn];
			reg = _regul[prn];
			phi_map[orn] = mea + (reg * org_reg_factor);
		}
	}
	return phi_map;
}


ParChangeSummarizer::ParChangeSummarizer(ParameterEnsemble *_base_pe_ptr, FileManager *_file_manager_ptr, OutputFileWriter* _output_file_writer_ptr)
	//base_pe_ptr(_base_pe), file_manager_ptr(_file_manager)
{
	base_pe_ptr = _base_pe_ptr;
	file_manager_ptr = _file_manager_ptr;
	output_file_writer_ptr = _output_file_writer_ptr;
	init_moments = base_pe_ptr->get_moment_maps();
	ParameterGroupInfo gi = base_pe_ptr->get_pest_scenario().get_base_group_info();
	string group;
	for (auto &n : base_pe_ptr->get_var_names())
	{
		group = gi.get_group_name(n);
		if (pargp2par_map.find(group) == pargp2par_map.end())
			pargp2par_map[group] = set<string>();	
		pargp2par_map[group].emplace(n);
	}
}


void ParChangeSummarizer::summarize(ParameterEnsemble &pe, int iiter, string filename)
{

	stringstream ss;
	ofstream &frec = file_manager_ptr->rec_ofstream();
	ss << endl << "   ---  Parameter Group Change Summmary  ---    " << endl;
	ss << "   (compared to the initial ensemble using active realizations)" << endl;
	cout << ss.str();
	frec << ss.str();
	ss.str("");
	ss << setw(15) << "group" << setw(12) << "mean change" << setw(12) << "std change" << setw(18) << "num at/near bnds" << setw(16) << "% at/near bnds" << endl;
	cout << ss.str();
	frec << ss.str();	
	update(pe);
	vector<string> grp_names;// = base_pe_ptr->get_pest_scenario().get_ctl_ordered_par_group_names();
	vector<pair<double, string>> mean_pairs;
	for (auto m : mean_change)
	{
		mean_pairs.push_back(pair<double, string>(abs(m.second), m.first));
	}
	sort(mean_pairs.begin(), mean_pairs.end());
	//for (auto m : mean_pairs)
	for (int i = mean_pairs.size() - 1; i >= 0; i--)
	{
		grp_names.push_back(mean_pairs[i].second);
	}

	int i = 0;
	for (auto &grp_name : grp_names)
	{
		double mean_diff = mean_change[grp_name];
		double std_diff = std_change[grp_name];
		int num_out = num_at_bounds[grp_name];
		int percent_out = percent_at_bounds[grp_name];
		ss.str("");
		ss << setw(15) << pest_utils::lower_cp(grp_name) << setw(12) << mean_diff * 100.0 << setw(12) << std_diff * 100.0 << setw(18);
		ss << num_out << setw(16) << setprecision(2) << percent_out << endl;
		if (i < 15)
			cout << ss.str();
		frec << ss.str();
		i++;
	}

	ss.str("");
	ss << "    Note: parameter change summary sorted according to abs 'mean change'." << endl;
	cout << ss.str();
	frec << ss.str();
	if (grp_names.size() > 15)
	{
		ss.str("");
		ss << "    Note: Only the first 15 parameter groups shown, see rec file for full listing" << endl;
		cout << ss.str();
	}

	cout << endl;
	frec << endl;

	if (filename.size() > 0)
		write_to_csv(filename);

}

void ParChangeSummarizer::write_to_csv(string& filename)
{
	ofstream f(filename);
	if (f.bad())
		throw runtime_error("ParChangeSummarizer::write_to_csv() error opening file " + filename);

	f << "group,mean_change,std_change,num_at_near_bounds,percent_at_near_bounds" << endl;
	for (auto grp_name : base_pe_ptr->get_pest_scenario_ptr()->get_ctl_ordered_par_group_names())
	{
		f << pest_utils::lower_cp(grp_name) << "," << mean_change[grp_name] << "," << std_change[grp_name] << ",";
		f << num_at_bounds[grp_name] << "," << percent_at_bounds[grp_name] << endl;
	}
	f.close();
	file_manager_ptr->rec_ofstream() << "...saved parameter change summary to " << filename << endl;
	cout << "...saved parameter change summary to " << filename << endl;

}


void ParChangeSummarizer::update(ParameterEnsemble& pe)
{
	mean_change.clear();
	std_change.clear();
	num_at_bounds.clear();
	percent_at_bounds.clear();
	pair<map<string, double>, map<string, double>> moments = pe.get_moment_maps();
	init_moments = base_pe_ptr->get_moment_maps(pe.get_real_names());
	
	double mean_diff = 0.0, std_diff = 0.0;
	double dsize, value1, value2, v;
	vector<string> pnames = pe.get_var_names();
	Parameters lb = pe.get_pest_scenario_ptr()->get_ctl_parameter_info().get_low_bnd(pnames);
	pe.get_pest_scenario_ptr()->get_base_par_tran_seq().active_ctl2numeric_ip(lb);
	Parameters ub = pe.get_pest_scenario_ptr()->get_ctl_parameter_info().get_up_bnd(pnames);
	pe.get_pest_scenario_ptr()->get_base_par_tran_seq().active_ctl2numeric_ip(ub);
	vector<string> grp_names = base_pe_ptr->get_pest_scenario().get_ctl_ordered_par_group_names();
	map<string, int> idx_map;
	for (int i = 0; i < pnames.size(); i++)
		idx_map[pnames[i]] = i;
	int num_out, num_pars;
	int num_reals = pe.get_real_names().size();
	Eigen::ArrayXd arr;
	for (auto& grp_name : grp_names)
	{
		mean_diff = 0.0, std_diff = 0.0;
		num_pars = pargp2par_map[grp_name].size();
		num_out = 0;
		for (auto& par_name : pargp2par_map[grp_name])
		{
			arr = pe.get_eigen_ptr()->col(idx_map[par_name]).array();
			for (int i = 0; i < num_reals; i++)
			{
				v = arr[i];
				if ((v > (ub[par_name] * 1.01)) || (v < (lb[par_name] * 0.99)))
					num_out++;
			}
			value1 = init_moments.first[par_name];
			value2 = value1 - moments.first[par_name];
			if ((value1 != 0.0) && (value2 != 0.0))
				mean_diff += value2 / value1;
			value1 = init_moments.second[par_name];
			value2 = value1 - moments.second[par_name];
			if ((value1 != 0.0) && (value2 != 0.0))
				std_diff += value2 / value1;
		}
		dsize = double(init_moments.first.size());
		if (mean_diff != 0.0)
			mean_diff = mean_diff / dsize;
		if (std_diff != 0.0)
			std_diff = std_diff / dsize;

		double percent_out = 0;
		if (num_pars > 0)
			percent_out = double(num_out) / double(num_pars * num_reals) * 100;

		mean_change[grp_name] = mean_diff;
		std_change[grp_name] = std_diff;
		num_at_bounds[grp_name] = num_out;
		percent_at_bounds[grp_name] = percent_out;

	}
}


void save_base_real_par_rei(Pest& pest_scenario, ParameterEnsemble& pe, ObservationEnsemble& oe,
	OutputFileWriter& output_file_writer, FileManager& file_manager, int iter)
{
	stringstream ss;
	map<string, int> vmap = pe.get_real_map();
	if (vmap.find(BASE_REAL_NAME) != vmap.end())
	{
		ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();
		Parameters pars;
		pars.update(pe.get_var_names(), eigenvec_2_stlvec(pe.get_real_vector(BASE_REAL_NAME)));
		if (pe.get_trans_status() == ParameterEnsemble::transStatus::NUM)
			pts.numeric2ctl_ip(pars);
		// save parameters to .par file
		if (iter >= 0)
			ss << iter << ".";
		ss << "base.par";
		output_file_writer.write_par(file_manager.open_ofile_ext(ss.str()), pars, *(pts.get_offset_ptr()),
			*(pts.get_scale_ptr()));
		file_manager.close_file("par");

		vmap = oe.get_real_map();
		if (vmap.find(BASE_REAL_NAME) == vmap.end())
		{
			//message(2, "unable to find 'BASE' realization in obs ensemble for saving .base.rei file, continuing...");
		}
		else
		{
			Observations obs;
			obs.update(oe.get_var_names(), eigenvec_2_stlvec(oe.get_real_vector(BASE_REAL_NAME)));
			ObjectiveFunc obj_func(&(pest_scenario.get_ctl_observations()), &(pest_scenario.get_ctl_observation_info()), &(pest_scenario.get_prior_info()));
			// save new residuals to .rei file
			ss.str("");
			if (iter >= 0)
				ss << iter << ".";
			ss << "base.rei";
			output_file_writer.write_rei(file_manager.open_ofile_ext(ss.str()), iter,
				pest_scenario.get_ctl_observations(), obs, obj_func, pars);
		}
	}

}


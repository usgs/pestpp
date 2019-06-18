#include <random>
#include <iomanip>
#include <unordered_set>
#include <iterator>
#include "Ensemble.h"
#include "RestartController.h"
#include "utilities.h"
#include "ParamTransformSeq.h"
#include "ObjectiveFunc.h"
#include "RedSVD-h.h"
#include "covariance.h"
#include "PerformanceLog.h"
#include "system_variables.h"
#include "Localizer.h"

bool Localizer::initialize(PerformanceLog *performance_log)
{
	stringstream ss;
	how == How::OBSERVATIONS; //set this for the case with no localization
	string how_str = pest_scenario_ptr->get_pestpp_options().get_ies_localize_how();
	if (how_str[0] == 'P')
	{
		how = How::PARAMETERS;
	}
	else if (how_str[0] == 'O')
	{
		how = How::OBSERVATIONS;
	}
	else
	{
		throw runtime_error("Localizer error: 'ies_localize_how' must start with 'P' (pars) or 'O' (obs) not " + how_str[0]);
	}
	string filename = pest_scenario_ptr->get_pestpp_options().get_ies_localizer();
	autoadaloc = pest_scenario_ptr->get_pestpp_options().get_ies_autoadaloc();
	sigma_dist = pest_scenario_ptr->get_pestpp_options().get_ies_autoadaloc_sigma_dist();
	use = true;
	if ((filename.size() == 0) && (!autoadaloc))
	{
		use = false;
		return use;
	}
	//use = true;

	if (filename.size() == 0)
		return use;
	performance_log->log_event("loading localizer matrix from file " + filename);
	mat.from_file(filename);
	
	performance_log->log_event("processing localizer matrix");
	process_mat(performance_log);
	if (autoadaloc)
	{
		//string how = pest_scenario_ptr->get_pestpp_options().get_ies_localize_how();
		if (how != How::PARAMETERS)
			throw runtime_error("using a localizer matrix and autoadaloc requires ies_localize_how == 'PARAMETERS'");
		for (auto i : localizer_map)
		{
			set<string> oset(i.second.first.begin(), i.second.first.end());
			for (auto pname : i.second.second)
			{
				listed_obs[pname] = oset;
			}
		}

	}
	return use;
}


void Localizer::process_mat(PerformanceLog *performance_log)
{
	stringstream ss;
	

	//error checking and building up container of names
	vector<string> names = pest_scenario_ptr->get_ctl_ordered_adj_par_names();
	set<string> par_names(names.begin(), names.end());
	names = pest_scenario_ptr->get_ctl_ordered_nz_obs_names();
	set<string> obs_names(names.begin(), names.end());

	map<string, vector<string>> pargp_map;
	ParameterGroupInfo *pi = pest_scenario_ptr->get_base_group_info_ptr();

	for (auto &pg : pest_scenario_ptr->get_ctl_ordered_par_group_names())
	{

		names.clear();
		for (auto &p : par_names)
			if (pi->get_group_name(p) == pg)
				names.push_back(p);
		pargp_map[pg] = names;
	}

	map<string, vector<string>> obgnme_map;
	for (auto &og : pest_scenario_ptr->get_ctl_ordered_obs_group_names())
	{
		ObservationInfo *oi = pest_scenario_ptr->get_observation_info_ptr();
		names.clear();
		for (auto &o : obs_names)
			if (oi->get_group(o) == og)
				names.push_back(o);
		obgnme_map[og] = names;
	}

	vector<string> missing, dups, not_allowed;
	vector<vector<string>> obs_map;
	set<string> dup_check;

	//for (auto &o : mat.get_row_names())
	string o;
	vector<string> row_names = mat.get_row_names();
	for (int i=0;i<mat.nrow();i++)
	{
		o = row_names[i];
		if (obs_names.find(o) != obs_names.end())
		{
			obs2row_map[o] = i;
			obs_map.push_back(vector<string>{o});
			if (dup_check.find(o) != dup_check.end())
				dups.push_back(o);
			dup_check.emplace(o);
		}
		else if (obgnme_map.find(o) != obgnme_map.end())
		{
			obs_map.push_back(obgnme_map[o]);
			if (obgnme_map[o].size() == 0)
				throw runtime_error("Localizer::process_mat() error: listed observation group '" + o + "' has no non-zero weight observations");
			for (auto &oo : obgnme_map[o])
			{
				obs2row_map[oo] = i;
				if (dup_check.find(oo) != dup_check.end())
					dups.push_back(oo);
				dup_check.emplace(oo);
				if (obs_names.find(oo) == obs_names.end())
					not_allowed.push_back(oo);

			}
		}
		else
			missing.push_back(o);
	}
	if (not_allowed.size() > 0)
	{
		ss.str("");
		ss << "Localizer::process_mat() error: the following obs names were identified through obs groups but are not in the non-zero weight obs names: ";
		for (auto &oo : not_allowed)
			ss << oo << ",";
		throw runtime_error(ss.str());
	}
		
	if (missing.size() > 0)
	{
		ss << "Localizer::process_mat() error:  the following rows in " << filename << " were not found in the non-zero-weight observation names or observation group names: ";
		for (auto &m : missing)
			ss << m << ',';
		performance_log->log_event("error:" + ss.str());
		throw runtime_error(ss.str());
	}
	if (dups.size() > 0)
	{
		ss << "Localizer::process_mat() error:  the following observations were listed more than once (possibly through an obs group): ";
		for (auto & d : dups)
			ss << d << ',';
		performance_log->log_event("error:" + ss.str());
		throw runtime_error(ss.str());
	}


	vector<string> col_names = mat.get_col_names();
	vector<vector<string>> par_map;
	dup_check.clear();
	string p;
	
	//for (auto &p : mat.get_col_names())
	for (int i=0;i<mat.ncol();++i)
	{
		p = col_names[i];
		if (par_names.find(p) != par_names.end())
		{
			par2col_map[p] = i;
			par_map.push_back(vector<string>{p});
			if (dup_check.find(p) != dup_check.end())
				dups.push_back(p);
			dup_check.emplace(p);
		}
		else if (pargp_map.find(p) != pargp_map.end())
		{
			par_map.push_back(pargp_map[p]);
			if (pargp_map[p].size() == 0)
				throw runtime_error("Localizer::process_mat() error:  listed parameter group '" + p + "' has no adjustable parameters");
			for (auto &pp : pargp_map[p])
			{
				par2col_map[pp] = i;
				if (dup_check.find(pp) != dup_check.end())
					dups.push_back(pp);
				dup_check.emplace(pp);
				if (par_names.find(pp) == par_names.end())
					not_allowed.push_back(pp);

			}
		}
		else
			missing.push_back(p);
	}
	if (not_allowed.size() > 0)
	{
		ss.str("");
		ss << "Localizer::process_mat() error: the following par names were identified through par groups but are not in the adj par names: ";
		for (auto &pp : not_allowed)
			ss << pp << ",";
		throw runtime_error(ss.str());

	}
	if (missing.size() > 0)
	{
		ss << "Localizer::process_mat() error: the following cols in " << filename << " were not found in the active parameter names or parameter group names: ";
		for (auto &m : missing)
			ss << m << ',';
		performance_log->log_event("error:" + ss.str());
		throw runtime_error(ss.str());
	}
	if (dups.size() > 0)
	{
		ss << "Localizer::process_mat() error: the following parameters were listed more than once (possibly through a par group): ";
		for (auto & d : dups)
			ss << d << ',';
		performance_log->log_event("error:" + ss.str());
		throw runtime_error(ss.str());
	}


	//map all the nz locations in the matrix
	map<int, vector<int>> idx_map;
	vector<string> vobs, vpar;

	//int i, j;
	if (how == How::PARAMETERS)
	{
		vector<string> col_names = mat.get_col_names();
		for (int k = 0; k < mat.e_ptr()->outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(*mat.e_ptr(), k); it; ++it)
			{
				if (idx_map.find(it.col()) == idx_map.end())
					idx_map[it.col()] = vector<int>{ it.row() };
				else
					idx_map[it.col()].push_back(it.row());
				//std::cout << "(" << it.row() << ","; // row index
				//std::cout << it.col() << ")\t"; // col index (here it is equal to k)

			}
		}
		//populate the localizer map
		for (auto &idx : idx_map)
		{
			vpar = par_map[idx.first];
			vobs.clear();
			for (auto &i : idx.second)
			{
				vobs.insert(vobs.end(), obs_map[i].begin(), obs_map[i].end());
			}
			pair<vector<string>, vector<string>> p(vobs, vpar);
			//localizer_map.push_back(p);
			localizer_map[col_names[idx.first]] = p;
		}

	}
	else
	{
		vector<string> row_names = mat.get_row_names();
		for (int k = 0; k < mat.e_ptr()->outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(*mat.e_ptr(), k); it; ++it)
			{
				if (idx_map.find(it.row()) == idx_map.end())
					idx_map[it.row()] = vector<int>{ it.col() };
				else
					idx_map[it.row()].push_back(it.col());
				//std::cout << "(" << it.row() << ","; // row index
				//std::cout << it.col() << ")\t"; // col index (here it is equal to k)

			}
		}
		//populate the localizer map
		for (auto &idx : idx_map)
		{
			vobs = obs_map[idx.first];
			vpar.clear();
			for (auto &i : idx.second)
			{
				vpar.insert(vpar.end(), par_map[i].begin(), par_map[i].end());
			}
			pair<vector<string>, vector<string>> p(vobs, vpar);
			//localizer_map.push_back(p);
			localizer_map[row_names[idx.first]] = p;
		}

	}

	
}


void Localizer::report(ofstream &f_rec)
{
	vector<string> zeros;
	vector<string> col_names = mat.get_col_names();

	double sum;
	for (int i = 0; i < mat.e_ptr()->cols(); i++)
	{
		sum = mat.e_ptr()->col(i).sum();
		if (sum == 0.0)
			zeros.push_back(col_names[i]);
	}
	if (zeros.size() > 0)
	{
		cout << "Note: " << zeros.size() << " parameters/parameter groups in the localizer have no non-zero entries";
		if (pest_scenario_ptr->get_pestpp_options().get_ies_verbose_level() > 1)
		{
			cout << ", see.rec file for listing" << endl;
			f_rec << "Note: the following parameters/parmaeter groups have no non-zero entries in the localizer meaning they will not be adjusted:" << endl;
			int i = 0;
			for (auto z : zeros)
			{
				f_rec << z << " ";
				i++;
				if (i > 10)
				{
					f_rec << endl;
					i = 0;
				}
			}
			f_rec << endl;
		}
		else
		{
			cout << endl;
		}
	}
	zeros.clear();
	vector<string> row_names = mat.get_row_names();
	for (int i = 0; i < mat.e_ptr()->rows(); i++)
	{
		sum = mat.e_ptr()->row(i).sum();
		if (sum == 0.0)
			zeros.push_back(row_names[i]);
	}
	if (zeros.size() > 0)
	{
		cout << "Note: " << zeros.size() << " observations/observation groups in the localizer have no non-zero entries, see .rec file for listing" << endl;

		f_rec << "WARNING: the following observations/observation groups have no non-zero entries in the localizer meaning their residuals cannot be reduced:" << endl;
		int i = 0;
		for (auto z : zeros)
		{
			f_rec << z << " ";
			i++;
			if (i > 10)
			{
				f_rec << endl;
				i = 0;
			}
		}
		f_rec << endl;

	}

}

unordered_map<string, pair<vector<string>, vector<string>>> Localizer::get_localizer_map(int iter, ObservationEnsemble &oe, ParameterEnsemble &pe, PerformanceLog *performance_log)
{
	if (!autoadaloc)
		return localizer_map;

	cout << "...starting automatic adaptive localization calculations" << endl;
	stringstream ss;
	performance_log->log_event("autoadaloc: checking for compat with par and obs ensembles");
	if (pe.shape().first != oe.shape().first)
	{
		ss.str("");
		ss << "Localizer::get_localizer_map() Error: for autoadaloc, pe must have same number of reals (" << pe.shape().first << ") as oe (" << oe.shape().first << ")";
		performance_log->log_event(ss.str());
		throw runtime_error(ss.str());
	}

	
	vector<string> par_names = pe.get_pest_scenario_ptr()->get_ctl_ordered_adj_par_names(), obs_names = pe.get_pest_scenario_ptr()->get_ctl_ordered_nz_obs_names();

	performance_log->log_event("autoadaloc: calculating correlation coefficients");
	Eigen::MatrixXd pe_diff = pe.get_eigen_anomalies(vector<string>(), par_names), oe_diff = oe.get_eigen_anomalies(vector<string>(), obs_names);
	double scale = 1.0 / double(pe.shape().first - 1);
	Eigen::ArrayXd par_std = ((pe_diff.array().square().colwise().sum()) * scale).cwiseSqrt();
	Eigen::ArrayXd obs_std = ((oe_diff.array().square().colwise().sum()) * scale).cwiseSqrt();

	int npar = par_names.size(), nobs = obs_names.size(), nreals = pe.shape().first;


	//construct the permutation matrix
	/*Eigen::Matrix<int, Eigen::Dynamic, 1> temp(nreals);
	temp[0] = nreals - 1;
	for (int i = 0; i < nreals - 1; i++)
	{
		temp[i + 1] = i;
	}
	Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(temp);*/

	int ies_verbose = pe.get_pest_scenario_ptr()->get_pestpp_options().get_ies_verbose_level();

	string pst_filename = pe.get_pest_scenario_ptr()->get_pst_filename();
	ofstream f_out;
	if (ies_verbose > 1)
	{
		ss.str("");
		ss << pst_filename.substr(0, pst_filename.size() - 4) << "." << iter << ".autoadaloc.csv";
		string filename = ss.str();
		f_out = ofstream(filename);
		if (!f_out.good())
			throw runtime_error("autoadaloc error opening filename " + filename + " for verbose output");
		f_out << "obsnme,parnme,correlation_coeff,background_mean,background_stdev,threshold,kept";
		for (int i = 0; i < nreals - 1; i++)
			f_out << "," << i;
		f_out << endl;
	}

	//here we go...
	vector<int> par_indices;
	for (int jpar = 0; jpar < npar; jpar++)
		par_indices.push_back(jpar);
	vector<Eigen::Triplet<double>> triplets;
	AutoAdaLocThread worker(performance_log, &f_out, iter, ies_verbose, npar, nobs, par_indices, pe_diff, oe_diff, par_std, obs_std, par_names, obs_names, triplets,sigma_dist,listed_obs);

	int num_threads = pe.get_pest_scenario_ptr()->get_pestpp_options().get_ies_num_threads();

	if (num_threads < 1)
		worker.work(0);
	else
	{
		Eigen::setNbThreads(1);
		vector<thread> threads;
		vector<exception_ptr> exception_ptrs;
		performance_log->log_event("launching autoadaloc threads");

		for (int i = 0; i < num_threads; i++)
		{
			exception_ptrs.push_back(exception_ptr());
		}

		for (int i = 0; i < num_threads; i++)
		{
			//threads.push_back(thread(&LocalUpgradeThread::work, &worker, i, iter, cur_lam));

			threads.push_back(thread(aal_upgrade_thread_function, i, std::ref(worker), std::ref(exception_ptrs[i])));

		}
		performance_log->log_event("waiting to join threads");
		//for (auto &t : threads)
		//	t.join();
		for (int i = 0; i < num_threads; ++i)
		{
			if (exception_ptrs[i])
			{
				try
				{
					rethrow_exception(exception_ptrs[i]);
				}
				catch (const std::exception& e)
				{
					ss.str("");
					ss << "thread " << i << "raised an exception: " << e.what();
					throw runtime_error(ss.str());
				}
			}
			threads[i].join();
		}
		performance_log->log_event("threaded localized upgrade calculation done");

	}

	if (triplets.size() == 0)
	{
		ss.str("");
		ss << "autoadaloc error: no non-zero entries in localization matrix, can not continue";
		throw runtime_error(ss.str());
	}

	if (ies_verbose > 1)
	{
		f_out.close();

		performance_log->log_event("instantiating thresholded correlation coefficient matrix from triplets");
		mat.from_triplets(obs_names, par_names, triplets);
		ss.str("");
		ss << pst_filename.substr(0, pst_filename.size() - 4) << "." << iter << ".autoadaloc.tCC";
		string filename = ss.str();
		if (pe.get_pest_scenario_ptr()->get_pestpp_options().get_ies_save_binary())
		{
			filename = filename + ".jcb";
			performance_log->log_event("saving thresholded CC matrix in binary format to " + filename);
			mat.to_binary_new(filename);
		}
		else
		{
			filename = filename + ".mat";
			performance_log->log_event("saving thresholded CC matrix in ASCII format to " + filename);
			mat.to_ascii(filename);
		}
	}


	//convert triplet values to 1's
	//and check for pars and obs with zero entries
	map<int, int> par_count, obs_count;
	for (int i = 0; i < par_names.size(); i++)
		par_count[i] = 0;
	for (int i = 0; i < obs_names.size(); i++)
		obs_count[i] = 0;
	vector<Eigen::Triplet<double>> ones;
	for (auto &t : triplets)
	{
		ones.push_back(Eigen::Triplet<double>(t.row(), t.col(), abs(t.value())));
		par_count[t.col()]++;
		obs_count[t.row()]++;
	}
	triplets.clear();
	vector<int> zeros;
	for (auto p : par_count)
		if (p.second == 0)
			zeros.push_back(p.first);
	if (zeros.size() > 0)
		cout << "Note: " << zeros.size() << " parameters have no nonzero entries in autoadaloc localizer" << endl;

	zeros.clear();
	for (auto o : obs_count)
		if (o.second == 0)
			zeros.push_back(o.first);

	if (zeros.size() > 0)
		cout << "Note: " << zeros.size() << " observations have no nonzero entries in autoadaloc localizer" << endl;





	performance_log->log_event("instantiating localizer matrix from triplets");
	mat.from_triplets(obs_names, par_names, ones);

	ss.str("");
	ss << "autoadaloc matrix constructed with " << mat.e_ptr()->nonZeros() << " non-zero elements";
	performance_log->log_event(ss.str());
	process_mat(performance_log);
	cout << "automatic adaptive localization calculations done" << endl;
	return localizer_map;
	//cout << endl;



	//Eigen::Array<double, 1, Eigen::Dynamic> std_dev = ((oe.get_eigen_ptr()->rowwise() - oe.get_eigen_ptr()->colwise().mean()).array().square().colwise().sum() / (oe.get_eigen_ptr()->rows - 1)).sqrt();

}



void aal_upgrade_thread_function(int id, AutoAdaLocThread &worker, exception_ptr &eptr)
{
	try
	{
		worker.work(id);
	}
	catch (...)
	{
		eptr = current_exception();
	}

	return;
}



AutoAdaLocThread::AutoAdaLocThread(PerformanceLog *_performance_log, ofstream *_f_out, int _iter, int _ies_verbose, int _npar, int _nobs, vector<int> &_par_indices,
	Eigen::MatrixXd &_pe_diff, Eigen::MatrixXd &_oe_diff, Eigen::ArrayXd &_par_std, Eigen::ArrayXd &_obs_std, vector<string> &_par_names, vector<string> &_obs_names,
	vector<Eigen::Triplet<double>> &_triplets, double _sigma_dist,map<string,set<string>> &_list_obs): pe_diff(_pe_diff), oe_diff(_oe_diff), par_indices(_par_indices), par_std(_par_std), obs_std(_obs_std),par_names(_par_names),
	obs_names(_obs_names),triplets(_triplets), list_obs(_list_obs)
{
	iter = _iter;
	npar = _npar;
	nobs = _nobs;
	ies_verbose = _ies_verbose;
	performance_log = _performance_log;
	f_out = _f_out;
	sigma_dist = _sigma_dist;
	
	for (int i = 0; i < obs_names.size(); i++)
		idx2obs[i] = obs_names[i];

}

void AutoAdaLocThread::work(int thread_id)
{

	stringstream ss;
	int jpar, pcount = 0, nreals = pe_diff.rows();
	double cc, bg_cc, bg_mean, bg_std, thres, t;
	double sign;
	double scale = 1.0 / double(nreals - 1);
	double pstd, ostd;
	//vector<Eigen::Triplet<double>> triplets, ones;
	Eigen::VectorXd par_ss, obs_ss, obs_ss_shift;

	par_ss.resize(nreals);
	obs_ss.resize(nreals);
	obs_ss_shift.resize(nreals);
	Eigen::ArrayXd bg_cc_vec(nreals - 1);
	set<string> sobs;
	unique_lock<mutex> par_indices_guard(par_indices_lock, defer_lock);
	unique_lock<mutex> oe_diff_gaurd(oe_diff_lock, defer_lock);
	unique_lock<mutex> f_out_guard(f_out_lock, defer_lock);
	unique_lock<mutex> triplets_guard(triplets_lock, defer_lock);
	unique_lock<mutex>pfm_guard(pfm_lock, defer_lock);
	unique_lock<mutex>par_names_guard(par_names_lock, defer_lock);
	bool use_list_obs = true;
	while (true)
	{

		while (true)
		{
			if (par_indices_guard.try_lock())
			{
				if (par_indices.size() == 0)
				{
					ss.str("");
					ss << "autoadaloc thread: " << thread_id << " processed " << pcount << " parameters ";
					if (ies_verbose > 1)
					{
						cout << ss.str() << endl;
					}
					while (true)
					{
						if (pfm_guard.try_lock())
						{
							performance_log->log_event(ss.str());
							pfm_guard.unlock();
							break;
						}
					}
					par_indices_guard.unlock();
					return;
				}
				if (par_indices.size() % 10000 == 0)
				{
					ss.str("");
					ss << "autoadaloc iter " << iter << " progress: " << par_indices.size() << " of " << npar << " parameters done";
					while (true)
					{
						if (pfm_guard.try_lock())
						{
							performance_log->log_event(ss.str());
							pfm_guard.unlock();
							break;
						}
					}
					if (ies_verbose > 1)
						cout << ss.str() << endl;
				}
				
				jpar = par_indices[par_indices.size() - 1];
				par_indices.pop_back();
				if (par_std[jpar] == 0.0)
				{
					par_indices_guard.unlock();
					continue;
				}
				par_ss = pe_diff.col(jpar) * (1.0 / par_std[jpar]);
				
				if (list_obs.size() > 0)
				{
					sobs = list_obs[par_names[jpar]];
					use_list_obs = true;
				}
				else
				{
					use_list_obs = false;
				}
				pcount++;
				par_indices_guard.unlock();
				break;
			}
		}

		
		string oname;
		bool no_obs = true;
		for (int iobs = 0; iobs < nobs; iobs++)
		{
			while (true)
			{
				/*if (oe_diff_gaurd.try_lock())
				{
					if (obs_std[iobs] == 0.0)
					{
						oe_diff_gaurd.unlock();
						continue;
					}
					string oname = obs_names[iobs];
					if ((sobs.size() > 0) && (sobs.find(oname) == sobs.end()))
					{
						oe_diff_gaurd.unlock();
						continue;
					}
					obs_ss = oe_diff.col(iobs) * (1.0 / obs_std[iobs]);
					oe_diff_gaurd.unlock();
					break;*/
				if (oe_diff_gaurd.try_lock())
				{
					
					obs_ss = oe_diff.col(iobs);
					oname = obs_names[iobs];
					oe_diff_gaurd.unlock();
					break;
				}

			}
			if (obs_std[iobs] == 0.0) 
			{
				continue;
			}
			
			if ((use_list_obs) && (sobs.size() == 0))
				continue;

			if ((sobs.size() > 0) && (sobs.find(oname) == sobs.end()))
			{
				continue;
			}
			obs_ss = obs_ss * (1.0 / obs_std[iobs]);
			cc = (par_ss.transpose() * obs_ss)[0] * scale;
			obs_ss_shift = 1.0 * obs_ss; //force a copy
			for (int ireal = 0; ireal < nreals - 1; ireal++)
			{

				//obs_ss_shift.transpose() = obs_ss_shift.transpose() * perm;
				//circular shift
				t = obs_ss_shift[nreals - 1];
				for (int i = nreals - 1; i > 0; i--)
					obs_ss_shift[i] = obs_ss_shift[i - 1];
				obs_ss_shift[0] = t;

				//bg_cc = (par_ss.transpose() * obs_ss_shift)[0] * scale;
				bg_cc_vec[ireal] = (par_ss.transpose() * obs_ss_shift)[0] * scale;
				//cout << ireal << " " << par_names[jpar] << " " << obs_names[iobs] << " " << cc << " " << bg_cc << endl;
			}

			//cout << par_names[jpar] << " " << obs_names[iobs] << " " << cc << " " << bg_cc << endl; 
			(cc < 0.0) ? sign = -1. : sign = 1.;

			bg_mean = bg_cc_vec.mean();
			bg_std = sqrt((bg_cc_vec - bg_mean).pow(2).sum() / (nreals - 1));
			thres = bg_mean + (sign * sigma_dist * bg_std);
			if (ies_verbose > 1)
			{
				
				while (true)
				{
					if (f_out_guard.try_lock())
					{
						*f_out << obs_names[iobs] << "," << par_names[jpar] << "," << cc << "," << bg_mean << "," << bg_std << "," << thres << "," << (((sign * cc) - (sign * thres)) > 0.0);
						for (int i = 0; i < nreals - 1; i++)
							*f_out << "," << bg_cc_vec[i];
						*f_out << endl;
						f_out_guard.unlock();
						break;

					}
				}
				
			}
			if (((sign * cc) - (sign * thres)) > 0.0)
			{
				//cout << par_names[jpar] << " " << obs_names[iobs] << " " << cc << " " << bg_mean << " " << bg_std << " " << thres << " kept " << endl;
				while (true)
				{
					if (triplets_guard.try_lock())
					{
						triplets.push_back(Eigen::Triplet<double>(iobs, jpar, cc));
						triplets_guard.unlock();
						break;
					}
				}
				no_obs = false;
				
			}
			

		}
		if (no_obs)
		{
			string pname;
			while (true)
			{
				if (par_names_guard.try_lock())
				{
					pname = par_names[jpar];
					par_names_guard.unlock();
					break;
				}
			}
			ss.str("");

			ss << "autoadaloc warning: parameter " << pname << " is completely localized -it maps to no observations";
			while (true)
			{
				if (pfm_guard.try_lock())
				{
					performance_log->log_event(ss.str());
					pfm_guard.unlock();
					break;
				}
			}
		}
	}
}


Eigen::MatrixXd Localizer::get_localizing_obs_hadamard_matrix(int num_reals, string col_name, vector<string> &obs_names)
{

	vector<double> values;
	vector<string> mat_cols = mat.get_col_names();
	vector<string>::iterator it = find(mat_cols.begin(), mat_cols.end(), col_name);
	if (it == mat_cols.end())
		throw runtime_error("Localizer::get_localizing_obs_hadamard_matrix() error: col_name not found in localizer matrix: " + col_name);
	int idx = it - mat_cols.begin();
	Eigen::VectorXd mat_vec = mat.e_ptr()->col(idx);
	int col_idx;
	Eigen::MatrixXd loc(obs_names.size(), num_reals);
	for (int i=0;i<obs_names.size();i++)
	{
		col_idx = obs2row_map[obs_names[i]];
		loc.row(i).setConstant(mat_vec[col_idx]);
	}
	return loc;

}


Eigen::MatrixXd Localizer::get_localizing_par_hadamard_matrix(int num_reals, string row_name, vector<string> &par_names)
{

	vector<double> values;
	vector<string> mat_rows = mat.get_row_names();
	vector<string>::iterator it = find(mat_rows.begin(), mat_rows.end(), row_name);
	if (it == mat_rows.end())
		throw runtime_error("Localizer::get_localizing_par_hadamard_matrix() error: row_name not found in localizer matrix: " + row_name);
	int idx = it - mat_rows.begin();
	Eigen::VectorXd mat_vec = mat.e_ptr()->row(idx);
	int col_idx;
	Eigen::MatrixXd loc(par_names.size(), num_reals);
	for (int i = 0; i<par_names.size(); i++)
	{
		col_idx = par2col_map[par_names[i]];
		loc.row(i).setConstant(mat_vec[col_idx]);
	}
	return loc;

}
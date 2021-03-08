#include <random>
#include <map>
#include <iomanip>
#include <mutex>
#include <thread>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "FileManager.h"
#include "Ensemble.h"
#include "RestartController.h"
#include "utilities.h"
#include "Ensemble.h"
#include "ObjectiveFunc.h"
#include "covariance.h"
#include "RedSVD-h.h"
#include "SVDPackage.h"
#include "eigen_tools.h"
#include "EnsembleMethodUtils.h"
#include "Localizer.h"
#include "Pest.h"
#include "PerformanceLog.h"
#include "RunStorage.h"
#include "RunManagerAbstract.h"




EnsembleSolver::EnsembleSolver(PerformanceLog* _performance_log, FileManager& _file_manager, Pest& _pest_scenario, ParameterEnsemble& _pe,
	ObservationEnsemble& _oe, ObservationEnsemble& _base_oe, Localizer& _localizer, Covariance& _parcov, Eigen::MatrixXd& _Am, L2PhiHandler& _ph,
	bool _use_localizer, int _iter, vector<string>& _act_par_names, vector<string>& _act_obs_names) :
	file_manager(_file_manager), pest_scenario(_pest_scenario), pe(_pe), oe(_oe), base_oe(_base_oe), localizer(_localizer), 
	parcov(_parcov), Am(_Am), ph(_ph), act_par_names(_act_par_names),act_obs_names(_act_obs_names)

{
	performance_log = _performance_log;
	use_localizer = _use_localizer;
	iter = _iter;
	verbose_level = pest_scenario.get_pestpp_options().get_ies_verbose_level();
	//prep the fast look par cov info
	message(0,"preparing fast-look containers for threaded localization solve");
	parcov_inv_map.clear();
	parcov_inv_map.reserve(pe.shape().second);
	Eigen::VectorXd parcov_inv;// = parcov.get(par_names).inv().e_ptr()->toDense().cwiseSqrt().asDiagonal();
	if (!parcov.isdiagonal())
	{
		parcov_inv = parcov.get_matrix().diagonal();
	}
	else
	{
		message(2,"extracting diagonal from prior parameter covariance matrix");
		Covariance parcov_diag;
		parcov_diag.from_diagonal(parcov);
		parcov_inv = parcov_diag.get_matrix().diagonal();
	}
	parcov_inv = parcov_inv.cwiseSqrt().cwiseInverse();
	vector<string> par_names = pe.get_var_names();
	for (int i = 0; i < parcov_inv.size(); i++)
		parcov_inv_map[par_names[i]] = parcov_inv[i];

	//prep the fast lookup weights info
	weight_map.clear();
	weight_map.reserve(oe.shape().second);
	vector<string> obs_names = oe.get_var_names();
	double w;
	for (int i = 0; i < obs_names.size(); i++)
	{
		w = pest_scenario.get_observation_info_ptr()->get_weight(obs_names[i]);
		//don't want to filter on weight here - might be changing weights, etc...
		//if (w == 0.0)
		//	continue;
		weight_map[obs_names[i]] = w;
	}

	//prep a fast look for obs, par and resid matrices - stored as column vectors in a map
	par_resid_map.clear();
	par_resid_map.reserve(par_names.size());
	par_diff_map.clear();
	par_diff_map.reserve(par_names.size());
	Am_map.clear();
	Am_map.reserve(par_names.size());
	obs_resid_map.clear();
	obs_resid_map.reserve(obs_names.size());
	obs_diff_map.clear();
	obs_diff_map.reserve(obs_names.size());
	obs_err_map.clear();
	obs_err_map.reserve(obs_names.size());

	//check for the 'center_on' real - it may have been dropped...
	string center_on = pest_scenario.get_pestpp_options().get_ies_center_on();
	if (center_on.size() > 0)
	{
		if (center_on != MEDIAN_CENTER_ON_NAME)
		{
			vector<string> real_names = oe.get_real_names();
			if (find(real_names.begin(), real_names.end(), center_on) == real_names.end())
			{
				message(1,"Warning: 'ies_center_on' real not found in obs en, reverting to mean...");
				center_on = "";
			}
			real_names = pe.get_real_names();
			if ((center_on.size() > 0) && find(real_names.begin(), real_names.end(), center_on) == real_names.end())
			{
				message(1,"Warning: 'ies_center_on' real not found in par en, reverting to mean...");
				center_on = "";
			}
		}
	}

	Eigen::MatrixXd mat = ph.get_obs_resid_subset(oe);
	for (int i = 0; i < obs_names.size(); i++)
	{
		obs_resid_map[obs_names[i]] = mat.col(i);
	}
	mat = oe.get_eigen_anomalies(center_on);
	for (int i = 0; i < obs_names.size(); i++)
	{
		obs_diff_map[obs_names[i]] = mat.col(i);
	}
	mat = base_oe.get_eigen(vector<string>(), obs_names);
	Observations ctl_obs = pest_scenario.get_ctl_observations();
	for (int i = 0; i < obs_names.size(); i++)
	{
		obs_err_map[obs_names[i]] = mat.col(i).array() - ctl_obs.get_rec(obs_names[i]);
	}
	mat = ph.get_par_resid_subset(pe);
	for (int i = 0; i < par_names.size(); i++)
	{
		par_resid_map[par_names[i]] = mat.col(i);
	}
	mat = pe.get_eigen_anomalies(center_on);
	for (int i = 0; i < par_names.size(); i++)
	{
		par_diff_map[par_names[i]] = mat.col(i);
	}
	if (!pest_scenario.get_pestpp_options().get_ies_use_approx())
	{
		for (int i = 0; i < par_names.size(); i++)
		{
			Am_map[par_names[i]] = Am.row(i);
		}
	}
	mat.resize(0, 0);
}

template<typename T, typename A>
void EnsembleSolver::message(int level, const string& _message, vector<T, A> _extras, bool echo)
{
	stringstream ss;
	if (level == 0)
		ss << endl << "  ---  ";
	else if (level == 1)
		ss << "...";
	ss << _message;
	if (_extras.size() > 0)
	{

		for (auto& e : _extras)
			ss << e << " , ";

	}
	if (level == 0)
		ss << "  ---  ";
	if ((echo) && ((verbose_level >= 2) || (level < 2)))
		cout << ss.str() << endl;
	file_manager.rec_ofstream() << ss.str() << endl;
	performance_log->log_event(_message);

}

void upgrade_thread_function(int id, int iter, double cur_lam, bool use_glm_form, vector<string> par_names,
								vector<string> obs_names, UpgradeThread& worker, exception_ptr& eptr)
{
	try
	{
		worker.work(id, iter, cur_lam, use_glm_form, par_names, obs_names);
	}
	catch (...)
	{
		eptr = current_exception();
	}

	return;
}


void EnsembleSolver::solve(int num_threads, double cur_lam, bool use_glm_form, ParameterEnsemble& pe_upgrade, unordered_map<string, pair<vector<string>, vector<string>>>& loc_map)
{
	pe_upgrade.set_zeros();
	stringstream ss;
	Localizer::How _how = localizer.get_how();
	Localizer::LocTyp loctyp = localizer.get_loctyp();
	bool use_cov_loc = true;
	if (loctyp == Localizer::LocTyp::LOCALANALYSIS)
		use_cov_loc = false;
	//LocalAnalysisUpgradeThread worker(performance_log, par_resid_map, par_diff_map, obs_resid_map, obs_diff_map,obs_err_map,
	//	localizer, parcov_inv_map, weight_map, pe_upgrade, loc_map, Am_map, _how);
	UpgradeThread* ut_ptr;
	if (!use_cov_loc)
		ut_ptr = new LocalAnalysisUpgradeThread(performance_log, par_resid_map, par_diff_map, obs_resid_map, obs_diff_map, obs_err_map,
			localizer, parcov_inv_map, weight_map, pe_upgrade, loc_map, Am_map, _how);
	else
		ut_ptr = new CovLocalizationUpgradeThread(performance_log, par_resid_map, par_diff_map, obs_resid_map, obs_diff_map, obs_err_map,
			localizer, parcov_inv_map, weight_map, pe_upgrade, loc_map, Am_map, _how);
	if ((num_threads < 1) || (loc_map.size() == 1))
		//if (num_threads < 1)
	{
		//worker.work(0, iter, cur_lam, use_glm_form, act_par_names, act_obs_names);
		ut_ptr->work(0, iter, cur_lam, use_glm_form, act_par_names, act_obs_names);
	}
	else
	{
		Eigen::setNbThreads(1);
		vector<thread> threads;
		vector<exception_ptr> exception_ptrs;
		message(2, "launching threads");

		for (int i = 0; i < num_threads; i++)
		{
			exception_ptrs.push_back(exception_ptr());
		}

		for (int i = 0; i < num_threads; i++)
		{
			//threads.push_back(thread(&LocalUpgradeThread::work, &worker, i, iter, cur_lam));

			//threads.push_back(thread(localanalysis_upgrade_thread_function, i, iter, cur_lam, use_glm_form, act_par_names, act_obs_names, std::ref(worker), std::ref(exception_ptrs[i])));
			threads.push_back(thread(upgrade_thread_function, i, iter, cur_lam, use_glm_form, act_par_names, 
				act_obs_names, std::ref(*ut_ptr), std::ref(exception_ptrs[i])));



		}
		message(2, "waiting to join threads");
		//for (auto &t : threads)
		//	t.join();
		ss.str("");
		int num_exp = 0;

		for (int i = 0; i < num_threads; ++i)
		{
			bool found = false;
			if (exception_ptrs[i])
			{
				found = true;
				num_exp++;
				try
				{
					rethrow_exception(exception_ptrs[i]);
				}
				catch (const std::exception& e)
				{
					//ss.str("");
					ss << " thread " << i << "raised an exception: " << e.what();
					//throw runtime_error(ss.str());
				}
				catch (...)
				{
					//ss.str("");
					ss << " thread " << i << "raised an exception";
					//throw runtime_error(ss.str());
				}
			}
			threads[i].join();
			if ((exception_ptrs[i]) && (!found))
			{
				num_exp++;
				try
				{
					rethrow_exception(exception_ptrs[i]);
				}
				catch (const std::exception& e)
				{
					//ss.str("");
					ss << " thread " << i << "raised an exception: " << e.what();
					//throw runtime_error(ss.str());
				}
				catch (...)
				{
					//ss.str("");
					ss << " thread " << i << "raised an exception: ";
					//throw runtime_error(ss.str());
				}
			}
		}
		if (num_exp > 0)
		{
			throw runtime_error(ss.str());
		}
		delete ut_ptr;
		message(2, "threaded localized upgrade calculation done");
	}
}


void EnsembleSolver::message(int level, const string& _message)
{
	message(level, _message, vector<string>());
}

template<typename T>
void EnsembleSolver::message(int level, const string& _message, T extra)
{
	stringstream ss;
	ss << _message << " " << extra;
	string s = ss.str();
	message(level, s);
}


UpgradeThread::UpgradeThread(PerformanceLog* _performance_log, unordered_map<string, 
	Eigen::VectorXd>& _par_resid_map, unordered_map<string, Eigen::VectorXd>& _par_diff_map, 
	unordered_map<string, Eigen::VectorXd>& _obs_resid_map, unordered_map<string, 
	Eigen::VectorXd>& _obs_diff_map, unordered_map<string, Eigen::VectorXd>& _obs_err_map, 
	Localizer& _localizer, unordered_map<string, double>& _parcov_inv_map, 
	unordered_map<string, double>& _weight_map, ParameterEnsemble& _pe_upgrade, 
	unordered_map<string, pair<vector<string>, vector<string>>>& _cases, 
	unordered_map<string, Eigen::VectorXd>& _Am_map, Localizer::How& _how):
	par_resid_map(_par_resid_map),par_diff_map(_par_diff_map), obs_resid_map(_obs_resid_map), 
	obs_diff_map(_obs_diff_map), obs_err_map(_obs_err_map), localizer(_localizer),
    pe_upgrade(_pe_upgrade), cases(_cases), parcov_inv_map(_parcov_inv_map), 
	weight_map(_weight_map), Am_map(_Am_map)
	{
		performance_log = _performance_log;
		how = _how;
		count = 0;

		for (auto& c : cases)
		{
			keys.push_back(c.first);
		}
		//sort(keys.begin(), keys.end());
		total = keys.size();
		//random_shuffle(keys.begin(), keys.end());

	
}



void CovLocalizationUpgradeThread::work(int thread_id, int iter, double cur_lam, bool use_glm_form, vector<string> par_names,
									vector<string> obs_names)
{

	//declare these helpers in here so they are thread safe...
		class local_utils
	{
	public:
		static Eigen::DiagonalMatrix<double, Eigen::Dynamic> get_matrix_from_map(vector<string>& names, unordered_map<string, double>& dmap)
		{
			Eigen::VectorXd vec(names.size());
			int i = 0;
			for (auto name : names)
			{
				vec[i] = dmap.at(name);
				++i;
			}
			Eigen::DiagonalMatrix<double, Eigen::Dynamic> m = vec.asDiagonal();
			return m;
		}
		static Eigen::MatrixXd get_matrix_from_map(int num_reals, vector<string>& names, unordered_map<string, Eigen::VectorXd>& emap)
		{
			Eigen::MatrixXd mat(num_reals, names.size());
			mat.setZero();

			for (int j = 0; j < names.size(); j++)
			{
				mat.col(j) = emap[names[j]];
			}

			return mat;
		}
		static void save_mat(int verbose_level, int tid, int iter, int t_count, string prefix, Eigen::MatrixXd& mat)
		{
			if (verbose_level < 2)
				return;

			if (verbose_level < 3)
				return;
			//cout << "thread: " << tid << ", " << t_count << ", " << prefix << " rows:cols" << mat.rows() << ":" << mat.cols() << endl;
			stringstream ss;

			ss << "thread_" << tid << ".count_ " << t_count << ".iter_" << iter << "." << prefix << ".dat";
			string fname = ss.str();
			ofstream f(fname);
			if (!f.good())
				cout << "error getting ofstream " << fname << endl;
			else
			{

				try
				{
					f << mat << endl;
					f.close();
				}
				catch (...)
				{
					cout << "error saving matrix " << fname << endl;
				}
			}
		}
	};

	stringstream ss;

	unique_lock<mutex> ctrl_guard(ctrl_lock, defer_lock);
	int maxsing, num_reals, verbose_level, pcount = 0, t_count=0;
	double eigthresh;
	bool use_approx;
	bool use_prior_scaling;
	bool use_localizer = false;
	bool loc_by_obs = true;

	while (true)
	{
		if (ctrl_guard.try_lock())
		{
			
			maxsing = pe_upgrade.get_pest_scenario_ptr()->get_svd_info().maxsing;
			eigthresh = pe_upgrade.get_pest_scenario_ptr()->get_svd_info().eigthresh;
			use_approx = pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_ies_use_approx();
			use_prior_scaling = pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_ies_use_prior_scaling();
			num_reals = pe_upgrade.shape().first;
			verbose_level = pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_ies_verbose_level();
			ctrl_guard.unlock();
			//if (pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_ies_localize_how()[0] == 'P')
			if (how == Localizer::How::PARAMETERS)
				loc_by_obs = false;
			else
			{
				throw runtime_error("Covariance localization only supporte for localization by parameters...");
			}
			break;
		}
	}
	ofstream f_thread;
	if (verbose_level > 2)
	{
		ss.str("");
		ss << "thread_" << thread_id << "part_map.csv";
		f_thread.open(ss.str());
		ss.str("");
	}
	Eigen::MatrixXd par_resid, par_diff, Am;
	Eigen::MatrixXd obs_resid, obs_diff, obs_err, loc;
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> weights, parcov_inv;
	vector<string> case_par_names, case_obs_names;
	string key;


	//these locks are used to control (thread-safe) access to the fast look up containers
	unique_lock<mutex> next_guard(next_lock, defer_lock);
	unique_lock<mutex> obs_diff_guard(obs_diff_lock, defer_lock);
	unique_lock<mutex> obs_resid_guard(obs_resid_lock, defer_lock);
	unique_lock<mutex> obs_err_guard(obs_err_lock, defer_lock);
	unique_lock<mutex> par_diff_guard(par_diff_lock, defer_lock);
	unique_lock<mutex> par_resid_guard(par_resid_lock, defer_lock);
	unique_lock<mutex> loc_guard(loc_lock, defer_lock);
	unique_lock<mutex> weight_guard(weight_lock, defer_lock);
	unique_lock<mutex> parcov_guard(parcov_lock, defer_lock);
	unique_lock<mutex> am_guard(am_lock, defer_lock);

	//reset all the solution parts
	par_resid.resize(0, 0);
	par_diff.resize(0, 0);
	obs_resid.resize(0, 0);
	obs_diff.resize(0, 0);
	obs_err.resize(0, 0);
	loc.resize(0, 0);
	Am.resize(0, 0);
	weights.resize(0);
	parcov_inv.resize(0);
	Am.resize(0, 0);


	//loop until this thread gets access to all the containers it needs to solve with 
	//which in this solution, is all active pars and obs...
	
	while (true)
	{
		//if all the solution pieces are filled, break out and solve!
		if (((use_approx) || (par_resid.rows() > 0)) &&
			(weights.size() > 0) &&
			(parcov_inv.size() > 0) &&
			(par_diff.rows() > 0) &&
			(obs_resid.rows() > 0) &&
			(obs_err.rows() > 0) &&
			(obs_diff.rows() > 0) &&
			((!use_localizer) || (loc.rows() > 0)) &&
			((use_approx) || (Am.rows() > 0)))
			break;


		//get access to the obs_diff container
		if ((obs_diff.rows() == 0) && (obs_diff_guard.try_lock()))
		{
			//piggy back here for thread safety
			//if (pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_svd_pack() == PestppOptions::SVD_PACK::PROPACK)
			//	use_propack = true;
			obs_diff = local_utils::get_matrix_from_map(num_reals, obs_names, obs_diff_map);
			obs_diff_guard.unlock();
		}

		//get access to the residual container
		if ((obs_resid.rows() == 0) && (obs_resid_guard.try_lock()))
		{
			obs_resid = local_utils::get_matrix_from_map(num_reals, obs_names, obs_resid_map);
			obs_resid_guard.unlock();
		}

		//get access to the obs noise container
		if ((obs_err.rows() == 0) && (obs_err_guard.try_lock()))
		{
			obs_err = local_utils::get_matrix_from_map(num_reals, obs_names, obs_err_map);
			obs_err_guard.unlock();
		}

		//get access to the par diff container
		if ((par_diff.rows() == 0) && (par_diff_guard.try_lock()))
		{
			par_diff = local_utils::get_matrix_from_map(num_reals, par_names, par_diff_map);
			par_diff_guard.unlock();
		}

		//get access to the par residual container
		if ((par_resid.rows() == 0) && (par_resid_guard.try_lock()))
		{
			par_resid = local_utils::get_matrix_from_map(num_reals, par_names, par_resid_map);
			par_resid_guard.unlock();
		}

		//get access to the obs weights container
		if ((weights.rows() == 0) && (weight_guard.try_lock()))
		{
			weights = local_utils::get_matrix_from_map(obs_names, weight_map);
			weight_guard.unlock();
		}

		//get access to the inverse prior parcov
		if ((parcov_inv.rows() == 0) && (parcov_guard.try_lock()))
		{
			parcov_inv = local_utils::get_matrix_from_map(par_names, parcov_inv_map);
			parcov_guard.unlock();
		}

		//if needed, get access to the Am container - needed in the full glm solution
		if ((!use_approx) && (Am.rows() == 0) && (am_guard.try_lock()))
		{
			//Am = local_utils::get_matrix_from_map(num_reals, par_names, Am_map).transpose();
			int am_cols = Am_map[par_names[0]].size();
			Am.resize(par_names.size(), am_cols);
			Am.setZero();

			for (int j = 0; j < par_names.size(); j++)
			{
				Am.row(j) = Am_map[par_names[j]];
			}
			am_guard.unlock();
		}
	}


	par_diff.transposeInPlace();
	obs_diff.transposeInPlace();
	obs_resid.transposeInPlace();
	par_resid.transposeInPlace();
	obs_err.transposeInPlace();


	//container to quickly look indices
	map<string, int> par2col_map;
	for (int i = 0; i < par_names.size(); i++)
		par2col_map[par_names[i]] = i;

	map<string, int> obs2row_map;
	for (int i = 0; i < obs_names.size(); i++)
		obs2row_map[obs_names[i]] = i;


	//form the scaled obs resid matrix
	local_utils::save_mat(verbose_level, thread_id, iter, t_count, "obs_resid", obs_resid);
	//Eigen::MatrixXd scaled_residual = weights * obs_resid;

	//form the (optionally) scaled par resid matrix
	local_utils::save_mat(verbose_level, thread_id, iter, t_count, "par_resid", par_resid);
	
	ss.str("");
	double scale = (1.0 / (sqrt(double(num_reals - 1))));

	local_utils::save_mat(verbose_level, thread_id, iter, t_count, "obs_diff", obs_diff);

	//calculate some full solution components first...

	Eigen::MatrixXd ivec, upgrade_1, s, s2, V, Ut, t;
	Eigen::MatrixXd upgrade_2;
	Eigen::VectorXd loc_vec;
	
	upgrade_2.resize(num_reals, pe_upgrade.shape().second);
	upgrade_2.setZero();

	if (!use_glm_form)
	{
		obs_diff = scale * obs_diff;// (H-Hm)/sqrt(N-1)		
		par_diff = scale * par_diff;// (K-Km)/sqrt(N-1)		
		obs_err = scale * obs_err; //  (E-Em)/sqrt(N-1)	

		SVD_REDSVD rsvd;
		Eigen::MatrixXd C = obs_diff + (cur_lam * obs_err); // curr_lam is the inflation factor 		
		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "C", C);
		rsvd.solve_ip(C, s, Ut, V, eigthresh, maxsing);
		Ut.transposeInPlace();
		V.resize(0, 0);
		C.resize(0, 0);
		obs_err.resize(0, 0);

		s2 = s.asDiagonal().inverse();
		for (int i = 0; i < s.size(); i++)
		{
			if (s(i) < 1e-50)
			{
				s2(i, i) = 0;
			}
		}

		t = obs_diff.transpose() * Ut.transpose() * s2 * Ut;
		Ut.resize(0, 0);
		obs_diff.resize(0,0);
	}


	//----------------------------------
	//glm solution
	//----------------------------------
	else
	{
		
		obs_diff = scale * (weights * obs_diff);
		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "par_diff", par_diff);
		if (use_prior_scaling)
			par_diff = scale * parcov_inv * par_diff;
		else
			par_diff = scale * par_diff;
		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "scaled_par_diff", par_diff);
		SVD_REDSVD rsvd;
		rsvd.solve_ip(obs_diff, s, Ut, V, eigthresh, maxsing);

		Ut.transposeInPlace();
		obs_diff.resize(0, 0);
		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "Ut", Ut);
		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "s", s);
		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "V", V);

		Eigen::MatrixXd s2 = s.cwiseProduct(s);

		ivec = ((Eigen::VectorXd::Ones(s2.size()) * (cur_lam + 1.0)) + s2).asDiagonal().inverse();
		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "ivec", ivec);

		obs_resid = weights * obs_resid;
		t = V * s.asDiagonal() * ivec * Ut;
		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "t", t);

		if ((!use_approx) && (iter > 1))
		{
			if (use_prior_scaling)
			{
				par_resid = parcov_inv * par_resid;
			}
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "Am", Am);
			Eigen::MatrixXd x4 = Am.transpose() * par_resid;
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X4", x4);

			par_resid.resize(0, 0);

			Eigen::MatrixXd x5 = Am * x4;
			x4.resize(0, 0);
			Am.resize(0, 0);

			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X5", x5);
			Eigen::MatrixXd x6 = par_diff.transpose() * x5;
			x5.resize(0, 0);

			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X6", x6);
			Eigen::MatrixXd x7 = V * ivec * V.transpose() * x6;
			x6.resize(0, 0);

			if (use_prior_scaling)
			{
				upgrade_2 = -1.0 * parcov_inv * par_diff * x7;
			}
			else
			{
				upgrade_2 = -1.0 * (par_diff * x7);
			}
			x7.resize(0, 0);

			//add upgrade_2 piece 
			//upgrade_1 = upgrade_1 + upgrade_2.transpose();
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "upgrade_2", upgrade_2);
			//upgrade_2.resize(0, 0);
		}
	}

	//This is the main thread loop - it continues until all upgrade pieces have been completed
	map<string, Eigen::VectorXd> loc_map;
	//solve for each case par_name
	Eigen::VectorXd pt = parcov_inv.diagonal();
	string name;
	while (true)
	{	
		//clear the pieces used last time...
		key = "";
		loc_map.clear();

		//get new pieces...
		while (true)
		{
			if ((key != "") && (loc_map.size() > 0))
				break;
			if (next_guard.try_lock())
			{
				//if all the pieces have been completed, return
				if (count == keys.size())
				{
					if (verbose_level > 1)
					{
						cout << "upgrade thread: " << thread_id << " processed " << pcount << " upgrade parts" << endl;
					}
					if (f_thread.good())
						f_thread.close();
					return;
				}
				key = keys[count];
				pair<vector<string>, vector<string>> p = cases.at(key);
				//we can potentially optimize the speed of this loc type by changing how many pars are
				//passed in each "case" so herein, we support a generic number of pars per case...
				case_par_names = p.second;
				//In this solution, we ignore case obs names since we are using the full set of obs for the solution...
				case_obs_names = p.first;
				
				if (count % 1000 == 0)
				{
					ss.str("");
					ss << "upgrade thread progress: " << count << " of " << total << " parts done";
					if (verbose_level > 1)
						cout << ss.str() << endl;
					performance_log->log_event(ss.str());
				}
				count++;
				t_count = count;
				pcount++;
				next_guard.unlock();
				
			}
			//get access to the localizer
			if ((key != "") && (loc_map.size() == 0) && (loc_guard.try_lock()))
			{
				//get the nobs-length localizing vector for each case par name
				for (auto& par_name : case_par_names)
				{
					loc_map[par_name] = localizer.get_obs_hadamard_vector(par_name, obs_names);
				}
				loc_guard.unlock();
			}
		}

		if (verbose_level > 2)
		{
			f_thread << t_count << "," << iter;
			for (auto name : case_par_names)
				f_thread << "," << name;
			for (auto name : case_obs_names)
				f_thread << "," << name;
			f_thread << endl;
		}
		upgrade_1.resize(num_reals,case_par_names.size());
		upgrade_1.setZero();

		for (int i = 0; i < case_par_names.size(); i++)
		{
			name = case_par_names[i];
			loc_vec = loc_map[name];
			Eigen::VectorXd par_vec = par_diff.row(par2col_map[name]) * t;
			//apply the localizer
			par_vec = par_vec.cwiseProduct(loc_vec);
			if (!use_glm_form)
				par_vec = -1.0 * par_vec.transpose() * obs_resid;
			else
			{
				par_vec = -1.0 * par_vec.transpose() * obs_resid;
				//if (use_prior_scaling)
				//	par_vec *= pt[i];
			}
			upgrade_1.col(i) += par_vec.transpose();
			//add the par change part for the full glm solution
			if ((!use_approx) && (iter > 1))
			{
				//update_2 is transposed relative to upgrade_1
				upgrade_1.col(i) += upgrade_2.row(par2col_map[name]);
			}

		}
		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "upgrade_1", upgrade_1);

		//put this piece of the upgrade vector in
		unique_lock<mutex> put_guard(put_lock, defer_lock);
		while (true)
		{
			if (put_guard.try_lock())
			{
				pe_upgrade.add_2_cols_ip(case_par_names, upgrade_1);
				put_guard.unlock();
				break;
			}
		}
	}
}


void LocalAnalysisUpgradeThread::work(int thread_id, int iter, double cur_lam, bool use_glm_form,
									vector<string> par_names, vector<string> obs_names)
{

	//declare these helpers in here so they are thread safe...
	class local_utils
	{
	public:
		static Eigen::DiagonalMatrix<double, Eigen::Dynamic> get_matrix_from_map(vector<string>& names, unordered_map<string, double>& dmap)
		{
			Eigen::VectorXd vec(names.size());
			int i = 0;
			for (auto name : names)
			{
				vec[i] = dmap.at(name);
				++i;
			}
			Eigen::DiagonalMatrix<double, Eigen::Dynamic> m = vec.asDiagonal();
			return m;
		}
		static Eigen::MatrixXd get_matrix_from_map(int num_reals, vector<string>& names, unordered_map<string, Eigen::VectorXd>& emap)
		{
			Eigen::MatrixXd mat(num_reals, names.size());
			mat.setZero();

			for (int j = 0; j < names.size(); j++)
			{
				mat.col(j) = emap[names[j]];
			}

			return mat;
		}
		static void save_mat(int verbose_level, int tid, int iter, int t_count, string prefix, Eigen::MatrixXd& mat)
		{
			if (verbose_level < 2)
				return;

			if (verbose_level < 3)
				return;
			//cout << "thread: " << tid << ", " << t_count << ", " << prefix << " rows:cols" << mat.rows() << ":" << mat.cols() << endl;
			stringstream ss;

			ss << "thread_" << tid << ".count_ " << t_count << ".iter_" << iter << "." << prefix << ".dat";
			string fname = ss.str();
			ofstream f(fname);
			if (!f.good())
				cout << "error getting ofstream " << fname << endl;
			else
			{

				try
				{
					f << mat << endl;
					f.close();
				}
				catch (...)
				{
					cout << "error saving matrix " << fname << endl;
				}
			}
		}
	};

	stringstream ss;

	unique_lock<mutex> ctrl_guard(ctrl_lock, defer_lock);
	int maxsing, num_reals, verbose_level, pcount = 0, t_count;
	double eigthresh;
	bool use_approx;
	bool use_prior_scaling;
	bool use_localizer = false;
	bool loc_by_obs = true;

	while (true)
	{
		if (ctrl_guard.try_lock())
		{
			maxsing = pe_upgrade.get_pest_scenario_ptr()->get_svd_info().maxsing;
			eigthresh = pe_upgrade.get_pest_scenario_ptr()->get_svd_info().eigthresh;
			use_approx = pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_ies_use_approx();
			use_prior_scaling = pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_ies_use_prior_scaling();
			num_reals = pe_upgrade.shape().first;
			verbose_level = pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_ies_verbose_level();
			ctrl_guard.unlock();
			//if (pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_ies_localize_how()[0] == 'P')
			if (how == Localizer::How::PARAMETERS)
				loc_by_obs = false;
			break;
		}
	}
	ofstream f_thread;
	if (verbose_level > 2)
	{
		ss.str("");
		ss << "thread_" << thread_id << "part_map.csv";
		f_thread.open(ss.str());
		ss.str("");
	}
	Eigen::MatrixXd par_resid, par_diff, Am;
	Eigen::MatrixXd obs_resid, obs_diff, obs_err, loc;
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> weights, parcov_inv;
	string key;


	//these locks are used to control (thread-safe) access to the fast look up containers
	unique_lock<mutex> next_guard(next_lock, defer_lock);
	unique_lock<mutex> obs_diff_guard(obs_diff_lock, defer_lock);
	unique_lock<mutex> obs_resid_guard(obs_resid_lock, defer_lock);
	unique_lock<mutex> obs_err_guard(obs_err_lock, defer_lock);
	unique_lock<mutex> par_diff_guard(par_diff_lock, defer_lock);
	unique_lock<mutex> par_resid_guard(par_resid_lock, defer_lock);
	unique_lock<mutex> loc_guard(loc_lock, defer_lock);
	unique_lock<mutex> weight_guard(weight_lock, defer_lock);
	unique_lock<mutex> parcov_guard(parcov_lock, defer_lock);
	unique_lock<mutex> am_guard(am_lock, defer_lock);
	unique_lock<mutex> put_guard(put_lock, defer_lock);

	//This is the main thread loop - it continues until all upgrade pieces have been completed
	while (true)
	{
		//First get a new case of par and obs names to solve with
		par_names.clear();
		obs_names.clear();
		use_localizer = false;
		
		while (true)
		{
			if (next_guard.try_lock())
			{
				//if all the pieces have been completed, return
				if (count == keys.size())
				{
					if (verbose_level > 1)
					{
						cout << "upgrade thread: " << thread_id << " processed " << pcount << " upgrade parts" << endl;
					}
					if (f_thread.good())
						f_thread.close();
					return;
				}
				key = keys[count];
				pair<vector<string>, vector<string>> p = cases.at(key);
				par_names = p.second;
				obs_names = p.first;
				if (localizer.get_use())
				{
					/*if ((loc_by_obs) && (par_names.size() == 1) && (k == par_names[0]))
						use_localizer = true;
					else if ((!loc_by_obs) && (obs_names.size() == 1) && (k == obs_names[0]))
						use_localizer = true;*/
					//if ((loc_by_obs) && (obs_names.size() == 1) && (k == obs_names[0]))	
					//else if ((!loc_by_obs) && (par_names.size() == 1) && (k == par_names[0]))
					//	use_localizer = true;
					use_localizer = true;
				}
				if (count % 1000 == 0)
				{
					ss.str("");
					ss << "upgrade thread progress: " << count << " of " << total << " parts done";
					if (verbose_level > 1)
						cout << ss.str() << endl;
					performance_log->log_event(ss.str());
				}
				count++;
				t_count = count;
				pcount++;
				next_guard.unlock();
				break;
			}
		}

		if (verbose_level > 2)
		{
			f_thread << t_count << "," << iter;
			for (auto name : par_names)
				f_thread << "," << name;
			for (auto name : obs_names)
				f_thread << "," << name;
			f_thread << endl;
		}

		//reset all the solution parts
		par_resid.resize(0, 0);
		par_diff.resize(0, 0);
		obs_resid.resize(0, 0);
		obs_diff.resize(0, 0);
		obs_err.resize(0, 0);
		loc.resize(0, 0);
		Am.resize(0, 0);
		weights.resize(0);
		parcov_inv.resize(0);
		Am.resize(0, 0);

		
		//now loop until this thread gets access to all the containers it needs to solve with
		while (true)
		{
			//if all the solution pieces are filled, break out and solve!
			if (((use_approx) || (par_resid.rows() > 0)) &&
				(weights.size() > 0) &&
				(parcov_inv.size() > 0) &&
				(par_diff.rows() > 0) &&
				(obs_resid.rows() > 0) &&
				(obs_err.rows() > 0) &&
				(obs_diff.rows() > 0) &&
				((!use_localizer) || (loc.rows() > 0)) &&
				((use_approx) || (Am.rows() > 0)))
				break;
			
			//get access to the localizer
			if ((use_localizer) && (loc.rows() == 0) && (loc_guard.try_lock()))
			{
				//get a matrix that is either the shape of par diff or obs diff
				if (loc_by_obs)
				{
					//loc = localizer.get_localizing_par_hadamard_matrix(num_reals, obs_names[0], par_names);
					loc = localizer.get_pardiff_hadamard_matrix(num_reals, key, par_names);
				}

				else
				{
					//loc = localizer.get_localizing_obs_hadamard_matrix(num_reals, par_names[0], obs_names);
					loc = localizer.get_obsdiff_hadamard_matrix(num_reals, key, obs_names);
				}
				loc_guard.unlock();
			}

			//get access to the obs_diff container
			if ((obs_diff.rows() == 0) && (obs_diff_guard.try_lock()))
			{
				//piggy back here for thread safety
				//if (pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_svd_pack() == PestppOptions::SVD_PACK::PROPACK)
				//	use_propack = true;
				obs_diff = local_utils::get_matrix_from_map(num_reals, obs_names, obs_diff_map);
				obs_diff_guard.unlock();
			}

			//get access to the residual container
			if ((obs_resid.rows() == 0) && (obs_resid_guard.try_lock()))
			{
				obs_resid = local_utils::get_matrix_from_map(num_reals, obs_names, obs_resid_map);
				obs_resid_guard.unlock();
			}

			//get access to the obs noise container
			if ((obs_err.rows() == 0) && (obs_err_guard.try_lock()))
			{
				obs_err = local_utils::get_matrix_from_map(num_reals, obs_names, obs_err_map);
				obs_err_guard.unlock();
			}

			//get access to the par diff container
			if ((par_diff.rows() == 0) && (par_diff_guard.try_lock()))
			{
				par_diff = local_utils::get_matrix_from_map(num_reals, par_names, par_diff_map);
				par_diff_guard.unlock();
			}

			//get access to the par residual container
			if ((par_resid.rows() == 0) && (par_resid_guard.try_lock()))
			{
				par_resid = local_utils::get_matrix_from_map(num_reals, par_names, par_resid_map);
				par_resid_guard.unlock();
			}

			//get access to the obs weights container
			if ((weights.rows() == 0) && (weight_guard.try_lock()))
			{
				weights = local_utils::get_matrix_from_map(obs_names, weight_map);
				weight_guard.unlock();
			}

			//get access to the inverse prior parcov
			if ((parcov_inv.rows() == 0) && (parcov_guard.try_lock()))
			{
				parcov_inv = local_utils::get_matrix_from_map(par_names, parcov_inv_map);
				parcov_guard.unlock();
			}

			//if needed, get access to the Am container - needed in the full glm solution
			if ((!use_approx) && (Am.rows() == 0) && (am_guard.try_lock()))
			{
				//Am = local_utils::get_matrix_from_map(num_reals, par_names, Am_map).transpose();
				int am_cols = Am_map[par_names[0]].size();
				Am.resize(par_names.size(), am_cols);
				Am.setZero();

				for (int j = 0; j < par_names.size(); j++)
				{
					Am.row(j) = Am_map[par_names[j]];
				}
				am_guard.unlock();
			}
		}
		

		par_diff.transposeInPlace();
		obs_diff.transposeInPlace();
		obs_resid.transposeInPlace();
		par_resid.transposeInPlace();
		obs_err.transposeInPlace();

		//form the scaled obs resid matrix
		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "obs_resid", obs_resid);
		//Eigen::MatrixXd scaled_residual = weights * obs_resid;

		//form the (optionally) scaled par resid matrix
		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "par_resid", par_resid);
		

		stringstream ss;

		double scale = (1.0 / (sqrt(double(num_reals - 1))));

		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "obs_diff", obs_diff);

		//apply the localizer here...
		if (use_localizer)
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "loc", loc);
		if (use_localizer)
		{
			if (loc_by_obs)
				par_diff = par_diff.cwiseProduct(loc);
			else
				obs_diff = obs_diff.cwiseProduct(loc);
		}

		Eigen::MatrixXd ivec, upgrade_1, s, s2, V, Ut;

		//----------------------------------
		//es-mda solution
		//----------------------------------
		if (!use_glm_form)
		{
			obs_diff = scale * obs_diff;// (H-Hm)/sqrt(N-1)		
			par_diff = scale * par_diff;// (K-Km)/sqrt(N-1)		
			obs_err = scale * obs_err; //  (E-Em)/sqrt(N-1)	
			
			SVD_REDSVD rsvd;
			Eigen::MatrixXd C = obs_diff + (cur_lam * obs_err); // curr_lam is the inflation factor 		
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "C", C);
			rsvd.solve_ip(C, s, Ut, V, eigthresh, maxsing);
			Ut.transposeInPlace();
			V.resize(0, 0);
			C.resize(0, 0);
			obs_err.resize(0, 0);

			s2 = s.asDiagonal().inverse();
			for (int i = 0; i < s.size(); i++)
			{
				if (s(i) < 1e-50)
				{
					s2(i, i) = 0;
				}
			}
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "s2", s2);
			Eigen::MatrixXd X1 = s2 * Ut;
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X1", X1);
			X1 = X1 * obs_resid;
			obs_resid.resize(0, 0);
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X1_obs_resid", X1);
			X1 = Ut.transpose() * X1;
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X1_Ut", X1);
			Ut.resize(0, 0);
			X1 = obs_diff.transpose() * X1;
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X1_obs_diff", X1);
			obs_diff.resize(0, 0);
			upgrade_1 = -1.0 * par_diff * X1;
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "upgrade_1", upgrade_1);
			upgrade_1.transposeInPlace();
			
		}


		//----------------------------------
		//glm solution
		//----------------------------------
		else
		{
			obs_resid = weights * obs_resid;
			obs_diff = scale * (weights * obs_diff);
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "par_diff", par_diff);
			if (use_prior_scaling)
				par_diff = scale * parcov_inv * par_diff;
			else
				par_diff = scale * par_diff;
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "scaled_par_diff", par_diff);
			SVD_REDSVD rsvd;
			rsvd.solve_ip(obs_diff, s, Ut, V, eigthresh, maxsing);

			Ut.transposeInPlace();
			obs_diff.resize(0, 0);
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "Ut", Ut);
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "s", s);
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "V", V);

			Eigen::MatrixXd s2 = s.cwiseProduct(s);

			ivec = ((Eigen::VectorXd::Ones(s2.size()) * (cur_lam + 1.0)) + s2).asDiagonal().inverse();
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "ivec", ivec);

			Eigen::MatrixXd X1 = Ut * obs_resid;
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X1", X1);
			
			Eigen::MatrixXd X2 = ivec * X1;
			X1.resize(0, 0);
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X2", X2);
			
			Eigen::MatrixXd X3 = V * s.asDiagonal() * X2;
			X2.resize(0, 0);
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X3", X3);
			upgrade_1 = -1.0 * par_diff * X3;
			
			if (use_prior_scaling)
			{
				//upgrade_1 = parcov_inv * upgrade_1;
			}
			
			upgrade_1.transposeInPlace();
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "upgrade_1", upgrade_1);
			X3.resize(0, 0);

			Eigen::MatrixXd upgrade_2;
			if ((!use_approx) && (iter > 1))
			{
				if (use_prior_scaling)
				{
					par_resid = parcov_inv * par_resid;
				}
			
				local_utils::save_mat(verbose_level, thread_id, iter, t_count, "Am", Am);
				Eigen::MatrixXd x4 = Am.transpose() * par_resid;
				local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X4", x4);

				par_resid.resize(0, 0);

				Eigen::MatrixXd x5 = Am * x4;
				x4.resize(0, 0);
				Am.resize(0, 0);

				local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X5", x5);
				Eigen::MatrixXd x6 = par_diff.transpose() * x5;
				x5.resize(0, 0);

				local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X6", x6);
				Eigen::MatrixXd x7 = V * ivec * V.transpose() * x6;
				x6.resize(0, 0);

				if (use_prior_scaling)
				{
					upgrade_2 = -1.0 * parcov_inv * par_diff * x7;
				}
				else
				{
					upgrade_2 = -1.0 * (par_diff * x7);
				}
				x7.resize(0, 0);

				upgrade_1 = upgrade_1 + upgrade_2.transpose();
				local_utils::save_mat(verbose_level, thread_id, iter, t_count, "upgrade_2", upgrade_2);
				upgrade_2.resize(0, 0);

			}
		}
		
		while (true)
		{
			if (put_guard.try_lock())
			{
				pe_upgrade.add_2_cols_ip(par_names, upgrade_1);
				put_guard.unlock();
				break;
			}
		}
	}

}





L2PhiHandler::L2PhiHandler(Pest *_pest_scenario, FileManager *_file_manager,
	ObservationEnsemble *_oe_base, ParameterEnsemble *_pe_base,
	Covariance *_parcov, bool should_prep_csv, string _tag)
{
	pest_scenario = _pest_scenario;
	file_manager = _file_manager;
	oe_base = _oe_base;
	pe_base = _pe_base;
	tag = _tag;
	
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
		prepare_csv(file_manager->open_ofile_ext(tag+"phi.actual.csv"), oreal_names);
		prepare_csv(file_manager->open_ofile_ext(tag+"phi.meas.csv"), oreal_names);
		prepare_csv(file_manager->open_ofile_ext(tag+"phi.composite.csv"), oreal_names);
		prepare_csv(file_manager->open_ofile_ext(tag+"phi.regul.csv"), preal_names);
		prepare_group_csv(file_manager->open_ofile_ext(tag+"phi.group.csv"));
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
	if (pest_scenario->get_pestpp_options().get_save_binary())
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
	if (pest_scenario->get_pestpp_options().get_save_binary())
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

map<string, double>* L2PhiHandler::get_phi_map_ptr(L2PhiHandler::phiType pt)
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

map<string, double> L2PhiHandler::get_phi_map(L2PhiHandler::phiType pt)
{
	switch (pt)
	{
	case L2PhiHandler::phiType::ACTUAL:
		return actual;
	case L2PhiHandler::phiType::COMPOSITE:
		return composite;
	case L2PhiHandler::phiType::MEAS:
		return meas;
	case L2PhiHandler::phiType::REGUL:
		return regul;
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
	map<string, double>* phi_map = get_phi_map_ptr(pt);
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
	map<string, double>* phi_map = get_phi_map_ptr(pt);
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
	map<string, double>* phi_map = get_phi_map_ptr(pt);
	map<string, double>::iterator pi = phi_map->begin(), end = phi_map->end();
	for (; pi != end; ++pi)
		mx = (pi->second > mx) ? pi->second : mx;
	return mx;
}

double L2PhiHandler::get_min(phiType pt)
{
	double mn = 1.0e+30;
	map<string, double>* phi_map = get_phi_map_ptr(pt);
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
	write_csv(iter_num, total_runs, file_manager->get_ofstream(tag+"phi.actual.csv"), phiType::ACTUAL,oreal_names);
	write_csv(iter_num, total_runs, file_manager->get_ofstream(tag+"phi.meas.csv"), phiType::MEAS, oreal_names);
	if (pest_scenario->get_pestpp_options().get_ies_reg_factor() != 0.0)
	{
		write_csv(iter_num, total_runs, file_manager->get_ofstream(tag+"phi.regul.csv"), phiType::REGUL, preal_names);	
	}
	write_csv(iter_num, total_runs, file_manager->get_ofstream(tag+"phi.composite.csv"), phiType::COMPOSITE, oreal_names);
	if (write_group)
		write_group_csv(iter_num, total_runs, file_manager->get_ofstream(tag+"phi.group.csv"));
}

void L2PhiHandler::write_group(int iter_num, int total_runs, vector<double> extra)
{
	write_group_csv(iter_num, total_runs, file_manager->get_ofstream(tag+"phi.group.csv"),extra);
}

void L2PhiHandler::write_csv(int iter_num, int total_runs, ofstream &csv, phiType pt, vector<string> &names)
{
	map<string, double>* phi_map = get_phi_map_ptr(pt);
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
	for (int i = 0; i < names.size(); i++)
	{
		
		if ((_meas[names[i]] > bad_phi) || (_meas[names[i]] > mean + (std * bad_phi_sigma)))
		{
			if (names[i] == BASE_REAL_NAME)
				cout << "...not dropping 'base' real even though phi is 'bad'" << endl;
			else
				idxs.push_back(i);
		}	
	}
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

	Eigen::MatrixXd resid = get_obs_resid(oe);// this is (Ho - Hs)
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
	map<string, double> mean_map, std_map;
	base_pe_ptr->fill_moment_maps(mean_map, std_map);
	init_moments = pair<map<string, double>, map<string, double>>(mean_map, std_map);
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
	int mxlen = 0;
	for (auto& g : grp_names)
		mxlen = max(mxlen, (int)g.size());

	stringstream ss;
	ofstream &frec = file_manager_ptr->rec_ofstream();
	ss << endl << "   ---  Parameter Group Change Summmary  ---    " << endl;
	ss << "   (compared to the initial ensemble using active realizations)" << endl;
	cout << ss.str();
	frec << ss.str();
	ss.str("");
	ss << setw(mxlen) << "group" << setw(12) << "mean change" << setw(12) << "std change" << setw(18) << "num at/near bnds" << setw(16) << "% at/near bnds" << endl;
	cout << ss.str();
	frec << ss.str();	
	
	
	int i = 0;
	for (auto &grp_name : grp_names)
	{
		double mean_diff = mean_change[grp_name];
		double std_diff = std_change[grp_name];
		int num_out = num_at_bounds[grp_name];
		int percent_out = percent_at_bounds[grp_name];
		ss.str("");
		ss << setw(mxlen) << pest_utils::lower_cp(grp_name) << setw(12) << mean_diff * 100.0 << setw(12) << std_diff * 100.0 << setw(18);
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
		f << pest_utils::lower_cp(grp_name) << "," << mean_change[grp_name]*100.0 << "," << std_change[grp_name]*100.0 << ",";
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
	//pair<map<string, double>, map<string, double>> moments = pe.get_moment_maps();
	//init_moments = base_pe_ptr->get_moment_maps(pe.get_real_names());
	map<string, double> mean_map, std_map;
	base_pe_ptr->fill_moment_maps(mean_map, std_map);
	init_moments = pair<map<string, double>, map<string, double>>(mean_map, std_map);
	mean_map.clear();
	std_map.clear();
	pe.fill_moment_maps(mean_map, std_map);
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
			value2 = value1 - mean_map[par_name];
			if ((value1 != 0.0) && (value2 != 0.0))
				mean_diff += value2 / value1;
			value1 = init_moments.second[par_name];
			value2 = value1 - std_map[par_name];
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


pair<Parameters, Observations> save_real_par_rei(Pest& pest_scenario, ParameterEnsemble& pe, ObservationEnsemble& oe,
	OutputFileWriter& output_file_writer, FileManager& file_manager, int iter,string tag)
{
	stringstream ss;
	map<string, int> vmap = pe.get_real_map();
	Parameters pars;
	Observations obs;
	if (vmap.find(tag) != vmap.end())
	{
		ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();
		Parameters pars;
		pars.update(pe.get_var_names(), eigenvec_2_stlvec(pe.get_real_vector(tag)));
		if (pe.get_trans_status() == ParameterEnsemble::transStatus::NUM)
			pts.numeric2ctl_ip(pars);
		// save parameters to .par file
		if (iter >= 0)
			ss << iter << ".";
		ss << pest_utils::lower_cp(tag) << ".par";
		output_file_writer.write_par(file_manager.open_ofile_ext(ss.str()), pars, *(pts.get_offset_ptr()),
			*(pts.get_scale_ptr()));
		file_manager.close_file("par");

		vmap = oe.get_real_map();
		if (vmap.find(tag) == vmap.end())
		{
			//message(2, "unable to find 'BASE' realization in obs ensemble for saving .base.rei file, continuing...");
		}
		else
		{
			Observations obs;
			obs.update(oe.get_var_names(), eigenvec_2_stlvec(oe.get_real_vector(tag)));
			ObjectiveFunc obj_func(&(pest_scenario.get_ctl_observations()), &(pest_scenario.get_ctl_observation_info()), &(pest_scenario.get_prior_info()));
			// save new residuals to .rei file
			ss.str("");
			if (iter >= 0)
				ss << iter << ".";
			ss << pest_utils::lower_cp(tag) <<  ".rei";
			output_file_writer.write_rei(file_manager.open_ofile_ext(ss.str()), iter,
				pest_scenario.get_ctl_observations(), obs, obj_func, pars);
		}
		cout << "saved par and rei files for realization " << tag << " for iteration " << iter << endl;
	}
	
	return pair<Parameters, Observations>(pars, obs);

}


vector<int> run_ensemble_util(PerformanceLog* performance_log, ofstream& frec,ParameterEnsemble& _pe, ObservationEnsemble& _oe, 
	RunManagerAbstract* run_mgr_ptr, bool check_pe_consistency, const vector<int>& real_idxs, int da_cycle)
{
	stringstream ss;
	ss << "queuing " << _pe.shape().first << " runs";
	performance_log->log_event(ss.str());
	run_mgr_ptr->reinitialize();
	map<int, int> real_run_ids;
	try
	{
		real_run_ids = _pe.add_runs(run_mgr_ptr, real_idxs,da_cycle);
	}
	catch (const exception& e)
	{
		ss.str("");
		ss << "run_ensemble() error queueing runs: " << e.what();
		throw runtime_error(ss.str());
	}
	catch (...)
	{
		throw runtime_error(string("run_ensemble() error queueing runs"));
	}
	performance_log->log_event("making runs");
	try
	{
		run_mgr_ptr->run();
	}
	catch (const exception& e)
	{
		ss.str("");
		ss << "error running ensemble: " << e.what();
		throw runtime_error(ss.str());
	}
	catch (...)
	{
		throw runtime_error(string("error running ensemble"));
	}

	performance_log->log_event("processing runs");
	if (real_idxs.size() > 0)
	{
		_oe.keep_rows(real_idxs);
	}
	vector<int> failed_real_indices;
	ParameterEnsemble run_mgr_pe;
	try
	{
		failed_real_indices = _oe.update_from_runs(real_run_ids, run_mgr_ptr, run_mgr_pe);
	} 
	catch (const exception& e)
	{
		ss.str("");
		ss << "error processing runs: " << e.what();
		throw runtime_error(ss.str());
	}
	catch (...)
	{
		throw runtime_error(string("error processing runs"));
	}
	
	if (failed_real_indices.size() > 0)
	{
		ss.str("");
		vector<string> par_real_names = _pe.get_real_names();
		vector<string> obs_real_names = _oe.get_real_names();
		ss << "the following par:obs realization runs failed: ";
		for (auto& i : failed_real_indices)
		{
			ss << par_real_names[i] << ":" << obs_real_names[i] << ',';
		}
		performance_log->log_event(ss.str());
		performance_log->log_event("dropping failed realizations");
		_pe.drop_rows(failed_real_indices);
		_oe.drop_rows(failed_real_indices);
	}
	if (check_pe_consistency)
	{
		if (_pe.shape().first != run_mgr_pe.shape().first)
		{
			ss.str("");
			ss << "error checking pe consistency: run_mgr_pe has different number of rows than _pe: " << run_mgr_pe.shape().first << " vs " << _pe.shape().first;
			throw runtime_error(ss.str());
		}
		run_mgr_pe.transform_ip(ParameterEnsemble::transStatus::NUM);
		if (_pe.shape().second != run_mgr_pe.shape().second)
		{
			ss.str("");
			ss << "error checking pe consistency: run_mgr_pe has different number of columns than _pe: " << run_mgr_pe.shape().second << " vs " << _pe.shape().second;
			throw runtime_error(ss.str());
		}

		
		Eigen::VectorXd org_real, new_real, diff;
		double percent_diff_max;
		int max_idx;
		vector<string> real_names = _pe.get_real_names(), par_names = _pe.get_var_names();
		ss.str("");
		frec << "checking for consistency in parameter ensemble" << endl;
		frec << "   full listing of consistency check results:" << endl;
		vector<string> failing_reals;
		for (auto rri : real_run_ids)
		{
			org_real = _pe.get_real_vector(rri.first);
			new_real = run_mgr_pe.get_real_vector(rri.first);
			diff = (org_real.array() - new_real.array()).cwiseAbs();
			diff = 100.0 * diff.cwiseQuotient(org_real);
			percent_diff_max = diff.maxCoeff(&max_idx);
			frec << "   realization: " << real_names[rri.first] << ", run_id: " << rri.second << ", max % diff: " << percent_diff_max << " at parameter: " << par_names[max_idx] << endl;
			if (percent_diff_max > 2.0)
			{
				failing_reals.push_back(real_names[rri.first]);
			}
		}
		if (failing_reals.size() > 0)
		{
			cout << "parameter ensemble consisteny check failed for " << failing_reals.size() << ", see .rec file for listing" << endl;
			frec << "ERROR: the following realizations failed consistency check:" << endl;
			for (auto fr : failing_reals)
				frec << fr << ",";
			frec << endl;
			throw runtime_error("parameter ensemble consistency check failed, see .rec file");
			
			

		}
	}
	return failed_real_indices;
}


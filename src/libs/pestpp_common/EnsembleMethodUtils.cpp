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
#include "utilities.h"




EnsembleSolver::EnsembleSolver(PerformanceLog* _performance_log, FileManager& _file_manager, Pest& _pest_scenario, ParameterEnsemble& _pe,
	ObservationEnsemble& _oe, ObservationEnsemble& _base_oe, ObservationEnsemble& _weights, Localizer& _localizer, Covariance& _parcov, Eigen::MatrixXd& _Am, L2PhiHandler& _ph,
	bool _use_localizer, int _iter, vector<string>& _act_par_names, vector<string>& _act_obs_names, double _reg_factor) :
	file_manager(_file_manager), pest_scenario(_pest_scenario), pe(_pe), oe(_oe), base_oe(_base_oe), weights(_weights), localizer(_localizer),
	parcov(_parcov), Am(_Am), ph(_ph), act_par_names(_act_par_names),act_obs_names(_act_obs_names), reg_factor(_reg_factor) {
    performance_log = _performance_log;
    use_localizer = _use_localizer;
    iter = _iter;
    verbose_level = pest_scenario.get_pestpp_options().get_ies_verbose_level();
    //prep the fast look par cov info

    if (localizer.get_use()) {
        message(1, "preparing fast-lookup containers for threaded localization solve");
        initialize_for_localized_solve();
    }

}

MmNeighborThread::MmNeighborThread(unordered_map<string,Eigen::VectorXd>& _real_vec_map,
                                   unordered_map<string,vector<int>>& _mm_real_idx_map,
                                   unordered_map<string,pair<vector<string>,vector<string>>>& _mm_real_name_map,
                                   unordered_map<string,unordered_map<string,double>>& _neighbor_phi_map,
                                   unordered_map<string,unordered_map<string,double>>& _neighbor_pardist_map):
real_vec_map(_real_vec_map),mm_real_idx_map(_mm_real_idx_map),mm_real_name_map(_mm_real_name_map),
neighbor_phi_map(_neighbor_phi_map),neighbor_pardist_map(_neighbor_pardist_map)
{
    count = 0;
    total = 0;
    for (int i=0;i<real_vec_map.size();i++)
        indexes.push_back(i);
}

void mm_neighbor_thread_function(int id, int verbose_level, double mm_alpha, map<string,map<string,double>> weight_phi_map,
                                 vector<string> preal_names, vector<string> oreal_names,Eigen::SparseMatrix<double> parcov_inv,
                                 map<string,int> real_map,
                                 MmNeighborThread& worker, exception_ptr& eptr) {
    try {
        worker.work(id,verbose_level,mm_alpha,weight_phi_map,preal_names,oreal_names,real_map,parcov_inv);
    }
    catch (...) {
        eptr = current_exception();
    }
    return;
}

void MmNeighborThread::work(int tid, int verbose_level, double mm_alpha, map<string,map<string,double>> weight_phi_map, vector<string> preal_names,
                            vector<string> oreal_names,map<string,int> real_map,Eigen::SparseMatrix<double> parcov_inv)
{
    unique_lock<mutex> next_guard(next_lock, defer_lock);
    unique_lock<mutex> results_guard(results_lock, defer_lock);


    stringstream ss;
    //util functions
    typedef std::function<bool(std::pair<std::string, double>, std::pair<std::string, double>)> Comparator;
    // Defining a lambda function to compare two pairs. It will compare two pairs using second field
    Comparator compFunctor = [](std::pair<std::string, double> elem1, std::pair<std::string, double> elem2) {
        return elem1.second < elem2.second;
    };
    typedef std::set<std::pair<std::string, double>, Comparator> sortedset;
    typedef std::set<std::pair<std::string, double>, Comparator>::iterator sortedset_iter;

    //need to guard these
    int num_reals,idx;
    //map<string,int> real_map;
    num_reals = real_vec_map.size();

    int subset_size = (int) (((double) num_reals) * mm_alpha);
    if (subset_size <= 0)
        throw runtime_error("solve_multimodal error: subsize size zero or negative");

    string real_name, oreal_name;

    map<string, double> euclid_par_dist;
    Eigen::VectorXd real, diff;
    double edist;

    vector<int> real_idxs;
    vector<string> pe_real_names_case,oe_real_names_case;
    string prname, orname;
    map<string, double> phi_map;
    Eigen::MatrixXd dmat;
    unordered_map<string,double> nphi,npardist;
    while (true) {

        while (true) {
            if (next_guard.try_lock()) {
                //if all the pieces have been completed, return
                if (count == indexes.size()) {

                    return;
                }

                idx = count;
                count++;
                next_guard.unlock();
                break;
            }
        }

        real_name = preal_names[idx];
        oreal_name = oreal_names[idx];
        euclid_par_dist.clear();
        phi_map = weight_phi_map.at(oreal_name);

        real = real_vec_map.at(real_name);


        //scale the phi values from 0 to 1

        double mx = -1.0e+300;
        for (auto &p : phi_map) {
            if (p.second > mx)
                mx = p.second;
        }
        for (auto &p : phi_map) {
            p.second /= mx;
        }

        //flip to map to par realization names
        map<string, double> par_phi_map;
        for (int ii = 0; ii < oreal_names.size(); ii++) {
            par_phi_map[preal_names[ii]] = phi_map.at(oreal_names[ii]);
        }
        phi_map.clear();

        for (int ii=0;ii<num_reals;ii++)
        {
            if (preal_names[ii] == real_name)
                continue;
            //diff = dmat.row(ii);
            //diff = pe.get_eigen_ptr()->row(ii) - real;
            diff = real_vec_map[preal_names[ii]] - real;
            edist = diff.transpose() * parcov_inv * diff;
            euclid_par_dist[preal_names[ii]] = edist;
        }

        //scale the par space distance from 0 to 1
        mx = -1e300;
        for (auto &e : euclid_par_dist) {
            if (e.second > mx)
                mx = e.second;
        }
        if (mx <= 0.0) {
            ss.str("");
            ss << "multimodal solve error: maximum par diff for realization '" << real_name << "' not valid";
            throw runtime_error(ss.str());
        } else {
            for (auto &e : euclid_par_dist) {
                e.second /= mx;
            }
        }

        //calc composite score
        map<string, double> composite_score;
        for (auto &rname : preal_names) {
            if (rname == real_name)
                continue;
            composite_score[rname] = euclid_par_dist.at(rname) + par_phi_map.at(rname);
        }
        //sort by score
        sortedset fitness_sorted(composite_score.begin(), composite_score.end(), compFunctor);

        //fill the containers
        real_idxs.clear();
        oe_real_names_case.clear();
        pe_real_names_case.clear();
        real_idxs.push_back(real_map.at(real_name));
        pe_real_names_case.push_back(real_name);
        oe_real_names_case.push_back(oreal_name);
        nphi.clear();
        npardist.clear();

        int iii = 0;
        for (sortedset_iter ii = fitness_sorted.begin(); ii != fitness_sorted.end(); ++ii) {
            int iidx = real_map.at(ii->first);
            real_idxs.push_back(iidx);
            pe_real_names_case.push_back(ii->first);
            oe_real_names_case.push_back(oreal_names[iidx]);
            nphi[preal_names[iidx]] = par_phi_map.at(ii->first);
            npardist[preal_names[iidx]] = euclid_par_dist.at(ii->first);
            iii++;
            if (iii >= subset_size)//plus one to count 'real_name'
                break;

        }
        if (real_idxs.size() != subset_size+1)
        {
            throw runtime_error("multimodal_solve: real_idxs.size() != subset_size");
        }

        //need to guard these
        while (true)
        {
            if (results_lock.try_lock())
            {
                mm_real_idx_map[real_name] = real_idxs;
                mm_real_name_map[real_name] = make_pair(pe_real_names_case,oe_real_names_case);
                neighbor_phi_map[real_name] = nphi;
                neighbor_pardist_map[real_name] = npardist;

                results_lock.unlock();
                break;
            }
        }
    }
}


void EnsembleSolver::update_multimodal_components(const double mm_alpha) {
    mm_real_idx_map.clear();
    mm_q_vec_map.clear();
    mm_real_name_map.clear();
    stringstream ss;
    int num_threads = pest_scenario.get_pestpp_options().get_ies_num_threads();
    Eigen::SparseMatrix<double> parcov_inv = parcov.inv().get_matrix();
    Eigen::MatrixXd wmat = weights.get_eigen(vector<string>(),act_obs_names);
    performance_log->log_event("getting phi vectors for all weights");
    map<string,map<string,double>> weight_phi_map = ph.get_meas_phi_weight_ensemble(oe,weights);
    //int verbose_level = pest_scenario.get_pestpp_options().get_ies_verbose_level();
    if (num_threads > 1)
    {
        performance_log->log_event("starting multithreaded MM neighbor calcs");
        ss.str("");
        Eigen::setNbThreads(1);
        vector<thread> threads;
        vector<exception_ptr> exception_ptrs;
        vector<string> preal_names = pe.get_real_names();
        vector<string> oreal_names = oe.get_real_names();
        unordered_map<string,Eigen::VectorXd> real_vec_map;
        unordered_map<string,unordered_map<string,double>> neighbor_phi_map, neighbor_pardist_map;
        mm_q_vec_map.clear();
        for (int i=0;i<preal_names.size();i++) {
            real_vec_map[preal_names[i]] = pe.get_eigen_ptr()->row(i);
            mm_q_vec_map[preal_names[i]] = wmat.row(i);
        }
        map<string,int> real_map = pe.get_real_map();

        MmNeighborThread *ut_ptr = new MmNeighborThread(real_vec_map,mm_real_idx_map,mm_real_name_map,neighbor_phi_map,neighbor_pardist_map);
        for (int i = 0; i < num_threads; i++) {
            exception_ptrs.push_back(exception_ptr());
        }
        for (int i = 0; i < num_threads; i++) {
            threads.push_back(thread(mm_neighbor_thread_function, i, verbose_level, mm_alpha, weight_phi_map,preal_names,
                                     oreal_names, parcov_inv, real_map, std::ref(*ut_ptr), std::ref(exception_ptrs[i])));
        }

        int num_exp = 0;

        for (int i = 0; i < num_threads; ++i) {
            bool found = false;
            if (exception_ptrs[i]) {
                found = true;
                num_exp++;
                try {
                    rethrow_exception(exception_ptrs[i]);
                }
                catch (const std::exception &e) {
                    //ss.str("");
                    ss << " thread " << i << "raised an exception: " << e.what();
                    //throw runtime_error(ss.str());
                }
                catch (...) {
                    //ss.str("");
                    ss << " thread " << i << "raised an exception";

                }
            }
            threads[i].join();
            if ((exception_ptrs[i]) && (!found)) {
                num_exp++;
                try {
                    rethrow_exception(exception_ptrs[i]);
                }
                catch (const std::exception &e) {
                    //ss.str("");
                    ss << " thread " << i << "raised an exception: " << e.what();

                }
                catch (...) {
                    //ss.str("");
                    ss << " thread " << i << "raised an exception: ";

                }
            }
        }
        if (num_exp > 0) {
            throw runtime_error(ss.str());
        }
        delete ut_ptr;

        if (pest_scenario.get_pestpp_options().get_ies_verbose_level() > 1) {
            int subset_size = (int) (((double) pe.shape().first) * mm_alpha);
            ofstream csv;
            ss.str("");
            ss << file_manager.get_base_filename() << "." << iter << ".mm.info.csv";
            csv.open(ss.str());
            csv << "pe_real_name,oe_real_name";
            for (int j = 0; j < subset_size; j++) {
                csv << ",par_real_neighbor_" << j << ",phi,pdiff";
            }
            csv << endl;
            string prname, orname;
            //for (auto& names : mm_real_name_map)
            for (auto real_name : preal_names)
            {
                auto names = mm_real_name_map.at(real_name);
                csv << real_name << "," << mm_real_name_map.at(real_name).second[0];
                for (int i=0;i<names.first.size();i++)
                {
                    prname = names.first[i];
                    orname = names.second[i];
                    if (prname == real_name)
                        continue;
                    csv  << "," << prname << "," << neighbor_phi_map.at(real_name).at(prname) << ",";
                    csv << neighbor_pardist_map.at(real_name).at(prname);
                }
                csv << endl;
            }
            csv.close();
        }
        return;
    }



    //util functions
    typedef std::function<bool(std::pair<std::string, double>, std::pair<std::string, double>)> Comparator;
    // Defining a lambda function to compare two pairs. It will compare two pairs using second field
    Comparator compFunctor = [](std::pair<std::string, double> elem1, std::pair<std::string, double> elem2) {
        return elem1.second < elem2.second;
    };
    typedef std::set<std::pair<std::string, double>, Comparator> sortedset;
    typedef std::set<std::pair<std::string, double>, Comparator>::iterator sortedset_iter;

    int subset_size = (int) (((double) pe.shape().first) * mm_alpha);
    if (subset_size <= 0)
        throw runtime_error("solve_multimodal error: subsize size zero or negative");
    ss.str("");
    ss << "calculating multimodal neighborhoods with " << subset_size << " realizations";
    performance_log->log_event(ss.str());

    vector<string> real_names = pe.get_real_names(), oreal_names = oe.get_real_names(),preal_names;
    string real_name, oreal_name;

    map<string, double> euclid_par_dist;
    Eigen::VectorXd real, diff;

    double edist;
    map<string, int> real_map = pe.get_real_map();
    oe.update_var_map();
    map<string, int> ovar_map = oe.get_var_map();

    vector<int> real_idxs;
    vector<string> pe_real_names_case,oe_real_names_case;
    //Eigen::MatrixXd* real_ptr = pe_upgrade.get_eigen_ptr_4_mod();

    ofstream csv;
    if (pest_scenario.get_pestpp_options().get_ies_verbose_level() > 1) {
        ss.str("");
        ss << file_manager.get_base_filename() << "." << iter << "." << ".mm.info.csv";
        csv.open(ss.str());
        csv << "pe_real_name,oe_real_name";
        for (int j = 0; j < subset_size; j++) {
            csv << ",par_real_neighbor_" << j << ",phi,pdiff";
        }
        csv << endl;

    }

    string prname, orname;


    Eigen::VectorXd q_vec;
    map<string, double> phi_map;
    for (int i = 0; i < pe.shape().first; i++) {
        real_name = real_names[i];

        performance_log->log_event("calculating multimodal neighborhood for realization " + real_name);

        euclid_par_dist.clear();
        real = pe.get_real_vector(real_name);

        //get the phi map for this realization using the weights vector for this realization
        oreal_name = oreal_names[i];
        performance_log->log_event("...getting weights");
        q_vec = wmat.row(i);
        oreal_name = oreal_names[i];
        phi_map = weight_phi_map.at(oreal_name);
        double mx = -1.0e+300;
        for (auto &p : phi_map) {
            if (p.second > mx)
                mx = p.second;
        }
        for (auto &p : phi_map) {
            p.second /= mx;
        }
        //flip to map to par realization names
        map<string, double> par_phi_map;
        oreal_names = oe.get_real_names();
        preal_names = pe.get_real_names();
        for (int ii = 0; ii < oreal_names.size(); ii++) {
            par_phi_map[preal_names[ii]] = phi_map.at(oreal_names[ii]);
        }
        phi_map.clear();

        //cout << real << endl;
        performance_log->log_event("...getting par distances");
        for (auto &rname : real_names) {
            if (rname == real_name)
                continue;
            diff = real - pe.get_real_vector(rname);
            edist = diff.transpose() * parcov_inv * diff;
            euclid_par_dist[rname] = edist;

        }

        mx = -1e300;
        for (auto &e : euclid_par_dist) {
            if (e.second > mx)
                mx = e.second;
        }
        if (mx <= 0.0) {
            ss.str("");
            ss << "multimodal solve error: maximum par diff for realization '" << real_name << "' not valid";
            message(0, ss.str());
        } else {
            for (auto &e : euclid_par_dist) {
                e.second /= mx;
            }
        }

        map<string, double> composite_score;
        for (auto &rname : real_names) {
            if (rname == real_name)
                continue;
            composite_score[rname] = euclid_par_dist.at(rname) + par_phi_map.at(rname);
        }
        performance_log->log_event("...sorting composite score");
        sortedset fitness_sorted(composite_score.begin(), composite_score.end(), compFunctor);
        set<string> dups;
        if (fitness_sorted.size() != composite_score.size())
        {
            performance_log->log_event("duplicates in composite score...dealing");
            map<double,int> counts;
            for (auto& score : composite_score)
            {
                if (counts.find(score.second) == counts.end())
                    counts[score.second] = 0;
                counts[score.second]++;
            }
            for (auto& score : composite_score)
            {
                if (counts[score.second] > 1)
                {
                    bool isin = false;
                    for (auto &ii: fitness_sorted) {
                        if (ii.first == score.first) {
                            isin = true;
                            break;
                        }
                    }

                    if (!isin) {
                        dups.insert(score.first);
                    }
                }

            }

        }
        real_idxs.clear();
        oe_real_names_case.clear();
        pe_real_names_case.clear();
        real_idxs.push_back(real_map.at(real_name));
        pe_real_names_case.push_back(real_name);
        oe_real_names_case.push_back(oreal_name);
        int iii = 0;
        for (sortedset_iter ii = fitness_sorted.begin(); ii != fitness_sorted.end(); ++ii) {
            int idx = real_map.at(ii->first);
            real_idxs.push_back(idx);
            pe_real_names_case.push_back(ii->first);
            oe_real_names_case.push_back(oreal_names[idx]);
            iii++;
            if (iii >= subset_size)//plus one to count 'real_name'
                break;
            //cout << real_name << "," << ii->first << "," << ii->second << endl;

        }
        for (auto& dup: dups)
        {
            real_idxs.push_back(real_map.at(dup));
        }
        if (real_idxs.size() != subset_size+1)
        {

            throw runtime_error("multimodal_solve: real_idxs.size() != subset_size");

        }
        if (pest_scenario.get_pestpp_options().get_ies_verbose_level() > 1) {
            csv << real_name;
        //for (auto &rname : pe_real_names_case) {
        for (int ii=0;ii<pe_real_names_case.size();ii++)
        {
            prname = pe_real_names_case[ii];
            orname = oe_real_names_case[ii];
            if (prname == real_name)
                continue;
            csv << "," << prname << ","  << par_phi_map.at(prname) << "," << euclid_par_dist.at(prname);
            }

            csv << endl;
        }
        mm_real_idx_map[real_name] = real_idxs;
        mm_q_vec_map[real_name] = q_vec;
        mm_real_name_map[real_name] = make_pair(pe_real_names_case,oe_real_names_case);
    }
    csv.close();
}

void EnsembleSolver::initialize_for_mm_solve()
{
    vector<string> pe_real_names = pe.get_real_names();
    vector<string> oe_real_names = oe.get_real_names();
    vector<string> t;
    vector<string> obs_names = oe.get_var_names();


    //prep a fast look for obs, par and resid matrices - stored as column vectors in a map
    par_resid_map.clear();
    par_resid_map.reserve(pe_real_names.size());
    par_diff_map.clear();
    par_diff_map.reserve(pe_real_names.size());

    obs_resid_map.clear();
    obs_resid_map.reserve(oe_real_names.size());
    obs_diff_map.clear();
    obs_diff_map.reserve(oe_real_names.size());
    obs_err_map.clear();
    obs_err_map.reserve(oe_real_names.size());

    //check for the 'center_on' real - it may have been dropped...

    Eigen::MatrixXd mat = ph.get_obs_resid_subset(oe,true,oe_real_names);
    for (int i = 0; i < oe_real_names.size(); i++)
    {
        obs_resid_map[oe_real_names[i]] = mat.row(i);
    }

    //oe should be in the right row and col order...
    for (int i = 0; i < oe_real_names.size(); i++)
    {
        obs_diff_map[oe_real_names[i]] = oe.get_eigen_ptr()->row(i);
    }

    mat = base_oe.get_eigen(oe_real_names, obs_names);
    mat.transposeInPlace();
    Observations ctl_obs = pest_scenario.get_ctl_observations();
    Eigen::ArrayXd ovals = ctl_obs.get_data_eigen_vec(obs_names).array();
    for (int i = 0; i < oe_real_names.size(); i++)
    {
        obs_err_map[oe_real_names[i]] = mat.col(i).array() - ovals;
    }

    mat = ph.get_par_resid_subset(pe,pe_real_names);
    for (int i = 0; i < pe_real_names.size(); i++)
    {
        par_resid_map[pe_real_names[i]] = mat.row(i);
    }

    //vector<string> par_names = pe.get_var_names();
    for (int i = 0; i < pe_real_names.size(); i++)
    {
        par_diff_map[pe_real_names[i]] = pe.get_eigen_ptr()->row(i);
    }

    mat.resize(0, 0);
}


void EnsembleSolver::initialize_for_localized_solve(string center_on, vector<int> real_idxs)
{
    vector<string> pe_real_names = pe.get_real_names();
    vector<string> oe_real_names = oe.get_real_names();
    vector<string> t;
    if (center_on.size() == 0)
    {
        center_on = pest_scenario.get_pestpp_options().get_ies_center_on();
    }
    if (real_idxs.size() == 0)
    {
        for (int i=0;i<pe.shape().first;i++)
            real_idxs.push_back(i);
    }
    else
    {
        for (auto ri : real_idxs)
        {
            t.push_back(pe_real_names[ri]);
        }
        pe_real_names = t;
        t.clear();
        for (auto ri: real_idxs)
        {
            t.push_back(oe_real_names[ri]);
        }
        oe_real_names = t;
        t.clear();
    }

	parcov_inv_map.clear();
	parcov_inv_map.reserve(pe.shape().second);
	Eigen::VectorXd parcov_inv;// = parcov.get(par_names).inv().e_ptr()->toDense().cwiseSqrt().asDiagonal();
	if (!parcov.isdiagonal())
	{
		parcov_inv = parcov.get_matrix().diagonal();
	}
	else
	{
		//message(2,"extracting diagonal from prior parameter covariance matrix");
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

	Eigen::MatrixXd mat = ph.get_obs_resid_subset(oe,true,oe_real_names);
	for (int i = 0; i < obs_names.size(); i++)
	{
		obs_resid_map[obs_names[i]] = mat.col(i);
	}
	mat = oe.get_eigen_anomalies(oe_real_names, t, center_on);
	for (int i = 0; i < obs_names.size(); i++)
	{
		obs_diff_map[obs_names[i]] = mat.col(i);
	}
	mat = base_oe.get_eigen(oe_real_names, obs_names);
	Observations ctl_obs = pest_scenario.get_ctl_observations();
	for (int i = 0; i < obs_names.size(); i++)
	{
		obs_err_map[obs_names[i]] = mat.col(i).array() - ctl_obs.get_rec(obs_names[i]);
	}
	mat = ph.get_par_resid_subset(pe,pe_real_names);
	for (int i = 0; i < par_names.size(); i++)
	{
		par_resid_map[par_names[i]] = mat.col(i);
	}
	mat = pe.get_eigen_anomalies(pe_real_names,t,center_on);
	for (int i = 0; i < par_names.size(); i++)
	{
		par_diff_map[par_names[i]] = mat.col(i);
	}
	if ((!pest_scenario.get_pestpp_options().get_ies_use_approx() && (Am.rows() > 0)))
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

//void work(int thread_id, int iter, double cur_lam, bool use_glm_form, Eigen::VectorXd parcov_inv_vec, Eigen::MatrixXd Am);
void upgrade_thread_function_mm(int id, int iter, double cur_lam, bool use_glm_form, Eigen::VectorXd parcov_inv_vec,
                             Eigen::MatrixXd Am, MmUpgradeThread& worker, exception_ptr& eptr) {
    //worker.work(id, iter, cur_lam, use_glm_form, parcov_inv_vec, Am);
    try {
        worker.work(id, iter, cur_lam, use_glm_form, parcov_inv_vec, Am);
    }
    catch (...) {
        eptr = current_exception();
    }

    return;
}


void EnsembleSolver::solve_multimodal(int num_threads, double cur_lam, bool use_glm_form, ParameterEnsemble& pe_upgrade,
                                      unordered_map<string, pair<vector<string>, vector<string>>>& loc_map,
                                      double mm_alpha) {

    stringstream ss;
    if ((!localizer.get_use()) && (num_threads > 1))
    {
        initialize_for_mm_solve();
        Eigen::setNbThreads(1);
        vector<thread> threads;
        vector<exception_ptr> exception_ptrs;
        message(2, "launching threads");

        MmUpgradeThread* ut_ptr = new MmUpgradeThread(performance_log, par_resid_map, par_diff_map, obs_resid_map, obs_diff_map, obs_err_map,
                                                                        mm_q_vec_map, pe_upgrade,mm_real_name_map,reg_factor);

        Eigen::VectorXd parcov_inv_vec = 1. / parcov.e_ptr()->diagonal().array();
        for (int i = 0; i < num_threads; i++)
        {
            exception_ptrs.push_back(exception_ptr());
        }

        for (int i = 0; i < num_threads; i++)
        {
            threads.push_back(thread(upgrade_thread_function_mm, i, iter, cur_lam, use_glm_form, parcov_inv_vec,
                                     Am, std::ref(*ut_ptr), std::ref(exception_ptrs[i])));



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
        message(1, "upgrade calculation done");

        return;
    }
    vector<string> real_names = pe.get_real_names(),oreal_names = oe.get_real_names(),upgrade_real_names;
    string real_name;
    vector<int> real_idxs;
    Eigen::VectorXd q_vec,real;
    for (int i=0;i<pe.shape().first;i++) {
        real_name = real_names[i];
        real_idxs.clear();
        upgrade_real_names.clear();
        real_idxs = mm_real_idx_map.at(real_name);
        q_vec = mm_q_vec_map.at(real_name);
        for (auto idx : real_idxs)
            upgrade_real_names.push_back(real_names[idx]);

        ParameterEnsemble pe_real(&pest_scenario, pe.get_rand_gen_ptr());
        pe_real.reserve(upgrade_real_names, pe.get_var_names());

        if (!localizer.get_use()) {
            performance_log->log_event("calculating nonlocalized multimodal upgrade for " + real_name);
            //get weights here
            //Eigen::DiagonalMatrix<double, Eigen::Dynamic> weights(q_vec);
            nonlocalized_solve(cur_lam, use_glm_form, pe_real, real_name, real_idxs, q_vec);
            //solve(num_threads,cur_lam,use_glm_form,pe_real,loc_map);
        } else {
            performance_log->log_event("re-initializing solver for localized multimodal upgrade for " + real_name);
            initialize_for_localized_solve(string(), real_idxs);
            //reset the weight map
            weight_map.clear();
            weight_map.reserve(q_vec.size());
            int iii = 0;
            for (auto &name : act_obs_names) {
                weight_map[name] = q_vec[iii];
                iii++;
            }
            performance_log->log_event("calculating localized multimodal upgrade for " + real_name);

            solve(num_threads, cur_lam, use_glm_form, pe_real, loc_map);
        }


        real = pe_real.get_real_vector(real_name);
        //cout << real << endl;
        pe_upgrade.update_real_ip(real_name, real);
    }
}

void EnsembleSolver::nonlocalized_solve(double cur_lam,bool use_glm_form, ParameterEnsemble& pe_upgrade,
                        string center_on, vector<int> real_idxs,Eigen::VectorXd q_vec)
{
    pe_upgrade.set_zeros();
    if (real_idxs.size() == 0)
    {
        for (int i = 0;i<pe_upgrade.shape().first;i++)
            real_idxs.push_back(i);
    }
    vector<string> names = oe.get_real_names();
    vector<string> oe_real_names;
    for (auto i : real_idxs)
        oe_real_names.push_back(names[i]);
    vector<string> pe_real_names = pe_upgrade.get_real_names();
    string obs_center_on;
    if ((center_on.size() > 0) && (center_on != MEDIAN_CENTER_ON_NAME)) {
        for (int i = 0; i < pe_real_names.size(); i++) {
            if (pe_real_names[i] == center_on) {
                obs_center_on = oe_real_names[i];
                break;
            }
        }
        if (obs_center_on.size() == 0) {
            throw runtime_error(
                    "EnsembleSolver::nonlocalized_solve: obs_center_on not found for center_on: " + center_on);
        }
    }
    Eigen::MatrixXd obs_resid = ph.get_obs_resid_subset(oe,true,oe_real_names);

    Eigen::MatrixXd obs_diff = oe.get_eigen_anomalies(oe_real_names, act_obs_names, obs_center_on);

    Eigen::MatrixXd obs_err = base_oe.get_eigen(oe_real_names, act_obs_names);
    Observations ctl_obs = pest_scenario.get_ctl_observations();
    for (int j = 0; j < act_obs_names.size(); j++)
    {
        obs_err.col(j) = obs_err.col(j).array() - ctl_obs.get_rec(act_obs_names[j]);
    }
    Eigen::MatrixXd par_resid = ph.get_par_resid_subset(pe,pe_real_names);

    Eigen::MatrixXd par_diff = pe.get_eigen_anomalies(pe_real_names,act_par_names,center_on);
    int maxsing = pest_scenario.get_svd_info().maxsing;
    bool use_prior_scaling = pest_scenario.get_pestpp_options().get_ies_use_prior_scaling();
    bool use_approx = pest_scenario.get_pestpp_options().get_ies_use_approx();
    double eigthresh = pest_scenario.get_svd_info().eigthresh;
    if (q_vec.size() == 0)
    {
        q_vec.resize(act_obs_names.size());
        if ((center_on.size() > 0) && (center_on != MEDIAN_CENTER_ON_NAME))
        {

            map<string, int> oreal_map = oe.get_real_map();
            int ireal = oreal_map.at(obs_center_on);
            map<string, int> ovar_map = oe.get_var_map();
            int ii = 0;
            for (auto &aon : act_obs_names)
            {
                q_vec[ii] = weights.get_eigen_ptr()->coeff(ireal, ovar_map.at(aon));
                ii++;
            }
        }
        else
        {
            int ii = 0;
            ObservationInfo *oi_ptr = pest_scenario.get_observation_info_ptr();
            for (auto &aon : act_obs_names) {
                q_vec[ii] = oi_ptr->get_weight(aon);
                ii++;
            }
        }
    }
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> local_weights(q_vec);

    Eigen::VectorXd parcov_inv_vec;// = parcov.get(par_names).inv().e_ptr()->toDense().cwiseSqrt().asDiagonal();
    if (parcov.get_row_names().size() != act_par_names.size())
    {
        throw runtime_error("parcov not aligned with act par names");
    }
    if (parcov.isdiagonal())
    {
        parcov_inv_vec = parcov.get_matrix().diagonal();
    }
    else
    {
        Covariance parcov_diag;
        parcov_diag.from_diagonal(parcov);
        parcov_inv_vec = parcov_diag.get_matrix().diagonal();
    }
    parcov_inv_vec = parcov_inv_vec.cwiseSqrt().cwiseInverse();
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> parcov_inv(parcov_inv_vec);
    Eigen::MatrixXd upgrade_1;
    par_diff.transposeInPlace();
    obs_diff.transposeInPlace();
    obs_resid.transposeInPlace();
    par_resid.transposeInPlace();
    obs_err.transposeInPlace();
    UpgradeThread::ensemble_solution(iter,verbose_level,maxsing,0,0,use_prior_scaling,use_approx,use_glm_form,cur_lam,eigthresh,par_resid,
                      par_diff,Am,obs_resid,obs_diff,upgrade_1,obs_err,local_weights,parcov_inv, act_obs_names,act_par_names,reg_factor);
    pe_upgrade.add_2_cols_ip(act_par_names, upgrade_1);


}

void EnsembleSolver::solve(int num_threads, double cur_lam, bool use_glm_form, ParameterEnsemble& pe_upgrade, unordered_map<string, pair<vector<string>, vector<string>>>& loc_map)
{


    if (!localizer.get_use())
    {
        string center_on = pest_scenario.get_pestpp_options().get_ies_center_on();
        nonlocalized_solve(cur_lam,use_glm_form,pe_upgrade,center_on);
        return;
    }
	pe_upgrade.set_zeros();
	stringstream ss;
	Localizer::How _how = localizer.get_how();
	Localizer::LocTyp loctyp = localizer.get_loctyp();
	bool use_cov_loc = true;
	if (loctyp == Localizer::LocTyp::COVARIANCE) {
        throw runtime_error("EnsembleSolver::solve(): 'covariance' localization is deprecated");
    }
	//LocalAnalysisUpgradeThread worker(performance_log, par_resid_map, par_diff_map, obs_resid_map, obs_diff_map,obs_err_map,
	//	localizer, parcov_inv_map, weight_map, pe_upgrade, loc_map, Am_map, _how);
	UpgradeThread* ut_ptr;
    ut_ptr = new LocalAnalysisUpgradeThread(performance_log, par_resid_map, par_diff_map, obs_resid_map,
                                            obs_diff_map, obs_err_map,
                                            localizer, parcov_inv_map, weight_map, pe_upgrade, loc_map, Am_map,
                                            _how, reg_factor);
    performance_log->log_event("using local analysis upgrade thread");

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
		    performance_log->log_event(ss.str());
			throw runtime_error(ss.str());
		}
		delete ut_ptr;
		message(1, "upgrade calculation done");
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


UpgradeThread::UpgradeThread(PerformanceLog* _performance_log, unordered_map<string, Eigen::VectorXd>& _par_resid_map, unordered_map<string, Eigen::VectorXd>& _par_diff_map,
                             unordered_map<string, Eigen::VectorXd>& _obs_resid_map, unordered_map<string, Eigen::VectorXd>& _obs_diff_map, unordered_map<string, Eigen::VectorXd>& _obs_err_map,
                             Localizer& _localizer, unordered_map<string, double>& _parcov_inv_map,
                             unordered_map<string, double>& _weight_map, ParameterEnsemble& _pe_upgrade,
                             unordered_map<string, pair<vector<string>, vector<string>>>& _cases,
                             unordered_map<string, Eigen::VectorXd>& _Am_map, Localizer::How& _how, double _reg_factor):
	par_resid_map(_par_resid_map),par_diff_map(_par_diff_map), obs_resid_map(_obs_resid_map), 
	obs_diff_map(_obs_diff_map), obs_err_map(_obs_err_map), localizer(_localizer),
    pe_upgrade(_pe_upgrade), cases(_cases), parcov_inv_map(_parcov_inv_map), 
	weight_map(_weight_map), Am_map(_Am_map), reg_factor(_reg_factor)
	{
		performance_log = _performance_log;
		how = _how;
		count = 0;

		for (auto& c : cases)
		{
			keys.push_back(c.first);
		}
		total = keys.size();
	
}

void UpgradeThread::ensemble_solution(const int iter, const int verbose_level,const int maxsing,  const int thread_id,
                              const int t_count, const bool use_prior_scaling,const bool use_approx, const bool use_glm,
                              const double cur_lam,const double eigthresh, Eigen::MatrixXd& par_resid, Eigen::MatrixXd& par_diff,
                              const Eigen::MatrixXd& Am, Eigen::MatrixXd& obs_resid,Eigen::MatrixXd& obs_diff, Eigen::MatrixXd& upgrade_1,
                              Eigen::MatrixXd& obs_err, const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& weights,
                              const Eigen::DiagonalMatrix<double, Eigen::Dynamic>& parcov_inv,
                              const vector<string>& act_obs_names,const vector<string>& act_par_names,double _reg_factor)
{
    class local_utils
    {
    public:

        static void save_names(int verbose_level, int tid, int iter, int t_count, string prefix, const vector<string>& names)
        {
            if (verbose_level < 2)
                return;

            if (verbose_level < 3)
                return;
            stringstream ss;
            ss << "thread_" << tid << ".count_ " << t_count << ".iter_" << iter << "." << prefix << ".dat";
            string fname = ss.str();
            ofstream f(fname);
            if (!f.good())
                cout << "error getting ofstream " << fname << endl;
            else
            {
                for (const auto& name : names)
                {
                    f << name << endl;
                }

            }
            f.close();
            return;
        }

        static void save_mat(int verbose_level, int tid, int iter, int t_count, string prefix, const Eigen::MatrixXd& mat)
        {
            if (verbose_level < 2)
                return;

            if (verbose_level < 3)
                return;
            //cout << "thread: " << tid << ", " << t_count << ", " << prefix << " rows:cols" << mat.rows() << ":" << mat.cols() << endl;
            stringstream ss;

            ss << "thread_" << tid << ".count_" << t_count << ".iter_" << iter << "." << prefix << ".dat";
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

    local_utils::save_names(verbose_level,thread_id, iter, t_count,"act_obs_names", act_obs_names);
    local_utils::save_names(verbose_level,thread_id, iter, t_count,"act_par_names", act_par_names);

    stringstream ss;
    if (!use_glm)
    {
        if (true)
        {
            // Low rank Cee. Section 14.3.2 Evenson Book
            obs_err = obs_err.colwise() - obs_err.rowwise().mean();
            obs_err = sqrt(cur_lam) * obs_err;
            Eigen::MatrixXd s0, V0, U0, s0_i;
            SVD_REDSVD rsvd;
            rsvd.solve_ip(obs_diff, s0, U0, V0, eigthresh, maxsing);
            s0_i = s0.asDiagonal().inverse();
            Eigen::MatrixXd X0 = U0.transpose() * obs_err;
            X0 = s0_i * X0;
            Eigen::MatrixXd s1, V1, U1, s1_2, s1_2i;
            rsvd.solve_ip(X0, s1, U1, V1, 0, maxsing);

            s1_2 = s1.cwiseProduct(s1);
            s1_2i = (Eigen::VectorXd::Ones(s1_2.size()) + s1_2).asDiagonal().inverse();
            Eigen::MatrixXd X1 = s0_i * U1;
            X1 = U0 * X1;

            Eigen::MatrixXd X4 = s1_2i * X1.transpose();
            Eigen::MatrixXd X2 = X4 * obs_resid;
            Eigen::MatrixXd X3 = X1 * X2;

            X3 = obs_diff.transpose() * X3;
            upgrade_1 = -1 * par_diff * X3;
            local_utils::save_mat(verbose_level, thread_id, iter, t_count, "upgrade_1", upgrade_1);
            upgrade_1.transposeInPlace();

        }
        else
        {
            // Use when ensemble size is larger than number of observation.
            // This is a good option when number of observation is small
            local_utils::save_mat(verbose_level, thread_id, iter, t_count, "obs_err", obs_err);
            obs_err = obs_err.colwise() - obs_err.rowwise().mean();

            Eigen::MatrixXd s2_, s, V, U, cum_sum;
            SVD_REDSVD rsvd;
            Eigen::MatrixXd C;
            C = obs_diff + (sqrt(cur_lam) * obs_err); // curr_lam is the inflation factor
            local_utils::save_mat(verbose_level, thread_id, iter, t_count, "C", C);
            Eigen::VectorXd s2;


            rsvd.solve_ip(C, s, U, V, eigthresh, maxsing);
            local_utils::save_mat(verbose_level, thread_id, iter, t_count, "s", s);
            local_utils::save_mat(verbose_level, thread_id, iter, t_count, "U", U);
            local_utils::save_mat(verbose_level, thread_id, iter, t_count, "V", V);
            s2 = s.cwiseProduct(s);
            s2_ = s2.asDiagonal().inverse();
            local_utils::save_mat(verbose_level, thread_id, iter, t_count, "inv_s2_", s2_);


            Eigen::MatrixXd X1 = s2_ * U.transpose();
            local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X1", X1);

            X1 = X1 * obs_resid;
            local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X1_obs_resid", X1);

            X1 = U * X1;
            local_utils::save_mat(verbose_level, thread_id, iter, t_count, "U_X1", X1);

            X1 = obs_diff.transpose() * X1;
            local_utils::save_mat(verbose_level, thread_id, iter, t_count, "obs_diff_X1", X1);

            upgrade_1 = -1 * par_diff * X1;
            local_utils::save_mat(verbose_level, thread_id, iter, t_count, "upgrade_1", upgrade_1);
            upgrade_1.transposeInPlace();
        }

    }
    else
    {

        Eigen::MatrixXd ivec, s, s2, V, Ut, d_dash;
        string key;
        Eigen::MatrixXd wtemp = weights.diagonal().matrix();
        local_utils::save_mat(verbose_level, thread_id, iter, t_count, "weights", wtemp);
        local_utils::save_mat(verbose_level, thread_id, iter, t_count, "obs_resid", obs_resid);
        obs_resid = weights * obs_resid;
        local_utils::save_mat(verbose_level, thread_id, iter, t_count, "scaled_obs_resid", obs_resid);

        int num_reals = par_resid.cols();
        double scale = (1.0 / (sqrt(double(num_reals - 1))));
        local_utils::save_mat(verbose_level, thread_id, iter, t_count, "obs_diff", obs_diff);
        obs_diff = scale * (weights * obs_diff);
        local_utils::save_mat(verbose_level, thread_id, iter, t_count, "scaled_obs_diff", obs_diff);
        local_utils::save_mat(verbose_level, thread_id, iter, t_count, "par_diff", par_diff);
        if (verbose_level > 1) {
            if (parcov_inv.size() < 10000) {

                Eigen::MatrixXd temp = parcov_inv.toDenseMatrix();
                local_utils::save_mat(verbose_level, thread_id, iter, t_count, "parcov_inv", temp);
            }
            else
            {
                Eigen::MatrixXd temp = parcov_inv.diagonal().matrix();
                local_utils::save_mat(verbose_level, thread_id, iter, t_count, "parcov_inv", temp);
            }

        }
        if (use_prior_scaling)

            par_diff = scale * parcov_inv * par_diff;
        else
            par_diff = scale * par_diff;
        local_utils::save_mat(verbose_level, thread_id, iter, t_count, "scaled_par_diff", par_diff);
        SVD_REDSVD rsvd;

        rsvd.solve_ip(obs_diff, s, Ut, V, eigthresh, maxsing);

        Ut.transposeInPlace();
        //obs_diff.resize(0, 0);
        local_utils::save_mat(verbose_level, thread_id, iter, t_count, "Ut", Ut);
        local_utils::save_mat(verbose_level, thread_id, iter, t_count, "s", s);
        local_utils::save_mat(verbose_level, thread_id, iter, t_count, "V", V);

        s2 = s.cwiseProduct(s);

        ivec = ((Eigen::VectorXd::Ones(s2.size()) * (cur_lam + 1.0)) + s2).asDiagonal().inverse();
        local_utils::save_mat(verbose_level, thread_id, iter, t_count, "ivec", ivec);

        Eigen::MatrixXd X1 = Ut * obs_resid;
        local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X1", X1);

        Eigen::MatrixXd X2 = ivec * X1;
        //X1.resize(0, 0);
        local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X2", X2);

        Eigen::MatrixXd X3 = V * s.asDiagonal() * X2;
        //X2.resize(0, 0);
        local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X3", X3);
        upgrade_1 = -1.0 * par_diff * X3;

        if (use_prior_scaling) {
            //upgrade_1 = parcov_inv * upgrade_1;
        }

        upgrade_1.transposeInPlace();
        local_utils::save_mat(verbose_level, thread_id, iter, t_count, "upgrade_1", upgrade_1);
        //X3.resize(0, 0);

        Eigen::MatrixXd upgrade_2;
        if ((!use_approx) && (iter > 1)) {
            if (use_prior_scaling) {
                par_resid = parcov_inv * par_resid;
            }

            local_utils::save_mat(verbose_level, thread_id, iter, t_count, "Am", Am);
            Eigen::MatrixXd x4 = Am.transpose() * par_resid;
            local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X4", x4);

            par_resid.resize(0, 0);

            Eigen::MatrixXd x5 = Am * x4;
            //x4.resize(0, 0);
            //Am.resize(0, 0);

            local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X5", x5);
            Eigen::MatrixXd x6 = par_diff.transpose() * x5;
            //x5.resize(0, 0);

            local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X6", x6);
            Eigen::MatrixXd x7 = V * ivec * V.transpose() * x6;
            //x6.resize(0, 0);

            if (use_prior_scaling) {
                upgrade_2 = -1.0 * parcov_inv * par_diff * x7;
            } else {
                upgrade_2 = -1.0 * (par_diff * x7);
            }
            //x7.resize(0, 0);
            if (_reg_factor >= 0)
            {
                upgrade_1 = upgrade_1 + (_reg_factor * upgrade_2.transpose());
            }
            else {
                upgrade_1 = upgrade_1 + upgrade_2.transpose();
            }

            local_utils::save_mat(verbose_level, thread_id, iter, t_count, "upgrade_2", upgrade_2);
            //upgrade_2.resize(0, 0);

        }
    }

}


void MmUpgradeThread::work(int thread_id, int iter, double cur_lam, bool use_glm_form, Eigen::VectorXd parcov_inv_vec, Eigen::MatrixXd Am)
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
        static Eigen::MatrixXd get_matrix_from_map(vector<string>& row_names, unordered_map<string, Eigen::VectorXd>& emap)
        {
            int num_vars = emap.at(row_names[0]).size();
            Eigen::MatrixXd mat(row_names.size(),num_vars);
            mat.setZero();

            for (int i = 0; i < row_names.size(); i++)
            {
                mat.row(i) = emap.at(row_names[i]);
            }

            return mat;
        }
        static void save_mat(int verbose_level, int tid, int iter, int t_count, string prefix, Eigen::MatrixXd& mat)
        {
            if (verbose_level < 2)
                return;

            if (verbose_level < 3)
                return;
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
    int maxsing, verbose_level, pcount = 0, t_count;
    double eigthresh;
    bool use_approx;
    bool use_prior_scaling;

    while (true)
    {
        if (ctrl_guard.try_lock()) {
            maxsing = pe_upgrade.get_pest_scenario_ptr()->get_svd_info().maxsing;
            eigthresh = pe_upgrade.get_pest_scenario_ptr()->get_svd_info().eigthresh;
            use_approx = pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_ies_use_approx();
            use_prior_scaling = pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_ies_use_prior_scaling();
            verbose_level = pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_ies_verbose_level();
            ctrl_guard.unlock();
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
    Eigen::MatrixXd par_resid, par_diff;
    Eigen::MatrixXd obs_resid, obs_diff, obs_err, loc;
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> weights;
    string key;
    Eigen::VectorXd row_vec;


    //these locks are used to control (thread-safe) access to the fast look up containers
    unique_lock<mutex> next_guard(next_lock, defer_lock);
    unique_lock<mutex> obs_diff_guard(obs_diff_lock, defer_lock);
    unique_lock<mutex> obs_resid_guard(obs_resid_lock, defer_lock);
    unique_lock<mutex> obs_err_guard(obs_err_lock, defer_lock);
    unique_lock<mutex> par_diff_guard(par_diff_lock, defer_lock);
    unique_lock<mutex> par_resid_guard(par_resid_lock, defer_lock);
    unique_lock<mutex> weight_guard(weight_lock, defer_lock);
    unique_lock<mutex> parcov_guard(parcov_lock, defer_lock);
    unique_lock<mutex> put_guard(put_lock, defer_lock);

    vector<string> pe_real_names,oe_real_names;
    int num_reals;

    //form parcov_inv
    Eigen::DiagonalMatrix<double,Eigen::Dynamic> parcov_inv = parcov_inv_vec.asDiagonal();

    //This is the main thread loop - it continues until all upgrade pieces have been completed
    while (true)
    {
        pe_real_names.clear();
        oe_real_names.clear();
        while (true)
        {
            if (next_guard.try_lock())
            {
                //if all the pieces have been completed, return
                if (count == keys.size())
                {
                    if (verbose_level > 3)
                    {
                        cout << "upgrade thread: " << thread_id << " processed " << pcount << " upgrade parts" << endl;
                    }
                    if (f_thread.good())
                        f_thread.close();
                    return;
                }
                key = keys[count];
                pair<vector<string>,vector<string>> p = cases.at(key);
                pe_real_names = p.first;
                oe_real_names = p.second;
                if ((count % 100 == 0) && (count > 0))
                {
                    ss.str("");
                    ss << "upgrade thread progress: " << count << " of " << total << " parts done";
                    if (verbose_level > 3)
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
            for (auto& name : pe_real_names)
                f_thread << "," << name;
            for (auto& name : oe_real_names)
                f_thread << "," << name;
            f_thread << endl;
        }

        //reset all the solution parts
        par_resid.resize(0, 0);
        par_diff.resize(0, 0);
        obs_resid.resize(0, 0);
        obs_diff.resize(0, 0);
        obs_err.resize(0, 0);
        weights.resize(0);

        num_reals = pe_real_names.size();

        obs_diff = local_utils::get_matrix_from_map(oe_real_names, obs_diff_map);
        obs_resid = local_utils::get_matrix_from_map(oe_real_names, obs_resid_map);
        obs_err = local_utils::get_matrix_from_map(oe_real_names, obs_err_map);
        par_diff = local_utils::get_matrix_from_map(pe_real_names, par_diff_map);
        par_resid = local_utils::get_matrix_from_map(pe_real_names, par_resid_map);
        weights = weight_map.at(key).asDiagonal();

        //now form the deviations - both pars and obs
        row_vec = par_diff.row(0);
        for (int i=0;i<par_diff.rows();i++)
            par_diff.row(i) -= row_vec;
        row_vec = obs_diff.row(0);
        for (int i=0;i<obs_diff.rows();i++)
            obs_diff.row(i) -= row_vec;

        par_diff.transposeInPlace();
        obs_diff.transposeInPlace();
        obs_resid.transposeInPlace();
        par_resid.transposeInPlace();
        obs_err.transposeInPlace();

        stringstream ss;

        double scale = (1.0 / (sqrt(double(num_reals - 1))));

        //local_utils::save_mat(verbose_level, thread_id, iter, t_count, "obs_diff", obs_diff);


        Eigen::MatrixXd upgrade_1;

        //Eigen::MatrixXd upgrade_1;
        vector<string> empty_obs_names,empty_par_names;
        UpgradeThread::ensemble_solution(iter,verbose_level,maxsing,thread_id,t_count, use_prior_scaling,use_approx,use_glm_form,cur_lam,eigthresh,par_resid,par_diff,Am,obs_resid,
                          obs_diff,upgrade_1,obs_err,weights,parcov_inv,empty_obs_names,empty_par_names,reg_factor);


        //assuming that the fist row is the realization we are after...
        row_vec = upgrade_1.row(0);

        while (true)
        {
            if (put_guard.try_lock())
            {
                //pe_upgrade.add_2_cols_ip(par_names, upgrade_1);
                pe_upgrade.add_2_row_ip(key,row_vec);
                put_guard.unlock();
                break;
            }
        }
    }

}


MmUpgradeThread::MmUpgradeThread(PerformanceLog* _performance_log, unordered_map<string, Eigen::VectorXd>& _par_resid_map,
                                 unordered_map<string, Eigen::VectorXd>& _par_diff_map,
                                 unordered_map<string, Eigen::VectorXd>& _obs_resid_map, unordered_map<string, Eigen::VectorXd>& _obs_diff_map,
                                 unordered_map<string, Eigen::VectorXd>& _obs_err_map,
                                 unordered_map<string, Eigen::VectorXd>& _weight_map, ParameterEnsemble& _pe_upgrade,
                                 unordered_map<string, pair<vector<string>, vector<string>>>& _cases, double _reg_factor):
        par_resid_map(_par_resid_map),par_diff_map(_par_diff_map), obs_resid_map(_obs_resid_map),
        obs_diff_map(_obs_diff_map), obs_err_map(_obs_err_map),
        pe_upgrade(_pe_upgrade), cases(_cases),
        weight_map(_weight_map),
        reg_factor(_reg_factor)
{
    performance_log = _performance_log;
    count = 0;

    for (auto& c : cases)
    {
        keys.push_back(c.first);
    }
    //sort(keys.begin(), keys.end());
    total = keys.size();
    //random_shuffle(keys.begin(), keys.end());

}




void LocalAnalysisUpgradeThread::work(int thread_id, int iter, double cur_lam, bool use_glm_form,
									vector<string> par_names, vector<string> obs_names)
{

	//declare these helpers in here so they are thread safe...
	class local_utils
	{
	public:
        static void save_names(int verbose_level, int tid, int iter, int t_count, string prefix, vector<string>& names)
        {
            if (verbose_level < 2)
                return;

            if (verbose_level < 3)
                return;
            stringstream ss;

            ss << "thread_" << tid << ".count_ " << t_count << ".iter_" << iter << "." << prefix << ".dat";
            string fname = ss.str();
            ofstream f(fname);
            if (!f.good())
                cout << "error getting ofstream " << fname << endl;
            else
            {
                for (const auto& name : names)
                {
                    f << name << endl;
                }

            }
            f.close();
            return;
        }

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
				mat.col(j) = emap.at(names[j]);
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
		ss << "thread_" << thread_id << ".part_map.csv";
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

	set<string> sorg_par_names(par_names.begin(),par_names.end());
	set<string> sorg_obs_names(obs_names.begin(),obs_names.end());


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
					if (verbose_level > 3)
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

					use_localizer = true;
				}
				if ((count % 1000 == 0) && (count > 0))
				{
					ss.str("");
					ss << "upgrade thread progress: " << count << " of " << total << " parts done";
					if (verbose_level > 3)
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

		

		if ((use_localizer) && (loc.rows() == 0)) {
            while (true) {
                //get access to the localizer
                if (loc_guard.try_lock()) {
                    //get a matrix that is either the shape of par diff or obs diff
                    if (loc_by_obs) {
                        //loc = localizer.get_localizing_par_hadamard_matrix(num_reals, obs_names[0], par_names);
                        loc = localizer.get_pardiff_hadamard_matrix(num_reals, key, par_names);
                    } else {
                        //loc = localizer.get_localizing_obs_hadamard_matrix(num_reals, par_names[0], obs_names);
                        loc = localizer.get_obsdiff_hadamard_matrix(num_reals, key, obs_names);
                    }
                    loc_guard.unlock();
                    break;
                }

            }
        }

        obs_diff = local_utils::get_matrix_from_map(num_reals, obs_names, obs_diff_map);
        obs_resid = local_utils::get_matrix_from_map(num_reals, obs_names, obs_resid_map);
        obs_err = local_utils::get_matrix_from_map(num_reals, obs_names, obs_err_map);
        par_diff = local_utils::get_matrix_from_map(num_reals, par_names, par_diff_map);
        par_resid = local_utils::get_matrix_from_map(num_reals, par_names, par_resid_map);
        weights = local_utils::get_matrix_from_map(obs_names, weight_map);
        parcov_inv = local_utils::get_matrix_from_map(par_names, parcov_inv_map);

        if ((!use_approx) && (Am.rows() == 0))
        {
            //Am = local_utils::get_matrix_from_map(num_reals, par_names, Am_map).transpose();
            int am_cols = Am_map[par_names[0]].size();
            Am.resize(par_names.size(), am_cols);
            Am.setZero();

            for (int j = 0; j < par_names.size(); j++)
            {
                Am.row(j) = Am_map[par_names[j]];
            }
        }

        Eigen::MatrixXd ivec, upgrade_1, s, s2, V, Ut, d_dash;
        par_diff.transposeInPlace();
		obs_diff.transposeInPlace();
		obs_resid.transposeInPlace();
		par_resid.transposeInPlace();
		obs_err.transposeInPlace();

        UpgradeThread::ensemble_solution(iter, verbose_level,maxsing,  thread_id,
                                      t_count, use_prior_scaling,use_approx, use_glm_form,
                                      cur_lam,eigthresh, par_resid, par_diff,
                                      Am, obs_resid,obs_diff, upgrade_1,
                                      obs_err, weights,
                                      parcov_inv,
                                      obs_names,par_names,reg_factor);

		
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
    map<string,double> extfile_gt = pest_scenario->get_ext_file_double_map("observation data external","greater_than");
    map<string,double> extfile_lt = pest_scenario->get_ext_file_double_map("observation data external","less_than");

    for (auto &oname : pest_scenario->get_ctl_ordered_obs_names())
	{
		og = oi->get_group(oname);
		weight = oi->get_weight(oname);
		if (weight == 0)
			continue;
		if ((og.compare(0, 2, "L_") == 0) ||
		    (og.compare(0, 4, "LESS")==0) ||
                (og.compare(0,2,">@") == 0))
		{
			lt_obs_names.push_back(oname);
		}
		else if ((og.compare(0, 2, "G_")==0) ||
		        (og.compare(0, 7, "GREATER")==0) ||
                (og.compare(0,2,"<@") == 0))
		{
			gt_obs_names.push_back(oname);
		}
        else
        {
            if ((extfile_gt.find(oname) != extfile_gt.end()) &&
                    (extfile_lt.find(oname) != extfile_lt.end()))
            {
                double_obs_bounds[oname] = pair<double,double>(extfile_gt.at(oname),extfile_lt.at(oname));
            }
            else if (extfile_gt.find(oname) != extfile_gt.end())
            {
                gt_obs_bounds[oname] = extfile_gt.at(oname);
            }
            else if (extfile_lt.find(oname) != extfile_lt.end())
            {
                lt_obs_bounds[oname] = extfile_lt.at(oname);
            }
        }
    }

	//save the org reg factor and org q vector
	org_reg_factor = pest_scenario->get_pestpp_options().get_ies_reg_factor();
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
        //prepare_csv(file_manager->open_ofile_ext(tag+"phi.noise.csv"), oreal_names);
		prepare_csv(file_manager->open_ofile_ext(tag+"phi.composite.csv"), oreal_names);
		prepare_csv(file_manager->open_ofile_ext(tag+"phi.regul.csv"), preal_names);
		prepare_group_csv(file_manager->open_ofile_ext(tag+"phi.group.csv"));
        prepare_lambda_csv(file_manager->open_ofile_ext(tag+"phi.lambda.csv"));

	}

}

Eigen::MatrixXd L2PhiHandler::get_obs_resid(ObservationEnsemble &oe, bool apply_ineq)
{
	vector<string> names = oe_base->get_var_names();
    Eigen::MatrixXd oe_vals = oe.get_eigen(vector<string>(),names);
	Eigen::MatrixXd resid = oe_vals -
		oe_base->get_eigen(oe.get_real_names(), vector<string>());
	
	if (apply_ineq)
		apply_ineq_constraints(resid,oe_vals, names);
	return resid;
}

map<string,map<string,double>> L2PhiHandler::get_swr_real_map(ObservationEnsemble& oe, ObservationEnsemble& weights,phiType ptype)
{
    Eigen::MatrixXd resid;
    if (ptype == phiType::MEAS)
    {
      resid = get_obs_resid(oe,true);
    }
    else if (ptype == phiType::ACTUAL)
    {
        resid = get_actual_obs_resid(oe);
    }
    else
    {
        throw runtime_error("L2PhiHandler::get_swr_real_map() unsupported ptype");
    }

    Eigen::MatrixXd wmat = weights.get_eigen(oe.get_real_names(),oe_base->get_var_names());
    resid = resid.array() * wmat.array();
    resid = resid.array() * resid.array();
    map<string,map<string,double>> swr_map;
    vector<string> var_names = oe_base->get_var_names();
    vector<string> real_names = oe.get_real_names();
    for (int i=0;i<real_names.size();i++) {
        swr_map[real_names[i]] = map<string,double>();
        for (int j = 0; j < var_names.size(); j++) {
            swr_map[real_names[i]][var_names[j]] = resid(i, j);
        }
    }
    return swr_map;
}

map<string,double> L2PhiHandler::get_swr_map(ObservationEnsemble& oe, string real_name,phiType ptype)
{

    Eigen::MatrixXd resid;
    if (ptype == phiType::ACTUAL)
    {
       resid = get_actual_obs_resid(oe);
    }
    else if (ptype == phiType::MEAS)
    {
        resid = get_obs_resid(oe,true);
    }
    else
    {
        throw runtime_error("L2PhiHandler::get_swr_map() unsupported ptype");
    }


    Eigen::VectorXd q = get_q_vector();
    for (int i=0;i<resid.rows();i++) {
        resid.row(i) = resid.row(i).array() * q.array().transpose();
        resid.row(i) = resid.row(i).array() * resid.row(i).array();
    }
    Eigen::VectorXd mean_vec;
    map<string,double> swr_map;
    vector<string> var_names = oe_base->get_var_names();
    if (real_name.size() > 0)
    {
        map<string,int> real_map = oe.get_real_map();
        if (real_map.find(real_name) != real_map.end())
        {
            mean_vec = resid.row(real_map[real_name]);
        }
        else
        {
            mean_vec = resid.colwise().mean();
        }
    }
    else
    {
        mean_vec = resid.colwise().mean();
    }

    for (int i=0;i<mean_vec.size();i++)
    {
        swr_map[var_names[i]] = mean_vec[i];
    }
    return swr_map;
}


Eigen::MatrixXd L2PhiHandler::get_obs_resid_subset(ObservationEnsemble &oe, bool apply_ineq, vector<string> real_names)
{
    if (real_names.size() == 0)
    {
        real_names = oe.get_real_names();
    }
	vector<string> names = oe.get_var_names();
    Eigen::MatrixXd oe_vals = oe.get_eigen(real_names,vector<string>());
	Eigen::MatrixXd resid = oe_vals - oe_base->get_eigen(real_names, names);
	if (apply_ineq)
		apply_ineq_constraints(resid, oe_vals,names);
	return resid;
}

Eigen::MatrixXd L2PhiHandler::get_par_resid(ParameterEnsemble &pe)
{
	Eigen::MatrixXd resid = pe.get_eigen(vector<string>(), pe_base->get_var_names()) -
		pe_base->get_eigen(pe.get_real_names(), vector<string>());
	return resid;
}

Eigen::MatrixXd L2PhiHandler::get_par_resid_subset(ParameterEnsemble &pe, vector<string> real_names)
{
    if (real_names.size() == 0)
    {
        real_names = pe.get_real_names();
    }
	Eigen::MatrixXd resid = pe.get_eigen(real_names,vector<string>()) - pe_base->get_eigen(real_names,pe.get_var_names());
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
	apply_ineq_constraints(resid, oe_vals, act_obs_names);
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
    ObservationInfo oinfo = pest_scenario->get_ctl_observation_info();
    num_conflict_group.clear();
    for (auto& group : pest_scenario->get_ctl_ordered_obs_group_names())
    {
        num_conflict_group[group] = 0;
    }
    vector<string> in_conflict = detect_simulation_data_conflict(oe,"");
    string group;
    for (auto& ic : in_conflict)
    {
        group = oinfo.get_group(ic);
        num_conflict_group[group]++;
    }
	//build up obs group and par group idx maps for group reporting
	obs_group_idx_map.clear();
	vector<string> nnz_obs = oe_base->get_var_names();

	vector<int> idx;
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

        for (auto& r : reg_map)
        {
            regul[r.first] = r.second.sum();
            par_group_phi_map[r.first] = get_par_group_contrib(r.second);
        }
        composite.clear();
        composite = calc_composite(meas, regul);
	}
	else
    {
	    composite = meas;
    }
	
	actual.clear();
	for (auto &pv : calc_actual(oe, q))
	{
		actual[pv.first] = pv.second.sum();
		obs_group_phi_map[pv.first] = get_obs_group_contrib(pv.second);
	}

}

void L2PhiHandler::update(ObservationEnsemble & oe, ParameterEnsemble & pe, ObservationEnsemble& weights)
{
    ObservationInfo oinfo = pest_scenario->get_ctl_observation_info();
    num_conflict_group.clear();
    for (auto& group : pest_scenario->get_ctl_ordered_obs_group_names())
    {
        num_conflict_group[group] = 0;
    }
    vector<string> in_conflict = detect_simulation_data_conflict(oe,"");
    string group;
    for (auto& ic : in_conflict)
    {
        group = oinfo.get_group(ic);
        num_conflict_group[group]++;
    }
    //build up obs group and par group idx maps for group reporting
    obs_group_idx_map.clear();
    vector<string> nnz_obs = oe_base->get_var_names();
    vector<int> idx;
    for (auto& og : pest_scenario->get_ctl_ordered_obs_group_names())
        obs_group_idx_map[og] = vector<int>();

    for (int i = 0; i < nnz_obs.size(); i++)
    {
        obs_group_idx_map[oinfo.get_group(nnz_obs[i])].push_back(i);
    }

    //update the various phi component vectors
    meas.clear();
    obs_group_phi_map.clear();
    map<string, Eigen::VectorXd> meas_map = calc_meas(oe, weights);
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

        for (auto& r : reg_map)
        {
            regul[r.first] = r.second.sum();
            par_group_phi_map[r.first] = get_par_group_contrib(r.second);
        }
        composite.clear();
        composite = calc_composite(meas, regul);
    }
    else
    {
        composite = meas;
    }

    actual.clear();
    for (auto &pv : calc_actual(oe, weights))
    {
        actual[pv.first] = pv.second.sum();
        obs_group_phi_map[pv.first] = get_obs_group_contrib(pv.second);
    }

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
    case L2PhiHandler::phiType::NOISE:
        return &noise;
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
    case L2PhiHandler::phiType::NOISE:
        return noise;
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

double L2PhiHandler::calc_median(const std::vector<double>& values) {
	size_t const size = values.size();
	if (size % 2 != 0) {
		// odd number of elements: median is the middle element
		return values[size / 2];
	} else {
		// even number of elements: median is the average of the two middle elements
		return (values[(size / 2) - 1] + values[size / 2]) / 2.0;
	}
}
double L2PhiHandler::calc_iqr_thresh(map<string, double> *phi_map, double bad_phi_sigma) {
	// get the phi map
	map<string, double>::iterator pi = phi_map->begin(), end = phi_map->end();

	// put the values into a vector
	std::vector<double> values;
	for (; pi != end; ++pi)
		values.push_back(pi->second);
	//for (const auto& pair : pi) {
	//	values.push_back(pair.second);
	//}
	// sort 'em
	std::sort(values.begin(), values.end());
	// calc the median depending on vector size
	double median = calc_median(values);

	// get the total size
	int const n = values.size();

	// q1: median of the lower half)
	std::vector<double> lower_half;
	for (int i = 0; i < n / 2; ++i) {
		lower_half.push_back(values[i]);
	}
	const double q1 = calc_median(lower_half);

	// q3: median of the upper half)
	std::vector<double> upper_half;
	// Handle odd/even cases for upper half starting point
	int upper_half_start_index = (n % 2 == 0) ? n / 2 : n / 2 + 1;
	for (int i = upper_half_start_index; i < n; ++i) {
		upper_half.push_back(values[i]);
	}
	const double q3 = calc_median(upper_half);

	// calculate iqr
	const double iqr =  q3 - q1;

	// return the threshold
	if (bad_phi_sigma<0)
		bad_phi_sigma = -1*bad_phi_sigma;
	return q3+bad_phi_sigma*iqr;
}


double L2PhiHandler::get_representative_phi(phiType pt)
{
    //if (pest_scenario->get_pestpp_options().get_ies_n_iter_reinflate() < 0)
    bool use_min = false;
//    for (auto& fac : pest_scenario->get_pestpp_options().get_ies_n_iter_reinflate())
//    {
//        if (fac < 0)
//        {
//            use_min = true;
//            break;
//        }
//    }
    if (use_min)
    {
        return get_min(pt);
    }

    else {
        return get_mean(pt);
    }
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
	map<string, double>* phi_map = get_phi_map_ptr(pt);
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
	double mn = numeric_limits<double>::max();//1.0e+30;
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
    case L2PhiHandler::phiType::NOISE:
        typ = "noise";
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

bool cmp_pair(pair<string,double>& first, pair<string,double>& second)
{
    return first.second > second.second;
}

vector<string> L2PhiHandler::detect_simulation_data_conflict(ObservationEnsemble& _oe, string csv_tag) {
    vector<string> in_conflict;

    ofstream pdccsv;
    if (csv_tag.size() > 0)
        pdccsv.open(file_manager->get_base_filename() + csv_tag);

    double smin, smax, omin, omax, smin_stat, smax_stat, omin_stat, omax_stat;
    map<string, int> smap, omap;
    vector<string> snames = _oe.get_var_names();
    vector<string> onames = oe_base->get_var_names();
    vector<string> temp = get_lt_obs_names();
    set<string> ineq_lt(temp.begin(), temp.end());
    //set<string>::iterator end = ineq.end();
    temp = get_gt_obs_names();
    set<string> ineq_gt(temp.begin(), temp.end());
    temp.resize(0);

    for (int i = 0; i < snames.size(); i++) {
        smap[snames[i]] = i;
    }
    for (int i = 0; i < onames.size(); i++) {
        omap[onames[i]] = i;
    }
    int sidx, oidx;
    bool use_stat_dist = true;
    if (pest_scenario->get_pestpp_options().get_ies_pdc_sigma_distance() <= 0.0)
        use_stat_dist = false;

    double smn, sstd, omn, ostd, dist;
    double sd = abs(pest_scenario->get_pestpp_options().get_ies_pdc_sigma_distance());
    int oe_nr = _oe.shape().first;
    int oe_base_nr = oe_base->shape().first;
    Eigen::VectorXd t;

    if (csv_tag.size() > 0)
    {
        pdccsv << "name,obs_mean,obs_std,obs_min,obs_max,obs_stat_min,obs_stat_max,sim_mean,sim_std,sim_min,sim_max,sim_stat_min,sim_stat_max,distance";
        pdccsv << endl;
    }

    for (auto oname: pest_scenario->get_ctl_ordered_nz_obs_names()) {
        //if (ineq.find(oname) != end)
        //	continue;
        sidx = smap[oname];
        oidx = omap[oname];
        smin = _oe.get_eigen_ptr()->col(sidx).minCoeff();
        omin = oe_base->get_eigen_ptr()->col(oidx).minCoeff();
        smax = _oe.get_eigen_ptr()->col(sidx).maxCoeff();
        omax = oe_base->get_eigen_ptr()->col(oidx).maxCoeff();
        t = _oe.get_eigen_ptr()->col(sidx);
        smn = t.mean();
        sstd = std::sqrt((t.array() - smn).square().sum() / (oe_nr - 1));
        smin_stat = smn - (sd * sstd);
        smax_stat = smn + (sd * sstd);
        t = oe_base->get_eigen_ptr()->col(oidx);
        omn = t.mean();
        ostd = std::sqrt((t.array() - omn).square().sum() / (oe_base_nr - 1));
        omin_stat = omn - (sd * ostd);
        omax_stat = omn + (sd * ostd);
        bool conflicted = false;
        if (use_stat_dist) {
            if (ineq_lt.find(oname) != ineq_lt.end()) {
                if (smin_stat > omax_stat)
                    conflicted = true;
            } else if (ineq_gt.find(oname) != ineq_gt.end()) {
                if (smax_stat < omin_stat)
                    conflicted = true;
            } else if ((smin_stat > omax_stat) || (smax_stat < omin_stat)) {
                conflicted = true;
            }
        } else {
            if (ineq_lt.find(oname) != ineq_lt.end()) {
                if (smin > omax)
                    conflicted = true;
            } else if (ineq_gt.find(oname) != ineq_gt.end()) {
                if (smax < omin)
                    conflicted = true;
            } else if ((smin > omax) || (smax < omin)) {
                conflicted = true;
            }
        }
        if (conflicted) {
            in_conflict.push_back(oname);
            if (csv_tag.size() > 0) {
                dist = max((smin - omax), (omin - smax));

                pdccsv << pest_utils::lower_cp(oname) << "," << omn << "," << ostd << "," << omin << "," << omax << ","
                       << omin_stat << ","
                       << omax_stat;
                pdccsv << "," << smn << "," << sstd << "," << smin << "," << smax << "," << smin_stat << ","
                       << smax_stat << "," << dist << endl;
            }
        }
    }
    pdccsv.close();
    return in_conflict;
}

void L2PhiHandler::report_group(bool echo) {

    map<string,double> mn_map;
    map<string,double> std_map;
    map<string,double> pmn_map;
    map<string,double> pstd_map;
    map<string,double> mx_map;
    map<string,double> mmn_map;
    map<string,double> pmx_map;
    map<string,double> pmmn_map;

    set<string> snzgroups;
    ObservationInfo* oi_ptr = pest_scenario->get_observation_info_ptr();
    for (auto& o : pest_scenario->get_ctl_ordered_nz_obs_names())
    {
        if (snzgroups.find(oi_ptr->get_group(o)) == snzgroups.end())
            snzgroups.emplace(oi_ptr->get_group(o));
    }

    if (snzgroups.size() == 0)
        return;

    double tot = 0, ptot = 0;
    double v = 0,pv = 0;
    int c = 0;
    for (auto& g : snzgroups)
    {
        tot = 0;
        ptot = 0;
        v = 0;
        pv = 0;
        c = 0;
        mx_map[g] = -1e+300;
        mmn_map[g] = 1e+300;
        pmx_map[g] = -1e+300;
        pmmn_map[g] = 1e+300;
        for (auto& o : obs_group_phi_map)
        {
            if (actual.find(o.first) == actual.end())
            {
                continue;
            }
            tot = tot + o.second[g];
            if (actual[o.first] > 0)
                ptot = ptot + (o.second[g]/actual[o.first]);
            mx_map[g] = max(mx_map[g],o.second[g]);
            mmn_map[g] = min(mmn_map[g],o.second[g]);
            pmx_map[g] = max(pmx_map[g],(o.second[g]/actual[o.first]));
            pmmn_map[g] = min(pmmn_map[g],(o.second[g]/actual[o.first]));
            c++;
        }
        tot = tot / (double)c;
        if (ptot > 0)
            ptot = ptot/(double)c;
        mn_map[g] = tot;
        pmn_map[g] = ptot;

        for (auto& o : obs_group_phi_map)
        {
            if (actual.find(o.first) == actual.end())
            {
                continue;
            }
            v = v + (pow(o.second[g] - tot,2));
            if (actual[o.first] > 0)
                pv = pv + (pow(ptot - (o.second[g]/actual[o.first]),2));
        }
        if (v != 0)
            std_map[g] = sqrt(v/(double)(c-1));
        else
            std_map[g] = 0.0;
        if (pv != 0)
            pstd_map[g] = sqrt(pv/(double)(c-1));
        else
            pstd_map[g] = 0.0;
    }
    ofstream& f = file_manager->rec_ofstream();
    int len = 1;
    for (auto& o : obs_group_phi_map)
    {
        for (auto& oo : o.second)
        {
            //if (oo.second == 0.0)
            if (snzgroups.find(oo.first) == snzgroups.end())
                continue;
            len = max(len,(int)oo.first.size());
        }
        break;
    }
    len++;

    string s;
    stringstream ss;
    ss.str("");
    ss << "  ---  observation group phi summary ---  " << endl;
    ss << "       (computed using 'actual' phi)" << endl;
    ss << "           (sorted by mean phi)" << endl;
    ss << left << setw(len) << "group" << right << setw(7) << "count" << setw(11) << "nconflict" << setw(10) << "mean" << setw(10) << "std";
    ss << setw(10) << "min" << setw(10) << "max";
    ss << setw(10) << "percent" << setw(10) << "std" << endl; //<< setw(10) << "min " << setw(10) << "max " << endl;
    f << ss.str();
    if (echo)
        cout << ss.str();

    vector<string> nzgroups(snzgroups.begin(),snzgroups.end());
    //sort(nzgroups.begin(),nzgroups.end());
    vector<pair<string,double>> pairs;
    for (auto& it : mn_map)
        pairs.push_back(it);


    sort(pairs.begin(),pairs.end(),cmp_pair);

    c = 0;
    int nzc = 0;
    string g;
    for (auto& pair : pairs)
    {
        g = pair.first;
        nzc = 0;
        for (auto& n : oi_ptr->observations)
        {
            if ((n.second.weight > 0) && (n.second.group == g))
            {
                nzc++;
            }
        }
        ss.str("");
        ss << left << setw(len) << pest_utils::lower_cp(g) << " ";
        ss << right << setw(6)  << setprecision(0) << nzc << " ";
        ss << right << setw(10) << num_conflict_group[g] << " ";

        ss << right << setw(9) << setprecision(3) << mn_map[g] << " ";
        ss << setw(9) << setprecision(3) << std_map[g] << " ";
        ss << setw(9) << setprecision(3) << mmn_map[g] << " ";
        ss << setw(9) << setprecision(3) << mx_map[g] << " ";

        ss << setw(9) << setprecision(3) << 100. * pmn_map[g] << " ";
        ss << setw(9) << setprecision(3) << 100. * pstd_map[g];
        //ss << setw(15) << 100. * pmmn_map[g];
        //ss << setw(15) << 100. * pmx_map[g];
        ss << endl;

        f << ss.str();
        if ((echo) && (c <= 10))
            cout << ss.str();
        c++;
    }
    f << "    Note: 'percent' is the percentage of the actual phi for each realization." << endl << endl;
    //if (echo)
    //    cout << "    Note: 'percent' is the percentage of the actual phi for each realization." << endl;
    if ((echo) && (c > 10))
        cout << "    Note: only the first 10 highest mean phi groups are listed. See rec file for full listing" << endl;
    if (echo)
        cout << endl;
}

void L2PhiHandler::report(bool echo, bool group_report)
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
	


    if (group_report)
    {
        report_group(echo);
    }

    f << endl << endl;
    f.flush();

}




void L2PhiHandler::write(int iter_num, int total_runs, bool write_group)
{
	write_csv(iter_num, total_runs, file_manager->get_ofstream(tag+"phi.actual.csv"), phiType::ACTUAL,oreal_names);
	write_csv(iter_num, total_runs, file_manager->get_ofstream(tag+"phi.meas.csv"), phiType::MEAS, oreal_names);
    //write_csv(iter_num, total_runs, file_manager->get_ofstream(tag+"phi.noise.csv"), phiType::NOISE, oreal_names);

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

void L2PhiHandler::write_lambda(int iteration,int num_reals,double current_lambda,double current_comp_mean_phi,
                                double current_comp_std_phi,double lambda_mult,
                                double lambda, double comp_mean_phi,double comp_std_phi)
{
    ofstream& csv = file_manager->get_ofstream(tag+"phi.lambda.csv");
    csv << iteration << ',' << num_reals << ',' << current_lambda << ',' << current_comp_mean_phi << ',',
    csv << current_comp_std_phi << ',' << lambda_mult << ',';
    csv << lambda << ',' << comp_mean_phi << ',' << comp_std_phi << endl;
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

void L2PhiHandler::prepare_lambda_csv(ofstream &csv)
{
    csv << "iteration,num_reals,current_lambda,current_comp_mean_phi,current_comp_std_phi,lambda_scale_fac,lambda,comp_mean_phi,comp_std_phi" << endl;
}



vector<int> L2PhiHandler::get_idxs_greater_than(double bad_phi, double bad_phi_sigma, ObservationEnsemble &oe, ObservationEnsemble& weights)
{
    //todo: handle weights ensemble here...
	map<string, double> _meas;
	//Eigen::VectorXd q = get_q_vector();
	//for (auto &pv : calc_meas(oe, q))
    for (auto &pv : calc_meas(oe, weights))
        _meas[pv.first] = pv.second.sum();
	double mean = calc_mean(&_meas);
	double std = calc_std(&_meas);
	vector<int> idxs;
	vector<string> names = oe.get_real_names();
	double bad_thres = std::numeric_limits<double>::max();
	if (bad_phi_sigma < std::numeric_limits<double>::max()) {
        if (bad_phi_sigma > 0) {
            bad_thres = min(std::numeric_limits<double>::max(), mean + (std * bad_phi_sigma));

            ofstream &frec = file_manager->rec_ofstream();
            frec << "...bad_phi_sigma " << bad_phi_sigma << " and mean " << mean << " yields bad phi threshold of "
                 << bad_thres << endl;
        } else {
            const double qval = -1. * bad_phi_sigma;
            bad_thres = min(std::numeric_limits<double>::max(),calc_iqr_thresh(&_meas, bad_phi_sigma));//Q3 + bad_phi_sigma * IQR
            ofstream &frec = file_manager->rec_ofstream();
            frec << "...bad_phi_sigma interquartile factor " << qval << " yields bad phi threshold of " << bad_thres << endl;
        	if (qval > 3.) {
        		frec << "WARNING:...bad_phi_sigma interquartile factor " << qval << " greater than 3.0 may be too large... " << endl;

        	}

        }
    }

	for (int i = 0; i < names.size(); i++) {

        if ((_meas[names[i]] > bad_phi) || (_meas[names[i]] > bad_thres)) {
            if (names[i] == BASE_REAL_NAME)
                cout << "...not dropping 'base' real even though phi is 'bad'" << endl;
            else
                idxs.push_back(i);
        }
    }
	return idxs;
}

void upgrade_thread_function_phi(int id, Eigen::MatrixXd& wmat, Eigen::MatrixXd& resid, vector<string> oe_real_names, map<string,map<string,double>>& phi_map, PhiThread& worker, exception_ptr& eptr) {
    try {
        worker.work(id, wmat, resid, oe_real_names, phi_map);
    }
    catch (...) {
        eptr = current_exception();
    }

    return;
}


PhiThread::PhiThread(vector<string> _oe_real_names)
: keys(_oe_real_names)
{
    count = 0;
    total = 0;

}

void PhiThread::work(int thread_id, Eigen::MatrixXd& wmat, Eigen::MatrixXd& resid, vector<string>& oe_real_names, map<string,map<string,double>>& phi_map) {
    stringstream ss;

    //these locks are used to control (thread-safe) access to the fast look up containers
    unique_lock<mutex> next_guard(next_lock, defer_lock);
    unique_lock<mutex> phi_map_guard(phi_map_lock, defer_lock);


    //This is the main thread loop - it continues until all upgrade pieces have been completed
    int idx;
    Eigen::MatrixXd rresid;
    map<string,double> pmap;
    string rname;
    Eigen::VectorXd phi;
    while (true) {

        while (true) {
            if (next_guard.try_lock()) {
                //if all the pieces have been completed, return
                if (count == keys.size()) {

                    return;
                }

                idx = count;
                count++;
                next_guard.unlock();
                break;
            }
        }

        rresid = resid.array().rowwise() * wmat.row(idx).array();
        rresid = rresid.array().cwiseProduct(rresid.array());
        phi = rresid.rowwise().sum();
        pmap.clear();
        for (int ii = 0; ii<resid.rows(); ii++)
        {
            pmap[oe_real_names[ii]] = phi(ii);
        }

        while (true) {
            if (phi_map_guard.try_lock()) {
                phi_map[oe_real_names[idx]] = pmap;
                phi_map_guard.unlock();
                break;
            }
        }

    }
}


map<string,map<string,double>> L2PhiHandler::get_meas_phi_weight_ensemble(ObservationEnsemble& oe, ObservationEnsemble& weights)
{
    assert(oe.shape().first == weights.shape().first);
    map<string,map<string,double>> phi_map;
    vector<string> base_real_names = oe_base->get_real_names(), oe_real_names = oe.get_real_names();
    Eigen::MatrixXd resid = get_obs_resid(oe);// this is (Ho - Hs)
    Eigen::MatrixXd wmat = weights.get_eigen(vector<string>(),oe.get_var_names());
    assert(oe_real_names.size() == resid.rows());
    int num_threads = pest_scenario->get_pestpp_options().get_ies_num_threads();
    stringstream ss;
    if ((num_threads > 1)) {

        Eigen::setNbThreads(1);
        vector<thread> threads;
        vector<exception_ptr> exception_ptrs;
        PhiThread *ut_ptr = new PhiThread(oe_real_names);
        for (int i = 0; i < num_threads; i++) {
            exception_ptrs.push_back(exception_ptr());
        }
        for (int i = 0; i < num_threads; i++) {
            threads.push_back(thread(upgrade_thread_function_phi, i, std::ref(wmat),std::ref(resid),oe_real_names,std::ref(phi_map),std::ref(*ut_ptr), std::ref(exception_ptrs[i])));
        }

        int num_exp = 0;

        for (int i = 0; i < num_threads; ++i) {
            bool found = false;
            if (exception_ptrs[i]) {
                found = true;
                num_exp++;
                try {
                    rethrow_exception(exception_ptrs[i]);
                }
                catch (const std::exception &e) {
                    //ss.str("");
                    ss << " thread " << i << "raised an exception: " << e.what();
                    //throw runtime_error(ss.str());
                }
                catch (...) {
                    //ss.str("");
                    ss << " thread " << i << "raised an exception";
                    //throw runtime_error(ss.str());
                }
            }
            threads[i].join();
            if ((exception_ptrs[i]) && (!found)) {
                num_exp++;
                try {
                    rethrow_exception(exception_ptrs[i]);
                }
                catch (const std::exception &e) {
                    ss << " thread " << i << "raised an exception: " << e.what();

                }
                catch (...) {
                    ss << " thread " << i << "raised an exception: ";
                }
            }
        }
        if (num_exp > 0) {
            throw runtime_error(ss.str());
        }
        delete ut_ptr;
        return phi_map;
    }

    set<string> bset(base_real_names.begin(),base_real_names.end());
    set<string>::iterator end = bset.end();
    string rname;

    map<string,double> pmap;
    Eigen::MatrixXd rresid;
    for (int i=0;i<oe.shape().first;i++)
    {
        rresid = resid.array().rowwise() * wmat.row(i).array();
        rresid = rresid.array().cwiseProduct(rresid.array());
        pmap.clear();
        for (int ii = 0; ii<resid.rows(); ii++)
        {
            rname = oe_real_names[ii];

            pmap[rname] = rresid.row(ii).sum();
        }
        phi_map[oe_real_names[i]] = pmap;
    }
    return phi_map;
}




map<string, Eigen::VectorXd> L2PhiHandler::calc_meas(ObservationEnsemble & oe, Eigen::VectorXd &q_vec)
{
	map<string, Eigen::VectorXd> phi_map;
	Eigen::VectorXd oe_base_vec, oe_vec, diff, w_vec;
	vector<string> act_obs_names = pest_scenario->get_ctl_ordered_nz_obs_names();
	vector<string> base_real_names = oe_base->get_real_names(), oe_real_names = oe.get_real_names();
	set<string> bset(base_real_names.begin(),base_real_names.end());
	set<string>::iterator end = bset.end();
	string rname;
	if (act_obs_names.size() == 0)
	{
		for (auto name : oe.get_real_names())
			phi_map[name] = Eigen::VectorXd();
		return phi_map;
	}
	Eigen::MatrixXd resid = get_obs_resid(oe);// this is (Ho - Hs)
    assert(oe_real_names.size() == resid.rows());
    resid = resid.array().rowwise() * q_vec.transpose().array();
    resid = resid.array().cwiseProduct(resid.array());
	for (int i = 0; i<resid.rows(); i++)
	{
		rname = oe_real_names[i];
		if (bset.find(rname) == end)
			continue;
	    phi_map[rname] = resid.row(i);
	}
	return phi_map;
}


map<string, Eigen::VectorXd> L2PhiHandler::calc_meas(ObservationEnsemble & oe, ObservationEnsemble& weights)
{
    map<string, Eigen::VectorXd> phi_map;
    Eigen::VectorXd oe_base_vec, oe_vec, diff, w_vec;
    vector<string> act_obs_names = pest_scenario->get_ctl_ordered_nz_obs_names();

    if (act_obs_names.size() == 0)
    {
        for (auto name : oe.get_real_names())
            phi_map[name] = Eigen::VectorXd();
        return phi_map;
    }
    vector<string> base_real_names = oe_base->get_real_names(), oe_real_names = oe.get_real_names();
    set<string> bset(base_real_names.begin(),base_real_names.end());
    set<string>::iterator end = bset.end();
    string rname;
    Eigen::MatrixXd resid = get_obs_resid(oe);// this is (Ho - Hs)
    Eigen::MatrixXd wmat = weights.get_eigen(oe.get_real_names(),oe_base->get_var_names());
    assert(oe_real_names.size() == resid.rows());
    resid = resid.array() * wmat.array();
    resid = resid.array().cwiseProduct(resid.array());
    for (int i = 0; i<resid.rows(); i++)
    {
        rname = oe_real_names[i];
        if (bset.find(rname) == end)
            continue;
        phi_map[rname] = resid.row(i);
    }
    return phi_map;
}


map<string, Eigen::VectorXd> L2PhiHandler::calc_regul(ParameterEnsemble & pe)
{
	map<string, Eigen::VectorXd> phi_map;
	vector<string> real_names = pe.get_real_names();
	pe_base->transform_ip(ParameterEnsemble::transStatus::NUM);
	pe.transform_ip(ParameterEnsemble::transStatus::NUM);
	Eigen::MatrixXd resid = get_par_resid(pe);

	resid = resid.array().rowwise() * parcov_inv_diag.transpose().array();
    resid = resid.array().cwiseProduct(resid.array());

	for (int i = 0; i < real_names.size(); i++)
	{
        phi_map[real_names[i]] = resid.row(i);
	}
	return phi_map;
}


vector<string> L2PhiHandler::get_violating_realizations(ObservationEnsemble& oe, const vector<string>& viol_obs_names)
{
    stringstream ss;
    vector<string> viol_real_names;
    if (viol_obs_names.size() == 0)
        return viol_real_names;

    Eigen::MatrixXd resid = get_actual_obs_resid(oe);
    map<string,int> vmap;
    vector<string> nz_onames = oe_base->get_var_names();
    for (int i=0;i<nz_onames.size();i++)
        vmap[nz_onames[i]] = i;
    map<string,int>::iterator end = vmap.end();
    vector<string> missing;
    ObservationInfo* oi = pest_scenario->get_observation_info_ptr();
    vector<string> nz_viol_obs_names;
    map<string,int> nz_viol_vmap;
    for (auto& name : viol_obs_names)
    {
        if (vmap.find(name) == end)
            missing.push_back(name);
        else
        {
            if (oi->get_weight(name) > 0)
            {
                nz_viol_obs_names.push_back(name);
                nz_viol_vmap[name] = vmap.at(name);
            }
        }
    }

    Eigen::VectorXd real;
    double sum;
    vector<string> rnames = oe.get_real_names();
    for (int i=0;i<resid.rows();i++)
    {
        real = resid.row(i).cwiseAbs();
        sum = 0.0;
        for (auto& v : nz_viol_vmap)
        {
            sum += real(v.second);
        }
        if (sum > 1.0e-7)
        {
            viol_real_names.push_back(rnames[i]);
        }
    }
    return viol_real_names;
}


void L2PhiHandler::apply_ineq_constraints(Eigen::MatrixXd &resid, Eigen::MatrixXd &sim_vals,vector<string> &names)
{
	if ((lt_obs_names.empty()) && (gt_obs_names.empty()) &&
            (lt_obs_bounds.empty()) && (gt_obs_bounds.empty()) &&
            (double_obs_bounds.empty()))
	    return;
	assert(names.size() == resid.cols());

	map<string, double> lt_vals,gt_vals;
	Observations obs = pest_scenario->get_ctl_observations();
	for (auto &n : lt_obs_names)
		lt_vals[n] = obs.get_rec(n);
	for (auto &n : gt_obs_names)
		gt_vals[n] = obs.get_rec(n);
	if ((lt_vals.empty()) && (gt_vals.empty()) &&
       (lt_obs_bounds.empty()) && (gt_obs_bounds.empty()) &&
            (double_obs_bounds.empty()))
		return;
	map<string, int> idxs;
	for (int i = 0; i < names.size(); i++)
		idxs[names[i]] = i;
	int idx;
	double val,val2;
	Eigen::VectorXd col, scol;

    for (auto const iv : double_obs_bounds)
    {
        idx = idxs[iv.first];
        col = resid.col(idx);
        scol = sim_vals.col(idx);
        val = iv.second.first;
        val2 = iv.second.second;

        for (int i = 0; i < resid.rows(); i++)
            col(i) = ((scol(i) > val) && (scol(i) < val2)) ? 0.0 : col(i);
        resid.col(idx) = col;
    }

	for (auto const iv : lt_vals)
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


    for (auto const iv : lt_obs_bounds)
    {
        idx = idxs[iv.first];
        col = resid.col(idx);
        scol = sim_vals.col(idx);
        val = iv.second;
        //cout << resid.col(idx) << endl;
        for (int i = 0; i < resid.rows(); i++)
            col(i) = (scol(i) < val) ? 0.0 : col(i);
        //cout << resid.col(idx) << endl;
        //cout << col << endl << endl;
        resid.col(idx) = col;
        //cout << resid.col(idx) << endl;
    }


    //Eigen::MatrixXd temp = resid;
	for (auto const iv : gt_vals)
	{
		idx = idxs[iv.first];
		col = resid.col(idx);
		val = iv.second;
		for (int i = 0; i < resid.rows(); i++)
			col(i) = (col(i) > 0.0) ? 0.0 : col(i);
		resid.col(idx) = col;
	}

    for (auto const iv : gt_obs_bounds)
    {
        idx = idxs[iv.first];
        col = resid.col(idx);
        scol = sim_vals.col(idx);
        val = iv.second;
        for (int i = 0; i < resid.rows(); i++)
            col(i) = (scol(i) > val) ? 0.0 : col(i);
        resid.col(idx) = col;
    }
}



map<string, Eigen::VectorXd> L2PhiHandler::calc_actual(ObservationEnsemble & oe, Eigen::VectorXd &q_vec)
{
	map<string, Eigen::VectorXd> phi_map;
	Eigen::MatrixXd resid = get_actual_obs_resid(oe);
	vector<string> base_real_names = oe_base->get_real_names(), oe_real_names = oe.get_real_names();
	set<string> bset(base_real_names.begin(),base_real_names.end());
	set<string>::iterator end = bset.end();
	string rname;

	resid = resid.array().rowwise() * q_vec.transpose().array();
    resid = resid.array().cwiseProduct(resid.array());
	for (int i = 0; i<oe.shape().first; i++)
	{
		rname = oe_real_names[i];
		if (bset.find(rname) == end)
			continue;
		phi_map[rname] = resid.row(i);
	}
	return phi_map;
}

map<string, Eigen::VectorXd> L2PhiHandler::calc_actual(ObservationEnsemble & oe, ObservationEnsemble& weights)
{
    map<string, Eigen::VectorXd> phi_map;
    Eigen::MatrixXd resid = get_actual_obs_resid(oe);
    Eigen::MatrixXd wmat = weights.get_eigen(oe.get_real_names(),oe_base->get_var_names());
    vector<string> base_real_names = oe_base->get_real_names(), oe_real_names = oe.get_real_names();
    set<string> bset(base_real_names.begin(),base_real_names.end());
    set<string>::iterator end = bset.end();
    string rname;

    resid = resid.array() * wmat.array();
    resid = resid.array().cwiseProduct(resid.array());
    for (int i = 0; i<oe.shape().first; i++)
    {
        rname = oe_real_names[i];
        if (bset.find(rname) == end)
            continue;
        phi_map[rname] = resid.row(i);
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
			mea = _meas.at(orn);
			reg = _regul.at(prn);
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


void ParChangeSummarizer::summarize(ParameterEnsemble &pe, string filename)
{
	update(pe);
    map<string,int> excess_std_reduction = get_npar_per_group_with_excess_std_reduction(pe);
	vector<string> grp_names;// = base_pe_ptr->get_pest_scenario().get_ctl_ordered_par_group_names();
	vector<pair<double, string>> pairs;
    //for (auto m : mean_change)
	for (auto& m : percent_at_lbound)
    {
		//mean_pairs.push_back(pair<double, string>(abs(m.second), m.first));
        pairs.push_back(pair<int,string>(m.second+percent_at_ubound.at(m.first),m.first));
	}
	sort(pairs.begin(), pairs.end());
	//for (auto m : mean_pairs)
	for (int i = pairs.size() - 1; i >= 0; i--)
	{
		grp_names.push_back(pairs[i].second);
	}
	int mxlen = 7;
	for (auto& g : grp_names)
		mxlen = max(mxlen, (int)g.size());

	stringstream ss;
	ofstream &frec = file_manager_ptr->rec_ofstream();
	ss << endl << "   ---  parameter group change summary  ---    " << endl;
	cout << ss.str();
	frec << ss.str();
	ss.str("");
	ss << setw(mxlen) << left << "group" << right << setw(6) << "count" << setw(10) << right << "mean chg" << setw(9) << "std chg";
	ss << setw(11) << "n at ubnd" << setw(11) << "% at ubnd";
    ss << setw(11) << "n at lbnd" << setw(11) << "% at lbnd";
	ss << setw(12) << "n std decr" << endl;
	cout << ss.str();
	frec << ss.str();	
	
	
	int i = 0;
	int num_out;
	double percent_out;
	for (auto &grp_name : grp_names)
	{
		double mean_diff = mean_change[grp_name];
		double std_diff = std_change[grp_name];
        if (pargp2par_map[grp_name].size() == 0)
            continue;
		ss.str("");
		ss << setw(mxlen) << left << pest_utils::lower_cp(grp_name) << right << setw(6) << pargp2par_map[grp_name].size();
        ss << setw(10) << setprecision(4) << right << mean_diff * 100.0;
		ss << setw(9) << setprecision(4) << std_diff * 100.0;
        num_out = num_at_ubound[grp_name];
        percent_out = percent_at_ubound[grp_name];
		ss << setw(11) << num_out << setw(11) << setprecision(4) << percent_out;
        num_out = num_at_lbound[grp_name];
        percent_out = percent_at_lbound[grp_name];
        ss << setw(11) << num_out << setw(11) << setprecision(4) << percent_out;
		//ss << setw(10) << setprecision(4) << init_cv[grp_name] << setw(10) << curr_cv[grp_name] << setw(10) << setprecision(4) << endl;
		ss << setw(12) << excess_std_reduction[grp_name] << endl;
        if ((i < 15) && (pargp2par_map[grp_name].size() > 0)) {
            cout << ss.str();
            i++;
        }
		frec << ss.str();
	}

	ss.str("");
	ss << "    Note: parameter change summary sorted according to percent at bounds." << endl;
	//ss << "          'n CV decr' is the number of parameters with current CV less " << cv_dec_threshold*100.0 << "% of the initial CV" << endl;
	ss << "    Note: the parameter change statistics implicitly include the effect of " << endl;
	ss << "          realizations that have failed or have been dropped." << endl;
    ss << "    Note: the 'n std decr' is the number of parameters with current" << endl;
    ss << "          std less 5% of their initial std." << endl;

	//cout << ss.str();
	frec << ss.str();
	if (grp_names.size() > 15)
	{
		ss.str("");
		ss << "    Note: Only the first 15 adjustable parameter groups shown, see rec file for full listing" << endl;
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

	f << "group,mean_change,std_change,num_at_near_lbound,percent_at_near_lbound,num_at_near_ubound,percent_at_near_ubound,initial_cv,current_cv" << endl;
	f << setprecision(20);
	for (auto grp_name : base_pe_ptr->get_pest_scenario_ptr()->get_ctl_ordered_par_group_names())
	{
		f << pest_utils::lower_cp(grp_name) << "," << mean_change[grp_name]*100.0 << "," << std_change[grp_name]*100.0 << ",";
		f << num_at_lbound[grp_name] << "," << percent_at_lbound[grp_name];
        f << "," << num_at_ubound[grp_name] << "," << percent_at_ubound[grp_name];
		f <<","<< init_cv[grp_name] << "," << curr_cv[grp_name] <<  endl;
	}
	f.close();
	file_manager_ptr->rec_ofstream() << "...saved parameter change summary to " << filename << endl;
	cout << "...saved parameter change summary to " << filename << endl;

}

map<string,int> ParChangeSummarizer::get_npar_per_group_with_excess_std_reduction(ParameterEnsemble& _pe, double thresh)
{
    stringstream ss;
    map<string,double> pr_mn,pr_std;
    base_pe_ptr->fill_moment_maps(pr_mn,pr_std);
    map<string,double> pt_mn,pt_std;
    _pe.fill_moment_maps(pt_mn,pt_std);

    double ratio;
    map<string,int> results;

    ParameterInfo* pi = _pe.get_pest_scenario_ptr()->get_ctl_parameter_info_ptr_4_mod();
    string group;
    for (auto& pname : _pe.get_pest_scenario_ptr()->get_ctl_parameters())
    {
        group = pi->get_parameter_rec_ptr(pname.first)->group;
        if (results.find(group) == results.end())
        {
            results[group] = 0;
        }
    }

    double demon = 0.0;
    for (auto& std : pt_std)
    {
        demon = pr_std.at(std.first);
        if (demon <= 0.0)
        {
            ratio = 1.0;
        }
        else {
            ratio = std.second / demon;
        }
        if (!isnormal(ratio) && (ratio != 0.0))
        {
            ss.str("");
            ss << "L2PhiHandler: posterior-prior std ratio not normal for par '" << std.first << "': " << ratio;
            throw runtime_error(ss.str());
        }
        ratio = 1.0 - ratio;
        group = pi->get_parameter_rec_ptr(std.first)->group;

        if (ratio >= thresh)
        {
            results[group]++;
        }

    }
    return results;


}


void ParChangeSummarizer:: update(ParameterEnsemble& pe)
{
	mean_change.clear();
	std_change.clear();
	num_at_ubound.clear();
	percent_at_ubound.clear();
    num_at_lbound.clear();
    percent_at_lbound.clear();
	init_cv.clear();
	curr_cv.clear();
    map<string, double> mean_map, std_map;
	base_pe_ptr->fill_moment_maps(mean_map, std_map);
	init_moments = pair<map<string, double>, map<string, double>>(mean_map, std_map);
	mean_map.clear();
	std_map.clear();
	pe.fill_moment_maps(mean_map, std_map);
	double mean_diff = 0.0, std_diff = 0.0;
	double icv = 0.0, ccv = 0.0;
	double iicv = 0.0, cccv = 0.0;

	int ndec_cv = 0, dsize=0;
	double  value1, value2, v;
	vector<string> pnames = pe.get_var_names();
	Parameters lb = pe.get_pest_scenario_ptr()->get_ctl_parameter_info().get_low_bnd(pnames);
	pe.get_pest_scenario_ptr()->get_base_par_tran_seq().active_ctl2numeric_ip(lb);
	Parameters ub = pe.get_pest_scenario_ptr()->get_ctl_parameter_info().get_up_bnd(pnames);
	pe.get_pest_scenario_ptr()->get_base_par_tran_seq().active_ctl2numeric_ip(ub);
	//vector<string> grp_names = base_pe_ptr->get_pest_scenario().get_ctl_ordered_par_group_names();
	vector<string> grp_names = pe.get_pest_scenario_ptr()->get_ctl_ordered_par_group_names();
	map<string, int> idx_map;
	for (int i = 0; i < pnames.size(); i++)
		idx_map[pnames[i]] = i;
	int num_out_u,num_out_l, num_pars;
	int num_reals = pe.get_real_names().size();
	Eigen::ArrayXd arr;
	int mean_size, std_size,icv_size,ccv_size;
	for (auto& grp_name : grp_names)
	{
		mean_diff = 0.0, std_diff = 0.0;
		icv = 0.0; ccv = 0.0;
		icv_size = 0; ccv_size = 0;
		ndec_cv = 0;
		num_pars = pargp2par_map[grp_name].size();
		num_out_u = 0;
		num_out_l = 0;
		mean_size = 0;
		std_size = 0;
		for (auto& par_name : pargp2par_map[grp_name])
		{
			arr = pe.get_eigen_ptr()->col(idx_map[par_name]).array();
			for (int i = 0; i < num_reals; i++)
			{
				v = arr[i];
				if (v > (ub[par_name] - (abs(ub[par_name]) * 0.01)))
					num_out_u++;
				else if (v < (lb[par_name] + (abs(lb[par_name]) * 0.01)))
				    num_out_l++;
			}
			value1 = init_moments.first[par_name];
			value2 = value1 - mean_map[par_name];
			//if ((value1 != 0.0) && (value2 != 0.0))
			if (value1 != 0.0)
			{
				mean_diff += abs(value2) / abs(value1);
				mean_size++;
			}
			value1 = init_moments.second[par_name];
			value2 = value1 - std_map[par_name];
			if (value1 != 0.0)
			{
				std_diff += abs(value2) / value1;
				std_size++;
			}
			//cv calcs
			value1 = init_moments.first[par_name];
			value2 = init_moments.second[par_name];
			iicv = 0.0;
			if (value1 != 0.0)
			{
				iicv = value2 / abs(value1);
				icv += iicv;
				icv_size++;
			}
			value1 = mean_map[par_name];
			value2 = std_map[par_name];
			cccv = 0.0;
			if (value1 != 0.0)
			{
				cccv = value2 / abs(value1);
				ccv += cccv;
				ccv_size++;
			}
			if ((cccv != 0.0) && (iicv != 0.0))
			{
				value1 = cccv / iicv;
				if (value1 < cv_dec_threshold)
					ndec_cv++;
			}
		}
		
		if (mean_diff != 0.0)
		{	
			mean_diff = mean_diff / double(mean_size);
		}
		if (std_diff != 0.0)
		{
			std_diff = std_diff / double(std_size);
		}
		if (icv != 0.0)
		{
			icv = icv / double(icv_size);
		}
		if (ccv != 0.0)
		{
			ccv = ccv / double(ccv_size);
		}
			

		double percent_out = 0;


		mean_change[grp_name] = mean_diff;
		std_change[grp_name] = std_diff;
		num_at_lbound[grp_name] = num_out_l;
        num_at_ubound[grp_name] = num_out_u;
        percent_out = 0.0;
		if (num_pars > 0)
            percent_out = double(num_out_l) / double(num_pars * num_reals) * 100.;
		percent_at_lbound[grp_name] = percent_out;
        percent_out = 0.0;
        if (num_pars > 0)
            percent_out = double(num_out_u) / double(num_pars * num_reals) * 100.;
        percent_at_ubound[grp_name] = percent_out;
        curr_cv[grp_name] = ccv;
		init_cv[grp_name] = icv;
	}
}


pair<Parameters, Observations> save_real_par_rei(Pest& pest_scenario, ParameterEnsemble& pe, ObservationEnsemble& oe,
    OutputFileWriter& output_file_writer, FileManager& file_manager,
    int iter, string tag, int cycle,map<string,double> base_weights)
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
		if (pe.get_trans_status() == ParameterEnsemble::transStatus::NUM) {
            pts.numeric2ctl_ip(pars);
        }
        vector<string> frnames = pe.get_fixed_info().get_real_names();
        if (find(frnames.begin(),frnames.end(),tag) != frnames.end()) {
            map<string, double> fmap = pe.get_fixed_info().get_real_fixed_values(tag);

            for (auto &item: fmap) {
                pars.update_rec(item.first, item.second);
            }
        }
		if (cycle != NetPackage::NULL_DA_CYCLE)
			ss << cycle << ".";
		if (iter >= 0)
			ss << iter << ".";
		
		ss << pest_utils::lower_cp(tag) << ".par";
		output_file_writer.write_par(file_manager.open_ofile_ext(ss.str()), pars, *(pts.get_offset_ptr()),
			*(pts.get_scale_ptr()));
		//file_manager.close_file("par");
		file_manager.close_all_files_ending_with("par");
		vmap = oe.get_real_map();
		if (vmap.find(tag) == vmap.end())
		{
			//message(2, "unable to find 'BASE' realization in obs ensemble for saving .base.rei file, continuing...");
		}
		else
		{
			Observations obs = pest_scenario.get_ctl_observations();
			obs.update_without_clear(oe.get_var_names(), eigenvec_2_stlvec(oe.get_real_vector(tag)));
			ObjectiveFunc obj_func(&(pest_scenario.get_ctl_observations()), &(pest_scenario.get_ctl_observation_info()), &(pest_scenario.get_prior_info()));
			// save new residuals to .rei file
			ss.str("");
			if (cycle != NetPackage::NULL_DA_CYCLE)
				ss << cycle << ".";
			if (iter >= 0)
				ss << iter << ".";
			ss << pest_utils::lower_cp(tag) <<  ".rei";
			output_file_writer.write_rei(file_manager.open_ofile_ext(ss.str()), iter,
				pest_scenario.get_ctl_observations(), obs, obj_func, pars,base_weights);
			file_manager.close_all_files_ending_with("rei");
		}
		cout << "saved par and rei files for realization " << tag;
		if (iter >= 0)
		    cout << " for iteration " << iter;
		cout << endl;

	}
	
	return pair<Parameters, Observations>(pars, obs);

}


vector<int> run_ensemble_util(PerformanceLog* performance_log, ofstream& frec,ParameterEnsemble& _pe, ObservationEnsemble& _oe, 
	RunManagerAbstract* run_mgr_ptr, bool check_pe_consistency, const vector<int>& real_idxs, int da_cycle, string additional_tag)
{
	stringstream ss;
	ss << "queuing " << _pe.shape().first << " runs";
	performance_log->log_event(ss.str());
	run_mgr_ptr->reinitialize();
	map<int, int> real_run_ids;
	try
	{
		real_run_ids = _pe.add_runs(run_mgr_ptr, real_idxs,da_cycle,additional_tag);
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
		performance_log->log_event(ss.str());
		throw runtime_error(ss.str());
	}
	catch (...)
	{
        performance_log->log_event("error running ensemble");
		throw runtime_error(string("error running ensemble"));
	}

	performance_log->log_event("processing runs");
	if (real_idxs.size() > 0)
	{
		_oe.keep_rows(real_idxs);
	}
	vector<int> failed_real_indices;
	ParameterEnsemble run_mgr_pe = _pe.zeros_like(0);
	try
	{
		failed_real_indices = _oe.update_from_runs(real_run_ids, run_mgr_ptr, run_mgr_pe);
	} 
	catch (const exception& e)
	{
		ss.str("");
		ss << "error processing runs: " << e.what();
        performance_log->log_event(ss.str());
		throw runtime_error(ss.str());
	}
	catch (...)
	{
        performance_log->log_event("error processing runs");
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
			cout << "parameter ensemble consistency check failed for " << failing_reals.size() << ", see .rec file for listing" << endl;
			frec << "ERROR: the following realizations failed consistency check:" << endl;
			for (auto fr : failing_reals)
				frec << fr << ",";
			frec << endl;
			throw runtime_error("parameter ensemble consistency check failed, see .rec file");
			
		}
	}
	return failed_real_indices;
}

EnsembleMethod::EnsembleMethod(Pest& _pest_scenario, FileManager& _file_manager,
	OutputFileWriter& _output_file_writer, PerformanceLog* _performance_log,
	RunManagerAbstract* _run_mgr_ptr, string _alg_tag) : alg_tag(_alg_tag), pest_scenario(_pest_scenario),
	file_manager(_file_manager),output_file_writer(_output_file_writer), performance_log(_performance_log),
	run_mgr_ptr(_run_mgr_ptr), pe(&_pest_scenario), verbose_level(1)
{
	rand_gen = std::mt19937(pest_scenario.get_pestpp_options().get_random_seed());

	subset_rand_gen = std::mt19937(pest_scenario.get_pestpp_options().get_random_seed());
	pe.set_pest_scenario(&pest_scenario);
	oe.set_pest_scenario_ptr(&pest_scenario);
	pe.set_rand_gen(&rand_gen);
	oe.set_rand_gen(&rand_gen);
	oe_base.set_pest_scenario_ptr(&pest_scenario);
	oe_base.set_rand_gen(&rand_gen);
    weights.set_pest_scenario_ptr(&pest_scenario);
    weights.set_rand_gen(&rand_gen);

    localizer.set_pest_scenario(&pest_scenario);
	pp_args = pest_scenario.get_pestpp_options().get_passed_args();
	act_obs_names = pest_scenario.get_ctl_ordered_nz_obs_names();
	act_par_names = pest_scenario.get_ctl_ordered_adj_par_names();
	iter = 0;
    warn_min_reals = 10;
    error_min_reals = 2;
    org_obs_info = pest_scenario.get_ctl_observation_info_copy();

}

void EnsembleMethod::sanity_checks()
{
    PestppOptions* ppo = pest_scenario.get_pestpp_options_ptr();
    vector<string> errors;
    vector<string> warnings;
    stringstream ss;
    string par_csv = ppo->get_ies_par_csv();
    string obs_csv = ppo->get_ies_obs_csv();
    string restart_obs = ppo->get_ies_obs_restart_csv();
    string restart_par = ppo->get_ies_par_restart_csv();

    if (pest_scenario.get_pestpp_options().get_ies_use_mda() && (pest_scenario.get_pestpp_options().get_ies_loc_type()[0] == 'C'))
    {
        errors.push_back("Covariance-based localization not supported with MDA solver");
    }

//    if (pest_scenario.get_control_info().noptmax > 10)
//	{
//		warnings.push_back("noptmax > 10, don't expect anything meaningful from the results!");
//	}

	else if (pest_scenario.get_control_info().noptmax > 15)
	{
		warnings.push_back("noptmax > 15, this is a lot of iterations for an ensemble method...");
	}

    if (pest_scenario.get_control_info().pestmode == ControlInfo::PestMode::REGUL)
    {
        warnings.push_back("'pestmode' == 'regularization', in pestpp-ies, this is controlled with the ++ies_reg_factor argument, resetting to 'estimation'");
        //throw_ies_error("'pestmode' == 'regularization', please reset to 'estimation'");
    }
    else if (pest_scenario.get_control_info().pestmode == ControlInfo::PestMode::UNKNOWN)
    {
        warnings.push_back("unrecognized 'pestmode', using 'estimation'");
    }


    if (pest_scenario.get_ctl_ordered_pi_names().size() > 0)
    {
        warnings.push_back("prior information equations not supported in ensemble methods, ignoring...");
    }
    double acc_phi = ppo->get_ies_accept_phi_fac();
    if (acc_phi < 1.0)
        warnings.push_back("ies_accept_phi_fac < 1.0, not good!");
    if (acc_phi > 10.0)
        warnings.push_back("ies_accept_phi_fac > 10.0, this is prob too big, typical values 1.05 to 1.3");

    double lam_inc = ppo->get_ies_lambda_inc_fac();
    if (lam_inc < 1.0)
        errors.push_back("ies_lambda_inc_fac < 1.0, nope! how can lambda increase if this is less than 1.0???");

    double lam_dec = ppo->get_ies_lambda_dec_fac();
    if (lam_dec > 1.0)
        errors.push_back("ies_lambda_dec_fac > 1.0, nope!  how can lambda decrease if this is greater than 1.0???");

    if ((ppo->get_ies_par_csv().size() == 0) && (ppo->get_ies_use_empirical_prior()))
    {
        warnings.push_back("no point in using an empirical prior if we are drawing the par ensemble...resetting ies_use_empirical_prior to false");
        ppo->set_ies_use_empirical_prior(false);
    }
    if ((par_csv.size() == 0) && (restart_par.size() > 0))
        errors.push_back("ies_par_en is empty but ies_restart_par_en is not - how can this work?");
    if ((restart_par.size() > 0) && (restart_obs.size() == 0))
        errors.push_back("use of ies_restart_par_en requires ies_restart_obs_en");
    if ((par_csv.size() == 0) && (restart_obs.size() > 0))
        errors.push_back("ies_par_en is empty but ies_restart_obs_en is not - how can this work?");
    if ((obs_csv.size() == 0) && (ppo->get_ies_weights_csv().size() > 0))
        errors.push_back("ies_obs_en is empty but ies_weight_en is not");
    if (ppo->get_ies_bad_phi() <= 0.0)
        errors.push_back("ies_bad_phi <= 0.0, really?");
    if ((ppo->get_ies_num_reals() < error_min_reals) && (par_csv.size() == 0))
    {
        ss.str("");
        ss << "ies_num_reals < " << error_min_reals << ", this is redic!";
        warnings.push_back(ss.str());
        //ppo->set_ies_num_reals(warn_min_reals);
    }
    if ((ppo->get_ies_num_reals() < warn_min_reals) && (par_csv.size() == 0))
    {
        ss.str("");
        ss << "ies_num_reals < " << warn_min_reals << ", this is prob too few";
        warnings.push_back(ss.str());
    }
    //if (ppo->get_ies_reg_factor() < 0.0)
    //    errors.push_back("ies_reg_factor < 0.0 - WRONG!");
    //if (ppo->get_ies_reg_factor() > 1.0)
    //	errors.push_back("ies_reg_factor > 1.0 - nope");
    if ((par_csv.size() == 0) && (ppo->get_ies_subset_size() < 10000000) && (ppo->get_ies_num_reals() < ppo->get_ies_subset_size() * 2))
        warnings.push_back("ies_num_reals < 2*ies_subset_size: you not gaining that much using subset here");
    //if ((ppo->get_ies_subset_size() < 100000001) && (ppo->get_ies_lam_mults().size() == 1))
    //{
    //	warnings.push_back("only one lambda mult to test, no point in using a subset");
    //	//ppo->set_ies_subset_size(100000000);
    //}

    string how = pest_scenario.get_pestpp_options().get_ies_subset_how();
    if ((how != "FIRST") && (how != "LAST") && (how != "RANDOM") && (how != "PHI_BASED"))
    {
        ss.str("");
        ss << "'subset_how' is '" << how << "' but should be 'FIRST','LAST','RANDOM','PHI_BASED'";
        errors.push_back(ss.str());
    }

    if ((ppo->get_ies_verbose_level() < 0) || (ppo->get_ies_verbose_level() > 3))
    {
        warnings.push_back("ies_verbose_level must be between 0 and 3, resetting to 3");
        ppo->set_ies_verbose_level(3);
    }
    if ((ppo->get_ies_no_noise()) && (ppo->get_obscov_filename().size() > 0))
    {
        ss.str("");
        ss << "ies_no_noise is true but obscov file supplied - these two are not compatible";
        errors.push_back(ss.str());
    }

    if ((ppo->get_obscov_filename().size() > 0) && (ppo->get_ies_drop_conflicts()))
    {
        ss.str("");
        ss << "use of a full obscov with ies_drop_conflicts is not currently supported";
        errors.push_back(ss.str());
    }
    if ((ppo->get_ies_obs_csv().size() > 0) && (ppo->get_ies_no_noise()))
    {
        ss.str("");
        ss << "ies_no_noise can't be used with an ies_observation_ensemble";
        errors.push_back(ss.str());
    }

    if ((ppo->get_ies_save_rescov()) && (pest_scenario.get_ctl_ordered_nz_obs_names().size() > 60000))
    {
        errors.push_back("'ies_save_rescov' requires too much memory for greater than 60,000 observations");
    }
    else if ((ppo->get_ies_save_rescov()) && (pest_scenario.get_ctl_ordered_nz_obs_names().size() > 30000))
    {
        warnings.push_back("'ies_save_rescov' with more than 30,000 observations will likely produce allocation errors, be prepared!");
    }
    if (ppo->get_ies_weights_csv().size() > 0)
    {
        warnings.push_back("ies_weight_ensemble is highly experimental - user beware");
        if (ppo->get_ies_phi_fractions_file().size() > 0)
        {
            warnings.push_back("weight ensemble used with internal weight adjustment is even more experimental!");
        }
    }
    if (ppo->get_ies_multimodal_alpha() > 1.0)
    {
        errors.push_back("multimodal alpha > 1.0");
    }
    if (ppo->get_ies_multimodal_alpha() < 0.001)
    {
        errors.push_back("multimodal alpha < 0.001");
    }

    for (auto& fac : ppo->get_ies_reinflate_factor())
    {
        if (fac <= 0.0)
        {
            warnings.push_back("reinflation factor <= 0.0");
        }
        else if (fac > 1.0)
        {
            warnings.push_back("reinflation factor > 1.0..user beware!");
        }
    }


    if (warnings.size() > 0)
    {
        message(0, "sanity_check warnings");
        for (auto &w : warnings)
            message(1, w);
        message(1, "continuing initialization...");
    }
    if (errors.size() > 0)
    {
        message(0, "sanity_check errors - uh oh");
        for (auto &e : errors)
            message(1, e);
        throw_em_error(string("sanity_check() found some problems - please review rec file"));
    }

}


void EnsembleMethod::throw_em_error(string message)
{
	performance_log->log_event(alg_tag + " error: " + message);
		cout << endl << "   ************   " << endl << "    " << alg_tag << " error: " << message << endl << endl;
		file_manager.rec_ofstream() << endl << "   ************   " << endl << "    " << alg_tag << " error: " << message << endl << endl;
		file_manager.close_file("rec");
		performance_log->~PerformanceLog();
		throw runtime_error(alg_tag + " error: " + message);
	
}

bool EnsembleMethod::should_terminate(int current_n_iter_mean)
{
	//todo: use ies accept fac here?
	double phiredstp = pest_scenario.get_control_info().phiredstp;
	int nphistp = pest_scenario.get_control_info().nphistp;
	int nphinored = pest_scenario.get_control_info().nphinored;
	bool phiredstp_sat = false, nphinored_sat = false, consec_sat = false;
	double phi, ratio;
	int count = 0;
	int nphired = 0;
	//best_mean_phis = vector<double>{ 1.0,0.8,0.81,0.755,1.1,0.75,0.75,1.2 };



	/*if ((!consec_sat )&& (best_mean_phis.size() == 0))
		return false;*/
	message(0, "phi-based termination criteria check");
	message(1, "phiredstp: ", phiredstp);
	message(1, "nphistp: ", nphistp);
	message(1, "nphinored (also used for consecutive bad lambda cycles): ", nphinored);
	//int n_mean_iter = pest_scenario.get_pestpp_options().get_ies_n_iter_reinflate();
    vector<double>::iterator begin_idx = best_mean_phis.begin();
    //if ((current_n_iter_mean > 0) && (best_mean_phis.size() > current_n_iter_mean))
    //    begin_idx = best_mean_phis.end() - (current_n_iter_mean+1); //bc of prior phi and then adding the mean shift to the list
    if (best_mean_phis.size() > 0)
	{

		vector<double>::iterator idx = min_element(begin_idx, best_mean_phis.end());
        nphired = (best_mean_phis.end() - idx) - 1;
		best_phi_yet = best_mean_phis[idx - best_mean_phis.begin()];// *pest_scenario.get_pestpp_options().get_ies_accept_phi_fac();
		message(1, "best mean phi sequence: ", best_mean_phis);
		message(1, "best phi yet: ", best_phi_yet);
	}
	message(1, "number of consecutive bad lambda testing cycles: ", consec_bad_lambda_cycles);
	if (consec_bad_lambda_cycles >= nphinored)
	{
		message(1, "number of consecutive bad lambda testing cycles > nphinored");
		consec_sat = true;
	}

    int i = 0;
    for (auto& phi : best_mean_phis)
	{
		ratio = (phi - best_phi_yet) / phi;
    	//if ((i>=(iter - current_n_iter_mean)) && (ratio <= phiredstp))
        if (ratio <= phiredstp)
			count++;
        i++;
	}
	message(1, "number of iterations satisfying phiredstp criteria: ", count);
	if (count >= nphistp)
	{
		message(1, "number iterations satisfying phiredstp criteria > nphistp");
		phiredstp_sat = true;
	}

	message(1, "number of iterations since best yet mean phi: ", nphired);
	if (nphired >= nphinored)
	{
		message(1, "number of iterations since best yet mean phi > nphinored");
		nphinored_sat = true;
	}
	if (best_mean_phis[best_mean_phis.size() - 1] == 0.0)
	{
		message(1, "mean phi is zero, all done");
		return true;
	}

	if ((nphinored_sat) || (phiredstp_sat) || (consec_sat))
	{
		message(1, "phi-based termination criteria satisfied, all done");
		return true;
	}
	int q = pest_utils::quit_file_found();
	if ((q == 1) || (q == 2))
    {
	    message(1,"'pest.stp' found, quitting");
	    return true;
    }
	else if (q == 4)
    {
	    message(0,"pest.stp found with '4'.  run mgr has returned control, removing file.");
	    if (!pest_utils::try_remove_quit_file())
        {
	        message(0,"error removing pest.stp file, bad times ahead...");
        }
    }
	return false;
}


vector<ObservationEnsemble> EnsembleMethod::run_lambda_ensembles(vector<ParameterEnsemble>& pe_lams, vector<double>& lam_vals, 
	vector<double>& scale_vals, int cycle, vector<int>& pe_subset_idxs, vector<int>& oe_subset_idxs)
{
	ofstream& frec = file_manager.rec_ofstream();
	stringstream ss;
	ss << "queuing " << pe_lams.size() << " ensembles";
	performance_log->log_event(ss.str());
	run_mgr_ptr->reinitialize();
	vector<string> names = pe_lams[0].get_real_names();
	ss.str("");
	for (auto i : pe_subset_idxs)
		ss << i << ":" << names[i] << ", ";
	message(1, "subset idx:pe real name: ", ss.str());
	ss.str("");
	names = oe.get_real_names();
	for (auto i : oe_subset_idxs)
		ss << i << ":" << names[i] << ", ";
	message(1, "subset idx:oe real name: ", ss.str());

	//set_subset_idx(pe_lams[0].shape().first);
	vector<map<int, int>> real_run_ids_vec;

	string additional_tag;
	for (int i=0;i<pe_lams.size();i++)
	{
	    ss.str("");
	    ss << " lambda:" << lam_vals[i] << " scale_fac:" << scale_vals[i] << " iteration:" << iter;
	    additional_tag = ss.str();
		try
		{
			real_run_ids_vec.push_back(pe_lams[i].add_runs(run_mgr_ptr, pe_subset_idxs,cycle,additional_tag));
		}
		catch (const exception& e)
		{
			stringstream ss;
			ss << "run_ensemble() error queueing runs: " << e.what();
			throw_em_error(ss.str());
		}
		catch (...)
		{
			throw_em_error(string("run_ensembles() error queueing runs"));
		}
	}
	performance_log->log_event("making runs");
	try
	{

		run_mgr_ptr->run();
	}
	catch (const exception& e)
	{
		stringstream ss;
		ss << "error running ensembles: " << e.what();
		throw_em_error(ss.str());
	}
	catch (...)
	{
		throw_em_error(string("error running ensembles"));
	}

	performance_log->log_event("processing runs");
	vector<int> failed_real_indices;
	vector<ObservationEnsemble> obs_lams;

	map<int, int> real_run_ids;
	for (int i = 0; i < pe_lams.size(); i++)
	{
		ObservationEnsemble _oe = oe;//copy
		vector<double> rep_vals{ lam_vals[i],scale_vals[i] };
		real_run_ids = real_run_ids_vec[i];

		if ((use_subset) && ((_oe.shape().first > pe_lams[i].shape().first) || (pe_subset_idxs.size() < pe_lams[i].shape().first)))
		{
			_oe.keep_rows(oe_subset_idxs);
			int ireal = 0;
			map<int, int> temp;
			for (auto& rri : real_run_ids)
			{
				temp[ireal] = rri.second;
				ireal++;
			}

			real_run_ids = temp;
		}

		try
		{
			failed_real_indices = _oe.update_from_runs(real_run_ids, run_mgr_ptr);
		}
		catch (const exception& e)
		{
			stringstream ss;
			ss << "error processing runs for lambda,scale: " << lam_vals[i] << ',' << scale_vals[i] << ':' << e.what();
			throw_em_error(ss.str());
		}
		catch (...)
		{
			stringstream ss;
			ss << "error processing runs for lambda,scale: " << lam_vals[i] << ',' << scale_vals[i];
			throw_em_error(ss.str());
		}

		if (pest_scenario.get_pestpp_options().get_ies_debug_fail_subset())
		{
			stringstream ss;
			if ((pe_lams.size() > 1) && (i == pe_lams.size() - 1))
			{
				ss << "'ies_debug_fail_subset' is true, failing all realizations for inflation factor " << lam_vals[i] << ", backtrack factor " << scale_vals[i];
				message(0, ss.str());
				for (int j=0;j<real_run_ids.size();j++)
					failed_real_indices.push_back(j);
			}
			else
			{
				ss << "'ies_debug_fail_subset' is true, failing last realization for inflation factor " << lam_vals[i] << ", backtrack factor " << scale_vals[i];
				message(0, ss.str());
				failed_real_indices.push_back(real_run_ids.size() - 1);
			}
		}
				
			

		if (failed_real_indices.size() > 0)
		{
			stringstream ss;
			vector<string> par_real_names = pe_lams[i].get_real_names();
			vector<string> obs_real_names = oe.get_real_names();
			vector<string> failed_par_names, failed_obs_names;
			string oname, pname;
			ss << "the following par:obs realization runs failed for lambda,scale " << lam_vals[i] << ',' << scale_vals[i] << "-->";
			for (auto& i : failed_real_indices)
			{
				pname = par_real_names[pe_subset_idxs[i]];
				oname = obs_real_names[oe_subset_idxs[i]];
				failed_par_names.push_back(pname);
				failed_obs_names.push_back(oname);
				ss << pname << ":" << oname << ',';
			}
			string s = ss.str();
			message(1, s);
			if (failed_real_indices.size() == _oe.shape().first)
			{
				message(0, "WARNING: all realizations failed for lambda, scale :", rep_vals);
				_oe = ObservationEnsemble(&pest_scenario);

			}
			else
			{
				performance_log->log_event("dropping failed realizations");

				_oe.drop_rows(failed_obs_names);
				pe_lams[i].drop_rows(failed_par_names);
			}

		}
		obs_lams.push_back(_oe);
	}
	return obs_lams;
}

pair<string,string> EnsembleMethod::save_ensembles(string tag, int cycle, ParameterEnsemble& _pe, ObservationEnsemble& _oe)
{
	stringstream ss;
	ss << file_manager.get_base_filename();
	if (cycle != NetPackage::NULL_DA_CYCLE)
		ss << "." << cycle;
	if (tag.size() > 0)
		ss << "." << tag;
	if (iter == -1)
		ss << "." << "prior" << ".obs";
	else if (iter == -2)
		ss << "." << "mean" << ".obs";
	else
		ss << "." << iter << ".obs";
    if (pest_scenario.get_pestpp_options().get_save_dense())
    {
        ss << dense_file_ext;
        _oe.to_dense(ss.str());
    }
	else if (pest_scenario.get_pestpp_options().get_save_binary())
	{
		ss << ".jcb";
		_oe.to_binary(ss.str());
	}
	else
	{
		ss << ".csv";
		_oe.to_csv(ss.str());
	}
	string oname = ss.str();
	ss.str("");
	ss << file_manager.get_base_filename();
	if (cycle != NetPackage::NULL_DA_CYCLE)
		ss << "." << cycle;
	if (tag.size() > 0)
		ss << "." << tag;
	if (iter == -1)
		ss << "." << "prior" << ".par";
	else if (iter == -2)
		ss << "." << "mean" << ".par";
	else
		ss << "." << iter << ".par";
    if (pest_scenario.get_pestpp_options().get_save_dense())
    {
        ss << dense_file_ext;
        _pe.to_dense(ss.str());
    }
    else if (pest_scenario.get_pestpp_options().get_save_binary())
	{
		ss << ".jcb";
		_pe.to_binary(ss.str());
	}
	else
	{
		ss << ".csv";
		_pe.to_csv(ss.str());
	}
	string pname = ss.str();
	return pair<string, string>(pname, oname);
}

void EnsembleMethod::report_and_save(int cycle)
{
	ofstream& frec = file_manager.rec_ofstream();
	frec << endl << "  ---  " << alg_tag << " iteration " << iter << " report  ---  " << endl;
	frec << "   number of active realizations:  " << pe.shape().first << endl;
	frec << "   number of model runs:           " << run_mgr_ptr->get_total_runs() << endl;

	cout << endl << "  ---  " << alg_tag << " iteration " << iter << " report  ---  " << endl;
	cout << "   number of active realizations:   " << pe.shape().first << endl;
	cout << "   number of model runs:            " << run_mgr_ptr->get_total_runs() << endl;

	pair<string, string> names = save_ensembles(string(), cycle, pe, oe);
	frec << "      current obs ensemble saved to " << names.second << endl;
	cout << "      current obs ensemble saved to " << names.second << endl;
	frec << "      current par ensemble saved to " << names.first << endl;
	cout << "      current par ensemble saved to " << names.first << endl;

    map<string,double> base_weights = weights.get_real_map(BASE_REAL_NAME, true);
	save_real_par_rei(pest_scenario, pe, oe, output_file_writer, file_manager, iter, BASE_REAL_NAME, cycle,base_weights);
	save_real_par_rei(pest_scenario, pe, oe, output_file_writer, file_manager, -1, BASE_REAL_NAME, cycle,base_weights);
}


void EnsembleMethod::save_mat(string prefix, Eigen::MatrixXd& mat)
{
	stringstream ss;
	ss << iter << '.' << prefix;
	try
	{
		ofstream& f = file_manager.open_ofile_ext(ss.str());
		f << mat << endl;
		f.close();
		file_manager.close_file(ss.str());
	}
	catch (...)
	{
		message(1, "error saving matrix", ss.str());
	}
}

vector<int> EnsembleMethod::run_ensemble(ParameterEnsemble& _pe,
	ObservationEnsemble& _oe, const vector<int>& real_idxs, int cycle)
{
	stringstream ss;
	ss.str("");
	ss << " iteration:" << iter;
	vector<int> failed_real_indices;
	try
	{
		failed_real_indices = run_ensemble_util(performance_log, file_manager.rec_ofstream(), 
			_pe, _oe, run_mgr_ptr, 
			pest_scenario.get_pestpp_options().get_debug_check_par_en_consistency(), 
			real_idxs, cycle,ss.str());
	}
	catch (const exception& e)
	{
		ss.str("");
		ss << "run_ensemble() error: " << e.what();
		throw_em_error(ss.str());
	}
	return failed_real_indices;
}


void EnsembleMethod::initialize(int cycle, bool run, bool use_existing)
{
	message(0, "initializing");
	pp_args = pest_scenario.get_pestpp_options().get_passed_args();

	act_obs_names = pest_scenario.get_ctl_ordered_nz_obs_names();
	act_par_names = pest_scenario.get_ctl_ordered_adj_par_names();

	stringstream ss;

    map<string,vector<string>> group_map,group_to_obs_map;
    map<string,map<string,double>> phi_fracs_by_real;
    vector<string> index;
    check_and_fill_phi_factors(group_to_obs_map,group_map,phi_fracs_by_real,index,false);

	if (pest_scenario.get_control_info().noptmax == 0)
	{
        if (pest_scenario.get_pestpp_options().get_debug_parse_only()) {
            return;
        }

        ParamTransformSeq pts = pe.get_par_transform();
		
		ParameterEnsemble _pe(&pest_scenario, &rand_gen);
		if (cycle == NetPackage::NULL_DA_CYCLE)
		{
			message(0, "'noptmax'=0, running control file parameter values and quitting");
			_pe.reserve(vector<string>(), pest_scenario.get_ctl_ordered_adj_par_names());
			_pe.set_trans_status(ParameterEnsemble::transStatus::NUM);
			Parameters pars = pest_scenario.get_ctl_parameters();
			pts.ctl2numeric_ip(pars);
			_pe.append(BASE_REAL_NAME, pars);
		}
		else
		{
			map<string, int> rmap = pe.get_real_map();
			if (rmap.find(BASE_REAL_NAME) == rmap.end())
			{
				message(0, "'noptmax'=0, running current cycle mean parameter values");
				message(1, "calculating mean parameter values");
				Parameters pars = pest_scenario.get_ctl_parameters();
				pars.update(pe.get_var_names(), pe.get_mean_stl_var_vector());
				_pe.reserve(vector<string>(), pe.get_var_names());
				_pe.set_trans_status(pe.get_trans_status());
				_pe.append("mean", pars);
			}
			else
			{
				message(0, "'noptmax'=0, running current cycle base parameter values");
				_pe.reserve(vector<string>(), pest_scenario.get_ctl_ordered_adj_par_names());
				_pe.set_trans_status(ParameterEnsemble::transStatus::NUM);
				Parameters pars = pest_scenario.get_ctl_parameters();
				pars.update(pe.get_var_names(), pe.get_real_vector(BASE_REAL_NAME));
				_pe.set_trans_status(pe.get_trans_status());
				_pe.append(BASE_REAL_NAME, pars);
			}
		}
		ss.str("");
		ss << file_manager.get_base_filename();
		if (cycle != NetPackage::NULL_DA_CYCLE)
			ss << "." << cycle;
		ss << "." << pest_utils::lower_cp(BASE_REAL_NAME) << ".par.csv";
		string par_csv = ss.str();
		message(1, "saving control file parameter ensemble to ", par_csv);
		//_pe.to_csv(par_csv);
		pe_base = _pe;
		pe_base.reorder(vector<string>(), act_par_names);
		ObservationEnsemble _oe(&pest_scenario, &rand_gen);
		_oe.reserve(vector<string>(), pest_scenario.get_ctl_ordered_obs_names());
		_oe.append(BASE_REAL_NAME, pest_scenario.get_ctl_observations());
		oe_base = _oe;
		oe_base.reorder(vector<string>(), act_obs_names);
		initialize_parcov();
		//initialize the phi handler
		ph = L2PhiHandler(&pest_scenario, &file_manager, &oe_base, &pe_base, &parcov);
		if (ph.get_lt_obs_names().size() > 0)
		{
			message(1, "less_than inequality defined for observations: ", ph.get_lt_obs_names().size());
		}
		if (ph.get_gt_obs_names().size())
		{
			message(1, "greater_than inequality defined for observations: ", ph.get_gt_obs_names().size());
		}
        map<string,double> t;
        t = ph.get_lt_obs_bounds();
        if (!t.empty())
        {
            ss.str("");
            ss << "less_than inequality defined through 'less_than' data for observations:" << endl;
            for (const auto it : t)
            {
                ss << it.first << "," << it.second << endl;
            }
            ss << endl;
            message(1,ss.str());
        }
        t = ph.get_gt_obs_bounds();
        if (!t.empty())
        {
            ss.str("");
            ss << "greater_than inequality defined through 'greater_than' data for observations:" << endl;
            for (const auto it : t)
            {
                ss << it.first << "," << it.second << endl;
            }
            ss << endl;
            message(1,ss.str());
        }
        t.clear();
        map<string,pair<double,double>> tt = ph.get_double_obs_bounds();
        if (!tt.empty())
        {
            ss.str("");
            ss << "double inequality defined through 'greater_than' and 'less_than' data for " << tt.size() << " observations" << endl;
            ss.str("");
            ss << "double inequality defined through 'greater_than' and 'less_than' data for observations:" << endl;
            for (const auto it : tt)
            {
                ss << it.first << "," << it.second.first << " to " << it.second.second << endl;
            }
            ss << endl;
            message(1,ss.str(),false);

        }


		message(1, "running control file parameter values");

		vector<int> failed_idxs = run_ensemble(_pe, _oe,vector<int>(),cycle);
		if (failed_idxs.size() != 0)
		{
			message(0, "control file parameter value run failed...bummer");
			throw_em_error("control file parameter value run failed");
		}

		ss.str("");
		ss << file_manager.get_base_filename();
		if (cycle != NetPackage::NULL_DA_CYCLE)
			ss << "." << cycle;
		ss << "." << pest_utils::lower_cp(BASE_REAL_NAME) << ".obs.csv";
		string obs_csv = ss.str();
		message(1, "saving results from control file parameter value run to ", obs_csv);
		_oe.to_csv(obs_csv);

		ph.update(_oe, _pe);
		message(0, "control file parameter phi report:");
		ph.report(true);
		ph.write(0, 1);
		save_real_par_rei(pest_scenario, _pe, _oe, output_file_writer, file_manager, -1, BASE_REAL_NAME, cycle);
		//transfer_dynamic_state_from_oe_to_initial_pe(_pe, _oe);
		pe = _pe;
		oe = _oe;

		return;
	}

	//set some defaults
	PestppOptions* ppo = pest_scenario.get_pestpp_options_ptr();

	if (ppo->get_ies_use_mda())
	{
	    message(1, "using multiple-data-assimilation algorithm");
		int noptmax = pest_scenario.get_control_info().noptmax;
		if (noptmax > 0)
		{

            if (ppo->get_ies_no_noise())
            {
                throw_em_error("'no noise'is not compatible with 'use_mda' as this solution relies on noise draws");
            }
		}
	}
	else
	{
		message(1, "using glm algorithm");
	}

	verbose_level = pest_scenario.get_pestpp_options_ptr()->get_ies_verbose_level();

	message(1, "using REDSVD for truncated svd solve");
	message(1, "maxsing:", pest_scenario.get_svd_info().maxsing);
	message(1, "eigthresh: ", pest_scenario.get_svd_info().eigthresh);
	if (localizer.is_initialized())
	{
		message(2, "using previously initialized localizer");
		use_localizer = localizer.get_use();
	}
	else
	{
		message(1, "initializing localizer");
		bool forgive_missing = pest_scenario.get_pestpp_options().get_ies_localizer_forgive_missing();
		ofstream& frec = file_manager.rec_ofstream();
		use_localizer = localizer.initialize(performance_log, frec, forgive_missing);
	}
	num_threads = pest_scenario.get_pestpp_options().get_ies_num_threads();
	if (!use_localizer)
		message(1, "not using localization");
	else
	{
		if (localizer.get_autoadaloc())
		{
			message(1, "using automatic adaptive localization");
			message(2, "with autoadaloc_sigma_dist ", ppo->get_ies_autoadaloc_sigma_dist());
		}
		if (localizer.get_filename().size() > 0)
		{
			message(1, "using localization matrix " + localizer.get_filename());
			localizer.report(file_manager.rec_ofstream());
		}
	}
	if ((use_localizer) && (!localizer.get_autoadaloc()))
	{
		ss.str("");
		ss << "using localized solution with " << localizer.get_num_upgrade_steps() << " sequential upgrade steps";
		message(1, ss.str());
		ss.str("");
		if (num_threads > 0)
		{
			ss.str("");
			ss << "using multithreaded localization calculation with " << num_threads << " threads";
			message(1, ss.str());

		}
		if (localizer.get_how() == Localizer::How::OBSERVATIONS)
			message(1, "localizing by observations");
		else
			message(1, "localizing by parameters");
	}
	iter = 0;
	//ofstream &frec = file_manager.rec_ofstream();
	last_best_mean = 1.0E+30;
	last_best_std = 1.0e+30;
	lambda_max = 1.0E+30;
	lambda_min = 1.0E-30;

	consec_bad_lambda_cycles = 0;
    reinflate_to_minphi_real = false;
    bool use_min = false;
    for (auto& fac : pest_scenario.get_pestpp_options().get_ies_n_iter_reinflate())
    {
        if (fac < 0)
        {
            use_min = true;
            //fac *= -1;
        }
    }

    if (use_min)
    {
        message(2,"n_iter_reinflate < 0, using min-phi real for re-inflation, resetting n_iter_reinflate to positive");
        reinflate_to_minphi_real = true;
        //pest_scenario.get_pestpp_options_ptr()->set_ies_n_iter_reinflate(-1 * pest_scenario.get_pestpp_options().get_ies_n_iter_reinflate());
    }
	lam_mults = pest_scenario.get_pestpp_options().get_ies_lam_mults();
	if (lam_mults.size() == 0)
		lam_mults.push_back(1.0);
	message(1, "using lambda multipliers: ", lam_mults);
	vector<double> scale_facs = pest_scenario.get_pestpp_options().get_lambda_scale_vec();
	message(1, "using lambda scaling factors: ", scale_facs);
	double acc_fac = pest_scenario.get_pestpp_options().get_ies_accept_phi_fac();
	message(1, "acceptable phi factor: ", acc_fac);
	double inc_fac = pest_scenario.get_pestpp_options().get_ies_lambda_inc_fac();
	message(1, "lambda increase factor: ", inc_fac);
	double dec_fac = pest_scenario.get_pestpp_options().get_ies_lambda_dec_fac();
	message(1, "lambda decrease factor: ", dec_fac);
	message(1, "max run fail: ", ppo->get_max_run_fail());

    prep_drop_violations();

	//todo: add sanity checks for weights en
	//like conflict with da weight cycle table
	sanity_checks();

	bool echo = false;
	if (verbose_level > 1)
		echo = true;

	initialize_parcov();
	initialize_obscov();

	int subset_size = pest_scenario.get_pestpp_options().get_ies_subset_size();

    reg_factor = pest_scenario.get_pestpp_options().get_ies_reg_factor();
    if (reg_factor < 0)
    {
        reg_factor *= -1.0;
        message(1, "using reg_factor in upgrade calculations with full solution: ", reg_factor);
        //reset the passed reg factor to 0.0 so that the phi handler wont try to use it
        pest_scenario.get_pestpp_options_ptr()->set_ies_reg_factor(0.0);

        if (!pest_scenario.get_pestpp_options().get_ies_use_approx())
        {
            message(1, "WARNING: negative reg_factor passed, implying a full solution, resetting 'ies_use_approx' to false ");
            pest_scenario.get_pestpp_options_ptr()->set_ies_use_approx(false);
        }
    }

	double bad_phi = pest_scenario.get_pestpp_options().get_ies_bad_phi();
	if (bad_phi < std::numeric_limits<double>::max())
		message(1, "using bad_phi: ", bad_phi);

	int num_reals = pest_scenario.get_pestpp_options().get_ies_num_reals();

	if ((pe.shape().first > 0) && (cycle != NetPackage::NULL_DA_CYCLE))
	{
		message(2, "using previously initialized parameter ensemble");
		pe_drawn = false;
	}
	else
		pe_drawn = initialize_pe(parcov);
	
	if (pest_scenario.get_pestpp_options().get_ies_use_prior_scaling())
	{
		message(1, "forming inverse sqrt of prior parameter covariance matrix");

		if (parcov.isdiagonal())
			parcov_inv_sqrt = parcov.inv(echo).get_matrix().diagonal().cwiseSqrt().asDiagonal();
		else
		{
			message(1, "first extracting diagonal from prior parameter covariance matrix");
			Covariance parcov_diag;
			parcov_diag.from_diagonal(parcov);
			parcov_inv_sqrt = parcov_diag.inv(echo).get_matrix().diagonal().cwiseSqrt().asDiagonal();
		}
	}
	else {
		message(1, "not using prior parameter covariance matrix scaling");
	}

    //check to see if any explicit obs noise options are set
    bool reset_to_nonoise = true;
	set<string> passed = ppo->get_passed_args();
    if (passed.find("IES_NO_NOISE") != passed.end())
    {
        reset_to_nonoise = false;
    }
    else if (!ppo->get_ies_obs_csv().empty())
        reset_to_nonoise = false;
    else if (!ppo->get_obscov_filename().empty())
        reset_to_nonoise = false;
    else if (ppo->get_ies_use_mda())
        reset_to_nonoise = false;
    else
    {
        map<string, double> obs_std = pest_scenario.get_ext_file_double_map("observation data external", "standard_deviation");
        if (obs_std.size() > 0)
            reset_to_nonoise = false;
    }

    if (reset_to_nonoise)
    {
        ss.str("");
        ss << "NOTE: no obs-noise-specific options have been passed, resetting to `ies_no_noise` to true";
        message(0,ss.str());
        ppo->set_ies_no_noise(true);
    }


	if ((oe_base.shape().first > 0) && (cycle != NetPackage::NULL_DA_CYCLE)) {
        message(2, "using previously initialized observation (simulated output) ensemble");
        oe_drawn = false;
    }
	else
		oe_drawn = initialize_oe(obscov);

	string center_on = ppo->get_ies_center_on();

	try
	{
		pe.check_for_dups();
	}
	catch (const exception& e)
	{
		string message = e.what();
		throw_em_error("error in parameter ensemble: " + message);
	}
	try
	{
		oe_base.check_for_dups();
	}
	catch (const exception& e)
	{
		string message = e.what();
		throw_em_error("error in observation ensemble: " + message);
	}

	if (pe.shape().first != oe_base.shape().first)
	{
		//the special case where par en < obs en and all par reals are found in obs en...

		if (pe.shape().first < oe_base.shape().first)
		{
			vector<string> oe_names = oe_base.get_real_names();
			set<string> oset(oe_names.begin(), oe_names.end());
			vector<string> missing;
			for (auto n : pe.get_real_names())
				if (oset.find(n) == oset.end())
				{
					missing.push_back(n);
				}
			if (missing.size() == 0)
			{
				ss.str("");
				ss << "par en has " << pe.shape().first << " realizations, compared to " << oe_base.shape().first << " obs+noise realizations";
				message(1, ss.str());
				message(1, " the realization names are compatible");
				message(1, "re-indexing obs+noise en to align with par en...");
                if (verbose_level > 1) {
                    cout << "oe names: " << endl;
                    for (auto &name : oe_names)
                        cout << name << endl;
                    oe_names = pe.get_real_names();
                    cout << endl << "pe names: " << endl;
                    for (auto &name : oe_names)
                        cout << name << endl;
                }

                    oe_base.reorder(pe.get_real_names(), vector<string>());

                if (verbose_level > 1)
                {
                    oe_names = oe_base.get_real_names();
                    cout << "new oe names: " << endl;
                    for (auto &name : oe_names)
                        cout << name << endl;
                    cout << endl;
                }
			}
			else
			{
				ss.str("");
				ss << "the following par en real names were not found in the obs en: ";
				for (auto m : missing)
				{
					ss << m << ",";
				}
				throw_em_error(ss.str());
			}
		}
		else
		{
			ss.str("");
			ss << "parameter ensemble rows (" << pe.shape().first << ") not equal to observation ensemble rows (";
			ss << oe_base.shape().first << ")";
			throw_em_error(ss.str());
		}
	}

	string obs_restart_csv = pest_scenario.get_pestpp_options().get_ies_obs_restart_csv();
	string par_restart_csv = pest_scenario.get_pestpp_options().get_ies_par_restart_csv();


    initialize_weights();
    ss.str("");
    if (pest_scenario.get_pestpp_options().get_save_dense())
    {
        ss << file_manager.get_base_filename();
        if (cycle != NetPackage::NULL_DA_CYCLE)
            ss << "." << cycle;
        ss << ".weights" << dense_file_ext;
        weights.to_dense(ss.str());
    }
    else if (pest_scenario.get_pestpp_options().get_save_binary())
    {
        ss << file_manager.get_base_filename();
        if (cycle != NetPackage::NULL_DA_CYCLE)
            ss << "." << cycle;
        ss << ".weights.jcb";
        weights.to_binary(ss.str());
    }
    else
    {
        ss << file_manager.get_base_filename();
        if (cycle != NetPackage::NULL_DA_CYCLE)
            ss << "." << cycle;
        ss << ".weights.csv";
        weights.to_csv(ss.str());
    }
    message(1, "saved weight ensemble to ", ss.str());
    message(2, "checking for denormal values in weights ensemble");
    weights.check_for_normal("weights ensemble");
    // if no restart and a pe was passed and no oe was passed, reset here before adding base
	if (pest_scenario.get_pestpp_options().get_ies_include_base()) {
        if (pp_args.find("IES_RESTART_OBS_EN") != pp_args.end()) {
            message(1,
                    "Warning: even though `ies_include_base` is true, you passed a restart obs en, not adding 'base' realization...");
        } else {
            add_bases();
        }
    }
	if ((obs_restart_csv.size() == 0) && (!pe_drawn) && (oe_drawn))
	{
		vector<string> rnames = pe.get_real_names();
		oe_base.set_real_names(rnames,true);
		weights.set_real_names(rnames,true);
		message(2, "resetting obs + noise ensemble real names to parameter ensemble real names");
	}

	//need this here for Am calcs...
	//message(1, "transforming parameter ensemble to numeric");
	pe.transform_ip(ParameterEnsemble::transStatus::NUM);

	//now we check to see if we need to try to align the par and obs en
	//this would only be needed if either of these were not drawn
	if (!pe_drawn || !oe_drawn)
	{
		bool aligned = pe.try_align_other_rows(performance_log, oe_base);
		if (aligned)
		{
			message(2, "observation ensemble reordered to align rows with parameter ensemble");
		}
	}

	//just check to see if common real names are found but are not in the same location
	map<string, int> pe_map = pe.get_real_map(), oe_map = oe_base.get_real_map();
	vector<string> misaligned;
	for (auto item : pe_map)
	{
		if (oe_map.find(item.first) == oe_map.end())
			continue;
		if (item.second != oe_map[item.first])
			misaligned.push_back(item.first);
	}
	if (misaligned.size() > 0)
	{
		message(1, "WARNING: common realization names shared between the parameter and observation + noise ensembles but they are not in the same row locations, see .rec file for listing");
		ofstream& frec = file_manager.rec_ofstream();
		frec << endl << "WARNING: the following " << misaligned.size() << " realization names are shared between the ";
		frec << "parameter and observation ensembles but they are not in the same row locations:" << endl;
		for (auto ma : misaligned)
			frec << ma << endl;
	}

    int num_par_elements = pest_scenario.get_ctl_parameters().size();
    int num_obs_elements = pest_scenario.get_ctl_observations().size();
    if ((num_par_elements > 1000000) || (num_obs_elements > 1000000))
    {
        if ((pest_scenario.get_pestpp_options().get_save_binary()) && (!pest_scenario.get_pestpp_options().get_save_dense()))
        {
            message(0,"WARNING: npar and/or nobs > 1e6, you are close to going out-of-range for jcb format.  Switching to dense format but using '.jcb' file extension");
            pest_scenario.get_pestpp_options_ptr()->set_save_dense(true);
            dense_file_ext = "jcb";
        }
    }

    if ((!pest_scenario.get_pestpp_options().get_save_dense()) && (pp_args.find("IES_ORDERED_BINARY") == pp_args.end()) && (pp_args.find("DA_ORDERED_BINARY") == pp_args.end()))
	{

		//if ((pe.shape().second > 100000) && (pest_scenario.get_pestpp_options().get_save_binary()))
        if (((num_par_elements > 100000) || (num_obs_elements > 100000)) && (pest_scenario.get_pestpp_options().get_save_binary()))
		{
			message(1, "'ies_ordered_binary' was not passed, but 'ies_save_binary' is true and npar > 100,000, switching to unordered binary...");
			pest_scenario.get_pestpp_options_ptr()->set_ies_ordered_binary(false);
		}
	}

	if ((pp_args.find("IES_UPGRADES_IN_MEMORY") == pp_args.end()) && (pp_args.find("DA_UPGRADES_IN_MEMORY") == pp_args.end()))
	{
		if (num_par_elements > 100000)
		{
			message(1, "'ies_upgrades_in_memory' was not passed, but npar > 100,000, switching ies_upgrades_in_memory to false...");
			pest_scenario.get_pestpp_options_ptr()->set_ies_upgrades_in_memory(false);
		}
	}


	message(2, "checking for denormal values in pe");
	pe.check_for_normal("initial transformed parameter ensemble");
	ss.str("");

    if (pest_scenario.get_pestpp_options().get_save_dense())
    {
        ss << file_manager.get_base_filename();
        if (cycle != NetPackage::NULL_DA_CYCLE)
            ss << "." << cycle;
        ss << ".0.par" << dense_file_ext;
        pe.to_dense_unordered(ss.str());
    }
	else if (pest_scenario.get_pestpp_options().get_save_binary())
	{
		ss << file_manager.get_base_filename();
		if (cycle != NetPackage::NULL_DA_CYCLE)
			ss << "." << cycle;
		ss << ".0.par.jcb";
		pe.to_binary(ss.str());
	}
	else
	{
		ss << file_manager.get_base_filename();
		if (cycle != NetPackage::NULL_DA_CYCLE)
			ss << "." << cycle;
		ss << ".0.par.csv";
		pe.to_csv(ss.str());
	}
	message(1, "saved initial parameter ensemble to ", ss.str());
	message(2, "checking for denormal values in obs + noise ensemble");
	oe_base.check_for_normal("obs+noise observation ensemble");
	ss.str("");
    if (pest_scenario.get_pestpp_options().get_save_dense())
    {
        ss << file_manager.get_base_filename();
        if (cycle != NetPackage::NULL_DA_CYCLE)
            ss << "." << cycle;
        ss << ".obs+noise" + dense_file_ext;
        oe_base.to_dense(ss.str());
    }
	else if (pest_scenario.get_pestpp_options().get_save_binary())
	{
		ss << file_manager.get_base_filename();
		if (cycle != NetPackage::NULL_DA_CYCLE)
			ss << "." << cycle;
		ss << ".obs+noise.jcb";
		oe_base.to_binary(ss.str());
	}
	else
	{
		ss << file_manager.get_base_filename();
		if (cycle != NetPackage::NULL_DA_CYCLE)
			ss << "." << cycle;
		ss << ".obs+noise.csv";
		oe_base.to_csv(ss.str());
	}
	message(1, "saved obs+noise observation ensemble (obsval + noise realizations) to ", ss.str());

    ss.str("");

	if (subset_size > pe.shape().first)
	{
		use_subset = false;
	}
	else
	{
	    if (subset_size < 0)
        {
	        message(1, "using subset in lambda testing, percentage of realizations used in subset testing: ", -1. * subset_size);
        }
	    else
        {
            message(1, "using subset in lambda testing, number of realizations used in subset testing: ", subset_size);
        }
		string how = pest_scenario.get_pestpp_options().get_ies_subset_how();
		message(1, "subset how: ", how);
		use_subset = true;
	}

	oe_org_real_names = oe.get_real_names();
	pe_org_real_names = pe.get_real_names();

	if ((oe.shape().first > 0) && (cycle != NetPackage::NULL_DA_CYCLE))
	{
		message(2, "using previously initialized observation + noise ensemble");
	}
	else
		oe = oe_base; //copy
	
	//reorder this for later...
	if (act_obs_names.size() > 0)
		oe_base.reorder(vector<string>(), act_obs_names);

	pe_base = pe; //copy
	//reorder this for later
	if (pe.shape().second == act_par_names.size()) {
		pe_base.reorder(vector<string>(), act_par_names);
	}

	//the hard way to restart
	if (obs_restart_csv.size() > 0)
		initialize_restart();
    //we need this for the prior mean shifting
    weights_base = weights;

    if (pest_scenario.get_control_info().noptmax == -2)
    {
        if (pest_scenario.get_pestpp_options().get_debug_parse_only()) {
            return;
        }
        string rname = ppo->get_ies_run_realname();
        string org_rname = rname;
        Parameters pars;
        if (!rname.empty())
        {
            map<string,int> rmap = pe.get_real_map();
            if (rmap.find(rname) == rmap.end())
            {
                throw_em_error("'ies_run_realname' argument supplied '"+rname+"' but that realization is not in the parameter ensemble");
            }
            pars.update(pe.get_var_names(),pe.get_real_vector(rname));
            rname = "real"+rname;

        }
        else {
            rname = "mean";
            org_rname = "mean";
            message(0, "'noptmax'=-2, running mean parameter ensemble values and quitting");
            message(1, "calculating mean parameter values");
            pars.update(pe.get_var_names(), pe.get_mean_stl_var_vector());


            if (pe.get_fixed_info().get_map_size() > 0) {
                ss.str("");
                ss << "WARNING: 'fixed' parameter realizations provided but ctrl " << endl;
                ss << "         file parameter values are being used for 'fixed' parameters" << endl;
                ss << "         in the mean parameter value run." << endl;
                message(0, ss.str());
            }
        }


        ParamTransformSeq pts = pe.get_par_transform();

        ParameterEnsemble _pe(&pest_scenario, &rand_gen);
        _pe.reserve(vector<string>(), pe.get_var_names());
        _pe.set_trans_status(pe.get_trans_status());

        if (rname != "mean") {
            _pe.set_fixed_info(pe.get_fixed_info());
        }
        _pe.append(org_rname, pars);
        ss.str("");
        ss << file_manager.get_base_filename();
        if (cycle != NetPackage::NULL_DA_CYCLE)
            ss << "." << cycle;
        ss << "." <<  rname << ".par.csv";
        string par_csv = ss.str();
        message(1, "saving "+ rname +" parameter values to ", par_csv);
        _pe.to_csv(par_csv);
        pe_base = _pe;
        pe_base.reorder(vector<string>(), act_par_names);
        ObservationEnsemble _oe(&pest_scenario, &rand_gen);
        _oe.reserve(vector<string>(), oe_base.get_var_names());
        _oe.reserve(vector<string>(), pest_scenario.get_ctl_ordered_obs_names());
        _oe.append(org_rname, pest_scenario.get_ctl_observations());
        oe_base = _oe;
        oe_base.reorder(vector<string>(), act_obs_names);
        //initialize the phi handler
        ph = L2PhiHandler(&pest_scenario, &file_manager, &oe_base, &pe_base, &parcov);
        if (ph.get_lt_obs_names().size() > 0)
        {
            message(1, "less_than inequality defined for observations: ", ph.get_lt_obs_names().size());
        }
        if (ph.get_gt_obs_names().size())
        {
            message(1, "greater_than inequality defined for observations: ", ph.get_gt_obs_names().size());
        }
        map<string,double> t;
        t = ph.get_lt_obs_bounds();
        if (!t.empty())
        {
            ss.str("");
            ss << "less_than inequality defined through 'less_than' data for observations:" << endl;
            for (const auto it : t)
            {
                ss << it.first << "," << it.second << endl;
            }
            ss << endl;
            message(1,ss.str());
        }
        t = ph.get_gt_obs_bounds();
        if (!t.empty())
        {
            ss.str("");
            ss << "greater_than inequality defined through 'greater_than' data for observations:" << endl;
            for (const auto it : t)
            {
                ss << it.first << "," << it.second << endl;
            }
            ss << endl;
            message(1,ss.str());
        }
        t.clear();
        map<string,pair<double,double>> tt = ph.get_double_obs_bounds();
        if (!tt.empty())
        {
            ss.str("");
            ss << "double inequality defined through 'greater_than' and 'less_than' data for observations:" << endl;
            for (const auto it : tt)
            {
                ss << it.first << "," << it.second.first << " to " << it.second.second << endl;
            }
            ss << endl;
            message(1,ss.str());
        }


        message(1, "running " + rname + " parameter values");

        vector<int> failed_idxs = run_ensemble(_pe, _oe,vector<int>(),cycle);
        if (failed_idxs.size() != 0)
        {
            message(0, rname+" parameter value run failed...bummer");
            throw_em_error(rname+" parameter value run failed");
        }
        ss.str("");
        ss << file_manager.get_base_filename();
        if (cycle != NetPackage::NULL_DA_CYCLE)
            ss << "." << cycle;
        ss << "." << rname << ".obs.csv";
        string obs_csv = ss.str();
        message(1, "saving results from mean parameter value run to ", obs_csv);
        _oe.to_csv(obs_csv);

        ph.update(_oe, _pe);
        message(0,"realization "+ rname+" phi report:");
        ph.report(true);
        ph.write(0, 1);
        save_real_par_rei(pest_scenario, _pe, _oe, output_file_writer, file_manager, -1, org_rname, cycle);
        //transfer_dynamic_state_from_oe_to_initial_pe(_pe, _oe);
        pe = _pe;
        oe = _oe;
        return;
    }

	if (!run)
		return;

	//check for center on 
	if (center_on.size() > 0)
	{
		ss.str("");
		if (center_on == MEDIAN_CENTER_ON_NAME)
		{
			ss << "centering on ensemble median value";
			message(1, ss.str());
			pe.get_eigen_anomalies(center_on);
		}
		else
		{
			ss << "centering on realization: '" << center_on << "' ";
			message(1, ss.str());
			vector<string> names = pe.get_real_names();
			if (find(names.begin(), names.end(), center_on) == names.end())
				throw_em_error("'ies_center_on' realization not found in par en: " + center_on);
			names = oe.get_real_names();
			if (find(names.begin(), names.end(), center_on) == names.end())
				throw_em_error("'ies_center_on' realization not found in obs en: " + center_on);
		}
	}
	else
		message(1, "centering on ensemble mean vector");
    if (pest_scenario.get_pestpp_options().get_debug_parse_only()) {
        return;
    }

	int reinflate_num_reals = pest_scenario.get_pestpp_options().get_ies_reinflate_num_reals()[0];
	if ((abs(reinflate_num_reals) > 0) && (abs(reinflate_num_reals) < pe_base.shape().first))
	{
		ss.str("");
		ss << "'ies_reinflate_num_reals[0] < current ensemble size, truncating ensemble to " << reinflate_num_reals << " realizations";
		message(1,ss.str());
		vector<string> tpar,tobs;
		vector<string> pebase_real_names = pe.get_real_names();
		vector<string> oebase_real_names = oe.get_real_names();
		bool has_base = false;
		for (auto& n : pe.get_real_names())
		{
			if (n == BASE_REAL_NAME)
			{
				has_base = true;
				break;
			}
		}
		bool found_base = false;
		for (int i=0;i<abs(reinflate_num_reals);i++)
		{
			tpar.push_back(pebase_real_names[i]);
			tobs.push_back(oebase_real_names[i]);
			if (pebase_real_names[i] == BASE_REAL_NAME)
			{
				found_base = true;
			}
		}
		pebase_real_names = tpar;
		oebase_real_names = tobs;
		if ((has_base) && (!found_base))
		{
			pebase_real_names[pebase_real_names.size() -1] = BASE_REAL_NAME;
			oebase_real_names[oebase_real_names.size() -1] = BASE_REAL_NAME;
		}
		pe.keep_rows(pebase_real_names);
		oe.keep_rows(oebase_real_names);

		ss.str("");

		if (pest_scenario.get_pestpp_options().get_save_dense())
		{
			ss << file_manager.get_base_filename();
			if (cycle != NetPackage::NULL_DA_CYCLE)
				ss << "." << cycle;
			ss << ".0.par" << dense_file_ext;
			pe.to_dense_unordered(ss.str());
		}
		else if (pest_scenario.get_pestpp_options().get_save_binary())
		{
			ss << file_manager.get_base_filename();
			if (cycle != NetPackage::NULL_DA_CYCLE)
				ss << "." << cycle;
			ss << ".0.par.jcb";
			pe.to_binary(ss.str());
		}
		else
		{
			ss << file_manager.get_base_filename();
			if (cycle != NetPackage::NULL_DA_CYCLE)
				ss << "." << cycle;
			ss << ".0.par.csv";
			pe.to_csv(ss.str());
		}
		message(1, "saved initial truncated parameter ensemble to ", ss.str());
	}



    //ok, now run the prior ensemble - after checking for center_on
	//in case something is wrong with center_on
	if ((obs_restart_csv.size() == 0) && (!use_existing)) {
        performance_log->log_event("running initial ensemble");
        message(1, "running initial ensemble of size", oe.shape().first);
        //get these here so we can use them later for reporting since pe and oe have fails dropped...
        vector<string> pnames = pe.get_real_names(),onames = oe.get_real_names();
        vector<int> failed = run_ensemble(pe, oe, vector<int>(), cycle);
        if (pe.shape().first == 0)
            throw_em_error("all realizations failed during initial evaluation");
        if (pest_scenario.get_pestpp_options().get_ies_debug_fail_remainder())
        {
            failed.push_back(0);
        }
        if (failed.size() > 0)
        {
            ss.str("");
            ss << "the following " << failed.size() << " par:obs realization runs failed during evaluation of the initial parameter ensemble:" << endl;
            for (auto ifail : failed)
            {
                ss << pnames[ifail] << ":" << onames[ifail];
                if (ifail + 1 % 5 > 0)
                    ss << endl;
            }
            message(0,ss.str());

        }

		pe.transform_ip(ParameterEnsemble::transStatus::NUM);
	}

	ss.str("");
    if (pest_scenario.get_pestpp_options().get_save_dense())
    {
        ss << file_manager.get_base_filename();
        if (cycle != NetPackage::NULL_DA_CYCLE)
            ss << "." << cycle;
        ss << ".0.obs" << dense_file_ext;
        oe.to_dense(ss.str());
    }
	else if (pest_scenario.get_pestpp_options().get_save_binary())
	{
		ss << file_manager.get_base_filename();
		if (cycle != NetPackage::NULL_DA_CYCLE)
			ss << "." << cycle;
		ss << ".0.obs.jcb";
		oe.to_binary(ss.str());
	}
	else
	{
		ss << file_manager.get_base_filename();
		if (cycle != NetPackage::NULL_DA_CYCLE)
			ss << "." << cycle;
		ss << ".0.obs.csv";
		oe.to_csv(ss.str());
	}
	message(1, "saved initial obs ensemble to", ss.str());
	
	//save the 0th iter par and rei and well as the untagged par and rei
	save_real_par_rei(pest_scenario, pe, oe, output_file_writer, file_manager, iter, BASE_REAL_NAME, cycle);
	save_real_par_rei(pest_scenario, pe, oe, output_file_writer, file_manager, -1, BASE_REAL_NAME, cycle);


	performance_log->log_event("calc pre-drop phi");
	//initialize the phi handler
	ph = L2PhiHandler(&pest_scenario, &file_manager, &oe_base, &pe_base, &parcov);

	if (ph.get_lt_obs_names().size() > 0)
	{
		message(1, "less_than inequality defined for observations: ", ph.get_lt_obs_names().size());
	}
	if (ph.get_gt_obs_names().size())
	{
		message(1, "greater_than inequality defined for observations: ", ph.get_gt_obs_names().size());
	}
    map<string,double> t;
    t = ph.get_lt_obs_bounds();
    if (!t.empty())
    {
        ss.str("");
        if (t.size() < 50) {
            ss << "less_than inequality defined through 'less_than' data for observations:" << endl;
            for (const auto it: t) {
                ss << it.first << "," << it.second << endl;
            }
            ss << endl;
        }
        else
        {
            ss << "less_than inequality defined through 'less_than' data for " << t.size() << " observations" << endl;
        }
        message(1,ss.str());
    }
    t = ph.get_gt_obs_bounds();
    if (!t.empty())
    {
        ss.str("");
        if (t.size() < 50) {
            ss << "greater_than inequality defined through 'greater_than' data for observations:" << endl;
            for (const auto it: t) {
                ss << it.first << "," << it.second << endl;
            }
            ss << endl;
        }
        message(1,ss.str());
    }
    else
    {
        ss << "greater_than inequality defined through 'greater_than' data for " << t.size() << " observations" << endl;

    }
    t.clear();
    map<string,pair<double,double>> tt = ph.get_double_obs_bounds();
    if (!tt.empty())
    {
        ss.str("");
        if (tt.size() < 50) {
            ss << "double inequality defined through 'greater_than' and 'less_than' data for observations:" << endl;
            for (const auto it: tt) {
                ss << it.first << "," << it.second.first << " to " << it.second.second << endl;
            }
            ss << endl;
        }
        else
        {
            ss << "double inequality defined through 'greater_than' and 'less_than' data for " << tt.size() << " observations" << endl;

        }
        message(1,ss.str());
    }

	ph.update(oe, pe, weights);
	message(0, "pre-drop initial phi summary");
	ph.report(true);
	
	pcs = ParChangeSummarizer(&pe_base, &file_manager, &output_file_writer);
	message(1,"checking for prior-data conflict");
    vector<string> in_conflict = ph.detect_simulation_data_conflict(oe,".0.pdc.csv");
	if (in_conflict.size() > 0)
	{
		ss.str("");
		ss << "WARNING: " << in_conflict.size() << " non-zero weighted observations are in conflict";
		ss << " with the prior simulated ensemble." << endl;
		message(0, ss.str());

		cout << "...see " << file_manager.get_base_filename() << ".0.pdc.csv" << " for listing of conflicted observations" << endl << endl;
		ofstream& frec = file_manager.rec_ofstream();
        frec << "...see " << file_manager.get_base_filename() << ".0.pdc.csv" << " for listing of conflicted observations" << endl << endl;

		if (!ppo->get_ies_drop_conflicts())
		{
			ss.str("");
			ss << "  WARNING: Prior-data conflict detected.  Continuing with IES parameter" << endl;
			ss << "           adjustment will likely result in parameter and forecast bias." << endl;
			//ss << "           Consider using 'ies_drop_conflicts' as a quick fix.";
			message(1, ss.str());
		}
		else
		{
            ss.str("");
            ss << "  WARNING: 'ies_drop_conflicts' is not an ideal approach.  A better" << endl;
            ss << "           approach would be to dig in and figure out why these observations " << endl;
			ss << "           are in conflict." << endl;
            //check that all obs are in conflict
			message(1, "dropping conflicted observations");

			if (in_conflict.size() == act_obs_names.size()) {
                if (cycle == NetPackage::NULL_DA_CYCLE) {
                    throw_em_error("all non-zero weighted observations in conflict state, cannot continue");
                } else {
                    message(0, "all non-zero weighted observations in conflict state, continuing to next cycle");
                    zero_weight_obs(in_conflict, false, false);
                    ph.update(oe, pe, weights);
                    return;
                }
            }
			if (violation_obs.size() > 0)
            {
			    set<string> sviol_obs(violation_obs.begin(),violation_obs.end());
			    set<string>::iterator end = sviol_obs.end();
			    vector<string> temp,temp2;
			    for (auto& obs : in_conflict)
                {
			        if (sviol_obs.find(obs) == end)
                    {
			            temp.push_back(obs);
                    }
			        else
			        {
			            temp2.push_back(obs);
			        }
                }
                if (temp2.size() > 0)
                {
                    ss.str("");
                    ss << "WARNING: the following observations are 'in conflict' but" << endl;
                    ss << "         also contain 'drop_violation' conditions, " << endl;
                    ss << "         so they are not being treated as conflicts...user beware..." << endl;
                    message(1,ss.str(),temp2);
                    in_conflict = temp;
                }
            }
            if (in_conflict.size() > 0) {
                zero_weight_obs(in_conflict);
                if (ppo->get_ies_localizer().size() > 0) {
                    message(1, "updating localizer");
                    if (localizer.get_use())
                        localizer.get_orgmat_ptr()->clear_names();
                    ofstream &frec = file_manager.rec_ofstream();
                    use_localizer = localizer.initialize(performance_log, frec, true);
                }
            }
		}
	}
    if ((ppo->get_ies_obs_restart_csv().size() > 0) && (ppo->get_ies_par_restart_csv().size() > 0))
    {
        message(1,"not dropping any realizations during initialization b/c restart obs and restart par ensembles passed");
    }
    else {
        if (pest_scenario.get_control_info().noptmax < 0)
        {
            message(1,"not dropping any realizations during initialization b/c noptmax<0");
        }
        else
        {
            drop_bad_reals(pe, oe);
        }
    }
	if (oe.shape().first == 0)
	{
		throw_em_error(string("all realizations dropped as 'bad'"));
	}
	if (oe.shape().first < error_min_reals)
	{
		message(0, "too few active realizations:", oe.shape().first);
		message(1, "need at least ", error_min_reals);
		throw_em_error(string("too few active realizations, cannot continue"));
	}
	if (oe.shape().first < warn_min_reals)
	{
		ss.str("");
		ss << "WARNING: less than " << warn_min_reals << " active realizations...might not be enough";
		string s = ss.str();
		message(0, s);
	}

	adjust_weights(false);

    //need this here since during reinflation we reset weight ensemble and might rebalance weights
    if (ppo->get_ies_phi_fractions_file().size() > 0)
    {
        weights_base = weights;
        org_obs_info = pest_scenario.get_ctl_observation_info_copy();
    }
    if ((ppo->get_ies_phi_fractions_file().size() > 0) ||
        (ppo->get_obscov_filename().size() > 0) ||
        (in_conflict.size() > 0))
    {
        ss.str("");
        if (pest_scenario.get_pestpp_options().get_save_dense())
        {
            ss << file_manager.get_base_filename();
            ss << ".adjusted.weights" << dense_file_ext;
            weights.to_dense(ss.str());
        }
        else if (pest_scenario.get_pestpp_options().get_save_binary())
        {
            ss << file_manager.get_base_filename();
            ss << ".adjusted.weights.jcb";
            weights.to_binary(ss.str());
        }
        else
        {
            ss << file_manager.get_base_filename();
            ss << ".adjusted.weights.csv";
            weights.to_csv(ss.str());
        }
        message(1, "saved adjusted weight ensemble to ", ss.str());
    }
    if ((ppo->get_ies_phi_fractions_file().size() > 0) ||
        (ppo->get_obscov_filename().size() > 0) ||
        (in_conflict.size() > 0))
    {
        map<string,double> base_weights = weights.get_real_map(BASE_REAL_NAME,true);
        string filename = file_manager.get_base_filename() + ".adjusted.obs_data.csv";
        ofstream f_obs(filename);
        if (f_obs.bad())
            throw_em_error("error opening: " + filename);
        output_file_writer.scenario_obs_csv(f_obs, base_weights);
        f_obs.close();


        save_real_par_rei(pest_scenario, pe, oe, output_file_writer, file_manager, iter, BASE_REAL_NAME, cycle, base_weights);
        save_real_par_rei(pest_scenario, pe, oe, output_file_writer, file_manager, -1, BASE_REAL_NAME, cycle, base_weights);


    }


    performance_log->log_event("calc initial phi");
	ph.update(oe, pe, weights);
	message(0, "initial phi summary");
	ph.report(true);
	ph.write(0, run_mgr_ptr->get_total_runs());
	if (ppo->get_ies_save_rescov())
		ph.save_residual_cov(oe, 0);

    best_mean_phis.push_back(ph.get_representative_phi(L2PhiHandler::phiType::COMPOSITE));
    last_best_mean = ph.get_representative_phi(L2PhiHandler::phiType::COMPOSITE);

    if (!pest_scenario.get_pestpp_options().get_ies_use_approx())
	{
		message(1, "using full (MAP) update solution");

	}

	last_best_std = ph.get_std(L2PhiHandler::phiType::COMPOSITE);
	last_best_lam = pest_scenario.get_pestpp_options().get_ies_init_lam();
	bool continue_anyway = false;
	if ((pest_scenario.get_control_info().noptmax > 0) && (act_obs_names.size() > 0))
	{
	    int num_ineq = ph.get_gt_obs_names().size() + ph.get_lt_obs_names().size();
	    int nnz_obs = pest_scenario.get_ctl_ordered_nz_obs_names().size();
		if (ph.get_mean(L2PhiHandler::phiType::ACTUAL) < 1.0e-10)
		{
		    if ((cycle != NetPackage::NULL_DA_CYCLE) && (num_ineq == nnz_obs))
            {
		        message(0,"initial actual phi mean too low but only inequality obs are being used, continuing...");
                continue_anyway = true;
            }
		    else
			    throw_em_error("initial actual phi mean too low, something is wrong...or you have the perfect model that already fits the data shockingly well");
		}
		if (ph.get_std(L2PhiHandler::phiType::ACTUAL) < 1.0e-10)
		{
            if ((cycle != NetPackage::NULL_DA_CYCLE) && (num_ineq == nnz_obs))
            {
                message(0,"initial actual phi stdev is very low but only inequality obs are being used, continuing...");
                continue_anyway = true;
            }
			else
                message(0,"WARNING: initial actual phi stdev is very low, something is probably wrong...");
		}
		if (last_best_mean < 1.0e-10)
		{
            if ((cycle != NetPackage::NULL_DA_CYCLE) && (num_ineq == nnz_obs))
            {
                message(0,"initial composite phi mean too low but only inequality obs are being used, continuing...");
                continue_anyway = true;
            }
			else
                throw_em_error("initial composite phi mean too low, something is wrong...");
		}
		if (last_best_std < 1.0e-10)
		{
            if ((cycle != NetPackage::NULL_DA_CYCLE) && (num_ineq == nnz_obs))
            {
                message(0,"initial composite phi stdev is ver low but only inequality obs are being used, continuing...");
                continue_anyway = true;
            }
			else
            message(0, "WARNING: initial composite phi stdev is very low, something is probably wrong...");
		}
	}

	if (last_best_lam <= 0.0)
	{
        if (continue_anyway)
        {
            last_best_lam = 1000;//?
        }
        else {
            last_best_lam = get_lambda();
        }
	}

	if (cycle != NetPackage::NULL_DA_CYCLE)
	{
		pair<string, string> names = save_ensembles(string(),cycle, pe, oe);
		message(1, "saved cycle prior obs ensemble to", names.second);
		message(1, "saved cycle prior par ensemble to", names.first);
	}
    if (act_obs_names.size() > 0) {
        message(1, "current lambda:", last_best_lam);
    }
	pcs.summarize(pe);

    message(0, "initialization complete");
}

double EnsembleMethod::get_lambda()
{
    //double x = last_best_mean / (2.0 * double(oe.shape().second));
    stringstream ss;
    double lam_val = 0.0;
    double org_val = last_best_lam;
    double x = last_best_mean / (2.0 * double(pest_scenario.get_ctl_ordered_nz_obs_names().size()));
    lam_val = pow(10.0, (floor(log10(x))));
    if (lam_val < 1.0e-10) {
        message(1, "initial lambda estimation from phi failed, using 10,000");
        lam_val = 10000;
    }
    if (org_val < 0.0) {
        org_val *= -1.0;
        ss.str("");
        ss << "scaling phi-based initial lambda: " << last_best_lam
           << ", by user-supplied (negative) initial lambda: " << org_val;
        message(1, ss.str());
        lam_val *= org_val;
    }
    return lam_val;

}

void EnsembleMethod::prep_drop_violations()
{
    violation_obs.clear();
    string section = "observation data external";
    string viol_col = "drop_violations";
    map<string, string> viol_map = pest_scenario.get_ext_file_string_map(section,viol_col);
    if (viol_map.size() == 0)
    {
        return;
    }

    string sval;
    for (auto& v : viol_map)
    {
        sval = pest_utils::strip_cp(pest_utils::lower_cp(v.second));
        if (sval == "true")
        {
            violation_obs.push_back(v.first);
        }
    }
    stringstream ss;
    ss << violation_obs.size() << " 'drop_violations' observations detected, see rec file for listing";
    message(1,ss.str());
    ofstream& f_rec = file_manager.rec_ofstream();
    f_rec << endl << "The following observations are being used as 'drop_violations' constraints:" << endl;

    for (auto& o : violation_obs)
    {
        f_rec << o << endl;
    }
    f_rec << "Note: any 'drop_violation' observations that are not inequalities will be treated as " << endl;
    f_rec << "      equality constraints such that any realization that does not match the observation" << endl;
    f_rec << "      value (numerically) exactly will be dropped. User beware..." << endl << endl;

}


void EnsembleMethod::check_and_fill_phi_factors(map<string,vector<string>>& group_to_obs_map,map<string,vector<string>>& group_map,map<string,map<string,double>>& phi_fracs_by_real,
                                                vector<string>& index, bool check_reals)
{
    stringstream ss;
    string fname = pest_scenario.get_pestpp_options().get_ies_phi_fractions_file();
    if (fname.size() == 0)
    {
        return;
    }
    ss.str("");
    ss << "checking phi factors in file " << fname;
    message(1,ss.str());
    ;
    if (pest_scenario.get_pestpp_options().get_ies_phi_factors_by_real())
    {
        performance_log->log_event("reading 'ies_phi_factors_file for each realization': "+fname);
        phi_fracs_by_real = pest_utils::read_csv_to_nested_map(fname, index);
        performance_log->log_event("checking phi factors for each realizations': "+fname);
        //make sure all factors are strictly positive
        vector<pair<string,string>> problems;
        for (auto& row : phi_fracs_by_real)
        {
            for (auto& entry : row.second)
            {
                if (entry.second == -999.)
                {
                    ss.str("");
                    ss << "phi factor tag '" << row.first << "' has -999 value, indicating that weights for this tag should not be adjusted";
                    message(1,ss.str());

                }
                else if (entry.second <= 0.0)
                {
                    problems.push_back(make_pair(row.first,entry.first));
                }
            }
        }
        if (problems.size() > 0)
        {
            ss.str("");
            ss << "the following row/col positions have phi factors <= 0.0 :";
            for (auto& p : problems)
            {
                ss << p.first << ":" << p.second << ",";
            }
            throw_em_error(ss.str());
        }
        if (check_reals)
        {
            map<string,int> oe_base_real_map = oe_base.get_real_map();
            //base is in the obs+noise en
            if (oe_base_real_map.find(BASE_REAL_NAME) != oe_base_real_map.end())
            {
                //base not in the phi fracs
                if (phi_fracs_by_real.find(BASE_REAL_NAME) == phi_fracs_by_real.end())
                {
                    message(2,"'"+BASE_REAL_NAME+" not in phi factors but in obs+noise ensemble, leaving base realization weights unaltered");

                }
            }
            if (pest_scenario.get_pestpp_options().get_ies_obs_csv().size() == 0)
            {
                if (oe_base_real_map.size() > phi_fracs_by_real.size())
                    throw_em_error("too few phi frac realizations passed");
                map<string,map<string,double>> temp;
                vector<string> oe_base_real_names = oe_base.get_real_names();
                for (int i=0;i<oe_base_real_map.size();i++)
                {
                    temp[oe_base_real_names[i]] = phi_fracs_by_real[index[i]];
                }
                phi_fracs_by_real = temp;

            }
            else {
                vector<string> missing;
                map<string, map<string, double>>::iterator end = phi_fracs_by_real.end();
                for (auto &oereal : oe_base_real_map) {
                    if (phi_fracs_by_real.find(oereal.first) == end) {
                        if (oereal.first != BASE_REAL_NAME)
                            missing.push_back(oereal.first);
                    }
                }
                if (missing.size() > 0)
                {
                    ss.str("");
                    ss << "error in phi factors by realization: the obs+noise ensemble was passed but the realization names ";
                    ss << "in this ensemble do not align with the realziation names in the phi factor realizations";
                    throw_em_error(ss.str());
                }
            }




        }
        //now we need to somehow align the rows of the phi factor table with the realization names
        //first check if they are full coincident

        //if not, check if the obs+noise ensemble was passed, if not, then we can just assuming the user
        // doesn't care about real names
        //but we need to make sure this jives with the base real option
        //but we have to make sure we have at least enough rows in the phi factor table - warn if more than needed,
        //  error if not enough

    }
    else
    {
        performance_log->log_event("reading 'ies_phi_factors_file for all realizations': "+fname);
        map<string,double> phi_fracs = pest_utils::read_twocol_ascii_to_map(fname);
        phi_fracs_by_real[BASE_REAL_NAME] = phi_fracs;
    }


    ObservationInfo* oi = pest_scenario.get_observation_info_ptr();
    set<string> nzgroups;

    for (auto& oname : pest_scenario.get_ctl_ordered_obs_names())
    {
        if (oi->get_weight(oname) == 0)
            continue;
        if (nzgroups.find(oi->get_group(oname)) == nzgroups.end())
        {
            nzgroups.insert(oi->get_group(oname));
            group_to_obs_map[oi->get_group(oname)]= vector<string>();
        }
        group_to_obs_map[oi->get_group(oname)].push_back(oname);

    }

    //check that each nzgroup is not found more than once using the tags in the file
    map<string,vector<string>> rev_group_map;
    for (auto& g : nzgroups)
        rev_group_map[g] = vector<string>();


    map<string,double> phi_fracs;
    for (auto& entry : phi_fracs_by_real) {
        phi_fracs = entry.second;
        vector<string> in_groups;
        for (auto &pf : phi_fracs) {
            in_groups.clear();
            for (auto &g : nzgroups) {

                if (g.find(pest_utils::upper_cp(pf.first)) != string::npos) {
                    rev_group_map.at(g).push_back(pf.first);
                    in_groups.push_back(g);
                }
            }
            group_map[pf.first] = in_groups;
            if (in_groups.size() == 0) {
                message(1, "WARNING: no non-zero obs groups found for tag '" + pf.first + "'");
                continue;
            }

            ss.str("");
            ss << "file tag '" << pf.first << " maps to groups ";
            for (auto &g : group_map.at(pf.first))
                ss << g << ",";
            message(2, ss.str());
            if (pf.second == -999.0)
            {
                ss.str("");
                ss << "phi factor tag '" << pf.first << "' has -999 value, indicating that weights for this tag should not be adjusted";
                message(1,ss.str());
            }
            else if (pf.second <= 0.0) {
                ss.str("");
                ss << "adjust_weights(): phi factor '" << pf.first << "' less or equal 0.0 - this is not allowed";
                throw_em_error(ss.str());
            }
        }

        ss.str("");
        ss << "Errors in phi factors file: ";
        bool has_errors = false;
        for (auto &rg : rev_group_map) {
            if (rg.second.size() == 0) {
                ss << ", group '" << rg.first << "' not identified with any tags " << endl;
                has_errors = true;
            } else if (rg.second.size() > 1) {
                ss << ", group '" << rg.first << "' mapped to multiple tags: ";
                for (auto &tag : rg.second)
                    ss << " " << tag;
                ss << endl;
                has_errors = true;
            }
        }
        if (has_errors) {
            throw_em_error(ss.str());
        }
        //we can break after the first iterate since we have checked that if by reals, all the
        // tags are coherent...
        break;
    }

}

void EnsembleMethod::adjust_weights(bool save) {

    stringstream ss;
    string fname = pest_scenario.get_pestpp_options().get_ies_phi_fractions_file();
    if (fname.size() == 0) {
        return;
    }


    ss.str("");
    ss << "adjusting weights using phi factors in file " << fname;
    message(0, ss.str());

    ph.update(oe, pe, weights);
    message(0, "pre-weight-adjustment initial phi summary");
    ph.report(true);
    //performance_log->log_event("reading 'ies_phi_factors_file': "+fname);
    //map<string,double> phi_fracs = pest_utils::read_twocol_ascii_to_map(fname);
    map<string, map<string, double>> phi_fracs_by_real;
    vector<string> index;
    map<string, vector<string>> group_map, group_to_obs_map;
    check_and_fill_phi_factors(group_to_obs_map, group_map, phi_fracs_by_real, index, true);

    if (phi_fracs_by_real.size() == 1) {
        map<string, double> phi_fracs = phi_fracs_by_real.begin()->second;
        adjust_weights_single(group_to_obs_map, group_map, phi_fracs);
    } else {
        adjust_weights_by_real(group_to_obs_map, group_map, phi_fracs_by_real, index);

    }

    if (save) {
        ss.str("");
        ss << file_manager.get_base_filename() << "." << iter;
        if (pest_scenario.get_pestpp_options().get_save_dense()) {
            ss << ".adjusted.weights" << dense_file_ext;
            weights.to_dense(ss.str());
        } else if (pest_scenario.get_pestpp_options().get_save_binary()) {
            ss << ".adjusted.weights.jcb";
            weights.to_binary(ss.str());
        } else {
            ss << ".adjusted.weights.csv";
            weights.to_csv(ss.str());
        }
        message(1, "saved adjusted weight ensemble to ", ss.str());

        map<string, double> base_weights = weights.get_real_map(BASE_REAL_NAME, true);
        ss.str("");
        ss << file_manager.get_base_filename()  << "." << iter << ".adjusted.obs_data.csv";
        ofstream f_obs(ss.str());
        if (f_obs.bad())
            throw_em_error("error opening: " + ss.str());
        output_file_writer.scenario_obs_csv(f_obs, base_weights);
        f_obs.close();
    }
}


void EnsembleMethod::adjust_weights_by_real(map<string,vector<string>>& group_to_obs_map, map<string,vector<string>>& group_map,
                            map<string,map<string,double>>& phi_fracs_by_real,vector<string> index)
{
    stringstream ss;
    ss.str("");
    map<string,map<string,double>> actual_swr_map = ph.get_swr_real_map(oe, weights);
    map<string,double> phi_fracs;
    map<string,double> current_phi_fracs;
    map<string,double> init_group_phis;
    map<string,map<string,double>> real_init_group_phis,real_adj_group_phis;
    map<string,int> weight_real_map = weights.get_real_map();
    map<string,int> weight_var_map = weights.get_var_map();
    double total = 0;
    double scale_fac = 0;
    double sub_total = 0;
    double cur_mean_phi = 0;
    bool do_adjust = true;
    for (auto& swr_map : actual_swr_map)
    {
        if ((swr_map.first == BASE_REAL_NAME) && (phi_fracs_by_real.find(swr_map.first) == phi_fracs_by_real.end()))
        {
            //just use any phi fracs here because they are not used
            phi_fracs = phi_fracs_by_real.begin()->second;
            do_adjust = false;
        }
        else
        {
            phi_fracs = phi_fracs_by_real.at(swr_map.first);
            do_adjust = true;
        }
        cur_mean_phi = 0.0;
        for (auto& p : swr_map.second)
            cur_mean_phi += p.second;
        current_phi_fracs = swr_map.second;
        init_group_phis.clear();
        total = 0;
        sub_total=0;
        scale_fac = 0;
        for (auto& pf: phi_fracs) {
            if (pf.second == -999)
            {
                ss.str("");
                ss << "NOTE: adjust_weights(): tag " << pf.first << " for realization " << swr_map.first << " has a factor value of -999, skipping this tag";
                message(1,ss.str());
                continue;
            }
            total = 0;
            for (auto &g : group_map.at(pf.first)) {
                sub_total = 0;
                for (auto oname : group_to_obs_map.at(g)) {
                    total += swr_map.second.at(oname);
                    sub_total += swr_map.second.at(oname);
                }
                init_group_phis[g] = sub_total;
            }

            if (total == 0) {
                ss.str("");
                ss << "WARNING: adjust_weights(): tag " << pf.first << " has 0.0 phi";
                message(1, ss.str());
                continue;
            }
            current_phi_fracs[pf.first] = total / cur_mean_phi;
            ss.str("");
            ss << "realization " << swr_map.first << ", file tag '" << pf.first << "' original phi (factor): "
               << total << " (" << current_phi_fracs[pf.first] << ")";
            message(2, ss.str());
            if (!do_adjust)
                scale_fac = 1.0;
            else
                scale_fac = sqrt((cur_mean_phi * pf.second) / total);
            for (auto &g : group_map.at(pf.first)) {

                for (auto oname : group_to_obs_map.at(g)) {
                    //cout << swr_map.first << "," << g << "," << oname << endl;
                    weights.get_eigen_ptr_4_mod()->row(weight_real_map.at(swr_map.first))(weight_var_map.at(oname)) *= scale_fac;
                }
            }
        }

        real_init_group_phis[swr_map.first] = init_group_phis;

    }

    map<string,map<string,double>> adj_swr_map = ph.get_swr_real_map(oe, weights);
    for (auto& swr_map : adj_swr_map) {
        if ((swr_map.first == BASE_REAL_NAME) && (phi_fracs_by_real.find(swr_map.first) == phi_fracs_by_real.end()))
        {
            //just use any phi fracs here because they are not used
            phi_fracs = phi_fracs_by_real.begin()->second;
            do_adjust = false;
        }
        else
        {
            phi_fracs = phi_fracs_by_real.at(swr_map.first);
            do_adjust = true;
        }
        cur_mean_phi = 0.0;
        for (auto &p : swr_map.second)
            cur_mean_phi += p.second;
        current_phi_fracs = swr_map.second;
        init_group_phis.clear();
        total = 0;
        sub_total = 0;
        scale_fac = 0;
        for (auto &pf: phi_fracs) {
            total = 0;
            for (auto &g : group_map.at(pf.first)) {
                sub_total = 0;
                for (auto oname : group_to_obs_map.at(g)) {
                    total += swr_map.second.at(oname);
                    sub_total += swr_map.second.at(oname);
                }
                init_group_phis[g] = sub_total;
            }
            current_phi_fracs[pf.first] = total / cur_mean_phi;
            ss.str("");
            ss << "realization " << swr_map.first << ", file tag '" << pf.first << "' adjusted phi (factor): "
               << total << " (" << current_phi_fracs[pf.first] << ")";
            message(2, ss.str());

        }
        real_adj_group_phis[swr_map.first] = init_group_phis;

    }
    if (verbose_level > 2)
    {
        ss.str("");
        ss << file_manager.get_base_filename() << "." << iter << ".obsgroupadj.summary.csv";

        ofstream f(ss.str());
        if (!f.good())
            throw_em_error("error opening weight adjustment summary file "+ss.str());

        f << "realization,group,initial_phi,adjusted_phi" << endl;
        for (auto& grp_entry : real_init_group_phis) {
            for (auto &g : init_group_phis) {
                f << grp_entry.first << "," << g.first << "," << g.second << "," << real_adj_group_phis.at(grp_entry.first).at(g.first) << endl;
            }
        }
        f.close();
        message(2,"saved obs group weight adjustment summary to "+ss.str());
    }

}

void EnsembleMethod::adjust_weights_single(map<string,vector<string>>& group_to_obs_map, map<string,vector<string>>& group_map,
                                           map<string,double>& phi_fracs)
{
    ObservationInfo *oi = pest_scenario.get_observation_info_ptr();
    stringstream ss;
    ss.str("");
    string center_on = pest_scenario.get_pestpp_options().get_ies_center_on();
    if (center_on.size() != 0)
        message(1,"using realization "+center_on+" residuals for reweighting (if available");
    else
        message(1,"using mean residuals for reweighting");
    L2PhiHandler::phiType ptype = L2PhiHandler::phiType::MEAS;
    map<string,double> mean_swr_map = ph.get_swr_map(oe, center_on);
    double cur_mean_phi = 0;
    for (auto& sw : mean_swr_map)
        cur_mean_phi += sw.second;

    message(1,"original mean phi: ",cur_mean_phi);
    //if the current is really low, just return and the traps in initialize() will catch it.
    if (cur_mean_phi < 1.0e-30)
    {
        performance_log->log_event("mean phi too low - returning");
        return;
    }
    //cur_mean_phi = nzobs_obs_fac * cur_mean_phi;
    map<string,double> current_phi_fracs;
    map<string,double> init_group_phis,adj_group_phis;
    map<string,int> weight_var_map = weights.get_var_map();
    double total = 0;
    double scale_fac = 0;
    double sub_total = 0;
    for (auto& pf: phi_fracs)
    {
        if (pf.second == -999)
        {
            ss.str("");
            ss << "NOTE: adjust_weights(): tag " << pf.first << " has a factor value of -999, skipping this tag";
            message(1,ss.str());
            continue;
        }
        total = 0;
        for (auto& g : group_map.at(pf.first))
        {
            sub_total = 0;
            for (auto oname : group_to_obs_map.at(g))
            {
                total += mean_swr_map.at(oname);
                sub_total += mean_swr_map.at(oname);
            }
            init_group_phis[g] = sub_total;
        }

        if (total == 0)
        {
            ss.str("");
            ss << "WARNING: adjust_weights(): tag " << pf.first << " has 0.0 phi";
            message(1,ss.str());
            continue;
        }
        current_phi_fracs[pf.first] = total / cur_mean_phi;
        ss.str("");
        ss << "file tag '" << pf.first << "' original mean phi (factor): " << total << " (" << current_phi_fracs[pf.first] << ")";
        message(1,ss.str());
        scale_fac = sqrt((cur_mean_phi * pf.second) / total);
        for (auto& g : group_map.at(pf.first))
        {
            for (auto oname : group_to_obs_map.at(g))
            {
                oi->set_weight(oname,oi->get_weight(oname) * scale_fac );
                weights.get_eigen_ptr_4_mod()->col(weight_var_map.at(oname)) *= scale_fac;
            }
        }
    }
    ph.update(oe,pe);
    mean_swr_map = ph.get_swr_map(oe, center_on);
    for (auto& pf: phi_fracs) {
        total = 0;
        for (auto &g : group_map.at(pf.first)) {
            sub_total = 0;
            for (auto oname : group_to_obs_map.at(g)) {
                total += mean_swr_map.at(oname);
                sub_total += mean_swr_map.at(oname);
            }
            adj_group_phis[g] = sub_total;

        }
        current_phi_fracs[pf.first] = total / cur_mean_phi;
        ss.str("");
        ss << "file tag '" << pf.first << "' adjusted mean phi (factor): " << total << " ("
           << total / cur_mean_phi << ")";
        message(1, ss.str());
    }
    if (verbose_level > 2)
    {
        ss.str("");
        ss << file_manager.get_base_filename() << "." << iter << ".obsgroupadj.summary.csv";

        ofstream f(ss.str());
        if (!f.good())
            throw_em_error("error opening weight adjustment summary file "+ss.str());

        f << "group,initial_phi,adjusted_phi" << endl;
        for (auto& g : init_group_phis)
        {
            f << g.first << "," << g.second << "," << adj_group_phis.at(g.first) << endl;
        }
        f.close();
        message(2,"saved obs group weight adjustment summary to "+ss.str());
    }
}

void EnsembleMethod::transfer_dynamic_state_from_oe_to_initial_pe(ParameterEnsemble& _pe, ObservationEnsemble& _oe)
{
	if (obs_dyn_state_names.size() > 0)
	{
		map<string, int> par2col_map;
		for (int i = 0; i < par_dyn_state_names.size(); i++)
			par2col_map[par_dyn_state_names[i]] = i;

		ParameterEnsemble::transStatus org_status = _pe.get_trans_status();
		ParamTransformSeq bts = pest_scenario.get_base_par_tran_seq();
		Eigen::MatrixXd mat = _oe.get_eigen(vector<string>(), obs_dyn_state_names);
		stringstream ss;
		ss << "transferring " << obs_dyn_state_names.size() << " dynamic states from obs to par ensemble for " << _pe.shape().first << " realizations";
		message(1, ss.str());
		if (org_status == ParameterEnsemble::transStatus::NUM)
		{
			//get a vec of adj states
			vector<string> adj_par_dyn_state_names = pest_scenario.get_ctl_ordered_adj_par_names();
			set<string> sadj_pars(adj_par_dyn_state_names.begin(), adj_par_dyn_state_names.end());
			adj_par_dyn_state_names.clear();
			set<string>::iterator end = sadj_pars.end();
			for (auto& p : par_dyn_state_names)
				if (sadj_pars.find(p) != end)
					adj_par_dyn_state_names.push_back(p);
			for (int i = 0; i < mat.rows(); i++)
			{
				Parameters pars = pest_scenario.get_ctl_parameters();
				pars.update(par_dyn_state_names, mat.row(i));
				bts.ctl2numeric_ip(pars);
				for (auto& p : adj_par_dyn_state_names)
					mat(i, par2col_map[p]) = pars.get_rec(p);
			}
		}
		_pe.replace_col_vals_and_fixed(par_dyn_state_names, mat);
	}
}

void EnsembleMethod::transfer_dynamic_state_from_oe_to_final_pe(ParameterEnsemble& _pe, ObservationEnsemble& _oe)
{
    if (final2init_par_state_names.size() > 0) {
        stringstream ss;
        ss << "transferring " << final2init_par_state_names.size();
        ss << " simulated states from obs ensemble to final dynamic state par ensemble for ";
        ss << _pe.shape().first << " realizations, see rec file for listing";
        message(1, ss.str());

        //first get a flipped par-to-par map
        map<string,string> init_to_final_par;
        for (auto& m : final2init_par_state_names)
        {
            init_to_final_par[m.second] = m.first;
        }
        //now build a list of final par state names ordered against obs state names
        //obs_dyn_state_names in sync'd with par_dyn_state_names
        vector<string> final_par_dyn_state_names;
		vector<string> final_obs_dyn_state_names;

		for (int i = 0; i < obs_dyn_state_names.size(); i++)
		{
			if (init_to_final_par.find(par_dyn_state_names[i]) != init_to_final_par.end())
			{
				final_par_dyn_state_names.push_back(init_to_final_par.at(par_dyn_state_names[i]));
				final_obs_dyn_state_names.push_back(obs_dyn_state_names[i]);

			}
		}

        map<string,int> par2col_map;
        for (int i=0;i<final_par_dyn_state_names.size();i++)
            par2col_map[final_par_dyn_state_names[i]] = i;

        ParameterEnsemble::transStatus org_status = _pe.get_trans_status();
        ParamTransformSeq bts = pest_scenario.get_base_par_tran_seq();
        Eigen::MatrixXd mat = _oe.get_eigen(vector<string>(), final_obs_dyn_state_names);

        if (org_status == ParameterEnsemble::transStatus::NUM)
        {
            //get a vec of adj states
            vector<string> adj_par_dyn_state_names = pest_scenario.get_ctl_ordered_adj_par_names();
            set<string> sadj_pars(adj_par_dyn_state_names.begin(), adj_par_dyn_state_names.end());
            adj_par_dyn_state_names.clear();
            set<string>::iterator end = sadj_pars.end();
            for (auto& p : final_par_dyn_state_names)
                if (sadj_pars.find(p) != end)
                    adj_par_dyn_state_names.push_back(p);
            for (int i = 0; i < mat.rows(); i++)
            {
                Parameters pars = pest_scenario.get_ctl_parameters();
                pars.update(final_par_dyn_state_names, mat.row(i));
                bts.ctl2numeric_ip(pars);
                for (auto& p : adj_par_dyn_state_names)
                    mat(i, par2col_map[p]) = pars.get_rec(p);
            }
        }
        _pe.replace_col_vals_and_fixed(final_par_dyn_state_names, mat);

    }
}

void EnsembleMethod::transfer_par_dynamic_state_final_to_initial_ip(ParameterEnsemble& _pe)
{
    if (final2init_par_state_names.size() > 0) {
        stringstream ss;
        ss << "transferring " << final2init_par_state_names.size() << " final dynamic par states to initial dynamic par states for "
           << _pe.shape().first << " realizations, see rec file for listing";
        message(1, ss.str());

        _pe.update_var_map();
        vector<string> final_names,init_names;
        for (auto &sm : final2init_par_state_names) {
            final_names.push_back(sm.first);
            init_names.push_back(sm.second);
        }
        _pe.transform_ip(ParameterEnsemble::transStatus::CTL);
        Eigen::MatrixXd vals = _pe.get_eigen(vector<string>(),final_names);
        _pe.replace_col_vals_and_fixed(init_names,vals);
        _pe.transform_ip(ParameterEnsemble::transStatus::NUM);

    }
}

void EnsembleMethod::initialize_dynamic_states(bool rec_report)
{
	stringstream ss;
	// find a list of dynamic states
	obs_dyn_state_names.clear();
	par_dyn_state_names.clear();
	final2init_par_state_names.clear();
	vector<string> obs_names = pest_scenario.get_ctl_observations().get_keys();
	vector<string> par_names = pest_scenario.get_ctl_parameters().get_keys();
	set<string> spar_names(par_names.begin(), par_names.end());
	set<string>::iterator end = spar_names.end();
	par_names.clear();
	//for (int i = 0; i < obs_names.size(); i++)	
	double w;
	for (auto& name : obs_names)
	{

		if (spar_names.find(name) != end)
		{
			obs_dyn_state_names.push_back(name);
			par_dyn_state_names.push_back(name);
		}

	}
	if (obs_dyn_state_names.size() > 0)
	{
		ss.str("");
		ss << obs_dyn_state_names.size() << " dynamic states identified through shared-names";
		message(1, ss.str());
	}
	map<string, string> state_map = pest_scenario.get_ext_file_string_map("observation data external", "state_par_link");
	if (state_map.size() > 0)
	{

		set<string> pstates(par_dyn_state_names.begin(), par_dyn_state_names.end());
		set<string>::iterator send = pstates.end();
		vector<string> t = pest_scenario.get_ctl_ordered_par_names();
		set<string> pnames(t.begin(), t.end());
		t = pest_scenario.get_ctl_ordered_obs_names();
		set<string> onames(t.begin(), t.end());
		t.clear();
		set<string>::iterator pend = pnames.end();
		set<string>::iterator oend = onames.end();

		vector<string> dups, missing;
		int c = 0;
		for (auto& sm : state_map)
		{
			if (onames.find(sm.first) == oend)
			{
				continue;
			}
			if (pstates.find(sm.second) != send)
			{
				dups.push_back(sm.second);
			}
			else if (pnames.find(sm.second) == pend)
			{
				missing.push_back(sm.second);
			}
			else
			{
                obs_dyn_state_names.push_back(sm.first);
                par_dyn_state_names.push_back(sm.second);
                c++;
			}
		}
		if (dups.size() > 0)
		{
			ss.str("");
			ss << "the following state parameters nominated thru obs data linking " << endl;
			ss << "    were already tagged as 'states' by identically named observations:" << endl;
			for (auto& d : dups)
				ss << d << ",";
			throw_em_error(ss.str());
		}
		if (missing.size() > 0)
		{
            ss.str("");
			ss << "the following parameters nominated thru obs data linking " << endl;
			ss << "    were not found in par data section:" << endl;
			for (auto& m : missing)
				ss << m << ",";
			throw_em_error(ss.str());
		}
		if (c > 0)
		{
			ss.str("");
			ss << c << " dynamic states identified through observation data 'state_par_link'";
			message(1, ss.str());
		}
	}

    map<string, string> par2par_state_map = pest_scenario.get_ext_file_string_map("parameter data external", "state_par_link");
    if (par2par_state_map.size() > 0)
    {
        set<string> pstates(par_dyn_state_names.begin(), par_dyn_state_names.end());
        set<string>::iterator send = pstates.end();
        vector<string> t = pest_scenario.get_ctl_ordered_adj_par_names();
        set<string> adj_pnames(t.begin(), t.end());
        set<string>::iterator adj_end = adj_pnames.end();
        vector<string> dups, missing,already,not_adj;
        set<string> named;
        set<string>::iterator par_send = spar_names.end();
        int c = 0;
        for (auto& sm : par2par_state_map)
        {
            if ((spar_names.find(sm.first) == par_send) && (spar_names.find(sm.second) == par_send))
            {
                continue;
            }
            //if the linking par isn't in current par state names, that's a problem
            if (pstates.find(sm.second) == send)
            {
                missing.push_back(sm.second);
            }
            //if the par is in the current par state names, that's a problem
            else if (pstates.find(sm.first) != send)
            {
                already.push_back(sm.first);
            }
            else if (named.find(sm.second) != named.end())
            {
                dups.push_back(sm.second);
            }
            else
            {
                final2init_par_state_names[sm.first] = sm.second;
                named.emplace(sm.second);
                c++;
                if (adj_pnames.find(sm.first) == adj_end)
                {
                    not_adj.push_back(sm.first);
                }
            }

        }
        if (dups.size() > 0)
        {
            ss.str("");
            ss << "the following final state parameters identified with non-null parameter data 'state_par_link' " << endl;
            ss << "    were listed more than once:" << endl;
            for (auto& d : dups)
                ss << d << ",";
            throw_em_error(ss.str());
        }
        if (missing.size() > 0)
        {
            ss.str("");
            ss << "the following initial state parameters nominated thru parameter data 'state_par_link' " << endl;
            ss << "    were not found in obs-to-par dynamic state linkage:" << endl;
            for (auto& m : missing)
                ss << m << ",";
            throw_em_error(ss.str());
        }
        if (already.size() > 0)
        {
            ss.str("");
            ss << "the following initial state parameters nominated thru parameter data 'state_par_link' " << endl;
            ss << "    were listed more than once in parameter data 'state_par_link':" << endl;
            for (auto& a : already)
                ss << a << ",";
            throw_em_error(ss.str());

        }
        if ((!pest_scenario.get_pestpp_options().get_da_use_simulated_states()) && (not_adj.size() > 0))
        {
            ss.str("");
            ss << "not using simulated states and the following final parameter states are not adjustable:" << endl;
            for (auto& p : not_adj)
                ss << p << ",";
            throw_em_error(ss.str());
        }

        ss.str("");
        ss << c << " final-to-initial parameter states identified thru parameter data 'state_par_link' entries";
        message(1,ss.str());

    }
    if ((rec_report) && (obs_dyn_state_names.size() < 10000)) {
        ofstream &frec = file_manager.rec_ofstream();
        frec << "...observation-to-parameter state mapping: " << endl;
        for (int i = 0; i < obs_dyn_state_names.size(); i++) {
            frec << "observation state '" << obs_dyn_state_names[i] << "' maps to parameter initial state '"
                 << par_dyn_state_names[i] << "'" << endl;
        }
    }
    if ((rec_report) && (final2init_par_state_names.size() < 10000)) {
        ofstream &frec = file_manager.rec_ofstream();
        frec << endl << "...final-to-initial parameter state mapping: " << endl;
        for (auto &sm : final2init_par_state_names) {
            frec << "final state '" << sm.first << "' maps to initial state '" << sm.second << "'" << endl;
        }
    }
}

bool EnsembleMethod::solve_glm(int cycle)
{
	vector<double> lambdas = lam_mults;
	message(1, "current lambda: ", last_best_lam);
	for (auto& m : lambdas)
		m *= last_best_lam;

	vector<double> scale_facs = pest_scenario.get_pestpp_options().get_lambda_scale_vec();
	return solve(false, lambdas, scale_facs,cycle);

}

bool EnsembleMethod::solve_mda(bool last_iter,int cycle)
{

	//this function should cover the case where noptmax = 1 (vanilla) also...
	vector<double> mda_facs, scaled_mda_facs;
	mda_facs.push_back(pest_scenario.get_pestpp_options().get_ies_mda_init_fac());
	double tot, ttot;
	if (iter == 1)	{
		
		tot = 1.0 / mda_facs[0];
		double dec_fac = pest_scenario.get_pestpp_options().get_ies_mda_dec_fac();// get_ies_mda_dec_fac();
		
		ttot = 0.0;
		for (int i = 1; i < pest_scenario.get_control_info().noptmax; i++)
		{
			mda_facs.push_back(mda_facs[i - 1] * dec_fac);
			tot += (1.0 / mda_facs[i]);
		}		
		scaled_mda_facs.clear();
		for (auto& mda_fac : mda_facs)
		{
			mda_fac = (1.0 / mda_fac) / tot;
			scaled_mda_facs.push_back(1.0 / mda_fac);
			ttot += mda_fac;
		}
		mda_lambdas = scaled_mda_facs;
	}	
	else 
	{
		if (last_iter) // trim unused lambdas after the last iteration
		{
			mda_lambdas.erase(mda_lambdas.begin()+iter, mda_lambdas.end());
		}
		int ii = 1;
		double tot_fac1, tot_fac2, resid_fac;
		tot_fac1 = 0;
		tot_fac2 = 0;
		scaled_mda_facs.clear();		
		mda_lambdas[iter - 1] = last_best_lam;
		tot = 0;
		ttot = 0;
		for (auto& mda_fac : mda_lambdas)
		{		
			tot += (1.0 / mda_fac);
			if (ii<iter)
				tot_fac1 += (1.0 / mda_fac);
			else
				tot_fac2 += (1.0 / mda_fac);
			ii += 1;			
		}
		
		ii = 1;
		for (auto& mda_fac : mda_lambdas)
		{
			if (ii < iter)
			{
				ii += 1;
				ttot += 1.0/mda_fac;
				scaled_mda_facs.push_back(mda_fac);
				continue;
			}
			mda_fac = (1.0-tot_fac1)*(1.0 / mda_fac) / tot_fac2;
			scaled_mda_facs.push_back(1.0 / mda_fac);
			ttot += mda_fac;
			ii += 1;
		}
		mda_lambdas = scaled_mda_facs;

	}
	
	mda_facs = vector<double>{ mda_lambdas[iter-1] };
	vector<double> scale_facs = pest_scenario.get_pestpp_options().get_lambda_scale_vec();
	return solve(true, mda_facs, scale_facs,cycle);
	
}

bool EnsembleMethod::solve(bool use_mda, vector<double> inflation_factors, vector<double> backtrack_factors, int cycle)
{
	stringstream ss;
	ofstream& frec = file_manager.rec_ofstream();
	if (pe.shape().first < error_min_reals)
	{
		message(0, "too few active realizations:", oe.shape().first);
		message(1, "need at least ", error_min_reals);
		throw_em_error(string("too few active realizations, cannot continue"));
	}
	if (pe.shape().first < warn_min_reals)
	{
		ss.str("");
		ss << "WARNING: less than " << warn_min_reals << " active realizations...might not be enough";
		string s = ss.str();
		message(1, s);
	}

	int local_subset_size = pest_scenario.get_pestpp_options().get_ies_subset_size();
	if (local_subset_size < 0)
    {
	    ss.str("");

	    local_subset_size = (int)((double)pe.shape().first) * ((-1. * (double)local_subset_size) / 100.);

        ss << "subset defined as a percentage of ensemble size, using " << local_subset_size;
        ss << " realizations for subset" << endl;
        message(2,ss.str());
        if (local_subset_size < 4)
        {
            ss.str("");
            ss << "percentage-based subset size too small, increasing to 4" << endl;
            local_subset_size = 4;
            message(2,ss.str());
        }
    }

	if ((use_subset) && (local_subset_size > pe.shape().first))
	{
		ss.str("");
		ss << "subset size (" << local_subset_size << ") greater than ensemble size (" << pe.shape().first << ")";
		frec << "  ---  " << ss.str() << endl;
		cout << "  ---  " << ss.str() << endl;
		frec << "  ...reducing subset size to " << pe.shape().first << endl;
		cout << "  ...reducing subset size to " << pe.shape().first << endl;
		local_subset_size = pe.shape().first;
	}

	else if ((inflation_factors.size() == 1) && (backtrack_factors.size() == 1))
	{
		ss.str("");
		ss << "only testing one inflation/lambda factor and one scale/backtrack factor, not using subset";
		frec << "  ---  " << ss.str() << endl;
		cout << "  ---  " << ss.str() << endl;
		local_subset_size = pe.shape().first;
	}

	pe.transform_ip(ParameterEnsemble::transStatus::NUM);

	vector<ParameterEnsemble> pe_lams;
	vector<double> lam_vals, scale_vals;
	//update all the fast-lookup structures
	performance_log->log_event("reordering variables in pe");
	pe.reorder(vector<string>(), act_par_names);
	oe.update_var_map();
	pe.update_var_map();
	parcov.update_sets();
	obscov.update_sets();

	//make sure weight ensemble is in sync with oe
	if (weights.shape().first != oe.shape().first)
    {
	    if (weights.shape().first < oe.shape().first)
        {
            throw_em_error("weight ensemble has less realizations that obs ensemble");
        }
	    ss.str("");
	    ss << "resizing weight ensemble from " << weights.shape().first << " realizations to " << oe.shape().first << " realizations";
	    performance_log->log_event(ss.str());
	    weights.keep_rows(oe.get_real_names());
    }

	//build up this container here and then reuse it for each lambda later...
	unordered_map<string, pair<vector<string>, vector<string>>> loc_map;
	if (use_localizer)
	{

		loc_map = localizer.get_localanalysis_case_map(iter, act_obs_names, act_par_names, oe, pe, performance_log, frec);
	}
	else
	{
		pair<vector<string>, vector<string>> p(act_obs_names, act_par_names);
		loc_map["all"] = p;
	}
	if (loc_map.size() == 0)
	{
		throw_em_error("EnsembleMethod::solve() internal error: loc_map is empty");
	}
	//get this once and reuse it for each lambda
	Eigen::MatrixXd Am;
	if ((!use_mda) && (!pest_scenario.get_pestpp_options().get_ies_use_approx()))
	{
		Am = get_Am(pe.get_real_names(), pe.get_var_names());
	}

	vector<int> subset_idxs = get_subset_idxs(pe.shape().first, local_subset_size);
	vector<string> pe_filenames;

	performance_log->log_event("preparing EnsembleSolver");
	ParameterEnsemble pe_upgrade(pe.get_pest_scenario_ptr(), &rand_gen, pe.get_eigen(vector<string>(), act_par_names, false), pe.get_real_names(), act_par_names);
	pe_upgrade.set_zeros();
	pe_upgrade.set_trans_status(pe.get_trans_status());
	ObservationEnsemble oe_upgrade(oe.get_pest_scenario_ptr(), &rand_gen, oe.get_eigen(vector<string>(), act_obs_names, false), oe.get_real_names(), act_obs_names);
    EnsembleSolver es(performance_log, file_manager, pest_scenario, pe, oe_upgrade, oe_base, weights, localizer, parcov, Am, ph,
		use_localizer, iter, act_par_names, act_obs_names, reg_factor);
    double mm_alpha = pest_scenario.get_pestpp_options().get_ies_multimodal_alpha();
    if (mm_alpha != 1.0)
    {
        es.update_multimodal_components(mm_alpha);
    }


    //solve for each factor
    for (auto& cur_lam : inflation_factors)
	{
		ss.str("");
		if (!use_mda)
			message(1, "starting calcs for glm factor", cur_lam);
		else
			message(1, "starting calcs for mda factor", cur_lam);
		message(2, "see .log file for more details");

		pe_upgrade = ParameterEnsemble(pe.get_pest_scenario_ptr(), &rand_gen, pe.get_eigen(vector<string>(), act_par_names, false), pe.get_real_names(), act_par_names);
        pe_upgrade.set_zeros();
		pe_upgrade.set_trans_status(pe.get_trans_status());

		if (mm_alpha != 1.0)
        {
            message(1,"multimodal solve for inflation factor ",cur_lam);
            es.solve_multimodal(num_threads, cur_lam, !use_mda, pe_upgrade, loc_map, mm_alpha);
        }
		else{
            es.solve(num_threads, cur_lam, !use_mda, pe_upgrade, loc_map);
		}

		map<string, double> norm_map;
		for (auto sf : backtrack_factors)
		{
			ParameterEnsemble pe_lam_scale = pe;
			pe_lam_scale.set_eigen(*pe_lam_scale.get_eigen_ptr() + (*pe_upgrade.get_eigen_ptr() * sf));
			if (pest_scenario.get_pestpp_options().get_ies_enforce_bounds())
			{
				if (pest_scenario.get_pestpp_options().get_ies_enforce_chglim())
					norm_map = pe_lam_scale.enforce_change_limits_and_bounds(performance_log, pe);
				else
					norm_map = pe_lam_scale.enforce_bounds(performance_log, false);

				ss.str("");
				ss << " lambda " << cur_lam << ", scale factor " << sf;
				norm_map_report(norm_map, ss.str());
			}
			ss.str("");
			ss << file_manager.get_base_filename() << "." << iter << "." << cur_lam << ".lambda." << sf << ".scale.par";

			if ((!pest_scenario.get_pestpp_options().get_ies_upgrades_in_memory()) && (subset_idxs.size() < pe.shape().first) && ((inflation_factors.size() > 1) || (backtrack_factors.size() > 1)))
			{
				pe_lam_scale.to_dense(ss.str() + ".bin");
				pe_filenames.push_back(ss.str() + ".bin");
				pe_lam_scale.keep_rows(subset_idxs,true);
				message(1,"upgrade ensemble saved to " + ss.str() + ".bin");

			}
			else if (!pest_scenario.get_pestpp_options().get_ies_upgrades_in_memory())
			{
				message(1, "even though 'ies_upgrades_in_memory' is 'false', there is no benefit to using this option because either you aren't testing multiple upgrades or you aren't using a subset");
			}

			pe_lams.push_back(pe_lam_scale);
			lam_vals.push_back(cur_lam);
			scale_vals.push_back(sf);
			if (!pest_scenario.get_pestpp_options().get_ies_save_lambda_en())
				continue;

            if (pest_scenario.get_pestpp_options().get_save_dense())
            {
                pe_lam_scale.to_dense_unordered(ss.str()+dense_file_ext);
            }
			else if (pest_scenario.get_pestpp_options().get_save_binary())
			{
				pe_lam_scale.to_binary(ss.str()+".jcb");
			}
			else
			{
				pe_lam_scale.to_csv(ss.str()+".csv");
			}
			frec << "lambda, scale value " << cur_lam << ',' << sf << " par ensemble saved to " << ss.str() << endl;

		}

		ss.str("");
		message(1, "finished calcs for:", cur_lam);

	}

	if (pest_scenario.get_pestpp_options().get_ies_debug_upgrade_only())
	{
		message(0, "ies_debug_upgrade_only is true, exiting");
		if (pe_filenames.size() > 0)
		{
			message(1, "first attempting to load upgrade ensemble " + pe_filenames[0]);
			performance_log->log_event("loading dense pe upgrade ensemble");
			ParameterEnsemble remaining_pe_lam(&pest_scenario);
			remaining_pe_lam.from_binary(pe_filenames[0]);
			remaining_pe_lam.transform_ip(ParameterEnsemble::transStatus::NUM);
		}
		throw_em_error("ies_debug_upgrade_only is true, exiting");
	}

	//the case for all state estimation and non-iterative
	//only one upgrade lambda and all pars are states
	if (((pe_lams.size() == 1) && (par_dyn_state_names.size() == pe_lams[0].shape().second) && pest_scenario.get_control_info().noptmax == 1))
	{
		message(1, "non-iterative state-estimation detected, not evaluating update ensemble for current cycle");
		pe = pe_lams[0];
		//move the estimated states to the oe, which will then later be transferred back to the pe
		ph.update(oe, pe, weights);
        double best_mean = ph.get_representative_phi(L2PhiHandler::phiType::COMPOSITE);
		double best_std = ph.get_std(L2PhiHandler::phiType::COMPOSITE);
		best_mean_phis.push_back(best_mean);

		return true;
	}

	vector<map<int, int>> real_run_ids_lams;
	int best_idx = -1;
	double best_mean = 1.0e+300, best_std = 1.0e+300; // todo (Ayman): read those from input
	double mean, std;

	message(0, "running upgrade ensembles");
	vector<ObservationEnsemble> oe_lams;
	
	//if we are saving upgrades to disk
	if (pe_filenames.size() > 0)
	{
		vector<int> temp;
		for (int i = 0; i < subset_idxs.size(); i++)
			temp.push_back(i);
		oe_lams = run_lambda_ensembles(pe_lams, lam_vals, scale_vals, cycle, temp, subset_idxs);
	}
	else
 		oe_lams = run_lambda_ensembles(pe_lams, lam_vals, scale_vals, cycle, subset_idxs,subset_idxs);

	message(0, "evaluating upgrade ensembles");
	message(1, "last mean: ", last_best_mean);
	message(1, "last stdev: ", last_best_std);

    int q = pest_utils::quit_file_found();
    if ((q == 1) || (q == 2))
    {
        message(1,"'pest.stp' found, quitting");
    }
    else if (q == 4)
    {
        message(0,"pest.stp found with '4'.  run mgr has returned control, removing file.");
        if (!pest_utils::try_remove_quit_file())
        {
            message(0,"error removing pest.stp file, bad times ahead...");
        }
    }

    double acc_fac = pest_scenario.get_pestpp_options().get_ies_accept_phi_fac();
    double lam_inc = pest_scenario.get_pestpp_options().get_ies_lambda_inc_fac();
    double lam_dec = pest_scenario.get_pestpp_options().get_ies_lambda_dec_fac();
	ObservationEnsemble oe_lam_best(&pest_scenario);
	bool echo = false;
	if (verbose_level > 1)
		echo = true;
	for (int i = 0; i < pe_lams.size(); i++)
	{
		if (oe_lams[i].shape().first == 0)
			continue;
		vector<double> vals({ lam_vals[i],scale_vals[i] });
		if (pest_scenario.get_pestpp_options().get_ies_save_lambda_en())

		{
			ss.str("");
			ss << file_manager.get_base_filename() << "." << iter << "." << lam_vals[i] << ".lambda." << scale_vals[i] << ".scale.obs";

            if (pest_scenario.get_pestpp_options().get_save_dense())
            {
                ss << dense_file_ext;
                oe_lams[i].to_dense(ss.str());
            }
			else if (pest_scenario.get_pestpp_options().get_save_binary())
			{
				ss << ".jcb";
				oe_lams[i].to_binary(ss.str());
			}
			else
			{
				ss << ".csv";
				oe_lams[i].to_csv(ss.str());
			}
			frec << "lambda, scale value " << lam_vals[i] << ',' << scale_vals[i] << " obs ensemble saved to " << ss.str() << endl;

		}
        drop_bad_reals(pe_lams[i], oe_lams[i], subset_idxs);
		if (oe_lams[i].shape().first == 0)
		{
			message(1, "all realizations dropped as 'bad' for lambda, scale fac ", vals);
			continue;
		}

		ph.update(oe_lams[i], pe_lams[i], weights);

		message(0, "phi summary for lambda, scale fac:", vals, echo);
		ph.report(echo);

        mean = ph.get_representative_phi(L2PhiHandler::phiType::COMPOSITE);

		std = ph.get_std(L2PhiHandler::phiType::COMPOSITE);
        ph.write_lambda(iter,oe_lams[i].shape().first,last_best_lam,last_best_mean,
                        last_best_std,
                        scale_vals[i],lam_vals[i],mean,std);
		if (mean < best_mean)
		{
			oe_lam_best = oe_lams[i];
			best_mean = mean;
			best_std = std;
			best_idx = i;
		}
	}
	if (best_idx == -1)
	{
		message(0, "WARNING:  unsuccessful upgrade testing, multiplying lambda by 10000.0 and returning to upgrade calculations");
		last_best_lam *= 10000.0;
		return false;

	}

    if ((best_idx != -1) && (use_subset) && (local_subset_size < pe.shape().first)) {
        double acc_phi = last_best_mean * acc_fac;

        if (pest_scenario.get_pestpp_options().get_ies_debug_high_subset_phi()) {
            cout << "ies_debug_high_subset_phi active" << endl;
            frec << "ies_debug_high_subset_phi active" << endl;
            best_mean = acc_phi + 1.0;
        }

        if ((best_mean > acc_phi)) {

            ss.str("");
            ss << "best subset mean phi  (" << best_mean << ") greater than acceptable phi : " << acc_phi;
            string m = ss.str();
            message(0, m);
            ph.update(oe, pe, weights);
            best_mean_phis.push_back(ph.get_representative_phi(L2PhiHandler::phiType::COMPOSITE));
            if (!use_mda) {
                message(1, "updating realizations with reduced phi");
                update_reals_by_phi(pe_lams[best_idx], oe_lams[best_idx], subset_idxs);
            }

            ph.update(oe, pe, weights);
            //re-check phi
            double new_best_mean;
            new_best_mean = ph.get_representative_phi(L2PhiHandler::phiType::COMPOSITE);
            if (new_best_mean < best_mean) {
                best_mean = new_best_mean;
                best_std = ph.get_std(L2PhiHandler::phiType::COMPOSITE);
                //replace the last entry in the best mean phi tracker
                best_mean_phis[best_mean_phis.size() - 1] = best_mean;
                message(1, "current best mean phi (after updating reduced-phi reals): ", best_mean);
                if (best_mean < last_best_mean * acc_fac) {
                    if (best_std < last_best_std * acc_fac) {
                        double new_lam = lam_vals[best_idx] * lam_dec;
                        new_lam = (new_lam < lambda_min) ? lambda_min : new_lam;
                        message(0, "partial update improved phi stats, updating lambda to ", new_lam);
                        last_best_lam = new_lam;
                    } else {
                        message(0, "not updating lambda (standard deviation reduction criteria not met)");
                    }
                    last_best_std = best_std;
                }
            } else {
                double new_lam = last_best_lam * lam_inc;
                new_lam = (new_lam > lambda_max) ? lambda_max : new_lam;
                last_best_lam = new_lam;
                message(1,
                        "abandoning current upgrade ensembles, returning to upgrade calculations and increasing lambda to ",
                        new_lam);

            }
            message(1, "returning to upgrade calculations...");

            return false;
        }

        //release the memory of the unneeded pe_lams
        for (int i = 0; i < pe_lams.size(); i++) {
            if (i == best_idx)
                continue;
            pe_lams[i] = pe.zeros_like(0);
        }
        //need to work out which par and obs en real names to run - some may have failed during subset testing...
        ObservationEnsemble remaining_oe_lam = oe;//copy

        ParameterEnsemble remaining_pe_lam = pe_lams[best_idx];
        remaining_pe_lam.set_fixed_info(pe.get_fixed_info());

        if (pe_filenames.size() > 0) {
            performance_log->log_event(
                    "'ies_upgrades_in_memory' is 'false', loading 'best' parameter ensemble from file '" +
                    pe_filenames[best_idx] + "'");
            vector<string> missing;
            if (oe_lams[best_idx].shape().first != remaining_pe_lam.shape().first) {
                set<string> ssub_names;
                for (auto real_name: oe_lams[best_idx].get_real_names())
                    ssub_names.emplace(real_name);
                vector<string> oreal_names = oe.get_real_names();
                vector<string> preal_names = pe.get_real_names();

                for (auto idx: subset_idxs)
                    if (ssub_names.find(oreal_names[idx]) == ssub_names.end())
                        missing.push_back(preal_names[idx]);
                if (missing.size() > 0)
                    pe_lams[best_idx].drop_rows(missing);
            }

            remaining_pe_lam.from_binary(pe_filenames[best_idx]);
            remaining_pe_lam.transform_ip(ParameterEnsemble::transStatus::NUM);
            //remove any failed runs from subset testing
            if (missing.size() > 0)
                remaining_pe_lam.drop_rows(missing);
            remove_external_pe_filenames(pe_filenames);

        }


        vector<string> pe_keep_names, oe_keep_names;
        vector<string> pe_names = pe.get_real_names(), oe_names = oe.get_real_names();

        vector<string> org_pe_idxs, org_oe_idxs;
        set<string> ssub;
        for (auto &i: subset_idxs)
            ssub.emplace(pe_names[i]);
        for (int i = 0; i < pe_names.size(); i++)
            if (ssub.find(pe_names[i]) == ssub.end()) {
                pe_keep_names.push_back(pe_names[i]);
                //oe_keep_names.push_back(oe_names[i]);
            }
        ssub.clear();
        for (auto &i: subset_idxs)
            ssub.emplace(oe_names[i]);
        for (int i = 0; i < oe_names.size(); i++)
            if (ssub.find(oe_names[i]) == ssub.end()) {
                oe_keep_names.push_back(oe_names[i]);
            }
        message(0, "phi summary for best lambda, scale fac: ",
                vector<double>({lam_vals[best_idx], scale_vals[best_idx]}));
        ph.update(oe_lams[best_idx], pe_lams[best_idx], weights);
        ph.report(true, false);
        message(0, "running remaining realizations for best lambda, scale:",
                vector<double>({lam_vals[best_idx], scale_vals[best_idx]}));

        //pe_keep_names and oe_keep_names are names of the remaining reals to eval
        performance_log->log_event("dropping subset idxs from remaining_pe_lam");
        remaining_pe_lam.keep_rows(pe_keep_names);
        performance_log->log_event("dropping subset idxs from remaining_oe_lam");
        remaining_oe_lam.keep_rows(oe_keep_names);
        //save these names for later
        org_pe_idxs = remaining_pe_lam.get_real_names();
        org_oe_idxs = remaining_oe_lam.get_real_names();
        ///run
        vector<int> fails = run_ensemble(remaining_pe_lam, remaining_oe_lam, vector<int>(), cycle);
        int q = pest_utils::quit_file_found();
        if ((q == 1) || (q == 2)) {
            message(1, "'pest.stp' found, quitting");
            return true;
        }
        else if (q == 4)
        {
            message(0,"pest.stp found with '4'.  run mgr has returned control, removing file.");
            if (!pest_utils::try_remove_quit_file())
            {
                message(0,"error removing pest.stp file, bad times ahead...");
            }
        }
        //for testing
        if (pest_scenario.get_pestpp_options().get_ies_debug_fail_remainder()) {
            ss.str("");
            ss << "ies_debug_fail_remainder is True, failing  par:obs realization " << org_pe_idxs[0] << ":"
               << org_oe_idxs[0];
            message(0, ss.str());
            fails.push_back(0);
        }

        //if any of the remaining runs failed
        if (fails.size() == org_pe_idxs.size())
        {
            double new_lam = last_best_lam * lam_inc;
            new_lam = (new_lam > lambda_max) ? lambda_max : new_lam;
            last_best_lam = new_lam;
            message(1, "all remaining realizations failed...something is prob wrong, returning to upgrade calculations and increasing lambda to ",
                    new_lam);
            return false;
        }
		if (fails.size() > 0)
		{

			vector<string> new_pe_idxs, new_oe_idxs;
			vector<int>::iterator start = fails.begin(), end = fails.end();
			stringstream ss;
			ss << "the following par:obs realizations failed during evaluation of the remaining ensemble: ";
			for (int i = 0; i < org_pe_idxs.size(); i++)
				if (find(start, end, i) == end)
				{
					new_pe_idxs.push_back(org_pe_idxs[i]);
					new_oe_idxs.push_back(org_oe_idxs[i]);
				}
				else
				{
					ss << org_pe_idxs[i] << ":" << org_oe_idxs[i] << " , ";
				}
			string s = ss.str();
			message(1, s);
			remaining_oe_lam.keep_rows(new_oe_idxs);
			remaining_pe_lam.keep_rows(new_pe_idxs);

		}
		//drop the remaining runs from the par en then append the remaining par runs (in case some failed)
		performance_log->log_event("assembling ensembles");
		pe_lams[best_idx].drop_rows(pe_keep_names);
		pe_lams[best_idx].append_other_rows(remaining_pe_lam);
		pe_lams[best_idx].set_fixed_info(pe.get_fixed_info());
		//append the remaining obs en
		oe_lam_best.append_other_rows(remaining_oe_lam);
		assert(pe_lams[best_idx].shape().first == oe_lam_best.shape().first);
        drop_bad_reals(pe_lams[best_idx], oe_lam_best);
		if (oe_lam_best.shape().first == 0)
		{
			//throw_em_error(string("all realization dropped after finishing subset runs...something might be wrong..."));
            double new_lam = last_best_lam * lam_inc;
            new_lam = (new_lam > lambda_max) ? lambda_max : new_lam;
            last_best_lam = new_lam;
            message(1, "all realization dropped after finishing subset runs...something is prob wrong, returning to upgrade calculations and increasing lambda to ",
                    new_lam);
            return false;
		}
		performance_log->log_event("updating phi");
		ph.update(oe_lam_best, pe_lams[best_idx], weights);

        best_mean = ph.get_representative_phi(L2PhiHandler::phiType::COMPOSITE);
		best_std = ph.get_std(L2PhiHandler::phiType::COMPOSITE);
		message(1, "phi summary for entire ensemble using lambda,scale_fac ", vector<double>({ lam_vals[best_idx],scale_vals[best_idx] }));
		ph.report(true, false);
	}
	else
	{
		ph.update(oe_lam_best, pe_lams[best_idx], weights);
        best_mean = ph.get_representative_phi(L2PhiHandler::phiType::COMPOSITE);
		best_std = ph.get_std(L2PhiHandler::phiType::COMPOSITE);
	}

	ph.update(oe_lam_best, pe_lams[best_idx], weights);
    best_mean = ph.get_representative_phi(L2PhiHandler::phiType::COMPOSITE);
	best_std = ph.get_std(L2PhiHandler::phiType::COMPOSITE);
	message(1, "last best mean phi * acceptable phi factor: ", last_best_mean * acc_fac);
	message(1, "current best mean phi: ", best_mean);

	if (pest_scenario.get_pestpp_options().get_ies_debug_high_upgrade_phi())
	{
		cout << "ies_debug_high_upgrade_phi active" << endl;
        frec << "ies_debug_high_upgrade_phi active" << endl;
		best_mean = (last_best_mean * acc_fac) + 1.0;
	}

	//track this here for phi-based termination check
	best_mean_phis.push_back(best_mean);


	if ((best_mean < last_best_mean * acc_fac))
	{
		message(0, "updating parameter ensemble");
		performance_log->log_event("updating parameter ensemble");
		last_best_mean = best_mean;
        if (pest_scenario.get_pestpp_options().get_ies_updatebyreals())
        {

            if (pe.shape().first > pe_lams[best_idx].shape().first)
            {
                performance_log->log_event("aligning rows of pe with pe_lams[best_idx]");
                pe.keep_rows(pe_lams[best_idx].get_real_names());
            }
            if (oe.shape().first > oe_lam_best.shape().first) {
                performance_log->log_event("aligning rows of oe with oe_lam_best");
                oe.keep_rows(oe_lam_best.get_real_names());
            }
            message(0, "only updating realizations with reduced phi");
            update_reals_by_phi(pe_lams[best_idx], oe_lam_best);
        }
        else {
            pe = pe_lams[best_idx];
            oe = oe_lam_best;
        }
		if (best_std < last_best_std * acc_fac)
		{
			double new_lam = lam_vals[best_idx] * lam_dec;
			new_lam = (new_lam < lambda_min) ? lambda_min : new_lam;
			message(0, "updating lambda to ", new_lam);
			last_best_lam = new_lam;
		}
		else
		{
			message(0, "not updating lambda (standard deviation reduction criteria not met)");
		}
		last_best_std = best_std;
	}

	else
	{
		//message(0, "not updating parameter ensemble");
		if ((!use_mda) && (!pest_scenario.get_pestpp_options().get_ies_updatebyreals()))
		{
			message(0, "only updating realizations with reduced phi");
			update_reals_by_phi(pe_lams[best_idx], oe_lam_best);
		}
		ph.update(oe, pe, weights);
		//re-check phi
		double new_best_mean;
        new_best_mean = ph.get_representative_phi(L2PhiHandler::phiType::COMPOSITE);
		if (new_best_mean < best_mean)
		{
			best_mean = new_best_mean;

			best_std = ph.get_std(L2PhiHandler::phiType::COMPOSITE);
			//replace the last entry in the best mean phi tracker
			best_mean_phis[best_mean_phis.size() - 1] = best_mean;
			message(1, "current best mean phi (after updating reduced-phi reals): ", best_mean);
			if (best_mean < last_best_mean * acc_fac)
			{
				if (best_std < last_best_std * acc_fac)
				{
					double new_lam = lam_vals[best_idx] * lam_dec;
					new_lam = (new_lam < lambda_min) ? lambda_min : new_lam;
					message(0, "updating lambda to ", new_lam);
					last_best_lam = new_lam;
				}
				else
				{
					message(0, "not updating lambda (standard deviation reduction criteria not met)");
				}
				last_best_std = best_std;
			}
			else
			{
				double new_lam = last_best_lam * lam_inc;
				new_lam = (new_lam > lambda_max) ? lambda_max : new_lam;
				message(0, "increasing lambda to: ", new_lam);
				last_best_lam = new_lam;
			}
		}
        else {
            double new_lam = last_best_lam * lam_inc;
            new_lam = (new_lam > lambda_max) ? lambda_max : new_lam;
            message(0, "increasing lambda to: ", new_lam);
            last_best_lam = new_lam;
        }
        save_ensembles("rejected",cycle,pe_lams[best_idx],oe_lam_best);
	}

	return true;
}

void EnsembleMethod::reinflate_par_ensemble(double reinflate_factor,int reinflate_num_reals){

    string min_phi_name = "";
    //find the min phi real...
    map<string,double> pmap = ph.get_phi_map(L2PhiHandler::phiType::ACTUAL);
    double min_phi = numeric_limits<double>::max();

    for (auto& item : pmap)
    {
        if (item.second < min_phi)
        {
            min_phi = item.second;
            min_phi_name = item.first;
        }
    }
    vector<string> real_names = oe.get_real_names();
    int idx = distance(real_names.begin(),find(real_names.begin(),real_names.end(),min_phi_name));
    if (idx == real_names.size())
    {
        throw_em_error("min phi realization name not found in obs ensemble");
    }
    real_names = pe.get_real_names();
    min_phi_name = real_names[idx];

    message(0,"resetting current parameter ensemble to prior ensemble with current ensemble mean");
    message(1,"reinflation factor:",reinflate_factor);
    performance_log->log_event("getting prior parameter ensemble mean-centered anomalies");
    vector<string> pebase_real_names = pe_base.get_real_names();
	vector<string> oebase_real_names = oe_base.get_real_names();

    if ((reinflate_num_reals != 0) && (abs(reinflate_num_reals) <= pe_base.shape().first))
    {
        vector<string> tpar,tobs;
        bool has_base = false;
        for (auto& n : pebase_real_names)
        {
            if (n == BASE_REAL_NAME)
            {
                has_base = true;
                break;
            }
        }
        bool found_base = false;
        for (int i=0;i<abs(reinflate_num_reals);i++)
        {
            tpar.push_back(pebase_real_names[i]);
        	tobs.push_back(oebase_real_names[i]);
            if (pebase_real_names[i] == BASE_REAL_NAME)
            {
                found_base = true;
            }
        }
        pebase_real_names = tpar;
		oebase_real_names = tobs;
        if ((has_base) && (!found_base))
        {
            pebase_real_names[pebase_real_names.size() -1] = BASE_REAL_NAME;
        	oebase_real_names[oebase_real_names.size() -1] = BASE_REAL_NAME;
        }

    }
	Eigen::MatrixXd anoms;
	if (reinflate_num_reals < 0) {
		message(1,"draw new parameter ensemble from current ensemble of size ",abs(reinflate_num_reals));
		anoms = pe.get_eigen_anomalies();
		anoms *= 1.0/sqrt(static_cast<double>(anoms.rows()));
		Eigen::MatrixXd draws(anoms.rows(),pebase_real_names.size());
		draws.setZero();
		performance_log->log_event("making standard normal draws");

		for (int i = 0; i < pebase_real_names.size(); i++)
		{
			for (int j = 0; j < anoms.rows(); j++)
			{
				draws(j, i) = draw_standard_normal(rand_gen);
			}
		}
		performance_log->log_event("making standard normal draws");
		anoms = (anoms.transpose() * draws).transpose();
		draws.resize(0,0);
		if (abs(reinflate_factor) < 1.0) {
			performance_log->log_event("adding scaled prior anomalies to new ensemble");
			Eigen::MatrixXd pranoms = pe_base.get_eigen_anomalies(pebase_real_names, pe.get_var_names(),
				pest_scenario.get_pestpp_options().get_ies_center_on());
			vector<int> seq = uniform_int_draws(pebase_real_names.size(),0,pebase_real_names.size()-1,rand_gen);
			for (int i=0;i<pebase_real_names.size();i++) {
				if (pebase_real_names[i] == BASE_REAL_NAME) {
					continue;
				}
				anoms.row(i).array() += (abs(reinflate_factor) * pranoms.row(seq[i]).array());
			}
		}
	}
	else {
		anoms = pe_base.get_eigen_anomalies(pebase_real_names, pe.get_var_names(), pest_scenario.get_pestpp_options().get_ies_center_on());
		anoms = anoms * abs(reinflate_factor);
	}

    performance_log->log_event("getting current parameter ensemble mean vector");
    vector<double> mean_vec = pe.get_mean_stl_var_vector();
    Eigen::VectorXd offset(mean_vec.size());
    for (int i=0;i<mean_vec.size();i++)
        offset[i] = mean_vec[i];
    if (reinflate_to_minphi_real)
    {
        offset = pe.get_real_vector(min_phi_name);
        message(2,"using min-phi realization for offset");
    }

    performance_log->log_event("adding offset to anomalies");
    for (int i = 0; i < mean_vec.size(); i++)
    {
        anoms.col(i) = anoms.col(i).array() + offset[i];
    }
    performance_log->log_event("forming new parameter ensemble of mean-shifted realizations");
    ParameterEnsemble new_pe = ParameterEnsemble(&pest_scenario,&rand_gen,anoms,pebase_real_names,pe.get_var_names());

    new_pe.set_trans_status(pe.get_trans_status());
    new_pe.set_fixed_info(pe.get_fixed_info());
    if (pest_scenario.get_pestpp_options().get_ies_enforce_bounds()) {
        new_pe.enforce_bounds(performance_log, false);
    }

    message(0,"running new reinflated realizations: ",new_pe.shape().first);
    stringstream ss;
    ss.str("");
    ss << "iteration:" << iter;
    vector<int> temp;
    ofstream& frec = file_manager.rec_ofstream();
    oe.reserve(oebase_real_names,oe.get_var_names());
    weights = weights_base;
    for (auto& oname : oe.get_var_names())
    {
        pest_scenario.get_observation_info_ptr()->get_observation_rec_ptr_4_mod(oname)->weight = org_obs_info.get_weight(oname);
    }
    run_ensemble_util(performance_log,frec,new_pe,oe,run_mgr_ptr,false,temp,NetPackage::NULL_DA_CYCLE, ss.str());
    pe = new_pe;
    new_pe = ParameterEnsemble();
    report_and_save(NetPackage::NULL_DA_CYCLE);
    ph.update(oe,pe,weights);
    message(0,"reinflation phi report:");


    best_mean_phis.push_back(ph.get_representative_phi(L2PhiHandler::phiType::COMPOSITE));
    last_best_mean = ph.get_representative_phi(L2PhiHandler::phiType::COMPOSITE);

    last_best_std = ph.get_std(L2PhiHandler::phiType::COMPOSITE);

    ph.report(true,true);
    ph.write(iter, run_mgr_ptr->get_total_runs());
    ss.str("");
    ss << file_manager.get_base_filename() << "." << iter << ".reinflate.pcs.csv";
    pcs.summarize(pe, ss.str());
    last_best_lam = pest_scenario.get_pestpp_options().get_ies_init_lam();
    double phi_lam = get_lambda();
    last_best_lam = phi_lam;
    message(1,"iter = ies_n_iter_reinflate, resetting lambda to ",last_best_lam);
    consec_bad_lambda_cycles = 0;

}

void EnsembleMethod::remove_external_pe_filenames(vector<string>& pe_filenames)
{
    for (auto& pe_filename : pe_filenames)
    {
        performance_log->log_event("removing upgrade ensemble '" + pe_filename + "'");
        try
        {
            remove(pe_filename.c_str());
        }
        catch (exception& e)
        {
            message(2, "error removing upgrade ensemble: '" + pe_filename + "' :" + e.what());
        }
        catch (...)
        {
            message(2, "error removing upgrade ensemble: '" + pe_filename);
        }

    }
}


void EnsembleMethod::message(int level, const string& _message, vector<string> _extras, bool echo)
{
	stringstream ss;
	if (level == 0)
		ss << endl << "  ---  ";
	else if (level == 1)
		ss << "...";
	ss << _message;
	if (_extras.size() > 0)
	{

        int ee = 0;
        for (auto& e : _extras)
        {
            ss << e << " , ";
            if (ee > 5)
            {
                ss << endl << "     ";
                ee = 0;
            }
            ee++;
        }

	}
	if (level == 0)
		ss << "  ---  ";
	if ((echo) && ((verbose_level >= 2) || (level < 2)))
		cout << ss.str() << endl;
	file_manager.rec_ofstream() << ss.str() << endl;
	performance_log->log_event(_message);

}


void EnsembleMethod::message(int level, const string& _message, vector<int> _extras, bool echo)
{
	stringstream ss;
	if (level == 0)
		ss << endl << "  ---  ";
	else if (level == 1)
		ss << "...";
	ss << _message;
	if (_extras.size() > 0)
	{

        int ee = 0;
        for (auto& e : _extras)
        {
            ss << e << " , ";
            if (ee > 5)
            {
                ss << endl << "     ";
                ee = 0;
            }
            ee++;
        }

	}
	if (level == 0)
		ss << "  ---  ";
	if ((echo) && ((verbose_level >= 2) || (level < 2)))
		cout << ss.str() << endl;
	file_manager.rec_ofstream() << ss.str() << endl;
	performance_log->log_event(_message);

}

void EnsembleMethod::message(int level, const string& _message, vector<double> _extras, bool echo)
{
	stringstream ss;
	if (level == 0)
		ss << endl << "  ---  ";
	else if (level == 1)
		ss << "...";
	ss << _message;

	if (_extras.size() > 0)
	{
        int ee = 0;
		for (auto& e : _extras)
        {
            ss << e << " , ";
            if (ee > 5)
            {
                ss << endl << "     ";
                ee = 0;
            }
            ee++;
        }


	}
	if (level == 0)
		ss << "  ---  ";
	if ((echo) && ((verbose_level >= 2) || (level < 2)))
		cout << ss.str() << endl;
	file_manager.rec_ofstream() << ss.str() << endl;
	performance_log->log_event(_message);

}

void EnsembleMethod::message(int level, const string& _message)
{
	message(level, _message, vector<string>());
}

//template<typename T>

void EnsembleMethod::message(int level, const string& _message, string extra)
{
	stringstream ss;
	ss << _message << " " << extra;
	string s = ss.str();
	message(level, s);
}

void EnsembleMethod::message(int level, const string& _message, int extra)
{
	stringstream ss;
	ss << _message << " " << extra;
	string s = ss.str();
	message(level, s);
}

void EnsembleMethod::message(int level, const string& _message, double extra)
{
	stringstream ss;
	ss << _message << " " << extra;
	string s = ss.str();
	message(level, s);
}

void EnsembleMethod::message(int level, const string& _message, size_t extra)
{
	stringstream ss;
	ss << _message << " " << extra;
	string s = ss.str();
	message(level, s);
}


bool EnsembleMethod::initialize_pe(Covariance& cov)
{
	stringstream ss;
	int num_reals = pest_scenario.get_pestpp_options().get_ies_num_reals();
	string par_csv = pest_scenario.get_pestpp_options().get_ies_par_csv();
	bool drawn = false;
	if (par_csv.size() == 0)
	{
		ofstream& frec = file_manager.rec_ofstream();
		message(1, "drawing parameter realizations: ", num_reals);
		map<string, double> par_means = pest_scenario.get_ext_file_double_map("parameter data external", "mean");
		Parameters draw_par = pest_scenario.get_ctl_parameters();
		if (par_means.size() > 0)
		{

			frec << "Note: the following parameters contain 'mean' value information that will be used in place of " << endl;
			frec << "      the 'parval1' values as mean values during ensemble generation" << endl;
			double lb, ub;
			for (auto par_mean : par_means)
			{
				if (draw_par.find(par_mean.first) != draw_par.end())
				{
					lb = pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(par_mean.first)->lbnd;
					ub = pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(par_mean.first)->ubnd;
					if (par_mean.second < lb)
					{
						frec << "Warning: 'mean' value for parameter " << par_mean.first << " less than lower bound, using 'parval1'";
					}
					else if (par_mean.second > ub)
					{
						frec << "Warning: 'mean' value for parameter " << par_mean.first << " greater than upper bound, using 'parval1'";
					}
					else
					{
						draw_par[par_mean.first] = par_mean.second;
						frec << par_mean.first << " " << par_mean.second << endl;
					}

				}

			}
		}

		if ((pest_scenario.get_pestpp_options().get_ies_enforce_bounds()) && (pest_scenario.get_pestpp_options().get_ies_enforce_chglim()))
		{
			ss.str("");
			ss << "WARNING: 'ies_enforce_chglim' is true, so bounds are being enforced on the " << endl;
			ss << "          initial parameter ensemble by shrinking each realization towards the " << endl;
			ss << "          parameter values in the pest control file. If any parameters are at / near " << endl;
			ss << "          bounds in the control file this will cause many realizations to be nearly " << endl;
			ss << "          identical to the parameter values in the control file" << endl;
			message(1, ss.str());

		}

		map<string, double> norm_map = pe.draw(num_reals, draw_par, cov, performance_log, pest_scenario.get_pestpp_options().get_ies_verbose_level(), file_manager.rec_ofstream());
		norm_map_report(norm_map, "initial parameter realizations");
		drawn = true;
		
	}
	else
	{
		string par_ext = pest_utils::lower_cp(par_csv).substr(par_csv.size() - 3, par_csv.size());
		performance_log->log_event("processing par ensemble " + par_csv);
		if (par_ext.compare("csv") == 0)
		{
			message(1, "loading par ensemble from csv file", par_csv);
			try
			{
				pe.from_csv(par_csv);
			}
			catch (const exception& e)
			{
				ss << "error processing par csv: " << e.what();
				throw_em_error(ss.str());
			}
			catch (...)
			{
				throw_em_error(string("error processing par csv"));
			}
		}
		else if ((par_ext.compare("jcb") == 0) || (par_ext.compare("jco") == 0) || (par_ext.compare("bin")==0))
		{
			message(1, "loading par ensemble from binary file", par_csv);
			try
			{
				pe.from_binary(par_csv);
			}
			catch (const exception& e)
			{
				ss << "error processing par jcb: " << e.what();
				throw_em_error(ss.str());
			}
			catch (...)
			{
				throw_em_error(string("error processing par jcb"));
			}
		}
		else
		{
			ss << "unrecognized par ensemble extension " << par_ext << ", looking for csv, jcb, jco, or bin";
			throw_em_error(ss.str());
		}

		pe.transform_ip(ParameterEnsemble::transStatus::NUM);

		if (pp_args.find("IES_NUM_REALS") != pp_args.end())
		{
			int num_reals = pest_scenario.get_pestpp_options().get_ies_num_reals();

			if (num_reals < pe.shape().first)
			{
				message(1, "ies_num_reals arg passed, truncated parameter ensemble to ", num_reals);
				vector<string> keep_names, real_names = pe.get_real_names();
				for (int i = 0; i < num_reals; i++)
				{
					keep_names.push_back(real_names[i]);
				}
				pe.keep_rows(keep_names);
			}
		}


		if (pest_scenario.get_pestpp_options().get_ies_use_empirical_prior())
		{
			message(1, "initializing prior parameter covariance matrix from ensemble (using diagonal matrix)");
			parcov = pe.get_diagonal_cov_matrix();
			if (pest_scenario.get_pestpp_options().get_ies_verbose_level() > 1)
			{


				if (pest_scenario.get_pestpp_options().get_save_binary())
				{
					string filename = file_manager.get_base_filename() + ".prior.jcb";
					message(1, "saving empirical parameter covariance matrix to binary file: ", filename);
					parcov.to_binary(filename);
				}
				else
				{
					string filename = file_manager.get_base_filename() + ".prior.cov";
					message(1, "saving empirical parameter covariance matrix to ASCII file: ", filename);
					parcov.to_ascii(filename);
				}
			}
		}
		if (pest_scenario.get_pestpp_options().get_ies_enforce_bounds())
		{
			if (pest_scenario.get_pestpp_options().get_ies_obs_restart_csv().size() > 0)
				message(1, "Warning: even though ies_enforce_bounds is true, a restart obs en was passed, so bounds will not be enforced on the initial par en");
			else
			{

				//dont use the shrinking bounds - too many chances for zero length
				map<string, double> norm_map = pe.enforce_bounds(performance_log, false);
				norm_map_report(norm_map, "initial parameter");


			}

		}

	}
	
	return drawn;

}

void EnsembleMethod::norm_map_report(map<string, double>& norm_map, string tag, double thres)
{
	stringstream ss;
	ss.str("");
	ss << "WARNING: the following " << tag << " realizations have been shrunk to less than 10% of their original length: " << endl;
	int count = 0;
	for (auto& n : norm_map)
		if (n.second < thres)
		{
			ss << n.first << ": " << n.second * 100.0 << " % original length" << endl;
			count++;
		}
	if (count > 0)
	{
		file_manager.rec_ofstream() << ss.str() << endl;
		ss.str("");
		ss << "WARNING: " << count << " " << tag << " have been shrunk to less than 10% of their original length, see .rec file for listing" << endl;
		message(1, ss.str());
	}

}

void EnsembleMethod::add_bases()
{
    stringstream ss;
	//check that 'base' isn't already in ensemble
	vector<string> rnames = pe.get_real_names();
	bool inpar = false;
	if (find(rnames.begin(), rnames.end(), BASE_REAL_NAME) != rnames.end())
	{
		message(1, "'base' realization already in parameter ensemble, ignoring '++ies_include_base'");
		inpar = true;
	}
	else
	{
		message(1, "adding 'base' parameter values to ensemble");
		Parameters pars = pest_scenario.get_ctl_parameters();
		pe.get_par_transform().active_ctl2numeric_ip(pars);
		vector<int> drop{ pe.shape().first - 1 };
		pe.drop_rows(drop,true);
		pe.append(BASE_REAL_NAME, pars);
	}

	//check that 'base' isn't already in ensemble
	rnames = oe_base.get_real_names();
	vector<string> wrnames = weights.get_real_names();
	bool in_weight = false;
    if (find(wrnames.begin(),wrnames.end(),BASE_REAL_NAME) != wrnames.end()) {
        in_weight = true;
    }
	if (find(rnames.begin(), rnames.end(), BASE_REAL_NAME) != rnames.end())
	{
		message(1, "'base' realization already in observation ensemble, ignoring 'include_base'");
		if (in_weight)
        {
		    int odist = distance(rnames.begin(),find(rnames.begin(), rnames.end(), BASE_REAL_NAME));
            int wdist = distance(wrnames.begin(),find(wrnames.begin(), wrnames.end(), BASE_REAL_NAME));
            if (odist != wdist)
            {
                ss.str("");
                ss << "'base' realization found at different location in weight ensemble (" << wdist;
                ss << ") than in obs ensemble (" << odist << ")";
                throw_em_error(ss.str());
            }
        }
	}
	else
	{
		Observations obs = pest_scenario.get_ctl_observations();
		Observations wobs;
        wobs.clear();
        ObservationInfo* oi = pest_scenario.get_observation_info_ptr();
        for (auto& oname : weights.get_var_names())
        {
            wobs[oname] = oi->get_weight(oname);
        }

		if (inpar)
		{
			vector<string> prnames = pe.get_real_names();

			int idx = find(prnames.begin(), prnames.end(), BASE_REAL_NAME) - prnames.begin();
			//cout << idx << "," << rnames.size() << endl;
			string oreal = rnames[idx];
			stringstream ss;
			ss << "warning: 'base' realization in par ensenmble but not in obs ensemble," << endl;
			ss << "         replacing obs realization '" << oreal << "' with 'base'";
			string mess = ss.str();
			message(1, mess);
			vector<string> drop;
			drop.push_back(oreal);
			oe_base.drop_rows(drop,true);
			oe_base.append(BASE_REAL_NAME, obs);
			//rnames.insert(rnames.begin() + idx, string(base_name));
			rnames[idx] = BASE_REAL_NAME;
			oe_base.reorder(rnames, vector<string>());
			if (!in_weight)
            {
                string oreal = wrnames[idx];
                ss.str("");
                ss << "warning: 'base' realization in par ensemble but not in weight ensemble," << endl;
                ss << "         replacing weight realization '" << oreal << "' with 'base' weights";
                string mess = ss.str();
                message(1, mess);
                vector<string> drop;
                drop.push_back(oreal);
                weights.drop_rows(drop,true);

                weights.append(BASE_REAL_NAME, wobs);
                //rnames.insert(rnames.begin() + idx, string(base_name));
                wrnames[idx] = BASE_REAL_NAME;
                weights.reorder(wrnames, vector<string>());
            }

		}
		else
		{
			message(1, "adding 'base' observation values to ensemble");
			vector<int> drop{ oe_base.shape().first - 1 };
			oe_base.drop_rows(drop);
			oe_base.append(BASE_REAL_NAME, obs);
            message(1, "adding 'base' weight values to weight ensemble");
            vector<int> wdrop{ weights.shape().first - 1 };
            weights.drop_rows(wdrop);
            weights.append(BASE_REAL_NAME, wobs);
		}
	}
}

bool EnsembleMethod::initialize_weights()
{
	stringstream ss;
    if (act_obs_names.size() == 0)
    {
        int i = 0;
        Eigen::VectorXd wvec(oe_base.shape().second);
        for (auto& n : oe_base.get_var_names())
        {
            wvec[i] = 0.0;
            i++;
        }
        int num_reals = oe_base.shape().first;
        message(1, "setting weights ensemble to zeros - no non-zero weighted obs found");
        weights.reserve(oe_base.get_real_names(), oe_base.get_var_names());
        weights.get_eigen_ptr_4_mod()->setZero();
        for (int i = 0; i < weights.shape().first; i++) {
            weights.get_eigen_ptr_4_mod()->row(i) = wvec;
        }
        return true;

    }


	string weight_csv = pest_scenario.get_pestpp_options().get_ies_weights_csv();
	bool drawn = false;
    Eigen::VectorXd wvec(act_obs_names.size());

    int i = 0;
    for (auto& n : act_obs_names)
    {
        wvec[i] = pest_scenario.get_observation_info_ptr()->get_weight(n);
        i++;
    }

	if (weight_csv.size() == 0) {
        int num_reals = oe_base.shape().first;
        message(1, "setting weights ensemble from control file weights");
        weights.reserve(oe_base.get_real_names(), act_obs_names);
        weights.get_eigen_ptr_4_mod()->setZero();
        for (int i = 0; i < weights.shape().first; i++) {
            weights.get_eigen_ptr_4_mod()->row(i) = wvec;
        }
        drawn = true;
    }
	else
	{
		string obs_ext = pest_utils::lower_cp(weight_csv).substr(weight_csv.size() - 3, weight_csv.size());
		performance_log->log_event("processing weights ensemble " + weight_csv);
		if (obs_ext.compare("csv") == 0)
		{
			message(1, "loading weights ensemble from csv file", weight_csv);
			try
			{
				weights.from_csv(weight_csv);
			}
			catch (const exception& e)
			{
				ss << "error processing weights csv: " << e.what();
				throw_em_error(ss.str());
			}
			catch (...)
			{
				throw_em_error(string("error processing weights csv"));
			}
		}
		else if ((obs_ext.compare("jcb") == 0) || (obs_ext.compare("jco") == 0) || (obs_ext.compare("bin")==0))
		{
			message(1, "loading weights ensemble from binary file", weight_csv);
			try
			{
				weights.from_binary(weight_csv);
			}
			catch (const exception& e)
			{
				stringstream ss;
				ss << "error processing weights binary file: " << e.what();
				throw_em_error(ss.str());
			}
			catch (...)
			{
				throw_em_error(string("error processing weights binary file"));
			}
		}
		else
		{
			ss << "unrecognized weights ensemble extension " << obs_ext << ", looking for csv, jcb, jco, or bin";
			throw_em_error(ss.str());
		}
		//make sure all oe realizations are in weights
		vector<string> missing, keep, real_names = oe_base.get_real_names();
		map<string,int> rmap = weights.get_real_map();
		map<string,int>::iterator end = rmap.end();
		for (auto& rname : real_names)
        {
		    if (rmap.find(rname) == end)
		        missing.push_back(rname);
        }
		if (missing.size() >0)
        {
		    if ((missing.size() == 1) && (missing[0] == BASE_REAL_NAME) && (pest_scenario.get_pestpp_options().get_ies_include_base()))
            {
		        message(2,"'base' realization missing from weight ensemble, replacing last realizations with 'base'");
                Eigen::VectorXd wvec_full(weights.get_var_names().size());
                i = 0;
                for (auto& n : weights.get_var_names()) {
                    wvec_full[i] = pest_scenario.get_observation_info_ptr()->get_weight(n);
                    i++;
                }
                vector<int> drop{ weights.shape().first - 1 };
                weights.drop_rows(drop);
                weights.append(BASE_REAL_NAME, wvec_full);
            }
		    else {
                ss.str("");
                ss << "weights ensemble missing obs ensemble realizations: " << endl << "   ";
                for (auto &m : missing) {
                    ss << m << ",";
                }
                throw_em_error(ss.str());
            }
        }
		real_names = weights.get_real_names();
		rmap = oe_base.get_real_map();
		end = rmap.end();
		for (auto& rname : real_names)
        {
		    if (rmap.find(rname) != end)
		        keep.push_back(rname);
        }
		weights.keep_rows(keep);
	}
	return drawn;

}


bool EnsembleMethod::initialize_oe(Covariance& cov)
{
    stringstream ss;
    int num_reals = pe.shape().first;
    string obs_csv = pest_scenario.get_pestpp_options().get_ies_obs_csv();
    bool drawn = false;
    if (obs_csv.size() == 0)
    {
        if ((pest_scenario.get_pestpp_options().get_ies_no_noise()) || (act_obs_names.size() == 0))
        {
            message(1, "initializing no-noise observation ensemble of : ", num_reals);
            oe_base.initialize_without_noise(num_reals);

        }
        else
        {
            message(1, "drawing observation noise realizations: ", num_reals);
            oe_base.draw(num_reals, cov, performance_log, pest_scenario.get_pestpp_options().get_ies_verbose_level(), file_manager.rec_ofstream());

        }
        drawn = true;
    }
    else
    {
        string obs_ext = pest_utils::lower_cp(obs_csv).substr(obs_csv.size() - 3, obs_csv.size());
        performance_log->log_event("processing obs ensemble" + obs_csv);
        if (obs_ext.compare("csv") == 0)
        {
            message(1, "loading obs ensemble from csv file", obs_csv);
            try
            {
                oe_base.from_csv(obs_csv);
            }
            catch (const exception& e)
            {
                ss << "error processing obs csv: " << e.what();
                throw_em_error(ss.str());
            }
            catch (...)
            {
                throw_em_error(string("error processing obs csv"));
            }
        }
        else if ((obs_ext.compare("jcb") == 0) || (obs_ext.compare("jco") == 0) || (obs_ext.compare("bin")==0))
        {
            message(1, "loading obs ensemble from binary file", obs_csv);
            try
            {
                oe_base.from_binary(obs_csv);
            }
            catch (const exception& e)
            {
                stringstream ss;
                ss << "error processing obs binary file: " << e.what();
                throw_em_error(ss.str());
            }
            catch (...)
            {
                throw_em_error(string("error processing obs binary file"));
            }
        }
        else
        {
            ss << "unrecognized obs ensemble extension " << obs_ext << ", looking for csv, jcb, jco, bin";
            throw_em_error(ss.str());
        }
        if (pp_args.find("IES_NUM_REALS") != pp_args.end())
        {
            int num_reals = pest_scenario.get_pestpp_options().get_ies_num_reals();

            if (num_reals < oe_base.shape().first)
            {
                message(1, "ies_num_reals arg passed, truncated observation ensemble to ", num_reals);
                vector<string> keep_names, real_names = oe_base.get_real_names();
                for (int i = 0; i < num_reals; i++)
                {
                    keep_names.push_back(real_names[i]);
                }
                oe_base.keep_rows(keep_names);
            }
        }
    }
    return drawn;

}


void EnsembleMethod::initialize_restart()
{
	stringstream ss;
	string obs_restart_csv = pest_scenario.get_pestpp_options().get_ies_obs_restart_csv();
	string par_restart_csv = pest_scenario.get_pestpp_options().get_ies_par_restart_csv();

	//performance_log->log_event("restart with existing obs ensemble: " + obs_restart_csv);
	message(1, "restarting with existing obs ensemble", obs_restart_csv);
	string obs_ext = pest_utils::lower_cp(obs_restart_csv).substr(obs_restart_csv.size() - 3, obs_restart_csv.size());
	if (obs_ext.compare("csv") == 0)
	{
		message(1, "loading restart obs ensemble from csv file", obs_restart_csv);
		try
		{
			oe.from_csv(obs_restart_csv);
		}
		catch (const exception& e)
		{
			ss << "error processing restart obs csv: " << e.what();
			throw_em_error(ss.str());
		}
		catch (...)
		{
			throw_em_error(string("error processing restart obs csv"));
		}
	}
	else if ((obs_ext.compare("jcb") == 0) || (obs_ext.compare("jco") == 0) || (obs_ext.compare("bin")==0))
	{
		message(1, "loading restart obs ensemble from binary file", obs_restart_csv);
		try
		{
			oe.from_binary(obs_restart_csv);
		}
		catch (const exception& e)
		{
			ss << "error processing restart obs binary file: " << e.what();
			throw_em_error(ss.str());
		}
		catch (...)
		{
			throw_em_error(string("error processing restart obs binary file"));
		}
	}
	else
	{
		ss << "unrecognized restart obs ensemble extension " << obs_ext << ", looking for csv, jcb, jco, bin";
		throw_em_error(ss.str());
	}
	if (par_restart_csv.size() > 0)
	{
		string par_ext = pest_utils::lower_cp(par_restart_csv).substr(par_restart_csv.size() - 3, par_restart_csv.size());
		if (par_ext.compare("csv") == 0)
		{
			message(1, "loading restart par ensemble from csv file", par_restart_csv);
			try
			{
				pe.from_csv(par_restart_csv);
			}
			catch (const exception& e)
			{
				ss << "error processing restart par csv: " << e.what();
				throw_em_error(ss.str());
			}
			catch (...)
			{
				throw_em_error(string("error processing restart par csv"));
			}
		}
		else if ((par_ext.compare("jcb") == 0) || (par_ext.compare("jco") == 0) || (par_ext.compare("bin")==0))
		{
			message(1, "loading restart par ensemble from binary file", par_restart_csv);
			try
			{
				pe.from_binary(par_restart_csv);
			}
			catch (const exception& e)
			{
				ss << "error processing restart par binary file: " << e.what();
				throw_em_error(ss.str());
			}
			catch (...)
			{
				throw_em_error(string("error processing restart par binary file"));
			}
		}
		else
		{
			ss << "unrecognized restart par ensemble extension " << par_ext << ", looking for csv, jcb, jco, or bin";
			throw_em_error(ss.str());
		}
		if (pe.shape().first != oe.shape().first)
		{
			ss.str("");
			ss << "restart par en has " << pe.shape().first << " realizations but restart obs en has " << oe.shape().first;
			throw_em_error(ss.str());
		}

		//check that restart pe is in sync with pe_base
		vector<string> pe_real_names = pe.get_real_names(), pe_base_real_names = pe_base.get_real_names();
		vector<string>::const_iterator start, end;
		vector<string> missing;
		start = pe_base_real_names.begin();
		end = pe_base_real_names.end();
		for (auto& rname : pe_real_names)
			if (find(start, end, rname) == end)
				missing.push_back(rname);
		if (missing.size() > 0)
		{
			ss << "the following realization names were found in the restart par en but not in the 'base' par en:";
			for (auto& m : missing)
				ss << m << ",";
			throw_em_error(ss.str());
		}
		pe.transform_ip(ParameterEnsemble::transStatus::NUM);
	}

	if (pp_args.find("IES_NUM_REALS") != pp_args.end())
	{
		int num_reals = pest_scenario.get_pestpp_options().get_ies_num_reals();

		if (num_reals < oe.shape().first)
		{
			message(1, "ies_num_reals arg passed, truncated restart obs ensemble to ", num_reals);
			vector<string> keep_names, real_names = oe.get_real_names();
			for (int i = 0; i < num_reals; i++)
			{
				keep_names.push_back(real_names[i]);
			}
			oe.keep_rows(keep_names);
		}
	}

	//check that restart oe is in sync with oe_base
	vector<string> oe_real_names = oe.get_real_names(), oe_base_real_names = oe_base.get_real_names();
	vector<string>::const_iterator start, end;
	vector<string> missing;
	start = oe_base_real_names.begin();
	end = oe_base_real_names.end();
	if (verbose_level > 3)
	    cout << "restart oe real names: " << endl;
	for (auto& rname : oe_real_names)
	{
	    if (verbose_level > 3)
		    cout << rname << endl;
		if (find(start, end, rname) == end)
			missing.push_back(rname);
	}
	if (missing.size() > 0)
	{
		//the special case where the base real is what is missing...
		if ((missing.size() == 1) && (missing[0] == BASE_REAL_NAME))
		{
			//check that the base real is in the par en - restart_par_en should be accounted for by now
			int base_par_idx = -1;
			vector<string> pe_real_names = pe.get_real_names(), pe_base_real_names = pe_base.get_real_names();
			for (int i = 0; i < pe_base.shape().first; i++)
			{
				if (pe_base_real_names[i] == BASE_REAL_NAME)
				{
					base_par_idx = i;
					break;
				}
			}
			if (base_par_idx != -1)
			{

				ss.str("");
				ss << "WARNING: replacing obs+noise obs en realization '" << oe_base_real_names[base_par_idx] << "' with 'base' (noise free) values to match par en 'base' location";
				message(1, ss.str());
				Observations obs = pest_scenario.get_ctl_observations();
				oe_base.replace(base_par_idx, obs, BASE_REAL_NAME);
				ss.str("");
                if (pest_scenario.get_pestpp_options().get_save_dense())
                {
                    ss << file_manager.get_base_filename() << ".obs+noise"+dense_file_ext;
                    oe_base.to_dense(ss.str());
                }
				else if (pest_scenario.get_pestpp_options().get_save_binary())
				{
					ss << file_manager.get_base_filename() << ".obs+noise.jcb";
					oe_base.to_binary(ss.str());
				}
				else
				{
					ss << file_manager.get_base_filename() << ".obs+noise.csv";
					oe_base.to_csv(ss.str());
				}
				message(1, "re-saved obs+noise observation ensemble (obsval+noise) to ", ss.str());

				ss.str("");
				oe_base_real_names = weights.get_real_names();
                ss << "WARNING: replacing weight realization '" << oe_base_real_names[base_par_idx] << "' with 'base' weight values to match par en 'base' location";
                message(1, ss.str());
                Observations wobs;
                ObservationInfo* oi = pest_scenario.get_observation_info_ptr();
                for (auto& name : weights.get_var_names())
                {
                    wobs[name] = oi->get_weight(name);
                }
                weights.replace(base_par_idx, wobs, BASE_REAL_NAME);
                ss.str("");
                if (pest_scenario.get_pestpp_options().get_save_dense())
                {
                    ss << file_manager.get_base_filename() << ".weights.bin";
                    weights.to_binary(ss.str());
                }
                else if (pest_scenario.get_pestpp_options().get_save_binary())
                {
                    ss << file_manager.get_base_filename() << ".weights.jcb";
                    weights.to_binary(ss.str());
                }
                else
                {
                    ss << file_manager.get_base_filename() << ".weights.csv";
                    weights.to_csv(ss.str());
                }
                message(1, "re-saved weight ensemble to ", ss.str());
			}
			else
			{
				ss << "the 'base' realization was not found consistently in the restart obs en, par en, and/or the obs+noise en";
				throw_em_error(ss.str());
			}

		}
		else
		{
			ss << "the following realization names were found in the restart obs en but not in the 'obs+noise' obs en:";
			for (auto& m : missing)
				ss << m << ",";
			throw_em_error(ss.str());
		}
	}

	if (oe.shape().first != pe.shape().first)
	{
		//check if all oe names are found in par en, if so, we can reorder and proceed.  otherwise, die
		missing.clear();
		vector<string> pe_real_names = pe.get_real_names();
		for (auto& oname : oe_real_names)
		{
			if (find(pe_real_names.begin(), pe_real_names.end(), oname) == pe_real_names.end())
				missing.push_back(oname);
		}

		if (missing.size() > 0)
		{
			ss << "number of reals differ between restart obs en (" << oe.shape().first << ") and par en (" << pe.shape().first << ")";
			ss << " and realization names could not be aligned:";
			for (auto& m : missing)
				ss << m << ",";
			throw_em_error(ss.str());
		}

		message(2, "reordering pe to align with restart obs en, num reals: ", oe_real_names.size());
		try
		{
			pe.reorder(oe_real_names, vector<string>());
		}
		catch (exception& e)
		{
			ss << "error reordering pe with restart oe:" << e.what();
			throw_em_error(ss.str());
		}
		catch (...)
		{
			throw_em_error(string("error reordering pe with restart oe"));
		}

	}

	//if (oe.shape().first < oe_base.shape().first) //maybe some runs failed...
	if (oe.shape().first <= oe_base.shape().first)
	{
		message(2, "reordering oe_base to align with restart obs en,num reals:", oe_real_names.size());
		if ((oe_drawn) && (oe_base.shape().first == oe_real_names.size()))
		{
            if (pest_scenario.get_pestpp_options().get_ies_n_iter_reinflate().size() > 0)
                message(2,"reinflation is active, not reordering oe_base");
            else {
                oe_base.set_real_names(oe_real_names);
                weights.set_real_names(oe_real_names);
            }
		}
		else
		{
            if (pest_scenario.get_pestpp_options().get_ies_n_iter_reinflate().size() > 0)
                message(2,"reinflation is active, not reordering oe_base or weight ensembles");
			else {


                try {
                    oe_base.reorder(oe_real_names, vector<string>());
                }
                catch (exception &e) {
                    ss << "error reordering oe_base with restart oe:" << e.what();
                    throw_em_error(ss.str());
                }
                catch (...) {
                    throw_em_error(string("error reordering oe_base with restart oe"));
                }

                try {
                    weights.reorder(oe_real_names, vector<string>());
                }
                catch (exception &e) {
                    ss << "error reordering weights with restart oe:" << e.what();
                    throw_em_error(ss.str());
                }
                catch (...) {
                    throw_em_error(string("error reordering weights with restart oe"));
                }
            }


		}
		//if (par_restart_csv.size() > 0)
		if (pest_scenario.get_pestpp_options().get_ies_n_iter_reinflate().size() == 0)
		{
			vector<string> pe_real_names = pe.get_real_names();
			message(2, "reordering pe_base to align with restart par en,num reals:", pe_real_names.size());
			try
			{
				pe_base.reorder(pe_real_names, vector<string>());
			}
			catch (exception& e)
			{
				ss << "error reordering pe_base with restart pe:" << e.what();
				throw_em_error(ss.str());
			}
			catch (...)
			{
				throw_em_error(string("error reordering pe_base with restart pe"));
			}
		}
        else
        {
            message(2,"reinflation is active, not reordering pe_base ensemble");
            pe.set_fixed_info(pe_base.get_fixed_info());

        }


	}
	else if (oe.shape().first > oe_base.shape().first) //something is wrong
	{
		ss << "oe read from file has too many rows: " << oe.shape().first << " compared to obs+noise oe: " << oe_base.shape().first;
		throw_em_error(ss.str());
	}

}


void EnsembleMethod::initialize_parcov()
{
	stringstream ss;
	message(1, "initializing prior parameter covariance matrix");
	performance_log->log_event("initializing parcov");

	if (pest_scenario.get_pestpp_options().get_ies_use_empirical_prior())
		return;
	string how = parcov.try_from(pest_scenario, file_manager);
	message(1, "parcov loaded ", how);
	//if (parcov.e_ptr()->rows() > 0)
	parcov = parcov.get(act_par_names);

}


void EnsembleMethod::initialize_obscov()
{
	if (act_obs_names.size() == 0)
		return;
	message(1, "initializing observation noise covariance matrix");
	string obscov_filename = pest_scenario.get_pestpp_options().get_obscov_filename();

	string how = obscov.try_from(pest_scenario, file_manager, false, true);
	message(1, "obscov loaded ", how);
	if (obscov_filename.size() > 0)
	{
		vector<string> cov_names = obscov.get_col_names();
		vector<string> nz_obs_names = pest_scenario.get_ctl_ordered_nz_obs_names();
		if (cov_names.size() < nz_obs_names.size())
		{
			//this means we need to drop some nz obs
			set<string> scov(cov_names.begin(), cov_names.end());
			set<string>::iterator end = scov.end();
			vector<string> drop;
			for (auto name : nz_obs_names)
				if (scov.find(name) == end)
					drop.push_back(name);
			ofstream& frec = file_manager.rec_ofstream();
			frec << "Note: zero-weighting the following " << drop.size() << " observations that are not in the obscov:" << endl;
			int i = 0;
			for (auto n : drop)
			{
				frec << ',' << n;
				i++;
				if (i > 10)
				{
					frec << endl;
					i = 0;
				}
			}
			frec << endl;
			zero_weight_obs(drop, false, true);

		}
		message(1, "resetting weights based on obscov diagonal");
		Eigen::VectorXd diag = obscov.e_ptr()->diagonal();
		ObservationInfo* oi = pest_scenario.get_observation_info_ptr();
		for (int i = 0; i < diag.size(); i++)
		{
			oi->set_weight(cov_names[i], min(1.0 / sqrt(diag[i]), oi->get_weight(cov_names[i])));
		}
	}
	else 
		obscov = obscov.get(act_obs_names);

}
void EnsembleMethod::zero_weight_obs(vector<string>& obs_to_zero_weight, bool update_obscov, bool update_oe_base)
{
	//drop from act_obs_names
	vector<string> t;
	set<string> sdrop(obs_to_zero_weight.begin(), obs_to_zero_weight.end());
	for (auto oname : act_obs_names)
		if (sdrop.find(oname) == sdrop.end())
			t.push_back(oname);
	act_obs_names = t;
	int org_nnz_obs = pest_scenario.get_ctl_ordered_nz_obs_names().size();

	//update obscov
	if (update_obscov)
		obscov.drop(obs_to_zero_weight);

	//drop from oe_base
	if (update_oe_base)
		oe_base.drop_cols(obs_to_zero_weight);
    ObservationInfo* oi = pest_scenario.get_observation_info_ptr();
    for (auto n : obs_to_zero_weight)
    {
        oi->set_weight(n, 0.0);
    }
	//shouldn't need to update localizer since we dropping not adding
	//updating weights in control file
	if (weights.shape().first == 0)
	{
		message(2, "weight ensemble not yet initialized, skipping zero weighting");
	}
	else
	{
		map<string, int> weight_var_map = weights.get_var_map();
		for (auto n : obs_to_zero_weight)
		{
			weights.get_eigen_ptr_4_mod()->col(weight_var_map.at(n)).setZero();
		}
	}
    act_obs_names = pest_scenario.get_ctl_ordered_nz_obs_names();
	stringstream ss;
	ss << "number of non-zero weighted observations reduced from " << org_nnz_obs;
	ss << " to " << act_obs_names.size() << endl;
	message(1, ss.str());
}


Eigen::MatrixXd EnsembleMethod::get_Am(const vector<string>& real_names, const vector<string>& par_names)
{

	double scale = (1.0 / (sqrt(double(real_names.size() - 1))));
	Eigen::MatrixXd par_diff = scale * pe_base.get_eigen_anomalies(real_names, par_names, pest_scenario.get_pestpp_options().get_ies_center_on());
	par_diff.transposeInPlace();
	if (verbose_level > 1)
	{
		cout << "prior_par_diff shape: " << par_diff.rows() << ',' << par_diff.cols() << endl;
		if (verbose_level > 2)
			save_mat("prior_par_diff.dat", par_diff);
	}

	Eigen::MatrixXd ivec, upgrade_1, s, V, U, st;
	SVD_REDSVD rsvd;
	//SVD_EIGEN rsvd;
	rsvd.set_performance_log(performance_log);

	rsvd.solve_ip(par_diff, s, U, V, pest_scenario.get_svd_info().eigthresh, pest_scenario.get_svd_info().maxsing);
	par_diff.resize(0, 0);
	Eigen::MatrixXd temp = s.asDiagonal().inverse();
	Eigen::MatrixXd Am = U * temp;
	return Am;
}

void EnsembleMethod::drop_bad_reals(ParameterEnsemble& _pe, ObservationEnsemble& _oe, vector<int> subset_idxs)
{
    stringstream ss;
	//don't use this assert because _pe maybe full size, but _oe might be subset size
	bool is_subset = false;
	if (subset_idxs.size() > 0)
		is_subset = true;
	if (!is_subset)
		if (_pe.shape().first != _oe.shape().first)
			throw_em_error("EnsembleMethod::drop_bad_reals() error: _pe != _oe and not subset");

	double bad_phi = pest_scenario.get_pestpp_options().get_ies_bad_phi();
	double bad_phi_sigma = pest_scenario.get_pestpp_options().get_ies_bad_phi_sigma();
	vector<int> idxs = ph.get_idxs_greater_than(bad_phi, bad_phi_sigma, _oe, weights);
    vector<string> viol_reals = ph.get_violating_realizations(_oe,violation_obs);
    map<string,int> rmap = _oe.get_real_map();
    int idx;
    if (viol_reals.size() > 0) {
        ss.str("");
        ss << viol_reals.size() << " realizations meet 'drop_violations' conditions";
        message(2,ss.str());
        for (auto &v : viol_reals) {
            if (v == BASE_REAL_NAME) {
                message(3, "not dropping 'base' real even though it meets 'drop_violations' conditions");

            } else
            {
                idx = rmap.at(v);
                if (find(idxs.begin(), idxs.end(), idx) == idxs.end()) {
                    idxs.push_back(idx);
                }
            }
        }
    }

	if (pest_scenario.get_pestpp_options().get_ies_debug_bad_phi())
		idxs.push_back(0);

	if (idxs.size() > 0)
	{

		message(0, "dropping realizations as bad: ", idxs.size());

		vector<string> par_real_names = _pe.get_real_names(), obs_real_names = _oe.get_real_names();

		string pname;
		string oname;

		int pidx;
		vector<string> full_onames, full_pnames;
		// if a subset drop, then use the full oe index, otherwise, just use _oe index

		full_onames = oe.get_real_names();
		full_pnames = pe.get_real_names();
		vector<string> pdrop, odrop;
		for (auto i : idxs)
		{
			oname = obs_real_names[i];

			if (is_subset)
			{
				pidx = find(full_onames.begin(), full_onames.end(), oname) - full_onames.begin();
				if (find(subset_idxs.begin(), subset_idxs.end(), pidx) == subset_idxs.end())
				{
					ss.str("");
					ss << "drop_bad_reals() error: idx " << pidx << " not found in subset_idxs";
					throw_em_error(ss.str());
				}
				pname = full_pnames[pidx];
			}
			else
			{
				pidx = i;
				pname = par_real_names[pidx];
			}
			ss << pname << " : " << obs_real_names[i] << " , ";
			pdrop.push_back(pname);
			odrop.push_back(obs_real_names[i]);
		}

		string s = "dropping par:obs realizations: " + ss.str();
		message(1, s);
		try
		{
			_pe.drop_rows(pdrop);
			_oe.drop_rows(odrop);
		}
		catch (const exception& e)
		{
			stringstream ss;
			ss << "drop_bad_reals() error : " << e.what();
			throw_em_error(ss.str());
		}
		catch (...)
		{
			throw_em_error(string("drop_bad_reals() error"));
		}
	}
}

void EnsembleMethod::update_reals_by_phi(ParameterEnsemble& _pe, ObservationEnsemble& _oe, vector<int> subset_idxs)
{

	vector<string> oe_names = _oe.get_real_names();
	vector<string> pe_names = _pe.get_real_names();
	vector<string> oe_base_names = oe.get_real_names();
	vector<string> pe_base_names = pe.get_real_names();

	map<string, int> oe_name_to_idx;
	map<int, string> pe_idx_to_name;
	if (subset_idxs.size() > 0)
	{
		for (int i = 0; i < subset_idxs.size(); i++)
		{
			oe_name_to_idx[oe_base_names[subset_idxs[i]]] = i;
			pe_idx_to_name[i] = pe_base_names[subset_idxs[i]];
		}
	}
	else
	{
		for (int i = 0; i < oe_names.size(); i++)
			oe_name_to_idx[oe_names[i]] = i;
		for (int i = 0; i < pe_names.size(); i++)
			pe_idx_to_name[i] = pe_names[i];
	}
	//store map of current phi values
	ph.update(oe, pe, weights);
	L2PhiHandler::phiType pt = L2PhiHandler::phiType::COMPOSITE;
	map<string, double>* phi_map = ph.get_phi_map_ptr(pt);
	map<string, double> cur_phi_map;
	for (auto p : *phi_map)
		cur_phi_map[p.first] = p.second;

	//now get a phi map of the new phi values
	ph.update(_oe, _pe, weights);
	phi_map = ph.get_phi_map_ptr(pt);

	double acc_fac = pest_scenario.get_pestpp_options().get_ies_accept_phi_fac();
	double cur_phi, new_phi;
	string oname, pname;
	Eigen::VectorXd real;
	stringstream ss;
	for (int i = 0; i < _oe.shape().first; i++)
	{
		oname = oe_names[i];
		new_phi = phi_map->at(oname);
		cur_phi = cur_phi_map.at(oname);
		if (new_phi < cur_phi * acc_fac)
		{
			pname = pe_idx_to_name.at(oe_name_to_idx.at(oname));
			if (find(pe_names.begin(), pe_names.end(), pname) == pe_names.end())
				throw runtime_error("EnsembeMethod::update_reals_by_phi() error: pname not in pe_names: " + pname);
			ss.str("");
			ss << "updating pe:oe real =" << pname << ":" << oname << ", current phi: new phi  =" << cur_phi << ":" << new_phi;
			message(2, ss.str());

			real = _pe.get_real_vector(pname);
			pe.update_real_ip(pname, real);
			real = _oe.get_real_vector(oname);
			oe.update_real_ip(oname, real);
		}
	}
	ph.update(oe, pe, weights);

}


vector<int> EnsembleMethod::get_subset_idxs(int size, int nreal_subset)
{
	vector<int> subset_idxs;
	if ((!use_subset) || (nreal_subset >= size))
	{
		for (int i = 0; i < size; i++)
			subset_idxs.push_back(i);
		return subset_idxs;
	}
	vector<string> pe_names = pe.get_real_names();

	vector<string>::iterator bidx = find(pe_names.begin(), pe_names.end(), BASE_REAL_NAME);
	if (bidx != pe_names.end())
	{

		subset_idxs.push_back(bidx - pe_names.begin());
	}
	//int size = pe.shape().first;
	string how = pest_scenario.get_pestpp_options().get_ies_subset_how();
	if (how == "FIRST")
	{
		for (int i = 0; i < size; i++)
		{
			if (subset_idxs.size() >= nreal_subset)
				break;
			if (find(subset_idxs.begin(), subset_idxs.end(), i) != subset_idxs.end())
				continue;

			subset_idxs.push_back(i);

		}

	}
	else if (how == "LAST")
	{

		for (int i = size - 1; i >= 0; i--)
		{
			if (subset_idxs.size() >= nreal_subset)
				break;
			if (find(subset_idxs.begin(), subset_idxs.end(), i) != subset_idxs.end())
				continue;

			subset_idxs.push_back(i);

		}

	}

	else if (how == "RANDOM")
    {
		//std::uniform_int_distribution<int> uni(0, size - 1);
		int idx;
		for (int i = 0; i < 1000000000; i++)
		{
			if (subset_idxs.size() >= nreal_subset)
				break;
			//idx = uni(subset_rand_gen);
            idx = uniform_int_draws(1,0,size-1,rand_gen)[0];
            if ((idx < 0) || (idx > (size-1)))
            {
                continue;
            }
			if (find(subset_idxs.begin(), subset_idxs.end(), idx) != subset_idxs.end())
				continue;
			subset_idxs.push_back(idx);
		}
		if (subset_idxs.size() != nreal_subset)
			throw_em_error("max iterations exceeded when trying to find random subset idxs");

	}
	else if (how == "PHI_BASED")
	{
		//sidx needs to be index of realization, not realization number
		vector<pair<double, int>> phis;
		//vector<int> sidx;
		int step;
		int idx;
		L2PhiHandler::phiType pt = L2PhiHandler::phiType::COMPOSITE;
		map<string, double>* phi_map = ph.get_phi_map_ptr(pt);
		map<string, double>::iterator pi = phi_map->begin(), end = phi_map->end();

		int i = 0;
		for (; pi != end; ++pi)
		{
			phis.push_back(make_pair(pi->second, i)); //phival,idx?
			++i;
		}
		sort(phis.begin(), phis.end());

		//include idx for lowest and highest phi reals
		if (subset_idxs.size() < nreal_subset)
		{
			for (auto phi : phis)
			{
				if (find(subset_idxs.begin(), subset_idxs.end(), phi.second) == subset_idxs.end())
				{
					subset_idxs.push_back(phi.second);
					break;
				}
			}
		}
		if (subset_idxs.size() < nreal_subset)
		{
			for (int i = phis.size() - 1; i >= 0; i--)
			{
				if (find(subset_idxs.begin(), subset_idxs.end(), phis[i].second) == subset_idxs.end())
				{
					subset_idxs.push_back(phis[i].second);
					break;
				}
			}
		}


		step = (phis.size() - 1) / nreal_subset;

		for (i = 1; i < nreal_subset; ++i)
		{
			//add higher phis first
			idx = phis.size() - (i * step);
			if ((subset_idxs.size() < nreal_subset) && (find(subset_idxs.begin(), subset_idxs.end(), phis[idx].second) == subset_idxs.end()))
			{
				subset_idxs.push_back(phis[idx].second);

			}
		}
	}
	else
	{
		throw_em_error("unknown 'subset_how'");
	}

	sort(subset_idxs.begin(), subset_idxs.end());
	return subset_idxs;

}

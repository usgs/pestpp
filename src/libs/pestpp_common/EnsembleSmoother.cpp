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
#include "EnsembleMethodUtils.h"




void IterEnsembleSmoother::throw_ies_error(string message)
{
	EnsembleMethod::throw_em_error(message);
}

void IterEnsembleSmoother::sanity_checks()
{
	EnsembleMethod::sanity_checks();
}



void IterEnsembleSmoother::iterate_2_solution()
{
	stringstream ss;
	ofstream &frec = file_manager.rec_ofstream();
	
	bool accept;
	vector<int> n_iter_mean = pest_scenario.get_pestpp_options().get_ies_n_iter_mean();

    int iters_since_reinflate = 0;
    int n_iter_mean_idx = 0;
    int current_n_iter_mean = n_iter_mean[n_iter_mean_idx];
    int solution_iter = 0;
	for (int i = 0; i < pest_scenario.get_control_info().noptmax; i++)
	{
		iter++;
        solution_iter++;
		message(0, "starting solve for iteration:", iter);
		ss.str("");
		ss << "starting solve for iteration: " << iter;
		performance_log->log_event(ss.str());
		if (pest_scenario.get_pestpp_options().get_ies_use_mda())
        {
		    accept = solve_mda(false);
        }
		else {
            accept = solve_glm();
        }
		report_and_save(NetPackage::NULL_DA_CYCLE);
		ph.update(oe,pe, weights);
		last_best_mean = ph.get_representative_phi(L2PhiHandler::phiType::COMPOSITE);
		last_best_std = ph.get_std(L2PhiHandler::phiType::COMPOSITE);
		ph.report(true);
        ss.str("");
        ss << "." << iter << ".pdc.csv";
        ph.detect_simulation_data_conflict(oe,ss.str());
		ph.write(iter, run_mgr_ptr->get_total_runs());
		if (pest_scenario.get_pestpp_options().get_ies_save_rescov())
			ph.save_residual_cov(oe,iter);
		ss.str("");
		ss << file_manager.get_base_filename() << "." << iter << ".pcs.csv";
		pcs.summarize(pe,ss.str());
		if (accept)
			consec_bad_lambda_cycles = 0;
		else
			consec_bad_lambda_cycles++;

		//if ((n_iter_mean > 0) && (solution_iter % n_iter_mean == 0))
        if ((current_n_iter_mean > 0) && (iters_since_reinflate >= current_n_iter_mean))
        {
            iter++;
            reset_par_ensemble_to_prior_mean();
            iters_since_reinflate = 0;
            if (n_iter_mean.size() > n_iter_mean_idx)
            {
                n_iter_mean_idx++;
                current_n_iter_mean = n_iter_mean[n_iter_mean_idx];
            }
        }

		if (should_terminate(current_n_iter_mean))
        {
		    //if (iter > pest_scenario.get_pestpp_options().get_ies_n_iter_mean()) {
            if (iter > current_n_iter_mean) {
                break;
            }
		    else{
		        message(1,"continuing iterations to satisfy ies_n_iter_mean");
		    }
        }

	}
}


void IterEnsembleSmoother::finalize()
{

}

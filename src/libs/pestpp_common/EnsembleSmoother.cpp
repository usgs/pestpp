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
	vector<int> n_iter_reinflate = pest_scenario.get_pestpp_options().get_ies_n_iter_reinflate();
    vector<double> reinflate_factor = pest_scenario.get_pestpp_options().get_ies_reinflate_factor();
    vector<int> reinflate_num_reals = pest_scenario.get_pestpp_options().get_ies_reinflate_num_reals();

    int iters_since_reinflate = 0;
    int n_iter_reinflate_idx = 0;
    int current_n_iter_reinflate = abs(n_iter_reinflate[n_iter_reinflate_idx]);
    double current_reinflate_factor = reinflate_factor[n_iter_reinflate_idx];
    int current_num_reals = reinflate_num_reals[n_iter_reinflate_idx];
	if (reinflate_num_reals.size() > n_iter_reinflate_idx+1) {
		current_num_reals = reinflate_num_reals[n_iter_reinflate_idx+1];
	}
    int solution_iter = 0;
    int q;
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
        q = pest_utils::quit_file_found();
        if ((q == 1) || (q == 2)) {
            message(1, "'pest.stp' found, quitting");
            return;
        }
		if (accept)
			consec_bad_lambda_cycles = 0;
		else
			consec_bad_lambda_cycles++;

		//if ((n_iter_reinflate > 0) && (solution_iter % n_iter_reinflate == 0))
        iters_since_reinflate++;
        if ((current_n_iter_reinflate != 0) && (iters_since_reinflate >= current_n_iter_reinflate))
        {
            //do this so that we get a phi sequence report
            should_terminate(current_n_iter_reinflate);
            message(2,"incrementing iteration count for reinflation cycle");
            iter++;

            reinflate_par_ensemble(current_reinflate_factor,current_num_reals);
            //adjust_weights(true);
            iters_since_reinflate = 0;
            n_iter_reinflate_idx++;
            if (reinflate_factor.size() > n_iter_reinflate_idx)
            {
                current_reinflate_factor = reinflate_factor[n_iter_reinflate_idx];
            }
            if (n_iter_reinflate.size() > n_iter_reinflate_idx) {
	            current_n_iter_reinflate = abs(n_iter_reinflate[n_iter_reinflate_idx]);
            }
        	if (reinflate_num_reals.size() > n_iter_reinflate_idx+1) {
        		current_num_reals = reinflate_num_reals[n_iter_reinflate_idx+1];
        	}
            //now report again to get the new phi sequence after reinflation
            should_terminate(current_n_iter_reinflate);

        }

		else if (should_terminate(current_n_iter_reinflate))
        {
		    //if (iter > pest_scenario.get_pestpp_options().get_ies_n_iter_reinflate()) {
            if (current_n_iter_reinflate == 0) {
                break;
            }
		    else{
		        message(1,"continuing iterations because reinflation is in use");
		    }
        }
        else if (solution_iter >= pest_scenario.get_control_info().noptmax){
            message(1,"solution iterations >= noptmax, all done");
            break;
        }

	}
}


void IterEnsembleSmoother::finalize()
{

}

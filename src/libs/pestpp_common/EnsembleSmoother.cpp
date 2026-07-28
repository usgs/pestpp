/**
 * @file EnsembleSmoother.cpp
 * @brief Implementation of EnsembleSmoother.
 */
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




/**
 * @brief Throw ies error.
 *
 * @param message Description.
 */
void IterEnsembleSmoother::throw_ies_error(string message)
{
	EnsembleMethod::throw_em_error(message);
}

/**
 * @brief Sanity checks.
 */
void IterEnsembleSmoother::sanity_checks()
{
	EnsembleMethod::sanity_checks();
}



/**
 * @brief Iterate 2 solution.
 */
void IterEnsembleSmoother::iterate_2_solution()
{
	stringstream ss;
	ofstream &frec = file_manager.rec_ofstream();
	
	bool accept;
	// when/how much/how many to reinflate - the walk over the three parallel option vectors
	// lives in the schedule so an API caller driving its own loop uses the same rules
	ReinflationSchedule reinflation(pest_scenario);
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
		    accept = (solve_mda(false) != UpgradeStatus::REJECTED_RETRY);
        }
		else {
            accept = (solve_glm() != UpgradeStatus::REJECTED_RETRY);
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
        reinflation.tick();
        if (reinflation.due())
        {
            //do this so that we get a phi sequence report
            should_terminate(reinflation.get_n_iter());
            message(2,"incrementing iteration count for reinflation cycle");
            iter++;

            reinflate_par_ensemble(reinflation.get_factor(),reinflation.get_num_reals());
            //adjust_weights(true);
            reinflation.advance();
            //now report again to get the new phi sequence after reinflation
            should_terminate(reinflation.get_n_iter());

        }

		else if (should_terminate(reinflation.get_n_iter()))
        {
		    //if (iter > pest_scenario.get_pestpp_options().get_ies_n_iter_reinflate()) {
            if (!reinflation.is_active()) {
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

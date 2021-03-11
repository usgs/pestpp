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
	PestppOptions* ppo = pest_scenario.get_pestpp_options_ptr();
	vector<string> errors;
	vector<string> warnings;
	stringstream ss;
	string par_csv = ppo->get_ies_par_csv();
	string obs_csv = ppo->get_ies_obs_csv();
	string restart_obs = ppo->get_ies_obs_restart_csv();
	string restart_par = ppo->get_ies_par_restart_csv();


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
		warnings.push_back("prior information equations not supported in pestpp-ies, ignoring...");
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
	if (ppo->get_ies_bad_phi() <= 0.0)
		errors.push_back("ies_bad_phi <= 0.0, really?");
	if ((ppo->get_ies_num_reals() < error_min_reals) && (par_csv.size() == 0))
	{
		ss.str("");
		ss << "ies_num_reals < " << error_min_reals << ", this is redic, increaing to " << warn_min_reals;
		warnings.push_back(ss.str());
		ppo->set_ies_num_reals(warn_min_reals);
	}
	if ((ppo->get_ies_num_reals() < warn_min_reals) && (par_csv.size() == 0))
	{
		ss.str("");
		ss << "ies_num_reals < " << warn_min_reals << ", this is prob too few";
		warnings.push_back(ss.str());
	}
	if (ppo->get_ies_reg_factor() < 0.0)
		errors.push_back("ies_reg_factor < 0.0 - WRONG!");
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
		throw_ies_error(string("sanity_check() found some problems - please review rec file"));
	}
}



void IterEnsembleSmoother::iterate_2_solution()
{
	stringstream ss;
	ofstream &frec = file_manager.rec_ofstream();
	
	bool accept;
	for (int i = 0; i < pest_scenario.get_control_info().noptmax; i++)
	{
		iter++;
		message(0, "starting solve for iteration:", iter);
		ss.str("");
		ss << "starting solve for iteration: " << iter;
		performance_log->log_event(ss.str());
		accept = solve_glm();
		report_and_save(NetPackage::NULL_DA_CYCLE);
		ph.update(oe,pe);
		last_best_mean = ph.get_mean(L2PhiHandler::phiType::COMPOSITE);
		last_best_std = ph.get_std(L2PhiHandler::phiType::COMPOSITE);
		ph.report(true);
		ph.write(iter, run_mgr_ptr->get_total_runs());
		if (pest_scenario.get_pestpp_options().get_ies_save_rescov())
			ph.save_residual_cov(oe,iter);
		ss.str("");
		ss << file_manager.get_base_filename() << "." << iter << ".pcs.csv";
		pcs.summarize(pe,iter,ss.str());
			
			
		if (accept)
			consec_bad_lambda_cycles = 0;
		else
			consec_bad_lambda_cycles++;

		if (should_terminate())
			break;
	}
}


void IterEnsembleSmoother::finalize()
{

}

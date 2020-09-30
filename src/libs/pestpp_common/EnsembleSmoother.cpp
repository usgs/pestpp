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



IterEnsembleSmoother::IterEnsembleSmoother(Pest &_pest_scenario, FileManager &_file_manager,
	OutputFileWriter &_output_file_writer, PerformanceLog *_performance_log,
	RunManagerAbstract* _run_mgr_ptr) : pest_scenario(_pest_scenario), file_manager(_file_manager),
	output_file_writer(_output_file_writer), performance_log(_performance_log),
	run_mgr_ptr(_run_mgr_ptr)
{
	rand_gen = std::mt19937(pest_scenario.get_pestpp_options().get_random_seed());
	subset_rand_gen = std::mt19937(pest_scenario.get_pestpp_options().get_random_seed());
	pe.set_pest_scenario(&pest_scenario);
	oe.set_pest_scenario(&pest_scenario);
	pe.set_rand_gen(&rand_gen);
	oe.set_rand_gen(&rand_gen);
	localizer.set_pest_scenario(&pest_scenario);
	
}

void IterEnsembleSmoother::throw_ies_error(string message)
{
	performance_log->log_event("IterEnsembleSmoother error: " + message);
	cout << endl << "   ************   " << endl << "    IterEnsembleSmoother error: " << message << endl << endl;
	file_manager.rec_ofstream() << endl << "   ************   " << endl << "    IterEnsembleSmoother error: " << message << endl << endl;
	file_manager.close_file("rec");
	performance_log->~PerformanceLog();
	throw runtime_error("IterEnsembleSmoother error: " + message);
}

bool IterEnsembleSmoother::initialize_pe(Covariance &cov)
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
		pe.draw(num_reals, draw_par,cov, performance_log, pest_scenario.get_pestpp_options().get_ies_verbose_level(), file_manager.rec_ofstream());
		// stringstream ss;
		// ss << file_manager.get_base_filename() << ".0.par.csv";
		// message(1, "saving initial parameter ensemble to ", ss.str());
		// pe.to_csv(ss.str());
		drawn = true;
	}
	else
	{
		string par_ext = pest_utils::lower_cp(par_csv).substr(par_csv.size() - 3, par_csv.size());
		performance_log->log_event("processing par csv " + par_csv);
		if (par_ext.compare("csv") == 0)
		{
			message(1, "loading par ensemble from csv file", par_csv);
			try
			{
				pe.from_csv(par_csv);
			}
			catch (const exception &e)
			{
				ss << "error processing par csv: " << e.what();
				throw_ies_error(ss.str());
			}
			catch (...)
			{
				throw_ies_error(string("error processing par csv"));
			}
		}
		else if ((par_ext.compare("jcb") == 0) || (par_ext.compare("jco") == 0))
		{
			message(1, "loading par ensemble from binary file", par_csv);
			try
			{
				pe.from_binary(par_csv);
			}
			catch (const exception &e)
			{
				ss << "error processing par jcb: " << e.what();
				throw_ies_error(ss.str());
			}
			catch (...)
			{
				throw_ies_error(string("error processing par jcb"));
			}
		}
		else
		{
			ss << "unrecognized par csv extension " << par_ext << ", looking for csv, jcb, or jco";
			throw_ies_error(ss.str());
		}

		pe.transform_ip(ParameterEnsemble::transStatus::NUM);
		
		if (pp_args.find("IES_NUM_REALS") != pp_args.end())
		{
			int num_reals = pest_scenario.get_pestpp_options().get_ies_num_reals();
			/*if (pest_scenario.get_pestpp_options().get_ies_include_base())
			{
				message(1, "Note: increasing num_reals by 1 to account for 'base' realization in existing par ensemble");
				num_reals++;
			}*/
			if (num_reals < pe.shape().first)
			{
				message(1,"ies_num_reals arg passed, truncated parameter ensemble to ",num_reals);
				vector<string> keep_names,real_names=pe.get_real_names();
				for (int i=0;i<num_reals;i++)
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

				if (pest_scenario.get_pestpp_options().get_ies_save_binary())
				{
					string filename = file_manager.get_base_filename() + ".prior.jcb";
					message(1, "saving emprirical parameter covariance matrix to binary file: ", filename);
					parcov.to_binary(filename);
				}
				else
				{
					string filename = file_manager.get_base_filename() + ".prior.cov";
					message(1, "saving emprirical parameter covariance matrix to ASCII file: ", filename);
					parcov.to_ascii(filename);
				}
			}
		}
		if (pest_scenario.get_pestpp_options().get_ies_enforce_bounds())
		{
			if (pest_scenario.get_pestpp_options().get_ies_obs_restart_csv().size() > 0)
				message(1, "Warning: even though ies_enforce_bounds is true, a restart obs en was passed, so bounds will not be enforced on the initial par en");
			else
				pe.enforce_limits(performance_log, pest_scenario.get_pestpp_options().get_ies_enforce_chglim());
		}

	}
	return drawn;

}

void IterEnsembleSmoother::add_bases()
{
	//check that 'base' isn't already in ensemble
	vector<string> rnames = pe.get_real_names();
	bool inpar = false;
	if (find(rnames.begin(), rnames.end(), base_name) != rnames.end())
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
		pe.drop_rows(drop);
		pe.append(base_name, pars);
	}

	//check that 'base' isn't already in ensemble
	rnames = oe.get_real_names();
	if (find(rnames.begin(), rnames.end(), base_name) != rnames.end())
	{
		message(1, "'base' realization already in observation ensemble, ignoring '++ies_include_base'");
	}
	else
	{
		Observations obs = pest_scenario.get_ctl_observations();
		if (inpar)
		{
			vector<string> prnames = pe.get_real_names();

			int idx = find(prnames.begin(), prnames.end(), base_name) - prnames.begin();
			//cout << idx << "," << rnames.size() << endl;
			string oreal = rnames[idx];
			stringstream ss;
			ss << "warning: 'base' realization in par ensenmble but not in obs ensemble," << endl;
			ss << "         replacing obs realization '" << oreal << "' with 'base'";
			string mess = ss.str();
			message(1, mess);
			vector<string> drop;
			drop.push_back(oreal);
			oe.drop_rows(drop);
			oe.append(base_name, obs);
			//rnames.insert(rnames.begin() + idx, string(base_name));
			rnames[idx] = base_name;
			oe.reorder(rnames, vector<string>());
		}
		else
		{
			message(1, "adding 'base' observation values to ensemble");
			vector<int> drop{ oe.shape().first - 1 };
			oe.drop_rows(drop);
			oe.append(base_name, obs);
		}
	}
}

bool IterEnsembleSmoother::initialize_oe(Covariance &cov)
{
	stringstream ss;
	int num_reals = pe.shape().first;
	string obs_csv = pest_scenario.get_pestpp_options().get_ies_obs_csv();
	bool drawn = false;
	if (obs_csv.size() == 0)
	{
		if (pest_scenario.get_pestpp_options().get_ies_no_noise())
		{
			message(1, "initializing no-noise observation ensemble of : ", num_reals);
			oe.initialize_without_noise(num_reals);
		
		}
		else
		{
			message(1, "drawing observation noise realizations: ", num_reals);
			oe.draw(num_reals, cov, performance_log, pest_scenario.get_pestpp_options().get_ies_verbose_level(), file_manager.rec_ofstream());
			
		}
		drawn = true;
	}
	else
	{
		string obs_ext = pest_utils::lower_cp(obs_csv).substr(obs_csv.size() - 3, obs_csv.size());
		performance_log->log_event("processing obs csv " + obs_csv);
		if (obs_ext.compare("csv") == 0)
		{
			message(1, "loading obs ensemble from csv file", obs_csv);
			try
			{
				oe.from_csv(obs_csv);
			}
			catch (const exception &e)
			{
				ss << "error processing obs csv: " << e.what();
				throw_ies_error(ss.str());
			}
			catch (...)
			{
				throw_ies_error(string("error processing obs csv"));
			}
		}
		else if ((obs_ext.compare("jcb") == 0) || (obs_ext.compare("jco") == 0))
		{
			message(1, "loading obs ensemble from binary file", obs_csv);
			try
			{
				oe.from_binary(obs_csv);
			}
			catch (const exception &e)
			{
				stringstream ss;
				ss << "error processing obs binary file: " << e.what();
				throw_ies_error(ss.str());
			}
			catch (...)
			{
				throw_ies_error(string("error processing obs binary file"));
			}
		}
		else
		{
			ss << "unrecognized obs ensemble extension " << obs_ext << ", looking for csv, jcb, or jco";
			throw_ies_error(ss.str());
		}
		if (pp_args.find("IES_NUM_REALS") != pp_args.end())
		{
			int num_reals = pest_scenario.get_pestpp_options().get_ies_num_reals();
			/*if (pest_scenario.get_pestpp_options().get_ies_include_base())
			{
				message(1, "Note: increasing num_reals by 1 to account for 'base' realization in existing obs ensemble");
				num_reals++;
			}*/
			if (num_reals < oe.shape().first)
			{
				message(1,"ies_num_reals arg passed, truncated observation ensemble to ",num_reals);
				vector<string> keep_names,real_names=oe.get_real_names();
				for (int i=0;i<num_reals;i++)
				{
					keep_names.push_back(real_names[i]);
				}
				oe.keep_rows(keep_names);
			}
		}
	}
	return drawn;

}

//template<typename T, typename A>
//void IterEnsembleSmoother::message(int level, char* _message, vector<T, A> _extras)
//{
//	string s(_message);
//	message(level, s, _extras);
//}

//void IterEnsembleSmoother::message(int level, char* _message)
//{
//	string s(_message);
//	message(level, s);
//}

//template<typename T>
//void IterEnsembleSmoother::message(int level, char* _message, T extra)
//{
//	string s(_message);
//	message(level, s, extra);
//
//}

template<typename T, typename A>
void IterEnsembleSmoother::message(int level, const string &_message, vector<T, A> _extras, bool echo)
{
	stringstream ss;
	if (level == 0)
		ss << endl << "  ---  ";
	else if (level == 1)
		ss << "...";
	ss << _message;
	if (_extras.size() > 0)
	{

		for (auto &e : _extras)
			ss << e << " , ";

	}
	if (level == 0)
		ss << "  ---  ";
	if ((echo) && ((verbose_level >= 2) || (level < 2)))
		cout << ss.str() << endl;
	file_manager.rec_ofstream() <<ss.str() << endl;
	performance_log->log_event(_message);

}

void IterEnsembleSmoother::message(int level, const string &_message)
{
	message(level, _message, vector<string>());
}

template<typename T>
void IterEnsembleSmoother::message(int level, const string &_message, T extra)
{
	stringstream ss;
	ss << _message << " " << extra;
	string s = ss.str();
	message(level, s);
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
	
	/*if (!ppo->get_ies_use_prior_scaling())
	{
		warnings.push_back("not using prior scaling - this is really a dev option, you should always use prior scaling...");
	}*/

	/*if ((ppo->get_ies_num_threads() > 0) && (!use_localizer))
	{
		warnings.push_back("'ies_num_threads > 0 but no localization, resetting 'ies_num_threads' = 0");
		ppo->set_ies_num_threads(-1);
		num_threads = -1;
	}
*/


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
	//cout << endl << endl;
}

void IterEnsembleSmoother::initialize_restart()
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
		catch (const exception &e)
		{
			ss << "error processing restart obs csv: " << e.what();
			throw_ies_error(ss.str());
		}
		catch (...)
		{
			throw_ies_error(string("error processing restart obs csv"));
		}
	}
	else if ((obs_ext.compare("jcb") == 0) || (obs_ext.compare("jco") == 0))
	{
		message(1, "loading restart obs ensemble from binary file", obs_restart_csv);
		try
		{
			oe.from_binary(obs_restart_csv);
		}
		catch (const exception &e)
		{
			ss << "error processing restart obs binary file: " << e.what();
			throw_ies_error(ss.str());
		}
		catch (...)
		{
			throw_ies_error(string("error processing restart obs binary file"));
		}
	}
	else
	{
		ss << "unrecognized restart obs ensemble extension " << obs_ext << ", looking for csv, jcb, or jco";
		throw_ies_error(ss.str());
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
			catch (const exception &e)
			{
				ss << "error processing restart par csv: " << e.what();
				throw_ies_error(ss.str());
			}
			catch (...)
			{
				throw_ies_error(string("error processing restart par csv"));
			}
		}
		else if ((par_ext.compare("jcb") == 0) || (par_ext.compare("jco") == 0))
		{
			message(1, "loading restart par ensemble from binary file", par_restart_csv);
			try
			{
				pe.from_binary(par_restart_csv);
			}
			catch (const exception &e)
			{
				ss << "error processing restart par binary file: " << e.what();
				throw_ies_error(ss.str());
			}
			catch (...)
			{
				throw_ies_error(string("error processing restart par binary file"));
			}
		}
		else
		{
			ss << "unrecognized restart par ensemble extension " << par_ext << ", looking for csv, jcb, or jco";
			throw_ies_error(ss.str());
		}
		if (pe.shape().first != oe.shape().first)
		{
			ss.str("");
			ss << "restart par en has " << pe.shape().first << " realizations but restart obs en has " << oe.shape().first;
			throw_ies_error(ss.str());
		}

		//check that restart pe is in sync with pe_base
		vector<string> pe_real_names = pe.get_real_names(), pe_base_real_names = pe_base.get_real_names();
		vector<string>::const_iterator start, end;
		vector<string> missing;
		start = pe_base_real_names.begin();
		end = pe_base_real_names.end();
		for (auto &rname : pe_real_names)
			if (find(start, end, rname) == end)
				missing.push_back(rname);
		if (missing.size() > 0)
		{
			ss << "the following realization names were found in the restart par en but not in the 'base' par en:";
			for (auto &m : missing)
				ss << m << ",";
			throw_ies_error(ss.str());
		}
		pe.transform_ip(ParameterEnsemble::transStatus::NUM);
	}

	if (pp_args.find("IES_NUM_REALS") != pp_args.end())
	{
		int num_reals = pest_scenario.get_pestpp_options().get_ies_num_reals();
		/*if (pest_scenario.get_pestpp_options().get_ies_include_base())
		{
			message(1, "Note: increasing num_reals by 1 to account for 'base' realization in existing obs restart ensemble");
			num_reals++;
		}*/
		if (num_reals < oe.shape().first)
		{
			message(1, "ies_num_reals arg passed, truncated restart obs ensemble to ", num_reals);
			vector<string> keep_names, real_names = oe.get_real_names();
			for (int i = 0; i<num_reals; i++)
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
	for (auto &rname : oe_real_names)
		if (find(start, end, rname) == end)
			missing.push_back(rname);
	if (missing.size() > 0)
	{
		//the special case where the base real is what is missing...
		if ((missing.size() == 1) && (missing[0] == "BASE"))
		{
			//check that the base real is in the par en - restart_par_en should be accounted for by now
			int base_par_idx = -1;
			vector<string> pe_real_names = pe.get_real_names(), pe_base_real_names = pe_base.get_real_names();
			for (int i = 0; i < pe_base.shape().first; i++)
			{
				if (pe_base_real_names[i] == "BASE")
				{
					base_par_idx = i;
					break;
				}
			}
			if (base_par_idx != -1)
			{
				ss.str("");
				ss << "WARNING: replacing base obs en realization '" << oe_base_real_names[base_par_idx] << "' with 'base' (noise free) values to match par en 'base' location";
				message(2, ss.str());
				Observations obs = pest_scenario.get_ctl_observations();
				oe_base.replace(base_par_idx, obs, "BASE");
			}
			else
			{
				ss << "the 'base' realization was not found in the restart obs en and also not found in the par en";
				throw_ies_error(ss.str());
			}

		}
		else
		{
			ss << "the following realization names were found in the restart obs en but not in the 'base' obs en:";
			for (auto &m : missing)
				ss << m << ",";
			throw_ies_error(ss.str());
		}

	}

	if (oe.shape().first != pe.shape().first)
	{
		//check if all oe names are found in par en, if so, we can reorder and proceed.  otherwise, die
		missing.clear();
		vector<string> pe_real_names = pe.get_real_names();
		for (auto &oname : oe_real_names)
		{
			if (find(pe_real_names.begin(), pe_real_names.end(), oname) == pe_real_names.end())
				missing.push_back(oname);
		}

		if (missing.size() > 0)
		{
			ss << "number of reals differ between restart obs en (" << oe.shape().first << ") and par en (" << pe.shape().first << ")";
			ss << " and realization names could not be aligned:";
			for (auto &m : missing)
				ss << m << ",";
			throw_ies_error(ss.str());
		}

		message(2, "reordering pe to align with restart obs en, num reals: ", oe_real_names.size());
		try
		{
			pe.reorder(oe_real_names, vector<string>());
		}
		catch (exception &e)
		{
			ss << "error reordering pe with restart oe:" << e.what();
			throw_ies_error(ss.str());
		}
		catch (...)
		{
			throw_ies_error(string("error reordering pe with restart oe"));
		}

	}

	//if (oe.shape().first < oe_base.shape().first) //maybe some runs failed...
	if (oe.shape().first <= oe_base.shape().first)
	{
		//find which realizations are missing and reorder oe_base, pe and pe_base

		//message(1, "shape mismatch detected with restart obs ensemble...checking for compatibility");
		
		/*vector<string> pe_real_names;
		start = oe_base_real_names.begin();
		end = oe_base_real_names.end();
		vector<string>::const_iterator it;
		int iit;
		for (int i = 0; i < oe.shape().first; i++)
		{
			it = find(start, end, oe_real_names[i]);
			if (it != end)
			{
				iit = it - start;
				pe_real_names.push_back(pe_org_real_names[iit]);
			}
		}*/
		message(2, "reordering oe_base to align with restart obs en,num reals:", oe_real_names.size());
		if ((oe_drawn) && (oe_base.shape().first == oe_real_names.size()))
		{
			oe_base.set_real_names(oe_real_names);
		}
		else
		{
			try
			{
				oe_base.reorder(oe_real_names, vector<string>());
			}
			catch (exception& e)
			{
				ss << "error reordering oe_base with restart oe:" << e.what();
				throw_ies_error(ss.str());
			}
			catch (...)
			{
				throw_ies_error(string("error reordering oe_base with restart oe"));
			}
		}
		//if (par_restart_csv.size() > 0)
		if (true)
		{
			vector<string> pe_real_names = pe.get_real_names();
			message(2, "reordering pe_base to align with restart par en,num reals:", pe_real_names.size());
			try
			{
				pe_base.reorder(pe_real_names, vector<string>());
			}
			catch (exception &e)
			{
				ss << "error reordering pe_base with restart pe:" << e.what();
				throw_ies_error(ss.str());
			}
			catch (...)
			{
				throw_ies_error(string("error reordering pe_base with restart pe"));
			}
		}


		/*try
		{
			pe.reorder(pe_real_names, vector<string>());
		}
		catch (exception &e)
		{
			ss << "error reordering pe with restart oe:" << e.what();
			throw_ies_error(ss.str());
		}
		catch (...)
		{
			throw_ies_error(string("error reordering pe with restart oe"));
		}*/
	}
	else if (oe.shape().first > oe_base.shape().first) //something is wrong
	{
		ss << "restart oe has too many rows: " << oe.shape().first << " compared to oe_base: " << oe_base.shape().first;
		throw_ies_error(ss.str());
	}
}


void IterEnsembleSmoother::initialize_parcov()
{
	stringstream ss;
	performance_log->log_event("initializing parcov");

	if (pest_scenario.get_pestpp_options().get_ies_use_empirical_prior())
		return;
	string how = parcov.try_from(pest_scenario, file_manager);
	message(1, "parcov loaded ", how);
	//if (parcov.e_ptr()->rows() > 0)
	parcov = parcov.get(act_par_names);

}


void IterEnsembleSmoother::initialize_obscov()
{
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
	{
		obscov = obscov.get(act_obs_names);
	}
	
}


void IterEnsembleSmoother::initialize()
{	
	message(0, "initializing");
	pp_args = pest_scenario.get_pestpp_options().get_passed_args();

	act_obs_names = pest_scenario.get_ctl_ordered_nz_obs_names();
	act_par_names = pest_scenario.get_ctl_ordered_adj_par_names();

	stringstream ss;

	if (pest_scenario.get_control_info().noptmax == 0)
	{
		message(0, "'noptmax'=0, running control file parameter values and quitting");
		
		Parameters pars = pest_scenario.get_ctl_parameters();
		ParamTransformSeq pts = pe.get_par_transform();

		ParameterEnsemble _pe(&pest_scenario, &rand_gen);
		_pe.reserve(vector<string>(), pest_scenario.get_ctl_ordered_par_names());
		_pe.set_trans_status(ParameterEnsemble::transStatus::CTL);
		_pe.append("BASE", pars);
		string par_csv = file_manager.get_base_filename() + ".par.csv";
		//message(1, "saving parameter values to ", par_csv);
		//_pe.to_csv(par_csv);
		pe_base = _pe;
		pe_base.reorder(vector<string>(), act_par_names);
		ObservationEnsemble _oe(&pest_scenario, &rand_gen);
		_oe.reserve(vector<string>(), pest_scenario.get_ctl_ordered_obs_names());
		_oe.append("BASE", pest_scenario.get_ctl_observations());
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
		message(1, "running control file parameter values");

		vector<int> failed_idxs = run_ensemble(_pe, _oe);
		if (failed_idxs.size() != 0)
		{
			message(0, "control file parameter value run failed...bummer");
			throw_ies_error("control file parameter value run failed");
		}
		string obs_csv = file_manager.get_base_filename() + ".obs.csv";
		message(1, "saving results from control file parameter value run to ", obs_csv);
		_oe.to_csv(obs_csv);

		ph.update(_oe, _pe);
		message(0, "control file parameter phi report:");
		ph.report(true);
		ph.write(0, 1);
		save_base_real_par_rei(pest_scenario, _pe, _oe, output_file_writer, file_manager, -1);			
		return;
	}

	//set some defaults
	PestppOptions *ppo = pest_scenario.get_pestpp_options_ptr();

	verbose_level = pest_scenario.get_pestpp_options_ptr()->get_ies_verbose_level();
	if (pest_scenario.get_n_adj_par() >= 1e6)
	{
		message(0, "You are a god among mere mortals!");
	}

	message(1, "using REDSVD for truncated svd solve");
	message(1, "maxsing:", pest_scenario.get_svd_info().maxsing);
	message(1, "eigthresh: ", pest_scenario.get_svd_info().eigthresh);

	message(1, "initializing localizer");
	use_localizer = localizer.initialize(performance_log);
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
			message(1, "localizing by obseravtions");
		else
			message(1, "localizing by parameters");
	}
	iter = 0;
	//ofstream &frec = file_manager.rec_ofstream();
	last_best_mean = 1.0E+30;
	last_best_std = 1.0e+30;
	lambda_max = 1.0E+30;
	lambda_min = 1.0E-30;
	warn_min_reals = 10;
	error_min_reals = 2;
	consec_bad_lambda_cycles = 0;

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

	sanity_checks();

	bool echo = false;
	if (verbose_level > 1)
		echo = true;

	initialize_parcov();
	initialize_obscov();

	subset_size = pest_scenario.get_pestpp_options().get_ies_subset_size();
	reg_factor = pest_scenario.get_pestpp_options().get_ies_reg_factor();
	message(1,"using reg_factor: ", reg_factor);
	double bad_phi = pest_scenario.get_pestpp_options().get_ies_bad_phi();
	if (bad_phi < 1.0e+30)
		message(1, "using bad_phi: ", bad_phi);

	int num_reals = pest_scenario.get_pestpp_options().get_ies_num_reals();

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

	oe_drawn = initialize_oe(obscov);
	string center_on = ppo->get_ies_center_on();
	
	try
	{
		pe.check_for_dups();
	}
	catch (const exception &e)
	{
		string message = e.what();
		throw_ies_error("error in parameter ensemble: " + message);
	}

	try
	{
		oe.check_for_dups();
	}
	catch (const exception &e)
	{
		string message = e.what();
		throw_ies_error("error in observation ensemble: " + message);
	}

	if (pe.shape().first != oe.shape().first)
	{
		//the special case where par en < obs en and all par reals are found in obs en...

		if (pe.shape().first < oe.shape().first)
		{
			vector<string> oe_names = oe.get_real_names();
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
				ss << "par en has " << pe.shape().first << " realizations, compared to " << oe.shape().first << " obs realizations";
				message(1, ss.str());
				message(1," the realization names are compatible");
				message(1, "re-indexing obs en to align with par en...");

				oe.reorder(pe.get_real_names(), vector<string>());
			}
			else
			{
				ss.str("");
				ss << "the following par en real names were not found in the obs en: ";
				for (auto m : missing)
				{
					ss << m << ",";
				}
				throw_ies_error(ss.str());
			}
		}
		else
		{
			ss.str("");
			ss << "parameter ensemble rows (" << pe.shape().first << ") not equal to observation ensemble rows (" << oe.shape().first << ")";
			throw_ies_error(ss.str());
		}
	}


	//need this here for Am calcs...
	//message(1, "transforming parameter ensemble to numeric");
	pe.transform_ip(ParameterEnsemble::transStatus::NUM);

	if (pest_scenario.get_pestpp_options().get_ies_include_base())
		if (pp_args.find("IES_RESTART_OBS_EN") != pp_args.end())
		{
			message(1, "Warning: even though `ies_include_base` is true, you passed a restart obs en, not adding 'base' realization...");
		}
		else
			add_bases();


	//now we check to see if we need to try to align the par and obs en
	//this would only be needed if either of these were not drawn
	if (!pe_drawn || !oe_drawn)
	{
		bool aligned = pe.try_align_other_rows(performance_log, oe);
		if (aligned)
		{
			message(2, "observation ensemble reordered to align rows with parameter ensemble");
		}
	}

	//just check to see if common real names are found but are not in the same location
	map<string, int> pe_map = pe.get_real_map(), oe_map = oe.get_real_map();
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
		message(1, "WARNING: common realization names shared between the parameter and observation ensembles but they are not in the same row locations, see .rec file for listing");
		ofstream& frec = file_manager.rec_ofstream();
		frec << endl <<  "WARNING: the following " << misaligned.size() << " realization names are shared between the parameter and observation ensembles but they are not in the same row locations:" << endl;
		for (auto ma : misaligned)
			frec << ma << endl;
	}



	message(2, "checking for denormal values in pe");
	pe.check_for_normal("initial transformed parameter ensemble");
	ss.str("");
	if (pest_scenario.get_pestpp_options().get_ies_save_binary())
	{
		ss << file_manager.get_base_filename() << ".0.par.jcb";
		pe.to_binary(ss.str());
	}
	else
	{
		ss << file_manager.get_base_filename() << ".0.par.csv";
		pe.to_csv(ss.str());
	}
	message(1, "saved initial parameter ensemble to ", ss.str());
	message(2, "checking for denormal values in base oe");
	oe.check_for_normal("obs+noise observation ensemble");
	ss.str("");
	if (pest_scenario.get_pestpp_options().get_ies_save_binary())
	{
		ss << file_manager.get_base_filename() << ".obs+noise.jcb";
		oe.to_binary(ss.str());
	}
	else
	{
		ss << file_manager.get_base_filename() << ".obs+noise.csv";
		oe.to_csv(ss.str());
	}
	message(1, "saved obs+noise observation ensemble (obsval+noise) to ", ss.str());


	if (center_on.size() > 0)
	{
		ss.str("");
		ss << "centering on realization: '" << center_on << "' ";
		message(1, ss.str());
		vector<string> names = pe.get_real_names();
		if (find(names.begin(), names.end(), center_on) == names.end())
			throw_ies_error("'ies_center_on' realization not found in par en: " + center_on);
		names = oe.get_real_names();
		if (find(names.begin(), names.end(), center_on) == names.end())
			throw_ies_error("'ies_center_on' realization not found in obs en: " + center_on);
	}
	else
		message(1, "centering on ensemble mean vector");

	if (pest_scenario.get_control_info().noptmax == -2)
	{
		message(0, "'noptmax'=-2, running mean parameter ensemble values and quitting");
		message(1, "calculating mean parameter values");
		Parameters pars;
		vector<double> mv = pe.get_mean_stl_var_vector();
		pars.update(pe.get_var_names(), pe.get_mean_stl_var_vector());
		ParamTransformSeq pts = pe.get_par_transform();

		ParameterEnsemble _pe(&pest_scenario, &rand_gen);
		_pe.reserve(vector<string>(), pe.get_var_names());
		_pe.set_trans_status(pe.get_trans_status());
		_pe.append("mean", pars);
		string par_csv = file_manager.get_base_filename() + ".mean.par.csv";
		message(1, "saving mean parameter values to ", par_csv);
		_pe.to_csv(par_csv);
		pe_base = _pe;
		pe_base.reorder(vector<string>(), act_par_names);
		ObservationEnsemble _oe(&pest_scenario, &rand_gen);
		_oe.reserve(vector<string>(), oe.get_var_names());
		_oe.append("mean", pest_scenario.get_ctl_observations());
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
		message(1, "running mean parameter values");

		vector<int> failed_idxs = run_ensemble(_pe, _oe);
		if (failed_idxs.size() != 0)
		{
			message(0, "mean parmeter value run failed...bummer");
			return;
		}
		string obs_csv = file_manager.get_base_filename() + ".mean.obs.csv";
		message(1, "saving results from mean parameter value run to ", obs_csv);
		_oe.to_csv(obs_csv);

		ph.update(_oe, _pe);
		message(0, "mean parameter phi report:");
		ph.report(true);
		ph.write(0, 1);

		return;
	}

	if (subset_size > pe.shape().first)
	{
		use_subset = false;
	}
	else
	{
		message(1, "using subset in lambda testing, number of realizations used in subset testing: ", subset_size);
		string how = pest_scenario.get_pestpp_options().get_ies_subset_how();
		message(1, "subset how: ", how);
		use_subset = true;
	}

	oe_org_real_names = oe.get_real_names();
	pe_org_real_names = pe.get_real_names();
	string obs_restart_csv = pest_scenario.get_pestpp_options().get_ies_obs_restart_csv();
	string par_restart_csv = pest_scenario.get_pestpp_options().get_ies_par_restart_csv();

	oe_base = oe; //copy
	//reorder this for later...
	oe_base.reorder(vector<string>(), act_obs_names);

	pe_base = pe; //copy
	//reorder this for later
	pe_base.reorder(vector<string>(), act_par_names);


	//the hard way to restart
	if (obs_restart_csv.size() > 0)
		initialize_restart();
	
	//no restart
	else
	{
		performance_log->log_event("running initial ensemble");
		message(1, "running initial ensemble of size", oe.shape().first);
		vector<int> failed = run_ensemble(pe, oe);
		if (pe.shape().first == 0)
			throw_ies_error("all realizations failed during initial evaluation");

		pe.transform_ip(ParameterEnsemble::transStatus::NUM);
	}
	
	ss.str("");
	if (pest_scenario.get_pestpp_options().get_ies_save_binary())
	{
		ss << file_manager.get_base_filename() << ".0.obs.jcb";
		oe.to_binary(ss.str());
	}
	else
	{
		ss << file_manager.get_base_filename() << ".0.obs.csv";
		oe.to_csv(ss.str());
	}
	message(1, "saved initial obs ensemble to", ss.str());

	//save the 0th iter par and rei and well as the untagged par and rei
	save_base_real_par_rei(pest_scenario, pe, oe, output_file_writer, file_manager, iter);
	save_base_real_par_rei(pest_scenario, pe, oe, output_file_writer, file_manager, -1);


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

	ph.update(oe, pe);
	message(0, "pre-drop initial phi summary");
	ph.report(true);
	drop_bad_phi(pe, oe);
	if (oe.shape().first == 0)
	{
		throw_ies_error(string("all realizations dropped as 'bad'"));
	}
	if (oe.shape().first <= error_min_reals)
	{
		message(0, "too few active realizations:", oe.shape().first);
		message(1, "need at least ", error_min_reals);
		throw_ies_error(string("too few active realizations, cannot continue"));
	}
	if (oe.shape().first < warn_min_reals)
	{
		ss.str("");
		ss << "WARNING: less than " << warn_min_reals << " active realizations...might not be enough";
		string s = ss.str();
		message(0, s);
	}

	pcs = ParChangeSummarizer(&pe_base, &file_manager,&output_file_writer);
	vector<string> in_conflict = detect_prior_data_conflict();
	if (in_conflict.size() > 0)
	{
		ss.str("");
		ss << "WARNING: " << in_conflict.size() << " non-zero weighted observations are in conflict";
		ss << " with the prior simulated ensemble." << endl;
		message(0, ss.str());

		cout << "...see rec file or " << file_manager.get_base_filename() << ".pdc.csv" << "for listing of conflicted observations" << endl << endl;
		ofstream& frec = file_manager.rec_ofstream();
		frec << endl << "...conflicted observations: " << endl;
		for (auto oname : in_conflict)
		{
			frec << oname << endl;
		}
		
		if (!ppo->get_ies_drop_conflicts())
		{
			ss.str("");
			ss << "  WARNING: Prior-data conflict detected.  Continuing with IES parameter" << endl;
			ss << "           adjustment will likely result in parameter and forecast bias." << endl;
			ss << "           Consider using 'ies_drop_conflicts' as a quick fix.";
			message(1, ss.str());
		}
		else
		{

			//check that all obs are in conflict
			message(1, "dropping conflicted observations");
			if (in_conflict.size() == oe.shape().second)
			{
				throw_ies_error("all non-zero weighted observations in conflict state, cannot continue");
			}
			zero_weight_obs(in_conflict);
			if (ppo->get_ies_localizer().size() > 0)
			{
				message(1, "updating localizer");
					use_localizer = localizer.initialize(performance_log,true);
			}

		}
		string filename = file_manager.get_base_filename() + ".adjusted.obs_data.csv";
		ofstream f_obs(filename);
		if (f_obs.bad())
			throw_ies_error("error opening: " + filename);
		output_file_writer.scenario_obs_csv(f_obs);
		f_obs.close();
		message(1, "updated observation data information written to file ", filename);
	}
	else if (ppo->get_obscov_filename().size() > 0)
	{
		string filename = file_manager.get_base_filename() + ".adjusted.obs_data.csv";
		ofstream f_obs(filename);
		if (f_obs.bad())
			throw_ies_error("error opening: " + filename);
		output_file_writer.scenario_obs_csv(f_obs);
		f_obs.close();
		message(1, "updated observation data information written to file ", filename);
	}
	performance_log->log_event("calc initial phi");
	ph.update(oe, pe);
	message(0, "initial phi summary");
	ph.report(true);
	ph.write(0, run_mgr_ptr->get_total_runs());
	if (ppo->get_ies_save_rescov())
		ph.save_residual_cov(oe, 0);
	best_mean_phis.push_back(ph.get_mean(L2PhiHandler::phiType::COMPOSITE));
	if (!pest_scenario.get_pestpp_options().get_ies_use_approx())
	{
		message(1, "using full (MAP) update solution");

	}

	if (ph.get_mean(L2PhiHandler::phiType::ACTUAL) < 1.0e-10)
	{	
		throw_ies_error("initial actual phi mean too low, something is wrong...or you have the perfect model that already fits the data shockingly well");	
	}
	if (ph.get_std(L2PhiHandler::phiType::ACTUAL) < 1.0e-10)
	{
		throw_ies_error("initial actual phi stdev too low, something is wrong...");
	}


	last_best_mean = ph.get_mean(L2PhiHandler::phiType::COMPOSITE);
	last_best_std = ph.get_std(L2PhiHandler::phiType::COMPOSITE);
	last_best_lam = pest_scenario.get_pestpp_options().get_ies_init_lam();
	if (last_best_mean < 1.0e-10)
	{
		throw_ies_error("initial composite phi mean too low, something is wrong...");
	}
	if (last_best_std < 1.0e-10)
	{
		throw_ies_error("initial composite phi stdev too low, something is wrong...");
	}
	if (last_best_lam <= 0.0)
	{
		//double x = last_best_mean / (2.0 * double(oe.shape().second));
		double x = last_best_mean / (2.0 * double(pest_scenario.get_ctl_ordered_nz_obs_names().size()));
		last_best_lam = pow(10.0, (floor(log10(x))));
		if (last_best_lam < 1.0e-10)
		{
			message(1, "initial lambda estimation from phi failed, using 10,000");
			last_best_lam = 10000;
		}
	}
	message(1, "current lambda:", last_best_lam);
	message(0, "initialization complete");
}

void IterEnsembleSmoother::zero_weight_obs(vector<string>& obs_to_zero_weight, bool update_obscov, bool update_oe_base)
{
	//drop from act_obs_names
	vector<string> t;
	set<string> sdrop(obs_to_zero_weight.begin(), obs_to_zero_weight.end());
	for (auto oname : act_obs_names)
		if (sdrop.find(oname) == sdrop.end())
			t.push_back(oname);
	act_obs_names = t;

	//update obscov
	if (update_obscov)
		obscov.drop(obs_to_zero_weight);

	//drop from oe_base
	if (update_oe_base)
		oe_base.drop_cols(obs_to_zero_weight);
	
	//shouldnt need to update localizer since we dropping not adding
	//updating weights in control file

	ObservationInfo* oi = pest_scenario.get_observation_info_ptr();
	int org_nnz_obs = pest_scenario.get_ctl_ordered_nz_obs_names().size();
	for (auto n : obs_to_zero_weight)
	{
		oi->set_weight(n, 0.0);
	}

	stringstream ss;
	ss << "number of non-zero weighted observations reduced from " << org_nnz_obs;
	ss << " to " << pest_scenario.get_ctl_ordered_nz_obs_names().size() << endl;
	message(1, ss.str());
}

vector<string> IterEnsembleSmoother::detect_prior_data_conflict()
{
	message(1, "checking for prior-data conflict...");
	//for now, just really simple metric - checking for overlap
	// write out conflicted obs and some related info to a csv file
	ofstream pdccsv(file_manager.get_base_filename() + ".pdc.csv");
	
	vector<string> in_conflict;
	double smin, smax, omin, omax,smin_stat, smax_stat, omin_stat,omax_stat;
	map<string, int> smap, omap;
	vector<string> snames = oe.get_var_names();
	vector<string> onames = oe_base.get_var_names();
	vector<string> temp = ph.get_lt_obs_names();
	set<string> ineq(temp.begin(), temp.end());
	set<string>::iterator end = ineq.end();
	temp = ph.get_gt_obs_names();
	ineq.insert(temp.begin(), temp.end());
	temp.resize(0);
	
	for (int i = 0; i < snames.size(); i++)
	{
		smap[snames[i]] = i;
	}
	for (int i = 0; i < onames.size(); i++)
	{
		omap[onames[i]] = i;
	}
	int sidx, oidx;
	bool use_stat_dist = true;
	if (pest_scenario.get_pestpp_options().get_ies_pdc_sigma_distance() <= 0.0)
		use_stat_dist = false;
	
	double smn, sstd, omn, ostd,dist;
	double sd = pest_scenario.get_pestpp_options().get_ies_pdc_sigma_distance();
	int oe_nr = oe.shape().first;
	int oe_base_nr = oe_base.shape().first;
	Eigen::VectorXd t;
	pdccsv << "name,obs_mean,obs_std,obs_min,obs_max,obs_stat_min,obs_stat_max,sim_mean,sim_std,sim_min,sim_max,sim_stat_min,sim_stat_max,distance" << endl;
	for (auto oname : pest_scenario.get_ctl_ordered_nz_obs_names())
	{
		if (ineq.find(oname) != end)
			continue;
		sidx = smap[oname];
		oidx = omap[oname];
		smin = oe.get_eigen_ptr()->col(sidx).minCoeff();
		omin = oe_base.get_eigen_ptr()->col(oidx).minCoeff();
		smax = oe.get_eigen_ptr()->col(sidx).maxCoeff();
		omax = oe_base.get_eigen_ptr()->col(oidx).maxCoeff();
		t = oe.get_eigen_ptr()->col(sidx);
		smn = t.mean();
		sstd = std::sqrt((t.array() - smn).square().sum() / (oe_nr-1));
		smin_stat = smn - (sd * sstd);
		smax_stat = smn + (sd * sstd);
		t = oe_base.get_eigen_ptr()->col(oidx);
		omn = t.mean();
		ostd = std::sqrt((t.array() - omn).square().sum() / (oe_base_nr-1));
		omin_stat = omn - (sd * ostd);
		omax_stat = omn + (sd * ostd);	

		if (use_stat_dist)
		{
			if ((smin_stat > omax_stat) || (smax_stat < omin_stat))
			{
				in_conflict.push_back(oname);
				dist = max((smin_stat - omax_stat), (omin_stat - smax_stat));
				pdccsv << oname << "," << omn << "," << ostd << "," << omin << "," << omax << "," << omin_stat << "," << omax_stat;
				pdccsv << "," << smn << "," << sstd << "," << smin << "," << smax << "," << smin_stat << "," << smax_stat << "," << dist << endl;
			}
		}
		else
		{
			if ((smin > omax) || (smax < omin))
			{
				in_conflict.push_back(oname);
				dist = max((smin - omax), (omin - smax));
				pdccsv << oname << "," << omn << "," << ostd << "," << omin << "," << omax << "," << omin_stat << "," << omax_stat;
				pdccsv << "," << smn << "," << sstd << "," << smin << "," << smax << "," << smin_stat << "," << smax_stat << "," << dist << endl;
			}
		}
	}
	
	return in_conflict;
}

Eigen::MatrixXd IterEnsembleSmoother::get_Am(const vector<string> &real_names, const vector<string> &par_names)
{

	double scale = (1.0 / (sqrt(double(real_names.size() - 1))));
	Eigen::MatrixXd par_diff = scale * pe_base.get_eigen_anomalies(real_names,par_names,pest_scenario.get_pestpp_options().get_ies_center_on());
	par_diff.transposeInPlace();
	if (verbose_level > 1)
	{
		cout << "prior_par_diff: " << par_diff.rows() << ',' << par_diff.cols() << endl;
		if (verbose_level > 2)
			save_mat("prior_par_diff.dat", par_diff);
	}

	Eigen::MatrixXd ivec, upgrade_1, s, V, U, st;
	SVD_REDSVD rsvd;
	//SVD_EIGEN rsvd;
	rsvd.set_performance_log(performance_log);

	rsvd.solve_ip(par_diff, s, U, V, pest_scenario.get_svd_info().eigthresh, pest_scenario.get_svd_info().maxsing);
	par_diff.resize(0, 0);
	Eigen::MatrixXd temp = s.asDiagonal();
	Eigen::MatrixXd Am = U * temp;
	return Am;
}

void IterEnsembleSmoother::drop_bad_phi(ParameterEnsemble &_pe, ObservationEnsemble &_oe, bool is_subset)
{
	//don't use this assert because _pe maybe full size, but _oe might be subset size
	if (!is_subset)
		if (_pe.shape().first != _oe.shape().first)
			throw_ies_error("IterEnsembleSmoother::drop_bad_phi() error: _pe != _oe and not subset");
		
	double bad_phi = pest_scenario.get_pestpp_options().get_ies_bad_phi();
	double bad_phi_sigma = pest_scenario.get_pestpp_options().get_ies_bad_phi_sigma();
	vector<int> idxs = ph.get_idxs_greater_than(bad_phi,bad_phi_sigma, _oe);

	if (pest_scenario.get_pestpp_options().get_ies_debug_bad_phi())
		idxs.push_back(0);
	
	if (idxs.size() > 0)
	{

		message(0, "dropping realizations as bad: ", idxs.size());

		vector<string> par_real_names = _pe.get_real_names(), obs_real_names = _oe.get_real_names();
		stringstream ss;
		string pname;
		string oname;

		int pidx;
		vector<string> full_onames, full_pnames;
		// if a subset drop, then use the full oe index, otherwise, just use _oe index
		/*if (_oe.shape().first != _pe.shape().first)
		{
			full_onames = oe.get_real_names();
		}
		else
		{
			full_onames = _oe.get_real_names();
		}*/
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
					ss << "drop_bad_phi() error: idx " << pidx << " not found in subset_idxs";
					throw_ies_error(ss.str());
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

		string s = "dropping par:obs realizations: "+ ss.str();
		message(1, s);
		try
		{
			_pe.drop_rows(pdrop);
			_oe.drop_rows(odrop);
		}
		catch (const exception &e)
		{
			stringstream ss;
			ss << "drop_bad_phi() error : " << e.what();
			throw_ies_error(ss.str());
		}
		catch (...)
		{
			throw_ies_error(string("drop_bad_phi() error"));
		}
	}
}

void IterEnsembleSmoother::save_mat(string prefix, Eigen::MatrixXd &mat)
{
	stringstream ss;
	ss << iter << '.' << prefix;
	try
	{
		ofstream &f = file_manager.open_ofile_ext(ss.str());
		f << mat << endl;
		f.close();
		file_manager.close_file(ss.str());
	}
	catch (...)
	{
		message(1, "error saving matrix", ss.str());
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
		accept = solve_new();
		report_and_save();
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

bool IterEnsembleSmoother::should_terminate()
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
	if (best_mean_phis.size() > 0)
	{
		vector<double>::iterator idx = min_element(best_mean_phis.begin(), best_mean_phis.end());
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

	for (auto &phi : best_mean_phis)
	{
		ratio = (phi - best_phi_yet) / phi;
		if (ratio <= phiredstp)
			count++;
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

	if ((nphinored_sat) || (phiredstp_sat) || (consec_sat))
	{
		message(1, "phi-based termination criteria satisfied, all done");
		return true;
	}
	return false;
}




LocalUpgradeThread::LocalUpgradeThread(PerformanceLog *_performance_log, unordered_map<string, Eigen::VectorXd> &_par_resid_map, unordered_map<string, Eigen::VectorXd> &_par_diff_map,
	unordered_map<string, Eigen::VectorXd> &_obs_resid_map, unordered_map<string, Eigen::VectorXd> &_obs_diff_map,
	Localizer &_localizer, unordered_map<string, double> &_parcov_inv_map, unordered_map<string, double> &_weight_map,
	ParameterEnsemble &_pe_upgrade, unordered_map<string,pair<vector<string>,vector<string>>> &_cases,
	unordered_map<string, Eigen::VectorXd> &_Am_map, Localizer::How &_how): par_resid_map(_par_resid_map),
	par_diff_map(_par_diff_map), obs_resid_map(_obs_resid_map),obs_diff_map(_obs_diff_map),localizer(_localizer),
	pe_upgrade(_pe_upgrade),cases(_cases), parcov_inv_map(_parcov_inv_map), weight_map(_weight_map), Am_map(_Am_map)
{
	performance_log = _performance_log;
	how = _how;
	parcov_inv_map = _parcov_inv_map;
	weight_map = _weight_map;
	count = 0;

	for (auto &c : cases)
	{
		keys.push_back(c.first);
	}
	//sort(keys.begin(), keys.end());
	total = keys.size();
	//random_shuffle(keys.begin(), keys.end());

}


void LocalUpgradeThread::work(int thread_id, int iter, double cur_lam)
{
	class local_utils
	{
	public:
		static Eigen::DiagonalMatrix<double, Eigen::Dynamic> get_matrix_from_map(vector<string> &names, unordered_map<string, double> &dmap)
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
		static Eigen::MatrixXd get_matrix_from_map(int num_reals, vector<string> &names, unordered_map<string, Eigen::VectorXd> &emap)
		{
			Eigen::MatrixXd mat(num_reals, names.size());
			mat.setZero();

			for (int j = 0; j < names.size(); j++)
			{
				mat.col(j) = emap[names[j]];
			}

			return mat;
		}
		static void save_mat(int verbose_level, int tid, int iter, int t_count, string prefix, Eigen::MatrixXd &mat)
		{
			if (verbose_level < 2)
				return;
			
			if (verbose_level < 3)
				return;
			//cout << "thread: " << tid << ", " << t_count << ", " << prefix << " rows:cols" << mat.rows() << ":" << mat.cols() << endl;
			stringstream ss;
			
			ss << "thread_"<< tid << ".count_ " << t_count << ".iter_" << iter << "." << prefix << ".dat";
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
	int maxsing, num_reals, verbose_level, pcount=0,t_count;
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
	Eigen::MatrixXd obs_resid, obs_diff, loc;
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> weights, parcov_inv;
	vector<string> par_names, obs_names;
	while (true)
	{
		unique_lock<mutex> next_guard(next_lock, defer_lock);
		par_names.clear();
		obs_names.clear();
		use_localizer = false;
		//the end condition
		//the end condition
		while (true)
		{	
			if (next_guard.try_lock())
			{
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
				string k = keys[count];
				pair<vector<string>, vector<string>> p = cases.at(k);
				par_names = p.second;
				obs_names = p.first;
				if (localizer.get_use())
				{
					if ((loc_by_obs) && (par_names.size() == 1) && (k == par_names[0]))
						use_localizer = true;
					else if ((!loc_by_obs) && (obs_names.size() == 1) && (k == obs_names[0]))
					{
						use_localizer = true;
						//loc_by_obs = false;
					}
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

		par_resid.resize(0, 0);
		par_diff.resize(0, 0);
		obs_resid.resize(0, 0);
		obs_diff.resize(0, 0);
		loc.resize(0, 0);
		Am.resize(0, 0);
		weights.resize(0);
		parcov_inv.resize(0);
		Am.resize(0, 0);

		unique_lock<mutex> obs_diff_guard(obs_diff_lock, defer_lock);
		unique_lock<mutex> obs_resid_guard(obs_resid_lock, defer_lock);
		unique_lock<mutex> par_diff_guard(par_diff_lock, defer_lock);
		unique_lock<mutex> par_resid_guard(par_resid_lock, defer_lock);
		unique_lock<mutex> loc_guard(loc_lock, defer_lock);
		unique_lock<mutex> weight_guard(weight_lock, defer_lock);
		unique_lock<mutex> parcov_guard(parcov_lock, defer_lock);
		unique_lock<mutex> am_guard(am_lock, defer_lock);

		while (true)
		{
			if (((use_approx) || (par_resid.rows() > 0)) &&
				(weights.size() > 0) &&
				(parcov_inv.size() > 0) &&
				(par_diff.rows() > 0) &&
				(obs_resid.rows() > 0) &&
				(obs_diff.rows() > 0) &&
				((!use_localizer) || (loc.rows() > 0)) && 
				((use_approx) || (Am.rows() > 0)))
				break;
			if ((use_localizer) && (loc.rows() == 0) && (loc_guard.try_lock()))
			{
				if (loc_by_obs)
					loc = localizer.get_localizing_par_hadamard_matrix(num_reals, obs_names[0], par_names);
				else
					loc = localizer.get_localizing_obs_hadamard_matrix(num_reals, par_names[0], obs_names);
				loc_guard.unlock();
			}
			if ((obs_diff.rows() == 0) && (obs_diff_guard.try_lock()))
			{
				//piggy back here for thread safety
				//if (pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_svd_pack() == PestppOptions::SVD_PACK::PROPACK)
				//	use_propack = true;
				obs_diff = local_utils::get_matrix_from_map(num_reals, obs_names, obs_diff_map);
				obs_diff_guard.unlock();
			}
			if ((obs_resid.rows() == 0) && (obs_resid_guard.try_lock()))
			{
				obs_resid = local_utils::get_matrix_from_map(num_reals, obs_names, obs_resid_map);
				obs_resid_guard.unlock();
			}
			if ((par_diff.rows() == 0) && (par_diff_guard.try_lock()))
			{
				par_diff = local_utils::get_matrix_from_map(num_reals, par_names, par_diff_map);
				par_diff_guard.unlock();
			}
			if ((par_resid.rows() == 0) && (par_resid_guard.try_lock()))
			{
				par_resid = local_utils::get_matrix_from_map(num_reals, par_names, par_resid_map);
				par_resid_guard.unlock();
			}
			if ((weights.rows() == 0) && (weight_guard.try_lock()))
			{
				weights = local_utils::get_matrix_from_map(obs_names, weight_map);
				weight_guard.unlock();
			}
			if ((parcov_inv.rows() == 0) && (parcov_guard.try_lock()))
			{
				parcov_inv = local_utils::get_matrix_from_map(par_names, parcov_inv_map);
				parcov_guard.unlock();
			}
			if ((!use_approx) && (Am.rows() == 0) && (am_guard.try_lock()))
			{
				//Am = local_utils::get_matrix_from_map(num_reals, par_names, Am_map).transpose();
				int am_cols = Am_map[par_names[0]].size();
				Am.resize(par_names.size(),am_cols);
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

		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "obs_resid", obs_resid);
		Eigen::MatrixXd scaled_residual = weights * obs_resid;
		
		

		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "par_resid", par_resid);
		Eigen::MatrixXd scaled_par_resid;
		if ((!use_approx) && (iter > 1))
		{
			if (use_prior_scaling)
			{
				scaled_par_resid = parcov_inv * par_resid;
			}
			else
			{
				scaled_par_resid = par_resid;
			}
		}

		stringstream ss;

		double scale = (1.0 / (sqrt(double(num_reals - 1))));
		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "obs_diff", obs_diff);

		if (use_localizer)
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "loc", loc);
		if (use_localizer)
		{
			if (loc_by_obs)
				par_diff = par_diff.cwiseProduct(loc);
			else	
				obs_diff = obs_diff.cwiseProduct(loc);

		}
		
		obs_diff = scale * (weights * obs_diff);
		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "par_diff", par_diff);
		if (use_prior_scaling)
			par_diff = scale * parcov_inv * par_diff;
		else
			par_diff = scale * par_diff;


		//performance_log->log_event("SVD of obs diff");
		Eigen::MatrixXd ivec, upgrade_1, s, V, Ut;
		
		
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
		Eigen::MatrixXd X1 = Ut * scaled_residual;
		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X1", X1);
		Eigen::MatrixXd X2 = ivec * X1;
		X1.resize(0, 0);

		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X2", X2);
		Eigen::MatrixXd X3 = V * s.asDiagonal() * X2;
		X2.resize(0, 0);

		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X3", X3);
		if (use_prior_scaling)
		{
			upgrade_1 = -1.0 * parcov_inv * par_diff * X3;
		}
		else
		{
			upgrade_1 = -1.0 * par_diff * X3;
		}
		upgrade_1.transposeInPlace();
		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "upgrade_1",upgrade_1);
		X3.resize(0, 0);


		Eigen::MatrixXd upgrade_2;
		if ((!use_approx) && (iter > 1))
		{
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "Am",Am);
			Eigen::MatrixXd x4 = Am.transpose() * scaled_par_resid;
			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X4", x4);
			
			par_resid.resize(0, 0);

			Eigen::MatrixXd x5 = Am * x4;
			x4.resize(0, 0);
			Am.resize(0, 0);

			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X5", x5);
			Eigen::MatrixXd x6 = par_diff.transpose() * x5;
			x5.resize(0, 0);

			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X6", x6);
			Eigen::MatrixXd x7 = V * ivec *V.transpose() * x6;
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
		
		unique_lock<mutex> put_guard(put_lock, defer_lock);
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



void upgrade_thread_function(int id, int iter,double cur_lam, LocalUpgradeThread &worker, exception_ptr &eptr)
{
	try
	{
		worker.work(id, iter, cur_lam);
	}
	catch (...)
	{
		eptr = current_exception();
	}
	
	return;
}


ParameterEnsemble IterEnsembleSmoother::calc_localized_upgrade_threaded(double cur_lam, unordered_map<string, pair<vector<string>, vector<string>>> &loc_map)
{
	stringstream ss;
	
	ObservationEnsemble oe_upgrade(oe.get_pest_scenario_ptr(), &rand_gen, oe.get_eigen(vector<string>(), act_obs_names, false), oe.get_real_names(), act_obs_names);
	ParameterEnsemble pe_upgrade(pe.get_pest_scenario_ptr(),&rand_gen, pe.get_eigen(vector<string>(), act_par_names, false), pe.get_real_names(), act_par_names);
	
	//this copy of the localizer map will be consumed by the worker threads
	//unordered_map<string, pair<vector<string>, vector<string>>> loc_map;
	if (use_localizer)
	{
		
		//loc_map = localizer.get_localizer_map(iter, oe, pe, performance_log);
		//localizer.report(file_manager.rec_ofstream());
	}
	else
	{
		pair<vector<string>, vector<string>> p(act_obs_names, act_par_names);
		loc_map["all"] = p;
	}
	
	//prep the fast look par cov info
	message(2, "preparing fast-look containers for threaded localization solve");
	unordered_map<string, double> parcov_inv_map;
	parcov_inv_map.reserve(pe_upgrade.shape().second);
	Eigen::VectorXd parcov_inv;// = parcov.get(par_names).inv().e_ptr()->toDense().cwiseSqrt().asDiagonal();
	if (!parcov.isdiagonal())
	{
		//parcov_inv = parcov.inv().get_matrix().diagonal().cwiseSqrt();
		parcov_inv = parcov.get_matrix().diagonal();
		
	}
	else
	{
		message(2, "extracting diagonal from prior parameter covariance matrix");
		Covariance parcov_diag;
		parcov_diag.from_diagonal(parcov);
		//parcov_inv = parcov_diag.inv().get_matrix().diagonal().cwiseSqrt();
		parcov_inv = parcov_diag.get_matrix().diagonal();
		
	}
	parcov_inv = parcov_inv.cwiseSqrt().cwiseInverse();

	vector<string> par_names = pe_upgrade.get_var_names();
	for (int i = 0; i < parcov_inv.size(); i++)
		parcov_inv_map[par_names[i]] = parcov_inv[i];


	//prep the fast lookup weights info
	unordered_map<string, double> weight_map;
	weight_map.reserve(oe_upgrade.shape().second);
	vector<string> obs_names = oe_upgrade.get_var_names();
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
	unordered_map<string, Eigen::VectorXd> par_resid_map, obs_resid_map, Am_map;
	unordered_map<string, Eigen::VectorXd> par_diff_map, obs_diff_map;
	par_resid_map.reserve(par_names.size());
	par_diff_map.reserve(par_names.size());
	Am_map.reserve(par_names.size());
	obs_resid_map.reserve(obs_names.size());
	obs_diff_map.reserve(obs_names.size());
	
	//check for the 'center_on' real - it may have been dropped...
	string center_on = pest_scenario.get_pestpp_options().get_ies_center_on();
	if (center_on.size() > 0)
	{
		vector<string> real_names = oe_upgrade.get_real_names();
		if (find(real_names.begin(), real_names.end(), center_on) == real_names.end())
		{
			message(0, "Warning: 'ies_center_on' real not found in obs en, reverting to mean...");
			center_on = "";
		}
		real_names = pe_upgrade.get_real_names();
		if ((center_on.size() > 0) && find(real_names.begin(), real_names.end(), center_on) == real_names.end())
		{
			message(0, "Warning: 'ies_center_on' real not found in par en, reverting to mean...");
			center_on = "";
		}
	}

	Eigen::MatrixXd mat = ph.get_obs_resid_subset(oe_upgrade);
	for (int i = 0; i < obs_names.size(); i++)
	{
		obs_resid_map[obs_names[i]] = mat.col(i);
	}
	mat = oe_upgrade.get_eigen_anomalies(center_on);
	for (int i = 0; i < obs_names.size(); i++)
	{
		obs_diff_map[obs_names[i]] = mat.col(i);
	}
	mat = ph.get_par_resid_subset(pe_upgrade);
	for (int i = 0; i <par_names.size(); i++)
	{
		par_resid_map[par_names[i]] = mat.col(i);
	}	
	mat = pe_upgrade.get_eigen_anomalies(center_on);
	for (int i = 0; i <par_names.size(); i++)
	{
		par_diff_map[par_names[i]] = mat.col(i);
	}
	if (!pest_scenario.get_pestpp_options().get_ies_use_approx())
	{
		mat = get_Am(pe_upgrade.get_real_names(), pe_upgrade.get_var_names());
		for (int i = 0; i < par_names.size(); i++)
		{
			Am_map[par_names[i]] = mat.row(i);
		}
	}
	mat.resize(0, 0);
	// clear the upgrade ensemble
	pe_upgrade.set_zeros();
	Localizer::How _how = localizer.get_how();
	LocalUpgradeThread worker(performance_log, par_resid_map, par_diff_map, obs_resid_map, obs_diff_map,
		localizer, parcov_inv_map, weight_map, pe_upgrade, loc_map, Am_map, _how);

	if ((num_threads < 1) || (loc_map.size() == 1))
	//if (num_threads < 1)
	{
		worker.work(0, iter, cur_lam);
	}
	else
	{
		Eigen::setNbThreads(1);
		vector<thread> threads;
		vector<exception_ptr> exception_ptrs;
		message(2, "launching threads");
		
		for (int i = 0; i <num_threads; i++)
		{
			exception_ptrs.push_back(exception_ptr());
		}

		for (int i = 0; i < num_threads; i++)
		{
			//threads.push_back(thread(&LocalUpgradeThread::work, &worker, i, iter, cur_lam));
			
			threads.push_back(thread(upgrade_thread_function, i, iter, cur_lam, std::ref(worker),std::ref( exception_ptrs[i])));
			
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
		message(2, "threaded localized upgrade calculation done");
	}
	
	return pe_upgrade;
}


//ParameterEnsemble IterEnsembleSmoother::calc_localized_upgrade(double cur_lam)
//{
//	stringstream ss;
//	int i = 0;
//	int lsize = localizer.get_localizer_map().size();
//	ParameterEnsemble pe_upgrade = pe.zeros_like();
//	for (auto local_pair : localizer.get_localizer_map())
//	{
//		ss.str("");
//		ss << "localized upgrade part " << i + 1 << " of " << lsize;
//		message(2, ss.str());
//		ParameterEnsemble pe_local;
//		if (localizer.get_how() == Localizer::How::PARAMETERS)
//		{
//			pe_local = calc_upgrade(local_pair.second.first, local_pair.second.second, cur_lam, pe.shape().first, local_pair.first);
//		}
//		else
//		{
//			pe_local = calc_upgrade(local_pair.second.first, local_pair.second.second, cur_lam, pe.shape().first);
//		}
//
//		pe_upgrade.add_2_cols_ip(pe_local);
//		i++;
//	}
//	return pe_upgrade;
//
//}

void IterEnsembleSmoother::update_reals_by_phi(ParameterEnsemble &_pe, ObservationEnsemble &_oe)
{

	vector<string> oe_names = _oe.get_real_names();
	vector<string> pe_names = _pe.get_real_names();
	vector<string> oe_base_names = oe.get_real_names();
	vector<string> pe_base_names = pe.get_real_names();

	//if (pe_names.size() != oe_base_names.size())
	//	throw runtime_error("IterEnsembleSmoother::update_reals_by_phi() error: pe_names != oe_base_names");
	map<string, int> oe_name_to_idx;
	map<int,string> pe_idx_to_name;

	for (int i = 0; i < oe_base_names.size(); i++)
		oe_name_to_idx[oe_base_names[i]] = i;
	
	for (int i = 0; i < pe_base_names.size(); i++)
		pe_idx_to_name[i] = pe_base_names[i];
	//store map of current phi values
	ph.update(oe, pe);
	L2PhiHandler::phiType pt = L2PhiHandler::phiType::COMPOSITE;
	map<string, double> *phi_map = ph.get_phi_map(pt);
	map<string, double> cur_phi_map;
	for (auto p : *phi_map)
		cur_phi_map[p.first] = p.second;

	//now get a phi map of the new phi values
	ph.update(_oe, _pe);
	phi_map = ph.get_phi_map(pt);
	
	double acc_fac = pest_scenario.get_pestpp_options().get_ies_accept_phi_fac();
	double cur_phi, new_phi;
	string oname, pname;
	Eigen::VectorXd real;
	stringstream ss;
	for (int i=0;i<_oe.shape().first;i++)
	{
		oname = oe_names[i];
		new_phi = phi_map->at(oname);
		cur_phi = cur_phi_map.at(oname);
		if (new_phi < cur_phi * acc_fac)
		{
			//pname = pe_names[i];
			//pname = pe_names[oe_name_to_idx[oname]];
			pname = pe_idx_to_name[oe_name_to_idx[oname]];
			if (find(pe_names.begin(), pe_names.end(), pname) == pe_names.end())
				throw runtime_error("IterEnsembleSmoother::update_reals_by_phi() error: pname not in pe_names: " + pname);
			ss.str("");
			ss << "updating pe:oe real =" << pname << ":" << oname << ", current phi: new phi  =" << cur_phi << ":" << new_phi;
			message(3, ss.str());
			
			real = _pe.get_real_vector(pname);
			pe.update_real_ip(pname, real);
			real = _oe.get_real_vector(oname);
			oe.update_real_ip(oname, real);
		}
	}
	ph.update(oe, pe);

}

bool IterEnsembleSmoother::solve_new()
{
	stringstream ss;
	ofstream &frec = file_manager.rec_ofstream();
	if (pe.shape().first <= error_min_reals)
	{
		message(0, "too few active realizations:", oe.shape().first);
		message(1, "need at least ", error_min_reals);
		throw_ies_error(string("too few active realizations, cannot continue"));
	}
	if (pe.shape().first < warn_min_reals)
	{
		ss.str("");
		ss << "WARNING: less than " << warn_min_reals << " active realizations...might not be enough";
		string s = ss.str();
		message(1, s);
	}

	if ((use_subset) && (subset_size > pe.shape().first))
	{
		ss.str("");
		ss << "++ies_subset size (" << subset_size << ") greater than ensemble size (" << pe.shape().first << ")";
		frec << "  ---  " << ss.str() << endl;
		cout << "  ---  " << ss.str() << endl;
		frec << "  ...reducing ++ies_subset_size to " << pe.shape().first << endl;
		cout << "  ...reducing ++ies_subset_size to " << pe.shape().first << endl;
		subset_size = pe.shape().first;
	}

	pe.transform_ip(ParameterEnsemble::transStatus::NUM);

	vector<ParameterEnsemble> pe_lams;
	vector<double> lam_vals, scale_vals;
	//update all the fast-lookup structures
	oe.update_var_map();
	pe.update_var_map();
	parcov.update_sets();
	obscov.update_sets();
	unordered_map<string, pair<vector<string>, vector<string>>> loc_map;
	if (use_localizer)
	{

		loc_map = localizer.get_localizer_map(iter, oe, pe, performance_log);
		//localizer.report(file_manager.rec_ofstream());
	}
	else
	{
		pair<vector<string>, vector<string>> p(act_obs_names, act_par_names);
		loc_map["all"] = p;
	}


	for (auto &lam_mult : lam_mults)
	{
		ss.str("");
		double cur_lam = last_best_lam * lam_mult;
		ss << "starting calcs for lambda" << cur_lam;
		message(1, "starting lambda calcs for lambda", cur_lam);
		message(2, "see .log file for more details");
		ParameterEnsemble pe_upgrade;
		
		pe_upgrade = calc_localized_upgrade_threaded(cur_lam, loc_map);

		for (auto sf : pest_scenario.get_pestpp_options().get_lambda_scale_vec())
		{

			ParameterEnsemble pe_lam_scale = pe;
			pe_lam_scale.set_eigen(*pe_lam_scale.get_eigen_ptr() + (*pe_upgrade.get_eigen_ptr() * sf));
			if (pest_scenario.get_pestpp_options().get_ies_enforce_bounds())
			{
				pe_lam_scale.enforce_limits(performance_log, pest_scenario.get_pestpp_options().get_ies_enforce_chglim());
			}

			pe_lams.push_back(pe_lam_scale);
			lam_vals.push_back(cur_lam);
			scale_vals.push_back(sf);
			if (!pest_scenario.get_pestpp_options().get_ies_save_lambda_en())
				continue;
			ss.str("");
			ss << file_manager.get_base_filename() << "." << iter << "." << cur_lam << ".lambda." << sf << ".scale.par";

			if (pest_scenario.get_pestpp_options().get_ies_save_binary())
			{
				ss << ".jcb";
				pe_lam_scale.to_binary(ss.str());
			}
			else
			{
				ss << ".csv";
				pe_lam_scale.to_csv(ss.str());
			}
			frec << "lambda, scale value " << cur_lam << ',' << sf << " pars saved to " << ss.str() << endl;

		}

		ss.str("");
		message(1, "finished calcs for lambda:", cur_lam);

	}
	
	if (pest_scenario.get_pestpp_options().get_ies_debug_upgrade_only())
	{
		message(0, "ies_debug_upgrade_only is true, exiting");
		throw_ies_error("ies_debug_upgrade_only is true, exiting");
	}

	vector<map<int, int>> real_run_ids_lams;
	int best_idx = -1;
	double best_mean = 1.0e+30, best_std = 1.0e+30;
	double mean, std;

	message(0, "running lambda ensembles");
	vector<ObservationEnsemble> oe_lams = run_lambda_ensembles(pe_lams, lam_vals, scale_vals);

	message(0, "evaluting lambda ensembles");
	message(1, "last mean: ", last_best_mean);
	message(1, "last stdev: ", last_best_std);

	ObservationEnsemble oe_lam_best;
	bool echo = false;
	if (verbose_level > 1)
		echo = true;
	for (int i = 0; i<pe_lams.size(); i++)
	{
		if (oe_lams[i].shape().first == 0)
			continue;
		vector<double> vals({ lam_vals[i],scale_vals[i] });
		if (pest_scenario.get_pestpp_options().get_ies_save_lambda_en())
			
		{
			ss.str("");
			ss << file_manager.get_base_filename() << "." << iter << "." << lam_vals[i] << ".lambda." << scale_vals[i] << ".scale.obs";

			if (pest_scenario.get_pestpp_options().get_ies_save_binary())
			{
				ss << ".jcb";
				oe_lams[i].to_binary(ss.str());
			}
			else
			{
				ss << ".csv";
				oe_lams[i].to_csv(ss.str());
			}
			frec << "lambda, scale value " << lam_vals[i] << ',' << scale_vals[i] << " pars saved to " << ss.str() << endl;

		}
		drop_bad_phi(pe_lams[i], oe_lams[i], true);
		if (oe_lams[i].shape().first == 0)
		{
			message(1, "all realizations dropped as 'bad' for lambda, scale fac ", vals);
			continue;
		}
		
		ph.update(oe_lams[i], pe_lams[i]);
		
		message(0, "phi summary for lambda, scale fac:", vals,echo);
		ph.report(echo);
		
		mean = ph.get_mean(L2PhiHandler::phiType::COMPOSITE);
		std = ph.get_std(L2PhiHandler::phiType::COMPOSITE);
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
		message(0, "WARNING:  unsuccessful lambda testing, resetting lambda to 10000.0");
		last_best_lam = 10000.0;
		return false;

	}
	double acc_fac = pest_scenario.get_pestpp_options().get_ies_accept_phi_fac();
	double lam_inc = pest_scenario.get_pestpp_options().get_ies_lambda_inc_fac();
	double lam_dec = pest_scenario.get_pestpp_options().get_ies_lambda_dec_fac();
	

	//subset stuff here
	if ((best_idx != -1) && (use_subset) && (subset_size < pe.shape().first))
	{

		double acc_phi = last_best_mean * acc_fac;
		
		if (pest_scenario.get_pestpp_options().get_ies_debug_high_subset_phi())
		{
			cout << "ies_debug_high_subset_phi active" << endl;
			best_mean = acc_phi + 1.0;
		}

		if (best_mean > acc_phi)
		{
			//ph.update(oe_lams[best_idx],pe_lams[best_idx]);
			
			double new_lam = last_best_lam * lam_inc;
			new_lam = (new_lam > lambda_max) ? lambda_max : new_lam;
			last_best_lam = new_lam;
			ss.str("");
			ss << "best subset mean phi  (" << best_mean << ") greater than acceptable phi : " << acc_phi;
			string m = ss.str();
			message(0, m);
			message(1, "abandoning current lambda ensembles, increasing lambda to ", new_lam);
			message(1, "updating realizations with reduced phi");
			update_reals_by_phi(pe_lams[best_idx], oe_lams[best_idx]);
			message(1, "returing to lambda calculations...");
			return false;
		}

		//release the memory of the unneeded pe_lams
		for (int i = 0; i < pe_lams.size(); i++)
		{
			if (i == best_idx)
				continue;
			pe_lams[i] = ParameterEnsemble();
		}
		//need to work out which par and obs en real names to run - some may have failed during subset testing...
		ObservationEnsemble remaining_oe_lam = oe;//copy
		ParameterEnsemble remaining_pe_lam = pe_lams[best_idx];
		vector<string> pe_keep_names, oe_keep_names;
		vector<string> pe_names = pe.get_real_names(), oe_names = oe.get_real_names();

		vector<string> org_pe_idxs,org_oe_idxs;
		set<string> ssub;
		for (auto &i : subset_idxs)
			ssub.emplace(pe_names[i]);
		for (int i=0;i<pe_names.size();i++)
			if (ssub.find(pe_names[i]) == ssub.end())
			{
				pe_keep_names.push_back(pe_names[i]);
				//oe_keep_names.push_back(oe_names[i]);
			}
		ssub.clear();
		for (auto &i : subset_idxs)
			ssub.emplace(oe_names[i]);
		for (int i = 0; i<oe_names.size(); i++)
			if (ssub.find(oe_names[i]) == ssub.end())
			{
				oe_keep_names.push_back(oe_names[i]);
			}
		message(0, "phi summary for best lambda, scale fac: ", vector<double>({ lam_vals[best_idx],scale_vals[best_idx] }));
		ph.update(oe_lams[best_idx], pe_lams[best_idx]);
		ph.report(true);
		message(0, "running remaining realizations for best lambda, scale:", vector<double>({ lam_vals[best_idx],scale_vals[best_idx] }));

		//pe_keep_names and oe_keep_names are names of the remaining reals to eval
		performance_log->log_event("dropping subset idxs from remaining_pe_lam");
		remaining_pe_lam.keep_rows(pe_keep_names);
		performance_log->log_event("dropping subset idxs from remaining_oe_lam");
		remaining_oe_lam.keep_rows(oe_keep_names);
		//save these names for later
		org_pe_idxs = remaining_pe_lam.get_real_names();
		org_oe_idxs = remaining_oe_lam.get_real_names();
		///run
		vector<int> fails = run_ensemble(remaining_pe_lam, remaining_oe_lam);

		//for testing
		if (pest_scenario.get_pestpp_options().get_ies_debug_fail_remainder())
			fails.push_back(0);

		//if any of the remaining runs failed
		if (fails.size() == org_pe_idxs.size())
			throw_ies_error(string("all remaining realizations failed...something is prob wrong"));
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
		//append the remaining obs en
		oe_lam_best.append_other_rows(remaining_oe_lam);
		assert(pe_lams[best_idx].shape().first == oe_lam_best.shape().first);
		drop_bad_phi(pe_lams[best_idx], oe_lam_best);
		if (oe_lam_best.shape().first == 0)
		{
			throw_ies_error(string("all realization dropped after finishing subset runs...something might be wrong..."));
		}
		performance_log->log_event("updating phi");
		ph.update(oe_lam_best, pe_lams[best_idx]);
		best_mean = ph.get_mean(L2PhiHandler::phiType::COMPOSITE);
		best_std = ph.get_std(L2PhiHandler::phiType::COMPOSITE);
		message(1, "phi summary for entire ensemble using lambda,scale_fac ", vector<double>({ lam_vals[best_idx],scale_vals[best_idx] }));
		ph.report(true);
	}
	else
	{
		ph.update(oe_lam_best, pe_lams[best_idx]);
		best_mean = ph.get_mean(L2PhiHandler::phiType::COMPOSITE);
		best_std = ph.get_std(L2PhiHandler::phiType::COMPOSITE);
		
	}

	ph.update(oe_lam_best, pe_lams[best_idx]);
	best_mean = ph.get_mean(L2PhiHandler::phiType::COMPOSITE);
	best_std = ph.get_std(L2PhiHandler::phiType::COMPOSITE);
	message(1, "last best mean phi * acceptable phi factor: ", last_best_mean * acc_fac);
	message(1, "current best mean phi: ", best_mean);

	if (pest_scenario.get_pestpp_options().get_ies_debug_high_upgrade_phi())
	{
		cout << "ies_debug_high_upgrade_phi active" << endl;
		best_mean = (last_best_mean * acc_fac) + 1.0;
	}

	//track this here for phi-based termination check
	best_mean_phis.push_back(best_mean);


	if (best_mean < last_best_mean * acc_fac)
	{
		message(0, "updating parameter ensemble");
		performance_log->log_event("updating parameter ensemble");
		last_best_mean = best_mean;

		pe = pe_lams[best_idx];
		oe = oe_lam_best;
		if (best_std < last_best_std * acc_fac)
		{
			double new_lam = lam_vals[best_idx] * lam_dec;
			new_lam = (new_lam < lambda_min) ? lambda_min : new_lam;
			message(0, "updating lambda to ", new_lam);
			last_best_lam = new_lam;
		}
		else
		{
			message(0, "not updating lambda");
		}
		last_best_std = best_std;
	}

	else 
	{
		//message(0, "not updating parameter ensemble");
		message(0, "only updating realizations with reduced phi");
		update_reals_by_phi(pe_lams[best_idx], oe_lam_best);
		ph.update(oe, pe);
		double new_lam = last_best_lam * lam_inc;
		new_lam = (new_lam > lambda_max) ? lambda_max : new_lam;
		message(0, "incresing lambda to: ", new_lam);
		last_best_lam = new_lam;
	}
	//report_and_save();
	return true;
}


void IterEnsembleSmoother::report_and_save()
{
	ofstream &frec = file_manager.rec_ofstream();
	frec << endl << "  ---  IterEnsembleSmoother iteration " << iter << " report  ---  " << endl;
	frec << "   number of active realizations:  " << pe.shape().first << endl;
	frec << "   number of model runs:           " << run_mgr_ptr->get_total_runs() << endl;

	cout << endl << "  ---  IterEnsembleSmoother iteration " << iter << " report  ---  " << endl;
	cout << "   number of active realizations:   " << pe.shape().first << endl;
	cout << "   number of model runs:            " << run_mgr_ptr->get_total_runs() << endl;

	stringstream ss;
	if (pest_scenario.get_pestpp_options().get_ies_save_binary())
	{
		ss << file_manager.get_base_filename() << "." << iter << ".obs.jcb";
		oe.to_binary(ss.str());
	}
	else
	{
		ss << file_manager.get_base_filename() << "." << iter << ".obs.csv";
		oe.to_csv(ss.str());
	}
	frec << "      current obs ensemble saved to " << ss.str() << endl;
	cout << "      current obs ensemble saved to " << ss.str() << endl;
	ss.str("");
	if (pest_scenario.get_pestpp_options().get_ies_save_binary())
	{
		ss << file_manager.get_base_filename() << "." << iter << ".par.jcb";
		pe.to_binary(ss.str());
	}
	else
	{
		ss << file_manager.get_base_filename() << "." << iter << ".par.csv";
		pe.to_csv(ss.str());
	}
	save_base_real_par_rei(pest_scenario, pe, oe, output_file_writer, file_manager, iter);
	save_base_real_par_rei(pest_scenario, pe, oe, output_file_writer, file_manager, -1);
	//ss << file_manager.get_base_filename() << "." << iter << ".par.csv";
	//pe.to_csv(ss.str());
	frec << "      current par ensemble saved to " << ss.str() << endl;
	cout << "      current par ensemble saved to " << ss.str() << endl;

	

}


void IterEnsembleSmoother::set_subset_idx(int size)
{
	//map<int,int> subset_idx_map;
	subset_idxs.clear();
	int nreal_subset = pest_scenario.get_pestpp_options().get_ies_subset_size();
	if ((!use_subset) || (nreal_subset >= size))
	{
		for (int i = 0; i < size; i++)
			subset_idxs.push_back(i);
		return;
	}
	vector<string> pe_names = pe.get_real_names();

	vector<string>::iterator bidx = find(pe_names.begin(), pe_names.end(), base_name);
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

		for (int i = size-1; i >= 0; i--)
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
		std::uniform_int_distribution<int> uni(0, size-1);
		int idx;
		for (int i = 0; i < 10000000; i++)
		{
			if (subset_idxs.size() >= nreal_subset)
				break;
			idx = uni(subset_rand_gen);
			if (find(subset_idxs.begin(), subset_idxs.end(), idx) != subset_idxs.end())
				continue;
			subset_idxs.push_back(idx);
		}
		if (subset_idxs.size() != nreal_subset)
			throw_ies_error("max iterations exceeded when trying to find random subset idxs");

	}
	else if (how == "PHI_BASED")
	{
		//sidx needs to be index of realization, not realization number
		vector<pair<double, int>> phis;
		//vector<int> sidx;
		int step;
		int idx;
		L2PhiHandler::phiType pt = L2PhiHandler::phiType::COMPOSITE;
		map<string, double>* phi_map = ph.get_phi_map(pt);
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


		step = (phis.size()-1) / nreal_subset;
		//cout << step << endl;
		//cout << (phis.size() - 1) << endl;
		for (i = 1; i < nreal_subset; ++i)
		{
			//add higher phis first
			idx = phis.size() - (i * step);
			if ((subset_idxs.size() < nreal_subset) && (find(subset_idxs.begin(), subset_idxs.end(), phis[idx].second) == subset_idxs.end()))
			{
				subset_idxs.push_back(phis[idx].second);
				//cout << i << endl;
				//cout << idx << endl;
				//cout << phis[idx].first << endl;
				//cout << phis[idx].second << endl;
			}
		}
	}
	else
	{
		//throw runtime_error("unkonwn 'subset_how'");
		throw_ies_error("unknown 'subset_how'");
	}
	stringstream ss;
	for (auto i : subset_idxs)
		ss << i << ":" << pe_names[i] << ", ";
	message(1,"subset idx:pe real name: ",ss.str());
	return;
	//return subset_idx_map;
}

vector<ObservationEnsemble> IterEnsembleSmoother::run_lambda_ensembles(vector<ParameterEnsemble> &pe_lams, vector<double> &lam_vals, vector<double> &scale_vals)
{
	ofstream &frec = file_manager.rec_ofstream();
	stringstream ss;
	ss << "queuing " << pe_lams.size() << " ensembles";
	performance_log->log_event(ss.str());
	run_mgr_ptr->reinitialize();
	
	set_subset_idx(pe_lams[0].shape().first);
	vector<map<int, int>> real_run_ids_vec;
	//ParameterEnsemble pe_lam;
	//for (int i=0;i<pe_lams.size();i++)
	for (auto &pe_lam : pe_lams)
	{
		try
		{
			real_run_ids_vec.push_back(pe_lam.add_runs(run_mgr_ptr,subset_idxs));
		}
		catch (const exception &e)
		{
			stringstream ss;
			ss << "run_ensemble() error queueing runs: " << e.what();
			throw_ies_error(ss.str());
		}
		catch (...)
		{
			throw_ies_error(string("run_ensembles() error queueing runs"));
		}
	}
	performance_log->log_event("making runs");
	try
	{

		run_mgr_ptr->run();
	}
	catch (const exception &e)
	{
		stringstream ss;
		ss << "error running ensembles: " << e.what();
		throw_ies_error(ss.str());
	}
	catch (...)
	{
		throw_ies_error(string("error running ensembles"));
	}

	performance_log->log_event("processing runs");
	vector<int> failed_real_indices;
	vector<ObservationEnsemble> obs_lams;
	//ObservationEnsemble _oe = oe;//copy
	//if (subset_size < pe_lams[0].shape().first)
	//	_oe.keep_rows(subset_idxs);
	map<int, int> real_run_ids;
	//for (auto &real_run_ids : real_run_ids_vec)
	for (int i=0;i<pe_lams.size();i++)
	{
		ObservationEnsemble _oe = oe;//copy
		vector<double> rep_vals{ lam_vals[i],scale_vals[i] };
		real_run_ids = real_run_ids_vec[i];
		//if using subset, reset the real_idx in real_run_ids to be just simple counter
		//if (subset_size < pe_lams[0].shape().first)
		if ((use_subset) && (subset_size < pe_lams[i].shape().first))
		{
			_oe.keep_rows(subset_idxs);
			int ireal = 0;
			map<int, int> temp;
			for (auto &rri : real_run_ids)
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
		catch (const exception &e)
		{
			stringstream ss;
			ss << "error processing runs for lambda,scale: " << lam_vals[i] << ',' << scale_vals[i] << ':' << e.what();
			throw_ies_error(ss.str());
		}
		catch (...)
		{
			stringstream ss;
			ss << "error processing runs for lambda,scale: " << lam_vals[i] << ',' << scale_vals[i];
			throw_ies_error(ss.str());
		}

		if (pest_scenario.get_pestpp_options().get_ies_debug_fail_subset())
			failed_real_indices.push_back(real_run_ids.size()-1);

		if (failed_real_indices.size() > 0)
		{
			stringstream ss;
			vector<string> par_real_names = pe.get_real_names();
			vector<string> obs_real_names = oe.get_real_names();
			vector<string> failed_par_names, failed_obs_names;
			string oname, pname;
			ss << "the following par:obs realization runs failed for lambda,scale " << lam_vals[i] << ',' << scale_vals[i] << "-->";
			for (auto &i : failed_real_indices)
			{
				pname = par_real_names[subset_idxs[i]];
				oname = obs_real_names[subset_idxs[i]];
				failed_par_names.push_back(pname);
				failed_obs_names.push_back(oname);
				ss << pname << ":" << oname << ',';
			}
			string s = ss.str();
			message(1,s);
			if (failed_real_indices.size() == _oe.shape().first)
			{
				message(0, "WARNING: all realizations failed for lambda, scale :", rep_vals);
				_oe = ObservationEnsemble();

			}
			else
			{
				performance_log->log_event("dropping failed realizations");
				//_oe.drop_rows(failed_real_indices);
				//pe_lams[i].drop_rows(failed_real_indices);
				_oe.drop_rows(failed_obs_names);
				pe_lams[i].drop_rows(failed_par_names);
			}

		}
		obs_lams.push_back(_oe);
	}
	return obs_lams;
}


vector<int> IterEnsembleSmoother::run_ensemble(ParameterEnsemble &_pe, ObservationEnsemble &_oe, const vector<int> &real_idxs)
{
	stringstream ss;
	ss << "queuing " << _pe.shape().first << " runs";
	performance_log->log_event(ss.str());
	run_mgr_ptr->reinitialize();
	map<int, int> real_run_ids;
	try
	{
		real_run_ids = _pe.add_runs(run_mgr_ptr,real_idxs);
	}
	catch (const exception &e)
	{
		stringstream ss;
		ss << "run_ensemble() error queueing runs: " << e.what();
		throw_ies_error(ss.str());
	}
	catch (...)
	{
		throw_ies_error(string("run_ensemble() error queueing runs"));
	}
	performance_log->log_event("making runs");
	try
	{
		run_mgr_ptr->run();
	}
	catch (const exception &e)
	{
		stringstream ss;
		ss << "error running ensemble: " << e.what();
		throw_ies_error(ss.str());
	}
	catch (...)
	{
		throw_ies_error(string("error running ensemble"));
	}

	performance_log->log_event("processing runs");
	if (real_idxs.size() > 0)
	{
		_oe.keep_rows(real_idxs);
	}
	vector<int> failed_real_indices;
	try
	{
		failed_real_indices = _oe.update_from_runs(real_run_ids,run_mgr_ptr);
	}
	catch (const exception &e)
	{
		stringstream ss;
		ss << "error processing runs: " << e.what();
		throw_ies_error(ss.str());
	}
	catch (...)
	{
		throw_ies_error(string("error processing runs"));
	}
	//for testing
	//failed_real_indices.push_back(0);

	if (failed_real_indices.size() > 0)
	{
		stringstream ss;
		vector<string> par_real_names = _pe.get_real_names();
		vector<string> obs_real_names = _oe.get_real_names();
		ss << "the following par:obs realization runs failed: ";
		for (auto &i : failed_real_indices)
		{
			ss << par_real_names[i] << ":" << obs_real_names[i] << ',';
		}
		performance_log->log_event(ss.str());
		message(1, "failed realizations: ", failed_real_indices.size());
		string s = ss.str();
		message(1, s);
		performance_log->log_event("dropping failed realizations");
		_pe.drop_rows(failed_real_indices);
		_oe.drop_rows(failed_real_indices);
	}
	return failed_real_indices;
}


void IterEnsembleSmoother::finalize()
{

}

#include <random>
#include <map>
#include <iomanip>
#include <mutex>
#include <thread>
#include <math.h>
#include <algorithm>
#include "Ensemble.h"
#include "RestartController.h"
#include "utilities.h"
#include "Ensemble.h"
#include "DataAssimilator.h"
#include "ObjectiveFunc.h"
#include "covariance.h"
#include "RedSVD-h.h"
#include "SVDPackage.h"
#include "DataAssimilator.h"
#include "EnsembleMethodUtils.h"
//#include <Eigen/SVD>
//using namespace Eigen;



//DataAssimilator::DataAssimilator(Pest& _pest_scenario, FileManager& _file_manager,
//	OutputFileWriter& _output_file_writer, PerformanceLog* _performance_log,
//	RunManagerAbstract* _run_mgr_ptr) : pest_scenario(_pest_scenario), file_manager(_file_manager),
//	output_file_writer(_output_file_writer), performance_log(_performance_log),
//	run_mgr_ptr(_run_mgr_ptr)
//{
//	rand_gen = std::mt19937(pest_scenario.get_pestpp_options().get_random_seed());
//
//	subset_rand_gen = std::mt19937(pest_scenario.get_pestpp_options().get_random_seed());
//
//	pe.set_pest_scenario(&pest_scenario);
//	pe.set_rand_gen(&rand_gen);
//	oe.set_pest_scenario(&pest_scenario);
//	oe.set_rand_gen(&rand_gen);
//	weights.set_pest_scenario(&pest_scenario);
//	localizer.set_pest_scenario(&pest_scenario);
//	icycle = 0;
//	pp_args = pest_scenario.get_pestpp_options().get_passed_args();
//	act_obs_names = pest_scenario.get_ctl_ordered_nz_obs_names();
//	act_par_names = pest_scenario.get_ctl_ordered_adj_par_names();
//
//	
//}

//void DataAssimilator::throw_da_error(string message)
//{
//	performance_log->log_event("DataAssimilator error: " + message);
//	cout << endl << "   ************   " << endl << "    DataAssimilator error: " << message << endl << endl;
//	file_manager.rec_ofstream() << endl << "   ************   " << endl << "    DataAssimilator error: " << message << endl << endl;
//	file_manager.close_file("rec");
//	performance_log->~PerformanceLog();
//	throw runtime_error("DataAssimilator error: " + message);
//}

bool DataAssimilator::initialize_pe(Covariance& cov)
{
	stringstream ss;
	int num_reals;
	
	string par_csv;
	
	par_csv = pest_scenario.get_pestpp_options().get_da_par_csv();
	num_reals = pest_scenario.get_pestpp_options().get_da_num_reals();

	bool drawn = false;
	if (par_csv.size() == 0)
	{
		message(1, "drawing parameter realizations: ", num_reals);
		pe.draw(num_reals, pest_scenario.get_ctl_parameters(), cov, performance_log, pest_scenario.get_pestpp_options().get_ies_verbose_level(),file_manager.rec_ofstream());
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
			catch (const exception & e)
			{
				ss << "error processing par csv: " << e.what();
				throw_em_error(ss.str());
			}
			catch (...)
			{
				throw_em_error(string("error processing par csv"));
			}
		}
		else if ((par_ext.compare("jcb") == 0) || (par_ext.compare("jco") == 0))
		{
			message(1, "loading par ensemble from binary file", par_csv);
			try
			{
				pe.from_binary(par_csv);
			}
			catch (const exception & e)
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
			ss << "unrecognized par csv extension " << par_ext << ", looking for csv, jcb, or jco";
			throw_em_error(ss.str());
		}

		pe.transform_ip(ParameterEnsemble::transStatus::NUM);
		
		if (pp_args.find("DA_NUM_REALS") != pp_args.end())
		{
			int num_reals = pest_scenario.get_pestpp_options().get_da_num_reals();			
	
			if (num_reals < pe.shape().first)
			{
				message(1, "da_num_reals arg passed, truncated parameter ensemble to ", num_reals);
				vector<string> keep_names, real_names = pe.get_real_names();
				for (int i = 0; i < num_reals; i++)
				{
					keep_names.push_back(real_names[i]);
				}
				pe.keep_rows(keep_names);
			}
		}
		
		// todo: Ayman ?? what empirical prior means
		if (pest_scenario.get_pestpp_options().get_ies_use_empirical_prior())
		{
			message(1, "initializing prior parameter covariance matrix from ensemble (using diagonal matrix)");
			parcov = pe.get_diagonal_cov_matrix();
			if (pest_scenario.get_pestpp_options().get_ies_verbose_level() > 1)
			{

				if (pest_scenario.get_pestpp_options().get_save_binary())
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

		//todo: apply this for da
		if (pest_scenario.get_pestpp_options().get_ies_enforce_bounds())
		{
			if (pest_scenario.get_pestpp_options().get_ies_obs_restart_csv().size() > 0)
				message(1, "Warning: even though ies_enforce_bounds is true, a restart obs en was passed, so bounds will not be enforced on the initial par en");
			else
				pe.enforce_bounds(performance_log, false);
		}

	}
	return drawn;

}

//void DataAssimilator::add_bases()
//{
//	//check that 'base' isn't already in ensemble
//	vector<string> rnames = pe.get_real_names();
//	bool inpar = false;
//	if (find(rnames.begin(), rnames.end(), BASE_REAL_NAME) != rnames.end())
//	{
//		message(1, "'base' realization already in parameter ensemble, ignoring '++ies_include_base'");
//		inpar = true;
//	}
//	else
//	{
//		message(1, "adding 'base' parameter values to ensemble");
//		Parameters pars = pest_scenario.get_ctl_parameters();
//		pe.get_par_transform().active_ctl2numeric_ip(pars);
//		vector<int> drop{ pe.shape().first - 1 };
//		pe.drop_rows(drop);
//		pe.append(BASE_REAL_NAME, pars);
//	}
//
//	//check that 'base' isn't already in ensemble
//	rnames = oe.get_real_names();
//	if (find(rnames.begin(), rnames.end(), BASE_REAL_NAME) != rnames.end())
//	{
//		message(1, "'base' realization already in observation ensemble, ignoring '++ies_include_base'");
//	}
//	else
//	{
//		Observations obs = pest_scenario.get_ctl_observations();
//		if (inpar)
//		{
//			vector<string> prnames = pe.get_real_names();
//
//			int idx = find(prnames.begin(), prnames.end(), BASE_REAL_NAME) - prnames.begin();
//			//cout << idx << "," << rnames.size() << endl;
//			string oreal = rnames[idx];
//			stringstream ss;
//			ss << "warning: 'base' realization in par ensenmble but not in obs ensemble," << endl;
//			ss << "         replacing obs realization '" << oreal << "' with 'base'";
//			string mess = ss.str();
//			message(1, mess);
//			vector<string> drop;
//			drop.push_back(oreal);
//			oe.drop_rows(drop);
//			oe.append(BASE_REAL_NAME, obs);
//			//rnames.insert(rnames.begin() + idx, string(base_name));
//			rnames[idx] = BASE_REAL_NAME;
//			oe.reorder(rnames, vector<string>());
//		}
//		else
//		{
//			message(1, "adding 'base' observation values to ensemble");
//			vector<int> drop{ oe.shape().first - 1 };
//			oe.drop_rows(drop);
//			oe.append(BASE_REAL_NAME, obs);
//		}
//	}
//
//	//check that 'base' isn't already in ensemble
//	rnames = weights.get_real_names();
//	if (rnames.size() == 0)
//		return;
//	if (find(rnames.begin(), rnames.end(), BASE_REAL_NAME) != rnames.end())
//	{
//		message(1, "'base' realization already in weights ensemble, ignoring '++ies_include_base'");
//	}
//	else
//	{
//		//Observations obs = pest_scenario.get_ctl_observations();
//		ObservationInfo oinfo = pest_scenario.get_ctl_observation_info();
//		Eigen::VectorXd q;
//		//vector<string> act_obs_names = pest_scenario->get_ctl_ordered_nz_obs_names();
//		vector<string> vnames = weights.get_var_names();
//		q.resize(vnames.size());
//		double w;
//		for (int i = 0; i < vnames.size(); i++)
//		{
//			q(i) = oinfo.get_weight(vnames[i]);
//		}
//
//		Observations wobs(vnames, q);
//		if (inpar)
//		{
//			vector<string> prnames = pe.get_real_names();
//
//			int idx = find(prnames.begin(), prnames.end(), BASE_REAL_NAME) - prnames.begin();
//			//cout << idx << "," << rnames.size() << endl;
//			string oreal = rnames[idx];
//			stringstream ss;
//			ss << "warning: 'base' realization in par ensenmble but not in weights ensemble," << endl;
//			ss << "         replacing weights realization '" << oreal << "' with 'base'";
//			string mess = ss.str();
//			message(1, mess);
//			vector<string> drop;
//			drop.push_back(oreal);
//			weights.drop_rows(drop);
//			weights.append(BASE_REAL_NAME, wobs);
//			//rnames.insert(rnames.begin() + idx, string(base_name));
//			rnames[idx] = BASE_REAL_NAME;
//			weights.reorder(rnames, vector<string>());
//		}
//		else
//		{
//			message(1, "adding 'base' weight values to weights");
//
//
//			weights.append(BASE_REAL_NAME, wobs);
//		}
//	}
//}

bool DataAssimilator::initialize_oe(Covariance& cov)
{
	//if there are no active obs, then just reserve a generic oe and return
	if (act_obs_names.size() == 0)
	{
		oe.reserve(pe.get_real_names(), pest_scenario.get_ctl_ordered_obs_names());
		return true;
	}


	stringstream ss;
	int num_reals = pe.shape().first;

	string obs_csv;
	obs_csv = da_ctl_params.get_svalue("DA_OBSERVATION_ENSEMBLE");
	
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
			oe.draw(num_reals, cov, performance_log, pest_scenario.get_pestpp_options().get_ies_verbose_level(),file_manager.rec_ofstream());

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
			catch (const exception & e)
			{
				ss << "error processing obs csv: " << e.what();
				throw_em_error(ss.str());
			}
			catch (...)
			{
				throw_em_error(string("error processing obs csv"));
			}
		}
		else if ((obs_ext.compare("jcb") == 0) || (obs_ext.compare("jco") == 0))
		{
			message(1, "loading obs ensemble from binary file", obs_csv);
			try
			{
				oe.from_binary(obs_csv);
			}
			catch (const exception & e)
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
			ss << "unrecognized obs ensemble extension " << obs_ext << ", looing for csv, jcb, or jco";
			throw_em_error(ss.str());
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
				message(1, "ies_num_reals arg passed, truncated observation ensemble to ", num_reals);
				vector<string> keep_names, real_names = oe.get_real_names();
				for (int i = 0; i < num_reals; i++)
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
//void DataAssimilator::message(int level, char* _message, vector<T, A> _extras)
//{
//	string s(_message);
//	message(level, s, _extras);
//}

//void DataAssimilator::message(int level, char* _message)
//{
//	string s(_message);
//	message(level, s);
//}

//template<typename T>
//void DataAssimilator::message(int level, char* _message, T extra)
//{
//	string s(_message);
//	message(level, s, extra);
//
//}

template<typename T, typename A>
void DataAssimilator::message(int level, const string& _message, vector<T, A> _extras, bool echo)
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
	performance_log->log_event(ss.str());

}

void DataAssimilator::message(int level, const string& _message)
{
	message(level, _message, vector<string>());
}

template<typename T>
void DataAssimilator::message(int level, const string& _message, T extra)
{
	stringstream ss;
	ss << _message << " " << extra;
	string s = ss.str();
	message(level, s);
}

void DataAssimilator::sanity_checks()
{
	PestppOptions* ppo = pest_scenario.get_pestpp_options_ptr();
	vector<string> errors;
	vector<string> warnings;
	stringstream ss;
	

	if (pest_scenario.get_control_info().pestmode == ControlInfo::PestMode::UNKNOWN)
	{
		warnings.push_back("unrecognized 'pestmode', using 'estimation'");
	}


	if (pest_scenario.get_ctl_ordered_pi_names().size() > 0)
	{
		warnings.push_back("prior information equations not supported in pestpp-da, ignoring...");
	}
	
	if ((da_type != "VANILLA") && (da_type != "ITERATIVE") && (da_type != "MDA"))
	{
		ss.str("");
		ss << " DA_TYPE is '" << da_type << "' but must be 'VANILLA','ITERATIVE', or'MDA'";
		errors.push_back(ss.str());
	}

	if (warnings.size() > 0)
	{
		message(0, "sanity_check warnings");
		for (auto& w : warnings)
			message(1, w);
		message(1, "continuing initialization...");
	}
	if (errors.size() > 0)
	{
		message(0, "sanity_check errors - uh oh");
		for (auto& e : errors)
			message(1, e);
		throw_em_error(string("sanity_check() found some problems - please review rec file"));
	}
	//cout << endl << endl;
}

void DataAssimilator::initialize_restart()
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
		catch (const exception & e)
		{
			ss << "error processing restart obs csv: " << e.what();
			throw_em_error(ss.str());
		}
		catch (...)
		{
			throw_em_error(string("error processing restart obs csv"));
		}
	}
	else if ((obs_ext.compare("jcb") == 0) || (obs_ext.compare("jco") == 0))
	{
		message(1, "loading restart obs ensemble from binary file", obs_restart_csv);
		try
		{
			oe.from_binary(obs_restart_csv);
		}
		catch (const exception & e)
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
		ss << "unrecognized restart obs ensemble extension " << obs_ext << ", looking for csv, jcb, or jco";
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
			catch (const exception & e)
			{
				ss << "error processing restart par csv: " << e.what();
				throw_em_error(ss.str());
			}
			catch (...)
			{
				throw_em_error(string("error processing restart par csv"));
			}
		}
		else if ((par_ext.compare("jcb") == 0) || (par_ext.compare("jco") == 0))
		{
			message(1, "loading restart par ensemble from binary file", par_restart_csv);
			try
			{
				pe.from_binary(par_restart_csv);
			}
			catch (const exception & e)
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
			ss << "unrecognized restart par ensemble extension " << par_ext << ", looking for csv, jcb, or jco";
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
	for (auto& rname : oe_real_names)
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
				throw_em_error(ss.str());
			}

		}
		else
		{
			ss << "the following realization names were found in the restart obs en but not in the 'base' obs en:";
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
		catch (exception & e)
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
			catch (exception & e)
			{
				ss << "error reordering oe_base with restart oe:" << e.what();
				throw_em_error(ss.str());
			}
			catch (...)
			{
				throw_em_error(string("error reordering oe_base with restart oe"));
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
			catch (exception & e)
			{
				ss << "error reordering pe_base with restart pe:" << e.what();
				throw_em_error(ss.str());
			}
			catch (...)
			{
				throw_em_error(string("error reordering pe_base with restart pe"));
			}
		}


		/*try
		{
			pe.reorder(pe_real_names, vector<string>());
		}
		catch (exception &e)
		{
			ss << "error reordering pe with restart oe:" << e.what();
			throw_em_error(ss.str());
		}
		catch (...)
		{
			throw_em_error(string("error reordering pe with restart oe"));
		}*/
	}
	else if (oe.shape().first > oe_base.shape().first) //something is wrong
	{
		ss << "restart oe has too many rows: " << oe.shape().first << " compared to oe_base: " << oe_base.shape().first;
		throw_em_error(ss.str());
	}
}



void DataAssimilator::initialize_parcov()
{
	stringstream ss;
	performance_log->log_event("initializing parcov");

	if (pest_scenario.get_pestpp_options().get_ies_use_empirical_prior())
		return;
	string how = parcov.try_from(pest_scenario, file_manager);
	message(1, "parcov loaded ", how);
	//if (parcov.e_ptr()->rows() > 0)
	if (act_par_names.size() > 0)
	{
		parcov = parcov.get(act_par_names);
	}

}


void DataAssimilator::initialize_obscov()
{
	if (act_obs_names.size() == 0)
		message(1, "no non-zero weighted observations for cycle ", icycle);
	else
	{
		message(1, "initializing observation noise covariance matrix");
		string obscov_filename = pest_scenario.get_pestpp_options().get_obscov_filename();

		string how = obscov.try_from(pest_scenario, file_manager, false);
		message(1, "obscov loaded ", how);
		obscov = obscov.get(act_obs_names);
	}
}


//void DataAssimilator::initialize(int _icycle)
//{
//	icycle = _icycle;
//	message(0, "initializing cycle", icycle);
//	stringstream ss;
//	ofstream& f_rec = file_manager.rec_ofstream();
//
//	CtlPar_container da_ctl_params = pest_scenario.get_pestpp_options().da_ctl_params;
//
//	pp_args = pest_scenario.get_pestpp_options().get_passed_args();
//	act_obs_names = pest_scenario.get_ctl_ordered_nz_obs_names();
//	act_par_names = pest_scenario.get_ctl_ordered_adj_par_names();
//
//
//	// run the model one time using the initial parameters values. 
//	if (pest_scenario.get_control_info().noptmax == 0)
//	{
//		message(0, "'noptmax'=0, running control file parameter values and quitting");
//
//		Parameters pars = pest_scenario.get_ctl_parameters();
//		ParamTransformSeq pts = pe.get_par_transform();
//
//		ParameterEnsemble _pe(&pest_scenario, &rand_gen);
//		_pe.reserve(vector<string>(), pest_scenario.get_ctl_ordered_par_names());
//		_pe.set_trans_status(ParameterEnsemble::transStatus::CTL);
//		_pe.append("BASE", pars);
//		string par_csv = file_manager.get_base_filename() + ".par.csv";
//		pe_base = _pe;
//		pe_base.reorder(vector<string>(), act_par_names);
//		ObservationEnsemble _oe(&pest_scenario, &rand_gen);
//		_oe.reserve(vector<string>(), pest_scenario.get_ctl_ordered_obs_names());
//		_oe.append("BASE", pest_scenario.get_ctl_observations());
//		oe_base = _oe;
//		oe_base.reorder(vector<string>(), act_obs_names);
//		//initialize the phi handler
//		ph = L2PhiHandler(&pest_scenario, &file_manager, &oe_base, &pe_base, &parcov);
//		if (ph.get_lt_obs_names().size() > 0)
//		{
//			message(1, "less_than inequality defined for observations: ", ph.get_lt_obs_names().size());
//		}
//		if (ph.get_gt_obs_names().size())
//		{
//			message(1, "greater_than inequality defined for observations: ", ph.get_gt_obs_names().size());
//		}
//		message(1, "running control file parameter values");
//
//		vector<int> failed_idxs = run_ensemble(_pe, _oe);
//		if (failed_idxs.size() != 0)
//		{
//			message(0, "control file parameter value run failed...bummer");
//			throw_em_error("control file parameter value run failed");
//		}
//		string obs_csv = file_manager.get_base_filename() + ".obs.csv";
//		message(1, "saving results from control file parameter value run to ", obs_csv);
//		_oe.to_csv(obs_csv);
//
//		ph.update(_oe, _pe);
//		message(0, "control file parameter phi report:");
//		ph.report(true);
//		ph.write(0, 1);
//		ObjectiveFunc obj_func(&(pest_scenario.get_ctl_observations()), &(pest_scenario.get_ctl_observation_info()), &(pest_scenario.get_prior_info()));
//		Observations obs;
//		Eigen::VectorXd v = _oe.get_real_vector("BASE");
//		vector<double> vv;
//		vv.resize(v.size());
//		Eigen::VectorXd::Map(&vv[0], v.size()) = v;
//		obs.update(_oe.get_var_names(), vv);
//
//		// save parameters to .par file
//		output_file_writer.write_par(file_manager.open_ofile_ext("base.par"), pars, *(pts.get_offset_ptr()),
//			*(pts.get_scale_ptr()));
//		file_manager.close_file("par");
//
//		// save new residuals to .rei file
//		output_file_writer.write_rei(file_manager.open_ofile_ext("base.rei"), 0,
//			pest_scenario.get_ctl_observations(), obs, obj_func, pars);
//
//		return;
//	}
//
//	// TODO: (Ayman) Maybe here we can add option to make one Monte Carlo forward run for all realizations
//
//	//set some defaults
//	PestppOptions* ppo = pest_scenario.get_pestpp_options_ptr();
//
//	// TODO: (Ayman) We need to develop a da verbose level, for now ies is being used
//	verbose_level = pest_scenario.get_pestpp_options_ptr()->get_ies_verbose_level();
//	if (pest_scenario.get_n_adj_par() >= 1e6)
//	{
//		message(0, "You are a god among mere mortals!");
//	}
//
//	/*PestppOptions::SVD_PACK svd = ppo->get_svd_pack();
//	if (svd == PestppOptions::SVD_PACK::PROPACK)
//	{
//		message(1, "using PROPACK for truncated svd solve");
//	}
//	else
//	{*/
//	message(1, "using REDSVD for truncated svd solve");
//	//}
//	message(1, "maxsing:", pest_scenario.get_svd_info().maxsing);
//	message(1, "eigthresh: ", pest_scenario.get_svd_info().eigthresh);
//
//	//if ((ppo->get_ies_localizer().size() > 0) & (ppo->get_ies_autoadaloc()))
//	//{
//	//	throw_em_error("use of localization matrix and autoadaloc not supported...yet!");
//	//}
//	message(1, "initializing localizer--");
//	use_localizer = localizer.initialize(performance_log);
//	num_threads = pest_scenario.get_pestpp_options().get_ies_num_threads();
//	if (!use_localizer)
//		message(1, "not using localization");
//	else
//	{
//		if (localizer.get_autoadaloc())
//		{
//			message(1, "using automatic adaptive localization");
//			message(2, "with autoadaloc_sigma_dist ", ppo->get_ies_autoadaloc_sigma_dist());
//		}
//		if (localizer.get_filename().size() > 0)
//		{
//			message(1, "using localization matrix " + localizer.get_filename());
//			localizer.report(file_manager.rec_ofstream());
//		}
//	}
//	if ((use_localizer) && (!localizer.get_autoadaloc()))
//	{
//
//		ss.str("");
//		ss << "using localized solution with " << localizer.get_num_upgrade_steps() << " sequential upgrade steps";
//		message(1, ss.str());
//		ss.str("");
//		if (num_threads > 0)
//		{
//			ss.str("");
//			ss << "using multithreaded localization calculation with " << num_threads << " threads";
//			message(1, ss.str());
//
//		}
//		if (localizer.get_how() == Localizer::How::OBSERVATIONS)
//			message(1, "localizing by obseravtions");
//		else
//			message(1, "localizing by parameters");
//	}
//
//
//
//	iter = 0;
//	//ofstream &frec = file_manager.rec_ofstream();
//	last_best_mean = 1.0E+30;
//	last_best_std = 1.0e+30;
//	lambda_max = 1.0E+30;
//	lambda_min = 1.0E-30;
//	warn_min_reals = 10;
//	error_min_reals = 2;
//	consec_bad_lambda_cycles = 0;
//
//	if (use_ies)
//	{	
//		lam_mults = pest_scenario.get_pestpp_options().get_ies_lam_mults();
//		if (lam_mults.size() == 0)
//			lam_mults.push_back(1.0);
//		message(1, "using lambda multipliers: ", lam_mults);
//		vector<double> scale_facs = pest_scenario.get_pestpp_options().get_lambda_scale_vec();
//		message(1, "using lambda scaling factors: ", scale_facs);
//		double acc_fac = pest_scenario.get_pestpp_options().get_ies_accept_phi_fac();
//		message(1, "acceptable phi factor: ", acc_fac);
//		double inc_fac = pest_scenario.get_pestpp_options().get_ies_lambda_inc_fac();
//		message(1, "lambda increase factor: ", inc_fac);
//		double dec_fac = pest_scenario.get_pestpp_options().get_ies_lambda_dec_fac();
//		message(1, "lambda decrease factor: ", dec_fac);
//		message(1, "max run fail: ", ppo->get_max_run_fail());
//	}
//	else
//	{
//		// report DA paramters use
//		lam_mults = da_ctl_params.get_vvalue("DA_LAMBDA_MULTS");
//
//
//	}
//
//	sanity_checks();
//
//
//	bool echo = false;
//	if (verbose_level > 1)
//		echo = true;
//
//	// ies_use_empirical --- TODO (Ayman: check if we need to add if-statement for use only in ies
//	initialize_parcov();	
//	initialize_obscov();
//
//	subset_size = pest_scenario.get_pestpp_options().get_ies_subset_size();
//	reg_factor = pest_scenario.get_pestpp_options().get_ies_reg_factor();
//	message(1, "using reg_factor: ", reg_factor);
//	double bad_phi = pest_scenario.get_pestpp_options().get_ies_bad_phi();
//	if (bad_phi < 1.0e+30)
//		message(1, "using bad_phi: ", bad_phi);
//
//	int num_reals;
//	if (use_ies)
//		num_reals = pest_scenario.get_pestpp_options().get_ies_num_reals();
//	else
//		num_reals = pest_scenario.get_pestpp_options().get_da_num_reals();
//
//	//if (icycle == 0) 
//	//{
//	//	// Get Par realizations either from external file or by drawing internally
//	//	pe_drawn = initialize_pe(parcov);	
//	//	
//	//}
//
//	if (pest_scenario.get_pestpp_options().get_ies_use_prior_scaling())
//	{
//		message(1, "forming inverse sqrt of prior parameter covariance matrix");
//
//		if (parcov.isdiagonal())
//			parcov_inv_sqrt = parcov.inv(echo).get_matrix().diagonal().cwiseSqrt().asDiagonal();
//		else
//		{
//			message(1, "first extracting diagonal from prior parameter covariance matrix");
//			Covariance parcov_diag;
//			parcov_diag.from_diagonal(parcov);
//			parcov_inv_sqrt = parcov_diag.inv(echo).get_matrix().diagonal().cwiseSqrt().asDiagonal();
//		}
//	}
//	else 
//	{
//		message(1, "not using prior parameter covariance matrix scaling");
//	}
//
//	oe_drawn = initialize_oe(obscov);
//	string center_on = ppo->get_ies_center_on();
//
//	if (icycle == 0)
//	{
//		try
//		{
//			pe.check_for_dups();
//		}
//		catch (const exception & e)
//		{
//			string message = e.what();
//			throw_em_error("error in parameter ensemble: " + message);
//		}
//
//	}
//
//	try
//	{
//		oe.check_for_dups();
//	}
//	catch (const exception & e)
//	{
//		string message = e.what();
//		throw_em_error("error in observation ensemble: " + message);
//	}
//
//	if (pe.shape().first != oe.shape().first)
//	{
//		//the special case where par en < obs en and all par reals are found in obs en...
//
//		if (pe.shape().first < oe.shape().first)
//		{
//			vector<string> oe_names = oe.get_real_names();
//			set<string> oset(oe_names.begin(), oe_names.end());
//			vector<string> missing;
//			for (auto n : pe.get_real_names())
//				if (oset.find(n) == oset.end())
//				{
//					missing.push_back(n);
//				}
//			if (missing.size() == 0)
//			{
//				ss.str("");
//				ss << "par en has " << pe.shape().first << " realizations, compared to " << oe.shape().first << " obs realizations";
//				message(1, ss.str());
//				message(1, " the realization names are compatible");
//				message(1, "re-indexing obs en to align with par en...");
//
//				oe.reorder(pe.get_real_names(), vector<string>());
//			}
//			else
//			{
//				ss.str("");
//				ss << "the following par en real names were not found in the obs en: ";
//				for (auto m : missing)
//				{
//					ss << m << ",";
//				}
//				throw_em_error(ss.str());
//
//			}
//		}
//		else
//		{
//			ss.str("");
//			ss << "parameter ensemble rows (" << pe.shape().first << ") not equal to observation ensemble rows (" << oe.shape().first << ")";
//			throw_em_error(ss.str());
//		}
//	}
//
//		
//	//need this here for Am calcs...
//	//message(1, "transforming parameter ensemble to numeric");
//	//if (icycle == 0)
//		// parameters need to be transformed once ?
//		//pe.transform_ip(ParameterEnsemble::transStatus::NUM);
//
//	if (use_ies)
//	{
//		if (pest_scenario.get_pestpp_options().get_ies_include_base())
//			if (pp_args.find("IES_RESTART_OBS_EN") != pp_args.end())
//			{
//				message(1, "Warning: even though `ies_include_base` is true, you passed a restart obs en, not adding 'base' realization...");
//			}
//			else
//				add_bases();
//	}
//	else
//	{	
//		
//		if (da_ctl_params.get_bvalue("DA_ADD_BASE"))
//			if (pp_args.find("DA_RESTART_OBS_EN") != pp_args.end())  // todo: think about how this restart works in da?
//			{
//				message(1, "Warning: even though `da_include_base` is true, you passed a restart obs en, not adding 'base' realization...");
//			}
//			else
//				add_bases();
//
//	}
//	message(2, "checking for denormal values in pe");
//	//if (icycle == 0)
//	pe.check_for_normal("initial transformed parameter ensemble");
//	/* Ayman commented this and moved it after we update dynamic states exist in pe
//	ss.str("");
//	if (pest_scenario.get_pestpp_options().get_ies_save_binary())
//	{
//		ss << file_manager.get_base_filename() << "_cycle_" << icycle << ".par.jcb";
//		pe.to_binary(ss.str());
//	}
//	else
//	{
//		ss << file_manager.get_base_filename() << "_cycle_" << icycle << ".par.csv";
//		pe.to_csv(ss.str());
//	}
//	*/
//
//
//	//message(1, "saved initial parameter ensemble to ", ss.str());
//	if (act_obs_names.size() > 0)
//	{
//		message(2, "checking for denormal values in base oe");
//		oe.check_for_normal("base observation ensemble");
//		ss.str("");
//		if (pest_scenario.get_pestpp_options().get_ies_save_binary())
//		{
//			ss << file_manager.get_base_filename() << ".base.obs.jcb";
//			oe.to_binary(ss.str());
//		}
//		else
//		{
//			ss << file_manager.get_base_filename() << ".base.obs.csv";
//			oe.to_csv(ss.str());
//		}
//		message(1, "saved base observation ensemble (obsval+noise) to ", ss.str());
//	}
//
//	if (center_on.size() > 0)
//	{
//		ss.str("");
//		ss << "centering on realization: '" << center_on << "' ";
//		message(1, ss.str());
//		vector<string> names = pe.get_real_names();
//		if (find(names.begin(), names.end(), center_on) == names.end())
//			throw_em_error("'ies_center_on' realization not found in par en: " + center_on);
//		names = oe.get_real_names();
//		if (find(names.begin(), names.end(), center_on) == names.end())
//			throw_em_error("'ies_center_on' realization not found in obs en: " + center_on);
//	}
//	else
//		message(1, "centering on ensemble mean vector");
//
//	if (pest_scenario.get_control_info().noptmax == -2)
//	{
//		message(0, "'noptmax'=-2, running mean parameter ensemble values and quitting");
//		message(1, "calculating mean parameter values");
//		Parameters pars;
//		vector<double> mv = pe.get_mean_stl_vector();
//		pars.update(pe.get_var_names(), pe.get_mean_stl_vector());
//		ParamTransformSeq pts = pe.get_par_transform();
//
//		ParameterEnsemble _pe(&pest_scenario, &rand_gen);
//		_pe.reserve(vector<string>(), pe.get_var_names());
//		_pe.set_trans_status(pe.get_trans_status());
//		_pe.append("mean", pars);
//		string par_csv = file_manager.get_base_filename() + ".mean.par.csv";
//		message(1, "saving mean parameter values to ", par_csv);
//		_pe.to_csv(par_csv);
//		pe_base = _pe;
//		pe_base.reorder(vector<string>(), act_par_names);
//		ObservationEnsemble _oe(&pest_scenario, &rand_gen);
//		_oe.reserve(vector<string>(), oe.get_var_names());
//		_oe.append("mean", pest_scenario.get_ctl_observations());
//		oe_base = _oe;
//		oe_base.reorder(vector<string>(), act_obs_names);
//		//initialize the phi handler
//		ph = L2PhiHandler(&pest_scenario, &file_manager, &oe_base, &pe_base, &parcov);
//		if (ph.get_lt_obs_names().size() > 0)
//		{
//			message(1, "less_than inequality defined for observations: ", ph.get_lt_obs_names().size());
//		}
//		if (ph.get_gt_obs_names().size())
//		{
//			message(1, "greater_than inequality defined for observations: ", ph.get_gt_obs_names().size());
//		}
//		message(1, "running mean parameter values");
//
//		vector<int> failed_idxs = run_ensemble(_pe, _oe);
//		if (failed_idxs.size() != 0)
//		{
//			message(0, "mean parameter value run failed...bummer");
//			return;
//		}
//		string obs_csv = file_manager.get_base_filename() + ".mean.obs.csv";
//		message(1, "saving results from mean parameter value run to ", obs_csv);
//		_oe.to_csv(obs_csv);
//
//		ph.update(_oe, _pe);
//		message(0, "mean parameter phi report:");
//		ph.report();
//
//		return;
//	}
//
//	if (subset_size > pe.shape().first)
//	{
//		use_subset = false;
//	}
//	else
//	{
//		message(1, "using subset in lambda testing, number of realizations used in subset testing: ", subset_size);
//		string how = pest_scenario.get_pestpp_options().get_ies_subset_how();
//		message(1, "subset how: ", how);
//		use_subset = true;
//	}
//
//	oe_org_real_names = oe.get_real_names();
//	pe_org_real_names = pe.get_real_names();
//	string obs_restart_csv = pest_scenario.get_pestpp_options().get_ies_obs_restart_csv();
//	string par_restart_csv = pest_scenario.get_pestpp_options().get_ies_par_restart_csv();
//
//	oe_base = oe; //copy
//	//reorder this for later...
//	oe_base.reorder(vector<string>(), act_obs_names);
//
//
//	pe_base = pe; //copy
//	//reorder this for later
//	pe_base.reorder(vector<string>(), act_par_names);
//	
//	//the hard way to restart
//	if ((act_obs_names.size() > 0) && (obs_restart_csv.size() > 0))
//		initialize_restart();
//
//	//no restart
//	else
//	{
//		performance_log->log_event("running initial ensemble");
//		message(1, "running initial ensemble of size", oe.shape().first);
//		vector<int> failed = run_ensemble(pe, oe);
//
//		if (pe.shape().first == 0)
//			throw_em_error("all realizations failed during initial evaluation");
//		if (failed.size() > 0)
//			oe_base.drop_rows(failed);
//		
//
//		//string obs_csv = file_manager.get_base_filename() + ".0.obs.csv";
//		//message(1, "saving results of initial ensemble run to", obs_csv);
//		//oe.to_csv(obs_csv);
//		/*if (icycle == 0)
//		{
//			pe.transform_ip(ParameterEnsemble::transStatus::NUM);
//		}*/
//	}
//
//	// extract dynamic states forecast from oe and add them to pe.
//	//States are flaged as ones that have zero weights and exist in both oe and pe
//	vector <string> dyn_states;
//	dyn_states = get_dynamic_states();	
//	add_dynamic_state_to_pe();
//	ss.str("");
//	if (pest_scenario.get_pestpp_options().get_ies_save_binary())
//	{
//		ss << file_manager.get_base_filename() << "_cycle_" << icycle << ".0.par.jcb";
//		pe.to_binary(ss.str());
//	}
//	else
//	{
//		ss << file_manager.get_base_filename() << "_cycle_" << icycle << ".par.csv";
//		pe.to_csv(ss.str());
//	}
//	
//	ss.str("");
//	if (pest_scenario.get_pestpp_options().get_ies_save_binary())
//	{
//		ss << file_manager.get_base_filename() << ".0.obs.jcb";
//		oe.to_binary(ss.str());
//	}
//	else
//	{
//		ss << file_manager.get_base_filename() << ".0.obs.csv";
//		oe.to_csv(ss.str());
//	}
//	message(1, "saved initial obs ensemble to", ss.str());
//
//
//	performance_log->log_event("calc pre-drop phi");
//	//initialize the phi handler
//	ph = L2PhiHandler(&pest_scenario, &file_manager, &oe_base, &pe_base, &parcov);
//
//	if (ph.get_lt_obs_names().size() > 0) // todo: (Ayman) how to handle those obs in DA? 
//	{
//		message(1, "less_than inequality defined for observations: ", ph.get_lt_obs_names().size());
//	}
//	if (ph.get_gt_obs_names().size())
//	{
//		message(1, "greater_than inequality defined for observations: ", ph.get_gt_obs_names().size());
//	}
//
//	ph.update(oe, pe);
//	message(0, "pre-drop initial phi summary");
//	ph.report();
//	drop_bad_phi(pe, oe);
//	if (oe.shape().first == 0)
//	{
//		throw_em_error(string("all realizations dropped as 'bad'"));
//	}
//	if (oe.shape().first <= error_min_reals)
//	{
//		message(0, "too few active realizations:", oe.shape().first);
//		message(1, "need at least ", error_min_reals);
//		throw_em_error(string("too few active realizations, cannot continue"));
//	}
//	if (oe.shape().first < warn_min_reals)
//	{
//		ss.str("");
//		ss << "WARNING: less than " << warn_min_reals << " active realizations...might not be enough";
//		string s = ss.str();
//		message(0, s);
//	}
//
//
//	pcs = ParChangeSummarizer(&pe_base, &file_manager, &output_file_writer);
//	vector<string> in_conflict = detect_prior_data_conflict();
//	if (in_conflict.size() > 0)
//	{
//		ss.str("");
//		ss << "WARNING: " << in_conflict.size() << " non-zero weighted observations are in conflict";
//		ss << " with the prior simulated ensemble." << endl;
//		message(0, ss.str());
//	}
//
//	cout << "...see rec file for listing of conflicted observations" << endl << endl;
//	ofstream& frec = file_manager.rec_ofstream();
//	frec << endl << "...conflicted observations: " << endl;
//	for (auto oname : in_conflict)
//	{
//		frec << oname << endl;
//	}
//	if (!ppo->get_ies_drop_conflicts())
//	{
//		ss.str("");
//		ss << "  Continuing with data assimilation will likely result ";
//		ss << " in parameter bias and, ultimately, forecast bias";
//		message(1, ss.str());
//	}
//	else
//	{
//
//		//check that all obs are in conflict
//		message(1, "dropping conflicted observations");
//		if (in_conflict.size() == oe.shape().second)
//		{
//			throw_em_error("all non-zero weighted observations in conflict state, cannot continue");
//		}
//		//drop from act_obs_names
//		vector<string> t;
//		set<string> sconflict(in_conflict.begin(), in_conflict.end());
//		for (auto oname : act_obs_names)
//			if (sconflict.find(oname) == sconflict.end())
//				t.push_back(oname);
//		act_obs_names = t;
//
//		//update obscov
//		obscov.drop(in_conflict);
//
//		//drop from oe_base
//		oe_base.drop_cols(in_conflict);
//		//shouldnt need to update localizer since we dropping not adding
//		//updating weights in control file
//
//		ObservationInfo* oi = pest_scenario.get_observation_info_ptr();
//		int org_nnz_obs = pest_scenario.get_ctl_ordered_nz_obs_names().size();
//		for (auto n : in_conflict)
//		{
//			oi->set_weight(n, 0.0);
//		}
//
//		stringstream ss;
//		ss << "number of non-zero weighted observations reduced from " << org_nnz_obs;
//		ss << " to " << pest_scenario.get_ctl_ordered_nz_obs_names().size() << endl;
//		message(1, ss.str());
//	}
//	performance_log->log_event("calc initial phi");
//	ph.update(oe, pe);
//	message(0, "initial phi summary");
//	ph.report();
//	ph.write(0, run_mgr_ptr->get_total_runs());
//	best_mean_phis.push_back(ph.get_mean(L2PhiHandler::phiType::COMPOSITE));
//	if (!pest_scenario.get_pestpp_options().get_ies_use_approx())
//	{
//		message(1, "using full (MAP) update solution");
//
//	}
//
//	last_best_mean = ph.get_mean(L2PhiHandler::phiType::COMPOSITE);
//	last_best_std = ph.get_std(L2PhiHandler::phiType::COMPOSITE);
//	last_best_lam = pest_scenario.get_pestpp_options().get_ies_init_lam();
//	if (last_best_lam <= 0.0)
//	{
//		//double x = last_best_mean / (2.0 * double(oe.shape().second));
//		double x = last_best_mean / (2.0 * double(pest_scenario.get_ctl_ordered_nz_obs_names().size()));
//		last_best_lam = pow(10.0, (floor(log10(x))));
//	}
//	message(1, "current lambda:", last_best_lam);
//	message(0, "initialization complete");
//}

void DataAssimilator::forward_run_noptmax_0(int icycle)
{
	// Dry Run. No parameter update occurs. Run initial values and move to the next time cycle 
	vector <string> dyn_states;
	stringstream ss;

	message(0, "'noptmax'=0, running control file parameter values and quitting");

	Parameters pars = pest_scenario.get_ctl_parameters();
	ParamTransformSeq pts = pe.get_par_transform();

	ParameterEnsemble _pe(&pest_scenario, &rand_gen);
	_pe.reserve(vector<string>(), pest_scenario.get_ctl_ordered_par_names());
	_pe.set_trans_status(ParameterEnsemble::transStatus::CTL);
	_pe.append(BASE_REAL_NAME, pars);
	if (icycle == 0)
	{
		pe = _pe; // pe is initialized in the first cycle
	}
	else
	{
		_pe = pe;
	}
	string par_csv = file_manager.get_base_filename() + ".par.csv";
	pe_base = _pe;
	pe_base.reorder(vector<string>(), act_par_names);
	ObservationEnsemble _oe(&pest_scenario, &rand_gen);
	_oe.reserve(vector<string>(), pest_scenario.get_ctl_ordered_obs_names());
	_oe.append(BASE_REAL_NAME, pest_scenario.get_ctl_observations());
	oe_base = _oe;
	oe = _oe; // 
			  // 
	oe_base.reorder(vector<string>(), act_obs_names);

	//initialize the phi handler
	ss.str("");
	ss << icycle << ".";
	ph = L2PhiHandler(&pest_scenario, &file_manager, &oe_base, &pe_base, &parcov, true, ss.str());
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
		throw_em_error("control file parameter value run failed");
	}
	string obs_csv = file_manager.get_base_filename() + ".obs.csv";
	message(1, "saving results from control file parameter value run to ", obs_csv);
	_oe.to_csv(obs_csv);

	ph.update(_oe, _pe);
	message(0, "control file parameter phi report:");
	ph.report(true);
	ph.write(0, 1);
	ObjectiveFunc obj_func(&(pest_scenario.get_ctl_observations()), &(pest_scenario.get_ctl_observation_info()), &(pest_scenario.get_prior_info()));
	Observations obs;
	Eigen::VectorXd v = _oe.get_real_vector(BASE_REAL_NAME);
	vector<double> vv;
	vv.resize(v.size());
	Eigen::VectorXd::Map(&vv[0], v.size()) = v;
	obs.update(_oe.get_var_names(), vv);

	// save parameters to .par file
	output_file_writer.write_par(file_manager.open_ofile_ext("base.par"), pars, *(pts.get_offset_ptr()),
		*(pts.get_scale_ptr()));
	file_manager.close_file("par");

	// save new residuals to .rei file
	output_file_writer.write_rei(file_manager.open_ofile_ext("base.rei"), 0,
		pest_scenario.get_ctl_observations(), obs, obj_func, pars);

	// for models y(t+1) = g(x, y(t)), extract y(t+1) and use it as initial state for next time cycle
	//dyn_states = get_dynamic_states();
	//vector<string> real_names = _oe.get_real_names();
	Eigen::MatrixXd obs_i;
	if (obs_dyn_state_names.size() > 0) {
		message(1, "update initial dynamic states for next time cycle. Dynamic states  matrix size is", dyn_states.size());

		Eigen::MatrixXd mat = _oe.get_eigen(vector<string>(), obs_dyn_state_names);
		obs_i = mat.replicate(pe.shape().first, 1);
		pe.replace_col_vals(par_dyn_state_names, obs_i);

	}
	oe = _oe;
	return;
}


void DataAssimilator::initialize(int _icycle)
{
	icycle = _icycle;
	message(0, "initializing cycle", icycle);
	stringstream ss;
	ofstream& f_rec = file_manager.rec_ofstream();
	vector <string> dyn_states;

	da_ctl_params = pest_scenario.get_pestpp_options().da_ctl_params;
	pp_args = pest_scenario.get_pestpp_options().get_passed_args();
	act_obs_names = pest_scenario.get_ctl_ordered_nz_obs_names();
	act_par_names = pest_scenario.get_ctl_ordered_adj_par_names();

	// da type that is lower case is ok, I should move this somewhere else
	transform(da_type.begin(), da_type.end(), da_type.begin(), ::toupper);
	
	initialize_dynamic_states();

	// run the model one time using the initial parameters values. 
	if (pest_scenario.get_control_info().noptmax == 0)
	{
		forward_run_noptmax_0(_icycle);
		return;
	}
	
	//set some defaults
	PestppOptions* ppo = pest_scenario.get_pestpp_options_ptr();

	// TODO: (Ayman) We need to develop a da verbose level, for now ies is being used
	//verbose_level = pest_scenario.get_pestpp_options_ptr()->get_ies_verbose_level();
	verbose_level = da_ctl_params.get_ivalue("DA_VERBOSE_LEVEL");
	if (pest_scenario.get_n_adj_par() >= 1e6)
	{
		message(0, "Too many parameters ... !");
	}
	message(1, "using REDSVD for truncated svd solve");
	message(1, "maxsing:", pest_scenario.get_svd_info().maxsing);
	message(1, "eigthresh: ", pest_scenario.get_svd_info().eigthresh);

	message(1, "initializing localizer--");
	use_localizer = localizer.initialize(performance_log);	

	num_threads = da_ctl_params.get_ivalue("DA_NUM_THREADS");
	if (!use_localizer)
		message(1, "not using localization");
	else
	{
		if (localizer.get_autoadaloc())
		{
			message(1, "using automatic adaptive localization");
			message(2, "with autoadaloc_sigma_dist ", da_ctl_params.get_dvalue("DA_AUTOADALOC_SIGMA_DIST"));
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
	last_best_mean = 1.0E+30;
	last_best_std = 1.0e+30;
	lambda_max = 1.0E+30;
	lambda_min = 1.0E-30;
	warn_min_reals = 10;
	error_min_reals = 2;
	consec_bad_lambda_cycles = 0;

	
	// report DA paramters use
	lam_mults = da_ctl_params.get_vvalue("DA_INFLATION_MULT");
	infl_facs = da_ctl_params.get_vvalue("DA_INFLATION_FAC");
	

	sanity_checks();

	bool echo = false;
	if (verbose_level > 1)
		echo = true;

	// ies_use_empirical --- TODO (Ayman: check if we need to add if-statement for use only in ies
	initialize_parcov();
	initialize_obscov();

	
	subset_size = da_ctl_params.get_ivalue("DA_SUBSET_SIZE");
	reg_factor = pest_scenario.get_pestpp_options().get_ies_reg_factor(); //Todo: how this fit in da
	message(1, "using reg_factor: ", reg_factor);
	double bad_phi = pest_scenario.get_pestpp_options().get_ies_bad_phi();
	if (bad_phi < 1.0e+30)
		message(1, "using bad_phi: ", bad_phi);
	


	int num_reals = pest_scenario.get_pestpp_options().get_da_num_reals();

	
	pe_drawn = true;
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

	//this is being done globally now
	/*if (icycle == 0)
	{
		try
		{
			pe.check_for_dups();
		}
		catch (const exception & e)
		{
			string message = e.what();
			throw_em_error("error in parameter ensemble: " + message);
		}

	}*/

	try
	{
		oe.check_for_dups();
	}
	catch (const exception & e)
	{
		string message = e.what();
		throw_em_error("error in observation ensemble: " + message);
	}

	if ((act_obs_names.size() > 0) && (pe.shape().first != oe.shape().first))
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
				message(1, " the realization names are compatible");
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
				throw_em_error(ss.str());

			}
		}
		else
		{
			ss.str("");
			ss << "parameter ensemble rows (" << pe.shape().first << ") not equal to observation ensemble rows (" << oe.shape().first << ")";
			throw_em_error(ss.str());
		}
	}


	//need this here for Am calcs...
	//message(1, "transforming parameter ensemble to numeric");
	//if (icycle == 0)
		// parameters need to be transformed once ?
	pe.transform_ip(ParameterEnsemble::transStatus::NUM);

	if (da_ctl_params.get_bvalue("DA_ADD_BASE"))
	{
		if (pp_args.find("DA_RESTART_OBS_EN") != pp_args.end())  // todo: think about how this restart works in da?
		{
			message(1, "Warning: even though `da_include_base` is true, you passed a restart obs en, not adding 'base' realization...");
		}
		else
			add_bases();
	}

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
		frec << endl << "WARNING: the following " << misaligned.size() << " realization names are shared between the parameter and observation ensembles but they are not in the same row locations:" << endl;
		for (auto ma : misaligned)
			frec << ma << endl;
	}

	message(2, "checking for denormal values in pe");
	//if (icycle == 0)
	pe.check_for_normal("initial transformed parameter ensemble");

	//message(1, "saved initial parameter ensemble to ", ss.str());
	message(2, "checking for denormal values in base oe");
	oe.check_for_normal("base observation ensemble");
	ss.str("");
	if (pest_scenario.get_pestpp_options().get_save_binary())
	{
		ss << file_manager.get_base_filename() << ".base.obs.jcb";
		oe.to_binary(ss.str());
	}
	else
	{
		ss << file_manager.get_base_filename() << ".base.obs.csv";
		oe.to_csv(ss.str());
	}
	message(1, "saved base observation ensemble (obsval+noise) to ", ss.str());

	if (center_on.size() > 0)
	{
		ss.str("");
		ss << "centering on realization: '" << center_on << "' ";
		message(1, ss.str());
		vector<string> names = pe.get_real_names();
		if (find(names.begin(), names.end(), center_on) == names.end())
			throw_em_error("'ies_center_on' realization not found in par en: " + center_on);
		names = oe.get_real_names();
		if (find(names.begin(), names.end(), center_on) == names.end())
			throw_em_error("'ies_center_on' realization not found in obs en: " + center_on);
	}
	else
		message(1, "centering on ensemble mean vector");

	// Run the ensemble mean and move to the next cycle
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
			message(0, "mean parameter value run failed...bummer");
			return;
		}
		string obs_csv = file_manager.get_base_filename() + ".mean.obs.csv";
		message(1, "saving results from mean parameter value run to ", obs_csv);
		_oe.to_csv(obs_csv);

		ph.update(_oe, _pe);
		message(0, "mean parameter phi report:");
		ph.report();

		// for models y(t+1) = g(x, y(t)), extract y(t+1) and as initial state for next time cycle
		//dyn_states = get_dynamic_states();
		if (obs_dyn_state_names.size() != 0)
		{
			message(1, "update initial dynamic states for next time cycle. Dynamic states  matrix size is", dyn_states.size());
		}
		//add_dynamic_state_to_pe();
		//vector<string> real_names = _oe.get_real_names();
		Eigen::MatrixXd obs_i;

		if (obs_dyn_state_names.size() > 0) {
			Eigen::MatrixXd mat = _oe.get_eigen(vector<string>(), obs_dyn_state_names);
			obs_i = mat.replicate(pe.shape().first, 1);
			pe.replace_col_vals(par_dyn_state_names, obs_i);

		}

		return;
	}

	if (subset_size > pe.shape().first)
	{
		use_subset = false;
	}
	else
	{
		message(1, "using subset in lambda testing, number of realizations used in subset testing: ", subset_size);
		string how = da_ctl_params.get_svalue("DA_SUBSET_HOW");
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

	// ==================================================
	// ==================== forward run =================
	// ==================================================

	if (obs_restart_csv.size() > 0)
		//the hard way to restart
		initialize_restart();

	//no restart
	else
	{
		performance_log->log_event("running ensemble for data assimilation cycle No." + icycle);
		message(1, "runing ensemble of size", pe.shape().first);

		vector<int> failed = run_ensemble(pe, oe, vector<int>(), icycle);

		if (pe.shape().first == 0)
			throw_em_error("all realizations failed during initial evaluation");

		if (failed.size() > 0)
			oe_base.drop_rows(failed);

		//if (icycle == 0)
		//{
			pe.transform_ip(ParameterEnsemble::transStatus::NUM);
		//}
	}

	// ======================================================================
	// extract dynamic states forecast from oe and add them to pe.
	// for a model y(t+1) = g(x,y(t)), this code use simulated y(t+1) as input for the next time
	// cycle

	//States are flaged as ones that have zero weights and exist in both oe and pe

	//dyn_states = get_dynamic_states();
	
	if (obs_dyn_state_names.size() != 0)
	{
		message(1, "add dynamic states to forecast ensemble. Dynamic states  matrix size is ", dyn_states.size());
	}
	add_dynamic_state_to_pe();
	// ==

	ss.str("");
	da_save_ensemble_pe("_forecast_cycle_", ".par");
	da_save_ensemble_oe("_forecast_cycle_", ".obs");

	performance_log->log_event("calc pre-drop phi");
	//initialize the phi handler
	ss.str("");
	ss << icycle << ".";
	ph = L2PhiHandler(&pest_scenario, &file_manager, &oe_base, &pe_base, &parcov, true, ss.str());

	if (ph.get_lt_obs_names().size() > 0) // todo: (Ayman) how to handle those obs in DA? 
	{
		message(1, "less_than inequality defined for observations: ", ph.get_lt_obs_names().size());
	}
	if (ph.get_gt_obs_names().size())
	{
		message(1, "greater_than inequality defined for observations: ", ph.get_gt_obs_names().size());
	}

	ph.update(oe, pe);
	message(0, "pre-drop initial phi summary");
	ph.report();
	drop_bad_phi(pe, oe);
	if (oe.shape().first == 0)
	{
		throw_em_error(string("all realizations dropped as 'bad'"));
	}
	if (oe.shape().first <= error_min_reals)
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


	pcs = ParChangeSummarizer(&pe_base, &file_manager, &output_file_writer);
	vector<string> in_conflict = detect_prior_data_conflict();
	if (in_conflict.size() > 0)
	{
		ss.str("");
		ss << "WARNING: " << in_conflict.size() << " non-zero weighted observations are in conflict";
		ss << " with the prior simulated ensemble." << endl;
		message(0, ss.str());
	}

	cout << "...see rec file for listing of conflicted observations" << endl << endl;
	ofstream& frec = file_manager.rec_ofstream();
	frec << endl << "...conflicted observations: " << endl;
	for (auto oname : in_conflict)
	{
		frec << oname << endl;
	}
	if (!ppo->get_ies_drop_conflicts())
	{
		ss.str("");
		ss << "  Continuing with data assimilation will likely result ";
		ss << " in parameter bias and, ultimately, forecast bias";
		message(1, ss.str());
	}
	else
	{

		//check that all obs are in conflict
		message(1, "dropping conflicted observations");
		if (in_conflict.size() == oe.shape().second)
		{
			throw_em_error("all non-zero weighted observations in conflict state, cannot continue");
		}
		//drop from act_obs_names
		vector<string> t;
		set<string> sconflict(in_conflict.begin(), in_conflict.end());
		for (auto oname : act_obs_names)
			if (sconflict.find(oname) == sconflict.end())
				t.push_back(oname);
		act_obs_names = t;

		//update obscov
		obscov.drop(in_conflict);

		//drop from oe_base
		oe_base.drop_cols(in_conflict);
		//shouldnt need to update localizer since we dropping not adding
		//updating weights in control file

		ObservationInfo* oi = pest_scenario.get_observation_info_ptr();
		int org_nnz_obs = pest_scenario.get_ctl_ordered_nz_obs_names().size();
		for (auto n : in_conflict)
		{
			oi->set_weight(n, 0.0);
		}

		stringstream ss;
		ss << "number of non-zero weighted observations reduced from " << org_nnz_obs;
		ss << " to " << pest_scenario.get_ctl_ordered_nz_obs_names().size() << endl;
		message(1, ss.str());
	}
	performance_log->log_event("calc initial phi");
	ph.update(oe, pe);
	message(0, "initial phi summary");
	ph.report();
	ph.write(0, run_mgr_ptr->get_total_runs());
	best_mean_phis.push_back(ph.get_mean(L2PhiHandler::phiType::COMPOSITE));
	if (!pest_scenario.get_pestpp_options().get_ies_use_approx())
	{
		message(1, "using full (MAP) update solution");

	}

	if (ph.get_mean(L2PhiHandler::phiType::ACTUAL) < 1.0e-10)
	{
		throw_em_error("initial actual phi mean too low, something is wrong...or you have the perfect model that already fits the data shockingly well");
	}
	if (ph.get_std(L2PhiHandler::phiType::ACTUAL) < 1.0e-10)
	{
		throw_em_error("initial actual phi stdev too low, something is wrong...");
	}

	last_best_mean = ph.get_mean(L2PhiHandler::phiType::COMPOSITE);
	last_best_std = ph.get_std(L2PhiHandler::phiType::COMPOSITE);
	last_best_lam = da_ctl_params.get_dvalue("DA_INITIAL_INF_FAC");

	if (last_best_mean < 1.0e-10)
	{
		throw_em_error("initial composite phi mean too low, something is wrong...");
	}
	if (last_best_std < 1.0e-10)
	{
		throw_em_error("initial composite phi stdev too low, something is wrong...");
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

void DataAssimilator:: da_save_ensemble_pe(string fprefix, string dtyp)
{
	// todo: update for da
	stringstream ss;
	ss.str("");
	bool save_binary;
	save_binary = da_ctl_params.get_bvalue("DA_SAVE_BINARY");
	
	if (save_binary)
	{
		ss << file_manager.get_base_filename() << fprefix << icycle << dtyp <<".jcb";
		pe.to_binary(ss.str());
	}
	else
	{
		ss << file_manager.get_base_filename() << fprefix << icycle << dtyp <<".csv";
		pe.to_csv(ss.str());
	}
	//message(1, "saved forecast (parameter and state) ensemble to ", ss.str());
}
void DataAssimilator::da_save_ensemble_oe(string fprefix, string dtyp)
{
	// todo: update for da
	stringstream ss;
	ss.str("");
	bool save_binary;
	save_binary = da_ctl_params.get_bvalue("DA_SAVE_BINARY");
	if (save_binary)
	{
		ss << file_manager.get_base_filename() << fprefix << icycle << dtyp << ".jcb";
		oe.to_binary(ss.str());
	}
	else
	{
		ss << file_manager.get_base_filename() << fprefix << icycle << dtyp << ".csv";
		oe.to_csv(ss.str());
	}
	//message(1, "saved forecast (parameter and state) ensemble to ", ss.str());
}
void DataAssimilator::add_dynamic_state_to_pe()
{
	
	//vector<string> real_names = oe.get_real_names();
	ParameterEnsemble::transStatus org_status = pe.get_trans_status();
	ParamTransformSeq bts = pest_scenario.get_base_par_tran_seq();
	if (obs_dyn_state_names.size() > 0)	{	
			
		Eigen::MatrixXd mat = oe.get_eigen(vector<string>(), obs_dyn_state_names);
		if (org_status == ParameterEnsemble::transStatus::NUM)
		{
			for (int i = 0; i < mat.rows(); i++)
			{
				Parameters pars = pest_scenario.get_ctl_parameters();
				
				pars.update(par_dyn_state_names, mat.row(i));
				bts.ctl2numeric_ip(pars);
				mat.row(i) = pars.get_data_eigen_vec(par_dyn_state_names);
			}
		}
		pe.replace_col_vals(par_dyn_state_names, mat);
		
	}
	
}
void DataAssimilator::initialize_dynamic_states()
{
	stringstream ss;
	// find a list of dynamic states
	//vector <string> dyn_states_names;	
	obs_dyn_state_names.clear();
	par_dyn_state_names.clear();
	vector<string> obs_names = oe.get_var_names(); // todo: get obs names from a different source. 
	vector<string> par_names = pe.get_var_names();
	set<string> spar_names(par_names.begin(), par_names.end());
	set<string>::iterator end = spar_names.end();
	par_names.clear();
	//for (int i = 0; i < obs_names.size(); i++)	
	double w;
	for (auto& name : obs_names)
	{
		
		w = pest_scenario.get_observation_info_ptr()->get_weight(name);		
		//if ((w == 0) && (spar_names.find(name) != end))
		if (spar_names.find(name) != end)
		{
			obs_dyn_state_names.push_back(name);
			par_dyn_state_names.push_back(name);
		}

	}
	if (obs_dyn_state_names.size() > 0)
	{
		ss.str("");
		ss << obs_dyn_state_names.size() << " non-zero weighted dynamic states identified through shared-names";
		message(1, ss.str());
	}
	map<string, string> state_map = pest_scenario.get_ext_file_string_map("observation data external", "state_par_link");
	if (state_map.size() > 0)
	{
		
		set<string> pstates(par_dyn_state_names.begin(), par_dyn_state_names.end());
		set<string>::iterator send = pstates.end();
		vector<string> t = pest_scenario.get_ctl_ordered_par_names();
		set<string> pnames(t.begin(), t.end());
		t.clear();
		set<string>::iterator pend = pnames.end();
		vector<string> dups,missing;
		int c = 0;
		for (auto& sm : state_map)
		{
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
				w = pest_scenario.get_observation_info_ptr()->get_weight(sm.first);
				//if (w == 0)
				{
					obs_dyn_state_names.push_back(sm.first);
					par_dyn_state_names.push_back(sm.second);
					c++;
				}
			}
		}
		if (dups.size() > 0)
		{
			stringstream ss;
			ss << "the following state parameters nominated thru obs data linking " << endl;
			ss << "    were already tagged as 'states' by identically named observations:" << endl;
			for (auto& d : dups)
				ss << d << ",";		
			throw_em_error(ss.str());
		}
		if (missing.size() > 0)
		{
			stringstream ss;
			ss << "the following parameters nominated thru obs data linking " << endl;
			ss << "    were not found in par data section:" << endl;
			for (auto& m : missing)
				ss << m << ",";
			throw_em_error(ss.str());
		}
		if (c > 0)
		{
			ss.str("");
			ss << c << " non-zero weighted dynamic states identified through 'state_par_link'";
			message(1, ss.str());
		}
	}
}

vector<string> DataAssimilator::detect_prior_data_conflict()
{
	message(1, "checking for prior-data conflict...");
	//for now, just really simple metric - checking for overlap
	vector<string> in_conflict;
	double smin, smax, omin, omax;
	map<string, int> smap, omap;
	vector<string> snames = oe.get_var_names();
	vector<string> onames = oe_base.get_var_names();

	for (int i = 0; i < snames.size(); i++)
	{
		smap[snames[i]] = i;
	}
	for (int i = 0; i < onames.size(); i++)
	{
		omap[onames[i]] = i;
	}
	int sidx, oidx; // Ayman: this is to check if the model predictions does not envelop obs values. 
				    // If yes, we still can run KF, but it means that our prior is not wide enough or the model is biased!
	for (auto oname : pest_scenario.get_ctl_ordered_nz_obs_names())
	{
		sidx = smap[oname];
		oidx = omap[oname];
		smin = oe.get_eigen_ptr()->col(sidx).minCoeff();
		omin = oe_base.get_eigen_ptr()->col(oidx).minCoeff();
		smax = oe.get_eigen_ptr()->col(sidx).maxCoeff();
		omax = oe_base.get_eigen_ptr()->col(oidx).maxCoeff();
		if ((smin > omax) || (smax < omin))
			in_conflict.push_back(oname);
	}
	return in_conflict;
}

//Eigen::MatrixXd DataAssimilator::get_Am(const vector<string>& real_names, const vector<string>& par_names)
//{
//
//	double scale = (1.0 / (sqrt(double(real_names.size() - 1))));
//	Eigen::MatrixXd par_diff = scale * pe_base.get_eigen_anomalies(real_names, par_names, pest_scenario.get_pestpp_options().get_ies_center_on());
//	par_diff.transposeInPlace();
//	if (verbose_level > 1)
//	{
//		cout << "prior_par_diff: " << par_diff.rows() << ',' << par_diff.cols() << endl;
//		if (verbose_level > 2)
//			save_mat("prior_par_diff.dat", par_diff);
//	}
//
//	Eigen::MatrixXd ivec, upgrade_1, s, V, U, st;
//	SVD_REDSVD rsvd;
//	//SVD_EIGEN rsvd;
//	rsvd.set_performance_log(performance_log);
//
//	rsvd.solve_ip(par_diff, s, U, V, pest_scenario.get_svd_info().eigthresh, pest_scenario.get_svd_info().maxsing);
//	par_diff.resize(0, 0);
//	Eigen::MatrixXd temp = s.asDiagonal();
//	Eigen::MatrixXd Am = U * temp;
//	return Am;
//}

//void DataAssimilator::drop_bad_phi(ParameterEnsemble& _pe, ObservationEnsemble& _oe, bool is_subset)
//{
//	//don't use this assert because _pe maybe full size, but _oe might be subset size
//	if (!is_subset)
//		if (_pe.shape().first != _oe.shape().first)
//			throw_em_error("DataAssimilator::drop_bad_phi() error: _pe != _oe and not subset");
//
//	double bad_phi = pest_scenario.get_pestpp_options().get_ies_bad_phi();
//	double bad_phi_sigma = pest_scenario.get_pestpp_options().get_ies_bad_phi_sigma();
//	vector<int> idxs = ph.get_idxs_greater_than(bad_phi, bad_phi_sigma, _oe);
//
//	if (pest_scenario.get_pestpp_options().get_ies_debug_bad_phi())
//		idxs.push_back(0);
//
//	if (idxs.size() > 0)
//	{
//
//		message(0, "dropping realizations as bad: ", idxs.size());
//
//		vector<string> par_real_names = _pe.get_real_names(), obs_real_names = _oe.get_real_names();
//		stringstream ss;
//		string pname;
//		string oname;
//
//		int pidx;
//		vector<string> full_onames, full_pnames;
//		// if a subset drop, then use the full oe index, otherwise, just use _oe index
//		/*if (_oe.shape().first != _pe.shape().first)
//		{
//			full_onames = oe.get_real_names();
//		}
//		else
//		{
//			full_onames = _oe.get_real_names();
//		}*/
//		full_onames = oe.get_real_names();
//		full_pnames = pe.get_real_names();
//		vector<string> pdrop, odrop;
//		for (auto i : idxs)
//		{
//			oname = obs_real_names[i];
//
//			if (is_subset)
//			{
//				pidx = find(full_onames.begin(), full_onames.end(), oname) - full_onames.begin();
//				if (find(subset_idxs.begin(), subset_idxs.end(), pidx) == subset_idxs.end())
//				{
//					ss.str("");
//					ss << "drop_bad_phi() error: idx " << pidx << " not found in subset_idxs";
//					throw_em_error(ss.str());
//				}
//				pname = full_pnames[pidx];
//			}
//			else
//			{
//				pidx = i;
//				pname = par_real_names[pidx];
//			}
//			ss << pname << " : " << obs_real_names[i] << " , ";
//			pdrop.push_back(pname);
//			odrop.push_back(obs_real_names[i]);
//		}
//
//		string s = "dropping par:obs realizations: " + ss.str();
//		message(1, s);
//		try
//		{
//			_pe.drop_rows(pdrop);
//			_oe.drop_rows(odrop);
//		}
//		catch (const exception & e)
//		{
//			stringstream ss;
//			ss << "drop_bad_phi() error : " << e.what();
//			throw_em_error(ss.str());
//		}
//		catch (...)
//		{
//			throw_em_error(string("drop_bad_phi() error"));
//		}
//	}
//}

//void DataAssimilator::save_mat(string prefix, Eigen::MatrixXd& mat)
//{
//	stringstream ss;
//	ss << iter << '.' << prefix;
//	try
//	{
//		ofstream& f = file_manager.open_ofile_ext(ss.str());
//		f << mat << endl;
//		f.close();
//		file_manager.close_file(ss.str());
//	}
//	catch (...)
//	{
//		message(1, "error saving matrix", ss.str());
//	}
//}

void DataAssimilator::da_upate()
{
	stringstream ss;
	ofstream& frec = file_manager.rec_ofstream();

	string da_method = da_ctl_params.get_svalue("DA_TYPE");
	message(0, "Assimilation Method :", da_method);
	int user_noptmax = pest_scenario.get_control_info().noptmax;
	
	if (da_type == "MDA")
	{
		// make sure that noptmax and number of inflations are consistent. always noptmax will be respected. 
		if (user_noptmax > infl_facs.size())
		{
			int beta_diff = user_noptmax - infl_facs.size();
			for (int ibeta = 0; ibeta < beta_diff; ibeta++)
			{
				infl_facs.push_back(infl_facs.back());
			}
		}
		else if (user_noptmax < infl_facs.size())
		{
			int beta_diff = infl_facs.size() - user_noptmax;
			vector<double> v2 = std::vector<double>(infl_facs.begin(), infl_facs.end() - beta_diff);
			infl_facs = v2;
		}
		solution_iterations = infl_facs.size();
	}
	else if (da_type == "VANILLA")
		solution_iterations = 1;
	else
		solution_iterations = pest_scenario.get_control_info().noptmax;
	bool accept;
	for (int i = 0; i < solution_iterations; i++)
	{
		iter++;
		message(0, "starting solve for iteration:", iter);
		ss << "starting solve for iteration: " << iter;
		performance_log->log_event(ss.str());
		accept = solve_new_da();
		report_and_save();
		ph.update(oe, pe);
		last_best_mean = ph.get_mean(L2PhiHandler::phiType::COMPOSITE);
		last_best_std = ph.get_std(L2PhiHandler::phiType::COMPOSITE);
		ph.report(true);
		ph.write(iter, run_mgr_ptr->get_total_runs());
		if (pest_scenario.get_pestpp_options().get_ies_save_rescov())
			ph.save_residual_cov(oe, iter);
		pcs.summarize(pe, iter);


		if (accept)
			consec_bad_lambda_cycles = 0;
		else
			consec_bad_lambda_cycles++;

		if (should_terminate())		
			break;

	}
	//todo return posterior here
	da_save_ensemble_pe("_analysis_cycle_", ".par");
	da_save_ensemble_oe("_analysis_cycle_", ".obs");
	if (obs_dyn_state_names.size() > 0)
		update_starting_state();


}

void DataAssimilator::kf_upate()
{
	stringstream ss;
	ofstream& frec = file_manager.rec_ofstream();

	bool accept;
	for (int i = 0; i < pest_scenario.get_control_info().noptmax; i++)
	{
		iter++;
		message(0, "starting solve for iteration:", iter);
		ss << "starting solve for iteration: " << iter;
		performance_log->log_event(ss.str());
		accept = solve_new_da();
		report_and_save();
		ph.update(oe, pe);
		last_best_mean = ph.get_mean(L2PhiHandler::phiType::COMPOSITE);
		last_best_std = ph.get_std(L2PhiHandler::phiType::COMPOSITE);
		ph.report(true);
		ph.write(iter, run_mgr_ptr->get_total_runs());
		if (pest_scenario.get_pestpp_options().get_ies_save_rescov())
			ph.save_residual_cov(oe, iter);
		pcs.summarize(pe, iter);


		if (accept)
			consec_bad_lambda_cycles = 0;
		else
			consec_bad_lambda_cycles++;

		if (should_terminate())
			break;
	}
	
}

//bool DataAssimilator::should_terminate()
//{
//	//todo: use ies accept fac here?
//	double phiredstp = pest_scenario.get_control_info().phiredstp;
//	int nphistp = pest_scenario.get_control_info().nphistp;
//	int nphinored = pest_scenario.get_control_info().nphinored;
//	bool phiredstp_sat = false, nphinored_sat = false, consec_sat = false;
//	double phi, ratio;
//	int count = 0;
//	int nphired = 0;
//	//best_mean_phis = vector<double>{ 1.0,0.8,0.81,0.755,1.1,0.75,0.75,1.2 };
//
//
//
//	/*if ((!consec_sat )&& (best_mean_phis.size() == 0))
//		return false;*/
//	message(0, "phi-based termination criteria check");
//	message(1, "phiredstp: ", phiredstp);
//	message(1, "nphistp: ", nphistp);
//	message(1, "nphinored (also used for consecutive bad lambda cycles): ", nphinored);
//	if (best_mean_phis.size() > 0)
//	{
//		vector<double>::iterator idx = min_element(best_mean_phis.begin(), best_mean_phis.end());
//		nphired = (best_mean_phis.end() - idx) - 1;
//		best_phi_yet = best_mean_phis[idx - best_mean_phis.begin()];// *pest_scenario.get_pestpp_options().get_ies_accept_phi_fac();
//		message(1, "best mean phi sequence: ", best_mean_phis);
//		message(1, "best phi yet: ", best_phi_yet);
//	}
//	message(1, "number of consecutive bad lambda testing cycles: ", consec_bad_lambda_cycles);
//	if (consec_bad_lambda_cycles >= nphinored)
//	{
//		message(1, "number of consecutive bad lambda testing cycles > nphinored");
//		consec_sat = true;
//	}
//
//	for (auto& phi : best_mean_phis)
//	{
//		ratio = (phi - best_phi_yet) / phi;
//		if (ratio <= phiredstp)
//			count++;
//	}
//	message(1, "number of iterations satisfying phiredstp criteria: ", count);
//	if (count >= nphistp)
//	{
//		message(1, "number iterations satisfying phiredstp criteria > nphistp");
//		phiredstp_sat = true;
//	}
//
//	message(1, "number of iterations since best yet mean phi: ", nphired);
//	if (nphired >= nphinored)
//	{
//		message(1, "number of iterations since best yet mean phi > nphinored");
//		nphinored_sat = true;
//	}
//
//	if ((nphinored_sat) || (phiredstp_sat) || (consec_sat))
//	{
//		message(1, "phi-based termination criteria satisfied, all done");
//		return true;
//	}
//	return false;
//}
//



LocalUpgradeThread_da::LocalUpgradeThread_da(PerformanceLog* _performance_log, unordered_map<string, Eigen::VectorXd>& _par_resid_map, unordered_map<string, Eigen::VectorXd>& _par_diff_map,
	unordered_map<string, Eigen::VectorXd>& _obs_resid_map, unordered_map<string, Eigen::VectorXd>& _obs_diff_map,
	Localizer& _localizer, unordered_map<string, double>& _parcov_inv_map, unordered_map<string, double>& _weight_map,
	ParameterEnsemble& _pe_upgrade, unordered_map<string, pair<vector<string>, vector<string>>>& _cases,
	unordered_map<string, Eigen::VectorXd>& _Am_map, Localizer::How& _how, unordered_map<string, Eigen::VectorXd>& _obs_err_map) : par_resid_map(_par_resid_map),
	par_diff_map(_par_diff_map), obs_resid_map(_obs_resid_map), obs_diff_map(_obs_diff_map), localizer(_localizer),
	pe_upgrade(_pe_upgrade), cases(_cases), parcov_inv_map(_parcov_inv_map), weight_map(_weight_map), Am_map(_Am_map), obs_err_map(_obs_err_map)
{
	performance_log = _performance_log;
	how = _how;
	parcov_inv_map = _parcov_inv_map;
	weight_map = _weight_map;
	//obs_err_map = _obs_err_map;
	count = 0;

	for (auto& c : cases)
	{
		keys.push_back(c.first);
	}
	//sort(keys.begin(), keys.end());
	total = keys.size();
	//random_shuffle(keys.begin(), keys.end());

}


//void LocalUpgradeThread_da::work(int thread_id, int iter, double cur_lam)
//{
//	class local_utils
//	{
//	public:
//		static Eigen::DiagonalMatrix<double, Eigen::Dynamic> get_matrix_from_map(vector<string>& names, unordered_map<string, double>& dmap)
//		{
//			Eigen::VectorXd vec(names.size());
//			int i = 0;
//			for (auto name : names)
//			{
//				vec[i] = dmap.at(name);
//				++i;
//			}
//			Eigen::DiagonalMatrix<double, Eigen::Dynamic> m = vec.asDiagonal();
//			return m;
//		}
//		static Eigen::MatrixXd get_matrix_from_map(int num_reals, vector<string>& names, unordered_map<string, Eigen::VectorXd>& emap)
//		{
//			Eigen::MatrixXd mat(num_reals, names.size());
//			mat.setZero();
//
//			for (int j = 0; j < names.size(); j++)
//			{
//				mat.col(j) = emap[names[j]];
//			}
//
//			return mat;
//		}
//		static void save_mat(int verbose_level, int tid, int iter, int t_count, string prefix, Eigen::MatrixXd& mat)
//		{
//			if (verbose_level < 2)
//				return;
//
//			if (verbose_level < 3)
//				return;
//			//cout << "thread: " << tid << ", " << t_count << ", " << prefix << " rows:cols" << mat.rows() << ":" << mat.cols() << endl;
//			stringstream ss;
//
//			ss << "thread_" << tid << ".count_ " << t_count << ".iter_" << iter << "." << prefix << ".dat";
//			string fname = ss.str();
//			ofstream f(fname);
//			if (!f.good())
//				cout << "error getting ofstream " << fname << endl;
//			else
//			{
//
//				try
//				{
//					f << mat << endl;
//					f.close();
//				}
//				catch (...)
//				{
//					cout << "error saving matrix " << fname << endl;
//				}
//			}
//		}
//	};
//
//	stringstream ss;
//
//
//	unique_lock<mutex> ctrl_guard(ctrl_lock, defer_lock);
//	int maxsing, num_reals, verbose_level, pcount = 0;
//	int t_count;
//	double eigthresh;
//	bool use_approx;
//	bool use_prior_scaling;
//	bool use_localizer = false;
//	bool loc_by_obs = true;
//
//	while (true)
//	{
//		if (ctrl_guard.try_lock())
//		{
//			maxsing = pe_upgrade.get_pest_scenario_ptr()->get_svd_info().maxsing;
//			eigthresh = pe_upgrade.get_pest_scenario_ptr()->get_svd_info().eigthresh;
//			use_approx = pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_ies_use_approx();
//			use_prior_scaling = pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_ies_use_prior_scaling();
//			num_reals = pe_upgrade.shape().first;
//			verbose_level = pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_ies_verbose_level();
//			ctrl_guard.unlock();
//			//if (pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_ies_localize_how()[0] == 'P')
//			if (how == Localizer::How::PARAMETERS)
//				loc_by_obs = false;
//			break;
//		}
//	}
//	ofstream f_thread;
//	if (verbose_level > 2)
//	{
//		ss.str("");
//		ss << "thread_" << thread_id << "part_map.csv";
//		f_thread.open(ss.str());
//		ss.str("");
//	}
//	Eigen::MatrixXd par_resid, par_diff, Am;
//	Eigen::MatrixXd obs_resid, obs_diff, loc, obs_err;
//	Eigen::DiagonalMatrix<double, Eigen::Dynamic> weights, parcov_inv;
//	vector<string> par_names, obs_names;
//	while (true)
//	{
//		unique_lock<mutex> next_guard(next_lock, defer_lock);
//		par_names.clear();
//		obs_names.clear();
//		use_localizer = false;
//		//the end condition
//		//the end condition
//		while (true)
//		{
//			if (next_guard.try_lock())
//			{
//				if (count == keys.size())
//				{
//					if (verbose_level > 1)
//					{
//						cout << "upgrade thread: " << thread_id << " processed " << pcount << " upgrade parts" << endl;
//					}
//					if (f_thread.good())
//						f_thread.close();
//					return;
//				}
//				string k = keys[count];
//				pair<vector<string>, vector<string>> p = cases.at(k);
//				par_names = p.second;
//				obs_names = p.first;
//				if (localizer.get_use())
//				{
//					if ((loc_by_obs) && (par_names.size() == 1) && (k == par_names[0]))
//						use_localizer = true;
//					else if ((!loc_by_obs) && (obs_names.size() == 1) && (k == obs_names[0]))
//					{
//						use_localizer = true;
//						//loc_by_obs = false;
//					}
//				}
//				if (count % 1000 == 0)
//				{
//					ss.str("");
//					ss << "upgrade thread progress: " << count << " of " << total << " parts done";
//					if (verbose_level > 1)
//						cout << ss.str() << endl;
//					performance_log->log_event(ss.str());
//				}
//				count++;
//				t_count = count;
//				pcount++;
//				next_guard.unlock();
//				break;
//			}
//		}
//
//
//		if (verbose_level > 2)
//		{
//			f_thread << t_count << "," << iter;
//			for (auto name : par_names)
//				f_thread << "," << name;
//			for (auto name : obs_names)
//				f_thread << "," << name;
//			f_thread << endl;
//		}
//
//		par_resid.resize(0, 0);
//		par_diff.resize(0, 0);
//		obs_resid.resize(0, 0);
//		obs_diff.resize(0, 0);
//		loc.resize(0, 0);
//		obs_err.resize(0, 0);
//		Am.resize(0, 0);
//		weights.resize(0);
//		parcov_inv.resize(0);
//		Am.resize(0, 0);
//
//		unique_lock<mutex> obs_diff_guard(obs_diff_lock, defer_lock);
//		unique_lock<mutex> obs_resid_guard(obs_resid_lock, defer_lock);
//		unique_lock<mutex> obs_err_guard(obs_err_lock, defer_lock);
//		unique_lock<mutex> par_diff_guard(par_diff_lock, defer_lock);
//		unique_lock<mutex> par_resid_guard(par_resid_lock, defer_lock);
//		unique_lock<mutex> loc_guard(loc_lock, defer_lock);
//		unique_lock<mutex> weight_guard(weight_lock, defer_lock);
//		unique_lock<mutex> parcov_guard(parcov_lock, defer_lock);
//		unique_lock<mutex> am_guard(am_lock, defer_lock);
//
//		while (true)
//		{
//			if (((use_approx) || (par_resid.rows() > 0)) &&
//				(weights.size() > 0) &&
//				(parcov_inv.size() > 0) &&
//				(par_diff.rows() > 0) &&
//				(obs_resid.rows() > 0) &&
//				(obs_diff.rows() > 0) &&
//				((!use_localizer) || (loc.rows() > 0)) &&
//				((use_approx) || (Am.rows() > 0)))
//				break;
//			if ((use_localizer) && (loc.rows() == 0) && (loc_guard.try_lock()))
//			{
//				if (loc_by_obs)
//					loc = localizer.get_localizing_par_hadamard_matrix(num_reals, obs_names[0], par_names);
//				else
//					loc = localizer.get_localizing_obs_hadamard_matrix(num_reals, par_names[0], obs_names);
//				loc_guard.unlock();
//			}
//			if ((obs_diff.rows() == 0) && (obs_diff_guard.try_lock()))
//			{
//				//piggy back here for thread safety
//				//if (pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_svd_pack() == PestppOptions::SVD_PACK::PROPACK)
//				//	use_propack = true;
//				obs_diff = local_utils::get_matrix_from_map(num_reals, obs_names, obs_diff_map);
//				obs_diff_guard.unlock();
//			}
//			if ((obs_resid.rows() == 0) && (obs_resid_guard.try_lock()))
//			{
//				obs_resid = local_utils::get_matrix_from_map(num_reals, obs_names, obs_resid_map);
//				obs_resid_guard.unlock();
//			}
//			if ((obs_err.rows() == 0) && (obs_err_guard.try_lock()))
//			{
//				obs_err = local_utils::get_matrix_from_map(num_reals, obs_names, obs_err_map);
//				obs_err_guard.unlock();
//			}
//			if ((par_diff.rows() == 0) && (par_diff_guard.try_lock()))
//			{
//				par_diff = local_utils::get_matrix_from_map(num_reals, par_names, par_diff_map);
//				par_diff_guard.unlock();
//			}
//			if ((par_resid.rows() == 0) && (par_resid_guard.try_lock()))
//			{
//				par_resid = local_utils::get_matrix_from_map(num_reals, par_names, par_resid_map);
//				par_resid_guard.unlock();
//			}
//			if ((weights.rows() == 0) && (weight_guard.try_lock()))
//			{
//				weights = local_utils::get_matrix_from_map(obs_names, weight_map);
//				weight_guard.unlock();
//			}
//			if ((parcov_inv.rows() == 0) && (parcov_guard.try_lock()))
//			{
//				parcov_inv = local_utils::get_matrix_from_map(par_names, parcov_inv_map);
//				parcov_guard.unlock();
//			}
//			if ((!use_approx) && (Am.rows() == 0) && (am_guard.try_lock()))
//			{
//				//Am = local_utils::get_matrix_from_map(num_reals, par_names, Am_map).transpose();
//				int am_cols = Am_map[par_names[0]].size();
//				Am.resize(par_names.size(), am_cols);
//				Am.setZero();
//
//				for (int j = 0; j < par_names.size(); j++)
//				{
//					Am.row(j) = Am_map[par_names[j]];
//				}
//				am_guard.unlock();
//			}
//		}
//
//		par_diff.transposeInPlace();
//		obs_diff.transposeInPlace();
//		obs_resid.transposeInPlace();
//		obs_err.transposeInPlace();
//		par_resid.transposeInPlace();
//
//		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "obs_resid", obs_resid);
//		Eigen::MatrixXd scaled_residual = weights * obs_resid;
//
//
//
//		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "par_resid", par_resid);
//		Eigen::MatrixXd scaled_par_resid;
//		if ((!use_approx) && (iter > 1))
//		{
//			if (use_prior_scaling)
//			{
//				scaled_par_resid = parcov_inv * par_resid;
//			}
//			else
//			{
//				scaled_par_resid = par_resid;
//			}
//		}
//
//		stringstream ss;
//
//		double scale = (1.0 / (sqrt(double(num_reals - 1))));
//		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "obs_diff", obs_diff);
//
//		if (use_localizer)
//			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "loc", loc);
//		if (use_localizer)
//		{
//			if (loc_by_obs)
//				par_diff = par_diff.cwiseProduct(loc);
//			else
//				obs_diff = obs_diff.cwiseProduct(loc);
//
//		}
//
//		obs_diff = scale * (weights * obs_diff);
//		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "par_diff", par_diff);
//		if (use_prior_scaling)
//			par_diff = scale * parcov_inv * par_diff;
//		else
//			par_diff = scale * par_diff;
//
//
//		//performance_log->log_event("SVD of obs diff");
//		Eigen::MatrixXd ivec, upgrade_1, s, V, Ut;
//
//
//		SVD_REDSVD rsvd;
//		rsvd.solve_ip(obs_diff, s, Ut, V, eigthresh, maxsing);
//
//		Ut.transposeInPlace();
//		obs_diff.resize(0, 0);
//		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "Ut", Ut);
//		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "s", s);
//		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "V", V);
//
//		Eigen::MatrixXd s2 = s.cwiseProduct(s);
//
//		ivec = ((Eigen::VectorXd::Ones(s2.size()) * cur_lam) + s2).asDiagonal().inverse();
//		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "ivec", ivec);
//		Eigen::MatrixXd X1 = Ut * scaled_residual;
//		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X1", X1);
//		Eigen::MatrixXd X2 = ivec * X1;
//		X1.resize(0, 0);
//
//		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X2", X2);
//		Eigen::MatrixXd X3 = V * s.asDiagonal() * X2;
//		X2.resize(0, 0);
//
//		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X3", X3);
//		if (use_prior_scaling)
//		{
//			upgrade_1 = -1.0 * parcov_inv * par_diff * X3;
//		}
//		else
//		{
//			upgrade_1 = -1.0 * par_diff * X3;
//		}
//		upgrade_1.transposeInPlace();
//		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "upgrade_1", upgrade_1);
//		X3.resize(0, 0);
//
//
//		Eigen::MatrixXd upgrade_2;
//		if ((!use_approx) && (iter > 1))
//		{
//			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "Am", Am);
//			Eigen::MatrixXd x4 = Am.transpose() * scaled_par_resid;
//			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X4", x4);
//
//			par_resid.resize(0, 0);
//
//			Eigen::MatrixXd x5 = Am * x4;
//			x4.resize(0, 0);
//			Am.resize(0, 0);
//
//			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X5", x5);
//			Eigen::MatrixXd x6 = par_diff.transpose() * x5;
//			x5.resize(0, 0);
//
//			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X6", x6);
//			Eigen::MatrixXd x7 = V * ivec * V.transpose() * x6;
//			x6.resize(0, 0);
//
//			if (use_prior_scaling)
//			{
//				upgrade_2 = -1.0 * parcov_inv * par_diff * x7;
//			}
//			else
//			{
//				upgrade_2 = -1.0 * (par_diff * x7);
//			}
//			x7.resize(0, 0);
//
//			upgrade_1 = upgrade_1 + upgrade_2.transpose();
//			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "upgrade_2", upgrade_2);
//			upgrade_2.resize(0, 0);
//
//		}
//
//		unique_lock<mutex> put_guard(put_lock, defer_lock);
//		while (true)
//		{
//			if (put_guard.try_lock())
//			{
//				pe_upgrade.add_2_cols_ip(par_names, upgrade_1);
//				put_guard.unlock();
//				break;
//			}
//		}
//	}
//
//}
//
//void LocalUpgradeThread_da::work2da(int thread_id, int iter, double cur_lam)
//{
//	class local_utils
//	{
//	public:
//		static Eigen::DiagonalMatrix<double, Eigen::Dynamic> get_matrix_from_map(vector<string>& names, unordered_map<string, double>& dmap)
//		{
//			Eigen::VectorXd vec(names.size());
//			int i = 0;
//			for (auto name : names)
//			{
//				vec[i] = dmap.at(name);
//				++i;
//			}
//			Eigen::DiagonalMatrix<double, Eigen::Dynamic> m = vec.asDiagonal();
//			return m;
//		}
//		static Eigen::MatrixXd get_matrix_from_map(int num_reals, vector<string>& names, unordered_map<string, Eigen::VectorXd>& emap)
//		{
//			Eigen::MatrixXd mat(num_reals, names.size());
//			mat.setZero();
//
//			for (int j = 0; j < names.size(); j++)
//			{
//				mat.col(j) = emap[names[j]];
//			}
//
//			return mat;
//		}
//		static void save_mat(int verbose_level, int tid, int iter, int t_count, string prefix, Eigen::MatrixXd& mat)
//		{
//			if (verbose_level < 2)
//				return;
//
//			if (verbose_level < 3)
//				return;
//			//cout << "thread: " << tid << ", " << t_count << ", " << prefix << " rows:cols" << mat.rows() << ":" << mat.cols() << endl;
//			stringstream ss;
//
//			ss << "thread_" << tid << ".count_ " << t_count << ".iter_" << iter << "." << prefix << ".dat";
//			string fname = ss.str();
//			ofstream f(fname);
//			if (!f.good())
//				cout << "error getting ofstream " << fname << endl;
//			else
//			{
//
//				try
//				{
//					f << mat << endl;
//					f.close();
//				}
//				catch (...)
//				{
//					cout << "error saving matrix " << fname << endl;
//				}
//			}
//		}
//	};
//
//	stringstream ss;
//
//
//	unique_lock<mutex> ctrl_guard(ctrl_lock, defer_lock);
//	int maxsing, num_reals, verbose_level, pcount = 0;
//	int t_count;
//	double eigthresh;
//	bool use_approx;
//	bool use_prior_scaling;
//	bool use_localizer = false;
//	bool loc_by_obs = true;
//
//	while (true)
//	{
//		if (ctrl_guard.try_lock())
//		{
//			maxsing = pe_upgrade.get_pest_scenario_ptr()->get_svd_info().maxsing;
//			eigthresh = pe_upgrade.get_pest_scenario_ptr()->get_svd_info().eigthresh;
//			use_approx = pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_ies_use_approx();
//			use_prior_scaling = pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_ies_use_prior_scaling();
//			num_reals = pe_upgrade.shape().first;
//			verbose_level = pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_ies_verbose_level();
//			ctrl_guard.unlock();
//			//if (pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_ies_localize_how()[0] == 'P')
//			if (how == Localizer::How::PARAMETERS)
//				loc_by_obs = false;
//			break;
//		}
//	}
//	ofstream f_thread;
//	if (verbose_level > 2)
//	{
//		ss.str("");
//		ss << "thread_" << thread_id << "part_map.csv";
//		f_thread.open(ss.str());
//		ss.str("");
//	}
//	Eigen::MatrixXd par_resid, par_diff, Am;
//	Eigen::MatrixXd obs_resid, obs_diff, loc, obs_err;
//	Eigen::DiagonalMatrix<double, Eigen::Dynamic> weights, parcov_inv;
//	vector<string> par_names, obs_names;
//	while (true)
//	{
//		unique_lock<mutex> next_guard(next_lock, defer_lock);
//		par_names.clear();
//		obs_names.clear();
//		use_localizer = false;
//		//the end condition
//		//the end condition
//		while (true)
//		{
//			if (next_guard.try_lock())
//			{
//				if (count == keys.size())
//				{
//					if (verbose_level > 1)
//					{
//						cout << "upgrade thread: " << thread_id << " processed " << pcount << " upgrade parts" << endl;
//					}
//					if (f_thread.good())
//						f_thread.close();
//					return;
//				}
//				string k = keys[count];
//				pair<vector<string>, vector<string>> p = cases.at(k);
//				par_names = p.second;
//				obs_names = p.first;
//				if (localizer.get_use())
//				{
//					if ((loc_by_obs) && (par_names.size() == 1) && (k == par_names[0]))
//						use_localizer = true;
//					else if ((!loc_by_obs) && (obs_names.size() == 1) && (k == obs_names[0]))
//					{
//						use_localizer = true;
//						//loc_by_obs = false;
//					}
//				}
//				if (count % 1000 == 0)
//				{
//					ss.str("");
//					ss << "upgrade thread progress: " << count << " of " << total << " parts done";
//					if (verbose_level > 1)
//						cout << ss.str() << endl;
//					performance_log->log_event(ss.str());
//				}
//				count++;
//				t_count = count;
//				pcount++;
//				next_guard.unlock();
//				break;
//			}
//		}
//
//
//		if (verbose_level > 2)
//		{
//			f_thread << t_count << "," << iter;
//			for (auto name : par_names)
//				f_thread << "," << name;
//			for (auto name : obs_names)
//				f_thread << "," << name;
//			f_thread << endl;
//		}
//
//		par_resid.resize(0, 0);
//		par_diff.resize(0, 0);
//		obs_resid.resize(0, 0);
//		obs_diff.resize(0, 0);
//		loc.resize(0, 0);
//		obs_err.resize(0, 0);
//		Am.resize(0, 0);
//		weights.resize(0);
//		parcov_inv.resize(0);
//		Am.resize(0, 0);
//
//		unique_lock<mutex> obs_diff_guard(obs_diff_lock, defer_lock);
//		unique_lock<mutex> obs_resid_guard(obs_resid_lock, defer_lock);
//		unique_lock<mutex> obs_err_guard(obs_err_lock, defer_lock);
//		unique_lock<mutex> par_diff_guard(par_diff_lock, defer_lock);
//		unique_lock<mutex> par_resid_guard(par_resid_lock, defer_lock);
//		unique_lock<mutex> loc_guard(loc_lock, defer_lock);
//		unique_lock<mutex> weight_guard(weight_lock, defer_lock);
//		unique_lock<mutex> parcov_guard(parcov_lock, defer_lock);
//		unique_lock<mutex> am_guard(am_lock, defer_lock);
//
//		while (true)
//		{
//			if (((use_approx) || (par_resid.rows() > 0)) &&
//				(weights.size() > 0) &&
//				(parcov_inv.size() > 0) &&
//				(par_diff.rows() > 0) &&
//				(obs_resid.rows() > 0) &&
//				(obs_diff.rows() > 0) &&
//				((!use_localizer) || (loc.rows() > 0)) &&
//				((use_approx) || (Am.rows() > 0)))
//				break;
//			if ((use_localizer) && (loc.rows() == 0) && (loc_guard.try_lock()))
//			{
//				if (loc_by_obs)
//					loc = localizer.get_localizing_par_hadamard_matrix(num_reals, obs_names[0], par_names);
//				else
//					loc = localizer.get_localizing_obs_hadamard_matrix(num_reals, par_names[0], obs_names);
//				loc_guard.unlock();
//			}
//			if ((obs_diff.rows() == 0) && (obs_diff_guard.try_lock()))
//			{
//				//piggy back here for thread safety
//				//if (pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_svd_pack() == PestppOptions::SVD_PACK::PROPACK)
//				//	use_propack = true;
//				obs_diff = local_utils::get_matrix_from_map(num_reals, obs_names, obs_diff_map);
//				obs_diff_guard.unlock();
//			}
//			if ((obs_resid.rows() == 0) && (obs_resid_guard.try_lock()))
//			{
//				obs_resid = local_utils::get_matrix_from_map(num_reals, obs_names, obs_resid_map);
//				obs_resid_guard.unlock();
//			}
//			if ((obs_err.rows() == 0) && (obs_err_guard.try_lock()))
//			{
//				obs_err = local_utils::get_matrix_from_map(num_reals, obs_names, obs_err_map);
//				obs_err_guard.unlock();
//			}
//			if ((par_diff.rows() == 0) && (par_diff_guard.try_lock()))
//			{
//				par_diff = local_utils::get_matrix_from_map(num_reals, par_names, par_diff_map);
//				par_diff_guard.unlock();
//			}
//			if ((par_resid.rows() == 0) && (par_resid_guard.try_lock()))
//			{
//				par_resid = local_utils::get_matrix_from_map(num_reals, par_names, par_resid_map);
//				par_resid_guard.unlock();
//			}
//			if ((weights.rows() == 0) && (weight_guard.try_lock()))
//			{
//				weights = local_utils::get_matrix_from_map(obs_names, weight_map);
//				weight_guard.unlock();
//			}
//			if ((parcov_inv.rows() == 0) && (parcov_guard.try_lock()))
//			{
//				parcov_inv = local_utils::get_matrix_from_map(par_names, parcov_inv_map);
//				parcov_guard.unlock();
//			}
//			if ((!use_approx) && (Am.rows() == 0) && (am_guard.try_lock()))
//			{
//				//Am = local_utils::get_matrix_from_map(num_reals, par_names, Am_map).transpose();
//				int am_cols = Am_map[par_names[0]].size();
//				Am.resize(par_names.size(), am_cols);
//				Am.setZero();
//
//				for (int j = 0; j < par_names.size(); j++)
//				{
//					Am.row(j) = Am_map[par_names[j]];
//				}
//				am_guard.unlock();
//			}
//		}
//		Eigen::MatrixXd upgrade_1;
//		par_diff.transposeInPlace();
//		obs_diff.transposeInPlace();
//		obs_resid.transposeInPlace();
//		obs_err.transposeInPlace();
//		par_resid.transposeInPlace();
//
//		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "obs_resid", obs_resid);
//		Eigen::MatrixXd scaled_residual = weights * obs_resid;
//
//		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "par_resid", par_resid);
//		Eigen::MatrixXd scaled_par_resid;
//		if ((!use_approx) && (iter > 1))
//		{
//			if (use_prior_scaling)
//			{
//				scaled_par_resid = parcov_inv * par_resid;
//			}
//			else
//			{
//				scaled_par_resid = par_resid;
//			}
//		}
//
//		stringstream ss;
//
//		double scale = (1.0 / (sqrt(double(num_reals - 1))));
//		local_utils::save_mat(verbose_level, thread_id, iter, t_count, "obs_diff", obs_diff);
//
//		if (use_localizer)
//			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "loc", loc);
//		if (use_localizer)
//		{
//			if (loc_by_obs)
//				par_diff = par_diff.cwiseProduct(loc);
//			else
//				obs_diff = obs_diff.cwiseProduct(loc);
//
//		}
//
//		/// ****************************************
//		obs_diff = scale * obs_diff;// (H-Hm)/sqrt(N-1)		
//		par_diff = scale * par_diff;// (K-Km)/sqrt(N-1)		
//		obs_err = scale * obs_err; //  (E-Em)/sqrt(N-1)	
//		obs_resid = (-1.0) * obs_resid; // D-H  , the negative one is becuase obs_resid is calculated as H-D ???
//
//		Eigen::MatrixXd s2_, s, V, U, cum_sum;
//		SVD_REDSVD rsvd;
//		Eigen::MatrixXd C;
//		C = obs_diff + (cur_lam * obs_err); // curr_lam is the inflation factor 		
//		Eigen::VectorXd s2;
//		
//		if (false) // check this 
//		{
//			rsvd.solve_ip(C, s, U, V, eigthresh, maxsing);
//		}
//		else
//		{
//			
//			Eigen::JacobiSVD<Eigen::MatrixXd> svd(C, Eigen::ComputeThinU | Eigen::ComputeThinV);
//			U = svd.matrixU();
//			s = svd.singularValues();
//
//			int ssize = s.size();
//			double threshold_frac;
//			s2 = s.cwiseProduct(s);
//			double s_sum = s.sum();
//			cum_sum = s;
//			for (int i = 0; i < ssize; i++)
//			{
//
//				threshold_frac = s(i) / s_sum;
//				if (threshold_frac < eigthresh)
//					s2(i) = 0;
//			}
//
//		}
//		s2_ = s2.asDiagonal().inverse();
//		for (int i = 0; i < s2.size(); i++)
//		{
//			if (s2(i) < 1e-50)
//			{
//				s2_(i, i) = 0;
//			}
//		}
//		Eigen::MatrixXd X1 = s2_ * U.transpose();
//		
//		X1 = X1 * obs_resid;
//		
//		X1 = U * X1;
//		
//		X1 = obs_diff.transpose() * X1;
//		
//		upgrade_1 = par_diff * X1;
//		upgrade_1.transposeInPlace();
//				
//
//		////******************************************
//		if (0)
//		{
//			if (true) //Evensen's solution
//			{
//				obs_diff = obs_diff + cur_lam * obs_err;
//			}
//			else
//			{
//				// do nothing
//			}
//
//			obs_diff = scale * (weights * obs_diff);
//			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "par_diff", par_diff);
//			if (use_prior_scaling)
//				par_diff = scale * parcov_inv * par_diff;
//			else
//				par_diff = scale * par_diff;
//
//
//			//performance_log->log_event("SVD of obs diff");
//			Eigen::MatrixXd ivec, upgrade_1, s, V, Ut;
//
//
//			SVD_REDSVD rsvd;
//			rsvd.solve_ip(obs_diff, s, Ut, V, eigthresh, maxsing);
//
//			Ut.transposeInPlace();
//			obs_diff.resize(0, 0);// this is becuase obs_diff = USVt
//			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "Ut", Ut);
//			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "s", s);
//			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "V", V);
//
//			Eigen::MatrixXd s2 = s.cwiseProduct(s);
//
//			if (true) // Evensen's solution
//			{
//				ivec = s2.asDiagonal().inverse();
//			}
//			else
//			{
//				ivec = ((Eigen::VectorXd::Ones(s2.size()) * cur_lam) + s2).asDiagonal().inverse();
//			}
//
//			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "ivec", ivec);
//			Eigen::MatrixXd X1 = Ut * scaled_residual;
//			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X1", X1);
//			Eigen::MatrixXd X2 = ivec * X1;
//			X1.resize(0, 0);
//
//			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X2", X2);
//			Eigen::MatrixXd X3 = V * s.asDiagonal() * X2;
//			X2.resize(0, 0);
//
//			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X3", X3);
//			if (use_prior_scaling)
//			{
//				upgrade_1 = -1.0 * parcov_inv * par_diff * X3;
//			}
//			else
//			{
//				upgrade_1 = -1.0 * par_diff * X3;
//			}
//			upgrade_1.transposeInPlace();
//			local_utils::save_mat(verbose_level, thread_id, iter, t_count, "upgrade_1", upgrade_1);
//			X3.resize(0, 0);
//
//
//			Eigen::MatrixXd upgrade_2;
//			if ((!use_approx) && (iter > 1))
//			{
//				local_utils::save_mat(verbose_level, thread_id, iter, t_count, "Am", Am);
//				Eigen::MatrixXd x4 = Am.transpose() * scaled_par_resid;
//				local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X4", x4);
//
//				par_resid.resize(0, 0);
//
//				Eigen::MatrixXd x5 = Am * x4;
//				x4.resize(0, 0);
//				Am.resize(0, 0);
//
//				local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X5", x5);
//				Eigen::MatrixXd x6 = par_diff.transpose() * x5;
//				x5.resize(0, 0);
//
//				local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X6", x6);
//				Eigen::MatrixXd x7 = V * ivec * V.transpose() * x6;
//				x6.resize(0, 0);
//
//				if (use_prior_scaling)
//				{
//					upgrade_2 = -1.0 * parcov_inv * par_diff * x7;
//				}
//				else
//				{
//					upgrade_2 = -1.0 * (par_diff * x7);
//				}
//				x7.resize(0, 0);
//
//				upgrade_1 = upgrade_1 + upgrade_2.transpose();
//				local_utils::save_mat(verbose_level, thread_id, iter, t_count, "upgrade_2", upgrade_2);
//				upgrade_2.resize(0, 0);
//
//			}
//		}
//
//		unique_lock<mutex> put_guard(put_lock, defer_lock);
//		while (true)
//		{
//			if (put_guard.try_lock())
//			{
//				pe_upgrade.add_2_cols_ip(par_names, upgrade_1);
//				put_guard.unlock();
//				break;
//			}
//		}
//	}
//
//}


//void upgrade_thread_function(int id, int iter, double cur_lam, LocalUpgradeThread_da& worker, exception_ptr& eptr)
//{
//	try
//	{
//		worker.work(id, iter, cur_lam);
//	}
//	catch (...)
//	{
//		eptr = current_exception();
//	}
//
//	return;
//}
//
//
//void upgrade_thread_function_da(int id, int iter, double cur_lam, LocalUpgradeThread_da& worker, exception_ptr& eptr)
//{
//	try
//	{
//		worker.work2da(id, iter, cur_lam);
//	}
//	catch (...)
//	{
//		eptr = current_exception();
//	}
//
//	return;
//}


//ParameterEnsemble DataAssimilator::calc_localized_upgrade_threaded(double cur_lam, unordered_map<string, pair<vector<string>, vector<string>>>& loc_map)
//{
//	stringstream ss;
//
//	ObservationEnsemble oe_upgrade(oe.get_pest_scenario_ptr(), &rand_gen, oe.get_eigen(vector<string>(), act_obs_names, false), oe.get_real_names(), act_obs_names);
//	ParameterEnsemble pe_upgrade(pe.get_pest_scenario_ptr(), &rand_gen, pe.get_eigen(vector<string>(), act_par_names, false), pe.get_real_names(), act_par_names);
//
//	//this copy of the localizer map will be consumed by the worker threads
//	//unordered_map<string, pair<vector<string>, vector<string>>> loc_map;
//	if (use_localizer)
//	{
//
//		//loc_map = localizer.get_localizer_map(iter, oe, pe, performance_log);
//		//localizer.report(file_manager.rec_ofstream());
//	}
//	else
//	{
//		pair<vector<string>, vector<string>> p(act_obs_names, act_par_names);
//		loc_map["all"] = p;
//	}
//
//	//prep the fast look par cov info
//	message(2, "preparing fast-look containers for threaded localization solve");
//	unordered_map<string, double> parcov_inv_map;
//	parcov_inv_map.reserve(pe_upgrade.shape().second);
//	Eigen::VectorXd parcov_inv;// = parcov.get(par_names).inv().e_ptr()->toDense().cwiseSqrt().asDiagonal();
//	if (!parcov.isdiagonal())
//	{
//		//parcov_inv = parcov.inv().get_matrix().diagonal().cwiseSqrt();
//		parcov_inv = parcov.get_matrix().diagonal();
//
//	}
//	else
//	{
//		message(2, "extracting diagonal from prior parameter covariance matrix");
//		Covariance parcov_diag;
//		parcov_diag.from_diagonal(parcov);
//		//parcov_inv = parcov_diag.inv().get_matrix().diagonal().cwiseSqrt();
//		parcov_inv = parcov_diag.get_matrix().diagonal();
//
//	}
//	parcov_inv = parcov_inv.cwiseSqrt().cwiseInverse();
//
//	vector<string> par_names = pe_upgrade.get_var_names();
//	for (int i = 0; i < parcov_inv.size(); i++)
//		parcov_inv_map[par_names[i]] = parcov_inv[i];
//
//
//	//prep the fast lookup weights info
//	unordered_map<string, double> weight_map;
//	weight_map.reserve(oe_upgrade.shape().second);
//	vector<string> obs_names = oe_upgrade.get_var_names();
//
//	Eigen::MatrixXd meas_errors = oe_base.get_eigen(oe_base.get_real_names(), obs_names);
//	unordered_map<string, Eigen::VectorXd> obs_err_map;
//	obs_err_map.reserve(obs_names.size());
//	Observations obs_value = pest_scenario.get_ctl_observations();
//
//	double w;
//	for (int i = 0; i < obs_names.size(); i++)
//	{
//		w = pest_scenario.get_observation_info_ptr()->get_weight(obs_names[i]);
//		//don't want to filter on weight here - might be changing weights, etc...
//		//if (w == 0.0)
//		//	continue;
//		weight_map[obs_names[i]] = w;
//
//		meas_errors.col(i).array() = meas_errors.col(i).array() - obs_value.get_rec(obs_names[i]);// ;
//		obs_err_map[obs_names[i]] = meas_errors.col(i);
//	}
//
//	//prep a fast look for obs, par and resid matrices - stored as column vectors in a map
//	unordered_map<string, Eigen::VectorXd> par_resid_map, obs_resid_map, Am_map;
//	unordered_map<string, Eigen::VectorXd> par_diff_map, obs_diff_map;
//	par_resid_map.reserve(par_names.size());
//	par_diff_map.reserve(par_names.size());
//	Am_map.reserve(par_names.size());
//	obs_resid_map.reserve(obs_names.size());
//	obs_diff_map.reserve(obs_names.size());
//
//	//check for the 'center_on' real - it may have been dropped...
//	string center_on = pest_scenario.get_pestpp_options().get_ies_center_on();
//	if (center_on.size() > 0)
//	{
//		vector<string> real_names = oe_upgrade.get_real_names();
//		if (find(real_names.begin(), real_names.end(), center_on) == real_names.end())
//		{
//			message(0, "Warning: 'ies_center_on' real not found in obs en, reverting to mean...");
//			center_on = "";
//		}
//		real_names = pe_upgrade.get_real_names();
//		if ((center_on.size() > 0) && find(real_names.begin(), real_names.end(), center_on) == real_names.end())
//		{
//			message(0, "Warning: 'ies_center_on' real not found in par en, reverting to mean...");
//			center_on = "";
//		}
//	}
//
//	Eigen::MatrixXd mat = ph.get_obs_resid_subset(oe_upgrade);
//	for (int i = 0; i < obs_names.size(); i++)
//	{
//		obs_resid_map[obs_names[i]] = mat.col(i);
//	}
//	mat = oe_upgrade.get_eigen_anomalies(center_on);
//	for (int i = 0; i < obs_names.size(); i++)
//	{
//		obs_diff_map[obs_names[i]] = mat.col(i);
//	}
//	mat = ph.get_par_resid_subset(pe_upgrade);
//	for (int i = 0; i < par_names.size(); i++)
//	{
//		par_resid_map[par_names[i]] = mat.col(i);
//	}
//	mat = pe_upgrade.get_eigen_anomalies(center_on);
//	for (int i = 0; i < par_names.size(); i++)
//	{
//		par_diff_map[par_names[i]] = mat.col(i);
//	}
//	if (!pest_scenario.get_pestpp_options().get_ies_use_approx())
//	{
//		mat = get_Am(pe_upgrade.get_real_names(), pe_upgrade.get_var_names());
//		for (int i = 0; i < par_names.size(); i++)
//		{
//			Am_map[par_names[i]] = mat.row(i);
//		}
//	}
//	mat.resize(0, 0);
//	// clear the upgrade ensemble
//	pe_upgrade.set_zeros();
//	Localizer::How _how = localizer.get_how();
//	LocalUpgradeThread_da worker(performance_log, par_resid_map, par_diff_map, obs_resid_map, obs_diff_map,
//		localizer, parcov_inv_map, weight_map, pe_upgrade, loc_map, Am_map, _how, obs_err_map);
//
//	if ((num_threads < 1) || (loc_map.size() == 1))
//		//if (num_threads < 1)
//	{
//		if (true) // use da 
//		{
//			worker.work2da(0, iter, cur_lam);
//		}
//		else
//		{
//			worker.work(0, iter, cur_lam);
//		}
//	}
//	else
//	{
//		Eigen::setNbThreads(1);
//		vector<thread> threads;
//		vector<exception_ptr> exception_ptrs;
//		message(2, "launching threads");
//
//		for (int i = 0; i < num_threads; i++)
//		{
//			exception_ptrs.push_back(exception_ptr());
//		}
//
//		for (int i = 0; i < num_threads; i++)
//		{
//			//threads.push_back(thread(&LocalUpgradeThread::work, &worker, i, iter, cur_lam));
//			if (true) // use da
//			{
//				threads.push_back(thread(upgrade_thread_function_da, i, iter, cur_lam, std::ref(worker), std::ref(exception_ptrs[i])));
//
//			}
//			else
//			{
//				threads.push_back(thread(upgrade_thread_function, i, iter, cur_lam, std::ref(worker), std::ref(exception_ptrs[i])));
//			}
//		}
//		message(2, "waiting to join threads");
//		//for (auto &t : threads)
//		//	t.join();
//		ss.str("");
//		int num_exp = 0;
//
//		for (int i = 0; i < num_threads; ++i)
//		{
//			bool found = false;
//			if (exception_ptrs[i])
//			{
//				found = true;
//				num_exp++;
//				try
//				{
//					rethrow_exception(exception_ptrs[i]);
//				}
//				catch (const std::exception & e)
//				{
//					//ss.str("");
//					ss << " thread " << i << "raised an exception: " << e.what();
//					//throw runtime_error(ss.str());
//				}
//				catch (...)
//				{
//					//ss.str("");
//					ss << " thread " << i << "raised an exception";
//					//throw runtime_error(ss.str());
//				}
//			}
//			threads[i].join();
//			if ((exception_ptrs[i]) && (!found))
//			{
//				num_exp++;
//				try
//				{
//					rethrow_exception(exception_ptrs[i]);
//				}
//				catch (const std::exception & e)
//				{
//					//ss.str("");
//					ss << " thread " << i << "raised an exception: " << e.what();
//					//throw runtime_error(ss.str());
//				}
//				catch (...)
//				{
//					//ss.str("");
//					ss << " thread " << i << "raised an exception: ";
//					//throw runtime_error(ss.str());
//				}
//			}
//		}
//		if (num_exp > 0)
//		{
//			throw runtime_error(ss.str());
//		}
//		message(2, "threaded localized upgrade calculation done");
//	}
//
//	return pe_upgrade;
//}



//ParameterEnsemble DataAssimilator::calc_localized_upgrade(double cur_lam)
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

//void DataAssimilator::update_reals_by_phi(ParameterEnsemble& _pe, ObservationEnsemble& _oe)
//{
//
//	vector<string> oe_names = _oe.get_real_names();
//	vector<string> pe_names = _pe.get_real_names();
//	vector<string> oe_base_names = oe.get_real_names();
//	vector<string> pe_base_names = pe.get_real_names();
//
//	//if (pe_names.size() != oe_base_names.size())
//	//	throw runtime_error("DataAssimilator::update_reals_by_phi() error: pe_names != oe_base_names");
//	map<string, int> oe_name_to_idx;
//	map<int, string> pe_idx_to_name;
//
//	for (int i = 0; i < oe_base_names.size(); i++)
//		oe_name_to_idx[oe_base_names[i]] = i;
//
//	for (int i = 0; i < pe_base_names.size(); i++)
//		pe_idx_to_name[i] = pe_base_names[i];
//	//store map of current phi values
//	ph.update(oe, pe);
//	L2PhiHandler::phiType pt = L2PhiHandler::phiType::COMPOSITE;
//	map<string, double>* phi_map = ph.get_phi_map_ptr(pt);
//	map<string, double> cur_phi_map;
//	for (auto p : *phi_map)
//		cur_phi_map[p.first] = p.second;
//
//	//now get a phi map of the new phi values
//	ph.update(_oe, _pe);
//	phi_map = ph.get_phi_map_ptr(pt);
//
//	double acc_fac = pest_scenario.get_pestpp_options().get_ies_accept_phi_fac();
//	double cur_phi, new_phi;
//	string oname, pname;
//	Eigen::VectorXd real;
//	stringstream ss;
//	for (int i = 0; i < _oe.shape().first; i++)
//	{
//		oname = oe_names[i];
//		new_phi = phi_map->at(oname);
//		cur_phi = cur_phi_map.at(oname);
//		if (new_phi < cur_phi * acc_fac)
//		{
//			//pname = pe_names[i];
//			//pname = pe_names[oe_name_to_idx[oname]];
//			pname = pe_idx_to_name[oe_name_to_idx[oname]];
//			if (find(pe_names.begin(), pe_names.end(), pname) == pe_names.end())
//				throw runtime_error("DataAssimilator::update_reals_by_phi() error: pname not in pe_names: " + pname);
//			ss.str("");
//			ss << "updating pe:oe real =" << pname << ":" << oname << ", current phi: new phi  =" << cur_phi << ":" << new_phi;
//			message(3, ss.str());
//
//			real = _pe.get_real_vector(pname);
//			pe.update_real_ip(pname, real);
//			real = _oe.get_real_vector(oname);
//			oe.update_real_ip(oname, real);
//		}
//	}
//	ph.update(oe, pe);
//
//}

ParameterEnsemble DataAssimilator::calc_kf_upgrade(double cur_lam, unordered_map<string, pair<vector<string>, vector<string>>>& loc_map)
{
	stringstream ss;

	ObservationEnsemble oe_upgrade(oe.get_pest_scenario_ptr(), &rand_gen, oe.get_eigen(vector<string>(), act_obs_names, false), oe.get_real_names(), act_obs_names);
	ParameterEnsemble pe_upgrade(pe.get_pest_scenario_ptr(), &rand_gen, pe.get_eigen(vector<string>(), act_par_names, false), pe.get_real_names(), act_par_names);

	//this copy of the localizer map will be consumed by the worker threads
	//unordered_map<string, pair<vector<string>, vector<string>>> loc_map;
	if (use_localizer)
	{

		//loc_map = localizer.get_localizer_map(iter, oe, pe, performance_log);
		//localizer.report(file_manager.rec_ofstream());
	}
	else
	{
		// Ayman noticed that p is empty when no localizer is used
		//pair<vector<string>, vector<string>> p(act_obs_names, act_par_names);
		//loc_map["all"] = p;
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
	Observations obs_value = pest_scenario.get_ctl_observations();
	vector<string> names = oe.get_var_names();	
	Eigen::MatrixXd meas_errors = oe_base.get_eigen(oe_base.get_real_names(), obs_names);// -obs_value.get_rec(obs_names);
	unordered_map<string, Eigen::VectorXd> obs_err_map;
	obs_err_map.reserve(obs_names.size());
	
	double w;
	for (int i = 0; i < obs_names.size(); i++)
	{
		w = pest_scenario.get_observation_info_ptr()->get_weight(obs_names[i]);
		meas_errors.col(i).array() = meas_errors.col(i).array() - obs_value.get_rec(obs_names[i]);// ;
		obs_err_map[obs_names[i]] = meas_errors.col(i);		
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

	Eigen::MatrixXd mat = ph.get_obs_resid_subset(oe_upgrade); // oe - oe_base --->(H-D)
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
	for (int i = 0; i < par_names.size(); i++)
	{
		par_resid_map[par_names[i]] = mat.col(i);
	}
	mat = pe_upgrade.get_eigen_anomalies(center_on);
	for (int i = 0; i < par_names.size(); i++)
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



	kf_work(performance_log, par_resid_map, par_diff_map, obs_resid_map, obs_diff_map, obs_err_map,
		localizer, parcov_inv_map, weight_map, pe_upgrade, cur_lam, loc_map, Am_map, _how);
	/*
	
	*/

	return pe_upgrade;
}



bool DataAssimilator::solve_new_da()
{
	stringstream ss;
	ofstream& frec = file_manager.rec_ofstream();
	if (pe.shape().first <= error_min_reals)
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

	 // check subset size	
	if ((use_subset) && (subset_size > pe.shape().first))
	{
		ss.str("");
		ss << "++da_subset size (" << subset_size << ") greater than ensemble size (" << pe.shape().first << ")";
		frec << "  ---  " << ss.str() << endl;
		cout << "  ---  " << ss.str() << endl;
		frec << "  ...reducing subset_size to " << pe.shape().first << endl;
		cout << "  ...reducing subset_size to " << pe.shape().first << endl;
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
		loc_map = localizer.get_localanalysis_case_map(iter, oe, pe, performance_log);
		//localizer.report(file_manager.rec_ofstream());
	}
	else
	{
		pair<vector<string>, vector<string>> p(act_obs_names, act_par_names);
		loc_map["all"] = p;
	}

	
	if ((da_type == "MDA") && ((iter == 1)))
	{   
		// Multiple data assimilation (MDA) only. Make sure that sum(1/infl_fac) = 1		
		double inv_infl_sum = 0;
		for (auto& inf_fac : infl_facs)
		{
			inv_infl_sum = inv_infl_sum + (1.0 / inf_fac);
		}
		if (abs(inv_infl_sum - 1.0) < 1e-5)
		{
			// ok close enough to 1, maybe msg later 
		}
		else
		{			
			message(0, "WARNING:  Summation of reciprocals of inflation factors should be 1. Now it is ", inv_infl_sum);
			message(1, "Inflation factors are rescaled to correct this issue...", " ");
		}
		int fac_ii = 0;
		//ss.str("");
		for (auto& inf_fac : infl_facs)
		{			 
			infl_facs[fac_ii] = infl_facs[fac_ii] * inv_infl_sum;
			fac_ii = fac_ii + 1;		
		}
	}

	// from now on, beta's has different meaning in different contexts
	vector<double> betas;
	betas.clear();
	if (da_type == "MDA") // beta is the inflation factor at each cycle
	{		
		betas.push_back(infl_facs[iter-1]);
	}
	else if (da_type == "ITERATIVE")
	{
		betas = lam_mults;
	}
	else if (da_type == "VANILLA")
	{
		betas.push_back(1.0);
	}

	Eigen::MatrixXd Am;
	if (!pest_scenario.get_pestpp_options().get_ies_use_approx())
	{
		Am = get_Am(pe.get_real_names(), pe.get_var_names());
	}

	// =========== test all lam and scales =======
	for (auto& beta : betas) 
		// NOTE : beta is inflation factor when da_type = MDA and is inflation factor 
		//multiplier whehn da_type = ITERATIVE, and is 1 for Vanilla
	{
		ss.str("");
		double cur_lam;

		//ParameterEnsemble pe_upgrade;
		bool use_mda = false;
		if (da_type == "MDA")
		{
			cur_lam = beta;
			use_mda = true;
		}
		else if (da_type == "ITERATIVE")
		{
			cur_lam = 1+beta * last_best_lam;
		}
		else if (da_type == "VANILLA")
		{
			cur_lam = beta;
		}
		ss << "starting calcs for inflation factors" << cur_lam;
		message(1, "starting inflation calcs for inflation factor ", cur_lam);
		message(2, "see .log file for more details");

		// Here is the update .....
		/*if (0)
		{
			pe_upgrade = calc_kf_upgrade(cur_lam, loc_map);
		}
		else
		{
			pe_upgrade = calc_localized_upgrade_threaded(cur_lam, loc_map);
		}*/
		ParameterEnsemble pe_upgrade(pe.get_pest_scenario_ptr(), &rand_gen, pe.get_eigen(vector<string>(), act_par_names, false), pe.get_real_names(), act_par_names);
		pe_upgrade.set_trans_status(pe.get_trans_status());
		ObservationEnsemble oe_upgrade(oe.get_pest_scenario_ptr(), &rand_gen, oe.get_eigen(vector<string>(), act_obs_names, false), oe.get_real_names(), act_obs_names);

		EnsembleSolver es(performance_log, file_manager, pest_scenario,pe_upgrade,oe_upgrade,oe_base,localizer,parcov,
			Am,ph,use_localizer,iter,act_par_names,act_obs_names);
		es.solve(num_threads, cur_lam, !use_mda, pe_upgrade, loc_map);


		vector<double> scale_facs;
		scale_facs = da_ctl_params.get_vvalue("DA_SCALE_FAC");
		
		if ((da_type == "MDA") || (da_type == "VANILLA"))
		{
			// mda makes no changes to scaling factors; maybe in the future, users  should choose this
			scale_facs.clear();
			scale_facs.push_back(1.0);
		}

		for (auto sf : scale_facs)
		{

			ParameterEnsemble pe_lam_scale = pe;
			pe_lam_scale.set_eigen(*pe_lam_scale.get_eigen_ptr() + (*pe_upgrade.get_eigen_ptr() * sf));
			//todo: what chglim enforcement?
			if (da_ctl_params.get_bvalue("DA_ENFORCE_BOUNDS"))
			{				
				pe_lam_scale.enforce_bounds(performance_log, false);
			}

			pe_lams.push_back(pe_lam_scale); // This is a list of all pe's computed using different lam and scale
			lam_vals.push_back(cur_lam);
			scale_vals.push_back(sf);
			if (!pest_scenario.get_pestpp_options().get_ies_save_lambda_en()) // todo: add da specific option for this
				continue;
			ss.str("");
			ss << file_manager.get_base_filename() << "." << iter << "." << cur_lam << ".lambda." << sf << ".scale.par";
			
			if (pest_scenario.get_pestpp_options().get_save_binary())
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

		}// end scales loop

		ss.str("");
		message(1, "finished calcs for lambda:", cur_lam);

	}// end lams loop

	if (pest_scenario.get_pestpp_options().get_ies_debug_upgrade_only())
	{
		message(0, "ies_debug_upgrade_only is true, exiting");
		throw_em_error("ies_debug_upgrade_only is true, exiting");
	}

	//  =============== forward runs to test lams and scales ==================

	vector<map<int, int>> real_run_ids_lams;
	int best_idx = -1;
	double best_mean = 1.0e+30, best_std = 1.0e+30;
	double mean, std;
	vector<ObservationEnsemble> oe_lams;
	vector<ParameterEnsemble> posterior_dyn_states;
	bool evaluate_update = da_ctl_params.get_bvalue("DA_EVALUATE_UPDATE");

	if (da_type == "ITERATIVE")
	{
		evaluate_update = true;
	}

	if (da_type == "MDA")
	{
		if (iter < solution_iterations) 
			evaluate_update = true;					
	}

		
	if (evaluate_update)
	{
		message(0, "running inflation ensembles");
		if (par_dyn_state_names.size() > 0)
		{
			// before forward runing remove the dynamic states temporarly 
			// and replace them with prior dynamic states
			posterior_dyn_states = temp_remove_dyn_state(pe_lams);
		}
		oe_lams = run_lambda_ensembles(pe_lams, lam_vals, scale_vals, icycle);
	}
	/*
	if (dyn_states_names.size() > 0)
	{
		// after forward running put the dynamic states back into pe
		if (da_type == " VANILLA") {
			return_post_dyn_state(pe_lams, posterior_dyn_states);
		}
	}
	*/
		
	
	
	// ============= end of forward runs =================
	message(0, "evaluting inflation factor ensembles");
	message(1, "last mean: ", last_best_mean);
	message(1, "last stdev: ", last_best_std);

	ObservationEnsemble oe_lam_best;
	bool echo = false;
	if (verbose_level > 1)
		echo = true;

	//*******************************************************************************************
	// ===== loop over all subsets to find best inflation factor and best scaling factor ======
	if (evaluate_update)
	{
		for (int i = 0; i < pe_lams.size(); i++)
		{
			if (oe_lams[i].shape().first == 0)
				continue;
			vector<double> vals({ lam_vals[i],scale_vals[i] });
			if (pest_scenario.get_pestpp_options().get_ies_save_lambda_en())

			{
				ss.str("");
				ss << file_manager.get_base_filename() << "." << iter << "." << lam_vals[i] << ".lambda." << scale_vals[i] << ".scale.obs";

				if (pest_scenario.get_pestpp_options().get_save_binary())
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
			message(0, "phi summary for inflation, scale fac:", vals, echo);
			
			ph.report(true);


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
	}
	// ======== end of looping over subsets ==============
	if (da_type == "VANILLA")
	{
		pe = pe_lams[0];
		pe_post = posterior_dyn_states[0];
		if (evaluate_update)
		{
			oe = oe_lams[0];
		}
		return true;

	}
	
	if (da_type == "MDA") 
	{			
		if (iter == infl_facs.size()) // I thinks is also not needed
		{
			if (evaluate_update)
			if (par_dyn_state_names.size() > 0)
			{
				message(1, "Last MDA cycle; update dynamic forecast ");
				return_post_dyn_state(pe_lams, posterior_dyn_states);
			}			
		}
		pe = pe_lams[0];
		//save best posterior ...
		if (posterior_dyn_states.size() > 0)
			pe_post = posterior_dyn_states[0];
		if (evaluate_update)
			oe = oe_lams[0];
		return true;
	}
	


	if (best_idx == -1) // here break to do another iteration becuase prior phi is better than posterior phi
	{
		message(0, "WARNING:  unsuccessful lambda testing, resetting lambda to 10000.0");
		last_best_lam = 10000.0;
		return false;
	}

	//save best posterior ...
	pe_post = posterior_dyn_states[best_idx];
	posterior_dyn_states.clear();

	//use_ies = true;
	if (false)
	{
		pe = pe_lams[best_idx];
		oe = oe_lam_best;
		return true;
	}
	// ============= if you are here, then update is not good enough =======
	double acc_fac;
	double lam_inc;
	double lam_dec;
	acc_fac = da_ctl_params.get_dvalue("DA_ACCEPT_PHI_FAC");
	lam_inc = da_ctl_params.get_dvalue("DA_LAMBDA_INC_FAC");
	lam_dec = da_ctl_params.get_dvalue("DA_LAMBDA_DEC_FAC");

	// from this point till the end, it is only for the iterative mode
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

			vector<string> org_pe_idxs, org_oe_idxs;
			set<string> ssub;
			for (auto& i : subset_idxs)
				ssub.emplace(pe_names[i]);
			for (int i = 0; i < pe_names.size(); i++)
				if (ssub.find(pe_names[i]) == ssub.end())
				{
					pe_keep_names.push_back(pe_names[i]);
					//oe_keep_names.push_back(oe_names[i]);
				}
			ssub.clear();
			for (auto& i : subset_idxs)
				ssub.emplace(oe_names[i]);
			for (int i = 0; i < oe_names.size(); i++)
				if (ssub.find(oe_names[i]) == ssub.end())
				{
					oe_keep_names.push_back(oe_names[i]);
				}
			message(0, "phi summary for best lambda, scale fac: ", vector<double>({ lam_vals[best_idx],scale_vals[best_idx] }));
			ph.update(oe_lams[best_idx], pe_lams[best_idx]);
			ph.report(true);
			message(0, "running remaining realizations for best inflation factor, scale:", vector<double>({ lam_vals[best_idx],scale_vals[best_idx] }));

			//pe_keep_names and oe_keep_names are names of the remaining reals to eval
			performance_log->log_event("dropping subset idxs from remaining_pe_lam");
			remaining_pe_lam.keep_rows(pe_keep_names);
			performance_log->log_event("dropping subset idxs from remaining_oe_lam");
			remaining_oe_lam.keep_rows(oe_keep_names);
			//save these names for later
			org_pe_idxs = remaining_pe_lam.get_real_names();
			org_oe_idxs = remaining_oe_lam.get_real_names();
			///run
			//vector<int> fails = run_ensemble(remaining_pe_lam, remaining_oe_lam);

			/// ***************** Run remaining ****************************
			vector<ParameterEnsemble> pe_lams_remain;
			vector<ParameterEnsemble> posterior_dyn_states; //todo: remove, not used!
			pe_lams_remain.clear();
			pe_lams_remain.push_back(remaining_pe_lam);
			// TODO: we do not need this here becuase we have already posteriors from all the ensembles 
			if (true)
				if (par_dyn_state_names.size() > 0)
				{
					// before forward runing remove the dynamic states temporarly 
					// and replace them with prior dynamic states
					posterior_dyn_states = temp_remove_dyn_state(pe_lams_remain);
				}

			//oe_lams = run_lambda_ensembles(pe_lams, lam_vals, scale_vals);
			vector<int> fails = run_ensemble(pe_lams_remain[0], remaining_oe_lam, vector<int>(), icycle);
			if (false) // only no need for this here
				if (par_dyn_state_names.size() > 0)
				{
					// after forward running put the dynamic states back into pe
					return_post_dyn_state(pe_lams_remain, posterior_dyn_states);
				}
			remaining_pe_lam = pe_lams_remain[0];
			//clear posterior matrix
			posterior_dyn_states.clear();
			pe_lams_remain.clear();
			
			///*************************************************************

			//for testing
			if (pest_scenario.get_pestpp_options().get_ies_debug_fail_remainder())
				fails.push_back(0);

			//if any of the remaining runs failed
			if (fails.size() == org_pe_idxs.size())
				throw_em_error(string("all remaining realizations failed...something is prob wrong"));
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
				throw_em_error(string("all realization dropped after finishing subset runs...something might be wrong..."));
			}
			performance_log->log_event("updating phi");
			ph.update(oe_lam_best, pe_lams[best_idx]);
			best_mean = ph.get_mean(L2PhiHandler::phiType::COMPOSITE);
			best_std = ph.get_std(L2PhiHandler::phiType::COMPOSITE);
			message(1, "phi summary for entire ensemble using lambda,scale_fac ", vector<double>({ lam_vals[best_idx],scale_vals[best_idx] }));
			ph.report(true);
		}
	else //no subsets
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
	{   // new phi is decreasing or slightly above best!! 
		message(0, "updating parameter ensemble");
		performance_log->log_event("updating parameter ensemble");
		last_best_mean = best_mean;

		pe = pe_lams[best_idx];
		oe = oe_lam_best;
		if (best_std < last_best_std * acc_fac)
		{
			double new_lam = lam_vals[best_idx] * lam_dec;
			new_lam = (new_lam < lambda_min) ? lambda_min : new_lam;
			message(0, "updating inflation factor to ", new_lam);
			last_best_lam = new_lam;
		}
		else
		{
			message(0, "not updating inflation factor");
		}
		last_best_std = best_std;
	}

	else
	{ //  ************ phi is not reduced!!!
		//message(0, "not updating parameter ensemble");
		message(0, "only updating realizations with reduced phi");
		update_reals_by_phi(pe_lams[best_idx], oe_lam_best);
		ph.update(oe, pe);
		double new_lam = last_best_lam * lam_inc;
		new_lam = (new_lam > lambda_max) ? lambda_max : new_lam;
		message(0, "incresing inflation factor to: ", new_lam);
		last_best_lam = new_lam;
	}
	//report_and_save();
	return true;
}
void DataAssimilator::update_starting_state()
{
	
	vector<string> real_names, var_names;	
	real_names = pe.get_real_names();
	var_names = pe_post.get_var_names();
	pe_post.keep_rows(real_names);	
	pe.replace_col_vals(var_names, pe_post.get_eigen());
	
	
}
void DataAssimilator::return_post_dyn_state(vector<ParameterEnsemble>& pe_lams, vector<ParameterEnsemble> posterior_dyn_states)
{
	vector<string> real_names_; // it is not used!

	int ireal = 0;
	//real_names_ = pe_base.get_real_names();	

	for (auto _pe : pe_lams)

	{
		vector<string> real_names, var_names;
		ParameterEnsemble Ens = posterior_dyn_states[ireal];
		real_names = pe_lams[ireal].get_real_names();
		var_names = Ens.get_var_names();
		Ens.keep_rows(real_names);

		//real_names_ = _pe.get_real_names();		
		pe_lams[ireal].replace_col_vals(var_names, Ens.get_eigen());
		ireal = ireal + 1;
	}

}

vector<ParameterEnsemble> DataAssimilator::temp_remove_dyn_state(vector<ParameterEnsemble>& pe_lams)
{
		vector<string> real_names_;
		vector<ParameterEnsemble> poterior_dyn_states;
		int ireal = 0;
		real_names_ = pe_base.get_real_names();
		Eigen::MatrixXd matPrior;
		vector<string> static_parms_names;
		ParameterEnsemble local_base_pe;
		ParameterEnsemble current_pe;

		static_parms_names = pe_base.get_var_names();
		set<string> dyn_nam(par_dyn_state_names.begin(),par_dyn_state_names.end());		

		// get a vector of non-dynamic parameters
		auto pred = [&dyn_nam](const std::string& key) ->bool
		{
			return dyn_nam.find(key) != dyn_nam.end();
		};
		static_parms_names.erase(remove_if(static_parms_names.begin(), static_parms_names.end(), pred), static_parms_names.end());
		
		for (auto _pe : pe_lams)

		{
			real_names_ = _pe.get_real_names();			
			current_pe = _pe;
			current_pe.drop_cols(static_parms_names);			
			poterior_dyn_states.push_back(current_pe);			
			local_base_pe = pe_base;
			local_base_pe.keep_rows(real_names_);
			matPrior = local_base_pe.get_eigen(real_names_, par_dyn_state_names);
			pe_lams[ireal].replace_col_vals(par_dyn_state_names, matPrior);
			ireal = ireal + 1;
		}		
		return poterior_dyn_states;
	
}


void DataAssimilator::eig2csv(string name, Eigen::MatrixXd matrix)
{
	ofstream file(name.c_str());

	for (int i = 0; i < matrix.rows(); i++) {
		for (int j = 0; j < matrix.cols(); j++) {
			string str = std::to_string(matrix(i, j));
			if (j + 1 == matrix.cols()) {
				file << str;
			}
			else {
				file << str << ',';
			}
		}
		file << '\n';
	}
}
ParameterEnsemble DataAssimilator::kf_work(PerformanceLog* performance_log, unordered_map<string, Eigen::VectorXd>& par_resid_map,
	unordered_map<string, Eigen::VectorXd>& par_diff_map, 	unordered_map<string, Eigen::VectorXd>& obs_resid_map,
	unordered_map<string, Eigen::VectorXd>& obs_diff_map, unordered_map<string, Eigen::VectorXd>& obs_err_map, Localizer& localizer,
	unordered_map<string, double>& parcov_inv_map, 	unordered_map<string, double>& weight_map,
	ParameterEnsemble& pe_upgrade, double cur_lam,	unordered_map<string, pair<vector<string>, vector<string>>>& loc_map,
	unordered_map<string, Eigen::VectorXd>& Am_map, Localizer::How& how)
	
{

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


	//unique_lock<mutex> ctrl_guard(ctrl_lock, defer_lock);
	int maxsing, num_reals, verbose_level, pcount = 0;
	double eigthresh;
	bool use_approx;
	bool use_prior_scaling;
	bool use_localizer = false;
	bool loc_by_obs = true;

	maxsing = pe_upgrade.get_pest_scenario_ptr()->get_svd_info().maxsing;
	eigthresh = pe_upgrade.get_pest_scenario_ptr()->get_svd_info().eigthresh;
	use_approx = pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_ies_use_approx();
	use_prior_scaling = pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_ies_use_prior_scaling();
	num_reals = pe_upgrade.shape().first;
	verbose_level = pe_upgrade.get_pest_scenario_ptr()->get_pestpp_options().get_ies_verbose_level();
	
	if (how == Localizer::How::PARAMETERS)
		loc_by_obs = false;
	//break;
	
	ofstream f_thread;
	
	Eigen::MatrixXd par_resid, par_diff, Am;
	Eigen::MatrixXd obs_resid, obs_diff, obs_err, loc;
	Eigen::DiagonalMatrix<double, Eigen::Dynamic> weights, parcov_inv;
	vector<string> par_names, obs_names;
	while (true)
	{
		//unique_lock<mutex> next_guard(next_lock, defer_lock);
		par_names.clear();
		obs_names.clear();
		use_localizer = false;
		//the end condition
		//the end condition

		if (verbose_level > 2)
		{
			f_thread << iter;
			for (auto name : par_names)
				f_thread << "," << name;
			for (auto name : obs_names)
				f_thread << "," << name;
			f_thread << endl;
		}
		obs_names = loc_map.at("all").first;
		par_names = loc_map.at("all").second;
		par_resid.resize(0, 0);
		par_diff.resize(0, 0);
		obs_resid.resize(0, 0);
		obs_diff.resize(0, 0);
		loc.resize(0, 0);
		Am.resize(0, 0);
		weights.resize(0);
		parcov_inv.resize(0);
		Am.resize(0, 0);

		obs_diff = local_utils::get_matrix_from_map(num_reals, obs_names, obs_diff_map);
		par_diff = local_utils::get_matrix_from_map(num_reals, par_names, par_diff_map);
		par_resid = local_utils::get_matrix_from_map(num_reals, par_names, par_resid_map);
		weights = local_utils::get_matrix_from_map(obs_names, weight_map);
		obs_resid = local_utils::get_matrix_from_map(num_reals, obs_names, obs_resid_map);
		obs_err = local_utils::get_matrix_from_map(num_reals, obs_names, obs_err_map);
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


		// *******************************************************
		par_diff.transposeInPlace();
		obs_diff.transposeInPlace();
		obs_resid.transposeInPlace();
		par_resid.transposeInPlace();
		obs_err.transposeInPlace();


		//local_utils::save_mat(verbose_level, thread_id, iter, t_count, "obs_resid", obs_resid);
		Eigen::MatrixXd scaled_residual = weights * obs_resid;

		//local_utils::save_mat(verbose_level, thread_id, iter, t_count, "par_resid", par_resid);
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

		if (use_localizer) // todo: check localizer for kf
		{
			if (loc_by_obs)
				par_diff = par_diff.cwiseProduct(loc);
			else
				obs_diff = obs_diff.cwiseProduct(loc);

		}

		obs_diff = scale * obs_diff;// (H-Hm)/sqrt(N-1)		
		par_diff = scale * par_diff;// (K-Km)/sqrt(N-1)		
		obs_err = scale * obs_err; //  (E-Em)/sqrt(N-1)	
		obs_resid = (-1.0) * obs_resid; // D-H  , the negative one is becuase obs_resid is calculated as H-D ???
		if (true) // just debuging
		{
			eig2csv("kdash.csv", par_diff);
			eig2csv("hdash.csv", obs_diff);
			eig2csv("H_D.csv", obs_resid);
			eig2csv("noise.csv", obs_err);
		}


		Eigen::MatrixXd s2_, upgrade_1, s, V, U, cum_sum;


		SVD_REDSVD rsvd;
		Eigen::MatrixXd C;
		C = obs_diff + (cur_lam*obs_err); // curr_lam is the inflation factor 
		eig2csv("C.csv", C);
		Eigen::VectorXd s2;
		//eigthresh = 0;
		if (false)
		{
			rsvd.solve_ip(C, s, U, V, eigthresh, maxsing);
		}
		else
		{
		
			//Eigen::BDCSVD<Eigen::MatrixXd> svd(C, Eigen::ComputeThinU | Eigen::ComputeThinV);
			Eigen::JacobiSVD<Eigen::MatrixXd> svd(C, Eigen::ComputeThinU | Eigen::ComputeThinV);
			U = svd.matrixU();
			s = svd.singularValues();
			
			int ssize = s.size();			
			double threshold_frac;
			s2 = s.cwiseProduct(s);
			double s_sum = s.sum();
			cum_sum = s;
			for (int i=0; i<ssize; i++)
			{			
				
				threshold_frac = s(i) / s_sum;
				if (threshold_frac < eigthresh)
					s2(i) = 0;				
			}		

	    }
		eig2csv("s.csv", s);
		
		
		
		//obs_diff.resize(0, 0);
		//local_utils::save_mat(verbose_level, thread_id, iter, t_count, "Ut", Ut);
		//local_utils::save_mat(verbose_level, thread_id, iter, t_count, "s", s);
		//local_utils::save_mat(verbose_level, thread_id, iter, t_count, "V", V);

		

		// old // ivec = ((Eigen::VectorXd::Ones(s2.size()) * (cur_lam + 1.0)) + s2).asDiagonal().inverse();

		s2_ = s2.asDiagonal().inverse();
		for (int i=0; i < s2.size(); i++)
		{
			if (s2(i)<1e-50)
			{
				s2_(i, i) = 0;
			}
		}
		eig2csv("s2_.csv", s2_);
		//local_utils::save_mat(verbose_level, thread_id, iter, t_count, "ivec", ivec);
		//old// Eigen::MatrixXd X1 = Ut * scaled_residual;
		Eigen::MatrixXd X1 = s2_ * U.transpose();
		eig2csv("X1.csv", X1);
		X1 = X1 * obs_resid; 
		eig2csv("X2.csv", X1);
		X1 = U * X1;
		eig2csv("X3.csv", X1);
		X1 = obs_diff.transpose() * X1;
		eig2csv("X4.csv", X1);
		upgrade_1 = par_diff * X1;
		eig2csv("upgrade.csv", upgrade_1);
		pe_upgrade.add_2_cols_ip(par_names, upgrade_1.transpose());
		return pe_upgrade;
		break;
		
		//X1 = 1;
		//local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X1", X1);
		//Eigen::MatrixXd X2 = ivec * X1;
		//X1.resize(0, 0);

		//local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X2", X2);
		//Eigen::MatrixXd X3 = V * s.asDiagonal() * X2;
		//X2.resize(0, 0);

		//local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X3", X3);
		/*if (use_prior_scaling)
		{
			upgrade_1 = -1.0 * parcov_inv * par_diff * X3;
		}
		else
		{
			upgrade_1 = -1.0 * par_diff * X3;
		}*/
		//upgrade_1.transposeInPlace();
		//local_utils::save_mat(verbose_level, thread_id, iter, t_count, "upgrade_1", upgrade_1);
		//X3.resize(0, 0);


		/*Eigen::MatrixXd upgrade_2;
		if ((!use_approx) && (iter > 1))
		{
			//local_utils::save_mat(verbose_level, thread_id, iter, t_count, "Am", Am);
			Eigen::MatrixXd x4 = Am.transpose() * scaled_par_resid;
			//local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X4", x4);

			par_resid.resize(0, 0);

			Eigen::MatrixXd x5 = Am * x4;
			x4.resize(0, 0);
			Am.resize(0, 0);

			//local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X5", x5);
			Eigen::MatrixXd x6 = par_diff.transpose() * x5;
			x5.resize(0, 0);

			//local_utils::save_mat(verbose_level, thread_id, iter, t_count, "X6", x6);
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
			//local_utils::save_mat(verbose_level, thread_id, iter, t_count, "upgrade_2", upgrade_2);
			upgrade_2.resize(0, 0);

		} */

		//unique_lock<mutex> put_guard(put_lock, defer_lock);
		/*
		while (true)
		{
			if (put_guard.try_lock())
			{
				pe_upgrade.add_2_cols_ip(par_names, upgrade_1);
				put_guard.unlock();
				break;
			}
		}
		*/
	}

}



//void DataAssimilator::report_and_save()
//{
//	ofstream& frec = file_manager.rec_ofstream();
//	frec << endl << "  ---  DataAssimilator iteration " << iter << " report  ---  " << endl;
//	frec << "   number of active realizations:  " << pe.shape().first << endl;
//	frec << "   number of model runs:           " << run_mgr_ptr->get_total_runs() << endl;
//
//	cout << endl << "  ---  DataAssimilator iteration " << iter << " report  ---  " << endl;
//	cout << "   number of active realizations:   " << pe.shape().first << endl;
//	cout << "   number of model runs:            " << run_mgr_ptr->get_total_runs() << endl;
//
//	stringstream ss;
//	if (pest_scenario.get_pestpp_options().get_save_binary())
//	{
//		ss << file_manager.get_base_filename() << "." << iter << ".obs.jcb";
//		oe.to_binary(ss.str());
//	}
//	else
//	{
//		ss << file_manager.get_base_filename() << "." << iter << ".obs.csv";
//		oe.to_csv(ss.str());
//	}
//	frec << "      current obs ensemble saved to " << ss.str() << endl;
//	cout << "      current obs ensemble saved to " << ss.str() << endl;
//	ss.str("");
//	if (pest_scenario.get_pestpp_options().get_save_binary())
//	{
//		ss << file_manager.get_base_filename() << "." << iter << ".par.jcb";
//		pe.to_binary(ss.str());
//	}
//	else
//	{
//		ss << file_manager.get_base_filename() << "." << iter << ".par.csv";
//		pe.to_csv(ss.str());
//	}
//	//ss << file_manager.get_base_filename() << "." << iter << ".par.csv";
//	//pe.to_csv(ss.str());
//	frec << "      current par ensemble saved to " << ss.str() << endl;
//	cout << "      current par ensemble saved to " << ss.str() << endl;
//
//
//
//}
//

void DataAssimilator::set_subset_idx(int size)
{
	//map<int,int> subset_idx_map;
	subset_idxs.clear();
	int nreal_subset;
	
	//CtlPar_container da_ctl_params = pest_scenario.get_pestpp_options().da_ctl_params;
	if ((da_type == "MDA")|| (da_type == "VANILLA"))
	{
		int nnreal = pe.get_real_names().size();
		nreal_subset = nnreal;
		use_subset = true;
	}
	else if (da_type == "ITERATIVE")
	{
		nreal_subset = da_ctl_params.get_ivalue("DA_SUBSET_SIZE");
		if (nreal_subset <= 0)
		{
			int nnreal = pe.get_real_names().size();
			nreal_subset = nnreal;
		}
		use_subset = true;
	}
	
	if ((!use_subset) || (nreal_subset >= size))
	{
		for (int i = 0; i < size; i++)
			subset_idxs.push_back(i);
		return;
	}
	vector<string> pe_names = pe.get_real_names();

	vector<string>::iterator bidx = find(pe_names.begin(), pe_names.end(), BASE_REAL_NAME);
	if (bidx != pe_names.end())
	{

		subset_idxs.push_back(bidx - pe_names.begin());
	}
	//int size = pe.shape().first;
	string how = pest_scenario.get_pestpp_options().get_ies_subset_how();
	how = da_ctl_params.get_svalue("DA_SUBSET_HOW");

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
		std::uniform_int_distribution<int> uni(0, size - 1);
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
		for (pi; pi != end; ++pi)
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
		throw_em_error("unknown 'subset_how'");
	}
	stringstream ss;
	for (auto i : subset_idxs)
		ss << i << ":" << pe_names[i] << ", ";
	message(1, "subset idx:pe real name: ", ss.str());
	return;
	//return subset_idx_map;
}

//vector<ObservationEnsemble> DataAssimilator::run_lambda_ensembles(vector<ParameterEnsemble>& pe_lams, vector<double>& lam_vals, vector<double>& scale_vals)
//{
//	ofstream& frec = file_manager.rec_ofstream();
//	stringstream ss;
//	ss << "queuing " << pe_lams.size() << " ensembles";
//	performance_log->log_event(ss.str());
//	run_mgr_ptr->reinitialize();
//
//	set_subset_idx(pe_lams[0].shape().first); // get realization index to test tuning parameters
//	vector<map<int, int>> real_run_ids_vec;
//	//ParameterEnsemble pe_lam;
//	//for (int i=0;i<pe_lams.size();i++)
//	for (auto& pe_lam : pe_lams)
//	{
//		try
//		{
//			real_run_ids_vec.push_back(pe_lam.add_runs(run_mgr_ptr, subset_idxs,icycle));
//		}
//		catch (const exception & e)
//		{
//			stringstream ss;
//			ss << "run_ensemble() error queueing runs: " << e.what();
//			throw_em_error(ss.str());
//		}
//		catch (...)
//		{
//			throw_em_error(string("run_ensembles() error queueing runs"));
//		}
//	}
//	performance_log->log_event("making runs");
//	try
//	{
//
//		run_mgr_ptr->run();
//	}
//	catch (const exception & e)
//	{
//		stringstream ss;
//		ss << "error running ensembles: " << e.what();
//		throw_em_error(ss.str());
//	}
//	catch (...)
//	{
//		throw_em_error(string("error running ensembles"));
//	}
//
//	performance_log->log_event("processing runs");
//	vector<int> failed_real_indices;
//	vector<ObservationEnsemble> obs_lams;
//	//ObservationEnsemble _oe = oe;//copy
//	//if (subset_size < pe_lams[0].shape().first)
//	//	_oe.keep_rows(subset_idxs);
//	map<int, int> real_run_ids;
//	//for (auto &real_run_ids : real_run_ids_vec)
//	for (int i = 0; i < pe_lams.size(); i++)
//	{
//		ObservationEnsemble _oe = oe;//copy
//		vector<double> rep_vals{ lam_vals[i],scale_vals[i] };
//		real_run_ids = real_run_ids_vec[i];
//		//if using subset, reset the real_idx in real_run_ids to be just simple counter
//		if (subset_size < pe_lams[0].shape().first)
//		{
//			_oe.keep_rows(subset_idxs);
//			int ireal = 0;
//			map<int, int> temp;
//			for (auto& rri : real_run_ids)
//			{
//				temp[ireal] = rri.second;
//				ireal++;
//			}
//
//			real_run_ids = temp;
//		}
//
//		try
//		{
//			failed_real_indices = _oe.update_from_runs(real_run_ids, run_mgr_ptr);
//		}
//		catch (const exception & e)
//		{
//			stringstream ss;
//			ss << "error processing runs for lambda,scale: " << lam_vals[i] << ',' << scale_vals[i] << ':' << e.what();
//			throw_em_error(ss.str());
//		}
//		catch (...)
//		{
//			stringstream ss;
//			ss << "error processing runs for lambda,scale: " << lam_vals[i] << ',' << scale_vals[i];
//			throw_em_error(ss.str());
//		}
//
//		if (pest_scenario.get_pestpp_options().get_ies_debug_fail_subset())
//			failed_real_indices.push_back(real_run_ids.size() - 1);
//
//		if (failed_real_indices.size() > 0)
//		{
//			stringstream ss;
//			vector<string> par_real_names = pe.get_real_names();
//			vector<string> obs_real_names = oe.get_real_names();
//			vector<string> failed_par_names, failed_obs_names;
//			string oname, pname;
//			ss << "the following par:obs realization runs failed for lambda,scale " << lam_vals[i] << ',' << scale_vals[i] << "-->";
//			for (auto& i : failed_real_indices)
//			{
//				pname = par_real_names[subset_idxs[i]];
//				oname = obs_real_names[subset_idxs[i]];
//				failed_par_names.push_back(pname);
//				failed_obs_names.push_back(oname);
//				ss << pname << ":" << oname << ',';
//			}
//			string s = ss.str();
//			message(1, s);
//			if (failed_real_indices.size() == _oe.shape().first)
//			{
//				message(0, "WARNING: all realizations failed for lambda, scale :", rep_vals);
//				_oe = ObservationEnsemble();
//
//			}
//			else
//			{
//				performance_log->log_event("dropping failed realizations");
//				//_oe.drop_rows(failed_real_indices);
//				//pe_lams[i].drop_rows(failed_real_indices);
//				_oe.drop_rows(failed_obs_names);
//				pe_lams[i].drop_rows(failed_par_names);
//			}
//
//		}
//		obs_lams.push_back(_oe);
//	}
//	return obs_lams;
//}
//
//
//vector<int> DataAssimilator::run_ensemble(ParameterEnsemble& _pe, ObservationEnsemble& _oe, const vector<int>& real_idxs)
//{
//	stringstream ss;
//	ss << "queuing " << _pe.shape().first << " runs";
//	performance_log->log_event(ss.str());
//	run_mgr_ptr->reinitialize();
//	map<int, int> real_run_ids;
//	try
//	{
//		real_run_ids = _pe.add_runs(run_mgr_ptr, real_idxs,icycle);
//	}
//	catch (const exception & e)
//	{
//		stringstream ss;
//		ss << "run_ensemble() error queueing runs: " << e.what();
//		throw_em_error(ss.str());
//	}
//	catch (...)
//	{
//		throw_em_error(string("run_ensemble() error queueing runs"));
//	}
//	performance_log->log_event("making runs");
//	try
//	{
//		run_mgr_ptr->run();
//	}
//	catch (const exception & e)
//	{
//		stringstream ss;
//		ss << "error running ensemble: " << e.what();
//		throw_em_error(ss.str());
//	}
//	catch (...)
//	{
//		throw_em_error(string("error running ensemble"));
//	}
//
//	performance_log->log_event("processing runs");
//	if (real_idxs.size() > 0)
//	{
//		_oe.keep_rows(real_idxs);
//	}
//	vector<int> failed_real_indices;
//	try
//	{
//		failed_real_indices = _oe.update_from_runs(real_run_ids, run_mgr_ptr);
//	}
//	catch (const exception & e)
//	{
//		stringstream ss;
//		ss << "error processing runs: " << e.what();
//		throw_em_error(ss.str());
//	}
//	catch (...)
//	{
//		throw_em_error(string("error processing runs"));
//	}
//	//for testing
//	//failed_real_indices.push_back(0);
//
//	if (failed_real_indices.size() > 0)
//	{
//		stringstream ss;
//		vector<string> par_real_names = _pe.get_real_names();
//		vector<string> obs_real_names = _oe.get_real_names();
//		ss << "the following par:obs realization runs failed: ";
//		for (auto& i : failed_real_indices)
//		{
//			ss << par_real_names[i] << ":" << obs_real_names[i] << ',';
//		}
//		performance_log->log_event(ss.str());
//		message(1, "failed realizations: ", failed_real_indices.size());
//		string s = ss.str();
//		message(1, s);
//		performance_log->log_event("dropping failed realizations");
//		_pe.drop_rows(failed_real_indices);
//		_oe.drop_rows(failed_real_indices);
//	}
//	return failed_real_indices;
//}
//

void DataAssimilator::finalize()
{

}

map<int, map<string, double>> process_da_obs_cycle_table(Pest& pest_scenario, vector <int>& ncycles_in_tables, ofstream& fout_rec, set<string>& obs_in_tbl)
{
	//process da par cycle table
	string filename = pest_scenario.get_pestpp_options().get_da_obs_cycle_table();
	map<int, map<string, double>> obs_cycle_info;
	if (filename.size() > 0)
	{
		fout_rec << "processing 'DA_OBSEERVATION_CYCLE_TABLE' file " << filename;
		pest_utils::ExternalCtlFile cycle_table(filename);
		cycle_table.read_file(fout_rec);
		vector<string> col_names = cycle_table.get_col_names();
		fout_rec << "...using the first column ('" << col_names[0] << "') as observation names" << endl;
		vector<string> onames = pest_scenario.get_ctl_ordered_obs_names();
		set<string> obs_names(onames.begin(), onames.end());
		onames = cycle_table.get_col_string_vector(col_names[0]);
		vector<string> missing, zeroweight, notneg, tbl_obs_names;
		ObservationInfo* oi = pest_scenario.get_observation_info_ptr();
		for (auto oname : onames)
		{
			pest_utils::upper_ip(oname);
			if (obs_names.find(oname) == obs_names.end())
			{
				missing.push_back(oname);
				continue;
			}
			else if (oi->get_weight(oname) == 0.0)
			{
				zeroweight.push_back(oname);
				continue;
			}
			else if (oi->get_observation_rec_ptr(oname)->cycle != -1)
			{
				notneg.push_back(oname);
				continue;
			}
			else
			{
				tbl_obs_names.push_back(oname);
			}
		}
		if (missing.size() > 0)
		{
			fout_rec << "ERROR: The following observations in DA_OBSERVATION_CYCLE_TABLE are not the control file:" << endl;
			for (auto m : missing)
			{
				fout_rec << m << endl;
			}
			throw runtime_error("DA_OBSERVATION_CYCLE_TABLE contains missing observations, see rec file for listing");
		}
		if (zeroweight.size() > 0)
		{
			fout_rec << "ERROR: The following observations in DA_OBSERVATION_CYCLE_TABLE have zero weight:" << endl;
			for (auto p : zeroweight)
			{
				fout_rec << p << endl;
			}
			throw runtime_error("DA_OBSERVATION_CYCLE_TABLE contains zero-weighted observations, see rec file for listing");
		}
		if (notneg.size() > 0)
		{
			fout_rec << "ERROR: The following OBSERVATIONS in DA_OBSERVATION_CYCLE_TABLE do not have cycle=-1:" << endl;
			for (auto p : notneg)
			{
				fout_rec << p << endl;
			}
			throw runtime_error("DA_OBSERVATION_CYCLE_TABLE contains observations with cycle!=-1, see rec file for listing");
		}
		//process the remaining columns - these should be cycle numbers
		obs_in_tbl.insert(tbl_obs_names.begin(), tbl_obs_names.end());
		string col_name;
		vector<string> cycle_vals;
		int cycle;
		double val;
		bool parse_fail = false;
		for (int i = 1; i < col_names.size(); i++)
		{
			col_name = col_names[i];

			try
			{
				cycle = stoi(col_name);
			}
			catch (...)
			{
				fout_rec << "ERROR: could not parse DA_OBSERVATION_CYCLE_TABLE column '" << col_name << "' to integer" << endl;
				throw runtime_error("ERROR parsing DA_OBSERVATIONS_CYCLE_TABLE column " + col_name + " to integer");
			}
			if (obs_cycle_info.find(cycle) != obs_cycle_info.end())
			{
				throw runtime_error("ERROR: DA_OBSERVATION_CYCLE_TABLE cycle column '" + col_name + "' listed more than once");
			}
			cycle_vals = cycle_table.get_col_string_vector(col_name);
			map<string, double> cycle_map;
			for (int i = 0; i < cycle_vals.size(); i++)
			{
				try
				{
					val = stod(cycle_vals[i]);
				}
				catch (...)
				{
					fout_rec << "WARNING: error parsing '" << cycle_vals[i] << "' for observation " << tbl_obs_names[i] << " in cycle " << cycle << ", continuing..." << endl;
					parse_fail = true;
					continue;
				}
				cycle_map[tbl_obs_names[i]] = val;
				ncycles_in_tables.push_back(cycle);
				ncycles_in_tables.push_back(cycle);

				// make sure cycle numbers are sorted and unique
				sort(ncycles_in_tables.begin(), ncycles_in_tables.end());
				ncycles_in_tables.erase(unique(ncycles_in_tables.begin(),
					ncycles_in_tables.end()),
					ncycles_in_tables.end());
			}
			obs_cycle_info[cycle] = cycle_map;
		}
		if (parse_fail)
		{
			cout << "WARNING: error parsing at least one cycle-based observation value" << endl;
			cout << "         from DA_OBSERVATION_CYCLE_TABLE, see rec file for listing." << endl;
		}


	}
	return obs_cycle_info;


}

map<int, map<string, double>> process_da_weight_cycle_table(Pest& pest_scenario, vector <int>& ncycles_in_tables, ofstream& fout_rec, set<string>& obs_in_tbl)
{
	//process da par cycle table
	string filename = pest_scenario.get_pestpp_options().get_da_weight_cycle_table();
	map<int, map<string, double>> weight_cycle_info;
	if (filename.size() > 0)
	{
		fout_rec << "processing 'DA_WEIGHT_CYCLE_TABLE' file " << filename;
		pest_utils::ExternalCtlFile cycle_table(filename);
		cycle_table.read_file(fout_rec);
		vector<string> col_names = cycle_table.get_col_names();
		fout_rec << "...using the first column ('" << col_names[0] << "') as observation names" << endl;
		vector<string> onames = pest_scenario.get_ctl_ordered_obs_names();
		set<string> obs_names(onames.begin(), onames.end());
		onames = cycle_table.get_col_string_vector(col_names[0]);
		vector<string> missing, zeroweight, notneg, tbl_obs_names;
		ObservationInfo* oi = pest_scenario.get_observation_info_ptr();
		for (auto oname : onames)
		{
			pest_utils::upper_ip(oname);
			if (obs_names.find(oname) == obs_names.end())
			{
				missing.push_back(oname);
				continue;
			}
			else if (oi->get_observation_rec_ptr(oname)->cycle != -1)
			{
				notneg.push_back(oname);
				continue;
			}
			else
			{
				tbl_obs_names.push_back(oname);
			}
		}
		if (missing.size() > 0)
		{
			fout_rec << "ERROR: The following observations in DA_WEIGHT_CYCLE_TABLE are not the control file:" << endl;
			for (auto m : missing)
			{
				fout_rec << m << endl;
			}
			throw runtime_error("DA_WEIGHT_CYCLE_TABLE contains missing observations, see rec file for listing");
		}
		
		if (notneg.size() > 0)
		{
			fout_rec << "ERROR: The following OBSERVATIONS in DA_WEIGHT_CYCLE_TABLE do not have cycle=-1:" << endl;
			for (auto p : notneg)
			{
				fout_rec << p << endl;
			}
			throw runtime_error("DA_WEIGHT_CYCLE_TABLE contains observations with cycle!=-1, see rec file for listing");
		}
		//process the remaining columns - these should be cycle numbers
		obs_in_tbl.insert(tbl_obs_names.begin(), tbl_obs_names.end());
		string col_name;
		vector<string> cycle_vals;
		int cycle;
		double val;
		bool parse_fail = false;
		for (int i = 1; i < col_names.size(); i++)
		{
			col_name = col_names[i];

			try
			{
				cycle = stoi(col_name);
			}
			catch (...)
			{
				fout_rec << "ERROR: could not parse DA_WEIGHT_CYCLE_TABLE column '" << col_name << "' to integer" << endl;
				throw runtime_error("ERROR parsing DA_WEIGHT_CYCLE_TABLE column " + col_name + " to integer");
			}
			if (weight_cycle_info.find(cycle) != weight_cycle_info.end())
			{
				throw runtime_error("ERROR: DA_WEIGHT_CYCLE_TABLE cycle column '" + col_name + "' listed more than once");
			}
			cycle_vals = cycle_table.get_col_string_vector(col_name);
			map<string, double> cycle_map;
			for (int i = 0; i < cycle_vals.size(); i++)
			{
				try
				{
					val = stod(cycle_vals[i]);
				}
				catch (...)
				{
					fout_rec << "WARNING: error parsing weight '" << cycle_vals[i] << "' for observation " << tbl_obs_names[i] << " in cycle " << cycle << ", continuing..." << endl;
					parse_fail = true;
					continue;
				}
				cycle_map[tbl_obs_names[i]] = val;
			}
			weight_cycle_info[cycle] = cycle_map;
			ncycles_in_tables.push_back(cycle);

			// make sure cycle numbers are sorted and unique
			sort(ncycles_in_tables.begin(), ncycles_in_tables.end());
			ncycles_in_tables.erase(unique(ncycles_in_tables.begin(),
								ncycles_in_tables.end()),
								ncycles_in_tables.end());
		}
		if (parse_fail)
		{
			cout << "WARNING: error parsing at least one cycle-based observation value" << endl;
			cout << "         from DA_WEIGHT_CYCLE_TABLE, see rec file for listing." << endl;
		}


	}
	return weight_cycle_info;


}

map<int, map<string, double>> process_da_par_cycle_table(Pest& pest_scenario, vector <int>& ncycles_in_tables, ofstream& fout_rec)
{
	//process da par cycle table
	string filename = pest_scenario.get_pestpp_options().get_da_par_cycle_table();
	map<int, map<string, double>> par_cycle_info;
	if (filename.size() > 0)
	{
		fout_rec << "processing 'DA_PARAMETER_CYCLE_TABLE' file " << filename;
		pest_utils::ExternalCtlFile cycle_table(filename);
		cycle_table.read_file(fout_rec);
		vector<string> col_names = cycle_table.get_col_names();
		fout_rec << "...using the first column ('" << col_names[0] << "') as parameter names" << endl;
		vector<string> pnames = pest_scenario.get_ctl_ordered_par_names();
		set<string> par_names(pnames.begin(), pnames.end());
		pnames = cycle_table.get_col_string_vector(col_names[0]);
		ParameterInfo pi = pest_scenario.get_ctl_parameter_info();
		vector<string> missing, notfixed, notneg, tbl_par_names;

		for (auto pname : pnames)
		{
			pest_utils::upper_ip(pname);
			if (par_names.find(pname) == par_names.end())
			{
				missing.push_back(pname);
				continue;
			}
			else if (pi.get_parameter_rec_ptr(pname)->tranform_type != ParameterRec::TRAN_TYPE::FIXED)
			{
				notfixed.push_back(pname);
				continue;
			}
			else if (pi.get_parameter_rec_ptr(pname)->cycle != -1)
			{
				notneg.push_back(pname);
				continue;
			}
			else
			{
				tbl_par_names.push_back(pname);
			}
		}
		if (missing.size() > 0)
		{
			fout_rec << "ERROR: The following parameters in DA_PARAMETER_CYCLE_TABLE are not the control file:" << endl;
			for (auto p : missing)
			{
				fout_rec << p << endl;
			}
			throw runtime_error("DA_PARAMTER_CYCLE_TABLE contains missing parameters, see rec file for listing");
		}
		if (notfixed.size() > 0)
		{
			fout_rec << "ERROR: The following parameters in DA_PARAMETER_CYCLE_TABLE are not 'fixed':" << endl;
			for (auto p : notfixed)
			{
				fout_rec << p << endl;
			}
			throw runtime_error("DA_PARAMTER_CYCLE_TABLE contains non-fixed parameters, see rec file for listing");
		}
		if (notneg.size() > 0)
		{
			fout_rec << "ERROR: The following parameters in DA_PARAMETER_CYCLE_TABLE do not have cycle=-1:" << endl;
			for (auto p : notneg)
			{
				fout_rec << p << endl;
			}
			throw runtime_error("DA_PARAMTER_CYCLE_TABLE contains parameters with cycle!=-1, see rec file for listing");
		}
		//process the remaining columns - these should be cycle numbers
		string col_name;
		vector<string> cycle_vals;
		int cycle;
		double val;
		bool parse_fail = false;
		for (int i = 1; i < col_names.size(); i++)
		{
			col_name = col_names[i];

			try
			{
				cycle = stoi(col_name);
			}
			catch (...)
			{
				fout_rec << "ERROR: could not parse DA_PARAMETER_CYCLE_TABLE column '" << col_name << "' to integer" << endl;
				throw runtime_error("ERROR parsing DA_PARAMETER CYCLE TABLE column " + col_name + " to integer");
			}
			if (par_cycle_info.find(cycle) != par_cycle_info.end())
			{
				throw runtime_error("ERROR: DA_PARAMETER_CYCLE_TABLE cycle column '" + col_name + "' listed more than once");
			}
			cycle_vals = cycle_table.get_col_string_vector(col_name);
			map<string, double> cycle_map;
			for (int i = 0; i < cycle_vals.size(); i++)
			{
				try
				{
					val = stod(cycle_vals[i]);
				}
				catch (...)
				{
					fout_rec << "WARNING: error parsing '" << cycle_vals[i] << "' for parameter " << tbl_par_names[i] << " in cycle " << cycle << ", continuing..." << endl;
					parse_fail = true;
					continue;
				}
				cycle_map[tbl_par_names[i]] = val;
			}
			par_cycle_info[cycle] = cycle_map;
			ncycles_in_tables.push_back(cycle);

			// make sure cycle numbers are sorted and unique
			sort(ncycles_in_tables.begin(), ncycles_in_tables.end());
			ncycles_in_tables.erase(unique(ncycles_in_tables.begin(), 
											ncycles_in_tables.end()),
											ncycles_in_tables.end());


		}
		if (parse_fail)
		{
			cout << "WARNING: error parsing at least one cycle-based parameter value" << endl;
			cout << "         from DA_PARAMETER_CYCLE_TABLE, see rec file for listing." << endl;
		}


	}
	return par_cycle_info;
}


void write_global_phi_info(int icycle, ofstream& f_phi, DataAssimilator& da, vector<string>& init_real_names)
{
	int iter = da.get_iter();
	f_phi << icycle << "," << iter << "," << da.get_phi_handler().get_mean(L2PhiHandler::phiType::ACTUAL);
	f_phi << "," << da.get_phi_handler().get_std(L2PhiHandler::phiType::ACTUAL);
	f_phi << "," << da.get_phi_handler().get_min(L2PhiHandler::phiType::ACTUAL);
	f_phi << "," << da.get_phi_handler().get_max(L2PhiHandler::phiType::ACTUAL);
	map<string, double> final_actual = da.get_phi_handler().get_phi_map(L2PhiHandler::phiType::ACTUAL);
	for (auto rname : init_real_names)
	{
		if (final_actual.find(rname) != final_actual.end())
			f_phi << "," << final_actual.at(rname);
		else
			f_phi << ",";
	}
	f_phi << endl;
}

void generate_global_ensembles(DataAssimilator& da, ofstream& fout_rec, ParameterEnsemble& curr_pe, ObservationEnsemble& curr_oe)
{
	//generate a parent ensemble which includes all parameters across all cycles
	cout << "...preparing global parameter ensemble for all parameters across all cycles" << endl;
	fout_rec << "...preparing global parameter ensemble for all parameters across all cycles" << endl;
	Pest pest_scenario = da.get_pest_scenario();
	da.initialize_parcov();
	da.initialize_pe(*da.get_parcov_ptr());
	curr_pe = da.get_pe();
	Parameters pars = pest_scenario.get_ctl_parameters();
	ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();
	if (curr_pe.get_trans_status() != ParameterEnsemble::transStatus::NUM)
		throw runtime_error("parameter ensemble not in numeric transform status");
	pts.ctl2numeric_ip(pars);
	Eigen::VectorXd v = pars.get_data_eigen_vec(pest_scenario.get_ctl_ordered_adj_par_names());
	if (pest_scenario.get_control_info().noptmax == 0)
	{
		curr_pe.reserve(vector<string>{BASE_REAL_NAME}, pest_scenario.get_ctl_ordered_adj_par_names());
		curr_pe.update_real_ip(BASE_REAL_NAME, v);
	}
	else if (pest_scenario.get_pestpp_options().get_ies_include_base())
	{
		cout << "...replacing last realization with 'base' realization in global parameter ensemble" << endl;
		fout_rec << "...replacing last realization with 'base' realization in global parameter ensemble" << endl;
		curr_pe.replace(curr_pe.shape().first - 1, pars, BASE_REAL_NAME);
	}

	cout << "...preparing global observation ensemble for all observations across all cycles" << endl;
	fout_rec << "...preparing global observation ensemble for all observations across all cycles" << endl;
	mt19937 rand_gen = da.get_rand_gen();
	curr_oe.set_pest_scenario(&pest_scenario);
	curr_oe.set_rand_gen(&rand_gen);
	curr_oe.reserve(curr_pe.get_real_names(), pest_scenario.get_ctl_ordered_obs_names());

	try
	{
		curr_pe.check_for_dups();
	}
	catch (const exception& e)
	{
		string message = e.what();
		throw runtime_error("error in parameter ensemble: " + message);
	}
	//todo: add the base realization here so that it is present thru out the cycles below

	cout << "...global parameter ensemble has " << curr_pe.shape().first << " rows and " << curr_pe.shape().second << " columns" << endl;
	fout_rec << "...global parameter ensemble has " << curr_pe.shape().first << " rows and " << curr_pe.shape().second << " columns" << endl;
	curr_pe.to_csv(da.get_file_manager().get_base_filename() + ".global.prior.pe.csv");
	cout << "...global observation ensemble has " << curr_oe.shape().first << " rows and " << curr_oe.shape().second << " columns" << endl;
	fout_rec << "...global observation ensemble has " << curr_oe.shape().first << " rows and " << curr_oe.shape().second << " columns" << endl;

	return;
}

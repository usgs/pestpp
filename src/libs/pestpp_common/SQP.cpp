#include <random>
#include <map>
#include <iomanip>
#include <mutex>
#include <thread>
#include "Ensemble.h"
#include "RestartController.h"
#include "utilities.h"
#include "Ensemble.h"
#include "SQP.h"
#include "ObjectiveFunc.h"
#include "covariance.h"
#include "RedSVD-h.h"
#include "SVDPackage.h"
#include "eigen_tools.h"
#include "EnsembleMethodUtils.h"
#include "constraints.h"



SeqQuadProgram::SeqQuadProgram(Pest &_pest_scenario, FileManager &_file_manager,
	OutputFileWriter &_output_file_writer, PerformanceLog *_performance_log,
	RunManagerAbstract* _run_mgr_ptr) : pest_scenario(_pest_scenario), file_manager(_file_manager),
	output_file_writer(_output_file_writer), performance_log(_performance_log),
	run_mgr_ptr(_run_mgr_ptr), 
	constraints(_pest_scenario, &_file_manager, _output_file_writer, *_performance_log),
	jco(_file_manager,_output_file_writer)
{
	rand_gen = std::mt19937(pest_scenario.get_pestpp_options().get_random_seed());
	subset_rand_gen = std::mt19937(pest_scenario.get_pestpp_options().get_random_seed());
	dv.set_pest_scenario(&pest_scenario);
	oe.set_pest_scenario(&pest_scenario);
	dv.set_rand_gen(&rand_gen);
	oe.set_rand_gen(&rand_gen);
	
}

void SeqQuadProgram::throw_sqp_error(string message)
{
	performance_log->log_event("SeqQuadProgram error: " + message);
	cout << endl << "   ************   " << endl << "    SeqQuadProgram error: " << message << endl << endl;
	file_manager.rec_ofstream() << endl << "   ************   " << endl << "    SeqQuadProgram error: " << message << endl << endl;
	file_manager.close_file("rec");
	performance_log->~PerformanceLog();
	throw runtime_error("SeqQuadProgram error: " + message);
}

//SeqQuadProgram::apply_draw_mult()  // right place?
//{
//	PestppOptions* ppo = pest_scenario.get_pestpp_options_ptr();
//	float dm = ppo->get_sqp_dv_draw_mult(); // TODO add as ppo arg
//	cov * dm
//	return cov 
//}

bool SeqQuadProgram::initialize_dv(Covariance &cov)
{
	stringstream ss;
	int num_reals = pest_scenario.get_pestpp_options().get_sqp_num_reals();
	string dv_file = pest_scenario.get_pestpp_options().get_sqp_dv_en();
	bool drawn = false;
	if (dv_file.size() == 0)
	{
		//only draw for dv names
		Covariance dv_cov = cov.get(dv_names);
		message(1, "dv_cov:", dv_cov);  // tmp
		// apply draw_mult (++ arg?)
		// dv_cov = apply_draw_mult();
		ofstream& frec = file_manager.rec_ofstream();
		message(1, "drawing decision variable realizations: ", num_reals);
		map<string, double> par_means = pest_scenario.get_ext_file_double_map("parameter data external", "mean");
		Parameters draw_par = pest_scenario.get_ctl_parameters().get_subset(dv_names.begin(),dv_names.end());
		if (par_means.size() > 0)
		{
			
			frec << "Note: the following decision variables contain 'mean' value information that will be used in place of " << endl;
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
						frec << "Warning: 'mean' value for decision variable " << par_mean.first << " less than lower bound, using 'parval1'";
					}
					else if (par_mean.second > ub)
					{
						frec << "Warning: 'mean' value for decision variable " << par_mean.first << " greater than upper bound, using 'parval1'";
					}
					else
					{
						draw_par[par_mean.first] = par_mean.second;
						frec << par_mean.first << " " << par_mean.second << endl;
					}
					
				}

			}
		}
		dv.draw(num_reals, draw_par, dv_cov, performance_log, pest_scenario.get_pestpp_options().get_ies_verbose_level(), file_manager.rec_ofstream());
		drawn = true;
	}
	else
	{
		string par_ext = pest_utils::lower_cp(dv_file).substr(dv_file.size() - 3, dv_file.size());
		performance_log->log_event("processing par csv " + dv_file);
		if (par_ext.compare("csv") == 0)
		{
			message(1, "loading dv ensemble from csv file", dv_file);
			try
			{
				dv.from_csv(dv_file);
			}
			catch (const exception &e)
			{
				ss << "error processing dv csv file: " << e.what();
				throw_sqp_error(ss.str());
			}
			catch (...)
			{
				throw_sqp_error(string("error processing dv csv file"));
			}
		}
		else if ((par_ext.compare("jcb") == 0) || (par_ext.compare("jco") == 0))
		{
			message(1, "loading dv ensemble from binary file", dv_file);
			try
			{
				dv.from_binary(dv_file);
			}
			catch (const exception &e)
			{
				ss << "error processing dv jcb file: " << e.what();
				throw_sqp_error(ss.str());
			}
			catch (...)
			{
				throw_sqp_error(string("error processing dv jcb file"));
			}
		}
		else
		{
			ss << "unrecognized dv ensemble file extension " << par_ext << ", looking for csv, jcb, or jco";
			throw_sqp_error(ss.str());
		}

		dv.transform_ip(ParameterEnsemble::transStatus::NUM);
		
		if (pp_args.find("SQP_NUM_REALS") != pp_args.end())
		{
			int num_reals = pest_scenario.get_pestpp_options().get_sqp_num_reals();
			/*if (pest_scenario.get_pestpp_options().get_ies_include_base())
			{
				message(1, "Note: increasing num_reals by 1 to account for 'base' realization in existing par ensemble");
				num_reals++;
			}*/
			if (num_reals < dv.shape().first)
			{
				message(1,"ies_num_reals arg passed, truncated parameter ensemble to ",num_reals);
				vector<string> keep_names,real_names=dv.get_real_names();
				for (int i=0;i<num_reals;i++)
				{
					keep_names.push_back(real_names[i]);
				}
				dv.keep_rows(keep_names);
			}
		}
		

		//TODO: sqp version of this arg?
		if (pest_scenario.get_pestpp_options().get_ies_enforce_bounds())
		{
			if (pest_scenario.get_pestpp_options().get_ies_obs_restart_csv().size() > 0)
				message(1, "Warning: even though ies_enforce_bounds is true, a restart obs en was passed, so bounds will not be enforced on the initial par en");
			else
				dv.enforce_limits(performance_log, pest_scenario.get_pestpp_options().get_ies_enforce_chglim());
		}

	}

	if (dv_names.size() < pest_scenario.get_ctl_ordered_adj_par_names().size())
	{
		performance_log->log_event("filling non-decision-variable columns with control file values");
		Parameters ctl_num_pars = pest_scenario.get_ctl_parameters();
		pest_scenario.get_base_par_tran_seq().ctl2numeric_ip(ctl_num_pars);
		vector<string> ctl_adj_par_names = pest_scenario.get_ctl_ordered_adj_par_names();
		Eigen::VectorXd ctl_vals = ctl_num_pars.get_data_eigen_vec(ctl_adj_par_names);
		Eigen::MatrixXd temp(dv.shape().first, ctl_vals.size());
		for (int i = 0; i < temp.rows(); i++)
			temp.row(i) = ctl_vals;
		ParameterEnsemble dv_full(&pest_scenario, &rand_gen, temp, dv.get_real_names(), ctl_adj_par_names);
		dv_full.set_trans_status(ParameterEnsemble::transStatus::NUM);
		dv.update_var_map();
		for (auto d : dv.get_var_map())
		{
			Eigen::VectorXd col = dv.get_eigen_ptr()->col(d.second);
			dv_full.replace_col(d.first, col);
		}
		dv = dv_full;
	}

	return drawn;

}

void SeqQuadProgram::add_bases()
{
	//check that 'base' isn't already in ensemble
	vector<string> rnames = dv.get_real_names();
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
		dv.get_par_transform().active_ctl2numeric_ip(pars);
		vector<int> drop{ dv.shape().first - 1 };
		dv.drop_rows(drop);
		dv.append(base_name, pars);
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
			vector<string> prnames = dv.get_real_names();

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


//template<typename T, typename A>
//void SeqQuadProgram::message(int level, char* _message, vector<T, A> _extras)
//{
//	string s(_message);
//	message(level, s, _extras);
//}

//void SeqQuadProgram::message(int level, char* _message)
//{
//	string s(_message);
//	message(level, s);
//}

//template<typename T>
//void SeqQuadProgram::message(int level, char* _message, T extra)
//{
//	string s(_message);
//	message(level, s, extra);
//
//}

template<typename T, typename A>
void SeqQuadProgram::message(int level, const string &_message, vector<T, A> _extras, bool echo)
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

void SeqQuadProgram::message(int level, const string &_message)
{
	message(level, _message, vector<string>());
}

template<typename T>
void SeqQuadProgram::message(int level, const string &_message, T extra)
{
	stringstream ss;
	ss << _message << " " << extra;
	string s = ss.str();
	message(level, s);
}

void SeqQuadProgram::sanity_checks()
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
		warnings.push_back("'pestmode' == 'regularization', in pestpp-sqp, this has no meaning...");
		//throw_sqp_error("'pestmode' == 'regularization', please reset to 'estimation'");
	}
	else if (pest_scenario.get_control_info().pestmode == ControlInfo::PestMode::UNKNOWN)
	{
		warnings.push_back("unrecognized 'pestmode', using 'estimation'");
	}

	
	
	
	
	if ((ppo->get_sqp_num_reals() < warn_min_reals) && (par_csv.size() == 0))
	{
		ss.str("");
		ss << "ies_num_reals < " << warn_min_reals << ", this is prob too few";
		warnings.push_back(ss.str());
	}
	

	string how = pest_scenario.get_pestpp_options().get_ies_subset_how();
	if ((how != "FIRST") && (how != "LAST") && (how != "RANDOM") && (how != "PHI_BASED"))
	{
		ss.str("");
		ss << "'subset_how' is '" << how << "' but should be 'FIRST','LAST','RANDOM','PHI_BASED'";
		errors.push_back(ss.str());
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
		throw_sqp_error(string("sanity_check() found some problems - please review rec file"));
	}
	//cout << endl << endl;
}

void SeqQuadProgram::initialize_objfunc()
{
	//initialize the objective function
	obj_func_str = pest_scenario.get_pestpp_options().get_opt_obj_func();
	obj_sense = (pest_scenario.get_pestpp_options().get_opt_direction() == 1) ? "minimize" : "maximize";

	ofstream& f_rec = file_manager.rec_ofstream();


	//check if the obj_str is an observation
	use_obj_obs = false;
	if (pest_scenario.get_ctl_observations().find(obj_func_str) != pest_scenario.get_ctl_observations().end())
	{
		use_obj_obs = true;
		obj_obs = obj_func_str;
		//check
		vector<string> cnames = constraints.get_obs_constraint_names();
		set<string> names(cnames.begin(), cnames.end());
		if (names.find(obj_obs) != names.end())
		{
			throw_sqp_error("objective function obs is a constraint, #sad");
		}
		names.clear();
		cnames = constraints.get_nz_obs_names();
		names.insert(cnames.begin(), cnames.end());
		if (names.find(obj_obs) != names.end())
		{
			throw_sqp_error("objective function obs has non-zero weight and chance constraints are active");
		}
	}

	else
	{
		if (obj_func_str.size() == 0)
		{
			f_rec << " warning: no ++opt_objective_function-->forming a generic objective function (1.0 coef for each decision var)" << endl;
			for (auto& name : dv_names)
				obj_func_coef_map[name] = 1.0;
		}

		//or if it is a prior info equation
		else if (pest_scenario.get_prior_info().find(obj_func_str) != pest_scenario.get_prior_info().end())
		{
			obj_func_coef_map = pest_scenario.get_prior_info().get_pi_rec_ptr(obj_func_str).get_atom_factors();
		}
		else
		{
			//check if this obj_str is a filename
			ifstream if_obj(obj_func_str);
			if (!if_obj.good())
				throw_sqp_error("unrecognized ++opt_objective_function arg (tried file name, obs name, prior info name): " + obj_func_str);
			else
				obj_func_coef_map = pest_utils::read_twocol_ascii_to_map(obj_func_str);
		}


		//check that all obj_coefs are decsision vars
		vector<string> missing_vars;
		set<string> s_dv_names(dv_names.begin(), dv_names.end());
		for (auto& coef : obj_func_coef_map)
			if (s_dv_names.find(coef.first) == s_dv_names.end())
				missing_vars.push_back(coef.first);
		if (missing_vars.size() > 0)
		{
			stringstream ss;
			ss << "the following objective function components are not decision variables: ";
			for (auto m : missing_vars)
			{
				ss << m << ",";
				
			}
			throw_sqp_error(ss.str());
		}
	}
}


bool SeqQuadProgram::initialize_restart()
{
	stringstream ss;
	string obs_restart_csv = pest_scenario.get_pestpp_options().get_sqp_obs_restart_en();
	if (obs_restart_csv.size() == 0)
	{
		oe.initialize_without_noise(dv.shape().first);
		vector<string> real_names = dv.get_real_names();
		oe.set_real_names(real_names);
		return true;
	}
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
			throw_sqp_error(ss.str());
		}
		catch (...)
		{
			throw_sqp_error(string("error processing restart obs csv"));
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
			throw_sqp_error(ss.str());
		}
		catch (...)
		{
			throw_sqp_error(string("error processing restart obs binary file"));
		}
	}
	else
	{
		ss << "unrecognized restart obs ensemble extension " << obs_ext << ", looking for csv, jcb, or jco";
		throw_sqp_error(ss.str());
	}
	

	if (pp_args.find("SQP_NUM_REALS") != pp_args.end())
	{
		int num_reals = pest_scenario.get_pestpp_options().get_ies_num_reals();
		/*if (pest_scenario.get_pestpp_options().get_ies_include_base())
		{
			message(1, "Note: increasing num_reals by 1 to account for 'base' realization in existing obs restart ensemble");
			num_reals++;
		}*/
		if (num_reals < oe.shape().first)
		{
			message(1, "sqp_num_reals arg passed, truncated restart obs ensemble to ", num_reals);
			vector<string> keep_names, real_names = oe.get_real_names();
			for (int i = 0; i<num_reals; i++)
			{
				keep_names.push_back(real_names[i]);
			}
			oe.keep_rows(keep_names);
		}
	}

	
	if (oe.shape().first != dv.shape().first)
	{
		//check if all oe names are found in par en, if so, we can reorder and proceed.  otherwise, die
		vector<string> missing;
		vector<string> oe_real_names = oe.get_real_names();
		vector<string> pe_real_names = dv.get_real_names();
		for (auto &oname : oe_real_names)
		{
			if (find(pe_real_names.begin(), pe_real_names.end(), oname) == pe_real_names.end())
				missing.push_back(oname);
		}

		if (missing.size() > 0)
		{
			ss << "number of reals differ between restart obs en (" << oe.shape().first << ") and par en (" << dv.shape().first << ")";
			ss << " and realization names could not be aligned:";
			for (auto &m : missing)
				ss << m << ",";
			throw_sqp_error(ss.str());
		}

		message(2, "reordering dv to align with restart obs en, num reals: ", oe_real_names.size());
		try
		{
			dv.reorder(oe_real_names, vector<string>());
		}
		catch (exception &e)
		{
			ss << "error reordering dv with restart oe:" << e.what();
			throw_sqp_error(ss.str());
		}
		catch (...)
		{
			throw_sqp_error(string("error reordering dv with restart oe"));
		}

	}

	return false;
}


void SeqQuadProgram::initialize_parcov()
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


void SeqQuadProgram::initialize_obscov()
{
	message(1, "initializing observation noise covariance matrix");
	string obscov_filename = pest_scenario.get_pestpp_options().get_obscov_filename();

	string how = obscov.try_from(pest_scenario, file_manager, false, true);
	message(1, "obscov loaded ", how);
	
}


void SeqQuadProgram::initialize()
{	
	message(0, "initializing");
	pp_args = pest_scenario.get_pestpp_options().get_passed_args();

	iter = 1;

	act_obs_names = pest_scenario.get_ctl_ordered_nz_obs_names();
	act_par_names = pest_scenario.get_ctl_ordered_adj_par_names();

	stringstream ss;
	//set some defaults
	PestppOptions* ppo = pest_scenario.get_pestpp_options_ptr();

	if (pp_args.find("PAR_SIGMA_RANGE") == pp_args.end())
	{
		message(1, "resetting par_sigma_range to 20.0");
		ppo->set_par_sigma_range(20.0);
	}


	//reset the par bound PI augmentation since that option is just for simplex
	ppo->set_opt_include_bnd_pi(false);

	//process dec var args
	vector<string> dec_var_groups = ppo->get_opt_dec_var_groups();
	if (dec_var_groups.size() != 0)
	{
		//first make sure all the groups are actually listed in the control file
		vector<string> missing;
		vector<string> pst_groups = pest_scenario.get_ctl_ordered_par_group_names();
		vector<string>::iterator end = pst_groups.end();
		vector<string>::iterator start = pst_groups.begin();
		for (auto grp : dec_var_groups)
			if (find(start, end, grp) == end)
				missing.push_back(grp);
		if (missing.size() > 0)
		{
			ss.str("");
			ss << "the following ++opt_dec_var_groups were not found: ";
			for (auto m : missing)
				ss << m << ",";
			throw_sqp_error(ss.str());
		}


		//find the parameter in the dec var groups
		ParameterGroupInfo pinfo = pest_scenario.get_base_group_info();
		string group;
		end = dec_var_groups.end();
		start = dec_var_groups.begin();
		for (auto& par_name : pest_scenario.get_ctl_ordered_par_names())
		{
			group = pinfo.get_group_name(par_name);
			if (find(start, end, group) != end)
			{
				dv_names.push_back(par_name);

			}
		}

		if (dv_names.size() == 0)
		{
			ss.str("");
			ss << "no adjustable decision variables found in supplied dec var groups : ";
			for (auto g : dec_var_groups)
			{
				ss << g << ",";
			}
			throw_sqp_error(ss.str());
		}
		ss.str("");
		ss << "'opt_dec_var_groups' passed, using " << dv_names.size() << " adjustable parameters as decision variables";
		message(2, ss.str());
		ofstream& frec = file_manager.rec_ofstream();
		frec << "decision variables:" << endl;
		int icol = 0;
		for (auto dv_name : dv_names)
		{
			frec << dv_name << " ";
			icol++;
			if (icol == 10)
			{
				frec << endl;
				icol = 0;
			}
		}
	}

	//otherwise, just use all adjustable parameters as dec vars
	else
	{
		message(2, "using all adjustable parameters as decision variables: ", act_par_names.size());
		dv_names = act_par_names;
	}

	constraints.initialize(dv_names, numeric_limits<double>::max());
	constraints.initial_report();
	//some risk-based stuff here
	string chance_points = ppo->get_opt_chance_points();
	if (chance_points == "ALL")
	{
		//evaluate the chance constraints at every individual, very costly, but most robust
		//throw_sqp_error("'opt_chance_points' == 'all' not implemented");
		chancepoints = chancePoints::ALL;
	}

	else if (chance_points == "SINGLE")
	{
		//evaluate the chance constraints only at the population member nearest the optimal tradeoff.
		//much cheaper, but assumes linear coupling
		chancepoints = chancePoints::SINGLE;

	}
	else
	{
		ss.str("");
		ss << "unrecognized 'opt_chance_points' value :" << chance_points << ", should be 'all' or 'single'";
		throw_sqp_error(ss.str());
	}

	iter = 0;

	if (pest_scenario.get_control_info().noptmax == 0)
	{
		message(0, "'noptmax'=0, running control file parameter values and quitting");
		
		Parameters pars = pest_scenario.get_ctl_parameters();
		ParamTransformSeq pts = dv.get_par_transform();

		ParameterEnsemble _pe(&pest_scenario, &rand_gen);
		_pe.reserve(vector<string>(), pest_scenario.get_ctl_ordered_par_names());
		_pe.set_trans_status(ParameterEnsemble::transStatus::CTL);
		_pe.append("BASE", pars);
		string par_csv = file_manager.get_base_filename() + ".par.csv";
		//message(1, "saving parameter values to ", par_csv);
		//_pe.to_csv(par_csv);
		dv_base = _pe;
		dv_base.reorder(vector<string>(), act_par_names);
		ObservationEnsemble _oe(&pest_scenario, &rand_gen);
		_oe.reserve(vector<string>(), pest_scenario.get_ctl_ordered_obs_names());
		_oe.append("BASE", pest_scenario.get_ctl_observations());
		oe_base = _oe;
		oe_base.reorder(vector<string>(), act_obs_names);
		//initialize the phi handler
		ph = L2PhiHandler(&pest_scenario, &file_manager, &oe_base, &dv_base, &parcov);
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
			throw_sqp_error("control file parameter value run failed");
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

	
	//ofstream &frec = file_manager.rec_ofstream();
	last_best_mean = 1.0E+30;
	last_best_std = 1.0e+30;
	
	warn_min_reals = 10;
	error_min_reals = 2;
	
	vector<double> scale_facs = pest_scenario.get_pestpp_options().get_lambda_scale_vec();
	message(1, "using scaling factors: ", scale_facs);
	
	message(1, "max run fail: ", ppo->get_max_run_fail());

	//TODO: update sanity checks for SQP context
	//check that if using fd, chance points == single
	sanity_checks();

	bool echo = false;
	if (verbose_level > 1)
		echo = true;

	initialize_parcov();
	//I dont think SQP needs an obscov? 
	//initialize_obscov();


	

	//todo: add the FD approach and initialize it here
	use_ensemble_grad = false;
	if (ppo->get_sqp_num_reals() > 0)
		use_ensemble_grad = true;
	if (use_ensemble_grad)
	{
		prep_4_ensemble_grad();
	}
	else
	{
		prep_4_fd_grad();
	}
		
	//these will be the ones we track...
	current_ctl_dv_values = pest_scenario.get_ctl_parameters();
	current_obs = pest_scenario.get_ctl_observations();
	

	message(0, "initialization complete");
}



void SeqQuadProgram::prep_4_fd_grad()
{
	stringstream ss;
	message(1, "using finite-difference approximation to gradient (Jacobian)");
	string base_jco = pest_scenario.get_pestpp_options().get_basejac_filename();
	if (base_jco.size() > 0)
	{
		message(1, "loading existing base jacobian " + base_jco);
		jco.read(base_jco);
		//todo: error trapping to make sure all the needed rows and cols are found
		vector<string> vnames = jco.get_base_numeric_par_names();
		set<string> snames(vnames.begin(), vnames.end());
		vnames.clear();
		for (auto& dv_name : dv_names)
			if (snames.find(dv_name) == snames.end())
				vnames.push_back(dv_name);
		if (vnames.size() > 0)
		{
			ss.str("");
			ss << "existing jacobian missing the following decision variables:" << endl;
			for (auto m : vnames)
				ss << vnames << endl;
			throw_sqp_error(ss.str());
		}
		snames.clear();
		vnames = constraints.get_obs_constraint_names();
		snames.insert(vnames.begin(), vnames.end());
		vnames.clear();
		for (auto name : jco.get_sim_obs_names())
			if (snames.find(name) == snames.end())
				vnames.push_back(name);

		if (vnames.size() > 0)
		{
			ss.str("");
			ss << "existing jacobian missing the following obs constraints:" << endl;
			for (auto m : vnames)
				ss << vnames << endl;
			throw_sqp_error(ss.str());
		}
	}
	else
	{
		
		//todo: handle hotstart_resfile here...
		bool init_obs = true;

		run_jacobian(current_ctl_dv_values, current_obs, init_obs);
		

	}
}

void SeqQuadProgram::run_jacobian(Parameters& _current_ctl_dv_vals, Observations& _current_obs, bool init_obs)
{
	stringstream ss;
	ParamTransformSeq par_trans = pest_scenario.get_base_par_tran_seq();
	ParameterGroupInfo pgi = pest_scenario.get_base_group_info();
	Parameters current_pars = pest_scenario.get_ctl_parameters();
	PriorInformation pi = pest_scenario.get_prior_info();
	current_pars.update_without_clear(dv_names,_current_ctl_dv_vals.get_data_eigen_vec(dv_names));

	set<string> out_of_bounds;
	ss.str("");
	ss << "queuing " << dv_names.size() << " finite difference runs";
	message(2, ss.str());
	bool success = jco.build_runs(current_pars, _current_obs, dv_names, par_trans,
		pest_scenario.get_base_group_info(), pest_scenario.get_ctl_parameter_info(),
		*run_mgr_ptr, out_of_bounds, false, init_obs);
	if (!success)
		throw_sqp_error("error building jacobian runs for FD grad");
	//todo: think about freezind dec vars that go out of bounds? - yuck!
	if (out_of_bounds.size() > 0)
	{
		ss.str("");
		ss << "the following decision variable are out of bounds: " << endl;
		for (auto& o : out_of_bounds)
			ss << o << ",";
		throw_sqp_error(ss.str());
	}
	queue_chance_runs();
	message(2, "starting finite difference gradient pertubation runs");
	jco.make_runs(*run_mgr_ptr);
	
	success = jco.process_runs(par_trans,pgi,*run_mgr_ptr,pi,false,false);
	if (!success)
	{
		throw_sqp_error("error processing finite difference gradient pertubation runs");
	}
	constraints.process_runs(run_mgr_ptr, iter);
	if (init_obs)
	{
		run_mgr_ptr->get_run(0, current_pars, _current_obs, false);

	}

}

void SeqQuadProgram::prep_4_ensemble_grad()
{
	stringstream ss;
	message(1, "using ensemble approximation to gradient (EnOpt)");
	
	//I think a bad phi option has use in SQP?
	/*double bad_phi = pest_scenario.get_pestpp_options().get_ies_bad_phi();
	if (bad_phi < 1.0e+30)
		message(1, "using bad_phi: ", bad_phi);*/

	int num_reals = pest_scenario.get_pestpp_options().get_sqp_num_reals();

	dv_drawn = initialize_dv(parcov);

	oe_drawn = initialize_restart();

	try
	{
		dv.check_for_dups();
	}
	catch (const exception& e)
	{
		string message = e.what();
		throw_sqp_error("error in dv ensemble: " + message);
	}

	try
	{
		oe.check_for_dups();
	}
	catch (const exception& e)
	{
		string message = e.what();
		throw_sqp_error("error in observation ensemble: " + message);
	}

	if (dv.shape().first != oe.shape().first)
	{
		//the special case where par en < obs en and all par reals are found in obs en...

		if (dv.shape().first < oe.shape().first)
		{
			vector<string> oe_names = oe.get_real_names();
			set<string> oset(oe_names.begin(), oe_names.end());
			vector<string> missing;
			for (auto n : dv.get_real_names())
				if (oset.find(n) == oset.end())
				{
					missing.push_back(n);
				}
			if (missing.size() == 0)
			{
				ss.str("");
				ss << "dv en has " << dv.shape().first << " realizations, compared to " << oe.shape().first << " obs realizations";
				message(1, ss.str());
				message(1, " the realization names are compatible");
				message(1, "re-indexing obs en to align with dv en...");

				oe.reorder(dv.get_real_names(), vector<string>());
			}
			else
			{
				ss.str("");
				ss << "the following dv en real names were not found in the obs en: ";
				for (auto m : missing)
				{
					ss << m << ",";
				}
				throw_sqp_error(ss.str());
			}
		}
		else
		{
			ss.str("");
			ss << "dv ensemble rows (" << dv.shape().first << ") not equal to observation ensemble rows (" << oe.shape().first << ")";
			throw_sqp_error(ss.str());
		}
	}


	//need this here for Am calcs...
	//message(1, "transforming parameter ensemble to numeric");
	dv.transform_ip(ParameterEnsemble::transStatus::NUM);


	//TODO: think about what adding the base would do for SQP?
	/*if (pest_scenario.get_pestpp_options().get_ies_include_base())
		if (pp_args.find("SQP_RESTART_OBS_EN") != pp_args.end())
		{
			message(1, "Warning: even though `sqp_include_base` is true, you passed a restart obs en, not adding 'base' realization...");
		}
		else
			add_bases();*/


			//now we check to see if we need to try to align the par and obs en
			//this would only be needed if either of these were not drawn
	if (!dv_drawn || !oe_drawn)
	{
		bool aligned = dv.try_align_other_rows(performance_log, oe);
		if (aligned)
		{
			message(2, "observation ensemble reordered to align rows with dv ensemble");
		}
	}

	//just check to see if common real names are found but are not in the same location
	map<string, int> pe_map = dv.get_real_map(), oe_map = oe.get_real_map();
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
		message(1, "WARNING: common realization names shared between the dv and observation ensembles but they are not in the same row locations, see .rec file for listing");
		ofstream& frec = file_manager.rec_ofstream();
		frec << endl << "WARNING: the following " << misaligned.size() << " realization names are shared between the dv and observation ensembles but they are not in the same row locations:" << endl;
		for (auto ma : misaligned)
			frec << ma << endl;
	}



	message(2, "checking for denormal values in dv");
	dv.check_for_normal("initial transformed dv ensemble");
	ss.str("");
	//TODO: setup an sqp save bin flag?  or piggy back?
	if (pest_scenario.get_pestpp_options().get_ies_save_binary())
	{
		ss << file_manager.get_base_filename() << ".0.par.jcb";
		dv.to_binary(ss.str());
	}
	else
	{
		ss << file_manager.get_base_filename() << ".0.par.csv";
		dv.to_csv(ss.str());
	}
	message(1, "saved initial dv ensemble to ", ss.str());
	message(2, "checking for denormal values in base oe");
	oe.check_for_normal("observation ensemble");
	ss.str("");


	/*if (pest_scenario.get_pestpp_options().get_ies_save_binary())
	{
		ss << file_manager.get_base_filename() << ".obs.jcb";
		oe.to_binary(ss.str());
	}
	else
	{
		ss << file_manager.get_base_filename() << ".obs.csv";
		oe.to_csv(ss.str());
	}
	message(1, "saved initial observation ensemble to ", ss.str());*/
	message(1, "centering on ensemble mean vector");

	if (pest_scenario.get_control_info().noptmax == -2)
	{
		message(0, "'noptmax'=-2, running mean dv ensemble values and quitting");
		message(1, "calculating mean dv values");
		Parameters pars;
		vector<double> mv = dv.get_mean_stl_var_vector();
		pars.update(dv.get_var_names(), dv.get_mean_stl_var_vector());
		ParamTransformSeq pts = dv.get_par_transform();

		ParameterEnsemble _pe(&pest_scenario, &rand_gen);
		_pe.reserve(vector<string>(), dv.get_var_names());
		_pe.set_trans_status(dv.get_trans_status());
		_pe.append("mean", pars);
		string par_csv = file_manager.get_base_filename() + ".mean.par.csv";
		message(1, "saving mean dv values to ", par_csv);
		_pe.to_csv(par_csv);
		dv_base = _pe;
		dv_base.reorder(vector<string>(), act_par_names);
		ObservationEnsemble _oe(&pest_scenario, &rand_gen);
		_oe.reserve(vector<string>(), oe.get_var_names());
		_oe.append("mean", pest_scenario.get_ctl_observations());
		oe_base = _oe;
		oe_base.reorder(vector<string>(), act_obs_names);
		//initialize the phi handler
		ph = L2PhiHandler(&pest_scenario, &file_manager, &oe_base, &dv_base, &parcov);
		if (ph.get_lt_obs_names().size() > 0)
		{
			message(1, "less_than inequality defined for observations: ", ph.get_lt_obs_names().size());
		}
		if (ph.get_gt_obs_names().size())
		{
			message(1, "greater_than inequality defined for observations: ", ph.get_gt_obs_names().size());
		}
		message(1, "running mean dv values");

		vector<int> failed_idxs = run_ensemble(_pe, _oe);
		if (failed_idxs.size() != 0)
		{
			message(0, "mean dv value run failed...bummer");
			return;
		}
		string obs_csv = file_manager.get_base_filename() + ".mean.obs.csv";
		message(1, "saving results from mean dv value run to ", obs_csv);
		_oe.to_csv(obs_csv);


		//TODO: rather than report l2 phi, report actual obj func and feas status
		ph.update(_oe, _pe);
		message(0, "mean parameter phi report:");
		ph.report(true);
		ph.write(0, 1);

		return;
	}

	oe_org_real_names = oe.get_real_names();
	pe_org_real_names = dv.get_real_names();
	string obs_restart_csv = pest_scenario.get_pestpp_options().get_ies_obs_restart_csv();
	string par_restart_csv = pest_scenario.get_pestpp_options().get_ies_par_restart_csv();

	//TODO: I think the base_oe should just be a "no noise" obs ensemble?
	//or do we even need it?  Or can we use this to do a chen-oliver style 
	//robust opt?
	oe_base = oe; //copy

	//reorder this for later...
	oe_base.reorder(vector<string>(), act_obs_names);

	dv_base = dv; //copy
	//reorder this for later
	dv_base.reorder(vector<string>(), act_par_names);

	//no restart
	if (oe_drawn)
	{
		performance_log->log_event("running initial ensemble");
		message(1, "running initial ensemble of size", oe.shape().first);
		vector<int> failed = run_ensemble(dv, oe);
		if (dv.shape().first == 0)
			throw_sqp_error("all realizations failed during initial evaluation");

		dv.transform_ip(ParameterEnsemble::transStatus::NUM);
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
	save_base_real_par_rei(pest_scenario, dv, oe, output_file_writer, file_manager, iter);
	save_base_real_par_rei(pest_scenario, dv, oe, output_file_writer, file_manager, -1);


	performance_log->log_event("calc pre-drop phi");
	//todo: I think we probably want to use the opt_objective_function arg
	//to define the objective function.  Support and obs or a pi.  
	//what about an external file like pestpp-opt with just coefficients?


	//initialize the phi handler
	//todo: can we also use the l2 norm of residuals to ident 
	//"bad" realizations?

	/*ph = L2PhiHandler(&pest_scenario, &file_manager, &oe_base, &dv_base, &parcov);

	if (ph.get_lt_obs_names().size() > 0)
	{
		message(1, "less_than inequality defined for observations: ", ph.get_lt_obs_names().size());
	}
	if (ph.get_gt_obs_names().size())
	{
		message(1, "greater_than inequality defined for observations: ", ph.get_gt_obs_names().size());
	}

	ph.update(oe, dv);
	message(0, "pre-drop initial phi summary");
	ph.report(true);
	drop_bad_phi(dv, oe);
	if (oe.shape().first == 0)
	{
		throw_sqp_error(string("all realizations dropped as 'bad'"));
	}
	if (oe.shape().first <= error_min_reals)
	{
		message(0, "too few active realizations:", oe.shape().first);
		message(1, "need at least ", error_min_reals);
		throw_sqp_error(string("too few active realizations, cannot continue"));
	}
	if (oe.shape().first < warn_min_reals)
	{
		ss.str("");
		ss << "WARNING: less than " << warn_min_reals << " active realizations...might not be enough";
		string s = ss.str();
		message(0, s);
	}*/

	pcs = ParChangeSummarizer(&dv_base, &file_manager, &output_file_writer);

}

void SeqQuadProgram::drop_bad_phi(ParameterEnsemble &_pe, ObservationEnsemble &_oe, bool is_subset)
{
	//don't use this assert because _pe maybe full size, but _oe might be subset size
	if (!is_subset)
		if (_pe.shape().first != _oe.shape().first)
			throw_sqp_error("SeqQuadProgram::drop_bad_phi() error: _pe != _oe and not subset");
		
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
		full_pnames = dv.get_real_names();
		vector<string> pdrop, odrop;
		for (auto i : idxs)
		{
			oname = obs_real_names[i];
			
			/*if (is_subset)
			{
				pidx = find(full_onames.begin(), full_onames.end(), oname) - full_onames.begin();
				if (find(subset_idxs.begin(), subset_idxs.end(), pidx) == subset_idxs.end())
				{
					ss.str("");
					ss << "drop_bad_phi() error: idx " << pidx << " not found in subset_idxs";
					throw_sqp_error(ss.str());
				}
				pname = full_pnames[pidx];
			}*/
			//else
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
			throw_sqp_error(ss.str());
		}
		catch (...)
		{
			throw_sqp_error(string("drop_bad_phi() error"));
		}
	}
}

void SeqQuadProgram::save_mat(string prefix, Eigen::MatrixXd &mat)
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


void SeqQuadProgram::iterate_2_solution()
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
		ph.update(oe,dv);
		last_best_mean = ph.get_mean(L2PhiHandler::phiType::COMPOSITE);
		last_best_std = ph.get_std(L2PhiHandler::phiType::COMPOSITE);
		ph.report(true);
		ph.write(iter, run_mgr_ptr->get_total_runs());
		if (pest_scenario.get_pestpp_options().get_ies_save_rescov())
			ph.save_residual_cov(oe,iter);
		ss.str("");
		ss << file_manager.get_base_filename() << "." << iter << ".pcs.csv";
		pcs.summarize(dv,iter,ss.str());
			
		if (should_terminate())
			break;
	}
}

bool SeqQuadProgram::should_terminate()
{
	//todo: work out some SQP-based criteria here...
	//double phiredstp = pest_scenario.get_control_info().phiredstp;
	//int nphistp = pest_scenario.get_control_info().nphistp;
	//int nphinored = pest_scenario.get_control_info().nphinored;
	//bool phiredstp_sat = false, nphinored_sat = false, consec_sat = false;
	//double phi, ratio;
	//int count = 0;
	//int nphired = 0;

	//message(0, "phi-based termination criteria check");
	//message(1, "phiredstp: ", phiredstp);
	//message(1, "nphistp: ", nphistp);
	//message(1, "nphinored (also used for consecutive bad lambda cycles): ", nphinored);
	//if (best_mean_phis.size() > 0)
	//{
	//	vector<double>::iterator idx = min_element(best_mean_phis.begin(), best_mean_phis.end());
	//	nphired = (best_mean_phis.end() - idx) - 1;
	//	best_phi_yet = best_mean_phis[idx - best_mean_phis.begin()];// *pest_scenario.get_pestpp_options().get_ies_accept_phi_fac();
	//	message(1, "best mean phi sequence: ", best_mean_phis);
	//	message(1, "best phi yet: ", best_phi_yet);
	//}
	//message(1, "number of consecutive bad lambda testing cycles: ", consec_bad_lambda_cycles);
	//if (consec_bad_lambda_cycles >= nphinored)
	//{
	//	message(1, "number of consecutive bad lambda testing cycles > nphinored");
	//	consec_sat = true;
	//}

	//for (auto &phi : best_mean_phis)
	//{
	//	ratio = (phi - best_phi_yet) / phi;
	//	if (ratio <= phiredstp)
	//		count++;
	//}
	//message(1, "number of iterations satisfying phiredstp criteria: ", count);
	//if (count >= nphistp)
	//{
	//	message(1, "number iterations satisfying phiredstp criteria > nphistp");
	//	phiredstp_sat = true;
	//}

	//message(1, "number of iterations since best yet mean phi: ", nphired);
	//if (nphired >= nphinored)
	//{
	//	message(1, "number of iterations since best yet mean phi > nphinored");
	//	nphinored_sat = true;
	//}

	//if ((nphinored_sat) || (phiredstp_sat) || (consec_sat))
	//{
	//	message(1, "phi-based termination criteria satisfied, all done");
	//	return true;
	//}
	return false;
}







void SeqQuadProgram::update_reals_by_phi(ParameterEnsemble &_pe, ObservationEnsemble &_oe)
{

	vector<string> oe_names = _oe.get_real_names();
	vector<string> pe_names = _pe.get_real_names();
	vector<string> oe_base_names = oe.get_real_names();
	vector<string> pe_base_names = dv.get_real_names();

	//if (pe_names.size() != oe_base_names.size())
	//	throw runtime_error("SeqQuadProgram::update_reals_by_phi() error: pe_names != oe_base_names");
	map<string, int> oe_name_to_idx;
	map<int,string> pe_idx_to_name;

	for (int i = 0; i < oe_base_names.size(); i++)
		oe_name_to_idx[oe_base_names[i]] = i;
	
	for (int i = 0; i < pe_base_names.size(); i++)
		pe_idx_to_name[i] = pe_base_names[i];
	//store map of current phi values
	ph.update(oe, dv);
	L2PhiHandler::phiType pt = L2PhiHandler::phiType::COMPOSITE;
	map<string, double> phi_map = ph.get_phi_map(pt);
	map<string, double> cur_phi_map;
	for (auto p : phi_map)
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
		new_phi = phi_map.at(oname);
		cur_phi = cur_phi_map.at(oname);
		if (new_phi < cur_phi * acc_fac)
		{
			//pname = pe_names[i];
			//pname = pe_names[oe_name_to_idx[oname]];
			pname = pe_idx_to_name[oe_name_to_idx[oname]];
			if (find(pe_names.begin(), pe_names.end(), pname) == pe_names.end())
				throw runtime_error("SeqQuadProgram::update_reals_by_phi() error: pname not in pe_names: " + pname);
			ss.str("");
			ss << "updating dv:oe real =" << pname << ":" << oname << ", current phi: new phi  =" << cur_phi << ":" << new_phi;
			message(3, ss.str());
			
			real = _pe.get_real_vector(pname);
			dv.update_real_ip(pname, real);
			real = _oe.get_real_vector(oname);
			oe.update_real_ip(oname, real);
		}
	}
	ph.update(oe, dv);

}

Eigen::VectorXd SeqQuadProgram::calc_gradient_vector(const Parameters& _current_dv_values)
{
	Eigen::VectorXd grad(dv_names.size());
	string center_on = pest_scenario.get_pestpp_options().get_ies_center_on();
	if (use_ensemble_grad)
	{
		//ensemble stuff here
		if (use_obj_obs)
		{
			
			//todo: need to make sqp specific or refactor ies arg names
			//this is really just a vector of obj function anoms around the center_on point
			Eigen::MatrixXd obj_anoms = oe.get_eigen_anomalies(vector<string>(), vector<string>{obj_func_str}, center_on);
			//this is a matrix of dv anoms around the center_on point
			Eigen::MatrixXd dv_anoms = oe.get_eigen_anomalies(vector<string>(), dv_names, center_on);
			
			// calc the dec var ensemble (approx) covariance vector (e.g., eq (8) of Dehdari and Oliver 2012 SPE)
			for (int i = 0; i < dv_anoms.rows(); i++)
			{
				//assuming obj_anoms is being treated as a matrix - could swap it to a vector above
				dv_anoms.row(i).array() *= obj_anoms.row(i)[0];
			}
			Eigen::VectorXd dv_cross_cov = dv_anoms.rowwise().sum();
			//assumes that parcov has been pseudo inv ip'd before now...probably best to pass in a ref to an already inv'd parcov since we 
			//will be using adaptive cov matrix
			dv_cross_cov = *parcov.e_ptr() * dv_cross_cov;

			//parcov.pseudo_inv_ip()
			// calc the dec var-phi ensemble (approx) cross-covariance vector (e.g., eq (9) of Dehdari and Oliver 2012 SPE)

			// calc pseudo inv of ensemble dec var covariance vector

			// calc gradient: (outer) product of pseudo inv of dec var cov vector and dec var-phi cross-cov vector
			
			//todo: actually calc the grad!

			throw_sqp_error("obs-based obj for ensembles not implemented");

			//todo: localize the gradient

		}
		//pi base obj, need representative dv values using the "center_on" arg
		//represent the mean/median/base - that is, derived from the "center_on" arg
		//todo: for now, just using mean dv values
		else
		{
			//if not center_on arg, use the mean dv values
			if (center_on.size() == 0)
			{
				pair<map<string, double>, map<string, double>> mm = dv.get_moment_maps();
				for (int i = 0; i < dv_names.size(); i++)
				{
					grad[i] = obj_func_coef_map[dv_names[i]] * mm.first[dv_names[i]];
				}
			}
			else
			{
				
					grad = dv.get_real_vector(pest_scenario.get_pestpp_options().get_ies_center_on());
			}
		}
	}
	else
	{
		//obs-based obj
		if (use_obj_obs)
		{
			//just a jco row
			vector<string> obj_name_vec{obj_func_str};
			grad = jco.get_matrix(obj_name_vec, dv_names);
		}
		//pi based obj
		else
		{
			for (int i=0;i<dv_names.size();i++)
			{
				grad[i] = obj_func_coef_map[dv_names[i]] * _current_dv_values.get_rec(dv_names[i]);
			}
		}
	}
	return grad;
}


Parameters SeqQuadProgram::fancy_solve_routine(double scale_val, const Parameters& _current_dv_values)
{
	Parameters candidate;
	Eigen::VectorXd grad = calc_gradient_vector(_current_dv_values);
	return candidate;
}

bool SeqQuadProgram::solve_new()
{
	stringstream ss;
	ofstream& frec = file_manager.rec_ofstream();
	if (dv.shape().first <= error_min_reals)
	{
		message(0, "too few active realizations:", oe.shape().first);
		message(1, "need at least ", error_min_reals);
		throw_sqp_error(string("too few active realizations, cannot continue"));
	}
	if (dv.shape().first < warn_min_reals)
	{
		ss.str("");
		ss << "WARNING: less than " << warn_min_reals << " active realizations...might not be enough";
		string s = ss.str();
		message(1, s);
	}

	Parameters _current_num_dv_values = current_ctl_dv_values;  // make copy
	ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();
	pts.ctl2numeric_ip(_current_num_dv_values);  // covert to numeric format
	_current_num_dv_values = _current_num_dv_values.get_subset(dv_names.begin(), dv_names.end());  // just dec var operation
	if (use_ensemble_grad)
	{
		dv.transform_ip(ParameterEnsemble::transStatus::NUM);
		vector<double> mean_vec = dv.get_mean_stl_var_vector();
		_current_num_dv_values.update(dv_names, mean_vec);
	}

	// compute gradient
	//if (use_ensemble_grad)
	//{
	//	if (LBFGS) &(num_it > 2)& (constraints False); // constraint False as need phi_grad for Lagrangian
	//	{
	//		ss.str("");
	//		ss << "(re)use grad from Wolfe testing during upgrade evaluations last iteration";
	//		string s = ss.str();
	//		message(1, s);
	//		throw_sqp_error("TODO");
	//	}
	//	else
	//	{
	//		ss.str("");
	//		ss << "compute dec var en cov vector and dec var en-phi en cross-covariance vector";
	//		string s = ss.str();
	//		message(1, s);
	//		throw_sqp_error("TODO");

	//		// CMA implementation to go here

	//		ss.str("");
	//		ss << "compute gradient vector through product of pseudo inv of dec var and cross-cov";
	//		string s = ss.str();
	//		message(1, s);
	//		throw_sqp_error("TODO");

	//		ss.str("");
	//		ss << "compute ensemble approx to (active) constraint jacobian";
	//		string s = ss.str();
	//		message(1, s);
	//		throw_sqp_error("TODO");
	//	}
	//}
	//else // finite differences for gradient
	//{
	//	if (LBFGS) &(num_it > 2);
	//	{
	//		ss.str("");
	//		ss << "(re)use grad from Wolfe testing during upgrade evaluations last iteration";
	//		string s = ss.str();
	//		message(1, s);
	//		throw_sqp_error("TODO");
	//	}
	//	else // standard BFGS
	//	{
	//		ss.str("");
	//		ss << "populate jco and partition into phi grad vector and (active) constraint jco";
	//		string s = ss.str();
	//		message(1, s);
	//		throw_sqp_error("TODO");
	//	}
	//}

	//// compute search direction
	//if (constraints is True) &(len(working_set) > 0);
	//{
	//	ss.str("");
	//	ss << "solve (E)QP sub-problem via active set method for search direction";
	//	string s = ss.str();
	//	message(1, s);
	//	//self.search_d, self.lagrang_mults = self._solve_eqp(qp_solve_method=self.qp_solve_method)
	//	//search direction tests contained therein
	//	throw_sqp_error("TODO");
	//}
	//else if ((constraints is False) | ((constraints is True) & (len(self.working_set) == 0)));
	//{
	//	if (LBFGS);
	//	{
	//		ss.str("");
	//		ss << "use limited-memory BFGS algorithm to solve for search direction in unconstrained case";
	//		string s = ss.str();
	//		message(1, s);
	//		//self.search_d = self._LBFGS_hess_update(memory = memory, self_scale = hess_self_scaling)
	//		throw_sqp_error("TODO");
	//	}
	//	else  //BFGS
	//	{
	//		ss.str("");
	//		ss << "use standard BFGS algorithm to solve for search direction in unconstrained case";
	//		string s = ss.str();
	//		message(1, s);
	//		//self.search_d = self._LBFGS_hess_update(memory = memory, self_scale = hess_self_scaling)
	//		throw_sqp_error("TODO");
	//	}

	//	// now first attempt at updating active set
	//	if (constraints is True);  //irrespective of shape of working set
	//	{
	//		ss.str("");
	//		ss << "updating active set of constraints; first attempt";
	//		string s = ss.str();
	//		message(1, s);
	//		//alpha, next_it = self._active_set_method(first_pass = True)
	//		throw_sqp_error("TODO");
	//	}
	//
	//// undertake some search direction-related tests here, e.g., point down-hill
	//// and rectify if bad!
	//
	//}

	// attempt update of active set again, now where constraints are active 
	//if ((constraints is True) & (len(working_set) > 0));
	//{
	//	ss.str("");
	//	ss << "updating active set of constraints; first attempt";
	//	string s = ss.str();
	//	message(1, s);
	//	//alpha, next_it = self._active_set_method(first_pass = True)
	//	throw_sqp_error("TODO");
	//}

	Parameters dv_candidate;
	ParameterEnsemble dv_candidates(&pest_scenario,&rand_gen);
	dv_candidates.set_trans_status(ParameterEnsemble::transStatus::NUM);
	
	
	//backtracking/line search factors
	//TODO: make this a ++ arg or tunable or something clever
	vector<double> scale_vals{ 0.1,0.5,1.0 };
	vector<string> real_names;
	
	for (auto sv : scale_vals)
	{
		ss.str("");
		ss << "dv_cand_sv:" << setw(8) << setprecision(3) << sv;
		real_names.push_back(ss.str());
	}
	dv_candidates.reserve(real_names, dv_names);

	for (int i=0;i<scale_vals.size();i++)
	{
		double scale_val = scale_vals[i];
		ss.str("");
		ss << "starting calcs for scaling factor" << scale_val;
		message(1, "starting lambda calcs for scaling factor", scale_val);
		message(2, "see .log file for more details");
		
		dv_candidate = fancy_solve_routine(scale_val, _current_num_dv_values);
		Eigen::VectorXd vec = dv_candidate.get_data_eigen_vec(dv_names);
		dv_candidates.update_real_ip(real_names[i], vec);
		//TODO: add sqp option to save candidates

		ss.str("");
		message(1, "finished calcs for scaling factor:", scale_val);

	}
	//// alpha trialing
	//if (search_dir != 0 or next_it == False);  // or break out of iteration here if True
	//{
	//	if (constraints);
	//	{
	//		ss.str("");
	//		ss << "specifying 'base' step length based on active set solve";
	//		string s = ss.str();
	//		message(1, s);
	//		//step = alpha
	//		throw_sqp_error("TODO");
	//	}
	//	else
	//	{
	//		ss.str("");
	//		ss << "specifying 'base' step length based on bound-related heuristics";
	//		string s = ss.str();
	//		message(1, s);
	//		//step = alpha
	//		throw_sqp_error("TODO");
	//	}

	//	// compute step length and evaluate
	//	ss.str("");
	//	ss << "compute trial step lengths and associated mean dec var vectors";
	//	string s = ss.str();
	//	message(1, s);
	//	throw_sqp_error("TODO");

	//	// some bound handling here
	//	// and drop subsequent (larger) alphas

	//	if (use_ensemble_grad)
	//	{
	//		ss.str("");
	//		ss << "draw dec var en at candidate mean vector";
	//		string s = ss.str();
	//		message(1, s);
	//		throw_sqp_error("TODO");

	//		ss.str("");
	//		ss << "evaluate ensembles at trial alphas - efficiently";
	//		string s = ss.str();
	//		message(1, s);
	//		throw_sqp_error("TODO");
	//	}
	//	else  // finite differences
	//	{
	//		ss.str("");
	//		ss << "evaluate model at trial alphas - efficiently";
	//		string s = ss.str();
	//		message(1, s);
	//		throw_sqp_error("TODO");
	//	}

	//	if (constraints); // and if active set size is > 1?
	//	{
	//		ss.str("");
	//		ss << "adopt filtering method to handle constraints and select best step size";
	//		string s = ss.str();
	//		message(1, s);
	//		//self._filter, accept, c_viol = self._filter_constraint_eval(self.obsensemble_1, self._filter, step_size,
	//			//	biobj_weight = biobj_weight, biobj_transf = biobj_transf,
	//				//opt_direction = self.opt_direction)
	//		throw_sqp_error("TODO");

	//		if (accept);
	//		{
	//			// tracking

	//			ss.str("");
	//			ss << "check for blocking constraints; break if so";
	//			string s = ss.str();
	//			message(1, s);
	//			throw_sqp_error("TODO");
	//		}
	//	}
	//	else  // unconstrained
	//	{
	//		if ((LBFGS) & (iter_num > 1));  // BFGS comes later
	//		{
	//			ss.str("");
	//			ss << "check first (sufficiency) Wolfe condition for each alpha; next alpha if not";
	//			string s = ss.str();
	//			message(1, s);
	//			throw_sqp_error("TODO");

	//			ss.str("");
	//			ss << "evaluate grad at candidate step size";
	//			string s = ss.str();
	//			message(1, s);
	//			// can also do some re-use here
	//			throw_sqp_error("TODO");

	//			ss.str("");
	//			ss << "check second (curvature) Wolfe condition for each alpha; next alpha if not";
	//			string s = ss.str();
	//			message(1, s);
	//			throw_sqp_error("TODO");
	//		}
	//		// can compute curvature (y^Ts) here for reference and phi-curvature trade off

	//		if; //current alpha is best
	//		{
	//			// tracking
	//		}
	//	}
	//}
	
	if (pest_scenario.get_pestpp_options().get_ies_debug_upgrade_only())
	{
		message(0, "ies_debug_upgrade_only is true, exiting");
		throw_sqp_error("ies_debug_upgrade_only is true, exiting");
	}

	vector<map<int, int>> real_run_ids_lams;
	int best_idx = -1;
	double best_mean = 1.0e+30, best_std = 1.0e+30;
	double mean, std;

	message(0, "running candidate decision variable ensembles");
	ObservationEnsemble oe_candidates = run_candidate_ensemble(dv_candidates, scale_vals);

	//todo: decide which if any dv candidate to accept...

	//// deal with unsuccessful iteration

	//// quasi-Newton Hessian updating via BFGS
	//// only if combination of conditions satisfies (some of which are user-specified)
	//if (Hessian update or self scale) & (BFGS);
	//{
	//	ss.str("");
	//	ss << "update Hessian via standard quasi-Newton BFGS";
	//	string s = ss.str();
	//	message(1, s);
	//	throw_sqp_error("TODO");
	//}
	//else
	//{
	//	ss.str("");
	//	ss << "skipping Hessian scaling and updating";
	//	string s = ss.str();
	//	message(1, s);
	//}
	
	//report and save, etc...
	
	return true;
}


void SeqQuadProgram::report_and_save()
{
	ofstream &frec = file_manager.rec_ofstream();
	frec << endl << "  ---  SeqQuadProgram iteration " << iter << " report  ---  " << endl;
	frec << "   number of active realizations:  " << dv.shape().first << endl;
	frec << "   number of model runs:           " << run_mgr_ptr->get_total_runs() << endl;

	cout << endl << "  ---  SeqQuadProgram iteration " << iter << " report  ---  " << endl;
	cout << "   number of active realizations:   " << dv.shape().first << endl;
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
		dv.to_binary(ss.str());
	}
	else
	{
		ss << file_manager.get_base_filename() << "." << iter << ".par.csv";
		dv.to_csv(ss.str());
	}
	save_base_real_par_rei(pest_scenario, dv, oe, output_file_writer, file_manager, iter);
	save_base_real_par_rei(pest_scenario, dv, oe, output_file_writer, file_manager, -1);
	//ss << file_manager.get_base_filename() << "." << iter << ".par.csv";
	//dv.to_csv(ss.str());
	frec << "      current par ensemble saved to " << ss.str() << endl;
	cout << "      current par ensemble saved to " << ss.str() << endl;

	

}


//void SeqQuadProgram::set_subset_idx(int size)
//{
//	//map<int,int> subset_idx_map;
//	subset_idxs.clear();
//	int nreal_subset = pest_scenario.get_pestpp_options().get_ies_subset_size();
//	if ((!use_subset) || (nreal_subset >= size))
//	{
//		for (int i = 0; i < size; i++)
//			subset_idxs.push_back(i);
//		return;
//	}
//	vector<string> pe_names = dv.get_real_names();
//
//	vector<string>::iterator bidx = find(pe_names.begin(), pe_names.end(), base_name);
//	if (bidx != pe_names.end())
//	{
//
//		subset_idxs.push_back(bidx - pe_names.begin());
//	}
//	//int size = dv.shape().first;
//	string how = pest_scenario.get_pestpp_options().get_ies_subset_how();
//	if (how == "FIRST")
//	{
//		for (int i = 0; i < size; i++)
//		{
//			if (subset_idxs.size() >= nreal_subset)
//				break;
//			if (find(subset_idxs.begin(), subset_idxs.end(), i) != subset_idxs.end())
//				continue;
//
//			subset_idxs.push_back(i);
//
//		}
//
//	}
//	else if (how == "LAST")
//	{
//
//		for (int i = size-1; i >= 0; i--)
//		{
//			if (subset_idxs.size() >= nreal_subset)
//				break;
//			if (find(subset_idxs.begin(), subset_idxs.end(), i) != subset_idxs.end())
//				continue;
//
//			subset_idxs.push_back(i);
//
//		}
//
//	}
//
//	else if (how == "RANDOM")
//	{
//		std::uniform_int_distribution<int> uni(0, size-1);
//		int idx;
//		for (int i = 0; i < 10000000; i++)
//		{
//			if (subset_idxs.size() >= nreal_subset)
//				break;
//			idx = uni(subset_rand_gen);
//			if (find(subset_idxs.begin(), subset_idxs.end(), idx) != subset_idxs.end())
//				continue;
//			subset_idxs.push_back(idx);
//		}
//		if (subset_idxs.size() != nreal_subset)
//			throw_sqp_error("max iterations exceeded when trying to find random subset idxs");
//
//	}
//	else if (how == "PHI_BASED")
//	{
//		//sidx needs to be index of realization, not realization number
//		vector<pair<double, int>> phis;
//		//vector<int> sidx;
//		int step;
//		int idx;
//		L2PhiHandler::phiType pt = L2PhiHandler::phiType::COMPOSITE;
//		map<string, double> phi_map = ph.get_phi_map(pt);
//		map<string, double>::iterator pi = phi_map.begin(), end = phi_map.end();
//
//		int i = 0;
//		for (; pi != end; ++pi)
//		{
//			phis.push_back(make_pair(pi->second, i)); //phival,idx?
//			++i;
//		}
//		sort(phis.begin(), phis.end());
//
//		//include idx for lowest and highest phi reals
//		if (subset_idxs.size() < nreal_subset)
//		{
//			for (auto phi : phis)
//			{
//				if (find(subset_idxs.begin(), subset_idxs.end(), phi.second) == subset_idxs.end())
//				{
//					subset_idxs.push_back(phi.second);
//					break;
//				}
//			}
//		}
//		if (subset_idxs.size() < nreal_subset)
//		{
//			for (int i = phis.size() - 1; i >= 0; i--)
//			{
//				if (find(subset_idxs.begin(), subset_idxs.end(), phis[i].second) == subset_idxs.end())
//				{
//					subset_idxs.push_back(phis[i].second);
//					break;
//				}
//			}
//		}
//
//
//		step = (phis.size()-1) / nreal_subset;
//		//cout << step << endl;
//		//cout << (phis.size() - 1) << endl;
//		for (i = 1; i < nreal_subset; ++i)
//		{
//			//add higher phis first
//			idx = phis.size() - (i * step);
//			if ((subset_idxs.size() < nreal_subset) && (find(subset_idxs.begin(), subset_idxs.end(), phis[idx].second) == subset_idxs.end()))
//			{
//				subset_idxs.push_back(phis[idx].second);
//				//cout << i << endl;
//				//cout << idx << endl;
//				//cout << phis[idx].first << endl;
//				//cout << phis[idx].second << endl;
//			}
//		}
//	}
//	else
//	{
//		//throw runtime_error("unkonwn 'subset_how'");
//		throw_sqp_error("unknown 'subset_how'");
//	}
//	stringstream ss;
//	for (auto i : subset_idxs)
//		ss << i << ":" << pe_names[i] << ", ";
//	message(1,"subset idx:dv real name: ",ss.str());
//	return;
//	//return subset_idx_map;
//}

ObservationEnsemble SeqQuadProgram::run_candidate_ensemble(ParameterEnsemble& dv_candidates, vector<double> &scale_vals)
{
	ofstream &frec = file_manager.rec_ofstream();
	stringstream ss;
	ss << "queuing " << dv_candidates.shape().first << " candidate solutions";
	performance_log->log_event(ss.str());
	run_mgr_ptr->reinitialize();
	
	//set_subset_idx(dv_candidates[0].shape().first);
	map<int, int> real_run_ids;
	//ParameterEnsemble pe_lam;
	//for (int i=0;i<pe_lams.size();i++)
	try
	{
		real_run_ids = dv_candidates.add_runs(run_mgr_ptr);
	}
	catch (const exception &e)
	{
		stringstream ss;
		ss << "run_ensemble() error queueing runs: " << e.what();
		throw_sqp_error(ss.str());
	}
	catch (...)
	{
		throw_sqp_error(string("run_ensembles() error queueing runs"));
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
		throw_sqp_error(ss.str());
	}
	catch (...)
	{
		throw_sqp_error(string("error running ensembles"));
	}

	performance_log->log_event("processing runs");
	vector<int> failed_real_indices;
	
	ObservationEnsemble _oe = oe;//copy
	try
	{
		failed_real_indices = _oe.update_from_runs(real_run_ids, run_mgr_ptr);
	}
	catch (const exception &e)
	{
		stringstream ss;
		ss << "error processing dv candidate runs: " << e.what();
		throw_sqp_error(ss.str());
	}
	catch (...)
	{
		stringstream ss;
		ss << "error processing dv candidate runs";
		throw_sqp_error(ss.str());
	}

	if (pest_scenario.get_pestpp_options().get_ies_debug_fail_subset())
		failed_real_indices.push_back(real_run_ids.size()-1);

	if (failed_real_indices.size() > 0)
	{
		stringstream ss;
		vector<string> par_real_names = dv.get_real_names();
		vector<string> obs_real_names = oe.get_real_names();
		vector<string> failed_par_names, failed_obs_names;
		string oname, pname;
		ss << "the following dv candidate runs failed -->";
		for (auto& i : failed_real_indices)
		{
			pname = par_real_names[i];
			oname = obs_real_names[i];
			failed_par_names.push_back(pname);
			failed_obs_names.push_back(oname);
			ss << pname << ":" << oname << ',';
		}
		string s = ss.str();
		message(1, s);
		if (failed_real_indices.size() == _oe.shape().first)
		{
			message(0, "WARNING: all dv candidate runs failed");
			_oe = ObservationEnsemble();

		}
		else
		{
			performance_log->log_event("dropping failed realizations");
			//_oe.drop_rows(failed_real_indices);
			//pe_lams[i].drop_rows(failed_real_indices);
			_oe.drop_rows(failed_obs_names);
			dv_candidates.drop_rows(failed_par_names);
		}
	}

	
	return _oe;
}

void SeqQuadProgram::queue_chance_runs()
{
	/* queue up chance-related runs using the class attributes dp and op*/
	if (pest_scenario.get_control_info().noptmax == 0)
		return;
	stringstream ss;
	if (constraints.should_update_chance(iter))
	{
		//just use dp member nearest the mean dec var values
		dv.transform_ip(ParameterEnsemble::transStatus::NUM);
		vector<double> t = dv.get_mean_stl_var_vector();
		Eigen::VectorXd dv_mean = stlvec_2_eigenvec(t);
		t.resize(0);
			
		ss << "using mean decision variables for chance calculations";
		
		Parameters pars = pest_scenario.get_ctl_parameters();
		pest_scenario.get_base_par_tran_seq().ctl2numeric_ip(pars);
		pars.update_without_clear(dv.get_var_names(), dv_mean);
		Observations obs = pest_scenario.get_ctl_observations();
		pest_scenario.get_base_par_tran_seq().numeric2ctl_ip(pars);
		constraints.add_runs(iter, pars, obs, run_mgr_ptr);
	}
}



vector<int> SeqQuadProgram::run_ensemble(ParameterEnsemble &_pe, ObservationEnsemble &_oe, const vector<int> &real_idxs)
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
		throw_sqp_error(ss.str());
	}
	catch (...)
	{
		throw_sqp_error(string("run_ensemble() error queueing runs"));
	}

	queue_chance_runs();

	performance_log->log_event("making runs");
	try
	{
		run_mgr_ptr->run();
	}
	catch (const exception &e)
	{
		stringstream ss;
		ss << "error running ensemble: " << e.what();
		throw_sqp_error(ss.str());
	}
	catch (...)
	{
		throw_sqp_error(string("error running ensemble"));
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
		throw_sqp_error(ss.str());
	}
	catch (...)
	{
		throw_sqp_error(string("error processing runs"));
	}
	

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

	constraints.process_runs(run_mgr_ptr, iter);

	return failed_real_indices;
}


void SeqQuadProgram::finalize()
{

}

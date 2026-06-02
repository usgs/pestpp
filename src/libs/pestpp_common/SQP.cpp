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
#include "EnsembleSmoother.h"


bool SqpFilter::accept(double obj_val, double violation_val, double violation_padded, Parameters p, Observations o, string rname, int iter, bool keep)
{
	FilterRec candidate{ obj_val, violation_val,iter, p, o, rname, violation_padded };
	if (obj_viol_pairs.size() == 0)
	{
		obj_viol_pairs.insert(candidate);
		return true;
	}

	FilterRec candidate_with_tol = candidate;
	if (minimize)
		candidate_with_tol.obj_val *= (1 + obj_tol);
	else
		candidate_with_tol.obj_val *= (1 - obj_tol);
	candidate_with_tol.viol_val *= (1 + viol_tol);

	bool accept = true;
	for (auto& p : obj_viol_pairs)
		if (!first_partially_dominates_second(candidate_with_tol, p))
		{
			accept = false;
			break;
		}
	if ((keep) && (accept))
	{
		obj_viol_pairs.insert(candidate);
	}
	return accept;

}


/**
 * @brief First partially dominates second.
 *
 * @param first Description.
 * @param second Description.
 *
 * @return Description.
 */
bool SqpFilter::first_partially_dominates_second(const FilterRec& first, const FilterRec& second)
{
	if (minimize)
	{
		if ((first.obj_val < second.obj_val) || (first.viol_val < second.viol_val))
			return true;
		else
			return false;
	}
	else
	{
		if ((first.obj_val > second.obj_val) || (first.viol_val < second.viol_val))
			return true;
		else
			return false;
	}
}

/**
 * @brief First strictly dominates second.
 *
 * @param first Description.
 * @param second Description.
 *
 * @return Description.
 */
bool SqpFilter::first_strictly_dominates_second(const FilterRec& first, const FilterRec& second)
{
    if (minimize)
    {
        if ((first.obj_val < second.obj_val) && (first.viol_val < second.viol_val))
            return true;
        else
            return false;
    }
    else
    {
        if ((first.obj_val > second.obj_val) && (first.viol_val < second.viol_val))
            return true;
        else
            return false;
    }
}

/**
 * @brief Report.
 *
 * @param frec Description.
 * @param iter Description.
 */
void SqpFilter::report(ofstream& frec, int iter)
{
	frec << endl << "SQP filter members (" << obj_viol_pairs.size() << ") for iteration " << iter << ":" << endl;
	frec << string(30, '-') << endl;
	frec << left << setw(15) << "Obj" << setw(15) << "Violation" << endl;
	frec << string(30, '-') << endl; 

	vector<FilterRec> sorted_pairs(obj_viol_pairs.begin(), obj_viol_pairs.end());
	sort(sorted_pairs.begin(), sorted_pairs.end(),
		[](const FilterRec& a, const FilterRec& b) {
			return a.obj_val > b.obj_val;
		});

	double omin = 1.0e+300, omax = -1e+300, vmin = 1e+300, vmax = -1e+300;
	for (auto& fr : sorted_pairs)
	{
		frec << left << setw(15) << setprecision(6) << fr.obj_val
			<< setw(15) << setprecision(6) << fr.viol_val << endl;
		omin = min(fr.obj_val, omin);
		omax = max(fr.obj_val, omax);
		vmin = min(fr.viol_val, vmin);
		vmax = max(fr.viol_val, vmax);
	}
	frec << string(30, '-') << endl;

    stringstream ss;
    ss.str("");
    ss << "...filter summary with " << obj_viol_pairs.size() << " pairs for iteration " << iter << ":" << endl;
    ss << "          obj min: " <<  setw(10) << omin << endl;
    ss << "          obj max: " << setw(10) << omax << endl;
    ss << "    violation min: " << setw(10) << vmin << endl;
    ss << "    violation max: " << setw(10) << vmax << endl;
    ss << endl;

    frec << ss.str();
    cout << ss.str();

}

vector<FilterRec> SqpFilter::get_feasible_solutions(bool padded) const 
{
	vector<FilterRec> feasible;
	if (!padded)
	{
		for (const auto& rec : obj_viol_pairs)
		{
			if (rec.viol_val <= 1E-6)
				feasible.push_back(rec);
		}
	}
	else
	{
		for (const auto& rec : obj_viol_pairs)
		{
			if (rec.viol_padded <= 1E-6)
				feasible.push_back(rec);
		}
	}
	return feasible;
}

vector<FilterRec> SqpFilter::get_filter_members() const
{
	vector<FilterRec> filterset;
	for (const auto& rec : obj_viol_pairs)
		filterset.push_back(rec);
	
	return filterset;
}

bool SqpFilter::update(double obj_val, double violation_val, double violation_padded, Parameters p, Observations o, string rname, int iter)
{
	FilterRec candidate;
	candidate.obj_val = obj_val;
	candidate.viol_val = violation_val;
	candidate.iter = iter;
	candidate.dp_val = p;
	candidate.oe_val = o;
	candidate.real_name = rname;
	candidate.viol_padded = violation_padded;
	multiset<FilterRec> updated;
	obj_viol_pairs.insert(candidate);
	bool i_is_dominated = false;
	multiset<FilterRec>::iterator first = obj_viol_pairs.begin();
    multiset<FilterRec>::iterator second = obj_viol_pairs.begin();
	for (int i=0;i<obj_viol_pairs.size();i++)
    {

	    i_is_dominated = false;
	    second = obj_viol_pairs.begin();
	    for (int j=0;j<obj_viol_pairs.size();j++)
        {
			if (i != j)
			{
				if (first_strictly_dominates_second(*second, *first)) {
					i_is_dominated = true;
					break;
				}
			}
	        second++;
        }
	    if (!i_is_dominated)
        {
			bool exists = false;
			for (const auto& u : updated)
			{
				if (fabs(u.obj_val - first->obj_val) <= 1E-10 && fabs(u.viol_val - first->viol_val) <= 1E-10)
					exists = true;
			}
			if (!exists)
				updated.insert(*first);
        }
	    first++;

    }
	obj_viol_pairs = updated;

	return true;
 }

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
	oe.set_pest_scenario_ptr(&pest_scenario);
	dv.set_rand_gen(&rand_gen);
	oe.set_rand_gen(&rand_gen);

	
}

/**
 * @brief Throw sqp error.
 *
 * @param message Description.
 */
void SeqQuadProgram::throw_sqp_error(string message)
{
	performance_log->log_event("SeqQuadProgram error: " + message);
	cout << endl << "   ************   " << endl << "    SeqQuadProgram error: " << message << endl << endl;
	file_manager.rec_ofstream() << endl << "   ************   " << endl << "    SeqQuadProgram error: " << message << endl << endl;
	file_manager.close_file("rec");
	performance_log->~PerformanceLog();
	throw runtime_error("SeqQuadProgram error: " + message);
}

bool SeqQuadProgram::initialize_dv(Covariance &cov)
{
	stringstream ss;
	int num_reals = pest_scenario.get_pestpp_options().get_sqp_num_reals();
	string dv_file = pest_scenario.get_pestpp_options().get_sqp_dv_en();
	bool drawn = false;
	if (dv_file.size() == 0)
	{
		ofstream& frec = file_manager.rec_ofstream();
		ParameterEnsemble dv_ensemble(&pest_scenario, &rand_gen);
		if (dv_names.size() > 0)
		{
			Covariance dv_cov = cov.get(dv_names);
			map<string, double> par_means = pest_scenario.get_ext_file_double_map("PARAMETER DATA EXTERNAL", MEAN_REAL_NAME);
			Parameters draw_dv_par = pest_scenario.get_ctl_parameters().get_subset(dv_names.begin(), dv_names.end());

			if (par_means.size() > 0)
			{
				frec << "Note: the following decision variables contain 'mean' value information that will be used in place of " << endl;
				frec << "      the 'parval1' values as mean values during ensemble generation" << endl;
				double lb, ub;
				for (auto par_mean : par_means)
				{
					if (draw_dv_par.find(par_mean.first) != draw_dv_par.end())
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
							draw_dv_par[par_mean.first] = par_mean.second;
							frec << par_mean.first << " " << par_mean.second << endl;
						}
					}
				}
			}

			message(1, "drawing decision variable realizations: ", num_reals);
			dv_ensemble.draw(num_reals, draw_dv_par, dv_cov, performance_log, pest_scenario.get_pestpp_options().get_ies_verbose_level(), file_manager.rec_ofstream());
		}

		ParameterEnsemble uncertain_ensemble(&pest_scenario, &rand_gen);
		double opt_risk = pest_scenario.get_pestpp_options().get_opt_risk();
		double sqp_risk = pest_scenario.get_pestpp_options().get_sqp_risk();
		bool use_chance = (opt_risk != 0.5) || (sqp_risk != 0.5);
		if (adj_par_names.size() > 0 && use_chance)
		{
			Covariance unc_cov = uncertain_parcov.get(adj_par_names);
			map<string, double> par_means = pest_scenario.get_ext_file_double_map("PARAMETER DATA EXTERNAL", MEAN_REAL_NAME);
			Parameters draw_unc_par = pest_scenario.get_ctl_parameters().get_subset(adj_par_names.begin(), adj_par_names.end());

			if (par_means.size() > 0)
			{
				double lb, ub;
				for (auto par_mean : par_means)
				{
					if (draw_unc_par.find(par_mean.first) != draw_unc_par.end())
					{
						lb = pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(par_mean.first)->lbnd;
						ub = pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(par_mean.first)->ubnd;
						if (par_mean.second < lb)
						{
							frec << "Warning: 'mean' value for uncertain parameter " << par_mean.first << " less than lower bound, using 'parval1'";
						}
						else if (par_mean.second > ub)
						{
							frec << "Warning: 'mean' value for uncertain parameter " << par_mean.first << " greater than upper bound, using 'parval1'";
						}
						else
						{
							draw_unc_par[par_mean.first] = par_mean.second;
							frec << par_mean.first << " " << par_mean.second << endl;
						}
					}
				}
			}

			message(1, "drawing uncertain parameter realizations: ", num_reals);
			uncertain_ensemble.draw(num_reals, draw_unc_par, unc_cov, performance_log,
				pest_scenario.get_pestpp_options().get_ies_verbose_level(),
				file_manager.rec_ofstream());
		}

		if (dv_names.size() > 0 && adj_par_names.size() > 0 && use_chance)
		{
			vector<string> real_names = dv_ensemble.get_real_names();
			uncertain_ensemble.reorder(real_names, vector<string>());

			vector<string> all_par_names = dv_names;
			all_par_names.insert(all_par_names.end(), adj_par_names.begin(), adj_par_names.end());

			dv.reserve(real_names, all_par_names);
			dv.set_trans_status(dv_ensemble.get_trans_status());

			for (auto& real_name : real_names)
			{
				Eigen::VectorXd dv_vec = dv_ensemble.get_real_vector(real_name);
				Eigen::VectorXd unc_vec = uncertain_ensemble.get_real_vector(real_name);

				Eigen::VectorXd combined_vec(dv_vec.size() + unc_vec.size());
				combined_vec << dv_vec, unc_vec;
				dv.update_real_ip(real_name, combined_vec);

			}
		}
		else if (dv_names.size() > 0)
			dv = dv_ensemble;
		
		else if (adj_par_names.size() > 0)
			dv = uncertain_ensemble;

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
				dv.enforce_bounds(performance_log, false);
			    //throw_sqp_error("not implemented");
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

/**
 * @brief Add current as bases.
 *
 * @param _dv Description.
 * @param _oe Description.
 */
void SeqQuadProgram::add_current_as_bases(ParameterEnsemble& _dv, ObservationEnsemble& _oe)
{
	//check that 'base' isn't already in ensemble
	vector<string> rnames = _dv.get_real_names();

	bool inpar = false;
	if (find(rnames.begin(), rnames.end(), BASE_REAL_NAME) != rnames.end())
	{
		message(1, "'base' realization already in parameter ensemble, ignoring 'include_base'");
		inpar = true;
	}
	else
	{
		message(1, "adding 'base' parameter values to ensemble");
		Parameters pars = pest_scenario.get_ctl_parameters();
		pars.update_without_clear(dv_names,current_ctl_dv_values.get_data_vec(dv_names));
		_dv.get_par_transform().active_ctl2numeric_ip(pars);
		//BASE simply added, no dropping/replacement of last row -- BASE not counted in num_reals for StoSAG
		//vector<int> drop{ _dv.shape().first - 1 };
		//_dv.drop_rows(drop);
		_dv.append(BASE_REAL_NAME, pars);
	}

	rnames = _oe.get_real_names();
	if (find(rnames.begin(), rnames.end(), BASE_REAL_NAME) != rnames.end())
	{
		message(1, "'base' realization already in observation ensemble, ignoring 'include_base'");
	}
	else
	{
		Observations obs = pest_scenario.get_ctl_observations();
		if (inpar)
		{
			vector<string> prnames = _dv.get_real_names();

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
			_oe.drop_rows(drop);
			_oe.append(BASE_REAL_NAME, obs);
			//rnames.insert(rnames.begin() + idx, string(base_name));
			rnames[idx] = BASE_REAL_NAME;
			_oe.reorder(rnames, vector<string>());
		}
		else
		{
			message(1, "adding 'base' observation values to ensemble");
			/*vector<int> drop{ _oe.shape().first - 1 };
			_oe.drop_rows(drop);*/
			_oe.append(BASE_REAL_NAME, obs);
		}
	}
}

template<typename T, typename A>
/**
 * @brief Message.
 *
 * @param level Description.
 * @param _message Description.
 * @param _extras Description.
 * @param echo Description.
 */
void SeqQuadProgram::message(int level, const string &_message, vector<T, A> _extras, bool echo)
{
	stringstream ss;
	if (level == 0)
		ss << endl << "  ---  ";
	else if (level == 1)
		ss << "...";
	else if (level == 2)
		ss << "   ";
	ss << _message;
	if (_extras.size() > 0)
	{

		for (auto &e : _extras)
			ss << e << " , ";

	}
	if (level == 0)
		ss << "  ---  ";
	if ((echo) && ((verbose_level >= 2) || (level < 3)))
		cout << ss.str() << endl;
	file_manager.rec_ofstream() <<ss.str() << endl;
	performance_log->log_event(_message);

}

/**
 * @brief Message.
 *
 * @param level Description.
 * @param _message Description.
 */
void SeqQuadProgram::message(int level, const string &_message)
{
	message(level, _message, vector<string>());
}

template<typename T>
/**
 * @brief Message.
 *
 * @param level Description.
 * @param _message Description.
 * @param extra Description.
 */
void SeqQuadProgram::message(int level, const string &_message, T extra)
{
	stringstream ss;
	ss << _message << " " << extra;
	string s = ss.str();
	message(level, s);
}

/**
 * @brief Sanity checks.
 */
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
	}
	else if (pest_scenario.get_control_info().pestmode == ControlInfo::PestMode::UNKNOWN)
	{
		warnings.push_back("unrecognized 'pestmode', using 'estimation'");
	}
	if ((use_ensemble_grad) && (ppo->get_sqp_num_reals() < warn_min_reals) && (par_csv.size() == 0))
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
}

/**
 * @brief Initialize objfunc.
 */
void SeqQuadProgram::initialize_objfunc()
{
	//initialize the objective function
	obj_func_str = pest_scenario.get_pestpp_options().get_opt_obj_func();
	obj_sense = (pest_scenario.get_pestpp_options().get_opt_direction() == 1) ? "minimize" : "maximize";

	ofstream& f_rec = file_manager.rec_ofstream();

	//check if the obj_str is an observation
	use_obj_obs = false;
	use_obj_pi = false;
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
		message(1, "using observation '" + obj_func_str + "' as the objective function");

		string obs_group = pest_scenario.get_ctl_observation_info().get_group(obj_func_str);
		pair<Constraints::ConstraintSense, string> sense = constraints.get_sense_from_group_name(obs_group);
		if (sense.first == Constraints::ConstraintSense::greater_than)
		{
			obj_sense = "maximize";
			message(1, "observation group '" + obs_group + "' indicates maximize objective");
		}
		else if (sense.first == Constraints::ConstraintSense::less_than)
		{
			obj_sense = "minimize";
			message(1, "observation group '" + obs_group + "' indicates minimize objective");
		}
	}

	else
	{
		if (obj_func_str.size() == 0)
		{
			
			message(0, " warning: no ++opt_objective_function-->forming a generic objective function (1.0 coef for each decision var)");
			ParameterInfo pi = pest_scenario.get_ctl_parameter_info();
			for (auto& name : dv_names)
			{
				if (pi.get_parameter_rec_ptr(name)->tranform_type != ParameterRec::TRAN_TYPE::NONE)
				{
					throw_sqp_error("only 'none' type decision variable transform supported for generic obj function");
				}
				obj_func_coef_map[name] = 1.0;
			}
				
		}

		//or if it is a prior info equation
		else if (pest_scenario.get_prior_info().find(obj_func_str) != pest_scenario.get_prior_info().end())
		{
			message(1, "using prior information equation '" + obj_func_str + "' as the objective function");
			obj_func_coef_map = pest_scenario.get_prior_info().get_pi_rec(obj_func_str).get_atom_factors();
			use_obj_pi = true;

			string pi_group = pest_scenario.get_prior_info().get_pi_rec(obj_func_str).get_group();
			pair<Constraints::ConstraintSense, string> sense = constraints.get_sense_from_group_name(pi_group);
			if(sense.first == Constraints::ConstraintSense::greater_than)
			{
				obj_sense = "maximize";
				message(1, "prior info equation group '" + pi_group + "' indicates maximize objective");
			}
			else if (sense.first == Constraints::ConstraintSense::less_than)
			{
				obj_sense = "minimize";
				message(1, "prior info equation group '" + pi_group + "' indicates minimize objective");
			}
				
		}

		else
		{
			
			//check if this obj_str is a filename
			ifstream if_obj(obj_func_str);
			if (!if_obj.good())
				throw_sqp_error("unrecognized ++opt_objective_function arg (tried file name, obs name, prior info name): " + obj_func_str);
			else
			{
				message(1, "loading objective function coefficients from ascii file ", obj_func_str);
				obj_func_coef_map = pest_utils::read_twocol_ascii_to_map(obj_func_str);
				ParameterInfo pi = pest_scenario.get_ctl_parameter_info();
				for (auto& name : dv_names)
				{
					if (pi.get_parameter_rec_ptr(name)->tranform_type != ParameterRec::TRAN_TYPE::NONE)
					{
						throw_sqp_error("only 'none' type decision variable transform supported for external file obj function");
					}
				}
			}
		}


		//check that all obj_coefs are decision vars
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


/**
 * @brief Initialize restart.
 *
 * @return Description.
 */
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


/**
 * @brief Initialize parcov.
 */
void SeqQuadProgram::initialize_parcov()
{
	stringstream ss;
	performance_log->log_event("initializing parcov");

	if (pest_scenario.get_pestpp_options().get_ies_use_empirical_prior())
		return;

	string cov_fname = pest_scenario.get_pestpp_options().get_parcov_filename();
	ofstream& frec = file_manager.rec_ofstream();

	if (!cov_fname.empty())
	{
		Covariance full_parcov;
		string how = full_parcov.try_from(pest_scenario, file_manager);
		message(1, "parcov loaded ", how);

		if (dv_names.size() > 0)
		{
			parcov = full_parcov.get(dv_names);
			ss.str("");
			ss << "cov for " << dv_names.size() << " decision variables";
			message(1, ss.str());
		}

		if (adj_par_names.size() > 0)
		{
			uncertain_parcov = full_parcov.get(adj_par_names);
			ss.str("");
			ss << "cov for " << adj_par_names.size() << " uncertain parameters";
			message(1, ss.str());
		}
	}
	else
	{
		map<string, double> par_std = pest_scenario.get_ext_file_double_map("parameter data external", "standard_deviation");
		const ParameterInfo& par_info = pest_scenario.get_ctl_parameter_info();

		double default_sigma_range = pest_scenario.get_pestpp_options().get_par_sigma_range();
		if (dv_names.size() > 0)
		{
			parcov.from_parameter_bounds(frec, dv_names, par_info, par_std, default_sigma_range);
			ss.str("");
			ss << "created cov for " << dv_names.size() << " decision variables using par_sigma_range = " << default_sigma_range;
			message(1, ss.str());
		}

		double uncpar_sigma_range = 4.0;
		if (adj_par_names.size() > 0)
		{
			uncertain_parcov.from_parameter_bounds(frec, adj_par_names, par_info, par_std, uncpar_sigma_range);
			ss.str("");
			ss << "created cov for " << adj_par_names.size() << " uncertain parameters using par_sigma_range = " << uncpar_sigma_range;
			message(1, ss.str());
		}
	}

}

/**
 * @brief Initialize.
 */
void SeqQuadProgram::initialize()
{	
	message(0, "initializing");
	pp_args = pest_scenario.get_pestpp_options().get_passed_args();

	iter = 1;

    //safety clear
	best_phis.clear();
	best_violations.clear();
	grad_vector_map.clear();
	constraint_mat_en.clear();
	cnames_en.clear();
	search_d_en.clear();
	lm_en.clear();
	current_obj_en.clear();
	constraint_jco_en.clear();
	hessian_en.clear();
	obj_map.clear();
	total_viol_map.clear();
	step_length_map.clear();
	ls_parent_map.clear();
	sv_lineage_map.clear();
	cname_sf_map.clear();
	selected_dv_indices.clear();
	unselected_dv_indices.clear();

	act_obs_names = pest_scenario.get_ctl_ordered_nz_obs_names();
	act_par_names = pest_scenario.get_ctl_ordered_adj_par_names();
	MAX_CONSEC_INFEAS_IES = pest_scenario.get_pestpp_options().get_sqp_max_consec_infeas_ies();
	SF_DEC_FAC = pest_scenario.get_pestpp_options().get_sqp_scale_down_factor();

	stringstream ss;
	PestppOptions* ppo = pest_scenario.get_pestpp_options_ptr();

	int subset_size = pest_scenario.get_pestpp_options().get_sqp_subset_size();
	if (subset_size >= 0)
		use_subset = true;
	else 
		use_subset = false;
	
	if (pp_args.find("PAR_SIGMA_RANGE") == pp_args.end())
	{
		message(1, "resetting par_sigma_range to 20.0");
		ppo->set_par_sigma_range(20.0);
	}
	ppo->set_opt_include_bnd_pi(false);

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

		adj_par_names.clear();
		set<string> dv_set(dv_names.begin(), dv_names.end());
		for (auto& par_name : act_par_names)
		{
			if (dv_set.find(par_name) == dv_set.end())
			{
				adj_par_names.push_back(par_name);
			}
		}
		if (adj_par_names.size() > 0)
		{
			ss.str("");
			ss << "identified " << adj_par_names.size() << " uncertain model parameters";
			message(2, ss.str());
		}
	}
	else
	{
		message(2, "using all adjustable parameters as decision variables: ", act_par_names.size());
		dv_names = act_par_names;
	}

	diagonal_scaling = Eigen::VectorXd::Ones(dv_names.size());
	constraints.initialize(dv_names, numeric_limits<double>::max());
	constraints.initial_report();
    initialize_objfunc();
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

	constraint_sense = constraints.get_constraint_sense();
	iter = 0;

	if (pest_scenario.get_control_info().noptmax == 0)
	{
		message(0, "'noptmax'=0, running control file parameter values and quitting");
		
		current_ctl_dv_values = pest_scenario.get_ctl_parameters();
		ParamTransformSeq pts = dv.get_par_transform();

		ParameterEnsemble _pe(&pest_scenario, &rand_gen);
		_pe.reserve(vector<string>(), pest_scenario.get_ctl_ordered_par_names());
		_pe.set_trans_status(ParameterEnsemble::transStatus::CTL);
		_pe.append(BASE_REAL_NAME, current_ctl_dv_values);
		string par_csv = file_manager.get_base_filename() + ".par.csv";
		//message(1, "saving parameter values to ", par_csv);
		//_pe.to_csv(par_csv);
		dv_base = _pe;
		dv_base.reorder(vector<string>(), act_par_names);
		ObservationEnsemble _oe(&pest_scenario, &rand_gen);
		_oe.reserve(vector<string>(), pest_scenario.get_ctl_ordered_obs_names());
		_oe.append(BASE_REAL_NAME, pest_scenario.get_ctl_observations());
		oe_base = _oe;
		oe_base.reorder(vector<string>(), act_obs_names);
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
		Eigen::VectorXd o = _oe.get_real_vector(BASE_REAL_NAME);
		current_obs = pest_scenario.get_ctl_observations();
		current_obs.update_without_clear(_oe.get_var_names(), o);
		save_real_par_rei(pest_scenario, _pe, _oe, output_file_writer, file_manager, -1, BASE_REAL_NAME);	
		constraints.sqp_report(0,current_ctl_dv_values, current_obs);
		return;
	}

	message(1, "using the following upgrade vector scale (e.g. 'line search') values:", ppo->get_sqp_alpha_mults());
	
	last_best = 1.0E+30;
	last_viol = 0.0;
	
	warn_min_reals = 10;
	error_min_reals = 2;
	
	message(1, "max run fail: ", ppo->get_max_run_fail());

	use_ensemble_grad = false;
	if ((ppo->get_sqp_num_reals() > 0) || (ppo->get_sqp_dv_en().size() > 0))
	{
		use_ensemble_grad = true;
		sampling_tracking_initialized = false;
	}
	sanity_checks();

	bool echo = false;
	if (verbose_level > 1)
		echo = true;

	initialize_parcov();
	if (use_cma && ppo->get_sqp_num_reals() > 0)
	{
		cma = CovMatAdap(&pest_scenario, &rand_gen, &file_manager);
		cma.initialize(dv_names.size(), ppo->get_sqp_num_reals());
		cma.set_covariance(parcov.get_matrix());
	}
	else if (use_cma && (ppo->get_sqp_num_reals() <= 0) && (ppo->get_sqp_dv_en().size() > 0))
	{
		message(1, "WARNING: CMA requires sqp_num_reals > 0, disabling CMA for finite-difference mode");
		use_cma = false;
	}
	current_ctl_dv_values = pest_scenario.get_ctl_parameters();
	current_ctl_dv_values = pest_scenario.get_ctl_parameters();
	current_obs = pest_scenario.get_ctl_observations();

	if (use_ensemble_grad)
	{
		prep_4_ensemble_grad();
	}
	else
	{
		prep_4_fd_grad();
	}

	sqp_risk = pest_scenario.get_pestpp_options().get_sqp_risk();
	if (constraints.get_use_chance())
	{
		constraints.presolve_chance_report(iter, current_obs, &oe, sqp_risk, true, "initial chance constraint report");
	}
	working_set_tol = pest_scenario.get_pestpp_options().get_sqp_working_set_tol();

	message(2, "calculating initial objective function gradient");
	current_grad_vector = calc_gradient_vector(current_ctl_dv_values);
	grad_vector_map[0] = current_grad_vector;
	
	last_best = get_obj_value(current_ctl_dv_values, current_obs);
	last_viol = constraints.get_sum_of_violations(current_ctl_dv_values, current_obs);
	if (last_viol < 1.0E-12)
		found_feasible = true;
	ss.str("");
	ss << "Initial phi, infeasibility values: " << last_best << " , " << last_viol;
	message(0, ss.str());
	constraints.sqp_report(iter, current_ctl_dv_values, current_obs);
	best_phis.push_back(last_best);
    best_violations.push_back(last_viol);

    double vpad, v = constraints.get_sum_of_violations(current_ctl_dv_values, current_obs);
	if (v < pest_scenario.get_pestpp_options().get_sqp_viol_pad())
		vpad = 0.0;
	else
		vpad = v;
	filter = SqpFilter((obj_sense == "minimize") ? true : false);
	filter.set_tol(pest_scenario.get_pestpp_options().get_sqp_filter_tol());
	filter.update(last_best, v, vpad, current_ctl_dv_values, current_obs, BASE_REAL_NAME, 0);

	if ((v > 0.0) && !use_ensemble_grad)
	{
	    message(0,"initial solution infeasible, seeking feasible solution");
		seek_feasible();
	}
	
	if (pest_scenario.get_pestpp_options().get_sqp_use_ensemble_approx_hessian())
	{
		message(2, "calculating initial objective function hessian");
		hessian = calc_objective_hessian();
		Eigen::MatrixXd hessian_dense = hessian.e_ptr()->toDense();
	}
	else
	{
		message(2, "initializing hessian matrix with identity");
		Eigen::SparseMatrix<double> h(dv_names.size(), dv_names.size());
		h.setIdentity();
		hessian = Covariance(dv_names, h);
	}
	if (pest_scenario.get_pestpp_options().get_sqp_debug_hessian())
	{
		ss.str("");
		ss << "hessian (iter 0):" << endl;
		ss << hessian.e_ptr()->toDense() << endl;
		ofstream& frec = file_manager.rec_ofstream();
		frec << ss.str() << endl;
	}

	message(0, "initialization complete");
}

/**
 * @brief Save current dv obs.
 */
void SeqQuadProgram::save_current_dv_obs()
{
    stringstream ss;
    ss.str("");
    ss << file_manager.get_base_filename() << "." << iter << "." << BASE_REAL_NAME << ".par";
    string par_name = ss.str();
    pest_utils::lower_ip(par_name);
    ofstream of(par_name);
    if (of.bad())
    {
        throw_sqp_error("error opening par file"+par_name);
    }
    const TranOffset& toff = *pest_scenario.get_base_par_tran_seq().get_offset_ptr();
    const TranScale& tsc = *pest_scenario.get_base_par_tran_seq().get_scale_ptr();
    output_file_writer.write_par(of,current_ctl_dv_values,toff,tsc);
    of.close();
    ObjectiveFunc obj_func(&(pest_scenario.get_ctl_observations()), &(pest_scenario.get_ctl_observation_info()), &(pest_scenario.get_prior_info()));
    ss.str("");
    ss << iter << "." << BASE_REAL_NAME << ".rei";
    string rei_name = ss.str();
    pest_utils::lower_ip(rei_name);
    ofstream& ofr = file_manager.open_ofile_ext(rei_name);
    output_file_writer.write_rei(ofr, iter,
                                 pest_scenario.get_ctl_observations(), current_obs, obj_func, current_ctl_dv_values);
    file_manager.close_all_files_ending_with("rei");

}

/**
 * @brief Prep 4 fd grad.
 */
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
		vnames = jco.get_sim_obs_names(); 
		snames.insert(vnames.begin(), vnames.end());
		vnames.clear();
		for (auto name : constraints.get_obs_constraint_names())
			if (snames.find(name) == snames.end())
				vnames.push_back(name);

		if (vnames.size() > 0)
		{
			for (auto m : vnames)
			ss.str("");
			ss << "existing jacobian missing the following obs constraints:" << endl;
				ss << vnames << endl;
			throw_sqp_error(ss.str());
		}
		string res_filename = pest_scenario.get_pestpp_options().get_hotstart_resfile();
		if (res_filename.size() == 0)
		{
			//make the initial base run
			cout << "  ---  running the model once with initial decision variables  ---  " << endl;
			ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();
			int run_id = run_mgr_ptr->add_run(pts.ctl2model_cp(current_ctl_dv_values));
			queue_chance_runs();

			run_mgr_ptr->run();
			bool success = run_mgr_ptr->get_run(run_id, current_ctl_dv_values, current_obs);
			if (!success)
				throw_sqp_error("initial (base) run with initial decision vars failed...cannot continue");
			pts.model2ctl_ip(current_ctl_dv_values);
			constraints.process_runs(run_mgr_ptr,iter);
		}
		else
		{
			stringstream message;
			message << "  reading  residual file " << res_filename << " for hot-start...";
			cout << message.str();
			file_manager.rec_ofstream() << message.str();
			for (auto& oname : pest_scenario.get_ctl_ordered_obs_names())
				current_obs[oname] = -1.0e+30;
			pest_utils::read_res(res_filename, current_obs);
			file_manager.rec_ofstream() << "done" << endl;
			cout << "done" << endl;
			if (constraints.get_use_chance())
			{
				queue_chance_runs();
				run_mgr_ptr->run();
				constraints.process_runs(run_mgr_ptr, iter);
			}
		}
			
	}
	else
	{
		//todo: handle hotstart_resfile here...
		bool init_obs = true;
		run_jacobian(current_ctl_dv_values, current_obs, init_obs);
	}
	jco.save("0.jcb");
	message(1, "saved initial jacobian to " + file_manager.get_base_filename() + ".0.jcb");
	save_current_dv_obs();
}

/**
 * @brief Run jacobian.
 *
 * @param _current_ctl_dv_vals Description.
 * @param _current_obs Description.
 * @param init_obs Description.
 */
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
		*run_mgr_ptr, out_of_bounds, false, init_obs,true);
	if (!success)
		throw_sqp_error("error building jacobian runs for FD grad");

	if (out_of_bounds.size() > 0)
	{
		ss.str("");
		ss << "the following decision variable are out of bounds: " << endl;
		for (auto& o : out_of_bounds)
			ss << o << ",";
		throw_sqp_error(ss.str());
	}

	queue_chance_runs();
	message(2, "starting finite difference gradient perturbation runs");
	jco.make_runs(*run_mgr_ptr);

	success = jco.process_runs(par_trans,pgi,*run_mgr_ptr,pi,false,false);
	if (!success)
	{
		throw_sqp_error("error processing finite difference gradient perturbation runs");
	}

	if (init_obs)
	{
		run_mgr_ptr->get_run(0, current_pars, _current_obs, false);
	}
}

/**
 * @brief Make gradient runs.
 *
 * @param _current_dv_vals Description.
 * @param _current_obs Description.
 */
void SeqQuadProgram::make_gradient_runs(Parameters& _current_dv_vals, Observations& _current_obs)
{
	stringstream ss;
	if (use_ensemble_grad)
	{
		
		ParameterEnsemble _dv(&pest_scenario, &rand_gen);
		ofstream& frec = file_manager.rec_ofstream();

		message(1, "generating new dv and pars ensemble at current best phi via CMA");
		if (reset)
		{
			message(1, "resetting Hessian to identity");
			Eigen::SparseMatrix<double> h(dv_names.size(), dv_names.size());
			h.setIdentity();
			hessian = Covariance(dv_names, h);

			cma.reinflate_C(1.0, true, pest_scenario.get_pestpp_options().get_sqp_max_reinflation_cond_num());
			cma.clear_archives();

			_dv = cma.generate_population(current_ctl_dv_values, dv);
			reset = false;
		}
		else if (reset_corr)
		{
			message(1, "dropping covariance elements in C");
			cma.reinflate_C(1.0, true, pest_scenario.get_pestpp_options().get_sqp_max_reinflation_cond_num());
			_dv = cma.generate_population(current_ctl_dv_values, dv);
			reset_corr = false;				
		}
		else if (seek_ies)
		{
			seek_ies = false;
			cma.reinflate_C(pest_scenario.get_pestpp_options().get_sqp_cma_reinflation_factor(), false, pest_scenario.get_pestpp_options().get_sqp_max_reinflation_cond_num());
			_dv = cma.generate_population(current_ctl_dv_values, dv);
		}
		else
			_dv = cma.generate_population(current_ctl_dv_values, dv);

		if (pest_scenario.get_pestpp_options().get_sqp_debug_cma())
		{
			ss.str("");
			ss << endl << "CMA approximated covariance: " << endl << cma.get_covariance_matrix() << endl << endl;
			frec << ss.str();

			message(0, "CMA metrics summary");
			message(2, cma.get_cma_update_summary());
		}
		
		ObservationEnsemble _oe(&pest_scenario, &rand_gen);
		_oe.reserve(_dv.get_real_names(), constraints.get_obs_constraint_names());
		add_current_as_bases(_dv, _oe);

		Eigen::VectorXd base_dv_vec;
		Eigen::VectorXd base_obs_vec;
		bool removed_base = false;
		if ((iter > 0) && (_dv.shape().first > 1))
		{
			auto real_map = _dv.get_real_map();
			auto base_itr = real_map.find(BASE_REAL_NAME);
			if (base_itr != real_map.end())
			{
				int base_row = base_itr->second;
				base_dv_vec = _dv.get_real_vector(base_row);
				base_obs_vec = _current_obs.get_data_eigen_vec(_oe.get_var_names());
				std::vector<int> drop_idx{ base_row };
				_dv.drop_rows(drop_idx);
				if (base_row < _oe.shape().first)
					_oe.drop_rows(drop_idx);
				removed_base = true;
				message(2, "skipping 'base' realization when queueing dv ensemble");
			}
		}

		message(1, "running new dv ensemble");
		run_ensemble(_dv, _oe);


		if (removed_base)
		{
			_dv.append(BASE_REAL_NAME, base_dv_vec);
			if (base_obs_vec.size() == _oe.shape().second)
				_oe.append(BASE_REAL_NAME, base_obs_vec);
			else
			{
				Observations base_obs = _current_obs;
				_oe.append(BASE_REAL_NAME, base_obs);
			}
		}
		dv = _dv;
		oe = _oe;
	}
	else
	{
		ss.str("");
		ss << iter << ".jcb";
		message(1, "running jacobian for gradient");
		run_jacobian(_current_dv_vals, _current_obs, false);
		jco.save(ss.str());
		message(1, "saved jacobian to " + file_manager.get_base_filename() + "." + ss.str());
	}
}

/**
 * @brief Prep 4 ensemble grad.
 */
void SeqQuadProgram::prep_4_ensemble_grad()
{
	stringstream ss;
	message(1, "using stochastic gradient approximation (StoSAG)");

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
	dv.transform_ip(ParameterEnsemble::transStatus::NUM);

    if (pp_args.find("SQP_RESTART_OBS_EN") != pp_args.end())
    {
        message(1, "Warning: even though `sqp_include_base` is true, you passed a restart obs en, not adding 'base' realization...");
    }
    else
        add_current_as_bases(dv, oe);

	if (!dv_drawn || !oe_drawn)
	{
		bool aligned = dv.try_align_other_rows(performance_log, oe);
		if (aligned)
		{
			message(2, "observation ensemble reordered to align rows with dv ensemble");
		}
	}

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

	if (pest_scenario.get_pestpp_options().get_save_binary())
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

	message(1, "centering on 'base' realization");

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
		_pe.append(MEAN_REAL_NAME, pars);
		string par_csv = file_manager.get_base_filename() + ".mean.par.csv";
		message(1, "saving mean dv values to ", par_csv);
		_pe.to_csv(par_csv);
		dv_base = _pe;
		dv_base.reorder(vector<string>(), act_par_names);
		ObservationEnsemble _oe(&pest_scenario, &rand_gen);
		_oe.reserve(vector<string>(), oe.get_var_names());
		_oe.append(MEAN_REAL_NAME, pest_scenario.get_ctl_observations());
		oe_base = _oe;
		oe_base.reorder(vector<string>(), act_obs_names);

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


        Eigen::VectorXd o = _oe.get_real_vector(MEAN_REAL_NAME);
        current_obs = pest_scenario.get_ctl_observations();
        current_obs.update_without_clear(_oe.get_var_names(), o);
        save_real_par_rei(pest_scenario, _pe, _oe, output_file_writer, file_manager, -1, "mean");
        constraints.sqp_report(0,current_ctl_dv_values, current_obs);

		return;
	}

	oe_org_real_names = oe.get_real_names();
	pe_org_real_names = dv.get_real_names();
	string obs_restart_csv = pest_scenario.get_pestpp_options().get_ies_obs_restart_csv();
	string par_restart_csv = pest_scenario.get_pestpp_options().get_ies_par_restart_csv();

	oe_base = oe;
	oe_base.reorder(vector<string>(), act_obs_names);

	dv_base = dv;
	dv_base.reorder(vector<string>(), act_par_names);

	if (oe_drawn)
	{
		performance_log->log_event("running initial ensemble");
		message(1, "running initial ensemble of size", oe.shape().first);
		vector<int> failed = run_ensemble(dv, oe);
		if (dv.shape().first == 0)
			throw_sqp_error("all realizations failed during initial evaluation");

		dv.transform_ip(ParameterEnsemble::transStatus::NUM);
	}

	vector<string> pinames = pest_scenario.get_ctl_ordered_pi_names();
	ObservationEnsemble combined_oe = combine_obs_and_pi(oe, dv);

	ss.str("");
	if (pest_scenario.get_pestpp_options().get_save_binary())
	{
		ss << file_manager.get_base_filename() << ".0.obs.jcb";
		combined_oe.to_binary(ss.str());
	}
	else
	{
		ss << file_manager.get_base_filename() << ".0.obs.csv";
		combined_oe.to_csv(ss.str());
	}
	message(1, "saved initial obs ensemble to", ss.str());

	save_real_par_rei(pest_scenario, dv, oe, output_file_writer, file_manager, iter);
	save_real_par_rei(pest_scenario, dv, oe, output_file_writer, file_manager, -1);

	pcs = ParChangeSummarizer(&dv_base, &file_manager, &output_file_writer);

	dv.transform_ip(ParameterEnsemble::transStatus::NUM);

	//commenting out for now -- for StoSAG, we do not need the mean of dv and oe draws...we use the current dv/oe values as is
	//perhaps add an option later? e.g.., for non StoSAG grad approx?
	/*vector<double> vals = dv.get_mean_stl_var_vector();
	vector<string> names = dv.get_var_names();
	current_ctl_dv_values.update(names, vals);
	ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();
	pts.numeric2ctl_ip(current_ctl_dv_values);
	vals = oe.get_mean_stl_var_vector();
	names = oe.get_var_names();
	current_obs.update(names, vals);*/
	current_ctl_dv_values.update(dv.get_var_names(), dv.get_real_vector(BASE_REAL_NAME));
	current_obs.update(oe.get_var_names(), oe.get_real_vector(BASE_REAL_NAME));

}



/**
 * @brief Save mat.
 *
 * @param prefix Description.
 * @param mat Description.
 */
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

bool SeqQuadProgram::try_modify_hessian() 
{
	Eigen::MatrixXd H = *hessian.e_ptr();
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(H);
	if (eigensolver.info() != Eigen::Success)
		return false;
	
	Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
	Eigen::MatrixXd eigenvectors = eigensolver.eigenvectors();

	double min_eig = eigenvalues.minCoeff();
	//jwhite: should this be a user option?
	const double min_allowed_eig = 1e-3; 

	if (min_eig >= min_allowed_eig) 
		return true; 

	//Algorithm 3.3, p. 51 Nocedal and Wright
	double tau = 2 * abs(min_eig) + min_allowed_eig;
	Eigen::MatrixXd modified_H = H + tau * Eigen::MatrixXd::Identity(H.rows(), H.cols()); 

	hessian = Covariance(dv_names, modified_H.sparseView());

	message(1, "Modified Hessian to ensure positive definiteness. tau = ", tau);
	return true;
}

bool SeqQuadProgram::hessian_update_sr1(Eigen::VectorXd s_k, Eigen::VectorXd y_k, Covariance old_hessian)
{
	message(1, "starting SR1 hessian update for iteration ", iter);

	const double eps = 1e-10;
	const double max_condition = 1e8;
	const double sr1_threshold = 1e-8; // Threshold for denominator in SR1 update

	// Check if step or gradient difference is too small
	if (s_k.norm() < eps || y_k.norm() < eps)
	{
		message(1, "skipping SR1 update - step or gradient difference too small");
		return false;
	}

	// Get current Hessian matrix
	Eigen::MatrixXd H = *old_hessian.e_ptr();

	// Calculate SR1 update components
	// Eq. 6.24 in Nocedal and Wright, p. 144
	Eigen::VectorXd Hs = H * s_k;
	Eigen::VectorXd y_minus_Hs = y_k - Hs;
	double denominator = y_minus_Hs.dot(s_k);

	// Check if SR1 update is numerically safe
	if (abs(denominator) < sr1_threshold * y_minus_Hs.norm() * s_k.norm())
	{
		message(1, "skipping SR1 update - denominator too small for numerical stability");
		return false;
	}

	// Apply SR1 update formula: H_{k+1} = H_k + (y_k - H_k*s_k)(y_k - H_k*s_k)^T / ((y_k - H_k*s_k)^T * s_k)
	// Eq. 6.24 in Nocedal and Wright, p. 144
	Eigen::MatrixXd H_new = H + (y_minus_Hs * y_minus_Hs.transpose()) / denominator;

	// Check condition number for stability
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(H_new);

	// SR1 can produce indefinite matrices, which is actually okay for SQP
	// But we should check if the condition number is reasonable
	Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
	double max_eig = eigenvalues.cwiseAbs().maxCoeff();
	double min_eig = eigenvalues.cwiseAbs().minCoeff();

	if (min_eig < eps)
	{
		message(1, "warning: very small eigenvalues detected in SR1 update");
		// We can add a small regularization to avoid numerical issues
		double reg = eps - min_eig;
		H_new += reg * Eigen::MatrixXd::Identity(H_new.rows(), H_new.cols());
	}

	double cond = max_eig / min_eig;
	if (cond > max_condition)
	{
		message(1, "warning: condition number too large in SR1 update: ", cond);
		// Apply more aggressive regularization or scaling
		H_new = 0.5 * (H_new + H); // Dampen the update
	}

	hessian = Covariance(dv_names, H_new.sparseView());
	message(2, "SR1 Hessian update complete");
	return true;

}

Covariance SeqQuadProgram::calc_objective_hessian()
{
	message(1, "starting ensemble hessian approximation for iteration ", iter);

	if (!use_ensemble_grad)
	{
		message(1, "Ensemble Hessian requires ensemble gradient mode - skipping");
		return Covariance();
	}

	
	if (dv.shape().first < 2)
	{
		message(1, "insufficient ensemble size to approximate Hessian - need at least 2 realizations");
		return Covariance();
	}

	performance_log->log_event("computing approximate Hessian from ensemble covariance");
	Eigen::MatrixXd dv_anoms = dv.get_eigen_anomalies(vector<string>(), dv_names, BASE_REAL_NAME);
	Eigen::MatrixXd dv_cov_matrix = 1.0 / (dv.shape().first - 1.0) * (dv_anoms.transpose() * dv_anoms);

	Eigen::MatrixXd obj_anoms(dv.shape().first, 1);
	if (use_obj_obs) 
	{
		obj_anoms = oe.get_eigen_anomalies(vector<string>(), vector<string>{obj_func_str}, BASE_REAL_NAME);
	}
	else
	{
		dv.update_var_map();
		map<string, int> vmap = dv.get_var_map();
		Eigen::VectorXd real;
		double oval;
		int i = 0;
		for (auto& real_name : dv.get_real_names())
		{
			oval = 0;
			real = dv.get_real_vector(real_name);
			for (auto& dv_name : dv_names)
			{
				oval += obj_func_coef_map.at(dv_name) * real(vmap.at(dv_name));
			}
			obj_anoms(i, 0) = oval;
			i++;
		}
		obj_anoms.array() -= obj_anoms.mean();
	}

	// Compute objective variance for scaling
	double obj_variance = obj_anoms.squaredNorm() / (dv.shape().first - 1.0);

	// use inverse of decvar cov scaled by obj var
	// H ~ (1/σ^2_obj) * C_dv^(-1)
	// this assumes the Hessian is proportional to the inv cov
	Eigen::MatrixXd s, V, U;
	SVD_REDSVD rsvd;
	rsvd.set_performance_log(performance_log);
	rsvd.solve_ip(dv_cov_matrix, s, U, V, pest_scenario.get_svd_info().eigthresh, pest_scenario.get_svd_info().maxsing);

	
	Eigen::MatrixXd dv_cov_inv = V * s.asDiagonal().inverse() * U.transpose();
	double scale_factor = 1.0;
	if (obj_variance > 1E-10)
	{
		scale_factor = 1.0 / obj_variance;
	}
	else
	{
		message(1, "warning: very small objective variance, using unscaled Hessian");
	}

	Eigen::MatrixXd H_stosag = scale_factor * dv_cov_inv;

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(H_stosag);
	if (eigensolver.info() != Eigen::Success)
	{
		message(1, "eigenvalue decomposition failed for StoSAG Hessian");
		return Covariance();
	}

	Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
	double min_eig = eigenvalues.minCoeff();
	const double min_allowed_eig = 1e-3;

	if (min_eig < min_allowed_eig)
	{
		double tau = 2 * abs(min_eig) + min_allowed_eig;
		H_stosag += tau * Eigen::MatrixXd::Identity(H_stosag.rows(), H_stosag.cols());
		message(1, "Modified approximate Hessian to ensure positive definiteness. tau = ", tau);
	}

	Covariance obj_hessian = Covariance(dv_names, H_stosag.sparseView());

	double cond = eigenvalues.maxCoeff() / eigenvalues.minCoeff();
	message(1, "Approximate Hessian condition number: ", cond);

	return obj_hessian;
}

bool SeqQuadProgram::hessian_update_bfgs(Eigen::VectorXd s_k, Eigen::VectorXd y_k, Covariance old_hessian)
{
	stringstream ss;
	ss.str("");

	message(2, "");
	message(1, "starting BFGS hessian update for iteration ", iter);

	const double eps = 1E-5;
	const double damping_factor = pest_scenario.get_pestpp_options().get_sqp_powell_damping_factor();
	const double max_scale = 1E6;  
	const double min_s_dot_y = 1E-8;  
	const double max_hessian_norm = 1E10;

	if (s_k.norm() < eps || y_k.norm() < eps)
	{
		ss << "skipping BFGS update - step or gradient difference too small" << endl;
		message(1, "skipping BFGS update - step or gradient difference too small");
		performance_log->log_event(ss.str());
		return false;
	}

	if (seek_ies)
	{
		bool stosag_success = false;

		if (pest_scenario.get_pestpp_options().get_sqp_use_ensemble_approx_hessian())
		{
			if (use_ensemble_grad && dv.shape().first >= 2)
			{
				Covariance obj_hessian = calc_objective_hessian();
				if (obj_hessian.get_col_names().empty())
				{
					stosag_success = true;
					ss.str("");
					ss << "successfully computed ensemble-approximated Hessian" << endl;
					message(1, ss.str());
				}
				else
				{
					ss.str("");
					ss << "Hessian approximation failed, falling back to identity" << endl;
					message(1, ss.str());
				}
			}
		}

		if (!stosag_success)
		{
			message(1, "resetting Hessian to identity");
			Eigen::SparseMatrix<double> h(dv_names.size(), dv_names.size());
			h.setIdentity();
			hessian = Covariance(dv_names, h);
		}
		ss.str("");
		ss << "BFGS Hessian update skipped (seek_ies mode)" << endl;
		message(2, "BFGS Hessian update skipped (seek_ies mode)");
		performance_log->log_event(ss.str());
		return true;
	}

	// Check curvature condition and apply Powell's damping if needed
	double s_dot_y = y_k.dot(s_k); //Eq. 6.7 Nocedal and Wright, p. 137
	if (s_dot_y <= eps * s_k.norm() * y_k.norm())
	{
		ss << "applying Powell's damping to maintain positive definiteness" << endl;
		Eigen::VectorXd Hs = (*old_hessian.e_ptr()) * s_k;
		double s_dot_Hs = s_k.dot(Hs); //skTBksk rhs of Eq 18.14 in Nocedal and Wright, p. 537

		// Powell's damping formula
		double theta = 1.0;
		if (s_dot_y < damping_factor * s_dot_Hs)
			theta = (1.0 - damping_factor) * s_dot_Hs / (s_dot_Hs - s_dot_y); //lhs of Eq. 18.15 in Nocedal and Wright, p. 537

		// Modify y_k with damping
		y_k = theta * y_k + (1.0 - theta) * Hs; //r_k before Eq. 18.15 in Nocedal and Wright, p. 537
		s_dot_y = y_k.dot(s_k);  // Recalculate with damped y_k

		if (s_dot_y < min_s_dot_y * s_k.norm() * y_k.norm())
		{
			ss.str("");
			ss << "skipping BFGS update - s^T y too small even after damping" << endl;
			message(1, ss.str());
			return false;
		}
	}

	// BFGS Update formula with scaling
	Eigen::MatrixXd H = regularize_hessian(*old_hessian.e_ptr(), "BFGS update");
	Eigen::MatrixXd H_new = H;

	// Initial scaling factor (Nocedal & Wright scaling)
	Eigen::MatrixXd H_identity = Eigen::MatrixXd::Identity(H.rows(), H.cols());
	bool is_identity = H.isApprox(H_identity, 1E-6);
	if (is_identity || iter == 1)
	{
		double scale = s_dot_y / (y_k.squaredNorm());
		
		if (scale > max_scale || scale < 1.0 / max_scale)
		{
			scale = (scale > max_scale) ? max_scale : 1.0 / max_scale;
		}

		H_new *= scale;
		ss << "applying initial scaling factor: " << scale << endl;
		performance_log->log_event(ss.str());
	}

	// First term: H_k*s_k*s_k^T*H_k from Eq. 6.19, Nocedal and Wright, p. 140
	Eigen::VectorXd Hs = H_new * s_k;
	double s_dot_Hs = s_k.dot(Hs); 
	if (abs(s_dot_Hs) < eps) {
		ss << "skipping BFGS update - s^T H s too small" << endl;
		performance_log->log_event(ss.str());
		return false;
	}

	// First term: subtract (H_k s_k s_k^T H_k) / (s_k^T H_k s_k)
	H_new -= (Hs * Hs.transpose()) / s_dot_Hs;

	// Second term: add (y_k y_k^T) / (y_k^T s_k) from Eq. 6.19, Nocedal and Wright, p. 140
	H_new += (y_k * y_k.transpose()) / s_dot_y;

	double hessian_norm = H_new.norm();

	if (hessian_norm > max_hessian_norm) 
	{
		ss.str("");
		ss << "BFGS update created very large Hessian norm: " << hessian_norm << endl;
		message(1, ss.str());
	}

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_H(H_new);
	double cond_H = eig_H.eigenvalues().maxCoeff() / eig_H.eigenvalues().minCoeff();
	const double max_condition_warning = pest_scenario.get_pestpp_options().get_sqp_hess_max_cond_num() * 1E-2;
	if (cond_H > max_condition_warning) 
	{
		ss.str("");
		ss << "WARNING: BFGS Hessian condition number large: " << cond_H << endl;
		message(2, ss.str());
	}
	
	hessian = Covariance(dv_names, H_new.sparseView());
	ss.str("");
	ss << "BFGS Hessian update complete" << endl;
	message(2, "BFGS Hessian update complete");

	return true;
}

bool SeqQuadProgram::update_hessian(string how)
{
	stringstream ss;
	ss.str("");

	if (!pest_scenario.get_pestpp_options().get_sqp_update_hessian())
	{
		message(2, "hessian_update is false...");
		return false;
	}

	if (pest_scenario.get_pestpp_options().get_sqp_hessian_update_method() == "STOSAG")
		return false;

	Covariance old_hessian = hessian;

	Eigen::VectorXd prev_grad = grad_vector_map.at(iter-1).get_data_eigen_vec(dv_names);
	Eigen::VectorXd curr_grad = grad_vector_map.at(iter).get_data_eigen_vec(dv_names);
	
	//compute gradient difference (y_k) and step (s_k): Eq. 18.13 Nocedal and Wright, p. 536; Eq. 15.25 Andrei, p. 529 
	Eigen::VectorXd y_k = curr_grad - prev_grad; 
	//Eigen::VectorXd s_k = current_ctl_dv_values.get_data_eigen_vec(dv_names) - prev_ctl_dv_values.get_data_eigen_vec(dv_names);

	//check if there's an active constraint for the current dv then compute constraint jco and update y_k
	vector<string> prev_cnames, curr_cnames;
	prev_cnames = prev_constraint_mat.get_row_names();
	current_constraint_mat = get_constraint_mat(current_ctl_dv_values, current_obs, (working_set_tol)).first;
	curr_cnames = current_constraint_mat.get_row_names();
	
	set<string> all_constraint_names;
	for (const auto& name : prev_cnames)
		all_constraint_names.insert(name);
	for (const auto& name : curr_cnames)
		all_constraint_names.insert(name);

	if (!all_constraint_names.empty()) {
		Eigen::MatrixXd curr_full_jco = Eigen::MatrixXd::Zero(all_constraint_names.size(), dv_names.size());
		Eigen::MatrixXd prev_full_jco = Eigen::MatrixXd::Zero(all_constraint_names.size(), dv_names.size());
		Eigen::VectorXd curr_full_lambda = Eigen::VectorXd::Zero(all_constraint_names.size());
		Eigen::VectorXd prev_full_lambda = Eigen::VectorXd::Zero(all_constraint_names.size());

		map<string, int> constraint_to_row;
		int row_idx = 0;
		for (const auto& name : all_constraint_names)
			constraint_to_row[name] = row_idx++;

		if (!curr_cnames.empty()) 
		{
			Eigen::MatrixXd current_jco = current_constraint_mat.e_ptr()->toDense();

			if (current_jco.cols() != dv_names.size())
			{
				vector<string> jco_var_names = current_constraint_mat.get_col_names();
				vector<int> dv_idxs;
				for (int j = 0; j < dv_names.size(); j++)
				{
					auto it = find(jco_var_names.begin(), jco_var_names.end(), dv_names[j]);
					if (it != jco_var_names.end())
					{
						dv_idxs.push_back(distance(jco_var_names.begin(), it));
					}
					else
					{
						throw_sqp_error("DV '" + dv_names[j] + "' not found in constraint jacobian");
					}
				}

				Eigen::MatrixXd dv_jco(current_jco.rows(), dv_names.size());
				for (int i = 0; i < current_jco.rows(); i++)
				{
					for (int j = 0; j < dv_idxs.size(); j++)
					{
						dv_jco(i, j) = current_jco(i, dv_idxs[j]);
					}
				}
				current_jco = dv_jco;
			}

			for (int i = 0; i < curr_cnames.size(); i++) 
			{
				int row = constraint_to_row[curr_cnames[i]];
				curr_full_jco.row(row) = current_jco.row(i);

				if (constraint_sense[curr_cnames[i]] == "less_than") 
				{
					curr_full_jco.row(row) *= -1;
				}
			}

			//Eq. 18.21 from Nocedal and Wright, p. 539, approx new lambda w/o computing new Hessian
			Eigen::BDCSVD<Eigen::MatrixXd> svd_AAT(curr_full_jco * curr_full_jco.transpose(), Eigen::ComputeThinU | Eigen::ComputeThinV);
			curr_full_lambda = svd_AAT.solve(curr_full_jco * curr_grad);
		}

		if (!prev_cnames.empty()) 
		{
			Eigen::MatrixXd previous_jco = prev_constraint_mat.e_ptr()->toDense();
			if (cnames_en.find(selected_ls_parent) == cnames_en.end()) 
			{
				ss.str("");
				ss << "update_hessian error: selected_ls_parent '" << selected_ls_parent << "' not found in cnames_en";
				throw_sqp_error(ss.str());
			}

			if (lm_en.find(selected_ls_parent) == lm_en.end())
			{
				ss.str("");
				ss << "update_hessian error: selected_ls_parent '" << selected_ls_parent << "' not found in lm_en";
				throw_sqp_error(ss.str());
			}
			
			for (int i = 0; i < prev_cnames.size(); i++) 
			{
				int row = constraint_to_row[prev_cnames[i]];
				prev_full_jco.row(row) = previous_jco.row(i);

				auto lambda_it = find(cnames_en[selected_ls_parent].begin(), cnames_en[selected_ls_parent].end(), prev_cnames[i]);
				if (lambda_it != cnames_en[selected_ls_parent].end()) 
				{
					int lambda_idx = distance(cnames_en[selected_ls_parent].begin(), lambda_it);
					if (lambda_idx < lm_en[selected_ls_parent].size()) 
						prev_full_lambda(row) = lm_en[selected_ls_parent](lambda_idx);
				}

				if (constraint_sense[prev_cnames[i]] == "less_than") 
					prev_full_jco.row(row) *= -1;
			}
		}

		y_k -= curr_full_jco.transpose() * curr_full_lambda - prev_full_jco.transpose() * prev_full_lambda;

		ss << "applied full Lagrangian correction with separate lambdas" << endl;
		performance_log->log_event(ss.str());

	}
	else 
	{
		ss << "no active constraints in either iteration, using objective gradient difference only" << endl;
		performance_log->log_event(ss.str());
	}

	if (step_k.size() != dv_names.size()) 
	{
		ss.str("");
		ss << "step_k size mismatch: expected " << dv_names.size() << ", got " << step_k.size() << ", skipping BFGS update";
		message(1, ss.str());
		return false;
	}

	if (how == "BFGS")
		return hessian_update_bfgs(step_k, y_k, old_hessian);
	else if (how == "SR1")
		return hessian_update_sr1(step_k, y_k, old_hessian);
	else
		throw_sqp_error("unknown hessian update method: " + how);
	return false;
}

void SeqQuadProgram::update_scaling(const Eigen::VectorXd& step, const Eigen::VectorXd& grad) {
	if (iter % scaling_update_frequency != 0) return;
	Eigen::VectorXd rel_step = (step.array().abs() / (grad.array().abs() + 1e-10)).matrix();

	double avg_step = rel_step.mean();
	if (avg_step > 1e-10) {
		rel_step /= avg_step;
	}

	for (int i = 0; i < dv_names.size(); i++) {
		if (rel_step(i) < 0.5) {
			diagonal_scaling(i) *= (1.0 - adaptation_rate);
		}
		else if (rel_step(i) > 2.0) {
			diagonal_scaling(i) *= (1.0 + adaptation_rate);
		}

		diagonal_scaling(i) = max(0.1, min(diagonal_scaling(i), 1000.0));
	}

	message(2, "Updated adaptive Hessian scaling");
}

Eigen::MatrixXd SeqQuadProgram::regularize_hessian(const Eigen::MatrixXd& H, const string& context)
{
	stringstream ss;

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_H(H);
	double min_eig_H = eig_H.eigenvalues().minCoeff();
	const double min_allowed_eig = 1E-6;

	const double max_condition_critical = pest_scenario.get_pestpp_options().get_sqp_hess_max_cond_num();
	const double max_condition_action = max_condition_critical * 1E-1;
	const double max_condition_warning = max_condition_action * 1E-1;

	Eigen::MatrixXd H_reg = H;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_H_reg = eig_H;

	bool report_regul = true;

	if (min_eig_H < 0) {
		double tau = abs(min_eig_H) + min_allowed_eig;
		H_reg += tau * Eigen::MatrixXd::Identity(H.rows(), H.cols());
		ss.str("");
		ss << "Hessian has negative eigenvalue (" << min_eig_H << "), applying shift tau = " << tau;
		if (!context.empty())
			ss << " (" << context << ")";
		message(1, ss.str());

		eig_H_reg.compute(H_reg);
		min_eig_H = eig_H_reg.eigenvalues().minCoeff();
	}


	double max_eig = eig_H_reg.eigenvalues().maxCoeff();
	double cond_H = max_eig / min_eig_H;

	if (min_eig_H < min_allowed_eig) 
	{
		double delta = min_allowed_eig - min_eig_H;
		H_reg += delta * Eigen::MatrixXd::Identity(H.rows(), H.cols());
		ss.str("");
		ss << "Hessian minimum eigenvalue too small after shift, adding delta = " << delta;
		if (!context.empty())
			ss << " (" << context << ")";
		message(2, ss.str());

		eig_H_reg.compute(H_reg);
		min_eig_H = eig_H_reg.eigenvalues().minCoeff();
		max_eig = eig_H_reg.eigenvalues().maxCoeff();
		cond_H = max_eig / min_eig_H;
	}


	if (cond_H > max_condition_critical) 
	{
		ss.str("");
		ss << "Hessian condition number extremely large";
		if (!context.empty())
			ss << " (" << context << ")";
		ss << ": " << cond_H;
		message(2, ss.str());

		bool use_stosag = false;
		Eigen::MatrixXd H_stosag;

		if (use_ensemble_grad && dv.shape().first >= 2 && (pest_scenario.get_pestpp_options().get_sqp_hessian_update_method() == "STOSAG"))
		{
			
			Covariance obj_hessian = calc_objective_hessian();
			if (!obj_hessian.get_col_names().empty())
			{
				H_stosag = obj_hessian.e_ptr()->toDense();
				Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_stosag(H_stosag);
				if (eig_stosag.info() == Eigen::Success)
				{
					Eigen::VectorXd eigvals_stosag = eig_stosag.eigenvalues();
					double cond_stosag = eigvals_stosag.maxCoeff() / eigvals_stosag.minCoeff();
					double min_eig_stosag = eigvals_stosag.minCoeff();

					if (cond_stosag < max_condition_critical && min_eig_stosag >= min_allowed_eig)
					{
						H_reg = H_stosag;
						use_stosag = true;
						ss.str("");
						ss << "using StoSAG-approximated Hessian, condition number: " << cond_stosag;
						if (!context.empty())
							ss << " (" << context << ")";
						message(1, ss.str());
					}
					else
					{
						ss.str("");
						ss << "StoSAG Hessian also ill-conditioned, condition number: " << cond_stosag;
						if (!context.empty())
							ss << " (" << context << ")";
						message(2, ss.str());
					}
				}
				else
				{
					ss.str("");
					ss << "StoSAG Hessian eigenvalue decomposition failed";
					if (!context.empty())
						ss << " (" << context << ")";
					message(2, ss.str());
				}
			}
			else
			{
				ss.str("");
				ss << "StoSAG Hessian computation failed, using identity";
				if (!context.empty())
					ss << " (" << context << ")";
				message(2, ss.str());
			}
		}

		if (!use_stosag)
		{
			H_reg = Eigen::MatrixXd::Identity(H.rows(), H.cols());
			ss.str("");
			ss << "using identity Hessian";
			if (!context.empty())
				ss << " (" << context << ")";
			message(2, ss.str());
		}

	}
	else if (cond_H > max_condition_action) 
	{
		ss.str("");
		ss << "Hessian condition number very large";

		if (!context.empty()) 
			ss << " (" << context << ")";
		ss << ": " << cond_H << ", applying strong regularization";
		message(2, ss.str());
		double reg_factor = 1.0 / sqrt(cond_H);
		H_reg += reg_factor * Eigen::MatrixXd::Identity(H.rows(), H.cols());
	}
	else if (cond_H > max_condition_warning) 
	{
		if (min_eig_H < min_allowed_eig) {
			double delta = min_allowed_eig - min_eig_H + 1e-6;
			H_reg += delta * Eigen::MatrixXd::Identity(H.rows(), H.cols());
			ss.str("");
			ss << "applied regularization to Hessian";
			if (!context.empty()) 
				ss << " (" << context << ")";
			ss << ", delta = " << delta;
			message(2, ss.str());
		}
	}
	else if (min_eig_H < min_allowed_eig) {
		double delta = min_allowed_eig - min_eig_H + 1e-6;
		H_reg += delta * Eigen::MatrixXd::Identity(H.rows(), H.cols());
		ss.str("");
		ss << "applied regularization to Hessian";
		if (!context.empty()) 
			ss << " (" << context << ")";
		ss << ", delta = " << delta;
		message(2, ss.str());
	}
	else
		report_regul = false;

	if (report_regul && pest_scenario.get_pestpp_options().get_sqp_debug_hessian())
	{
		ss.str("");
		ss << "   approx hessian after regularization " << "(" << context << "): " << endl;
		ss << "   " << H_reg << endl;
		ofstream& frec = file_manager.rec_ofstream();
		frec << ss.str() << endl;
	}

	used_hessian = Covariance(dv_names, H_reg.sparseView());
	return H_reg;
}

bool SeqQuadProgram::isfullrank(const Eigen::MatrixXd& mat)
{
	if (mat.rows() == 0 || mat.cols() == 0) 
		return true;

	Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(mat.transpose());
	int rank = qr.rank();
	return (rank == min(mat.rows(), mat.cols()));
}

/**
 * @brief Iterate 2 solution.
 */
void SeqQuadProgram::iterate_2_solution()
{
	stringstream ss;
	ofstream &frec = file_manager.rec_ofstream();

	bool accept = false;
	n_consec_infeas = 0;
	for (int i = 0; i < pest_scenario.get_control_info().noptmax; i++)
	{ 		
		iter++;
		if (use_ensemble_grad)
			accept = solve_new_ensemble();
		else {
			throw_sqp_error("only ensemble gradient currently implemented");
			//accept = solve_new();
		}
		if (accept)
		{
			n_consec_infeas = 0;
			if (!use_ensemble_grad)
				working_set_tol = max(0.05, working_set_tol * 0.5);
		}
		else
		{
			n_consec_infeas++;
			if (n_consec_infeas > MAX_CONSEC_INFEAS_IES)
			{
				seek_ies = true;
				ss.str("");
				ss << "number of consecutive infeasible iterations > " << MAX_CONSEC_INFEAS_IES << ", switching to IES to seek feasibility";
				message(0, ss.str());
				seek_feasible();
				n_consec_infeas = 0;
			}
		}
        
		if (use_cma && !seek_ies)
		{
			message(1, "updating CMA with approximate gradient");
			cma.update(prev_ctl_dv_values, current_ctl_dv_values, iter);

			if (pest_scenario.get_pestpp_options().get_sqp_save_cov_every() > 0)
			{
				if (iter % pest_scenario.get_pestpp_options().get_sqp_save_cov_every() == 0)
				{
					ss.str("");
					ss << file_manager.get_base_filename() << "." << iter << ".CMA";
					Covariance cma_cov(dv_names, cma.get_covariance_matrix().sparseView());
					if (pest_scenario.get_pestpp_options().get_save_binary()) {
						ss << ".jcb";
						cma_cov.to_binary_new(ss.str());
					}
					else {
						ss << ".cov";
						cma_cov.to_ascii(ss.str());
					}

				}
				string fname = ss.str();
				ss.str("");
				ss << "CMA covariance matrix for iteration " << iter << " saved to: " << fname;
				message(1, ss.str());
			}
		}

		make_gradient_runs(current_ctl_dv_values, current_obs);

		if (use_ensemble_grad)
			report_and_save_ensemble(dv, oe);
		else
			save_current_dv_obs();

		grad_vector_map[iter] = calc_gradient_vector(current_ctl_dv_values);
		current_grad_vector = grad_vector_map.at(iter);

		
		if (pest_scenario.get_pestpp_options().get_sqp_hessian_update_method() == "STOSAG")
		{
			message(2, "calculating objective function hessian");
			hessian = calc_objective_hessian();
		}

        constraints.sqp_report(iter, current_ctl_dv_values, current_obs, true);

        if (use_ensemble_grad)
        {
            ss.str("");
            ss << file_manager.get_base_filename() << "." << iter << ".pcs.csv";
            pcs.summarize(dv, ss.str());
        }

		if (should_terminate() && last_viol < 1E-10)
			break;

		if (pest_scenario.get_pestpp_options().get_sqp_update_hessian())
		{
			if (pest_scenario.get_pestpp_options().get_sqp_reset_hessian_every() > 0)
			{
				if ((iter + 1) % pest_scenario.get_pestpp_options().get_sqp_reset_hessian_every() == 0)
				{
					message(1, "resetting hessian to identity");
					Eigen::SparseMatrix<double> h(dv_names.size(), dv_names.size());
					h.setIdentity();
					hessian = Covariance(dv_names, h);
				}
			}
			else if (iter > 0 && is_good_search)
				update_hessian(pest_scenario.get_pestpp_options().get_sqp_hessian_update_method());
			else
				message(1, "skipping hessian update");

			Eigen::MatrixXd hessian_dense = hessian.e_ptr()->toDense();
			
			if (pest_scenario.get_pestpp_options().get_sqp_debug_hessian())
			{
				ss.str("");
				ss << "hessian (iter " << iter << "):" << endl;
				ss << hessian_dense << endl;
				ofstream& frec = file_manager.rec_ofstream();
				frec << ss.str() << endl;
			}
		}
	}
}

/**
 * @brief Should terminate.
 *
 * @return Description.
 */
bool SeqQuadProgram::should_terminate()
{
    stringstream ss;
	
    double phiredstp = pest_scenario.get_control_info().phiredstp;
    int nphistp = pest_scenario.get_control_info().nphistp;
    int nphinored = pest_scenario.get_control_info().nphinored;
    bool phiredstp_sat = false, nphinored_sat = false, consec_sat = false, stationary = false;
    double phi, ratio, viol;
    int count = 0;
    int nphired = 0;
    best_phi_yet = 1.0e+300;
    int best_idx_yet = -1;
    for (int i=0;i<best_phis.size();i++)
    {
		if ((best_phis[i] < best_phi_yet) && (best_violations[i] == 0))
		{
			best_phi_yet = best_phis[i];
			best_violation_yet = best_violations[i];
		}
		
    }

	for (int i = 0; i < best_phis.size(); i++)
	{
		if ((best_phis[i] <= best_phi_yet) && (best_violations[i] == 0))
			nphired++;
		else
			nphired = 0;

	}

    ss.str("");
    ss << "best phi, total violation progress" << endl;
	ss << "   " << left << setw(5) << "iter"
		          << right << setw(12) << "phi"
		          << setw(15) << "violation" << endl;

	ss << "   " << left << setw(5) << "iter"
		<< right << setw(12) << "phi"
		<< setw(15) << "violation" << endl;

	//save phi and violation data to CSV file
	ss.str("");
	ss << file_manager.get_base_filename() << ".phi_viol.summary.csv";
	string csv_filename = ss.str();
	ofstream csv_file(csv_filename);
	if (!csv_file.good())
	{
		stringstream err_ss;
		err_ss << "error opening csv file for writing: " << csv_filename;
		throw_sqp_error(err_ss.str());
	}

	csv_file << "iter,phi,violation" << endl;
	for (int i = 0; i < best_phis.size(); i++)
	{
		csv_file << i << "," << fixed << setprecision(10) << best_phis[i]
			<< "," << fixed << setprecision(10) << best_violations[i] << endl;
	}
	csv_file.close();
	message(1, "saved phi and violation summary to: ", csv_filename);

	ss.str("");
	ss << " --- phi and violation sequence --- " << endl;
    for (int i=0;i<best_phis.size();i++)
    {
        phi = best_phis[i];
        viol = best_violations[i];

		if (viol == 0.0)
		{
			ratio = (phi - best_phi_yet) / phi;
			if (ratio <= phiredstp)
				count++;
		}
		ss << "   " << left << setw(5) << i
			          << right << setw(12) << fixed << setprecision(3) << phi
			          << setw(15) << fixed << setprecision(6) << viol << endl;
    }
     cout << ss.str() << endl;
	file_manager.rec_ofstream() << ss.str() << endl;

    message(0, "phi-based termination criteria check");
    message(2, "phiredstp: ", phiredstp);
    message(2, "nphistp: ", nphistp);
    message(2, "nphinored: ", nphinored);
    message(2, "best phi yet: ", best_phi_yet);
    message(2, "number of consecutive infeasible solutions: ",n_consec_infeas);

    message(2, "number of iterations satisfying phiredstp criteria: ", count);
    if (count >= nphistp)
    {
        message(1, "number iterations satisfying phiredstp criteria > nphistp");
        phiredstp_sat = true;
    }

    message(2, "number of iterations since best yet mean phi: ", nphired);
    if (nphired >= nphinored)
    {
        message(1, "number of iterations since best yet mean phi > nphinored");
        nphinored_sat = true;
    }
    if (best_phis[best_phis.size() - 1] < numeric_limits<double>::denorm_min())
    {
        message(1, "phi is zero, all done");
        return true;
    }

	const double grad_tol = 1e-6;
	if (current_grad_vector.get_data_eigen_vec(dv_names).norm() < grad_tol)
		stationary = true;

    if ((nphinored_sat) || (phiredstp_sat) || (consec_sat))
    {
        message(1, "phi-based termination criteria satisfied, all done");
        return true;
    }
	else if (stationary)
	{
		message(1, "stationary point detected, all done");
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


/**
 * @brief Calc gradient vector from coeffs.
 *
 * @param _current_dv_values Description.
 *
 * @return Description.
 */
Eigen::VectorXd SeqQuadProgram::calc_gradient_vector_from_coeffs(const Parameters& _current_dv_values)
{
	Eigen::VectorXd grad(dv_names.size());
	//first calc the current obj function value
	double current_obj_val = 0.0;
	for (auto& dv : dv_names)
	{
		current_obj_val += obj_func_coef_map.at(dv) * _current_dv_values.get_rec(dv);
	}
	//now perturb each dec var and re calc
	//just use a plain ole perturb here since we dont
	//case
	double pert = 1.1;
	double pert_val;
	double pert_obj_val, derv, dv_val;
	int i = 0;
	for (auto& dv : dv_names)
	{
		dv_val = _current_dv_values.get_rec(dv);
		Parameters pert_dv_values = _current_dv_values;
		if (dv_val != 0.0)
			pert_val = dv_val * pert;
		else
			pert_val = dv_val + pert;
		pert_dv_values.update_rec(dv,pert_val);
		pert_obj_val = 0.0;
		for (auto& dv_name : dv_names) 
		{
			pert_obj_val += obj_func_coef_map.at(dv_name) * pert_dv_values.get_rec(dv_name);
		}
		derv = (current_obj_val - pert_obj_val) / (dv_val - pert_val);
		grad[i] = derv;
		i++;
	}
	return grad;
}


/**
 * @brief Calc gradient vector.
 *
 * @param _current_dv_values Description.
 * @param _center_on Description.
 *
 * @return Description.
 */
Parameters SeqQuadProgram::calc_gradient_vector(const Parameters& _current_dv_values, string _center_on)
{
	stringstream ss;
	Eigen::VectorXd grad(dv_names.size());
	ofstream& frec = file_manager.rec_ofstream();

	string center_on = pest_scenario.get_pestpp_options().get_ies_center_on();
	if (!_center_on.empty())
		center_on = _center_on;

	if (use_ensemble_grad)
	{
		// compute sample dec var cov matrix and its pseudo inverse
		// see eq (8) of Dehdari and Oliver 2012 SPE and Fonseca et al 2015 SPE
		performance_log->log_event("form dv cov matrix");
		Eigen::MatrixXd dv_anoms = dv.get_eigen_anomalies(vector<string>(), dv_names, BASE_REAL_NAME);
		Eigen::MatrixXd dv_cov_matrix = 1.0 / (dv.shape().first - 1.0) * (dv_anoms.transpose() * dv_anoms);

		performance_log->log_event("svd of dv cov matrix");
		Eigen::MatrixXd s, V, U, st;
		SVD_REDSVD rsvd;
		rsvd.set_performance_log(performance_log);
		rsvd.solve_ip(dv_cov_matrix, s, U, V, pest_scenario.get_svd_info().eigthresh, pest_scenario.get_svd_info().maxsing);
		Eigen::MatrixXd dv_cov_pseudoinv = V * s.asDiagonal().inverse() * U.transpose();

		// Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(dv_cov_matrix);
		// Eigen::MatrixXd dv_cov_pseudoinv = cod.pseudoInverse();



		//objective function anomalies
		//Note: These objective values already reflect the effect of uncertain parameters
		//because each real includes both dec vars and uncertain params
		performance_log->log_event("filling obj fxn anomalies");

		Eigen::MatrixXd obj_anoms(dv.shape().first, 1);
		if (use_obj_obs) {
			obj_anoms = oe.get_eigen_anomalies(vector<string>(), vector<string>{obj_func_str},BASE_REAL_NAME);
		}
		else
		{
			dv.update_var_map();
			map<string,int> vmap = dv.get_var_map();
			Eigen::VectorXd real;
			double oval;
			int i =0;
			for (auto& real_name: dv.get_real_names())
			{
				oval = 0;
				real = dv.get_real_vector(real_name);
				for (auto& dv_name : dv_names)
				{
					oval += obj_func_coef_map.at(dv_name) * real(vmap.at(dv_name));
				}
				obj_anoms(i, 0) = oval;
				i++;
			}
			obj_anoms.array() -= obj_anoms.mean();
		}

		//cross-covariance between decvar and obj
		//this cross-cov naturally accounts for the effect of uncertain params
		//bc the obj values in obj_anoms already reflect their variability
		// see eq (9) of Dehdari and Oliver 2012 SPE and Fonseca et al 2015 SPE
		Eigen::VectorXd cross_cov_vector = 1.0 / (dv.shape().first - 1.0) * (dv_anoms.transpose() * obj_anoms);

		// now compute grad vector: Chen et al. (2009) and Fonseca et al. (2015) 
		grad = dv_cov_pseudoinv * cross_cov_vector;
	}
	else
	{
		if (use_obj_obs)
		{
			vector<string> obj_name_vec{ obj_func_str };
			Eigen::MatrixXd t = jco.get_matrix(obj_name_vec, dv_names);
			grad = t.row(0);
		}
		else
		{
			grad = calc_gradient_vector_from_coeffs(_current_dv_values);
		}
	}
	Parameters pgrad = _current_dv_values;
	pgrad.update_without_clear(dv_names, grad);

	if (pest_scenario.get_pestpp_options().get_sqp_debug_stosag_grad())
	{
		ss.str("");
		ss << endl << "StoSAG-calculated gradient: " << endl << grad << endl;
		frec << ss.str();
	}

	return pgrad;
}

pair<Eigen::VectorXd, Eigen::VectorXd> SeqQuadProgram::_kkt_null_space(const Eigen::MatrixXd& G, Eigen::MatrixXd& _constraint_jco, Eigen::VectorXd& constraint_diff, Eigen::VectorXd& curved_grad, vector<string>* _cnames)
{
	stringstream ss;
	ofstream& frec = file_manager.rec_ofstream();

	ss.str("");
	ss << "starting KKT null space solve";
	performance_log->log_event(ss.str());

	Eigen::MatrixXd scaled_constraint_jco = _constraint_jco;
	Eigen::VectorXd scaled_constraint_diff = constraint_diff;

	Eigen::VectorXd row_scales(scaled_constraint_jco.rows());
	const double min_scale = 1E-12;  

	for (int i = 0; i < scaled_constraint_jco.rows(); i++)
	{
		double row_norm = scaled_constraint_jco.row(i).norm();
		if (row_norm > min_scale)
			row_scales[i] = 1.0 / row_norm;
		else
			row_scales[i] = 1.0;  
	}

	for (int i = 0; i < scaled_constraint_jco.rows(); i++)
	{
		scaled_constraint_jco.row(i) *= row_scales[i];
		scaled_constraint_diff[i] *= row_scales[i];
	}

	ss.str("");
	ss << "   applied constraint Jacobian scaling:";
	if (scaled_constraint_jco.rows() > 0)
	{
		ss << " min_scale=" << row_scales.minCoeff();
		ss << " max_scale=" << row_scales.maxCoeff();
		ss << " condition_improvement=" << (row_scales.maxCoeff() / row_scales.minCoeff());
	}
	ss << endl;
	frec << ss.str();
	performance_log->log_event(ss.str());


	Eigen::VectorXd search_d, lm;

	int m = scaled_constraint_jco.rows();
	int n = scaled_constraint_jco.cols();

	if (m == 0)
	{
		Eigen::MatrixXd G_reg = regularize_hessian(G, "null space unconstrained");
		Eigen::LDLT<Eigen::MatrixXd> ldlt_H(G_reg);
		if (ldlt_H.info() != Eigen::Success || !ldlt_H.isPositive())
			search_d = -curved_grad;
		else
			search_d = ldlt_H.solve(-curved_grad);
		return pair<Eigen::VectorXd, Eigen::VectorXd>(search_d, Eigen::VectorXd::Zero(0));
	}

	Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(scaled_constraint_jco.transpose());

	int rank = qr.rank();

	if (rank < m)
	{
		ss.str("");
		ss << "   WARNING: Constraint matrix rank-deficient after scaling. Rank = " << rank
			<< ", Rows = " << m << ". Using QR to identify independent constraints." << endl;
		frec << ss.str();
		performance_log->log_event(ss.str());

		Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm = qr.colsPermutation();
		Eigen::VectorXi perm_indices = perm.indices();

		if (_cnames != nullptr)
		{
			vector<string> reduced_cnames;
			Eigen::MatrixXd reduced_constr_jco(rank, scaled_constraint_jco.cols());
			Eigen::VectorXd reduced_constraint_diff(rank);
			Eigen::VectorXd reduced_row_scales(rank);

			for (int i = 0; i < rank; i++)
			{
				int orig_idx = perm_indices(i);
				reduced_cnames.push_back((*_cnames)[orig_idx]);
				reduced_constr_jco.row(i) = scaled_constraint_jco.row(orig_idx);
				reduced_constraint_diff(i) = scaled_constraint_diff(orig_idx);
				reduced_row_scales(i) = row_scales[orig_idx];
			}

			scaled_constraint_jco = reduced_constr_jco;
			scaled_constraint_diff = reduced_constraint_diff;
			row_scales = reduced_row_scales;
			*_cnames = reduced_cnames;
			m = rank;

			// Recompute QR with reduced matrix
			qr.compute(scaled_constraint_jco.transpose());
			rank = qr.rank();
		}
		else
		{
			ss.str("");
			ss << "Constraint matrix rank-deficient after scaling. Rank = " << rank << ", Rows = " << m;
			throw_sqp_error(ss.str());
		}
	}

	// Get Q matrix from QR: A^T = Q * R
	Eigen::MatrixXd Q = qr.householderQ();

	// Y: first 'rank' columns of Q span range(A^T) = column space of A^T
	Eigen::MatrixXd Y = Q.leftCols(rank);

	// Z: remaining columns of Q span null(A) 
	Eigen::MatrixXd Z;
	if (Q.cols() > rank)
		Z = Q.rightCols(Q.cols() - rank);
	else
		Z = Eigen::MatrixXd::Zero(n, 0);

	Eigen::MatrixXd G_reg = regularize_hessian(G, "null space solve");
	const double max_condition_warning = pest_scenario.get_pestpp_options().get_sqp_hess_max_cond_num() * 1E-2;
	const double min_allowed_eig = 1E-6;

	// solve p_range_space
	Eigen::VectorXd p_y, rhs;
	Eigen::MatrixXd AY = scaled_constraint_jco * Y; // A*Y is m×m, should be full rank

	// Check condition number of A*Y
	Eigen::JacobiSVD<Eigen::MatrixXd> svd_AY(AY);
	double cond_AY = svd_AY.singularValues()(0) / svd_AY.singularValues()(svd_AY.singularValues().size() - 1);
	if (cond_AY > max_condition_warning)
	{
		ss.str("");
		ss << "   WARNING: A*Y condition number very large: " << cond_AY << endl;
		frec << ss.str();
		performance_log->log_event(ss.str());
	}

	Eigen::LDLT<Eigen::MatrixXd> ldlt_AY(AY);
	if (ldlt_AY.info() != Eigen::Success)
	{
		ss.str("");
		ss << "   WARNING: LDLT failed for A*Y, falling back to SVD" << endl;
		frec << ss.str();
		performance_log->log_event(ss.str());

		Eigen::BDCSVD<Eigen::MatrixXd> svd_AY(AY, Eigen::ComputeThinU | Eigen::ComputeThinV);
		p_y = svd_AY.solve(scaled_constraint_diff);
	}
	else
	{
		p_y = ldlt_AY.solve(scaled_constraint_diff); //from Eq 18.19a pp. 538 Nocedal and Wright
	}

	// Check magnitude of p_y
	const double max_p_norm = 1E+6;
	ss.str("");
	if (p_y.norm() > max_p_norm)
	{
		ss.str("");
		ss << "   WARNING: p_y norm too large: " << p_y.norm() << ", scaling down to " << max_p_norm << endl;
		frec << ss.str();
		performance_log->log_event(ss.str());
		p_y *= (max_p_norm / p_y.norm());
	}

	// solve p_null_space
	bool simplified_null_space_approach = false; //TODO: add as ++arg; too much to be an arg?
	Eigen::VectorXd p_z;

	if (Z.cols() > 0)
	{
		Eigen::MatrixXd red_hess = Z.transpose() * G_reg * Z;

		// Check condition number of reduced Hessian
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_red(red_hess);
		double cond_red = eig_red.eigenvalues().maxCoeff() / eig_red.eigenvalues().minCoeff();
		if (cond_red > max_condition_warning)
		{
			ss.str("");
			ss << "   WARNING: Reduced Hessian condition number very large: " << cond_red << endl;
			frec << ss.str();
			performance_log->log_event(ss.str());

			double min_eig_red = eig_red.eigenvalues().minCoeff();
			if (min_eig_red < min_allowed_eig)
			{
				double delta_red = min_allowed_eig - min_eig_red + 1e-6;
				red_hess += delta_red * Eigen::MatrixXd::Identity(red_hess.rows(), red_hess.cols());
				ss.str("");
				ss << "   applied regularization to reduced Hessian, delta = " << delta_red << endl;
				frec << ss.str();
				performance_log->log_event(ss.str());
			}
		}

		//check if positive definite
		Eigen::LDLT<Eigen::MatrixXd> ldlt(red_hess);
		if (ldlt.info() != Eigen::Success || !ldlt.isPositive())
			throw_sqp_error("Reduced Hessian Z^T G Z is not positive definite");

		bool cholesky = false; //TODO: add as ++arg; too much to be an arg?
		if (simplified_null_space_approach)
		{
			ss.str("");
			ss << "   using simplified approach in KKT null space solve..." << endl;
			frec << ss.str();
			performance_log->log_event(ss.str());
			rhs = -Z.transpose() * curved_grad;
			// simplify by removing cross term (or ``partial hessian'') matrix (zTgy), which is approp when approximating hessian (zTgz) (as p_y goes to zero faster than p_z)
			if (cholesky)
			{
				throw_sqp_error("cholesky decomp for null space KKT solve not implemented");
			}
			else
			{
				p_z = ldlt.solve(rhs); //From Eq 18.23 pp. 539 Nocedal and Wright
			}
		}
		else
		{
			Eigen::MatrixXd ZtGY;
			Eigen::VectorXd rhs;
			ZtGY = Z.transpose() * G * Y;
			rhs = (-1. * ZtGY * p_y) - (Z.transpose() * curved_grad); //Eq 18.19b pp. 538 Nocedal and Wright
			if (cholesky)
			{
				//TODO: need to test this
				Eigen::LLT<Eigen::MatrixXd> llt(red_hess);
				if (llt.info() != Eigen::Success)
				{
					throw_sqp_error("Cholesky decomposition failed - matrix not positive definite");
				}
				p_z = llt.solve(rhs);
			}
			else
			{
				p_z = ldlt.solve(rhs); //From Eq 18.19b pp. 538 Nocedal and Wright
			}
		}

		ss.str("");
		if (p_z.norm() > max_p_norm)
		{
			ss << "   WARNING: p_z norm too large: " << p_z.norm() << ", scaling down" << endl;
			frec << ss.str();
			performance_log->log_event(ss.str());
			p_z *= (max_p_norm / p_z.norm());
		}
	}
	else
	{
		p_z = Eigen::VectorXd::Zero(0);
	}

	if (Z.cols() > 0)
		search_d = Y * p_y + Z * p_z; // Eq. 18.18 p. 539 Nocedal and Wright 
	else
		search_d = Y * p_y;

	if (search_d.size() != curved_grad.size())
	{
		throw_sqp_error("search direction vector computation error (in null space KKT solve method)!");
	}

	// compute lagrangian multipliers
	if (simplified_null_space_approach)
	{
		Eigen::BDCSVD<Eigen::MatrixXd> svd_AAT(scaled_constraint_jco * scaled_constraint_jco.transpose(), Eigen::ComputeThinU | Eigen::ComputeThinV);
		lm = svd_AAT.solve(scaled_constraint_jco * curved_grad);
	}
	else
	{
		// Nocedal and Wright pg. 457 and 538
		rhs = Y.transpose() * (curved_grad + G * search_d);
		Eigen::MatrixXd AY_transpose = (scaled_constraint_jco * Y).transpose(); // (A*Y)^T, also m×m
		Eigen::LDLT<Eigen::MatrixXd> ldlt_AYT(AY_transpose);
		if (ldlt_AYT.info() != Eigen::Success)
		{
			Eigen::BDCSVD<Eigen::MatrixXd> svd_AYT(AY_transpose, Eigen::ComputeThinU | Eigen::ComputeThinV);
			lm = svd_AYT.solve(rhs);
		}
		else
		{
			lm = ldlt_AYT.solve(rhs);
		}
	}

	for (int i = 0; i < lm.size(); i++)
		lm[i] /= row_scales[i];

	return pair<Eigen::VectorXd, Eigen::VectorXd>(search_d, lm);
}

pair<Eigen::VectorXd, Eigen::VectorXd> SeqQuadProgram::_kkt_direct(const Eigen::MatrixXd& G, Eigen::MatrixXd& _constraint_jco, Eigen::VectorXd& constraint_diff, Eigen::VectorXd& curved_grad, vector<string>* cnames)
{
	stringstream ss;
	ofstream& frec = file_manager.rec_ofstream();

	ss << "starting KKT direct solve";
	performance_log->log_event(ss.str());

	if (cnames == nullptr || cnames->empty()) 
	{
		Eigen::MatrixXd G_reg = regularize_hessian(G, "direct unconstrained");
		Eigen::LDLT<Eigen::MatrixXd> ldlt(G_reg);
		if (ldlt.info() != Eigen::Success || !ldlt.isPositive()) 
		{
			message(1, "WARNING: LDLT failed for unconstrained case in direct method, using steepest descent");
			Eigen::VectorXd search_d = -curved_grad;
			return pair<Eigen::VectorXd, Eigen::VectorXd>(search_d, Eigen::VectorXd());
		}
		Eigen::VectorXd search_d = ldlt.solve(-curved_grad);
		return pair<Eigen::VectorXd, Eigen::VectorXd>(search_d, Eigen::VectorXd());
	}

	int n = dv_names.size();
	int m = cnames->size();

	if (_constraint_jco.rows() != m || _constraint_jco.cols() != n)
	{
		ss.str("");
		ss << "ERROR: Constraint Jacobian dimension mismatch. Expected " << m << "x" << n << ", got " << _constraint_jco.rows() << "x" << _constraint_jco.cols() << endl;
		throw_sqp_error(ss.str());
	}

	if (constraint_diff.size() != m)
	{
		ss.str("");
		ss << "ERROR: constraint_diff size mismatch. Expected " << m << ", got " << constraint_diff.size() << endl;
		throw_sqp_error(ss.str());
	}

	Eigen::MatrixXd G_reg = regularize_hessian(G, "direct solve");

	//KKT matrix according to Nocedal & Wright Eq. 16.4:
	Eigen::MatrixXd kkt_matrix(n + m, n + m);
	kkt_matrix.topLeftCorner(n, n) = G_reg;
	kkt_matrix.topRightCorner(n, m) = _constraint_jco.transpose();
	kkt_matrix.bottomLeftCorner(m, n) = _constraint_jco;
	kkt_matrix.bottomRightCorner(m, m) = Eigen::MatrixXd::Zero(m, m);

	// form the rhs vector
	// RHS = [-∇f; -c] where c is constraint_diff
	Eigen::VectorXd rhs(n + m);
	rhs.head(n) = -curved_grad;  
	rhs.tail(m) = -constraint_diff; 

	// Solve the KKT system
	// Use LDLT factorization for symmetric indefinite systems (KKT matrix is indefinite)
	Eigen::LDLT<Eigen::MatrixXd> ldlt(kkt_matrix);
	Eigen::VectorXd x;

	if (ldlt.info() != Eigen::Success)
	{
		message(1, "WARNING: LDLT factorization failed for KKT system, falling back to SVD");
		frec << "   WARNING: LDLT factorization failed for KKT system, using SVD" << endl;

		Eigen::BDCSVD<Eigen::MatrixXd> svd(kkt_matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
		x = svd.solve(rhs);
	}
	else
	{
		x = ldlt.solve(rhs);
	}

	// Extract solution components
	Eigen::VectorXd search_d = x.head(n);
	Eigen::VectorXd lm = x.tail(m);

	// Compute residual and error (shared code path)
	Eigen::VectorXd residual = kkt_matrix * x - rhs;
	double kkt_error = residual.norm() / (1.0 + rhs.norm());

	if (kkt_error > 1E-6)
	{
		ss.str("");
		ss << "   WARNING: KKT system solution has high residual: " << kkt_error << endl;
		frec << ss.str();
		message(1, "WARNING: KKT system solution has high residual: ", kkt_error);
	}
	else if (verbose_level >= 2)
	{
		ss.str("");
		ss << "   KKT direct solve completed. Residual: " << kkt_error << endl;
		frec << ss.str();
		message(2, "KKT direct solve completed. Residual: ", kkt_error);
	}

	return pair<Eigen::VectorXd, Eigen::VectorXd>(search_d, lm);
}

pair<Mat, bool> SeqQuadProgram::get_constraint_mat(Parameters& _dv_vals, Observations& _obs_vals, double working_set_tol, const Eigen::VectorXd* lm, vector<string> curr_ws)
{
	if (use_ensemble_grad) 
	{
		Parameters decvar = _dv_vals.get_subset(dv_names.begin(), dv_names.end());

		if (constraints.get_use_chance() && (pest_scenario.get_pestpp_options().get_sqp_risk() != 0.5))
		{
			Observations base_obs = _obs_vals;
			if (oe.shape().first > 0)
			{
				vector<string> obs_names = oe.get_var_names();
				Eigen::VectorXd base_vec = oe.get_real_vector(BASE_REAL_NAME);
				base_obs.update_without_clear(obs_names, base_vec);
			}

			Observations shifted_obs = constraints.get_chance_shifted_constraints(base_obs, oe, sqp_risk);
			vector<string> constraint_names = constraints.get_obs_constraint_names();
			for (auto& name : constraint_names)
			{
				if (shifted_obs.find(name) != shifted_obs.end())
				{
					_obs_vals.update_rec(name, shifted_obs.get_rec(name));
				}
			}
		}

		return constraints.get_working_set_constraint_matrix(decvar, _obs_vals, dv, oe, true, lm, curr_ws, (working_set_tol));
	}
	else
	{
		message(2, "getting working set constraint matrix");
		return constraints.get_working_set_constraint_matrix(_dv_vals, _obs_vals, jco, true, lm, (working_set_tol));
	}
}

bool SeqQuadProgram::trust_region_step(Parameters& current_dv_values, Eigen::VectorXd grad)
{
	prev_ctl_dv_values = trial_ctl_dv_values; //saving a copy for BFGS later

	// Algorithm 4.1, pp. 68-69 Nocedal and Wright
	double current_obj = get_obj_value(trial_ctl_dv_values, trial_obs);

	// Solve the trust region subproblem to get step
	Eigen::VectorXd p = solve_trust_region_subproblem_dogleg(hessian.get_matrix(), grad, trust_radius);
	double step_norm = p.norm();

	// Create trial point
	Parameters trial_dv_values = trial_ctl_dv_values;
	Eigen::VectorXd trial_vec = current_dv_values.get_data_eigen_vec(dv_names);
	trial_vec += p;
	trial_dv_values.update_without_clear(dv_names, trial_vec);

	// Enforce bounds on trial point
	ParameterEnsemble trial_pe(&pest_scenario, &rand_gen);
	trial_pe.reserve({ "trial" }, dv_names);
	trial_pe.update_real_ip("trial", trial_vec);
	trial_pe.enforce_bounds(performance_log, false);

	// Get bounded trial vector
	trial_vec = trial_pe.get_real_vector("trial");
	trial_dv_values.update_without_clear(dv_names, trial_vec);

	// Recalculate actual step after bounds enforcement
	p = trial_vec - current_dv_values.get_data_eigen_vec(dv_names);
	step_norm = p.norm();

	// Run the model for the trial point
	ObservationEnsemble trial_oe = run_candidate_ensemble(trial_pe);

	// If model run failed, reject step and reduce trust region
	if (trial_oe.shape().first == 0) {
		trust_radius = max(trust_radius_min, 0.5 * trust_radius);
		stringstream ss;
		ss << "model run failed for trust region step, reducing radius to " << trust_radius;
		message(1, ss.str());
		return false;
	}

	// Get trial observations
	trial_obs = pest_scenario.get_ctl_observations();
	Eigen::VectorXd obs_vec = trial_oe.get_real_vector("trial");
	trial_obs.update(trial_oe.get_var_names(), obs_vec);

	cout << "trial_dv: " << trial_pe.get_eigen() << endl;

	// Calculate actual and predicted reduction
	double actual_reduction = compute_actual_reduction(trial_dv_values, trial_obs);
	double predicted_reduction = compute_predicted_reduction(p, grad);

	// Define trust region parameters
	const double eta1 = 0.25;  // Threshold for accepting step
	const double eta2 = 0.75;  // Threshold for very successful step

	// Calculate ratio of actual to predicted reduction
	double rho = 0.0;
	if (abs(predicted_reduction) > 1e-10) {
		rho = actual_reduction / predicted_reduction;
	}

	// Update trust region radius and accept/reject step based on ratio
	if (rho < eta1) {
		// Unsuccessful step - reduce trust region
		trust_radius = max(trust_radius_min, 0.25 * trust_radius);
		stringstream ss;
		ss << "rejected trust region step, rho=" << rho
			<< ", new radius=" << trust_radius;
		message(1, ss.str());

		trial_ctl_dv_values = trial_dv_values;
		Parameters new_grad = calc_gradient_vector(trial_ctl_dv_values);
		current_grad_vector = new_grad;
		make_gradient_runs(trial_ctl_dv_values, trial_obs);
		//update_hessian_and_grad_vector();

		if (trust_radius == trust_radius_min)
			return true;
		else
			return false;
	}
	else {
		// Accept step
		current_obs = trial_obs;
		current_ctl_dv_values = trial_dv_values;

		// Update trust region radius based on step quality
		if (rho > eta2 && step_norm > 0.8 * trust_radius) {
			// Very successful step and we're near boundary - increase radius
			trust_radius = min(2.0 * trust_radius, trust_radius_max);
			stringstream ss;
			ss << "very successful trust region step with rho=" << rho
				<< ", increasing radius to " << trust_radius;
			message(1, ss.str());
		}
		else {
			// Moderately successful step - keep same radius
			stringstream ss;
			ss << "accepted trust region step with rho=" << rho
				<< ", keeping radius=" << trust_radius;
			message(1, ss.str());
		}

		// Update objective history for non-monotone criteria
		double new_obj = get_obj_value(current_ctl_dv_values, current_obs);
		previous_obj_values.push_back(new_obj);
		if (previous_obj_values.size() > memory_length) {
			previous_obj_values.erase(previous_obj_values.begin());
		}

		return true;
	}
}


FilterRec SeqQuadProgram::trust_region_step(Eigen::VectorXd& grad, map<string, double> current_obj_ens, map<string, vector<string>>& cnames_en,
	map<string, Eigen::MatrixXd>& constraint_jco_en, ParameterEnsemble* dvs_subset, bool recalc)
{
	stringstream ss;
	if (!recalc)
	{
		ss << "performing line search on ensemble subset:";
		for (auto& n : dvs_subset->get_real_names())
			ss << " " << n << ",";
		ss << " BASE" << endl;
		message(1, ss.str());
	}
	else
	{
		ss << "performing line search on realizations:";
		const auto& n = dvs_subset->get_real_names();
		for (size_t i = 0; i < n.size(); i++)
		{
			ss << ' ' << n[i];
			if (i + 1 != n.size())
				ss << ',';
		}
		message(1, ss.str());
	}

	ParameterEnsemble dv_candidates(&pest_scenario, &rand_gen);
	dv_candidates.set_trans_status(ParameterEnsemble::transStatus::NUM);

	map<double, map<string, string>> radius_real_map;
	vector<string> real_names;
	vector<double> radius_scale_vals;
	for (auto sf : pest_scenario.get_pestpp_options().get_sqp_alpha_mults())
	{
		if (sf > 0)
			radius_scale_vals.push_back(sf * BASE_SCALE_FACTOR);
	}

	if (SOLVE_EACH_REAL)
	{
		for (double radius_scale : radius_scale_vals)
		{
			for (const auto& real_name : dv.get_real_names())
			{
				ss.str("");
				ss << "cand_tr" << real_name << "_rs:" << left << setw(8) << setprecision(3) << radius_scale;
				real_names.push_back(ss.str());
			}
		}
	}
	else
	{
		if (dvs_subset != nullptr)
		{
			ParameterEnsemble d;
			const auto& names = dvs_subset->get_real_names();
			if (find(names.begin(), names.end(), BASE_REAL_NAME) == names.end() && !recalc)
			{
				d.reserve(vector<string>{ BASE_REAL_NAME }, dv_names);
				d.add_2_row_ip(BASE_REAL_NAME, dv.get_real_vector(BASE_REAL_NAME));
				dvs_subset->append_other_rows(d);
			}

			for (const auto& real_name : dvs_subset->get_real_names())
			{
				for (double radius_scale : radius_scale_vals)
				{
					ss.str("");
					ss << "cand_tr" << real_name << "_rs:" << left << setw(8) << setprecision(3) << radius_scale;
					real_names.push_back(ss.str());
					radius_real_map[radius_scale][real_name] = ss.str();
				}
			}
		}
		else
			throw_sqp_error("use_ensemble_grad is true but subset dv ensemble is null");
	}
	dv_candidates.reserve(real_names, dv_names);

	vector<double> used_radius_scales;
	map<string, double> real_rs_map;
	Eigen::MatrixXd H = hessian.get_matrix();

	for (size_t i = 0; i < radius_scale_vals.size(); ++i)
	{
		double radius_scale = radius_scale_vals[i];
		message(1, "starting trust region calcs for radius scale", radius_scale);

		double candidate_trust_radius = trust_radius * radius_scale;
		candidate_trust_radius = max(trust_radius_min, min(candidate_trust_radius, trust_radius_max));

		if (dvs_subset != nullptr)
		{
			map<string, Eigen::VectorXd> candidate_steps;
			vector<string> short_steps;

			for (const auto& real_name : dvs_subset->get_real_names())
			{
				Eigen::VectorXd real_grad = grad;

				Eigen::MatrixXd A_real = Eigen::MatrixXd::Zero(0, H.cols());
				auto a_it = constraint_jco_en.find(real_name);
				if (a_it != constraint_jco_en.end())
					A_real = a_it->second;

				Eigen::VectorXd step = solve_constrained_trust_region_step(H, real_grad, A_real, candidate_trust_radius);
				Eigen::VectorXd current_dv_vec = dvs_subset->get_real_vector(real_name);

				ParameterInfo par_info = pest_scenario.get_ctl_parameter_info();
				Parameters lbnd = par_info.get_low_bnd(dv_names);
				Parameters ubnd = par_info.get_up_bnd(dv_names);
				ParamTransformSeq par_transform = pest_scenario.get_base_par_tran_seq();
				par_transform.ctl2numeric_ip(lbnd);
				par_transform.ctl2numeric_ip(ubnd);

				Eigen::VectorXd alpha = Eigen::VectorXd::Ones(dv_names.size());
				for (int j = 0; j < dv_names.size(); ++j)
				{
					double x = current_dv_vec(j);
					double d = step(j);
					double lb = lbnd[dv_names[j]];
					double ub = ubnd[dv_names[j]];

					if (abs(d) < 1e-10)
						continue;

					if (d > 0)
						alpha(j) = min(alpha(j), (ub - x) / d);
					else if (d < 0)
						alpha(j) = min(alpha(j), (x - lb) / (-d));
				}
				double min_alpha = alpha.minCoeff();
				if (min_alpha < 1.0)
					step *= min_alpha;

				candidate_steps[real_name] = step;
				if (step.squaredNorm() < 1.0e-10)
					short_steps.push_back(real_name);
			}

			if (!short_steps.empty())
			{
				ss.str("");
				ss << " very short trust-region steps for radius scale " << radius_scale << " realizations:" << short_steps;
				message(1, ss.str());
			}

			for (const auto& real_name : dvs_subset->get_real_names())
			{
				Eigen::VectorXd dv_upgrade = dvs_subset->get_real_vector(real_name);
				dv_upgrade += candidate_steps.at(real_name);

				string candidate_rname = radius_real_map.at(radius_scale).at(real_name);
				dv_candidates.update_real_ip(candidate_rname, dv_upgrade);
				real_rs_map[candidate_rname] = radius_scale;
				used_radius_scales.push_back(radius_scale);
			}
		}
		message(1, "finished trust region calcs for radius scale:", radius_scale);
	}

	if (pest_scenario.get_pestpp_options().get_ies_debug_upgrade_only())
	{
		message(0, "ies_debug_upgrade_only is true, exiting");
		throw_sqp_error("ies_debug_upgrade_only is true, exiting");
	}

	dv_candidates.enforce_bounds(performance_log, false);
	ss.str("");
	ss << file_manager.get_base_filename() << "." << iter << ".dv_candidates.csv";

	if (!recalc)
		dv_to_save = dv_candidates;
	else
	{
		auto save = dv_to_save.get_real_map();
		for (const auto& rname : dv_candidates.get_real_names())
		{
			Eigen::VectorXd row = dv_candidates.get_real_vector(rname);
			if (save.find(rname) != save.end())
			{
				dv_to_save.update_real_ip(rname, row);
			}
			else
			{
				throw_sqp_error("active set recalc has new candidate not in original candidate set:" + rname);
				//dv_to_save.append(rname, row);
				//save[rname] = dv_to_save.get_real_map().at(rname);
			}
		}
	}

	dv_candidates.to_csv(ss.str());
	message(1, "saved trust region candidate ensemble to:", ss.str());

	Eigen::VectorXd v1, v2;
	double d;
	vector<string> drop;
	set<int> jvals;
	for (int i = 0; i < dv_candidates.shape().first; ++i)
	{
		v1 = dv_candidates.get_real_vector(i);
		for (int j = i + 1; j < dv_candidates.shape().first; ++j)
		{
			v2 = dv_candidates.get_real_vector(j);
			d = (v2 - v1).squaredNorm();
			if ((abs(d) < 1E-12) && (jvals.find(j) == jvals.end()))
			{
				message(1, "duplicate trust-region candidates:",
					vector<string>{ dv_candidates.get_real_names()[i],
					dv_candidates.get_real_names()[j] });
				drop.push_back(dv_candidates.get_real_names()[j]);
				jvals.emplace(j);
			}
		}
	}
	if (!drop.empty())
	{
		dv_candidates.drop_rows(drop, true);
		used_radius_scales.clear();
		for (const auto& real_name : dv_candidates.get_real_names())
			used_radius_scales.push_back(real_rs_map.at(real_name));
	}

	message(0, "running candidate dv/pars");
	ObservationEnsemble oe_candidates = run_candidate_ensemble(dv_candidates);

	if (!recalc)
		oe_to_save = oe_candidates;
	else
	{
		auto save = oe_to_save.get_real_map();
		for (const auto& rname : oe_candidates.get_real_names())
		{
			Eigen::VectorXd row = oe_candidates.get_real_vector(rname);
			if (save.find(rname) != save.end())
			{
				oe_to_save.update_real_ip(rname, row);
			}
			else
			{
				throw_sqp_error("active set recalc has new candidate not in original candidate set:" + rname);
				//oe_to_save.append(rname, row);
				//save[rname] = oe_to_save.get_real_map().at(rname);
			}
		}
	}

	ObservationEnsemble combined_oe_candidates = combine_obs_and_pi(oe_to_save, dv_to_save);
	ss.str("");
	ss << file_manager.get_base_filename() << "." << iter << ".oe_candidates.csv";
	combined_oe_candidates.to_csv(ss.str());
	message(1, "saved trust region candidate ensemble obs to:", ss.str());

	ParameterEnsemble dv_all = dv;
	ObservationEnsemble oe_all = oe;
	dv_all.append_other_rows(dv_candidates);
	oe_all.append_other_rows(oe_candidates);

	bool success;
	FilterRec result = pick_upgrade_and_update_current(dv_all, oe_all, cma_reset_archive, true);
	if (result.viol_val > 0.0)
		success = false;
	else
		success = true;

	if (success)
	{
		Eigen::VectorXd current_vec = prev_ctl_dv_values.get_data_eigen_vec(dv_names);
		Eigen::VectorXd step_taken =
			current_ctl_dv_values.get_data_eigen_vec(dv_names) - current_vec;

		double predicted = compute_predicted_reduction(step_taken, grad);
		double actual = get_obj_value(prev_ctl_dv_values, prev_ctl_dv_values == current_ctl_dv_values ? current_obs : current_obs)
			- get_obj_value(current_ctl_dv_values, current_obs);

		double rho = 0.0;
		if (abs(predicted) > 1e-10)
			rho = actual / predicted;

		const double eta1 = 0.25;
		const double eta2 = 0.75;

		if (rho > eta2 && step_taken.norm() > 0.8 * trust_radius)
		{
			trust_radius = min(2.0 * trust_radius, trust_radius_max);
			ss.str("");
			ss << "very successful trust region step with rho=" << rho << ", increasing radius to " << trust_radius;
			message(1, ss.str());
		}
		else if (rho > eta1)
		{
			ss.str("");
			ss << "accepted trust region step with rho=" << rho << ", keeping radius=" << trust_radius;
			message(1, ss.str());
		}
		else
		{
			trust_radius = max(trust_radius_min, 0.25 * trust_radius);
			ss.str("");
			ss << "trust region step accepted by filter but rho=" << rho << ", reducing radius to " << trust_radius;
			message(1, ss.str());
		}
	}
	else
	{
		trust_radius = max(trust_radius_min, 0.5 * trust_radius);
		message(1, "no trust region candidate accepted by filter, reducing radius to", trust_radius);
	}
	return result;
	
}


Eigen::VectorXd SeqQuadProgram::solve_trust_region_subproblem_dogleg(const Eigen::MatrixXd& B, const Eigen::VectorXd& g, double radius)
{
	// Dogleg method for trust region subproblem (Algorithm 4.3 in Nocedal & Wright)
	int n = g.size();

	// Step 1: Compute the Cauchy point (steepest descent direction)
	double gBg = g.transpose() * B * g;
	double alpha;

	if (gBg <= 0) {
		// If Hessian is not positive definite, use a simple scaling
		alpha = radius / g.norm();
	}
	else {
		alpha = pow(g.norm(), 3) / (radius * gBg);

		// Limit step to trust region boundary if needed
		Eigen::VectorXd p_sd = -alpha * g;
		if (p_sd.norm() > radius) {
			alpha = radius / g.norm();
		}
	}

	Eigen::VectorXd p_cauchy = -alpha * g;

	// If Cauchy point is at the boundary, return it
	if (p_cauchy.norm() >= radius) {
		return p_cauchy;
	}

	// Step 2: Compute the Newton point (full step)
	Eigen::VectorXd p_newton;

	Eigen::LDLT<Eigen::MatrixXd> ldlt(B);
	if (ldlt.info() == Eigen::Success) {
		p_newton = -ldlt.solve(g);
	}
	else {
		// If decomposition fails, use a regularized version
		double lambda = 1e-6;
		Eigen::MatrixXd B_reg = B + lambda * Eigen::MatrixXd::Identity(n, n);
		Eigen::LDLT<Eigen::MatrixXd> ldlt_reg(B_reg);
		p_newton = -ldlt_reg.solve(g);
	}

	// If Newton point is within trust region, return it
	if (p_newton.norm() <= radius) {
		return p_newton;
	}

	// Step 3: Compute the dogleg path intersection with trust region boundary
	// Find tau where ||p_cauchy + tau * (p_newton - p_cauchy)|| = radius
	Eigen::VectorXd d = p_newton - p_cauchy;

	// Solve quadratic equation: ||p_cauchy + tau*d||^2 = radius^2
	double a = d.squaredNorm();
	double b = 2 * p_cauchy.dot(d);
	double c = p_cauchy.squaredNorm() - radius * radius;

	// Quadratic formula: tau = (-b + sqrt(b^2 - 4ac)) / (2a)
	double discriminant = b * b - 4 * a * c;
	double tau = (-b + sqrt(discriminant)) / (2 * a);

	// Compute the dogleg point
	Eigen::VectorXd p_dogleg = p_cauchy + tau * d;

	return p_dogleg;
}

double SeqQuadProgram::compute_actual_reduction(Parameters& trial_dv_values, Observations& trial_obs)
{
	double current_obj = get_obj_value(current_ctl_dv_values, current_obs);
	double trial_obj = get_obj_value(trial_dv_values, trial_obs);
	return current_obj - trial_obj;
}

double SeqQuadProgram::compute_predicted_reduction(const Eigen::VectorXd& step,
	const Eigen::VectorXd& grad)
{
	//predicted reduction using quadratic model Eq 4.2, p. 68 Nocedal and Wright
	double linear_term = -grad.dot(step);
	double quadratic_term = -0.5 * step.dot(hessian.get_matrix() * step);
	return linear_term + quadratic_term;
}

FilterRec SeqQuadProgram::line_search(map<string, Eigen::VectorXd>& search_d_map, Eigen::VectorXd& grad, map<string, double> current_obj_ens, ParameterEnsemble* dvs_subset, bool recalc)
{
	stringstream ss;
	const auto& dv_real_names = dv.get_real_names();
	vector <string> subset_real_names = dvs_subset->get_real_names();

	if (!recalc)
	{
		ss << "performing line search on ensemble subset:";
		for (auto& n : subset_real_names)
			ss << " " << n << ",";
		ss << " BASE" << endl;
		message(1, ss.str());
	}
	else
	{
		ss << "performing line search on realizations:";
		const auto& n = subset_real_names;
		for (size_t i = 0; i < n.size(); i++)
		{
			ss << ' ' << n[i];
			if (i + 1 != n.size())
				ss << ',';
		}
		message(1, ss.str());
	}
	
	ParameterEnsemble dv_candidates(&pest_scenario, &rand_gen);
	dv_candidates.set_trans_status(ParameterEnsemble::transStatus::NUM);

	
	vector<string> real_names;
	vector<double> scale_vals;
	
	/*for (auto& sf : pest_scenario.get_pestpp_options().get_sqp_alpha_mults())
	{
		scale_vals.push_back(sf * BASE_SCALE_FACTOR);
	}*/

	vector<double> alpha_mults = pest_scenario.get_pestpp_options().get_sqp_alpha_mults();
	sort(alpha_mults.begin(), alpha_mults.end());

	for (auto it = alpha_mults.begin(); it != alpha_mults.end(); it++)
	{
		if (it == alpha_mults.begin())
		{
			scale_vals.push_back(*it * BASE_SCALE_FACTOR);
		}
		else
		{
			double prev = *(it - 1);
			double curr = *it;
			scale_vals.push_back(prev + (curr - prev) * BASE_SCALE_FACTOR);
		}
	}

	sv_lineage_map.clear();
	if ((use_ensemble_grad) && (SOLVE_EACH_REAL))
	{
		for (auto sv : scale_vals)
		{
			for (auto& rname : dv_real_names)
			{
				ss.str("");
				ss << "cand_" << rname << "_sv:" << left << setw(8) << setprecision(3) << sv;
				real_names.push_back(ss.str()); 
				ls_parent_map[ss.str()] = rname;
			}
		}

	}
	else if (use_ensemble_grad)
	{
		if (dvs_subset != nullptr) 
		{
			ParameterEnsemble d;
			map<string,int> rmap = dvs_subset->get_real_map();
			if ((rmap.find(BASE_REAL_NAME) == rmap.end()) && !recalc)
			{
				if (subset_real_names.size() == 0)
				{
					dvs_subset->reserve(vector<string>{ BASE_REAL_NAME }, dv_names);
					Eigen::VectorXd real_vec = dv.get_real_vector(BASE_REAL_NAME);
					dvs_subset->update_real_ip(BASE_REAL_NAME, real_vec);
					subset_real_names.push_back(BASE_REAL_NAME);
				}
				else
				{
					d.reserve(vector<string>{ BASE_REAL_NAME }, dv_names);
					d.add_2_row_ip(BASE_REAL_NAME, dv.get_real_vector(BASE_REAL_NAME));
					dvs_subset->append_other_rows(d);
					subset_real_names.push_back(BASE_REAL_NAME);
				}
			}
			
			for (auto& rname : subset_real_names) 
			{
				for (auto sv : scale_vals) 
				{
					ss.str("");
					ss << "cand_" << rname << "_sv:" << left << setw(8) << setprecision(3) << sv;
					real_names.push_back(ss.str());
					sv_lineage_map[sv][rname] = ss.str();
					ls_parent_map[ss.str()] = rname;
				}
			}
		}
		else
			throw_sqp_error("use_ensemble_grad is true but subset dv ensemble is null");
	}
	dv_candidates.reserve(real_names, dv_names);

	vector<double> used_scale_vals;
	cname_sf_map.clear();
	step_length_map.clear();

	for (int i = 0;i < scale_vals.size();i++)
	{
		double scale_val = scale_vals[i];
		ss.str("");
		ss << "starting calcs for scaling factor" << scale_val;
		message(1, "starting lambda calcs for scaling factor", scale_val);

		map<string, Eigen::VectorXd> scaled_sdir_map;
		vector<string> short_upgrades;
		for (auto& sd : search_d_map)
		{
			scaled_sdir_map[sd.first] = sd.second * scale_val;
			if (scaled_sdir_map[sd.first].squaredNorm() < 1.0E-10)
				short_upgrades.push_back(sd.first);
		}
		if (short_upgrades.size() > 0)
		{
			ss.str("");
			ss << "there are very short upgrades for scale value " << scale_val << " for realizations: " << short_upgrades;
			message(1, ss.str());
		}
		
		for (const auto& real_name : dvs_subset->get_real_names())
		{
			Eigen::VectorXd dv_upgrade = dvs_subset->get_real_vector(real_name);
			dv_upgrade += scaled_sdir_map[real_name];			
			
			string cand_rname = sv_lineage_map[scale_val][real_name];
			dv_candidates.update_real_ip(cand_rname, dv_upgrade);
			step_length_map[cand_rname] = scaled_sdir_map[real_name];
			cname_sf_map[cand_rname] = scale_val;
			used_scale_vals.push_back(scale_val);
		}
			
		ss.str("");
		message(1, "finished calcs for scaling factor:", scale_val);

	}

	if (pest_scenario.get_pestpp_options().get_sqp_enforce_bounds())
		dv_candidates.enforce_bounds(performance_log, false);

	ss.str("");
	ss << file_manager.get_base_filename() << "." << iter << ".dv_candidates.csv";
	
	const vector<string> cand_real_names = dv_candidates.get_real_names();
	if (!recalc)
		dv_to_save = dv_candidates;
	else
	{
		auto save = dv_to_save.get_real_map();
		for (const auto& rname : cand_real_names)
		{
			Eigen::VectorXd row = dv_candidates.get_real_vector(rname);
			if (save.find(rname) != save.end())
			{
				dv_to_save.update_real_ip(rname, row);
			}
			else
			{
				throw_sqp_error("active set recalc has new candidate not in original candidate set:" + rname);
			}
		}
	}

	Eigen::VectorXd v1, v2;
	double d;
	vector<string> drop;
	set<int> jvals;
	for (int i = 0;i < dv_candidates.shape().first; i++)
	{
		v1 = dv_candidates.get_real_vector(i);
		for (int j = i + 1;j < dv_candidates.shape().first; j++) {
			v2 = dv_candidates.get_real_vector(j);
			d = (v2 - v1).squaredNorm();
			if ((abs(d) < 1E-12) && (jvals.find(j) == jvals.end())) {
				message(1, "duplicate candidates:", vector<string>{cand_real_names[i], cand_real_names[j]});
				drop.push_back(cand_real_names[j]);
				jvals.emplace(j);
			}
		}

	}
	if (drop.size() > 0)
	{
		message(1, "dropping the following duplicate candidates: ", drop);
		dv_candidates.drop_rows(drop, true);
		used_scale_vals.clear();
		for (auto& rname : cand_real_names)
		{
			used_scale_vals.push_back(cname_sf_map.at(rname));
		}
	}
	message(0, "running candidate dv/pars");
	ObservationEnsemble oe_candidates = run_candidate_ensemble(dv_candidates);

	if (!recalc)
		oe_to_save = oe_candidates;
	else
	{
		auto save = oe_to_save.get_real_map();
		for (const auto& rname : oe_candidates.get_real_names())
		{
			Eigen::VectorXd row = oe_candidates.get_real_vector(rname);
			if (save.find(rname) != save.end() || find(drop.begin(), drop.end(), rname) != drop.end())
			{
				oe_to_save.update_real_ip(rname, row);
			}
			else
			{
				oe_to_save.append(rname, row);
				save[rname] = oe_to_save.get_real_map().at(rname);
			}
		}
	}

	const int num_intermediate_pts = pest_scenario.get_pestpp_options().get_sqp_num_refined_search_pts();
	if (!recalc && scale_vals.size() > 2 && num_intermediate_pts > 0)
	{
		message(1, "performing adaptive refinement around promising candidates");

		map<string, double> omap = get_obj_map(dv, oe);
		map<string, double> omap_cand = get_obj_map(dv_candidates, oe_candidates);
		map<string, map<string, double>> vmap = constraints.get_ensemble_violations_map(dv, oe, 0.0, true);
		map<string, map<string, double>> vmap_cand = constraints.get_ensemble_violations_map(dv_candidates, oe_candidates, 0.0, true);

		map<string, pair<string, double>> best_candidates; 
		for (const auto& rname : subset_real_names)
		{
			double best_obj = omap[rname];
			double best_viol = 0.0;
			for (const auto& v : vmap[rname])
				best_viol += v.second;
			string best_cand = rname;

			for (const auto& cand_name : cand_real_names)
			{
				if (ls_parent_map[cand_name] != rname)
					continue;

				double obj_val = omap_cand[cand_name];
				double viol_sum = 0.0;
				for (const auto& v : vmap_cand[cand_name])
					viol_sum += v.second;

				bool is_better = false;
				if (viol_sum <= 1E-6)
				{
					if (best_viol > 1E-6)
						is_better = true;
					else if (obj_sense == "minimize" && obj_val < best_obj)
						is_better = true;
					else if (obj_sense == "maximize" && obj_val > best_obj)
						is_better = true;
				}
				else if (viol_sum < best_viol)
				{
					is_better = true;
				}

				if (is_better)
				{
					best_obj = obj_val;
					best_viol = viol_sum;
					best_cand = cand_name;
				}
			}

			if (find(subset_real_names.begin(), subset_real_names.end(), best_cand) != subset_real_names.end())
				best_candidates[rname] = { best_cand, 0.0 };
			else
				best_candidates[rname] = { best_cand, cname_sf_map[best_cand] };
		}

		ParameterEnsemble dv_intermediate(&pest_scenario, &rand_gen);
		dv_intermediate.set_trans_status(ParameterEnsemble::transStatus::NUM);
		vector<string> intermediate_cand_names;
		scale_vals.push_back(0.0);
		sort(scale_vals.begin(), scale_vals.end());

		for (const auto& bc : best_candidates)
		{
			string parent_name = bc.first;
			double best_scale = bc.second.second;

			auto scale_it = find(scale_vals.begin(), scale_vals.end(), best_scale);
			if (scale_it == scale_vals.end())
				continue;

			size_t best_idx = distance(scale_vals.begin(), scale_it);
			if (best_idx == 0)
				generate_intermediate_candidates(parent_name, -scale_vals[1], 0.0, num_intermediate_pts, dvs_subset, search_d_map, dv_intermediate, intermediate_cand_names);
			else if (best_idx == scale_vals.size() - 1)
				generate_intermediate_candidates(parent_name, scale_vals[best_idx], best_scale * (1.0 + BASE_SCALE_FACTOR), num_intermediate_pts, dvs_subset, search_d_map, dv_intermediate, intermediate_cand_names);
			
			if (best_idx > 0) 
				generate_intermediate_candidates(parent_name, scale_vals[best_idx - 1], best_scale,	num_intermediate_pts, dvs_subset, search_d_map, dv_intermediate, intermediate_cand_names);
			
			if (best_idx < scale_vals.size() - 1)
				generate_intermediate_candidates(parent_name, scale_vals[best_idx + 1], best_scale, num_intermediate_pts, dvs_subset, search_d_map, dv_intermediate, intermediate_cand_names);
			
		}

		if (!intermediate_cand_names.empty())
		{
			ss.str("");
			ss << "evaluating " << intermediate_cand_names.size() << " intermediate candidates";
			message(1, ss.str());

			if (pest_scenario.get_pestpp_options().get_sqp_enforce_bounds())
				dv_intermediate.enforce_bounds(performance_log, false);

			ObservationEnsemble oe_intermediate = run_candidate_ensemble(dv_intermediate);

			dv_candidates.append_other_rows(dv_intermediate);
			oe_candidates.append_other_rows(oe_intermediate);

			if (!recalc)
			{
				dv_to_save.append_other_rows(dv_intermediate, true);
				oe_to_save.append_other_rows(oe_intermediate, true);
			}
		}
	}

	ss.str("");
	ss << file_manager.get_base_filename() << "." << iter << ".dv_candidates.csv";
	dv_to_save.to_csv(ss.str());
	message(1, "saved candidate decvar/parameter ensemble to: ", ss.str());

	ObservationEnsemble combined_oe_candidates = combine_obs_and_pi(oe_to_save, dv_to_save);
	ss.str("");
	ss << file_manager.get_base_filename() << "." << iter << ".oe_candidates.csv";
	combined_oe_candidates.to_csv(ss.str());
	message(1, "saved candidate ensemble obs to: ", ss.str());

	ParameterEnsemble dv_all = dv;
	ObservationEnsemble oe_all = oe;
	dv_all.append_other_rows(dv_candidates);
	oe_all.append_other_rows(oe_candidates);
	
	return pick_upgrade_and_update_current(dv_all , oe_all, cma_reset_archive, true, dvs_subset, recalc);
}

void SeqQuadProgram::generate_intermediate_candidates(const string& parent_name, double start_scale, double end_scale, int num_points, 
ParameterEnsemble* dvs_subset, const map<string, Eigen::VectorXd>& search_d_map, ParameterEnsemble& dv_intermediate, vector<string>& intermediate_cand_names)
{
	if (num_points <= 0 || fabs(end_scale - start_scale) < 1E-10)
		return;

	Eigen::VectorXd parent_dv = dvs_subset->get_real_vector(parent_name);
	Eigen::VectorXd search_dir = search_d_map.at(parent_name);
	double scale_range = end_scale - start_scale;

	for (int i = 1; i <= num_points; i++)
	{
		double alpha = static_cast<double>(i) / (num_points + 1.0);
		double intermediate_scale = start_scale + alpha * scale_range;

		stringstream ss_ref;
		ss_ref << "cand_" << parent_name << "_sv:" << left << setw(8) << setprecision(6) << intermediate_scale;
		string intermediate_name = ss_ref.str();

		Eigen::VectorXd intermediate_dv = parent_dv + search_dir * intermediate_scale;

		if (intermediate_cand_names.empty())
		{
			dv_intermediate.reserve(vector<string>{intermediate_name}, dv_names);
			dv_intermediate.update_real_ip(intermediate_name, intermediate_dv);
		}
		else
			dv_intermediate.append(intermediate_name, intermediate_dv);

		intermediate_cand_names.push_back(intermediate_name);
		cname_sf_map[intermediate_name] = intermediate_scale;
		ls_parent_map[intermediate_name] = parent_name;
		step_length_map[intermediate_name] = search_dir * intermediate_scale;
		sv_lineage_map[intermediate_scale][parent_name] = intermediate_name;
	}
}

Ensemble SeqQuadProgram::get_pi_ensemble(ParameterEnsemble& _dv, vector<string>& pinames)
{
	vector<string> parnames = _dv.get_var_names();
	vector<string> realnames = _dv.get_real_names();
	PriorInformation pi_info = pest_scenario.get_prior_info();
	
	Ensemble pioe(&pest_scenario);
	pioe.reserve(realnames, pinames);
	Eigen::MatrixXd piX(realnames.size(), pinames.size());

	Parameters pars = pest_scenario.get_ctl_parameters();
	ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();

	int i = 0;
	for (auto& real : realnames)
	{
		Eigen::VectorXd parvals = _dv.get_real_vector(real);
		pars.update_without_clear(parnames, parvals);

		if (_dv.get_trans_status() == ParameterEnsemble::transStatus::NUM)
			pts.numeric2ctl_ip(pars);

		int j = 0;
		for (auto& piname : pinames)
		{
			PriorInformationRec pi_rec = pi_info.get_pi_rec(piname);
			pair<double, double> pi_sim_resid = pi_rec.calc_sim_and_resid(pars);
			piX(i, j) = pi_sim_resid.first; 
			j++;
		}
		i++;
	}

	pioe.from_eigen_mat(piX, realnames, pinames);
	return pioe;
}


ObservationEnsemble SeqQuadProgram::combine_obs_and_pi(ObservationEnsemble& _oe, ParameterEnsemble& _pe)
{
	vector<string> pinames = pest_scenario.get_ctl_ordered_pi_names();
	ObservationEnsemble combined_oe = _oe;

	if (pinames.size() > 0)
	{
		Ensemble pi_oe = get_pi_ensemble(_pe, pinames);

		Eigen::MatrixXd obs_mat = _oe.get_eigen();
		Eigen::MatrixXd pi_mat = pi_oe.get_eigen(_oe.get_real_names(),pinames);

		Eigen::MatrixXd combined_mat(obs_mat.rows(), obs_mat.cols() + pi_mat.cols());
		combined_mat << obs_mat, pi_mat;

		vector<string> combined_var_names = _oe.get_var_names();
		vector<string> pi_var_names = pi_oe.get_var_names();
		combined_var_names.insert(combined_var_names.end(), pi_var_names.begin(), pi_var_names.end());

		combined_oe = ObservationEnsemble(&pest_scenario, &rand_gen);
		combined_oe.from_eigen_mat(combined_mat, _oe.get_real_names(), combined_var_names);
	}

	return combined_oe;
}

pair<Eigen::VectorXd, Eigen::VectorXd> SeqQuadProgram::calc_search_direction_vector(Parameters& _current_dv_values, Observations& _current_obs_values, Eigen::VectorXd& grad_vector, Eigen::MatrixXd* _constraint_jco, vector<string>* _cnames)
{
	Eigen::VectorXd search_d, lm;
	vector<string> Cnames = _cnames != nullptr ? *_cnames : this->cnames;
	Eigen::MatrixXd constr_jco = _constraint_jco != nullptr ? *_constraint_jco : base_constraint_jco;
	pair<Eigen::VectorXd, Eigen::VectorXd> x;
	stringstream ss;
	ofstream& frec = file_manager.rec_ofstream();

	if (Cnames.size() > 0)
	{
		Eigen::VectorXd constraint_diff = constraints.get_obs_resid_constraint_vectors(_current_dv_values, _current_obs_values, Cnames).second;
		for (int i = 0;i < Cnames.size();i++) {
			if (constraint_sense[Cnames[i]] == "less_than")
			{
				constr_jco.row(i) *= -1;
				constraint_diff[i] *= -1;
			}
		}

		if (constr_jco.rows() > 0)
		{	
			Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(constr_jco.transpose());
			int effective_rank = qr.rank();
			int full_rank = min(constr_jco.rows(), constr_jco.cols());

			if (effective_rank < full_rank)
			{
				frec << "   WARNING: constraint_jco is not full rank. Identifying independent constraints.";
				performance_log->log_event("   WARNING: constraint_jco is not full rank. Identifying independent constraints.");

				Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm = qr.colsPermutation();
				Eigen::VectorXi perm_indices = perm.indices();

				vector<string> reduced_Cnames;
				Eigen::MatrixXd reduced_constr_jco(effective_rank, constr_jco.cols());
				Eigen::VectorXd reduced_constraint_diff(effective_rank);

				for (int i = 0; i < effective_rank; i++)
				{
					int orig_idx = perm_indices(i);
					reduced_Cnames.push_back(Cnames[orig_idx]);
					reduced_constr_jco.row(i) = constr_jco.row(orig_idx);
					reduced_constraint_diff(i) = constraint_diff(orig_idx);
				}

				constr_jco = reduced_constr_jco;
				if (_constraint_jco != nullptr)
				{
					*_constraint_jco = constr_jco;
				}
				constraint_diff = reduced_constraint_diff;
				Cnames = reduced_Cnames;

				if (_cnames != nullptr)
					*_cnames = reduced_Cnames;
				

				ss.str("");
				ss << endl << "   reduced constraint matrix from rank " << constr_jco.rows() + (Cnames.size() - effective_rank)
					<< " to rank " << effective_rank << ". Kept constraints: ";
				for (size_t i = 0; i < reduced_Cnames.size(); i++)
				{
					if (i > 0) ss << ", ";
					ss << reduced_Cnames[i];
				}
				ss << endl;
				frec << ss.str();
				performance_log->log_event(ss.str());
			}
		}
		
		if ((constraint_diff.array().abs() > filter.get_viol_tol()).any()) 
		{
			ss.str("");
			ss << "   not on constraint but working set not empty" << endl;
			ss << "   current constraint diffs: " << constraint_diff.transpose() << endl;
			performance_log->log_event(ss.str());
			frec << ss.str();
		}
		const Eigen::MatrixXd& G = *hessian.e_ptr();

		string sqp_solve_method = pest_scenario.get_pestpp_options().get_sqp_solve_method(); 
		if (sqp_solve_method == "NULL" || sqp_solve_method == "NULL_SPACE")
		{
			x = _kkt_null_space(G, constr_jco, constraint_diff, grad_vector, _cnames);
			search_d = x.first;
			lm = x.second;
		}
		else if (sqp_solve_method == "DIRECT")
		{
			x = _kkt_direct(G, constr_jco, constraint_diff, grad_vector, _cnames);
			search_d = x.first;
			lm = x.second;
		}
		else // if "schur", "cg", ...
		{
			throw_sqp_error("sqp_solve_method not implemented");
		}
	}
	else  // solve unconstrained QP subproblem
	{
		Eigen::MatrixXd H_reg = regularize_hessian(*hessian.e_ptr(), "unconstrained");
		Eigen::LDLT<Eigen::MatrixXd> ldlt_H(H_reg);

		if (ldlt_H.info() != Eigen::Success || !ldlt_H.isPositive()) {
			message(1, "WARNING: LDLT failed for unconstrained case, using steepest descent");
			search_d = -grad_vector;
		}
		else {
			search_d = ldlt_H.solve(-grad_vector);
			double dir_dot_grad = search_d.dot(grad_vector);
			if (dir_dot_grad > 0) {
				message(2, "WARNING: search direction not descent, using steepest descent");
				search_d = -grad_vector;
			}
		}

		lm = Eigen::VectorXd::Zero(0);
	}
	return pair<Eigen::VectorXd, Eigen::VectorXd> (search_d, lm);
}

bool SeqQuadProgram::recalc_search_direction_vector(const string& rname, Parameters& dv_vals, Observations& obs_vals, Eigen::VectorXd& grad)
{
	stringstream ss;
	ofstream& frec = file_manager.rec_ofstream();

	if ((lm_en[rname].array() < 0).any())
	{
		constraint_mat_en[rname] = get_constraint_mat(dv_vals, obs_vals, working_set_tol, &lm_en[rname], cnames_en[rname]);
		if (constraint_mat_en[rname].first.get_row_names() != cnames_en[rname])
		{
			vector<string> prev_cnames = cnames_en[rname];
			cnames_en[rname] = constraint_mat_en[rname].first.get_row_names();
			set<string> snames(cnames_en[rname].begin(),cnames_en[rname].end());
			vector<string> dropped_cnames;
			for (auto c : prev_cnames)
			{
				if (snames.find(c) == snames.end())
					dropped_cnames.push_back(c);
			}

			ss.str("");
			ss << "   recalculating active set for realization " << rname << endl;
			frec << ss.str();
			performance_log->log_event(ss.str());

			Eigen::VectorXd prior_sd = search_d_en[rname];
			double prior_search_d_norm = search_d_en[rname].norm();

			constraint_jco_en[rname] = constraint_mat_en[rname].first.e_ptr()->toDense();

			Covariance backup_hessian = hessian;
			hessian = hessian_en[rname];
			used_hessian = Covariance();

			pair<Eigen::VectorXd, Eigen::VectorXd> x = calc_search_direction_vector(dv_vals, obs_vals, grad, &constraint_jco_en[rname], &cnames_en[rname]);
			search_d_en[rname] = x.first;
			lm_en[rname] = x.second;

			if (!used_hessian.get_col_names().empty())
				hessian_en[rname] = used_hessian;

			hessian = backup_hessian;

			Eigen::VectorXd unscaled_search_d = search_d_en[rname];
			if (pest_scenario.get_pestpp_options().get_sqp_rescale_search_dir())
			{
				ParameterInfo par_info = pest_scenario.get_ctl_parameter_info();
				Parameters lbnd = par_info.get_low_bnd(dv_names);
				Parameters ubnd = par_info.get_up_bnd(dv_names);
				ParamTransformSeq par_transform = pest_scenario.get_base_par_tran_seq();
				par_transform.ctl2numeric_ip(lbnd);
				par_transform.ctl2numeric_ip(ubnd);

				double rangesq = 0.0;
				for (int i = 0; i < dv_names.size(); i++)
				{
					double lb = lbnd[dv_names[i]];
					double ub = ubnd[dv_names[i]];
					rangesq += pow(ub - lb, 2);
				}
				rangesq = pow(rangesq, 0.5);
				if (search_d_en[rname].norm() > rangesq)
				{
					search_d_en[rname] = rangesq * search_d_en[rname] / search_d_en[rname].norm();
				}
			}

			Eigen::VectorXd delta_sd = search_d_en[rname] - prior_sd;

			ss.str("");
			if (dropped_cnames.size() == 0)
				ss << "   constraint dropped: NONE" << endl;
			else
			{
				ss << "   constraint dropped:";
				for (auto d : dropped_cnames)
					ss << "  " << d;
				ss << endl;
				if (cnames_en[rname].size() == 0)
				{
					ss << "   new working set: EMPTY" << endl;
					ss << "   new unscaled step length: " << unscaled_search_d.norm() << endl;
					ss << "   new scaled step length: " << search_d_en[rname].norm() << endl;
					ss << "   new search direction: " << search_d_en[rname].transpose() << endl;
					
				}
				else
				{
					ss << "   new working set:" << endl;
					int i = 0;
					for (auto& r : cnames_en[rname])
					{
						ss << "         " << r << " (lm = " << lm_en[rname][i] << ")" << endl;
						i++;
					}
					ss << "   new unscaled step length: " << unscaled_search_d.norm() << endl;
					ss << "   new scaled step length: " << search_d_en[rname].norm() << endl;
					ss << "   new search direction: " << search_d_en[rname].transpose() << endl;
				}

			}
			ss << "   summary of search direction change: " << endl;
			ss << "      max abs decvar step length change: " << delta_sd.cwiseAbs().maxCoeff() << endl;
			ss << "      min abs decvar step length change: " << delta_sd.cwiseAbs().minCoeff() << endl;
			ss << "      change in unscaled step length: " << search_d_en[rname].norm() - prior_search_d_norm << endl;

			frec << ss.str() << endl;

			return (lm_en[rname].array() < 0).any();
		}
		else
		{
			throw_sqp_error("negative multipliers exist but working set didnt change for " + rname);
			return false;
		}
	}
	else
		return false;
}

bool SeqQuadProgram::solve_new_ensemble()
{
	stringstream ss;
	ofstream& frec = file_manager.rec_ofstream();

	if (dv.shape().first <= error_min_reals)
	{
		message(0, "too few active realizations:", oe.shape().first);
		message(1, "need more than ", error_min_reals);
		throw_sqp_error(string("too few active realizations, cannot continue"));
	}
	else if (dv.shape().first < warn_min_reals)
	{
		ss.str("");
		ss << "WARNING: less than " << warn_min_reals << " active realizations...might not be enough";
		string s = ss.str();
		message(1, s);
	}

	prev_ctl_dv_values = current_ctl_dv_values;

	message(0, "starting solve for iteration:", iter);
	cout << "...calculating search direction for each realization (see .rec file for more details)" << endl;

	performance_log->log_event("reordering variables in dv");
	dv.reorder(vector<string>(), dv_names);
	dv.transform_ip(ParameterEnsemble::transStatus::NUM);
	ParameterEnsemble _dvs = dv;
	_dvs.drop_rows(vector<string>{BASE_REAL_NAME}, true);

	int local_subset_size = pest_scenario.get_pestpp_options().get_sqp_subset_size();
	if (local_subset_size < 0)
	{
		ss.str("");

		local_subset_size = (int)((double)_dvs.shape().first) * ((-1. * (double)local_subset_size) / 100.);

		ss << "subset defined as a percentage of ensemble size, using " << local_subset_size;
		ss << " realizations for subset" << endl;
		message(2, ss.str());
		if (local_subset_size < 4)
		{
			ss.str("");
			ss << "percentage-based subset size too small, increasing to 4" << endl;
			local_subset_size = 4;
			message(2, ss.str());
		}
	}
	if ((use_subset) && (local_subset_size > _dvs.shape().first))
	{
		ss.str("");
		ss << "subset size (" << local_subset_size << ") greater than ensemble size (" << _dvs.shape().first << ")";
		frec << "  ---  " << ss.str() << endl;
		cout << "  ---  " << ss.str() << endl;
		frec << "  ...reducing subset size to " << _dvs.shape().first << endl;
		cout << "  ...reducing subset size to " << _dvs.shape().first << endl;
		local_subset_size = _dvs.shape().first;
	}
	else if (pest_scenario.get_pestpp_options().get_sqp_alpha_mults().size() == 1)
	{
		ss.str("");
		ss << "only testing one scale factor, not using subset";
		frec << "  ---  " << ss.str() << endl;
		cout << "  ---  " << ss.str() << endl;
		local_subset_size = _dvs.shape().first;
	}
		
	if (!sampling_tracking_initialized)
	{
		unselected_dv_indices.clear();
		selected_dv_indices.clear();
		for (int i = 0; i < _dvs.shape().first; i++)
			unselected_dv_indices.insert(i);
		
		sampling_tracking_initialized = true;
	}

	if (unselected_dv_indices.empty() || local_subset_size > unselected_dv_indices.size())
	{
		unselected_dv_indices.clear();
		selected_dv_indices.clear();
		for (int i = 0; i < _dvs.shape().first; i++)
			unselected_dv_indices.insert(i);
	}

	vector<int> subset_idxs = get_subset_idxs(_dvs.shape().first, local_subset_size);
	vector<string> subset_real_names;
	subset_real_names.reserve(subset_idxs.size());
	for (int idx : subset_idxs)
	{
		unselected_dv_indices.erase(idx);
		selected_dv_indices.insert(idx);
		subset_real_names.push_back(_dvs.get_real_names()[idx]);
	}

	ParameterEnsemble _drawn_dvs = _dvs;
	_drawn_dvs.keep_rows(subset_real_names, true);

	Parameters dv_vals = current_ctl_dv_values;
	Observations obs_vals = current_obs;

	Eigen::VectorXd grad = current_grad_vector.get_data_eigen_vec(dv_names);
	vector<string> drawn_real_names = _drawn_dvs.get_real_names();
	drawn_real_names.push_back(BASE_REAL_NAME);

	double rangesq = -1.0;
	if (pest_scenario.get_pestpp_options().get_sqp_rescale_search_dir())
	{
		ParameterInfo par_info = pest_scenario.get_ctl_parameter_info();
		Parameters lbnd = par_info.get_low_bnd(dv_names);
		Parameters ubnd = par_info.get_up_bnd(dv_names);
		ParamTransformSeq par_transform = pest_scenario.get_base_par_tran_seq();
		par_transform.ctl2numeric_ip(lbnd);
		par_transform.ctl2numeric_ip(ubnd);

		rangesq = 0.0;
		for (int i = 0; i < dv_names.size(); i++)
		{
			double lb = lbnd[dv_names[i]];
			double ub = ubnd[dv_names[i]];
			rangesq += pow(ub - lb, 2);
		}
		rangesq = pow(rangesq, 0.5);
	}

	hessian_en.clear();
	for (auto d : drawn_real_names)
	{
		ss.str("");
		ss << "...calculating search direction for realization " << d << endl;
		frec << ss.str();
		performance_log->log_event(ss.str());

		Eigen::VectorXd real_dv_vec = dv.get_real_vector(d);
		dv_vals.update_without_clear(dv_names, real_dv_vec);
		Eigen::VectorXd real_obs_vec = oe.get_real_vector(d);
		obs_vals.update_without_clear(oe.get_var_names(), real_obs_vec);
		
		constraint_mat_en[d] = get_constraint_mat(dv_vals, obs_vals, working_set_tol);
		cnames_en[d] = constraint_mat_en[d].first.get_row_names();
		constraint_jco_en[d] = constraint_mat_en[d].first.e_ptr()->toDense();
		current_obj_en[d] = get_obj_value(dv_vals, obs_vals);

		if (hessian_en.find(d) == hessian_en.end())
			hessian_en[d] = hessian;
		
		Covariance backup_hessian = hessian;
		hessian = hessian_en[d];
		used_hessian = Covariance();

		pair<Eigen::VectorXd, Eigen::VectorXd> x = calc_search_direction_vector(dv_vals, obs_vals, grad, &constraint_jco_en[d], &cnames_en[d]);
		search_d_en[d] = x.first;
		lm_en[d] = x.second;

		if (!used_hessian.get_col_names().empty())
			hessian_en[d] = used_hessian;
		
		hessian = backup_hessian;

		Eigen::VectorXd unscaled_search_d = search_d_en[d];
		double dir_norm = search_d_en[d].norm();
		if (rangesq > 0.0 && dir_norm > rangesq)
		{
			search_d_en[d] = rangesq * search_d_en[d] / dir_norm;
		}

		ss.str("");
		if (cnames_en[d].size() == 0)
		{
			ss << "   working set: EMPTY" << endl;
			ss << "   unscaled step length: " << unscaled_search_d.norm() << endl;
			ss << "   scaled step length: " << search_d_en[d].norm() << endl;
			ss << "   search direction: " << search_d_en[d].transpose() << endl;
		}
		else
		{
			ss << "   working set:" << endl;
			int i = 0;
			for (auto c : cnames_en[d])
			{
				ss << "         " << c << " (lm = " << lm_en[d][i] << ")" << endl;
				i++;
			}
			ss << "   unscaled step length: " << unscaled_search_d.norm() << endl;
			ss << "   scaled step length: " << search_d_en[d].norm() << endl;
			ss << "   search direction: " << search_d_en[d].transpose() << endl;
		}
		frec << ss.str() << endl;

		while (true)
		{
			bool changed = recalc_search_direction_vector(d, dv_vals, obs_vals, grad);
			if (!changed)
				break;
		}
	}

	FilterRec search = run_search_routine(grad, &_drawn_dvs);

	//needed for bfgs hessian update
	selected_ls_child = search.real_name;
	selected_ls_parent = ls_parent_map[selected_ls_child];
	step_k = step_length_map[search.real_name];
	if (step_k.size() != 0 && search.iter == iter)
	{
		is_good_search = true;
		Eigen::VectorXd parent_dv_vec = dv.get_real_vector(selected_ls_parent);
		prev_ctl_dv_values.update_without_clear(dv_names, parent_dv_vec);
		constraint_jco = constraint_jco_en[selected_ls_parent];
		current_constraint_mat = constraint_mat_en[selected_ls_parent].first;
		cnames = cnames_en[selected_ls_parent];
		prev_constraint_mat = current_constraint_mat;
		lambda = lm_en[selected_ls_parent];
		if (hessian_en.find(selected_ls_parent) != hessian_en.end())
		{
			hessian = hessian_en[selected_ls_parent];
		}
		BASE_SCALE_FACTOR = 1.0;
	}
	else
	{
		is_good_search = false;
		BASE_SCALE_FACTOR = max(1E-4, BASE_SCALE_FACTOR * SF_DEC_FAC);
	}
	message(1, "new base scale factor: ", BASE_SCALE_FACTOR);

	if (search.viol_val == 0.0)
	{
		is_base_infeas = false;
		cma_reset_archive = true;
	}
	else
	{
		is_base_infeas = true;
		cma_reset_archive = false;
	}

	return (search.viol_val == 0.0);
}

FilterRec SeqQuadProgram::run_search_routine(Eigen::VectorXd& grad, ParameterEnsemble* drawn_dvs, bool recalc)
{
	if (recalc)
		recalc_attempt++;
	else
		recalc_attempt = 0;

	string search_method = pest_scenario.get_pestpp_options().get_sqp_search_method();
	if (search_method == "LINE" || search_method == "LINE_SEARCH" || search_method == "LS")
		return line_search(search_d_en, grad, current_obj_en, drawn_dvs, recalc);
	else if (search_method == "TRUST_REGION" || search_method == "TRUST" || search_method == "TR")
		return trust_region_step(grad, current_obj_en, cnames_en, constraint_jco_en, drawn_dvs, recalc);
	else
		throw_sqp_error("search_method not recognized");
	return FilterRec();
}

/**
 * @brief Seek feasible.
 *
 * @return Description.
 */
bool SeqQuadProgram::seek_feasible()
{
	stringstream ss;
	message(1, "seeking feasibility with iterative ensemble smoother solution");
	Pest ies_pest_scenario;
	string pst_filename = pest_scenario.get_pst_filename();
	ifstream fin(pest_scenario.get_pst_filename());
	ies_pest_scenario.process_ctl_file(fin, pst_filename);
	set<string>snames(dv_names.begin(), dv_names.end());
	set<string>::iterator send = snames.end();
	ParameterInfo* pi = ies_pest_scenario.get_ctl_parameter_info_ptr_4_mod();
	ParamTransformSeq pts = ies_pest_scenario.get_base_par_tran_seq_4_mod();
	TranFixed* tf_ptr = pts.get_fixed_ptr_4_mod();

	Parameters& ctl_pars = ies_pest_scenario.get_ctl_parameters_4_mod();

	for (auto& name : ies_pest_scenario.get_ctl_ordered_par_names())
	{
		if (snames.find(name) == send)
		{
			if (pi->get_parameter_rec_ptr_4_mod(name)->tranform_type != ParameterRec::TRAN_TYPE::FIXED)
			{
				pi->get_parameter_rec_ptr_4_mod(name)->tranform_type = ParameterRec::TRAN_TYPE::FIXED;
				tf_ptr->insert(name, ctl_pars.get_rec(name));
			}
		}
		else
			ctl_pars.update_rec(name,current_ctl_dv_values.get_rec(name));

	}
	snames.clear();
	vector<string> names = constraints.get_obs_constraint_names();
	if (names.size() == 0)
		throw_sqp_error("SQP::seek_feasible() error: no obs-based constraints found");
	snames.insert(names.begin(), names.end());
	send = snames.end();
	if (snames.find(obj_obs) != send)
	{
		snames.erase(obj_obs);
	}
	ObservationInfo* oi = ies_pest_scenario.get_observation_info_ptr();

	Observations shifted = pest_scenario.get_ctl_observations();
	if (constraints.get_use_chance())
    {
	    shifted = constraints.get_chance_shifted_constraints(current_obs);
    }
	Observations& ctl_obs = ies_pest_scenario.get_ctl_observations_4_mod();
	map<string,double> viol_map = constraints.get_unsatified_obs_constraints(current_obs,filter.get_viol_tol(), true);
	for (auto& name : ies_pest_scenario.get_ctl_ordered_obs_names())
	{
		if (snames.find(name) == send)
		{
			oi->get_observation_rec_ptr_4_mod(name)->weight = 0.0;

		}
		else
		{
		    if (viol_map.find(name) != viol_map.end()) {
                ctl_obs.update_rec(name, shifted.get_rec(name));
                oi->get_observation_rec_ptr_4_mod(name)->group = "__eqconstraint__" + name;
            }
		}
	}

	snames = ies_pest_scenario.get_pestpp_options().get_passed_args();
	if (snames.find("IES_BAD_PHI_SIGMA") == snames.end())
    {
	    ies_pest_scenario.get_pestpp_options_ptr()->set_ies_bad_phi_sigma(1.25);
    }

    if (snames.find("IES_LAMBBDA_MULTS") == snames.end())
    {
        ies_pest_scenario.get_pestpp_options_ptr()->set_ies_lam_mults(vector<double>{0.1,1.0,10});
    }

    if (snames.find("LAMBBDA_SCALE_FAC") == snames.end())
    {
        ies_pest_scenario.get_pestpp_options_ptr()->set_lambda_scale_vec(vector<double>{0.5,1.0});
    }
    if (snames.find("IES_NUM_REALS") == snames.end()) {
        ies_pest_scenario.get_pestpp_options_ptr()->set_ies_num_reals(
                max(max(pest_scenario.get_pestpp_options().get_sqp_num_reals(), (int)(constraints.num_constraints()*1.1)),30));
    }
    if (snames.find("IES_SUBSET_SIZE") == snames.end()) {
        ies_pest_scenario.get_pestpp_options_ptr()->set_ies_subset_size(-5);
    }
    ies_pest_scenario.get_pestpp_options_ptr()->set_ies_no_noise(true);
	ies_pest_scenario.get_pestpp_options_ptr()->set_ies_obs_csv("");
    ies_pest_scenario.get_pestpp_options_ptr()->set_ies_obs_restart_csv("");
    ies_pest_scenario.get_pestpp_options_ptr()->set_ies_par_csv("");
    ies_pest_scenario.get_pestpp_options_ptr()->set_ies_par_restart_csv("");
    ies_pest_scenario.get_control_info_4_mod().noptmax = pest_scenario.get_pestpp_options().get_sqp_seek_feas_max_iter();
    ss.str("");
    string org_base = file_manager.get_base_filename();
	string safe_base = org_base;
	size_t pos = safe_base.find_last_of("/\\");
	if (pos != string::npos)
		safe_base = safe_base.substr(pos + 1);
	for (char& ch : safe_base)
	{
		if (ch == '/' || ch == '\\')
			ch = '_';
	}
	ss << "feas_ies_" << iter << "_" << safe_base;

    file_manager.set_base_filename(ss.str());
    IterEnsembleSmoother ies(ies_pest_scenario, file_manager, output_file_writer, performance_log, run_mgr_ptr);
    if (use_ensemble_grad) {

        ies.set_pe(dv);
        ies.set_oe(oe);
        ies.set_noise_oe(oe_base);
        ies.initialize(iter,true,true);
    }
    else{
        ies.initialize();
    }

	ies.iterate_2_solution();
    file_manager.set_base_filename(org_base);
	
	ParameterEnsemble* ies_pe_ptr = ies.get_pe_ptr();
	ObservationEnsemble* ies_oe_ptr = ies.get_oe_ptr();
	vector<string> oreal_names = ies_oe_ptr->get_real_names();
	map<string,double> aphi_map = ies.get_phi_handler().get_phi_map(L2PhiHandler::phiType::ACTUAL);

	ies_pe_ptr->transform_ip(ParameterEnsemble::transStatus::CTL);
	names = ies_pe_ptr->get_var_names();

	Eigen::VectorXd cdv = current_ctl_dv_values.get_data_eigen_vec(dv_names);
	double mndiff = 1.0e+300;
	int mndiff_idx = -1;
	for (int i = 0; i < ies_pe_ptr->shape().first; i++)
	{
		if (aphi_map[oreal_names[i]] < mndiff)
		{
			mndiff = aphi_map[oreal_names[i]];
			mndiff_idx = i;
		}
	}
	ss.str("");
	ss << "updating current decision variable values with realization " << ies_pe_ptr->get_real_names()[mndiff_idx];
	ss << ", with minimum weighted constraint phi of " << mndiff;
	message(1,ss.str());
	cdv = ies_pe_ptr->get_real_vector(mndiff_idx);
	current_ctl_dv_values.update_without_clear(names, cdv);

	cdv = ies.get_oe().get_real_vector(mndiff_idx);
	names = ies.get_oe().get_var_names();
	current_obs.update(names, cdv);
	constraints.sqp_report(iter, current_ctl_dv_values, current_obs, true, "post feasible seek");

	last_best = get_obj_value(current_ctl_dv_values, current_obs);
	last_viol = constraints.get_sum_of_violations(current_ctl_dv_values, current_obs);
	best_phis[best_phis.size()-1] = last_best;
	best_violations[best_violations.size() -1] = last_viol;
	message(1, "finished seeking feasible, reset best phi,infeasible value to ", vector<double>{last_best,last_viol});

	message(1, "updating CMA archives with IES results for seeking feasibility");
	map<string, double> obj_map = get_obj_map(*ies_pe_ptr, *ies_oe_ptr);
	map<string, double> total_viol_map;
	map<string, map<string, double>> violations_nominal = constraints.get_ensemble_violations_map(*ies_pe_ptr, *ies_oe_ptr, 0.0, true);
	for (auto& real_name : oreal_names)
	{
		double infeas_sum_nom = 0.0;
		for (auto& v : violations_nominal[real_name])
			infeas_sum_nom += v.second;
		total_viol_map[real_name] = infeas_sum_nom;

		double obj_val = obj_map[real_name];
		double viol_val = total_viol_map[real_name];

		Parameters p = pest_scenario.get_ctl_parameters();
		Eigen::VectorXd pv = ies_pe_ptr->get_real_vector(real_name);
		p.update_without_clear(dv_names, pv);

		Observations o = pest_scenario.get_ctl_observations();
		Eigen::VectorXd ov = ies_oe_ptr->get_real_vector(real_name);
		o.update_without_clear(ies_oe_ptr->get_var_names(), ov);

		double vpad = viol_val;
		if (viol_val < pest_scenario.get_pestpp_options().get_sqp_viol_pad())
			vpad = 0.0;
		filter.update(obj_val, viol_val, vpad, p, o, real_name, -iter);  
	}

	if ((iter > 0) && (use_cma))
	{
		if (last_viol < filter.get_viol_tol())
			cma.update_archives(*ies_pe_ptr, obj_map, total_viol_map, to_string(-iter), true);
		else
			cma.update_archives(*ies_pe_ptr, obj_map, total_viol_map, to_string(-iter), true);
	}

	return false;
}

Eigen::VectorXd SeqQuadProgram::solve_constrained_trust_region_step(const Eigen::MatrixXd& B, const Eigen::VectorXd& g, const Eigen::MatrixXd& A, double radius)
{
	if (A.rows() == 0)
		return solve_trust_region_subproblem_dogleg(B, g, radius);

	Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(A);
	Eigen::MatrixXd Z = cod.matrixZ();

	if (Z.cols() == 0)
		return Eigen::VectorXd::Zero(B.rows());

	Eigen::MatrixXd B_reduced = Z.transpose() * B * Z;
	Eigen::VectorXd g_reduced = Z.transpose() * g;

	Eigen::VectorXd s = solve_trust_region_subproblem_dogleg(B_reduced, g_reduced, radius);
	return Z * s;
}

/**
 * @brief Get obj value.
 *
 * @param _current_ctl_dv_vals Description.
 * @param _current_obs Description.
 *
 * @return Description.
 */
double SeqQuadProgram::get_obj_value(Parameters& _current_ctl_dv_vals, Observations& _current_obs)
{
	double v = 0;
	if (use_obj_obs)
	{
		v =  _current_obs.get_rec(obj_func_str);
	}
	else
	{
		if (use_obj_pi)
		{
			PriorInformationRec pi = pest_scenario.get_prior_info_ptr()->get_pi_rec(obj_func_str);
			v = pi.calc_sim_and_resid(_current_ctl_dv_vals).first;

		}
		else
		{
			Parameters pars = _current_ctl_dv_vals;
			ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();
			pts.ctl2numeric_ip(pars);
			for (auto& dv_name : dv_names)
				v += obj_func_coef_map[dv_name] * pars.get_rec(dv_name);
		}
	}
	return v;
}

/**
 * @brief Get obj map.
 *
 * @param _dv Description.
 * @param _oe Description.
 *
 * @return Description.
 */
map<string, double> SeqQuadProgram::get_obj_map(ParameterEnsemble& _dv, ObservationEnsemble& _oe)
{
	Eigen::VectorXd obj_vec = get_obj_vector(_dv, _oe);
	vector<string> real_names = _dv.get_real_names();
	map<string, double> obj_map;
	for (int i = 0; i < real_names.size(); i++)
		obj_map[real_names[i]] = obj_vec[i];

	return obj_map;


}

/**
 * @brief Get obj vector.
 *
 * @param _dv Description.
 * @param _oe Description.
 *
 * @return Description.
 */
Eigen::VectorXd SeqQuadProgram::get_obj_vector(ParameterEnsemble& _dv, ObservationEnsemble& _oe)
{
	Eigen::VectorXd obj_vec(_dv.shape().first);
	if (use_obj_obs)
	{
		obj_vec = _oe.get_var_vector(obj_func_str);
	}
	else
	{
		_dv.transform_ip(ParameterEnsemble::transStatus::NUM);
		Parameters pars = pest_scenario.get_ctl_parameters();
		ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();
		pts.ctl2numeric_ip(pars);
		Eigen::VectorXd real;
		vector<string> vnames = _dv.get_var_names();
		double v;
		for (int i = 0; i < _dv.shape().first; i++)
		{
			real = _dv.get_real_vector(i);
			pars.update_without_clear(vnames, real);
			pts.numeric2ctl_ip(pars);
			v = get_obj_value(pars, current_obs);
			obj_vec[i] = v;
			pts.ctl2numeric_ip(pars);
		}
	}
	return obj_vec;
}

tuple<FilterRec, SqpFilter> SeqQuadProgram::pick_from_filter(ParameterEnsemble& dv_candidates, ObservationEnsemble& _oe, bool recalc)
{
	stringstream ss;
	ofstream& frec = file_manager.rec_ofstream();
	Eigen::VectorXd obj_vec = get_obj_vector(dv_candidates, _oe);
	double oext, oviol = 0.0, nviol = 0.0;

	vector<string> real_names = dv_candidates.get_real_names();

	map<string, Parameters> par_map;
	map<string, Observations> obs_map;
	for (auto& d : real_names)
	{
		Parameters p = current_ctl_dv_values;
		Eigen::VectorXd t = dv_candidates.get_real_vector(d);
		p.update_without_clear(dv_names, t);
		par_map[d] = p;

		Observations o = current_obs;
		t = _oe.get_real_vector(d);
		o.update_without_clear(_oe.get_var_names(), t);
		obs_map[d] = o;
	}
	double viol_pad = pest_scenario.get_pestpp_options().get_sqp_viol_pad();
	map<string, map<string, double>> violations, violations_nominal;
	if (constraints.get_use_chance() && (sqp_risk != 0.5))
	{
		violations = constraints.get_ensemble_violations_map(dv_candidates, _oe, viol_pad, true, &_oe, sqp_risk);
		violations_nominal = constraints.get_ensemble_violations_map(dv_candidates, _oe, 0.0, true, &_oe, sqp_risk);
	}
	else
	{
		violations = constraints.get_ensemble_violations_map(dv_candidates, _oe, viol_pad, true);
		violations_nominal = constraints.get_ensemble_violations_map(dv_candidates, _oe, 0.0, true);
	}

	vector<string> onames = _oe.get_var_names();
	vector<string> vnames = dv_candidates.get_var_names();
	obj_map = get_obj_map(dv_candidates, _oe);
	
	vector<double> infeas_vec, nviol_vec, feas_dist_map;
	vector<int> accept_idxs, feas_idxs;

	message(0, "current best phi:", last_best);
	ss.str("");
	ss << "evaluating " << obj_vec.size() << " candidate realizations and updating filter";
	message(1, ss.str());

	ss.str("");
	ss << string(80, '-') << endl;
	ss << left << setw(25) << "real_name"
		<< setw(15) << "phi"
		<< setw(20) << "infeasibility"
		<< setw(15) << "accept/reject" << endl;
	ss << string(80, '-') << endl;

	SqpFilter candidate_filter;
	if (!recalc)
		candidate_filter = filter;


	for (int i = 0; i < obj_vec.size(); i++)
	{
		ss << left << setw(25) << real_names[i]
			<< setw(15) << obj_vec[i];

		double infeas_sum = 0.0, infeas_sum_nom = 0.0, feas_dist = 0;
		for (auto& v : violations[real_names[i]])
			infeas_sum += v.second;
		infeas_vec.push_back(infeas_sum);
		for (auto& v : violations_nominal[real_names[i]])
			infeas_sum_nom += v.second;
		nviol_vec.push_back(infeas_sum_nom);
		total_viol_map[real_names[i]] = infeas_sum_nom;
		bool filter_accept = candidate_filter.accept(obj_vec[i], infeas_sum_nom, infeas_sum, par_map[real_names[i]], obs_map[real_names[i]], real_names[i], iter, true);
		string accept_reject = filter_accept ? "accept" : "reject";

		ss << setw(20) << infeas_sum_nom
			<< setw(15) << accept_reject
			<< endl;

		if (filter_accept)
			accept_idxs.push_back(i);
	}

	ss << string(80, '-') << endl;
	frec << ss.str();
	if (accept_idxs.size() > 0)
		message(1, "number of candidate realizations passing filter: ", accept_idxs.size());
	else
		message(1, "no realizations passed the filter");
	candidate_filter.report(frec, iter);

	FilterRec selected;
	vector<FilterRec> feasible = candidate_filter.get_feasible_solutions();

	if (obj_sense == "minimize")
		oext = 1E+30;
	else
		oext = -1E+30;

	if (feasible.size() > 0)
	{
		for (const auto& f : feasible)
		{
			if ((obj_sense == "minimize") && (f.obj_val < oext))
			{
				oext = f.obj_val;
				selected = f;
			}
			else if ((obj_sense == "maximize") && (f.obj_val > oext))
			{
				oext = f.obj_val;
				selected = f;
			}
		}
	}
	else
	{
		vector<FilterRec> filterset = candidate_filter.get_filter_members();
		oext = 1E+30;
		for (const auto f : filterset)
		{
			if (f.viol_val < oext)
			{
				oext = f.viol_val;
				selected = f;
			}
		}
	}

	return { selected, candidate_filter };
}

FilterRec SeqQuadProgram::pick_upgrade_and_update_current(ParameterEnsemble& dv_candidates, ObservationEnsemble& _oe, bool cma_reset_arc, bool report, ParameterEnsemble* dvs_subset, bool recalc)
{
	stringstream ss;
	ofstream& frec = file_manager.rec_ofstream();

	vector<string> onames = _oe.get_var_names();
	vector<string> vnames = dv_candidates.get_var_names();
	
	auto pick = pick_from_filter(dv_candidates, _oe, recalc);
	FilterRec selected = get<0>(pick);
	SqpFilter candidate_filter = get<1>(pick);

	if (recalc)
		return selected;

	if (selected.obj_val == last_best && selected.viol_val == last_viol)
	{

		if (!found_feasible)
		{
			message(1, "still seeking feasible solution...choosing next least infeasible filter member");

			vector<FilterRec> filterset = candidate_filter.get_filter_members();
			sort(filterset.begin(), filterset.end(),
				[](const FilterRec& a, const FilterRec& b)
				{
					return a.viol_val < b.viol_val;
				});

			for (const auto& fr : filterset)
			{
				if (fr.viol_val != last_viol)
				{
					selected = fr;
					break;
				}
			}
		}
		else
		{
			message(1, "still no better solution...allowing some violation");
			vector<FilterRec> filterset = candidate_filter.get_feasible_solutions(true);
			bool found = false;
			if (filterset.size() > 0)
			{
				if (obj_sense == "minimize")
				{
					sort(filterset.begin(), filterset.end(),
						[](const FilterRec& a, const FilterRec& b)
						{
							return a.obj_val < b.obj_val;
						});
				}
				else
				{
					sort(filterset.begin(), filterset.end(),
						[](const FilterRec& a, const FilterRec& b)
						{
							return a.obj_val > b.obj_val;
						});
				}

				auto it = filterset.begin();

				while (it != filterset.end())
				{
					bool better_obj = false;
					if (obj_sense == "minimize")
					{
						better_obj = (it->obj_val < last_best);
					}
					else
					{
						better_obj = (it->obj_val > last_best);
					}

					if (better_obj && it->obj_val != last_best && find(best_phis.begin(), best_phis.end(), it->obj_val) == best_phis.end())
					{
						selected = *it;
						found = true;
						break;
					}
					it++;
				}
			}

			if (!found)
			{
				filterset = candidate_filter.get_filter_members();

				sort(filterset.begin(), filterset.end(),
					[](const FilterRec& a, const FilterRec& b)
					{
						return a.viol_val < b.viol_val;
					});

				auto it = filterset.begin();

				while (it != filterset.end())
				{
					bool better_obj = false;
					if (obj_sense == "minimize")
					{
						better_obj = (it->obj_val < last_best);
					}
					else
					{
						better_obj = (it->obj_val > last_best);
					}

					if (better_obj && it->obj_val != last_best && find(best_phis.begin(), best_phis.end(), it->obj_val) == best_phis.end())
					{
						selected = *it;
						found = true;
						break;
					}
					it++;
				}
			}

			if (!found)
			{
				message(1, "WARNING: no unselected candidate found with better objective than last_best");
			}

			if (selected.obj_val == last_best && selected.viol_val == last_viol && candidate_filter.get_filter_members().size() > 0)
			{
				message(1, "resorting to global filter by merit");
				selected = pick_from_filter_by_merit(candidate_filter);
				reset = true;
			}
		}
	}
	
	if (selected.viol_val > 0.0 || find(best_phis.begin(), best_phis.end(), selected.obj_val) == best_phis.end())
	{
		if (pest_scenario.get_pestpp_options().get_sqp_cma_parent_num() == 0)
		{
			int curr_parent_num = cma.get_parent_num();
			int ratio = pest_scenario.get_pestpp_options().get_sqp_num_reals() / curr_parent_num + 1;
			int new_parent_num = max(5, pest_scenario.get_pestpp_options().get_sqp_num_reals() / ratio);
			cma.set_parent_num(new_parent_num);
			
		}
		message(1, "updating CMA archive of size: ", cma.get_parent_num());
		cma.update_archives(dv_candidates, obj_map, total_viol_map, to_string(iter), true);
	}
	else
	{
		
		if (pest_scenario.get_pestpp_options().get_sqp_cma_parent_num() == 0)
		{
			cma.set_parent_num(pest_scenario.get_pestpp_options().get_sqp_num_reals() / 4);
		}
		message(1, "updating CMA archive of size: ", cma.get_parent_num());
		cma.update_archives(dv_candidates, obj_map, total_viol_map, to_string(iter), true);
	}

	bool is_violated = (selected.viol_val >= 1E-10);
	bool is_recycled = (find(best_phis.begin(), best_phis.end(), selected.obj_val) != best_phis.end());

	filter.update(selected.obj_val, selected.viol_val, selected.viol_padded, selected.dp_val, selected.oe_val, selected.real_name, iter);
	filter.report(frec, iter);
	current_ctl_dv_values = selected.dp_val;
	last_best = selected.obj_val;
	last_viol = selected.viol_val;
	if (last_viol < 1E-12)
		found_feasible = true;

	best_phis.push_back(last_best);
	best_violations.push_back(last_viol);

	ss.str("");
	ss << "updating BASE with realization: " << selected.iter << "|" << selected.real_name << endl;
	ss << "   phi: " << last_best << ", violation: " << last_viol;
	message(1, ss.str());

	cnames_base = cnames_en[selected.real_name];
	lm_base = lm_en[selected.real_name];
	Parameters p;
	Observations o;

	Eigen::VectorXd t = selected.dp_val.get_data_eigen_vec(vnames);
	p.update_without_clear(vnames, t);

	t = selected.oe_val.get_data_eigen_vec(onames);
	o.update_without_clear(onames, t);
	constraints.sqp_report(iter, p, o, true, selected.real_name);

	for (auto& d : vnames)
		current_ctl_dv_values[d] = p.get_rec(d);
	current_obs.update_without_clear(_oe.get_var_names(), o.get_data_eigen_vec(_oe.get_var_names()));

	return selected;
}

FilterRec SeqQuadProgram::pick_from_filter_by_merit(SqpFilter _filtered)
{
	vector<FilterRec> filterset = _filtered.get_filter_members();

	double vmin = numeric_limits<double>::infinity(), vmax = -numeric_limits<double>::infinity();
	for (const auto& fr : filterset) {
		vmin = min(vmin, fr.viol_val);
		vmax = max(vmax, fr.viol_val);
	}
	double mu = (1 + vmax) / (1 + vmin);
	double vref = max(vmax - vmin, 1e-30);

	auto merit = [&](const FilterRec& fr) {
		const double f_signed = (obj_sense == "minimize") ? fr.obj_val : -fr.obj_val;
		return f_signed + mu * (fr.viol_val / vref);
		};

	const FilterRec* best = &filterset[0];
	double best_m = merit(filterset[0]);
	for (size_t i = 1; i < filterset.size(); ++i)
	{
		double m = merit(filterset[i]);
		if (m < best_m)
		{
			best_m = m;
			best = &filterset[i];
		}
	}
	return *best;
}

void SeqQuadProgram::report_and_save_ensemble(ParameterEnsemble& _dv, ObservationEnsemble& _oe)
{
	map<string, map<string, double>> violations_nominal = constraints.get_ensemble_violations_map(_dv, _oe, 0.0, true);
	vector<string> real_names = _dv.get_real_names();
	double viol_tol = filter.get_viol_tol();
	int active_feas = 0;

	for (const auto& real_name : real_names)
	{
		double infeas_sum = 0.0;
		for (const auto& v : violations_nominal[real_name])
			infeas_sum += v.second;

		if (infeas_sum < viol_tol)
			active_feas++;
	}

	ofstream& frec = file_manager.rec_ofstream();
	frec << endl << "  ---  SeqQuadProgram iteration " << iter << " report  ---  " << endl;
	frec << "   number of active realizations:  " << _dv.shape().first << endl;
	frec << "   number of model runs:           " << run_mgr_ptr->get_total_runs() << endl;
	frec << "   number of active feasible solutions (% of active):   " << active_feas << "(" << 100 * active_feas /_dv.shape().first << "%)" << endl;

	cout << endl << "  ---  SeqQuadProgram iteration " << iter << " report  ---  " << endl;
	cout << "   number of active realizations:   " << _dv.shape().first << endl;
	cout << "   number of model runs:            " << run_mgr_ptr->get_total_runs() << endl;
	cout << "   number of active feasible solutions (% of active):    " << active_feas << "(" << 100 * active_feas / _dv.shape().first << "%)" << endl;
	save(_dv, _oe);
}

/**
 * @brief Save.
 *
 * @param _dv Description.
 * @param _oe Description.
 * @param save_base Description.
 */
void SeqQuadProgram::save(ParameterEnsemble& _dv, ObservationEnsemble& _oe, bool save_base)
{
	ofstream& frec = file_manager.rec_ofstream();
	stringstream ss;
	if (pest_scenario.get_pestpp_options().get_save_binary())
	{
		ss << file_manager.get_base_filename() << "." << iter << ".obs.jcb";
		_oe.to_binary(ss.str());
	}
	else
	{
		ss << file_manager.get_base_filename() << "." << iter << ".obs.csv";
		_oe.to_csv(ss.str());
	}
	frec << "      obs ensemble saved to " << ss.str() << endl;
	cout << "      obs ensemble saved to " << ss.str() << endl;
	ss.str("");
	if (pest_scenario.get_pestpp_options().get_save_binary())
	{
		ss << file_manager.get_base_filename() << "." << iter << ".par.jcb";
		_dv.to_binary(ss.str());
	}
	else
	{
		ss << file_manager.get_base_filename() << "." << iter << ".par.csv";
		_dv.to_csv(ss.str());
	}
	if (save_base)
	{
		save_real_par_rei(pest_scenario, _dv, _oe, output_file_writer, file_manager, iter);
		save_real_par_rei(pest_scenario, _dv, _oe, output_file_writer, file_manager, -1);
	}
	//ss << file_manager.get_base_filename() << "." << iter << ".par.csv";
	//dv.to_csv(ss.str());
	frec << "      par ensemble saved to " << ss.str() << endl;
	cout << "      par ensemble saved to " << ss.str() << endl;

	

}

ObservationEnsemble SeqQuadProgram::run_candidate_ensemble(ParameterEnsemble& dv_candidates)
{
	run_mgr_ptr->reinitialize();
	ofstream &frec = file_manager.rec_ofstream();
	stringstream ss;
	ss << "queuing " << dv_candidates.shape().first << " candidate solutions";
	performance_log->log_event(ss.str());
	run_mgr_ptr->reinitialize();
	
	map<int, int> real_run_ids;
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
	
	ObservationEnsemble _oe(&pest_scenario, &rand_gen);
	_oe.reserve(dv_candidates.get_real_names(), pest_scenario.get_ctl_ordered_obs_names());

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
		vector<string> par_real_names = dv_candidates.get_real_names();
		vector<string> obs_real_names = _oe.get_real_names();
		vector<string> failed_par_names, failed_obs_names;
		string oname, pname;
		ss << "the following dv candidate runs failed -->";
		for (auto& i : failed_real_indices)
		{
			if (i >= 0 && i < par_real_names.size() && i < obs_real_names.size())
			{
				pname = par_real_names[i];
				oname = obs_real_names[i];

				bool par_exists = (find(par_real_names.begin(), par_real_names.end(), pname) != par_real_names.end());
				bool obs_exists = (find(obs_real_names.begin(), obs_real_names.end(), oname) != obs_real_names.end());

				if (par_exists && obs_exists)
				{
					failed_par_names.push_back(pname);
					failed_obs_names.push_back(oname);
					ss << pname << ":" << oname << ',';
				}
				else
				{
					message(1, "WARNING: failed index " + to_string(i) + " name mismatch, skipping drop");
				}
			}
			else
			{
				message(1, "WARNING: failed index " + to_string(i) + " out of bounds (par_size=" +
					to_string(par_real_names.size()) + ", obs_size=" + to_string(obs_real_names.size()) + "), skipping");
			}
		}
		string s = ss.str();
		message(1, s);
		if (failed_par_names.size() > 0 && failed_obs_names.size() > 0)
		{
			if (failed_real_indices.size() == _oe.shape().first)
			{
				message(0, "WARNING: all dv candidate runs failed");
				_oe = ObservationEnsemble(&pest_scenario);
			}
			else
			{
				if (failed_par_names.size() > 0 && failed_obs_names.size() > 0 &&
					failed_par_names.size() == failed_obs_names.size())
				{
					performance_log->log_event("dropping failed realizations");
					try
					{
						_oe.drop_rows(failed_obs_names);
						dv_candidates.drop_rows(failed_par_names);
					}
					catch (const exception& e)
					{
						message(0, "ERROR dropping failed runs: " + string(e.what()));
						throw_sqp_error("error dropping failed realizations: " + string(e.what()));
					}
					catch (...)
					{
						message(0, "ERROR dropping failed runs: unknown exception");
						throw_sqp_error("error dropping failed realizations: unknown exception");
					}
				}
				else
				{
					message(1, "WARNING: cannot drop failed runs - name vectors empty or mismatched (par_size=" +
						to_string(failed_par_names.size()) + ", obs_size=" + to_string(failed_obs_names.size()) + ")");
				}
			}
		}
		else
		{
			performance_log->log_event("dropping failed realizations");
			//_oe.drop_rows(failed_real_indices);
			//pe_lams[i].drop_rows(failed_real_indices);
			_oe.drop_rows(failed_obs_names);
			dv_candidates.drop_rows(failed_par_names);
			//update scale_vals 
			/*vector<double> new_scale_vals;
			for (int i = 0; i < real_names.size(); i++)
				if (find(failed_real_indices.begin(), failed_real_indices.end(), i) == failed_real_indices.end())
					new_scale_vals.push_back(scale_vals[i]);
			scale_vals = new_scale_vals;*/
		}
	}
	
	
	return _oe;
}

/**
 * @brief Queue chance runs.
 */
void SeqQuadProgram::queue_chance_runs()
{
	/* queue up chance-related runs using the class attributes dp and op*/
	if (pest_scenario.get_control_info().noptmax == 0)
		return;
	stringstream ss;
	if (constraints.should_update_chance(iter))
	{
		if (use_ensemble_grad)
		{
			if (chancepoints == chancePoints::ALL)
			{
				message(1, "queueing up chance runs using nested chance points");
				constraints.add_runs(iter, dv, current_obs, run_mgr_ptr);
			}
			else
			{
				//just use dp member nearest the mean dec var values
				dv.transform_ip(ParameterEnsemble::transStatus::NUM);
				vector<double> t = dv.get_mean_stl_var_vector();
				Eigen::VectorXd dv_mean = stlvec_2_eigenvec(t);
				t.resize(0);
				ss << "queueing up chance runs using mean decision variables";
				message(1, ss.str());
				Parameters pars = pest_scenario.get_ctl_parameters();
				pest_scenario.get_base_par_tran_seq().ctl2numeric_ip(pars);
				pars.update_without_clear(dv.get_var_names(), dv_mean);
				Observations obs = pest_scenario.get_ctl_observations();
				pest_scenario.get_base_par_tran_seq().numeric2ctl_ip(pars);
				constraints.add_runs(iter, pars, obs, run_mgr_ptr);
			}
		}
		else
		{
			message(1, "queuing chance runs");
			constraints.add_runs(iter, current_ctl_dv_values, current_obs, run_mgr_ptr);
		}
	}
}



/**
 * @brief Run ensemble.
 *
 * @param _pe Description.
 * @param _oe Description.
 * @param real_idxs Description.
 *
 * @return Description.
 */
vector<int> SeqQuadProgram::run_ensemble(ParameterEnsemble &_pe, ObservationEnsemble &_oe, const vector<int> &real_idxs)
{
	run_mgr_ptr->reinitialize();
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
	_oe.reserve(_pe.get_real_names(), pest_scenario.get_ctl_ordered_obs_names());
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

	if (pest_scenario.get_pestpp_options().get_ies_debug_fail_remainder()) {
		failed_real_indices.push_back(0);
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


/**
 * @brief Finalize.
 */
void SeqQuadProgram::finalize()
{

}

vector<int> SeqQuadProgram::get_subset_idxs(int size, int nreal_subset)
{
	vector<int> subset_idxs;
	if ((!use_subset) || (nreal_subset >= size))
	{
		for (int i = 0; i < size; i++)
			subset_idxs.push_back(i);
		return subset_idxs;
	}

	vector<string> dv_real_names = dv.get_real_names();

	map<string, map<string, double>> violations_nominal;
	if (constraints.get_use_chance() && (sqp_risk != 0.5))
	{
		violations_nominal = constraints.get_ensemble_violations_map(dv, oe, 0.0, true, &oe, sqp_risk);
	}
	else
	{
		violations_nominal = constraints.get_ensemble_violations_map(dv, oe, 0.0, true);
	}

	vector<pair<int, double>> idx_vsum_pairs;
	for (size_t i = 0; i < size; ++i)
	{
		const string& rname = dv_real_names[i];
		if (rname == BASE_REAL_NAME)
			continue;

		double vsum = 0;
		auto viol_it = violations_nominal.find(rname);
		if (viol_it != violations_nominal.end())
		{
			for (const auto& v : viol_it->second)
				vsum += v.second;
		}

		idx_vsum_pairs.push_back({ static_cast<int>(i), vsum });
	}

	sort(idx_vsum_pairs.begin(), idx_vsum_pairs.end(),
		[](const pair<int, double>& a, const pair<int, double>& b) {
			return a.second < b.second;
		});

	set<int> selected_set;

	for (const auto& pair : idx_vsum_pairs)
	{
		if (subset_idxs.size() >= nreal_subset)
			break;

		if (selected_set.find(pair.first) == selected_set.end())
		{
			subset_idxs.push_back(pair.first);
			selected_set.insert(pair.first);
		}
	}

	sort(subset_idxs.begin(), subset_idxs.end());
	return subset_idxs;

}

void CovMatAdap::initialize(int n_params, int _num_reals)
{
	lambda = _num_reals;
	clear_archives();
	if (pest_scenario_ptr->get_pestpp_options().get_sqp_cma_parent_num() <= 0)
		mu = lambda / 4;
	else
		mu = pest_scenario_ptr->get_pestpp_options().get_sqp_cma_parent_num();

	sigma = 1.0;

	m = Eigen::VectorXd::Zero(n_params);
	C = Eigen::MatrixXd::Identity(n_params, n_params);

	pc = Eigen::VectorXd::Zero(n_params);
	ps = Eigen::VectorXd::Zero(n_params);

	B = Eigen::MatrixXd::Identity(n_params, n_params);
	D = Eigen::VectorXd::Ones(n_params);

	weights.resize(max(0,_num_reals));
	for (int i = 0; i < mu; i++) {
		weights[i] = log(mu + 0.5) - log(i + 1.0);
	}

	for (int i = mu; i < lambda; i++) {
		weights[i] = 0.0;
	}

	double sum_positive_weights = 0.0;
	for (int i = 0; i < mu; i++) {
		sum_positive_weights += weights[i];
	}

	for (int i = 0; i < mu; i++) {
		weights[i] /= sum_positive_weights;
	}

	double sum_squared_weights = 0.0;
	for (int i = 0; i < mu; i++) {
		sum_squared_weights += weights[i] * weights[i];
	}
	mu_eff = 1.0 / sum_squared_weights;

	if (pest_scenario_ptr->get_pestpp_options().get_sqp_cma_c1() != -1)
		c_1 = pest_scenario_ptr->get_pestpp_options().get_sqp_cma_c1();
	else
		c_1 = 2.0 / ((n_params + 1.3) * (n_params + 1.3) + mu_eff);

	if (pest_scenario_ptr->get_pestpp_options().get_sqp_cma_cmu() != -1)
		c_mu = pest_scenario_ptr->get_pestpp_options().get_sqp_cma_cmu();
	else
		c_mu = min(1 - c_1, 2.0 / (n_params + sqrt(2.0)) + min(1.0, 2.0 * mu_eff / (n_params + 2.0)) * (1.0 / (n_params + 2.0)));
	
	if (pest_scenario_ptr->get_pestpp_options().get_sqp_cma_cc() != -1)
		c_c = pest_scenario_ptr->get_pestpp_options().get_sqp_cma_cc();
	else
		c_c = (4.0 + mu_eff / n_params) / (n_params + 4.0 + 2.0 * mu_eff / n_params);
	
	c_sigma = (mu_eff + 2.0) / (n_params + mu_eff + 5.0);
	d_sigma = 1.0 + 2.0 * max(0.0, sqrt((mu_eff - 1.0) / (n_params + 1.0)) - 1.0) + c_sigma;
	chi_n = sqrt(n_params) * (1.0 - 1.0 / (4.0 * n_params) + 1.0 / (21.0 * n_params * n_params));

}

void CovMatAdap::update(Parameters prev_m, Parameters curr_m, int iter) 
{
	ofstream& frec = file_manager->rec_ofstream();
	CovMetrics metrics_prior = compute_cov_metrics();
	if (iter == 1)
		metrics_init = metrics_prior;

	ParameterEnsemble U = sorted_dp_archive;
	vector<string> par_names = U.get_var_names();
	m = curr_m.get_data_eigen_vec(par_names);
	vector<string> mu_best_real_names;
	int i;

	if (U.shape().first != 0)
	{
		i = 0;
		for (auto r : U.get_real_names())
		{
			i++;
			if (i > mu)
				break;
			mu_best_real_names.push_back(r);

		}
		U.keep_rows(mu_best_real_names);
	}

	if (U.shape().first != 0)
	{
		if (pest_scenario_ptr->get_pestpp_options().get_sqp_cma_stepsize_control())
		{
			Eigen::VectorXd y_w = (m - prev_m.get_data_eigen_vec(par_names)) / sigma;
			Eigen::VectorXd ww = D.cwiseInverse().cwiseSqrt().asDiagonal() * (B.transpose() * y_w);
			ps = (1.0 - c_sigma) * ps + sqrt(c_sigma * (2.0 - c_sigma) * mu_eff) * (B * ww);
			sigma = sigma * exp((c_sigma / d_sigma) * (ps.norm() / chi_n - 1.0));
		}

		Eigen::MatrixXd U_anoms = U.get_eigen_anomalies(vector<string>(), vector<string>()) / sigma;

		Eigen::MatrixXd rank_mu_update = Eigen::MatrixXd::Zero(C.rows(), C.cols());
		for (int i = 0; i < mu; i++)
		{
			Eigen::VectorXd y_i = U_anoms.row(i).transpose();
			rank_mu_update += weights[i] * (y_i * y_i.transpose());
		}
		rank_mu_update = c_mu * rank_mu_update;

		pc = (1 - c_c) * pc + sqrt(c_c * (2 - c_c) * mu_eff) * (m - prev_m.get_data_eigen_vec(par_names)) / sigma;
		Eigen::MatrixXd rank_one_update = c_1 * (pc * pc.transpose());

		C = (1.0 - c_1 - c_mu) * C + rank_mu_update + rank_one_update;

		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(C);
		B = eigensolver.eigenvectors();
		D = eigensolver.eigenvalues();

		for (int i = 0; i < D.size(); i++)
		{
			D(i) = max(D(i), D.maxCoeff() * 1E-12);
		}
	}
	else
	{
		frec << "...nothing to learn for covariance here...skipping CMA update" << endl;
		cout << "...nothing to learn for covariance here...skipping CMA update" << endl;
	}

	CovMetrics metrics_post = compute_cov_metrics();
	const double max_cond_num = pest_scenario_ptr->get_pestpp_options().get_sqp_max_reinflation_cond_num();

	if (metrics_post.condition_number > max_cond_num)
	{

		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig_C(C);
		Eigen::VectorXd eigenvals = eig_C.eigenvalues();
		double lambda_max = eigenvals.maxCoeff();
		double lambda_min = eigenvals.minCoeff();

		double min_eig_floor = lambda_max / max_cond_num;

		double delta = 0.0;
		if (lambda_min < min_eig_floor) 
		{
			delta = min_eig_floor - lambda_min + 1E-12;
			C += delta * Eigen::MatrixXd::Identity(C.rows(), C.cols());

			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(C);
			B = eigensolver.eigenvectors();
			D = eigensolver.eigenvalues();
		}

		ofstream& frec = file_manager->rec_ofstream();
		frec << "...WARNING: ensemble is shrinking too much. Regularizing covariance by delta: " << delta << endl;
		cout << "...WARNING: ensemble is shrinking too much. Regularizing covariance by delta: " << delta << endl;

		metrics_post = compute_cov_metrics();
	}
	cma_update_summary = report_cma_metrics(metrics_prior, metrics_post, iter); 
}

void CovMatAdap::reinflate_C(double reinflation_factor, bool reset_corr, double max_cond_num)
{
	if (reinflation_factor < 0.0)
	{
		CovMetrics metrics = compute_cov_metrics();
		double factor = metrics.determinant / pow(10.0, floor(log10(fabs(metrics.determinant))));
		reinflate_C(factor, reset_corr, max_cond_num);
	}
	else
	{
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(C);
		B = eigensolver.eigenvectors();
		D = eigensolver.eigenvalues();

		if (reset_corr)
		{
			Eigen::MatrixXd C_diag = C.diagonal().asDiagonal();
			C = C_diag * reinflation_factor;

			eigensolver.compute(C);
			B = eigensolver.eigenvectors();
			D = eigensolver.eigenvalues();
		}
		else if (max_cond_num > 0.0)
		{
			std::vector<std::pair<double, int>> eigval_idx;
			for (int i = 0; i < D.size(); i++)
			{
				eigval_idx.push_back(std::make_pair(D(i), i));
			}
			std::sort(eigval_idx.begin(), eigval_idx.end(), std::greater<std::pair<double, int>>());

			double lambda_max = eigval_idx[0].first;
			double lambda_min = eigval_idx.back().first;
			double current_cond_num = lambda_max / lambda_min;

			for (int i = 0; i < D.size(); i++)
			{
				D(i) *= reinflation_factor;
			}
		
			if (current_cond_num > max_cond_num)
			{
				double new_lambda_max = lambda_max * reinflation_factor;
				double target_lambda_min = new_lambda_max / max_cond_num;

					double min_eigenvalue_floor = target_lambda_min;

				for (int i = 0; i < D.size(); i++)
				{
						double scaled_eigenvalue = D(i);

						if (scaled_eigenvalue < min_eigenvalue_floor)
							D(i) = 0.5 * (scaled_eigenvalue + min_eigenvalue_floor);

						D(i) = max(D(i), 1e-12);
			}

			Eigen::MatrixXd D_matrix = D.asDiagonal();
					C = B * D_matrix * B.transpose();
					C = (C + C.transpose()) / 2.0;
			}
			else
				{
					Eigen::MatrixXd D_matrix = D.asDiagonal();
					C = B * D_matrix * B.transpose();
					C = (C + C.transpose()) / 2.0;
				}
		}
		else
		{
				C *= reinflation_factor;

			eigensolver.compute(C);
			B = eigensolver.eigenvectors();
			D = eigensolver.eigenvalues();
		}

		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> final_eigensolver(C);
		B = final_eigensolver.eigenvectors();
		D = final_eigensolver.eigenvalues();

		for (int i = 0; i < D.size(); i++)
		{
				D(i) = max(D(i), 1e-12);
		}

		Eigen::MatrixXd D_final = D.asDiagonal();
		C = B * D_final * B.transpose();
		C = (C + C.transpose()) / 2.0;
	}
}
void CovMatAdap::clear_archives()
{
	sorted_obj_map.clear();
	sorted_dp_archive = ParameterEnsemble(pest_scenario_ptr, rand_gen_ptr);
}

void CovMatAdap::update_archives(const ParameterEnsemble& pe, map<string, double> obj_map, map<string, double> viol_map, string tag, bool clear)
{
	map<string, double> unique_obj_map;
	set<double> seen_values;
	ParameterEnsemble curr_pe = pe;
	
	if (clear) 
	{
		clear_archives();
		sorted_dp_archive = curr_pe.zeros_like(0);
		unique_obj_map = obj_map;
	}
	else
	{
		if (sorted_dp_archive.shape().second == 0) 
			sorted_dp_archive = curr_pe.zeros_like(0);

		for (auto& o : obj_map)
		{
			bool is_duplicate = false;
			for (auto& s : sorted_obj_map)
			{
				if (fabs(o.second - s.second) < 1E-10)
				{
					is_duplicate = true;
					break;
				}
			}

			if (!is_duplicate)
				unique_obj_map[o.first] = o.second;
		}
	}

	for (auto o : unique_obj_map)
	{
		sorted_obj_map[tag + "|" + o.first] = o.second;
		sorted_dp_archive.append(tag + "|" + o.first, curr_pe.get_real_vector(o.first));
	}

	if (sorted_obj_map.size() > 0)
	{
		vector<pair<string, double>> sorted_obj_map_vec(sorted_obj_map.begin(), sorted_obj_map.end());
		sort(sorted_obj_map_vec.begin(), sorted_obj_map_vec.end(),
			[](const pair<string, double>& a, const pair<string, double>& b) {
				return a.second < b.second;
			});

		vector<string> sorted_names_from_obj;
		sorted_obj_map.clear();
		int i = 0;
		for (const auto& pair : sorted_obj_map_vec)
		{
			i++;
			if (i > lambda)
				break;
			sorted_names_from_obj.push_back(pair.first);
			sorted_obj_map[pair.first] = pair.second;

		}
		sorted_dp_archive.keep_rows(sorted_names_from_obj, true);
		sorted_dp_archive.reorder(sorted_names_from_obj, curr_pe.get_var_names(), true);
	}

	if (sorted_obj_map.size() == 0)
		throw runtime_error("no members in sorted obj maps after CovMatAdap::update_archives()");

}

CovMatAdap::CovMetrics CovMatAdap::compute_cov_metrics() const
{
	CovMetrics metrics;

	metrics.trace = C.trace();
	metrics.determinant = C.determinant();
	metrics.frobenius_norm = C.norm();
	metrics.max_eigenvalue = D.maxCoeff();
	metrics.min_eigenvalue = D.minCoeff();
	metrics.condition_number = metrics.max_eigenvalue / (metrics.min_eigenvalue + 1E-10);

	return metrics;
}

string CovMatAdap::report_cma_metrics(const CovMetrics& prior, const CovMetrics& post, int iter)
{

	trace_ratio = post.trace / (prior.trace + 1E-30);
	det_ratio = post.determinant / (prior.determinant + 1E-30);
	frobenius_ratio = post.frobenius_norm / (prior.frobenius_norm + 1E-30);
	max_eigenval_ratio = post.max_eigenvalue / (prior.max_eigenvalue + 1E-30);

	int criteria_sat = 0;
	double tol = 0.05;
	trace_ratio_0 = post.trace / (metrics_init.trace + 1E-30);
	if (trace_ratio_0 < tol)
		criteria_sat++;
	det_ratio_0 = post.determinant / (metrics_init.determinant + 1E-30);
	if (det_ratio_0 < tol)
		criteria_sat++;
	frobenius_ratio_0 = post.frobenius_norm / (metrics_init.frobenius_norm + 1E-30);
	if (frobenius_ratio_0 < tol)
		criteria_sat++;
	max_eigenval_ratio_0 = post.max_eigenvalue / (metrics_init.max_eigenvalue + 1E-30);
	if (max_eigenval_ratio_0 < tol)
		criteria_sat++;

	stringstream ss;
	ss << left; 
	ss << setw(20) << "  " 
		<< setw(14) << "Prior"
		<< setw(14) << "Post"
		<< setw(10) << "Po/Pr"
		<< setw(10) << "Po/Pr_0" << endl;

	ss << scientific << setprecision(4);
	ss << "   Trace           " << setw(12) << prior.trace
		<< "   " << setw(12) << post.trace
		<< "   " << fixed << setprecision(6) 
		<< setw(12) << trace_ratio << setw(10) << trace_ratio_0 << endl;

	ss << scientific << setprecision(4);
	ss << "   Determinant     " << setw(12) << prior.determinant
		<< "   " << setw(12) << post.determinant
		<< "   " << fixed << setprecision(6) 
		<< setw(12) << det_ratio << setw(10) << det_ratio_0 << endl;

	ss << scientific << setprecision(4);
	ss << "   Frobenius Norm  " << setw(12) << prior.frobenius_norm
		<< "   " << setw(12) << post.frobenius_norm
		<< "   " << fixed << setprecision(6) 
		<< setw(12) << frobenius_ratio << setw(10) << frobenius_ratio_0 << endl;

	ss << scientific << setprecision(4);
	ss << "   Max Eigenvalue  " << setw(12) << prior.max_eigenvalue
		<< "   " << setw(12) << post.max_eigenvalue
		<< "   " << fixed << setprecision(6) 
		<< setw(12) << max_eigenval_ratio << setw(10) << max_eigenval_ratio_0 << endl;

	ss << scientific << setprecision(4);
	ss << "   Condition No.   " << setw(12) << prior.condition_number
		<< "   " << setw(12) << post.condition_number
		<< "   " << setw(12) << " " << setw(10) << " " << endl;
	ss << "   number of cov-based criteria satisfied: " << criteria_sat << endl << endl;

	ss << left << setw(20) << "   CMA parameter"
		<<setw(18) << "Value" << endl;
	ss << scientific << setprecision(6);
	ss << setw(20) << "      c_c" << fixed << setprecision(6) << c_c << endl;
	ss << setw(20) << "      c_1" << fixed << setprecision(6) << c_1 << endl;
	ss << setw(20) << "      c_mu" << fixed << setprecision(6) << c_mu << endl;
	ss << setw(20) << "      c_sigma" << fixed << setprecision(6) << c_sigma << endl;
	ss << setw(20) << "      d_sigma" << fixed << setprecision(6) << d_sigma << endl;
	ss << setw(20) << "      chi_n" << fixed << setprecision(6) << chi_n << endl;
	ss << setw(20) << "      sigma" << fixed << setprecision(6) << sigma << endl;
	ss << setw(20) << "      lambda" << fixed << setprecision(6) << lambda << endl;
	ss << setw(20) << "      mu" << fixed << setprecision(6) << mu << endl;
	ss << setw(20) << "      mu_eff" << fixed << setprecision(6) << mu_eff << endl;
	return ss.str();

}

bool CovMatAdap::should_terminate()
{
	if (trace_ratio_0 < 0.05 && det_ratio_0 < 0.05 && frobenius_ratio_0 < 0.05 && max_eigenval_ratio_0 < 0.05)
		return true;
	else
		return false;

}

ParameterEnsemble CovMatAdap::generate_population(Parameters& _curr_m, ParameterEnsemble _dv) 
{	
	vector<string> rnames;
	vector<string> parnames = _dv.get_var_names();
	ParameterEnsemble new_reals = _dv;

	rnames = new_reals.get_real_names();
	auto it = find(rnames.begin(), rnames.end(), BASE_REAL_NAME);
	if (it != rnames.end())
		rnames.erase(it);
	new_reals.keep_rows(rnames, true);

	const ParameterInfo& par_info = pest_scenario_ptr->get_ctl_parameter_info();
	
	vector<string> target_rnames = new_reals.get_generic_real_names(lambda);
	vector<string> missing_rnames;
	set<string> current_rnames_set(rnames.begin(), rnames.end());
	for (const auto& target_name : target_rnames)
	{
		if (current_rnames_set.find(target_name) == current_rnames_set.end())
		{
			missing_rnames.push_back(target_name);
		}
	}

	for (const auto& missing_name : missing_rnames)
	{
		Eigen::VectorXd x(parnames.size());
		int draws = 1;
		//jwhite: is this big enough considering we might be at/near bounds
		//in many cases, esp in high dimensions?
		const int max_draws = 1000;
		bool found = false;

		while (!found)
		{
			if (draws > max_draws) 
			{
				for (int j = 0; j < parnames.size(); ++j) 
				{
					const ParameterRec* par_rec = par_info.get_parameter_rec_ptr(parnames[j]);
					double lbnd = par_rec->lbnd;
					double ubnd = par_rec->ubnd;

					if (par_rec->tranform_type == ParameterRec::TRAN_TYPE::LOG) 
					{
						lbnd = log10(lbnd);
						ubnd = log10(ubnd);
					}

					x(j) = max(lbnd, min(ubnd, x(j)));
				}
				break;
			}

			Eigen::VectorXd z(m.size());
			for (int j = 0; j < m.size(); j++) 
			{
				z(j) = draw_standard_normal(*rand_gen_ptr);
			}

			Eigen::VectorXd y = B * (D.cwiseSqrt().asDiagonal() * z);
			x = _curr_m.get_data_eigen_vec(parnames) + sigma * y;

			found = true;
			for (int j = 0; j < parnames.size(); j++) 
			{
				const ParameterRec* par_rec = par_info.get_parameter_rec_ptr(parnames[j]);
				double lbnd = par_rec->lbnd;
				double ubnd = par_rec->ubnd;

				if (par_rec->tranform_type == ParameterRec::TRAN_TYPE::LOG) 
				{
					lbnd = log10(lbnd);
					ubnd = log10(ubnd);
				}

				if (x(j) < lbnd || x(j) > ubnd) {
					found = false;
					break;
				}
			}
			draws++;
		}

		new_reals.append(missing_name, x);
		rnames.push_back(missing_name);
	}

	if (rnames.size() != target_rnames.size() || rnames != target_rnames)
	{
		set<string> current_set(rnames.begin(), rnames.end());
		for (const auto& target_name : target_rnames)
		{
			if (current_set.find(target_name) == current_set.end())
			{
				stringstream ss;
				ss << "generate_population() error: target name '" << target_name
					<< "' not found after restoration. Current size: " << rnames.size()
					<< ", target size: " << lambda;
				throw runtime_error(ss.str());
			}
		}

		new_reals.reorder(target_rnames, vector<string>());
		rnames = target_rnames;
	}

	if (rnames.size() != lambda)
	{
		stringstream ss;
		ss << "generate_population() error: rnames.size()=" << rnames.size()
			<< " != lambda=" << lambda << " after restoration";
		throw runtime_error(ss.str());
	}

	
	if (new_reals.shape().first != lambda)
	{
		stringstream ss;
		ss << "generate_population() error: new_reals.shape().first="
			<< new_reals.shape().first << " != lambda=" << lambda;
		throw runtime_error(ss.str());
	}


	for (int i = 0; i < lambda; i++) 
	{
		if (i >= rnames.size())
		{
			stringstream ss;
			ss << "generate_population() error: index " << i
				<< " >= rnames.size()=" << rnames.size();
			throw runtime_error(ss.str());
		}

		Eigen::VectorXd x(parnames.size());
		int draws = 1;
		//jwhite: is this big enough considering we might be at/near bounds
		//in many cases, esp in high dimensions?
		const int max_draws = 1000;

		bool found = false;
		while (!found) 
		{
			if (draws > max_draws) {
				for (int j = 0; j < parnames.size(); ++j) {
					const ParameterRec* par_rec = par_info.get_parameter_rec_ptr(parnames[j]);
					double lbnd = par_rec->lbnd;
					double ubnd = par_rec->ubnd;

					if (par_rec->tranform_type == ParameterRec::TRAN_TYPE::LOG) {
						lbnd = log10(lbnd);
						ubnd = log10(ubnd);
					}

					x(j) = max(lbnd, min(ubnd, x(j)));
				}
				break;
			}

			Eigen::VectorXd z(m.size());
			for (int j = 0; j < m.size(); ++j) {
				z(j) = draw_standard_normal(*rand_gen_ptr);
			}

			Eigen::VectorXd y = B * (D.cwiseSqrt().asDiagonal() * z);
			x = _curr_m.get_data_eigen_vec(parnames) + sigma * y;

			found = true;;
			for (int j = 0; j < parnames.size(); ++j) {
				const ParameterRec* par_rec = par_info.get_parameter_rec_ptr(parnames[j]);
				double lbnd = par_rec->lbnd;
				double ubnd = par_rec->ubnd;

				if (par_rec->tranform_type == ParameterRec::TRAN_TYPE::LOG) {
					lbnd = log10(lbnd);
					ubnd = log10(ubnd);
				}

				if (x(j) < lbnd || x(j) > ubnd) {
					found = false;
					break;
				}
			}
			draws++;
		}
		new_reals.update_real_ip(rnames[i], x);
		
	}

	vector<string> all_adj_par_names = pest_scenario_ptr->get_ctl_ordered_adj_par_names();
	if (all_adj_par_names.size() > parnames.size())
	{
		Parameters all_pars = pest_scenario_ptr->get_ctl_parameters();
		ParamTransformSeq pts = pest_scenario_ptr->get_base_par_tran_seq();

		if (new_reals.get_trans_status() == ParameterEnsemble::transStatus::NUM)
			pts.ctl2numeric_ip(all_pars);

		ParameterEnsemble expanded_reals(pest_scenario_ptr, rand_gen_ptr);
		expanded_reals.reserve(rnames, all_adj_par_names);
		expanded_reals.set_trans_status(new_reals.get_trans_status());

		for (const auto& rname : rnames) 
		{
			Parameters real_pars = all_pars;

			real_pars.update_without_clear(parnames, new_reals.get_real_vector(rname));

			Eigen::VectorXd full_vec = real_pars.get_data_eigen_vec(all_adj_par_names);
			expanded_reals.update_real_ip(rname, full_vec);
		}

		return expanded_reals;
	}

	return new_reals;
}
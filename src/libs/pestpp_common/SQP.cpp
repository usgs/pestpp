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


bool SqpFilter::accept(double obj_val, double violation_val, int iter, double alpha,bool keep)
{
	FilterRec candidate{ obj_val, violation_val,iter,alpha };
	if (obj_viol_pairs.size() == 0)
	{
		obj_viol_pairs.insert(candidate);
		return true;
	}
	//I think its cheaper to combine the tols with the candidate, rather adding them to every 
	//existing pair...
	if (minimize)
		candidate.obj_val *= (1 + obj_tol);
	else
		candidate.obj_val *= (1 - obj_tol);
	candidate.viol_val *= (1 + viol_tol);
	
	bool accept = true;
	for (auto& p : obj_viol_pairs)
		if (!first_partially_dominates_second(candidate, p))
		{
			accept = false;
			break;
		}
	if ((keep) && (accept))
    {
	    //cout << "obj:" << obj_val << ", viol:" << violation_val << ", alpha:" << alpha << endl;
	    obj_viol_pairs.insert(candidate);
    }
	return accept;
}


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

void SqpFilter::report(ofstream& frec, int iter)
{
    frec << "...SQP filter members (" << obj_viol_pairs.size() <<") for iteration " << iter << ":" << endl << "    obj, violation" << endl;
    double omin = 1.0e+300,omax = -1e+300,vmin = 1e+300,vmax = -1e+300;
    for (auto& fr : obj_viol_pairs)
    {
        frec << setw(6) << setprecision(3) << fr.obj_val << "," << fr.viol_val << endl;
        omin = min(fr.obj_val,omin);
        omax = max(fr.obj_val,omax);
        vmin = min(fr.viol_val,vmin);
        vmax = max(fr.viol_val,vmax);
    }
    stringstream ss;
    ss.str("");
    ss << endl << "... filter summary with " << obj_viol_pairs.size() << " pairs for iteration " << iter << ":" << endl;
    ss << "         obj min: " <<  setw(10) << omin << endl;
    ss << "         obj max: " << setw(10) << omax << endl;
    ss << "   violation min: " << setw(10) << vmin << endl;
    ss << "   violation max: " << setw(10) << vmax << endl;
    ss << endl;

    frec << ss.str();
    cout << ss.str();

}

bool SqpFilter::update(double obj_val, double violation_val, int iter, double alpha)
{
    //check if this candidate is nondom
	//bool acc = accept(obj_val, violation_val,iter, alpha);
	//if (!acc)
	//	return false;
	FilterRec candidate;
	candidate.obj_val = obj_val;
	candidate.viol_val = violation_val;
	candidate.iter = iter;
	candidate.alpha = alpha;
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
	        if (i == j)
	            continue;
	        if (first_strictly_dominates_second(*first,*second)) {
                i_is_dominated = true;
                break;
            }
	        second++;
        }
	    if (!i_is_dominated)
        {
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
				//dv.enforce_bounds(performance_log, pest_scenario.get_pestpp_options().get_ies_enforce_chglim());
			    message(1, "TODO");
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
		vector<int> drop{ _dv.shape().first - 1 };
		_dv.drop_rows(drop);
		_dv.append(BASE_REAL_NAME, pars);
	}

	//check that 'base' isn't already in ensemble
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
			vector<int> drop{ _oe.shape().first - 1 };
			_oe.drop_rows(drop);
			_oe.append(BASE_REAL_NAME, obs);
		}
	}
}

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

	iter = 0;

	if (pest_scenario.get_control_info().noptmax == 0)
	{
		message(0, "'noptmax'=0, running control file parameter values and quitting");
		
		current_ctl_dv_values = pest_scenario.get_ctl_parameters();
		ParamTransformSeq pts = dv.get_par_transform();

		ParameterEnsemble _pe(&pest_scenario, &rand_gen);
		_pe.reserve(vector<string>(), pest_scenario.get_ctl_ordered_par_names());
		_pe.set_trans_status(ParameterEnsemble::transStatus::CTL);
		_pe.append("BASE", current_ctl_dv_values);
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


	message(1, "using the following upgrade vector scale (e.g. 'line search') values:", ppo->get_sqp_scale_facs());
	
	//ofstream &frec = file_manager.rec_ofstream();
	last_best = 1.0E+30;
	last_viol = 0.0;
	
	warn_min_reals = 10;
	error_min_reals = 2;
	
	//vector<double> scale_facs = pest_scenario.get_pestpp_options().get_lambda_scale_vec();
	//message(1, "using scaling factors: ", scale_facs);
	set<string> passed = ppo->get_passed_args();
	if (passed.find("SQP_SCALE_FACS") == passed.end())
	{
	    if ((use_ensemble_grad) && (SOLVE_EACH_REAL))
        {
	        message(1,"'sqp_scale_facs' not passed, using ensemble gradient, and solving each real, resetting scale facs");
	        vector<double> new_scale_facs{0.000001,0.00001,0.0005,0.01,.1};
	        message(1,"new sqp_scale_facs",new_scale_facs);
	        ppo->set_sqp_scale_facs(new_scale_facs);
        };
	}

	
	message(1, "max run fail: ", ppo->get_max_run_fail());

	//TODO: update sanity checks for SQP context
	//check that if using fd, chance points == single
	use_ensemble_grad = false;
	if (ppo->get_sqp_num_reals() > 0)
		use_ensemble_grad = true;
	sanity_checks();

	bool echo = false;
	if (verbose_level > 1)
		echo = true;

	initialize_parcov();

	//these will be the ones we track...
	//this means the initial dv vals in the control file will be the "center" of the enopt ensemble
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

	

	constraints.sqp_report(iter, current_ctl_dv_values, current_obs);
	if (constraints.get_use_chance())
	{
		constraints.presolve_chance_report(iter, current_obs, true, "initial chance constraint report");
	}

	//set the initial grad vector
	message(2, "calculating initial objective function gradient");
	current_grad_vector = calc_gradient_vector(current_ctl_dv_values);
	grad_vector_map[0] = current_grad_vector;
	//todo: save and report on initial gradient - make some checks would be useful?
	
	last_best = get_obj_value(current_ctl_dv_values, current_obs);
	last_viol = constraints.get_sum_of_violations(current_ctl_dv_values, current_obs);
	message(0, "Initial phi value, infeasible value:", vector<double>{last_best,last_viol});
	best_phis.push_back(last_best);
    best_violations.push_back(last_viol);



    double v = constraints.get_sum_of_violations(current_ctl_dv_values, current_obs);
	filter.update(last_best, v, 0, -1.0);

	if (v > 0.0)
	{
	    message(0,"initial solution infeasible, seeking feasible solution");
		seek_feasible();
	}
	
	message(2, "initializing hessian matrix with identity");
	Eigen::SparseMatrix<double> h(dv_names.size(), dv_names.size());
	h.setIdentity();
	hessian = Covariance(dv_names, h);
	message(0, "initialization complete");
}

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
			ss.str("");
			ss << "existing jacobian missing the following obs constraints:" << endl;
			for (auto m : vnames)
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
	//todo: think about freezind dec vars that go out of bounds? - yuck!
	if (out_of_bounds.size() > 0)
	{
		ss.str("");
		ss << "the following decision variable are out of bounds: " << endl;
		for (auto& o : out_of_bounds)
			ss << o << ",";
		throw_sqp_error(ss.str());
	}
	//todo: mod queue chance runs for FD grad
	queue_chance_runs();
	message(2, "starting finite difference gradient perturbation runs");
	jco.make_runs(*run_mgr_ptr);
	
	success = jco.process_runs(par_trans,pgi,*run_mgr_ptr,pi,false,false);
	if (!success)
	{
		throw_sqp_error("error processing finite difference gradient perturbation runs");
	}
	//constraints.process_runs(run_mgr_ptr, iter);
	if (init_obs)
	{
		run_mgr_ptr->get_run(0, current_pars, _current_obs, false);
	}
}

void SeqQuadProgram::make_gradient_runs(Parameters& _current_dv_vals, Observations& _current_obs)
{
	stringstream ss;
	if (use_ensemble_grad)
	{
		message(1, "generating new dv ensemble at current best location");
		//draw new dv ensemble using - assuming parcov has been updated by now...
		ParameterEnsemble _dv(&pest_scenario, &rand_gen);
		Parameters dv_par = _current_dv_vals.get_subset(dv_names.begin(),dv_names.end());
		ofstream& frec = file_manager.rec_ofstream();
		_dv.draw(pest_scenario.get_pestpp_options().get_sqp_num_reals(), dv_par, parcov, performance_log, 0, frec);

		//todo: save _dv here in case something bad happens...
		ObservationEnsemble _oe(&pest_scenario, &rand_gen);
		_oe.reserve(_dv.get_real_names(), constraints.get_obs_constraint_names());
        add_current_as_bases(_dv, _oe);
		message(1, "running new dv ensemble");
		run_ensemble(_dv, _oe);
		save(_dv, _oe);
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

	//message(1, "transforming parameter ensemble to numeric");
	dv.transform_ip(ParameterEnsemble::transStatus::NUM);


	//TODO: think about what adding the base would do for SQP
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


        Eigen::VectorXd o = _oe.get_real_vector("mean");
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
	if (pest_scenario.get_pestpp_options().get_save_binary())
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
	save_real_par_rei(pest_scenario, dv, oe, output_file_writer, file_manager, iter);
	save_real_par_rei(pest_scenario, dv, oe, output_file_writer, file_manager, -1);



	pcs = ParChangeSummarizer(&dv_base, &file_manager, &output_file_writer);

	dv.transform_ip(ParameterEnsemble::transStatus::NUM);
	vector<double> vals = dv.get_mean_stl_var_vector();
	vector<string> names = dv.get_var_names();
	current_ctl_dv_values.update(names, vals);
	ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();
	pts.numeric2ctl_ip(current_ctl_dv_values);
	vals = oe.get_mean_stl_var_vector();
	names = oe.get_var_names();
	current_obs.update(names, vals);


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

bool SeqQuadProgram::update_hessian_and_grad_vector()
{
	Parameters new_grad = calc_gradient_vector(current_ctl_dv_values);
	if (!pest_scenario.get_pestpp_options().get_sqp_update_hessian())
	{
		message(2, "hessian_update is false...");

		current_grad_vector = new_grad;
		return false;
	}
	//fancy shit here...
	message(1, "starting hessian update for iteration ", iter);
	
	Eigen::VectorXd old_grad = current_grad_vector.get_data_eigen_vec(dv_names);


	//update
	current_grad_vector = new_grad;
	//if accepted, return true
	return true;
}

void SeqQuadProgram::iterate_2_solution()
{
	stringstream ss;
	ofstream &frec = file_manager.rec_ofstream();
	
	bool accept;
	n_consec_infeas = 0;
	for (int i = 0; i < pest_scenario.get_control_info().noptmax; i++)
	{
		iter++;
		message(0, "starting solve for iteration:", iter);
		ss.str("");
		ss << "starting solve for iteration: " << iter;
		performance_log->log_event(ss.str());
 		
		//solve for and test candidates
		accept = solve_new();

        //save some stuff...
        if (use_ensemble_grad)
            report_and_save_ensemble();
        else
        {
            //save par and res files for this iteration
            save_current_dv_obs();
        }
        constraints.sqp_report(iter, current_ctl_dv_values, current_obs, true);

        //report dec var change stats - only for ensemble form
        if (use_ensemble_grad)
        {
            ss.str("");
            ss << file_manager.get_base_filename() << "." << iter << ".pcs.csv";
            pcs.summarize(dv, ss.str());
        }

        //store the grad vector used for this iteration
        grad_vector_map[iter] = current_grad_vector;

		//check to break here before making more runs
		if (should_terminate())
			break;

		if (n_consec_infeas > MAX_CONSEC_INFEAS_IES)
        {
		    ss.str("");
		    ss << "number of consecutive infeasible iterations > " << MAX_CONSEC_INFEAS << ", switching to IES to seek feasibility";
		    message(0,ss.str());
		    seek_feasible();
		    n_consec_infeas = 0;
        }

        //update the underlying runs
        make_gradient_runs(current_ctl_dv_values,current_obs);
		//update the hessian
		update_hessian_and_grad_vector();

		//todo: report constraint stats
		//a la constraints.mou_report();



	}
}

bool SeqQuadProgram::should_terminate()
{
    stringstream ss;
    //todo: use ies accept fac here?
    double phiredstp = pest_scenario.get_control_info().phiredstp;
    int nphistp = pest_scenario.get_control_info().nphistp;
    int nphinored = MAX_CONSEC_PHIINC;
    bool phiredstp_sat = false, nphinored_sat = false, consec_sat = false;
    double phi, ratio, infeas;
    int count = 0;
    int nphired = 0;
    best_phi_yet = 1.0e+300;
    int best_idx_yet = -1;
    for (int i=0;i<best_phis.size();i++)
    {
        if (best_phis[i]<=best_phi_yet)
        {
            best_phi_yet = best_phis[i];
            best_violation_yet = best_violations[i];
            best_idx_yet = i;
        }
    }
    if (best_idx_yet == -1)
    {
        throw_sqp_error("something is wrong in shouuld_terminate()");
    }
    nphired = best_phis.size() - best_idx_yet;



    //todo: save and write out the current phi grad vector (maybe save all of them???)
    ss.str("");
    ss << "best phi,infeas sequence:" << endl;
    //ss << "       ";
    int ii = 0;
    for (int i=0;i<best_phis.size();i++)
    //for (auto phi : best_phis)
    {
        phi = best_phis[i];
        infeas = best_violations[i];
        ss << "    " << setw(5) << setprecision(4) << right << phi << "," << setw(5) << setprecision(4) << left << infeas << endl;
        ii++;
//        if (ii % 6 == 0) {
//            ss << endl;
//            ss << "       ";
//        }

    }
    ss << endl;
    message(0, ss.str());

    /*if ((!consec_sat )&& (best_mean_phis.size() == 0))
        return false;*/
    message(0, "phi-based termination criteria check");
    message(1, "phiredstp: ", phiredstp);
    message(1, "nphistp: ", nphistp);
    message(1, "nphinored: ", nphinored);
    message(1, "best phi yet: ", best_phi_yet);
    message(1,"number of consecutive infeasible solutions: ",n_consec_infeas);
    for (auto& phi : best_phis)
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
    if (best_phis[best_phis.size() - 1] == 0.0)
    {
        message(1, "phi is zero, all done");
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
		for (auto& ddv_val : pert_dv_values)
		{
			pert_obj_val += obj_func_coef_map[ddv_val.first] * ddv_val.second;
		}
		derv = (current_obj_val - pert_obj_val) / (dv_val - pert_val);
		grad[i] = derv;
		i++;
	}
	return grad;
}


Parameters SeqQuadProgram::calc_gradient_vector(const Parameters& _current_dv_values, string _center_on)
{
	stringstream ss;
	Eigen::VectorXd grad(dv_names.size());
    //TODO: should this be optional?

	string center_on = pest_scenario.get_pestpp_options().get_ies_center_on();
	if (!_center_on.empty())
	    center_on = _center_on;
	
	//if don't already have or if already have and exit
	//	if (LBFGS) &(num_it > 2)& (constraints False); // constraint False as need phi_grad for Lagrangian
	//	{
	//		ss.str("");
	//		ss << "(re)use grad from Wolfe testing during upgrade evaluations last iteration";
	//		string s = ss.str();
	//		message(1, s);
	//		throw_sqp_error("TODO");
	//	}
	
	if (use_ensemble_grad)
	{
		//ensemble stuff here
		//if (use_obj_obs)
		{
			// compute sample dec var cov matrix and its pseudo inverse
			// see eq (8) of Dehdari and Oliver 2012 SPE and Fonseca et al 2015 SPE
			// TODO: so can pseudo inverse: Covariance dv_cov_matrix; 
			//Eigen::MatrixXd dv_cov_matrix;
			//Eigen::MatrixXd parcov_inv;
			// start by computing mean-shifted dec var ensemble
			Eigen::MatrixXd dv_anoms = dv.get_eigen_anomalies(vector<string>(), dv_names, center_on);  // need this for both cov and cross-cov
			if (dv.shape().first > 1000)  // until we encounter
			{
				// lower rank - diag elements only
				throw_sqp_error("TODO: use dv.get_diagonal_cov matrix()? need to check for consistency if so"); 
				//parcov = dv.get_diagonal_cov_matrix();  // check ok to instantiate empricially here  // pass helpful center_on here too
				//Eigen::VectorXd parcov_inv;
				//Covariance parcov_diag;
				//parcov_diag.from_diagonal(parcov);
				//parcov_inv = parcov_diag.get_matrix().diagonal();
				//parcov_inv = parcov_inv.cwiseInverse();  // equivalent to pseudo inv?
			}
			//dv_cov_matrix = 1.0 / (dv.shape().first - 1.0) * (dv_anoms.transpose() * dv_anoms);
			//message(1, "dv_cov:", dv_cov_matrix);
			//parcov_inv = dv_cov_matrix.cwiseInverse();  // check equivalence to pseudo inv

//			if (pest_scenario.get_pestpp_options().get_ies_use_empirical_prior()) {
//                //the second return matrix should be shrunk optimally to be nonsingular...but who knows!
//                Covariance shrunk_cov = dv.get_empirical_cov_matrices(&file_manager).second;
//                shrunk_cov.inv_ip();
//                parcov_inv = shrunk_cov.e_ptr()->toDense();
//            }
//			else
//            {
//			    parcov_inv = parcov.inv().e_ptr()->toDense();
//            }
			//cout << "parcov inv: " << endl << parcov_inv << endl;
			//TODO: Matt to check consistency being sample cov forms
            //message(1, "empirical parcov inv:", parcov_inv);  // tmp
			
			// try pseudo_inv_ip()
			//Covariance x;
			//x = Covariance(dv_names, dv_cov_matrix.sparseView(), Covariance::MatType::SPARSE);
			//x.pseudo_inv_ip(pest_scenario.get_svd_info().eigthresh, pest_scenario.get_svd_info().maxsing);
			//message(1, "pseudo inv:", x);  // tmp

			// CMA implementation to go here

			// compute dec var-phi cross-cov vector
			// see eq (9) of Dehdari and Oliver 2012 SPE and Fonseca et al 2015 SPE
			// start by computing mean-shifted obj function ensemble
            Eigen::MatrixXd ivec, upgrade_1, s, V, U, st;
            SVD_REDSVD rsvd;
            //SVD_EIGEN rsvd;
            rsvd.set_performance_log(performance_log);
            //dv_anoms.transposeInPlace();

            rsvd.solve_ip(dv_anoms, s, U, V, pest_scenario.get_svd_info().eigthresh, pest_scenario.get_svd_info().maxsing);
            Eigen::MatrixXd dv_anoms_pseudoinv = V * s.asDiagonal().inverse() * U.transpose();
            Eigen::MatrixXd obj_anoms(dv.shape().first,1);
            if (use_obj_obs) {
                obj_anoms = oe.get_eigen_anomalies(vector<string>(), vector<string>{obj_func_str},center_on);
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
                    for (auto& dv :dv_names)
                    {
                        oval += obj_func_coef_map.at(dv) * real(vmap.at(dv));
                    }
                    obj_anoms(i,0) = oval;
                    i++;
                }
                obj_anoms.array() -= obj_anoms.mean();

            }
			Eigen::MatrixXd cross_cov_vector;  // or Eigen::VectorXd?
			//cross_cov_vector = 1.0 / (dv.shape().first - 1.0) * (dv_anoms.transpose() * obj_anoms);
			//cout << "dv-obj_cross_cov:" << endl << cross_cov_vector << endl;
			
			// now compute grad vector
			// this is a matrix-vector product; the matrix being the pseudo inv of diag empirical dec var cov matrix and the vector being the dec var-phi cross-cov vector\
			// see, e.g., Chen et al. (2009) and Fonseca et al. (2015) 
			//grad = parcov_inv * cross_cov_vector;//*parcov.e_ptr() * cross_cov_vector;
			grad = 1.0 / ((double)dv.shape().first - 1.0) * (dv_anoms_pseudoinv * obj_anoms);
			// if (constraints)
			//{
			//	ss.str("");
			//	ss << "compute ensemble approx to (active) constraint jacobian";
			//	string s = ss.str();
			//	message(1, s);
			//	throw_sqp_error("TODO");
			//}

             //throw_sqp_error("obs-based obj for ensembles not implemented");

             //todo: localize the gradient here - fun times

		}
		//pi base obj, need representative dv values using the "center_on" arg
		//represent the mean/median/base - that is, derived from the "center_on" arg
		//todo: for now, just using mean dv values
//		else
//		{
//			//if not center_on arg, use the mean dv values
//			//if (center_on.size() == 0)
//			//{
//			//	//pair<map<string, double>, map<string, double>> mm = dv.get_moment_maps();
//			//	for (int i = 0; i < dv_names.size(); i++)
//			//	{
//			//		grad[i] = obj_func_coef_map[dv_names[i]];// * mm.first[dv_names[i]];
//			//	}
//			//}
//			//else
//			//{
//
//			//	grad = dv.get_real_vector(pest_scenario.get_pestpp_options().get_ies_center_on());
//			//}
//			//
//			//I think we should just eval the gradient around the current dv values
//			grad = calc_gradient_vector_from_coeffs(_current_dv_values);
//		}
			
	}
	else
	{
		//obs-based obj
		if (use_obj_obs)
		{
			//just a jco row
			vector<string> obj_name_vec{ obj_func_str };
			Eigen::MatrixXd t = jco.get_matrix(obj_name_vec, dv_names);
			grad = t.row(0);
		}
		//pi based obj
		else
		{
			grad = calc_gradient_vector_from_coeffs(_current_dv_values);
		}
	}
	Parameters pgrad = _current_dv_values;
	pgrad.update_without_clear(dv_names, grad);
	return pgrad;
}

//Eigen::MatrixXd SeqQuadProgram::update_hessian()
//{
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
//}

//pair<Eigen::VectorXd, Eigen::VectorXd> SeqQuadProgram::_solve_eqp()
//{
//	//return pair<>
//}


pair<Eigen::VectorXd, Eigen::VectorXd> SeqQuadProgram::_kkt_null_space(Eigen::MatrixXd& G, Eigen::MatrixXd& constraint_jco, Eigen::VectorXd& constraint_diff, Eigen::VectorXd& curved_grad, vector<string>& cnames)
{

	Eigen::VectorXd search_d, lm;
	
	// check: A full rank
	// check: reduced hessian ZTGZ is non pos def

	// compute orthog bases of A
	// prob put the beblow in a standalone function too
	// ------
	bool use_qr = false;  // unsure if needed with randomized SVD option
	Eigen::MatrixXd Z;
	// check A has more rows than cols // this should have been caught before this point
	//if (use_qr)
	//{
	//	throw_sqp_error("QR decomposition for orthog basis matrix computation not implemented");
	//}
	//else
	//{
	// rSVD
	Eigen::VectorXd x;
	Eigen::MatrixXd V, U, S_, s;
	SVD_REDSVD rsvd;
	//SVD_EIGEN rsvd;
	message(1, "using randomized SVD to compute basis matrices of constraint JCO for null space KKT solve", constraint_jco);
	rsvd.set_performance_log(performance_log);

	message(1, "A before", constraint_jco);
	rsvd.solve_ip(constraint_jco, s, U, V, pest_scenario.get_svd_info().eigthresh, pest_scenario.get_svd_info().maxsing);
	message(1, "A after", constraint_jco);
	S_ = s.asDiagonal();  // check truncation done automatically
	message(1, "singular values of A matrix", S_);  // tmp
	message(1, "V", V);  // tmp
	//Z = V.transpose().conjugate();
	Z = V.conjugate();  // tmp
	message(1, "Z", Z);  // tmp

	Eigen::MatrixXd Y = constraint_jco.transpose();
	message(1, "Y", Y);  // tmp
	//}
	// ------

	// solve p_range_space
	Eigen::VectorXd p_y, rhs;
	Eigen::MatrixXd coeff;
	coeff = constraint_jco * Y;
	message(1, "coeff matrix for p_range_space component", coeff);  // tmp
	rhs = (-1. * constraint_diff);
	message(1, "rhs", rhs);
	cout << "starting SVD..." << endl;
	rsvd.solve_ip(coeff, s, U, V, pest_scenario.get_svd_info().eigthresh, pest_scenario.get_svd_info().maxsing);
	cout << "...done..." << endl;
	S_ = s.asDiagonal().inverse();
	cout << "s:" << S_.rows() << "," << S_.cols() << ", U:" << U.rows() << "," << U.cols() << ", V:" << V.rows() << "," << V.cols() << endl;
	cout << endl << endl;

	coeff = V * S_ * U.transpose();

	message(1, "coeff inv", coeff);
	p_y = coeff * rhs;
	message(1, "p_y", p_y);  // tmp
	// todo assert here for p_y != 0 if sum_viol == 0 (this is the component that rectifies constraint viol)

	// solve p_null_space
	bool simplified_null_space_approach = false; // see Nocedal and Wright (2006) pg 538-9
	message(1, "hess", G);
	Eigen::MatrixXd red_hess = G;
	red_hess = Z.transpose() * red_hess * Z;
	message(1, "red hess", red_hess);
	bool cholesky = false;
	Eigen::VectorXd p_z;
	//try following, else ZTGZ is not pos-def
	//if LBFGS
	if (Z.size() > 0)
	{
		// check use of grads or curved grads below (want second order)?
		if (simplified_null_space_approach)
		{
			message(1, "using simplified approach in KKT null space solve...");
			//rhs = -1. * Z.transpose() * curved_grad;
			// simplify by removing cross term (or ``partial hessian'') matrix (zTgy), which is approp when approximating hessian (zTgz) (as p_y goes to zero faster than p_z)
			if (cholesky)
			{
				// todo throw_sqp_error("cholesky decomp for null space KKT solve not implemented");
			}
			else
			{
				//throw_sqp_error("todo");
				//p_z = solve: red_hess, rhs
			}
		}
		else
		{
			message(1, "carrying out KKT null space solve...");
			Eigen::MatrixXd X1, coeff;
			Eigen::VectorXd rhs;
			coeff = red_hess;
			X1 = Z.transpose() * G * Y;
			message(1, "larger cross-term matrix", X1);
			rhs = (-1. * X1 * p_y) - (Z.transpose() * curved_grad);
			if (cholesky)
			{
				
				//throw_sqp_error("cholesky decomp for null space KKT solve not implemented");  // todo JDub I noticed there is an built in cholesky decomposition? 
				//l = cholesky(red_hess);
				//rhs2 = solve: l, rhs;
				//p_z = solve: l.transpose(), rhs2;
			    //jwhite: I was tinkering with adding cholesky as a cov matrix method but it old as!
				//jwhite: I just hacked this in - not tested!
				Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
				Eigen::VectorXd rhs2 = solver.compute(red_hess.sparseView()).solve(rhs);
				p_z = solver.solve(rhs2);
				
			}
			else
			{
				// try straight inverse here
				// todo rSVD here instead?
				p_z = coeff.inverse() * rhs;
				message(1, "p_z", p_z);  // tmp

			}
		}
	}
	else
	{
		throw_sqp_error("zero null space dimensionality wrt active constraints");  // tmp
	}
	//else not sure can be done

	// combine to make total direction
	message(1, "combining range and null space components of search direction");  // tmp
	search_d = Y * p_y + Z * p_z;
    //message(1, "SD", search_d);  // tmp

	if (search_d.size() != curved_grad.size())
	{
		throw_sqp_error("search direction vector computation error (in null space KKT solve method)!");
	}

	
	// compute lagrangian multipliers
	//if LBFGS
	message(1, "computing lagrangian multipliers...");  // tmp
	if (simplified_null_space_approach)
	{
		// simplify by dropping dependency of lm on hess (considered appropr given p converges to zero whereas grad does not..; Nocedal and Wright pg. 539)
		// throw_sqp_error("reduced hessian in KKT null space solve not implemented");
		// todo
		// lm = solve: constraint_jco * constraint_jco.transpose(), constraint_jco * curved_grad
	}
	else
	{
		// Nocedal and Wright pg. 457 and 538
		rhs = Y.transpose() * curved_grad;
		coeff = (constraint_jco * Y).transpose();
		lm = coeff.inverse() * rhs;
		message(1, "lm", lm);  // tmp

	}
	//else not sure can be done


	return pair<Eigen::VectorXd, Eigen::VectorXd>(search_d, lm);
}


pair<Eigen::VectorXd, Eigen::VectorXd> SeqQuadProgram::_kkt_direct(Eigen::MatrixXd& G, Eigen::MatrixXd& constraint_jco, Eigen::VectorXd& constraint_diff, Eigen::VectorXd& curved_grad, vector<string>& cnames)
{
	
	//check A full rank

	// forming system to be solved - this is filth but it works..
	Eigen::MatrixXd coeff_u(dv_names.size(), dv_names.size() + cnames.size());  // todo only in WS
	coeff_u << G, constraint_jco.transpose();
	Eigen::MatrixXd coeff_l(cnames.size(), dv_names.size() + cnames.size());  // todo only in WS
	coeff_l << constraint_jco, Eigen::MatrixXd::Zero(cnames.size(), cnames.size());
	Eigen::MatrixXd coeff(dv_names.size() + cnames.size(), dv_names.size() + cnames.size());  // todo only in WS
	coeff << coeff_u, coeff_l;
	message(1, "coeff", coeff);  // tmp

	Eigen::VectorXd rhs(curved_grad.size() + constraint_diff.size());
	rhs << curved_grad, constraint_diff;  // << vec1, vec2;
	message(1, "rhs", rhs);  // tmp

	Eigen::VectorXd x;
	Eigen::MatrixXd V, U, S_, s;
	SVD_REDSVD rsvd;
	//SVD_EIGEN rsvd;
	rsvd.set_performance_log(performance_log);

	rsvd.solve_ip(coeff, s, U, V, pest_scenario.get_svd_info().eigthresh, pest_scenario.get_svd_info().maxsing);
	S_ = s.asDiagonal();
	message(1, "singular values of KKT matrix", S_);  // tmp

	// an old friend!
	x = V * S_.inverse() * U.transpose() * rhs;
	message(1, "solution vector of steps and lagrange mults", x);  // tmp

	Eigen::VectorXd search_d, lm;
	search_d = x.head(dv_names.size());  // add rigor here or at least asserts to ensure operating on correct elements
	lm = x.tail(x.size() - dv_names.size());  // add rigor here or at least asserts to ensure operating on correct elements

	
	return pair<Eigen::VectorXd, Eigen::VectorXd> (search_d, lm);
}


pair<Eigen::VectorXd, Eigen::VectorXd> SeqQuadProgram::calc_search_direction_vector(const Parameters& _current_dv_values, Eigen::VectorXd& grad_vector)
{

	Eigen::VectorXd search_d, lm;
	pair<Eigen::VectorXd, Eigen::VectorXd> x;
	

	//message(1, "hessian:", hessian);  // tmp
    Mat constraint_mat;
    if (use_ensemble_grad) {
        message(1, "getting ensemble-based working set constraint matrix");
        constraint_mat = constraints.get_working_set_constraint_matrix(current_ctl_dv_values, current_obs, dv, oe,true);
    }
    else
    {
        message(2, "getting working set constraint matrix");
        constraint_mat = constraints.get_working_set_constraint_matrix(current_ctl_dv_values, current_obs, jco, true);
    }
	//todo:probably need to check if constraint_mat has any nonzeros?
	vector<string> cnames = constraint_mat.get_row_names();

	if (cnames.size() > 0)  // solve constrained QP subproblem
	{
		// solve (E)QP sub-problem via active set method for search direction

		// todo put in _solve_eqp() func - self.search_d, self.lagrang_mults = self._solve_eqp(qp_solve_method=self.qp_solve_method)

		// Direct QP (KKT system) solve method herein; assumes convexity (pos-def-ness) and that A has full row rank.
		// Direct method here is just for demon purposes - will require the other methods offered here given that the KKT system will be indefinite.
		// Refer to (16.5, 18.9) of Nocedal and Wright (2006)

		// collate the pieces

		// constraint jco
		
		//todo: something better than just echoing working set to screen and rec...
		message(0, "current working set:", cnames);
		Eigen::MatrixXd constraint_jco = constraint_mat.e_ptr()->toDense();  // or would you pref to slice and dice each time - this won't get too big but want to avoid replicates  // and/or make Jacobian obj?

		//constraint_jco = jco.get_matrix(cnames, dv_names);
		message(1, "A:", constraint_jco);  // tmp
		// add check here that A is full rank; warn that linearly dependent will be removed via factorization

		// constraint diff (h = Ax - b)
		// todo only for constraints in WS only
		Eigen::VectorXd Ax, b;
		Eigen::VectorXd constraint_diff(cnames.size());
		
		pair<Eigen::VectorXd, Eigen::VectorXd> p = constraints.get_obs_resid_constraint_vectors(current_ctl_dv_values, current_obs, cnames);
		b = p.first;
		constraint_diff = p.second;
		message(1, "constraint diff:", constraint_diff);  // tmp
		
		// throw error here if not all on/near constraint
		if ((constraint_diff.array().abs() > filter.get_viol_tol()).any())  // todo make some level of forgiveness with a tolerance parameter here
		{
			//throw_sqp_error("not on constraint");  // better to pick this up elsewhere (before) anyway
			cout << "constraint diff vector: " << constraint_diff.array() << endl;
			message(0,"WARNING: not on constraint but working set not empty, continuing...");
		}

		// some transforms for solve
		Eigen::MatrixXd G; // or //Eigen::SparseMatrix<double> G;?
		//will hessian ever be singular?  If not, then just use inv method...

		G = *hessian.inv().e_ptr() * 2.;  // TODO: double check and give ref

		Eigen::VectorXd c;
		c = grad_vector + G * _current_dv_values.get_data_eigen_vec(dv_names);  // TODO: check not just grad (see both) and check sign too...


		string eqp_solve_method; // probably too heavy to be a ++arg
		eqp_solve_method = "null_space";
		if (eqp_solve_method == "null_space")
		{
			x = _kkt_null_space(G, constraint_jco, constraint_diff, c, cnames);
			search_d = x.first;
			lm = x.second;
			message(1, "sd:", search_d);  // tmp
			message(1, "lm:", lm);  // tmp
		}
		else if (eqp_solve_method == "direct")
		{
			x = _kkt_direct(G, constraint_jco, constraint_diff, c, cnames);
			search_d = x.first;
			lm = x.second;
			message(1, "sd:", search_d);  // tmp
			message(1, "lm:", lm);  // tmp

		}
		else // if "schur", "cg", ...
		{
			throw_sqp_error("eqp_solve_method not implemented");
		}

		//todo check sign of search_d; think the way the system is formulated, it solves for -dx


		// put the following in constraints.reduce_working_set(bool& unsuccessful) func
		// call here constraints.reduce_working_set(false);
		//----------
		// see alg (16.3) of Nocedal and Wright (2006)

		message(1, "current working set:", cnames);  // tmp
		message(1, "lagrangian multipliers:", lm);  // tmp

		// check whether to `stop and drop` a single constraint from working set due to p ~= 0 or convergence in p (or unsuccessful iteration)
		// note sign of search_d doesn't matter here
		bool converged = false;  // todo approx step convergence test here? e.g. current search_d is within 5% of last itn and WS is same as last itn
		bool unsuccessful = false;  // todo JDub pass unsuccessful arg to constraints.reduce_working_set() if accept is false 
		// can also check convergence with orthogonality of p_y / p_z > 1.0 if using null_space KKT solve method.. later

		bool search_d_approx_zero = false;
		double tol = 0.01; //todo: move this to a ++arg
		for (int i=0;i<search_d.size();i++)
			if (abs(search_d[i]) < tol)
			{
				search_d_approx_zero = true;
				break;
			}
		
		//search_d_approx_zero = true;  // tmp

		if (search_d_approx_zero | converged | unsuccessful)
		{
			vector<string> working_set_ineq_names = constraints.get_working_set_ineq_names(cnames);  // lm sign only iterpret-able for ineq constraints in working set
			message(1, "ineq constraints in working set:", working_set_ineq_names);  // tmp

			double lm_min = 0.0;
			int idx;
			for (int i = 0; i < cnames.size(); i++)
			{
				// todo JDub - need another if here such that we get the idx and lm_min for ineq constraints only
				//jwhite: this is slow if heaps (1000+) in working_set_ineq_names - can swap to set<string> to be faster
				if (find(working_set_ineq_names.begin(), working_set_ineq_names.end(), cnames[i]) == working_set_ineq_names.end())
					continue;
				if ((lm[i] < 0.0) && (lm[i] < lm_min))
				{
					idx = i;
					lm_min = lm[i];
				}
			}

			if (lm_min == 0.0)
			{
				throw_sqp_error("optimal soln detected at solve EQP step (lagrangian multiplier for all ineq constraints in working set is non-neg)...");
			}
			message(1, "lm_min:", lm_min);  // tmp
			message(1, "lm_min_idx:", cnames[idx]);  // tmp

			string to_drop = cnames[idx];
			cnames.erase(cnames.begin() + idx);
			message(1, "cnames after dropping attempting:", cnames);  // tmp

			// todo JDub now we need to skip alpha testing for this iter and recalc search_d without this constraint... i moved to next iter in proto because it was cheap, probably want to go back to start of search_d computation?
			// jwhite: will this process happen multiple times?  If so, we probably need to wrap this process in a "while true" loop
			// mjk: as we discussed, the while true would be awesome!
			//}
		}
		else  // progress with this iter
		{
			// todo determine clever (base) alpha
			double alpha = 1.;
		}

		// return alpha, go_to
		//----------


	}
	else  // solve unconstrained QP subproblem
	{
		message(1, "constraint working set is empty, problem is currently unconstrained...");
		search_d = *hessian.e_ptr() * grad_vector;


//		if (LBFGS)
//		ss.str("");
//		ss << "use LBFGS algorithm to solve for search direction in unconstrained case";
//		string s = ss.str();
//		message(1, s);
//		//self.search_d = self._LBFGS_hess_update(memory = memory, self_scale = hess_self_scaling)
//		throw_sqp_error("TODO");
	}

	// sign of direction
	if (obj_sense == "minimize")
		search_d *= -1.;

	message(1, "sd:", search_d.transpose());  // tmp

	return pair<Eigen::VectorXd, Eigen::VectorXd> (search_d, lm);  // lm will be empty if non-constrained solve
}

Eigen::VectorXd SeqQuadProgram::fancy_solve_routine(const Parameters& _current_dv_num_values, const Parameters& _grad_vector)
{

	// grad vector computation
	Eigen::VectorXd grad = _grad_vector.get_data_eigen_vec(dv_names);

	// search direction computation
	//Eigen::VectorXd search_d = calc_search_direction_vector(_current_dv_num_values, grad);  
	pair<Eigen::VectorXd, Eigen::VectorXd> x = calc_search_direction_vector(_current_dv_num_values, grad);

	// todo undertake search direction-related tests, e.g., point down-hill
	// and check if constraints in working set cause zero search_d (and go to next iteration if so)
	//	if (constraints is True);  // irrespective of shape of working set
	//	{
	//		ss.str("");
	//		ss << "check if constraints in working set cause zero search_d";
	//		string s = ss.str();
	//		message(1, s);
	//		//alpha, next_it = self._active_set_method(first_pass = True)
	//		throw_sqp_error("TODO");
	//	}


	return x.first;  // search_d;
}

bool SeqQuadProgram::solve_new()
{
	stringstream ss;
	ofstream& frec = file_manager.rec_ofstream();
	if ((use_ensemble_grad) && (dv.shape().first <= error_min_reals))
	{
		message(0, "too few active realizations:", oe.shape().first);
		message(1, "need more than ", error_min_reals);
		throw_sqp_error(string("too few active realizations, cannot continue"));
	}
	else if ((use_ensemble_grad) && (dv.shape().first < warn_min_reals))
	{
		ss.str("");
		ss << "WARNING: less than " << warn_min_reals << " active realizations...might not be enough";
		string s = ss.str();
		message(1, s);
	}

	Parameters _current_num_dv_values = current_ctl_dv_values;  // make copy
	ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();
	pts.ctl2numeric_ip(_current_num_dv_values);  // convert to numeric format
	_current_num_dv_values = _current_num_dv_values.get_subset(dv_names.begin(), dv_names.end());  // just dec var operation
	
	if (use_ensemble_grad)
	{
		dv.transform_ip(ParameterEnsemble::transStatus::NUM);
		vector<double> mean_vec = dv.get_mean_stl_var_vector();
		_current_num_dv_values.update(dv_names, mean_vec);
	}

	//Parameters dv_num_candidate;
	ParameterEnsemble dv_candidates(&pest_scenario,&rand_gen);
	dv_candidates.set_trans_status(ParameterEnsemble::transStatus::NUM);
	
	
	// search direction computation
	Eigen::VectorXd search_d;
	search_d = fancy_solve_routine(_current_num_dv_values,current_grad_vector);

	
	//backtracking/line search factors
	//TODO: make this a ++ arg or tunable or something clever
		//if (constraints);
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
	
	vector<string> real_names;
    vector<double> scale_vals;
	scale_vals.clear();
	for (auto& sf : pest_scenario.get_pestpp_options().get_sqp_scale_facs())
    {
	    scale_vals.push_back(sf * BASE_SCALE_FACTOR);
    }

    if ((use_ensemble_grad) && (SOLVE_EACH_REAL))
    {
        for (auto sv : scale_vals)
        {
            for (auto& real_name : dv.get_real_names())
             {
                ss.str("");
                ss << "dv_cand_" << real_name << "_sv:" << left << setw(8) << setprecision(3) << sv;
                real_names.push_back(ss.str());
            }
        }

    }
    else {
        for (auto sv : scale_vals) {
            ss.str("");
            ss << "dv_cand_sv:" << left << setw(8) << setprecision(3) << sv;
            real_names.push_back(ss.str());
        }
    }
	dv_candidates.reserve(real_names, dv_names);
    int ii = 0;
    vector<double> used_scale_vals;
    map<string,double> real_sf_map;
	for (int i=0;i<scale_vals.size();i++)
	{
		double scale_val = scale_vals[i];
		ss.str("");
		ss << "starting calcs for scaling factor" << scale_val;
		message(1, "starting lambda calcs for scaling factor", scale_val);
		message(2, "see .log file for more details");

		if ((use_ensemble_grad) && (SOLVE_EACH_REAL))
        {
		    Parameters real_grad,num_candidate=current_ctl_dv_values;
		    pest_scenario.get_base_par_tran_seq().ctl2numeric_ip(num_candidate);
		    Eigen::VectorXd real,scale_search_d,cvals;
		    dv.transform_ip(ParameterEnsemble::transStatus::NUM);

		    for (auto& real_name : dv.get_real_names())
            {
		        real_grad = calc_gradient_vector(current_ctl_dv_values,real_name);

                real = dv.get_real_vector(real_name);
                num_candidate.update_without_clear(dv.get_var_names(),real);
                search_d = fancy_solve_routine(num_candidate,real_grad);
                scale_search_d = search_d * scale_val;
                cvals = num_candidate.get_data_eigen_vec(dv_names);
                cvals.array() += scale_search_d.array();
                //num_candidate.update_without_clear(dv_names, cvals);
                //Eigen::VectorXd vec = num_candidate.get_data_eigen_vec(dv_names);
                dv_candidates.update_real_ip(real_names[ii], cvals);
                real_sf_map[real_names[ii]] = scale_val;
                ii++;
                used_scale_vals.push_back(scale_val);


            }
        }
		else {
            //dv_num_candidate = fancy_solve_routine(scale_val, _current_num_dv_values);
            Parameters num_candidate = _current_num_dv_values;

            Eigen::VectorXd scale_search_d = search_d * scale_val;
            if (scale_search_d.squaredNorm() < 1.0 - 10)
                message(1, "very short upgrade for scale value", scale_val);

            Eigen::VectorXd cvals = num_candidate.get_data_eigen_vec(dv_names);

            cvals.array() += scale_search_d.array();
            num_candidate.update_without_clear(dv_names, cvals);

            //Eigen::VectorXd vec = dv_num_candidate.get_data_eigen_vec(dv_names);
            Eigen::VectorXd vec = num_candidate.get_data_eigen_vec(dv_names);
            dv_candidates.update_real_ip(real_names[i], vec);
            used_scale_vals.push_back(scale_val);
            real_sf_map[real_names[i]] = scale_val;
        }

		ss.str("");
		message(1, "finished calcs for scaling factor:", scale_val);

	}


//// check that we need to have the gradient at candidate, e.g., for Hessian purposes
// if (use_ensemble_grad)
//	{
//		ss.str("");
//		ss << "draw dec var en at candidate mean vector";
//		string s = ss.str();
//		message(1, s);
//		throw_sqp_error("TODO");

//		ss.str("");
//		ss << "evaluate ensembles at trial alphas";
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
	//			// can also do some reuse here
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

	//enforce bounds on candidates - TODO: report the shrinkage summary that enforce_bounds returns
	dv_candidates.enforce_bounds(performance_log,false);
	ss.str("");
	ss << file_manager.get_base_filename() << "." << iter << ".dv_candidates.csv";
	dv_candidates.to_csv(ss.str());

	//check for duplicate candidates
	Eigen::VectorXd v1,v2;
	double d;
	vector<string> drop;
	set<int> jvals;
	for (int i=0;i<dv_candidates.shape().first;i++)
    {
	    v1 = dv_candidates.get_real_vector(i);
	    for (int j=i+1;j<dv_candidates.shape().first;j++) {
            v2 = (dv_candidates.get_real_vector(j) - v1).array() / v1.array().cwiseAbs();
            d = v2.transpose() * v2;
            if ((abs(d) < 1e-7) && (jvals.find(j) == jvals.end())) {
                message(1, "duplicate candidates:", vector<string>{real_names[i], real_names[j]});
                drop.push_back(real_names[j]);
                jvals.emplace(j);
            }
        }

    }
	if (drop.size() > 0)
    {
	    message(1,"dropping the following duplicate candidates: ",drop);
	    dv_candidates.drop_rows(drop,true);
	    used_scale_vals.clear();
	    for (auto& real_name : dv_candidates.get_real_names())
        {
	        used_scale_vals.push_back(real_sf_map.at(real_name));
        }

    }


	message(0, "running candidate decision variable batch");
	vector<double> passed_scale_vals = scale_vals;	
	//passed_scale_vals will get amended in this function based on run fails...

	ObservationEnsemble oe_candidates = run_candidate_ensemble(dv_candidates);
    ss.str("");
    ss << file_manager.get_base_filename() << "." << iter << ".oe_candidates.csv";
    oe_candidates.to_csv(ss.str());
	//todo: decide which if any dv candidate to accept...
	map<string,double> sf_map;
	for (int i=0;i<real_names.size();i++)
    {
	    sf_map[real_names[i]] = used_scale_vals[i];
    }
	bool success = pick_candidate_and_update_current(dv_candidates, oe_candidates,sf_map);
	if (!success)
	{
		//// deal with unsuccessful iteration

		// call constraints.reduce_working_set(true) here
		
		return false;
		
	}
	
	return true;  // reporting and saving done next
}

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
	map<string,double> viol_map = constraints.get_unsatified_obs_constraints(current_obs,filter.get_viol_tol());
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
    ies_pest_scenario.get_control_info_4_mod().noptmax = 3; //TODO: make this an option some how?
    ss.str("");
    string org_base = file_manager.get_base_filename();
    ss << "feas_ies_" << iter << "_" << org_base;
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
	//what to do here? maybe we need to eval the kkt conditions to pick a new point that maintains the hessian?
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
	    //cout << "real:" << oreal_names[i] << ", phi: " << aphi_map[oreal_names[i]] << endl;
		//Eigen::VectorXd real = ies_pe_ptr->get_eigen_ptr()->row(i);
		//Eigen::VectorXd d = real - cdv;
		//double diff = (d.array() * d.array()).sum();
		//if (diff < mndiff)
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
	//update current obs
	cdv = ies.get_oe().get_real_vector(mndiff_idx);
	names = ies.get_oe().get_var_names();
	current_obs.update(names, cdv);
	constraints.sqp_report(iter, current_ctl_dv_values, current_obs, true, "post feasible seek");
	//todo: probably more algorithmic things here...
	last_best = get_obj_value(current_ctl_dv_values, current_obs);
	last_viol = constraints.get_sum_of_violations(current_ctl_dv_values, current_obs);
	best_phis[best_phis.size()-1] = last_best;
	best_violations[best_violations.size() -1] = last_viol;
	message(1, "finished seeking feasible, reset best phi,infeasible value to ", vector<double>{last_best,last_viol});
	return false;
}


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

map<string, double> SeqQuadProgram::get_obj_map(ParameterEnsemble& _dv, ObservationEnsemble& _oe)
{
	Eigen::VectorXd obj_vec = get_obj_vector(_dv, _oe);
	vector<string> real_names = _dv.get_real_names();
	map<string, double> obj_map;
	for (int i = 0; i < real_names.size(); i++)
		obj_map[real_names[i]] = obj_vec[i];

	return obj_map;


}

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
			//pts.ctl2numeric_ip(pars);
			real = _dv.get_real_vector(i);
			pars.update_without_clear(vnames, real);
			pts.numeric2ctl_ip(pars);
			v = get_obj_value(pars, current_obs); //shouldn't be using current obs since this is dv-based obj
			obj_vec[i] = v;
			pts.ctl2numeric_ip(pars);
		}
	}
	return obj_vec;
}

bool SeqQuadProgram::pick_candidate_and_update_current(ParameterEnsemble& dv_candidates, ObservationEnsemble& _oe, map<string,double>& sf_map)
{
	//decide!
	message(0, " current best phi:", last_best);
	stringstream ss;
	Eigen::VectorXd obj_vec = get_obj_vector(dv_candidates, _oe);
	double oext,oviol;
	if (obj_sense == "minimize")
		oext = numeric_limits<double>::max();
	else
		oext = numeric_limits<double>::min();
	int idx = -1;
	vector<string> real_names = dv_candidates.get_real_names();
	map<string, double> obj_map = get_obj_map(dv_candidates, _oe);
	//todo make sure chances have been applied before now...
	map<string, map<string, double>> violations = constraints.get_ensemble_violations_map(dv_candidates,_oe,filter.get_viol_tol(),true);
	Parameters cand_dv_values = current_ctl_dv_values;
	Observations cand_obs_values = current_obs;
	Eigen::VectorXd t;
	vector<string> onames = _oe.get_var_names();
	ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();
	bool filter_accept;
	string tag;
	vector<double> infeas_vec;
	vector<int> accept_idxs;
	bool accept;
	for (int i = 0; i < obj_vec.size(); i++)
	{
		ss.str("");
		//ss << "scale factor " << setprecision(4) << scale_vals[i];
		ss << "candidate: " << real_names[i] << ", scale factor: " << sf_map.at(real_names[i]);
		tag = ss.str();
		ss.str("");
		ss << "candidate: " << tag << " phi: " << obj_vec[i];
		double infeas_sum = 0.0;
		for (auto& v : violations[real_names[i]])
		{
			infeas_sum += v.second;
		}
		ss << " infeasibilty total: " << infeas_sum << ", ";
		infeas_vec.push_back(infeas_sum);
		filter_accept = filter.accept(obj_vec[i], infeas_sum,iter,sf_map.at(real_names[i]),false);
		if (filter_accept)
			ss << " filter accepted ";
		else
			ss << " filter rejected ";
		message(1, ss.str());
        if (filter_accept)
        {
            if ((best_violation_yet > 1e-7) && (infeas_vec[i] > (best_violation_yet * 2.0)))
            {
                ss << ", infeasibility exceeds previous best infeasibility threshold";
            }
            else{
                accept_idxs.push_back(i);
            }
            accept_idxs.push_back(i);
//            idx = i;
//            oext = obj_vec[i];
//            // now update the filter recs to remove any dominated pairs
//            filter.update(obj_vec[i], infeas_sum, iter, alpha_vals[i]);
//            accept = true;
        }
		//report the constraint info for this candidate
		t = dv_candidates.get_real_vector(real_names[i]);
		cand_dv_values = current_ctl_dv_values;
		cand_dv_values.update_without_clear(dv_names, t);
		pts.numeric2ctl_ip(cand_dv_values);
		t = _oe.get_real_vector(real_names[i]);
		cand_obs_values.update_without_clear(onames, t);
		constraints.sqp_report(iter, cand_dv_values, cand_obs_values, false,tag);
	}

    if (accept_idxs.size() > 0)
    {
        accept = true;
        //since all of the these passed the filter, choose the one with lowest phi
        double min_obj = 1.0e+300;
        ss.str("");
        ss << "number of scale factors passing filter:" << accept_idxs.size();
        message(1,ss.str());
        for (auto iidx : accept_idxs)
        {
            if (obj_vec[iidx] < min_obj)
            {
                min_obj = obj_vec[iidx];
                idx = iidx;
                oext = obj_vec[iidx];
                oviol = infeas_vec[iidx];
            }
        }
        filter.update(min_obj,infeas_vec[idx],iter,sf_map.at(real_names.at(idx)));
    }
	else
    {
	    message(0,"filter failed, checking for feasible solutions....");
	    double viol_tol = filter.get_viol_tol();
        for (int i=0;i<infeas_vec.size();i++)
        {
            if (infeas_vec[i] < viol_tol)
            {
                if ((obj_sense=="minimize") && (obj_vec[i] < oext))
                {
                    idx = i;
                    oext = obj_vec[i];
                    oviol = infeas_vec[i];
                }
                else if ((obj_sense == "maximize") && (obj_vec[i] > oext))
                {
                    idx = i;
                    oext = obj_vec[i];
                    oviol = infeas_vec[i];
                }
            }
        }
        if (idx == -1)
        {
            message(0,"no feasible solutions, choosing lowest constraint violation...");
            //now what?
            //how about least violation - probably going to hand off to ies now anyway...
            double viol_min = 1e+300;
            for (int i=0;i<infeas_vec.size();i++)
            {
                if (infeas_vec[i] < viol_min)
                {
                    viol_min = infeas_vec[i];
                    idx = i;
                    oext = obj_vec[i];
                    oviol = infeas_vec[i];
                }
            }
        }

    }
	if (idx == -1)
	    throw_sqp_error("shits busted");
    message(0, "best phi and infeas this iteration: ", vector<double>{oext,oviol});
    t = dv_candidates.get_real_vector(real_names[idx]);
    cand_dv_values = current_ctl_dv_values;
    cand_dv_values.update_without_clear(dv_names, t);
    pts.numeric2ctl_ip(cand_dv_values);
    t = _oe.get_real_vector(real_names[idx]);
    cand_obs_values.update_without_clear(onames, t);
    ss.str("");
    ss << "best candidate (scale factor: " << setprecision(4) << sf_map.at(real_names[idx]) << ", phi: " << oext << ", infeas: " << oviol << ")";
    constraints.sqp_report(iter, cand_dv_values, cand_obs_values, true,ss.str());
    filter.report(file_manager.rec_ofstream(),iter);
	//TODO: need more thinking here - do we accept only if filter accepts?  I think so....
	//if (((obj_sense == "minimize") && (oext < last_best)) || ((obj_sense == "maximize") && (oext > last_best)))
	if (accept)
	{
		//todo:update current_dv and current_obs
		message(0, "accepting upgrade", real_names[idx]);
		t = dv_candidates.get_real_vector(idx);
		vector<string> vnames = dv_candidates.get_var_names();
		Parameters p;
		p.update_without_clear(vnames, t);
		
		pts.numeric2ctl_ip(p);
		for (auto& d : dv_names)
			current_ctl_dv_values[d] = p[d];
		t = _oe.get_real_vector(idx);
		current_obs.update_without_clear(onames, t);
		last_best = oext;
		last_viol = oviol;
		message(0, "new best phi and infeas:", vector<double>{last_best,last_viol});
        best_phis.push_back(oext);
        best_violations.push_back(oviol);

		// todo add constraint (largest violating constraint not already in working set) to working set
		// is this the right place to do this? after accepting a particular candidate? 
		// also can we adapt alpha_mult based on subset? using concept of blocking constraint here?
		// take diff between vector of strings of constraints in working set and constraints with non-zero violation (return constraint idx from filter?)

		//if no filter-accepted solutions and we are in violation...
		if (infeas_vec[idx] > filter.get_viol_tol())
        {
		    n_consec_infeas++;
        }
		else {
            if (use_ensemble_grad) {
                double new_par_sigma = pest_scenario.get_pestpp_options().get_par_sigma_range();
                new_par_sigma = new_par_sigma * (PAR_SIGMA_INC_FAC);
                new_par_sigma = min(new_par_sigma, par_sigma_max);

                message(1, "increasing par_sigma_range to", new_par_sigma);
                message(1, "regenerating parcov");
                pest_scenario.get_pestpp_options_ptr()->set_par_sigma_range(new_par_sigma);
                parcov.try_from(pest_scenario, file_manager);
                cout << parcov << endl;
            }
            BASE_SCALE_FACTOR = BASE_SCALE_FACTOR * SF_DEC_FAC;
            message(0, "new base scale factor", BASE_SCALE_FACTOR);
        }

        return true;
		
	}
	else
	{
		message(0, "not accepting upgrade #sad");
        best_phis.push_back(last_best);
        best_violations.push_back(last_viol);
        if (infeas_vec[idx] > filter.get_viol_tol())
        {
            n_consec_infeas++;
        }
        if (use_ensemble_grad) {
            double new_par_sigma = pest_scenario.get_pestpp_options().get_par_sigma_range();
            new_par_sigma = new_par_sigma * PAR_SIGMA_DEC_FAC;
            new_par_sigma = max(new_par_sigma, par_sigma_min);
            message(1, "decreasing par_sigma_range to", new_par_sigma);
            message(1, "regenerating parcov");
            parcov.try_from(pest_scenario, file_manager);
            pest_scenario.get_pestpp_options_ptr()->set_par_sigma_range(new_par_sigma);
            cout << parcov << endl;
        }
        BASE_SCALE_FACTOR = BASE_SCALE_FACTOR * SF_INC_FAC;
        message(0, "new base scale factor", BASE_SCALE_FACTOR);
		return false;
	}
}

void SeqQuadProgram::report_and_save_ensemble()
{
	if (use_ensemble_grad)
		report_and_save_ensemble(dv, oe);

}
void SeqQuadProgram::report_and_save_ensemble(ParameterEnsemble& _dv, ObservationEnsemble& _oe)
{
	ofstream& frec = file_manager.rec_ofstream();
	frec << endl << "  ---  SeqQuadProgram iteration " << iter << " report  ---  " << endl;
	frec << "   number of active realizations:  " << _dv.shape().first << endl;
	frec << "   number of model runs:           " << run_mgr_ptr->get_total_runs() << endl;

	cout << endl << "  ---  SeqQuadProgram iteration " << iter << " report  ---  " << endl;
	cout << "   number of active realizations:   " << _dv.shape().first << endl;
	cout << "   number of model runs:            " << run_mgr_ptr->get_total_runs() << endl;
	save(_dv, _oe);
}

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
//		//throw runtime_error("unknown 'subset_how'");
//		throw_sqp_error("unknown 'subset_how'");
//	}
//	stringstream ss;
//	for (auto i : subset_idxs)
//		ss << i << ":" << pe_names[i] << ", ";
//	message(1,"subset idx:dv real name: ",ss.str());
//	return;
//	//return subset_idx_map;
//}

ObservationEnsemble SeqQuadProgram::run_candidate_ensemble(ParameterEnsemble& dv_candidates)
{
	run_mgr_ptr->reinitialize();
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
			_oe = ObservationEnsemble(&pest_scenario);

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

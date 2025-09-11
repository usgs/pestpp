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
	cout << endl << parcov << endl;
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

	message(1, "using the following upgrade vector scale (e.g. 'line search') values:", ppo->get_sqp_alpha_mults());
	
	//ofstream &frec = file_manager.rec_ofstream();
	last_best = 1.0E+30;
	last_viol = 0.0;
	
	warn_min_reals = 10;
	error_min_reals = 2;
	
	//vector<double> scale_facs = pest_scenario.get_pestpp_options().get_lambda_scale_vec();
	//message(1, "using scaling factors: ", scale_facs);
	set<string> passed = ppo->get_passed_args();
	if (passed.find("sqp_alpha_mults") == passed.end())
	{
	    if ((use_ensemble_grad) && (SOLVE_EACH_REAL))
        {
	        message(1,"'sqp_alpha_mults' not passed, using ensemble gradient, and solving each real, resetting scale facs");
	        vector<double> new_scale_facs{0.000001,0.00001,0.0005,0.01,.1};
	        message(1,"new sqp_alpha_mults",new_scale_facs);
	        ppo->set_sqp_alpha_mults(new_scale_facs);
        };
	}

	
	message(1, "max run fail: ", ppo->get_max_run_fail());

	//TODO: update sanity checks for SQP context
	//check that if using fd, chance points == single
	use_ensemble_grad = false;
	if (ppo->get_sqp_num_reals() > 0)
	{
		use_ensemble_grad = true;
		sampling_tracking_initialized = false;
	}
	sanity_checks();

	bool echo = false;
	if (verbose_level > 1)
		echo = true;

	initialize_parcov();
	if (use_cmaes)
	{
		cmaes = CovMatAdapES(&pest_scenario, &rand_gen);
		
		cmaes.initialize(dv_names.size(), ppo->get_sqp_num_reals());
		cmaes.set_covariance(parcov.get_matrix());
	}
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

	working_set_tol = pest_scenario.get_pestpp_options().get_sqp_working_set_tol();

	use_localization = true; //TODO: add ++args

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
	filter.set_tol(pest_scenario.get_pestpp_options().get_sqp_filter_tol());
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
		
		ParameterEnsemble _dv(&pest_scenario, &rand_gen);
		ofstream& frec = file_manager.rec_ofstream();

		if (use_cmaes)
		{
			message(1, "generating new dv ensemble at current best phi via CMA-ES");
			Covariance cmaes_cov = cmaes.get_covariance_matrix(dv_names);

			_dv = cmaes.generate_population(current_ctl_dv_values, dv);
			//_dv.draw(pest_scenario.get_pestpp_options().get_sqp_num_reals(), _current_dv_vals, cmaes_cov, performance_log, 0, frec);
		}
		else
		{
			message(1, "generating new dv ensemble at current best phi");
			_dv.draw(pest_scenario.get_pestpp_options().get_sqp_num_reals(), _current_dv_vals, parcov, performance_log, 0, frec);
		}

		cout << "dv: " << endl << _dv.get_eigen() << endl;
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
	message(1, "using stochastic gradient approximation (StoSAG)");
	
	//I think a bad phi option has use in SQP?
	/*double bad_phi = pest_scenario.get_pestpp_options().get_ies_bad_phi();
	if (bad_phi < 1.0e+30)
		message(1, "using bad_phi: ", bad_phi);*/

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
	current_ctl_dv_values.update(dv.get_var_names(), dv.get_real_vector("BASE"));
	current_obs.update(oe.get_var_names(), oe.get_real_vector("BASE"));

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

bool SeqQuadProgram::try_modify_hessian() {
	// Get current Hessian matrix
	Eigen::MatrixXd H = *hessian.e_ptr();

	// Compute eigendecomposition
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(H);
	if (eigensolver.info() != Eigen::Success) {
		return false;
	}

	// Get eigenvalues and eigenvectors
	Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
	Eigen::MatrixXd eigenvectors = eigensolver.eigenvectors();

	// Check if modification is needed
	double min_eig = eigenvalues.minCoeff();
	const double min_allowed_eig = 1e-3;  // Minimum allowed eigenvalue

	if (min_eig >= min_allowed_eig) {
		return true; // No modification needed
	}

	//Algorithm 3.3, p. 51 Nocedal and Wright
	double tau = 2 * abs(min_eig) + min_allowed_eig;
	Eigen::MatrixXd modified_H = H + tau * Eigen::MatrixXd::Identity(H.rows(), H.cols()); 

	// Update the Hessian
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
	if (std::abs(denominator) < sr1_threshold * y_minus_Hs.norm() * s_k.norm())
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

bool SeqQuadProgram::hessian_update_bfgs(Eigen::VectorXd s_k, Eigen::VectorXd y_k, Covariance old_hessian)
{
	stringstream ss;
	//ofstream& frec = file_manager.rec_ofstream();
	ss.str("");

	ss << "starting BFGS hessian update for iteration " << iter << endl;
	message(1, ss.str());
	// quasi-Newton Hessian updating via BFGS
	// TODO: check if there are conditions to satisfy (some of which are user-specified)

	const double eps = 1e-10;
	const double max_condition = 1e8;
	const double damping_factor = 0.2;

	if (s_k.norm() < eps || y_k.norm() < eps)
	{
		ss << "skipping BFGS update - step or gradient difference too small" << endl;
		message(1, "skipping BFGS update - step or gradient difference too small");
		performance_log->log_event(ss.str());
		return false;
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
	}

	// BFGS Update formula with scaling
	Eigen::MatrixXd H = *old_hessian.e_ptr();
	Eigen::MatrixXd H_new = H;

	// Initial scaling factor (Nocedal & Wright scaling)
	if (iter == 1)
	{
		double scale = s_dot_y / (y_k.squaredNorm());
		H_new *= scale;
		ss << "applying initial scaling factor: " << scale << endl;
	}

	// First term: H_k*s_k*s_k^T*H_k from Eq. 6.19, Nocedal and Wright, p. 140
	Eigen::VectorXd Hs = H_new * s_k;
	H_new -= (Hs * s_k.transpose() * H_new) / (s_k.dot(Hs));

	// Second term: y_k*y_k^T/(y_k^T*s_k) from Eq. 6.19, Nocedal and Wright, p. 140
	H_new += (y_k * y_k.transpose()) / s_dot_y;

	// Check condition number per Nocedal and Wright p. 117. This is required to maintain stability
	// Condition number is the ration between the max eigenvalue and min eigenvalue
	// Too high condition number means that the matrix is close to being singular
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(H_new);
	double cond = eigensolver.eigenvalues().maxCoeff() /
		eigensolver.eigenvalues().minCoeff();

	if (cond > max_condition)
	{
		ss << "warning: condition number too large: " << cond << endl;
		message(1, "warning: condition number too large: ", cond);
		// Apply regularization
		double min_eig = eigensolver.eigenvalues().minCoeff();
		if (min_eig < eps)
		{
			double reg = eps - min_eig;
			H_new += reg * Eigen::MatrixXd::Identity(H_new.rows(), H_new.cols());
			ss << "applying regularization: " << reg << endl;
			message(2, "applying regularization: ", reg);
		}
	}

	hessian = Covariance(dv_names, H_new.sparseView());
	ss << "BFGS Hessian update complete" << endl;
	message(2, "BFGS Hessian update complete");
	performance_log->log_event(ss.str());

	return true;
}

bool SeqQuadProgram::update_hessian()
{
	stringstream ss;
	ss.str("");

	if (!pest_scenario.get_pestpp_options().get_sqp_update_hessian())
	{
		message(2, "hessian_update is false...");
		return false;
	}

	if (n_consec_failures >= max_consec_failures)
	{
		message(2, "filter rejected too many times, resetting hessian to identity matrix");
		return false;
	}
	
	Covariance old_hessian = hessian;
	
	Eigen::VectorXd prev_grad = grad_vector_map[iter-1].get_data_eigen_vec(dv_names);
	Eigen::VectorXd curr_grad = grad_vector_map[iter].get_data_eigen_vec(dv_names);
	
	//compute gradient difference (y_k) and step (s_k): Eq. 18.13 Nocedal and Wright, p. 536; Eq. 15.25 Andrei, p. 529 
	Eigen::VectorXd y_k = curr_grad - prev_grad; 
	Eigen::VectorXd s_k = current_ctl_dv_values.get_data_eigen_vec(dv_names) - prev_ctl_dv_values.get_data_eigen_vec(dv_names);

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

		if (!curr_cnames.empty()) {
			Eigen::MatrixXd current_jco = current_constraint_mat.e_ptr()->toDense();
			for (int i = 0; i < curr_cnames.size(); i++) {
				int row = constraint_to_row[curr_cnames[i]];
				curr_full_jco.row(row) = current_jco.row(i);

				if (constraint_sense[curr_cnames[i]] == "less_than") {
					curr_full_jco.row(row) *= -1;
				}
			}

			//Eq. 18.21 from Nocedal and Wright, p. 539, approx new lambda w/o computing new Hessian
			Eigen::BDCSVD<Eigen::MatrixXd> svd_AAT(curr_full_jco * curr_full_jco.transpose(), Eigen::ComputeThinU | Eigen::ComputeThinV);
			curr_full_lambda = svd_AAT.solve(curr_full_jco * curr_grad);
		}

		if (!prev_cnames.empty()) {
			Eigen::MatrixXd previous_jco = prev_constraint_mat.e_ptr()->toDense();
			for (int i = 0; i < prev_cnames.size(); i++) {
				int row = constraint_to_row[prev_cnames[i]];
				prev_full_jco.row(row) = previous_jco.row(i);

				if (i < lambda.size()) {
					prev_full_lambda(row) = lambda(i);
				}

				if (constraint_sense[prev_cnames[i]] == "less_than") {
					prev_full_jco.row(row) *= -1;
				}
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

	
	//return hessian_update_sr1(s_k, y_k, old_hessian);
	return hessian_update_bfgs(s_k, y_k, old_hessian);

	cout << "new hessian: " << endl << hessian << endl;
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

		diagonal_scaling(i) = std::max(0.1, std::min(diagonal_scaling(i), 1000.0));
	}

	message(2, "Updated adaptive Hessian scaling");
}

bool SeqQuadProgram::isfullrank(const Eigen::MatrixXd& mat)
{
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::VectorXd singularValues = svd.singularValues();
	double tol = 1e-5;
	int rank = (singularValues.array() > tol).count();
	if (rank == min(mat.rows(), mat.cols()))
		return true;
	else
		return false;
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
 		
		double ws = working_set_tol;

		if (use_ensemble_grad)
			accept = solve_new_ensemble();
		else
			accept = solve_new();
		//accept = solve_new();
		if (accept)
			working_set_tol = max(0.05, ws * 0.5);


        if (use_ensemble_grad)
            report_and_save_ensemble();
        else
            save_current_dv_obs();
        
		if (use_cmaes)
		{
			message(1, "updating CMA-ES with approximate gradient");
			
			cmaes.update(prev_ctl_dv_values, current_ctl_dv_values);
		}
		make_gradient_runs(current_ctl_dv_values, current_obs);

		grad_vector_map[iter] = calc_gradient_vector(current_ctl_dv_values);
		current_grad_vector = grad_vector_map[iter];

        constraints.sqp_report(iter, current_ctl_dv_values, current_obs, true);

        if (use_ensemble_grad)
        {
            ss.str("");
            ss << file_manager.get_base_filename() << "." << iter << ".pcs.csv";
            pcs.summarize(dv, ss.str());
        }

		if (should_terminate())
			break;

		if (n_consec_infeas > MAX_CONSEC_INFEAS_IES || !accept)
        {
		    ss.str("");
		    ss << "number of consecutive infeasible iterations > " << MAX_CONSEC_INFEAS << ", switching to IES to seek feasibility";
		    message(0,ss.str());
		    seek_feasible();
		    n_consec_infeas = 0;
        }

		if (pest_scenario.get_pestpp_options().get_sqp_update_hessian())
			update_hessian();
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
        if (best_phis[i]<best_phi_yet)
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

    ss.str("");
    ss << "best phi,infeas sequence:" << endl;
    int ii = 0;
    for (int i=0;i<best_phis.size();i++)
    {
        phi = best_phis[i];
        infeas = best_violations[i];
        ss << "    " << setw(5) << setprecision(4) << right << phi << "," << setw(5) << setprecision(4) << left << infeas << endl;
        ii++;
    }
    ss << endl;
    message(0, ss.str());

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
	if (converged)
	{
		message(1, "optimal solution detected at solve EQP step (lagrangian multiplier for all ineq constraints in working set is non-neg)");
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
	ofstream& frec = file_manager.rec_ofstream();

	string center_on = pest_scenario.get_pestpp_options().get_ies_center_on();
	if (!_center_on.empty())
	    center_on = _center_on;
	
	if (use_ensemble_grad)
	{
		//ensemble stuff here
		//if (use_obj_obs)
		{
			// compute sample dec var cov matrix and its pseudo inverse
			// see eq (8) of Dehdari and Oliver 2012 SPE and Fonseca et al 2015 SPE
			// TODO: so can pseudo inverse: Covariance dv_cov_matrix; 
			//Eigen::MatrixXd parcov_inv;
			// start by computing mean-shifted dec var ensemble
			Eigen::MatrixXd dv_anoms = dv.get_eigen_anomalies(vector<string>(), dv_names, "BASE"); 
			//dv_anoms.conservativeResize(dv_anoms.rows() - 1, dv_anoms.cols());
			/*if (dv.shape().first > 1000)  // until we encounter
			{
				// lower rank - diag elements only
				throw_sqp_error("TODO: use dv.get_diagonal_cov matrix()? need to check for consistency if so"); 
				//parcov = dv.get_diagonal_cov_matrix();  // check ok to instantiate empricially here  // pass helpful center_on here too
				//Eigen::VectorXd parcov_inv;
				//Covariance parcov_diag;
				//parcov_diag.from_diagonal(parcov);
				//parcov_inv = parcov_diag.get_matrix().diagonal();
				//parcov_inv = parcov_inv.cwiseInverse();  // equivalent to pseudo inv?
			}*/
			Eigen::MatrixXd dv_cov_matrix = 1.0 / (dv.shape().first - 1.0) * (dv_anoms.transpose() * dv_anoms);
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
            Eigen::MatrixXd s, V, U, st;
            SVD_REDSVD rsvd;
            //SVD_EIGEN rsvd;
            rsvd.set_performance_log(performance_log);
            //dv_anoms.transposeInPlace();
            rsvd.solve_ip(dv_cov_matrix, s, U, V, pest_scenario.get_svd_info().eigthresh, pest_scenario.get_svd_info().maxsing);
			Eigen::MatrixXd dv_cov_pseudoinv = V * s.asDiagonal().inverse() * U.transpose();
			Eigen::MatrixXd obj_anoms(dv.shape().first,1);
            if (use_obj_obs) {
                obj_anoms = oe.get_eigen_anomalies(vector<string>(), vector<string>{obj_func_str},"BASE");
				//obj_anoms.conservativeResize(obj_anoms.rows() - 1, obj_anoms.cols());
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
			Eigen::VectorXd cross_cov_vector = 1.0 / (dv.shape().first - 1.0) * (dv_anoms.transpose() * obj_anoms);

			// now compute grad vector
			// this is a matrix-vector product; the matrix being the pseudo inv of diag empirical dec var cov matrix and the vector being the dec var-phi cross-cov vector\
			// see, e.g., Chen et al. (2009) and Fonseca et al. (2015) 
			grad = dv_cov_pseudoinv * cross_cov_vector;

			//if (use_localization)
			//{
			//	Eigen::VectorXd loc_weights;

			//	loc_weights = compute_adaptive_localization_weights(dv_names, dv_anoms);
			//	{
			//		stringstream ss;
			//		ss << "Applied localization. Weights: " << loc_weights.transpose();
			//		frec << ss.str() << endl;
			//	}
			//	grad = grad.cwiseProduct(loc_weights);

			//}

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
	ss << "approximate gradient" << endl << grad << endl;
	frec << ss.str();

	return pgrad;
}

pair<Eigen::VectorXd, Eigen::VectorXd> SeqQuadProgram::_kkt_null_space(Eigen::MatrixXd& G, Eigen::MatrixXd& _constraint_jco, Eigen::VectorXd& constraint_diff, Eigen::VectorXd& curved_grad)
{
	stringstream ss;
	ss.str("");
	ss << "starting KKT null space solve";
	performance_log->log_event(ss.str());

	Eigen::VectorXd search_d, lm;
	
	Eigen::VectorXd x;
	Eigen::MatrixXd V, U, S_, s;
	SVD_REDSVD rsvd;
	ss.str("");
	ss << "using randomized SVD to compute basis matrices of constraint JCO for null space KKT solve " << _constraint_jco;
	performance_log->log_event(ss.str());

	Eigen::BDCSVD<Eigen::MatrixXd> svd_A(_constraint_jco, Eigen::DecompositionOptions::ComputeFullU | 
														Eigen::DecompositionOptions::ComputeFullV);
	s = svd_A.singularValues();
	U = svd_A.matrixU(); 
	V = svd_A.matrixV(); 

	double rank_tol = pest_scenario.get_svd_info().eigthresh;
	int rank = (s.array() > rank_tol).count();

	if (rank < _constraint_jco.rows())
	{
		stringstream ss;
		ss << "Constraint matrix does not have full row rank. Rank = " << rank << ", Rows = " << _constraint_jco.rows();
		throw_sqp_error(ss.str());
	}

	Eigen::MatrixXd Y = V.leftCols(rank);
	Eigen::MatrixXd Z;
	if (V.cols() > rank)
		Z = V.rightCols(V.cols() - rank);
	else
		Z = Eigen::MatrixXd::Zero(V.rows(), 0);
	

	// solve p_range_space
	Eigen::VectorXd p_y, rhs;
	Eigen::MatrixXd coeff = _constraint_jco * Y; //AY Eq. 16.18 Nocedal and Wright, pp. 457
	Eigen::BDCSVD<Eigen::MatrixXd> AY(coeff, Eigen::ComputeThinU | Eigen::ComputeThinV);
	p_y = AY.solve(-constraint_diff); //from Eq 18.19a pp. 538 Nocedal and Wright

	// todo assert here for p_y != 0 if sum_viol == 0 (this is the component that rectifies constraint viol)

	// solve p_null_space
	bool simplified_null_space_approach = false; //TODO: add as ++arg; too much to be an arg?
	Eigen::VectorXd p_z;

	if (Z.size() > 0)
	{
		Eigen::MatrixXd red_hess = Z.transpose() * G * Z; 

		//check if positive definite
		Eigen::LDLT<Eigen::MatrixXd> ldlt(red_hess);
		if (ldlt.info() != Eigen::Success || !ldlt.isPositive())
			throw_sqp_error("Reduced Hessian Z^T G Z is not positive definite");

		bool cholesky = false; //TODO: add as ++arg; too much to be an arg?
		if (simplified_null_space_approach)
		{
			message(1, "using simplified approach in KKT null space solve...");
			rhs = - Z.transpose() * curved_grad;
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
		Eigen::BDCSVD<Eigen::MatrixXd> svd_AAT(_constraint_jco* _constraint_jco.transpose(),Eigen::ComputeThinU | Eigen::ComputeThinV);
		lm = svd_AAT.solve(_constraint_jco * curved_grad);
	}
	else
	{
		// Nocedal and Wright pg. 457 and 538
		rhs = Y.transpose() * (curved_grad + G * search_d); 
		coeff = (_constraint_jco * Y).transpose();
		Eigen::BDCSVD<Eigen::MatrixXd> AY(coeff,Eigen::ComputeThinU | Eigen::ComputeThinV);
		lm = AY.solve(rhs);
	}
	return pair<Eigen::VectorXd, Eigen::VectorXd>(search_d, lm);
}

pair<Eigen::VectorXd, Eigen::VectorXd> SeqQuadProgram::_kkt_direct(Eigen::MatrixXd& G, Eigen::MatrixXd& _constraint_jco, Eigen::VectorXd& constraint_diff, Eigen::VectorXd& curved_grad, vector<string>& cnames)
{
	// 1. Check if we can use a more efficient approach for special cases
	if (cnames.empty()) {
		// No constraints - just solve the unconstrained problem
		Eigen::VectorXd search_d = G.ldlt().solve(-curved_grad);
		return pair<Eigen::VectorXd, Eigen::VectorXd>(search_d, Eigen::VectorXd());
	}

	// 2. Form the KKT matrix with proper regularization
	int n = dv_names.size();
	int m = cnames.size();

	// Apply regularization to G if needed to ensure positive definiteness
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(G);
	double min_eig = eigensolver.eigenvalues().minCoeff();
	const double min_allowed_eig = 1e-6;

	Eigen::MatrixXd G_reg = G;
	if (min_eig < min_allowed_eig) {
		double delta = min_allowed_eig - min_eig + 1e-6;
		G_reg += delta * Eigen::MatrixXd::Identity(n, n);
		message(1, "Applied regularization to Hessian, delta = ", delta);
	}

	// Form the KKT matrix with proper scaling
	// Following Nocedal & Wright Eq. 16.62
	double constraint_scaling = 1.0;
	if (G_reg.norm() > 1e-8) {
		constraint_scaling = std::sqrt(G_reg.norm());
	}

	Eigen::MatrixXd scaled_constraint_jco = constraint_scaling * _constraint_jco;

	Eigen::MatrixXd kkt_matrix(n + m, n + m);
	kkt_matrix.topLeftCorner(n, n) = G_reg;
	kkt_matrix.topRightCorner(n, m) = scaled_constraint_jco.transpose();
	kkt_matrix.bottomLeftCorner(m, n) = scaled_constraint_jco;
	kkt_matrix.bottomRightCorner(m, m) = Eigen::MatrixXd::Zero(m, m);

	// 3. Form the right-hand side
	// Following Nocedal & Wright Eq. 16.4
	Eigen::VectorXd rhs(n + m);
	rhs.head(n) = -curved_grad;
	rhs.tail(m) = -constraint_diff * constraint_scaling;

	// 4. Solve the KKT system using an appropriate method
	// For stability, use LDLT factorization with pivoting
	Eigen::LDLT<Eigen::MatrixXd> ldlt(kkt_matrix);

	// Check if factorization succeeded
	if (ldlt.info() != Eigen::Success) {
		message(1, "LDLT factorization failed, falling back to SVD");

		// Use Eigen's built-in SVD instead of custom SVD_REDSVD
		Eigen::BDCSVD<Eigen::MatrixXd> svd(kkt_matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);

		// Solve the system using SVD
		Eigen::VectorXd x = svd.solve(rhs);

		Eigen::VectorXd search_d = x.head(n);
		Eigen::VectorXd lm = x.tail(m) / constraint_scaling; // Rescale back

		return pair<Eigen::VectorXd, Eigen::VectorXd>(search_d, lm);
	}

	// Solve using LDLT
	Eigen::VectorXd x = ldlt.solve(rhs);

	// 5. Extract the solution components
	Eigen::VectorXd search_d = x.head(n);
	Eigen::VectorXd lm = x.tail(m) / constraint_scaling; // Rescale back

	// 6. Verify the solution
	double kkt_error = (kkt_matrix * x - rhs).norm() / (1.0 + rhs.norm());
	if (kkt_error > 1e-6) {
		message(1, "Warning: KKT system solution has high residual: ", kkt_error);
	}

	message(1, "KKT system solved with residual: ", kkt_error);
	message(2, "search_d: ", search_d.transpose());
	message(2, "lagrange multipliers: ", lm.transpose());

	return pair<Eigen::VectorXd, Eigen::VectorXd>(search_d, lm);
}
//
//pair<Eigen::VectorXd, Eigen::VectorXd> SeqQuadProgram::_kkt_direct(Eigen::MatrixXd& G, Eigen::MatrixXd& constraint_jco, Eigen::VectorXd& constraint_diff, Eigen::VectorXd& curved_grad, vector<string>& cnames)
//{
//	
//	//check A full rank
//
//	// forming system to be solved - this is filth but it works..
//	Eigen::MatrixXd coeff_u(dv_names.size(), dv_names.size() + cnames.size());  // todo only in WS
//	coeff_u << G, constraint_jco.transpose();
//	Eigen::MatrixXd coeff_l(cnames.size(), dv_names.size() + cnames.size());  // todo only in WS
//	coeff_l << constraint_jco, Eigen::MatrixXd::Zero(cnames.size(), cnames.size());
//	Eigen::MatrixXd coeff(dv_names.size() + cnames.size(), dv_names.size() + cnames.size());  // todo only in WS
//	coeff << coeff_u, coeff_l;
//	message(1, "coeff", coeff);  // tmp
//
//	Eigen::VectorXd rhs(curved_grad.size() + constraint_diff.size());
//	rhs << curved_grad, constraint_diff;  // << vec1, vec2;
//	message(1, "rhs", rhs);  // tmp
//
//	Eigen::VectorXd x;
//	Eigen::MatrixXd V, U, S_, s;
//	SVD_REDSVD rsvd;
//	//SVD_EIGEN rsvd;
//	rsvd.set_performance_log(performance_log);
//
//	rsvd.solve_ip(coeff, s, U, V, pest_scenario.get_svd_info().eigthresh, pest_scenario.get_svd_info().maxsing);
//	S_ = s.asDiagonal();
//	message(1, "singular values of KKT matrix", S_);  // tmp
//
//	// an old friend!
//	x = V * S_.inverse() * U.transpose() * rhs;
//	message(1, "solution vector of steps and lagrange mults", x);  // tmp
//
//	Eigen::VectorXd search_d, lm;
//	search_d = x.head(dv_names.size());  // add rigor here or at least asserts to ensure operating on correct elements
//	lm = x.tail(x.size() - dv_names.size());  // add rigor here or at least asserts to ensure operating on correct elements
//
//	
//	return pair<Eigen::VectorXd, Eigen::VectorXd> (search_d, lm);
//}

pair<Mat, bool> SeqQuadProgram::get_constraint_mat(Parameters& _dv_vals, Observations& _obs_vals, double working_set_tol, const Eigen::VectorXd* lm)
{
	if (use_ensemble_grad) {
		return constraints.get_working_set_constraint_matrix(_dv_vals, _obs_vals, dv, oe, true, lm, (working_set_tol));
	}
	else
	{
		message(2, "getting working set constraint matrix");
		return constraints.get_working_set_constraint_matrix(_dv_vals, _obs_vals, jco, true, lm, (working_set_tol));
	}
}

Eigen::VectorXd SeqQuadProgram::solve_trust_region_subproblem(const Eigen::MatrixXd& B,	const Eigen::VectorXd& g, double radius)
{
	//Steihaug-Point CG method algorithm 7.2, p. 171 Nocedal and Wright
	int n = g.size();
	Eigen::VectorXd p = Eigen::VectorXd::Zero(n);
	Eigen::VectorXd r = g;  // r = g + Bp = g + 0 initially
	Eigen::VectorXd d = -r;

	cout << endl << B << endl;
	double g_norm = g.norm();
	double tol = min(0.1, sqrt(g_norm)) * g_norm;

	for (int i = 0; i < n; i++)
	{
		Eigen::VectorXd Bd = B * d;
		double dBd = d.dot(Bd);

		if (dBd <= 0)
		{
			return compute_boundary_solution(p, d, radius);
		}

		double alpha = r.dot(r) / dBd;
		Eigen::VectorXd p_new = p + alpha * d;

		if (p_new.norm() >= radius)
		{
			return compute_boundary_solution(p, d, radius);
		}

		Eigen::VectorXd r_new = r + alpha * Bd;

		if (r_new.norm() < tol)
		{
			return p_new;
		}

		double beta = r_new.dot(r_new) / r.dot(r);
		d = -r_new + beta * d;
		p = p_new;
		r = r_new;
	}

	return p;
}

Eigen::VectorXd SeqQuadProgram::compute_boundary_solution(const Eigen::VectorXd& p,
	const Eigen::VectorXd& d, double radius)
{
	//compute tau where ||p + tau*d|| = radius
	double a = d.squaredNorm();
	double b = 2 * p.dot(d);
	double c = p.squaredNorm() - radius * radius;

	//solve: a*tau^2 + b*tau + c = 0
	double discriminant = b * b - 4 * a * c;
	double tau = (-b + sqrt(discriminant)) / (2 * a);

	return p + tau * d;
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

	// Use LDLT decomposition for stability
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
	// Compute predicted reduction using quadratic model Eq 4.2, p. 68 Nocedal and Wright
	double linear_term = -grad.dot(step);
	double quadratic_term = -0.5 * step.dot(hessian.get_matrix() * step);
	return linear_term + quadratic_term;
}

Eigen::VectorXd SeqQuadProgram::compute_adaptive_localization_weights(const vector<string>& dv_names, const Eigen::MatrixXd& dv_anoms)
{
	Eigen::VectorXd weights(dv_names.size());

	Eigen::VectorXd par_std = dv_anoms.colwise().norm() / sqrt(dv_anoms.rows() - 1);

	for (int i = 0; i < dv_names.size(); i++) {
		double uncertainty = par_std[i];
		weights[i] = 1.0 / (1.0 + uncertainty);
	}

	return weights;
}

bool SeqQuadProgram::iterative_partial_step(const string& _blocking_constraint)
{
	Eigen::VectorXd curr_dv_vec = current_ctl_dv_values.get_data_eigen_vec(dv_names);
	Eigen::VectorXd infeas_dv_vec = infeas_cand_dv_values.get_data_eigen_vec(dv_names);
	Eigen::VectorXd full_step =  infeas_dv_vec - curr_dv_vec;

	double alpha_low = 0.0, alpha_high = 1.0, alpha = 0.5;
	const double tol = 1e-4;
	const int max_iter = 8; 
	int iter_count = 0;

	Parameters test_dv_values;
	Observations test_obs;
	bool found_boundary = false;

	while ((alpha_high - alpha_low > tol) && (iter_count <= max_iter)) {
		alpha = (alpha_low + alpha_high) / 2.0;
		test_dv_values = current_ctl_dv_values;
		Eigen::VectorXd test_vec = curr_dv_vec + alpha * full_step;
		test_dv_values.update_without_clear(dv_names, test_vec);

		ParameterEnsemble test_pe(&pest_scenario, &rand_gen);
		test_pe.set_trans_status(ParameterEnsemble::transStatus::NUM);
		test_pe.reserve({ "test_point" }, dv_names);
		test_pe.update_real_ip("test_point", test_vec);
		test_pe.enforce_bounds(performance_log, false);

		test_vec = test_pe.get_real_vector("test_point");
		test_dv_values.update_without_clear(dv_names, test_vec);

		ObservationEnsemble test_oe = run_candidate_ensemble(test_pe);

		if (test_oe.shape().first == 0) {
			message(1, "model run failed at alpha =", alpha);
			alpha_high = alpha;
			continue;
		}
		test_obs = current_obs;
		test_obs.update_without_clear(test_oe.get_var_names(), test_oe.get_real_vector("test_point"));

		map<string, double> test_violations = constraints.get_unsatified_obs_constraints(test_obs, 0.0);
		bool is_violated = (test_violations.find(_blocking_constraint) != test_violations.end());

		stringstream ss;
		ss << "iteration " << iter_count << ", alpha = " << alpha
			<< ", constraint " << (is_violated ? "violated" : "satisfied");
		message(1, ss.str());

		//update search interval
		if (is_violated) {
			alpha_high = alpha;
		}
		else {
			found_boundary = true;
			alpha_low = alpha;

			current_ctl_dv_values = test_dv_values;
			current_obs = test_obs;
		}

		iter_count++;
	}

	if (found_boundary) {
		message(1, "found approximate constraint boundary at alpha =", alpha_low);

		if (find(cnames.begin(), cnames.end(), _blocking_constraint) == cnames.end()) {
			message(1, "adding blocking constraint to working set:", _blocking_constraint);
			cnames.push_back(_blocking_constraint);
			pair<Mat, bool> new_constraint_mat = get_constraint_mat(current_ctl_dv_values, current_obs);
			current_constraint_mat = new_constraint_mat.first;
		}

		constraints.sqp_report(iter, current_ctl_dv_values, current_obs, true,
			"iterative partial step result");

		double obj_val = get_obj_value(current_ctl_dv_values, current_obs);
		double viol_val = constraints.get_sum_of_violations(current_ctl_dv_values, current_obs);
		last_best = obj_val;
		last_viol = viol_val;
		best_phis.push_back(obj_val);
		best_violations.push_back(viol_val);

		message(0, "new best phi and infeas:", vector<double>{last_best, last_viol});
		return true;
	}
	else {
		message(1, "could not find feasible point near constraint boundary");
		return false;
	}
}

bool SeqQuadProgram::line_search(Eigen::VectorXd& search_d, const Parameters& _current_dv_values, Eigen::VectorXd& grad)
{
	double initial_obj = get_obj_value(current_ctl_dv_values, current_obs);
	double initial_slope = grad.dot(search_d);

	//check if we're at a stationary point
	const double grad_tol = 1e-6;  //can be made a class member
	if (grad.norm() < grad_tol) {
		message(1, "Possible stationary point detected - gradient norm below tolerance");
		return false;
	}

	//handle non-descent direction
	if (initial_slope >= 0.1)
	{
		message(1, "Warning: search direction is not a descent direction");
		Covariance original_hessian = hessian;

		// First try: Modify Hessian to make it positive definite
		bool modified_success = try_modify_hessian();

		if (modified_success) {
			// Recompute search direction with modified Hessian
			pair<Eigen::VectorXd, Eigen::VectorXd> x = calc_search_direction_vector(current_ctl_dv_values, current_obs, grad, &base_constraint_jco);
			search_d = x.first;
			initial_slope = grad.dot(search_d);
		}

		// If modification fails or still gives non-descent direction
		if (!modified_success || initial_slope >= 0) {
			message(1, "Resetting Hessian to scaled identity matrix");
			Eigen::SparseMatrix<double> h(dv_names.size(), dv_names.size());
			h.setIdentity();
			update_scaling(search_d, grad);
			for (int i = 0; i < dv_names.size(); i++) {
				h.coeffRef(i, i) *= diagonal_scaling(i);
			}

			Covariance identity_hessian(dv_names, h);
			hessian = identity_hessian;
			pair<Eigen::VectorXd, Eigen::VectorXd> x = calc_search_direction_vector(current_ctl_dv_values, current_obs, grad, &base_constraint_jco);
			search_d = x.first;
			initial_slope = grad.dot(search_d);

			// If still not descent, use pure steepest descent and restore original Hessian
			if (initial_slope >= 0) {
				message(1, "Using pure steepest descent and restoring original Hessian");
				search_d = -grad;
				initial_slope = grad.dot(search_d);

				hessian = original_hessian;
				cout << endl << "hessian" << endl << hessian << endl;
			}
		}
	}


	stringstream ss;
	ParameterEnsemble dv_candidates(&pest_scenario, &rand_gen);
	dv_candidates.set_trans_status(ParameterEnsemble::transStatus::NUM);

	vector<string> real_names;
	vector<double> scale_vals;
	scale_vals.clear();
	for (auto& sf : pest_scenario.get_pestpp_options().get_sqp_alpha_mults())
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
	map<string, double> real_sf_map;
	for (int i = 0; i < scale_vals.size(); i++)
	{
		double scale_val = scale_vals[i];
		ss.str("");
		ss << "starting calcs for scaling factor" << scale_val;
		message(1, "starting lambda calcs for scaling factor", scale_val);
		message(2, "see .log file for more details");


		Parameters num_candidate = _current_dv_values;

		Eigen::VectorXd scale_search_d = search_d * scale_val;
		if (scale_search_d.squaredNorm() < 1.0 - 10)
			message(1, "very short upgrade for scale value", scale_val);

		Eigen::VectorXd cvals = num_candidate.get_data_eigen_vec(dv_names);
		cvals.array() += (scale_search_d / search_d.norm()).array();
		num_candidate.update_without_clear(dv_names, cvals);

		Eigen::VectorXd vec = num_candidate.get_data_eigen_vec(dv_names);
		dv_candidates.update_real_ip(real_names[i], vec);
		used_scale_vals.push_back(scale_val);
		real_sf_map[real_names[i]] = scale_val;
		

		ss.str("");
		message(1, "finished calcs for scaling factor:", scale_val);

	}

	if (pest_scenario.get_pestpp_options().get_ies_debug_upgrade_only())
	{
		message(0, "ies_debug_upgrade_only is true, exiting");
		throw_sqp_error("ies_debug_upgrade_only is true, exiting");
	}

	//enforce bounds on candidates - TODO: report the shrinkage summary that enforce_bounds returns
	dv_candidates.enforce_bounds(performance_log, false);
	ss.str("");
	ss << file_manager.get_base_filename() << "." << iter << ".dv_candidates.csv";
	dv_candidates.to_csv(ss.str());

	Eigen::VectorXd v1, v2;
	double d;
	vector<string> drop;
	set<int> jvals;
	for (int i = 0; i < dv_candidates.shape().first; i++)
	{
		v1 = dv_candidates.get_real_vector(i);
		for (int j = i + 1; j < dv_candidates.shape().first; j++) {
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
		message(1, "dropping the following duplicate candidates: ", drop);
		dv_candidates.drop_rows(drop, true);
		used_scale_vals.clear();
		for (auto& real_name : dv_candidates.get_real_names())
		{
			used_scale_vals.push_back(real_sf_map.at(real_name));
		}

	}

	message(0, "running candidate decision variable batch");
	vector<double> passed_scale_vals = scale_vals;

	ObservationEnsemble oe_candidates = run_candidate_ensemble(dv_candidates);
	ss.str("");
	ss << file_manager.get_base_filename() << "." << iter << ".oe_candidates.csv";
	oe_candidates.to_csv(ss.str());

	map<string, double> sf_map;
	for (int i = 0; i < dv_candidates.get_real_names().size(); i++)
	{
		sf_map[dv_candidates.get_real_names()[i]] = used_scale_vals[i];
	}
	return pick_candidate_and_update_current(dv_candidates, oe_candidates, sf_map);
}

bool SeqQuadProgram::line_search(map<string, Eigen::VectorXd>& search_d, Eigen::VectorXd& grad, ParameterEnsemble* dvs_subset, bool first_pick)
{
	stringstream ss;
	ss.str("");
	ss << "performing line search on ensemble subset:";
	for (auto& n : dvs_subset->get_real_names())
		ss << " " << n << ",";
	ss << " BASE" << endl;
	message(1, ss.str());

	const double grad_tol = 1e-6;
	if (grad.norm() < grad_tol) {
		message(1, "Possible stationary point detected - gradient norm below tolerance");
		return false;
	}

	/*double initial_obj = get_obj_value(current_ctl_dv_values, current_obs);
	double initial_slope = grad.dot(search_d);*/

	//TODO: REVISIT FOR ENSEMBLE handle non-descent direction
	//if (initial_slope >= 0.1) 
	//{
	//	message(1, "Warning: search direction is not a descent direction");
	//	Covariance original_hessian = hessian;
	//	
	//	// First try: Modify Hessian to make it positive definite
	//	bool modified_success = try_modify_hessian();

	//	if (modified_success) {
	//		// Recompute search direction with modified Hessian
	//		pair<Eigen::VectorXd, Eigen::VectorXd> x = calc_search_direction_vector(_current_dv_values, grad);
	//		search_d = x.first;
	//		initial_slope = grad.dot(search_d);
	//	}

	//	// If modification fails or still gives non-descent direction
	//	if (!modified_success || initial_slope >= 0) {
	//		message(1, "Resetting Hessian to scaled identity matrix");
	//		Eigen::SparseMatrix<double> h(dv_names.size(), dv_names.size());
	//		h.setIdentity();
	//		update_scaling(search_d, grad);
	//		for (int i = 0; i < dv_names.size(); i++) {
	//			h.coeffRef(i, i) *= diagonal_scaling(i);
	//		}
	//		
	//		Covariance identity_hessian(dv_names, h);
	//		hessian = identity_hessian;
	//		pair<Eigen::VectorXd, Eigen::VectorXd> x = calc_search_direction_vector(_current_dv_values, grad);
	//		search_d = x.first;
	//		initial_slope = grad.dot(search_d);

	//		// If still not descent, use pure steepest descent and restore original Hessian
	//		if (initial_slope >= 0) {
	//			message(1, "Using pure steepest descent and restoring original Hessian");
	//			search_d = -grad;
	//			initial_slope = grad.dot(search_d);

	//			hessian = original_hessian;
	//			cout << endl << "hessian" << endl << hessian << endl;
	//		}
	//	}
	//}

	ParameterEnsemble dv_candidates(&pest_scenario, &rand_gen);
	dv_candidates.set_trans_status(ParameterEnsemble::transStatus::NUM);

	map<double, map<string, string>> sv_real_map;
	map<string, Eigen::VectorXd> search_dirs;
	vector<string> real_names;
	vector<double> scale_vals;
	
	scale_vals.clear();
	for (auto& sf : pest_scenario.get_pestpp_options().get_sqp_alpha_mults())
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
	else if (use_ensemble_grad)
	{
		if (dvs_subset != nullptr) 
		{
			ParameterEnsemble d;
			const auto& names = dvs_subset->get_real_names();
			if (std::find(names.begin(), names.end(), "BASE") == names.end())
			{
				d.reserve(vector<string>{ "BASE" }, pest_scenario.get_ctl_ordered_par_names());
				d.add_2_row_ip("BASE", dv.get_real_vector("BASE"));
				dvs_subset->append_other_rows(d);
			}
			
			for (auto& real_name : dvs_subset->get_real_names()) 
			{
				for (auto sv : scale_vals) 
				{
					ss.str("");
					ss << "dv_cand_" << real_name << "_sv:" << left << setw(8) << setprecision(3) << sv;
					real_names.push_back(ss.str());
					sv_real_map[sv][real_name] = ss.str();
				}
				search_dirs[real_name] = search_d.at(real_name);
			}
		}
		else
			throw_sqp_error("use_ensemble_grad is true but subset dv ensemble is null");
	}
	dv_candidates.reserve(real_names, dv_names);

	vector<double> used_scale_vals;
	map<string, double> real_sf_map;
	for (int i = 0;i < scale_vals.size();i++)
	{
		double scale_val = scale_vals[i];
		ss.str("");
		ss << "starting calcs for scaling factor" << scale_val;
		message(1, "starting lambda calcs for scaling factor", scale_val);
		message(2, "see .log file for more details");

		map<string, Eigen::VectorXd> scale_search_d;
		vector<string> short_upgrades;
		for (auto& sd : search_dirs)
		{
			scale_search_d[sd.first] = sd.second * scale_val / sd.second.norm();
			if (scale_search_d[sd.first].squaredNorm() < 1.0 - 10)
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
			dv_upgrade += scale_search_d[real_name];			
			
			string upgrade_rname = sv_real_map[scale_val][real_name];
			dv_candidates.update_real_ip(upgrade_rname, dv_upgrade);
			real_sf_map[upgrade_rname] = scale_val;
			used_scale_vals.push_back(scale_val);
		}
			
		ss.str("");
		message(1, "finished calcs for scaling factor:", scale_val);

	}

	if (pest_scenario.get_pestpp_options().get_ies_debug_upgrade_only())
	{
		message(0, "ies_debug_upgrade_only is true, exiting");
		throw_sqp_error("ies_debug_upgrade_only is true, exiting");
	}

	//enforce bounds on candidates - TODO: report the shrinkage summary that enforce_bounds returns
	dv_candidates.enforce_bounds(performance_log, false);
	ss.str("");
	ss << file_manager.get_base_filename() << "." << iter << ".dv_candidates.csv";
	dv_candidates.to_csv(ss.str());

	Eigen::VectorXd v1, v2;
	double d;
	vector<string> drop;
	set<int> jvals;
	for (int i = 0;i < dv_candidates.shape().first;i++)
	{
		v1 = dv_candidates.get_real_vector(i);
		for (int j = i + 1;j < dv_candidates.shape().first;j++) {
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
		message(1, "dropping the following duplicate candidates: ", drop);
		dv_candidates.drop_rows(drop, true);
		used_scale_vals.clear();
		for (auto& real_name : dv_candidates.get_real_names())
		{
			used_scale_vals.push_back(real_sf_map.at(real_name));
		}

	}
	message(0, "running candidate decision variable batch");
	ObservationEnsemble oe_candidates = run_candidate_ensemble(dv_candidates);
	ss.str("");
	ss << file_manager.get_base_filename() << "." << iter << ".oe_candidates.csv";
	oe_candidates.to_csv(ss.str());

	return pick_candidate_and_update_current(dv_candidates, oe_candidates, real_sf_map, first_pick);
}

bool SeqQuadProgram::check_wolfe_conditions(Parameters& trial_dv_values, Observations& trial_obs, const Eigen::VectorXd& search_d,
	const Eigen::VectorXd& grad, double scale, double initial_obj, double initial_slope)
{
	//Algorithm 3.5, pp. 60-61 Nocedal and Wright
	double trial_obj = get_obj_value(trial_dv_values, trial_obs);

	//Armijo condition (sufficient decrease)
	if (trial_obj > initial_obj + c1 * scale * initial_slope) //Eq. 3.6a pp. 34, Nocedal and Wright
		return false;

	//get new gradient at trial point
	Parameters trial_grad = calc_gradient_vector(trial_dv_values);
	Eigen::VectorXd trial_grad_vec = trial_grad.get_data_eigen_vec(dv_names);

	//curvature condition
	double trial_slope = trial_grad_vec.dot(search_d);
	if (abs(trial_slope) > c2 * abs(initial_slope)) //Eq. 3.6b pp. 34, Nocedal and Wright
		return false;

	return true;
}

double SeqQuadProgram::get_reference_obj()
{
	if (previous_obj_values.empty()) {
		return get_obj_value(current_ctl_dv_values, current_obs);
	}
	return *max_element(previous_obj_values.begin(), previous_obj_values.end());
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
		//message(0, "current working set:", Cnames);

		Eigen::VectorXd constraint_diff(Cnames.size());
		pair<Eigen::VectorXd, Eigen::VectorXd> p = constraints.get_obs_resid_constraint_vectors(_current_dv_values, _current_obs_values, Cnames);
		constraint_diff = p.second;

		for (int i = 0;i < Cnames.size();i++) {
			if (constraint_sense[Cnames[i]] == "less_than")
			{
				if (use_ensemble_grad)
					constr_jco.row(i) *= -1;
				constraint_diff[i] *= -1;
			}
		}

		if ((constr_jco.rows() > 0) && (!isfullrank(constr_jco)))
        {
			message(0, "WARNING: constraint_jco is not full rank. Using complete orthogonal decomposition.");
            //todo: swap this to redsvd and use control file truncation limits
			Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(constr_jco);
			double threshold = cod.threshold();
			int effective_rank = cod.rank();

			Eigen::MatrixXd U = cod.matrixQ();
			Eigen::MatrixXd V = cod.matrixZ();
			Eigen::MatrixXd T = cod.pseudoInverse();

			constr_jco = U.leftCols(effective_rank);
		}
		
		if ((constraint_diff.array().abs() > filter.get_viol_tol()).any())  // todo make some level of forgiveness with a tolerance parameter here
		{
			ss.str("");
			ss << "WARNING: not on constraint but working set not empty" << endl;
			ss << "current constraint diffs: " << constraint_diff.transpose() << endl;
			ss << "continuing...";
			performance_log->log_event(ss.str());
			frec << ss.str() << endl;
		}
		Eigen::MatrixXd G = *hessian.e_ptr();  // initialize the hessian as identity matrix (Dhedari et al 2012)
		Eigen::VectorXd c;
		//c is the rhs of Eq 18.20, p. 538 in Nocedal and Wright

		string eqp_solve_method; // probably too heavy to be a ++arg
		eqp_solve_method = "null_space";
		if (eqp_solve_method == "null_space")
		{
			x = _kkt_null_space(G, constr_jco, constraint_diff, grad_vector);
			search_d = x.first;
			lm = x.second;
		}
		else if (eqp_solve_method == "direct")
		{
			x = _kkt_direct(G, constr_jco, constraint_diff, c, Cnames);
			search_d = x.first;
			lm = x.second;
		}
		else // if "schur", "cg", ...
		{
			throw_sqp_error("eqp_solve_method not implemented");
		}
		lambda = lm; //TODO: refactor for ensemble
	}
	else  // solve unconstrained QP subproblem
	{
		Eigen::MatrixXd H = *hessian.e_ptr();
		search_d = H.ldlt().solve(-grad_vector); //Eq 7.9, Nocedal and Wright p. 169
		double dir_dot_grad = search_d.dot(grad_vector);
		if (dir_dot_grad > 0) {
			search_d = -grad_vector;
		}
		lm = Eigen::VectorXd::Zero(0);

	}
	return pair<Eigen::VectorXd, Eigen::VectorXd> (search_d, lm);
}

Eigen::VectorXd SeqQuadProgram::fancy_solve_routine(const Parameters& _current_dv_num_values, const Parameters& _grad_vector)
{

	// grad vector computation
	Eigen::VectorXd grad = _grad_vector.get_data_eigen_vec(dv_names);

	// search direction computation; RQM: this is for the QP subproblem assuming inqeuality constraints are active
	//Eigen::VectorXd search_d = calc_search_direction_vector(_current_dv_num_values, grad);  
	//pair<Eigen::VectorXd, Eigen::VectorXd> x = calc_search_direction_vector(_current_dv_num_values, grad);

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


	//return x.first;  // search_d;
	return grad;
}

bool SeqQuadProgram::solve_new()
{
	//stringstream ss;
	//ofstream& frec = file_manager.rec_ofstream();
	//if ((use_ensemble_grad) && (dv.shape().first <= error_min_reals))
	//{
	//	message(0, "too few active realizations:", oe.shape().first);
	//	message(1, "need more than ", error_min_reals);
	//	throw_sqp_error(string("too few active realizations, cannot continue"));
	//}
	//else if ((use_ensemble_grad) && (dv.shape().first < warn_min_reals))
	//{
	//	ss.str("");
	//	ss << "WARNING: less than " << warn_min_reals << " active realizations...might not be enough";
	//	string s = ss.str();
	//	message(1, s);
	//}

	//Parameters _current_num_dv_values = current_ctl_dv_values;  // make copy
	//ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();
	//pts.ctl2numeric_ip(_current_num_dv_values); 
	//_current_num_dv_values = _current_num_dv_values.get_subset(dv_names.begin(), dv_names.end()); 

	//pair<Mat, bool> constraint_mat = get_constraint_mat(current_ctl_dv_values, current_obs,(working_set_tol));
	//current_constraint_mat = constraint_mat.first;
	//cnames = constraint_mat.first.get_row_names();

	////copy for BFGS later
	//prev_ctl_dv_values = current_ctl_dv_values; 
	//prev_constraint_mat = current_constraint_mat;

	//Eigen::VectorXd search_d, lm;
	//Eigen::VectorXd grad = current_grad_vector.get_data_eigen_vec(dv_names);
	//bool successful = false;
	//Covariance old_hessian = hessian;

	//int line_search_attempts = 0;
	//while (!successful && line_search_attempts < max_line_search_attempts)
	//{
	//	constraint_jco = constraint_mat.first.e_ptr()->toDense();
	//	infeas_cand_obs.clear();
	//	infeas_cand_dv_values.clear();
	//	pair<Eigen::VectorXd, Eigen::VectorXd> x = calc_search_direction_vector(_current_num_dv_values, grad);
	//	search_d = x.first;
	//	lm = x.second;

	//	message(1, "constraint_jco:", constraint_jco); // tmp
	//	message(1, "sd:", search_d.transpose());  // tmp
	//	message(1, "sd_norm:", search_d.norm()); //tmp
	//	message(1, "lm:", lm); //tmp

	//	if (cnames.size() > 0)
	//	{
	//		//Algorithm 16.3 in Nocedal and Wright, pp. 472-473
	//		bool search_d_approx_zero = false;
	//		double tol = 0.0001;  // should be a carefully chosen tolerance
	//		if (search_d.norm() < tol) {
	//			search_d_approx_zero = true;
	//		}

	//		if (search_d_approx_zero)
	//		{
	//			pair<Mat, bool> constraint_mat = get_constraint_mat(current_ctl_dv_values, current_obs, (working_set_tol), &lm);

	//			//if all multipliers non-negative, we're at optimal solution
	//			if (constraint_mat.second) {
	//				message(1, "optimal solution found - all Lagrange multipliers non-negative");
	//				converged = true;
	//				return true;
	//			}
	//		}

	//		if ((lm.array() > 0).any())
	//		{
	//			pair<Mat, bool> constraint_mat = get_constraint_mat(current_ctl_dv_values, current_obs, (working_set_tol), &lm);
	//			if (constraint_mat.first.get_row_names() != cnames) {
	//				current_constraint_mat = constraint_mat.first;
	//				cnames = constraint_mat.first.get_row_names();
	//				message(1, "constraints dropped due to negative multipliers:", cnames);
	//				continue;
	//			}
	//		}
	//	}

	//	//trial_ctl_dv_values = current_ctl_dv_values;
	//	//trial_obs = current_obs;
	//	//successful = trust_region_step(current_ctl_dv_values, grad); //should consider switching to trust region at some point?
	//	is_blocking_constraint = false;
	//	successful = line_search(search_d, _current_num_dv_values, grad);
	//	string blocking_constraint = "";
	//	if (successful)
	//	{
	//		if (is_blocking_constraint)
	//		{
	//			map<string, double> violations;
	//			if (infeas_cand_obs.size() != 0)
	//				violations = constraints.get_unsatified_obs_constraints(infeas_cand_obs, 0.0);

	//			if (!violations.empty())
	//			{
	//				double max_violation = -1.0;
	//				for (const auto& v : violations)
	//				{
	//					if (v.second > max_violation)
	//					{
	//						max_violation = v.second;
	//						blocking_constraint = v.first;
	//					}
	//				}

	//				if (find(cnames.begin(), cnames.end(), blocking_constraint) == cnames.end())
	//				{
	//					message(1, "adding blocking constraint to working set:", blocking_constraint);
	//					cnames.push_back(blocking_constraint);
	//					pair<Mat, bool> new_constraint_mat = get_constraint_mat(current_ctl_dv_values, current_obs);
	//					current_constraint_mat = new_constraint_mat.first;
	//				}
	//				else
	//					blocking_constraint = "";
	//			}
	//		}
	//	}

	//	
	//	if (blocking_constraint != "")
	//	{
	//		successful = false;
	//		message(1, "performing binary search for constraint boundary with working set:", cnames);
	//		if (pest_scenario.get_pestpp_options().get_sqp_solve_partial_step())
	//		{
	//			if (iterative_partial_step(blocking_constraint))
	//			{
	//				constraint_mat = get_constraint_mat(current_ctl_dv_values, current_obs, (working_set_tol));
	//				current_constraint_mat = constraint_mat.first;
	//				successful = true;
	//				BASE_SCALE_FACTOR *= SF_INC_FAC;
	//				break;
	//			}
	//			else
	//				throw_sqp_error("Something is wrong with iterative partial step...");
	//		}
	//	}

	//	if (successful)
	//	{
	//		if ((search_d.norm() < 0.1) && (((lm.array() < 0).all()) || (lm.size() != 0)))
	//			BASE_SCALE_FACTOR /= SF_INC_FAC;
	//		else
	//			BASE_SCALE_FACTOR *= SF_INC_FAC;
	//	}
	//	else
	//	{
	//		n_consec_failures++;
	//		line_search_attempts++;

	//		if (use_ensemble_grad)
	//			BASE_SCALE_FACTOR /= SF_INC_FAC;
	//		else
	//			BASE_SCALE_FACTOR *= SF_DEC_FAC;

	//		if (n_consec_failures >= max_consec_failures)
	//		{
	//			Eigen::SparseMatrix<double> h(dv_names.size(), dv_names.size());
	//			h.setIdentity();
	//			update_scaling(search_d, grad);
	//			for (int i = 0; i < dv_names.size(); i++) {
	//				h.coeffRef(i, i) *= diagonal_scaling(i);
	//			}

	//			hessian = Covariance(dv_names, h);
	//			n_consec_failures = 0;
	//			BASE_SCALE_FACTOR *= SF_INC_FAC;
	//		}
	//			if (n_consec_failures >= max_consec_failures)
	//				break;
	//	}

	//}
	//message(0, "new base scale factor: ", BASE_SCALE_FACTOR);
	//return successful;

	return true;
}

bool SeqQuadProgram::solve_new_ensemble()
{
	stringstream ss;
	ofstream& frec = file_manager.rec_ofstream();

	prev_ctl_dv_values = current_ctl_dv_values;
	pair<bool, pair<Parameters, Observations>> first_pick = pick_upgrade_and_update_current(dv, oe);
	if (first_pick.first)
	{
		current_ctl_dv_values = first_pick.second.first;
		current_obs = first_pick.second.second;
	}

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

	performance_log->log_event("reordering variables in dv");
	dv.reorder(vector<string>(), dv_names);
	dv.transform_ip(ParameterEnsemble::transStatus::NUM);
	ParameterEnsemble _dvs = dv;
	_dvs.drop_rows(vector<string>{"BASE"}, true);

	//use subset to determine optimal search direction
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

	ParameterEnsemble _avail_dvs = _dvs, _drawn_dvs(&pest_scenario, &rand_gen);
		
	if (!sampling_tracking_initialized)
	{
		unselected_dv_indices.clear();
		selected_dv_indices.clear();
		for (int i = 0; i < _dvs.shape().first; i++)
			unselected_dv_indices.insert(i);
		
		sampling_tracking_initialized = true;
	}

	string how = pest_scenario.get_pestpp_options().get_ies_subset_how();

	if (unselected_dv_indices.empty() || local_subset_size > unselected_dv_indices.size())
	{
		unselected_dv_indices.clear();
		selected_dv_indices.clear();
		for (int i = 0; i < _dvs.shape().first; i++)
			unselected_dv_indices.insert(i);
	}

	vector<int> subset_idxs = get_subset_idxs(_dvs.shape().first, local_subset_size);

	for (int idx : subset_idxs)
	{
		unselected_dv_indices.erase(idx);
		selected_dv_indices.insert(idx);

		string par_name = _dvs.get_real_names()[idx];

		ParameterEnsemble t;
		t.reserve(vector<string>{ par_name }, pest_scenario.get_ctl_ordered_par_names());
		t.add_2_row_ip(par_name, _dvs.get_real_vector(par_name));

		if (_drawn_dvs.shape().first == 0) 
			_drawn_dvs = t;
		else
			_drawn_dvs.append_other_rows(t, true);
		
		_avail_dvs.drop_rows({ par_name }, true);
	}

	Covariance old_hessian = hessian;
	map <string, pair<Mat, bool>> constraint_mat_en;
	map<string, vector<string>> cnames_en;
	map<string, Eigen::VectorXd> search_d_en, lm_en;
	map<string, Eigen::MatrixXd> constraint_jco_en;
	Parameters dv_vals = current_ctl_dv_values;
	Observations obs_vals = current_obs;

	Eigen::VectorXd grad = current_grad_vector.get_data_eigen_vec(dv_names);

	for (auto d : dv.get_real_names())
	{
		Eigen::VectorXd real_dv_vec = dv.get_real_vector(d);
		dv_vals.update_without_clear(dv_names, real_dv_vec);
		Eigen::VectorXd real_obs_vec = oe.get_real_vector(d);
		obs_vals.update_without_clear(oe.get_var_names(), real_obs_vec);
		
		constraint_mat_en[d] = get_constraint_mat(dv_vals, obs_vals, working_set_tol);
		Mat current_cmat = constraint_mat_en[d].first;
		cnames_en[d] = constraint_mat_en[d].first.get_row_names();
		constraint_jco_en[d] = constraint_mat_en[d].first.e_ptr()->toDense();

		ss.str("");
		ss << "starting search direction calcs for realization " << d << " with working set: " << cnames_en[d];
		performance_log->log_event(ss.str());
		frec << ss.str() << endl;

		pair<Eigen::VectorXd, Eigen::VectorXd> x = calc_search_direction_vector(dv_vals, obs_vals, grad, &constraint_jco_en[d], &cnames_en[d]);
		search_d_en[d] = x.first;
		lm_en[d] = x.second;

		if ((lm_en[d].array() > 0).any())
		{
			constraint_mat_en[d] = get_constraint_mat(dv_vals, obs_vals, working_set_tol, &lm_en[d]);
			if (constraint_mat_en[d].first.get_row_names() != cnames_en[d]) 
			{
				//current_constraint_mat = constraint_mat.first;
				vector<string> prev_cnames = cnames_en[d];
				cnames_en[d] = constraint_mat_en[d].first.get_row_names();
				vector<string> dropped_cnames;
				for (auto c : prev_cnames)
				{
					if (find(cnames_en[d].begin(), cnames_en[d].end(), c) == cnames_en[d].end())
						dropped_cnames.push_back(c);
				}
				ss.str("");
				ss << "constraints dropped for realization " << d << " due to non-negative multipliers: " << dropped_cnames;
				ss << "recalculating search_d with new working set: " << cnames_en[d];
				frec << ss.str();

				constraint_jco_en[d] = constraint_mat_en[d].first.e_ptr()->toDense();
				pair<Eigen::VectorXd, Eigen::VectorXd> x = calc_search_direction_vector(dv_vals, obs_vals, grad, &constraint_jco_en[d] , &cnames_en[d]);
				search_d_en[d] = x.first;
				lm_en[d] = x.second;
			}
		}
	}

	//message(1, "constraint_jco:", constraint_jco); // tmp
	//message(1, "sd:", search_d_en["BASE"].transpose());  // tmp
	//message(1, "sd_norm:", search_d.norm()); //tmp
	//message(1, "lm:", lm_en["BASE"]); //tmp


	//pair<Mat, bool> constraint_mat = get_constraint_mat(current_ctl_dv_values, current_obs, (working_set_tol));
	constraint_jco = constraint_mat_en["BASE"].first.e_ptr()->toDense();
	current_constraint_mat = constraint_mat_en["BASE"].first;
	cnames = constraint_mat_en["BASE"].first.get_row_names();

	
	prev_constraint_mat = current_constraint_mat;

	Eigen::VectorXd search_d, lm;
	bool successful = false;
	int line_search_attempts = 0;
	while (!successful && line_search_attempts < max_line_search_attempts)
	{
		infeas_cand_obs.clear();
		infeas_cand_dv_values.clear();

		is_blocking_constraint = false;
		successful = line_search(search_d_en, grad, &_drawn_dvs, first_pick.first);
		string blocking_constraint = "";
		/*if (successful)
		{
			if (is_blocking_constraint && !use_ensemble_grad)
			{
				map<string, double> violations;
				if (infeas_cand_obs.size() != 0)
					violations = constraints.get_unsatified_obs_constraints(infeas_cand_obs, 0.0);

				if (!violations.empty())
				{
					double max_violation = -1.0;
					for (const auto& v : violations)
					{
						if (v.second > max_violation)
						{
							max_violation = v.second;
							blocking_constraint = v.first;
						}
					}

					if (find(cnames.begin(), cnames.end(), blocking_constraint) == cnames.end())
					{
						message(1, "adding blocking constraint to working set:", blocking_constraint);
						cnames.push_back(blocking_constraint);
						pair<Mat, bool> new_constraint_mat = get_constraint_mat(current_ctl_dv_values, current_obs);
						current_constraint_mat = new_constraint_mat.first;
					}
					else
						blocking_constraint = "";
				}
			}
		}*/


		//if (blocking_constraint != "" && !use_ensemble_grad)
		//{
		//	successful = false;
		//	message(1, "performing binary search for constraint boundary with working set:", cnames);
		//	if (pest_scenario.get_pestpp_options().get_sqp_solve_partial_step())
		//	{
		//		if (iterative_partial_step(blocking_constraint))
		//		{
		//			constraint_mat = get_constraint_mat(current_ctl_dv_values, current_obs, (working_set_tol));
		//			current_constraint_mat = constraint_mat.first;
		//			successful = true;
		//			BASE_SCALE_FACTOR *= SF_INC_FAC;
		//			break;
		//		}
		//		else
		//			throw_sqp_error("Something is wrong with iterative partial step...");
		//	}
		//}


		if (successful)
		{
			if ((search_d.norm() < 0.1) && (((lm.array() < 0).all()) || (lm.size() != 0)))
				BASE_SCALE_FACTOR /= SF_INC_FAC;
			else
				BASE_SCALE_FACTOR *= SF_INC_FAC;
		}
		else
		{
			if (first_pick.first)
			{
				ss.str("");
				ss << "line search failed. centgering new ensemble at best member of the current ensemble: ";
				message(1, ss.str());
				
				successful = true;
				if ((search_d.norm() < 0.1) && (((lm.array() < 0).all()) || (lm.size() != 0)))
					BASE_SCALE_FACTOR /= SF_INC_FAC;
				else
					BASE_SCALE_FACTOR *= SF_INC_FAC;

				break;
			}
			if (use_cmaes)
			{
				message(1, "line search failed. updating covariance from current ensemble and regenerating ensemble via CMA-ES...");
				break;
			}
			n_consec_failures++;
			line_search_attempts++;

			if (use_ensemble_grad)
				BASE_SCALE_FACTOR /= SF_INC_FAC;
			else
				BASE_SCALE_FACTOR *= SF_DEC_FAC;

			if (n_consec_failures >= max_consec_failures)
			{
				Eigen::SparseMatrix<double> h(dv_names.size(), dv_names.size());
				h.setIdentity();
				update_scaling(search_d_en["BASE"], grad);
				for (int i = 0; i < dv_names.size(); i++) {
					h.coeffRef(i, i) *= diagonal_scaling(i);
				}

				hessian = Covariance(dv_names, h);
				n_consec_failures = 0;
				BASE_SCALE_FACTOR *= SF_INC_FAC;
			}
			if (n_consec_failures >= max_consec_failures)
				break;
		}

	}
	message(0, "new base scale factor: ", BASE_SCALE_FACTOR);
	return successful;
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

pair<bool, pair<Parameters, Observations>> SeqQuadProgram::pick_upgrade_and_update_current(ParameterEnsemble& dv_candidates, ObservationEnsemble& _oe)
{
	message(0, " current best phi:", last_best);
	stringstream ss;
	Eigen::VectorXd obj_vec = get_obj_vector(dv_candidates, _oe);
	double oext, oviol = 0.0, nviol = 0.0;
	if (obj_sense == "minimize")
		oext = numeric_limits<double>::max();
	else
		oext = numeric_limits<double>::min();
	int idx = -1;
	vector<string> real_names = dv_candidates.get_real_names();
	map<string, double> total_viol_map;
	map<string, double> obj_map = get_obj_map(dv_candidates, _oe);
	map<string, map<string, double>> violations = constraints.get_ensemble_violations_map(dv_candidates, _oe, filter.get_viol_tol(), true);
	map<string, map<string, double>> violations_nominal = constraints.get_ensemble_violations_map(dv_candidates, _oe, 0.0, true);
	map<string, map<string, double>> feasible_distance = constraints.get_ensemble_violations_map(dv_candidates, _oe, -1, true);
	Parameters cand_dv_values = current_ctl_dv_values;
	Observations cand_obs_values = current_obs;
	Eigen::VectorXd t;
	vector<string> onames = _oe.get_var_names();
	ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();
	bool filter_accept;
	string tag;
	vector<double> infeas_vec, nviol_vec, feas_dist_map;
	vector<int> accept_idxs;
	int num_feas_filter_accepts = 0, num_infeas_filter_accepts = 0;
	bool accept = false;

	for (int i = 0; i < obj_vec.size(); i++)
	{
		ss.str("");
		ss << "candidate: " << real_names[i];
		tag = ss.str();
		ss.str("");
		ss << tag << " phi: " << obj_vec[i];
		double infeas_sum = 0.0, infeas_sum_nom = 0.0, feas_dist = 0;
		for (auto& v : violations[real_names[i]])
			infeas_sum += v.second;
		ss << " infeasibility total: " << infeas_sum << ", ";
		infeas_vec.push_back(infeas_sum);
		for (auto& v : violations_nominal[real_names[i]])
			infeas_sum_nom += v.second;
		nviol_vec.push_back(infeas_sum_nom);
		total_viol_map[real_names[i]] = infeas_sum_nom;
		for (auto& f : feasible_distance[real_names[i]])
			feas_dist -= f.second;
		feas_dist_map.push_back(feas_dist);
		filter_accept = filter.accept(obj_vec[i], infeas_sum, iter);
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
			else {
				accept_idxs.push_back(i);
			}
		}

		if (obj_vec[i] < oext)
		{
			idx = i;
			oext = obj_vec[idx];
			oviol = infeas_vec[idx];
		}

		t = dv_candidates.get_real_vector(real_names[i]);
		cand_dv_values = current_ctl_dv_values;
		cand_dv_values.update_without_clear(dv_names, t);
		pts.numeric2ctl_ip(cand_dv_values);
		t = _oe.get_real_vector(real_names[i]);
		cand_obs_values.update_without_clear(onames, t);
		constraints.sqp_report(iter, cand_dv_values, cand_obs_values, false, tag);
	}

	cmaes.update_archives(dv_candidates, obj_map, total_viol_map, iter, true);

	Parameters p;
	Observations o;
	if (accept_idxs.size() > 0)
	{
		ss.str("");
		ss << "number of realizations passing first filter:" << accept_idxs.size();
		message(1, ss.str());
		message(1, "accepting realization ", real_names[idx]);

		t = dv_candidates.get_real_vector(idx);
		vector<string> vnames = dv_candidates.get_var_names();

		p.update_without_clear(vnames, t);

		t = _oe.get_real_vector(idx);

		o.update_without_clear(onames, t);
		//current_obs.update_without_clear(onames, t);
		last_best = oext;
		last_viol = oviol;
		message(0, "new best phi and infeas:", vector<double>{last_best, last_viol});
		//best_phis.push_back(oext);
		//best_violations.push_back(oviol);
		//if (last_viol == 0)
		//	best_feas_phis.push_back(oext);

		filter.update(oext, infeas_vec[idx], iter);

		return pair <bool, pair<Parameters, Observations>>(true, pair<Parameters, Observations>(p, o));
	}
	else
	{
		message(1, "no realization passed first filter. continuing to line search...");
		return pair <bool, pair<Parameters, Observations>>(false, pair<Parameters, Observations>(p, o));
	}
}

bool SeqQuadProgram::pick_candidate_and_update_current(ParameterEnsemble& dv_candidates, ObservationEnsemble& _oe, map<string,double>& sf_map, bool first_pick)
{
	message(0, "selecting candidate via filter. current best phi:", last_best);

	stringstream ss;
	Eigen::VectorXd obj_vec = get_obj_vector(dv_candidates, _oe);
	double oext,oviol = 0.0, nviol = 0.0, lviol = numeric_limits<double>::max();
	if (obj_sense == "minimize")
		oext = numeric_limits<double>::max();
	else
		oext = numeric_limits<double>::min();
	int idx = -1, jdx = -1;
	vector<string> real_names = dv_candidates.get_real_names();
	map<string, double> total_viol_map;
	map<string, double> obj_map = get_obj_map(dv_candidates, _oe);
	map<string, map<string, double>> violations = constraints.get_ensemble_violations_map(dv_candidates,_oe,filter.get_viol_tol(),true);
	map<string, map<string, double>> violations_nominal = constraints.get_ensemble_violations_map(dv_candidates, _oe, 0.0, true);
	map<string, map<string, double>> feasible_distance = constraints.get_ensemble_violations_map(dv_candidates, _oe, -1, true);
	Parameters cand_dv_values = current_ctl_dv_values;
	Observations cand_obs_values = current_obs;
	Eigen::VectorXd t;
	vector<string> onames = _oe.get_var_names();
	ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();
	bool filter_accept;
	string tag;
	vector<double> infeas_vec, nviol_vec, feas_dist_map;
	vector<int> accept_idxs;
	int num_feas_filter_accepts = 0, num_infeas_filter_accepts = 0;
	bool accept = false;

	for (int i = 0; i < obj_vec.size(); i++)
	{
		ss.str("");
		ss << "candidate: " << real_names[i] << ", scale factor: " << sf_map.at(real_names[i]);
		tag = ss.str();
		ss.str("");
		ss << tag << " phi: " << obj_vec[i];
		double infeas_sum = 0.0, infeas_sum_nom = 0.0, feas_dist = 0;
		for (auto& v : violations[real_names[i]])
			infeas_sum += v.second;
		ss << " infeasibility total: " << infeas_sum << ", ";
		infeas_vec.push_back(infeas_sum);
		for (auto& v : violations_nominal[real_names[i]])
			infeas_sum_nom += v.second;
		nviol_vec.push_back(infeas_sum_nom);
		total_viol_map[real_names[i]] = infeas_sum_nom;
		for (auto& f : feasible_distance[real_names[i]])
			feas_dist -= f.second;
		feas_dist_map.push_back(feas_dist);
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
        }

		if ((!use_ensemble_grad) && (obj_vec[i] < oext * 1.25))
		{
			idx = i;
			oext = obj_vec[idx];
			if (infeas_sum > 0)
				is_blocking_constraint = true;
		}

		t = dv_candidates.get_real_vector(real_names[i]);
		cand_dv_values = current_ctl_dv_values;
		cand_dv_values.update_without_clear(dv_names, t);
		pts.numeric2ctl_ip(cand_dv_values);
		t = _oe.get_real_vector(real_names[i]);
		cand_obs_values.update_without_clear(onames, t);
		constraints.sqp_report(iter, cand_dv_values, cand_obs_values, false,tag);
	}

	cmaes.update_archives(dv_candidates, obj_map, total_viol_map, iter, false);

	if (obj_sense == "minimize")
		oext = numeric_limits<double>::max();
	else
		oext = numeric_limits<double>::min();
    if (accept_idxs.size() > 0)
    {
        accept = true;
        ss.str("");
        ss << "number of scale factors passing filter:" << accept_idxs.size();
        message(1,ss.str());
		double delta_feas = -999;
		for (auto iidx : accept_idxs)
        {
			if ((accept_idxs.size() > 1) && (iidx > 0))
				delta_feas = feas_dist_map[iidx] * feas_dist_map[iidx - 1];
			else
				delta_feas = feas_dist_map[iidx];

			if (obj_vec[iidx] < oext)
            {
				if ((nviol_vec[iidx] <= 1E-6) && (delta_feas >= 0))
				{
					idx = iidx;
					oext = obj_vec[iidx];
					oviol = infeas_vec[iidx];
					nviol = nviol_vec[iidx];
				}
				if (infeas_vec[iidx] <= 1E-6)
					num_feas_filter_accepts++;
            }
			if (nviol_vec[iidx] > 1E-6)
			{
				if (nviol_vec[iidx] < lviol)
				{
					jdx = iidx;
					lviol = nviol_vec[iidx];
				}
			}
        }

		if (jdx == -1)
		{
			for (int i = idx; i < obj_vec.size(); i++)
			{
				if (nviol_vec[i] > 1E-6)
				{
					if (nviol_vec[i] < lviol)
					{
						jdx = i;
						lviol = nviol_vec[i];
					}
				}
			}
		}
		
		if (idx != -1)
			filter.update(oext,infeas_vec[idx],iter,sf_map.at(real_names.at(idx)));
    }
	else
    {
	    message(0,"filter failed, checking for feasible solutions....");
	    double viol_tol = filter.get_viol_tol();
		nviol_vec.clear();
		double infeas_sum_nom = 0.0;
		for (int i = 0; i < obj_vec.size(); i++)
		{
			violations_nominal = constraints.get_ensemble_violations_map(dv_candidates, _oe, -1.0, true);
			for (auto& v : violations_nominal[real_names[i]])
				infeas_sum_nom += v.second;
			nviol_vec.push_back(infeas_sum_nom);
		}

        for (int i=0;i<nviol_vec.size();i++)
        {
            if (nviol_vec[i] < viol_tol)
            {
                if ((obj_sense=="minimize") && (obj_vec[i] < oext))
                {
                    idx = i;
                    oext = obj_vec[i];
                    oviol = infeas_vec[i];
					nviol = nviol_vec[i];
					num_feas_filter_accepts++;
                }
                else if ((obj_sense == "maximize") && (obj_vec[i] > oext))
                {
                    idx = i;
                    oext = obj_vec[i];
                    oviol = infeas_vec[i];
					nviol = nviol_vec[i];
					num_feas_filter_accepts++;
                }
            }
        }

		if (idx != -1)
			filter.update(oext, infeas_vec[idx], iter, sf_map.at(real_names.at(idx)));
		/*if (idx == -1)
		{
			message(0, "no feasible solutions, choosing lowest constraint violation...");
			double viol_min = 1e+300;
			for (int i = 0;i < infeas_vec.size();i++)
			{
				if (nviol_vec[i] < viol_min)
				{
					viol_min = infeas_vec[i];
					idx = i;
					jdx = i;
					oext = obj_vec[i];
					oviol = infeas_vec[i];
				}
			}
		}
		else
			accept = true;*/

		//message(0, "no feasible solutions, choosing lowest constraint violation...");
		//double viol_min = 1e+300;
		//for (int i = 0;i < infeas_vec.size();i++)
		//{
		//	if (nviol_vec[i] < viol_min)
		//	{
		//		viol_min = infeas_vec[i];
		//		idx = i;
		//		jdx = i;
		//		oext = obj_vec[i];
		//		oviol = infeas_vec[i];
		//		num_feas_filter_accepts++;
		//	}
		//}
		//if (idx != -1)
		//	accept = true;

    }
	
	if (jdx != -1)
	{
		infeas_cand_dv_values = pest_scenario.get_ctl_parameters();
		t = dv_candidates.get_real_vector(jdx);
		vector<string> vnames = dv_candidates.get_var_names();
		Parameters p;
		p.update_without_clear(vnames, t);
		pts.numeric2ctl_ip(p);

		for (auto& d : dv_names)
			infeas_cand_dv_values[d] = p[d];
		
		infeas_cand_obs = pest_scenario.get_ctl_observations();
		t = _oe.get_real_vector(jdx);
		infeas_cand_obs.update_without_clear(onames, t);
	}

	if (accept)
	{
		t = dv_candidates.get_real_vector(real_names[idx]);
		cand_dv_values = current_ctl_dv_values;
		cand_dv_values.update_without_clear(dv_names, t);
		pts.numeric2ctl_ip(cand_dv_values);
		t = _oe.get_real_vector(real_names[idx]);
		cand_obs_values.update_without_clear(onames, t);

		ss.str("");
		ss << "best candidate (scale factor: " << setprecision(4) << sf_map.at(real_names[idx]) << ", phi: " << oext << ", infeas: " << oviol << ")";
		constraints.sqp_report(iter, cand_dv_values, cand_obs_values, true, ss.str());
		filter.report(file_manager.rec_ofstream(), iter);

		if (num_feas_filter_accepts == 0)
		{
			BASE_SCALE_FACTOR = BASE_SCALE_FACTOR * SF_DEC_FAC;
			return false;
		}
		if (infeas_vec[idx] > 1E-6)
			n_consec_infeas++;
		else {
			n_consec_infeas = 0;
			if (use_ensemble_grad) {
				
				double new_par_sigma = pest_scenario.get_pestpp_options().get_par_sigma_range();
				new_par_sigma = new_par_sigma * pow(PAR_SIGMA_INC_FAC, 1 / (static_cast<double>(iter)));
				new_par_sigma = min(new_par_sigma, par_sigma_max);

				message(1, "increasing par_sigma_range to", new_par_sigma);
				message(1, "regenerating parcov");
				pest_scenario.get_pestpp_options_ptr()->set_par_sigma_range(new_par_sigma);
				parcov.try_from(pest_scenario, file_manager);
			}
		}

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
		if (last_viol == 0)
			best_feas_phis.push_back(oext);

		return true;
		
	}
	/*else if (!first_pick)
	{
		if (idx != -1)
		{
			t = dv_candidates.get_real_vector(real_names[idx]);
			cand_dv_values = current_ctl_dv_values;
			cand_dv_values.update_without_clear(dv_names, t);
			pts.numeric2ctl_ip(cand_dv_values);
			t = _oe.get_real_vector(real_names[idx]);
			cand_obs_values.update_without_clear(onames, t);
			ss.str("");
			ss << "best candidate (scale factor: " << setprecision(4) << sf_map.at(real_names[idx]) << ", phi: " << oext << ", infeas: " << infeas_vec[idx] << ")";
			constraints.sqp_report(iter, cand_dv_values, cand_obs_values, false, ss.str());

			message(0, "temporarily accepting upgrade");
			t = dv_candidates.get_real_vector(idx);
			vector<string> vnames = dv_candidates.get_var_names();
			Parameters p;
			p.update_without_clear(vnames, t);

			pts.numeric2ctl_ip(p);

			for (auto& d : dv_names)
				current_ctl_dv_values[d] = p[d];
			t = _oe.get_real_vector(idx);
			current_obs.update_without_clear(onames, t);

			double new_par_sigma = pest_scenario.get_pestpp_options().get_par_sigma_range();
			new_par_sigma = new_par_sigma * PAR_SIGMA_DEC_FAC;
			new_par_sigma = max(new_par_sigma, par_sigma_min);
			message(1, "decreasing par_sigma_range to", new_par_sigma);
			message(1, "regenerating parcov");
			parcov.try_from(pest_scenario, file_manager);
			pest_scenario.get_pestpp_options_ptr()->set_par_sigma_range(new_par_sigma);

		}
		
		return true;
	}*/
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
            cout << "parcov: " << endl << parcov << endl;
        }
        //BASE_SCALE_FACTOR = BASE_SCALE_FACTOR * SF_DEC_FAC;
        
		return false;
	}
}

bool SeqQuadProgram::pick_partial_step(ParameterEnsemble& dv_candidates, ObservationEnsemble& _oe, map<string, double>& sf_map)
{
	SqpFilter partial_step_filter = filter;
	partial_step_filter.set_tol(0.0);
	message(0, " current best phi:", last_best);
	stringstream ss;
	Eigen::VectorXd obj_vec = get_obj_vector(dv_candidates, _oe);
	double oext, oviol = 0.0;
	if (obj_sense == "minimize")
		oext = numeric_limits<double>::max();
	else
		oext = numeric_limits<double>::min();
	int idx = -1;
	vector<string> real_names = dv_candidates.get_real_names();
	map<string, double> obj_map = get_obj_map(dv_candidates, _oe);
	map<string, map<string, double>> violations = constraints.get_ensemble_violations_map(dv_candidates, _oe, partial_step_filter.get_viol_tol(), true);
	Parameters cand_dv_values = current_ctl_dv_values;
	Observations cand_obs_values = current_obs;
	Eigen::VectorXd t;
	vector<string> onames = _oe.get_var_names();
	ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();
	bool filter_accept;
	string tag;
	vector<double> infeas_vec;
	vector<int> accept_idxs;
	bool accept = false;
	for (int i = 0; i < obj_vec.size(); i++)
	{
		ss.str("");
		ss << "partial step: " << real_names[i] << ", scale factor: " << sf_map.at(real_names[i]);
		tag = ss.str();
		ss.str("");
		ss << tag << " phi: " << obj_vec[i];
		double infeas_sum = 0.0;
		for (auto& v : violations[real_names[i]])
			infeas_sum += v.second;
		ss << " infeasibilty total: " << infeas_sum << ", ";
		infeas_vec.push_back(infeas_sum);
		filter_accept = partial_step_filter.accept(obj_vec[i], infeas_sum, iter, sf_map.at(real_names[i]), false);
		if (filter_accept)
		{
			ss << " filter accepted ";
			accept_idxs.push_back(i);
		}
		else
			ss << " filter rejected ";
		message(1, ss.str());
		if (filter_accept)
			
			
		t = dv_candidates.get_real_vector(real_names[i]);
		cand_dv_values = current_ctl_dv_values;
		cand_dv_values.update_without_clear(dv_names, t);
		pts.numeric2ctl_ip(cand_dv_values);
		t = _oe.get_real_vector(real_names[i]);
		cand_obs_values.update_without_clear(onames, t);
		constraints.sqp_report(iter, cand_dv_values, cand_obs_values, false, tag);
	}

	if (accept_idxs.size() > 0)
	{
		accept = true;
		ss.str("");
		ss << "number of partial steps passing filter:" << accept_idxs.size();
		message(1, ss.str());

		for (auto iidx : accept_idxs)
		{
			if (obj_vec[iidx] < last_best)
			{
				if (infeas_vec[iidx] <= 1E-6)
				{
					idx = iidx;
					oext = obj_vec[iidx];
					oviol = infeas_vec[iidx];
				}
			}
		}
		partial_step_filter.update(oext, infeas_vec[idx], iter, sf_map.at(real_names.at(idx)));
	}

	if (accept)
	{
		message(0, "best phi and infeas this iteration: ", vector<double>{oext, oviol});
		t = dv_candidates.get_real_vector(real_names[idx]);
		cand_dv_values = current_ctl_dv_values;
		cand_dv_values.update_without_clear(dv_names, t);
		pts.numeric2ctl_ip(cand_dv_values);
		t = _oe.get_real_vector(real_names[idx]);
		cand_obs_values.update_without_clear(onames, t);
		ss.str("");
		ss << "partial step closest to constraint (scale factor: " << setprecision(4) << sf_map.at(real_names[idx]) << ", phi: " << oext << ", infeas: " << oviol << ")";
		constraints.sqp_report(iter, cand_dv_values, cand_obs_values, true, ss.str());

		message(0, "accepting partial step", real_names[idx]);
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
		message(0, "new best phi and infeas:", vector<double>{last_best, last_viol});
		best_phis.push_back(oext);
		best_violations.push_back(oviol);
		if (last_viol == 0)
			best_feas_phis.push_back(oext);
		return true;
	}
	else
	{
		message(0, "not accepting partial step. Current position is close enough to the constraint boundary.");
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

vector<int> SeqQuadProgram::get_subset_idxs(int size, int nreal_subset)
{
	vector<int> subset_idxs;
	if ((!use_subset) || (nreal_subset >= size))
	{
		for (int i = 0; i < size; i++)
			subset_idxs.push_back(i);
		return subset_idxs;
	}

	vector<string>::iterator bidx = find(dv_names.begin(), dv_names.end(), BASE_REAL_NAME);
	if (bidx != dv_names.end())
	{

		subset_idxs.push_back(bidx - dv_names.begin());
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
		std::uniform_int_distribution<int> uni(0, size - 1);
		int idx;
		for (int i = 0; i < 1000000000; i++)
		{
			if (subset_idxs.size() >= nreal_subset)
				break;
			idx = uni(subset_rand_gen);
			if (find(subset_idxs.begin(), subset_idxs.end(), idx) != subset_idxs.end())
				continue;
			subset_idxs.push_back(idx);
		}
		if (subset_idxs.size() != nreal_subset)
			throw_sqp_error("max iterations exceeded when trying to find random subset idxs");

	}
	else
	{
		throw_sqp_error("unknown 'subset_how'");
	}

	sort(subset_idxs.begin(), subset_idxs.end());
	return subset_idxs;

}

void CovMatAdapES::initialize(int n_params, int _num_reals)
{
	lambda = _num_reals;
	mu = lambda / 2;
	sigma = 1.0;

	m = Eigen::VectorXd::Zero(n_params);
	C = Eigen::MatrixXd::Identity(n_params, n_params);

	pc = Eigen::VectorXd::Zero(n_params);
	ps = Eigen::VectorXd::Zero(n_params);

	B = Eigen::MatrixXd::Identity(n_params, n_params);
	D = Eigen::VectorXd::Ones(n_params);

	weights.resize(mu);
	double sum_weights = 0.0;
	for (int i = 0; i < mu; i++) {
		weights[i] = log(mu + 0.5) - log(i + 1.0);
		sum_weights += weights[i];
	}

	for (int i = 0; i < mu; i++) {
		weights[i] /= sum_weights;
	}

	c_c = 0.5;
	c_m = 1.0;
	c_mu = pest_scenario_ptr->get_pestpp_options().get_sqp_cma_cmu();
	c_1 = pest_scenario_ptr->get_pestpp_options().get_sqp_cma_c1();
	/*c_1 = 2.0 / pow(static_cast<double>(n_params), 2);
	c_mu = min(mu/ pow(static_cast<double>(n_params), 2), 1 - c_1);*/
}

void CovMatAdapES::update(Parameters prev_m, Parameters curr_m) 
{
	m = curr_m.get_data_eigen_vec(curr_m.get_keys());
	vector<string> par_names = feas_dp_archive.get_var_names();
	ParameterEnsemble U = feas_dp_archive;
	vector<string> mu_best_real_names;
	int i = 0;
	for (auto r : U.get_real_names())
	{
		i++;
		if (i > mu)
			break;
		mu_best_real_names.push_back(r);

	}
	U.keep_rows(mu_best_real_names);
	Eigen::MatrixXd Uu_anoms = U.get_eigen_anomalies(vector<string>(), par_names) / sigma;
	//Uu_anoms.conservativeResize(Uu_anoms.rows() - 1, Uu_anoms.cols());
	
	Eigen::MatrixXd rank_mu_update = c_mu * Uu_anoms.transpose() * Uu_anoms / mu;
	
	pc = (1 - c_c) * pc + sqrt(c_c * (2 - c_c) * mu) * (curr_m.get_data_eigen_vec(par_names) - prev_m.get_data_eigen_vec(par_names)) / (c_m * sigma);
	Eigen::MatrixXd rank_one_update = c_1 * (pc * pc.transpose());

	C = (1.0 - c_1 - c_mu) * C + rank_mu_update + rank_one_update;
	cout << "CMA-ES approximated covariance: " << endl << C << endl;

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(C);
	B = eigensolver.eigenvectors();
	D = eigensolver.eigenvalues();

	for (int i = 0; i < D.size(); i++) {
		D(i) = max(D(i), 1e-12);
	}
}

void CovMatAdapES::update_archives(const ParameterEnsemble& pe, map<string, double> obj_map, map<string, double> viol_map, int iter, bool clear)
{
	vector<string> curr_names;
	
	if (clear) 
	{
		sorted_obj_map.clear();
		sorted_viol_map.clear();

		feas_dp_archive = pe;
		feas_dp_archive.keep_rows(vector<string>(), true);
		infeas_dp_archive = feas_dp_archive;
	}

	ParameterEnsemble curr_pe = pe;
	for (auto o : obj_map) 
	{
		curr_names.push_back(o.first);
		if (viol_map[o.first] == 0)
		{
			sorted_obj_map[to_string(iter) + "|" + o.first] = o.second;

			Parameters pars = pest_scenario_ptr->get_ctl_parameters();
			pars.update_without_clear(curr_pe.get_var_names(), curr_pe.get_real_vector(o.first));
			pe.get_par_transform().active_ctl2numeric_ip(pars);
			feas_dp_archive.append(to_string(iter) + "|" + o.first, pars);
		}
		else
		{
			sorted_viol_map[to_string(iter) + "|" + o.first] = viol_map[o.first];
			
			Parameters pars = pest_scenario_ptr->get_ctl_parameters();
			pars.update_without_clear(curr_pe.get_var_names(), curr_pe.get_real_vector(o.first));
			pe.get_par_transform().active_ctl2numeric_ip(pars);
			infeas_dp_archive.append(to_string(iter) + "|" + o.first, pars);
		}
	}

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
	feas_dp_archive.keep_rows(sorted_names_from_obj, true);
	feas_dp_archive.reorder(sorted_names_from_obj, curr_pe.get_var_names(), true);

	
	if (sorted_viol_map.size() > 0)
	{
		vector<pair<string, double>> sorted_viol_map_vec(sorted_viol_map.begin(), sorted_viol_map.end());
		sort(sorted_viol_map_vec.begin(), sorted_viol_map_vec.end(),
			[](const pair<string, double>& a, const pair<string, double>& b) {
				return a.second < b.second;
			});
		vector<string> sorted_names_from_viol;
		sorted_viol_map.clear();
		int i = 0;
		for (const auto& pair : sorted_viol_map_vec)
		{
			i++;
			if (i > lambda)
				break;
			sorted_names_from_viol.push_back(pair.first);
			sorted_viol_map[pair.first] = pair.second;
		}
		infeas_dp_archive.keep_rows(sorted_names_from_viol, true);
		infeas_dp_archive.reorder(sorted_names_from_viol, curr_pe.get_var_names(), true);
	}
}

ParameterEnsemble CovMatAdapES::generate_population(Parameters& _curr_m, ParameterEnsemble _dv) 
{	
	vector<string> rnames;
	vector<string> parnames = _dv.get_var_names();
	ParameterEnsemble new_reals = _dv;
	rnames = new_reals.get_real_names();
	auto it = find(rnames.begin(), rnames.end(), BASE_REAL_NAME);
	if (it != rnames.end()) {
		rnames.erase(it);
	}
	new_reals.keep_rows(rnames, true);

	const ParameterInfo& par_info = pest_scenario_ptr->get_ctl_parameter_info();

	for (int i = 0; i < lambda; i++) {
		Eigen::VectorXd x;
		int draws = 1;
		const int max_draws = 1000;

		while (true) 
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
			x = m + sigma * y;

			for (int j = 0; j < parnames.size(); ++j) {
				const ParameterRec* par_rec = par_info.get_parameter_rec_ptr(parnames[j]);
				double lbnd = par_rec->lbnd;
				double ubnd = par_rec->ubnd;

				if (par_rec->tranform_type == ParameterRec::TRAN_TYPE::LOG) {
					lbnd = log10(lbnd);
					ubnd = log10(ubnd);
				}

				if (x(j) < lbnd || x(j) > ubnd) {
					break;
				}
			}
			draws++;
		}
		new_reals.update_real_ip(rnames[i], x);
		
	}
	return new_reals;
}

Covariance CovMatAdapES::get_covariance_matrix(const vector<string>& par_names)
{
	vector<string> pnames = par_names;
	Covariance cov(pnames);

	// Get parameter bounds and scale the covariance matrix
	const ParameterInfo& par_info = pest_scenario_ptr->get_ctl_parameter_info();
	double sigma_range = pest_scenario_ptr->get_pestpp_options().get_par_sigma_range();

	vector<Eigen::Triplet<double>> triplets;
	for (int i = 0; i < C.rows(); ++i) {
		for (int j = 0; j < C.cols(); ++j) {
			if (C(i, j) != 0.0) {
				triplets.push_back(Eigen::Triplet<double>(i, j, C(i, j)));
			}
		}
	}

	Eigen::SparseMatrix<double> sparse_cov(par_names.size(), par_names.size());
	sparse_cov.setFromTriplets(triplets.begin(), triplets.end());
	cov = Covariance(par_names, sparse_cov);
	cout << "parcov: " << endl << cov << endl; //tmp
	return cov;
}
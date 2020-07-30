#include <random>
#include <iomanip>
#include <iterator>
#include "MOEA.h"
#include "Ensemble.h"
#include "RunManagerAbstract.h"
#include "ModelRunPP.h"
#include "RestartController.h"
#include "EnsembleMethodUtils.h"
#include "constraints.h"

using namespace std;

ParetoObjectives::ParetoObjectives(Pest& _pest_scenario, FileManager& _file_manager, PerformanceLog* _performance_log)
	: pest_scenario(_pest_scenario),file_manager(_file_manager),performance_log(_performance_log)

{

}
vector<string> ParetoObjectives::pareto_dominance_sort(ObservationEnsemble& op)
{
	//TODO: I think the nsga-II or spea-II sorts into levels, but this is just a dummy function for now
	stringstream ss;
	ss << "ParetoObjectives::pareto_dominance_sort() for " << op.shape().first << " population members";
	performance_log->log_event(ss.str());
	
	int nondom = 0;
	ss.str("");
	ss << nondom << " non-dominated members found in " << op.shape().first << " population";
	performance_log->log_event(ss.str());
	file_manager.rec_ofstream() << "...pareto dominance sort yielded " << nondom << "non-dominated members in population of " << op.shape().first << " members";
	return vector<string>();
}



MOEA::MOEA(Pest &_pest_scenario, FileManager &_file_manager, OutputFileWriter &_output_file_writer, 
	PerformanceLog *_performance_log, RunManagerAbstract* _run_mgr_ptr)
	: pest_scenario(_pest_scenario), file_manager(_file_manager),
	output_file_writer(_output_file_writer), performance_log(_performance_log),
	run_mgr_ptr(_run_mgr_ptr), constraints(_pest_scenario, &_file_manager, _output_file_writer, *_performance_log),
	objectives(_pest_scenario,_file_manager,_performance_log)
	
{
	rand_gen = std::mt19937(pest_scenario.get_pestpp_options().get_random_seed());
	dp.set_rand_gen(&rand_gen);
	dp.set_pest_scenario(&pest_scenario);
	op.set_rand_gen(&rand_gen);
	op.set_pest_scenario(&pest_scenario);
}


template<typename T, typename A>
void MOEA::message(int level, const string& _message, vector<T, A> _extras, bool echo)
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
	if (echo)// && ((verbose_level >= 2) || (level < 2)))
		cout << ss.str() << endl;
	file_manager.rec_ofstream() << ss.str() << endl;
	performance_log->log_event(_message);

}

void MOEA::message(int level, const string& _message)
{
	message(level, _message, vector<string>());
}

template<typename T>
void MOEA::message(int level, const string& _message, T extra)
{
	stringstream ss;
	ss << _message << " " << extra;
	string s = ss.str();
	message(level, s);
}

void MOEA::throw_moea_error(const string& message)
{
	performance_log->log_event("MOEA error: " + message);
	cout << endl << "   ************   " << endl << "    MOEAerror: " << message << endl << endl;
	file_manager.rec_ofstream() << endl << "   ************   " << endl << "    MOEA error: " << message << endl << endl;
	file_manager.close_file("rec");
	performance_log->~PerformanceLog();
	throw runtime_error("MOEA error: " + message);
}

void MOEA::sanity_checks()
{
	PestppOptions* ppo = pest_scenario.get_pestpp_options_ptr();
	vector<string> errors;
	vector<string> warnings;
	stringstream ss;

	if ((population_dv_file.size() == 0) && (population_obs_restart_file.size() > 0))
		errors.push_back("population_dv_file is empty but population_obs_restart_file is not - how can this work?");
	
	if ((ppo->get_mou_population_size() < error_min_members) && (population_dv_file.size() == 0))
	{
		ss.str("");
		ss << "population_size < " << error_min_members << ", this is redic, increaing to " << warn_min_members;
		warnings.push_back(ss.str());
		ppo->set_ies_num_reals(warn_min_members);
	}
	if ((ppo->get_mou_population_size() < warn_min_members) && (population_dv_file.size() == 0))
	{
		ss.str("");
		ss << "mou_population_size < " << warn_min_members << ", this is prob too few";
		warnings.push_back(ss.str());
	}
	if (ppo->get_ies_reg_factor() < 0.0)
		errors.push_back("ies_reg_factor < 0.0 - WRONG!");
	
	
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
		throw_moea_error(string("sanity_check() found some problems - please review rec file"));
	}
	//cout << endl << endl;
}


vector<int> MOEA::run_population(ParameterEnsemble& _pe, ObservationEnsemble& _oe, const vector<int>& real_idxs)
{
	message(1, "running population of size ", _pe.shape().first);
	stringstream ss;
	ss << "queuing " << _pe.shape().first << " runs";
	performance_log->log_event(ss.str());
	run_mgr_ptr->reinitialize();
	map<int, int> real_run_ids;
	try
	{
		real_run_ids = _pe.add_runs(run_mgr_ptr, real_idxs);
	}
	catch (const exception& e)
	{
		stringstream ss;
		ss << "run_ensemble() error queueing runs: " << e.what();
		throw_moea_error(ss.str());
	}
	catch (...)
	{
		throw_moea_error(string("run_ensemble() error queueing runs"));
	}
	performance_log->log_event("making runs");
	try
	{
		run_mgr_ptr->run();
	}
	catch (const exception& e)
	{
		stringstream ss;
		ss << "error running ensemble: " << e.what();
		throw_moea_error(ss.str());
	}
	catch (...)
	{
		throw_moea_error(string("error running ensemble"));
	}

	performance_log->log_event("processing runs");
	if (real_idxs.size() > 0)
	{
		_oe.keep_rows(real_idxs);
	}
	vector<int> failed_real_indices;
	try
	{
		failed_real_indices = _oe.update_from_runs(real_run_ids, run_mgr_ptr);
	}
	catch (const exception& e)
	{
		stringstream ss;
		ss << "error processing runs: " << e.what();
		throw_moea_error(ss.str());
	}
	catch (...)
	{
		throw_moea_error(string("error processing runs"));
	}
	//for testing
	//failed_real_indices.push_back(0);

	if (failed_real_indices.size() > 0)
	{
		stringstream ss;
		vector<string> par_real_names = _pe.get_real_names();
		vector<string> obs_real_names = _oe.get_real_names();
		ss << "the following par:obs realization runs failed: ";
		for (auto& i : failed_real_indices)
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

void MOEA::finalize()
{

}

void MOEA::initialize()
{
	message(0, "initializing MOEA process");
	
	pp_args = pest_scenario.get_pestpp_options().get_passed_args();

	act_obs_names = pest_scenario.get_ctl_ordered_nz_obs_names();
	act_par_names = pest_scenario.get_ctl_ordered_adj_par_names();

	//TODO: add arg to specify dv names rather than just use all adj names - can copy from pestpp-opt
	dv_names = act_par_names;

	stringstream ss;

	if (pest_scenario.get_control_info().noptmax == 0)
	{
		message(0, "'noptmax'=0, running control file parameter values and quitting");

		Parameters pars = pest_scenario.get_ctl_parameters();
		ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();

		ParameterEnsemble _pe(&pest_scenario, &rand_gen);
		_pe.reserve(vector<string>(), pest_scenario.get_ctl_ordered_par_names());
		_pe.set_trans_status(ParameterEnsemble::transStatus::CTL);
		_pe.append("BASE", pars);
		string par_csv = file_manager.get_base_filename() + ".par.csv";
		//message(1, "saving parameter values to ", par_csv);
		//_pe.to_csv(par_csv);
		ParameterEnsemble pe_base = _pe;
		pe_base.reorder(vector<string>(), act_par_names);
		ObservationEnsemble _oe(&pest_scenario, &rand_gen);
		_oe.reserve(vector<string>(), pest_scenario.get_ctl_ordered_obs_names());
		_oe.append("BASE", pest_scenario.get_ctl_observations());
		ObservationEnsemble oe_base = _oe;
		oe_base.reorder(vector<string>(), act_obs_names);
		//initialize the phi handler
		Covariance parcov;
		parcov.from_parameter_bounds(pest_scenario, file_manager.rec_ofstream());
		L2PhiHandler ph(&pest_scenario, &file_manager, &oe_base, &pe_base, &parcov);
		if (ph.get_lt_obs_names().size() > 0)
		{
			message(1, "less_than inequality defined for observations: ", ph.get_lt_obs_names().size());
		}
		if (ph.get_gt_obs_names().size())
		{
			message(1, "greater_than inequality defined for observations: ", ph.get_gt_obs_names().size());
		}
		message(1, "running control file parameter values");

		vector<int> failed_idxs = run_population(_pe, _oe);
		if (failed_idxs.size() != 0)
		{
			message(0, "control file parameter value run failed...bummer");
			throw_moea_error(string("control file parameter value run failed"));
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
	PestppOptions* ppo = pest_scenario.get_pestpp_options_ptr();

	
	iter = 0;
	member_count = 0;
	
	warn_min_members = 20;
	error_min_members = 4;
	
	message(1, "max run fail: ", ppo->get_max_run_fail());

	//TODO: how do we get info about the objective directions?  
	//one option would be for objective to use the 
	//constraint naming convention so that "less_than" would be minimize
	//and greater_than would be maximize. This would also make the risk-based objectives easy as...
	
	//TODO: deal with prior info objectives
	Constraints::ConstraintSense gt = Constraints::ConstraintSense::greater_than, lt = Constraints::ConstraintSense::less_than;
	pair<Constraints::ConstraintSense, string> sense;
	map<string, string> sense_map;
	vector<string> onames = pest_scenario.get_ctl_ordered_nz_obs_names();
	if (obj_names.size() == 0)
	{
		for (auto oname : onames)
		{
			sense = Constraints::get_sense_from_group_name(pest_scenario.get_ctl_observation_info().get_group(oname));
			if (sense.first == gt)
			{
				obj_names.push_back(oname);
				sense_map[oname] = "maximize";
			}
			else if (sense.first == lt)
			{
				obj_names.push_back(oname);
				sense_map[oname] = "minimize";
			}
		}

		message(1, "'mou_objectives' not passed, using all nonzero weighted obs that use the proper obs group naming convention");
	}
	else
	{
		vector<string> onames = pest_scenario.get_ctl_ordered_nz_obs_names();
		set<string> oset(onames.begin(), onames.end());
		onames.clear();
		vector<string> missing,keep,err_sense;
		for (auto obj_name : obj_names)
		{
			if (oset.find(obj_name) == oset.end())
				missing.push_back(obj_name);
			else
			{
				sense = Constraints::get_sense_from_group_name(pest_scenario.get_ctl_observation_info().get_group(obj_name));
				if ((sense.first != gt) && (sense.first != lt))
					err_sense.push_back(obj_name);
				else
				{
					if (sense.first == gt)
					{
						keep.push_back(obj_name);
						sense_map[obj_name] = "maximize";
					}
					else if (sense.first == lt)
					{
						keep.push_back(obj_name);
						sense_map[obj_name] = "minimize";
					}
					
				}
			}
		}
		if (err_sense.size() > 0)
		{
			ss.str("");
			ss << "the following non-zero weighted 'mou_objectives' do not have the correct obs group naming convention (needed to identify objective direction):";
			for (auto e : err_sense)
				ss << e << ";";
			throw_moea_error(ss.str());
		}
		if (keep.size() == 0)
		{
			throw_moea_error("none of the supplied 'mou_objectives' were found in the zero-weighted observations");
		}
		
		else if (missing.size() > 0)
		{
			ss.str("");
			ss << "WARNING: the following mou_objectives were not found in the zero-weighted observations: ";
			for (auto m : missing)
				ss << m << ",";
			message(1, ss.str());

		}
		obj_names = keep;


	}
	ss.str("");
	ss << "...using the following observations as objectives: " << endl;
	for (auto s : sense_map)
	{
		ss << setw(30) << s.first << "   " << s.second << endl;
	}
	file_manager.rec_ofstream() << ss.str();
	cout << ss.str();

	if (obj_names.size() > 5)
		message(1, "WARNING: more than 5 objectives, this is pushing the limits!");

	sanity_checks();


	int num_members = pest_scenario.get_pestpp_options().get_mou_population_size();
	population_dv_file = ppo->get_mou_dv_population_file();
	population_obs_restart_file = ppo->get_mou_obs_population_restart_file();
	
	initialize_dv_population();
	
	initialize_obs_restart_population();
	
	try
	{
		dp.check_for_dups();
	}
	catch (const exception& e)
	{
		string message = e.what();
		throw_moea_error("error in dv population: " + message);
	}

	try
	{
		op.check_for_dups();
	}
	catch (const exception& e)
	{
		string message = e.what();
		throw_moea_error("error in obs population: " + message);
	}


	//we are restarting
	if (population_obs_restart_file.size() > 0)
	{
	
		//since mou reqs strict linking of realization names, let's see if we can find an intersection set 
		vector<string> temp = dp.get_real_names();
		set<string> dvnames(temp.begin(), temp.end());
		temp = op.get_real_names();
		set<string> obsnames(temp.begin(), temp.end());
		set<string> common;
		set_intersection(dvnames.begin(), dvnames.end(), obsnames.begin(), obsnames.end(),std::inserter(common,common.end()));
		
		// all members are common to both dp and op
		if (common.size() == dp.shape().first)
		{
			op.reorder(dp.get_real_names(), vector<string>());
		}
		
		//otherwise some members are not common
		else
		{
			ss.str("");
			ss << "WARNING: only " << common.size() << " members are common between the dv population and obs restart population.";
			message(1, ss.str());
			if (common.size() < error_min_members)
			{
				throw_moea_error("too few members to continue");
			}
			
			message(2,"aligning dv and obs populations");
			temp.clear();
			temp.resize(common.size());
			copy(common.begin(), common.end(), temp.begin());
			sort(temp.begin(), temp.end());
			dp.reorder(temp, vector<string>(), true);
			op.reorder(temp, vector<string>(), true);
			message(2, "dv population size: ", dp.shape().first);
			message(2, "obs population size", op.shape().first);
			message(2, "checking for denormal values in dv population");
			dp.check_for_normal("initial transformed dv population");
			message(2, "checking for denormal values in obs restart population");
			op.check_for_normal("restart obs population");
		}
		
	}
	else
	{
		message(2, "checking for denormal values in dv population");
		dp.check_for_normal("initial transformed dv population");
		//save the initial population once here
		ss.str("");
		ss << file_manager.get_base_filename() << ".0." << dv_pop_file_tag << ".csv";
		dp.to_csv(ss.str());
		message(1, "saved initial dv population to ", ss.str());
		performance_log->log_event("running initial ensemble");
		message(1, "running initial ensemble of size", dp.shape().first);
		vector<int> failed = run_population(dp, op);
		if (dp.shape().first == 0)
			throw_moea_error(string("all members failed during initial population evaluation"));
		dp.transform_ip(ParameterEnsemble::transStatus::NUM);
	}
	ss.str("");
	ss << file_manager.get_base_filename() << ".0." << obs_pop_file_tag << ".csv";
	op.to_csv(ss.str());
	message(1, "saved observation population to ", ss.str());

	//save the initial dv population again in case runs failed or members were dropped as part of restart
	ss.str("");
	ss << file_manager.get_base_filename() << ".0." << dv_pop_file_tag << ".csv";
	dp.to_csv(ss.str());
	message(1, "saved initial dv population to ", ss.str());

	//TODO: think about a bad phi (or phis) for MOEA
	
	if (op.shape().first <= error_min_members)
	{
		message(0, "too few population members:", op.shape().first);
		message(1, "need at least ", error_min_members);
		throw_moea_error(string("too few active population members, cannot continue"));
	}
	if (op.shape().first < warn_min_members)
	{
		ss.str("");
		ss << "WARNING: less than " << warn_min_members << " active population members...might not be enough";
		string s = ss.str();
		message(0, s);
	}

	
	message(0, "initialization complete");
}


void MOEA::iterate_to_solution()
{


}


bool MOEA::initialize_dv_population()
{
	stringstream ss;
	int num_members = pest_scenario.get_pestpp_options().get_mou_population_size();
	string dv_filename = pest_scenario.get_pestpp_options().get_mou_dv_population_file();
	bool drawn = false;
	if (dv_filename.size() == 0)
	{
		ofstream& frec = file_manager.rec_ofstream();
		message(1, "drawing initial dv population of size: ", num_members);
		Parameters draw_par = pest_scenario.get_ctl_parameters();
		
		dp.draw_uniform(num_members, dv_names, performance_log, 1, file_manager.rec_ofstream());
		drawn = true;
	}
	else
	{
		string par_ext = pest_utils::lower_cp(dv_filename).substr(dv_filename.size() - 3, dv_filename.size());
		performance_log->log_event("processing dv population file " + dv_filename);
		if (par_ext.compare("csv") == 0)
		{
			message(1, "loading dv population from csv file", dv_filename);
			try
			{
				dp.from_csv(dv_filename);
			}
			catch (const exception& e)
			{
				ss << "error processing dv population file: " << e.what();
				throw_moea_error(ss.str());
			}
			catch (...)
			{
				throw_moea_error(string("error processing dv population file"));
			}
		}
		else if ((par_ext.compare("jcb") == 0) || (par_ext.compare("jco") == 0))
		{
			message(1, "loading dv population from binary file", dv_filename);
			try
			{
				dp.from_binary(dv_filename);
			}
			catch (const exception& e)
			{
				ss << "error processing binary file: " << e.what();
				throw_moea_error(ss.str());
			}
			catch (...)
			{
				throw_moea_error(string("error processing binary file"));
			}
		}
		else
		{
			ss << "unrecognized dv population file extension " << par_ext << ", looking for csv, jcb, or jco";
			throw_moea_error(ss.str());
		}

		dp.transform_ip(ParameterEnsemble::transStatus::NUM);
		ss.str("");
		ss << "dv population with " << dp.shape().first << " members read from '" << dv_filename << "'" << endl;
		message(1, ss.str());

		if (dp.shape().first == 0)
		{
			throw_moea_error("zero members found in dv population file");
		}

		if (pp_args.find("MOU_POPULATION_SIZE") != pp_args.end())
		{
			int num_members = pest_scenario.get_pestpp_options().get_mou_population_size();
			
			if (num_members < dp.shape().first)
			{
				message(1, "'mou_population_size' arg passed, truncated dv population to ", num_members);
				vector<string> keep_names, real_names = dp.get_real_names();
				for (int i = 0; i < num_members; i++)
				{
					keep_names.push_back(real_names[i]);
				}
				dp.keep_rows(keep_names);
			}
		}
		
	}
	return drawn;

}


void MOEA::initialize_obs_restart_population()
{
	string obs_filename = pest_scenario.get_pestpp_options().get_mou_obs_population_restart_file();
	if (obs_filename.size() == 0)
	{
		op.reserve(dp.get_real_names(), pest_scenario.get_ctl_ordered_obs_names());
		return;
	}
	stringstream ss;
	//throw_moea_error(string("restart obs population not implemented"));
	string par_ext = pest_utils::lower_cp(obs_filename).substr(obs_filename.size() - 3, obs_filename.size());
	performance_log->log_event("processing obs population file " + obs_filename);
	if (par_ext.compare("csv") == 0)
	{
		message(1, "loading obs population from csv file", obs_filename);
		try
		{
			op.from_csv(obs_filename);
		}
		catch (const exception& e)
		{
			ss << "error processing obs population file: " << e.what();
			throw_moea_error(ss.str());
		}
		catch (...)
		{
			throw_moea_error(string("error processing obs population file"));
		}
	}
	else if ((par_ext.compare("jcb") == 0) || (par_ext.compare("jco") == 0))
	{
		message(1, "loading obs population from binary file", obs_filename);
		try
		{
			op.from_binary(obs_filename);
		}
		catch (const exception& e)
		{
			ss << "error processing obs population binary file: " << e.what();
			throw_moea_error(ss.str());
		}
		catch (...)
		{
			throw_moea_error(string("error processing obs population binary file"));
		}
	}
	else
	{
		ss << "unrecognized obs population restart file extension " << par_ext << ", looking for csv, jcb, or jco";
		throw_moea_error(ss.str());
	}

	ss.str("");
	ss << "obs population with " << op.shape().first << " members read from '" << obs_filename << "'" << endl;
	message(1, ss.str());
	if (op.shape().first == 0)
	{
		throw_moea_error("zero members found in obs population restart file");
	}


}

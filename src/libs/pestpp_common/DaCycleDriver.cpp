#include "DaCycleDriver.h"

#include <algorithm>
#include <iostream>
#include <set>
#include <sstream>

#include "EnsembleMethodUtils.h"
#include "RunManagerSerial.h"
#include "utilities.h"

using namespace std;

DaCycleDriver::DaCycleDriver(Pest& _parent_scenario, FileManager& _file_manager,
	PerformanceLog* _performance_log, RunManagerAbstract* _run_mgr_ptr,
	RunManagerKind _rm_kind, bool _restart_flag)
	: parent_scenario(_parent_scenario), file_manager(_file_manager),
	performance_log(_performance_log), run_mgr_ptr(_run_mgr_ptr), rm_kind(_rm_kind),
	restart_flag(_restart_flag), start_cycle(0), cycle_index(-1), current_cycle(-1),
	curr_pe(&_parent_scenario), curr_oe(&_parent_scenario), curr_noise(&_parent_scenario),
	initialized(false), verbose_level(0), stop_requested(false), owns_cycle_run_mgr(false),
	use_existing(false), cycle_nnz_obs(0),
	cycle_curr_pe(&_parent_scenario), cycle_curr_oe(&_parent_scenario)
{
	pathname = ".";
}

void DaCycleDriver::set_pathname(const string& _pathname)
{
	pathname = _pathname;
}

int DaCycleDriver::get_n_cycles_remaining() const
{
	if (cycle_index < 0)
		return (int)assimilation_cycles.size();
	int remaining = (int)assimilation_cycles.size() - cycle_index;
	return (remaining < 0) ? 0 : remaining;
}

void DaCycleDriver::initialize()
{
	ofstream& frec = file_manager.rec_ofstream();

	// the cycle list comes from the assimilation/dci tables in the control file. we work it out
	// here instead of making the caller do it, because the noptmax schedule below is keyed on it
	// and if the two disagree nobody would ever notice.
	verbose_level = parent_scenario.get_pestpp_options().get_ies_verbose_level();

	// read the cycle tables before the cycle list - get_assim_dci_cycles() uses the cycles the
	// tables mention, so a table naming a cycle that isnt in the control file is how a cycle
	// gets into the sequence in the first place.
	par_cycle_info = process_da_par_cycle_table(parent_scenario, cycles_in_tables, frec);
	obs_cycle_info = process_da_obs_cycle_table(parent_scenario, cycles_in_tables, frec, obs_in_tbl);
	// use a separate set here on purpose. obs_in_tbl decides which obs a cycle zero-weights, so
	// if the weight table could add names to it we would end up zeroing obs that the obs table
	// never mentioned
	set<string> weights_in_tbl;
	weight_cycle_info = process_da_weight_cycle_table(parent_scenario, cycles_in_tables, frec,
		weights_in_tbl);

	assimilation_cycles = parent_scenario.get_assim_dci_cycles(frec, cycles_in_tables);
	sort(assimilation_cycles.begin(), assimilation_cycles.end());
	if (assimilation_cycles.size() == 0)
		throw runtime_error("no assimilation cycles found");

	// a throwaway DataAssimilator on the parent pest object, just to build the global ensembles
	// and work out the schedule. each cycle makes its own in begin_cycle() - this one never
	// assimilates anything.
	parent_ofw.reset(new OutputFileWriter(file_manager, parent_scenario, restart_flag));
	DataAssimilator bootstrap(parent_scenario, file_manager, *parent_ofw, performance_log,
		run_mgr_ptr);
	bootstrap.initialize_dynamic_states();
	bootstrap.sanity_checks();

	generate_global_ensembles(bootstrap, frec, curr_pe, curr_oe, curr_noise);
	// the base realization is a per-cycle thing - the global ensembles just carry realizations
	parent_scenario.get_pestpp_options_ptr()->set_ies_include_base(false);
	noptmax_schedule = bootstrap.initialize_noptmax_schedule(assimilation_cycles);

	string phi_file = file_manager.get_base_filename() + ".global.phi.actual.csv";
	f_phi.open(phi_file);
	if (f_phi.bad())
		throw runtime_error("error opening " + phi_file + " for writing");
	f_phi << "cycle,iteration,mean,standard_deviation,min,max";
	init_real_names = curr_oe.get_real_names();
	for (auto& real : init_real_names)
		f_phi << "," << pest_utils::lower_cp(real);
	f_phi << endl;

	cout << "initializing 'global' localizer" << endl;
	frec << "initializing 'global' localizer" << endl;
	global_loc.reset(new Localizer(&parent_scenario));
	bool forgive_missing =
		parent_scenario.get_pestpp_options().get_ies_localizer_forgive_missing();
	global_loc->initialize(performance_log, frec, forgive_missing);

	pcs.reset(new ParChangeSummarizer(&curr_pe, &file_manager, parent_ofw.get()));

	start_cycle = parent_scenario.get_pestpp_options().get_da_hotstart_cycle();
	int max_cycle = assimilation_cycles[assimilation_cycles.size() - 1];
	if (start_cycle > max_cycle)
	{
		stringstream ss;
		ss << "'da_hotstart_cycle' (" << start_cycle << ") greater than max cycle (" << max_cycle << ")";
		throw runtime_error(ss.str());
	}

	cycle_index = 0;
	current_cycle = -1;
	initialized = true;
}

void DaCycleDriver::finalize()
{
	if (f_phi.is_open())
		f_phi.flush();
	// flush but dont close. closing kills the stream, and then any later write - including one
	// from a destructor while an exception is unwinding - takes the whole process down
}

bool DaCycleDriver::begin_cycle()
{
	if (!initialized)
		throw runtime_error("DaCycleDriver::begin_cycle(): initialize() first");
	if (stop_requested)
		return false;

	ofstream& frec = file_manager.rec_ofstream();
	stringstream ss;

	// handle ++da_hotstart_cycle in here instead of making the caller do it, so somebody
	// stepping through cycles by hand cant accidentally assimilate one the exe would have
	// skipped
	while ((cycle_index < (int)assimilation_cycles.size()) &&
		(assimilation_cycles[cycle_index] < start_cycle))
	{
		cout << "fast-forwarding past cycle " << assimilation_cycles[cycle_index] << endl;
		cycle_index++;
	}
	if (cycle_index >= (int)assimilation_cycles.size())
		return false;

	current_cycle = assimilation_cycles[cycle_index];
	const int icycle = current_cycle;

	cout << endl << endl;
	cout << " =======================================" << endl;
	cout << " >>>> processing cycle " << icycle << endl;
	cout << " =======================================" << endl << endl;

	frec << endl;
	frec << " =======================================" << endl;
	frec << " >>>> processing cycle " << icycle << endl;
	frec << " =======================================" << endl;

	performance_log->log_event("instantiating child pest object");

	child_scenario.reset(new Pest(parent_scenario.get_child_pest(icycle)));
	Pest& childPest = *child_scenario;
	childPest.get_control_info_4_mod().noptmax = noptmax_schedule[icycle];
	child_ofw.reset(new OutputFileWriter(file_manager, childPest, restart_flag));
	OutputFileWriter& output_file_writer = *child_ofw;

	cout << "checking inputs...";
	childPest.check_inputs(frec, false, true, icycle);
	cout << "done" << endl;

	// the run manager for this cycle. panther gets re-initialized against the cycle's par/obs so
	// the workers and master agree. serial and external get rebuilt as a RunManagerSerial
	// because the model exec info comes from the child pest object, which is different each
	// cycle.
	//
	// external getting a serial manager looks like a mistake but it isnt. pestpp-da /e builds a
	// RunManagerExternal for the parent, then swaps in a serial one each cycle that actually
	// runs the model - without it the cycles dont run at all. i tried "just keep the caller's
	// manager" thinking it was a copy-paste bug and basic_xplat_g07_test went from total_runs
	// 140 to 0. no cycle ran a single model call and the run still finished, reporting the
	// initial values like they were results.
	//
	// other keeps the caller's manager - that is the c api case, where the session owns it.
	owns_cycle_run_mgr = false;
	if (rm_kind == RunManagerKind::PANTHER_MASTER)
	{
		performance_log->log_event("reinitializing panther master");
		run_mgr_ptr->initialize(childPest.get_ctl_parameters(), childPest.get_ctl_observations());
	}
	else if ((rm_kind == RunManagerKind::SERIAL) || (rm_kind == RunManagerKind::EXTERNAL))
	{
		performance_log->log_event("starting basic model IO error checking");
		cout << "checking model IO files...";
		childPest.check_io(frec);
		performance_log->log_event("finished basic model IO error checking");
		cout << "done" << endl;
		const ModelExecInfo& exi = childPest.get_model_exec_info();
		run_mgr_ptr = new RunManagerSerial(exi.comline_vec,
			exi.tplfile_vec, exi.inpfile_vec, exi.insfile_vec, exi.outfile_vec,
			file_manager.build_filename("rns"), pathname,
			childPest.get_pestpp_options().get_max_run_fail(),
			childPest.get_pestpp_options().get_fill_tpl_zeros(),
			childPest.get_pestpp_options().get_additional_ins_delimiters(),
			parent_scenario.get_pestpp_options().get_num_tpl_ins_threads(),
			parent_scenario.get_pestpp_options().get_tpl_force_decimal());
		owns_cycle_run_mgr = true;
	}

	ParamTransformSeq& base_trans_seq = childPest.get_base_par_tran_seq_4_mod();

	//check for entries in the par cycle table
	if (par_cycle_info.find(icycle) != par_cycle_info.end())
	{
		performance_log->log_event("updating pars using da par cycle table info");
		map<string, double> cycle_map = par_cycle_info[icycle];
		for (const auto& item : cycle_map)
		{
			base_trans_seq.get_fixed_ptr_4_mod()->insert(item.first, item.second);
			childPest.get_ctl_parameters_4_mod().update_rec(item.first, item.second);
		}
	}
	//check for entries in the obs cycle table
	if (obs_cycle_info.find(icycle) != obs_cycle_info.end())
	{
		performance_log->log_event("updating obs using da obs cycle table info");
		map<string, double> cycle_map = obs_cycle_info[icycle];
		map<string, double> weight_cycle_map;
		if (weight_cycle_info.find(icycle) != weight_cycle_info.end())
			weight_cycle_map = weight_cycle_info[icycle];
		ObservationInfo oi = parent_scenario.get_ctl_observation_info();
		childPest.set_observation_info(oi);
		for (auto tbl_obs_name : obs_in_tbl)
		{
			//if this obs is not in this cycle, give it a zero weight
			if (cycle_map.find(tbl_obs_name) == cycle_map.end())
			{
				oi.set_weight(tbl_obs_name, 0.0);
			}
			else
			{
				childPest.get_ctl_observations_4_mod().update_rec(tbl_obs_name,
					cycle_map.at(tbl_obs_name));
				//check if this obs is in this cycle's weight info
				if (weight_cycle_map.find(tbl_obs_name) != weight_cycle_map.end())
					oi.set_weight(tbl_obs_name, weight_cycle_map.at(tbl_obs_name));
			}
		}
		childPest.set_observation_info(oi);
	}

	if (verbose_level > 1) // todo: add da verbose
	{
		frec << "...verbose level > 1, cycle " << icycle << " Pest Interface Summary:" << endl;
		output_file_writer.scenario_io_report(frec);
		output_file_writer.scenario_pargroup_report(frec);
		output_file_writer.scenario_par_report(frec);
		output_file_writer.scenario_obs_report(frec);
	}

	Parameters par1 = childPest.get_ctl_parameters();
	base_trans_seq.ctl2numeric_ip(par1);
	base_trans_seq.numeric2model_ip(par1);
	ParameterInfo pi = parent_scenario.get_ctl_parameter_info();

	int nadj_par = 0;
	for (const auto& par : par1)
	{
		if (cycle_in_range(icycle, pi.get_parameter_rec_ptr(par.first)->dci))
		{
			if ((pi.get_parameter_rec_ptr(par.first)->tranform_type != ParameterRec::TRAN_TYPE::FIXED) &&
				(pi.get_parameter_rec_ptr(par.first)->tranform_type != ParameterRec::TRAN_TYPE::TIED))
				nadj_par++;
		}
	}

	cout << "...number of parameters in cycle " << icycle << ": " << childPest.get_ctl_parameters().size() << endl;
	frec << "...number of parameters in cycle " << icycle << ": " << childPest.get_ctl_parameters().size() << endl;
	cout << "...number of adjustable parameters in cycle " << icycle << ": " << nadj_par << endl;
	frec << "...number of adjustable parameters in cycle " << icycle << ": " << nadj_par << endl;

	ObservationInfo oi = childPest.get_ctl_observation_info();
	cycle_nnz_obs = 0;
	for (const auto& o : childPest.get_ctl_observations())
	{
		if (oi.get_observation_rec_ptr(o.first)->weight != 0.0)
			cycle_nnz_obs++;
	}

	cout << "...number of observations in cycle " << icycle << ": " << childPest.get_ctl_observations().size() << endl;
	frec << "...number of observations in cycle " << icycle << ": " << childPest.get_ctl_observations().size() << endl;
	cout << "...number of non-zero weighted observations in cycle " << icycle << ": " << cycle_nnz_obs << endl;
	frec << "...number of non-zero weighted observations in cycle " << icycle << ": " << cycle_nnz_obs << endl;
	cout << "...number of template files in cycle " << icycle << ": " << childPest.get_model_exec_info().tplfile_vec.size() << endl;
	frec << "...number of template files in cycle " << icycle << ": " << childPest.get_model_exec_info().tplfile_vec.size() << endl;
	cout << "...number of instruction files in cycle " << icycle << ": " << childPest.get_model_exec_info().insfile_vec.size() << endl;
	frec << "...number of instruction files in cycle " << icycle << ": " << childPest.get_model_exec_info().insfile_vec.size() << endl;

	Parameters cur_ctl_parameters = childPest.get_ctl_parameters();
	//Allocates Space for Run Manager.  This initializes the model parameter names and observations names.
	//make sure we use the vector-based initializer so that the pars and obs are in order on the
	//workers - PantherAgent uses this same strategy (child pest with cycle number, then sorted par and
	//obs names)
	vector<string> par_names = base_trans_seq.ctl2model_cp(cur_ctl_parameters).get_keys();
	sort(par_names.begin(), par_names.end());
	vector<string> obs_names = childPest.get_ctl_observations().get_keys();
	sort(obs_names.begin(), obs_names.end());
	run_mgr_ptr->set_save_all_runs(parent_scenario.get_pestpp_options().get_save_all_runs());
	run_mgr_ptr->initialize(par_names, obs_names);

	performance_log->log_event("instantiating DA instance");
	da.reset(new DataAssimilator(childPest, file_manager, output_file_writer, performance_log,
		run_mgr_ptr));
	da->initialize_dynamic_states();
	std::mt19937 rand_gen = da->get_rand_gen();

	vector<string> act_par_names = childPest.get_ctl_ordered_adj_par_names();
	performance_log->log_event("instantiating cycle parameter ensemble instance");
	cycle_curr_pe = ParameterEnsemble(&childPest, &rand_gen,
		curr_pe.get_eigen(vector<string>(), act_par_names), curr_pe.get_real_names(), act_par_names);
	cycle_curr_pe.set_trans_status(curr_pe.get_trans_status());
	cycle_curr_pe.set_fixed_info(curr_pe.get_fixed_info());
	if (par_cycle_info.find(icycle) != par_cycle_info.end())
	{
		map<string, double> cycle_info = par_cycle_info.at(icycle);
		cycle_curr_pe.get_fixed_info().update_par_values(cycle_info);
	}
	da->set_pe(cycle_curr_pe);

	obs_names = childPest.get_ctl_ordered_obs_names();
	cycle_curr_oe = ObservationEnsemble(&childPest, &rand_gen,
		curr_oe.get_eigen(vector<string>(), obs_names), curr_oe.get_real_names(), obs_names);
	da->set_oe(cycle_curr_oe);

	if (cycle_nnz_obs > 0)
	{
		obs_names = childPest.get_ctl_ordered_nz_obs_names();
		ObservationEnsemble cycle_curr_noise(&childPest, &rand_gen,
			curr_noise.get_eigen(cycle_curr_oe.get_real_names(), obs_names),
			cycle_curr_oe.get_real_names(), obs_names);
		//correct for obs cycle table
		if ((parent_scenario.get_pestpp_options().get_da_weight_cycle_table().size() > 0) &&
			(!parent_scenario.get_pestpp_options().get_ies_no_noise()))
		{
			cout << "...'da_weight_cycle_table' passed, redrawing noise realizations" << endl;
			frec << "...'da_weight_cycle_table' passed, redrawing noise realizations" << endl;
			da->initialize_obscov();
			cycle_curr_noise.draw(cycle_curr_oe.shape().first, *da->get_obscov_ptr(),
				performance_log, verbose_level, frec);
			vector<string> names = cycle_curr_oe.get_real_names();
			cycle_curr_noise.set_real_names(names);
		}
		else if (parent_scenario.get_pestpp_options().get_da_obs_cycle_table().size() > 0)
		{
			cout << "...'da_obs_cycle_table' passed, translating noise realizations" << endl;
			frec << "...'da_obs_cycle_table' passed, translating noise realizations" << endl;
			cycle_curr_noise.update_var_map();
			Observations org_obs = parent_scenario.get_ctl_observations();
			map<string, int> vmap = cycle_curr_noise.get_var_map();
			for (const auto& o : childPest.get_ctl_observations())
			{
				if (oi.get_observation_rec_ptr(o.first)->weight != 0.0)
				{
					cycle_curr_noise.get_eigen_ptr_4_mod()->col(vmap[o.first]).array() +=
						o.second - org_obs[o.first];
				}
			}
		}
		da->set_noise_oe(cycle_curr_noise);
		if (verbose_level >= 3)
			cycle_curr_noise.to_csv("cycle_curr_noise.csv");
	}
	else
	{
		da->set_noise_oe(cycle_curr_oe);
	}
	da->set_localizer(*global_loc);

	//check if we can use the previous outputs.
	//if there are no dynamic states, then we can possibly reuse sim outputs
	use_existing = false;
	if (da->get_par_dyn_state_names().size() == 0)
	{
		// if npar is the same...
		if (prev_cycle_ctl_pars.size() == cur_ctl_parameters.size())
		{
			//if nobs is the same
			if (prev_cycle_obs.size() == childPest.get_ctl_observations().size())
			{
				//check that all pars have the same values
				Parameters::iterator pend = cur_ctl_parameters.end();
				bool pars_same = true;
				for (auto& p : prev_cycle_ctl_pars)
				{
					if (cur_ctl_parameters.find(p.first) == pend)
					{
						pars_same = false;
						break;
					}
					if (cur_ctl_parameters.get_rec(p.first) != p.second)
					{
						pars_same = false;
						break;
					}
				}
				if (pars_same)
				{
					//check that all obs have the same values
					Observations cur_cycle_obs = childPest.get_ctl_observations();
					bool obs_same = true;
					Observations::iterator oend = cur_cycle_obs.end();
					for (auto& o : prev_cycle_obs)
					{
						if (cur_cycle_obs.find(o.first) == oend)
						{
							obs_same = false;
							break;
						}
						if (cur_cycle_obs.get_rec(o.first) != o.second)
						{
							obs_same = false;
							break;
						}
					}
					//if we get to here, it should be safe to reuse the sim outputs from the
					//previous cycle...
					if (obs_same)
						use_existing = true;
				}
			}
		}
	}

	if (use_existing)
	{
		ss.str("");
		ss << "...parameters and observations are consistent with previous cycle, reusing existing simulated outputs" << endl;
		cout << ss.str();
		frec << ss.str();
	}
	output_file_writer.scenario_io_report(frec);
	if (verbose_level > 1)
	{
		output_file_writer.scenario_pargroup_report(frec);
		output_file_writer.scenario_par_report(frec);
		output_file_writer.scenario_obs_report(frec);
	}

	return true;
}

void DaCycleDriver::drive_cycle()
{
	if (!da)
		throw runtime_error("DaCycleDriver::drive_cycle(): begin_cycle() first");

	ofstream& frec = file_manager.rec_ofstream();
	stringstream ss;
	const int icycle = current_cycle;

	da->initialize(icycle, true, use_existing);

	write_global_phi_info(icycle, f_phi, *da, init_real_names);

	da->transfer_dynamic_state_from_oe_to_final_pe(*da->get_pe_ptr(), *da->get_oe_ptr());

	if (child_scenario->get_ctl_ordered_nz_obs_names().size() > 0)
	{
		if (parent_scenario.get_control_info().noptmax > 0)
		{
			if ((da->get_phi_handler().get_mean(L2PhiHandler::phiType::ACTUAL) == 0) ||
				(da->get_phi_handler().get_mean(L2PhiHandler::phiType::MEAS) == 0))
			{
				ss.str("");
				ss << "...Note:current mean actual and/or measurement phi is too low for solution, continuing..." << endl;
				frec << ss.str();
				cout << ss.str();
			}
			else
			{
				da->da_update(icycle);
				ss.str("");
				ss << file_manager.get_base_filename() << ".global." << icycle << "."
					<< da->get_iter() << ".pcs.csv";
				pcs->summarize(*da->get_pe_ptr(), ss.str());
			}
		}
		write_global_phi_info(icycle, f_phi, *da, init_real_names);
	}
	else
	{
		cout << "...Note: no non-zero-weighted observations in cycle " << icycle << ", continuing..." << endl;
		frec << "...Note: no non-zero-weighted observations in cycle " << icycle << ", continuing..." << endl;
		performance_log->log_event("no non-zero-weighted observations in current cycle");
	}
}

void DaCycleDriver::end_cycle()
{
	if (!da)
		throw runtime_error("DaCycleDriver::end_cycle(): begin_cycle() first");

	ofstream& frec = file_manager.rec_ofstream();
	stringstream ss;
	const int icycle = current_cycle;

	//replace all the pars used in this cycle in the parent parameter ensemble
	performance_log->log_event("updating global ensemble with cycle ensemble columns");
	cycle_curr_pe = da->get_pe();
	drop_missing_realizations(cycle_curr_pe);
	cycle_curr_pe.transform_ip(curr_pe.get_trans_status());
	curr_pe.replace_col_vals(cycle_curr_pe.get_var_names(), *cycle_curr_pe.get_eigen_ptr());

	cycle_curr_oe = da->get_oe();
	drop_missing_realizations(cycle_curr_oe);
	curr_oe.replace_col_vals(cycle_curr_oe.get_var_names(), *cycle_curr_oe.get_eigen_ptr());

	ss.str("");
	ss << ".global." << icycle << ".oe.csv";
	curr_oe.to_csv(file_manager.get_base_filename() + ss.str());
	ss.str("");
	ss << ".global." << icycle << ".pe.csv";
	curr_pe.to_csv(file_manager.get_base_filename() + ss.str());

	file_manager.close_all_files_containing(".phi.");

	//transfer the best (current) simulated final states to the initial states pars in the pe for
	//the next cycle
	if (parent_scenario.get_pestpp_options().get_da_use_simulated_states())
		da->transfer_dynamic_state_from_oe_to_initial_pe(curr_pe, curr_oe);
	else
		da->transfer_par_dynamic_state_final_to_initial_ip(curr_pe);

	// clean the cycle up before deciding whether to stop, so stopping leaves things in the same
	// state as a normal cycle ending - get_da() null, run manager freed. main() could just break
	// out of its loop and let scope handle it, but we cant do that here.
	bool stop_now = false;
	if (icycle >= parent_scenario.get_pestpp_options().get_da_stop_cycle())
	{
		cout << "'da_stop_cycle' criteria met" << endl;
		frec << "'da_stop_cycle' criteria met" << endl;
		stop_now = true;
	}
	else
	{
		prev_cycle_ctl_pars = Parameters(child_scenario->get_ctl_parameters());
		prev_cycle_obs = Observations(child_scenario->get_ctl_observations());

		int q = pest_utils::quit_file_found();
		if ((q == 1) || (q == 2))
		{
			cout << "'pest.stp' found, quitting" << endl;
			frec << "'pest.stp' found, quitting" << endl;
			stop_now = true;
		}
		else if (q == 4)
		{
			cout << "...pest.stp found with '4'.  run mgr has returned control, removing file." << endl;
			frec << "...pest.stp found with '4'.  run mgr has returned control, removing file." << endl;
			if (!pest_utils::try_remove_quit_file())
			{
				cout << "...error removing pest.stp file, bad times ahead..." << endl;
				frec << "...error removing pest.stp file, bad times ahead..." << endl;
			}
		}
	}

	da.reset();
	if (owns_cycle_run_mgr)
	{
		delete run_mgr_ptr;
		run_mgr_ptr = nullptr;
		owns_cycle_run_mgr = false;
	}

	stop_requested = stop_now;
	cycle_index++;
}

/// any realizations a cycle lost have to come out of the global ensemble too, otherwise the
/// next cycle builds its view from rows the last cycle doesnt have anymore. that shows up as
/// garbage values instead of an error, which is worse.
void DaCycleDriver::drop_missing_realizations(ParameterEnsemble& cycle_pe)
{
	if (curr_pe.shape().first <= cycle_pe.shape().first)
		return;
	ofstream& frec = file_manager.rec_ofstream();
	vector<string> missing;
	vector<string> t = cycle_pe.get_real_names();
	set<string> ccp_reals(t.begin(), t.end());
	for (auto r : curr_pe.get_real_names())
		if (ccp_reals.find(r) == ccp_reals.end())
			missing.push_back(r);
	frec << "...dropping the following " << missing.size() << " realizations from global parameter ensemble:" << endl;
	cout << "...dropping " << missing.size() << " realizations from global parameter ensemble, see rec file for listing" << endl;
	int i = 0;
	for (auto m : missing)
	{
		frec << m << ",";
		if (i > 10)
		{
			frec << endl;
			i = 0;
		}
	}
	curr_pe.drop_rows(missing);
	curr_pe.reorder(t, vector<string>());
}

void DaCycleDriver::drop_missing_realizations(ObservationEnsemble& cycle_oe)
{
	if (curr_oe.shape().first <= cycle_oe.shape().first)
		return;
	ofstream& frec = file_manager.rec_ofstream();
	vector<string> missing;
	vector<string> t = cycle_oe.get_real_names();
	set<string> ccp_reals(t.begin(), t.end());
	for (auto r : curr_oe.get_real_names())
		if (ccp_reals.find(r) == ccp_reals.end())
			missing.push_back(r);
	frec << "...dropping the following " << missing.size() << " realizations from global observation ensemble:" << endl;
	cout << "...dropping " << missing.size() << " realizations from global observation ensemble, see rec file for listing" << endl;
	int i = 0;
	for (auto m : missing)
	{
		frec << m << ",";
		if (i > 10)
		{
			frec << endl;
			i = 0;
		}
	}
	curr_oe.drop_rows(missing);
	curr_oe.reorder(t, vector<string>());
}

void DaCycleDriver::run_all_cycles()
{
	if (!initialized)
		initialize();
	while (begin_cycle())
	{
		drive_cycle();
		end_cycle();
	}
	finalize();
}

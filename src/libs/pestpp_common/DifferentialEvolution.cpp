#include <random>
#include <iomanip>
#include "DifferentialEvolution.h"
#include "RunManagerAbstract.h"
#include "ModelRunPP.h"
#include "RestartController.h"

mt19937_64 DifferentialEvolution::rand_engine = mt19937_64(1);

const string DifferentialEvolution::solver_type_name = "differential_evolution";

DifferentialEvolution::DifferentialEvolution(Pest &_pest_scenario,
	FileManager &_file_manager, ObjectiveFunc *_obj_func_ptr,
	const ParamTransformSeq &_par_transform, OutputFileWriter &_output_file_writer,
	PerformanceLog *_performance_log, unsigned int seed)
	: file_manager(_file_manager), obj_func_ptr(_obj_func_ptr), par_transform(_par_transform),
	output_file_writer(_output_file_writer), performance_log(_performance_log),
	gen_1(_file_manager.build_filename("de1"))
{
	// initialize random number generator
	rand_engine.seed(seed);

	// DE only works for one to one transformations
	if (!par_transform.is_one_to_one())
	{
		throw PestError("Error: Differential Evolution only supports one to one transformations.  Please insure the SVDA transformation is turned off.");
	}

	par_list = _pest_scenario.get_ctl_ordered_adj_par_names();
	Parameters inti_pars = _pest_scenario.get_ctl_parameters();
	for (const auto &i : inti_pars)
	{
		const string &p_name = i.first;
		const ParameterRec *p_info = _pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(p_name);
		max_numeric_pars[p_name] = p_info->ubnd;
		min_numeric_pars[p_name] = p_info->lbnd;
	}
	par_transform.ctl2numeric_ip(max_numeric_pars);
	par_transform.ctl2numeric_ip(min_numeric_pars);
}

void DifferentialEvolution::solve(RunManagerAbstract &run_manager,
	RestartController &restart_controller,
	int max_gen, double f, double cr, bool dither_f, ModelRun &cur_run)
{
	ostream &os = file_manager.rec_ofstream();
	ostream &fout_restart = file_manager.get_ofstream("rst");

	for (int iter = 0; iter < max_gen && best_phi > std::numeric_limits<double>::min(); ++iter)
	{
		RestartController::write_start_iteration(fout_restart, solver_type_name, iter+1, iter+1);

		Parameters tmp_pars;
		Observations tmp_obs;

		// write header for iteration
		cout << endl;
		output_file_writer.iteration_report(cout, iter+1, run_manager.get_total_runs(), "differential evolution");
		os << endl;
		output_file_writer.iteration_report(os, iter + 1, run_manager.get_total_runs(), "differential evolution");
		// write initial phi report for this iteration
		bool run_target_ok = gen_1.get_run(best_run_idx, tmp_pars, tmp_obs);
		par_transform.model2ctl_ip(tmp_pars);
		PhiData phi_data  = obj_func_ptr->phi_report(tmp_obs, tmp_pars, DynamicRegularization::get_unit_reg_instance());
		output_file_writer.phi_report(cout, iter, run_manager.get_nruns(), phi_data, DynamicRegularization::get_unit_reg_instance().get_weight(), false);
		output_file_writer.phi_report(os, iter, run_manager.get_nruns(), phi_data, DynamicRegularization::get_unit_reg_instance().get_weight(), false);

		run_manager.reinitialize();
		mutation(run_manager, f, dither_f, cr);
		RestartController::write_upgrade_runs_built(fout_restart);
		// make trial vector model runs
		cout << endl;
		cout << "  performing trial vector model runs... ";
		cout.flush();
		run_manager.run();
		os << endl;
		best_run_idx = recombination(run_manager);
		os << endl;

		run_target_ok = gen_1.get_run(best_run_idx, tmp_pars, tmp_obs);
		par_transform.model2ctl_ip(tmp_pars);
		// write parameter file for this iteration
		output_file_writer.write_par(file_manager.open_ofile_ext("par"), tmp_pars, *(par_transform.get_offset_ptr()),
			*(par_transform.get_scale_ptr()));
		file_manager.close_file("par");
		// write final phi report for this iteration
		phi_data = obj_func_ptr->phi_report(tmp_obs, tmp_pars, DynamicRegularization::get_unit_reg_instance());
		output_file_writer.phi_report(cout, iter, run_manager.get_nruns(), phi_data, DynamicRegularization::get_unit_reg_instance().get_weight(), true);
		cout << endl;
		output_file_writer.phi_report(os, iter, run_manager.get_nruns(), phi_data, DynamicRegularization::get_unit_reg_instance().get_weight(), true);
		os << endl;
	}
}


void DifferentialEvolution::initialize_population(RunManagerAbstract &run_manager, int d)
{
	ostream &fout_restart = file_manager.get_ofstream("rst");
	int iter = 0;

	RestartController::write_start_iteration(fout_restart, solver_type_name, iter, iter);
	Parameters numeric_pars;
	for (int i = 0; i < d; ++i)
	{
		numeric_pars.clear();
		initialize_vector(numeric_pars);
		par_transform.numeric2model_ip(numeric_pars);
		run_manager.add_run(numeric_pars);
	}
	RestartController::write_upgrade_runs_built(fout_restart);
	// make initial population vector model runs
	cout << endl;
	cout << "  performing initial population model runs... ";
	cout.flush();
	run_manager.run();
	gen_1.copy(run_manager.get_runstorage_ref());

	// get the best_run to track phi
	int r_status;
	int n_par = par_list.size();
	best_phi = std::numeric_limits<double>::max();

	Parameters tmp_pars;
	Observations tmp_obs;
	ModelRun tmp_run(obj_func_ptr);
	for (int i_run = 0; i_run < d; ++i_run)
	{
		bool r_status = gen_1.get_run(i_run, tmp_pars, tmp_obs);
		if (r_status)
		{
			par_transform.model2ctl_ip(tmp_pars);
			tmp_run.update_ctl(tmp_pars, tmp_obs);
			double tmp_phi = tmp_run.get_phi(DynamicRegularization::get_unit_reg_instance());
			if (tmp_phi < best_phi)
			{
				best_phi = tmp_phi;
				best_run_idx = i_run;
			}
		}
	}

}

void DifferentialEvolution::initialize_vector(Parameters &numeric_pars)
{
	std::uniform_real_distribution<double> distribution(0.0, 1.0);

	for (const auto &i : par_list)
	{
		double p_val;
		const string &p_name = i;
		double p_min = min_numeric_pars[p_name];
		double p_max = max_numeric_pars[p_name];
		double rnum = rand_engine();

		// use uniform distribution to initialize parameters
		p_val = p_min + distribution(rand_engine) * (p_max - p_min);
		numeric_pars.insert(p_name, p_val);
	}
}

void DifferentialEvolution::mutation(RunManagerAbstract &run_manager, double f, bool dither_f, double cr)
{
	int d = gen_1.get_nruns();
	int r_status;
	int n_par = par_list.size();
	std::uniform_int_distribution<int> uni_par(0, n_par-1);
	std::uniform_real_distribution<double> cr_prob(0.0, 1.0);
	std::uniform_real_distribution<double> dither_f_prob(.5, 1.0);
	vector<int> successful_run_ids;
	// generate a vector of successful runs
	for (int i_run = 0; i_run < d; ++i_run)
	{
		r_status = gen_1.get_run_status(i_run);
		if (r_status > 0)
		{
			successful_run_ids.push_back(i_run);
		}
	}
	int d_ok = successful_run_ids.size();
	std::uniform_int_distribution<int> uni_run_ok(0, d_ok-1);

	Parameters xa;
	Parameters xb;
	Parameters xc;
	Parameters x_trial;
	for (int i_run = 0; i_run < d; ++i_run)
	{
		int xa_id = successful_run_ids[uni_run_ok(rand_engine)];
		int xb_id = successful_run_ids[uni_run_ok(rand_engine)];
		while (xa_id == xb_id)
		{
			xb_id = successful_run_ids[uni_run_ok(rand_engine)];
		}
		int xc_id = successful_run_ids[uni_run_ok(rand_engine)];
		//initialize trail vector with the target vector
		gen_1.get_parameters(i_run, x_trial);
		gen_1.get_parameters(xa_id, xa);
		gen_1.get_parameters(xb_id, xb);
		gen_1.get_parameters(xc_id, xb);
		par_transform.model2numeric_ip(xa);
		par_transform.model2numeric_ip(xb);
		par_transform.model2numeric_ip(xc);
		int par_id_chg = uni_par(rand_engine);
		for (int idx=0; idx<n_par; ++idx)
		{
			const string &ipar = par_list[idx];
			double a = xa[ipar];
			double b = xb[ipar];
			double c = xc[ipar];
			double delta = a - b;
			double tmp_f = f;
			if (dither_f)
			{
				tmp_f = dither_f_prob(rand_engine);
			}
			double c_p = c + tmp_f * (delta);
			double p_min = min_numeric_pars[ipar];
			double p_max = max_numeric_pars[ipar];
			// reflect purturbation if parameter is outside it's bounds
			while (c_p < p_min || c_p > p_max)
			{
				if (c_p < p_min)
				{
					c_p += p_min - c_p;
				}
				if (c_p > p_max)
				{
					c_p -= c_p - p_max;
				}
			}
			// do cross over
			//the trial vector was initialized with the target vector
			// so the parameters from the parent don't need to be set
			double rand_cr = cr_prob(rand_engine);
			if (rand_cr > cr || idx == par_id_chg)
			{
				x_trial[ipar] = c_p;
			}
		}
		par_transform.numeric2model_ip(x_trial);
		run_manager.add_run(x_trial);
	}
}

int DifferentialEvolution::recombination(RunManagerAbstract &run_manager)
{
	ostream &os = file_manager.rec_ofstream();

	int best_run_idx = 0;
	ModelRun run_target(obj_func_ptr);
	ModelRun run_canidate(obj_func_ptr);
	Parameters tmp_pars_targ;
	Observations tmp_obs_targ;
	Parameters tmp_pars_can;
	Observations tmp_obs_can;

	best_phi = std::numeric_limits<double>::max();

	int d = gen_1.get_nruns();
	int n_good_runs_targ = gen_1.get_num_good_runs();
	int n_good_runs_can = run_manager.get_num_good_runs();
	double phi_sum_targ = 0.0;
	double phi_sum_can = 0.0;
	double phi_sum_new = 0.0;
	double phi_max_targ = 0.0;
	double phi_max_can = 0.0;
	double phi_max_new = 0.0;
	double phi_min_targ = std::numeric_limits<double>::max();
	double phi_min_can = std::numeric_limits<double>::max();
	double phi_min_new = std::numeric_limits<double>::max();

	os << "  population phi values:" << endl;
	os << "              parent         candidate" << endl;
	os << "    id        phi            phi" << endl;
	os << "    ----     ---------      ---------" << endl;
	for (int i_run = 0; i_run < d; ++i_run)
	{
		double new_phi = std::numeric_limits<double>::max();
		bool run_target_ok = gen_1.get_run(i_run, tmp_pars_targ, tmp_obs_targ);
		bool  run_canidate_ok = run_manager.get_run(i_run, tmp_pars_can, tmp_obs_can);
		if (!run_canidate_ok && !run_target_ok)
		{
			//keep current target
			new_phi = std::numeric_limits<double>::max();
			++failed_runs_old;
			++failed_runs_new;
			os << "    " << left << setw(10) << i_run << setw(15) << "N/A" << setw(15) << "N/A";
			os << endl;

		}
		else if (!run_canidate_ok)
		{
			//keep current target
			++failed_runs_new;
			par_transform.model2ctl_ip(tmp_pars_targ);
			// compute phi
			run_target.update_ctl(tmp_pars_targ, tmp_obs_targ);
			double phi_target = run_target.get_phi(DynamicRegularization::get_unit_reg_instance());
			phi_sum_targ += phi_target;
			phi_sum_new += phi_target;
			phi_min_targ = min(phi_min_targ, phi_target);
			phi_max_targ = max(phi_max_targ, phi_target);
			os << "    " << left << setw(10) << i_run << setw(15) << phi_target << setw(15) << "N/A";
			os << endl;
		}
		else if (!run_target_ok)
		{
			gen_1.update_run(i_run, tmp_pars_can, tmp_obs_can);
			// compute phi
			par_transform.model2ctl_ip(tmp_pars_can);
			run_canidate.update_ctl(tmp_pars_can, tmp_obs_can);
			new_phi = run_canidate.get_phi(DynamicRegularization::get_unit_reg_instance());
			++failed_runs_old;
			phi_sum_can += new_phi;
			phi_sum_new += new_phi;
			phi_min_can = min(phi_min_targ, new_phi);
			phi_max_can = max(phi_max_targ, new_phi);
			os << "    " << left << setw(10) << i_run << setw(15) << "N/A" << setw(15) << new_phi;
			os << endl;
		}
		else
		{
			// process target parameters and observations
			par_transform.model2ctl_ip(tmp_pars_targ);
			run_target.update_ctl(tmp_pars_targ, tmp_obs_targ);
			double phi_target = run_target.get_phi(DynamicRegularization::get_unit_reg_instance());
			//process candidate parameters and observations
			par_transform.model2ctl_ip(tmp_pars_can);
			run_canidate.update_ctl(tmp_pars_can, tmp_obs_can);
			double phi_canidate = run_canidate.get_phi(DynamicRegularization::get_unit_reg_instance());
			new_phi = min(phi_target, phi_canidate);
			os << "    " << left << setw(10) << i_run;
			os << setw(15) << phi_target;
			os << setw(15) << phi_canidate;
			os << endl;
			if (phi_canidate < phi_target)
			{
				gen_1.update_run(i_run, tmp_pars_can, tmp_obs_can);
			}
			phi_sum_targ += phi_target;
			phi_sum_can += phi_canidate;
			phi_sum_new += new_phi;
			phi_min_can = min(phi_min_can, phi_canidate);
			phi_max_can = max(phi_max_can, phi_canidate);
			phi_min_targ = min(phi_min_targ, phi_target);
			phi_max_targ = max(phi_max_targ, phi_target);

		}
		phi_min_new = min(phi_min_new, new_phi);
		phi_max_new = max(phi_max_new, new_phi);
		if (new_phi < best_phi)
		{
			best_phi = new_phi;
			best_run_idx = i_run;
		}
	}
	double phi_avg_targ = std::numeric_limits<double>::max();
	double phi_avg_can = std::numeric_limits<double>::max();
	double phi_avg_new = std::numeric_limits<double>::max();
	int n_good_runs_new = gen_1.get_num_good_runs();
	if (n_good_runs_targ > 0)
	{
		phi_avg_targ = phi_sum_targ / n_good_runs_targ;
	}
	if (n_good_runs_can > 0)
	{
		phi_avg_can = phi_sum_can / n_good_runs_can;
	}
	if (n_good_runs_new > 0)
	{
		phi_avg_new = phi_sum_new / n_good_runs_new;
	}

	write_run_summary(cout, n_good_runs_targ, phi_avg_targ, phi_min_targ, phi_max_targ,
		n_good_runs_can, phi_avg_can, phi_min_can, phi_max_can,
		n_good_runs_new, phi_avg_new, phi_min_new, phi_max_new);
	cout << endl;
	os << endl;
	write_run_summary(os, n_good_runs_targ, phi_avg_targ, phi_min_targ, phi_max_targ,
		n_good_runs_can, phi_avg_can, phi_min_can, phi_max_can,
		n_good_runs_new, phi_avg_new, phi_min_new, phi_max_new);
	os << endl;

	return best_run_idx;
}

void DifferentialEvolution::write_run_summary(std::ostream &os,
	int nrun_par, double avg_par, double min_par, double max_par,
	int nrun_can, double avg_can, double min_can, double max_can,
	int nrun_child, double avg_child, double min_child, double max_child)
{
	os << "  summary of differential evolution iteration:" << endl;
	os << "    group      num runs     avg phi        min phi        max phi" << endl;
	os << "    --------   --------    ---------      ---------      --------- " << endl;

		os << left << setw(17) << "    parent";
		os << setw(11) << nrun_par;
		os << setw(15) << avg_par;
		os << setw(15) << min_par;
		os << setw(15) << max_par;
		os << endl;
		os << left << setw(17) << "    candidate";
		os << setw(11) << nrun_can;
		os << setw(15) << avg_can;
		os << setw(15) << min_can;
		os << setw(15) << max_can;
		os << endl;
		os << left << setw(17) << "    child";
		os << setw(11) << nrun_child;
		os << setw(15) << avg_child;
		os << setw(15) << min_child;
		os << setw(15) << max_child;
		os << endl;
		//cout << "  summary:" << endl;
		//cout << "    parents  : runs=" << n_good_runs_targ << ";  phi(avg=" << phi_avg_targ << "; min=" << phi_min_targ << "; max=" << phi_max_targ << endl;
		//cout << "    candidates: runs=" << n_good_runs_can << ";  phi(avg=" << phi_avg_can << "; min=" << phi_min_can << "; max=" << phi_max_can << endl;
		//cout << "    children : runs=" << n_good_runs_new << ";  phi(avg=" << phi_avg_new << "; min=" << phi_min_new << "; max=" << phi_max_new << endl;

}
DifferentialEvolution::~DifferentialEvolution()
{
}

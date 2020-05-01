#include <random>
#include <iomanip>
#include "MOEA.h"
#include "Ensemble.h"
#include "RunManagerAbstract.h"
#include "ModelRunPP.h"
#include "RestartController.h"
#include <iostream>
#include <numeric>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>

using namespace mou;
using namespace std;

//ZAK: this should be done in the initialize() method
mt19937_64 MOEA::rand_engine = mt19937_64(1);

const string MOEA::solver_type_name = "MOEA";

MOEA::MOEA(Pest &_pest_scenario, FileManager &_file_manager,
	Objectives *_objs_ptr, Constraints *_cons_ptr, const ParamTransformSeq &_par_transform,
	OutputFileWriter &_output_file_writer, PerformanceLog *_performance_log, RunManagerAbstract* _run_mgr_ptr)
	: pest_scenario(_pest_scenario), file_manager(_file_manager), par_transform(_par_transform),
	output_file_writer(_output_file_writer), performance_log(_performance_log),
	run_mgr_ptr(_run_mgr_ptr),
	gen_1(_file_manager.build_filename("zak")),
	nreal(-1),
	nbin(-1),
	nobj(-1),
	ncon(-1),
	popsize(-1),
	ngen(-1),
	nreport(1),
	pcross_real(-1),
	pcross_bin(-1),
	pmut_real(-1),
	pmut_bin(-1),
	eta_c(-1),
	eta_m(-1),
	epsilon_c(EPS),
	nbits(0),
	limits_realvar(0),
	limits_binvar(0),
	function(0),
	popFunction(0),
	reportFunction(0),
	// choice(0),
	// obj1(0),
	// obj2(0),
	// obj3(0),
	// angle1(0),
	// angle2(0),
	backupFilename("nsga2_backup_pop.data"),
	nbinmut(0),
	nrealmut(0),
	nbincross(0),
	nrealcross(0),
	bitlength(0),
	parent_pop(0),
	child_pop(0),
	mixed_pop(0),
	crowd_obj(true)
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

void MOEA::solve(RunManagerAbstract &run_manager,
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
		output_file_writer.iteration_report(cout, iter+1, run_manager.get_total_runs(), "differntial evolution");
		os << endl;
		output_file_writer.iteration_report(os, iter + 1, run_manager.get_total_runs(), "differntial evolution");
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


void MOEA::initialize(RunManagerAbstract &run_manager, int d) throw (MOEAexception)
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
	// make innitial population vector model runs
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
	ModelRun tmp_run(objs_ptr);
	// copied from NSGA2.cpp
	cout << "Initializing NSGA-II v0.2.1\n"
		<< "Checking configuration" << endl;

	if (nreal < 0)
		throw MOEAexception("Invalid number of real variables");
	if (nbin < 0)
		throw MOEAexception("Invalid number of binary variables");
	if (nreal == 0 && nbin == 0)
		throw MOEAexception("Zero real and binary variables");
	if (nobj < 1)
		throw MOEAexception("Invalid number of objective functions");
	if (ncon < 0)
		throw MOEAexception("Invalid number of constraints");
	if (popsize < 4 || (popsize % 4) != 0)
		throw MOEAexception("Invalid size of population");
	if (pcross_real < 0.0 || pcross_real>1.0)
		throw MOEAexception("Invalid probability of real crossover");
	if (pmut_real < 0.0 || pmut_real>1.0)
		throw MOEAexception("Invalid probability of real mutation");
	if (pcross_bin < 0.0 || pcross_bin>1.0)
		throw MOEAexception("Invalid probability of binary crossover");
	if (pmut_bin < 0.0 || pmut_bin>1.0)
		throw MOEAexception("Invalid probability of binary mutation");
	if (eta_c <= 0)
		throw MOEAexception("Invalid distribution index for crossover");
	if (eta_m <= 0)
		throw MOEAexception("Invalid distribution index for mutation");
	if (ngen < 1)
		throw MOEAexception("Invalid number of generations");
	if (nbin != 0 && nbits.size() == 0)
		throw MOEAexception("Invalid number of bits for binary variables");
	if (limits_realvar.size() != nreal)
		throw MOEAexception("Invalid number of real variable limits");
	if (limits_binvar.size() != nbin)
		throw MOEAexception("Invalid number of binary variable limits");
	if (function == 0)
		throw MOEAexception("Evaluation function not defined");

	init_streams();
	report_parameters(fpt5);

	nbinmut = 0;
	nrealmut = 0;
	nbincross = 0;
	nrealcross = 0;
	bitlength = std::accumulate(nbits.begin(), nbits.end(), 0);

	parent_pop = new population(popsize,
		nreal,
		nbin,
		ncon,
		nbits,
		limits_realvar,
		limits_binvar,
		nobj,
		pmut_real,
		pmut_bin,
		eta_m,
		epsilon_c,
		function);
	child_pop = new population(popsize,
		nreal,
		nbin,
		ncon,
		nbits,
		limits_realvar,
		limits_binvar,
		nobj,
		pmut_real,
		pmut_bin,
		eta_m,
		epsilon_c,
		function);
	mixed_pop = new population(popsize * 2,
		nreal,
		nbin,
		ncon,
		nbits,
		limits_realvar,
		limits_binvar,
		nobj,
		pmut_real,
		pmut_bin,
		eta_m,
		epsilon_c,
		function);

	if (popFunction) {
		parent_pop->set_popfunction(popFunction);
		child_pop->set_popfunction(popFunction);
		mixed_pop->set_popfunction(popFunction);
	}

	parent_pop->crowd_obj = crowd_obj;
	child_pop->crowd_obj = crowd_obj;
	mixed_pop->crowd_obj = crowd_obj;

	//randomize();

	bool fromBackup = load_backup();
	if (!fromBackup) {
		parent_pop->initialize();
		cout << "Initialization done, now performing first generation" << endl;

		parent_pop->decode();
		parent_pop->custom_evaluate();
		parent_pop->fast_nds();
		parent_pop->crowding_distance_all();

		t = 1;
	}
	else {
		cout << "Initialization made from backup file" << endl;
	}

	custom_report(*parent_pop);

	report_pop(*parent_pop, fpt1);
	fpt4 << "# gen = " << t << '\n';
	report_pop(*parent_pop, fpt4);

	fpt1.flush();
	fpt4.flush();
	fpt5.flush();
	// end copy
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

void MOEA::initialize_vector(Parameters &numeric_pars)
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

void MOEA::mutation(RunManagerAbstract &run_manager, double f, bool dither_f, double cr)
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
		//initialize trail vector with the traget vector
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

int MOEA::recombination(RunManagerAbstract &run_manager)
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
	os << "              parent         canidate" << endl;
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
			//process canidate parameters and observations
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

void MOEA::write_run_summary(std::ostream &os,
	int nrun_par, double avg_par, double min_par, double max_par,
	int nrun_can, double avg_can, double min_can, double max_can,
	int nrun_child, double avg_child, double min_child, double max_child)
{
	os << "  summary of differntial evolution iteration:" << endl;
	os << "    group      num runs     avg phi        min phi        max phi" << endl;
	os << "    --------   --------    ---------      ---------      --------- " << endl;

		os << left << setw(17) << "    parent";
		os << setw(11) << nrun_par;
		os << setw(15) << avg_par;
		os << setw(15) << min_par;
		os << setw(15) << max_par;
		os << endl;
		os << left << setw(17) << "    canidate";
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
		//cout << "    canidates: runs=" << n_good_runs_can << ";  phi(avg=" << phi_avg_can << "; min=" << phi_min_can << "; max=" << phi_max_can << endl;
		//cout << "    children : runs=" << n_good_runs_new << ";  phi(avg=" << phi_avg_new << "; min=" << phi_min_new << "; max=" << phi_max_new << endl;
}

void MOEA::advance() {

	cout << "Advancing to generation " << t + 1 << endl;

	std::pair<int, int> res;

	// create next population Qt
	selection(*parent_pop, *child_pop);
	res = child_pop->mutate();
	child_pop->generation = t + 1;
	child_pop->decode();
	child_pop->custom_evaluate();

	// mutation book-keeping
	nrealmut += res.first;
	nbinmut += res.second;

	// fpt4 << "#Child pop\n";
	// report_pop(*child_pop,fpt4);

	// create population Rt = Pt U Qt
	mixed_pop->merge(*parent_pop, *child_pop);
	mixed_pop->generation = t + 1;

	// fpt4 << "#Mixed\n";
	// report_pop(*mixed_pop, fpt4);

	mixed_pop->fast_nds();
	//mixed_pop->crowding_distance_all();

	// fpt4 << "#Mixed nfs\n";
	// report_pop(*mixed_pop, fpt4);


	// Pt+1 = empty
	parent_pop->ind.clear();

	int i = 0;
	// until |Pt+1| + |Fi| <= N, i.e. until parent population is filled
	while (parent_pop->size() + mixed_pop->front[i].size() < popsize) {
		std::vector<int>& Fi = mixed_pop->front[i];
		mixed_pop->crowding_distance(i);           // calculate crowding in Fi
		for (int j = 0; j < Fi.size(); ++j)        // Pt+1 = Pt+1 U Fi
			parent_pop->ind.push_back(mixed_pop->ind[Fi[j]]);
		i += 1;
	}

	mixed_pop->crowding_distance(i);           // calculate crowding in Fi
	std::sort(mixed_pop->front[i].begin(),
		mixed_pop->front[i].end(),
		sort_n(*mixed_pop));// sort remaining front using <n

	const int extra = popsize - parent_pop->size();
	for (int j = 0; j < extra; ++j) // Pt+1 = Pt+1 U Fi[1:N-|Pt+1|]
		parent_pop->ind.push_back(mixed_pop->ind[mixed_pop->front[i][j]]);

	t += 1;

	// if (popFunction) {
	//   (*popFunction)(*parent_pop);
	// }

	parent_pop->generation = t;
	custom_report(*parent_pop);

	if (t % nreport == 0) {
		fpt4 << "# gen = " << t << '\n';
		report_pop(*parent_pop, fpt4);
		fpt4.flush();
	}

	// save a backup
	save_backup();
}

void MOEA::evolve() {
	while (t < ngen)
		advance();
	report_pop(*parent_pop, fpt2);
}

MOEA::~MOEA()
{
	if (parent_pop) {
		delete parent_pop;
		parent_pop = 0;
	}
	if (child_pop) {
		delete child_pop;
		child_pop = 0;
	}
	if (mixed_pop) {
		delete mixed_pop;
		mixed_pop = 0;
	}
}

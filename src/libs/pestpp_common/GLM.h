#ifndef GLM_H_
#define GLM_H_

#include <map>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "covariance.h"
#include "FileManager.h"
#include "Jacobian_1to1.h"
#include "ModelRunPP.h"
#include "ObjectiveFunc.h"
#include "OutputFileWriter.h"
#include "PerformanceLog.h"
#include "Pest.h"
#include "RestartController.h"
#include "RunManagerAbstract.h"
#include "SVDSolver.h"
#include "TerminationController.h"

/// pestpp-glm's top-level loop, as a class - the peer of MOEA, SeqQuadProgram and
/// EnsembleMethod.
///
/// The algorithm is not here. SVDSolver already owns the Gauss-Levenberg-Marquardt maths, and
/// this does not touch it. What lives here is everything that used to sit loose in
/// pest++.cpp's main(): building the solver and its collaborators, the noptmax==0 and
/// noptmax<0 special cases, the restart branching, the iteration loop, the termination
/// summary and the posterior FOSM block. main() is left with argument parsing and scenario
/// construction, the same as the other tools' mains.
///
/// The shape is deliberately identical to the peers - initialize() / iterate_2_solution() /
/// finalize(), with initialize() also available split in two - because that is the shape the
/// C ABI's ToolAdapter expects. glm has no adapter yet; this is what makes one writable
/// without touching the algorithm again.
///
/// Two things follow from being callable through a library rather than only from a main():
///
///  1. NOTHING here calls exit(). The loose code did, in the noptmax==0 path and in the
///     iteration handler, which under the API would take the host process down with it.
///     Failures throw, and the caller decides.
///  2. State that an iteration needs is held as members, not locals, so a caller can drive
///     one iteration at a time and look at what happened in between.
class GLM
{
public:
	/// The first five arguments are exactly the peers' constructor, so an adapter can build
	/// this the way it builds the others. Restart is glm's own and is decided from the
	/// command line before any tool exists, so it arrives separately and is optional -
	/// a library caller passes nothing and gets a session with no restart.
	GLM(Pest& _pest_scenario, FileManager& _file_manager,
		OutputFileWriter& _output_file_writer, PerformanceLog* _performance_log,
		RunManagerAbstract* _run_mgr_ptr,
		RestartController* _restart_ctl = nullptr, bool _restart_flag = false,
		bool _save_restart_rec_header = true);

	/// Build the solver, load parcov, prime the run manager, and handle noptmax==0 (which is
	/// a single forward run and a report, not an iteration).
	void initialize();

	/// initialize(), split so a caller can own the initial run. Returns the number of runs to
	/// service before initialize_finish() - 1 for the noptmax==0 forward run, 0 otherwise,
	/// since every other path computes its first batch inside the first iteration.
	int  initialize_prepare();
	void initialize_finish();

	/// Iterate until the termination controller says stop. The in-tree composition of
	/// solve_iteration().
	void iterate_2_solution();

	/// One iteration - the loop body of iterate_2_solution(), including the restart branching
	/// and the termination check. Returns whether the loop should continue, so a caller
	/// driving its own loop runs the same steps in the same order rather than an
	/// approximation of them.
	bool solve_iteration();

	/// Termination summary, final parameter/residual reporting, and the posterior FOSM
	/// calculations. Safe to call after a failed iteration: it reports what there is.
	void finalize();

	/// Report to the record file and throw. The record stream is FLUSHED, never closed:
	/// closing it destroys the ofstream, and a later write - including one from a destructor
	/// during unwinding - then terminates the process. That is not hypothetical; it is what
	/// the equivalent in sqp, mou and ies used to do.
	void throw_glm_error(const std::string& message);

	// -- what a caller wants to look at between iterations ---------------------------------
	const ModelRun& get_optimum_run() const { return optimum_run; }
	const ModelRun& get_current_run() const { return cur_run; }
	TerminationController& get_termination_ctl() { return termination_ctl; }
	Jacobian& get_base_jacobian() { return *base_jacobian_ptr; }
	/// Whether the Jacobian exists yet. It is built in initialize_prepare(), so a caller that
	/// asks before initialize() gets a clear answer instead of a dereferenced null.
	bool has_base_jacobian() const { return base_jacobian_ptr != nullptr; }

	/// The parameter vector the next iteration will work from, in CTL space, and the best one
	/// found so far. Both live in their ModelRun, so these read through rather than copy.
	/// Setting the current one is the point of exposing it: a caller can move the starting
	/// position between iterations.
	const Parameters& get_current_parameters() const { return cur_run.get_ctl_pars(); }
	const Parameters& get_optimum_parameters() const { return optimum_run.get_ctl_pars(); }
	void set_current_parameters(const Parameters& pars) { cur_run.set_ctl_parameters(pars); }
	int get_iteration_number() { return termination_ctl.get_iteration_number(); }

private:
	Pest& pest_scenario;
	FileManager& file_manager;
	OutputFileWriter& output_file_writer;
	PerformanceLog* performance_log;
	RunManagerAbstract* run_mgr_ptr;

	/// Not owned. Null when there is no restart, which is every library caller.
	RestartController* restart_ctl;
	bool restart_flag;
	bool save_restart_rec_header;

	ObjectiveFunc obj_func;
	std::mt19937 rand_gen;
	Covariance parcov;
	TerminationController termination_ctl;

	/// unique_ptr rather than a bare new/delete pair: the loose version leaked the Jacobian
	/// on every error path out of the iteration loop.
	std::unique_ptr<Jacobian> base_jacobian_ptr;
	std::unique_ptr<SVDSolver> base_svd_ptr;

	ModelRun cur_run;
	ModelRun optimum_run;
	Parameters cur_ctl_parameters;
	int noptmax;
	bool initialized;

	/// noptmax==0: one forward run with the initial parameters, reported and then done.
	void run_initial_parameters_only();
	/// The posterior Schur-complement block, and the wall of explanatory text that goes with it.
	void write_fosm();
};

#endif /* GLM_H_ */

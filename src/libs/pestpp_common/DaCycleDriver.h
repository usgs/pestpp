#ifndef DA_CYCLE_DRIVER_H_
#define DA_CYCLE_DRIVER_H_

#include <map>
#include <set>
#include <memory>
#include <string>
#include <vector>
#include <fstream>

#include "DataAssimilator.h"
#include "Ensemble.h"
#include "EnsembleMethodUtils.h"
#include "FileManager.h"
#include "Localizer.h"
#include "OutputFileWriter.h"
#include "PerformanceLog.h"
#include "Pest.h"
#include "RunManagerAbstract.h"

/// Drives pestpp-da's SEQUENCE of assimilation cycles.
///
/// A DataAssimilator handles ONE cycle. What turns that into data assimilation is the sequence:
/// each cycle gets its own child scenario, its own run-manager parameter/observation set and its
/// own DataAssimilator, and the POSTERIOR ensemble of cycle N becomes the PRIOR of cycle N+1 -
/// with failed realizations dropped from the global ensembles along the way so the two stay
/// aligned. That machinery lived in pestpp-da.cpp's main(), which is why the C ABI could only
/// ever run the first cycle and a multi-cycle control file silently produced a one-cycle answer.
///
/// It lives here so main() and the API drive the SAME sequence. That is not a stylistic
/// preference: every precondition left in a main() this month - set_default_dynreg(),
/// prep_glm_files(), write_par_iter(), the regularization-mode guards - turned into a bug the
/// moment the library gained a second entry point, because a check only one caller performs is
/// not a check.
class DaCycleDriver
{
public:
	/// Which run manager is in play. Only two behaviours differ per cycle, so this takes an
	/// enum rather than CmdLine - the driver has no business knowing about argv.
	enum class RunManagerKind { SERIAL, PANTHER_MASTER, EXTERNAL, OTHER };

	DaCycleDriver(Pest& _parent_scenario, FileManager& _file_manager,
		PerformanceLog* _performance_log, RunManagerAbstract* _run_mgr_ptr,
		RunManagerKind _rm_kind = RunManagerKind::OTHER, bool _restart_flag = false);

	/// Where a per-cycle RunManagerSerial should run the model. Defaults to ".", which is what
	/// every in-tree caller uses; set it before initialize() if the caller ran from elsewhere.
	void set_pathname(const std::string& _pathname);

	/// Build the global ensembles, resolve the cycle list and the per-cycle noptmax schedule,
	/// and initialize the localizer. After this get_cycles() is meaningful.
	void initialize();

	/// Every assimilation cycle in the control file, in the order they will run.
	const std::vector<int>& get_cycles() const { return assimilation_cycles; }
	/// The cycle currently set up, or -1 before the first begin_cycle().
	int get_current_cycle() const { return current_cycle; }
	/// How many cycles remain, including the current one.
	int get_n_cycles_remaining() const;

	/// Set up the next cycle: child scenario, run manager, a fresh DataAssimilator, and the
	/// carried ensembles restricted to this cycle's names. False when the sequence is done.
	/// Cycles below ++da_start_cycle are skipped here rather than by the caller.
	bool begin_cycle();

	/// The cycle's tool, valid between begin_cycle() and end_cycle(). Null outside that.
	/// This is what an API caller drives - initialize(), solve, and so on - so the per-cycle
	/// surface is the ordinary DataAssimilator one and needs nothing new.
	DataAssimilator* get_da() { return da.get(); }

	/// The cycle's OWN scenario, valid between begin_cycle() and end_cycle(). A cycle has its
	/// own parameter and observation sets - that is what a cycle IS - so a caller asking "what
	/// parameters am I looking at" must be answered from here, not from the parent.
	Pest* get_child_scenario() { return child_scenario.get(); }

	/// Run the cycle that begin_cycle() set up: initialize the tool, assimilate, report phi.
	/// Split from begin_cycle() so a caller can instead take get_da() and drive the cycle
	/// itself - which is the whole point of exposing the sequence rather than just running it.
	void drive_cycle();

	/// Harvest this cycle's posterior back into the global ensembles: drop realizations the
	/// cycle lost, realign, and write the cycle's columns back. This is the step that makes it
	/// assimilation rather than a loop, and the reason a caller cannot just run cycles itself.
	void end_cycle();

	/// initialize() then begin/drive/end for every cycle - the in-tree composition, so the
	/// executable and a caller stepping cycles run the same code.
	void run_all_cycles();

	void finalize();

	// -- the global ensembles, carried across cycles -----------------------------------------
	ParameterEnsemble& get_global_pe() { return curr_pe; }
	ObservationEnsemble& get_global_oe() { return curr_oe; }
	ObservationEnsemble& get_global_noise() { return curr_noise; }

private:
	Pest& parent_scenario;
	FileManager& file_manager;
	PerformanceLog* performance_log;
	RunManagerAbstract* run_mgr_ptr;
	RunManagerKind rm_kind;
	bool restart_flag;

	std::vector<int> assimilation_cycles;
	std::map<int, int> noptmax_schedule;
	std::vector<std::string> init_real_names;
	int start_cycle;
	int cycle_index;      ///< index into assimilation_cycles of the cycle set up now
	int current_cycle;

	/// The ensembles that survive across cycles: cycle N's posterior is cycle N+1's prior.
	ParameterEnsemble curr_pe;
	ObservationEnsemble curr_oe;
	ObservationEnsemble curr_noise;

	/// Per-cycle objects, rebuilt by begin_cycle() and torn down by end_cycle().
	std::unique_ptr<Pest> child_scenario;
	std::unique_ptr<OutputFileWriter> child_ofw;
	std::unique_ptr<DataAssimilator> da;

	std::unique_ptr<Localizer> global_loc;
	bool initialized;

	/// The cycle tables, read once by initialize(). A cycle absent from a table simply has no
	/// entry - the tables are sparse, and obs_in_tbl is what lets a cycle zero the weight of an
	/// observation the table lists for OTHER cycles but not this one.
	std::vector<int> cycles_in_tables;
	std::map<int, std::map<std::string, double>> par_cycle_info;
	std::map<int, std::map<std::string, double>> obs_cycle_info;
	std::map<int, std::map<std::string, double>> weight_cycle_info;
	std::set<std::string> obs_in_tbl;

	/// The parent scenario's writer. A MEMBER, not a local in initialize(): ParChangeSummarizer
	/// keeps an OutputFileWriter* and would otherwise point at a destroyed stack object - the
	/// same defect that made sequentialLP crash the first time it was driven through the ABI.
	std::unique_ptr<OutputFileWriter> parent_ofw;
	/// Built once curr_pe exists, because it holds a pointer to it.
	std::unique_ptr<ParChangeSummarizer> pcs;

	std::string pathname;   ///< passed to a per-cycle RunManagerSerial
	int verbose_level;

	/// Set by end_cycle() when da_stop_cycle or pest.stp says to stop. begin_cycle() reports
	/// false on the next call, so the stop looks like the sequence ending rather than a break
	/// out of somebody's loop - which is what makes the API and the executable agree.
	bool stop_requested;

	/// True when begin_cycle() created the run manager for this cycle and end_cycle() must
	/// destroy it. Only the SERIAL path does; every other kind reuses the caller's.
	bool owns_cycle_run_mgr;
	/// Whether the cycle's simulated outputs can be reused. Computed in begin_cycle() because it
	/// compares against the PREVIOUS cycle, consumed by drive_cycle().
	bool use_existing;
	int cycle_nnz_obs;

	/// Live only between begin_cycle() and end_cycle(): the cycle's view of the global
	/// ensembles. end_cycle() reads them back out, which is why they outlive begin_cycle()'s
	/// scope rather than being locals as they were in main().
	ParameterEnsemble cycle_curr_pe;
	ObservationEnsemble cycle_curr_oe;
	/// Carried across cycles so a cycle can see what the previous one ended at.
	Parameters prev_cycle_ctl_pars;
	Observations prev_cycle_obs;
	std::ofstream f_phi;   ///< the global phi log, opened once by initialize()

	void drop_missing_realizations(ParameterEnsemble& cycle_pe);
	void drop_missing_realizations(ObservationEnsemble& cycle_oe);
};

#endif /* DA_CYCLE_DRIVER_H_ */

#ifndef DA_CYCLE_DRIVER_H_
#define DA_CYCLE_DRIVER_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "DataAssimilator.h"
#include "Ensemble.h"
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

	void drop_missing_realizations(ParameterEnsemble& cycle_pe);
	void drop_missing_realizations(ObservationEnsemble& cycle_oe);
};

#endif /* DA_CYCLE_DRIVER_H_ */

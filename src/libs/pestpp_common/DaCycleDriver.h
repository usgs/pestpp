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

/// runs pestpp-da's sequence of assimilation cycles.
///
/// a DataAssimilator only does one cycle. what makes it data assimilation is the sequence - each
/// cycle gets its own child pest object, its own par/obs set for the run manager and its own
/// DataAssimilator, and the posterior ensemble from cycle N becomes the prior for cycle N+1.
/// any realizations that fail along the way get dropped from the global ensembles too so the
/// two stay lined up. all of that used to live in main() in pestpp-da.cpp, which is why the c
/// api could only ever run the first cycle - a control file with several cycles would quietly
/// give you a one cycle answer.
///
/// it lives here so main() and the api run the same code. we kept getting bit by this: every
/// setup step left in a main() this month - set_default_dynreg(), prep_glm_files(),
/// write_par_iter(), the regularization mode checks - turned into a bug as soon as the library
/// got a second entry point. a check that only one caller does isnt really a check.
class DaCycleDriver
{
public:
	/// which run manager is being used. only a couple of things differ per cycle, so this takes
	/// an enum instead of CmdLine - no reason for this class to know anything about argv.
	enum class RunManagerKind { SERIAL, PANTHER_MASTER, EXTERNAL, OTHER };

	DaCycleDriver(Pest& _parent_scenario, FileManager& _file_manager,
		PerformanceLog* _performance_log, RunManagerAbstract* _run_mgr_ptr,
		RunManagerKind _rm_kind = RunManagerKind::OTHER, bool _restart_flag = false);

	/// where a per-cycle RunManagerSerial should run the model. defaults to ".", which is what
	/// everything in the tree uses. set it before initialize() if you are running from
	/// somewhere else.
	void set_pathname(const std::string& _pathname);

	/// build the global ensembles, work out the cycle list and the per-cycle noptmax schedule,
	/// and set up the localizer. after this get_cycles() means something.
	void initialize();

	/// every assimilation cycle in the control file, in the order they run.
	const std::vector<int>& get_cycles() const { return assimilation_cycles; }
	/// the cycle that is set up now, or -1 before the first begin_cycle().
	int get_current_cycle() const { return current_cycle; }
	/// how many cycles are left, counting the current one.
	int get_n_cycles_remaining() const;

	/// set up the next cycle - child pest object, run manager, a new DataAssimilator, and the
	/// carried ensembles cut down to this cycle's names. returns false when there are no cycles
	/// left. cycles before ++da_start_cycle get skipped in here instead of by the caller.
	bool begin_cycle();

	/// the cycle's tool. only valid between begin_cycle() and end_cycle(), null otherwise.
	/// this is what an api caller drives - initialize(), solve and so on - so driving a cycle
	/// is just the normal DataAssimilator calls and doesnt need anything new.
	DataAssimilator* get_da() { return da.get(); }

	/// the cycle's own pest object, valid between begin_cycle() and end_cycle(). each cycle has
	/// its own parameters and observations - that is what a cycle is - so if you want to know
	/// what parameters you are looking at, the answer comes from here, not the parent.
	Pest* get_child_scenario() { return child_scenario.get(); }

	/// run the cycle that begin_cycle() set up - initialize the tool, assimilate, report phi.
	/// kept separate from begin_cycle() so you can grab get_da() and drive the cycle yourself
	/// instead, which is the whole point of exposing the sequence at all.
	void drive_cycle();

	/// pull this cycle's posterior back into the global ensembles - drop any realizations the
	/// cycle lost, line things back up, and write the cycle's columns back. this is the step
	/// that makes it assimilation instead of just a loop, and why you cant simply run the
	/// cycles yourself.
	void end_cycle();

	/// initialize() then begin/drive/end for every cycle. this is what the exe calls, so the exe
	/// and somebody stepping cycles by hand run the same code.
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

	/// the ensembles that carry across cycles - cycle N's posterior is cycle N+1's prior.
	ParameterEnsemble curr_pe;
	ObservationEnsemble curr_oe;
	ObservationEnsemble curr_noise;

	/// per-cycle stuff, rebuilt by begin_cycle() and thrown away by end_cycle().
	std::unique_ptr<Pest> child_scenario;
	std::unique_ptr<OutputFileWriter> child_ofw;
	std::unique_ptr<DataAssimilator> da;

	std::unique_ptr<Localizer> global_loc;
	bool initialized;

	/// the cycle tables, read once in initialize(). the tables are sparse - a cycle that isnt in
	/// a table just has no entry. obs_in_tbl is what lets a cycle zero the weight on an obs that
	/// the table lists for other cycles but not this one.
	std::vector<int> cycles_in_tables;
	std::map<int, std::map<std::string, double>> par_cycle_info;
	std::map<int, std::map<std::string, double>> obs_cycle_info;
	std::map<int, std::map<std::string, double>> weight_cycle_info;
	std::set<std::string> obs_in_tbl;

	/// the parent pest object's writer. this is a member instead of a local in initialize()
	/// because ParChangeSummarizer keeps an OutputFileWriter* - if it were a local it would end
	/// up pointing at something that got destroyed. same thing that made sequentialLP crash the
	/// first time it got driven through the api.
	std::unique_ptr<OutputFileWriter> parent_ofw;
	/// built after curr_pe exists, since it holds a pointer to it.
	std::unique_ptr<ParChangeSummarizer> pcs;

	std::string pathname;   ///< handed to a per-cycle RunManagerSerial
	int verbose_level;

	/// set by end_cycle() when da_stop_cycle or pest.stp says to quit. the next begin_cycle()
	/// returns false, so stopping looks the same as running out of cycles instead of somebody
	/// breaking out of a loop. that is what keeps the api and the exe doing the same thing.
	bool stop_requested;

	/// true when begin_cycle() made the run manager for this cycle and end_cycle() has to delete
	/// it. only the serial path does that - everything else reuses the caller's.
	bool owns_cycle_run_mgr;
	/// whether we can reuse the simulated outputs from last cycle. worked out in begin_cycle()
	/// since it compares against the previous cycle, then used by drive_cycle().
	bool use_existing;
	int cycle_nnz_obs;

	/// only alive between begin_cycle() and end_cycle() - this cycle's view of the global
	/// ensembles. end_cycle() reads them back out, which is why they are members here instead of
	/// locals like they were in main().
	ParameterEnsemble cycle_curr_pe;
	ObservationEnsemble cycle_curr_oe;
	/// carried across cycles so a cycle can see where the last one ended up.
	Parameters prev_cycle_ctl_pars;
	Observations prev_cycle_obs;
	std::ofstream f_phi;   ///< the global phi file, opened once in initialize()

	void drop_missing_realizations(ParameterEnsemble& cycle_pe);
	void drop_missing_realizations(ObservationEnsemble& cycle_oe);
};

#endif /* DA_CYCLE_DRIVER_H_ */

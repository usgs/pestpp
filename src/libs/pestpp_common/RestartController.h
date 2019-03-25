#ifndef RESTART_CTL_H_
#define RESTART_CTL_H_

#include <string>
#include <fstream>
#include <iostream>

class TerminationController;
class SVDSolver;
class Parameters;

class RestartController
{
public:
	enum class RestartOption {NONE, REUSE_JACOBIAN, RESUME_NEW_ITERATION, RESUME_JACOBIAN_RUNS, RESUME_UPGRADE_RUNS};
	enum class IterationType{BASE, SUPER};
	RestartController(void);
	static void write_start_failed_super(std::ostream &fout);
	static void write_start_iteration(std::ostream &fout, const std::string &solver_type, int _iter_num, int _global_iter_num);
	static void write_start_parameters_updated(std::ostream &fout, const std::string &parameter_filename);
	static void write_finish_parameters_updated(std::ostream &fout, const std::string &parameter_filename);
	static void write_jac_runs_built(std::ostream &fout);
	static void write_upgrade_runs_built(std::ostream &fout);
	static void write_iteration_complete(std::ostream &fout);
	RestartOption get_restart_option() const { return restart_option; }
	RestartOption& get_restart_option() { return restart_option; }
	IterationType get_iteration_type() { return iteration_type; }
	void process_rst_file(std::ifstream &fin);
	void update_termination_ctl(TerminationController &term_ctl);
	Parameters get_restart_parameters(const std::string &restart_par_file, const std::string &prev_par_file);
	~RestartController(void);
private:
	enum class PARAMETER_STATE {INIT_PAR, RESTART_PAR, PREV_PAR};
	PARAMETER_STATE parameter_state;
	int global_iter_no;
	int local_iter_no;
	RestartOption restart_option;
	IterationType iteration_type;
	int nopt_count;
	int nphinored_count;
	int nrelpar_count;
	std::vector<double> lowest_phi;
	std::string best_par_file;
};

#endif /* RESTART_CTL_H_ */

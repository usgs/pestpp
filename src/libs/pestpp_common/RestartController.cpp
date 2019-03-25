#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "RestartController.h"
#include "TerminationController.h"
#include "utilities.h"
#include "SVDASolver.h"

using namespace std;
using namespace::pest_utils;

RestartController::RestartController(void)
	: global_iter_no(0), local_iter_no(0), restart_option(RestartOption::NONE),
	iteration_type(IterationType::BASE), nopt_count(0), nphinored_count(0),
	nrelpar_count(), parameter_state(PARAMETER_STATE::INIT_PAR), best_par_file("")
{
}

void RestartController::write_start_iteration(ostream &fout, const string &solver_type, int _iter_num, int _global_iter_num)
{
	fout << "start_iteration " << _iter_num << "  " << _global_iter_num << "  " << solver_type << endl;

}

void RestartController::write_start_failed_super(ostream &fout)
{
	fout << "failed_super" << endl;

}

void RestartController::write_start_parameters_updated(ostream &fout, const string &parameter_filename)
{
	fout << "parameter_file_save_started " << parameter_filename << endl;
}


void RestartController::write_finish_parameters_updated(ostream &fout, const string &parameter_filename)
{
	fout << "parameter_file_save_finished " << parameter_filename << endl;
}

void RestartController::write_jac_runs_built(ostream &fout)
{
	fout << "jacobian_model_runs_built" << endl;
}

void RestartController::write_upgrade_runs_built(ostream &fout)
{
	fout << "upgrade_model_runs_built" << endl;
}

void RestartController::write_iteration_complete(ostream &fout)
{
	fout << "iteration_complete" << endl;
}

void RestartController::process_rst_file(std::ifstream &fin)
{
	string line;
	vector<string> tokens;
	while (getline(fin, line))
	{
		tokens.clear();
		tokenize(line, tokens);

		if (tokens.empty())
		{
		}
		else if (tokens[0] == "start_iteration")
		{
			convert_ip(tokens[1], local_iter_no);
			convert_ip(tokens[2], global_iter_no);
			if (tokens[3] == "svd_base_par")
			{
				iteration_type = IterationType::BASE;
			}
			else if (tokens[3] == "svda_super_par")
			{
				iteration_type = IterationType::SUPER;
			}
			restart_option = RestartOption::RESUME_NEW_ITERATION;
		}
		else if (tokens[0] == "parameter_file_save_started")
		{
			if (parameter_state != PARAMETER_STATE::INIT_PAR)
			{
				parameter_state = PARAMETER_STATE::PREV_PAR;
			}
		}
		else if (tokens[0] == "parameter_file_save_finished")
		{
			parameter_state = PARAMETER_STATE::RESTART_PAR;
		}
		else if (tokens[0] == "termination_info_1")
		{
			//tokens[1] is noptmax;
			convert_ip(tokens[2], nopt_count);
			//tokens[3] is nphinored);
			convert_ip(tokens[4], nphinored_count);
			//tokens[5] is term_ctl.nrelpar);

		}
		else if (tokens[0] == "termination_info_2")
		{
			convert_ip(tokens[1], nrelpar_count);
			//tokens[2] is nphistp;
			//tokens[3] is phiredstp;
			//tokens[4] is relparstp;

		}
		else if (tokens[0] == "termination_info_3")
		{
			lowest_phi.clear();
			for (size_t i = 1; i<tokens.size(); ++i)
			{
				double val = convert_cp<double>(tokens[i]);
				lowest_phi.push_back(val);
			}
		}
		else if (tokens[0] == "jacobian_model_runs_built")
		{
			restart_option = RestartOption::RESUME_JACOBIAN_RUNS;
		}
		else if (tokens[0] == "upgrade_model_runs_built")
		{
			restart_option = RestartOption::RESUME_UPGRADE_RUNS;
		}
		else if (tokens[0] == "failed_super")
		{
			restart_option = RestartOption::RESUME_NEW_ITERATION;
			iteration_type = IterationType::BASE;
		}
	}
}

Parameters RestartController::get_restart_parameters(const string &restart_par_file, const string &prev_par_file)
{
	Parameters new_pars;
	map<string, double> offset;
	map<string, double> scale;
	if (parameter_state == PARAMETER_STATE::INIT_PAR)
	{}
	else if (parameter_state == PARAMETER_STATE::RESTART_PAR)
	{
		ifstream fin(restart_par_file);
		new_pars.read_par_file(fin, offset, scale);
	}
	else if (parameter_state == PARAMETER_STATE::PREV_PAR)
	{
		ifstream fin(prev_par_file);
		new_pars.read_par_file(fin, offset, scale);
	}
	return new_pars;
}


void RestartController::update_termination_ctl(TerminationController &term_ctl)
{
	term_ctl.nopt_count = nopt_count;
	term_ctl.nphinored_count = nphinored_count;
	term_ctl.nrelpar_count = nrelpar_count;
	term_ctl.lowest_phi = lowest_phi;
}

RestartController::~RestartController(void)
{
}

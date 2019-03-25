#include <algorithm>
#include "RunManagerCWrapper.h"
#include "utilities.h"
#include "RunManagerPanther.h"
#include "RunManagerSerial.h"
#include "RunManagerGenie.h"
#include "pest_error.h"

typedef class RunManagerAbstract RunManagerAbstract;

static ofstream fout_run_manager_log_file;

using namespace pest_utils;

extern "C"
{

RunManager* rmic_create_serial(char **comline, int comline_array_len,
	char **tpl, int tpl_array_len,
	char **inp, int inp_array_len,
	char **ins, int ins_array_len,
	char **out, int out_array_len,
	char *storfile,
	char *rundir,
	int n_max_fail)
{
	RunManager *run_manager_ptr = nullptr;
	vector<string> comline_vec(comline, comline+comline_array_len);
	vector<string> tpl_vec(tpl, tpl+tpl_array_len);
	vector<string> inp_vec(inp, inp+inp_array_len);
	vector<string> ins_vec(ins, ins+ins_array_len);
	vector<string> out_vec(out, out+out_array_len);
	run_manager_ptr = new RunManagerSerial(comline_vec, tpl_vec, inp_vec, ins_vec,
		out_vec, storfile, rundir, n_max_fail);
	return run_manager_ptr;
}


RunManager* rmic_create_panther(
	char *storfile,
	char *port,
	char *info_filename,
	int n_max_fail,
	double overdue_reched_fac,
	double overdue_giveup_fac,
	double overdue_giveup_minutes)
{
	RunManager *run_manager_ptr = nullptr;
	fout_run_manager_log_file.open(info_filename);
	run_manager_ptr = new RunManagerPanther(storfile, port,
		fout_run_manager_log_file, n_max_fail, overdue_reched_fac,
		overdue_giveup_fac, overdue_giveup_minutes);
	return run_manager_ptr;
}

RunManager* rmic_create_genie(char **comline, int comline_array_len,
	char **tpl, int tpl_array_len,
	char **inp, int inp_array_len,
	char **ins, int ins_array_len,
	char **out, int out_array_len,
	char *storfile,
	char *genie_tag)
{
	RunManager *run_manager_ptr = nullptr;
	vector<string> comline_vec(comline, comline+comline_array_len);
	vector<string> tpl_vec(tpl, tpl+tpl_array_len);
	vector<string> inp_vec(inp, inp+inp_array_len);
	vector<string> ins_vec(ins, ins+ins_array_len);
	vector<string> out_vec(out, out+out_array_len);
	run_manager_ptr = new RunManagerGenie(comline_vec, tpl_vec, inp_vec, ins_vec, out_vec, storfile, genie_tag);
	return run_manager_ptr;
}

int rmic_initialize(RunManager *run_manager_ptr,
	char **pname, int pname_array_len,
	char **oname, int oname_array_len)

{
	int err = 0;
	try {
		vector<string> pname_vec(pname, pname+pname_array_len);
		vector<string> oname_vec(oname, oname+oname_array_len);
		run_manager_ptr->initialize(pname_vec, oname_vec);
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}

int rmic_reinitialize(RunManager *run_manager_ptr)
{
    int err = 0;
	try {
		run_manager_ptr->reinitialize();
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}


int rmic_add_run(RunManager *run_manager_ptr, double *parameter_data, int npar, int *id)
{
	int err = 0;
	try {
		vector<double> data(parameter_data, parameter_data+npar);
		*id = run_manager_ptr->add_run(data);
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}



int rmic_run(RunManager *run_manager_ptr)
{
	int err = 0;
	try {
		run_manager_ptr->run();
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}

int rmic_run_until_(RunManager *run_manager_ptr, int *condition, int no_ops, double time_sec, int *return_cond)
{
	int err = 0;
	RunManagerAbstract::RUN_UNTIL_COND enum_input_cond;
	RunManagerAbstract::RUN_UNTIL_COND enum_return_cond;
	enum_input_cond = static_cast<RunManagerAbstract::RUN_UNTIL_COND>(*condition);
	try {
		enum_return_cond = run_manager_ptr->run_until(enum_input_cond, no_ops, time_sec);
	}
	catch (...)
	{
		err = 1;
	}
	*return_cond = static_cast<int>(enum_return_cond);
	return err;
}


int rmic_get_run(RunManager *run_manager_ptr, int run_id, double *parameter_data, int npar, double *obs_data, int nobs)
{
	int err = 1;
	try
	{
	   bool success;
	   success = run_manager_ptr->get_run(run_id, parameter_data, npar, obs_data, nobs);
	   if (success) err = 0;
	}
	catch(PestIndexError ex) {
		cerr << ex.what() << endl;
		err = 1;
	}
	return err;
}

int rmic_get_num_failed_runs(RunManager *run_manager_ptr, int *nfail)
{
	int err = 0;
	*nfail = -999;
	   try
	{
        const std::set<int> &fail_set = run_manager_ptr->get_failed_run_ids();
        *nfail = fail_set.size();
	}
    catch(PestIndexError ex)
	{
		cerr << ex.what() << endl;
		err = 1;
	}
	return err;
}

int rmic_get_failed_runs_alloc(RunManager *run_manager_ptr, int *run_id_array, int *nfail)
{
    int err = 0;
    int n_failed;
    try
	{
        const std::set<int> &fail_set = run_manager_ptr->get_failed_run_ids();
        n_failed = fail_set.size();
        run_id_array = new int[n_failed];
        std::copy_n(fail_set.begin(), n_failed, run_id_array);
	*nfail = n_failed;
	}
    catch(PestIndexError ex)
	{
		cerr << ex.what() << endl;
		err = 1;
	}
	return err;
}

int rmic_get_failed_runs_n(RunManager *run_manager_ptr, int *run_id_array, int nfail)
{
    int err = 0;
    try
	{
		std::fill_n(run_id_array, nfail, -999);
        const std::set<int> &fail_set = run_manager_ptr->get_failed_run_ids();
        int n_failed_tmp = min(int(fail_set.size()), nfail);
        std::copy_n(fail_set.begin(), n_failed_tmp, run_id_array);
	}
    catch(PestIndexError ex)
	{
		cerr << ex.what() << endl;
		err = 1;
	}
	return err;
}



int rmic_get_nruns(RunManager *run_manager_ptr, int *nruns)
{
    int err = 0;
    try {
        *nruns = run_manager_ptr->get_nruns();
	}
    catch(PestIndexError ex)
	{
		cerr << ex.what() << endl;
		err = 1;
	}
	return err;
}

int rmic_get_total_runs(RunManager *run_manager_ptr, int *total_runs)
{
    int err = 0;
    try {
        *total_runs = run_manager_ptr->get_total_runs();
	}
    catch(PestIndexError ex)
	{
		cerr << ex.what() << endl;
		err = 1;
	}
	return err;
}


int rmic_delete(RunManager *run_manager_ptr)
{
	int err = 0;
	try {
		delete run_manager_ptr;
	}
		catch(...)
	{
		err = 1;
	}
	return err;
}

}

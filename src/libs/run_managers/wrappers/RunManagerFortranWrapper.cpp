#include "RunManagerFortranWrapper.h"
#include "utilities.h"
#include "RunManagerPanther.h"
#include "RunManagerSerial.h"
#include "RunManagerGenie.h"
#include "pest_error.h"

typedef class RunManagerAbstract RunManagerAbstract;

static RunManagerAbstract *_run_manager_ptr_ = nullptr;
static ofstream fout_run_manager_log_file;

using namespace pest_utils;

extern "C"
{

int rmif_create_serial_(char *f_comline, int  *comline_str_len, int *comline_array_len,
	char *f_tpl, int  *tpl_str_len, int *tpl_array_len,
	char *f_inp, int  *inp_str_len, int *inp_array_len,
	char *f_ins, int  *ins_str_len, int *ins_array_len,
	char *f_out, int  *out_str_len, int *out_array_len,
	char *f_storfile, int *storfile_len,
	char *f_rundir, int *rundir_len, int *n_max_fail)

{
	vector<string> comline_vec =  fortran_str_array_2_vec(f_comline, *comline_str_len, *comline_array_len);
	vector<string> tpl_vec =  fortran_str_array_2_vec(f_tpl, *tpl_str_len, *tpl_array_len);
	vector<string> inp_vec =  fortran_str_array_2_vec(f_inp, *inp_str_len, *inp_array_len);
	vector<string> ins_vec =  fortran_str_array_2_vec(f_ins, *ins_str_len, *ins_array_len);
	vector<string> out_vec =  fortran_str_array_2_vec(f_out, *out_str_len, *out_array_len);
	string storfile =  fortran_str_2_string(f_storfile, *storfile_len);
	string rundir =  fortran_str_2_string(f_rundir, *rundir_len);
	int err = 0;
	try {
		_run_manager_ptr_ = new RunManagerSerial(comline_vec, tpl_vec, inp_vec, ins_vec,
			out_vec, storfile, rundir, *n_max_fail);
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}


int rmif_create_panther_(
	char *f_storfile, int *storfile_len,
	char *f_port, int *f_port_len,
	char *f_info_filename, int *info_filename_len, int *n_max_fail,
	double *overdue_reched_fac, double *overdue_giveup_fac,
	double *overdue_giveup_minutes)
{
	int err = 0;
	try {
		string storfile =  fortran_str_2_string(f_storfile, *storfile_len);
		string port =  fortran_str_2_string(f_port, *f_port_len);
		string info_filename =  fortran_str_2_string(f_info_filename, *info_filename_len);
		fout_run_manager_log_file.open(info_filename);
		_run_manager_ptr_ = new RunManagerPanther(storfile,
			port, fout_run_manager_log_file, *n_max_fail, *overdue_reched_fac, *overdue_giveup_fac,
			*overdue_giveup_minutes);
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}

int rmif_create_genie_(char *f_comline, int  *comline_str_len, int *comline_array_len,
	char *f_tpl, int  *tpl_str_len, int *tpl_array_len,
	char *f_inp, int  *inp_str_len, int *inp_array_len,
	char *f_ins, int  *ins_str_len, int *ins_array_len,
	char *f_out, int  *out_str_len, int *out_array_len,
	char *f_storfile, int *storfile_len,
	char *f_host, int *f_host_len,
	char *f_genie_tag, int *genie_tag_len)
{
	int err = 0;
	try {
		vector<string> comline_vec =  fortran_str_array_2_vec(f_comline, *comline_str_len, *comline_array_len);
		vector<string> tpl_vec =  fortran_str_array_2_vec(f_tpl, *tpl_str_len, *tpl_array_len);
		vector<string> inp_vec =  fortran_str_array_2_vec(f_inp, *inp_str_len, *inp_array_len);
		vector<string> ins_vec =  fortran_str_array_2_vec(f_ins, *ins_str_len, *ins_array_len);
		vector<string> out_vec =  fortran_str_array_2_vec(f_out, *out_str_len, *out_array_len);
		string storfile =  fortran_str_2_string(f_storfile, *storfile_len);
		string host =  fortran_str_2_string(f_host, *f_host_len);
		string genie_tag =  fortran_str_2_string(f_genie_tag, *genie_tag_len);
		_run_manager_ptr_ = new RunManagerGenie(comline_vec, tpl_vec, inp_vec, ins_vec, out_vec, storfile, genie_tag);
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}


int rmif_add_run_(double *parameter_data, int *npar, int *id)
{
	int err = 0;
	try {
		vector<double> data(parameter_data, parameter_data+*npar);
		*id = _run_manager_ptr_->add_run(data);
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}

int rmif_add_run_with_info_(double *parameter_data, int *npar, int *id,
	char *f_info_txt, int  *info_txt_len, double *info_value)
{
	int err = 0;
	try {
		string info_txt = fortran_str_2_string(f_info_txt, *info_txt_len);
		vector<double> data(parameter_data, parameter_data + *npar);
		*id = _run_manager_ptr_->add_run(data, info_txt, *info_value);
	}
	catch (...)
	{
		err = 1;
	}
	return err;
}

int rmif_initialize_(char *f_pname, int  *pname_str_len, int *pname_array_len,
				 char *f_oname, int  *oname_str_len, int *oname_array_len)

{
	int err = 0;
	try {
		vector<string> pname_vec =  fortran_str_array_2_vec(f_pname, *pname_str_len, *pname_array_len);
		vector<string> oname_vec =  fortran_str_array_2_vec(f_oname, *oname_str_len, *oname_array_len);
		_run_manager_ptr_->initialize(pname_vec, oname_vec);
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}

int rmif_initialize_restart_(char *f_storfile, int *storfile_len)
{
	int err = 0;
	try {
		string storfile =  fortran_str_2_string(f_storfile, *storfile_len);
		_run_manager_ptr_->initialize_restart(storfile);
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}

int rmif_reinitialize_()
{
    int err = 0;
	try {
	    _run_manager_ptr_->reinitialize();
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}

int rmif_run_()
{
	int err = 0;
	try {
		_run_manager_ptr_->run();
	}
	catch(...)
	{
		err = 1;
	}
	return err;
}

int rmif_run_until_(int *condition, int *no_ops, double *time_sec, int *return_cond)
{
	int err = 0;
	RunManagerAbstract::RUN_UNTIL_COND enum_input_cond;
	RunManagerAbstract::RUN_UNTIL_COND enum_return_cond;
	enum_input_cond = static_cast<RunManagerAbstract::RUN_UNTIL_COND>(*condition);
	try {
		enum_return_cond = _run_manager_ptr_->run_until(enum_input_cond, *no_ops, *time_sec);
	}
	catch (...)
	{
		err = 1;
	}
	*return_cond = static_cast<int>(enum_return_cond);
	return err;
}
int rmif_get_run_(int *run_id, double *parameter_data, int *npar, double *obs_data, int *nobs)
{
	int err = 1;
	bool success = false;
	size_t n_par = *npar;
	size_t n_obs = *nobs;
	try {
		success = _run_manager_ptr_->get_run(*run_id, parameter_data, n_par, obs_data, n_obs);
		if (success) err = 0;
	}
	catch(PestIndexError ex) {
		cerr << ex.what() << endl;
		err = 1;
	}
	return err;
}


int rmif_get_run_with_info_(int *run_id, double *parameter_data, int *npar, double *obs_data, int *nobs,
	char *f_info_txt, int  *info_txt_len, double *info_value)
{
	int err = 1;
	bool success = false;
	size_t n_par = *npar;
	size_t n_obs = *nobs;
	try {
		string info_txt;
		success = _run_manager_ptr_->get_run(*run_id, parameter_data, n_par, obs_data, n_obs, info_txt, *info_value);
		string_to_fortran_char(info_txt, f_info_txt, *info_txt_len);
		if (success) err = 0;
	}
	catch (PestIndexError ex) {
		cerr << ex.what() << endl;
		err = 1;
	}
	return err;

}



int rmif_delete_()
{
	int err = 0;
	try {
		delete _run_manager_ptr_;
	}
		catch(...)
	{
		err = 1;
	}
	return err;
}

int rmif_get_num_failed_runs_(int *nfail)
{
    int err = 0;
	*nfail = -999;
	try
	{
        const std::set<int> &fail_set = _run_manager_ptr_->get_failed_run_ids();
        *nfail = fail_set.size();
	}
    catch(PestIndexError ex)
	{
		cerr << ex.what() << endl;
		err = 1;
	}
	return err;
}

int rmif_get_failed_run_ids_(int *run_id_array, int *len_run_id_array)
{
    int err = 0;
    try
	{
		std::fill_n(run_id_array, *len_run_id_array, -999);
        const std::set<int> &fail_set = _run_manager_ptr_->get_failed_run_ids();
        size_t n_failed_tmp = min(fail_set.size(), size_t(*len_run_id_array));
        std::copy_n(fail_set.begin(), n_failed_tmp, run_id_array);
	}
    catch(PestIndexError ex)
	{
		cerr << ex.what() << endl;
		err = 1;
	}
	return err;
}

int rmif_get_num_total_runs_(int *nruns)
{
    int err = 0;
    try {
        *nruns = _run_manager_ptr_->get_total_runs();
	}
    catch(PestIndexError ex)
	{
		cerr << ex.what() << endl;
		err = 1;
	}
	return err;
}

}

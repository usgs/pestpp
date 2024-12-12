#ifndef RUNMANAGER_C_WRAP_H_
#define RUNMANAGER_C_WRAP_H_
#include "RunManagerAbstract.h"
#include "config_os.h"

typedef struct RunManagerAbstract RunManager;
extern "C"
{
#ifdef OS_WIN
extern __declspec(dllexport)
#endif
RunManager* rmic_create_serial(char **comline, int comline_array_len,
	char **tpl, int tpl_array_len,
	char **inp, int inp_array_len,
	char **ins, int ins_array_len,
	char **out, int out_array_len,
	char *storfile,
	char *rundir,
	int n_max_fail);

#ifdef OS_WIN
extern __declspec(dllexport)
#endif
RunManager* rmic_create_panther(char *storfile,
	char *port,
	char *info_filename,
	int n_max_fail,
	double overdue_reched_fac, double overdue_giveup_fac,
	double overdue_giveup_minutes);

#ifdef OS_WIN
extern __declspec(dllexport)
#endif
int rmic_initialize(RunManager *run_manager_ptr,
	char **pname, int pname_array_len,
	char **oname, int oname_array_len);

#ifdef OS_WIN
extern __declspec(dllexport)
#endif
int rmic_reinitialize(RunManager *run_manager_ptr);

#ifdef OS_WIN
extern __declspec(dllexport)
#endif
int rmic_add_run(RunManager *run_manager_ptr, double *parameter_data, int npar, int *id);

#ifdef OS_WIN
extern __declspec(dllexport)
#endif
int rmic_run(RunManager *run_manager_ptr);


#ifdef OS_WIN
extern __declspec(dllexport)
#endif
int rmic_run_until(RunManager *run_manager_ptr, int condition, int n_nops , double sec, int * return_cond);


#ifdef OS_WIN
extern __declspec(dllexport)
#endif
int rmic_get_run(RunManager *run_manager_ptr, int run_id, double *parameter_data, int npar, double *obs_data, int nobs);


#ifdef OS_WIN
extern __declspec(dllexport)
#endif
int rmic_get_num_failed_runs(RunManager *run_manager_ptr, int *nfail);

//*************************************************************************************
//******************************** IMPORTANT ******************************************
//The calling program is responsible for freeing the memory associated with run_id_array
//after calling this function by involking delete[] run_id_array
//*************************************************************************************
#ifdef OS_WIN
extern __declspec(dllexport)
#endif
int rmic_get_failed_runs_alloc(RunManager *run_manager_ptr, int *run_id_array, int *nfail);

#ifdef OS_WIN
extern __declspec(dllexport)
#endif
int rmic_get_failed_runs_n(RunManager *run_manager_ptr, int *run_id_array, int nfail);

#ifdef OS_WIN
extern __declspec(dllexport)
#endif
int rmic_get_nruns(RunManager *run_manager_ptr, int *nruns);

#ifdef OS_WIN
extern __declspec(dllexport)
#endif
int rmic_get_total_runs(RunManager *run_manager_ptr, int *total_runs);

#ifdef OS_WIN
extern __declspec(dllexport)
#endif
int rmic_delete(RunManager *run_manager_ptr);


}
#endif //RUNMANAGER_C_WRAP_H_

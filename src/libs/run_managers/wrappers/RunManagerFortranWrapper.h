#include "RunManagerAbstract.h"

#ifndef RUNMANAGER_FORTRAN_WRAP_H_
#define RUNMANAGER_FORTRAN_WRAP_H_

#ifdef EXTERNAL_UPPER
/* E.g. Windows Intel uses different name conventions */
#define rmif_create_serial_ RMIF_CREATE_SERIAL
#define rmif_create_panther_ RMIF_CREATE_PANTHER
#define rmif_add_run_ RMIF_ADD_RUN
#define rmif_add_run_with_info_ RMIF_ADD_RUN_WITH_INFO
#define rmif_initialize_ RMIF_INITIALIZE
#define rmif_initialize_restart_ RMIF_INITIALIZE_RESTART
#define rmif_reinitialize_ RMIF_REINITIALIZE
#define rmif_run_ RMIF_RUN
#define rmif_run_until_ RMIF_RUN_UNTIL
#define rmif_get_run_ RMIF_GET_RUN
#define rmif_get_run_with_info_ RMIF_GET_RUN_WITH_INFO
#define rmif_delete_ RMIF_DELETE
#define rmif_get_failed_run_ids_ RMIF_GET_FAILED_RUN_IDS
#define rmif_get_num_total_runs_ RMIF_GET_NUM_TOTAL_RUNS

#endif //EXTERNAL_UPPER

#endif //RUNMANAGER_FORTRAN_WRAP_H_

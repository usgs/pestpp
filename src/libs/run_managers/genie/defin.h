/*
 Copyright 2012, S.S. Papadopulos & Assoc. Inc.

    This file is part of GENIE.

    The GENIE is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GENIE is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GENIE.  If not, see<http://www.gnu.org/licenses/>.
*/

// cmuffels, Genie version 1, 2012

#ifndef DEFIN_H
#define DEFIN_H

// MISC -------------------------------------------------
#define NOT_SET 0
#define MAX_LINE_SIZE 300
#define FORTRAN_STRING_WIDTH 50
#define CXN_RETRY_MAX 120
#define CXN_RETRY_TIME 5
#define MIN_PORT 1024
#define MAX_PORT 65535

// MODEL IO ---------------------------------------------
#define TYPE_INSTRUCTION 0
#define TYPE_TEMPLATE 1

// SOCKETS ----------------------------------------------
#define MAX_QUEUE_SIZE 100
#define COMM_WAIT_TIME 10
#define CONNECTION_LOST 0
#define SELECT_TIMEOUT 0

// MESSAGE TYPES ----------------------------------------
#define COMMAND 0
#define STATUS_UPDATE 1
#define CONNECTION_TYPE 2
#define CONNECTION_NAME 3
#define CONNECTION_SPEED 4
#define RUN 5
#define OBS 6
#define SLAVE_COUNT 7

// SLAVE STATUS FLAGS -----------------------------------
#define WAITING_FOR_RUN 0
#define RUN_FAILED 2
#define READY_TO_START 3
#define RUN_COMPLETE 4
#define SEND_NSLAVE 5

// COMMANDS ---------------------------------------------
#define DO_NOTHING 1
#define COMMENCE_RUNS 2
#define COMMENCE_RUN 3
#define KILL_RUN 4
#define END_CONNECTION 17
#define KILL_ALL 6
#define KILL_SLAVE 7
#define TERMINATE_GMAN 8
#define RETURN_NSLAVE 9

// RUN/EXEC COMPONENTS
#define DIMENSIONS 1
#define DATA 2

// CONNECTION TYPES -------------------------------------
#define TYPE_NOT_SET 0
#define TYPE_FACE 1
#define TYPE_SLAVE 2
#define TYPE_MODEL 3

// CLIENT -----------------------------------------------
#define CLIENT_LOST 1
#define CLIENT_WAITING 2
#define CLIENT_WORKING 3

// THREAD -----------------------------------------------
#define THREAD_ERROR 0

// RUN PROCESS
#define PROCESS_COMPLETE 0

// MANAGED RUN STATUS
#define MRUN_WAITING 0
#define MRUN_INPROGRESS 1
#define MRUN_COMPLETE 2

// MESSAGES ---------------------------------------------
#define ERROR_MSG_TYPE -2
#define ERROR_RECV_DATA -1
#define ERROR_RECV_HEADER 0
#define MSG_HEADER_RECVD 1
#define MSG_NONE 2

// MISC ERRORS ------------------------------------------
#define NORMAL_TERMINATION 1
#define ERROR_TERMINATION 0
#define ERROR_GET_SPEED -100

// WINDOWS
#define ERROR_WSASTARTUP -1
#define ERROR_WSACLEANUP -2

// NODE
#define ERROR_NODE_INIT -3
#define ERROR_NODE_START -4
#define ERROR_NODE_SENDNAME -5
#define ERROR_NODE_SENDTYPE -6
#define ERROR_NODE_SENDSTATUS -7
#define ERROR_NODE_SENDCOMMAND -8
#define ERROR_NODE_SENDSPEED -9
#define ERROR_NODE_SENDOBS -10

// THREAD
#define ERROR_THREAD_START -30
#define ERROR_THREAD_CHKEXIT -31
#define ERROR_THREAD_EXEC -32

// HANDLE
#define ERROR_EVENT_INIT -40
#define ERROR_EVENT_WAIT -41
#define ERROR_EVENT_SET -42
#define ERROR_EVENT_RELEASE -43

// SOCKET
#define ERROR_SOCKET_WAIT -50
#define ERROR_SOCKET_ACCEPT -51
#define ERROR_SOCKET_PARSE -52
#define ERROR_SOCKET_PORT -53
#define ERROR_SOCKET_IP -54

// RUN
#define ERROR_RUN_SEND -60
#define ERROR_RUN_RECEIVE -61
#define ERROR_DIMENSIONS -62
#define ERROR_RUN_START -63
#define ERROR_RUN_WRITE_IN -64
#define ERROR_RUN_READ_OUT -65
#define ERROR_RUN_KILL -66

// TPL FILE
#define ERROR_TPL_EOF -70
#define ERROR_TPL_OPEN -71
#define ERROR_TPL_INOPEN -72
#define ERROR_TPL_MARKER -73
#define ERROR_TPL_READ -74
#define ERROR_TPL_MAP -75

// RUN MANAGER
#define ERROR_RM_START -80

// RESTART FILES
#define ERROR_RESTART_OPEN -90
#define ERROR_RESTART_WRITE -91

#endif /* DEFIN_H */

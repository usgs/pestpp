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

#ifndef THREAD_H
#define THREAD_H

#include <process.h>
#define _WINSOCKAPI_
#include <windows.h>

#include "defin.h"

//****************************************************************************
class THREAD
//****************************************************************************
//CTM Dec 2010

/*
*/

{

  public:

      HANDLE handle;
      HANDLE h_thEnd;
      unsigned int threadid;
      bool *iset;
      //PROCESS_INFORMATION pinfo;

      int initialize();
      void start();
      void set_primary(HANDLE &_primary);
      void set_event_comm(HANDLE &_event);
      void set_mutex_comm(HANDLE &_event);
      HANDLE get_parent_end();

  protected:

      // handles belonging to this thread
      HANDLE h_thComm;
      HANDLE m_thComm;

      // parent version of handles
      HANDLE primary;
      HANDLE ph_thEnd;
      HANDLE ph_thComm;
      HANDLE pm_thComm;

      //STARTUPINFO sinfo;

      virtual int execute()
                  {return (0);}

//    ------------------------------------------------------------------------
      static unsigned _stdcall ThreadEntryPoint(void *pThis)
//    ------------------------------------------------------------------------
      /*
            This routine is the entry point for the thread. The windows
            "create" thread method requires such a routine for execution
            as a class.
      */
      {
        THREAD*pthX=(THREAD*)pThis;
        return pthX->execute();
      }

};

#endif //THREAD_H

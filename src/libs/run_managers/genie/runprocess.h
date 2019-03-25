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

#ifndef RUNPROCESS_H
#define RUNPROCESS_H

#include <process.h>
#define _WINSOCKAPI_
#include <windows.h>
#include <tlhelp32.h>
#include <time.h>

#include "defin.h"
#include "executable.h"

//****************************************************************************
class RUNPROCESS
//****************************************************************************
//CTM Jan 2011

/*
      run process class
*/

{

  private:

      STARTUPINFO sinfo;
      EXECUTABLE *exec;


  public:

      RUNPROCESS();
      ~RUNPROCESS();

      bool *iset,*console_visible;

      PROCESS_INFORMATION pinfo;

      void initialize(bool *console);
      void set_executable(EXECUTABLE *_executable);
      int start();
//      bool __fastcall KillProcessTree(DWORD myprocID,DWORD dwTimeout);  // routine from Stackoverflow
      int __fastcall KillProcessTree(DWORD myprocID);  // routine from Stackoverflow

};

#endif //RUNPROCESS_H

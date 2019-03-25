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

#include <iostream>

#include "thread.h"

using namespace std;

//****************************************************************************
int THREAD::initialize()
//****************************************************************************
//CTM Dec 2010

/*
      Initialize thread necessities
*/

{

  // Create all access security rights for thread
  SECURITY_DESCRIPTOR sd={0};
  SECURITY_ATTRIBUTES saAttr;

  InitializeSecurityDescriptor(&sd,SECURITY_DESCRIPTOR_REVISION);
  SetSecurityDescriptorDacl(&sd,true,NULL,false);
  saAttr.nLength=sizeof(SECURITY_ATTRIBUTES);
  saAttr.lpSecurityDescriptor=&sd;
  saAttr.bInheritHandle=true;

  try
  {

    // get current process handle
    primary=GetCurrentProcess();

    // create kill event for thread
    ph_thEnd=CreateEvent(&saAttr,true,false,NULL);
    if(ph_thEnd==NULL) throw ERROR_EVENT_INIT;

    // create thread
    handle=(HANDLE)_beginthreadex(&saAttr,
                                  0,
                                  THREAD::ThreadEntryPoint,
                                  this,
                                  CREATE_SUSPENDED,
                                  &threadid
                                  );
    if(handle==0) throw -1;
    return 1;
  }
  catch (...)
  {
    cout<<"\n\n Error initializing thread."<<'\n';
    return 0;
  }

}

//****************************************************************************
HANDLE THREAD::get_parent_end()
//****************************************************************************
//CTM Mar 2011

/*
      Return primary thread's end event handle
*/

{
  return ph_thEnd;
}

//****************************************************************************
void THREAD::set_event_comm(HANDLE &_event)
//****************************************************************************
//CTM Mar 2011

/*
      Set event by which this thread can signal primary
*/

{
  ph_thComm=_event;
}

//****************************************************************************
void THREAD::set_mutex_comm(HANDLE &_event)
//****************************************************************************
//CTM Mar 2011

/*
      Set mutex by which this thread can share memory with primary thread
*/

{
  pm_thComm=_event;
}

//****************************************************************************
void THREAD::set_primary(HANDLE &_primary)
//****************************************************************************
//CTM Mar 2011

/*
      Set handle of primary thread
*/

{
  primary=_primary;
}

//****************************************************************************
void THREAD::start()
//****************************************************************************
//CTM Dec 2010

/*
*/

{
  ResumeThread(handle);
}

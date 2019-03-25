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

#include "runprocess.h"

int childstartedafterparent(SYSTEMTIME &p_time,SYSTEMTIME &c_time);
//void wait(int seconds);

using namespace std;

//****************************************************************************
RUNPROCESS::RUNPROCESS()
//****************************************************************************
//CTM Sep 2011

/*
      Constructor
*/

{
  iset=new bool;
  *iset=false;
  console_visible=new bool;
  *console_visible=false;
}

//****************************************************************************
RUNPROCESS::~RUNPROCESS()
//****************************************************************************
//CTM Sep 2011

/*
      Destructor
*/

{
  delete iset;
  delete console_visible;
}

//****************************************************************************
void RUNPROCESS::initialize(bool *console)
//****************************************************************************
//CTM Jan 2011
/*
      Routine to initialize process properties
*/

{
  *iset=true;
  memset(&sinfo,0,sizeof(sinfo));
  memset(&pinfo,0,sizeof(pinfo));
  sinfo.cb=sizeof(sinfo);
  console_visible=console;

}

//****************************************************************************
void RUNPROCESS::set_executable(EXECUTABLE *_executable)
//****************************************************************************
//CTM Mar 2011
/*
      Routine to set the executable associated with run
*/

{

  exec=_executable;

}

//****************************************************************************
int RUNPROCESS::start()
//****************************************************************************
//CTM Jan 2011
/*
      Routine to initialize process properties
*/

{

  string cmdline;

  cmdline.assign(exec->name + " " + exec->args);

  // Create all access security rights for process
  SECURITY_DESCRIPTOR sd={0};
  SECURITY_ATTRIBUTES saAttr;

  InitializeSecurityDescriptor(&sd,SECURITY_DESCRIPTOR_REVISION);
  SetSecurityDescriptorDacl(&sd,true,NULL,false);
  saAttr.nLength=sizeof(SECURITY_ATTRIBUTES);
  saAttr.lpSecurityDescriptor=&sd;
  saAttr.bInheritHandle=true;

//  string sexec=exec.path+exec.name;
  //cout<<"Executable name  : "<<exec->name.c_str()<<'\n';
  //cout<<"Arguments        : "<<exec->args.c_str()<<'\n';
  //cout<<"Arguments        : "<<cmdline.c_str()<<'\n';

  try
  {
    if(*console_visible)
    {
      if(!CreateProcess(NULL,               // run executable
                        (char*)cmdline.c_str(),    // run command arguments
                        &saAttr,
                        &saAttr,
                        true,
                        CREATE_NEW_CONSOLE,
                        NULL,
                        NULL,
                        &sinfo,
                        &pinfo)) throw -1;
    }
    else
    {
      if(!CreateProcess(NULL,               // run executable
                        (char*)cmdline.c_str(),    // run command arguments
                        &saAttr,
                        &saAttr,
                        true,
                        CREATE_NO_WINDOW,
                        NULL,
                        NULL,
                        &sinfo,
                        &pinfo)) throw -1;
    }
    return 1;
  }
  catch (...)
  {
    cout<<"\n\n Error creating process to execute run."<<'\n';
    cout<<" Error Code: "<<GetLastError()<<'\n';
    return 0;
  }

}

//****************************************************************************
int __fastcall RUNPROCESS::KillProcessTree(DWORD myprocID)
//****************************************************************************
/*
    Modified routine from Stackoverflow to kill a windows process tree
*/
{

  bool bCont=true;
  HANDLE hSnap,hChildProc,hProc;
  PROCESSENTRY32 pe;
  FILETIME p_time,c_time,e_time,k_time,u_time;
  SYSTEMTIME p_stime,c_stime;

  hChildProc=NULL;

  memset(&pe, 0, sizeof(PROCESSENTRY32));
  pe.dwSize = sizeof(PROCESSENTRY32);

  //DWORD dwTimeout=COMM_WAIT_TIME*CLOCKS_PER_SEC;

  cout<<"Killing processes associated with PID: "<<myprocID<<'\n';

  try
  {

    hSnap=CreateToolhelp32Snapshot(TH32CS_SNAPPROCESS, 0);
    if(hSnap==INVALID_HANDLE_VALUE) throw -1;

    if(!Process32First(hSnap,&pe)) throw -2;


    // kill child processes
    while(bCont)
    {
      if(pe.th32ParentProcessID==myprocID)
      {
        // get creation time of parent process
        cout<<"open parent process"<<'\n';
        hProc=OpenProcess(PROCESS_QUERY_INFORMATION,FALSE,myprocID);
        if(hProc==NULL) throw -11;
        cout<<"get creation time of parent process"<<'\n';
        if(!GetProcessTimes(hProc,&p_time,&e_time,&k_time,&u_time)) throw -8;
        cout<<"KILL PID: "<<pe.th32ProcessID<<'\n';
        if(!KillProcessTree(pe.th32ProcessID)) throw -3;
        cout<<"get handle to child process"<<'\n';
        hChildProc=OpenProcess(PROCESS_QUERY_INFORMATION, FALSE, pe.th32ProcessID);
        if(hChildProc==NULL) throw -4;
        // get creation time of parent process
        if(!GetProcessTimes(hChildProc,&c_time,&e_time,&k_time,&u_time)) throw -9;
        //convert filetime to system time for both the parent and child
        if(!FileTimeToSystemTime(&p_time,&p_stime)) throw -10;
        if(!FileTimeToSystemTime(&c_time,&c_stime)) throw -10;
        //if creation time of child process is after parent process then kill it
        if(childstartedafterparent(p_stime,c_stime))
        {
          cout<<"force end of child."<<'\n';
          //wait(iwait);
          hChildProc=OpenProcess(PROCESS_TERMINATE, FALSE, pe.th32ProcessID);
          //if(!TerminateProcess(hChildProc,1)) throw -5;
          if(TerminateProcess(hChildProc,1))
            CloseHandle(hChildProc);
        }
      }
      bCont=Process32Next(hSnap,&pe);
    }

    // kill the main process
    cout<<"terminate the main process."<<'\n';
    hProc=OpenProcess(PROCESS_TERMINATE, FALSE, myprocID);
    if(hChildProc==NULL) throw -6;
    cout<<"force termination of main process."<<'\n';
    if(!TerminateProcess(hProc,1)) throw -7;
    CloseHandle(hProc);

    return 1;

  }
  catch(int err)
  {
    switch(err)
    {
    case -1:
      cout<<"   *** Error taking 'snapshot' of Windows processes."<<'\n';
      cout<<"      GetLastError: "<<GetLastError()<<'\n';
      break;
    case -2:
      cout<<"   *** Windows process 'snapshot' is empty."<<'\n';
      break;
    case -3:
      cout<<"   *** Error cutting child process branch."<<'\n';
      break;
    case -4:
      cout<<"   *** Error opening handle to this process."<<'\n';
      break;
    case -5:
      cout<<"   *** Error forcing early termination of child process."<<'\n';
      cout<<"      GetLastError: "<<GetLastError()<<'\n';
      break;
    case -6:
      cout<<"   *** Error opening handle to parent process."<<'\n';
      break;
    case -7:
      cout<<"   *** Error forcing early termination of parent process."<<'\n';
      cout<<"      GetLastError: "<<GetLastError()<<'\n';
      break;
    case -8:
      cout<<"   *** Error getting creation time of parent process."<<'\n';
      cout<<"      GetLastError: "<<GetLastError()<<'\n';
      break;
    case -9:
      cout<<"   *** Error getting creation time of child process."<<'\n';
      cout<<"      GetLastError: "<<GetLastError()<<'\n';
      break;
    case -10:
      cout<<"   *** Error converting FILETIME to SYSTEMTIME."<<'\n';
      break;
    case -11:
      cout<<"   *** Error opening parent process to get creation time."<<'\n';
      break;
    }
    cout<<"\n   Error terminating process tree."<<'\n';
    return 0;
  }
}

//****************************************************************************
int childstartedafterparent(SYSTEMTIME &p_time,SYSTEMTIME &c_time)
//****************************************************************************
// CTM Apr 2011
/*
    Routine to determine if child process started after parent
*/
{
  if(c_time.wYear>p_time.wYear)
    return 1;
  else if(c_time.wYear<p_time.wYear)
    return 0;
  else
  {
    if(c_time.wMonth>p_time.wMonth)
      return 1;
    else if(c_time.wMonth<p_time.wMonth)
      return 0;
    else
    {
      if(c_time.wDay>p_time.wDay)
        return 1;
      else if(c_time.wDay<p_time.wDay)
        return 0;
      else
      {
        if(c_time.wHour>p_time.wHour)
          return 1;
        else if(c_time.wHour<p_time.wHour)
          return 0;
        else
        {
          if(c_time.wMinute>p_time.wMinute)
            return 1;
          else if(c_time.wMinute<p_time.wMinute)
            return 0;
          else
          {
            if(c_time.wSecond>p_time.wSecond)
              return 1;
            else if(c_time.wSecond<p_time.wSecond)
              return 0;
            else
            {
              if(c_time.wMilliseconds>=p_time.wMilliseconds)
                return 1;
              else
                return 0;
            }
          }
        }
      }
    }
  }

}

////******************************************************************************************************
//static void wait(int seconds)
////******************************************************************************************************
////from http://www.cplusplus.com/reference/clibrary/ctime/clock/
//{
//  clock_t endwait;
//  endwait = clock () + seconds * CLOCKS_PER_SEC ;
//  while (clock() < endwait) {}
//}

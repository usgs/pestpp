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


// The extern C block should be uncommented for compilation with Fortran
extern "C"
{
  int GENIE_KILL_GMAN(char*,char*);

  int GENIE_GET_NSLAVE(char*,char*);

  int GENIE_INTERFACE(int*,
                      int*,
                      char*,
                      int*,
                      int*,
                      char*,
                      char*,
                      double*,
                      double*,
                      int*,
                      int*,
                      char*,
                      char*,
                      char*,
                      char*,
                      char*,
                      char*,
                      int*);
}

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <list>
#include <vector>

#include <winsock2.h>
#include <ws2tcpip.h>
#define _WINSOCKAPI_
#include <windows.h>

#include "defin.h"
#include "modelrun.h"
#include "message.h"
#include "receivor.h"
#include "client.h"
#include "modelresult.h"
#include "genie_misc_routines.cpp"

using namespace std;

//****************************************************************************
int GENIE_KILL_GMAN(char *id,char *host)
//****************************************************************************

  // C.Muffels, S.S. Papadopulos & Assoc. Inc.
  // W.Shreuder, Principia Mathematica

  // Routine to connect to GMAN

  // Version 1.0.0.0  Apr 2012

{

  //--- COMMANDLINE ARGS FOR CLIENT INITIALIZATION
  int argc=7;
  char *argv[7]={"XXXXXXIDXXXXXX",
                 "/name"," XXXXXXXXIDXXXXXXX",
                 "/type"," model",
                 "/host","111.111.111.111:11111"};
  //--- Socket
  WSADATA wsaData;
  //--- Client
  CLIENT cxn;

  // INITIALIZE WINSOCK
  if(WSAStartup(WINSOCK_VERSION,&wsaData)!=0)
  {
    cout<<"\n\n    Winsock function 'WSAStartup' failed with error: "<<WSAGetLastError()<<'\n';
    return ERROR_WSASTARTUP;
  }


  // UPDATE "COMMANDLINE"
  argv[0]=id;
  argv[2]=id;
  argv[6]=host;

  // INITIALIZE CLIENT AND CONNECT TO RUN MANAGER
  cout<<"\n   Connecting "<<id<<" to GENIE RUN MANAGER."<<'\n';
  if(!cxn.initialize(argc,argv))
  {
    cout<<"\n\n   Error initializing interface as HOST."<<'\n';
    return ERROR_NODE_INIT;
  }
  if(!cxn.start())
  {
    cout<<"   *** Error starting interface CLIENT."<<'\n';
    return ERROR_NODE_START;
  }

  // INTRODUCE SELF TO RUN MANAGER
  cout<<"   Sending introductory messages to GENIE RUN MANAGER."<<'\n';
  if(!cxn.sendname())
  {
    cout<<"   *** Error sending HOST ID to GENIE RUN MANAGER."<<'\n';
    return ERROR_NODE_SENDNAME;
  }
  if(!cxn.sendtype())
  {
    cout<<"   *** Error sending HOST type (model) to GENIE RUN MANAGER."<<'\n';
    return ERROR_NODE_SENDTYPE;
  }

  // SEND TERMINATE COMMAND
  if(!cxn.sendcommand(TERMINATE_GMAN))
  {
    cout<<"   *** Error sending TERMINATE COMMAND to GENIE RUN MANAGER."<<'\n';
    return ERROR_NODE_SENDCOMMAND;
  }

  // DISCONNECT SLAVE
  cxn.stop();
  if(WSACleanup()==SOCKET_ERROR)
  {
    cout<<"\n\n   Winsock function 'WSACleanup' failed with error: "<<WSAGetLastError()<<'\n';
    return ERROR_WSACLEANUP;
  }

  // NORMAL TERMINATION
  cout<<"\n\n   NORMAL TERMINATION OF GENIE."<<'\n';
  return NORMAL_TERMINATION;

}

//****************************************************************************
int GENIE_GET_NSLAVE(char *id,char *host)
//****************************************************************************

  // C.Muffels, S.S. Papadopulos & Assoc. Inc.
  // W.Shreuder, Principia Mathematica

  // Routine to connect to GMAN and get the number of connected slaves

  // Version 1.0.0.0  July 2011

{

  //--- COMMANDLINE ARGS FOR CLIENT INITIALIZATION
  int argc=7;
  char *argv[7]={"XXXXXXIDXXXXXX",
                 "/name"," XXXXXXXXIDXXXXXXX",
                 "/type"," model",
                 "/host","111.111.111.111:11111"};
  //--- Socket
  WSADATA wsaData;
  //--- Client
  CLIENT cxn;
  int msg_wait_header;
  int nslave;
  MESSAGE *message;

  // INITIALIZE WINSOCK
  if(WSAStartup(WINSOCK_VERSION,&wsaData)!=0)
  {
    cout<<"\n\n    Winsock function 'WSAStartup' failed with error: "<<WSAGetLastError()<<'\n';
    return ERROR_WSASTARTUP;
  }


  // UPDATE "COMMANDLINE"
  argv[0]=id;
  argv[2]=id;
  argv[6]=host;

  // INITIALIZE CLIENT AND CONNECT TO RUN MANAGER
  cout<<"\n   Connecting "<<id<<" to GENIE RUN MANAGER."<<'\n';
  if(!cxn.initialize(argc,argv))
  {
    cout<<"\n\n   Error initializing interface as HOST."<<'\n';
    return ERROR_NODE_INIT;
  }
  if(!cxn.start())
  {
    cout<<"   *** Error starting interface CLIENT."<<'\n';
    return ERROR_NODE_START;
  }

  // INTRODUCE SELF TO RUN MANAGER
  cout<<"   Sending introductory messages to GENIE RUN MANAGER."<<'\n';
  if(!cxn.sendtype())
  {
    cout<<"   *** Error sending HOST type (model) to GENIE RUN MANAGER."<<'\n';
    return ERROR_NODE_SENDTYPE;
  }
  if(!cxn.sendname())
  {
    cout<<"   *** Error sending HOST ID to GENIE RUN MANAGER."<<'\n';
    return ERROR_NODE_SENDNAME;
  }

  // SEND COMMAND TO GET NUMBER OF CONNECTED SLAVES
  if(!cxn.sendcommand(RETURN_NSLAVE))
  {
    cout<<"   *** Error sending RETURN_NSLAVE to GENIE RUN MANAGER."<<'\n';
    return ERROR_NODE_SENDCOMMAND;
  }
  else
  {
    // successfully sent command - process response
    message=new MESSAGE;
    msg_wait_header=message->receiveheader(cxn.sock);
    switch(msg_wait_header)
    {
    case ERROR_RECV_HEADER:
      delete message;
      return ERROR_RECV_HEADER;
      break;
    case MSG_NONE:
      delete message;
      break;
    case MSG_HEADER_RECVD:
      cout<<"Have HEADER..."<<'\n';
      message->alloc();
      if(!message->receivedata(cxn.sock))
      {
        message->dealloc();
        delete message;
        return ERROR_RECV_DATA;
      }
      // process message
      cout<<"Have DATA..."<<'\n';
      switch(message->header.type)
      {
      case SLAVE_COUNT:
        message->data->copyto(0,(char*)&nslave,message->datasize);
        cout<<"NSLAVE: "<<nslave<<'\n';
        break;
      default:
        nslave=ERROR_MSG_TYPE;
        break;
      }
      message->dealloc();
      delete message;
      break;
    }
  }

  // DISCONNECT SLAVE
  cxn.stop();
  if(WSACleanup()==SOCKET_ERROR)
  {
    cout<<"\n\n   Winsock function 'WSACleanup' failed with error: "<<WSAGetLastError()<<'\n';
    return ERROR_WSACLEANUP;
  }

  // NORMAL TERMINATION
  cout<<"\n\n   NORMAL TERMINATION OF GENIE."<<'\n';
  return nslave;

}

//****************************************************************************
int GENIE_INTERFACE(int *nrun,
                    int *nexec,
                    char *_execnams,
                    int *npar,
                    int *nobs,
                    char *_apar,
                    char *_aobs,
                    double *pval,
                    double *oval,
                    int *ntpl,
                    int *nins,
                    char *_tplfle,
                    char *_infle,
                    char *_insfle,
                    char *_oufle,
                    char *host,
                    char *id,
                    int *ifail)
//****************************************************************************

  // C.Muffels, S.S. Papadopulos & Assoc. Inc.
  // W.Shreuder, Principia Mathematica

  // Interface between calling program and the Genie run manager

  // Version 1.0.0.0  Mar. 2011

  // IFAIL is an NRUN sized array with 0 for a successful model run and a 1 for failed

{

  void dealloc(int nrun,string *aobs,string *apar,string *execnams,
               string *insfle,string *infle,string *tplfle,
               string *oufle,char *ctmp,bool *runcomplete,
               MODEL_RUN *_run,list<MODEL_RUN*> *runs);

  //--- MISC
  //size_t nr;
  size_t n,ndone=0,itrk;
  bool *runcomplete,chkcomplete;
  ofstream outfile;

  list<MODEL_RESULT*> *results;
  list<MODEL_RESULT*>::iterator result;

  //--- External functions
  string trim(string &_string);
  string chartostring(char *from);

  //--- Genie Specific Variables
  string *aobs,*apar,*execnams;
  string *insfle,*infle,*tplfle,*oufle;

  //--- Char to string conversion
  char *ctmp;
  stringstream aobs_,apar_,execnams_;
  stringstream infle_,oufle_;
  stringstream tplfle_,insfle_;

  //--- Client
  CLIENT cxn;

  //--- Threads
  RECEIVOR *thread;
  HANDLE h_thComm;
  HANDLE m_thRuns;
  SECURITY_DESCRIPTOR sd={0};
  SECURITY_ATTRIBUTES saAttr;
  DWORD thExitCode;

  //--- RUN COLLECTION
  MODEL_RUN *_run;
  list<MODEL_RUN*> *runs;
  list<MODEL_RUN*>::iterator run;
  list<MODEL_RUN*>::iterator run_chk;

  results=new list<MODEL_RESULT*>;

  //--- COMMANDLINE ARGS FOR CLIENT INITIALIZATION
  int argc=7;
  char *argv[7]={"XXXXXXIDXXXXXX",
                 "/name"," XXXXXXXXIDXXXXXXX",
                 "/type"," model",
                 "/host","111.111.111.111:11111"};

  cout<<"\n GENIE Interface version 1.0.0.0"<<'\n';

  //--- Socket
  WSADATA wsaData;

  // INITIALIZE WINSOCK
  if(WSAStartup(WINSOCK_VERSION,&wsaData)!=0) {
    cout<<"\n\n    Winsock function 'WSAStartup' failed with error: "<<WSAGetLastError()<<'\n';
    return ERROR_WSASTARTUP;
  }

  // Initialize security requirements for events/mutexes
  InitializeSecurityDescriptor(&sd,SECURITY_DESCRIPTOR_REVISION);
  SetSecurityDescriptorDacl(&sd,true,NULL,false);
  saAttr.nLength=sizeof(SECURITY_ATTRIBUTES);
  saAttr.lpSecurityDescriptor=&sd;
  saAttr.bInheritHandle=true;

  // Open log file
  outfile.open("genie_interface.log",ios_base::app);
  outfile<<"\n\n\n New call for runs --->\n\n";

  // open synchronization objects
  outfile<<"    Creating synchronization objects...";
  h_thComm=CreateEvent(&saAttr,false,false,"ThreadComm");
  if(h_thComm==NULL)
  {
    cout<<"\n\n   Error creating thread COMMUNICATIONS event."<<'\n';
    return ERROR_EVENT_INIT;
  }
  //m_thRuns=new HANDLE;
  m_thRuns=CreateMutex(&saAttr,false,"RunsQueue");
  if(m_thRuns==NULL)
  {
    cout<<"\n\n   Error creating COMPLETE RUNS mutex."<<'\n';
    return ERROR_EVENT_INIT;
  }
  outfile<<" Complete.\n";

  // UPDATE "COMMANDLINE"
  argv[0]=id;
  argv[2]=id;
  argv[6]=host;

  // INITIALIZE CLIENT AND CONNECT TO RUN MANAGER
  outfile<<"    Connecting "<<id<<" to GENIE RUN MANAGER...";
  if(!cxn.initialize(argc,argv))
  {
    cout<<"\n\n   Error initializing interface as HOST."<<'\n';
    outfile<<"\n\n   Error initializing interface as HOST."<<'\n';
    return ERROR_NODE_INIT;
  }
  if(!cxn.start())
  {
    cout<<"   *** Error starting interface CLIENT."<<'\n';
    outfile<<"   *** Error starting interface CLIENT."<<'\n';
    return ERROR_NODE_START;
  }
  outfile<<" Complete.\n";
  outfile.flush();

  // INTRODUCE SELF TO RUN MANAGER
  outfile<<"    Sending introductory messages to GENIE RUN MANAGER..."<<'\n';
  outfile.flush();
  outfile<<"       NAME...";
  outfile.flush();
  if(!cxn.sendname())
  {
    cout<<"   *** Error sending HOST ID to GENIE RUN MANAGER."<<'\n';
    return ERROR_NODE_SENDNAME;
  }
  outfile<<" Complete."<<'\n';
  outfile.flush();
  outfile<<"       TYPE...";
  if(!cxn.sendtype())
  {
    cout<<"   *** Error sending HOST type (model) to GENIE RUN MANAGER."<<'\n';
    outfile<<"   *** Error sending HOST type (model) to GENIE RUN MANAGER."<<'\n';
    outfile.flush();
    return ERROR_NODE_SENDTYPE;
  }
  outfile<<" Complete."<<'\n';
  outfile.flush();

  // ASSIGN CONNECTION TO RECEIVOR THREAD TO HANDLE ANY FURTHER MESSAGES
  outfile<<"    Starting RECEIVOR thread to handle communication with GMAN...";
  outfile.flush();
  thread=new RECEIVOR;
  if(!thread->initialize())
  {
    cout<<"\n\n   Error starting RECEIVOR thread to manage communication with GENIE RUN MANAGER."<<'\n';
    outfile<<"\n\n   Error starting RECEIVOR thread to manage communication with GENIE RUN MANAGER."<<'\n';
    outfile.flush();
    return ERROR_THREAD_START;
  }
  thread->set_slave(&cxn);
  thread->set_event_comm(h_thComm);
  thread->set_mutex_comm(m_thRuns);
  thread->set_queue_results(results);
  thread->start();
  outfile<<" Complete"<<'\n';
  outfile.flush();

  // PUT LIST OF EXECUTABLES INTO ARRAY
  outfile<<"    Preparing executable(s)...";
  outfile.flush();
  ctmp=new char[MAX_LINE_SIZE];
  execnams=new string[*nexec];
  execnams_<<_execnams;
  for(n=0;n<*nexec;n++)
  {
    // parse char array according to appended "|"
    execnams_.getline(ctmp,MAX_LINE_SIZE,'|');
    execnams[n]=trim(chartostring(ctmp));
  }
  outfile<<"Complete."<<'\n';
  outfile.flush();

  outfile<<"    Preparing arrays..."<<'\n';
  // PUT LIST OF OBS NAMES INTO ARRAY
  outfile<<"       Observation names... ";
  outfile.flush();
  aobs=new string[*nobs];
  aobs_<<_aobs;
  for(n=0;n<*nobs;n++) {
    aobs_>>aobs[n];
  }
  outfile<<"Complete."<<'\n';
  outfile.flush();

  // PUT LIST OF PAR NAMES INTO ARRAY
  outfile<<"       Parameter names... ";
  outfile.flush();
  apar=new string[*npar];
  apar_<<_apar;
  for(n=0;n<*npar;n++)
    apar_>>apar[n];
  outfile<<"Complete."<<'\n';
  outfile.flush();

  // PUT LIST OF INPUT FILES INTO ARRAY
  outfile<<"       TPL files... ";
  outfile.flush();
  infle=new string[*ntpl];
  tplfle=new string[*ntpl];
  infle_<<_infle;
  tplfle_<<_tplfle;
  for(n=0;n<*ntpl;n++)
  {
    tplfle_>>tplfle[n];
    infle_>>infle[n];
  }
  outfile<<"Complete."<<'\n';
  outfile.flush();

  // PUT LIST OF OUTPUT FILES INTO ARRAY
  outfile<<"       INS files... ";
  outfile.flush();
  insfle=new string[*nins];
  oufle=new string[*nins];
  insfle_<<_insfle;
  oufle_<<_oufle;
  for(n=0;n<*nins;n++)
  {
    insfle_>>insfle[n];
    oufle_>>oufle[n];
  }
  outfile<<"Complete.\n";
  outfile.flush();

  // SETUP RUN COLLECTION
  outfile<<"    Preparing run collection...";
  outfile.flush();
  runs=new list<MODEL_RUN*>;
  if(*nrun==1)
    runcomplete=new bool;
  else
    runcomplete=new bool[*nrun];
  for(n=0;n<*nrun;n++)
  {
    _run=new MODEL_RUN();
    _run->npar=new size_t;
    *_run->npar=*npar;
    _run->nobs=new size_t;
    *_run->nobs=*nobs;
    _run->ntpl=new size_t;
    *_run->ntpl=*ntpl;
    _run->nins=new size_t;
    *_run->nins=*nins;
    _run->nexec=new size_t;
    *_run->nexec=*nexec;
    _run->id=new size_t;
    *_run->id=n+1;
    _run->set_exec(execnams);
    _run->parnams=apar;
    _run->obsnams=aobs;
    _run->parvals=&pval[n*(*npar)];
    _run->obsvals=&oval[n*(*nobs)];
    _run->tplfiles=tplfle;
    _run->infiles=infle;
    _run->insfiles=insfle;
    _run->outfiles=oufle;
    runs->push_back(_run);
    if(*nrun==1) {
      *runcomplete=false;
    }
    else {
      runcomplete[n]=false;
    }
  }
  outfile<<" Complete.\n";
/*
INCOMPLETE AT THIS TIME
  // write summary of first 100 runs to log file
  outfile<<"\n Run summary (first 100...)";
  outfile<<'\n';
  outfile<<setw(25)<<"Parnam";
  for(int i=0;i<*npar) {
    outfile<<setw(12)<<"Run"<<i+1;
    if(i==99) {
      break; }
  }
  for(int i=0;i<*npar) {
    for(run=runs.begin();run!=runs.end();++run) {
      i++;
      outfile<<setw(12)<<(*run)->parnams[i-1]<<i+1;
    }
    if(i==99) {
      break; }
  }
  outfile<<'\n';
  outfile.flush();
*/
  outfile<<"    Ready to execute "<<runs->size()<<" runs."<<'\n';
  outfile.flush();

  /*
  // SEND EACH RUN TO RUN MANAGER
  outfile<<"    Sending RUNS to GMAN...";
  for(run=runs->begin();run!=runs->end();++run)
  {
    if(!(*run)->send(cxn.sock))
    {
      cout<<"\n   *** Error sending run: "<<(*run)->id<<" to GENIE RUN MANAGER."<<'\n';
      outfile<<"\n   *** Error sending run: "<<(*run)->id<<" to GENIE RUN MANAGER."<<'\n';
      dealloc(*nrun,aobs,apar,execnams,insfle,infle,tplfle,oufle,ctmp,runcomplete,_run,runs);
      return ERROR_RUN_SEND;
    }
  }
  outfile<<" Complete."<<'\n';
  */

  // WAIT FOR MESSAGES FROM RUN MANAGER
  itrk=0;
  cout<<"\n    Waiting for runs to complete ---> (ID of completed runs follow)\n    ";

  //nr=(size_t)*nrun;
  outfile<<"    Sending runs to GMAN and waiting for runs to complete ---> (ID of completed runs follow)\n";
  run=runs->begin();
  while(1)
  {

    // SEND EACH RUN TO RUN MANAGER
    if(run!=runs->end())
    {
      if(!(*run)->send(cxn.sock))
      {
        cout<<"\n   *** Error sending run: "<<*(*run)->id<<" to GENIE RUN MANAGER."<<'\n';
        outfile<<"\n   *** Error sending run: "<<*(*run)->id<<" to GENIE RUN MANAGER."<<'\n';
        dealloc(*nrun,aobs,apar,execnams,insfle,infle,tplfle,oufle,ctmp,runcomplete,_run,runs);
        return ERROR_RUN_SEND;
      }
      ++run;
    }

    if(WaitForSingleObject(m_thRuns,0)==WAIT_OBJECT_0)
    {
      for(result=results->begin();result!=results->end();++result)
      {
        if(*(*result)->id<1||*(*result)->id>*nrun)
        {
          cout<<"   *** Error in completed RUN ID."<<'\n';
          (*result)->dealloc();
          delete (*result);
          (*result)=NULL;
          dealloc(*nrun,aobs,apar,execnams,insfle,infle,tplfle,oufle,ctmp,runcomplete,_run,runs);
          return -1;
        }
        // check if run is already reported as being complete
        if(*nrun==1)
          chkcomplete=*runcomplete;
        else
          chkcomplete=runcomplete[*(*result)->id-1];
        // update observation and parameter values in run arrays
        if(!chkcomplete)
        {
          // NEW RUN FAILURE  ----------------------------
          if(*npar>=999999) {
            // run failed
            ifail[*(*result)->id-1]=1;
          }
          else {
            // run succeeded
            ifail[*(*result)->id-1]=0;
            for(n=0;n<*npar;n++)
              pval[(*(*result)->id-1)*(*npar)+n]=(*result)->parvals[n];
          }
          // ---------------------------------------------
          for(n=0;n<*nobs;n++)
            oval[(*(*result)->id-1)*(*nobs)+n]=(*result)->obsvals[n];
          if(*nrun==1)
            *runcomplete=true;
          else
            runcomplete[*(*result)->id-1]=true;
          ndone++;
          itrk++;
          cout<<setw(6)<<*(*result)->id;
          outfile<<setw(6)<<*(*result)->id;
          outfile.flush();
          if(itrk==20)
          {
            outfile<<'\n';
            cout<<'\n';
            itrk=0;
          }
        }
        (*result)->dealloc();
        delete *result;
        *result=NULL;
      }
      results->clear();
      *thread->haveresult=false;
      if(!ReleaseMutex(m_thRuns)) {
        cout<<" Unable to release RUNS mutex."<<'\n';
        break;
      }
    }
    if(ndone==runs->size()) break;
  }
  cout<<"\n\n";
  outfile<<"\n    All runs are complete."<<'\n';
  outfile.flush();

  // LET THE RUN MANAGER KNOW WE ARE DONE WITH IT FOR NOW
  cxn.sendcommand(END_CONNECTION);

  // TERMINATE RECEIVOR THREAD
  outfile<<"    Kill RECEIVOR thread...";
  if(!SetEvent(thread->h_thEnd))
  {
    cout<<"\n\n   Unable to raise 'kill' flag to terminate RECEIVOR thread."<<'\n';
    outfile<<"\n\n   Unable to raise 'kill' flag to terminate RECEIVOR thread."<<'\n';
    dealloc(*nrun,aobs,apar,execnams,insfle,infle,tplfle,oufle,ctmp,runcomplete,_run,runs);
    return ERROR_EVENT_SET;
  }
  switch(WaitForSingleObject(thread->handle,INFINITE))
  {
  case WAIT_OBJECT_0:
    outfile<<" OK...";
    break;
  default:
    outfile<<" ERROR...";
    break;
  }
  if(!GetExitCodeThread(thread->handle,&thExitCode))
    cout<<" Error getting thread exit code"<<'\n';
  else
    cout<<thExitCode<<'\n';
  if(!CloseHandle(thread->handle))
    cout<<" Error closing thread handle."<<'\n';
  delete thread;
  outfile<<"    Complete.\n";

  // DEALLOCATE MEMORY
  outfile<<"    Dealloc memory...";
  dealloc(*nrun,aobs,apar,execnams,insfle,infle,tplfle,oufle,ctmp,runcomplete,_run,runs);
  outfile<<"    Complete.\n";

  // DISCONNECT SLAVE
  outfile<<"    Disconnect host...";
  cxn.stop();
  outfile<<"    Complete.\n";

  // CLOSE EVENT HANDLES
  CloseHandle(h_thComm);
  CloseHandle(m_thRuns);

  if(WSACleanup()==SOCKET_ERROR)
  {
    cout<<"\n\n   Winsock function 'WSACleanup' failed with error: "<<WSAGetLastError()<<'\n';
    outfile<<"\n\n   Winsock function 'WSACleanup' failed with error: "<<WSAGetLastError()<<'\n';
    return ERROR_WSACLEANUP;
  }

  // NORMAL TERMINATION
  outfile.flush();
  outfile<<"    NORMAL TERMINATION OF GENIE INTERFACE."<<'\n';
  outfile.close();
  cout<<"    NORMAL TERMINATION OF GENIE INTERFACE."<<'\n';
  return NORMAL_TERMINATION;

}

//*****************************************************************************
void dealloc(int nrun,
             string *aobs,
             string *apar,
             string *execnams,
             string *insfle,
             string *infle,
             string *tplfle,
             string *oufle,
             char *ctmp,
             bool *runcomplete,
             MODEL_RUN *_run,
             list<MODEL_RUN*> *runs)
//*****************************************************************************
//CTM, Apr 2011

/*
      Deallocate memory
*/

{

  list<MODEL_RUN*>::iterator run;

  // CLEAR LISTS
  //cout<<"Clear run lists"<<'\n';
  runs->clear();

  // DELETE ALLOCATED MEMORY
  //cout<<"Delete allocated memory"<<'\n';
  delete runs;
  delete _run;
  delete[] aobs;
  delete[] apar;
  delete[] tplfle;
  delete[] infle;
  delete[] insfle;
  delete[] oufle;
  delete ctmp;
  if(nrun==1)
    delete runcomplete;
  else
    delete[] runcomplete;

  // RESET POINTERS
  runs=NULL;
  _run=NULL;
  runcomplete=NULL;
  aobs=NULL;
  apar=NULL;
  execnams=NULL;
  tplfle=NULL;
  infle=NULL;
  insfle=NULL;
  oufle=NULL;
  ctmp=NULL;

}

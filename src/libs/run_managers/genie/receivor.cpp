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
#include <sstream>

#include "receivor.h"
#include "header.h"

using namespace std;

//****************************************************************************
RECEIVOR::RECEIVOR()
//****************************************************************************
//CTM Sep 2011

/*
      Constructor
*/

{
  haveresult=new bool;
  *haveresult=false;
  iset=new bool;
  *iset=false;
}

//****************************************************************************
RECEIVOR::~RECEIVOR()
//****************************************************************************
//CTM Sep 2011

/*
      Constructor
*/

{
  delete haveresult;
  delete iset;
  //CloseHandle(primary);
  //CloseHandle(ph_thEnd);
  //CloseHandle(ph_thComm);
  //CloseHandle(pm_thComm);
}

//****************************************************************************
void RECEIVOR::set_slave(CLIENT *_client)
//****************************************************************************
//CTM Dec 2010

/*
      Set slave associated with thread
*/

{
  slave=_client;
}

//****************************************************************************
MODEL_RESULT* RECEIVOR::get_result()
//****************************************************************************
//CTM Apr 2011

/*
      Get result associated with thread
*/

{
  return result;
}

//****************************************************************************
MODEL_RUN* RECEIVOR::get_run()
//****************************************************************************
//CTM Apr 2011

/*
      Get result associated with thread
*/

{
  return run;
}

//****************************************************************************
void RECEIVOR::clear_run()
//****************************************************************************
//CTM Apr 2011

/*
      Get result associated with thread
*/

{
  run=NULL;
}

//****************************************************************************
void RECEIVOR::set_run(MODEL_RUN *_run)
//****************************************************************************
//CTM Dec 2010

/*
      Set run associated with thread
*/

{
  run=_run;
}

//****************************************************************************
void RECEIVOR::set_queue_runs(list<MODEL_RUN*> *_runs)
//****************************************************************************
//CTM Dec 2010

/*
      Set new run shared with primary thread
*/

{
  runs=_runs;
}

//****************************************************************************
void RECEIVOR::set_queue_results(list<MODEL_RESULT*> *_results)
//****************************************************************************
//CTM Dec 2010

/*
      Set new run shared with primary thread
*/

{
  results=_results;
}

//****************************************************************************
int RECEIVOR::execute()
//****************************************************************************
//CTM Dec 2010

/*
      Receive and Process Messages
*/

{

  unsigned int event_wait_result;
  int msg_wait_header;

  endthread=false;
  *haveresult=false;
  *iset=true;

  try
  {
#ifdef VERBOSE
    cout<<" *** RECEIVOR object: routine 'execute'... "<<'\n';
#endif

    // duplicate handles of primary thread
#ifdef VERBOSE
    cout<<" x RECEIVOR: Duplicate handles etc. from primary thread."<<'\n';
#endif
    if(!DuplicateHandle(primary,ph_thComm,GetCurrentProcess(),&h_thComm,
                        0,false,DUPLICATE_SAME_ACCESS)) throw -2;
    if(!DuplicateHandle(primary,pm_thComm,GetCurrentProcess(),&m_thComm,
                        0,false,DUPLICATE_SAME_ACCESS)) throw -3;
    if(!DuplicateHandle(primary,ph_thEnd,GetCurrentProcess(),&h_thEnd,
                        0,false,DUPLICATE_SAME_ACCESS)) throw -4;

    while(1)
    {

      // check whether to terminate thread
#ifdef VERBOSE
    cout<<" x RECEIVOR: Check whether to terminate."<<'\n';
#endif
      if(endthread)
      {
#ifdef VERBOSE
    cout<<" x RECEIVOR: Yes."<<'\n';
#endif
        // Close handles
        //CloseHandle(h_thEnd);
        //CloseHandle(h_thComm);
        //CloseHandle(m_thComm);
        if(!DuplicateHandle(GetCurrentProcess(),h_thComm,GetCurrentProcess(),NULL,
                            0,false,DUPLICATE_CLOSE_SOURCE)) throw -2;
        if(!DuplicateHandle(GetCurrentProcess(),m_thComm,GetCurrentProcess(),NULL,
                            0,false,DUPLICATE_CLOSE_SOURCE)) throw -3;
        if(!DuplicateHandle(GetCurrentProcess(),h_thEnd,GetCurrentProcess(),NULL,
                            0,false,DUPLICATE_CLOSE_SOURCE)) throw -4;
        break;
      }

      // check whether terminate command received
      event_wait_result=WaitForSingleObject(h_thEnd,COMM_WAIT_TIME*CLOCKS_PER_SEC/100);
      if(event_wait_result==WAIT_OBJECT_0)
      {
        // Close handles
        //CloseHandle(h_thEnd);
        //CloseHandle(h_thComm);
        //CloseHandle(m_thComm);
        if(!DuplicateHandle(GetCurrentProcess(),h_thComm,GetCurrentProcess(),NULL,
                            0,false,DUPLICATE_CLOSE_SOURCE)) throw -2;
        if(!DuplicateHandle(GetCurrentProcess(),m_thComm,GetCurrentProcess(),NULL,
                            0,false,DUPLICATE_CLOSE_SOURCE)) throw -3;
        if(!DuplicateHandle(GetCurrentProcess(),h_thEnd,GetCurrentProcess(),NULL,
                            0,false,DUPLICATE_CLOSE_SOURCE)) throw -4;
        break;
      }

      // receive message
#ifdef VERBOSE
    cout<<" x RECEIVOR: Getting message..."<<'\n';
#endif
      message=new MESSAGE;
      msg_wait_header=message->receiveheader(slave->sock);
      switch(msg_wait_header)
      {
      case ERROR_RECV_HEADER:
#ifdef VERBOSE
    cout<<" x RECEIVOR: Error receiving message..."<<'\n';
#endif
        delete message;
#ifdef GSLAVE
        // if this is a slave machine try to reconnect to HOST
        if(!slave->start())
          throw ERROR_RECV_HEADER;
#else
        throw ERROR_RECV_HEADER;
#endif
        break;
      case MSG_NONE:
#ifdef VERBOSE
    cout<<" x RECEIVOR: No message to get..."<<'\n';
#endif
        delete message;
        break;
      case MSG_HEADER_RECVD:
#ifdef VERBOSE
        cout<<" x RECEIVOR: Getting message..."<<'\n';
#endif
        message->alloc();
        if(!message->receivedata(slave->sock))
        {
          message->dealloc();
          delete message;
#ifdef GSLAVE
          // if this is a slave machine try to reconnect to HOST
          if(!slave->start())
            throw ERROR_RECV_DATA;
#else
          throw ERROR_RECV_DATA;
#endif
        }
        // process message
#ifdef VERBOSE
        cout<<" x RECEIVOR: Processing message..."<<'\n';
#endif
        switch(message->header.type)
        {
        case CONNECTION_TYPE:
#ifdef VERBOSE
        cout<<" x RECEIVOR: Connection type."<<'\n';
#endif
          process_type();
          break;
        case CONNECTION_SPEED:
#ifdef VERBOSE
        cout<<" x RECEIVOR: Connection speed indicator."<<'\n';
#endif
          process_speed();
          break;
        case CONNECTION_NAME:
#ifdef VERBOSE
        cout<<" x RECEIVOR: Connection name."<<'\n';
#endif
          process_name();
          break;
        case STATUS_UPDATE:
#ifdef VERBOSE
        cout<<" x RECEIVOR: Status update."<<'\n';
#endif
          process_status();
          break;
        case COMMAND:
#ifdef VERBOSE
        cout<<" x RECEIVOR: Command."<<'\n';
#endif
          process_command();
          break;
        case RUN:
#ifdef VERBOSE
        cout<<" x RECEIVOR: Run."<<'\n';
#endif
          process_run();
          break;
        case OBS:
#ifdef VERBOSE
        cout<<" x RECEIVOR: Results."<<'\n';
#endif
          process_obs();
          break;
        }
        message->dealloc();
        delete message;
        break;
      }

    }
    return 1;
  }

  catch(int err)
  {
    slave->status=CLIENT_LOST;
    switch(err)
    {
    case -2:
      cout<<"\n\n   Error duplicating misc. communication flag."<<'\n';
      break;
    case -3:
      cout<<"\n\n   Error duplicating misc. shared access flag."<<'\n';
      break;
    case -4:
      cout<<"\n\n   Error duplicating termination flag."<<'\n';
      break;
    case ERROR_RECV_HEADER:
      cout<<" *** Communication lost with Slave: "<<slave->name<<'\n';
      break;
    case ERROR_RECV_DATA:
      cout<<" *** Communication lost with Slave: "<<slave->name<<'\n';
      break;
    case ERROR_EVENT_SET:
      cout<<"\n\n Error raising EVENT to signal PRIMARY thread."<<'\n';
      break;
    default:
      cout<<"Unexpected error receiving message...."<<'\n';
    }
    return 0;
  }
}

//****************************************************************************
void RECEIVOR::process_type()
//****************************************************************************
//CTM Mar 2011

/*
      Process TYPE message
*/

{
  message->data->copyto(0,(char*)&slave->type,message->datasize);
//  cout<<"received type."<<'\n';
}


//****************************************************************************
void RECEIVOR::process_speed()
//****************************************************************************
//CTM Mar 2011

/*
      Process SPEED or RATING message
*/

{
  message->data->copyto(0,(char*)&slave->rating,message->datasize);
}

//****************************************************************************
void RECEIVOR::process_name()
//****************************************************************************
//CTM Mar 2011

/*
      Process NAME message
*/

{
  string rv_chartostring(char *from);
  slave->name=rv_chartostring(message->data->buf);
}

//****************************************************************************
void RECEIVOR::process_command()
//****************************************************************************
//CTM Mar 2011

/*
      Process COMMAND message
*/

{

  int irec;

  message->data->copyto(0,(char*)&irec,message->datasize);

  switch(irec)
  {
  case COMMENCE_RUN:
    //cout<<" COMMAND: COMMENCE RUN."<<'\n';
    //start the run
    break;
  case KILL_RUN:
    //cout<<" COMMAND: KILL RUN."<<'\n';
    // terminate the process executing the run
    slave->status=KILL_RUN;
    break;
  case KILL_SLAVE:
    //cout<<" COMMAND: KILL SLAVE."<<'\n';
    // terminate the process executing the run
    slave->status=KILL_SLAVE;
    endthread=true;
    break;
  case KILL_ALL:
    //cout<<" COMMAND: KILL ALL."<<'\n';
    // terminate the process executing the run
    slave->status=KILL_ALL;
    break;
  case DO_NOTHING:
    break;
  case END_CONNECTION:
    //cout<<"message received to end connection."<<'\n';
    slave->stop();
    endthread=true;
    break;
  case RETURN_NSLAVE:
    slave->status=SEND_NSLAVE;
    break;
  case TERMINATE_GMAN:
    //cout<<"message received to terminate GMAN."<<'\n';
    endthread=true;
    slave->terminate=true;
    break;
  }

}

//****************************************************************************
void RECEIVOR::process_status()
//****************************************************************************
//CTM Mar 2011

/*
      Process STATUS message
*/

{

  int irec;

  message->data->copyto(0,(char*)&irec,message->datasize);

  switch(irec)
  {
  case WAITING_FOR_RUN:
    cout<<" STATUS: WAITING_FOR_RUN."<<'\n';
    break;
  case READY_TO_START:
    cout<<" STATUS: READY_TO_START."<<'\n';
    break;
  case RUN_FAILED:
    cout<<"\n\n Run Failed on Slave: "<<slave->name<<'\n';
    break;
  }

}

//****************************************************************************
int RECEIVOR::process_run()
//****************************************************************************
//CTM Mar 2011

/*
      Process RUN message
*/

{

  unsigned int mutex_wait_result;

  try
  {

    run=new MODEL_RUN;
    run->init();
    if(!run->set_frm_msg(message)) throw -1;
    mutex_wait_result=WaitForSingleObject(m_thComm,INFINITE);
    //cout<<"\n receivor has runs mutex."<<'\n';
    if(mutex_wait_result!=WAIT_OBJECT_0) throw ERROR_EVENT_WAIT;
    runs->push_back(run);
    if(!ReleaseMutex(m_thComm)) throw ERROR_EVENT_RELEASE;
    //cout<<"\n receivor released runs mutex."<<'\n';
    return 1;

  }

  catch(...)
  {
    cout<<"\n\n   Error receiving run!!"<<'\n';
    return 0;
  }

}

//****************************************************************************
int RECEIVOR::process_obs()
//****************************************************************************
//CTM Mar 2011

/*
      Process OBS message
*/

{

  unsigned int mutex_wait_result;

  try
  {

    result=new MODEL_RESULT;
    result->set_frm_msg(message);
    //cout<<"\n receivor is waiting for results mutex..."<<'\n';
    mutex_wait_result=WaitForSingleObject(m_thComm,INFINITE);
    if(mutex_wait_result!=WAIT_OBJECT_0) throw ERROR_EVENT_WAIT;
    //cout<<"\n receivor has results mutex."<<'\n';
    *haveresult=true;
    results->push_back(result);
    if(!ReleaseMutex(m_thComm)) throw ERROR_EVENT_RELEASE;
    //cout<<"\n receivor released results mutex."<<'\n';
    slave->status=RUN_COMPLETE;

    return 1;

  }

  catch(int ierr)
  {
    switch (ierr)
	{
	case ERROR_EVENT_WAIT:
		cout<<"\n\n Error waiting for array mutex while receiving run from slave."<<'\n';
		break;
	case ERROR_EVENT_RELEASE:
		cout<<"\n\n Error releasing array mutex while receiving run from slave."<<'\n';
		break;
	}
    return 0;
  }

}

//*****************************************************************************
string rv_chartostring(char *from)
//*****************************************************************************
//CTM, Nov 2010

/*
Routine to convert char to string
*/

{
  stringstream ss;

  ss<<from;
  return ss.str();
}

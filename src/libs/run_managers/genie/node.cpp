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
#include <set>

#include <time.h>

#include "node.h"
#include "header.h"
#include "defin.h"

using namespace std;

//string chartostring(char *from);

//****************************************************************************
void NODE::stop()
//****************************************************************************
//CTM Dec 2010

/*
      Routine to terminate client connection
*/

{

  // clear and delete timers
  run_start->clear();
  run_end->clear();
  delete run_start;
  delete run_end;

  shutdown(sock,2); //shutdown send and receives on socket
  closesocket(sock);

}

//****************************************************************************
void NODE::init_timers()
//****************************************************************************
//CTM Apr 2011

/*
      Routine to initialize run timers
*/

{

  run_start=new list<clock_t>;
  run_end=new list<clock_t>;

}

//****************************************************************************
void NODE::settype(int _type)
//****************************************************************************
//CTM Jan 2011

/*
      Routine to set node type
*/

{

  type=_type;

}

//****************************************************************************
int NODE::makenonblocking()
//****************************************************************************
//CTM Jan 2011

/*
      Routine to set node socket to non-blocking
*/

{

  int ichk;
  unsigned long imode=1;

  try
  {
    ichk=ioctlsocket(sock,FIONBIO, &imode);
    if(ichk==SOCKET_ERROR) throw;
    return 1;
  }
  catch(...)
  {
    cout<<"\n\n Failed to make SOCKET non-blocking."<<'\n';
    return 0;
  }

}

//****************************************************************************
void NODE::setprocess(RUNPROCESS *_process)
//****************************************************************************
//CTM Jan 2011

/*
      Routine to set client socket and add to file descriptor set
*/

{

  myprocess=_process;

}


//****************************************************************************
void NODE::endprocess()
//****************************************************************************
//CTM Jan 2011

/*
      Routine to delete process
*/

{

  delete myprocess;

}

//**********************************************************************
int NODE::sendcommand(int _command)
//**********************************************************************
//CTM, Jan 2011

/*
      Send command to socket
*/

{

  MESSAGE message;

  try
  {
    // initialize message
    //message.setmsgsize(sizeof(_command)+sizeof(message.header));
    message.setmsgsize(0+sizeof(message.header));
    //message.setdatasize(sizeof(_command));
    message.setdatasize(0);
    message.alloc();
    // initialize the message header
    //message.header.type=COMMAND;
    message.header.type=_command;
    message.header.bytesize=message.datasize;
    //message.header.nbytes=1;  // end of line
    message.header.nbytes=0;  // end of line
    message.header.compression=0;
    // copy type to message
    //message.data->copyfrom(0,(char*)&_command,message.datasize);
//    cout<<" MESSAGE: COMMAND: "<<(int)*message.data->buf<<'\n';
    // send message to host
    if(!message.sendme(sock)) throw -1;
    message.dealloc();
    return 1;
  }
  catch (...)
  {
    cout<<"\n\n Error sending COMMAND to slave."<<'\n';
    message.dealloc();
    return 0;
  }

}

//**********************************************************************
int NODE::sendstatus(int _status)
//**********************************************************************
//CTM, Jan 2011

/*
      Send status to socket
*/

{

  MESSAGE message;

  try
  {
    // initialize message
    message.setmsgsize(sizeof(_status)+sizeof(message.header));
    message.setdatasize(sizeof(_status));
    message.alloc();
    // initialize the message header
    message.header.type=STATUS_UPDATE;
    message.header.bytesize=message.datasize;
    message.header.nbytes=1;  // end of line
    message.header.compression=0;
    // copy type to message
    message.data->copyfrom(0,(char*)&_status,message.datasize);
//    cout<<" MESSAGE: STATUS: "<<(int)*message.data->buf<<'\n';
    // send message to host
    if(!message.sendme(sock)) throw -1;
    message.dealloc();
    return 1;
  }
  catch (...)
  {
    cout<<"\n\n Error sending STATUS to host."<<'\n';
    message.dealloc();
    return 0;
  }

}

//**********************************************************************
int NODE::sendtype()
//**********************************************************************
//CTM, Jan 2011

/*
      Send type to socket
*/

{

  MESSAGE message;

  try
  {
    // initialize message
    message.setmsgsize(sizeof(type)+sizeof(message.header));
    message.setdatasize(sizeof(type));
    message.alloc();
    // initialize the message header
    message.header.type=CONNECTION_TYPE;
    message.header.bytesize=message.datasize;
    message.header.nbytes=1;  // end of line
    message.header.compression=0;
    // copy type to message
    message.data->copyfrom(0,(char*)&type,message.datasize);
//    cout<<" MESSAGE: TYPE: "<<(int)*message.data->buf<<'\n';
    // send message to host
    if(!message.sendme(sock)) throw -1;
    message.dealloc();
    return 1;
  }
  catch (...)
  {
    cout<<"\n\n Error sending TYPE to host."<<'\n';
    message.dealloc();
    return 0;
  }

}

//**********************************************************************
int NODE::sendspeed()
//**********************************************************************
//CTM, Mar 2011

/*
      Send speed rating
*/

{

  MESSAGE message;

  try
  {
    // initialize message
    message.setmsgsize(sizeof(rating)+sizeof(message.header));
    message.setdatasize(sizeof(rating));
    message.alloc();
    // initialize the message header
    message.header.type=CONNECTION_SPEED;
    message.header.bytesize=message.datasize;
    message.header.nbytes=1;  // end of line
    message.header.compression=0;
    // copy type to message
    message.data->copyfrom(0,(char*)&rating,message.datasize);
    // send message to host
    if(!message.sendme(sock)) throw -1;
    message.dealloc();
    return 1;
  }
  catch (...)
  {
    cout<<"\n\n Error sending SPEED to host."<<'\n';
    message.dealloc();
    return 0;
  }

}

//**********************************************************************
int NODE::sendnslave(int nslave)
//**********************************************************************
//CTM, Mar 2011

/*
      Send nslave
*/

{

  MESSAGE message;

  try
  {
    // initialize message
    message.setmsgsize(sizeof(nslave)+sizeof(message.header));
    message.setdatasize(sizeof(nslave));
    message.alloc();
    // initialize the message header
    message.header.type=SLAVE_COUNT;
    message.header.bytesize=message.datasize;
    message.header.nbytes=1;  // end of line
    message.header.compression=0;
    // copy type to message
    message.data->copyfrom(0,(char*)&nslave,message.datasize);
    // send message to host
    if(!message.sendme(sock)) throw -1;
    message.dealloc();
    return 1;
  }
  catch (...)
  {
    cout<<"\n\n Error sending NSLAVE to host."<<'\n';
    message.dealloc();
    return 0;
  }

}

//**********************************************************************
int NODE::sendname()
//**********************************************************************
//CTM, Jan 2011

/*
      Send name to socket
*/

{

  MESSAGE message;

  try
  {
    // initialize message (+1 for end of line)
    message.setmsgsize(name.size()+sizeof(message.header)+1);
    message.setdatasize(name.size()+1);
    message.alloc();
    // initialize the message header
    message.header.type=CONNECTION_NAME;
    message.header.bytesize=1;
    message.header.nbytes=message.datasize;  // end of line
    message.header.compression=0;
    // copy type to message
    message.data->copyfrom(0,name.c_str(),message.datasize);
//    cout<<" MESSAGE: NAME: "<<message.data->buf<<'\n';
    // send message to host
    if(!message.sendme(sock)) throw -1;
    message.dealloc();
    return 1;
  }
  catch (...)
  {
    cout<<"\n\n Error sending NAME to host."<<'\n';
    message.dealloc();
    return 0;
  }

}
//
////**********************************************************************
//int NODE::sendrun()
////**********************************************************************
////CTM, Feb 2011
//
///*
//      Send RUN to socket
//*/
//
//{
//
//  try
//  {
////    if(!run->sendme(sock)) throw -1;
//    return 1;
//  }
//  catch (int err)
//  {
//    if(err==-1)
//      cout<<"\n\n   Error sending RUN to socket."<<'\n';
//    return 0;
//  }
//
//}

//****************************************************************************
double NODE::median_runtime()
//****************************************************************************
// CTM Apr 2011

/*
      Routine to calculate the median time (seconds) this thread is taking
      to send/execute/return a model run
*/
{

  //double t;
  unsigned int i;

  list<clock_t>::iterator is,ie;

  multiset<double> duration;
  multiset<double>::iterator dur;

  // if no runs are complete then return 1 over the rating
  if(run_end->size()==0)
    return 1./rating*CLOCKS_PER_SEC;

  // calculate time to complete runs
  is=run_start->begin();
  for(ie=run_end->begin();ie!=run_end->end();++ie)
  {
    duration.insert((*ie-*is)/CLOCKS_PER_SEC);
    ++is;
  }

  // if only one run completed it is the median
  if(duration.size()==1)
    return *duration.begin();

  // set is already sorted so just return median value
  dur=duration.begin();
  for(i=0;i<run_end->size()/2-1;i++)
  {
    duration.erase(dur);
    ++dur;
  }
  if(run_end->size()%2==0)
    return (*dur+*(dur++))/2.;
  else
    return *(dur++);

}

//****************************************************************************
double NODE::expected_timeleft()
//****************************************************************************
// CTM Apr 2011

/*
      Routine to estimate the amount of time (seconds) it will take this
      slave to complete its current run
*/
{

  double current_time,median_time;

  // if no runs are complete then return 1 over the rating
  if(run_end->size()==0)
    return 1./rating*CLOCKS_PER_SEC;

  // get the current time and the median time to complete a run
  median_time=median_runtime();
  if(median_time*rating/CLOCKS_PER_SEC==rating)
    current_time=0.;
  else
    current_time=(clock()-run_start->back())/CLOCKS_PER_SEC;
  return median_time-current_time;

}

////*****************************************************************************
//static string chartostring(char *from)
////*****************************************************************************
////CTM, Nov 2010
//
///*
//Routine to convert char to string
//*/
//
//{
//  stringstream ss;
//
//  ss<<from;
//  return ss.str();
//}

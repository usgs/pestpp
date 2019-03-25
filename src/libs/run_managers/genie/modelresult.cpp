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

#include "modelresult.h"
#include "message.h"

using namespace std;

//string chartostring(char *from);

//****************************************************************************
MODEL_RESULT::MODEL_RESULT()
//****************************************************************************
//CTM Oct 2011

/*
      Constructor
*/

{
  id=new size_t;
  npar=new size_t;
  nobs=new size_t;
  *id=0;
  *npar=0;
  *nobs=0;
}

//****************************************************************************
MODEL_RESULT::~MODEL_RESULT()
//****************************************************************************
//CTM Oct 2011

/*
      Constructor
*/

{
  delete npar;
  npar=NULL;
  delete nobs;
  nobs=NULL;
  delete id;
  id=NULL;
}

//****************************************************************************
int MODEL_RESULT::send(SOCKET &_socket)
//****************************************************************************
//CTM Feb 2011

/*
      Send components of model results to socket
*/

{

  size_t n;
  size_t datasize=0,copied=0;

  MESSAGE *message;

  try
  {

#ifdef VERBOSE
    cout<<" *** MODEL_RESULT object: routine 'send'... "<<'\n';
    cout<<" - Result "<<*id<<'\n';
    cout<<" - Sending RESULT to socket "<<_socket<<'\n';
#endif

    // set message size ------------------------------------------------------

    // size of run dimensions
#ifdef VERBOSE
    cout<<" - Set dimensions."<<'\n';
#endif
    datasize+=sizeof(*id);
    datasize+=sizeof(*npar)+
              sizeof(*nobs);

    // size of parameter values
#ifdef VERBOSE
    cout<<" - Prepare parameter values."<<'\n';
#endif
    for(n=0;n<*npar;n++)
      datasize+=sizeof(parvals[n]);

    // size of observation values
#ifdef VERBOSE
    cout<<" - Prepare observation values."<<'\n';
#endif
    for(n=0;n<*nobs;n++)
      datasize+=sizeof(obsvals[n]);

    // initialize message
#ifdef VERBOSE
    cout<<" - Initialize message."<<'\n';
#endif
    message=new MESSAGE;
    message->setmsgsize(sizeof(message->header)+datasize);
    message->setdatasize(datasize);
    message->alloc();

    // set header
#ifdef VERBOSE
    cout<<" - Set message header."<<'\n';
#endif
    message->header.type=OBS;
    message->header.bytesize=1;
    message->header.nbytes=datasize;
    message->header.compression=0;

    // set message buffer ----------------------------------------------------
#ifdef VERBOSE
    cout<<" - Set message buffer."<<'\n';
#endif

    // run dimensions
    message->data->copyfrom(copied,(char*)id,sizeof(*id));
    copied+=sizeof(*id);
    message->data->copyfrom(copied,(char*)npar,sizeof(*npar));
    copied+=sizeof(*npar);
    message->data->copyfrom(copied,(char*)nobs,sizeof(*nobs));
    copied+=sizeof(*nobs);

    // parameter values
    for(n=0;n<*npar;n++)
    {
      message->data->copyfrom(copied,(char*)&parvals[n],sizeof(parvals[n]));
      copied+=sizeof(parvals[n]);
//      cout<<parnams[n]<<" --> "<<parvals[n]<<'\n';
    }

    // observation values
    for(n=0;n<*nobs;n++)
    {
      message->data->copyfrom(copied,(char*)&obsvals[n],sizeof(obsvals[n]));
      copied+=sizeof(obsvals[n]);
    }

    // ensure dimensions match
    if(copied!=datasize)
    {
      message->dealloc();
      delete message;
      throw ERROR_DIMENSIONS;
    }

    // send message
#ifdef VERBOSE
    cout<<" - Send message."<<'\n';
#endif
    if(!message->sendme(_socket))
    {
      message->dealloc();
      delete message;
      throw ERROR_RUN_SEND;
    }

    // clear message
    message->dealloc();
    delete message;

    return 1;
  }

  catch(int err)
  {

    if(err==ERROR_RUN_SEND)
      cout<<"\n\n   Error sending RESULTS to client."<<'\n';
    else if(err==ERROR_DIMENSIONS)
    {
      cout<<"\n\n   Error creating RESULTS message."<<'\n';
      cout<<"   *** Data and message sizes are inconsistent."<<'\n';
    }

    return 0;
  }

}

//****************************************************************************
int MODEL_RESULT::set_frm_msg(MESSAGE *message)
//****************************************************************************
//CTM Feb 2011

/*
      Decompose message into necessary run parts
*/

{

  size_t n;
  size_t copied=0;

  try
  {

#ifdef VERBOSE
    cout<<" *** MODEL_RESULT object: routine 'set_frm_msg'... "<<'\n';
    cout<<" - Setting from message..."<<'\n';
#endif

    //init();
#ifdef VERBOSE
    cout<<" - Dimensions"<<'\n';
#endif
    message->data->copyto(copied,(char*)id,sizeof(*id));
    copied+=sizeof(*id);
#ifdef VERBOSE
    cout<<" + Run "<<*id<<'\n';
#endif
    message->data->copyto(copied,(char*)npar,sizeof(*npar));
    copied+=sizeof(*npar);
    message->data->copyto(copied,(char*)nobs,sizeof(*nobs));
    copied+=sizeof(*nobs);
    alloc();

    // Parameter Values
#ifdef VERBOSE
    cout<<" - Parameter values"<<'\n';
#endif
    for(n=0;n<*npar;n++)
    {
      message->data->copyto(copied,(char*)&parvals[n],sizeof(parvals[n]));
      copied+=sizeof(parvals[n]);
      //cout<<"PAR VALUE: "<<parvals[n]<<'\n';
    }

    // Observation Values
    for(n=0;n<*nobs;n++)
    {
      message->data->copyto(copied,(char*)&obsvals[n],sizeof(obsvals[n]));
      copied+=sizeof(obsvals[n]);
      //cout<<"OBS VALUE: "<<obsvals[n]<<'\n';
    }

    return 1;
  }

  catch(...)
  {
    cout<<"\n\n   Error setting RESULTS from message."<<'\n';
    return 0;
  }

}

//****************************************************************************
void MODEL_RESULT::set_frm_run(MODEL_RUN *run)
//****************************************************************************
//CTM Feb 2011

/*

*/

{

  *npar=*run->npar;
  *nobs=*run->nobs;
  *id=*run->id;
  alloc();
  parvals=run->parvals;
  obsvals=run->obsvals;

  //cout<<*npar<<", "<<*nobs<<", "<<*id<<'\n';

}

//****************************************************************************
void MODEL_RESULT::alloc() //size_t &_npar,size_t &_nobs,int &_id)
//****************************************************************************
//CTM Mar 2011

/*
      Initialize and allocate run conmponents
*/

{

  // allocate space for the different arrays
  parvals=new double[*npar];
  obsvals=new double[*nobs];

}

//****************************************************************************
void MODEL_RESULT::dealloc()
//****************************************************************************
//CTM Mar 2011

/*
      Deallocate run components
*/

{

  delete[] parvals;
  delete[] obsvals;

  //delete npar;
  //delete nobs;
  //delete id;

}

////****************************************************************************
//void MODEL_RESULT::init()
////****************************************************************************
////CTM Mar 2011
//
///*
//      Initialize and allocate run conmponents
//*/
//
//{
//
//  npar=new size_t;
//  nobs=new size_t;
//  id=new int;
//
//}

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

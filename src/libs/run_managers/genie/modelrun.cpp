/*s
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

#include <stdlib.h>

#include "modelrun.h"
#include "message.h"

using namespace std;

extern "C"
{

  void WRTTPL(int *,
              char *,
              char *,
              int *,
              char *,
              double *,
              int *);

  void READINS(int *,
               char *,
               char *,
               int *,
               char *,
               double *,
               int *);

  //void WRTSIG(int *,
  //            double *,
  //            char *,
  //            int *,
  //            int *,
  //            double *,
  //            int *);

}

//****************************************************************************
MODEL_RUN::MODEL_RUN()
//****************************************************************************
//CTM Oct 2011

/*
      Constructor
*/

{
  id=new size_t;
  *id=0;
}

#ifdef GSLAVE
//****************************************************************************
int MODEL_RUN::write_input()
//****************************************************************************
//CTM Mar 2011

/*
      Write model input files
*/

{

  size_t i,j;
  size_t nw=FORTRAN_STRING_WIDTH;

  int ii,jj,ifail;

  char *_tplfiles;
  char *_infiles;
  char *_parnams;

  //char space=' ';

  // initialize arrays
  _tplfiles=new char[*ntpl*nw];
  _infiles=new char[*ntpl*nw];
  _parnams=new char[*npar*nw];

  try
  {

#ifdef VERBOSE
    cout<<" *** MODEL_RUN object: routine 'write_input'... "<<'\n';
    cout<<" - Run "<<*id<<'\n';
#endif

    // Convert arrays of string to char for tpl and input files
    for(i=0;i<*ntpl;i++)
    {
      strcpy(&_tplfiles[i*nw],tplfiles[i].c_str());
      for(j=tplfiles[i].length();j<nw;j++)
        _tplfiles[i*nw+j]=' ';
      strcpy(&_infiles[i*nw],infiles[i].c_str());
      for(j=infiles[i].length();j<nw;j++)
        _infiles[i*nw+j]=' ';
    }

    // Convert string to arrays of chars for parameter name
    for(i=0;i<*npar;i++)
    {
      strcpy(&_parnams[i*nw],parnams[i].c_str());
      for(j=parnams[i].length();j<nw;j++)
        _parnams[i*nw+j]=' ';
    }

    // Call FORTRAN WRTTPL to write the files
    ii=(int)*ntpl;
    jj=(int)*npar;
    WRTTPL(&ii,_tplfiles,_infiles,&jj,_parnams,parvals,&ifail);
    if(ifail==0)
      throw NORMAL_TERMINATION;
    else
      throw -1;

  }

  catch(int err)
  {

    // Deallocate memory
    delete[] _tplfiles;
    delete[] _infiles;
    delete[] _parnams;

    switch(err)
    {
    case NORMAL_TERMINATION:
#ifdef VERBOSE
      cout<<" - Complete: No errors."<<'\n';
#endif
      return 1;
      break;
    default:
#ifdef VERBOSE
      cout<<" - Complete: With errors."<<'\n';
#endif
      cout<<" + Error in TPL file(s)"<<'\n';
      return 0;
      break;
    }

  }

}
//****************************************************************************
int MODEL_RUN::read_output()
//****************************************************************************
//CTM Apr 2011

/*
      Read model output files
*/

{

  size_t i,j;
  size_t nw=FORTRAN_STRING_WIDTH;

  int ii,jj,ifail;

  char *_insfiles;
  char *_outfiles;
  char *_obsnams;

  //char space=' ';

  // initialize arrays
  _insfiles=new char[*nins*nw];
  _outfiles=new char[*nins*nw];
  _obsnams=new char[*nobs*nw];

  try
  {

#ifdef VERBOSE
    cout<<" *** MODEL_RUN object: routine 'read_output'... "<<'\n';
    cout<<" - Run "<<*id<<'\n';
#endif

    // Convert arrays of string to char for ins and input files
#ifdef VERBOSE
    cout<<" - Processing instruction and output file names."<<'\n';
#endif
    for(i=0;i<*nins;i++)
    {
      //cout<<insfiles[i]<<'\n';
      strcpy(&_insfiles[i*nw],insfiles[i].c_str());
      for(j=insfiles[i].length();j<nw;j++)
        _insfiles[i*nw+j]=' ';
      strcpy(&_outfiles[i*nw],outfiles[i].c_str());
      for(j=outfiles[i].length();j<nw;j++)
        _outfiles[i*nw+j]=' ';
    }

    // Convert string to arrays of chars for parameter name
#ifdef VERBOSE
      cout<<" - Processing parameter names."<<'\n';
#endif
    for(i=0;i<*nobs;i++)
    {
      strcpy(&_obsnams[i*nw],obsnams[i].c_str());
      for(j=obsnams[i].length();j<nw;j++)
        _obsnams[i*nw+j]=' ';
    }

    // Call FORTRAN WRTINS to write the files
#ifdef VERBOSE
    cout<<" - Calling FORTRAN routine READINS..."<<'\n';
#endif
    //cout<<*nins<<'\n';
    //cout<<_insfiles<<'\n';
    //cout<<_outfiles<<'\n';
    //cout<<*nobs<<'\n';
    //cout<<_obsnams<<'\n';
    //cout<<obsvals<<'\n';
    ii=(int)*nins;
    jj=(int)*nobs;
    //cout<<ii<<", "<<jj<<'\n';
    READINS(&ii,_insfiles,_outfiles,&jj,_obsnams,obsvals,&ifail);
    if(ifail==0)
      throw NORMAL_TERMINATION;
    else
      throw -1;

  }

  catch(int err)
  {

    // Deallocate memory
    delete[] _insfiles;
    delete[] _outfiles;
    delete[] _obsnams;

    switch(err)
    {
    case NORMAL_TERMINATION:
#ifdef VERBOSE
      cout<<" - Complete: No errors."<<'\n';
#endif
      return 1;
      break;
    default:
#ifdef VERBOSE
      cout<<" - Complete: With errors."<<'\n';
#endif
      cout<<" + Error in INS file(s)"<<'\n';
      return 0;
      break;
    }

  }

}

//****************************************************************************
int MODEL_RUN::delete_output()
//****************************************************************************
//CTM May 2011

/*
      Delete any model output files from current folder
*/

{

  size_t n;
  string fnam;

#ifdef VERBOSE
    cout<<" *** MODEL_RUN object: routine 'delete_output'... "<<'\n';
    cout<<" - Run "<<*id<<'\n';
    cout<<" - Deleting RUN output files."<<'\n';
#endif

  for(n=0;n<*nins;n++)
  {
    fnam="del "+outfiles[n]+" >nul";
    system(fnam.c_str());
  }
#ifdef VERBOSE
    cout<<" - Complete."<<'\n';
#endif

  return 1;

}
#endif

//****************************************************************************
void MODEL_RUN::init()
//****************************************************************************
//CTM Mar 2011

/*
      Initialize and allocate run conmponents
*/

{

#ifdef VERBOSE
    cout<<" *** MODEL_RUN object: routine 'init'... "<<'\n';
    cout<<" - Run id not assigned yet."<<'\n';
    cout<<" - Initializing RUN."<<'\n';
#endif

  npar=new size_t;
  nobs=new size_t;
  ntpl=new size_t;
  nins=new size_t;
  nexec=new size_t;
  id=new size_t;

#ifdef VERBOSE
    cout<<" - Complete."<<'\n';
#endif

}

//****************************************************************************
void MODEL_RUN::alloc(size_t &_npar,size_t &_nobs,
                      size_t &_ntpl,size_t &_nins,
                      size_t &_nexec, size_t &_id)
//****************************************************************************
//CTM Mar 2011

/*
      Initialize and allocate run conmponents
*/

{

  size_t n;

  npar=&_npar;
  nobs=&_nobs;
  ntpl=&_ntpl;
  nins=&_nins;
  nexec=&_nexec;
  id=&_id;

#ifdef VERBOSE
    cout<<" *** MODEL_RUN object: routine 'alloc'... "<<'\n';
    cout<<" - Run "<<*id<<'\n';
    cout<<" - Allocate arrays for RUN."<<'\n';
#endif

  // allocate space for the different arrays

  exec=new EXECUTABLE[*nexec];
  for(n=0;n<*nexec;n++)
  {
    exec[n].name='\0';
    exec[n].args='\0';
    exec[n].path='\0';
  }

  parnams=new string[*npar];
  parvals=new double[*npar];

  obsnams=new string[*nobs];
  obsvals=new double[*nobs];

  tplfiles=new string[*ntpl];
  infiles=new string[*ntpl];
  insfiles=new string[*nins];
  outfiles=new string[*nins];

#ifdef VERBOSE
    cout<<" - Complete."<<'\n';
#endif

}

//****************************************************************************
void MODEL_RUN::dealloc()
//****************************************************************************
//CTM Mar 2011

/*
      Deallocate run components
*/

{

#ifdef VERBOSE
    cout<<" *** MODEL_RUN object: routine 'dealloc'... "<<'\n';
    cout<<" - Run "<<*id<<'\n';
    cout<<" - Deallocate arrays for RUN."<<'\n';
#endif

  delete[] exec;

  delete[] parnams;
  delete[] parvals;

  delete[] obsnams;
  delete[] obsvals;

  delete[] tplfiles;
  delete[] infiles;
  delete[] insfiles;
  delete[] outfiles;

  delete npar;
  delete nobs;
  delete ntpl;
  delete nins;
  delete nexec;
  delete id;

#ifdef VERBOSE
    cout<<" - Complete."<<'\n';
#endif

}

//****************************************************************************
int MODEL_RUN::set_exec(string *exec_ary)
//****************************************************************************
//CTM Mar 2011

/*
      Set the name and arguments for each executable
*/

{

  size_t n;
  string stmp;
  stringstream _execnams,_args;

  exec=new EXECUTABLE[*nexec];

  try
  {

#ifdef VERBOSE
    cout<<" *** MODEL_RUN object: routine 'set_exec'... "<<'\n';
    cout<<" - Run "<<*id<<'\n';
    cout<<" - Setting name and arguments for each executable."<<'\n';
#endif

    for(n=0;n<*nexec;n++)
    {
      exec[n].name='\0';
      exec[n].args='\0';
      exec[n].path='\0';
      _execnams<<exec_ary[n].c_str();
      if(_execnams.eof()) throw -1;
      _execnams>>exec[n].name;
      while(!_execnams.eof())
      {
        _execnams>>stmp;
        _args<<stmp;
      }
      exec[n].args=_args.str();
    }
    return 1;
  }

  catch(int err)
  {
    if(err==-1) {
      cout<<" - Complete: with error."<<'\n';
      cout<<" * Empty model executable."<<'\n';
    }
#ifdef VERBOSE
    else {
      cout<<" - Complete: No error."<<'\n';
    }
#endif
    return 0;
  }

}

//****************************************************************************
int MODEL_RUN::send(SOCKET &_socket)
//****************************************************************************
//CTM Feb 2011

/*
      Send components of model run to socket
*/

{

  size_t n;
  size_t datasize=0,copied=0;
  size_t itmp;

  MESSAGE *message;

  try
  {

#ifdef VERBOSE
    cout<<" *** MODEL_RUN object: routine 'send'... "<<'\n';
    cout<<" - Run "<<*id<<'\n';
    cout<<" - Sending RUN to socket "<<_socket<<'\n';
#endif

    // set message size ------------------------------------------------------

    // size of run dimensions
#ifdef VERBOSE
    cout<<" - Set dimensions."<<'\n';
#endif
    datasize+=sizeof(*id);
    datasize+=sizeof(*nexec);
    datasize+=sizeof(*npar)+
              sizeof(*nobs)+
              sizeof(*ntpl)+
              sizeof(*nins);

    // size of executable parts
#ifdef VERBOSE
    cout<<" - Prepare executable info."<<'\n';
#endif
    for(n=0;n<*nexec;n++)
    {
      datasize+=sizeof(exec[n].name.size());
      if(exec[n].name.size()>1) datasize+=exec[n].name.size()+1;
      datasize+=sizeof(exec[n].args.size());
      if(exec[n].args.size()>1) datasize+=exec[n].args.size()+1;
      datasize+=sizeof(exec[n].path.size());
      if(exec[n].path.size()>1) datasize+=exec[n].path.size()+1;
    }

    // size of parameter names
#ifdef VERBOSE
    cout<<" - Prepare parameter names."<<'\n';
#endif
    for(n=0;n<*npar;n++)
    {
      datasize+=sizeof(parnams[n].size());
      datasize+=parnams[n].size()+1;
    }

    // size of parameter values
#ifdef VERBOSE
    cout<<" - Prepare parameter values."<<'\n';
#endif
    for(n=0;n<*npar;n++) {
      datasize+=sizeof(parvals[n]);
    }

    // size of observation names
#ifdef VERBOSE
    cout<<" - Prepare observation names."<<'\n';
#endif
    for(n=0;n<*nobs;n++)
    {
      datasize+=sizeof(obsnams[n].size());
      datasize+=obsnams[n].size()+1;
    }

    // size of observation values
#ifdef VERBOSE
    cout<<" - Prepare observation values."<<'\n';
#endif
    for(n=0;n<*nobs;n++) {
      datasize+=sizeof(obsvals[n]);
    }

    // size of tpl file names
#ifdef VERBOSE
    cout<<" - Prepare template file info."<<'\n';
#endif
    for(n=0;n<*ntpl;n++)
    {
      datasize+=sizeof(tplfiles[n].size());
      datasize+=tplfiles[n].size()+1;
      datasize+=sizeof(infiles[n].size());
      datasize+=infiles[n].size()+1;
    }

    // size of ins file names
#ifdef VERBOSE
    cout<<" - Prepare instruction file info."<<'\n';
#endif
    for(n=0;n<*nins;n++)
    {
      datasize+=sizeof(insfiles[n].size());
      datasize+=insfiles[n].size()+1;
      datasize+=sizeof(outfiles[n].size());
      datasize+=outfiles[n].size()+1;
    }

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
    message->header.type=RUN;
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
    message->data->copyfrom(copied,(char*)nexec,sizeof(*nexec));
    copied+=sizeof(*nexec);
    message->data->copyfrom(copied,(char*)npar,sizeof(*npar));
    copied+=sizeof(*npar);
    message->data->copyfrom(copied,(char*)nobs,sizeof(*nobs));
    copied+=sizeof(*nobs);
    message->data->copyfrom(copied,(char*)ntpl,sizeof(*ntpl));
    copied+=sizeof(*ntpl);
    message->data->copyfrom(copied,(char*)nins,sizeof(*nins));
    copied+=sizeof(*nins);

    // executable properties
    for(n=0;n<*nexec;n++)
    {
      // name
      if(exec[n].name.size()<=1)
        itmp=0;
      else
        itmp=exec[n].name.size()+1;
      message->data->copyfrom(copied,(char*)&itmp,sizeof(exec[n].name.size()));
      copied+=sizeof(exec[n].name.size());
      if(itmp>0)
      {
        message->data->copyfrom(copied,exec[n].name.c_str(),exec[n].name.size()+1);
        copied+=exec[n].name.size()+1;
      }
      // arguments
      if(exec[n].args.size()<=1) {
        itmp=0; }
      else {
        itmp=exec[n].args.size()+1; }
      message->data->copyfrom(copied,(char*)&itmp,sizeof(exec[n].args.size()));
      copied+=sizeof(exec[n].args.size());
      if(itmp>0)
      {
//        cout<<"ARGS: "<<exec[n].args<<'\n';
//        cout<<"ARGS: "<<exec[n].args.c_str()<<'\n';
        message->data->copyfrom(copied,exec[n].args.c_str(),exec[n].args.size()+1);
        copied+=exec[n].args.size()+1;
      }
      // path
      if(exec[n].path.size()<=1) {
        itmp=0; }
      else {
        itmp=exec[n].path.size()+1; }
      message->data->copyfrom(copied,(char*)&itmp,sizeof(exec[n].path.size()));
      copied+=sizeof(exec[n].path.size());
      if(itmp>0)
      {
        message->data->copyfrom(copied,exec[n].path.c_str(),exec[n].path.size()+1);
        copied+=exec[n].path.size()+1;
      }
    }

    // parameter names
    for(n=0;n<*npar;n++)
    {
      itmp=parnams[n].size()+1;
      message->data->copyfrom(copied,(char*)&itmp,sizeof(parnams[n].size()));
      copied+=sizeof(parnams[n].size());
      message->data->copyfrom(copied,parnams[n].c_str(),parnams[n].size()+1);
      copied+=parnams[n].size()+1;
    }

    // parameter values
    for(n=0;n<*npar;n++)
    {
      message->data->copyfrom(copied,(char*)&parvals[n],sizeof(parvals[n]));
      copied+=sizeof(parvals[n]);
//      cout<<parnams[n]<<" --> "<<parvals[n]<<'\n';
    }

    // observation names
    for(n=0;n<*nobs;n++)
    {
      itmp=obsnams[n].size()+1;
      message->data->copyfrom(copied,(char*)&itmp,sizeof(obsnams[n].size()));
      copied+=sizeof(obsnams[n].size());
      message->data->copyfrom(copied,obsnams[n].c_str(),obsnams[n].size()+1);
      copied+=obsnams[n].size()+1;
    }

    // observation values
    for(n=0;n<*nobs;n++)
    {
      message->data->copyfrom(copied,(char*)&obsvals[n],sizeof(obsvals[n]));
      copied+=sizeof(obsvals[n]);
    }

    // template file names
    for(n=0;n<*ntpl;n++)
    {
      itmp=tplfiles[n].size()+1;
      message->data->copyfrom(copied,(char*)&itmp,sizeof(tplfiles[n].size()));
      copied+=sizeof(tplfiles[n].size());
      message->data->copyfrom(copied,tplfiles[n].c_str(),tplfiles[n].size()+1);
      copied+=tplfiles[n].size()+1;
    }

    // model input file names
    for(n=0;n<*ntpl;n++)
    {
      itmp=infiles[n].size()+1;
      message->data->copyfrom(copied,(char*)&itmp,sizeof(infiles[n].size()));
      copied+=sizeof(infiles[n].size());
      message->data->copyfrom(copied,infiles[n].c_str(),infiles[n].size()+1);
      copied+=infiles[n].size()+1;
    }

    // instruction file names
    for(n=0;n<*nins;n++)
    {
      itmp=insfiles[n].size()+1;
      message->data->copyfrom(copied,(char*)&itmp,sizeof(insfiles[n].size()));
      copied+=sizeof(insfiles[n].size());
      message->data->copyfrom(copied,insfiles[n].c_str(),insfiles[n].size()+1);
      copied+=insfiles[n].size()+1;
    }

    // model output file names
    for(n=0;n<*nins;n++)
    {
      itmp=outfiles[n].size()+1;
      message->data->copyfrom(copied,(char*)&itmp,sizeof(outfiles[n].size()));
      copied+=sizeof(outfiles[n].size());
      message->data->copyfrom(copied,outfiles[n].c_str(),outfiles[n].size()+1);
      copied+=outfiles[n].size()+1;
    }

    // ensure dimensions match
#ifdef VERBOSE
    cout<<" COPIED SIZE: "<<copied<<'\n';
    cout<<" DATA SIZE: "<<datasize<<'\n';
#endif
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

    //cout<<" Number of bytes in message: "<<sizeof(*message)<<'\n';
    //cout<<" Data size: "<<datasize<<'\n';
    //cout<<" Header size: "<<sizeof(message->header)<<'\n';
    //cout<<" Buffer size: "<<sizeof(*message->data)<<'\n';

    // clear message
    message->dealloc();
    delete message;

    return 1;
  }

  catch(int err)
  {

    if(err==ERROR_RUN_SEND)
      cout<<"\n\n   Error sending RUN to client."<<'\n';
    else if(err==ERROR_DIMENSIONS)
    {
      cout<<"\n\n   Error creating RUN message."<<'\n';
      cout<<"   *** Data and message sizes are inconsistent."<<'\n';
    }

    return 0;
  }

}

//****************************************************************************
int MODEL_RUN::save(fstream &_file)
//****************************************************************************
//CTM May 2011

/*
      Save model run to binary file _file
*/

{

  size_t n;
  size_t itmp;
  size_t one=1,zero=0;

  try
  {

    // run dimensions
    _file.write((char*)id,sizeof(*id));
    _file.write((char*)nexec,sizeof(*nexec));
    _file.write((char*)npar,sizeof(*npar));
    _file.write((char*)nobs,sizeof(*nobs));
    _file.write((char*)ntpl,sizeof(*ntpl));
    _file.write((char*)nins,sizeof(*nins));

    // executable properties
    for(n=0;n<*nexec;n++)
    {
      // name
      if(exec[n].name.size()<=1)
        itmp=0;
      else
        itmp=exec[n].name.size()+1;
      _file.write((char*)&itmp,sizeof(exec[n].name.size()));
      if(itmp>0)
        _file.write(exec[n].name.c_str(),(streamsize)(exec[n].name.size()+one));
      // arguments
      if(exec[n].args.size()<=1)
        itmp=0;
      else
        itmp=exec[n].args.size()+1;
      _file.write((char*)&itmp,sizeof(exec[n].args.size()));
      if(itmp>0)
        _file.write(exec[n].args.c_str(),(streamsize)(exec[n].args.size()+one));
      // path
      if(exec[n].path.size()<=one)
        itmp=zero;
      else
        itmp=exec[n].path.size()+1;
      _file.write((char*)&itmp,sizeof(exec[n].args.size()));
      if(itmp>zero)
        _file.write(exec[n].path.c_str(),(streamsize)(exec[n].path.size()+one));
    }

    // parameter names
    for(n=0;n<*npar;n++)
    {
      itmp=parnams[n].size()+1;
      _file.write((char*)&itmp,sizeof(parnams[n].size()));
      _file.write(parnams[n].c_str(),(streamsize)(parnams[n].size()+one));
    }

    // parameter values
    for(n=0;n<*npar;n++)
      _file.write((char*)&parvals[n],sizeof(parvals[n]));

    // observation names
    for(n=0;n<*nobs;n++)
    {
      itmp=obsnams[n].size()+1;
      _file.write((char*)&itmp,sizeof(obsnams[n].size()));
      _file.write(obsnams[n].c_str(),(streamsize)(obsnams[n].size()+one));
    }

    // observation values
    for(n=0;n<*nobs;n++)
      _file.write((char*)&obsvals[n],sizeof(obsvals[n]));

    // template file names
    for(n=0;n<*ntpl;n++)
    {
      itmp=tplfiles[n].size()+1;
      _file.write((char*)&itmp,sizeof(tplfiles[n].size()));
      _file.write(tplfiles[n].c_str(),(streamsize)(tplfiles[n].size()+one));
    }

    // model input file names
    for(n=0;n<*ntpl;n++)
    {
      itmp=infiles[n].size()+1;
      _file.write((char*)&itmp,sizeof(infiles[n].size()));
      _file.write(infiles[n].c_str(),(streamsize)(infiles[n].size()+one));
    }

    // instruction file names
    for(n=0;n<*nins;n++)
    {
      itmp=insfiles[n].size()+1;
      _file.write((char*)&itmp,sizeof(insfiles[n].size()));
      _file.write(insfiles[n].c_str(),(streamsize)(insfiles[n].size()+one));
    }

    // model output file names
    for(n=0;n<*nins;n++)
    {
      itmp=outfiles[n].size()+1;
      _file.write((char*)&itmp,sizeof(outfiles[n].size()));
      _file.write(outfiles[n].c_str(),(streamsize)(outfiles[n].size()+one));
    }

    return 1;
  }

  catch(int err)
  {

    if(err==ERROR_DIMENSIONS)
    {
      cout<<"\n\n   Error saving RUN message to file."<<'\n';
      cout<<"   *** Data and message sizes are inconsistent."<<'\n';
    }

    return 0;
  }

}

//****************************************************************************
int MODEL_RUN::set_frm_msg(MESSAGE *message)
//****************************************************************************
//CTM Feb 2011

/*
      Decompose message into necessary run parts
*/

{

  size_t n;
  size_t copied=0;
  size_t nsiz;
  string mr_chartostring(char *from);

  try
  {

#ifdef VERBOSE
    cout<<" *** MODEL_RUN object: routine 'set_frm_msg'... "<<'\n';
    cout<<" - Setting from message..."<<'\n';
#endif

#ifdef VERBOSE
    cout<<" - Dimensions"<<'\n';
#endif
    message->data->copyto(copied,(char*)id,sizeof(*id));
    copied+=sizeof(*id);
    cout<<" - Run ID: "<<*id<<'\n';
    //cout<<*id<<'\n';
    message->data->copyto(copied,(char*)nexec,sizeof(*nexec));
    copied+=sizeof(*nexec);
    //cout<<*nexec<<'\n';
    message->data->copyto(copied,(char*)npar,sizeof(*npar));
    copied+=sizeof(*npar);
    //cout<<*npar<<'\n';
    message->data->copyto(copied,(char*)nobs,sizeof(*nobs));
    copied+=sizeof(*nobs);
    //cout<<*nobs<<'\n';
    message->data->copyto(copied,(char*)ntpl,sizeof(*ntpl));
    copied+=sizeof(*ntpl);
    //cout<<*ntpl<<'\n';
    message->data->copyto(copied,(char*)nins,sizeof(*nins));
    copied+=sizeof(*nins);
    //cout<<*nins<<'\n';
    alloc(*npar,*nobs,*ntpl,*nins,*nexec,*id);

#ifdef VERBOSE
    cout<<" - Executables"<<'\n';
#endif
    // Executable properties
    for(n=0;n<*nexec;n++)
    {
      // Executable Name
      message->data->copyto(copied,(char*)&nsiz,sizeof(nsiz));
      copied+=sizeof(nsiz);
      if(nsiz>0)
      {
        char *ctmp=new char[nsiz];
        message->data->copyto(copied,ctmp,nsiz);
        copied+=nsiz;
        exec[n].name=mr_chartostring(ctmp);
        //cout<<"Exec Name: "<<exec[n].name<<'\n';
        delete ctmp;
      }
      // Executable Arguments
      message->data->copyto(copied,(char*)&nsiz,sizeof(nsiz));
      copied+=sizeof(nsiz);
      if(nsiz>0)
      {
        char *ctmp=new char[nsiz];
        message->data->copyto(copied,ctmp,nsiz);
        copied+=nsiz;
        exec[n].args=mr_chartostring(ctmp);
        //cout<<"Exec Args: "<<exec[n].args<<'\n';
        delete ctmp;
      }
      // Executable Path
      message->data->copyto(copied,(char*)&nsiz,sizeof(nsiz));
      copied+=sizeof(nsiz);
      if(nsiz>0)
      {
        char *ctmp=new char[nsiz];
        message->data->copyto(copied,ctmp,nsiz);
        copied+=nsiz;
        exec[n].path=mr_chartostring(ctmp);
        //cout<<"Exec Path: "<<exec[n].path<<'\n';
        delete ctmp;
      }
    }

    // Parameter Names
#ifdef VERBOSE
    cout<<" - Parameter names"<<'\n';
#endif
    for(n=0;n<*npar;n++)
    {
      message->data->copyto(copied,(char*)&nsiz,sizeof(nsiz));
      copied+=sizeof(nsiz);
      char *ctmp=new char[nsiz];
      message->data->copyto(copied,ctmp,nsiz);
      copied+=nsiz;
      parnams[n]=mr_chartostring(ctmp);
      delete ctmp;
    }

    // Parameter Values
#ifdef VERBOSE
    cout<<" - Parameter values"<<'\n';
#endif
    for(n=0;n<*npar;n++)
    {
      message->data->copyto(copied,(char*)&parvals[n],sizeof(parvals[n]));
      copied+=sizeof(parvals[n]);
      //cout<<parnams[n]<<" --> "<<parvals[n]<<'\n';
    }

    // Observation Names
#ifdef VERBOSE
    cout<<" - Observation names"<<'\n';
#endif
    for(n=0;n<*nobs;n++)
    {
      message->data->copyto(copied,(char*)&nsiz,sizeof(nsiz));
      copied+=sizeof(nsiz);
      char *ctmp=new char[nsiz];
      message->data->copyto(copied,ctmp,nsiz);
      copied+=nsiz;
      obsnams[n]=mr_chartostring(ctmp);
      delete ctmp;
    }

    // Observation Values
#ifdef VERBOSE
    cout<<" - Observation values"<<'\n';
#endif
    for(n=0;n<*nobs;n++)
    {
      message->data->copyto(copied,(char*)&obsvals[n],sizeof(obsvals[n]));
      copied+=sizeof(obsvals[n]);
      //cout<<obsnams[n]<<" --> "<<obsvals[n]<<'\n';
    }

    // Template File Names
#ifdef VERBOSE
    cout<<" - Template file names"<<'\n';
#endif
    for(n=0;n<*ntpl;n++)
    {
      message->data->copyto(copied,(char*)&nsiz,sizeof(nsiz));
      copied+=sizeof(nsiz);
      char *ctmp=new char[nsiz];
      message->data->copyto(copied,ctmp,nsiz);
      copied+=nsiz;
      tplfiles[n]=mr_chartostring(ctmp);
      delete ctmp;
    }

    // Model Input File Names
#ifdef VERBOSE
    cout<<" - Run input file names"<<'\n';
#endif
    for(n=0;n<*ntpl;n++)
    {
      message->data->copyto(copied,(char*)&nsiz,sizeof(nsiz));
      copied+=sizeof(nsiz);
      char *ctmp=new char[nsiz];
      message->data->copyto(copied,ctmp,nsiz);
      copied+=nsiz;
      infiles[n]=mr_chartostring(ctmp);
      //cout<<tplfiles[n]<<" --> "<<infiles[n]<<'\n';
      delete ctmp;
    }

    // Instruction File Names
#ifdef VERBOSE
    cout<<" - Instruction file names"<<'\n';
#endif
    for(n=0;n<*nins;n++)
    {
      message->data->copyto(copied,(char*)&nsiz,sizeof(nsiz));
      copied+=sizeof(nsiz);
      char *ctmp=new char[nsiz];
      message->data->copyto(copied,ctmp,nsiz);
      copied+=nsiz;
      insfiles[n]=mr_chartostring(ctmp);
      delete ctmp;
    }

    // Model Output File Names
#ifdef VERBOSE
    cout<<" - Run output file names"<<'\n';
#endif
    for(n=0;n<*nins;n++)
    {
      message->data->copyto(copied,(char*)&nsiz,sizeof(nsiz));
      copied+=sizeof(nsiz);
      char *ctmp=new char[nsiz];
      message->data->copyto(copied,ctmp,nsiz);
      copied+=nsiz;
      outfiles[n]=mr_chartostring(ctmp);
      //cout<<insfiles[n]<<" --> "<<outfiles[n]<<'\n';
      delete ctmp;
    }

    return 1;
  }

  catch(...)
  {
    cout<<"\n\n   Error setting RUN from message."<<'\n';
    return 0;
  }

}

//*****************************************************************************
string mr_chartostring(char *from)
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

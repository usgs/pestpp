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

#include "client.h"
#include "commandline.h"
#include "socket_utilities.h"
#include "defin.h"

using namespace std;

typedef struct addrinfo addr;

//****************************************************************************
int CLIENT::initialize(int &_argc,
                       char *_argv[])
//****************************************************************************
//CTM Dec 2010
/*
      Routine to initialize client properties
*/

{

  int ierr;
  COMMANDLINE cmd;
  SOCKET_UTILITIES sock_util;
  string s_temp;

  run_start=new list<clock_t>;
  run_end=new list<clock_t>;

  try
  {

    //pass commandline arguments to cmd
    cmd.setarg(_argc,_argv);
    //get the name of this node
    name=cmd.getarg("/name");
    if(name.compare("NOT_FOUND")==0)
    {
      cout<<"\n\n Client name not specified. Set to 'NOT_PROVIDED'."<<'\n';
      name="NOT_PROVIDED";
    }

    //get the type of this node
    s_temp=cmd.getarg("/type");
    if(s_temp.compare("NOT_FOUND")==0)
      s_temp.assign("slave");
    s_temp=lowercase(s_temp);
    if(s_temp.compare("model")==0)
      type=TYPE_MODEL;
    else if(s_temp.compare("slave")==0)
      type=TYPE_SLAVE;
    else if(s_temp.compare("face")==0)
      type=TYPE_FACE;
    else
      throw -5;

    //get the host socket
    host=cmd.getarg("/host");
    if(host.compare("NOT_FOUND")==0) throw -1;
    if(!sock_util.parsesocket(host,hostIP,hostPORT)) throw -4;
    hostStrPORT=sock_util.strport(hostPORT);

    // check if restart option is invoked
    s_temp=cmd.getarg("/restart");
    if(s_temp.compare("NOT_FOUND")==0)
      restart=false;
    else
      restart=true;
    //cout<<"Restart: |"<<s_temp<<"|"<<'\n';

    //initialize socket (ip and port)
    s_temp=cmd.getarg("/sock");
    if(s_temp.compare("NOT_FOUND")==0)
    {
      //check for ip
      myIP=cmd.getarg("/ip");
      if(myIP.compare("NOT_FOUND")==0)
      {
        myIP=sock_util.ip_automatic();
        if(myIP.compare("NOT_SET")==0) throw -2;
      }
      //check for port
      myStrPORT=cmd.getarg("/port");
      if(myStrPORT=="NOT_FOUND")
      {
        myPORT=sock_util.port_automatic();
        if(myPORT==NOT_SET) throw -3;
        myStrPORT=sock_util.strport(myPORT);
      }
      else
        myPORT=sock_util.uintport(myStrPORT);
    }
    else
    {
      ierr=sock_util.parsesocket(s_temp,myIP,myPORT);
      myStrPORT=sock_util.strport(myPORT);
      if(ierr<0) throw -4;
    }
    cout<<" - CLIENT using IP: "<<myIP<<":"<<myPORT<<'\n';

    //initialize additional variables
    rating=-1.;
    status=WAITING_FOR_RUN;
    terminate=false;

    return 1;
  }
  catch (int err)
  {

    delete run_start;
    delete run_end;

    if(err==-1)
    {
      cout<<"\n\n *** Unable to determine HOST SOCKET."<<'\n';
      cout<<"     Specify HOST using '/host' switch."<<'\n';
    }
    if(err==-2)
    {
      cout<<"\n\n *** Unable to determine IP."<<'\n';
      cout<<"     Please verify internet connection or, specify IP using '/ip' switch."<<'\n';
    }
    else if(err==-3)
    {
      cout<<"\n\n *** Unable to determine PORT."<<'\n';
      cout<<"     Please verify internet connection or, specify PORT using '/port' switch."<<'\n';
    }
    else if(err==-4)
      cout<<"\n\n *** Error parsing supplied socket for IP and PORT."<<'\n';
    if(err==-5)
    {
      cout<<"\n\n *** Unable to determine slave TYPE."<<'\n';
      cout<<"     Specify TYPE using '/type' switch."<<'\n';
      cout<<"     Avaliable Types:"<<'\n';
      cout<<"        'MODEL' - Provide runs to run manager"<<'\n';
      cout<<"        'SLAVE' - Slave"<<'\n';
      cout<<"        'FACE' - Interface program"<<'\n';
    }

    cout<<"\n     See manual for more details."<<'\n';
    return 0;
  }

}

//****************************************************************************
int CLIENT::start()
//****************************************************************************
// CTM Dec 2010
// CTM Apr 2011 - Added retry functionality
/*
      Routine to open and connect CLIENT to HOST
*/

{

  //routine specific variables
  int ntry;
  int flag=1;
  addr *sa,hint;
  sockaddr_in serverINF;

  void wait(int seconds);

  //parse socket into port and ip strings
  try
  {
    //initialize hint structure
    memset(&hint,0,sizeof(hint));
    hint.ai_family=AF_INET;
    hint.ai_socktype=SOCK_STREAM;
    hint.ai_protocol=IPPROTO_TCP;
    hint.ai_flags=AI_PASSIVE;

    //initialize sa and convert PORT to address type
    if(getaddrinfo(myIP.c_str(),myStrPORT.c_str(),&hint,&sa)!=0) throw -1;
    //open socket
    sock=socket(sa->ai_family,sa->ai_socktype,sa->ai_protocol);
    if(sock==INVALID_SOCKET) throw -2;
    //make socket non-blocking
    //if(!makenonblocking()) throw -4;
    // connect to host
    if(setsockopt(sock,SOL_SOCKET,SO_KEEPALIVE,(char *)&flag,sizeof(int))!=0) throw -3;
    int nTimeout = 5000; // 5 seconds
    setsockopt(sock,SOL_SOCKET,SO_RCVTIMEO,(const char*)&nTimeout,sizeof(int));
    serverINF.sin_family=AF_INET;
    serverINF.sin_addr.s_addr=inet_addr(hostIP.c_str());
    serverINF.sin_port=htons(hostPORT);
    for(ntry=0;ntry<CXN_RETRY_MAX;ntry++)
    {
      cout<<" - Connection to GMAN attempt: "<<ntry+1<<'\n';
      if(connect(sock,(SOCKADDR *)&serverINF,sizeof(serverINF))==0) {
        break;
      }
      if(WSAGetLastError()==10056||WSAGetLastError()==10035) {
        break;
      }
      closesocket(sock);
      wait(CXN_RETRY_TIME);
    }
    if(ntry>=CXN_RETRY_MAX) throw -3;

    cout<<" - CLIENT connected to GMAN on socket: "<<sock<<'\n';

    return 1;

  }
  catch (int err)
  {
    if(err==-1)
      cout<<"\n\n Error converting IP and PORT to ADDRESS."<<'\n';
    if(err==-2)
      cout<<"\n\n Error opening SOCKET."<<'\n';
    if(err==-3)
      cout<<"\n\n Error connecting to HOST."<<'\n';
    if(err==-4)
      cout<<" *** Error setting SOCKET option."<<'\n';

    return 0;
  }
}

//****************************************************************************
void CLIENT::setsocket(SOCKET _socket)
//****************************************************************************
//CTM Dec 2010

/*
      Routine to set client socket and add to file descriptor set
*/

{

  sock=_socket;

}

//*****************************************************************************
string CLIENT::lowercase(string s)
//*****************************************************************************
//CTM, Nov 2010

/*
Routine to convert string to lower case
*/

{
  for(unsigned int i=0;i<s.size();i++) {
    s[i]=tolower(s[i]);
  }
  return s;
}

//******************************************************************************************************
void wait(int seconds)
//from http://www.cplusplus.com/reference/clibrary/ctime/clock/
{
  clock_t endwait;
  endwait = clock () + seconds * CLOCKS_PER_SEC ;
  while (clock() < endwait) {}
}

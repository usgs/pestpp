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
#include <string>
#include <time.h>
#include <winsock2.h>
#include <ws2tcpip.h>

#include "socket_utilities.h"

using namespace std;

typedef struct addrinfo addr;

//****************************************************************************
string SOCKET_UTILITIES::ip_automatic()
//****************************************************************************
//CTM Dec 2010
/*
      Automatically gets local IP address
*/

{

  //char buffer[512];
  SOCKET sock;
  string s1,s2;

  string su_chartostring(char *from);
  string su_trim(string *s);

  try
  {
    sock=socket(AF_INET,SOCK_DGRAM,0);
    if(sock==SOCKET_ERROR) throw -1;
    const char* kGoogleDnsIp = "8.8.8.8";
    unsigned short kDnsPort = 53;
    struct sockaddr_in serv;
    memset(&serv, 0, sizeof(serv));
    serv.sin_family = AF_INET;
    serv.sin_addr.s_addr = inet_addr(kGoogleDnsIp);
    serv.sin_port = htons(kDnsPort);
    int err = connect(sock, (const sockaddr*) &serv, sizeof(serv));
    if(err==-1) throw -1;
    sockaddr_in name;
    socklen_t namelen = sizeof(name);
    err = getsockname(sock, (sockaddr*) &name, &namelen);
    if(err==-1) throw -1;
//#if (_WIN32_WINNT >= 0x600)
//    const char* p = inet_ntop(AF_INET, &name.sin_addr, buffer, 512);
//    return trim(chartostring(buffer));
//#else
    char *p = inet_ntoa(name.sin_addr);
    closesocket(sock);
    s1=su_chartostring(p);
    s2=su_trim(&s1);
    return s2;
//#endif
  }
  catch (...)
  {
    cout<<"\n\n Error determining IP automatically."<<'\n';
    closesocket(sock);
    return "NOT_SET";
  }

}

//*****************************************************************************
unsigned int SOCKET_UTILITIES::port_automatic()
//*****************************************************************************

/*
    Routine to generate a random PORT number
*/

{

  int ntry=0;
  time_t seconds;
  unsigned int port;
  int i;

  // initialize
  seconds=time(NULL);
  srand((unsigned int) seconds);
  rand();

  // get random port number
  while(1)
  {
    ntry++;
    i=rand();
    port=i%(MAX_PORT-MIN_PORT+1)+MIN_PORT;
    if(port_check(port)||ntry==5) break;
  }
  if(ntry>=5) port=NOT_SET;

  return port;

}

//*****************************************************************************
bool SOCKET_UTILITIES::port_check(unsigned int &_port)
//*****************************************************************************
//CTM, Nov 2010

/*
    Routine to check whether supplied port is valid
*/

{

  int ierr;
  SOCKET socketTEMP;
  sockaddr_in serverINF;

  //ensure port is between 1024 and 65535
  if(_port>65535||_port<1024)
    return false;

  //verify socket can be connected via this port
  socketTEMP=socket(AF_INET,SOCK_STREAM,IPPROTO_TCP);
  ierr=WSAGetLastError();
  if(socketTEMP==INVALID_SOCKET)
    return false;

  //Bind Socket to Port
  serverINF.sin_family=AF_INET;
  serverINF.sin_addr.s_addr=INADDR_ANY;
  serverINF.sin_port=htons(_port);
  ierr=bind(socketTEMP,(SOCKADDR*)&serverINF,sizeof(serverINF));
  if(ierr==SOCKET_ERROR)
    return false;
  closesocket(socketTEMP);

  return true;

}

//*****************************************************************************
unsigned int SOCKET_UTILITIES::uintport(string &from)
//*****************************************************************************
//CTM, Nov 2010

/*
      Routine to convert string to integer
*/

{
  stringstream ss;
  int to;

  ss<<from;
  ss>>to;
  return to;
}

//*****************************************************************************
string SOCKET_UTILITIES::strport(unsigned int &from)
//*****************************************************************************
//CTM, Nov 2010

/*
    Routine to convert unsigned int to string
*/

{
  stringstream ss;

  ss<<from;
  return ss.str();
}

//*****************************************************************************
int SOCKET_UTILITIES::parsesocket(string &_socket,
                                  string &_ip,
                                  unsigned int &_port)
//*****************************************************************************
//CTM, Nov 2010

/*
    Routine to parse a socket in the form ip:port into strings of ip and port
*/

{

  string s_port;

  size_t i=_socket.find_first_of(":");
  if(i<0) return 0;
  //split string on ':'
  _ip=_socket.substr(0,i);
  s_port=_socket.substr(i+1,_socket.size());
  _port=uintport(s_port);
  return 1;
}

//*****************************************************************************
string su_chartostring(char *from)
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

//*****************************************************************************
string su_trim(string *s)
//*****************************************************************************
//CTM, Nov 2010

/*
    Routine to trim leading and trailing blank spaces from a string
*/

{
  size_t istart=s->find_first_not_of(" /t");
  size_t iend=s->find_last_not_of(" /t");
  if(istart==string::npos) {
    return "\0";
  }
  else {
    return s->substr(istart,iend-istart+1);
  }
}

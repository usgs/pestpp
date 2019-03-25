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

#ifndef MESSAGE_H
#define MESSAGE_H

#include <winsock2.h>
#include <ws2tcpip.h>

#include "defin.h"
#include "header.h"
#include "buffer.h"

//****************************************************************************
class MESSAGE
//****************************************************************************
{

  public:

      size_t msgsize;
      size_t datasize;
      HEADER header;
      BUFFER *data;

      void alloc();
      void dealloc();
      void setmsgsize(size_t _size);
      void setdatasize(size_t _size);
      int waitformessage(SOCKET &_socket,long sec,long usec);
      int receiveheader(SOCKET &_socket);
      int receivedata(SOCKET &_socket);
      int sendme(SOCKET &_socket);

};

#endif /* MESSAGE_H */

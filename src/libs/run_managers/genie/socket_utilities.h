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

#ifndef SOCKUTILS_H
#define SOCKUTILS_H

#include <string>
#include <windows.h>

#include "defin.h"

//****************************************************************************
class SOCKET_UTILITIES
//****************************************************************************
{
  public:

      std::string ip_automatic();
      unsigned int port_automatic();
      unsigned int uintport(std::string &from);
      std::string strport(unsigned int &from);
      int parsesocket(std::string &_socket,
                      std::string &_ip,
                      unsigned int &_port);
      bool port_check(unsigned int &_port);

  private:


};

#endif /* SOCKUTILS_H */

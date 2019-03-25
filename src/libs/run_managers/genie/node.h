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

#ifndef NODE_H
#define NODE_H

#include <string>
#include <list>
#include <winsock2.h>
#include <ws2tcpip.h>
#include <time.h>

#include "defin.h"
#include "message.h"
#include "runprocess.h"
#include "modelrun.h"

//****************************************************************************
class NODE
//****************************************************************************
{

  private:

  public:

      bool terminate,restart;
      int type,status;
      double interval;
      double rating;
      std::string name;
      SOCKET sock;
      std::list<clock_t> *run_start;
      std::list<clock_t> *run_end;

      virtual int initialize(int &_argc,
                             char *_argv[])
                  {return (0);}
      virtual int start()
                  {return (0);}
      void init_timers();
      void stop();
      int makenonblocking();
      void settype(int _type);
      void setprocess(RUNPROCESS *_process);
      void endprocess();
      int sendtype();
      int sendspeed();
      int sendname();
      int sendcommand(int _command);
      int sendstatus(int _status);
      int sendnslave(int nslave);

      double median_runtime();
      double expected_timeleft();

  protected:

      RUNPROCESS *myprocess;
      std::string myIP,myStrPORT;
      unsigned int myPORT;

};

#endif /* NODE_H */

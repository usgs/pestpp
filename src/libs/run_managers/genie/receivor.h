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

#ifndef RECEIVOR_H
#define RECEIVOR_H

#include <list>
//#include <vector>

#include "thread.h"
#include "client.h"
#include "message.h"
#include "modelrun.h"
#include "modelresult.h"

//****************************************************************************
class RECEIVOR : public THREAD
//****************************************************************************
//CTM Dec 2010

/*
      Thread class to receive messages on a dedicated socket
*/

{

  public:

      RECEIVOR();
      ~RECEIVOR();

      CLIENT *slave;
      bool *haveresult;

      int execute();
      void set_slave(CLIENT *_client);
      void set_run(MODEL_RUN *_run);
      void set_queue_runs(std::list<MODEL_RUN*> *_runs);
      void set_queue_results(std::list<MODEL_RESULT*> *_results);
      MODEL_RESULT *get_result();
      MODEL_RUN *get_run();
      void clear_run();

  private:

      MESSAGE *message;
      MODEL_RUN *run;
      std::list<MODEL_RUN*> *runs;
      MODEL_RESULT *result;
      std::list<MODEL_RESULT*> *results;
      bool endthread;

      void process_type();
      void process_speed();
      void process_name();
      void process_command();
      void process_status();
      int process_run();
      int process_obs();

};

#endif // RECEIVOR_H

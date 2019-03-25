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

#ifndef MODELRESULT_H
#define MODELRESULT_H

#include <string>

#include <winsock2.h>
#include <ws2tcpip.h>

#include "executable.h"
#include "defin.h"
#include "message.h"
#include "modelrun.h"

//****************************************************************************
class MODEL_RESULT
//****************************************************************************
{

  public:

    size_t *npar,*nobs;
    size_t *id;
    double *parvals,*obsvals;

    MODEL_RESULT();
    ~MODEL_RESULT();

    //void init();
    void alloc(); //size_t &_npar,size_t &_nobs,int &_id);
    void dealloc();
    int set_frm_msg(MESSAGE *message);
    void set_frm_run(MODEL_RUN *run);
    int send(SOCKET &_socket);

  private:

};

#endif /* MODELRESULT_H */

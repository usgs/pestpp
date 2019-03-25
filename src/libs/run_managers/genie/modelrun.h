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

#ifndef MODELRUN_H
#define MODELRUN_H

#include <string>
#include <fstream>

#include <winsock2.h>
#include <ws2tcpip.h>

#include "executable.h"
#include "defin.h"
#include "message.h"

//****************************************************************************
class MODEL_RUN
//****************************************************************************
{

  public:

    // properties

    size_t *npar,*nobs,*ntpl,*nins,*nexec;
    size_t *id;
    std::string *tplfiles,*insfiles;
    std::string *infiles,*outfiles;
    std::string *parnams,*obsnams;
    double *parvals,*obsvals;
    EXECUTABLE *exec;

    // methods

    MODEL_RUN();

    void init();
    void alloc(size_t &_npar,size_t &_nobs,
               size_t &_ntpl,size_t &_nins,
               size_t &_nexec, size_t &_id);
    void dealloc();

    int set_exec(std::string *exec_ary);

    int send(SOCKET &_socket);
    int set_frm_msg(MESSAGE *message);
    int save(std::fstream &_file);

#ifdef GSLAVE
    int write_input();
    int read_output();
    int delete_output();
#endif

  private:

    //std::string chartostring(char *from);


};

#endif /* MODELRUN_H */

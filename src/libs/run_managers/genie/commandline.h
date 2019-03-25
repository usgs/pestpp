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

#ifndef COMMANDLINE_H
#define COMMANDLINE_H

#include <list>
#include <string>

//*****************************************************************************
class COMMANDLINE
//*****************************************************************************
{

  public:

      void setarg(int _argc,char *_argv[]);
      std::string getarg(std::string _flag);

  private:

      std::list<std::string> args;
      std::list<std::string> pars;

};

#endif /* COMMANDLINE_H */

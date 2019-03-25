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

#ifndef EXECUTABLE_H
#define EXECUTABLE_H

#include <string>

//****************************************************************************
class EXECUTABLE
//****************************************************************************
//CTM Jan 2011

/*
      run executable class
*/

{

  private:

      //std::string appendseparator(std::string &_string);
      //std::string prependspace(std::string &_string);
      //std::string chartostring(char _char);

  public:

      std::string path;
      std::string name;
      std::string args;

      int initialize(int &_argc,char *_argv[]);

      //std::string doubleseparators(std::string &_string);

};

#endif //EXECUTABLE_H

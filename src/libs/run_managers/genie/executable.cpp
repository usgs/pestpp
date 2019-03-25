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

#include "executable.h"
#include "commandline.h"

using namespace std;

//****************************************************************************
int EXECUTABLE::initialize(int &_argc,char *_argv[])
//****************************************************************************
//CTM Jan 2011
/*
      Routine to initialize run executable properties
*/

{

  COMMANDLINE cmd;

  string appendseparator(string &_string);
  //string doubleseparators(string &_string);
  string prependspace(string &_string);

  try
  {

    // pass commandline arguments to cmd
    cmd.setarg(_argc,_argv);

    // get the path to the executable
    path=cmd.getarg("/path");
    if(name.compare("NOT_FOUND")==0)
      path='\0';
    else
    {
      path=appendseparator(path);
      //path=doubleseparators(path);
    }
    //cout<<path<<'\n';

    // get any arguments
    args=cmd.getarg("/args");
    if(args.compare("NOT_FOUND")==0)
      args='\0';
    else
      args=prependspace(args);
    //cout<<"|"<<args<<"|"<<'\n';

    // get name of executable - fail if not provided
    name=cmd.getarg("/exec");
    if(name.compare("NOT_FOUND")==0)
    {
      name='\0';
      //throw -1;
    }
    //cout<<"|"<<name<<"|"<<'\n';

    return 0;
  }

  catch(...)
  {
    cout<<"\n\n Name of run executable not specified.";
    cout<<"Specify name using '/exec' flag."<<'\n';
    return 0;
  }
}

//****************************************************************************
string appendseparator(string &_string)
//****************************************************************************
//CTM Jan 2011
/*
      Routine to add folder separators to end of path if not present
*/

{

  string lastchar;

  lastchar=_string.at(_string.size()-1);
  if(lastchar.compare("\\")!=0)
    _string.push_back('\\');
  return _string;
}

//****************************************************************************
string doubleseparators(string &_string)
//****************************************************************************
//CTM Jan 2011
/*
      Routine to double folder separators (\ to \\)
*/

{

  string::iterator _char;
  string ex_chartostring(char _char);

  _char=_string.begin();
  while(_char!=_string.end())
  {
    if(ex_chartostring(*_char).compare("\\")==0)
    {
      _char=_string.insert(_char,'\\');
      _char++;
    }
    _char++;
  }

  return _string;
}

//****************************************************************************
string prependspace(string &_string)
//****************************************************************************
//CTM Jan 2011
/*
      Routine to add folder separators to end of path if not present
*/

{

  string ex_chartostring(char _char);

  if(ex_chartostring(*_string.begin()).compare(" ")!=0)
    _string.insert(_string.begin(),' ');

  return _string;
}

//*****************************************************************************
string ex_chartostring(char _char)
//*****************************************************************************
//CTM, Nov 2010

/*
      Routine to convert char to string
*/

{
  stringstream ss;

  ss<<_char;
  return ss.str();
}

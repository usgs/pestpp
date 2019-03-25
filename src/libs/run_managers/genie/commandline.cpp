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
#include <cmath>
#include <algorithm>

#include "commandline.h"

using namespace std;

//****************************************************************************
void COMMANDLINE::setarg(int _argc,char *_argv[])
//****************************************************************************
//CTM Dec 2010
/*
      Routine to set the command line arguments
*/

{

  string firstchar;
  string _arg;
  string _par;
  bool havearg=false;
  //list<string>::iterator _a;
  //list<string>::iterator _p;

  string cl_chartostring(char *_char);
  string cl_lowercase(string s);
  string cl_trim(string &_string);

  // populate lists of arguments and corresponding parameters
  for(int i=1;i<_argc;i++)
  {
    _arg=cl_chartostring(_argv[i]);
    firstchar=_arg.at(0);
    if(firstchar.compare("/")==0)
    {
      if(havearg) pars.push_back(cl_trim(_par));
      _par.clear();
      args.push_back(cl_lowercase(_arg));
      havearg=true;
    }
    else
      _par+=" "+_arg;
  }
  pars.push_back(cl_trim(_par));

  //_a=args.begin();
  //_p=pars.begin();
  //while(_a!=args.end())
  //{
  //  cout<<*_a<<", "<<*_p<<'\n';
  //  _a++;
  //  _p++;
  //}

}

//****************************************************************************
string COMMANDLINE::getarg(string _flag)
//****************************************************************************
//CTM Jan 2011

/*
      Routine to return the command line prompt for a particular flag
*/

{

  list<string>::iterator _a;
  list<string>::iterator _p;

  _a=args.begin();
  _p=pars.begin();
  while(_a!=args.end())
  {
    if(_flag.compare(*_a)==0) return *_p;
    _a++;
    _p++;
  }
  return "NOT_FOUND";

}

//*****************************************************************************
string cl_chartostring(char *_char)
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

//*****************************************************************************
string cl_lowercase(string s)
//*****************************************************************************
//CTM, Nov 2010

/*
Routine to convert string to lower case
*/

{
  for(unsigned int i=0;i<s.size();i++) {
    s[i]=tolower(s[i]);
  }
  return s;
}

//*****************************************************************************
string cl_trim(string &_string)
//*****************************************************************************
//CTM, Nov 2010

/*
      Routine to trim leading and trailing spaces from a string
*/

{


  size_t is,ie;
  string _space=" ";
  string _tab="\t";

  is=max(_string.find_first_not_of(_space),
         _string.find_first_not_of(_tab));
  ie=max(_string.find_last_not_of(_space),
         _string.find_last_not_of(_tab));

  if(is==string::npos)
    return "\0";
  else
    return _string.substr(is,ie-is+1);

}

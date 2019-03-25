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

#include <string>
#include <sstream>
#include <algorithm>

using namespace std;

//*****************************************************************************
string trim(string &_string)
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

//*****************************************************************************
string chartostring(char *from)
//*****************************************************************************
//CTM, Nov 2010

/*
Routine to convert char to string
*/

{
  stringstream ss;

  ss<<from;
  return ss.str();
}

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

#ifndef BUFFER_H
#define BUFFER_H

//#include <winsock2.h>

//****************************************************************************
class BUFFER
//****************************************************************************
{

public:

      int size;
      char *buf;

      void alloc(size_t _size);
      void dealloc();
      void copyfrom(size_t _start, const char *_src, size_t _size);
      void copyto(size_t _start, char *_dst, size_t _size);

};

#endif /* BUFFER_H */

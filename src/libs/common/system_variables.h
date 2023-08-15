/*


    This file is part of PEST++.

    PEST++ is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PEST++ is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PEST++.  If not, see<http://www.gnu.org/licenses/>.
*/
#ifndef SYSTEM_VARIABLES_H_
#define SYSTEM_VARIABLES_H_

#include "config_os.h"
#include <string>

class OperSys
{
public:
	static const std::string DIR_SEP;
	static const std::string COMMAND_LINE_APPEND;
	void string2pathname(std::string &s);
	static std::string getcwd();
	static void chdir(const char *str);
	static char *gets_s(char *str, size_t len);
	static bool double_is_invalid(double x);
};

#ifdef OS_WIN
#include <Windows.h>
PROCESS_INFORMATION start(std::string &cmd_string);
#endif
#ifdef OS_LINUX
int start(std::string &cmd_string);
#endif




#endif /* SYSTEM_VARIABLES_H_ */

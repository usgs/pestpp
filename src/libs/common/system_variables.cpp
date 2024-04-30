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

#include <cstdlib>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <iostream>
#include "system_variables.h"

#ifdef OS_WIN
 #include <direct.h>
#endif

#ifdef OS_LINUX
#include "stdio.h"
#include <fcntl.h>
#include <unistd.h>

#endif

#ifdef OS_WIN
const std::string OperSys::DIR_SEP = "\\";
const std::string OperSys::COMMAND_LINE_APPEND = " ; ";
#endif

#ifdef OS_LINUX
const std::string OperSys::DIR_SEP = "/";
const std::string OperSys::COMMAND_LINE_APPEND = " & ";
#endif

using namespace std;

void OperSys::string2pathname(string &s)
{
	size_t i;
	size_t len(s.size());
	string dir_chars("/\\");
	stringstream new_s;

	for (i=0; i<len; ++i) {
		if (s[i] == '/') {new_s << DIR_SEP;}
		else {new_s << s[i];}
		if (i!=len-1 && (s[i] == '\\' || s[i]=='/')) {
			new_s << DIR_SEP;
		}
	}
	s = new_s.str();
}

string OperSys::getcwd()
{
    #ifdef OS_WIN
	char *buffer;
	buffer = _getcwd( NULL, 0 );
	string cwd(buffer);
    free(buffer);
	return cwd;
    #endif
    #ifdef OS_LINUX
        char *buffer;
	buffer = ::getcwd( NULL, 0 );
	string cwd(buffer);
        free(buffer);
	return cwd;
    #endif
}

void OperSys::chdir(const char *str)
{
   #ifdef OS_WIN
      _chdir(str);
   #endif
   #ifdef OS_LINUX
      //chdir(str);
   #endif
}

char* OperSys::gets_s(char *str, size_t len)
{
 #ifdef OS_WIN
  return ::gets_s(str, len);
 #endif
 #ifdef OS_LINUX
  return gets_s(str, len);
 #endif

}

bool OperSys::double_is_invalid(double x)
{
	bool test = false;
#if defined __INTEL_COMPILER && defined __APPLE__
    test = (isnan(x) || isinf(x));
#endif
#if defined __INTEL_COMPILER && defined OS_LINUX && !defined __APPLE__
     test = (::isnan(x) || ::isinf(x));
#endif
#if defined OS_WIN
     test = (std::isnan(x) || !std::isfinite(x));
#endif
	return test;
}


#ifdef OS_WIN
PROCESS_INFORMATION start(string &cmd_string)
{
	char* cmd_line = _strdup(cmd_string.c_str());
	STARTUPINFO si;
	PROCESS_INFORMATION pi;
	ZeroMemory(&si, sizeof(si));
	ZeroMemory(&pi, sizeof(pi));
	if (!CreateProcess(NULL, cmd_line, NULL, NULL, false, 0, NULL, NULL, &si, &pi))
	{
		std::string cmd_string(cmd_line);
		throw std::runtime_error("CreateProcess() failed for command: " + cmd_string);
	}
	delete cmd_line;
	return pi;
}
#endif


#ifdef OS_LINUX
int start(string &cmd_string)
{
	//split cmd_string on whitespaces
	stringstream cmd_ss(cmd_string);
	string cmd;
	vector<string> cmds;
	while (getline(cmd_ss, cmd,' '))
	{
		cmds.push_back(cmd);
	}

	//create the strurture for execv
	vector<char const*> arg_v;
	for (size_t icmd = 0; icmd<cmds.size(); ++icmd)
	{
		arg_v.push_back(cmds[icmd].data());
	}
	//char * const*argv = new char* [cmds.size()+1];
	//for (size_t icmd=0; icmd<cmds.size(); ++icmd)
	//{
	//  argv[icmd] = cmds[icmd].data();
	//}
	//argv[cmds.size() + 1] = NULL; //last arg must be NULL
   
	arg_v.push_back(NULL);
	pid_t pid = fork();
	if (pid == 0)
	{
		setpgid(0, 0);
//        int fd = open("stdout.dat", O_CREAT);
//        dup2(fd, 1);
//        std::cout << "file descrp " << fd << std::endl;
		int success = execvp(arg_v[0], const_cast<char* const*>(&(arg_v[0])));
		if (success == -1)
		{
			throw std::runtime_error("execv() failed for command: " + cmd_string);
		}

	}
	else
	{
		setpgid(pid, pid);
	}
	return pid;
}

#endif



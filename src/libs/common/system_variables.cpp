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
/**
 * @file system_variables.cpp
 * @brief Implementation of system_variables.
 */


#include <cstdlib>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <sstream>
#include <cmath>
#include <vector>
#include <iostream>
#include <filesystem>
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

// OS_MAC has to be tested before OS_LINUX: config_os.h defines BOTH on a mac, so an
// #ifdef OS_LINUX branch catches macs too and would send them looking for /proc, which
// is not there.
#if defined(OS_WIN)
// Windows.h arrives through system_variables.h
#elif defined(OS_MAC)
#include <sys/sysctl.h>
#include <mach/mach.h>
#else
#include <fstream>
#endif

#if !defined(OS_WIN) && !defined(OS_MAC)
/** The single number in a one-line file, or -1 if it is not there or is not a number.
 *
 * The cgroup files this reads hold the word "max" when no limit is set, which is why a failed
 * conversion has to be an ordinary answer here rather than an error.
 */
static long long read_one_number(const string& path)
{
	ifstream f(path);
	if (!f.good())
		return -1;
	string tok;
	if (!(f >> tok))
		return -1;
	try
	{
		return stoll(tok);
	}
	catch (...)
	{
		return -1;
	}
}
#endif

/**
 * @brief Physical memory on this machine, total and available, in bytes.
 */
SysMemory get_system_memory()
{
	SysMemory m;
	m.valid = false;
	m.total_bytes = 0;
	m.available_bytes = 0;

#if defined(OS_WIN)
	MEMORYSTATUSEX s;
	s.dwLength = sizeof(s);
	if (GlobalMemoryStatusEx(&s))
	{
		m.total_bytes = (long long)s.ullTotalPhys;
		m.available_bytes = (long long)s.ullAvailPhys;
		m.valid = true;
	}

#elif defined(OS_MAC)
	uint64_t total = 0;
	size_t len = sizeof(total);
	if (sysctlbyname("hw.memsize", &total, &len, nullptr, 0) == 0)
	{
		m.total_bytes = (long long)total;
		m.valid = true;
	}
	// free pages alone would badly understate this. the system parks memory in inactive and
	// purgeable pages that it hands straight back when something asks, so those are available
	// in every sense that matters here.
	vm_size_t page_size = 0;
	vm_statistics64_data_t vmstat;
	mach_msg_type_number_t count = HOST_VM_INFO64_COUNT;
	mach_port_t self = mach_host_self();
	if ((host_page_size(self, &page_size) == KERN_SUCCESS) &&
	    (host_statistics64(self, HOST_VM_INFO64, (host_info64_t)&vmstat, &count) == KERN_SUCCESS))
	{
		m.available_bytes = ((long long)vmstat.free_count + (long long)vmstat.inactive_count +
		                     (long long)vmstat.purgeable_count) * (long long)page_size;
	}

#else
	ifstream f("/proc/meminfo");
	if (f.good())
	{
		// read whole lines rather than streaming key/value/unit: not every line carries a
		// unit, and one that does not would knock a field-by-field read out of step for
		// everything after it
		string line;
		long long mem_free = -1;
		while (getline(f, line))
		{
			istringstream ss(line);
			string key;
			long long value = 0;
			if (!(ss >> key >> value))
				continue;
			if (key == "MemTotal:")
				m.total_bytes = value * 1024;
			else if (key == "MemAvailable:")
				m.available_bytes = value * 1024;
			else if (key == "MemFree:")
				mem_free = value * 1024;
		}
		// MemAvailable is absent on kernels older than 3.14. MemFree is a poorer answer - it
		// ignores cache that would be handed back - but it beats reporting nothing at all.
		if ((m.available_bytes == 0) && (mem_free > 0))
			m.available_bytes = mem_free;
		m.valid = (m.total_bytes > 0);
	}

	// a container or a scheduler caps memory below what the hardware has, and /proc/meminfo
	// knows nothing about it - it describes the whole host. without this an agent inside a
	// 4gb container cheerfully reports the machine's 512gb, then gets killed for using 5.
	long long limit = read_one_number("/sys/fs/cgroup/memory.max");                 // cgroup v2
	long long usage = read_one_number("/sys/fs/cgroup/memory.current");
	if (limit < 0)
	{
		limit = read_one_number("/sys/fs/cgroup/memory/memory.limit_in_bytes");     // cgroup v1
		usage = read_one_number("/sys/fs/cgroup/memory/memory.usage_in_bytes");
	}
	// v1 with no limit set does not say so, it reports a huge sentinel instead, so a "limit"
	// at or above what the machine has is not a real cap and is ignored
	if ((limit > 0) && ((m.total_bytes == 0) || (limit < m.total_bytes)))
	{
		m.total_bytes = limit;
		if (usage >= 0)
			m.available_bytes = (limit > usage) ? (limit - usage) : 0;
		m.valid = true;
	}
#endif

	return m;
}

/**
 * @brief Disk space on the filesystem holding a path, total and available, in bytes.
 */
SysStorage get_system_storage(const string& path)
{
	SysStorage s;
	s.valid = false;
	s.total_bytes = 0;
	s.available_bytes = 0;

	// the one part of this that needs no per-os code at all. the error_code form on purpose:
	// the throwing one raises if the path has gone away, and a directory disappearing under a
	// run is exactly the moment this gets asked
	error_code ec;
	filesystem::space_info si = filesystem::space(path, ec);
	if (ec)
		return s;

	// static_cast because these come back as uintmax_t, and a filesystem that cannot answer
	// reports -1 cast to unsigned - which would arrive as an enormous positive number
	if ((si.capacity == static_cast<uintmax_t>(-1)) ||
	    (si.available == static_cast<uintmax_t>(-1)))
		return s;

	s.total_bytes = (long long)si.capacity;
	s.available_bytes = (long long)si.available;
	s.valid = (s.total_bytes > 0);
	return s;
}

/**
 * @brief String2pathname.
 *
 * @param s Description.
 */
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

/**
 * @brief Getcwd.
 *
 * @return Description.
 */
string OperSys::getcwd()
{
	// NOTE: deliberately NOT std::filesystem::current_path().string(). On windows that
	// returns UTF-8 bytes, while _getcwd returns active-code-page bytes. The non-ascii cwd
	// guard in CmdLine (utilities.cpp) compares narrow bytes, so switching encodings changes
	// when it fires - it made the panther agent abort at startup in a non-ascii directory
	// instead of connecting to the master. Keep the platform-native narrow calls.
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

/**
 * @brief Chdir.
 *
 * @param str Description.
 */
void OperSys::chdir(const char *str)
{
   #ifdef OS_WIN
      _chdir(str);
   #endif
   #ifdef OS_LINUX
      //chdir(str);
   #endif
}

/**
 * @brief Gets s.
 *
 * @param str Description.
 * @param len Description.
 *
 * @return Description.
 */
char* OperSys::gets_s(char *str, size_t len)
{
 #ifdef OS_WIN
  return ::gets_s(str, len);
 #endif
 #ifdef OS_LINUX
  return gets_s(str, len);
 #endif

}

/**
 * @brief Double is invalid.
 *
 * @param x Description.
 *
 * @return Description.
 */
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
/**
 * @brief Start.
 *
 * @param cmd_string Description.
 *
 * @return Description.
 */
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



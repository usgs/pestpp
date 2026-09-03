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
	static bool double_is_invalid(double x);
};

/** What the machine has, and what it would actually give you right now, in bytes.
 *
 * `available` is the number worth acting on. It is not "free" - a healthy machine keeps very
 * little memory untouched, because anything spare is holding cache it can drop the moment
 * something asks. Free memory on a busy machine reads near zero while gigabytes are there for
 * the taking, so acting on it would have you back off when there was no need to.
 *
 * `valid` is false when none of it could be read. The numbers are zero then and mean nothing -
 * check it before using them, rather than treating zero as "no memory".
 */
struct SysMemory
{
	bool valid;
	long long total_bytes;
	long long available_bytes;
};

/** Physical memory on this machine, or the cap it is being held to.
 *
 * Under a container or a scheduler like slurm the cap is what matters, not what the hardware
 * has - see the comments in the linux branch.
 */
SysMemory get_system_memory();

/** Disk space on the filesystem holding a given path, in bytes.
 *
 * `available` is what this process could actually write. It is not the same as free: filesystems
 * hold some back for root, so free is the larger and more optimistic number, and a non-privileged
 * agent that believed it would run out early.
 *
 * Ask about the working directory rather than the machine. An agent's directory is often on a
 * different mount from `/`, and it is the one that fills up.
 */
struct SysStorage
{
	bool valid;
	long long total_bytes;
	long long available_bytes;
};

SysStorage get_system_storage(const std::string& path);

/** What THIS PROCESS is using, in bytes. Not the machine - this process.
 *
 * `current_bytes` is the resident set size: physical memory held right now.
 *
 * `peak_bytes` is the high-water mark since the process started, and it is the one that matters
 * here. pest++ peaks INSIDE the threaded upgrade solve, where every thread materialises dense
 * matrices for its own case - not when the ensembles load. Sampling current usage between
 * iterations therefore reports a comfortable number right up until the run dies. The kernel
 * records the peak whether or not anyone was watching.
 *
 * Separate from SysMemory on purpose: on a host running eight agents the machine figures are
 * shared between them and these are not.
 *
 * `valid` is false when nothing could be read. The numbers are zero then and mean nothing -
 * check it rather than reading zero as "no memory in use".
 */
struct SysProcessMemory
{
	bool valid;
	long long current_bytes;
	long long peak_bytes;
};

SysProcessMemory get_process_memory();

#ifdef OS_WIN
#include <Windows.h>
PROCESS_INFORMATION start(std::string &cmd_string);
#endif
#ifdef OS_LINUX
int start(std::string &cmd_string);
#endif




#endif /* SYSTEM_VARIABLES_H_ */

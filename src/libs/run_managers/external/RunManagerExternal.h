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
#ifndef _RUNMANAGEEXTERNAL_H_
#define _RUNMANAGEEXTERNAL_H_

#include "RunManagerAbstract.h"
#include <string>
#include <vector>


class RunManagerExternal :
	public RunManagerAbstract
{
public:
	RunManagerExternal(std::vector<std::string> _comline_vec,
		std::vector<std::string> _tplfile_vec, std::vector<std::string> _inpfile_vec,
		std::vector<std::string> _insfile_vec, std::vector<std::string> _outfile_vec,
		std::string &stor_filename,
		int _max_n_failure = 1,int _sleep_ms=10);
	virtual ~RunManagerExternal();
	virtual void run();
private:
	int sleep_ms;
};

#endif /* _RUNMANAGEEXTERNAL_H_ */

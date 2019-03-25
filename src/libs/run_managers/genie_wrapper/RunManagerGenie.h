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
#ifndef RUNMANAGERGENIE_H
#define RUNMANAGERGENIE_H

#include<string>
#include <vector>
#include <set>
#include "RunManagerAbstract.h"

class Parameters;
class Observations;
class ModelRun;
class ModelExecInfo;

class RunManagerGenie : public RunManagerAbstract
{
public:
	static const int LEN_PARAMETER_NAME;
	static const int LEN_OBSERVATION_NAME;
	RunManagerGenie(const std::vector<std::string> _comline_vec,
		const std::vector<std::string> _tplfile_vec, const std::vector<std::string> _inpfile_vec,
		const std::vector<std::string> _insfile_vec, const std::vector<std::string> _outfile_vec,
		const std::string &stor_filename, const std::string &_host, const std::string &_id="PPEST");
	virtual void run();
	virtual ~RunManagerGenie(void);
protected:
	std::string id;
	std::string host;
};

#endif /* RUNMANAGERGENIE_H */

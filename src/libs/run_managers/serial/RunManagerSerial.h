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
#ifndef RUNMANAGERSERIAL_H
#define RUNMANAGERSERIAL_H

#include "RunManagerAbstract.h"
#include <string>
#include "model_interface.h"

class RunManagerSerial : public RunManagerAbstract
{
public:
	RunManagerSerial(const std::vector<std::string> _comline_vec,
		const std::vector<std::string> _tplfile_vec, const std::vector<std::string> _inpfile_vec,
		const std::vector<std::string> _insfile_vec, const std::vector<std::string> _outfile_vec,
		const std::string &stor_filename, const std::string &run_dir, int _max_run_fail=1,
		bool fill_tpl_zeros=false, string additional_ins_delimiters="", int _num_threads=1,
		bool tpl_force_decimal=false, bool should_echo=true);
	virtual void run();
	~RunManagerSerial(void);
private:
	ModelInterface mi;
	std::string run_dir;


    void run_async(pest_utils::thread_flag* terminate, pest_utils::thread_flag* finished,
                   exception_ptr& run_exception,
                   Parameters* pars, Observations* obs);
    void run(Parameters* pars, Observations* obs);
};

#endif /* RUNMANAGERSERIAL_H */

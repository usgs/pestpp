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
	// Where the batch has got to, so the inner poll loop can report during a run that is
	// still executing. Without this the observer would only hear between runs, and a slow
	// model - the case that most needs a progress bar - would look frozen for its whole
	// duration.
	int prog_nruns = 0, prog_done = 0, prog_failed = 0, prog_cur_run = -1;
	std::chrono::system_clock::time_point prog_start;
	RunProgress make_progress() const
	{
		RunProgress p;
		p.n_total = prog_nruns;
		p.n_completed = prog_done;
		p.n_failed = prog_failed;
		p.n_running = (prog_cur_run >= 0) ? 1 : 0;   // serial runs one at a time
		p.run_id = prog_cur_run;
		std::chrono::duration<double> d = std::chrono::system_clock::now() - prog_start;
		p.elapsed_sec = d.count();
		return p;
	}

	ModelInterface mi;
	std::string run_dir;


    void run_async(pest_utils::thread_flag* terminate, pest_utils::thread_flag* finished,
                   exception_ptr& run_exception,
                   Parameters* pars, Observations* obs);
    void run(Parameters* pars, Observations* obs);
};

#endif /* RUNMANAGERSERIAL_H */

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

#ifndef PANTHERAGENT_H_
#define PANTHERAGENT_H_

#include "network_wrapper.h"
#include <iostream>
#include <fstream>
#include <iostream>
#include <fstream>
#include <memory>
#include "utilities.h"
#include "pest_error.h"
#include "network_package.h"
#include "Transformable.h"
#include "model_interface.h"

class PANTHERAgent{
public:
	PANTHERAgent(ofstream &_frec);
	void init_network(const std::string &host, const std::string &port);
	void start(const std::string &host, const std::string &port);
	~PANTHERAgent();

	pair<int,string> recv_message(NetPackage &net_pack, struct timeval *tv=NULL);
	pair<int,string> recv_message(NetPackage &net_pack, long  timeout_seconds, long  timeout_microsecs = 0);
	pair<int,string> send_message(NetPackage &net_pack, const void *data=NULL, unsigned long data_len=0);
	std::pair<NetPackage::PackType,std::string> run_model(Parameters &pars, Observations &obs, NetPackage &net_pack);
	void process_ctl_file(const string &ctl_filename);
private:
	ofstream& frec;
	int sockfd;
	int fdmax;
	double run_time;
	int poll_interval_seconds;
	int max_time_without_master_ping_seconds;
	bool restart_on_error;
	int current_da_cycle;

#ifdef _DEBUG
	static const int max_recv_fails = 100;
	static const int max_send_fails = 100;
#else
	static const int max_recv_fails = 1000;
	static const int max_send_fails = 1000;
#endif
	static const int recv_timeout_secs = 10;
	bool terminate;
	fd_set master;
	/*std::vector<std::string> comline_vec;
	std::vector<std::string> tplfile_vec;
	std::vector<std::string> inpfile_vec;
	std::vector<std::string> insfile_vec;
	std::vector<std::string> outfile_vec;
	std::vector<std::string> obs_name_vec;
	std::vector<std::string> par_name_vec;*/

	void start_impl(const std::string &host, const std::string &port);

	ModelInterface mi;
	void run_async(pest_utils::thread_flag* terminate, pest_utils::thread_flag* finished,
		exception_ptr& run_exception,
		Parameters* pars, Observations* obs);

	void terminate_or_restart(int error_code) const;

	//Observations ctl_obs;
	//Parameters ctl_pars;
	Pest pest_scenario;

	void report(const string& _message, bool to_cout);

	void transfer_files(const vector<string>& tfiles, int group, int run_id, string& desc, string tag);

};


class PANTHERAgentRestartError final : public std::runtime_error {
public:
	PANTHERAgentRestartError(const string &_message="") : runtime_error(_message), message(_message) {}
	virtual ~PANTHERAgentRestartError() throw () {}
	void add_front(const string &s) {message = s + message;}
	void add_back(const string &s) {message += s;}
	void raise() {throw *this;}
	virtual const char* what() const throw()
	{
		return message.c_str();
	}
protected:
	string message;
};


#endif /* PANTHERAGENT_H_ */

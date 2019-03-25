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

#ifndef PANTHERSLAVE_H_
#define PANTHERSLAVE_H_

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

class PANTHERSlave{
public:
	PANTHERSlave();
	void init_network(const std::string &host, const std::string &port);
	void start(const std::string &host, const std::string &port);
	~PANTHERSlave();
	void run();
	int recv_message(NetPackage &net_pack, struct timeval *tv=NULL);
	int recv_message(NetPackage &net_pack, long  timeout_seconds, long  timeout_microsecs = 0);
	int send_message(NetPackage &net_pack, const void *data=NULL, unsigned long data_len=0);
	NetPackage::PackType run_model(Parameters &pars, Observations &obs, NetPackage &net_pack);
	//int run_model(Parameters &pars, Observations &obs);
	std::string tpl_err_msg(int i);
	std::string ins_err_msg(int i);
	void check_io();
	//void listener(pest_utils::thread_flag* terminate, pest_utils::thread_flag* finished);
	void listener();
	void process_ctl_file(const string &ctl_filename);
	void process_panther_ctl_file(const string &ctl_filename);
private:
	int sockfd;
	int fdmax;
	double run_time;
	int poll_interval_seconds;
#ifdef _DEBUG
	static const int max_recv_fails = 100;
	static const int max_send_fails = 100;
#else
	static const int max_recv_fails = 1000;
	static const int max_send_fails = 1000;
#endif
	static const int recv_timeout_secs = 1;
	bool terminate;
	fd_set master;
	std::vector<std::string> comline_vec;
	std::vector<std::string> tplfile_vec;
	std::vector<std::string> inpfile_vec;
	std::vector<std::string> insfile_vec;
	std::vector<std::string> outfile_vec;
	std::vector<std::string> obs_name_vec;
	std::vector<std::string> par_name_vec;

	ModelInterface mi;
	void run_async(pest_utils::thread_flag* terminate, pest_utils::thread_flag* finished,
		pest_utils::thread_exceptions *shared_execptions,
		Parameters* pars, Observations* obs);

};

#endif /* PANTHERSLAVE_H_ */

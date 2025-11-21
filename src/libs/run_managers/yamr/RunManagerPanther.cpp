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

#include "network_wrapper.h"
#include "RunManagerPanther.h"
#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <list>
#include <iterator>
#include <cassert>
#include <cstring>
#include <map>
#include <deque>
#include <utility>
#include <algorithm>
#include "network_wrapper.h"
#include "network_package.h"
#include "Transformable.h"
#include "utilities.h"
#include "Serialization.h"


using namespace std;
using namespace pest_utils;

const int RunManagerPanther::BACKLOG = 1000;
const int RunManagerPanther::MAX_FAILED_PINGS = 60;
const int RunManagerPanther::N_PINGS_UNRESPONSIVE = 3;
const int RunManagerPanther::MIN_PING_INTERVAL_SECS = 60;				// Ping each slave at most once every minute
const int RunManagerPanther::MAX_PING_INTERVAL_SECS = 120;				// Ping each slave at least once every 2 minutes
const int RunManagerPanther::MAX_CONCURRENT_RUNS_LOWER_LIMIT = 1;
const int RunManagerPanther::IDLE_THREAD_SIGNAL_TIMEOUT_SECS = 10;  // Allow up to 10s for the run_idle_async() thread to acknowledge signals (pause idling, terminate)
const double RunManagerPanther::MIN_AVGRUNMINS_FOR_KILL = 0.01; //minimum avg runtime to try to kill and/or resched runs
//const int RunManagerPanther::MILLISECONDS_BETWEEN_ECHOS = 10;
//const int RunManagerPanther::TIMEOUT_MILLISECONDS = 10;

AgentInfoRec::AgentInfoRec(int _socket_fd)
{
	socket_fd = _socket_fd;
	name_info_vec = w_getnameinfo_vec(_socket_fd);
	run_id = UNKNOWN_ID;
	group_id = UNKNOWN_ID;
	state = AgentInfoRec::State::NEW;
	work_dir = "";
	linpack_time = std::chrono::hours(-500);
	run_time = std::chrono::hours(-500);
	start_time = std::chrono::system_clock::now();
	last_ping_time = std::chrono::system_clock::now();
	ping = false;
	failed_pings = 0;
	failed_runs = 0;
	state_strings = vector<string>({ "NEW", "CWD_REQ", "CWD_RCV", "NAMES_SENT", "LINPACK_REQ", "LINPACK_RCV", "WAITING", "ACTIVE",
	"KILLED", "KILLED_FAILED", "COMPLETE" });
	state = State::NEW;

}

bool AgentInfoRec::CompareTimes::operator() (const AgentInfoRec &a, const AgentInfoRec &b)
{
	bool ret = false;
	if (a.run_time > std::chrono::milliseconds(0) && b.run_time > std::chrono::milliseconds(0))
	{
		ret = (a.run_time < b.run_time);
	}
	else if (a.linpack_time > std::chrono::milliseconds(0) && b.linpack_time > std::chrono::milliseconds(0))
	{
		ret = (a.linpack_time < b.linpack_time);
	}
	return ret;
}


int AgentInfoRec::get_socket_fd() const
{
	return socket_fd;
}

void AgentInfoRec::set_socket_fd(int _socket_fd)
{
	socket_fd = _socket_fd;
}

string AgentInfoRec::get_hostname()const
{
	return name_info_vec[0];
}

string AgentInfoRec::get_port()const
{
	return name_info_vec[1];
}

string AgentInfoRec::get_socket_name()const
{
	return name_info_vec[0] + ":" + name_info_vec[1];
}

int AgentInfoRec::get_run_id() const
{
	return run_id;
}

void AgentInfoRec::set_run_id(int _run_id)
{
	run_id = _run_id;
}

int AgentInfoRec::get_group_id() const
{
	return group_id;
}

void AgentInfoRec::set_group_id(int _group_id)
{
	group_id = _group_id;
}

AgentInfoRec::State AgentInfoRec::get_state() const
{
	return state;
}

void AgentInfoRec::set_state(const State &_state)
{
	if (_state == AgentInfoRec::State::ACTIVE)
	{
		throw PestError("AgentInfo::set_state: run_id and group_id must be supplied when state it set to active");
	}
	state = _state;
}

void AgentInfoRec::set_state(const State &_state, int _run_id, int _group_id)
{
	state = _state;
	run_id = _run_id;
	group_id = _group_id;
}

void AgentInfoRec::set_work_dir(const std::string &_work_dir)
{
	work_dir = _work_dir;
}

string AgentInfoRec::get_work_dir() const
{
	return work_dir;
}

void AgentInfoRec::start_timer()
{
	start_time = std::chrono::system_clock::now();
}

void AgentInfoRec::end_run()
{
	auto dt = std::chrono::system_clock::now() - start_time;
	if (run_time > std::chrono::hours(0))
	{
		run_time = run_time + dt;
		run_time /= 2;
	}
	else
	{
		run_time = dt;
	}
}

void AgentInfoRec::end_linpack()
{
	linpack_time = std::chrono::system_clock::now() - start_time;
}

double AgentInfoRec::get_duration_sec() const
{
	chrono::system_clock::duration dt = chrono::system_clock::now() - start_time;
	return (double)std::chrono::duration_cast<std::chrono::milliseconds>(dt).count() / 1000.0;
}

double AgentInfoRec::get_duration_minute() const
{
	return get_duration_sec() / 60.0;
}

double AgentInfoRec::get_runtime_sec() const
{
	return(double)std::chrono::duration_cast<std::chrono::milliseconds>(run_time).count() / 1000.0;
}

double AgentInfoRec::get_runtime_minute() const
{
	double run_minutes = std::chrono::duration_cast<std::chrono::milliseconds>(run_time).count() / 60000.0;
	return run_minutes;
}

double AgentInfoRec::get_runtime() const
{
	return double(run_time.count());
}

double AgentInfoRec::get_linpack_time() const
{
	return double(linpack_time.count());
}


void AgentInfoRec::reset_failed_pings()
{
	failed_pings = 0;
}

int AgentInfoRec::add_failed_ping()
{
	failed_pings++;
	return failed_pings;
}

int AgentInfoRec::add_failed_run()
{
	failed_runs++;
	return failed_runs;
}

void AgentInfoRec::set_ping(bool val)
{
	ping = val;
	//a success response
	if (!val) reset_failed_pings();
	//sending a request
	else reset_last_ping_time();
}

bool AgentInfoRec::get_ping() const
{
	return ping;
}

int AgentInfoRec::get_failed_pings() const
{
	return failed_pings;
}

void AgentInfoRec::reset_last_ping_time()
{
	last_ping_time = chrono::system_clock::now();
}

int AgentInfoRec::seconds_since_last_ping_time() const
{
	return chrono::duration_cast<std::chrono::seconds>
		(chrono::system_clock::now() - last_ping_time).count();
}


RunManagerPanther::RunManagerPanther(const string& stor_filename, const string& _port, ofstream& _f_rmr, int _max_n_failure,
	double _overdue_reched_fac, double _overdue_giveup_fac, double _overdue_giveup_minutes, bool _should_echo, const vector<string>& par_names,
	const vector<string>& obs_names,int _timeout_milliseconds,int _echo_interval_milliseconds, bool _persistent_workers)

	: RunManagerAbstract(vector<string>(), vector<string>(), vector<string>(),
		vector<string>(), vector<string>(), stor_filename, _max_n_failure),
	overdue_reched_fac(_overdue_reched_fac), overdue_giveup_fac(_overdue_giveup_fac),
	port(_port), f_rmr(_f_rmr), n_no_ops(0), overdue_giveup_minutes(_overdue_giveup_minutes),
	terminate_idle_thread(false), currently_idle(true), idling(false), idle_thread_finished(false),
	idle_thread(nullptr), should_echo(_should_echo),nftx(0),timeout_milliseconds(_timeout_milliseconds),
    echo_interval_milliseconds(_echo_interval_milliseconds),persistent_workers(_persistent_workers)
{

	const char * t =
    "\n\n:~-._                                                 _.-~:\n"
    ": :.~^o._        ________---------________        _.o^~.:.:\n"
    " : ::.`?88booo~~~.::::::::...::::::::::::..~~oood88P'.::.:\n"
    " :  ::: `?88P .:::....         ........:::::. ?88P' :::. :\n"
    "  :  :::. `? .::.            . ...........:::. P' .:::. :\n"
    "   :  :::   ... ..  ...       .. .::::......::.   :::. :\n"
    "   `  :' .... ..  .:::::.     . ..:::::::....:::.  `: .'\n"
    "    :..    ____:::::::::.  . . ....:::::::::____  ... :\n"
    "   :... `:~    ^~-:::::..  .........:::::-~^    ~::.::::\n"
    "   `.::. `\\   (8)  \\b:::..::.:.:::::::d/  (8)   /'.::::'\n"
    "    ::::.  ~-._v    |b.::::::::::::::d|    v_.-~..:::::\n"
    "    `.:::::... ~~^?888b..:::::::::::d888P^~...::::::::'\n"
    "     `.::::::::::....~~~ .:::::::::~~~:::::::::::::::'\n"
    "      `..:::::::::::   .   ....::::    ::::::::::::,'\n"
    "        `. .:::::::    .      .::::.    ::::::::'.'\n"
    "          `._ .:::    .        :::::.    :::::_.'\n"
    "             `-. :    .        :::::      :,-'\n"
    "                :.   :___     .:::___   .::\n"
    "      ..--~~~~--:+::. ~~^?b..:::dP^~~.::++:--~~~~--..\n"
    "        ___....--`+:::.    `~8~'    .:::+'--....___\n"
    "      ~~   __..---`_=:: ___gd8bg___ :==_'---..__   ~~\n"
    "       -~~~  _.--~~`-.~~~~~~~~~~~~~~~,-' ~~--._ ~~~-\n\n";
             



	cout << t << endl;
	cout << "               starting PANTHER master..." << endl << endl;
	max_concurrent_runs = max(MAX_CONCURRENT_RUNS_LOWER_LIMIT, _max_n_failure);
	w_init();
	std::pair<int, string> status;
	struct addrinfo hints;
	struct addrinfo* servinfo;
	memset(&hints, 0, sizeof hints);
	//Use this for IPv4 and IPv6
	//hints.ai_family = AF_UNSPEC;
	//Use this just for IPv4;
	hints.ai_family = AF_INET;
	hints.ai_socktype = SOCK_STREAM;
	hints.ai_flags = AI_PASSIVE;

	status = w_getaddrinfo(NULL, port.c_str(), &hints, &servinfo);
	if (status.first != 0)
	{
		cout << "ERROR: getaddrinfo returned non-zero: " << status.second << endl;
		throw(PestError("ERROR: getaddrinfo returned non-zero: " + status.second));
	}

	w_print_servinfo(servinfo, cout);
	cout << endl;
	//make socket, bind and listen
	addrinfo* connect_addr = w_bind_first_avl(servinfo, listener);
	if (connect_addr == nullptr)
	{
		stringstream err_str;
		err_str << "Error: port \"" << port << "\n is busy.  Can not bind port" << endl;
		throw(PestError(err_str.str()));
	}
	else {

	}
	w_listen(listener, BACKLOG);
	f_rmr << endl;
	cout << "PANTHER master listening on socket: " << w_get_addrinfo_string(connect_addr) << endl;
	f_rmr << "PANTHER master listening on socket:" << w_get_addrinfo_string(connect_addr) << endl;
	//free servinfo
	freeaddrinfo(servinfo);
	fdmax = listener;
	FD_ZERO(&master);
	FD_SET(listener, &master);
	par_names_to_check_worker = par_names;
	obs_names_to_check_worker = obs_names;
	mgr_type = RUN_MGR_TYPE::PANTHER;
}

int RunManagerPanther::get_n_concurrent(int run_id)
{
	auto range_pair = active_runid_to_iterset_map.equal_range(run_id);
	int n = 0;
	for (auto &i = range_pair.first; i != range_pair.second; ++i)
	{
		if (i->second->get_state() == AgentInfoRec::State::ACTIVE)
		{
			++n;
		}
	}
	return n;
}

list<AgentInfoRec>::iterator RunManagerPanther::get_active_run_iter(int socket)
{
	auto iter = socket_to_iter_map.find(socket);

	if (iter != socket_to_iter_map.end())
	{
		return socket_to_iter_map.find(socket)->second;
	}
	else
	{
		return agent_info_set.end();
	}
}


void RunManagerPanther::initialize(const Parameters &model_pars, const Observations &obs, const string &_filename)
{
	RunManagerAbstract::initialize(model_pars, obs, _filename);
	cur_group_id = NetPackage::get_new_group_id();
	
}

void RunManagerPanther::initialize_restart(const std::string &_filename)
{
	file_stor.init_restart(_filename);
	free_memory();
	vector<int> waiting_run_id_vec = get_outstanding_run_ids();
	for (int &id : waiting_run_id_vec)
	{
		waiting_runs.push_back(id);
	}
}

void RunManagerPanther::reinitialize(const std::string &_filename)
{
	free_memory();
	RunManagerAbstract::reinitialize(_filename);
	cur_group_id = NetPackage::get_new_group_id();
}

void  RunManagerPanther::free_memory()
{
	waiting_runs.clear();
	model_runs_done = 0;
	failure_map.clear();
	active_runid_to_iterset_map.clear();
}

int RunManagerPanther::add_run(const Parameters &model_pars, const string &info_txt, double info_value)
{
	int run_id = file_stor.add_run(model_pars, info_txt, info_value);
	waiting_runs.push_back(run_id);
	return run_id;
}

int RunManagerPanther::add_run(const std::vector<double> &model_pars, const string &info_txt, double info_value)
{
	int run_id = file_stor.add_run(model_pars, info_txt, info_value);
	waiting_runs.push_back(run_id);
	return run_id;
}

int RunManagerPanther::add_run(const Eigen::VectorXd &model_pars, const string &info_txt, double info_value)
{
	int run_id = file_stor.add_run(model_pars, info_txt, info_value);
	waiting_runs.push_back(run_id);
	return run_id;
}

void RunManagerPanther::update_run(int run_id, const Parameters &pars, const Observations &obs)
{

	file_stor.update_run(run_id, pars, obs);
	// erase any waiting runs with this id
	for (auto it_run = waiting_runs.begin(); it_run != waiting_runs.end();)
	{
		if (*it_run == run_id)
		{
			it_run = waiting_runs.erase(it_run);
		}
		else
		{
			++it_run;
		}
	}
	// kill any active runs with this id
	kill_runs(run_id, false, "run not required");
}

void RunManagerPanther::run()
{
	run_until(RUN_UNTIL_COND::NORMAL);
}

RunManagerAbstract::RUN_UNTIL_COND RunManagerPanther::run_until(RUN_UNTIL_COND condition, int max_no_ops, double max_time_sec)
{
	RUN_UNTIL_COND terminate_reason = RUN_UNTIL_COND::NORMAL;
	stringstream message;
	NetPackage net_pack;

	// Pause idle pinging thread
	pause_idle();
	stringstream ss;
	model_runs_done = 0;
	model_runs_failed = 0;
	model_runs_timed_out = 0;
	bytes_transferred = 0;
	files_transferred = 0;
	failure_map.clear();
	active_runid_to_iterset_map.clear();
	open_file_trans_streams.clear();
	open_file_socket_map.clear();
	int num_runs = waiting_runs.size();
	cout << "    running model " << num_runs << " times" << endl;
	f_rmr << "running model " << num_runs << " times" << endl;
	cout << "    starting at " << pest_utils::get_time_string() << endl;
	if (agent_info_set.size() == 0) // first entry is the listener, slave appears after this
	{
		cout << endl << "    waiting for agents to appear..." << endl << endl;
		//f_rmr << endl << "    waiting for agents to appear..." << endl << endl;
		report("waiting for agents to appear", false);
	}
	else
	{
		for (auto &si : agent_info_set)
			si.reset_runtime();
		ss.str("");
		ss << "    " << agent_info_set.size() << " agents ready";
		report(ss.str(), false);
		cout << ss.str() << endl << endl;
		
	}
	cout << endl;
	f_rmr << endl;
	if (should_echo)
	{
		cout << "PANTHER progress" << endl;
		cout << "   avg = average model run time in minutes" << endl;
		cout << "   runs(C = completed | F = failed | T = timed out)" << endl;
		cout << "   agents(R = running | W = waiting | U = unavailable)" << endl;
		cout << "--------------------------------------------------------------------------------" << endl;
	}
	else
	{
		cout << "'panther_echo' is 'false', running in silent mode - see rmr file for details" << endl;
	}

	std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();
    last_echo_time = std::chrono::system_clock::now();

	while (!all_runs_complete() && terminate_reason == RUN_UNTIL_COND::NORMAL)
	{
        int q = pest_utils::quit_file_found();
        if ((q == 1) || (q == 2) || (q == 4))
        {
            cout << "'pest.stp' found" << endl;
            kill_all_active_runs();

        }
		echo();
		init_agents();
		//schedule runs on available nodes
		schedule_runs();
		echo();
		// get and process incoming messages
		if (!listen())
		{
			++n_no_ops;
		}
		else
		{
			n_no_ops = 0;
		}
        echo();
		if (ping())
		{
			n_no_ops = 0;
		}
        echo();
		if ((condition == RUN_UNTIL_COND::NO_OPS || condition == RUN_UNTIL_COND::NO_OPS_OR_TIME) && n_no_ops >= max_no_ops)
		{
			terminate_reason = RUN_UNTIL_COND::NO_OPS;
		}

		if ((condition == RUN_UNTIL_COND::TIME || condition == RUN_UNTIL_COND::NO_OPS_OR_TIME) && get_duration_sec(start_time) >= max_time_sec)
		{
			terminate_reason = RUN_UNTIL_COND::TIME;
		}

	}
    w_sleep(get_current_sleep_timeout_milliseconds(timeout_milliseconds));
	n_no_ops = 0;
    while (true)
    {
        if (!listen())
        {
            ++n_no_ops;
        }
        else
        {
            n_no_ops = 0;
        }

        cout << get_time_string_short() << " remaining file transfers: " << open_file_trans_streams.size() << "\r" << flush;
        if (ping())
        {
            n_no_ops = 0;
        }
        if (agent_info_set.size() == 0)
        {
            ss.str("");
            ss << "lost comms with all agents, closing all remaining open file transfers";
            report(ss.str(),true);
            for (auto& m : open_file_socket_map) {
                string fname = m.second;
                pair<map<string, ofstream *>::iterator, bool> ret = open_file_trans_streams.insert(
                        pair<string, ofstream *>(fname, new ofstream));
                ofstream &out = *ret.first->second;
                int file_size = out.tellp();
                out.flush();
                out.close();
                open_file_trans_streams.erase(ret.first);
                ss.str("");
                ss << "closed master_file:" << fname << " for file transfer, transferred bytes:" << file_size;
                report(ss.str(), false);
                files_transferred += 1;

            }
            open_file_socket_map.clear();

        }
        if (open_file_trans_streams.size() == 0)
        {
            if (!listen())
                break;
        }
        int q = pest_utils::quit_file_found();
        if ((q == 1) || (q == 2) || (q == 4))
        {
            cout << "'pest.stp' found" << endl;
            kill_all_active_runs();
            ss.str("");
            ss << "'pest.stp' found, closing all remaining open file transfers";
            report(ss.str(),true);
            for (auto& m : open_file_socket_map) {
                string fname = m.second;
                pair<map<string, ofstream *>::iterator, bool> ret = open_file_trans_streams.insert(
                        pair<string, ofstream *>(fname, new ofstream));
                ofstream &out = *ret.first->second;
                int file_size = out.tellp();
                out.flush();
                out.close();
                open_file_trans_streams.erase(ret.first);
                ss.str("");
                ss << "closed master_file:" << fname << " for file transfer, transferred bytes:" << file_size;
                report(ss.str(), false);
                files_transferred += 1;

            }

        }


    }

    if (terminate_reason == RUN_UNTIL_COND::NORMAL)
	{
		echo();
		total_runs += model_runs_done;
		//kill any remaining active runs
		kill_all_active_runs();
		echo();
		cout << endl << endl;
		message.str("");
		message << "   " << model_runs_done << " runs complete :  " << get_num_failed_runs() << " runs failed" << endl;
		double duration = pest_utils::get_duration_sec(start_time) / 60.0;
		message << "   " << setprecision(3) << get_global_runtime_minute() << " avg run time (min) : ";
		message << setprecision(3) << duration << " run mgr time (min)" << endl;
		message << "   " << agent_info_set.size() << " agents connected" << endl;
		if (bytes_transferred > 0)
        {

            message << "   " << files_transferred << " files transferred" << endl;
            message << "   " << bytes_transferred << " bytes transferred" << endl;

        }
		
		cout << message.str() << endl << endl;
		f_rmr << endl << "---------------------" << endl << message.str() << endl;

		//Removed because it was preventing the restart from functioning properly
		//if (model_runs_done == 0)
		//	throw PestError("no runs completed successfully");

		set<int> fids = get_failed_run_ids();
		if (fids.size() > 0)
		{
			f_rmr << "  failed run_ids and (attempts):";
			for (auto fid : fids)
				f_rmr << " " << fid << "(" << failure_map.count(fid) << ")";
		}
		f_rmr << endl << endl;
			

		if (init_sim.size() == 0)
		{
			vector<double> pars;
			int status = file_stor.get_run(0, pars, init_sim);
		}
	}





    // Resume idle pinging thread

    if (!persistent_workers) {
        stringstream ss;

        std::list<list<AgentInfoRec>::iterator> free_agent_list = get_free_agent_list();
        list<AgentInfoRec>::iterator it_agent, iter_e;
        for (int i=0;i<free_agent_list.size();i++) {
            for (it_agent = agent_info_set.begin(), iter_e = agent_info_set.end();
                 it_agent != iter_e; ++it_agent) {
                AgentInfoRec::State state = it_agent->get_state();
                if (state == AgentInfoRec::State::WAITING) {
                    ss.str("");
                    ss << "using non-persistent agents, closed connection to agent: " << it_agent->get_socket_name()
                       << ", number of agents: " << socket_to_iter_map.size();
                    close_agent(it_agent);


                    report(ss.str(), false);
                    break;
                }
            }
        }
    }
    resume_idle();

	return terminate_reason;
}

bool RunManagerPanther::ping(pest_utils::thread_flag* terminate/* = nullptr*/)
{
	bool ping_sent = false;
	for (auto &i : socket_to_iter_map)
	{
		if (terminate && terminate->get())
		{
			break;
		}

		if (ping(i.first))
		{
			ping_sent = true;
	}
	}
	return ping_sent;
}

void RunManagerPanther::run_idle_async()
{
	try
	{
		// Initialise idle thread state flags
		idle_thread_finished.set(false);
		idling.set(false);

		// Continue to accept new connections and ping agents while the manager is alive and not processing any runs
		while (!terminate_idle_thread.get())
		{
			// Don't do anything unless the run manager is currently idle (in between calls to run())
			if(!currently_idle.get())
			{
				idling.set(false);

				// Sleep 1s to avoid spinlock
				w_sleep(get_current_sleep_timeout_milliseconds(timeout_milliseconds));
				continue;
			}

			// Confirm that the thread has started idling
			idling.set(true);

			// Initialise any new agents
			init_agents(&terminate_idle_thread);

			if(terminate_idle_thread.get())
			{
				break;
			}

			// Receive any waiting data from agents
			listen(&terminate_idle_thread);

			if(terminate_idle_thread.get())
			{
				break;
			}

			// Ping agents periodically to ensure connections are not dropped
			ping(&terminate_idle_thread);
		}
	}
	catch(...)
	{
		try
		{
			idle_thread_finished.set(true);
			idling.set(false);
		}
		catch(...)
		{
			// Don't hide original exception if something goes wrong with flag lock and set
		}

		throw;
	}

	idle_thread_finished.set(true);
	idling.set(false);
}

void RunManagerPanther::start_run_idle_async()
{
	if(idle_thread)
	{
		return;
	}

	// Initialise idle thread flags
	idling.set(false);
	idle_thread_finished.set(false);
	terminate_idle_thread.set(false);
	currently_idle.set(true);

	// Start thread
	idle_thread = new thread(&RunManagerPanther::run_idle_async, this);
	//idle_thread_raii = new pest_utils::thread_RAII(*idle_thread);

	report("Started idle ping thread.", false);
}

void RunManagerPanther::end_run_idle_async()
{
	if(!idle_thread)
	{
		return;
	}

	// Signal thread to terminate
	terminate_idle_thread.set(true);

	// Wait for termination
	bool timed_out = false;
	std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();
	while (!idle_thread_finished.get())
	{
		// Give up if we've spent too much time waiting already
		auto wait_time = std::chrono::system_clock::now() - start_time;
		if (wait_time > std::chrono::seconds(IDLE_THREAD_SIGNAL_TIMEOUT_SECS))
		{
			timed_out = true;
			break;
		}
		
		// Sleep to avoid spinlock
		w_sleep(get_current_sleep_timeout_milliseconds(timeout_milliseconds));
	}

	report("Stopped idle ping thread, as Panther manager is shutting down.", false);
	if(timed_out)
	{
		report("Warning: timed out waiting for acknowledgement of signal from idle thread.", false);
	}

    idle_thread->join();
	// Clean up
	//delete idle_thread_raii;
	//delete idle_thread;
	//idle_thread_raii = nullptr;
	idle_thread = nullptr;
}

void RunManagerPanther::pause_idle()
{
	if(!idle_thread)
	{
		return;
	}

	// Signal thread to stop idle pings
	currently_idle.set(false);

	// Wait for the signal to be acknowledged
	bool timed_out = false;
	auto start_time = std::chrono::system_clock::now();
	while (idling.get() && !idle_thread_finished.get())
	{
		// Give up if we've spent too much time waiting already
		auto waited = std::chrono::system_clock::now() - start_time;
		if (waited > std::chrono::seconds(IDLE_THREAD_SIGNAL_TIMEOUT_SECS))
		{
			timed_out = true;
			break;
		}
		
		// Sleep to avoid spinlock
		w_sleep(get_current_sleep_timeout_milliseconds(timeout_milliseconds));
	}

	report("Panther idle ping thread paused prior to scheduling runs.", false);
	if(timed_out)
	{
		report("Warning: timed out waiting for acknowledgement of signal from idle thread.", false);
	}
    //delete idle_thread;
}

void RunManagerPanther::resume_idle()
{
	// Start up the thread if it has not already been started
	start_run_idle_async();

	// Signal thread to continue idle pings
	currently_idle.set(true);

	// Don't bother waiting for acknowledgement here, as none of the management code relies on it; we can happily go off and do other processing while the thread gets around to resuming idle pings
	report("Panther idle ping thread resumed.", false);
}

int RunManagerPanther::get_current_sleep_timeout_milliseconds(const int org_timeout_milliseconds)
{
    double avg = get_global_runtime_minute() / 1000.0;
    if (org_timeout_milliseconds > 0)
    {
        return org_timeout_milliseconds;
    }
    double timeout;
    if (avg == 0.0) {
        timeout = 500;
    }
    else {
        timeout = avg / 0.1;
        timeout = max<double>(timeout, 10);
        timeout = min<double>(timeout, 500);
    }

    //cout << timeout;
    return timeout;
}

bool RunManagerPanther::ping(int i_sock)
{
	bool ping_sent = false;
	list<AgentInfoRec>::iterator agent_info_iter = socket_to_iter_map.at(i_sock);
	AgentInfoRec::State state = agent_info_iter->get_state();
	if (state != AgentInfoRec::State::WAITING
		&& state != AgentInfoRec::State::ACTIVE
		&& state != AgentInfoRec::State::KILLED
		&& state != AgentInfoRec::State::KILLED_FAILED)
	{
		return ping_sent;
	}

	string sock_hostname = agent_info_iter->get_hostname();
	fd_set read_fds = master;
	//if the agent hasn't communicated since the last ping request
	if ((!FD_ISSET(i_sock, &read_fds)) && agent_info_iter->get_ping())
	{
		int fails = agent_info_iter->add_failed_ping();
		report("failed to receive ping response from agent: " + sock_hostname + "$" + agent_info_iter->get_work_dir(), false);
		if (fails >= MAX_FAILED_PINGS)
		{
			ping_sent = true;
			report("max failed ping communications since last successful run form agent:" + sock_hostname + "$" + agent_info_iter->get_work_dir() + "  -> terminating", false);
			close_agent(i_sock);
			return ping_sent;
		}
	}
	//check if it is time to ping again...
	double duration = (double)agent_info_iter->seconds_since_last_ping_time();
	double ping_time = min(double(MAX_PING_INTERVAL_SECS), max(double(MIN_PING_INTERVAL_SECS), agent_info_iter->get_runtime_sec()));
	if (duration >= ping_time)
	{
		ping_sent = true;
		const char* data = "\0";

		NetPackage net_pack(NetPackage::PackType::PING, 0, 0, "");
		pair<int,string> err = net_pack.send(i_sock, data, 0);
		
		if (err.first <= 0)
		{
			int fails = agent_info_iter->add_failed_ping();
			report("failed to send ping request to agent:" + sock_hostname + "$" + agent_info_iter->get_work_dir() + ": " + err.second, false);
			if (fails >= MAX_FAILED_PINGS)
			{
				report("max failed ping communications since last successful run for agent:" + sock_hostname + "$" + agent_info_iter->get_work_dir() + "  -> terminating", true);
				close_agent(i_sock);
				return ping_sent;
			}
		}
		else agent_info_iter->set_ping(true);
			//report("ping sent to agent:" + sock_hostname + "$" + agent_info_iter->get_work_dir(), false);

#ifdef _DEBUG
		report("ping sent to agent:" + sock_hostname + "$" + agent_info_iter->get_work_dir(), false);
#endif
	}
	return ping_sent;
}


bool RunManagerPanther::listen(pest_utils::thread_flag* terminate/* = nullptr*/)
{
	bool got_message = false;
	struct sockaddr_storage remote_addr;
	fd_set read_fds; // temp file descriptor list for select()
	socklen_t addr_len;
	timeval tv;
	tv.tv_sec = 0;
	tv.tv_usec = 250;
	read_fds = master; // copy it
	if (w_select(fdmax+1, &read_fds, NULL, NULL, &tv) == -1)
	{
		// there are no slaves available.  W need to keep listening until at least one appears
		got_message = true;
		return got_message;
	}
	// run through the existing connections looking for data to read
	for (int i = 0; i <= fdmax; i++)
	{
		// Stop early if we're requested to terminate
	 	if(terminate && terminate->get())
	 	{
	 		break;
	 	}

		if (FD_ISSET(i, &read_fds)) { // we got one!!
			got_message = true;
			if (i == listener)  // handle new connections
			{
				int newfd;
				addr_len = sizeof remote_addr;
				newfd = w_accept(listener,(struct sockaddr *)&remote_addr, &addr_len);
				if (newfd == -1) {}
				else
				{
					add_agent(newfd);
				}
			}
			else  // handle data from a client
			{
				//set the ping flag since the slave sent something back
				list<AgentInfoRec>::iterator iter = socket_to_iter_map.at(i);
				iter->set_ping(false);
				process_message(i);
			} // END handle data from client
		} // END got new incoming connection
	} // END looping through file descriptors
	return got_message;
}

void RunManagerPanther::close_agents()
{
	/*for (int i = 0; i <= fdmax; i++)
	{
		list<SlaveInfoRec>::iterator slave_info_iter = socket_to_iter_map.at(i);
		if (slave_info_iter != slave_info_set.end())
			close_slave(slave_info_iter);
	}*/
	while (socket_to_iter_map.size() > 0)
	{
		listen();
		vector<int> sock_nums;
		for (auto &si : socket_to_iter_map)
			sock_nums.push_back(si.first);
		for (auto si : sock_nums)
			close_agent(si);

		w_sleep(get_current_sleep_timeout_milliseconds(timeout_milliseconds));

	}
}

void RunManagerPanther::close_agent(int i_sock)
{
	list<AgentInfoRec>::iterator agent_info_iter = socket_to_iter_map.at(i_sock);
	close_agent(agent_info_iter);
}

void RunManagerPanther::close_agent(list<AgentInfoRec>::iterator agent_info_iter)
{
	int i_sock = agent_info_iter->get_socket_fd();
	int run_id = agent_info_iter->get_run_id();
	AgentInfoRec::State state = agent_info_iter->get_state();

	string socket_name = agent_info_iter->get_socket_name();
	w_close(i_sock); // bye!
	FD_CLR(i_sock, &master); // remove from master set
	// remove run from active_runid_to_iterset_map
	unschedule_run(agent_info_iter);

	// check if this run needs to be returned to the waiting queue
	int n_concurr = get_n_concurrent(run_id);
	if (run_id != AgentInfoRec::UNKNOWN_ID &&  agent_info_iter->get_state() == AgentInfoRec::State::ACTIVE && n_concurr == 0)
	{
		waiting_runs.push_front(run_id);
	}

	agent_info_set.erase(agent_info_iter);
	socket_to_iter_map.erase(i_sock);
	if (open_file_socket_map.find(i_sock) != open_file_socket_map.end())
    {
	    string fname = open_file_socket_map.at(i_sock);
        pair<map<string ,ofstream*>::iterator, bool> ret = open_file_trans_streams.insert(pair<string,ofstream*>(fname,new ofstream));
        ofstream& out = *ret.first->second;
        int file_size = out.tellp();
        out.flush();
        out.close();
        open_file_trans_streams.erase(ret.first);
        stringstream ss;
        ss.str("");
        ss << "lost comms with agent, closed file:" << fname << " bytes:" << file_size << "  transferred";

        report(ss.str(),false);
        files_transferred += 1;
        open_file_socket_map.erase(i_sock);
    }


	stringstream ss;
	ss << "closed connection to agent: " << socket_name << ", number of agents: " << socket_to_iter_map.size();
	report(ss.str(), false);
}


void RunManagerPanther::schedule_runs()
{
	NetPackage net_pack;

	std::list<list<AgentInfoRec>::iterator> free_agent_list = get_free_agent_list();
	int n_responsive_agents = get_n_responsive_agents();
	//first try to schedule waiting runs
	for (auto it_run = waiting_runs.begin(); !free_agent_list.empty() && it_run != waiting_runs.end();)
	{
		int success = schedule_run(*it_run, free_agent_list, n_responsive_agents);
		if (success >= 0)
		{
			it_run = waiting_runs.erase(it_run);
		}
		else
		{
			++it_run;
		}
	}

	//check for overdue runs if there are no runs waiting to be processed
	if (n_no_ops > 0)
	{
		try
		{
			double duration, avg_runtime;
			double global_avg_runtime = get_global_runtime_minute();
			bool should_schedule = false;

			list<AgentInfoRec>::iterator it_agent, iter_e;
			for (it_agent = agent_info_set.begin(), iter_e = agent_info_set.end();
				it_agent != iter_e; ++it_agent)
			{
				AgentInfoRec::State state = it_agent->get_state();
				if (state == AgentInfoRec::State::ACTIVE)
				{
					should_schedule = false;
					int run_id = it_agent->get_run_id();
					int act_sock_id = it_agent->get_socket_fd();
					int n_concur = get_n_concurrent(run_id);

					duration = it_agent->get_duration_minute();
					avg_runtime = it_agent->get_runtime_minute();
					if (avg_runtime <= 0) avg_runtime = global_avg_runtime;
					if (avg_runtime <= 0) avg_runtime = 1.0E+300;
					vector<int> overdue_kill_runs_vec = get_overdue_runs_over_kill_threshold(run_id);

					if (failure_map.count(run_id) + overdue_kill_runs_vec.size() >= max_n_failure)
					{
						// kill the overdue runs
						//kill_runs(run_id, true, "overdue");
						stringstream ss;
						ss << "overdue. duration:" << duration << ", avg:" << avg_runtime;
						//kill_run(it_slave, ss.str());
						kill_runs(run_id, true, ss.str());
						should_schedule = false;
						//update_run_failed(run_id, it_slave->get_socket_fd());
						model_runs_timed_out += overdue_kill_runs_vec.size();
					}
					else if (overdue_kill_runs_vec.size() > max_concurrent_runs)
					{
						// kill the overdue runs
						kill_runs(run_id, true, "overdue");
						// reschedule runs as we still haven't reach the max failure threshold
						// and there are not concurrent runs for this id because we just killed all of them
						should_schedule = true;
						model_runs_timed_out += overdue_kill_runs_vec.size();
					}

					else if (((duration > overdue_giveup_minutes) || ((duration > avg_runtime*overdue_giveup_fac) &&
                                                                      (avg_runtime > MIN_AVGRUNMINS_FOR_KILL)))
						&& free_agent_list.empty())
					{
						// If there are no free slaves kill the overdue ones
						// This is necessary to keep runs with small numbers of slaves behaving
						stringstream ss;
						ss << "overdue. duration:" << duration << ", avg:" << avg_runtime;
						kill_run(it_agent, ss.str());
						update_run_failed(run_id, it_agent->get_socket_fd());

						if (failure_map.count(run_id) + overdue_kill_runs_vec.size() <= max_n_failure)
						{
							should_schedule = true;
						}
						model_runs_timed_out += 1;
					}

					else if ((duration > avg_runtime*overdue_reched_fac) && (avg_runtime > MIN_AVGRUNMINS_FOR_KILL))
					{
						//check how many concurrent runs are going
						if (n_concur < max_concurrent_runs) should_schedule = true;
						else should_schedule = false;
					}

					if ((!free_agent_list.empty()) && should_schedule)
					{
						string host_name = it_agent->get_hostname();
						stringstream ss;
						ss << "rescheduling overdue run " << run_id << " (" << duration << "|" <<
							avg_runtime << " minutes) on: " << host_name << "$" <<
							it_agent->get_work_dir();
						report(ss.str(), false);
						int success = schedule_run(run_id, free_agent_list, n_responsive_agents);
						n_concur = get_n_concurrent(run_id);
						if (success >= 0)
						{
							stringstream ss;
							ss << " run_id:" << run_id << " concurrent:" << n_concur;
							report(ss.str(), false);
						}
						else
						{
							stringstream ss;
							ss << "failed to schedule concurrent run for run id:" << run_id;
							report(ss.str(), false);
						}
					}
					n_concur = get_n_concurrent(run_id);
					if (n_concur == 0 && should_schedule)
					{
						waiting_runs.push_front(run_id);
					}
				}
			}
		}
		catch (exception &e)
		{
			cout << "exception trying to find overdue runs: " << endl << e.what() << endl;
		}
	}
    if ((!persistent_workers) && (waiting_runs.size() < free_agent_list.size()))
    {
        stringstream ss;
        free_agent_list = get_free_agent_list();
        int num_to_close = free_agent_list.size() - waiting_runs.size();
        list<AgentInfoRec>::iterator it_agent, iter_e;

        int closed = 0;

        for (int i=0;i<num_to_close;i++)
        {
            for (it_agent = agent_info_set.begin(), iter_e = agent_info_set.end();
                 it_agent != iter_e; ++it_agent) {
                AgentInfoRec::State state = it_agent->get_state();
                if (state == AgentInfoRec::State::WAITING) {
                    ss.str("");
                    ss << "using non-persistent agents, closed connection to agent: " << it_agent->get_socket_name()
                       << ", number of agents: " << socket_to_iter_map.size();
                    close_agent(it_agent);
                    closed++;


                    report(ss.str(), false);
                    break;
                }
            }
        }

    }
}

int RunManagerPanther::schedule_run(int run_id, std::list<list<AgentInfoRec>::iterator> &free_agent_list, int n_responsive_agents)
{
	int scheduled = -1;
	auto it_agent = free_agent_list.end(); // iterator to current socket
	int n_concurrent = get_n_concurrent(run_id);

	if (run_finished(run_id))
	{
		// run already completed on different node.  Do nothing
		scheduled = 0;
	}
	else if (failure_map.count(run_id) >= max_n_failure)
	{
		//if this run has already failed the max number of times, do nothing
		scheduled = 0;
	}
	else if (failure_map.count(run_id) == 0)// || failure_map.count(run_id) >= slave_fd.size())
	{
		// schedule a run on a slave
		it_agent = free_agent_list.begin();
		scheduled = -1;
	}
	else if (failure_map.count(run_id) + n_concurrent >= n_responsive_agents)
	{
		// enough enough slaves to make all failed runs on different slaves
		// schedule a run on a slave
		it_agent = free_agent_list.begin();
		scheduled = -1;
	}
	else if (failure_map.count(run_id) > 0)
	{
		for (it_agent = free_agent_list.begin(); it_agent != free_agent_list.end(); ++it_agent)
		{
			int socket_fd = (*it_agent)->get_socket_fd();
			auto fail_iter_pair = failure_map.equal_range(run_id);

			auto i = fail_iter_pair.first;
			for (i = fail_iter_pair.first;
				i != fail_iter_pair.second && i->second != socket_fd;
				++i) {
			}
			if (i == fail_iter_pair.second)  // This is slave has not previously failed on this run
			{
				// This run has not previously failed on this slave
				// Schedule run on it_sock
				break;
			}
		}
	}
	if (it_agent != free_agent_list.end())
	{
		int socket_fd = (*it_agent)->get_socket_fd();
		vector<char> data = file_stor.get_serial_pars(run_id);
		string info_txt;
		double info_val;
		int rstat;
		file_stor.get_info(run_id, rstat, info_txt, info_val);
		string host_name = (*it_agent)->get_hostname();

		//  info_txt = "sending run to " + host_name + ":" + (*it_agent)->get_work_dir() + " at " + pest_utils::get_time_string();
		NetPackage net_pack(NetPackage::PackType::START_RUN, cur_group_id, run_id, info_txt);
		pair<int,string> err = net_pack.send(socket_fd, &data[0], data.size());
		if (err.first > 0)
		{
			(*it_agent)->set_state(AgentInfoRec::State::ACTIVE, run_id, cur_group_id);
			//report("changed agent " + host_name + ":" + (*it_agent)->get_work_dir() + " to 'active'",false);
			//start run timer
			(*it_agent)->start_timer();
			//reset the last ping time so we don't ping immediately after run is started
			(*it_agent)->reset_last_ping_time();
			active_runid_to_iterset_map.insert(make_pair(run_id, *it_agent));
			stringstream ss;
			ss << "Sending run_id:" << run_id << " to:" << host_name << "$" << (*it_agent)->get_work_dir() <<
				" group_id:" << cur_group_id << " " << info_txt << " concurrent:" << get_n_concurrent(run_id) << ")";
			report(ss.str(), false);
			free_agent_list.erase(it_agent);
			scheduled = 1;
		}
		else
		{
			stringstream ss;
			ss << "error sending run_id:" << run_id << " to:" << host_name << "$" << (*it_agent)->get_work_dir() + " error:" + err.second;
			report(ss.str(), false);
		}
	}
	return scheduled;  // 1 = run scheduled; -1 failed to schedule run; 0 run not needed
}


void RunManagerPanther::echo()
{
	if (!should_echo)
		return;
    std::chrono::system_clock::time_point now = chrono::system_clock::now();
    if (chrono::duration_cast<std::chrono::milliseconds> ( now- last_echo_time).count() < echo_interval_milliseconds)
        return;
    last_echo_time = now;
	map<string, int> stats_map = get_agent_stats();
	cout << get_time_string_short() << " mn:" << setw(5) << setprecision(2) << left << get_global_runtime_minute()  << " runs("
	     << "C" << setw(5) << left << model_runs_done
		<< "|F" << setw(5) << left << model_runs_failed
		<< "|T" << setw(5) << left << model_runs_timed_out << ") agents("
		<< "R" << setw(4) << left << stats_map["run"]
		<< "|W" << setw(4) << left << stats_map["wait"]
		<< "|U" << setw(4) << left << stats_map["unavailable"] << ") " << setw(3) << left << open_file_trans_streams.size() << "\r" << flush;
}

void RunManagerPanther::report(std::string message,bool to_cout)
{
	string t_str = pest_utils::get_time_string();
	f_rmr << t_str << "->" << message << endl;
	if (to_cout) cout << endl << t_str << "->" << message << endl;
}

void RunManagerPanther::process_message(int i_sock)
{
	NetPackage net_pack;
	pair<int,string>  err;
	list<AgentInfoRec>::iterator agent_info_iter = socket_to_iter_map.at(i_sock);

	string host_name = agent_info_iter->get_hostname();
	string port_name = agent_info_iter->get_port();
	string socket_name = agent_info_iter->get_socket_name();
	err = net_pack.recv(i_sock);
	if( err.first <=0) // error or lost connection
	{
		if (err.first  == -2) {
			report("received corrupt message from agent:" + host_name + "$" + agent_info_iter->get_work_dir() + " error:" + err.second, false);
		}
		else if (err.first < 0) {
			report("receive failed from agent:" + host_name + "$" + agent_info_iter->get_work_dir() + " error:" + err.second, false);
		}
		else {
			report("lost connection to agent:" + host_name + "$" + agent_info_iter->get_work_dir() + " error:" + err.second, false);
		}
		close_agent(i_sock);
	}
	else if (net_pack.get_type() == NetPackage::PackType::TERMINATE)
	{
		report("agent exiting:" + host_name + "$" + agent_info_iter->get_work_dir() + " " + net_pack.get_info_txt(), false);
		close_agent(i_sock);
	}
	else if (net_pack.get_type() == NetPackage::PackType::CORRUPT_MESG)
	{
		report("agent reporting corrupt message:" + host_name + "$" + agent_info_iter->get_work_dir() + " " + net_pack.get_info_txt(), false);
		close_agent(i_sock);
	}
	else if (net_pack.get_type() == NetPackage::PackType::RUNDIR)
	{
		bool good_work_dir = NetPackage::check_string(net_pack.get_data(), 0, net_pack.get_data().size());
		if (good_work_dir)
		{
			string work_dir = NetPackage::extract_string(net_pack.get_data(), 0, net_pack.get_data().size());
			stringstream ss;
			ss << "initializing new agent connection from:" << agent_info_iter->get_hostname() << "$" << work_dir << " socket:" << socket_name << " number of agents:" << socket_to_iter_map.size();
			report(ss.str(), false);
			agent_info_iter->set_work_dir(work_dir);
			agent_info_iter->set_state(AgentInfoRec::State::CWD_RCV);
		}
		else
		{
			report("received corrupt run directory from agent:" + host_name + " - terminating agent", false);
			close_agent(i_sock);
		}
	}
	else if (net_pack.get_type() == NetPackage::PackType::LINPACK)
	{
		agent_info_iter->end_linpack();
		agent_info_iter->set_state(AgentInfoRec::State::LINPACK_RCV);
		stringstream ss;
		ss << "new agent ready:" << agent_info_iter->get_hostname() << "$" << agent_info_iter->get_work_dir() << " socket:" << agent_info_iter->get_socket_name() ;
		report(ss.str(), false);
	}
	else if (net_pack.get_type() == NetPackage::PackType::READY)
	{
		// ready message received from agent
		agent_info_iter->set_state(AgentInfoRec::State::WAITING);
	}

	else if ( (net_pack.get_type() == NetPackage::PackType::RUN_FINISHED
		|| net_pack.get_type() == NetPackage::PackType::RUN_FAILED
		|| net_pack.get_type() == NetPackage::PackType::RUN_KILLED)
			&& net_pack.get_group_id() != cur_group_id)
	{
		// this is an old run that did not finish on time
		// just ignore it
		int run_id = net_pack.get_run_id();
		int group_id = net_pack.get_group_id();
		stringstream ss;
		ss << " run_id:" << run_id << " received from unexpected group_id:" << group_id << " should be group: " << cur_group_id;
		ss << " from:" << agent_info_iter->get_hostname() << "$" << agent_info_iter->get_work_dir() << "...ignoring";
		report(ss.str(), false);
		//throw PestError(ss.str());
	}
	else if (net_pack.get_type() == NetPackage::PackType::RUN_FINISHED)
	{
		int run_id = net_pack.get_run_id();
		int group_id = net_pack.get_group_id();

		//check if this run already completed on another node
		if (run_finished(run_id))
		{
			stringstream ss;
			ss << "Previously completed run_id:" << run_id << " finished on:" << host_name << "$" << agent_info_iter->get_work_dir() <<
				"  run time:" << agent_info_iter->get_runtime_minute() << " min group_id:" << group_id <<
				" " << net_pack.get_info_txt() << " concurrent:" << get_n_concurrent(run_id);
			report(ss.str(), false);
		}
		else
		{
			// keep track of model run time
			agent_info_iter->end_run();
			stringstream ss;
			ss << " run_id:" << run_id << " received from:" << host_name << "$" << agent_info_iter->get_work_dir() <<
				"  run_time:" << agent_info_iter->get_runtime_minute() << "  avg_run_time:" << get_global_runtime_minute() << " group_id:" << group_id <<
				" " << net_pack.get_info_txt() << " concurrent:" << get_n_concurrent(run_id);
			report(ss.str(), false);
			process_model_run(i_sock, net_pack);
			ss.str("");
			ss << "run " << run_id << " processed";
			report(ss.str(), false);
		}



	}
	else if (net_pack.get_type() == NetPackage::PackType::RUN_FAILED)
	{
		int run_id = net_pack.get_run_id();
		int group_id = net_pack.get_group_id();
		int n_concur = get_n_concurrent(run_id);
		stringstream ss;

		if (!run_finished(run_id))
		{
			ss << " run_id:" << run_id << " failed on agent:" << host_name << "$" << agent_info_iter->get_work_dir() << " group_id:" <<
			group_id << " " << net_pack.get_info_txt() << " concurrent:" << n_concur;
			string netpack_message = net_pack.get_info_txt();
			if (netpack_message.size() > 0)
			{
				ss << " worker_message:" << netpack_message;
			}
			report(ss.str(), false);
			model_runs_failed++;
			update_run_failed(run_id, i_sock);
			auto it = get_active_run_iter(i_sock);
			unschedule_run(it);
			n_concur = get_n_concurrent(run_id);
			if (n_concur == 0 && (failure_map.count(run_id) < max_n_failure))
			{
				//put model run back into the waiting queue
				waiting_runs.push_front(run_id);
			}
		}
	}
	else if (net_pack.get_type() == NetPackage::PackType::DEBUG_FAIL_FREEZE)
	{
		stringstream ss;
		ss << "Frozen agent:" << host_name << "$" << agent_info_iter->get_work_dir() << " is frozen because of panther_debug_freeze_on_fail = true...";
		report(ss.str(), true);
		ss.str("");
		ss <<"closing connection to frozen agent:" << host_name << "$" << agent_info_iter->get_work_dir() << ", note: the agent will continue in a frozen state";
		report(ss.str(), true);
        close_agent(agent_info_iter);

	}
	else if (net_pack.get_type() == NetPackage::PackType::RUN_KILLED)
	{
		int run_id = net_pack.get_run_id();
		int group_id = net_pack.get_group_id();
		int n_concur = get_n_concurrent(run_id);
		auto it = get_active_run_iter(i_sock);
		unschedule_run(it);
		stringstream ss;
		ss << " run_id:" << run_id << " killed on agent:" << host_name << "$" << agent_info_iter->get_work_dir()
		<< " " << net_pack.get_info_txt() << " concurrent:" << n_concur;
		report(ss.str(), false);
	}
	else if (net_pack.get_type() == NetPackage::PackType::PING)
	{
		//report("ping received from agent:" + host_name + "$" + agent_info_iter->get_work_dir(), false);
#ifdef _DEBUG
		report("ping received from agent:" + host_name + "$" + agent_info_iter->get_work_dir(), false);
#endif
	}
	else if (net_pack.get_type() == NetPackage::PackType::IO_ERROR)
	{
		//string err(net_pack.get_data().begin(),net_pack.get_data().end());
		report("error in model IO files on agent:" + host_name + "$" + agent_info_iter->get_work_dir() + " socket:" + net_pack.get_info_txt() + " - terminating agent. ", true);
		close_agent(i_sock);
	}
	else if (net_pack.get_type() == NetPackage::PackType::START_FILE_WRKR2MSTR) {
        stringstream ss;
        ss.str("");
        pair<string, string> fnames = get_recv_filenames(net_pack, host_name,agent_info_iter->get_work_dir());
        if ((fnames.first.size() == 0) || (fnames.second.size() == 0)) {
            //do something here
        } 
		else {
            if (open_file_trans_streams.find(fnames.second) != open_file_trans_streams.end()) {
                ss.str("");
                ss << "file transfer error from agent:" << host_name << "$" << agent_info_iter->get_work_dir() << " agent_file:" << fnames.first << " master_file:" << fnames.second
                   << " - already open, can't reopen, something is wrong";
                report(ss.str(), true);
            } else {
                pair<map<string, ofstream *>::iterator, bool> ret = open_file_trans_streams.insert(
                        pair<string, ofstream *>(fnames.second, new ofstream));
                ofstream &out = *ret.first->second;
                out.open(fnames.second.c_str(), ios::binary);
                if (out.bad())
                {
                    ss.str("");
                    ss << " error opening master_file:" << fnames.second << " for writing" << endl;
                    report(ss.str(),true);
                }
                ss.str("");
                string agent_dir = host_name + "$" + agent_info_iter->get_work_dir();
                if (agent_dir.find(" ") != string::npos)
                {
                    agent_dir = "\"" + agent_dir + "\"";
                }

                ss << "file transfer request from agent:" << agent_dir;
                ss << " from agent_file:" << fnames.first << " to master_file:" << fnames.second << " " << net_pack.get_info_txt();

                report(ss.str(),false);
                open_file_socket_map.insert(make_pair(i_sock,fnames.second));
            }
        }
    }

    else if (net_pack.get_type() == NetPackage::PackType::CONT_FILE_WRKR2MSTR)
    {
        stringstream ss;
        ss.str("");
        pair<string,string> fnames = get_recv_filenames(net_pack, host_name,agent_info_iter->get_work_dir());
        if ((fnames.first.size() == 0) || (fnames.second.size() == 0))
        {
            //do something here
        }
        else
        {
            if (open_file_trans_streams.find(fnames.second) == open_file_trans_streams.end())
            {
                ss << "file transfer error from agent:" << host_name << "$" << agent_info_iter->get_work_dir() << " requesting agent_file:" << fnames.first << " to master_file:" << fnames.second << " - not open yet, can't continue transfer, something is wrong";
                report(ss.str(),true);
            }
            else
            {
                //pair<map<string ,ofstream*>::iterator, bool> ret = open_file_trans_streams.insert(pair<string,ofstream*>(fnames.second,new ofstream));
                //ofstream& out = *ret.first->second;
                ofstream& out = *open_file_trans_streams.at(fnames.second);
                vector<int8_t> ibuf = net_pack.get_data();
                //cout << reinterpret_cast<char*>(ibuf.data()) << endl;
				if (out.bad())
				{
					ss.str("");
					ss << " error writing to master_file:" << fnames.second ;
					report(ss.str(), true);
				}
                out.write(reinterpret_cast<char*>(ibuf.data()),ibuf.size());
				if (out.bad())
				{
					ss.str("");
					ss << " error writing to master_file:" << fnames.second;
					report(ss.str(), true);
				}
                out.flush();
                bytes_transferred += ibuf.size();
                ss.str("");
                
				int file_size = out.tellp();
                string agent_dir = host_name + "$" + agent_info_iter->get_work_dir();
                if (agent_dir.find(" ") != string::npos)
                {
                    agent_dir = "\"" + agent_dir + "\"";
                }
				ss << " agent:" << agent_dir << " transferred bytes:" << ibuf.size() << " to master_file:" << fnames.second << " total_bytes: " << file_size;
				report(ss.str(), false);

            }
        }
    }

    else if (net_pack.get_type() == NetPackage::PackType::FINISH_FILE_WRKR2MSTR)
    {
        stringstream ss;
        ss.str("");
        pair<string,string> fnames = get_recv_filenames(net_pack, host_name,agent_info_iter->get_work_dir());
        if ((fnames.first.size() == 0) || (fnames.second.size() == 0))
        {
            //do something here
        }
        else
        {
            if (open_file_trans_streams.find(fnames.second) == open_file_trans_streams.end())
            {
                string agent_dir = host_name + "$" + agent_info_iter->get_work_dir();
                if (agent_dir.find(" ") != string::npos)
                {
                    agent_dir = "\"" + agent_dir + "\"";
                }
                ss << "file transfer error from agent:" << agent_dir << " agent_file:" << fnames.first << "' master_file:" << fnames.second << " - not open yet, can't close, something is wrong";
                report(ss.str(),true);
            }
            else
            {
                //pair<map<string ,ofstream*>::iterator, bool> ret = open_file_trans_streams.insert(pair<string,ofstream*>(fnames.second,new ofstream));
                //ofstream& out = *ret.first->second;
                ofstream& out = *open_file_trans_streams.at(fnames.second);
                int file_size = out.tellp();
                out.flush();
                out.close();
                delete open_file_trans_streams.at(fnames.second);
                open_file_trans_streams.erase(fnames.second);
                string agent_dir = host_name + "$" + agent_info_iter->get_work_dir();
                if (agent_dir.find(" ") != string::npos)
                {
                    agent_dir = "\"" + agent_dir + "\"";
                }
                ss.str("");
                ss << "closed master_file:" << fnames.second << " for file transfer of agent_file:" << fnames.first;
                ss << " from agent:" << agent_dir << " transferred bytes:" << file_size;
                report(ss.str(),false);
                files_transferred += 1;
                open_file_socket_map.erase(i_sock);
            }
        }
    }


    else
	{
		report("received unsupported message from agent:" , false);
		net_pack.print_header(f_rmr);
		//save results from model run
	}
}

pair<string,string> RunManagerPanther::get_recv_filenames(NetPackage& net_pack, string hostname, string working_dir)
{
    //sanitize hostname and working_dir
    replace(hostname.begin(),hostname.end(),'/','-');
    replace(hostname.begin(),hostname.end(),'\\','-');
    replace(hostname.begin(),hostname.end(),'.','-');
    replace(hostname.begin(),hostname.end(),'(','-');
    replace(hostname.begin(),hostname.end(),')','-');
    replace(hostname.begin(),hostname.end(),':','-');
    replace(hostname.begin(),hostname.end(),' ','-');
    replace(hostname.begin(),hostname.end(),'\t','-');

    replace(working_dir.begin(),working_dir.end(),'/','-');
    replace(working_dir.begin(),working_dir.end(),'\\','-');
    replace(working_dir.begin(),working_dir.end(),'.','-');
    replace(working_dir.begin(),working_dir.end(),'(','-');
    replace(working_dir.begin(),working_dir.end(),')','-');
    replace(working_dir.begin(),working_dir.end(),':','-');
    replace(working_dir.begin(),working_dir.end(),' ','-');
    replace(working_dir.begin(),working_dir.end(),'\t','-');
    stringstream ss;
    ss.str("");
    ss << "hostname=" << hostname;
    ss << ".agentdir=" << working_dir;
    string info_txt = net_pack.get_info_txt();
    vector<string> tokens;
    tokenize(info_txt,tokens);

    ss << ".runid=" << net_pack.get_run_id() << ".groupid=" << net_pack.get_group_id();
    string agent_filename_token = "";
    string agent_filename = "";
    string master_filename = "";
    string new_master_filename = "";
    for (auto& token : tokens) {
        if (token.find("agent_filename:") != string::npos)
            agent_filename_token = token;
        else
            ss << "." << token;
    }
    if (agent_filename_token.size() == 0)
    {
        ss.str("");
        ss << "generic_ftx_" << nftx << ".dat";
        master_filename = ss.str();
        ss.str("");
        ss << "WARNING: could not find 'agent_filename' tag in net_pack info txt: '" << info_txt;
        ss << " from agent:" << hostname << "$" << working_dir << " -  unable form file transfer names";
        report (ss.str(), true);
        new_master_filename = master_filename;
    }
    else {
        ss << "." << agent_filename_token;
        tokens.clear();
        tokenize(agent_filename_token, tokens, ":");
        if (tokens.size() != 2) {
            report("WARNING: invalid 'agent_filename' token found:"+agent_filename_token,true);
            agent_filename = tokens[1];
            master_filename = ss.str();

        }
        else {
            agent_filename = tokens[1];
            master_filename = ss.str();
        }
        ss.str("");
        ss << "ftx_" << nftx << "." << agent_filename;
        new_master_filename = ss.str();
		
    }
	replace(master_filename.begin(), master_filename.end(), ':', '=');
    replace(new_master_filename.begin(),new_master_filename.end(),'/','-');
    replace(new_master_filename.begin(),new_master_filename.end(),'\\','-');
    replace(new_master_filename.begin(),new_master_filename.end(),'(','-');
    replace(new_master_filename.begin(),new_master_filename.end(),')','-');
    replace(new_master_filename.begin(),new_master_filename.end(),':','-');
    replace(new_master_filename.begin(),new_master_filename.end(),' ','-');
    replace(new_master_filename.begin(),new_master_filename.end(),'\t','-');

    if (org_new_master_fxt_map.find(master_filename) == org_new_master_fxt_map.end())
    {
        nftx++;
        org_new_master_fxt_map[master_filename] = new_master_filename;
    }
    else
    {
        new_master_filename = org_new_master_fxt_map.at(master_filename);
    }
    ss.str("");
    ss << " org_master_file:" << master_filename << " new master_file:" << new_master_filename;
    report(ss.str(),false);
    if (new_master_filename.size() > 255)
    {
        ss.str("");
        ss << "WARNING: master_filename:" << new_master_filename << " size > 255" << endl;
        report(ss.str(),true);
    }


    return make_pair(agent_filename,new_master_filename);
}


bool RunManagerPanther::process_model_run(int sock_id, NetPackage &net_pack)
{
	list<AgentInfoRec>::iterator agent_info_iter = socket_to_iter_map.at(sock_id);
	bool use_run = false;
	int run_id = net_pack.get_run_id();

	//check if another instance of this model run has already completed
	if (!run_finished(run_id))
	{
		Parameters pars;
		Observations obs;
		double run_time = 0;
		Serialization::unserialize(net_pack.get_data(), pars, get_par_name_vec(), obs, get_obs_name_vec(), run_time);
		file_stor.update_run(run_id, pars, obs);
		agent_info_iter->set_state(AgentInfoRec::State::COMPLETE);
		//slave_info_iter->set_state(SlaveInfoRec::State::WAITING);
		use_run = true;
		model_runs_done++;

	}
	else
	{
		stringstream ss;
		ss << " run_id:" << run_id << " already finished";
		report(ss.str(), false);
	}
	// remove currently completed run from the active list
	auto it = get_active_run_iter(sock_id);
	unschedule_run(it);
	kill_runs(run_id, false, "completed on alternative node");
	return use_run;
}

void RunManagerPanther::kill_run(list<AgentInfoRec>::iterator agent_info_iter, const string &reason)
{
	int socket_id = agent_info_iter->get_socket_fd();
	AgentInfoRec::State state = agent_info_iter->get_state();
	unschedule_run(agent_info_iter);
	if (socket_id && (state == AgentInfoRec::State::ACTIVE || state == AgentInfoRec::State::KILLED_FAILED))
	{
		int run_id = agent_info_iter->get_run_id();
		agent_info_iter->set_state(AgentInfoRec::State::KILLED);
		//schedule run to be killed
		string host_name = agent_info_iter->get_hostname();
		stringstream ss;
		ss << "sending kill request. reason:" << reason << " run_id:" << run_id;
		ss<< ",  num_previous_fails:" << failure_map.count(run_id) << " agent:" << host_name << "$" << agent_info_iter->get_work_dir();
		report(ss.str(), false);
		NetPackage net_pack(NetPackage::PackType::REQ_KILL, 0, 0, "");
		char data = '\0';
		pair<int,string> err = net_pack.send(socket_id, &data, sizeof(data));
		if (err.first == 1)
		{
			agent_info_iter->set_state(AgentInfoRec::State::KILLED);
		}
		else
		{
			report("error sending kill request to agent:" + host_name + "$" +
				agent_info_iter->get_work_dir() + " error:" + err.second, true);
			agent_info_iter->set_state(AgentInfoRec::State::KILLED_FAILED);
		}
	}
}


void RunManagerPanther::kill_runs(int run_id, bool update_failure_map, const string &reason)
{
	auto range_pair = active_runid_to_iterset_map.equal_range(run_id);
	//runs with this id are not needed so kill them
	list<list<AgentInfoRec>::iterator> kill_list;

	for (auto b = range_pair.first; b != range_pair.second; ++b)
	{
		list<AgentInfoRec>::iterator agent_info_iter = (*b).second;
		kill_list.push_back(agent_info_iter);
	}
	for (auto &iter : kill_list)
	{
		kill_run(iter, reason);
		if (update_failure_map) update_run_failed(run_id, iter->get_socket_fd());
	}
}


void RunManagerPanther::kill_all_active_runs()
{
	list<list<AgentInfoRec>::iterator> iter_list;
	list<AgentInfoRec>::iterator iter_b, iter_e;
	bool active_runs = true;
	for (int n_tries = 0; active_runs && n_tries >= 100; ++n_tries)
	{
		init_agents();
		active_runs = false;
		for (iter_b = agent_info_set.begin(), iter_e = agent_info_set.end();
			iter_b != iter_e; ++iter_b)
		{
			int socket_id = iter_b->get_socket_fd();
			AgentInfoRec::State state = iter_b->get_state();
			if (socket_id && (state == AgentInfoRec::State::ACTIVE || state == AgentInfoRec::State::KILLED_FAILED))
			{
				active_runs = true;
				kill_run(iter_b, "completed run group");
			}
		}
		listen();
	}
}

 void RunManagerPanther::init_agents(pest_utils::thread_flag* terminate/* = nullptr*/)
 {
	 for (auto &i_agent : agent_info_set)
	 {
	 	// Stop early if we're requested to terminate
	 	if(terminate && terminate->get())
	 	{
	 		break;
	 	}

		int i_sock = i_agent.get_socket_fd();
		AgentInfoRec::State cur_state = i_agent.get_state();
		if (cur_state == AgentInfoRec::State::NEW)
		{
			NetPackage net_pack(NetPackage::PackType::REQ_RUNDIR, 0, 0, "");
			char data = '\0';
			pair<int,string> err = net_pack.send(i_sock, &data, sizeof(data));
			if (err.first > 0)
			{
				i_agent.set_state(AgentInfoRec::State::CWD_REQ);
			}
		}
		else if (cur_state == AgentInfoRec::State::CWD_RCV)
		{
			// send parameter and observation names
			NetPackage net_pack(NetPackage::PackType::PAR_NAMES, 0, 0, "");
			vector<int8_t> data;
			vector<string> tmp_vec;
			// send parameter names
			if (par_names_to_check_worker.size() > 0)
				tmp_vec = par_names_to_check_worker;
			else
				tmp_vec = file_stor.get_par_name_vec();
			data = Serialization::serialize(tmp_vec);
			pair<int,string> err_par = net_pack.send(i_sock, &data[0], data.size());
			
			//send observation names
			net_pack = NetPackage(NetPackage::PackType::OBS_NAMES, 0, 0, "");
			// send parameter names
			if (obs_names_to_check_worker.size() > 0)
				tmp_vec = obs_names_to_check_worker;
			else
				tmp_vec = file_stor.get_obs_name_vec();
			data = Serialization::serialize(tmp_vec);
			pair<int,string> err_obs = net_pack.send(i_sock, &data[0], data.size());

			if (err_par.first > 0 && err_obs.first > 0)
			{
				i_agent.set_state(AgentInfoRec::State::NAMES_SENT);
			}
			else if (err_par.first <= 0)
			{
				report("Error sending par names to agent:" + i_agent.get_hostname() + "$" + i_agent.get_work_dir() + " error:" + err_par.second,false);
			}
			else if (err_obs.first <= 0)
			{
				report("Error sending obs names to agent:" + i_agent.get_hostname() + "$" + i_agent.get_work_dir() + " error:" + err_obs.second, false);
			}
		}
		else if (cur_state == AgentInfoRec::State::NAMES_SENT)
		{
			NetPackage net_pack(NetPackage::PackType::REQ_LINPACK, 0, 0, "");
			char data = '\0';
			pair<int,string> err = net_pack.send(i_sock, &data, sizeof(data));
			if (err.first  > 0)
			{
				i_agent.set_state(AgentInfoRec::State::LINPACK_REQ);
				i_agent.start_timer();
			}
			else
			{
				report("error sending linpack request to agent:" + i_agent.get_hostname() + "$" + i_agent.get_work_dir() + " error:" + err.second,false);
			}
		}
		else if (cur_state == AgentInfoRec::State::LINPACK_RCV)
		{
			i_agent.set_state(AgentInfoRec::State::WAITING);
		}
		/*else if (cur_state == AgentInfoRec::State::COMPLETE)
		{
			i_agent.set_state(AgentInfoRec::State::WAITING);
		}*/
	}
 }

 vector<int> RunManagerPanther::get_overdue_runs_over_kill_threshold(int run_id)
 {
	 vector<int> sock_id_vec;
	 auto range_pair = active_runid_to_iterset_map.equal_range(run_id);


	 bool just_quit = false;
     int q = pest_utils::quit_file_found();
     if ((q == 1) || (q == 2) || (q == 4)) {
	     just_quit = true;
     }

	 double duration;
	 for (auto &i = range_pair.first; i != range_pair.second; ++i)
	 {
		 if (i->second->get_state() == AgentInfoRec::State::ACTIVE)
		 {
			 double avg_runtime = i->second->get_runtime_minute();
			 if (avg_runtime <= 0) avg_runtime = get_global_runtime_minute();;
			 if (avg_runtime <= 0) avg_runtime = 1.0E+10;
			 duration = i->second->get_duration_minute();
			 if ((just_quit) || (duration > overdue_giveup_minutes) ||
                     ((duration >= avg_runtime*overdue_giveup_fac) &&
                     (avg_runtime > MIN_AVGRUNMINS_FOR_KILL)))
			 {
				 sock_id_vec.push_back(i->second->get_socket_fd());
			 }
		 }
	 }
	 return sock_id_vec;
 }

 bool RunManagerPanther::all_runs_complete()
 {
     int q = pest_utils::quit_file_found();
     if ((q == 1) || (q == 2) || (q == 4))
    {

        kill_all_active_runs();
        for (auto run_id : waiting_runs)
        {
            file_stor.update_run_failed(run_id);
        }
        waiting_runs.clear();
        cout << endl << "'pest.stp' found, all remaining runs marked as fails. " << endl << endl;
        return true;
    }
	 // check for run in the waitng queue
	 if (!waiting_runs.empty())
	 {
		 return false;
	 }
	 // check for active runs
	 for (auto it_active = active_runid_to_iterset_map.begin(); it_active != active_runid_to_iterset_map.end(); ++it_active)
	 {
		 if (it_active->second->get_state() == AgentInfoRec::State::ACTIVE)
		 {
			 return false;
		 }
	 }
	 return true;
 }


 list<AgentInfoRec>::iterator RunManagerPanther::add_agent(int sock_id)
 {
	 stringstream ss;
	 ss << "new connection from:" << w_getnameinfo_string(sock_id);
	 report(ss.str(), false);
	 FD_SET(sock_id, &master); // add to master set
	 if (sock_id > fdmax) { // keep track of the max
		 fdmax = sock_id;
	 }

	 //list<SlaveInfoRec>::iterator
	agent_info_set.push_back(AgentInfoRec(sock_id));
	list<AgentInfoRec>::iterator iter = std::prev(agent_info_set.end());
	socket_to_iter_map[sock_id] = iter;
	return iter;
 }

 double RunManagerPanther::get_global_runtime_minute() const
 {
	 double global_runtime = 0;
	 double temp = 0;
	 int count = 0;
	 for (auto &si : agent_info_set)
	 {
		 temp = si.get_runtime_minute();
		 if (temp > 0)
		 {
			 count++;
			 global_runtime += temp;
		 }
	 }
	 if (count == 0)
		 return 0.0;
	 return global_runtime / (double)count;
 }

 void RunManagerPanther::unschedule_run(list<AgentInfoRec>::iterator agent_info_iter)
 {
	 int run_id = agent_info_iter->get_run_id();
	 auto range_pair = active_runid_to_iterset_map.equal_range(run_id);

	 for (auto iter = range_pair.first; iter != range_pair.second;)
	 {
		 if (iter->second == agent_info_iter)
		 {
			 iter = active_runid_to_iterset_map.erase(iter);
			 return;
		 }
		 else
		 {
			 ++iter;
		 }
	 }
 }

 list<list<AgentInfoRec>::iterator> RunManagerPanther::get_free_agent_list()
 {
	 list<list<AgentInfoRec>::iterator> iter_list;
	 list<AgentInfoRec>::iterator iter_b, iter_e;
	 stringstream ss;
	 for (iter_b = agent_info_set.begin(), iter_e = agent_info_set.end();
		 iter_b != iter_e; ++iter_b)
	 {
		 AgentInfoRec::State cur_state = iter_b->get_state();
		 /*ss.str("");
		 ss << iter_b->get_hostname() << ":" << iter_b->get_work_dir() << "," << iter_b->state_strings[static_cast<int>(cur_state)];
		 report(ss.str(), false);*/
		 
		 if (cur_state == AgentInfoRec::State::WAITING)
		 {
			 iter_list.push_back(iter_b);
		 }
	 }
	 return iter_list;
 }

 map<string, int> RunManagerPanther::get_agent_stats()
 {
	 map<string, int> stats_map;
	 list<AgentInfoRec>::iterator iter_b, iter_e;
	 int n_active = 0;
	 int n_waiting = 0;
	 int n_unavailable = 0;
	 for (iter_b = agent_info_set.begin(), iter_e = agent_info_set.end();
		 iter_b != iter_e; ++iter_b)
	 {
		 AgentInfoRec::State cur_state = iter_b->get_state();
		 if (cur_state == AgentInfoRec::State::WAITING)
		 {
			 ++n_waiting;
		 }
		 else if (cur_state == AgentInfoRec::State::ACTIVE)
		 {
			 ++n_active;
		 }
		 else
		 {
			 ++n_unavailable;
		 }
	 }
	 stats_map["wait"] = n_waiting;
	 stats_map["run"] = n_active;
	 stats_map["unavailable"] = n_unavailable;
	 stats_map["total"] = agent_info_set.size();
	 return stats_map;
 }

 int RunManagerPanther::get_n_unique_failures()
 {
	 set<int> run_id_set;
	 for (const auto &i_fail : failure_map)
	 {
		 run_id_set.insert(i_fail.first);
	 }
	 return run_id_set.size();
 }

 int RunManagerPanther::get_n_responsive_agents()
 {
	 int n = 0;
	 for (const auto &i : agent_info_set)
	 {
		 if (i.get_failed_pings() < N_PINGS_UNRESPONSIVE) ++n;
	 }
	 return n;
 }


 void RunManagerPanther::update_run_failed(int run_id, int socket_fd)
 {
	 file_stor.update_run_failed(run_id);
	 failure_map.insert(make_pair(run_id, socket_fd));
	 list<AgentInfoRec>::iterator agent_info_iter = socket_to_iter_map.at(socket_fd);
	 agent_info_iter->add_failed_run();
 }

 void RunManagerPanther::update_run_failed(int run_id)
 {
	 // must call void RunManagerPANTHER::update_run_failed(int run_id, int socket_fd) instead
	 throw(PestError("Error: Unsupported function call  RunManagerPANTHER::update_run_failed(int run_id)"  ));
 }

RunManagerPanther::~RunManagerPanther(void)
{
	// Shut down idle agent management thread
	end_run_idle_async();

	//close sockets and cleanup
	int err;
	err = w_close(listener);
	FD_CLR(listener, &master);
	// this is needed to ensure that the first slave closes properly
	w_sleep(get_current_sleep_timeout_milliseconds(timeout_milliseconds));
	for (int i = 0; i <= fdmax; i++)
	{
		if (FD_ISSET(i, &master))
		{
			NetPackage netpack(NetPackage::PackType::TERMINATE, 0, 0,"");
			char data;
			netpack.send(i, &data, 0);
			err = w_close(i);
			FD_CLR(i, &master);
		}
	}
	w_cleanup();
    delete idle_thread;

}

RunManagerYAMRCondor::RunManagerYAMRCondor(const std::string & stor_filename,
	const std::string & port, std::ofstream & _f_rmr, int _max_n_failure,
	double overdue_reched_fac, double overdue_giveup_fac, double overdue_giveup_minutes, string _condor_submit_file): RunManagerPanther(stor_filename,
		port,_f_rmr,_max_n_failure,overdue_reched_fac,overdue_giveup_fac, overdue_giveup_minutes)
{
	submit_file = _condor_submit_file;
	parse_submit_file();
}

void RunManagerYAMRCondor::run()
{
	int cluster = submit();
	cout << " on condor cluster " << cluster << endl;
	RunManagerPanther::run();
	cleanup(cluster);

}

void RunManagerYAMRCondor::write_submit_file()
{
	ofstream f_out("temp.sub");
	if (!f_out.good())
		throw runtime_error("error opening temp.sub for writing");
	for (auto &line : submit_lines)
		f_out << line << endl;
	int n_q = min(max_condor_queue, get_n_waiting_runs());
	cout << "queueing "  << n_q << " agents " ;
	f_out << "queue " << n_q << endl;

}

int RunManagerYAMRCondor::get_cluster()
{
	string line, lower_line;
	int cluster = -999;
	string tag = "submitted to cluster";
	vector<string> tokens;

	//first check for err from condor_submit
	vector<string> err_lines;
	ifstream f_err("cs_temp.stderr");
	if (f_err.good())
	{
		while (getline(f_err, line))
		{
			err_lines.push_back(line);
		}

		if (err_lines.size() > 0)
		{
			stringstream ss;
			for (auto &l : err_lines)
				ss << l << endl;
			throw runtime_error("condor_submit issued error: " + ss.str());
		}
	}
	f_err.close();

	ifstream f_in("cs_temp.stdout");
	if (!f_in.good())
		throw runtime_error("error opening cs_temp.stdout to read cluster info from condor_submit");
	double temp;
	while (getline(f_in, line))
	{
		if (line.find(tag) != string::npos)
		{
			pest_utils::tokenize(pest_utils::strip_cp(line), tokens);
			try
			{
				pest_utils::convert_ip(tokens[tokens.size() - 1], temp);
				cluster = int(temp);
			}
			catch (exception &e)
			{
				throw runtime_error("error parsing '" + tokens[tokens.size() - 1] + "' to int on line: " + line);
			}
			break;
		}
	}
	f_in.close();
	if (cluster == -999)
	{
		throw runtime_error("cluster number not found in cs_temp.stdout");
	}
	return cluster;
}

int RunManagerYAMRCondor::submit()
{
	write_submit_file();
	stringstream ss;
	system("condor_submit temp.sub 1>cs_temp.stdout 2>cs_temp.stderr");
	return get_cluster();
}

void RunManagerYAMRCondor::cleanup(int cluster)
{
	RunManagerPanther::close_agents();
	stringstream ss;
	ss << "condor_rm " << cluster << " 1>cr_temp.stdout 2>cr_temp.stderr";
	system(ss.str().c_str());
	w_sleep(1000);
	ss.str(string());
	ss << "condor_rm " << cluster << " -forcex 1>cr_temp.stdout 2>cr_temp.stderr";
	w_sleep(1000);
	system(ss.str().c_str());
	RunManagerPanther::close_agents();
	cout << "   all agents freed " << endl << endl;
}

void RunManagerYAMRCondor::parse_submit_file()
{
	ifstream f_in(submit_file);
	if (!f_in.good())
		throw runtime_error("error opening submit file '" + submit_file + "' for reading");
	string line,lower_line,q_line;
	string q_tag = "queue";
	vector<string> tokens;
	while (getline(f_in, line))
	{
		pest_utils::strip_ip(line);
		//check if this line starts with 'queue'
		lower_line = pest_utils::lower_cp(line);
		if (lower_line.compare(0, q_tag.size(), q_tag) == 0)
		{
			q_line = line;
			pest_utils::tokenize(line, tokens);
		}
		else
			submit_lines.push_back(line);
	}
	f_in.close();
	if (tokens.size() == 0)
		throw runtime_error("'queue' line not found in submit file " + submit_file);
	else
	{
		try
		{
			pest_utils::convert_ip(tokens[1], max_condor_queue);
		}
		catch (exception &e)
		{
			runtime_error("error converting '" + tokens[2] + "' from line '" + q_line + "' to int");
		}
	}
}

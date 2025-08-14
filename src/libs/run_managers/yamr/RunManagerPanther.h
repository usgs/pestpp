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
#ifndef RUNMANAGERPANTHER_H
#define RUNMANAGERPANTHER_H
#include "network_wrapper.h"
#include <string>
#include <set>
#include <deque>
#include <unordered_map>
#include <chrono>
#include <list>
#include <thread>
#include "network_wrapper.h"
#include "network_package.h"
#include "RunManagerAbstract.h"
#include "RunStorage.h"
#include "utilities.h"

class AgentInfoRec {
public:
	static const int UNKNOWN_ID = -9999;
	enum class State { NEW, CWD_REQ, CWD_RCV, NAMES_SENT, LINPACK_REQ, LINPACK_RCV, WAITING, ACTIVE, KILLED, KILLED_FAILED, COMPLETE};
	AgentInfoRec(int _socket_fd);
	std::vector<std::string> state_strings;
	int get_socket_fd() const;
	string get_hostname()const;
	string get_port()const;
	string get_socket_name()const;
	void set_socket_fd(int _socket_fd);
	int get_run_id() const;
	void set_run_id(int _run_id);
	int get_group_id() const;
	void set_group_id(int _group_id);
	State get_state() const;
	void set_state(const State &_state);
	void set_state(const State &_state, int run_id, int group_id);
	void set_work_dir(const std::string & wkd);
	std::string get_work_dir() const;
	void start_timer();
	void end_run();
	void end_linpack();
	double get_runtime() const;
	double get_duration_sec() const;
	double get_duration_minute() const;
	double get_runtime_sec() const;
	double get_runtime_minute() const;
	double get_linpack_time() const;
	int add_failed_ping();
	int add_failed_run();
	int get_failed_runs() const { return failed_runs; }
	void set_ping(bool val);
	bool get_ping() const;
	int get_failed_pings() const;
	void reset_failed_pings();
	void reset_last_ping_time();
	void reset_runtime() { run_time = std::chrono::system_clock::duration::zero(); }
	int seconds_since_last_ping_time() const;
	~AgentInfoRec(){}
private:
	int socket_fd;
	int run_id;
	int group_id;
	bool ping;
	int failed_pings;
	int failed_runs;
	State state;
	std::chrono::system_clock::duration linpack_time;
	std::chrono::system_clock::duration run_time;
	std::chrono::system_clock::time_point start_time;
	std::chrono::system_clock::time_point last_ping_time;
	std::string work_dir;
	std::vector<string> name_info_vec;
public:
	class CompareTimes
	{
	public:
		CompareTimes() {}
		bool operator() (const AgentInfoRec &a, const AgentInfoRec &b);
	};
};


class RunManagerPanther : public RunManagerAbstract
{
public:
	RunManagerPanther(const std::string &stor_filename, const std::string &port, std::ofstream &_f_rmr, int _max_n_failure,
		double overdue_reched_fac, double overdue_giveup_fac, double overdue_giveup_minutes, bool _should_echo = true, const vector<string>& par_names=vector<string>(),
		const vector<string>& obs_names=vector<string>(),int _timeout_milliseconds=10,int _echo_interval_milliseconds=10, bool _persistent_workers=true);

	virtual void initialize(const Parameters &model_pars, const Observations &obs, const std::string &_filename = std::string(""));
	virtual void initialize_restart(const std::string &_filename);
	virtual void reinitialize(const std::string &_filename = std::string(""));
	virtual void free_memory();
	virtual int add_run(const Parameters &model_pars, const std::string &info_txt="", double info_value=RunStorage::no_data);
	virtual int add_run(const std::vector<double> &model_pars, const std::string &info_txt="", double info_valuee=RunStorage::no_data);
	virtual int add_run(const Eigen::VectorXd &model_pars, const std::string &info_txt="", double info_valuee=RunStorage::no_data);
	virtual void update_run(int run_id, const Parameters &pars, const Observations &obs);
	virtual void run();
	virtual RunManagerAbstract::RUN_UNTIL_COND run_until(RUN_UNTIL_COND condition, int n_nops = 0, double sec = 0.0);
	~RunManagerPanther(void);
	int get_n_waiting_runs() { return waiting_runs.size(); }
	void close_agents();

private:
	std::string port;
	static const int BACKLOG;
	static const int MAX_FAILED_PINGS;
	static const int N_PINGS_UNRESPONSIVE;
	static const int MIN_PING_INTERVAL_SECS;
	static const int MAX_PING_INTERVAL_SECS;
	static const int MAX_CONCURRENT_RUNS_LOWER_LIMIT;
	static const int IDLE_THREAD_SIGNAL_TIMEOUT_SECS;
    static const double MIN_AVGRUNMINS_FOR_KILL;
    //static const int MILLISECONDS_BETWEEN_ECHOS;
    //static const int TIMEOUT_MILLISECONDS;
    int echo_interval_milliseconds;
    int timeout_milliseconds;
	double overdue_reched_fac;
	double overdue_giveup_fac;
	double overdue_giveup_minutes;
	int max_concurrent_runs;
	int n_no_ops;  //number of consecutive times tcp/ip has looked for slave communciations and not found any
	int listener;
	int fdmax;
	int model_runs_done;
	int model_runs_failed;
	int model_runs_timed_out;
	long long bytes_transferred;
	int files_transferred;
	bool should_echo;
    std::chrono::system_clock::time_point last_echo_time;
	int nftx;
	fd_set master; // master file descriptor list
	list<AgentInfoRec> agent_info_set;
	map<int, list<AgentInfoRec>::iterator> socket_to_iter_map;
	multimap<int, list<AgentInfoRec>::iterator> active_runid_to_iterset_map;
	std::deque<int> waiting_runs;
	std::unordered_multimap<int, int> failure_map;
	pest_utils::thread_flag terminate_idle_thread;
	pest_utils::thread_flag currently_idle;
	pest_utils::thread_flag idling;
	pest_utils::thread_flag idle_thread_finished;
	thread* idle_thread;
	map<string,ofstream*> open_file_trans_streams;
	map<int,string> open_file_socket_map;
	//pest_utils::thread_RAII* idle_thread_raii;
    bool persistent_workers;

	int schedule_run(int run_id, std::list<list<AgentInfoRec>::iterator> &free_agent_list, int n_responsive_agents);
	void unschedule_run(list<AgentInfoRec>::iterator agent_info_iter);
	void kill_run(list<AgentInfoRec>::iterator agent_info_iter, const std::string &reason="UNKNOWN");
	void kill_runs(int run_id, bool update_failure_map, const std::string &reason = "UNKNOWN");
	void kill_all_active_runs();
	void close_agent(int i_sock);
	void close_agent(list<AgentInfoRec>::iterator agent_info_iter);

	void run_idle_async();
	void start_run_idle_async();
	void end_run_idle_async();
	void pause_idle();
	void resume_idle();
    int get_current_sleep_timeout_milliseconds(const int org_timeout_milliseconds);

    std::ofstream &f_rmr;
	bool listen(pest_utils::thread_flag* terminate = nullptr);
	bool process_model_run(int sock_id, NetPackage &net_pack);
	void process_message(int i);
	void schedule_runs();
	void init_agents(pest_utils::thread_flag* terminate = nullptr);
	list<AgentInfoRec>::iterator add_agent(int sock_id);
	//void erase_agent(int sock_id);
	bool ping(int i_sock);
	bool ping(pest_utils::thread_flag* terminate = nullptr);
	void report(std::string message,bool to_cout);
	
	/*string get_time_string();
	string get_time_string_short();*/
	void echo();
	vector<int> get_overdue_runs_over_kill_threshold(int run_id);
	bool all_runs_complete();
	list<AgentInfoRec>::iterator get_active_run_iter(int socket);
	std::list<std::list<AgentInfoRec>::iterator> get_free_agent_list();
	double get_global_runtime_minute() const;
	int get_n_concurrent(int run_id);
	int get_n_unique_failures();
	int get_n_responsive_agents();
	virtual void update_run_failed(int run_id, int socket_fd);
	virtual void update_run_failed(int run_id);
	map<string, int> get_agent_stats();
	vector<string> par_names_to_check_worker;
	vector<string> obs_names_to_check_worker;
    pair<string,string> get_recv_filenames(NetPackage& net_pack, string host_name, string working_dir);
    map<string,string> org_new_master_fxt_map;
};

class RunManagerYAMRCondor : public RunManagerPanther
{
public:
	RunManagerYAMRCondor(const std::string &stor_filename, const std::string &port, std::ofstream &_f_rmr, int _max_n_failure,
		double overdue_reched_fac, double overdue_giveup_fac, double overdue_giveup_minutes, string _condor_submit_file);
	virtual void run();

private:
	int max_condor_queue;
	vector<string> submit_lines;
	void parse_submit_file();
	int get_cluster();
	string submit_file;
	void write_submit_file();
	int submit();
	void cleanup(int cluster);

};

#endif /* RUNMANAGERPANTHER_H */

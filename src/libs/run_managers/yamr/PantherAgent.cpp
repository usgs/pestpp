#include "PantherAgent.h"
#include "utilities.h"
#include "Serialization.h"
#include "system_variables.h"
#include <cassert>
#include <cstring>
#include <algorithm>
#include <thread>
#include "system_variables.h"
#include "utilities.h"
#include <regex>
#include "network_package.h"
#include "OutputFileWriter.h"
#include "Pest.h"
#include "config_os.h"
#include "pest_data_structs.h"

using namespace pest_utils;

int  linpack_wrap(void);

PANTHERAgent::PANTHERAgent(ofstream &_frec)
	: frec(_frec),
	  max_time_without_master_ping_seconds(300),
	  restart_on_error(false),
	  current_da_cycle(NetPackage::NULL_DA_CYCLE)
{
}


void PANTHERAgent::init_network(const string &host, const string &port)
{
	report("initializing network connection", true);
	w_init();

	stringstream ss;
	pair<int,string> status;
	struct addrinfo hints;
	struct addrinfo *servinfo;
	//cout << "setting hints" << endl;
	memset(&hints, 0, sizeof hints);
	//Use this for IPv4 and IPv6
	//hints.ai_family = AF_UNSPEC;
	//Use this just for IPv4;
	hints.ai_family = AF_INET;
	hints.ai_socktype = SOCK_STREAM;
	hints.ai_flags = AI_PASSIVE;
	//cout << "atttemping w_getaddrinfo" << endl;
	status = w_getaddrinfo(host.c_str(), port.c_str(), &hints, &servinfo);
	if (status.first != 0)
	{
		ss.str("");
		ss << "ERROR: getaddrinfo returned non-zero: " << status.second;
		report(ss.str(), true);
		throw(PestError("ERROR: getaddrinfo returned non-zero: " + status.second));
	}
	w_print_servinfo(servinfo, cout);
	cout << endl;
	// connect
	ss.str("");
	ss << "PANTHER Agent will poll for master connection every " << poll_interval_seconds << " seconds" << endl;
	report(ss.str(), true);
	addrinfo* connect_addr = nullptr;
	while  (connect_addr == nullptr)
	{

		connect_addr = w_connect_first_avl(servinfo, sockfd);
		if (connect_addr == nullptr) {
			report("failed to connect to master", true);
			w_sleep(poll_interval_seconds * 1000);

		}

	}
	ss.str("");

	ss << "connection to master succeeded on socket: " << w_get_addrinfo_string(connect_addr) << endl << endl;
	report(ss.str(), true);
	freeaddrinfo(servinfo);

	fdmax = sockfd;
	FD_ZERO(&master);
	FD_SET(sockfd, &master);
	// send run directory to master
}


PANTHERAgent::~PANTHERAgent()
{
	w_close(sockfd);
	w_cleanup();
}


void PANTHERAgent::process_ctl_file(const string &ctl_filename)
{
	string version = PESTPP_VERSION;
	frec << "panther agent starting..." << endl;
	frec << "using control file: \"" << ctl_filename << "\"" << endl << endl;
	frec << endl << endl << "version: " << version << endl;
	frec << "binary compiled on " << __DATE__ << " at " << __TIME__ << endl << endl;
	frec << "in directory: \"" << OperSys::getcwd() << "\"" << endl << endl;

	cout << "panther agent starting..." << endl;
	cout << "using control file: \"" << ctl_filename << "\"" << endl << endl;
	cout << "in directory: \"" << OperSys::getcwd() << "\"" << endl << endl;

	report("processing control file", true);

	ifstream fin;
	long lnum;
	long sec_begin_lnum;
	long sec_lnum;
	string section("");
	string line;
	string line_upper;
	vector<string> tokens;

	int num_par;
	int num_tpl_file;
	std::vector<std::string> pestpp_lines;

	fin.open(ctl_filename);
	if (!fin)
	{
		report("ERROR: PANTHER agent unable to open pest control file : " + ctl_filename, true);
		throw PestError("PANTHER agent unable to open pest control file: " + ctl_filename);
	}
    pest_scenario.set_default_dynreg();
	pest_scenario.process_ctl_file(fin,ctl_filename,frec);
	if ((pest_scenario.get_ctl_parameters().size() > 250000) || (pest_scenario.get_ctl_observations().size() > 250000))
	{
		set<string> pargs = pest_scenario.get_pestpp_options().get_passed_args();
		if (pargs.find("CHECK_TPLINS") == pargs.end())
		{
			pest_scenario.get_pestpp_options_ptr()->set_check_tplins(false);
			cout << "large problem detected, resetting 'check_tplins' to false" << endl;
 		}
	}
	report("checking model interface files", true);
	pest_scenario.check_io(frec);
	pest_scenario.release_unused_for_agent();
	poll_interval_seconds = pest_scenario.get_pestpp_options().get_worker_poll_interval();

	mi = ModelInterface(pest_scenario.get_model_exec_info().tplfile_vec, 
		pest_scenario.get_model_exec_info().inpfile_vec,
		pest_scenario.get_model_exec_info().insfile_vec,
		pest_scenario.get_model_exec_info().outfile_vec,
		pest_scenario.get_model_exec_info().comline_vec);
	mi.set_additional_ins_delimiters(pest_scenario.get_pestpp_options().get_additional_ins_delimiters());
	mi.set_fill_tpl_zeros(pest_scenario.get_pestpp_options().get_fill_tpl_zeros());
	mi.set_tpl_force_decimal(pest_scenario.get_pestpp_options().get_tpl_force_decimal());
	mi.set_num_threads(pest_scenario.get_pestpp_options().get_num_tpl_ins_threads());
    mi.set_sleep_ms(5);
	restart_on_error = pest_scenario.get_pestpp_options().get_panther_agent_restart_on_error();
	max_time_without_master_ping_seconds = pest_scenario.get_pestpp_options().get_panther_agent_no_ping_timeout_secs();
	FileManager fm("panther_agent");
	OutputFileWriter of(fm, pest_scenario);
	of.scenario_report(frec);
	//pest_scenario.clear_ext_files();
}

pair<int,string> PANTHERAgent::recv_message(NetPackage &net_pack, struct timeval *tv)
{
	fd_set read_fds;
	std::pair<int, string> err;
	err.first = -1;
	int recv_fails = 0;
	stringstream ss;
	while (recv_fails < max_recv_fails && err.first != 1)
	{
		read_fds = master; // copy master
		int result = w_select(fdmax + 1, &read_fds, NULL, NULL, tv);
		if (result == -1)
		{
			ss.str("");

			ss << "fatal network error while receiving messages. ERROR: select() failure";
			report(ss.str(), true);
			return pair<int,string> (-990, "fatal network error while receiving messages. ERROR: select() failure");
		}
		if (result == 0)
		{
			// no messages available for reading
			if (tv == NULL)
			{
				ss.str("");
				ss << "fatal network error while receiving messages. ERROR: blocking select() call failure";
				report(ss.str(), true);
				return pair<int, string>(-990, "fatal network error while receiving messages. ERROR: blocking select() call failure");
			}
			else
			{
				return pair<int, string>(2, err.second);
			}
		}
		for (int i = 0; i <= fdmax; i++) {
			if (FD_ISSET(i, &read_fds)) { // got message to read
				err = net_pack.recv(i); // error or lost connection
				if (err.first == -2) {
					vector<string> sock_name = w_getnameinfo_vec(i);
					ss.str("");
					ss << "received corrupt message from master: " << sock_name[0] << ":" << sock_name[1] << ": " << err.second << endl;
					report(ss.str(), true);
					//w_close(i); // bye!
					//FD_CLR(i, &master); // remove from master set
					err.first = -999;
					return err;
				}
				else if (err.first < 0) {
					recv_fails++;
					vector<string> sock_name = w_getnameinfo_vec(i);
					ss.str("");
					ss << "receive from master failed: " << sock_name[0] << ":" << sock_name[1] << endl;
					report(ss.str(), true);
					err.second = "receive from master failed";
					err.first = -1;
				}
				else if(err.first == 0) {
					vector<string> sock_name = w_getnameinfo_vec(i);
					ss.str("");
					ss << "lost connection to master: " << sock_name[0] << ":" << sock_name[1] << endl;
					report(ss.str(), true);
					w_close(i); // bye!
					FD_CLR(i, &master); // remove from master set
					err.first = -999;
					err.second = "lost connection to master";
					return err;
				}
				else
				{
					// received data sored in net_pack return to calling routine to process it
					err.first = 1;
					err.second = "successful receive from master";
					return err;
				}
			}
		}
	}
	ss.str("");
	ss << "recv from master failed " << max_recv_fails << " times, exiting..." << endl;
	report(ss.str(), true);
	return err;
	// returns -1  receive error
	//         -990  error in call to select()
	//         -991  connection closed
	//          1  message received
	//          2  no message received
}

pair<int,string> PANTHERAgent::recv_message(NetPackage &net_pack, long  timeout_seconds, long  timeout_microsecs)
{
	pair<int, string> err;
	err.first = -1;
	int result = 0;
	struct timeval tv;
	tv.tv_sec = timeout_seconds;
	tv.tv_usec = timeout_microsecs;
	err = recv_message(net_pack, &tv);
	return err;
}

void PANTHERAgent::transfer_files(const vector<string>& tfiles, int group, int run_id, string& desc, string tag)
{
    int bytes_read;

    stringstream ss;
    NetPackage pack;
    for (auto& filename :tfiles) {
        if (!check_exist_in(filename))
        {
            ss.str("");
            ss << "file " << filename << " does not exists";
            report(ss.str(), true);
            continue;
        }

        ifstream in;
        in.open(filename.c_str(), ifstream::binary);
        if (in.bad()) {
            ss.str("");
            ss << "error opening file " << filename << " for reading";
            report(ss.str(), true);
            in.close();
            continue;
        }
        string filename_desc = desc + " agent_filename:"+filename + " " + tag;

        pack = NetPackage(NetPackage::PackType::START_FILE_WRKR2MSTR,group,run_id,filename_desc);
        send_message(pack);
        report("starting file transfer of '" + filename + "' for info '" + filename_desc + "' ",true);
        int total_size = 0;
        pack = NetPackage(NetPackage::PackType::CONT_FILE_WRKR2MSTR,group,run_id,filename_desc);
        char buf[NetPackage::FILE_TRANS_BUF_SIZE]={'\0'};
        in.seekg(0,in.end);

        int file_size = in.tellg();
        in.seekg(0,in.beg);
        if (file_size > NetPackage::FILE_TRANS_BUF_SIZE) {
            while ((total_size < file_size) && (in.read(buf, sizeof(buf)))) {
                //cout << filename << ": " << buf << endl;
                send_message(pack, buf, sizeof(buf));
                total_size = total_size + sizeof(buf);
                if ((file_size - total_size) < NetPackage::FILE_TRANS_BUF_SIZE)
                    break;
                
            }
        }
        if (total_size < file_size) {
            char buf[NetPackage::FILE_TRANS_BUF_SIZE]={'\0'};
            in.read(buf,file_size-total_size);
            //cout << filename << "fsize: " << file_size << ", total_size: " << total_size << ", final: " << string(buf,file_size-total_size) << endl;
            send_message(pack, buf, file_size-total_size);
            total_size = total_size + (file_size-total_size);
        }
        
        pack = NetPackage(NetPackage::PackType::FINISH_FILE_WRKR2MSTR,group,run_id,filename_desc);
        send_message(pack);
        ss.str("");
        ss << "sent " << total_size << " bytes for file '" << filename << "', file size: " << file_size;
        report(ss.str(),true);
        in.close();
    }
}


pair<int,string> PANTHERAgent::send_message(NetPackage &net_pack, const void *data, unsigned long data_len)
{
	pair<int,string> err;
	int n;
	stringstream ss;
	for (err.first = -1, n = 0; err.first != 1 && n < max_send_fails; ++n)
	{
		err = net_pack.send(sockfd, data, data_len);
		if (err.first <= 0)
		{
			ss.str("");
			ss << "failed to send to master: " << err.second << ", trying again..." << endl;
			report(ss.str(), true);
		}
	}
	if (n >= max_send_fails)
	{
		ss.str("");
		ss << "send to master failed " << max_send_fails << " times, giving up..." << endl;
		report(ss.str(), true);
	}
	return err;
}


std::pair<NetPackage::PackType,std::string> PANTHERAgent::run_model(Parameters &pars, Observations &obs, NetPackage &net_pack)
{
	NetPackage::PackType final_run_status = NetPackage::PackType::RUN_FAILED;
	stringstream smessage;
	bool done = false;
	pair<int, string> err;
	err.first = 0;

	//if (!mi.get_initialized())
	//{
	//	//initialize the model interface
	//	mi.initialize(tplfile_vec, inpfile_vec, insfile_vec,
	//		outfile_vec, comline_vec, par_name_vec, obs_name_vec);
	//}

	thread_flag f_terminate(false);
	thread_flag f_finished(false);
	exception_ptr run_exception = nullptr;
	stringstream ss;
    vector<string> par_name_vec;
    vector<double> par_values;
    for (auto &i : pars)
    {
        par_name_vec.push_back(i.first);
        par_values.push_back(i.second);
    }

    vector<double> obs_vec;
    thread run_thread(&PANTHERAgent::run_async, this, &f_terminate, &f_finished, std::ref(run_exception),
                      &pars, &obs);

    try
	{

		while (true)
		{

			if (run_exception)
			{
                try {
                    rethrow_exception(run_exception);
                }
                catch (exception& e) {
                    ss.str("");
                    ss << "exception raised by run thread: " << std::endl;
                    ss << e.what() << std::endl;
                    report(ss.str(), true);
                    //don't break here, need to check one last time for incoming messages
                    done = true;
                }
			}
			//check if the runner thread has finished
			if (f_finished.get())
			{
				ss.str("");
				ss << "received finished signal from run thread " << std::endl;
				report(ss.str(), true);
				//don't break here, need to check one last time for incoming messages
				done = true;
			}
			//this call includes a "sleep" for the timeout
			err = recv_message(net_pack, 0, 100000);
			if (err.first < 0)
			{
				ss.str("");
				ss << "error receiving message from master: " << err.second << endl;
				report(ss.str(), true);
				f_terminate.set(true);
				terminate_or_restart(-1);
			}
			//timeout on recv
			else if (err.first == 2)
			{
			}
			else if (net_pack.get_type() == NetPackage::PackType::PING)
			{
				//cout << "ping request received...";
				net_pack.reset(NetPackage::PackType::PING, 0, 0, "");
				const char* data = "\0";
				report("sending ping response to master",false);
				err = send_message(net_pack, &data, 0);
				if (err.first != 1)
				{
					ss.str("");
					ss << "Error sending ping response to master: " << err.second << "...quitting" << endl;
					report(ss.str(), true);
					f_terminate.set(true);
					terminate_or_restart(-1);

					smessage << "Error sending ping response to master...quitting";
				}
				//cout << "ping response sent" << endl;
			}
			else if (net_pack.get_type() == NetPackage::PackType::REQ_KILL)
			{
				ss.str("");
				ss << "received kill request signal from master, ";
				ss << "sending terminate signal to run thread" << endl;
				report(ss.str(), true);
				f_terminate.set(true);
				final_run_status = NetPackage::PackType::RUN_KILLED;
				smessage << "received kill request signal from master";
				break;
			}
			else if (net_pack.get_type() == NetPackage::PackType::TERMINATE)
			{
				ss.str("");
				ss << "received terminate signal from master, ";
				ss << "sending terminate signal to run thread" << endl;
				report(ss.str(), true);
				f_terminate.set(true);
				terminate = true;
				final_run_status = NetPackage::PackType::TERMINATE;
				break;
			}
			else
			{
				ss.str("");
				ss << "Received unsupported message from master, only PING REQ_KILL or TERMINATE can be sent during model run, not: ";
				ss << net_pack.pack_strings[static_cast<int>(net_pack.get_type())] <<  " run_id:" << net_pack.get_run_id();
				report(ss.str(), true);
				f_terminate.set(true);
				final_run_status = NetPackage::PackType::CORRUPT_MESG;
				smessage << "Received unsupported message from master, only PING REQ_KILL or TERMINATE can be sent during model run, not:";
				smessage << net_pack.pack_strings[static_cast<int>(net_pack.get_type())] << " run_id:" << net_pack.get_run_id();
				//terminate_or_restart(-1);
				break;
			}
			if (done)
			{
				break;
			}
		}

		if (!f_terminate.get())
		{
			final_run_status = NetPackage::PackType::RUN_FINISHED;
		}
        run_thread.join();
	}

	catch(const PANTHERAgentRestartError&)
	{
		// Rethrow for start() method to handle
		throw;
	}
	catch(const std::exception& ex)
	{
		ss.str("");
	
		ss << "error(s) thrown during async run: " << ex.what() << " ";
		ss << "Aborting model run";
		report(ss.str(), true);
		smessage << ex.what();
		final_run_status = NetPackage::PackType::RUN_FAILED;

		
	}
	catch(...)
	{
		ss.str("");
 		ss << "   Error running model, aborting model run";
		report(ss.str(), true);
		final_run_status = NetPackage::PackType::RUN_FAILED;
	}

    if (run_exception)
    {
        try {
            rethrow_exception(run_exception);
        }
        catch (exception& e) {
            ss.str("");
            ss << "exception raised by run thread: " << std::endl;
            ss << e.what() << std::endl;
            report(ss.str(), true);
            final_run_status = NetPackage::PackType::RUN_FAILED;
        }
    }

	//sleep here just to give the os a chance to cleanup any remaining file handles
	w_sleep(poll_interval_seconds * 1000);
	return pair<NetPackage::PackType,std::string> (final_run_status,smessage.str());
}


void PANTHERAgent::run_async(pest_utils::thread_flag* terminate, pest_utils::thread_flag* finished, exception_ptr& run_exception,
	Parameters* pars, Observations* obs)
{
    mi.set_sleep_ms(5);
	mi.run(terminate,finished,run_exception, pars, obs);
}


void PANTHERAgent::start(const string &host, const string &port)
{
	stringstream ss;
	if(restart_on_error)
	{
		ss.str("");

		ss << "PANTHER worker will restart on any communication error.";
		report(ss.str(), true);
	}

	do
	{
		try
		{
			// Start the agent
			start_impl(host, port);

			// If we make it this far, there was no error, so we do not need to restart
			return;
		}
		catch(const PANTHERAgentRestartError& ex)
		{
			// A fatal comms error occurred; wait a bit and then restart
			this_thread::sleep_for(chrono::seconds(5));
			ss.str("");
			
			ss << "Restarting PANTHER worker...";
			report(ss.str(), true);
		}
	} while(restart_on_error);
}


void PANTHERAgent::start_impl(const string &host, const string &port)
{
	NetPackage net_pack;
	Observations obs = pest_scenario.get_ctl_observations();
	Parameters pars;
	vector<int8_t> serialized_data;
	pair<int,string> err;
	vector<string> par_name_vec, obs_name_vec;
	stringstream ss;
	//class attribute - can be modified in run_model()
	terminate = false;
	init_network(host, port);

	std::chrono::system_clock::time_point last_ping_time = chrono::system_clock::now();
	while (!terminate)
	{
		//get message from master
		err = recv_message(net_pack, recv_timeout_secs, 0);

		// Refresh ping timer
		if (err.first != 2)
		{
			last_ping_time = chrono::system_clock::now();
		}

		if (err.first == -999)
		{
			ss.str("");
			ss << "error receiving message from master: " << err.second << " , terminating";
			report(ss.str(), true);
			//terminate = true;
			net_pack.reset(NetPackage::PackType::CORRUPT_MESG, 0, 0, "recv security message error");
			char data;
			err = send_message(net_pack, &data, 0);
			//if (err != 1)
			//{
			//	terminate_or_restart(-1);
			//}

			terminate_or_restart(-1);
			}
			
		if (err.first < 0)
		{
			ss.str("");
			ss << "error receiving message from master: " << err.second << " terminating, header follows: ";
			net_pack.print_header(ss);
			report(ss.str(), true);
			//terminate = true;
			terminate_or_restart(-1);
		}
		else if (err.first == 2)
		{
			// Timeout on socket receive
			// Optionally: die if no data received from master for a long time (e.g. 5 minutes)
			// Set max_time_without_master_ping_seconds <= 0 to disable and wait forever
			auto time_without_master_ping_seconds = chrono::duration_cast<std::chrono::seconds>(chrono::system_clock::now() - last_ping_time).count();
			if (max_time_without_master_ping_seconds > 0 && time_without_master_ping_seconds > max_time_without_master_ping_seconds)
			{
				time_t rawtime = std::chrono::system_clock::to_time_t(last_ping_time);
				struct tm* timeinfo;
				char buffer[80];
				time(&rawtime);
				timeinfo = localtime(&rawtime);
				strftime(buffer, 80, "%m/%d/%y %H:%M:%S", timeinfo);
				string t_str(buffer);
				ss.str("");
				ss << "no ping received from master in the last " << max_time_without_master_ping_seconds << " seconds, last ping time: " <<  t_str << ",  terminating";
				report(ss.str(), true);
				net_pack.reset(NetPackage::PackType::CORRUPT_MESG, 0, 0, "ping overdue from master, exiting");
				send_message(net_pack);
				//terminate = true;
				terminate_or_restart(-1);
			}
		}
		else if(net_pack.get_type() == NetPackage::PackType::REQ_RUNDIR)
		{
			// Send Master the local run directory.  This information is only used by the master
			// for reporting purposes
			report("responding to REQ_RUNDIR", true);
			net_pack.reset(NetPackage::PackType::RUNDIR, 0, 0,"");
			string cwd =  OperSys::getcwd();
			err = send_message(net_pack, cwd.c_str(), cwd.size());
			if (err.first != 1)
			{
				ss.str("");
				ss << "error sending RUNDIR message to master: " << err.second << ", terminating";
				report(ss.str(), true);
				terminate_or_restart(-1);
			}
		}
		else if (net_pack.get_type() == NetPackage::PackType::PAR_NAMES)
		{
			report("received PAR_NAMES", true);
			//Don't check first8 bytes as these contain an integer which stores the size of the data.
			bool safe_data = NetPackage::check_string(net_pack.get_data(), 0, net_pack.get_data().size());
			if (!safe_data)
			{
				ss.str("");
				ss << "received corrupt parameter name packet from master,";
				ss << "terminating execution ..." << endl << endl;
				report(ss.str(), true);
				net_pack.reset(NetPackage::PackType::CORRUPT_MESG, 0, 0, "");
				char data;
				pair<int,string> np_err = send_message(net_pack, &data, 0);
				terminate_or_restart(-1);
			}
			Serialization::unserialize(net_pack.get_data(), par_name_vec);
			//make sure all par names are found in the scenario
			//vector<string> vnames = pest_scenario.get_ctl_ordered_par_names();
			vector<string> vnames = pest_scenario.get_ctl_parameters().get_keys();
			set<string> snames(vnames.begin(), vnames.end());
			vnames.clear();
			for (auto& pname : par_name_vec)
			{
				pest_utils::upper_ip(pname);
				if (snames.find(pname) == snames.end())
					vnames.push_back(pname);
			}
			if (vnames.size() > 0)
			{
				ss.str("");
				ss << vnames.size() << " par names not found in agent ctl file: ";
				report(ss.str(), true);
				
				for (auto vname : vnames)
				{
					ss << vname << ",";
					frec << vname << endl;
				}
				net_pack.reset(NetPackage::PackType::IO_ERROR, 0, 0, ss.str());
				err = send_message(net_pack);
				w_close(sockfd);
				w_cleanup();
				exit(-1);
			}	
			snames.clear();
			snames.insert(par_name_vec.begin(), par_name_vec.end());
			for (auto name : pest_scenario.get_ctl_parameters().get_keys())
			{
				if (snames.find(name) == snames.end())
					vnames.push_back(name);
			}
			if (vnames.size() > 0)
			{
				ss.str("");
				ss << "ERROR:" << vnames.size() << " extra par names found in agent ctl file: ";
				report(ss.str(), true);
				for (auto vname : vnames)
				{
					ss << vname << ",";
					frec << vname << endl;
				}
				net_pack.reset(NetPackage::PackType::IO_ERROR, 0, 0, ss.str());
				err = send_message(net_pack);
				w_close(sockfd);
				w_cleanup();
				exit(-1);
			}

		}
		else if (net_pack.get_type() == NetPackage::PackType::OBS_NAMES)
		{
			report("received OBS_NAMES", true);
			//Don't check first8 bytes as these contain an integer which stores the size of the data.
			bool safe_data = NetPackage::check_string(net_pack.get_data(), 0, net_pack.get_data().size());
			if (!safe_data)
			{
				ss.str("");
				ss << "received corrupt observation name packet from master," ;
				ss << "received corrupt observation name packet from master" << endl;
				report(ss.str(), true);
				net_pack.reset(NetPackage::PackType::CORRUPT_MESG, 0, 0, "");
				char data;
				pair<int,string> np_err = send_message(net_pack, &data, 0);
				terminate_or_restart(-1);
			}
			Serialization::unserialize(net_pack.get_data(), obs_name_vec);
			//make sure all par names are found in the scenario
			//vector<string> vnames = pest_scenario.get_ctl_ordered_obs_names();
            vector<string> vnames = pest_scenario.get_ctl_observations().get_keys();
			set<string> snames(vnames.begin(), vnames.end());
			vnames.clear();
			for (auto& oname : obs_name_vec)
			{
				pest_utils::upper_ip(oname);
				if (snames.find(oname) == snames.end())
					vnames.push_back(oname);
			}
			if (vnames.size() > 0)
			{
				ss.str("");
				ss << "ERROR: " << vnames.size() << " obs names not found in agent ctl file: ";
				report(ss.str(), true);
				for (auto vname : vnames)
				{
					ss << vname << ",";
					frec << vname << endl;
				}
				for (auto vname : vnames)
					ss << vname << ",";
				net_pack.reset(NetPackage::PackType::IO_ERROR, 0, 0, ss.str());
				err = send_message(net_pack);
				w_close(sockfd);
				w_cleanup();
				exit(-1);
			}
			snames.clear();
			snames.insert(obs_name_vec.begin(), obs_name_vec.end());
			for (auto name : pest_scenario.get_ctl_observations().get_keys())
			{
				if (snames.find(name) == snames.end())
					vnames.push_back(name);
			}
			if (vnames.size() > 0)
			{
				ss.str("");
				ss << "ERROR: " << vnames.size() << " extra obs names found in agent ctl file: ";
				report(ss.str(), true);
				for (auto vname : vnames)
				{
					ss << vname << ",";
					frec << vname << endl;
				}
				net_pack.reset(NetPackage::PackType::IO_ERROR, 0, 0, ss.str());
				err = send_message(net_pack);
				w_close(sockfd);
				w_cleanup();
				exit(-1);
			}
		}
		else if(net_pack.get_type() == NetPackage::PackType::REQ_LINPACK)
		{
			report("received REQ_LINPACK",true);
			//linpack_wrap();
			net_pack.reset(NetPackage::PackType::LINPACK, 0, 0,"");
			char data;
			err = send_message(net_pack, &data, 0);
			if (err.first != 1)
			{
				ss.str("");
				ss << "error sending LINPACK message to master: " << err.second << ", terminating" << endl;
				report(ss.str(), true);
				terminate_or_restart(-1);
			}
		}
		else if(net_pack.get_type() == NetPackage::PackType::START_RUN)
		{
			
			int group_id = net_pack.get_group_id();
			int run_id = net_pack.get_run_id();
			string info_txt = net_pack.get_info_txt();
			pest_utils::upper_ip(info_txt);
			int da_cycle = NetPackage::NULL_DA_CYCLE;
			if (info_txt.find("DA_CYCLE:") != string::npos)
			{
				frec << "Note: 'DA_CYCLE' information passed in START_RUN command" << endl;
				frec << "      info txt for group_id:run_id " << group_id << ":" << run_id << endl;
				cout << "Note: 'DA_CYCLE' information passed in START_RUN command" << endl;
				cout << "      info txt for group_id:run_id " << group_id << ":" << run_id << endl;
				vector<string> tokens,ttokens;
				pest_utils::tokenize(info_txt, tokens, " ");
				
				for (auto token : tokens)
				{
					if (token.find(":") != string::npos)
					{
						pest_utils::tokenize(token, ttokens, ":");
						if (ttokens[0] == "DA_CYCLE")
						{
							if (ttokens[1].size() > 0)
							{
								string s_cycle = ttokens[1];
								try
								{
									da_cycle = stoi(s_cycle);
								}
								catch (...)
								{
									frec << "WARNING: error casting '" + ttokens[1] + "' to int for da_cycle...continuing" << endl;
									frec << "WARNING: error casting '" + ttokens[1] + "' to int for da_cycle...continuing" << endl;

								}
							}
						}
					}
				}
				if (da_cycle != NetPackage::NULL_DA_CYCLE)
				{
					if (da_cycle != current_da_cycle)
					{
						try
						{
							pest_scenario.assign_da_cycles(frec);
							Pest childPest = pest_scenario.get_child_pest(da_cycle);
							const ParamTransformSeq& base_trans_seq = childPest.get_base_par_tran_seq();
							Parameters cur_ctl_parameters = childPest.get_ctl_parameters();
							vector<string> par_names = base_trans_seq.ctl2model_cp(cur_ctl_parameters).get_keys();
							sort(par_names.begin(), par_names.end());
							vector<string> obs_names = childPest.get_ctl_observations().get_keys();
							sort(obs_names.begin(), obs_names.end());
							par_name_vec = par_names;
							obs_name_vec = obs_names;
							mi = ModelInterface(childPest.get_tplfile_vec(), childPest.get_inpfile_vec(),
								childPest.get_insfile_vec(), childPest.get_outfile_vec(), childPest.get_comline_vec());
							obs = childPest.get_ctl_observations();
							stringstream ss;
							ss << "Updated interface components for DA_CYCLE " << da_cycle << " as follows: " << endl;
							report(ss.str(), true);
							int i = 0;
							ss.str("");
							ss << "parameter names:" << endl;
							for (int i = 0; i < par_name_vec.size(); i++)
							{
								ss << par_name_vec[i] << " ";
								if (i % 10 == 0)
									ss << endl;
							}
							frec << ss.str() << endl;
							//cout << ss.str() << endl;
							
							cout << par_name_vec.size() << " parameters in current cycle, see rec file for listing" << endl;
							ss.str("");
							ss << endl << "observation names:" << endl;
							for (int i = 0; i < obs_name_vec.size(); i++)
							{
								ss << obs_name_vec[i] << " ";
								if (i % 10 == 0)
									ss << endl;
							}
							frec << ss.str() << endl;
							//cout << ss.str() << endl;
							cout << obs_name_vec.size() << " observations in current cycle, see rec file for listing" << endl;

							ss.str("");
							ss << endl << "tpl:in file names:" << endl;
							vector<string> tpl_vec = childPest.get_tplfile_vec();
							vector<string> in_vec = childPest.get_inpfile_vec();
							for (int i = 0; i < tpl_vec.size(); i++)
							{
								ss << tpl_vec[i] << ":" << in_vec[i] << " ";
								if (i % 5 == 0)
									ss << endl;
							}
							frec << ss.str() << endl;
							//cout << ss.str() << endl;
							cout << tpl_vec.size() << " template files in current cycle, see rec file for listing" << endl;

							ss.str("");
							ss << endl << "ins:out file names:" << endl;
							vector<string> ins_vec = childPest.get_insfile_vec();
							vector<string> out_vec = childPest.get_outfile_vec();
							for (int i = 0; i < ins_vec.size(); i++)
							{
								ss << ins_vec[i] << ":" << out_vec[i] << " ";
								if (i % 5 == 0)
									ss << endl;
							}
							frec << ss.str() << endl << endl;
							//cout << ss.str() << endl << endl;
							cout << ins_vec.size() << " instruction files in current cycle, see rec file for listing" << endl;

							current_da_cycle = da_cycle;
						}
						catch (exception& e)
						{
							stringstream ss;
							ss << "ERROR: could not process 'DA_CYCLE' " << da_cycle << ": " << e.what();
							frec << ss.str() << endl;
							cout << ss.str() << endl;
							net_pack.reset(NetPackage::PackType::RUN_FAILED, group_id, run_id, ss.str());
							char data;
							err = send_message(net_pack, &data, 0);
							terminate = true;
							continue;
						}
						catch (...)
						{
							stringstream ss;
							ss << "ERROR: could not process 'DA_CYCLE' " << da_cycle;
							frec << ss.str() << endl;
							cout << ss.str() << endl;
							net_pack.reset(NetPackage::PackType::RUN_FAILED, group_id, run_id, ss.str());
							char data;
							err = send_message(net_pack, &data, 0);
							terminate = true;
							continue;
						}
					}
					else
					{
						stringstream ss;
						ss << "reusing 'DA_CYCLE' " << da_cycle;
						frec << ss.str() << endl;
						cout << ss.str() << endl;

					}
				}
			}
			
			//do this after we handle a cycle change so that par_name_vec is updated
			Serialization::unserialize(net_pack.get_data(), pars, par_name_vec);

			/*frec << "parameters for run_id: " << run_id << ", group_id: " << group_id;
			if (da_cycle != NetPackage::NULL_DA_CYCLE)
				frec << ", da_cycle: " << da_cycle;
			frec << endl << "name value " << endl;

			for (auto name : par_name_vec)
			{
				frec << name << " " << pars.get_rec(name) << endl;
			}*/

			// run model
			if (pest_scenario.get_pestpp_options().get_panther_debug_loop())
			{
				ss.str("");
				ss << "PANTHER_DEBUG_LOOP = true, returning ctl obs values";
				
				report(ss.str(), true);
				
				serialized_data = Serialization::serialize(pars, par_name_vec, obs, obs_name_vec, run_time);
				ss.str("");
				double rd = ((double)rand() / (double)RAND_MAX);
				if (rd < 0.1)
				{
					ss << "debug loop returning failed run for: " << run_id << "," << group_id;
					net_pack.reset(NetPackage::PackType::RUN_FAILED, group_id, run_id, ss.str());
				}
				else
				{
					ss << "debug loop returning ctl obs for run_id, group_id: " << run_id << "," << group_id;
					net_pack.reset(NetPackage::PackType::RUN_FINISHED, group_id, run_id, ss.str());
				}
				err = send_message(net_pack, serialized_data.data(), serialized_data.size());
				if (err.first != 1)
				{
					ss.str("");
					ss << "error sending RUN_FINISHED message to master: " << err.second << ", terminating";
					report(ss.str(), true);
					terminate_or_restart(-1);
				}
				
				ss.str("");
				ss << "results of run_id " << run_id << " sent successfully";
				report(ss.str(), true);
				ss.str("");
				ss << "sending ready signal to master";
				report(ss.str(), true);
				net_pack.reset(NetPackage::PackType::READY, 0, 0, "lets do it");
				char data;
				err = send_message(net_pack, &data, 0);
				if (err.first != 1)
				{
					ss.str("");
					ss << "error sending READY message to master: " << err.second << ", terminating";
					report(ss.str(), true);
					terminate_or_restart(-1);
				}
				continue;

			}
			ss.str("");
			ss << "received parameters ( group_id=" << group_id << ", run_id=" << run_id << ", info_txt=" << info_txt << " ), ";
			ss << "starting model run..." << endl;
			report(ss.str(), true);

			try
			{
				ofstream fout("run.info");
				if (fout.good())
				{
					fout << "run_id, " << run_id << endl;
					fout << "group_id, " << group_id << endl;
					fout << "info_txt," << info_txt << endl;
				}
				fout.close();
			}
			catch (...)
			{
				
			}

			std::chrono::system_clock::time_point start_time = chrono::system_clock::now();
			pair<NetPackage::PackType,std::string> final_run_status = run_model(pars, obs, net_pack);
            double run_time = pest_utils::get_duration_sec(start_time);
			if (final_run_status.first == NetPackage::PackType::RUN_FINISHED)
			{

				//send model results back
				ss.str("");
				ss << "run complete, ";
				ss << "sending results to master (group_id=" << group_id << " , run_id=" << run_id << " , info_txt=" << info_txt <<" )...";
				ss << "run took: " << run_time << " seconds";
				report(ss.str(), true);
				ss.str("");
				ss << " worker_time:" << run_time/60.0;
				string message = info_txt + " " + final_run_status.second + ss.str();
				serialized_data = Serialization::serialize(pars, par_name_vec, obs, obs_name_vec, run_time);
				net_pack.reset(NetPackage::PackType::RUN_FINISHED, group_id, run_id, message);
				err = send_message(net_pack, serialized_data.data(), serialized_data.size());
				if (err.first != 1)
				{
					ss.str("");
					ss << "error sending RUN_FINISHED message to master: " << err.second << " , terminating";
					report(ss.str(), true);
					terminate_or_restart(-1);
				}
				ss.str("");
				ss << "results of run_id " << run_id << " sent successfully";
				report(ss.str(), true);
                transfer_files(pest_scenario.get_pestpp_options().get_panther_transfer_on_finish(), group_id,
                               run_id,info_txt, "run_status:completed");

            }
			else if (final_run_status.first == NetPackage::PackType::RUN_FAILED)
			{
				ss.str("");
				ss << "run failed for run_id: " << run_id << " " << info_txt << "  " << final_run_status.second;
				report(ss.str(), true);
				ss.str("");
				ss << "group_id:" << group_id << " run_id:" << run_id << " " << info_txt << " " << final_run_status.second;
				net_pack.reset(NetPackage::PackType::RUN_FAILED, group_id, run_id,ss.str());
				char data;
				err = send_message(net_pack, &data, 0);
				if (err.first != 1)
				{
					ss.str("");
					ss << "error sending RUN_FAILED message to master: " << err.second << " terminating";
					report(ss.str(), true);
					terminate_or_restart(-1);
				}
                transfer_files(pest_scenario.get_pestpp_options().get_panther_transfer_on_fail(), group_id,
                               run_id,info_txt, "run_status:failed");
				if (pest_scenario.get_pestpp_options().get_panther_debug_fail_freeze())
				{
					ss.str("");
					ss << "debug_panther_fail_freeze = true, entering frozen state...";
					report(ss.str(), true);
					net_pack.reset(NetPackage::PackType::DEBUG_FAIL_FREEZE, group_id, run_id, final_run_status.second);
					char data;
					err = send_message(net_pack, &data, 0);
					if (err.first != 1)
					{
						ss.str("");
						ss << "error sending DEBUG_FAIL_FREEZE message to master: " << err.second << "...freezing anyway";
						report(ss.str(), true);
					}
					while (true)
					{
						ss.str("");
						ss << "frozen";
						report(ss.str(), true);
						w_sleep(10000);
//						if (quit_file_found()) {
//                            report("pest.stp file found, resetting panther_agent_freeze_on_fail and continuing...",
//                                   true);
//                            pest_scenario.get_pestpp_options_ptr()->set_panther_debug_fail_freeze(false);
//                            break;
//                        }
					}
				}
			}
			else if (final_run_status.first == NetPackage::PackType::RUN_KILLED)
			{
				ss.str("");
				ss << "run_id:" << run_id << " " << info_txt << " killed";
				report(ss.str(), true);
				net_pack.reset(NetPackage::PackType::RUN_KILLED, group_id, run_id, final_run_status.second);
				char data;
				err = send_message(net_pack, &data, 0);
				if (err.first != 1)
				{
					ss.str("");
					ss << "error sending RUN_KILLED message to master: " << err.second << " terminating";
					report(ss.str(), true);
					terminate_or_restart(-1);
				}
                transfer_files(pest_scenario.get_pestpp_options().get_panther_transfer_on_fail(), group_id,
                               run_id,info_txt, "run_status:failed");
			}

			else if (final_run_status.first == NetPackage::PackType::CORRUPT_MESG)
			{
				ss << "corrupt/incorrect message received from master:" << final_run_status.second << " " << info_txt << " quitting for safety";
				net_pack.reset(NetPackage::PackType::RUN_KILLED, group_id, run_id, ss.str());
				char data;
				err = send_message(net_pack, &data, 0);
				if (err.first != 1)
				{
					ss.str("");
					ss << "error sending CORRUPT_MESG message to master:" << err.second << " terminating";
					report(ss.str(), true);
					terminate_or_restart(-1);
				}
				report(ss.str(), true);
				terminate_or_restart(-1);

			}
			else if (final_run_status.first == NetPackage::PackType::TERMINATE)
			{
				ss.str("");
				ss << "run preempted by termination requested";
				report(ss.str(), true);
				terminate = true;
			}

			if (!terminate)
			{
				// Send READY Message to master
				ss.str("");
				ss << "sending ready signal to master";
				report(ss.str(), true);
				net_pack.reset(NetPackage::PackType::READY, 0, 0, final_run_status.second);
				char data;
				err = send_message(net_pack, &data, 0);
				if (err.first != 1)
				{
					ss.str("");
					ss << "error sending READY message to master:" << err.second << " terminating";
					report(ss.str(), true);
					terminate_or_restart(-1);
				}
			}
		}
		else if (net_pack.get_type() == NetPackage::PackType::TERMINATE)
		{
			ss.str("");
			ss << "terminate requested" << endl;
			report(ss.str(), true);
			terminate = true;
		}
		else if (net_pack.get_type() == NetPackage::PackType::REQ_KILL)
		{
			ss.str("");
			ss << "received kill request from master. run already finished";
			report(ss.str(), true);
			
		}
		else if (net_pack.get_type() == NetPackage::PackType::PING)
		{
			ss.str("");
			ss << "ping request received...";
			report(ss.str(), true);
			net_pack.reset(NetPackage::PackType::PING, 0, 0, "");
			const char* data = "\0";
			err = send_message(net_pack, &data, 0);
			if (err.first != 1)
			{
				ss.str("");
				ss << "error sending PING message to master:" << err.second << " terminating";
				report(ss.str(), true);
				terminate_or_restart(-1);
			}
		
			report("ping response sent",true);
		}
		else
		{
			ss.str("");

			ss << "received unsupported messaged type: " << int(net_pack.get_type()) << ": " << net_pack.get_info_txt();
			report(ss.str(), true);
		}
		//w_sleep(100);
		this_thread::sleep_for(chrono::milliseconds(100));
	}
}


void PANTHERAgent::terminate_or_restart(int error_code) const
{
	
	w_sleep(poll_interval_seconds * 10000);
	if(!restart_on_error)
	{
		exit(error_code);
	}
	
	// Cleanup sockets and throw exception to signal restart
	w_close(sockfd);
	w_cleanup();
	throw PANTHERAgentRestartError("");
}

void PANTHERAgent::report(const string& message, bool to_cout)
{
	string t_str = pest_utils::get_time_string();
	frec << t_str << "->" << message << endl;
	if (to_cout) cout << endl << t_str << "->" << message << endl;
}
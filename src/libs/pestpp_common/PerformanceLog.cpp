#include <iomanip>
#include <iostream>
#include <string>
#include <algorithm>
#include <sstream>
#include <cstring>
#include "PerformanceLog.h"
#include "config_os.h"

using namespace std;
using std::chrono::system_clock;

void PerformanceLog::writetime(stringstream &os, time_t tc) {
	// alternative to put_time iomanip
        // as put_time is not implemented in gcc4.8
	locale loc;
	const time_put<char>& tp = use_facet<time_put<char>>(loc);
	const char *pat = "%F %T";
	tp.put(os,os,' ',localtime(&tc),pat,pat+strlen(pat));
}

PerformanceLog::PerformanceLog(ofstream &_fout)
: fout(_fout)
{
	prev_time = system_clock::now();
	//fout << "PEST++ performance logger started at:  " << time_to_string(prev_time) << endl;
	fout << "time,elapsed_seconds,message" << endl;
	log_event("PEST++ performance logger started");

}

void PerformanceLog::log_event(const string &message)
{
	system_clock::time_point time_now = system_clock::now();
	string elapsed_str = elapsed_time_to_string(time_now, prev_time);
	fout << time_to_string(time_now) << "," << elapsed_str << "," << message << endl;
	prev_time = time_now;
	fout.flush();
}

string PerformanceLog::time_to_string(const std::chrono::system_clock::time_point &tmp_time)
{
	stringstream time_str;
	auto tmp_time_c = system_clock::to_time_t(tmp_time);
	#ifdef OS_LINUX
	writetime(time_str, tmp_time_c);
	#endif
	#ifdef OS_WIN
	time_str <<  put_time(std::localtime(&tmp_time_c), "%X");
	#endif
	return time_str.str();
}

string PerformanceLog::elapsed_time_to_string(std::chrono::system_clock::time_point &current_time, std::chrono::system_clock::time_point &prev_time)
{
	ostringstream ss;
	std::chrono::duration<double, std::chrono::seconds::period> double_delta_t(current_time - prev_time);
	auto delta_t = current_time - prev_time;
	//double_delta_t = std::chrono::duration_cast<std::chrono::seconds>(delta_t);
	ss << double_delta_t.count();
	return ss.str();
	/*if (delta_t < std::chrono::milliseconds(1))
	{
		str << std::chrono::duration_cast<std::chrono::microseconds>(delta_t).count() << "us";
	}
	else if (delta_t < std::chrono::seconds(1))
	{
		str << std::chrono::duration_cast<std::chrono::milliseconds>(delta_t).count() << "ms";
	}
	else if (delta_t < std::chrono::minutes(1))
	{
		str << std::chrono::duration_cast<std::chrono::seconds>(delta_t).count() << "sec";
	}
	else if (delta_t < std::chrono::hours(1))
	{
		str << std::chrono::duration_cast<std::chrono::minutes>(delta_t).count() << "min";
	}
	else
	{
		str << std::chrono::duration_cast<std::chrono::hours>(delta_t).count() << "hr";
	}
	return str.str();*/
}


PerformanceLog::~PerformanceLog()
{
}

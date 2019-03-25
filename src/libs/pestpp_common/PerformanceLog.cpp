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
        // as put_time is not implimented in gcc4.8
	locale loc;
	const time_put<char>& tp = use_facet<time_put<char>>(loc);
	const char *pat = "%F %T";
	tp.put(os,os,' ',localtime(&tc),pat,pat+strlen(pat));
}


PerformanceLog::PerformanceLog(ofstream &_fout)
: fout(_fout), indent_size(2), indent_level(0)
{
	indent_level_prev = 0;
	prev_time = system_clock::now();
	fout << "PEST++ performance logger started at:  " << time_to_string(prev_time) << endl;
}

void PerformanceLog::log_blank_lines(int n)
{
	for (int i = 0; i < n; ++i) fout << endl;
}

void PerformanceLog::add_indent(int n)
{
	indent_level += n;
}
void PerformanceLog::log_event(const string &message, int delta_indent, const std::string &tag)
{
	system_clock::time_point time_now = system_clock::now();
	if (!tag.empty())
	{
		tagged_events[tag] = time_now;
	}

	indent_level = max(0, indent_level+delta_indent);
	fout << string(indent_prev() + 8, ' ') << "( time = " << time_to_string(time_now) << ",  elapsed time = " <<
		elapsed_time_to_string(time_now, prev_time) << " )" << endl;
	fout << string(indent(), ' ') << message << endl;
	prev_time = time_now;
	indent_level_prev = indent_level;
	fout.flush();
}

void PerformanceLog::log_summary(const string &message, const std::string &end_tag, const std::string &begin_tag, int delta_indent)
{
	indent_level = max(0, indent_level + delta_indent);
	fout << string(indent(), ' ') << message << endl;
	fout << string(indent()+8, ' ') << "( elapsed time = " << elapsed_time_to_string(tagged_events[end_tag], tagged_events[begin_tag]) << " )" << endl;
	indent_level_prev = indent_level;
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
	ostringstream str;
	auto delta_t = current_time - prev_time;
	if (delta_t < std::chrono::milliseconds(1))
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
	return str.str();
}


PerformanceLog::~PerformanceLog()
{
}

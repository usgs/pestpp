#ifndef PERFORMANCE_LOG_H_
#define PERFORMANCE_LOG_H_

#include <iostream>
#include <fstream>
#include <chrono>
#include <map>

class PerformanceLog
{
public:
	PerformanceLog(std::ofstream &_fout);
	void log_event(const std::string &message, int delta_indent=0, const std::string &tag="");
	void log_summary(const std::string &message, const std::string &end_tag, const std::string &begin_tag, int delta_indent = 0);
	void log_blank_lines(int n = 1);
	void add_indent(int n = 1);
	~PerformanceLog();
private:
	std::ofstream &fout;
	std::chrono::system_clock::time_point prev_time;
	int indent_level_prev;
	int indent_level;
	int indent_size;
	std::map<std::string, std::chrono::system_clock::time_point> tagged_events;
	int indent() { return indent_level * indent_size; }
	int indent_prev() { return indent_level_prev * indent_size; }
	std::string time_to_string(const std::chrono::system_clock::time_point &tmp_time);
	std::string elapsed_time_to_string(std::chrono::system_clock::time_point &current_time, std::chrono::system_clock::time_point &prev_time);
	void writetime(std::stringstream &os, time_t tc);
};

#endif //PERFORMANCE_LOG_H_

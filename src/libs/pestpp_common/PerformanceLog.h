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
	void log_event(const std::string &message);
	~PerformanceLog();
private:
	std::ofstream &fout;
	std::chrono::system_clock::time_point prev_time;
	std::map<std::string, std::chrono::system_clock::time_point> tagged_events;
	std::string time_to_string(const std::chrono::system_clock::time_point &tmp_time);
	std::string elapsed_time_to_string(std::chrono::system_clock::time_point &current_time, std::chrono::system_clock::time_point &prev_time);
	void writetime(std::stringstream &os, time_t tc);
};

#endif //PERFORMANCE_LOG_H_

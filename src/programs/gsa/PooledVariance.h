#ifndef POOLED_VAR_H_
#define POOLED_VAR_H_

#include <string>
#include <iostream>
#include <fstream>
#include <map>
class PooledVariance
{
public:
	PooledVariance();
	void process_pva_file(std::ifstream &fin);
	~PooledVariance();
private:
	std::map< std::string, std::string> obs_to_pool_grp;
};

#endif /* POOLED_VAR_H_ */

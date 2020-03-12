#ifndef CONSTRAINTS_H_
#define CONSTRAINTS_H_

#include <string>
#include <sstream>
#include <vector>
#include <random>
#include<Eigen/Sparse>

#include "Pest.h"
#include "logger.h"
#include "FileManager.h"
#include "OutputFileWriter.h"
#include "RunManagerAbstract.h"
#include "Transformable.h"
#include "covariance.h"
#include "Jacobian_1to1.h"

using namespace std;

class Constraints
{
public:
	Constraints(Pest* _pest_scenario_ptr, FileManager* _file_manager_ptr, OutputFileWriter* _of_wr_ptr) :
		pest_scenario_ptr(_pest_scenario_ptr), file_manager_ptr(_file_manager_ptr), of_wr_ptr(_of_wr_ptr) {;}
	Constraints() { ; }
	virtual map<string, double> get_residual_map(Observations& sim);
	virtual map<string, double> get_chance_offsets(Observations& sim);
	virtual void add_runs(RunManagerAbstract* run_mgr_ptr);
	void report(Observations& sim);
	void chance_report(Observations& sim);

private:
	Pest* pest_scenario_ptr;
	FileManager* file_manager_ptr;
	OutputFileWriter* of_wr_ptr;
};

class PiConstraints : public Constraints
{
public:
	PiConstraints() { ; }

private:

};

class FosmConstraints : public Constraints
{
public:
	FosmConstraints() { ; }
private:
	Covariance parcov;
	Jacobian_1to1 jco;
};

class StackConstraints : public Constraints
{
public:
	StackConstraints() { ; }
private:

		
};

#endif

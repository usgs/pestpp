#ifndef ENSEMBLESMOOTHER_H_
#define ENSEMBLESMOOTHER_H_

#include <map>
#include <random>
#include <mutex>
#include <thread>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "FileManager.h"
#include "ObjectiveFunc.h"
#include "OutputFileWriter.h"
#include "PerformanceLog.h"
#include "RunStorage.h"
#include "covariance.h"
#include "RunManagerAbstract.h"
#include "ObjectiveFunc.h"
#include "Localizer.h"
#include "EnsembleMethodUtils.h"




class IterEnsembleSmoother: public EnsembleMethod
{
public:

	using EnsembleMethod::EnsembleMethod;
	
	void iterate_2_solution();
	void finalize();
	void throw_ies_error(string message);

private:
	
	void sanity_checks();

};

#endif

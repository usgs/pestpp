#include <vector>
#include <map>
#include <string>

std::vector<double> vec_array_prod(const std::vector<double> &data_vec1, const std::vector<double> &data_vec2, double missing_val);

double vec_mean(const std::vector<double> &data_vec);
double vec_mean_missing_data(const std::vector<double> &data_vec, double missing_val);

std::pair<double, size_t> sum_of_prod_missing_data(const std::vector<double> &x_vec, const std::vector<double> &y_vec, double missing_val);
double sobol_u_missing_data(const std::vector<double> &x_vec, const std::vector<double> &y_vec, double missing_val);

class RunningStats
{
	//Based on paper "Computing the standard deviation efficiently" by Mark Hoemmen, Berkeley
public:
	RunningStats() : n(0), mk(0.0), mk_abs(0.0), qk(0.0) {}
	virtual void reset();
	virtual void add(double sample);
	virtual void add(const std::vector<double> &sample);
	virtual double comp_var() const;
	virtual double comp_sigma()const;
	virtual double comp_mean()const;
	virtual double comp_abs_mean()const;
	virtual long comp_nsamples()const;
	virtual ~RunningStats() {};
private:
	long n;
	double mk;
	double mk_abs;
	double qk;
};

class RunningStatsMissingData : public RunningStats
{
	//Based on paper "Computing the standard deviation efficiently" by Mark Hoemmen, Berkeley
public:
	RunningStatsMissingData(double _missing_value) : RunningStats(), missing_value(_missing_value){}
	virtual void add(double sample);
	virtual void add(const std::vector<double> &sample);
	~RunningStatsMissingData() {};
private:
	double missing_value;
};

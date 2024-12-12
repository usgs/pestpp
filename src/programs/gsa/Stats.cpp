#include <string>
#include <vector>
#include <cmath>
#include <cassert>
#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "Stats.h"

using namespace std;

double vec_mean(const vector<double> &data_vec)
{
	RunningStats stats;
	stats.add(data_vec);
	return stats.comp_mean();
}

double vec_mean_missing_data(const vector<double> &data_vec, double missing_val)
{
	RunningStatsMissingData stats(missing_val);
	stats.add(data_vec);
	return stats.comp_mean();
}

vector<double> vec_array_prod(const vector<double> &data_vec1, const vector<double> &data_vec2, double missing_val)
{
	vector<double> prod_vec = data_vec1;
	assert(data_vec1.size() == data_vec2.size());

	size_t len = data_vec1.size();

	for (size_t i = 0; i < len; ++i)
	{
		if (data_vec1[i] == missing_val || data_vec1[2] == missing_val)
		{
			prod_vec[i] = missing_val;
		}
		else
		{
			prod_vec[i] = data_vec1[i] * data_vec2[i];
		}
	}
	return prod_vec;
}

pair<double, size_t> sum_of_prod_missing_data(const std::vector<double> &x_vec, const std::vector<double> &y_vec, double missing_val)
{
	double m_xy_k = 0;
	double x_k, y_k, xy_k;

	assert(x_vec.size() == y_vec.size());
	size_t n = x_vec.size();
	size_t n_actual = 0;

	int k;
	//compute the mean value of the product of x and y
	for (k = 1; k <= n; ++k)
	{
		x_k = x_vec[k - 1];
		y_k = y_vec[k - 1];
		if (x_k != missing_val && y_k != missing_val)
		{
			xy_k = x_k * y_k;
			m_xy_k = xy_k;
			++n_actual;
			break;
		}
	}
	for (k += 1; k <= n; ++k)
	{
		x_k = x_vec[k - 1];
		y_k = y_vec[k - 1];

		if (x_k != missing_val && x_k != missing_val)
		{
			xy_k = x_k * y_k;
			m_xy_k += (xy_k - m_xy_k) / k;
			++n_actual;
		}
	}
	// m_xy_k is the mean of x * y
	return make_pair(m_xy_k * n_actual, n_actual);
}

double sobol_u_missing_data(const std::vector<double> &x_vec, const std::vector<double> &y_vec, double missing_val)
{
	pair<double, size_t> data = sum_of_prod_missing_data(x_vec, y_vec, missing_val);
	double u = data.first / (data.second - 1.0);
	return u;
}

double si_saltelli_numer(const std::vector<double>& y_a, const std::vector<double>& y_b, const std::vector<double>& y_c, double missing_val)
{
	assert (y_a.size() == y_b.size());
	assert (y_a.size() == y_c.size());
	int valid_count = 0;
	double a, b, c;
	double t = 0;
	for (int i = 0; i < y_a.size(); i++)
	{
		a = y_a[i], b = y_b[i], c = y_c[i];
		if ((a == missing_val) ||
			(b == missing_val) ||
			(c == missing_val))
			continue;
		valid_count++;
		t = t + (b * (c - a));

	}
	t = t / valid_count;

	return t;

}


double sti_saltelli_numer(const std::vector<double>& y_a, const std::vector<double>& y_c, double missing_val)
{
	assert(y_a.size() == y_c.size());
	int valid_count = 0;
	double a, b, c;
	double t = 0;
	for (int i = 0; i < y_a.size(); i++)
	{
		a = y_a[i], c = y_c[i];
		if ((a == missing_val) ||
			(c == missing_val))
			continue;
		valid_count++;
		t = t + pow((a - c), 2);

	}
	t = t / (2.0 * valid_count);
	return t;
}

void RunningStats::reset()
{
	n = 0;
	mk = 0.0;
	mk_abs = 0.0;
	qk = 0.0;
}

void RunningStats::add(double sample)
{
	++n;
	if (n==1)
	{
		mk = sample;
		mk_abs = abs(sample);
		qk = 0.0;
	}
	else
	{
		qk += (n - 1) * pow((sample - mk), 2.0) / n; // must happen before mk is updated
		mk += (sample - mk) / n;
		mk_abs += (abs(sample) - mk_abs) / n;
	}
}

void RunningStats::add(const std::vector<double> &sample)
{
	for (auto &s : sample)
	{
	   add(s);
	}
}

double RunningStats::comp_var() const
{
	return qk / (n-1);
}

double RunningStats::comp_sigma() const
{
	return sqrt(comp_var());
}

double RunningStats::comp_mean() const
{
return mk;
}

double RunningStats::comp_abs_mean() const
{
	return mk_abs;
}

long RunningStats::comp_nsamples() const
{
	return n;
}

void RunningStatsMissingData::add(double sample)
{
	if (sample != missing_value)
	{
		RunningStats::add(sample);
	}
}

void RunningStatsMissingData::add(const std::vector<double> &sample)
{
	for (auto &s : sample)
	{
		add(s);
	}
}

#include <iostream>
#include <list>
#include <iomanip>
#include "GsaAbstractBase.h"
#include "utilities.h"
#include "FileManager.h"
#include "ParamTransformSeq.h"

using namespace std;
using namespace pest_utils;

const double GsaAbstractBase::MISSING_DATA = -9999e50;
mt19937_64 GsaAbstractBase::rand_engine = mt19937_64(1);

GsaAbstractBase::GsaAbstractBase(Pest &_pest_scenario,
		FileManager &_file_manager, ObjectiveFunc *_obj_func_ptr,
		const ParamTransformSeq &_par_transform, PARAM_DIST _par_dist, unsigned int _seed)
	: file_manager_ptr(&_file_manager), obj_func_ptr(_obj_func_ptr), base_partran_seq_ptr(0),
	par_dist(_par_dist), seed(_seed)
{
//	// GSA methods only support one to one transformations
//	if (!_par_transform.is_one_to_one())
//	{
//		throw PestError("Error: GSA methods only support one to one transformations.  Please insure the SVDA transformation is turned off.");
//	}

	//copy transformations and delete all active_ctl to numeric transformation except the log trans
	ofstream &frec = file_manager_ptr->rec_ofstream();

	pest_scenario_ptr = &_pest_scenario;

	gsa_parm_tran_seq.deep_copy(_par_transform);
	gsa_parm_tran_seq.clear_tranSeq_active_ctl2numeric();
	gsa_parm_tran_seq.push_back_active_ctl2numeric(_par_transform.get_log10_ptr()->clone());
	base_partran_seq_ptr = &gsa_parm_tran_seq;
	Parameters inti_pars = _pest_scenario.get_ctl_parameters();
	if (_pest_scenario.get_pestpp_options().get_enforce_tied_bounds())
	{
		pair<Parameters, Parameters> ppair = _pest_scenario.get_effective_ctl_lower_upper_bnd(inti_pars);
		max_numeric_pars = ppair.second;
		min_numeric_pars = ppair.first;


		frec << endl << "  effective parameter bounds (accounting for tied parameters)" << endl;
		frec << setw(50) << left << "parnme" << setw(25) << "initial value" << setw(25) << "control file lower bound" << setw(25) << "effective lower bound" << setw(25) << "control file upper bound" << setw(25) << "effective upper bound" << endl;
		ParameterInfo pi = _pest_scenario.get_ctl_parameter_info();
		for (auto &p : _pest_scenario.get_ctl_ordered_par_names())
			frec << setw(50) << left << setprecision(7) << scientific << p << setw(25) << inti_pars.get_rec(p) << setw(25) << pi.get_parameter_rec_ptr(p)->lbnd << setw(25) << min_numeric_pars.get_rec(p) << setw(25) << pi.get_parameter_rec_ptr(p)->ubnd << setw(25) << max_numeric_pars.get_rec(p) << setw(25) << "" << endl;

	
	}
	else
	{
		for (const auto &i : inti_pars)
		{
			const string &p_name = i.first;
			const ParameterRec *p_info = _pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(p_name);
			max_numeric_pars[p_name] = p_info->ubnd;
			min_numeric_pars[p_name] = p_info->lbnd;
		}
	}

	
	base_partran_seq_ptr->ctl2numeric_ip(max_numeric_pars);
	base_partran_seq_ptr->ctl2numeric_ip(min_numeric_pars);
	for (const auto i : _pest_scenario.get_ctl_ordered_adj_par_names())
	{
		if (max_numeric_pars.find(i) != max_numeric_pars.end())
		{
			adj_par_name_vec.push_back(i);
		}
	}
	obs_name_vec = pest_scenario_ptr->get_ctl_ordered_obs_names();
}


GsaAbstractBase::~GsaAbstractBase(void)
{
}

map<string, string>  GsaAbstractBase::process_gsa_file(ifstream &fin, FileManager &file_manager)
{
	map<string, string> arg_map;
	string line;
	int lnum = 0;
	try {
		while(getline(fin, line))
		{
			++lnum;
			strip_ip(line);
			if (line[0] == '#')
			{
			}
			else
			{
				string line_upper = upper_cp(line);
				parce_line(line_upper, arg_map);
			}
		}
	}
	catch (PestConversionError &e) {
		std::stringstream out;
		out << "Error parsing \"" << file_manager.build_filename("gsa") << "\" on line number " << lnum << endl;
		out << e.what() << endl;
		e.add_front(out.str());
		e.raise();
	}
	return arg_map;
}


map<string, double> GsaAbstractBase::calc_parameter_norm_std_dev()
{
	map<string, double> std_dev_map;
	for (const auto &ipar : adj_par_name_vec)
	{
		double lower = min_numeric_pars[ipar];
		double upper = max_numeric_pars[ipar];
		std_dev_map[ipar] = (upper - lower) / 4.0;
	}
	return std_dev_map;
}

map<string, double> GsaAbstractBase::calc_parameter_unif_std_dev()
{
	map<string, double> std_dev_map;
	for (const auto &ipar : adj_par_name_vec)
	{
		double lower = min_numeric_pars[ipar];
		double upper = max_numeric_pars[ipar];
		std_dev_map[ipar] = (upper - lower) / sqrt(12.0);
	}
	return std_dev_map;
}

void GsaAbstractBase::parce_line(const string &line, map<string, string> &arg_map)
{
	string key;
	string value;
	size_t found = line.find_first_of("#");
	if (found == string::npos) {
		found = line.length();
	}
	string tmp_line = line.substr(0, found);
	strip_ip(tmp_line, "both", "\t\n\r+ ");
	upper_ip(tmp_line);
	list<string> tokens;
	tokenize(tmp_line, tokens, "\t\n\r() ");
	while(tokens.size() >=2)
	{
		string method;
		key = tokens.front();
		tokens.pop_front();
		value = tokens.front();
		tokens.pop_front();
		if (key=="METHOD"){
		}
		else if (key == "RAND_SEED"){
		}
		else if (key=="MORRIS_R"){
		}
		else if (key=="MORRIS_P"){
		}
		else if (key=="SOBOL_SAMPLES"){
		}
		else if (key == "MORRIS_POOLED_OBS"){
		}
		else if (key == "MORRIS_OBS_SEN"){
		}
		else if (key == "MORRIS_DELTA"){
		}
		else if (key == "SOBOL_PAR_DIST"){
		}

		else {
			throw PestParsingError(line, "Invalid key word \"" + key +"\"");
		}
		arg_map[key] = value;
	}
}


vector<double>  GsaAbstractBase::calc_interval_midpoints(int n_interval, double min, double max)
{
	assert(max > min);
	double del_interval = (max - min) / n_interval;
	vector<double> interval_mid_pts;
	interval_mid_pts.reserve(n_interval);
	for (int i=0; i<n_interval; ++i)
	{
		interval_mid_pts.push_back( (.5 + i) * del_interval );
	}
	return interval_mid_pts;
}


double GsaAbstractBase::ltqnorm(double p)
{
	/*
	 * Lower tail quantile for standard normal distribution function.
	 *
	 * This function returns an approximation of the inverse cumulative
	 * standard normal distribution function.  I.e., given P, it returns
	 * an approximation to the X satisfying P = Pr{Z <= X} where Z is a
	 * random variable from the standard normal distribution.
	 *
	 * The algorithm uses a minimax approximation by rational functions
	 * and the result has a relative error whose absolute value is less
	 * than 1.15e-9.
	 *
	 * Author:      Peter John Acklam
	 * Time-stamp:  2002-06-09 18:45:44 +0200
	 * E-mail:      jacklam@math.uio.no
	 * WWW URL:     http://www.math.uio.no/~jacklam
	 *
	 * C implementation adapted from Peter's Perl version
	 */
	/* Coefficients in rational approximations. */
	static const double a[] =
	{
		-3.969683028665376e+01,
		 2.209460984245205e+02,
		-2.759285104469687e+02,
		 1.383577518672690e+02,
		-3.066479806614716e+01,
		 2.506628277459239e+00
	};

	static const double b[] =
	{
		-5.447609879822406e+01,
		 1.615858368580409e+02,
		-1.556989798598866e+02,
		 6.680131188771972e+01,
		-1.328068155288572e+01
	};

	static const double c[] =
	{
		-7.784894002430293e-03,
		-3.223964580411365e-01,
		-2.400758277161838e+00,
		-2.549732539343734e+00,
		 4.374664141464968e+00,
		 2.938163982698783e+00
	};

	static const double d[] =
	{
		7.784695709041462e-03,
		3.224671290700398e-01,
		2.445134137142996e+00,
		3.754408661907416e+00
	};
	const double LOW = 0.02425;
	const double HIGH =  0.97575;

	double q, r;

	errno = 0;

	if (p < 0 || p > 1)
	{
		errno = EDOM;
		return 0.0;
	}
	else if (p == 0)
	{
		errno = ERANGE;
		return -HUGE_VAL /* minus "infinity" */;
	}
	else if (p == 1)
	{
		errno = ERANGE;
		return HUGE_VAL /* "infinity" */;
	}
	else if (p < LOW)
	{
		/* Rational approximation for lower region */
		q = sqrt(-2*log(p));
		return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	}
	else if (p > HIGH)
	{
		/* Rational approximation for upper region */
		q  = sqrt(-2*log(1-p));
		return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	}
	else
	{
		/* Rational approximation for central region */
			q = p - 0.5;
			r = q*q;
		return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
			(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
	}
}

ostream& operator<<(ostream& out, const vector<double> &rhs)
{
	int n=0;
	for(const auto &i : rhs)
	{
		out << "[" << n++ << "] = " << i << endl;
	}
	return out;
}

ostream& operator<<(ostream& out, const pair<double, double> &pr)
{
	out << "(" << pr.first << ", " << pr.second << ")";
	return out;
}

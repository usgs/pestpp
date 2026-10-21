/*


    This file is part of PEST++.

    PEST++ is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PEST++ is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PEST++.  If not, see<http://www.gnu.org/licenses/>.
*/
/**
 * @file PriorInformation.cpp
 * @brief Implementation of PriorInformation.
 */

#include "PriorInformation.h"
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "utilities.h"
#include "Transformable.h"

using namespace::std;
using namespace pest_utils;

/**
 * @brief P i atom.
 *
 * @param _par_name Description.
 * @param _log_trans Description.
 * @param _factor Description.
 *
 * @return Description.
 */
PIAtom::PIAtom(const string &_par_name, bool _log_trans, double _factor)
	: par_name(_par_name), log_trans(_log_trans), factor(_factor)
{

}

PriorInformationRec::PriorInformationRec(double _pival, double _weight,
	const string &_group, const vector<PIAtom> _pi_atoms)
	: pival(_pival), weight(_weight), group(_group), pi_atoms(_pi_atoms)
{
}

/**
 * @brief Prior information rec.
 *
 * @param rhs Description.
 *
 * @return Description.
 */
PriorInformationRec::PriorInformationRec(const PriorInformationRec &rhs)
{
	*this = rhs;
}

const PriorInformationRec& PriorInformationRec::operator=(const PriorInformationRec &rhs)
{
	pi_atoms = rhs.pi_atoms;
	pival = rhs.pival;
	weight = rhs.weight;
	group = rhs.group;
	return *this;
}

/**
 * @brief Get atom factors.
 *
 * @return Description.
 */
map<string, double> PriorInformationRec::get_atom_factors()
{
	map<string, double> piatom_factors;
	for (auto &piatom : pi_atoms)
		piatom_factors[piatom.par_name] = piatom.factor;
	return piatom_factors;
}

/**
 * @brief Calc residual.
 *
 * @param pars Description.
 *
 * @return Description.
 */
double PriorInformationRec::calc_residual(const Parameters &pars) const
{
	double sim_value = 0;
	double par_value;
	Parameters::const_iterator p_iter;
	for (vector<PIAtom>::const_iterator b=pi_atoms.begin(), e=pi_atoms.end();
		b!=e; ++b) {
		par_value = pars.get_rec((*b).par_name);
		if ((*b).log_trans) {
			sim_value += (*b).factor * log10(par_value);
		}
		else {
			sim_value += (*b).factor * par_value;
		}
	}
	return (sim_value - pival);
}


/**
 * @brief Calc sim and resid.
 *
 * @param pars Description.
 *
 * @return Description.
 */
pair<double,double> PriorInformationRec::calc_sim_and_resid(const Parameters &pars) const
{
	double sim_value = 0;
	double par_value;
	Parameters::const_iterator p_iter;
	for (vector<PIAtom>::const_iterator b = pi_atoms.begin(), e = pi_atoms.end();
		b != e; ++b) {
		par_value = pars.get_rec((*b).par_name);
		if ((*b).log_trans) {
			sim_value += (*b).factor * log10(par_value);
		}
		else {
			sim_value += (*b).factor * par_value;
		}
	}
	return pair<double,double>(sim_value,(sim_value - pival));
}

/**
 * @brief Is regularization.
 *
 * @return Description.
 */
bool PriorInformationRec::is_regularization() const
{
	bool is_reg = false;
	int found = upper_cp(group).find("REGUL");
	if (found == 0) is_reg=true;
	return is_reg;
}

/**
 * @brief Destructor for .
 */
PriorInformationRec::~PriorInformationRec(void)
{
}


/**
 * @brief Overloaded operator << operator.
 *
 * @param out Description.
 * @param rhs Description.
 *
 * @return Description.
 */
ostream& operator<< (ostream& out, const PriorInformationRec &rhs)
{
	struct Local {
		static string GetSign(double val) {
			if (val < 0) return "-";
			return "+";
		}
	};

	for(vector<PIAtom>::const_iterator b=rhs.pi_atoms.begin(), e=rhs.pi_atoms.end();
		b!=e; ++b) {
			out << " " << Local::GetSign((*b).factor);
			out << " " << abs((*b).factor);
			out << " * ";
			if ((*b).log_trans)
				out << "LOG(" << (*b).par_name << ")";
			else
				out << (*b).par_name;
			
	}
	out << " = " << rhs.pival;
	out << "   " << rhs.weight;
	out << "   " << rhs.group;
	out << endl;
	return out;
}

/**
 * @brief Add record.
 *
 * @param name Description.
 * @param pi_rec_ptr Description.
 */
void PriorInformation::AddRecord(const string &name, const PriorInformationRec* pi_rec_ptr)
{
	PriorInformationRec pi_rec(pi_rec_ptr->get_obs_value(), pi_rec_ptr->get_weight(), pi_rec_ptr->get_group(), pi_rec_ptr->get_atoms());
	prior_info_map[name] = pi_rec;

}

/**
 * @brief Add record.
 *
 * @param pi_line Description.
 *
 * @return Description.
 */
pair<string, string> PriorInformation::AddRecord(const string& pi_line)
{
	vector<string> tokens;
	tokenize(upper_cp(strip_cp(pi_line)), tokens);
	string par_name;
	double pifac;
	double pival;
	double weight;
	bool log_trans;
	string group;
	bool plus_sign;
	vector<PIAtom> pi_atoms;
	int n_tokens = tokens.size();

	group = tokens[n_tokens - 1];
	convert_ip(tokens[n_tokens - 2], weight);
	convert_ip(tokens[n_tokens - 3], pival);

    if ((weight != 0) && (!isnormal(weight)))
        throw runtime_error("denormal weight '"+tokens[n_tokens - 2]+"' for PI equation "+tokens[0]);
    if ((pival != 0) && (!isnormal(pival)))
        throw runtime_error("denormal pival '"+tokens[n_tokens - 3]+"' for PI equation "+tokens[0]);


    string prior_info_name = tokens[0];
	// process prior information equation
	plus_sign = true;
	for (int i = 1; i < n_tokens - 4; i += 4)
	{
		convert_ip(tokens[i], pifac);
		// token i+1 = "*"
		par_name = tokens[i + 2];
		if (par_name.find("LOG") == 0) // parameter name
		{
			log_trans = true;
			par_name = par_name.substr(4, par_name.size() - 5); // remove "log(" from front and ")" from rear of tokens[i+1]
		}
		else {
			log_trans = false;
		}
		if (!plus_sign)
			pifac *= -1.0;

		// add term to vector storing prior information equation
		pi_atoms.push_back(PIAtom(par_name, log_trans, pifac));


		if (tokens[i + 3] == "+") {
			plus_sign = true;
		}
		else {
			plus_sign = false;
		}
	}
	prior_info_map[prior_info_name] = PriorInformationRec(pival, weight, group, pi_atoms);
	return pair<string, string>(prior_info_name, group);
}


/**
 * @brief Add record.
 *
 * @param tokens Description.
 *
 * @return Description.
 */
pair<string, string> PriorInformation::AddRecord(const vector<string> tokens)
{
	//here we expect tokens to have the entire pi equation as a string, so we need to split it up and insert it into tokens
	vector<string> eq_tokens;// , full_tokens;
	tokenize(tokens[1], eq_tokens);
	//full_tokens = tokens;
	//vector<string>::iterator it = full_tokens.begin() + 1;// distance(full_tokens.begin(), 1);
	//full_tokens.insert(it, eq_tokens.begin(), eq_tokens.end());
	//it = full_tokens.end() - 3;
	//full_tokens.erase(it);
	string par_name;
	double pifac;
	double pival;
	double weight;
	bool log_trans;
	string group;
	bool plus_sign;
	vector<PIAtom> pi_atoms;
	int n_tokens = tokens.size();
	int n_tokens_eq = eq_tokens.size();

	group = tokens[n_tokens-1];
	convert_ip(tokens[n_tokens-2], weight);
	convert_ip(eq_tokens[n_tokens_eq-1], pival);

	string prior_info_name = tokens[0];
	// process prior information equation
	plus_sign = true;
	for (int i=0; i< n_tokens_eq-4; i+=4)
	{
		convert_ip(eq_tokens[i], pifac);
		// token i+1 = "*"
		par_name = eq_tokens[i+2];
		if (par_name.find("LOG") == 0) // parameter name
		{
			log_trans = true;
			par_name = par_name.substr(4, par_name.size() - 5); // remove "log(" from front and ")" from rear of tokens[i+1]
		}
		else {
			log_trans = false;
		}
		if (!plus_sign)
			pifac *= -1.0;

		// add term to vector storing prior information equation
		pi_atoms.push_back(PIAtom(par_name, log_trans, pifac));


		if (eq_tokens[i+3] == "+") {
			plus_sign = true;
		}
		else {
			plus_sign = false;
		}
	}
	prior_info_map[prior_info_name] = PriorInformationRec(pival, weight, group, pi_atoms);
	return pair<string, string>(prior_info_name, group);
}

/**
 * @brief Get keys.
 *
 * @return Description.
 */
vector<string> PriorInformation::get_keys() const
{
	vector<string> ret_val;

	for(std::map<std::string, PriorInformationRec>::const_iterator b=prior_info_map.begin(), e=prior_info_map.end();
		b!=e; ++b)
	{
	ret_val.push_back((*b).first);
	}
	return ret_val;
}

/**
 * @brief Get nnz pi.
 *
 * @return Description.
 */
int PriorInformation::get_nnz_pi() const
{
	int nnz = 0;
	for (const auto& pi : prior_info_map)
	{
		if (pi.second.get_weight() > 0.0)
			nnz++;
	}
	return nnz;
}

PriorInformation::IrlsStats PriorInformation::irls_reweight(const Parameters &pars, double eps, map<string, double> &w0, int max_dev)
{
	//the l1 penalty by reweighting.  pest applies weights as (w*r)^2 in phi, so to make the
	//penalty behave like |r| we need w^2 * r^2 = w0^2 * |r|, which means w = w0/sqrt(|r|).
	//(w = w0/|r| would give a constant penalty w0^2, which isnt l1 at all).  the floor eps
	//keeps a residual that lands on zero from making the weight infinite and welding that
	//equation shut for the rest of the run - it can still move if the data want it to.  w0
	//is the control file weight, stored once on the first pass so the reweighting doesnt
	//compound from iteration to iteration.  only regularization group equations are touched
	//and zero-weight ones stay zero.  the residual is whatever calc_residual gives, so for
	//log transformed pars it is in log10 units and eps has to be too.
	//
	//max_dev >= 0 turns eps into a lower bound and solves for the floor instead: the floor is
	//raised to the (max_dev+1)-th largest |residual|, so at most max_dev equations sit above it
	//and get the l1 treatment while the rest share one quadratic weight.  that is the whole
	//"solve for eps" - with the residuals fixed at the accepted parameters the floor that leaves
	//n equations active is just an order statistic, no search needed.  eps passed in still
	//applies as the minimum so the floor can never chase shrinking residuals down to zero.
	IrlsStats stats;
	vector<double> facs;
	if (w0.empty())
	{
		for (auto &pi : prior_info_map)
			w0[pi.first] = pi.second.get_weight();
	}
	struct Eq { PriorInformationRec* rec; double r; double w0; };
	vector<Eq> eligible;
	for (auto &pi : prior_info_map)
	{
		if (!pi.second.is_regularization())
			continue;
		map<string, double>::const_iterator it = w0.find(pi.first);
		if ((it == w0.end()) || (it->second <= 0.0))
			continue;
		double r = abs(pi.second.calc_residual(pars));
		if (!std::isfinite(r))
			continue;
		eligible.push_back(Eq{&pi.second, r, it->second});
	}
	if ((max_dev >= 0) && (max_dev < (int)eligible.size()))
	{
		vector<double> rs;
		for (auto &e : eligible)
			rs.push_back(e.r);
		sort(rs.begin(), rs.end(), greater<double>());
		//rs[max_dev] is the (max_dev+1)-th largest.  "active" below means strictly above the
		//floor, so this equation and everything tied with it land at the floor
		eps = max(eps, rs[max_dev]);
	}
	stats.eps = eps;
	for (auto &e : eligible)
	{
		if (e.r < eps)
			stats.n_floor++;
		if (e.r > eps)
			stats.n_active++;
		double fac = 1.0 / sqrt(max(e.r, eps));
		e.rec->set_weight(e.w0 * fac);
		facs.push_back(fac);
	}
	stats.n = facs.size();
	if (stats.n > 0)
	{
		sort(facs.begin(), facs.end());
		stats.fmin = facs.front();
		stats.fmax = facs.back();
		stats.fmed = facs[facs.size() / 2];
	}
	return stats;
}

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
#ifndef PRIORINFORMATION_H_
#define PRIORINFORMATION_H_


#include <vector>
#include <string>
#include <iostream>
#include <map>

class Parameters;

class PIAtom
{
public:
	PIAtom(const std::string &par_name, bool log_trans, double factor);
	std::string par_name;
	bool log_trans;
	double factor;
};

class PriorInformationRec
{
public:
	friend std::ostream& operator<< (std::ostream& out, const PriorInformationRec &rhs);
	PriorInformationRec(const PriorInformationRec &rhs);
	PriorInformationRec(double _pival=0.0, double _weight=0.0, const std::string &_group="",
		const std::vector<PIAtom> _pi_atoms = std::vector<PIAtom>());
	const PriorInformationRec& operator=(const PriorInformationRec &rhs);
	double calc_residual(const Parameters &pars) const;
	std::pair<double, double> calc_sim_and_resid(const Parameters &pars) const;
	bool is_regularization() const;
	double get_weight()const {return weight;}
	void set_weight(double _weight) { weight = _weight; }
	double get_obs_value()const {return pival;}
	const std::string& get_group() const{return group;}
	const std::string *get_group_ptr() const{return &group;}
	std::map<std::string, double> get_atom_factors();
	std::vector<PIAtom> get_atoms()const { return pi_atoms; }
	~PriorInformationRec(void);
private:
	std::vector<PIAtom> pi_atoms;
	double pival;
	double weight;
	std::string group;
};


class PriorInformation
{
public:
	typedef std::map<std::string, PriorInformationRec>::iterator iterator;
	typedef std::map<std::string, PriorInformationRec>::const_iterator const_iterator;
	PriorInformation() {}
	~PriorInformation() {}
	std::pair<std::string, std::string> AddRecord(const std::string &pi_line);
	std::pair<std::string, std::string> AddRecord(const std::vector<std::string> tokens);
	void AddRecord(const std::string &name, const PriorInformationRec* pi_rec_ptr);
	PriorInformation::iterator begin(){return prior_info_map.begin();}
	PriorInformation::const_iterator begin() const {return prior_info_map.begin();}
	PriorInformation::iterator end() {return prior_info_map.end();}
	PriorInformation::const_iterator end() const {return prior_info_map.end();}
	PriorInformation::const_iterator find(const std::string &key) const {return  prior_info_map.find(key);}
	PriorInformation::iterator find(const std::string &key) {return  prior_info_map.find(key);}
	PriorInformationRec get_pi_rec(std::string name)const { return prior_info_map.at(name); }
	size_t size() const {return prior_info_map.size();}
	int get_nnz_pi() const;
	std::vector<std::string> get_keys() const;
	void clear(){prior_info_map.clear();}
	//what one irls reweighting pass did, for the rec file
	struct IrlsStats
	{
		int n = 0;
		int n_floor = 0;
		int n_active = 0;   //equations with |residual| above the floor - the deviations let through
		double eps = 0.0;   //the floor actually used, raised above the one passed in when max_dev asks
		double fmin = 0.0;
		double fmed = 0.0;
		double fmax = 0.0;
	};
	//max_dev >= 0 raises the floor to the value that leaves at most max_dev equations above it
	IrlsStats irls_reweight(const Parameters &pars, double eps, std::map<std::string, double> &w0, int max_dev = -1);
private:
	std::map<std::string, PriorInformationRec> prior_info_map;
};


std::ostream& operator<< (std::ostream& out, const PriorInformationRec &rhs);

#endif //PRIORINFORMATION_H_

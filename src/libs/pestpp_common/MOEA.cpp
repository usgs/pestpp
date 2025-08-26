#include <random>
#include <iomanip>
#include <iterator>
#include <map>
#include "MOEA.h"
#include "Ensemble.h"
#include "RunManagerAbstract.h"
#include "ModelRunPP.h"
#include "RestartController.h"
#include "EnsembleMethodUtils.h"
#include "constraints.h"
#include "eigen_tools.h"
#include "RedSVD-h.h"


using namespace std;

//util functions
typedef std::function<bool(std::pair<std::string, double>, std::pair<std::string, double>)> Comparator;
// Defining a lambda function to compare two pairs. It will compare two pairs using second field
Comparator compFunctor = [](std::pair<std::string, double> elem1, std::pair<std::string, double> elem2)
{
	return elem1.second < elem2.second;
};

typedef std::set<std::pair<std::string, double>, Comparator> sortedset;



ParetoObjectives::ParetoObjectives(Pest& _pest_scenario, FileManager& _file_manager,
	PerformanceLog* _performance_log):
	pest_scenario(_pest_scenario),file_manager(_file_manager),
	performance_log(_performance_log)

{

}

map<string, map<string, double>> ParetoObjectives::get_member_struct(ObservationEnsemble& op, ParameterEnsemble& dp)
{
	map<string, map<string, double>> _member_struct;

	//map<string, map<double, string>> obj_struct;
	vector<string> real_names = op.get_real_names();
	Eigen::VectorXd obj_vals, obj_sd_vals;
	map<string, map<string, double>> temp;
	for (auto obj_name : *obs_obj_names_ptr)
	{
		obj_vals = op.get_eigen(vector<string>(), vector<string>{obj_name});

		//if this is a max obj, just flip the values here
		if (obj_dir_mult_ptr->find(obj_name) != obj_dir_mult_ptr->end())
		{
			obj_vals *= obj_dir_mult_ptr->at(obj_name);
		}
		map<double, string> obj_map;
		map<string, double> t;
		for (int i = 0; i < real_names.size(); i++)
		{
			obj_map[obj_vals[i]] = real_names[i];
			t[real_names[i]] = obj_vals[i];
		}
		//obj_struct[obj_name] = obj_map;
		temp[obj_name] = t;

	}

	//add variance info for obj values to member_struct
	if (prob_pareto) {
		map<double, string> obj_sd_map;
		map<string, double> t;
		
		for (auto obj_sd_name : *obs_obj_sd_names_ptr)
		{
			obj_sd_vals = op.get_eigen(vector<string>(), vector<string>{obj_sd_name});
					

			for (int i = 0; i < real_names.size(); i++)
			{
				//obj_sd_map[obj_sd_vals[i]] = real_names[i];
				t[real_names[i]] = obj_sd_vals[i];
			}
			temp[obj_sd_name] = t;

		}
		for (auto obj_sd_name : *pi_obj_sd_names_ptr)
		{
			obj_sd_vals = op.get_eigen(vector<string>(), vector<string>{obj_sd_name});

			for (int i = 0; i < real_names.size(); i++)
			{
				//obj_sd_map[obj_sd_vals[i]] = real_names[i];
				t[real_names[i]] = obj_sd_vals[i];
			}
			temp[obj_sd_name] = t;

		}
	}
	

	//map<string, map<string, double>> member_struct;

	for (auto real_name : real_names)
	{
		map<string, double> obj_map;
		for (auto t : temp)
		{
			//obj_map[t.second[real_name]] = t.first;
			obj_map[t.first] = t.second[real_name];
		}

		_member_struct[real_name] = obj_map;
	}
	temp.clear();

	



	//add any prior info obj values to member_struct
	if (pi_obj_names_ptr->size() > 0)
	{
		PriorInformation* pi_ptr = pest_scenario.get_prior_info_ptr();
		
		pair<double, double> pi_res_sim;
		ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();
		vector<string> pnames = dp.get_var_names();
		for (auto real_name : real_names)
		{
			Parameters pars = pest_scenario.get_ctl_parameters(); 
			pars.update_without_clear(pnames, dp.get_real_vector(real_name));
			if (dp.get_trans_status() == ParameterEnsemble::transStatus::NUM)
				pts.numeric2ctl_ip(pars);
			for (auto obj_name : *pi_obj_names_ptr)
			{
				pi_res_sim = pi_ptr->get_pi_rec(obj_name).calc_sim_and_resid(pars);
				//account for dir mult here
				_member_struct[real_name][obj_name] = pi_res_sim.first * obj_dir_mult_ptr->at(obj_name);
			}
		}
	}
	return _member_struct;

}

bool ParetoObjectives::compare_two(string& first, string& second, MouEnvType envtyp)
{
    if (infeas.find(first) != infeas.end())
	{	//if both are infeas, select the solution that is less infeasible
		if (infeas.find(second) != infeas.end())
		{
			if (infeas[first] > infeas[second])
				return false;
			else
				return true;
		}
		return false;
	}
	if (infeas.find(second) != infeas.end())
		return true;
	if (envtyp == MouEnvType::NSGA)
		return compare_two_nsga(first, second);
	else if (envtyp == MouEnvType::SPEA)
		return compare_two_spea(first, second);
	else
		throw runtime_error("ParetoObjectives::compare_two(): unrecognized envtyp");
}


bool ParetoObjectives::compare_two_spea(string& first, string& second)
{
	if (spea2_constrained_fitness_map.at(first) < spea2_constrained_fitness_map.at(second))
		return true;
	else
		return false;
}

bool ParetoObjectives::compare_two_nsga(string& first, string& second)
{
	if (member_front_map.at(first) < member_front_map.at(second))
		return true;
	if (member_front_map.at(first) == member_front_map.at(second))
	{
		if (crowd_map.at(first) > crowd_map.at(second))
			return true;
		else
			return false;
	}
	return false;
}

void ParetoObjectives::drop_duplicates(map<string, map<string, double>>& _member_struct)
{
	performance_log->log_event("checking for duplicate solutions");
	duplicates.clear();
	

	//performance_log->log_event("checking for duplicate solutions");
	vector<string> names;
	map<string, double> solution_p,solution_q;
	string name_p,name_q;
	for (auto m : _member_struct)
		names.push_back(m.first);

	for (int i=0;i<names.size();i++)
	{
		name_p = names[i];
		solution_p = _member_struct[names[i]];
		//for (auto solution_q : member_struct)
		for (int j=i+1;j<names.size();j++)
		{
			name_q = names[j];
			solution_q = _member_struct[name_q];
			//if (solution_p.first == solution_q.first)
			//	continue;
			if ((first_equals_second(solution_p, solution_q)) && (duplicates.find(name_q) == duplicates.end()))
			{
				duplicates.emplace(name_p);
			}
			
		}
	}
	if (duplicates.size() > 0)
	{
		stringstream ss;
		ss << "WARNING: " << duplicates.size() << " duplicate solutions found, dropping (see rec file for listing)";
		performance_log->log_event(ss.str());
		cout << ss.str() << endl;
		ofstream& frec = file_manager.rec_ofstream();
		frec << "WARNING: " << duplicates.size() << " duplicate solutions found:" << endl;
		int i = 0;
		for (auto d : duplicates)
		{
			frec << d << " ";
			i++;
			if (i > 5)
			{
				frec << endl;
				i = 0;
			}
		}
		frec << endl;
		for (auto d : duplicates)
			member_struct.erase(d);	
	}
}

map<string, double> ParetoObjectives::get_mopso_fitness(vector<string> members, ObservationEnsemble& op, ParameterEnsemble& dp)
{
	map<string, map<string, double>> _member_struct = get_member_struct(op, dp);
	return get_mopso_fitness(members, _member_struct);
}

map<string, double> ParetoObjectives::get_mopso_fitness(vector<string> members, map<string, map<string, double>>& _member_struct)
{
	double alpha = pest_scenario.get_pestpp_options().get_mou_pso_alpha();
	stringstream ss;
	if (alpha == 0)
	{
		//double maxarchivesize = pest_scenario.get_pestpp_options().get_mou_max_archive_size();
		double maxarchivesize = pest_scenario.get_pestpp_options().get_mou_max_archive_size();
		double pfull = members.size() / maxarchivesize;
		double rramp = pest_scenario.get_pestpp_options().get_mou_pso_rramp();
		double rfit = pest_scenario.get_pestpp_options().get_mou_pso_rfit();

		if (rramp == -5e+02)
			throw runtime_error("PSO alpha is zero");
		if (rramp == 0.0)
			throw runtime_error("PSO RRAMP is zero");

		alpha = 1 + (exp(rramp * pfull) - 1.0) / (exp(rramp) - 1) * (rfit - 1.0);

		stringstream ss;
		ss.str("");
		ss << "Computing fitness using alpha = " << alpha;
		performance_log->log_event(ss.str());
	}

	map<string, double> fitness;

	if (prob_pareto)
	{
		map<string, double> cluster_crowding = get_cluster_crowding_fitness(members);


		//normalize cd
		double mx = -1.0e+30;
		double mn = 1.0e+30;
		for (auto& cd : cluster_crowding)
		{
			if ((cd.second > mx) && (cd.second != CROWDING_EXTREME))
				mx = cd.second;
			else if (members.size() == 2)
				mx = 0.0;

			if ((cd.second < mn) && (cd.second != CROWDING_EXTREME))
				mn = cd.second;
			else if (members.size() == 2)
				mn = 0.0;
		}
		if (mx < 0.0)
		{
	        ss.str("");
	        ss << "WARNING: pso gbest solution max crowding distance == 0.0" << endl;
	        file_manager.rec_ofstream() << ss.str();
	        cout << ss.str();
	        mx = 0.0;
    	}

		for (auto& cd : cluster_crowding) {
			if (cd.second == CROWDING_EXTREME)
			{
				cd.second = 1;
			}
			else if (mx >= 0.0) {
				cd.second = pow(1 - (cd.second - mn) / (mx - mn + 1), alpha);
			}
			else {
				cd.second = pow(0.5, alpha);
			}
		}

		fitness = cluster_crowding;
	}
	else
	{
		stringstream ss;
		map<string, double> crowd_dist = get_cuboid_crowding_distance(members);
		sortedset crowd_sorted(crowd_dist.begin(), crowd_dist.end(), compFunctor);
		//normalize cd
		double mx = -1.0e+30;
		for (auto& cd : crowd_dist)
			if ((cd.second != CROWDING_EXTREME) && (cd.second > mx))
				mx = cd.second;
			else if (members.size() == 2)
				mx = cd.second;
			else if (crowd_sorted.size() == 1)
				mx = cd.second;
		if (mx < 0.0)
		{
	        ss.str("");
	        ss << "WARNING: pso gbest solution max crowding distance == 0.0" << endl;
	        file_manager.rec_ofstream() << ss.str();
	        cout << ss.str();
	        mx = 0.0;
    	}


		for (auto& cd : crowd_dist) {
			if (cd.second == CROWDING_EXTREME) {
				cd.second = 1.0;
			}
			else if (mx != 0.0) {
				cd.second = pow(cd.second / mx, alpha);
			}
			else {
				cd.second = pow(0.5, alpha);
			}
		}

		fitness = crowd_dist;

	}

	return fitness;
}

void ParetoObjectives::get_spea2_archive_names_to_keep(int num_members, vector<string>& keep, const ObservationEnsemble& op, const ParameterEnsemble& dp)
{
	ParameterEnsemble temp_dp = dp;
	ObservationEnsemble temp_op = op;
	string rm;
	map<string, double> kdist = get_spea2_kth_nn_crowding_distance(temp_op, temp_dp);
	while (keep.size() > num_members)
	{
		temp_dp.keep_rows(keep);
		temp_op.keep_rows(keep);
		kdist = get_spea2_kth_nn_crowding_distance(temp_op, temp_dp);
		sortedset fit_sorted(kdist.begin(), kdist.end(), compFunctor);
		rm = fit_sorted.begin()->first;
		keep.erase(remove(keep.begin(), keep.end(), rm), keep.end());
		//cout << keep.size() << endl;
	}
}


void ParetoObjectives::update(ObservationEnsemble& op, ParameterEnsemble& dp, Constraints* constraints_ptr)
{
	stringstream ss;
	ss << "ParetoObjectives::update() for  " << op.shape().first << " population members";
	performance_log->log_event(ss.str());
	performance_log->log_event("preparing fast-lookup containers");
	ofstream& frec = file_manager.rec_ofstream();

	//update the member struct container
	member_struct = get_member_struct(op, dp);

	drop_duplicates(member_struct);

	if (member_struct.size() == 0)
		throw runtime_error("ParetoObjectives error: member_struct is empty");

	vector<string> real_names;
	for (auto& m : member_struct)
		real_names.push_back(m.first);
	infeas_ordered.clear();
	infeas.clear();
	bool all_infeas = false;
	//TODO: work out if there are actually any constraints (instead of all objectives)
	bool check_constraints = false;
	if (constraints_ptr)
	{
		
		vector<string> t = *obs_obj_names_ptr;
		set<string> obj_names(t.begin(), t.end());
		for (auto name : constraints_ptr->get_obs_constraint_names())
		{
			if (obj_names.find(name) == obj_names.end())
			{
				check_constraints = true;
				break;
			}
		}

		if (!check_constraints)
		{
			obj_names.clear();
			t = *pi_obj_names_ptr;
			obj_names.insert(t.begin(), t.end());
			for (auto name : constraints_ptr->get_pi_constraint_names())
			{
				if (obj_names.find(name) == obj_names.end())
				{
					check_constraints = true;
					break;
				}
			}

		}
		if (check_constraints) {
            performance_log->log_event("feasible sorting");
            //sort into feasible and infeasible
            map<string, double> violations;
            double vsum;
            Observations obs = pest_scenario.get_ctl_observations();
            Parameters pars = pest_scenario.get_ctl_parameters();
            vector<string> onames = op.get_var_names(), pnames = dp.get_var_names();
            set<string> obs_obj_set(obs_obj_names_ptr->begin(), obs_obj_names_ptr->end());
			set<string> obs_obj_sd_set(obs_obj_sd_names_ptr->begin(), obs_obj_sd_names_ptr->end());
            set<string> pi_obj_set(pi_obj_names_ptr->begin(), pi_obj_names_ptr->end());
            set<string>::iterator end;
            ObservationInfo *oi = pest_scenario.get_observation_info_ptr();
            PriorInformation *pi = pest_scenario.get_prior_info_ptr();
            feas_member_struct.clear();

            for (auto real_name : real_names) {
                vsum = 0.0;
                obs.update_without_clear(onames, op.get_real_vector(real_name));
                //TODO: add a constraint tol ++ arg and use it here
                // the 'false' arg is to not apply risk shifting to the satisfaction calcs since
                // 'op' has already been shifted
                violations = constraints_ptr->get_unsatified_obs_constraints(obs, 0.0, false);
                end = obs_obj_set.end();
                for (auto v : violations) {
                    if (obs_obj_set.find(v.first) == end)
                        vsum += pow(v.second * oi->get_weight(v.first), 2);
                }
                pars.update_without_clear(pnames, dp.get_real_vector(real_name));
                violations = constraints_ptr->get_unsatified_pi_constraints(pars, 0.0);
                end = pi_obj_set.end();
                for (auto v : violations) {
                    if (pi_obj_set.find(v.first) == end)
                        vsum += pow(v.second * pi->get_pi_rec(v.first).get_weight(), 2);
                }
                if (vsum > 0.0) {
                    infeas[real_name] = vsum;
                } else
                    feas_member_struct[real_name] = member_struct[real_name];
            }
            /*if (feas_member_struct.size() == 0)
            {
                ss.str("");
                ss << "WARNING: all members are infeasible" << endl;
                frec << ss.str();
                cout << ss.str();
                all_infeas = true;
            }
            else
            {
                ss.str("");
                ss << feas_member_struct.size() << " feasible solutions" << endl;
                frec << ss.str();
                cout << ss.str();
                member_struct = feas_member_struct;
            }*/

            //sort the infeasible members by violation
            vector<pair<string, double>> infeas_vec;
            for (auto inf : infeas)
                infeas_vec.push_back(inf);

            std::sort(infeas_vec.begin(), infeas_vec.end(),
                      compFunctor);

            for (auto inf : infeas_vec)
                infeas_ordered.push_back(inf.first);
        }
		else
        {
            feas_member_struct = member_struct;
        }
	}


	performance_log->log_event("calculating fitness");
	pair<map<string, double>, map<string, double>> fitness_maps = get_spea2_fitness(member_struct);
	spea2_constrained_fitness_map = fitness_maps.first;
	spea2_unconstrained_fitness_map = fitness_maps.second;

	//remove infeasible solutions - after spea fitness calcs
	if (check_constraints)
	{
		if (infeas.size() == member_struct.size())
		{
			ss.str("");
			ss << "WARNING: all members are infeasible" << endl;
			frec << ss.str();
			cout << ss.str();
			all_infeas = true;
		}
		else
		{
			ss.str("");
			ss << feas_member_struct.size() << " feasible solutions" << endl;
			frec << ss.str();
			cout << ss.str();
			member_struct = feas_member_struct;
		}
	}
	performance_log->log_event("pareto front sorting");
	
	if (all_infeas)
	{
		front_map.clear();
		for (int i = 0; i < infeas_ordered.size(); i++)
			front_map[i].push_back(infeas_ordered[i]);	
	}
		
	else
		front_map = sort_members_by_dominance_into_fronts(member_struct);

	return;
}



map<string, double> ParetoObjectives::get_spea2_fitness(int generation, ObservationEnsemble& op, ParameterEnsemble& dp, Constraints* constraints_ptr, bool report, string sum_tag)
{
	stringstream ss;
	ofstream& frec = file_manager.rec_ofstream();
	ss << "ParetoObjectives::get_spea2_fitness() for " << op.shape().first << " population members";
	performance_log->log_event(ss.str());
	update(op, dp, constraints_ptr);

	if (member_struct.size() == 0)
		throw runtime_error("ParetoObjectives error: member_struct is empty");

	if (sum_tag.size() > 0)
	{
		vector<string> real_names = op.get_real_names();
		write_pareto_summary(sum_tag, generation, op, dp,constraints_ptr);
	}
	return spea2_constrained_fitness_map;
}

pair<vector<string>, vector<string>> ParetoObjectives::get_nsga2_pareto_dominance(int generation, ObservationEnsemble& op,
	ParameterEnsemble& dp, Constraints* constraints_ptr, bool sort_ppd, bool report, string sum_tag)
{
	ppd_sort = sort_ppd;
	iter = generation;
	//prob_pareto = sort_ppd;
	stringstream ss;
	ofstream& frec = file_manager.rec_ofstream();
	ss << "ParetoObjectives::get_nsga2_pareto_dominance() for " << op.shape().first << " population members";
	performance_log->log_event(ss.str());
	update(op, dp, constraints_ptr);

	if (member_struct.size() == 0)
		throw runtime_error("ParetoObjectives::get_nsga2_pareto_dominance() error: member_struct is empty");
	

	if (obs_obj_names_ptr->size() + pi_obj_names_ptr->size() > 1)
	{
		frec << "...pareto dominance sort yielded " << front_map.size() << " domination fronts" << endl;

		for (auto front : front_map)
		{
			if (front.second.size() == 0)
			{
				ss.str("");
				ss << "ParetoObjectives::get_nsga2_pareto_dominance() error: front " << front.first << " has no members";
				performance_log->log_event(ss.str());
				throw runtime_error(ss.str());
			}
			if (front.first < 10)
			    frec << front.second.size() << " in the front " << front.first << endl;
			else if (front.first == 10)
			    frec << "(not reporting remaining front membership)" << endl;
		}
	}

	vector<string> nondom_crowd_ordered,dom_crowd_ordered;
	vector<string> crowd_ordered_front;
	crowd_map.clear();
	member_front_map.clear();
	for (auto front : front_map)
	{	
		for (auto m : front.second)
			member_front_map[m] = front.first;
		//TODO: Deb says we only need to worry about crowding sort if not all
		//members of the front are going to be retained.  For now, just sorting all fronts...
		
		if (front.second.size() == 1)
		{
			crowd_ordered_front = front.second;
			crowd_map[front.second[0]] = -999.0;
			fitness_map[front.second[0]] = -999;
			if (sort_ppd)
			{
				for (auto& obj_name : *obs_obj_names_ptr)
				{
					if (member_struct[front.second[0]].at(ppd_obj_to_sd_ptr->at(obj_name)) < min_sd.at(obj_name) - FLOAT_EPSILON)
						member_struct[front.second[0]][ppd_obj_to_sd_ptr->at(obj_name) + "_SYN"] = min_sd.at(obj_name);
					else
						member_struct[front.second[0]][ppd_obj_to_sd_ptr->at(obj_name) + "_SYN"] = member_struct[front.second[0]].at(ppd_obj_to_sd_ptr->at(obj_name));
					//member_struct[front.second[0]][ppd_obj_to_sd_ptr->at(obj_name) + "_SYN"] = min_sd.at(obj_name);
				}
			}
		}

		else
		{
			crowd_ordered_front = sort_members_by_crowding_distance(front.first, front.second, crowd_map, member_struct);
		}

		if (front.first == 1)
			for (auto front_member : crowd_ordered_front)
				nondom_crowd_ordered.push_back(front_member);
		else
			for (auto front_member : crowd_ordered_front)
				dom_crowd_ordered.push_back(front_member);
	}

	//now add the infeasible members
	//if there is at least one feasible nondom solution, then add the infeasible ones to dom solutions
	bool all_infeas = true;
	if (infeas.size() < op.shape().first - duplicates.size())
	{

		if (nondom_crowd_ordered.size() > 0)
			for (auto inf : infeas_ordered)
				dom_crowd_ordered.push_back(inf);
		else
			for (auto inf : infeas_ordered)
				nondom_crowd_ordered.push_back(inf);
		all_infeas = false;
	}

	if (all_infeas)
	{
		if (infeas.size() != nondom_crowd_ordered.size() + dom_crowd_ordered.size())
		{
			ss.str("");
			ss << "ParetoObjectives::get_nsga2_pareto_dominance() internal error: final sorted population size: " <<
				nondom_crowd_ordered.size() + dom_crowd_ordered.size() << " != member_struct size: " << member_struct.size();
			cout << ss.str();
			throw runtime_error(ss.str());
		}
	}
	
	else if (infeas.size() + member_struct.size() != nondom_crowd_ordered.size() + dom_crowd_ordered.size())
	{
		ss.str("");
		ss << "ParetoObjectives::get_nsga2_pareto_dominance() internal error: final sorted population size: " <<
			nondom_crowd_ordered.size() + dom_crowd_ordered.size() << " != member_struct size: " << member_struct.size();
		cout << ss.str();
		throw runtime_error(ss.str());
	}
	

	if (sum_tag.size() > 0)
	{
		vector<string> real_names = op.get_real_names();
		write_pareto_summary(sum_tag, generation, op, dp, constraints_ptr);
	}	
	return pair<vector<string>, vector<string>>(nondom_crowd_ordered, dom_crowd_ordered);
}

void ParetoObjectives::write_pareto_summary(string& sum_tag, int generation, ObservationEnsemble& op, ParameterEnsemble& dp, Constraints* constr_ptr)
{
	//update(op, dp, constr_ptr);
	if (!file_manager.check_ofile_tag_exists(sum_tag))
		file_manager.open_ofile_ext(sum_tag);
	ofstream& sum = file_manager.get_ofstream(sum_tag);
	for (auto member : op.get_real_names())
	{
		if (member_struct.find(member) == member_struct.end())
			continue;
		sum << generation << "," << pest_utils::lower_cp(member);
		for (auto obj : *obs_obj_names_ptr)
		{
			sum << "," << obj_dir_mult_ptr->at(obj) * member_struct[member][obj];
		}
		for (auto obj : *pi_obj_names_ptr)
		{
			sum << "," << obj_dir_mult_ptr->at(obj) * member_struct[member][obj];
		}
		if (prob_pareto)
		{
			for (auto objsd : *obs_obj_sd_names_ptr)
			{
				sum << "," << member_struct[member].at(objsd);
			}
			for (auto objsd : *obs_obj_sd_names_ptr)
			{
				sum << "," << member_struct[member].at(objsd + "_SYN");
			}
		}
		sum << "," << member_front_map[member];
		sum << "," << crowd_map[member];
		if (prob_pareto)
			sum << "," << nn_map[member];
		
		sum << "," << spea2_constrained_fitness_map[member];
		sum << "," << spea2_unconstrained_fitness_map[member];
		if (infeas.find(member) != infeas.end())
		{
			sum << "," << 0 << "," << infeas[member];
		}
		else
		{
			sum << "," << 1 << "," << -999;
		}
		sum << endl;
	}
	cout << "...wrote pareto summary to " << file_manager.get_base_filename() << "." << sum_tag << endl;
	file_manager.rec_ofstream() << "...wrote pareto summary to " << file_manager.get_base_filename() << "." << sum_tag << endl;

}

void ParetoObjectives::prep_pareto_summary_file(string summary_tag)
{
	file_manager.open_ofile_ext(summary_tag);
	ofstream& sum = file_manager.get_ofstream(summary_tag);
	sum << "generation,member";
	for (auto obj : *obs_obj_names_ptr)
		sum << "," << pest_utils::lower_cp(obj);
	for (auto obj : *pi_obj_names_ptr)
		sum << "," << pest_utils::lower_cp(obj);
	if (prob_pareto)
	{
		for (auto objsd : *obs_obj_sd_names_ptr)
			sum << "," << pest_utils::lower_cp(objsd);
		for (auto objsd : *obs_obj_sd_names_ptr)
			sum << "," << pest_utils::lower_cp(objsd + "_SYN");
		sum << ",nsga2_front,nsga2_crowding_distance,nn_count,spea2_unconstrained_fitness,spea2_constrained_fitness,is_feasible,feasible_distance" << endl;
	}
	else
		sum << ",nsga2_front,nsga2_crowding_distance,spea2_unconstrained_fitness,spea2_constrained_fitness,is_feasible,feasible_distance" << endl;
}

map<string, double> ParetoObjectives::get_spea2_kth_nn_crowding_distance(ObservationEnsemble& op, ParameterEnsemble& dp)
{
	//this updates the complicated map-based structure that stores the member names: obj_names:value nested pairs
	map<string, map<string, double>> _member_struct = get_member_struct(op, dp);
	//just reuse the same routine used for the pareto dominance sort...
	return get_spea2_kth_nn_crowding_distance(_member_struct);
}


map<string, double> ParetoObjectives::get_spea2_kth_nn_crowding_distance(map<string, map<string, double>>& _member_struct)
{
	vector<string> members;
	for (auto m : _member_struct)
		members.push_back(m.first);
	return get_spea2_kth_nn_crowding_distance(members, _member_struct);
}

map<string, double> ParetoObjectives::get_spea2_kth_nn_crowding_distance(vector<string>& members, map<string, map<string, double>>& _member_struct)
{
	map<string, double> crowd_distance_map;
	if (members.size() < 2)
	{
		//throw runtime_error("ParetoObjectives::get_kth_nn_crowding_distance(): less than 2 members");
		crowd_distance_map[members[0]] = 0.0;
		return crowd_distance_map;
	}
	int k = ceil(sqrt(members.size())) - 2; //minus one for zero-based indexing and minus one since the distances will be one less than the num of members
	if (k > members.size())
		k = members.size() - 1;
	map<string, Eigen::VectorXd> obj_member_map;
	
	string m = members[0];
	vector<string> obj_names;
	for (auto obj_map : _member_struct[m])
	{
		obj_names.push_back(obj_map.first);
	}
	Eigen::VectorXd ovals(obj_names.size());
	for (auto member : members)
	{
		crowd_distance_map[member] = 1.0;
		for (int i = 0; i < obj_names.size(); i++)
		{
			ovals[i] = _member_struct[member][obj_names[i]];
		}
		obj_member_map[member] = ovals;
	}

	vector<double> distances;
	double distance;
	for (int i=0;i<members.size();i++)
	{
		distances.clear();
		ovals = obj_member_map[members[i]];
		for (int j = 0; j < members.size(); j++)
		{
			if (j == i)
				continue;
			distance = (ovals - obj_member_map[members[j]]).squaredNorm();
			distances.push_back(distance);
		}
		sort(distances.begin(), distances.end());
		crowd_distance_map[members[i]] = sqrt(distances[k]);	
	}
	return crowd_distance_map;
}


map<string, double> ParetoObjectives::get_cuboid_crowding_distance(ObservationEnsemble& op, ParameterEnsemble& dp)
{
	//this updates the complicated map-based structure that stores the member names: obj_names:value nested pairs
	map<string, map<string, double>> _member_struct = get_member_struct(op, dp);
	//just reuse the same routine used for the pareto dominance sort...
	return get_cuboid_crowding_distance(_member_struct);
}

map<string, double> ParetoObjectives::get_cuboid_crowding_distance(map<string, map<string, double>>& _member_struct)
{
	vector<string> members;
	for (auto m : _member_struct)
		members.push_back(m.first);
	return get_cuboid_crowding_distance(members, _member_struct);
}

map<string, double> ParetoObjectives::get_cuboid_crowding_distance(vector<string>& members)
{
	return get_cuboid_crowding_distance(members, member_struct);
}

map<string, double> ParetoObjectives::get_cuboid_crowding_distance(vector<string>& members, map<string, map<string, double>>& _member_struct)
{

	map<string, map<string, double>> obj_member_map;
	map<string, double> crowd_distance_map;
	string m = members[0];
	vector<string> obj_names;
	//for (auto obj_map : _member_struct[m])
	//{
	//	obj_member_map[obj_map.first] = map<string, double>();
	//	obj_names.push_back(obj_map.first);
	//}



	for (auto member : members)
	{
		crowd_distance_map[member] = 0.0;
		/*for (auto obj_map : _member_struct[member])
			obj_member_map[obj_map.first][member] = obj_map.second;*/

		for (auto obj_map : *obs_obj_names_ptr) //need to make sure the SDs are not included
			obj_member_map[obj_map][member] = _member_struct[member][obj_map];
	}

	//map<double,string>::iterator start, end;
	map<string, double> omap;
	double obj_range;

	for (auto obj_map : obj_member_map)
	{
		omap = obj_map.second;
		//note: for members with identical distances, only the first one gets into the 
		//sorted set but this is ok since we initialized the distance map with zeros
		//for all members, so it works out...
		sortedset crowd_sorted(omap.begin(), omap.end(), compFunctor);

		sortedset::iterator start = crowd_sorted.begin(), last = prev(crowd_sorted.end(), 1);

		obj_range = last->second - start->second;

		//the obj extrema - makes sure they are retained 
		crowd_distance_map[start->first] = CROWDING_EXTREME;
		crowd_distance_map[last->first] = CROWDING_EXTREME;
		if (crowd_sorted.size() == 3)
		{
			sortedset::iterator it = next(start, 1);
			crowd_distance_map[it->first] = crowd_distance_map[it->first] + ((last->second - start->second) / obj_range);

		}
		else if (crowd_sorted.size() > 3)
		{
			//need iterators to start and stop one off from the edges
			start = next(crowd_sorted.begin(), 1);
			last = prev(crowd_sorted.end(), 1);

			sortedset::iterator it = start;

			sortedset::iterator inext, iprev;
			for (; it != last; ++it)
			{
				iprev = prev(it, 1);
				inext = next(it, 1);
				crowd_distance_map[it->first] = crowd_distance_map[it->first] + ((inext->second - iprev->second) / obj_range);
			}
		}
	}
	return crowd_distance_map;
}

vector<double> ParetoObjectives::get_euclidean_distance(map<string, double> first, map<string, double> second)
{
	vector<double> euclidean_dist{ 0, 0 };

	for (auto& obj : *obs_obj_names_ptr)
		euclidean_dist.at(0) += pow(first.at(obj) - second.at(obj), 2);

	if (prob_pareto)
	{
		for (auto& obj : *obs_obj_names_ptr) {
            euclidean_dist.at(1) += 4 * pow(first.at(obj) - second.at(obj), 2) *
                                    (pow(first.at(ppd_obj_to_sd_ptr->at(obj)), 2) + pow(second.at(ppd_obj_to_sd_ptr->at(obj)), 2));
        }
            for (auto& objsd : *obs_obj_sd_names_ptr)
		{
			euclidean_dist.at(0) += pow(first.at(objsd), 2) + pow(second.at(objsd), 2);
			euclidean_dist.at(1) += 2 * pow(pow(first.at(objsd), 2) + pow(second.at(objsd), 2), 2);
		}

	}

	return euclidean_dist;
}

double ParetoObjectives::get_euclidean_fitness(double E, double V)
{
	//double beta = pest_scenario.get_pestpp_options().get_mou_fit_beta();
	double val;

	/*if (beta < 0)
		val = pow(E / exp(pow(V, 0.5)),0.5);
	else
		val = pow(E / (beta * pow(V, 0.5) + 1),0.5);*/
	
	return val;
}

map<string, double> ParetoObjectives::get_cluster_crowding_fitness(vector<string>& members)
{
	return get_cluster_crowding_fitness(members, member_struct);
}

map<string, double> ParetoObjectives::get_cluster_crowding_fitness(vector<string>& members, map<string, map<string, double>>& _member_struct)
{

	map<string, map<string, double>> obj_member_map;
	map<string, double> nondomprob_map, fit_map;
	
	vector<string> obj_names;

	for (auto obj_map : *obs_obj_names_ptr)
		obj_names.push_back(obj_map);

	for (auto member : members)
	{
		nondomprob_map[member] = 0.0;
		fit_map[member] = 0.0;

		for (auto obj_map : *obs_obj_names_ptr) //need to make sure the SDs are not included
			obj_member_map[obj_map][member] = _member_struct[member][obj_map];
	}

	//map<double,string>::iterator start, end;
	map<string, double> omap;
	double fitness;
	vector<double> nonuniq_obj;

	for (auto obj_map : obj_member_map)
	{
		omap = obj_map.second;
		//note: for members with identical distances, only the first one gets into the 
		//sorted set but this is ok since we initialized the distance map with zeros
		//for all members, so it works out...
		sortedset crowd_sorted(omap.begin(), omap.end(), compFunctor);

		sortedset::iterator start = crowd_sorted.begin(), last = prev(crowd_sorted.end(), 1);

		if (members.size() <= pest_scenario.get_pestpp_options().get_mou_max_archive_size())
			min_sd[obj_map.first] = (last->second - start->second) / (members.size());
		else
			min_sd[obj_map.first] = (last->second - start->second) / (pest_scenario.get_pestpp_options().get_mou_max_archive_size());

		for (auto m : members)
		{

            if (_member_struct.at(m).at(ppd_obj_to_sd_ptr->at(obj_map.first)) < min_sd.at(obj_map.first) - FLOAT_EPSILON)
                _member_struct[m][ppd_obj_to_sd_ptr->at(obj_map.first)+"_SYN"] = min_sd.at(obj_map.first);
            else
                _member_struct[m][ppd_obj_to_sd_ptr->at(obj_map.first)+"_SYN"] = _member_struct.at(m).at(ppd_obj_to_sd_ptr->at(obj_map.first));
		}

		nonuniq_obj.clear();

		if (omap.size() != crowd_sorted.size())
		{
			for (auto o : omap)
			{
				int count = 0;
				for (auto cd : crowd_sorted)
					if ((cd.second == o.second) && (cd.first != o.first))
						count++;

				if (count > 0)
					nonuniq_obj.push_back(o.second);
			}
		}

		if (nonuniq_obj.size() != 0)
		{
			for (auto m : members)
			{
				if (find(nonuniq_obj.begin(), nonuniq_obj.end(), _member_struct[m][obj_map.first]) != nonuniq_obj.end())
					fit_map[m] = -1;
			}
		}

		if (iter > 0)
		{
			map<string, map<string, double>> lower_extreme_candidates, upper_extreme_candidates;
			map<string, double> curr;
			vector<string> all_extreme_set;

			for (auto m : incumbent_front_extreme)
				curr[m.first] = incumbent_front_extreme[m.first][obj_map.first];

			sortedset curr_sorted(curr.begin(), curr.end(), compFunctor);
			double lb = incumbent_front_extreme[curr_sorted.begin()->first][obj_map.first];
			//double ub = incumbent_front_extreme[prev(curr_sorted.end(), 1)->first][obj_map.first];

			//get the lower and upper objective bounds of the incumbent true front
			for (auto m : members)
			{
				if (_member_struct[m][obj_map.first] < lb)
					lower_extreme_candidates[m] = _member_struct[m];
			}

			//assign the lower extreme in current population
			if (lower_extreme_candidates.size() == 0) //if size is 0, the extreme is an incumbent front member
				fit_map[start->first] = CROWDING_EXTREME;

			else //use ei to assign the end member
			{
				double mx = 0, ei;
				string endmem;
				for (auto l : lower_extreme_candidates)
				{
					ei = get_ei(_member_struct[l.first], obj_map.first, lb);
					if (ei > mx)
					{
						mx = ei;
						endmem = l.first;
					}
				}

				if (mx == 0)
					fit_map[start->first] = CROWDING_EXTREME;
				else
				{
					if (find(nonuniq_obj.begin(), nonuniq_obj.end(), lower_extreme_candidates[endmem][obj_map.first]) != nonuniq_obj.end())
					{
						map<string, double> ext_mems_pd;
						for (auto l : lower_extreme_candidates)
						{
							if (lower_extreme_candidates[l.first][obj_map.first] != lower_extreme_candidates[endmem][obj_map.first])
								continue;

							double pd = 1;
							for (auto m : lower_extreme_candidates)
							{
								if (l.first != m.first)
									pd *= dominance_probability(_member_struct[l.first], _member_struct[m.first]);
							}
							ext_mems_pd[l.first] = pd;
						}

						double mx = 0;
						string extreme_member_name;
						for (auto em : ext_mems_pd)
						{
							if (em.second > mx)
							{
								mx = em.second;
								extreme_member_name = em.first;
							}
						}
						for (auto em : ext_mems_pd)
						{
							if (em.first == extreme_member_name)
								fit_map[em.first] = CROWDING_EXTREME;

						}
					}
					else
						fit_map[endmem] = CROWDING_EXTREME;
				}
			}
		}
		else
			fit_map[start->first] = CROWDING_EXTREME;
	}

	//crowding distance calculation for non extreme members;
	int nn = 1, nn_max = pest_scenario.get_pestpp_options().get_mou_max_nn_search();
		
	double gamma1 = pest_scenario.get_pestpp_options().get_mou_fit_gamma();
	if (gamma1 == 0)
	{
		double epsilon = pest_scenario.get_pestpp_options().get_mou_fit_epsilon();
		gamma1 = 1 - ppd_beta - epsilon;
	}

	double gamma2 = pow (gamma1, 2);
	double pd, dp, nn_count;
	for (auto m : members)
	{		
		nn_count = 0;
		for (auto n : members)
		{
			if (m != n)
			{
				pd = dominance_prob_adhoc(_member_struct[n], _member_struct[m]);
				dp = dominance_prob_adhoc(_member_struct[m], _member_struct[n]);
				if ((pd > gamma1) && (dp > gamma2))
					nn_count += 1.0;
			}
		}

		nn_map[m] = nn_count;

		if (fit_map[m] != CROWDING_EXTREME)
			fit_map[m] = nn_count;
		
	}

	return fit_map;
}

vector<string> ParetoObjectives::sort_members_by_crowding_distance(int front, vector<string>& members, map<string,double>& crowd_map, 
	map<string, map<string, double>>& _member_struct)
{

	pair<map<string, double>, map<string, double>> euclidean_maps;
	map<string, double> expected_dist_map;
	map<string, double> var_dist_map;
	map<string, double> fit_map;

	if (prob_pareto)
	{
		fit_map = get_mopso_fitness(members, _member_struct);
	}
	else
		fit_map = get_cuboid_crowding_distance(members, _member_struct);

	vector <pair<string, double>> cs_vec;
	for (auto cd : fit_map)
	{
		cs_vec.push_back(cd);
		/*pair <string, double> pnd {cd.first, probnondom_map[cd.first]};
		cs_vec.push_back(pnd);*/
		crowd_map[cd.first] = cd.second;

		/*expected_crowd_map[cd.first] = expected_dist_map[cd.first];
		var_crowd_map[cd.first] = var_dist_map[cd.first];*/
		fitness_map[cd.first] = fit_map[cd.first];
	}

	std::sort(cs_vec.begin(), cs_vec.end(),
		compFunctor);

	reverse(cs_vec.begin(), cs_vec.end());

	vector<string> crowd_ordered;
	for (auto cs : cs_vec)
		crowd_ordered.push_back(cs.first);

	//TODO: check here that all solutions made it thru the crowd distance sorting
	if (crowd_ordered.size() != members.size())
		throw runtime_error("ParetoObjectives::sort_members_by_crowding_distance() error: final sort size != initial size");
	return crowd_ordered;
}

pair<map<string, double>, map<string, double>> ParetoObjectives::get_spea2_fitness(map<string, map<string, double>>& _member_struct)
{
	performance_log->log_event("get_spea_fitness");
	map<string, double> kdist = get_spea2_kth_nn_crowding_distance(_member_struct);
	map<string, vector<string>> solutions_dominated_map;
	map<string, int> num_dominating_map;
	int dom;
	fill_domination_containers(_member_struct, solutions_dominated_map, num_dominating_map,true);
	map<string, double> _unconstrained_fitness_map;
	map<string, double> _fitness_map;
	double max_infeas = -1.0e+30;
	for (auto& in : infeas)
		max_infeas = max(max_infeas, in.second);

	for (auto& sol_map : solutions_dominated_map)
	{
		dom = 0;
		for (auto& sol : sol_map.second)
			dom = dom + num_dominating_map[sol];
		_unconstrained_fitness_map[sol_map.first] = (double)dom + (1.0/(kdist[sol_map.first] + 2.0)); //convert the distance to density
		_fitness_map[sol_map.first] = _unconstrained_fitness_map[sol_map.first];

		//include scaled infeasibility sum in fitness...
		if (infeas.find(sol_map.first) != infeas.end())
			_fitness_map[sol_map.first] += (infeas[sol_map.first] / (max_infeas + 2.0)); // +2.0 to make sure it is less than 0.5, similar to distance penalty
	}
	return pair<map<string, double>, map<string, double>>(_fitness_map,_unconstrained_fitness_map);

}

void ParetoObjectives::fill_domination_containers(map<string, map<string, double>>& _member_struct, map<string,
	vector<string>>& solutions_dominated_map, map<string, int>& num_dominating_map, bool dup_as_dom)
{
	solutions_dominated_map.clear();
	num_dominating_map.clear();
	int domination_counter;
	vector<string> solutions_dominated, first_front;
	performance_log->log_event("fill domination containers");
		for (auto solution_p : _member_struct)
		{
			domination_counter = 0;
			solutions_dominated.clear();
			for (auto solution_q : _member_struct)
			{
				if (solution_p.first == solution_q.first) //string compare real name
					continue;

				//if the solutions are identical...
				if (first_equals_second(solution_p.second, solution_q.second))
				{
					if (dup_as_dom)
						domination_counter++;
					else
						throw runtime_error("ParetoObjectives::fill_domination_containers(): solution '" + solution_p.first + "' and '" + solution_q.first + "' are identical");
				}
				else if (first_dominates_second(solution_p.second, solution_q.second))
				{
					solutions_dominated.push_back(solution_q.first);
				}
				else if (first_dominates_second(solution_q.second, solution_p.second))
				{
					domination_counter++;
				}
			}
			num_dominating_map[solution_p.first] = domination_counter;
			solutions_dominated_map[solution_p.first] = solutions_dominated;
		}
}

map<int,vector<string>> ParetoObjectives::sort_members_by_dominance_into_fronts(map<string, map<string, double>>& _member_struct)
{
	//following fast non-dom alg in Deb
	performance_log->log_event("starting 'fast non-dom sort");
	//map<string,map<string,double>> Sp, F1;
	map<string, vector<string>> solutions_dominated_map;
	map<string, int> num_dominating_map;
	fill_domination_containers(_member_struct, solutions_dominated_map, num_dominating_map); 
	vector<string> solutions_dominated, first_front;
	performance_log->log_event("finding first front");
	for (auto  num_dom : num_dominating_map)
	{	
		//solution_p is in the first front
		if (num_dom.second == 0)
		{
			first_front.push_back( num_dom.first);
		}
	}
	performance_log->log_event("sorting remaining fronts");
	int i = 1;
	int nq;
    vector<string> q_front;
	map<int, vector<string>> front_map;
	front_map[1] = first_front;
	vector<string> front = first_front;
    vector<string> sorted = first_front;
	int num_front_solutions = front.size();

	while (true)
	{
		q_front.clear();
		for (auto solution_p : front)
		{
			solutions_dominated = solutions_dominated_map[solution_p];
			for (auto solution_q : solutions_dominated)
			{
				/*if (!first_dominates_second(_member_struct[solution_q], _member_struct[solution_p])) num_dominating_map[solution_q] = 0;
				else num_dominating_map[solution_q]--;
				*/
				num_dominating_map[solution_q]--;
				if (num_dominating_map[solution_q] <= 0)
					q_front.push_back(solution_q);
			}
		}
		if (q_front.size() == 0)
			break;
		i++;
		front_map[i] = q_front;

		num_front_solutions += q_front.size();
		if (num_front_solutions > _member_struct.size())
        {
		    cout << "note: nsga-ii front sorting: number of visited solutions " << num_front_solutions << " >= number of members " << _member_struct.size() << endl;
		    cout << "q_front:" <<endl;
		    for (auto& f : q_front)
		        cout << "  " << f << endl;
            cout << "front:" <<endl;
            for (auto& f : front)
                cout << "  " << f << endl;
            throw runtime_error("error in nsga-ii front sorting: number of visited solutions > number of members");

            break;
        }
        front = q_front;
//		for (auto& sol : front)
//        {
//		    sorted.push_back(sol);
//        }
//
//		//do we really need more than 10 front for multiobjectives - I dont think so...
//		if ((i > 10) && ((*obs_obj_names_ptr + *pi_obj_names_ptr) > 1))
//        {
//		    set<string> ssorted(sorted.begin(),sorted.end());
//		    set<string>::iterator end = ssorted.end();
//		    q_front.clear();
//		    for (auto& member : _member_struct)
//            {
//		        if (ssorted.find(member.first) == end)
//                {
//		            q_front.push_back(member.first);
//                }
//            }
//		    i++;
//		    front_map[i] = q_front;
//		    break;
//        }


	}
	
	if (num_front_solutions != _member_struct.size())
	{
		stringstream ss;
		ss << "ERROR: ParetoObjectives::sort_members_by_dominance_into_fronts(): number of solutions in fronts (";
		ss << num_front_solutions << ") != member_struct.size() (" << _member_struct.size() << "," << endl;
		file_manager.rec_ofstream() << ss.str();
		cout << ss.str();
		throw runtime_error(ss.str());
	}

	return front_map;
}

//compute probability of dominance
double ParetoObjectives::dominance_probability(map<string, double>& first, map<string, double>& second)
{
	double prob_dom = 1;

	for (auto obj_name : *obj_names_ptr)
	{
		prob_dom *= (std_norm_df(0, first.at(obj_name) - second.at(obj_name), sqrt(pow(first.at(ppd_obj_to_sd_ptr->at(obj_name)),2) + pow(second.at(ppd_obj_to_sd_ptr->at(obj_name)),2)), true));
	}

	return prob_dom;
}

double ParetoObjectives::dominance_prob_adhoc(map<string, double>& first, map<string, double>& second)
{
	map<string, double> f = first, s = second;
	for (auto obj_name : *obj_names_ptr)
	{

        if (f.at(ppd_obj_to_sd_ptr->at(obj_name)) < min_sd.at(obj_name) - FLOAT_EPSILON)
            f.at(ppd_obj_to_sd_ptr->at(obj_name)) = min_sd.at(obj_name);


        if (s.at(ppd_obj_to_sd_ptr->at(obj_name)) < min_sd.at(obj_name) - FLOAT_EPSILON)
            s.at(ppd_obj_to_sd_ptr->at(obj_name)) = min_sd.at(obj_name);
	}

	double pd = dominance_probability(f, s);

	return pd;
}

bool ParetoObjectives::first_equals_second(map<string, double>& first, map<string, double>& second)
{
	for (auto f : first)
	{
		if (abs(f.second - second[f.first]) >= FLOAT_EPSILON)
			return false;
	}
	return true;
}

bool ParetoObjectives::first_dominates_second(map<string, double>& first, map<string, double>& second)
{

	if (ppd_sort)
	{
		double pd = dominance_probability(first, second);

		if (pd < ppd_beta - FLOAT_EPSILON) {
			return false;
		}
		else
			return true;
	}
	else
	{
		for (auto f : first)
		{
			if (f.second > second[f.first] + FLOAT_EPSILON)
				return false;
		}
		return true;
	}

}

double ParetoObjectives::std_norm_df(double x, double mu, double sd, bool cumdf)
{
	double Z, val;
	const double sqrt_2 = sqrt(2.0);
	const double inv_sqrt_2pi = 0.3989422804;

	if (sd == 0)
		sd = 1E-30;

	Z = (x - mu) / sd;

	if (cumdf)
		val = 0.5 * (1 + erf((x - mu)/(sqrt_2*sd)));
	else
		val = inv_sqrt_2pi*exp(-0.5 * Z * Z);

	return val;
}

double ParetoObjectives::psi_function(double aa, double bb, double mu, double sd)
{
	double a = std_norm_df(bb, mu, sd, false);
	double b = std_norm_df(bb, mu, sd, true);
	double psi = sd * std_norm_df(bb, mu, sd, false) + (aa - mu) * std_norm_df(bb, mu, sd, true);
	return psi;
}

//hypervolume works for two-objective problems for now
void ParetoObjectives::set_hypervolume_partitions(map<string, map<string, double>> _hv_pts)
{
	stringstream ss;
	ofstream& frec = file_manager.rec_ofstream();
	ss << "ParetoObjectives::set_hypervolume_partitions() for outer pareto archive members";
	performance_log->log_event(ss.str());
		
	map<string, map<string, double>> hv_parts;
	map<int, vector<double>> hypervolume;
	map<string, double> hv_partition;
	vector<double> hpv;
	double hv_extreme = pest_scenario.get_pestpp_options().get_mou_hypervolume_extreme();
	
	//initialize reference values
	int mult = 1;
	vector <string> ref_tags {"r_0", "rfty"};
	for (string reftags : ref_tags)
	{
		for (auto obj_map : *obj_names_ptr)
		{
			hv_parts[reftags][obj_map] = hv_extreme * mult;
			mult *= -1;
		}
		mult *= -1;

		hv_partition[reftags] = hv_parts[reftags][obj_names_ptr->at(1)];
	}
	
	//set partition boundaries: rectangular strips along obj 2
	for (auto member : _hv_pts)
	{
		for (auto obj_map : *obj_names_ptr)
		{
			hv_parts[member.first][obj_map] = _hv_pts[member.first][obj_map]; //get only objective values, leave SDs, for easy referencing later
		}
		hv_partition[member.first] = _hv_pts[member.first][obj_names_ptr->at(1)];
	}

	//order points by increasing obj 2 values
	sortedset hv_parts_sorted(hv_partition.begin(), hv_partition.end(), compFunctor); 

	int i = 0;
	for (auto hv : hv_parts_sorted)
	{
		hpv.clear();
		for (auto obj_map : *obj_names_ptr)
		{
			hpv.push_back(hv_parts[hv.first][obj_map]);
		}
		hypervolume[i] = hpv;
		i++;
	}
	
	hypervolume_partitions = hypervolume;

	//get the extreme points of the incumbent front to be used for identifying extreme points of pareto cloud thru ehvi
	map<string, double> obj_map;
	for (auto pts : _hv_pts)
		obj_map[pts.first] = _hv_pts[pts.first][obj_names_ptr->at(1)];

	sortedset sortedpts(obj_map.begin(), obj_map.end(), compFunctor);

	sortedset::iterator start = sortedpts.begin(), end = prev(sortedpts.end(), 1);
	incumbent_front_extreme[start->first] = _hv_pts[start->first];
	incumbent_front_extreme[end->first] = _hv_pts[end->first];
}

//this works only for two objectives following the method of Yang et al (2019)
void ParetoObjectives::get_ehvi(ObservationEnsemble& op, ParameterEnsemble& dp)
{
	map<string, map<string, double>> _member_struct = get_member_struct(op, dp);
	string member;
	double ehvi_m;

	for (auto m : _member_struct)
	{
		member = m.first;
		ehvi_m = get_ehvi(member, _member_struct);
		ehvi_member_map[member] = ehvi_m;
	}
	
}

double ParetoObjectives::get_ei(map<string, double> phi, string obj, double curr_opt)
{
    double stdnorm = std_norm_df(curr_opt, phi.at(obj), phi.at(ppd_obj_to_sd_ptr->at(obj)),false);
    double ei = (curr_opt - phi[obj]) * std_norm_df(curr_opt, phi.at(obj), phi.at(ppd_obj_to_sd_ptr->at(obj)), true) + phi.at(ppd_obj_to_sd_ptr->at(obj)) * stdnorm;
	return ei;
}

map<string, double> ParetoObjectives::get_ehvi(vector<string>& members)
{
	map<string, double> ehvi_map;
	double ehvi_mem;

	for (auto m : members)
	{
		ehvi_mem = get_ehvi(m, member_struct);
		ehvi_map[m] = ehvi_mem;
	}

	return ehvi_map;
}

double ParetoObjectives::get_ehvi(string& member, map<string, map<string, double>>& _member_struct)
{
	map<int, vector<double>> hv_i = hypervolume_partitions;
	stringstream ss;
	double t1, t2, p1, p2, p3, ehvi=0;
	vector<double> obj, obj_sd;

	obj.clear();
	obj_sd.clear();
	for (auto& obj_map : *obj_names_ptr)
	{
		obj.push_back(_member_struct.at(member).at(obj_map));
		obj_sd.push_back(_member_struct.at(member).at(ppd_obj_to_sd_ptr->at(obj_map)));
	}

	t1 = 0;
	t2 = 0;
	ehvi = 0;

	map<int, vector<double>>::iterator it = next(hv_i.begin(), 1), iprev;
	for (; it != hv_i.end(); it++)
	{
		iprev = prev(it, 1);

		p1 = hv_i[iprev->first][0] - hv_i[it->first][0];
		p2 = std_norm_df(hv_i[it->first][0], obj.at(0), obj_sd.at(0), true);
		p3 = psi_function(hv_i[it->first][1], hv_i[it->first][1], obj.at(1), obj_sd.at(1));
		t1 += p1 * p2 * p3;

		p1 = psi_function(hv_i[iprev->first][0], hv_i[iprev->first][0], obj.at(0), obj_sd.at(0));
		p2 = psi_function(hv_i[iprev->first][0], hv_i[it->first][0], obj.at(0), obj_sd.at(0));
		p3 = psi_function(hv_i[it->first][1], hv_i[it->first][1], obj.at(1), obj_sd.at(1));
		t2 += (p1 - p2) * p3;
	}

	ehvi = t1 + t2;

	if (ehvi < -FLOAT_EPSILON) //Sometimes the value is only a little bit negative. Perhaps, due to the approximation of std normal. This happened only few times, though, but when it does, temporarily set the value to 0. Will revisit this later.
	{
		ss.str("");
		ss << "WARNING: EHVI of " << member << " is negative = " << ehvi << ". Resetting to 0.0.";
		performance_log->log_event(ss.str());
        cout << ss.str() << endl;
		ehvi = 0;
	}

	//if (ehvi <= -0.1) //If it is way too negative, something must be really wrong.
	//{
	//	ss.str("");
	//	ss << "EHVI of " << member << " is negative: " << ehvi;
	//	performance_log->log_event(ss.str());
	//	throw runtime_error(ss.str());
	//}

	return ehvi;
}


MOEA::MOEA(Pest &_pest_scenario, FileManager &_file_manager, OutputFileWriter &_output_file_writer, 
	PerformanceLog *_performance_log, RunManagerAbstract* _run_mgr_ptr)
	: pest_scenario(_pest_scenario), file_manager(_file_manager),
	output_file_writer(_output_file_writer), performance_log(_performance_log),
	run_mgr_ptr(_run_mgr_ptr), constraints(_pest_scenario, &_file_manager, _output_file_writer, *_performance_log),
	objectives(_pest_scenario,_file_manager,_performance_log)
	
{
	rand_gen = std::mt19937(pest_scenario.get_pestpp_options().get_random_seed());
	dp.set_rand_gen(&rand_gen);
	dp.set_pest_scenario(&pest_scenario);
	op.set_rand_gen(&rand_gen);
	op.set_pest_scenario_ptr(&pest_scenario);
	dp_archive.set_rand_gen(&rand_gen);
	dp_archive.set_pest_scenario(&pest_scenario);
	op_archive.set_rand_gen(&rand_gen);
	op_archive.set_pest_scenario_ptr(&pest_scenario);
}


template<typename T, typename A>
void MOEA::message(int level, const string& _message, vector<T, A> _extras, bool echo)
{
	stringstream ss;
	if (level == 0)
		ss << endl << "  ---  ";
	else if (level == 1)
		ss << "...";
	ss << _message;
	if (_extras.size() > 0)
	{

		for (auto& e : _extras)
			ss << e << " , ";

	}
	if (level == 0)
		ss << "  ---  ";
	if (echo)// && ((verbose_level >= 2) || (level < 2)))
		cout << ss.str() << endl;
	file_manager.rec_ofstream() << ss.str() << endl;
	performance_log->log_event(_message);

}

void MOEA::message(int level, const string& _message)
{
	message(level, _message, vector<string>());
}

template<typename T>
void MOEA::message(int level, const string& _message, T extra)
{
	stringstream ss;
	ss << _message << " " << extra;
	string s = ss.str();
	message(level, s);
}

void MOEA::throw_moea_error(const string& message)
{
	performance_log->log_event("MOEA error: " + message);
	cout << endl << "   ************   " << endl << "    MOEA error: " << message << endl << endl;
	file_manager.rec_ofstream() << endl << "   ************   " << endl << "    MOEA error: " << message << endl << endl;
	file_manager.close_file("rec");
	performance_log->~PerformanceLog();
	throw runtime_error("MOEA error: " + message);
}

int MOEA::get_max_len_obj_name()
{
	int max_len = 20;
	for (auto obs_obj : obs_obj_names)
	{
		max_len = max(int(obs_obj.size()), max_len);
	}
	for (auto pi_obj : pi_obj_names)
	{
		max_len = max(int(pi_obj.size()), max_len);
	}
	return max_len;
}

map<string, map<string, double>> MOEA::decvar_report(ParameterEnsemble& _dp)
{

    stringstream ss;
    ss.str("");
    int max_len = 19;
    for (auto& dv : dv_names)
    {
        max_len = max(max_len,(int)dv.size());
    }
    ss << left << setw(max_len) << "decision variable " << right << setw(10) << "ubnd " << setw(10) << "lbnd " << setw(12) << "mean " << setw(12);
    ss << "stdev " << setw(12) << "min " << setw(12) << "max " << endl;
    map<string,double> meanmap,stdmap;
    _dp.fill_moment_maps(meanmap,stdmap);
    Eigen::VectorXd vec;
    const ParameterRec* prec;
    map<string,map<string,double>> summary;
    map<string,double> sum;
    for (auto& dv : dv_names)
    {
        sum.clear();
        prec = pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(dv);
        ss << left << setw(max_len) << dv << " ";
        ss << right << setw(10) << prec->ubnd << " ";
        ss << right << setw(10) << prec->lbnd << " ";
        ss << right << setw(12) << meanmap[dv] << " ";
        ss << right << setw(12) << stdmap[dv] << " ";
        vec = _dp.get_var_vector(dv);
        ss << setw(12) << vec.minCoeff() << " ";
        ss << setw(12) << vec.maxCoeff();
        ss << endl;
        sum["mean"] = meanmap[dv];
        sum["std"] = stdmap[dv];
        sum["min"] = vec.minCoeff();
        sum["max"] = vec.maxCoeff();
        summary[dv] = sum;
    }
    ss << endl;
    file_manager.rec_ofstream() << ss.str();
    cout << ss.str();
    return summary;
}

map<string, map<string, double>> MOEA::decvar_change_report(map<string, map<string, double>>& current_dv_summary)
{
    map<string, map<string, double>> change_summary;
    if (previous_dv_summary.size() == 0)
        return change_summary;
    double change, percent_change;

    stringstream ss;
    int max_len = 19;
    for (auto& n : dv_names)
        max_len = max(max_len,(int)n.size());

    ss << left << setw(max_len) << "decision variable " << right << setw(11) << "mean change ";
    ss << setw(11) << "% change ";
    ss << setw(11) << "max change " << setw(11) << "% change ";
    ss << setw(11) << "min change " << setw(11) << "% change " << endl;

    vector<string> tags{ "mean","max","min" };
    for (auto dv : dv_names)
    {
        change_summary[dv] = map<string, double>();
        if (previous_dv_summary.find(dv) == previous_dv_summary.end())
            throw_moea_error("decvar_change_report() error: dv '" + dv + "' not in previous summary");
        if (current_dv_summary.find(dv) == current_dv_summary.end())
            throw_moea_error("decvar_change_report() error: dv '" + dv + "' not in current summary");
        ss << left << setw(max_len) << dv << " ";
        for (auto tag : tags)
        {

            change = previous_dv_summary[dv][tag] - current_dv_summary[dv][tag];
            if ((previous_dv_summary[dv][tag] <= 0.0) || (change == 0.0))
                percent_change = 0.0;
            else
                percent_change = 100.0 * (change / previous_dv_summary[dv][tag]);

            ss << right << setw(11) << change << " ";
            ss << setw(11) << percent_change << " ";
            change_summary[dv][tag] = change;
            change_summary[dv][tag+"_percent"] = percent_change;
        }
        ss << endl;
    }

    file_manager.rec_ofstream() << ss.str() << endl;
    cout << ss.str() << endl;
    return change_summary;
}

map<string, map<string, double>> MOEA::obj_func_report(ParameterEnsemble& _dp, ObservationEnsemble& _op)
{
	map<string, map<string, double>> summary = get_obj_func_summary_stats(_dp, _op);
	stringstream frec;
	//frec << endl << "  ---  Objective Function Summary  ---  " << endl;

	int max_len = get_max_len_obj_name();
	string dir;
	frec << left << setw(max_len) << "objective function" << right << setw(10) << "direction" << setw(13) << "mean";
	frec << setw(13) << "std dev" << setw(13) << "min" << setw(13) << "max" << setw(13) << "knee" << endl;
    frec << endl << setprecision(7);
	for (auto obs_obj : obs_obj_names)
	{

		frec << left << setw(max_len) << obs_obj;
		dir = "minimize";
		if (obj_dir_mult[obs_obj] == -1)
			dir = "maximize";

		frec << right << setw(10) << dir;
		frec << right << setw(13) << summary[obs_obj]["mean"];
		frec << setw(13) << summary[obs_obj]["std"];
		frec << setw(13) << summary[obs_obj]["min"];
		frec << setw(13) << summary[obs_obj]["max"] ;
		frec << setw(13) << summary[obs_obj]["knee"];
		frec << endl;
	}

	
	for (auto pi_obj : pi_obj_names)
	{
		frec << left << setw(max_len) << pi_obj;
		dir = "minimize";
		if (obj_dir_mult[pi_obj] == -1)
			dir = "maximize";
		frec << right << setw(10) << dir;
		frec << right << setw(13) << summary[pi_obj]["mean"];
		frec << setw(13) << summary[pi_obj]["std"];
		frec << setw(13) << summary[pi_obj]["min"];
		frec << setw(13) << summary[pi_obj]["max"];
        frec << setw(13) << summary[pi_obj]["knee"];
		frec << endl;
	}



	frec << endl;
	file_manager.rec_ofstream() << frec.str();
	cout << frec.str();
	return summary;
}

map<string, map<string, double>> MOEA::obj_func_change_report(map<string, map<string, double>>& current_obj_summary)
{
	map<string, map<string, double>> change_summary;
	if (previous_obj_summary.size() == 0)
		return change_summary;
	double change, percent_change;
	
	stringstream ss;
	int max_len = get_max_len_obj_name();
	ss << left << setw(max_len) << "objective function" << right << setw(12) << "mean change";
	ss << setw(12) << "% change";
	ss << setw(12) << "max change" << setw(12) << "% change";
	ss << setw(12) << "min change" << setw(12) << "% change";
    ss << setw(12) << "knee change" << setw(12) << "% change";
	ss << endl;

	vector<string> tags{ "mean","max","min","knee" };
	for (auto obs_obj : obs_obj_names)
	{
		change_summary[obs_obj] = map<string, double>();
		if (previous_obj_summary.find(obs_obj) == previous_obj_summary.end())
			throw_moea_error("obj_func_change_report() error: obs obj '" + obs_obj + "' not in previous summary");
		if (current_obj_summary.find(obs_obj) == current_obj_summary.end())
			throw_moea_error("obj_func_change_report() error: obs obj '" + obs_obj + "' not in current summary");
		ss << left << setw(max_len) << obs_obj;
		for (auto tag : tags)
		{
			
			change = previous_obj_summary[obs_obj][tag] - current_obj_summary[obs_obj][tag];
			if ((previous_obj_summary[obs_obj][tag] <= 0.0) || (change == 0.0))
				percent_change = 0.0;
			else
				percent_change = 100.0 * (change / previous_obj_summary[obs_obj][tag]);
			
			ss << right << setw(12) << change;
			ss << setw(12) << percent_change;
			change_summary[obs_obj][tag] = change;
			change_summary[obs_obj][tag+"_percent"] = percent_change;
		}
		ss << endl;
	}
	for (auto pi_obj : pi_obj_names)
	{
		change_summary[pi_obj] = map<string, double>();
		if (previous_obj_summary.find(pi_obj) == previous_obj_summary.end())
			throw_moea_error("obj_func_change_report() error: pi obj '" + pi_obj + "' not in previous summary");
		if (current_obj_summary.find(pi_obj) == current_obj_summary.end())
			throw_moea_error("obj_func_change_report() error: pi obj '" + pi_obj + "' not in current summary");
		ss << left << setw(max_len) << pi_obj;
		for (auto tag : tags)
		{
			change = previous_obj_summary[pi_obj][tag] - current_obj_summary[pi_obj][tag];
			if ((previous_obj_summary[pi_obj][tag] <= 0.0) || (change == 0.0))
				percent_change = 0.0;
			else
				percent_change = 100.0 * (change / previous_obj_summary[pi_obj][tag]);
			
			ss << right << setw(12) << change;
			ss << setw(12) << percent_change;
			change_summary[pi_obj][tag] = change;
			change_summary[pi_obj][tag + "_percent"] = percent_change;
		}
		ss << endl;
	}
	file_manager.rec_ofstream() << ss.str() << endl;
	cout << ss.str() << endl;
	return change_summary;
}

map<string, map<string, double>> MOEA::get_obj_func_summary_stats(ParameterEnsemble& _dp, ObservationEnsemble& _op)
{

	//pair<map<string, double>, map<string, double>> mm = _op.get_moment_maps();
	map<string, double> mean_map, std_map;
	_op.fill_moment_maps(mean_map, std_map);
	_op.update_var_map();
	map<string, int> var_map = _op.get_var_map();
	map<string, double> sum;
	map<string, map<string, double>> summary_stats;
	string opt_member_name;
	pair<Parameters,Observations> p = get_optimal_solution(_dp,_op,opt_member_name);
	for (auto obs_obj : obs_obj_names)
	{
		sum.clear();
		sum["mean"] = mean_map[obs_obj];
		sum["std"] = std_map[obs_obj];
		sum["min"] = _op.get_eigen_ptr()->col(var_map[obs_obj]).minCoeff();
		sum["max"] = _op.get_eigen_ptr()->col(var_map[obs_obj]).maxCoeff();
		sum["knee"] = p.second[obs_obj];
		summary_stats[obs_obj] = sum;
	}
	_dp.update_var_map();
	var_map = _dp.get_var_map();
	vector<string> dp_names = _dp.get_var_names();
	Parameters pars = pest_scenario.get_ctl_parameters();
	ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();
	pts.ctl2numeric_ip(pars);
	_dp.transform_ip(ParameterEnsemble::transStatus::NUM);
	Eigen::VectorXd vec;
	pair<double, double> sim_res;
	map<string, vector<double>> pi_vals;
	PriorInformation* prior_info_ptr = pest_scenario.get_prior_info_ptr();
	for (auto pi_obj : pi_obj_names)
		pi_vals[pi_obj] = vector<double>();
	for (int i = 0; i < _dp.shape().first; i++)
	{
		vec = _dp.get_eigen_ptr()->row(i);
		pts.ctl2numeric_ip(pars);
		pars.update_without_clear(dp_names, vec);
		pts.numeric2ctl_ip(pars);
		
		for (auto pi_obj : pi_obj_names)
		{
			sim_res = prior_info_ptr->get_pi_rec(pi_obj).calc_sim_and_resid(pars);
			pi_vals[pi_obj].push_back(sim_res.first);
		}
	}

	for (auto pi_obj : pi_obj_names)
	{
        sim_res = prior_info_ptr->get_pi_rec(pi_obj).calc_sim_and_resid(p.first);
		vec = stlvec_2_eigenvec(pi_vals[pi_obj]);
		sum.clear();
		sum["mean"] = vec.mean();
		sum["std"] = sqrt((vec.array() - vec.mean()).pow(2).sum() / (vec.size() - 1));
		sum["min"] = vec.minCoeff();
		sum["max"] = vec.maxCoeff();
		sum["knee"] = sim_res.first;
		summary_stats[pi_obj] = sum;
	}

//	//calculate relative hyper volumes
//	//first form the ideal solution vector
//	Eigen::VectorXd ideal(summary_stats.size());
//	int i=0;
//	for (auto& oname : obs_obj_names)
//    {
//	    if (obj_dir_mult[oname] == 1)
//	        ideal[i] = summary_stats.at(oname).at("min");
//	    else
//            ideal[i] = summary_stats.at(oname).at("max");
//	    i++;
//    }
//	for (auto& pname : pi_obj_names)
//    {
//        if (obj_dir_mult[pname] == 1)
//            ideal[i] = summary_stats.at(pname).at("min");
//        else
//            ideal[i] = summary_stats.at(pname).at("max");
//        i++;
//    }
	

	return summary_stats;
}


void MOEA::sanity_checks()
{
	PestppOptions* ppo = pest_scenario.get_pestpp_options_ptr();
	vector<string> errors;
	vector<string> warnings;
	stringstream ss;

	if ((population_dv_file.size() == 0) && (population_obs_restart_file.size() > 0))
		errors.push_back("population_dv_file is empty but population_obs_restart_file is not - how can this work?");
	
	if ((ppo->get_mou_population_size() < error_min_members) && (population_dv_file.size() == 0))
	{
		ss.str("");
		ss << "population_size < " << error_min_members << ", this is redic, increasing to " << warn_min_members;
		warnings.push_back(ss.str());
		ppo->set_ies_num_reals(warn_min_members);
	}
	if ((ppo->get_mou_population_size() < warn_min_members) && (population_dv_file.size() == 0))
	{
		ss.str("");
		ss << "mou_population_size < " << warn_min_members << ", this is prob too few";
		warnings.push_back(ss.str());
	}
	
	bool use_pso = false;
	for (auto gt : gen_types)
		if (gt == MouGenType::PSO)
			use_pso = true;
	if (use_pso)
		if (gen_types.size() > 1)
			errors.push_back("'PSO' generator cannot be used with other generators...");

	
	if (warnings.size() > 0)
	{
		message(0, "sanity_check warnings");
		for (auto& w : warnings)
			message(1, w);
		message(1, "continuing initialization...");
	}
	if (errors.size() > 0)
	{
		message(0, "sanity_check errors - uh oh");
		for (auto& e : errors)
			message(1, e);
		throw_moea_error(string("sanity_check() found some problems - please review rec file"));
	}

	

	//cout << endl << endl;
}


void MOEA::update_archive_nsga(ObservationEnsemble& _op, ParameterEnsemble& _dp)
{
	message(2, "updating archive");
	stringstream ss;
	if (op_archive.shape().first != dp_archive.shape().first)
	{
		ss.str("");
		ss << "MOEA::update_archive_nsga(): op_archive members " << op_archive.shape().first << " != dp_archive members " << dp_archive.shape().first;
		throw_moea_error(ss.str());
	}
	//if this is a population reset because of trying to reuse chances, then we need to reset the archive now also...
    if (should_use_multigen())
    {
        message(2,"resetting archive for multi-generational population");
        dp_archive = _dp;
        op_archive = _op;
        DomPair dompair = objectives.get_nsga2_pareto_dominance(iter, op_archive, dp_archive, &constraints, prob_pareto, true,
                                                                ARC_SUM_TAG);
        dp_archive.keep_rows(dompair.first);
        op_archive.keep_rows(dompair.first);

    }
    else {
        //check that members of _op aren't in the archive already
        vector<string> keep, temp = op_archive.get_real_names();
        set<string> archive_members(temp.begin(), temp.end());
        for (auto &member : _op.get_real_names()) {
            if (archive_members.find(member) == archive_members.end())
                keep.push_back(member);
        }
        if (keep.size() == 0) {
            message(2, "all nondominated members in already in archive");
            return;
        }

        ss.str("");
        ss << "adding " << keep.size() << " non-dominated members to archive";
        message(2, ss.str());
        Eigen::MatrixXd other = _op.get_eigen(keep, vector<string>());
        op_archive.append_other_rows(keep, other);
        other = _dp.get_eigen(keep, vector<string>());
        dp_archive.append_other_rows(keep, other);
        other.resize(0, 0);
        message(2, "pareto dominance sorting archive of size", op_archive.shape().first);
        DomPair dompair = objectives.get_nsga2_pareto_dominance(iter, op_archive, dp_archive, &constraints, prob_pareto, true, ARC_SUM_TAG);

        ss.str("");
        ss << "resizing archive from " << op_archive.shape().first << " to " << dompair.first.size()
           << " current non-dominated solutions";
        message(2, ss.str());
        op_archive.keep_rows(dompair.first);
        dp_archive.keep_rows(dompair.first);
    }

	if (op_archive.shape().first > archive_size)
	{
	    vector<string> keep;
		ss.str("");
		ss << "trimming archive size from " << op_archive.shape().first << " to max archive size " << archive_size;
		message(2, ss.str());
		vector<string> members = op_archive.get_real_names();
		keep.clear();
		for (int i = 0; i < archive_size; i++)
			keep.push_back(members[i]);
		op_archive.keep_rows(keep);
		dp_archive.keep_rows(keep);
	}

	//report only trimmed archive with updated fitness value
	objectives.get_nsga2_pareto_dominance(iter, op_archive, dp_archive, &constraints, prob_pareto, true, ARC_TRIM_SUM_TAG);

	dp_archive.reset_org_real_names();
	op_archive.reset_org_real_names();

	save_populations(dp_archive, op_archive, "archive");
}

void MOEA::update_archive_spea(ObservationEnsemble& _op, ParameterEnsemble& _dp)
{
    int current_archive_size = population_schedule[iter] * 2;
	message(2, "updating archive");
	stringstream ss;
	if (op_archive.shape().first != dp_archive.shape().first)
	{
		ss.str("");
		ss << "MOEA::update_archive_spea(): op_archive members " << op_archive.shape().first << " != dp_archive members " << dp_archive.shape().first;
		throw_moea_error(ss.str());
	}
    if (should_use_multigen())
    {
        message(1,"resetting archive for multigenerational population");
        dp_archive = _dp;
        op_archive = _op;
        vector<string> keep = dp_archive.get_real_names();
        if (op_archive.shape().first > current_archive_size)
        {
            objectives.get_spea2_archive_names_to_keep(current_archive_size, keep, op_archive, dp_archive);
            op_archive.keep_rows(keep);
            dp_archive.keep_rows(keep);
        }
        save_populations(dp_archive, op_archive, "archive");
        return;
    }

    //check that members of _op aren't in the archive already
    vector<string> keep, temp = op_archive.get_real_names();
    set<string> archive_members(temp.begin(), temp.end());
    for (auto &member : _op.get_real_names()) {
        if (archive_members.find(member) == archive_members.end())
            keep.push_back(member);
    }
    if (keep.size() == 0) {
        message(2, "all members in already in archive");
        return;
    }

    ss.str("");
    ss << "adding " << keep.size() << " members to archive";
    message(2, ss.str());
    Eigen::MatrixXd other = _op.get_eigen(keep, vector<string>());
    op_archive.append_other_rows(keep, other);
    other = _dp.get_eigen(keep, vector<string>());
    dp_archive.append_other_rows(keep, other);
    other.resize(0, 0);
    message(2, "spea fitness calculation for archive of size ", op_archive.shape().first);
    map<string, double> fit = objectives.get_spea2_fitness(iter, op_archive, dp_archive, &constraints, true,
                                                           ARC_SUM_TAG);


	keep.clear();
	for (auto& f : fit)
		if (f.second < 1.0)
			keep.push_back(f.first);

	ss.str("");
	ss << "resizing archive from " << op_archive.shape().first << " to " << keep.size() << " current non-dominated solutions";
	message(2, ss.str());
	op_archive.keep_rows(keep);
	dp_archive.keep_rows(keep);
	if (keep.size() > 0)
	{
		if (op_archive.shape().first > current_archive_size)
		{
			objectives.get_spea2_archive_names_to_keep(current_archive_size, keep, op_archive, dp_archive);
			op_archive.keep_rows(keep);
			dp_archive.keep_rows(keep);
		}
        dp_archive.reset_org_real_names();
        op_archive.reset_org_real_names();
        objectives.get_nsga2_pareto_dominance(iter, op_archive, dp_archive, &constraints, prob_pareto, true,
                                                                ARC_SUM_TAG);
		save_populations(dp_archive, op_archive, "archive");
	}
}

void MOEA::queue_resample_runs(ParameterEnsemble& _dp)
{
	//insert outer iter scripts
}


void MOEA::queue_chance_runs(ParameterEnsemble& _dp)
{
	/* queue up chance-related runs using the class attributes dp and op*/
	stringstream ss;
	if (constraints.should_update_chance(iter))
	{
		dp.transform_ip(ParameterEnsemble::transStatus::NUM);
		Parameters pars = pest_scenario.get_ctl_parameters();
		Parameters opars;
		pest_scenario.get_base_par_tran_seq().ctl2numeric_ip(pars);
		Observations obs = pest_scenario.get_ctl_observations();
		//if this is the first iter and no restart
        ss.str("");
        ss << "queuing chance runs for generation " << iter;
		message(1,ss.str());
		if (chancepoints == chancePoints::SINGLE)
		{
			//dont use the _dp, use the class attr dp and op here
			//because they are in sync. _dp hasn't been run yet...
			string opt_member;
			Parameters::iterator end = pars.end();
			pair<Parameters, Observations> po_pair = get_optimal_solution(dp, op, opt_member);
            opars = po_pair.first;
            pest_scenario.get_base_par_tran_seq().ctl2numeric_ip(opars);
			for (auto& item : opars)
				if (pars.find(item.first) != end)
					pars.update_rec(item.first,item.second);

			obs = po_pair.second;
			constraints.add_runs(iter, pars, obs, run_mgr_ptr);
		}
		else if (chancepoints == chancePoints::ALL)
		{
			constraints.add_runs(iter, _dp, obs, run_mgr_ptr);
		}
		 
		else
		{
			throw_moea_error("internal error: trying to call unsupported chancePoints type");
		}

	}
}

vector<int> MOEA::run_population(ParameterEnsemble& _dp, ObservationEnsemble& _op, bool allow_chance)
{
	run_mgr_ptr->reinitialize();
	//queue up any chance related runs
	if (allow_chance)
		queue_chance_runs(_dp);

	//queue up outer iter runs
	/*if (iter % pest_scenario.get_pestpp_options().get_mou_resample_every() == 0)
		queue_resample_runs(_dp);*/

	message(1, "running population of size ", _dp.shape().first);
	stringstream ss;
	ss << "queuing " << _dp.shape().first << " runs";
	performance_log->log_event(ss.str());
	//run_mgr_ptr->reinitialize();
	map<int, int> real_run_ids;
	try
	{
		real_run_ids = _dp.add_runs(run_mgr_ptr);
	}
	catch (const exception& e)
	{
		stringstream ss;
		ss << "run_ensemble() error queueing runs: " << e.what();
		throw_moea_error(ss.str());
	}
	catch (...)
	{
		throw_moea_error(string("run_ensemble() error queueing runs"));
	}
	performance_log->log_event("making runs");
	try
	{
		run_mgr_ptr->run();
	}
	catch (const exception& e)
	{
		stringstream ss;
		ss << "error running ensemble: " << e.what();
		throw_moea_error(ss.str());
	}
	catch (...)
	{
		throw_moea_error(string("error running ensemble"));
	}

	performance_log->log_event("processing runs");
	
	vector<int> failed_real_indices;
	try
	{
		failed_real_indices = _op.update_from_runs(real_run_ids, run_mgr_ptr);
	}
	catch (const exception& e)
	{
		stringstream ss;
		ss << "error processing runs: " << e.what();
		throw_moea_error(ss.str());
	}
	catch (...)
	{
		throw_moea_error(string("error processing runs"));
	}
	//for testing
	if (pest_scenario.get_pestpp_options().get_ies_debug_fail_subset())
		failed_real_indices.push_back(0);

	if (failed_real_indices.size() > 0)
	{
		stringstream ss;
		vector<string> par_real_names = _dp.get_real_names();
		vector<string> obs_real_names = _op.get_real_names();
		ss << "the following par:obs realization runs failed: ";
		for (auto& i : failed_real_indices)
		{
			ss << par_real_names[i] << ":" << obs_real_names[i] << ',';
		}
		performance_log->log_event(ss.str());
		message(1, "failed realizations: ", failed_real_indices.size());
		string s = ss.str();
		message(1, s);
		performance_log->log_event("dropping failed realizations");
		_dp.drop_rows(failed_real_indices);
		_op.drop_rows(failed_real_indices);
	}

	//do this here in case something is wrong, we know sooner than later
	if (allow_chance)
	{
		constraints.process_runs(run_mgr_ptr, iter);
		constraints.update_chance_offsets();
	}
	return failed_real_indices;
}

ObservationEnsemble MOEA::get_chance_shifted_op(ParameterEnsemble& _dp, ObservationEnsemble& _op, string& opt_member)
{
	if (risk_obj)
		return constraints.get_chance_shifted_constraints(_dp, _op, iter, RISK_NAME, opt_member);
	else
		return constraints.get_chance_shifted_constraints(_dp, _op, iter, string(), opt_member);
}

void MOEA::finalize()
{

}

void MOEA::initialize()
{
	stringstream ss;
	ofstream& frec = file_manager.rec_ofstream();
	message(0, "initializing MOEA process");
	
	pp_args = pest_scenario.get_pestpp_options().get_passed_args();

	act_obs_names = pest_scenario.get_ctl_ordered_nz_obs_names();
	act_par_names = pest_scenario.get_ctl_ordered_adj_par_names();

	//define these here to make sure the iter loop later behaves
	iter = 0;
	restart_iter_offset = 0;
	member_count = 0;

	warn_min_members = 20;
	error_min_members = 4;

	initialize_population_schedule();

	//set some defaults
	PestppOptions* ppo = pest_scenario.get_pestpp_options_ptr();

	string env = ppo->get_mou_env_selector();
	if (env == "NSGA")
	{
		envtype = MouEnvType::NSGA;
		prob_pareto = false;
		objectives.set_prob_pareto(prob_pareto);
		message(1, "using 'nsga2' env selector");
	}
	else if (env == "NSGA_PPD")
	{
		envtype = MouEnvType::NSGA;
		prob_pareto = true;
		objectives.set_ppd_beta();
		objectives.set_prob_pareto(prob_pareto);
		message(1, "using 'nsga2_ppd' env selector");
	}
	else if (env == "SPEA")
	{
		envtype = MouEnvType::SPEA;
		message(1, "using 'spea2' env selector");
	}
	else
		throw_moea_error("'mou_env_selector' type not recognized: " + env + ", should be 'NSGA' or 'SPEA'");

	string mate = ppo->get_mou_mating_selector();
	if (mate == "RANDOM")
	{
		mattype = MouMateType::RANDOM;
		message(1, "using random mating pool selector");
	}	
	else if (mate == "TOURNAMENT")
	{
		mattype = MouMateType::TOURNAMENT;
		message(1, "using binary tournament mating pool selector");
	}
	else
		throw_moea_error("'mou_mating_selector' type not recognized: " + mate + ", should be 'RANDOM' or 'TOURNAMENT'");

	save_every = ppo->get_mou_save_population_every();
	if (save_every <= 0)
	{
		message(1, "'mou_save_population_every' less than/equal to zero, not saving generation-specific populations (and archives)");
	}
	else
	{
		message(1, "saving generation specific populations and archives every nth generation", save_every);
	}


	//reset the par bound PI augmentation since that option is just for simplex
	ppo->set_opt_include_bnd_pi(false);

	//process dec var args
	vector<string> dec_var_groups = ppo->get_opt_dec_var_groups();
	if (dec_var_groups.size() != 0)
	{
		//first make sure all the groups are actually listed in the control file
		vector<string> missing;
		vector<string> pst_groups = pest_scenario.get_ctl_ordered_par_group_names();
		vector<string>::iterator end = pst_groups.end();
		vector<string>::iterator start = pst_groups.begin();
		for (auto grp : dec_var_groups)
			if (find(start, end, grp) == end)
				missing.push_back(grp);
		if (missing.size() > 0)
		{
			ss.str("");
			ss << "the following ++opt_dec_var_groups were not found: ";
			for (auto m : missing)
				ss << m << ",";
			throw_moea_error(ss.str());
		}


		//find the parameter in the dec var groups
		ParameterGroupInfo pinfo = pest_scenario.get_base_group_info();
		string group;
		end = dec_var_groups.end();
		start = dec_var_groups.begin();
		for (auto& par_name : pest_scenario.get_ctl_ordered_adj_par_names())
		{
			group = pinfo.get_group_name(par_name);
			if (find(start, end, group) != end)
			{
				dv_names.push_back(par_name);

			}
		}

		if (dv_names.size() == 0)
		{
			ss.str("");
			ss << "no adjustable decision variables found in supplied dec var groups : ";
			for (auto g : dec_var_groups)
			{
				ss << g << ",";
			}
			throw_moea_error(ss.str());
		}
		ss.str("");
		ss << "'opt_dec_var_groups' passed, using " << dv_names.size() << " adjustable parameters as decision variables";
		message(2, ss.str());
		ofstream& frec = file_manager.rec_ofstream();
		frec << "decision variables:" << endl;
		int icol = 0;
		for (auto dv_name : dv_names)
		{
			frec << dv_name << " ";
			icol++;
			if (icol == 10)
			{
				frec << endl;
				icol = 0;
			}
		}
		frec << endl;
	}
	//otherwise, just use all adjustable parameters as dec vars
	else
	{
		message(2, "using all adjustable parameters as decision variables: ", act_par_names.size());
		dv_names = act_par_names;
	}


	message(1, "number of decision variables: ", dv_names.size());
	message(1, "max run fail: ", ppo->get_max_run_fail());


	//some risk-based stuff here
	string chance_points = ppo->get_opt_chance_points();
	if (chance_points == "ALL")
	{
		//evaluate the chance constraints at every individual, very costly, but most robust
		chancepoints = chancePoints::ALL;
		message(1, "'opt_chance_points' = ALL, evaluating chance at all population members");
	}
	
	else if (chance_points == "SINGLE")
	{
		//evaluate the chance constraints only at the population member nearest the optimal tradeoff.
		//much cheaper, but assumes linear coupling
		chancepoints = chancePoints::SINGLE;
		message(1, "'opt_chance_points' = SINGLE, evaluating chance at representative point");
	}
	else
	{
		ss.str("");
		ss << "unrecognized 'opt_chance_points' value :" << chance_points << ", should be 'all' or 'single'";
		throw_moea_error(ss.str());
	}


	//process objectives
	Constraints::ConstraintSense gt = Constraints::ConstraintSense::greater_than, lt = Constraints::ConstraintSense::less_than;
	pair<Constraints::ConstraintSense, string> sense;
	map<string, string> obj_sense_map;
	vector<string> onames = pest_scenario.get_ctl_ordered_nz_obs_names();
	vector<string> passed_obj_names = ppo->get_mou_objectives();
	if (passed_obj_names.size() == 0)
	{
		for (auto oname : onames)
		{
			sense = Constraints::get_sense_from_group_name(pest_scenario.get_ctl_observation_info().get_group(oname));
			if (sense.first == gt)
			{
				obs_obj_names.push_back(oname);
				obj_sense_map[oname] = "maximize";
				obj_dir_mult[oname] = -1.0;
			}
			else if (sense.first == lt)
			{
				obs_obj_names.push_back(oname);
				obj_sense_map[oname] = "minimize";
				obj_dir_mult[oname] = 1.0;
			}

		}

		onames = pest_scenario.get_ctl_ordered_pi_names();
		for (auto oname : onames)
		{
			if (pest_scenario.get_prior_info().get_pi_rec(oname).get_weight() == 0.0)
				continue;
			sense = Constraints::get_sense_from_group_name(pest_scenario.get_prior_info().get_pi_rec(oname).get_group());
			if (sense.first == gt)
			{
				pi_obj_names.push_back(oname);
				obj_sense_map[oname] = "maximize";
				obj_dir_mult[oname] = -1.0;
			}
			else if (sense.first == lt)
			{
				pi_obj_names.push_back(oname);
				obj_sense_map[oname] = "minimize";
				obj_dir_mult[oname] = 1.0;
			}
		}

		message(1, "'mou_objectives' not passed, using all nonzero weighted obs and prior info eqs that use the proper obs group naming convention");
	}
	else
	{
		vector<string> onames = pest_scenario.get_ctl_ordered_nz_obs_names();
		set<string> oset(onames.begin(), onames.end());
		onames = pest_scenario.get_ctl_ordered_pi_names();
		set<string> pinames(onames.begin(), onames.end());
		onames.clear();
		vector<string> missing,keep_obs, keep_pi,err_sense,keep_obs_sd, keep_pi_sd;
        map<string,string> obslink = pest_scenario.get_ext_file_string_map("observation data external","link_to");
		for (auto obj_name : passed_obj_names)
		{
			if ((oset.find(obj_name) == oset.end()) && (pinames.find(obj_name) == pinames.end()))
				missing.push_back(obj_name);
			else if (oset.find(obj_name) != oset.end())
			{
				sense = Constraints::get_sense_from_group_name(pest_scenario.get_ctl_observation_info().get_group(obj_name));
				if ((sense.first != gt) && (sense.first != lt))
					err_sense.push_back(obj_name);
				else
				{
					if (sense.first == gt)
					{
						keep_obs.push_back(obj_name);
						obj_sense_map[obj_name] = "maximize";
						obj_dir_mult[obj_name] = -1.0;
					}
					else if (sense.first == lt)
					{
						keep_obs.push_back(obj_name);
						obj_sense_map[obj_name] = "minimize";
						obj_dir_mult[obj_name] = 1.0;
					}
					obs_obj_names.push_back(obj_name);
				}
				if (prob_pareto)
                {
                    string sdobs = obj_name + "_SD";
                    bool found = true;
                    if (oset.find(sdobs) == oset.end())
                    {
                        found = false;
                        if (obslink.find(obj_name) != obslink.end())
                        {
                            sdobs = obslink.at(obj_name);
                            if (oset.find(sdobs) != oset.end())
                            {
                                found = true;
                            }
                        }
                        if (!found) {
                            ss.str("");
                            ss << "PPD is active but objective '" << obj_name
                               << "' needs the corresponding standard deviation observation: '" << sdobs << "'";
                            throw_moea_error(ss.str());
                        }
                    }
                    message(1,"found PPD standard deviation observation: '"+sdobs+"' for objective: '"+obj_name+"'");

                    keep_obs_sd.push_back(sdobs);
                    if (ppd_obj_to_sd.find(obj_name) != ppd_obj_to_sd.end())
                    {
                        ss.str("");
                        ss << "objective '" << obj_name << "' already in ppd_obj_to_sd map";
                        throw_moea_error(ss.str());
                    }
                    ppd_obj_to_sd[obj_name] = sdobs;
                }
			}
			//else if (oset.find(obj_name+"_sd") != oset.end()) //find the corresponding sd observations
			//{
			//	keep_obs_sd.push_back(obj_name + "_sd");
			//}
			else
			{
				sense = Constraints::get_sense_from_group_name(pest_scenario.get_prior_info().get_pi_rec(obj_name).get_group());
				if ((sense.first != gt) && (sense.first != lt))
					err_sense.push_back(obj_name);
				else
				{
					if (sense.first == gt)
					{
						keep_pi.push_back(obj_name);
						obj_sense_map[obj_name] = "maximize";
						obj_dir_mult[obj_name] = -1.0;
					}
					else if (sense.first == lt)
					{
						keep_pi.push_back(obj_name);
						obj_sense_map[obj_name] = "minimize";
						obj_dir_mult[obj_name] = 1.0;
					}
					pi_obj_names.push_back(obj_name);
				}
				if (prob_pareto)
                {
                    keep_pi_sd.push_back(obj_name + "_SD");
                }
			}
		}
		if (err_sense.size() > 0)
		{
			ss.str("");
			ss << "the following non-zero weighted 'mou_objectives' do not have the correct obs group naming convention (needed to identify objective direction):";
			for (auto e : err_sense)
				ss << e << ";";
			throw_moea_error(ss.str());
		}
		if ((keep_obs.size() == 0) && (keep_pi.size() == 0))
		{
			throw_moea_error("none of the supplied observation and/or prior info 'mou_objectives' were found in the non-zero-weighted observations");
		}
		
		else if (missing.size() > 0)
		{
			ss.str("");
			ss << "the following mou_objectives were not found in the non-zero-weighted observations or prior info eqs: ";
			for (auto m : missing)
				ss << m << ",";
			//message(1, ss.str());
            throw_moea_error(ss.str());

		}
		obj_names = passed_obj_names;
		obs_obj_names = keep_obs;
		obs_obj_sd_names = keep_obs_sd;
		pi_obj_names = keep_pi;
		pi_obj_sd_names = keep_pi_sd;
	}

	ss.str("");
	ss << "...using the following observations as objectives: " << endl;
	for (auto name : obs_obj_names)
	{
		ss << setw(30) << name << "   " << obj_sense_map[name] << endl;
	}
	file_manager.rec_ofstream() << ss.str();
	cout << ss.str();

	risk_obj = pest_scenario.get_pestpp_options().get_mou_risk_obj();
	
	if (pi_obj_names.size() > 0)
	{
		ss.str("");
		ss << "...using the following prior info eqs as objectives: " << endl;
		for (auto name : pi_obj_names)
		{
			ss << setw(30) << name << "   " << obj_sense_map[name] << endl;
		}
		file_manager.rec_ofstream() << ss.str();
		cout << ss.str();

	}
	
	if (risk_obj)
	{
		set<string> snames(act_par_names.begin(), act_par_names.end());
		if (snames.find(RISK_NAME) == snames.end())
			throw_moea_error("couldn't find '" + RISK_NAME + "' in adj par names for risk objective");
		if (find(dv_names.begin(), dv_names.end(), RISK_NAME) == dv_names.end())
		{
			//throw_moea_error(RISK_NAME + " not found in decision variable names");
			message(1, "adding '" + RISK_NAME + "' to decision variable names");
			dv_names.push_back(RISK_NAME);
		}
		if (find(pi_obj_names.begin(), pi_obj_names.end(), RISK_NAME) == pi_obj_names.end())
		{
			//throw_moea_error(RISK_NAME + " not found in prior information objective names");
			PriorInformation* pi_ptr = pest_scenario.get_prior_info_ptr();
			ParameterInfo par_info = pest_scenario.get_ctl_parameter_info();
			ss.str("");
			ss << RISK_NAME << " 1.0 * " << RISK_NAME << " = 0.5 1.0 greater_than";
			pi_ptr->AddRecord(pest_utils::upper_cp(ss.str()));
			message(1, "added prior information objective for '" + RISK_NAME + "': ", ss.str());
			pi_obj_names.push_back(RISK_NAME);
			obj_dir_mult[RISK_NAME] = -1;
		}
		ss.str("");
		ss << "'mou_risk_objective' is true, using " << RISK_NAME << " decision variable as risk in chance calcs";
		message(1, ss.str());

		//reset bounds of the risk parameter
		double b = pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(RISK_NAME)->lbnd;
		pest_scenario.get_ctl_parameter_info_ptr_4_mod()->get_parameter_rec_ptr_4_mod(RISK_NAME)->lbnd = max(b,0.01);
		b = pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(RISK_NAME)->ubnd;
		pest_scenario.get_ctl_parameter_info_ptr_4_mod()->get_parameter_rec_ptr_4_mod(RISK_NAME)->ubnd = min(b, 0.99);
		//set this just to make sure everything gets initialized right
		pest_scenario.get_pestpp_options_ptr()->set_opt_risk(0.95);


	}
	n_adaptive_dvs = 0;
	set<string> s_dv_names(dv_names.begin(), dv_names.end());
	if (s_dv_names.find(DE_F_NAME) != s_dv_names.end())
	{
		message(1, "self-adaptive diffevol 'F' value found, resetting bounds to min:0.5 max:1.0");
		double b = pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(DE_F_NAME)->lbnd;
		pest_scenario.get_ctl_parameter_info_ptr_4_mod()->get_parameter_rec_ptr_4_mod(DE_F_NAME)->lbnd = max(b, 0.5);
		b = pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(DE_F_NAME)->ubnd;
		pest_scenario.get_ctl_parameter_info_ptr_4_mod()->get_parameter_rec_ptr_4_mod(DE_F_NAME)->ubnd = min(b, 1.0);
		n_adaptive_dvs++;
	}

	if (s_dv_names.find(CR_NAME) != s_dv_names.end())
	{
		message(1, "self-adaptive crossover probability value found, resetting bounds to min:0.8, max:1.0");
		double b = pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(CR_NAME)->lbnd;
		pest_scenario.get_ctl_parameter_info_ptr_4_mod()->get_parameter_rec_ptr_4_mod(CR_NAME)->lbnd = max(b, 0.8);
		b = pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(CR_NAME)->ubnd;
		pest_scenario.get_ctl_parameter_info_ptr_4_mod()->get_parameter_rec_ptr_4_mod(CR_NAME)->ubnd = min(b, 1.0);
		n_adaptive_dvs++;
	}

	if (s_dv_names.find(MR_NAME) != s_dv_names.end())
	{
		message(1, "self-adaptive mutation probability value found, resetting bounds to min:0.01, max:0.3");
		double b = pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(MR_NAME)->lbnd;
		pest_scenario.get_ctl_parameter_info_ptr_4_mod()->get_parameter_rec_ptr_4_mod(MR_NAME)->lbnd = max(b, 0.8);
		b = pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(MR_NAME)->ubnd;
		pest_scenario.get_ctl_parameter_info_ptr_4_mod()->get_parameter_rec_ptr_4_mod(MR_NAME)->ubnd = min(b, 1.0);
		n_adaptive_dvs++;
	}
	constraints.initialize(dv_names, numeric_limits<double>::max());

	constraints.initial_report();

	if (pest_scenario.get_control_info().noptmax == 0)
	{
		message(0, "'noptmax'=0, running control file parameter values and quitting");

		Parameters pars = pest_scenario.get_ctl_parameters();
		ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();
		pts.ctl2numeric_ip(pars);
		ParameterEnsemble _pe(&pest_scenario, &rand_gen);
		_pe.reserve(vector<string>(), pest_scenario.get_ctl_ordered_adj_par_names());
		_pe.set_trans_status(ParameterEnsemble::transStatus::NUM);
		_pe.append(BASE_REAL_NAME, pars);
		string par_csv = file_manager.get_base_filename() + ".par.csv";
		//message(1, "saving parameter values to ", par_csv);
		//_pe.to_csv(par_csv);
		ParameterEnsemble pe_base = _pe;
		pe_base.reorder(vector<string>(), act_par_names);
		ObservationEnsemble _oe(&pest_scenario, &rand_gen);
		_oe.reserve(vector<string>(), pest_scenario.get_ctl_ordered_obs_names());
		_oe.append(BASE_REAL_NAME, pest_scenario.get_ctl_observations());
		ObservationEnsemble oe_base = _oe;
		oe_base.reorder(vector<string>(), act_obs_names);
		//initialize the phi handler
		Covariance parcov;
		parcov.from_parameter_bounds(pest_scenario, file_manager.rec_ofstream());
		L2PhiHandler ph(&pest_scenario, &file_manager, &oe_base, &pe_base, &parcov);
		if (ph.get_lt_obs_names().size() > 0)
		{
			message(1, "less_than inequality defined for observations: ", ph.get_lt_obs_names().size());
		}
		if (ph.get_gt_obs_names().size())
		{
			message(1, "greater_than inequality defined for observations: ", ph.get_gt_obs_names().size());
		}
		message(1, "running control file parameter values");

		vector<int> failed_idxs = run_population(_pe, _oe, false);
		if (failed_idxs.size() != 0)
		{
			message(0, "control file parameter value run failed...bummer");
			throw_moea_error(string("control file parameter value run failed"));
		}
		string obs_csv = file_manager.get_base_filename() + ".obs.csv";
		message(1, "saving results from control file parameter value run to ", obs_csv);
		_oe.to_csv(obs_csv);

		ph.update(_oe, _pe);
		message(0, "control file parameter phi report:");
		ph.report(true);
		ph.write(0, 1);
		save_real_par_rei(pest_scenario, _pe, _oe, output_file_writer, file_manager, -1, BASE_REAL_NAME);

		message(0, "control file parameter objective function summary: ");
		obj_func_report(_pe, _oe);

		vector<string> names = _oe.get_var_names();
		Observations obs(names, _oe.get_real_vector(BASE_REAL_NAME));
		names = _pe.get_var_names();
		pars.update(names, eigenvec_2_stlvec(_pe.get_real_vector(BASE_REAL_NAME)));
		
		constraints.mou_report(0, pars, obs, obs_obj_names, pi_obj_names, false);
		return;
	}
	
	string mou_generator = pest_utils::upper_cp(pest_scenario.get_pestpp_options().get_mou_generator());
	vector<string> tokens;
	int pop_size = pest_scenario.get_pestpp_options().get_mou_population_size();
	pest_utils::tokenize(mou_generator, tokens, ",");
	for (auto token : tokens)
	{
		if (token == "DE")
		{
			gen_types.push_back(MouGenType::DE);
			message(1, "using differential evolution generator");
		}
		else if (token == "SBX")
		{
			gen_types.push_back(MouGenType::SBX);
			message(1, "using simulated binary cross over generator");
		}
		else if (token == "PM")
		{
			gen_types.push_back(MouGenType::PM);
			message(1, "using polynomial mutation generator");
		}
		else if (token == "PSO")
		{
			gen_types.push_back(MouGenType::PSO);
			message(1, "using particle swarm generator");
			inertia_info = pest_scenario.get_pestpp_options().get_mou_pso_inertia();
			curr_omega = inertia_info[0];
			pso_dv_bound_handling = pest_scenario.get_pestpp_options().get_mou_pso_dv_bound_handling();
		}
        else if (token == "SIMPLEX")
        {
            gen_types.push_back(MouGenType::SMP);
            message(1, "using simplex generator");
        }
        else
		{
			throw_moea_error("unrecognized generator type '" + token + "', should be in {'DE','SBX','PM','PSO'}");
		}
	}
	//TODO: report constraints being applied

	sanity_checks();

	file_manager.open_ofile_ext(lineage_tag);
	ofstream& lin = file_manager.get_ofstream(lineage_tag);
	lin << "child,parent_1,parent_2,parent_3" << endl;

	//initialize the constraints using ctl file pars and obs
	//throughout the process, we can update these pars and obs
	//to control where in dec var space the stack/fosm estimates are
	//calculated
	//effective_constraint_pars = pest_scenario.get_ctl_parameters();
	//effective_constraint_obs = pest_scenario.get_ctl_observations();
	

	//int num_members = pest_scenario.get_pestpp_options().get_mou_population_size();
	population_dv_file = ppo->get_mou_dv_population_file();
	population_obs_restart_file = ppo->get_mou_obs_population_restart_file();
	
	initialize_dv_population();
	
	initialize_obs_restart_population();
	
	try
	{
		dp.check_for_dups();
	}
	catch (const exception& e)
	{
		string message = e.what();
		throw_moea_error("error in dv population: " + message);
	}

	try
	{
		op.check_for_dups();
	}
	catch (const exception& e)
	{
		string message = e.what();
		throw_moea_error("error in obs population: " + message);
	}

	if (pest_scenario.get_control_info().noptmax == -2)
	{
		message(0, "'noptmax'=-2, running mean parameter/decision variable ensemble values and quitting");
		message(1, "calculating mean parameter values");
		
		Parameters pars;
		vector<double> mv = dp.get_mean_stl_var_vector();
		pars.update(dp.get_var_names(), dp.get_mean_stl_var_vector());
		ParamTransformSeq pts = dp.get_par_transform();

		ParameterEnsemble _pe(&pest_scenario, &rand_gen);
		_pe.reserve(vector<string>(), dp.get_var_names());
		_pe.set_trans_status(dp.get_trans_status());
		_pe.append("mean", pars);
		string par_csv = file_manager.get_base_filename() + ".mean.par.csv";
		message(1, "saving mean parameter/decision variable values to ", par_csv);
		_pe.to_csv(par_csv);
		ParameterEnsemble pe_base = _pe;
		pe_base.reorder(vector<string>(), act_par_names);
		ObservationEnsemble _oe(&pest_scenario, &rand_gen);
		_oe.reserve(vector<string>(), op.get_var_names());
		_oe.append("mean", pest_scenario.get_ctl_observations());
		ObservationEnsemble oe_base = _oe;
		oe_base.reorder(vector<string>(), act_obs_names);
		//initialize the phi handler
		Covariance parcov;
		parcov.from_parameter_bounds(pest_scenario, file_manager.rec_ofstream());

		L2PhiHandler ph(&pest_scenario, &file_manager, &oe_base, &pe_base, &parcov);
		if (ph.get_lt_obs_names().size() > 0)
		{
			message(1, "less_than inequality defined for observations: ", ph.get_lt_obs_names().size());
		}
		if (ph.get_gt_obs_names().size())
		{
			message(1, "greater_than inequality defined for observations: ", ph.get_gt_obs_names().size());
		}
		message(1, "running mean values");

		vector<int> failed_idxs = run_population(_pe, _oe,false);
		if (failed_idxs.size() != 0)
		{
			message(0, "mean value run failed...bummer");
			return;
		}
		string obs_csv = file_manager.get_base_filename() + ".mean.obs.csv";
		message(1, "saving results from mean value run to ", obs_csv);
		_oe.to_csv(obs_csv);

		ph.update(_oe, _pe);
		message(0, "mean phi report:");
		ph.report(true);
		ph.write(0, 1);
		save_real_par_rei(pest_scenario, _pe, _oe, output_file_writer, file_manager, -1, "mean");
		message(0, "mean objective function summary: ");
		obj_func_report(_pe, _oe);

		vector<string> names = _oe.get_var_names();
		Observations obs(names, _oe.get_real_vector("mean"));
		names = _pe.get_var_names();
		pars.update(names, eigenvec_2_stlvec(_pe.get_real_vector("mean")));

		constraints.mou_report(0, pars, obs, obs_obj_names, pi_obj_names, false);
		return;
	}


	

	//we are restarting
	if (population_obs_restart_file.size() > 0)
	{
		if (constraints.get_use_chance())
		{
			//this can be done, but we need to make sure the appropriate chance restart
			//args were supplied: base_jacobian or obs_stack
			throw_moea_error("chance constraints/objectives not yet supported with restart");
		}
	
		//since mou reqs strict linking of realization names, let's see if we can find an intersection set 
		vector<string> temp = dp.get_real_names();
		set<string> dvnames(temp.begin(), temp.end());
		temp = op.get_real_names();
		set<string> obsnames(temp.begin(), temp.end());
		set<string> common;
		set_intersection(dvnames.begin(), dvnames.end(), obsnames.begin(), obsnames.end(),std::inserter(common,common.end()));
		
		// all members are common to both dp and op
		if (common.size() == dp.shape().first)
		{
			op.reorder(dp.get_real_names(), vector<string>());
		}
		
		//otherwise some members are not common
		else
		{
			ss.str("");
			ss << "WARNING: only " << common.size() << " members are common between the dv population and obs restart population.";
			message(1, ss.str());
			if (common.size() < error_min_members)
			{
			    if (pest_scenario.get_control_info().noptmax > 0)
				    throw_moea_error("too few members to continue");
			    else
                {
                    ss.str("");
                    ss << "WARNING: very few population members..." << endl;
                    message(0,ss.str());
                }

			}
			
			message(2,"aligning dv and obs populations");
			temp.clear();
			temp.resize(common.size());
			copy(common.begin(), common.end(), temp.begin());
			sort(temp.begin(), temp.end());
			dp.reorder(temp, vector<string>(), true);
			op.reorder(temp, vector<string>(), true);
			message(2, "dv population size: ", dp.shape().first);
			message(2, "obs population size", op.shape().first);
			message(2, "checking for denormal values in dv population");
			dp.check_for_normal("initial transformed dv population");
			message(2, "checking for denormal values in obs restart population");
			op.check_for_normal("restart obs population");
		}
		//TODO: make any risk runs that need to be done here or do we assume the 
		//restart population has already been shifted?

		//TODO: save both sim and sim+chance observation populations
		
	}
	else
	{
		message(2, "checking for denormal values in dv population");
		dp.check_for_normal("initial transformed dv population");
		//save the initial population once here
		ss.str("");
		ss << file_manager.get_base_filename() << ".0." << dv_pop_file_tag;
		if (pest_scenario.get_pestpp_options().get_save_binary())
		{
            if (pest_scenario.get_pestpp_options().get_save_dense())
            {
                ss << ".bin";
                dp.to_dense(ss.str());
            }
            else {
                ss << ".jcb";
                dp.to_binary(ss.str());
            }

		}
		else
		{
			ss << ".csv";
			dp.to_csv(ss.str());
		}
		
		message(1, " saved initial dv population to ", ss.str());
		performance_log->log_event("running initial population");
		message(1, "running initial population of size", dp.shape().first);
	
		vector<int> failed = run_population(dp, op, true);
		if (dp.shape().first == 0)
			throw_moea_error(string("all members failed during initial population evaluation"));
		
		dp.transform_ip(ParameterEnsemble::transStatus::NUM);
	}
	ss.str("");
	ss << file_manager.get_base_filename() << ".0." << obs_pop_file_tag;
	if (pest_scenario.get_pestpp_options().get_save_binary())
	{
        if (pest_scenario.get_pestpp_options().get_save_dense())
        {
            ss << ".bin";
            op.to_dense(ss.str());
        }
        else {
            ss << ".jcb";
            op.to_binary(ss.str());
        }
	}
	else
	{
		ss << ".csv";
		op.to_csv(ss.str());
	}
	message(1, " saved observation population to ", ss.str());

    message(0, "initial population decision variable summary:");
    previous_dv_summary = decvar_report(dp);

	message(0, "initial population objective function summary:");
	previous_obj_summary = obj_func_report(dp, op);

	if (constraints.get_use_chance())
	{
        string sum = constraints.mou_population_observation_constraint_summary(0,op,"pre-shift",obs_obj_names);
        frec << sum << endl;
        cout << sum << endl;
	    string opt_member;
		ObservationEnsemble shifted_op = get_chance_shifted_op(dp, op, opt_member);
		ss.str("");
		ss << file_manager.get_base_filename() << ".0.chance." << obs_pop_file_tag;
		if (pest_scenario.get_pestpp_options().get_save_binary())
		{
            if (pest_scenario.get_pestpp_options().get_save_dense())
            {
                ss << ".bin";
                shifted_op.to_dense(ss.str());

            }
            else {
                ss << ".jcb";
                shifted_op.to_binary(ss.str());
            }
		}
		else
		{
			ss << ".csv";
			shifted_op.to_csv(ss.str());
		}
		message(1, " saved chance-shifted observation population to ", ss.str());

		op = shifted_op;
		message(0, "chance-shifted initial population objective function summary");
		previous_obj_summary = obj_func_report(dp, op);
	}

	//save the initial dv population again in case runs failed or members were dropped as part of restart
	ss.str("");
	ss << file_manager.get_base_filename() << ".0." << dv_pop_file_tag;
	if (pest_scenario.get_pestpp_options().get_save_binary())
	{
        if (pest_scenario.get_pestpp_options().get_save_dense())
        {
            ss << ".bin";
            dp.to_dense(ss.str());
        }
        else {
            ss << ".jcb";
            dp.to_binary(ss.str());
        }
	}
	else
	{
		ss << ".csv";
		dp.to_csv(ss.str());
	}
	message(1, " saved initial dv population to ", ss.str());

	//TODO: think about a bad phi (or phis) for MOEA

	if (op.shape().first < error_min_members)
	{
		message(0, "too few population members:", op.shape().first);
		message(1, "need at least ", error_min_members);
		if (pest_scenario.get_control_info().noptmax > 0)
		    throw_moea_error(string("too few active population members, cannot continue"));
	}
	if (op.shape().first < warn_min_members)
	{
		ss.str("");
		ss << "WARNING: less than " << warn_min_members << " active population members...might not be enough";
		string s = ss.str();
		message(0, s);
	}


	//do an initial pareto dominance sort
	message(1, "performing initial pareto dominance sort");
	objectives.set_pointers(obj_names, obs_obj_names, obs_obj_sd_names, pi_obj_names, pi_obj_sd_names, obj_dir_mult,ppd_obj_to_sd);
	archive_size = ppo->get_mou_max_archive_size();
	vector<string> keep;
	if (envtype == MouEnvType::NSGA)
	{

		DomPair dompair = objectives.get_nsga2_pareto_dominance(iter, op, dp, &constraints, false, true, POP_SUM_TAG);

		//drop any duplicates
		keep.clear();
		for (auto nondom : dompair.first)
		{
			keep.push_back(nondom);
		}
		for (auto nondom : dompair.second)
		{
			keep.push_back(nondom);
		}
		if (keep.size() == 0)
		{
			throw_moea_error("initial sorting yielded zero valid solutions");
		}
		dp.keep_rows(keep);
		op.keep_rows(keep);

		//initialize op and dp archives
		op_archive = ObservationEnsemble(&pest_scenario, &rand_gen,
			op.get_eigen(dompair.first, vector<string>()), dompair.first, op.get_var_names());
		dp_archive = ParameterEnsemble(&pest_scenario, &rand_gen,
			dp.get_eigen(dompair.first, vector<string>()), dompair.first, dp.get_var_names());
		dp_archive.set_trans_status(dp.get_trans_status());
		ss << "initialized archives with " << dompair.first.size() << " nondominated members";
		message(2, ss.str());


		//this causes the initial archive pareto summary file to be written
		objectives.get_nsga2_pareto_dominance(iter, op_archive, dp_archive, &constraints, false, true, ARC_SUM_TAG);
		objectives.get_nsga2_pareto_dominance(iter, op_archive, dp_archive, &constraints, false, true, ARC_TRIM_SUM_TAG);

		//set hypervolume partitions of nondom solutions from previous outer iteration
		if (prob_pareto)
		{
			map<string, map<string, double>> hv_pts;
			vector<string> tokens;

			string outer_repo_obs_filename = pest_scenario.get_pestpp_options().get_mou_outer_repo_obs_file();
			if (outer_repo_obs_filename != "")
			{
				message(1, "loading outer repository obs from csv file", outer_repo_obs_filename);
				try
				{
					ifstream csv(outer_repo_obs_filename);
					string line;
					getline(csv, line);

					while (getline(csv, line))
					{
						pest_utils::strip_ip(line);
						tokens.clear();
						pest_utils::tokenize(line, tokens, ",", false);
						map<string, double> vals;
						int i = 1;
						for (auto obj : obs_obj_names)
						{
							vals[obj] = stod(tokens[i]);
							i++;
						}
						hv_pts[tokens[0]] = vals;
					}

				}
				catch (const exception& e)
				{
					ss << "error processing outer repository obs file: " << e.what();
					throw_moea_error(ss.str());
				}
				catch (...)
				{
					throw_moea_error(string("error processing outer repository obs file"));
				}
			}
			else
			{
				message(1, "using the initial population for hypervolume partitioning");
				stringstream ss;
				ofstream& frec = file_manager.rec_ofstream();
				ss << "ParetoObjectives::get_hypervolume() for " << op.shape().first << " archive members";
				performance_log->log_event(ss.str());

				hv_pts = objectives.get_members(op_archive, dp_archive);
			}
			objectives.set_hypervolume_partitions(hv_pts);
			//objectives.get_ehvi(op, dp);
		}
			
				
	}
	else if (envtype == MouEnvType::SPEA)
	{
		map<string, double> fit = objectives.get_spea2_fitness(iter, op, dp, &constraints, true, POP_SUM_TAG);
		keep.clear();
		for (auto& rname : dp.get_real_names())
        {
		    if (fit.find(rname) != fit.end())
		        keep.push_back(rname);
        }
		dp.keep_rows(keep);
		op.keep_rows(keep);
		op_archive = op; //copy
		dp_archive = dp; //copy
		vector<string> keep;
		for (auto& f : fit)
			if (f.second < 1.0)
				keep.push_back(f.first);
		op_archive.keep_rows(keep);
		dp_archive.keep_rows(keep);

		if (keep.size() == 0)
		{
			message(2, "initial archive empty - no feasible non-dominated solutions...");
		}
		else
		{
			ss.str("");
			ss << "initialized archives with " << keep.size() << " nondominated members";
			message(2, ss.str());
			//dont set this here b/c it is set dynamically each gen using population schedule
			//archive_size = num_members * 2;

			//this causes the initial archive pareto summary file to be written
			objectives.get_spea2_fitness(iter, op_archive, dp_archive, &constraints, true, ARC_SUM_TAG);
		}

	}
    save_populations(dp_archive, op_archive, "archive");


    string opt_member;
	if (constraints.get_use_chance())
	{
		ofstream& f_rec = file_manager.rec_ofstream();

		pair<Parameters, Observations> po_pair = get_optimal_solution(dp, op, opt_member);
		constraints.presolve_chance_report(iter, po_pair.second, true," --- initial chance constraint summary (calculated at optimal/mean decision variable point) --- ");
	}

	constraints.mou_report(0,dp, op, obs_obj_names,pi_obj_names);

    initialize_pso();

    par_sim_map.clear();
    obs_sim_map.clear();
    update_sim_maps(dp,op);
    ss.str("");
    ss << "number of initial feasible solutions: " << objectives.get_num_feasible();
    message(1,ss.str());
	message(0, "initialization complete");
}

void MOEA::update_sim_maps(ParameterEnsemble& _dp, ObservationEnsemble& _op)
{
    map<string,int> rmap = _dp.get_real_map();
    for (auto& ridx : rmap)
    {
        par_sim_map[ridx.first] = _dp.get_eigen_ptr()->row(ridx.second);
    }
    rmap = _op.get_real_map();
    for (auto& ridx : rmap)
    {
        obs_sim_map[ridx.first] = _op.get_eigen_ptr()->row(ridx.second);
    }

}

ParameterEnsemble MOEA::get_initial_pso_velocities(int num_members) {
	ParameterEnsemble _pso_velocity = dp.zeros_like(num_members);
	Parameters lb = pest_scenario.get_ctl_parameter_info().get_low_bnd(dv_names);
	Parameters ub = pest_scenario.get_ctl_parameter_info().get_up_bnd(dv_names);
	_pso_velocity.get_par_transform().ctl2numeric_ip(lb);
	_pso_velocity.get_par_transform().ctl2numeric_ip(ub);
	Parameters dist = ub - lb;

	pso_vmax.clear();
	double vmax_scale_factor = pest_scenario.get_pestpp_options().get_mou_pso_vmax_factor();
	for (auto& dv_name : dv_names) 
		pso_vmax[dv_name] = dist[dv_name] * vmax_scale_factor;
	

	double init_vel_scale_fac = 0.5;
	for (auto& dv_name : dv_names)
	{
		vector<double> vals = uniform_draws(num_members, -dist[dv_name] * init_vel_scale_fac, dist[dv_name] * init_vel_scale_fac, rand_gen);
		Eigen::VectorXd real = stlvec_2_eigenvec(vals);
		_pso_velocity.replace_col(dv_name, real);
	}
	
    return _pso_velocity;
}

void MOEA::initialize_pso()
{
	pso_velocity = get_initial_pso_velocities(dp.shape().first);
    update_pso_velocity_map(pso_velocity);
	pso_pbest_dp = ParameterEnsemble(&pest_scenario, &rand_gen, dp.get_eigen(), dp.get_real_names(), dp.get_var_names());
	pso_pbest_dp.set_trans_status(dp.get_trans_status());
	pso_pbest_op = ObservationEnsemble(&pest_scenario, &rand_gen, op.get_eigen(), op.get_real_names(), op.get_var_names());

}

void MOEA::update_pso_velocity_map(ParameterEnsemble& _pso_velocity)
{
    map<string,int> rmap = _pso_velocity.get_real_map();
    for (auto& ridx : rmap)
    {
        pso_velocity_map[ridx.first] = _pso_velocity.get_eigen_ptr()->row(ridx.second);
    }
}

pair<Parameters, Observations> MOEA::get_optimal_solution(ParameterEnsemble& _dp, ObservationEnsemble& _op, string& opt_member_name)
{
	Parameters pars;
	Observations obs;
	stringstream ss;
	bool use_mean = false;
	if (((iter == 0) && (population_obs_restart_file.size() == 0)) || (obs_obj_names.size() == 0))
		use_mean = true;

	if (use_mean)
	{
		//just use dp member nearest the mean dec var values

		vector<double> t = _dp.get_mean_stl_var_vector();
		Eigen::VectorXd dp_mean = stlvec_2_eigenvec(t);
		t.resize(0);
		int idx_min = -1;
		double dist, dist_min = numeric_limits<double>::max();
		for (int i = 0; i < _dp.shape().first; i++)
		{
			//dist = dp_mean.dot(_dp.get_eigen_ptr()->row(i));
            dist = (dp_mean - _dp.get_eigen_ptr()->row(i).transpose()).squaredNorm();
			if (dist < dist_min)
			{
				idx_min = i;
				dist_min = dist;
			}
		}
		if (idx_min == -1)
			throw_moea_error("couldn't find nearest mean point");
		string min_member = _dp.get_real_names()[idx_min];
		if (dist_min > 0.0) dist_min = sqrt(dist_min);
		ss.str("");
		ss << "using member " << min_member << " as nearest-to-mean single point" << endl;
		ss << "    with distance of " << dist_min << " from mean of decision variable population";
		message(2, ss.str());
		opt_member_name = min_member;

		pars.update_without_clear(_dp.get_var_names(), _dp.get_real_vector(min_member));
		obs.update_without_clear(_op.get_var_names(), _op.get_real_vector(min_member));
		pest_scenario.get_base_par_tran_seq().numeric2ctl_ip(pars);
	}
	else
	{
		//calculate the optimal tradeoff point from the current op
		//dont worry about pi-based obj since they aren't chance-based
		message(2, "seeking optimal trade-off point for single 'optimal' chance point runs");
		vector<double> obj_extrema;
		Eigen::VectorXd obj_vec; 

		for (auto obj_name : obs_obj_names)
		{
			//if this is a max obj
			if ((obj_dir_mult.find(obj_name) != obj_dir_mult.end()) &&
				(obj_dir_mult[obj_name] == -1.0))
				obj_extrema.push_back(_op.get_var_vector(obj_name).maxCoeff());
			else
				obj_extrema.push_back(_op.get_var_vector(obj_name).minCoeff());
		}

		Eigen::VectorXd opt_vec = stlvec_2_eigenvec(obj_extrema);

		//find the member nearest the optimal tradeoff
		int opt_idx = -1;
		Eigen::MatrixXd obj_op = _op.get_eigen(vector<string>(), obs_obj_names);
		double dist, opt_dist = numeric_limits<double>::max();
		for (int i = 0; i < _op.shape().first; i++)
		{
			dist = (opt_vec - obj_op.row(i).transpose()).squaredNorm();
			if (dist < opt_dist)
			{
				opt_idx = i;
				opt_dist = dist;
			}
		}
		if (opt_idx == -1)
			throw_moea_error("couldn't find nearest optimal point");
		string opt_member = _op.get_real_names()[opt_idx];
		if (opt_dist > 0.0) opt_dist = sqrt(opt_dist);
		ss.str("");
		ss << "using member " << opt_member << " as single, 'optimal' point" << endl;
		ss << "   with distance of " << opt_dist << " from optimal trade - off";
		message(2, ss.str());
		pars.update_without_clear(_dp.get_var_names(), _dp.get_real_vector(opt_member));
		obs.update_without_clear(_op.get_var_names(), _op.get_real_vector(opt_member));
		opt_member_name=opt_member;
        pest_scenario.get_base_par_tran_seq().numeric2ctl_ip(pars);
	}
	return pair<Parameters, Observations>(pars, obs);
}


ParameterEnsemble MOEA::generate_population()
{
	//int total_new_members = pest_scenario.get_pestpp_options().get_mou_population_size();
    int total_new_members = population_schedule.at(iter);
	//add new members for any missing
	//total_new_members += (total_new_members - dp.shape().first);
	int new_members_per_gen = int(total_new_members / gen_types.size());
	ParameterEnsemble new_pop(&pest_scenario, &rand_gen);
	new_pop.set_trans_status(ParameterEnsemble::transStatus::NUM);
	objectives.get_nsga2_pareto_dominance(iter, op, dp, &constraints, prob_pareto, false);
	for (auto gen_type : gen_types)
	{
		ParameterEnsemble p(&pest_scenario);
		if (gen_type == MouGenType::DE)
		{
			p = generate_diffevol_population(new_members_per_gen, dp);
		}

		else if (gen_type == MouGenType::SBX)
		{
			p = generate_sbx_population(new_members_per_gen, dp);
		}
		else if (gen_type == MouGenType::PM)
		{
			p = generate_pm_population(new_members_per_gen, dp);
		}
		else if (gen_type == MouGenType::PSO)
		{
			p = generate_pso_population(new_members_per_gen, dp);
		}
        else if (gen_type == MouGenType::SMP)
        {
            p = generate_simplex_population(new_members_per_gen, dp, op);
        }
        else
			throw_moea_error("unrecognized mou generator");
		if (new_pop.shape().first == 0)
			new_pop = p;
		else
			new_pop.append_other_rows(p);
	}

	if ((pest_scenario.get_pestpp_options().get_mou_shuffle_fixed_pars()) && (new_pop.get_fixed_info().get_map_size() > 0))
    {

	    vector<string> real_names = new_pop.get_real_names();
        vector<string> fixed_names = new_pop.get_fixed_info().get_fixed_names();

        vector<int> ireals;
	    for (int i=0;i<real_names.size();i++)
	        ireals.push_back((i));
        shuffle(ireals.begin(),ireals.end(),rand_gen);
        string name1,name2;
        Eigen::MatrixXd new_fixed(real_names.size(),fixed_names.size());
        new_fixed.setZero();
        vector<double> t;
        for (int i=0;i<real_names.size();i++)
        {
            name1 = real_names[i];
            name2 = real_names[ireals[i]];
            t = new_pop.get_fixed_info().get_real_fixed_values(name2,fixed_names);
            new_fixed.row(i) = Eigen::Map<Eigen::VectorXd>(&t[0],t.size());
            //cout << name1 << ", " << new_fixed.row(i) << endl;
        }
        new_pop.get_fixed_info().update_realizations(fixed_names,real_names,new_fixed);
    }
	

	return new_pop;
}

void MOEA::fill_populations_from_maps(ParameterEnsemble& new_dp, ObservationEnsemble& new_op )
{
    vector<string> rnames;
    for (auto& entry : par_sim_map)
    {
        rnames.push_back(entry.first);
    }
    new_dp.reserve(rnames,dp.get_var_names());

    map<string,int> rmap = new_dp.get_real_map();
    for (auto& ridx : rmap)
    {
        new_dp.get_eigen_ptr_4_mod()->row(ridx.second) = par_sim_map.at(ridx.first);
    }
    new_op.reserve(rnames,op.get_var_names());
    rmap = new_op.get_real_map();
    for (auto& ridx : rmap)
    {
        new_op.get_eigen_ptr_4_mod()->row(ridx.second) = obs_sim_map.at(ridx.first);
    }

}

void MOEA::iterate_to_solution()
{
	iter = 1;
	int num_members = pest_scenario.get_pestpp_options().get_mou_population_size();
	vector<string> keep;
	stringstream ss;
	map<string, map<string, double>> summary;
	int noptmax = pest_scenario.get_control_info().noptmax;
	while(iter <= noptmax)
	{
		message(0, "starting generation ", iter);

		if (dp.shape().first < error_min_members)
        {
            throw_moea_error("too few members to continue");
        }
		if (dp.shape().first < warn_min_members)
        {
		    message(0,"WARNING: very few members in current population...");
        }

		//generate offspring
		ParameterEnsemble new_dp = generate_population();
		
		//run offspring thru the model while also running risk runs, possibly at many points in dec var space	
		ObservationEnsemble new_op(&pest_scenario, &rand_gen);
		new_op.reserve(new_dp.get_real_names(), op.get_var_names());
		run_population(new_dp, new_op, true);

		save_populations(new_dp, new_op);
        //update_sim_maps(new_dp,new_op);


        if (pest_scenario.get_pestpp_options().get_mou_use_multigen())
        {
            update_sim_maps(new_dp,new_op);
        }

        if (should_use_multigen())
        {
            message(1,"using multi-generational population in dominance sorting");
            fill_populations_from_maps(new_dp,new_op);
        }

		if (constraints.get_use_chance())
		{

            if (pest_scenario.get_pestpp_options().get_mou_verbose_level() > 2) {
                ss.str("");
                ss << "all.pre-shift";
                save_populations(new_dp, new_op, ss.str());
            }


            string csum = constraints.mou_population_observation_constraint_summary(iter,new_op,"pre-shift",obs_obj_names);
		    cout << csum;
		    file_manager.rec_ofstream() << csum;
		    string opt_member;
			pair<Parameters,Observations> po = get_optimal_solution(dp, op, opt_member);
			constraints.presolve_chance_report(iter, po.second,true, "chance constraint summary (calculated at optimal/mean decision variable point)");
			ObservationEnsemble new_op_shifted = get_chance_shifted_op(new_dp, new_op,opt_member);
			save_populations(new_dp, new_op_shifted,"chance");
			new_op = new_op_shifted;
		}
        else if (!should_use_multigen()) {
            //append offspring dp and (risk-shifted) op to make new dp and op containers

            new_dp.append_other_rows(dp);
            if (dp.get_fixed_info().get_map_size() > 0)
            {
                map<string,map<string,double>> fi = dp.get_fixed_info().get_fixed_info_map();
                new_dp.get_fixed_info().add_realizations(fi);
            }
            new_op.append_other_rows(op);
        }

        if (find(gen_types.begin(),gen_types.end(),MouGenType::PSO) != gen_types.end()) {
            update_pso_pbest(new_dp, new_op);
        }

		if (envtype == MouEnvType::NSGA)
		{
			message(1, "pareto dominance sorting combined parent-child populations of size ", new_dp.shape().first);
			DomPair dompair = objectives.get_nsga2_pareto_dominance(iter, new_op, new_dp, &constraints, prob_pareto, true, POP_SUM_TAG);
			//the ordering must not be by fitness

            if (should_use_multigen()) {
                message(2,"keeping all feasible nondom members from multi-generational population");
                keep = dompair.first;
            }
            else {
                keep.clear();
                for (auto nondom : dompair.first) {
                    if (keep.size() >= num_members)
                        break;
                    keep.push_back(nondom);
                }
            }

			if (keep.size() > 0)
			{
				//update the archive of nondom members
				ParameterEnsemble new_dp_nondom = new_dp;
				new_dp_nondom.keep_rows(keep);
				ObservationEnsemble new_op_nondom = new_op;
				new_op_nondom.keep_rows(keep);
				update_archive_nsga(new_op_nondom, new_dp_nondom);
			}

			//now fill out the rest of keep with dom solutions
			for (auto dom : dompair.second)
			{
				if (keep.size() >= num_members)
					break;
				keep.push_back(dom);
			}

			message(1, "resizing current populations to ", keep.size());
			new_dp.keep_rows(keep);
			new_op.keep_rows(keep);
			dp = new_dp;
			op = new_op;



		}

		else if (envtype == MouEnvType::SPEA)
		{
			map<string, double> fit = objectives.get_spea2_fitness(iter, new_op, new_dp, &constraints, true, POP_SUM_TAG);
			//first find all members with fitness less than 1 (nondom)
            if (should_use_multigen())
            {
                message(2,"keeping all feasible nondom members from multigenerational population");
                keep = new_dp.get_real_names();
            }
            else {
                keep.clear();
                for (auto member : new_dp.get_real_names()) {
                    if (fit[member] < 1.0)
                        keep.push_back(member);
                }
                ss.str("");
                ss << keep.size() << " non-dominated members (spea2 fitness less than 1.0)";
                cout << ss.str() << endl;
                file_manager.rec_ofstream() << ss.str() << endl;
                if (keep.size() < num_members) {
                    //fill with members of increasing fitness value
                    sortedset fit_sorted(fit.begin(), fit.end(), compFunctor);
                    sortedset::iterator it = next(fit_sorted.begin(), keep.size());
                    for (; it != fit_sorted.end(); ++it) {
                        keep.push_back(it->first);
                        if (keep.size() == num_members)
                            break;
                    }
                    //cout << keep.size() << endl;;
                }
            }
			if (keep.size() > num_members)
			{
				objectives.get_spea2_archive_names_to_keep(num_members, keep, new_op, new_dp);
			}
			message(1, "resizing current populations to ", keep.size());
			new_dp.keep_rows(keep);
			new_op.keep_rows(keep);
			dp = new_dp;
			op = new_op;
			update_archive_spea(op, dp);
		}

		else
		{
			throw_moea_error("unrecognized 'mou_env'");
		}
        ss.str("");
        ss << "generation " << iter << " decision variable summary:";
        message(0, ss.str());
        summary = decvar_report(dp);

        ss.str("");
        ss << "generation " << iter << " decision variable change summary:";
        message(0, ss.str());
        decvar_change_report(summary);
        previous_dv_summary = summary;

		ss.str("");
		ss << "generation " << iter << " objective function summary:";
		message(0, ss.str());
		summary = obj_func_report(dp, op);
		ss.str("");
		ss << "generation " << iter << " objective function change summary:";
		message(0, ss.str());
		obj_func_change_report(summary);
		previous_obj_summary = summary;
		constraints.mou_report(iter, new_dp, new_op, obs_obj_names, pi_obj_names, true);
        ss.str("");
        ss << "number of feasible solutions at the end of generation " << iter << ": " << objectives.get_num_feasible();
        message(1,ss.str());
		iter++;
        int q = pest_utils::quit_file_found();
        if ((q == 1) || (q == 2))
        {
		    message(0,"'pest.stp' found, quitting");
            message(1,"force-saving current populations");
            save_populations(dp,op,"",true);
            save_populations(dp_archive, op_archive, "archive",true);
		    break;
        }
        else if (q == 4) {
            message(0,"pest.stp found with '4'.  run mgr has returned control, removing file.");

            if (!pest_utils::try_remove_quit_file()) {
                message(0,"error removing pest.stp file, bad times ahead...");
            }
        }
	}

}

bool MOEA::should_use_multigen() {

    if (pest_scenario.get_pestpp_options().get_mou_use_multigen()) {
        return true;
    }
    //if ((constraints.should_update_chance(iter)) && (pest_scenario.get_pestpp_options().get_opt_recalc_fosm_every() != 1))
    //{
    //    return true;
    //}
    return false;
}

void MOEA::initialize_population_schedule()
{
    stringstream ss;
    population_schedule.clear();
    int num_members = pest_scenario.get_pestpp_options().get_mou_population_size();
    string fname = pest_scenario.get_pestpp_options().get_mou_population_schedule();
    string line;
    vector<string> tokens;
    int lcount = 0, gen,psize;
    for (int i=0;i<max(1,pest_scenario.get_control_info().noptmax+1);i++)
        population_schedule[i] = num_members;
    if (fname.size() > 0)
    {
        message(2,"reading population schedule from file '"+fname+"'");
        ifstream in(fname);
        if (in.bad())
        {
            throw_moea_error("error opening mou_population_schedule file '"+fname+"'");
        }
        while (getline(in,line))
        {
            lcount++;
            tokens.clear();
            pest_utils::tokenize(line,tokens,"\t ,");
            if (tokens.size() < 2)
            {
                ss.str("");
                ss << "mou_population_schedule file '" << fname << "' line " << lcount << " needs at least two entries";
                throw_moea_error(ss.str());
            }
            try
            {
                gen = stoi(tokens[0]);
            }
            catch (...)
            {
                ss.str("");
                ss << "error casting '" << tokens[0] << "' to generation integer on line " << lcount << " in mou_schedule_file";
                throw_moea_error(ss.str());
            }
            try
            {
                psize = stoi(tokens[1]);
            }
            catch (...)
            {
                ss.str("");
                ss << "error casting '" << tokens[1] << "' to population size integer on line " << lcount << " in mou_schedule_file";
                throw_moea_error(ss.str());
            }
            if (psize < error_min_members) {
                ss.str("");
                ss << "population size " << psize << " on line " << lcount << " in mou_schedule_file less than min "
                   << error_min_members;
                throw_moea_error(ss.str());
            }
            population_schedule[gen] = psize;
        }
        in.close();
    }

    ofstream& frec = file_manager.rec_ofstream();
    frec << "...population schedule: generation,population size:" << endl;
    for (int i=0;i<pest_scenario.get_control_info().noptmax;i++)
    {
        frec << "...   " << i << ", " << population_schedule.at(i) << endl;
    }
}

bool MOEA::initialize_dv_population()
{
	stringstream ss;
	//int num_members = pest_scenario.get_pestpp_options().get_mou_population_size();
	int num_members = population_schedule[0];
	string dv_filename = pest_scenario.get_pestpp_options().get_mou_dv_population_file();
	bool drawn = false;
	if (dv_filename.size() == 0)
	{
		ofstream& frec = file_manager.rec_ofstream();
		message(1, "drawing initial dv population of size: ", num_members);
		Parameters draw_par = pest_scenario.get_ctl_parameters();
		
		dp.draw_uniform(num_members, dv_names, performance_log, 1, file_manager.rec_ofstream());
		vector<string> real_names;
		
		for (int i = 0; i < num_members; i++)
		{
			ss.str("");
			ss << "GEN=0_MEMBER=" << i;
			real_names.push_back(ss.str());
			member_count++;
		}
		dp.set_real_names(real_names, true);
	
		drawn = true;
	}
	else
	{
		string par_ext = pest_utils::lower_cp(dv_filename).substr(dv_filename.size() - 3, dv_filename.size());
		performance_log->log_event("processing dv population file " + dv_filename);
		if (par_ext.compare("csv") == 0)
		{
			message(1, "loading dv population from csv file", dv_filename);
			try
			{
				dp.from_csv(dv_filename);
			}
			catch (const exception& e)
			{
				ss << "error processing dv population file: " << e.what();
				throw_moea_error(ss.str());
			}
			catch (...)
			{
				throw_moea_error(string("error processing dv population file"));
			}
		}
		else if ((par_ext.compare("jcb") == 0) || (par_ext.compare("jco") == 0) || (par_ext.compare("bin") == 0))
		{
			message(1, "loading dv population from binary file", dv_filename);
			try
			{
				dp.from_binary(dv_filename);
			}
			catch (const exception& e)
			{
				ss << "error processing binary file: " << e.what();
				throw_moea_error(ss.str());
			}
			catch (...)
			{
				throw_moea_error(string("error processing binary file"));
			}
		}
		else
		{
			ss << "unrecognized dv population file extension " << par_ext << ", looking for csv, jcb, bin, or jco";
			throw_moea_error(ss.str());
		}

		dp.transform_ip(ParameterEnsemble::transStatus::NUM);
		ss.str("");
		ss << "dv population with " << dp.shape().first << " members read from '" << dv_filename << "'" << endl;
		message(1, ss.str());

		if (dp.shape().first == 0)
		{
			throw_moea_error("zero members found in dv population file");
		}

		if (pp_args.find("MOU_POPULATION_SIZE") != pp_args.end())
		{
			int num_members = pest_scenario.get_pestpp_options().get_mou_population_size();
			
			if (num_members < dp.shape().first)
			{
				message(1, "'mou_population_size' arg passed, truncated dv population to ", num_members);
				vector<string> keep_names, real_names = dp.get_real_names();
				for (int i = 0; i < num_members; i++)
				{
					keep_names.push_back(real_names[i]);
				}
				dp.keep_rows(keep_names);
			}
		}
		//this is just to make sure we dont accidentally generate duplicate real nanmes later...
		int gen,max_gen = 0;
		vector<string> tokens;
		string t;
		for (auto name : dp.get_real_names())
		{
			tokens.clear();
			if (name.find("GEN=") != string::npos)
			{
				pest_utils::tokenize(name, tokens, "=");
				t = tokens[1];
				tokens.clear();
				pest_utils::tokenize(t, tokens, "_");
				if (tokens.size() > 0)
				{
					t = tokens[0];
					try
					{
						gen = stoi(t);
					}
					catch (...)
					{
						continue;
					}
					max_gen = max(max_gen, gen);
				}
			}
		}
		if (max_gen > 0)
		{
			message(1, "previous generation numbers detected in dv population, fast forwarding generation counter to", max_gen);
			restart_iter_offset = max_gen;
		}
	}
	return drawn;


}


void MOEA::initialize_obs_restart_population()
{
	string obs_filename = pest_scenario.get_pestpp_options().get_mou_obs_population_restart_file();
	if (obs_filename.size() == 0)
	{
		op.reserve(dp.get_real_names(), pest_scenario.get_ctl_ordered_obs_names());
		return;
	}
	stringstream ss;
	//throw_moea_error(string("restart obs population not implemented"));
	string par_ext = pest_utils::lower_cp(obs_filename).substr(obs_filename.size() - 3, obs_filename.size());
	performance_log->log_event("processing obs population file " + obs_filename);
	if (par_ext.compare("csv") == 0)
	{
		message(1, "loading obs population from csv file", obs_filename);
		try
		{
			op.from_csv(obs_filename);
		}
		catch (const exception& e)
		{
		ss << "error processing obs population file: " << e.what();
		throw_moea_error(ss.str());
		}
		catch (...)
		{
			throw_moea_error(string("error processing obs population file"));
		}
	}
	else if ((par_ext.compare("jcb") == 0) || (par_ext.compare("jco") == 0) || (par_ext.compare("bin") == 0))
	{
	message(1, "loading obs population from binary file", obs_filename);
	try
	{
		op.from_binary(obs_filename);
	}
	catch (const exception& e)
	{
		ss << "error processing obs population binary file: " << e.what();
		throw_moea_error(ss.str());
	}
	catch (...)
	{
		throw_moea_error(string("error processing obs population binary file"));
	}
	}
	else
	{
	ss << "unrecognized obs population restart file extension " << par_ext << ", looking for csv, jcb, bin,  or jco";
	throw_moea_error(ss.str());
	}

	ss.str("");
	ss << "obs population with " << op.shape().first << " members read from '" << obs_filename << "'" << endl;
	message(1, ss.str());
	if (op.shape().first == 0)
	{
		throw_moea_error("zero members found in obs population restart file");
	}

}



void MOEA::update_pso_pbest(ParameterEnsemble& _dp, ObservationEnsemble& _op)
{
	ParameterEnsemble tdp = _dp;
	ObservationEnsemble top = _op;
	objectives.update(top, tdp, &constraints);
	objectives.get_nsga2_pareto_dominance(-999, _op, _dp, &constraints, prob_pareto, false);
	Eigen::VectorXd real;
	string f, s;
	vector<string> names = _dp.get_real_names();
	set<string> snames(names.begin(), names.end());
	names.clear();
	bool new_dom_old;
	//vector<string> new_pbest_names;
	//pso_pbest_dp = _dp;
	//pso_pbest_op = _op;
	set<string> duplicates = objectives.get_duplicates();

	for (auto lm : current_pso_lineage_map)
	{
		
		f = lm.second; s = lm.first;
		if ((snames.find(f) == snames.end()) || (snames.find(s) == snames.end()))
		{
			names.push_back(f);
		}
		else if ((duplicates.find(f) != duplicates.end()) || (duplicates.find(s) != duplicates.end()))
		{
			names.push_back(f);
		}
		else
		{
			new_dom_old = objectives.compare_two(f, s, envtype);
			if (!new_dom_old)
			{
				real = pso_pbest_dp.get_real_vector(s);
				tdp.update_real_ip(f, real);
				real = pso_pbest_op.get_real_vector(s);
				top.update_real_ip(f, real);
			}
			//new_pbest_names.push_back(lm.first);

		}
	}
	pso_pbest_dp = tdp;
	pso_pbest_op = top;
}

pair<ParameterEnsemble, ParameterEnsemble> MOEA::get_updated_pso_velocity(ParameterEnsemble& _dp, vector<string>& gbest_solutions)
{
	double cog_const, social_const;
	stringstream ss;
	vector<double> cog_const_range = pest_scenario.get_pestpp_options().get_mou_pso_cognitive_const();
	if (cog_const_range.size() == 1)
		cog_const = cog_const_range[0];
	else if (cog_const_range.size() == 2)
	{
		cog_const = cog_const_range[0] + (cog_const_range[1] - cog_const_range[0]) * (iter / pest_scenario.get_control_info().noptmax);
		message(1, "computing pso velocity using cognitive const: ", cog_const);
	}
	else
		throw_moea_error("invalid cognitive const range");

	vector<double> social_const_range = pest_scenario.get_pestpp_options().get_mou_pso_social_const();
	if (social_const_range.size() == 1)
		social_const = social_const_range[0];
	else if (social_const_range.size() == 2)
	{
		social_const = social_const_range[0] + (social_const_range[1] - social_const_range[0]) * (iter / pest_scenario.get_control_info().noptmax);
		message(1, "computing pso velocity using social const: ", social_const);
	}
	else
		throw_moea_error("invalid social const range");

	int num_dv = _dp.shape().second;
	vector<double> r;
	Eigen::VectorXd rand1, rand2, cur_real, p_best, g_best, new_par_vel, new_par_dval, cur_vel, inertia_comp, social_comp, cog_comp;
	pso_pbest_dp.transform_ip(_dp.get_trans_status());
	dp_archive.set_trans_status(_dp.get_trans_status());
	Eigen::MatrixXd new_vel(_dp.shape().first, _dp.shape().second), child_dv(_dp.shape().first, _dp.shape().second);
	string real_name;
	vector<string> real_names = pso_velocity.get_real_names();
	set<string> snames(real_names.begin(), real_names.end());

	double omega;
	if (((iter - 1) <= inertia_info[2]) && (inertia_info[2] != 0))
	{
		omega = inertia_info[0] + (inertia_info[1] - inertia_info[0]) * ((iter - 1) / inertia_info[2]);
		curr_omega = omega;
		message(1, "computing pso velocity using inertia weight: ", omega);
	}
	else
		omega = curr_omega;

	Parameters lb = pest_scenario.get_ctl_parameter_info().get_low_bnd(dv_names);
	Parameters ub = pest_scenario.get_ctl_parameter_info().get_up_bnd(dv_names);

	real_names = _dp.get_real_names();
	for (int i = 0; i < _dp.shape().first; i++)
	{
		real_name = real_names[i];
		if (snames.find(real_name) != snames.end())
			cur_vel = pso_velocity.get_real_vector(real_name);
		else
		{
			//cur_vel = pso_velocity.get_real_vector(current_pso_lineage_map.at(real_name));
			cur_vel = pso_velocity_map.at(real_name);
		}

		r = uniform_draws(num_dv, 0.0, 1.0, rand_gen);
		rand1 = stlvec_2_eigenvec(r);
		r = uniform_draws(num_dv, 0.0, 1.0, rand_gen);
		rand2 = stlvec_2_eigenvec(r);
		cur_real = _dp.get_real_vector(real_name);
		p_best = pso_pbest_dp.get_real_vector(real_name);
		g_best = dp_archive.get_real_vector(gbest_solutions[i]);

		inertia_comp = omega * cur_vel.array();
		cog_comp = cog_const * rand1.array() * (p_best.array() - cur_real.array());
		social_comp = social_const * rand2.array() * (g_best.array() - cur_real.array());

		new_par_vel = inertia_comp + cog_comp + social_comp;
		new_par_dval = cur_real.array() + new_par_vel.array();

		if (pso_dv_bound_handling == "CLAMP") //this replaces the original "reset/basic"
		{
			for (int j = 0; j < dv_names.size(); j++)
			{

				double lb_val = lb[dv_names[j]];
				double ub_val = ub[dv_names[j]];
				double new_dv, cur_dv;
				new_dv = new_par_dval[j] < lb_val - FLOAT_EPSILON ? lb_val : new_par_dval[j];
				new_dv = new_par_dval[j] > ub_val + FLOAT_EPSILON ? ub_val : new_par_dval[j];
				new_par_dval[j] = new_dv;
			}
			
		}
		else
		{
			double new_dv, cur_dv;
			for (int j = 0; j < dv_names.size(); j++)
			{
				double lb_val = lb[dv_names[j]];
				double ub_val = ub[dv_names[j]];

				double vmax = pso_vmax[dv_names[j]];
				if (new_par_vel[j] > vmax + FLOAT_EPSILON) {
					new_par_vel[j] = vmax;
				}
				else if (new_par_vel[j] < -vmax - FLOAT_EPSILON) {
					new_par_vel[j] = -vmax;
				}

				new_dv = cur_real[j] + new_par_vel[j];
				new_par_dval[j] = new_dv;
				double curr_vel = new_par_vel[j];
				int draws = 0;
				while (true)
				{

					if (!((new_dv <= ub_val + FLOAT_EPSILON) && (new_dv >= lb_val - FLOAT_EPSILON)))
					{
						double wiggle_room = 0;
						if ((new_dv > ub_val + FLOAT_EPSILON))
							wiggle_room = ub_val - cur_real[j];
						else if ((new_dv < lb_val - FLOAT_EPSILON))
							wiggle_room = lb_val - cur_real[j];
						else
							throw_moea_error("invalid dv value in pso velocity calculation");

						draws++;
						if (draws > 1000)
						{
							ss << "problem with perturbing member: " << real_name << endl;
							ss << "at dv: " << dv_names[j] << endl;
							ss << setprecision(17) << fixed
								<< "wiggle room: " << wiggle_room << endl
								<< "inertia component: " << inertia_comp[j] << endl
								<< "cognitive component: " << cog_comp[j] << endl
								<< "social component: " << social_comp[j] << endl
								<< "current dv: " << cur_real[j] << endl
								<< "new dv: " << new_dv << endl
								<< "pbest: " << p_best[j] << endl
								<< "gbest: " << g_best[j] << endl;
							ofstream& frec = file_manager.rec_ofstream();
							frec << ss.str();
							throw_moea_error("infinite loop in pso velocity calculation (see rec file for details)");
						}

						//Adam's recursive perturbation algo to seek new feasible dv
						if (pso_dv_bound_handling == "REPERTURB" || "HYBRID")
						{
							double curr_dv = new_dv;

							vector<double> r1 = uniform_draws(1, 0.0, 1.0, rand_gen);
							vector<double> r2 = uniform_draws(1, 0.0, 1.0, rand_gen);
							vector<double> r3 = uniform_draws(1, 0.0, 1.0, rand_gen);

							//do clamping sometimes -- recommended for MOO; straight up REPERTURBATION generally performs better for SOO
							if ((2 * abs(wiggle_room) / (ub_val - lb_val)) < r3[0] - FLOAT_EPSILON &&
								(pso_dv_bound_handling == "HYBRID"))
							{
								new_dv = new_par_dval[j] < lb_val - FLOAT_EPSILON ? lb_val : new_par_dval[j];
								new_dv = new_par_dval[j] > ub_val + FLOAT_EPSILON ? ub_val : new_par_dval[j];
								new_par_dval[j] = new_dv;
								break;
							}

							inertia_comp[j] = omega * curr_vel;
							cog_comp[j] = cog_const * r1[0] * (p_best[j] - curr_dv);
							social_comp[j] = social_const * r2[0] * (g_best[j] - curr_dv);

							curr_vel = inertia_comp[j] + cog_comp[j] + social_comp[j];

							double vmax = pso_vmax[dv_names[j]];
							if (curr_vel > vmax + FLOAT_EPSILON) {
								curr_vel = vmax;
							}
							else if (curr_vel < -vmax - FLOAT_EPSILON) {
								curr_vel = -vmax;
							}
							new_dv = curr_dv + curr_vel;
							new_par_dval[j] = new_dv;
							//new_par_vel[j] = curr_vel;
							new_par_vel[j] = new_dv - cur_real[j];
						}
						else
							throw_moea_error("invalid pso_dv_bound_handling option. Choose between REPERTURB, CLAMP, or HYBRID");
					}
					else
						break;
				}
			}
		}
		new_vel.row(i) = new_par_vel;
		child_dv.row(i) = new_par_dval;
	}
	return pair<ParameterEnsemble, ParameterEnsemble>(ParameterEnsemble(&pest_scenario, &rand_gen, new_vel, _dp.get_real_names(), _dp.get_var_names()), 
		ParameterEnsemble(&pest_scenario, &rand_gen, child_dv, _dp.get_real_names(), _dp.get_var_names()));
}



vector<string> MOEA::get_pso_gbest_solutions(int num_reals, ParameterEnsemble& _dp, ObservationEnsemble& _op)
{
	stringstream ss;
	DomPair dompair = objectives.get_nsga2_pareto_dominance(-999, _op, _dp, &constraints, prob_pareto, false);
	vector<string> nondom_solutions = dompair.first;
	vector<string> gbest_solutions;
	double alpha = pest_scenario.get_pestpp_options().get_mou_pso_alpha();
    int num_objs = pi_obj_names.size()+obs_obj_names.size();

	//if no non dom solutions, then use the dominated ones...
	if (nondom_solutions.size() == 0)
	{
        ss.str("");
        ss << "WARNING: no nondom solutions for pso gbest calculation, using dominated solutions" << endl;
		nondom_solutions = dompair.second;
	}
    //todo: should we warn for nondom > 1 and objs == 1?
	else if ((nondom_solutions.size() == 1) && (num_objs > 1))
	{
	    ss.str("");
	    ss << "WARNING: only one nondom solution for pso gbest calculation" << endl;
	    file_manager.rec_ofstream() << ss.str();
	    cout << ss.str();
		for (int i = 0; i < num_reals; i++)
			gbest_solutions.push_back(nondom_solutions[0]);
		return gbest_solutions;
	}
	

	map<string, double> fitness = objectives.get_mopso_fitness(nondom_solutions, _op, _dp);
	vector<string> working;
	string candidate;
	int count = 0;
	vector < double> r;
	bool found;
	double size = nondom_solutions.size();
	for (int i = 0; i < num_reals; i++)
	{
		count = 0;
		found = false;

		
		
		while (true)
		{
			working = nondom_solutions;
			shuffle(working.begin(), working.end(), rand_gen);
			r = uniform_draws(nondom_solutions.size(), 0.0, 1.0, rand_gen);
			for (int i = 0; i < r.size(); i++)
				if (fitness[working[i]] >= r[i] - FLOAT_EPSILON)
				{
					candidate = working[i];
					found = true;
					break;
				}

			if (found)
				break;
			count++;
			if (count > 1000000) {
                throw_moea_error("MOEA::get_pso_gbest_solutions() seems to be stuck in a infinite loop....");
            }
		}
		gbest_solutions.push_back(candidate);
	}
	return gbest_solutions;
}

ParameterEnsemble MOEA::generate_pso_population(int num_members, ParameterEnsemble& _dp)
{
    //generate this first before pso resets the objectives member map...
    ParameterEnsemble temp(&pest_scenario, _dp.get_rand_gen_ptr());

    if (num_members > _dp.shape().first)
    {
        int num_reals = num_members - _dp.shape().first;
        message(1,"augmenting PSO population DE population of size",num_reals);
        temp = generate_diffevol_population(num_reals, _dp);

    }
    message(1, "generating PSO population of size", num_members);
	vector<string> gbest_solutions = get_pso_gbest_solutions(_dp.shape().first, dp_archive, op_archive);
	pair<ParameterEnsemble, ParameterEnsemble> new_gen = get_updated_pso_velocity(_dp, gbest_solutions);
	ParameterEnsemble cur_velocity = new_gen.first;
	ParameterEnsemble new_dp = new_gen.second;

    if (temp.shape().first > 0) {
        new_dp.append_other_rows(temp);

        pso_pbest_dp.append_other_rows(temp);
        vector<string> real_names = temp.get_real_names();
        ParameterEnsemble ptemp = get_initial_pso_velocities(temp.shape().first);
        ptemp.set_real_names(real_names);
        cur_velocity.append_other_rows(ptemp);
    }


    current_pso_lineage_map.clear();
	string new_name;
	vector<string> new_names;
	map<string,string> primary_parent_map;
    ofstream& lin = file_manager.get_ofstream(lineage_tag);
	for (auto real_name : new_dp.get_real_names())
	{
		new_name = get_new_member_name("pso");
		current_pso_lineage_map[real_name] = new_name;
		new_names.push_back(new_name);
		primary_parent_map[new_name] = real_name;
		lin << new_name << "," << real_name << ",," << endl;

	}
	cur_velocity.set_real_names(new_names);
	pso_velocity = cur_velocity;
    update_pso_velocity_map(pso_velocity);
	new_dp.set_real_names(new_names);
	new_dp.set_trans_status(_dp.get_trans_status());
	new_dp.enforce_bounds(performance_log,false);
	new_dp.check_for_normal("new pso population");
    if (num_members < new_dp.shape().first)
    {
        vector<string> keep;
        for (int i=0;i<num_members;i++)
            keep.push_back(new_names[i]);
        new_dp.keep_rows(keep,true);
        pso_velocity.keep_rows(keep,true);
    }
    if (_dp.get_fixed_info().get_map_size() > 0) {
        vector<string> fi_fixed_names = _dp.get_fixed_info().get_fixed_names();
        new_dp.get_fixed_info().set_fixed_names(fi_fixed_names);
        map<string, double> fi;
        vector<string> rnames = _dp.get_fixed_info().get_real_names();
        set<string> sdp_rnames(rnames.begin(),rnames.end());
        set<string>::iterator dpend = sdp_rnames.end();
        rnames = temp.get_fixed_info().get_real_names();
        set<string> stemp_rnames(rnames.begin(),rnames.end());
        set<string>::iterator tend = stemp_rnames.end();
        for (auto &p : primary_parent_map) {
            if (sdp_rnames.find(p.second) != dpend)
                fi = _dp.get_fixed_info().get_real_fixed_values(p.second);
            else if (stemp_rnames.find(p.second) != tend)
                fi = temp.get_fixed_info().get_real_fixed_values(p.second);
            else
                throw_moea_error("fixed info for existing realization '"+p.second+"' not found");
            new_dp.get_fixed_info().add_realization(p.first, fi);
        }
    }
	return new_dp;
}



ParameterEnsemble MOEA::simplex_cceua_kn(ParameterEnsemble s, int k, int optbounds)
{
	//C++ implementation of the cceua algorithm, Duan et al. (1992) with the addition of k worst points
	// 	   and n steps along the reflection path


	//TODO get npt from s.
	int nopt = 30; //number of variables in the model, in an ideal situations has nopt realizations, handle size of s in generate_simplex_population
	int nps = nopt + 1; // number of members in a simplex

	//TODO get parameters and fitness from s
	Eigen::MatrixXd svals(nps, nopt); //PARAMETERS 
	Eigen::VectorXd sfvals(nps);     //OBJECTIVE FUNCTION for members of the simplex

	//TOERASE, FILL WITH RANDOM NUMBERS FOR NOW
	for (int i = 0; i < nps; i++)
	{
		for (int j = 0; j < nopt; j++)
			svals(i, j) = uniform_draws(1, 0.0, 1.0, rand_gen)[0];
		sfvals(i) = uniform_draws(1, 0.0, 1.0, rand_gen)[0];
	}

	//TODO GET bl bu from s
	Eigen::VectorXd bl(nopt); 
	Eigen::VectorXd bu(nopt);
	for (int i = 0; i < nopt; i++)
	{
		bu(i) = 1.0; //zdt1 example
		bl(i) = 0.0; //zdt1 example

	}
	//Create vector with n steps from reflection [1, 1-1/n, 1-2/n, ...1-(n-1)/n]
	//Examples:                             n=1, [1]
	//                                      n=4, [1, 1-1/4, 1-2/4, 1-3/4]
	//TODO DECIDE TO INCLUDE one or more contraction points right the way or under some circumstance
	//A contraction point could use -1+2/n or something similar.
	vector<double> alpha_d_vec = pest_scenario.get_pestpp_options().get_mou_simplex_factors();
	int nsteps = alpha_d_vec.size();

	//initialize structures for kth worst values
	Eigen::MatrixXd skw(k, nopt);
	Eigen::VectorXd sfkw(k);

	//Separate the k worst points
	for (int ik = 0; ik < k; ik++)
	{
		skw.row(ik) = svals.row(svals.rows() - ik -1);
		sfkw(ik) = sfvals(svals.rows() - ik -1);
	}

	//Loop through the k worst points and n steps reflections
	Eigen::MatrixXd snewkn(k * alpha_d_vec.size(), nopt);
	int inew = 0;
	for (int ik = 0; ik < k; ik++)
	{
		// Compute the centroid of the simplex excluding the selected kth worst point
		Eigen::MatrixXd svalsek(nps - 1, nopt);
		int j = 0;
		for (int ikk = 0; ikk < nps; ikk++)
		{
			if (ikk != ik)
			{
				svalsek.row(j) = svals.row(ikk);
				j++;
			}

		}
		Eigen::VectorXd ce = svalsek.rowwise().mean();

		//Query reflection/contration points stored in vector of reflection/contraction  points
		for (int ia = 0; ia < alpha_d_vec.size(); ia++)
		{
			Eigen::VectorXd valskrow = skw.row(ik);
			Eigen::VectorXd delta = ce - valskrow;
			double alpha;
			alpha = alpha_d_vec[ia];
			Eigen::VectorXd delta_a = delta * alpha;
			Eigen::VectorXd ce_delta_a = ce + delta_a;
						
			//Check if is outside the bounds :
			int ibound = 0;
			
			Eigen::VectorXd sl = ce_delta_a - bl;

			if ((sl.array() < 0.0).count() > 0)
				ibound = 1;

			sl = bu - ce_delta_a;
			if ((sl.array() < 0.0).count() > 0)
				ibound = 2;

			if (ibound >= 1)
				//TODO BRING TO THE BOUND INSTEAD OF RANDOM.
				switch (optbounds){
					case 1:
						//RANDOM, ORIGINAL SCE
						ce_delta_a = bl.array() + uniform_draws(1, 0.0, 1.0, rand_gen)[0] * (bu.array() - bl.array()); //TODO CHECK RECIPE
						break;
					case 2:
						//ENFORCE BOUNDS CODE
						for (int j = 0; j < ce_delta_a.size();j++ )
						{
							if (ce_delta_a(j) > bu(j))
								ce_delta_a(j) = bu(j);
							if (ce_delta_a(j) < bl(j))
								ce_delta_a(j) = bl(j);
						}						
						break;
					case 3:
						//RANDOM, ONLY FOR THE PARAMETER BEYOND THE BOUND
						break;
				}
			snewkn.row(inew) = ce_delta_a;
			inew++; //index used to fill k x n matrix

		}

	}

	//TOD0 assign snewk back in s

	return s;//TODO fnew, icall
} 

ParameterEnsemble MOEA::generate_simplex_population(int num_members, ParameterEnsemble& _dp, ObservationEnsemble& _op)
{
    message(1, "generating simplex population of size", num_members);
    //for now just using nsga selector. TODO: work in spea2 selector if requested
    DomPair t = objectives.get_nsga2_pareto_dominance(-999,_op,_dp,&constraints,prob_pareto,false);
    vector<string> fitness;
    for (auto& tt : t.first)
        fitness.push_back(tt);
    for (auto& tt : t.second)
        fitness.push_back(tt);


	int num_reflect = min(pest_scenario.get_pestpp_options().get_mou_simplex_reflections(),_dp.shape().first-1);
	vector<double> step_factors = pest_scenario.get_pestpp_options().get_mou_simplex_factors();
    num_members = num_reflect * step_factors.size();

    Eigen::MatrixXd new_reals(num_members, _dp.shape().second);
    new_reals.setZero();
    _dp.transform_ip(ParameterEnsemble::transStatus::NUM);
    vector<string> new_member_names;

    ofstream& lin = file_manager.get_ofstream(lineage_tag);
    vector<string> real_names = _dp.get_real_names();
    string new_name,current_name;
    _dp.update_var_map();
    map<string, int> var_map = _dp.get_var_map();
    string dv_name;
    int ii;
    int i_last;
    Parameters lbnd = pest_scenario.get_ctl_parameter_info().get_low_bnd(dv_names);
    Parameters ubnd = pest_scenario.get_ctl_parameter_info().get_up_bnd(dv_names);
    ParamTransformSeq bts = pest_scenario.get_base_par_tran_seq();
    bts.ctl2numeric_ip(lbnd);
    bts.ctl2numeric_ip(ubnd);

    vector<string> best;
    for (int i=0;i<fitness.size()-num_reflect;i++)
        best.push_back(fitness[i]);

    //get the centroid in dec var space of the best individuals (excluding the ones that will be reflected)
    Eigen::VectorXd best_centroid = _dp.get_eigen(best,vector<string>()).colwise().mean();

    //to hold the current and newly reflected individuals
    Eigen::VectorXd current,reflected;

    stringstream ss;


    int i_newreal = 0;
    map<string,string> primary_parent_map;
    for (int k=0;k<num_reflect;k++)
    {
        //fitness is sorted from best to worst, so work from the end
        current_name = fitness[fitness.size()-k-1];
        //get the current worst individual
        current = _dp.get_real_vector(current_name);
        //reflect
        reflected = best_centroid - current;

        //now generate the step factor reflections
        for (auto& f : step_factors)
        {
            new_reals.row(i_newreal) = reflected * f;
            //get a new name for this individual
            ss.str("");
            ss << "simplex_k" << k << "_f" << setprecision(2) << f;
            new_name = get_new_member_name(ss.str());
            new_member_names.push_back(new_name);
            i_newreal++;
            //write the lineage info
            lin << new_name << "," << current_name << endl;
            primary_parent_map[new_name] = current_name;
        }
    }
    ParameterEnsemble new_dp(&pest_scenario, &rand_gen, new_reals, new_member_names, _dp.get_var_names());
    new_dp.set_trans_status(ParameterEnsemble::transStatus::NUM);
    if (pest_scenario.get_pestpp_options().get_mou_simplex_mutation())
    {
        gauss_mutation_ip(new_dp);
    }
    new_dp.enforce_bounds(performance_log, false);
    if (_dp.get_fixed_info().get_map_size() > 0) {
        vector<string> fi_fixed_names = _dp.get_fixed_info().get_fixed_names();
        new_dp.get_fixed_info().set_fixed_names(fi_fixed_names);
        map<string, double> fi;
        for (auto &p : primary_parent_map) {
            fi = _dp.get_fixed_info().get_real_fixed_values(p.second);
            new_dp.get_fixed_info().add_realization(p.first, fi);
        }
    }
    return new_dp;
}


ParameterEnsemble MOEA::generate_diffevol_population(int num_members, ParameterEnsemble& _dp)
{
	/* adaptive idea:  Front. Built Environ., 09 July 2020 | https://doi.org/10.3389/fbuil.2020.00102
	A Comparative Study of Differential Evolution Variants in Constrained Structural Optimization
	Manolis Georgioudakis1 and Vagelis Plevris2**/
	message(1, "generating diffevol population of size", num_members);
	vector<int> r_int_vec;
	for (int i = 0; i < dv_names.size() - n_adaptive_dvs; i++)
		r_int_vec.push_back(i);


	Eigen::VectorXd y, x, diff;
	double org_crossover = pest_scenario.get_pestpp_options().get_mou_crossover_probability();
	double r;
	
	vector<double> cr_vals;
	Eigen::MatrixXd new_reals(num_members, _dp.shape().second);
	new_reals.setZero();
	_dp.transform_ip(ParameterEnsemble::transStatus::NUM);
	vector<string> new_member_names;
	
	//since _dp might contain both dev vars and pars, we want to 
	//make sure we are only fooling with dec vars
	//the var_map lets us map between dv name and col index
	ofstream& lin = file_manager.get_ofstream(lineage_tag);
	vector<string> real_names = _dp.get_real_names();
	string new_name;
	_dp.update_var_map();
	map<string, int> var_map = _dp.get_var_map();
	string dv_name;
	int ii;
	int i_last;
	vector<int> selected;
	double org_f = pest_scenario.get_pestpp_options().get_mou_de_f();
	double f = org_f, crossover = org_crossover;
	bool adaptive_f = false;
	bool adaptive_cr = false;
	if (var_map.find(DE_F_NAME) != var_map.end())
		adaptive_f = true;
	if (var_map.find(CR_NAME) != var_map.end())
		adaptive_cr = true;
	Parameters lbnd = pest_scenario.get_ctl_parameter_info().get_low_bnd(dv_names);
	Parameters ubnd = pest_scenario.get_ctl_parameter_info().get_up_bnd(dv_names);
	ParamTransformSeq bts = pest_scenario.get_base_par_tran_seq();
	bts.ctl2numeric_ip(lbnd);
	bts.ctl2numeric_ip(ubnd);
	int tries = 0;
	int i = 0;
	map<string,string> primary_parent_map;
	while (i < num_members)
    {
		selected = selection(4, _dp, mattype);
		if (adaptive_f)
		{
			if (uniform_draws(1, 0.0, 1.0, rand_gen)[0] > 0.1)
				f = _dp.get_eigen_ptr()->row(selected[0])[var_map[DE_F_NAME]];
			else
				f = org_f;
		}
		if (adaptive_cr)
		{
			if (uniform_draws(1, 0.0, 1.0, rand_gen)[0] > 0.1)
				crossover = _dp.get_eigen_ptr()->row(selected[0])[var_map[CR_NAME]];
			else
				crossover = org_crossover;
		}
		
		//differential vector
		diff = _dp.get_eigen_ptr()->row(selected[0]) + (f * (_dp.get_eigen_ptr()->row(selected[1]) - _dp.get_eigen_ptr()->row(selected[2])));

		//current member if in range, otherwise, select randomly
		if (i < _dp.shape().first)
			x = _dp.get_eigen_ptr()->row(i);
		else
			//this risks "inbreeding" but maybe that's good?!
			x = _dp.get_eigen_ptr()->row(selected[3]);
		//copy to preserve non-dec var values;
		y = x; 
		//random cross over probs - one per decision variable
		cr_vals = uniform_draws(_dp.shape().second, 0.0, 1.0, rand_gen);
		
		//get the R values for this member;
		shuffle(r_int_vec.begin(), r_int_vec.end(), rand_gen);
		r = r_int_vec[0];

		//only change dec vars
		double dsum = 0.0;
		for(int idv=0;idv<dv_names.size();idv++)
		{
			ii = var_map[dv_names[idv]];
			if ((cr_vals[ii] < crossover) || (ii == r))
			{
				y[ii] = diff[ii];
				y[ii] = min(y[ii], ubnd[dv_names[idv]]);
				y[ii] = max(y[ii], lbnd[dv_names[idv]]);

				
			}
		}

		tries++;
		if (tries > 1000000)
			throw_moea_error("diffevol generator appears to be stuck in an infinite loop...");

		double n = (x - y).squaredNorm();
		if (n < epsilon)
			continue;
		new_name = get_new_member_name("de");
		new_member_names.push_back(new_name);
        primary_parent_map[new_name] = real_names[selected[0]];
		lin << new_name;
		for (auto idx : selected)
			lin << "," << real_names[idx];
		lin << endl;
		new_reals.row(i) = y;
		i++;
	}

	ParameterEnsemble new_dp(&pest_scenario, &rand_gen, new_reals, new_member_names, _dp.get_var_names());
	
	new_dp.set_trans_status(ParameterEnsemble::transStatus::NUM);
	new_dp.enforce_bounds(performance_log, false);
	//transfer any fixed par info
	if (_dp.get_fixed_info().get_map_size() > 0)
    {
	    vector<string> fi_fixed_names = _dp.get_fixed_info().get_fixed_names();
	    new_dp.get_fixed_info().set_fixed_names(fi_fixed_names);
	    map<string,double> fi;
	    for (auto& p : primary_parent_map)
        {
	        fi = _dp.get_fixed_info().get_real_fixed_values(p.second);
	        new_dp.get_fixed_info().add_realization(p.first,fi);
        }
    }
	return new_dp;
}

ParameterEnsemble MOEA::generate_pm_population(int num_members, ParameterEnsemble& _dp)
{
	message(1, "generating PM population of size", num_members);
	_dp.transform_ip(ParameterEnsemble::transStatus::NUM);
	Parameters lbnd = pest_scenario.get_ctl_parameter_info().get_low_bnd(dv_names);
	Parameters ubnd = pest_scenario.get_ctl_parameter_info().get_up_bnd(dv_names);
	pest_scenario.get_base_par_tran_seq().ctl2numeric_ip(lbnd);
	pest_scenario.get_base_par_tran_seq().ctl2numeric_ip(ubnd);
	vector<int> selected;
	Eigen::VectorXd child, parent;
	double org_mut_prob = pest_scenario.get_pestpp_options().get_mou_mutation_probability();
	if (org_mut_prob < -0.0)
		org_mut_prob = 1.0 / (double)_dp.shape().second;
	
	double disrupt_prob = 0.2;
	vector<string> _dv_names = _dp.get_var_names();
	vector<string> real_names = _dp.get_real_names();
	Eigen::MatrixXd new_reals(num_members, _dp.shape().second);
	new_reals.setZero();
	int tries = 0;
	int imember = 0;
	ofstream& lin = file_manager.get_ofstream(lineage_tag);
	vector<string> new_member_names;
	string new_name;

	_dp.update_var_map();
	map<string, int> var_map = _dp.get_var_map();
	bool adaptive_mr = false;
	if (var_map.find(MR_NAME) != var_map.end())
		adaptive_mr = true;
	double mut_prob = org_mut_prob;
	map<string,string> primary_parent_map;
	while (imember < num_members)
	{
		selected = selection(1, _dp, mattype);
		parent = _dp.get_real_vector(selected[0]);
		if (adaptive_mr)
			mut_prob = parent[var_map[MR_NAME]];
		child = hybrid_pm(parent, mut_prob, disrupt_prob, _dv_names, lbnd, ubnd);
		if ((parent - child).squaredNorm() < epsilon)
			continue;
		
		tries++;
		if (tries > 10000000)
			throw_moea_error("hybrid polynomial mutation appears to be stuck in an infinite loop...");
		
		new_reals.row(imember) = child;
		imember++;
		new_name = get_new_member_name("pm");
		new_member_names.push_back(new_name);
		lin << new_name;
		for (auto idx : selected)
			lin << "," << real_names[idx];
		lin << endl;
		primary_parent_map[new_name] = real_names[selected[0]];
	}

	ParameterEnsemble tmp_dp(&pest_scenario, &rand_gen, new_reals, new_member_names, _dp.get_var_names());
	tmp_dp.set_trans_status(ParameterEnsemble::transStatus::NUM);
	tmp_dp.enforce_bounds(performance_log,false);
    if (_dp.get_fixed_info().get_map_size() > 0) {
        vector<string> fi_fixed_names = _dp.get_fixed_info().get_fixed_names();
        tmp_dp.get_fixed_info().set_fixed_names(fi_fixed_names);
        map<string, double> fi;
        for (auto &p : primary_parent_map) {
            fi = _dp.get_fixed_info().get_real_fixed_values(p.second);
            tmp_dp.get_fixed_info().add_realization(p.first, fi);
        }
    }

	return tmp_dp;
}

vector<int> MOEA::selection(int num_to_select, ParameterEnsemble& _dp, MouMateType& _mattype)
{
	int i_member = 0, p1_idx,p2_idx;
	vector<int> member_count, working_count, selected, r_int_vec;
	vector<double> rnds;
	for (int i = 0; i < _dp.shape().first; i++)
		member_count.push_back(i);
	vector<string> real_names = _dp.get_real_names();
	for (int i = 0; i < _dp.shape().first; i++)
		r_int_vec.push_back(i);
	set<int> selected_members;
	string s1, s2;
	int tries = 0;
	while (selected_members.size() < num_to_select)
	{
		working_count = member_count;//copy member count index to working count
		// sampled with replacement
		shuffle(working_count.begin(), working_count.end(), rand_gen); //randomly shuffle working count
		//just take the first two since this should change each time thru
		p1_idx = working_count[0];
		if (selected_members.find(p1_idx) != selected_members.end())
			continue;
		p2_idx = working_count[1];
		if (selected_members.find(p2_idx) != selected_members.end())
			continue;
		s1 = real_names[p1_idx];
		s2 = real_names[p2_idx];
		if (_mattype==MouMateType::TOURNAMENT)
		{
			if (objectives.compare_two(s1, s2, envtype))
				selected_members.emplace(p1_idx);
			else
				selected_members.emplace(p2_idx);
		}
		else if (_mattype == MouMateType::RANDOM)
		{
			selected_members.emplace(p1_idx);
			if (selected_members.size() == num_to_select)
				break;
			selected_members.emplace(p2_idx);
		}
		else
		{
			throw_moea_error("selector error: unrecognized MouMatingType, should be 'random' or 'tournament'");
		}
		tries++;
		if (tries > 1000000000)
			throw_moea_error("selection process appears to be stuck in an infinite loop...");
	}
	vector<int> members(selected_members.begin(), selected_members.end());
	return members;
}

ParameterEnsemble MOEA::generate_sbx_population(int num_members, ParameterEnsemble& _dp)
{
	message(1, "generating SBX population of size", num_members);

	_dp.transform_ip(ParameterEnsemble::transStatus::NUM);

	vector<int> r_int_vec;
	vector<double> rnds;
	
	for (int i = 0; i < dv_names.size(); i++)
		r_int_vec.push_back(i);

	double org_crossover_probability = pest_scenario.get_pestpp_options().get_mou_crossover_probability();
	double crossover_distribution_index = 10.0;
	double crossover = org_crossover_probability;
	int i_member = 0;
	int p1_idx, p2_idx;
	Eigen::MatrixXd new_reals(num_members, _dp.shape().second);
	new_reals.setZero();
	pair<Eigen::VectorXd, Eigen::VectorXd> children;
	vector<string> new_member_names;
	vector<int> selected;
	ofstream& lin = file_manager.get_ofstream(lineage_tag);
	vector<string> real_names = _dp.get_real_names();
	vector<string> _dv_names = _dp.get_var_names();
	string new_name;

	Parameters lbnd = pest_scenario.get_ctl_parameter_info().get_low_bnd(dv_names);
	Parameters ubnd = pest_scenario.get_ctl_parameter_info().get_up_bnd(dv_names);
	pest_scenario.get_base_par_tran_seq().ctl2numeric_ip(lbnd);
	pest_scenario.get_base_par_tran_seq().ctl2numeric_ip(ubnd);
	Eigen::VectorXd parent1, parent2;
	int tries = 0;
	_dp.update_var_map();
	map<string, int> var_map = _dp.get_var_map();
	bool adaptive_cr = false;
	if (var_map.find(CR_NAME) != var_map.end())
		adaptive_cr = true;
	map<string,string> primary_parent_map;
	while (i_member < num_members)
	{
		selected = selection(2, _dp, mattype);
		p1_idx = selected[0];
		p2_idx = selected[1];

		//generate two children thru cross over
		parent1 = _dp.get_real_vector(p1_idx);
		parent2 = _dp.get_real_vector(p2_idx);
		if ((parent1 - parent2).squaredNorm() < epsilon)
			continue;
		if (adaptive_cr)
		{
			if (uniform_draws(1, 0.0, 1.0, rand_gen)[0] > 0.1)
			{
				if (uniform_draws(1, 0.0, 1.0, rand_gen)[0] > 0.5)
					crossover = parent1[var_map[CR_NAME]];
				else
					crossover = parent2[var_map[CR_NAME]];
			}
			else
				crossover = org_crossover_probability;
		}
		children = sbx_new(crossover, crossover_distribution_index, parent1, parent2,_dv_names, lbnd,ubnd);

		//put the two children into the child population
		new_reals.row(i_member) = children.first;
		new_name = get_new_member_name("sbx");
		lin << new_name << "," << real_names[p1_idx] << "," << real_names[p2_idx] << endl;
		new_member_names.push_back(new_name);
		//cout << i_member << "," << p1_idx << "," << p2_idx << new_names[new_names.size() - 1] << endl;
		i_member++;
		if (i_member >= num_members)
			break;
		new_reals.row(i_member) = children.second;
		new_name = get_new_member_name("sbx");
		lin << new_name << "," << real_names[p1_idx] << "," << real_names[p2_idx] << endl;
		new_member_names.push_back(new_name);
		primary_parent_map[new_name] = real_names[p1_idx];
		//cout << i_member << "," << p1_idx << "," << p2_idx << new_names[new_names.size() -1] << endl;
		i_member++;
		tries++;
		if (tries > 1000000)
			throw_moea_error("sbx appears to be stuck in an infinite loop");
	}
		
	ParameterEnsemble tmp_dp(&pest_scenario, &rand_gen, new_reals, new_member_names, _dp.get_var_names());
	tmp_dp.set_trans_status(ParameterEnsemble::transStatus::NUM);

	//if (find(gen_types.begin(),gen_types.end(),MouGenType::PM) == gen_types.end())
	gauss_mutation_ip(tmp_dp);
	//generate_pm_population(tmp_dp.shape().first, tmp_dp);
	
	tmp_dp.enforce_bounds(performance_log,false);
    if (_dp.get_fixed_info().get_map_size() > 0) {
        vector<string> fi_fixed_names = _dp.get_fixed_info().get_fixed_names();
        tmp_dp.get_fixed_info().set_fixed_names(fi_fixed_names);
        map<string, double> fi;
        for (auto &p : primary_parent_map) {
            fi = _dp.get_fixed_info().get_real_fixed_values(p.second);
            tmp_dp.get_fixed_info().add_realization(p.first, fi);
        }
    }
	return tmp_dp;
}

Ensemble MOEA::save_pi_constraints(ParameterEnsemble &_dp, vector<string> &pinames)
{
	//check if there are prior info equations
	//if (constraints.num_pi_constraints() > 0)
	//	{
		// get the parameter and realization names
		vector<string> parnames;
		parnames = _dp.get_var_names();
		vector<string> realnames;
		realnames = _dp.get_real_names();

		// empty Parameters and parvals vector for populating later
		Parameters pars;
		Eigen::VectorXd parvals;

		// pi obj for later
		PriorInformation constraints_pi = pest_scenario.get_prior_info();

		//Instantiate an empty Ensemble class with reserved row and col names
		Ensemble pioe(&pest_scenario);
		pioe.reserve(realnames,pinames);
		//pioe.update_var_map();
		Eigen::MatrixXd piX(realnames.size(),pinames.size());

		//piX.resize(realnames.size(),pinames.size());

		int i = 0;
		for (auto &real : realnames)
			{
			// get Parameters for each realization
			parvals = _dp.get_real_vector(real);
			pars.update_without_clear(parnames,parvals);

			Eigen::VectorXd pivals;
			pivals.resize(pinames.size());
			int j = 0;
			//loop over pi_constraint_names
			for (auto &piname : pinames)
				{
				PriorInformationRec pi_rec = constraints_pi.get_pi_rec(piname);

				pair<double,double> pi_sim_resid = pi_rec.calc_sim_and_resid(pars);
				// update the matrix
				piX(i,j) = pi_sim_resid.first;
				j++;
				}
			i++;
			}
	pioe.update_var_map();
	pioe.from_eigen_mat(piX,realnames,pinames);

	//}
	pioe.update_var_map();
	return pioe;
}


void MOEA::save_populations(ParameterEnsemble& _dp, ObservationEnsemble& _op, string tag, bool force_save)
{
	
	stringstream ss;
	string fname;
    _dp.reset_org_real_names();
	ss << file_manager.get_base_filename();
	if (tag.size() > 0)
	{
		ss << "." << tag;
	}
	ss << "." << dv_pop_file_tag;
	if (pest_scenario.get_pestpp_options().get_save_binary())
	{
        if (pest_scenario.get_pestpp_options().get_save_dense())
        {
            ss << ".bin";
            _dp.to_dense(ss.str());
        }
        else {
            ss << ".jcb";
            _dp.to_binary(ss.str());
        }
	}
	else
	{
		ss << ".csv";
		_dp.to_csv(ss.str());
	}
	string name = ss.str();
	ss.str("");
	ss << "saved decision variable population of size " << _dp.shape().first << " X " << _dp.shape().second << " to '" << name << "'";
	message(1, ss.str());
	ss.str("");
	if (((save_every > 0) && (iter % save_every == 0)) || (iter == pest_scenario.get_control_info().noptmax) || (force_save))
	{
		ss << file_manager.get_base_filename() << "." << iter;
		if (tag.size() > 0)
		{
			ss << "." << tag;
		}
		ss << "." << dv_pop_file_tag;
		if (pest_scenario.get_pestpp_options().get_save_binary()) {
            if (pest_scenario.get_pestpp_options().get_save_dense()) {
                ss << ".bin";
                _dp.to_dense(ss.str());
            } else {
                ss << ".jcb";
                _dp.to_binary(ss.str());
            }
        }
		else
		{
			ss << ".csv";
			_dp.to_csv(ss.str());
		}
		string name = ss.str();
		ss.str("");
		ss << "saved generation-specific decision variable population of size " << _dp.shape().first << " X " << _dp.shape().second << " to '" << name << "'";
		message(1, ss.str());
	}

	// record prior information constraint obs names
	vector<string> pinames;
	pinames = constraints.get_pi_constraint_names();
	// remove __risk__ from pinames
	pinames.erase(std::remove(pinames.begin(), pinames.end(), "_RISK_"), pinames.end());

	if (pinames.size() >0)
	{
		Ensemble _dpi= save_pi_constraints(_dp, pinames);

		ss.str("");
		ss << file_manager.get_base_filename();
		if (tag.size() > 0)
		{
			ss << "." << tag;
		}

		ss << "." << pi_pop_file_tag;
		if (pest_scenario.get_pestpp_options().get_save_binary())
		{
			if (pest_scenario.get_pestpp_options().get_save_dense())
			{
				ss << ".bin";
				_dpi.to_dense(ss.str());
			}
			else {
				ss << ".jcb";
				_dpi.to_binary(ss.str());
			}
		}
		else
		{
			ss << ".csv";
			_dpi.to_csv(ss.str());
		}
		string name = ss.str();
		ss.str("");
		ss << "saved prior information population of size " << _dpi.shape().first << " X " << _dpi.shape().second << " to '" << name << "'";
		message(1, ss.str());
		ss.str("");
		if (((save_every > 0) && (iter % save_every == 0)) || (iter == pest_scenario.get_control_info().noptmax) || (force_save))
		{
			ss << file_manager.get_base_filename() << "." << iter;
			if (tag.size() > 0)
			{
				ss << "." << tag;
			}
			ss << "." << pi_pop_file_tag;
			if (pest_scenario.get_pestpp_options().get_save_binary()) {
				if (pest_scenario.get_pestpp_options().get_save_dense()) {
					ss << ".bin";
					_dpi.to_dense(ss.str());
				} else {
					ss << ".jcb";
					_dpi.to_binary(ss.str());
				}
			}
			else
			{
				ss << ".csv";
				_dpi.to_csv(ss.str());
			}
			string name = ss.str();
			ss.str("");
			ss << "saved generation-specific prior information population of size " << _dpi.shape().first << " X " << _dpi.shape().second << " to '" << name << "'";
			message(1, ss.str());
		}
	}


	ss.str("");
	ss << file_manager.get_base_filename();
	if (tag.size() > 0)
	{
		ss << "." << tag;
	}
	ss << "." << obs_pop_file_tag;
	if (pest_scenario.get_pestpp_options().get_save_binary())
	{
        if (pest_scenario.get_pestpp_options().get_save_dense())
        {
            ss << ".bin";
            _op.to_dense(ss.str());
        }
        else {
            ss << ".jcb";
            _op.to_binary(ss.str());
        }
	}
	else
	{
		ss << ".csv";
		_op.to_csv(ss.str());
	}
	name = ss.str();
	ss.str("");
	ss << "saved observation population of size " << _op.shape().first << " X " << _op.shape().second << " to '" << name << "'";
	message(1, ss.str());

	if ((save_every > 0) && (iter % save_every == 0))
	{
		ss.str("");
		ss << file_manager.get_base_filename() << "." << iter;
		if (tag.size() > 0)
		{
			ss << "." << tag;
		}
		ss << "." << obs_pop_file_tag;
		if (pest_scenario.get_pestpp_options().get_save_binary())
		{
            if (pest_scenario.get_pestpp_options().get_save_dense())
            {
                ss << ".bin";
                _op.to_dense(ss.str());
            }
            else {
                ss << ".jcb";
                _op.to_binary(ss.str());
            }
		}
		else
		{
			ss << ".csv";
			_op.to_csv(ss.str());
		}
		name = ss.str();
		ss.str("");
		ss << "saved generation-specific observation population of size " << _op.shape().first << " X " << _op.shape().second << " to '" << name << "'";
		message(1, ss.str());
	}

}

string MOEA::get_new_member_name(string tag)
{
	stringstream ss;
	ss << "gen=" << iter+restart_iter_offset << "_member=" << member_count;
	if (tag.size() > 0)
	{
		ss << "_" << tag;
	}
	member_count++;
	return pest_utils::upper_cp(ss.str());
}

pair<Eigen::VectorXd, Eigen::VectorXd> MOEA::sbx(double probability, double di, int idx1, int idx2)
{
	int i;
	//double rnd1, rnd2, rnd3, rnd4;
	vector<double> rnds;
	double y1, y2, yL, yu;
	double c1, c2;
	double alpha, beta, betaq;
	double valueX1, valueX2, valueX1initial, valueX2initial;
	const double EPS = 1.0e-14;
	// get parents from dp
	Eigen::VectorXd x1 = dp.get_eigen_ptr()->row(idx1); // parent #1
	Eigen::VectorXd x2 = dp.get_eigen_ptr()->row(idx2); // parent #2
	// initialize offspring
	Eigen::VectorXd offs1(x1.size()); // ofspring #1
	Eigen::VectorXd offs2(x2.size()); // offspring #2

	//ZQ: get upper and lower bounds for each variable
	vector<string> var_names = dp.get_var_names();
	string vname;
	Parameters lbnd = pest_scenario.get_ctl_parameter_info().get_low_bnd(var_names);
	Parameters ubnd = pest_scenario.get_ctl_parameter_info().get_up_bnd(var_names);
	pest_scenario.get_base_par_tran_seq().ctl2numeric_ip(lbnd);
	pest_scenario.get_base_par_tran_seq().ctl2numeric_ip(ubnd);

	offs1.setZero();
	offs2.setZero();

	int n_var = dv_names.size();
	//can't set all rnds outside of loop or all vars will be treated the same
	rnds = uniform_draws(4, 0.0, 1.0, rand_gen);

	//if (rnds[0] <= probability) 
	if (true)
	{
		for (i = 0; i < n_var; i++) {
			valueX1 = x1[i];
			valueX2 = x2[i];
			valueX1initial = valueX1;
			valueX2initial = valueX2;
			vname = var_names[i];
			// cout << "valueX1 "<<valueX1<<" valueX2 "<<valueX2<<std::endl;
			if (rnds[1] <= 0.5) {
				if (fabs(valueX1 - valueX2) > EPS) {

					if (valueX1 < valueX2) {
						y1 = valueX1;
						y2 = valueX2;
					}
					else {
						y1 = valueX2;
						y2 = valueX1;
					} // if                       

					yL = lbnd[vname];
					yu = ubnd[vname];

					//cout << yL << " " << yu << endl;
					//cout <<"yL "<< yL << " yu "<<yu<<endl;

					beta = 1.0 + (2.0 * (y1 - yL) / (y2 - y1));
					//cout <<"beta "<< beta<<endl;
					alpha = 2.0 - pow(beta, -(di + 1.0));
					//cout <<"alpha "<< beta<<endl;
					//cout <<"distributionIndex_ "<< distributionIndex_<<endl;
					if (rnds[2] <= (1.0 / alpha)) {
						betaq = pow((rnds[2] * alpha), (1.0 / (di + 1.0)));
					}
					else {
						betaq = pow((1.0 / (2.0 - rnds[2] * alpha)), (1.0 / (di + 1.0)));
					} // if
					//cout <<"betaq "<< betaq<<endl;
					//cout <<"rand "<< rand<<endl;
					c1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
					beta = 1.0 + (2.0 * (yu - y2) / (y2 - y1));
					alpha = 2.0 - pow(beta, -(di + 1.0));
					//cout <<"beta2 "<< beta<<endl;
					//cout <<"alpha2 "<< alpha<<endl;
					if (rnds[2] <= (1.0 / alpha)) {
						betaq = pow((rnds[2] * alpha), (1.0 / (di + 1.0)));
					}
					else {
						betaq = pow((1.0 / (2.0 - rnds[2] * alpha)), (1.0 / (di + 1.0)));
					} // if
					//cout <<"betaq2 "<< betaq<<endl;
					c2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1));

					if (c1 < yL)
						c1 = yL;

					if (c2 < yL)
						c2 = yL;

					if (c1 > yu)
						c1 = yu;

					if (c2 > yu)
						c2 = yu;

					if (rnds[3] <= 0.5)
					{
						if (std::isnan(c2) || std::isnan(c1)) {
							c2 = valueX1initial;
							c1 = valueX2initial;
						}
						offs1[i] = c2;
						offs2[i] = c1;
					}
					else {
						if (std::isnan(c2) || std::isnan(c1)) { // sometimes, c1 or c2 can be nan, so in this case use the values from the parents
							c1 = valueX1initial;
							c2 = valueX2initial;
						}
						offs1[i] = c1;
						offs2[i] = c2;
					}
				}
				else {
					offs1[i] = valueX1;
					offs2[i] = valueX2;
				}
			}
			else {
				offs1[i] = valueX2;
				offs2[i] = valueX1;
			}
		}
	}
	pair<Eigen::VectorXd, Eigen::VectorXd> vec_pair(offs1, offs2);
	return vec_pair;
}


Eigen::VectorXd MOEA::hybrid_pm(Eigen::VectorXd& parent, double mutation_probability, double disrupt_probabilty,
	vector<string>& _dv_names, Parameters& lbnd, Parameters& ubnd)
{
	/*
	On the Disruption-level of Polynomial Mutation for Evolutionary Multi-objective Optimisation Algorithms.

    December 2009Computing and Informatics 29(5):783-800

    SourceDBLP

    Mohammad HamdanMohammad Hamdan
	
	*/
	stringstream ss;
	vector<double> rnds;
	Eigen::VectorXd child = parent; // parent #1
	double delta1, delta2, delta, deltaq1,deltaq2,deltaq;
	string vname;
	double nm = 20;
	for (int i = 0; i < child.size(); i++)
	{
		vname = _dv_names[i];
		rnds = uniform_draws(3, 0.0, 1.0, rand_gen);
		if (rnds[0] < mutation_probability)
		{
			delta1 = (child[i] - lbnd[vname]) / (ubnd[vname] - lbnd[vname]);
			delta2 = (ubnd[vname] - child[i]) / (ubnd[vname] - lbnd[vname]);
			if (rnds[1] > disrupt_probabilty)
				delta = min(delta1, delta2);
			else
			{
				if (rnds[2] <= 0.5)
					delta = delta1;
				else
					delta = delta2;
			}
			deltaq1 = pow((2.0 * rnds[2]) + ((1.0 - (2.0 * rnds[2])) * pow(1 - delta, nm + 1.0)), (1.0 / (nm + 1))) - 1.0;
			deltaq2 = 1.0 - pow((2.0 * (1 - rnds[2])) + (2.0 * (rnds[2] - 0.5) * pow(1 - delta, nm + 1.0)), 1.0 / (nm + 1));
			if (rnds[2] <= 0.5)
			{
				deltaq = deltaq1;
			}
			else
			{
				deltaq = deltaq2;
			}
			child[i] = child[i] + (deltaq * (ubnd[vname] - lbnd[vname]));

		}
	}
	return child;
}

//pair<Eigen::VectorXd, Eigen::VectorXd> MOEA::sbx_new(double crossover_probability, double di, Eigen::VectorXd& parent1, 
//	Eigen::VectorXd parent2, vector<string>& _dv_names, Parameters& lbnd, Parameters& ubnd)
//{
//
//	stringstream ss;
//	int i;
//	//double rnd1, rnd2, rnd3, rnd4;
//	vector<double> rnds;
//	double lt, ut;
//	double p1, p2, c1, c2;
//	// get parents from dp
//	Eigen::VectorXd child1 = parent1; // parent #1
//	Eigen::VectorXd child2 = parent2; // parent #2
//	string vname;
//	pair<double, double> betas;
//
//	int n_var = _dv_names.size();
//	//can't set all rnds outside of loop or all vars will be treated the same
//	rnds = uniform_draws(4, 0.0, 1.0, rand_gen);
//
//	int tries = 0;
//	double abs_diff;
//	while (true)
//	{
//		
//		child1 = parent1; 
//		child2 = parent2; 
//		for (i = 0; i < n_var; i++)
//		{
//			rnds = uniform_draws(1, 0.0, 1.0, rand_gen);
//			p1 = parent1[i];
//			p2 = parent2[i];
//			vname = _dv_names[i];
//			if (rnds[0] <= crossover_probability)
//			{
//				abs_diff = abs(p1 - p2);
//				if (abs_diff < epsilon)
//					abs_diff = 1e-10;
//				lt = (p1 + p2 + (2 * lbnd[vname]))/abs_diff;
//				ut = ((2. * ubnd[vname]) - p1 - p2) / abs_diff;
//				lt = max(lt, 1.0);
//				ut = max(ut, 1.0);
//				/*if (lt < 0)
//					throw_moea_error("sbx error: lower transform bound less than zero");
//				if (ut < 0)
//				*/	//throw_moea_error("sbx error: upper transform bound less than zero");
//				betas = get_betas(lt, ut, di);
//				c1 = 0.5 * ((p1 + p2) - betas.first * abs_diff);
//				c2 = 0.5 * ((p1 + p2) - betas.second * abs_diff);
//				if (isnan(c1))
//				{
//					ss.str("");
//					ss << "sbx error: denormal value generated for " << vname << ", beta1: " << betas.first << ", lt:" << lt << ", ut: " << ut << ", v1:" << p1 << ", v2: " << p2;
//					throw_moea_error(ss.str());
//				}
//				if (isnan(c2))
//				{
//					ss.str("");
//					ss << "sbx error: denormal value generated for " << vname << ", beta2: " << betas.second << ", lt:" << lt << ", ut: " << ut << ", v1:" << p1 << ", v2: " << p2;
//					throw_moea_error(ss.str());
//				}
//				child1[i] = c1;
//				child2[i] = c2;
// 			}
//		}
//		tries++;
//		if (tries > 10000000)
//			throw_moea_error("sbx generation process appears to be stuck in an infinite loop...");
//		//hack alert - we dont wanna waste time with identical solutions, so we will keep trying 
//		//until both child1 and child2 are different from parents - evolution!
//		if ((parent1 - child1).squaredNorm() < epsilon)
//			continue;
//		else if ((parent2 - child2).squaredNorm() < epsilon)
//			continue;
//		else
//			break;
//			
//	}
//	
//	return pair<Eigen::VectorXd, Eigen::VectorXd>(child1, child2);
//
//}

pair<Eigen::VectorXd, Eigen::VectorXd> MOEA::sbx_new(double crossover_probability, double di, Eigen::VectorXd& parent1,
	Eigen::VectorXd parent2, vector<string>& _dv_names, Parameters& lbnd, Parameters& ubnd)
{
	/* https://gist.github.com/Tiagoperes/1779d5f1c89bae0cfdb87b1960bba36d */

	stringstream ss;
	vector<double> rnds;
	double p1, p2, c1, c2;
	// get parents from dp
	Eigen::VectorXd child1 = parent1; // parent #1
	Eigen::VectorXd child2 = parent2; // parent #2
	string vname;

	int n_var = _dv_names.size();
	
	int tries = 0;
	double alpha, beta, betaq;
	while (true)
	{

		child1 = parent1;
		child2 = parent2;
		for (int i = 0; i < n_var; i++)
		{
			rnds = uniform_draws(4, 0.0, 1.0, rand_gen);
			p1 = parent1[i];
			p2 = parent2[i];
			vname = _dv_names[i];
			//this uses crossover like de instead of like standard sbx but we need to generate
			//a full child population for testing...
			if (rnds[0] <= crossover_probability)
			{			
				if (abs(p1 - p2) > epsilon)
				{
					get_sbx_child_values(p1, p2, lbnd[vname], ubnd[vname], di, rnds[1], c1, c2);
					child1[i] = c1;
					child2[i] = c2;
				}
			}
		}
		tries++;
		if (tries > 10000000)
			throw_moea_error("sbx generation process appears to be stuck in an infinite loop...");
		//hack alert - we dont wanna waste time with identical solutions, so we will keep trying 
		//until both child1 and child2 are different from parents - evolution!
		if ((parent1 - child1).squaredNorm() < epsilon)
			continue;
		else if ((parent2 - child2).squaredNorm() < epsilon)
			continue;
		else
			break;

	}
	return pair<Eigen::VectorXd, Eigen::VectorXd>(child1, child2);

}

void MOEA::get_sbx_child_values(const double& p1, const double& p2, const double& lbnd, const double& ubnd, const double& eta, double& rnd, double& c1, double& c2)
{
	double y1 = min(p1, p2);
	double y2 = max(p1, p2);
	double beta = 1.0 + (2.0 * (y1 - lbnd) / (y2 - y1));
	double alpha = 2.0 - pow(beta, -1.0 * (eta + 1));
	double betaq;
	
	if (rnd <= (1.0 / alpha))
	{
		betaq = pow((rnd * alpha), (1.0 / (eta + 1.0)));
	}
	else
	{
		betaq = pow((1.0 / (2.0 - rnd * alpha)), (1.0 / (eta + 1.0)));
	}
	c1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
	c1 = max(lbnd, c1);
	c1 = min(ubnd, c1);

	beta = 1.0 + (2.0 * (ubnd - y2) / (y2 - y1));
	alpha = 2.0 - pow(beta, -1.0 * (eta + 1.0));
	if (rnd <= (1.0 / alpha))
	{
		betaq = pow((rnd * alpha), (1.0 / (eta + 1.0)));
	}
	else
	{
		betaq = pow((1.0 / (2.0 - (rnd * alpha))), (1.0 / (eta + 1.0)));
	}
	c2 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
	c2 = max(lbnd, c2);
	c2 = min(ubnd, c2);


}
pair<double,double> MOEA::get_betas(double v1, double v2, double distribution_index)
{
	stringstream ss;
	bool are_close = false;
	if (abs(v1 - v2) <= epsilon)
		are_close = true;
	double p1,p2,u1,u2,b1,b2;	
	vector<double> rands = uniform_draws(2, 0.0, 1.0,rand_gen);
	if (are_close)
		p1 = 1., p2 = 1.;
	else
	{
		p1 = 1. - (1. / (2. * pow(v1, distribution_index)));
		p2 = 1. - (1. / (2. * pow(v2, distribution_index)));
	}
	u1 = p1 * rands[0];
	u2 = p2 * rands[1];
	if (u1 <= 0.5)
	{
		b1 = pow((2. * u1), (1. / distribution_index + 1));
	}
	else
	{
		b1 = pow((1. / (2. - (2. * u1))), (1.0 / (distribution_index + 1)));
	}
	if (u2 <= 0.5)
	{
		b2 = pow((2. * u2), (1. / distribution_index + 1));
	}
	else
	{
		b2 = pow((1. / (2. - (2. * u2))), (1.0 / (distribution_index + 1)));
	}

	if (isnan(b1))
	{
		ss.str("");
		ss << "sbx error: denormal beta1 generated: v1:" << v1 << ", v2:" << v2 << ", rand: "<< rands[0] << ", p1:" << p1 << ", p2: " << p2 << ", u1:" << u1 << ", u2: " << u2;
		throw_moea_error(ss.str());
	}
	if (isnan(b2))
	{
		ss.str("");
		ss << "sbx error: denormal beta1 generated: v1:" << v1 << ", v2:" << v2 << ", rand: " << rands[1] << ", p1:" << p1 << ", p2: " << p2 << ", u1:" << u1 << ", u2: " << u2;
		throw_moea_error(ss.str());
	}

	return pair<double, double>(b1, b2);

}

void MOEA::gauss_mutation_ip(ParameterEnsemble& _dp)
{
	/* 
	Information Sciences, Vol. 133/3-4, pp. 229-247 (2001.5)
	Search Space Boundary Extension Method inReal-Coded Genetic Algorithms
	Shigeyoshi Tsutsui* and David E. Goldberg**
	*/

	double org_mutation_probability = pest_scenario.get_pestpp_options().get_mou_mutation_probability();
	if (org_mutation_probability < 0.0)
		org_mutation_probability = 1.0 / _dp.shape().second;

	_dp.transform_ip(ParameterEnsemble::transStatus::NUM);
	_dp.update_var_map();
	map<string, int> var_map = _dp.get_var_map();
	bool adaptive_mr = false;
	if (var_map.find(MR_NAME) != var_map.end())
		adaptive_mr = true;
	vector<string> var_names =_dp.get_var_names();
	string vname;
	Parameters lbnd = pest_scenario.get_ctl_parameter_info().get_low_bnd(var_names);
	Parameters ubnd = pest_scenario.get_ctl_parameter_info().get_up_bnd(var_names);
	pest_scenario.get_base_par_tran_seq().ctl2numeric_ip(lbnd);
	pest_scenario.get_base_par_tran_seq().ctl2numeric_ip(ubnd);
	Eigen::VectorXd indiv;
	vector<double> rnds;
	double mut_val;
	Eigen::MatrixXd reals = _dp.get_eigen();
	double mutation_probability = org_mutation_probability;
	for (int i = 0; i < _dp.shape().first; i++)
	{
		indiv = reals.row(i);
		if (adaptive_mr)
		{
			if (uniform_draws(1, 0.0, 1.0, rand_gen)[0] > 0.1)
				mutation_probability = org_mutation_probability;
			else
				mutation_probability = org_mutation_probability;

		}
		for (int var = 0; var < var_names.size(); var++)
		{
			vname = var_names[var];
			rnds = uniform_draws(2, 0.0, 1.0, rand_gen);
			if (rnds[0] <= mutation_probability)
			{
				mut_val = draw_standard_normal(rand_gen) * (ubnd[vname] - lbnd[vname]) / 10.0;
				indiv[var] += mut_val;
			}
		} 
		reals.row(i) = indiv;	
	}
	_dp.set_eigen(reals);
}




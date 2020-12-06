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

using namespace std;

//util functions



ParetoObjectives::ParetoObjectives(Pest& _pest_scenario, FileManager& _file_manager,
	PerformanceLog* _performance_log):
	pest_scenario(_pest_scenario),file_manager(_file_manager),
	performance_log(_performance_log)

{

}

void ParetoObjectives::update_member_struct(ObservationEnsemble& op, ParameterEnsemble& dp)
{
	member_struct.clear();

	//map<string, map<double, string>> obj_struct;
	vector<string> real_names = op.get_real_names();
	Eigen::VectorXd obj_vals;
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


	//map<string, map<string, double>> member_struct;

	for (auto real_name : real_names)
	{
		map<string, double> obj_map;
		for (auto t : temp)
		{
			//obj_map[t.second[real_name]] = t.first;
			obj_map[t.first] = t.second[real_name];
		}

		member_struct[real_name] = obj_map;
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
				pi_res_sim = pi_ptr->get_pi_rec_ptr(obj_name).calc_sim_and_resid(pars);
				//account for dir mult here
				member_struct[real_name][obj_name] = pi_res_sim.first * obj_dir_mult_ptr->at(obj_name);
			}
		}
	}

}


void ParetoObjectives::drop_duplicates(ObservationEnsemble& op, ParameterEnsemble& dp)
{
	performance_log->log_event("checking for duplicate solutions");
	set<string> duplicates;
	//for (auto solution_p : member_struct)
	//{
	//	for (auto solution_q : member_struct)
	//	{
	//		if (solution_p.first == solution_q.first)
	//			continue;
	//		if ((first_equals_second(solution_p.second, solution_q.second)) && (duplicates.find(solution_q.first) == duplicates.end()))
	//		{
	//			duplicates.emplace(solution_p.first);
	//		}
	//	}
	//}
	//if (duplicates.size() > 0)
	//{
	//	stringstream ss;
	//	ss << "WARNING: " << duplicates.size() << " duplicate solutions found, dropping (see rec file for listing)";
	//	performance_log->log_event(ss.str());
	//	cout << ss.str() << endl;
	//	ofstream& frec = file_manager.rec_ofstream();
	//	frec << "WARNING: " << duplicates.size() << " duplicate solutions found:" << endl;
	//	int i = 0;
	//	for (auto d : duplicates)
	//	{
	//		frec << d << " ";
	//		i++;
	//		if (i > 5)
	//		{
	//			frec << endl;
	//			i = 0;
	//		}
	//	}
	//	/*vector<string> d(duplicates.begin(), duplicates.end());
	//	op.drop_rows(d);
	//	dp.drop_rows(d);
	//	performance_log->log_event("updating member struct after dropping duplicates");
	//	update_member_struct(op, dp);*/
	//}

	//duplicates.clear(); 

	//performance_log->log_event("checking for duplicate solutions");
	vector<string> names;
	map<string, double> solution_p,solution_q;
	string name_p,name_q;
	for (auto m : member_struct)
		names.push_back(m.first);

	//for (auto solution_p : member_struct)
	for (int i=0;i<names.size();i++)
	{
		name_p = names[i];
		solution_p = member_struct[names[i]];
		//for (auto solution_q : member_struct)
		for (int j=i+1;j<names.size();j++)
		{
			name_q = names[j];
			solution_q = member_struct[name_q];
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
		vector<string> d(duplicates.begin(), duplicates.end());
		op.drop_rows(d);
		dp.drop_rows(d);
		performance_log->log_event("updating member struct after dropping duplicates");
		update_member_struct(op, dp);
	}


}

pair<vector<string>, vector<string>> ParetoObjectives::pareto_dominance_sort(int generation, ObservationEnsemble& op,
	ParameterEnsemble& dp, Constraints* constraints_ptr, bool report, string sum_tag)
{
	stringstream ss;
	ss << "ParetoObjectives::pareto_dominance_sort() for " << op.shape().first << " population members";
	performance_log->log_event(ss.str());
	performance_log->log_event("preparing fast-lookup containers");


	ofstream& frec = file_manager.rec_ofstream();
	//TODO: check for a single objective and deal appropriately

	//update the member struct container
	update_member_struct(op, dp);



	if (member_struct.size() == 0)
		throw runtime_error("ParetoObjectives error: member_struct is empty");

	//check for and drop duplictes
	drop_duplicates(op, dp);

	

	vector<string> real_names = op.get_real_names();
	vector<string> infeas_ordered;
	map<string, double> infeas;
	bool all_infeas = false;
	//TODO: work out if there are actually any constraints (instead of all objectives)
	

	if (constraints_ptr)
	{
		bool check_constraints = false;
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
		if (check_constraints)
		{
			performance_log->log_event("feasible sorting");
			//sort into feasible and infeasible 
			map<string, double> violations;
			double vsum;
			Observations obs = pest_scenario.get_ctl_observations();
			Parameters pars = pest_scenario.get_ctl_parameters();
			vector<string> onames = op.get_var_names(), pnames = dp.get_var_names();
			set<string> obs_obj_set(obs_obj_names_ptr->begin(), obs_obj_names_ptr->end());
			set<string> pi_obj_set(pi_obj_names_ptr->begin(), pi_obj_names_ptr->end());
			map<string, map<string, double>> feas_member_struct;
			for (auto real_name : real_names)
			{
				vsum = 0.0;
				obs.update_without_clear(onames, op.get_real_vector(real_name));
				//TODO: add a constraint tol ++ arg and use it here
				// the 'false' arg is to not apply risk shifting to the satisfaction calcs since
				// 'op' has already been shifted
				violations = constraints_ptr->get_unsatified_obs_constraints(obs, 0.0, false);
				for (auto v : violations)
				{
					if (obs_obj_set.find(v.first) == obs_obj_set.end())
						vsum += v.second;
				}
				pars.update_without_clear(pnames, dp.get_real_vector(real_name));
				violations = constraints_ptr->get_unsatified_pi_constraints(pars, 0.0);
				for (auto v : violations)
				{
					if (pi_obj_set.find(v.first) == pi_obj_set.end())
						vsum += v.second;
				}
				if (vsum > 0.0)
				{
					infeas[real_name] = vsum;
				}
				else
					feas_member_struct[real_name] = member_struct[real_name];
			}
			if (feas_member_struct.size() == 0)
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

			//sort the infeasible members by violation
			vector <pair<string, double>> infeas_vec;
			for (auto inf : infeas)
				infeas_vec.push_back(inf);

			std::sort(infeas_vec.begin(), infeas_vec.end(),
				compFunctor);

			for (auto inf : infeas_vec)
				infeas_ordered.push_back(inf.first);
		}
	}

	performance_log->log_event("pareto sorting");
	map<int, vector<string>> front_map = sort_members_by_dominance_into_fronts(member_struct);
	

	if (obs_obj_names_ptr->size() + pi_obj_names_ptr->size() > 1)
	{
		frec << "...pareto dominance sort yielded " << front_map.size() << " domination fronts" << endl;
		for (auto front : front_map)
		{
			if (front.second.size() == 0)
			{
				ss.str("");
				ss << "ParetoObjectives::pareto_dominance_sort() error: front " << front.first << " has no members";
				performance_log->log_event(ss.str());
				throw runtime_error(ss.str());
			}
			frec << front.second.size() << " in the front " << front.first << endl;
		}
	}
	
	vector<string> nondom_crowd_ordered,dom_crowd_ordered;
	vector<string> crowd_ordered_front;
	map<string, double> crowd_map;
	map<string, int> member_front_map;
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
		}

		else
		{
			crowd_ordered_front = sort_members_by_crowding_distance(front.second, crowd_map);
		}

		if (front.first == 1)
			for (auto front_member : crowd_ordered_front)
				nondom_crowd_ordered.push_back(front_member);
		else
			for (auto front_member : crowd_ordered_front)
				dom_crowd_ordered.push_back(front_member);
	}

	//now add the infeasible members
	//if there is atleast one feasible nondom solution, then add the infeasible ones to dom solutions
	if (!all_infeas)
	{
		if (nondom_crowd_ordered.size() > 0)
			for (auto inf : infeas_ordered)
				dom_crowd_ordered.push_back(inf);
		else
			for (auto inf : infeas_ordered)
				nondom_crowd_ordered.push_back(inf);
	}

	if (op.shape().first != nondom_crowd_ordered.size() + dom_crowd_ordered.size())
	{
		ss.str("");
		ss << "ParetoObjectives::pareto_dominance_sort() internal error: final sorted population size: " << 
			nondom_crowd_ordered.size() + dom_crowd_ordered.size() << " != initial population size: " << op.shape().first;
		cout << ss.str();
		throw runtime_error(ss.str());
	}

	//TODO: some kind of reporting here.  Probably to the rec and maybe a csv file too...
	if (sum_tag.size() > 0)
	{

		write_pareto_summary(sum_tag, generation, real_names, member_front_map, crowd_map, infeas);
	}
		
	return pair<vector<string>, vector<string>>(nondom_crowd_ordered, dom_crowd_ordered);

}

void ParetoObjectives::write_pareto_summary(string& sum_tag, int generation, vector<string>& member_names, 
	map<string,int>& front_map, map<string,double>& crowd_map, map<string,double>& infeas_map)
{
	ofstream& sum = file_manager.get_ofstream(sum_tag);
	for (auto& member : member_names)
	{
		
		sum << generation << "," << member;
		for (auto obj : *obs_obj_names_ptr)
		{
			sum << "," << obj_dir_mult_ptr->at(obj) * member_struct[member][obj];
		}
		for (auto obj : *pi_obj_names_ptr)
		{
			sum << "," << obj_dir_mult_ptr->at(obj) * member_struct[member][obj];
		}
		sum << "," << front_map[member];
		sum << "," << crowd_map[member];
		if (infeas_map.find(member) != infeas_map.end())
		{
			sum << "," << 0 << "," << infeas_map[member];
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
	sum << ",front,crowding_distance,is_feasible,feasible_distance" << endl;

}

map<string, double> ParetoObjectives::get_crowding_distance(ObservationEnsemble& op, ParameterEnsemble& dp)
{
	//this updates the complicated map-based structure that stores the member names: obj_names:value nested pairs
	update_member_struct(op, dp);
	//just reuse the same routine used for the pareto dominance sort...
	return get_crowding_distance();
}

map<string, double> ParetoObjectives::get_hypervolume(ObservationEnsemble& op, ParameterEnsemble& dp)
{
	//this updates the complicated map-based structure that stores the member names: obj_names:value nested pairs
	update_member_struct(op, dp);

	map<string, double> hypervol_map;

	//do something cool here...
	


	return hypervol_map;

}

map<string, double> ParetoObjectives::get_crowding_distance()
{
	vector<string> members;
	for (auto m : member_struct)
		members.push_back(m.first);
	return get_crowding_distance(members);
}

map<string, double> ParetoObjectives::get_crowding_distance(vector<string>& members)
{
	
	map<string, map<string, double>> obj_member_map;
	map<string, double> crowd_distance_map;
	string m = members[0];
	vector<string> obj_names;
	for (auto obj_map : member_struct[m])
	{
		obj_member_map[obj_map.first] = map<string, double>();
		obj_names.push_back(obj_map.first);
	}

	for (auto member : members)
	{
		crowd_distance_map[member] = 0.0;
		for (auto obj_map : member_struct[member])
			obj_member_map[obj_map.first][member] = obj_map.second;

	}

	//map<double,string>::iterator start, end;
	map<string, double> omap;
	double obj_range;
	typedef std::set<std::pair<std::string, double>, Comparator> crowdset;
	for (auto obj_map : obj_member_map)
	{
		omap = obj_map.second;
		//note: for members with identical distances, only the first one gets into the 
		//sorted set but this is ok since we initialized the distance map with zeros
		//for all members, so it works out...
		crowdset crowd_sorted(omap.begin(), omap.end(), compFunctor);

		crowdset::iterator start = crowd_sorted.begin(), last = prev(crowd_sorted.end(), 1);

		obj_range = last->second - start->second;

		//the obj extrema - makes sure they are retained 
		crowd_distance_map[start->first] = 1.0e+30;
		crowd_distance_map[last->first] = 1.0e+30;
		if (crowd_sorted.size() == 3)
		{
			crowdset::iterator it = start;
			next(it, 1);
			crowd_distance_map[it->first] = crowd_distance_map[it->first] + ((last->second - start->second) / obj_range);

		}
		else if (crowd_sorted.size() > 3)
		{
			//need iterators to start and stop one off from the edges
			start = next(crowd_sorted.begin(), 1);
			last = prev(crowd_sorted.end(), 2);

			crowdset::iterator it = start;

			crowdset::iterator inext, iprev;
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


vector<string> ParetoObjectives::sort_members_by_crowding_distance(vector<string>& members, map<string,double>& crowd_map)
{

	map<string, double> crowd_distance_map = get_crowding_distance(members);
	

	vector <pair<string, double>> cs_vec;
	for (auto cd : crowd_distance_map)
	{
		cs_vec.push_back(cd);
		crowd_map[cd.first] = cd.second;
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

map<int,vector<string>> ParetoObjectives::sort_members_by_dominance_into_fronts(map<string, map<string, double>>& member_struct)
{
	//following fast non-dom alg in Deb
	performance_log->log_event("starting 'fast non-dom sort");
	//map<string,map<string,double>> Sp, F1;
	map<string, vector<string>> solutions_dominated_map;
	int domination_counter;
	map<string, int> num_dominating_map;
	vector<string> solutions_dominated, first_front;
	performance_log->log_event("finding first front");
	for (auto solution_p : member_struct)
	{
		domination_counter = 0;
		solutions_dominated.clear();
		for (auto solution_q : member_struct)
		{
			if (solution_p.first == solution_q.first) //string compare real name
				continue;

			//if the solutions are identical...
			if (first_equals_second(solution_p.second, solution_q.second))
			{
				throw runtime_error("ParetoObjectives::sort_members_by_dominance_into_fronts(): solution '" + solution_p.first + "' and '" + solution_q.first + "' are identical");
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
		//solution_p is in the first front
		if (domination_counter == 0)
		{
			first_front.push_back(solution_p.first);
		}
		num_dominating_map[solution_p.first] = domination_counter;
		solutions_dominated_map[solution_p.first] = solutions_dominated;

	}
	performance_log->log_event("sorting remaining fronts");
	int i = 1;
	int nq;
    vector<string> q_front;
	map<int, vector<string>> front_map;
	front_map[1] = first_front;
	vector<string> front = first_front;

	int num_front_solutions = front.size();

	while (true)
	{
		
		q_front.clear();
		for (auto solution_p : front)
		{
			solutions_dominated = solutions_dominated_map[solution_p];
			for (auto solution_q : solutions_dominated)
			{
				num_dominating_map[solution_q]--;
				if (num_dominating_map[solution_q] <= 0)
					q_front.push_back(solution_q);
			}
		}
		if (q_front.size() == 0)
			break;
		i++;
		front_map[i] = q_front;
		
		front = q_front;

		num_front_solutions += front.size();
	}
	
	if (num_front_solutions != member_struct.size())
	{
		stringstream ss;
		ss << "ERROR: ParetoObjectives::sort_members_by_dominance_into_fronts(): number of solutions in fronts (";
		ss << num_front_solutions << ") != member_stuct.size() (" << member_struct.size() << endl;
		file_manager.rec_ofstream() << ss.str();
		cout << ss.str();
		throw runtime_error(ss.str());
	}

	return front_map;
}

bool ParetoObjectives::first_equals_second(map<string, double>& first, map<string, double>& second)
{
	for (auto f : first)
	{
		if (f.second != second[f.first])
			return false;
	}
	return true;
}

bool ParetoObjectives::first_dominates_second(map<string,double>& first, map<string,double>& second)
{
	for (auto f: first)
	{
		if (f.second > second[f.first])
			return false;
	}
	return true;
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
	op.set_pest_scenario(&pest_scenario);
	dp_archive.set_rand_gen(&rand_gen);
	dp_archive.set_pest_scenario(&pest_scenario);
	op_archive.set_rand_gen(&rand_gen);
	op_archive.set_pest_scenario(&pest_scenario);
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
	cout << endl << "   ************   " << endl << "    MOEAerror: " << message << endl << endl;
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

map<string, map<string, double>> MOEA::obj_func_report(ParameterEnsemble& _dp, ObservationEnsemble& _op)
{
	map<string, map<string, double>> summary = get_obj_func_summary_stats(_dp, _op);
	stringstream frec;
	//frec << endl << "  ---  Objective Function Summary  ---  " << endl;

	int max_len = get_max_len_obj_name();
	string dir;
	frec << left << setw(max_len) << "objective function" << right << setw(10) << "direction" << setw(10) << "mean" << setw(20) << "standard devation" << setw(12) << "min" << setw(12) << "max" << endl;
	for (auto obs_obj : obs_obj_names)
	{

		frec << left << setw(max_len) << obs_obj;
		dir = "minimize";
		if (obj_dir_mult[obs_obj] == -1)
			dir = "maximize";
		frec << right << setw(10) << dir;
		frec << right << setw(10) << summary[obs_obj]["mean"];
		frec << setw(20) << summary[obs_obj]["std"];
		frec << setw(12) << summary[obs_obj]["min"];
		frec << setw(12) << summary[obs_obj]["max"] << endl;
	}

	
	for (auto pi_obj : pi_obj_names)
	{
		frec << left << setw(max_len) << pi_obj;
		dir = "minimize";
		if (obj_dir_mult[pi_obj] == -1)
			dir = "maximize";
		frec << right << setw(10) << dir;
		frec << right << setw(10) << summary[pi_obj]["mean"];
		frec << setw(20) << summary[pi_obj]["std"];
		frec << setw(12) << summary[pi_obj]["min"];
		frec << setw(12) << summary[pi_obj]["max"] << endl;
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
	ss << left << setw(max_len) << "objective function" << right << setw(11) << "mean change";
	ss << setw(11) << "% change";
	ss << setw(11) << "max change" << setw(11) << "% change";
	ss << setw(11) << "min change" << setw(11) << "% change" << endl;

	vector<string> tags{ "mean","max","min" };
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
			
			ss << right << setw(11) << change;
			ss << setw(11) << percent_change;
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
			
			ss << right << setw(11) << change;
			ss << setw(11) << percent_change;
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

	pair<map<string, double>, map<string, double>> mm = _op.get_moment_maps();
	_op.update_var_map();
	map<string, int> var_map = _op.get_var_map();
	map<string, double> sum;
	map<string, map<string, double>> summary_stats;
	for (auto obs_obj : obs_obj_names)
	{
		sum.clear();
		sum["mean"] = mm.first[obs_obj];
		sum["std"] = mm.second[obs_obj];
		sum["min"] = _op.get_eigen_ptr()->col(var_map[obs_obj]).minCoeff();
		sum["max"] = _op.get_eigen_ptr()->col(var_map[obs_obj]).maxCoeff();
		summary_stats[obs_obj] = sum;
	}

	mm = _dp.get_moment_maps();
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
			sim_res = prior_info_ptr->get_pi_rec_ptr(pi_obj).calc_sim_and_resid(pars);
			pi_vals[pi_obj].push_back(sim_res.first);
		}
	}

	for (auto pi_obj : pi_obj_names)
	{
		vec = stlvec_2_eigenvec(pi_vals[pi_obj]);
		sum.clear();
		sum["mean"] = vec.mean();
		sum["std"] = sqrt((vec.array() - vec.mean()).pow(2).sum() / (vec.size() - 1));
		sum["min"] = vec.minCoeff();
		sum["max"] = vec.maxCoeff();
		summary_stats[pi_obj] = sum;
	}
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
		ss << "population_size < " << error_min_members << ", this is redic, increaing to " << warn_min_members;
		warnings.push_back(ss.str());
		ppo->set_ies_num_reals(warn_min_members);
	}
	if ((ppo->get_mou_population_size() < warn_min_members) && (population_dv_file.size() == 0))
	{
		ss.str("");
		ss << "mou_population_size < " << warn_min_members << ", this is prob too few";
		warnings.push_back(ss.str());
	}
	if (ppo->get_ies_reg_factor() < 0.0)
		errors.push_back("ies_reg_factor < 0.0 - WRONG!");
	
	
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


void MOEA::update_archive(ObservationEnsemble& _op, ParameterEnsemble& _dp)
{
	message(2, "updating archive");
	stringstream ss;
	if (op_archive.shape().first != dp_archive.shape().first)
	{
		ss.str("");
		ss << "MOEA::update_archive(): op_archive members " << op_archive.shape().first << " != dp_archive members " << dp_archive.shape().first;
		throw_moea_error(ss.str());
	}

	//check that members of _op arent in the archive already
	vector<string> keep, temp = op.get_real_names();
	set<string> archive_members(temp.begin(), temp.end());
	for (auto& member : _op.get_real_names())
	{
		if (archive_members.find(member) == archive_members.end())
			keep.push_back(member);
	}
	if (keep.size() == 0)
	{
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
	DomPair dompair = objectives.pareto_dominance_sort(iter,op_archive, dp_archive,&constraints,true,ARC_SUM_TAG);
	
	ss.str("");
	ss << "resizing archive from " << op_archive.shape().first << " to " << dompair.first.size() << " current non-dominated solutions";
	message(2, ss.str());
	op_archive.keep_rows(dompair.first);
	dp_archive.keep_rows(dompair.first);

	if (op_archive.shape().first > archive_size)
	{
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

	save_populations(dp_archive, op_archive, "archive");
}


void MOEA::queue_chance_runs(ParameterEnsemble& _dp)
{
	/* queue up chance-related runs using the class attributes dp and op*/
	stringstream ss;
	if (constraints.should_update_chance(iter))
	{
		dp.transform_ip(ParameterEnsemble::transStatus::NUM);
		Parameters pars = pest_scenario.get_ctl_parameters();
		pest_scenario.get_base_par_tran_seq().ctl2numeric_ip(pars);
		Observations obs = pest_scenario.get_ctl_observations();
		//if this is the first iter and no restart
		
		if (chancepoints == chancePoints::SINGLE)
		{
			pair<Parameters, Observations> po_pair = get_optimal_solution(_dp, op);
			pest_scenario.get_base_par_tran_seq().numeric2ctl_ip(pars);
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
	//queue up any chance related runs
	if (allow_chance)
		queue_chance_runs(_dp);

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
	//failed_real_indices.push_back(0);

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

ObservationEnsemble MOEA::get_chance_shifted_op(ParameterEnsemble& _dp, ObservationEnsemble& _op)
{
	if (risk_obj)
		return constraints.get_chance_shifted_constraints(_dp, _op, RISK_NAME);
	else
		return constraints.get_chance_shifted_constraints(_dp, _op);
}

void MOEA::finalize()
{

}

void MOEA::initialize()
{
	stringstream ss;
	message(0, "initializing MOEA process");
	
	pp_args = pest_scenario.get_pestpp_options().get_passed_args();

	act_obs_names = pest_scenario.get_ctl_ordered_nz_obs_names();
	act_par_names = pest_scenario.get_ctl_ordered_adj_par_names();

	//define these here to make sure the iter loop later behaves
	iter = 0;
	member_count = 0;

	warn_min_members = 20;
	error_min_members = 4;
	

	
	//set some defaults
	PestppOptions* ppo = pest_scenario.get_pestpp_options_ptr();

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

	message(1, "max run fail: ", ppo->get_max_run_fail());


	//some risk-based stuff here
	string chance_points = ppo->get_opt_chance_points();
	if (chance_points == "ALL")
	{
		//evaluate the chance constraints at every individual, very costly, but most robust
		//throw_moea_error("'mou_chance_points' == 'all' not implemented");
		chancepoints = chancePoints::ALL;
	}
	
	else if (chance_points == "SINGLE")
	{
		//evaluate the chance constraints only at the population member nearest the optimal tradeoff.
		//much cheaper, but assumes linear coupling
		chancepoints = chancePoints::SINGLE;
	}
	else
	{
		ss.str("");
		ss << "unrecognized 'mou_chance_points' value :" << chance_points << ", should be 'all' or 'single'";
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
			if (pest_scenario.get_prior_info().get_pi_rec_ptr(oname).get_weight() == 0.0)
				continue;
			sense = Constraints::get_sense_from_group_name(pest_scenario.get_prior_info().get_pi_rec_ptr(oname).get_group());
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
		vector<string> missing,keep_obs, keep_pi,err_sense;
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
			}
			else
			{
				sense = Constraints::get_sense_from_group_name(pest_scenario.get_prior_info().get_pi_rec_ptr(obj_name).get_group());
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
			ss << "WARNING: the following mou_objectives were not found in the non-zero-weighted observations or prior info eqs: ";
			for (auto m : missing)
				ss << m << ",";
			message(1, ss.str());

		}
		obs_obj_names = keep_obs;
		pi_obj_names = keep_pi;
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
	
	if (obs_obj_names.size() + pi_obj_names.size() > 5)
		message(1, "WARNING: more than 5 objectives, this is pushing the limits!");


	if (risk_obj)
	{
		/*file_manager.rec_ofstream() << ss.str();
		cout << ss.str();
		dv_names.push_back(RISK_NAME);
		pi_obj_names.push_back(RISK_NAME);
		PIAtom pia(RISK_NAME, false, 1.0);
		PriorInformationRec pir(0.0, 1.0, "greater_than", vector<PIAtom>{pia});
		pest_scenario.get_prior_info_ptr()->AddRecord(RISK_NAME, &pir);*/
		if (find(dv_names.begin(), dv_names.end(), RISK_NAME) == dv_names.end())
			throw_moea_error(RISK_NAME + " not found in decision variable names");
		if (find(pi_obj_names.begin(), pi_obj_names.end(), RISK_NAME) == pi_obj_names.end())
			throw_moea_error(RISK_NAME + " not found in prior information objective names");
		ss.str("");
		ss << "'mou_risk_objective' is true, using " << RISK_NAME << " decision variable as risk in chance calcs" << endl;
		file_manager.rec_ofstream() << ss.str();
		cout << ss.str();

		//reset bounds of the risk parameter
		double b = pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(RISK_NAME)->lbnd;
		pest_scenario.get_ctl_parameter_info_ptr_4_mod()->get_parameter_rec_ptr_4_mod(RISK_NAME)->lbnd = max(b,0.01);
		b = pest_scenario.get_ctl_parameter_info().get_parameter_rec_ptr(RISK_NAME)->ubnd;
		pest_scenario.get_ctl_parameter_info_ptr_4_mod()->get_parameter_rec_ptr_4_mod(RISK_NAME)->ubnd = min(b, 0.99);
		//set this just to make sure everuything gets initialized right
		pest_scenario.get_pestpp_options_ptr()->set_opt_risk(0.95);


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
		_pe.append("BASE", pars);
		string par_csv = file_manager.get_base_filename() + ".par.csv";
		//message(1, "saving parameter values to ", par_csv);
		//_pe.to_csv(par_csv);
		ParameterEnsemble pe_base = _pe;
		pe_base.reorder(vector<string>(), act_par_names);
		ObservationEnsemble _oe(&pest_scenario, &rand_gen);
		_oe.reserve(vector<string>(), pest_scenario.get_ctl_ordered_obs_names());
		_oe.append("BASE", pest_scenario.get_ctl_observations());
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
		save_base_real_par_rei(pest_scenario, _pe, _oe, output_file_writer, file_manager, -1);

		message(0, "control file parameter objective function summary: ");
		obj_func_report(_pe, _oe);

		vector<string> names = _oe.get_var_names();
		Observations obs(names, _oe.get_real_vector(BASE_REAL_NAME));
		names = _pe.get_var_names();
		pars.update(names, eigenvec_2_stlvec(_pe.get_real_vector(BASE_REAL_NAME)));
		
		constraints.mou_report(0, pars, obs, obs_obj_names, pi_obj_names);
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
		}
		else if (token == "SBX")
		{
			gen_types.push_back(MouGenType::SBX);
		}
		else
		{
			throw_moea_error("unrecognized generator type '" + token + "', should be in {'DE','SBX'}");
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
	

	int num_members = pest_scenario.get_pestpp_options().get_mou_population_size();
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


	//we are restarting
	if (population_obs_restart_file.size() > 0)
	{
		if (constraints.get_use_chance())
		{
			//this can be done, but we need to make sure the appropriate chance restart
			//args were supplied: base_jacobian or obs_stack
			throw_moea_error("chance constraints not yet supported with restart");
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
				throw_moea_error("too few members to continue");
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
		ss << file_manager.get_base_filename() << ".0." << dv_pop_file_tag << ".csv";
		dp.to_csv(ss.str());
		message(1, "saved initial dv population to ", ss.str());
		performance_log->log_event("running initial population");
		message(1, "running initial population of size", dp.shape().first);
	
		vector<int> failed = run_population(dp, op, true);
		if (dp.shape().first == 0)
			throw_moea_error(string("all members failed during initial population evaluation"));
		
		dp.transform_ip(ParameterEnsemble::transStatus::NUM);
	}
	ss.str("");
	ss << file_manager.get_base_filename() << ".0." << obs_pop_file_tag << ".csv";
	op.to_csv(ss.str());
	message(1, "saved observation population to ", ss.str());

	message(0, "initial population objective function summary:");
	previous_obj_summary = obj_func_report(dp, op);

	if (constraints.get_use_chance())
	{
		ObservationEnsemble shifted_op = get_chance_shifted_op(dp, op);
		ss.str("");
		ss << file_manager.get_base_filename() << ".0." << obs_pop_file_tag << ".chance.csv";
		shifted_op.to_csv(ss.str());
		message(1, "saved chance-shifted observation population to ", ss.str());

		op = shifted_op;
		message(0, "chance-shifted initial population objective function summary");
		previous_obj_summary = obj_func_report(dp, op);
	}

	//save the initial dv population again in case runs failed or members were dropped as part of restart
	ss.str("");
	ss << file_manager.get_base_filename() << ".0." << dv_pop_file_tag << ".csv";
	dp.to_csv(ss.str());
	message(1, "saved initial dv population to ", ss.str());

	//TODO: think about a bad phi (or phis) for MOEA
	
	if (op.shape().first <= error_min_members)
	{
		message(0, "too few population members:", op.shape().first);
		message(1, "need at least ", error_min_members);
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
	objectives.set_pointers(obs_obj_names, pi_obj_names, obj_dir_mult);
	DomPair dompair = objectives.pareto_dominance_sort(iter, op, dp,&constraints,true,POP_SUM_TAG);
	
	//initialize op and dp archives
	op_archive = ObservationEnsemble(&pest_scenario, &rand_gen, 
		op.get_eigen(dompair.first, vector<string>()),dompair.first,op.get_var_names());
	
	dp_archive = ParameterEnsemble(&pest_scenario, &rand_gen,
		dp.get_eigen(dompair.first, vector<string>()), dompair.first, dp.get_var_names());
	ss.str("");
	ss << "initialized archives with " << dompair.first.size() << " nondominated members";
	message(2, ss.str());
	archive_size = ppo->get_mou_max_archive_size();

	//ad hoc archive update test
	/*vector<string> temp = op.get_real_names();
	for (auto& m : temp)
		m = "XXX" + m;
	op.set_real_names(temp);
	op.set_eigen(*op.get_eigen_ptr() * 2.0);
	dp.set_real_names(temp);
	update_archive(op, dp);*/
	if (constraints.get_use_chance())
	{
		ofstream& f_rec = file_manager.rec_ofstream();
		pair<Parameters, Observations> po_pair = get_optimal_solution(dp, op);
		constraints.presolve_chance_report(iter, po_pair.second, true,"initial chance constraint summary (calculated at optimal/mean decision variable point)");
	}

	constraints.mou_report(0,dp, op, obs_obj_names,pi_obj_names);

	message(0, "initialization complete");
}


pair<Parameters, Observations> MOEA::get_optimal_solution(ParameterEnsemble& _dp, ObservationEnsemble& _oe)
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

		vector<double> t = dp.get_mean_stl_var_vector();
		Eigen::VectorXd dp_mean = stlvec_2_eigenvec(t);
		t.resize(0);
		int idx_min;
		double dist, dist_min = numeric_limits<double>::max();
		for (int i = 0; i < dp.shape().first; i++)
		{
			dist = dp_mean.dot(dp.get_eigen_ptr()->row(i));
			if (dist < dist_min)
			{
				idx_min = i;
				dist_min = dist;
			}
		}
		string min_member = dp.get_real_names()[idx_min];
		ss.str("");
		ss << "using member " << min_member << " as nearest-to-mean single point" << endl;
		ss << "    with distance of " << dist_min << " from mean of decision variable population";
		message(2, ss.str());

		pars.update_without_clear(dp.get_var_names(), dp.get_real_vector(min_member));
		obs.update_without_clear(op.get_var_names(), op.get_real_vector(min_member));
		pest_scenario.get_base_par_tran_seq().numeric2ctl_ip(pars);
	}
	else
	{
		//calculate the optimal tradeoff point from the current op
		//dont worry about pi-based obj since they arent chance-based
		message(2, "seeking optimal trade-off point for single 'optimal' chance point runs");
		vector<double> obj_extrema;
		Eigen::VectorXd obj_vec; 

		for (auto obj_name : obs_obj_names)
		{
			//if this is a max obj
			if ((obj_dir_mult.find(obj_name) != obj_dir_mult.end()) &&
				(obj_dir_mult[obj_name] == -1.0))
				obj_extrema.push_back(op.get_var_vector(obj_name).maxCoeff());
			else
				obj_extrema.push_back(op.get_var_vector(obj_name).minCoeff());
		}

		Eigen::VectorXd opt_vec = stlvec_2_eigenvec(obj_extrema);

		//find the member nearest the optimal tradeoff
		int opt_idx = -1;
		Eigen::MatrixXd obj_op = op.get_eigen(vector<string>(), obs_obj_names);
		double dist, opt_dist = numeric_limits<double>::max();
		for (int i = 0; i < op.shape().first; i++)
		{
			dist = (opt_vec - obj_op.row(i).transpose()).squaredNorm();
			if (dist < opt_dist)
			{
				opt_idx = i;
				opt_dist = dist;
			}
		}
		string opt_member = op.get_real_names()[opt_idx];
		ss.str("");
		ss << "using member " << opt_member << " as single, 'optimal' point" << endl;
		ss << "   with distance of " << opt_dist << " from optimal trade - off";
		message(2, ss.str());
		pars.update_without_clear(dp.get_var_names(), dp.get_real_vector(opt_member));
		obs.update_without_clear(op.get_var_names(), op.get_real_vector(opt_member));
	}

	return pair<Parameters, Observations>(pars, obs);
}


ParameterEnsemble MOEA::generate_population()
{
	int total_new_members = pest_scenario.get_pestpp_options().get_mou_population_size();
	//add new members for any missing
	total_new_members += (total_new_members - dp.shape().first);
	int new_members_per_gen = int(total_new_members / gen_types.size());
	ParameterEnsemble new_pop(&pest_scenario, &rand_gen);
	new_pop.set_trans_status(ParameterEnsemble::transStatus::NUM);
	for (auto gen_type : gen_types)
	{
		ParameterEnsemble p;
		if (gen_type == MouGenType::DE)
		{
			p = generate_diffevol_population(new_members_per_gen, dp);
		}

		else if (gen_type == MouGenType::SBX)
		{
			p = generate_sbx_population(new_members_per_gen, dp);
		}
		else
			throw_moea_error("unrecognized mou generator");
		if (new_pop.shape().first == 0)
			new_pop = p;
		else
			new_pop.append_other_rows(p);
	}
	return new_pop;
}

void MOEA::iterate_to_solution()
{
	iter = 1;
	int num_members = pest_scenario.get_pestpp_options().get_mou_population_size();
	vector<string> keep;
	stringstream ss;
	map<string, map<string, double>> summary;
	while(iter <= pest_scenario.get_control_info().noptmax)
	{
		message(0, "starting iteration ", iter);

		//generate offspring
		ParameterEnsemble new_dp = generate_population();
		
		//run offspring thru the model while also running risk runs, possibly at many points in dec var space	
		ObservationEnsemble new_op(&pest_scenario, &rand_gen);
		new_op.reserve(new_dp.get_real_names(), op.get_var_names());
		run_population(new_dp, new_op, true);
		
		

		//and risk-shift
			
		//TODO: save new_dp, new_op and new_op_shifted?
		save_populations(new_dp, new_op);
		if (constraints.get_use_chance())
		{
			pair<Parameters,Observations> po = get_optimal_solution(dp, op);
			constraints.presolve_chance_report(iter, po.second,true, "chance constraint summary (calculated at optimal/mean decision variable point)");
			ObservationEnsemble new_op_shifted = get_chance_shifted_op(new_dp, new_op);
			save_populations(new_dp, new_op_shifted,"chance");
			new_op = new_op_shifted;
		}

		//append offspring dp and (risk-shifted) op to make new dp and op containers
		new_dp.append_other_rows(dp);
		new_op.append_other_rows(op);

		//sort according to pareto dominance, crowding distance, and, optionally, feasibility

		message(1, "pareto dominance sorting combined parent-child populations of size ", new_dp.shape().first);

		DomPair dompair = objectives.pareto_dominance_sort(iter, new_op, new_dp,&constraints,true,POP_SUM_TAG);

		//drop shitty members
		//TODO: this is just a cheap hack, prob something more meaningful to be done...
		keep.clear();
		for (auto nondom : dompair.first)
		{
			if (keep.size() >= num_members)
				break;
			keep.push_back(nondom);
			
		}

		if (keep.size() > 0)
		{
			//update the archive of nondom members
			ParameterEnsemble new_dp_nondom = new_dp;
			new_dp_nondom.keep_rows(keep);
			ObservationEnsemble new_op_nondom = new_op;
			new_op_nondom.keep_rows(keep);
			update_archive(new_op_nondom, new_dp_nondom);
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

		
		ss.str("");
		ss << "iteration " << iter << " objective function summary:";
		message(0, ss.str());
		summary = obj_func_report(dp, op);
		ss.str("");
		ss << "iteration " << iter << " objective function change summary:";
		message(0, ss.str());
		obj_func_change_report(summary);
		previous_obj_summary = summary;
		constraints.mou_report(iter, new_dp, new_op, obs_obj_names, pi_obj_names, true);


		iter++;
	}

}


bool MOEA::initialize_dv_population()
{
	stringstream ss;
	int num_members = pest_scenario.get_pestpp_options().get_mou_population_size();
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
			ss << "original_member=" << i;
			real_names.push_back(ss.str());
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
		else if ((par_ext.compare("jcb") == 0) || (par_ext.compare("jco") == 0))
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
			ss << "unrecognized dv population file extension " << par_ext << ", looking for csv, jcb, or jco";
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
	else if ((par_ext.compare("jcb") == 0) || (par_ext.compare("jco") == 0))
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
		ss << "unrecognized obs population restart file extension " << par_ext << ", looking for csv, jcb, or jco";
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

ParameterEnsemble MOEA::generate_diffevol_population(int num_members, ParameterEnsemble& _dp)
{
	message(1, "generating diffevol population of size", num_members);
	vector<int> member_count,working_count, selected, r_int_vec;
	for (int i = 0; i < _dp.shape().first; i++)
		member_count.push_back(i);

	for (int i = 0; i < dv_names.size(); i++)
		r_int_vec.push_back(i);


	Eigen::VectorXd y, x, diff;
	double F = pest_scenario.get_pestpp_options().get_de_f();
	double CR = pest_scenario.get_pestpp_options().get_de_cr();
	double R;
	
	vector<double> cr_vals;
	Eigen::MatrixXd new_reals(num_members, _dp.shape().second);
	new_reals.setZero();
	_dp.transform_ip(ParameterEnsemble::transStatus::NUM);
	vector<string> new_member_names;
	
	//since _dp migth contain both dev vars and pars, we want to 
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
	for (int i = 0; i < num_members; i++)
	{
		working_count = member_count;
		shuffle(working_count.begin(), working_count.end(), rand_gen);
		selected.clear();

		i_last = 0;
		while (true)
		{
			i_last += 1;
			if (i_last > working_count.size())
				throw_moea_error("MOEA::generate_diffevol_population(): internal error seeking random population members for differntial");
			if (i_last == i)
				continue;
			selected.push_back(working_count[i_last]);
			if (selected.size() == 3)
				break;
		}

		//differential vector
		diff = _dp.get_eigen_ptr()->row(selected[0]) + (F * (_dp.get_eigen_ptr()->row(selected[1]) - _dp.get_eigen_ptr()->row(selected[2])));

		//current member if in range, otherwise, select randomly
		if (i < _dp.shape().first)
			x = _dp.get_eigen_ptr()->row(i);
		else
			//this risks "inbreeding" but maybe thats good?!
			x = _dp.get_eigen_ptr()->row(working_count[i_last]);
		//copy to perserve non-dec var values;
		y = x; 
		//random cross over probs - one per decision variable
		cr_vals = uniform_draws(_dp.shape().second, 0.0, 1.0, rand_gen);
		
		//get the R values for this member;
		shuffle(r_int_vec.begin(), r_int_vec.end(), rand_gen);
		R = r_int_vec[0];

		//only change dec vars
		for(int idv=0;idv<dv_names.size();idv++)
		{
			ii = var_map[dv_names[idv]];
			if ((cr_vals[ii] < CR) || (ii == R))
			{
				y[ii] = diff[ii];
			}
		}
		new_name = get_new_member_name("diffevol");
		new_member_names.push_back(new_name);
		lin << new_name;
		for (auto idx : selected)
			lin << "," << real_names[idx];
		lin << endl;
		new_reals.row(i) = y;
	}

	ParameterEnsemble new_dp(&pest_scenario, &rand_gen, new_reals, new_member_names, _dp.get_var_names());
	
	new_dp.set_trans_status(ParameterEnsemble::transStatus::NUM);
	new_dp.enforce_limits(performance_log, false);

	return new_dp;
}

ParameterEnsemble MOEA::generate_sbx_population(int num_members, ParameterEnsemble& _dp)
{
	message(1, "generating SBX population of size", num_members);

	_dp.transform_ip(ParameterEnsemble::transStatus::NUM);

	vector<int> member_count, working_count, selected, r_int_vec;
	vector<double> rnds;
	for (int i = 0; i < _dp.shape().first; i++)
		member_count.push_back(i);

	for (int i = 0; i < dv_names.size(); i++)
		r_int_vec.push_back(i);

	//TODO: move the operators to classes
	//TODO: add algorithm parameters to pp args
	// crossover
	double crossover_probability = 0.9;
	double crossover_distribution_index = 10.0;
	int i_member = 0;
	int p1_idx, p2_idx;
	Eigen::MatrixXd new_reals(num_members, _dp.shape().second);
	new_reals.setZero();
	pair<Eigen::VectorXd, Eigen::VectorXd> children;
	vector<string> new_names;
	ofstream& lin = file_manager.get_ofstream(lineage_tag);
	vector<string> real_names = _dp.get_real_names();
	string new_name;
	while (i_member < num_members)
	
	{

		//randomly select two parents - this is just a temp routine, something better needed...
		working_count = member_count;//copy member count index to working count
		// sampled with or without replacement? TODO: need to figure this out
		shuffle(working_count.begin(), working_count.end(), rand_gen); //randomly shuffle working count
		//just take the first two since this should change each time thru
		p1_idx = working_count[0];
		p2_idx = working_count[1];

		//generate two children thru cross over
		children = sbx(crossover_probability, crossover_distribution_index, p1_idx, p2_idx);

		//put the two children into the child population
		new_reals.row(i_member) = children.first;
		new_name = get_new_member_name("sbx");
		lin << new_name << "," << real_names[p1_idx] << "," << real_names[p2_idx] << endl;
		new_names.push_back(new_name);
		//cout << i_member << "," << p1_idx << "," << p2_idx << new_names[new_names.size() - 1] << endl;
		i_member++;
		if (i_member >= num_members)
			break;
		new_reals.row(i_member) = children.second;
		new_name = get_new_member_name("sbx");
		lin << new_name << "," << real_names[p1_idx] << "," << real_names[p2_idx] << endl;
		new_names.push_back(new_name);
		//cout << i_member << "," << p1_idx << "," << p2_idx << new_names[new_names.size() -1] << endl;
		i_member++;

	}
		
	ParameterEnsemble tmp_dp(&pest_scenario, &rand_gen, new_reals, new_names, _dp.get_var_names());
	tmp_dp.set_trans_status(ParameterEnsemble::transStatus::NUM);
	//tmp_dp.to_csv("temp_cross.csv");
	//mutation
	double mutation_probability = 1.0 / pest_scenario.get_n_adj_par();
	double mutation_distribution_index = 20.0;
	mutate_ip(mutation_probability, mutation_distribution_index, tmp_dp);

	//tmp_dp.to_csv("temp_mut.csv");


	//TODO: return parameter ensemble for the next generation
	return tmp_dp;
}

void MOEA::save_populations(ParameterEnsemble& dp, ObservationEnsemble& op, string tag)
{
	stringstream ss;
	string fname;
	ss << file_manager.get_base_filename() << "." << iter ;
	if (tag.size() > 0)
	{
		ss << "." << tag;
	}
	ss << "." << dv_pop_file_tag << ".csv";
	fname = ss.str();
	dp.to_csv(fname);
	ss.str("");
	ss << "saved decision variable population of size " << dp.shape().first << " X " << dp.shape().second << " to '" << fname << "'";
	message(1, ss.str());
	
	ss.str("");
	ss << file_manager.get_base_filename() << "." << iter;
	if (tag.size() > 0)
	{
		ss << "." << tag;
	}
	ss << "." << obs_pop_file_tag << ".csv";
	fname = ss.str();
	op.to_csv(fname);
	ss.str("");
	ss << "saved observation population of size " << op.shape().first << " X " << op.shape().second << " to '" << fname << "'";
	message(1, ss.str());

}

string MOEA::get_new_member_name(string tag)
{
	stringstream ss;
	ss << "gen=" << iter << "_member=" << member_count;
	if (tag.size() > 0)
	{
		ss << "_" << tag;
	}
	member_count++;
	return ss.str();
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

	if (rnds[0] <= probability) {
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

void MOEA::mutate_ip(double probability, double eta_m, ParameterEnsemble& temp_dp)
{
	vector<double> rnds;
	double delta1, delta2, mut_pow, deltaq;
	double y, yl, yu, val, xy, initalval;
	// for each decision variable, randomly decide to mutate or not
	// for each individual
	vector<string> var_names = temp_dp.get_var_names();
	string vname;
	Parameters lbnd = pest_scenario.get_ctl_parameter_info().get_low_bnd(var_names);
	Parameters ubnd = pest_scenario.get_ctl_parameter_info().get_up_bnd(var_names);
	pest_scenario.get_base_par_tran_seq().ctl2numeric_ip(lbnd);
	pest_scenario.get_base_par_tran_seq().ctl2numeric_ip(ubnd);
	for (int i = 0; i < temp_dp.shape().first; i++)
	{
		Eigen::VectorXd indiv = temp_dp.get_eigen_ptr()->row(i);
		for (int var = 0; var < var_names.size(); var++)
		{

			vname = var_names[var];
			rnds = uniform_draws(2, 0.0, 1.0, rand_gen);
			if (rnds[0] <= probability)
			{
				// ZQ need to access the value for the current parameter
				y = indiv[var];
				initalval = y;
				yl = lbnd[vname];
				yu = ubnd[vname];
				delta1 = (y - yl) / (yu - yl);
				delta2 = (yu - y) / (yu - yl);
				mut_pow = 1.0 / (eta_m + 1.0);
				if (rnds[1] <= 0.5) {
					xy = 1.0 - delta1;
					val = 2.0 * rnds[1] + (1.0 - 2.0 * rnds[1]) * (pow(xy, (eta_m + 1.0)));
					deltaq = pow(val, mut_pow) - 1.0;
				}
				else {
					xy = 1.0 - delta2;
					val = 2.0 * (1.0 - rnds[1]) + 2.0 * (rnds[1] - 0.5) * (pow(xy, (eta_m + 1.0)));
					deltaq = 1.0 - (pow(val, mut_pow));
				}
				y = y + deltaq * (yu - yl);
				if (y < yl)
					y = yl;
				if (y > yu)
					y = yu;

				if (std::isnan(y)) // y can be nan result from the pow
					indiv[var] = initalval;
				else
					indiv[var] = y;
			}
		} // for
		//get the ctl file parameters but only the ones in temp_dp
		Parameters temp = pest_scenario.get_ctl_parameters().get_subset(var_names.begin(), var_names.end());
		//update the values in temp Parameters with the indiv mutated values
		temp.update_without_clear(var_names, indiv);
		//replave the ith row in temp_dp
		temp_dp.replace(i, temp);
	}
}



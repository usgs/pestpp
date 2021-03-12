#include <random>
#include <map>
#include <iomanip>
#include <mutex>
#include <thread>
#include <math.h>
#include <algorithm>
#include "Ensemble.h"
#include "RestartController.h"
#include "utilities.h"
#include "Ensemble.h"
#include "DataAssimilator.h"
#include "ObjectiveFunc.h"
#include "covariance.h"
#include "RedSVD-h.h"
#include "SVDPackage.h"
#include "DataAssimilator.h"
#include "EnsembleMethodUtils.h"


void DataAssimilator::sanity_checks()
{
	PestppOptions* ppo = pest_scenario.get_pestpp_options_ptr();
	vector<string> errors;
	vector<string> warnings;
	stringstream ss;
	

	if (pest_scenario.get_control_info().pestmode == ControlInfo::PestMode::UNKNOWN)
	{
		warnings.push_back("unrecognized 'pestmode', using 'estimation'");
	}


	if (pest_scenario.get_ctl_ordered_pi_names().size() > 0)
	{
		warnings.push_back("prior information equations not supported in pestpp-da, ignoring...");
	}


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
		throw_em_error(string("sanity_check() found some problems - please review rec file"));
	}
	//cout << endl << endl;
}


void DataAssimilator::da_update(int cycle)
{
	stringstream ss;
	ofstream& frec = file_manager.rec_ofstream();

	//string da_method = da_ctl_params.get_svalue("DA_TYPE");
	//message(0, "Assimilation Method :", da_method);
	int user_noptmax = pest_scenario.get_control_info().noptmax;

	bool accept;
	for (int i = 0; i < pest_scenario.get_control_info().noptmax; i++)
	{
		iter++;
		message(0, "starting solve for iteration:", iter);
		ss << "starting solve for iteration: " << iter;
		performance_log->log_event(ss.str());
		//accept = solve_new_da();
		if (pest_scenario.get_pestpp_options().get_ies_use_mda())
			accept = solve_mda();
		else
			accept = solve_glm();
		report_and_save(cycle);
		ph.update(oe, pe);
		last_best_mean = ph.get_mean(L2PhiHandler::phiType::COMPOSITE);
		last_best_std = ph.get_std(L2PhiHandler::phiType::COMPOSITE);
		ph.report(true);
		ph.write(iter, run_mgr_ptr->get_total_runs());
		if (pest_scenario.get_pestpp_options().get_ies_save_rescov())
			ph.save_residual_cov(oe, iter);
		ss.str("");
		ss << file_manager.get_base_filename() << "." << cycle << "." << iter << ".pcs.csv";
		pcs.summarize(pe, ss.str());
		
		if (accept)
			consec_bad_lambda_cycles = 0;
		else
			consec_bad_lambda_cycles++;

		if (should_terminate())		
			break;
	}
}





void DataAssimilator::eig2csv(string name, Eigen::MatrixXd matrix)
{
	ofstream file(name.c_str());

	for (int i = 0; i < matrix.rows(); i++) {
		for (int j = 0; j < matrix.cols(); j++) {
			string str = std::to_string(matrix(i, j));
			if (j + 1 == matrix.cols()) {
				file << str;
			}
			else {
				file << str << ',';
			}
		}
		file << '\n';
	}
}


void DataAssimilator::finalize()
{

}

map<int, map<string, double>> process_da_obs_cycle_table(Pest& pest_scenario, vector <int>& ncycles_in_tables, ofstream& fout_rec, set<string>& obs_in_tbl)
{
	//process da par cycle table
	string filename = pest_scenario.get_pestpp_options().get_da_obs_cycle_table();
	map<int, map<string, double>> obs_cycle_info;
	if (filename.size() > 0)
	{
		fout_rec << "processing 'DA_OBSEERVATION_CYCLE_TABLE' file " << filename;
		pest_utils::ExternalCtlFile cycle_table(filename);
		cycle_table.read_file(fout_rec);
		vector<string> col_names = cycle_table.get_col_names();
		fout_rec << "...using the first column ('" << col_names[0] << "') as observation names" << endl;
		vector<string> onames = pest_scenario.get_ctl_ordered_obs_names();
		set<string> obs_names(onames.begin(), onames.end());
		onames = cycle_table.get_col_string_vector(col_names[0]);
		vector<string> missing, zeroweight, notneg, tbl_obs_names;
		ObservationInfo* oi = pest_scenario.get_observation_info_ptr();
		for (auto oname : onames)
		{
			pest_utils::upper_ip(oname);
			if (obs_names.find(oname) == obs_names.end())
			{
				missing.push_back(oname);
				continue;
			}
			else if (oi->get_weight(oname) == 0.0)
			{
				zeroweight.push_back(oname);
				continue;
			}
			else if (oi->get_observation_rec_ptr(oname)->cycle != -1)
			{
				notneg.push_back(oname);
				continue;
			}
			else
			{
				tbl_obs_names.push_back(oname);
			}
		}
		if (missing.size() > 0)
		{
			fout_rec << "ERROR: The following observations in DA_OBSERVATION_CYCLE_TABLE are not the control file:" << endl;
			for (auto m : missing)
			{
				fout_rec << m << endl;
			}
			throw runtime_error("DA_OBSERVATION_CYCLE_TABLE contains missing observations, see rec file for listing");
		}
		if (zeroweight.size() > 0)
		{
			fout_rec << "ERROR: The following observations in DA_OBSERVATION_CYCLE_TABLE have zero weight:" << endl;
			for (auto p : zeroweight)
			{
				fout_rec << p << endl;
			}
			throw runtime_error("DA_OBSERVATION_CYCLE_TABLE contains zero-weighted observations, see rec file for listing");
		}
		if (notneg.size() > 0)
		{
			fout_rec << "ERROR: The following OBSERVATIONS in DA_OBSERVATION_CYCLE_TABLE do not have cycle=-1:" << endl;
			for (auto p : notneg)
			{
				fout_rec << p << endl;
			}
			throw runtime_error("DA_OBSERVATION_CYCLE_TABLE contains observations with cycle!=-1, see rec file for listing");
		}
		//process the remaining columns - these should be cycle numbers
		obs_in_tbl.insert(tbl_obs_names.begin(), tbl_obs_names.end());
		string col_name;
		vector<string> cycle_vals;
		int cycle;
		double val;
		bool parse_fail = false;
		for (int i = 1; i < col_names.size(); i++)
		{
			col_name = col_names[i];

			try
			{
				cycle = stoi(col_name);
			}
			catch (...)
			{
				fout_rec << "ERROR: could not parse DA_OBSERVATION_CYCLE_TABLE column '" << col_name << "' to integer" << endl;
				throw runtime_error("ERROR parsing DA_OBSERVATIONS_CYCLE_TABLE column " + col_name + " to integer");
			}
			if (obs_cycle_info.find(cycle) != obs_cycle_info.end())
			{
				throw runtime_error("ERROR: DA_OBSERVATION_CYCLE_TABLE cycle column '" + col_name + "' listed more than once");
			}
			cycle_vals = cycle_table.get_col_string_vector(col_name);
			map<string, double> cycle_map;
			for (int i = 0; i < cycle_vals.size(); i++)
			{
				try
				{
					val = stod(cycle_vals[i]);
				}
				catch (...)
				{
					fout_rec << "WARNING: error parsing '" << cycle_vals[i] << "' for observation " << tbl_obs_names[i] << " in cycle " << cycle << ", continuing..." << endl;
					parse_fail = true;
					continue;
				}
				cycle_map[tbl_obs_names[i]] = val;
				ncycles_in_tables.push_back(cycle);
				ncycles_in_tables.push_back(cycle);

				// make sure cycle numbers are sorted and unique
				sort(ncycles_in_tables.begin(), ncycles_in_tables.end());
				ncycles_in_tables.erase(unique(ncycles_in_tables.begin(),
					ncycles_in_tables.end()),
					ncycles_in_tables.end());
			}
			obs_cycle_info[cycle] = cycle_map;
		}
		if (parse_fail)
		{
			cout << "WARNING: error parsing at least one cycle-based observation value" << endl;
			cout << "         from DA_OBSERVATION_CYCLE_TABLE, see rec file for listing." << endl;
		}


	}
	return obs_cycle_info;


}

map<int, map<string, double>> process_da_weight_cycle_table(Pest& pest_scenario, vector <int>& ncycles_in_tables, ofstream& fout_rec, set<string>& obs_in_tbl)
{
	//process da par cycle table
	string filename = pest_scenario.get_pestpp_options().get_da_weight_cycle_table();
	map<int, map<string, double>> weight_cycle_info;
	if (filename.size() > 0)
	{
		fout_rec << "processing 'DA_WEIGHT_CYCLE_TABLE' file " << filename;
		pest_utils::ExternalCtlFile cycle_table(filename);
		cycle_table.read_file(fout_rec);
		vector<string> col_names = cycle_table.get_col_names();
		fout_rec << "...using the first column ('" << col_names[0] << "') as observation names" << endl;
		vector<string> onames = pest_scenario.get_ctl_ordered_obs_names();
		set<string> obs_names(onames.begin(), onames.end());
		onames = cycle_table.get_col_string_vector(col_names[0]);
		vector<string> missing, zeroweight, notneg, tbl_obs_names;
		ObservationInfo* oi = pest_scenario.get_observation_info_ptr();
		for (auto oname : onames)
		{
			pest_utils::upper_ip(oname);
			if (obs_names.find(oname) == obs_names.end())
			{
				missing.push_back(oname);
				continue;
			}
			else if (oi->get_observation_rec_ptr(oname)->cycle != -1)
			{
				notneg.push_back(oname);
				continue;
			}
			else
			{
				tbl_obs_names.push_back(oname);
			}
		}
		if (missing.size() > 0)
		{
			fout_rec << "ERROR: The following observations in DA_WEIGHT_CYCLE_TABLE are not the control file:" << endl;
			for (auto m : missing)
			{
				fout_rec << m << endl;
			}
			throw runtime_error("DA_WEIGHT_CYCLE_TABLE contains missing observations, see rec file for listing");
		}
		
		if (notneg.size() > 0)
		{
			fout_rec << "ERROR: The following OBSERVATIONS in DA_WEIGHT_CYCLE_TABLE do not have cycle=-1:" << endl;
			for (auto p : notneg)
			{
				fout_rec << p << endl;
			}
			throw runtime_error("DA_WEIGHT_CYCLE_TABLE contains observations with cycle!=-1, see rec file for listing");
		}
		//process the remaining columns - these should be cycle numbers
		obs_in_tbl.insert(tbl_obs_names.begin(), tbl_obs_names.end());
		string col_name;
		vector<string> cycle_vals;
		int cycle;
		double val;
		bool parse_fail = false;
		for (int i = 1; i < col_names.size(); i++)
		{
			col_name = col_names[i];

			try
			{
				cycle = stoi(col_name);
			}
			catch (...)
			{
				fout_rec << "ERROR: could not parse DA_WEIGHT_CYCLE_TABLE column '" << col_name << "' to integer" << endl;
				throw runtime_error("ERROR parsing DA_WEIGHT_CYCLE_TABLE column " + col_name + " to integer");
			}
			if (weight_cycle_info.find(cycle) != weight_cycle_info.end())
			{
				throw runtime_error("ERROR: DA_WEIGHT_CYCLE_TABLE cycle column '" + col_name + "' listed more than once");
			}
			cycle_vals = cycle_table.get_col_string_vector(col_name);
			map<string, double> cycle_map;
			for (int i = 0; i < cycle_vals.size(); i++)
			{
				try
				{
					val = stod(cycle_vals[i]);
				}
				catch (...)
				{
					fout_rec << "WARNING: error parsing weight '" << cycle_vals[i] << "' for observation " << tbl_obs_names[i] << " in cycle " << cycle << ", continuing..." << endl;
					parse_fail = true;
					continue;
				}
				cycle_map[tbl_obs_names[i]] = val;
			}
			weight_cycle_info[cycle] = cycle_map;
			ncycles_in_tables.push_back(cycle);

			// make sure cycle numbers are sorted and unique
			sort(ncycles_in_tables.begin(), ncycles_in_tables.end());
			ncycles_in_tables.erase(unique(ncycles_in_tables.begin(),
								ncycles_in_tables.end()),
								ncycles_in_tables.end());
		}
		if (parse_fail)
		{
			cout << "WARNING: error parsing at least one cycle-based observation value" << endl;
			cout << "         from DA_WEIGHT_CYCLE_TABLE, see rec file for listing." << endl;
		}


	}
	return weight_cycle_info;


}

map<int, map<string, double>> process_da_par_cycle_table(Pest& pest_scenario, vector <int>& ncycles_in_tables, ofstream& fout_rec)
{
	//process da par cycle table
	string filename = pest_scenario.get_pestpp_options().get_da_par_cycle_table();
	map<int, map<string, double>> par_cycle_info;
	if (filename.size() > 0)
	{
		fout_rec << "processing 'DA_PARAMETER_CYCLE_TABLE' file " << filename;
		pest_utils::ExternalCtlFile cycle_table(filename);
		cycle_table.read_file(fout_rec);
		vector<string> col_names = cycle_table.get_col_names();
		fout_rec << "...using the first column ('" << col_names[0] << "') as parameter names" << endl;
		vector<string> pnames = pest_scenario.get_ctl_ordered_par_names();
		set<string> par_names(pnames.begin(), pnames.end());
		pnames = cycle_table.get_col_string_vector(col_names[0]);
		ParameterInfo pi = pest_scenario.get_ctl_parameter_info();
		vector<string> missing, notfixed, notneg, tbl_par_names;

		for (auto pname : pnames)
		{
			pest_utils::upper_ip(pname);
			if (par_names.find(pname) == par_names.end())
			{
				missing.push_back(pname);
				continue;
			}
			else if (pi.get_parameter_rec_ptr(pname)->tranform_type != ParameterRec::TRAN_TYPE::FIXED)
			{
				notfixed.push_back(pname);
				continue;
			}
			else if (pi.get_parameter_rec_ptr(pname)->cycle != -1)
			{
				notneg.push_back(pname);
				continue;
			}
			else
			{
				tbl_par_names.push_back(pname);
			}
		}
		if (missing.size() > 0)
		{
			fout_rec << "ERROR: The following parameters in DA_PARAMETER_CYCLE_TABLE are not the control file:" << endl;
			for (auto p : missing)
			{
				fout_rec << p << endl;
			}
			throw runtime_error("DA_PARAMTER_CYCLE_TABLE contains missing parameters, see rec file for listing");
		}
		if (notfixed.size() > 0)
		{
			fout_rec << "ERROR: The following parameters in DA_PARAMETER_CYCLE_TABLE are not 'fixed':" << endl;
			for (auto p : notfixed)
			{
				fout_rec << p << endl;
			}
			throw runtime_error("DA_PARAMTER_CYCLE_TABLE contains non-fixed parameters, see rec file for listing");
		}
		if (notneg.size() > 0)
		{
			fout_rec << "ERROR: The following parameters in DA_PARAMETER_CYCLE_TABLE do not have cycle=-1:" << endl;
			for (auto p : notneg)
			{
				fout_rec << p << endl;
			}
			throw runtime_error("DA_PARAMTER_CYCLE_TABLE contains parameters with cycle!=-1, see rec file for listing");
		}
		//process the remaining columns - these should be cycle numbers
		string col_name;
		vector<string> cycle_vals;
		int cycle;
		double val;
		bool parse_fail = false;
		for (int i = 1; i < col_names.size(); i++)
		{
			col_name = col_names[i];

			try
			{
				cycle = stoi(col_name);
			}
			catch (...)
			{
				fout_rec << "ERROR: could not parse DA_PARAMETER_CYCLE_TABLE column '" << col_name << "' to integer" << endl;
				throw runtime_error("ERROR parsing DA_PARAMETER CYCLE TABLE column " + col_name + " to integer");
			}
			if (par_cycle_info.find(cycle) != par_cycle_info.end())
			{
				throw runtime_error("ERROR: DA_PARAMETER_CYCLE_TABLE cycle column '" + col_name + "' listed more than once");
			}
			cycle_vals = cycle_table.get_col_string_vector(col_name);
			map<string, double> cycle_map;
			for (int i = 0; i < cycle_vals.size(); i++)
			{
				try
				{
					val = stod(cycle_vals[i]);
				}
				catch (...)
				{
					fout_rec << "WARNING: error parsing '" << cycle_vals[i] << "' for parameter " << tbl_par_names[i] << " in cycle " << cycle << ", continuing..." << endl;
					parse_fail = true;
					continue;
				}
				cycle_map[tbl_par_names[i]] = val;
			}
			par_cycle_info[cycle] = cycle_map;
			ncycles_in_tables.push_back(cycle);

			// make sure cycle numbers are sorted and unique
			sort(ncycles_in_tables.begin(), ncycles_in_tables.end());
			ncycles_in_tables.erase(unique(ncycles_in_tables.begin(), 
											ncycles_in_tables.end()),
											ncycles_in_tables.end());


		}
		if (parse_fail)
		{
			cout << "WARNING: error parsing at least one cycle-based parameter value" << endl;
			cout << "         from DA_PARAMETER_CYCLE_TABLE, see rec file for listing." << endl;
		}


	}
	return par_cycle_info;
}


void write_global_phi_info(int icycle, ofstream& f_phi, DataAssimilator& da, vector<string>& init_real_names)
{
	int iter = da.get_iter();
	f_phi << icycle << "," << iter << "," << da.get_phi_handler().get_mean(L2PhiHandler::phiType::ACTUAL);
	f_phi << "," << da.get_phi_handler().get_std(L2PhiHandler::phiType::ACTUAL);
	f_phi << "," << da.get_phi_handler().get_min(L2PhiHandler::phiType::ACTUAL);
	f_phi << "," << da.get_phi_handler().get_max(L2PhiHandler::phiType::ACTUAL);
	map<string, double> final_actual = da.get_phi_handler().get_phi_map(L2PhiHandler::phiType::ACTUAL);
	for (auto rname : init_real_names)
	{
		if (final_actual.find(rname) != final_actual.end())
			f_phi << "," << final_actual.at(rname);
		else
			f_phi << ",";
	}
	f_phi << endl;
}

void generate_global_ensembles(DataAssimilator& da, ofstream& fout_rec, ParameterEnsemble& curr_pe, ObservationEnsemble& curr_oe)
{
	//generate a parent ensemble which includes all parameters across all cycles
	cout << "...preparing global parameter ensemble for all parameters across all cycles" << endl;
	fout_rec << "...preparing global parameter ensemble for all parameters across all cycles" << endl;
	Pest pest_scenario = da.get_pest_scenario();
	da.initialize_parcov();
	da.initialize_pe(*da.get_parcov_ptr());
	curr_pe = da.get_pe();
	Parameters pars = pest_scenario.get_ctl_parameters();
	ParamTransformSeq pts = pest_scenario.get_base_par_tran_seq();
	if (curr_pe.get_trans_status() != ParameterEnsemble::transStatus::NUM)
		throw runtime_error("parameter ensemble not in numeric transform status");
	pts.ctl2numeric_ip(pars);
	Eigen::VectorXd v = pars.get_data_eigen_vec(pest_scenario.get_ctl_ordered_adj_par_names());
	if (pest_scenario.get_control_info().noptmax == 0)
	{
		curr_pe.reserve(vector<string>{BASE_REAL_NAME}, pest_scenario.get_ctl_ordered_adj_par_names());
		curr_pe.update_real_ip(BASE_REAL_NAME, v);
	}
	else if (pest_scenario.get_pestpp_options().get_ies_include_base())
	{
		cout << "...replacing last realization with 'base' realization in global parameter ensemble" << endl;
		fout_rec << "...replacing last realization with 'base' realization in global parameter ensemble" << endl;
		curr_pe.replace(curr_pe.shape().first - 1, pars, BASE_REAL_NAME);
	}

	cout << "...preparing global observation ensemble for all observations across all cycles" << endl;
	fout_rec << "...preparing global observation ensemble for all observations across all cycles" << endl;
	mt19937 rand_gen = da.get_rand_gen();
	curr_oe.set_pest_scenario(&pest_scenario);
	curr_oe.set_rand_gen(&rand_gen);
	curr_oe.reserve(curr_pe.get_real_names(), pest_scenario.get_ctl_ordered_obs_names());

	try
	{
		curr_pe.check_for_dups();
	}
	catch (const exception& e)
	{
		string message = e.what();
		throw runtime_error("error in parameter ensemble: " + message);
	}
	//todo: add the base realization here so that it is present thru out the cycles below

	cout << "...global parameter ensemble has " << curr_pe.shape().first << " rows and " << curr_pe.shape().second << " columns" << endl;
	fout_rec << "...global parameter ensemble has " << curr_pe.shape().first << " rows and " << curr_pe.shape().second << " columns" << endl;
	curr_pe.to_csv(da.get_file_manager().get_base_filename() + ".global.prior.pe.csv");
	cout << "...global observation ensemble has " << curr_oe.shape().first << " rows and " << curr_oe.shape().second << " columns" << endl;
	fout_rec << "...global observation ensemble has " << curr_oe.shape().first << " rows and " << curr_oe.shape().second << " columns" << endl;

	return;
}

#include <vector>
#include <string>
#include <sstream>
#include "Pest.h"
#include "utilities.h"
#include "eigen_tools.h"
#include "covariance.h"
#include "linear_analysis.h"
#include <iomanip>
#include "OutputFileWriter.h"
#include "ModelRunPP.h"
#include "Ensemble.h"
#include <iterator>
#include "PerformanceLog.h"
#include "EnsembleSmoother.h"

vector<string> get_common(vector<string> v1, vector<string> v2)
{
	vector<string> common;
	for (auto e1 : v1)
	for (auto e1 : v1)
		if (find(v2.begin(), v2.end(), e1) != v2.end()) common.push_back(e1);
	return common;
}

map<string, double> get_obj_comps(string &filename)
{
	ifstream ifile(filename);
	string line;

	vector<string> tokens;
	map<string, double> obj_comps;
	if (!ifile.good()) throw runtime_error("linear_analysis::get_obj_comps() error opening file: " + filename);
	string ext = filename.substr(filename.find_last_of(".") + 1);
	pest_utils::upper_ip(ext);
	if ((ext == "REI") || (ext == "RES"))
	{
		double resid, weight;
		for (;;)
		{
			if (!getline(ifile, line)) throw runtime_error("lienar_analysis::get_obj_comps() error: EOF while looking for 'NAME'");
			pest_utils::upper_ip(line);
			if (line.find("NAME") != string::npos)
			{
				for (;;)
				{
					if (!getline(ifile, line)) break;
					pest_utils::upper_ip(line);
					tokens.clear();
					pest_utils::tokenize(line, tokens);
					pest_utils::convert_ip(tokens[4], resid);
					pest_utils::convert_ip(tokens[5], weight);
					if (obj_comps.find(tokens[1]) != obj_comps.end())obj_comps[tokens[1]] += pow(resid * weight, 2);
					else obj_comps[tokens[1]] = pow(resid * weight, 2);
				}
				break;
			}
		}
	}

	else if (ext == "REC")
	{
		throw runtime_error("linear_analysis::get_obj_comps() .rec not implemented");
	}
	else if (ext == "IOBJ")
	{
		throw runtime_error("linear_analysis::get_obj_comps() .iobj not implemented");
	}

	else throw runtime_error("linear_analysis::get_obj_comps() error: unrecognized file type: " + filename + " must be .rei, or .res");
	ifile.close();
	return obj_comps;

}


map<string, int> get_nnz_group(Pest &pest_scenario)
{
	vector<string> ogrp_names = pest_scenario.get_ctl_ordered_obs_group_names();
	map<string, string> pst_grps = pest_scenario.get_observation_groups();
	ObservationInfo obs_info = pest_scenario.get_ctl_observation_info();
	const ObservationRec* obs_rec;
	map<string, int> grp_nnz;
	for (auto &ogrp : ogrp_names) grp_nnz[ogrp] = 0;

	for (auto &ogrp : pst_grps)
	{
		obs_rec = obs_info.get_observation_rec_ptr(ogrp.first);
		if (obs_rec->weight > 0.0) grp_nnz[ogrp.second]++;
	}
	return grp_nnz;
}

ObservationInfo normalize_weights_by_residual(Pest &pest_scenario, PhiData obj_comps)
{
	ObservationInfo obs_info = pest_scenario.get_ctl_observation_info();
	map<string, string> pst_grps = pest_scenario.get_observation_groups();
	vector<string> ogrp_names = pest_scenario.get_ctl_ordered_obs_group_names();
	map<string, int> grp_nnz = get_nnz_group(pest_scenario);

	const ObservationRec* obs_rec;
	double weight;
	double exp_obj = 0.0;
	if (exp_obj > 0.0)
	{
		if (pest_scenario.get_regul_scheme_ptr())
		{
			cout << " WARNING:  can't use EXPECTED_OBJ option with dynamic" << endl;
			cout << "           regularization - using PHIMACCEPT instead." << endl;
		}
		double tot_obj = obj_comps.meas;
		if (tot_obj > 0.0)
		{
			double mult = exp_obj / tot_obj;
			for (auto &oc : obj_comps.group_phi)
			{
				obj_comps.group_phi[oc.first] = mult * oc.second;
			}
		}
	}
	//if using regularization, we need to check if scaling is needed -> is phimaccept been satisfied
	/*if ((pest_scenario.get_regul_scheme_ptr()->get_use_dynamic_reg()))
	{
		double phimlim = pest_scenario.get_regul_scheme_ptr()->get_phimlim();
		double phimaccept = pest_scenario.get_regul_scheme_ptr()->get_phimaccept();
		double tot_obj = obj_comps.meas;
		if (tot_obj > phimaccept)
		{
			double mult = phimaccept / tot_obj;
			for (auto &oc : obj_comps.group_phi)
				obj_comps.group_phi[oc.first] = mult * oc.second;
		}
		else
		{
			for (auto &oc : obj_comps.group_phi)
				obj_comps.group_phi[oc.first] = 0.0;
		}
	}
*/
	for (auto &ogrp : pst_grps)
	{
		if (obj_comps.group_phi[ogrp.second] <= numeric_limits<double>::min())
			continue;
		obs_rec = obs_info.get_observation_rec_ptr(ogrp.first);
		if (obs_rec->weight > 0.0)
		{
			weight = obs_info.get_observation_rec_ptr(ogrp.first)->weight *
				sqrt(((double)grp_nnz[ogrp.second]) / obj_comps.group_phi[ogrp.second]);
			if (weight <= numeric_limits<double>::min())
				weight = 0.0;
			else if (weight >= numeric_limits<double>::max())
				weight = 1.0e+30;
			obs_info.set_weight(ogrp.first, weight);
		}
	}

	return obs_info;
}

ObservationInfo normalize_weights_by_residual(Pest &pest_scenario, Observations &sim)
{
	ObservationInfo obs_info(pest_scenario.get_ctl_observation_info());
	Observations obs = pest_scenario.get_ctl_observations();

	const ObservationRec* obs_rec;
	double weight,swr,new_weight;

	for (auto &oname : pest_scenario.get_ctl_ordered_nz_obs_names())
	{
		weight = obs_info.get_observation_rec_ptr(oname)->weight;
		swr = pow(((obs[oname] - sim[oname]) * weight), 2);
		new_weight = weight * sqrt(1.0 / swr);
			
		if (new_weight > weight)
			new_weight = weight;
		else if (new_weight <= numeric_limits<double>::min())
			new_weight = 0.0;
		else if (new_weight >= numeric_limits<double>::max())
			new_weight = 1.0e+30;
		obs_info.set_weight(oname, new_weight);
		
	}
	return obs_info;
}

ObservationInfo LinearAnalysis::glm_iter_fosm(ModelRun& optimum_run, OutputFileWriter& output_file_writer, int iter,
	RunManagerAbstract* run_mgr_ptr)
{
	ofstream& fout_rec = file_manager.rec_ofstream();

	//ofstream& pfm = file_manager.get_ofstream("log");
	stringstream ss;
	ss << "starting FOSM uncertainty analyses for iteration " << iter;
	pfm.log_event(ss.str());
	fout_rec << ss.str() << endl << endl;

	/*if (base_jacobian_ptr->get_base_numeric_par_names().size() == 0)
	{
		cout << "WARNING: no parameters in base jacobian, can't calculate uncertainty with FOSM" << endl;
		fout_rec << "WARNING: no parameters in base jacobian, can't calculate uncertainty with FOSM" << endl;
		return 0;
	}*/
	if (jacobian.get_col_names().size() == 0)
	{
		throw runtime_error("LinearAnalysis::glml_iter_fosm() error: no parameters in jacobian");
	}

	//instance of a Mat for the jco
	/*Mat j(base_jacobian_ptr->get_sim_obs_names(), base_jacobian_ptr->get_base_numeric_par_names(),
		base_jacobian_ptr->get_matrix_ptr());*/

		//get a new obs info instance that accounts for residual phi (and expected objection value if passed)
		// and report new weights to the rec file
	fout_rec << endl;
	ObservationInfo reweight;
	Observations sim = optimum_run.get_obs();
	reweight = normalize_weights_by_residual(pest_scenario, sim);
	
	/*fout_rec << "Note: The observation covariance matrix has been constructed from " << endl;
	fout_rec << "      weights listed in the pest control file that have been scaled by " << endl;
	fout_rec << "      by the final residuals to account for " << endl;
	fout_rec << "      the level of measurement noise implied by the original weights so" << endl;
	fout_rec << "      the total objective function is equal to the number of  " << endl;
	fout_rec << "      non-zero weighted observations." << endl;
	fout_rec << endl;*/

	ss.str("");
	if (iter != -999)
		ss << file_manager.get_base_filename() << "." << iter << ".fosm_reweight.rei";
	else
		ss << file_manager.get_base_filename() << ".fosm_reweight.rei";
	string reres_filename = ss.str();
	ofstream reres_of(reres_filename);

	Observations obs = pest_scenario.get_ctl_observations();
	output_file_writer.obs_report(reres_of, obs, sim, reweight);
	//fout_rec << "Scaled observation weights used to form observation noise covariance matrix written to residual file '" << reres_filename << "'" << endl << endl;

	//reset the obscov using the scaled residuals - pass empty prior info names so none are included
	pfm.log_event("loading obscov");
	map<string, double> obs_std = pest_scenario.get_ext_file_double_map("observation data external", "standard_deviation");
	obscov.from_observation_weights(file_manager.rec_ofstream(), obscov.get_col_names(), reweight, vector<string>(), 
		pest_scenario.get_prior_info_ptr(),obs_std);

	//reset the parcov since pars are being frozen and unfrozen iter by iter
	pfm.log_event("loading parcov");
	const string parcov_filename = pest_scenario.get_pestpp_options().get_parcov_filename();
	parcov = Covariance();
	bool parcov_success = false;
	if (parcov_filename.size() > 0)
	{
		try
		{
			parcov.try_from(pest_scenario, file_manager, true);
			parcov_success = true;
		}
		catch (exception &e)
		{
			pfm.log_event("unable to load parcov from file: " + parcov_filename + ", reverting to parameter bounds "+e.what());
		}
	}
	if (!parcov_success)
	{
		try
		{
			parcov.from_parameter_bounds(pest_scenario, file_manager.rec_ofstream());
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::glm_fosm_iter() error setting parcov from parameter bounds:" + string(e.what()));
		}
	}

	//if needed, set the predictive sensitivity vectors
	vector<string> pred_names = pest_scenario.get_pestpp_options().get_prediction_names();

	//if no preds, check for zero-weighted obs to use
	if (pred_names.size() == 0)
	{

		for (auto& oname : pest_scenario.get_ctl_ordered_obs_names())
		{
			if (pest_scenario.get_ctl_observation_info().get_weight(oname) == 0.0)
			{
				pred_names.push_back(oname);
			}
		}
		if (pred_names.size() > 0)
		{
			//cout << "Note: since no forecast/predictions were passed, using " << pred_names.size() << " zero-weighted obs as forecasts" << endl;
			//fout_rec << "Note: since no forecast/predictions were passed, using " << pred_names.size() << " zero-weighted obs as forecasts" << endl;

		}
	}

	//make sure prediction weights are zero
	else
	{
		for (auto& pname : pred_names)
		{
			if (pest_scenario.get_ctl_observation_info().get_weight(pname) != 0.0)
			{
				//cout << endl << "WARNING: prediction: " << pname << " has a non-zero weight" << endl << endl;
				fout_rec << endl << "WARNING: prediction: " << pname << " has a non-zero weight" << endl << endl;
			}
		}
	}
	if (pred_names.size() > 0)
		pfm.log_event("setting predictions");
		set_predictions(pred_names, true);

	//drop all 'regul' obs and equations
	//no longer need to call this since the PI is not being added to the obscov during 
	//reconstruction with scaled weights
	//drop_prior_information(pest_scenario);

	//write the posterior covariance matrix
	ss.str("");
	if (iter != -999)
		ss << file_manager.get_base_filename() << "." << iter + 1 << ".post.cov";
	else
		ss << file_manager.get_base_filename() << ".post.cov";
	string postcov_filename = ss.str();
	pfm.log_event("parameter fosm calcs");
	posterior_parameter_ptr()->to_ascii(postcov_filename);
	fout_rec << "posterior parameter covariance matrix written to file '" + postcov_filename +
		"'" << endl << endl;

	//write a parameter prior and posterior summary to the rec file
	const ParamTransformSeq trans = pest_scenario.get_base_par_tran_seq();
	Parameters pars = pest_scenario.get_ctl_parameters();
	ss.str("");
	if (iter != -999)
		ss << file_manager.get_base_filename() << "." << iter << ".par.usum.csv";
	else
		ss << file_manager.get_base_filename() << ".par.usum.csv";
	string parsum_filename = ss.str();
	write_par_credible_range(fout_rec, parsum_filename, pest_scenario.get_ctl_parameter_info(),
		trans.active_ctl2numeric_cp(pest_scenario.get_ctl_parameters()),
		trans.active_ctl2numeric_cp(optimum_run.get_ctl_pars()),
		pest_scenario.get_ctl_ordered_par_names());
	fout_rec << "the above parameter uncertainty summary was written to file '" + parsum_filename +
		"'" << endl << endl;


	//if predictions were defined, write a prior and posterior summary to the rec file
	if (pred_names.size() > 0)
	{
		pfm.log_event("forecast FOSM calcs");
		map<string, pair<double, double>> init_final_pred_values;
		double ival, fval;
		vector<string> run_mgr_obs_names = run_mgr_ptr->get_obs_name_vec();
		vector<double> run_mgr_obs_vals = run_mgr_ptr->get_init_sim();
		for (auto& pred_name : pred_names)
		{
			fval = optimum_run.get_obs().get_rec(pred_name);
			if (run_mgr_ptr->get_init_sim().size() > 0)
			{
				int idx = std::distance(run_mgr_obs_names.begin(), find(run_mgr_obs_names.begin(),
					run_mgr_obs_names.end(), pred_name));
				ival = run_mgr_obs_vals[idx];
			}
			else
			{
				//cout << "WARNING: initial simulation results not available, falling back to optimum run outputs for prior forecast mean" << endl;
				// << "WARNING: initial simulation results not available for FOSM calcs, falling back to optimum run outputs for prior forecast mean" << endl;
				pfm.log_event("no initial sim results, using current run instead: " + pred_name);
				ival = fval;
			}

			init_final_pred_values[pred_name] = pair<double, double>(ival, fval);
		}

		ss.str("");
		if (iter != -999)
			ss << file_manager.get_base_filename() << "." << iter << ".pred.usum.csv";
		else
			ss << file_manager.get_base_filename() << ".pred.usum.csv";
		string predsum_filename = ss.str();
		write_pred_credible_range(fout_rec, predsum_filename, init_final_pred_values);
		fout_rec << "Note : the above forecast uncertainty summary was written to file '" + predsum_filename +
			"'" << endl << endl;
	}
	return reweight;
}

pair<ParameterEnsemble,map<int,int>> LinearAnalysis::draw_fosm_reals(RunManagerAbstract* run_mgr_ptr, int iter, 
	ModelRun& optimum_run)
{
	set<string> args = pest_scenario.get_pestpp_options().get_passed_args();
	ofstream& fout_rec = file_manager.rec_ofstream();
	map<int, int> run_map;
	Covariance cov = posterior_parameter_matrix();
	//check for missing adjustable pars
	/*vector<string> adj_names = pest_scenario.get_ctl_ordered_adj_par_names();
	set<string> sadj_names(adj_names.begin(), adj_names.end());
	vector<string> missing;
	for (auto name : cov.get_row_names())
	{
		if (sadj_names.find(name) == sadj_names.end())
			missing.push_back(name);
	}
	if (missing.size() > 0)
	{
		for (auto m : missing)
		{
			pest_scenario.
		}
	}*/
	ParameterEnsemble pe(&pest_scenario,rand_gen_ptr);
	
	if (pest_scenario.get_pestpp_options().get_glm_num_reals() > 0)
	{
		pfm.log_event("drawing, saving and queuing FOSM parameter realizations");
		bool binary = pest_scenario.get_pestpp_options().get_save_binary();
		int num_reals = pest_scenario.get_pestpp_options().get_glm_num_reals();
		pe.draw(num_reals, optimum_run.get_ctl_pars(), cov, &pfm, 2, file_manager.rec_ofstream());
		stringstream ss;
		ss.str("");
		if (iter != -999)
			ss << file_manager.get_base_filename() << "." << iter + 1 << ".post.paren";
		else
			ss << file_manager.get_base_filename() << ".post.paren";
		if (binary)
		{
			pe.to_binary(ss.str() + ".jcb");
			cout << "...posterior parameter ensemble saved to " << ss.str() << ".jcb" << endl;
			file_manager.rec_ofstream() << "...posterior parameter ensemble saved to " << ss.str() << ".jcb" << endl;
		}
		else
		{
			pe.to_csv(ss.str() + ".csv");
			cout << "...posterior parameter ensemble saved to " << ss.str() << ".csv" << endl;
			file_manager.rec_ofstream() << "...posterior parameter ensemble saved to " << ss.str() << ".csv" << endl;
		}
		
		pfm.log_event("queueing realizations");
		run_map = pe.add_runs(run_mgr_ptr);
		
		/*run_mgr_ptr->run();
		
		ObservationEnsemble oe(&pest_scenario);
		Covariance obscov = get_obscov();
		oe.draw(num_reals, obscov, &pfm, 1);
		oe.update_from_runs(run_map, run_mgr_ptr);
		ss.str("");
		if (iter != -999)
			ss << file_manager.get_base_filename() << "." << iter << ".post.obsen";
		else
			ss << file_manager.get_base_filename() << ".post.obsen";
		if (binary)
			oe.to_binary(ss.str() + ".jcb");
		else
			oe.to_csv(ss.str() + ".csv");*/
	}
	//fout_rec << "  ---  finished uncertainty analysis calculations  ---  " << endl << endl << endl;
	return pair<ParameterEnsemble,map<int,int>>(pe,run_map);
}


pair<ObservationEnsemble,map<string,double>> LinearAnalysis::process_fosm_reals(RunManagerAbstract* run_mgr_ptr, pair<ParameterEnsemble,map<int, int>>& fosm_real_info, int iter,
														double last_best_phi)
{
	int num_reals = pest_scenario.get_pestpp_options().get_glm_num_reals();
	bool binary = pest_scenario.get_pestpp_options().get_save_binary();
	pfm.log_event("processing FOSM realization runs");
	ObservationEnsemble oe(&pest_scenario, rand_gen_ptr);
	if (num_reals <= 0)
		return pair<ObservationEnsemble, map<string, double>>(oe, map<string, double>());
	//Covariance obscov = get_obscov();
	//oe.draw(num_reals, obscov, &pfm, 1);
	oe.reserve(oe.get_generic_real_names(num_reals), pest_scenario.get_ctl_ordered_obs_names());
	oe.update_from_runs(fosm_real_info.second, run_mgr_ptr);
	if (pest_scenario.get_pestpp_options().get_glm_debug_real_fail())
	{
		vector<string> drop;
		drop.push_back(oe.get_real_names()[0]);
		oe.drop_rows(drop);
	}
	stringstream ss;
	ss.str("");
	if (iter != -999)
		ss << file_manager.get_base_filename() << "." << iter + 1 << ".post.obsen";
	else
		ss << file_manager.get_base_filename() << ".post.obsen";
	if (binary)
	{
		oe.to_binary(ss.str() + ".jcb");
		cout << "...posterior observation ensemble saved to " << ss.str() << ".jcb" << endl;
		file_manager.rec_ofstream() << "...posterior observation ensemble saved to " << ss.str() << ".jcb" << endl;
	}
	else
	{
		oe.to_csv(ss.str() + ".csv");
		cout << "...posterior observation ensemble saved to " << ss.str() << ".csv" << endl;
		file_manager.rec_ofstream() << "...posterior observation ensemble saved to " << ss.str() << ".csv" << endl;
	}
		
	map<string, double> t;
	if (oe.shape().first > 0)
	{
		ofstream& os = file_manager.rec_ofstream();
		ParameterEnsemble pe = fosm_real_info.first;
		double reg_fac = 0.0;
		L2PhiHandler ph(&pest_scenario, &file_manager, &oe, &pe, get_parcov_ptr(),false);
		ph.update(oe, pe);
		L2PhiHandler::phiType pt = L2PhiHandler::phiType::ACTUAL;
		map<string,double>* phi_map = ph.get_phi_map_ptr(pt);
		t = *phi_map;
		os << endl << "  FOSM-based Monte Carlo phi summary:" << endl;
		os << setw(15) << "realization" << setw(20) << "phi" << endl;
		//cout << endl << "  FOSM-based Monte Carlo phi summary:" << endl;
		//cout << setw(15) << "realization" << setw(20) << "phi" << endl;
		for (auto real : oe.get_real_names())
		{
			os << setw(20) << real << setw(20) << phi_map->at(real) << " (" << (phi_map->at(real) / last_best_phi) * 100. << "% of starting phi)" << endl;
			//cout << setw(20) << real << setw(20) << phi_map->at(real) << " (" << (phi_map->at(real) / last_best_phi) * 100. << "% of starting phi)" << endl;
		}
	}
	return pair<ObservationEnsemble,map<string,double>>(oe,t);
}


void LinearAnalysis::throw_error(const string &message)
{
	pfm.log_event("Error in LinearAnalysis:" + message);
	throw runtime_error(message);
}


void LinearAnalysis::load_pst(Pest &pest_scenario, const string &pst_filename)
{
	ifstream ipst(pst_filename);
	if (!ipst.good())
		throw_error("linear_analysis::load_pst() error opening pest control file: " + pst_filename);
	try
	{
		pest_scenario.process_ctl_file(ipst, pst_filename);
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::load_pst() error processing pest control file : " + pst_filename + " : " +string(e.what()));
	}
	ipst.close();
}


void LinearAnalysis::load_jco(Mat& jco, const string& jco_filename)
{
	pfm.log_event("LinearAnalysis::load_jco");
	if (!pest_utils::check_exist_in(jco_filename))
		throw_error("linear_analysis::load_jco() error: jco_filename does not exist");
	string ext = jco_filename.substr(jco_filename.find_last_of(".") + 1);
	pest_utils::upper_ip(ext);
	if ((ext == "JCO") || (ext == "JCB"))
	{
		try
		{
			jco.from_binary(jco_filename);
		}
		catch (exception& e)
		{
			throw_error("linear_analysis::load_jco() error loading jco from binary: " + string(e.what()));
		}
	}
	else
	{
		try
		{
			jco.from_ascii(jco_filename);
		}
		catch (exception& e)
		{
			throw_error("linear_analysis::load_jco() error loading jco from ascii: " + string(e.what()));
		}
	}
}


void LinearAnalysis::load_parcov(const string &parcov_filename)
{
	pfm.log_event("Linear_Analysis::load_parcov from "+parcov_filename);
	if (!pest_utils::check_exist_in(parcov_filename))
		throw_error("linear_analysis::load_parcov() error: parcov_filename does not exist");
	string ext = parcov_filename.substr(parcov_filename.find_last_of(".") + 1);
	vector<string> empty;
	pest_utils::upper_ip(ext);
	if (ext == "UNC")
	{
		try
		{
			parcov.from_uncertainty_file(parcov_filename, empty);
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::load_parcov() error loading parcov from UNC file:" + parcov_filename + " : " + string(e.what()));
		}
	}
	else if (ext == "PST")
	{
		try
		{
			parcov.from_parameter_bounds(parcov_filename, file_manager.rec_ofstream());
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::load_parcov() error loading parcov from PST file:" + parcov_filename + " : " + string(e.what()));
		}
	}
	else
	{
		try
		{
			parcov.from_ascii(parcov_filename);
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::load_parcov() error loading parcov from ASCII file:" + parcov_filename + " : " + string(e.what()));
		}
	}
	
}


void LinearAnalysis::load_obscov(const string &obscov_filename)
{
	pfm.log_event("LinearAnalysis::load_obscov from "+obscov_filename);
	if (!pest_utils::check_exist_in(obscov_filename))
		throw_error("linear_analysis::load_obscov() error: obscov_filename does not exist");
	string ext = obscov_filename.substr(obscov_filename.find_last_of(".") + 1);
	pest_utils::upper_ip(ext);
	vector<string> empty;
	if (ext == "UNC")
	{
		try
		{
			obscov.from_uncertainty_file(obscov_filename,empty);
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::load_obscov() error loading obscov from UNC file:" + obscov_filename + " : " + string(e.what()));
		}
	}
	else if (ext == "PST")
	{
		try
		{
			obscov.from_observation_weights(obscov_filename, file_manager.rec_ofstream());
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::load_obscov() error loading obscov from PST file:" + obscov_filename + " : " + string(e.what()));
		}
	}
	else
	{
		try
		{
			obscov.from_ascii(obscov_filename);
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::load_obscov() error loading obscov from ASCII file:" + obscov_filename + " : " + string(e.what()));
		}
	}
}


LinearAnalysis::LinearAnalysis(Mat &_jacobian, Pest &_pest_scenario, FileManager& _file_manager, PerformanceLog &_pfm, Covariance& _parcov,
	std::mt19937* _rand_gen_ptr): 
	pest_scenario(_pest_scenario),file_manager(_file_manager),
			jacobian(_jacobian), pfm(_pfm),parcov(_parcov), rand_gen_ptr(_rand_gen_ptr)
{
	bool parcov_success = false;
	//if (jacobian.nrow() != pest_scenario.get_nonregul_obs().size())
	//	jacobian = jacobian.get(pest_scenario.get_ctl_ordered_obs_names(), pest_scenario.get_ctl_ordered_adj_par_names());

	/*const string parcov_filename = pest_scenario.get_pestpp_options().get_parcov_filename();
	if (parcov_filename.size() > 0)
	{
		try
		{
			parcov.try_from(pest_scenario, file_manager, true);
			parcov_success = true;
		}
		catch (exception &e)
		{
			pfm.log_event("unable to load parcov from file: " + parcov_filename + ", reverting to parameter bounds "+e.what());
		}
	}
	if (!parcov_success)
	{
		try
		{
			parcov.from_parameter_bounds(pest_scenario);
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::linear_analysis() error setting parcov from parameter bounds:" + string(e.what()));
		}
	}*/
	const string obscov_filename = pest_scenario.get_pestpp_options().get_obscov_filename();
	bool obscov_success = false;
	if (obscov_filename.size() > 0)
	{
		try
		{
			obscov.try_from(pest_scenario, file_manager, false);
			obscov_success = true;
		}
		catch (exception &e)
		{
			pfm.log_event("unable to load obscov from file: " + obscov_filename + ", reverting to observation weights " + e.what());
		}
	}
	if (!obscov_success)
	{
		try
		{
			obscov.from_observation_weights(pest_scenario, file_manager.rec_ofstream());
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::linear_analysis() error setting obscov from observation weights:" + string(e.what()));

		}
	}
	R_sv = -999, G_sv = -999, ImR_sv = -999, V1_sv = -999;
}

void  LinearAnalysis::set_parcov(Mat& _parcov)
{
	parcov = _parcov;
}


void LinearAnalysis::align()
{
	pfm.log_event("LinearAnalysis::align");
	vector<string> common;
	if (jacobian.get_col_names() != parcov.get_col_names())
	{
		try
		{
       		vector<string> tmp_vec = jacobian.get_col_names();
			parcov = parcov.get(tmp_vec);
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::align() error getting aligned parcov: " + string(e.what()));
		}
		catch (...)
		{
			throw_error("linear_analysis::align() error getting aligned parcov");
		}
	}
	if (jacobian.get_row_names() != obscov.get_col_names())
	{
		try
		{
		  vector<string> tmp_vec = jacobian.get_row_names();
		  obscov = obscov.get(tmp_vec);
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::align() error getting aligned obscov: " + string(e.what()));
		}
		catch (...)
		{
			throw_error("linear_analysis::align() error getting aligned parcov");
		}
	}
	for (auto &p : predictions)
	{
		if (*jacobian.cn_ptr() != *p.second.rn_ptr())
		{
			Mat new_pred;
			try
			{
				new_pred = p.second.get(*jacobian.cn_ptr(), *p.second.rn_ptr());
			}
			catch (exception &e)
			{
				throw_error("linear_analysis::align() error getting aligned prediction " + p.first + " : " + string(e.what()));
			}
			predictions[p.first] = new_pred;
		}
	}
}


map<string, double> LinearAnalysis::prior_parameter_variance()
{
	map<string, double> results;
	for (auto &pname : parcov.get_col_names())
		results[pname] = prior_parameter_variance(pname);
	return results;
}


double LinearAnalysis::prior_parameter_variance(string &par_name)
{
	//pfm.log_event("prior_parameter_variance");
	int ipar = find(parcov.rn_ptr()->begin(), parcov.rn_ptr()->end(), par_name) - parcov.rn_ptr()->begin();
	if (ipar == parcov.nrow())
		throw_error("linear_analysis::prior_parameter_variance() error: parameter: " + par_name + " not found");
	const Eigen::SparseMatrix<double>* ptr = parcov.e_ptr();
	double val;
	try
	{
		val = ptr->diagonal()[ipar];
	}
	catch (exception &e)
	{
		stringstream ss;
		ss << ipar;
		throw_error("linear_analysis::prior_parameter_variance() error accessing parameter variance at index " + ss.str() + " : " + string(e.what()));
	}
	return val;
}


map<string, double> LinearAnalysis::posterior_parameter_variance()
{
	map<string, double> results;
	for (auto &pname : parcov.get_col_names())
		results[pname] = posterior_parameter_variance(pname);
	return results;
}


double LinearAnalysis::posterior_parameter_variance(string &par_name)
{
	//pfm.log_event("posterior_parameter_variance");
	if (posterior.nrow() == 0) calc_posterior();
	int ipar = find(posterior.rn_ptr()->begin(), posterior.rn_ptr()->end(), par_name) - posterior.rn_ptr()->begin();
	if (ipar == posterior.nrow())
		throw_error("linear_analysis::posterior_parameter_variance() error: parameter: " + par_name + " not found");
	const Eigen::SparseMatrix<double>* ptr = posterior.e_ptr();
	double val;
	try
	{
		val = ptr->diagonal()[ipar];
	}
	catch (exception &e)
	{
		stringstream ss;
		ss << ipar;
		throw_error("linear_analysis::posterior_parameter_variance() error accessing parameter variance at index " + ss.str() + " : " + string(e.what()));
	}
	return val;
}


Mat LinearAnalysis::posterior_parameter_matrix()
{
	if (posterior.nrow() == 0) calc_posterior();
	return posterior;
}

Mat* LinearAnalysis::posterior_parameter_ptr()
{
	if (posterior.nrow() == 0) calc_posterior();
	Mat* ptr = &posterior;
	return ptr;
}

Covariance LinearAnalysis::posterior_parameter_covariance_matrix()
{
	if (posterior.nrow() == 0) calc_posterior();
	return posterior;
}


double LinearAnalysis::prior_prediction_variance(string &pred_name)
{
	pest_utils::upper_ip(pred_name);
	map<string, Mat>::iterator p_iter = predictions.find(pred_name);
	if (p_iter == predictions.end())
		throw_error("linear_analysis::prior_pred_variance() error: pred:" + pred_name + " not found in predictions");

	if (p_iter->second.e_ptr()->nonZeros() == 0)
		return 0.0;
	double val;
	try
	{
		Eigen::SparseMatrix<double> result = (*p_iter->second.transpose().e_ptr() * *parcov.e_ptr() * *p_iter->second.e_ptr());
		val = result.valuePtr()[0];
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::prior_prediction_variance() error calculating variance :" + string(e.what()));
	}
	return val;
}

map<string, double> LinearAnalysis::prior_prediction_variance()
{
	pfm.log_event("LinearAnalysis::prior_prediction_variance");
	map<string, double> result;
	for (auto &pred : predictions)
	{
		string pname(pred.first);
		result[pname] = prior_prediction_variance(pname);
	}
	return result;
}

double LinearAnalysis::posterior_prediction_variance(string &pred_name)
{
	pest_utils::upper_ip(pred_name);
	map<string, Mat>::iterator p_iter = predictions.find(pred_name);
	if (p_iter == predictions.end())
		throw_error("linear_analysis::prior_pred_variance() error: pred:" + pred_name + " not found in predictions");
	if (p_iter->second.e_ptr()->nonZeros() == 0)
		return 0.0;
	if (posterior.nrow() == 0) calc_posterior();
	double val;
	try
	{

		Eigen::SparseMatrix<double> result = (*p_iter->second.transpose().e_ptr() * *posterior.e_ptr() * *p_iter->second.e_ptr());
		val = result.valuePtr()[0];
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::posterior_prediction_variance() error calculating variance : " + string(e.what()));
	}
	return val;

}

map<string, double> LinearAnalysis::posterior_prediction_variance()
{
	pfm.log_event("LinearAnalysis::prior_prediction_variance");
	map<string, double> result;
	for (auto &pred : predictions)
	{
		string pname(pred.first);
		result[pname] = posterior_prediction_variance(pname);
	}
	return result;
}

void LinearAnalysis::calc_posterior()
{
	pfm.log_event("LinearAnalysis::calc_posterior");
	try
	{
		align();
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::calc_posterior() error in align() : " + string(e.what()));
	}

	pfm.log_event("LinearAnalysis::prior_prediction_variance()::invert obscov");
	try
	{
		obscov.inv_ip(pfm);
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::calc_posterior() error inverting obscov : " + string(e.what()));
	}
	

	try
	{
		pfm.log_event("LinearAnalysis::calc_posterior() form JtQJ");
		Covariance JtQJ(*parcov.rn_ptr(), (*jacobian.transpose().e_ptr() * *obscov.e_ptr() *
			*jacobian.e_ptr()));
		
		pfm.log_event("LinearAnalysis::calc_posterior() invert prior parcov");
		Covariance parcov_inv = parcov.inv();

		pfm.log_event("LinearAnalysis::calc_posterior() form posterior parcov");
		posterior = Covariance(*parcov.rn_ptr(), (*JtQJ.e_ptr() + *parcov_inv.e_ptr()));
		
		pfm.log_event("LinearAnalysis::calc_posterior() invert posterior parcov");
		posterior.inv_ip(pfm);
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::calc_posterior() error calculating posterior : " + string(e.what()));
	}

}


void LinearAnalysis::set_predictions(vector<string> preds, bool forgive)
{
	pfm.log_event("set_predictions");
	const vector<string>* obs_names = jacobian.rn_ptr();
	set<string> oset(obs_names->begin(), obs_names->end());
	for (auto pred : preds)
	{
		pest_utils::upper_ip(pred);
		/*if (predictions.find(pred) != predictions.end())
		{
			throw_error("linear_analysis::set_predictions() error: pred:" + pred + " already in predictions");
		}*/
			
		//if (find(obs_names->begin(), obs_names->end(), pred) != obs_names->end())
		if (oset.find(pred) != oset.end())
		{
			
			Mat mpred;
			pfm.log_event("extracting prediction " + pred + " from jacobian");
			try
			{
				mpred = jacobian.extract(pred, vector<string>());
			}
			catch (exception &e)
			{
				throw_error("linear_analysis::set_predictions() error extracting prediction " + pred + " : " + string(e.what()));
			}
			if (mpred.e_ptr()->nonZeros() == 0)
			{
				pfm.log_event("Prediction " + pred + " has no non-zero entries in jacobian row/");
				/*cerr << endl << "WARNING: Prediction " + pred + " has no non-zero entries in jacobian. " << endl;
				cerr << "         This mean that the adjustable parameters have no effect on " << endl;
				cerr << "         prediction " + pred + ".  The uncertainty for this prediction " << endl;
				cerr << "         is essentially infinite." << endl << endl;*/
			}
			mpred.transpose_ip();
			predictions[pred] = mpred;
		}
		else
		{
			if (!pest_utils::check_exist_in(pred))
			{
				//if the pred is not in the jco rows and its not a file
				//then...if forgive, just continue, otherwise throw
				//forgive is so the glm iter fosm can be done repeatedly
				if ((forgive) && (predictions.find(pred) != predictions.end()))
					continue;
				else
					throw_error("linear_analysis::set_predictions() error: pred: " + pred + " not found in jco rows and is not an accessible file");

			}
				Mat mpred;
			pfm.log_event("loading prediction " + pred + " from ASCII file");
			try
			{
				mpred.from_ascii(pred);
			}
			catch (exception &e)
			{
				throw_error("linear_analysis::set_predictions() error loading prediction " + pred + " from ASCII file :" + string(e.what()));
			}
			pfm.log_event("loading prediction " + pred + " from ASCII file");
			if (mpred.ncol() != 1)
			{
				if (mpred.nrow() == 1)
				{
					mpred.transpose_ip();
				}
				else
				{
					throw_error("linear_analysis::set_predictions() error: pred: " + pred + "must be shape (1,npar)");
				}
			}
			//if the pred vector is not completely aligned with the JCO
			if (*mpred.rn_ptr() != *jacobian.cn_ptr())
			{
				vector<string> missing;
				const vector<string> *mpred_par_names = mpred.rn_ptr();
				for (auto jco_par : *jacobian.cn_ptr())
				{
					if (find(mpred_par_names->begin(), mpred_par_names->end(), jco_par) == mpred_par_names->end())
					{
						missing.push_back(jco_par);
					}
				}
				if (missing.size() > 0)
				{
					stringstream ss;
					for (auto m : missing) ss << m << ",";
					throw_error("linear_analysis::set_predictions() error: parameters missing from pred: " + pred + " : " + ss.str());
				}
				try
				{
					mpred = mpred.get(*jacobian.cn_ptr(), *mpred.cn_ptr());
				}
				catch (exception &e)
				{
					throw_error("linear_analysis::set_predictions() error getting/realigning prediction " + pred + " : " + string(e.what()));
				}
			}
			string pname = mpred.get_col_names()[0];
			//if (predictions.find(pname) != predictions.end())
			//	throw_error("linear_analysis::set_predictions() error: pred:" + pred + " already in predictions");
			if (mpred.e_ptr()->nonZeros() == 0)
			{
				pfm.log_event("Prediction " + pred + " has no non-zero entries in jacobian row/");
				/*cerr << endl << "WARNING: Prediction " + pred + " has no non-zero entries in jacobian. " << endl;
				cerr << "         This mean that the adjustable parameters have no effect on " << endl;
				cerr << "         prediction " + pred + ".  The uncertainty for this prediction " << endl;
				cerr << "         is essential infinite." << endl << endl;*/
			}
			predictions[pname] = mpred;
		}
	}
	
}




map<string, double> LinearAnalysis::like_preds(double val)
{
	map<string, double> result;
	for (auto &pred : predictions)
		result[pred.first] = val;
	return result;
}


void LinearAnalysis::write_par_credible_range(ofstream &fout, string sum_filename, ParameterInfo parinfo,
	Parameters init_pars, Parameters opt_pars, vector<string> ordered_names)
{
    int name_len = 20;
    for (auto &pname : ordered_names)
        name_len = max((int)pname.size(),name_len);
	
	fout << "current parameter uncertainty summary: " << endl << endl;
	fout << setw(name_len) << left << " name" << setw(20) << right << "prior_mean" << setw(20) << "prior_stdev" ;
	fout << setw(20) << "prior_lower_bound" << setw(20) << "prior_upper_bound";
	fout << setw(20) << "post_mean" << setw(20) << "post_stdev";
	fout << setw(20) << "post_lower_bound" << setw(20) << "post_upper_bound" << endl;

	ofstream sout(sum_filename);
	sout << "name,prior_mean,prior_stdev,prior_lower_bound,prior_upper_bound,";
	sout << "post_mean,post_stdev,post_lower_bound,post_upper_bound" << endl;

	map<string, double> prior_vars = prior_parameter_variance();
	map<string, double> post_vars = posterior_parameter_variance();
	vector<string> missing;
	double value,stdev;
	pair<double, double> range;
	for (auto &pname : ordered_names)
	{
		//if (find(jacobian.cn_ptr()->begin(), jacobian.cn_ptr()->end(), pname) == jacobian.cn_ptr()->end())
		if (prior_vars.find(pname) == prior_vars.end())
			missing.push_back(pname);
		else
		{
			//prior
			value = init_pars.get_rec(pname);
			//if (parinfo.get_parameter_rec_ptr(pname)->tranform_type == ParameterRec::TRAN_TYPE::LOG)
			//	value = log10(value);
			//range = get_range(value, prior_vars[pname], parinfo.get_parameter_rec_ptr(pname)->tranform_type);
			stdev = sqrt(prior_vars.at(pname));

			fout << setw(name_len+1) << left << pest_utils::lower_cp(pname) << setw(20) << right << value << setw(20) << stdev << setw(20) <<
				value - (2.0*stdev) << setw(20) << value + (2.0*stdev);
			sout << pest_utils::lower_cp(pname) << "," << value << "," << stdev << "," <<
				value - (2.0*stdev) << "," << value + (2.0*stdev);

			//posterior
			value = opt_pars.get_rec(pname);
			stdev = sqrt(post_vars.at(pname));
			//range = get_range(value, post_vars[pname], parinfo.get_parameter_rec_ptr(pname)->tranform_type);

			fout << setw(20) << value << setw(20) << stdev << setw(20) <<
				value - (2.0*stdev) << setw(20) << value + (2.0*stdev) << endl;
			sout << "," << value << "," << stdev << "," <<
				value - (2.0*stdev) << "," << value + (2.0*stdev) << endl;
		}
	}
	/*if (missing.size() > 0)
	{
		fout << endl;
		fout << "WARNING: the following parameters were not found in the final " << endl;
		fout << "      base parameter jacobian and were subsequently not included " << endl;
		fout << "      in the uncertainty analysis calculations.  This may lead to " << endl;
		fout << "      NON-CONSERVATIVE uncertainty estimates. Please " << endl;
		fout << "      consider including all parameters:" << endl;
		int i = 0;
		for (auto &m : missing)
		{
			fout << setw(20) << m;
			i++;
			if (i == 5)
			{
				fout << endl;
				i = 0;
			}
		}
		fout << endl;*/
	//}
	/*fout << endl;
	fout << "Note: Upper and lower uncertainty bounds reported above are " << endl;
	fout << "      calculated as: <prior,post>_mean +/- (2.0 * <prior,post>_stdev). " << endl;
	fout << "      For log-transformed parameters, the mean, stdev and range are reported " << endl;
	fout << "      with respect to the log of the parameter value. " << endl << endl;*/
}


void LinearAnalysis::write_pred_credible_range(ofstream &fout, string sum_filename,
	map<string,pair<double,double>> init_final_pred_values)
{
    int name_len = 20;
    for (auto &pred : predictions)
        name_len = max((int)pred.first.size(),name_len);
	//fout << endl << "----------------------------------------" << endl;
	fout << "current forecast uncertainty summary: " << endl << endl;
	fout << setw(name_len) << left << " name" << setw(20) << right << "prior_mean";
	fout << setw(20) << "prior_stdev" << setw(20) << "prior_lower_bound";
	fout << setw(20) << "prior_upper_bound" << setw(20) << "post_mean";
	fout << setw(20) << "post_stdev" << setw(20) << "post_lower_bound";
	fout << setw(20) << "post_upper_bound" << endl;

	ofstream sout(sum_filename);
	sout << "name,prior_mean,prior_stdev,prior_lower_bound,prior_upper_bound,";
	sout << "post_mean,post_stdev,post_lower_bound,post_upper_bound" << endl;

	map<string, double> prior_vars = prior_prediction_variance();
	map<string, double> post_vars = posterior_prediction_variance();
	double val, stdev, lower, upper;
	for (auto &pred : predictions)
	{
		val = init_final_pred_values[pred.first].first;
		stdev = sqrt(prior_vars[pred.first]);
		lower = val - (2.0 * stdev);
		upper = val + (2.0 * stdev);
		fout << setw(name_len+1) << left << pest_utils::lower_cp(pred.first) << setw(20) << right << val << setw(20) << stdev;
		fout << setw(20) << lower << setw(20) << upper;
		sout << pest_utils::lower_cp(pred.first) << "," << val << "," << stdev;
		sout << "," << lower << "," << upper;

		val = init_final_pred_values[pred.first].second;
		stdev = sqrt(post_vars[pred.first]);
		lower = val - (2.0 * stdev);
		upper = val + (2.0 * stdev);
		fout << setw(20) << val << setw(20) << stdev;
		fout << setw(20) << lower << setw(20) << upper << endl;
		sout << "," << val << "," << stdev;
		sout << "," << lower << "," << upper << endl;
	}

	/*fout << endl << endl;
	fout << "Note: predictive sensitivity vectors for both prior and " << endl;
	fout << "      posterior uncertainty calculations are the same " << endl;
	fout << "      and were extracted from final base parameter jacobian." << endl;
	fout << "      The upper and lower uncertainty bounds were calculated " << endl;
	fout << "      as: <prior,post>_mean +/- (2.0*<prior,post>_stdev" << endl;*/
	sout.close();
}


void LinearAnalysis::drop_prior_information(const Pest &pest_scenario)
{
	vector<string> pi_names;
	vector<string> obs_names = jacobian.get_row_names();
	const vector<string> nonregul = pest_scenario.get_nonregul_obs();
	if (nonregul.size() < pest_scenario.get_ctl_ordered_obs_names().size())
		for (auto &oname : pest_scenario.get_ctl_ordered_obs_names())
			if ((find(nonregul.begin(), nonregul.end(), oname) == nonregul.end()) &&
				(find(obs_names.begin(),obs_names.end(),oname) != obs_names.end()))
				pi_names.push_back(oname);
	for (auto &oname : pest_scenario.get_ctl_ordered_pi_names())
	{
		if (find(obs_names.begin(), obs_names.end(), oname) != obs_names.end())
			pi_names.push_back(oname);
	}
	//pi_names.insert(pi_names.end(), pest_scenario.get_ctl_ordered_pi_names().begin(),
	//	pest_scenario.get_ctl_ordered_pi_names().end());
	if (pi_names.size() > 0)
	{
		jacobian.drop_rows(pi_names);
	}
	pi_names.clear();
	obs_names.clear();
	obs_names = obscov.get_row_names();
	if (nonregul.size() < pest_scenario.get_ctl_ordered_obs_names().size())
		for (auto &oname : pest_scenario.get_ctl_ordered_obs_names())
			if ((find(nonregul.begin(), nonregul.end(), oname) == nonregul.end()) &&
				(find(obs_names.begin(), obs_names.end(), oname) != obs_names.end()))
				pi_names.push_back(oname);
	for (auto &oname : pest_scenario.get_ctl_ordered_pi_names())
	{
		if (find(obs_names.begin(), obs_names.end(), oname) != obs_names.end())
			pi_names.push_back(oname);
	}
	if (pi_names.size() > 0)
	{
		obscov.drop_rows(pi_names);
	}
}


LinearAnalysis::~LinearAnalysis()
{
}

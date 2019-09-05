#include <vector>
#include <string>
#include <sstream>
#include "Pest.h"
#include "utilities.h"
#include "eigen_tools.h"
#include "covariance.h"
#include "linear_analysis.h"
#include <iomanip>


vector<string> get_common(vector<string> v1, vector<string> v2)
{
	vector<string> common;
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
	ObservationInfo obs_info = pest_scenario.get_ctl_observation_info();
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





void linear_analysis::throw_error(const string &message)
{
	log->error(message);
	throw runtime_error(message);
}


void linear_analysis::load_pst(Pest &pest_scenario, const string &pst_filename)
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


void linear_analysis::load_jco(Mat &jco, const string &jco_filename)
{
	log->log("load_jco");
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
		catch (exception &e)
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
		catch (exception &e)
		{
			throw_error("linear_analysis::load_jco() error loading jco from ascii: " + string(e.what()));
		}
	}
	log->log("load_jco");
}


void linear_analysis::load_parcov(const string &parcov_filename)
{
	log->log("load_parcov from "+parcov_filename);
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
			parcov.from_parameter_bounds(parcov_filename);
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
	log->log("load_parcov from " + parcov_filename);
}


void linear_analysis::load_obscov(const string &obscov_filename)
{
	log->log("load_obscov from "+obscov_filename);
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
			obscov.from_observation_weights(obscov_filename);
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
	log->log("load_obscov from " + obscov_filename);
}


linear_analysis::linear_analysis(Mat _jacobian, Mat _parcov, Mat _obscov, map<string, Mat> _predictions, Logger* _log)
{
	jacobian = _jacobian; parcov = _parcov; obscov = _obscov, predictions = _predictions;
	log = _log;
}


//linear_analysis::linear_analysis(string &jco_filename,Logger* _log)
//{
//	log = _log;
//	string pst_filename = jco_filename;
//	pst_filename.replace(jco_filename.size() - 4, jco_filename.size() - 1, ".pst");
//	if (pest_utils::check_exist_in(pst_filename))
//	{
//		pest_utils::upper_ip(jco_filename);
//		Mat jco;
//		Pest pest_scenario;
//		try
//		{
//			load_jco(jco, jco_filename);
//		}
//		catch (exception &e)
//		{
//			throw_error("linear_analysis::linear_analysis() error loading jco: " + string(e.what()));
//		}
//		try
//		{
//			load_pst(pest_scenario, pst_filename);
//		}
//		catch (exception &e)
//		{
//			throw_error("linear_analysis::linear_analysis() error loading pst: " + string(e.what()));
//		}
//		linear_analysis(jco, pest_scenario);
//	}
//	else
//	{
//		try
//		{
//			load_jco(jacobian, jco_filename);
//		}
//		catch (exception &e)
//		{
//			throw_error("linear_analysis::linear_analysis() error loading jco: " + string(e.what()));
//		}
//	}
//	R_sv = -999, G_sv = -999, ImR_sv = -999, V1_sv = -999;
//}


//linear_analysis::linear_analysis(string &jco_filename, string &parcov_filename, string &obscov_filename,Logger* _log)
//{
//	log = _log;
//	try
//	{
//		load_jco(jacobian, jco_filename);
//	}
//	catch (exception &e)
//	{
//		throw_error("linear_analysis::linear_analysis() error loading jco: " + string(e.what()));
//	}
//	try
//	{
//		load_parcov(parcov_filename);
//	}
//	catch (exception &e)
//	{
//		throw_error("linear_analysis::linear_analysis() error loading parcov: " + string(e.what()));
//	}
//	try
//	{
//		load_obscov(obscov_filename);
//	}
//	catch (exception &e)
//	{
//		throw_error("linear_analysis::linear_analysis() error loading obscov: " + string(e.what()));
//	}
//	R_sv = -999, G_sv = -999, ImR_sv = -999, V1_sv = -999;
//}


//linear_analysis::linear_analysis(Mat _jacobian, Pest pest_scenario,Logger* _log)
//{
//	log = _log;
//	jacobian = _jacobian;
//	try
//	{
//		parcov.from_parameter_bounds(pest_scenario);
//	}
//	catch (exception &e)
//	{
//		throw_error("linear_analysis::linear_analysis() error setting parcov from parameter bounds:" + string(e.what()));
//	}
//	try
//	{
//		obscov.from_observation_weights(pest_scenario);
//	}
//	catch (exception &e)
//	{
//		throw_error("linear_analysis::linear_analysis() error setting obscov from observation weights:" + string(e.what()));
//	}
//	R_sv = -999, G_sv = -999, ImR_sv = -999,V1_sv=-999;
//}

linear_analysis::linear_analysis(Mat &_jacobian, Pest &pest_scenario, FileManager &file_manager, Logger* _log)
{
	bool parcov_success = false;
	log = _log;
	jacobian = _jacobian;
	//if (jacobian.nrow() != pest_scenario.get_nonregul_obs().size())
	//	jacobian = jacobian.get(pest_scenario.get_ctl_ordered_obs_names(), pest_scenario.get_ctl_ordered_adj_par_names());
	const string parcov_filename = pest_scenario.get_pestpp_options().get_parcov_filename();
	if (parcov_filename.size() > 0)
	{
		try
		{
			parcov.try_from(pest_scenario, file_manager, true);
			parcov_success = true;
		}
		catch (exception &e)
		{
			log->warning("unable to load parcov from file: " + parcov_filename + ", reverting to parameter bounds "+e.what());
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
	}
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
			log->warning("unable to load obscov from file: " + obscov_filename + ", reverting to observation weights " + e.what());
		}
	}
	if (!obscov_success)
	{
		try
		{
			obscov.from_observation_weights(pest_scenario);
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::linear_analysis() error setting obscov from observation weights:" + string(e.what()));

		}
	}
	R_sv = -999, G_sv = -999, ImR_sv = -999, V1_sv = -999;
}

void  linear_analysis::set_parcov(Mat &_parcov)
{
	parcov = _parcov;
	//check that everything is kosher
	//vector<string> missing;
	//bool aligned = true;
	//if (jacobian.get_col_names().size() != _parcov.get_row_names.size())
}

//linear_analysis::linear_analysis(Mat* _jacobian, Pest* pest_scenario, Mat* _obscov, Logger* _log)
//{
//	bool parcov_success = false;
//	log = _log;
//	jacobian = *_jacobian;
//	obscov = *_obscov;
//	const string parcov_filename = pest_scenario->get_pestpp_options().get_parcov_filename();
//	if (parcov_filename.size() > 0)
//	{
//		try
//		{
//			load_parcov(parcov_filename);
//			parcov_success = true;
//		}
//		catch (exception &e)
//		{
//			log->warning("unable to load parcov from file: " + parcov_filename + ", reverting to parameter bounds " + e.what());
//			cout << "WARNING: unable to load parcov from file : " << parcov_filename << ", reverting to parameter bounds" << endl;
//		}
//	}
//	if (!parcov_success)
//	{
//		try
//		{
//			parcov.from_parameter_bounds(*pest_scenario);
//		}
//		catch (exception &e)
//		{
//			throw_error("linear_analysis::linear_analysis() error setting parcov from parameter bounds:" + string(e.what()));
//		}
//	}
//	R_sv = -999, G_sv = -999, ImR_sv = -999, V1_sv = -999;
//
//	}

void linear_analysis::align()
{
	log->log("align");
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
	log->log("align");
}

map<string, double> linear_analysis::worth(vector<string> &obs_names)
{
	if (predictions.size() == 0)
		throw_error("linear_analysis::worth() error: no predictions are set");

	log->log("worth");
	try
	{
		align();
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::worth() error in align() : " + string(e.what()));
	}
	for (int i = 0; i != obs_names.size(); i++)
		pest_utils::upper_ip(obs_names[i]);

	//check the inputs
	vector<string> errors;
	const vector<string>* jobs_names = jacobian.rn_ptr();
	for (auto &oname : obs_names)
	{
		if (find(jobs_names->begin(), jobs_names->end(), oname) == jobs_names->end())
			errors.push_back("obs not found in jacobian: "+oname);
	}
	if (errors.size() > 0)
	{
		stringstream ss;
		for (auto &e : errors)
			ss << e << ',';
		throw_error("linear_analysis::worth() errors: " + ss.str());
	}

	//get a list of remaining obs
	vector<string> keep_obs_names;
	for (auto &oname : *jobs_names)
	{
		if (find(obs_names.begin(), obs_names.end(), oname) == obs_names.end())
		{
			keep_obs_names.push_back(oname);
		}
	}
	if (keep_obs_names.size() == 0)
		throw_error("linear_analysis::worth() error at least one observation must not be in the worth group");
	//get a new linear analysis object with only the keep names
	Covariance* parcov_ptr = &parcov;
	map<string,Mat>* pred_ptr = &predictions;
	linear_analysis keep;
	try
	{
		keep = linear_analysis(jacobian.get(keep_obs_names, *jacobian.cn_ptr()), *parcov_ptr, obscov.get(keep_obs_names), *pred_ptr);
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::worth() error getting 'keep' linear analysis object : "+string(e.what()));
	}

	//calculate posterior with all obs and with only keep obs
	map<string, double> org_post, keep_post;
	try
	{
		org_post = posterior_prediction_variance();
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::worth() error calculating posterior variance for original problem : "+string(e.what()));
	}
	try
	{
		keep_post = keep.posterior_prediction_variance();
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::worth() error calculating posterior variance for 'keep' problem : "+string(e.what()));
	}

	map<string, double> results;
	for (auto &pred : predictions)
		results[pred.first] = keep_post[pred.first] - org_post[pred.first];
	log->log("worth");
	return results;
}


map<string, pair<double, double>> linear_analysis::contribution(vector<string> &cond_par_names)
{
	if (predictions.size() == 0)
		throw_error("linear_analysis::contribution() error: no predictions are set");
	log->log("contribution");
	try
	{
		align();
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::contribution() error in align(): " + string(e.what()));
	}

	for (int i = 0; i != cond_par_names.size(); i++)
		pest_utils::upper_ip(cond_par_names[i]);
	//check the inputs
	vector<string> errors;
	const vector<string>* jpar_names = jacobian.cn_ptr();
	for (auto par_name : cond_par_names)
	{
		if (find(jpar_names->begin(), jpar_names->end(), par_name) == jpar_names->end())
			errors.push_back("par not found in jacobian: " + par_name);
	}
	if (errors.size() > 0)
	{
		stringstream ss;
		for (auto e : errors)
			ss << e << ',';
		throw_error("linear_analysis::contribution() errors: " + ss.str());
	}

	//the parameters that will remain uncertain
	vector<string> keep_par_names;
	for (auto pname : *parcov.rn_ptr())
	{
		if (find(cond_par_names.begin(), cond_par_names.end(), pname) == cond_par_names.end())
			keep_par_names.push_back(pname);
	}
	if (keep_par_names.size() == 0)
		throw_error("linear_analysis::contribution() error atleast one parameter must not be in the contribution group");


	map<string,Mat> cond_preds;
	for (auto &pred : predictions)
	{
		try
		{
			cond_preds[pred.first] = pred.second.get(keep_par_names, *pred.second.cn_ptr());
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::contribution() error getting conditional prediction: " + pred.first + " : " + string(e.what()));
		}
	}

	//get a new linear analysis object - we can use the same obscov b/c it doesn't change
	Covariance* obscov_ptr = &obscov;

	linear_analysis cond_la;
	try
	{
		cond_la = linear_analysis(jacobian.get(*jacobian.rn_ptr(), keep_par_names), condition_on(keep_par_names, cond_par_names), *obscov_ptr, cond_preds);
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::contribution() error getting 'conditional' linear analysis object : " + string(e.what()));
	}

	//calculate the prior and posterior with all parameters
	map<string, double> org_prior, org_post;
	try
	{
		org_prior = prior_prediction_variance();
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::contribution() error calculating org prior variance :" + string(e.what()));
	}
	try
	{
		org_post = posterior_prediction_variance();
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::contribution() error calculating org posterior variance :" + string(e.what()));
	}


	//calculate the prior and posterior with some conditioned parameters
	map<string, double> cond_prior, cond_post;
	try
	{
		cond_prior = cond_la.prior_prediction_variance();
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::contribution() error calculating cond prior variance :" + string(e.what()));
	}
	try
	{
		cond_post = cond_la.posterior_prediction_variance();
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::contribution() error calculating cond posterior variance :" + string(e.what()));
	}

	//calculate the reduction in prior and posterior uncertainty for each prediction
	map<string, pair<double, double>> results;
	for (auto pred : predictions)
	{
		pair<double, double> reduction(org_prior[pred.first] - cond_prior[pred.first], org_post[pred.first] - cond_post[pred.first]);
		results[pred.first] = reduction;
	}
	return results;
	log->log("contribution");
}


Covariance linear_analysis::condition_on(vector<string> &keep_par_names, vector<string> &cond_par_names)
{
	log->log("condition_on");
	//C11
	Covariance new_parcov;
	try
	{
		new_parcov = parcov.get(keep_par_names);
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::condition_on() error getting C11 matrix :" + string(e.what()));
	}
	//if parcov is diagonal, then there is no conditioning
	if (parcov.get_mattype() == Mat::MatType::DIAGONAL)
	{
		return new_parcov;
	}
	//the portion of parcov that is becoming "known": C22
	Covariance cond_parcov;
	try
	{
		cond_parcov = parcov.get(cond_par_names);
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::condition_on() error getting C22 matrix :" + string(e.what()));
	}
	//C22^-1
	try
	{
		cond_parcov.inv_ip();
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::condition_on() error inverting C22 matrix :" + string(e.what()));
	}
	//C12
	Mat upper_off_diag;
	try
	{
		upper_off_diag = parcov.get(keep_par_names, cond_par_names);
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::condition_on() error getting C12 matrix :" + string(e.what()));
	}
	//C11 - (C12*C22^-1*C21)
	Covariance new_cond_parcov;
	try
	{
		Eigen::SparseMatrix<double> mcond = *new_parcov.e_ptr() - (*upper_off_diag.e_ptr() * *cond_parcov.e_ptr() *
			*upper_off_diag.transpose().e_ptr());
		new_cond_parcov = Covariance(keep_par_names, mcond);
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::condition_on() error calculating conditional parcov :" + string(e.what()));
	}
	log->log("condition_on");
	return new_cond_parcov;
}


map<string, double> linear_analysis::prior_parameter_variance()
{
	map<string, double> results;
	for (auto &pname : parcov.get_col_names())
		results[pname] = prior_parameter_variance(pname);
	return results;
}


double linear_analysis::prior_parameter_variance(string &par_name)
{
	log->log("prior_parameter_variance");
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
	log->log("prior_parameter_variance");
	return val;
}


map<string, double> linear_analysis::posterior_parameter_variance()
{
	map<string, double> results;
	for (auto &pname : parcov.get_col_names())
		results[pname] = posterior_parameter_variance(pname);
	return results;
}


double linear_analysis::posterior_parameter_variance(string &par_name)
{
	log->log("posterior_parameter_variance");
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
	log->log("posterior_parameter_variance");
	return val;
}


Mat linear_analysis::posterior_parameter_matrix()
{
	if (posterior.nrow() == 0) calc_posterior();
	return posterior;
}

Mat* linear_analysis::posterior_parameter_ptr()
{
	if (posterior.nrow() == 0) calc_posterior();
	Mat* ptr = &posterior;
	return ptr;
}

Covariance linear_analysis::posterior_parameter_covariance_matrix()
{
	if (posterior.nrow() == 0) calc_posterior();
	return posterior;
}


double linear_analysis::prior_prediction_variance(string &pred_name)
{
	pest_utils::upper_ip(pred_name);
	map<string, Mat>::iterator p_iter = predictions.find(pred_name);
	if (p_iter == predictions.end())
		throw_error("linear_analysis::prior_pred_variance() error: pred:" + pred_name + " not found in predicitons");

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

map<string, double> linear_analysis::prior_prediction_variance()
{
	log->log("prior_prediction_variance");
	map<string, double> result;
	for (auto &pred : predictions)
	{
		string pname(pred.first);
		result[pname] = prior_prediction_variance(pname);
	}
	log->log("prior_prediction_variance");
	return result;
}

double linear_analysis::posterior_prediction_variance(string &pred_name)
{
	pest_utils::upper_ip(pred_name);
	map<string, Mat>::iterator p_iter = predictions.find(pred_name);
	if (p_iter == predictions.end())
		throw_error("linear_analysis::prior_pred_variance() error: pred:" + pred_name + " not found in predicitons");
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

map<string, double> linear_analysis::posterior_prediction_variance()
{
	log->log("posterior_prediction_variance");
	map<string, double> result;
	for (auto &pred : predictions)
	{
		string pname(pred.first);
		result[pname] = posterior_prediction_variance(pname);
	}
	log->log("posterior_prediction_variance");
	return result;
}



void linear_analysis::calc_posterior()
{
	log->log("calc_posterior");
	try
	{
		align();
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::calc_posterior() error in align() : " + string(e.what()));
	}

	log->log("invert obscov");
	try
	{
		obscov.inv_ip(log);
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::calc_posterior() error inverting obscov : " + string(e.what()));
	}
	log->log("invert obscov");

	try
	{
		log->log("form JtQJ");
		Covariance JtQJ(*parcov.rn_ptr(), (*jacobian.transpose().e_ptr() * *obscov.e_ptr() *
			*jacobian.e_ptr()));
		log->log("form JtQJ");

		log->log("invert parcov");
		Covariance parcov_inv = parcov.inv();
		log->log("invert parcov");

		log->log("form posterior");
		posterior = Covariance(*parcov.rn_ptr(), (*JtQJ.e_ptr() + *parcov_inv.e_ptr()));
		log->log("form posterior");

		log->log("invert posterior");
		posterior.inv_ip(log);
		log->log("invert posterior");
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::calc_posterior() error calculating posterior : " + string(e.what()));
	}
	log->log("calc_posterior");
}


void linear_analysis::set_predictions(Mat* preds)
{

}

void linear_analysis::set_predictions(vector<string> preds)
{
	log->log("set_predictions");
	const vector<string>* obs_names = jacobian.rn_ptr();
	for (auto pred : preds)
	{
		pest_utils::upper_ip(pred);
		if (predictions.find(pred) != predictions.end())
			throw_error("linear_analysis::set_predictions() error: pred:" + pred + " already in predictions");
		if (find(obs_names->begin(), obs_names->end(), pred) != obs_names->end())
		{
			
			Mat mpred;
			log->log("extracting prediction " + pred + " from jacobian");
			try
			{
				mpred = jacobian.extract(pred, vector<string>());
			}
			catch (exception &e)
			{
				throw_error("linear_analysis::set_predictions() error extracting prediction " + pred + " : " + string(e.what()));
			}
			log->log("extracting prediction " + pred + " from jacobian");
			if (mpred.e_ptr()->nonZeros() == 0)
			{
				log->warning("Prediction " + pred + " has no non-zero entries in jacobian row/");
				cerr << endl << "WARNING: Prediction " + pred + " has no non-zero entries in jacobian. " << endl;
				cerr << "         This mean that the adjustable parameters have no effect on " << endl;
				cerr << "         prediction " + pred + ".  The uncertainty for this prediction " << endl;
				cerr << "         is essentially infinite." << endl << endl;
			}
			mpred.transpose_ip();
			predictions[pred] = mpred;
		}
		else
		{
			if (!pest_utils::check_exist_in(pred))
				throw_error("linear_analysis::set_predictions() error: pred: " + pred + " not found in jco rows and is not an accessible file");
			Mat mpred;
			log->log("loading prediction " + pred + " from ASCII file");
			try
			{
				mpred.from_ascii(pred);
			}
			catch (exception &e)
			{
				throw_error("linear_analysis::set_predictions() error loading prediction " + pred + " from ASCII file :" + string(e.what()));
			}
			log->log("loading prediction " + pred + " from ASCII file");
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
			if (predictions.find(pname) != predictions.end())
				throw_error("linear_analysis::set_predictions() error: pred:" + pred + " already in predictions");
			if (mpred.e_ptr()->nonZeros() == 0)
			{
				log->warning("Prediction " + pred + " has no non-zero entries in jacobian row/");
				cerr << endl << "WARNING: Prediction " + pred + " has no non-zero entries in jacobian. " << endl;
				cerr << "         This mean that the adjustable parameters have no effect on " << endl;
				cerr << "         prediction " + pred + ".  The uncertainty for this prediction " << endl;
				cerr << "         is essential infinite." << endl << endl;
			}
			predictions[pname] = mpred;
		}
	}
	log->log("set_predictions");
}

void linear_analysis::set_predictions(vector<Mat> preds)
{
	throw_error("linear_analysis::set_predictions() not implemented vector<Mat> args");
}

void linear_analysis::svd()
{
	log->log("svd");
	//Eigen::JacobiSVD<Eigen::MatrixXd> svd_fac(*get_normal_ptr()->eptr(), Eigen::DecompositionOptions::ComputeFullU | Eigen::DecompositionOptions::ComputeFullV);
	//Eigen::JacobiSVD<Eigen::MatrixXd> svd_fac(*get_normal_ptr()->eptr(), Eigen::DecompositionOptions::ComputeFullV);
	log->log("svd - factorization");
	Eigen::JacobiSVD<Eigen::MatrixXd> svd_fac;
	try
	{
		svd_fac = Eigen::JacobiSVD<Eigen::MatrixXd>(*get_normal_ptr()->e_ptr(), Eigen::DecompositionOptions::ComputeThinV);
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::svd() error calculating svd: " + string(e.what()));
	}
	log->log("svd - factorization");
	vector<string> sv_names,rs_names;
	stringstream ss;
	for (int i = 0; i != jacobian.cn_ptr()->size(); i++)
	{
		ss.str(string());
		ss.clear();
		ss << "singular_value_" << i;
		sv_names.push_back(ss.str());
		ss.str(string());
		ss.clear();
		ss << "right_singular_vector_" << i;
		rs_names.push_back(ss.str());
	}
	try
	{
		V = Mat(jacobian.get_col_names(), rs_names, svd_fac.matrixV().sparseView());
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::svd() error setting V :" + string(e.what()));
	}
	try
	{
		S = Covariance(sv_names, eigenvec_2_diagsparse(svd_fac.singularValues()));
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::svd() error setting S :" + string(e.what()));
	}

	log->log("svd");
}


Mat* linear_analysis::get_R_ptr(int sv)
{
	if (R_sv != sv)
		build_R(sv);
	Mat* ptr = &R;
	return ptr;

}

void linear_analysis::build_R(int sv)
{
	log->log("build_R");
	if (V.nrow() == 0)
	{
		try
		{
			svd();
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::build_R() error in svd() :" + string(e.what()));
		}
	}
	if (sv == 0)
	{
		try
		{

			Eigen::SparseMatrix<double> Z(V.nrow(), V.ncol());
			Z.setZero();
			R = Mat(V.get_row_names(), V.get_col_names(), Z);
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::build_R() error setting R to zero :" + string(e.what()));
		}
	}
	else if (sv > jacobian.ncol())
	{
		R = V.identity();
	}
	else
	{
		vector<string> cnames;
		vector<string> base_cnames = *V.cn_ptr();
		for (int i = 0; i < sv; i++)
			cnames.push_back(base_cnames[i]);
		log->log("build_R - MM");
		try
		{
			R = Mat(V.get_row_names(), V.get_row_names(), *get_V1_ptr(sv)->e_ptr() * get_V1_ptr(sv)->e_ptr()->transpose());
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::build_R() error calculating R : " + string(e.what()));
		}
		log->log("build_R - MM");
	}
	R_sv = sv;
	log->log("build_R");
}

Mat* linear_analysis::get_V1_ptr(int sv)
{
	if (sv != V1_sv)
		build_V1(sv);
	Mat* ptr = &V1;
	return ptr;
}

void linear_analysis::build_V1(int sv)
{
	log->log("build_V1");
	if (V.nrow() == 0)
	{
		try
		{
			svd();
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::build_V1() error in svd() :" + string(e.what()));
		}
	}
	vector<string> cnames;
	vector<string> base_cnames = *V.cn_ptr();
	for (int i = 0; i < sv; i++)
		cnames.push_back(base_cnames[i]);
	try
	{
		//V1 = Mat(V.get_row_names(), cnames, V.eptr()->leftCols(sv));
		V1 = V.leftCols(sv);
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::build_V1() error forming V1 matrix:" + string(e.what()));
	}
	V1_sv = sv;
	log->log("build_V1");
}

Mat* linear_analysis::get_G_ptr(int sv)
{
	if (G_sv != sv)
		build_G(sv);
	Mat* ptr = &G;
	return ptr;
}

void linear_analysis::build_G(int sv)
{
	log->log("build_G");
	if (V.nrow() == 0)
	{
		try
		{
			svd();
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::build_G() error in svd() :" + string(e.what()));
		}
	}
	Eigen::SparseMatrix<double> s_inv;
	try
	{
		s_inv = S.inv().e_ptr()->topLeftCorner(sv, sv);
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::build_G() error building s_inv matrix : " + string(e.what()));
	}
	log->log("build_G - MMM");
	try
	{
		Eigen::SparseMatrix<double> g = *get_V1_ptr(sv)->e_ptr() * s_inv * get_V1_ptr(sv)->e_ptr()->transpose() *
			jacobian.e_ptr()->transpose() * *obscov.e_ptr();
		log->log("build_G - MMM");
		G = Mat(jacobian.get_col_names(), jacobian.get_row_names(), g);
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::build_G() error calculating G : " + string(e.what()));
	}
	G_sv = sv;
	log->log("build_G");
}

Mat* linear_analysis::get_ImR_ptr(int sv)
{
	if (ImR_sv != sv)
		build_ImR(sv);
	Mat* ptr = &ImR;
	return ptr;
}

void linear_analysis::build_ImR(int sv)
{
	log->log("build_ImR");
	if (sv == 0)
		ImR = parcov.identity();
	else
	{
		if (V.nrow() == 0)
		{
			try
			{
				svd();
			}
			catch (exception &e)
			{
				throw_error("linear_analysis::build_G() error in svd() :" + string(e.what()));
			}
		}
		Mat V2;
		try
		{
			V2 = V.rightCols(V.ncol() - sv);
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::build_ImR() error forming V2 matrix " + string(e.what()));
		}
		log->log("build_ImR - MM");
		try
		{
			ImR = Mat(V2.get_row_names(), V2.get_row_names(), *V2.e_ptr() * V2.e_ptr()->transpose());
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::build_ImR() error calculating ImR :" + string(e.what()));
		}
		log->log("build_ImR - MM");
	}
	ImR_sv = sv;
	log->log("build_ImR");
}


Mat * linear_analysis::get_normal_ptr()
{
	if (normal.nrow() == 0)
		build_normal();
	Mat* ptr = &normal;
	return ptr;
}

void linear_analysis::build_normal()
{
	log->log("build_normal");
	try
	{
		align();
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::build_normal() error in align() : " + string(e.what()));
	}
	log->log("build_normal - MMM");
	try
	{
		normal = Mat(jacobian.get_col_names(), jacobian.get_col_names(), jacobian.e_ptr()->transpose() * *obscov.e_ptr() *
			*jacobian.e_ptr());
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::build_normal() error calculating normal matrix :" + string(e.what()));
	}
	log->log("build_normal - MMM");
	log->log("build_normal");
}

Covariance linear_analysis::first_parameter(int sv)
{
	stringstream ss;
	ss << sv;
	string sv_str = ss.str();

	if (sv >= jacobian.ncol())
		return parcov.zero();
	else
	{
		try
		{
			log->log("first_parameter @" + sv_str);
			Eigen::SparseMatrix<double> first = *get_ImR_ptr(sv)->e_ptr() * *parcov.e_ptr() * *get_ImR_ptr(sv)->e_ptr();
			Covariance f(parcov.get_col_names(), first);
			log->log("first_parameter @" + sv_str);
			return f;
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::first_parameter() error @" +sv_str + " : " + string(e.what()));
		}
	}

}

Covariance linear_analysis::second_parameter(int sv)
{
	stringstream ss;
	ss << sv;
	string sv_str = ss.str();

	if (sv == 0)
		return parcov.zero();
	else if (sv >= jacobian.ncol())
		return parcov.diagonal(1.0E+20);
	else
	{
		try
		{
			log->log("second_parameter @" + sv_str);
			Eigen::SparseMatrix<double> second = *get_G_ptr(sv)->e_ptr() * *obscov.e_ptr() * get_G_ptr(sv)->e_ptr()->transpose();
			Covariance s(parcov.get_col_names(), second);
			log->log("second_parameter @" + sv_str);
			return s;
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::second_parameter() error @"+sv_str+ " : " + string(e.what()));
		}

	}

}

Covariance linear_analysis::third_parameter(int sv)
{
	stringstream ss;
	ss << sv;
	string sv_str = ss.str();
	log->log("third_parameter @"+sv_str);
	if (sv == 0)
		return parcov.zero();
	else if (sv >= jacobian.ncol())
		return parcov.diagonal(1.0E+20);
	else
	{
		Eigen::SparseMatrix<double> GZo, third;
		try
		{
			GZo = *get_G_ptr(sv)->e_ptr() * *omitted_jacobian.e_ptr();
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::third_parameter() error calculating GZo matrix @"+sv_str+ " : " + string(e.what()));
		}
		try
		{
			third = GZo * *omitted_parcov.e_ptr() * GZo.transpose();
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::third_parameter() error calculating thrid matrix@"+sv_str+" : " + string(e.what()));
		}
		return Covariance(parcov.get_col_names(), third);
	}
	log->log("third_parameter @" + sv_str);
}

map<string, double> linear_analysis::first_prediction(int sv)
{
	stringstream ss;
	ss << sv;
	string sv_str = ss.str();
	map<string, double> result;
	if (sv >= jacobian.ncol())
		return like_preds(0.0);
	else
	{
		log->log("first_prediction @" + sv_str);
		Covariance first;
		try
		{
			first = first_parameter(sv);
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::first_prediction() error calculation first matrix : " + string(e.what()));
		}

		for (auto &p : predictions)
		{
			try
			{
				result[p.first] = (p.second.e_ptr()->transpose() * *first.e_ptr() * *p.second.e_ptr()).eval().valuePtr()[0];
			}
			catch (exception &e)
			{
				stringstream ss;
				ss << sv;
				throw_error("linear_analysis::first_prediction() error calculating result for prediction " + p.first + " for sv " + ss.str()+ " : " + string(e.what()));
			}
		}
		log->log("first_prediction @" + sv_str);
		return result;
	}
}

map<string, double> linear_analysis::second_prediction(int sv)
{
	stringstream ss;
	ss << sv;
	string sv_str = ss.str();
	map<string, double> result;
	if (sv == 0)
		return like_preds(0.0);
	else if (sv >= jacobian.ncol())
		return like_preds(1.0E+20);
	else
	{
		log->log("second_prediction @" + sv_str);
		Covariance second;
		try
		{
			second = second_parameter(sv);
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::second_prediction() error calculation second matrix : " + string(e.what()));
		}
		for (auto &p : predictions)
		{
			try
			{
				result[p.first] = (p.second.e_ptr()->transpose() * *second.e_ptr() * *p.second.e_ptr()).eval().valuePtr()[0];
			}
			catch (exception &e)
			{
				throw_error("linear_analysis::second_prediction() error calculating result for prediction " + p.first + " : " + string(e.what()));
			}
		}
		log->log("second_prediction @" + sv_str);
		return result;
	}
}

map<string, double> linear_analysis::third_prediction(int sv)
{
	stringstream ss;
	ss << sv;
	string sv_str = ss.str();
	map<string, double> result;
	if (sv >= jacobian.ncol())
		return like_preds(1.0E+20);
	else if (sv == 0)
	{
		log->log("third_prediction @" + sv_str);
		for (auto &opred : omitted_predictions)
		{
			try
			{
				result[opred.first] = (opred.second.e_ptr()->transpose() * *omitted_parcov.e_ptr() * *opred.second.e_ptr()).eval().valuePtr()[0];
			}
			catch (exception &e)
			{
				throw_error("linear_analysis::third_prediction() error calculating result for prediction " + opred.first + " : " + string(e.what()));
			}
		}
		log->log("third_prediction @" + sv_str);
		return result;
	}
	else
	{
		log->log("third_prediction with 'p' @" + sv_str);
		for (auto &pred : predictions)
		{

			Eigen::SparseMatrix<double, Eigen::RowMajor> p;
			try
			{
				p = (pred.second.e_ptr()->transpose() * *get_G_ptr(sv)->e_ptr() * *omitted_jacobian.e_ptr()).eval();
			}
			catch (exception &e)
			{
				throw_error("linear_analysis::third_prediction() error calculating p vector for prediction " + pred.first + " : " + string(e.what()));
			}
			try
			{
				p = (p - omitted_predictions[pred.first].e_ptr()->transpose()).eval().transpose();
			}
			catch (exception &e)
			{
				throw_error("linear_analysis::third_prediction() error differencing p vector for prediction " + pred.first + " : " + string(e.what()));
			}
			try
			{
				result[pred.first] = (p.transpose() * *omitted_parcov.e_ptr() * p).eval().valuePtr()[0];
			}
			catch (exception &e)
			{
				throw_error("linear_analysis::third_prediction() error calculating result for prediction " + pred.first + " : " + string(e.what()));
			}
		}
		log->log("third_prediction with 'p' @" + sv_str);
		return result;
	}
}

void linear_analysis::extract_omitted(vector<string> &omitted_par_names)
{
	log->log("extract_omitted");
	try
	{
		align();
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::extract_omitted() error in align() : " + string(e.what()));
	}
	if (omitted_par_names.size() == 0)
		throw_error("linear_analysis::extract_omitted() error: omitted par name vector empty");
	const vector<string>* jpar_names = jacobian.cn_ptr();
	vector<string> missing;
	for (auto &o : omitted_par_names)
	{
		pest_utils::upper_ip(o);
		if (find(jpar_names->begin(), jpar_names->end(), o) == jpar_names->end())
			missing.push_back(o);
	}

	if (missing.size() > 0)
	{
		stringstream ss;
		for (auto m : missing)
			ss << m << ",";
		throw_error("linear_analysis::extract_omitted() error some omitted par names were not found: " + ss.str());
	}
	try
	{
		omitted_jacobian = jacobian.extract(vector<string>(), omitted_par_names);
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::extract_omitted() error extracting omitted jco from jco : " + string(e.what()));
	}
	if (omitted_jacobian.e_ptr()->nonZeros() == 0)
	{
		log->warning("omitted Jacobian has no non-zero entries.");
		cerr << endl << "WARNING: omitted parameter Jacobian has no non-zero entries in jacobian. " << endl;
		cerr << "         This mean that the observations are not sensitive to the " << endl;
		cerr << "         omitted parameters." << endl;
	}
	try
	{
		omitted_parcov = parcov.extract(omitted_par_names);
	}
	catch (exception &e)
	{
		throw_error("linear_analysis::extract_omitted() error extracting omitted parcov from parcov : " + string(e.what()));
	}
	for (auto &p : predictions)
	{
		try
		{
			Mat opred = p.second.extract(omitted_par_names, vector<string>());
			if (opred.e_ptr()->nonZeros() == 0)
			{
				log->warning("omitted prediction "+p.first +" has no non-zero entries.");
				cerr << endl << "WARNING: omitted prediction "+p.first+" has no non-zero entries." << endl;
				cerr << "         This mean that this prediction is not sensitive to the " << endl;
				cerr << "         omitted parameters." << endl;
			}
			omitted_predictions[p.first] = opred;
		}
		catch (exception &e)
		{
			throw_error("linear_analysis::extract_omitted() error extracting omitted prediction from prediction "+p.first+ " : " + string(e.what()));
		}

	}
	log->log("extract_omitted");
}

void linear_analysis::extract_omitted(string &omitted_par_name)
{
  vector<string> tmp_vec;
  tmp_vec.push_back(omitted_par_name);
  extract_omitted(tmp_vec);
}


map<string, double> linear_analysis::like_preds(double val)
{
	map<string, double> result;
	for (auto &pred : predictions)
		result[pred.first] = val;
	return result;
}


void linear_analysis::write_par_credible_range(ofstream &fout, string sum_filename, ParameterInfo parinfo,
	Parameters init_pars, Parameters opt_pars, vector<string> ordered_names)
{
	fout << endl << "---------------------------------------" << endl;
	fout << "---- parameter uncertainty summary ----" << endl;
	fout << "---------------------------------------" << endl << endl << endl;
	fout << setw(20) << "name" << setw(20) << "prior_mean" << setw(20) << "prior_stdev" ;
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
		if (find(jacobian.cn_ptr()->begin(), jacobian.cn_ptr()->end(), pname) == jacobian.cn_ptr()->end())
			missing.push_back(pname);
		else
		{
			//prior
			value = init_pars.get_rec(pname);
			//if (parinfo.get_parameter_rec_ptr(pname)->tranform_type == ParameterRec::TRAN_TYPE::LOG)
			//	value = log10(value);
			//range = get_range(value, prior_vars[pname], parinfo.get_parameter_rec_ptr(pname)->tranform_type);
			stdev = sqrt(prior_vars[pname]);

			fout << setw(20) << pname << setw(20) << value << setw(20) << stdev << setw(20) <<
				value - (2.0*stdev) << setw(20) << value + (2.0*stdev);
			sout << pname << "," << value << "," << stdev << "," <<
				value - (2.0*stdev) << "," << value + (2.0*stdev);

			//posterior
			value = opt_pars.get_rec(pname);
			stdev = sqrt(post_vars[pname]);
			//range = get_range(value, post_vars[pname], parinfo.get_parameter_rec_ptr(pname)->tranform_type);

			fout << setw(20) << value << setw(20) << stdev << setw(20) <<
				value - (2.0*stdev) << setw(20) << value + (2.0*stdev) << endl;
			sout << "," << value << "," << stdev << "," <<
				value - (2.0*stdev) << "," << value + (2.0*stdev) << endl;
		}
	}
	if (missing.size() > 0)
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
		fout << endl;
	}
	fout << endl;
	fout << "Note: Upper and lower uncertainty bounds reported above are " << endl;
	fout << "      calculated as: <prior,post>_mean +/- (2.0 * <prior,post>_stdev). " << endl;
	fout << "      For log-transformed parameters, the mean, stdev and range are reported " << endl;
	fout << "      with respect to the log of the parameter value. " << endl << endl;
}

pair<double, double> linear_analysis::get_range(double value, double variance, const ParameterRec::TRAN_TYPE &tt)
{
	double stdev = 0.0, lower = 0.0 , upper = 0.0 ,lvalue = 0.0 ;
	if (variance > numeric_limits<double>::min())
		stdev = sqrt(variance);
	//stdev = sqrt(variance);
	if (tt == ParameterRec::TRAN_TYPE::LOG)
	{
		lvalue = log10(value);
		lower = pow(10, (lvalue - (2.0*stdev)));
		upper = pow(10, (lvalue + (2.0*stdev)));
	}
	else if (tt == ParameterRec::TRAN_TYPE::NONE)
	{
		lower = value - (2.0 * stdev);
		upper = value + (2.0 * stdev);
	}
	else
	{
		stringstream ss;
		ss << "linear_analysis::get_range() unsupported trans type " << static_cast<int>(tt);
		throw PestError(ss.str());
	}
	return pair<double, double>(lower, upper);
}

void linear_analysis::write_pred_credible_range(ofstream &fout, string sum_filename,
	map<string,pair<double,double>> init_final_pred_values)
{
	fout << endl << "----------------------------------------" << endl;
	fout << "---- prediction uncertainty summary ----" << endl;
	fout << "----------------------------------------" << endl << endl << endl;
	fout << setw(20) << "prediction_name" << setw(20) << "prior_mean";
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
		fout << setw(20) << pred.first << setw(20) << val << setw(20) << stdev;
		fout << setw(20) << lower << setw(20) << upper;
		sout << pred.first << "," << val << "," << stdev;
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

	fout << endl << endl;
	fout << "Note: predictive sensitivity vectors for both prior and " << endl;
	fout << "      posterior uncertainty calculations are the same " << endl;
	fout << "      and were extracted from final base parameter jacobian." << endl;
	fout << "      The upper and lower uncertainty bounds were calculated " << endl;
	fout << "      as: <prior,post>_mean +/- (2.0*<prior,post>_stdev" << endl;
	sout.close();
}


void linear_analysis::drop_prior_information(const Pest &pest_scenario)
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


linear_analysis::~linear_analysis()
{
}

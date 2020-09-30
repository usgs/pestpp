#include <random>
#include <iomanip>
#include <unordered_set>
#include <iterator>
#include <limits>
#include <cstddef>
#include "Ensemble.h"
#include "RestartController.h"
#include "utilities.h"
#include "ParamTransformSeq.h"
#include "ObjectiveFunc.h"
#include "RedSVD-h.h"
#include "covariance.h"
#include "PerformanceLog.h"
#include "system_variables.h"
#include "pest_data_structs.h"
#include "eigen_tools.h"

Ensemble::Ensemble(Pest *_pest_scenario_ptr, std::mt19937* _rand_gen_ptr): pest_scenario_ptr(_pest_scenario_ptr),
rand_gen_ptr(_rand_gen_ptr)
{
}

bool Ensemble::try_align_other_rows(PerformanceLog* performance_log, Ensemble& other)
{
	map<string, int> other_real_map = other.get_real_map(), this_real_map = get_real_map();
	bool need_to_align = false;
	for (auto item : this_real_map)
	{
		//if even one realization name is not common, we are done here
		if (other_real_map.find(item.first) == other_real_map.end())
			return false;
		//if this realization name is not at the same location then we do need to reorder
		if (item.second != other_real_map[item.first])
			need_to_align = true;
	}
	//if we made it to here and need_to_align is false, then two ensembles are already row aligned
	if (!need_to_align)
		return false;
	performance_log->log_event("Ensemble::try_align_other_rows(): reorderig other ensemble");
	other.reorder(real_names, vector<string>(),true);
	return true;
}

void Ensemble::check_for_dups()
{
	vector<string> dups;
	set<string> names;
	for (auto &n : var_names)
	{
		if (names.find(n) != names.end())
			dups.push_back(n);
		names.insert(n);
	}

	names.clear();
	for (auto &n : real_names)
	{
		if (names.find(n) != names.end())
			dups.push_back(n);
		names.insert(n);
	}

	if (dups.size() > 0)
	{
		throw_ensemble_error("duplicate var/real names in ensemble: ", dups);
	}

}

void Ensemble::reserve(vector<string> _real_names, vector<string> _var_names)
{
	reals.resize(_real_names.size(), _var_names.size());
	var_names = _var_names;
	real_names = _real_names;
	org_real_names = real_names;
}

Ensemble Ensemble::zero_like()
{
	Eigen::MatrixXd new_reals = Eigen::MatrixXd::Zero(real_names.size(), var_names.size());
	Ensemble new_en(pest_scenario_ptr, rand_gen_ptr);
	new_en.from_eigen_mat(new_reals, real_names, var_names);
	return new_en;
}


void Ensemble::add_2_cols_ip(const vector<string> &other_var_names, const Eigen::MatrixXd &mat)
{
	
	if (shape().first != mat.rows())
		throw_ensemble_error("Ensemble::add_2_cols_ip(): first dimensions don't match");
	
	map<string, int> this_varmap, other_varmap;
	for (int i = 0; i < var_names.size(); i++)
		this_varmap[var_names[i]] = i;
	
	vector<string> missing;
	set<string> svnames(var_names.begin(), var_names.end());
	set<string>::iterator end = svnames.end();
	for (int i = 0; i < other_var_names.size(); i++)
	{
		if (svnames.find(other_var_names[i]) == end)
			missing.push_back(other_var_names[i]);
		other_varmap[other_var_names[i]] = i;
	}
	if (missing.size() > 0)
		throw_ensemble_error("Ensemble::add_2_cols_ip(): the following var names in other were not found", missing);
	for (auto &ovm : other_varmap)
	{
		reals.col(this_varmap[ovm.first]) += mat.col(ovm.second);
	}
}


void Ensemble::add_2_cols_ip(Ensemble &other)
{
	//add values to (a subset of the) columns of reals
	if (shape().first != other.shape().first)
	throw_ensemble_error("Ensemble::add_2_cols_ip(): first dimensions don't match");
	vector<string> other_real_names = other.get_real_names();
	vector<string> mismatch;
	for (int i = 0; i < shape().first; i++)
	{
		if (other_real_names[i] != real_names[i])
			mismatch.push_back(real_names[i]);

	}
	if (mismatch.size() > 0)
	throw_ensemble_error("the following real_names don't match other", mismatch);
	vector<string> other_var_names = other.get_var_names();
	add_2_cols_ip(other_var_names, other.get_eigen());

}


void draw_thread_function(int id, DrawThread &worker, int num_reals, int ies_verbose, map<string, int> idx_map, map<string, double> std_map, exception_ptr &eptr)
{
	try
	{
		worker.work(id, num_reals, ies_verbose, idx_map, std_map);
	}
	catch (...)
	{
		eptr = current_exception();
	}

	return;
}


void Ensemble::draw(int num_reals, Covariance cov, Transformable &tran, const vector<string> &draw_names,
	const map<string, vector<string>> &grouper, PerformanceLog *plog, int level)
{
	//draw names should be "active" var_names (nonzero weight obs and not fixed/tied pars)
	//just a quick sanity check...
	//if ((draw_names.size() > 50000) && (!cov.isdiagonal()))
	//	cout << "  ---  Ensemble::draw() warning: non-diagonal cov used to draw for lots of variables...this might run out of memory..." << endl << endl;

	//matrix to hold the standard normal draws
	Eigen::MatrixXd draws(num_reals, draw_names.size());

	draws.setZero();
	

	//make sure the cov is aligned
	if (cov.get_col_names() != draw_names)
		cov = cov.get(draw_names);

	//make standard normal draws
	plog->log_event("making standard normal draws");
	//RedSVD::sample_gaussian(draws);
	for (int i = 0; i < num_reals; i++)
	{
		for (int j = 0; j < draw_names.size(); j++)
		{
			draws(i, j) = draw_standard_normal(*rand_gen_ptr);
		}
	}

	if (level > 2)
	{
		ofstream f("standard_normal_draws.dat");
		f << draws << endl;
		f.close();
	}
	//Eigen::MatrixXd draws_temp = draws;

	Eigen::VectorXd std = cov.e_ptr()->diagonal().cwiseSqrt();
	map<string, double> std_map;
	for (int i = 0; i < std.size(); i++)
	{
		std_map[var_names[i]] = std(i);
	}
	//if diagonal cov, then scale by std
	if (cov.isdiagonal())
	{
		plog->log_event("scaling by std");

		for (int j = 0; j < draw_names.size(); j++)
		{
			//cout << var_names[j] << " , " << std(j) << endl;
			draws.col(j) *= std(j);
		}
	}
	//if not diagonal, eigen decomp of cov then project the standard normal draws
	else
	{

		if (grouper.size() > 0)
		{
			stringstream ss;
			ss.str("");
			ss << "...drawing by group" << endl;
			plog->log_event(ss.str());
			cout << ss.str();
			map<string, int> idx_map;
			
			for (int i = 0; i < var_names.size(); i++)
				idx_map[var_names[i]] = i;
			vector<string> group_keys;
			for (auto gi : grouper)
				group_keys.push_back(gi.first);
			DrawThread worker(plog, cov, &draws, group_keys, grouper);
			int num_threads = pest_scenario_ptr->get_pestpp_options().get_ies_num_threads();
			if ((num_threads <= 0) || (group_keys.size() == 1))
			{
				worker.work(0, num_reals, level, idx_map, std_map);
			}
			else
			{
				Eigen::setNbThreads(1);
				vector<thread> threads;
				vector<exception_ptr> exception_ptrs;
				plog->log_event("launching draw threads");

				for (int i = 0; i < num_threads; i++)
				{
					exception_ptrs.push_back(exception_ptr());
				}

				for (int i = 0; i < num_threads; i++)
				{
					//threads.push_back(thread(&LocalUpgradeThread::work, &worker, i, iter, cur_lam));

					threads.push_back(thread(draw_thread_function, i, std::ref(worker), num_reals, level, idx_map, std_map,std::ref(exception_ptrs[i])));

				}
				plog->log_event("waiting to join threads");
				//for (auto &t : threads)
				//	t.join();
				for (int i = 0; i < num_threads; ++i)
				{
					if (exception_ptrs[i])
					{
						try
						{
							rethrow_exception(exception_ptrs[i]);
						}
						catch (const std::exception& e)
						{
							ss.str("");
							ss << "thread " << i << "raised an exception: " << e.what();
							throw runtime_error(ss.str());
						}
					}
					threads[i].join();
					if (exception_ptrs[i])
					{
						try
						{
							rethrow_exception(exception_ptrs[i]);
						}
						catch (const std::exception& e)
						{
							ss.str("");
							ss << "thread " << i << "raised an exception: " << e.what();
							throw runtime_error(ss.str());
						}
					}
				}
				plog->log_event("threaded draws done");

			}
			

			//Covariance gcov;
			//
			//RedSVD::RedSymEigen<Eigen::SparseMatrix<double>> eig;

			//Eigen::MatrixXd proj;
			////map<string, int> idx_map;
			//vector<int> idx;
			//for (int i = 0; i < var_names.size(); i++)
			//	idx_map[var_names[i]] = i;
			//for (auto &gi : grouper)
			//{
			//	ss.str("");
			//	ss << "...processing " << gi.first << " with " << gi.second.size() << " elements" << endl;
			//	cout << ss.str();
			//	plog->log_event(ss.str());
			//	if (gi.second.size() == 0)
			//		throw_ensemble_error("Ensemble::draw(): no elements found for group", gi.second);
			//	if (gi.second.size() == 1)
			//	{
			//		ss.str("");
			//		ss << "only one element in group " << gi.first << ", scaling by std";
			//		plog->log_event(ss.str());
			//		int j = idx_map[gi.second[0]];
			//		draws.col(j) *= std(j);
			//		continue;
			//	}

			//	gcov = cov.get(gi.second);
			//	if (level > 2)
			//	{
			//		gcov.to_ascii(gi.first + "_cov.dat");
			//	}

			//	idx.clear();
			//	for (auto n : gi.second)
			//		idx.push_back(idx_map[n]);
			//	if (idx.size() != idx[idx.size() - 1] - idx[0] + 1)
			//		throw_ensemble_error("Ensemble:: draw() error in full cov group draw: idx out of order");

			//	//cout << var_names[idx[0]] << "," << var_names[idx.size() - idx[0]] << endl;
			//	//cout << idx[0] << ',' << idx[idx.size()-1] << ',' <<  idx.size() << " , " << idx[idx.size()-1] - idx[0] <<  endl;

			//	double fac = gcov.e_ptr()->diagonal().minCoeff();
			//	ss.str("");
			//	ss << "min variance for group " << gi.first << ": " << fac;
			//	plog->log_event(ss.str());
			//	Eigen::MatrixXd block = draws.block(0, idx[0], num_reals, idx.size());
			//	ss.str("");
			//	ss << "Randomized Eigen decomposition of full cov for " << gi.second.size() << " element matrix" << endl;
			//	plog->log_event(ss.str());
			//	eig.compute(*gcov.e_ptr() * (1.0/fac), gi.second.size());
			//	//RedSVD::RedSVD<Eigen::SparseMatrix<double>> svd;
			//	//svd.compute(*gcov.e_ptr(),gcov.get_col_names().size());// , gi.second.size());
			//	//cout << svd.singularValues() << endl;
			//	//Eigen::JacobiSVD<Eigen::MatrixXd> svd(gcov.e_ptr()->toDense(), Eigen::ComputeFullU);
			//	//svd.computeU();

			//	proj = (eig.eigenvectors() * (fac *eig.eigenvalues()).cwiseSqrt().asDiagonal());
			//	//proj = (svd.matrixU() * svd.singularValues().asDiagonal());

			//	if (level > 2)
			//	{
			//		ofstream f(gi.first + "_evec.dat");
			//		f << eig.eigenvectors() << endl;
			//		//f << svd.matrixU() << endl;
			//		f.close();
			//		ofstream ff(gi.first + "_sqrt_evals.dat");
			//		ff << (fac * eig.eigenvalues()).cwiseSqrt() << endl;
			//		//ff << svd.singularValues() << endl;
			//		ff.close();
			//		ofstream fff(gi.first+"_proj.dat");
			//		fff << proj << endl;
			//		fff.close();
			//	}
			//	//cout << "block " << block.rows() << " , " << block.cols() << endl;
			//	//cout << " proj " << proj.rows() << " , " << proj.cols() << endl;
			//	plog->log_event("projecting group block");
			//	draws.block(0, idx[0], num_reals, idx.size()) = (proj * block.transpose()).transpose();

			//}
			
		}
		else
		{
			int ncomps = draw_names.size();
			stringstream ss;
			ss << "Randomized Eigen decomposition of full cov using " << ncomps << " components";
			plog->log_event(ss.str());
			double fac = cov.e_ptr()->diagonal().minCoeff();
			ss.str("");
			ss << "min variance: " << fac;
			RedSVD::RedSymEigen<Eigen::SparseMatrix<double>> eig(*cov.e_ptr()* (1.0 / fac));
			//dirty trick alert: apply the abs to make sure all eigen values are positive - nasty!
			Eigen::MatrixXd proj = (eig.eigenvectors() * (fac * eig.eigenvalues()).cwiseAbs().cwiseSqrt().asDiagonal());

			if (level > 2)
			{
				ofstream f("cov_eigenvectors.dat");
				f << eig.eigenvectors() << endl;
				f.close();
				f.open("cov_sqrt_evals.dat");
				f << (fac * eig.eigenvalues()).cwiseAbs().cwiseSqrt() << endl;
				f.close();
				f.open("cov_projection_matrix.dat");
				f << proj << endl;
				f.close();
			}

			//project each standard normal draw in place
			plog->log_event("projecting realizations");
			/*for (int i = 0; i < num_reals; i++)
			{
				draws.row(i) = proj * draws.row(i).transpose();
			}*/
			draws = proj * draws.transpose();
			draws.transposeInPlace();
		}
	}

	//check for invalid values
	plog->log_event("checking realization for invalid values");
	bool found_invalid = false;
	int iv = 0;
	for (int j = 0; j < draw_names.size(); j++)
	{
		iv = 0;
		for (int i = 0; i < num_reals; i++)
		{


			if (OperSys::double_is_invalid(draws(i, j)))
			{
				found_invalid = true;
				iv++;
			}

		}
		if ((level>2) && (iv > 0))
		{
			cout << iv+1 << " invalid values found for " << draw_names[j] << endl;
		}
	}



	real_names = get_generic_real_names(num_reals);

	org_real_names = real_names;
	//add the mean values - using the Transformable instance (initial par value or observed value)
	plog->log_event("resizing reals matrix");
	reals.resize(num_reals, var_names.size());
	reals.setZero(); // zero-weighted obs and fixed/tied pars get zero values here.
	plog->log_event("filling reals matrix and adding mean values");
	vector<string>::const_iterator start = draw_names.begin(), end=draw_names.end(), name;
	set<string> dset(draw_names.begin(), draw_names.end());
	map<string, int> dmap;
	for (int i = 0; i < draw_names.size(); i++)
		dmap[draw_names[i]] = i;
	for (int j = 0; j < var_names.size(); j++)
	{
		//int jj;
		//name = find(start, end, var_names[j]);
		//if (name != end)
		//{
		//	jj = name - start;
		//	reals.col(j) = draws.col(jj).array() + tran.get_rec(var_names[j]);
		//}
		if (dset.find(var_names[j]) != dset.end())
		{
			//Eigen::MatrixXd temp = draws.col(dmap[var_names[j]]);
			//double dtemp = tran.get_rec(var_names[j]);
			reals.col(j) = draws.col(dmap[var_names[j]]).array() + tran.get_rec(var_names[j]);
		}
	}
	if (found_invalid)
	{
		to_csv("trouble.csv");
		throw_ensemble_error("invalid values in realization draws - trouble.csv written");
	}
}

vector<string> Ensemble::get_generic_real_names(int num_reals)
{
	//form some realization names
	vector<string> rnames;
	stringstream ss;
	for (int i = 0; i < num_reals; i++)
	{
		ss.str("");
		ss << i;
		rnames.push_back(ss.str());
	}
	return rnames;
}

pair<Covariance,Covariance> Ensemble::get_empirical_cov_matrices(FileManager* file_manager_ptr)
{
	//Eigen::MatrixXd rmat = get_obs_resid(oe, false); //dont apply ineq constraints
	//ObservationEnsemble(Pest *_pest_scenario_ptr, Eigen::MatrixXd _reals, vector<string> _real_names, vector<string> _var_names);
	Eigen::MatrixXd ercov = reals.transpose() * reals;
	Covariance rcov(var_names, ercov.sparseView());
	Eigen::MatrixXd anom = get_eigen_anomalies();
	Eigen::VectorXd wij;
	double num_reals = static_cast<double>(shape().first);
	double wij_sum = 0;
	double demon = 0;
	for (int i = 0; i < ercov.rows(); i++)
	{
		for (int j = 0; j < ercov.cols(); j++)
		{
			if (i == j)
				continue;
			demon = demon + (ercov(i, j) * ercov(i, j));
		}
	}
	for (int i = 0; i < shape().second; i++)
	{
		for (int j = 0; j < shape().second; j++)
		{
			if (i == j)
				continue;
			wij = anom.col(i).cwiseProduct(anom.col(j));
			wij_sum = wij_sum + (wij.array() - wij.mean()).square().sum();

		}
	}
	double scale = (num_reals / ((num_reals - 1.) * (num_reals - 1.) * (num_reals - 1.))) * wij_sum;
	scale = scale / demon;
	cout << "optimal residual covariance matrix shrinkage factor: " << scale << endl;
	file_manager_ptr->rec_ofstream() << "optimal residual covariance matrix shrinkage factor : " << scale << endl;
		
	Covariance rcov_diag;
	rcov_diag.from_diagonal(rcov);
	Eigen::MatrixXd shrunk = rcov_diag.e_ptr()->toDense();
	Eigen::MatrixXd t = (rcov_diag.e_ptr()->toDense().array() * scale) + (rcov.e_ptr()->toDense().array() * (1. - scale));
	Covariance rcov_shrunk(rcov.get_row_names(), t.sparseView());

	return pair<Covariance,Covariance> (rcov,rcov_shrunk);
}

Covariance Ensemble::get_diagonal_cov_matrix()
{
	//build an empirical diagonal covariance matrix from the realizations


	Eigen::MatrixXd meandiff = get_eigen_anomalies();
	double var;
	vector<Eigen::Triplet<double>> triplets;
	double num_reals = double(reals.rows());
	for (int j = 0; j < var_names.size(); j++)
	{
		//calc variance for this var_name
		var = (reals.col(j).cwiseProduct(reals.col(j))).sum() / num_reals;
		//if (var == 0.0)
		if (var <=1.0e-30)
		{
			stringstream ss;
			ss << "invalid variance for ensemble variable: " << var_names[j] << ": " << var;
			throw_ensemble_error(ss.str());

		}
		triplets.push_back(Eigen::Triplet<double>(j, j, var));
	}

	//allocate and fill a spare Eigen matrix
	Eigen::SparseMatrix<double> mat;
	mat.conservativeResize(triplets.size(), triplets.size());
	mat.setFromTriplets(triplets.begin(), triplets.end());

	return Covariance(var_names, mat, Covariance::MatType::DIAGONAL);
}

Eigen::MatrixXd Ensemble::get_eigen_anomalies(string on_real)
{
	return get_eigen_anomalies(vector<string>(),vector<string>(),on_real);
}

Eigen::MatrixXd Ensemble::get_eigen_anomalies(const vector<string> &_real_names, const vector<string> &_var_names, string on_real)
{
	//get a matrix this is the differences of var_names  realized values from the mean realized value

	//the new Eigen matrix
	Eigen::MatrixXd _reals;
	if ((_real_names.size() == 0) && (_var_names.size() == 0))
		_reals = reals;
	else
		_reals = get_eigen(_real_names, _var_names);


	map<int, double> center_on_map;
	if (on_real.size() > 0)
	{
		vector<string>::iterator it = find(real_names.begin(), real_names.end(), on_real);
		if (it == real_names.end())
			throw runtime_error("Ensemble::get_eigen_mean_diff() error: 'on_real' not found: " + on_real);
		int idx = distance(real_names.begin(), it);
		for (int i = 0; i <_reals.cols(); i++)
			center_on_map[i] = _reals(idx, i);
	}
	else
	{
		for (int j = 0; j < _reals.cols(); j++)
			center_on_map[j] = _reals.col(j).mean();
	}

	//process each var name
	//double mean;
	int s = _reals.rows();
	for (int j = 0; j < _reals.cols(); j++)
	{
		//mean = _reals.col(j).mean();
		//Eigen::MatrixXd temp = _reals.col(j);
		_reals.col(j) = _reals.col(j) - (Eigen::VectorXd::Ones(s) * center_on_map[j]);
		//temp = _reals.col(j);
		//cout << temp << endl;
	}
	return _reals;
}

vector<double> Ensemble::get_mean_stl_var_vector()
{
	vector<double> mean_vec;
	mean_vec.reserve(var_names.size());
	for (int j = 0; j < reals.cols(); j++)
	{
		mean_vec.push_back(reals.col(j).mean());
	}
	return mean_vec;
}

pair<map<string, double>, map<string, double>>  Ensemble::get_moment_maps(const vector<string> &_real_names)
{
	Eigen::VectorXd mean, std;
	if (_real_names.size() == 0)
	{
		mean = reals.colwise().mean();
		Eigen::MatrixXd mean_diff = get_eigen_anomalies();
		std = mean_diff.array().pow(2).colwise().sum().sqrt();
	}
	else
	{
		mean = get_eigen(_real_names,vector<string>()).colwise().mean();
		Eigen::MatrixXd mean_diff = get_eigen_anomalies(_real_names,vector<string>());
		std = mean_diff.array().pow(2).colwise().sum().sqrt();
	}
	map<string, double> mean_map, std_map;
	string name;
	for (int i = 0; i < reals.cols(); i++)
	{
		name = var_names[i];
		mean_map[name] = mean[i];
		std_map[name] = std[i];
	}
	return pair<map<string, double>, map<string, double>>(mean_map,std_map);
}

void Ensemble::replace_col(string var_name, Eigen::VectorXd& vec)
{
	update_var_map();
	if (var_map.find(var_name) == var_map.end())
		throw_ensemble_error("replace_col(): var_name not found: " + var_name);
	if (vec.size() != reals.rows())
	{
		stringstream ss;
		ss << "replace_col(): vec of length " << vec.size() << " not aligned with reals first dimen " << reals.rows();
		throw_ensemble_error(ss.str());
	}
	reals.col(var_map[var_name]) = vec;
}

//Ensemble Ensemble::get_mean()
//{
//	Ensemble new_en(pest_scenario_ptr);
//	new_en.app
//
//}

void Ensemble::from_eigen_mat(Eigen::MatrixXd _reals, const vector<string> &_real_names, const vector<string> &_var_names)
{
	//create a new Ensemble from components
	if (_reals.rows() != _real_names.size())
		throw_ensemble_error("Ensemble.from_eigen_mat() rows != real_names.size");
	if (_reals.cols() != _var_names.size())
		throw_ensemble_error("Ensemble.from_eigen_mat() cols != var_names.size");
	reals = _reals;
	var_names = _var_names;
	real_names = _real_names;
	org_real_names = real_names;
}

void Ensemble::set_eigen(Eigen::MatrixXd _reals)
{
	///reset the reals matrix attribute
	if (_reals.rows() != real_names.size())
		throw_ensemble_error("Ensemble.set_reals() rows != real_names.size");
	if (_reals.cols() != var_names.size())
		throw_ensemble_error("Ensemble.set_reals() cols != var_names.size");
	reals = _reals;
}


void Ensemble::reorder(const vector<string> &_real_names, const vector<string> &_var_names, bool update_org_real_names)
{
	//reorder inplace
	reals = get_eigen(_real_names, _var_names);
	if (_var_names.size() != 0)
		var_names = _var_names;
	if (_real_names.size() != 0)
	{
		real_names = _real_names;
		if (update_org_real_names)
			org_real_names = _real_names;
	}
}

void Ensemble::drop_rows(const vector<int> &row_idxs)
{
	vector<int>::const_iterator start = row_idxs.begin(), end = row_idxs.end();
	vector<string> keep_names;
	for (int ireal = 0; ireal < reals.rows(); ireal++)
		if (find(start, end, ireal) == end)
			keep_names.push_back(real_names[ireal]);
	if (keep_names.size() == 0)
		reals = Eigen::MatrixXd();
	else
		reals = get_eigen(keep_names, vector<string>());
	real_names = keep_names;
}

void Ensemble::drop_rows(const vector<string> &drop_names)
{
	vector<string> keep_names;
	vector<string>::const_iterator start = drop_names.begin(), end = drop_names.end();
	for (auto &n : real_names)
		if (find(start, end, n) == end)
			keep_names.push_back(n);
	if (keep_names.size() == 0)
		reals = Eigen::MatrixXd();
	else
		reals = get_eigen(keep_names, vector<string>());
	real_names = keep_names;

}

void Ensemble::drop_cols(const vector<string>& drop_names)
{
	vector<string> keep_names;
	vector<string>::const_iterator start = drop_names.begin(), end = drop_names.end();
	for (auto& n : var_names)
		if (find(start, end, n) == end)
			keep_names.push_back(n);
	if (keep_names.size() == 0)
		reals = Eigen::MatrixXd();
	else
		reals = get_eigen(vector<string>(), keep_names);
	var_names = keep_names;

}



void Ensemble::keep_rows(const vector<int> &row_idxs)
{
	vector<int>::const_iterator start = row_idxs.begin(), end = row_idxs.end();
	vector<string> keep_names;
	for (int ireal = 0; ireal < reals.rows(); ireal++)
		if (find(start, end, ireal) != end)
			keep_names.push_back(real_names[ireal]);
	reals = get_eigen(keep_names, vector<string>());
	real_names = keep_names;
}

void Ensemble::keep_rows(const vector<string> &keep_names)
{
	//make sure all names are in real_names
	vector<string>::const_iterator start = real_names.begin(), end = real_names.end();
	vector<string> missing;
	for (auto &n : keep_names)
		if (find(start, end, n) == end)
			missing.push_back(n);
	if (missing.size() > 0)
		throw_ensemble_error("Ensemble::keep_rows() error: the following real names not found: ", missing);
	reals = get_eigen(keep_names, vector<string>());
	real_names = keep_names;
}


Eigen::MatrixXd Ensemble::get_eigen(vector<string> row_names, vector<string> col_names, bool update_vmap)
{
	//get a dense eigen matrix from reals by row and col names
	vector<string> missing_rows,missing_cols;
	vector<string>::iterator iter, start = real_names.begin(), end = real_names.end();
	vector<int> row_idxs, col_idxs;

	//check for missing
	vector<string> missing;

	if (row_names.size() > 0)
	{
		map<string, int> real_map;
		for (int i = 0; i < real_names.size(); i++)
			real_map[real_names[i]] = i;
		set<string> real_set(real_names.begin(), end = real_names.end());
		set<string>::iterator end = real_set.end();
		for (auto &name : row_names)
		{
			if (real_set.find(name) == end)
				missing.push_back(name);
			row_idxs.push_back(real_map[name]);
		}
		if (missing.size() > 0)
			throw_ensemble_error("Ensemble.get_eigen() error: the following realization names were not found:", missing);
	}
	if (col_names.size() > 0)
	{
		//map<string, int> var_map;
		//for (int i = 0; i < var_names.size(); i++)
		//	var_map[var_names[i]] = i;
		if (update_vmap)
			update_var_map();

		//set<string> var_set(var_names.begin(), var_names.end());
		map<string,int>::iterator end = var_map.end();
		for (auto &name : col_names)
		{
			if (var_map.find(name) == end)
				missing.push_back(name);
			col_idxs.push_back(var_map[name]);
		}
		if (missing.size() > 0)
			throw_ensemble_error("Ensemble.get_eigen() error: the following variable names were not found:", missing);
	}


	Eigen::MatrixXd mat;

	// only mess with columns, keep rows the same
	if (row_names.size() == 0)
	{
		if (missing_cols.size() > 0)
			throw_ensemble_error("Ensemble.get_eigen() the following col_names not found:", missing_cols);
		mat.resize(real_names.size(), col_names.size());
		int j = 0;
		for (auto &jj : col_idxs)
		{
			mat.col(j) = reals.col(jj);
			j++;
		}
		return mat;
	}

	// only mess with rows, keep cols the same
	if (col_names.size() == 0)
	{
		if (missing_rows.size() > 0)
			throw_ensemble_error("Ensemble.get_eigen() the following row_names not found:", missing_rows);
		mat.resize(row_names.size(), var_names.size());
		int i = 0;
		for (auto &ii : row_idxs)
		{
			mat.row(i) = reals.row(ii);
			i++;
		}
		return mat;
	}

	//the slow one where we are rearranging rows and cols
	if (missing_rows.size() > 0)
		throw_ensemble_error("Ensemble.get_eigen() the following row_names not found:", missing_rows);

	if (missing_cols.size() > 0)
		throw_ensemble_error("Ensemble.get_eigen() the following col_names not found:", missing_cols);
	mat.resize(row_names.size(), col_names.size());
	int i=0, j=0;
	for (auto &ii : row_idxs)
	{
		j = 0;
		for (auto &jj : col_idxs)
		{
			//cout << i << ',' << j << ';' << ii << ',' << jj << endl;
			mat(i, j) = reals(ii, jj);
			j++;
		}
		i++;
	}
	return mat;

}


void Ensemble::to_csv(string file_name)
{
	///write the ensemble to a csv file
	ofstream csv(file_name);
	if (!csv.good())
	{
		throw_ensemble_error("Ensemble.to_csv() error opening csv file " + file_name + " for writing");
	}
	if (pest_scenario_ptr->get_pestpp_options().get_ies_csv_by_reals())
	{
		to_csv_by_reals(csv);
	}
	else
	{
		to_csv_by_vars(csv);
	}
}

void Ensemble::to_csv_by_vars(ofstream &csv)
{
	csv << "var_name";
	
	//for (int ireal = 0; ireal < reals.rows(); ireal++)
	int ireal = 0;
	map<string, int> real_map;
	for (int i = 0; i < real_names.size(); i++)
		real_map[real_names[i]] = i;
	map<string, int>::iterator end = real_map.end();
	for (auto rname : org_real_names)
	{
		if (real_map.find(rname) == end)
			continue;

		ireal = real_map[rname];
		csv << ',' << pest_utils::lower_cp(real_names[ireal]);
	}
	csv << endl;
	for (int ivar = 0; ivar < reals.cols(); ivar++)
	{
		csv << pest_utils::lower_cp(var_names[ivar]);
		for (auto rname : org_real_names)
		{
			if (real_map.find(rname) == end)
				continue;
			ireal = real_map[rname];
			csv << ',' << reals.block(ireal, ivar, 1, 1);
		}
		csv << endl;
	}
}

void Ensemble::to_csv_by_reals(ofstream &csv)
{
	csv << "real_name";
	for (auto &vname : var_names)
		csv << ',' << pest_utils::lower_cp(vname);
	csv << endl;
	//for (int ireal = 0; ireal < reals.rows(); ireal++)
	int ireal = 0;
	map<string, int> real_map;
	for (int i = 0; i < real_names.size(); i++)
		real_map[real_names[i]] = i;
	map<string, int>::iterator end = real_map.end();	
	for (auto rname : org_real_names)
	{
		if (real_map.find(rname) == end)
			continue;
		ireal = real_map[rname];
		csv << pest_utils::lower_cp(real_names[ireal]);
		for (int ivar=0; ivar < reals.cols(); ivar++)
		{
			csv << ',' << reals.block(ireal, ivar,1,1);
		}
		csv << endl;
	}

}

const vector<string> Ensemble::get_real_names(vector<int> &indices)
{
	//returns a subset of real_names by int index
	vector<string> names;
	for (auto &i : indices)
	{
		names.push_back(real_names[i]);
	}
	return names;
}

Eigen::VectorXd Ensemble::get_real_vector(int ireal)
{
	//get a row vector from reals by index
	if (ireal >= shape().first)
	{
		stringstream ss;
		ss << "Ensemble::get_real_vector() : ireal (" << ireal << ") >= reals.shape[0] (" << ireal << ")";
		throw_ensemble_error(ss.str());
	}
	return reals.row(ireal);
}


Eigen::VectorXd Ensemble::get_var_vector(const string& var_name)
{
	update_var_map();
	if (var_map.find(var_name) == var_map.end())
		throw_ensemble_error("Ensemble::get_var_vector(): var_name not found: " + var_name);
	return reals.col(var_map[var_name]);

}

Eigen::VectorXd Ensemble::get_real_vector(const string &real_name)
{
	//get a row vector from reals by realization name
	int idx = find(real_names.begin(), real_names.end(), real_name) - real_names.begin();
	if (idx >= real_names.size())
	{
		stringstream ss;
		ss << "Ensemble::get_real_vector() real_name '" << real_name << "' not found";
		throw_ensemble_error(ss.str());
	}
	return get_real_vector(idx);
}

void Ensemble::update_real_ip(string & rname, Eigen::VectorXd & real)
{
	vector<string>::iterator idx = find(real_names.begin(), real_names.end(), rname);
	if (idx == real_names.end())
	{
		throw_ensemble_error("Ensemble::update_real_ip() real name not found:" + rname);
	}
	//assume real is in order with reals along columns (var names)
	if (real.rows() != reals.cols())
		throw runtime_error("Ensemble::update_real_ip() real has wrong number of entires compared to reals");
	int i = idx - real_names.begin();
	reals.row(i) = real;

}

void Ensemble::update_var_map()
{
	var_map.clear();
	for (int i = 0; i < var_names.size(); i++)
		var_map[var_names[i]] = i;
}

void Ensemble::throw_ensemble_error(string message, vector<string> vec)
{
	stringstream ss;
	ss << ' ';
	for (auto &v : vec)
		ss << v << ',';
	throw_ensemble_error(message + ss.str());
}

void Ensemble::throw_ensemble_error(string message)
{
	string full_message = "Ensemble Error: " + message;
	cout << endl << endl << full_message << endl << endl;
	//cerr << endl << endl << full_message << endl << endl;
	throw runtime_error(full_message);
}

void Ensemble::set_real_names(vector<string>& _real_names)
{
	if (_real_names.size() != real_names.size())
	{
		stringstream ss;
		ss << " set_real_names() _real_names.size(): " << _real_names.size() << " != real_names.size(): " << real_names.size();
		throw_ensemble_error(ss.str());

	}
	real_names = _real_names;
}

Ensemble::~Ensemble()
{
}

map<string, int> Ensemble::get_real_map()
{

	map<string, int> real_map;
	for (int i = 0; i < real_names.size(); i++)
		real_map[real_names[i]] = i;
	return real_map;
}

//Ensemble& Ensemble::operator=(const Ensemble& other)
//{
//	if (this != &other)
//	{
//		pest_scenario_ptr = other.pest_scenario_ptr;
//		rand_gen_ptr = other.rand_gen_ptr;
//	}
//	return *this;
//	
//}

pair<map<string,int>, map<string, int>> Ensemble::prepare_csv(const vector<string> &names, ifstream &csv, bool forgive)
{
	//prepare the input csv for reading checks for compatibility with var_names, forgives extra names in csv
	if (!csv.good())
	{
		throw runtime_error("ifstream not good");
	}

	//process the header
	//any missing header labels will be marked to ignore those columns later
	string line;
	vector<string> header_tokens;
	if (!getline(csv, line))
		throw runtime_error("error reading header (first) line from csv file :");
	pest_utils::strip_ip(line);
	pest_utils::upper_ip(line);
	pest_utils::tokenize(line, header_tokens, ",", false);
	
	//read the index labels
	vector<string> index_tokens,tokens;
	int lcount = 1, nerr = 0;
	stringstream ss;
	
	while (getline(csv, line))
	{
		pest_utils::strip_ip(line);
		pest_utils::upper_ip(line);
		tokens.clear();
		pest_utils::tokenize(line, tokens, ",", false);
		if (header_tokens.size() != tokens.size())
		{
			ss << "wrong number of items on line " << lcount << ", expecting " << header_tokens.size() << " but found " << tokens.size() << endl;
			nerr++;
		}
		index_tokens.push_back(tokens[0]);
		lcount++;
	}

	if (nerr > 0)
	{
		string err = ss.str();
		throw_ensemble_error(err);
	}
	
	unordered_set<string> hset;
	bool csv_by_reals = pest_scenario_ptr->get_pestpp_options().get_ies_csv_by_reals();
	if (csv_by_reals)
		hset = unordered_set<string>(header_tokens.begin(), header_tokens.end());
	else
		hset = unordered_set<string>(index_tokens.begin(), index_tokens.end());

	// check for parameter names that in the pest control file but that are missing from the csv file
	
	vector<string> missing_names;
	unordered_set<string>::iterator end = hset.end();
	string name;
	for (auto &name : names)
		if (hset.find(name) == end)
			missing_names.push_back(name);
	/*if (pest_scenario_ptr->get_pestpp_options().get_ies_csv_by_reals())
	{
		if (missing_names.size() > 0)
		{
			stringstream ss;
			ss << " the following names were not found in the csv file header:" << endl;
			for (auto &n : missing_names) ss << n << endl;
			if (!forgive)
				throw runtime_error(ss.str());
			else
				cout << ss.str() << endl << "continuing anyway..." << endl;
		}
	}*/
	if (missing_names.size() > 0)
	{
		stringstream ss;
		if (csv_by_reals)
			ss << " the following names were not found in the csv file header:" << endl;
		else
			ss << " the following names were not found in the first column of the csv file:" << endl;
		for (auto& n : missing_names) ss << n << endl;
		if (!forgive)
			throw runtime_error(ss.str());
		else
			cout << ss.str() << endl << "continuing anyway..." << endl;
	}

	vector<string> header_names;
	map<string, int> header_info,index_info;
	hset.clear();
	hset = unordered_set<string>(names.begin(), names.end());
	end = hset.end();
	if (csv_by_reals)
	{
		for (int i = 0; i < header_tokens.size(); i++)
		{
			if (hset.find(header_tokens[i]) != end)
			{
				header_info[header_tokens[i]] = i;
			}
		}
		for (int i = 0; i < index_tokens.size(); i++)
			index_info[index_tokens[i]] = i;
		real_names = index_tokens;
	}

	else
	{
		for (int i = 0; i < index_tokens.size(); i++)
		{
			if (hset.find(index_tokens[i]) != end)
			{
				index_info[index_tokens[i]] = i;
			}
		}
		//skip the index label
		for (int i = 1; i < header_tokens.size(); i++)
		{
			header_info[header_tokens[i]] = i;
			real_names.push_back(header_tokens[i]);
		}

			
	}
	
	org_real_names = real_names;
	return pair<map<string,int>,map<string,int>> (header_info,index_info);

}

void Ensemble::extend_cols(Eigen::MatrixXd &_reals, const vector<string> &_var_names)
{
	//add new columns to reals
	vector<string>::iterator start = var_names.begin(), end = var_names.end();
	vector<string> missing;
	for (auto &vname : _var_names)
		if (find(start, end, vname) == end)
			missing.push_back(vname);
	if (missing.size() > 0)
		throw_ensemble_error("add_to_cols(): the following var_names not found: ", missing);
	if (_var_names.size() != _reals.cols())
		throw_ensemble_error("add_to_cols(): _reals.cols() != _var_names.size()");
	if (_reals.rows() != reals.rows())
		throw_ensemble_error("add_to_cols(): _reals.rows() != reals.rows()");
	int i_this;
	string vname;
	for (int i = 0; i < _var_names.size(); i++)
	{
		vname = _var_names[i];
		i_this = find(start, end, vname) - start;
		reals.col(i_this) += _reals.col(i);
	}
}

void Ensemble::append_other_rows(const vector<string>& _real_names, Eigen::MatrixXd& _reals)
{

	//append rows to the end of reals
	if (_reals.cols() != shape().second)
		throw_ensemble_error("append_other_rows(): different number of columns in _reals");
	vector<string> probs;
	set<string> rnames(real_names.begin(), real_names.end());
	set<string>::iterator end = rnames.end();
	for (auto& rname : _real_names)
		//if (find(start, end, rname) != end)
		if (rnames.find(rname) != end)
			probs.push_back(rname);
	if (probs.size() > 0)
		throw_ensemble_error("append_other_rows(): the following _real_names are also in this::real_names: ", probs);
	vector<string> new_real_names = real_names;
	for (auto& rname : _real_names)
	{
		new_real_names.push_back(rname);
		//org_real_names.push_back(rname);
	}
	reals.conservativeResize(new_real_names.size(), var_names.size());

	int iother = 0;
	for (int i = real_names.size(); i < new_real_names.size(); i++)
	{
		reals.row(i) = _reals.row(iother);
		iother++;
	}
	real_names = new_real_names;
}

void Ensemble::append_other_rows(Ensemble &other)
{
	//append rows to the end of reals
	if (other.shape().second != shape().second)
		throw_ensemble_error("append_other_rows(): different number of var_names in other");
	vector<string> probs;
	set<string> vnames(var_names.begin(), var_names.end());
	//vector<string>::iterator start = var_names.begin(), end = var_names.end();
	set<string>::iterator end = vnames.end();
	for (auto &vname : other.get_var_names())
		//if (find(start, end, vname) == end)
		if (vnames.find(vname) == end)
			probs.push_back(vname);
	if (probs.size() > 0)
		throw_ensemble_error("append_other_rows(): the following other::var_names not in this::var_names: ", probs);
	//start = real_names.begin();
	//end = real_names.end();
	set<string> rnames(real_names.begin(), real_names.end());
	end = rnames.end();
	for (auto &rname : other.get_real_names())
		//if (find(start, end, rname) != end)
		if (rnames.find(rname) != end)
			probs.push_back(rname);
	if (probs.size() > 0)
		throw_ensemble_error("append_other_rows(): the following other::real_names are also in this::real_names: ", probs);
	vector<string> new_real_names = real_names;
	for (auto &rname : other.get_real_names())
	{
		new_real_names.push_back(rname);
		//org_real_names.push_back(rname);
	}
	reals.conservativeResize(new_real_names.size(), var_names.size());

	int iother = 0;
	for (int i = real_names.size(); i < new_real_names.size(); i++)
	{
		reals.row(i) = other.get_real_vector(iother);
		iother++;
	}
	real_names = new_real_names;
}

void Ensemble::append(string real_name, const Transformable &trans)
{
	stringstream ss;
	//make sure this real_name isn't ready used
	if (find(real_names.begin(), real_names.end(), real_name) != real_names.end())
	{
		ss << "Ensemble::append() error: real_name '" << real_name << "' already in real_names";
		throw_ensemble_error(ss.str());
	}
	//make sure all var_names are found
	vector<string> keys = trans.get_keys();
	//set<string> tset(keys.begin(), keys.end());
	vector<string> missing;
	Transformable::const_iterator end = trans.end();
	for (auto &n : var_names)
		//if (tset.find(n) == tset.end())
		if (trans.find(n) == end)
			missing.push_back(n);
	if (missing.size() > 0)
	{
		ss.str("");
		ss << "Ensemble::append() error: the following var_names not found: ";
		for (auto &m : missing)
			ss << m << " , ";
		throw_ensemble_error(ss.str());
	}

	reals.conservativeResize(real_names.size() + 1, var_names.size());
	reals.row(real_names.size()) = trans.get_data_eigen_vec(var_names);
	real_names.push_back(real_name);
	if (find(org_real_names.begin(),org_real_names.end(),real_name) == org_real_names.end())
		org_real_names.push_back(real_name);
}

void Ensemble::replace(int idx, const Transformable &trans, string real_name)
{
	stringstream ss;
	
	//make sure all var_names are found
	vector<string> keys = trans.get_keys();
	set<string> tset(keys.begin(), keys.end());
	vector<string> missing;
	for (auto &n : var_names)
		if (tset.find(n) == tset.end())
			missing.push_back(n);
	if (missing.size() > 0)
	{
		ss.str("");
		ss << "Ensemble::append() error: the following var_names not found: ";
		for (auto &m : missing)
			ss << m << " , ";
		throw_ensemble_error(ss.str());
	}

	reals.row(idx) = trans.get_data_eigen_vec(var_names);
	if (real_name.size() > 0)
	{
		//make sure this real_name isn't ready used
		if (find(real_names.begin(), real_names.end(), real_name) != real_names.end())
		{
			ss << "Ensemble::replace() error: real_name '" << real_name << "' already in real_names";
			throw_ensemble_error(ss.str());
		}
		real_names[idx] = real_name;
		org_real_names[idx] = real_name;
	}
}


void Ensemble::to_binary_old(string file_name,bool transposed)
{
	ofstream fout(file_name, ios::binary);
	if (!fout.good())
	{
		throw runtime_error("error opening file for binary ensemble:" + file_name);
	}
	int n_var = var_names.size();
	int n_real =real_names.size();
	int n;
	int tmp;
	double data;
	char par_name[12];
	char obs_name[20];

	// write header
	if (transposed)
	{
		tmp = -n_real;
		fout.write((char*)&tmp, sizeof(tmp));
		tmp = -n_var;
		fout.write((char*)&tmp, sizeof(tmp));

	}
	else
	{
		tmp = -n_var;
		fout.write((char*)&tmp, sizeof(tmp));
		tmp = -n_real;
		fout.write((char*)&tmp, sizeof(tmp));
	}
	

	//write number nonzero elements in jacobian (includes prior information)
	n = reals.size();
	fout.write((char*)&n, sizeof(n));

	//write matrix
	n = 0;
	//map<string, double>::const_iterator found_pi_par;
	//map<string, double>::const_iterator not_found_pi_par;
	//icount = row_idxs + 1 + col_idxs * self.shape[0]
	if (transposed)
	{
		for (int irow = 0; irow < n_var; ++irow)
		{
			for (int jcol = 0; jcol < n_real; ++jcol)
			{
				n = irow + 1 + (jcol * n_var);
				data = reals(jcol, irow);
				fout.write((char*) &(n), sizeof(n));
				fout.write((char*) &(data), sizeof(data));
			}
		}
	}
	else
	{ 
		for (int irow = 0; irow < n_real; ++irow)
		{
			for (int jcol = 0; jcol < n_var; ++jcol)
			{
				n = irow + 1 + (jcol * n_real);
				data = reals(irow, jcol);
				fout.write((char*) &(n), sizeof(n));
				fout.write((char*) &(data), sizeof(data));
			}
		}
	}
	//save parameter names

	if (transposed)
	{
		//save observation and Prior information names
		for (vector<string>::const_iterator b = real_names.begin(), e = real_names.end();
			b != e; ++b) {
			string l = pest_utils::lower_cp(*b);
			pest_utils::string_to_fortran_char(l, par_name, 12);
			fout.write(par_name, 12);
		}
		
		for (vector<string>::const_iterator b = var_names.begin(), e = var_names.end();
			b != e; ++b) {
			string l = pest_utils::lower_cp(*b);
			pest_utils::string_to_fortran_char(l, obs_name, 20);
			fout.write(obs_name, 20);
		}

	}
	else
	{
		for (vector<string>::const_iterator b = var_names.begin(), e = var_names.end();
			b != e; ++b) {
			string l = pest_utils::lower_cp(*b);
			pest_utils::string_to_fortran_char(l, par_name, 12);
			fout.write(par_name, 12);
		}

		//save observation and Prior information names
		for (vector<string>::const_iterator b = real_names.begin(), e = real_names.end();
			b != e; ++b) {
			string l = pest_utils::lower_cp(*b);
			pest_utils::string_to_fortran_char(l, obs_name, 20);
			fout.write(obs_name, 20);
		}
	}
	//save observation names (part 2 prior information)
	fout.close();
}

void Ensemble::check_for_normal(string context)
{
	stringstream ss;
	ss << "realization,variable,value" << endl;
	bool nn_found = false;
	for (int i = 0; i < reals.rows(); i++)
		for (int j = 0; j < reals.cols(); j++)
			if (!isnormal(reals(i, j)) && (reals(i, j) != 0.0))
			{
				ss << real_names[i] << "," << var_names[j] << "," << reals(i,j) << endl;
				nn_found = true;
			}
	if (nn_found)
	{
		ofstream of("ensemble_not_normal.csv");
		of << ss.str();
		ss.str("");
		ss << "Ensemble::check_for_normal() - " << context << " - not normal values found, see file 'ensemble_not_normal.csv'";
		throw_ensemble_error(ss.str());
	}
}

void Ensemble::to_binary(string file_name, bool transposed)
{
	ofstream fout(file_name, ios::binary);
	if (!fout.good())
	{
		throw runtime_error("error opening file for binary ensemble:" + file_name);
	}
	int n_var = var_names.size();
	int n_real = real_names.size();
	int n;
	int tmp;
	double data;
	char par_name[200];
	char obs_name[200];

	// write header
	
	tmp = n_var;
	fout.write((char*)&tmp, sizeof(tmp));
	tmp = n_real;
	fout.write((char*)&tmp, sizeof(tmp));
	
	//write number nonzero elements in jacobian (includes prior information)
	n = reals.size();
	fout.write((char*)&n, sizeof(n));

	//write matrix
	n = 0;
	//map<string, double>::const_iterator found_pi_par;
	//map<string, double>::const_iterator not_found_pi_par;
	//icount = row_idxs + 1 + col_idxs * self.shape[0]
	
	for (int irow = 0; irow < n_real; ++irow)
	{
		for (int jcol = 0; jcol < n_var; ++jcol)
		{
			
			data = reals(irow, jcol);
			fout.write((char*) &(irow), sizeof(irow));
			fout.write((char*) &(jcol), sizeof(jcol));
			fout.write((char*) &(data), sizeof(data));
		}
	}
	
	//save parameter names=
	for (vector<string>::const_iterator b = var_names.begin(), e = var_names.end();
		b != e; ++b) {
		string l = pest_utils::lower_cp(*b);
		pest_utils::string_to_fortran_char(l, par_name, 200);
		fout.write(par_name, 200);
	}

	//save observation and Prior information names
	for (vector<string>::const_iterator b = real_names.begin(), e = real_names.end();
		b != e; ++b) {
		string l = pest_utils::lower_cp(*b);
		pest_utils::string_to_fortran_char(l, obs_name, 200);
		fout.write(obs_name, 200);
	}
	
	//save observation names (part 2 prior information)
	fout.close();
}

map<string, int> Ensemble::from_binary(string file_name, vector<string> &names, bool transposed)
{
	var_names.clear();
	real_names.clear();
	reals.resize(0, 0);
	bool is_new_format = pest_utils::read_binary(file_name, real_names, var_names, reals);
	if ((!is_new_format) && (transposed))
	{
		vector<string> temp = real_names;
		real_names = var_names;
		var_names = temp;
		temp.clear();
		reals.transposeInPlace();
	}

	map<string, int> header_info;
	for (int i = 0; i < var_names.size(); i++)
		header_info[var_names.at(i)] = i;
	org_real_names = real_names;
	return header_info;
}


map<string,int> Ensemble::from_binary_old(string file_name, vector<string> &names, bool transposed)
{
	//load an ensemble from a binary jco-type file.  if transposed=true, reals is transposed and row/col names are swapped for var/real names.
	//needed to store observation ensembles in binary since obs names are 20 chars and par names are 12 chars
	var_names.clear();
	real_names.clear();
	ifstream in;
	in.open(file_name.c_str(), ifstream::binary);
	if (!in.good())
		throw_ensemble_error("Ensemble::from_binary(): error opening binary matrix file: " + file_name);
	int n_col;
	int n_nonzero;
	int n_row;
	int i, j;
	unsigned int n;
	double data;
	
	map<string, int> header_info;

	// read header
	in.read((char*)&n_col, sizeof(n_col));
	in.read((char*)&n_row, sizeof(n_row));
	if (n_col > 0) //throw runtime_error("Ensemble:::from_binary() error: binary matrix file " + file_name + " was produced by deprecated version of PEST");
	{
		char col_name[200];
		char row_name[200];
		n_col = n_col;
		n_row = n_row;
		in.read((char*)&n_nonzero, sizeof(n_nonzero));

		// record current position in file
		streampos begin_sen_pos = in.tellg();

		//advance to col names section
		in.seekg(n_nonzero*(sizeof(double) + sizeof(int) + sizeof(int)), ios_base::cur);

		//read col names
		vector<string>* col_names = &var_names;
		vector<string>* row_names = &real_names;
		

		for (int i_rec = 0; i_rec < n_col; ++i_rec)
		{
			in.read(col_name, 200);
			string temp_col = string(col_name, 200);
			pest_utils::strip_ip(temp_col);
			pest_utils::upper_ip(temp_col);
			col_names->push_back(temp_col);
		}
		//read row names
		for (int i_rec = 0; i_rec < n_row; ++i_rec)
		{
			in.read(row_name, 200);
			string temp_row = pest_utils::strip_cp(string(row_name, 200));
			pest_utils::upper_ip(temp_row);
			row_names->push_back(temp_row);
		}

		//make sure that var_names is compatible with names
		if (var_names.size() != names.size())
		{
			set<string> vset(var_names.begin(), var_names.end());
			set<string> nset(names.begin(), names.end());
			vector<string> diff;
			set_symmetric_difference(vset.begin(), vset.end(), nset.begin(), nset.end(), std::back_inserter(diff));
			throw_ensemble_error("the following names are common between the var names in the binary file and the var names expected", diff);
		}

		//return to data section of file
		in.seekg(begin_sen_pos, ios_base::beg);

		reals.resize(n_row, n_col);
		reals.setZero();
		for (int i_rec = 0; i_rec < n_nonzero; ++i_rec)
		{
			in.read((char*)&(i), sizeof(i));
			in.read((char*)&(j), sizeof(j));
			in.read((char*)&(data), sizeof(data));
			reals(i, j) = data;
		}
		in.close();
		
		set<string> vset(var_names.begin(), var_names.end());
		set<string> nset(names.begin(), names.end());
		vector<string> diff;
		set_symmetric_difference(vset.begin(), vset.end(), nset.begin(), nset.end(), std::back_inserter(diff));
		if (diff.size() > 0)
			throw_ensemble_error("the following names are common between the var names in the binary file and the var names expected", diff);
		
		map<string, int> header_info;
		for (int i = 0; i < col_names->size(); i++)
			header_info[col_names->at(i)] = i;
	}
	else
	{
		char col_name[12];
		char row_name[20];
		n_col = -n_col;
		n_row = -n_row;
		if ((n_col > 1e15) || (n_row > 1e15))
			throw_ensemble_error("Ensemble::from_binary() n_col or n_row > 1e15, something is prob wrong");
		in.read((char*)&n_nonzero, sizeof(n_nonzero));

		// record current position in file
		streampos begin_sen_pos = in.tellg();

		//advance to col names section
		in.seekg(n_nonzero*(sizeof(double) + sizeof(int)), ios_base::cur);

		//read col names
		vector<string>* col_names = &var_names;
		vector<string>* row_names = &real_names;
		if (transposed)
		{
			col_names = &real_names;
			row_names = &var_names;
		}

		for (int i_rec = 0; i_rec < n_col; ++i_rec)
		{
			in.read(col_name, 12);
			string temp_col = string(col_name, 12);
			pest_utils::strip_ip(temp_col);
			pest_utils::upper_ip(temp_col);
			col_names->push_back(temp_col);
		}
		//read row names
		for (int i_rec = 0; i_rec < n_row; ++i_rec)
		{
			in.read(row_name, 20);
			string temp_row = pest_utils::strip_cp(string(row_name, 20));
			pest_utils::upper_ip(temp_row);
			row_names->push_back(temp_row);
		}

		//make sure that var_names is compatible with names
		if (var_names.size() != names.size())
		{
			set<string> vset(var_names.begin(), var_names.end());
			set<string> nset(names.begin(), names.end());
			vector<string> diff;
			set_symmetric_difference(vset.begin(), vset.end(), nset.begin(), nset.end(), std::back_inserter(diff));
			throw_ensemble_error("the following names are common between the var names in the binary file and the var names expected", diff);
		}

		//return to sensitivity section of file
		in.seekg(begin_sen_pos, ios_base::beg);

		reals.resize(n_row, n_col);
		reals.setZero();
		for (int i_rec = 0; i_rec < n_nonzero; ++i_rec)
		{
			in.read((char*)&(n), sizeof(n));
			--n;
			in.read((char*)&(data), sizeof(data));
			j = int(n / (n_row)); // column index
			i = (n - n_row * j) % n_row;  //row index
			reals(i, j) = data;
		}
		if (transposed)
			reals.transposeInPlace();
		in.close();

		map<string, int> header_info;
		for (int i = 0; i < col_names->size(); i++)
			header_info[col_names->at(i)] = i;
	}
	return header_info;

}



void Ensemble::read_csv_by_reals(int num_reals,ifstream &csv, map<string,int> &header_info, map<string,int> &index_info)
{
	//read a csv file to an Ensmeble
	int lcount = 1;
	//vector<vector<double>> vectors;
	double val;
	string line;
	vector<string> tokens;
	string real_id;

	//real_names.clear();
	reals.resize(num_reals, var_names.size());
	reals.setZero();
	int irow = 0;

	map<string, int> var_map;
	for (int i = 0; i < var_names.size(); i++)
		var_map[var_names[i]] = i;

	while (getline(csv, line))
	{
		pest_utils::strip_ip(line);
		tokens.clear();
		pest_utils::tokenize(line, tokens, ",", false);
		if (tokens[tokens.size() - 1].size() == 0)
			tokens.pop_back();

		try
		{
			pest_utils::convert_ip(tokens[0], real_id);
		}
		catch (exception &e)
		{
			stringstream ss;
			ss << "error converting token '" << tokens[0] << "' to <string> real_name on line " << lcount << ": " << line << endl << e.what();
			throw runtime_error(ss.str());
		}
		//real_names.push_back(real_id);
		

		for (auto &hi : header_info)
		{
			try
			{
				//val = pest_utils::convert_cp<double>(tokens[hi.second]);
				val = stod(tokens[hi.second]);
			}
			catch (exception &e)
			{
				stringstream ss;
				ss << "error converting token '" << tokens[hi.second] << "' to double for " << hi.first << " on line " << lcount << " : " << e.what();
				throw runtime_error(ss.str());
			}
			reals(irow, var_map[hi.first]) = val;
		}
		lcount++;
		irow++;

	}
	if (lcount-1 != num_reals)
		throw runtime_error("different number of reals found");
	
}


void Ensemble::read_csv_by_vars(int num_reals, ifstream &csv, map<string, int> &header_info, map<string, int> &index_info)
{
	//read a csv file to an Ensmeble
	int lcount = 1;
	//vector<vector<double>> vectors;
	double val;
	string line;
	vector<string> tokens;
	string var_name;

	//real_names.clear();
	reals.resize(num_reals, var_names.size());
	reals.setZero();
	int var_idx;

	map<string, int> var_map;
	for (int i = 0; i < var_names.size(); i++)
		var_map[var_names[i]] = i;

	while (getline(csv, line))
	{
		
		pest_utils::strip_ip(line);
		pest_utils::upper_ip(line);
		tokens.clear();
		pest_utils::tokenize(line, tokens, ",", false);
		if (tokens[tokens.size() - 1].size() == 0)
			tokens.pop_back();

		try
		{
			pest_utils::convert_ip(tokens[0], var_name);
		}
		catch (exception &e)
		{
			stringstream ss;
			ss << "error converting token '" << tokens[0] << "' to <string> var_name on line " << lcount << ": " << line << endl << e.what();
			throw runtime_error(ss.str());
		}
		//real_names.push_back(real_id);
		var_idx = var_map[var_name];
		for (auto &hi : header_info)
		{
			try
			{
				val = pest_utils::convert_cp<double>(tokens[hi.second]);
			}
			catch (exception &e)
			{
				stringstream ss;
				ss << "error converting token '" << tokens[hi.second] << "' to double for " << hi.first << " on line " << lcount << " : " << e.what();
				throw runtime_error(ss.str());
			}
			reals(hi.second-1, var_idx) = val;
		}
		lcount++;
		

	}
	
}




ParameterEnsemble::ParameterEnsemble(Pest *_pest_scenario_ptr, std::mt19937* rand_gen_ptr):Ensemble(_pest_scenario_ptr, rand_gen_ptr)
{
	par_transform = pest_scenario_ptr->get_base_par_tran_seq();
	tstat = transStatus::CTL;
}

ParameterEnsemble::ParameterEnsemble(Pest *_pest_scenario_ptr, std::mt19937* _rand_gen_ptr, Eigen::MatrixXd _reals,
	vector<string> _real_names, vector<string> _var_names):Ensemble(_pest_scenario_ptr, _rand_gen_ptr)
{
	par_transform = pest_scenario_ptr->get_base_par_tran_seq();
	if (_reals.rows() != _real_names.size())
		throw_ensemble_error("ParameterEnsemble() _reals.rows() != _real_names.size()");
	if (_reals.cols() != _var_names.size())
		throw_ensemble_error("ParameterEnsemble() _reals.cols() != _var_names.size()");
	reals = _reals;
	var_names = _var_names;
	real_names = _real_names;
	tstat = transStatus::CTL;
}

void ParameterEnsemble::draw_uniform(int num_reals, vector<string> par_names, PerformanceLog* plog, int level, ofstream& frec)
{
	var_names = par_names;
	ParamTransformSeq ptrans = pest_scenario_ptr->get_base_par_tran_seq();
	Parameters lb = pest_scenario_ptr->get_ctl_parameter_info().get_low_bnd(par_names);
	Parameters ub = pest_scenario_ptr->get_ctl_parameter_info().get_up_bnd(par_names);
	ptrans.ctl2numeric_ip(lb);
	ptrans.ctl2numeric_ip(ub);
	Parameters dist = ub - lb;
	reals.resize(num_reals, var_names.size());
	reals.setZero();

	vector<double> vals;
	for (int j=0;j<par_names.size();j++)
	{
		vals = uniform_draws(num_reals, lb[par_names[j]], ub[par_names[j]], *rand_gen_ptr);
		reals.col(j) = stlvec_2_eigenvec(vals);
	}

	stringstream ss;
	real_names.clear();
	for (int i = 0; i < num_reals; i++)
	{
		ss.str("");
		ss << "member_" << i;
		real_names.push_back(ss.str());
	}
	org_real_names = real_names;
	tstat = ParameterEnsemble::transStatus::NUM;
}

void ParameterEnsemble::draw(int num_reals, Parameters par, Covariance &cov, PerformanceLog *plog, int level, ofstream& frec)
{
	///draw a parameter ensemble
	var_names = pest_scenario_ptr->get_ctl_ordered_adj_par_names(); //only draw for adjustable pars
	vector<string> cov_names = cov.get_col_names();
	set<string> scov_names(cov_names.begin(), cov_names.end());
	var_names.clear();
	for (auto name : pest_scenario_ptr->get_ctl_ordered_adj_par_names())
	{
		if (scov_names.find(name) != scov_names.end())
			var_names.push_back(name);
	}
	
	//Parameters par = pest_scenario_ptr->get_ctl_parameters();
	par_transform.active_ctl2numeric_ip(par);//removes fixed/tied pars
	tstat = transStatus::NUM;
	ParameterGroupInfo pgi = pest_scenario_ptr->get_base_group_info();
	//vector<string> group_names = pgi.get_group_names();
	vector<string> group_names = pest_scenario_ptr->get_ctl_ordered_par_group_names();
	vector<string> vars_in_group,sorted_var_names;
	map<string, vector<string>> grouper;
	sorted_var_names.reserve(var_names.size());
	bool same = true;
	if (pest_scenario_ptr->get_pestpp_options().get_ies_group_draws())
	{
		for (auto &group : group_names)
		{
			vars_in_group.clear();
			for (auto name : var_names)
			{
				if (pgi.get_group_rec_ptr(name)->name == group)
					vars_in_group.push_back(name);
			}
			if (vars_in_group.size() == 0)
				continue;
			sort(vars_in_group.begin(), vars_in_group.end());
			sorted_var_names.insert(sorted_var_names.end(), vars_in_group.begin(), vars_in_group.end());

			grouper[group] = vars_in_group;
		}

		//check
		if (var_names.size() != sorted_var_names.size())
			throw_ensemble_error("sorted par names not equal to org par names");
		for (int i = 0; i < var_names.size(); i++)
			if (var_names[i] != sorted_var_names[i])
			{
				same = false;
				break;
			}
		if (!same)
		{
			plog->log_event("parameters not grouped by parameter groups, reordering par ensemble");
			cout << "parameters not grouped by parameter groups, reordering par ensemble" << endl;
			var_names = sorted_var_names;
		}
	}
	Ensemble::draw(num_reals, cov, par, var_names, grouper, plog, level);
	/*map<string, int> header_info;
	for (int i = 0; i < var_names.size(); i++)
		header_info[var_names[i]] = i;
	
	ParameterInfo pi = pest_scenario_ptr->get_ctl_parameter_info();
	ParameterRec::TRAN_TYPE ft = ParameterRec::TRAN_TYPE::FIXED;
	for (auto name : pest_scenario_ptr->get_ctl_ordered_par_names() )
	{
		if (pi.get_parameter_rec_ptr(name)->tranform_type == ft)
		{
			fixed_names.push_back(name);
		}
	}*/
	//fill_fixed(header_info);
	//save_fixed();
	if (!same)
	{

		reorder(vector<string>(), pest_scenario_ptr->get_ctl_ordered_adj_par_names());
	}
	if (pest_scenario_ptr->get_pestpp_options().get_ies_enforce_bounds())
		//dont enforce parchglim-style here - if a pars initial value is near its bound, lots of realizations 
		//will be near zero length.  
		enforce_limits(plog, false);


}

void ParameterEnsemble::set_pest_scenario(Pest *_pest_scenario)
{
	pest_scenario_ptr = _pest_scenario;
	par_transform = pest_scenario_ptr->get_base_par_tran_seq();
}

void ParameterEnsemble::set_zeros()
{
	reals.setZero();
}

ParameterEnsemble ParameterEnsemble::zeros_like()
{

	Eigen::MatrixXd new_reals = Eigen::MatrixXd::Zero(real_names.size(), var_names.size());
	
	ParameterEnsemble new_en(pest_scenario_ptr, rand_gen_ptr);
	new_en.from_eigen_mat(new_reals, real_names, var_names);
	return new_en;


}

map<int,int> ParameterEnsemble::add_runs(RunManagerAbstract *run_mgr_ptr,const vector<int> &real_idxs)
{
	//add runs to the run manager using int indices
	map<int,int> real_run_ids;

	//get the pars and transform to be in sync with ensemble trans status
	Parameters pars = pest_scenario_ptr->get_ctl_parameters();
	if (tstat == transStatus::NUM)
	{
		par_transform.active_ctl2numeric_ip(pars);
	}
	else if (tstat == transStatus::MODEL)
	{
		par_transform.active_ctl2model_ip(pars);
	}
	
	Parameters pars_real = pars;
	Eigen::VectorXd evec;
	vector<double> svec;
	int run_id;
	vector<string> run_real_names;
	if (real_idxs.size() > 0)
		for (auto i : real_idxs)
			run_real_names.push_back(real_names[i]);
	else
		run_real_names = real_names;
	int idx;
	map<string, int> rmap;

	for (int i = 0; i < real_names.size(); i++)
		rmap[real_names[i]] = i;
	vector<string> nn;
	for (auto &rname : run_real_names)
	{
		//idx = find(real_names.begin(), real_names.end(), rname) - real_names.begin();
		idx = rmap[rname];
		//Eigen::VectorXd rvector = get_real_vector(idx);
		pars_real = pars;
		pars_real.update_without_clear(var_names, get_real_vector(idx));
		//make sure the pars are in the right trans status
		if (tstat == ParameterEnsemble::transStatus::CTL)
			par_transform.active_ctl2model_ip(pars_real);
		else if (tstat == ParameterEnsemble::transStatus::NUM)
			par_transform.numeric2model_ip(pars_real);
		replace_fixed(rname, pars_real);
		nn = pars_real.get_notnormal_keys();
		if (nn.size() > 0)
		{
			stringstream ss;
			ss << "ParameterEnsemble:: add_runs() error: denormal values for realization " << rname << " : ";
			for (auto n : nn)
				ss << n << ",";
			throw_ensemble_error(ss.str());
		}
		run_id = run_mgr_ptr->add_run(pars_real,"realization: "+rname);
		real_run_ids[idx]  = run_id;
	}
	return real_run_ids;
}

void ParameterEnsemble::from_eigen_mat(Eigen::MatrixXd mat, const vector<string> &_real_names, const vector<string> &_var_names, ParameterEnsemble::transStatus _tstat)
{
	//create a par ensemble from components
	vector<string> missing;
	/*
	vector<string>::const_iterator start = _var_names.begin();
	vector<string>::const_iterator end = _var_names.end();*/

	vector<string> vnames = pest_scenario_ptr->get_ctl_ordered_par_names();
	set<string> vset(vnames.begin(), vnames.end());
	set<string>::iterator end = vset.end();

	/*for (auto &name : pest_scenario_ptr->get_ctl_ordered_par_names())
		if (find(start, end, name) == end)
			missing.push_back(name);*/
	for (auto &name : _var_names)
	{
		if (vset.find(name) == end)
			missing.push_back(name);
	}
	if (missing.size() > 0)
		throw_ensemble_error("ParameterEnsemble.from_eigen_mat() the following par names not found: ", missing);
	Ensemble::from_eigen_mat(mat, _real_names, _var_names);
	tstat = _tstat;
}


void ParameterEnsemble::from_binary(string file_name)
{
	//overload for ensemble::from_binary - just need to set tstat
	vector<string> names = pest_scenario_ptr->get_ctl_ordered_par_names();
	map<string,int> header_info = Ensemble::from_binary(file_name, names, false);
	ParameterInfo pi = pest_scenario_ptr->get_ctl_parameter_info();
	ParameterRec::TRAN_TYPE ft = ParameterRec::TRAN_TYPE::FIXED;
	for (auto &name : var_names)
	{
		if (pi.get_parameter_rec_ptr(name)->tranform_type == ft)
		{
			fixed_names.push_back(name);
		}
	}
	fill_fixed(header_info);
	save_fixed();
	tstat = transStatus::CTL;

}

//ParameterEnsemble ParameterEnsemble::get_new(const vector<string> &_real_names, const vector<string> &_var_names)
//{
//	
//	ParameterEnsemble new_pe(pest_scenario_ptr);
//	new_pe.set_var_names(_var_names);
//	new_pe.set_real_names(_real_names);
//	new_pe.set_eigen(get_eigen(_real_names, _var_names))
//
//	//return ParameterEnsmeble(get_eigen(_real_names,_var_names),)
//}

void ParameterEnsemble::from_csv(string file_name)
{
	ifstream csv(file_name);
	if (!csv.good())
		throw runtime_error("error opening parameter csv " + file_name + " for reading");
	bool csv_by_reals = pest_scenario_ptr->get_pestpp_options().get_ies_csv_by_reals();
	//var_names = pest_scenario_ptr->get_ctl_ordered_adj_par_names();
	var_names = pest_scenario_ptr->get_ctl_ordered_par_names();
	//map<string,int>header_info = prepare_csv(var_names, csv, true);
	pair<map<string, int>, map<string, int>> p = prepare_csv(var_names, csv, true);
	map<string, int> header_info = p.first, index_info = p.second;

	//blast through the file to get number of reals
	/*string line;
	int num_reals = 0;
	while (getline(csv, line))
		num_reals++;
	csv.close();*/
	int num_reals;
	if (csv_by_reals)
		num_reals = index_info.size();
	else
		num_reals = header_info.size();

	//make sure all adjustable parameters are present
	vector<string> missing;
	map<string, int> *map_ptr;
	if (csv_by_reals)
		map_ptr = &header_info;
	else
		map_ptr = &index_info;
	map<string, int>::iterator end = map_ptr->end();

	for (auto p : pest_scenario_ptr->get_ctl_ordered_adj_par_names())
		if (map_ptr->find(p) == end)
			missing.push_back(p);
	if (missing.size() > 0)
		throw_ensemble_error("ParameterEnsemble.from_csv() error: the following adjustable pars not in csv:",missing);

	csv.close();
	csv.open(file_name);
	if (!csv.good())
		throw runtime_error("error re-opening parameter csv " + file_name + " for reading");
	string line;
	getline(csv, line);
	if (csv_by_reals)
		Ensemble::read_csv_by_reals(num_reals, csv, header_info, index_info);
	else
		Ensemble::read_csv_by_vars(num_reals, csv, header_info, index_info);
	ParameterInfo pi = pest_scenario_ptr->get_ctl_parameter_info();
	ParameterRec::TRAN_TYPE ft = ParameterRec::TRAN_TYPE::FIXED;
	for (auto &name : var_names)
	{
		if (pi.get_parameter_rec_ptr(name)->tranform_type == ft)
		{
			fixed_names.push_back(name);
		}
	}
	fill_fixed(header_info);
	save_fixed();

	tstat = transStatus::CTL;
	org_real_names = real_names;

}

void ParameterEnsemble::fill_fixed(const map<string, int> &header_info)
{
	if (fixed_names.size() == 0)
		return;
	map<string, int> var_map;
	for (int i = 0; i < var_names.size(); i++)
		var_map[var_names[i]] = i;

	Parameters pars = pest_scenario_ptr->get_ctl_parameters();
	int c = 0;
	for (auto &fname : fixed_names)
	{
		if (header_info.find(fname) == header_info.end())
		{
			Eigen::VectorXd vec(reals.rows());
			vec.setOnes();
			//reals.col(var_map[fname]) = reals.col(var_map[fname]) + (pars[fname] * vec);
			reals.col(var_map[fname]) = (pars[fname] * vec);
			c++;
		}
	}
	if (c > 0)
		cout << "filled " << c << " fixed pars not listed in user-supplied par ensemble with `parval1` values from control file" << endl;


}

void ParameterEnsemble::save_fixed()
{
	
	if (fixed_names.size() == 0)
		return;
	
	Eigen::MatrixXd fixed_reals = get_eigen(vector<string>(), fixed_names);

	for (int i = 0; i < real_names.size(); i++)
	{
		for (int j = 0; j < fixed_names.size(); j++)
		{
			pair<string, string> key(real_names[i], fixed_names[j]);
			fixed_map[key] = fixed_reals(i, j);

		}
	}
	// add the "base" if its not in the real names already
	if (find(real_names.begin(), real_names.end(), base_name) == real_names.end())
	{
		Parameters pars = pest_scenario_ptr->get_ctl_parameters();
		for (auto fname : fixed_names)
		{
			pair<string, string> key(base_name, fname);
			fixed_map[key] = pars[fname];
		}
	}
}

void ParameterEnsemble::enforce_limits(PerformanceLog* plog, bool enforce_chglim)
{
	
	if (tstat != ParameterEnsemble::transStatus::NUM)
	{
		throw_ensemble_error("pe.enforce_bounds() tstat != NUM not implemented");

	}
	//fancy vector shrinking style enforcement
	if (enforce_chglim)
	{
		Parameters base = pest_scenario_ptr->get_ctl_parameters();
		for (int i = 0; i < reals.rows(); i++)
		{

			Parameters real(var_names, reals.row(i));

			pest_scenario_ptr->enforce_par_limits(plog, real, base, true, true);
			//pest_scenario_ptr->enforce_par_limits(real, base, true, false);
			reals.row(i) = real.get_data_eigen_vec(var_names);
		}
	}
	//reset parameters to be inbounds - very crude
	else
	{

		double l, u, v;
		if (pest_scenario_ptr->get_pestpp_options().get_enforce_tied_bounds())

		{
			Parameters real = pest_scenario_ptr->get_ctl_parameters();//.get_subset(var_names.begin(), var_names.end());
			//pest_scenario_ptr->get_base_par_tran_seq().ctl2numeric_ip(real);
			pair<Parameters, Parameters> ppar;
			ppar = pest_scenario_ptr->get_effective_ctl_lower_upper_bnd(real);
			for (int i = 0; i < reals.rows(); i++)
			{
				real.update_without_clear(var_names, reals.row(i));
				Parameters real_ctl = pest_scenario_ptr->get_base_par_tran_seq().numeric2ctl_cp(real);

				//cout << "";
				for (auto n : var_names)
				{
					v = real_ctl[n];
					l = ppar.first[n];
					u = ppar.second[n];
					v = v < l ? l : v;
					v = v > u ? u : v;
					real_ctl.update_rec(n, v);
				}
				reals.row(i) = pest_scenario_ptr->get_base_par_tran_seq().ctl2numeric_cp(real_ctl).get_data_eigen_vec(var_names);

			}
		}
		else
		{
			//this works except maybe not when tied pars are present
			ParameterInfo pinfo = pest_scenario_ptr->get_ctl_parameter_info();
			Parameters lower = pest_scenario_ptr->get_ctl_parameter_info().get_low_bnd(var_names);
			Parameters upper = pest_scenario_ptr->get_ctl_parameter_info().get_up_bnd(var_names);
			par_transform.ctl2numeric_ip(lower);
			par_transform.ctl2numeric_ip(upper);
			Eigen::VectorXd col;

			for (int j = 0; j < reals.cols(); j++)
			{
				l = lower[var_names[j]];
				u = upper[var_names[j]];
				col = reals.col(j);
				for (int i = 0; i < reals.rows(); i++)
				{
					v = col(i);
					v = v < l ? l : v;
					col(i) = v > u ? u : v;
				}
				reals.col(j) = col;
			}
		}
		

	}
}

void ParameterEnsemble::to_binary(string file_name)
{


	ofstream fout(file_name, ios::binary);
	if (!fout.good())
	{
		throw runtime_error("error opening file for binary parameter ensemble:" + file_name);
	}

	//vector<string> vnames = var_names;
	vector<string> vnames = pest_scenario_ptr->get_ctl_ordered_par_names();
	//vnames.insert(vnames.end(), fixed_names.begin(), fixed_names.end());
	ParameterInfo pi = pest_scenario_ptr->get_ctl_parameter_info();
	ParameterRec::TRAN_TYPE ft = ParameterRec::TRAN_TYPE::FIXED;
	ParameterRec::TRAN_TYPE tt = ParameterRec::TRAN_TYPE::TIED;

	vector<string> f_names;
	for (auto &name : pest_scenario_ptr->get_ctl_ordered_par_names())
	{
		if ((pi.get_parameter_rec_ptr(name)->tranform_type == ft) || 
			(pi.get_parameter_rec_ptr(name)->tranform_type == tt))
		{
			f_names.push_back(name);
		}
	}
	//vnames.insert(vnames.end(),f_names.begin(),f_names.end());
	

	int n_var = vnames.size();
	int n_real = real_names.size();
	int n;
	int tmp;
	double data;
	char par_name[200];
	char obs_name[200];

	// write header
	tmp = n_var;
	fout.write((char*)&tmp, sizeof(tmp));
	tmp = n_real;
	fout.write((char*)&tmp, sizeof(tmp));

	//write number nonzero elements in jacobian (includes prior information)

	n = reals.size() + (f_names.size() * shape().first);
	fout.write((char*)&n, sizeof(n));

	//write matrix
	n = 0;
	//map<string, double>::const_iterator found_pi_par;
	//map<string, double>::const_iterator not_found_pi_par;
	//icount = row_idxs + 1 + col_idxs * self.shape[0]
	Parameters pars;
	for (int irow = 0; irow<n_real; ++irow)
	{
		pars.update_without_clear(var_names, reals.row(irow));
		if (tstat == transStatus::MODEL)
			par_transform.model2ctl_ip(pars);
		else if (tstat == transStatus::NUM)
			par_transform.numeric2ctl_ip(pars);
		replace_fixed(real_names[irow], pars);

		for (int jcol = 0; jcol<n_var; ++jcol)
		{
			n = irow;
			fout.write((char*) &(n), sizeof(n));
			n = jcol;
			fout.write((char*) &(n), sizeof(n));
			data = pars[vnames[jcol]];
			fout.write((char*) &(data), sizeof(data));
		}
		/*int jcol = n_var;
		for (auto &fname : fixed_names)
		{
			n = irow + 1 + (jcol * n_real);
			data = fixed_map.at(pair<string, string>(real_names[irow], fname));
			fout.write((char*) &(n), sizeof(n));
			fout.write((char*) &(data), sizeof(data));
			jcol++;
		}*/
	}
	//save parameter names
	for (vector<string>::const_iterator b = vnames.begin(), e = vnames.end();
		b != e; ++b) {
		string l = pest_utils::lower_cp(*b);
		pest_utils::string_to_fortran_char(l, par_name, 200);
		fout.write(par_name, 200);
	}
	/*for (auto fname : fixed_names)
	{
		pest_utils::string_to_fortran_char(pest_utils::lower_cp(fname), par_name, 12);
		fout.write(par_name, 12);
	}*/

	//save observation and Prior information names
	for (vector<string>::const_iterator b = real_names.begin(), e = real_names.end();
		b != e; ++b) {
		string l = pest_utils::lower_cp(*b);
		pest_utils::string_to_fortran_char(l, obs_name, 200);
		fout.write(obs_name, 200);
	}
	//save observation names (part 2 prior information)
	fout.close();
}


void ParameterEnsemble::to_csv(string file_name)
{
	//write the par ensemble to csv file - transformed back to CTL status
	
	ofstream csv(file_name);
	if (!csv.good())
	{
		throw_ensemble_error("ParameterEnsemble.to_csv() error opening csv file " + file_name + " for writing");
	}
	if (pest_scenario_ptr->get_pestpp_options().get_ies_csv_by_reals())
		to_csv_by_reals(csv);
	else
		to_csv_by_vars(csv);
}


void ParameterEnsemble::to_csv_by_vars(ofstream &csv)
{
	vector<string> names = pest_scenario_ptr->get_ctl_ordered_par_names();
	csv << "var_name";
	int ireal = 0;
	map<string, int> real_map;
	for (int i = 0; i < real_names.size(); i++)
		real_map[real_names[i]] = i;
	map<string, int>::iterator end = real_map.end();
	for (auto rname : org_real_names)
	{
		if (real_map.find(rname) == end)
			continue;
		ireal = real_map[rname];
		csv << ',' << pest_utils::lower_cp(real_names[ireal]);
	}
	csv << endl;

	//get the pars and transform to be in sync with ensemble trans status
	Parameters pars = pest_scenario_ptr->get_ctl_parameters();
	if (tstat == transStatus::NUM)
	{
		par_transform.active_ctl2numeric_ip(pars);
	}
	else if (tstat == transStatus::MODEL)
	{
		par_transform.active_ctl2model_ip(pars);
	}

	//for (int ireal = 0; ireal < reals.rows(); ireal++)
	map<string, Parameters> tpar_map;
	for (auto rname : org_real_names)
	{
		if (real_map.find(rname) == end)
			continue;
		ireal = real_map[rname];
		pars.update_without_clear(var_names, reals.row(ireal));
		if (tstat == transStatus::MODEL)
			par_transform.model2ctl_ip(pars);
		else if (tstat == transStatus::NUM)
			par_transform.numeric2ctl_ip(pars);
		replace_fixed(real_names[ireal], pars);
		tpar_map[rname] = pars;
	}
	for (auto &name : names)
	{
		csv << pest_utils::lower_cp(name);
		for (auto rname : org_real_names)
		{
			if (real_map.find(rname) == end)
				continue;
			Parameters *tpars = &tpar_map[rname];
			csv << ',' << tpars->get_rec(name);
		}
		csv << endl;
	}
}

void ParameterEnsemble::to_csv_by_reals(ofstream &csv)
{
	vector<string> names = pest_scenario_ptr->get_ctl_ordered_par_names();
	csv << "real_name";
	for (auto &vname : names)
		csv << ',' << pest_utils::lower_cp(vname);
	csv << endl;

	//get the pars and transform to be in sync with ensemble trans status
	Parameters pars = pest_scenario_ptr->get_ctl_parameters();
	if (tstat == transStatus::NUM)
	{
		par_transform.active_ctl2numeric_ip(pars);
	}
	else if (tstat == transStatus::MODEL)
	{
		par_transform.active_ctl2model_ip(pars);
	}

	//for (int ireal = 0; ireal < reals.rows(); ireal++)
	int ireal = 0;
	map<string, int> real_map;
	for (int i = 0; i < real_names.size(); i++)
		real_map[real_names[i]] = i;
	map<string, int>::iterator end = real_map.end();
	vector<string> rnames = org_real_names;
	if (rnames.size() == 0)
		rnames = real_names;
	for (auto rname : rnames)
	{
		if (real_map.find(rname) == end)
			continue;
		ireal = real_map[rname];
		csv << pest_utils::lower_cp(real_names[ireal]);
		//evec = reals.row(ireal);
		//svec.assign(evec.data(), evec.data() + evec.size());
		pars.update_without_clear(var_names, reals.row(ireal));
		if (tstat == transStatus::MODEL)
			par_transform.model2ctl_ip(pars);
		else if (tstat == transStatus::NUM)
			par_transform.numeric2ctl_ip(pars);
		replace_fixed(real_names[ireal],pars);
		for (auto &name : names)
			csv << ',' << pars[name];
		csv << endl;
	}
}

void ParameterEnsemble::replace_fixed(string real_name,Parameters &pars)
{
	if (fixed_names.size() > 0)
	{
		for (auto fname : fixed_names)
		{
			pair<string, string> key(real_name, fname);
			double val = fixed_map.at(key);
			pars.update_rec(fname, val);
		}

	}

}

void ParameterEnsemble::transform_ip(transStatus to_tstat)
{
	//transform the ensemble in place
	if (to_tstat == tstat)
		return;
	if ((to_tstat == transStatus::NUM) && (tstat == transStatus::CTL))
	{
		Parameters pars = pest_scenario_ptr->get_ctl_parameters();
		vector<string> adj_par_names = pest_scenario_ptr->get_ctl_ordered_adj_par_names();
		Eigen::MatrixXd new_reals = Eigen::MatrixXd(shape().first, adj_par_names.size());
		for (int ireal = 0; ireal < reals.rows(); ireal++)
		{
			pars.update_without_clear(var_names, reals.row(ireal));
			par_transform.ctl2numeric_ip(pars);
			assert(pars.size() == adj_par_names.size());
			new_reals.row(ireal) = pars.get_data_eigen_vec(adj_par_names);
		}
		reals = new_reals;
		var_names = adj_par_names;
		tstat = to_tstat;

	}
	else
		throw_ensemble_error("ParameterEnsemble::transform_ip() only CTL to NUM implemented");

}
Covariance ParameterEnsemble::get_diagonal_cov_matrix()
{
	transform_ip(transStatus::NUM);
	return Ensemble::get_diagonal_cov_matrix();
}

ObservationEnsemble::ObservationEnsemble(Pest *_pest_scenario_ptr, std::mt19937* _rand_gen_ptr): Ensemble(_pest_scenario_ptr, _rand_gen_ptr)
{
}

ObservationEnsemble::ObservationEnsemble(Pest *_pest_scenario_ptr, std::mt19937* _rand_gen_ptr, Eigen::MatrixXd _reals,
	vector<string> _real_names, vector<string> _var_names) : Ensemble(_pest_scenario_ptr, _rand_gen_ptr)
{
	if (_reals.rows() != _real_names.size())
		throw_ensemble_error("ParameterEnsemble() _reals.rows() != _real_names.size()");
	if (_reals.cols() != _var_names.size())
		throw_ensemble_error("ParameterEnsemble() _reals.cols() != _var_names.size()");
	reals = _reals;
	var_names = _var_names;
	real_names = _real_names;
}


void ObservationEnsemble::initialize_without_noise(int num_reals)
{
	var_names = pest_scenario_ptr->get_ctl_ordered_obs_names();
	Observations obs = pest_scenario_ptr->get_ctl_observations();
	Eigen::VectorXd obsvals = obs.get_data_eigen_vec(var_names);
	reals.resize(num_reals, var_names.size());
	for (int i = 0; i < num_reals; i++)
	{
		reals.row(i) = obsvals;
	}
	real_names.clear();
	stringstream ss;
	for (int i = 0; i < num_reals; i++)
	{
		ss.str("");
		ss << i;
		real_names.push_back(ss.str());
	}
	org_real_names = real_names;

}

void ObservationEnsemble::draw(int num_reals, Covariance &cov, PerformanceLog *plog, int level, ofstream& frec)
{
	//draw an obs ensemble using only nz obs names
	var_names = pest_scenario_ptr->get_ctl_ordered_nz_obs_names();
	Observations obs = pest_scenario_ptr->get_ctl_observations();
	ObservationInfo oi = pest_scenario_ptr->get_ctl_observation_info();
	map<string, vector<string>> grouper;
	vector<string> ogroups = pest_scenario_ptr->get_ctl_ordered_obs_group_names();
	if (pest_scenario_ptr->get_pestpp_options().get_ies_group_draws())
	{
		for (auto group : ogroups)
		{
			for (auto &oname : var_names)
			{
				if (oi.get_group(oname) == group)
					grouper[group].push_back(oname);
			}
		}
	}
	Ensemble::draw(num_reals, cov, obs, pest_scenario_ptr->get_ctl_ordered_nz_obs_names(), grouper, plog, level);

	//apply any bounds that were supplied
	map<string, double> lower_bnd = pest_scenario_ptr->get_ext_file_double_map("observation data external", "lower_bound");
	map<string, double> upper_bnd = pest_scenario_ptr->get_ext_file_double_map("observation data external", "upper_bound");
	set<string> snames(var_names.begin(), var_names.end());
	double v, lb, ub;
	Eigen::VectorXd col;
	string var_name;
	if ((lower_bnd.size() > 0) || (upper_bnd.size() > 0))
	{
		frec << "Note: the following observations contain 'lower_bound' and/or 'upper_bound' information that will be" << endl;
		frec << "      enforced on the additive noise realizations: " << endl;
		for (int j = 0; j < reals.cols(); j++)
		{
			var_name = var_names[j];
			lb = std::numeric_limits<double>::lowest();
			ub = 1.79769e+308;//std::numeric_limits<double>::max();
			if (lower_bnd.find(var_name) != lower_bnd.end())
			{
				lb = lower_bnd[var_name];
			}
			if (upper_bnd.find(var_name) != upper_bnd.end())
			{
				ub = upper_bnd[var_name];
			}
			frec << var_name << " " << lb << " " << ub << endl;
			col = reals.col(j);

			for (int i = 0; i < reals.rows(); i++)
			{
				v = col(i);
				v = v < lb ? lb : v;
				col(i) = v > ub ? ub : v;
			}
			reals.col(j) = col;
			
		}
	}
	


	//now fill in all the zero-weighted obs
	Eigen::MatrixXd drawn = reals;
	vector<string> full_var_names = pest_scenario_ptr->get_ctl_ordered_obs_names();
	reals.resize(num_reals, full_var_names.size());
	reals.setConstant(-1.0e30);
	map<string, int> vmap;
	for (int i = 0; i < var_names.size(); i++)
		vmap[var_names[i]] = i;
	map<string, int>::iterator end = vmap.end();
	//for (auto n : full_var_names)
	string n;
	double val;
	for (int i = 0; i < full_var_names.size(); i++)
	{
		n = full_var_names[i];
		if (vmap.find(n) != end)
			reals.col(i) = drawn.col(vmap[n]);
		else
		{
			val = obs.get_rec(n);
			reals.col(i).setConstant(val);
		}
	}
	var_names = full_var_names;
	
}

void ObservationEnsemble::update_from_obs(int row_idx, Observations &obs)
{
	//update a row in reals from an int id
	if (row_idx >= real_names.size())
		throw_ensemble_error("ObservtionEnsemble.update_from_obs() obs_idx out of range");
	//Eigen::VectorXd temp = obs.get_data_eigen_vec(var_names);
	//cout << reals.row(row_idx) << endl << endl;
	//cout << temp << endl;
	reals.row(row_idx) = obs.get_data_eigen_vec(var_names);
	//cout << reals.row(row_idx) << endl;
	return;
}

void ObservationEnsemble::update_from_obs(string real_name, Observations &obs)
{
	//update a row in reals from a string id
	vector<string>::iterator real,start = real_names.begin(), end = real_names.end();
	real = find(start, end, real_name);
	if (real == end)
	{
		throw_ensemble_error("real_name not in real_names:" + real_name);
	}
	update_from_obs(real - start, obs);
}

vector<int> ObservationEnsemble::update_from_runs(map<int,int> &real_run_ids, RunManagerAbstract *run_mgr_ptr)
{
	//update the obs ensemble in place from the run manager
	set<int> failed_runs = run_mgr_ptr->get_failed_run_ids();
	vector<int> failed_real_idxs;
	Parameters pars = pest_scenario_ptr->get_ctl_parameters();
	Observations obs = pest_scenario_ptr->get_ctl_observations();
	string real_name;
	int ireal = 0;
	for (auto &real_run_id : real_run_ids)
	{
		if (failed_runs.find(real_run_id.second) != failed_runs.end())
		{
			failed_real_idxs.push_back(real_run_id.first);
			//failed_real_idxs.push_back(ireal);
		}

		else
		{
			run_mgr_ptr->get_run(real_run_id.second, pars, obs);
			//real_name = real_names[real_run_id.first];
			update_from_obs(real_run_id.first, obs);
			//update_from_obs(ireal, obs);
		}
	}
	return failed_real_idxs;
}

void ObservationEnsemble::from_binary(string file_name)
{
	//load obs en from binary jco-type file
	vector<string> names = pest_scenario_ptr->get_ctl_ordered_obs_names();
	Ensemble::from_binary(file_name, names, true);
}

void ObservationEnsemble::from_csv(string file_name)
{
	//load the obs en from a csv file
	var_names = pest_scenario_ptr->get_ctl_ordered_obs_names();
	bool csv_by_reals = pest_scenario_ptr->get_pestpp_options().get_ies_csv_by_reals();
	ifstream csv(file_name);
	if (!csv.good())
		throw runtime_error("error opening observation csv " + file_name + " for reading");
	pair<map<string,int>, map<string, int>> p = prepare_csv(pest_scenario_ptr->get_ctl_ordered_nz_obs_names(), csv, false);

	map<string, int> header_info = p.first, index_info = p.second;
	//blast through the file to get number of reals
	string line;
	int num_reals;
	if (csv_by_reals)
		num_reals = index_info.size();
	else
		num_reals = header_info.size();
	
	csv.close();
	csv.open(file_name);
	if (!csv.good())
		throw runtime_error("error re-opening observation csv " + file_name + " for reading");
	getline(csv, line);
	//Ensemble::read_csv(num_reals, csv,header_info,index_info);
	if (csv_by_reals)
		Ensemble::read_csv_by_reals(num_reals, csv, header_info, index_info);
	else
		Ensemble::read_csv_by_vars(num_reals, csv, header_info, index_info);
}

void ObservationEnsemble::from_eigen_mat(Eigen::MatrixXd mat, const vector<string> &_real_names, const vector<string> &_var_names)
{
	//reset obs en from components
	vector<string> missing;
	vector<string>::const_iterator start = _var_names.begin();
	vector<string>::const_iterator end = _var_names.end();

	for (auto &name : pest_scenario_ptr->get_ctl_ordered_obs_names())
		if (find(start, end, name) == end)
			missing.push_back(name);
	if (missing.size() > 0)
		throw_ensemble_error("ObservationEnsemble.from_eigen_mat() the following obs names not found: ", missing);
	Ensemble::from_eigen_mat(mat, _real_names, _var_names);
}

DrawThread::DrawThread(PerformanceLog * _performance_log, Covariance & _cov,
	Eigen::MatrixXd *_draws_ptr, vector<string> &_group_keys, const map<string, vector<string>> &_grouper) : cov(_cov),
	group_keys(_group_keys), grouper(_grouper)
{
	//idx_map = _idx_map;
	//std_map = _std_map;
	//std_map = _std_map;
	performance_log = _performance_log;
	draws_ptr = _draws_ptr;
}


void DrawThread::work(int thread_id, int num_reals, int ies_verbose, map<string, int> idx_map, map<string,double> std_map)
{
	stringstream ss;

	unique_lock<mutex> cov_guard(cov_lock,defer_lock);
	unique_lock<mutex> draw_guard(draw_lock,defer_lock);
	//unique_lock<mutex> grouper_guard(grouper_lock,defer_lock);
	unique_lock<mutex> pfm_guard(pfm_lock,defer_lock);
	unique_lock<mutex> key_guard(key_lock, defer_lock);
	int count = 0;
	string group;
	vector<string> names;
	Covariance gcov;
	vector<int> idx;
	Eigen::MatrixXd block, proj;
	RedSVD::RedSymEigen<Eigen::SparseMatrix<double>> eig;
	while (true)
	{
		//get the key for more work, or return if all work done
		while (true)
		{
			if (key_guard.try_lock())
			{
				if (group_keys.size() == 0)
				{
					ss.str("");
					ss << "draw thread: " << thread_id << " processed " << count << " groups";
					if (ies_verbose > 1)
					{
						cout << ss.str() << endl;

					}
					while (true)
					{
						if (pfm_guard.try_lock())
						{
							performance_log->log_event(ss.str());
							pfm_guard.unlock();
							break;
						}
					}
					key_guard.unlock();
					return;
				}
			}
			group = group_keys[group_keys.size() - 1];
			group_keys.pop_back();
			names = grouper[group];
			if (names.size() == 0)
			{
				ss.str("");
				ss << "no entries for grouper key:" << group;
				while (true)
				{
					if (pfm_guard.try_lock())
					{
						performance_log->log_event(ss.str());
						pfm_guard.unlock();
						break;
					}
				}
				key_guard.unlock();
				continue;
			}
			key_guard.unlock();
			break;
		}
		
		if (ies_verbose > 1)
		{
			ss.str("");
			ss << "...processing " << group << " with " <<names.size() << " elements" << endl;
			cout << ss.str();
			while (true)
			{
				if (pfm_guard.try_lock())
				{
					performance_log->log_event(ss.str());
					pfm_guard.unlock();
					break;
				}
			}
		}

		//if there is only one par in the group
		if (names.size() == 1)
		{
			ss.str("");
			ss << "thread: " << thread_id << " - only one element in group " <<group << ", scaling by std";
			while (true)
			{
				if (pfm_guard.try_lock())
				{
					performance_log->log_event(ss.str());
					pfm_guard.unlock();
					break;
				}
			}
			int j = idx_map[names[0]];
			while (true)
			{
				if (draw_guard.try_lock())
				{
					draws_ptr->col(j) *= std_map[names[0]];
					draw_guard.unlock();
					break;
				}
			}
			
			continue;
		}

		//get a sub cov
		while (true)
		{
			if (cov_guard.try_lock())
			{
				gcov = cov.get(names);
				cov_guard.unlock();
				break;
			}
		}
		
		if (ies_verbose > 2)
		{
			gcov.to_ascii(group + "_cov.dat");
		}

		idx.clear();
		for (auto n : names)
			idx.push_back(idx_map[n]);

		if (idx.size() != idx[idx.size() - 1] - idx[0] + 1)
		{
			ss.str("");
			ss << "thread: " << thread_id << " - DrawThread error: idx out of order for group: " << group;
			while (true)
			{
				if (pfm_guard.try_lock())
				{
					performance_log->log_event(ss.str());
					pfm_guard.unlock();
					break;
				}
			}
			throw runtime_error(ss.str());
		}


		double fac = gcov.e_ptr()->diagonal().minCoeff();
		ss.str("");
		ss << "thread: " << thread_id <<  " - min variance for group " << group << ": " << fac;
		while (true)
		{
			if (pfm_guard.try_lock())
			{
				performance_log->log_event(ss.str());
				pfm_guard.unlock();
				break;
			}
		}
		/*while (true)
		{
			if (draw_guard.try_lock())
			{
				block = draws_ptr->block(0, idx[0], num_reals, idx.size());
				draw_guard.unlock();
				break;
			}
		}*/
		
		ss.str("");
		ss << "thread: " << thread_id <<  " - Randomized Eigen decomposition of full cov for " << names.size() << " element matrix" << endl;
		while (true)
		{
			if (pfm_guard.try_lock())
			{
				performance_log->log_event(ss.str());
				pfm_guard.unlock();
				break;
			}
		}
		eig.compute(*gcov.e_ptr() * (1.0 / fac), names.size());
		//RedSVD::RedSVD<Eigen::SparseMatrix<double>> svd;
		//svd.compute(*gcov.e_ptr(),gcov.get_col_names().size());// , gi.second.size());
		//cout << svd.singularValues() << endl;
		//Eigen::JacobiSVD<Eigen::MatrixXd> svd(gcov.e_ptr()->toDense(), Eigen::ComputeFullU);
		//svd.computeU();

		proj = (eig.eigenvectors() * (fac *eig.eigenvalues()).cwiseSqrt().asDiagonal());
		//proj = (svd.matrixU() * svd.singularValues().asDiagonal());

		if (ies_verbose > 2)
		{
			ofstream f(group + "_evec.dat");
			f << eig.eigenvectors() << endl;
			//f << svd.matrixU() << endl;
			f.close();
			ofstream ff(group + "_sqrt_evals.dat");
			ff << (fac * eig.eigenvalues()).cwiseSqrt() << endl;
			//ff << svd.singularValues() << endl;
			ff.close();
			ofstream fff(group + "_proj.dat");
			fff << proj << endl;
			fff.close();
		}
		//cout << "block " << block.rows() << " , " << block.cols() << endl;
		//cout << " proj " << proj.rows() << " , " << proj.cols() << endl;
		while (true)
		{
			if (draw_guard.try_lock())
			{
				block = draws_ptr->block(0, idx[0], num_reals, idx.size());
				draws_ptr->block(0, idx[0], num_reals, idx.size()) = (proj * block.transpose()).transpose();
				draw_guard.unlock();
				break;
			}
		}
		

	}

}
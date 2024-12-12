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

#include "SVDPackage.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "RedSVD-h.h"


using namespace Eigen;

SVDPackage::SVDPackage(std::string _descritpion, int _n_max_sing, double _eign_thres) : description(_descritpion), n_max_sing(_n_max_sing), eign_thres(_eign_thres)
{
	performance_log = nullptr;
}


void SVDPackage::set_performance_log(PerformanceLog *_performance_log)
{
	performance_log = _performance_log;
}

void SVDPackage::set_max_sing(int _n_max_sing) {
	n_max_sing = _n_max_sing;
}

void SVDPackage::set_eign_thres(double _eign_thres)
{
	eign_thres = _eign_thres;
}

int SVDPackage::get_max_sing() {
	return n_max_sing;
}

double SVDPackage::get_eign_thres()
{
	return eign_thres;
}

void SVD_REDSVD::solve_ip(Eigen::SparseMatrix<double>& A, Eigen::VectorXd &Sigma, Eigen::SparseMatrix<double>& U,
	Eigen::SparseMatrix<double>& Vt, Eigen::VectorXd &Sigma_trunc)
{
	solve_ip(A, Sigma, U, Vt, Sigma_trunc, eign_thres);
}

void SVD_REDSVD::solve_ip(Eigen::MatrixXd& A, Eigen::MatrixXd &Sigma, Eigen::MatrixXd& U,
	Eigen::MatrixXd& V, double _eigen_thres, int _max_sing)
{
	if (performance_log)
		performance_log->log_event("starting REDSVD");

	RedSVD::RedSVD<MatrixXd> red_svd(A);
	
	U = red_svd.matrixU();
	V = red_svd.matrixV();
	VectorXd Sigma_full = red_svd.singularValues();

	int kmax = (Sigma_full.size() < n_max_sing) ? Sigma_full.size() : n_max_sing;
	int num_sing_used = 0;
	double eig_ratio;
	for (int i_sing = 0; i_sing < kmax; ++i_sing)
	{
		eig_ratio = Sigma_full[i_sing] / Sigma_full[0];
		if ((eig_ratio > _eigen_thres) && (i_sing <= _max_sing))
		{
			++num_sing_used;
		}
		else
		{
			break;
		}
	}
	std::stringstream ss;
	ss << "trimming REDSVD components to " << num_sing_used << "elements";
	if (performance_log)
		performance_log->log_event(ss.str());
	//std::cout << Sigma_full << std::endl;
	Sigma = Sigma_full.head(num_sing_used);
	//std::cout << Sigma << std::endl;
	//Sigma_trunc = Sigma_full.tail(Sigma_full.size() - num_sing_used);
	if (num_sing_used == 1)
	{
		Eigen::VectorXd temp = V.col(0);
		V.resize(V.rows(),1);
		V.col(0) = temp;
		temp = U.col(0);
		U.resize(U.rows(),1);
		U.col(0) = temp;
	}
	else
	{
		//std::cout << U.col(0);
		Eigen::MatrixXd temp = V.leftCols(num_sing_used);
		V = temp;
		temp = U.leftCols(num_sing_used);
		U = temp;
		//std::cout << U.col(0);


	}
	if (performance_log)
		performance_log->log_event("done REDSVD");
	return;
}

void SVD_REDSVD::solve_ip(Eigen::MatrixXd& A, Eigen::MatrixXd &Sigma, Eigen::MatrixXd& U,
	Eigen::MatrixXd& V)
{
	if (performance_log)
		performance_log->log_event("starting REDSVD");

	RedSVD::RedSVD<MatrixXd> red_svd(A);
	if (performance_log)
		performance_log->log_event("retrieving REDSVD components");
	U = red_svd.matrixU();
	V = red_svd.matrixV();
	Sigma = red_svd.singularValues();
	if (performance_log)
		performance_log->log_event("done REDSVD");
	return;
}


void SVD_REDSVD::solve_ip(Eigen::SparseMatrix<double>& A, Eigen::VectorXd &Sigma, Eigen::SparseMatrix<double>& U,
	Eigen::SparseMatrix<double>& VT, Eigen::VectorXd &Sigma_trunc, double _eigen_thres)
{
	if (performance_log)
		performance_log->log_event("starting REDSVD");

	RedSVD::RedSVD<MatrixXd> red_svd(A,n_max_sing);
	if (performance_log)
		performance_log->log_event("retrieving REDSVD components");
	U = red_svd.matrixU().sparseView();
	VT = red_svd.matrixV().transpose().sparseView();
	VectorXd Sigma_full = red_svd.singularValues();

	int kmax = (Sigma_full.size() < n_max_sing) ? Sigma_full.size() : n_max_sing;
	int num_sing_used = 0;
	double eig_ratio;
	for (int i_sing = 0; i_sing < kmax; ++i_sing)
	{
		eig_ratio = Sigma_full[i_sing] / Sigma_full[0];
		if (eig_ratio > _eigen_thres)
		{
			++num_sing_used;
		}
		else
		{
			break;
		}
	}
	std::stringstream ss;
	ss << "trimming REDSVD components to " << num_sing_used << "elements";
	if (performance_log)
		performance_log->log_event(ss.str());

	Sigma = Sigma_full.head(num_sing_used);
	Sigma_trunc = Sigma_full.tail(Sigma_full.size() - num_sing_used);

	VT = Eigen::SparseMatrix<double>(VT.topRows(num_sing_used));
	U = Eigen::SparseMatrix<double>(U.leftCols(num_sing_used));
	if (performance_log)
		performance_log->log_event("done REDSVD");
	return;

}


void SVD_EIGEN::solve_ip(Eigen::SparseMatrix<double>& A, Eigen::VectorXd &Sigma, Eigen::SparseMatrix<double>& U,
	Eigen::SparseMatrix<double>& VT, Eigen::VectorXd &Sigma_trunc)
{
	solve_ip(A, Sigma, U, VT, Sigma_trunc, eign_thres);
}


void SVD_EIGEN::solve_ip(Eigen::MatrixXd& A, Eigen::MatrixXd &Sigma, Eigen::MatrixXd& U,
	Eigen::MatrixXd& V, double _eigen_thres, double _max_sing)
{

	JacobiSVD<MatrixXd> svd_fac(A, ComputeFullU | ComputeFullV);

	performance_log->log_event("starting REDSVD");

	//RedSVD::RedSVD<MatrixXd> red_svd(A);
	performance_log->log_event("retrieving REDSVD components");
	U = svd_fac.matrixU();
	V = svd_fac.matrixV();
	VectorXd Sigma_full = svd_fac.singularValues();

	int kmax = (Sigma_full.size() < n_max_sing) ? Sigma_full.size() : n_max_sing;
	int num_sing_used = 0;
	double eig_ratio;
	for (int i_sing = 0; i_sing < kmax; ++i_sing)
	{
		eig_ratio = Sigma_full[i_sing] / Sigma_full[0];
		if ((eig_ratio > _eigen_thres) && (i_sing <= _max_sing))
		{
			++num_sing_used;
		}
		else
		{
			break;
		}
	}
	std::stringstream ss;
	ss << "trimming REDSVD components to " << num_sing_used << "elements";
	performance_log->log_event(ss.str());
	//std::cout << Sigma_full << std::endl;
	Sigma = Sigma_full.head(num_sing_used);
	//std::cout << Sigma << std::endl;
	//Sigma_trunc = Sigma_full.tail(Sigma_full.size() - num_sing_used);
	if (num_sing_used == 1)
	{
		Eigen::VectorXd temp = V.col(0);
		V.resize(V.rows(), 1);
		V.col(0) = temp;
		temp = U.col(0);
		U.resize(U.rows(), 1);
		U.col(0) = temp;
	}
	else
	{
		//std::cout << U.col(0);
		Eigen::MatrixXd temp = V.leftCols(num_sing_used);
		V = temp;
		temp = U.leftCols(num_sing_used);
		U = temp;
		//std::cout << U.col(0);


	}
	performance_log->log_event("done REDSVD");
	return;
}

void SVD_EIGEN::solve_ip(Eigen::SparseMatrix<double>& A, Eigen::VectorXd &Sigma, Eigen::SparseMatrix<double>& U,
	Eigen::SparseMatrix<double>& Vt, Eigen::VectorXd &Sigma_trunc,  double _eigen_thres)
{

	JacobiSVD<MatrixXd> svd_fac(A,  ComputeThinU |  ComputeThinV);
	VectorXd Sigma_full = svd_fac.singularValues();
	U = svd_fac.matrixU().sparseView();
	Vt = svd_fac.matrixV().transpose().sparseView();

	//Compute number of singular values to be used in the solution
	int num_sing_used = 0;
	double eig_ratio;

	int kmax = (Sigma_full.size() < n_max_sing) ? Sigma_full.size() : n_max_sing;
	for (int i_sing = 0; i_sing < kmax; ++i_sing)
	{
		eig_ratio = Sigma_full[i_sing] / Sigma_full[0];
		if (eig_ratio > _eigen_thres)
		{
			++num_sing_used;
		}
		else
		{
			break;
		}
	}
	//Trim the Matrices based on the number of singular values to be used
	Sigma = Sigma_full.head(num_sing_used);
	Sigma_trunc = Sigma_full.tail(Sigma_full.size() - num_sing_used);

	Vt = Eigen::SparseMatrix<double>(Vt.topRows(num_sing_used));
	U = Eigen::SparseMatrix<double>(U.leftCols(num_sing_used));
}

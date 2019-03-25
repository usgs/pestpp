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

#ifndef SVDPACKAGE_H_
#define SVDPACKAGE_H_
#include <string>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include "PerformanceLog.h"

class SVDPackage
{
public:
	SVDPackage(std::string _descritpion="undefined", int _n_max_sing=1000, double _eign_thres=1.0e-7);
	SVDPackage(const SVDPackage &rhs) : description(rhs.description), n_max_sing(rhs.n_max_sing), eign_thres(rhs.eign_thres), performance_log(rhs.performance_log) {}
	virtual void solve_ip(Eigen::SparseMatrix<double>& A, Eigen::VectorXd &Sigma, Eigen::SparseMatrix<double> & U, Eigen::SparseMatrix<double>& VT, Eigen::VectorXd &Sigma_trunc) = 0;
	virtual void solve_ip(Eigen::SparseMatrix<double>& A, Eigen::VectorXd &Sigma, Eigen::SparseMatrix<double> & U, Eigen::SparseMatrix<double>& VT, Eigen::VectorXd &Sigma_trunc, double _eigen_thres) = 0;
	virtual void set_max_sing(int _n_max_sing);
	virtual int get_max_sing();
	virtual void set_eign_thres(double _eign_thres);
	virtual double get_eign_thres();
	virtual void set_performance_log(PerformanceLog *_performance_log);
	virtual ~SVDPackage(void){};
	const std::string description;
	virtual SVDPackage *clone() const {return 0;}
protected:
	int n_max_sing;
	double eign_thres;
	PerformanceLog *performance_log;
};

class SVD_EIGEN : public SVDPackage
{
public:
	SVD_EIGEN(int _n_max_sing = 1000, double _eign_thres = 1.0e-7) : SVDPackage("Eigen JacobiSVD", _n_max_sing, _eign_thres)  {}
	SVD_EIGEN(const SVD_EIGEN &rhs) : SVDPackage(rhs) {}
	virtual void solve_ip(Eigen::SparseMatrix<double>& A, Eigen::VectorXd &Sigma, Eigen::SparseMatrix<double>& U, Eigen::SparseMatrix<double>& VT, Eigen::VectorXd &Sigma_trunc);
	virtual void solve_ip(Eigen::SparseMatrix<double>& A, Eigen::VectorXd &Sigma, Eigen::SparseMatrix<double>& U, Eigen::SparseMatrix<double>& VT, Eigen::VectorXd &Sigma_trunc, double _eigen_thres);
	virtual SVD_EIGEN *clone() const{ return new SVD_EIGEN(*this); }
	virtual void solve_ip(Eigen::MatrixXd& A, Eigen::MatrixXd &Sigma, Eigen::MatrixXd& U, Eigen::MatrixXd& VT, double _eigen_thres, double _max_sing);

	virtual ~SVD_EIGEN(void) {}
};


class SVD_REDSVD : public SVDPackage
{
public:
	SVD_REDSVD(int _n_max_sing = 1000, double _eign_thres = 1.0e-7) : SVDPackage("RedSVD", _n_max_sing, _eign_thres) {}
	SVD_REDSVD(const SVD_REDSVD &rhs) : SVDPackage(rhs) {}
	virtual void solve_ip(Eigen::SparseMatrix<double>& A, Eigen::VectorXd &Sigma, Eigen::SparseMatrix<double>& U,
		Eigen::SparseMatrix<double>& VT, Eigen::VectorXd &Sigma_trunc);
	virtual void solve_ip(Eigen::SparseMatrix<double>& A, Eigen::VectorXd &Sigma, Eigen::SparseMatrix<double>& U,
		Eigen::SparseMatrix<double>& VT, Eigen::VectorXd &Sigma_trunc, double _eigen_thres);
	virtual void solve_ip(Eigen::MatrixXd& A, Eigen::MatrixXd &Sigma, Eigen::MatrixXd& U,
		Eigen::MatrixXd& V, double _eigen_thres, int _max_sing);
	virtual SVD_REDSVD *clone() const { return new SVD_REDSVD(*this); }
	virtual void solve_ip(Eigen::MatrixXd& A, Eigen::MatrixXd &Sigma, Eigen::MatrixXd& U,
		Eigen::MatrixXd& V);
	virtual ~SVD_REDSVD(void) {}
};

#endif //SVDPACKAGE_H_

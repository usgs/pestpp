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
#ifndef SVD_PROPACK_H_
#define SVD_PROPACK_H_

#include "SVDPackage.h"
#include<Eigen/Dense>
#include<Eigen/Sparse>

class SVD_PROPACK  : public SVDPackage
{
public:
	SVD_PROPACK(int _n_max_sing = 1000, double _eign_thres = 1.0e-7);
	void solve_ip(Eigen::SparseMatrix<double>& A, Eigen::VectorXd &Sigma, Eigen::SparseMatrix<double>& U, Eigen::SparseMatrix<double>& Vt, Eigen::VectorXd &Sigma_trunc);
	void solve_ip(Eigen::SparseMatrix<double>& A, Eigen::VectorXd &Sigma, Eigen::SparseMatrix<double> &U, Eigen::SparseMatrix<double>& Vt, Eigen::VectorXd &Sigma_trunc, double _eigen_thres);
	void solve_ip(Eigen::MatrixXd& A, Eigen::MatrixXd &Sigma, Eigen::MatrixXd& U, Eigen::MatrixXd& V, double _eigen_thres, int _max_sing);
	void test();
	~SVD_PROPACK(void);
private:
};

#endif /*SVD_PROPACK_H_*/

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
#include <vector>
#include <iostream>
#include <algorithm>
#include "SVD_PROPACK.h"
#include "config_os.h"

using namespace std;
using namespace Eigen;

extern "C" {
  double DEF_DLAMCH(char*);
//  void DEF_DLANBPRO_SPARCE(int *m, int *n, int *k0, int *k, double *U,
//			   int *ldu, double *V, int *ldv, double *B, int *ldb,
//		   double *rnorm, double *doption, int *ioption, double *work,
//			   int *iwork, double *dparm, int *iparm, int *ierr);
  void DEF_DLANSVD(char *jobu, char *jobv,int *m,int *n,int *k, int *kmax, double *U, int *ldu, double *Sigma, double *bnd,
		   double *V, int *ldv, double *tolin, double *work, int *lwork, int *iwork,int *liwork, double *doption, int *ioption, int *info,
		   double *dparm, int *iparm, long *jobu_len, long *jobv_len);
}

SVD_PROPACK::SVD_PROPACK(int _n_max_sing, double _eign_thres) : SVDPackage("PROPACK", _n_max_sing, _eign_thres)

{
}

void SVD_PROPACK::solve_ip(Eigen::SparseMatrix<double>& A, VectorXd &Sigma, Eigen::SparseMatrix<double> &U, Eigen::SparseMatrix<double>& Vt, VectorXd &Sigma_trunc)
{
	solve_ip(A, Sigma, U, Vt, Sigma_trunc, eign_thres);
}

void SVD_PROPACK::solve_ip(Eigen::SparseMatrix<double>& A, VectorXd &Sigma, Eigen::SparseMatrix<double> &U, Eigen::SparseMatrix<double>& Vt, VectorXd &Sigma_trunc, double _eigen_thres)
{
	class local_utils {
		public:
		static void init_array(double data[], int len, double value = 0.0)
		{
			for(int i=0; i<len; ++i)
			{
				data[i] = value;
			}
		}
	};
	int m_rows = A.rows();
	int n_cols = A.cols();
	int n_nonzero;
	int k = min(m_rows, n_cols);
	k = min(n_max_sing, k);
        int kmax = k;
	int ioption[] = {0, 1};
	char eps_char = 'e';
	double eps = DEF_DLAMCH(&eps_char);
	double doption[] = {0.0, sqrt(eps), pow(eps, 3.0/4.0), 0.0};
	int ierr = 0;

	//count number of nonzero entries
	n_nonzero = A.nonZeros();

	// Allocate and initialize arrays
	double *dparm = new double[n_nonzero];
	int *iparm = new int[2*n_nonzero+1];
	double *tmp_u = new double[m_rows*(kmax+1)];
	double *tmp_sigma = new double[k];
	double *tmp_bnd = new double[k];
	double *tmp_v = new double[n_cols*kmax];
	int nb = 1;
	int lwork =  m_rows + n_cols + 9*kmax + 5*kmax*kmax + 4 + max(3*kmax*kmax+4*kmax+4, nb*max(m_rows,n_cols));
	if (performance_log)
	{
		stringstream info_str;
		info_str << "allocating lwork; size = " << lwork << endl;
		performance_log->log_event(info_str.str());
	}
	double *tmp_work = new double[lwork];
	//cerr << "lwork allocated" << endl;
	int liwork = 8*kmax;
	if (performance_log)
	{
		stringstream info_str;
		info_str << "allocating tmp_iwork; size = " << liwork << endl;
		performance_log->log_event(info_str.str());
	}
	int *tmp_iwork = new int[liwork];
	iparm[0] = n_nonzero;
	int n=0;
	double data;
	for (int icol=0; icol<A.outerSize(); ++icol)
	{
		for (SparseMatrix<double>::InnerIterator it(A, icol); it; ++it)
		{
			data = it.value();
			dparm[n] = data;
			++n;
			iparm[n] = it.row()+1;
			iparm[n_nonzero+n] = it.col()+1;
		}
	}
	local_utils::init_array(tmp_u, m_rows*(kmax+1), 0.0);
	local_utils::init_array(tmp_sigma, k, 0.0);
	local_utils::init_array(tmp_bnd, k, 0.0);;
	local_utils::init_array(tmp_v, n_cols*kmax, 0.0);
	local_utils::init_array(tmp_work,lwork, 0.0);
	char jobu = 'Y';
    char jobv = 'Y';
	long jobu_len = 1;
	long jobv_len = 1;
	int ld_tmpu = m_rows;
	int ld_tmpv = n_cols;
	int ld_tmpb = kmax;

	// Compute singular values and vectors
	int info=0;
	double tolin = 1.0E-4 ;
	if (performance_log)
	{
		performance_log->log_event("calling DEF_DLANSVD");
	}
	DEF_DLANSVD(&jobu, &jobv, &m_rows, &n_cols, &k, &kmax, tmp_u, &ld_tmpu, tmp_sigma, tmp_bnd,
		   tmp_v, &ld_tmpv, &tolin, tmp_work, &lwork, tmp_iwork, &liwork, doption, ioption, &info,
		   dparm,iparm, &jobu_len, &jobv_len);

	if (performance_log)
	{
		performance_log->log_event("updating Vt, S and U matrices");
	}
	//Compute number of singular values to be used in the solution
	int n_sing_used = 0;
	int actual_sing = kmax;
	if (info > 0)
	  actual_sing = info;
	for (int i_sing = 0; i_sing < actual_sing; ++i_sing)
	{
		double eig_ratio = tmp_sigma[i_sing] / tmp_sigma[0];
		if (eig_ratio > _eigen_thres)
		{
			++n_sing_used;
		}
		else
		{
			break;
		}
	}
	// Update Sigma
	Sigma.resize(n_sing_used);
	Sigma.setConstant(0.0);
	for (int i_sing = 0; i_sing<n_sing_used; ++i_sing)
	{
		Sigma(i_sing) = tmp_sigma[i_sing];
	}

	Sigma_trunc.resize(kmax - n_sing_used);
	Sigma_trunc.setConstant(0.0);
	for (int i_sing = n_sing_used; i_sing<kmax; ++i_sing)
	{
		Sigma_trunc(i_sing - n_sing_used) = tmp_sigma[i_sing];
	}

	std::vector<Eigen::Triplet<double> > triplet_list;
	triplet_list.reserve(n_sing_used);
	// Update U
	for (int i_sing = 0; i_sing<n_sing_used; ++i_sing)
	{
		for(int irow=0; irow<m_rows; ++ irow)
		{
			if (tmp_u[i_sing*m_rows+irow] != 0)
			{
				triplet_list.push_back(Eigen::Triplet<double>(irow,i_sing, tmp_u[i_sing*m_rows+irow]));
			}
		}
	}
	U.resize(m_rows, n_sing_used);
	U.setZero();
	U.setFromTriplets(triplet_list.begin(), triplet_list.end());

	triplet_list.clear();
	// Update Vt
	for (int i_sing = 0; i_sing<n_sing_used; ++i_sing)
	{
		for(int icol=0; icol<n_cols; ++ icol)
		{
			if (tmp_v[i_sing*n_cols+icol] != 0)
			{
				triplet_list.push_back(Eigen::Triplet<double>(i_sing, icol, tmp_v[i_sing*n_cols+icol]));
			}
		}
	}
	Vt.resize(n_sing_used, n_cols);
	Vt.setZero();
	Vt.setFromTriplets(triplet_list.begin(), triplet_list.end());

	delete [] dparm;
	delete [] iparm;
	delete [] tmp_u;
	delete [] tmp_sigma;
	delete [] tmp_bnd;
	delete [] tmp_v;
	delete [] tmp_work;
	delete [] tmp_iwork;

}


void SVD_PROPACK::solve_ip(Eigen::MatrixXd& A, Eigen::MatrixXd &Sigma, Eigen::MatrixXd& U, Eigen::MatrixXd& V, double _eigen_thres, int _max_sing)
{
	class local_utils {
	public:
		static void init_array(double data[], int len, double value = 0.0)
		{
			for (int i = 0; i<len; ++i)
			{
				data[i] = value;
			}
		}
	};
	int m_rows = A.rows();
	int n_cols = A.cols();
	int n_nonzero;
	int k = min(m_rows, n_cols);
	k = min(n_max_sing, k);
	int kmax = k;
	int ioption[] = { 0, 1 };
	char eps_char = 'e';
	double eps = DEF_DLAMCH(&eps_char);
	double doption[] = { 0.0, sqrt(eps), pow(eps, 3.0 / 4.0), 0.0 };
	int ierr = 0;

	//count number of nonzero entries
	//n_nonzero = A.nonZeros();
	n_nonzero = m_rows * n_cols;

	// Allocate and initialize arrays
	double *dparm = new double[n_nonzero];
	int *iparm = new int[2 * n_nonzero + 1];
	double *tmp_u = new double[m_rows*(kmax + 1)];
	double *tmp_sigma = new double[k];
	double *tmp_bnd = new double[k];
	double *tmp_v = new double[n_cols*kmax];
	int nb = 1;
	int lwork = m_rows + n_cols + 9 * kmax + 5 * kmax*kmax + 4 + max(3 * kmax*kmax + 4 * kmax + 4, nb*max(m_rows, n_cols));
	if (performance_log)
	{
		stringstream info_str;
		info_str << "allocating lwork; size = " << lwork << endl;
		performance_log->log_event(info_str.str());
	}
	double *tmp_work = new double[lwork];
	//cerr << "lwork allocated" << endl;
	int liwork = 8 * kmax;
	if (performance_log)
	{
		stringstream info_str;
		info_str << "allocating tmp_iwork; size = " << liwork << endl;
		performance_log->log_event(info_str.str());
	}
	int *tmp_iwork = new int[liwork];
	iparm[0] = n_nonzero;
	int n = 0;
	double data;
	//for (int icol = 0; icol<A.outerSize(); ++icol)
	for (int icol = 0; icol < n_cols; ++icol)
	{
		//for (SparseMatrix<double>::InnerIterator it(A, icol); it; ++it)
		for (int jrow=0;jrow< m_rows;++jrow)
		{
			dparm[n] = A(jrow,icol);
			++n;
			iparm[n] = jrow + 1;
			iparm[n_nonzero + n] = icol + 1;
		}
	}
	local_utils::init_array(tmp_u, m_rows*(kmax + 1), 0.0);
	local_utils::init_array(tmp_sigma, k, 0.0);
	local_utils::init_array(tmp_bnd, k, 0.0);;
	local_utils::init_array(tmp_v, n_cols*kmax, 0.0);
	local_utils::init_array(tmp_work, lwork, 0.0);
	char jobu = 'Y';
	char jobv = 'Y';
	long jobu_len = 1;
	long jobv_len = 1;
	int ld_tmpu = m_rows;
	int ld_tmpv = n_cols;
	int ld_tmpb = kmax;

	// Compute singular values and vectors
	int info = 0;
	double tolin = 1.0E-4;
	if (performance_log)
	{
		performance_log->log_event("calling DEF_DLANSVD");
	}
	DEF_DLANSVD(&jobu, &jobv, &m_rows, &n_cols, &k, &kmax, tmp_u, &ld_tmpu, tmp_sigma, tmp_bnd,
		tmp_v, &ld_tmpv, &tolin, tmp_work, &lwork, tmp_iwork, &liwork, doption, ioption, &info,
		dparm, iparm, &jobu_len, &jobv_len);

	if (performance_log)
	{
		performance_log->log_event("updating Vt, S and U matrices");
	}
	//Compute number of singular values to be used in the solution
	int n_sing_used = 0;
	int actual_sing = kmax;
	if (info > 0)
		actual_sing = info;
	for (int i_sing = 0; i_sing < actual_sing; ++i_sing)
	{
		double eig_ratio = tmp_sigma[i_sing] / tmp_sigma[0];
		if (eig_ratio > _eigen_thres)
		{
			++n_sing_used;
		}
		else
		{
			break;
		}
	}
	// Update Sigma
	Sigma.resize(n_sing_used, n_sing_used);
	Sigma.setConstant(0.0);
	for (int i_sing = 0; i_sing<n_sing_used; ++i_sing)
	{
		Sigma(i_sing,i_sing) = tmp_sigma[i_sing];
	}

	/*Sigma_trunc.resize(kmax - n_sing_used);
	Sigma_trunc.setConstant(0.0);
	for (int i_sing = n_sing_used; i_sing<kmax; ++i_sing)
	{
		Sigma_trunc(i_sing - n_sing_used) = tmp_sigma[i_sing];
	}*/
	
	// Update U
	U.resize(m_rows, n_sing_used);
	U.setZero();
	for (int i_sing = 0; i_sing<n_sing_used; ++i_sing)
	{
		for (int irow = 0; irow<m_rows; ++irow)
		{
			U(irow, i_sing) = tmp_u[i_sing*m_rows + irow];
		}
	}
	
	// Update Vt
	V.resize(n_sing_used, n_cols);
	V.setZero();
	for (int i_sing = 0; i_sing<n_sing_used; ++i_sing)
	{
		for (int icol = 0; icol<n_cols; ++icol)
		{
			V(i_sing, icol) = tmp_v[i_sing*n_cols + icol];
		}
	}
	V.transposeInPlace();

	delete[] dparm;
	delete[] iparm;
	delete[] tmp_u;
	delete[] tmp_sigma;
	delete[] tmp_bnd;
	delete[] tmp_v;
	delete[] tmp_work;
	delete[] tmp_iwork;
}


SVD_PROPACK::~SVD_PROPACK(void)
{
}



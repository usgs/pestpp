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
#include "eigen_tools.h"
#include <Eigen/Dense>
#include "Transformable.h"
#include <vector>
#include <set>
#include <algorithm>
#include <iomanip>
#include <cstdio>
#include <cstdint>


using namespace Eigen;
using namespace std;


void get_MatrixXd_row_abs_max(const MatrixXd &m, int row, int *max_col, double *max_val)
{
	size_t nrows = m.rows();
	size_t ncols = m.cols();
	assert(row <= nrows);

	*max_col = -999;
	*max_val = 0;
	for (size_t icol=0; icol<ncols; ++icol)
	{
		if(abs(m(row,icol)) > abs(*max_val))
		{
			*max_col = icol;
			*max_val = m(row,icol);
		}
	}
}

VectorXd stlvec_2_eigenvec(const std::vector<double> &stl_vec)
{
	size_t len = stl_vec.size();
	VectorXd la_vec(len);
	for (size_t i=0; i<len; ++i)
	{
		la_vec(i) = stl_vec[i];
	}
	return la_vec;
}

vector<double> eigenvec_2_stlvec(const VectorXd &eigen_vec)
{
	size_t len = eigen_vec.size();
	vector<double> stl_vec;
	stl_vec.reserve(len);
	for (size_t i = 0; i<len; ++i)
	{
		stl_vec.push_back(eigen_vec(i));
	}
	return stl_vec;
}


Eigen::VectorXd transformable_2_eigen_vec(const Transformable &data, vector<string> oredered_names)
{
	size_t len = oredered_names.size();
	VectorXd new_vec(len);
	auto data_end = data.end();

	for (size_t i = 0; i<len; ++i)
	{
		const auto &it = data.find(oredered_names[i]);
		assert(it != data_end);
		new_vec(i) = it->second;
	}
	return new_vec;
}


void matrix_del_rows_cols(Eigen::SparseMatrix<double> &mat, const vector<size_t> &id_vec, bool perm_rows, bool perm_cols)
{
	if (!id_vec.empty())
	{
		size_t ncols = mat.cols();
		set<int> del_set(id_vec.begin(), id_vec.end());

		std::vector<Eigen::Triplet<int> > triplet_list;

		// add rows to be retained to the beginning of the  permutation matrix
		int icol_new = 0;
		int n_save = 0;
		for (int icol_old = 0; icol_old<ncols; ++icol_old)
		{
			if (del_set.find(icol_old) == del_set.end())
			{
				triplet_list.push_back(Eigen::Triplet<int>(icol_old, icol_new, 1));
				++icol_new;
				++n_save;
			}
		}
		// add rows to be deleted to end to permutation matrix
		for (int icol_old :  id_vec)
		{
			triplet_list.push_back(Eigen::Triplet<int>(icol_old, icol_new, 1));
			++icol_new;
		}
		Eigen::SparseMatrix<double> perm_sparse_matrix(ncols, ncols);
		perm_sparse_matrix.setZero();
		perm_sparse_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
		Eigen::SparseMatrix<double> new_matrix;
		if ((perm_rows) && (perm_cols))
		{
			new_matrix = perm_sparse_matrix.transpose() * mat * perm_sparse_matrix;
			mat = new_matrix.topLeftCorner(n_save, n_save);
		}
		else if (perm_cols)
		{
			new_matrix = mat * perm_sparse_matrix;
			mat = new_matrix.leftCols(n_save);
		}
		else if (perm_rows)
		{
			new_matrix = perm_sparse_matrix.transpose() * mat;
			mat = new_matrix.topRows(n_save);
		}

	}
}

Eigen::SparseMatrix<double> get_diag_matrix(const Eigen::SparseMatrix<double> &mat)
{
	auto  ncols = mat.cols();
	auto  nrows = mat.rows();

	assert(ncols == nrows);
	if (nrows != ncols)
	{
		PestError("Error in get_diag_matrix: Can not return the diagonal of a non-square matrix");
	}
	VectorXd diag_vec = mat.diagonal();
	std::vector<Eigen::Triplet<double> > triplet_list;
	triplet_list.reserve(nrows);
	for (size_t i = 0; i<nrows; ++i)
	{
		triplet_list.push_back(Eigen::Triplet<double>(i, i, diag_vec(i)));
	}

	Eigen::SparseMatrix<double> diag_matrix;
	diag_matrix.resize(nrows, ncols);
	diag_matrix.setZero();
	diag_matrix.setFromTriplets(triplet_list.begin(), triplet_list.end());
	return diag_matrix;
}

void print(const MatrixXd &mat, ostream & fout)
{
	size_t nrows = mat.rows();
	size_t ncols = mat.cols();

	for (size_t i=0; i<nrows; ++i)
	{
		for (size_t j=0; j<ncols; ++j)
		{
			fout << mat(i,j);
			if (j < ncols-1) {fout << ", ";}
		}
		fout << endl;
	}

}

void print(const MatrixXd &mat, ostream & fout, int n_per_line)
{
	size_t nrows = mat.rows();
	size_t ncols = mat.cols();

	for (size_t i=0; i<nrows; ++i)
	{
		for (size_t j=0; j<ncols; ++j)
		{
			fout << setw(15) << setiosflags(ios::right) << mat(i,j);
			if ((j+1) % (n_per_line) == 0 || j+1==ncols)
			{
				fout << '\n';
			}
		}
		fout << endl;
	}

}

void print(const VectorXd &vec, ostream & fout, int n_per_line)
{
	size_t n = vec.size();

	for (size_t i=0; i<n; ++i)
	{
		fout << showpoint << setw(15) << setiosflags(ios::right) << vec(i);
			if ((i+1) % (n_per_line) == 0 || i+1==n)
		{
			fout << '\n';
		}

	}
}

bool save_triplets_bin(const SparseMatrix<double> &mat, ostream &fout)
{
	int32_t xyn[3] = { static_cast<int32_t>(mat.rows()), static_cast<int32_t>(mat.cols()), static_cast<int32_t>(mat.nonZeros()) };
	fout.write((char*)xyn, sizeof(xyn));

	for (int k = 0; k < mat.outerSize(); ++k)
	{
		SparseMatrix<double>::InnerIterator it(mat, k);
		for (; it; ++it)
		{
			int32_t rc[2] = { static_cast<int32_t>(it.row()), static_cast<int32_t>(it.col() )};
			fout.write((char*)rc, sizeof(rc));
			double v = it.value();
			fout.write((char*)&v, sizeof(v));
		}
	}
	return true;
}

bool save_vector_bin(const VectorXd &vec, ostream &fout)
{
	int32_t size = vec.size();
	//vector<double> buf = egienvec_2_stlvec(vec);
	fout.write((char*)&size, sizeof(size));
	fout.write((char*)vec.data(), sizeof(double)*size);
	return true;
}

bool load_vector_bin(VectorXd &vec, istream &fin)
{
	int32_t size = 0;
	fin.read((char*)&size, sizeof(size));
	vec.resize(size);
	fin.read((char*)vec.data(), sizeof(double)*size);
	return true;
}

bool load_triplets_bin(SparseMatrix<double> &a, istream &fin)
{
	int32_t xyn[3];
	fin.read((char*)xyn, sizeof(xyn));
	a.resize(xyn[0], xyn[1]);
	vector<Triplet<double>> trips(xyn[2]);

	for (int k = 0; k < trips.size(); ++k)
	{
		int32_t rc[2];
		fin.read((char*)rc, sizeof(rc));
		double v;
		fin.read((char*)&v, sizeof(v));

		trips[k] = Triplet<double>(rc[0], rc[1], v);
	}
	a.setFromTriplets(trips.begin(), trips.end());
	return true;
}

Eigen::SparseMatrix<double> eigenvec_2_diagsparse(Eigen::VectorXd vec)
{
	vector<Eigen::Triplet<double>> triplet_list;
	for (int i = 0; i != vec.size(); i++)
		triplet_list.push_back(Eigen::Triplet<double>(i, i, vec[i]));
	Eigen::SparseMatrix<double> mat(vec.size(), vec.size());
	mat.setZero();
	mat.setFromTriplets(triplet_list.begin(), triplet_list.end());
	return mat;



}

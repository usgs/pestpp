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
#ifndef QSQRT_MATRIX_H_
#define QSQRT_MATRIX_H_

#include <vector>
#include <Eigen/Dense>
#include<Eigen/Sparse>

class ObservationInfo;
class Observations;
class PriorInformation;
class Parameters;
class DynamicRegularization;

using namespace std;

class QSqrtMatrix
{
public:
	QSqrtMatrix(){};
	QSqrtMatrix(const ObservationInfo *obs_info_ptr, const PriorInformation *prior_info_ptr);
	Eigen::SparseMatrix<double> get_sparse_matrix(const vector<string> &obs_names, const DynamicRegularization &regul, bool get_square = false) const;
	~QSqrtMatrix(void);
private:
	const ObservationInfo *obs_info_ptr;
	const PriorInformation *prior_info_ptr;
};

#endif /* QSQRT_MATRIX_H_ */

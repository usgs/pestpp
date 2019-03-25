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
#include <utility>
#include <cmath>
#include <Eigen/Dense>
#include <algorithm>
#include <vector>
#include "QSqrtMatrix.h"
#include "Transformable.h"
#include "PriorInformation.h"
#include "pest_data_structs.h"
#include "Regularization.h"

using namespace std;
using namespace Eigen;

QSqrtMatrix::QSqrtMatrix(const ObservationInfo *_obs_info_ptr, const PriorInformation *_prior_info_ptr)
: obs_info_ptr(_obs_info_ptr), prior_info_ptr(_prior_info_ptr)
{
}

Eigen::SparseMatrix<double> QSqrtMatrix::get_sparse_matrix(const vector<string> &obs_names, const DynamicRegularization &regul, bool get_sqaure) const
{

	Eigen::SparseMatrix<double> weights(obs_names.size(), obs_names.size());
	unordered_map<string, ObservationRec>::const_iterator found_obsinfo_iter;
	unordered_map<string, ObservationRec>::const_iterator non_found_obsinfo_iter = obs_info_ptr->observations.end();
	PriorInformation::const_iterator found_prior_info;
	PriorInformation::const_iterator not_found_prior_info = prior_info_ptr->end();
	vector<string>::const_iterator b;
	vector<string>::const_iterator e;
	int i = 0;
	double weight = 0;
	double tikhonov_weight = 1.0;
	const string *group = nullptr;
	// PEST convention is the the weights are 1/standard deviation but the regularization weight is the square
	// of this
	bool use_regul = regul.get_use_dynamic_reg();
	if (use_regul) tikhonov_weight = sqrt(regul.get_weight());
	std::vector<Eigen::Triplet<double> > triplet_list;
	for (i = 0, b = obs_names.begin(), e = obs_names.end(); b != e; ++b, ++i)
	{
		found_obsinfo_iter = obs_info_ptr->observations.find(*b);
		found_prior_info = prior_info_ptr->find(*b);
		// This section handles Observations
		if (found_obsinfo_iter != non_found_obsinfo_iter)
		{
			group = &((*found_obsinfo_iter).second.group);
			weight = (*found_obsinfo_iter).second.weight;
			bool is_reg_grp = ObservationGroupRec::is_regularization(*group);
			if (use_regul && is_reg_grp)
			{
				if (regul.get_adj_grp_weights() && is_reg_grp)
				{
					double grp_factor = regul.get_grp_weight_fact(*group);
					weight *= grp_factor;
				}
				weight *= tikhonov_weight;
			}
			if (get_sqaure) weight = weight * weight;
			triplet_list.push_back(Eigen::Triplet<double>(i, i, weight));
		}
		// This section handles Prior Information
		else if (found_prior_info != not_found_prior_info)
		{
			group = &((*found_prior_info).second.get_group());
			weight = (*found_prior_info).second.get_weight();
			bool is_reg_grp = (*found_prior_info).second.is_regularization();
			if (use_regul && is_reg_grp)
			{
				if (regul.get_adj_grp_weights())
				{
					double grp_factor = regul.get_grp_weight_fact(*group);
					weight *= grp_factor;
				}
				weight *= sqrt(regul.get_weight());
			}
			if (get_sqaure) weight = weight * weight;
			triplet_list.push_back(Eigen::Triplet<double>(i, i, weight));
		}
		else {
			assert(true);  //observation not in standard observations or prior information
		}
	}
	weights.resize(obs_names.size(), obs_names.size());
	weights.setZero();
	weights.setFromTriplets(triplet_list.begin(), triplet_list.end());
	return weights;
}

QSqrtMatrix::~QSqrtMatrix(void)
{
}

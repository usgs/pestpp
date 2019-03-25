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

#ifndef SERIALIZE_H_
#define SERIALIZE_H_

#include <vector>
#include <climits>
#include <string>

class Transformable;
class Parameters;
class Observations;

class Serialization
{
public:
	static std::vector<int8_t> serialize(int64_t data);
	static std::vector<int8_t> serialize(const Transformable &tr_data);
	static std::vector<int8_t> serialize(const std::vector<const Transformable*> tr_vec);
	static std::vector<int8_t> serialize(const std::vector<Transformable*> &tr_vec);
	static std::vector<int8_t> serialize(const Parameters &pars, const Observations &obs);
	static std::vector<int8_t> serialize(const Parameters &pars, const std::vector<std::string> &par_names_vec, const Observations &obs, const std::vector<std::string> &obs_names_vec, double run_time);
	static std::vector<int8_t> serialize(const std::vector<std::string> &string_vec);
	static std::vector<int8_t> serialize(const std::vector<std::vector<std::string> const*> &string_vec_vec);
	static unsigned long unserialize(const std::vector<int8_t> &ser_data, int64_t &data, unsigned long start_loc = 0);
	static unsigned long unserialize(const std::vector<int8_t> &ser_data, Transformable &tr_data, unsigned long start_loc = 0);
	static unsigned long unserialize(const std::vector<int8_t> &ser_data, std::vector<Transformable*> &tr_vec, unsigned long start_loc = 0);
	static unsigned long unserialize(const std::vector<int8_t> &ser_data, Parameters &pars, Observations &obs, unsigned long start_loc = 0);
	static unsigned long unserialize(const std::vector<int8_t> &ser_data, std::vector<std::string> &string_vec, unsigned long start_loc = 0, unsigned long max_read_bytes = ULONG_MAX);
	static unsigned long unserialize(const std::vector<int8_t> &ser_data, Transformable &items, const std::vector<std::string> &names_vec, unsigned long start_loc = 0);
	static unsigned long unserialize(const std::vector<int8_t> &ser_data, Parameters &pars, const std::vector<std::string> &par_names, Observations &obs, const std::vector<std::string> &obs_names, double &run_time);
private:
};


#endif /* SERIALIZE_H_ */





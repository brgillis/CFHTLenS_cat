/**********************************************************************\
 @file count_fitting_functor.cpp
 ------------------

 Source file for count_fitting_functor.

 **********************************************************************

 Copyright (C) 2014  Bryan R. Gillis

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

 \**********************************************************************/

#include <cassert>
#include <cstdlib>

#include "brg/global.h"

#include "brg/file_access/ascii_table_map.hpp"
#include "brg/math/functor/functor.hpp"
#include "brg/math/misc_math.hpp"
#include "brg/physics/units/unit_obj.h"

#include "count_fitting_functor.h"

void count_fitting_functor::_load() const
{
	if(_loaded_) return;
	auto map = brgastro::load_table_map<double>(_filename_);

	const auto count_label = "count";
	const auto bin_lower_label = "mag_bin_lower";

	size_t num_bins = map[count_label].size();
	_mag_bin_counts_.reserve(num_bins-2);
	_mag_bin_limits_.reserve(num_bins-1);

	// Fill in all but the overflow bins

	for(size_t i=1; i<num_bins-1; ++i)
	{
		_mag_bin_counts_.push_back(map[count_label][i]);
		_mag_bin_limits_.push_back(map[bin_lower_label][i]);
	}
	_mag_bin_limits_.push_back(map[bin_lower_label][num_bins-1]);

	_loaded_ = true;
}

std::vector<BRG_UNITS> count_fitting_functor::operator()( const std::vector<BRG_UNITS> & in_params,
		const bool silent ) const
{
	assert(_f_!=NULL);
	_load();
	_f_->set_params(in_params);
	double chi_sq = 0;
	for(size_t i=0; i<_mag_bin_counts_.size(); ++i)
	{
		double mid = (_mag_bin_limits_[i]+_mag_bin_limits_[i+1])/2;
		BRG_UNITS size = field_size()*(_mag_bin_limits_[i+1]-_mag_bin_limits_[i]);
		double error = _mag_bin_counts_[i]/std::sqrt(_mag_bin_counts_[i]-1);
		chi_sq += brgastro::square( ((*_f_)(mid,silent)*size-_mag_bin_counts_[i]) / error);
	}
	return std::vector<BRG_UNITS>(1,chi_sq);
}

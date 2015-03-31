/**********************************************************************\
 @file lensing_fitting_bin.hpp
 ------------------

 TODO <Insert file description here>

 **********************************************************************

 Copyright (C) 2015  Bryan R. Gillis

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

// body file: lensing_fitting_bin.cpp
#ifndef _BRG_LENSING_FITTING_BIN_HPP_INCLUDED_
#define _BRG_LENSING_FITTING_BIN_HPP_INCLUDED_

#include <cassert>
#include <vector>

#include "brg/container/labeled_array.hpp"
#include "brg/math/misc_math.hpp"

#include "common.h"
#include "fitting_bin.hpp"

/**
 *
 */
class lensing_fitting_bin
{
private:

	// Bin boundaries
	flt_type _z_min_, _z_max_;
	flt_type _m_min_, _m_max_;

	// R bins
	std::vector<fitting_bin> shear_fitting_bins;
	std::vector<fitting_bin> magf_fitting_bins;

public:

	// Constructors and destructors
#if(1)
	lensing_fitting_bin()
	{
	}
	lensing_fitting_bin(const brgastro::labeled_array<flt_type>::row_reference & row)
	{
		_z_min_ = row.at_label(z_min_label);
		_z_max_ = row.at_label(z_max_label);
		_m_min_ = row.at_label(m_min_label);
		_m_max_ = row.at_label(m_max_label);

		insert(row);
	}

	virtual ~lensing_fitting_bin()
	{
	}
#endif // Constructors and destructors

	// Set/insert data and bounds
#if(1)
	void set_bounds(const brgastro::labeled_array<flt_type>::row_reference & row)
	{
		_z_min_ = row.at_label(z_min_label);
		_z_max_ = row.at_label(z_max_label);
		_m_min_ = row.at_label(m_min_label);
		_m_max_ = row.at_label(m_max_label);

		return;
	}
	template< typename Tz1, typename Tz2, typename Tm1, typename Tm2 >
	void set_bounds( Tz1 && new_z_min, Tz2 && new_z_max, Tm1 && new_m_min, Tm2 && new_m_max )
	{
		_z_min_ = std::forward<Tz1>(new_z_min);
		_z_max_ = std::forward<Tz1>(new_z_max);
		_m_min_ = std::forward<Tz1>(new_m_min);
		_m_max_ = std::forward<Tz1>(new_m_max);

		return;
	}

	void insert(const brgastro::labeled_array<flt_type>::row_reference & row)
	{
		assert(_z_min_==row.at_label(z_min_label));
		assert(_z_max_==row.at_label(z_max_label));
		assert(_m_min_==row.at_label(m_min_label));
		assert(_m_max_==row.at_label(m_max_label));

		shear_fitting_bins.push_back(fitting_bin(row.at_label(shear_R_mean_label),
												 row.at_label(dS_t_mean_label),
												 row.at_label(dS_t_stderr_label)));
		magf_fitting_bins.push_back(fitting_bin(row.at_label(magf_R_mean_label),
												row.at_label(Sigma_mean_label),
												row.at_label(Sigma_stderr_label)));

		return;
	}
#endif // Set/insert data and bounds

	// Get bounds
#if(1)

	const flt_type & z_min() const { return _z_min_; }
	const flt_type & z_max() const { return _z_max_; }
	const flt_type & m_min() const { return _m_min_; }
	const flt_type & m_max() const { return _m_max_; }

	const flt_type z_mid() const { return (_z_min_+_z_max_)/2;}
	const flt_type m_mid() const { return (_m_min_+_m_max_)/2;}

#endif

	// Clear data
#if(1)
	void clear()
	{
		shear_fitting_bins.clear();
		magf_fitting_bins.clear();
	}
#endif // Clear data

	// Get Chi-squared for a given function
#if(1)
	template< typename func >
	flt_type get_shear_chi_squared(const func & f) const
	{
		flt_type chi_squared = 0;

		for( const auto & bin : shear_fitting_bins )
		{
			chi_squared += brgastro::square( (bin.y_mean-f(bin.x))/bin.y_stderr );
		}

		return chi_squared;
	}

	template< typename func >
	flt_type get_magf_chi_squared(const func & f) const
	{
		flt_type chi_squared = 0;

		for( const auto & bin : magf_fitting_bins )
		{
			chi_squared += brgastro::square( (bin.y_mean-f(bin.x))/bin.y_stderr );
		}

		return chi_squared;
	}
#endif // Get Chi-squared for a given function

};

#endif // _BRG_LENSING_FITTING_BIN_HPP_INCLUDED_

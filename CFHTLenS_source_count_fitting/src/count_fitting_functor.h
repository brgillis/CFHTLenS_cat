/**********************************************************************\
 @file count_fitting_functor.hpp
 ------------------

 TODO <Insert file description here>

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

// body file: count_fitting_functor.cpp

#ifndef _BRG_COUNT_FITTING_FUNCTOR_H_INCLUDED_
#define _BRG_COUNT_FITTING_FUNCTOR_H_INCLUDED_

#include "brg/global.h"

#include "brg/math/functor/functor.hpp"
#include "brg/physics/units/unit_conversions.hpp"
#include "brg/physics/units/unit_obj.h"

/**
 *
 */
class count_fitting_functor: public brgastro::functor<std::vector<BRG_UNITS>,BRG_UNITS> {
private:
	brgastro::functor<double,std::vector<BRG_UNITS>> *_f_;
	std::string _filename_;
	double _z_bin_size_;

	mutable std::vector<double> _mag_bin_limits_;
	mutable std::vector<double> _mag_bin_counts_;
	mutable bool _loaded_;

	void _load() const;

public:
	count_fitting_functor()
	: functor(1),
	  _f_(NULL),
	  _filename_(""),
	  _z_bin_size_(1),
	  _loaded_(false)

	{
	}
	count_fitting_functor(brgastro::functor<double,std::vector<BRG_UNITS>> *init_f,
			std::string init_filename="",
			BRG_UNITS init_field_size=1,
			double init_z_bin_size=1)
	: functor(init_field_size),
	  _f_(init_f),
	  _filename_(init_filename),
	  _z_bin_size_(init_z_bin_size),
	  _loaded_(false)
	{
	}
	virtual ~count_fitting_functor()
	{
	}

	// Setters
#if (1)

	void set_function(brgastro::functor<double,std::vector<BRG_UNITS>> *new_f)
	{
		_f_ = new_f;
	}
	void set_filename(std::string new_filename)
	{
		_filename_ = new_filename;
	}
	void set_field_size(CONST_BRG_UNITS_REF new_field_size)
	{
		set_params(new_field_size);
	}
	void set_z_bin_size_size(double new_z_bin_size)
	{
		_z_bin_size_ = new_z_bin_size;
	}

#endif

	// Accessors
#if (1)
	const brgastro::functor<double,std::vector<BRG_UNITS>> *function()
	{
		return _f_;
	}
	const std::string & filename()
	{
		return _filename_;
	}
	CONST_BRG_UNITS_REF field_size() const
	{
		return params();
	}
	double z_bin_size() const
	{
		return _z_bin_size_;
	}
#endif

	// Function method

	std::vector<BRG_UNITS> operator()( const std::vector<BRG_UNITS> & in_params,
			const bool silent = false ) const;
};

#endif // _BRG_COUNT_FITTING_FUNCTOR_H_INCLUDED_

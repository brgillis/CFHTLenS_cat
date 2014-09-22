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

	mutable std::vector<double> _mag_bin_limits_;
	mutable std::vector<double> _mag_bin_counts_;
	mutable bool _loaded_;

	void _load() const;

public:
	count_fitting_functor()
	: _f_(NULL),
	  _filename_(""),
	  _loaded_(false),
	  functor(1)
	{
	}
	count_fitting_functor(brgastro::functor<double,std::vector<BRG_UNITS>> *init_f,
			std::string init_filename="",
			BRG_UNITS init_field_size=1)
	: _f_(init_f),
	  _filename_(init_filename),
	  _loaded_(false),
	  functor(init_field_size)
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
#endif

	// Function method

	std::vector<BRG_UNITS> operator()( const std::vector<BRG_UNITS> & in_params,
			const bool silent = false ) const;
};

#endif // _BRG_COUNT_FITTING_FUNCTOR_H_INCLUDED_

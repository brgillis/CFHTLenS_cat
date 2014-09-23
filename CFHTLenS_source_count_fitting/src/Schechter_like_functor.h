/**********************************************************************\
 @file Schechter_like_functor.h
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

// body file: Schechter_like_functor.cpp

#ifndef _BRG_SCHECHTER_LIKE_FUNCTOR_H_INCLUDED_
#define _BRG_SCHECHTER_LIKE_FUNCTOR_H_INCLUDED_

#include <cassert>

#include "brg/math/functor/functor.hpp"

/**
 *
 */
class Schechter_like_functor: public brgastro::functor<double,std::vector<BRG_UNITS>> {
private:
	const size_t _num_params_ = 7;

public:
	Schechter_like_functor()
	{
	}
	Schechter_like_functor(const std::vector<BRG_UNITS> & init_params )
	: functor(init_params)
	{
	}
	virtual ~Schechter_like_functor()
	{
	}

	// Params accessors
#if (1)

	CONST_BRG_UNITS_REF N_scale() const
	{
		assert(params().size()==_num_params_);
		return params()[0];
	}
	CONST_BRG_UNITS_REF m_star() const
	{
		assert(params().size()==_num_params_);
		return params()[1];
	}
	CONST_BRG_UNITS_REF alpha() const
	{
		assert(params().size()==_num_params_);
		return params()[2];
	}
	CONST_BRG_UNITS_REF mag_lower_lim_sharpness() const
	{
		assert(params().size()==_num_params_);
		return params()[3];
	}
	CONST_BRG_UNITS_REF mag23_jump() const
	{
		assert(params().size()==_num_params_);
		return params()[4];
	}
	CONST_BRG_UNITS_REF mag_upper_lim() const
	{
		assert(params().size()==_num_params_);
		return params()[5];
	}
	CONST_BRG_UNITS_REF mag_upper_lim_sharpness() const
	{
		assert(params().size()==_num_params_);
		return params()[6];
	}

#endif

	// Function method

	double operator()( const double in_params,
			const bool silent = false ) const;

};

#endif // _BRG_SCHECHTER_LIKE_FUNCTOR_H_INCLUDED_

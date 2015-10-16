/**********************************************************************\
 @file exp_quadratic_functor.hpp
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

// body file: exp_quadratic_functor.cpp

#ifndef _BRG_EXP_QUADRATIC_FUNCTOR_H_INCLUDED_
#define _BRG_EXP_QUADRATIC_FUNCTOR_H_INCLUDED_

#include <cassert>
#include <utility>

#include "IceBRG_main/container/tuple.hpp"
#include "IceBRG_main/math/functor/functor.hpp"
#include "IceBRG_main/units/units.hpp"

namespace IceBRG {

typedef custom_unit_type<0,0,0,-2,0> inverse_square_angle;

/**
 *
 */
class exp_quadratic_functor: public functor<flt_type,
	inverse_square_angle,
	tuple<inverse_square_angle,flt_type,flt_type>> {
private:
	const size_t _num_params_ = 3;

public:
	typedef flt_type fin_type;
	typedef inverse_square_angle fout_type;
	typedef tuple<inverse_square_angle,flt_type,flt_type> params_type;

	typedef functor<flt_type, inverse_square_angle, params_type> base_type;

	exp_quadratic_functor()
	{
	}
	template< typename T >
	exp_quadratic_functor(T && init_params )
	: base_type(std::forward<T>(init_params))
	{
	}
	virtual ~exp_quadratic_functor()
	{
	}

	// Params accessors
#if (1)

	inverse_square_angle N_scale() const
	{
		return params().get<0>();
	}
	flt_type beta_0() const
	{
		return params().get<1>();
	}
	flt_type d_beta() const
	{
		return params().get<2>();
	}

#endif

	// Function method

	fout_type operator()( const fin_type & mag ) const
	{
		fin_type m = mag-23;
		fout_type result = N_scale() * std::exp(beta_0()*m + d_beta()*m*m);

		return result;
	}

};

} // namespace IceBRG

#endif // _BRG_EXP_QUADRATIC_FUNCTOR_H_INCLUDED_

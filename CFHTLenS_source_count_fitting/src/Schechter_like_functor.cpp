/**********************************************************************\
 @file Schechter_like_functor.cpp
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

#include <cmath>
#include <cstdlib>

#include "Schechter_like_functor.h"


double Schechter_like_functor::operator()( const double in_param,
			const bool silent = false ) const
{
	double x = std::pow(10,0.4*(m_star()-in_param));
	return N_scale()*std::pow(x,alpha())*std::exp(-x);
}

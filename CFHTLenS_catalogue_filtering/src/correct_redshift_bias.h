/**********************************************************************\
 @file correct_redshift_bias.h
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

// body file: correct_redshift_bias.cpp

#ifndef _BRG_SRC_CORRECT_REDSHIFT_BIAS_H_INCLUDED_
#define _BRG_SRC_CORRECT_REDSHIFT_BIAS_H_INCLUDED_

#include <string>
#include <vector>

#include "brg/math/interpolator/interpolator.h"

struct redshift_calibrator_wrapper
{
	static brgastro::interpolator redshift_calibration_interpolator;
};

void correct_redshift_bias( std::vector<std::string> & redshifts);


#endif // _BRG_SRC_CORRECT_REDSHIFT_BIAS_H_INCLUDED_

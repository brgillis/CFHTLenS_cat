/**********************************************************************\
 @file gg_lensing_config.h
 ------------------

 An object to store configuration settings for the gg_lensing program.

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

// body file: gg_lensing_config.cpp

#ifndef _BRG_GG_LENSING_CONFIG_H_INCLUDED_
#define _BRG_GG_LENSING_CONFIG_H_INCLUDED_

#include <cassert>
#include <cstdlib>

#include "brg/global.h"

#include "brg/physics/units/unit_conversions.hpp"
#include "brg/physics/units/unit_obj.h"

/**
 *
 */
struct gg_lensing_config {

	BRG_DISTANCE R_min, R_max, R_step;
	BRG_MASS m_min, m_max, m_step;
	double z_min, z_max, z_step;
	double mag_min, mag_max, mag_step;

	bool R_log, m_log, z_log, mag_log;
	size_t R_bins, m_bins, z_bins, mag_bins;

	double z_buffer;

	gg_lensing_config( const int argc=0, const char *argv[]=nullptr );
	virtual ~gg_lensing_config()
	{
	}
};

#endif // _BRG_GG_LENSING_CONFIG_H_INCLUDED_

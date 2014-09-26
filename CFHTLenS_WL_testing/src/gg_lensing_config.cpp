/**********************************************************************\
 @file gg_lensing_config.cpp
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

#include <limits>

#include "gg_lensing_config.h"

gg_lensing_config::gg_lensing_config( const int argc, const char *argv[] )
{
	using namespace brgastro::unitconv;

	if(argc==1)
	{
		R_min = 10*kpctom;
		R_max = 2000*kpctom;
		R_step = 10*kpctom;
		R_log = true;
		R_bins = 20;

		m_min = 1e9*Msuntokg;
		m_max = 1e10*Msuntokg;
		m_step = 1e11*Msuntokg;
		m_log = true;
		m_bins = 1;

		z_min = 0.2;
		z_max = 1.1;
		z_step = 0.9;
		z_log = true;
		z_bins = 1;

		mag_min = -std::numeric_limits<double>::infinity();
		mag_max = 25;
		mag_step = std::numeric_limits<double>::infinity();
		mag_log = false;
		mag_bins = 1;

		z_buffer = 0.1;

		return;
	}
	else
	{
		// TODO Construct from arguments passed to the program
		assert(false);
		return;
	}
}

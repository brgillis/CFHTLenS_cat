/**********************************************************************\
 @file pass_configs_to_binner.cpp
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


#include "brg/physics/lensing/pair_binner.h"

#include "gg_lensing_config.h"
#include "pass_configs_to_binner.h"


brgastro::pair_binner pass_configs_to_binner(gg_lensing_config config)
{
	brgastro::pair_binner binner(config.R_min,config.R_max,config.R_step,
			config.m_min, config.m_max, config.m_step,
			config.z_min, config.z_max, config.z_step,
			config.mag_min, config.mag_max, config.mag_step);

	if(config.R_log)
	{
		binner.set_log_R_limits(config.R_min,config.R_max,config.R_bins);
	}
	if(config.m_log)
	{
		binner.set_log_m_limits(config.m_min,config.m_max,config.m_bins);
	}
	if(config.z_log)
	{
		binner.set_log_z_limits(config.z_min,config.z_max,config.z_bins);
	}
	if(config.mag_log)
	{
		binner.set_log_mag_limits(config.mag_min,config.mag_max,config.mag_bins);
	}
	return binner;
}

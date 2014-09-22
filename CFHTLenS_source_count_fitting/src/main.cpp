/**********************************************************************\
 @file main.cpp
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

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/density.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/weighted_density.hpp>
#include <boost/lexical_cast.hpp>

#include "brg/file_access/ascii_table_map.hpp"
#include "brg/file_access/open_file.hpp"
#include "brg/file_access/table_typedefs.hpp"
#include "brg/math/misc_math.hpp"
#include "brg/physics/astro.h"
#include "brg/physics/lensing/pair_binner.h"
#include "brg/physics/lensing/pair_bins_summary.h"
#include "brg/physics/lensing/source_galaxy.h"
#include "brg/physics/sky_obj/galaxy.h"
#include "brg/vector/limit_vector.hpp"
#include "brg/vector/manipulations.hpp"

// Magic values
std::string fields_directory = "/disk2/brg/git/CFHTLenS_cat/Data/";
std::string output_table_root = fields_directory + "magnitude_hist_z";
std::string output_table_tail = ".dat";
unsigned int zlo = 0;
unsigned int zstep = 10;
unsigned int zhi = 390;

int main( const int argc, const char *argv[] )
{

	// Set up the redshift bins
	std::vector<double> z_bin_limits = brgastro::make_limit_vector<double>(zlo*0.01,zhi*0.01,zstep*0.01);


	return 0;
}

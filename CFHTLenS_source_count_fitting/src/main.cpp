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
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/lexical_cast.hpp>

#include "brg/file_access/ascii_table_map.hpp"
#include "brg/math/misc_math.hpp"
#include "brg/math/solvers/solvers.hpp"
#include "brg/vector/elementwise_functions.hpp"

#include "count_fitting_functor.h"
#include "Schechter_like_functor.h"

// Magic values
const std::string fields_directory = "/disk2/brg/git/CFHTLenS_cat/Data/";
const std::string count_table_root = fields_directory + "magnitude_hist_z";
const std::string count_table_tail = ".dat";
const std::string output_filename = fields_directory + "count_fitting_results.dat";
const unsigned int zlo = 20;
const unsigned int zstep = 10;
const unsigned int zhi = 390;
const double z_bin_size = zstep*0.01;

const BRG_UNITS field_size(130.98*brgastro::square(brgastro::unitconv::degtorad));

int main( const int argc, const char *argv[] )
{
	// General set-up
	Schechter_like_functor estimator;

	std::vector<BRG_UNITS> mins, maxes, steps, inits;

	// Bounds for N_scale
	mins.push_back(0);
	maxes.push_back(1e5/brgastro::square(brgastro::unitconv::degtorad));
	steps.push_back(1e3/brgastro::square(brgastro::unitconv::degtorad));
	inits.push_back(3e3/brgastro::square(brgastro::unitconv::degtorad));

	// Bounds for m_star
	mins.push_back(15);
	maxes.push_back(30);
	steps.push_back(0.1);
	inits.push_back(18.5);

	// Bounds for alpha
	mins.push_back(-10);
	maxes.push_back(2);
	steps.push_back(0.05);
	inits.push_back(-1);

	// Bounds for mag_lower_lim_sharpness
	mins.push_back(0);
	maxes.push_back(2);
	steps.push_back(0.1);
	inits.push_back(1);

	// Bounds for mag23_jump
	mins.push_back(0);
	maxes.push_back(1e5/brgastro::square(brgastro::unitconv::degtorad));
	steps.push_back(1e3/brgastro::square(brgastro::unitconv::degtorad));
	inits.push_back(2e3/brgastro::square(brgastro::unitconv::degtorad));

	// Bounds for mag_upper_lim
	mins.push_back(24);
	maxes.push_back(30);
	steps.push_back(0.1);
	inits.push_back(25);

	// Bounds for mag_upper_lim_sharpness
	mins.push_back(0);
	maxes.push_back(10);
	steps.push_back(0.1);
	inits.push_back(1);

	// Result map
	brgastro::table_map_t<double> result_map;

	for(unsigned int z100=zlo; z100<=zhi; z100+=zstep)
	{
		try
		{
			std::string filename = count_table_root + boost::lexical_cast<std::string>(z100)
					+ count_table_tail;
			count_fitting_functor fitter(&estimator,filename,field_size,z_bin_size);

			std::vector<BRG_UNITS> result_in_params = brgastro::solve_MCMC(&fitter,
					inits, mins, maxes, steps, 1000000, 25000);

			std::cout << "Best params for z=" << z100/100. << ":\t"
					<< result_in_params.at(0)*brgastro::square(brgastro::unitconv::degtorad) << " "
					<< result_in_params.at(1) << " "
					<< result_in_params.at(2) << " "
					<< result_in_params.at(3) << " "
					<< result_in_params.at(4)*brgastro::square(brgastro::unitconv::degtorad) << " "
					<< result_in_params.at(5) << " "
					<< result_in_params.at(6) << "\t";

			std::cout << "Chi^2: " << fitter(result_in_params).at(0) << std::endl;

			inits = result_in_params; // For next bin

			// Add this to the result map
			result_map["z_mid"].push_back(0.01*(z100+zstep/2));
			result_map["N_scale"].push_back(result_in_params.at(0));
			result_map["m_star_lower"].push_back(result_in_params.at(1));
			result_map["alpha"].push_back(result_in_params.at(2));
			result_map["lower_cutoff_sharpness"].push_back(result_in_params.at(3));
			result_map["mag23_jump"].push_back(result_in_params.at(4));
			result_map["m_star_upper"].push_back(result_in_params.at(5));
			result_map["upper_cutoff_sharpness"].push_back(result_in_params.at(6));
		}
		catch(const std::exception &e)
		{
			std::cerr << e.what();
		}
	}

	brgastro::print_table_map(output_filename,result_map);
	brgastro::print_table_map(std::cout,result_map);

	return 0;
}

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
#include "brg/units/units.hpp"
#include "brg/vector/elementwise_functions.hpp"

#include "count_fitting_functor.hpp"
#include "Schechter_like_functor.h"

// Magic values
const std::string fields_directory = "/disk2/brg/git/CFHTLenS_cat/Data/";
const std::string count_table_root = fields_directory + "magnitude_hist_z";
const std::string count_table_tail = ".dat";
const std::string output_filename = fields_directory + "count_fitting_results.dat";
const unsigned int zlo = 20;
const unsigned int zstep = 4;
const unsigned int zhi = 196;
const double z_bin_size = zstep*0.01;

using namespace brgastro;

const custom_unit_type<0,0,0,2,0> field_size(130.98*square((M_PI/180 * radian)));

int main( const int argc, const char *argv[] )
{
	// General set-up
	typedef Schechter_like_functor estimator_type;

	Schechter_like_functor estimator;

	typedef estimator_type::params_type params_type;

	const params_type mins( 0./square(unitconv::degtorad*rad), // N_scale
				      	    15., // m_star
							-10., // alpha
							0., // mag_lower_lim_sharpness
							0./square(unitconv::degtorad*rad), // mag23_jump
							24., // mag_upper_lim
							0.); // mag_upper_lim_sharpness

	const params_type maxes( 1.e5/square(unitconv::degtorad*rad), // N_scale
      	    				 30., // m_star
							 2., // alpha
							 2., // mag_lower_lim_sharpness
							 0./square(unitconv::degtorad*rad), // mag23_jump
							 30., // mag_upper_lim
							 10.); // mag_upper_lim_sharpness

	const params_type steps( 1.e3/square(unitconv::degtorad*rad), // N_scale
      	    				 0.1, // m_star
							 0.05, // alpha
							 0.1, // mag_lower_lim_sharpness
							 0./square(unitconv::degtorad*rad), // mag23_jump
							 0.1, // mag_upper_lim
							 0.1); // mag_upper_lim_sharpness

	params_type inits( 3.e3/square(unitconv::degtorad*rad), // N_scale
      	    				 18.5, // m_star
							 -1., // alpha
							 1., // mag_lower_lim_sharpness
							 0./square(unitconv::degtorad*rad), // mag23_jump
							 25., // mag_upper_lim
							 1.); // mag_upper_lim_sharpness

	// Result map
	table_map_t<double> result_map;

	for(unsigned int z100=zlo; z100<=zhi; z100+=zstep)
	{
		std::string filename = count_table_root + boost::lexical_cast<std::string>(z100)
				+ count_table_tail;
		count_fitting_functor<decltype(estimator)> fitter(&estimator,filename,field_size,z_bin_size);

		try
		{
			params_type result_in_params = solve_MCMC(fitter,
					inits, mins, maxes, steps, 1000000, 25000);

			std::cout << "Best params for z=" << z100/100. << ":\t"
					<< value_of(result_in_params.get<0>()*square(unitconv::degtorad)) << " "
					<< result_in_params.get<1>() << " "
					<< result_in_params.get<2>() << " "
					<< result_in_params.get<3>() << " "
					<< value_of(result_in_params.get<4>()*square(unitconv::degtorad)) << " "
					<< result_in_params.get<5>() << " "
					<< result_in_params.get<6>() << "\t";

			std::cout << "Chi^2: " << fitter(result_in_params) << std::endl;

			inits = result_in_params; // For next bin

			// Add this to the result map
			result_map["z_mid"].push_back(0.01*(z100+zstep/2));
			result_map["N_scale"].push_back(value_of(result_in_params.get<0>()));
			result_map["m_star_lower"].push_back(result_in_params.get<1>());
			result_map["alpha"].push_back(result_in_params.get<2>());
			result_map["lower_cutoff_sharpness"].push_back(result_in_params.get<3>());
			result_map["mag23_jump"].push_back(value_of(result_in_params.get<4>()));
			result_map["m_star_upper"].push_back(result_in_params.get<5>());
			result_map["upper_cutoff_sharpness"].push_back(result_in_params.get<6>());
		}
		catch(const std::exception &e)
		{
			// Most likely no counts in this bin, so use the inits (to match last good bin), but
			// with N_scale and mag23_jump set to zero

			result_map["z_mid"].push_back(0.01*(z100+zstep/2));
			result_map["N_scale"].push_back(0.);
			result_map["m_star_lower"].push_back(inits.get<1>());
			result_map["alpha"].push_back(inits.get<2>());
			result_map["lower_cutoff_sharpness"].push_back(inits.get<3>());
			result_map["mag23_jump"].push_back(0.);
			result_map["m_star_upper"].push_back(inits.get<5>());
			result_map["upper_cutoff_sharpness"].push_back(inits.get<6>());
		}
	}

	print_table_map(output_filename,result_map);
	print_table_map(std::cout,result_map);

	return 0;
}

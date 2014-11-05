/**********************************************************************\
 @file main_make_mock_cats.cpp
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

#include <cassert>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <boost/serialization/vector.hpp>

#include "brg/file_access/binary_archive.hpp"
#include "brg/file_access/open_file.hpp"
#include "brg/file_access/ascii_table_map.hpp"
#include "brg/math/misc_math.hpp"
#include "brg/math/random/random_functions.h"
#include "brg/physics/units/unit_conversions.hpp"

#include "get_ra_dec.h"
#include "num_good_pixels.hpp"

// Magic values
const std::string data_directory = "/disk2/brg/git/CFHTLenS_cat/Data/";
const std::string mask_directory = data_directory + "masks/";
const std::string field_directory = data_directory + "filtered_tables/";
const std::string fields_list = data_directory + "fields_list.txt";
const std::string pixel_map_root = "_lens_good_pixels.bin";
const std::string lens_output_root = "_mock_lens.dat";
const std::string source_output_root = "_mock_source.dat";

constexpr unsigned num_lenses_to_generate = 300000;
constexpr unsigned num_sources_to_generate = 500000;
constexpr double min_lens_z = 0.2;
constexpr double max_lens_z = 0.8;
constexpr double min_source_z = 0.2;
constexpr double max_source_z = 2.0;
constexpr double min_good_frac = 0.01;

int main( const int argc, const char *argv[] )
{
	// General setup
#ifdef _OPENMP
	omp_set_num_threads(5);
#endif
	
	// Open and read in the fields list
	std::ifstream fi;
	brgastro::open_file_input(fi,fields_list);

	std::vector<std::string> field_names;
	std::string temp_field_name;

	while(fi>>temp_field_name)
	{
		field_names.push_back(temp_field_name);
	}

	size_t num_fields = field_names.size();

	// Loop over all fields
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for(unsigned i=0;i<num_fields;++i)
	{
		std::string & field_name = field_names[i];
		std::string field_name_root = field_name.substr(0,6);

		std::string pixel_map_name = field_directory + field_name_root + pixel_map_root;
		std::string lens_output_file_name = field_directory + field_name_root + lens_output_root;
		std::string source_output_file_name = field_directory + field_name_root + source_output_root;

		const std::vector<std::vector<bool>> good_pixels =
				brgastro::binary_load_vector<std::vector<std::vector<bool>>>(pixel_map_name);

		const size_t ncol = good_pixels.size();
		assert(ncol>0);
		const size_t nrow = good_pixels[0].size();

		assert(static_cast<double>(num_good_pixels(good_pixels))/(ncol*nrow)>min_good_frac);

		std::vector<std::pair<double,double>> galaxy_positions;

		// For each galaxy we want to generate, get good random pixel positions for it
		for(unsigned j=0; j<num_lenses_to_generate; ++j)
		{
			bool bad_result=true;

			// Generate a random position inside the mask
			while(bad_result)
			{
				const double x = brgastro::drand(0,ncol)-0.5;
				const double y = brgastro::drand(0,ncol)-0.5;

				const unsigned xp = brgastro::round_int(x);
				const unsigned yp = brgastro::round_int(y);

				if(good_pixels[xp][yp])
				{
					bad_result = false;
					galaxy_positions.push_back(std::make_pair(x,y));
				}
			}
		}

		auto galaxy_sky_positions =
				get_ra_dec(field_name,galaxy_positions);

		// Create the output table, and then add each galaxy to it
		brgastro::table_map_t<double> field_output_table;

		for( size_t j=0; j<num_lenses_to_generate; ++j )
		{
			field_output_table["SeqNr"].push_back(j);
			field_output_table["ra_radians"].push_back(galaxy_sky_positions[j].first.first*brgastro::unitconv::degtorad);
			field_output_table["dec_radians"].push_back(galaxy_sky_positions[j].first.second*brgastro::unitconv::degtorad);
			field_output_table["Xpos"].push_back(galaxy_sky_positions[j].second.first);
			field_output_table["Ypos"].push_back(galaxy_sky_positions[j].second.second);
			field_output_table["Z_B"].push_back(brgastro::drand(min_lens_z,max_lens_z));
			field_output_table["T_B"].push_back(1);
			field_output_table["ODDS"].push_back(1);
			field_output_table["CHI_SQUARED_BPZ"].push_back(1);
			field_output_table["Mstel_kg"].push_back(1e40);
			field_output_table["Mstel_lo_kg"].push_back(0.5e40);
			field_output_table["Mstel_hi_kg"].push_back(2e40);
			field_output_table["MAG_i"].push_back(23);
			field_output_table["MAGERR_i"].push_back(1);
			field_output_table["EXTINCTION_i"].push_back(0);
			field_output_table["MAG_r"].push_back(23);
			field_output_table["MAGERR_r"].push_back(1);
			field_output_table["EXTINCTION_r"].push_back(0);
		}

		// Output the table
		#ifdef _OPENMP
		#pragma omp critical(output_mock_lenses)
		#endif
		{
			brgastro::print_table_map(lens_output_file_name,field_output_table);
			std::cout << "Finished generating " << lens_output_file_name << ".\n";
		}

		// Repeat for sources

		galaxy_positions.clear();

		// For each galaxy we want to generate, get good random pixel positions for it
		for(unsigned j=0; j<num_sources_to_generate; ++j)
		{
			bool bad_result=true;

			// Generate a random position inside the mask
			while(bad_result)
			{
				const double x = brgastro::drand(0,ncol)-0.5;
				const double y = brgastro::drand(0,ncol)-0.5;

				const unsigned xp = brgastro::round_int(x);
				const unsigned yp = brgastro::round_int(y);

				if(good_pixels[xp][yp])
				{
					bad_result = false;
					galaxy_positions.push_back(std::make_pair(x,y));
				}
			}
		}

		galaxy_sky_positions =
				get_ra_dec(field_name,galaxy_positions);

		// Create the output table, and then add each galaxy to it
		field_output_table.clear();

		for( size_t j=0; j<num_sources_to_generate; ++j )
		{
			field_output_table["SeqNr"].push_back(j);
			field_output_table["ra_radians"].push_back(galaxy_sky_positions[j].first.first*brgastro::unitconv::degtorad);
			field_output_table["dec_radians"].push_back(galaxy_sky_positions[j].first.second*brgastro::unitconv::degtorad);
			field_output_table["Xpos"].push_back(galaxy_sky_positions[j].second.first);
			field_output_table["Ypos"].push_back(galaxy_sky_positions[j].second.second);
			field_output_table["Z_B"].push_back(brgastro::drand(min_lens_z,max_lens_z));
			field_output_table["T_B"].push_back(1);
			field_output_table["ODDS"].push_back(1);
			field_output_table["e1"].push_back(0);
			field_output_table["e2"].push_back(0);
			field_output_table["weight"].push_back(1);
			field_output_table["m"].push_back(0);
			field_output_table["c2"].push_back(0);
			field_output_table["Mstel_kg"].push_back(1e40);
			field_output_table["Mstel_lo_kg"].push_back(0.5e40);
			field_output_table["Mstel_hi_kg"].push_back(2e40);
			field_output_table["MAG_i"].push_back(23);
			field_output_table["MAGERR_i"].push_back(1);
			field_output_table["EXTINCTION_i"].push_back(0);
			field_output_table["MAG_r"].push_back(23);
			field_output_table["MAGERR_r"].push_back(1);
			field_output_table["EXTINCTION_r"].push_back(0);
		}

		// Output the table
		#ifdef _OPENMP
		#pragma omp critical(output_mock_sources)
		#endif
		{
			brgastro::print_table_map(source_output_file_name,field_output_table);
			std::cout << "Finished generating " << source_output_file_name << ".\n";
		}
	}

	return 0;
}

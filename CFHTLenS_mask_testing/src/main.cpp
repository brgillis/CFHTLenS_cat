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

#include <cmath>
#include <iostream>
#include <fstream>
#include <functional>
#include <sstream>
#include <string>
#include <unordered_map>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <CCfits/CCfits>

#include "brg/file_access/ascii_table_map.hpp"
#include "brg/file_access/open_file.hpp"
#include "brg/math/misc_math.hpp"
#include "brg/vector/limit_vector.hpp"

#include "get_good_positions.hpp"
#include "is_good_position.hpp"
#include "load_pixel_table.h"

// Magic values
std::string data_directory = "/disk2/brg/git/CFHTLenS_cat/Data/";
std::string mask_directory = data_directory + "masks/";
std::string fields_list = data_directory + "fields_list.txt";

/*
 * Pixel is valid galaxy position if surrounded by unmasked pixels
 *
 * Max physical separation is 2000 kpc.
 * At redshift 0.2 this is 601.5"
 * Pixel scale is 0.186"/px
 * Maximum pixel separation we care about is thus 3235 px.
 */

constexpr double min_kpc_sep=0;
constexpr double max_kpc_sep=2000;
constexpr double kpc_sep_step=10;

constexpr unsigned min_px_sep=0;
constexpr unsigned max_px_sep=3240;
constexpr unsigned px_sep_step=10;



constexpr unsigned short sampling_factor=2;

static_assert(sampling_factor>0,"Subsample must be positive.");

int main( const int argc, const char *argv[] )
{

	// General setup
#ifdef _OPENMP
	omp_set_num_threads(5);
#endif

	// Set up separation limits vector
	std::vector<unsigned> sep_limits = brgastro::make_limit_vector<unsigned>(min_px_sep,max_px_sep,px_sep_step);

	const unsigned num_bins = sep_limits.size()-1;

	// Open and read in the fields list
	std::ifstream fi;
	brgastro::open_file_input(fi,fields_list);

	std::vector<std::string> field_names;
	std::string temp_field_name;

	while(fi>>temp_field_name)
	{
		field_names.push_back(temp_field_name);
	}

	// Set up map for each field's separation hists
	brgastro::table_map_t<float> result;

	// Add a column to the result table which gives pixel limits
	std::vector<float> bin_mids(sep_limits.size()-1);
	for(unsigned i=0; i<num_bins; ++i)
	{
		bin_mids[i] = (sep_limits[i]+sep_limits[i+1])/2;
	}
	result["bin_px_mid"] = bin_mids;

#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for(unsigned i=0;i<field_names.size();++i)
	{
		std::string & field_name = field_names[i];
		std::string field_name_root = field_name.substr(0,6);

		// Get the input file name
		std::stringstream ss("");
		ss << mask_directory << field_name << ".fits.fz";
		std::string input_file_name = ss.str();

		std::vector<std::vector<bool>> good_positions;

		try
		{
			good_positions = get_good_positions(load_pixel_table(input_file_name));
		}
		catch(const CCfits::FitsError &e)
		{
			std::cerr << "Cannot open file " << input_file_name << "!\n"
					<< e.message() << std::endl;
			continue;
		}
		catch(const std::runtime_error &e)
		{
			std::cerr << "Cannot process file " << input_file_name << "!\n"
					<< e.what() << std::endl;
			continue;
		}

#ifdef _OPENMP
		#pragma omp critical(mask_status_update)
#endif
		{
			std::cout << "Mask data loaded for " << field_name_root << ".\n";
		}

		// We'll use that table for possible source positions (since we care about where
		// they could be, and not just where they are). For lenses, we can save time
		// by using only the actual positions of galaxies.

		// Get the lens table file name
		ss.str("");
		ss << data_directory << "filtered_tables/" << field_name_root << "_lens.dat";
		std::string lens_file_name = ss.str();

		// Load the lens table
		const auto lens_table_map = brgastro::load_table_map<double>(lens_file_name);

		// Set up vectors to store result data
		std::vector<unsigned> total_px_per_bin(num_bins,0);
		std::vector<unsigned> good_px_per_bin(num_bins,0);

		// Loop through lens positions now to calculate data
		unsigned num_lenses = lens_table_map.at("Xpos").size();
		for(size_t lens_i = 0; lens_i<num_lenses; ++lens_i)
		{
			const float lens_x = lens_table_map.at("Xpos")[lens_i];
			const float lens_y = lens_table_map.at("Ypos")[lens_i];

			const int lens_xp = static_cast<int>(lens_x);
			const int lens_yp = static_cast<int>(lens_y);

			const int lens_xp_min = lens_xp-max_px_sep-1;
			const int lens_yp_min = lens_yp-max_px_sep-1;

			const int lens_xp_max = lens_xp+max_px_sep+2;
			const int lens_yp_max = lens_yp+max_px_sep+2;

			unsigned short x_start_offset=0;
			unsigned short y_start_offset=0;

			if(sampling_factor>1)
			{
				x_start_offset = std::floor(sampling_factor*rand());
				y_start_offset = std::floor(sampling_factor*rand());
			}

			// Loop over possible source positions
			for(int i=lens_xp_min+x_start_offset; i<lens_xp_max; i+=sampling_factor)
			{
				// Loop over possible source positions
				for(int j=lens_yp_min+y_start_offset; j<lens_yp_max; j+=sampling_factor)
				{
					const unsigned d = brgastro::dist2d(lens_x-i,lens_y-j);
					if(d>max_px_sep) continue;
					unsigned bin_i = brgastro::get_bin_index(d,sep_limits);

					++total_px_per_bin[bin_i];
					if(is_good_position(i,j,good_positions))
					{
						++good_px_per_bin[bin_i];
					}
				}
			}
		}

		// Set up result vector for this field
		std::vector<float> field_res(num_bins);
		for(unsigned i=0; i<num_bins; ++i)
		{
			field_res[i] = static_cast<float>(good_px_per_bin[i])/static_cast<float>(total_px_per_bin[i]);
		}

		// Add this to the full result map

#ifdef _OPENMP
		#pragma omp critical(add_mask_test_results_to_map)
#endif
		{
			result[field_name_root] = field_res;
			std::cout << "Finished processing field " << field_name_root << "!\n";
		}
	}

	// Get the output file name
	std::stringstream ss("");
	ss << mask_directory << "position_corr.dat";
	std::string output_name = ss.str();

	brgastro::print_table_map(output_name,result);

	std::cout << "Done!\n";

	return 0;
}


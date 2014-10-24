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

#include <boost/lexical_cast.hpp>
#include <boost/serialization/vector.hpp>

#include <CCfits/CCfits>

#include "brg/file_access/ascii_table_map.hpp"
#include "brg/file_access/binary_archive.hpp"
#include "brg/file_access/open_file.hpp"
#include "brg/math/misc_math.hpp"
#include "brg/math/random/random_functions.h"
#include "brg/physics/astro.h"
#include "brg/physics/units/unit_conversions.hpp"
#include "brg/vector/limit_vector.hpp"

#include "get_good_positions.hpp"
#include "is_good_position.hpp"
#include "load_pixel_table.h"

// Magic values
const std::string data_directory = "/disk2/brg/git/CFHTLenS_cat/Data/";
const std::string mask_directory = data_directory + "masks/";
const std::string field_directory = data_directory + "filtered_tables/";
const std::string fields_list = data_directory + "fields_list.txt";
const std::string lens_output_root = "_lens_mask_frac.dat";
const std::string lens_pixel_map_root = "_lens_good_pixels.bin";

/*
 * Max physical separation is 2000 kpc.
 * At redshift 0.2 this is 601.5"
 * Pixel scale is 0.186"/px
 * Maximum pixel separation we care about is thus 3235 px.
 */

constexpr float min_kpc_sep=0;
constexpr float max_kpc_sep=2000;
constexpr float kpc_sep_step=10;

constexpr double rad_per_px = 0.186*brgastro::unitconv::asectorad;

constexpr unsigned max_px_sep=3240;

constexpr unsigned short sampling_factor=4;

static_assert(sampling_factor>0,"Subsample must be positive.");

int main( const int argc, const char *argv[] )
{

	// General setup
#ifdef _OPENMP
	omp_set_num_threads(5);
#endif

	// Set up separation limits vector
	std::vector<float> sep_limits =
			brgastro::make_limit_vector<float>(min_kpc_sep,max_kpc_sep,kpc_sep_step);

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

	// Set up map for each field's average separation hists
	brgastro::table_map_t<float> result;

	// Add a column to the result table which gives pixel limits
	std::vector<float> bin_mids(sep_limits.size()-1);
	for(unsigned i=0; i<num_bins; ++i)
	{
		bin_mids[i] = (sep_limits[i]+sep_limits[i+1])/2;
	}
	result["bin_mid_kpc"] = bin_mids;

	// Set up the field sizes table
	brgastro::table_map_t<double> field_sizes;

	size_t num_fields = field_names.size();

#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for(unsigned i=0;i<num_fields;++i)
	{
		std::string & field_name = field_names[i];
		std::string field_name_root = field_name.substr(0,6);

		// Get the input file name
		std::stringstream ss("");
		ss << mask_directory << field_name << ".fits.fz";
		std::string input_file_name = ss.str();

		std::vector<std::vector<bool>> good_pixels;
		std::vector<std::vector<bool>> good_positions;

		try
		{
			good_pixels = load_pixel_table(input_file_name);
			good_positions = get_good_positions(good_pixels);
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
		ss << field_directory << field_name_root << "_lens.dat";
		std::string lens_file_name = ss.str();

		// Get the lens output file name
		ss.str("");
		ss << field_directory << field_name_root << lens_output_root;
		std::string lens_output_file_name = ss.str();

		// Get the lens pixel map file name
		ss.str("");
		ss << field_directory << field_name_root << lens_pixel_map_root;
		std::string lens_pixel_map_file_name = ss.str();

		// Load the lens table
		const auto lens_table_map = brgastro::load_table_map<double>(lens_file_name);

		// Set up vectors to store result data
		std::vector<unsigned> total_px_per_bin(num_bins,0);
		std::vector<unsigned> good_px_per_bin(num_bins,0);

		// Set up a map for each lens's info
		brgastro::table_map_t<float> field_result;

		// Add column to this table for bin middles
		field_result["bin_mid_kpc"] = bin_mids;

		// Loop through lens positions now to calculate data
		unsigned num_lenses = lens_table_map.at("Xpos").size();
		for(size_t lens_i = 0; lens_i<num_lenses; ++lens_i)
		{
			// Get the lens's position on the image
			const unsigned lens_xp = brgastro::round_int(lens_table_map.at("Xpos")[lens_i]);
			const unsigned lens_yp = brgastro::round_int(lens_table_map.at("Ypos")[lens_i]);

			// Get the lens's redshift and calculate related information
			const float lens_z = lens_table_map.at("Z_B")[lens_i];

			const double pxfd_fact = brgastro::afd(brgastro::unitconv::kpctom,lens_z)/rad_per_px;

			const int lens_max_px_sep = pxfd_fact*max_kpc_sep;

			// Initialize vectors for the lens's masked fraction
			std::vector<unsigned> lens_good_px_per_bin(num_bins,0);
			std::vector<unsigned> lens_total_px_per_bin(num_bins,0);

			auto increment_bin = [&] (unsigned bin_i, bool good)
			{
				++total_px_per_bin[bin_i];
				++lens_total_px_per_bin[bin_i];
				if(good)
				{
					++good_px_per_bin[bin_i];
					++lens_good_px_per_bin[bin_i];
				}
			};

			// Loop over possible source positions

			// Start with those along horizontal and vertical directions
			for(int i=0; i<lens_max_px_sep; i+=sampling_factor)
			{
				const float d = i;
				unsigned bin_i = brgastro::get_bin_index<float>(d/pxfd_fact,sep_limits);

				increment_bin(bin_i,is_good_position(lens_xp+i,lens_yp,good_pixels));
				increment_bin(bin_i,is_good_position(lens_xp-i,lens_yp,good_pixels));
				increment_bin(bin_i,is_good_position(lens_xp,lens_yp+i,good_pixels));
				increment_bin(bin_i,is_good_position(lens_xp,lens_yp-i,good_pixels));
			}

			// Now do diagonals
			for(int i=0; i<lens_max_px_sep/std::sqrt(2.); i+=sampling_factor)
			{
				const float d = i*std::sqrt(2.);
				unsigned bin_i = brgastro::get_bin_index<float>(d/pxfd_fact,sep_limits);

				increment_bin(bin_i,is_good_position(lens_xp+i,lens_yp+i,good_pixels));
				increment_bin(bin_i,is_good_position(lens_xp-i,lens_yp+i,good_pixels));
				increment_bin(bin_i,is_good_position(lens_xp+i,lens_yp-i,good_pixels));
				increment_bin(bin_i,is_good_position(lens_xp-i,lens_yp-i,good_pixels));
			}

			// And now all other positions
			for(int i=2; i<lens_max_px_sep; i+=sampling_factor)
			{
				// Loop over possible source positions
				for(int j=1; j<i; j+=sampling_factor)
				{
					const float d = brgastro::dist2d(i,j);
					if(d>lens_max_px_sep) continue;
					unsigned bin_i = brgastro::get_bin_index<float>(d/pxfd_fact,sep_limits);

					increment_bin(bin_i,is_good_position(lens_xp+i,lens_yp+j,good_pixels));
					increment_bin(bin_i,is_good_position(lens_xp-i,lens_yp+j,good_pixels));
					increment_bin(bin_i,is_good_position(lens_xp+i,lens_yp-j,good_pixels));
					increment_bin(bin_i,is_good_position(lens_xp-i,lens_yp-j,good_pixels));

					increment_bin(bin_i,is_good_position(lens_xp+j,lens_yp+i,good_pixels));
					increment_bin(bin_i,is_good_position(lens_xp-j,lens_yp+i,good_pixels));
					increment_bin(bin_i,is_good_position(lens_xp+j,lens_yp-i,good_pixels));
					increment_bin(bin_i,is_good_position(lens_xp-j,lens_yp-i,good_pixels));


				}
			}

			// Set up result vector for this lens
			std::vector<float> lens_res(num_bins);
			for(unsigned i=0; i<num_bins; ++i)
			{
				if(static_cast<float>(lens_total_px_per_bin[i])==0)
					lens_res[i] = 1;
				else
					lens_res[i] = static_cast<float>(lens_good_px_per_bin[i])/
						static_cast<float>(lens_total_px_per_bin[i]);
			}

			// Add this to the field table map
			field_result[boost::lexical_cast<std::string>(lens_table_map.at("SeqNr")[lens_i])] =
					lens_res;
		}

		// Set up result vector for this field
		std::vector<float> field_res(num_bins);
		for(unsigned i=0; i<num_bins; ++i)
		{
			if(static_cast<float>(total_px_per_bin[i])==0)
				field_res[i] = 1;
			else
				field_res[i] = static_cast<float>(good_px_per_bin[i])/
				static_cast<float>(total_px_per_bin[i]);
		}

		// Get the good size of this field now
		const double good_size_upper = brgastro::square(rad_per_px)*num_good_pixels(good_pixels);
		const double good_size_lower = brgastro::square(rad_per_px)*num_good_pixels(good_positions);
		std::vector<double> size_measures;
		size_measures.push_back(good_size_lower);
		size_measures.push_back(good_size_upper);

		// Add this to the full result map and output this table

#ifdef _OPENMP
		#pragma omp critical(add_mask_test_results_to_map)
#endif
		{
			result[field_name_root] = field_res;
			field_sizes[field_name_root] = size_measures;

			//brgastro::print_table_map(lens_output_file_name,field_result);

			// Archive the good pixel map
			brgastro::binary_save(lens_pixel_map_file_name,good_pixels);

			std::cout << "Finished processing field " << field_name_root << "!\n";
		}
	}

	// Get the output file name
	std::stringstream ss("");
	ss << mask_directory << "position_corr.dat";
	std::string output_name = ss.str();

	brgastro::print_table_map(output_name,result);

	ss.str("");
	ss << mask_directory << "field_sizes.dat";
	std::string field_sizes_name = ss.str();

	brgastro::print_table_map(field_sizes_name,field_sizes);

	std::cout << "Done!\n";

	return 0;
}


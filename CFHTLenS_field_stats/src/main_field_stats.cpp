/**********************************************************************\
 @file main_field_stats.cpp
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


#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "brg/container/labeled_array.hpp"
#include "brg/file_access/ascii_table_map.hpp"
#include "brg/file_access/binary_archive.hpp"
#include "brg/file_access/open_file.hpp"
#include "brg/math/misc_math.hpp"
#include "brg/vector/limit_vector.hpp"

#include "brg_physics/units/unit_conversions.hpp"

// Magic values
const std::string data_directory = "/disk2/brg/git/CFHTLenS_cat/Data/";
const std::string field_directory = data_directory + "filtered_tables/";
const std::string fields_list = data_directory + "fields_list.txt";
const std::string lens_suffix = "_lens.dat";
const std::string source_suffix = "_source.dat";
const std::string pixel_map_suffix = "_lens_good_pixels.bin";

const std::string output_filename = data_directory + "field_stats.dat";

constexpr double z_bin_min = 0.2;
constexpr double z_bin_max = 1.3;
constexpr unsigned z_bins = 1;

constexpr double rad_per_px = 0.185965*brgastro::unitconv::asectorad;

inline size_t num_good_pixels(const std::vector<std::vector<bool>> & input)
{
	size_t result(0);

	unsigned ncol = input.size();
	assert(ncol>0);
	unsigned nrow = input[0].size();

	for(unsigned i=0; i<ncol; ++i)
	{
		for(unsigned j=0; j<nrow; ++j)
		{
			if(input[i][j]) ++result;
		}
	}

	return result;
}

/**
 * TODO (description)
 *
 * @param TODO (params to be passed at command-line)
 * @return
 */
int main( const int argc, const char *argv[] )
{
	
#ifdef _OPENMP
	omp_set_num_threads(5);
#endif

	brgastro::limit_vector<double> z_limits(z_bin_min,z_bin_max,z_bins);

	// Set up output data table
	brgastro::labeled_array<double> output_table;
	std::vector<std::string> output_header;

	std::stringstream ss;

	output_header.push_back("field_index");

	for(unsigned i=0; i<z_limits.num_bins(); ++i)
	{
		ss.str("");
		ss << "lens_count_z_" << z_limits.lower_limit(i);
		output_header.push_back(ss.str());
		ss.str("");
		ss << "source_count_z_" << z_limits.lower_limit(i);
		output_header.push_back(ss.str());
		ss.str("");
		ss << "lens_dens_z_" << z_limits.lower_limit(i);
		output_header.push_back(ss.str());
		ss.str("");
		ss << "source_dens_z_" << z_limits.lower_limit(i);
		output_header.push_back(ss.str());
	}

	output_table.set_labels(output_header);

	// Open and read in the fields list
	std::ifstream fi;
	brgastro::open_file_input(fi,fields_list);

	std::vector<std::string> field_names;
	std::string temp_field_name;

	while(fi>>temp_field_name)
	{
		field_names.push_back(temp_field_name);
	}

	unsigned num_fields = field_names.size();
	unsigned num_completed_fields = 0;

	//num_fields = 1;

	#ifdef _OPENMP
	#pragma omp parallel for schedule(dynamic)
	#endif
	for(unsigned field_i=0;field_i<num_fields;++field_i)
	{
		const std::string & field_name = field_names[field_i];
		const std::string field_name_root = field_name.substr(0,6);

		std::stringstream ss;

		ss.str("");
		ss << field_directory << field_name_root << pixel_map_suffix;
		const std::string pixel_map_file_name = ss.str();

		// Get the size of this field

		std::vector<std::vector<bool>> good_pixels;

		good_pixels = brgastro::binary_load_vector<std::vector<std::vector<bool>>>(pixel_map_file_name);
		// Get the good size of this field now
		size_t num_good = num_good_pixels(good_pixels);
		const double field_size = brgastro::square(rad_per_px)*num_good;

		// Get the lens file name
		ss.str("");
		ss << field_directory << field_name_root << lens_suffix;
		const std::string lens_file_name = ss.str();

		// Get the source file name
		ss.str("");
		ss << field_directory << field_name_root << source_suffix;
		const std::string source_file_name = ss.str();

		brgastro::labeled_array<double> lens_table;
		brgastro::labeled_array<double> source_table;

		lens_table.load(lens_file_name);
		source_table.load(source_file_name);

		// Loop over lenses and bin them by redshift
		Eigen::ArrayXd lens_count = Eigen::ArrayXd::Zero(z_bins+2);

		for(const auto & lens : lens_table.rows())
		{
			lens_count[z_limits.get_bin_index(lens.at_label("Z_B"))] += 1;
		}

		// Loop over sources and bin them by redshift
		Eigen::ArrayXd source_count = Eigen::ArrayXd::Zero(z_bins+2);

		for(const auto & source : source_table.rows())
		{
			source_count[z_limits.get_bin_index(source.at_label("Z_B"))] += 1;
		}

		Eigen::ArrayXd lens_dens = lens_count/field_size;
		Eigen::ArrayXd source_dens = lens_count/field_size;

		// Add this data to the output table

		std::vector<double> output_data;

		output_data.push_back(field_i);

		for(unsigned i=0; i<z_limits.num_bins(); ++i)
		{
			output_data.push_back(lens_count(i));
			output_data.push_back(source_count(i));
			output_data.push_back(lens_dens(i));
			output_data.push_back(source_dens(i));
		}

		#ifdef _OPENMP
		#pragma omp critical(update_field_stats_output)
		#endif
		{
			output_table.insert_row(output_data);
			std::cout << "Field " << field_name_root << " (#" <<
					++num_completed_fields << "/" << num_fields << ") complete.\n";
		}
	}

	output_table.save(output_filename);

	return 0;
}
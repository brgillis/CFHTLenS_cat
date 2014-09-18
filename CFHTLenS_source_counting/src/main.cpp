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
#include <vector>

#include <boost/lexical_cast.hpp>

#include "brg/file_access/open_file.hpp"
#include "brg/file_access/ascii_table_map.hpp"
#include "brg/physics/astro.h"
#include "brg/physics/lensing/pair_binner.h"
#include "brg/physics/lensing/pair_bins_summary.h"
#include "brg/physics/lensing/source_galaxy.h"
#include "brg/physics/sky_obj/galaxy.h"
#include "brg/vector/limit_vector.hpp"
#include "brg/vector/manipulations.hpp"

// Magic values
std::string fields_directory = "/disk2/brg/git/CFHTLenS_cat/Data/";
std::string fields_list = fields_directory + "fields_list.txt";
std::string output_table = fields_directory + "magnitude_counts.dat";

int main( const int argc, const char *argv[] )
{

	// Open and read in the fields list
	std::ifstream fi;
	brgastro::open_file_input(fi,fields_list);

	// Set up the redshift bins
	std::vector<double> z_bin_limits = brgastro::make_limit_vector<double>(0,2.0,0.1);
	std::vector<std::vector<double>> z_bin_magnitudes(z_bin_limits.size()-1);

	std::vector<std::string> field_names;
	std::string field_name;

	while(fi>>field_name)
	{
		field_names.push_back(field_name);
	}

	fi.close();

	// Load each field in turn and process it

	size_t num_fields = field_names.size();
	size_t num_processed = 0;

	for(size_t field_i=0;field_i<num_fields;++field_i)
	{
		std::string field_name_root = field_names[field_i].substr(0,6);

		try
		{
			// Get the source file names
			std::stringstream ss("");
			ss << fields_directory << "filtered_tables/" << field_name_root << "_source.dat";
			std::string source_input_name = ss.str();

			// Load in sources
			const brgastro::table_map_t<double> source_map(brgastro::load_table_map<double>(source_input_name));
			size_t num_sources = source_map.begin()->second.size();
			for(size_t i=0; i<num_sources; ++i)
			{
				z_bin_magnitudes[brgastro::get_bin_index(source_map.at("Z_B").at(i),
						z_bin_limits)].push_back(source_map.at("MAG_i").at(i));
			}
		}
		catch (const std::exception &e)
		{
			std::cerr << "Error processing field " << field_name_root << " (#" <<
					++num_processed << "/" << num_fields << ")!\n"
					<< e.what();
			continue;
		}
		std::cout << "Field " << field_name_root << " (#" <<
				++num_processed << "/" << num_fields << ") complete.\n";

	}

	// Pad the data with zeros so it can be printed
	z_bin_magnitudes = brgastro::pad(z_bin_magnitudes,0.);

	// Set up the header for the output table
	std::vector<std::string> header;
	for(auto it=z_bin_limits.begin();it!=z_bin_limits.end();++it)
	{
		header.push_back(boost::lexical_cast<std::string>(*it));
	}
	header.pop_back();

	brgastro::print_table(output_table,z_bin_magnitudes,header);

	return 0;
}

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
std::string fields_list = fields_directory + "fields_list.txt";
std::string output_table_root = fields_directory + "magnitude_hist_z";
std::string output_table_tail = ".dat";

int main( const int argc, const char *argv[] )
{

	// Open and read in the fields list
	std::ifstream fi;
	brgastro::open_file_input(fi,fields_list);

	// Set up the redshift bins
	std::vector<double> z_bin_limits = brgastro::make_limit_vector<double>(0.2,4.0,0.1);
	typedef boost::accumulators::accumulator_set<double,
			boost::accumulators::stats<
				boost::accumulators::tag::density >> hist_accum;
	typedef boost::accumulators::accumulator_set<double,
			boost::accumulators::stats<
				boost::accumulators::tag::density >,
					double > weighted_hist_accum;
	std::vector<hist_accum> z_bin_hists(z_bin_limits.size()-1,
			hist_accum(boost::accumulators::density_cache_size=10000,
					boost::accumulators::density_num_bins=100));
	std::vector<weighted_hist_accum> weighted_z_bin_hists(z_bin_limits.size()-1,
			weighted_hist_accum(boost::accumulators::density_cache_size=10000,
					boost::accumulators::density_num_bins=100));
	std::vector<unsigned long> z_bin_counts(z_bin_limits.size()-1,0);

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
				size_t z_i = brgastro::get_bin_index(source_map.at("Z_B").at(i),
						z_bin_limits);
				double mag = source_map.at("MAG_i").at(i);
				double weight = source_map.at("weight").at(i);
				z_bin_hists[z_i](mag);
				weighted_z_bin_hists[z_i](mag, boost::accumulators::weight =
						weight);
				++z_bin_counts[z_i];
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

	// Set up a histogram table for each redshift bin

	// Common header for all tables
	brgastro::header_t header;
	header.push_back("mag_bin_lower");
	header.push_back("weighted_count");

	auto z_hist_it = z_bin_hists.begin();
	auto z_w_hist_it = weighted_z_bin_hists.begin();
	auto z_count_it = z_bin_counts.begin();
	for(auto z_it=z_bin_limits.begin();z_it!=z_bin_limits.end(); ++z_it, ++z_hist_it, ++z_w_hist_it, ++z_count_it)
	{
		// Get the name for the table we'll output to
		std::string z_label = boost::lexical_cast<std::string>(brgastro::round_int(100 * *z_it));
		std::string output_file_name = output_table_root + z_label + output_table_tail;

		// Set up the data table, and make sure it has the needed columns
		brgastro::table_map_t<double> data;
		data["mag_bin_lower"];
		data["count"];
		data["weighted_count"];

		// Add each bin to the table
		auto hist = boost::accumulators::density(*z_hist_it);
		auto w_hist = boost::accumulators::weighted_density(*z_w_hist_it);
		auto w_h_it = w_hist.begin();
		for(auto h_it = hist.begin(); h_it!=hist.end(); ++h_it, ++w_h_it)
		{
			data["mag_bin_lower"].push_back(h_it->first);
			data["count"].push_back(*z_count_it * h_it->second);
			data["weighted_count"].push_back(*z_count_it * w_h_it->second);
		}

		brgastro::print_table_map(output_file_name,data);

	}

	return 0;
}

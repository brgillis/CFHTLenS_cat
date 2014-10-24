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

#include <boost/lexical_cast.hpp>

#include "brg/external/sgsmooth.h"
#include "brg/file_access/ascii_table_map.hpp"
#include "brg/file_access/open_file.hpp"
#include "brg/file_access/table_typedefs.hpp"
#include "brg/math/misc_math.hpp"
#include "brg/physics/astro.h"
#include "brg/physics/lensing/magnification/mag_global_values.h"
#include "brg/physics/lensing/pair_binner.h"
#include "brg/physics/lensing/pair_bins_summary.h"
#include "brg/physics/lensing/source_galaxy.h"
#include "brg/physics/sky_obj/galaxy.h"
#include "brg/vector/elementwise_functions.hpp"
#include "brg/vector/limit_vector.hpp"
#include "brg/vector/manipulations.hpp"

// Magic values
const std::string fields_directory = "/disk2/brg/git/CFHTLenS_cat/Data/";
const std::string fields_list = fields_directory + "fields_list.txt";
const std::string output_table_root = fields_directory + "magnitude_hist_z";
const std::string g_output_table_root = fields_directory + "magnitude_hist_gz";
const std::string output_table_tail = ".dat";
const std::string field_size_file_name = "/disk2/brg/git/CFHTLenS_cat/Data/masks/field_sizes.dat";

constexpr short unsigned sg_window = 5;
constexpr short unsigned sg_deg = 3;

constexpr short unsigned cache_size = 1000;
constexpr short unsigned num_mag_bins = 10;

int main( const int argc, const char *argv[] )
{

	// Open and read in the fields list
	std::ifstream fi;
	brgastro::open_file_input(fi,fields_list);

	// Set up the redshift bins
	std::vector<double> z_bin_limits = brgastro::make_limit_vector<double>(brgastro::mag_z_min,
			brgastro::mag_z_max,brgastro::mag_z_step);

	std::vector<double> mag_bin_limits = brgastro::make_limit_vector<double>(brgastro::mag_m_min,
			brgastro::mag_m_max,brgastro::mag_m_step);

	size_t num_z_bins = mag_bin_limits.size()-1;
	size_t num_mag_bins = mag_bin_limits.size()-1;

	typedef std::vector<unsigned long> int_hist;
	typedef std::vector<long double> double_hist;

	// Hists for exactly in bin
	std::vector<int_hist> z_bin_hists(num_z_bins,
			int_hist(num_mag_bins,0));
	std::vector<double_hist> weighted_z_bin_hists(num_z_bins,
			double_hist(num_mag_bins,0));
	int_hist z_bin_counts(num_z_bins,0);

	// Hists for in this bin or greater
	std::vector<int_hist> gz_bin_hists(num_z_bins,
			int_hist(num_mag_bins,0));
	std::vector<double_hist> weighted_gz_bin_hists(num_z_bins,
			double_hist(num_mag_bins,0));
	int_hist gz_bin_counts(num_z_bins,0);

	std::vector<std::string> field_names;
	std::string field_name;

	std::vector<double> field_sizes;
	const auto field_size_table = brgastro::load_table_map<double>(field_size_file_name);

	while(fi>>field_name)
	{
		field_names.push_back(field_name);
	}

	fi.close();

	// Load each field in turn and process it

	size_t num_fields = field_names.size();
	size_t num_processed = 0;

	//num_fields = 1;

	for(size_t field_i=0;field_i<num_fields;++field_i)
	{
		std::string field_name_root = field_names[field_i].substr(0,6);

		try
		{
			// Get the field size
			const auto size_measurements = field_size_table.at(field_name_root);
			//const double field_size = brgastro::mean(size_measurements);
			//const double field_size = size_measurements[0];
			const double field_size = size_measurements[1];
			field_sizes.push_back(field_size);

			// Get the source file names
			std::stringstream ss("");
			ss << fields_directory << "filtered_tables/" << field_name_root << "_source.dat";
			std::string source_input_name = ss.str();

			// Load in sources
			//const brgastro::table_map_t<double> source_map(brgastro::load_table_map<double>(source_input_name));
			brgastro::table_map_t<double> source_map;
			source_map = (brgastro::load_table_map<double>(source_input_name));
			size_t num_sources = source_map.begin()->second.size();
			for(size_t i=0; i<num_sources; ++i)
			{
				const auto & source_z = source_map.at("Z_B").at(i);
				size_t z_i = brgastro::get_bin_index(source_z,
						z_bin_limits);
				const double & mag = source_map.at("MAG_r").at(i);
				const double & weight = source_map.at("weight").at(i);

				if(brgastro::outside_limits(mag,mag_bin_limits)) continue;

				size_t mag_i = brgastro::get_bin_index(mag,
						mag_bin_limits);

				// Add to the specific bin
				++z_bin_hists[z_i][mag_i];
				weighted_z_bin_hists[z_i][mag_i]+=weight;
				++z_bin_counts[z_i];

				// Add to the cumulative bins
				for(size_t j=0; j<=z_i; ++j)
				{
					++gz_bin_hists[j][mag_i];
					weighted_gz_bin_hists[j][mag_i]+=weight;
					++gz_bin_counts[j];
				}
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
#ifndef NDEBUG
		break; // debug line
#endif

	}

	// Get the total survey size
	const double survey_size = brgastro::sum(field_sizes);

	// Set up and print histogram tables

	// Specific redshift tables

	auto z_hist_it = z_bin_hists.begin();
	auto z_w_hist_it = weighted_z_bin_hists.begin();
	auto z_count_it = z_bin_counts.begin();
	for(auto z_it=z_bin_limits.begin();z_hist_it!=z_bin_hists.end(); ++z_it, ++z_hist_it, ++z_w_hist_it, ++z_count_it)
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
		auto h_it = z_hist_it->begin();
		auto w_h_it = z_w_hist_it->begin();
		auto mag_it = mag_bin_limits.begin();
		for(; h_it!=z_hist_it->end(); ++h_it, ++w_h_it, ++mag_it)
		{
			data["mag_bin_lower"].push_back(*mag_it);
			data["count"].push_back(*h_it / survey_size);
			data["weighted_count"].push_back(*w_h_it / survey_size);
		}

		double mag_step = data["mag_bin_lower"].at(1)-data["mag_bin_lower"].at(0);

		// Replace any zeros with 0.1 so we can take the log safely
		auto rep_func = [] (const double & v) {if(v>0) return v; else return 0.1;};

		// Now calculate the smoothed count
		auto log_count = brgastro::divide(brgastro::log(brgastro::apply(rep_func,data["count"])),
				std::log(10.));
		auto smoothed_log = brgastro::sg_smooth(log_count,sg_window,sg_deg);
		data["smoothed_count"] = brgastro::pow(10.,smoothed_log);

		// Check for a systematic over or underestimation from smoothing
		double smooth_correction_factor = brgastro::sum(data["count"])/
				brgastro::sum(data["smoothed_count"]);

		data["smoothed_count"] = brgastro::multiply(data["smoothed_count"],smooth_correction_factor);

		data["smoothed_alpha"] = brgastro::add(brgastro::multiply(
				brgastro::sg_derivative(smoothed_log,sg_window,sg_deg),
				2.5/mag_step),0);

		brgastro::print_table_map(output_file_name,data);

	}

	// Cumulative redshift tables

	z_hist_it = gz_bin_hists.begin();
	z_w_hist_it = weighted_gz_bin_hists.begin();
	z_count_it = gz_bin_counts.begin();
	for(auto z_it=z_bin_limits.begin();z_hist_it!=gz_bin_hists.end(); ++z_it, ++z_hist_it, ++z_w_hist_it, ++z_count_it)
	{
		// Get the name for the table we'll output to
		std::string z_label = boost::lexical_cast<std::string>(brgastro::round_int(100 * *z_it));
		std::string output_file_name = g_output_table_root + z_label + output_table_tail;

		// Set up the data table, and make sure it has the needed columns
		brgastro::table_map_t<double> data;
		data["mag_bin_lower"];
		data["count"];
		data["weighted_count"];

		// Add each bin to the table
		auto h_it = z_hist_it->begin();
		auto w_h_it = z_w_hist_it->begin();
		auto mag_it = mag_bin_limits.begin();
		for(; h_it!=z_hist_it->end(); ++h_it, ++w_h_it, ++mag_it)
		{
			data["mag_bin_lower"].push_back(*mag_it);
			data["count"].push_back(*h_it / survey_size);
			data["weighted_count"].push_back(*w_h_it / survey_size);
		}

		double mag_step = data["mag_bin_lower"].at(1)-data["mag_bin_lower"].at(0);

		// Replace any zeros with 0.1 so we can take the log safely
		auto rep_func = [] (const double & v) {if(v>0) return v; else return 0.1;};

		// Now calculate the smoothed count
		auto log_count = brgastro::divide(brgastro::log(brgastro::apply(rep_func,data["count"])),
				std::log(10.));
		auto smoothed_log = brgastro::sg_smooth(log_count,sg_window,sg_deg);
		data["smoothed_count"] = brgastro::pow(10.,smoothed_log);

		// Check for a systematic over or underestimation from smoothing
		double smooth_correction_factor = brgastro::sum(data["count"])/
				brgastro::sum(data["smoothed_count"]);

		data["smoothed_count"] = brgastro::multiply(data["smoothed_count"],smooth_correction_factor);

		data["smoothed_alpha"] = brgastro::add(brgastro::multiply(
				brgastro::sg_derivative(smoothed_log,sg_window,sg_deg),
				2.5/mag_step),0);

		brgastro::print_table_map(output_file_name,data);

	}
	return 0;
}

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

#include "brg/file_access/open_file.hpp"
#include "brg/file_access/ascii_table_map.hpp"
#include "brg/physics/astro.h"
#include "brg/physics/lensing/magnification/expected_count_cache.h"
#include "brg/physics/lensing/magnification/expected_count_derivative_cache.h"
#include "brg/physics/lensing/magnification/mag_signal_integral_cache.h"
#include "brg/physics/lensing/magnification/mag_weight_integral_cache.h"
#include "brg/physics/lensing/magnification/mag_global_values.h"
#include "brg/physics/lensing/pair_binner.h"
#include "brg/physics/lensing/pair_bins_summary.h"
#include "brg/physics/lensing/source_galaxy.h"
#include "brg/physics/sky_obj/galaxy.h"
#include "brg/physics/units/unitconv_map.hpp"

#include "gg_lensing_config.h"
#include "pass_configs_to_binner.h"

// Magic values
std::string fields_directory = "/disk2/brg/git/CFHTLenS_cat/Data/";
std::string fields_list = fields_directory + "fields_list.txt";
std::string output_table = fields_directory + "gg_lensing_signal.dat";
std::string output_data = fields_directory + "gg_lensing_data.dat";

std::string expected_count_cache_output_file = fields_directory + "ex_count_cache.dat";
std::string expected_count_derivative_cache_output_file = fields_directory + "alpha_cache.dat";
std::string mag_signal_integral_cache_output_file = fields_directory + "mag_sig_integral_cache.dat";
std::string mag_weight_integral_cache_output_file = fields_directory + "mag_W_integral_cache.dat";

int main( const int argc, const char *argv[] )
{
	// Set up the caches before we get to the parallel section, so they can be calculated in parallel
	brgastro::expected_count_cache().print(expected_count_cache_output_file);
	brgastro::expected_count_derivative_cache().print(expected_count_derivative_cache_output_file);
	brgastro::mag_weight_integral_cache().print(mag_weight_integral_cache_output_file);
	brgastro::mag_signal_integral_cache().print(mag_signal_integral_cache_output_file);

	constexpr size_t batch_size = 1;
	const gg_lensing_config config(argc,argv);

	// Set up the bins summary
	brgastro::pair_bins_summary bins_summary(pass_configs_to_binner(config));

	// Check if we're using preloaded or saved data
	if(config.use_precalculated_data)
	{
		bins_summary.load(config.precalculated_data_filename);
	}
	else
	{

		// Open and read in the fields list
		std::ifstream fi;
		brgastro::open_file_input(fi,fields_list);

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

		#pragma omp parallel for
		for(size_t field_i=0;field_i<num_fields;++field_i)
		{
			brgastro::pair_bins_summary field_bins_summary(pass_configs_to_binner(config));
			brgastro::pair_binner lens_binner(pass_configs_to_binner(config));
			std::string field_name_root = field_names[field_i].substr(0,6);

			try
			{
				// Get the lens and source file names
				std::stringstream ss("");
				ss << fields_directory << "filtered_tables/" << field_name_root << "_lens.dat";
				std::string lens_input_name = ss.str();

				ss.str("");
				ss << fields_directory << "filtered_tables/" << field_name_root << "_source.dat";
				std::string source_input_name = ss.str();

				// Set up vectors
				std::vector<brgastro::galaxy> lens_galaxies;
				std::vector<brgastro::source_galaxy> source_galaxies;

				// Load in lenses
				const brgastro::table_map_t<double> lens_map(brgastro::load_table_map<double>(lens_input_name));
				size_t num_lenses = lens_map.begin()->second.size();
				for(size_t i=0; i<num_lenses; ++i)
				{
					brgastro::galaxy lens;
					lens.set_z(lens_map.at("Z_B").at(i));
					lens.set_ra(lens_map.at("ALPHA_J2000").at(i));
					lens.set_dec(lens_map.at("DELTA_J2000").at(i));
					lens.stellar_mass = lens_map.at("LP_log10_SM_MED").at(i);
					lens.imag = lens_map.at("MAG_i").at(i);
					lens.set_index(i);

					lens_galaxies.push_back(lens);
				}

				// Load in sources
				const brgastro::table_map_t<double> source_map(brgastro::load_table_map<double>(source_input_name));
				size_t num_sources = source_map.begin()->second.size();
				for(size_t i=0; i<num_sources; ++i)
				{
					source_galaxies.push_back(
							brgastro::source_galaxy(source_map.at("ALPHA_J2000").at(i),
									source_map.at("DELTA_J2000").at(i),
									source_map.at("Z_B").at(i),
									source_map.at("e1").at(i), source_map.at("e2").at(i), 0,
									source_map.at("LP_log10_SM_MED").at(i), source_map.at("MAG_i").at(i)));
					source_galaxies.back().set_weight(source_map.at("weight").at(i));
				}

				// Find lens-source pairs and add them to the binner
				for(size_t lens_i=0; lens_i<num_lenses; ++lens_i)
				{
					brgastro::galaxy & lens = lens_galaxies[lens_i];

					// Check if the lens fits somewhere within the binner's limits
					if(!lens_binner.binnable(lens)) continue;

					lens_binner.add_lens_id(lens.index(),lens.m(),lens.z(),lens.mag());

					BRG_ANGLE max_angle_sep = brgastro::afd(config.R_max,lens.z());

					for(size_t source_i=0; source_i<num_sources; ++source_i)
					{
						brgastro::source_galaxy & source = source_galaxies[source_i];

						// Check that the lens is sufficiently in front of the source
						if(lens.z() >= source.z() - config.z_buffer) continue;

						// Check that the source is within our magnitude limitations
						if(source.mag() > brgastro::mag_m_max) continue;
						if(source.mag() < brgastro::mag_m_min) continue;

						// Check against maximum angular separation in ra and dec simply first for speed
						auto ddec = std::fabs(lens.dec()-source.dec());
						if(ddec>max_angle_sep) continue;
						double cosdec = std::cos(lens.dec());
						auto dra = std::fabs(lens.ra()-source.ra())*cosdec;
						if(dra>max_angle_sep) continue;

						BRG_DISTANCE R = brgastro::dfa(brgastro::dist2d(dra,ddec),lens.z());

						if(R <= config.R_max)
						{
							lens_binner.add_pair(brgastro::lens_source_pair(&lens,&source));
						}
					}

					if((lens_i+1) % batch_size == 0)
					{
						// Add this binner to the field summary and empty it
						field_bins_summary += lens_binner;
						lens_binner.empty();
					}
				}
				field_bins_summary += lens_binner;
			}
			catch (const std::exception &e)
			{

				#pragma omp critical(CFHTLenS_gg_lens_combine_fields)
				{
					std::cerr << "Error processing field " << field_name_root << " (#" <<
							++num_processed << "/" << num_fields << ")!\n"
							<< e.what();
				}
				continue;

			}

			#pragma omp critical(CFHTLenS_gg_lens_combine_fields)
			{
				try
				{
					bins_summary += field_bins_summary;
					std::cout << "Field " << field_name_root << " (#" <<
							++num_processed << "/" << num_fields << ") complete.\n";
				}
				catch (const std::exception &e)
				{
					std::cerr << "Error combining data from field " << field_name_root << " (#" <<
							++num_processed << "/" << num_fields << ")!\n"
							<< e.what();
				}
			}

		}

		// Save the data so we can load it in again without needing to rerun
		bins_summary.save(output_data);
	}

	// Set up the units we want data to be output in
	brgastro::unitconv_map u_map;

	u_map["R_min"] = u_map["R_max"] = u_map["R_mean"] = brgastro::unitconv::kpctom;
	u_map["m_min"] = u_map["m_max"] = u_map["m_mean"] = brgastro::unitconv::Msuntokg;
	u_map["dS_t_mean"] = u_map["dS_t_stddev"] = u_map["dS_t_stderr"] =
			u_map["dS_x_mean"] = u_map["dS_x_stddev"] = u_map["dS_x_stderr"] =
					brgastro::unitconv::Msuntokg/brgastro::square(brgastro::unitconv::pctom);

	bins_summary.print_bin_data(output_table,u_map);
	bins_summary.print_bin_data(std::cout,u_map);

	return 0;
}

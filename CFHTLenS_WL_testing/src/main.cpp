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
#include <brg/file_access/ascii_table_map.hpp>
#include <brg/file_access/open_file.hpp>

#include "brg/vector/elementwise_functions.hpp"
#include "brg/vector/limit_vector.hpp"

#include "brg_lensing/magnification/expected_count_cache.h"
#include "brg_lensing/magnification/expected_count_derivative_cache.h"
#include "brg_lensing/magnification/mag_calibration_cache.h"
#include "brg_lensing/magnification/mag_signal_integral_cache.h"
#include "brg_lensing/magnification/mag_weight_integral_cache.h"
#include "brg_lensing/magnification/mag_global_values.h"
#include "brg_lensing/pair_binner.h"
#include "brg_lensing/pair_bins_summary.h"
#include "brg_lensing/source_galaxy.h"

#include "brg_physics/astro.h"
#include "brg_physics/sky_obj/galaxy.h"
#include "brg_physics/units/unitconv_map.hpp"

#include "gg_lensing_config.h"
#include "pass_configs_to_binner.h"

#undef USE_CALIBRATION_LENSES

#define USE_MOCK_LENSES
#undef USE_MOCK_SOURCES

// Magic values
std::string fields_directory = "/disk2/brg/git/CFHTLenS_cat/Data/";
std::string fields_list = fields_directory + "fields_list.txt";

#ifdef USE_CALIBRATION_LENSES

std::string output_table = fields_directory + "gg_calibration_lensing_signal.dat";
std::string output_data = fields_directory + "gg_calibration_lensing_data.dat";

std::string lens_root = "_calibration_lens.dat";
std::string lens_unmasked_frac_root = "_calibration_lens_mask_frac.dat";

#else // #ifdef USE_CALIBRATION_LENSES

#ifdef USE_MOCK_LENSES

std::string output_table = fields_directory + "gg_mock_lensing_signal.dat";
std::string output_data = fields_directory + "gg_mock_lensing_data.dat";

std::string lens_root = "_small_mock_lens.dat";
std::string lens_unmasked_frac_root = "_small_mock_lens_mask_frac.dat";

#else // #ifdef USE_MOCK_LENSES

std::string output_table = fields_directory + "gg_lensing_signal.dat";
std::string output_data = fields_directory + "gg_lensing_data.dat";

std::string lens_root = "_lens.dat";
std::string lens_unmasked_frac_root = "_lens_mask_frac.dat";

#endif // #ifdef USE_MOCK_LENSES // #else

#endif // #ifdef USE_CALIBRATION_LENSES // #else

#ifdef USE_MOCK_SOURCES
std::string source_root = "_small_mock_source.dat";
#else
std::string source_root = "_source.dat";
#endif

std::string expected_count_cache_output_file = fields_directory + "ex_count_cache.dat";
std::string expected_count_derivative_cache_output_file = fields_directory + "alpha_cache.dat";
std::string mag_signal_integral_cache_output_file = fields_directory + "mag_sig_integral_cache.dat";
std::string mag_weight_integral_cache_output_file = fields_directory + "mag_W_integral_cache.dat";
std::string mag_calibration_cache_output_file = fields_directory + "mag_calibration_cache.dat";

constexpr double mag_fudge_shift = 0;

int main( const int argc, const char *argv[] )
{
	// Set up the caches before we get to the parallel section, so they can be calculated in parallel
	brgastro::expected_count_cache().print(expected_count_cache_output_file);
	brgastro::expected_count_derivative_cache().print(expected_count_derivative_cache_output_file);
	brgastro::mag_weight_integral_cache().print(mag_weight_integral_cache_output_file);
	brgastro::mag_signal_integral_cache().print(mag_signal_integral_cache_output_file);
#ifndef USE_CALIBRATION_LENSES
	brgastro::mag_calibration_cache().print(mag_calibration_cache_output_file);
#endif // #ifndef USE_CALIBRATION_LENSES

	constexpr size_t batch_size = 1000000;
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

		//num_fields = 1;

		#ifdef _OPENMP
		#pragma omp parallel for schedule(dynamic)
		#endif
		for(size_t field_i=0;field_i<num_fields;++field_i)
		{
			brgastro::pair_bins_summary field_bins_summary(pass_configs_to_binner(config));
			brgastro::pair_binner lens_binner(pass_configs_to_binner(config));
			std::string field_name_root = field_names[field_i].substr(0,6);

			try
			{
				// Get the lens and source file names
				std::stringstream ss("");
				ss << fields_directory << "filtered_tables/" << field_name_root << lens_root;
				std::string lens_input_name = ss.str();

				ss.str("");
				ss << fields_directory << "filtered_tables/" << field_name_root << lens_unmasked_frac_root;
				std::string lens_unmasked_name = ss.str();

				ss.str("");
				ss << fields_directory << "filtered_tables/" << field_name_root << source_root;
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
					lens.set_ra(lens_map.at("ra_radians").at(i));
					lens.set_dec(lens_map.at("dec_radians").at(i));
					lens.stellar_mass = lens_map.at("Mstel_kg").at(i);
					lens.imag = lens_map.at("MAG_i").at(i);
					lens.set_index(lens_map.at("SeqNr").at(i));

					// Check if the lens fits somewhere within the binner's limits
					if(!lens_binner.binnable(lens)) continue;

					lens_galaxies.push_back(std::move(lens));
				}

				// Load the masked fraction table
				const brgastro::table_map_t<double> lens_unmasked_frac_map(
						brgastro::load_table_map<double>(lens_unmasked_name));
				brgastro::limit_vector<double> unmasked_frac_bin_limits;

				unmasked_frac_bin_limits.reconstruct_from_bin_mids(brgastro::multiply(lens_unmasked_frac_map.at("bin_mid_kpc"),
																						brgastro::unitconv::kpctom));

				// Load in sources
				const brgastro::table_map_t<double> source_map(brgastro::load_table_map<double>(source_input_name));
				size_t num_sources = source_map.begin()->second.size();

				for(size_t i=0; i<num_sources; ++i)
				{
					brgastro::source_galaxy source(source_map.at("ra_radians").at(i),
										source_map.at("dec_radians").at(i),
										source_map.at("Z_B").at(i),
										source_map.at("e1").at(i), source_map.at("e2").at(i), 0,
										source_map.at("Mstel_kg").at(i), source_map.at("MAG_r").at(i)+mag_fudge_shift);

					source.set_weight(source_map.at("weight").at(i));
					source.set_index(source_map.at("SeqNr").at(i));

					// Check that this source is valid for either shear or magnification

					const auto & mag = source.mag();
					const auto & z = source.z();
					const auto & shear_weight = source.weight();

					if( (shear_weight > 0) || // If it's valid for shear
						((mag>=brgastro::mag_m_min) && (mag<brgastro::mag_m_max) &&
						 (z>=brgastro::mag_z_min) && (z<brgastro::mag_z_max)) ) // or it's valid for magnification
					{
						source_galaxies.push_back( std::move(source) ); // Add it to the list of sources
					}
				}

				// Find lens-source pairs and add them to the binner
				size_t lens_i = 0;
				for(const auto & lens : lens_galaxies)
				{

					// Get the unmasked fraction for this lens and set it
					lens_binner.set_unmasked_fractions( unmasked_frac_bin_limits,
							lens_unmasked_frac_map.at(boost::lexical_cast<std::string>(lens.index())));

					lens_binner.add_lens_id(lens.index(),lens.m(),lens.z(),lens.mag());

					const BRG_ANGLE max_angle_sep = brgastro::afd(config.R_max,lens.z());

					// And loop over all sources

					//unsigned counter = 0; //!!

					for(const auto & source : source_galaxies)
					{

						// Check that the lens is sufficiently in front of the source
						if(lens.z() > source.z() - config.z_buffer + std::numeric_limits<double>::epsilon()) continue;

//						if(((source.mag()>=brgastro::mag_m_min) && (source.mag()<brgastro::mag_m_max) &&
//						 (source.z()>=1.14) && (source.z()<brgastro::mag_z_max)) ) ++counter; //!!

						// Check against maximum angular separation in ra and dec simply first for speed
						auto ddec = std::fabs(lens.dec()-source.dec());
						if(ddec>max_angle_sep) continue;
						double cosdec = std::cos((lens.dec()+source.dec())/2);
						auto dra = std::fabs(lens.ra()-source.ra())*cosdec;
						if(dra>max_angle_sep) continue;

						double da = brgastro::dist2d(dra,ddec);

						if(da <= max_angle_sep)
						{
							lens_binner.add_pair(brgastro::lens_source_pair(&lens,&source));
						}

					}
					//std::cout << counter << std::endl; //!!

					if((++lens_i) % batch_size == 0)
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

#ifdef _OPENMP
				#pragma omp critical(CFHTLenS_gg_lens_combine_fields)
				{
					std::cerr << "Error processing field " << field_name_root << " (#" <<
							++num_processed << "/" << num_fields << ")!\n"
							<< e.what() << std::endl;
				}
				continue;
#else
				throw;
#endif

			}

#ifdef _OPENMP
			#pragma omp critical(CFHTLenS_gg_lens_combine_fields)
#endif
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

	u_map["R_min"] = u_map["R_max"] = u_map["shear_R_mean"] = u_map["magf_R_mean"] =
		brgastro::unitconv::kpctom;
	u_map["m_min"] = u_map["m_max"] = u_map["shear_lens_m_mean"] = u_map["magf_lens_m_mean"] =
		brgastro::unitconv::Msuntokg;
	u_map["dS_t_mean"] = u_map["dS_t_stddev"] = u_map["dS_t_stderr"] =
		u_map["dS_x_mean"] = u_map["dS_x_stddev"] = u_map["dS_x_stderr"] =
		u_map["model_dS_t"] = u_map["Sigma"] = u_map["Sigma_stderr"] =
		u_map["model_Sigma"] = u_map["shear_Sigma_crit"] = u_map["magf_Sigma_crit"] =
			brgastro::unitconv::Msuntokg/brgastro::square(brgastro::unitconv::pctom);
	u_map["magf_area"] = brgastro::square(brgastro::unitconv::asectorad);

	bins_summary.print_bin_data(output_table,u_map);
	bins_summary.print_bin_data(std::cout,u_map);

	return 0;
}

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
#include <tuple>
#include <vector>

#include <boost/lexical_cast.hpp>

#include <Eigen/Core>

#include "brg/container/coerce.hpp"
#include "brg/file_access/open_file.hpp"
#include "brg/file_access/ascii_table_map.hpp"
#include "brg/vector/limit_vector.hpp"
#include "brg_physics/units/unit_conversions.hpp"
#include "brg_lensing/magnification/mag_global_values.h"

#define SMALL_MOCKS

#define LENSING_STYLE

#ifdef LENSING_STYLE
#include "brg_physics/lensing_correlation_function_estimator.h"
#else
#include "brg_physics/correlation_function_estimator.h"
#endif

// Magic values
std::string fields_directory = "/disk2/brg/git/CFHTLenS_cat/Data/";
std::string fields_list = fields_directory + "fields_list.txt";

const std::string lens_root = "_lens.dat";
const std::string source_root = "_source.dat";

#ifdef SMALL_MOCKS
#ifdef LENSING_STYLE
std::string output_name = fields_directory + "lensing_corr_funcs_quick.dat";
#else
std::string output_name = fields_directory + "auto_corr_funcs_quick.dat";
#endif
const std::string lens_mock_root = "_small_mock_lens.dat";
const std::string source_mock_root = "_small_mock_source.dat";
#else
#ifdef LENSING_STYLE
std::string output_name = fields_directory + "lensing_corr_funcs.dat";
#else
std::string output_name = fields_directory + "auto_corr_funcs.dat";
#endif
const std::string lens_mock_root = "_mock_lens.dat";
const std::string source_mock_root = "_mock_source.dat";
#endif

constexpr double lens_z_min = 0.5;
constexpr double lens_z_max = 0.6;
constexpr double lens_m_min = 1e10;
constexpr double lens_m_max = 1e11;

#ifdef LENSING_STYLE
constexpr double source_z_min = brgastro::mag_z_min;
constexpr double source_z_max = brgastro::mag_z_max;
constexpr double z_buffer = 0.3;
#else
constexpr double source_z_min = lens_z_min;
constexpr double source_z_max = lens_z_max;
#endif

constexpr double r_min = 0.00001616; // About 10 kpc at redshift 0.2
constexpr double r_max = 0.00291; // About 2000 kpc at redshift 0.2
constexpr size_t r_steps = 200;

int main( const int argc, const char *argv[] )
{
	// Set up global data
	brgastro::limit_vector<double> r_limits(r_min,r_max,r_steps);
	Eigen::ArrayXd monopole_corr_func_sum = Eigen::ArrayXd::Zero(r_steps);
	Eigen::ArrayXd dipole_1_corr_func_sum = Eigen::ArrayXd::Zero(r_steps);
	Eigen::ArrayXd dipole_2_corr_func_sum = Eigen::ArrayXd::Zero(r_steps);
	Eigen::ArrayXd quadrupole_1_corr_func_sum = Eigen::ArrayXd::Zero(r_steps);
	Eigen::ArrayXd quadrupole_2_corr_func_sum = Eigen::ArrayXd::Zero(r_steps);
	Eigen::ArrayXd octopole_1_corr_func_sum = Eigen::ArrayXd::Zero(r_steps);
	Eigen::ArrayXd octopole_2_corr_func_sum = Eigen::ArrayXd::Zero(r_steps);

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
		std::string field_name_root = field_names[field_i].substr(0,6);

		Eigen::ArrayXd monopole_corr_func = Eigen::ArrayXd::Zero(r_steps);
		Eigen::ArrayXd dipole_1_corr_func = Eigen::ArrayXd::Zero(r_steps);
		Eigen::ArrayXd dipole_2_corr_func = Eigen::ArrayXd::Zero(r_steps);
		Eigen::ArrayXd quadrupole_1_corr_func = Eigen::ArrayXd::Zero(r_steps);
		Eigen::ArrayXd quadrupole_2_corr_func = Eigen::ArrayXd::Zero(r_steps);
		Eigen::ArrayXd octopole_1_corr_func = Eigen::ArrayXd::Zero(r_steps);
		Eigen::ArrayXd octopole_2_corr_func = Eigen::ArrayXd::Zero(r_steps);

		try
		{
			// Get the lens and source file names
			std::stringstream ss("");
			ss << fields_directory << "filtered_tables/" << field_name_root << lens_root;
			std::string lens_input_name = ss.str();

			ss.str("");
			ss << fields_directory << "filtered_tables/" << field_name_root << source_root;
			std::string source_input_name = ss.str();

			ss.str("");
			ss << fields_directory << "filtered_tables/" << field_name_root << lens_mock_root;
			std::string mock_lens_input_name = ss.str();

			ss.str("");
			ss << fields_directory << "filtered_tables/" << field_name_root << source_mock_root;
			std::string mock_source_input_name = ss.str();

			// Set up vectors
			#ifdef LENSING_STYLE
			std::vector<std::tuple<double,double,double>> lens_positions;
			std::vector<std::tuple<double,double,double>> source_positions;
			std::vector<std::tuple<double,double,double>> mock_lens_positions;
			std::vector<std::tuple<double,double,double>> mock_source_positions;
			#else
			std::vector<std::pair<double,double>> lens_positions;
			std::vector<std::pair<double,double>> source_positions;
			std::vector<std::pair<double,double>> mock_lens_positions;
			std::vector<std::pair<double,double>> mock_source_positions;
			#endif

			// Load in lenses
			{
				const brgastro::table_map_t<double> lens_map(brgastro::load_table_map<double>(lens_input_name));
				size_t num_lenses = lens_map.begin()->second.size();
				for(size_t i=0; i<num_lenses; ++i)
				{
					const double & z = lens_map.at("Z_B").at(i);
					if((z<lens_z_min)||(z>lens_z_max)) continue;
					const double & m = lens_map.at("Mstel_kg").at(i)*brgastro::unitconv::kgtoMsun;
					if((m<lens_m_min)||(m>lens_m_max)) continue;

//					const double & T = lens_map.at("T_B").at(i);
//					if((T<brgastro::mag_lens_T_min)||(T>brgastro::mag_lens_T_max)) continue;

					#ifdef LENSING_STYLE
					lens_positions.push_back(std::tie(
							lens_map.at("ra_radians").at(i),
							lens_map.at("dec_radians").at(i),
							lens_map.at("Z_B").at(i)));
					#else
					lens_positions.push_back(std::make_pair(
							lens_map.at("ra_radians").at(i),
							lens_map.at("dec_radians").at(i)));
					#endif
				}
			}

			// Load in sources
			{
				const brgastro::table_map_t<double> source_map(brgastro::load_table_map<double>(source_input_name));
				size_t num_sources = source_map.begin()->second.size();
				for(size_t i=0; i<num_sources; ++i)
				{
					const double & z = source_map.at("Z_B").at(i);
					if((z<source_z_min)||(z>source_z_max)) continue;
					const double & mag = source_map.at("MAG_r").at(i);
					if((mag<brgastro::mag_m_min)||(mag>brgastro::mag_m_max)) continue;

					#ifdef LENSING_STYLE
					source_positions.push_back(std::tie(
						source_map.at("ra_radians").at(i),
						source_map.at("dec_radians").at(i),
						source_map.at("Z_B").at(i)));
					#else
					source_positions.push_back(std::make_pair(
						source_map.at("ra_radians").at(i),
						source_map.at("dec_radians").at(i)));
					#endif
				}
			}

			// Load in mock lenses
			{
				const brgastro::table_map_t<double> mock_lens_map(brgastro::load_table_map<double>(mock_lens_input_name));
				size_t num_mock_lenses = mock_lens_map.begin()->second.size();
				for(size_t i=0; i<num_mock_lenses; ++i)
				{
					const double & z = mock_lens_map.at("Z_B").at(i);
					if((z<lens_z_min)||(z>lens_z_max)) continue;

//					const double & T = mock_lens_map.at("T_B").at(i);
//					if((T<brgastro::mag_lens_T_min)||(T>brgastro::mag_lens_T_max)) continue;

					#ifdef LENSING_STYLE
					mock_lens_positions.push_back(std::tie(
						mock_lens_map.at("ra_radians").at(i),
						mock_lens_map.at("dec_radians").at(i),
						mock_lens_map.at("Z_B").at(i)));
					#else
					mock_lens_positions.push_back(std::make_pair(
						mock_lens_map.at("ra_radians").at(i),
						mock_lens_map.at("dec_radians").at(i)));
					#endif
				}
			}

			// Load in mock sources
			{
				const brgastro::table_map_t<double> mock_source_map(brgastro::load_table_map<double>(mock_source_input_name));
				size_t num_mock_sources = mock_source_map.begin()->second.size();
				for(size_t i=0; i<num_mock_sources; ++i)
				{
					const double & z = mock_source_map.at("Z_B").at(i);
					if((z<source_z_min)||(z>source_z_max)) continue;

					const double & mag = mock_source_map.at("MAG_r").at(i);
					if((mag<brgastro::mag_m_min)||(mag>brgastro::mag_m_max)) continue;

					#ifdef LENSING_STYLE
					mock_source_positions.push_back(std::tie(
						mock_source_map.at("ra_radians").at(i),
						mock_source_map.at("dec_radians").at(i),
						mock_source_map.at("Z_B").at(i)));
					#else
					mock_source_positions.push_back(std::make_pair(
						mock_source_map.at("ra_radians").at(i),
						mock_source_map.at("dec_radians").at(i)));
					#endif
				}
			}

			// Set up the correlation function estimator
			#ifdef LENSING_STYLE
			brgastro::lensing_correlation_function_estimator estimator(r_limits,lens_positions,
					source_positions,mock_lens_positions,mock_source_positions,z_buffer);
			#else
			brgastro::correlation_function_estimator estimator(r_limits,lens_positions,
					source_positions,mock_lens_positions,mock_source_positions);
			#endif

			// Get the correlation functions

			monopole_corr_func = estimator.calculate();
			dipole_1_corr_func = estimator.calculate_dipole(0);
			dipole_2_corr_func = estimator.calculate_dipole(0.5);
			quadrupole_1_corr_func = estimator.calculate_quadrupole(0);
			quadrupole_2_corr_func = estimator.calculate_quadrupole(0.5);
			octopole_1_corr_func = estimator.calculate_octopole(0);
			octopole_2_corr_func = estimator.calculate_octopole(0.5);
		}
		catch (const std::exception &e)
		{

			#ifdef _OPENMP
			#pragma omp critical(CFHTLenS_corr_func_combine_fields)
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
				monopole_corr_func_sum += monopole_corr_func;
				dipole_1_corr_func_sum += dipole_1_corr_func;
				dipole_2_corr_func_sum += dipole_2_corr_func;
				quadrupole_1_corr_func_sum += quadrupole_1_corr_func;
				quadrupole_2_corr_func_sum += quadrupole_2_corr_func;
				octopole_1_corr_func_sum += octopole_1_corr_func;
				octopole_2_corr_func_sum += octopole_2_corr_func;
				std::cout << "Field " << field_name_root << " (#" <<
						++num_processed << "/" << num_fields << ") complete!\n";
			}
			catch (const std::exception &e)
			{
				std::cerr << "Error combining data from field " << field_name_root << " (#" <<
						++num_processed << "/" << num_fields << ")!\n"
						<< e.what();
			}
		}

	}

	// Set up the output table
	brgastro::table_map_t<double> output_table;
	output_table["r_bin_mid_radians"] = r_limits.get_bin_mids();

	// Divide the sums by the number of fields to get the means
	output_table["mp_eps"] = brgastro::coerce<std::vector<double>>(static_cast< Eigen::ArrayXd >(
			monopole_corr_func_sum/static_cast<double>(num_fields)));
	output_table["dp_1_eps"] = brgastro::coerce<std::vector<double>>(static_cast< Eigen::ArrayXd >(
			dipole_1_corr_func_sum/static_cast<double>(num_fields)));
	output_table["dp_2_eps"] = brgastro::coerce<std::vector<double>>(static_cast< Eigen::ArrayXd >(
			dipole_2_corr_func_sum/static_cast<double>(num_fields)));
	output_table["qp_1_eps"] = brgastro::coerce<std::vector<double>>(static_cast< Eigen::ArrayXd >(
			quadrupole_1_corr_func_sum/static_cast<double>(num_fields)));
	output_table["qp_2_eps"] = brgastro::coerce<std::vector<double>>(static_cast< Eigen::ArrayXd >(
			quadrupole_2_corr_func_sum/static_cast<double>(num_fields)));
	output_table["op_1_eps"] = brgastro::coerce<std::vector<double>>(static_cast< Eigen::ArrayXd >(
			octopole_1_corr_func_sum/static_cast<double>(num_fields)));
	output_table["op_2_eps"] = brgastro::coerce<std::vector<double>>(static_cast< Eigen::ArrayXd >(
			octopole_2_corr_func_sum/static_cast<double>(num_fields)));

	brgastro::print_table_map(output_name,output_table);


	return 0;
}

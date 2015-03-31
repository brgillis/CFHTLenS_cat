/**********************************************************************\
 @file main_lensing_model_fitter.cpp
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

#include <cmath>
#include <iostream>
#include <string>
#include <utility>

#include <boost/bimap.hpp>

#include "brg/container/labeled_array.hpp"
#include "brg/math/solvers/solvers.hpp"

#include "brg_lensing/lensing_tNFW_profile.h"

#include "brg_physics/units/unit_conversions.hpp"
#include "brg_physics/units/unitconv_map.hpp"

#include "common.h"

#include "get_errors_on_fit.hpp"
#include "lensing_fitting_bin.hpp"

// Magic values
constexpr flt_type init_frac_in_groups = 0.1;
constexpr flt_type min_frac_in_groups = 0.;
constexpr flt_type max_frac_in_groups = 1.;

constexpr flt_type max_group_mass = 1e16*brgastro::unitconv::Msuntokg;

constexpr flt_type max_Sigma_offset = 50 * brgastro::unitconv::Msuntokg/(brgastro::unitconv::pctom*brgastro::unitconv::pctom);
constexpr flt_type min_Sigma_offset = -max_Sigma_offset;

constexpr int_type MCMC_max_steps = 10000;
constexpr int_type MCMC_annealing_period = MCMC_max_steps/10;

/**
 * TODO (description)
 *
 * @param TODO (params to be passed at command-line)
 * @return
 */
int main( const int argc, const char *argv[] )
{
	// Check that we got a filename in the command line
	if(argc<=2)
	{
		std::cerr << "ERROR: Filename of lensing data to be fit must be passed at command line.\n";
		return 1;
	}

	// Load in the lensing data
	const std::string input_filename = argv[1];
	brgastro::labeled_array<flt_type> lensing_data(input_filename);

	std::string fitting_results_filename = argv[2];

	bool output_new_models = false;
	std::string output_filename;
	if(argc>=4)
	{
		output_new_models = true;
		output_filename = argv[3];
	}

	// Convert the units of the array to standard set
	brgastro::unitconv_map u_map;

	u_map[R_min_label] = u_map[R_max_label] = u_map[shear_R_mean_label] = u_map[magf_R_mean_label] =
		brgastro::unitconv::kpctom;
	u_map[m_min_label] = u_map[m_max_label] = u_map[shear_lens_m_mean_label] = u_map[magf_lens_m_mean_label] =
		brgastro::unitconv::Msuntokg;
	u_map[dS_t_mean_label] = u_map[dS_t_std_label] = u_map[dS_t_stderr_label] =
		u_map[dS_x_mean_label] = u_map[dS_x_std_label] = u_map[dS_x_stderr_label] =
		u_map[model_dS_t_label] = u_map[Sigma_mean_label] = u_map[Sigma_stderr_label] =
		u_map[model_Sigma_label] = u_map[shear_Sigma_crit_label] = u_map[magf_Sigma_crit_label] =
			brgastro::unitconv::Msuntokg/brgastro::square(brgastro::unitconv::pctom);
	u_map[magf_area_label] = brgastro::square(brgastro::unitconv::asectorad);

	lensing_data.apply_inverse_unitconvs(u_map);

	std::cout << "Loaded data from " << input_filename << ".\n";

	typedef std::pair<flt_type,flt_type> bound_type;

	typedef boost::bimap< bound_type, int_type > map_type;

	map_type z_bounds;
	map_type m_bounds;

	// Fill up the bound sets
	for( const auto & bin : lensing_data.rows())
	{
		z_bounds.left.insert(map_type::left_value_type(bound_type(bin.at_label(z_min_label),bin.at_label(z_max_label)),z_bounds.size()));
		m_bounds.left.insert(map_type::left_value_type(bound_type(bin.at_label(m_min_label),bin.at_label(m_max_label)),m_bounds.size()));
	}

	int num_z_bins = z_bounds.size();
	int num_m_bins = m_bounds.size();

	// Output the number of z and m bins we've found
	std::cout << "Found " << num_z_bins << " redshift bins and " << num_m_bins << " mass bins, for "
		<< "a total of " << num_z_bins*num_m_bins << " bins.\n";

	std::vector<std::vector<lensing_fitting_bin>> fitting_bins(num_z_bins,std::vector<lensing_fitting_bin>(num_m_bins));

	// Set the proper bounds for each bin
	for(int z_i=0; z_i<num_z_bins; ++z_i)
	{
		const auto & zb = z_bounds.right.at(z_i);
		for(int m_i=0; m_i<num_m_bins; ++m_i)
		{
			const auto & mb = m_bounds.right.at(m_i);
			fitting_bins.at(z_i).at(m_i).set_bounds(zb.first,zb.second,mb.first,mb.second);
		}
	}

	// Go through each bin again and add it to the proper set of bins
	for( const auto & row : lensing_data.rows())
	{
		const bound_type zb(row.at_label(z_min_label),row.at_label(z_max_label));
		const bound_type mb(row.at_label(m_min_label),row.at_label(m_max_label));

		const int_type & z_i = z_bounds.left.at(zb);
		const int_type & m_i = m_bounds.left.at(mb);

		fitting_bins.at(z_i).at(m_i).insert(row);
	}

	// Set up the labels for the results labeled_array
	brgastro::labeled_array<flt_type> results;

	std::vector<std::string> result_header;
	result_header.push_back("z_min");
	result_header.push_back("m_min");
	result_header.push_back("shear_Chi_squared");
	result_header.push_back("shear_sat_m_best");
	result_header.push_back("shear_sat_m_err");
	result_header.push_back("shear_group_m_best");
	result_header.push_back("shear_group_m_err");
	result_header.push_back("shear_sat_frac_best");
	result_header.push_back("shear_sat_frac_err");
	result_header.push_back("magf_Chi_squared");
	result_header.push_back("magf_sat_m_best");
	result_header.push_back("magf_sat_m_err");
	result_header.push_back("magf_group_m_best");
	result_header.push_back("magf_group_m_err");
	result_header.push_back("magf_sat_frac_best");
	result_header.push_back("magf_sat_frac_err");
	result_header.push_back("magf_Sigma_offset_best");
	result_header.push_back("magf_Sigma_offset_err");
	result_header.push_back("overall_Chi_squared");
	result_header.push_back("overall_sat_m_best");
	result_header.push_back("overall_sat_m_err");
	result_header.push_back("overall_group_m_best");
	result_header.push_back("overall_group_m_err");
	result_header.push_back("overall_sat_frac_best");
	result_header.push_back("overall_sat_frac_err");
	result_header.push_back("overall_Sigma_offset_best");
	result_header.push_back("overall_Sigma_offset_err");

	results.set_labels(std::move(result_header));

	// Set up a vector of the best fit results
	std::vector<std::vector<std::vector<std::vector<flt_type>>>> best_fit_params(num_z_bins,
																				  std::vector<std::vector<std::vector<flt_type>>>(num_m_bins));

	// Set up the models we'll use
	const auto dS_t_model = [] (const flt_type & sat_m, const flt_type & z, const flt_type & group_m, const flt_type & sat_frac,
		const flt_type & R) {
		flt_type one_halo_term = brgastro::lensing_tNFW_profile(sat_m,z).WLsig(R);
		flt_type offset_halo_term = sat_frac*brgastro::lensing_tNFW_profile(group_m,z).quick_group_WLsig(R);

		return one_halo_term+offset_halo_term;
	};
	const auto Sigma_model = [] (const flt_type & sat_m, const flt_type & z, const flt_type & group_m, const flt_type & sat_frac,
		const flt_type & Sigma_offset, const flt_type & R) {
		flt_type one_halo_term = brgastro::lensing_tNFW_profile(sat_m,z).Sigma(R);
		flt_type offset_halo_term = sat_frac*brgastro::lensing_tNFW_profile(group_m,z).quick_group_Sigma(R);

		return one_halo_term+offset_halo_term+Sigma_offset;
	};

	// Go through each fitting bin and fit the best model to it
	for( const auto & fitting_bin_row : fitting_bins )
	{
		for( const auto & fitting_bin : fitting_bin_row )
		{
			// Set up fitting functions

			// Shear fitting function
			const auto shear_fitting_function = [&dS_t_model, &fitting_bin] (const flt_type & sat_m, const flt_type & group_m, const flt_type & sat_frac)
			{

				const auto shear_chi_squared_func = [&dS_t_model, &sat_m, &group_m, &sat_frac, &fitting_bin] (const flt_type & R)
				{
					return dS_t_model(sat_m,fitting_bin.z_mid(),group_m,sat_frac,R);
				};
				return fitting_bin.get_shear_chi_squared(shear_chi_squared_func);
			};

			// Vector version wrapper
			const auto shear_fitting_vector_function = [&shear_fitting_function] (const std::vector<flt_type> & in)
			{
				assert(in.size()>=3);
				return std::vector<flt_type>(1,shear_fitting_function(std::pow(10.,in[0]),std::pow(10.,in[1]),in[2]));
			};

			const auto magf_fitting_function = [&Sigma_model, &fitting_bin] (const flt_type & sat_m, const flt_type & group_m, const flt_type & sat_frac, const flt_type & Sigma_offset)
			{
				const auto magf_chi_squared_func = [&Sigma_model, &sat_m, &group_m, &sat_frac, &Sigma_offset, &fitting_bin] (const flt_type & R)
				{
					return Sigma_model(sat_m,fitting_bin.z_mid(),group_m,sat_frac,Sigma_offset,R);
				};

				return fitting_bin.get_magf_chi_squared(magf_chi_squared_func);
			};

			// Vector version wrapper
			const auto magf_fitting_vector_function = [&magf_fitting_function] (const std::vector<flt_type> & in)
			{
				assert(in.size()>=4);
				return std::vector<flt_type>(1,magf_fitting_function(std::pow(10.,in[0]),std::pow(10.,in[1]),in[2],in[3]));
			};

			const auto fitting_function = [&shear_fitting_function, &magf_fitting_function]
										   (const flt_type & sat_m, const flt_type & group_m, const flt_type & sat_frac,
											   const flt_type & Sigma_offset)
			{
				flt_type shear_chi_squared = shear_fitting_function(sat_m,group_m,sat_frac);
				flt_type magf_chi_squared = magf_fitting_function(sat_m,group_m,sat_frac,Sigma_offset);

				return shear_chi_squared+magf_chi_squared;
			};

			const auto fitting_vector_function = [&fitting_function]
												  (const std::vector<flt_type> & in)
			{
				assert(in.size()>=4);
				return std::vector<flt_type>(1,fitting_function(std::pow(10.,in[0]),std::pow(10.,in[1]),in[2],in[3]));
			};

			// Set up boundaries and step values. Parameters are in order: sat_m, group_m, Sigma_offset
			std::vector<flt_type> init_in = { std::log10(50*fitting_bin.m_mid()),
				std::log10(10000*fitting_bin.m_mid()),
				init_frac_in_groups,
				0. };
			std::vector<flt_type> min_in = { std::log10(fitting_bin.m_mid()),
				std::log10(fitting_bin.m_mid()),
				min_frac_in_groups,
				min_Sigma_offset };
			std::vector<flt_type> max_in = { std::log10(200*fitting_bin.m_mid()),
				std::log10(max_group_mass),
				max_frac_in_groups,
				max_Sigma_offset };
			std::vector<flt_type> in_step = {  0.1,
				0.1,
				(max_frac_in_groups-min_frac_in_groups)/100.,
				(max_Sigma_offset-min_Sigma_offset)/100. };

			std::vector<flt_type> best_shear_in_params = brgastro::solve_MCMC(&shear_fitting_vector_function,init_in,
																			  min_in,max_in,in_step,
																			  MCMC_max_steps,MCMC_annealing_period);

			std::vector<flt_type> best_magf_in_params = brgastro::solve_MCMC(&magf_fitting_vector_function,init_in,
																			  min_in,max_in,in_step,
																			  MCMC_max_steps,MCMC_annealing_period);

			std::vector<flt_type> best_overall_in_params = brgastro::solve_MCMC(&fitting_vector_function,init_in,
																			  min_in,max_in,in_step,
																			  MCMC_max_steps,MCMC_annealing_period);

			flt_type best_shear_chi_squared = shear_fitting_vector_function(best_shear_in_params).at(0);
			flt_type best_magf_chi_squared = magf_fitting_vector_function(best_magf_in_params).at(0);
			flt_type best_overall_chi_squared = fitting_vector_function(best_overall_in_params).at(0);

			// Get the errors on each of these fits
			std::vector<flt_type> shear_errors = get_errors_on_in_params(shear_fitting_vector_function,best_shear_in_params,
																		 best_shear_chi_squared);
			std::vector<flt_type> magf_errors = get_errors_on_in_params(magf_fitting_vector_function,best_magf_in_params,
																		 best_magf_chi_squared);
			std::vector<flt_type> overall_errors = get_errors_on_in_params(fitting_vector_function,best_overall_in_params,
																		 best_overall_chi_squared);

			// Set the best-fitting sigma offset for the best_shear_in_params vector to that from the best_overall_in_params vector
			best_shear_in_params.at(3) = best_overall_in_params.at(3);

			std::vector<flt_type> new_row = { fitting_bin.z_min(), fitting_bin.m_min(),
				best_shear_chi_squared,
				best_shear_in_params.at(0), shear_errors.at(0),
				best_shear_in_params.at(1), shear_errors.at(1),
				best_shear_in_params.at(2), shear_errors.at(2),
				best_magf_chi_squared,
				best_magf_in_params.at(0), magf_errors.at(0),
				best_magf_in_params.at(1), magf_errors.at(1),
				best_magf_in_params.at(2), magf_errors.at(2),
				best_magf_in_params.at(3), magf_errors.at(3),
				best_overall_chi_squared,
				best_overall_in_params.at(0), overall_errors.at(0),
				best_overall_in_params.at(1), overall_errors.at(1),
				best_overall_in_params.at(2), overall_errors.at(2),
				best_overall_in_params.at(3), overall_errors.at(3) };
			results.insert_row(std::move(new_row));

			std::cout << "Finished processing bin " << fitting_bin.z_min() << ", " << fitting_bin.m_min() << ".\n";

			// Save this if we plan to output models later

			if(output_new_models)
			{
				// Get the index for this one's z and m bounds
				const bound_type zb(fitting_bin.z_min(),fitting_bin.z_max());
				const bound_type mb(fitting_bin.m_min(),fitting_bin.m_max());

				const int_type & z_i = z_bounds.left.at(zb);
				const int_type & m_i = m_bounds.left.at(mb);

				best_fit_params.at(z_i).at(m_i).reserve(3);
				best_fit_params.at(z_i).at(m_i).push_back(best_shear_in_params);
				best_fit_params.at(z_i).at(m_i).push_back(best_magf_in_params);
				best_fit_params.at(z_i).at(m_i).push_back(best_overall_in_params);
			}
		}
	}

	results.save(std::cout,true);
	results.save(fitting_results_filename,true);

	// Output the models if desired
	if(output_new_models)
	{
		// Set up the desired columns first
		lensing_data[best_fit_shear_model_dS_t_label];
		lensing_data[best_fit_shear_model_Sigma_label];
		lensing_data[best_fit_magf_model_dS_t_label];
		lensing_data[best_fit_magf_model_Sigma_label];
		lensing_data[best_fit_overall_model_dS_t_label];
		lensing_data[best_fit_overall_model_Sigma_label];

		// Loop through rows of the table now
		for(auto row : lensing_data.rows())
		{
			// Get the index for this one's z and m bounds
			const bound_type zb(row.at_label(z_min_label),row.at_label(z_max_label));
			const bound_type mb(row.at_label(m_min_label),row.at_label(m_max_label));

			const int_type & z_i = z_bounds.left.at(zb);
			const int_type & m_i = m_bounds.left.at(mb);

			const std::vector<flt_type> & best_shear_fit_params = best_fit_params.at(z_i).at(m_i).at(0);
			const std::vector<flt_type> & best_magf_fit_params = best_fit_params.at(z_i).at(m_i).at(1);
			const std::vector<flt_type> & best_overall_fit_params = best_fit_params.at(z_i).at(m_i).at(2);

			row.at_label(best_fit_shear_model_dS_t_label) = dS_t_model(std::pow(10.,best_shear_fit_params.at(0)),
																	   (zb.first+zb.second)/2.,
																	   std::pow(10.,best_shear_fit_params.at(1)),
																	   best_shear_fit_params.at(2),
																	   row.at_label(shear_R_mean_label));
			row.at_label(best_fit_shear_model_Sigma_label) = Sigma_model(std::pow(10.,best_shear_fit_params.at(0)),
																		 (zb.first+zb.second)/2.,
																		 std::pow(10.,best_shear_fit_params.at(1)),
																		 best_shear_fit_params.at(2),
																	     best_shear_fit_params.at(3),
																		 row.at_label(magf_R_mean_label));

			row.at_label(best_fit_magf_model_dS_t_label) = dS_t_model(std::pow(10.,best_magf_fit_params.at(0)),
																	   (zb.first+zb.second)/2.,
																	   std::pow(10.,best_magf_fit_params.at(1)),
																	   best_magf_fit_params.at(2),
																	   row.at_label(shear_R_mean_label));
			row.at_label(best_fit_magf_model_Sigma_label) = Sigma_model(std::pow(10.,best_magf_fit_params.at(0)),
																		 (zb.first+zb.second)/2.,
																		 std::pow(10.,best_magf_fit_params.at(1)),
																		 best_magf_fit_params.at(2),
																	     best_magf_fit_params.at(3),
																		 row.at_label(magf_R_mean_label));

			row.at_label(best_fit_overall_model_dS_t_label) = dS_t_model(std::pow(10.,best_overall_fit_params.at(0)),
																	   (zb.first+zb.second)/2.,
																	   std::pow(10.,best_overall_fit_params.at(1)),
																	   best_overall_fit_params.at(2),
																	   row.at_label(shear_R_mean_label));
			row.at_label(best_fit_overall_model_Sigma_label) = Sigma_model(std::pow(10.,best_overall_fit_params.at(0)),
																		 (zb.first+zb.second)/2.,
																		 std::pow(10.,best_overall_fit_params.at(1)),
																		 best_overall_fit_params.at(2),
																	     best_overall_fit_params.at(3),
																		 row.at_label(magf_R_mean_label));

		}

		u_map[best_fit_shear_model_dS_t_label] = u_map[best_fit_shear_model_Sigma_label] =
			u_map[best_fit_magf_model_dS_t_label] = u_map[best_fit_magf_model_Sigma_label] =
				u_map[best_fit_overall_model_dS_t_label] = u_map[best_fit_overall_model_Sigma_label] = u_map[model_Sigma_label];

		lensing_data.apply_unitconvs(u_map);
		lensing_data.save(output_filename);
	}

	return 0;
}

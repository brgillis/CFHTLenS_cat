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
#include <Eigen/LU>
#include <iostream>
#include <string>
#include <sstream>
#include <utility>

#include <boost/algorithm/string/replace.hpp>
#include <boost/bimap.hpp>

#include "IceBRG_main/container/labeled_array.hpp"
#include "IceBRG_main/file_system.hpp"
#include "IceBRG_main/math/solvers/solver_classes.hpp"
#include "IceBRG_main/units/unit_conversions.hpp"
#include "IceBRG_main/units/unitconv_map.hpp"

#include "IceBRG_lensing/lensing_tNFW_profile.h"
#include "IceBRG_lensing/lensing_tNFW_caches.h"

#include "common.h"

#include "lensing_fitting_bin.hpp"

#define SKIP_CACHE_OUTPUT

using namespace IceBRG;

// Magic values
constexpr flt_type init_frac_in_groups = 0.1;
constexpr flt_type min_frac_in_groups = 0.;
constexpr flt_type max_frac_in_groups = 1.;

const mass_type max_group_mass = 1e16*unitconv::Msuntokg*kg;

constexpr flt_type max_Sigma_offset = 50 * unitconv::Msuntokg/(unitconv::pctom*unitconv::pctom);
constexpr flt_type min_Sigma_offset = -max_Sigma_offset;

constexpr int_type MCMC_max_steps = 10000;
constexpr int_type MCMC_annealing_period = MCMC_max_steps/10;
constexpr flt_type MCMC_annealing_factor = 4;
constexpr int_type MCMC_skip_factor = 10;

#ifdef NDEBUG
constexpr int_type store_points_factor = 10000;
#else
constexpr int_type store_points_factor = 1;
#endif

constexpr int_type num_params = 4;

#ifndef SKIP_CACHE_OUTPUT
constexpr const char * offset_Sigma_cache_name = "tNFW_offset_Sigma_cache.dat";
constexpr const char * group_Sigma_cache_name = "tNFW_group_Sigma_cache.dat";
constexpr const char * offset_sig_cache_name = "tNFW_offset_Delta_Sigma_cache.dat";
constexpr const char * group_sig_cache_name = "tNFW_group_Delta_Sigma_cache.dat";
#endif

typedef Eigen::Array<flt_type,num_params,1> array_type;
typedef MCMC_solver<flt_type,num_params> solver_type;

/**
 * TODO (description)
 *
 * @param TODO (params to be passed at command-line)
 * @return
 */
int_type main( const int_type argc, const char *argv[] )
{

	// Check that we got a filename in the command line
	if(argc<=2)
	{
		std::cerr << "ERROR: Filename of lensing data to be fit and output name must be passed at "
				<< "command line.\n";
		return 1;
	}

	// Load in the lensing data
	const std::string input_filename = argv[1];
	labeled_array<flt_type> lensing_data(input_filename);

	std::string fitting_results_filename = argv[2];

	bool output_new_models = false;
	std::string output_filename;
	if(argc>=4)
	{
		output_new_models = true;
		output_filename = argv[3];
	}

	flt_type fitting_R_min = 0;
	if(argc>=5)
	{
		fitting_R_min = boost::lexical_cast<flt_type>(argv[4])*unitconv::kpctom;
	}

	// Check if there's any request to store test points for any bins
	std::vector<std::pair<flt_type,flt_type>> bins_to_store_points_for;
	if(argc>5)
	{
		// Make sure there's an even number of extra arguments
		int_type num_extra_args = argc-5;
		if(num_extra_args % 2 != 0)
		{
			std::cerr << "ERROR: Bins for which test point data should be stored must be passed in "
					<< "pairs of z, m.";
			return 1;
		}
		int_type num_pairs = num_extra_args/2;

		// For each pair, store the point in the bins_to_store_points_for vector
		for(int i = 0; i < num_pairs; ++i)
		{
			flt_type z_point = std::atof(argv[5+2*i]);
			flt_type m_point = std::atof(argv[5+2*i+1]); // In log10(Msun)
			m_point = std::pow(10.,m_point)*unitconv::Msuntokg;
			bins_to_store_points_for.push_back( std::make_pair(z_point,m_point) );
		}
	}

	// Output caches
	#ifndef SKIP_CACHE_OUTPUT
	tNFW_offset_Sigma_cache().print(offset_Sigma_cache_name);
	tNFW_group_Sigma_cache().print(group_Sigma_cache_name);
	tNFW_offset_sig_cache().print(offset_sig_cache_name);
	tNFW_group_sig_cache().print(group_sig_cache_name);
	#endif

	// Convert the units of the array to standard set
	unitconv_map u_map;

	u_map[R_min_label] = u_map[R_max_label] = u_map[shear_R_mean_label] = u_map[magf_R_mean_label] =
		unitconv::kpctom;

	u_map[m_min_label] = u_map[m_max_label] = u_map[shear_lens_m_mean_label] =
			u_map[magf_lens_m_mean_label] = unitconv::Msuntokg;

	u_map[dS_t_mean_label] = u_map[dS_t_std_label] = u_map[dS_t_stderr_label] =
		u_map[dS_x_mean_label] = u_map[dS_x_std_label] = u_map[dS_x_stderr_label] =
		u_map[model_dS_t_label] = u_map[Sigma_mean_label] = u_map[Sigma_stderr_label] =
		u_map[model_Sigma_label] = u_map[shear_Sigma_crit_label] = u_map[magf_Sigma_crit_label] =
			unitconv::Msuntokg/square(unitconv::pctom);

	u_map[magf_area_label] = square(unitconv::asectorad);

	lensing_data.apply_inverse_unitconvs(u_map);

	std::cout << "Loaded data from " << input_filename << ".\n";

	typedef std::pair<flt_type,flt_type> bound_type;

	typedef boost::bimap< bound_type, int_type > map_type;

	map_type z_bounds;
	map_type m_bounds;

	// Fill up the bound sets
	for( const auto & bin : lensing_data.rows())
	{
		z_bounds.left.insert(map_type::left_value_type(bound_type(bin.at_label(z_min_label),
				bin.at_label(z_max_label)),z_bounds.size()));
		m_bounds.left.insert(map_type::left_value_type(bound_type(bin.at_label(m_min_label),
				bin.at_label(m_max_label)),m_bounds.size()));
	}

	int_type num_z_bins = z_bounds.size();
	int_type num_m_bins = m_bounds.size();

	// Output the number of z and m bins we've found
	std::cout << "Found " << num_z_bins << " redshift bins and " << num_m_bins << " mass bins, for "
		<< "a total of " << num_z_bins*num_m_bins << " bins.\n";

	std::vector<std::vector<lensing_fitting_bin>> fitting_bins(
			num_z_bins,std::vector<lensing_fitting_bin>(num_m_bins));

	// Set the proper bounds for each bin
	for(int_type z_i=0; z_i<num_z_bins; ++z_i)
	{
		const auto & zb = z_bounds.right.at(z_i);
		for(int_type m_i=0; m_i<num_m_bins; ++m_i)
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

		// Check that the row meets the minimum R requirement, and if so, add it to the bin
		if(row.at_label(R_min_label)>fitting_R_min)
			fitting_bins.at(z_i).at(m_i).insert(row);
	}

	// Set up the labels for the results labeled_array
	labeled_array<flt_type> results;

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
	result_header.push_back("shear_sqrt_det_covar");
	result_header.push_back("magf_Chi_squared");
	result_header.push_back("magf_sat_m_best");
	result_header.push_back("magf_sat_m_err");
	result_header.push_back("magf_group_m_best");
	result_header.push_back("magf_group_m_err");
	result_header.push_back("magf_sat_frac_best");
	result_header.push_back("magf_sat_frac_err");
	result_header.push_back("magf_Sigma_offset_best");
	result_header.push_back("magf_Sigma_offset_err");
	result_header.push_back("magf_sqrt_det_covar");
	result_header.push_back("overall_Chi_squared");
	result_header.push_back("overall_sat_m_best");
	result_header.push_back("overall_sat_m_err");
	result_header.push_back("overall_group_m_best");
	result_header.push_back("overall_group_m_err");
	result_header.push_back("overall_sat_frac_best");
	result_header.push_back("overall_sat_frac_err");
	result_header.push_back("overall_Sigma_offset_best");
	result_header.push_back("overall_Sigma_offset_err");
	result_header.push_back("overall_sqrt_det_covar");
	result_header.push_back("Sigma_crit");

	results.set_labels(std::move(result_header));

	// Set up a vector of the best fit results
	std::vector<std::vector<std::vector<std::vector<flt_type>>>> best_fit_params(num_z_bins,
		std::vector<std::vector<std::vector<flt_type>>>(num_m_bins));

	// Set up the models we'll use
	const auto dS_t_model = [] (const flt_type & sat_m, const flt_type & z, const flt_type & group_m,
			const flt_type & sat_frac, const flt_type & R) {
		flt_type one_halo_term = value_of(lensing_tNFW_profile(sat_m*kg,z).Delta_Sigma(R*m));
		flt_type offset_halo_term = value_of(sat_frac *
				lensing_tNFW_profile(group_m*kg,z).quick_group_Delta_Sigma(R*m));

		return one_halo_term+offset_halo_term;
	};
	const auto Sigma_model = [] (const flt_type & sat_m, const flt_type & z, const flt_type & group_m,
			const flt_type & sat_frac,
		const flt_type & Sigma_offset, const flt_type & R) {
		flt_type one_halo_term = value_of(lensing_tNFW_profile(sat_m*kg,z).Sigma(R*m));
		flt_type offset_halo_term = value_of(sat_frac *
				lensing_tNFW_profile(group_m*kg,z).quick_group_Sigma(R*m));

		return one_halo_term+offset_halo_term+Sigma_offset;
	};

	// Go through each fitting bin and fit the best model to it
	for( const auto & fitting_bin_row : fitting_bins )
	{
		for( const auto & fitting_bin : fitting_bin_row )
		{

			// Check if we want to store test points for this bin
			bool store_points = false;
			for( const auto & bin_location : bins_to_store_points_for )
			{
				// Check if this point is inside this bin
				if((fitting_bin.z_min() <= bin_location.first) and
				   (fitting_bin.z_max() > bin_location.first) and
				   (fitting_bin.m_min() <= bin_location.second) and
				   (fitting_bin.m_max() > bin_location.second))
				{
					store_points = true;
					break;
				}
			}

			// Set up fitting functions

			// Shear fitting function
			const auto shear_fitting_function = [&dS_t_model, &fitting_bin]
				(const flt_type & sat_m, const flt_type & group_m, const flt_type & sat_frac)
			{

				const auto shear_chi_squared_func = [&dS_t_model, &sat_m, &group_m, &sat_frac,
													 &fitting_bin] (const flt_type & R)
				{
					return dS_t_model(sat_m,fitting_bin.z_mid(),group_m,sat_frac,R);
				};
				return fitting_bin.get_shear_chi_squared(shear_chi_squared_func);
			};

			// Array version wrapper
			const auto shear_fitting_array_function = [&shear_fitting_function] (const array_type & in)
			{
				return shear_fitting_function(std::pow(10.,in[0]),std::pow(10.,in[1]),in[2]);
			};

			// Magf fitting function
			const auto magf_fitting_function = [&Sigma_model, &fitting_bin]
				(const flt_type & sat_m, const flt_type & group_m, const flt_type & sat_frac,
						const flt_type & Sigma_offset)
			{
				const auto magf_chi_squared_func = [&Sigma_model, &sat_m, &group_m, &sat_frac,
													&Sigma_offset, &fitting_bin] (const flt_type & R)
				{
					return Sigma_model(sat_m,fitting_bin.z_mid(),group_m,sat_frac,Sigma_offset,R);
				};

				return fitting_bin.get_magf_chi_squared(magf_chi_squared_func);
			};

			// Array version wrapper
			const auto magf_fitting_array_function = [&magf_fitting_function] (const array_type & in)
			{
				return magf_fitting_function(std::pow(10.,in[0]),std::pow(10.,in[1]),in[2],in[3]);
			};

			// Overall fitting function
			const auto fitting_function = [&shear_fitting_function, &magf_fitting_function]
										   (const flt_type & sat_m, const flt_type & group_m,
											const flt_type & sat_frac, const flt_type & Sigma_offset)
			{
				flt_type shear_chi_squared = shear_fitting_function(sat_m,group_m,sat_frac);
				flt_type magf_chi_squared = magf_fitting_function(sat_m,group_m,sat_frac,Sigma_offset);

				return shear_chi_squared+magf_chi_squared;
			};

			// Array version wrapper
			const auto fitting_array_function = [&fitting_function]
												  (const array_type & in)
			{
				assert(in.size()>=4);
				return fitting_function(std::pow(10.,in[0]),std::pow(10.,in[1]),in[2],in[3]);
			};

			// Set up boundaries and step values. Parameters are in order: sat_m, group_m, Sigma_offset
			array_type init_in = { std::log10(50*fitting_bin.m_mid()),
				std::log10(10000*fitting_bin.m_mid()),
				init_frac_in_groups,
				0. };
			array_type min_in = { std::log10(fitting_bin.m_mid()),
				std::log10(fitting_bin.m_mid()),
				min_frac_in_groups,
				min_Sigma_offset };
			array_type max_in = { std::log10(200*fitting_bin.m_mid()),
				std::log10(value_of(max_group_mass)),
				max_frac_in_groups,
				max_Sigma_offset };
			array_type in_step = {  0.1,
				0.1,
				(max_frac_in_groups-min_frac_in_groups)/100.,
				(max_Sigma_offset-min_Sigma_offset)/100. };

			solver_type shear_solver;
			solver_type magf_solver;
			solver_type overall_solver;

			// If we're storing points, use 100x as many points so it can be plotted nicely
			int steps_for_this_bin = MCMC_max_steps;
			if(store_points) steps_for_this_bin *= store_points_factor;

			array_type best_shear_in_params = shear_solver.solve(shear_fitting_array_function,init_in,
																  min_in,max_in,in_step,
																  steps_for_this_bin,steps_for_this_bin/4,
																  MCMC_annealing_factor,MCMC_skip_factor);

			array_type best_magf_in_params = magf_solver.solve(magf_fitting_array_function,init_in,
																  min_in,max_in,in_step,
																  steps_for_this_bin,steps_for_this_bin/4,
																  MCMC_annealing_factor,MCMC_skip_factor);

			array_type best_overall_in_params = overall_solver.solve(fitting_array_function,init_in,
																	  min_in,max_in,in_step,
																	  steps_for_this_bin,steps_for_this_bin/4,
																	  MCMC_annealing_factor,MCMC_skip_factor);

			flt_type best_shear_chi_squared = shear_fitting_array_function(best_shear_in_params);
			flt_type best_magf_chi_squared = magf_fitting_array_function(best_magf_in_params);
			flt_type best_overall_chi_squared = fitting_array_function(best_overall_in_params);

			// Get the errors on each of these fits
			array_type shear_errors = shear_solver.get_stddevs();
			array_type magf_errors = magf_solver.get_stddevs();
			array_type overall_errors = overall_solver.get_stddevs();

			// Get covariances
			flt_type shear_det_covar =
					std::sqrt(shear_solver.get_covars().block<num_params-1,num_params-1>(0,0).determinant());
			flt_type magf_det_covar =
					std::sqrt(magf_solver.get_covars().block<num_params-1,num_params-1>(0,0).determinant());
			flt_type overall_det_covar =
					std::sqrt(overall_solver.get_covars().block<num_params-1,num_params-1>(0,0).determinant());

			// Set the best-fitting sigma offset for the best_shear_in_params vector to that from
			// the best_overall_in_params vector
			best_shear_in_params[3] = best_overall_in_params[3];

			std::vector<flt_type> new_row = { fitting_bin.z_min(), fitting_bin.m_min(),
				best_shear_chi_squared,
				best_shear_in_params[0], shear_errors[0],
				best_shear_in_params[1], shear_errors[1],
				best_shear_in_params[2], shear_errors[2],
				shear_det_covar,
				best_magf_chi_squared,
				best_magf_in_params[0], magf_errors[0],
				best_magf_in_params[1], magf_errors[1],
				best_magf_in_params[2], magf_errors[2],
				best_magf_in_params[3], magf_errors[3],
				magf_det_covar,
				best_overall_chi_squared,
				best_overall_in_params[0], overall_errors[0],
				best_overall_in_params[1], overall_errors[1],
				best_overall_in_params[2], overall_errors[2],
				best_overall_in_params[3], overall_errors[3],
				overall_det_covar,
				fitting_bin.Sigma_crit() };
			results.insert_row(std::move(new_row));

			std::cout << "Finished processing bin " << fitting_bin.z_min() << ", "
					<< fitting_bin.m_min() << ".\n";

			// Save this if we plan to output models later

			if(output_new_models)
			{
				// Get the index for this one's z and m bounds
				const bound_type zb(fitting_bin.z_min(),fitting_bin.z_max());
				const bound_type mb(fitting_bin.m_min(),fitting_bin.m_max());

				const int_type & z_i = z_bounds.left.at(zb);
				const int_type & m_i = m_bounds.left.at(mb);

				best_fit_params.at(z_i).at(m_i).reserve(3);
				best_fit_params.at(z_i).at(m_i).push_back(coerce<std::vector<flt_type>>(
						best_shear_in_params));
				best_fit_params.at(z_i).at(m_i).push_back(coerce<std::vector<flt_type>>(
						best_magf_in_params));
				best_fit_params.at(z_i).at(m_i).push_back(coerce<std::vector<flt_type>>(
						best_overall_in_params));
			}

			// Store the generated points if desired
			if(store_points)
			{
				// Get the points from the solver
				const auto & shear_test_points = shear_solver.get_test_run_points().transpose();
				const auto & magf_test_points = magf_solver.get_test_run_points().transpose();
				const auto & overall_test_points = overall_solver.get_test_run_points().transpose();

				// Set up the header for the test points output files
				std::vector<std::string> test_point_header( { "sat_m",
				                                              "group_m",
														      "sat_frac",
														      "kappa_offset"} );

				// Make tables for each of the test points arrays
				labeled_array<flt_type> shear_test_points_table(shear_test_points,test_point_header);
				labeled_array<flt_type> magf_test_points_table(magf_test_points,test_point_header);
				labeled_array<flt_type> overall_test_points_table(overall_test_points,test_point_header);

				// Determine the filenames to use for saving these test points
				std::stringstream ss("");

				ss << std::setprecision(3) << "z-" << fitting_bin.z_min();
				std::string bin_z_label = ss.str();
				boost::replace_all(bin_z_label,".","p"); // Replace periods with p to simplify filenames

				ss.str("");
				ss << std::setprecision(3) << "m-" << std::log10(fitting_bin.m_min()*unitconv::kgtoMsun);
				std::string bin_m_label = ss.str();
				boost::replace_all(bin_m_label,".","p"); // Replace periods with p to simplify filenames

				std::string output_filename_root = get_base_filename(fitting_results_filename);

				std::string output_test_points_filename_base = output_filename_root + "_" + bin_z_label +
						"_" + bin_m_label;

				// Convert the units
				shear_test_points_table.at_label("sat_m").raw() -= std::log10(unitconv::Msuntokg);
				shear_test_points_table.at_label("group_m").raw() -= std::log10(unitconv::Msuntokg);
				shear_test_points_table.at_label("kappa_offset").raw() /= fitting_bin.Sigma_crit();

				magf_test_points_table.at_label("sat_m").raw() -= std::log10(unitconv::Msuntokg);
				magf_test_points_table.at_label("group_m").raw() -= std::log10(unitconv::Msuntokg);
				magf_test_points_table.at_label("kappa_offset").raw() /= fitting_bin.Sigma_crit();

				overall_test_points_table.at_label("sat_m").raw() -= std::log10(unitconv::Msuntokg);
				overall_test_points_table.at_label("group_m").raw() -= std::log10(unitconv::Msuntokg);
				overall_test_points_table.at_label("kappa_offset").raw() /= fitting_bin.Sigma_crit();


				// Save the tables
				shear_test_points_table.save(output_test_points_filename_base + "_shear_test_points.dat");
				magf_test_points_table.save(output_test_points_filename_base + "_magf_test_points.dat");
				overall_test_points_table.save(output_test_points_filename_base + "_overall_test_points.dat");
			}
		}
	}

	results.save(std::cout,true);
	results.save(fitting_results_filename,true);

	std::cout << "Results saved to " << fitting_results_filename << "." << std::endl;

	// Output the models if desired
	if(output_new_models)
	{
		// Set up the desired columns first
		lensing_data[best_fit_shear_model_dS_t_label];
		lensing_data[best_fit_shear_model_Sigma_label];
		lensing_data[best_fit_magf_model_dS_t_label];
		lensing_data[best_fit_magf_model_Sigma_label];
		lensing_data[best_fit_magf_Sigma_offset_label];
		lensing_data[best_fit_overall_model_dS_t_label];
		lensing_data[best_fit_overall_model_Sigma_label];
		lensing_data[best_fit_overall_Sigma_offset_label];

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

			row.at_label(best_fit_shear_model_dS_t_label) =
					dS_t_model(std::pow(10.,best_shear_fit_params.at(0)),
							   (zb.first+zb.second)/2.,
							   std::pow(10.,best_shear_fit_params.at(1)),
							   best_shear_fit_params.at(2),
							   row.at_label(shear_R_mean_label));
			row.at_label(best_fit_shear_model_Sigma_label) =
					Sigma_model(std::pow(10.,best_shear_fit_params.at(0)),
								 (zb.first+zb.second)/2.,
								 std::pow(10.,best_shear_fit_params.at(1)),
								 best_shear_fit_params.at(2),
								 best_shear_fit_params.at(3),
								 row.at_label(magf_R_mean_label));

			row.at_label(best_fit_magf_model_dS_t_label) =
					dS_t_model(std::pow(10.,best_magf_fit_params.at(0)),
							   (zb.first+zb.second)/2.,
							   std::pow(10.,best_magf_fit_params.at(1)),
							   best_magf_fit_params.at(2),
							   row.at_label(shear_R_mean_label));
			row.at_label(best_fit_magf_model_Sigma_label) =
					Sigma_model(std::pow(10.,best_magf_fit_params.at(0)),
								 (zb.first+zb.second)/2.,
								 std::pow(10.,best_magf_fit_params.at(1)),
								 best_magf_fit_params.at(2),
								 best_magf_fit_params.at(3),
								 row.at_label(magf_R_mean_label));
			row.at_label(best_fit_magf_Sigma_offset_label) = best_magf_fit_params.at(3);

			row.at_label(best_fit_overall_model_dS_t_label) =
					dS_t_model(std::pow(10.,best_overall_fit_params.at(0)),
							   (zb.first+zb.second)/2.,
							   std::pow(10.,best_overall_fit_params.at(1)),
							   best_overall_fit_params.at(2),
							   row.at_label(shear_R_mean_label));
			row.at_label(best_fit_overall_model_Sigma_label) =
					Sigma_model(std::pow(10.,best_overall_fit_params.at(0)),
								 (zb.first+zb.second)/2.,
								 std::pow(10.,best_overall_fit_params.at(1)),
								 best_overall_fit_params.at(2),
								 best_overall_fit_params.at(3),
								 row.at_label(magf_R_mean_label));
			row.at_label(best_fit_overall_Sigma_offset_label) = best_overall_fit_params.at(3);

		}

		u_map[best_fit_shear_model_dS_t_label] = u_map[best_fit_shear_model_Sigma_label] =
			u_map[best_fit_magf_model_dS_t_label] = u_map[best_fit_magf_model_Sigma_label] =
				u_map[best_fit_overall_model_dS_t_label] = u_map[best_fit_overall_model_Sigma_label] =
						u_map[model_Sigma_label];

		lensing_data.apply_unitconvs(u_map);
		lensing_data.save(output_filename);

		std::cout << "Data with best-fit models saved to " << output_filename << "." << std::endl;
	}

	return 0;
}

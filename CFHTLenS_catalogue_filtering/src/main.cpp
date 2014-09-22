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

#include "brg/file_access/ascii_table_map.hpp"
#include "brg/file_access/open_file.hpp"
#include "brg/physics/units/unit_conversions.hpp"

#include "get_filtered_indices.h"
#include "make_output_map.h"

// Magic values
std::string fields_directory = "/disk2/brg/git/CFHTLenS_cat/Data/";
std::string fields_list = fields_directory + "fields_list.txt";

int main( const int argc, const char *argv[] )
{

	// Open and read in the fields list
	std::ifstream fi;
	brgastro::open_file_input(fi,fields_list);

	std::string field_name;

	while(fi>>field_name)
	{
		std::string field_name_root = field_name.substr(0,6);

		// Get the input file name
		std::stringstream ss("");
		ss << fields_directory << "full_tables/" << field_name << ".dat";
		std::string input_file_name = ss.str();

		// Get the lens and source output file names
		ss.str("");
		ss << fields_directory << "filtered_tables/" << field_name_root << "_lens.dat";
		std::string lens_output_name = ss.str();

		ss.str("");
		ss << fields_directory << "filtered_tables/" << field_name_root << "_source.dat";
		std::string source_output_name = ss.str();

		// Load in the input file
		brgastro::table_map_t<std::string> table_map;
		try
		{
			table_map = brgastro::load_table_map<std::string>(input_file_name,false,std::string(""));
		}
		catch(const std::runtime_error &e)
		{
			std::cerr << "WARNING: " << e.what() << std::endl;
			continue;
		}

		// Rename the y magnitude column i if it's present
		move_y_column_to_i(table_map);

		// Lens file
		std::vector<size_t> filtered_indices(get_filtered_lenses(table_map));

		std::cout << "Generating " << lens_output_name << "... ";
		std::cout.flush();

		// Set up the header columns vector for the ones we want to output
		brgastro::header_t lens_header_columns;
		lens_header_columns.push_back("SeqNr");
		lens_header_columns.push_back("ALPHA_J2000");
		lens_header_columns.push_back("DELTA_J2000");
		lens_header_columns.push_back("Z_B");
		lens_header_columns.push_back("T_B");
		lens_header_columns.push_back("ODDS");
		lens_header_columns.push_back("LP_log10_SM_MED");
		lens_header_columns.push_back("LP_log10_SM_INF");
		lens_header_columns.push_back("LP_log10_SM_SUP");
		lens_header_columns.push_back("MAG_i");
		lens_header_columns.push_back("MAGERR_i");
		lens_header_columns.push_back("EXTINCTION_i");
		lens_header_columns.push_back("MAG_r");
		lens_header_columns.push_back("MAGERR_r");
		lens_header_columns.push_back("EXTINCTION_r");

		// Set up a map of the conversions to apply
		std::map<std::string,std::function<double(double)>> lens_conversions;
		auto deg_to_rad = [] (double theta) {return theta*brgastro::unitconv::degtorad;};
		auto l10_Msun_to_kg = [] (double l10_Msun)
				{return std::pow(10.,l10_Msun)*brgastro::unitconv::Msuntokg;};
		lens_conversions["ALPHA_J2000"] = deg_to_rad;
		lens_conversions["DELTA_J2000"] = deg_to_rad;
		lens_conversions["LP_log10_SM_MED"] = l10_Msun_to_kg;
		lens_conversions["LP_log10_SM_INF"] = l10_Msun_to_kg;
		lens_conversions["LP_log10_SM_SUP"] = l10_Msun_to_kg;

		brgastro::table_map_t<std::string> output_map(make_output_map(table_map,filtered_indices,lens_header_columns,
				lens_conversions));

		brgastro::print_table_map(lens_output_name,output_map);

		std::cout << "Done!\nGenerating " << source_output_name << "... ";
		std::cout.flush();

		// Source file

		filtered_indices = get_filtered_sources(table_map);

		// Set up the header columns vector for the ones we want to output
		brgastro::header_t source_header_columns;
		source_header_columns.push_back("SeqNr");
		source_header_columns.push_back("ALPHA_J2000");
		source_header_columns.push_back("DELTA_J2000");
		source_header_columns.push_back("Z_B");
		source_header_columns.push_back("T_B");
		source_header_columns.push_back("e1");
		source_header_columns.push_back("e2");
		source_header_columns.push_back("weight");
		source_header_columns.push_back("m");
		source_header_columns.push_back("c2");
		source_header_columns.push_back("LP_log10_SM_MED");
		source_header_columns.push_back("LP_log10_SM_INF");
		source_header_columns.push_back("LP_log10_SM_SUP");
		source_header_columns.push_back("MAG_i");
		source_header_columns.push_back("MAGERR_i");
		source_header_columns.push_back("EXTINCTION_i");
		source_header_columns.push_back("MAG_r");
		source_header_columns.push_back("MAGERR_r");
		source_header_columns.push_back("EXTINCTION_r");

		// Set up a map of the conversions to apply
		std::map<std::string,std::function<double(double)>> source_conversions;
		source_conversions["ALPHA_J2000"] = deg_to_rad;
		source_conversions["DELTA_J2000"] = deg_to_rad;
		source_conversions["LP_log10_SM_MED"] = l10_Msun_to_kg;
		source_conversions["LP_log10_SM_INF"] = l10_Msun_to_kg;
		source_conversions["LP_log10_SM_SUP"] = l10_Msun_to_kg;

		output_map = make_output_map(table_map,filtered_indices,source_header_columns,source_conversions);

		brgastro::print_table_map(source_output_name,output_map);

		std::cout << "Done!\n";

	}

	return 0;
}


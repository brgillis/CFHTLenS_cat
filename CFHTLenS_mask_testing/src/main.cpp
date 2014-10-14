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

#include "brg/file_access/open_file.hpp"

#include "load_pixel_table.h"

// Magic values
std::string data_directory = "/disk2/brg/git/CFHTLenS_cat/Data/";
std::string mask_directory = data_directory + "masks/";
std::string fields_list = data_directory + "fields_list.txt";

/*
 * Pixel is valid galaxy position if surrounded by unmasked pixels
 *
 * Max physical separation is 2000 kpc.
 * At redshift 0.2 this is 601.5"
 * Pixel scale is 0.186"/px
 * Maximum pixel separation we care about is thus 3235 px.
 */

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
		ss << mask_directory << field_name << ".fits.fz";
		std::string input_file_name = ss.str();

		// Get the output file names
		ss.str("");
		ss << mask_directory << field_name_root << "_corr.dat";
		std::string lens_output_name = ss.str();

		std::vector<std::vector<bool>> pixel_table = load_pixel_table(input_file_name);

	}

	return 0;
}


/**********************************************************************\
 @file gg_lensing_config.cpp
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
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "brg/file_access/open_file.hpp"
#include "brg/file_access/trim_comments.hpp"
#include "brg/lexical_cast.hpp"

#include "gg_lensing_config.h"

gg_lensing_config::gg_lensing_config( const int argc, const char *argv[] )
{
	using namespace brgastro::unitconv;

	if(argc==1)
	{
		use_precalculated_data = false;
		precalculated_data_filename = "";

		R_min = 10*kpctom;
		R_max = 2000*kpctom;
		R_step = 10*kpctom;
		R_log = true;
		R_bins = 100;

		m_min = 1e8*Msuntokg;
		m_max = 1e11*Msuntokg;
		m_step = 1e11*Msuntokg;
		m_log = true;
		m_bins = 3;

		z_min = 0.2;
		z_max = 1.1;
		z_step = 0.9;
		z_log = true;
		z_bins = 1;

		mag_min = -std::numeric_limits<double>::infinity();
		mag_max = 25;
		mag_step = std::numeric_limits<double>::infinity();
		mag_log = false;
		mag_bins = 1;

		z_buffer = 0.1;
	}
	else
	{
		// Open the config file
		std::string filename(argv[1]);
		std::ifstream fi;
		brgastro::open_file_input(fi,filename);

		// Set up a vector to store config values in
		std::vector<std::string> config_value_strings(num_config_params);

		auto get_config_value = [] (std::istream & fi)
		{
			std::string line_buffer, word_buffer;
			std::istringstream line_data_stream;

			brgastro::trim_comments_all_at_top(fi);

			do
			{
				std::getline( fi, line_buffer );
			} while(line_buffer.size()==0);
			line_data_stream.str(line_buffer);
			do
			{
				line_data_stream >> word_buffer;
			} while(line_data_stream);
			return word_buffer;
		};

		for(std::string & val : config_value_strings)
		{
			val = get_config_value(fi);
		}

		// Load in the values
		size_t i=0;
		use_precalculated_data = brgastro::bool_cast(config_value_strings.at(i++));
		precalculated_data_filename = config_value_strings.at(i++);

		R_min = brgastro::min_cast<double>(config_value_strings.at(i++))*kpctom;
		R_max = brgastro::max_cast<double>(config_value_strings.at(i++))*kpctom;
		R_step = brgastro::max_cast<double>(config_value_strings.at(i++))*kpctom;
		R_log = brgastro::bool_cast(config_value_strings.at(i++));
		R_bins = boost::lexical_cast<size_t>(config_value_strings.at(i++));

		m_min = brgastro::min_cast<double>(config_value_strings.at(i++))*Msuntokg;
		m_max = brgastro::max_cast<double>(config_value_strings.at(i++))*Msuntokg;
		m_step = brgastro::max_cast<double>(config_value_strings.at(i++))*Msuntokg;
		m_log = brgastro::bool_cast(config_value_strings.at(i++));
		m_bins = boost::lexical_cast<size_t>(config_value_strings.at(i++));

		z_min = brgastro::min_cast<double>(config_value_strings.at(i++));
		z_max = brgastro::max_cast<double>(config_value_strings.at(i++));
		z_step = brgastro::max_cast<double>(config_value_strings.at(i++));
		z_log = brgastro::bool_cast(config_value_strings.at(i++));
		z_bins = boost::lexical_cast<size_t>(config_value_strings.at(i++));

		mag_min = brgastro::min_cast<double>(config_value_strings.at(i++));
		mag_max = brgastro::max_cast<double>(config_value_strings.at(i++));
		mag_step = brgastro::max_cast<double>(config_value_strings.at(i++));
		mag_log = brgastro::bool_cast(config_value_strings.at(i++));
		mag_bins = boost::lexical_cast<size_t>(config_value_strings.at(i++));

		z_buffer = boost::lexical_cast<double>(config_value_strings.at(i++));

	}
}

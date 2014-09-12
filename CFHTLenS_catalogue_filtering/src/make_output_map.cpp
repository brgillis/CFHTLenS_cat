/**********************************************************************\
 @file make_output_map.cpp
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


#include <functional>
#include <map>
#include <string>
#include <vector>

#include <boost/lexical_cast.hpp>

#include "brg/file_access/table_typedefs.hpp"

#include "make_output_map.h"

brgastro::table_map_t<std::string> make_output_map(const brgastro::table_map_t<std::string> & map,
		const std::vector<size_t> & filtered_indices, const brgastro::header_t & header_columns,
		const std::map<std::string,std::function<double(double)>> & conversions)
{
	brgastro::table_map_t<std::string> result_map;

	// Loop over columns
	for(auto col_it = header_columns.begin(); col_it != header_columns.end(); ++col_it)
	{
		// Initialize column
		result_map[*col_it] = std::vector<std::string>();

		// Check if we'll apply a conversion to this column
		auto conv_it = conversions.find(*col_it);
		const bool apply_conversion = !(conv_it==conversions.end());

		auto ele_to_skip_next = filtered_indices.begin();

		// Add appropriate elements
		for(size_t i=0; i<map.at(*col_it).size(); ++i)
		{
			if(i==*ele_to_skip_next)
			{
				++ele_to_skip_next;
				if(ele_to_skip_next==filtered_indices.end())
					ele_to_skip_next = filtered_indices.begin();
			}
			else
			{
				if(!apply_conversion)
				{
					result_map[*col_it].push_back(map.at(*col_it)[i]);
				}
				else
				{
					double val = boost::lexical_cast<double>(map.at(*col_it)[i]);
					double new_val = conv_it->second(val);
					result_map[*col_it].push_back(boost::lexical_cast<std::string>(new_val));
				}
			}
		}
	}

	return result_map;
}



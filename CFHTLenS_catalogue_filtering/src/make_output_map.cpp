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



#include <map>
#include <string>
#include <vector>

#include "brg/boost_deep_includes/hold_any.hpp"
#include "brg/file_access/table_typedefs.hpp"

#include "make_output_map.h"

brgastro::table_map_t<std::string> make_output_map(const brgastro::table_map_t<boost::hold_any> & map,
		const std::vector<size_t> & filtered_indices, const brgastro::header_t & header_columns)
{
	brgastro::table_map_t<std::string> result_map;

	// Loop over columns
	for(auto col_it = header_columns.begin(); col_it != header_columns.end(); ++col_it)
	{
		// Initialize column
		result_map[*col_it] = std::vector<std::string>();

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
				result_map[*col_it].push_back(boost::any_cast<std::string>(map.at(*col_it)[i]));
			}
		}
	}

	return result_map;
}



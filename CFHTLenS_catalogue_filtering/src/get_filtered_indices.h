/**********************************************************************\
 @file get_filtered_indices.h
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


// body file: get_filtered_indices.cpp

#ifndef _BRG_GET_FILTERED_INDICES_H_INCLUDED_
#define _BRG_GET_FILTERED_INDICES_H_INCLUDED_

#include <vector>

#include "brg/boost_deep_includes/hold_any.hpp"
#include "brg/file_access/table_typedefs.hpp"

void move_y_column_to_i(brgastro::table_map_t<boost::hold_any> & map);

std::vector<size_t> get_filtered_objects(const brgastro::table_map_t<boost::hold_any> & map);

bool column_passes_lens_filter(const std::string & col_name,boost::hold_any value);
bool column_passes_source_filter(const std::string & col_name,boost::hold_any value);
bool column_passes_global_filter(const std::string & col_name,boost::hold_any value);

inline bool column_passes_filter(const std::string & col_name,boost::hold_any value)
{
	if(!column_passes_global_filter(col_name,value)) return false;
	return (column_passes_lens_filter(col_name,value)&&column_passes_source_filter(col_name,value));
}



#endif // _BRG_GET_FILTERED_INDICES_H_INCLUDED_

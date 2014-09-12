/**********************************************************************\
 @file get_filtered_indices.cpp
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

#include <limits>
#include <map>
#include <utility>
#include <vector>

#include <boost/lexical_cast.hpp>

#include "brg/file_access/table_typedefs.hpp"

#include "get_filtered_indices.h"

void move_y_column_to_i(brgastro::table_map_t<std::string> & map)
{
	// Check if the i columns exist, and if not, move in y columns
	if(map.find("MAG_i")==map.end())
	{
		map["MAG_i"] = std::move(map.at("MAG_y"));
		map.erase("MAG_y");
	}
	if(map.find("MAGERR_i")==map.end())
	{
		map["MAGERR_i"] = std::move(map.at("MAGERR_y"));
		map.erase("MAGERR_y");
	}
	if(map.find("EXTINCTION_i")==map.end())
	{
		map["EXTINCTION_i"] = std::move(map.at("EXTINCTION_y"));
		map.erase("EXTINCTION_y");
	}
}

std::vector<size_t> get_filtered_lenses(const brgastro::table_map_t<std::string> & map)
{
	// Determine the number of rows
	const size_t num_rows = map.begin()->second.size();

	std::vector<size_t> filtered_objects;
	filtered_objects.reserve(num_rows);

	for(size_t i=0; i<num_rows; ++i)
	{
		// Go over each column, and check if it passes the filter
		for(auto col_it = map.begin(); col_it != map.end(); ++col_it)
		{
			const std::string & key = col_it->first;
			if(!column_passes_lens_filter(key,col_it->second.at(i)))
			{
				filtered_objects.push_back(i);
				break;
			}
		}
	}

	return filtered_objects;
}

std::vector<size_t> get_filtered_sources(const brgastro::table_map_t<std::string> & map)
{
	// Determine the number of rows
	const size_t num_rows = map.begin()->second.size();

	std::vector<size_t> filtered_objects;
	filtered_objects.reserve(num_rows);

	for(size_t i=0; i<num_rows; ++i)
	{
		// Go over each column, and check if it passes the filter
		for(auto col_it = map.begin(); col_it != map.end(); ++col_it)
		{
			const std::string & key = col_it->first;
			if(!column_passes_lens_filter(key,col_it->second.at(i)))
			{
				filtered_objects.push_back(i);
				break;
			}
		}
	}

	return filtered_objects;
}

bool column_passes_lens_filter(const std::string & col_name,const std::string & value)
{
	if(!column_passes_global_filter(col_name,value)) return false;
	if(col_name=="Z_B")
	{
		double dval = boost::lexical_cast<double>(value);
		return (dval <= 0.8);
	}
	else
	{
		return true;
	}
}
bool column_passes_source_filter(const std::string & col_name,const std::string & value)
{
	if(!column_passes_global_filter(col_name,value)) return false;
	if(col_name=="MAG_i")
	{
		return boost::lexical_cast<double>(value) <= 24.7;
	}
	else
	{
		return true;
	}
}
bool column_passes_global_filter(const std::string & col_name,const std::string & value)
{
	if(col_name=="Z_B")
	{
		double dval = boost::lexical_cast<double>(value);
		return ((0.2 <= dval) && (dval <= 1.3));
	}
	else if(col_name=="ODSS")
	{
		double dval = boost::lexical_cast<double>(value);
		return (0.8 <= dval);
	}
	else
	{
		return true;
	}
}

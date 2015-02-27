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

#include "brg/container/table_typedefs.hpp"
#include "brg/math/misc_math.hpp"

#include "is_good_position.hpp"

#include "get_filtered_indices.h"

#undef WARN_MASK_MISMATCH

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

std::vector<size_t> get_bad_lenses(const brgastro::table_map_t<std::string> & map, const std::vector<std::vector<bool>> good_pixels)
{
	// Determine the number of rows
	const size_t num_rows = map.begin()->second.size();

	std::vector<size_t> bad_indices;
	bad_indices.reserve(num_rows);

	for(size_t i=0; i<num_rows; ++i)
	{
		// Check if this object is within the mask
		const unsigned xp = brgastro::round_int(boost::lexical_cast<double>(map.at("Xpos").at(i)));
		const unsigned yp = brgastro::round_int(boost::lexical_cast<double>(map.at("Ypos").at(i)));

		if(!is_good_position(xp,yp,good_pixels))
		{
			bad_indices.push_back(i);

#ifdef WARN_MASK_MISMATCH
			// Check that the mask value in the table agrees
			if(boost::lexical_cast<unsigned>(map.at("MASK").at(i))<=1)
			{
				std::cerr << "WARNING: Mask value mismatch for lens index " << i << " at position (" << xp << ", " << yp << ")." << std::endl;
				std::cerr << "Table's mask value is " << map.at("MASK").at(i) << ", but saved value is 'false'.\n";
			}
#endif

			continue;
		}
		else
		{
#ifdef WARN_MASK_MISMATCH
			// Check that the mask value in the table agrees
			if(boost::lexical_cast<unsigned>(map.at("MASK").at(i))>1)
			{
				std::cerr << "WARNING: Mask value mismatch for lens index " << i << " at position (" << xp << ", " << yp << ")." << std::endl;
				std::cerr << "Table's mask value is " << map.at("MASK").at(i) << ", but saved value is 'true'.\n";
			}
#endif
		}

		// Go over each column, and check if it passes the filter
		for(auto col_it = map.begin(); col_it != map.end(); ++col_it)
		{
			const std::string & key = col_it->first;
			if(!column_passes_lens_filter(key,col_it->second.at(i)))
			{
				bad_indices.push_back(i);
				break;
			}
		}
	}

	return bad_indices;
}

std::vector<size_t> get_bad_sources(const brgastro::table_map_t<std::string> & map, const std::vector<std::vector<bool>> good_pixels)
{
	// Determine the number of rows
	const size_t num_rows = map.begin()->second.size();

	std::vector<size_t> bad_indices;
	bad_indices.reserve(num_rows);

	for(size_t i=0; i<num_rows; ++i)
	{
		// Check if this object is within the mask
		const unsigned xp = brgastro::round_int(boost::lexical_cast<double>(map.at("Xpos").at(i)));
		const unsigned yp = brgastro::round_int(boost::lexical_cast<double>(map.at("Ypos").at(i)));

		if(!is_good_position(xp,yp,good_pixels))
		{
			bad_indices.push_back(i);

#ifdef WARN_MASK_MISMATCH
			// Check that the mask value in the table agrees
			if(boost::lexical_cast<unsigned>(map.at("MASK").at(i))<=1)
			{
				std::cerr << "WARNING: Mask value mismatch for source index " << i << " at position (" << xp << ", " << yp << ")." << std::endl;
				std::cerr << "Table's mask value is " << map.at("MASK").at(i) << ", but saved value is 'false'.\n";
			}
#endif

			continue;
		}
		else
		{
#ifdef WARN_MASK_MISMATCH
			// Check that the mask value in the table agrees
			if(boost::lexical_cast<unsigned>(map.at("MASK").at(i))>1)
			{
				std::cerr << "WARNING: Mask value mismatch for source index " << i << " at position (" << xp << ", " << yp << ")." << std::endl;
				std::cerr << "Table's mask value is " << map.at("MASK").at(i) << ", but saved value is 'true'.\n";
			}
#endif
		}

		// Go over each column, and check if it passes the filter
		for(auto col_it = map.begin(); col_it != map.end(); ++col_it)
		{
			const std::string & key = col_it->first;
			if(!column_passes_source_filter(key,col_it->second.at(i)))
			{
				bad_indices.push_back(i);
				break;
			}
		}
	}

	return bad_indices;
}

bool column_passes_lens_filter(const std::string & col_name,const std::string & value)
{
	if(!column_passes_global_filter(col_name,value)) return false;
	return true;
	if(col_name=="Z_B")
	{
		double dval = boost::lexical_cast<double>(value);
		return (dval <= 1.3) && (dval >= 0.2);
	}
	else
	{
		return true;
	}
}
bool column_passes_source_filter(const std::string & col_name,const std::string & value)
{
	if(!column_passes_global_filter(col_name,value)) return false;
	if(col_name=="Z_B")
	{
		return boost::lexical_cast<double>(value) >= 0.2;
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
		return (dval <= 4);
	}
	else if(col_name=="MAG_r")
	{
		double dval = std::fabs(boost::lexical_cast<double>(value));
		return (dval <= 24.7);
	}
	else if(col_name=="CHI_SQUARED_BPZ")
	{
		double dval = boost::lexical_cast<double>(value);
		return (dval <= 2);
	}
//	else if(col_name=="ODDS")
//	{
//		double dval = boost::lexical_cast<double>(value);
//		return (dval >= 0.8);
//	}
	else
	{
		return true;
	}
}

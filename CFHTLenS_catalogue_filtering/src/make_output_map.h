/**********************************************************************\
 @file make_output_map.h
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


// body file: make_output_map.cpp


#ifndef _BRG_MAKE_OUTPUT_MAP_H_INCLUDED_
#define _BRG_MAKE_OUTPUT_MAP_H_INCLUDED_

#include <functional>
#include <map>
#include <string>
#include <vector>

#include "brg/container/table_typedefs.hpp"

brgastro::table_map_t<std::string> make_output_map(const brgastro::table_map_t<std::string> & map,
		const std::vector<size_t> & filtered_indices, const brgastro::header_t & header_columns,
		const std::map<std::string,std::function<double(double)>> & conversions);


#endif // _BRG_MAKE_OUTPUT_MAP_H_INCLUDED_

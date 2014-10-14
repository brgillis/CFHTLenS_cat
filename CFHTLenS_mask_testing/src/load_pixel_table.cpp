/**********************************************************************\
 @file load_pixel_table.cpp
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

#include <iostream>
#include <memory>
#include <string>
#include <stdexcept>
#include <unistd.h>
#include <vector>
#include <wait.h>

#include <CCfits/CCfits>

#include "load_pixel_table.h"

std::string unpack_fits(const std::string & filename)
{
	constexpr size_t suffix_lengthp1=4;
	constexpr const char suffix[suffix_lengthp1] = ".fz";
	constexpr size_t suffix_length = suffix_lengthp1-1;
	const size_t pos = filename.rfind(suffix);

	if(pos==filename.size()-suffix_length)
	{
		std::string unpacked_name(filename);
		unpacked_name.erase(pos,suffix_length);

		pid_t pid;

		delete_file(unpacked_name);

		const char *unpacker("/disk2/brg/bin/funpack");

		// Spawn a subprocess to unpack the file
		switch (pid = fork())
		{
		case -1:
			throw std::runtime_error("Fork failed in unpack_fits");
		case 0:
			/* This is processed by the child */
			execl(unpacker, unpacker, filename.c_str(), nullptr);
			std::cerr << "Error executing unpacker in unpack_fits.\n";
			_exit(1);
		default:
			/* This is processed by the parent */
			int status=0;
			waitpid(0,&status,0);
			if(status)
				throw std::runtime_error("Error executing unpacker in unpack_fits.\n");
			break;
		}

		return unpacked_name;
	}
	else
	{
		// It has the wrong suffix
		return filename;
	}
}

void delete_file(const std::string & filename)
{
	pid_t pid;
	const char *deleter("/bin/rm");

	// Spawn a subprocess to delete the unpacked file
	switch (pid = fork())
	{
	case -1:
		throw std::runtime_error("Fork failed in delete_file");
	case 0:
		/* This is processed by the child */
		execl(deleter, deleter, "-f", filename.c_str(), nullptr);
		std::cerr << "Error executing rm in delete_file.\n";
		_exit(1);
	default:
		/* This is processed by the parent */
		int status=0;
		waitpid(0,&status,0);
		if(status)
			throw std::runtime_error("Error executing rm in delete_file.\n");
		break;
	}
}

std::vector<std::vector<bool>> load_pixel_table(const std::string & filename)
{
	using namespace CCfits;

#ifndef NDEBUG
	FITS::setVerboseMode(true);
#endif

	std::string unpacked_filename = unpack_fits(filename);

	std::unique_ptr<FITS> pInfile(new FITS(unpacked_filename,Read,true));

	PHDU& image = pInfile->pHDU();

	std::valarray<unsigned short> contents;
	// read all user-specifed, coordinate, and checksum keys in the image
	image.readAllKeys();

	image.read(contents);

	size_t ax1(image.axis(0));
	size_t ax2(image.axis(1));

	std::vector<std::vector<bool>> result(ax1,std::vector<bool>(ax2));

	for (size_t j = 0; j < ax2; ++j)
	{
		for (size_t i=0; i < ax1; ++i)
		{
			auto & v = contents[ax1*j+i];
			if(v==0)
			{
				result[i][j]=true;
			}
			else
			{
				result[i][j]=false;
			}
		}
	}

	pInfile.release();

	delete_file(unpacked_filename);

	return result;
}


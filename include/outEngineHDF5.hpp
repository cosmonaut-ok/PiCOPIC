/*
 * This file is part of the PiCoPiC distribution (https://github.com/cosmonaut-ok/PiCoPiC).
 * Copyright (c) 2020 Alexander Vynnyk.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _OUT_ENGINE_HDF5_HPP_
#define _OUT_ENGINE_HDF5_HPP_

#include <sys/resource.h>

#include "defines.hpp"
#include "loguru.hpp"
#include "msg.hpp"
#include "lib.hpp"
#include "outEngine.hpp"

#include "H5Cpp.h"

using namespace std;

class OutEngineHDF5 : public OutEngine
{
public:
  OutEngineHDF5 () {};
  OutEngineHDF5 (string a_path, string a_subpath, unsigned int a_shape, int *a_size,
                  bool a_append, bool a_compress, unsigned int a_compress_level);

  void write_rec(string a_name, Grid<double> data);
  void write_vec(string a_name, Grid<double> data);
  void write_dot(string a_name, Grid<double> data);
  void write_1d_vector(string a_name, vector<double> data);
  void write_metadata(string metadata);

private:
  string datafile_name;
};

#endif // end of _OUT_ENGINE_HDF5_HPP_

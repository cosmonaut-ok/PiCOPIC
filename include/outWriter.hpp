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

#ifndef _OUT_WRITER_HPP_
#define _OUT_WRITER_HPP_

#include <string>
#include <vector>

#include "algo/grid.hpp"
#include "algo/grid3d.hpp"
#include "timeSim.hpp"

#ifdef ENABLE_HDF5
#include "outEngine/outEngineHDF5.hpp"
#endif // ENABLE_HDF5

class OutWriter
{
private:
#ifdef ENABLE_HDF5
  OutEngineHDF5 engine;
#endif // end of ENABLE_HDF5

  TimeSim *time;
  unsigned short schedule;
  unsigned short shape; // 0 - cube; 1 - rec; 2 - vec; 3 - dot
  std::vector<short> position; // { a_start, b_start, c_start, ... a_end, b_end, c_end ... } for rectangles and cubes - diagonal points
  Grid<double> *values;
  std::string path;

public:
  OutWriter () {};
  OutWriter ( std::string _path, unsigned short _shape,
              std::vector<short> _position, std::vector<size_t> _engine_offset,
              unsigned short _schedule, bool _append, unsigned short _compress,
              TimeSim *_time, Grid<double> *_values );

#ifdef ENABLE_HDF5
  OutWriter ( HighFive::File* _file, std::string _path, unsigned short _shape,
              std::vector<short> _position, std::vector<size_t> _engine_offset,
              unsigned short _schedule, bool _append, unsigned short _compress,
              TimeSim *_time, Grid<double> *_values );
#endif // ENABLE_HDF5

  void operator()();

#ifdef ENABLE_HDF5
  HighFive::File *hdf5_file;
#endif
};

#endif // end of _OUT_WRITER_HPP_

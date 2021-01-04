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

#include <highfive/H5File.hpp>

#include "outEngine.hpp"

class OutEngineHDF5 : public OutEngine
{
public:
  OutEngineHDF5 () {};

  OutEngineHDF5 ( std::string _path,
                  std::vector<size_t> _size,
                  std::vector<size_t> _offset,
                  bool _append, unsigned short _compress );

  OutEngineHDF5 ( HighFive::File* _file, std::string _path,
                  std::vector<size_t> _size,
                  std::vector<size_t> _offset,
                  bool _append, unsigned short _compress );

  void create_dataset();
  void extend_dataset(size_t num); // extend dataset to number of slices
  void create_path(); // create group
  void write_metadata(picojson::value _metadata);

  void write_cub(size_t _slice, vector<vector<vector<double>>> data); // cube
  void write_rec(size_t _slice, vector<vector<double>> data);         // rectangle
  void write_vec(size_t _slice, vector<double> data);                 // vector
  void write_dot(size_t _slice, double data);                         // dot

// private:
  HighFive::File *data_file;
};

#endif // end of _OUT_ENGINE_HDF5_HPP_

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

  OutEngineHDF5 ( std::string _group,
                    std::vector<size_t> _offset,
                    bool _append, unsigned short _compress );

  OutEngineHDF5 ( HighFive::File* _file, std::string _group,
                    std::vector<size_t> _offset,
                    bool _append, unsigned short _compress );

  void create_dataset(std::string _name, vector<unsigned int> _dims);
  void create_path();
  void write_metadata(std::string _metadata);

  // void write_cub(string _name, Grid3D<double> data); // cube
  void write_rec(string _name, vector<vector<double>> data);   // rectangle
  void write_vec(string _name, vector<double> data);           // vector
  void write_dot(string _name, double data);                   // dot

// private:
  HighFive::File *data_file;
};

#endif // end of _OUT_ENGINE_HDF5_HPP_

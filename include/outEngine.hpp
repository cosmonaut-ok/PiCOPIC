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

#ifndef _OUT_ENGINE_HPP_
#define _OUT_ENGINE_HPP_

#include <string>
#include <vector>

#include "algo/grid.hpp"

class OutEngine
{
public:
  OutEngine () {};
  OutEngine ( std::string _group, std::vector<size_t> _offset,
                bool _append, unsigned short _compress )
  {
    group = _group;
    offset = _offset;
    append = _append;
    compress = _compress;
  };

  virtual void create_dataset(std::string _name, vector<unsigned int> _dims) = 0;
  virtual void create_path() = 0;
  virtual void write_metadata(std::string _metadata) = 0;

  // virtual void write_cub(string _name, Grid3D<double> data) = 0; // cube
  virtual void write_rec(string _name, vector<vector<double>> data) = 0;   // rectangle
  virtual void write_vec(string _name, vector<double> data) = 0;           // vector
  virtual void write_dot(string _name, double data) = 0;                   // dot

protected:
  std::string group;
  std::vector<size_t> offset;
  bool append;
  unsigned short compress;
};

#endif // end of _OUT_ENGINE_HPP_

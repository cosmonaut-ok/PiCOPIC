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

#include "picojson.h"

#include "algo/grid.hpp"

class OutEngine
{
public:
  OutEngine () {};
  OutEngine ( std::string _path, vector<size_t> _size,
              std::vector<size_t> _offset,
              bool _append, unsigned short _compress )
  {
    path = _path;
    size = _size;
    offset = _offset;
    append = _append;
    compress = _compress;
  };

  virtual void create_dataset() = 0; // extend dataset to number of slices
  virtual void extend_dataset(size_t num) = 0; // extend dataset to number of slices
  virtual void create_path() = 0;
  virtual void write_metadata(picojson::value _metadata) = 0;

  virtual void write_cub(size_t _slice, vector<vector<vector<double>>> data) = 0;   // rectangle
  virtual void write_rec(size_t _slice, vector<vector<double>> data) = 0;   // rectangle
  virtual void write_vec(size_t _slice, vector<double> data) = 0;           // vector
  virtual void write_dot(size_t _slice, double data) = 0;                   // dot

protected:
  std::string path;
  std::vector<size_t> size;
  std::vector<size_t> offset;
  bool append;
  unsigned short compress;
};

#endif // end of _OUT_ENGINE_HPP_

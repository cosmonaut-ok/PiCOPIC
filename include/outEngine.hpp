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
#include "algo/grid.hpp"

#include "defines.hpp"
#include "msg.hpp"

using namespace std;

class OutEngine
{
public:
  string path;
  string subpath;
  unsigned int shape;
  int size[4];
  bool append;
  bool compress;
  unsigned int compress_level;


public:
  OutEngine () {};
  OutEngine (string _path, string _subpath, unsigned int _shape, int *_size,
             bool _append, bool _compress, unsigned int _compress_level)
  {
    path = _path;
    subpath = _subpath;
    shape = _shape;
    append = _append;
    compress = _compress;
    compress_level = _compress_level;

    for (unsigned int i = 0; i < 4; ++i)
      size[i] = _size[i];
  };

  ~OutEngine () {};

  virtual void write_rec(string _name, Grid<double> data) = 0;
  virtual void write_vec(string _name, Grid<double> data) = 0;
  virtual void write_dot(string _name, Grid<double> data) = 0;
  virtual void write_1d_vector(string _name, vector<double> data) = 0;
  virtual void write_metadata(string _metadata) = 0;
};

#endif // end of _OUT_ENGINE_HPP_

/* 
 * This file is part of the PiCOPIC distribution (https://github.com/cosmonaut-ok/PiCOPIC).
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

#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h> /* floor */

#include "grid.hpp"

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
  OutEngine (string a_path, string a_subpath, unsigned int a_shape, int *a_size,
             bool a_append, bool a_compress, unsigned int a_compress_level)
  {
    path = a_path;
    subpath = a_subpath;
    shape = a_shape;
    append = a_append;
    compress = a_compress;
    compress_level = a_compress_level;

    for (unsigned int i = 0; i < 4; ++i)
      size[i] = a_size[i];
  };

  ~OutEngine () {};

  virtual void write_rec(string a_name, Grid<double> data) = 0;
  virtual void write_vec(string a_name, Grid<double> data) = 0;
  virtual void write_dot(string a_name, Grid<double> data) = 0;
  virtual void write_1d_vector(string a_name, vector<double> data) = 0;
  virtual void write_metadata(string metadata) = 0;
};

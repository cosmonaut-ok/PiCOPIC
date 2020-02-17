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

#include <typeinfo>

#include "outEnginePlain.hpp"

OutEnginePlain::OutEnginePlain (string a_path, string a_subpath, unsigned int a_shape, int *a_size,
                                bool a_append, bool a_compress, unsigned int a_compress_level)
  : OutEngine (a_path, a_subpath, a_shape, a_size, a_append, a_compress, a_compress_level)
{
  // make root directory
  lib::make_directory(path + "/" + subpath);

  metadata_file = "metadata.json";

  if (compress)
    LOG_WARN("compression is not supported by plaintext output engine");
}

void OutEnginePlain::write_rec(string a_name, Grid<double> data)
{
  std::ofstream::openmode omode;
  if (append)
    omode = ios::app;
  else
    omode = ios::trunc;

  ofstream out_val(path + "/" + subpath + "/" + a_name + ".dat", omode);

  for (int i = size[0]; i < size[2]; i++)
    for (int j = size[1]; j < size[3]; j++)
      out_val << data(i, j) << " ";

  out_val.close();
}

void OutEnginePlain::write_vec(string a_name, Grid<double> data)
{
  std::ofstream::openmode omode;
  if (append)
    omode = ios::app;
  else
    omode = ios::trunc;

  ofstream out_val(path + "/" + subpath + "/" + a_name + ".dat", omode);

  // vector by Z-component with fixed r (r_begin)
  if (size[1] == -1 && size[0] == -1 && size[3] == -1)
  {
    for (unsigned int i = 0; i < data.y_size; i++)
      out_val << data(size[2], i) << " ";
  }
  // vector by R-component with fixed z (z_begin)
  else if (size[0] == -1 && size[1] == -1 && size[2] == -1)
  {
    for (unsigned int i = 0; i < data.x_size; i++)
      out_val << data(i, size[3]) << " ";
  }
  else
    LOG_CRIT("Incorrect shape for vector output", 1);

  out_val.close();
}

void OutEnginePlain::write_dot(string a_name, Grid<double> data)
{
    std::ofstream::openmode omode;
  if (append)
    omode = ios::app;
  else
    omode = ios::trunc;

  ofstream out_val(path + "/" + subpath + "/" + a_name + ".dat", omode);

  out_val << data(size[2], size[3]) << endl;

  out_val.close();
}

void OutEnginePlain::write_1d_vector(string a_name, vector<double> data)
{
  std::ofstream::openmode omode;
  if (append)
    omode = ios::app;
  else
    omode = ios::trunc;

  ofstream out_val(path + "/" + subpath + "/" + a_name + ".dat", omode);

  for (auto i = data.begin(); i != data.end(); ++i)
    out_val << *i << " ";

  out_val.close();
}

void OutEnginePlain::write_metadata(string metadata)
{
  ofstream out_val(path + "/" + metadata_file, ios::out);
  out_val << metadata;
  out_val.close();
}

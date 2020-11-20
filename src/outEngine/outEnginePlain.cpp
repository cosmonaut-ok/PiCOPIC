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

#include "outEngine/outEnginePlain.hpp"

OutEnginePlain::OutEnginePlain (string _path, string _subpath, unsigned int _shape, int *_size,
                                bool _append, bool _compress, unsigned int _compress_level)
  : OutEngine (_path, _subpath, _shape, _size, _append, _compress, _compress_level)
{
  // make root directory
  algo::common::make_directory(path + "/" + subpath);

  metadata_file = "metadata.json";

  if (compress)
    LOG_S(WARNING) << "compression is not supported by plaintext output engine";
}

void OutEnginePlain::write_rec(string _name, Grid<double> data)
{
  std::ofstream::openmode omode;
  if (append)
    omode = ios::app;
  else
    omode = ios::trunc;

  ofstream out_val(path + "/" + subpath + "/" + _name + ".dat", omode);

  for (int i = size[0]; i < size[2]; i++)
    for (int j = size[1]; j < size[3]; j++)
      out_val << data(i, j) << " ";

  out_val.close();
}

void OutEnginePlain::write_vec(string _name, Grid<double> data)
{
  std::ofstream::openmode omode;
  if (append)
    omode = ios::app;
  else
    omode = ios::trunc;

  ofstream out_val(path + "/" + subpath + "/" + _name + ".dat", omode);

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
    LOG_S(FATAL) << "Incorrect shape for vector output";

  out_val.close();
}

void OutEnginePlain::write_dot(string _name, Grid<double> data)
{
    std::ofstream::openmode omode;
  if (append)
    omode = ios::app;
  else
    omode = ios::trunc;

  ofstream out_val(path + "/" + subpath + "/" + _name + ".dat", omode);

  out_val << data(size[2], size[3]) << endl;

  out_val.close();
}

void OutEnginePlain::write_1d_vector(string _name, vector<double> data)
{
  std::ofstream::openmode omode;
  if (append)
    omode = ios::app;
  else
    omode = ios::trunc;

  ofstream out_val(path + "/" + subpath + "/" + _name + ".dat", omode);

  for (auto i = data.begin(); i != data.end(); ++i)
    out_val << *i << " ";

  out_val.close();
}

void OutEnginePlain::write_metadata(string _metadata)
{
  ofstream out_val(path + "/" + metadata_file, ios::out);
  out_val << _metadata;
  out_val.close();
}

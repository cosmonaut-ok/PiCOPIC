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

#ifndef _OUT_ENGINE_PLAIN_HPP_
#define _OUT_ENGINE_PLAIN_HPP_

#include "outEngine.hpp"

using namespace std;

class OutEnginePlain : public OutEngine
{
public:
  OutEnginePlain () {};
  OutEnginePlain (string _path, string _name, unsigned int _shape, int *_size,
                  bool _append, bool _compress, unsigned int _compress_level);

  void write_rec(string _name, Grid<double> data);
  void write_vec(string _name, Grid<double> data);
  void write_dot(string _name, Grid<double> data);
  void write_1d_vector(string _name, vector<double> data);
  void write_metadata(string _metadata);

private:
  string metadata_file;
};
#endif // end of _OUT_ENGINE_PLAIN_HPP_

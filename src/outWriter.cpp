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

#include "outWriter.hpp"

using namespace std;

OutWriter::OutWriter ( string _path, unsigned short _shape,
                       vector<short> _position, vector<size_t> _engine_offset,
                       unsigned short _schedule, bool _append, unsigned short _compress,
                       TimeSim *_time, Grid<double> *_values )
  : time(_time)
{
  path = _path;
  schedule = _schedule;
  shape = _shape;
  position = _position;
  values = _values;
}

#ifdef ENABLE_HDF5
OutWriter::OutWriter ( HighFive::File* _file, string _path, unsigned short _shape,
                       vector<short> _position, vector<size_t> _engine_offset,
                       unsigned short _schedule, bool _append, unsigned short _compress,
                       TimeSim *_time, Grid<double> *_values )
  : OutWriter ( _path, _shape, _position, _engine_offset, _schedule,
                _append, _compress, _time, _values)
{
  engine = OutEngineHDF5 (_file, _path, {128, 512}, _engine_offset, _append, _compress);
  hdf5_file = _file;
  // engine.data_file = hdf5_file;
}
#endif // ENABLE_HDF5

void OutWriter::operator()()
{
  int current_time_step = ceil(time->current / time->step);
  int is_run = current_time_step % schedule;

  if (is_run == 0)
  {
    size_t slice = (size_t)(ceil(current_time_step / schedule));
    LOG_S(MAX) << "Launching " << path << "/" << slice;

    switch (shape)
    {
    case 0: // rectangle shape
    {
      vector<vector<double>> val (position[2] - position[0],
                                  vector<double> (position[3] - position[1], 0));

      for (short i = position[0]; i < position[2]; ++i)
        for (short j = position[1]; j < position[3]; ++j)
          val[i - position[0]][j - position[1]] = (*values)(i, j);

      engine.write_rec(slice, val);
      break;
    }
    case 1: // column shape
    {
      vector<double> val;
      short col_offset = position[3];
      for (unsigned int i = 0; i < values->x_size; ++i)
        val.push_back( (*values)(i, col_offset) );

      engine.write_vec(slice, val);
      break;
    }
    case 2: // row shape
    {
      vector<double> val;
      short row_offset = position[2];
      for (unsigned int i = 0; i < values->y_size; ++i)
        val.push_back( (*values)(row_offset, i) );

      engine.write_vec(slice, val);
      break;
    }
    case 3: // dot shape
    {
      engine.write_dot(slice, (*values)(position[2], position[3]));
      break;
    }
    }
  }
}

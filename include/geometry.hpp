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

#ifndef _GEOMETRY_HPP_
#define _GEOMETRY_HPP_

#include <vector>

#include "defines.hpp"
#include "constant.hpp"
#include "msg.hpp"

struct Geometry
{
public:
  std::vector<double> size;
  std::vector<double> global_size;
  std::vector<double> cell_size;

  std::vector<size_t> cell_dims; // offsets of all base edges
  std::vector<size_t> cell_amount;
  std::vector<size_t> global_cell_amount;

  std::vector<size_t> domains_amount;

  // bool walls[4]; // ex. [true, true, false, false]
  std::vector<bool> walls;

  //! comparative pml lengths on walls
  //! r=0, z=0, r=wall, z=wall
  std::vector<size_t> pml_size;
  std::vector<double> pml_sigma;

  // Geometry constructors
  Geometry ( std::vector<double> _size,
             std::vector<double> _global_size,
             std::vector<size_t> _cell_dims,
             std::vector<size_t> _pml_size,
             std::vector<double> _pml_sigma,
             std::vector<bool> _walls )
  {
    size = _size;
    global_size = _global_size;
    cell_dims = _cell_dims;

    // cell dimensions format is {a0, b0, c0, ... a1, b1, c1, ...}
    // so beginings and endings of the edge points are separated
    // by half of array length
    size_t half_length = (size_t)(cell_dims.size() / 2.);
    for (size_t i = 0; i < half_length; ++i)
      cell_amount.push_back(cell_dims[i+half_length] - cell_dims[i]);

    for (unsigned int i = 0; i < size.size(); ++i)
      cell_size.push_back(size[i] / cell_amount[i]);

    for (size_t i = 0; i < _global_size.size(); ++i)
      global_cell_amount.push_back(global_size[i] / cell_size[i]);

    walls = _walls;
    pml_size = _pml_size;
    pml_sigma = _pml_sigma;
  };

  Geometry ( std::vector<double> _size,
             std::vector<size_t> _cell_dims,
             std::vector<size_t> _pml_size,
             std::vector<double> _pml_sigma,
             std::vector<bool> _walls)
    : Geometry ( _size, _size, _cell_dims, _pml_size, _pml_sigma, _walls)
  {};

  Geometry ( std::vector<double> _size,
             std::vector<double> _global_size,
             std::vector<size_t> _cell_dims,
             std::vector<bool> _walls)
    : Geometry ( _size, _global_size, _cell_dims, {0, 0, 0, 0}, {0, 0}, _walls)
  {};

  Geometry ( std::vector<double> _size,
             std::vector<size_t> _cell_dims,
             std::vector<bool> _walls)
    : Geometry ( _size, _size, _cell_dims, {0, 0, 0, 0}, {0, 0}, _walls)
  {};

  Geometry() {};

  ~Geometry() {};

  double dist_x_to_end (size_t cell_num)
  {
    size_t cell_offset_rel_end = global_cell_amount[0] - cell_dims[0];

    return (double)(cell_offset_rel_end - cell_num) * cell_size[0];
  };

  double dist_y_to_end (size_t cell_num)
  {
    size_t cell_offset_rel_end = global_cell_amount[1] - cell_dims[1];

    return (double)(cell_offset_rel_end - cell_num) * cell_size[1];
  };

  // double dist_x_to_end (double pos)
  // {
  //   return global_size[0] - pos;
  // };

  // double dist_y_to_end (double pos)
  // {
  //   return global_size[1] - pos;
  // };

  // double cells_begin_to_x (size_t num)
  // {
  //   // how many cells from the begining of global simulation area
  //   // by 2nd dimension
  //   return cell_dims[0] + num;
  // };

  // double cells_begin_to_y (size_t num)
  // {
  //   // how many cells from the begining of global simulation area
  //   // by 2nd dimension
  //   return cell_dims[1] + num;
  // };

  // size_t cells_x_to_end (size_t num)
  // {
  //   // how many cells left to the end of global simulation area
  //   // by 1st dimension
  //   // NOTE: num is local for the domain
  //   size_t x_cell_abs = num + cell_dims[0]; // assume, coords are [x,y,z,...]
  //   return global_cell_amount[0] - x_cell_abs;
  // };

  // size_t cells_y_to_end (size_t num)
  // {
  //   // how many cells left to the end of global simulation area
  //   // by 2nd dimension
  //   // NOTE: num is local for the domain
  //   size_t y_cell_abs = num + cell_dims[1]; // assume, coords are [x,y,z,...]
  //   return global_cell_amount[1] - y_cell_abs;
  // };

  // size_t cell_num (double pos)
  // {
  //   // cell number by radius in cylindrical frame
  //   return (size_t)floor(pos / cell_size[0]);
  // };

  // size_t cell_num_local_x (double pos)
  // {
  //   // cell number by radius in cylindrical frame
  //   return (size_t)(floor(pos / cell_size[0]) - cell_dims[0]);
  // };

  // size_t cell_num_local_y (double pos)
  // {
  //   // cell number by radius in cylindrical frame
  //   return (size_t)(floor(pos / cell_size[0]) - cell_dims[1]);
  // };

  // double cell_vol_cyl(size_t num_r)
  // {
  //   // cell colume in cylindrical frame
  //   return cell_size[1] * cell_size[0] * cell_size[0] * (2. * num_r + 1);
  // };
};

// some geometry-related macros
//! \f$ ( \pi \times (dr * (i+0.5))^2 - \pi \times (dr * (i-0.5))^2 ) * dz \f$
#define CELL_VOLUME(i, dr, dz) constant::PI * (dz) * (dr) * (dr) * (2.0 * i + 1)

//! volume of the cylindrical ring (internal cylinder on r1 is cut out)
#define CYL_RNG_VOL(z, r1, r2) constant::PI * (z) * ((r2) * (r2) - (r1) * (r1))

//! volume of the cylinder
// #define CYL_VOL(z, r) PI * (z) * (r) * (r) / 4.
#define CYL_VOL(z, r) constant::PI * (z) * (r) * (r)

// #define PARTICLE_VOLUME(x,y) (PI * dz * dr * dr * 2.0 * i)

//! get cell number by 'radius'
#define CELL_NUMBER(position, dx) (int)floor((position) / (dx))

#endif // end of _GEOMETRY_HPP_

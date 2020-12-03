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
  std::vector<size_t> cell_dims;
  std::vector<size_t> cell_amount;
  std::vector<double> cell_size;

  // bool walls[4]; // ex. [true, true, false, false]
  std::vector<bool> walls;

  //! comparative pml lengths on walls
  //! r=0, z=0, r=wall, z=wall
  std::vector<double> pml_length;
  std::vector<double> pml_sigma;
  std::vector<size_t> domains_amount;

  // Geometry constructors
  Geometry ( std::vector<double> _size,
             std::vector<size_t> _cell_dims,
             std::vector<double> _pml_length,
             std::vector<double> _pml_sigma,
             std::vector<bool> _walls )
  {
    size = _size;
    cell_dims = _cell_dims;

    // cell dimensions format is {a0, b0, c0, ... a1, b1, c1, ...}
    // so beginings and endings of the edge points are separated
    // by half of array length
    size_t half_length = (size_t)(cell_dims.size() / 2.);
    for (size_t i = 0; i < half_length; ++i)
      cell_amount.push_back(cell_dims[i+half_length] - cell_dims[i]);

    for (unsigned int i = 0; i < size.size(); ++i)
      cell_size.push_back(size[i] / cell_amount[i]);

    walls = _walls;
    pml_length = _pml_length;
    pml_sigma = _pml_sigma;
  };

  Geometry ( std::vector<double> _size,
             std::vector<size_t> _cell_dims,
             std::vector<bool> _walls)
    : Geometry ( _size, _cell_dims, {0, 0, 0, 0}, {0, 0}, _walls)
  {};

  Geometry() {};

  ~Geometry() {};
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

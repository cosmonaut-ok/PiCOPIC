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

#include <cmath>
#include "constant.hpp"

using namespace constant;

class Geometry
{
public:
  double r_size;
  double z_size;

  int top_r_grid_number;
  int bottom_r_grid_number;
  int left_z_grid_number;
  int right_z_grid_number;

  int r_grid_amount;
  int z_grid_amount;

  double r_cell_size;
  double z_cell_size;

  bool walls[4]; // ex. [true, true, false, false]

  //! comparative pml lengths on walls
  //! r=0, z=0, r=wall, z=wall
  double pml_length[4] = {0, 0, 0, 0}; // [0, 0, 0, 0]
  double pml_sigma[2] = {0, 0}; // [0, 0]

  bool is_near_z_axis;

  unsigned int domains_by_r;
  unsigned int domains_by_z;

  Geometry(double rs, double zs,
           int bot_ngr, int top_ngr, int left_ngz, int right_ngz,
           double pml_l_z0,
           double pml_l_zwall, double pml_l_rwall, double pml_sigma1,
           double pml_sigma2,
           bool wall_r0, bool wall_z0, bool wall_rr, bool wall_zz);
  Geometry() {};

  ~Geometry(void);

private:
  void set_pml(double comparative_l_1, double comparative_l_2, double comparative_l_3,
               double sigma1, double sigma2);


};

//! *_OVERLAY is for shifting cell param from overlay to normal
//! *_OVERLAY_REVERSE is for shifting cell param from normal to overlay

// some geometry-related macros
//! \f$ ( \pi \times (dr * (i+0.5))^2 - \pi \times (dr * (i-0.5))^2 ) * dz \f$
#define CELL_VOLUME(i, dr, dz) PI * (dz) * (dr) * (dr) * (2.0 * i + 1)

//! volume of the cylindrical ring (internal cylinder on r1 is cut out)
#define CYL_RNG_VOL(z, r1, r2) PI * (z) * ((r2) * (r2) - (r1) * (r1))

//! volume of the cylinder
// #define CYL_VOL(z, r) PI * (z) * (r) * (r) / 4.
#define CYL_VOL(z, r) PI * (z) * (r) * (r)

// #define PARTICLE_VOLUME(x,y) (PI * dz * dr * dr * 2.0 * i)

//! get cell number by 'radius'
#define CELL_NUMBER(position, dx) (int)floor((position) / (dx))

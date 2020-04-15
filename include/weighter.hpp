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

#ifndef _WEIGHTER_HPP_
#define _WEIGHTER_HPP_

#include "geometry.hpp"
#include "grid.hpp"
#include "specieP.hpp"

template <typename T>
void weight_cylindrical(Geometry *geometry, Grid<T> *grid, double pos_r, double pos_z, T weight_value)
{
  T dr = (T)geometry->r_cell_size;
  T dz = (T)geometry->z_cell_size;
  int bottom_shift = geometry->bottom_r_grid_number;
  int left_shift = geometry->left_z_grid_number;
  T r1, r2, r3; // radiuses
  T dz1, dz2; // longitudes

  T rho_w_value, v_0, value;

  // finding number of i and k cell. example: dr = 0.5; r = 0.4; i =0
  unsigned int r_i = CELL_NUMBER(pos_r, dr);
  unsigned int z_k = CELL_NUMBER(pos_z, dz);
  unsigned int r_i_shift = r_i - bottom_shift;
  unsigned int z_k_shift = z_k - left_shift;

  // volumes
  T v_1 = CELL_VOLUME(r_i, dr, dz);
  T v_2 = CELL_VOLUME(r_i + 1, dr, dz);

  if (pos_r > dr)
  {
    r1 =  pos_r - 0.5 * dr;
    r2 = (r_i + 0.5) * dr;
    r3 = pos_r + 0.5 * dr;
    v_0 = 2. * PI * dz * dr * pos_r;
    dz1 = (z_k + 0.5) * dz - (pos_z - 0.5 * dz);
    dz2 = (pos_z + 0.5 * dz) - (z_k + 0.5) * dz;
    rho_w_value = weight_value / v_0;

    // weighting in ro[i][k] cell
    value = CYL_RNG_VOL(dz1, r1, r2) / v_1;
    grid->inc(r_i_shift, z_k_shift, rho_w_value * value);

    // weighting in ro[i + 1][k] cell
    value = CYL_RNG_VOL(dz1, r2, r3) / v_2;
    grid->inc(r_i_shift + 1, z_k_shift, rho_w_value * value);

    // weighting in ro[i][k + 1] cell
    value = CYL_RNG_VOL(dz2, r1, r2) / v_1;
    grid->inc(r_i_shift, z_k_shift+1, rho_w_value * value);

    // weighting in ro[i + 1][k + 1] cell
    value = CYL_RNG_VOL(dz2, r2, r3) / v_2;
    grid->inc(r_i_shift + 1, z_k_shift + 1, rho_w_value * value);
  }
  else if (pos_r <= dr / 2.)
  {
    r_i = 0;
    r1 =  0.;
    r2 = (r_i + 0.5) * dr;
    r3 = pos_r + 0.5 * dr;
    dz1 = (z_k + 0.5) * dz - (pos_z - 0.5 * dz);
    dz2 = (pos_z + 0.5 * dz) - (z_k + 0.5) * dz;
    v_0 = PI * dz * (2. * pos_r * pos_r + dr * dr / 2.);
    rho_w_value = weight_value / v_0;

    // weighting in ro[i][k] cell
    value = PI * dz1 * (dr * dr / 2. - pos_r * dr + pos_r * pos_r) / v_1;

    grid->inc(r_i_shift, z_k_shift, rho_w_value * value);

    // weighting in ro[i + 1][k] cell
    value = CYL_RNG_VOL(dz1, r2, r3) / v_2;
    grid->inc(r_i_shift + 1, z_k_shift, rho_w_value * value);

    // weighting in ro[i][k + 1] cell
    value = PI * dz2 * (dr * dr / 2. - pos_r * dr + pos_r * pos_r) / v_1;
    grid->inc(r_i_shift, z_k_shift + 1, rho_w_value * value);

    // weighting in ro[i + 1][k + 1] cell
    value = CYL_RNG_VOL(dz2, r2, r3) / v_2;
    grid->inc(r_i_shift + 1, z_k_shift + 1, rho_w_value * value);
  }
  else
  {
    r1 = pos_r - 0.5 * dr;
    r2 = (r_i + 0.5) * dr;
    r3 = pos_r + 0.5 * dr;
    dz1 = (z_k + 0.5) * dz - (pos_z - 0.5 * dz);
    dz2 = (pos_z + 0.5 * dz) - (z_k + 0.5) * dz;
    v_0 = 2. * PI * dz * dr * pos_r;
    v_1 = CYL_VOL(dz, dr);
    v_2 = CELL_VOLUME(r_i + 1, dr, dz);

    rho_w_value = weight_value / v_0;

    // weighting in ro[i][k] cell
    value = CYL_RNG_VOL(dz1, r1, r2) / v_1;
    grid->inc(r_i_shift, z_k_shift, rho_w_value * value);

    // weighting in ro[i + 1][k] cell
    value = CYL_RNG_VOL(dz1, r2, r3) / v_2;
    grid->inc(r_i_shift + 1, z_k_shift, rho_w_value * value);

    // weighting in ro[i][k + 1] cell
    value = CYL_RNG_VOL(dz2, r1, r2) / v_1;
    grid->inc(r_i_shift, z_k_shift + 1, rho_w_value * value);

    // weighting in ro[i + 1][k + 1] cell
    value = CYL_RNG_VOL(dz2, r2, r3) / v_2;
    grid->inc(r_i_shift + 1, z_k_shift + 1, rho_w_value * value);
  }
}

#endif // end of _WEIGHTER_HPP_

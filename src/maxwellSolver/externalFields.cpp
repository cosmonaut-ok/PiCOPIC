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

#include "maxwellSolver/externalFields.hpp"

using namespace constant;

ExternalFields::ExternalFields ( Geometry *_geometry, TimeSim *_time,
                                 Current *_current )
  : MaxwellSolver(_geometry, _time, _current)
{
  // emulate dielectric walls
  if (geometry->walls[0]) // r=0
    r_begin = 1;
  else
    r_begin = 0;

  if (geometry->walls[1]) // z=0
    z_begin = 1;
  else
    z_begin = 0;

  if (geometry->walls[2]) // r=r
    r_end = geometry->cell_amount[0] - 1;
  else
    r_end = geometry->cell_amount[0];

  if (geometry->walls[3]) // z=z
    z_end = geometry->cell_amount[1] - 1;
  else
    z_end = geometry->cell_amount[1];

#ifdef ENABLE_PML
    set_pml();
#endif // ENABLE_PML
}

void ExternalFields::set_pml()
{
  double cell_size_r = geometry->cell_size[0];
  double cell_size_z = geometry->cell_size[1];

  // defining lenght of sigma calculation ion left wall
  double lenght_sigma_left = cell_size_z * geometry->pml_size[1];
  double lenght_sigma_right = cell_size_z * geometry->pml_size[3];
  double lenght_sigma_extern = cell_size_r * geometry->pml_size[2];

  double radial_cell_offset = geometry->cell_dims[0] * geometry->cell_size[0];
  double longitudinal_cell_offset = geometry->cell_dims[1] * geometry->cell_size[1];

  double system_radius = geometry->global_size[0];
  double system_longitude =  geometry->global_size[1];

  double delta_sigma = geometry->pml_sigma[1] - geometry->pml_sigma[0];

  // r=r wall
  if (geometry->pml_size[2] != 0)
    for(int i = 0; i < geometry->cell_amount[0]; i++)
      for(int k = 0; k < geometry->cell_amount[1]; k++)
        if (geometry->dist_x_to_end(i) < lenght_sigma_extern)
          sigma.inc(i, k, geometry->pml_sigma[0]
                    + delta_sigma
                    / pow(lenght_sigma_extern, 2)
                    * pow(cell_size_r * (i + 1) + radial_cell_offset
                          - system_radius
                          + lenght_sigma_extern, 2));

  // sigma on z=0 and z=z walls
  if (geometry->pml_size[1] != 0)
    for(int k = 0; k < geometry->cell_amount[1]; k++)
      for(int i = 0; i < geometry->cell_amount[0]; i++)
        // z=0 wall
        if (cell_size_z * k + longitudinal_cell_offset < lenght_sigma_left)
          sigma.inc(i, k, geometry->pml_sigma[0]
                    + delta_sigma
                    / pow(lenght_sigma_left, 2)
                    * pow(lenght_sigma_left - cell_size_z * k - longitudinal_cell_offset, 2));

  // z=z wall
  if (geometry->pml_size[3] != 0)
    for(int k = 0; k < geometry->cell_amount[1]; k++)
      for(int i = 0; i < geometry->cell_amount[0]; i++)
        if (geometry->dist_y_to_end(k) < lenght_sigma_right)
          sigma.inc(i, k, geometry->pml_sigma[0]
                    + delta_sigma
                    / pow(lenght_sigma_right, 2)
                    * pow(cell_size_z * (k + 1) + longitudinal_cell_offset
                          - system_longitude
                          + lenght_sigma_right, 2));
}

void ExternalFields::calc_field_h()
{
  field_h.overlay_set(0);

  double dr = geometry->cell_size[0];
  double dz = geometry->cell_size[1];

  // H_r on outer wall (r=r)
  if (geometry->walls[2])
    for(int k = 0; k < geometry->cell_amount[1]; k++)
    {
      int i = geometry->cell_amount[0] - 1;

      field_h[1].set(i, k, 0);
    }

  // regular case
  for (int i = 0; i < geometry->cell_amount[0]; i++)
    for (int k = 0; k < geometry->cell_amount[1]; k++)
    {
      field_h[0].set(i, k, 0);
      field_h[1].set(i, k, 0);
      field_h[2].set(i, k, 0);
    }
}

void ExternalFields::calc_field_e()
{
  field_e.overlay_set(0);

  double dr = geometry->cell_size[0];
  double dz = geometry->cell_size[1];

  // E at the center axis (r=0) case
  if (geometry->walls[0]) // calculate only at the center axis (r=0)
    for (unsigned int k = z_begin; k < z_end; ++k)
    {
      int i = 0;
      double epsilonx2 = 2 * epsilon(i, k);

#ifdef ENABLE_PML
      double sigma_t = sigma(i, k) * time->step;
#else
      double sigma_t = 0;
#endif // ENABLE_PML

      field_e[0].set(i, k, 0);
      field_e[2].set(i, k, 0);
    }

  // E_z at the left wall (z=0) case
  if (geometry->walls[1]) // calculate only at the left wall (z=0)
    for (unsigned int i = r_begin; i < r_end; ++i)
    {
      int k = 0;
      double epsilonx2 = 2 * epsilon(i, k);

#ifdef ENABLE_PML
      double sigma_t = sigma(i, k) * time->step;
#else
      double sigma_t = 0;
#endif // ENABLE_PML

      field_e[2].set(i, k, 0);
    }

// regular case
  for (unsigned int i = r_begin; i < r_end; i++)
    for (unsigned int k = z_begin; k < z_end; k++)
    {
      double epsilonx2 = 2 * epsilon(i, k);

#ifdef ENABLE_PML
      double sigma_t = sigma(i, k) * time->step;
#else
      double sigma_t = 0;
#endif // ENABLE_PML

      double koef_e = (epsilonx2 - sigma_t) / (epsilonx2 + sigma_t);
      double koef_h = 2 * time->step / (epsilonx2 + sigma_t);

      field_e[0].set(i, k, 0);
      field_e[1].set(i, k, 0);
      field_e[2].set(i, k, 0);
    }
}

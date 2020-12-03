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

#include "maxwellSolver/maxwellSolverYee.hpp"

using namespace constant;

MaxwellSolverYee::MaxwellSolverYee(Geometry *_geometry, TimeSim *_time,
                                   vector<SpecieP *> _species_p, Current *_current)
  : MaxwellSolver(_geometry, _time, _species_p, _current)
{
  // set field_h_at_et to zero-state
  field_h_at_et = Grid3D<double> (geometry->cell_amount[0], geometry->cell_amount[1], 2);
  field_h_at_et = 0.;
  field_h_at_et.overlay_set(0.);

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

void MaxwellSolverYee::set_pml()
{
  // defining lenght of sigma calculation ion left wall
  double lenght_sigma_left = geometry->cell_size[1]
    * (floor(geometry->cell_amount[1] * geometry->pml_length[1]));

  // defining lenght of sigma calculation on right wall
  double lenght_sigma_right = geometry->cell_size[1]
    * (floor(geometry->cell_amount[1] * geometry->pml_length[3]));

  // defining lenght of sigma calculation on z-wall
  double lenght_sigma_extern = geometry->cell_size[0]
    * (floor(geometry->cell_amount[0] * geometry->pml_length[2]));

  double radial_shift = geometry->cell_dims[0] * geometry->cell_size[0];
  double longitudinal_shift = geometry->cell_dims[1] * geometry->cell_size[1];

  // r=r wall
  if (geometry->pml_length[2] != 0)
    for(int i = 0; i < geometry->cell_amount[0]; i++)
      for(int k = 0; k < geometry->cell_amount[1]; k++)
        if ((geometry->domains_amount[0] * geometry->cell_amount[0] - geometry->cell_dims[0] - i) * geometry->cell_size[0]
            < lenght_sigma_extern)
          sigma.inc(i, k, geometry->pml_sigma[0] +
                    (geometry->pml_sigma[1] - geometry->pml_sigma[0])
                    / pow(lenght_sigma_extern, 2)
                    * pow(geometry->cell_size[0] * (i + 1) + radial_shift
                          - geometry->size[0] * geometry->domains_amount[0]
                          + lenght_sigma_extern, 2));

  // sigma on z=0 and z=z walls
  if (geometry->pml_length[1] != 0)
    for(int k = 0; k < geometry->cell_amount[1]; k++)
      for(int i = 0; i < geometry->cell_amount[0]; i++)
        // z=0 wall
        if (geometry->cell_size[1] * k + longitudinal_shift < lenght_sigma_left)
          sigma.inc(i, k, geometry->pml_sigma[0]
                    + (geometry->pml_sigma[1] - geometry->pml_sigma[0])
                    / pow(lenght_sigma_left, 2)
                    * pow(lenght_sigma_left - geometry->cell_size[1] * k - longitudinal_shift, 2));

  // z=z wall
  if (geometry->pml_length[3] != 0)

    for(int k = 0; k < geometry->cell_amount[1]; k++)
      for(int i = 0; i < geometry->cell_amount[0]; i++)
        if ((geometry->domains_amount[1] * geometry->cell_amount[1] - geometry->cell_dims[1] - k) * geometry->cell_size[1]
            < lenght_sigma_right)
          sigma.inc(i, k, geometry->pml_sigma[0]
                    + (geometry->pml_sigma[1] - geometry->pml_sigma[0])
                    / pow(lenght_sigma_right, 2)
                    * pow(geometry->cell_size[1] * (k + 1) + longitudinal_shift
                          - geometry->size[1] * geometry->domains_amount[1]
                          + lenght_sigma_right, 2));
}

void MaxwellSolverYee::calc_field_h()
{
  field_h.overlay_set(0);
  field_h_at_et.overlay_set(0);

  double dr = geometry->cell_size[0];
  double dz = geometry->cell_size[1];

  // H_r on outer wall (r=r)
  if (geometry->walls[2])
    for(int k = 0; k < geometry->cell_amount[1]; k++)
    {
      int i = geometry->cell_amount[0] - 1;
      // alpha constant and delta_t production (to optimize calculations)
      double alpha_t = time->step
        * (field_e(1, i, k + 1) - field_e(1, i, k)) / (dz * MAGN_CONST);

      field_h[0].set(i, k, field_h_at_et(0, i, k) + alpha_t / 2);
      field_h_at_et[0].inc(i, k, alpha_t);
    }

  // regular case
  for (int i = 0; i < geometry->cell_amount[0]; i++)
    for (int k = 0; k < geometry->cell_amount[1]; k++)
    {
      double alpha_t_r = time->step
        * (field_e(1, i, k + 1) - field_e(1, i, k)) / (dz * MAGN_CONST);

      field_h[0].set(i, k, field_h_at_et(0, i, k) + alpha_t_r / 2);
      field_h_at_et[0].inc(i, k, alpha_t_r);

      double alpha_t_phi = time->step
        * ((field_e(2, i+1, k) - field_e(2, i, k)) / dr
           - (field_e(0, i, k+1) - field_e(0, i, k)) / dz)
        / MAGN_CONST;

      field_h[1].set(i, k, field_h_at_et(1, i, k) + alpha_t_phi / 2);
      field_h_at_et[1].inc(i, k, alpha_t_phi);

      double alpha_t_z = time->step
        * (
          (field_e(1, i+1, k) + field_e(1, i, k)) / (2. * dr * (i + 0.5 + geometry->cell_dims[0]))
          + (field_e(1, i+1, k) - field_e(1, i, k)) / dr)
        / MAGN_CONST;

      field_h[2].set(i, k, field_h_at_et(2, i, k) - alpha_t_z / 2);
      field_h_at_et[2].dec(i, k, alpha_t_z);
    }
}

void MaxwellSolverYee::calc_field_e()
{
  field_e.overlay_set(0);

  Grid3D<double> curr = current->current;

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

      double koef_e = (epsilonx2 - sigma_t) / (epsilonx2 + sigma_t);
      double koef_h =  2 * time->step / (epsilonx2 + sigma_t);

      field_e[0].m_a(i, k, koef_e);
      field_e[0].dec(i, k, (curr(0, i, k)
                            + (field_h_at_et(1, i, k)
                               - field_h_at_et(1, i, k-1)) / dz) * koef_h);

      field_e[2].m_a(i, k, koef_e);
      field_e[2].dec(i, k, (curr(2, i, k)
                            - field_h_at_et(1, i, k) * 4. / dr) * koef_h);
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

      double koef_e = (epsilonx2 - sigma_t) / (epsilonx2 + sigma_t);
      double koef_h =  2 * time->step / (epsilonx2 + sigma_t);

      field_e[2].m_a(i, k, koef_e);
      field_e[2].dec(i, k, (curr(2, i, k)
                            - (field_h_at_et(1, i, k) - field_h_at_et(1, i - 1, k)) / dr
                            - (field_h_at_et(1, i, k) + field_h_at_et(1, i-1, k))
                            / (2. * dr * (i + geometry->cell_dims[0])))
                     * koef_h);
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

      field_e[0].m_a(i, k, koef_e);
      field_e[0].dec(i, k, (curr(0, i, k)
                            + (field_h_at_et(1, i, k)
                               - field_h_at_et(1, i, k-1)) / dz) * koef_h);

      field_e[1].m_a(i, k, koef_e);
      field_e[1].dec(i, k, (curr(1, i, k)
                            - (field_h_at_et(0, i, k) - field_h_at_et(0, i, k - 1))
                            / dz + (field_h_at_et(2, i, k)
                                    - field_h_at_et(2, i - 1, k)) / dr) * koef_h);

      field_e[2].m_a(i, k, koef_e);
      field_e[2].dec(i, k, (curr(2, i, k)
                            - (field_h_at_et(1, i, k) - field_h_at_et(1, i - 1, k)) / dr
                            - (field_h_at_et(1, i, k)
                               + field_h_at_et(1, i-1, k))
                            / (2. * dr * (i + geometry->cell_dims[0]))) * koef_h);

      if ( isnan(field_e[0](i,k)) || isnan(field_e[1](i,k)) || isnan(field_e[2](i,k)) )
        LOG_S(FATAL) << "fld " << i << " " << k << " "
                     << field_e[0](i,k) << " " << field_e[1](i,k) << " " << field_e[2](i,k);
    }
}

vector3d<double> MaxwellSolverYee::get_field_h(double radius, double longitude)
//! function for magnetic field weighting
{
  vector3d<double> cmp;
  int i_r = 0; // number of particle i cell
  int k_z = 0; // number of particle k cell
  int i_r_shift = 0;
  int k_z_shift = 0;

  double dr = geometry->cell_size[0];
  double dz = geometry->cell_size[1];
  double r1, r2, r3; // temp variables for calculation
  double dz1, dz2; // temp var.: width of k and k+1 cell
  double vol_1 = 0; // volume of i cell; Q/V, V - volume of elementary cell
  double vol_2 = 0; // volume of i+1 cell;

  r1 = radius - 0.5 * dr;
  r3 = radius + 0.5 * dr;

  //// weighting of H_z
  // finding number of cell. example dr=0.5, radius = 0.7, i_r =0;!!
  i_r = CELL_NUMBER(radius - 0.5 * dr, dr);
  k_z = CELL_NUMBER(longitude, dz);
  i_r_shift = i_r - geometry->cell_dims[0];
  k_z_shift = k_z - geometry->cell_dims[1];
  // // TODO: workaround: sometimes it gives -1.
  // // Just get 0 cell if it happence
  // if (i_r < 0) i_r = 0;
  // if (k_z < 0) k_z = 0;

  vol_1 = CELL_VOLUME(i_r+1, dr, dz);
  vol_2 = CELL_VOLUME(i_r+3, dr, dz);
  dz1 = (k_z + 1) * dz - longitude;
  dz2 = longitude - k_z * dz;
  r2 = (i_r + 1) * dr;

  //weighting Hz[i][k]//
  cmp[2] += field_h(2, i_r_shift, k_z_shift) * CYL_RNG_VOL(dz1, r1, r2) / vol_1;

  //weighting Hz[i+1][k]//
  cmp[2] += field_h(2, i_r_shift + 1, k_z_shift) * CYL_RNG_VOL(dz1, r2, r3) / vol_2;

  //weighting Hz[i][k+1]//
  cmp[2] += field_h(2, i_r_shift, k_z_shift + 1) * CYL_RNG_VOL(dz2, r1, r2) / vol_1;

  //weighting Hz[i+1][k+1]//
  cmp[2] += field_h(2, i_r_shift + 1, k_z_shift + 1) * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

  //// weighting of Hr
  // finding number of cell. example dz=0.5, longitude = 0.7, z_k =0;!!
  i_r = CELL_NUMBER(radius, dr);
  k_z = CELL_NUMBER(longitude - 0.5 * dz, dz);
  i_r_shift = i_r - geometry->cell_dims[0];
  k_z_shift = k_z - geometry->cell_dims[1];
  // // TODO: workaround: sometimes it gives -1.
  // // Just get 0 cell if it happence
  // if (i_r < 0) i_r = 0;
  // if (k_z < 0) k_z = 0;

  vol_1 = CELL_VOLUME(i_r, dr, dz);
  vol_2 = CELL_VOLUME(i_r+2, dr, dz);

  r2 = (i_r + 0.5) * dr;
  dz1 = (k_z + 1.5) * dz - longitude;
  dz2 = longitude - (k_z + 0.5) * dz;

  //weighting Hr[i][k]//
  cmp[0] += field_h(0, i_r_shift, k_z_shift) * CYL_RNG_VOL(dz1, r1, r2) / vol_1;

  //weighting Hr[i+1][k]//
  cmp[0] += field_h(0, i_r_shift + 1, k_z_shift) * CYL_RNG_VOL(dz1, r2, r3) / vol_2;

  //weighting Hr[i][k+1]//
  cmp[0] += field_h(0, i_r_shift, k_z_shift + 1) * CYL_RNG_VOL(dz2, r1, r2) / vol_1;

  //weighting Hr[i+1][k+1]//
  cmp[0] += field_h(0, i_r_shift + 1, k_z_shift + 1) * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

  //// weighting of H_fi
  // finding number of cell. example dz=0.5, longitude = 0.7, z_k =0;
  i_r = CELL_NUMBER(radius - 0.5 * dr, dr);
  k_z = CELL_NUMBER(longitude - 0.5 * dz, dz);
  i_r_shift = i_r - geometry->cell_dims[0];
  k_z_shift = k_z - geometry->cell_dims[1];
  // // TODO: workaround: sometimes it gives -1.
  // // Just get 0 cell if it happence
  // if (i_r < 0) i_r = 0;
  // if (k_z < 0) k_z = 0;

  r2 = (i_r+1) * dr;
  vol_1 = CELL_VOLUME(i_r + 1, dr, dz);
  vol_2 = CELL_VOLUME(i_r + 3, dr, dz);
  dz1 = (k_z+1.5) * dz - longitude;
  dz2 = longitude - (k_z+0.5) * dz;

  //weighting Hphi[i][k]//
  cmp[1] += field_h(1, i_r_shift, k_z_shift) * CYL_RNG_VOL(dz1, r1, r2) / vol_1;

  //weighting Hphi[i+1][k]//
  cmp[1] += field_h(1, i_r_shift + 1, k_z_shift) * CYL_RNG_VOL(dz1, r2, r3) / vol_2;

  //weighting Hphi[i][k+1]//
  cmp[1] += field_h(1, i_r_shift, k_z_shift + 1) * CYL_RNG_VOL(dz2, r1, r2) / vol_1;

  //weighting Hphi[i+1][k+1]//
  cmp[1] += field_h(1, i_r_shift + 1, k_z_shift + 1) * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

  return cmp;
}

vector3d<double> MaxwellSolverYee::get_field_e(double radius, double longitude)
//! function for electric field weighting
{
  vector3d<double> cmp;

  int i_r = 0; // number of particle i cell
  int k_z = 0; // number of particle k cell
  int i_r_shift = 0; // number of particle i cell
  int k_z_shift = 0; // number of particle k cell

  double dr = geometry->cell_size[0];
  double dz = geometry->cell_size[1];
  double r1, r2, r3; // temp variables for calculation
  double dz1, dz2; // temp var.: width of k and k+1 cell
  double vol_1 = 0; // volume of i cell; Q/V, V - volume of elementary cell
  double vol_2 = 0; // volume of i+1 cell;

  r1 = radius - 0.5 * dr;
  r3 = radius + 0.5 * dr;

  // weighting of E_r
  // finding number of cell. example dr=0.5, radius = 0.7, i_r =0;!!
  i_r = CELL_NUMBER(radius - 0.5 * dr, dr);
  k_z = CELL_NUMBER(longitude, dz);
  i_r_shift = i_r - geometry->cell_dims[0];
  k_z_shift = k_z - geometry->cell_dims[1];
  // // TODO: workaround: sometimes it gives -1.
  // // Just get 0 cell if it happence
  // if (i_r < 0) i_r = 0;
  // if (k_z < 0) k_z = 0;

  vol_1 = CELL_VOLUME(i_r+1, dr, dz);
  vol_2 = CELL_VOLUME(i_r+3, dr, dz);
  dz1 = (k_z+1) * dz - longitude;
  dz2 = longitude - k_z * dz;
  r2 = (i_r + 1) * dr;
  //weighting Er[i][k]//
  cmp[0] += field_e(0, i_r_shift, k_z_shift) * CYL_RNG_VOL(dz1, r1, r2) / vol_1;

  //weighting Er[i+1][k]//
  cmp[0] += field_e(0, i_r_shift + 1, k_z_shift) * CYL_RNG_VOL(dz1, r2, r3) / vol_2;

  //weighting Er[i][k+1]//
  cmp[0] += field_e(0, i_r_shift, k_z_shift + 1) * CYL_RNG_VOL(dz2, r1, r2) / vol_1;

  //weighting Er[i+1][k+1]//
  cmp[0] += field_e(0, i_r_shift + 1, k_z_shift + 1) * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

  // weighting of E_z
  // finding number of cell. example dz=0.5, longitude = 0.7, z_k =0;!!
  i_r = CELL_NUMBER(radius, dr);
  k_z = CELL_NUMBER(longitude - 0.5 * dz, dz);
  i_r_shift = i_r - geometry->cell_dims[0];
  k_z_shift = k_z - geometry->cell_dims[1];

  vol_1 = CELL_VOLUME(i_r, dr, dz);
  vol_2 = CELL_VOLUME(i_r+2, dr, dz);
  r2 = (i_r + 0.5) * dr;

  dz1 = (k_z + 1.5) * dz - longitude;
  dz2 = longitude - (k_z + 0.5) * dz;

  // weighting Ez[i][k]
  cmp[2] += field_e(2, i_r_shift, k_z_shift) * CYL_RNG_VOL(dz1, r1, r2) / vol_1;

  // weighting Ez[i+1][k]
  cmp[2] += field_e(2, i_r_shift+1, k_z_shift) * CYL_RNG_VOL(dz1, r2, r3) / vol_2;

  // weighting Ez[i][k+1]
  cmp[2] += field_e(2, i_r_shift, k_z_shift+1) * CYL_RNG_VOL(dz2, r1, r2) / vol_1;

  // weighting Ez[i+1][k+1]
  cmp[2] += field_e(2, i_r_shift+1, k_z_shift+1) * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

  // weighting of E_fi
  // finding number of cell. example dz=0.5, longitude = 0.7, z_k =1;
  i_r = CELL_NUMBER(radius, dr);
  k_z = CELL_NUMBER(longitude, dz);
  i_r_shift = i_r - geometry->cell_dims[0];
  k_z_shift = k_z - geometry->cell_dims[1];
  // // TODO: workaround: sometimes it gives -1.
  // // Just get 0 cell if it happence
  // if (i_r < 0) i_r = 0;
  // if (k_z < 0) k_z = 0;

  vol_1 = CELL_VOLUME(i_r, dr, dz);
  vol_2 = CELL_VOLUME(i_r+2, dr, dz);
  r2 = (i_r + 0.5) * dr;

  dz1 = (k_z+1)*dz-longitude;
  dz2 = longitude - k_z * dz;

  // weighting Efi[i][k]
  cmp[1] += field_e(1, i_r_shift, k_z_shift) * CYL_RNG_VOL(dz1, r1, r2) / vol_1;

  // weighting Efi[i+1][k]
  cmp[1] += field_e(1, i_r_shift+1, k_z_shift) * CYL_RNG_VOL(dz1, r2, r3) / vol_2;

  // weighting Efi[i][k+1]
  cmp[1] += field_e(1, i_r_shift, k_z_shift+1) * CYL_RNG_VOL(dz2, r1, r2) / vol_1;

  // weighting Efi[i+1][k+1]
  cmp[2] += field_e(1, i_r_shift+1, k_z_shift+1) * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

  return cmp;
}

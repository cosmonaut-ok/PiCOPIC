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

#ifndef _MAXWELLSOLVER_HPP_
#define _MAXWELLSOLVER_HPP_

#include "defines.hpp"
#include "msg.hpp"

#include "math/vector3d.hpp"
#include "algo/grid3d.hpp"

#include "geometry.hpp"

#include "current.hpp"

using namespace std;

class Current;

class MaxwellSolver
{

public:
  Grid3D<double> field_e;
  Grid3D<double> field_h;

protected:
  Geometry *geometry;
  TimeSim *time;

  Grid<double> epsilon;
  Grid<double> sigma; // for PML
  Current *current;

public:
  MaxwellSolver ( void ) {};
  MaxwellSolver ( Geometry *_geometry, TimeSim *_time,
                  Current* _current )
    : geometry(_geometry), time(_time), current(_current)
  {
    // initialize fields to zero-state
    field_e = Grid3D<double> (geometry->cell_amount[0], geometry->cell_amount[1], 2);
    field_h = Grid3D<double> (geometry->cell_amount[0], geometry->cell_amount[1], 2);
    field_e = 0.;
    field_h = 0.;
    field_e.overlay_set(0.);
    field_h.overlay_set(0.);

    // initialize epsilon and sigma (for PML)
    epsilon = Grid<double> (geometry->cell_amount[0], geometry->cell_amount[1], 2);
    epsilon = constant::EPSILON0;
    epsilon.overlay_set(constant::EPSILON0);

    sigma = Grid<double> (geometry->cell_amount[0], geometry->cell_amount[1], 2);
    sigma = 0.;
    sigma.overlay_set(0.);
  };

  ~MaxwellSolver(void) {};

  vector3d<double> get_field_dummy(__attribute__((unused)) double radius, __attribute__((unused)) double longitude)
  {
    // just a dummy method for test/debug
    vector3d<double> cmp; // field components
    return cmp;
  };

  virtual void set_pml() = 0;
  virtual void calc_field_h() = 0;
  virtual void calc_field_e() = 0;
  // virtual vector3d<double> get_field_h(double radius, double longitude) = 0;
  // virtual vector3d<double> get_field_e(double radius, double longitude) = 0;

  virtual vector3d<double> get_field_h(double radius, double longitude)
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
  };

  virtual vector3d<double> get_field_e(double radius, double longitude)
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
  };
};

#endif // end of _MAXWELLSOLVER_HPP_

#pragma once

#include <vector>
#include "grid3d.hpp"
#include "geometry.hpp"
#include "timeSim.hpp"
#include "math/vector3d.hpp"
class Field
{
public:
  Grid3D<double> field;
  Geometry *geometry;
  TimeSim *time;

// private:
// vector<SpecieP> *species_p;

public:
  Field(void) {};
  Field(Geometry *geom, TimeSim *t) : geometry(geom), time(t)
  {
    field = Grid3D<double> (geometry->r_grid_amount+1, geometry->z_grid_amount+1);
    time = t;
    field.setall(0);
  };

  ~Field(void) {};

  void set_homogenous_field(double r_value, double phi_value, double z_value)
  {
    field[0].setall(r_value);
    field[1].setall(phi_value);
    field[2].setall(z_value);
  };

  vector3d<double> get_field(double radius, double longitude)
  {
    vector3d<double> cmp(0, 0, 0); // field components

    int i_r, k_z, i_r_shift, k_z_shift; // absolute and related numbers of cell

    double dr = geometry->r_cell_size;
    double dz = geometry->z_cell_size;
    double r1, r2, r3; // temp variables for calculation
    double dz1, dz2; // temp var.: width of k and k+1 cell
    // double er = 0;
    // double efi = 0;
    // double ez = 0;
    double vol_1 = 0; // volume of i cell; Q/V, V - volume of elementary cell
    double vol_2 = 0; // volume of i+1 cell;

    r1 = radius - 0.5 * dr;
    r3 = radius + 0.5 * dr;

    // weighting of E_r
    // finding number of cell. example dr=0.5, radius = 0.7, i_r =0;!!
    i_r = CELL_NUMBER(radius - 0.5 * dr, dr); // TODO: why we shifting to cell/2?
    k_z = CELL_NUMBER(longitude, dz);
    i_r_shift = i_r - geometry->bottom_r_grid_number;
    k_z_shift = k_z - geometry->left_z_grid_number;
    // TODO: workaround: it gives -1, because of shifting to cell/2
    // Just get 0 cell if it happens
    if (i_r < 0) i_r = 0;
    if (i_r_shift < 0) i_r_shift = 0;
    if (k_z < 0) k_z = 0;
    if (k_z_shift < 0) k_z_shift = 0;

    vol_1 = CELL_VOLUME(i_r+1, dr, dz);
    vol_2 = CELL_VOLUME(i_r+3, dr, dz);
    dz1 = (k_z+1) * dz - longitude;
    dz2 = longitude - k_z * dz;
    r2 = (i_r + 1) * dr;

    //weighting Er[i][k]//
    cmp[0] += field[0].get(i_r_shift, k_z_shift) * CYL_RNG_VOL(dz1, r1, r2) / vol_1;

    //weighting Er[i+1][k]//
    cmp[0] += field[0].get(i_r_shift + 1, k_z_shift) * CYL_RNG_VOL(dz1, r2, r3) / vol_2;

    //weighting Er[i][k+1]//
    cmp[0] += field[0].get(i_r_shift, k_z_shift + 1) * CYL_RNG_VOL(dz2, r1, r2) / vol_1;

    //weighting Er[i+1][k+1]//
    cmp[0] += field[0].get(i_r_shift + 1, k_z_shift + 1) * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

    // weighting of E_z
    // finding number of cell. example dz=0.5, longitude = 0.7, z_k =0;!!
    i_r = CELL_NUMBER(radius, dr);
    k_z = CELL_NUMBER(longitude - 0.5 * dz, dz); // TODO: why we shifting to cell/2?
    i_r_shift = i_r - geometry->bottom_r_grid_number;
    k_z_shift = k_z - geometry->left_z_grid_number;
    // TODO: workaround: it gives -1, because of shifting to cell/2
    // Just get 0 cell if it happence
    if (i_r < 0) i_r = 0;
    if (i_r_shift < 0) i_r_shift = 0;
    if (k_z < 0) k_z = 0;
    if (k_z_shift < 0) k_z_shift = 0;

    if (radius > dr)
      vol_1 = CELL_VOLUME(i_r, dr, dz);
    else
      vol_1 = CYL_VOL(dz, dr); // volume of first cell

    r2 = (i_r + 0.5) * dr;
    vol_2 = CELL_VOLUME(i_r+2, dr, dz);
    dz1 = (k_z + 1.5) * dz - longitude;
    dz2 = longitude - (k_z + 0.5) * dz;

    // weighting Ez[i][k]
    cmp[2] += field[2].get(i_r_shift, k_z_shift) * CYL_RNG_VOL(dz1, r1, r2) / vol_1;

    // weighting Ez[i+1][k]
    cmp[2] += field[2].get(i_r_shift+1, k_z_shift) * CYL_RNG_VOL(dz1, r2, r3) / vol_2;

    // weighting Ez[i][k+1]
    cmp[2] += field[2].get(i_r_shift, k_z_shift+1) * CYL_RNG_VOL(dz2, r1, r2) / vol_1;

    //weighting Ez[i+1][k+1]//
    cmp[2] += field[2].get(i_r_shift+1, k_z_shift+1) * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

    // weighting of E_fi
    // finding number of cell. example dz=0.5, longitude = 0.7, z_k =1;
    i_r = CELL_NUMBER(radius, dr);
    k_z = CELL_NUMBER(longitude, dz);
    i_r_shift = i_r - geometry->bottom_r_grid_number;
    k_z_shift = k_z - geometry->left_z_grid_number;

    if (i_r < 0) i_r = 0;
    if (i_r_shift < 0) i_r_shift = 0;
    if (k_z < 0) k_z = 0;
    if (k_z_shift < 0) k_z_shift = 0;

    if(radius>dr)
      vol_1 = CELL_VOLUME(i_r, dr, dz);
    else
      vol_1 = CYL_VOL(dz, dr); // volume of first cell

    r2 = (i_r + 0.5) * dr;
    vol_2 = CELL_VOLUME(i_r+2, dr, dz);
    dz1 = (k_z+1)*dz-longitude;
    dz2 = longitude - k_z * dz;

    // weighting Efi[i][k]
    cmp[1] += field[1].get(i_r_shift, k_z_shift) * CYL_RNG_VOL(dz1, r1, r2) / vol_1;

    // weighting Efi[i+1][k]
    cmp[1] += field[1].get(i_r_shift+1, k_z_shift) * CYL_RNG_VOL(dz1, r2, r3) / vol_2;

    // weighting Efi[i][k+1]
    cmp[1] += field[1].get(i_r_shift, k_z_shift+1) * CYL_RNG_VOL(dz2, r1, r2) / vol_1;

    // weighting Efi[i+1][k+1]
    cmp[2] += field[1].get(i_r_shift+1, k_z_shift+1) * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

    return cmp;
  };

  vector3d<double> get_field_dummy(__attribute__((unused)) double radius, __attribute__((unused)) double longitude)
  {
    // just dummy method for test/debug
    vector3d<double> cmp(0, 0, 0); // field components
    return cmp;
  }

};

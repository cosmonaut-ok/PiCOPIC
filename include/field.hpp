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

protected:
  // vector<SpecieP> *species_p;
  unsigned int r_begin;
  unsigned int z_begin;
  unsigned int r_end;
  unsigned int z_end;

public:
  Field(void) {};
  Field(Geometry *geom, TimeSim *t) : geometry(geom), time(t)
  {
    field = Grid3D<double> (geometry->r_grid_amount, geometry->z_grid_amount);
    time = t;
    field = 0;

    r_begin = 1; // FIXME: 0 gives segfault with -nan, -nan velocity
    z_begin = 0;
    r_end = geometry->r_grid_amount;
    z_end = geometry->z_grid_amount;

    // emulate dielectric walls
    if (geometry->walls[0]) // r=0
      r_begin = 1;

    if (geometry->walls[1]) // z=0
      z_begin = 1;

    if (geometry->walls[2]) // r=r
      r_end = geometry->r_grid_amount - 1;

    if (geometry->walls[3]) // z=z
      z_end = geometry->z_grid_amount - 1;
  };

  ~Field(void) {};

  void set_homogenous_field(double r_value, double phi_value, double z_value)
  {
    field[0] = r_value;
    field[1] = phi_value;
    field[2] = z_value;
  };

  vector3d<double> get_field_dummy(__attribute__((unused)) double radius, __attribute__((unused)) double longitude)
  {
    // just dummy method for test/debug
    vector3d<double> cmp(0, 0, 0); // field components
    return cmp;
  };
};

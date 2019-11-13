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
    //! The important feature is to leave so called 'overlay areas'
    //! around field and current grids. In this case it takes
    //! one grid cell before and after main grid area.
    //! e.g. overlay areas are grid[0:1][i], grid[i][0:1],
    //! grid[0:1][geometry->z_grid_amount-1:geometry->z_grid_amount-2]
    //! and grid[geometry->r_grid_amount-1:geometry->r_grid_amount-2][i]
    //! So, we need to shift all of the grid numbers to 2
    //! and count it from 1 to geometry->(r,z)_grid_amount
    //! instead of 0 to geometry->(r,z)_grid_amount

    field = Grid3D<double> (geometry->r_grid_amount + 2, geometry->z_grid_amount + 2);
    time = t;
    field = 0;
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

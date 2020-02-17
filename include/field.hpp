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
    //! The important feature is to leave so called 'overlay areas'
    //! around field and current grids. In this case it takes
    //! one grid cell before and after main grid area.
    //! e.g. overlay areas are grid[0:1][i], grid[i][0:1],
    //! grid[0:1][geometry->z_grid_amount-1:geometry->z_grid_amount-2]
    //! and grid[geometry->r_grid_amount-1:geometry->r_grid_amount-2][i]
    //! So, we need to shift all of the grid numbers to 2
    //! and count it from 1 to geometry->(r,z)_grid_amount
    //! instead of 0 to geometry->(r,z)_grid_amount

    field = Grid3D<double> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
    time = t;
    field = 0;
    field.overlay_set(0);
    
    // shift to use overlay area (2 before and 2 after main grid)
    r_begin = 0;
    z_begin = 0;
    r_end = geometry->r_grid_amount;
    z_end = geometry->z_grid_amount;
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

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

#include "constant.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include "lib.hpp"
#include "specieP.hpp"
#include "density.hpp"

class Temperature
{
public:
  Geometry *geometry;
  Grid<double> vel_full;
  Grid<double> vel_r;
  Grid<double> vel_phi;
  Grid<double> vel_z;
  Grid<double> temperature;
  vector<SpecieP *> species_p;

  Temperature(void) {};
  Temperature(Geometry *geom, vector<SpecieP *> species) : geometry(geom)
  {
    vel_full = Grid<double> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
    vel_r = Grid<double> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
    vel_phi = Grid<double> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
    vel_z = Grid<double> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
    temperature = Grid<double> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
    species_p = species;
  };
  
  ~Temperature(void) {};

  virtual void calc_temperature_cylindrical(string specie) = 0;
};

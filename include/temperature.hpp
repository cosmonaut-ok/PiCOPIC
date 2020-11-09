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

#ifndef _TEMPERATURE_HPP_
#define _TEMPERATURE_HPP_

#include <vector>

#include "constant.hpp"
#include "geometry.hpp"
#include "algo/grid.hpp"
#include "algo/common.hpp"
#include "specieP.hpp"
#include "density.hpp"

class Temperature
{
public:
  Geometry *geometry;
  Grid<double> p_abs;
  Grid<double> p_r;
  Grid<double> p_phi;
  Grid<double> p_z;
  Grid<double> tmpr;
  vector<SpecieP *> species_p;

  Temperature(void) {};
  Temperature(Geometry *geom, vector<SpecieP *> species) : geometry(geom)
  {
    p_abs = Grid<double> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
    p_r = Grid<double> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
    p_phi = Grid<double> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
    p_z = Grid<double> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
    tmpr = Grid<double> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
    species_p = species;
  };

  ~Temperature(void) {};

  virtual void operator()(string specie) = 0;
};

#endif // end of _TEMPERATURE_HPP_

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

#ifndef _DENSITY_HPP_
#define _DENSITY_HPP_

#include "defines.hpp"
#include "msg.hpp"
#include "algo/weighter.hpp"
#include "specieP.hpp"

class Density
{
public:
  Geometry *geometry;
  Grid<double> density;
  vector<SpecieP *> species_p;

  Density(void) {};
  Density(Geometry *geom, vector<SpecieP *> species) : geometry(geom)
  {
    density = Grid<double> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
    species_p = species;

    density = 0;
    density.overlay_set(0);
  };
  ~Density() {};

  virtual void operator()(string specie) = 0;
};

#endif // end of _DENSITY_HPP_

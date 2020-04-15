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

#ifndef _CURRENT_HPP_
#define _CURRENT_HPP_

#include <vector>

#include "geometry.hpp"
#include "grid3d.hpp"
#include "timeSim.hpp"
#include "specieP.hpp"

class SpecieP;
class Geometry;

class Current
{
public:
  TimeSim *time;
  Geometry *geometry;
  Grid3D<double> current;
  vector<SpecieP *> species_p;

  Current() {};
  Current(Geometry *geom, TimeSim *t, vector<SpecieP *> species) : geometry(geom), time(t)
  {
    current = Grid3D<double> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
    species_p = species;

    current = 0;
    current.overlay_set(0);
  };

  virtual void current_distribution() = 0;
};

#endif // end of _CURRENT_HPP_

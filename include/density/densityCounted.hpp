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

#ifndef _DENSITY_COUNTED_HPP_
#define _DENSITY_COUNTED_HPP_

#include "geometry.hpp"
#include "algo/grid.hpp"
#include "algo/weighter.hpp"
#include "density.hpp"
#include "specieP.hpp"

class DensityCounted : public Density
{
public:
  DensityCounted () {};

  DensityCounted (Geometry *geom, vector<SpecieP *> species) : Density(geom, species) {};
  ~DensityCounted () {};

  void calc_density_cylindrical(string specie);
};

#endif // end of _DENSITY_COUNTED_HPP_

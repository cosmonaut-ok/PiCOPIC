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

#ifndef _DENSITY_HPP_
#define _DENSITY_HPP_

#include "geometry.hpp"
#include "grid.hpp"
#include "weighter.hpp"
#include "specieP.hpp"

class Density
{
public:
  Geometry *geometry;
  Grid<double> density;
  vector<SpecieP *> species_p;
  
  Density(void) {};
  Density(Geometry *geom, vector<SpecieP *> species);
  ~Density() {};
  
  void calc_density_cylindrical(string specie);
};

#endif // end of _DENSITY_HPP_

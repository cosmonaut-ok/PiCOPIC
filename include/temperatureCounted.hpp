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

#ifndef _TEMPERATURE_COUNTED_HPP_
#define _TEMPERATURE_COUNTED_HPP_

#include "constant.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include "specieP.hpp"
#include "temperature.hpp"

class TemperatureCounted : public Temperature
{
public:
  Grid<double> count;

  TemperatureCounted () {};
  TemperatureCounted (Geometry *geom, vector<SpecieP *> species) : Temperature(geom, species)
  {
    count = Grid<double> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
  };
  ~TemperatureCounted () {};

  void calc_temperature_cylindrical(string specie);

private:
  void weight_temperature_cylindrical(string specie);
};
#endif // end of _TEMPERATURE_COUNTED_HPP_

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

#ifndef _TEMPERATURE_WEIGHTED_HPP_
#define _TEMPERATURE_WEIGHTED_HPP_

#include "constant.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include "weighter.hpp"
#include "specieP.hpp"
#include "densityWeighted.hpp"
#include "temperature.hpp"

class TemperatureWeighted : public Temperature
{
public:
  DensityWeighted density;

  TemperatureWeighted (void) {};
  TemperatureWeighted (Geometry *geom, vector<SpecieP *> species) : Temperature(geom, species)
  {
    density = DensityWeighted (geometry, species);
  };
  ~TemperatureWeighted () {};

  void calc_temperature_cylindrical(string specie);

private:
  void weight_temperature_cylindrical(string specie);
};
#endif // end of _TEMPERATURE_WEIGHTED_HPP_

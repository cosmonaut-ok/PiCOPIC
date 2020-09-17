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

#ifndef _COLLISIONS_HPP_
#define _COLLISIONS_HPP_

#include <vector>
#include <algorithm>    // std::min, std::random_shuffle
#include <typeinfo>
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <math.h>       // floor, asin

#include <string>

#include "defines.hpp"
#include "loguru.hpp"
#include "lib.hpp"
#include "geometry.hpp"
#include "specieP.hpp"

#include "temperatureCounted.hpp"
#include "densityCounted.hpp"

const double PROTON_MASS = EL_MASS * 1836;

class Collisions
{
public:

  TemperatureCounted temperature_el;
  TemperatureCounted temperature_ion;

  Collisions(void) {};
  Collisions(Geometry* _geometry, TimeSim *_time, vector <SpecieP *> _species_p);
  ~Collisions(void) {};

  void sort_to_cells();
  void random_sort();

  virtual void clear ();
  virtual void run ();

protected:
  vector <SpecieP *> species_p;
  Grid < vector< vector<double> * > > map_el2cell;
  Grid < vector< vector<double> * > > map_ion2cell;

  Grid < double > energy_tot_el;
  Grid < double > mass_tot_el;
  Grid3D < double > moment_tot_el;

  Grid < double > energy_tot_ion;
  Grid < double > mass_tot_ion;
  Grid3D < double > moment_tot_ion;

  Geometry* geometry;
  TimeSim *time;

  double get_el_density(int i, int j);
  double get_ion_density(int i, int j);
  double get_el_temperature(int i, int j);
  double get_ion_temperature(int i, int j);

  double geometry_cell_volume(int i);

  virtual void collect_weighted_params_tot_grid();
  virtual void collide () = 0;

};

#endif // end of _COLLISIONS_HPP_

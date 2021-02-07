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

#ifndef _COLLISIONS_HPP_
#define _COLLISIONS_HPP_

#include <vector>
#include <algorithm>    // std::min, std::random_shuffle

#include "defines.hpp"
#include "msg.hpp"

#include "algo/grid3d.hpp"
#include "phys/plasma.hpp"

#include "specieP.hpp"

class Collisions
{
public:
  double mass_el;
  double mass_ion;

  double charge_el;
  double charge_ion;

  Collisions(void) {};
  Collisions(Geometry* _geometry, TimeSim *_time, vector <SpecieP *> _species_p);
  ~Collisions(void) {};

  void sort_to_cells();
  void random_sort();

  virtual void clear ();
  virtual void operator()() = 0;

protected:
  vector <SpecieP *> species_p;
  SpecieP* specie_el;
  SpecieP* specie_ion;
  Grid < vector< Particle * > > map_el2cell;
  Grid < vector< Particle * > > map_ion2cell;

  Grid < double > energy_tot_el;
  Grid < double > amount_tot_el;
  Grid3D < double > moment_tot_el;

  Grid < double > energy_tot_ion;
  Grid < double > amount_tot_ion;
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

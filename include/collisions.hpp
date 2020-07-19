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

#include "specieP.hpp"

class Collisions
{
public:
  Collisions(void) {};
  Collisions(Geometry* _geometry, TimeSim *_time, vector <SpecieP *> _species_p);
  ~Collisions(void) {};

  // virtual void calc_collisions() = 0;
  void sort_to_cells();
  void clear();
  void random_sort();
  void collide_single(int i, int j, vector<double> &p1, vector<double> &p2);
  void collide();

protected:
  vector <SpecieP *> species_p;
  Grid < vector< vector<double> * > > map_el2cell;
  Grid < vector< vector<double> * > > map_ion2cell;
  Geometry* geometry;
  TimeSim *time;

  double get_el_density(int i, int j);
  double get_ion_density(int i, int j);
  double geometry_cell_volume(int i);

};

#endif // end of _COLLISIONS_HPP_

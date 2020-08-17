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

#ifndef _COLLISIONS_TANAKA_HPP_
#define _COLLISIONS_TANAKA_HPP_

#include <vector>
#include <algorithm>    // std::min, std::random_shuffle
#include <typeinfo>
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <math.h>       // floor, asin

#include <string>

#include "lib.hpp"
#include "phys/plasma.hpp"
#include "geometry.hpp"
#include "specieP.hpp"

#include "specieP.hpp"

#include "collisions.hpp"

class Collisions;

class CollisionsTanaka : public Collisions
{
public:
  CollisionsTanaka(void) {};
  CollisionsTanaka(Geometry* _geometry, TimeSim *_time, vector <SpecieP *> _species_p);
  ~CollisionsTanaka(void) {};

  // virtual void calc_collisions() = 0;
  void collide_single(int i, int j, double m_real_a, double m_real_b,
                      vector<double> &p1, vector<double> &p2);

  void run ();

protected:
  Grid < double > energy_tot_el;
  Grid < double > mass_tot_el;
  Grid3D < double > moment_tot_el;

  Grid < double > energy_tot_ion;
  Grid < double > mass_tot_ion;
  Grid3D < double > moment_tot_ion;

  void correct_velocities();

  void collide();
};

#endif // end of _COLLISIONS_TANAKA_HPP_

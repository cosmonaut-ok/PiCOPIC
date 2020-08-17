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

#ifndef _COLLISIONS_SENTOKU_M_HPP_
#define _COLLISIONS_SENTOKU_M_HPP_

#include <vector>
#include <algorithm>    // std::min, std::random_shuffle
#include <typeinfo>
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <math.h>       // floor, asin

#include <string>

#include "lib.hpp"
#include "geometry.hpp"
#include "specieP.hpp"

#include "specieP.hpp"

#include "collisions.hpp"

class Collisions;

class CollisionsSentokuM : public Collisions
{
public:
  CollisionsSentokuM(void) {};
  CollisionsSentokuM(Geometry* _geometry, TimeSim *_time, vector <SpecieP *> _species_p);
  ~CollisionsSentokuM(void) {};

  void collide_single(int i, int j, double m_real_a, double m_real_b,
                      vector<double> &p1, vector<double> &p2);

  void run ();

protected:
  void collide();
  double get_rel_p0 (double mass, double sq_vel);

double get_coulomb_logarithm (double m_a, double m_b,
                              double debye_length,
                              double v_rel);

double get_collision_freq (double e_a, double e_b,
                           double density_lowest,
                           double L,
                           double p_rel, double v_rel);
};

#endif // end of _COLLISIONS_SENTOKU_M_HPP_

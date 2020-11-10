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

#ifndef _COLLISIONS_SK98_HPP_
#define _COLLISIONS_SK98_HPP_

#include "collisions.hpp"

class CollisionsSK98 : public Collisions
{
public:
  CollisionsSK98(void) {};
  CollisionsSK98(Geometry* _geometry, TimeSim *_time, vector <SpecieP *> _species_p);
  ~CollisionsSK98(void) {};

  void collide_single(double m_real_a, double m_real_b,
		      double q_real_a, double q_real_b,
                      vector<double> &p1, vector<double> &p2,
                      double _density_a, double _density_b, double debye);

  void operator()();

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

#endif // end of _COLLISIONS_SK98_HPP_

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

#ifndef _CURRENT_VB_HPP_
#define _CURRENT_VB_HPP_

#include <vector>

#include "geometry.hpp"
#include "grid3d.hpp"
#include "timeSim.hpp"
#include "specieP.hpp"
#include "current.hpp"

#define SOME_SHIT_DENSITY_STRICT(q, r, dr, dz, delta_t) \
  (q) / (PI * 4. * (r) * (dz) * (dz) * (dr) * (delta_t));

#define SOME_SHIT_DENSITY_R(q, r, dr, dz, delta_t) \
  (q) / (2. * PI * (r) * (dr) * (dz)                \
         * (delta_t) * 2. * (dr));

#define SOME_SHIT_DENSITY_Z(q, r, dr, dz, delta_t) \
  (q) / (2. * PI * (r) * (dr) * (dz)                \
         * (delta_t) * (dz));

using namespace std;

class SpecieP;
class Geometry;

class CurrentVB : public Current
{
public:
  CurrentVB() {};
  CurrentVB(Geometry *geom, TimeSim *t, vector<SpecieP *> species) : Current(geom, t, species) {};

  void current_distribution();

private:
  void rz_current_distribution();
  void azimuthal_current_distribution();

  void simple_current_distribution (double radius_new,
                                    double longitude_new,
                                    double radius_old,
                                    double longitude_old,
                                    int i_n,
                                    int k_n,
                                    double p_charge);

  void strict_motion_distribution (double radius_new,
                                   double longitude_new,
                                   double radius_old,
                                   double longitude_old,
                                   double p_charge);
};

#endif // end of _CURRENT_VB_HPP_

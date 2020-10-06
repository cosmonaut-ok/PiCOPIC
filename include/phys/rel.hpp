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

#ifndef _REL_HPP_
#define _REL_HPP_

#include <math.h>
#include <algorithm>

#include "defines.hpp"
#include "constant.hpp"
#include "loguru.hpp"
#include "lib.hpp"
#include "math/vector3d.hpp"

namespace phys::rel
{
  double lorenz_factor (double sq_velocity);
  double lorenz_factor_inv (double sq_velocity);
  double energy (double mass, double velocity_2); // mass and velocity powered to 2
  double energy_m (double mass, double momentum_2);  // mass and momentum powered to 2

  double momentum_0 (double mass, vector3d<double> velocity);
  vector3d<double> momentum (double mass, vector3d<double> velocity);
}
#endif // end of _REL_HPP_

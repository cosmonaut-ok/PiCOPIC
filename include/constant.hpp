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

#ifndef _CONSTANT_HPP_
#define _CONSTANT_HPP_

#include <ctime>
#include "defines.hpp"

// WARNING!
// This header file used only to store
// physical and mathematical constants

namespace constant
{
  // some math constants and common numbers
  const double ONE_OVER_3 = 1. / 3.;
  const double TWO_OVER_3 = 2. / 3.;
  const double FOUR_OVER_3 = 4. / 3.;
  const double THREE_OVER_4 = 4. / 3.;
  // PI constant
  const double PI = 3.1415926535897932;
  // Vacuum permittivity (electric constant), F*m(e-1)
  const double EPSILON0 = 8.85E-12;
  // Electron mass, kg
  const double EL_MASS = 9.1E-31;
  // Electron charge, coulon
  const double EL_CHARGE = 1.6E-19;
  const double EL_CHARGE_INV = 1. / EL_CHARGE;
  // light speed in vacuum m/s
  const double LIGHT_VEL = 3.0E8;
  // define C^2 to decrease number of operations
  const double LIGHT_VEL_POW_2 = LIGHT_VEL * LIGHT_VEL;
  // Vacuum permeability (magnetic constant), m*kg*s(e-2)*A(e-2)
  const double MAGN_CONST = 1.26E-6;
  // Plank constant
  const double PLANK_CONST = 6.6261e-34;
  const double PLANK_BAR_CONST = PLANK_CONST / (2. * PI);
  // simulation start time
  const clock_t SIMULATION_START_TIME = std::time(nullptr);
  // minimal possible distance or velocity.
  // Smaller values should be rounded to zero
  const double MNZL = 1e-15; // Minimal Non-Zeroing Limit
}

// define some service constants
// use classical calculations, if velocity lower, than minimal
#ifdef REL_LIMIT
// define REL_LIMIT^2 to decrease number of operations
const double REL_LIMIT_POW_2 = REL_LIMIT * REL_LIMIT;
#endif

#endif // end of _CONSTANT_HPP_

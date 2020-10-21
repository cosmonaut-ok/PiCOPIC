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

#ifndef _MAXWELL_JUETTNER_HPP_
#define _MAXWELL_JUETTNER_HPP_

#include <cmath>
#include <ctime>
#include <cstdlib>

#include <iostream>

// IDRIS
#include <cstring>
// IDRIS

#include <vector>
#include "constant.hpp"
#include "math/rand.hpp"

namespace math::maxwell_juettner
{
  vector<double> maxwellJuettner(unsigned int npoints, double temperature);
}

#endif // end of _MAXWELL_JUETTNER_HPP_

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

#ifndef _RAND_HPP_
#define _RAND_HPP_

#include <random>

#include "defines.hpp"
#include "constant.hpp"
#include "msg.hpp"

using namespace std;

namespace math::random
{
  double uniform();
  double uniform1();
  double uniform2();
  double uniform_angle();
  double normal(double stddev);
  double random_reverse(double vel, int power);
}

#endif // end of _RAND_HPP_

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

#include "math/rand.hpp"

using namespace std;

namespace math::random
{
  std::random_device device;
  std::mt19937 gen(device());

  std::uniform_real_distribution<double> uniform_distribution(0., 1.);
  double uniform()
  {
    return uniform_distribution(gen);
  }

  std::uniform_real_distribution<double> uniform_distribution1(0., 1.-1e-11);
  double uniform1()
  {
    return uniform_distribution1(gen);
  }

  std::uniform_real_distribution<double> uniform_distribution2(-1., 1.);
  double uniform2()
  {
    return uniform_distribution2(gen);
  }

  double normal(double stddev)
  {
    std::normal_distribution<double> normal_distribution(0., stddev);
    return normal_distribution(gen);
  }

  // random_reverse from PDP2: pseudorandom number generator
  // could be useful for debug to get the same values in same
  // cases
  double random_reverse(double vel, int power)
  {
    int int_vel =(int) floor(vel);
    double ost = 0;
    double r = 0;
    int order = 1;
    while(int_vel >= 1)
    {
      ost = int_vel % power;

      r = r + ost * pow((double)power, (-order));

      int_vel = (int_vel - ost)/power;
      order = order+1;
    }
    return r;
  }
}

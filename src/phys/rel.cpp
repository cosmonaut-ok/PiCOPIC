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

#include "phys/rel.hpp"

using namespace constant;

namespace phys::rel
{
  double lorenz_factor (double sq_velocity)
  {
    //! get_gamma takes squared velocity,
    //! only because of some features of code
    //! and optimisation issues
    double gamma, beta;

    beta = sq_velocity / LIGHT_VEL_POW_2;

    if (beta > 1) // it's VERY BAD! Beta should not be more, than 1
      LOG_S(FATAL) << "(lorenz_factor): Lorentz factor aka gamma is complex. Velocity is: " << lib::sq_rt(sq_velocity);

    gamma = 1 / lib::sq_rt(1.0 - beta);

    if (isinf(gamma) == 1)
      { // avoid infinity values
	LOG_S(WARNING) << "(lorenz_factor): gamma (Lorenz factor) girects to infinity for velocity" << lib::sq_rt(sq_velocity);
	return 1e100; // just return some very big value
      }
    else
      return gamma;
  }

  double lorenz_factor_inv (double sq_velocity)
  {
    //! get_gamma_inv takes squared velocity,
    //! only because of some features of code
    //! and optimisation issues
    double gamma = pow(1.0 + sq_velocity / LIGHT_VEL_POW_2, -0.5);

    return gamma;
  }

  double momentum_0 (double mass, vector3d<double> velocity)
  // get 0th momentum in 4-momentum space, aka E/c
  {
    if (mass == 0) LOG_S(FATAL) << "mass must not be zero";
    double gamma = lorenz_factor(velocity.length2());

    return gamma * mass * LIGHT_VEL; // gamma m c^2 / c
  }

  vector3d<double> momentum (double mass, vector3d<double> velocity)
  {
    if (mass == 0) LOG_S(FATAL) << "mass must not be zero";
    double gamma = lorenz_factor(velocity.length2());

    return velocity * gamma * mass; // gamma m v
  }

  double energy (double mass, double velocity_2)
  {
    double e;
    if (velocity_2 > REL_LIMIT_POW_2)
    {
      double gamma = lorenz_factor(velocity_2);
      e = ( gamma - 1 ) * mass * LIGHT_VEL_POW_2;
    }
    else
      e = 0.5 * mass * velocity_2;

    return e;
  }

  double energy_m (double mass, double momentum_2)
  // ! get energy from momentum
  {
    double e;
    if (momentum_2 > REL_LIMIT_POW_2 * mass * mass)
    {
      double e_rest = mass * LIGHT_VEL_POW_2;
      e = lib::sq_rt(momentum_2 * LIGHT_VEL_POW_2 + pow(e_rest, 2)) - e_rest;
    }
    else
      e = 0.5 * momentum_2 / mass;

    return e;
  }
}

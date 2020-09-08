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

#ifndef _PLASMA_HPP_
#define _PLASMA_HPP_

#include <math.h>
#include <algorithm>

#include "defines.hpp"
#include "constant.hpp"
#include "loguru.hpp"
#include "lib.hpp"
#include "math/vector3d.hpp"
#include "phys/rel.hpp"

namespace phys::plasma
{
  double debye_length (double density_el, double density_ion,
                       double temperature_ion, double temperature_el);

  double plasma_frequency (double density);

  double coulomb_logarithm (double mass_a, double mass_b,
			    double debye_length, double v_rel);

  double collision_freqency (double e_a, double e_b,
                             double density_lowest,
                             double L,
                             double p_rel, double v_rel);
}
#endif // end of _PLASMA_HPP_

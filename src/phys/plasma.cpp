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

#include "phys/plasma.hpp"

using namespace constant;

namespace phys::plasma
{
  double debye_length (double density, double temperature)
  {
    // sqrt ( epsilon_0 * k_boltzmann / q_el^2 ) is 7400,
    // when density in m-3 and temperature in eV
    const double const1 = 7400;
    if (density <= 0) LOG_S(FATAL) << "(debye_length): density must be positive. Actucal value is: ``" << density << "''";
    return ( const1 * lib::sq_rt ( temperature / density ) );
  }

  double plasma_frequency (double density)
  {
    // const1 is sqrt ( q_e^2 / ( m_e * epsilon_0 ) )
    const double const1 = 56.3803;
    return ( const1 * lib::sq_rt ( density ) );
  }

  double coulomb_logarithm (double m_a, double m_b, double debye_length, double v_rel)
  {
    // get reduced mass
    double m_ab = m_a * m_b / (m_a + m_b);

    // FIXME: temporary simple formula
    return log ( debye_length * m_ab * v_rel / PLANK_BAR_CONST );
  }

  double collision_freqency (double e_a, double e_b,
                             double density,
                             double L,
                             double p_rel, double v_rel)
  {
    return 4 * constant::PI * pow(e_a * e_b, 2) * density * L / (p_rel * p_rel * v_rel);
  }


}

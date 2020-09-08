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
  double debye_length (double density_el, double density_ion,
                       double temperature_el, double temperature_ion)
  {
    // sqrt ( epsilon_0 * k_boltzmann * T(eV)->T(K) / q_el^2 ) is 7433.944,
    // when density in m^-3 and temperature in eV
    const double const1 = 7433.944;
    if (density_el <= 0) LOG_S(FATAL) << "(debye_length): electron density must be positive. Actucal value is: ``" << density_el << "''";
    if (density_ion <= 0) LOG_S(FATAL) << "(debye_length): ion density must be positive. Actucal value is: ``" << density_ion << "''";
    return ( const1 / lib::sq_rt ( density_el / temperature_el + density_ion / temperature_ion ) );
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
    return pow(e_a * e_b, 2) * density * L
      / (8 * constant::PI
         * pow(constant::EPSILON0, 2)
         * pow(p_rel, 2) * v_rel, 2);
  }
}

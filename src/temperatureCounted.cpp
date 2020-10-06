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

#include "temperatureCounted.hpp"

void TemperatureCounted::weight_temperature_cylindrical(string specie)
{
  for (auto ps = species_p.begin(); ps != species_p.end(); ++ps)
    if (specie.compare((**ps).name) == 0)
      for (auto i = (**ps).particles.begin(); i != (**ps).particles.end(); ++i)
      {
        // finding number of i and k cell. example: dr = 0.5; r = 0.4; i =0
        unsigned int r_i = P_CELL_R((**i));
        unsigned int z_k = P_CELL_Z((**i));
        unsigned int r_i_shift = r_i - geometry->bottom_r_grid_number;
        unsigned int z_k_shift = z_k - geometry->left_z_grid_number;

        double vel_r_single = P_VEL_R((**i));
        double vel_phi_single = P_VEL_PHI((**i));
        double vel_z_single = P_VEL_Z((**i));

        double vel_abs_single = lib::sq_rt(pow(vel_r_single, 2)
                                           + pow(vel_phi_single, 2)
                                           + pow(vel_z_single, 2));

        double m_weighted = (**ps).mass * P_WEIGHT((**i));

        double p_r_single = m_weighted * vel_r_single;
        double p_phi_single = m_weighted * vel_phi_single;
        double p_z_single = m_weighted * vel_z_single;
        double p_abs_single = m_weighted * vel_abs_single;

        if (vel_abs_single < REL_LIMIT)
        {
          double gamma = phys::rel::lorenz_factor(vel_abs_single * vel_abs_single);
          p_r *= gamma;
          p_phi *= gamma;
          p_z *= gamma;
          p_abs *= gamma;
        }

        p_r.inc(r_i_shift, z_k_shift, p_r_single);
        p_phi.inc(r_i_shift, z_k_shift, p_phi_single);
        p_z.inc(r_i_shift, z_k_shift, p_z_single);
        p_abs.inc(r_i_shift, z_k_shift, p_abs_single);
        count.inc(r_i_shift, z_k_shift, P_WEIGHT((**i)));
      }
}

void TemperatureCounted::calc_temperature_cylindrical(string specie)
{
  // get specie
  SpecieP *speciep;

  for (auto ps = species_p.begin(); ps != species_p.end(); ++ps)
    if (specie.compare((**ps).name) == 0)
      speciep = (*ps);

  // clear temperature grid
  tmpr = 0;
  tmpr.overlay_set(0);

  // clear velociry components grid
  p_r = 0;
  p_r.overlay_set(0);
  p_phi = 0;
  p_phi.overlay_set(0);
  p_z = 0;
  p_z.overlay_set(0);

  // clear full velocity grid
  p_abs = 0;
  p_abs.overlay_set(0);

  // clear count
  count = 0;
  count.overlay_set(0);

  weight_temperature_cylindrical(specie);

  for (int r = 0; r < geometry->r_grid_amount; r++)
    for (int z = 0; z < geometry->z_grid_amount; z++)
    {
      double p_vec_sum_2 =
        pow(p_r(r, z), 2)
        + pow(p_phi(r, z), 2)
        + pow(p_z(r, z), 2);

      double p_sc_sum_2 = pow(p_abs(r, z), 2) - p_vec_sum_2;

      // normalize to count and convert Joules to eV
      p_sc_sum_2 /= pow(count(r, z), 2);

      double energy = phys::rel::energy_m(speciep->mass, p_sc_sum_2);

      // convert Joules to eV
      energy *= EL_CHARGE_INV;

      tmpr.set(r, z, energy);
    }

#if defined TEMPERATURE_POSTPROC_BILINEAR
  Grid<double> temperature_src = tmpr;
  lib::bilinear_interpolation<Grid<double>>(temperature_src, tmpr);
#endif

}

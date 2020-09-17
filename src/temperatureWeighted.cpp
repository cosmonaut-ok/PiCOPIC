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

#include "temperatureWeighted.hpp"

void TemperatureWeighted::weight_temperature_cylindrical(string specie)
{
  for (auto ps = species_p.begin(); ps != species_p.end(); ++ps)
    if (specie.compare((**ps).name) == 0)
      for (auto i = (**ps).particles.begin(); i != (**ps).particles.end(); ++i)
      {
	double _mass = (**ps).mass * P_WEIGHT((**i));
        double p_vel_r = _mass * P_VEL_R((**i));
        double p_vel_phi = _mass * P_VEL_PHI((**i));
        double p_vel_z = _mass * P_VEL_Z((**i));
        double p_vel_full = lib::sq_rt(
          pow(p_vel_r, 2)
          + pow(p_vel_phi, 2)
          + pow(p_vel_z, 2)
          );

        weight_cylindrical<double>(geometry, &vel_r,
                                   P_POS_R((**i)),
                                   P_POS_Z((**i)),
                                   p_vel_r);

        weight_cylindrical<double>(geometry, &vel_phi,
                                   P_POS_R((**i)),
                                   P_POS_Z((**i)),
                                   p_vel_phi);

        weight_cylindrical<double>(geometry, &vel_z,
                                   P_POS_R((**i)),
                                   P_POS_Z((**i)),
                                   p_vel_z);

        weight_cylindrical<double>(geometry, &vel_full,
                                   P_POS_R((**i)),
                                   P_POS_Z((**i)),
                                   p_vel_full);
      }
}

void TemperatureWeighted::calc_temperature_cylindrical(string specie)
{
  // get specie
  SpecieP *speciep;

  for (auto ps = species_p.begin(); ps != species_p.end(); ++ps)
    if (specie.compare((**ps).name) == 0)
      speciep = (*ps);

  // calculate density initially
  density.density = 0;
  density.density.overlay_set(0);
  density.calc_density_cylindrical(specie);

  // clear temperature grid
  temperature = 0;
  temperature.overlay_set(0);

  // clear velociry components grid
  vel_r = 0;
  vel_r.overlay_set(0);
  vel_phi = 0;
  vel_phi.overlay_set(0);
  vel_z = 0;
  vel_z.overlay_set(0);

  // clear full velocity grid
  vel_full = 0;
  vel_full.overlay_set(0);

  // FIXME: temprary just copy full velocity
  // to temperature to get output
  // temperature.copy(vel_full);

  weight_temperature_cylindrical(specie);

  /////////////////////// FIXME: Legacy algorithm ////////////////////////////
  for (int r = 0; r < geometry->r_grid_amount; r++)
    for (int z = 0; z < geometry->z_grid_amount; z++)
    {
      double p_sum_2 = pow(vel_r(r, z), 2)
        + pow(vel_phi(r, z), 2) + pow(vel_z(r, z), 2);

      if (p_sum_2 < pow(REL_LIMIT, 2) * pow(speciep->mass, 2))
        temperature.set(r, z, (pow(vel_full(r, z), 2) - p_sum_2)
                        / (2. * speciep->mass * pow(density.density(r, z), 2)));
      else
      {
        double mc_2 = speciep->mass * LIGHT_VEL_POW_2;
        temperature.set(r, z, lib::sq_rt((pow(vel_full(r, z), 2) - p_sum_2)
                                         * LIGHT_VEL_POW_2
                                         / pow(density.density(r, z), 2)
                                         + pow(mc_2, 2)) - mc_2);
      }

      // convert joules to eV
      temperature.d_a(r, z, abs(speciep->charge));
    }

#if defined TEMPERATURE_POSTPROC_BILINEAR
  Grid<double> temperature_src = temperature;
  lib::bilinear_interpolation<Grid<double>>(temperature_src, temperature);
// #elif defined TEMPERATURE_POSTPROC_BICUBIC
//   lib::bicubic_interpolation(t_src, t, geom->n_grid_r, geom->n_grid_z);
#endif

}

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
        double dr = geometry->r_cell_size;
        double dz = geometry->z_cell_size;

        // finding number of i and k cell. example: dr = 0.5; r = 0.4; i =0
        unsigned int r_i = CELL_NUMBER(P_POS_R((**i)), dr);
        unsigned int z_k = CELL_NUMBER(P_POS_Z((**i)), dz);
        unsigned int r_i_shift = r_i - geometry->bottom_r_grid_number;
        unsigned int z_k_shift = z_k - geometry->left_z_grid_number;

        double p_vel_r = P_VEL_R((**i));
        double p_vel_phi = P_VEL_PHI((**i));
        double p_vel_z = P_VEL_Z((**i));

        double velocity = lib::sq_rt(pow(p_vel_r, 2) + pow(p_vel_phi, 2) + pow(p_vel_z, 2));

        vel_r.inc(r_i_shift, z_k_shift, p_vel_r);
        vel_phi.inc(r_i_shift, z_k_shift, p_vel_phi);
        vel_z.inc(r_i_shift, z_k_shift, p_vel_z);
        vel_full.inc(r_i_shift, z_k_shift, velocity);
        count.inc(r_i_shift, z_k_shift, 1);
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

  // clear count
  count = 0;
  count.overlay_set(0);

  // FIXME: temprary just copy full velocity
  // to temperature to get output
  // temperature.copy(vel_full);

  weight_temperature_cylindrical(specie);

  /////////////////////// FIXME: Legacy algorithm ////////////////////////////
  for (unsigned int r = 0; r < geometry->r_grid_amount; r++)
    for (unsigned int z = 0; z < geometry->z_grid_amount; z++)
    {
      double v_vec_sum_2 = vel_r(r, z) * vel_r(r, z)
        + vel_phi(r, z) * vel_phi(r, z) + vel_z(r, z) * vel_z(r, z);

      double v_sc_sum_2 = (vel_full(r, z) * vel_full(r, z) - v_vec_sum_2) / (count(r, z) * count(r, z));

      if (v_sc_sum_2 < pow(REL_LIMIT, 2))
        temperature.set(r, z, speciep->mass * v_sc_sum_2 / 2);
      else
      {
        double gamma = lib::get_gamma(v_sc_sum_2);
        temperature.set(r, z, speciep->mass * LIGHT_VEL_POW_2 * (gamma - 1));
      }

      // convert joules to eV
      temperature.d_a(r, z, abs(speciep->charge));
    }

// #if defined TEMPERATURE_POSTPROC_BILINEAR
//   lib::bilinear_interpolation(t_src, t, geom->n_grid_r, geom->n_grid_z);
// #elif defined TEMPERATURE_POSTPROC_BICUBIC
//   lib::bicubic_interpolation(t_src, t, geom->n_grid_r, geom->n_grid_z);
// #else
//   t = t_src;
// #endif

}

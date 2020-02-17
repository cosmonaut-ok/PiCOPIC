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
        double dr = geometry->r_cell_size;
        double dz = geometry->z_cell_size;
        double r1, r2, r3; // temp variables for calculation
        double dz1, dz2; // temp var.: width of k and k + 1 cell

        double v_0 = 0; // charge temperature Q/V, V - volume of particle
        double v_1 = 0; // volume of [i][k] cell
        double v_2 = 0; // volume of [i + 1][k] cell
        double value = 0;

        // finding number of i and k cell. example: dr = 0.5; r = 0.4; i =0
        unsigned int r_i = CELL_NUMBER(P_POS_R((**i)), dr);
        unsigned int z_k = CELL_NUMBER(P_POS_Z((**i)), dz);
        unsigned int r_i_shift = r_i - geometry->bottom_r_grid_number;
        unsigned int z_k_shift = z_k - geometry->left_z_grid_number;

        double p_vel_r = P_VEL_R((**i));
        double p_vel_phi = P_VEL_PHI((**i));
        double p_vel_z = P_VEL_Z((**i));

        double velocity = lib::sq_rt(pow(p_vel_r, 2) + pow(p_vel_phi, 2) + pow(p_vel_z, 2));

        double ro_vel;
        double ro_vel_r;
        double ro_vel_phi;
        double ro_vel_z;
        double ro_p;

        if (P_POS_R((**i)) > dr)
        {
          r1 =  P_POS_R((**i)) - 0.5 * dr;
          r2 = (r_i + 0.5) * dr;
          r3 = P_POS_R((**i)) + 0.5 * dr;
          v_0 = 2. * PI * dz * dr * P_POS_R((**i));

          ro_vel_r = P_MASS((**i)) * p_vel_r / v_0;
          ro_vel_phi = P_MASS((**i)) * p_vel_phi / v_0;
          ro_vel_z = P_MASS((**i)) * p_vel_z / v_0;
          ro_vel = P_MASS((**i)) * velocity / v_0;
          ro_p = P_MASS((**i)) / (**ps).mass / v_0;

          v_1 = CELL_VOLUME(r_i, dr, dz);
          v_2 = CELL_VOLUME(r_i + 1, dr, dz);
          dz1 = (z_k + 0.5) * dz - (P_POS_Z((**i)) - 0.5 * dz);
          dz2 = (P_POS_Z((**i)) + 0.5 * dz) - (z_k + 0.5) * dz;

          // weighting in ro[i][k] cell
          value = CYL_RNG_VOL(dz1, r1, r2) / v_1;
          vel_r.inc(r_i_shift, z_k_shift, ro_vel_r * value);
          vel_phi.inc(r_i_shift, z_k_shift, ro_vel_phi * value);
          vel_z.inc(r_i_shift, z_k_shift, ro_vel_z * value);
          vel_full.inc(r_i_shift, z_k_shift, ro_vel * value);

          // weighting in ro[i + 1][k] cell
          value = CYL_RNG_VOL(dz1, r2, r3) / v_2;
          vel_r.inc(r_i_shift + 1, z_k_shift, ro_vel_r * value);
          vel_phi.inc(r_i_shift + 1, z_k_shift, ro_vel_phi * value);
          vel_z.inc(r_i_shift + 1, z_k_shift, ro_vel_z * value);
          vel_full.inc(r_i_shift + 1, z_k_shift, ro_vel * value);

          // weighting in ro[i][k + 1] cell
          value = CYL_RNG_VOL(dz2, r1, r2) / v_1;
          vel_r.inc(r_i_shift, z_k_shift + 1, ro_vel_r * value);
          vel_phi.inc(r_i_shift, z_k_shift + 1, ro_vel_phi * value);
          vel_z.inc(r_i_shift, z_k_shift + 1, ro_vel_z * value);
          vel_full.inc(r_i_shift, z_k_shift + 1, ro_vel * value);

          // weighting in ro[i + 1][k + 1] cell
          value = CYL_RNG_VOL(dz2, r2, r3) / v_2;
          vel_r.inc(r_i_shift + 1, z_k_shift + 1, ro_vel_r * value);
          vel_phi.inc(r_i_shift + 1, z_k_shift + 1, ro_vel_phi * value);
          vel_z.inc(r_i_shift + 1, z_k_shift + 1, ro_vel_z * value);
          vel_full.inc(r_i_shift + 1, z_k_shift + 1, ro_vel * value);
        }
        else if (P_POS_R((**i)) <= dr / 2.)
        {
          r_i = 0;
          r1 =  0.;
          r2 = (r_i + 0.5) * dr;
          r3 = P_POS_R((**i)) + 0.5 * dr;
          dz1 = (z_k + 0.5) * dz - (P_POS_Z((**i)) - 0.5 * dz);
          dz2 = (P_POS_Z((**i)) + 0.5 * dz) - (z_k + 0.5) * dz;
          v_0 = PI * dz * (2. * P_POS_R((**i)) * P_POS_R((**i)) + dr * dr / 2.);

          ro_vel_r = P_MASS((**i)) * p_vel_r / v_0;
          ro_vel_phi = P_MASS((**i)) * p_vel_phi / v_0;
          ro_vel_z = P_MASS((**i)) * p_vel_z / v_0;
          ro_vel = P_MASS((**i)) * velocity / v_0;
          ro_p = P_MASS((**i)) / (**ps).mass / v_0;

          v_1 = CYL_VOL(dz, dr);
          v_2 = CELL_VOLUME(r_i + 1, dr, dz);

          // weighting in ro[i][k] cell
          value = PI * dz1 * (dr * dr / 2. - P_POS_R((**i)) * dr + P_POS_R((**i)) * P_POS_R((**i))) / v_1;
          vel_r.inc(r_i_shift, z_k_shift, ro_vel_r * value);
          vel_phi.inc(r_i_shift, z_k_shift, ro_vel_phi * value);
          vel_z.inc(r_i_shift, z_k_shift, ro_vel_z * value);
          vel_full.inc(r_i_shift, z_k_shift, ro_vel * value);

          // weighting in ro[i + 1][k] cell
          value = CYL_RNG_VOL(dz1, r2, r3) / v_2;
          vel_r.inc(r_i_shift + 1, z_k_shift, ro_vel_r * value);
          vel_phi.inc(r_i_shift + 1, z_k_shift, ro_vel_phi * value);
          vel_z.inc(r_i_shift + 1, z_k_shift, ro_vel_z * value);
          vel_full.inc(r_i_shift + 1, z_k_shift, ro_vel * value);

          // weighting in ro[i][k + 1] cell
          value = PI * dz2 * (dr * dr / 2. - P_POS_R((**i)) * dr + P_POS_R((**i)) * P_POS_R((**i))) / v_1;
          vel_r.inc(r_i_shift, z_k_shift + 1, ro_vel_r * value);
          vel_phi.inc(r_i_shift, z_k_shift + 1, ro_vel_phi * value);
          vel_z.inc(r_i_shift, z_k_shift + 1, ro_vel_z * value);
          vel_full.inc(r_i_shift, z_k_shift + 1, ro_vel * value);

          // weighting in ro[i + 1][k + 1] cell
          value = CYL_RNG_VOL(dz2, r2, r3) / v_2;
          vel_r.inc(r_i_shift + 1, z_k_shift + 1, ro_vel_r * value);
          vel_phi.inc(r_i_shift + 1, z_k_shift + 1, ro_vel_phi * value);
          vel_z.inc(r_i_shift + 1, z_k_shift + 1, ro_vel_z * value);
          vel_full.inc(r_i_shift + 1, z_k_shift + 1, ro_vel * value);
        }
        else
        {
          r1 = P_POS_R((**i)) - 0.5 * dr;
          r2 = (r_i + 0.5) * dr;
          r3 = P_POS_R((**i)) + 0.5 * dr;
          dz1 = (z_k + 0.5) * dz - (P_POS_Z((**i)) - 0.5 * dz);
          dz2 = (P_POS_Z((**i)) + 0.5 * dz) - (z_k + 0.5) * dz;
          v_0 = 2. * PI * dz * dr * P_POS_R((**i));

          ro_vel_r = P_MASS((**i)) * p_vel_r / v_0;
          ro_vel_phi = P_MASS((**i)) * p_vel_phi / v_0;
          ro_vel_z = P_MASS((**i)) * p_vel_z / v_0;
          ro_vel = P_MASS((**i)) * velocity / v_0;
          ro_p = P_MASS((**i)) / (**ps).mass / v_0;

          v_1 = CYL_VOL(dz, dr);
          v_2 = CELL_VOLUME(r_i + 1, dr, dz);

          // weighting in ro[i][k] cell
          value = CYL_RNG_VOL(dz1, r1, r2) / v_1;
          vel_r.inc(r_i_shift, z_k_shift, ro_vel_r * value);
          vel_phi.inc(r_i_shift, z_k_shift, ro_vel_phi * value);
          vel_z.inc(r_i_shift, z_k_shift, ro_vel_z * value);
          vel_full.inc(r_i_shift, z_k_shift, ro_vel * value);

          // weighting in ro[i + 1][k] cell
          value = CYL_RNG_VOL(dz1, r2, r3) / v_2;
          vel_r.inc(r_i_shift + 1, z_k_shift, ro_vel_r * value);
          vel_phi.inc(r_i_shift + 1, z_k_shift, ro_vel_phi * value);
          vel_z.inc(r_i_shift + 1, z_k_shift, ro_vel_z * value);
          vel_full.inc(r_i_shift + 1, z_k_shift, ro_vel * value);

          // weighting in ro[i][k + 1] cell
          value = CYL_RNG_VOL(dz2, r1, r2) / v_1;
          vel_r.inc(r_i_shift, z_k_shift + 1, ro_vel_r * value);
          vel_phi.inc(r_i_shift, z_k_shift + 1, ro_vel_phi * value);
          vel_z.inc(r_i_shift, z_k_shift + 1, ro_vel_z * value);
          vel_full.inc(r_i_shift, z_k_shift + 1, ro_vel * value);

          // weighting in ro[i + 1][k + 1] cell
          value = CYL_RNG_VOL(dz2, r2, r3) / v_2;
          vel_r.inc(r_i_shift + 1, z_k_shift + 1, ro_vel_r * value);
          vel_phi.inc(r_i_shift + 1, z_k_shift + 1, ro_vel_phi * value);
          vel_z.inc(r_i_shift + 1, z_k_shift + 1, ro_vel_z * value);
          vel_full.inc(r_i_shift + 1, z_k_shift + 1, ro_vel * value);
        }
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
  for (unsigned int r = 0; r < geometry->r_grid_amount; r++)
    for (unsigned int z = 0; z < geometry->z_grid_amount; z++)
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

// #if defined TEMPERATURE_POSTPROC_BILINEAR
//   lib::bilinear_interpolation(t_src, t, geom->n_grid_r, geom->n_grid_z);
// #elif defined TEMPERATURE_POSTPROC_BICUBIC
//   lib::bicubic_interpolation(t_src, t, geom->n_grid_r, geom->n_grid_z);
// #else
//   t = t_src;
// #endif

}

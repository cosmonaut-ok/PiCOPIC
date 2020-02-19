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

#include "density.hpp"

Density::Density(Geometry *geom, vector<SpecieP *> species) : geometry(geom)
{
  density = Grid<double> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
  species_p = species;

  density = 0;
  density.overlay_set(0);
}

void Density::calc_density_cylindrical(string specie)
{
  for (auto ps = species_p.begin(); ps != species_p.end(); ++ps)
    if (specie.compare((**ps).name) == 0)
      for (auto i = (**ps).particles.begin(); i != (**ps).particles.end(); ++i)
      {
        double dr = geometry->r_cell_size;
        double dz = geometry->z_cell_size;
        double r1, r2, r3; // temp variables for calculation
        double dz1, dz2; // temp var.: width of k and k + 1 cell

        double ro_p = 0; // charge density Q/V, V - volume of particle
        double v_0 = 0; // charge density Q/V, V - volume of particle
        double v_1 = 0; // volume of [i][k] cell
        double v_2 = 0; // volume of [i + 1][k] cell
        double value = 0;

        // finding number of i and k cell. example: dr = 0.5; r = 0.4; i =0
        unsigned int r_i = CELL_NUMBER(P_POS_R((**i)), dr);
        unsigned int z_k = CELL_NUMBER(P_POS_Z((**i)), dz);
        unsigned int r_i_shift = r_i - geometry->bottom_r_grid_number;
        unsigned int z_k_shift = z_k - geometry->left_z_grid_number;

        if (P_POS_R((**i)) > dr)
        {
          r1 =  P_POS_R((**i)) - 0.5 * dr;
          r2 = (r_i + 0.5) * dr;
          r3 = P_POS_R((**i)) + 0.5 * dr;
          v_0 = 2. * PI * dz * dr * P_POS_R((**i));

          ro_p = P_MASS((**i)) / EL_MASS / v_0;

          v_1 = CELL_VOLUME(r_i, dr, dz);
          v_2 = CELL_VOLUME(r_i + 1, dr, dz);
          dz1 = (z_k + 0.5) * dz - (P_POS_Z((**i)) - 0.5 * dz);
          dz2 = (P_POS_Z((**i)) + 0.5 * dz) - (z_k + 0.5) * dz;

          // weighting in ro[i][k] cell
          value = CYL_RNG_VOL(dz1, r1, r2) / v_1;
          density.inc(r_i_shift, z_k_shift, ro_p * value);

          // weighting in ro[i + 1][k] cell
          value = CYL_RNG_VOL(dz1, r2, r3) / v_2;
          density.inc(r_i_shift + 1, z_k_shift, ro_p * value);

          // weighting in ro[i][k + 1] cell
          value = CYL_RNG_VOL(dz2, r1, r2) / v_1;
          density.inc(r_i_shift, z_k_shift+1, ro_p * value);

          // weighting in ro[i + 1][k + 1] cell
          value = CYL_RNG_VOL(dz2, r2, r3) / v_2;
          density.inc(r_i_shift + 1, z_k_shift + 1, ro_p * value);
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

          ro_p = P_MASS((**i)) / EL_MASS / v_0;
          v_1 = CYL_VOL(dz, dr);
          v_2 = CELL_VOLUME(r_i + 1, dr, dz);

          // weighting in ro[i][k] cell
          value = PI * dz1 * (dr * dr / 2. - P_POS_R((**i)) * dr + P_POS_R((**i)) * P_POS_R((**i))) / v_1;
          density.inc(r_i_shift, z_k_shift, ro_p * value);

          // weighting in ro[i + 1][k] cell
          value = CYL_RNG_VOL(dz1, r2, r3) / v_2;
          density.inc(r_i_shift + 1, z_k_shift, ro_p * value);

          // weighting in ro[i][k + 1] cell
          value = PI * dz2 * (dr * dr / 2. - P_POS_R((**i)) * dr + P_POS_R((**i)) * P_POS_R((**i))) / v_1;
          density.inc(r_i_shift, z_k_shift + 1, ro_p* value);

          // weighting in ro[i + 1][k + 1] cell
          value = CYL_RNG_VOL(dz2, r2, r3) / v_2;
          density.inc(r_i_shift + 1, z_k_shift + 1, ro_p * value);
        }
        else
        {
          r1 = P_POS_R((**i)) - 0.5 * dr;
          r2 = (r_i + 0.5) * dr;
          r3 = P_POS_R((**i)) + 0.5 * dr;
          dz1 = (z_k + 0.5) * dz - (P_POS_Z((**i)) - 0.5 * dz);
          dz2 = (P_POS_Z((**i)) + 0.5 * dz) - (z_k + 0.5) * dz;
          v_0 = 2. * PI * dz * dr * P_POS_R((**i));

          ro_p = P_MASS((**i)) / EL_MASS / v_0;
          v_1 = CYL_VOL(dz, dr);
          v_2 = CELL_VOLUME(r_i + 1, dr, dz);

          // weighting in ro[i][k] cell
          value = CYL_RNG_VOL(dz1, r1, r2) / v_1;
          density.inc(r_i_shift, z_k_shift, ro_p * value);

          // weighting in ro[i + 1][k] cell
          value = CYL_RNG_VOL(dz1, r2, r3) / v_2;
          density.inc(r_i_shift + 1, z_k_shift, ro_p * value);

          // weighting in ro[i][k + 1] cell
          value = CYL_RNG_VOL(dz2, r1, r2) / v_1;
          density.inc(r_i_shift, z_k_shift + 1, ro_p * value);

          // weighting in ro[i + 1][k + 1] cell
          value = CYL_RNG_VOL(dz2, r2, r3) / v_2;
          density.inc(r_i_shift + 1, z_k_shift + 1, ro_p * value);
        }
      }
#if defined DENSITY_POSTPROC_BILINEAR
  Grid<double> density_src = density;
  lib::bilinear_interpolation<Grid<double>>(density_src, density);
#endif
}

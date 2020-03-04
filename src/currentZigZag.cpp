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

#include "currentZigZag.hpp"

void CurrentZigZag::current_distribution()
{
  current.overlay_set(0);

  double dr = geometry->r_cell_size;
  double dz = geometry->z_cell_size;
  double dt = time->step;

  for (auto ps = species_p.begin(); ps != species_p.end(); ++ps)
    for (auto i = (**ps).particles.begin(); i != (**ps).particles.end(); ++i)
    {
      //// Preparation

      // particle's position at t and \f$t + \Delta t \f$
      // TODO: is it P_POS and P_POS_OLD?
      double r_pos_old = P_POS_OLD_R((**i));
      double z_pos_old = P_POS_OLD_Z((**i));
      double r_pos_new = P_POS_R((**i));
      double z_pos_new = P_POS_Z((**i));

      double charge_over_dt = P_CHARGE((**i)) / dt;

      // finding number new and old cells
      int i_n = CELL_NUMBER(r_pos_new, dr);
      int i_o = CELL_NUMBER(r_pos_old, dr);
      int k_n = CELL_NUMBER(z_pos_new, dz);
      int k_o = CELL_NUMBER(z_pos_old, dz);

      //! shift also to take overlaying into account
      int i_o_shift = i_o - geometry->bottom_r_grid_number;
      int i_n_shift = i_n - geometry->bottom_r_grid_number;
      int k_o_shift = k_o - geometry->left_z_grid_number;
      int k_n_shift = k_n - geometry->left_z_grid_number;

      // find number of cell for middle point
      // between \f$ r_{old} \f$ and \f$ r_{new} \f$
      int i_middle = CELL_NUMBER((r_pos_new + r_pos_old) / 2., dr);
      double one_over_volume = 1 / (CELL_VOLUME(i_middle, dr, dz));

      //// Calculation

      //! \f$ F_{\phi} = \frac{Q_{prtl} * (\phi{new} + \phi_{new})}{\Delta t} \f$
      double F_phi = P_CHARGE((**i)) * P_VEL_PHI((**i));

      //! the formula is \f$ \frac{\phi_{old} + \phi_{new}}{2 \Delta \phi} - j_{old} \f$
      //! where \f$ j_{old} \f$ is 0 as well as \f$ j_{new} \f$ and
      //! \f$ \Delta \phi \f$ is full circle = \f$ 2 \pi r\f$,
      //! where \f$ r = (r_{position \; old} + r_{position \; new}) / 2\f$
      //! and (phi_pos_old + phi_pos_new) / 2 = \frac{v_{\phi} \Delta t}{2}
      double W_phi = P_VEL_PHI((**i)) * dt / (2. * PI * (r_pos_old + r_pos_new));

      double r_relay_pos = RELAY_POINT(i_o, i_n, r_pos_old, r_pos_new, dr);
      double z_relay_pos = RELAY_POINT(k_o, k_n, z_pos_old, z_pos_new, dz);

      double F_r1 = charge_over_dt * (r_relay_pos - r_pos_old);
      double F_r2 = charge_over_dt * (r_pos_new - r_relay_pos);
      double F_z1 = charge_over_dt * (z_relay_pos - z_pos_old);
      double F_z2 = charge_over_dt * (z_pos_new - z_relay_pos);

      double W_r1 = (r_pos_old + r_relay_pos) / 2 - i_o * dr;
      double W_r2 = (r_pos_new + r_relay_pos) / 2 - i_n * dr;
      double W_z1 = (z_pos_old + z_relay_pos) / 2 - k_o * dz;
      double W_z2 = (z_pos_new + z_relay_pos) / 2 - k_n * dz;

      double one_minus_Wr1 = 1 - W_r1;
      double one_minus_Wr2 = 1 - W_r2;
      double one_minus_Wz1 = 1 - W_z1;
      double one_minus_Wz2 = 1 - W_z2;
      double one_minus_Wphi = 1 - W_phi;

      double Jr_i1h_j1_k1 = one_over_volume * F_r1 * one_minus_Wphi * one_minus_Wz1;
      double Jr_i2h_j2_k2 = one_over_volume * F_r2 * one_minus_Wphi * one_minus_Wz2;
      double Jr_i1h_j1_1_k1 = one_over_volume * F_r1 * W_phi * one_minus_Wz1;
      double Jr_i2h_j2_1_k2 = one_over_volume * F_r2 * W_phi * one_minus_Wz2;
      double Jr_i1h_j1_k1_1 = one_over_volume * F_r1 * one_minus_Wphi * W_z1;
      double Jr_i2h_j2_k2_1 = one_over_volume * F_r2 * one_minus_Wphi * W_z2;
      double Jr_i1h_j1_1_k1_1 = one_over_volume * F_r1 * W_phi * W_z1;
      double Jr_i2h_j2_1_k2_1 = one_over_volume * F_r2 * W_phi * W_z2;

      // Let \f$ \phi_{relay} = \phi_{1} = 0
      // because geometry is axisymmetric
      double Jphi_i_jh_k = one_over_volume * F_phi * one_minus_Wr2 * one_minus_Wz2;
      double Jphi_i_jh_k_1 = one_over_volume * F_phi * one_minus_Wr2 * W_z1;
      double Jphi_i_1_jh_k = one_over_volume * F_phi * W_r2 * one_minus_Wz2;
      double Jphi_i_1_jh_k_1 = one_over_volume * F_phi * W_r2 * W_z2;

      double Jz_i1_j1_k1h = one_over_volume * F_z1 * one_minus_Wr1 * one_minus_Wphi;
      double Jz_i2_j2_k2h = one_over_volume * F_z2 * one_minus_Wr2 * one_minus_Wphi;
      double Jz_i1_j1_1_k1h = one_over_volume * F_z1 * one_minus_Wr1 * W_phi;
      double Jz_i2_j2_1_k2h = one_over_volume * F_z2 * one_minus_Wr2 * W_phi;
      double Jz_i1_1_j1_k1h = one_over_volume * F_z1 * W_r1 * one_minus_Wphi;
      double Jz_i2_1_j2_k2h = one_over_volume * F_z2 * W_r2 * one_minus_Wphi;
      double Jz_i1_1_j1_1_k1h = one_over_volume * F_z1 * W_r1 * W_phi;
      double Jz_i2_1_j2_1_k2h = one_over_volume * F_z2 * W_r2 * W_phi;

      // weight currents with calculated values
      current[0].inc(i_o_shift, k_o_shift, Jr_i1h_j1_k1 + Jr_i1h_j1_1_k1);
      current[0].inc(i_o_shift, k_o_shift+1, Jr_i1h_j1_k1_1 + Jr_i1h_j1_1_k1_1);
      current[0].inc(i_n_shift, k_n_shift, Jr_i2h_j2_k2 + Jr_i2h_j2_1_k2);
      current[0].inc(i_n_shift, k_n_shift+1, Jr_i2h_j2_k2_1 + Jr_i2h_j2_1_k2_1);

      // weight only to grid nodes, related to new position
      // because oldone is "zeroed", because of
      // 2.5D simplifications
      current[1].inc(i_n_shift, k_n_shift, Jphi_i_jh_k);
      current[1].inc(i_n_shift, k_n_shift+1, Jphi_i_jh_k_1);
      current[1].inc(i_n_shift+1, k_n_shift, Jphi_i_1_jh_k);
      current[1].inc(i_n_shift+1, k_n_shift+1, Jphi_i_1_jh_k_1);

      current[2].inc(i_o_shift, k_o_shift, Jz_i1_j1_k1h + Jz_i1_j1_1_k1h);
      current[2].inc(i_o_shift+1, k_o_shift, Jz_i1_1_j1_k1h + Jz_i1_1_j1_1_k1h);
      current[2].inc(i_n_shift, k_n_shift, Jz_i2_j2_k2h + Jz_i2_j2_1_k2h);
      current[2].inc(i_n_shift+1, k_n_shift, Jz_i2_1_j2_k2h + Jz_i2_1_j2_1_k2h);
    }
}

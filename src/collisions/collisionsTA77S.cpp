/*
 * This file is part of the PiCoPiC distribution (https://github.com/cosmonaut-ok/PiCoPiC).
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

#include "collisions/collisionsTA77S.hpp"

using namespace constant;

// say: TA77   - Takizuka, Abe; 1977; DOI: 10.1016/0021-9991(77)90099-7
//      TA77S18 - TA77S et al.; 2018; DOI: 10.1002/ctpp.201700121

using namespace std;

CollisionsTA77S::CollisionsTA77S (Geometry* _geometry, TimeSim *_time, vector <SpecieP *> _species_p) : Collisions ( _geometry, _time, _species_p)
{
}

void CollisionsTA77S::collide_single(double m_real_a, double m_real_b,
                                     double q_real_a, double q_real_b,
                                     Particle &pa, Particle &pb,
                                     double _density_a, double _density_b,
                                     double debye)
{
  // get required parameters
  double vr_a, vphi_a, vz_a, charge_a, mass_a, vr_b, vphi_b, vz_b,
    charge_b, mass_b, density_a, density_b;

  bool swap = false;

  // TA77S18: find weight ratio
  double w_ratio = P_WEIGHT(pa) / P_WEIGHT(pb);

  // a-particle should be lighter, than b-particle
  if (w_ratio <= 1)
  {
    vr_a = P_VEL_R(pa);
    vphi_a = P_VEL_PHI(pa);
    vz_a = P_VEL_Z(pa);
    charge_a = q_real_a;
    mass_a = m_real_a;
    density_a = _density_a;

    vr_b = P_VEL_R(pb);
    vphi_b = P_VEL_PHI(pb);
    vz_b = P_VEL_Z(pb);
    charge_b = q_real_b;
    mass_b = m_real_b;
    density_b = _density_b;
  }
  else
  {
    vr_a = P_VEL_R(pb);
    vphi_a = P_VEL_PHI(pb);
    vz_a = P_VEL_Z(pb);
    charge_a = q_real_b;
    mass_a = m_real_b;
    density_a = _density_b;

    vr_b = P_VEL_R(pa);
    vphi_b = P_VEL_PHI(pa);
    vz_b = P_VEL_Z(pa);
    charge_b = q_real_a;
    mass_b = m_real_a;
    density_b = _density_a;

    // swap particles when b-particle is lighter, than a-particle
    w_ratio = 1 / w_ratio;
    swap = true;
  }

  double density_lowest = min(density_a, density_b);
  double m_ab = mass_a * mass_b / (mass_a + mass_b);

  // relative velocity
  double ux = vr_a - vr_b;
  double uy = vphi_a - vphi_b;
  double uz = vz_a - vz_b;
  double u = algo::common::sq_rt(pow(ux, 2) + pow(uy, 2) + pow(uz, 2));

  // if ``u'' (relative velocity) is zero, particles can not collide
  if (u == 0) return;

  // get lambda Coulomb
  double lambda_coulomb = phys::plasma::coulomb_logarithm (mass_a, mass_b, debye, u);

  // TA77: calculate u perpendicular
  double u_p = algo::common::sq_rt(pow(ux, 2) + pow(uy, 2));

  // TA77: calculate scattering angles Theta and Phi
  // TA77 sub: find delta
  //// find variance of delta
  double variance_d = pow(charge_a, 2) * pow(charge_b, 2) * density_lowest
    * lambda_coulomb / (8 * constant::PI * pow(EPSILON0, 2) * pow(m_ab, 2) * pow(u, 3))
    * time->step;

  //// find standard deviation of delta
  double std_dev_d = algo::common::sq_rt(variance_d);
  //// find delta
  double delta = math::random::normal(std_dev_d);

  // TA77: find Theta angle
  double sin_Theta = 2 * delta / (1 + pow(delta, 2));
  double cos_Theta = 1 - 2 * pow(delta, 2) / (1 + pow(delta, 2));

  // TA77S18: find Theta_2 angle
  double sin_Theta2 = algo::common::sq_rt(
    w_ratio * pow(sin_Theta, 2)
    + w_ratio * ( 1 - w_ratio ) * pow(( 1 - cos_Theta), 2));
  double cos_Theta2 = 1 - w_ratio * (1 - cos_Theta);
  // TA77S18: sign of the sinus of Theta_angle
  // should be the same as sinus of Theta2_angle
  if (sin_Theta2 * sin_Theta < 0)
    sin_Theta2 = -sin_Theta2;

  // TA77: find Phi angle
  double Phi_angle = math::random::uniform_angle();
  double sin_Phi = sin(Phi_angle);
  double cos_Phi = cos(Phi_angle);

  // TA77S18: find Phi_2 angle
  //// find dzeta for Phi_2 angle
  double dzeta = math::random::uniform();
  double Phi2_angle = Phi_angle + 2 * constant::PI * ( 1 - algo::common::sq_rt(w_ratio)) * (1 - dzeta);
  double sin_Phi2 = sin(Phi2_angle);
  double cos_Phi2 = cos(Phi2_angle);

  // TA77: calculate delta u components
  double d_ux, d_uy, d_uz, d_ux2, d_uy2, d_uz2;
  // simplify calculations
  if (u_p == 0)
  {
    d_ux = u * sin_Theta * cos_Phi;
    d_uy = u * sin_Theta * sin_Phi;
    d_uz = -u * (1 - cos_Theta);

    d_ux2 = u * sin_Theta2 * cos_Phi2;
    d_uy2 = u * sin_Theta2 * sin_Phi2;
    d_uz2 = -u * (1 - cos_Theta2);
  }
  else
  {
    d_ux = (ux / u_p) * uz * sin_Theta * cos_Phi
      - (uy / u_p) * u * sin_Theta * sin_Phi
      - ux * (1 - cos_Theta);

    d_uy = (uy / u_p) * uz * sin_Theta * cos_Phi
      + (ux / u_p) * u * sin_Theta * sin_Phi
      - uy * (1 - cos_Theta);

    d_uz = u_p * sin_Theta * cos_Phi
      - uz * (1 - cos_Theta);

    ////

    d_ux2 = (ux / u_p) * uz * sin_Theta2 * cos_Phi2
      - (uy / u_p) * u * sin_Theta2 * sin_Phi2
      - ux * (1 - cos_Theta2);

    d_uy2 = (uy / u_p) * uz * sin_Theta2 * cos_Phi2
      + (ux / u_p) * u * sin_Theta2 * sin_Phi2
      - uy * (1 - cos_Theta2);

    d_uz2 = u_p * sin_Theta2 * cos_Phi2
      - uz * (1 - cos_Theta2);
  }

  double vr_a_new = vr_a + m_ab/mass_a * d_ux;
  double vphi_a_new = vphi_a + m_ab/mass_a * d_uy;
  double vz_a_new = vz_a + m_ab/mass_a * d_uz;

  double vr_b_new = vr_b + m_ab/mass_b * d_ux2;
  double vphi_b_new = vphi_b + m_ab/mass_b * d_uy2;
  double vz_b_new = vz_b + m_ab/mass_b * d_uz2;

  // set new velocity components
  if (swap)
  {
    P_VEL_R(pa) = vr_b_new;
    P_VEL_PHI(pa) = vphi_b_new;
    P_VEL_Z(pa) = vz_b_new;

    P_VEL_R(pb) = vr_a_new;
    P_VEL_PHI(pb) = vphi_a_new;
    P_VEL_Z(pb) = vz_a_new;
  }
  else
  {
    P_VEL_R(pa) = vr_a_new;
    P_VEL_PHI(pa) = vphi_a_new;
    P_VEL_Z(pa) = vz_a_new;

    P_VEL_R(pb) = vr_b_new;
    P_VEL_PHI(pb) = vphi_b_new;
    P_VEL_Z(pb) = vz_b_new;
  }

  if (vr_a_new > LIGHT_VEL || vphi_a_new > LIGHT_VEL || vz_a_new > LIGHT_VEL
      || vr_b_new > LIGHT_VEL || vphi_b_new > LIGHT_VEL || vz_b_new > LIGHT_VEL
      || vr_a_new * vr_a_new + vphi_a_new * vphi_a_new + vz_a_new * vz_a_new > LIGHT_VEL_POW_2 - 1e3
      || vr_b_new * vr_b_new + vphi_b_new * vphi_b_new + vz_b_new * vz_b_new > LIGHT_VEL_POW_2 - 1e3)
    LOG_S(ERROR) << "OH SHI! " << vr_a << " " << vphi_a << " " << vz_a << " => "
        << vr_a_new << " " << vphi_a_new << " " << vz_a_new << " -- "
        << charge_a << " " << mass_a << "; "
        << vr_b << " " << vphi_b << " " << vz_b << " => "
        << vr_b_new << " " << vphi_b_new << " " << vz_b_new << " -- "
        << charge_b << " " << mass_b;
}

void CollisionsTA77S::collide ()
{
  // pairing
  for (int i = 0; i < geometry->cell_amount[0]; ++i)
    for (int j = 0; j < geometry->cell_amount[1]; ++j)
    {
      unsigned int vec_size_ions = map_ion2cell(i, j).size();
      unsigned int vec_size_electrons = map_el2cell(i, j).size();

      // get temperatures, densities and debye length
      double temperature_el = get_el_temperature(i, j);
      double temperature_ion = get_ion_temperature(i, j);
      double density_el = get_el_density(i, j);
      double density_ion = get_ion_density(i, j);

      if (!isnormal(temperature_ion)) break;
      if (!isnormal(temperature_el)) break;
      if (!isnormal(density_ion)) break;
      if (!isnormal(density_el)) break;

      double debye = phys::plasma::debye_length(density_el, density_ion,
                                                temperature_el, temperature_ion);

      // TA77: case 1a
      // ions
      if (vec_size_ions % 2 == 0)
        for (unsigned int k = 0; k < vec_size_ions; k = k + 2)
          collide_single(mass_ion, mass_ion,
			 charge_ion, charge_ion,
                         (*map_ion2cell(i, j)[k]),
                         (*map_ion2cell(i, j)[k+1]),
                         density_ion, density_ion,
                         debye);
      // electrons
      if (vec_size_electrons % 2 == 0)
        for (unsigned int k = 0; k < vec_size_electrons; k = k + 2)
          collide_single(mass_el, mass_el,
			 charge_el, charge_el,
                         (*map_el2cell(i, j)[k]),
                         (*map_el2cell(i, j)[k+1]),
                         density_el, density_el,
                         debye);

      // TA77: case 1b
      // ions
      if (vec_size_ions % 2 != 0)
      {
        if (vec_size_ions >= 3)
        {
          // first 3 collisions in special way
          collide_single(mass_ion, mass_ion,
			 charge_ion, charge_ion,
                         (*map_ion2cell(i, j)[0]),
                         (*map_ion2cell(i, j)[1]),
                         density_ion, density_ion,
                         debye);
          collide_single(mass_ion, mass_ion,
			 charge_ion, charge_ion,
                         (*map_ion2cell(i, j)[1]),
                         (*map_ion2cell(i, j)[2]),
                         density_ion, density_ion,
                         debye);
          collide_single(mass_ion, mass_ion,
			 charge_ion, charge_ion,
                         (*map_ion2cell(i, j)[2]),
                         (*map_ion2cell(i, j)[0]),
                         density_ion, density_ion,
                         debye);
        }
        if (vec_size_ions >= 5)
          for (unsigned int k = 3; k < vec_size_ions; k = k + 2)
            collide_single(mass_ion, mass_ion,
			   charge_ion, charge_ion,
                           (*map_ion2cell(i, j)[k]),
                           (*map_ion2cell(i, j)[k+1]),
                           density_ion, density_ion,
                           debye);
      }
      // TA77: electrons, case 1b
      // electrons
      if (vec_size_electrons % 2 != 0)
      {
        if (vec_size_electrons >= 3)
        {
          // first 3 collisions in special way
          collide_single(mass_el, mass_el,
			 charge_el, charge_el,
                         (*map_el2cell(i, j)[0]),
                         (*map_el2cell(i, j)[1]),
                         density_el, density_el,
                         debye);
          collide_single(mass_el, mass_el,
			 charge_el, charge_el,
                         (*map_el2cell(i, j)[1]),
                         (*map_el2cell(i, j)[2]),
                         density_el, density_el,
                         debye);
          collide_single(mass_el, mass_el,
			 charge_el, charge_el,
                         (*map_el2cell(i, j)[2]),
                         (*map_el2cell(i, j)[0]),
                         density_el, density_el,
                         debye);
        }
        if (vec_size_electrons >= 5)
          for (unsigned int k = 3; k < vec_size_electrons; k = k + 2)
            collide_single(mass_el, mass_el,
			   charge_el, charge_el,
                           (*map_el2cell(i, j)[k]),
                           (*map_el2cell(i, j)[k+1]),
                           density_el, density_el,
                         debye);
      }

      // TA77: case 2a. electrons-ions
      if (vec_size_ions == vec_size_electrons)
        for (unsigned int k = 0; k < vec_size_electrons; ++k)
          collide_single(mass_el, mass_ion,
			 charge_el, charge_ion,
                         (*map_el2cell(i, j)[k]),
                         (*map_ion2cell(i, j)[k]),
                         density_el, density_ion,
                         debye);

      // TA77: case 2b. electrons-ions
      if (vec_size_ions > vec_size_electrons && vec_size_electrons > 0 && vec_size_ions > 0)
      {
        unsigned int c_i = floor( (float)vec_size_ions / (float)vec_size_electrons );
        double c_r = (float)vec_size_ions / (float)vec_size_electrons - c_i;

        unsigned int ions_1st_group = (c_i + 1) * c_r * vec_size_electrons;
        unsigned int els_1st_group = c_r * vec_size_electrons;

        unsigned int ions_2nd_group = c_i * (1 - c_r) * vec_size_electrons;
        // int els_2nd_group = (1 - c_r) * vec_size_electrons;

        // TA77: case 2b, 1st group, ions
        for (unsigned int fgi = 0; fgi < ions_1st_group; ++fgi)
        {
          int fge = floor(float(fgi) / float(c_i+1));
          collide_single(mass_el, mass_ion,
			 charge_el, charge_ion,
                         (*map_el2cell(i, j)[fge]),
                         (*map_ion2cell(i, j)[fgi]),
                         density_el, density_ion,
                         debye);
        }
        // TA77: case 2b, 2nd group, ions
        for (unsigned int fgi = 0; fgi < ions_2nd_group; ++fgi)
        {
          int fge = floor(float(fgi) / float(c_i));
          collide_single(mass_el, mass_ion,
			 charge_el, charge_ion,
                         (*map_el2cell(i, j)[fge+els_1st_group]),
                         (*map_ion2cell(i, j)[fgi+ions_1st_group]),
                         density_el, density_ion,
                         debye);
        }
      }
      // TA77: case 2b. electrons-ions
      if (vec_size_ions < vec_size_electrons && vec_size_electrons > 0 && vec_size_ions > 0)
      {
        double c_i = floor( (float)vec_size_electrons / (float)vec_size_ions );
        double c_r = (float)vec_size_electrons / (float)vec_size_ions - c_i;

        unsigned int els_1st_group = (c_i + 1) * c_r * vec_size_ions;
        unsigned int ions_1st_group = c_r * vec_size_ions;

        unsigned int els_2nd_group = c_i * (1 - c_r) * vec_size_ions;
        // int ions_2nd_group = (1 - c_r) * vec_size_ions;

        // TA77: case 2b, 1st group, electrons
        for (unsigned int fge = 0; fge < els_1st_group; ++fge)
        {
          int fgi = floor(float(fge) / float(c_i+1));
          collide_single(mass_el, mass_ion,
			 charge_el, charge_ion,
                         (*map_el2cell(i, j)[fge]),
                         (*map_ion2cell(i, j)[fgi]),
                         density_el, density_ion,
                         debye);
        }
        // TA77: case 2b, 2nd group, electrons
        for (unsigned int fge = 0; fge < els_2nd_group; ++fge)
        {
          int fgi = floor(float(fge) / float(c_i));
          collide_single(mass_el, mass_ion,
			 charge_el, charge_ion,
                         (*map_el2cell(i, j)[fge+els_1st_group]),
                         (*map_ion2cell(i, j)[fgi+ions_1st_group]),
                         density_el, density_ion,
                         debye);
        }
      }
    }
}

void CollisionsTA77S::correct_velocities()
{
  // TA77S18: correct velocities
  for (int i = 0; i < geometry->cell_amount[0]; ++i)
    for (int j = 0; j < geometry->cell_amount[1]; ++j)
    {
      unsigned int vec_size_ions = map_ion2cell(i, j).size();
      unsigned int vec_size_electrons = map_el2cell(i, j).size();

      // TA77S18: calculate delta V
      //// calculate summary moment components and total energy for t+delta_t
      double moment_new_r_ion = 0, moment_new_phi_ion = 0, moment_new_z_ion = 0,
        moment_new_r_el = 0, moment_new_phi_el = 0, moment_new_z_el = 0,
        E_tot_ion_new = 0, E_tot_el_new = 0;

      for (unsigned int k = 0; k < vec_size_ions; ++k)
      {
        double vr = P_VEL_R((*map_ion2cell(i, j)[k]));
        double vphi = P_VEL_PHI((*map_ion2cell(i, j)[k]));
        double vz = P_VEL_Z((*map_ion2cell(i, j)[k]));
        double v_sq = vr*vr + vphi*vphi + vz*vz;

        // mass of the macroparticle
        double mass_m = mass_ion * P_WEIGHT((*map_ion2cell(i, j)[k]));

        moment_new_r_ion += mass_m * vr;
        moment_new_phi_ion += mass_m * vphi;
        moment_new_z_ion += mass_m * vz;
        E_tot_ion_new += mass_m * v_sq / 2;
      }

      for (unsigned int k = 0; k < vec_size_electrons; ++k)
      {
        double vr = P_VEL_R((*map_el2cell(i, j)[k]));
        double vphi = P_VEL_PHI((*map_el2cell(i, j)[k]));
        double vz = P_VEL_Z((*map_el2cell(i, j)[k]));
        double v_sq = vr*vr + vphi*vphi + vz*vz;

        double mass_m = mass_el * P_WEIGHT((*map_el2cell(i, j)[k]));

        moment_new_r_el += mass_m * vr;
        moment_new_phi_el += mass_m * vphi;
        moment_new_z_el += mass_m * vz;
        E_tot_el_new += mass_m * v_sq / 2;
      }

      //// calculate delta V components and delta E total
      double mass_tot_ion = amount_tot_ion(i, j) * mass_ion;
      double mass_tot_el = amount_tot_el(i, j) * mass_el;
      double delta_V_r_ion = (moment_new_r_ion - moment_tot_ion[0](i, j)) / mass_tot_ion;
      double delta_V_phi_ion = (moment_new_phi_ion - moment_tot_ion[1](i, j)) / mass_tot_ion;
      double delta_V_z_ion = (moment_new_z_ion - moment_tot_ion[2](i, j)) / mass_tot_ion;
      double delta_V_r_el = (moment_new_r_el - moment_tot_el[0](i, j)) / mass_tot_el;
      double delta_V_phi_el = (moment_new_phi_el - moment_tot_el[1](i, j)) / mass_tot_el;
      double delta_V_z_el = (moment_new_z_el - moment_tot_el[2](i, j)) / mass_tot_el;

      //// calculate V_0 components
      double V_0_r_ion = moment_tot_ion[0](i, j) / mass_tot_ion;
      double V_0_phi_ion = moment_tot_ion[1](i, j) / mass_tot_ion;
      double V_0_z_ion = moment_tot_ion[2](i, j) / mass_tot_ion;
      double V_0_r_el = moment_tot_el[0](i, j) / mass_tot_el;
      double V_0_phi_el = moment_tot_el[1](i, j) / mass_tot_el;
      double V_0_z_el = moment_tot_el[2](i, j) / mass_tot_el;
      //// calculate delta E total
      double delta_E_tot_ion = E_tot_ion_new - energy_tot_ion(i, j);
      double delta_E_tot_el = E_tot_el_new - energy_tot_el(i, j);
      //// link E total to local vars
      double E_tot_ion = energy_tot_ion(i, j);
      double E_tot_el = energy_tot_el(i, j);
      double V_0_sq_ion = V_0_r_ion*V_0_r_ion + V_0_phi_ion*V_0_phi_ion + V_0_z_ion*V_0_z_ion;
      double V_0_sq_el = V_0_r_el*V_0_r_el + V_0_phi_el*V_0_phi_el + V_0_z_el*V_0_z_el;
      /// calculate alpha
      double alpha_ion =
        ( E_tot_ion - mass_tot_ion * V_0_sq_ion / 2 )
        / ( E_tot_ion
            + delta_E_tot_ion
            - mass_tot_ion
            * ( pow(V_0_r_ion + delta_V_r_ion, 2)
                + pow(V_0_phi_ion + delta_V_phi_ion, 2)
                + pow(V_0_z_ion + delta_V_z_ion, 2)
              )
            / 2
          );
      double alpha_el =
        ( E_tot_el - mass_tot_el * V_0_sq_el / 2 )
        / ( E_tot_el
            + delta_E_tot_el
            - mass_tot_el
            * ( pow(V_0_r_el + delta_V_r_el, 2)
                + pow(V_0_phi_el + delta_V_phi_el, 2)
                + pow(V_0_z_el + delta_V_z_el, 2)
              )
            / 2
          );

      // correct ion velocity
      if (isnormal(alpha_ion)) // FIXME: sometimes E_tot_ion = mass_tot_ion(i, j) * V_0_sq_ion / 2
      {
        // FIXME: quick and dirty workadound
        // caused negative values under squared root
        // and zeros
        if (alpha_ion < 0) alpha_ion = -alpha_ion;
        if (alpha_ion == 0) alpha_ion = 1;
        alpha_ion = algo::common::sq_rt(alpha_ion);

        for (unsigned int k = 0; k < vec_size_ions; ++k)
        {
          double vr = P_VEL_R((*map_ion2cell(i, j)[k]));
          double vphi = P_VEL_PHI((*map_ion2cell(i, j)[k]));
          double vz = P_VEL_Z((*map_ion2cell(i, j)[k]));

          double vr_corr = V_0_r_ion + alpha_ion * (vr - V_0_r_ion - delta_V_r_ion);
          double vphi_corr = V_0_phi_ion + alpha_ion * (vphi - V_0_phi_ion - delta_V_phi_ion);
          double vz_corr = V_0_z_ion + alpha_ion * (vz - V_0_z_ion - delta_V_z_ion);

          P_VEL_R((*map_ion2cell(i, j)[k])) = vr_corr;
          P_VEL_PHI((*map_ion2cell(i, j)[k])) = vphi_corr;
          P_VEL_Z((*map_ion2cell(i, j)[k])) = vz_corr;
        }
      }

      // correct electron velocity
      if (isnormal(alpha_el)) // FIXME: sometimes E_tot_el = mass_tot_el(i, j) * V_0_sq_el / 2
      {
        // FIXME: quick and dirty workadound
        // caused negative values under squared root
        // and zeros
        if (alpha_el < 0) alpha_el = -alpha_el;
        if (alpha_el == 0) alpha_el = 1;
        alpha_el = algo::common::sq_rt(alpha_el);

        for (unsigned int k = 0; k < vec_size_electrons; ++k)
        {
          double vr = P_VEL_R((*map_el2cell(i, j)[k]));
          double vphi = P_VEL_PHI((*map_el2cell(i, j)[k]));
          double vz = P_VEL_Z((*map_el2cell(i, j)[k]));

          double vr_corr = V_0_r_el + alpha_el * (vr - V_0_r_el - delta_V_r_el);
          double vphi_corr = V_0_phi_el + alpha_el * (vphi - V_0_phi_el - delta_V_phi_el);
          double vz_corr = V_0_z_el + alpha_el * (vz - V_0_z_el - delta_V_z_el);

          P_VEL_R((*map_el2cell(i, j)[k])) = vr_corr;
          P_VEL_PHI((*map_el2cell(i, j)[k])) = vphi_corr;
          P_VEL_Z((*map_el2cell(i, j)[k])) = vz_corr;
        }
      }
    }
}

void CollisionsTA77S::operator()()
{
  clear();
  sort_to_cells();
  random_sort();
  collect_weighted_params_tot_grid();
  collide();
  correct_velocities();
}

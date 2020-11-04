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

#include "collisionsP12.hpp"
#include <iostream>
// say: TA77   - Takizuka, Abe; 1977; DOI: 10.1016/0021-9991(77)90099-7
//      TA77S18 - P12 et al.; 2018; DOI: 10.1002/ctpp.201700121

using namespace std;

CollisionsP12::CollisionsP12 (Geometry* _geometry, TimeSim *_time, vector <SpecieP *> _species_p) : Collisions ( _geometry, _time, _species_p)
{}

void CollisionsP12::collide_single(double m_real_a, double m_real_b,
				   double q_real_a, double q_real_b,
                                   vector<double> &pa, vector<double> &pb,
                                   double _density_a, double _density_b,
                                   double _debye)
{
  // get required parameters
  double charge_a, mass_a, charge_b, mass_b, w_a, w_b, density_a, density_b;
  bool swap = false;

  // TA77S18: find weight ratio
  double w_ratio = P_WEIGHT(pa) / P_WEIGHT(pb);

  vector3d<double> v_a;
  vector3d<double> v_b;

  // a-particle should be lighter, than b-particle
  if (w_ratio <= 1)
  {
    v_a[0] = P_VEL_R(pa);
    v_a[1] = P_VEL_PHI(pa);
    v_a[2] = P_VEL_Z(pa);
    charge_a = q_real_a;
    mass_a = m_real_a;
    w_a = P_WEIGHT(pa);
    density_a = _density_a;

    v_b[0] = P_VEL_R(pb);
    v_b[1] = P_VEL_PHI(pb);
    v_b[2] = P_VEL_Z(pb);
    charge_b = q_real_b;
    mass_b = m_real_b;
    w_b = P_WEIGHT(pb);
    density_b = _density_b;
  }
  else
  {
    v_a[0] = P_VEL_R(pb);
    v_a[1] = P_VEL_PHI(pb);
    v_a[2] = P_VEL_Z(pb);
    charge_a = q_real_b;
    mass_a = m_real_b;
    w_a = P_WEIGHT(pb);
    density_a = _density_b;

    v_b[0] = P_VEL_R(pa);
    v_b[1] = P_VEL_PHI(pa);
    v_b[2] = P_VEL_Z(pa);
    charge_b = q_real_a;
    mass_b = m_real_a;
    w_b = P_WEIGHT(pa);
    density_b = _density_a;

    // swap particles when b-particle is lighter, than a-particle
    swap = true;
  }

  ////
  //// main calculation
  ////

  // do not collide particles with same velocities
  if (v_a == v_b) return;

  // find 4-momentum components for \alpha and \beta particles --- LAB frame
  vector3d<double> p_a, p_b;
  double p0_a, p0_b;

  double gamma_a = phys::rel::lorenz_factor(v_a.length2());
  double gamma_b = phys::rel::lorenz_factor(v_b.length2());

  p0_a = phys::rel::momentum_0(mass_a, v_a);
  p0_b = phys::rel::momentum_0(mass_b, v_b);
  p_a = phys::rel::momentum(mass_a, v_a);
  p_b = phys::rel::momentum(mass_b, v_b);

  // find \beta ( aka \frac{v}{c} ) --- center-of-mass frame
  vector3d<double> beta_cm;
  beta_cm = (p_a + p_b) / (p0_a + p0_b);

  // find lorenz-factor --- center-of-mass frame
  double v_cm_2 = beta_cm.length2() * LIGHT_VEL_POW_2;
  double gamma_cm = phys::rel::lorenz_factor(v_cm_2);

  // convert LAB frame's p^0 momentum (aka \gamma * m * c)
  // to center-of-mass frame's p^0
  double p0_a_cm = gamma_cm * ( p0_a - beta_cm.dot(p_a));
  double p0_b_cm = gamma_cm * ( p0_b - beta_cm.dot(p_b));

  // convert LAB frame's 3-momentum to center-of-mass frame's 3-momentum
  vector3d<double> p_a_cm, p_b_cm;
  p_a_cm = p_a;
  p_a_cm += beta_cm * beta_cm.dot(p_a) * ( gamma_cm - 1 ) / beta_cm.length2();
  p_a_cm -= beta_cm * gamma_cm * p0_a;

  p_b_cm[0] = -p_a_cm[0];
  p_b_cm[1] = -p_a_cm[1];
  p_b_cm[2] = -p_a_cm[2];

  // get gamma --- center-of-mass frame
  vector3d<double> v_cm = beta_cm * LIGHT_VEL;
  double gamma_a_cm = (1 - v_a.dot(v_cm) / LIGHT_VEL_POW_2) * gamma_a * gamma_cm;
  double gamma_b_cm = (1 - v_b.dot(v_cm) / LIGHT_VEL_POW_2) * gamma_b * gamma_cm;

  // get velocities --- center-of-mass frame
  vector3d<double> v_a_cm, v_b_cm;
  v_a_cm = p_a_cm / (mass_a * gamma_a_cm);
  v_b_cm = p_b_cm / (mass_b * gamma_b_cm);

  double v_rel_abs = ( v_a_cm.length() - v_b_cm.length() )
    / ( 1 - v_a_cm.dot(v_b_cm) / LIGHT_VEL_POW_2 );

  double gamma_rel = phys::rel::lorenz_factor(v_rel_abs * v_rel_abs);
  double p_rel_abs = gamma_rel * mass_b * v_rel_abs;

  // check if collision is possible
  if (p_rel_abs == 0) return;
  if (v_rel_abs < constant::MNZL) return;
  if (!isnormal(density_a)) return;
  if (!isnormal(density_b)) return;

  // take into account relativistic effects
  double debye;
  if (v_a.length2() > REL_LIMIT_POW_2 || v_b.length2() > REL_LIMIT_POW_2)
    debye = _debye * gamma_rel;
  else
    debye = _debye;

  // get coulomb logarithm
  double L_coulomb = phys::plasma::coulomb_logarithm (mass_a, mass_b, debye, v_rel_abs);

  // n weighted
  double density_w = density_a * density_b / min (density_a, density_b);
  double mass_ab = mass_a * mass_b / (mass_a + mass_b);

  // P12: calculate s_{1, 2} --- CM frame. FIXME: should be n_2 instead of density_highest
  double s12;
  if (v_a.length2() > REL_LIMIT_POW_2 || v_b.length2() > REL_LIMIT_POW_2)
  {
    s12 = density_w

      * time->step * L_coulomb * pow(charge_a * charge_b, 2)
      / ( 4 * constant::PI * pow(constant::EPSILON0 * LIGHT_VEL_POW_2, 2)
          * mass_a * gamma_a * mass_b * gamma_b)

      * gamma_cm * p_a.length()
      / (mass_a * gamma_a + mass_b * gamma_b)

      * pow(
        mass_a * gamma_a_cm * mass_b * gamma_b_cm * LIGHT_VEL_POW_2
        / p_a_cm.length2()
        + 1,
        2);
  }
  else
  {
    s12 = density_w *
      time->step * L_coulomb * pow(charge_a * charge_b, 2)
      / ( 4 * constant::PI * pow(constant::EPSILON0 * mass_ab, 2) * pow(v_rel_abs, 3) );
  }

  // Low temperature correction
  if (v_rel_abs < 1e6 // temperature should be less, than approx. 100 eV
      && density_w > 7e13) // density should be solid enough
  {
    double s12max = pow(FOUR_OVER_3 * constant::PI, ONE_OVER_3)
      * density_w * time->step
      * (mass_a + mass_b)
      / max(mass_a * pow(density_a, TWO_OVER_3),
            mass_b * pow(density_b, TWO_OVER_3))
      * v_rel_abs;

    s12 = min(s12, s12max);
  }

  double U = math::random::uniform();
  double cos_xi_cm;

  // P12: calculate scattering angle \xi^{*}
  if (s12 < 0.1)
    cos_xi_cm = 1 + s12 * log(U);
  else if (s12 < 3)
  {
    double A_inv = 0.00569578
      + 0.9560202 * s12
      - 0.508139 * pow(s12, 2)
      + 0.47913906 * pow(s12, 3)
      - 0.12788975 * pow(s12, 4)
      + 0.02389567 * pow(s12, 5);
    double A = 1./A_inv;
    cos_xi_cm = A_inv * log( exp(-A) + 2.* U * sinh(A) );
  }
  else if (s12 < 6.)
  {
    double A = 3. * exp( -s12 );
    double A_inv = 1./A;
    cos_xi_cm = A_inv * log( exp(-A) + 2.* U * sinh(A) );
  }
  else
    cos_xi_cm = 2 * U - 1;

  double sin_xi_cm = lib::sq_rt ( 1 - pow(cos_xi_cm, 2) );

  // find Phi angle
  double phi_cm = math::random::uniform_angle();
  double sin_phi_cm = sin(phi_cm);
  double cos_phi_cm = cos(phi_cm);

  vector3d<double> d_p, p_a_prime_cm, p_b_prime_cm;
  double p_a_cm_abs = p_a_cm.length();
  double p_cm_prp = lib::sq_rt(p_a_cm[0] * p_a_cm[0] + p_a_cm[1] * p_a_cm[1]);

  p_a_prime_cm[0] = p_a_cm[0] * p_a_cm[2] / p_cm_prp * sin_xi_cm * cos_phi_cm
    - p_a_cm[1] * p_a_cm_abs / p_cm_prp * sin_xi_cm * sin_phi_cm
    + p_a_cm[0] * cos_xi_cm;
  p_a_prime_cm[1] = p_a_cm[1] * p_a_cm[2] / p_cm_prp * sin_xi_cm * cos_phi_cm
    + p_a_cm[0] * p_a_cm_abs / p_cm_prp * sin_xi_cm * sin_phi_cm
    + p_a_cm[1] * cos_xi_cm;
  p_a_prime_cm[2] = -p_cm_prp * sin_xi_cm * cos_phi_cm + p_a_cm[2] * cos_xi_cm;

  // find deflection probability
  double U_defl = math::random::uniform();
  // double w_ab = min( w_a , w_b);

  vector3d<double> v_a_prime, v_b_prime;
  if (U_defl < w_a)
  {
    vector3d<double> p_a_prime;
    p_a_prime = p_a_prime_cm;
    p_a_prime += beta_cm * beta_cm.dot(p_a_prime_cm) * ( gamma_cm - 1 ) / beta_cm.length2();
    p_a_prime += beta_cm * gamma_cm * p0_a_cm;

    v_a_prime = p_a_prime / mass_a;
    double gamma_a_prime_inv = phys::rel::lorenz_factor_inv(v_a_prime.length2());
    v_a_prime *= gamma_a_prime_inv;
  }

  if (U_defl < w_b)
  {
    vector3d<double> p_b_prime;

    p_b_prime_cm[0] = -p_a_prime_cm[0];
    p_b_prime_cm[1] = -p_a_prime_cm[1];
    p_b_prime_cm[2] = -p_a_prime_cm[2];

    p_b_prime = p_b_prime_cm;
    p_b_prime += beta_cm * beta_cm.dot(p_b_prime_cm) * ( gamma_cm - 1 ) / beta_cm.length2();
    p_b_prime += beta_cm * gamma_cm * p0_b_cm;

    v_b_prime = p_b_prime / mass_b;
    double gamma_b_prime_inv = phys::rel::lorenz_factor_inv(v_b_prime.length2());
    v_b_prime *= gamma_b_prime_inv;
  }

  // collect mean values
  debye_mean += debye;
  s_mean += s12;
  L_mean += L_coulomb;
  ++ncol;

  ////
  //// end of main calculation
  ////

  // set new velocity components
  if (swap)
  {
    // deflect particles only with some probability
    // according to rejection scheme
    if (U_defl < w_a)
    {
      P_VEL_R(pa) = v_b_prime[0];
      P_VEL_PHI(pa) = v_b_prime[1];
      P_VEL_Z(pa) = v_b_prime[2];
    }

    if (U_defl < w_b)
    {
      P_VEL_R(pb) = v_a_prime[0];
      P_VEL_PHI(pb) = v_a_prime[1];
      P_VEL_Z(pb) = v_a_prime[2];
    }
  }
  else
  {
    if (U_defl < w_a)
    {
      P_VEL_R(pa) = v_a_prime[0];
      P_VEL_PHI(pa) = v_a_prime[1];
      P_VEL_Z(pa) = v_a_prime[2];
    }

    if (U_defl < w_b)
    {
      P_VEL_R(pb) = v_b_prime[0];
      P_VEL_PHI(pb) = v_b_prime[1];
      P_VEL_Z(pb) = v_b_prime[2];
    }
  }
}

void CollisionsP12::collide ()
{
  // reset all mean values
  debye_mean = 0;
  s_mean = 0;
  L_mean = 0;
  ncol = 0;

  // pairing
  for (int i = 0; i < geometry->r_grid_amount; ++i)
    for (int j = 0; j < geometry->z_grid_amount; ++j)
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

  s_mean /= ncol;
  if (s_mean > 1.2) // s_mean should not more, than 1 as angles should be small
  {
    // collect debug stuff
    L_mean /= ncol;
    debye_mean /= ncol;

    double d_el;
    double d_ion;
    double t_el;
    double t_ion;
    double g_counter;
    for (int i = 0; i < geometry->r_grid_amount; ++i)
      for (int j = 0; j < geometry->z_grid_amount; ++j)
      {
        // get temperatures, densities and debye length
        t_el += get_el_temperature(i, j);
        t_ion += get_ion_temperature(i, j);
        d_el += get_el_density(i, j);
        d_ion += get_ion_density(i, j);
        ++g_counter;
      }
    d_el /= g_counter;
    d_ion /= g_counter;
    t_el /= g_counter;
    t_ion /= g_counter;

    // print debug info
    LOG_S(WARNING) << "Value of s-parameter for P12 collisions is too large: " << s_mean;
    LOG_S(WARNING) << "\t other helpul info:";
    LOG_S(WARNING) << "\t\t domain number (r, z): "
                   << geometry->bottom_r_grid_number / geometry->r_grid_amount << "," << geometry->left_z_grid_number / geometry->z_grid_amount;
    LOG_S(WARNING) << "\t\t L coulomb mean: " << L_mean;
    LOG_S(WARNING) << "\t\t Debye length mean: " << L_mean;
    LOG_S(WARNING) << "\t\t Density mean (el, ion): " << d_el << "," << d_ion;
    LOG_S(WARNING) << "\t\t Temperature mean (el, ion): " << t_el << "," << t_ion;
  }
}

void CollisionsP12::run ()
{
  clear();
  sort_to_cells();
  random_sort();
  collect_weighted_params_tot_grid();
  collide();
  // correct_velocities();
}

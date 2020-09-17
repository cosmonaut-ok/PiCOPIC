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

#include "collisionsSK98.hpp"
#include <iostream>
// say: TA77   - Takizuka, Abe; 1977; DOI: 10.1016/0021-9991(77)90099-7
//      TA77S18 - SK98 et al.; 2018; DOI: 10.1002/ctpp.201700121

using namespace std;

CollisionsSK98::CollisionsSK98 (Geometry* _geometry, TimeSim *_time, vector <SpecieP *> _species_p) : Collisions ( _geometry, _time, _species_p)
{}

void CollisionsSK98::collide_single(double m_real_a, double m_real_b,
                                    vector<double> &pa, vector<double> &pb,
                                    double _density_a, double _density_b,
                                    double debye)
{
  // get required parameters
  double charge_a, mass_a, charge_b, mass_b, density_a, density_b; //, w_a, w_b;
  bool swap = false;

  // TA77S18: find weight ratio
  double w_ratio = P_MASS(pa) * m_real_b / (P_MASS(pb) * m_real_a);

  vector3d<double> v_a;
  vector3d<double> v_b;

  // a-particle should be lighter, than b-particle
  if (w_ratio <= 1)
  {
    v_a[0] = P_VEL_R(pa);
    v_a[1] = P_VEL_PHI(pa);
    v_a[2] = P_VEL_Z(pa);
    charge_a = P_CHARGE(pa);
    mass_a = P_MASS(pa);
    // w_a = mass_a / m_real_a;
    density_a = _density_a;

    v_b[0] = P_VEL_R(pb);
    v_b[1] = P_VEL_PHI(pb);
    v_b[2] = P_VEL_Z(pb);
    charge_b = P_CHARGE(pb);
    mass_b = P_MASS(pb);
    // w_b = mass_b / m_real_b;
    density_b = _density_b;
  }
  else
  {
    v_a[0] = P_VEL_R(pb);
    v_a[1] = P_VEL_PHI(pb);
    v_a[2] = P_VEL_Z(pb);
    charge_a = P_CHARGE(pb);
    mass_a = P_MASS(pb);
    // w_a = mass_a / m_real_b;
    density_a = _density_b;

    v_b[0] = P_VEL_R(pa);
    v_b[1] = P_VEL_PHI(pa);
    v_b[2] = P_VEL_Z(pa);
    charge_b = P_CHARGE(pa);
    mass_b = P_MASS(pa);
    // w_b = mass_b / m_real_a;
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
  // p_b_cm = p_b;
  // p_b_cm += beta_cm * beta_cm.dot(p_b) * ( gamma_cm - 1 ) / beta_cm.length2();
  // p_b_cm -= beta_cm * gamma_cm * p0_b;

  // get gamma --- center-of-mass frame
  double gamma_a_cm = p0_a_cm / mass_a / LIGHT_VEL;
  double gamma_b_cm = p0_b_cm / mass_b / LIGHT_VEL;

  // get velocities --- center-of-mass frame
  vector3d<double> v_a_cm, v_b_cm;
  v_a_cm = p_a_cm / (mass_a * gamma_a_cm);
  v_b_cm = p_b_cm / (mass_b * gamma_b_cm);

  vector3d<double> p_cm; // (0, 0, p_a_cm.length());
  p_cm = p_a_cm;

  double v_rel_abs = ( v_a_cm.length() - v_b_cm.length() )
    / ( 1 - v_a_cm.dot(v_b_cm) / LIGHT_VEL_POW_2 );

  double gamma_rel = phys::rel::lorenz_factor(v_rel_abs * v_rel_abs);
  double p_rel_abs = gamma_rel * mass_b * v_rel_abs;

  double density_lowest = min(density_a, density_b);

  // check if collision is possible
  if (p_rel_abs == 0) return;
  if (v_rel_abs < constant::MNZL) return;

  // get coulomb logarithm
  double L_coulomb = phys::plasma::coulomb_logarithm (mass_a, mass_b,
                                                      debye, v_rel_abs);

  // FIXME: figure out, why coulomb logarithm can acquire negative values
  if (L_coulomb <= 0) return;

  // get collision frequency
  double coll_freq = phys::plasma::collision_freqency ( charge_a, charge_b,
                                                        density_lowest,
                                                        L_coulomb,
                                                        p_rel_abs,
                                                        v_rel_abs );

  double variance_d = coll_freq * time->step;

  //// find standard deviation of delta
  double std_dev_d = lib::sq_rt(variance_d);
  //// find delta
  double delta = math::random::normal(std_dev_d);

  double sin_theta = 2 * delta / (1 + pow(delta, 2));
  double cos_theta = 1 - 2 * pow(delta, 2) / (1 + pow(delta, 2));

  // find theta angle for the COM frame
  double tg_theta_cm = sin_theta / (gamma_cm * ( cos_theta - lib::sq_rt(v_cm_2) / v_rel_abs ) );
  double sin_theta_cm = sin(atan(tg_theta_cm));
  double cos_theta_cm = cos(atan(tg_theta_cm));

  // find Phi angle
  double phi_angle = math::random::uniform_angle();
  double sin_phi = sin(phi_angle);
  double cos_phi = cos(phi_angle);

  vector3d<double> d_p;
  double p_cm_abs = p_cm.length();
  double p_cm_prp = lib::sq_rt(p_cm[0] * p_cm[0] + p_cm[1] * p_cm[1]);

  d_p[0] = p_cm[0] * p_cm[2] / p_cm_prp * sin_theta_cm * cos_phi
    - p_cm[1] * p_cm_abs / p_cm_prp * sin_theta_cm * sin_phi
    - p_cm[0] * ( 1 - cos_theta_cm);
  d_p[1] = p_cm[1] * p_cm[2] / p_cm_prp * sin_theta_cm * cos_phi
    + p_cm[0] * p_cm_abs / p_cm_prp * sin_theta_cm * sin_phi
    - p_cm[1] * ( 1 - cos_theta_cm);
  d_p[2] = -p_cm_prp * sin_theta_cm * cos_phi - p_cm[2] * ( 1 - cos_theta_cm );

  // get p_a_prime_cm and p_b_prime_cm --- COM frame
  vector3d<double> p_a_prime_cm, p_b_prime_cm;
  p_a_prime_cm = p_a_cm + d_p;
  p_b_prime_cm = p_b_cm - d_p;

  vector3d<double> p_a_prime, p_b_prime;
  p_a_prime = p_a_prime_cm;
  p_a_prime += beta_cm * beta_cm.dot(p_a_prime_cm) * ( gamma_cm - 1 ) / beta_cm.length2();
  p_a_prime += beta_cm * gamma_cm * p0_a_cm;

  p_b_prime = p_b_prime_cm;
  p_b_prime += beta_cm * beta_cm.dot(p_b_prime_cm) * ( gamma_cm - 1 ) / beta_cm.length2();
  p_b_prime += beta_cm * gamma_cm * p0_b_cm;

  // calculate velocities
  vector3d<double> v_a_prime, v_b_prime;
  v_a_prime = p_a_prime / mass_a;
  v_b_prime = p_b_prime / mass_b;

  double gamma_a_prime_inv = phys::rel::lorenz_factor_inv(v_a_prime.length2());
  double gamma_b_prime_inv = phys::rel::lorenz_factor_inv(v_b_prime.length2());

  v_a_prime *= gamma_a_prime_inv;
  v_b_prime *= gamma_b_prime_inv;

  ////
  //// end of main calculation
  ////

  // set new velocity components
  if (swap)
  {
    P_VEL_R(pa) = v_b_prime[0];
    P_VEL_PHI(pa) = v_b_prime[1];
    P_VEL_Z(pa) = v_b_prime[2];

    P_VEL_R(pb) = v_a_prime[0];
    P_VEL_PHI(pb) = v_a_prime[1];
    P_VEL_Z(pb) = v_a_prime[2];
  }
  else
  {
    P_VEL_R(pa) = v_a_prime[0];
    P_VEL_PHI(pa) = v_a_prime[1];
    P_VEL_Z(pa) = v_a_prime[2];

    P_VEL_R(pb) = v_b_prime[0];
    P_VEL_PHI(pb) = v_b_prime[1];
    P_VEL_Z(pb) = v_b_prime[2];
  }
}

void CollisionsSK98::collide ()
{
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
          collide_single(PROTON_MASS, PROTON_MASS,
                         (*map_ion2cell(i, j)[k]),
                         (*map_ion2cell(i, j)[k+1]),
                         density_ion, density_ion,
                         debye);
      // electrons
      if (vec_size_electrons % 2 == 0)
        for (unsigned int k = 0; k < vec_size_electrons; k = k + 2)
          collide_single(EL_MASS, EL_MASS,
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
          collide_single(PROTON_MASS, PROTON_MASS,
                         (*map_ion2cell(i, j)[0]),
                         (*map_ion2cell(i, j)[1]),
                         density_ion, density_ion,
                         debye);
          collide_single(PROTON_MASS, PROTON_MASS,
                         (*map_ion2cell(i, j)[1]),
                         (*map_ion2cell(i, j)[2]),
                         density_ion, density_ion,
                         debye);
          collide_single(PROTON_MASS, PROTON_MASS,
                         (*map_ion2cell(i, j)[2]),
                         (*map_ion2cell(i, j)[0]),
                         density_ion, density_ion,
                         debye);
        }
        if (vec_size_ions >= 5)
          for (unsigned int k = 3; k < vec_size_ions; k = k + 2)
            collide_single(PROTON_MASS, PROTON_MASS,
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
          collide_single(EL_MASS, EL_MASS,
                         (*map_el2cell(i, j)[0]),
                         (*map_el2cell(i, j)[1]),
                         density_el, density_el,
                         debye);
          collide_single(EL_MASS, EL_MASS,
                         (*map_el2cell(i, j)[1]),
                         (*map_el2cell(i, j)[2]),
                         density_el, density_el,
                         debye);
          collide_single(EL_MASS, EL_MASS,
                         (*map_el2cell(i, j)[2]),
                         (*map_el2cell(i, j)[0]),
                         density_el, density_el,
                         debye);
        }
        if (vec_size_electrons >= 5)
          for (unsigned int k = 3; k < vec_size_electrons; k = k + 2)
            collide_single(EL_MASS, EL_MASS,
                           (*map_el2cell(i, j)[k]),
                           (*map_el2cell(i, j)[k+1]),
                           density_el, density_el,
                         debye);
      }

      // TA77: case 2a. electrons-ions
      if (vec_size_ions == vec_size_electrons)
        for (unsigned int k = 0; k < vec_size_electrons; ++k)
          collide_single(EL_MASS, PROTON_MASS,
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
          collide_single(EL_MASS, PROTON_MASS,
                         (*map_el2cell(i, j)[fge]),
                         (*map_ion2cell(i, j)[fgi]),
                         density_el, density_ion,
                         debye);
        }
        // TA77: case 2b, 2nd group, ions
        for (unsigned int fgi = 0; fgi < ions_2nd_group; ++fgi)
        {
          int fge = floor(float(fgi) / float(c_i));
          collide_single(EL_MASS, PROTON_MASS,
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
          collide_single(EL_MASS, PROTON_MASS,
                         (*map_el2cell(i, j)[fge]),
                         (*map_ion2cell(i, j)[fgi]),
                         density_el, density_ion,
                         debye);
        }
        // TA77: case 2b, 2nd group, electrons
        for (unsigned int fge = 0; fge < els_2nd_group; ++fge)
        {
          int fgi = floor(float(fge) / float(c_i));
          collide_single(EL_MASS, PROTON_MASS,
                         (*map_el2cell(i, j)[fge+els_1st_group]),
                         (*map_ion2cell(i, j)[fgi+ions_1st_group]),
                         density_el, density_ion,
                         debye);
        }
      }
    }
}

void CollisionsSK98::run ()
{
  clear();
  sort_to_cells();
  random_sort();
  collect_weighted_params_tot_grid();
  collide();
  // correct_velocities();
}

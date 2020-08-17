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

#include "collisionsSentokuM.hpp"
#include <iostream>
// say: TA77   - Takizuka, Abe; 1977; DOI: 10.1016/0021-9991(77)90099-7
//      TITU18 - SentokuM et al.; 2018; DOI: 10.1002/ctpp.201700121

using namespace std;

CollisionsSentokuM::CollisionsSentokuM (Geometry* _geometry, TimeSim *_time, vector <SpecieP *> _species_p) : Collisions ( _geometry, _time, _species_p)
{
  //! TA77: make dummies to place particles there
  // map_el2cell = Grid<vector< vector<double> * >> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
  // map_ion2cell = Grid<vector< vector<double> * >> (geometry->r_grid_amount, geometry->z_grid_amount, 2);

}

void CollisionsSentokuM::collide_single(int i, int j, double m_real_a, double m_real_b,
                                vector<double> &pa, vector<double> &pb)
{
  // get required parameters
  double charge_a, mass_a, charge_b, mass_b;
  bool swap = false;

  // TITU18: find weight ratio
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

    v_b[0] = P_VEL_R(pb);
    v_b[1] = P_VEL_PHI(pb);
    v_b[2] = P_VEL_Z(pb);
    charge_b = P_CHARGE(pb);
    mass_b = P_MASS(pb);
  }
  else
  {
    v_a[0] = P_VEL_R(pb);
    v_a[1] = P_VEL_PHI(pb);
    v_a[2] = P_VEL_Z(pb);
    charge_a = P_CHARGE(pb);
    mass_a = P_MASS(pb);

    v_b[0] = P_VEL_R(pa);
    v_b[1] = P_VEL_PHI(pa);
    v_b[2] = P_VEL_Z(pa);
    charge_b = P_CHARGE(pa);
    mass_b = P_MASS(pa);

    // swap particles when b-particle is lighter, than a-particle
    w_ratio = 1 / w_ratio;
    swap = true;
  }

  //// main calculation

  // get p^0_a, p^0_b, \gamma_a and \gamma_b --- LAB frame
  double p0_a = phys::rel::momentum_0(mass_a, v_a);
  double p0_b = phys::rel::momentum_0(mass_a, v_b);

  // get p_a and p_b --- LAB frame
  vector3d<double> p_a;
  vector3d<double> p_b;
  p_a = phys::rel::momentum (mass_a, v_a);
  p_b = phys::rel::momentum (mass_b, v_b);

  // do not collide, if momentums are equal
  if (v_a == v_b) return;

  // get v of CM frame
  vector3d<double> v_cm(0,0,0);
  v_cm = (p_a + p_b) / (p0_a + p0_b);

  // get \gamma of CM frame
  double gamma_cm = phys::rel::lorenz_factor(v_cm.length2());

  // get p^0_a and p^0_b  --- CM frame
  double p0_a_cm = gamma_cm * ( p0_a - v_cm.dot(p_a) );
  double p0_b_cm = gamma_cm * ( p0_b - v_cm.dot(p_b) );

  // LOG_S(ERROR) << "v_cm: " << v_cm[0] << " " << v_cm[1] << " " << v_cm[2];
  // LOG_S(ERROR) << "p_a : " << p_a[0] << " " << p_a[1] << " " << p_a[2];
  // LOG_S(ERROR) << "p0_a, p0_a_cm : " << p0_a << " " << p0_a_cm;
  // LOG_S(ERROR);

  // P, P_a and P_b --- CM frame
  vector3d<double> p_cm;
  vector3d<double> p_a_cm;
  vector3d<double> p_b_cm;

  // get v_a and v_b  --- CM frame
  vector3d<double> v_a_cm;
  // find p_a in CM frame initially
  p_a_cm = p_a;
  p_a_cm += v_cm * v_cm.dot(p_a) * (gamma_cm - 1) / v_cm.length2();
  p_a_cm -= v_cm * p0_a * gamma_cm;
  p_cm = p_a_cm; // just set p_cm. FIXME: should it be like this?
  v_a_cm = p_a_cm / p0_a_cm; // v_a_cm = P_a_cm / (gamma_a_cm * m_a) -> P_a_cm / p0_a_cm

  vector3d<double> v_b_cm(0,0,0);
  p_b_cm = p_b;
  p_b_cm += v_cm * v_cm.dot(p_b) * (gamma_cm - 1) / v_cm.length2();
  p_b_cm -= v_cm * p0_b * gamma_cm;
  v_b_cm = p_b_cm / p0_b_cm;

  // get p and v --- OPR frame (alpha particle is in rest)
  vector3d<double> p_rel;
  vector3d<double> v_rel;
  p_rel = p_b_cm - p_a_cm;
  v_rel = v_a_cm - v_b_cm;
  v_rel /= 1 - v_a_cm.dot(v_b_cm) / LIGHT_VEL_POW_2;

  // check if collision is possible
  if (p_rel.length2() < constant::MNZL) return;
  if (v_rel.length2() < constant::MNZL) return;

  // get densities and electron temperature
  double temperature_el = get_el_temperature(i, j);
  double density_el = get_el_density(i, j);
  double density_lowest = min(get_el_density(i, j), get_ion_density(i, j));
  double density_highest = max(get_el_density(i, j), get_ion_density(i, j));

  // get debye length
  double debye = phys::plasma::debye_length(density_el, temperature_el);

  // get coulomb logarithm
  double L_coulomb = phys::plasma::coulomb_logarithm (mass_a, mass_b,
                                                      debye, v_rel.length());
  // TODO: figure out, why coulomb logarithm can acquire negative values
  if (L_coulomb < 1) L_coulomb = 1;

  // get collision frequency
  double coll_freq = phys::plasma::collision_freqency (charge_a, charge_b,
                                                       density_lowest,
                                                       L_coulomb,
                                                       p_rel.length(),
                                                       v_rel.length());

  // LOG_S(ERROR) << "debye: " << debye
  //              << " L Coulomb: " << L_coulomb
  //              << " Collisions freq: " << coll_freq;

  // get Theta angle --- OPR frame
  double variance_tg_theta_half = coll_freq * time->step;
  double tg_theta = 2 * math::random::normal(lib::sq_rt(variance_tg_theta_half));
  double sin_theta = sin(atan(tg_theta));
  double cos_theta = cos(atan(tg_theta));

  // LOG_S(ERROR) << "variance tg theta: " << variance_tg_theta_half
  //              << " tg theta: " << tg_theta
  //              << " sin theta: " << sin_theta
  //              << " cos theta: " << cos_theta;

  // get Phi angle --- CM frame
  double phi_cm = math::random::uniform_angle();
  double cos_phi_cm = cos(phi_cm);
  double sin_phi_cm = sin(phi_cm);

  // LOG_S(ERROR) << "Phi_CM: " << phi_cm
  //              << " cos Phi_CM: " << cos_phi_cm
  //              << " sin Phi_CM: " << sin_phi_cm;

  // get \beta and \beta_{cm}
  double beta = 0.01; // v_a.length() * v_b.length() / LIGHT_VEL_POW_2;
  double beta_cm = 0.02; //v_a_cm.length() * v_b_cm.length() / LIGHT_VEL_POW_2;

  // get Theta angle --- CM frame
  double tg_theta_cm = sin_theta / (gamma_cm * (cos_theta - beta/beta_cm));
  double sin_theta_cm = sin(atan(tg_theta_cm));
  double cos_theta_cm = cos(atan(tg_theta_cm));


  // LOG_S(ERROR) << "tg Theta_CM: " << tg_theta_cm
  //              << " cos Theta_CM: " << cos_theta_cm
  //              << " sin Theta_CM: " << sin_theta_cm;

  // get Theta angle --- rotation // FIXME: not P_x-P_y ?
  double cos_theta_r = p_cm[2] / p_cm.length();
  double sin_theta_r = lib::sq_rt(p_cm[0] * p_cm[0] + p_cm[1] * p_cm[1]) / p_cm.length();

  // LOG_S(ERROR) << "cos Theta_r: " << cos_theta_r
  //              << " sin Theta_r: " << sin_theta_r;

  // get Phi angle --- rotation
  double cos_phi_r = p_cm[0] / lib::sq_rt(p_cm[0] * p_cm[0] + p_cm[1] * p_cm[1]);
  double sin_phi_r = p_cm[1] / lib::sq_rt(p_cm[0] * p_cm[0] + p_cm[1] * p_cm[1]);

  // LOG_S(ERROR) << "cos Phi_r: " << cos_phi_r
  //              << " sin Phi_r: " << sin_phi_r;

  double p_cm_l = p_cm.length();
  vector3d<double> d_p(0,0,0);

  d_p[0] = p_cm_l * (
    cos_theta_r * cos_phi_r * sin_theta_cm * cos_phi_cm
    - cos_theta_r * sin_phi_r * sin_theta_cm * sin_phi_cm
    + sin_theta_r * sin_theta_cm) - p_cm[0];
  d_p[1] = p_cm_l * (
    sin_phi_r * sin_theta_cm * cos_phi_cm
    + cos_phi_r * sin_theta_cm * sin_phi_cm) - p_cm[1];
  d_p[2] = p_cm_l * (
    - sin_theta_r * cos_phi_r * sin_theta_cm * cos_phi_cm
    - sin_theta_r * sin_phi_r * sin_theta_cm * sin_phi_cm
    + cos_theta_r * cos_theta_cm) - p_cm[2];

  // LOG_S(ERROR) << "P_a    : " << p_a[0] << " " << p_a[1] << " " << p_b[2] << " " << p0_a;
  // LOG_S(ERROR) << "P_b    : " << p_b[0] << " " << p_b[1] << " " << p_b[2] << " " << p0_b;
  // LOG_S(ERROR) << "P_a_CM : " << p_a_cm[0] << " " << p_a_cm[1] << " " << p_b_cm[2] << " " << p0_a_cm;
  // LOG_S(ERROR) << "P_b_CM : " << p_b_cm[0] << " " << p_b_cm[1] << " " << p_b_cm[2] << " " << p0_b_cm;
  // LOG_S(ERROR) << "Delta P: " << d_p[0] << " " << d_p[1] << " " << d_p[2];

  vector3d<double> p_a_bar_cm;
  vector3d<double> p_b_bar_cm;

  p_a_bar_cm = p_a_cm + d_p;
  p_b_bar_cm = p_b_cm - d_p;

  vector3d<double> p_a_bar;
  vector3d<double> p_b_bar;

  p_a_bar = p_a_bar_cm;
  p_a_bar += v_cm * v_cm.dot(p_a_bar_cm) * ( ( gamma_cm - 1 ) / v_cm.length2() );
  p_a_bar += v_cm * gamma_cm * p0_a_cm;

  // LOG_S(ERROR);
  // LOG_S(ERROR) << "v_cm      : " << v_cm[0] << " " << v_cm[1] << " " << v_cm[2];
  // LOG_S(ERROR) << "p_a_cm    : " << p_a_cm[0] << " " << p_a_cm[1] << " " << p_a_cm[2];
  // LOG_S(ERROR) << "p_a_bar_cm: " << p_a_bar_cm[0] << " " << p_a_bar_cm[1] << " " << p_a_bar_cm[2];
  // LOG_S(ERROR) << "p_a_bar   : " << p_a_bar[0] << " " << p_a_bar[1] << " " << p_a_bar[2];
  // LOG_S(ERROR);

  p_b_bar = p_b_bar_cm;
  p_b_bar += v_cm * v_cm.dot(p_b_bar_cm) * ( ( gamma_cm - 1 ) / v_cm.length2() );
  p_b_bar += v_cm * gamma_cm * p0_b_cm;

  vector3d<double> v_a_bar;
  vector3d<double> v_b_bar;

  v_a_bar = p_a_bar;
  v_a_bar /= mass_a;
  v_a_bar *= phys::rel::lorenz_factor_inv(v_a_bar.length2());

  v_b_bar = p_b_bar;
  v_b_bar /= mass_b;
  v_b_bar *= phys::rel::lorenz_factor_inv(v_b_bar.length2());

  // LOG_S(ERROR) << "v_a -> v_b         : " << v_a[0] << " " << v_a[1] << " " << v_a[2] << "\t" << "\t" << v_b[0] << " " << v_b[1] << " " << v_b[2];
  // LOG_S(ERROR) << "v_a_bar -> v_b_bar : " << v_a_bar[0] << " " << v_a_bar[1] << " " << v_a_bar[2] << "\t" << "\t" << v_b_bar[0] << " " << v_b_bar[1] << " " << v_b_bar[2];

  // LOG_S(ERROR) << "p_a, p_b         : " << p_a.length() << ", " << p_b.length();
  // LOG_S(ERROR) << "p_a_bar, p_b_bar : " << p_a_bar.length() << ", " << p_b_bar.length();

  //// end of main calculation

  if ( !isnormal(v_a_bar[0])
       || !isnormal(v_a_bar[1])
       || !isnormal(v_a_bar[2])
       || !isnormal(v_b_bar[0])
       || !isnormal(v_b_bar[1])
       || !isnormal(v_b_bar[2])
       )
    {
    LOG_S(ERROR) << "Something went wrong during collide. Velocities a :"
     << L_coulomb << ","
     << p_rel.length() << ","
     << v_rel.length() << ","
     << coll_freq << ";"
     << sin_theta << ";" // / (gamma_cm * (cos_theta - beta/beta_cm)) << ";"
     << d_p[0] << ","
     << d_p[1] << ","
     << d_p[2] << ";"
     // << p_a_bar_cm[0] << ","
     // << p_a_bar_cm[1] << ","
     // << p_a_bar_cm[2] << ","
     // << p_b_bar_cm[0] << ","
     // << p_b_bar_cm[1] << ","
     // << p_b_bar_cm[2] << ";"
      ;
    }
  // 7.54897e+16,1,0,8.57476e-148,1;inf;1;-nan,-nan,-nan
  // charge_a, charge_b, density_lowest, L_coulomb, p_rel.length(), v_rel.length());
  // double variance_tg_theta_half = coll_freq * time->step;
  // double tg_theta = 2 * math::random::normal(lib::sq_rt(variance_tg_theta_half));

  // d_p[0] = p_cm_l * (
  //   cos_theta_r * cos_phi_r * sin_theta_cm * cos_phi_cm
  //   - cos_theta_r * sin_phi_r * sin_theta_cm * sin_phi_cm
  //   + sin_theta_r * sin_theta_cm) - p_cm[0];
  // d_p[1] = p_cm_l * (
  //   sin_phi_r * sin_theta_cm * cos_phi_cm
  //   + cos_phi_r * sin_theta_cm * sin_phi_cm) - p_cm[1];
  // d_p[2] = p_cm_l * (
  //   - sin_theta_r * cos_phi_r * sin_theta_cm * cos_phi_cm
  //   - sin_theta_r * sin_phi_r * sin_theta_cm * sin_phi_cm
  //   + cos_theta_r * cos_theta_cm) - p_cm[2];


  // set new velocity components
  if (swap)
  {
    P_VEL_R(pa) = p_b_bar[0];
    P_VEL_PHI(pa) = p_b_bar[1];
    P_VEL_Z(pa) = p_b_bar[2];

    P_VEL_R(pb) = p_a_bar[0];
    P_VEL_PHI(pb) = p_a_bar[1];
    P_VEL_Z(pb) = p_a_bar[2];
  }
  else
  {
    P_VEL_R(pa) = p_a_bar[0];
    P_VEL_PHI(pa) = p_a_bar[1];
    P_VEL_Z(pa) = p_a_bar[2];

    P_VEL_R(pb) = p_b_bar[0];
    P_VEL_PHI(pb) = p_b_bar[1];
    P_VEL_Z(pb) = p_b_bar[2];
  }
  // LOG_S(FATAL);
}

void CollisionsSentokuM::collide ()
{
  // pairing
  for (int i = 0; i < geometry->r_grid_amount; ++i)
    for (int j = 0; j < geometry->z_grid_amount; ++j)
    {
      unsigned int vec_size_ions = map_ion2cell(i, j).size();
      unsigned int vec_size_electrons = map_el2cell(i, j).size();

      // TA77: case 1a
      // ions
      if (vec_size_ions % 2 == 0)
        for (unsigned int k = 0; k < vec_size_ions; k = k + 2)
          collide_single(i, j, PROTON_MASS, PROTON_MASS,
                         (*map_ion2cell(i, j)[k]),
                         (*map_ion2cell(i, j)[k+1]));
      // electrons
      if (vec_size_electrons % 2 == 0)
        for (unsigned int k = 0; k < vec_size_electrons; k = k + 2)
          collide_single(i, j, EL_MASS, EL_MASS,
                         (*map_el2cell(i, j)[k]),
                         (*map_el2cell(i, j)[k+1]));

      // TA77: case 1b
      // ions
      if (vec_size_ions % 2 != 0)
      {
        if (vec_size_ions >= 3)
        {
          // first 3 collisions in special way
          collide_single(i, j, PROTON_MASS, PROTON_MASS,
                         (*map_ion2cell(i, j)[0]),
                         (*map_ion2cell(i, j)[1]));
          collide_single(i, j, PROTON_MASS, PROTON_MASS,
                         (*map_ion2cell(i, j)[1]),
                         (*map_ion2cell(i, j)[2]));
          collide_single(i, j, PROTON_MASS, PROTON_MASS,
                         (*map_ion2cell(i, j)[2]),
                         (*map_ion2cell(i, j)[0]));
        }
        if (vec_size_ions >= 5)
          for (unsigned int k = 3; k < vec_size_ions; k = k + 2)
            collide_single(i, j, PROTON_MASS, PROTON_MASS,
                           (*map_ion2cell(i, j)[k]),
                           (*map_ion2cell(i, j)[k+1]));
      }
      // TA77: electrons, case 1b
      // electrons
      if (vec_size_electrons % 2 != 0)
      {
        if (vec_size_electrons >= 3)
        {
          // first 3 collisions in special way
          collide_single(i, j, EL_MASS, EL_MASS,
                         (*map_el2cell(i, j)[0]),
                         (*map_el2cell(i, j)[1]));
          collide_single(i, j, EL_MASS, EL_MASS,
                         (*map_el2cell(i, j)[1]),
                         (*map_el2cell(i, j)[2]));
          collide_single(i, j, EL_MASS, EL_MASS,
                         (*map_el2cell(i, j)[2]),
                         (*map_el2cell(i, j)[0]));
        }
        if (vec_size_electrons >= 5)
          for (unsigned int k = 3; k < vec_size_electrons; k = k + 2)
            collide_single(i, j, EL_MASS, EL_MASS,
                           (*map_el2cell(i, j)[k]),
                           (*map_el2cell(i, j)[k+1]));
      }

      // TA77: case 2a. electrons-ions
      if (vec_size_ions == vec_size_electrons)
        for (unsigned int k = 0; k < vec_size_electrons; ++k)
          collide_single(i, j, EL_MASS, PROTON_MASS,
                         (*map_el2cell(i, j)[k]),
                         (*map_ion2cell(i, j)[k]));

      // TA77: case 2b. electrons-ions
      if (vec_size_ions > vec_size_electrons && vec_size_electrons > 0 && vec_size_ions > 0)
      {
        unsigned int c_i = floor( (float)vec_size_ions / (float)vec_size_electrons );
        double c_r = (float)vec_size_ions / (float)vec_size_electrons - c_i;

        int ions_1st_group = (c_i + 1) * c_r * vec_size_electrons;
        int els_1st_group = c_r * vec_size_electrons;

        int ions_2nd_group = c_i * (1 - c_r) * vec_size_electrons;
        // int els_2nd_group = (1 - c_r) * vec_size_electrons;

        // TA77: case 2b, 1st group, ions
        for (unsigned int fgi = 0; fgi < ions_1st_group; ++fgi)
        {
          unsigned int fge = floor(float(fgi) / float(c_i+1));
          collide_single(i, j, EL_MASS, PROTON_MASS,
                         (*map_el2cell(i, j)[fge]),
                         (*map_ion2cell(i, j)[fgi]));
        }
        // TA77: case 2b, 2nd group, ions
        for (unsigned int fgi = 0; fgi < ions_2nd_group; ++fgi)
        {
          unsigned int fge = floor(float(fgi) / float(c_i));
          collide_single(i, j, EL_MASS, PROTON_MASS,
                         (*map_el2cell(i, j)[fge+els_1st_group]),
                         (*map_ion2cell(i, j)[fgi+ions_1st_group]));
        }
      }
      // TA77: case 2b. electrons-ions
      if (vec_size_ions < vec_size_electrons && vec_size_electrons > 0 && vec_size_ions > 0)
      {
        double c_i = floor( (float)vec_size_electrons / (float)vec_size_ions );
        double c_r = (float)vec_size_electrons / (float)vec_size_ions - c_i;

        int els_1st_group = (c_i + 1) * c_r * vec_size_ions;
        int ions_1st_group = c_r * vec_size_ions;

        int els_2nd_group = c_i * (1 - c_r) * vec_size_ions;
        // int ions_2nd_group = (1 - c_r) * vec_size_ions;

        // TA77: case 2b, 1st group, electrons
        for (unsigned int fge = 0; fge < els_1st_group; ++fge)
        {
          unsigned int fgi = floor(float(fge) / float(c_i+1));
          // MSG("els 1st:  FGE " << fge << " FGI " << fgi << " VSI " << vec_size_ions << " VSE " << vec_size_electrons);
          collide_single(i, j, EL_MASS, PROTON_MASS,
                         (*map_el2cell(i, j)[fge]),
                         (*map_ion2cell(i, j)[fgi]));
        }
        // TA77: case 2b, 2nd group, electrons
        for (unsigned int fge = 0; fge < els_2nd_group; ++fge)
        {
          unsigned int fgi = floor(float(fge) / float(c_i));
          // MSG("els 2nd:  FGE " << fge << " FGI " << fgi << " VSI " << vec_size_ions << " VSE " << vec_size_electrons);
          collide_single(i, j, EL_MASS, PROTON_MASS,
                         (*map_el2cell(i, j)[fge+els_1st_group]),
                         (*map_ion2cell(i, j)[fgi+ions_1st_group]));
        }
      }
    }
}

// void CollisionsSentokuM::correct_velocities()
// {
//   // TITU18: correct velocities
//   for (int i = 0; i < geometry->r_grid_amount; ++i)
//     for (int j = 0; j < geometry->z_grid_amount; ++j)
//     {
//       unsigned int vec_size_ions = map_ion2cell(i, j).size();
//       unsigned int vec_size_electrons = map_el2cell(i, j).size();

//       // TITU18: calculate delta V
//       //// calculate summary moment components and total energy for t+delta_t
//       double moment_new_r_ion = 0, moment_new_phi_ion = 0, moment_new_z_ion = 0,
//         moment_new_r_el = 0, moment_new_phi_el = 0, moment_new_z_el = 0,
//         E_tot_ion_new = 0, E_tot_el_new = 0;

//       for (unsigned int k = 0; k < vec_size_ions; ++k)
//       {
//         double vr = P_VEL_R((*map_ion2cell(i, j)[k]));
//         double vphi = P_VEL_PHI((*map_ion2cell(i, j)[k]));
//         double vz = P_VEL_Z((*map_ion2cell(i, j)[k]));
//         double v_sq = vr*vr + vphi*vphi + vz*vz;

//         double mass = P_MASS((*map_ion2cell(i, j)[k]));
//         double weight = mass / PROTON_MASS;
//         double weighted_m = weight * mass;

//         moment_new_r_ion += weighted_m * vr;
//         moment_new_phi_ion += weighted_m * vphi;
//         moment_new_z_ion += weighted_m * vz;
//         E_tot_ion_new += weighted_m * v_sq / 2;
//       }

//       for (unsigned int k = 0; k < vec_size_electrons; ++k)
//       {
//         double vr = P_VEL_R((*map_el2cell(i, j)[k]));
//         double vphi = P_VEL_PHI((*map_el2cell(i, j)[k]));
//         double vz = P_VEL_Z((*map_el2cell(i, j)[k]));
//         double v_sq = vr*vr + vphi*vphi + vz*vz;

//         double mass = P_MASS((*map_el2cell(i, j)[k]));
//         double weight = mass / EL_MASS;
//         double weighted_m = weight * mass;

//         moment_new_r_el += weighted_m * vr;
//         moment_new_phi_el += weighted_m * vphi;
//         moment_new_z_el += weighted_m * vz;
//         E_tot_el_new += weighted_m * v_sq / 2;
//       }

//       //// calculate delta V components and delta E total
//       double delta_V_r_ion = (moment_new_r_ion - moment_tot_ion[0](i, j)) / mass_tot_ion(i, j);
//       double delta_V_phi_ion = (moment_new_phi_ion - moment_tot_ion[1](i, j)) / mass_tot_ion(i, j);
//       double delta_V_z_ion = (moment_new_z_ion - moment_tot_ion[2](i, j)) / mass_tot_ion(i, j);
//       double delta_V_r_el = (moment_new_r_el - moment_tot_el[0](i, j)) / mass_tot_el(i, j);
//       double delta_V_phi_el = (moment_new_phi_el - moment_tot_el[1](i, j)) / mass_tot_el(i, j);
//       double delta_V_z_el = (moment_new_z_el - moment_tot_el[2](i, j)) / mass_tot_el(i, j);

//       //// calculate V_0 components
//       double V_0_r_ion = moment_tot_ion[0](i, j) / mass_tot_ion(i, j);
//       double V_0_phi_ion = moment_tot_ion[1](i, j) / mass_tot_ion(i, j);
//       double V_0_z_ion = moment_tot_ion[2](i, j) / mass_tot_ion(i, j);
//       double V_0_r_el = moment_tot_el[0](i, j) / mass_tot_el(i, j);
//       double V_0_phi_el = moment_tot_el[1](i, j) / mass_tot_el(i, j);
//       double V_0_z_el = moment_tot_el[2](i, j) / mass_tot_el(i, j);
//       //// calculate delta E total
//       double delta_E_tot_ion = E_tot_ion_new - energy_tot_ion(i, j);
//       double delta_E_tot_el = E_tot_el_new - energy_tot_el(i, j);
//       //// link E total to local vars
//       double E_tot_ion = energy_tot_ion(i, j);
//       double E_tot_el = energy_tot_el(i, j);
//       double V_0_sq_ion = V_0_r_ion*V_0_r_ion + V_0_phi_ion*V_0_phi_ion + V_0_z_ion*V_0_z_ion;
//       double delta_V_sq_ion = delta_V_r_ion*delta_V_r_ion
//         + delta_V_phi_ion*delta_V_phi_ion
//         + delta_V_z_ion*delta_V_z_ion;
//       double V_0_sq_el = V_0_r_el*V_0_r_el + V_0_phi_el*V_0_phi_el + V_0_z_el*V_0_z_el;
//       double delta_V_sq_el = delta_V_r_el*delta_V_r_el
//         + delta_V_phi_el*delta_V_phi_el
//         + delta_V_z_el*delta_V_z_el;
//       /// calculate alpha
//       double alpha_ion =
//         ( E_tot_ion - mass_tot_ion(i, j) * V_0_sq_ion / 2 )
//         / ( E_tot_ion
//             + delta_E_tot_ion
//             - mass_tot_ion(i, j)
//             * ( pow(V_0_r_ion + delta_V_r_ion, 2)
//                 + pow(V_0_phi_ion + delta_V_phi_ion, 2)
//                 + pow(V_0_z_ion + delta_V_z_ion, 2)
//               )
//             / 2
//           );
//       double alpha_el =
//         ( E_tot_el - mass_tot_el(i, j) * V_0_sq_el / 2 )
//         / ( E_tot_el
//             + delta_E_tot_el
//             - mass_tot_el(i, j)
//             * ( pow(V_0_r_el + delta_V_r_el, 2)
//                 + pow(V_0_phi_el + delta_V_phi_el, 2)
//                 + pow(V_0_z_el + delta_V_z_el, 2)
//               )
//             / 2
//           );

//       // correct ion velocity
//       if (isnormal(alpha_ion)) // FIXME: sometimes E_tot_ion = mass_tot_ion(i, j) * V_0_sq_ion / 2
//       {
//         // FIXME: quick and dirty workadound
//         // caused negative values under squared root
//         // and zeros
//         if (alpha_ion < 0) alpha_ion = -alpha_ion;
//         if (alpha_ion == 0) alpha_ion = 1;
//         alpha_ion = lib::sq_rt(alpha_ion);

//         for (unsigned int k = 0; k < vec_size_ions; ++k)
//         {
//           double vr = P_VEL_R((*map_ion2cell(i, j)[k]));
//           double vphi = P_VEL_PHI((*map_ion2cell(i, j)[k]));
//           double vz = P_VEL_Z((*map_ion2cell(i, j)[k]));

//           double vr_corr = V_0_r_ion + alpha_ion * (vr - V_0_r_ion - delta_V_r_ion);
//           double vphi_corr = V_0_phi_ion + alpha_ion * (vphi - V_0_phi_ion - delta_V_phi_ion);
//           double vz_corr = V_0_z_ion + alpha_ion * (vz - V_0_z_ion - delta_V_z_ion);

//           P_VEL_R((*map_ion2cell(i, j)[k])) = vr_corr;
//           P_VEL_PHI((*map_ion2cell(i, j)[k])) = vphi_corr;
//           P_VEL_Z((*map_ion2cell(i, j)[k])) = vz_corr;
//         }
//       }

//       // correct electron velocity
//       if (isnormal(alpha_el)) // FIXME: sometimes E_tot_el = mass_tot_el(i, j) * V_0_sq_el / 2
//       {
//         // FIXME: quick and dirty workadound
//         // caused negative values under squared root
//         // and zeros
//         if (alpha_el < 0) alpha_el = -alpha_el;
//         if (alpha_el == 0) alpha_el = 1;
//         alpha_el = lib::sq_rt(alpha_el);

//         for (unsigned int k = 0; k < vec_size_electrons; ++k)
//         {
//           double vr = P_VEL_R((*map_el2cell(i, j)[k]));
//           double vphi = P_VEL_PHI((*map_el2cell(i, j)[k]));
//           double vz = P_VEL_Z((*map_el2cell(i, j)[k]));

//           double vr_corr = V_0_r_el + alpha_el * (vr - V_0_r_el - delta_V_r_el);
//           double vphi_corr = V_0_phi_el + alpha_el * (vphi - V_0_phi_el - delta_V_phi_el);
//           double vz_corr = V_0_z_el + alpha_el * (vz - V_0_z_el - delta_V_z_el);

//           P_VEL_R((*map_el2cell(i, j)[k])) = vr_corr;
//           P_VEL_PHI((*map_el2cell(i, j)[k])) = vphi_corr;
//           P_VEL_Z((*map_el2cell(i, j)[k])) = vz_corr;
//         }
//       }
//     }

// }

void CollisionsSentokuM::run ()
{
  clear();
  sort_to_cells();
  random_sort();
  collect_weighted_params_tot_grid();
  collide();
  // correct_velocities();
}

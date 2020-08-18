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

#include "collisionsTanaka.hpp"

// say: TA77   - Takizuka, Abe; 1977; DOI: 10.1016/0021-9991(77)90099-7
//      TITU18 - Tanaka et al.; 2018; DOI: 10.1002/ctpp.201700121

using namespace std;

CollisionsTanaka::CollisionsTanaka (Geometry* _geometry, TimeSim *_time, vector <SpecieP *> _species_p) : Collisions ( _geometry, _time, _species_p)
{
}

void CollisionsTanaka::collide_single(int i, int j, double m_real_a, double m_real_b,
                                vector<double> &pa, vector<double> &pb)
{
  // get required parameters
  double vr_a, vphi_a, vz_a, charge_a, mass_a, vr_b, vphi_b, vz_b,
    charge_b, mass_b;

  bool swap = false;

  // TITU18: find weight ratio
  double w_ratio = 1; // P_MASS(pa) * m_real_b / (P_MASS(pb) * m_real_a);

  // a-particle should be lighter, than b-particle
  if (w_ratio <= 1)
  {
    vr_a = P_VEL_R(pa);
    vphi_a = P_VEL_PHI(pa);
    vz_a = P_VEL_Z(pa);
    charge_a = P_CHARGE(pa);
    mass_a = P_MASS(pa);

    vr_b = P_VEL_R(pb);
    vphi_b = P_VEL_PHI(pb);
    vz_b = P_VEL_Z(pb);
    charge_b = P_CHARGE(pb);
    mass_b = P_MASS(pb);
  }
  else
  {
    vr_a = P_VEL_R(pb);
    vphi_a = P_VEL_PHI(pb);
    vz_a = P_VEL_Z(pb);
    charge_a = P_CHARGE(pb);
    mass_a = P_MASS(pb);

    vr_b = P_VEL_R(pa);
    vphi_b = P_VEL_PHI(pa);
    vz_b = P_VEL_Z(pa);
    charge_b = P_CHARGE(pa);
    mass_b = P_MASS(pa);

    // swap particles when b-particle is lighter, than a-particle
    w_ratio = 1 / w_ratio;
    swap = true;
  }

  double density_el = get_el_density(i, j);
  double density_ion = get_ion_density(i, j);
  double density_lowest = min(density_el, density_ion);
  double temperature_el = get_el_temperature(i, j);
  double m_ab = mass_a * mass_b / (mass_a + mass_b);

  // relative velocity
  double ux = vr_a - vr_b;
  double uy = vphi_a - vphi_b;
  double uz = vz_a - vz_b;
  double u = lib::sq_rt(pow(ux, 2) + pow(uy, 2) + pow(uz, 2));


  // if ``u'' (relative velocity) is zero, particles can not collide
  if (u == 0) return;
  // do not collide, if density or temperature is zero
  if (density_el == 0 || temperature_el == 0) return;

  // get lambda Coulomb
  double debye = phys::plasma::debye_length(density_el, temperature_el);
  double lambda_coulomb = phys::plasma::coulomb_logarithm (mass_a, mass_b, debye, u);

  // TA77: calculate u perpendicular
  double u_p = lib::sq_rt(pow(ux, 2) + pow(uy, 2));

  // TA77: calculate scattering angles Theta and Phi
  // TA77 sub: find delta
  //// find variance of delta
  double variance_d = pow(charge_a, 2) * pow(charge_b, 2) * density_lowest
    * lambda_coulomb / (8 * constant::PI * pow(EPSILON0, 2) * m_ab * pow(u, 3))
    * time->step;
  //// find standard deviation of delta
  double std_dev_d = lib::sq_rt(variance_d);
  //// find delta
  double delta = math::random::normal(std_dev_d);

  // TA77: find Theta angle
  double sin_Theta = 2 * delta / (1 + pow(delta, 2));
  double cos_Theta = 1 - 2 * pow(delta, 2) / (1 + pow(delta, 2));

  // TITU18: find Theta_2 angle
  double sin_Theta2 = lib::sq_rt(
    w_ratio * pow(sin_Theta, 2)
    + w_ratio * ( 1 - w_ratio ) * pow(( 1 - cos_Theta), 2));
  double cos_Theta2 = 1 - w_ratio * (1 - cos_Theta);
  // TITU18: sign of the sinus of Theta_angle
  // should be the same as sinus of Theta2_angle
  if (sin_Theta2 * sin_Theta < 0)
    sin_Theta2 = -sin_Theta2;

  // TA77: find Phi angle
  double Phi_angle = math::random::uniform_angle();
  double sin_Phi = sin(Phi_angle);
  double cos_Phi = cos(Phi_angle);

  // TITU18: find Phi_2 angle
  //// find dzeta for Phi_2 angle
  double dzeta = math::random::uniform();
  double Phi2_angle = Phi_angle + 2 * constant::PI * ( 1 - lib::sq_rt(w_ratio)) * (1 - dzeta);
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


  if ( !isnormal(vr_a_new)
       || !isnormal(vphi_a_new)
       || !isnormal(vz_a_new)
       || !isnormal(vr_b_new)
       || !isnormal(vphi_b_new)
       || !isnormal(vz_b_new)
    )
  {
    LOG_S(ERROR) << "Something went wrong during collide. Velocities a :"
            << vr_a_new << ","
            << vphi_a_new << ","
            << vz_a_new << ","
            << vr_b_new << ","
            << vphi_b_new << ","
            << vz_b_new << ";"
            << "Sinuses: "
            << sin_Theta << ","
            << sin_Theta2 << ","
            << sin_Phi << ","
            << sin_Phi2 << ";"
            << "Cosinuses: "
            << cos_Theta << ","
            << cos_Theta2 << ","
            << cos_Phi << ","
            << cos_Phi2 << ";"
            << "Variance of delta: "
            << m_ab << " "
            << pow(u, 3);
  }

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

void CollisionsTanaka::collide ()
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

void CollisionsTanaka::correct_velocities()
{
  // TITU18: correct velocities
  for (int i = 0; i < geometry->r_grid_amount; ++i)
    for (int j = 0; j < geometry->z_grid_amount; ++j)
    {
      unsigned int vec_size_ions = map_ion2cell(i, j).size();
      unsigned int vec_size_electrons = map_el2cell(i, j).size();

      // TITU18: calculate delta V
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

        double mass = P_MASS((*map_ion2cell(i, j)[k]));
        double weight = mass / PROTON_MASS;
        double weighted_m = weight * mass;

        moment_new_r_ion += weighted_m * vr;
        moment_new_phi_ion += weighted_m * vphi;
        moment_new_z_ion += weighted_m * vz;
        E_tot_ion_new += weighted_m * v_sq / 2;
      }

      for (unsigned int k = 0; k < vec_size_electrons; ++k)
      {
        double vr = P_VEL_R((*map_el2cell(i, j)[k]));
        double vphi = P_VEL_PHI((*map_el2cell(i, j)[k]));
        double vz = P_VEL_Z((*map_el2cell(i, j)[k]));
        double v_sq = vr*vr + vphi*vphi + vz*vz;

        double mass = P_MASS((*map_el2cell(i, j)[k]));
        double weight = mass / EL_MASS;
        double weighted_m = weight * mass;

        moment_new_r_el += weighted_m * vr;
        moment_new_phi_el += weighted_m * vphi;
        moment_new_z_el += weighted_m * vz;
        E_tot_el_new += weighted_m * v_sq / 2;
      }

      //// calculate delta V components and delta E total
      double delta_V_r_ion = (moment_new_r_ion - moment_tot_ion[0](i, j)) / mass_tot_ion(i, j);
      double delta_V_phi_ion = (moment_new_phi_ion - moment_tot_ion[1](i, j)) / mass_tot_ion(i, j);
      double delta_V_z_ion = (moment_new_z_ion - moment_tot_ion[2](i, j)) / mass_tot_ion(i, j);
      double delta_V_r_el = (moment_new_r_el - moment_tot_el[0](i, j)) / mass_tot_el(i, j);
      double delta_V_phi_el = (moment_new_phi_el - moment_tot_el[1](i, j)) / mass_tot_el(i, j);
      double delta_V_z_el = (moment_new_z_el - moment_tot_el[2](i, j)) / mass_tot_el(i, j);

      //// calculate V_0 components
      double V_0_r_ion = moment_tot_ion[0](i, j) / mass_tot_ion(i, j);
      double V_0_phi_ion = moment_tot_ion[1](i, j) / mass_tot_ion(i, j);
      double V_0_z_ion = moment_tot_ion[2](i, j) / mass_tot_ion(i, j);
      double V_0_r_el = moment_tot_el[0](i, j) / mass_tot_el(i, j);
      double V_0_phi_el = moment_tot_el[1](i, j) / mass_tot_el(i, j);
      double V_0_z_el = moment_tot_el[2](i, j) / mass_tot_el(i, j);
      //// calculate delta E total
      double delta_E_tot_ion = E_tot_ion_new - energy_tot_ion(i, j);
      double delta_E_tot_el = E_tot_el_new - energy_tot_el(i, j);
      //// link E total to local vars
      double E_tot_ion = energy_tot_ion(i, j);
      double E_tot_el = energy_tot_el(i, j);
      double V_0_sq_ion = V_0_r_ion*V_0_r_ion + V_0_phi_ion*V_0_phi_ion + V_0_z_ion*V_0_z_ion;
      double delta_V_sq_ion = delta_V_r_ion*delta_V_r_ion
        + delta_V_phi_ion*delta_V_phi_ion
        + delta_V_z_ion*delta_V_z_ion;
      double V_0_sq_el = V_0_r_el*V_0_r_el + V_0_phi_el*V_0_phi_el + V_0_z_el*V_0_z_el;
      double delta_V_sq_el = delta_V_r_el*delta_V_r_el
        + delta_V_phi_el*delta_V_phi_el
        + delta_V_z_el*delta_V_z_el;
      /// calculate alpha
      double alpha_ion =
        ( E_tot_ion - mass_tot_ion(i, j) * V_0_sq_ion / 2 )
        / ( E_tot_ion
            + delta_E_tot_ion
            - mass_tot_ion(i, j)
            * ( pow(V_0_r_ion + delta_V_r_ion, 2)
                + pow(V_0_phi_ion + delta_V_phi_ion, 2)
                + pow(V_0_z_ion + delta_V_z_ion, 2)
              )
            / 2
          );
      double alpha_el =
        ( E_tot_el - mass_tot_el(i, j) * V_0_sq_el / 2 )
        / ( E_tot_el
            + delta_E_tot_el
            - mass_tot_el(i, j)
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
        alpha_ion = lib::sq_rt(alpha_ion);

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
        alpha_el = lib::sq_rt(alpha_el);

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

void CollisionsTanaka::run ()
{
  clear();
  sort_to_cells();
  random_sort();
  collect_weighted_params_tot_grid();
  collide();
  // correct_velocities();
}

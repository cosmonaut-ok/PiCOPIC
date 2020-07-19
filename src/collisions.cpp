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

#include "collisions.hpp"

#include <vector>
#include <algorithm>    // std::min, std::random_shuffle
#include <typeinfo>
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <math.h>       // floor, asin

// #include <random>
// #include <iostream>
// #include <stdexcept>
// #include <algorithm>

// using std::mt19937_64;
// using std::random_device;
// using std::uniform_int_distribution;
// using std::vector;
// using std::cout;
// using std::cerr;
// using std::endl;
// using std::out_of_range;
#include <string>

#include "lib.hpp"
#include "geometry.hpp"
#include "specieP.hpp"

using namespace std;

Collisions::Collisions (Geometry* _geometry, TimeSim *_time, vector <SpecieP *> _species_p):geometry(_geometry),time(_time)
{
  species_p = _species_p;
  geometry = _geometry;

  //! TA77: make dummies to place particles there
  map_el2cell = Grid<vector< vector<double> * >> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
  map_ion2cell = Grid<vector< vector<double> * >> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
}

//! TA77: clear
void Collisions::clear()
{
  for (int i = 0; i < geometry->r_grid_amount; ++i)
    for (int j = 0; j < geometry->z_grid_amount; ++j)
    {
      map_el2cell(i, j).clear();
      map_ion2cell(i, j).clear();
    }
}

//! TA77: 5. group particles to cells
void Collisions::sort_to_cells()
{
  for (auto ps = species_p.begin(); ps != species_p.end(); ++ps)
    for (auto i = (**ps).particles.begin(); i != (**ps).particles.end(); ++i)
    {
      double dr = geometry->r_cell_size;
      double dz = geometry->z_cell_size;

      // finding number new and old cells
      int i_n = CELL_NUMBER(P_POS_R((**i)), dr);
      int k_n = CELL_NUMBER(P_POS_Z((**i)), dz);

      //! shift also to take overlaying into account
      int i_n_shift = i_n - geometry->bottom_r_grid_number;
      int k_n_shift = k_n - geometry->left_z_grid_number;

      // TODO: push_back particle to sorter
      string name = "ions";
      string cname = "Ions";
      if (name.compare((**ps).name) == 0 || cname.compare((**ps).name) == 0)
        map_ion2cell(i_n_shift, k_n_shift).push_back((*i));
      else
        map_el2cell(i_n_shift, k_n_shift).push_back((*i));
    }
}

//! TA77: 6. resort cell elements randomly to simply prepare the pairs
void Collisions::random_sort ()
{
  for (int i = 0; i < geometry->r_grid_amount; ++i)
    for (int j = 0; j < geometry->z_grid_amount; ++j)
    {
      std::srand ( unsigned ( std::time(0) ) );
      std::random_shuffle(map_el2cell(i, j).begin(), map_el2cell(i, j).end());
      std::random_shuffle(map_ion2cell(i, j).begin(), map_ion2cell(i, j).end());
    }
}

void Collisions::collide_single(int i, int j, vector<double> &pa, vector<double> &pb)
{
  // get required parameters
  // ion - b; electron - a
  double vr_a = P_VEL_R(pa);
  double vphi_a = P_VEL_PHI(pa);
  double vz_a = P_VEL_Z(pa);
  double charge_a = P_CHARGE(pa);
  double mass_a = P_MASS(pa);

  double vr_b = P_VEL_R(pb);
  double vphi_b = P_VEL_PHI(pb);
  double vz_b = P_VEL_Z(pb);
  double charge_b = P_CHARGE(pb);
  double mass_b = P_MASS(pb);

  double density_el = get_el_density(i, j);
  double density_ion = get_ion_density(i, j);
  double density_lowest = min(density_el, density_ion);
  double lambda_coulomb = 10; // TODO: clarify and/or set it as a constant
  double m_ab = mass_a * mass_b / (mass_a + mass_b);

  // relative velocity
  double ux = vr_a - vr_b;
  double uy = vphi_a - vphi_b;
  double uz = vz_a - vz_b;
  double u = lib::sq_rt(pow(ux, 2) + pow(uy, 2) + pow(uz, 2));

  // TA77: calculate u perpendicular
  double u_p = lib::sq_rt(pow(ux, 2) + pow(uy, 2));

  // TA77: calculate scattering angles Theta and Phi
  double n_lowest = 1e17;
  double lambda = 10;
  // TA77 sub: find delta
  // find variance of delta
  double variance_d = pow(charge_a, 2) * pow(charge_b, 2) * n_lowest * lambda
    / (8 * constant::PI * pow(EPSILON0, 2) * m_ab * pow(u, 3));
  // find standard deviation of delta
  double std_dev_d = lib::sq_rt(variance_d);
  // find delta
  double delta = math::random::normal(std_dev_d);

  double Theta_angle = 2 * atan(delta);
  double Phi_angle = math::random::uniform_angle();

  // MSG(Theta_angle << " " << Phi_angle);

  // TA77: calculate delta u components
  double d_ux, d_uy, d_uz;
  double sin_Theta = sin(Theta_angle);
  double cos_Theta = cos(Theta_angle);
  double sin_Phi = sin(Phi_angle);
  double cos_Phi = cos(Phi_angle);
  // simplify calculations
  if (u_p == 0)
  {
    d_ux = u * sin_Theta * cos_Phi;
    d_uy = u * sin_Theta * sin_Phi;
    d_uz = -u * (1 - cos_Theta);
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
  }

  double vr_a_new = vr_a + m_ab/mass_a * d_ux;
  double vphi_a_new = vphi_a + m_ab/mass_a * d_uy;
  double vz_a_new = vz_a + m_ab/mass_a * d_uz;

  double vr_b_new = vr_b + m_ab/mass_b * d_ux;
  double vphi_b_new = vphi_b + m_ab/mass_b * d_uy;
  double vz_b_new = vz_b + m_ab/mass_b * d_uz;

  // MSG("R_a:   " << vr_a << "=>" << vr_a_new);
  // MSG("PHI_a: " << vphi_a << "=>" << vphi_a_new);
  // MSG("Z_a:   " << vz_a << "=>" << vz_a_new);
  // MSG("R_b:   " << vr_b << "=>" << vr_b_new);
  // MSG("PHI_b: " << vphi_b << "=>" << vphi_b_new);
  // MSG("Z_b:   " << vz_b << "=>" << vz_b_new);

  // exit(1);

  // TA77: set new velocity components
  P_VEL_R(pa) = vr_a_new;
  P_VEL_PHI(pa) = vphi_a_new;
  P_VEL_Z(pa) = vz_a_new;

  P_VEL_R(pb) = vr_b_new;
  P_VEL_PHI(pb) = vphi_b_new;
  P_VEL_Z(pb) = vz_b_new;

  if (vr_a_new > LIGHT_VEL || vphi_a_new > LIGHT_VEL || vz_a_new > LIGHT_VEL
      || vr_b_new > LIGHT_VEL || vphi_b_new > LIGHT_VEL || vz_b_new > LIGHT_VEL)
    MSG("OH SHI! " << vr_a << " " << vphi_a << " " << vz_a << " => "
        << vr_a_new << " " << vphi_a_new << " " << vz_a_new << " -- "
        << charge_a << " " << mass_a << "; "
        << vr_b << " " << vphi_b << " " << vz_b << " => "
        << vr_b_new << " " << vphi_b_new << " " << vz_b_new << " -- "
        << charge_b << " " << mass_b);
}

void Collisions::collide ()
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
          collide_single(i, j,
                         (*map_ion2cell(i, j)[k]),
                         (*map_ion2cell(i, j)[k+1]));
      // electrons
      if (vec_size_electrons % 2 == 0)
        for (unsigned int k = 0; k < vec_size_electrons; k = k + 2)
          collide_single(i, j,
                         (*map_el2cell(i, j)[k]),
                         (*map_el2cell(i, j)[k+1]));

      // TA77: case 1b
      // ions
      if (vec_size_ions % 2 != 0)
      {
        if (vec_size_ions >= 3)
        {
          // first 3 collisions in special way
          collide_single(i, j,
                         (*map_ion2cell(i, j)[0]),
                         (*map_ion2cell(i, j)[1]));
          collide_single(i, j,
                         (*map_ion2cell(i, j)[1]),
                         (*map_ion2cell(i, j)[2]));
          collide_single(i, j,
                         (*map_ion2cell(i, j)[2]),
                         (*map_ion2cell(i, j)[0]));
        }
        if (vec_size_ions >= 5)
          for (unsigned int k = 3; k < vec_size_ions; k = k + 2)
            collide_single(i, j,
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
          collide_single(i, j,
                         (*map_el2cell(i, j)[0]),
                         (*map_el2cell(i, j)[1]));
          collide_single(i, j,
                         (*map_el2cell(i, j)[1]),
                         (*map_el2cell(i, j)[2]));
          collide_single(i, j,
                         (*map_el2cell(i, j)[2]),
                         (*map_el2cell(i, j)[0]));
        }
        if (vec_size_electrons >= 5)
          for (unsigned int k = 3; k < vec_size_electrons; k = k + 2)
            collide_single(i, j,
                           (*map_el2cell(i, j)[k]),
                           (*map_el2cell(i, j)[k+1]));
      }

      // TA77: case 2a. electrons-ions
      if (vec_size_ions == vec_size_electrons)
        for (unsigned int k = 0; k < vec_size_electrons; ++k)
          collide_single(i, j,
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
        int els_2nd_group = (1 - c_r) * vec_size_electrons;

        // TA77: case 2b, 1st group, ions
        for (unsigned int fgi = 0; fgi < ions_1st_group; ++fgi)
        {
          unsigned int fge = floor(float(fgi) / float(c_i+1));
          collide_single(i, j,
                         (*map_el2cell(i, j)[fge]),
                         (*map_ion2cell(i, j)[fgi]));
        }
        // TA77: case 2b, 2nd group, ions
        for (unsigned int fgi = 0; fgi < ions_2nd_group; ++fgi)
        {
          unsigned int fge = floor(float(fgi) / float(c_i));
          collide_single(i, j,
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
        int ions_2nd_group = (1 - c_r) * vec_size_ions;

        // TA77: case 2b, 1st group, electrons
        for (unsigned int fge = 0; fge < els_1st_group; ++fge)
        {
          unsigned int fgi = floor(float(fge) / float(c_i+1));
          // MSG("els 1st:  FGE " << fge << " FGI " << fgi << " VSI " << vec_size_ions << " VSE " << vec_size_electrons);
          collide_single(i, j,
                         (*map_el2cell(i, j)[fge]),
                         (*map_ion2cell(i, j)[fgi]));
        }
        // TA77: case 2b, 2nd group, electrons
        for (unsigned int fge = 0; fge < els_2nd_group; ++fge)
        {
          unsigned int fgi = floor(float(fge) / float(c_i));
          // MSG("els 2nd:  FGE " << fge << " FGI " << fgi << " VSI " << vec_size_ions << " VSE " << vec_size_electrons);
          collide_single(i, j,
                         (*map_el2cell(i, j)[fge+els_1st_group]),
                         (*map_ion2cell(i, j)[fgi+ions_1st_group]));
        }
      }
    }
}

double Collisions::get_el_density(int i, int j)
{
  int len = map_el2cell(i, j).size();
  double sum_mass;
  double cell_volume = geometry_cell_volume(i);

  // summary electron density in the cell
  for (unsigned int p = 0; p < len; ++p)
    sum_mass += P_MASS((*map_el2cell(i, j)[p]));

  return sum_mass / cell_volume;
}

double Collisions::get_ion_density(int i, int j)
{
  int len = map_ion2cell(i, j).size();
  double sum_mass;
  double cell_volume = geometry_cell_volume(i);

  // summary electron density in the cell
  for (unsigned int p = 0; p < len; ++p)
    sum_mass += P_MASS((*map_ion2cell(i, j)[p]));

  return sum_mass / cell_volume;
}

double Collisions::geometry_cell_volume(int i)
{
  double dr = geometry->r_cell_size;
  double dz = geometry->z_cell_size;
  double shift = geometry->bottom_r_grid_number;

  return CELL_VOLUME(i+shift, dr, dz);

}

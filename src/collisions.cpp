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

#include "collisions.hpp"

using namespace std;
using namespace constant;

Collisions::Collisions (Geometry* _geometry, TimeSim *_time, vector <SpecieP *> _species_p):geometry(_geometry),time(_time)
{
  species_p = _species_p;

  // find electron and ion masses and charges
  string name_el = "electrons";
  string cname_el = "Electrons";

  string name_ion = "ions";
  string cname_ion = "Ions";

  for (auto ps = species_p.begin(); ps != species_p.end(); ++ps)
  {
    if (name_el.compare((**ps).name) == 0 || cname_el.compare((**ps).name) == 0)
    {
      specie_el = *ps;
      mass_el = (**ps).mass;
      charge_el = (**ps).charge;
    }
    if (name_ion.compare((**ps).name) == 0 || cname_ion.compare((**ps).name) == 0)
    {
      specie_ion = *ps;
      mass_ion = (**ps).mass;
      charge_ion = (**ps).charge;
    }
  }
  map_el2cell = Grid<vector< Particle * >> (geometry->cell_amount[0], geometry->cell_amount[1], 2);
  map_ion2cell = Grid<vector< Particle * >> (geometry->cell_amount[0], geometry->cell_amount[1], 2);

  energy_tot_el = Grid<double> (geometry->cell_amount[0], geometry->cell_amount[1], 2);
  amount_tot_el = Grid<double> (geometry->cell_amount[0], geometry->cell_amount[1], 2);
  moment_tot_el = Grid3D<double> (geometry->cell_amount[0], geometry->cell_amount[1], 2);

  energy_tot_ion = Grid<double> (geometry->cell_amount[0], geometry->cell_amount[1], 2);
  amount_tot_ion = Grid<double> (geometry->cell_amount[0], geometry->cell_amount[1], 2);
  moment_tot_ion = Grid3D<double> (geometry->cell_amount[0], geometry->cell_amount[1], 2);
}

//! TA77: clear
void Collisions::clear()
{
  energy_tot_el = 0;
  energy_tot_el.overlay_set(0);
  amount_tot_el = 0;
  amount_tot_el.overlay_set(0);
  moment_tot_el = 0;
  moment_tot_el.overlay_set(0);

  energy_tot_ion = 0;
  energy_tot_ion.overlay_set(0);
  amount_tot_ion = 0;
  amount_tot_ion.overlay_set(0);
  moment_tot_ion = 0;
  moment_tot_ion.overlay_set(0);

  for (int i = 0; i < geometry->cell_amount[0]; ++i)
    for (int j = 0; j < geometry->cell_amount[1]; ++j)
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
      // finding number new and old cells
      int i_n = P_CELL_R((**i));
      int k_n = P_CELL_Z((**i));

      //! shift also to take overlaying into account
      int i_n_shift = i_n - geometry->cell_dims[0];
      int k_n_shift = k_n - geometry->cell_dims[1];

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
  for (int i = 0; i < geometry->cell_amount[0]; ++i)
    for (int j = 0; j < geometry->cell_amount[1]; ++j)
    {
      std::srand ( unsigned ( std::time(0) ) );
      std::random_shuffle(map_el2cell(i, j).begin(), map_el2cell(i, j).end());
      std::random_shuffle(map_ion2cell(i, j).begin(), map_ion2cell(i, j).end());
    }
}

double Collisions::get_el_density(int i, int j)
{
  unsigned int len = map_el2cell(i, j).size();
  double sum_amount = 0;
  double cell_volume = geometry_cell_volume(i);

  // summary electron density in the cell
  for (unsigned int p = 0; p < len; ++p)
    sum_amount += P_WEIGHT((*map_el2cell(i, j)[p]));

  return sum_amount / cell_volume;
}

double Collisions::get_ion_density(int i, int j)
{
  unsigned int len = map_ion2cell(i, j).size();
  double sum_amount = 0;
  double cell_volume = geometry_cell_volume(i);

  // summary electron density in the cell
  for (unsigned int p = 0; p < len; ++p)
    sum_amount += P_WEIGHT((*map_ion2cell(i, j)[p]));

  return sum_amount / cell_volume;
}

double Collisions::get_el_temperature(int i, int j)
{
  double e_tot = energy_tot_el(i, j);
  double amount = amount_tot_el(i, j);
  return e_tot / amount * EL_CHARGE_INV; // also, convert from J to eV
}

double Collisions::get_ion_temperature(int i, int j)
{
  double e_tot = energy_tot_ion(i, j);
  double amount = amount_tot_ion(i, j);

  return e_tot / amount * EL_CHARGE_INV;
}

double Collisions::geometry_cell_volume(int i)
{
  double dr = geometry->cell_size[0];
  double dz = geometry->cell_size[1];
  double shift = geometry->cell_dims[0];

  return CELL_VOLUME(i+shift, dr, dz);
}

void Collisions::collect_weighted_params_tot_grid ()
{
  for (int i = 0; i < geometry->cell_amount[0]; ++i)
    for (int j = 0; j < geometry->cell_amount[1]; ++j)
    {
      unsigned int vec_size_ions = map_ion2cell(i, j).size();
      unsigned int vec_size_electrons = map_el2cell(i, j).size();

      // ion weighting
      for (unsigned int k = 0; k < vec_size_ions; ++k)
      {
        double vr = P_VEL_R((*map_ion2cell(i, j)[k]));
        double vphi = P_VEL_PHI((*map_ion2cell(i, j)[k]));
        double vz = P_VEL_Z((*map_ion2cell(i, j)[k]));
        double v_sq = vr*vr + vphi*vphi + vz*vz;

        double weight = P_WEIGHT((*map_ion2cell(i, j)[k]));
        double weighted_m = weight * mass_ion;

        // increase sum of moment
        moment_tot_ion[0].inc(i, j, weighted_m * vr);
        moment_tot_ion[1].inc(i, j, weighted_m * vphi);
        moment_tot_ion[2].inc(i, j, weighted_m * vz);
        // increase sum os mass
        amount_tot_ion.inc(i, j, weight);
        // increase sum of moment
        energy_tot_ion.inc(i, j, weighted_m * v_sq / 2);
      }

      // electron weighting
      for (unsigned int k = 0; k < vec_size_electrons; ++k)
      {
        double vr = P_VEL_R((*map_el2cell(i, j)[k]));
        double vphi = P_VEL_PHI((*map_el2cell(i, j)[k]));
        double vz = P_VEL_Z((*map_el2cell(i, j)[k]));
        double v_sq = vr*vr + vphi*vphi + vz*vz;

        double weight = P_WEIGHT((*map_el2cell(i, j)[k]));
        double weighted_m = weight * mass_el;

        // increase sum of moment
        moment_tot_el[0].inc(i, j, weighted_m * vr);
        moment_tot_el[1].inc(i, j, weighted_m * vphi);
        moment_tot_el[2].inc(i, j, weighted_m * vz);
        // increase sum os mass
        amount_tot_el.inc(i, j, weight);
        // increase sum of moment
        energy_tot_el.inc(i, j, weighted_m * v_sq / 2);
      }
    }
  // calculate temperature of electrons and ions
  specie_el->calc_temperature();
  specie_ion->calc_temperature();
  // temperature_el("electrons");
  // temperature_ion("ions");
}

void Collisions::operator()()
{
  clear();
  sort_to_cells();
  random_sort();
  collect_weighted_params_tot_grid();
  collide();
}

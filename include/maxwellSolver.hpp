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

#ifndef _MAXWELLSOLVER_HPP_
#define _MAXWELLSOLVER_HPP_

#include "defines.hpp"
#include "msg.hpp"

#include "math/vector3d.hpp"

#include "current.hpp"
#include "specieP.hpp"

using namespace std;

class SpecieP;
class Current;

class MaxwellSolver
{

public:
  Grid3D<double> field_e;
  Grid3D<double> field_h;

protected:
  Geometry *geometry;
  TimeSim *time;
  vector<SpecieP *> species_p;

  Grid<double> epsilon;
  Grid<double> sigma; // for PML
  Current *current;

public:
  MaxwellSolver ( void ) {};
  MaxwellSolver ( Geometry *_geometry, TimeSim *_time,
                  vector<SpecieP *> _species_p, Current* _current)
    : geometry(_geometry), time(_time), current(_current)
  {
    // attach pointer to vector of particle species
    species_p = _species_p;

    // initialize fields to zero-state
    field_e = Grid3D<double> (geometry->cell_amount[0], geometry->cell_amount[1], 2);
    field_h = Grid3D<double> (geometry->cell_amount[0], geometry->cell_amount[1], 2);
    field_e = 0.;
    field_h = 0.;
    field_e.overlay_set(0.);
    field_h.overlay_set(0.);

    // initialize epsilon and sigma (for PML)
    epsilon = Grid<double> (geometry->cell_amount[0], geometry->cell_amount[1], 2);
    epsilon = constant::EPSILON0;
    epsilon.overlay_set(constant::EPSILON0);

    sigma = Grid<double> (geometry->cell_amount[0], geometry->cell_amount[1], 2);
    sigma = 0.;
    sigma.overlay_set(0.);
  };

  ~MaxwellSolver(void) {};

  vector3d<double> get_field_dummy(__attribute__((unused)) double radius, __attribute__((unused)) double longitude)
  {
    // just a dummy method for test/debug
    vector3d<double> cmp; // field components
    return cmp;
  };

  virtual void set_pml() = 0;
  virtual void calc_field_h() = 0;
  virtual void calc_field_e() = 0;
  virtual vector3d<double> get_field_h(double radius, double longitude) = 0;
  virtual vector3d<double> get_field_e(double radius, double longitude) = 0;
};

#endif // end of _MAXWELLSOLVER_HPP_

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

#ifndef _MAXWELLSOLVERYEE_HPP_
#define _MAXWELLSOLVERYEE_HPP_

#include "maxwellSolver.hpp"

class MaxwellSolverYee : public MaxwellSolver
{

public:
  Grid3D<double> field_h_at_et;

private:
  // internal variables to calculation for dielectric walls
  unsigned int r_begin;
  unsigned int z_begin;
  unsigned int r_end;
  unsigned int z_end;

public:
  MaxwellSolverYee ( void ) {};
  MaxwellSolverYee ( Geometry *_geometry, TimeSim *_time,
                     vector<SpecieP *> _species_p, Current *_current);

  ~MaxwellSolverYee(void) {};

  void set_pml();
  void calc_field_h();
  void calc_field_e();
  vector3d<double> get_field_h(double radius, double longitude);
  vector3d<double> get_field_e(double radius, double longitude);
};

#endif // end of _MAXWELLSOLVERYEE_HPP_

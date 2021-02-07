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

#ifndef _EXTERNALFIELDS_HPP_
#define _EXTERNALFIELDS_HPP_

#include "maxwellSolver.hpp"

class ExternalFields : public MaxwellSolver
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
  ExternalFields ( void ) {};
  ExternalFields ( Geometry *_geometry, TimeSim *_time,
                   Current *_current );

  ~ExternalFields(void) {};

  void set_pml();
  void calc_field_h();
  void calc_field_e();
};

#endif // end of _EXTERNALFIELDS_HPP_

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

#ifndef _SMB_HPP_
#define _SMB_HPP_

#ifdef _OPENMP
#include <omp.h>
// #else
// #define omp_get_thread_num() 0
#endif

#include "defines.hpp"
#include "constant.hpp"
#include "msg.hpp"
#include "algo/grid.hpp"
#include "domain.hpp"

#include "cfg.hpp"
#include "timeSim.hpp"

using namespace std;

#define BEAM_ID_START 1000

class SMB
{
public:
  Grid<Domain*> domains;
  unsigned int r_domains;
  unsigned int z_domains;

private:
  Geometry *geometry;
  Cfg *cfg;
  TimeSim *time;

public:
  SMB ( void ) {};
  SMB ( Cfg* _cfg, Geometry *_geometry, TimeSim *_time );

  ~SMB( void ) {};

// private:
  void particles_runaway_collector ();
  void current_overlay ();
  void field_h_overlay ();
  void field_e_overlay ();

  // public:
  void solve_maxvell();
  void solve_current();
  void advance_particles();
  void inject_beam();
  void distribute();
};
#endif // end of _SMB_HPP_

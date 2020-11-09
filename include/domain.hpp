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

#ifndef _DOMAIN_HPP_
#define _DOMAIN_HPP_

#include <vector>

#include <algorithm>

#include "algo/grid.hpp"
#include "algo/grid3d.hpp"
#include "geometry.hpp"
#include "timeSim.hpp"
#include "cfg.hpp"
#include "msg.hpp"

#include "specieP.hpp"

#if defined (PUSHER_BORIS_ADAPTIVE) || defined (PUSHER_BORIS_CLASSIC) || defined (PUSHER_BORIS_RELATIVISTIC)
#include "pusher/pusherBoris.hpp"
#elif PUSHER_VAY
#include "pusher/pusherVay.hpp"
#elif PUSHER_HIGUERA_CARY
#include "pusher/pusherHC.hpp"
#endif

#ifdef MAXWELL_SOLVER_YEE
#include "maxwellSolver/maxwellSolverYee.hpp"
#endif

#ifdef CCS_VILLASENOR_BUNEMAN
#include "current/currentVB.hpp"
#elif CCS_ZIGZAG
#include "current/currentZigZag.hpp"
#endif

#include "density/densityCharge.hpp"

#ifdef DENSITY_CALC_COUNTING
#include "density/densityCounted.hpp"
#elif DENSITY_CALC_WEIGHTING
#include "density/densityWeighted.hpp"
#endif

#ifdef TEMP_CALC_COUNTING
#include "temperature/temperatureCounted.hpp"
#elif TEMP_CALC_WEIGHTING
#include "temperature/temperatureWeighted.hpp"
#endif

#ifdef COLLISIONS
#ifdef COULOMB_COLLISIONS_TA77S
#include "collisions/collisionsTA77S.hpp"
#elif COULOMB_COLLISIONS_SK98
#include "collisions/collisionsSK98.hpp"
#elif COULOMB_COLLISIONS_P12
#include "collisions/collisionsP12.hpp"
#endif
#endif

using namespace std;

class Domain
{
public:
  vector<SpecieP *> species_p;
  Current *current;

#if defined (PUSHER_BORIS_ADAPTIVE) || defined (PUSHER_BORIS_CLASSIC) || defined (PUSHER_BORIS_RELATIVISTIC)
  PusherBoris *pusher;
#elif PUSHER_VAY
  PusherVay *pusher;
#elif PUSHER_HIGUERA_CARY
  PusherHC *pusher;
#endif

#ifdef MAXWELL_SOLVER_YEE
  MaxwellSolverYee *maxwell_solver;
#endif

  Temperature *temperature;
  Density *density;
  DensityCharge *charge;
#ifdef COLLISIONS
  Collisions *collisions;
#endif

  Geometry geometry;

  TimeSim *time; // simulation time object

public:
  Domain() {};
  Domain(Geometry _geometry, vector<SpecieP *> species_p, TimeSim* _time);

  // wrapper methods
  void distribute();
  void weight_density(string specie);
  void weight_temperature(string specie);
  void weight_charge(string specie);
  void push_particles();
  void weight_current();
  void update_particles_coords();
  // void weight_current_azimuthal();
  void reset_current();
  void reset_charge();
  void reset_field_e() {};
  void reset_field_h() {};
  void weight_field_h();
  void weight_field_e();
  void particles_back_velocity_to_rz();
  void particles_back_position_to_rz();
  void reflect();
  void manage_beam();
  void dump_particle_positions_to_old();
  void collide();
  void bind_cell_numbers();
};
#endif // end of _DOMAIN_HPP_

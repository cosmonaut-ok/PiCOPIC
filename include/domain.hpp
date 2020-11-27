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

#include "defines.hpp"
#include "msg.hpp"

#include "specieP.hpp"
#include "beamP.hpp"

#if defined(SWITCH_PUSHER_BORIS_ADAPTIVE) || defined(SWITCH_PUSHER_BORIS) || defined(SWITCH_PUSHER_BORIS_RELATIVISTIC)
#include "pusher/pusherBoris.hpp"
#elif defined(SWITCH_PUSHER_VAY
#include "pusher/pusherVay.hpp"
#elif defined(SWITCH_PUSHER_HC
#include "pusher/pusherHC.hpp"
#endif

#ifdef SWITCH_MAXWELL_SOLVER_YEE
#include "maxwellSolver/maxwellSolverYee.hpp"
#endif

#ifdef SWITCH_CCS_VB
#include "current/currentVB.hpp"
#elif defined(SWITCH_CCS_ZIGZAG)
#include "current/currentZigZag.hpp"
#endif

// #ifdef SWITCH_DENSITY_CALC_COUNTING
// #include "density/densityCounted.hpp"
// #elif defined(SWITCH_DENSITY_CALC_WEIGHTING)
// #include "density/densityWeighted.hpp"
// #endif

// #ifdef SWITCH_TEMP_CALC_COUNTING
// #include "temperature/temperatureCounted.hpp"
// #elif defined(WITH_TEMP_CALC_WEIGHTING)
// #include "temperature/temperatureWeighted.hpp"
// #endif

#ifdef ENABLE_COULOMB_COLLISIONS
#ifdef SWITCH_COULOMB_COLLISIONS_TA77S
#include "collisions/collisionsTA77S.hpp"
#elif defined(SWITCH_COULOMB_COLLISIONS_SK98)
#include "collisions/collisionsSK98.hpp"
#elif defined(SWITCH_COULOMB_COLLISIONS_P12
#include "collisions/collisionsP12.hpp"
#endif
#endif

// #include "density/densityCharge.hpp"

using namespace std;

class Domain
{
public:
  vector<SpecieP *> species_p;
  Current *current;

#if defined(SWITCH_PUSHER_BORIS_ADAPTIVE) || defined(SWITCH_PUSHER_BORIS) || defined(SWITCH_PUSHER_BORIS_RELATIVISTIC)
  PusherBoris *pusher;
#elif defined(SWITCH_PUSHER_VAY)
  PusherVay *pusher;
#elif defined(SWITCH_PUSHER_HC)
  PusherHC *pusher;
#endif

#ifdef SWITCH_MAXWELL_SOLVER_YEE
  MaxwellSolverYee *maxwell_solver;
#endif

  // Temperature *temperature;
  // Density *density;
  // DensityCharge *charge;
#ifdef ENABLE_COULOMB_COLLISIONS
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

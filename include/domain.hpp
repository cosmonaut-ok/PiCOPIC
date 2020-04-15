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

#ifndef _DOMAIN_HPP_
#define _DOMAIN_HPP_

#include <vector>

#include <algorithm>

#include "grid.hpp"
#include "grid3d.hpp"
#include "geometry.hpp"
#include "timeSim.hpp"
#include "cfg.hpp"
#include "msg.hpp"

#include "specieP.hpp"
#include "fieldE.hpp"
#include "fieldH.hpp"

#ifdef CCS_VILLASENOR_BUNEMAN
#include "currentVB.hpp"
#elif CCS_ZIGZAG
#include "currentZigZag.hpp"
#endif

#include "densityCharge.hpp"
#include "density.hpp"

#ifdef TEMP_CALC_COUNTING
#include "temperatureCounted.hpp"
#elif TEMP_CALC_WEIGHTING
#include "temperatureWeighted.hpp"
#endif

using namespace std;

class Domain
{
public:
  vector<SpecieP *> species_p;

  FieldE *field_e;
  FieldH *field_h;
  Current *current;
  Temperature *temperature;
  Density *density;
  DensityCharge *charge;

  // Temperature temperature;
  // DensityP density_particles;
  // DensityC density_charge;

  Geometry geometry;

  TimeSim *time_sim; // simulation time object

public:
  Domain();
  Domain(Geometry geom, vector<SpecieP *> species, TimeSim* time);

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
};
#endif // end of _DOMAIN_HPP_

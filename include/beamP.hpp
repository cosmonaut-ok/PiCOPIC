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

#ifndef _BEAM_P_HPP_
#define _BEAM_P_HPP_

#include "math/rand.hpp"
#include "specieP.hpp"

class SpecieP;

class BeamP : public SpecieP
{
public:
  double velocity;
  double duration;
  double density;
  double radius;
  double start_time;
  double finish_time;
  unsigned int bunches_amount;
  double bunches_distance;
  double bunch_length;
  double bunch_macro_amount;
  double macro_per_step_to_inject;
  double domain_radius;

public:
  BeamP (unsigned int id, // ID for every particles specie
         string p_name,
         double p_charge, // in electron charges
         double p_mass, // in electron masses
         unsigned int p_macro_amount,
         double p_start_time,
         double b_radius,
         double b_density, // CI units
         double b_amount,
         double b_length,
         double b_distance,
         double b_velocity, // electronvolts
         Geometry *geom,
         TimeSim *t );

  void inject();
  void reflect();

  // just dummy methods. Not used by BeamP
  void fullyfill_spatial_distribution() {};
  void velocity_distribution () {};
};

#endif // end of _BEAM_P_HPP_

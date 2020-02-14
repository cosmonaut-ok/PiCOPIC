#pragma once

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
  double area_radius;

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

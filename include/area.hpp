// #pragma once

#include <vector>

#include <algorithm>

#include "grid.hpp"
#include "grid3d.hpp"
#include "geometry.hpp"
#include "timeSim.hpp"
#include "cfg.hpp"
#include "msg.hpp"

#include "specieP.hpp"
#include "charge.hpp"
#include "fieldE.hpp"
#include "fieldH.hpp"
#include "current.hpp"
#include "temperature.hpp"
#include "density.hpp"

using namespace std;

class Area
{
public:
  vector<SpecieP *> species_p;
  Charge *charge;

  FieldE *field_e;
  FieldH *field_h;
  Current *current;
  Temperature *temperature;
  Density *density;

  // Temperature temperature;
  // DensityP density_particles;
  // DensityC density_charge;

  Geometry geometry;

  TimeSim *time_sim; // simulation time object

public:
  Area();
  Area(Geometry geom, vector<SpecieP *> species, TimeSim* time);

  // wrapper methods
  void distribute();
  void weight_charge();
  void push_particles();
  void weight_current();
  void update_particles_coords_at_half();
  void weight_current_azimuthal();
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

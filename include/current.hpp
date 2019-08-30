#pragma once

#include "geometry.hpp"
#include "timeSim.hpp"
#include "specieP.hpp"

#define SOME_SHIT_DENSITY_STRICT(q, r, dr, dz, delta_t) \
  (q) / (PI * 4. * (r) * (dz) * (dz) * (dr) * (delta_t));

#define SOME_SHIT_DENSITY_R(q, r, dr, dz, delta_t) \
  (q) / (2 * PI * (r) * (dr) * (dz)                \
         * (delta_t) * 2 * (dr));

#define SOME_SHIT_DENSITY_Z(q, r, dr, dz, delta_t) \
  (q) / (2 * PI * (r) * (dr) * (dz)                \
         * (delta_t) * (dz));

class SpecieP;
class Geometry;

class Current
{
public:
  TimeSim *time;
  Geometry *geometry;
  Grid3D<double> current;
  vector<SpecieP *> species_p;

  Current() {};
  Current(Geometry *geom, TimeSim *t, vector<SpecieP *> species);
  ~Current(void) {MSG("Current");};

  void current_distribution();
  void azimuthal_current_distribution();

private:
  void simple_current_distribution (double radius_new,
                                    double longitude_new,
                                    double radius_old,
                                    double longitude_old,
                                    int i_n,
                                    int k_n,
                                    double p_charge);

  void strict_motion_weighting (double radius_new,
                                double longitude_new,
                                double radius_old,
                                double longitude_old,
                                double p_charge);
};

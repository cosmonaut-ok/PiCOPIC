#pragma once

#include "geometry.hpp"
#include "grid.hpp"
#include "specieP.hpp"

class DensityCharge
{
public:
  Geometry *geometry;
  Grid<double> density;
  vector<SpecieP *> species_p;
  
  DensityCharge(void) {};
  DensityCharge(Geometry *geom, vector<SpecieP *> species);
  ~DensityCharge() {};
  
  void calc_density_cylindrical(string specie);
};

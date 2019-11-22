#pragma once

#include "geometry.hpp"
#include "grid.hpp"
#include "specieP.hpp"

class Density
{
public:
  Geometry *geometry;
  Grid<double> density;
  vector<SpecieP *> species_p;
  
  Density(void) {};
  Density(Geometry *geom, vector<SpecieP *> species);
  ~Density() {};
  
  void calc_density_cylindrical(string specie);
};

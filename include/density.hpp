#pragma once

#include "msg.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include "specieP.hpp"

class Density
{
public:
  Density() {};
  Density(Geometry *geom, vector<SpecieP *> species):geometry(geom)
  {
    species_p = species;
  };
  ~Density() {};

// private:
  Grid<double> density;
  Geometry *geometry;
  vector<SpecieP *> species_p;
};

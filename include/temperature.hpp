#pragma once

#include "msg.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include "specieP.hpp"

class Temperature
{
public:
  Temperature() {};
  Temperature(Geometry *geom, vector<SpecieP *> species):geometry(geom)
  {
    species_p = species;
  };
  ~Temperature() {};

// private:
  Grid<double> temperature;
  Geometry *geometry;
  vector<SpecieP *> species_p;
};

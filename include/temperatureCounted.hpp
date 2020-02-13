#pragma once

#include "constant.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include "specieP.hpp"
#include "temperature.hpp"

class TemperatureCounted : public Temperature
{
public:
  Grid<double> count;

  TemperatureCounted(void) {};
  TemperatureCounted(Geometry *geom, vector<SpecieP *> species) : Temperature(geom, species)
  {
    count = Grid<double> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
  };
  ~TemperatureCounted() {};

  void calc_temperature_cylindrical(string specie);

private:
  void weight_temperature_cylindrical(string specie);
};

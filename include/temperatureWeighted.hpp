#pragma once

#include "constant.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include "specieP.hpp"
#include "density.hpp"
#include "temperature.hpp"

class TemperatureWeighted : public Temperature
{
public:
  Density density;

  TemperatureWeighted(void) {};
  TemperatureWeighted(Geometry *geom, vector<SpecieP *> species) : Temperature(geom, species)
  {
    density = Density (geometry, species);
  };
  ~TemperatureWeighted() {};

  void calc_temperature_cylindrical(string specie);

private:
  void weight_temperature_cylindrical(string specie);
};

#pragma once

#include "constant.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include "specieP.hpp"
#include "density.hpp"

class Temperature
{
public:
  Geometry *geometry;
  Grid<double> vel_full;
  Grid<double> vel_r;
  Grid<double> vel_phi;
  Grid<double> vel_z;
  Grid<double> temperature;
  vector<SpecieP *> species_p;
  Density density;

  Temperature(void) {};
  Temperature(Geometry *geom, vector<SpecieP *> species);
  ~Temperature() {};

  void calc_temperature_cylindrical(string specie);

private:
  void weight_temperature_cylindrical(string specie);
};

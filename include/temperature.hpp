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

  Temperature(void) {};
  Temperature(Geometry *geom, vector<SpecieP *> species) : geometry(geom)
  {
    vel_full = Grid<double> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
    vel_r = Grid<double> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
    vel_phi = Grid<double> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
    vel_z = Grid<double> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
    temperature = Grid<double> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
    species_p = species;
  };
  
  ~Temperature(void) {};

  virtual void calc_temperature_cylindrical(string specie) = 0;
};

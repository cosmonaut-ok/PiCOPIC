#pragma once

#include "msg.hpp"
#include "geometry.hpp"
#include "grid.hpp"
#include "specieP.hpp"

class Charge
{
public:
  Charge() {};
  Charge(Geometry *geom, vector<SpecieP *> species);
  ~Charge() {};

  void weight_cylindrical();

  // double **get_rho() const;
  // void inc_rho(int i, int k, double value);

  // void reset_rho();

// private:
  Grid<double> density;
  Geometry *geometry;
  vector<SpecieP *> species_p;
};

#pragma once

#include <fstream>
#include "math/vector3d.hpp"

#include "geometry.hpp"
#include "field.hpp"
#include "specieP.hpp"
#include "current.hpp"
#include "fieldH.hpp"

class Current;
class SpecieP; // define dummy class to avoid an error while loading real links

class FieldE : public Field
{
public:
  vector<SpecieP *> species_p;
  Grid<double> epsilon;
  Grid<double> sigma;

  Current *current;
  FieldH *field_h;

  FieldE(Geometry *geom, TimeSim *t, vector<SpecieP *> species);
  FieldE() {};
  ~FieldE(void) {};

  void set_pml();
  void calc_field_cylindrical();

  // void calc_field(HField *h_field1, Time *time1, Current *current1);
  // void poisson_equation2(Geometry *geom1, ChargeDensity *ro1);
  // void cosine_ftansfrom(double **fi_ro, int lenght_n, int k);
  // void set_homogeneous_efield(double E_r, double E_phi, double E_z);
  // void set_fi_on_z();
  // void boundary_conditions();
  // double* get_field(double radius, double longitude);
  // bool test_poisson_equation(ChargeDensity *rho);
  // void tridiagonal_solve(const double *a,
  //                        const double *b,
  //                        double *c,
  //                        double *d,
  //                        double *x,
  //                        int n);
};

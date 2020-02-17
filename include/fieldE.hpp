/* 
 * This file is part of the PiCOPIC distribution (https://github.com/cosmonaut-ok/PiCOPIC).
 * Copyright (c) 2020 Alexander Vynnyk.
 * 
 * This program is free software: you can redistribute it and/or modify  
 * it under the terms of the GNU General Public License as published by  
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

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
  
  vector3d<double> get_field(double radius, double longitude);

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

/* 
 * This file is part of the PiCoPiC distribution (https://github.com/cosmonaut-ok/PiCoPiC).
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

#ifndef _FIELD_H_HPP_
#define _FIELD_H_HPP_

#include "field.hpp"
#include "fieldE.hpp"
#include "specieP.hpp"

class SpecieP; // define dummy class to avoid an error while loading real links
class FieldE;

class FieldH : public Field
{
public:
  Grid3D<double> field_at_et;
  vector<SpecieP *> species_p;
  FieldE* field_e;

public:
  FieldH(Geometry *geom, TimeSim *t, vector<SpecieP *> species);
  FieldH() {};
  ~FieldH(void) {};

  void calc_field_cylindrical();
  vector3d<double> get_field(double radius, double longitude);

  // HField(Geometry *geom1);
  // HField(void);
  // ~HField(void);
  // void calc_field(EField *e_field1, Time *time1);
  // void set_homogeneous_h(double E_r, double E_phi, double E_z);
  // double* get_field(double radius, double longitude);

  // double *get_1d_field_r();
  // double *get_1d_field_phi();
  // double *get_1d_field_z();
};

#endif // end of _FIELD_H_HPP_

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

#ifndef _FIELD_E_HPP_
#define _FIELD_E_HPP_

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
};
#endif // end of _FIELD_E_HPP_

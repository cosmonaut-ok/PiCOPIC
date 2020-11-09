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

#ifndef _CURRENT_ZIGZAG_
#define _CURRENT_ZIGZAG_

#include <vector>
#include <algorithm>

#include "geometry.hpp"
#include "algo/grid3d.hpp"
#include "timeSim.hpp"
#include "specieP.hpp"
#include "current.hpp"

using namespace std;

class SpecieP;
class Geometry;

#define RELAY_POINT(i1, i2, x1, x2, dx)                           \
  min( min((i1) * (dx), (i2) * (dx)) + (dx),                      \
       max( max((i1) * (dx), (i2) * (dx)), ((x1) + (x2)) / 2.));

class CurrentZigZag : public Current
{
public:
  CurrentZigZag() {};
  CurrentZigZag(Geometry *geom, TimeSim *t, vector<SpecieP *> species) : Current(geom, t, species) {};
  ~CurrentZigZag() {};

  void current_distribution();

private:

};
#endif // end of _CURRENT_ZIGZAG_

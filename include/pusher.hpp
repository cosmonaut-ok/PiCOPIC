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

#ifndef _PUSHER_HPP_
#define _PUSHER_HPP_

#include <vector>

#include "geometry.hpp"
#include "algo/grid3d.hpp"
#include "timeSim.hpp"
#include "specieP.hpp"

#include "maxwellSolver.hpp"

#ifdef MAXWELL_SOLVER_YEE
#include "maxwellSolver/maxwellSolverYee.hpp"
#endif

// class MaxwellSolver;
// class SpecieP;

class Pusher
{
protected:
  MaxwellSolver *maxwell_solver;
  vector<vector<double> * > particles;
  TimeSim *time;
  vector<SpecieP *> species_p;

public:
  Pusher() {};
  Pusher(MaxwellSolver *_maxwell_solver, vector<SpecieP *> _species_p, TimeSim *_time)
  {
    maxwell_solver = _maxwell_solver;
    species_p = _species_p;
    time = _time;
  };

  virtual void operator()() = 0;
};

#endif // end of _PUSHER_HPP_

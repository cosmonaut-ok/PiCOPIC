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

#include "defines.hpp"
#include "msg.hpp"
#include "specieP.hpp"

#ifdef SWITCH_MAXWELL_SOLVER_YEE
#include "maxwellSolver/maxwellSolverYee.hpp"
#endif

#ifdef ENABLE_EXTERNAL_FIELDS
#include "maxwellSolver/externalFields.hpp"
#endif

class Pusher
{
protected:
  MaxwellSolver *maxwell_solver;
#ifdef ENABLE_EXTERNAL_FIELDS
  ExternalFields *external_fields;
#endif
  vector<vector<double> * > particles;
  TimeSim *time;
  vector<SpecieP *> species_p;

public:
  Pusher() {};

#ifdef ENABLE_EXTERNAL_FIELDS
  Pusher(MaxwellSolver *_maxwell_solver, ExternalFields *_external_fields, vector<SpecieP *> _species_p, TimeSim *_time)
    : maxwell_solver(_maxwell_solver), external_fields(_external_fields), time(_time)
#else
  Pusher(MaxwellSolver *_maxwell_solver, vector<SpecieP *> _species_p, TimeSim *_time)
    : maxwell_solver(_maxwell_solver), time(_time)
#endif
  {
    species_p = _species_p;
  };

  virtual void operator()() = 0;
};

#endif // end of _PUSHER_HPP_

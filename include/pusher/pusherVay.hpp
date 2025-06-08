/*
 * This file is part of the PiCoPiC distribution
 * (https://github.com/cosmonaut-ok/PiCoPiC).
 * Copyright (c) 2020,2025 Alexander Vynnyk.
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

#ifndef _PUSHERVAY_HPP_
#define _PUSHERVAY_HPP_

#include "pusher.hpp"

class PusherVay : public Pusher
{
public:
#ifdef ENABLE_EXTERNAL_FIELDS
  PusherVay (MaxwellSolver* _maxwell_solver, ExternalFields* _external_fields, vector<SpecieP*> _species_p,
             TimeSim* _time)
    : Pusher(_maxwell_solver, _external_fields, _species_p, _time) {};
#else
  PusherVay (MaxwellSolver* _maxwell_solver, vector<SpecieP*> _species_p, TimeSim* _time)
    : Pusher(_maxwell_solver, _species_p, _time) {};
#endif

  void operator() ();
};

#endif // end of _PUSHERVAY_HPP_

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

#include "domain.hpp"

Domain::Domain(Geometry _geometry, vector<SpecieP *> _species_p, TimeSim* _time) : geometry(_geometry)
{
  //! Pass geometry, particle species configuration and time
  // TODO: is it ok to pass only configuration,
  // but not prepared particles specie w/o linking?

  time = _time;

  species_p = _species_p;

#ifdef SWITCH_CCS_VB
  current = new CurrentVB(&geometry, time, species_p);
#elif defined(SWITCH_CCS_ZIGZAG)
  current = new CurrentZigZag(&geometry, time, species_p);
#endif

#ifdef SWITCH_MAXWELL_SOLVER_YEE
  maxwell_solver = new MaxwellSolverYee(&geometry, time, current);
#endif

#ifdef ENABLE_EXTERNAL_FIELDS
  external_fields = new ExternalFields(&geometry, time, current);
#endif

#ifdef ENABLE_COULOMB_COLLISIONS
#ifdef SWITCH_COULOMB_COLLISIONS_TA77S
  collisions = new CollisionsTA77S(&geometry, time, species_p);
#elif defined(SWITCH_COULOMB_COLLISIONS_SK98)
  collisions = new CollisionsSK98(&geometry, time, species_p);
#elif defined(SWITCH_COULOMB_COLLISIONS_P12)
  collisions = new CollisionsP12(&geometry, time, species_p);
#endif // WITH_COULOMB_COLLISIONS
#endif // ENABLE_COULOMB_COLLISIONS

  // ! Linking classes:
#if defined(SWITCH_PUSHER_BORIS_ADAPTIVE) || defined(SWITCH_PUSHER_BORIS) || defined(SWITCH_PUSHER_BORIS_RELATIVISTIC)
#ifdef ENABLE_EXTERNAL_FIELDS
  pusher = new PusherBoris(maxwell_solver, external_fields, species_p, time);
#else
  pusher = new PusherBoris(maxwell_solver, species_p, time);
#endif
#elif defined(SWITCH_PUSHER_VAY)
#ifdef ENABLE_EXTERNAL_FIELDS
  pusher = new PusherVay(maxwell_solver, external_fields, species_p, time);
#else
  pusher = new PusherVay(maxwell_solver, species_p, time);
#endif
#elif defined(SWITCH_PUSHER_HC)
#ifdef ENABLE_EXTERNAL_FIELDS
  pusher = new PusherHC(maxwell_solver, external_fields, species_p, time);
#else
  pusher = new PusherHC(maxwell_solver, species_p, time);
#endif
#endif
}

void Domain::distribute()
{
  for (auto i = species_p.begin(); i != species_p.end(); i++)
  {
    (**i).fullyfill_spatial_distribution();
    (**i).bind_cell_numbers(); // calculate cell numbers at initial state
    (**i).velocity_distribution();
  }
}

void Domain::weight_density(string specie)
{
  SpecieP *speciep;

  for (auto ps = species_p.begin(); ps != species_p.end(); ++ps)
    if (specie.compare((**ps).name) == 0)
      speciep = (*ps);

  speciep->calc_density();
}

void Domain::weight_temperature(string specie)
{
  SpecieP *speciep;

  for (auto ps = species_p.begin(); ps != species_p.end(); ++ps)
    if (specie.compare((**ps).name) == 0)
      speciep = (*ps);

  speciep->calc_temperature();
}

void Domain::weight_field_h()
{
  maxwell_solver->calc_field_h();
#ifdef ENABLE_EXTERNAL_FIELDS
  external_fields->calc_field_h();
#endif
}

void Domain::weight_field_e()
{
  maxwell_solver->calc_field_e();
#ifdef ENABLE_EXTERNAL_FIELDS
  external_fields->calc_field_e();
#endif
}

void Domain::reset_current()
{
  current->current = 0;
  current->current.overlay_set(0);
}

void Domain::push_particles ()
{
  (*pusher)();
}

void Domain::update_particles_coords()
{
  // ! update particles coordinates
  for (auto i = species_p.begin(); i != species_p.end(); i++)
    (**i).mover_cylindrical();
}

void Domain::bind_cell_numbers()
{
  // ! bind cell r-number and z-number to each particle
  for (auto i = species_p.begin(); i != species_p.end(); i++)
    (**i).bind_cell_numbers();
}

void Domain::reflect()
{
  // ! update particles coordinates
  for (auto i = species_p.begin(); i != species_p.end(); i++)
    (*i)->reflect();
}


void Domain::particles_back_position_to_rz()
{
  for (auto i = species_p.begin(); i != species_p.end(); i++)
    (*i)->back_position_to_rz();
}

void Domain::particles_back_velocity_to_rz()
{
  for (auto i = species_p.begin(); i != species_p.end(); i++)
    (*i)->back_velocity_to_rz();
}

void Domain::manage_beam()
{
  // injecting bunch
  if (geometry.cell_dims[1] == 0) // inject only from left wall
    for (auto i = species_p.begin(); i != species_p.end(); i++)
      (**i).inject();
}

void Domain::weight_current()
{
  current->current_distribution();
}

void Domain::dump_particle_positions_to_old()
{
  for (auto i = species_p.begin(); i != species_p.end(); i++)
    (**i).dump_position_to_old();
}

#ifdef ENABLE_COULOMB_COLLISIONS
void Domain::collide()
{
  collisions->run();
}
#endif // ENABLE_COULOMB_COLLISIONS

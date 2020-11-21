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
  maxwell_solver = new MaxwellSolverYee(&geometry, time, species_p, current);
#endif

#ifdef SWITCH_TEMP_CALC_COUNTING
  temperature = new TemperatureCounted(&geometry, species_p);
#elif defined(SWITCH_TEMP_CALC_WEIGHTING)
  temperature = new TemperatureWeighted(&geometry, species_p);
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

#ifdef SWITCH_DENSITY_CALC_COUNTING
  density = new DensityCounted(&geometry, species_p);
#elif defined(SWITCH_DENSITY_CALC_WEIGHTING)
  density = new DensityWeighted(&geometry, species_p);
#endif

  charge = new DensityCharge(&geometry, species_p);

  // ! Linking classes:
  // ! * insert field pointers to each particles specie
  for (auto sp = species_p.begin(); sp != species_p.end(); ++sp)
  {
    (**sp).maxwell_solver = maxwell_solver;
  }

#if defined(SWITCH_PUSHER_BORIS_ADAPTIVE) || defined(SWITCH_PUSHER_BORIS) || defined(SWITCH_PUSHER_BORIS_RELATIVISTIC)
  pusher = new PusherBoris(maxwell_solver, species_p, time);
#elif defined(SWITCH_PUSHER_VAY)
  pusher = new PusherVay(maxwell_solver, species_p, time);
#elif defined(SWITCH_PUSHER_HC)
  pusher = new PusherHC(maxwell_solver, species_p, time);
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
  density->density = 0;
  density->density.overlay_set(0);
  (*density)(specie);
}

void Domain::weight_temperature(string specie)
{
  temperature->tmpr = 0;
  temperature->tmpr.overlay_set(0);
  (*temperature)(specie);
}

void Domain::weight_charge(string specie)
{
  charge->density = 0;
  charge->density.overlay_set(0);
  (*charge)(specie);
}

void Domain::weight_field_h()
{
  maxwell_solver->calc_field_h();
}

void Domain::weight_field_e()
{
  maxwell_solver->calc_field_e();
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
  if (geometry.left_z_grid_number == 0) // inject only from left wall
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

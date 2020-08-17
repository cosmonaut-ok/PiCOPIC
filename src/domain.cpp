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

#include "domain.hpp"

Domain::Domain(){}

Domain::Domain(Geometry geom, vector<SpecieP *> species, TimeSim* time):geometry(geom)
{
  //! Pass geometry, particle species configuration and time
  // TODO: is it ok to pass only configuration,
  // but not prepared particles specie w/o linking?

  time_sim = time;

  species_p = species;

  field_e = new FieldE(&geometry, time, species_p);
  field_h = new FieldH(&geometry, time, species_p);

#ifdef CCS_VILLASENOR_BUNEMAN
  current = new CurrentVB(&geometry, time, species_p);
#elif CCS_ZIGZAG
  current = new CurrentZigZag(&geometry, time, species_p);
#endif

#ifdef TEMP_CALC_COUNTING
  temperature = new TemperatureCounted(&geometry, species_p);
#elif TEMP_CALC_WEIGHTING
  temperature = new TemperatureWeighted(&geometry, species_p);
#endif

#ifdef COLLISIONS
#ifdef COULOMB_COLLISIONS_TANAKA
  collisions = new CollisionsTanaka(&geometry, time, species_p);
#elif COULOMB_COLLISIONS_SENTOKU_M
  collisions = new CollisionsSentokuM(&geometry, time, species_p);
#endif
#endif

#ifdef DENSITY_CALC_COUNTING
  density = new DensityCounted(&geometry, species_p);
#elif DENSITY_CALC_WEIGHTING
  density = new DensityWeighted(&geometry, species_p);
#endif

  charge = new DensityCharge(&geometry, species_p);

  // ! Linking classes:
  // ! * insert pointer to field_e into field_h
  field_h->field_e = field_e;
  // ! * insert pointer to current and field_h into field_e
  field_e->current = current;
  field_e->field_h = field_h;
  // ! * insert field pointers to each particles specie
  for (auto sp = species_p.begin(); sp != species_p.end(); ++sp)
  {
    (**sp).field_h = field_h;
    (**sp).field_e = field_e;
  }
}

void Domain::distribute()
{
  for (auto i = species_p.begin(); i != species_p.end(); i++)
  {
    (**i).fullyfill_spatial_distribution();
    (**i).velocity_distribution();
  }
}

void Domain::weight_density(string specie)
{
  density->density = 0;
  density->density.overlay_set(0);
  density->calc_density_cylindrical(specie);
}

void Domain::weight_temperature(string specie)
{
  temperature->temperature = 0;
  temperature->temperature.overlay_set(0);
  temperature->calc_temperature_cylindrical(specie);
}

void Domain::weight_charge(string specie)
{
  charge->density = 0;
  charge->density.overlay_set(0);
  charge->calc_density_cylindrical(specie);
}

void Domain::weight_field_h()
{
  field_h->calc_field_cylindrical();
}

void Domain::weight_field_e()
{
  field_e->calc_field_cylindrical();
}

void Domain::reset_current()
{
  current->current = 0;
}

// void Domain::reset_charge()
// {
//   charge->density = 0;
// }

void Domain::push_particles ()
{
  // ! update particles velocity
  for (auto i = species_p.begin(); i != species_p.end(); i++)
  {
#if defined PUSHER_BORIS_ADAPTIVE || defined PUSHER_BORIS_CLASSIC || defined PUSHER_BORIS_RELATIVISTIC
    (*i)->boris_pusher();
#elif defined PUSHER_VAY
    (*i)->vay_pusher();
#elif defined PUSHER_HIGUERA_CARY
    (*i)->hc_pusher();
#else
    LOG_S(FATAL) << "Undefined particles pusher used";
#endif
  }
}

void Domain::update_particles_coords()
{
  // ! update particles coordinates
  for (auto i = species_p.begin(); i != species_p.end(); i++)
    (*i)->mover_cylindrical();
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

#ifdef COLLISIONS
void Domain::collide()
{
  collisions->run();
}
#endif

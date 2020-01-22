#include "area.hpp"

Area::Area(){}

Area::Area(Geometry geom, vector<SpecieP *> species, TimeSim* time):geometry(geom)
{
  //! Pass geometry, particle species configuration and time
  // TODO: is it ok to pass only configuration,
  // but not prepared particles specie w/o linking?

  time_sim = time;

  species_p = species;

  field_e = new FieldE(&geometry, time, species_p);
  field_h = new FieldH(&geometry, time, species_p);

  current = new Current(&geometry, time, species_p);

  temperature = new Temperature(&geometry, species_p);
  density = new Density(&geometry, species_p);
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

void Area::distribute()
{
  for (auto i = species_p.begin(); i != species_p.end(); i++)
  {
    (**i).wakeup();
    (**i).fullyfill_spatial_distribution();
    (**i).thermal_velocity_distribution();
  }
}

void Area::weight_density(string specie)
{
  density->density = 0;
  density->density.overlay_set(0);
  density->calc_density_cylindrical(specie);
}

void Area::weight_temperature(string specie)
{
  temperature->temperature = 0;
  temperature->temperature.overlay_set(0);
  temperature->calc_temperature_cylindrical(specie);
}

void Area::weight_charge(string specie)
{
  charge->density = 0;
  charge->density.overlay_set(0);
  charge->calc_density_cylindrical(specie);
}

void Area::weight_field_h()
{
  field_h->calc_field_cylindrical();
}

void Area::weight_field_e()
{
  field_e->calc_field_cylindrical();
}

void Area::reset_current()
{
  current->current = 0;
}

// void Area::reset_charge()
// {
//   charge->density = 0;
// }

void Area::push_particles ()
{
  // ! update particles velocity
  for (auto i = species_p.begin(); i != species_p.end(); i++)
    (*i)->boris_pusher();
}

void Area::update_particles_coords_at_half()
{
  // ! update particles coordinates
  for (auto i = species_p.begin(); i != species_p.end(); i++)
    (*i)->half_step_mover_cylindrical();
}

void Area::reflect()
{
  // ! update particles coordinates
  for (auto i = species_p.begin(); i != species_p.end(); i++)
    (*i)->reflect();
}


void Area::particles_back_position_to_rz()
{
  for (auto i = species_p.begin(); i != species_p.end(); i++)
    (*i)->back_position_to_rz();
}

void Area::particles_back_velocity_to_rz()
{
  for (auto i = species_p.begin(); i != species_p.end(); i++)
    (*i)->back_velocity_to_rz();
}

void Area::manage_beam()
{
  // injecting bunch
  if (geometry.left_z_grid_number == 0) // inject only from left wall
    for (auto i = species_p.begin(); i != species_p.end(); i++)
      (**i).inject();
}

void Area::weight_current_azimuthal()
{
  current->azimuthal_current_distribution();
}

void Area::weight_current()
{
  current->current_distribution();
}

void Area::dump_particle_positions_to_old()
{
  for (auto i = species_p.begin(); i != species_p.end(); i++)
    (**i).dump_position_to_old();
}

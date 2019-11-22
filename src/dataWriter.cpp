#include "dataWriter.hpp"

DataWriter::DataWriter(string a_path, string a_component,
                       string a_specie, unsigned int a_shape,
                       int *a_size, unsigned int a_schedule,
                       bool a_compress,  unsigned int a_compress_level,
                       Geometry *a_geom, TimeSim *a_time, Grid<Area *> a_areas)
{
  // ! path - path to save the results
  // ! component - some parameter, that should be dumped:
  // !             charge, current(r, phi, z), field (E, H) (r, phi, z),
  // !             temperature, density, particles position or velocity
  // ! specie - particles specie to which writer should be applied
  // !          (only for specie-dependent components, like density)
  // ! shape - shape of writer: dot (3), column (col - 1), row (2), rectangle (rec - 0)
  // ! size - size of writer (could be used only part of size numbers - depends on shape):
  // !        dot - [r, z, 0, 0]
  // !        col - [r_begin, z, r_end, 0]
  // !        row - [0, z_begin, 0, z_end]
  // !        rec - [r_begin, z_begin, r_end, z_end]
  // ! compress - if output data should be compressed
  // ! compress_level - output data compress level
  // ! geom - pointer to global geometry object
  // ! time - pointer to global time object
  // ! areas - grid of pointers to areas

  path = a_path;
  component = a_component;
  specie = a_specie;
  if (a_shape < 4)
    shape = a_shape;
  else
    LOG_CRIT("Unknown DataWriter shape: " << a_shape, 1);

  for (unsigned int i = 0; i < 4; i++)
    size[i] = a_size[i];

  schedule = a_schedule;
  compress = a_compress;
  compress_level = a_compress_level;
  geometry = a_geom;
  time = a_time;
  areas = a_areas;

  out_data = Grid<double> (geometry->r_grid_amount, geometry->z_grid_amount, 0);
  out_data_plain = vector<double> (0);
  // string shape_name, size_name;

  name = component + "/";

  if (component.compare("temperature") == 0
      || component.compare("density") == 0
      || component.compare("position") == 0  // TODO: implement specie-dependent
      || component.compare("velocity") == 0
      || component.compare("p_charge") == 0  // TODO: implement specie-dependent
      || component.compare("p_mass") == 0) //        behavior for positions and velocities
  {
    name += specie;
    name += "/";
  }

  switch (shape)
  {
  case 0:
    name += "rec";
    name += "/";
    name += to_string(size[0]);
    name += "-";
    name += to_string(size[2]);
    name += "_";
    name += to_string(size[1]);
    name += "-";
    name += to_string(size[3]);
    break;
  case 1:
    name += "col";
    name += "/";
    name += to_string(size[1]);
    break;
  case 2:
    name += "row";
    name += "/";
    name += to_string(size[0]);
    break;
  case 3:
    name += "dot";
    name += "/";
    name += to_string(size[0]);
    name += "_";
    name += to_string(size[1]);
    break;
  }

#ifdef USE_HDF5
  engine = OutEnginePlain (path + "/" + name, shape, size, true, compress, compress_level);
#else
  engine = OutEnginePlain (path + "/" + name, shape, size, true, compress, compress_level);
#endif // USE_HDF5

}

void DataWriter::go()
{
  int current_time_step = (int)(time->current / time->step);
  int is_run = current_time_step % schedule;

  if (is_run == 0)
  {
    if (! DEBUG)
    {
      int print_header_step = 30;
      if (time->print_header_counter % print_header_step == 0)
      {
        MSG("+------------+-------------+-------------------------------------------------+-------+------------------+---------------------+");
        MSG(left << setw(13) << "|    Step"
            << left << setw(14) << "| Saved Frame"
            << left << setw(50) << "|                Dumping Probe Name"
            << left << setw(8) << "| Shape"
            << left << setw(19) << "| Model Time (sec)"
            << left << setw(22) << "| Simulation Duration"
            << left << "|"
          );
        MSG("+------------+-------------+-------------------------------------------------+-------+------------------+---------------------+");
      }
    }

    string dump_step = to_string((int)current_time_step / schedule);
    string shape_name;
    string component_name;

    if (component.compare("position") == 0 || component.compare("velocity") == 0)
    {
      LOG_DBG("Launch 1d-vector-shaped writer " << component << " at step " << dump_step);
      for (unsigned int i=0; i < 3; ++i)
      {
        string dump_step_name = dump_step + "_" + to_string(i);
        out_data_plain.resize(0);
        merge_particle_areas(component, i, specie);
        engine.write_1d_vector(dump_step_name, out_data_plain);
        shape_name = "rec";
        component_name = specie
          + "-" + component + ":" + to_string(size[0]) + "-"
          + to_string(size[2]) + "x"  + to_string(size[1]) + "-"
          + to_string(size[3]);
      };
    }
    else if (component.compare("p_mass") == 0
             || component.compare("p_charge") == 0)
    {
      LOG_DBG("Launch 1d-vector-shaped writer " << component << " at step " << dump_step);
      string dump_step_name = dump_step;
      out_data_plain.resize(0);
      merge_particle_areas(component, 0, specie);
      engine.write_1d_vector(dump_step_name, out_data_plain);
      shape_name = "rec";
      component_name = specie
        + "-" + component + ":" + to_string(size[0]) + "-"
        + to_string(size[2]) + "x"  + to_string(size[1]) + "-"
        + to_string(size[3]);
    }
    else
    {
      out_data = 0;
      merge_areas(component, specie);

      // shape: dot (3), column (col - 1), row (2), rectangle (rec - 0)
      switch(shape)
      {
      case 0:
        LOG_DBG("Launch rectangle-shaped writer " << component << " at step " << dump_step);
        engine.write_rec(dump_step, out_data);
        shape_name = "rec";
        component_name = component + ":" + to_string(size[0]) + "-"  + to_string(size[2]) + "x"  + to_string(size[1]) + "-"  + to_string(size[3]);
        break;
      case 1:
        LOG_DBG("Launch column-shaped writer " << component << " at step " << dump_step);
        engine.write_vec(dump_step, out_data);
        shape_name = "col";
        component_name = component + ":" + to_string(size[1]);
        break;
      case 2:
        LOG_DBG("Launch row-shaped writer " << component << " at step " << dump_step);
        engine.write_vec(dump_step, out_data);
        shape_name = "row";
        component_name = component + ":" + to_string(size[0]);
        break;
      case 3:
        LOG_DBG("Launch dot-shaped writer " << component << " at step " << dump_step);
        engine.write_dot(dump_step, out_data);
        shape_name = "dot";
        component_name = component + ":" + to_string(size[0]) + "x"  + to_string(size[1]);
        break;
      }
    }

    if (! DEBUG)
    {
      std::ostringstream time_conv;
      time_conv << time->current;
      std::string cur_time_str = time_conv.str();

      MSG(left << setw(13) << "| " + to_string(current_time_step)
          << left << setw(14) << "| " + dump_step
          << left << setw(50) << "| " + component_name
          << left << setw(8) << "| " + shape_name
          << left << setw(19) << "| " + cur_time_str
          << left << setw(22) << "| " + (string)lib::get_simulation_duration()
          << left << "|"
        );
    }

    ++time->print_header_counter;
  }
}

void DataWriter::merge_areas(string component, string specie)
{
  // prepare to merge areas
  for (unsigned int r = 0; r < geometry->areas_by_r; ++r)
    for (unsigned int z = 0; z < geometry->areas_by_z; ++z)
    {
      Area *sim_area = areas(r, z);

      if (component.compare("density") == 0)
        sim_area->weight_density(specie);
      else if (component.compare("temperature") == 0)
        sim_area->weight_temperature(specie);
    }

  for (unsigned int i=0; i < geometry->areas_by_r; i++)
    for (unsigned int j = 0; j < geometry->areas_by_z; j++)
    {
      Area *sim_area = areas(i, j);

      // update grid
      if (i < geometry->areas_by_r - 1)
      {
        Area *dst_area = areas(i+1, j);
        if (component.compare("density") == 0)
          sim_area->density->density.overlay_top(dst_area->density->density);
        else if (component.compare("temperature") == 0)
          sim_area->temperature->temperature.overlay_top(dst_area->temperature->temperature);
      }

      if (j < geometry->areas_by_z - 1)
      {
        Area *dst_area = areas(i, j + 1);
        if (component.compare("density") == 0)
          sim_area->density->density.overlay_right(dst_area->density->density);
        else if (component.compare("temperature") == 0)
          sim_area->temperature->temperature.overlay_right(dst_area->temperature->temperature);
      }

      // if (i < geometry_global->areas_by_r - 1 && j < geometry_global->areas_by_z - 1)
      // {
      //   Area *dst_area = areas(i + 1, j + 1);
      //   sim_area->current->current.overlay_top_right(dst_area->current->current);
      // }
    }

  // merge values
  for (unsigned int r = 0; r < geometry->areas_by_r; ++r)
    for (unsigned int z = 0; z < geometry->areas_by_z; ++z)
    {
      Grid<double> value;
      Area *sim_area = areas(r, z);

      if (component.compare("E_r") == 0)
        value = sim_area->field_e->field[0];
      else if (component.compare("E_phi") == 0)
        value = sim_area->field_e->field[1];
      else if (component.compare("E_z") == 0)
        value = sim_area->field_e->field[2];
      else if (component.compare("H_r") == 0)
        value = sim_area->field_h->field_at_et[0];
      else if (component.compare("H_phi") == 0)
        value = sim_area->field_h->field_at_et[1];
      else if (component.compare("H_z") == 0)
        value = sim_area->field_h->field_at_et[2];
      else if (component.compare("J_r") == 0)
        value = sim_area->current->current[0];
      else if (component.compare("J_phi") == 0)
        value = sim_area->current->current[1];
      else if (component.compare("J_z") == 0)
        value = sim_area->current->current[2];
      else if (component.compare("temperature") == 0)
        value = sim_area->temperature->temperature;
      else if (component.compare("density") == 0)
        value = sim_area->density->density;
      else
        LOG_ERR("Unknown DataWriter component ``" << component << "''");

      for (unsigned int i = 0; i < value.x_size; ++i) // shifting to avoid overlay areas
        for (unsigned int j = 0; j < value.y_size; ++j)
          out_data.set (i + r * sim_area->geometry.r_grid_amount,
                        j + z * sim_area->geometry.z_grid_amount,
                        value(i, j));
    }
}

void DataWriter::merge_particle_areas(string parameter,
                                      unsigned int component,
                                      string specie)
{
  Grid<double> value;
  if (shape == 0)
  {
    for (unsigned int r = 0; r < geometry->areas_by_r; ++r)
      for (unsigned int z = 0; z < geometry->areas_by_z; ++z)
      {
        Area *sim_area = areas(r, z);
        vector<SpecieP *> species = sim_area->species_p;

        for (auto ps = species.begin(); ps != species.end(); ++ps)
        {
          if (specie.compare((**ps).name) == 0)
          {
            vector< vector<double> * > pp = (**ps).particles;

            for (unsigned int p = 0; p < pp.size(); ++p)
            {
              if (parameter.compare("position") == 0)
                out_data_plain.push_back((*pp[p])[component]);
              else if (parameter.compare("velocity") == 0)
                out_data_plain.push_back((*pp[p])[component+6]); // positions(3), positions_old(3), velocity(3)
              else if (parameter.compare("p_charge") == 0)
                out_data_plain.push_back((*pp[p])[component+9]); // positions(3), positions_old(3), velocity(3)
              else if (parameter.compare("p_mass") == 0)
                out_data_plain.push_back((*pp[p])[component+10]); // positions(3), positions_old(3), velocity(3)
              else
                LOG_CRIT("Unknown DataWriter component ``" << parameter << "''", 1);
            }
          }
      }
      }
  }
  else
    LOG_CRIT("Incorrect DataWriter shape ``" << shape
             << "''(not allowed for particle sets) or size", 1);
}

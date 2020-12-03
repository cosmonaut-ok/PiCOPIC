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

#include "outController.hpp"

#ifdef ENABLE_HDF5
OutController::OutController ( HighFive::File* _file,
                               Geometry *_geometry, TimeSim *_time,
                               vector<probe> &_probes, SMB *_smb,
                               std::string _metadata )
  : geometry(_geometry), time(_time), smb(_smb)
#else
OutController::OutController ( Geometry *_geometry, TimeSim *_time,
                               vector<probe> &_probes, SMB *_smb,
                               std::string _metadata )
  : geometry(_geometry), time(_time), smb(_smb)
#endif
{
// write metadata
#ifdef ENABLE_HDF5
  hdf5_file = _file;
  OutEngineHDF5 engine (hdf5_file, _probes[0].path, {0,0}, {0,0}, true, false);
#endif // end of ENABLE_HDF5
  engine.write_metadata( _metadata );

  // initialize probes
  probes = _probes;

  // setup probes for every domain
  for (auto prb = probes.begin(); prb != probes.end(); ++prb)
  {
    vector<short> prb_size = {prb->r_start, prb->z_start, prb->r_end, prb->z_end};

    // initialize engine paths
#ifdef ENABLE_HDF5
    hdf5_file = _file;
    vector<size_t> hdf5_prb_size;

    switch (prb->shape)
    {
    case 0:
      hdf5_prb_size = {prb->r_end - prb->r_start, prb->z_end - prb->z_start};
      break;
    case 1:
      hdf5_prb_size = {geometry->cell_amount[0]};
      break;
    case 2:
      hdf5_prb_size = {geometry->cell_amount[1]};
      break;
    case 3:
      hdf5_prb_size = {1};
    }

    OutEngineHDF5 engine (hdf5_file, prb->path, hdf5_prb_size, {0,0}, true, false);
#endif // end of ENABLE_HDF5

    engine.create_dataset();

    // push probes to domains
    for (unsigned int r = 0; r < geometry->domains_amount[0]; ++r)
      for (unsigned int z = 0; z < geometry->domains_amount[1]; ++z)
      {
        Domain *dmn = smb->domains(r, z);

        vector<int> domain_size = { dmn->geometry.cell_dims[0],
                                    dmn->geometry.cell_dims[1],
                                    dmn->geometry.cell_dims[2],
                                    dmn->geometry.cell_dims[3] };

        // check if probe possition intersects the domain
        if (
          ( prb->shape == 0                      // condition for rec
            && prb_size[0] <= domain_size[2]
            && prb_size[1] <= domain_size[3]
            && prb_size[2] >= domain_size[0]
            && prb_size[3] >= domain_size[1] )
          || ( prb->shape == 1                    // condition for column
               && prb_size[3] >= domain_size[1]
               && prb_size[3] <= domain_size[3]
            )
          || ( prb->shape == 2                    // condition for row
               && prb_size[2] >= domain_size[0]
               && prb_size[2] <= domain_size[2]
            )
          || ( prb->shape == 3                      // condition for dot
               && prb_size[2] >= domain_size[0]
               && prb_size[2] <= domain_size[2]
               && prb_size[3] >= domain_size[1]
               && prb_size[3] <= domain_size[3] )
          )
        {
          // calculate effective probe size
          vector<short> eff_prb_size = {0,0,0,0};
          if (prb_size[0] < domain_size[0])
            eff_prb_size[0] = domain_size[0];
          else
            eff_prb_size[0] = prb_size[0];

          if (prb_size[1] < domain_size[1])
            eff_prb_size[1] = domain_size[1];
          else
            eff_prb_size[1] = prb_size[1];

          if (prb_size[2] > domain_size[2])
            eff_prb_size[2] = domain_size[2];
          else
            eff_prb_size[2] = prb_size[2];

          if (prb_size[3] > domain_size[3])
            eff_prb_size[3] = domain_size[3];
          else
            eff_prb_size[3] = prb_size[3];

          // shift to the begining
          eff_prb_size[0] -= dmn->geometry.cell_dims[0];
          eff_prb_size[2] -= dmn->geometry.cell_dims[0];
          eff_prb_size[1] -= dmn->geometry.cell_dims[1];
          eff_prb_size[3] -= dmn->geometry.cell_dims[1];

          // calculate effective engine offset
          vector<size_t> eff_engine_offset;

          if (prb->shape == 0 || prb->shape == 1)
          {
            if (dmn->geometry.cell_dims[0] <= prb_size[0])
              eff_engine_offset.push_back(0);
            else
              eff_engine_offset.push_back(dmn->geometry.cell_dims[0] - prb_size[0]);
          }
          if (prb->shape == 0 || prb->shape == 2)
          {
            if (dmn->geometry.cell_dims[1] <= prb_size[1])
              eff_engine_offset.push_back(0);
            else
              eff_engine_offset.push_back(dmn->geometry.cell_dims[1] - prb_size[1]);
          }
          Grid<double> *value;

          // map outWriter to the grid of values
          if (prb->component.compare("E_r") == 0)
            value = &(dmn->maxwell_solver->field_e[0]);
          else if (prb->component.compare("E_phi") == 0)
            value = &(dmn->maxwell_solver->field_e[1]);
          else if (prb->component.compare("E_z") == 0)
            value = &(dmn->maxwell_solver->field_e[2]);
          else if (prb->component.compare("H_r") == 0)
            value = &(dmn->maxwell_solver->field_h_at_et[0]);
          else if (prb->component.compare("H_phi") == 0)
            value = &(dmn->maxwell_solver->field_h_at_et[1]);
          else if (prb->component.compare("H_z") == 0)
            value = &(dmn->maxwell_solver->field_h_at_et[2]);
          else if (prb->component.compare("J_r") == 0)
            value = &(dmn->current->current[0]);
          else if (prb->component.compare("J_phi") == 0)
            value = &(dmn->current->current[1]);
          else if (prb->component.compare("J_z") == 0)
            value = &(dmn->current->current[2]);
          else if (prb->component.compare("temperature") == 0)
          {
            SpecieP *speciep;
            for (auto ps = dmn->species_p.begin(); ps != dmn->species_p.end(); ++ps)
              if (prb->specie.compare((**ps).name) == 0)
                speciep = (*ps);

            value = &(speciep->temperature_map);
          }
          else if (prb->component.compare("density") == 0)
          {
            SpecieP *speciep;
            for (auto ps = dmn->species_p.begin(); ps != dmn->species_p.end(); ++ps)
              if (prb->specie.compare((**ps).name) == 0)
                speciep = (*ps);

            value = &(speciep->density_map);
          }
          else
            LOG_S(ERROR) << "Unknown probe component ``" << prb->component << "''";

          // create and push out writer
          OutWriter writer (hdf5_file, prb->path, prb->shape,
                            eff_prb_size, eff_engine_offset,
                            prb->schedule, true, (unsigned short)0,
                            time, value);
          writer.hdf5_file = hdf5_file;

          dmn->out_writers.push_back(writer);
        }
      }
  }
}

void OutController::init_datasets()
{
  // create empty
  for (auto prb = probes.begin(); prb != probes.end(); ++prb)
  {
    vector<size_t> prb_size = {prb->r_start, prb->z_start, prb->r_end, prb->z_end};
    vector<size_t> offset = {0, 0};

#ifdef ENABLE_HDF5
    vector<size_t> hdf5_prb_size;

    switch (prb->shape)
    {
    case 0:
      hdf5_prb_size = {prb->r_end - prb->r_start, prb->z_end - prb->z_start};
      break;
    case 1:
      hdf5_prb_size = {geometry->cell_amount[0]};
      break;
    case 2:
      hdf5_prb_size = {geometry->cell_amount[1]};
      break;
    case 3:
      hdf5_prb_size = {1};
    }

    OutEngineHDF5 engine (hdf5_file, prb->path, hdf5_prb_size, offset, true, false);
#endif // end of ENABLE_HDF5

    int current_time_step = ceil(time->current / time->step);
    int is_run = current_time_step % prb->schedule;

    if (is_run == 0)
    {
      //// print header
#ifndef ENABLE_DEBUG
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
#endif // ENABLE_DEBUG

      string shape_name;
      size_t slices = (size_t)(ceil(current_time_step / prb->schedule));
      vector<unsigned int> probe_size;

      switch (prb->shape)
      {
      case 0:
        shape_name = "rec";
        probe_size.push_back(prb->r_end - prb->r_start);
        probe_size.push_back(prb->z_end - prb->z_start);
        break;
      case 1:
        shape_name = "col";
        probe_size.push_back(geometry->cell_amount[0]);
        break;
      case 2:
        shape_name = "row";
        probe_size.push_back(geometry->cell_amount[1]);
        break;
      case 3:
        shape_name = "dot";
        probe_size.push_back(1);
        break;
      }

      engine.extend_dataset(slices);

      ////
      //// calculate and overlay temperatures and densities before dump
      ///
      // calculate temperature and/or density before dump
      for (unsigned int r = 0; r < geometry->domains_amount[0]; ++r)
        for (unsigned int z = 0; z < geometry->domains_amount[1]; ++z)
        {
          Domain *sim_domain = smb->domains(r, z);
          SpecieP *speciep;

          for (auto ps = sim_domain->species_p.begin(); ps != sim_domain->species_p.end(); ++ps)
            if (prb->specie.compare((**ps).name) == 0)
              speciep = (*ps);

          if (prb->component.compare("density") == 0)
            speciep->calc_density();
          else if (prb->component.compare("temperature") == 0)
            speciep->calc_temperature();
        }

      // overlay domains before dump
      for (unsigned int i=0; i < geometry->domains_amount[0]; i++)
        for (unsigned int j = 0; j < geometry->domains_amount[1]; j++)
        {
          Domain *sim_domain = smb->domains(i, j);
          SpecieP *speciep;

          for (auto ps = sim_domain->species_p.begin(); ps != sim_domain->species_p.end(); ++ps)
            if (prb->specie.compare((**ps).name) == 0)
              speciep = (*ps);

          // update grid
          if (i < geometry->domains_amount[0] - 1)
          {
            Domain *dst_domain = smb->domains(i+1, j);
            SpecieP *speciep_dst;

            for (auto ps = dst_domain->species_p.begin(); ps != dst_domain->species_p.end(); ++ps)
              if (prb->specie.compare((**ps).name) == 0)
                speciep_dst = (*ps);

            if (prb->component.compare("density") == 0)
              speciep->density_map.overlay_x(speciep_dst->density_map);
            else if (prb->component.compare("temperature") == 0)
              speciep->temperature_map.overlay_x(speciep_dst->temperature_map);
          }

          if (j < geometry->domains_amount[1] - 1)
          {
            Domain *dst_domain = smb->domains(i, j + 1);
            SpecieP *speciep_dst;

            for (auto ps = dst_domain->species_p.begin(); ps != dst_domain->species_p.end(); ++ps)
              if (prb->specie.compare((**ps).name) == 0)
                speciep_dst = (*ps);

            if (prb->component.compare("density") == 0)
              speciep->density_map.overlay_y(speciep_dst->density_map);
            else if (prb->component.compare("temperature") == 0)
              speciep->temperature_map.overlay_y(speciep_dst->temperature_map);
          }

          // if (i < geometry_global->domains_amount[0] - 1 && j < geometry_global->domains_amount[1] - 1)
          // {
          //   Domain *dst_domain = smb->domains(i + 1, j + 1);
          //   sim_domain->current->current.overlay_xy(dst_domain->current->current);
          // }
        }

#ifndef DEBUG
      string dump_step = to_string((int)current_time_step / prb->schedule);
      std::ostringstream time_conv;
      time_conv << time->current;
      std::string cur_time_str = time_conv.str();

      MSG(left << setw(13) << "| " + to_string(current_time_step)
          << left << setw(14) << "| " + dump_step
          << left << setw(50) << "| " + prb->path
          << left << setw(8) << "| " + shape_name
          << left << setw(19) << "| " + cur_time_str
          << left << setw(22) << "| " + (string)algo::common::get_simulation_duration()
          << left << "|"
        );
#endif
      ++time->print_header_counter;
    }
  }
}

void OutController::operator()()
{
  init_datasets();

  for (unsigned int r = 0; r < geometry->domains_amount[0]; ++r)
    for (unsigned int z = 0; z < geometry->domains_amount[1]; ++z)
    {
      Domain *dmn = smb->domains(r, z);
      LOG_S(MAX) << "Launching writers for domain ``" << r << "," << z << "''";
      for (auto w = dmn->out_writers.begin(); w != dmn->out_writers.end(); ++w)
        (*w)();
    }
}
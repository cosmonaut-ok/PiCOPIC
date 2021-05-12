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
                               picojson::value _metadata,
                               bool _print_progress_table )
  : geometry(_geometry), time(_time), smb(_smb)
#else
OutController::OutController ( Geometry *_geometry, TimeSim *_time,
                               vector<probe> &_probes, SMB *_smb,
                               picojson::value _metadata,
                               bool _print_progress_table )
  : geometry(_geometry), time(_time), smb(_smb)
#endif
{
// write metadata
#ifdef ENABLE_HDF5
  hdf5_file = _file;
  OutEngineHDF5 engine (hdf5_file, _probes[0].path, {0,0}, {0,0}, true, false);
#endif // end of ENABLE_HDF5

  engine.write_metadata( _metadata );

  print_progress_table = _print_progress_table;

  // initialize probes
  probes = _probes;

  // setup probes for every domain
  for (auto prb = probes.begin(); prb != probes.end(); ++prb)
  {
    vector<size_t> prb_size = {prb->dims[0], prb->dims[1], prb->dims[2], prb->dims[3]};

    // initialize engine paths
#ifdef ENABLE_HDF5
    hdf5_file = _file;
    vector<size_t> hdf5_prb_size;

    switch (prb->shape)
    {
    case 0:
      hdf5_prb_size = {prb->dims[2] - prb->dims[0], prb->dims[3] - prb->dims[1]};
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

        vector<size_t> domain_size = { dmn->geometry.cell_dims[0],
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

          if (prb->shape == 0 || prb->shape == 1 || prb->shape == 3)
            eff_engine_offset.push_back(eff_prb_size[0] + dmn->geometry.cell_dims[0] - prb->dims[0]);

          if (prb->shape == 0 || prb->shape == 2 || prb->shape == 3)
            eff_engine_offset.push_back(eff_prb_size[1] + dmn->geometry.cell_dims[1] - prb->dims[1]);

          Grid<double> *value;

          // map outWriter to the grid of values
          if (prb->component.compare("E/r") == 0)
            value = &(dmn->maxwell_solver->field_e[0]);
          else if (prb->component.compare("E/phi") == 0)
            value = &(dmn->maxwell_solver->field_e[1]);
          else if (prb->component.compare("E/z") == 0)
            value = &(dmn->maxwell_solver->field_e[2]);
          else if (prb->component.compare("H/r") == 0)
            value = &(dmn->maxwell_solver->field_h_at_et[0]);
          else if (prb->component.compare("H/phi") == 0)
            value = &(dmn->maxwell_solver->field_h_at_et[1]);
          else if (prb->component.compare("H/z") == 0)
            value = &(dmn->maxwell_solver->field_h_at_et[2]);
          else if (prb->component.compare("J/r") == 0)
            value = &(dmn->current->current[0]);
          else if (prb->component.compare("J/phi") == 0)
            value = &(dmn->current->current[1]);
          else if (prb->component.compare("J/z") == 0)
            value = &(dmn->current->current[2]);
          else if (prb->component.compare("temperature") == 0)
          {
            SpecieP *speciep;
            bool specie_exists = false;
            for (auto ps = dmn->species_p.begin(); ps != dmn->species_p.end(); ++ps)
              if (prb->specie.compare((**ps).name) == 0)
              {
                speciep = (*ps);
                value = &(speciep->temperature_map);
                specie_exists = true;
              }

            if (! specie_exists)
            {
              LOG_S(WARNING) << "There is no particles specie ``"
                             << prb->specie
                             << "''. Probe dump is impossible";
              continue;
            }
          }
          else if (prb->component.compare("density") == 0)
          {
            SpecieP *speciep;
            bool specie_exists = false;

            for (auto ps = dmn->species_p.begin(); ps != dmn->species_p.end(); ++ps)
              if (prb->specie.compare((**ps).name) == 0)
              {
                speciep = (*ps);
                value = &(speciep->density_map);
                specie_exists = true;
              }

            if (!specie_exists)
            {
              LOG_S(WARNING) << "There is no particles specie ``"
                             << prb->specie
                             << "''. Probe dump is impossible";
              continue;
            }
          }
//////////////////////////////////////////////////////////////////////
          else if (prb->component.compare("energy") == 0)
          {
            SpecieP *speciep;
            bool specie_exists = false;

            for (auto ps = dmn->species_p.begin(); ps != dmn->species_p.end(); ++ps)
              if (prb->specie.compare((**ps).name) == 0)
              {
                speciep = (*ps);
                value = &(speciep->energy_map);
                specie_exists = true;
              }

            if (!specie_exists)
            {
              LOG_S(WARNING) << "There is no particles specie ``"
                             << prb->specie
                             << "''. Probe dump is impossible";
              continue;
            }
          }
          else if (prb->component.compare("momentum/r") == 0)
          {
            SpecieP *speciep;
            bool specie_exists = false;

            for (auto ps = dmn->species_p.begin(); ps != dmn->species_p.end(); ++ps)
              if (prb->specie.compare((**ps).name) == 0)
              {
                speciep = (*ps);
                value = &(speciep->momentum_map[0]);
                specie_exists = true;
              }

            if (!specie_exists)
            {
              LOG_S(WARNING) << "There is no particles specie ``"
                             << prb->specie
                             << "''. Probe dump is impossible";
              continue;
            }
          }
          else if (prb->component.compare("momentum/phi") == 0)
          {
            SpecieP *speciep;
            bool specie_exists = false;

            for (auto ps = dmn->species_p.begin(); ps != dmn->species_p.end(); ++ps)
              if (prb->specie.compare((**ps).name) == 0)
              {
                speciep = (*ps);
                value = &(speciep->momentum_map[1]);
                specie_exists = true;
              }

            if (!specie_exists)
            {
              LOG_S(WARNING) << "There is no particles specie ``"
                             << prb->specie
                             << "''. Probe dump is impossible";
              continue;
            }
          }
          else if (prb->component.compare("momentum/z") == 0)
          {
            SpecieP *speciep;
            bool specie_exists = false;

            for (auto ps = dmn->species_p.begin(); ps != dmn->species_p.end(); ++ps)
              if (prb->specie.compare((**ps).name) == 0)
              {
                speciep = (*ps);
                value = &(speciep->momentum_map[2]);
                specie_exists = true;
              }

            if (!specie_exists)
            {
              LOG_S(WARNING) << "There is no particles specie ``"
                             << prb->specie
                             << "''. Probe dump is impossible";
              continue;
            }
          }
//////////////////////////////////////////////////////////////////////
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
    vector<size_t> prb_size = {prb->dims[0], prb->dims[1], prb->dims[2], prb->dims[3]};
    vector<size_t> offset = {0, 0};

#ifdef ENABLE_HDF5
    vector<size_t> hdf5_prb_size;

    switch (prb->shape)
    {
    case 0:
      hdf5_prb_size = {prb->dims[2] - prb->dims[0], prb->dims[3] - prb->dims[1]};
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
      //// print header every 30 values
      if (print_progress_table)
        msg::print_header (30, time);

      string shape_name;
      size_t slices = (size_t)(ceil(current_time_step / prb->schedule));
      vector<unsigned int> probe_size;

      switch (prb->shape)
      {
      case 0:
        shape_name = "rec";
        probe_size.push_back(prb->dims[2] - prb->dims[0]);
        probe_size.push_back(prb->dims[3] - prb->dims[1]);
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
      ////
      // calculate temperature and/or density before dump
      for (unsigned int r = 0; r < geometry->domains_amount[0]; ++r)
        for (unsigned int z = 0; z < geometry->domains_amount[1]; ++z)
        {
          Domain *sim_domain = smb->domains(r, z);
          SpecieP *speciep;
          bool specie_exists = false;

          for (auto ps = sim_domain->species_p.begin(); ps != sim_domain->species_p.end(); ++ps)
            if (prb->specie.compare((**ps).name) == 0)
            {
              speciep = (*ps);
              specie_exists = true;
            }

          if (!specie_exists)
            continue; // skip if there is no specie for probe

          if (prb->component.compare("density") == 0)
            speciep->calc_density();
          else if (prb->component.compare("temperature") == 0)
            speciep->calc_temperature();
          else if (prb->component.compare("energy") == 0)
            speciep->calc_energy();
          else if (prb->component.compare("momentum/r") == 0
                   || prb->component.compare("momentum/phi") == 0
                   || prb->component.compare("momentum/z") == 0)
            speciep->calc_momentum();
        }

      // overlay domains before dump
      for (unsigned int i=0; i < geometry->domains_amount[0]; i++)
        for (unsigned int j = 0; j < geometry->domains_amount[1]; j++)
        {
          Domain *sim_domain = smb->domains(i, j);
          SpecieP *speciep;
          bool specie_exists = false;

          for (auto ps = sim_domain->species_p.begin(); ps != sim_domain->species_p.end(); ++ps)
            if (prb->specie.compare((**ps).name) == 0)
            {
              speciep = (*ps);
              specie_exists = true;
            }

          if (!specie_exists)
            continue; // skip if there is no specie for probe

          // update grid
          if (i < geometry->domains_amount[0] - 1)
          {
            Domain *dst_domain = smb->domains(i+1, j);
            SpecieP *speciep_dst;
            bool specie_exists = false;

            for (auto ps = dst_domain->species_p.begin(); ps != dst_domain->species_p.end(); ++ps)
              if (prb->specie.compare((**ps).name) == 0)
              {
                speciep_dst = (*ps);
                specie_exists = true;
              }
            if (!specie_exists)
              LOG_S(FATAL) << "There is no particles specie ``"
                           << prb->specie
                           << "'' in destination domain, while overlaying for probe dump. Exiting";

            if (prb->component.compare("density") == 0)
              speciep->density_map.overlay_x(speciep_dst->density_map);
            else if (prb->component.compare("temperature") == 0)
              speciep->temperature_map.overlay_x(speciep_dst->temperature_map);
            else if (prb->component.compare("energy") == 0)
              speciep->energy_map.overlay_x(speciep_dst->energy_map);
            else if (prb->component.compare("momentum/r") == 0
                     || prb->component.compare("momentum/phi") == 0
                     || prb->component.compare("momentum/z") == 0)
            {
              speciep->momentum_map[0].overlay_x(speciep_dst->momentum_map[0]);
              speciep->momentum_map[1].overlay_x(speciep_dst->momentum_map[1]);
              speciep->momentum_map[2].overlay_x(speciep_dst->momentum_map[2]);
            }
          }

          if (j < geometry->domains_amount[1] - 1)
          {
            Domain *dst_domain = smb->domains(i, j + 1);
            SpecieP *speciep_dst;
            bool specie_exists = false;

            for (auto ps = dst_domain->species_p.begin(); ps != dst_domain->species_p.end(); ++ps)
              if (prb->specie.compare((**ps).name) == 0)
              {
                speciep_dst = (*ps);
                specie_exists = true;
              }

            if (!specie_exists)
              LOG_S(FATAL) << "There is no particles specie ``"
                           << prb->specie
                           << "'' in destination domain, while overlaying for probe dump. Exiting";

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
      if (print_progress_table)
        msg::print_values (prb->path, shape_name, prb->schedule, time);
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

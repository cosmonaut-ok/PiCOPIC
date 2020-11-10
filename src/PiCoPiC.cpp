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

// enable openmp optional
#include <signal.h>
#include <unistd.h>
#include <vector>
#include <string>
#include <atomic>         // std::atomic, std::memory_order_relaxed

#ifdef _OPENMP
#include <omp.h>
// #else
// #define omp_get_thread_num() 0
#endif

#include "defines.hpp"
#include "msg.hpp"
#include "cfg.hpp"
#include "algo/common.hpp"
#include "algo/grid.hpp"

#include "domain.hpp"
#include "dataWriter.hpp"
#include "specieP.hpp"
#include "beamP.hpp"

using namespace std;

#define BEAM_ID_START 1000

#ifdef USE_HDF5
#include "H5Cpp.h"

std::atomic<bool> lock_(false);
std::atomic<bool> recreate_writers_(false);

string hdf5_filepath;

H5::H5File *file;
#endif // USE_HDF5

void signalHandler( int signum )
{
#ifdef USE_HDF5
  if (signum == SIGUSR1 && !lock_)
  {
    LOG_S(INFO) << "Pause calculation";
    LOG_S(INFO) << "Unlocking data file...";
    file->close();
    lock_.store(true, std::memory_order_relaxed);
    LOG_S(INFO) << "Unlocked";
    cout << "Waiting for ``USR2'' OS signal to continue" << flush;
  }
  else if (signum == SIGUSR2 && lock_) // reopen
  {
    cout << endl;
    LOG_S(INFO) << "Continue calculation";
    recreate_writers_.store(true, std::memory_order_relaxed);
  }
#endif // USE_HDF5

  // exit
  if ( signum == SIGINT
       || signum == SIGQUIT
       || signum == SIGPIPE
       || signum == SIGTERM
       || signum == SIGTSTP
       || signum == SIGABRT
       || signum == SIGBUS
       || signum == SIGFPE
       || signum == SIGILL
       || signum == SIGSEGV )
  {
#ifdef USE_HDF5
    file->close();
#endif // USE_HDF5
    LOG_S(ERROR) << "Signal ``" << signum << "'' received. Exiting";
    exit(signum);
  }
}

string parse_argv_get_config(int argc, char **argv)
{
  string filename;

  if (algo::common::cmd_option_exists(argv, argv+argc, "-h"))
  {
    cerr << "USAGE:" << endl << "  picopic [ --version | -f path/to/PiCoPiC.json ]" << endl;
    exit(1);
  }

  if (algo::common::cmd_option_exists(argv, argv+argc, "--version"))
  {
    cerr << PACKAGE_NAME << " " << PACKAGE_VERSION << endl;
    exit(0);
  }

  if (algo::common::cmd_option_exists(argv, argv+argc, "-f"))
  {
    filename = algo::common::get_cmd_option(argv, argv + argc, "-f");
  }
  else
  {
    filename = std::string(PACKAGE_NAME) + std::string(".json");
  }
  if (filename.empty())
    LOG_S(FATAL) << "Can not open configuration file ``" << filename << "''";

  return filename;
}

void particles_runaway_collector (Grid<Domain*> domains, Geometry *geometry_global)
{
  // ! collects particles, that runaways from their domains and moves it to
  // ! domain, corresponding to their actual position
  // ! also, erase particles, that run out of simulation domain
  unsigned int r_domains = domains.x_size;
  unsigned int z_domains = domains.y_size;
  int j_c = 0;
  int r_c = 0;
  for (unsigned int idx = 0; idx < 2; ++idx)
    for (unsigned int idy = 0; idy < 2; ++idy)
    {
#pragma omp parallel for
      for (unsigned int i = idx; i < r_domains; i+=2)
        for (unsigned int j = idy; j < z_domains; j+=2)
        {
          Domain *sim_domain = domains(i, j);

          for (auto ps = sim_domain->species_p.begin(); ps != sim_domain->species_p.end(); ++ps)
          {
            (**ps).particles.erase(
              std::remove_if(
                (**ps).particles.begin(), (**ps).particles.end(),
                [ &j_c, &r_c, &ps, &domains, &sim_domain, &i, &j,
                  &geometry_global ] ( vector <double> * & o )
                {
                  bool res = false;

                  // not unsigned, because it could be less, than zero
                  int r_cell = P_CELL_R((*o));
                  int z_cell = P_CELL_Z((*o));

                  unsigned int i_dst = (unsigned int)ceil(r_cell / sim_domain->geometry.r_grid_amount);
                  unsigned int j_dst = (unsigned int)ceil(z_cell / sim_domain->geometry.z_grid_amount);

                  if (r_cell < 0 || z_cell < 0)
                  {
                    LOG_S(ERROR) << "Particle's position is less, than 0. Position is: ["
                                 << P_VEL_R((*o)) << ", "
                                 << P_VEL_Z((*o)) << "]. Removing";
                    ++r_c;
                    res = true;
                  }

                  else if (r_cell >= geometry_global->r_grid_amount)
                  {
                    if ((**ps).id >= BEAM_ID_START)
                    {
                      LOG_S(MAX) << "Beam particle is out of simulation domain: ["
                                 << P_POS_R((*o)) << ", "
                                 << P_POS_Z((*o)) << "]. Removing";
                    }
                    else
                    {
                      LOG_S(ERROR) << "Particle's r-position is more, than geometry r-size: "
                                   << geometry_global->r_grid_amount
                                   << ". Position is: ["
                                   << P_POS_R((*o)) << ", "
                                   << P_POS_Z((*o)) << "]. Removing";
                    }
                    ++r_c;
                    res = true;
                  }

                  // remove out-of-simulation particles
                  else if (z_cell >= geometry_global->z_grid_amount)
                  {
                    if ((**ps).id >= BEAM_ID_START)
                    {
                      LOG_S(MAX) << "Beam particle is out of simulation domain: ["
                                 << P_POS_R((*o)) << ", "
                                 << P_POS_Z((*o)) << "]. Removing";
                    }
                    else
                    {
                      LOG_S(ERROR) << "Particle's z-position is more, than geometry z-size: "
                                   << geometry_global->z_grid_amount
                                   << ". Position is: ["
                                   << P_POS_R((*o)) << ", "
                                   << P_POS_Z((*o)) << "]. Removing";
                    }
                    ++r_c;
                    res = true;
                  }

                  // move particles between cells
                  else if (i_dst != i || j_dst != j) // check that destination domain is different, than source
                  {
                    ++j_c;

                    Domain *dst_domain = domains(i_dst, j_dst);
                    for (auto pd = dst_domain->species_p.begin(); pd != dst_domain->species_p.end(); ++pd)
                      if ((**pd).id == (**ps).id)
                      {
                        LOG_S(MAX) << "Particle with specie "
                                   << (**ps).id
                                   << " jump from domain "
                                   << i << "," << j
                                   << " to domain "
                                   << i_dst << "," << j_dst;
                        (**pd).particles.push_back(o);
                      }

                    res = true;
                  }
                  return res;
                }),
              (**ps).particles.end());
          }
        }
    }
  if (j_c > 0)
    LOG_S(MAX) << "Amount of particles to jump between domains: " << j_c;
  if (r_c > 0)
    LOG_S(MAX) << "Amount of particles to remove: " << r_c;
}

void current_overlay (Grid<Domain*> domains, Geometry *geometry_global)
{
  unsigned int r_domains = domains.x_size;
  unsigned int z_domains = domains.y_size;

  for (unsigned int idx = 0; idx < 2; ++idx)
    for (unsigned int idy = 0; idy < 2; ++idy)
    {
#pragma omp parallel for
      for (unsigned int i = idx; i < r_domains; i+=2)
        for (unsigned int j = idy; j < z_domains; j+=2)
        {
          Domain *sim_domain = domains(i, j);

          // update grid
          if (i < geometry_global->domains_by_r - 1)
          {
            Domain *dst_domain = domains(i+1, j);
            sim_domain->current->current.overlay_x(dst_domain->current->current);
          }

          if (j < geometry_global->domains_by_z - 1)
          {
            Domain *dst_domain = domains(i, j + 1);
            sim_domain->current->current.overlay_y(dst_domain->current->current);
          }

          // if (i < geometry_global->domains_by_r - 1 && j < geometry_global->domains_by_z - 1)
          // {
          //   Domain *dst_domain = domains(i + 1, j + 1);
          //   sim_domain->current->current.overlay_xy(dst_domain->current->current);
          // }
        }
    }
}

void field_h_overlay (Grid<Domain*> domains, Geometry *geometry_global)
{
  unsigned int r_domains = domains.x_size;
  unsigned int z_domains = domains.y_size;

  for (unsigned int idx = 0; idx < 2; ++idx)
    for (unsigned int idy = 0; idy < 2; ++idy)
    {
#pragma omp parallel for
      for (unsigned int i = idx; i < r_domains; i+=2)
        for (unsigned int j = idy; j < z_domains; j+=2)
        {
          Domain *sim_domain = domains(i, j);

          // update grid
          if (i < geometry_global->domains_by_r - 1)
          {
            Domain *dst_domain = domains(i+1, j);
            sim_domain->maxwell_solver->field_h.overlay_x(
              dst_domain->maxwell_solver->field_h
              );
            sim_domain->maxwell_solver->field_h_at_et.overlay_x(
              dst_domain->maxwell_solver->field_h_at_et
              );
          }

          if (j < geometry_global->domains_by_z - 1)
          {
            Domain *dst_domain = domains(i, j + 1);
            sim_domain->maxwell_solver->field_h.overlay_y(
              dst_domain->maxwell_solver->field_h);
            sim_domain->maxwell_solver->field_h_at_et.overlay_y(
              dst_domain->maxwell_solver->field_h_at_et);
          }

          // if (i < geometry_global->domains_by_r - 1 && j < geometry_global->domains_by_z - 1)
          // {
          //   Domain *dst_domain = domains(i + 1, j + 1);
          //   sim_domain->field_h->field.overlay_xy(dst_domain->field_h->field);
          //   sim_domain->field_h->field_at_et.overlay_xy(dst_domain->field_h->field_at_et);
          // }
        }
    }
}

void field_e_overlay (Grid<Domain*> domains, Geometry *geometry_global)
{
  unsigned int r_domains = domains.x_size;
  unsigned int z_domains = domains.y_size;

  for (unsigned int idx = 0; idx < 2; ++idx)
    for (unsigned int idy = 0; idy < 2; ++idy)
    {
#pragma omp parallel for
      for (unsigned int i = idx; i < r_domains; i+=2)
        for (unsigned int j = idy; j < z_domains; j+=2)
        {
          Domain *sim_domain = domains(i, j);

          // update grid
          if (i < geometry_global->domains_by_r - 1)
          {
            Domain *dst_domain = domains(i+1, j);
            sim_domain->maxwell_solver->field_e.overlay_x(
              dst_domain->maxwell_solver->field_e
              );
          }

          if (j < geometry_global->domains_by_z - 1)
          {
            Domain *dst_domain = domains(i, j + 1);
            sim_domain->maxwell_solver->field_e.overlay_y(
              dst_domain->maxwell_solver->field_e
              );
          }

          // if (i < geometry_global->domains_by_r - 1 && j < geometry_global->domains_by_z - 1)
          // {
          //   Domain *dst_domain = domains(i + 1, j + 1);
          //   sim_domain->field_e->field.overlay_xy(dst_domain->field_e->field);
          // }
        }
    }
}

int main(int argc, char **argv)
{
#ifdef USE_HDF5
  lock_.store(false, std::memory_order_relaxed);
  recreate_writers_.store(false, std::memory_order_relaxed);
#endif

  // handle signals
  signal(SIGUSR1, signalHandler);
  signal(SIGUSR2, signalHandler);
  signal(SIGTERM, signalHandler); // kill
  signal(SIGINT, signalHandler);  // control+c
  signal(SIGQUIT, signalHandler);
  signal(SIGPIPE, signalHandler);
  signal(SIGTSTP, signalHandler);
  signal(SIGABRT, signalHandler);
  signal(SIGBUS, signalHandler);
  signal(SIGFPE, signalHandler);
  signal(SIGILL, signalHandler);
  signal(SIGSEGV, signalHandler);


  // do not handle signals by Loguru
  static loguru::SignalOptions s_options = loguru::SignalOptions();
  s_options.unsafe_signal_handler = false;
  s_options.sigabrt = false;
  s_options.sigbus = false;
  s_options.sigfpe = false;
  s_options.sigill = false;
  s_options.sigint = false;
  s_options.sigsegv = false;
  s_options.sigterm = false;

  static loguru::Options options = loguru::Options();
  options.signals = s_options;

  loguru::init(argc, argv, options);

#ifdef _OPENMP
#ifdef OPENMP_DYNAMIC_THREADS
  omp_set_dynamic(1); // Explicitly enable dynamic teams
  LOG_S(MAX) << "Number of Calculation Processors Changing Dynamically";
#else
  int cores = omp_get_num_procs();
  LOG_S(MAX) << "Number of Calculation Processors: " << cores;
  omp_set_dynamic(0); // Explicitly disable dynamic teams
  omp_set_num_threads(cores); // Use 4 threads for all consecutive parallel regions
#endif
#else
  LOG_S(MAX) << "There is no Openmp Here";
#endif

  ////!
  ////! Program begin
  ////!
  LOG_S(INFO) << "Initialization";

  // string cfgname;
  string cfgname = parse_argv_get_config(argc, argv);

  LOG_S(MAX) << "Reading Configuration File ``" << cfgname << "''";
  Cfg cfg = Cfg(cfgname);

  Geometry* geometry_global = cfg.geometry;

  TimeSim* sim_time_clock = cfg.time;

  LOG_S(MAX) << "Initializing Geometry, Particle Species and Simulation Domains";

  Grid<Domain*> domains (geometry_global->domains_by_r, geometry_global->domains_by_z, 0);

  LOG_S(MAX) << "Initializing Data Paths";

#ifdef USE_HDF5
  std::string path = cfg.output_data->data_root;
  std::string datafile_name = "data.h5";

  // make root directory
  algo::common::make_directory(path);

  hdf5_filepath = path + "/" + datafile_name;

  try
  {
    H5::Exception::dontPrint();
    file = new H5::H5File(hdf5_filepath.c_str(), H5F_ACC_EXCL);
    // uncomment to truncate if exists
    // file = new H5::H5File(hdf5_filepath.c_str(), H5F_ACC_TRUNC);
  }
  catch (const H5::FileIException& error)
  {
    error.printErrorStack();
    LOG_S(FATAL) << "Can not create new file ``" << hdf5_filepath << "'', file already exists. Exiting";
  }
#endif // USE_HDF5


  vector<DataWriter> data_writers;

  for (auto i = cfg.probes.begin(); i != cfg.probes.end(); ++i)
  {
    int probe_size[4] = {i->r_start, i->z_start, i->r_end, i->z_end};

    DataWriter writer (cfg.output_data->data_root, i->component,
                       i->specie, i->shape, probe_size, i->schedule,
                       cfg.output_data->compress, cfg.output_data->compress_level,
                       geometry_global, sim_time_clock, domains, cfg.cfg2str());


#ifdef USE_HDF5
    writer.hdf5_file = file; // pass pointer to opened HDF5 file to outEngineHDF5's object
    writer.hdf5_init(cfg.cfg2str());
#endif // USE_HDF5

    data_writers.push_back(writer);
  }

  unsigned int r_domains = geometry_global->domains_by_r;
  unsigned int z_domains = geometry_global->domains_by_z;

  unsigned int p_id_counter = 0;
  unsigned int b_id_counter = BEAM_ID_START;

  for (unsigned int i=0; i < r_domains; i++)
    for (unsigned int j = 0; j < z_domains; j++)
    {
      //! init geometry
      bool wall_r0 = false;
      bool wall_rr = false;
      bool wall_z0 = false;
      bool wall_zz = false;

      double pml_l_z0 = 0;
      double pml_l_zwall = 0;
      double pml_l_rwall = 0;

      // set walls to domains
      if (i == 0)
        wall_r0 = true;
      if (i == r_domains - 1)
        wall_rr = true;
      if (j == 0)
        wall_z0 = true;
      if (j == z_domains - 1)
        wall_zz = true;

#ifdef USE_PML
      // set PML to domains
      if (geometry_global->r_size - geometry_global->r_size / r_domains * (i + 1)
          < geometry_global->pml_length[2])
        pml_l_rwall = geometry_global->pml_length[2];
      if (j * geometry_global->z_size / z_domains < geometry_global->pml_length[1])
        pml_l_z0 = geometry_global->pml_length[1];
      if (geometry_global->z_size - geometry_global->z_size / z_domains * (j + 1)
          < geometry_global->pml_length[3])
        pml_l_zwall = geometry_global->pml_length[3];
#endif

      unsigned int bot_r = (unsigned int)geometry_global->r_grid_amount * i / r_domains;
      unsigned int top_r = (unsigned int)geometry_global->r_grid_amount * (i + 1) / r_domains;

      unsigned int left_z = (unsigned int)geometry_global->z_grid_amount * j / z_domains;
      unsigned int right_z = (unsigned int)geometry_global->z_grid_amount * (j + 1) / z_domains;

#ifdef USE_PML
      Geometry *geom_domain = new Geometry (
        geometry_global->r_size / r_domains,
        geometry_global->z_size / z_domains,
        bot_r, top_r, left_z, right_z,
        pml_l_z0 * z_domains,    // multiplying is a workaround, because domain
        pml_l_zwall * z_domains, // doesn't know abount whole simulation domain size
        pml_l_rwall * r_domains, // aka geometry_global
        geometry_global->pml_sigma[0],
        geometry_global->pml_sigma[1],
        wall_r0,
        wall_z0,
        wall_rr,
        wall_zz
        );

      // WORKAROUND: // used just to set PML
      // for information about global geometry
      // WARNING! don't use it in local geometries!
      geom_domain->domains_by_r = r_domains;
      geom_domain->domains_by_z = z_domains;
      // /WORKAROUND
#else
      Geometry *geom_domain = new Geometry (
        geometry_global->r_size / r_domains,
        geometry_global->z_size / z_domains,
        bot_r, top_r, left_z, right_z,
        wall_r0,
        wall_z0,
        wall_rr,
        wall_zz
        );
#endif

      // init particle species
      vector<SpecieP *> species_p;

      p_id_counter = 0; // counter for particle specie IDs
      b_id_counter = 1000; // counter for particle specie IDs

      for (auto k = cfg.particle_species.begin(); k != cfg.particle_species.end(); ++k)
      {
        unsigned int grid_cell_macro_amount = (int)(k->macro_amount / r_domains / z_domains);

        double drho_by_dz = (k->right_density - k->left_density) / geometry_global->z_size;
        double ld_local = k->left_density + drho_by_dz * left_z * geom_domain->z_cell_size;
        double rd_local = k->left_density + drho_by_dz * right_z * geom_domain->z_cell_size;

        SpecieP *pps = new SpecieP (p_id_counter,
                                    k->name,
                                    k->charge, k->mass, grid_cell_macro_amount,
                                    ld_local, rd_local,
                                    k->temperature, geom_domain, sim_time_clock);
        species_p.push_back(pps);

        ++p_id_counter;
      };

      // init particle beams
      if (! cfg.particle_beams.empty())
        for (auto bm = cfg.particle_beams.begin(); bm != cfg.particle_beams.end(); ++bm)
        {
          if (bm->bunches_amount > 0)
          {
            BeamP *beam = new BeamP (b_id_counter, ((string)"beam_").append(bm->name),
                                     bm->charge, bm->mass, bm->macro_amount,
                                     bm->start_time, bm->bunch_radius, bm->density,
                                     bm->bunches_amount, bm->bunch_length,
                                     bm->bunches_distance, bm->velocity,
                                     geom_domain, sim_time_clock);

            species_p.push_back(beam);
          }
          ++b_id_counter;
        }

      Domain *sim_domain = new Domain(*geom_domain, species_p, sim_time_clock);
      domains.set(i, j, sim_domain);
    };

  LOG_S(INFO) << "Preparation to calculation";

#pragma omp parallel for collapse(2)
  for (unsigned int i=0; i < r_domains; i++)
    for (unsigned int j = 0; j < z_domains; j++)
    {
      Domain *sim_domain = domains(i, j);

      sim_domain->distribute(); // spatial and velocity distribution
    }

  //! Main calculation loop
  LOG_S(INFO) << "Launching calculation";

  // need to set p_id_counter to zero, because bunches are injecting dynamically
  // and we need mark it
  p_id_counter = 0;

  while (sim_time_clock->current <= sim_time_clock->end)
  {
    //! Steps:
    //! efield->calc_field
    //! == efield overlay ==
    //! 2. hfield->calc_field
    //! == hfield overlay ==
    //! 3. c_current->reset_j
    //! 1. init beam
    //! 6. push_particles
    //! 7. dump_position_to_old
    //! 8. half_step_pos
    //! 9. back_position_to_rz
    //! == runaway ==
    //! 10. azimuthal_current_distribution
    //! 11. half_step_pos
    //! 13. back_position_to_rz
    //! 12. reflection
    //! == runaway ==
    //! == current overlay ==
    //! 14. current_distribution
    //! 15. back_velocity_to_rz
    //! == current overlay ==

    LOG_S(MAX) << "Run calculation loop for timestamp: " << sim_time_clock->current;

//// inject beam
    if (! cfg.particle_beams.empty())
    {
#pragma omp parallel for collapse(2)
      for (unsigned int i=0; i < r_domains; i++)
        for (unsigned int j = 0; j < z_domains; j++)
        {
          Domain *sim_domain = domains(i, j);

          // ! 1. manage beam
          sim_domain->manage_beam();
        }
    }

//// solve maxwell equations
#pragma omp parallel for collapse(2)
    for (unsigned int i=0; i < r_domains; i++)
      for (unsigned int j = 0; j < z_domains; j++)
      {
        Domain *sim_domain = domains(i, j);

        // ! 5. Calculate electric field (E)
        sim_domain->weight_field_e(); // +
      }
    field_e_overlay(domains, geometry_global);

#pragma omp parallel for collapse(2)
    for (unsigned int i=0; i < r_domains; i++)
      for (unsigned int j = 0; j < z_domains; j++)
      {
        Domain *sim_domain = domains(i, j);

        // ! 2. Calculate magnetic field (H)
        sim_domain->weight_field_h(); // +
      }
    field_h_overlay(domains, geometry_global);

//// advance particles
#pragma omp parallel for collapse(2)
    for (unsigned int i=0; i < r_domains; i++)
      for (unsigned int j = 0; j < z_domains; j++)
      {
        Domain *sim_domain = domains(i, j);

        // ! 3. Calculate velocity
        sim_domain->push_particles(); // +

#ifdef COLLISIONS
        sim_domain->collide(); // collide before reflect
#endif

        sim_domain->dump_particle_positions_to_old(); // +
        sim_domain->update_particles_coords(); // + +reflect
        sim_domain->particles_back_position_to_rz(); // +

        sim_domain->reflect();

        sim_domain->particles_back_velocity_to_rz();

        sim_domain->bind_cell_numbers();
      }
    particles_runaway_collector(domains, geometry_global);

//// solve currents
#pragma omp parallel for collapse(2)
    for (unsigned int i=0; i < r_domains; i++)
      for (unsigned int j = 0; j < z_domains; j++)
      {
        Domain *sim_domain = domains(i, j);

        sim_domain->reset_current();
        sim_domain->weight_current();
      }
    current_overlay(domains, geometry_global);

//// output data
    // for (unsigned int i=0; i < data_writers.size(); ++i)
    for (auto i = data_writers.begin(); i != data_writers.end(); i++)
      (*i)();

    sim_time_clock->current += sim_time_clock->step;

//// check if the simulation pause/unpause requested
    // (signals USR1 for pause and USR2 for unpause
#ifdef USE_HDF5
    while (lock_)
    {
      if (recreate_writers_)
      {
        // reopen HDF file initially
        try
        {
	  LOG_S(INFO) << "Locking data file...";
          H5::Exception::dontPrint();
          file = new H5::H5File(hdf5_filepath.c_str(), H5F_ACC_RDWR);
	  LOG_S(INFO) << "Locked";
        }
        catch (const H5::FileIException&)
        {
	  LOG_S(FATAL) << "Can not open data file ``" << hdf5_filepath << "''. Not accessible, or corrupted";
        }

        // recreate all of the writers
        data_writers.clear();

        for (auto i = cfg.probes.begin(); i != cfg.probes.end(); ++i)
        {
          int probe_size[4] = {i->r_start, i->z_start, i->r_end, i->z_end};

          DataWriter writer (cfg.output_data->data_root, i->component,
                             i->specie, i->shape, probe_size, i->schedule,
                             cfg.output_data->compress, cfg.output_data->compress_level,
                             geometry_global, sim_time_clock, domains, cfg.cfg2str());

          writer.hdf5_file = file; // pass pointer to opened HDF5 file to outEngineHDF5's object
          writer.hdf5_init(cfg.cfg2str());

          data_writers.push_back(writer);
        }
	lock_.store(false, std::memory_order_relaxed);
	recreate_writers_.store(false, std::memory_order_relaxed);
      }
      else
      {
        cout << "." << flush;
      }
      sleep(1);
    }
#endif

  }

#ifdef USE_HDF5
  file->close();
#endif

#ifndef DEBUG
  MSG("+------------+-------------+-------------------------------------------------+-------+------------------+---------------------+");
#endif

  LOG_S(INFO) << "SIMULATION COMPLETE";

  return 0;
}

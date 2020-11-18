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

#include "SMB.hpp"
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

void signal_handler( int signum )
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

int main(int argc, char **argv)
{
#ifdef USE_HDF5
  lock_.store(false, std::memory_order_relaxed);
  recreate_writers_.store(false, std::memory_order_relaxed);
#endif

  // handle signals
  signal(SIGUSR1, signal_handler);
  signal(SIGUSR2, signal_handler);
  signal(SIGTERM, signal_handler); // kill
  signal(SIGINT, signal_handler);  // control+c
  signal(SIGQUIT, signal_handler);
  signal(SIGPIPE, signal_handler);
  signal(SIGTSTP, signal_handler);
  signal(SIGABRT, signal_handler);
  signal(SIGBUS, signal_handler);
  signal(SIGFPE, signal_handler);
  signal(SIGILL, signal_handler);
  signal(SIGSEGV, signal_handler);


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

// #ifdef _OPENMP
// #ifdef OPENMP_DYNAMIC_THREADS
//   omp_set_dynamic(1); // Explicitly enable dynamic teams
//   LOG_S(MAX) << "Number of Calculation Processors Changing Dynamically";
// #else
//   int cores = omp_get_num_procs();
//   LOG_S(MAX) << "Number of Calculation Processors: " << cores;
//   omp_set_dynamic(0); // Explicitly disable dynamic teams
//   omp_set_num_threads(cores); // Use 4 threads for all consecutive parallel regions
// #endif
// #else
//   LOG_S(MAX) << "There is no Openmp Here";
// #endif

  ////!
  ////! Program begin
  ////!
  LOG_S(INFO) << "Initialization";

  // string cfgname;
  string cfgname = parse_argv_get_config(argc, argv);

  LOG_S(MAX) << "Reading Configuration File ``" << cfgname << "''";
  Cfg cfg = Cfg(cfgname);

  LOG_S(MAX) << "Initializing Geometry, Time and Simulation Domains";
  Geometry* geometry_global = cfg.geometry;

  TimeSim* sim_time_clock = cfg.time;

  // define a shared memory block
  SMB shared_mem_blk ( &cfg, geometry_global, sim_time_clock);

  Grid<Domain *> domains = shared_mem_blk.domains;

  unsigned int r_domains = geometry_global->domains_by_r;
  unsigned int z_domains = geometry_global->domains_by_z;

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

  LOG_S(INFO) << "Preparation to calculation";

  shared_mem_blk.distribute();
// #pragma omp parallel for collapse(2)
//   for (unsigned int i=0; i < r_domains; i++)
//     for (unsigned int j = 0; j < z_domains; j++)
//     {
//       Domain *sim_domain = domains(i, j);

//       sim_domain->distribute(); // spatial and velocity distribution
//     }

  //! Main calculation loop
  LOG_S(INFO) << "Launching calculation";

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
    shared_mem_blk.inject_beam();

//// solve maxwell equations
    shared_mem_blk.solve_maxvell();

//// advance particles
    shared_mem_blk.advance_particles();

//// solve currents
    shared_mem_blk.solve_current();

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

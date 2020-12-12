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

#include "defines.hpp"

#ifdef ENABLE_MPI
#include <mpi.h>
using namespace MPI;
#endif // ENABLE_MPI

#include "msg.hpp"
#include "cfg.hpp"
#include "algo/common.hpp"
#include "algo/grid.hpp"

#include "SMB.hpp"
#include "outController.hpp"
#include "specieP.hpp"
#include "beamP.hpp"

using namespace std;

#ifdef ENABLE_HDF5
#include <highfive/H5File.hpp>

#define HDF5_DATA_FILE "data.h5"

std::atomic<bool> lock_(false);
std::atomic<bool> recreate_writers_(false);

string hdf5_filepath;

// using namespace HighFive;

HighFive::File *file;
#endif // ENABLE_HDF5

void signal_handler( int signum )
{
#ifdef ENABLE_HDF5
  if (signum == SIGUSR1 && !lock_)
  {
    LOG_S(INFO) << "Pause calculation";
    LOG_S(INFO) << "Unlocking data file...";
    delete file;
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
#endif // ENABLE_HDF5

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
#ifdef ENABLE_HDF5
    delete file;
#endif // ENABLE_HDF5

    const char * sigstr = strsignal(signum);
    std::cerr << "Signal ``" << strdup(sigstr) << " (" << signum << ")'' received. Exiting" << std::endl;
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
#ifdef ENABLE_HDF5
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
  bool print_progress_table = true;
#ifdef ENABLE_MPI
  try
  {
    // initialize MPI
    Init ();
    // COMM_WORLD.Set_errhandler (MPI::ERRORS_THROW_EXCEPTIONS);

    // get information about our world
    const int NPROCS = COMM_WORLD.Get_size ();
    const int ID = COMM_WORLD.Get_rank ();

    // number of tasks
    LOG_S(INFO) << "Number of MPI processes is: " << NPROCS;

    // size of each task
    const size_t SZ = 1024 * 1024;

    if (ID != 0)
      print_progress_table = false;

    Geometry* geometry_smb = geometry_global;

    // update macro amount per task
    cfg.macro_amount /= NPROCS;
    for (auto k = cfg.particle_species.begin(); k != cfg.particle_species.end(); ++k)
      k->macro_amount /= NPROCS;

    // update cell dimensions
    size_t smb_cell_amount = geometry_global->cell_amount[1] / NPROCS;
    double smb_z_size = geometry_global->size[1] / NPROCS;

    size_t smb_cell_begin = smb_cell_amount * ID;
    size_t smb_cell_end = smb_cell_amount * ( ID + 1 );

    geometry_smb->cell_amount[1] = smb_cell_amount;
    geometry_smb->cell_dims[1] = smb_cell_begin;
    geometry_smb->cell_dims[3] = smb_cell_end;
    geometry_smb->size[1] = smb_z_size;

    // remove walls for SMBs at requirement
    if (ID < NPROCS - 1)
      geometry_smb->walls[3] = false;
    if (ID > 0)
      geometry_smb->walls[1] = false;

    // define a shared memory block
    SMB shared_mem_blk ( &cfg, geometry_smb, sim_time_clock, ID, NPROCS);
#else
    // define a shared memory block
    SMB shared_mem_blk ( &cfg, geometry_global, sim_time_clock, 0, 1);
#endif // ENABLE_MPI

    LOG_S(MAX) << "Initializing Data Paths";

#ifdef ENABLE_HDF5
    // suppress traditional HDF5's flooding error output
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    std::string path = cfg.output_data->data_root;
    std::string datafile_name = HDF5_DATA_FILE;

    // make root directory
    algo::common::make_directory(path);

    hdf5_filepath = path + "/" + datafile_name;

    try
    {
#ifdef ENABLE_MPI
      file = new HighFive::File (
        hdf5_filepath.c_str(),
        HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Excl,
        HighFive::MPIOFileDriver(COMM_WORLD, MPI_INFO_NULL)
        );
#else
      file = new HighFive::File (
        hdf5_filepath.c_str(),
        HighFive::File::Create | HighFive::File::Excl
        );
#endif // ENABLE_MPI
    }
    catch (const HighFive::FileException& error)
    {
      LOG_S(FATAL) << "Can not create new file ``"
                   << hdf5_filepath << "''"
                   << endl << error.what();
    }
#endif // ENABLE_HDF5

#ifdef ENABLE_HDF5
    OutController out_controller ( file, geometry_global, sim_time_clock,
                                   cfg.probes, &shared_mem_blk, cfg.cfg2str(),
                                   print_progress_table );

    out_controller.hdf5_file = file;
#else
    OutController out_controller ( geometry_global, sim_time_clock,
                                   cfg.probes, &shared_mem_blk, cfg.cfg2str(),
                                   print_progress_table );
#endif

    LOG_S(INFO) << "Preparation to calculation";

    shared_mem_blk.distribute();

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
#ifdef ENABLE_MPI
      LOG_S(MAX) << "Launching writers for SMB with ID: ``" << ID
                 << "'' at ``" << sim_time_clock->current << "''";
#else
      LOG_S(MAX) << "Launching writers for SMB at ``"
                 << sim_time_clock->current << "''";
#endif // ENABLE_MPI
      out_controller();

      sim_time_clock->current += sim_time_clock->step;

//// check if the simulation pause/unpause requested
      // (signals USR1 for pause and USR2 for unpause
#ifdef ENABLE_HDF5
      while (lock_)
      {
        if (recreate_writers_)
        {
          // reopen HDF file initially
          try
          {
            LOG_S(INFO) << "Locking data file...";
            file = new HighFive::File(hdf5_filepath.c_str(), HighFive::File::ReadWrite | HighFive::File::Excl);
            LOG_S(INFO) << "Locked";
          }
          catch (const HighFive::Exception&)
          {
            LOG_S(FATAL) << "Can not open data file ``" << hdf5_filepath << "''. Not accessible, or corrupted";
          }

          // recreate all of the writers
          // data_writers.clear();

          // for (auto i = cfg.probes.begin(); i != cfg.probes.end(); ++i)
          // {
          //   int probe_size[4] = {i->dims[0], i->dims[1], i->dims[2], i->dims[3]};

          //   DataWriter writer (cfg.output_data->data_root, i->component,
          //                      i->specie, i->shape, probe_size, i->schedule,
          //                      cfg.output_data->compress, cfg.output_data->compress_level,
          //                      geometry_global, sim_time_clock, domains, cfg.cfg2str());

          //   writer.hdf5_file = file; // pass pointer to opened HDF5 file to outEngineHDF5's object
          //   writer.hdf5_init(cfg.cfg2str());

          //   data_writers.push_back(writer);
          // }
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

#ifdef ENABLE_HDF5
    delete file;
#endif
    if (print_progress_table)
      msg::print_final();

#ifdef ENABLE_MPI
    Finalize();
  }
  catch (MPI::Exception& e )
  {
    cerr << "MPI error: " << e.Get_error_string() << e.Get_error_code() << endl;
    MPI::COMM_WORLD.Abort (-1);
    return -1;
  }
#endif // ENABLE_MPI

  return 0;
}

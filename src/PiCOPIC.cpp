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

// enable openmp optional
#include <typeinfo>

#include "msg.hpp"

#ifdef _OPENMP
#include <omp.h>
// #else
// #define omp_get_thread_num() 0
#endif

// #ifdef USE_HDF5
// #include "ioHDF5.h"
// #endif

#include "defines.hpp"
#include "msg.hpp"
#include "cfg.hpp"
#include "lib.hpp"

#include <string>

#include "math/rand.hpp"
#include "math/maxwellJuttner.hpp"

#include "geometry.hpp"

#include "grid.hpp"

// #include "particles.hpp"

// #include "domain.hpp"
#include "timeSim.hpp"

#include "specieP.hpp"
#include "beamP.hpp"

#include "fieldE.hpp"
#include "specieP.hpp"

#include "dataWriter.hpp"

// #include "pBunch.hpp"

// #include "probePlain.hpp"

using namespace std;

#define BEAM_ID_START 1000

string parse_argv_get_config(int argc, char **argv)
{
  string filename;

  if (lib::cmd_option_exists(argv, argv+argc, "-h"))
  {
    cerr << "USAGE:" << endl << "  picopic [ --version | -f path/to/PiCOPIC.json ]" << endl;
    exit(1);
  }

  if (lib::cmd_option_exists(argv, argv+argc, "--version"))
  {
    cerr << PACKAGE_NAME << " " << PACKAGE_VERSION << endl;
    exit(0);
  }

  if (lib::cmd_option_exists(argv, argv+argc, "-f"))
  {
    filename = lib::get_cmd_option(argv, argv + argc, "-f");
    if (filename.empty())
    {
      cerr << "ERROR: configuration path is not specified" << endl;
      exit(1);
    }
  }
  else
    filename = std::string(PACKAGE_NAME) + std::string(".json");

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
              std::remove_if((**ps).particles.begin(), (**ps).particles.end(),
                             [&j_c, &r_c, &ps, &domains, &sim_domain, &i, &j, &geometry_global](vector <double> * & o)
                             {
                               bool res = false;

                               // not unsigned, because it could be less, than zero
                               int r_cell = CELL_NUMBER(P_POS_R((*o)), sim_domain->geometry.r_cell_size);
                               int z_cell = CELL_NUMBER(P_POS_Z((*o)), sim_domain->geometry.z_cell_size);

                               unsigned int i_dst = (unsigned int)ceil(r_cell / sim_domain->geometry.r_grid_amount);
                               unsigned int j_dst = (unsigned int)ceil(z_cell / sim_domain->geometry.z_grid_amount);

                               if (r_cell < 0 || z_cell < 0)
                               {
                                 LOG_ERR("Particle's position is less, than 0. Position is: ["
                                         << P_POS_R((*o)) << ", "
                                         << P_POS_Z((*o)) << "]. Removing");
                                 ++r_c;
                                 res = true;
                               }

                               if (r_cell >= geometry_global->r_grid_amount)
                               {
				 if ((**ps).id >= BEAM_ID_START)
                                 {
                                   LOG_DBG("Beam particle is out of simulation domain: ["
                                           << P_POS_R((*o)) << ", "
                                           << P_POS_Z((*o)) << "]. Removing");
                                 }
				 else
				 {
                                   LOG_ERR("Particle's r-position is more, than geometry r-size: "
                                           << geometry_global->r_grid_amount
                                           << ". Position is: ["
                                           << P_POS_R((*o)) << ", "
                                           << P_POS_Z((*o)) << "]. Removing");
				 }
                                 ++r_c;
                                 res = true;
                               }

                               // remove out-of-simulation particles
                               else if (z_cell >= geometry_global->z_grid_amount)
                               {
                                 if ((**ps).id >= BEAM_ID_START)
                                 {
                                   LOG_DBG("Beam particle is out of simulation domain: ["
                                           << P_POS_R((*o)) << ", "
                                           << P_POS_Z((*o)) << "]. Removing");
                                 }
                                 else
                                 {
                                   LOG_ERR("Particle's z-position is more, than geometry z-size: "
                                           << geometry_global->z_grid_amount
                                           << ". Position is: ["
                                           << P_POS_R((*o)) << ", "
                                           << P_POS_Z((*o)) << "]. Removing");
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
                                     LOG_DBG("Particle with specie "
                                             << (**ps).id
                                             << " jump from domain "
                                             << i << "," << j
                                             << " to domain "
                                             << i_dst << "," << j_dst);
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
  {
    LOG_DBG("Amount of particles to jump between domains: " << j_c);
  }
  if (r_c > 0)
  {
    LOG_DBG("Amount of particles to remove: " << r_c);
  }
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
            sim_domain->field_h->field.overlay_x(dst_domain->field_h->field);
            sim_domain->field_h->field_at_et.overlay_x(dst_domain->field_h->field_at_et);
          }

          if (j < geometry_global->domains_by_z - 1)
          {
            Domain *dst_domain = domains(i, j + 1);
            sim_domain->field_h->field.overlay_y(dst_domain->field_h->field);
            sim_domain->field_h->field_at_et.overlay_y(dst_domain->field_h->field_at_et);
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
            sim_domain->field_e->field.overlay_x(dst_domain->field_e->field);
          }

          if (j < geometry_global->domains_by_z - 1)
          {
            Domain *dst_domain = domains(i, j + 1);
            sim_domain->field_e->field.overlay_y(dst_domain->field_e->field);
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

#ifdef _OPENMP
#ifdef OPENMP_DYNAMIC_THREADS
  omp_set_dynamic(1); // Explicitly enable dynamic teams
  LOG_DBG("Number of Calculation Processors Changing Dynamically");
#else
  int cores = omp_get_num_procs();
  LOG_DBG("Number of Calculation Processors: " << cores);
  omp_set_dynamic(0); // Explicitly disable dynamic teams
  omp_set_num_threads(cores); // Use 4 threads for all consecutive parallel regions
#endif
#else
  LOG_DBG("There is no Openmp Here");
#endif

  ////!
  ////! Program begin
  ////!
  LOG_INFO("Initialization");

  // string cfgname;
  string cfgname = parse_argv_get_config(argc, argv);

  LOG_DBG("Reading Configuration File ``" << cfgname << "''");
  Cfg cfg = Cfg(cfgname.c_str());

  Geometry* geometry_global = cfg.geometry;

  TimeSim* sim_time_clock = cfg.time;

  LOG_DBG("Initializing Geometry, Particle Species and Simulation Domains");

  Grid<Domain*> domains (geometry_global->domains_by_r, geometry_global->domains_by_z, 0);

  LOG_DBG("Initializing Data Paths");

  vector<DataWriter> data_writers;

  for (auto i = cfg.probes.begin(); i != cfg.probes.end(); ++i)
  {
    int probe_size[4] = {i->r_start, i->z_start, i->r_end, i->z_end};

    DataWriter writer (cfg.output_data->data_root, i->component,
                       i->specie, i->shape, probe_size, i->schedule,
                       cfg.output_data->compress, cfg.output_data->compress_level,
                       geometry_global, sim_time_clock, domains, cfg.cfg2str());

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

      // set PML to domains
      if (geometry_global->r_size - geometry_global->r_size / r_domains * (i + 1)
          < geometry_global->pml_length[2])
        pml_l_rwall = geometry_global->pml_length[2];
      if (j * geometry_global->z_size / z_domains < geometry_global->pml_length[1])
        pml_l_z0 = geometry_global->pml_length[1];
      if (geometry_global->z_size - geometry_global->z_size / z_domains * (j + 1)
          < geometry_global->pml_length[3])
        pml_l_zwall = geometry_global->pml_length[3];

      unsigned int bot_r = (unsigned int)geometry_global->r_grid_amount * i / r_domains;
      unsigned int top_r = (unsigned int)geometry_global->r_grid_amount * (i + 1) / r_domains;

      unsigned int left_z = (unsigned int)geometry_global->z_grid_amount * j / z_domains;
      unsigned int right_z = (unsigned int)geometry_global->z_grid_amount * (j + 1) / z_domains;

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
      else
      {
        LOG_WARN("There is no particle beams");
      }
      Domain *sim_domain = new Domain(*geom_domain, species_p, sim_time_clock);
      domains.set(i, j, sim_domain);
    };

  LOG_INFO("Preparation to calculation");

#pragma omp parallel for collapse(2)
  for (unsigned int i=0; i < r_domains; i++)
    for (unsigned int j = 0; j < z_domains; j++)
    {
      Domain *sim_domain = domains(i, j);

      sim_domain->distribute(); // spatial and velocity distribution
    }

  LOG_DBG("Initializing Boundary Conditions (TBD)");
  // TODO: this is legacy PDP3 code. Need to update it
//   void LoadInitParam::init_boundary ()
// {
//   //! initialize boundaries and boundary conditions

//   // Maxwell initial conditions
//   BoundaryMaxwellConditions maxwell_rad(efield); // TODO: WTF?
//   maxwell_rad.specify_initial_field(params->geom,
//                                     params->boundary_maxwell_e_phi_upper,
//                                     params->boundary_maxwell_e_phi_left,
//                                     params->boundary_maxwell_e_phi_right);

//   if (params->boundary_conditions == 0)
//   {
//     p_list->charge_weighting(c_rho_new);

//     // Seems: https://en.wikipedia.org/wiki/Dirichlet_distribution
//     PoissonDirichlet dirih(params->geom);
//     dirih.poisson_solve(efield, c_rho_new);
//   }
// }

  //! Main calculation loop
  LOG_INFO("Launching calculation");

  // need to set p_id_counter to zero, because bunches are injecting dynamically
  // and we need mark it
  p_id_counter = 0;

  while (sim_time_clock->current <= sim_time_clock->end)
  {
    //! Steps:
    //! 1. init beam
    //! 2. hfield->calc_field
    //! == hfield overlay ==
    //! 3. c_current->reset_j
    //! 6. boris_pusher
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
    //! efield->calc_field
    //! == efield overlay ==

    LOG_DBG("Run calculation loop for timestamp: " << sim_time_clock->current);

#pragma omp parallel for collapse(2)
    for (unsigned int i=0; i < r_domains; i++)
      for (unsigned int j = 0; j < z_domains; j++)
      {
        Domain *sim_domain = domains(i, j);

        // ! 1. manage beam
        if (! cfg.particle_beams.empty())
          sim_domain->manage_beam();

        // ! 2. Calculate magnetic field (H)
        sim_domain->weight_field_h(); // +
      }
    field_h_overlay(domains, geometry_global);

#pragma omp parallel for collapse(2)
    for (unsigned int i=0; i < r_domains; i++)
      for (unsigned int j = 0; j < z_domains; j++)
      {
        Domain *sim_domain = domains(i, j);

        sim_domain->reset_current();
        // ! 3. Calculate velocity
        sim_domain->push_particles(); // +
        sim_domain->dump_particle_positions_to_old(); // +
        sim_domain->update_particles_coords(); // + +reflect
        sim_domain->particles_back_position_to_rz(); // +
        sim_domain->reflect();
      }
    particles_runaway_collector(domains, geometry_global);

#pragma omp parallel for collapse(2)
    for (unsigned int i=0; i < r_domains; i++)
      for (unsigned int j = 0; j < z_domains; j++)
      {
        Domain *sim_domain = domains(i, j);

        sim_domain->weight_current();
        sim_domain->particles_back_velocity_to_rz();
      }
    current_overlay(domains, geometry_global);

#pragma omp parallel for collapse(2)
    for (unsigned int i=0; i < r_domains; i++)
      for (unsigned int j = 0; j < z_domains; j++)
      {
        Domain *sim_domain = domains(i, j);

        // ! 5. Calculate electric field (E)
        sim_domain->weight_field_e(); // +
      }
    field_e_overlay(domains, geometry_global);

    // dump data
// #pragma omp parallel for
    for (unsigned int i=0; i < data_writers.size(); ++i)
      data_writers[i].go();

    sim_time_clock->current += sim_time_clock->step;
  }

  // this is only to finish pretty pringing of data_writers output
  if (! DEBUG)
  {
    MSG("+------------+-------------+-------------------------------------------------+-------+------------------+---------------------+");
  }

  LOG_INFO("SIMULATION COMPLETE");

  return 0;
}

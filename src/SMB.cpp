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

#include "SMB.hpp"

#ifdef ENABLE_MPI
#include <mpi.h>
// using namespace MPI;
#endif // ENABLE_MPI

using namespace std;

// SMB constructor
SMB::SMB ( Cfg* _cfg, Geometry *_geometry, TimeSim *_time,
           int _world_rank, int _world_size)
  : cfg(_cfg), time(_time), geometry(_geometry)
{
  LOG_S(MAX) << "Initializing Geometry, Particle Species and Simulation Domains";

  //
  // initialize domains grid
  //
  r_domains = geometry->domains_amount[0];
  z_domains = geometry->domains_amount[1];

  world_rank = _world_rank;
  world_size = _world_size;

  Grid<Domain*> _domains (r_domains, z_domains, 0);
  domains = _domains;

  //
  // initialize geometry
  //
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
        wall_r0 = geometry->walls[0];
      if (i == r_domains - 1)
        wall_rr = geometry->walls[2];
      if (j == 0)
        wall_z0 = geometry->walls[1];
      if (j == z_domains - 1)
        wall_zz = geometry->walls[3];

      // cell dimensions for domain
      size_t bot_r = (size_t)(geometry->cell_amount[0] * i / r_domains
                              + geometry->cell_dims[0]);
      size_t top_r = (size_t)(geometry->cell_amount[0] * (i + 1) / r_domains
                              + geometry->cell_dims[0]);

      size_t left_z = (size_t)(geometry->cell_amount[1] * j / z_domains
                               + geometry->cell_dims[1]);
      size_t right_z = (size_t)(geometry->cell_amount[1] * (j + 1) / z_domains
                                + geometry->cell_dims[1]);


#ifdef ENABLE_PML
      // set PML to domains
      if (geometry->global_cell_amount[0] - top_r < geometry->pml_size[2])
        pml_l_rwall = geometry->pml_size[2];

      if (geometry->global_cell_amount[1] - right_z < geometry->pml_size[3])
        pml_l_zwall = geometry->pml_size[3];

      if (left_z < geometry->pml_size[1])
        pml_l_z0 = geometry->pml_size[1];
#endif // ENABLE_PML

#ifdef ENABLE_PML
      Geometry *geom_domain = new Geometry (
        { geometry->size[0] / r_domains, geometry->size[1] / z_domains},
        { geometry->global_size[0], geometry->global_size[1] }, // global size
        { bot_r, left_z, top_r, right_z },
        { 0, pml_l_z0, pml_l_rwall, pml_l_zwall, },
        geometry->pml_sigma,
        { wall_r0, wall_z0, wall_rr, wall_zz }
        );

  // LOG_S(FATAL) << "BBB "
  //              << geom_domain->global_size[0] << " "
  //              << geom_domain->global_size[1] << " "
  //              << geometry->global_size[0] << " "
  //              << geometry->global_size[1] << " ";

      // WORKAROUND: // used just to set PML
      // for information about global geometry
      // WARNING! don't use it in local geometries!
      geom_domain->domains_amount.push_back(r_domains);
      geom_domain->domains_amount.push_back(z_domains);
      // /WORKAROUND
#else
      Geometry *geom_domain = new Geometry (
        { geometry->size[0] / r_domains,
          geometry->size[1] / z_domains },
        { bot_r, left_z, top_r, right_z },
        { wall_r0,
          wall_z0,
          wall_rr,
          wall_zz }
        );
#endif // ENABLE_PML

      //
      // initialize particle species and beams
      //
      unsigned int p_id_counter = 0;
      unsigned int b_id_counter = BEAM_ID_START;
      vector<SpecieP *> species_p;

      for (auto k = cfg->particle_species.begin(); k != cfg->particle_species.end(); ++k)
      {
        unsigned int grid_cell_macro_amount = (int)(k->macro_amount / r_domains / z_domains);

        double drho_by_dz = (k->right_density - k->left_density) / geometry->size[1];
        double ld_local = k->left_density + drho_by_dz * left_z * geom_domain->cell_size[1];
        double rd_local = k->left_density + drho_by_dz * right_z * geom_domain->cell_size[1];

        SpecieP *pps = new SpecieP (k->id,
                                    k->name,
                                    k->charge, k->mass, grid_cell_macro_amount,
                                    ld_local, rd_local,
                                    k->temperature, geom_domain, time);
        species_p.push_back(pps);

        ++p_id_counter;
      };

      // init particle beams
      if (! cfg->particle_beams.empty())
        for (auto bm = cfg->particle_beams.begin(); bm != cfg->particle_beams.end(); ++bm)
        {
          if (bm->bunches_amount > 0)
          {
            BeamP *beam = new BeamP (bm->id, bm->name,
                                     bm->charge, bm->mass, bm->macro_amount,
                                     bm->start_time, bm->bunch_radius, bm->density,
                                     bm->bunches_amount, bm->bunch_length,
                                     bm->bunches_distance, bm->velocity,
                                     geom_domain, time);

            species_p.push_back(beam);
          }
          ++b_id_counter;
        }

      Domain *sim_domain = new Domain(*geom_domain, species_p, time);
      domains.set(i, j, sim_domain);
    };

#ifdef _OPENMP
#ifdef ENABLE_OMP_DYNAMIC
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
}

void SMB::particles_runaway_collector ()
{
  // ! collects particles, that runaways from their __domains and moves it to
  // ! domain, corresponding to their actual position
  // ! also, erase particles, that run out of simulation domain

  // ! General description of send/recv particles to/from +/-1 SMB:
  // ! 1. form queues to send particles to +/-1 SMB (vectors of particle pointers)
  // ! 2. clear particles from domains (unpoint from vectors)
  // ! 3. send to +/-1 SMB
  // ! 4. recv from +/-1 SMB
  // ! 5. add new particles to local SMB's domains

  // create vectors to place particles,
  // scheduled to be sent to other SMBs

  vector< Particle * > queue_particles_plus;
  vector< Particle * > queue_particles_minus;

  int j_c = 0;
  int r_c = 0;

  Grid<Domain *> __domains = domains;
  Geometry *__geometry = geometry;

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
            (**ps).particles.erase (
              std::remove_if (
                (**ps).particles.begin(), (**ps).particles.end(),
                [ &j_c, &r_c, &ps, &__domains, &sim_domain,
                  &queue_particles_minus, &queue_particles_plus,
                  &i, &j, &__geometry ] ( Particle * & o )
                {
                  bool res = false;

                  // not unsigned, because it could be less, than zero
                  int r_cell = P_CELL_R((*o));
                  int z_cell = P_CELL_Z((*o));

                  // this is unshifted domain numbers (local for SMB)
                  unsigned int i_dst = (unsigned int)ceil (
                    ( r_cell - __geometry->cell_dims[0] ) // we should make cell numbers local for SMB
                    / sim_domain->geometry.cell_amount[0] );

                  unsigned int j_dst = (unsigned int)ceil (
                    ( z_cell - __geometry->cell_dims[1] ) // we should make cell numbers local for SMB
                    / sim_domain->geometry.cell_amount[1] );

                  if (r_cell < 0 || z_cell < 0)
                  {
                    LOG_S(ERROR) << "Particle's position is less, than 0. Position is: ["
                                 << P_VEL_R((*o)) << ", "
                                 << P_VEL_Z((*o)) << "]. Removing";
                    ++r_c;
                    res = true;
                  }
                  ////
                  //// check for conditions to send particle to previous/next SMB
                  ////
                  // send particle to previous SMB
#ifdef ENABLE_MPI
                  else if (z_cell < __geometry->cell_dims[1])
                  {
                    o->push_back((**ps).id);
                    queue_particles_minus.push_back(o);
                    res = true;
                  }
                  // send particle to next SMB
                  else if (z_cell >= __geometry->cell_dims[3])
                  {
                    o->push_back((**ps).id);
                    queue_particles_plus.push_back(o);
                    res = true;
                  }
#endif // ENABLE_MPI
                  else if (r_cell >= __geometry->cell_dims[2])
                  {
                    // ``beam_'' at the begining of the name
                    if ((**ps).name.find("beam_") == 0)
                    {
                      LOG_S(MAX) << "Beam particle is out of simulation domain: ["
                                 << P_POS_R((*o)) << ", "
                                 << P_POS_Z((*o)) << "]. Removing";
                    }
                    else
                    {
                      LOG_S(ERROR) << "Particle's r-position is more, than geometry r-size: "
                                   << __geometry->size[0]
                                   << ". Position is: ["
                                   << P_POS_R((*o)) << ", "
                                   << P_POS_Z((*o)) << "]. Removing";
                    }
                    ++r_c;
                    res = true;
                  }

                  // remove out-of-simulation particles
                  else if (z_cell >= __geometry->cell_dims[3])
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
                                   << __geometry->cell_size[1]
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

                    Domain *dst_domain = __domains(i_dst, j_dst);
                    for (auto pd = dst_domain->species_p.begin(); pd != dst_domain->species_p.end(); ++pd)
                      if ((**pd).id == (**ps).id)
                      {
                        LOG_S(MAX) << "Particle with specie ``"
                                   << (**ps).name
                                   << "'' jumps from domain ``"
                                   << i << "," << j
                                   << "'' to domain ``"
                                   << i_dst << "," << j_dst << "''";
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
    LOG_S(MAX) << "Amount of particles to jump between OMP domains: " << j_c;
  if (r_c > 0)
    LOG_S(MAX) << "Amount of particles to remove: " << r_c;

#ifdef ENABLE_MPI

  MPI_Barrier(MPI_COMM_WORLD); // FIXME: most likely it is not needed

  ////
  //// send particles to other domain
  ////

  unsigned int prtls_plus = queue_particles_plus.size();
  unsigned int prtls_minus = queue_particles_minus.size();

  if (world_rank < world_size - 1)
  {
    MPI_Send (
      /* data         = */ &prtls_plus,
      /* count        = */ 1,
      /* datatype     = */ MPI_UNSIGNED,
      /* destination  = */ world_rank + 1,
      /* tag          = */ 0,
      /* communicator = */ MPI_COMM_WORLD);

    LOG_S(MAX) << "Number of particles to be sent from MPI node ``"
               << world_rank
               << "'' to node ``" << world_rank + 1
               << "'' is ``" << prtls_plus << "''";

    if (prtls_plus > 0)
      for (unsigned prtl = 0; prtl < prtls_plus; ++prtl)
      {
        double dbuf[P_VEC_SIZE+1];
        for (unsigned int i = 0; i <= P_VEC_SIZE; ++i)
          dbuf[i] = (*queue_particles_plus[prtl])[i];

        MPI_Send (
          /* data         = */ dbuf,
          /* count        = */ P_VEC_SIZE+1,
          /* datatype     = */ MPI_DOUBLE,
          /* destination  = */ world_rank + 1,
          /* tag          = */ 0,
          /* communicator = */ MPI_COMM_WORLD);
      }
  }

  if (world_rank > 0)
  {
    MPI_Send (
      /* data         = */ &prtls_minus,
      /* count        = */ 1,
      /* datatype     = */ MPI_UNSIGNED,
      /* destination  = */ world_rank - 1,
      /* tag          = */ 0,
      /* communicator = */ MPI_COMM_WORLD);

    LOG_S(MAX) << "Number of particles to be sent from MPI node ``"
               << world_rank
               << "'' to node ``" << world_rank - 1
               << "'' is ``" << prtls_minus << "''";

    if (prtls_minus > 0)
      for (unsigned prtl = 0; prtl < prtls_minus; ++prtl)
      {
        double dbuf[P_VEC_SIZE+1];
        for (unsigned int i = 0; i <= P_VEC_SIZE; ++i)
          dbuf[i] = (*queue_particles_minus[prtl])[i];

        MPI_Send (
          /* data         = */ dbuf,
          /* count        = */ P_VEC_SIZE+1,
          /* datatype     = */ MPI_DOUBLE,
          /* destination  = */ world_rank - 1,
          /* tag          = */ 0,
          /* communicator = */ MPI_COMM_WORLD);
      }
  }

  ////
  //// receive particles from other MPI nodes
  ////
  if (world_rank > 0)
  {
    unsigned int howrcv; // how many particles to receive

    MPI_Recv (
      /* data         = */ &howrcv,
      /* count        = */ 1,
      /* datatype     = */ MPI_UNSIGNED,
      /* source       = */ world_rank - 1,
      /* tag          = */ 0,
      /* communicator = */ MPI_COMM_WORLD,
      /* status       = */ MPI_STATUS_IGNORE);

    LOG_S(MAX) << "Number of particles to be received from MPI node ``"
               << world_rank - 1
               << "'' to node ``" << world_rank
               << "'' is ``" << howrcv << "''";

    if (howrcv > 0)
      for (unsigned int i = 0; i < howrcv; ++i)
      {
        double dbuf[P_VEC_SIZE+1];

        MPI_Recv (
          /* data         = */ dbuf,
          /* count        = */ P_VEC_SIZE+1,
          /* datatype     = */ MPI_DOUBLE,
          /* source       = */ world_rank-1,
          /* tag          = */ 0,
          /* communicator = */ MPI_COMM_WORLD,
          /* status       = */ MPI_STATUS_IGNORE);

        // create vector and copy buffer there
        // vector<double> *n = new vector<double>(P_VEC_SIZE, 0);
        Particle *n = new Particle();
        for (unsigned int v=0; v < P_VEC_SIZE; ++v)
          (*n)[v] = dbuf[v];

        // get specie ID for particle
        unsigned int prtl_id = (unsigned int)dbuf[P_VEC_SIZE];

        // find proper domain for particle
        int r_cell = P_CELL_R((*n));
        int z_cell = P_CELL_Z((*n));

        // FIXME: this is a hardcode
        Domain *sim_domain = domains(0, 0);

        unsigned int i_dst = (unsigned int)ceil (
          ( r_cell - geometry->cell_dims[0] ) // we should make cell numbers local for SMB
          / sim_domain->geometry.cell_amount[0] );

        unsigned int j_dst = (unsigned int)ceil (
          ( z_cell - geometry->cell_dims[1] ) // we should make cell numbers local for SMB
          / sim_domain->geometry.cell_amount[1] );

        Domain *dst_domain = domains(i_dst, j_dst);

        // find proper specie for particle in domain
        for (auto sp = dst_domain->species_p.begin(); sp != dst_domain->species_p.end(); ++sp)
          if ((**sp).id == prtl_id)
            (**sp).particles.push_back(n);
      }
  }

  if (world_rank < world_size - 1)
  {
    unsigned int howrcv; // how many particles to receive

    MPI_Recv(
      /* data         = */ &howrcv,
      /* count        = */ 1,
      /* datatype     = */ MPI_UNSIGNED,
      /* source       = */ world_rank + 1,
      /* tag          = */ 0,
      /* communicator = */ MPI_COMM_WORLD,
      /* status       = */ MPI_STATUS_IGNORE);

    LOG_S(MAX) << "Number of particles to be received from MPI node "
               << world_rank + 1
               << " is " << howrcv;

    if (howrcv > 0)
      for (unsigned int i = 0; i < howrcv; ++i)
      {
        double dbuf[P_VEC_SIZE+1];

        MPI_Recv (
          /* data         = */ dbuf,
          /* count        = */ P_VEC_SIZE+1,
          /* datatype     = */ MPI_DOUBLE,
          /* source       = */ world_rank+1,
          /* tag          = */ 0,
          /* communicator = */ MPI_COMM_WORLD,
          /* status       = */ MPI_STATUS_IGNORE);

        // create vector and copy buffer there
        Particle *n = new Particle();
        // vector<double> *n = new vector<double>(P_VEC_SIZE, 0);
        for (unsigned int v=0; v < P_VEC_SIZE; ++v)
          (*n)[v] = dbuf[v];

        // get specie ID for particle
        unsigned int prtl_id = (unsigned int)dbuf[P_VEC_SIZE];

        // find proper domain for particle
        int r_cell = P_CELL_R((*n));
        int z_cell = P_CELL_Z((*n));

        // this is unshifted domain numbers (local for SMB)
        Domain *sim_domain = domains(0, 0);

        unsigned int i_dst = (unsigned int)ceil (
          ( r_cell - geometry->cell_dims[0] ) // we should make cell numbers local for SMB
          / sim_domain->geometry.cell_amount[0] );

        unsigned int j_dst = (unsigned int)ceil (
          ( z_cell - geometry->cell_dims[1] ) // we should make cell numbers local for SMB
          / sim_domain->geometry.cell_amount[1] );

        Domain *dst_domain = domains(i_dst, j_dst);

        // find proper specie for particle in domain
        for (auto sp = dst_domain->species_p.begin(); sp != dst_domain->species_p.end(); ++sp)
          if ((**sp).id == prtl_id)
            (**sp).particles.push_back(n);
      }
  }

  // send/receive and update domain species, then continue
  MPI_Barrier(MPI_COMM_WORLD);

  // clear temporary particle vectors after barrier
  for (auto q = queue_particles_minus.begin(); q != queue_particles_minus.end(); ++q)
    {
      (**q).clear();
      delete *q;
    }
  queue_particles_minus.clear();
  for (auto q = queue_particles_plus.begin(); q != queue_particles_plus.end(); ++q)
    {
      (**q).clear();
      delete *q;
    }
  queue_particles_plus.clear();
#endif // ENABLE_MPI
  ////
  ////
  ////
}

void SMB::current_overlay ()
{
  for (unsigned int idx = 0; idx < 2; ++idx)
    for (unsigned int idy = 0; idy < 2; ++idy)
    {
#pragma omp parallel for
      for (unsigned int i = idx; i < r_domains; i+=2)
        for (unsigned int j = idy; j < z_domains; j+=2)
        {
          Domain *sim_domain = domains(i, j);

          // update grid
          if (i < geometry->domains_amount[0] - 1)
          {
            Domain *dst_domain = domains(i+1, j);
            sim_domain->current->current.overlay_x(dst_domain->current->current);
          }

          if (j < geometry->domains_amount[1] - 1)
          {
            Domain *dst_domain = domains(i, j + 1);
            sim_domain->current->current.overlay_y(dst_domain->current->current);
          }

          // if (i < __geometry->domains_amount[0] - 1 && j < __geometry->domains_amount[1] - 1)
          // {
          //   Domain *dst_domain = __domains(i + 1, j + 1);
          //   sim_domain->current->current.overlay_xy(dst_domain->current->current);
          // }
        }
    }

#ifdef ENABLE_MPI
  ////
  //// send particles to other domain
  ////
  vector<MPICommOverlay> domains_plus;
  vector<MPICommOverlay> domains_minus;

  vector<MPICommOverlay> domains_plus_recv;
  vector<MPICommOverlay> domains_minus_recv;

  for (unsigned int i = 0; i < r_domains; ++i)
  {
    // send to +1 rank
    Domain *sim_domain_plus = domains(i, z_domains - 1);

    int domain_r_size = sim_domain_plus->current->current[0].x_size;
    int domain_r_real_size = sim_domain_plus->current->current[0].x_real_size;
    int domain_o_s = sim_domain_plus->current->current[0].o_s;
    int domain_z_real_size = sim_domain_plus->current->current[0].y_real_size;

    MPICommOverlay item_plus (domain_r_real_size, i, 0);
    MPICommOverlay item_plus_recv (domain_r_real_size, i, 0);

    for (int j = -domain_o_s; j < domain_r_size + domain_o_s; ++j)
    {
      int j_s = j + domain_o_s;
      int d_r_s = domain_z_real_size - domain_o_s * 2;

      item_plus.col_0_0[j_s] = sim_domain_plus->current->current[0](j, d_r_s  - 4);
      item_plus.col_0_1[j_s] = sim_domain_plus->current->current[1](j, d_r_s  - 4);
      item_plus.col_0_2[j_s] = sim_domain_plus->current->current[2](j, d_r_s  - 4);

      item_plus.col_1_0[j_s] = sim_domain_plus->current->current[0](j, d_r_s  - 3);
      item_plus.col_1_1[j_s] = sim_domain_plus->current->current[1](j, d_r_s  - 3);
      item_plus.col_1_2[j_s] = sim_domain_plus->current->current[2](j, d_r_s  - 3);

      item_plus.col_2_0[j_s] = sim_domain_plus->current->current[0](j, d_r_s - 2);
      item_plus.col_2_1[j_s] = sim_domain_plus->current->current[1](j, d_r_s - 2);
      item_plus.col_2_2[j_s] = sim_domain_plus->current->current[2](j, d_r_s - 2);

      item_plus.col_3_0[j_s] = sim_domain_plus->current->current[0](j, d_r_s - 1);
      item_plus.col_3_1[j_s] = sim_domain_plus->current->current[1](j, d_r_s - 1);
      item_plus.col_3_2[j_s] = sim_domain_plus->current->current[2](j, d_r_s - 1);

    }

    domains_plus.push_back( item_plus );
    domains_plus_recv.push_back( item_plus_recv );

    // send to -1 rank
    Domain *sim_domain_minus = domains(i, 0);

    domain_r_size = sim_domain_minus->current->current[0].x_size;
    domain_o_s = sim_domain_minus->current->current[0].o_s;

    MPICommOverlay item_minus (domain_r_size + 2 * domain_o_s, i, 0);
    MPICommOverlay item_minus_recv (domain_r_size + 2 * domain_o_s, i, 0);

    for (int j = -domain_o_s; j < domain_r_size + domain_o_s; ++j)
    {
      int j_s = j + domain_o_s;

      item_minus.col_0_0[j_s] = sim_domain_minus->current->current[0](j, 0 - domain_o_s);
      item_minus.col_0_1[j_s] = sim_domain_minus->current->current[1](j, 0 - domain_o_s);
      item_minus.col_0_2[j_s] = sim_domain_minus->current->current[2](j, 0 - domain_o_s);

      item_minus.col_1_0[j_s] = sim_domain_minus->current->current[0](j, 1 - domain_o_s);
      item_minus.col_1_1[j_s] = sim_domain_minus->current->current[1](j, 1 - domain_o_s);
      item_minus.col_1_2[j_s] = sim_domain_minus->current->current[2](j, 1 - domain_o_s);

      item_minus.col_2_0[j_s] = sim_domain_minus->current->current[0](j, 2 - domain_o_s);
      item_minus.col_2_1[j_s] = sim_domain_minus->current->current[1](j, 2 - domain_o_s);
      item_minus.col_2_2[j_s] = sim_domain_minus->current->current[2](j, 2 - domain_o_s);

      item_minus.col_3_0[j_s] = sim_domain_minus->current->current[0](j, 3 - domain_o_s);
      item_minus.col_3_1[j_s] = sim_domain_minus->current->current[1](j, 3 - domain_o_s);
      item_minus.col_3_2[j_s] = sim_domain_minus->current->current[2](j, 3 - domain_o_s);
    }

    domains_minus.push_back( item_minus );
    domains_minus_recv.push_back( item_minus_recv );
  }

  if (world_rank < world_size - 1)
    for (auto i = domains_plus.begin(); i != domains_plus.end(); ++i)
    {
      // LOG_S(WARNING) << "1 Sending from " << world_rank << " to " << world_rank + 1;
      MPI_Send(MPI_BOTTOM, 1, i->mpi_dtype(), world_rank + 1, 0, MPI_COMM_WORLD);
    }

  if (world_rank > 0)
     for (auto i = domains_minus_recv.begin(); i != domains_minus_recv.end(); ++i)
     {
       // LOG_S(WARNING) << "2 Receiving from " << world_rank - 1 << " at " << world_rank;
       MPI_Recv(MPI_BOTTOM, 1, i->mpi_dtype(), world_rank - 1, 0, MPI_COMM_WORLD, NULL);
     }

  if (world_rank > 0)
    for (auto i = domains_minus.begin(); i != domains_minus.end(); ++i)
    {
      // LOG_S(WARNING) << "3 Sending from " << world_rank << " to " << world_rank - 1;
      MPI_Send(MPI_BOTTOM, 1, i->mpi_dtype(), world_rank - 1, 0, MPI_COMM_WORLD);
    }

  if (world_rank < world_size - 1)
    for (auto i = domains_plus_recv.begin(); i != domains_plus_recv.end(); ++i)
    {
      // LOG_S(WARNING) << "4 Receiving from " << world_rank + 1 << " at " << world_rank;
      MPI_Recv(MPI_BOTTOM, 1, i->mpi_dtype(), world_rank + 1, 0, MPI_COMM_WORLD, NULL);
    }

  MPI_Barrier(MPI_COMM_WORLD);

  ////////////////////////////////////////////////////////////

  for (auto i = domains_minus_recv.begin(); i != domains_minus_recv.end(); ++i)
  {
    Domain *dst_domain_minus = domains(i->dst_domain[0], 0);

    int domain_r_size = dst_domain_minus->current->current[0].x_size;
    int domain_o_s = dst_domain_minus->current->current[0].o_s;
    int domain_z_size = dst_domain_minus->current->current[0].y_size;

    for (int j = -domain_o_s; j < domain_r_size + domain_o_s; ++j)
    {
      int j_s = j + domain_o_s;

      dst_domain_minus->current->current[0].inc(j, 0 - domain_o_s, i->col_0_0[j_s]);
      dst_domain_minus->current->current[1].inc(j, 0 - domain_o_s, i->col_0_1[j_s]);
      dst_domain_minus->current->current[2].inc(j, 0 - domain_o_s, i->col_0_2[j_s]);

      dst_domain_minus->current->current[0].inc(j, 1 - domain_o_s, i->col_1_0[j_s]);
      dst_domain_minus->current->current[1].inc(j, 1 - domain_o_s, i->col_1_1[j_s]);
      dst_domain_minus->current->current[2].inc(j, 1 - domain_o_s, i->col_1_2[j_s]);

      dst_domain_minus->current->current[0].inc(j, 2 - domain_o_s, i->col_2_0[j_s]);
      dst_domain_minus->current->current[1].inc(j, 2 - domain_o_s, i->col_2_1[j_s]);
      dst_domain_minus->current->current[2].inc(j, 2 - domain_o_s, i->col_2_2[j_s]);

      dst_domain_minus->current->current[0].inc(j, 3 - domain_o_s, i->col_3_0[j_s]);
      dst_domain_minus->current->current[1].inc(j, 3 - domain_o_s, i->col_3_1[j_s]);
      dst_domain_minus->current->current[2].inc(j, 3 - domain_o_s, i->col_3_2[j_s]);
    }
  }

  for (auto i = domains_plus_recv.begin(); i != domains_plus_recv.end(); ++i)
  {
    Domain *dst_domain_plus = domains(i->dst_domain[0], z_domains - 1);

    int domain_r_size = dst_domain_plus->current->current[0].x_size;
    int domain_o_s = dst_domain_plus->current->current[0].o_s;
    int domain_z_real_size = dst_domain_plus->current->current[0].y_real_size;

    for (unsigned int j = -domain_o_s; j < domain_r_size + domain_o_s; ++j)
    {
      int j_s = j + domain_o_s;
      int d_r_s = domain_z_real_size - domain_o_s * 2;

      dst_domain_plus->current->current[0].inc(j, d_r_s - 1, i->col_3_0[j_s]);
      dst_domain_plus->current->current[1].inc(j, d_r_s - 1, i->col_3_1[j_s]);
      dst_domain_plus->current->current[2].inc(j, d_r_s - 1, i->col_3_2[j_s]);

      dst_domain_plus->current->current[0].inc(j, d_r_s - 2, i->col_2_0[j_s]);
      dst_domain_plus->current->current[1].inc(j, d_r_s - 2, i->col_2_1[j_s]);
      dst_domain_plus->current->current[2].inc(j, d_r_s - 2, i->col_2_2[j_s]);

      dst_domain_plus->current->current[0].inc(j, d_r_s - 3, i->col_1_0[j_s]);
      dst_domain_plus->current->current[1].inc(j, d_r_s - 3, i->col_1_1[j_s]);
      dst_domain_plus->current->current[2].inc(j, d_r_s - 3, i->col_1_2[j_s]);

      dst_domain_plus->current->current[0].inc(j, d_r_s - 4, i->col_0_0[j_s]);
      dst_domain_plus->current->current[1].inc(j, d_r_s - 4, i->col_0_1[j_s]);
      dst_domain_plus->current->current[2].inc(j, d_r_s - 4, i->col_0_2[j_s]);
    }
  }
#endif // ENABLE_MPI
}

void SMB::field_h_overlay ()
{
  for (unsigned int idx = 0; idx < 2; ++idx)
    for (unsigned int idy = 0; idy < 2; ++idy)
    {
#pragma omp parallel for
      for (unsigned int i = idx; i < r_domains; i+=2)
        for (unsigned int j = idy; j < z_domains; j+=2)
        {
          Domain *sim_domain = domains(i, j);

          // update grid
          if (i < geometry->domains_amount[0] - 1)
          {
            Domain *dst_domain = domains(i+1, j);
            sim_domain->maxwell_solver->field_h.overlay_x(
              dst_domain->maxwell_solver->field_h
              );
            sim_domain->maxwell_solver->field_h_at_et.overlay_x(
              dst_domain->maxwell_solver->field_h_at_et
              );
          }

          if (j < geometry->domains_amount[1] - 1)
          {
            Domain *dst_domain = domains(i, j + 1);
            sim_domain->maxwell_solver->field_h.overlay_y(
              dst_domain->maxwell_solver->field_h);
            sim_domain->maxwell_solver->field_h_at_et.overlay_y(
              dst_domain->maxwell_solver->field_h_at_et);
          }

          // if (i < __geometry->domains_amount[0] - 1 && j < __geometry->domains_amount[1] - 1)
          // {
          //   Domain *dst_domain = __domains(i + 1, j + 1);
          //   sim_domain->field_h->field.overlay_xy(dst_domain->field_h->field);
          //   sim_domain->field_h->field_at_et.overlay_xy(dst_domain->field_h->field_at_et);
          // }
        }
    }

#ifdef ENABLE_MPI
  ////
  //// send particles to other domain
  ////
  vector<MPICommOverlay> domains_plus;
  vector<MPICommOverlay> domains_minus;

  vector<MPICommOverlay> domains_plus_recv;
  vector<MPICommOverlay> domains_minus_recv;

  for (unsigned int i = 0; i < r_domains; ++i)
  {
    // send to +1 rank
    Domain *sim_domain_plus = domains(i, z_domains - 1);

    int domain_r_size = sim_domain_plus->maxwell_solver->field_h_at_et[0].x_size;
    int domain_r_real_size = sim_domain_plus->maxwell_solver->field_h_at_et[0].x_real_size;
    int domain_o_s = sim_domain_plus->maxwell_solver->field_h_at_et[0].o_s;
    int domain_z_real_size = sim_domain_plus->maxwell_solver->field_h_at_et[0].y_real_size;

    MPICommOverlay item_plus (domain_r_real_size, i, 0);
    MPICommOverlay item_plus_recv (domain_r_real_size, i, 0);

    for (int j = -domain_o_s; j < domain_r_size + domain_o_s; ++j)
    {
      int j_s = j + domain_o_s;
      int d_r_s = domain_z_real_size - domain_o_s * 2;

      item_plus.col_0_0[j_s] = sim_domain_plus->maxwell_solver->field_h_at_et[0](j, d_r_s  - 4);
      item_plus.col_0_1[j_s] = sim_domain_plus->maxwell_solver->field_h_at_et[1](j, d_r_s  - 4);
      item_plus.col_0_2[j_s] = sim_domain_plus->maxwell_solver->field_h_at_et[2](j, d_r_s  - 4);

      item_plus.col_1_0[j_s] = sim_domain_plus->maxwell_solver->field_h_at_et[0](j, d_r_s  - 3);
      item_plus.col_1_1[j_s] = sim_domain_plus->maxwell_solver->field_h_at_et[1](j, d_r_s  - 3);
      item_plus.col_1_2[j_s] = sim_domain_plus->maxwell_solver->field_h_at_et[2](j, d_r_s  - 3);

      item_plus.col_2_0[j_s] = sim_domain_plus->maxwell_solver->field_h_at_et[0](j, d_r_s - 2);
      item_plus.col_2_1[j_s] = sim_domain_plus->maxwell_solver->field_h_at_et[1](j, d_r_s - 2);
      item_plus.col_2_2[j_s] = sim_domain_plus->maxwell_solver->field_h_at_et[2](j, d_r_s - 2);

      item_plus.col_3_0[j_s] = sim_domain_plus->maxwell_solver->field_h_at_et[0](j, d_r_s - 1);
      item_plus.col_3_1[j_s] = sim_domain_plus->maxwell_solver->field_h_at_et[1](j, d_r_s - 1);
      item_plus.col_3_2[j_s] = sim_domain_plus->maxwell_solver->field_h_at_et[2](j, d_r_s - 1);

    }

    domains_plus.push_back( item_plus );
    domains_plus_recv.push_back( item_plus_recv );

    // send to -1 rank
    Domain *sim_domain_minus = domains(i, 0);

    domain_r_size = sim_domain_minus->maxwell_solver->field_h_at_et[0].x_size;
    domain_o_s = sim_domain_minus->maxwell_solver->field_h_at_et[0].o_s;

    MPICommOverlay item_minus (domain_r_size + 2 * domain_o_s, i, 0);
    MPICommOverlay item_minus_recv (domain_r_size + 2 * domain_o_s, i, 0);

    for (int j = -domain_o_s; j < domain_r_size + domain_o_s; ++j)
    {
      int j_s = j + domain_o_s;

      item_minus.col_0_0[j_s] = sim_domain_minus->maxwell_solver->field_h_at_et[0](j, 0 - domain_o_s);
      item_minus.col_0_1[j_s] = sim_domain_minus->maxwell_solver->field_h_at_et[1](j, 0 - domain_o_s);
      item_minus.col_0_2[j_s] = sim_domain_minus->maxwell_solver->field_h_at_et[2](j, 0 - domain_o_s);

      item_minus.col_1_0[j_s] = sim_domain_minus->maxwell_solver->field_h_at_et[0](j, 1 - domain_o_s);
      item_minus.col_1_1[j_s] = sim_domain_minus->maxwell_solver->field_h_at_et[1](j, 1 - domain_o_s);
      item_minus.col_1_2[j_s] = sim_domain_minus->maxwell_solver->field_h_at_et[2](j, 1 - domain_o_s);

      item_minus.col_2_0[j_s] = sim_domain_minus->maxwell_solver->field_h_at_et[0](j, 2 - domain_o_s);
      item_minus.col_2_1[j_s] = sim_domain_minus->maxwell_solver->field_h_at_et[1](j, 2 - domain_o_s);
      item_minus.col_2_2[j_s] = sim_domain_minus->maxwell_solver->field_h_at_et[2](j, 2 - domain_o_s);

      item_minus.col_3_0[j_s] = sim_domain_minus->maxwell_solver->field_h_at_et[0](j, 3 - domain_o_s);
      item_minus.col_3_1[j_s] = sim_domain_minus->maxwell_solver->field_h_at_et[1](j, 3 - domain_o_s);
      item_minus.col_3_2[j_s] = sim_domain_minus->maxwell_solver->field_h_at_et[2](j, 3 - domain_o_s);
    }

    domains_minus.push_back( item_minus );
    domains_minus_recv.push_back( item_minus_recv );
  }

  if (world_rank < world_size - 1)
    for (auto i = domains_plus.begin(); i != domains_plus.end(); ++i)
    {
      // LOG_S(WARNING) << "1 Sending from " << world_rank << " to " << world_rank + 1;
      MPI_Send(MPI_BOTTOM, 1, i->mpi_dtype(), world_rank + 1, 0, MPI_COMM_WORLD);
    }

  if (world_rank > 0)
     for (auto i = domains_minus_recv.begin(); i != domains_minus_recv.end(); ++i)
     {
       // LOG_S(WARNING) << "2 Receiving from " << world_rank - 1 << " at " << world_rank;
       MPI_Recv(MPI_BOTTOM, 1, i->mpi_dtype(), world_rank - 1, 0, MPI_COMM_WORLD, NULL);
     }

  if (world_rank > 0)
    for (auto i = domains_minus.begin(); i != domains_minus.end(); ++i)
    {
      // LOG_S(WARNING) << "3 Sending from " << world_rank << " to " << world_rank - 1;
      MPI_Send(MPI_BOTTOM, 1, i->mpi_dtype(), world_rank - 1, 0, MPI_COMM_WORLD);
    }

  if (world_rank < world_size - 1)
    for (auto i = domains_plus_recv.begin(); i != domains_plus_recv.end(); ++i)
    {
      // LOG_S(WARNING) << "4 Receiving from " << world_rank + 1 << " at " << world_rank;
      MPI_Recv(MPI_BOTTOM, 1, i->mpi_dtype(), world_rank + 1, 0, MPI_COMM_WORLD, NULL);
    }

  MPI_Barrier(MPI_COMM_WORLD);

  ////////////////////////////////////////////////////////////

  for (auto i = domains_minus_recv.begin(); i != domains_minus_recv.end(); ++i)
  {
    Domain *dst_domain_minus = domains(i->dst_domain[0], 0);

    int domain_r_size = dst_domain_minus->maxwell_solver->field_h_at_et[0].x_size;
    int domain_o_s = dst_domain_minus->maxwell_solver->field_h_at_et[0].o_s;
    int domain_z_size = dst_domain_minus->maxwell_solver->field_h_at_et[0].y_size;

    for (int j = -domain_o_s; j < domain_r_size + domain_o_s; ++j)
    {
      int j_s = j + domain_o_s;

      dst_domain_minus->maxwell_solver->field_h_at_et[0].inc(j, 0 - domain_o_s, i->col_0_0[j_s]);
      dst_domain_minus->maxwell_solver->field_h_at_et[1].inc(j, 0 - domain_o_s, i->col_0_1[j_s]);
      dst_domain_minus->maxwell_solver->field_h_at_et[2].inc(j, 0 - domain_o_s, i->col_0_2[j_s]);

      dst_domain_minus->maxwell_solver->field_h_at_et[0].inc(j, 1 - domain_o_s, i->col_1_0[j_s]);
      dst_domain_minus->maxwell_solver->field_h_at_et[1].inc(j, 1 - domain_o_s, i->col_1_1[j_s]);
      dst_domain_minus->maxwell_solver->field_h_at_et[2].inc(j, 1 - domain_o_s, i->col_1_2[j_s]);

      dst_domain_minus->maxwell_solver->field_h_at_et[0].inc(j, 2 - domain_o_s, i->col_2_0[j_s]);
      dst_domain_minus->maxwell_solver->field_h_at_et[1].inc(j, 2 - domain_o_s, i->col_2_1[j_s]);
      dst_domain_minus->maxwell_solver->field_h_at_et[2].inc(j, 2 - domain_o_s, i->col_2_2[j_s]);

      dst_domain_minus->maxwell_solver->field_h_at_et[0].inc(j, 3 - domain_o_s, i->col_3_0[j_s]);
      dst_domain_minus->maxwell_solver->field_h_at_et[1].inc(j, 3 - domain_o_s, i->col_3_1[j_s]);
      dst_domain_minus->maxwell_solver->field_h_at_et[2].inc(j, 3 - domain_o_s, i->col_3_2[j_s]);
    }
  }

  for (auto i = domains_plus_recv.begin(); i != domains_plus_recv.end(); ++i)
  {
    Domain *dst_domain_plus = domains(i->dst_domain[0], z_domains - 1);

    int domain_r_size = dst_domain_plus->maxwell_solver->field_h_at_et[0].x_size;
    int domain_o_s = dst_domain_plus->maxwell_solver->field_h_at_et[0].o_s;
    int domain_z_real_size = dst_domain_plus->maxwell_solver->field_h_at_et[0].y_real_size;

    for (unsigned int j = -domain_o_s; j < domain_r_size + domain_o_s; ++j)
    {
      int j_s = j + domain_o_s;
      int d_r_s = domain_z_real_size - domain_o_s * 2;

      dst_domain_plus->maxwell_solver->field_h_at_et[0].inc(j, d_r_s - 1, i->col_3_0[j_s]);
      dst_domain_plus->maxwell_solver->field_h_at_et[1].inc(j, d_r_s - 1, i->col_3_1[j_s]);
      dst_domain_plus->maxwell_solver->field_h_at_et[2].inc(j, d_r_s - 1, i->col_3_2[j_s]);

      dst_domain_plus->maxwell_solver->field_h_at_et[0].inc(j, d_r_s - 2, i->col_2_0[j_s]);
      dst_domain_plus->maxwell_solver->field_h_at_et[1].inc(j, d_r_s - 2, i->col_2_1[j_s]);
      dst_domain_plus->maxwell_solver->field_h_at_et[2].inc(j, d_r_s - 2, i->col_2_2[j_s]);

      dst_domain_plus->maxwell_solver->field_h_at_et[0].inc(j, d_r_s - 3, i->col_1_0[j_s]);
      dst_domain_plus->maxwell_solver->field_h_at_et[1].inc(j, d_r_s - 3, i->col_1_1[j_s]);
      dst_domain_plus->maxwell_solver->field_h_at_et[2].inc(j, d_r_s - 3, i->col_1_2[j_s]);

      dst_domain_plus->maxwell_solver->field_h_at_et[0].inc(j, d_r_s - 4, i->col_0_0[j_s]);
      dst_domain_plus->maxwell_solver->field_h_at_et[1].inc(j, d_r_s - 4, i->col_0_1[j_s]);
      dst_domain_plus->maxwell_solver->field_h_at_et[2].inc(j, d_r_s - 4, i->col_0_2[j_s]);
    }
  }
#endif // ENABLE_MPI
}

void SMB::field_e_overlay ()
{
  for (unsigned int idx = 0; idx < 2; ++idx)
    for (unsigned int idy = 0; idy < 2; ++idy)
    {
#pragma omp parallel for
      for (unsigned int i = idx; i < r_domains; i+=2)
        for (unsigned int j = idy; j < z_domains; j+=2)
        {
          Domain *sim_domain = domains(i, j);

          // update grid
          if (i < geometry->domains_amount[0] - 1)
          {
            Domain *dst_domain = domains(i+1, j);
            sim_domain->maxwell_solver->field_e.overlay_x(
              dst_domain->maxwell_solver->field_e
              );
          }

          if (j < geometry->domains_amount[1] - 1)
          {
            Domain *dst_domain = domains(i, j + 1);
            sim_domain->maxwell_solver->field_e.overlay_y(
              dst_domain->maxwell_solver->field_e
              );
          }

          // if (i < __geometry->domains_amount[0] - 1 && j < __geometry->domains_amount[1] - 1)
          // {
          //   Domain *dst_domain = __domains(i + 1, j + 1);
          //   sim_domain->field_e->field.overlay_xy(dst_domain->field_e->field);
          // }
        }
    }

#ifdef ENABLE_MPI
  ////
  //// send particles to other domain
  ////
  vector<MPICommOverlay> domains_plus;
  vector<MPICommOverlay> domains_minus;

  vector<MPICommOverlay> domains_plus_recv;
  vector<MPICommOverlay> domains_minus_recv;

  for (unsigned int i = 0; i < r_domains; ++i)
  {
    // send to +1 rank
    Domain *sim_domain_plus = domains(i, z_domains - 1);

    int domain_r_size = sim_domain_plus->maxwell_solver->field_e[0].x_size;
    int domain_r_real_size = sim_domain_plus->maxwell_solver->field_e[0].x_real_size;
    int domain_o_s = sim_domain_plus->maxwell_solver->field_e[0].o_s;
    int domain_z_real_size = sim_domain_plus->maxwell_solver->field_e[0].y_real_size;

    MPICommOverlay item_plus (domain_r_real_size, i, 0);
    MPICommOverlay item_plus_recv (domain_r_real_size, i, 0);

    for (int j = -domain_o_s; j < domain_r_size + domain_o_s; ++j)
    {
      int j_s = j + domain_o_s;
      int d_r_s = domain_z_real_size - domain_o_s * 2;

      item_plus.col_0_0[j_s] = sim_domain_plus->maxwell_solver->field_e[0](j, d_r_s  - 4);
      item_plus.col_0_1[j_s] = sim_domain_plus->maxwell_solver->field_e[1](j, d_r_s  - 4);
      item_plus.col_0_2[j_s] = sim_domain_plus->maxwell_solver->field_e[2](j, d_r_s  - 4);

      item_plus.col_1_0[j_s] = sim_domain_plus->maxwell_solver->field_e[0](j, d_r_s  - 3);
      item_plus.col_1_1[j_s] = sim_domain_plus->maxwell_solver->field_e[1](j, d_r_s  - 3);
      item_plus.col_1_2[j_s] = sim_domain_plus->maxwell_solver->field_e[2](j, d_r_s  - 3);

      item_plus.col_2_0[j_s] = sim_domain_plus->maxwell_solver->field_e[0](j, d_r_s - 2);
      item_plus.col_2_1[j_s] = sim_domain_plus->maxwell_solver->field_e[1](j, d_r_s - 2);
      item_plus.col_2_2[j_s] = sim_domain_plus->maxwell_solver->field_e[2](j, d_r_s - 2);

      item_plus.col_3_0[j_s] = sim_domain_plus->maxwell_solver->field_e[0](j, d_r_s - 1);
      item_plus.col_3_1[j_s] = sim_domain_plus->maxwell_solver->field_e[1](j, d_r_s - 1);
      item_plus.col_3_2[j_s] = sim_domain_plus->maxwell_solver->field_e[2](j, d_r_s - 1);

    }

    domains_plus.push_back( item_plus );
    domains_plus_recv.push_back( item_plus_recv );

    // send to -1 rank
    Domain *sim_domain_minus = domains(i, 0);

    domain_r_size = sim_domain_minus->maxwell_solver->field_e[0].x_size;
    domain_o_s = sim_domain_minus->maxwell_solver->field_e[0].o_s;

    MPICommOverlay item_minus (domain_r_size + 2 * domain_o_s, i, 0);
    MPICommOverlay item_minus_recv (domain_r_size + 2 * domain_o_s, i, 0);

    for (int j = -domain_o_s; j < domain_r_size + domain_o_s; ++j)
    {
      int j_s = j + domain_o_s;

      item_minus.col_0_0[j_s] = sim_domain_minus->maxwell_solver->field_e[0](j, 0 - domain_o_s);
      item_minus.col_0_1[j_s] = sim_domain_minus->maxwell_solver->field_e[1](j, 0 - domain_o_s);
      item_minus.col_0_2[j_s] = sim_domain_minus->maxwell_solver->field_e[2](j, 0 - domain_o_s);

      item_minus.col_1_0[j_s] = sim_domain_minus->maxwell_solver->field_e[0](j, 1 - domain_o_s);
      item_minus.col_1_1[j_s] = sim_domain_minus->maxwell_solver->field_e[1](j, 1 - domain_o_s);
      item_minus.col_1_2[j_s] = sim_domain_minus->maxwell_solver->field_e[2](j, 1 - domain_o_s);

      item_minus.col_2_0[j_s] = sim_domain_minus->maxwell_solver->field_e[0](j, 2 - domain_o_s);
      item_minus.col_2_1[j_s] = sim_domain_minus->maxwell_solver->field_e[1](j, 2 - domain_o_s);
      item_minus.col_2_2[j_s] = sim_domain_minus->maxwell_solver->field_e[2](j, 2 - domain_o_s);

      item_minus.col_3_0[j_s] = sim_domain_minus->maxwell_solver->field_e[0](j, 3 - domain_o_s);
      item_minus.col_3_1[j_s] = sim_domain_minus->maxwell_solver->field_e[1](j, 3 - domain_o_s);
      item_minus.col_3_2[j_s] = sim_domain_minus->maxwell_solver->field_e[2](j, 3 - domain_o_s);
    }

    domains_minus.push_back( item_minus );
    domains_minus_recv.push_back( item_minus_recv );
  }

  if (world_rank < world_size - 1)
    for (auto i = domains_plus.begin(); i != domains_plus.end(); ++i)
    {
      // LOG_S(WARNING) << "1 Sending from " << world_rank << " to " << world_rank + 1;
      MPI_Send(MPI_BOTTOM, 1, i->mpi_dtype(), world_rank + 1, 0, MPI_COMM_WORLD);
    }

  if (world_rank > 0)
     for (auto i = domains_minus_recv.begin(); i != domains_minus_recv.end(); ++i)
     {
       // LOG_S(WARNING) << "2 Receiving from " << world_rank - 1 << " at " << world_rank;
       MPI_Recv(MPI_BOTTOM, 1, i->mpi_dtype(), world_rank - 1, 0, MPI_COMM_WORLD, NULL);
     }

  if (world_rank > 0)
    for (auto i = domains_minus.begin(); i != domains_minus.end(); ++i)
    {
      // LOG_S(WARNING) << "3 Sending from " << world_rank << " to " << world_rank - 1;
      MPI_Send(MPI_BOTTOM, 1, i->mpi_dtype(), world_rank - 1, 0, MPI_COMM_WORLD);
    }

  if (world_rank < world_size - 1)
    for (auto i = domains_plus_recv.begin(); i != domains_plus_recv.end(); ++i)
    {
      // LOG_S(WARNING) << "4 Receiving from " << world_rank + 1 << " at " << world_rank;
      MPI_Recv(MPI_BOTTOM, 1, i->mpi_dtype(), world_rank + 1, 0, MPI_COMM_WORLD, NULL);
    }

  MPI_Barrier(MPI_COMM_WORLD);

  ////////////////////////////////////////////////////////////

  for (auto i = domains_minus_recv.begin(); i != domains_minus_recv.end(); ++i)
  {
    Domain *dst_domain_minus = domains(i->dst_domain[0], 0);

    int domain_r_size = dst_domain_minus->maxwell_solver->field_e[0].x_size;
    int domain_o_s = dst_domain_minus->maxwell_solver->field_e[0].o_s;
    int domain_z_size = dst_domain_minus->maxwell_solver->field_e[0].y_size;

    for (int j = -domain_o_s; j < domain_r_size + domain_o_s; ++j)
    {
      int j_s = j + domain_o_s;

      dst_domain_minus->maxwell_solver->field_e[0].inc(j, 0 - domain_o_s, i->col_0_0[j_s]);
      dst_domain_minus->maxwell_solver->field_e[1].inc(j, 0 - domain_o_s, i->col_0_1[j_s]);
      dst_domain_minus->maxwell_solver->field_e[2].inc(j, 0 - domain_o_s, i->col_0_2[j_s]);

      dst_domain_minus->maxwell_solver->field_e[0].inc(j, 1 - domain_o_s, i->col_1_0[j_s]);
      dst_domain_minus->maxwell_solver->field_e[1].inc(j, 1 - domain_o_s, i->col_1_1[j_s]);
      dst_domain_minus->maxwell_solver->field_e[2].inc(j, 1 - domain_o_s, i->col_1_2[j_s]);

      dst_domain_minus->maxwell_solver->field_e[0].inc(j, 2 - domain_o_s, i->col_2_0[j_s]);
      dst_domain_minus->maxwell_solver->field_e[1].inc(j, 2 - domain_o_s, i->col_2_1[j_s]);
      dst_domain_minus->maxwell_solver->field_e[2].inc(j, 2 - domain_o_s, i->col_2_2[j_s]);

      dst_domain_minus->maxwell_solver->field_e[0].inc(j, 3 - domain_o_s, i->col_3_0[j_s]);
      dst_domain_minus->maxwell_solver->field_e[1].inc(j, 3 - domain_o_s, i->col_3_1[j_s]);
      dst_domain_minus->maxwell_solver->field_e[2].inc(j, 3 - domain_o_s, i->col_3_2[j_s]);
    }
  }

  for (auto i = domains_plus_recv.begin(); i != domains_plus_recv.end(); ++i)
  {
    Domain *dst_domain_plus = domains(i->dst_domain[0], z_domains - 1);

    int domain_r_size = dst_domain_plus->maxwell_solver->field_e[0].x_size;
    int domain_o_s = dst_domain_plus->maxwell_solver->field_e[0].o_s;
    int domain_z_real_size = dst_domain_plus->maxwell_solver->field_e[0].y_real_size;

    for (unsigned int j = -domain_o_s; j < domain_r_size + domain_o_s; ++j)
    {
      int j_s = j + domain_o_s;
      int d_r_s = domain_z_real_size - domain_o_s * 2;

      dst_domain_plus->maxwell_solver->field_e[0].inc(j, d_r_s - 1, i->col_3_0[j_s]);
      dst_domain_plus->maxwell_solver->field_e[1].inc(j, d_r_s - 1, i->col_3_1[j_s]);
      dst_domain_plus->maxwell_solver->field_e[2].inc(j, d_r_s - 1, i->col_3_2[j_s]);

      dst_domain_plus->maxwell_solver->field_e[0].inc(j, d_r_s - 2, i->col_2_0[j_s]);
      dst_domain_plus->maxwell_solver->field_e[1].inc(j, d_r_s - 2, i->col_2_1[j_s]);
      dst_domain_plus->maxwell_solver->field_e[2].inc(j, d_r_s - 2, i->col_2_2[j_s]);

      dst_domain_plus->maxwell_solver->field_e[0].inc(j, d_r_s - 3, i->col_1_0[j_s]);
      dst_domain_plus->maxwell_solver->field_e[1].inc(j, d_r_s - 3, i->col_1_1[j_s]);
      dst_domain_plus->maxwell_solver->field_e[2].inc(j, d_r_s - 3, i->col_1_2[j_s]);

      dst_domain_plus->maxwell_solver->field_e[0].inc(j, d_r_s - 4, i->col_0_0[j_s]);
      dst_domain_plus->maxwell_solver->field_e[1].inc(j, d_r_s - 4, i->col_0_1[j_s]);
      dst_domain_plus->maxwell_solver->field_e[2].inc(j, d_r_s - 4, i->col_0_2[j_s]);
    }
  }
#endif // ENABLE_MPI
}


void SMB::solve_maxvell()
{
#pragma omp parallel for collapse(2)
  for (unsigned int i=0; i < r_domains; i++)
    for (unsigned int j = 0; j < z_domains; j++)
    {
      Domain *sim_domain = domains(i, j);

      sim_domain->weight_field_e();
    }
  field_e_overlay();

#pragma omp parallel for collapse(2)
  for (unsigned int i=0; i < r_domains; i++)
    for (unsigned int j = 0; j < z_domains; j++)
    {
      Domain *sim_domain = domains(i, j);

      // ! 2. Calculate magnetic field (H)
      sim_domain->weight_field_h(); // +
    }
  field_h_overlay();
}

void SMB::solve_current()
{
#pragma omp parallel for collapse(2)
  for (unsigned int i=0; i < r_domains; i++)
    for (unsigned int j = 0; j < z_domains; j++)
    {
      Domain *sim_domain = domains(i, j);

      sim_domain->reset_current();
      sim_domain->weight_current();
    }
  current_overlay();
}

void SMB::advance_particles()
{
#pragma omp parallel for collapse(2)
  for (unsigned int i=0; i < r_domains; i++)
    for (unsigned int j = 0; j < z_domains; j++)
    {
      Domain *sim_domain = domains(i, j);

      // ! 3. Calculate velocity
      sim_domain->push_particles();

#ifdef ENABLE_COULOMB_COLLISIONS
      sim_domain->collide(); // collide before reflect
#endif // ENABLE_COULOMB_COLLISIONS

      sim_domain->dump_particle_positions_to_old();
      sim_domain->update_particles_coords();
      sim_domain->particles_back_position_to_rz();

      sim_domain->reflect();

      sim_domain->particles_back_velocity_to_rz();

      sim_domain->bind_cell_numbers();
    }
  particles_runaway_collector();
}

void SMB::inject_beam()
{
  if (! cfg->particle_beams.empty())
  {
#pragma omp parallel for collapse(2)
    for (unsigned int i = 0; i < r_domains; i++)
      for (unsigned int j = 0; j < z_domains; j++)
      {
        Domain *sim_domain = domains(i, j);

        // ! 1. manage beam
        sim_domain->manage_beam();
      }
  }
}

void SMB::distribute()
{
#pragma omp parallel for collapse(2)
  for (unsigned int i=0; i < r_domains; i++)
    for (unsigned int j = 0; j < z_domains; j++)
    {
      Domain *sim_domain = domains(i, j);

      sim_domain->distribute(); // spatial and velocity distribution
    }
}

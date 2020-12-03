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

// SMB constructor
SMB::SMB ( Cfg* _cfg, Geometry *_geometry, TimeSim *_time )
  : cfg(_cfg), time(_time), geometry(_geometry)
{
  LOG_S(MAX) << "Initializing Geometry, Particle Species and Simulation Domains";

  //
  // initialize domains grid
  //
  r_domains = geometry->domains_amount[0];
  z_domains = geometry->domains_amount[1];

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
        wall_r0 = true;
      if (i == r_domains - 1)
        wall_rr = true;
      if (j == 0)
        wall_z0 = true;
      if (j == z_domains - 1)
        wall_zz = true;

#ifdef ENABLE_PML
      // set PML to domains
      if (geometry->size[0] - geometry->size[0] / r_domains * (i + 1)
          < geometry->pml_length[2])
        pml_l_rwall = geometry->pml_length[2];
      if (j * geometry->size[1] / z_domains < geometry->pml_length[1])
        pml_l_z0 = geometry->pml_length[1];
      if (geometry->size[1] - geometry->size[1] / z_domains * (j + 1)
          < geometry->pml_length[3])
        pml_l_zwall = geometry->pml_length[3];
#endif // ENABLE_PML

      unsigned int bot_r = (unsigned int)geometry->cell_amount[0] * i / r_domains;
      unsigned int top_r = (unsigned int)geometry->cell_amount[0] * (i + 1) / r_domains;

      unsigned int left_z = (unsigned int)geometry->cell_amount[1] * j / z_domains;
      unsigned int right_z = (unsigned int)geometry->cell_amount[1] * (j + 1) / z_domains;

#ifdef ENABLE_PML
      Geometry *geom_domain = new Geometry (
        { geometry->size[0] / r_domains, geometry->size[1] / z_domains},
        { bot_r, left_z, top_r, right_z },
        { pml_l_z0 * z_domains,    // multiplying is a workaround, because domain
          pml_l_zwall * z_domains, // doesn't know abount whole simulation domain size
          pml_l_rwall * r_domains }, // aka geometry
        geometry->pml_sigma,
        { wall_r0,
          wall_z0,
          wall_rr,
          wall_zz }
        );

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

        SpecieP *pps = new SpecieP (p_id_counter,
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
            BeamP *beam = new BeamP (b_id_counter, ((string)"beam_").append(bm->name),
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
      // LOG_S(FATAL) << domains(1,1);
    };





////////////////////////////////////////////////////////////////////////////////
  unsigned int r_domains = domains.x_size;
  unsigned int z_domains = domains.y_size;

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
            (**ps).particles.erase(
              std::remove_if(
                (**ps).particles.begin(), (**ps).particles.end(),
                [ &j_c, &r_c, &ps, &__domains, &sim_domain, &i, &j,
                  &__geometry ] ( vector <double> * & o )
                {
                  bool res = false;

                  // not unsigned, because it could be less, than zero
                  int r_cell = P_CELL_R((*o));
                  int z_cell = P_CELL_Z((*o));

                  unsigned int i_dst = (unsigned int)ceil(r_cell / sim_domain->geometry.cell_amount[0]);
                  unsigned int j_dst = (unsigned int)ceil(z_cell / sim_domain->geometry.cell_amount[1]);

                  if (r_cell < 0 || z_cell < 0)
                  {
                    LOG_S(ERROR) << "Particle's position is less, than 0. Position is: ["
                                 << P_VEL_R((*o)) << ", "
                                 << P_VEL_Z((*o)) << "]. Removing";
                    ++r_c;
                    res = true;
                  }

                  else if (r_cell >= __geometry->cell_amount[0])
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
                                   << __geometry->cell_amount[0]
                                   << ". Position is: ["
                                   << P_POS_R((*o)) << ", "
                                   << P_POS_Z((*o)) << "]. Removing";
                    }
                    ++r_c;
                    res = true;
                  }

                  // remove out-of-simulation particles
                  else if (z_cell >= __geometry->cell_amount[1])
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
                                   << __geometry->cell_amount[1]
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

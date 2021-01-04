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

#include "beamP.hpp"

using namespace constant;

BeamP::BeamP (unsigned int id, // ID for every particles specie
              string p_name,
              double p_charge, // in electron charges
              double p_mass, // in electron masses
              unsigned int p_macro_amount,
              double p_start_time,
              double b_radius,
              double b_density, // CI units
              double b_amount,
              double b_length,
              double b_distance,
              double b_velocity, // electronvolts
              Geometry *geom,
              TimeSim *t )
  : SpecieP (id, p_name, p_charge, p_mass,
             p_macro_amount, b_density,
             b_density, 0, geom, t)
{
  current_bunch_number = 0;

  radius = b_radius;
  density =  b_density;
  velocity = b_velocity;
  start_time = p_start_time;
  bunches_amount = b_amount;
  bunch_length = b_length;
  bunches_distance = b_distance;
  finish_time = start_time
    + bunch_length / velocity * bunches_amount
    + bunches_distance / velocity * (bunches_amount - 1);


  // find bunch radius, compared to current domain
  if (geometry->cell_dims[1] > 0)
    domain_radius = 0;
  else if (radius > geometry->cell_dims[2] * geometry->cell_size[0])
    domain_radius = geometry->size[0];
  else if (radius > geometry->cell_dims[0] * geometry->cell_size[0])
    domain_radius = radius - geometry->cell_dims[0] * geometry->cell_size[0];
  else
    domain_radius = 0;

  double macro_normalizer = domain_radius / radius;

  bunch_macro_amount = macro_amount / bunches_amount * macro_normalizer;

  macro_per_step_to_inject = bunch_macro_amount * velocity * time->step / bunch_length;
}

void BeamP::inject()
{
  // check if limit of bunches is reached
  if (current_bunch_number < bunches_amount)
  {
    double bunch_start_time = start_time + (bunch_length + bunches_distance) / velocity * current_bunch_number;
    double bunch_finish_time = bunch_start_time + bunch_length / velocity;

    if (geometry->cell_dims[0] * geometry->cell_size[0] < radius // prevent injection for unused bunches
        && time->current >= bunch_start_time  // launch only
        && time->current < bunch_finish_time) // at injection time
    {
      double dr = geometry->cell_size[0];
      double dz = geometry->cell_size[1];
      double dl = velocity * time->step;
      double half_z_cell_size = geometry->cell_size[1] / 2.;

      // as volume of single macro is $2 \pi r dr dz$
      // than its proportional to $r$
      // so average value, respecting
      // to flat random particles distrribution,
      // should be just $2 \pi \frac{radius}{2} dr dz$
      // average macro volume is max_volume/2
      // so radius to get such volume is r_max/sqrt(2)
      // WARNING: used only for homogenous distribution
      double v_macro_avg = PI * radius * dr * dz;
      double N_total = PI * radius * radius * bunch_length * bunches_amount * density;
      double n_per_macro_avg = N_total / macro_amount;

      for (unsigned int i = 0; i < macro_per_step_to_inject; ++i)
      {
        // vector<double> *v = new vector<double>(P_VEC_SIZE, 0);
        Particle *v = new Particle();

        double rand_r = math::random::uniform();
        double rand_z = math::random::uniform();

        double pos_r = domain_radius * rand_r + dr / 2 + geometry->cell_dims[0] * geometry->cell_size[0];
        // 1. set pos_r
        P_POS_R((*v)) = pos_r;
        // 2. set pos_phi
        P_POS_PHI((*v)) = 0;
        // 3. set pos_z
        P_POS_Z((*v)) = dl * rand_z + half_z_cell_size;
        // 4. set vel_r
        P_VEL_R((*v)) = 0;
        // 5. set vel_phi
        P_VEL_PHI((*v)) = 0;
        // 6. set vel_z
        P_VEL_Z((*v)) = velocity;

        // coefitient of normalization
        double norm;
        if (pos_r > dr / 2)
          norm = 2 * PI * pos_r * dr * dz / v_macro_avg;
        else
          norm = PI * (pos_r * pos_r + pos_r * dr + dr * dr / 4.) * dz / v_macro_avg;

        // number of real particles per macroparticle
        double n_per_macro = n_per_macro_avg * norm;

        // 7. set charge
        P_WEIGHT((*v)) = n_per_macro;

        // set cell numbers
        int r_cell = CELL_NUMBER(P_POS_R((*v)), geometry->cell_size[0]);
        int z_cell = CELL_NUMBER(P_POS_R((*v)), geometry->cell_size[1]);

        P_CELL_R((*v)) = r_cell;
        P_CELL_Z((*v)) = z_cell;

        // push particle to particles beam vector
        particles.push_back(v);
        // particles_old.push_back(v_old);
      }
    }

    // increase current bunch number after whole bunch injection
    if (time->current > bunch_finish_time
        && geometry->cell_dims[1] == 0 // print message only for domain 0,0
        && geometry->cell_dims[0] == 0)
    {
      LOG_S(MAX) << "Bunch #" << current_bunch_number << " of beam ``" << name << "'' has been injected";
      ++current_bunch_number;
    }
  }
}

void BeamP::reflect ()
{
  double dr = geometry->cell_size[0];
  double half_dr = dr / 2.;

  // shift for converting local positions into global and back
  double r_shift = geometry->cell_dims[0] * dr;

  for (auto p = particles.begin(); p != particles.end(); ++p)
  {
    double pos_r = P_POS_R((**p)) - r_shift;

    if (pos_r < half_dr && geometry->walls[0])
    {
      P_POS_R((**p)) = dr - pos_r + r_shift;
      P_VEL_R((**p)) = - P_VEL_R((**p));
    }
  }
}

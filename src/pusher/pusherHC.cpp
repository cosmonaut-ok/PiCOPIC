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


#include "pusher/pusherHC.hpp"

void PusherHC::operator()()
{
  // !
  // ! Higuera-Cary pusher
  // !
  for (auto sp = species_p.begin(); sp != species_p.end(); ++sp)
  {
    double charge = (**sp).charge;
    double mass = (**sp).mass;

    for (auto p = (**sp).particles.begin(); p != (**sp).particles.end(); ++p)
    {
      vector3d<double> velocity(P_VEL_R((**p)), P_VEL_PHI((**p)), P_VEL_Z((**p)));
      vector3d<double> uplocity; // u prime
      vector3d<double> psm; // pxsm, pysm, pzsm
      vector3d<double> um; // pxsm, pysm, pzsm
      vector3d<double> up; // pxsm, pysm, pzsm

      // we don't care, if it is particle's or macroparticle's
      // charge over mass ratio
      double charge_over_2mass_dt = charge * time->step / (2 * mass);

      double pos_r = P_POS_R((**p));
      double pos_z = P_POS_Z((**p));

      vector3d<double> e = maxwell_solver->get_field_e(pos_r, pos_z);
      vector3d<double> b = maxwell_solver->get_field_h(pos_r, pos_z);
      vector3d<double> b2;
      vector3d<double> b_cross;

      double gamma, sq_vel, B2;

      // convert velocity to relativistic momentum
      sq_vel = velocity.length2();
      gamma = phys::rel::lorenz_factor(sq_vel);
      velocity *= gamma;

      //// enter main algo
      // init Half-acceleration in the electric field
      e *= charge_over_2mass_dt;
      psm = e;

      um = velocity;
      um += psm;
      // Intermediate gamma factor: only this part differs from the Boris scheme
      // Square Gamma factor from um
      double gfm2 = (1. + um.length2());

      b *= charge_over_2mass_dt * MAGN_CONST;
      B2 = b.length2();

      // Equivalent of 1/\gamma_{new} in the paper
      gamma = 1. / sqrt ( 0.5*( gfm2 - B2 +
                                sqrt ( pow( gfm2 - B2, 2 )
                                       + 4.0 * ( B2 + pow( b[0]*um[0]
                                                           + b[1]*um[1]
                                                           + b[2]*um[2], 2 ) ) ) ) );

      b *= gamma;
      b2 = b;
      b2 *= b;

      b_cross[0] = b[0]*b[1];
      b_cross[1] = b[1]*b[2];
      b_cross[2] = b[2]*b[0];
      double inv_det_B = 1.0/( 1.0+b2[0]+b2[1]+b2[0] );

      up[0] = ( ( 1.0+b2[0]-b2[1]-b2[2] ) * um[0] + 2. * ( b_cross[0]+b[2] )
                * um[1] + 2. * ( b_cross[2] - b[2] ) * um[2] ) * inv_det_B;
      up[1] = ( 2. * ( b_cross[0]-b[2] ) * um[0] + ( 1. - b2[0]+b2[1]-b2[2] )
                * um[1] + 2. * ( b_cross[1] + b[0] ) * um[2] ) * inv_det_B;
      up[2] = ( 2. * ( b_cross[2] + b[1] ) * um[0] + 2. * ( b_cross[1] - b[0] )
                * um[1] + ( 1. - b2[0]-b2[1]+b2[2] )* um[2] ) * inv_det_B;

      // finalize Half-acceleration in the electric field
      psm += up;

      //// exit main algo

      sq_vel = psm.length2();
      gamma = phys::rel::lorenz_factor_inv(sq_vel);
      psm *= gamma;

      P_VEL_R((**p)) = psm[0];
      P_VEL_PHI((**p)) = psm[1];
      P_VEL_Z((**p)) = psm[2];
    }
  }
}

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

#include "pusher/pusherVay.hpp"

using namespace constant;

void PusherVay::operator()()
{
  // !
  // ! Vay pusher
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

      // we don't care, if it is particle's or macroparticle's
      // charge over mass ratio
      double charge_over_2mass_dt = charge * time->step / (2 * mass);

      double pos_r = P_POS_R((**p));
      double pos_z = P_POS_Z((**p));

      vector3d<double> e = maxwell_solver->get_field_e(pos_r, pos_z);
      vector3d<double> b = maxwell_solver->get_field_h(pos_r, pos_z);

      double gamma, sq_vel, s, us2, alpha, B2;

      // convert velocity to relativistic momentum
      sq_vel = velocity.length2();
      gamma = phys::rel::lorenz_factor(sq_vel);
      velocity *= gamma;

      //
      // Part I: Computation of uprime
      //

      // Add Electric field
      uplocity = velocity;
      e *= 2. * charge_over_2mass_dt;
      uplocity += e;

      // Add magnetic field
      b *= charge_over_2mass_dt * MAGN_CONST;

      // Smilei: For unknown reason, this has to be computed again
      sq_vel = velocity.length2();
      gamma = phys::rel::lorenz_factor_inv(sq_vel);

      uplocity[0] += gamma * ( velocity[1] * b[2] - velocity[2] * b[1] );
      uplocity[1] += gamma * ( velocity[2] * b[0] - velocity[0] * b[2] );
      uplocity[2] += gamma * ( velocity[0] * b[1] - velocity[1] * b[0] );

      // alpha is gamma^2
      alpha = 1. + uplocity.length2();
      B2 = b.length2();

      //
      // Part II: Computation of Gamma^{i+1}
      //

      // s is sigma
      s = alpha - B2;
      // TODO: implement * operator for vector3d
      us2 = pow(uplocity.dot(b), 2);

      // alpha becomes 1/gamma^{i+1}
      alpha = 1. / sqrt( 0.5 * ( s + sqrt( s * s + 4. * ( B2 + us2 ) ) ) );

      b *= alpha;

      s = 1. / ( 1. + b.length2() );
      alpha = uplocity.dot(b);

      psm[0] = s * ( uplocity[0] + alpha*b[0] + b[2]*uplocity[1] - b[1]*uplocity[2] );
      psm[1] = s * ( uplocity[1] + alpha*b[1] + b[0]*uplocity[2] - b[2]*uplocity[0] );
      psm[2] = s * ( uplocity[2] + alpha*b[2] + b[1]*uplocity[0] - b[0]*uplocity[1] );

      sq_vel = psm.length2();
      gamma = phys::rel::lorenz_factor_inv(sq_vel);
      psm *= gamma;

      P_VEL_R((**p)) = psm[0];
      P_VEL_PHI((**p)) = psm[1];
      P_VEL_Z((**p)) = psm[2];
    }
  }
}

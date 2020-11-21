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


#include "pusher/pusherBoris.hpp"

void PusherBoris::operator()()
{
  // !
  // ! boris pusher
  // !
  for (auto sp = species_p.begin(); sp != species_p.end(); ++sp)
  {
    double charge = (**sp).charge;
    double mass = (**sp).mass;

    for (auto p = (**sp).particles.begin(); p != (**sp).particles.end(); ++p)
    {
      // define vars directly in loop, because of multithreading
      double charge_over_2mass_dt, const2, sq_velocity;
#ifdef SWITCH_PUSHER_BORIS_ADAPTIVE
      bool use_rel; // use relativistic calculations
#endif // SWITCH_PUSHER
#if defined(SWITCH_PUSHER_BORIS_ADAPTIVE) || defined(SWITCH_PUSHER_BORIS_RELATIVISTIC) // do not use gamma in non-relativistic boris pusher
      double gamma = 1;
#endif // SWITCH_PUSHER

      vector3d<double> velocity(P_VEL_R((**p)), P_VEL_PHI((**p)), P_VEL_Z((**p)));
      vector3d<double> vtmp;

      double pos_r = P_POS_R((**p));
      double pos_z = P_POS_Z((**p));

      // check if radius and longitude are correct
      if (isnan(pos_r) ||
          isinf(pos_r) != 0 ||
          isnan(pos_z) ||
          isinf(pos_z) != 0)
        LOG_S(FATAL) << "(boris_pusher): radius[" << pos_r
                     << "] or longitude[" << pos_z
                     << "] is not valid number. Can not continue.";

      vector3d<double> e = maxwell_solver->get_field_e(pos_r, pos_z);
      vector3d<double> b = maxwell_solver->get_field_h(pos_r, pos_z);

      charge_over_2mass_dt = charge * time->step / (2 * mass); // we just shortened particle weight and use only q/m relation

      e *= charge_over_2mass_dt;
      b *= charge_over_2mass_dt * MAGN_CONST;

      // ! 0. check, if we should use classical calculations.
      // ! Required to increase modeling speed
#ifdef SWITCH_PUSHER_BORIS_ADAPTIVE
      if (pow(velocity[0], 2) + pow(velocity[1], 2) + pow(velocity[2], 2) > REL_LIMIT_POW_2)
        use_rel = true;
#endif // SWITCH_PUSHER
      // ! 1. Multiplication by relativistic factor (only for relativistic case)
      // ! \f$ u_{n-\frac{1}{2}} = \gamma_{n-\frac{1}{2}} * v_{n-\frac{1}{2}} \f$
#ifdef SWITCH_PUSHER_BORIS_ADAPTIVE
      if (use_rel)
#endif // SWITCH_PUSHER
#if defined(SWITCH_PUSHER_BORIS_RELATIVISTIC) || defined(SWITCH_PUSHER_BORIS_ADAPTIVE)
      {
        sq_velocity = velocity.length2();

        gamma = phys::rel::lorenz_factor(sq_velocity);
        velocity *= gamma;
      }
#endif // SWITCH_PUSHER

      // ! 2. Half acceleration in the electric field
      // ! \f$ u'_n = u_{n-\frac{1}{2}} + \frac{q dt}{2 m E(n)} \f$
      // ! \f$ u'_n = u_{n-1 / 2} + \frac{q dt}{2 m E(n)} \f$
      velocity += e;

      // ! 3. Rotation in the magnetic field
      // ! \f$ u" = u' + \frac{2}{1 + B'^2}  [(u' + [u' \times B'(n)] ) \times B'(n)] \f$,
      // ! \f$ B'(n) = \frac{B(n) q dt}{2 m * \gamma_n} \f$
#ifdef SWITCH_PUSHER_BORIS_ADAPTIVE
      if (use_rel)
#endif // SWITCH_PUSHER
#if defined(SWITCH_PUSHER_BORIS_RELATIVISTIC) || defined(SWITCH_PUSHER_BORIS_ADAPTIVE)
      {
        sq_velocity = velocity.length2();
        gamma = phys::rel::lorenz_factor_inv(sq_velocity);
        b *= gamma;
      }
#endif // SWITCH_PUSHER
      // ! \f$ const2 = \frac{2}{1 + b_1^2 + b_2^2 + b_3^2} \f$
      const2 = 2. / (1. + b.length2());

      // set temporary velocity as old values
      // to calculate magnetic rotation
      vtmp = velocity;

      velocity[0] = vtmp[0] + const2 * (
        (vtmp[1] - vtmp[0] * b[2] + vtmp[2] * b[0]) * b[2]
        - (vtmp[2] + vtmp[0] * b[1] - vtmp[1] * b[0]) * b[1]
        );
      velocity[1] = vtmp[1] + const2 * (
        -(vtmp[0] + vtmp[1] * b[2] - vtmp[2] * b[1]) * b[2]
        + (vtmp[2] + vtmp[0] * b[1] - vtmp[1] * b[0]) * b[0]
        );
      velocity[2] = vtmp[2] + const2 * (
        (vtmp[0] + vtmp[1] * b[2] - vtmp[2] * b[1]) * b[1]
        - (vtmp[1] - vtmp[0] * b[2] + vtmp[2] * b[0]) * b[0]
        );

      // ! 4. Half acceleration in the electric field
      // ! \f$ u_{n + \frac{1}{2}} = u_n + \frac{q dt}{2 m E(n)} \f$
      velocity += e;

      // ! 5. Division by relativistic factor
#ifdef SWITCH_PUSHER_BORIS_ADAPTIVE
      if (use_rel)
#endif // SWITCH_PUSHER
#if defined(SWITCH_PUSHER_BORIS_RELATIVISTIC) || defined(SWITCH_PUSHER_BORIS_ADAPTIVE)
      {
        sq_velocity = velocity.length2();
        gamma = phys::rel::lorenz_factor_inv(sq_velocity);
        velocity *= gamma;
      }
#endif // SWITCH_PUSHER

      P_VEL_R((**p)) = velocity[0];
      P_VEL_PHI((**p)) = velocity[1];
      P_VEL_Z((**p)) = velocity[2];
    }
  }
}

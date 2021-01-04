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

#ifndef _SPECIE_P_HPP_
#define _SPECIE_P_HPP_

#include <math.h>
#include <string>
#include <vector>

#include "defines.hpp"
#include "msg.hpp"

#include "math/vector3d.hpp"
#include "math/rand.hpp"
#include "math/maxwellJuettner.hpp"
#include "phys/rel.hpp"
#include "algo/weighter.hpp"

#ifdef SWITCH_MAXWELL_SOLVER_YEE
#include "maxwellSolver/maxwellSolverYee.hpp"
#endif

#include "geometry.hpp"
#include "timeSim.hpp"

using namespace std;

//! define service constant "particle`s vector size"
#define P_VEC_SIZE 13

// getters from particle directly
#define P_POS_R(var) var.pos_r
#define P_POS_PHI(var) var.pos_phi
#define P_POS_Z(var) var.pos_z

#define P_POS_OLD_R(var) var.pos_old_r
#define P_POS_OLD_PHI(var) var.pos_old_phi
#define P_POS_OLD_Z(var) var.pos_old_z

#define P_VEL_R(var) var.vel_r
#define P_VEL_PHI(var) var.vel_phi
#define P_VEL_Z(var) var.vel_z

#define P_WEIGHT(var) var.weight

// service variables to correct cartesian to cylindrical geometry
#define P_SIN(var) var.sin
#define P_COS(var) algo::common::sq_rt(1 - var.sin * var.sin)

#define P_CELL_R(var) var.cell_r
#define P_CELL_Z(var) var.cell_z

#define P_MARK(var) var.mark

struct Particle
{
// getters from particle directly
  double pos_r;
  double pos_phi;
  double pos_z;

  double pos_old_r;
  double pos_old_phi;
  double pos_old_z;

  double vel_r;
  double vel_phi;
  double vel_z;

  double weight;
  double sin;

  size_t cell_r;
  size_t cell_z;
  size_t mark;

  Particle ()
  {
// getters from particle directly
    pos_r = 0;
    pos_phi = 0;
    pos_z = 0;

    pos_old_r = 0;
    pos_old_phi = 0;
    pos_old_z = 0;
    vel_r = 0;
    vel_phi = 0;
    vel_z = 0;
    weight = 0;
    sin = 0;
    cell_r = 0;
    cell_z = 0;
    mark = 0;
  }
};

class FieldE;
class FieldH;
class MaxwellSolver;

class SpecieP
{
public:
  unsigned int id;

  string name;

  unsigned int macro_amount;
  double density[2]; // [left, right]

  // The specie charge in electron charges
  double charge;
  // The specie  *mass
  double mass; // electron masses

  //! Array of particle properties
  //! in format:
  //! \f$ [[r_1, \phi_1, z_1, v_{r_1}, v_{\phi_1}, v_{z_1}, Q_1, m_1, \sin(\theta_r), \cos(\theta_r), is_alive], ...] \f$
  //! So, it is 2D array of particles
  //! which includes arrays of position,
  //! velocity, charge, mass, corrections
  //! and aliveness properties per particle
  //!
  //! You can call: particles[particle_number][component_number]
  //! or use macros to get required component:
  //!
  //! - P_POS_R(particles_variable, particle_number)
  //! - P_POS_PHI(particles_variable, particle_number)
  //! - P_POS_Z(particles_variable, particle_number)
  //! - P_VEL_R(particles_variable, particle_number)
  //! - P_VEL_PHI(particles_variable, particle_number)
  //! - P_VEL_Z(particles_variable, particle_number)
  //! - P_CHARGE(particles_variable, particle_number)
  //! - P_MASS(particles_variable, particle_number)
  //! - P_SIN(particles_variable, particle_number)
  //! - P_COS(particles_variable, particle_number)
  //! - P_ALIVE(particles_variable, particle_number)
  // vector< vector<double> * > particles;
  vector< Particle * > particles;

  Geometry *geometry;

  TimeSim *time;

  double temperature; // in electronvolts
  Grid<double> density_map;
  Grid<double> temperature_map;

  Grid<double> p_abs;
  Grid<double> p_r;
  Grid<double> p_phi;
  Grid<double> p_z;
  Grid<double> count;

  double start_time;
  double bunch_length;
  double bunches_distance;
  unsigned int current_bunch_number;
  int velocity;

public:
  SpecieP() {};

  SpecieP(unsigned int id, // ID for every particles specie
          string p_name,
          double p_charge, // in electron charges
          double p_mass, // in electron masses
          unsigned int p_macro_amount,
          double p_left_density, // CI units
          double p_right_density, // CI units
          double p_temperature, // electronvolts
          Geometry *geom,
          TimeSim *t
    );

  ~SpecieP();

  MaxwellSolver* maxwell_solver;



public:
  virtual void fullyfill_spatial_distribution (void);
  void linear_spatial_distribution (unsigned int int_cell_number,
                                    unsigned int ext_cell_number);

  void rectangular_spatial_distribution(unsigned int int_cell_number,
                                        unsigned int ext_cell_number,
                                        unsigned int left_cell_number,
                                        unsigned int right_cell_number);

  virtual void velocity_distribution ();
  virtual void inject();
  void inject_bunch();
  void boris_pusher();
  void vay_pusher();
  void hc_pusher();
  void mover_cylindrical();
  virtual void reflect();
  void back_position_to_rz();
  void back_velocity_to_rz();
  void dump_position_to_old();
  void bind_cell_numbers ();

  void calc_density ();
  void calc_temperature ();

protected:
  void rectangular_random_placement (unsigned int int_cell_number,
                                     unsigned int ext_cell_number,
                                     unsigned int left_cell_number,
                                     unsigned int right_cell_number);
  void rectangular_flat_placement (unsigned int int_cell_number,
                                   unsigned int ext_cell_number,
                                   unsigned int left_cell_number,
                                   unsigned int right_cell_number);
  void rectangular_regular_placement (unsigned int int_cell_number,
                                      unsigned int ext_cell_number,
                                      unsigned int left_cell_number,
                                      unsigned int right_cell_number);
  void rectangular_centered_placement (unsigned int int_cell_number,
                                       unsigned int ext_cell_number,
                                       unsigned int left_cell_number,
                                       unsigned int right_cell_number);
  void set_weight (unsigned int int_cell_number,
                   unsigned int ext_cell_number,
                   unsigned int left_cell_number,
                   unsigned int right_cell_number);

  void thermal_velocity_distribution ();
  void rectangular_velocity_distribution ();
  void eigen_velocity_distribution ();
  void eigen_directed_velocity_distribution (unsigned int dir); // 0,1,2 are for r, phi, z
};

#endif // end of _SPECIE_P_HPP_

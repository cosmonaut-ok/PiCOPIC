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

#ifndef _SPECIE_P_HPP_
#define _SPECIE_P_HPP_

#include "constant.hpp"
#include "phys/rel.hpp"
#include "lib.hpp"
#include "msg.hpp"
#include "grid.hpp"
#include "geometry.hpp"
#include "fieldE.hpp"
#include "fieldH.hpp"
#include "timeSim.hpp"

#include "math/rand.hpp"
#include "math/maxwellJuttner.hpp"
#include "math/vector3d.hpp"

//! define service constant "particle`s vector size"
#define P_VEC_SIZE 13

// getters from particles vector
#define PP_POS_R(var, p_num) (*var[p_num])[0]
#define PP_POS_PHI(var, p_num) (*var[p_num])[1]
#define PP_POS_Z(var, p_num) (*var[p_num])[2]

#define PP_POS_OLD_R(var, num) (*var[num])[3]
#define PP_POS_OLD_PHI(var, num) (*var[num])[4]
#define PP_POS_OLD_Z(var, num) (*var[num])[5]

#define PP_VEL_R(var, num) (*var[num])[6]
#define PP_VEL_PHI(var, num) (*var[num])[7]
#define PP_VEL_Z(var, num) (*var[num])[8]

#define PP_WEIGHT(var, num) (*var[num])[9]

// service variables to correct cartesian to cylindrical geometry
#define PP_SIN(var, num) (*var[num])[10]
#define PP_COS(var, num) lib::sq_rt( 1 - (*var[num])[10] * (*var[num])[10] )

#define PP_CELL_R(var, num) (*var[num])[11]
#define PP_CELL_Z(var, num) (*var[num])[12]

// getters from particle directly
#define P_POS_R(var) var[0]
#define P_POS_PHI(var) var[1]
#define P_POS_Z(var) var[2]

#define P_POS_OLD_R(var) var[3]
#define P_POS_OLD_PHI(var) var[4]
#define P_POS_OLD_Z(var) var[5]

#define P_VEL_R(var) var[6]
#define P_VEL_PHI(var) var[7]
#define P_VEL_Z(var) var[8]

#define P_WEIGHT(var) var[9]

// service variables to correct cartesian to cylindrical geometry
#define P_SIN(var) var[10]
#define P_COS(var) lib::sq_rt(1 - var[10] * var[10])

#define P_CELL_R(var) var[11]
#define P_CELL_Z(var) var[12]

class FieldE;
class FieldH;

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
  vector< vector<double> * > particles;
  vector< vector<double> * > particles_old;

  Geometry *geometry;

  TimeSim *time;

  double temperature; // in electronvolts

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

  // void set_spatial_distribution(){};
  // void set_velocity_distribution(){};
// private:
  FieldH* field_h;
  FieldE* field_e;

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

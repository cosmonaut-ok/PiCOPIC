#pragma once

#include "constant.hpp"
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

// define some service constants
// use classical calculations, if velocity lower, than minimal
#ifdef REL_LIMIT
// define REL_LIMIT^2 to decrease number of operations
  const double REL_LIMIT_POW_2 = pow (REL_LIMIT, 2);
#endif

// C^2 define c^2 to decrease number of operations
const double LIGHT_VEL_POW_2 = pow (LIGHT_VEL, 2);

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

#define PP_CHARGE(var, num) (*var[num])[9]
#define PP_MASS(var, num) (*var[num])[10]

// service variables to correct cartesian to cylindrical geometry
#define PP_SIN(var, num) (*var[num])[11]
#define PP_COS(var, num) (*var[num])[12]

// #define PP_ALIVE(var, num) (*var[num])[13]

// #define PP_JUMP(var, num) (*var[num])[14]

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

#define P_CHARGE(var) var[9]
#define P_MASS(var) var[10]

// service variables to correct cartesian to cylindrical geometry
#define P_SIN(var) var[11]
#define P_COS(var) var[12]

// #define P_ALIVE(var) var[13]

// #define P_JUMP(var) var[14]

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
  int current_bunch_number;
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
  virtual void fullyfill_spatial_distribution();
  void linear_spatial_distribution(unsigned int int_cell_number,
                                   unsigned int ext_cell_number);
  void rectangular_spatial_distribution(unsigned int int_cell_number,
                                        unsigned int ext_cell_number,
                                        unsigned int left_cell_number,
                                        unsigned int right_cell_number);
  virtual void thermal_velocity_distribution ();
  void rectangular_velocity_distribution ();
  void eigen_velocity_distribution ();
  void eigen_directed_velocity_distribution (unsigned int dir);
  virtual void wakeup();
  virtual void inject();
  void inject_bunch();
  void boris_pusher();
  void half_step_mover_cylindrical();
  virtual void reflect();
  void back_position_to_rz();
  void back_velocity_to_rz();

  virtual bool hes_dead_jim() { return false; }; // check if particles bunch is dead

  void dump_position_to_old();
};

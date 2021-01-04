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

#ifndef _CFG_HPP_
#define _CFG_HPP_

#include <vector>
#include <sstream> // std::stringstream
#include <fstream> // std::ifstream

#include "picojson.h"

#include "defines.hpp"
#include "msg.hpp"
#include "phys/plasma.hpp"

#include "geometry.hpp"
#include "timeSim.hpp"

#define PATH_DELIMITER "/"
#define RANGE_DELIMITER "-"
#define SPACE_DELIMITER "_"

#define SPECIE_HASH_SALT 1394
#define BEAM_HASH_SALT 18132

#define MIN_DOMAIN_GRID_AMOUNT 64

struct probe
{
  std::string path;
  std::string component;
  std::string specie;
  unsigned int shape;
  std::vector<size_t> dims;
  unsigned int schedule;
};

struct particle_specie
{
  unsigned int id;
  std::string name;
  unsigned int mass;
  int charge;
  double macro_amount;
  double left_density;
  double right_density;
  double temperature;
};

struct particle_beam : particle_specie
{
  double velocity;
  double density;
  int bunches_amount;
  double start_time;
  double bunches_distance;
  double bunch_length;
  double bunch_radius;
  unsigned int current_bunch_number; // to mark, how many bunches already injcted
};

struct save_data
{
  char *path;
  char *data_root;
  int fpf; // frames per file
  bool compress;
  int compress_level;
};

class Cfg
{
public:
  // Cfg(void);
  Cfg(const std::string json_file_name);
  ~Cfg(void);

  std::string cfg2str();
  picojson::value cfg2value();

public:
  char *log_file;

  bool use_hdf5 = false;

  double macro_amount; // amount of macroparticles in system
  unsigned int beam_macro_ratio;

  Geometry *geometry;

  /* <time> */
  TimeSim *time;

  /* output data params structure */
  save_data *output_data;

  vector<particle_specie> particle_species;
  vector<particle_beam> particle_beams;

  //! use "-particle-" infix for physical particles
  //! use "-macro-" infix for macro- or superparticles
  //! use "beam" for full set of electron bunches
  //! use "bunch" for sinble electron bunch (part of beam)

//   char *beam_name;
//   //! particles charge in beam
//   //! (mean, all bunches in the beam and charge/mass of all particles
//   //! in the bunch are equal.
//   int beam_particle_charge;
//   //! particles mass in beam
//   unsigned int beam_particle_mass;
//   //! number of particles in beam
//   unsigned int beam_number_bunches;
//   //! distance between bunches in beam
//   double beam_bunches_distance;
//   //! initial velocity of bunches in beam
//   //! (mean, it is equal for all bucnhes in beam)
//   double beam_initial_velocity;
//   //! number of macroparticles in bunch
//   unsigned int bunch_number_macro;
//   //! length of bunch
//   double bunch_lenght;
//   //! radius of bunch
//   double bunch_radius;
//   //! particles density in bunch
//   double bunch_density;

//   double boundary_maxwell_e_phi_upper;
//   double boundary_maxwell_e_phi_left;
//   double boundary_maxwell_e_phi_right;
//   int boundary_conditions;

//   char *dump_result_path;
//   char *dump_data_root;
//   unsigned int dump_data_interval;
//   unsigned int dump_frames_per_file;
//   unsigned int dump_system_state_interval;
//   bool dump_compress = false;
//   int dump_compress_level = 0;

//   // dump different kinds of data
//   bool dump_e_r = false;
//   bool dump_e_phi = false;
//   bool dump_e_z = false;
//   bool dump_h_r = false;
//   bool dump_h_phi = false;
//   bool dump_h_z = false;
//   bool dump_position = false;
//   bool dump_velocity = false;
//   bool dump_rho_beam = false;

//   // electrical current dump
//   bool dump_current_r = false;
//   bool dump_current_phi = false;
//   bool dump_current_z = false;

  // probes
  vector<probe> probes;

private:
  picojson::value json_data;

  void init_particles();
  void init_probes();
  void init_beam();
  void init_geometry();
  void init_pml();
  void init_time();
  void init_boundary();
  void init_output_data();
  void weight_macro_amount();
  bool method_limitations_check();
};
#endif // end of _CFG_HPP_

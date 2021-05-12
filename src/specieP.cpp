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

#include "specieP.hpp"
#include "geometry.hpp"

using namespace constant;

SpecieP::SpecieP (unsigned int p_id,
                  string p_name,
                  double p_charge, // in electron charges
                  double p_mass, // in electron masses
                  unsigned int p_macro_amount,
                  double p_left_density, // CI units
                  double p_right_density, // CI units
                  double p_temperature, // electronvolts
                  Geometry *geom,
                  TimeSim *t
  ):geometry(geom),time(t)
{
  id = p_id;

  name = p_name;

  // distribute the particles
  charge = p_charge * EL_CHARGE;
  mass = p_mass * EL_MASS;

/// Correct plasma macroparticles amount
/// to satisfy conditions of regular/centered
/// spatial distributions
#if defined(SWITCH_PLASMA_SPATIAL_REGULAR) || defined(SWITCH_PLASMA_SPATIAL_CENTERED)
  LOG_S(MAX) << "Correcting domain plasma particles macro amount to satisfy spatial distribution";
  double r_size = geometry->size[0];
  double z_size = geometry->size[1];
  double dr = geometry->cell_size[0];
  double dz = geometry->cell_size[1];

  macro_amount = algo::common::nearest_divide(p_macro_amount, r_size / dr * z_size / dz);
#else
  macro_amount = p_macro_amount;
#endif // WITH_PLASMA_SPATIAL

  density[0] = p_left_density;
  density[1] = p_right_density;

  temperature = p_temperature;

  density_map = Grid<double> (geometry->cell_amount[0], geometry->cell_amount[1], 2);
  temperature_map = Grid<double> (geometry->cell_amount[0], geometry->cell_amount[1], 2);
  energy_map = Grid<double> (geometry->cell_amount[0], geometry->cell_amount[1], 2);

  momentum_abs_map_w = Grid<double> (geometry->cell_amount[0], geometry->cell_amount[1], 2);
  momentum_map_w = Grid3D<double> (geometry->cell_amount[0], geometry->cell_amount[1], 2);
  momentum_map = Grid3D<double> (geometry->cell_amount[0], geometry->cell_amount[1], 2);
}

SpecieP::~SpecieP()
{
  for (auto i = particles.begin(); i != particles.end(); i++)
    delete *i;
  particles.clear();
}

void SpecieP::fullyfill_spatial_distribution()
// ! apply flat spatial macroparticles distribution
// ! with gradient, declared for object
// ! gradient achived with macroparticles
// ! mass differentiation
{
  rectangular_spatial_distribution(geometry->cell_dims[0],
                                   geometry->cell_dims[2],
                                   geometry->cell_dims[1],
                                   geometry->cell_dims[3]);
}

void SpecieP::linear_spatial_distribution(unsigned int int_cell_number,
                                          unsigned int ext_cell_number)
// ! spatial distribution for linear (for cartese coordinates)
// ! or tube-shaped macroparticles cloud
{
  rectangular_spatial_distribution(int_cell_number, ext_cell_number,
                                   0, geometry->size[1]);
}

void SpecieP::rectangular_spatial_distribution(unsigned int int_cell_number,
                                               unsigned int ext_cell_number,
                                               unsigned int left_cell_number,
                                               unsigned int right_cell_number)
// ! spatial distribution for rectangular (for cartese coordinates)
// ! or cylindrical ring-shaped macroparticles cloud
{
#ifdef SWITCH_PLASMA_SPATIAL_CENTERED
  rectangular_centered_placement (int_cell_number,
                                  ext_cell_number,
                                  left_cell_number,
                                  right_cell_number);
#elif defined(SWITCH_PLASMA_SPATIAL_FLAT) || defined(SWITCH_PLASMA_SPATIAL_RANDOM)
  rectangular_random_placement (int_cell_number,
                                ext_cell_number,
                                left_cell_number,
                                right_cell_number);
#elif defined(SWITCH_PLASMA_SPATIAL_REGULAR)
  rectangular_regular_placement (int_cell_number,
                                 ext_cell_number,
                                 left_cell_number,
                                 right_cell_number);
#endif // WITH_PLASMA_SPATIAL

  set_weight(int_cell_number, ext_cell_number, left_cell_number, right_cell_number);
}

void SpecieP::rectangular_random_placement (unsigned int int_cell_number,
                                            unsigned int ext_cell_number,
                                            unsigned int left_cell_number,
                                            unsigned int right_cell_number)
{
  //! macroparticles spatial placement
  //! for "random" and "flat" cases
  double dr = geometry->cell_size[0];
  double dz = geometry->cell_size[1];
  double dn = density[1] - density[0];
  double dl = density[0]; // dl is for "density left"

  double r_size = (ext_cell_number - int_cell_number) * dr;
  double z_size = (right_cell_number - left_cell_number) * dz;

  double rand_r, rand_z;

  // decrease size at r=r and z=z walls
  // this caused by particles formfactor
  if (geometry->walls[2]) r_size -= dr;
  if (geometry->walls[3]) z_size -= dz;

#ifdef SWITCH_PLASMA_SPATIAL_FLAT
  unsigned int macro_count = 0;
#endif

  for (unsigned int i = 0; i < macro_amount; i++)
  {
    Particle *n = new Particle();
    P_SPECIE_ID((*n)) = id;

#ifdef SWITCH_PLASMA_SPATIAL_RANDOM
    rand_r = math::random::uniform1();
    rand_z = math::random::uniform1();
#elif defined(SWITCH_PLASMA_SPATIAL_FLAT)
    rand_r = math::random::random_reverse(macro_count, 13);
    rand_z = math::random::random_reverse(macro_amount - 1 - macro_count, 11);
#endif // WITH_PLASMA_SPATIAL
    if (rand_r == 0) rand_r = MNZL;
    if (rand_z == 0) rand_z = MNZL;

    P_POS_R((*n)) = r_size * rand_r;

    P_POS_PHI((*n)) = 0;

    if (density[0] == density[1])
      P_POS_Z((*n)) = z_size * rand_z;
    else
      P_POS_Z((*n)) = z_size / dn
        * (algo::common::sq_rt(pow(dl, 2) + rand_z * (2 * dl * dn + pow(dn, 2))) - dl);

    P_POS_R((*n)) += int_cell_number * dr;
    P_POS_R((*n)) += dr / 2.;
    P_POS_Z((*n)) += left_cell_number * dz; // shift by z to respect geometry with domains
    P_POS_Z((*n)) += dz / 2.;

#ifdef SWITCH_PLASMA_SPATIAL_FLAT
    ++macro_count;
#endif

    particles.push_back(n);
  }
}

void SpecieP::rectangular_regular_placement (unsigned int int_cell_number,
                                             unsigned int ext_cell_number,
                                             unsigned int left_cell_number,
                                             unsigned int right_cell_number)
{
  //! macroparticles spatial placement
  //! for "regular" case

  double dr = geometry->cell_size[0];
  double dz = geometry->cell_size[1];

  double r_size = (ext_cell_number - int_cell_number) * dr;
  double z_size = (right_cell_number - left_cell_number) * dz;

  double r_macro = algo::common::sq_rt(macro_amount * r_size / z_size);
  double z_macro = algo::common::sq_rt(macro_amount * z_size / r_size);

  double r_macro_interval = r_size / r_macro;
  double z_macro_interval = z_size / z_macro;

  // decrease size at r=r and z=z walls
  // this caused by particles formfactor
  if (geometry->walls[2]) r_size -= dr;
  if (geometry->walls[3]) z_size -= dz;

  for (double i = MNZL; i <= r_size - r_macro_interval; i += r_macro_interval)
    for (double j = MNZL; j <= z_size - z_macro_interval; j += z_macro_interval)
    {
      Particle *n = new Particle();
      P_SPECIE_ID((*n)) = id;

      P_POS_R((*n)) = i;
      P_POS_Z((*n)) = j;

      P_POS_R((*n)) += int_cell_number * dr;
      P_POS_R((*n)) += dr / 2.;
      P_POS_Z((*n)) += left_cell_number * dz; // shift by z to respect geometry with domains
      P_POS_Z((*n)) += dz / 2.;

      particles.push_back(n);
    }
}

void SpecieP::rectangular_centered_placement (unsigned int int_cell_number,
                                              unsigned int ext_cell_number,
                                              unsigned int left_cell_number,
                                              unsigned int right_cell_number)
{
  //! macroparticles spatial placement
  //! for "centered" case

  double dr = geometry->cell_size[0];
  double dz = geometry->cell_size[1];
  double r_cells = (ext_cell_number - int_cell_number);
  double z_cells = (right_cell_number - left_cell_number);

  // decrease size at r=r and z=z walls
  // this caused by particles formfactor
  if (geometry->walls[2]) r_cells -= 1;
  if (geometry->walls[3]) z_cells -= 1;

  unsigned int macro_per_cell = floor(macro_amount / (r_cells * z_cells));

  for (unsigned int rc = 0; rc < r_cells; ++rc)
    for (unsigned int zc = 0; zc < z_cells; ++zc)
    {
      unsigned int macro_counter = 0;
      while (macro_counter < macro_per_cell)
      {
        double r_place = dr * rc + MNZL;
        double z_place = dz * zc + MNZL;

        Particle *n = new Particle();
	P_SPECIE_ID((*n)) = id;

        P_POS_R((*n)) = r_place;
        P_POS_Z((*n)) = z_place;

        P_POS_R((*n)) += int_cell_number * dr;
        P_POS_R((*n)) += dr / 2.;
        P_POS_Z((*n)) += left_cell_number * dz; // shift by z to respect geometry with domains
        P_POS_Z((*n)) += dz / 2.;

        particles.push_back(n);
        ++macro_counter;
      }
    }
}

void SpecieP::set_weight (unsigned int int_cell_number,
                          unsigned int ext_cell_number,
                          unsigned int left_cell_number,
                          unsigned int right_cell_number)
{
  double dr = geometry->cell_size[0];
  double dz = geometry->cell_size[1];
  double r_size = (ext_cell_number - int_cell_number) * dr;
  double z_size = (right_cell_number - left_cell_number) * dz;

  // decrease size at r=r and z=z walls
  // this caused by particles formfactor
  if (geometry->walls[2]) r_size -= dr;
  if (geometry->walls[3]) z_size -= dz;

  double v_sum = 0; // summary volume of all particles
  for (auto n = particles.begin(); n != particles.end(); ++n)
    if (P_POS_R((**n)) > dr / 2)
      v_sum += 2 * PI * P_POS_R((**n)) * dr * dz;
    else
      v_sum += PI * P_POS_R((**n)) * P_POS_R((**n)) * dz;

  double v_avg = v_sum / macro_amount;

  // set top and bottom radius cell numbers (decrease at r=r wall)
  double top_r_num = geometry->cell_dims[2];
  double bot_r_num = geometry->cell_dims[0];
  if (geometry->walls[2]) top_r_num -= 1;

  double N_total = (density[0] + density[1]) / 2
    * PI * dr * dr * (top_r_num * top_r_num - bot_r_num * bot_r_num) * z_size;

  double n_per_macro_avg = N_total / macro_amount;

  double norm;

  // for (unsigned int n = 0; n < macro_amount; n++)
  for (auto n = particles.begin(); n != particles.end(); ++n)
  {
    // coefitient of normalization
    if (P_POS_R((**n)) > dr / 2)
      norm = 2 * PI * P_POS_R((**n)) * dr * dz / v_avg;
    else
      norm = PI * P_POS_R((**n)) * P_POS_R((**n)) * dz / v_avg;

    // number of real particles per macroparticle
    double n_per_macro = n_per_macro_avg * norm;

    // set charge and mass of macroparticle
    P_WEIGHT((**n)) = n_per_macro;
  }
}

void SpecieP::velocity_distribution ()
{
#ifdef SWITCH_PLASMA_VELOCITY_THERMAL
  thermal_velocity_distribution();
#elif defined(SWITCH_PLASMA_VELOCITY_RECTANGULAR)
  rectangular_velocity_distribution();
#elif defined(SWITCH_PLASMA_VELOCITY_EIGEN)
  eigen_velocity_distribution();
#else
  LOG_S(FATAL) << "Plasma velocity distribution type unknown";
#endif
}

void SpecieP::thermal_velocity_distribution ()
{
  // Sample the energies in the MJ distribution
  vector<double> energies = math::maxwell_juettner::maxwellJuettner(macro_amount, temperature);

  unsigned int macro_count = 0;

  double two_over_mass = 2 * EL_CHARGE / mass; // '* EL_CHARGE' - convert eV to J
  double mc_inv = EL_CHARGE / (mass * LIGHT_VEL); // 'EL_CHARGE /' - convert eV to J

  // normalize to 1 / sqrt(2).
  // TODO: I don't know, why, but it returns temperature correct values
  const double norm = 0.7071067811865475;

  for (auto p = particles.begin(); p != particles.end(); ++p)
  {
    double therm_vel_cmp = energies[macro_count] * mc_inv;

    if (therm_vel_cmp > REL_LIMIT)
    {
      double gamma_inv = phys::rel::lorenz_factor_inv(pow(therm_vel_cmp, 2));
      therm_vel_cmp *= gamma_inv;
    }
    else
      therm_vel_cmp = algo::common::sq_rt(energies[macro_count] * two_over_mass);

    therm_vel_cmp *= norm;

    double rnd_0 = math::random::uniform2();
    double rnd_1 = math::random::uniform2();
    double rnd_2 = math::random::uniform2();

    P_VEL_R((**p)) = rnd_0 * therm_vel_cmp;
    P_VEL_PHI((**p)) = rnd_1 * therm_vel_cmp;
    P_VEL_Z((**p)) = rnd_2 * therm_vel_cmp;

    ++macro_count;
  }
}

void SpecieP::rectangular_velocity_distribution ()
{
  double two_over_mass = 2 * EL_CHARGE / mass; // '* EL_CHARGE' - convert eV to J
  double mc_inv = EL_CHARGE / (mass * LIGHT_VEL); // 'EL_CHARGE /' - convert eV to J

  double therm_vel_cmp = temperature * mc_inv;
  therm_vel_cmp *= 2; // because we want a medium value for temperature

  if (therm_vel_cmp > REL_LIMIT)
  {
    double gamma_inv = phys::rel::lorenz_factor_inv(pow(therm_vel_cmp, 2));
    therm_vel_cmp *= gamma_inv;
  }
  else
    therm_vel_cmp = algo::common::sq_rt(temperature * two_over_mass);

  for (auto p = particles.begin(); p != particles.end(); ++p)
  {
    double rnd_0, rnd_1, rnd_2;

    rnd_0 = math::random::uniform2();
    rnd_1 = math::random::uniform2();
    rnd_2 = math::random::uniform2();

    P_VEL_R((**p)) = rnd_0 * therm_vel_cmp;
    P_VEL_PHI((**p)) = rnd_1 * therm_vel_cmp;
    P_VEL_Z((**p)) = rnd_2 * therm_vel_cmp;
  }
}

void SpecieP::eigen_velocity_distribution ()
// ! singular velocity distribution
{
  double two_over_mass = 2 * EL_CHARGE / mass; // '* EL_CHARGE' - convert eV to J
  double mc_inv = EL_CHARGE / (mass * LIGHT_VEL); // 'EL_CHARGE /' - convert eV to J

  double therm_vel_cmp = temperature * mc_inv;
  therm_vel_cmp *= 2; // because we want a medium value for temperature

  if (therm_vel_cmp > REL_LIMIT)
  {
    double gamma_inv = phys::rel::lorenz_factor_inv(pow(therm_vel_cmp, 2));
    therm_vel_cmp *= gamma_inv;
  }
  else
    therm_vel_cmp = algo::common::sq_rt(temperature * two_over_mass);

  for (auto p = particles.begin(); p != particles.end(); ++p)
  {
    P_VEL_R((**p)) = therm_vel_cmp;
    P_VEL_PHI((**p)) = therm_vel_cmp;
    P_VEL_Z((**p)) = therm_vel_cmp;
  }
}

void SpecieP::eigen_directed_velocity_distribution (unsigned int dir)
// ! directed singular velocity distribution
// ! 0 is for "R" direction
// ! 1 is for "PHI" direction
// ! 2 is for "Z" direction
{
  double two_over_mass = 2 * EL_CHARGE / mass; // '* EL_CHARGE' - convert eV to J
  double mc_inv = EL_CHARGE / (mass * LIGHT_VEL); // 'EL_CHARGE /' - convert eV to J

  double therm_vel_cmp = temperature * mc_inv;

  if (therm_vel_cmp > REL_LIMIT)
  {
    double gamma_inv = phys::rel::lorenz_factor_inv(pow(therm_vel_cmp, 2));
    therm_vel_cmp *= gamma_inv;
  }
  else
    therm_vel_cmp = algo::common::sq_rt(temperature * two_over_mass);

  switch (dir)
  {
  case 0:
    for (auto p = particles.begin(); p != particles.end(); ++p)
      P_VEL_R((**p)) = therm_vel_cmp;
    break;
  case 1:
    for (auto p = particles.begin(); p != particles.end(); ++p)
      P_VEL_PHI((**p)) = therm_vel_cmp;
    break;
  case 2:
    for (auto p = particles.begin(); p != particles.end(); ++p)
      P_VEL_Z((**p)) = therm_vel_cmp;
    break;
  default:
    LOG_S(FATAL) << "Incorrect switch of rectangular directed velocity component: " << dir;
    break;
  }
}

void SpecieP::inject () {}

void SpecieP::mover_cylindrical()
{
  for (auto p = particles.begin(); p != particles.end(); ++p)
  {
    P_POS_R((**p)) = P_POS_R((**p)) + P_VEL_R((**p)) * time->step;
    //! we use "fake" rotation component to correct position from xy to rz pane
    P_POS_PHI((**p)) = P_POS_PHI((**p)) + P_VEL_PHI((**p)) * time->step;
    P_POS_Z((**p)) = P_POS_Z((**p)) + P_VEL_Z((**p)) * time->step;
  }
}

void SpecieP::reflect ()
{
  double dr = geometry->cell_size[0];
  double dz = geometry->cell_size[1];
  double half_dr = dr / 2.;
  double half_dz = dz / 2.;

  double radius_wall = geometry->size[0] - dr / 2.;
  double longitude_wall = geometry->size[1] - dz / 2.;

  double radius_wallX2 = radius_wall * 2.;
  double longitude_wallX2 = longitude_wall * 2.;

  // shift for converting local positions into global and back
  double r_shift = geometry->cell_dims[0] * dr;
  double z_shift = geometry->cell_dims[1] * dz;

  for (auto p = particles.begin(); p != particles.end(); ++p)
  {
    // set temporary position as it located in domain 0,0
    double pos_r = P_POS_R((**p)) - r_shift;
    double pos_z = P_POS_Z((**p)) - z_shift;
    double pos_old_r = P_POS_OLD_R((**p)) - r_shift;
    double pos_old_z = P_POS_OLD_Z((**p)) - z_shift;

    double pos_delta_r = pos_r - pos_old_r;
    double pos_delta_z = pos_z - pos_old_z;

    while ( // catch multiple reflections
      isnormal(pos_r) && isnormal(pos_z) && // ensure, that position components are not NANs
      (
        (pos_r > radius_wall && geometry->walls[2])
        || (pos_z > longitude_wall && geometry->walls[3])
        || (pos_r < half_dr && geometry->walls[0])
        || (pos_z < half_dz && geometry->walls[1])
        )
      )
    {
      if (pos_r > radius_wall && geometry->walls[2])
      {
        P_VEL_R((**p)) = - P_VEL_R((**p));
        double dr_small = radius_wall - pos_old_r;
        double dz_small = dr_small * pos_delta_z / pos_delta_r;

        // new set position
        pos_r = radius_wallX2 - pos_r;

        // update old_position and deltas
        pos_delta_r = pos_delta_r - dr_small;
        pos_delta_z = pos_delta_z - dz_small;
        pos_old_r = radius_wall;
        pos_old_z = pos_old_z - dz_small;
      }

      if (pos_z > longitude_wall && geometry->walls[3])
      {
        P_VEL_Z((**p)) = - P_VEL_Z((**p));

        double dz_small = longitude_wall - pos_old_z;
        double dr_small = dz_small * pos_delta_r / pos_delta_z;

        // new set position
        pos_z = longitude_wallX2 - pos_z;

        // update old_position and deltas
        pos_delta_r = pos_delta_r - dr_small;
        pos_delta_z = pos_delta_z - dz_small;
        pos_old_z = longitude_wall;
        pos_old_r = pos_old_r - dr_small;
      }

      if (pos_r < half_dr && geometry->walls[0])
      {
        P_VEL_R((**p)) = - P_VEL_R((**p));

        double dr_small = half_dr - pos_old_r;
        double dz_small = dr_small * pos_delta_z / pos_delta_r;

        // new set position
        pos_r = dr - pos_r;

        // update old_position and deltas
        pos_delta_r = pos_delta_r - dr_small;
        pos_delta_z = pos_delta_z - dz_small;
        pos_old_r = half_dr;
        pos_old_z = pos_old_z - dz_small;

        // fix change -small to +small
        if (pos_r < half_dr) pos_r = half_dr;
      }

      if (pos_z < half_dz && geometry->walls[1])
      {
        P_VEL_Z((**p)) = - P_VEL_Z((**p));

        double dz_small = half_dz - pos_old_z;
        double dr_small = dz_small * pos_delta_r / pos_delta_z;

        // new set position
        pos_z = dz - pos_z;

        // update old_position and deltas
        pos_delta_r = pos_delta_r - dr_small;
        pos_delta_z = pos_delta_z - dz_small;
        pos_old_z = dz;
        pos_old_r = pos_old_r - dr_small;

        // fix change -small to +small
        if (pos_z < half_dz) pos_z = half_dz;
      }
    }
    P_POS_R((**p)) = pos_r + r_shift;
    P_POS_Z((**p)) = pos_z + z_shift;
  }
}

void SpecieP::back_position_to_rz()
{
  // ! implementation of backing coodrinates to rz pane
  // ! taken from https: // www.particleincell.com / 2015 / rz-pic /
  for (auto p = particles.begin(); p != particles.end(); ++p)
  {
    double pos_r = P_POS_R((**p));
    double pos_phi = P_POS_PHI((**p));
    double r = algo::common::sq_rt(pos_r * pos_r + pos_phi * pos_phi);
    P_SIN((**p)) = pos_phi / r;
    P_POS_R((**p)) = r;
    P_POS_PHI((**p)) = 0;
  }
}

void SpecieP::back_velocity_to_rz()
{
  // ! implementation of backing coodrinates to rz pane
  // ! taken from https: // www.particleincell.com / 2015 / rz-pic /
  for (auto p = particles.begin(); p != particles.end(); ++p)
  {
    double sin = P_SIN((**p));
    double cos = P_COS((**p));
    double v_r = P_VEL_R((**p));
    double v_phi = P_VEL_PHI((**p));

    double u_2 = cos * v_r - sin * v_phi;
    double v_2 = sin * v_r + cos * v_phi;
    P_VEL_R((**p)) = u_2;
    P_VEL_PHI((**p)) = v_2;
  }
}

void SpecieP::dump_position_to_old()
{
  for (auto p = particles.begin(); p != particles.end(); ++p)
  {
    P_POS_OLD_R((**p)) = P_POS_R((**p));
    P_POS_OLD_PHI((**p)) = P_POS_PHI((**p));
    P_POS_OLD_Z((**p)) = P_POS_Z((**p));
  }
}

void SpecieP::bind_cell_numbers()
{
  for (auto p = particles.begin(); p != particles.end(); ++p)
  {
    // renumerate cells
    int r_cell = CELL_NUMBER(P_POS_R((**p)), geometry->cell_size[0]);
    int z_cell = CELL_NUMBER(P_POS_Z((**p)), geometry->cell_size[1]);

    P_CELL_R((**p)) = r_cell;
    P_CELL_Z((**p)) = z_cell;
  }
}

void SpecieP::calc_density()
{ // FIXME: it weights density in both cases
  // clear grid values
  density_map = 0;
  density_map.overlay_set(0);

  for (auto i = particles.begin(); i != particles.end(); i++)
    weight_cylindrical<double> ( geometry, &density_map,
                                 P_POS_R((**i)),
                                 P_POS_Z((**i)),
                                 P_WEIGHT((**i)));
}

void SpecieP::calc_temperature()
{
  // clear grid values
  temperature_map = 0;
  temperature_map.overlay_set(0);

  // calculate density initially
  calc_density();
  calc_momentum_w(); // calculates absolute momentum values also

  for (int r = 0; r < geometry->cell_amount[0]; r++)
    for (int z = 0; z < geometry->cell_amount[1]; z++)
    {
      double density_local = density_map(r, z);

      if (density_local == 0)
        temperature_map.set(r, z, 0);
      else
      {
        double p_vec_sum_2 =
          pow(momentum_map_w[0](r, z), 2)
          + pow(momentum_map_w[1](r, z), 2)
          + pow(momentum_map_w[2](r, z), 2);

        double p_sc_sum_2 = pow(momentum_abs_map_w(r, z), 2) - p_vec_sum_2;

        // normalize to count/density and convert Joules to eV
        p_sc_sum_2 /= pow(density_local, 2);

        double energy = phys::rel::energy_m(mass, p_sc_sum_2);

        // convert Joules to eV
        energy *= EL_CHARGE_INV;

        temperature_map.set(r, z, energy);
      }
    }
}

void SpecieP::calc_momentum_abs_w()
{ // clear grid values
  momentum_abs_map_w = 0;
  momentum_abs_map_w.overlay_set(0);

  for (auto i = particles.begin(); i != particles.end(); i++)
  {
    double v_abs = algo::common::sq_rt (
      pow(P_VEL_R((**i)), 2)
      + pow(P_VEL_PHI((**i)), 2)
      + pow(P_VEL_Z((**i)), 2)
      );

    double p_abs = mass * P_WEIGHT((**i)) * v_abs;

    if (v_abs < REL_LIMIT)
    {
      double gamma = phys::rel::lorenz_factor(pow(v_abs, 2));
      p_abs *= gamma;
    }

    weight_cylindrical<double>(geometry, &momentum_abs_map_w,
                               P_POS_R((**i)),
                               P_POS_Z((**i)),
                               p_abs);
  }
}

void SpecieP::calc_momentum_w()
{ // clear momentum grid
  momentum_map_w = 0;
  momentum_map_w.overlay_set(0);

  for (auto i = particles.begin(); i != particles.end(); i++)
  {
    double v_r = P_VEL_R((**i));
    double v_phi = P_VEL_PHI((**i));
    double v_z = P_VEL_Z((**i));

    double m_weighted = mass * P_WEIGHT((**i));

    double p_r = m_weighted * v_r;
    double p_phi = m_weighted * v_phi;
    double p_z = m_weighted * v_z;

    // check for relativistic case
    double v_abs_2 = pow(v_r, 2) + pow(v_phi, 2) + pow(v_z, 2);
    if (v_abs_2 < REL_LIMIT_POW_2)
    {
      double gamma = phys::rel::lorenz_factor(v_abs_2);
      p_r *= gamma;
      p_phi *= gamma;
      p_z *= gamma;
    }

    weight_cylindrical<double>(geometry, &(momentum_map_w[0]),
                               P_POS_R((**i)),
                               P_POS_Z((**i)),
                               p_r);

    weight_cylindrical<double>(geometry, &(momentum_map_w[1]),
                               P_POS_R((**i)),
                               P_POS_Z((**i)),
                               p_phi);

    weight_cylindrical<double>(geometry, &(momentum_map_w[2]),
                               P_POS_R((**i)),
                               P_POS_Z((**i)),
                               p_z);
  }
}

void SpecieP::calc_energy()
{
  // clear grid values
  energy_map = 0;
  energy_map.overlay_set(0);

  // calculate density initially
  calc_density();
  calc_momentum_abs_w();

  for (int r = 0; r < geometry->cell_amount[0]; r++)
    for (int z = 0; z < geometry->cell_amount[1]; z++)
    {
      double density_local = density_map(r, z);

      if (density_local == 0)
        energy_map.set(r, z, 0);
      else
      {
        double p_sc_sum_2 = pow(momentum_abs_map_w(r, z), 2);

        // normalize to count/density and convert Joules to eV
        p_sc_sum_2 /= pow(density_local, 2);

        double energy = phys::rel::energy_m(mass, p_sc_sum_2);

        // convert Joules to eV
        energy *= EL_CHARGE_INV;

        energy_map.set(r, z, energy);
      }
    }
}

void SpecieP::calc_momentum ()
{
  calc_density();
  calc_momentum_w(); // calculates absolute momentum values also

  momentum_map = 0;
  momentum_map.overlay_set(0);

  for (unsigned int i = 0; i < geometry->cell_amount[0]; ++i)
    for (unsigned int j = 0; j < geometry->cell_amount[1]; ++j)
      if (density_map(i, j) != 0)
      {
        momentum_map[0].set(i, j, momentum_map_w(0, i, j) / density_map(i, j));
        momentum_map[1].set(i, j, momentum_map_w(1, i, j) / density_map(i, j));
        momentum_map[2].set(i, j, momentum_map_w(2, i, j) / density_map(i, j));
      }
}

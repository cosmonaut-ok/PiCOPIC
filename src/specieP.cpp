#include "specieP.hpp"

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
#if defined (PLASMA_SPATIAL_REGULAR) || defined (PLASMA_SPATIAL_CENTERED)
  LOG_DBG("Correcting area plasma particles macro amount to satisfy spatial distribution");
  double r_size = geometry->r_size;
  double z_size = geometry->z_size;
  double dr = geometry->r_cell_size;
  double dz = geometry->z_cell_size;

  if (geometry->walls[2]) r_size -= geometry->r_cell_size;
  if (geometry->walls[3]) z_size -= geometry->z_cell_size;

  macro_amount = lib::nearest_divide(p_macro_amount, r_size / dr * z_size / dz);
#else
  macro_amount = p_macro_amount;
#endif

  density[0] = p_left_density;
  density[1] = p_right_density;

  temperature = p_temperature;
}

SpecieP::~SpecieP()
{
  for (auto i = particles.begin(); i != particles.end(); i++)
    (**i).clear();
  particles.clear();
}

void SpecieP::fullyfill_spatial_distribution()
// ! apply flat spatial macroparticles distribution
// ! with gradient, declared for object
// ! gradient achived with macroparticles
// ! mass differentiation
{
  rectangular_spatial_distribution(geometry->bottom_r_grid_number,
                                   geometry->top_r_grid_number,
                                   geometry->left_z_grid_number,
                                   geometry->right_z_grid_number);
}

void SpecieP::linear_spatial_distribution(unsigned int int_cell_number,
                                          unsigned int ext_cell_number)
// ! spatial distribution for linear (for cartese coordinates)
// ! or tube-shaped macroparticles cloud
{
  rectangular_spatial_distribution(int_cell_number, ext_cell_number,
                                   0, geometry->z_size);
}

void SpecieP::rectangular_spatial_distribution(unsigned int int_cell_number,
                                               unsigned int ext_cell_number,
                                               unsigned int left_cell_number,
                                               unsigned int right_cell_number)
// ! spatial distribution for rectangular (for cartese coordinates)
// ! or cylindrical ring-shaped macroparticles cloud
{
#if defined (PLASMA_SPATIAL_CENTERED)
  rectangular_centered_placement (int_cell_number,
                                  ext_cell_number,
                                  left_cell_number,
                                  right_cell_number);
#elif defined (PLASMA_SPATIAL_FLAT) || defined (PLASMA_SPATIAL_RANDOM)
  rectangular_random_placement (int_cell_number,
                                ext_cell_number,
                                left_cell_number,
                                right_cell_number);
#elif PLASMA_SPATIAL_REGULAR
  rectangular_regular_placement (int_cell_number,
                                 ext_cell_number,
                                 left_cell_number,
                                 right_cell_number);
#endif

  set_mass_charges(int_cell_number, ext_cell_number, left_cell_number, right_cell_number);
}

void SpecieP::rectangular_random_placement (unsigned int int_cell_number,
                                            unsigned int ext_cell_number,
                                            unsigned int left_cell_number,
                                            unsigned int right_cell_number)
{
  //! macroparticles spatial placement
  //! for "random" and "flat" cases
  double dr = geometry->r_cell_size;
  double dz = geometry->z_cell_size;
  double dn = density[1] - density[0];
  double dl = density[0]; // dl is for "density left"

  double r_size = (ext_cell_number - int_cell_number) * dr;
  double z_size = (right_cell_number - left_cell_number) * dz;

  double rand_r, rand_z;

  // decrease size at r=r and z=z walls
  // this caused by particles formfactor
  if (geometry->walls[2]) r_size -= dr;
  if (geometry->walls[3]) z_size -= dz;

#ifdef PLASMA_SPATIAL_FLAT
  unsigned int macro_count = 0;
#endif

  for (unsigned int i = 0; i < macro_amount; i++)
  {
    vector<double> *n = new vector<double>(13, 0);

#if defined PLASMA_SPATIAL_RANDOM
    rand_r = math::random::uniform1();
    rand_z = math::random::uniform1();
#elif defined PLASMA_SPATIAL_FLAT
    rand_r = math::random::random_reverse(macro_count, 13);
    rand_z = math::random::random_reverse(macro_amount - 1 - macro_count, 11);
#endif
    if (rand_r == 0) rand_r = MNZL;
    if (rand_z == 0) rand_z = MNZL;

    P_POS_R((*n)) = r_size * rand_r;

    P_POS_PHI((*n)) = 0;

    if (density[0] == density[1])
      P_POS_Z((*n)) = z_size * rand_z;
    else
      P_POS_Z((*n)) = z_size / dn
        * (lib::sq_rt(pow(dl, 2) + rand_z * (2 * dl * dn + pow(dn, 2))) - dl);

    P_POS_R((*n)) += int_cell_number * dr;
    P_POS_R((*n)) += dr / 2.;
    P_POS_Z((*n)) += left_cell_number * dz; // shift by z to respect geometry with areas
    P_POS_Z((*n)) += dz / 2.;

#ifdef PLASMA_SPATIAL_FLAT
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

  double dr = geometry->r_cell_size;
  double dz = geometry->z_cell_size;
  // double dn = density[1] - density[0];
  // double dl = density[0]; // dl is for "density left"

  double r_size = (ext_cell_number - int_cell_number) * dr;
  double z_size = (right_cell_number - left_cell_number) * dz;

  // decrease size at r=r and z=z walls
  // this caused by particles formfactor
  if (geometry->walls[2]) r_size -= dr;
  if (geometry->walls[3]) z_size -= dz;

  // +dr and +dz because, we counting particles from "0" to "r_size".
  // So, there are "plus 1" interval, related to particles.
  double macro_distance = lib::sq_rt((r_size + dr) * (z_size + dz) / macro_amount);
  // count macroparticles by r and z
  double r_macro_count = 0;
  double z_macro_count = 0;

  unsigned int macro_count = 0;
  while (true)
  {
    double r_place = r_macro_count * macro_distance + MNZL;
    double z_place = z_macro_count * macro_distance + MNZL;

    if (z_place > z_size) // FIXME: some f*cking magic. Should be z_size - macro_distance
      break;

    if (r_place > r_size - macro_distance)
    {
      r_macro_count = 0;
      ++z_macro_count;
    }
    else
      ++r_macro_count;
    ++macro_count;

    vector<double> *n = new vector<double>(13, 0);

    P_POS_R((*n)) = r_place;
    P_POS_Z((*n)) = z_place;

    P_POS_R((*n)) += int_cell_number * dr;
    P_POS_R((*n)) += dr / 2.;
    P_POS_Z((*n)) += left_cell_number * dz; // shift by z to respect geometry with areas
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
  LOG_CRIT("rectangular_centered_placement still not implemented", 1);
}

void SpecieP::set_mass_charges (unsigned int int_cell_number,
                                unsigned int ext_cell_number,
                                unsigned int left_cell_number,
                                unsigned int right_cell_number)
{
  double dr = geometry->r_cell_size;
  double dz = geometry->z_cell_size;
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
  double top_r_num = geometry->top_r_grid_number;
  double bot_r_num = geometry->bottom_r_grid_number;
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
    P_CHARGE((**n)) = charge * n_per_macro;
    P_MASS((**n)) = mass * n_per_macro;
  }
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void SpecieP::velocity_distribution ()
{
#ifdef PLASMA_VELOCITY_THERMAL
  thermal_velocity_distribution();
#elif PLASMA_VELOCITY_RECTANGULAR
  rectangular_velocity_distribution();
#elif PLASMA_VELOCITY_EIGEN
  eigen_velocity_distribution();
#else
  LOG_CRIT("Plasma velocity distribution type unknown", 1);
#endif
}

void SpecieP::thermal_velocity_distribution ()
{
  // double therm_vel = lib::sq_rt(2. * EL_CHARGE * temperature / mass);

  // Sample the energies in the MJ distribution
  vector<double> energies = math::maxwell_juttner::maxwellJuttner(macro_amount, temperature);

  unsigned int macro_count = 0;

  // normalize to $\frac{1}{\sqrt{2}}$
  double norm = 0.707;

  for (auto p = particles.begin(); p != particles.end(); ++p)
  {
    double therm_vel_el = lib::sq_rt(2 * EL_CHARGE * energies[macro_count] / mass);

    double rnd_0 = math::random::uniform2();
    double rnd_1 = math::random::uniform2();
    double rnd_2 = math::random::uniform2();

    P_VEL_R((**p)) = rnd_0 * therm_vel_el * norm;
    P_VEL_PHI((**p)) = rnd_1 * therm_vel_el * norm;
    P_VEL_Z((**p)) = rnd_2 * therm_vel_el * norm;

    // take into account relativistic factor
    // maxwellJuttner procedure returns relativistic momentums
    P_VEL_R((**p)) *= lib::get_gamma_inv(P_VEL_R((**p)) * P_VEL_R((**p)));
    P_VEL_PHI((**p)) *= lib::get_gamma_inv(P_VEL_PHI((**p)) * P_VEL_PHI((**p)));
    P_VEL_Z((**p)) *= lib::get_gamma_inv(P_VEL_Z((**p)) * P_VEL_Z((**p)));

    ++macro_count;
  }
}

void SpecieP::rectangular_velocity_distribution ()
{
  // therm velocity for singular velocity component
  double therm_vel_cmp = lib::sq_rt(2. * EL_CHARGE * temperature / mass / 3.);

  for (auto p = particles.begin(); p != particles.end(); ++p)
  {
    double rnd_0, rnd_1, rnd_2;

    rnd_0 = math::random::uniform2();
    rnd_1 = math::random::uniform2();
    rnd_2 = math::random::uniform2();

    P_VEL_R((**p)) = rnd_0 * therm_vel_cmp;
    P_VEL_PHI((**p)) = rnd_1 * therm_vel_cmp;
    P_VEL_Z((**p)) = rnd_2 * therm_vel_cmp;

    // take into account relativistic factor
    // maxwellJuttner procedure returns relativistic momentums
    P_VEL_R((**p)) *= lib::get_gamma_inv(P_VEL_R((**p)) * P_VEL_R((**p)));
    P_VEL_PHI((**p)) *= lib::get_gamma_inv(P_VEL_PHI((**p)) * P_VEL_PHI((**p)));
    P_VEL_Z((**p)) *= lib::get_gamma_inv(P_VEL_Z((**p)) * P_VEL_Z((**p)));
  }
}

void SpecieP::eigen_velocity_distribution ()
// ! singular velocity distribution
{
  double therm_vel_cmp = lib::sq_rt(2. * EL_CHARGE * temperature / mass / 3.);
  double gamma = lib::get_gamma_inv(therm_vel_cmp);
  therm_vel_cmp *= gamma;

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
  double therm_vel = lib::sq_rt(2. * EL_CHARGE * temperature / mass);
  double gamma = lib::get_gamma_inv(therm_vel);
  therm_vel *= gamma;

  switch (dir)
  {
  case 0:
    for (auto p = particles.begin(); p != particles.end(); ++p)
      P_VEL_R((**p)) = therm_vel;
    break;
  case 1:
    for (auto p = particles.begin(); p != particles.end(); ++p)
      P_VEL_PHI((**p)) = therm_vel;
    break;
  case 2:
    for (auto p = particles.begin(); p != particles.end(); ++p)
      P_VEL_Z((**p)) = therm_vel;
    break;
  default:
    LOG_CRIT("Incorrect switch of rectangular directed velocity component: " << dir, 1);
    break;
  }
}

void SpecieP::inject () {}

void SpecieP::boris_pusher()
{
  // !
  // ! boris pusher
  // !
  for (auto p = particles.begin(); p != particles.end(); ++p)
  {
    // define vars directly in loop, because of multithreading
    double const1, const2, sq_velocity;
#ifdef PUSHER_BORIS_ADAPTIVE
    bool use_rel; // use relativistic calculations
#endif
#ifndef PUSHER_BORIS_CLASSIC // do not use gamma in non-relativistic boris pusher
    double gamma = 1;
#endif

    vector3d<double> e(0, 0, 0);
    vector3d<double> b(0, 0, 0);
    vector3d<double> velocity(P_VEL_R((**p)), P_VEL_PHI((**p)), P_VEL_Z((**p)));
    vector3d<double> vtmp(0, 0, 0);

    double pos_r = P_POS_R((**p));
    double pos_z = P_POS_Z((**p));

    // check if radius and longitude are correct
    if (isnan(pos_r) ||
        isinf(pos_r) != 0 ||
        isnan(pos_z) ||
        isinf(pos_z) != 0)
      LOG_CRIT("(boris_pusher): radius[" << pos_r
               << "] or longitude[" << pos_z
               << "] is not valid number. Can not continue.", 1);

    e = field_e->get_field(pos_r, pos_z);
    b = field_h->get_field(pos_r, pos_z);

    // ! \f$ const1 = \frac{q t}{2 m} \f$
    // ! where \f$ q, m \f$ - particle charge and mass, \f$ t = \frac{\Delta t_{step}}{2} \f$
    const1 = P_CHARGE((**p)) * time->step / (2 * P_MASS((**p)));

    e *= const1;
    b *= const1 * MAGN_CONST;

    // ! 0. check, if we should use classical calculations.
    // ! Required to increase modeling speed
#ifdef PUSHER_BORIS_ADAPTIVE
    if (pow(velocity[0], 2) + pow(velocity[1], 2) + pow(velocity[2], 2) > REL_LIMIT_POW_2)
      use_rel = true;
#endif
    // ! 1. Multiplication by relativistic factor (only for relativistic case)
    // ! \f$ u_{n-\frac{1}{2}} = \gamma_{n-\frac{1}{2}} * v_{n-\frac{1}{2}} \f$
#ifdef PUSHER_BORIS_ADAPTIVE
    if (use_rel)
#endif
#if defined (PUSHER_BORIS_RELATIVISTIC) || defined (PUSHER_BORIS_ADAPTIVE)
    {
      sq_velocity = velocity.length2();

      gamma = lib::get_gamma(sq_velocity);
      velocity *= gamma;
    }
#endif

    // ! 2. Half acceleration in the electric field
    // ! \f$ u'_n = u_{n-\frac{1}{2}} + \frac{q dt}{2 m E(n)} \f$
    // ! \f$ u'_n = u_{n-1 / 2} + \frac{q dt}{2 m E(n)} \f$
    velocity += e;

    // ! 3. Rotation in the magnetic field
    // ! \f$ u" = u' + \frac{2}{1 + B'^2}  [(u' + [u' \times B'(n)] ) \times B'(n)] \f$,
    // ! \f$ B'(n) = \frac{B(n) q dt}{2 m * \gamma_n} \f$
#ifdef PUSHER_BORIS_ADAPTIVE
    if (use_rel)
#endif
#if defined (PUSHER_BORIS_RELATIVISTIC) || defined (PUSHER_BORIS_ADAPTIVE)
    {
      sq_velocity = velocity.length2();
      gamma = lib::get_gamma_inv(sq_velocity);
      b *= gamma;
    }
#endif
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
#ifdef PUSHER_BORIS_ADAPTIVE
    if (use_rel)
#endif
#if defined (PUSHER_BORIS_RELATIVISTIC) || defined (PUSHER_BORIS_ADAPTIVE)
    {
      sq_velocity = velocity.length2();
      gamma = lib::get_gamma_inv(sq_velocity);
      velocity *= gamma;
    }
#endif

    P_VEL_R((**p)) = velocity[0];
    P_VEL_PHI((**p)) = velocity[1];
    P_VEL_Z((**p)) = velocity[2];
  }
}

void SpecieP::vay_pusher()
{
  // !
  // ! Vay pusher
  // !
  // LOG_CRIT("Vay pusher not implemented", 1)

  for (auto p = particles.begin(); p != particles.end(); ++p)
  {
    vector3d<double> velocity(P_VEL_R((**p)), P_VEL_PHI((**p)), P_VEL_Z((**p)));
    vector3d<double> uplocity(0, 0, 0); // u prime
    vector3d<double> psm(0, 0, 0); // pxsm, pysm, pzsm

    double charge = P_CHARGE((**p));
    double mass = P_MASS((**p));

    // MSG(P_CHARGE((**p)) << " " << charge);
    double charge_over_2mass_dt = charge * time->step / (2 * mass);

    double pos_r = P_POS_R((**p));
    double pos_z = P_POS_Z((**p));

    vector3d<double> e = field_e->get_field(pos_r, pos_z);
    vector3d<double> b = field_h->get_field(pos_r, pos_z);

    double gamma, sq_vel, s, us2, alpha, B2;

    // convert velocity to relativistic momentum
    sq_vel = velocity.length2();
    gamma = lib::get_gamma(sq_vel);
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
    gamma = lib::get_gamma_inv(sq_vel);

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
    gamma = lib::get_gamma_inv(sq_vel);
    psm *= gamma;

    P_VEL_R((**p)) = psm[0];
    P_VEL_PHI((**p)) = psm[1];
    P_VEL_Z((**p)) = psm[2];
  }
}



void SpecieP::hc_pusher()
{
  // !
  // ! Higuera-Cary pusher
  // !

  for (auto p = particles.begin(); p != particles.end(); ++p)
  {
    vector3d<double> velocity(P_VEL_R((**p)), P_VEL_PHI((**p)), P_VEL_Z((**p)));
    vector3d<double> uplocity(0, 0, 0); // u prime
    vector3d<double> psm(0, 0, 0); // pxsm, pysm, pzsm
    vector3d<double> um(0, 0, 0); // pxsm, pysm, pzsm
    vector3d<double> up(0, 0, 0); // pxsm, pysm, pzsm

    double charge = P_CHARGE((**p));
    double mass = P_MASS((**p));

    // MSG(P_CHARGE((**p)) << " " << charge);
    double charge_over_2mass_dt = charge * time->step / (2 * mass);

    double pos_r = P_POS_R((**p));
    double pos_z = P_POS_Z((**p));

    vector3d<double> e = field_e->get_field(pos_r, pos_z);
    vector3d<double> b = field_h->get_field(pos_r, pos_z);
    vector3d<double> b2(0, 0, 0);
    vector3d<double> b_cross(0, 0, 0);

    double gamma, sq_vel, s, us2, alpha, B2;

    // convert velocity to relativistic momentum
    sq_vel = velocity.length2();
    gamma = lib::get_gamma(sq_vel);
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
    gamma = 1. / sqrt( 0.5*( gfm2 - B2 +
                             sqrt( pow( gfm2 - B2, 2 )
                                   + 4.0 * ( B2 + pow( b[0]*um[0]
                                                       + b[1]*um[1]
                                                       + b[2]*um[2], 2 )
                                     )
                               )
                         )
      );

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
    gamma = lib::get_gamma_inv(sq_vel);
    psm *= gamma;

    P_VEL_R((**p)) = psm[0];
    P_VEL_PHI((**p)) = psm[1];
    P_VEL_Z((**p)) = psm[2];
  }
}



void SpecieP::half_step_mover_cylindrical()
{
  double half_dt = time->step / 2.;

  for (auto p = particles.begin(); p != particles.end(); ++p)
  {
    // check if radius and longitude are correct
    double pos_r = P_POS_R((**p));
    double pos_phi = P_POS_PHI((**p));
    double pos_z = P_POS_Z((**p));

    // check if radius and longitude are correct
    if (isnan(pos_r) ||
        isinf(pos_r) != 0 ||
        isnan(pos_z) ||
        isinf(pos_z) != 0)
      LOG_CRIT("(half_step_mover_cylindrical): radius[" << pos_r
               << "] or longitude[" << pos_z
               << "] is not valid number", 1);

    P_POS_R((**p)) = pos_r + P_VEL_R((**p)) * half_dt;
    //! we use "fake" rotation component to correct position from xy to rz pane
    P_POS_PHI((**p)) = pos_phi + P_VEL_PHI((**p)) * half_dt;
    P_POS_Z((**p)) = pos_z + P_VEL_Z((**p)) * half_dt;
  }
}

void SpecieP::reflect ()
{
  double dr = geometry->r_cell_size;
  double dz = geometry->z_cell_size;
  double half_dr = dr / 2.;
  double half_dz = dz / 2.;

  double radius_wall = geometry->r_size - dr / 2.;
  double longitude_wall = geometry->z_size - dz / 2.;

  double radius_wallX2 = radius_wall * 2.;
  double longitude_wallX2 = longitude_wall * 2.;

  // shift for converting local positions into global and back
  double r_shift = geometry->bottom_r_grid_number * dr;
  double z_shift = geometry->left_z_grid_number * dz;

  for (auto p = particles.begin(); p != particles.end(); ++p)
  {
    double pos_r = P_POS_R((**p)) - r_shift;
    double pos_z = P_POS_Z((**p)) - z_shift;

    if (pos_r > radius_wall && geometry->walls[2])
    {
      P_POS_R((**p)) = radius_wallX2 - pos_r + r_shift;
      P_VEL_R((**p)) = - P_VEL_R((**p));
    }

    if (pos_z > longitude_wall && geometry->walls[3])
    {
      P_POS_Z((**p)) = longitude_wallX2 - pos_z + z_shift;
      P_VEL_Z((**p)) = - P_VEL_Z((**p));
    }

    if (pos_r < half_dr && geometry->walls[0])
    {
      P_POS_R((**p)) = dr - pos_r + r_shift;
      P_VEL_R((**p)) = - P_VEL_R((**p));
    }

    if (pos_z < half_dz && geometry->walls[1])
    {
      P_POS_Z((**p)) = dr - pos_z + z_shift;
      P_VEL_Z((**p)) = - P_VEL_Z((**p));
    }
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
    double r = lib::sq_rt(pos_r * pos_r + pos_phi * pos_phi);
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
    P_COS((**p)) = lib::sq_rt(1 - sin * sin);

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

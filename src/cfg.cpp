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

#include "cfg.hpp"
#include <typeinfo>

#define MIN_DOMAIN_GRID_AMOUNT 64

Cfg::Cfg(const char *json_file_name)
{
  //! read given json file and parse it
  stringstream ss;
  ifstream f;
  // unsigned int i;

  time = new TimeSim();
  output_data = new save_data();

  // Read Json file
  f.open(json_file_name, ios::binary);
  if (!f.is_open()) exit(1);
  ss << f.rdbuf();
  f.close();

  //! Parse Json data
  value v;
  ss >> v;
  string err = get_last_error();
  if(!err.empty())
    LOG_S(FATAL) << err;

  json_data = v;

  //! set logfile
  log_file = (char*)json_data.get<object>()["log_file"].get<string>().c_str();
  // set log level, file etc
#ifdef DEBUG
  loguru::add_file(log_file, loguru::Truncate, loguru::Verbosity_MAX);
  loguru::g_stderr_verbosity = loguru::Verbosity_MAX;
#else
  loguru::add_file(log_file, loguru::Truncate, loguru::Verbosity_INFO);
  loguru::g_stderr_verbosity = loguru::Verbosity_INFO;
#endif

  //! update json with additional items
  //! to inform about build options
  //! and package version
  // set package version to config metadata
  // json_data.get<object>()["package_version"] = value((string)PACKAGE_VERSION);
  object o;
  o["package_version"] = value((string)PACKAGE_VERSION);
#ifdef DEBUG
  o["debug"] = value(true);
#else
  o["debug"] = value(false);
#endif

#ifdef SPEEDUP
  o["ieee"] = value(false);
#else
  o["ieee"] = value(true);
#endif

#if defined(PLASMA_SPATIAL_CENTERED)
  o["plasma_spatial_distribution"] = value("centered");
#elif defined(PLASMA_SPATIAL_FLAT)
  o["plasma_spatial_distribution"] = value("flat");
#elif defined(PLASMA_SPATIAL_RANDOM)
  o["plasma_spatial_distribution"] = value("random");
#elif defined(PLASMA_SPATIAL_REGULAR)
  o["plasma_spatial_distribution"] = value("regular");
#endif

#if defined(PLASMA_VELOCITY_THERMAL)
  o["plasma_velocity_distribution"] = value("thermal");
#elif defined(PLASMA_VELOCITY_EIGEN)
  o["plasma_velocity_distribution"] = value("eigen");
#elif defined(PLASMA_VELOCITY_RECTANGULAR)
  o["plasma_velocity_distribution"] = value("rectangular");
#endif

#if defined(PUSHER_BORIS_ADAPTIVE)
  o["particles_pusher"] = value("boris adaptive");
#elif defined(PUSHER_BORIS_CLASSIC)
  o["particles_pusher"] = value("boris classic");
#elif defined(PUSHER_BORIS_RELATIVISTIC)
  o["particles_pusher"] = value("boris relativistic");
#elif defined(PUSHER_HIGUERA_CARY)
  o["particles_pusher"] = value("higuera-cary");
  #elif defined(PUSHER_VAY)
  o["particles_pusher"] = value("vay");
#endif

#if defined(CCS_ZIGZAG)
  o["charge_conservation"] = value("zigzag");
#elif defined(CCS_VILLASENOR_BUNEMAN)
  o["charge_conservation"] = value("villasenor-buneman");
#endif

#if defined(TEMP_CALC_COUNTING)
  o["temperature_calculation_algorithm"] = value("counting");
#elif defined(TEMP_CALC_WEIGHTING)
  o["temperature_calculation_algorithm"] = value("weighting");
#endif

#if defined(DENSITY_CALC_COUNTING)
  o["density_calculation_algorithm"] = value("counting");
#elif defined(DENSITY_CALC_WEIGHTING)
  o["density_calculation_algorithm"] = value("weighting");
#endif

  o["build_flags"] = value(CXXFLAGS);

  json_data.get<object>()["build_options"] = value(o);

  //! set amount of macroparticles
  macro_amount = json_data.get<object>()["macro_amount"].get<double>();

  //! make beam macroparticles X times smaller, than
  //! macroparticles of background plasma particle species
  //! required for short beams, which has really small volume
  //! related to whole simulation volume
  beam_macro_ratio = (unsigned int)json_data.get<object>()["beam_plasma_macro_size_ratio"].get<double>();

  //! initialize geometry, set PML
  init_geometry();

  //! initialize time
  init_time();

  //! initialize particle parameters
  init_particles();

  //! initialize particles beam parameters
  init_beam();

  //! weight macroparticles alignment
  //! for each particles specie and beam
  weight_macro_amount();

  init_output_data();
  init_probes();

  method_limitations_check();
}

Cfg::~Cfg()
{
}

void Cfg::init_particles()
{
  //! initialize particles charge, mass, number,
  //! density and temperature for all particle species

  object& json_root = json_data.get<object>();
  picojson::array& particle_species_list = json_root["particles"].get<picojson::array>();

  for (auto i=particle_species_list.begin(); i != particle_species_list.end(); ++i)
  {
    particle_specie p_s;

    object& o = i->get<object>();

    p_s.name = (char*)o["name"].get<string>().c_str();
    p_s.mass = (int)o["mass"].get<double>();
    p_s.charge = (int)o["charge"].get<double>();
    // p_s.left_density = o["density"].get<object>()["left"].get<double>();
    // p_s.right_density = o["density"].get<object>()["right"].get<double>();
    p_s.left_density = o["density"].get<double>(); // TODO: density profile still not supported, but planned
    p_s.right_density = o["density"].get<double>(); // TODO: density profile still not supported, but planned
    p_s.temperature = o["temperature"].get<double>();

    particle_species.push_back(p_s);
  }
}

void Cfg::init_probes ()
{
  //! initialize particles charge, mass, number,
  //! density and temperature for all particle species

  object& json_root = json_data.get<object>()["data"].get<object>();
  picojson::array& probes_list = json_root["probes"].get<picojson::array>();

  for (auto i=probes_list.begin(); i != probes_list.end(); ++i)
  {
    probe p_p;

    object& o = i->get<object>();

    p_p.component = o["component"].get<string>().c_str();

    try { p_p.specie = o["specie"].get<string>().c_str(); }
    catch (std::exception& e) { p_p.specie = "unknown"; }

     string p_shape = o["shape"].get<string>().c_str();

    if ( p_shape.compare("rec") == 0 )
    {
      p_p.shape = 0;
      p_p.r_start = (int)o["start"].get<object>()["r"].get<double>();
      p_p.r_end = (int)o["end"].get<object>()["r"].get<double>();
      p_p.z_start = (int)o["start"].get<object>()["z"].get<double>();
      p_p.z_end = (int)o["end"].get<object>()["z"].get<double>();
    }
    else if ( p_shape.compare("col") == 0 )
    {
      p_p.shape = 1;
      p_p.r_start = -1;
      p_p.r_end = -1;
      p_p.z_end = (int)o["z"].get<double>();
      p_p.z_start = -1;
    }
    else if ( p_shape.compare("row") == 0 )
    {
      p_p.shape = 2;
      p_p.r_end = (int)o["r"].get<double>();
      p_p.r_start = -1;
      p_p.z_start = -1;
      p_p.z_end = -1;
    }
    else if ( p_shape.compare("dot") == 0 )
    {
      p_p.shape = 3;
      p_p.r_end = (int)o["r"].get<double>();
      p_p.r_start = -1;
      p_p.z_end = (int)o["z"].get<double>();
      p_p.z_start = -1;
    }
    else
      LOG_S(FATAL) << "Unknown probe shape ``" << p_shape << "''";

    p_p.schedule = (int)o["schedule"].get<double>();

    if (p_p.r_start > p_p.r_end || p_p.z_start > p_p.z_end)
    {
      LOG_S(FATAL) << "Incorrect probe's " << p_p.component << "/" << p_shape << " shape: ["
               << p_p.r_start << "," << p_p.r_end << ","
               << p_p.z_start << "," << p_p.z_end << "]";
    }
    else if (p_p.r_end > geometry->r_grid_amount)
    {
      LOG_S(FATAL) << "Probe's " << p_p.component << "/" << p_shape << " radius is out of simulation domain: "
               << p_p.r_end << ". Must be less, than "
               << geometry->r_grid_amount;
    }
    else if (p_p.z_end > geometry->z_grid_amount)
    {
      LOG_S(FATAL) << "Probe's " << p_p.component << "/" << p_shape << " longitude is out of simulation domain: "
               << p_p.z_end << ". Must be less, than "
		   << geometry->z_grid_amount;
    }
    else
      probes.push_back(p_p);
  }
}

void Cfg::init_beam()
{
  //! initialize particles bunch data
  //! particles bunch should be injected to plasma

  object& json_root = json_data.get<object>();
  picojson::array& particle_beams_list = json_root["beams"].get<picojson::array>();

  if (! particle_beams_list.empty())
    for (auto i=particle_beams_list.begin(); i != particle_beams_list.end(); ++i)
      {
	particle_beam p_b;

	object& o = i->get<object>();
	object& b_o = o["bunch"].get<object>();

	p_b.name = (char*)o["name"].get<string>().c_str();
	p_b.mass = (int)o["mass"].get<double>();
	p_b.charge = (int)o["charge"].get<double>();
	p_b.start_time = o["start_time"].get<double>();
	p_b.velocity = o["velocity"].get<double>();
	p_b.bunches_amount = b_o["amount"].get<double>();
	p_b.bunches_distance = b_o["distance"].get<double>();
	p_b.bunch_length = b_o["length"].get<double>();
	p_b.bunch_radius = b_o["radius"].get<double>();
	p_b.density = b_o["density"].get<double>();
	// p_b.macro_amount = b_o["macro_amount"].get<double>();

	p_b.current_bunch_number = 0;

	particle_beams.push_back(p_b);
      }
}

// void Cfg::init_boundary ()
// {
//   //! initialize boundaries and boundary conditions
//   XMLElement *sub_root = try_first_child(xml_data, "boundary_maxwell_conditions");
//   //
//   boundary_maxwell_e_phi_upper = atof(try_first_child(sub_root, "e_fi_upper")->GetText());
//   boundary_maxwell_e_phi_left = atof(try_first_child(sub_root, "e_fi_left")->GetText());
//   boundary_maxwell_e_phi_right = atof(try_first_child(sub_root, "e_fi_right")->GetText());
//   boundary_conditions = atoi(try_first_child(xml_data, "boundary_conditions")->GetText());
// }

void Cfg::init_geometry ()
{
  //! initialize geometry
  object& json_root = json_data.get<object>()["geometry"].get<object>();
  object& json_root_size = json_root["size"].get<object>();
  object& json_root_grid = json_root["grid"].get<object>();

  double radius = json_root_size["radius"].get<double>();
  double longitude = json_root_size["longitude"].get<double>();

  int n_grid_r = (int)json_root_grid["radius"].get<double>();
  int n_grid_z = (int)json_root_grid["longitude"].get<double>();

  // init PML
#ifdef USE_PML
  object& json_root_pml = json_root["pml"].get<object>();
  double left_wall = json_root_pml["left_wall"].get<double>();
  double right_wall = json_root_pml["right_wall"].get<double>();
  double outer_wall = json_root_pml["outer_wall"].get<double>();
  double sigma_1 = json_root_pml["sigma_1"].get<double>();
  double sigma_2 = json_root_pml["sigma_2"].get<double>();

  geometry = new Geometry(radius, longitude, 0, n_grid_r, 0, n_grid_z,
                          left_wall, right_wall, outer_wall, sigma_1, sigma_2,
                          true, true, true, true); // init all walls
  #else
  geometry = new Geometry(radius, longitude, 0, n_grid_r, 0, n_grid_z,
                          true, true, true, true); // init all walls
#endif

  //
#ifdef SINGLETHREAD
  LOG_S(INFO) << "Singlethread mode. Reducing domains amount to 1";
  geometry->domains_by_r = 1;
  geometry->domains_by_z = 1;
#else
  geometry->domains_by_r = (int)json_root["domains_amount"].get<object>()["radius"].get<double>();
  geometry->domains_by_z = (int)json_root["domains_amount"].get<object>()["longitude"].get<double>();
#endif
}

void Cfg::init_time ()
{
  //! initialize time parameters
  object& json_root = json_data.get<object>()["time"].get<object>();

  time->start = json_root["start"].get<double>();
  time->end = json_root["end"].get<double>();
  time->step = json_root["step"].get<double>();
  time->current = 0;
}

void Cfg::init_output_data()
{
  //! initialize file saving parameters, like path to computed data files,
  //! path to system state data files max frames number, placed to one file etc.
  object& json_root = json_data.get<object>()["data"].get<object>();

  output_data->data_root = (char*)json_root["data_root"].get<string>().c_str();

  output_data->compress = json_root["compression"].get<object>()["use"].get<bool>();
  output_data->compress_level = (int)json_root["compression"].get<object>()["level"].get<double>();
}

void Cfg::weight_macro_amount()
//! calculate alignment of macroparticles amount
//! to each particles specie and each beam
{
  double norm_sum = 0;

// increase normalization sum for particle species (assume, that specie size times species amount)
  norm_sum += geometry->r_grid_amount * geometry->z_grid_amount * particle_species.size();

// increase normalization sum for particle beams
  for (auto b = particle_beams.begin(); b != particle_beams.end(); ++b)
  {
    double beam_r_grid_size = (*b).bunch_radius / geometry->r_cell_size;
    double beam_z_grid_size = (*b).bunch_length * (*b).bunches_amount / geometry->z_cell_size;

    norm_sum += beam_r_grid_size * beam_z_grid_size * beam_macro_ratio;
  }

  // calculate normalization coeffitient
  double norm = macro_amount / norm_sum;

  // message about pusher, used in system
#ifdef PUSHER_BORIS_CLASSIC
  LOG_S(INFO) << "Using classic Boris particles pusher";
#elif defined PUSHER_BORIS_ADAPTIVE
  LOG_S(INFO) << "Using adaptive Boris particles pusher";
#elif defined PUSHER_BORIS_RELATIVISTIC
  LOG_S(INFO) << "Using fully relativistic Boris particles pusher";
#elif defined PUSHER_VAY
  LOG_S(INFO) << "Using Vay particles pusher [10.1063/1.2837054]";
#elif defined PUSHER_HIGUERA_CARY
  LOG_S(INFO) << "Using Higuera-Cary particles pusher [10.1063/1.4979989]";
#else
  LOG_S(FATAL) << "Undefined particles pusher used";
#endif

  // message about charge conservation scheme, used in system
#if defined(CCS_ZIGZAG)
  LOG_S(INFO) << "Using ZigZag charge conservation scheme [10.1016/S0010-4655(03)00437-5]";
#elif defined(CCS_VILLASENOR_BUNEMAN)
  LOG_S(INFO) << "Using Villasenor-Buneman charge conservation scheme [10.1016/0010-4655(92)90169-Y]";
#endif

  // message about coulomb colisions
#ifdef COLLISIONS
#ifdef COULOMB_COLLISIONS_TA77S
  LOG_S(INFO) << "Using Takizuka and Abe Coulomb collisions scheme";
  LOG_S(INFO) << "\twith symmetric weighted particles correction [10.1002/ctpp.201700121]";
#elif COULOMB_COLLISIONS_SK98
  LOG_S(INFO) << "Using Sentoku and Kemp Coulomb collisions scheme [10.1143/JPSJ.67.4084, 10.1016/j.jcp.2008.03.043]";
#elif COULOMB_COLLISIONS_P12
  LOG_S(INFO) << "Using Perez et al. collisions scheme [10.1063/1.4742167]";
#endif
#endif

  // align macroparticles amount to particle species
  // with respect to normalization
#if defined (PLASMA_SPATIAL_REGULAR) || defined (PLASMA_SPATIAL_CENTERED)
  LOG_S(INFO) << "Using amount-dependent plama macroparticles spatial distribution";
  LOG_S(INFO) << "Amount of plasma macroparticles could be changed to satisfy spatial distribution requirement";
#endif
  for (auto p_s = particle_species.begin(); p_s != particle_species.end(); ++p_s)
  {
    (*p_s).macro_amount = (unsigned int)(geometry->r_grid_amount * geometry->z_grid_amount * norm);
    LOG_S(INFO) << "Macro amount for ``" << (*p_s).name << "'' speice is " << (*p_s).macro_amount;
  }

  // align macroparticles amount to particle beams
  // with respect to normalization
  for (auto b = particle_beams.begin(); b != particle_beams.end(); ++b)
  {
    double beam_r_grid_size = (*b).bunch_radius / geometry->r_cell_size;
    double beam_z_grid_size = (*b).bunch_length * (*b).bunches_amount / geometry->z_cell_size;

    (*b).macro_amount = (unsigned int)(beam_r_grid_size * beam_z_grid_size * norm * beam_macro_ratio);
    LOG_S(INFO) << "Macro amount for ``" << (*b).name << "'' beam is " << (*b).macro_amount;
  }
}

bool Cfg::method_limitations_check ()
{
  double full_macro_amount = 0;
  double electron_density = 0;
  double ion_density = 0;
  double electron_temperature = 0;
  double ion_temperature = 0;

  // grid limitations
  if ((fmod(geometry->domains_by_r, 2) != 0 && geometry->domains_by_r != 1) || (fmod(geometry->domains_by_z, 2) != 0 && geometry->domains_by_z != 1))
    {
      LOG_S(FATAL) << "Amount of domains should be multiple of 2, or 1";
    }

  // if (geometry->r_grid_amount / geometry->domains_by_r < MIN_DOMAIN_GRID_AMOUNT)
  //   {
  //     LOG_S(FATAL) << "Too small domain size by r. Must be ``" << MIN_DOMAIN_GRID_AMOUNT << "'' cells or more", 1);
  //   }

  // if (geometry->z_grid_amount / geometry->domains_by_z < MIN_DOMAIN_GRID_AMOUNT)
  //   {
  //     LOG_S(FATAL) << "Too small domain size by z. Must be ``" << MIN_DOMAIN_GRID_AMOUNT << "'' cells or more", 1);
  //   }

  // if (geometry->r_grid_amount / geometry->domains_by_r * geometry->z_grid_amount / geometry->domains_by_z < MIN_DOMAIN_GRID_AMOUNT * MIN_DOMAIN_GRID_AMOUNT)
  //   {
  //     LOG_S(FATAL) << "Too small domain size. Must be ``" << MIN_DOMAIN_GRID_AMOUNT * MIN_DOMAIN_GRID_AMOUNT << "'' cells or more", 1);
  //   }

// try to figure out, where is electrons
  for (auto i = particle_species.begin(); i != particle_species.end(); ++i)
  {
    full_macro_amount += (*i).macro_amount; // calculate total number of macroparticles

    if (strcmp((*i).name, "Electrons") == 0
        || strcmp((*i).name, "electrons") == 0
        || strcmp((*i).name, "Electron") == 0
        || strcmp((*i).name, "electron") == 0
        || (*i).mass == 1)
    {
      electron_density = ((*i).left_density + (*i).right_density) / 2;
      electron_temperature = (*i).temperature;
    }

    if (strcmp((*i).name, "Ions") == 0
        || strcmp((*i).name, "ions") == 0
        || strcmp((*i).name, "Ion") == 0
        || strcmp((*i).name, "ion") == 0)
    {
      ion_density = ((*i).left_density + (*i).right_density) / 2;
      ion_temperature = (*i).temperature;
    }
  }

  if (electron_density == 0)
  {
    LOG_S(FATAL) << "There is no electrons present in system";
  }

  double plasma_freq = phys::plasma::plasma_frequency ( electron_density );

  // TODO: WTF
  unsigned int time_multiplicator = 100;

  // particle should not fly more, than 1 cell per time step
  double max_time_c = min(geometry->r_cell_size / constant::LIGHT_VEL,
                          geometry->z_cell_size / constant::LIGHT_VEL);
  double max_time_wp = 2 / plasma_freq / time_multiplicator;
  double max_time = min(max_time_c, max_time_wp);

  if (time->step > max_time)
  {
    LOG_S(FATAL) << "Too large time step: ``"
             << time->step
             << " s.''. Should be less, than ``"
             << max_time
             << " s.''";
  }

  // cell size d{r,z} < \lambda debye

  double debye_length = phys::plasma::debye_length (electron_density, ion_density, ion_temperature, electron_temperature);
  unsigned int debye_multiplicator = 100;

  if ( electron_temperature != 0)
  {
    if (geometry->r_cell_size > debye_length * debye_multiplicator
	|| geometry->z_cell_size > debye_length * debye_multiplicator)
    {
      LOG_S(FATAL) << "Too large grid size: ``"
	       << geometry->r_cell_size << " x " << geometry->z_cell_size
	       << " m.''. Should be less, than ``"
	       << debye_length * debye_multiplicator
		   << " m.''";
    }

    // \f$ L >> R_{debye} \f$
    // if (geometry->r_size < debye_length * debye_multiplicator
    // 	|| geometry->z_size < debye_length * debye_multiplicator)
    // {
    //   LOG_S(FATAL) << "Too small system size: ``"
    // 	       << geometry->r_size << " x " << geometry->z_size
    // 	       << " m.''. Should be more, than ``"
    // 	       <<  debye_length * debye_multiplicator
    // 		   << " m.''";
    // }

    // \f$ L << N_{particles} * R_{debye} \f$
    if (geometry->r_size > debye_length * full_macro_amount * debye_multiplicator
	|| geometry->z_size > debye_length * full_macro_amount * debye_multiplicator)
    {
      LOG_S(FATAL) << "Too large system size: ``"
	       << geometry->r_size << " x " << geometry->z_size
	       << " m.''. Should be less, than ``"
	       <<  debye_length * full_macro_amount * debye_multiplicator
	       << " m.''";
    }
  }

  unsigned int grid_multiplicator = 10;

  if (geometry->r_grid_amount * geometry->r_grid_amount * grid_multiplicator > full_macro_amount)
  {
    LOG_S(WARNING) << "Too small summary number of macroparticles: ``"
             << full_macro_amount
             << "''. Should be more, than ``"
             <<  geometry->r_grid_amount * geometry->r_grid_amount * grid_multiplicator
             << "''. You could get not relevant results";
  }

  return true;
}

string Cfg::cfg2str()
{
  // make some generic preparations
  value json_data_for_dump = json_data;
  json_data_for_dump.get<object>()["data"].get<object>()["data_root"].set<string>(".");

#ifdef DEBUG
  string str = json_data_for_dump.serialize(true);
#else
  string str = json_data_for_dump.serialize();
#endif // DEBUG
  return str;
}

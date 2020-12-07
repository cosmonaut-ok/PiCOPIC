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

#include "cfg.hpp"

using namespace std;
using namespace picojson;

Cfg::Cfg(const std::string json_file_name)
{
  //! read given json file and parse it
  stringstream ss;
  ifstream f;
  // unsigned int i;

  time = new TimeSim();
  output_data = new save_data();

  // Read Json file
  f.open(json_file_name.c_str(), ios::binary);
  if (!f.is_open())
    LOG_S(FATAL) << "Can not read configuration file ``" << json_file_name.c_str() << "''";
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
#ifdef ENABLE_DEBUG
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
#ifdef ENABLE_DEBUG
  o["debug"] = value(true);
#else
  o["debug"] = value(false);
#endif

#ifndef ENABLE_IEEE
  o["ieee"] = value(false);
#else
  o["ieee"] = value(true);
#endif

#ifdef SWITCH_PLASMA_SPATIAL_CENTERED
  o["plasma_spatial_distribution"] = value("centered");
#elif defined(SWITCH_PLASMA_SPATIA_FLAT)
  o["plasma_spatial_distribution"] = value("flat");
#elif defined(SWITCH_PLASMA_SPATIAL_RANDOM)
  o["plasma_spatial_distribution"] = value("random");
#elif defined(SWITCH_PLASMA_SPATIAL_REGULAR)
  o["plasma_spatial_distribution"] = value("regular");
#endif

#ifdef SWITCH_PLASMA_VELOCITY_THERMAL
  o["plasma_velocity_distribution"] = value("thermal");
#elif defined(SWITCH_PLASMA_VELOCITY_EIGEN)
  o["plasma_velocity_distribution"] = value("eigen");
#elif defined(SWITCH_PLASMA_VELOCITY_RECTANGULAR)
  o["plasma_velocity_distribution"] = value("rectangular");
#endif

#ifdef SWITCH_PUSHER_BORIS_ADAPTIVE
  o["particles_pusher"] = value("boris adaptive");
#elif defined(SWITCH_PUSHER_BORIS)
  o["particles_pusher"] = value("boris classic");
#elif defined(SWITCH_PUSHER_BORIS_RELATIVISTIC)
  o["particles_pusher"] = value("boris relativistic");
#elif defined(SWITCH_PUSHER_HC)
  o["particles_pusher"] = value("higuera-cary");
#elif defined(SWITCH_PUSHER_VAY)
  o["particles_pusher"] = value("vay");
#endif

#ifdef SWITCH_CCS_ZIGZAG
  o["current_deposition"] = value("zigzag");
#elif defined(SWITCH_CCS_VB)
  o["current_deposition"] = value("villasenor-buneman");
#endif

#ifdef SWITCH_TEMP_CALC_COUNTING
  o["temperature_calculation_algorithm"] = value("counting");
#elif defined(SWITCH_TEMP_CALC_WEIGHTING)
  o["temperature_calculation_algorithm"] = value("weighting");
#endif

#ifdef SWITCH_DENSITY_CALC_COUNTING
  o["density_calculation_algorithm"] = value("counting");
#elif defined(SWITCH_DENSITY_CALC_WEIGHTING)
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

    string specie_name = o["name"].get<string>();

    unsigned int current_specie_id = algo::common::hash_from_string(specie_name, SPECIE_HASH_SALT);

    p_s.id = current_specie_id;
    p_s.name = specie_name;
    p_s.mass = (int)o["mass"].get<double>();
    p_s.charge = (int)o["charge"].get<double>();
    // p_s.left_density = o["density"].get<object>()["left"].get<double>();
    // p_s.right_density = o["density"].get<object>()["right"].get<double>();
    p_s.left_density = o["density"].get<double>(); // TODO: density profile still not supported, but planned
    p_s.right_density = o["density"].get<double>(); // TODO: density profile still not supported, but planned
    p_s.temperature = o["temperature"].get<double>();

    particle_species.push_back(p_s);

    ++current_specie_id;
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

    //// generate probe path
    p_p.path = p_p.component + PATH_DELIMITER;

    // add path part for specie-based probes
    if ( p_p.component.compare("temperature") == 0
         || p_p.component.compare("density") == 0 )
      p_p.path += p_p.specie + PATH_DELIMITER;

    if ( p_shape.compare("rec") == 0 )
    {
      p_p.shape = 0;
      p_p.dims.push_back((size_t)o["start"].get<object>()["r"].get<double>());
      p_p.dims.push_back((size_t)o["start"].get<object>()["z"].get<double>());
      p_p.dims.push_back((size_t)o["end"].get<object>()["r"].get<double>());
      p_p.dims.push_back((size_t)o["end"].get<object>()["z"].get<double>());

      // p_p.dims[2] = (int)o["end"].get<object>()["r"].get<double>();
      // p_p.dims[1] = (int)o["start"].get<object>()["z"].get<double>();
      // p_p.dims[3] = (int)o["end"].get<object>()["z"].get<double>();

      // add shape and size to probe path
      p_p.path += "rec";
      p_p.path += PATH_DELIMITER;
      p_p.path += to_string(p_p.dims[0]);
      p_p.path += RANGE_DELIMITER;
      p_p.path += to_string(p_p.dims[2]);
      p_p.path += SPACE_DELIMITER;
      p_p.path += to_string(p_p.dims[1]);
      p_p.path += RANGE_DELIMITER;
      p_p.path += to_string(p_p.dims[3]);
    }
    else if ( p_shape.compare("col") == 0 )
    {
      p_p.shape = 1;
      p_p.dims[0] = 0;
      p_p.dims[2] = 0;
      p_p.dims[3] = (int)o["z"].get<double>();
      p_p.dims[1] = 0;

      // add shape and size to probe path
      p_p.path += "col";
      p_p.path += PATH_DELIMITER;
      p_p.path += to_string(p_p.dims[3]);
    }
    else if ( p_shape.compare("row") == 0 )
    {
      p_p.shape = 2;
      p_p.dims[2] = (int)o["r"].get<double>();
      p_p.dims[0] = 0;
      p_p.dims[1] = 0;
      p_p.dims[3] = 0;

      // add shape and size to probe path
      p_p.path += "row";
      p_p.path += PATH_DELIMITER;
      p_p.path += to_string(p_p.dims[2]);
    }
    else if ( p_shape.compare("dot") == 0 )
    {
      p_p.shape = 3;
      p_p.dims[2] = (int)o["r"].get<double>();
      p_p.dims[0] = 0;
      p_p.dims[3] = (int)o["z"].get<double>();
      p_p.dims[1] = 0;

      // add shape and size to probe path
      p_p.path += "dot";
      p_p.path += PATH_DELIMITER;
      p_p.path += to_string(p_p.dims[2]);
      p_p.path += SPACE_DELIMITER;
      p_p.path += to_string(p_p.dims[3]);
    }
    else
      LOG_S(FATAL) << "Unknown probe shape ``" << p_shape << "''";

    p_p.schedule = (int)o["schedule"].get<double>();

    if (p_p.dims[0] > p_p.dims[2] || p_p.dims[1] > p_p.dims[3])
    {
      LOG_S(FATAL) << "Incorrect probe's " << p_p.component << "/" << p_shape << " shape: ["
                   << p_p.dims[0] << "," << p_p.dims[2] << ","
                   << p_p.dims[1] << "," << p_p.dims[3] << "]";
    }
    else if (p_p.dims[2] > geometry->cell_amount[0])
    {
      LOG_S(FATAL) << "Probe's " << p_p.component << "/" << p_shape << " radius is out of simulation domain: "
                   << p_p.dims[2] << ". Must be less, than "
                   << geometry->cell_amount[0];
    }
    else if (p_p.dims[3] > geometry->cell_amount[1])
    {
      LOG_S(FATAL) << "Probe's " << p_p.component << "/" << p_shape << " longitude is out of simulation domain: "
                   << p_p.dims[3] << ". Must be less, than "
                   << geometry->cell_amount[1];
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

      string beam_name = "beam_" + o["name"].get<string>();

      // ! calculate current beam ID from its name:
      unsigned int current_beam_id = algo::common::hash_from_string(beam_name, BEAM_HASH_SALT);

      p_b.id = current_beam_id;
      p_b.name = beam_name;
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

      // ++current_beam_id;
    }
  else
    LOG_S(INFO) << "There is no particle beams to inject";
}

void Cfg::init_geometry ()
{
  //! initialize geometry
  object& json_root = json_data.get<object>()["geometry"].get<object>();
  object& json_root_size = json_root["size"].get<object>();
  object& json_root_grid = json_root["grid"].get<object>();

  double radius = json_root_size["radius"].get<double>();
  double longitude = json_root_size["longitude"].get<double>();

  size_t n_grid_r = (size_t)json_root_grid["radius"].get<double>();
  size_t n_grid_z = (size_t)json_root_grid["longitude"].get<double>();

  vector<double> _size = {radius, longitude};
  vector<size_t> _dims = {0, 0, n_grid_r, n_grid_z};
  vector<bool> _walls = {true, true, true, true};

  // init PML
#ifdef ENABLE_PML
  object& json_root_pml = json_root["pml"].get<object>();
  double left_wall = json_root_pml["left_wall"].get<double>();
  double right_wall = json_root_pml["right_wall"].get<double>();
  double outer_wall = json_root_pml["outer_wall"].get<double>();
  double sigma_1 = json_root_pml["sigma_1"].get<double>();
  double sigma_2 = json_root_pml["sigma_2"].get<double>();

  vector<double> _pml = {0, left_wall, right_wall, outer_wall};
  vector<double> _sigma = {sigma_1, sigma_2};

  geometry = new Geometry ( _size, _dims, _pml, _sigma, _walls);
  #else
  geometry = new Geometry ( _size, _dims, _walls);
#endif
  vector<size_t> domains (2, 0);
  geometry->domains_amount = domains;

  //
#ifdef ENABLE_SINGLETHREAD
  LOG_S(INFO) << "Singlethread mode. Reducing domains amount to 1";
  geometry->domains_amount[0] = 1;
  geometry->domains_amount[1] = 1;
#else
  geometry->domains_amount[0] = (int)json_root["domains_amount"].get<object>()["radius"].get<double>();
  geometry->domains_amount[1] = (int)json_root["domains_amount"].get<object>()["longitude"].get<double>();
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
  norm_sum += geometry->cell_amount[0] * geometry->cell_amount[1] * particle_species.size();

// increase normalization sum for particle beams
  for (auto b = particle_beams.begin(); b != particle_beams.end(); ++b)
  {
    double beam_r_grid_size = (*b).bunch_radius / geometry->cell_size[0];
    double beam_z_grid_size = (*b).bunch_length * (*b).bunches_amount / geometry->cell_size[1];

    norm_sum += beam_r_grid_size * beam_z_grid_size * beam_macro_ratio;
  }

  // calculate normalization coeffitient
  double norm = macro_amount / norm_sum;

  // message about pusher, used in system
#ifdef SWITCH_PUSHER_BORIS
  LOG_S(INFO) << "Using classic Boris particles pusher";
#elif defined(SWITCH_PUSHER_BORIS_ADAPTIVE)
  LOG_S(INFO) << "Using adaptive Boris particles pusher";
#elif defined(SWITCH_PUSHER_BORIS_RELATIVISTIC)
  LOG_S(INFO) << "Using fully relativistic Boris particles pusher";
#elif defined(SWITCH_PUSHER_VAY)
  LOG_S(INFO) << "Using Vay particles pusher [10.1063/1.2837054]";
#elif defined(SWITCH_PUSHER_HC)
  LOG_S(INFO) << "Using Higuera-Cary particles pusher [10.1063/1.4979989]";
#else
  LOG_S(FATAL) << "Undefined particles pusher used";
#endif

  // message about charge conservation scheme, used in system
#ifdef SWITCH_CCS_ZIGZAG
  LOG_S(INFO) << "Using ZigZag charge conservation scheme [10.1016/S0010-4655(03)00437-5]";
#elif defined(SWITCH_CCS_VB)
  LOG_S(INFO) << "Using Villasenor-Buneman charge conservation scheme [10.1016/0010-4655(92)90169-Y]";
#endif

  // message about coulomb colisions
#ifdef ENABLE_COULOMB_COLLISIONS
  LOG_S(WARNING) << "Particle collisions are enabled. This is still development feature. Ensure, that you know, what you doing";
#ifdef SWITCH_COULOMB_COLLISIONS_TA77S
  LOG_S(INFO) << "Using Takizuka and Abe Coulomb collisions scheme";
  LOG_S(INFO) << "\twith symmetric weighted particles correction [10.1002/ctpp.201700121]";
#elif defined(SWITCH_COULOMB_COLLISIONS_SK98)
  LOG_S(INFO) << "Using Sentoku and Kemp Coulomb collisions scheme [10.1143/JPSJ.67.4084, 10.1016/j.jcp.2008.03.043]";
#elif defined(SWITCH_COULOMB_COLLISIONS_P12)
  LOG_S(INFO) << "Using Perez et al. collisions scheme [10.1063/1.4742167]";
#endif
#endif

#ifdef SWITCH_MAXWELL_SOLVER_YEE
  LOG_S(INFO) << "Using Yee (FDTD) maxwellian solver [10.1109/TAP.1966.1138693]";
#endif

  // align macroparticles amount to particle species
  // with respect to normalization
#if defined(SWITCH_PLASMA_SPATIAL_REGULAR) || defined(SWITCH_PLASMA_SPATIAL_CENTERED)
  LOG_S(INFO) << "Using amount-dependent plama macroparticles spatial distribution";
  LOG_S(INFO) << "Amount of plasma macroparticles could be changed to satisfy spatial distribution requirement";
#endif
  for (auto p_s = particle_species.begin(); p_s != particle_species.end(); ++p_s)
  {
    (*p_s).macro_amount = (unsigned int)(geometry->cell_amount[0] * geometry->cell_amount[1] * norm);
    LOG_S(INFO) << "Macro amount for ``" << (*p_s).name << "'' speice is " << (*p_s).macro_amount;
  }

  // align macroparticles amount to particle beams
  // with respect to normalization
  for (auto b = particle_beams.begin(); b != particle_beams.end(); ++b)
  {
    double beam_r_grid_size = (*b).bunch_radius / geometry->cell_size[0];
    double beam_z_grid_size = (*b).bunch_length * (*b).bunches_amount / geometry->cell_size[1];

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
  if ((fmod(geometry->domains_amount[0], 2) != 0 && geometry->domains_amount[0] != 1) || (fmod(geometry->domains_amount[1], 2) != 0 && geometry->domains_amount[1] != 1))
    {
      LOG_S(FATAL) << "Amount of domains should be multiple of 2, or 1";
    }

  // if (geometry->cell_amount[0] / geometry->domains_amount[0] < MIN_DOMAIN_GRID_AMOUNT)
  //   {
  //     LOG_S(FATAL) << "Too small domain size by r. Must be ``" << MIN_DOMAIN_GRID_AMOUNT << "'' cells or more", 1);
  //   }

  // if (geometry->cell_amount[1] / geometry->domains_amount[1] < MIN_DOMAIN_GRID_AMOUNT)
  //   {
  //     LOG_S(FATAL) << "Too small domain size by z. Must be ``" << MIN_DOMAIN_GRID_AMOUNT << "'' cells or more", 1);
  //   }

  // if (geometry->cell_amount[0] / geometry->domains_amount[0] * geometry->cell_amount[1] / geometry->domains_amount[1] < MIN_DOMAIN_GRID_AMOUNT * MIN_DOMAIN_GRID_AMOUNT)
  //   {
  //     LOG_S(FATAL) << "Too small domain size. Must be ``" << MIN_DOMAIN_GRID_AMOUNT * MIN_DOMAIN_GRID_AMOUNT << "'' cells or more", 1);
  //   }

// try to figure out, where is electrons
  for (auto i = particle_species.begin(); i != particle_species.end(); ++i)
  {
    full_macro_amount += (*i).macro_amount; // calculate total number of macroparticles

    if (i->name.compare("Electrons") == 0
        || i->name.compare("electrons") == 0
        || i->name.compare("Electron") == 0
        || i->name.compare("electron") == 0
        || i->mass == 1)
    {
      electron_density = ((*i).left_density + (*i).right_density) / 2;
      electron_temperature = (*i).temperature;
    }

    if (i->name.compare("Ions") == 0
        || i->name.compare("ions") == 0
        || i->name.compare("Ion") == 0
        || i->name.compare("ion") == 0)
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
  double max_time_c = min(geometry->cell_size[0] / constant::LIGHT_VEL,
                          geometry->cell_size[1] / constant::LIGHT_VEL);
  double max_time_wp = 2 / plasma_freq / time_multiplicator;
  double max_time = min(max_time_c, max_time_wp);

  if (time->step > max_time)
  {
    LOG_S(ERROR) << "Too large time step: ``"
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
    if (geometry->cell_size[0] > debye_length * debye_multiplicator
	|| geometry->cell_size[1] > debye_length * debye_multiplicator)
    {
      LOG_S(ERROR) << "Too large grid size: ``"
	       << geometry->cell_size[0] << " x " << geometry->cell_size[1]
	       << " m.''. Should be less, than ``"
	       << debye_length * debye_multiplicator
		   << " m.''";
    }

    // \f$ L >> R_{debye} \f$
    if (geometry->size[0] < debye_length * debye_multiplicator
    	|| geometry->size[1] < debye_length * debye_multiplicator)
    {
      LOG_S(ERROR) << "Too small system size: ``"
    	       << geometry->size[0] << " x " << geometry->size[1]
    	       << " m.''. Should be more, than ``"
    	       <<  debye_length * debye_multiplicator
    		   << " m.''";
    }

    // \f$ L << N_{particles} * R_{debye} \f$
    if (geometry->size[0] > debye_length * full_macro_amount * debye_multiplicator
	|| geometry->size[1] > debye_length * full_macro_amount * debye_multiplicator)
    {
      LOG_S(ERROR) << "Too large system size: ``"
	       << geometry->size[0] << " x " << geometry->size[1]
	       << " m.''. Should be less, than ``"
	       <<  debye_length * full_macro_amount * debye_multiplicator
	       << " m.''";
    }
  }

  unsigned int grid_multiplicator = 10;

  if (geometry->cell_amount[0] * geometry->cell_amount[0] * grid_multiplicator > full_macro_amount)
  {
    LOG_S(WARNING) << "Too small amount of macroparticles: ``"
             << full_macro_amount
             << "''. Should be more, than ``"
             <<  geometry->cell_amount[0] * geometry->cell_amount[0] * grid_multiplicator
             << "''. You could get not relevant results";
  }

  return true;
}

string Cfg::cfg2str()
{
  // make some generic preparations
  value json_data_for_dump = json_data;
  json_data_for_dump.get<object>()["data"].get<object>()["data_root"].set<string>(".");

#ifdef ENABLE_DEBUG
  string str = json_data_for_dump.serialize(true);
#else
  string str = json_data_for_dump.serialize();
#endif // ENABLE_DEBUG
  return str;
}

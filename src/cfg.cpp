/*
 * This file is part of the PiCoPiC distribution
 * (https://github.com/cosmonaut-ok/PiCoPiC).
 * Copyright (c) 2020,2025 Alexander Vynnyk.
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

Cfg::Cfg (const std::string json_file_name)
{
  //! read given json file and parse it
  stringstream ss;
  ifstream f;
  // unsigned int i;

  time = new TimeSim();
  output_data = new save_data();

  // Read Json file
  f.open(json_file_name.c_str(), ios::binary);
  if ( !f.is_open() )
  {
    spdlog::critical("Can not read configuration file ``{}''", json_file_name);
    exit(EXIT_FAILURE);
  }
  ss << f.rdbuf();
  f.close();

  //! Parse Json data
  value v;
  ss >> v;
  string err = get_last_error();
  if ( !err.empty() )
  {
    spdlog::critical(err);
    exit(EXIT_FAILURE);
  }

  json_data = v;

  //! setup logging
  //! set logfile
  log_file = (char*)json_data.get<object>()["log_file"].get<string>().c_str();
  spdlog::init_thread_pool(8192, 1);
  auto stdout_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
  auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file, true);
  // TODO: set pattern. e.g.: console_sink->set_pattern("[%T.%e] %^%-8l%$: %v");
  std::vector<spdlog::sink_ptr> sinks{stdout_sink, file_sink};
  auto logger = std::make_shared<spdlog::async_logger>("default", sinks.begin(), sinks.end(), spdlog::thread_pool(),
                                                       spdlog::async_overflow_policy::block);
  spdlog::register_logger(logger);
  spdlog::set_default_logger(logger);

  // set log level
#ifdef ENABLE_DEBUG
  logger->set_level(spdlog::level::trace);
#else
  logger->set_level(spdlog::level::info);
#endif

  //! update json with additional items
  //! to inform about build options
  //! and package version
  // set package version to config metadata
  // json_data.get<object>()["package_version"] =
  // value((string)PACKAGE_VERSION);
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
#elif defined(SWITCH_PLASMA_SPATIAL_FLAT)
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

#ifdef SWITCH_MAXWELL_SOLVER_YEE
  o["maxwell_solver"] = value("yee");
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

#ifdef ENABLE_EXTERNAL_FIELDS
  init_external_fields();
#endif

  method_limitations_check();
}

Cfg::~Cfg () {}

void Cfg::init_particles ()
{
  //! initialize particles charge, mass, number,
  //! density and temperature for all particle species

  object& json_root = json_data.get<object>();
  picojson::array& particle_species_list = json_root["particles"].get<picojson::array>();

  for ( auto i = particle_species_list.begin(); i != particle_species_list.end(); ++i )
  {
    particle_specie p_s;

    object& o = i->get<object>();

    string specie_name = o["name"].get<string>();

    for ( auto sp = particle_species.begin(); sp != particle_species.end(); ++sp )
      if ( sp->name.compare(specie_name) == 0 )
      {
        spdlog::critical("Particles specie name ``{}'' already exists. Names "
                         "should be unique",
                         specie_name);
        exit(EXIT_FAILURE);
      }
    unsigned short current_specie_id = algo::common::hash_from_string(specie_name, SPECIE_HASH_SALT);

    spdlog::debug("Specie ID is ``{}'' for ``{}''", current_specie_id, specie_name);

    p_s.id = current_specie_id;
    p_s.name = specie_name;
    p_s.mass = (int)o["mass"].get<double>();
    p_s.charge = (int)o["charge"].get<double>();
    p_s.left_density = o["density"].get<double>();  // TODO: density profile still not
                                                    // supported, but planned
    p_s.right_density = o["density"].get<double>(); // TODO: density profile still not
                                                    // supported, but planned
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

  for ( auto i = probes_list.begin(); i != probes_list.end(); ++i )
  {
    probe p_p;

    object& o = i->get<object>();

    p_p.component = o["component"].get<string>().c_str();

    try
    {
      p_p.specie = o["specie"].get<string>().c_str();
    }
    catch ( std::exception& e )
    {
      p_p.specie = "unknown";
    }

    string p_shape = o["shape"].get<string>().c_str();

    //// generate probe path
    p_p.path = p_p.component + PATH_DELIMITER;

    // add path part for specie-based probes
    if ( p_p.component.compare("temperature") == 0 || p_p.component.compare("energy") == 0 ||
         p_p.component.compare("momentum/r") == 0 || p_p.component.compare("momentum/phi") == 0 ||
         p_p.component.compare("momentum/z") == 0 || p_p.component.compare("density") == 0 )
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
    {
      spdlog::critical("Unknown probe shape ``{}''", p_shape);
      exit(EXIT_FAILURE);
    }

    p_p.schedule = (int)o["schedule"].get<double>();

    if ( p_p.dims[0] > p_p.dims[2] || p_p.dims[1] > p_p.dims[3] )
    {
      spdlog::critical("Incorrect probe's {}/{}, shape: [{},{},{},{}]", p_p.component, p_shape, p_p.dims[0],
                       p_p.dims[2], p_p.dims[1], p_p.dims[3]);
      exit(EXIT_FAILURE);
    }
    else if ( p_p.dims[2] > geometry->cell_amount[0] )
    {
      spdlog::critical("Probe's {}/{} radius is out of simulation domain: {}. "
                       "Must be less than {}",
                       p_p.component, p_shape, p_p.dims[2], geometry->cell_amount[0]);
      exit(EXIT_FAILURE);
    }
    else if ( p_p.dims[3] > geometry->cell_amount[1] )
    {
      spdlog::critical("Probe's {}/{} longitude is out of simulation domain: "
                       "{}. Must be less, than {}",
                       p_p.component, p_shape, p_p.dims[3], geometry->cell_amount[1]);
      exit(EXIT_FAILURE);
    }
    else
      probes.push_back(p_p);
  }
}

void Cfg::init_beam ()
{
  //! initialize particles bunch data
  //! particles bunch should be injected to plasma

  object& json_root = json_data.get<object>();
  picojson::array& particle_beams_list = json_root["beams"].get<picojson::array>();

  if ( !particle_beams_list.empty() )
    for ( auto i = particle_beams_list.begin(); i != particle_beams_list.end(); ++i )
    {
      particle_beam p_b;

      object& o = i->get<object>();
      object& b_o = o["bunch"].get<object>();

      string beam_name = "beam_" + o["name"].get<string>();

      for ( auto bm = particle_beams.begin(); bm != particle_beams.end(); ++bm )
        if ( bm->name.compare(beam_name) == 0 )
        {
          spdlog::critical("Beam name ``{}'' is already exists. Beam names "
                           "should be unique. Exiting",
                           o["name"].get<string>());
          exit(EXIT_FAILURE);
        }

      // ! calculate current beam ID from its name:
      unsigned int current_beam_id = algo::common::hash_from_string(beam_name, BEAM_HASH_SALT);

      spdlog::debug("Beam ID is ``{}'' for ``{}''", current_beam_id, beam_name);

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

      p_b.radial_profile = b_o["profile"].get<object>()["radial"].get<string>();
      p_b.longitudinal_profile = b_o["profile"].get<object>()["longitudinal"].get<string>();

      p_b.current_bunch_number = 0;

      particle_beams.push_back(p_b);

      // ++current_beam_id;
    }
  else
    spdlog::info("There is no particle beams to inject");
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

  vector<size_t> _pml = {0, 0, 0, 0};
  vector<double> _sigma = {sigma_1, sigma_2};

  _pml[1] = floor(left_wall * n_grid_z);
  _pml[2] = floor(outer_wall * n_grid_r);
  _pml[3] = floor(right_wall * n_grid_z);

  geometry = new Geometry(_size, _dims, _pml, _sigma, _walls);
#else
  geometry = new Geometry(_size, _dims, _walls);
#endif
  vector<size_t> domains(2, 0);
  geometry->domains_amount = domains;

  //
#ifdef ENABLE_SINGLETHREAD
  spdlog::info("Singlethread mode. Reducing domains amount to 1");
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

void Cfg::init_output_data ()
{
  //! initialize file saving parameters, like path to computed data files,
  //! path to system state data files max frames number, placed to one file etc.
  object& json_root = json_data.get<object>()["data"].get<object>();

  output_data->data_root = (char*)json_root["data_root"].get<string>().c_str();

  output_data->compress = json_root["compression"].get<object>()["use"].get<bool>();
  output_data->compress_level = (int)json_root["compression"].get<object>()["level"].get<double>();
}

#ifdef ENABLE_EXTERNAL_FIELDS
void Cfg::init_external_fields ()
{
  object& json_root = json_data.get<object>()["external_fields"].get<object>();

  // init electric field
  object& json_el_field = json_root["electric"].get<object>();
  external_field el_f;
  el_f.name = "electric";
  el_f.formfactor = json_el_field["formfactor"].get<string>().c_str();
  // create array of params for electric field
  value::array el_params = json_el_field["params"].get<value::array>();
  for ( value item : el_params )
    el_f.params.push_back(item.get<double>());

  // init magnetic field
  object& json_magn_field = json_root["magnetic"].get<object>();
  external_field magn_f;
  magn_f.name = "magnetic";
  magn_f.formfactor = json_magn_field["formfactor"].get<string>().c_str();
  // create array of params for magnetic field
  value::array magn_params = json_magn_field["params"].get<value::array>();
  for ( value item : magn_params )
    magn_f.params.push_back(item.get<double>());

  external_fields.push_back(el_f);
  external_fields.push_back(magn_f);

  spdlog::info("Using External Fields. Electric field formfactor: ``{}'', "
               "Magnetic field formfactor: ``{}''",
               external_fields[0].formfactor, external_fields[1].formfactor);
}
#endif

void Cfg::weight_macro_amount ()
//! calculate alignment of macroparticles amount
//! to each particles specie and each beam
{
  double norm_sum = 0;

  // increase normalization sum for particle species (assume, that specie size
  // times species amount)
  norm_sum += geometry->cell_amount[0] * geometry->cell_amount[1] * particle_species.size();

  // increase normalization sum for particle beams
  for ( auto b = particle_beams.begin(); b != particle_beams.end(); ++b )
  {
    double beam_r_grid_size = (*b).bunch_radius / geometry->cell_size[0];
    double beam_z_grid_size = (*b).bunch_length * (*b).bunches_amount / geometry->cell_size[1];

    norm_sum += beam_r_grid_size * beam_z_grid_size * beam_macro_ratio;
  }

  // calculate normalization coeffitient
  double norm = macro_amount / norm_sum;

  // message about pusher, used in system
#ifdef SWITCH_PUSHER_BORIS
  spdlog::info("Using classic Boris particles pusher");
#elif defined(SWITCH_PUSHER_BORIS_ADAPTIVE)
  spdlog::info("Using adaptive Boris particles pusher");
#elif defined(SWITCH_PUSHER_BORIS_RELATIVISTIC)
  spdlog::info("Using fully relativistic Boris particles pusher");
#elif defined(SWITCH_PUSHER_VAY)
  spdlog::info("Using Vay particles pusher [10.1063/1.2837054]");
#elif defined(SWITCH_PUSHER_HC)
  spdlog::info("Using Higuera-Cary particles pusher [10.1063/1.4979989]");
#else
  // spdlog::critical("Undefined particles pusher used");
  //       exit(EXIT_FAILURE);
#endif

  // message about charge conservation scheme, used in system
#ifdef SWITCH_CCS_ZIGZAG
  spdlog::info("Using ZigZag charge conservation scheme "
               "[10.1016/S0010-4655(03)00437-5]");
#elif defined(SWITCH_CCS_VB)
  spdlog::info("Using Villasenor-Buneman charge conservation scheme "
               "[10.1016/0010-4655(92)90169-Y]");
#endif

  // message about coulomb colisions
#ifdef ENABLE_COULOMB_COLLISIONS
  spdlog::warn("Particle collisions are enabled. This is still development "
               "feature. Ensure, that you know, what you doing");
#ifdef SWITCH_COULOMB_COLLISIONS_TA77S
  spdlog::info("Using Takizuka and Abe Coulomb collisions scheme");
  spdlog::info("\twith symmetric weighted particles correction "
               "[10.1002/ctpp.201700121]");
#elif defined(SWITCH_COULOMB_COLLISIONS_SK98)
  spdlog::info("Using Sentoku and Kemp Coulomb collisions scheme "
               "[10.1143/JPSJ.67.4084, 10.1016/j.jcp.2008.03.043]");
#elif defined(SWITCH_COULOMB_COLLISIONS_P12)
  spdlog::info("Using Perez et al. collisions scheme [10.1063/1.4742167]");
#endif
#endif

#ifdef SWITCH_MAXWELL_SOLVER_YEE
  spdlog::info("Using Yee (FDTD) maxwellian solver [10.1109/TAP.1966.1138693]");
#endif

  // align macroparticles amount to particle species
  // with respect to normalization
#if defined(SWITCH_PLASMA_SPATIAL_REGULAR) || defined(SWITCH_PLASMA_SPATIAL_CENTERED)
  spdlog::info("Using amount-dependent plama macroparticles spatial distribution");
  spdlog::info("Amount of plasma macroparticles could be changed to satisfy "
               "spatial distribution requirement");
#endif
  for ( auto p_s = particle_species.begin(); p_s != particle_species.end(); ++p_s )
  {
    (*p_s).macro_amount = (unsigned int)(geometry->cell_amount[0] * geometry->cell_amount[1] * norm);
    spdlog::info("Macro amount for ``{}'' specie is ``{}''", (*p_s).name, (*p_s).macro_amount);
  }

  // align macroparticles amount to particle beams
  // with respect to normalization
  for ( auto b = particle_beams.begin(); b != particle_beams.end(); ++b )
  {
    double beam_r_grid_size = (*b).bunch_radius / geometry->cell_size[0];
    double beam_z_grid_size = (*b).bunch_length * (*b).bunches_amount / geometry->cell_size[1];

    (*b).macro_amount = (unsigned int)(beam_r_grid_size * beam_z_grid_size * norm * beam_macro_ratio);
    spdlog::info("Macro amount for ``{}'' beam is ``{}''", (*b).name, (*b).macro_amount);
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
  if ( (fmod(geometry->domains_amount[0], 2) != 0 && geometry->domains_amount[0] != 1) ||
       (fmod(geometry->domains_amount[1], 2) != 0 && geometry->domains_amount[1] != 1) )
  {
    spdlog::critical("Amount of domains should be multiple of 2, or 1");
    exit(EXIT_FAILURE);
  }

  // if (geometry->cell_amount[0] / geometry->domains_amount[0] <
  // MIN_DOMAIN_GRID_AMOUNT)
  //   {
  //     spdlog::critical("Too small domain size by r. Must be ``{}'' cells or
  //     more", MIN_DOMAIN_GRID_AMOUNT); exit(EXIT_FAILURE);
  //   }

  // if (geometry->cell_amount[1] / geometry->domains_amount[1] <
  // MIN_DOMAIN_GRID_AMOUNT)
  //   {
  //     spdlog::critical("Too small domain size by z. Must be ``{}'' cells or
  //     more", MIN_DOMAIN_GRID_AMOUNT); exit(EXIT_FAILURE);
  //   }

  // if (geometry->cell_amount[0] / geometry->domains_amount[0] *
  // geometry->cell_amount[1] / geometry->domains_amount[1] <
  // MIN_DOMAIN_GRID_AMOUNT * MIN_DOMAIN_GRID_AMOUNT)
  //   {
  //     spdlog::critical("Too small domain size. Must be ``{}'' cells or more",
  //     MIN_DOMAIN_GRID_AMOUNT * MIN_DOMAIN_GRID_AMOUNT); exit(EXIT_FAILURE);
  //   }

  // try to figure out, where is electrons
  for ( auto i = particle_species.begin(); i != particle_species.end(); ++i )
  {
    full_macro_amount += (*i).macro_amount; // calculate total number of macroparticles

    if ( i->name.compare("Electrons") == 0 || i->name.compare("electrons") == 0 || i->name.compare("Electron") == 0 ||
         i->name.compare("electron") == 0 || i->mass == 1 )
    {
      electron_density = ((*i).left_density + (*i).right_density) / 2;
      electron_temperature = (*i).temperature;
    }

    if ( i->name.compare("Ions") == 0 || i->name.compare("ions") == 0 || i->name.compare("Ion") == 0 ||
         i->name.compare("ion") == 0 )
    {
      ion_density = ((*i).left_density + (*i).right_density) / 2;
      ion_temperature = (*i).temperature;
    }
  }

  if ( electron_density == 0 )
  {
    spdlog::critical("There is no electrons present in system");
    exit(EXIT_FAILURE);
  }

  double plasma_freq = phys::plasma::plasma_frequency(electron_density);

  // TODO: WTF
  unsigned int time_multiplicator = 100;

  // particle should not fly more, than 1 cell per time step
  double max_time_c = min(geometry->cell_size[0] / constant::LIGHT_VEL, geometry->cell_size[1] / constant::LIGHT_VEL);
  double max_time_wp = 2 / plasma_freq / time_multiplicator;
  double max_time = min(max_time_c, max_time_wp);

  if ( time->step > max_time )
    spdlog::error("Too large time step: ``{} s.''. Should be less, than ``{} s.''", time->step, max_time);

  // cell size d{r,z} < \lambda debye

  double debye_length =
    phys::plasma::debye_length(electron_density, ion_density, ion_temperature, electron_temperature);
  unsigned int debye_multiplicator = 100;

  if ( electron_temperature != 0 )
  {
    if ( geometry->cell_size[0] > debye_length * debye_multiplicator ||
         geometry->cell_size[1] > debye_length * debye_multiplicator )
      spdlog::error("Too large grid size: ``{}x{} m.''. Should be less, than ``{} m.''", geometry->cell_size[0],
                    geometry->cell_size[1], debye_length * debye_multiplicator);

    // \f$ L >> R_{debye} \f$
    if ( geometry->size[0] < debye_length * debye_multiplicator ||
         geometry->size[1] < debye_length * debye_multiplicator )
      spdlog::error("Too small system size: ``{}x{} m.''. Should be more, than ``{} m.''", geometry->size[0],
                    geometry->size[1], debye_length * debye_multiplicator);

    // \f$ L << N_{particles} * R_{debye} \f$
    if ( geometry->size[0] > debye_length * full_macro_amount * debye_multiplicator ||
         geometry->size[1] > debye_length * full_macro_amount * debye_multiplicator )
      spdlog::error("Too large system size: ``{}x{} m.''. Should be less, than ``{} m.''", geometry->size[0],
                    geometry->size[1], debye_length * full_macro_amount * debye_multiplicator);
  }

  unsigned int grid_multiplicator = 10;

  if ( geometry->cell_amount[0] * geometry->cell_amount[0] * grid_multiplicator > full_macro_amount )
  {
    spdlog::warn("Too small amount of macroparticles: ``{}''. Should be more, "
                 "than ``{}''. You could get not relevant results",
                 full_macro_amount, geometry->cell_amount[0] * geometry->cell_amount[0] * grid_multiplicator);
  }

  return true;
}

string Cfg::cfg2str ()
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

picojson::value Cfg::cfg2value () { return json_data; }

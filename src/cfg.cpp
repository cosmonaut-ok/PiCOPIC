#include "cfg.hpp"
#include <typeinfo>

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
    LOG_CRIT(err, 1);

  json_data = v;

  //! set debug option
  debug = json_data.get<object>()["debug"].get<bool>();

  //! check if use hdf5 backend
  use_hdf5 = json_data.get<object>()["hdf5"].get<bool>();

  //! check if use hdf5 backend
  macro_amount = json_data.get<object>()["macro_amount"].get<double>();

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
    // p_s.macro_amount = o["macro_amount"].get<double>();
    p_s.left_density = o["density"].get<object>()["left"].get<double>();
    p_s.right_density = o["density"].get<object>()["right"].get<double>();
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
      p_p.z_start = (int)o["z"].get<double>();
      p_p.z_end = -1;
    }
    else if ( p_shape.compare("row") == 0 )
    {
      p_p.shape = 2;
      p_p.r_start = (int)o["r"].get<double>();
      p_p.r_end = -1;
      p_p.z_start = -1;
      p_p.z_end = -1;
    }
    else if ( p_shape.compare("dot") == 0 )
    {
      p_p.shape = 3;
      p_p.r_start = (int)o["r"].get<double>();
      p_p.r_end = -1;
      p_p.z_start = (int)o["z"].get<double>();
      p_p.z_end = -1;
    }
    else
      LOG_CRIT("Unknown probe shape ``" << p_shape << "''", 1);

    p_p.schedule = (int)o["schedule"].get<double>();

    if (p_p.r_start > p_p.r_start || p_p.z_start > p_p.z_start)
    {
      LOG_CRIT("Incorrect probe's " << p_p.component << "/" << p_shape << " shape: ["
               << p_p.r_start << "," << p_p.r_end << ","
               << p_p.z_start << "," << p_p.z_end << "]", 1);
    }
    else if (p_p.r_end > geometry->r_grid_amount)
    {
      LOG_CRIT("Probe's " << p_p.component << "/" << p_shape << " radius is out of simulation area: "
               << p_p.r_end << ". Must be less, than "
               << geometry->r_grid_amount, 1);
    }
    else if (p_p.z_end > geometry->z_grid_amount)
    {
      LOG_CRIT("Probe's " << p_p.component << "/" << p_shape << " longitude is out of simulation area: "
               << p_p.z_end << ". Must be less, than "
               << geometry->z_grid_amount, 1);
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
  object& json_root_pml = json_data.get<object>()["pml"].get<object>();

  double radius = json_root["radius"].get<double>();
  double longitude = json_root["longitude"].get<double>();

  int n_grid_r = (int)json_root["grid_size"].get<object>()["radius"].get<double>();
  int n_grid_z = (int)json_root["grid_size"].get<object>()["longitude"].get<double>();

  // init PML
  double left_wall = json_root_pml["left_wall"].get<double>();
  double right_wall = json_root_pml["right_wall"].get<double>();
  double outer_wall = json_root_pml["outer_wall"].get<double>();
  double sigma_1 = json_root_pml["sigma_1"].get<double>();
  double sigma_2 = json_root_pml["sigma_2"].get<double>();

  geometry = new Geometry(radius, longitude, 0, n_grid_r, 0, n_grid_z,
                          left_wall, right_wall, outer_wall, sigma_1, sigma_2,
                          true, true, true, true); // init all walls

  //
  geometry->areas_by_r = (int)json_root["areas_amount"].get<object>()["radius"].get<double>();
  geometry->areas_by_z = (int)json_root["areas_amount"].get<object>()["longitude"].get<double>();
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

  // make beam macroparticles X times smaller, than
  // macroparticles of regular particle species
  const unsigned int beam_macro_const = 2;

// increase normalization sum for particle species (assume, that specie size times species amount)
  norm_sum += geometry->r_grid_amount * geometry->z_grid_amount * particle_species.size();

// increase normalization sum for particle beams
  for (auto b = particle_beams.begin(); b != particle_beams.end(); ++b)
  {
    double beam_r_grid_size = (*b).bunch_radius / geometry->r_cell_size;
    double beam_z_grid_size = (*b).bunch_length * (*b).bunches_amount / geometry->z_cell_size;

    norm_sum += beam_r_grid_size * beam_z_grid_size * beam_macro_const;
  }

  // calculate normalization coeffitient
  double norm = macro_amount / norm_sum;

  // align macroparticles amount to particle species
  // with respect to normalization
  for (auto p_s = particle_species.begin(); p_s != particle_species.end(); ++p_s)
  {
    (*p_s).macro_amount = (unsigned int)(geometry->r_grid_amount * geometry->z_grid_amount * norm);
    LOG_INFO("Macro amount for ``" << (*p_s).name << "'' speice is " << (*p_s).macro_amount);
  }

  // align macroparticles amount to particle beams
  // with respect to normalization
  for (auto b = particle_beams.begin(); b != particle_beams.end(); ++b)
  {
    double beam_r_grid_size = (*b).bunch_radius / geometry->r_cell_size;
    double beam_z_grid_size = (*b).bunch_length * (*b).bunches_amount / geometry->z_cell_size;

    (*b).macro_amount = (unsigned int)(beam_r_grid_size * beam_z_grid_size * norm * beam_macro_const);
    LOG_INFO("Macro amount for ``" << (*b).name << "'' beam is " << (*b).macro_amount);
  }
}

bool Cfg::method_limitations_check ()
{
  double full_macro_amount = 0;
  double electron_density = 0;
  double electron_temperature = 0;

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
  }

  if (electron_density == 0 || electron_temperature == 0)
  {
    LOG_CRIT("There is no electrons present in system", 1);
  }

  double plasma_freq = sqrt(electron_density * EL_CHARGE * EL_CHARGE / (EL_MASS * EPSILON0));

  // TODO: WTF
  unsigned int time_multiplicator = 100;

  // particle should not fly more, than 1 cell per time step
  double max_time_c = min(geometry->r_cell_size / constant::LIGHT_VEL,
                          geometry->z_cell_size / constant::LIGHT_VEL);
  double max_time_wp = 2 / plasma_freq / time_multiplicator;
  double max_time = min(max_time_c, max_time_wp);

  if (time->step > max_time)
  {
    LOG_CRIT("Too large time step: ``"
             << time->step
             << " s.''. Should be less, than ``"
             << max_time
             << " s.''", 1);
  }

  // cell size d{r,z} < \lambda debye

  // 7400 is a coefficient, when T is in electron volts and N is in \f$ m^-3 \f$
  double debye_length = 7400 * sqrt(electron_temperature / electron_density);
  unsigned int debye_multiplicator = 100;

  if (geometry->r_cell_size > debye_length * debye_multiplicator
      || geometry->z_cell_size > debye_length * debye_multiplicator)
  {
    LOG_CRIT("Too large grid size: ``"
             << geometry->r_cell_size << " x " << geometry->z_cell_size
             << " m.''. Should be less, than ``"
             << debye_length * debye_multiplicator
             << " m.''", 1);
  }

  // \f$ L >> R_{debye} \f$
  if (geometry->r_size < debye_length * debye_multiplicator
      || geometry->z_size < debye_length * debye_multiplicator)
  {
    LOG_CRIT("Too small system size: ``"
             << geometry->r_size << " x " << geometry->z_size
             << " m.''. Should be more, than ``"
             <<  debye_length * debye_multiplicator
             << " m.''", 1);
  }

  // \f$ L << N_{particles} * R_{debye} \f$
  if (geometry->r_size > debye_length * full_macro_amount * debye_multiplicator
      || geometry->z_size > debye_length * full_macro_amount * debye_multiplicator)
  {
    LOG_CRIT("Too large system size: ``"
             << geometry->r_size << " x " << geometry->z_size
             << " m.''. Should be less, than ``"
             <<  debye_length * full_macro_amount * debye_multiplicator
             << " m.''", 1);

  }

  unsigned int grid_multiplicator = 10;

  if (geometry->r_grid_amount * geometry->r_grid_amount * grid_multiplicator > full_macro_amount)
  {
    LOG_WARN("Too small summary number of macroparticles: ``"
             << full_macro_amount
             << "''. Should be more, than ``"
             <<  geometry->r_grid_amount * geometry->r_grid_amount * grid_multiplicator
             << "''. You could get not relevant results");
  }

  return true;
}

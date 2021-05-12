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

// https://portal.hdfgroup.org/pages/viewpage.action?pageId=50073884

#include "outEngine/outEngineHDF5.hpp"

using namespace std;
using namespace HighFive;

OutEngineHDF5::OutEngineHDF5 ( string _path,
                               vector<size_t> _size,
                               vector<size_t> _offset,
                               bool _append, unsigned short _compress )
  :OutEngine ( _path, _size, _offset, _append, _compress)
{
  true;
}

OutEngineHDF5::OutEngineHDF5 ( HighFive::File* _file, string _path,
                               vector<size_t> _size,
                               vector<size_t> _offset,
                               bool _append, unsigned short _compress )
  :OutEngineHDF5 ( _path, _size, _offset, _append, _compress)
{
  data_file = _file;

  if (compress)
    LOG_S(WARNING) << "compression is not supported by plaintext output engine";
}

// dummy methods just to show interface
void OutEngineHDF5::write_cub ( size_t _num, vector <vector <vector<double>>> data )
{
  LOG_S(ERROR) << "Not implemented";
}

void OutEngineHDF5::write_rec ( size_t _num, vector< vector<double>> data )
{
  try
  {
    // Create the dataset
    DataSet dataset = data_file->getDataSet ( path );

    vector<size_t> local_offset = {_num};
    for (auto i = offset.begin(); i != offset.end(); ++i)
      local_offset.push_back((*i));
    vector<size_t> local_size = {1, data.size(), data[0].size()};

    vector <vector <vector <double> > > local_data
      (1, std::vector <std::vector <double> >
       (data.size(), std::vector <double>
        (data[0].size(), 0.) ) );
    local_data[0] = data;

    dataset.select(local_offset, local_size).write( local_data );
  }
  catch (Exception& error)
  {
    LOG_S(ERROR) << error.what();
    LOG_S(FATAL) << "Can not write data to dataset ``/" << path << "''";
  }
}

void OutEngineHDF5::write_vec ( size_t _num, vector<double> data)
{
  try
  {
    DataSet dataset = data_file->getDataSet ( path );

    vector<size_t> local_offset = {_num};
    for (auto i = offset.begin(); i != offset.end(); ++i)
      local_offset.push_back((*i));
    vector<size_t> local_size = {1, data.size()};

    vector <vector <double> > local_data
      (1, vector <double>
       (data.size(), 0.) );
    local_data[0] = data;

    dataset.select(local_offset, local_size).write( local_data );
  }
  catch (Exception& error)
  {
    LOG_S(ERROR) << error.what();
    LOG_S(FATAL) << "Can not write data to dataset ``/" << path << "''";
  }
}

void OutEngineHDF5::write_dot(size_t _num, double data)
{
  try
  {
    DataSet dataset = data_file->getDataSet ( path );

    dataset.select({_num, 0}, {1, 1}).write( data );
  }
  catch (Exception& error)
  {
    LOG_S(ERROR) << error.what();
    LOG_S(FATAL) << "Can not write data to dataset ``/" << path << "''";
  }
}

void OutEngineHDF5::create_dataset()
{
  string group_name = path.substr(0, path.find_last_of("\\/"));

  try
  {
    LOG_S(MAX) << "Creating group ``" << group_name << "''";
    Group group_instance = data_file->createGroup('/' + group_name);
    LOG_S(MAX) << "Group ``" << group_name << "'' created";

    std::vector<size_t> dims_init = {1};
    std::vector<size_t> dims_final = {DataSpace::UNLIMITED};
    for (unsigned int i=0; i < size.size(); ++i)
    {
      dims_init.push_back(size[i]);
      dims_final.push_back(size[i]);
    }

    // Create a dataspace with initial shape and max shape
    DataSpace dataspace = DataSpace(dims_init, dims_final);

    // Use chunking
    DataSetCreateProps props;
    vector<hsize_t> chunk_dims;
    for (auto i = dims_init.begin(); i != dims_init.end(); ++i)
      chunk_dims.push_back((*i));

    props.add(Chunking(chunk_dims)); // std::vector<hsize_t>{1, 4, 5}));

    // Create the dataset
    DataSet dataset = data_file->createDataSet ( path, dataspace,
                                                 AtomicType<double>(), props );

    // // Create the dataset
    // DataSet dataset = data_file.createDataSet<double>(_path, dataspace);
    // exit(0);
  }
  catch (FileException& error)
  {
    LOG_S(ERROR) << "FileException: " << error.what();
  }

  catch (GroupException& error)
  {
    LOG_S(ERROR) << error.what();
    LOG_S(FATAL) << "Can not write data to dataset ``/" << path << "''";
  }

  // catch failure caused by the DataSet operations
  catch (DataSetException& error)
  {
    LOG_S(ERROR) << "DataSetException: " << error.what();
  }


  // try
  // {
  //   Group group_instance = data_file->getGroup('/' + group);

  //   std::vector<size_t> dims;
  //   for (unsigned int i=0; i < _dims.size(); ++i)
  //     dims.push_back((size_t)_dims[i]);

  //   // Create the dataset
  //   DataSet dataset = group_instance.createDataSet<double>(_name, DataSpace(dims));
  // }
  // catch (FileException& error)
  // {
  //   LOG_S(ERROR) << "FileException: " << error.what();
  // }

  // catch (GroupException& error)
  // {
  //   LOG_S(ERROR) << error.what();
  //   LOG_S(FATAL) << "Can not write dataset ``" << _name << "'' to group ``/" << group << "''";
  // }

  // // catch failure caused by the DataSet operations
  // catch (DataSetException& error)
  // {
  //   LOG_S(ERROR) << "DataSetException: " << error.what();
  // }
}

void OutEngineHDF5::extend_dataset(size_t _num)
{
  try
  {
    // Get the dataset
    DataSet dataset = data_file->getDataSet('/' + path);

    if (dataset.getDimensions()[0] > _num)
      return;
    else
    {
      vector<size_t> local_size = {_num+1};
      for (auto i = size.begin(); i != size.end(); ++i)
        local_size.push_back((*i));

      dataset.resize(local_size);
    }
  }
  catch (FileException& error)
  {
    LOG_S(ERROR) << "FileException: " << error.what();
  }

  catch (GroupException& error)
  {
    LOG_S(ERROR) << error.what();
    LOG_S(FATAL) << "Can not write data to dataset ``/" << path << "''";
  }

  // catch failure caused by the DataSet operations
  catch (DataSetException& error)
  {
    LOG_S(ERROR) << "DataSetException: " << error.what();
  }


  // try
  // {
  //   Group group_instance = data_file->getGroup('/' + group);

  //   std::vector<size_t> dims;
  //   for (unsigned int i=0; i < _dims.size(); ++i)
  //     dims.push_back((size_t)_dims[i]);

  //   // Create the dataset
  //   DataSet dataset = group_instance.createDataSet<double>(_name, DataSpace(dims));
  // }
  // catch (FileException& error)
  // {
  //   LOG_S(ERROR) << "FileException: " << error.what();
  // }

  // catch (GroupException& error)
  // {
  //   LOG_S(ERROR) << error.what();
  //   LOG_S(FATAL) << "Can not write dataset ``" << _name << "'' to group ``/" << group << "''";
  // }

  // // catch failure caused by the DataSet operations
  // catch (DataSetException& error)
  // {
  //   LOG_S(ERROR) << "DataSetException: " << error.what();
  // }
}

void OutEngineHDF5::write_metadata(picojson::value _metadata)
{
#ifdef ENABLE_DEBUG
  string metadata_str = _metadata.serialize(true);
#else
  string metadata_str = _metadata.serialize();
#endif // ENABLE_DEBUG

  try
  {
    Group group = data_file->getGroup("/");

    LOG_S(MAX) << "Writing metadata";

    Attribute attribute = group.createAttribute<std::string>("metadata", DataSpace::From(metadata_str));
    attribute.write(metadata_str);

    // software name, version, build flags and options
    std::string package_name = (string)PACKAGE_NAME;
    std::string package_version = (string)PACKAGE_VERSION;
    std::string package_build_flags = (string)CXXFLAGS;
    Attribute pkg_name_a = group
      .createAttribute<std::string>("software", DataSpace::From(package_name));
    Attribute pkg_ver_a = group
      .createAttribute<std::string>("softwareVersion",
                                    DataSpace::From(package_version));
    Attribute pkg_build_flags_a = group
      .createAttribute<std::string>("softwareBuildFlags",
                                    DataSpace::From(package_build_flags));

    pkg_name_a.write(package_name);
    pkg_ver_a.write(package_version);
    pkg_build_flags_a.write(package_build_flags);

    // softwareDependencies
    std::string package_deps = (string)PACKAGE_DEPS;
    Attribute pkg_deps_a = group
      .createAttribute<std::string>("softwareDependencies",
                                    DataSpace::From(package_deps));
    pkg_deps_a.write(package_deps);

    // debug
    bool is_debug_enabled = _metadata.get<picojson::object>()["build_options"]
      .get<picojson::object>()["debug"]
      .get<bool>();
    Attribute is_debug_enabled_a = group
      .createAttribute<bool>("softwareConfigDebug",
                                    DataSpace::From(is_debug_enabled));
    is_debug_enabled_a.write(is_debug_enabled);

    // accept IEEE
    bool is_ieee_enabled = _metadata.get<picojson::object>()["build_options"]
      .get<picojson::object>()["ieee"]
      .get<bool>();
    Attribute is_ieee_enabled_a = group
      .createAttribute<bool>("softwareConfigIEEE",
                             DataSpace::From(is_ieee_enabled));
    is_ieee_enabled_a.write(is_ieee_enabled);

    // current deposition
    std::string current_depos = _metadata.get<picojson::object>()["build_options"]
      .get<picojson::object>()["current_deposition"]
      .get<std::string>();
    Attribute current_depos_a = group
      .createAttribute<std::string>("softwareConfigCurrentDeposition",
                                    DataSpace::From(current_depos));
    current_depos_a.write(current_depos);

    // field solver
    std::string maxwell_solver = _metadata.get<picojson::object>()["build_options"]
      .get<picojson::object>()["maxwell_solver"]
      .get<std::string>();
    Attribute maxwell_solver_a = group
      .createAttribute<std::string>("softwareConfigFieldSolver",
                                    DataSpace::From(maxwell_solver));
    maxwell_solver_a.write(maxwell_solver);

    // pusher
    std::string prtls_pusher = _metadata.get<picojson::object>()["build_options"]
      .get<picojson::object>()["particles_pusher"]
      .get<std::string>();
    Attribute prtls_pusher_a = group
      .createAttribute<std::string>("softwareConfigParticlesPusher",
                                    DataSpace::From(prtls_pusher));
    prtls_pusher_a.write(prtls_pusher);

    // plasma_velocity_distribution
    std::string plasma_velocity_distribution = _metadata.get<picojson::object>()["build_options"]
      .get<picojson::object>()["plasma_velocity_distribution"]
      .get<std::string>();
    Attribute plasma_velocity_distribution_a = group
      .createAttribute<std::string>("softwareConfigPlasmaVelocityDistribution",
                                    DataSpace::From(plasma_velocity_distribution));
    plasma_velocity_distribution_a.write(plasma_velocity_distribution);

    // plasma_spatial_distribution
    std::string plasma_spatial_distribution = _metadata.get<picojson::object>()["build_options"]
      .get<picojson::object>()["plasma_spatial_distribution"]
      .get<std::string>();
    Attribute plasma_spatial_distribution_a = group
      .createAttribute<std::string>("softwareConfigPlasmaSpatialDistribution",
                                    DataSpace::From(plasma_spatial_distribution));
    plasma_spatial_distribution_a.write(plasma_spatial_distribution);

    // comment
    std::string comment = _metadata.get<picojson::object>()["comment"].get<std::string>();
    Attribute comment_a = group
      .createAttribute<std::string>("comment",
                                    DataSpace::From(comment));
    comment_a.write(comment);

  }
  catch(GroupException& error)
  {
    LOG_S(MAX) << "Group ``metadata'' already exists. Skipping";
  }

  catch(DataSetException& error)
  {
    LOG_S(MAX) << "metadata already exists. Skipping";
  }
}

void OutEngineHDF5::create_path()
{
  string group_name = path.substr(0, path.find_last_of("\\/"));
  try
  {
    LOG_S(MAX) << "Creating group ``" << group_name << "''";
    Group group_instance = data_file->createGroup('/' + group_name);
    LOG_S(MAX) << "Group ``" << group_name << "'' created";
  }
  catch (const GroupException&)
  {
    LOG_S(MAX) << "Group ``" << group_name << "'' already exists. Skipping";
  }
}

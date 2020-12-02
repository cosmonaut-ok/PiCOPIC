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

OutEngineHDF5::OutEngineHDF5 ( std::string _group,
                                   std::vector<size_t> _offset,
                                   bool _append, unsigned short _compress )
  :OutEngine ( _group, _offset, _append, _compress)
{
  true;
}

OutEngineHDF5::OutEngineHDF5 ( HighFive::File* _file, std::string _group,
				   std::vector<size_t> _offset,
				   bool _append, unsigned short _compress )
  :OutEngineHDF5 ( _group, _offset, _append, _compress)
{
  data_file = _file;
}

// dummy methods just to show interface
void OutEngineHDF5::write_rec ( string _name, vector< vector<double>> data )
{
  try
  {
    Group group_instance = data_file->getGroup('/' + group);

    // Create the dataset
    DataSet dataset = group_instance.getDataSet( _name );

    vector<size_t> size = {data.size(), data[0].size()};
    dataset.select(offset, size).write( data );
  }
  catch (Exception& error)
  {
    LOG_S(ERROR) << error.what();
    LOG_S(FATAL) << "Can not write dataset ``" << _name << "'' to group ``/" << group << "''";
  }
}

void OutEngineHDF5::write_vec(string _name, vector<double> data)
{
  try
  {
    Group group_instance = data_file->getGroup('/' + group);

    // Create the dataset
    DataSet dataset = group_instance.getDataSet( _name );

    dataset.select(offset, {data.size()}).write( data );
  }
  catch (Exception& error)
  {
    LOG_S(ERROR) << error.what();
    LOG_S(FATAL) << "Can not write dataset ``" << _name << "'' to group ``/" << group << "''";
  }
}

void OutEngineHDF5::write_dot(string _name, double data)
{
  try
  {
    Group group_instance = data_file->getGroup('/' + group);

    // Create the dataset
    DataSet dataset = group_instance.getDataSet( _name );
    dataset.write( data );
  }
  catch (Exception& error)
  {
    LOG_S(ERROR) << error.what();
    LOG_S(FATAL) << "Can not write dataset ``" << _name << "'' to group ``/" << group << "''";
  }
}

void OutEngineHDF5::create_dataset(std::string _name, vector<unsigned int> _dims)
{
  try
  {
    Group group_instance = data_file->getGroup('/' + group);

    std::vector<size_t> dims;
    for (unsigned int i=0; i < _dims.size(); ++i)
      dims.push_back((size_t)_dims[i]);

    // Create the dataset
    DataSet dataset = group_instance.createDataSet<double>(_name, DataSpace(dims));
  }
  catch (FileException& error)
  {
    LOG_S(ERROR) << "FileException: " << error.what();
  }

  catch (GroupException& error)
  {
    LOG_S(ERROR) << error.what();
    LOG_S(FATAL) << "Can not write dataset ``" << _name << "'' to group ``/" << group << "''";
  }

  // catch failure caused by the DataSet operations
  catch (DataSetException& error)
  {
    LOG_S(ERROR) << "DataSetException: " << error.what();
  }
}

void OutEngineHDF5::write_metadata(std::string _metadata)
{
  try
    {
      LOG_S(MAX) << "Creating group ``/metadata''";
      Group group = data_file->createGroup("/metadata");
      LOG_S(MAX) << "Group ``/metadata'' created";

      LOG_S(MAX) << "Writing metadata";

      Attribute attribute = group.createAttribute<std::string>("metadata", DataSpace::From(_metadata));

      attribute.write(_metadata);
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
  try
  {
    LOG_S(MAX) << "Creating group ``" << group << "''";
    Group g = data_file->createGroup(group);
    LOG_S(MAX) << "Group ``" << group << "'' created";
  }
  catch (const GroupException&)
  {
    LOG_S(MAX) << "Group ``" << group << "'' already exists. Skipping";
  }

  if (compress)
    LOG_S(WARNING) << "compression is not supported by plaintext output engine";
}

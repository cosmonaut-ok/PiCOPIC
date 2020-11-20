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

OutEngineHDF5::OutEngineHDF5 (File* _file, string _path, string _subpath,
                              unsigned int _shape, int *_size,
                              bool _append, bool _compress,
                              unsigned int _compress_level)
  : OutEngine (_path, _subpath, _shape, _size, _append, _compress, _compress_level)
{
  out_file = _file;

  // create group with subpath name recursively
  std::vector<std::string> group_array;
  const char delim = '/';
  algo::common::splitstr(subpath, delim, group_array);

  string group_name;

  for (auto i = group_array.begin(); i != group_array.end(); ++i)
  {
    group_name += delim;
    group_name += (*i);

    try
    {
      LOG_S(MAX) << "Creating group ``" << group_name << "''";
      Group g = out_file->createGroup(group_name.c_str());
      LOG_S(MAX) << "Group ``" << group_name << "'' created";
    }
    catch (const GroupException&)
    {
      LOG_S(MAX) << "Group ``" << group_name << "'' already exists. Skipping";
    }
  }

  if (compress)
    LOG_S(WARNING) << "compression is not supported by plaintext output engine";
}

void OutEngineHDF5::write_rec(string _name, Grid<double> data)
{
  // FIXME: workaround: set stack size to required value
  unsigned int size_x = size[2] - size[0];
  unsigned int size_y = size[3] - size[1];

  const rlim_t kStackSize = size_x * size_y * sizeof(double) * 1024;
  struct rlimit rl;
  int result;

  result = getrlimit(RLIMIT_STACK, &rl);
  if (result == 0)
    if (rl.rlim_cur < kStackSize)
    {
      rl.rlim_cur = kStackSize;
      result = setrlimit(RLIMIT_STACK, &rl);
    };

  try
  {
    Group group = out_file->getGroup('/' + subpath);

    std::vector<size_t> dims{size_x, size_y};

    // Create the dataset
    DataSet dataset = group.createDataSet<double>(_name, DataSpace(dims));

    vector<vector<double>> tmparr(size_x , vector<double> (size_y, 0));
    for (unsigned int i = 0; i < size_x; ++i)
      for (unsigned int j = 0; j < size_y; ++j)
        tmparr[i][j] = data(i+size[0], j+size[1]);

    dataset.write(tmparr);
  }
  catch (FileException& error)
  {
    LOG_S(ERROR) << "FileException: " << error.what();
  }

  catch (GroupException& error)
  {
    LOG_S(ERROR) << error.what();
    LOG_S(FATAL) << "Can not write dataset ``" << _name << "'' to group ``/" << subpath << "''";
  }

  // catch failure caused by the DataSet operations
  catch (DataSetException& error)
  {
    LOG_S(ERROR) << "DataSetException: " << error.what();
  }
}

void OutEngineHDF5::write_vec(string _name, Grid<double> data)
{
  try
  {
    Group group = out_file->getGroup('/' + subpath);

    std::vector<size_t> dims{2};

    // vector by Z-component with fixed r (r_begin)
    if (size[1] == -1 && size[0] == -1 && size[3] == -1)
      dims[0] = data.y_size;
    else if (size[0] == -1 && size[2] == -1 && size[1] == -1)
      dims[0] = data.x_size;
    else
      LOG_S(FATAL) << "Incorrect shape for vector output";

    // Create the dataset
    DataSet dataset = group.createDataSet<double>(_name, DataSpace(dims));

    vector<double> tmparr(dims[0]);
    for (hsize_t i = 0; i < dims[0]; ++i)
      if (dims[0] == data.y_size)
        tmparr[i] = data(size[2]-1, i);
      else if (dims[0] == data.x_size)
        tmparr[i] = data(i, size[3]-1);

    dataset.write(tmparr);
  }
  catch(FileException& error)
  {
    LOG_S(ERROR) << error.what();
  }

  catch (GroupException& error)
  {
    LOG_S(ERROR) << error.what();
    LOG_S(FATAL) << "Can not write dataset ``" << _name << "'' to group ``/" << subpath << "''";
  }
  // catch failure caused by the DataSet operations
  catch(DataSetException& error)
  {
    LOG_S(ERROR) << error.what();
  }
}

void OutEngineHDF5::write_dot(string _name, Grid<double> data)
{
  try
  {
    Group group = out_file->getGroup('/' + subpath);

    // Create dataspace for the dataset in the file.
    std::vector<size_t> dims{1};

    DataSet dataset = group.createDataSet<double>(_name, DataSpace(dims));

    dataset.write(data(size[2], size[3]));
  }
  catch(FileException& error)
  {
    LOG_S(ERROR) << error.what();
  }

  catch (GroupException& error)
  {
    LOG_S(ERROR) << error.what();
    LOG_S(FATAL) << "Can not write dataset ``" << _name << "'' to group ``/" << subpath << "''";
  }

  // catch failure caused by the DataSet operations
  catch(DataSetException& error)
  {
    LOG_S(ERROR) << error.what();
  }
}

void OutEngineHDF5::write_1d_vector(string _name, vector<double> data)
{
  double* arr = &data[0];

  for (unsigned int i = 0; i < data.size(); ++i)
    arr[i] = data[i];

  try
  {
    Group group = out_file->getGroup('/' + subpath);

    // Create dataspace for the dataset in the file.
    std::vector<size_t> dims{data.size()};

    DataSet dataset = group.createDataSet<double>(_name, DataSpace(dims));

    dataset.write(arr);
  }
  catch(FileException& error)
  {
    LOG_S(ERROR) << error.what();
  }

  catch(GroupException& error)
  {
    LOG_S(ERROR) << error.what();
    LOG_S(FATAL) << "Can not write dataset ``" << _name << "'' to group ``/" << subpath << "''";
  }

  // catch failure caused by the DataSet operations
  catch(DataSetException& error)
  {
    LOG_S(ERROR) << error.what();
  }
}

void OutEngineHDF5::write_metadata(string _metadata)
{
  try
  {
    LOG_S(MAX) << "Creating group ``/metadata''";
    Group group = out_file->createGroup("/metadata");
    LOG_S(MAX) << "Group ``/metadata'' created";

    LOG_S(MAX) << "Writing metadata";
    // DataSpace attr_dataspace = DataSpace(H5S_SCALAR);
    // StrType datatype(PredType::C_S1, _metadata.size());

    Attribute attribute = group.createAttribute<std::string>(
      "metadata", DataSpace::From(_metadata));

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

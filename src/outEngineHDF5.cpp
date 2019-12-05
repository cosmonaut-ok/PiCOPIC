// https://portal.hdfgroup.org/pages/viewpage.action?pageId=50073884

#include <typeinfo>

#include "outEngineHDF5.hpp"

using namespace H5;

OutEngineHDF5::OutEngineHDF5 (string a_path, string a_subpath, unsigned int a_shape, int *a_size,
                                bool a_append, bool a_compress, unsigned int a_compress_level)
  : OutEngine (a_path, a_subpath, a_shape, a_size, a_append, a_compress, a_compress_level)
{
  // make root directory
  lib::make_directory(path);

  Exception::dontPrint();

  // Exception::dontPrint();
  H5File file;
  string fname = path + "/data.h5";

  try
  {
    file = H5File(fname.c_str(), H5F_ACC_RDWR);
  }
  catch (const FileIException&)
  {
    file = H5File(fname.c_str(), H5F_ACC_TRUNC);
  }

  // H5File file(path + "/data.h5", H5F_ACC_TRUNC );

  // create group with subpath name recursively
	std::vector<std::string> group_array;
  const char delim = '/';
  lib::splitstr(subpath, delim, group_array);

  string group_name;
  for (auto i = group_array.begin(); i != group_array.end(); ++i)
  {
    group_name += delim;
    group_name += (*i);
    try
    {
      LOG_DBG("Creating group " << group_name);
      file.createGroup(group_name.c_str());
    }
    catch (const FileIException&)
    {
      LOG_DBG("Group " << group_name << " already exists. Skipping");
    }
  }

  file.close();

  if (compress)
    LOG_WARN("compression is not supported by plaintext output engine");
}

void OutEngineHDF5::write_rec(string a_name, Grid<double> data)
{
  try
  {
    Exception::dontPrint();

    H5File file(path + "/data.h5", H5F_ACC_RDWR );

    Group group(file.openGroup('/' + subpath));

    // Create property list for a dataset and set up fill values.
    int fillvalue = 0;   /* Fill value for the dataset */
    DSetCreatPropList plist;
    plist.setFillValue(PredType::NATIVE_INT, &fillvalue);

    // Create dataspace for the dataset in the file.
    unsigned int size_x = size[2] - size[0];
    unsigned int size_y = size[3] - size[1];
    hsize_t fdim[] = {size_x, size_y};
    DataSpace fspace( 2, fdim );

    DataSet* dataset = new DataSet(group.createDataSet(
                                     a_name,
                                     PredType::NATIVE_DOUBLE, fspace, plist));

    if (size_x == data.x_size && size_y == data.y_size)
      dataset->write(data.get_grid(), PredType::NATIVE_DOUBLE);
    else
    {
      // copy part of the grid to temporary,
      // if size less, than fullframe
      Grid<double> tmp(size_x, size_y, 0);
      for (int i = size[0]; i < size[2]; ++i)
        for (int j = size[1]; j < size[3]; ++j)
          tmp.set(i-size[0], i-size[1], data(i, j));
      dataset->write(tmp.get_grid(), PredType::NATIVE_DOUBLE);
    }

    file.close();
  }
  catch(FileIException error)
  {
    error.printErrorStack();
  }

  // catch failure caused by the DataSet operations
  catch(DataSetIException error)
  {
    error.printErrorStack();
  }
}

void OutEngineHDF5::write_vec(string a_name, Grid<double> data)
{
  MSG_FIXME("OutEngineHDF5::write_vec: is not implemented");
try
  {
    Exception::dontPrint();

    H5File file(path + "/data.h5", H5F_ACC_RDWR );

    Group group(file.openGroup('/' + subpath));

    hsize_t fdim[2];

    // Create property list for a dataset and set up fill values.
    int fillvalue = 0;   /* Fill value for the dataset */
    DSetCreatPropList plist;
    plist.setFillValue(PredType::NATIVE_INT, &fillvalue);

    // vector by Z-component with fixed r (r_begin)
    if (size[1] == -1 && size[0] == -1 && size[3] == -1)
    {
      fdim[0] = data.y_size;

      //   for (unsigned int i = 0; i < data.y_size; i++)
      //     out_val << data(size[0], i) << " ";
      // }
      // // vector by R-component with fixed z (z_begin)
    }
    else if (size[0] == -1 && size[2] == -1 && size[1] == -1)
    {
      fdim[0] = data.x_size;

      //   for (unsigned int i = 0; i < data.x_size; i++)
      //     out_val << data(i, size[1]) << " ";
    }
    else
      LOG_CRIT("Incorrect shape for vector output", 1);

    // Create dataspace for the dataset in the file.

    double arr[fdim[0]];

    int fdd = fdim[0] - 1;

    for (int i = 0; i < fdim[0]; ++i)
      if (fdim[0] == data.y_size)
        arr[i] = data(size[2]-1, i);
      else if (fdim[0] == data.x_size)
        arr[i] = data(i, size[3]-1);

    DataSpace fspace( 1, fdim );

    DataSet* dataset = new DataSet(group.createDataSet(
                                     a_name,
                                     PredType::NATIVE_DOUBLE, fspace, plist));

    dataset->write(arr, PredType::NATIVE_DOUBLE);

    file.close();
  }

  catch(FileIException error)
  {
    error.printErrorStack();
  }

  // catch failure caused by the DataSet operations
  catch(DataSetIException error)
  {
    error.printErrorStack();
  }


  // std::ofstream::openmode omode;
  // if (append)
  //   omode = ios::app;
  // else
  //   omode = ios::trunc;

  // ofstream out_val(path + "/" + subpath + "/" + a_name + ".dat", omode);

  // // vector by Z-component with fixed r (r_begin)
  // if (size[1] == -1 && size[2] == -1 && size[3] == -1)
  // {
  //   for (unsigned int i = 0; i < data.y_size; i++)
  //     out_val << data(size[0], i) << " ";
  // }
  // // vector by R-component with fixed z (z_begin)
  // else if (size[0] == -1 && size[2] == -1 && size[3] == -1)
  // {
  //   for (unsigned int i = 0; i < data.x_size; i++)
  //     out_val << data(i, size[1]) << " ";
  // }
  // else
  //   LOG_CRIT("Incorrect shape for vector output", 1);

  // out_val.close();
}

void OutEngineHDF5::write_dot(string a_name, Grid<double> data)
{
  try
  {
    Exception::dontPrint();

    H5File file(path + "/data.h5", H5F_ACC_RDWR );

    Group group(file.openGroup('/' + subpath));

    // Create property list for a dataset and set up fill values.
    int fillvalue = 0;   /* Fill value for the dataset */
    DSetCreatPropList plist;
    plist.setFillValue(PredType::NATIVE_INT, &fillvalue);

    // Create dataspace for the dataset in the file.
    hsize_t fdim[] = {1};
    DataSpace fspace( 1, fdim );

    DataSet* dataset = new DataSet(group.createDataSet(
                                     a_name,
                                     PredType::NATIVE_DOUBLE, fspace, plist));

    double ds[] = {data(size[2], size[3])};
    dataset->write(ds, PredType::NATIVE_DOUBLE);

    file.close();
  }
  catch(FileIException error)
  {
    error.printErrorStack();
  }

  // catch failure caused by the DataSet operations
  catch(DataSetIException error)
  {
    error.printErrorStack();
  }
}

void OutEngineHDF5::write_1d_vector(string a_name, vector<double> data)
{
  double* arr = &data[0];

  for (unsigned int i = 0; i < data.size(); ++i)
    arr[i] = data[i];

  try
  {
    // Exception::dontPrint();

    H5File file(path + "/data.h5", H5F_ACC_RDWR );

    Group group(file.openGroup('/' + subpath));

    // Create property list for a dataset and set up fill values.
    int fillvalue = 0;   /* Fill value for the dataset */
    DSetCreatPropList plist;
    plist.setFillValue(PredType::NATIVE_INT, &fillvalue);

    // Create dataspace for the dataset in the file.
    hsize_t fdim[] = {data.size()};
    DataSpace fspace( 1, fdim );

    DataSet* dataset = new DataSet(group.createDataSet(
                                     a_name,
                                     PredType::NATIVE_DOUBLE, fspace, plist));

    dataset->write(arr, PredType::NATIVE_DOUBLE);

    file.close();
  }
  catch(FileIException error)
  {
    error.printErrorStack();
  }

  // catch failure caused by the DataSet operations
  catch(DataSetIException error)
  {
    error.printErrorStack();
  }
}
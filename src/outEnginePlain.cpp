#include <typeinfo>

#include "outEnginePlain.hpp"

OutEnginePlain::OutEnginePlain (string a_path, unsigned int a_shape, int *a_size,
                                bool a_append, bool a_compress, unsigned int a_compress_level)
  : OutEngine (a_path, a_shape, a_size, a_append, a_compress, a_compress_level)
{
  // make root directory
  lib::make_directory(path);

  if (compress)
    LOG_WARN("compression is not supported by plaintext output engine");
}

void OutEnginePlain::write_rec(string a_name, Grid<double> data)
{
  std::ofstream::openmode omode;
  if (append)
    omode = ios::app;
  else
    omode = ios::trunc;

  ofstream out_val(path + "/" + a_name + ".dat", omode);

  for (int i = size[0]; i < size[2]; i++)
    for (int j = size[1]; j < size[3]; j++)
      out_val << data(i, j) << " ";

  out_val.close();
}

void OutEnginePlain::write_vec(string a_name, Grid<double> data)
{
  std::ofstream::openmode omode;
  if (append)
    omode = ios::app;
  else
    omode = ios::trunc;

  ofstream out_val(path + "/" + a_name + ".dat", omode);

  // vector by Z-component with fixed r (r_begin)
  if (size[1] == -1 && size[2] == -1 && size[3] == -1)
  {
    for (unsigned int i = 0; i < data.y_size; i++)
      out_val << data(size[0], i) << " ";
  }
  // vector by R-component with fixed z (z_begin)
  else if (size[0] == -1 && size[2] == -1 && size[3] == -1)
  {
    for (unsigned int i = 0; i < data.x_size; i++)
      out_val << data(i, size[1]) << " ";
  }
  else
    LOG_CRIT("Incorrect shape for vector output", 1);

  out_val.close();
}

void OutEnginePlain::write_dot(string a_name, Grid<double> data)
{
    std::ofstream::openmode omode;
  if (append)
    omode = ios::app;
  else
    omode = ios::trunc;

  ofstream out_val(path + "/" + a_name + ".dat", omode);

  out_val << data(size[0], size[1]) << endl;

  out_val.close();
}

void OutEnginePlain::write_1d_vector(string a_name, vector<double> data)
{
  std::ofstream::openmode omode;
  if (append)
    omode = ios::app;
  else
    omode = ios::trunc;

  ofstream out_val(path + "/" + a_name + ".dat", omode);

  for (auto i = data.begin(); i != data.end(); ++i)
    out_val << *i << " ";

  out_val.close();
}

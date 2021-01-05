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

#ifndef _LIB_HPP_
#define _LIB_HPP_

#include <math.h>
#include <vector>
#include <string>
#include <algorithm>  // std::transform
#include <sys/stat.h> // stat
#include <fstream>    // std::ifstream

#ifdef __SSE__
#include <pmmintrin.h>
#endif

#ifdef __AVX__
#include <immintrin.h> //AVX Intrinsic Functions
#endif

#ifdef _WIN32
#include <direct.h> // _mkdir
#endif

#include "defines.hpp"
#include "msg.hpp"

#include "algo/grid.hpp"
#include "constant.hpp"

using namespace std;

namespace algo::common
{
  double sqrt_recip (float x);
  bool to_bool (string str);
  std::string get_simulation_duration ();
  bool directory_exists (const std::string& path);
  bool make_directory (const std::string& path);
  double sq_rt (double x);
  char *get_cmd_option (char **begin, char **end, const std::string &option);
  bool cmd_option_exists (char **begin, char **end, const std::string &option);
  std::vector<double> read_file_to_double (const char *filename);
  void splitstr (std::string const &str, const char delim,
                 std::vector<std::string> &out);
  unsigned int nearest_divide (unsigned int number, double what);
  unsigned short hash_from_string(const std::string& str, unsigned short salt);

  // template is workaround, because c++ does not see Grid class template
  // from inside of custom namespace (for unknown reason)
  template<typename T>
  void bilinear_interpolation(T a, T b)
  {
#pragma omp parallel
    {
      // "cross" linear interpolation
#pragma omp for
      for (unsigned int i = 1; i < a.x_size; i++)
        for (unsigned int j = 1; j < a.y_size; j++)
          b.set(i, j, (a(i-1, j-1) + a(i, j-1) + a(i-1, j) + a(i, j)) / 4);

#pragma omp for
      // process border conditions
      for (unsigned int i = 1; i < a.x_size; i++)
        b.set(i, 0, (a(i-1, 0) + a(i, 0)) / 2);

#pragma omp for
      for (unsigned int j = 1; j < a.y_size; j++)
        b.set(0, j, (a(0, j-1) + a(0, j)) / 2);
    }
  }

  template<typename T>
  void bicubic_interpolation (T a, T b, unsigned int x_size, unsigned int y_size)
  {
    double b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16;

    int x = 0.5;
    int y = 0.5;

#pragma omp parallel for
    for (unsigned int i = 2; i < x_size - 2; i++)
      for (unsigned int j = 2; j < y_size - 2; j++)
      {
        b1 = (x-1) * (x-2) * (x+1) * (y-1) * (y-2) * (y+1) / 4;
        b2 = - x * (x+1) * (x-2) * (y-1) * (y-2) * (y+1) / 4;
        b3 = - (x-1) * (x-2) * (x+1) * y * (y-2) * (y+1) / 4;
        b4 = x * (x-2) * (x+1) * y * (y-2) * (y+1) / 4;
        b5 = - x * (x-1) * (x-2) * (y-1) * (y-2) * (y+1) / 12;
        b6 = - (x-1) * (x-2) * (x+1) * y * (y-2) * (y-1) / 12;
        b7 = x * (x-2) * (x-1) * y * (y-2) * (y+1) / 12;
        b8 = x * (x-2) * (x+1) * y * (y-2) * (y-1) / 12;
        b9 = x * (x-1) * (x+1) * (y-1) * (y-2) * (y+1) / 12;
        b10 = (x-1) * (x-2) * (x+1) * y * (y-1) * (y+1) / 12;
        b11 = x * (x-2) * (x-1) * y * (y-2) * (y-1) / 36;
        b12 = - x * (x-1) * (x+1) * y * (y-2) * (y+1) / 12;
        b13 = - x * (x-2) * (x+1) * y * (y-1) * (y+1) / 12;
        b14 = - x * (x-1) * (x+1) * y * (y-2) * (y-1) / 36;
        b15 = - x * (x-2) * (x-1) * y * (y-1) * (y+1) / 36;
        b16 = x * (x-1) * (x+1) * y * (y-1) * (y+1) / 36;

        b.set(i, j, b1 * a(i, j)
              + b2 * a(i, j+1)
              + b3 * a(i+1, j)
              + b4 * a(i+1, j+1)
              + b5 * a(i, j-1)
              + b6 * a(i-1, j)
              + b7 * a(i+1, j-1)
              + b8 * a(i-1, j+1)
              + b9 * a(i, j+2)
              + b10 * a(i+2, j)
              + b11 * a(i-1, j-1)
              + b12 * a(i+1, j+2)
              + b13 * a(i+2, j+1)
              + b14 * a(i-1, j+2)
              + b15 * a(i+2, j-1)
              + b16 * a(i+2, j+2)
          );
      }
    // TODO: make correct border processing
    for (unsigned int j = 0; j < a.y_size; j++)
    {
      b.set(0, j, 0);
      b.set(1, j, 0);
      b.set(a.x_size, j, 0);
      b.set(a.x_size, j, 0);
    }

    for (unsigned int i = 0; i < a.x_size; i++)
    {
      b.set(i, 0, 0);
      b.set(i, 1, 0);
      b.set(i, a.y_size, 0);
      b.set(i, a.y_size, 0);
    }
  }
}

#endif // end of _LIB_HPP_

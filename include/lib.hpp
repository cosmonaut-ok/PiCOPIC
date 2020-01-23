#pragma once

#include <sstream>
#include <iostream>
#include <string>
#include <sys/stat.h> // stat
#include <errno.h> // errno, ENOENT, EEXIST
#include <math.h>
#include <algorithm>
#include <ctime>

#include <fstream>
#include <vector>
#include <cstdlib>

#ifdef __SSE__
#include <pmmintrin.h>
#endif

#ifdef __AVX__
#include <immintrin.h>//AVX Intrinsic Functions
#endif

#ifdef _WIN32
#include <direct.h> // _mkdir
#endif

#include "constant.hpp"
#include "defines.hpp"
#include "msg.hpp"

using namespace std;

namespace lib
{
  double sqrt_recip(float x);
  bool to_bool(string str);
  double get_gamma (double velocity);
  double get_gamma_inv (double velocity);
  char* get_simulation_duration();
  bool directory_exists(const std::string& path);
  bool make_directory(const std::string& path);
  double sq_rt(double x);
  char *get_cmd_option(char **begin, char **end, const std::string &option);
  bool cmd_option_exists(char **begin, char **end, const std::string &option);
  std::vector<double> read_file_to_double(const char *filename);
  void splitstr(std::string const &str, const char delim,
                std::vector<std::string> &out);
  unsigned int nearest_divide (unsigned int number, double what);
}

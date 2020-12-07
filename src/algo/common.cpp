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

#include "algo/common.hpp"

using namespace constant;

namespace algo::common
{
  // C^2 define c^2 to decrease number of operations
  const double LIGHT_VEL_POW_2 = pow (LIGHT_VEL, 2);

#if __cplusplus >= 201103L
  using std::isnan;
#endif

#ifdef __AVX__
  double sqrt_recip(double x)
  {
    //! 1/sqrt(x), using AVX squared root procedure
    double d[4];
    __m256d c = _mm256_sqrt_pd(_mm256_set_pd(x, 0, 0, 0));

    _mm256_store_pd(d, c);

    return d[3];
  }
#endif

  double sq_rt(double x)
  {
#if defined (__AVX__)
    return sqrt_recip(x);
#else
    return sqrt(x);
#endif
  }

  bool to_bool(string str)
  {
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    std::istringstream is(str);
    bool b;
    is >> std::boolalpha >> b;
    return b;
  }

  std::string get_simulation_duration()
  // get the time, spent since simulation launched (in "d h m s" format)
  {
    double time_sec = std::time(nullptr) - SIMULATION_START_TIME;

    int d, hr, min;
    double sec;
    char* the_time = new char;

    int sec_in_min = 60;
    int sec_in_hr = 3600;
    int sec_in_day = 86400;

    if (time_sec > sec_in_min && time_sec < sec_in_hr)
    {
      min = (int)time_sec / sec_in_min;
      sec = time_sec - min * sec_in_min;
      sprintf(the_time, "%dm %.0fs", min, sec);
    }
    else if (time_sec > sec_in_hr && time_sec < sec_in_day)
    {
      hr = (int)time_sec / sec_in_hr;
      min = (int)((time_sec - hr * sec_in_hr) / sec_in_min);
      sec = time_sec - hr * sec_in_hr - min * sec_in_min;
      sprintf(the_time, "%dh %dm %.0fs", hr, min, sec);
    }
    else if (time_sec > sec_in_day)
    {
      d = (int)time_sec / sec_in_day;
      hr = (int)((time_sec - d * sec_in_day) / sec_in_hr);
      min = (int)((time_sec - d * sec_in_day - hr * sec_in_hr) / sec_in_min);
      sec = time_sec - d * sec_in_day - hr * sec_in_hr - min * sec_in_min;
      sprintf(the_time, "%dd %dh %dm %.0fs", d, hr, min, sec);
    }
    else
      sprintf(the_time, "%.0fs", time_sec);

    return (std::string)the_time;
  }

// Make directory and check if it exists
  bool directory_exists(const std::string& path)
  {
#ifdef _WIN32
    struct _stat info;
    if (_stat(path.c_str(), &info) != 0)
      return false;
    return (info.st_mode & _S_IFDIR) != 0;
#else
    struct stat info;
    if (stat(path.c_str(), &info) != 0)
      return false;
    return (info.st_mode & S_IFDIR) != 0;
#endif
  }

  bool make_directory(const std::string& path)
  {
    LOG_S(MAX) << "Creating directory " << path;

#ifdef _WIN32
    int ret = _mkdir(path.c_str());
#else
    mode_t mode = 0755;
    int ret = mkdir(path.c_str(), mode);
#endif
    if (ret == 0)
      return true;

    switch (errno)
    {
    case ENOENT:
      // parent didn't exist, try to create it
    {
      long unsigned int pos = path.find_last_of('/');
      if (pos == std::string::npos)
#ifdef _WIN32
        pos = path.find_last_of('\\');
      if (pos == std::string::npos)
#endif
        return false;
      if (!make_directory(path.substr(0, pos)))
        return false;
    }
    // now, try to create again
#ifdef _WIN32
    return 0 == _mkdir(path.c_str());
#else
    return 0 == mkdir(path.c_str(), mode);
#endif

    case EEXIST:
      // done!
      return directory_exists(path);

    default:
      return false;
    }
  }

  char *get_cmd_option(char **begin, char **end, const std::string &option)
  {
    char * *itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
      return *itr;
    return 0;
  }

  bool cmd_option_exists(char **begin, char **end, const std::string &option)
  {
    return std::find(begin, end, option) != end;
  }

  std::vector<double> read_file_to_double(const char *filename)
  {
    std::ifstream ifile(filename, std::ios::in);
    std::vector<double> scores;

    //check to see that the file was opened correctly:
    if (!ifile.is_open())
      LOG_S(FATAL) << "Can not open file ";

    double num = 0.0;
    //keep storing values from the text file so long as data exists:
    while (ifile >> num)
      scores.push_back(num);

    //verify that the scores were stored correctly:
    // for (int i = 0; i < scores.size(); ++i) {
    //   std::cout << scores[i] << std::endl;
    // }

    return scores;
  }

  void splitstr(std::string const &str, const char delim,
                std::vector<std::string> &out)
  {
    size_t start;
    size_t end = 0;

    while ((start = str.find_first_not_of(delim, end)) != std::string::npos)
    {
      end = str.find(delim, start);
      out.push_back(str.substr(start, end - start));
    }
  }

  unsigned int nearest_divide (unsigned int number, double what)
  {
    unsigned int res = number;
    while ( res / what != floor(res / what) )
      ++res;

    return res;
  }

  unsigned int hash_from_string (const std::string& str, unsigned int salt)
  {
    unsigned int ret = salt;
    for (std::size_t i = 0; i < str.length(); i++)
      ret ^= ((ret << 5) + str[i] + (ret >> 2));

    return (ret & 0x7FFFFFFF);
  }
}

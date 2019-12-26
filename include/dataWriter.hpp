#pragma once

#include <string>
#include <sstream>

#include "msg.hpp"
#include "lib.hpp"
#include "grid.hpp"
#include "timeSim.hpp"
#include "geometry.hpp"
#include "area.hpp"

#ifdef USE_HDF5
#include "outEngineHDF5.hpp"
#else
#include "outEnginePlain.hpp"
#endif // USE_HDF5


using namespace std;

class Area;

class DataWriter
{
public:
  string path;
  string name;
  string component;
  string specie;
  unsigned int shape;
  unsigned int schedule;
  int size[4];
  bool compress;
  unsigned int compress_level;

  Geometry *geometry;
  TimeSim *time;
  Grid<Area*> areas;

#ifdef USE_HDF5
  OutEngineHDF5 engine;
#else
  OutEnginePlain engine;
#endif // USE_HDF5

private:
  Grid<double> out_data;
  vector<double> out_data_plain;

public:
  DataWriter(string a_path, string a_component,
             string a_specie, unsigned int a_shape,
             int *a_size, unsigned int a_schedule,
             bool a_compress, unsigned int a_compress_level,
             Geometry *a_geom, TimeSim *a_time,
             Grid<Area *> a_areas, string a_metadata);

  void merge_areas(string component, string specie);
  void merge_particle_areas(string parameter, unsigned int component, string specie);
  void go();
};

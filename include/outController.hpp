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

#ifndef _OUT_CONTROLLER_HPP_
#define _OUT_CONTROLLER_HPP_

#include <string>
#include <highfive/H5File.hpp>

#include "defines.hpp"
#include "msg.hpp"

#include "geometry.hpp"
#include "timeSim.hpp"
#include "cfg.hpp"
#include "SMB.hpp"

#ifdef ENABLE_HDF5
#include "outEngine/outEngineHDF5.hpp"
#endif // ENABLE_HDF5

#include "outWriter.hpp"

class OutController
{
public:
  Geometry *geometry;
  TimeSim *time;
  vector<probe> probes;
  SMB *smb;

#ifdef ENABLE_HDF5
  HighFive::File *hdf5_file;
#endif // ENABLE_HDF5

public:
  OutController() {};
#ifdef ENABLE_HDF5
  OutController(HighFive::File *_file, Geometry *_geometry, TimeSim *_time,
                vector<probe> &_probes, SMB *_smb, std::string _metadata,
		bool _print_progress_table);
#else
  OutController(Geometry *_geometry, TimeSim *_time,
                vector<probe> &_probes, SMB *_smb, std::string _metadata,
		bool _print_progress_table);
#endif
  ~OutController() {};

  void operator()(); // launch writers on all domains of all SMBs

private:
  void init_datasets();
  bool print_progress_table;
};

#endif // end of _OUT_CONTROLLER_HPP_

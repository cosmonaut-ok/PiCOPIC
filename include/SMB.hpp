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

#ifndef _SMB_HPP_
#define _SMB_HPP_

#ifdef _OPENMP
#include <omp.h>
// #else
// #define omp_get_thread_num() 0
#endif // _OPENMP

#ifdef ENABLE_MPI
#include <mpi.h>
// using namespace MPI;
#endif // ENABLE_MPI

#include "defines.hpp"
#include "constant.hpp"
#include "msg.hpp"
#include "algo/grid.hpp"
#include "domain.hpp"

#include "cfg.hpp"
#include "timeSim.hpp"

#define BEAM_ID_START 1000

#ifdef ENABLE_MPI
#if SIZE_MAX == UCHAR_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
   #define MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
   #error "what is happening here?"
#endif

//// MPI comm overlay struct
struct MPICommOverlay
{
  size_t len;

  std::vector<double> col_0_0;
  std::vector<double> col_0_1;
  std::vector<double> col_0_2;
  std::vector<double> col_1_0;
  std::vector<double> col_1_1;
  std::vector<double> col_1_2;
  std::vector<double> col_2_0;
  std::vector<double> col_2_1;
  std::vector<double> col_2_2;
  std::vector<double> col_3_0;
  std::vector<double> col_3_1;
  std::vector<double> col_3_2;

  std::vector<int> dst_domain;

  MPICommOverlay() : dtype(MPI_DATATYPE_NULL) {}
  MPICommOverlay(size_t _len, unsigned int _r, unsigned int _z) : dtype(MPI_DATATYPE_NULL)
  {
    len = _len;

    dst_domain.push_back(_r);
    dst_domain.push_back(_z);

    for (unsigned int i = 0; i < len; ++i)
    {
      col_0_0.push_back(0);
      col_0_1.push_back(0);
      col_0_2.push_back(0);
      col_1_0.push_back(0);
      col_1_1.push_back(0);
      col_1_2.push_back(0);
      col_2_0.push_back(0);
      col_2_1.push_back(0);
      col_2_2.push_back(0);
      col_3_0.push_back(0);
      col_3_1.push_back(0);
      col_3_2.push_back(0);
    }
  };
  ~MPICommOverlay()
  {
    if (dtype != MPI_DATATYPE_NULL)
      MPI_Type_free(&dtype);
  };

  const MPI_Datatype mpi_dtype()
  {
    if (dtype == MPI_DATATYPE_NULL)
      make_dtype();
    return dtype;
  };

  void invalidate_dtype()
  {
    if (dtype != MPI_DATATYPE_NULL)
      MPI_Type_free(&dtype);
  }

  void make_dtype()
  {
    const int nblock = 14;
    const int block_count[nblock] =
      {
        1,
        (int)col_0_0.size(),
        (int)col_0_1.size(),
        (int)col_0_2.size(),
        (int)col_1_0.size(),
        (int)col_1_1.size(),
        (int)col_1_2.size(),
        (int)col_2_0.size(),
        (int)col_2_1.size(),
        (int)col_2_2.size(),
        (int)col_3_0.size(),
        (int)col_3_1.size(),
        (int)col_3_2.size(),
        (int)dst_domain.size()
      };

    MPI_Datatype block_type[nblock] =
      {
        MPI_INT,
        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
        MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
        MPI_INT
      };
    MPI_Aint offset[nblock];
    MPI_Get_address(&len, &offset[0]);
    MPI_Get_address(&col_0_0[0], &offset[1]);
    MPI_Get_address(&col_0_1[0], &offset[2]);
    MPI_Get_address(&col_0_2[0], &offset[3]);
    MPI_Get_address(&col_1_0[0], &offset[4]);
    MPI_Get_address(&col_1_1[0], &offset[5]);
    MPI_Get_address(&col_1_2[0], &offset[6]);
    MPI_Get_address(&col_2_0[0], &offset[7]);
    MPI_Get_address(&col_2_1[0], &offset[8]);
    MPI_Get_address(&col_2_2[0], &offset[9]);
    MPI_Get_address(&col_3_0[0], &offset[10]);
    MPI_Get_address(&col_3_1[0], &offset[11]);
    MPI_Get_address(&col_3_2[0], &offset[12]);
    MPI_Get_address(&dst_domain[0], &offset[13]);

    MPI_Type_create_struct(nblock, block_count, offset, block_type, &dtype);
    MPI_Type_commit(&dtype);
  };

private:
  MPI_Datatype dtype;

};
//// /MPI comm overlay struct
#endif // ENABLE_MPI

class SMB
{
public:
  Grid<Domain*> domains;
  unsigned int r_domains;
  unsigned int z_domains;

private:
  Geometry *geometry;
  Cfg *cfg;
  TimeSim *time;
  int world_rank;
  int world_size;

public:
  SMB ( void ) {};
  SMB ( Cfg* _cfg, Geometry *_geometry, TimeSim *_time,
        int _world_rank, int _world_size);

  ~SMB( void ) {};

// private:
  void particles_runaway_collector ();
  void current_overlay ();
  void field_h_overlay ();
  void field_e_overlay ();

  // public:
  void solve_maxvell();
  void solve_current();
  void advance_particles();
  void inject_beam();
  void distribute();
};
#endif // end of _SMB_HPP_

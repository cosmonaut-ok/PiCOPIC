#include <gtest/gtest.h>
#include "grid3d.hpp"

namespace {
  int deref(int * pint)
  {
    return *pint;
  }

#define VALUE 8.3

  // TEST(grid3d, get_grid)
  // {
  //   Grid<double> grid (10, 10, 3);
  //   grid = VALUE;
  //   double **g_grd = grid.get_grid();

  //   ASSERT_DOUBLE_EQ(g_grd[12][12], VALUE);
  //   ASSERT_NE(g_grd[15][15], VALUE);
  // }

  TEST(grid3d, constructor)
  {
    Grid3D<double> grid (10, 10, 2);
    grid = VALUE;

    ASSERT_EQ(grid[0](0, 0), VALUE);
    ASSERT_EQ(grid[2](9, 9), VALUE);
    ASSERT_NE(grid[1].get_grid()[13][13], VALUE);
    ASSERT_EXIT((deref(nullptr),grid[3](10, 10)),::testing::KilledBySignal(SIGSEGV),".*");
  }

//   TEST(grid3d, set)
//   {
//     Grid<double> grid (10, 10, 2);
//     grid = 0;
//     grid.set(3, 4, VALUE);

//     ASSERT_EQ(grid(3, 4), VALUE);
//     for (unsigned int i = 0; i < 10; ++i)
//       for (unsigned int j = 0; j < 10; ++j)
//         if (! (i == 3 && j == 4))
//           ASSERT_NE(grid(i, j), VALUE);
//   }

//   TEST(grid3d, inc)
//   {
//     Grid<double> grid (10, 10, 2);
//     grid = 0;
//     grid.inc(3, 4, VALUE);
//     grid.inc(3, 4, VALUE);

//     ASSERT_EQ(grid(3, 4), VALUE * 2);
//     for (unsigned int i = 0; i < 10; ++i)
//       for (unsigned int j = 0; j < 10; ++j)
//         if (! (i == 3 && j == 4))
//           ASSERT_NE(grid(i, j), VALUE * 2);
//   }

//   TEST(grid3d, dec)
//   {
//     Grid<double> grid (10, 10, 2);
//     grid = VALUE;
//     grid.dec(3, 4, VALUE);
//     grid.dec(3, 4, VALUE);

//     ASSERT_EQ(grid(3, 4), - VALUE);
//     for (unsigned int i = 0; i < 10; ++i)
//       for (unsigned int j = 0; j < 10; ++j)
//         if (! (i == 3 && j == 4))
//           ASSERT_NE(grid(i, j), - VALUE);
//   }


//   TEST(grid3d, m_a)
//   {
//     Grid<double> grid (10, 10, 2);
//     grid = VALUE;
//     grid.m_a(3, 4, VALUE);
//     grid.m_a(3, 4, VALUE);

// #define MUA 571.787

//     ASSERT_DOUBLE_EQ(grid(3, 4), MUA);
//     for (unsigned int i = 0; i < 10; ++i)
//       for (unsigned int j = 0; j < 10; ++j)
//         if (! (i == 3 && j == 4))
//           ASSERT_NE(grid(i, j), MUA);
//   }


//   TEST(grid3d, d_a)
//   {
//     Grid<double> grid (10, 10, 2);
//     grid = VALUE;
//     grid.d_a(3, 4, VALUE);
//     grid.d_a(3, 4, VALUE);

// #define DIVA 0.12048192771084336

//     ASSERT_DOUBLE_EQ(grid(3, 4), DIVA);
//     for (unsigned int i = 0; i < 10; ++i)
//       for (unsigned int j = 0; j < 10; ++j)
//         if (! (i == 3 && j == 4))
//           ASSERT_NE(grid(i, j), DIVA);
//   }

  TEST(grid3d, overlay_set)
  {
    for (unsigned int i_shift = 0; i_shift < 4; ++i_shift)
    {
      Grid3D<double> grid (10, 10, i_shift);
      grid = VALUE;
      grid.overlay_set(VALUE * 2);

      double **g_grd_0 = grid[0].get_grid();
      double **g_grd_1 = grid[1].get_grid();
      double **g_grd_2 = grid[2].get_grid();

      for (unsigned int i = 0; i < 10 + i_shift; ++i)
        for (unsigned int j = 0; j < 10 + i_shift; ++j)
          for (unsigned int k = 0; k < 2; ++k)
            if (i < i_shift || j < i_shift || i >= 10 + i_shift  || j >= 10 + i_shift)
            {
              ASSERT_EQ(g_grd_0[i][j], VALUE * 2);
              ASSERT_EQ(g_grd_1[i][j], VALUE * 2);
              ASSERT_EQ(g_grd_2[i][j], VALUE * 2);
            }
            else
            {
              ASSERT_EQ(g_grd_0[i][j], VALUE);
              ASSERT_EQ(g_grd_1[i][j], VALUE);
              ASSERT_EQ(g_grd_2[i][j], VALUE);
            }
    }
  }

  TEST(grid3d, overlay_x)
  {
    Grid3D<double> grid (10, 10, 3);
    Grid3D<double> grid2 (10, 10, 3);
    grid = VALUE;
    grid.overlay_set(VALUE);
    grid2 = 20.3;
    grid2.overlay_set(20.3);

    grid.overlay_x(grid2);

    double **g_grd_0 = grid[0].get_grid();
    double **g_grd_1 = grid[1].get_grid();
    double **g_grd_2 = grid[2].get_grid();

    double **g_grd2_0 = grid2[0].get_grid();
    double **g_grd2_1 = grid2[1].get_grid();
    double **g_grd2_2 = grid2[2].get_grid();

    ASSERT_EQ(g_grd_0[15][4], g_grd2_0[0][4]);
    ASSERT_EQ(g_grd_1[14][5], g_grd2_1[1][5]);
    ASSERT_EQ(g_grd_2[13][6], g_grd2_2[2][6]);

    ASSERT_NE(g_grd_0[0][4], g_grd2_0[15][4]);
    ASSERT_NE(g_grd_1[1][5], g_grd2_1[14][5]);
    ASSERT_NE(g_grd_2[2][6], g_grd2_2[13][6]);

    ASSERT_EQ(g_grd_0[15][4], 28.6);
    ASSERT_EQ(g_grd_1[14][5], 28.6);
    ASSERT_EQ(g_grd_2[13][6], 28.6);

    ASSERT_NE(g_grd2_0[15][4], 28.6);
    ASSERT_NE(g_grd2_1[14][5], 28.6);
    ASSERT_NE(g_grd2_2[13][6], 28.6);
  }

  TEST(grid3d, overlay_y)
  {
    Grid3D<double> grid (10, 10, 3);
    Grid3D<double> grid2 (10, 10, 3);
    grid = VALUE;
    grid.overlay_set(VALUE);
    grid2 = 20.3;
    grid2.overlay_set(20.3);

    grid.overlay_y(grid2);

    double **g_grd_0 = grid[0].get_grid();
    double **g_grd_1 = grid[1].get_grid();
    double **g_grd_2 = grid[2].get_grid();

    double **g_grd2_0 = grid2[0].get_grid();
    double **g_grd2_1 = grid2[1].get_grid();
    double **g_grd2_2 = grid2[2].get_grid();

    ASSERT_EQ(g_grd_0[4][15], g_grd2_0[4][0]);
    ASSERT_EQ(g_grd_1[5][14], g_grd2_1[5][1]);
    ASSERT_EQ(g_grd_2[6][13], g_grd2_2[6][2]);

    ASSERT_NE(g_grd_0[4][0], g_grd2_0[4][15]);
    ASSERT_NE(g_grd_1[5][1], g_grd2_1[5][14]);
    ASSERT_NE(g_grd_2[6][2], g_grd2_2[6][13]);

    ASSERT_EQ(g_grd_0[4][15], 28.6);
    ASSERT_EQ(g_grd_1[5][14], 28.6);
    ASSERT_EQ(g_grd_2[6][13], 28.6);

    ASSERT_NE(g_grd2_0[4][15], 28.6);
    ASSERT_NE(g_grd2_1[5][14], 28.6);
    ASSERT_NE(g_grd2_2[6][13], 28.6);
  }

  TEST(grid3d, overlay_xy)
  {
    Grid3D<double> grid (10, 10, 3);
    Grid3D<double> grid2 (10, 10, 3);
    grid = VALUE;
    grid.overlay_set(VALUE);
    grid2 = 20.3;
    grid2.overlay_set(20.3);

    grid.overlay_xy(grid2);

    double **g_grd_0 = grid[0].get_grid();
    double **g_grd_1 = grid[1].get_grid();
    double **g_grd_2 = grid[2].get_grid();

    double **g_grd2_0 = grid2[0].get_grid();
    double **g_grd2_1 = grid2[1].get_grid();
    double **g_grd2_2 = grid2[2].get_grid();

    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
      {
        ASSERT_EQ(g_grd_0[13+i][13+j], g_grd2_0[i][j]);
        ASSERT_EQ(g_grd_1[13+i][13+j], g_grd2_1[i][j]);
        ASSERT_EQ(g_grd_2[13+i][13+j], g_grd2_2[i][j]);

        ASSERT_NE(g_grd2_0[13+i][13+j], g_grd_0[i][j]);
        ASSERT_NE(g_grd2_1[13+i][13+j], g_grd_1[i][j]);
        ASSERT_NE(g_grd2_2[13+i][13+j], g_grd_2[i][j]);

        ASSERT_EQ(g_grd_0[13+i][13+j], 28.6);
        ASSERT_EQ(g_grd_1[13+i][13+j], 28.6);
        ASSERT_EQ(g_grd_2[13+i][13+j], 28.6);

        ASSERT_EQ(g_grd2_0[i][j], 28.6);
        ASSERT_EQ(g_grd2_1[i][j], 28.6);
        ASSERT_EQ(g_grd2_2[i][j], 28.6);
      }
  }

  TEST(grid3d, _operator_assign)
  {
    Grid3D<double> grid (10, 10, 3);
    grid = VALUE;

    for (unsigned int k = 0; k < 3; ++k)
      for (unsigned int i = 0; i < 10; ++i)
        for (unsigned int j = 0; j < 10; ++j)
          ASSERT_DOUBLE_EQ(grid(k, i, j), VALUE);
  }

  TEST(grid3d, _operator_parenthesis)
  {
    Grid3D<double> grid (10, 10, 3);
    grid = VALUE;

    ASSERT_DOUBLE_EQ(grid(0, 1, 2), VALUE);
    // ASSERT_EXIT(grid(3, 1, 2),::testing::ExitedWithCode(1),".*");
    ASSERT_EXIT((deref(nullptr),grid(0, 10, 10)), ::testing::KilledBySignal(SIGSEGV),".*");
  }

  TEST(grid3d, _operator_squared_braces)
  {
    Grid3D<double> grid (10, 10, 3);
    grid = VALUE;

    ASSERT_DOUBLE_EQ(grid[0](1, 2), VALUE);
    // ASSERT_EXIT(grid[3],::testing::ExitedWithCode(1),".*");
  }
}

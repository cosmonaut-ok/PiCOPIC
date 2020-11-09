#define LOGURU_WITH_STREAMS 1

#include <gtest/gtest.h>
#include "algo/grid.hpp"
  
namespace {
  int deref(int * pint)
  {
    return *pint;
  }

#define VALUE 8.3


  TEST(grid, get_grid)
  {
    Grid<double> grid (10, 10, 3);
    grid = VALUE;
    double **g_grd = grid.get_grid();

    ASSERT_DOUBLE_EQ(g_grd[12][12], VALUE);
    ASSERT_NE(g_grd[15][15], VALUE);
  }

  TEST(grid, constructor)
  {
    Grid<double> grid (10, 10, 2);
    grid = VALUE;

    ASSERT_EQ(grid(0, 0), VALUE);
    ASSERT_EQ(grid(9, 9), VALUE);
    ASSERT_NE(grid.get_grid()[13][13], VALUE);
    ASSERT_EXIT((deref(nullptr),grid(10, 10)),::testing::KilledBySignal(SIGSEGV),".*");
  }

  TEST(grid, set)
  {
    Grid<double> grid (10, 10, 2);
    grid = 0;
    grid.set(3, 4, VALUE);

    ASSERT_EQ(grid(3, 4), VALUE);
    for (unsigned int i = 0; i < 10; ++i)
      for (unsigned int j = 0; j < 10; ++j)
        if (i != 3 && j != 4)
          { ASSERT_NE(grid(i, j), VALUE); }
  }

  TEST(grid, inc)
  {
    Grid<double> grid (10, 10, 2);
    grid = 0;
    grid.inc(3, 4, VALUE);
    grid.inc(3, 4, VALUE);

    ASSERT_EQ(grid(3, 4), VALUE * 2);
    for (unsigned int i = 0; i < 10; ++i)
      for (unsigned int j = 0; j < 10; ++j)
        if (i != 3 && j != 4)
          { ASSERT_NE(grid(i, j), VALUE * 2); }
  }

  TEST(grid, dec)
  {
    Grid<double> grid (10, 10, 2);
    grid = VALUE;
    grid.dec(3, 4, VALUE);
    grid.dec(3, 4, VALUE);

    ASSERT_EQ(grid(3, 4), - VALUE);
    for (unsigned int i = 0; i < 10; ++i)
      for (unsigned int j = 0; j < 10; ++j)
        if (i != 3 && j != 4)
	  { ASSERT_NE(grid(i, j), - VALUE); }
  }


  TEST(grid, m_a)
  {
    Grid<double> grid (10, 10, 2);
    grid = VALUE;
    grid.m_a(3, 4, VALUE);
    grid.m_a(3, 4, VALUE);

#define MUA 571.787

    ASSERT_DOUBLE_EQ(grid(3, 4), MUA);
    for (unsigned int i = 0; i < 10; ++i)
      for (unsigned int j = 0; j < 10; ++j)
        if (i != 3 && j != 4)
          { ASSERT_NE(grid(i, j), MUA); }
  }


  TEST(grid, d_a)
  {
    Grid<double> grid (10, 10, 2);
    grid = VALUE;
    grid.d_a(3, 4, VALUE);
    grid.d_a(3, 4, VALUE);

#define DIVA 0.12048192771084336

    ASSERT_DOUBLE_EQ(grid(3, 4), DIVA);
    for (unsigned int i = 0; i < 10; ++i)
      for (unsigned int j = 0; j < 10; ++j)
        if (i != 3 && j != 4)
          { ASSERT_NE(grid(i, j), DIVA); }
  }

  TEST(grid, overlay_set)
  {
    for (unsigned int i_shift = 0; i_shift < 4; ++i_shift)
    {
      Grid<double> grid (10, 10, i_shift);
      grid = VALUE;
      grid.overlay_set(VALUE * 2);

      double **g_grd = grid.get_grid();

      for (unsigned int i = 0; i < 10 + i_shift; ++i)
        for (unsigned int j = 0; j < 10 + i_shift; ++j)
          if (i < i_shift || j < i_shift || i >= 10 + i_shift  || j >= 10 + i_shift)
            ASSERT_EQ(g_grd[i][j], VALUE * 2);
          else
            ASSERT_EQ(g_grd[i][j], VALUE);
    }
  }

  TEST(grid, overlay_x)
  {
    Grid<double> grid (10, 10, 3);
    Grid<double> grid2 (10, 10, 3);
    grid = VALUE;
    grid.overlay_set(VALUE);
    grid2 = 20.3;
    grid2.overlay_set(20.3);

    grid.overlay_x(grid2);

    double **g_grd = grid.get_grid();
    double **g_grd2 = grid2.get_grid();

    ASSERT_EQ(g_grd[15][4], g_grd2[0][4]);
    ASSERT_EQ(g_grd[14][5], g_grd2[1][5]);
    ASSERT_EQ(g_grd[13][6], g_grd2[2][6]);

    ASSERT_NE(g_grd[0][4], g_grd2[15][4]);
    ASSERT_NE(g_grd[1][5], g_grd2[14][5]);
    ASSERT_NE(g_grd[2][6], g_grd2[13][6]);

    ASSERT_EQ(g_grd[15][4], 28.6);
    ASSERT_EQ(g_grd[14][5], 28.6);
    ASSERT_EQ(g_grd[13][6], 28.6);

    ASSERT_NE(g_grd2[15][4], 28.6);
    ASSERT_NE(g_grd2[14][5], 28.6);
    ASSERT_NE(g_grd2[13][6], 28.6);
  }

  TEST(grid, overlay_y)
  {
    Grid<double> grid (10, 10, 3);
    Grid<double> grid2 (10, 10, 3);
    grid = VALUE;
    grid.overlay_set(VALUE);
    grid2 = 20.3;
    grid2.overlay_set(20.3);

    grid.overlay_y(grid2);

    double **g_grd = grid.get_grid();
    double **g_grd2 = grid2.get_grid();

    ASSERT_EQ(g_grd[4][15], g_grd2[4][0]);
    ASSERT_EQ(g_grd[5][14], g_grd2[5][1]);
    ASSERT_EQ(g_grd[6][13], g_grd2[6][2]);

    ASSERT_NE(g_grd[4][0], g_grd2[4][15]);
    ASSERT_NE(g_grd[5][1], g_grd2[5][14]);
    ASSERT_NE(g_grd[6][2], g_grd2[6][13]);

    ASSERT_EQ(g_grd[4][15], 28.6);
    ASSERT_EQ(g_grd[5][14], 28.6);
    ASSERT_EQ(g_grd[6][13], 28.6);

    ASSERT_NE(g_grd2[4][15], 28.6);
    ASSERT_NE(g_grd2[5][14], 28.6);
    ASSERT_NE(g_grd2[6][13], 28.6);
  }

  TEST(grid, overlay_xy)
  {
    Grid<double> grid (10, 10, 3);
    Grid<double> grid2 (10, 10, 3);
    grid = VALUE;
    grid.overlay_set(VALUE);
    grid2 = 20.3;
    grid2.overlay_set(20.3);

    grid.overlay_xy(grid2);

    double **g_grd = grid.get_grid();
    double **g_grd2 = grid2.get_grid();

    for (unsigned int i = 0; i < 3; ++i)
      for (unsigned int j = 0; j < 3; ++j)
      {
        ASSERT_EQ(g_grd[13+i][13+j], g_grd2[i][j]);
        ASSERT_NE(g_grd2[13+i][13+j], g_grd[i][j]);

        ASSERT_EQ(g_grd[13+i][13+j], 28.6);
        ASSERT_EQ(g_grd2[i][j], 28.6);
      }
  }

  TEST(grid, copy)
  {
    Grid<double> grid (10, 10, 3);
    Grid<double> grid2 (10, 10, 3);

    grid = VALUE;
    grid2.copy(grid);

    ASSERT_EQ(grid2(3, 4), VALUE);
  }

  TEST(grid, _operator_parenthesis)
  {
    Grid<double> grid (10, 10, 3);
    grid = VALUE;
    grid.set(3, 4, 8.1);

    for (unsigned int i = 0; i < 10; ++i)
      for (unsigned int j = 0; j < 10; ++j)
        if (i == 3 && j == 4)
          ASSERT_EQ(grid(i, j), 8.1);
        else
          ASSERT_EQ(grid(i, j), VALUE);
  }

  TEST(grid, _operator_inc)
  {
    Grid<double> grid (10, 10, 3);
    grid = VALUE;
    grid += 1.2;

    for (unsigned int i = 0; i < 10; ++i)
      for (unsigned int j = 0; j < 10; ++j)
        ASSERT_DOUBLE_EQ(grid(i, j), 9.5);
  }

  TEST(grid, _operator_dec)
  {
    Grid<double> grid (10, 10, 3);
    grid = VALUE;
    grid -= 1.2;

    for (unsigned int i = 0; i < 10; ++i)
      for (unsigned int j = 0; j < 10; ++j)
        ASSERT_DOUBLE_EQ(grid(i, j), 7.1);
  }

  TEST(grid, _operator_m_a)
  {
    Grid<double> grid (10, 10, 3);
    grid = VALUE;
    grid *= 1.2;

    for (unsigned int i = 0; i < 10; ++i)
      for (unsigned int j = 0; j < 10; ++j)
        ASSERT_DOUBLE_EQ(grid(i, j), 9.96);
  }

  TEST(grid, _operator_d_a)
  {
    Grid<double> grid (10, 10, 3);
    grid = VALUE;
    grid /= 1.2;

    for (unsigned int i = 0; i < 10; ++i)
      for (unsigned int j = 0; j < 10; ++j)
        ASSERT_DOUBLE_EQ(grid(i, j), 6.916666666666668);
  }
}

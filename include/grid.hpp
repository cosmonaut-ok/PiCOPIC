#pragma once

//! grid class implented to take into account
//! overlay areas. They takes first 2 elements
//! and last 2 elements of the grid

#include "msg.hpp"

template <class T>
class Grid
{
  T **grid;

private:
  int x_size;
  int y_size;

  int x_real_size;
  int y_real_size;
  int o_s; // overlay shift

public:
  Grid() {};
  Grid(int x_amount, int y_amount)
  {
    o_s = 2;

    x_size = x_amount;
    y_size = y_amount;
    x_real_size = x_amount + 2 * o_s; // first 2 and last 2 elements
    y_real_size = y_amount + 2 * o_s;

    grid = new T *[x_real_size];
    for (int i = 0; i < x_real_size; i++)
      grid[i] = new T[y_real_size];

    // set all values to 0
    for (int i = 0; i < x_real_size; ++i)
      for (int j = 0; j < y_real_size; ++j)
        grid[i][j] = 0;
  };

  void set(int x, int y, T value)
  {
    grid[x+o_s][y+o_s] = value;
  };

  void inc(int x, int y, T value)
  {
    grid[x+o_s][y+o_s] += value;
  };

  void dec(int x, int y, T value)
  {
    grid[x+o_s][y+o_s] -= value;
  };

  void d_a(int x, int y, T value)
  // Division assignment
  {
    grid[x+o_s][y+o_s] /= value;
  };

  void m_a(int x, int y, T value)
  // Multiplication assignment
  {
    grid[x+o_s][y+o_s] *= value;
  };

  int size_x()
  {
    return x_size;
  };

  int size_y()
  {
    return y_size;
  };

  // operators overloading
  T& operator() (int x, int y)
  {
    return grid[x+o_s][y+o_s];
  }

  // operatros for update all of the elements
  // of the grid
  Grid& operator= (T value)&
  {
    for (int i = 0; i < x_size; ++i)
      for (int j = 0; j < y_size; ++j)
        grid[i+o_s][j+o_s] = value;

    return *this;
  };

  Grid& operator+= (T value)&
  {
    for (int i = 0; i < x_size; ++i)
      for (int j = 0; j < y_size; ++j)
        grid[i+o_s][j+o_s] += value;

    return *this;
  };

  Grid& operator-= (T value)&
  {
    for (int i = 0; i < x_size; ++i)
      for (int j = 0; j < y_size; ++j)
        grid[i+o_s][j+o_s] -= value;

    return *this;
  };

  Grid& operator*= (T value)&
  {
    for (int i = 0; i < x_size; ++i)
      for (int j = 0; j < y_size; ++j)
        grid[i+o_s][j+o_s] *= value;

    return *this;
  };

  Grid& operator/= (T value)&
  {
    for (int i = 0; i < x_size; ++i)
      for (int j = 0; j < y_size; ++j)
        grid[i+o_s][j+o_s] /= value;

    return *this;
  };

  // overlay-related methods
  void overlay_reset()
  {
    for (int i = 0; i < x_real_size; ++i)
    {
      grid[i][0] = 0;
      grid[i][1] = 0;
      grid[i][y_real_size-1] = 0;
      grid[i][y_real_size-2] = 0;
    }

    for (int j = 0; j < y_real_size; ++j)
    {
      grid[0][j] = 0;
      grid[1][j] = 0;
      grid[x_real_size-1][j] = 0;
      grid[x_real_size-2][j] = 0;
    }
  };

  void overlay_right(Grid<T> rgrid)
  {
    if (x_size == rgrid.size_x())
      for (int i = 2; i < x_real_size - 2; ++i)
      {
        rgrid.inc(i, 1, grid[i][y_real_size-1]);
        rgrid.inc(i, 0, grid[i][y_real_size-2]);

        grid[i][y_real_size-4] += rgrid(i, -2);
        grid[i][y_real_size-3] += rgrid(i, -1);

        rgrid.set(i, -2, grid[i][y_real_size-4]);
        rgrid.set(i, -1, grid[i][y_real_size-3]);

        grid[i][y_real_size-2] = rgrid(i, 0);
        grid[i][y_real_size-1] = rgrid(i, 1);
      }
    else
    {
      LOG_CRIT("overlay_right: X sizes of left and right grid are not equal. Can not overlay", 1);
    }
  };

  void overlay_top(Grid<T> tgrid)
  {
    if (y_size == tgrid.size_y())
      for (int i = 2; i < y_real_size - 2; ++i)
      {
        tgrid.inc(1, i, grid[x_real_size-1][i]);
        tgrid.inc(0, i, grid[x_real_size-2][i]);

        grid[x_real_size-4][i] += tgrid(-2, i);
        grid[x_real_size-3][i] += tgrid(-1, i);

        tgrid.set(-2, i, grid[x_real_size-4][i]);
        tgrid.set(-1, i, grid[x_real_size-3][i]);

        grid[x_real_size-2][i] = tgrid(0, i);
        grid[x_real_size-1][i] = tgrid(1, i);
      }
    else
    {
      LOG_CRIT("overlay_top: Y sizes of bottom and top grid are not equal. Can not overlay", 1);
    }
  };

  void overlay_top_right(Grid<T> trgrid)
  {
    trgrid.inc(1, 1, grid[x_real_size-1][y_real_size-1]);
    trgrid.inc(1, 0, grid[x_real_size-2][y_real_size-1]);
    trgrid.inc(0, 1, grid[x_real_size-1][y_real_size-2]);
    trgrid.inc(0, 0, grid[x_real_size-2][y_real_size-2]);

    grid[x_size-3][y_real_size-3] += trgrid(1, 1);
    grid[x_size-3][y_real_size-4] += trgrid(1, 0);
    grid[x_size-4][y_real_size-3] += trgrid(0, 1);
    grid[x_size-4][y_real_size-4] += trgrid(0, 0);

    trgrid.set(-2, -2, grid[x_real_size-1][y_real_size-1]);
    trgrid.set(-2, -1, grid[x_real_size-1][y_real_size-2]);
    trgrid.set(-1, -2, grid[x_real_size-2][y_real_size-1]);
    trgrid.set(-1, -1, grid[x_real_size-2][y_real_size-2]);

    grid[x_real_size-1][y_real_size-1] += trgrid(1, 1);
    grid[x_real_size-1][y_real_size-2] += trgrid(1, 0);
    grid[x_real_size-2][y_real_size-1] += trgrid(0, 1);
    grid[x_real_size-2][y_real_size-2] += trgrid(0, 0);
  };

};

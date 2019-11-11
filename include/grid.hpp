#pragma once

#include "msg.hpp"

template <class T>
class Grid
{
  T **grid;

private:
  unsigned int x;
  unsigned int y;

public:
  Grid() {};
  Grid(unsigned int x_amount, unsigned int y_amount)
  {
    x = x_amount;
    y = y_amount;
    grid = new T *[x_amount];
    for (unsigned int i = 0; i < x_amount; i++)
      grid[i] = new T[y_amount];
  };

  void set(unsigned int x, unsigned int y, T value)
  {
    grid[x][y] = value;
  };

  void inc(unsigned int x, unsigned int y, T value)
  {
    grid[x][y] += value;
  };

  void dec(unsigned int x, unsigned int y, T value)
  {
    grid[x][y] -= value;
  };

  void d_a(unsigned int x, unsigned int y, T value)
  // Division assignment
  {
    grid[x][y] /= value;
  };

  void m_a(unsigned int x, unsigned int y, T value)
  // Multiplication assignment
  {
    grid[x][y] *= value;
  };

  unsigned int size_x()
  {
    return x;
  };

  unsigned int size_y()
  {
    return y;
  };

  void overlay_right(Grid<T> rgrid)
  {
    if (x == rgrid.size_x())
      for (unsigned int i = 0; i < x; ++i)
      {
        rgrid.inc(i, 0,  grid[i][y-2]);
        rgrid.inc(i, 1, grid[i][y-1]);

        // grid[i][y-2] = rgrid(i, 0);
        // grid[i][y-1] = rgrid(i, 1);
      }
    else
    {
      LOG_CRIT("overlay_right: X sizes of left and right grid are not equal. Can not overlay", 1);
    }
  };

  void overlay_top(Grid<T> tgrid)
  {
    if (y == tgrid.size_y())
      for (unsigned int i = 0; i < y; ++i)
      {
        tgrid.inc(0, i, grid[x-2][i]);
        tgrid.inc(1, i, grid[x-1][i]);

        // grid[x-2][i] = tgrid(0, i);
        // grid[x-1][i] = tgrid(1, i);
      }
    else
    {
      LOG_CRIT("overlay_top: Y sizes of bottom and top grid are not equal. Can not overlay", 1);
    }
  };

  void overlay_top_right(Grid<T> trgrid)
  {
    trgrid.inc(0, 0, grid[x-2][y-2]);
    trgrid.inc(0, 1, grid[x-2][y-1]);
    trgrid.inc(1, 0, grid[x-1][y-2]);
    trgrid.inc(1, 1, grid[x-1][y-1]);
  };

  // operators overloading
  T& operator() (unsigned int x, unsigned int y)
  {
    return grid[x][y];
  }

  // operatros for update all of the elements
  // of the grid
  Grid& operator= (T value)&
  {
    for (unsigned int i = 0; i < x; ++i)
      for (unsigned int j = 0; j < y; ++j)
        grid[i][j] = value;

    return *this;
  };

  Grid& operator+= (T value)&
  {
    for (unsigned int i = 0; i < x; ++i)
      for (unsigned int j = 0; j < y; ++j)
        grid[i][j] += value;

    return *this;
  };

  Grid& operator-= (T value)&
  {
    for (unsigned int i = 0; i < x; ++i)
      for (unsigned int j = 0; j < y; ++j)
        grid[i][j] -= value;

    return *this;
  };

  Grid& operator*= (T value)&
  {
    for (unsigned int i = 0; i < x; ++i)
      for (unsigned int j = 0; j < y; ++j)
        grid[i][j] *= value;

    return *this;
  };

  Grid& operator/= (T value)&
  {
    for (unsigned int i = 0; i < x; ++i)
      for (unsigned int j = 0; j < y; ++j)
        grid[i][j] /= value;

    return *this;
  };
};

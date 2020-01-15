#pragma once

#include "msg.hpp"

template <class T>
class Grid
{
  T **grid;

private:
  unsigned int o_s; // overlay shift

  unsigned int x_real_size;
  unsigned int y_real_size;

public:
  unsigned int x_size;
  unsigned int y_size;

public:
  Grid() {};
  Grid(unsigned int x_amount, unsigned int y_amount, unsigned int overlay_shift)
  {
    o_s = overlay_shift;
    x_real_size = x_amount + o_s * 2;
    y_real_size = y_amount + o_s * 2;

    x_size = x_amount;
    y_size = y_amount;

    grid = new T *[x_real_size];
    for (unsigned int i = 0; i < x_real_size; i++)
      grid[i] = new T[y_real_size];
  };

  void set(unsigned int x, unsigned int y, T value)
  {
    grid[x+o_s][y+o_s] = value;
  };

  void inc(unsigned int x, unsigned int y, T value)
  {
    grid[x+o_s][y+o_s] += value;
  };

  void dec(unsigned int x, unsigned int y, T value)
  {
    grid[x+o_s][y+o_s] -= value;
  };

  void d_a(unsigned int x, unsigned int y, T value)
  // Division assignment
  {
    grid[x+o_s][y+o_s] /= value;
  };

  void m_a(unsigned int x, unsigned int y, T value)
  // Multiplication assignment
  {
    grid[x+o_s][y+o_s] *= value;
  };

  void overlay_set(T value)
  {
    // unsigned int j = 0;
    for (unsigned int i = 0; i < x_real_size; ++i)
      for (unsigned int j = 0; j < y_real_size; ++j)
        if (i < o_s || j < o_s
            ||
            i >= x_real_size - o_s || j >= y_real_size - o_s
          )
          grid[i][j] = value;
  };

  void overlay_y(Grid<T> rhsgrid)
  {
    if (x_real_size == rhsgrid.x_real_size)
      for (unsigned int d = 0; d < o_s; ++d)
        for (unsigned int i = o_s; i < x_real_size - o_s; ++i)
        {
          rhsgrid.grid[i][2*o_s-d-1] += grid[i][y_real_size-1-d];
          grid[i][y_real_size-2*o_s+d] += rhsgrid.grid[i][d];

          rhsgrid.grid[i][d] = grid[i][y_real_size-2*o_s+d];
          grid[i][y_real_size-1-d] = rhsgrid.grid[i][2*o_s-1-d];
        }
    else
    {
      LOG_CRIT("overlay_x: X sizes of left and right grid are not equal. Can not overlay", 1);
    }
  };

  void overlay_x(Grid<T> rhsgrid)
  {
    if (y_real_size == rhsgrid.y_real_size)
      for (unsigned int d = 0; d < o_s; ++d)
        for (unsigned int i = o_s; i < y_real_size - o_s; ++i)
        {
          rhsgrid.grid[2*o_s-d-1][i] += grid[x_real_size-1-d][i];
          grid[x_real_size-2*o_s+d][i] += rhsgrid.grid[d][i];

          rhsgrid.grid[d][i] = grid[x_real_size-2*o_s+d][i];
          grid[x_real_size-1-d][i] = rhsgrid.grid[2*o_s-1-d][i];
        }
    else
    {
      LOG_CRIT("overlay_y: Y sizes of bottom and top grid are not equal. Can not overlay", 1);
    }
  };

  void overlay_xy(Grid<T> rhsgrid)
  {
    if (x_real_size == rhsgrid.x_real_size && y_real_size == rhsgrid.y_real_size)
      for (unsigned int d = 0; d < o_s; ++d)
        for (unsigned int e = 0; e < o_s; ++e)
        {
          rhsgrid.grid[2*o_s-d-1][2*o_s-e-1] += grid[x_real_size-1-d][x_real_size-1-e];
          grid[x_real_size-2*o_s+d][x_real_size-2*o_s+e] += rhsgrid.grid[d][e];

          rhsgrid.grid[d][e] = grid[y_real_size-2*o_s+d][y_real_size-2*o_s+e];
          grid[x_real_size-1-d][x_real_size-1-e] = rhsgrid.grid[2*o_s-1-d][2*o_s-1-e];
        }
    else
    {
      LOG_CRIT("overlay_xy: X or Y sizes of bottom-left and top-right grid are not equal. Can not overlay", 1);
    }
  };

  void copy(Grid<T> rhs)
  // fully copy rhs grid to current grid element-by-element
  {
    if (x_real_size == rhs.x_real_size && y_real_size == rhs.y_real_size)
      for (unsigned int i = 0; i < x_real_size; ++i)
        for (unsigned int j = 0; j < y_real_size; ++j)
          grid[i][j] = rhs.grid[i][j];
    else
      LOG_CRIT("overlay_xy: X or Y sizes of bottom-left and top-right grid are not equal. Can not overlay", 1);
  };

  // operators overloading
  T& operator() (unsigned int x, unsigned int y)
  {
    return grid[x+o_s][y+o_s];
  }

  // operatros for update all of the elements
  // of the grid
  Grid& operator= (T value)&
  {
    for (unsigned int i = 0; i < x_size; ++i)
      for (unsigned int j = 0; j < y_size; ++j)
        grid[i+o_s][j+o_s] = value;

    return *this;
  };

  Grid& operator+= (T value)&
  {
    for (unsigned int i = 0; i < x_size; ++i)
      for (unsigned int j = 0; j < y_size; ++j)
        grid[i+o_s][j+o_s] += value;

    return *this;
  };

  Grid& operator-= (T value)&
  {
    for (unsigned int i = 0; i < x_size; ++i)
      for (unsigned int j = 0; j < y_size; ++j)
        grid[i+o_s][j+o_s] -= value;

    return *this;
  };

  Grid& operator*= (T value)&
  {
    for (unsigned int i = 0; i < x_size; ++i)
      for (unsigned int j = 0; j < y_size; ++j)
        grid[i+o_s][j+o_s] *= value;

    return *this;
  };

  Grid& operator/= (T value)&
  {
    for (unsigned int i = 0; i < x_size; ++i)
      for (unsigned int j = 0; j < y_size; ++j)
        grid[i+o_s][j+o_s] /= value;

    return *this;
  };

  T** get_grid()
  {
    return grid;
  }
};

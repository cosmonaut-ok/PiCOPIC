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

  T get(unsigned int x, unsigned int y)
  {
    return grid[x][y];
  };

  void set(unsigned int x, unsigned int y, T value)
  {
    grid[x][y] = value;
  };

  void inc(unsigned int x, unsigned int y, T value)
  {
    grid[x][y] += value;
  };

  void setall(T value)
  {
    for (unsigned int i = 0; i < x; ++i)
      for (unsigned int j = 0; j < y; ++j)
        grid[i][j] = value;
  };

  int size_x()
  {
    return x;
  };

  int size_y()
  {
    return y;
  };
};

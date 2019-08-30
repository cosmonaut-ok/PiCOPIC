#pragma once

#include "msg.hpp"
#include "grid.hpp"

template<class T>
class Grid3D
{
public:
  Grid <T> r_component;
  Grid <T> phi_component;
  Grid <T> z_component;

public:
  Grid3D() {};
  Grid3D(unsigned int r_amount, unsigned int z_amount)
  {
    r_component = Grid<T> (r_amount, z_amount);
    phi_component = Grid<T> (r_amount, z_amount);
    z_component = Grid<T> (r_amount, z_amount);
  };

  void setall(T value)
  {
    r_component.setall(value);
    phi_component.setall(value);
    z_component.setall(value);
  };

  Grid<T>& operator[] (int index)
  {
    // Grid <T> ret;
    if (index == 0) return r_component;
    if (index == 1) return phi_component;
    if (index == 2) return z_component;
    LOG_CRIT("Grid3D: Out of index", 1);
  };

};

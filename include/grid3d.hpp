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
  Grid3D(unsigned int r_amount, unsigned int z_amount, unsigned int overlay_shift)
  {
    r_component = Grid<T> (r_amount, z_amount, overlay_shift);
    phi_component = Grid<T> (r_amount, z_amount, overlay_shift);
    z_component = Grid<T> (r_amount, z_amount, overlay_shift);
  };

  void overlay_x(Grid3D<T> rgrid)
  {
    r_component.overlay_x(rgrid.r_component);
    phi_component.overlay_x(rgrid.phi_component);
    z_component.overlay_x(rgrid.z_component);
  };

  void overlay_y(Grid3D<T> rgrid)
  {
    r_component.overlay_y(rgrid.r_component);
    phi_component.overlay_y(rgrid.phi_component);
    z_component.overlay_y(rgrid.z_component);
  };

  void overlay_xy(Grid3D<T> trgrid)
  {
    r_component.overlay_y(trgrid.r_component);
    phi_component.overlay_y(trgrid.phi_component);
    z_component.overlay_y(trgrid.z_component);
  };

  void overlay_set(T value)
  {
    r_component.overlay_set(value);
    phi_component.overlay_set(value);
    z_component.overlay_set(value);
  };

  Grid3D<T>& operator= (const T value)
  {
    r_component = value;
    phi_component = value;
    z_component = value;

    return *this;
  };

  Grid<T>& operator[] (int index)
  {
    // Grid <T> ret;
    if (index == 0) return r_component;
    if (index == 1) return phi_component;
    if (index == 2) return z_component;
    LOG_CRIT("Grid3D: Out of index", 1);
  };

  T& operator() (unsigned int x, unsigned int y, unsigned int z)
  {
    if (x == 0) return r_component(y, z);
    if (x == 1) return phi_component(y, z);
    if (x == 2) return z_component(y, z);
    LOG_CRIT("Grid3D: Out of index", 1);
  };
};

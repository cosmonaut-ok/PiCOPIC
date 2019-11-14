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

  void overlay_right(Grid3D<T> rgrid)
  {
    r_component.overlay_right(rgrid.r_component);
    phi_component.overlay_right(rgrid.phi_component);
    z_component.overlay_right(rgrid.z_component);
  };

  void overlay_top(Grid3D<T> rgrid)
  {
    r_component.overlay_top(rgrid.r_component);
    phi_component.overlay_top(rgrid.phi_component);
    z_component.overlay_top(rgrid.z_component);
  };

  void overlay_top_right(Grid3D<T> trgrid)
  {
    r_component.overlay_top(trgrid.r_component);
    phi_component.overlay_top(trgrid.phi_component);
    z_component.overlay_top(trgrid.z_component);
  };

  void overlay_reset()
  {
    r_component.overlay_reset();
    phi_component.overlay_reset();
    z_component.overlay_reset();
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

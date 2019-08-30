#pragma once
#include <ctime>

// WARNING!
// This header file used only to store
// physical and mathematical constants

namespace constant
{
  const double PI = 3.1415926535897932;
  // Vacuum permittivity (electric constant), F*m(e-1)
  const double EPSILON0 = 8.85E-12;
  // Electron mass, kg
  const double EL_MASS = 9.1E-31;
  // Electron charge, coulon
  const double EL_CHARGE = 1.6E-19;
  // light speed in vacuum m/s
  const double LIGHT_VEL = 3.0E8;
  // Vacuum permeability (magnetic constant), m*kg*s(e-2)*A(e-2)
  const double MAGN_CONST = 1.26E-6;
  // simulation start time
  const clock_t SIMULATION_START_TIME = std::time(nullptr);
  // minimal possible distance or velocity.
  // Smaller values should be rounded to zero
  const double MNZL = 1e-15; // Minimal Non-Zeroing Limit
}

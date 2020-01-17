#include "geometry.hpp"

// Geometry constructor
Geometry::Geometry (double rs, double zs,
                    int bot_ngr, int top_ngr, int left_ngz, int right_ngz,
                    double pml_l_z0, double pml_l_zwall,
                    double pml_l_rwall, double pml_sigma1,
                    double pml_sigma2,
                    bool wall_r0, bool wall_z0,
                    bool wall_rr,bool wall_zz)
{
  r_size = rs;
  z_size = zs;
  r_grid_amount = top_ngr - bot_ngr;
  z_grid_amount = right_ngz - left_ngz;

  top_r_grid_number = top_ngr;
  bottom_r_grid_number = bot_ngr;
  left_z_grid_number = left_ngz;
  right_z_grid_number = right_ngz;

  r_cell_size = r_size / r_grid_amount;
  z_cell_size = z_size / z_grid_amount;

  // initialize walls
  walls[0] = wall_r0;
  walls[1] = wall_z0;
  walls[2] = wall_rr;
  walls[3] = wall_zz;

  pml_length[0] = 0; // r=0
  pml_length[1] = pml_l_z0; // z=0
  pml_length[2] = pml_l_rwall; // r=wall
  pml_length[3] = pml_l_zwall; // z=wall

  pml_sigma[0] = pml_sigma1;
  pml_sigma[1] = pml_sigma2;
}

Geometry::~Geometry() {}

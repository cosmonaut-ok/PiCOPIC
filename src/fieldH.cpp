#include "fieldH.hpp"

FieldH::FieldH(Geometry *geom, TimeSim *t, vector<SpecieP *> species) : Field(geom, t)
{
  // ! we use ``field'' member variable to mark magnetic field at half of timestep
  // ! and ``field_at_et'' to mark field at round timesteps
  // ! (same timestamp, as for electrical field)
  // ! it useful to have single implementation of ``get_field'' method
  // ! for both classes - FieldE and FieldH
  field_at_et = Grid3D<double> (geometry->r_grid_amount, geometry->z_grid_amount, 2);
  species_p = species;

  field_at_et = 0;
  field_at_et.overlay_set(0);
}

// Field calculation
void FieldH::calc_field_cylindrical()
{
  field.overlay_set(0);
  field_at_et.overlay_set(0);

  Grid3D<double> el_field = field_e->field;
  double dr = geometry->r_cell_size;
  double dz = geometry->z_cell_size;

  // H_r on outer wall (r=r)
  if (geometry->walls[2])
    for(int k = 0; k < geometry->z_grid_amount; k++)
    {
      int i = geometry->r_grid_amount - 1;
      // alpha constant and delta_t production (to optimize calculations)
      double alpha_t = time->step
        * (el_field(1, i, k + 1) - el_field(1, i, k)) / (dz * MAGN_CONST);

      field[0].set(i, k, field_at_et(0, i, k) + alpha_t / 2);
      field_at_et[0].inc(i, k, alpha_t);
    }

  // regular case
  for (int i = 0; i < geometry->r_grid_amount; i++)
    for (int k = 0; k < geometry->z_grid_amount; k++)
    {
      double alpha_t_r = time->step
        * (el_field(1, i, k + 1) - el_field(1, i, k)) / (dz * MAGN_CONST);

      field[0].set(i, k, field_at_et(0, i, k) + alpha_t_r / 2);
      field_at_et[0].inc(i, k, alpha_t_r);

      double alpha_t_phi = time->step
        * ((el_field(2, i+1, k) - el_field(2, i, k)) / dr
           - (el_field(0, i, k+1) - el_field(0, i, k)) / dz)
        / MAGN_CONST;

      field[1].set(i, k, field_at_et(1, i, k) + alpha_t_phi / 2);
      field_at_et[1].inc(i, k, alpha_t_phi);

      double alpha_t_z = time->step
        * (
          (el_field(1, i+1, k) + el_field(1, i, k)) / (2. * dr * (i + 0.5 + geometry->bottom_r_grid_number))
          + (el_field(1, i+1, k) - el_field(1, i, k)) / dr)
        / MAGN_CONST;

      field[2].set(i, k, field_at_et(2, i, k) - alpha_t_z / 2);
      field_at_et[2].dec(i, k, alpha_t_z);
    }
}

vector3d<double> FieldH::get_field(double radius, double longitude)
//! function for magnetic field weighting
{
  vector3d<double> cmp(0., 0., 0.);
  int i_r = 0; // number of particle i cell
  int k_z = 0; // number of particle k cell
  int i_r_shift = 0;
  int k_z_shift = 0;

  double dr = geometry->r_cell_size;
  double dz = geometry->z_cell_size;
  double r1, r2, r3; // temp variables for calculation
  double dz1, dz2; // temp var.: width of k and k+1 cell
  double vol_1 = 0; // volume of i cell; Q/V, V - volume of elementary cell
  double vol_2 = 0; // volume of i+1 cell;

  r1 = radius - 0.5 * dr;
  r3 = radius + 0.5 * dr;

  //// weighting of H_z
  //finding number of cell. example dr=0.5, radius = 0.7, i_r =0;!!
  i_r = CELL_NUMBER(radius - 0.5 * dr, dr);
  k_z = CELL_NUMBER(longitude, dz);
  i_r_shift = i_r - geometry->bottom_r_grid_number;
  k_z_shift = k_z - geometry->left_z_grid_number;
  // TODO: workaround: sometimes it gives -1.
  // Just get 0 cell if it happence
  if (i_r < 0) i_r = 0;
  if (k_z < 0) k_z = 0;
  if (i_r_shift < 0) i_r_shift = 0;
  if (k_z_shift < 0) k_z_shift = 0;

  // FIXME: it can be more, than field.size_x - 2
  // for some unknown reason
  if (i_r_shift > field[0].x_size)
  {
    MSG_FIXME("fieldH::get_field: i_r_shift is more, than field[0].x_size. Applying workaround");
    i_r_shift = field[0].x_size;
  }

  if (k_z_shift > field[0].y_size)
  {
    MSG_FIXME("fieldH::get_field: k_z_shift is more, than field[0].y_size . Applying workaround");
    k_z_shift = field[0].y_size;
  }

  vol_1 = CELL_VOLUME(i_r+1, dr, dz);
  vol_2 = CELL_VOLUME(i_r+3, dr, dz);
  dz1 = (k_z + 1) * dz - longitude;
  dz2 = longitude - k_z * dz;
  r2 = (i_r + 1) * dr;

  //weighting Hz[i][k]//
  cmp[2] += field(2, i_r_shift, k_z_shift) * CYL_RNG_VOL(dz1, r1, r2) / vol_1;

  //weighting Hz[i+1][k]//
  cmp[2] += field(2, i_r_shift + 1, k_z_shift) * CYL_RNG_VOL(dz1, r2, r3) / vol_2;

  //weighting Hz[i][k+1]//
  cmp[2] += field(2, i_r_shift, k_z_shift + 1) * CYL_RNG_VOL(dz2, r1, r2) / vol_1;

  //weighting Hz[i+1][k+1]//
  cmp[2] += field(2, i_r_shift + 1, k_z_shift + 1) * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

  //// weighting of Hr
  // finding number of cell. example dz=0.5, longitude = 0.7, z_k =0;!!
  i_r = CELL_NUMBER(radius, dr);
  k_z = CELL_NUMBER(longitude - 0.5 * dz, dz);
  i_r_shift = i_r - geometry->bottom_r_grid_number;
  k_z_shift = k_z - geometry->left_z_grid_number;
  // TODO: workaround: sometimes it gives -1.
  // Just get 0 cell if it happence
  if (i_r < 0) i_r = 0;
  if (k_z < 0) k_z = 0;
  if (i_r_shift < 0) i_r_shift = 0;
  if (k_z_shift < 0) k_z_shift = 0;

  // FIXME: it can be more, than field.size_x - 2
  // for some unknown reason
  if (i_r_shift > field[0].x_size)
  {
    MSG_FIXME("fieldH::get_field: i_r_shift is more, than field[0].x_size. Applying workaround");
    i_r_shift = field[0].x_size;
  }

  if (k_z_shift > field[0].y_size)
  {
    MSG_FIXME("fieldH::get_field: k_z_shift is more, than field[0].y_size . Applying workaround");
    k_z_shift = field[0].y_size;
  }

  if(radius > dr)
    vol_1 = CELL_VOLUME(i_r, dr, dz);
  else
    vol_1 = CYL_VOL(dz, dr); // volume of first cell

  r2 = (i_r + 0.5) * dr;

  vol_2 = CELL_VOLUME(i_r+2, dr, dz);
  dz1 = (k_z + 1.5) * dz - longitude;
  dz2 = longitude - (k_z + 0.5) * dz;

  //weighting Hr[i][k]//
  cmp[0] += field(0, i_r_shift, k_z_shift) * CYL_RNG_VOL(dz1, r1, r2) / vol_1;

  //weighting Hr[i+1][k]//
  cmp[0] += field(0, i_r_shift + 1, k_z_shift) * CYL_RNG_VOL(dz1, r2, r3) / vol_2;

  //weighting Hr[i][k+1]//
  cmp[0] += field(0, i_r_shift, k_z_shift + 1) * CYL_RNG_VOL(dz2, r1, r2) / vol_1;

  //weighting Hr[i+1][k+1]//
  cmp[0] += field(0, i_r_shift + 1, k_z_shift + 1) * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

  //// weighting of H_fi
  // finding number of cell. example dz=0.5, longitude = 0.7, z_k =0;
  i_r = CELL_NUMBER(radius - 0.5 * dr, dr);
  k_z = CELL_NUMBER(longitude - 0.5 * dz, dz);
  i_r_shift = i_r - geometry->bottom_r_grid_number;
  k_z_shift = k_z - geometry->left_z_grid_number;
  // TODO: workaround: sometimes it gives -1.
  // Just get 0 cell if it happence
  if (i_r < 0) i_r = 0;
  if (k_z < 0) k_z = 0;
  if (i_r_shift < 0) i_r_shift = 0;
  if (k_z_shift < 0) k_z_shift = 0;

  // FIXME: it can be more, than field.size_x - 2
  // for some unknown reason
  if (i_r_shift > field[0].x_size)
  {
    MSG_FIXME("fieldH::get_field: i_r_shift is more, than field[0].x_size. Applying workaround");
    i_r_shift = field[0].x_size;
  }

  if (k_z_shift > field[0].y_size)
  {
    MSG_FIXME("fieldH::get_field: k_z_shift is more, than field[0].y_size . Applying workaround");
    k_z_shift = field[0].y_size;
  }

  r2 = (i_r+1) * dr;
  vol_1 = CELL_VOLUME(i_r + 1, dr, dz);
  vol_2 = CELL_VOLUME(i_r + 3, dr, dz);
  dz1 = (k_z+1.5) * dz - longitude;
  dz2 = longitude - (k_z+0.5) * dz;

  //weighting Hphi[i][k]//
  cmp[1] += field(1, i_r_shift, k_z_shift) * CYL_RNG_VOL(dz1, r1, r2) / vol_1;

  //weighting Hphi[i+1][k]//
  cmp[1] += field(1, i_r_shift + 1, k_z_shift) * CYL_RNG_VOL(dz1, r2, r3) / vol_2;

  //weighting Hphi[i][k+1]//
  cmp[1] += field(1, i_r_shift, k_z_shift + 1) * CYL_RNG_VOL(dz2, r1, r2) / vol_1;

  //weighting Hphi[i+1][k+1]//
  cmp[1] += field(1, i_r_shift + 1, k_z_shift + 1) * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

  return cmp;
}

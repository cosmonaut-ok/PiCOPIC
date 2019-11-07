#include "fieldH.hpp"

FieldH::FieldH(Geometry *geom, TimeSim *t, vector<SpecieP *> species) : Field(geom, t)
{
  // ! we use ``field'' member variable to mark magnetic field at half of timestep
  // ! and ``field_at_et'' to mark field at round timesteps
  // ! (same timestamp, as for electrical field)
  // ! it useful to have single implementation of ``get_field'' method
  // ! for both classes - FieldE and FieldH
  field_at_et = Grid3D<double> (geometry->r_grid_amount + 2, geometry->z_grid_amount + 2);
  species_p = species;

  field_at_et.setall(0);
}

// Field calculation
void FieldH::calc_field_cylindrical() // EField *e_field1, Time *time1
{
  Grid<double> e_r = field_e->field[0];
  Grid<double> e_phi = field_e->field[1];
  Grid<double> e_z = field_e->field[2];
  double dr = geometry->r_cell_size;
  double dz = geometry->z_cell_size;

  // double alpha;

  // Hr - last i value
#pragma omp parallel
  {
#pragma omp for
    for(int k = 0; k<(geometry->z_grid_amount - 1); k++)
    {
      int i=geometry->r_grid_amount - 1;
      // alpha constant and delta_t production (to optimize calculations)
      double alpha_t = time->step
        * (e_phi.get(i, k+1) - e_phi.get(i, k)) / (dz * MAGN_CONST);

      field[0].set(i, k, field_at_et[0].get(i, k) + alpha_t / 2);
      field_at_et[0].set(i, k, field_at_et[0].get(i, k) + alpha_t);
    }

// #pragma omp for
    for(int i = 0; i < (geometry->r_grid_amount - 1); i++)
      for(int k = 0; k < (geometry->z_grid_amount - 1); k++)
      {
        double alpha_t = time->step
          * (e_phi.get(i, k+1) - e_phi.get(i, k)) / (dz * MAGN_CONST);

        field[0].set(i, k, field_at_et[0].get(i, k) + alpha_t / 2);
        field_at_et[0].set(i, k, field_at_et[0].get(i, k) + alpha_t);
      }

// #pragma omp for
    for(int i = 0; i < (geometry->r_grid_amount - 1); i++)
      for(int k = 0; k < (geometry->z_grid_amount - 1); k++)
      {
        double alpha_t = time->step
          * ((e_z.get(i+1, k) - e_z.get(i, k)) / dr
             - (e_r.get(i, k+1) - e_r.get(i, k)) / dz)
          / MAGN_CONST;

        field[1].set(i, k, field_at_et[1].get(i, k) + alpha_t / 2);
        field_at_et[1].set(i, k, field_at_et[1].get(i, k) + alpha_t);
      }

// #pragma omp for
    for(int i = 0; i < (geometry->r_grid_amount - 1); i++)
      for(int k = 0; k < (geometry->z_grid_amount - 1); k++)
      {
        double alpha_t = time->step
          * (
            (e_phi.get(i+1, k) + e_phi.get(i, k)) / (2. * dr * (i + 0.5))
            + (e_phi.get(i+1, k) - e_phi.get(i, k)) / dr)
          / MAGN_CONST;

        field[2].set(i, k, field_at_et[2].get(i, k) - alpha_t / 2);
        field_at_et[2].set(i, k, field_at_et[2].get(i, k) - alpha_t);
      }
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

  vol_1 = CELL_VOLUME(i_r+1, dr, dz);
  vol_2 = CELL_VOLUME(i_r+3, dr, dz);
  dz1 = (k_z + 1) * dz - longitude;
  dz2 = longitude - k_z * dz;
  r2 = (i_r + 1) * dr;

  //weighting Hz[i][k]//
  cmp[2] += field[2].get(i_r_shift, k_z_shift) * CYL_RNG_VOL(dz1, r1, r2) / vol_1;
  
  //weighting Hz[i+1][k]//
  cmp[2] += field[2].get(i_r_shift + 1, k_z_shift) * CYL_RNG_VOL(dz1, r2, r3) / vol_2;
  
  //weighting Hz[i][k+1]//
  cmp[2] += field[2].get(i_r_shift, k_z_shift + 1) * CYL_RNG_VOL(dz2, r1, r2) / vol_1;
  
  //weighting Hz[i+1][k+1]//
  cmp[2] += field[2].get(i_r_shift + 1, k_z_shift + 1) * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

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

  if(radius>dr)
    vol_1 = CELL_VOLUME(i_r, dr, dz);
  else
    vol_1 = CYL_VOL(dz, dr); // volume of first cell

  r2 = (i_r + 0.5) * dr;

  vol_2 = CELL_VOLUME(i_r+2, dr, dz);
  dz1 = (k_z + 1.5) * dz - longitude;
  dz2 = longitude - (k_z + 0.5) * dz;

  //weighting Hr[i][k]//
  cmp[0] += field[0].get(i_r_shift, k_z_shift) * CYL_RNG_VOL(dz1, r1, r2) / vol_1;
  
  //weighting Hr[i+1][k]//
  cmp[0] += field[0].get(i_r_shift + 1, k_z_shift) * CYL_RNG_VOL(dz1, r2, r3) / vol_2;
  
  //weighting Hr[i][k+1]//
  cmp[0] += field[0].get(i_r_shift, k_z_shift + 1) * CYL_RNG_VOL(dz2, r1, r2) / vol_1;
  
  //weighting Hr[i+1][k+1]//
  cmp[0] += field[0].get(i_r_shift + 1, k_z_shift + 1) * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

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
  
  r2 = (i_r+1) * dr;
  vol_1 = CELL_VOLUME(i_r+1, dr, dz);
  vol_2 = CELL_VOLUME(i_r+3, dr, dz);
  dz1 = (k_z+1.5) * dz - longitude;
  dz2 = longitude - (k_z+0.5) * dz;

  //weighting Hphi[i][k]//
  cmp[1] += field[1].get(i_r_shift, k_z_shift) * CYL_RNG_VOL(dz1, r1, r2) / vol_1;
  
  //weighting Hphi[i+1][k]//
  cmp[1] += field[1].get(i_r_shift + 1, k_z_shift) * CYL_RNG_VOL(dz1, r2, r3) / vol_2;
  
  //weighting Hphi[i][k+1]//
  cmp[1] += field[1].get(i_r_shift, k_z_shift + 1) * CYL_RNG_VOL(dz2, r1, r2) / vol_1;
  
  //weighting Hphi[i+1][k+1]//
  cmp[1] += field[1].get(i_r_shift + 1, k_z_shift + 1) * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

  return cmp;
}

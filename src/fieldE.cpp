#include "fieldE.hpp"

FieldE::FieldE(Geometry *geom, TimeSim *t, vector<SpecieP *> species) : Field(geom, t)
{
  species_p = species;
  epsilon = Grid<double> (geometry->r_grid_amount + 4, geometry->z_grid_amount + 4);
  sigma = Grid<double> (geometry->r_grid_amount + 4, geometry->z_grid_amount + 4);

  epsilon = EPSILON0;
  sigma = 0;

  set_pml();
}

// PML
void FieldE::set_pml()
{
  // defining lenght of sigma calculation ion left wall
  double lenght_sigma_left = geometry->z_cell_size
    * (floor(geometry->z_grid_amount * geometry->pml_length[1]));

  // defining lenght of sigma calculation on right wall
  double lenght_sigma_right = geometry->z_cell_size
    * (floor(geometry->z_grid_amount * geometry->pml_length[3]));

  // defining lenght of sigma calculation on z-wall
  double lenght_sigma_extern = geometry->r_cell_size
    * (floor(geometry->r_grid_amount * geometry->pml_length[2]));

  double radial_shift = geometry->bottom_r_grid_number * geometry->r_cell_size;
  double longitudinal_shift = geometry->left_z_grid_number * geometry->z_cell_size;

  double global_r_grid_amount = geometry->r_grid_amount * geometry->areas_by_r;
  double global_z_grid_amount = geometry->z_grid_amount * geometry->areas_by_z;

  // r=r wall
  if (geometry->pml_length[2] != 0)
    for(int i = 0; i < geometry->r_grid_amount; i++)
      for(int k = 0; k < geometry->z_grid_amount; k++)
        if ((geometry->areas_by_r * geometry->r_grid_amount - geometry->bottom_r_grid_number - i) * geometry->r_cell_size
            < lenght_sigma_extern)
          sigma.inc(i, k, geometry->pml_sigma[0] +
                    (geometry->pml_sigma[1] - geometry->pml_sigma[0])
                    / pow(lenght_sigma_extern, 2)
                    * pow(geometry->r_cell_size * (i + 1) + radial_shift
                          - geometry->r_size * geometry->areas_by_r
                          + lenght_sigma_extern, 2));

  // sigma on z=0 and z=z walls
  if (geometry->pml_length[1] != 0)
    for(int k = 0; k < geometry->z_grid_amount; k++)
      for(int i = 0; i < geometry->r_grid_amount; i++)
        // z=0 wall
        if (geometry->z_cell_size * k + longitudinal_shift < lenght_sigma_left)
          sigma.inc(i, k, geometry->pml_sigma[0]
                    + (geometry->pml_sigma[1] - geometry->pml_sigma[0])
                    / pow(lenght_sigma_left, 2)
                    * pow(lenght_sigma_left - geometry->z_cell_size * k - longitudinal_shift, 2));

  // z=z wall
  if (geometry->pml_length[3] != 0)

    for(int k = 0; k < geometry->z_grid_amount; k++)
      for(int i = 0; i < geometry->r_grid_amount; i++)
        if ((geometry->areas_by_z * geometry->z_grid_amount - geometry->left_z_grid_number - k) * geometry->z_cell_size
            < lenght_sigma_right)
          sigma.inc(i, k, geometry->pml_sigma[0]
                    + (geometry->pml_sigma[1] - geometry->pml_sigma[0])
                    / pow(lenght_sigma_right, 2)
                    * pow(geometry->z_cell_size * (k + 1) + longitudinal_shift
                          - geometry->z_size * geometry->areas_by_z
                          + lenght_sigma_right, 2));
}

// Electric field calculation
void FieldE::calc_field_cylindrical()
{
  field.overlay_reset();

  Grid3D<double> curr = current->current;
  Grid3D<double> magn_fld = field_h->field_at_et;

  double dr = geometry->r_cell_size;
  double dz = geometry->z_cell_size;

  // E at the center axis (r=0) case
  if (geometry->walls[0]) // calculate only at the center axis (r=0)
    for(int k = z_begin; k < geometry->z_grid_amount; k++)
    {
      int i = r_begin;
      double epsilonx2 = 2 * epsilon(i, k);
      double sigma_t = sigma(i, k) * time->step;

      double koef_e = (epsilonx2 - sigma_t) / (epsilonx2 + sigma_t);
      double koef_h =  2 * time->step / (epsilonx2 + sigma_t);

      // E_r at the center axis (r=0)
      field[0].m_a(i, k, koef_e);
      field[0].dec(i, k, (curr(0, i, k) + (magn_fld(1, i, k) - magn_fld(1, i, k-1)) / dz) * koef_h);

      // E_z at the center on axis
      field[2].m_a(i, k, koef_e);
      field[2].dec(i, k, (curr(2, i, k) - magn_fld(1, i, k) / dr) * koef_h);
    }

// regular case
  for(int i = r_begin; i < r_end; i++) // TODO: it should be r_begin, instead of 1
    for(int k = z_begin; k < z_end; k++) // TODO: it should be z_begin, instead of 1
    {
      double epsilonx2 = 2 * epsilon(i, k);
      double sigma_t = sigma(i, k) * time->step;

      double koef_e = (epsilonx2 - sigma_t) / (epsilonx2 + sigma_t);
      double koef_h = 2 * time->step / (epsilonx2 + sigma_t);

      field[0].m_a(i, k, koef_e);
      field[0].dec(i, k, (curr(0, i, k) + (magn_fld(1, i, k) - magn_fld(1, i, k-1)) / dz) * koef_h);

      field[1].m_a(i, k, koef_e);
      field[1].dec(i, k, (curr(1, i, k) - (magn_fld(0, i, k) - magn_fld(0, i, k - 1))
                          / dz + (magn_fld(2, i, k) - magn_fld(2, i - 1, k)) / dr) * koef_h);

      field[2].m_a(i, k, koef_e);
      field[2].dec(i, k, (curr(2, i, k) - (magn_fld(1, i, k) - magn_fld(1, i - 1, k)) / dr
                          - (magn_fld(1, i, k) + magn_fld(1, i-1, k)) / (2. * dr * i)) * koef_h);
    }
}

vector3d<double> FieldE::get_field(double radius, double longitude)
//! function for electric field weighting
{
  vector3d<double> cmp(0., 0., 0.);

  int i_r = 0; // number of particle i cell
  int k_z = 0; // number of particle k cell
  int i_r_shift = 0; // number of particle i cell
  int k_z_shift = 0; // number of particle k cell

  double dr = geometry->r_cell_size;
  double dz = geometry->z_cell_size;
  double r1, r2, r3; // temp variables for calculation
  double dz1, dz2; // temp var.: width of k and k+1 cell
  double vol_1 = 0; // volume of i cell; Q/V, V - volume of elementary cell
  double vol_2 = 0; // volume of i+1 cell;

  r1 = radius - 0.5 * dr;
  r3 = radius + 0.5 * dr;

  // weighting of E_r
  // finding number of cell. example dr=0.5, radius = 0.7, i_r =0;!!
  i_r = CELL_NUMBER_OVERLAY(radius - 0.5 * dr, dr);
  k_z = CELL_NUMBER_OVERLAY(longitude, dz);
  i_r_shift = i_r - geometry->bottom_r_grid_number;
  k_z_shift = k_z - geometry->left_z_grid_number;
  // TODO: workaround: sometimes it gives -1.
  // Just get 0 cell if it happence
  if (i_r < 0) i_r = 0;
  if (k_z < 0) k_z = 0;
  if (i_r_shift < 0) i_r_shift = 0;
  if (k_z_shift < 0) k_z_shift = 0;

  // FIXME: it can be more, than current.size_x - 2
  // for some unknown reason
  if (i_r_shift > field[0].size_x() - 2)
  {
    MSG_FIXME("fieldE::get_field: i_r_shift is more, than field[0].size_x() - 2. Applying workaround");
    i_r_shift = field[0].size_x() - 2;
  }

  if (k_z_shift > field[0].size_y() - 2)
  {
    MSG_FIXME("fieldE::get_field: k_z_shift is more, than current[0].size_y() - 2 . Applying workaround");
      k_z_shift = field[0].size_y() - 2;
  }

  vol_1 = CELL_VOLUME(i_r+1, dr, dz);
  vol_2 = CELL_VOLUME(i_r+3, dr, dz);
  dz1 = (k_z+1) * dz - longitude;
  dz2 = longitude - k_z * dz;
  r2 = (i_r + 1) * dr;
  //weighting Er[i][k]//
  cmp[0] += field(0, i_r_shift, k_z_shift) * CYL_RNG_VOL(dz1, r1, r2) / vol_1;

  //weighting Er[i+1][k]//
  cmp[0] += field(0, i_r_shift + 1, k_z_shift) * CYL_RNG_VOL(dz1, r2, r3) / vol_2;

  //weighting Er[i][k+1]//
  cmp[0] += field(0, i_r_shift, k_z_shift + 1) * CYL_RNG_VOL(dz2, r1, r2) / vol_1;

  //weighting Er[i+1][k+1]//
  cmp[0] += field(0, i_r_shift + 1, k_z_shift + 1) * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

  // weighting of E_z
  // finding number of cell. example dz=0.5, longitude = 0.7, z_k =0;!!
  i_r = CELL_NUMBER_OVERLAY(radius, dr);
  k_z = CELL_NUMBER_OVERLAY(longitude - 0.5 * dz, dz);
  i_r_shift = i_r - geometry->bottom_r_grid_number;
  k_z_shift = k_z - geometry->left_z_grid_number;
  // TODO: workaround: sometimes it gives -1.
  // Just get 0 cell if it happence
  if (i_r < 0) i_r = 0;
  if (k_z < 0) k_z = 0;
  if (i_r_shift < 0) i_r_shift = 0;
  if (k_z_shift < 0) k_z_shift = 0;

  // FIXME: it can be more, than current.size_x - 2
  // for some unknown reason
  if (i_r_shift > field[0].size_x() - 2)
  {
    MSG_FIXME("fieldE::get_field: i_r_shift is more, than field[0].size_x() - 2. Applying workaround");
    i_r_shift = field[0].size_x() - 2;
  }

  if (k_z_shift > field[0].size_y() - 2)
  {
    MSG_FIXME("fieldE::get_field: k_z_shift is more, than current[0].size_y() - 2 . Applying workaround");
      k_z_shift = field[0].size_y() - 2;
  }

  if (radius > dr)
    vol_1 = CELL_VOLUME(i_r, dr, dz);
  else
    vol_1 = CYL_VOL(dz, dr); // volume of first cell

  r2 = (i_r + 0.5) * dr;
  vol_2 = CELL_VOLUME(i_r+2, dr, dz);
  dz1 = (k_z + 1.5) * dz - longitude;
  dz2 = longitude - (k_z + 0.5) * dz;

  // weighting Ez[i][k]
  cmp[2] += field(2, i_r_shift, k_z_shift) * CYL_RNG_VOL(dz1, r1, r2) / vol_1;

  // weighting Ez[i+1][k]
  cmp[2] += field(2, i_r_shift+1, k_z_shift) * CYL_RNG_VOL(dz1, r2, r3) / vol_2;

  // weighting Ez[i][k+1]
  cmp[2] += field(2, i_r_shift, k_z_shift+1) * CYL_RNG_VOL(dz2, r1, r2) / vol_1;

  //weighting Ez[i+1][k+1]//
  cmp[2] += field(2, i_r_shift+1, k_z_shift+1) * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

  // weighting of E_fi
  // finding number of cell. example dz=0.5, longitude = 0.7, z_k =1;

  i_r = CELL_NUMBER_OVERLAY(radius, dr);
  k_z = CELL_NUMBER_OVERLAY(longitude, dz);
  i_r_shift = i_r - geometry->bottom_r_grid_number;
  k_z_shift = k_z - geometry->left_z_grid_number;
  // TODO: workaround: sometimes it gives -1.
  // Just get 0 cell if it happence
  if (i_r < 0) i_r = 0;
  if (k_z < 0) k_z = 0;
  if (i_r_shift < 0) i_r_shift = 0;
  if (k_z_shift < 0) k_z_shift = 0;

  // FIXME: it can be more, than current.size_x - 2
  // for some unknown reason
  if (i_r_shift > field[0].size_x() - 2)
  {
    MSG_FIXME("fieldE::get_field: i_r_shift is more, than field[0].size_x() - 2. Applying workaround");
    i_r_shift = field[0].size_x() - 2;
  }

  if (k_z_shift > field[0].size_y() - 2)
  {
    MSG_FIXME("fieldE::get_field: k_z_shift is more, than current[0].size_y() - 2 . Applying workaround");
      k_z_shift = field[0].size_y() - 2;
  }

  if(radius>dr)
    vol_1 = CELL_VOLUME(i_r, dr, dz);
  else
    vol_1 = CYL_VOL(dz, dr); // volume of first cell

  r2 = (i_r + 0.5) * dr;
  vol_2 = CELL_VOLUME(i_r+2, dr, dz);
  dz1 = (k_z+1)*dz-longitude;
  dz2 = longitude - k_z * dz;

  // weighting Efi[i][k]
  cmp[1] += field(1, i_r_shift, k_z_shift) * CYL_RNG_VOL(dz1, r1, r2) / vol_1;

  // weighting Efi[i+1][k]
  cmp[1] += field(1, i_r_shift+1, k_z_shift) * CYL_RNG_VOL(dz1, r2, r3) / vol_2;

  // weighting Efi[i][k+1]
  cmp[1] += field(1, i_r_shift, k_z_shift+1) * CYL_RNG_VOL(dz2, r1, r2) / vol_1;

  // weighting Efi[i+1][k+1]
  cmp[2] += field(1, i_r_shift+1, k_z_shift+1) * CYL_RNG_VOL(dz2, r2, r3) / vol_2;

  return cmp;
}

#include "fieldH.hpp"

FieldH::FieldH(Geometry *geom, TimeSim *t, vector<SpecieP *> species) : Field(geom, t)
{
  // ! we use ``field'' member variable to mark magnetic field at half of timestep
  // ! and ``field_at_et'' to mark field at round timesteps
  // ! (same timestamp, as for electrical field)
  // ! it useful to have single implementation of ``get_field'' method
  // ! for both classes - FieldE and FieldH
  field_at_et = Grid3D<double> (geometry->r_grid_amount+1, geometry->z_grid_amount+1);
  species_p = species;

  field_at_et.setall(0);
}

void FieldH::calc_field_cylindrical()
{
  double dr = geometry->r_cell_size;
  double dz = geometry->z_cell_size;

// #pragma omp parallel for
  // calculate magnetic field r-component
  for(int i = 0; i < (geometry->r_grid_amount); i++)
    for(int k = 0; k < (geometry->z_grid_amount); k++)
    {
      double alpha_t = time->step
        * (field_at_et[1].get(i, k + 1) - field_at_et[1].get(i, k)) / (dz * MAGN_CONST);

      field[0].set(i, k, field_at_et[0].get(i, k) + alpha_t / 2.);
      field_at_et[0].set(i, k, field_at_et[0].get(i, k) + alpha_t);
    }

// calculate magnetic field phi-component
// #pragma omp parallel for
  for(int i = 0; i < (geometry->r_grid_amount); i++)
    for(int k = 0; k < (geometry->z_grid_amount); k++)
    {
      double alpha_t = time->step
        * ((field_at_et[2].get(i+1, k) - field_at_et[2].get(i, k)) / dr
           - (field_at_et[0].get(i, k+1) - field_at_et[0].get(i, k)) / dz)
        / MAGN_CONST;

      field[1].set(i, k, field_at_et[1].get(i, k) + alpha_t / 2.);
      field_at_et[1].set(i, k, field_at_et[1].get(i, k) + alpha_t);
    }

// calculate magnetic field z-component
// #pragma omp parallel for
  for(int i = 0; i < (geometry->r_grid_amount); i++)
    for(int k = 0; k < (geometry->z_grid_amount); k++)
    {
      double alpha_t = time->step
        * (
          (field_at_et[1].get(i+1, k) + field_at_et[1].get(i, k)) / (2. * dr * (i + 0.5))
          + (field_at_et[1].get(i+1, k) - field_at_et[1].get(i, k)) / dr)
        / MAGN_CONST;

      field[2].set(i, k, field_at_et[2].get(i, k) + alpha_t / 2.);
      field_at_et[2].set(i, k, field_at_et[2].get(i, k) + alpha_t);
    }
}

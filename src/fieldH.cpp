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

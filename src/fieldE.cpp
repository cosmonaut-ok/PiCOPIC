#include "fieldE.hpp"

FieldE::FieldE(Geometry *geom, TimeSim *t, vector<SpecieP *> species) : Field(geom, t)
{
  species_p = species;
  epsilon = Grid<double> (geometry->r_grid_amount+1, geometry->z_grid_amount+1);
  sigma = Grid<double> (geometry->r_grid_amount+1, geometry->z_grid_amount+1);

  epsilon.setall(EPSILON0);
  sigma.setall(0);

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
  Grid<double> j_r = current->current[0];
  Grid<double> j_phi = current->current[1];
  Grid<double> j_z = current->current[2];
  Grid<double> h_r = field_h->field_at_et[0];
  Grid<double> h_phi = field_h->field_at_et[1];
  Grid<double> h_z = field_h->field_at_et[2];
  double dr = geometry->r_cell_size;
  double dz = geometry->z_cell_size;

  // Er first[i] value
// #pragma omp parallel
  {
// #pragma omp for
    for(int k = 1; k < (geometry->z_grid_amount - 1); k++)
    {
      int i = 0;
      double epsilonx2 = 2 * epsilon.get(i, k);
      double sigma_t = sigma.get(i, k) * time->step;

      double koef_e = (epsilonx2 - sigma_t) / (epsilonx2 + sigma_t);

      double koef_h =  2 * time->step / (epsilonx2 + sigma_t);

      field[0].set(i, k, field[0].get(i, k) * koef_e
                   - (j_r.get(i, k)
                      + (h_phi.get(i, k) - h_phi.get(i, k-1)) / dz) * koef_h);
    }

    // Ez=on axis
// #pragma omp for
    for(int k = 0; k < (geometry->z_grid_amount - 1); k++)
    {
      int i = 0;
      double epsilonx2 = 2 * epsilon.get(i, k);
      double sigma_t = sigma.get(i, k) * time->step;

      double koef_e = (epsilonx2 - sigma_t) / (epsilonx2 + sigma_t);
      double koef_h = 2 * time->step / (epsilonx2 + sigma_t);

      field[2].set(i, k, field[2].get(i, k) * koef_e
                   - (j_z.get(i, k)
                      - 4. / dr * h_phi.get(i, k)) * koef_h);
    }

// #pragma omp for
    for(int i=1; i < (geometry->r_grid_amount - 1); i++)
      for(int k=1; k < (geometry->z_grid_amount - 1); k++)
      {
        double epsilonx2 = 2 * epsilon.get(i, k);
        double sigma_t = sigma.get(i, k) * time->step;

        double koef_e = (epsilonx2 - sigma_t) / (epsilonx2 + sigma_t);
        double koef_h = 2 * time->step / (epsilonx2 + sigma_t);

        field[0].set(i, k, field[0].get(i, k) * koef_e
                     - (j_r.get(i, k)
                        + (h_phi.get(i, k) - h_phi.get(i, k-1)) / dz) * koef_h);

        field[1].set(i, k, field[1].get(i, k) * koef_e
                     - (j_phi.get(i, k) - (h_r.get(i, k) - h_r.get(i, k-1))
                        / dz + (h_z.get(i, k) - h_z.get(i-1, k)) / dr) * koef_h);

        field[2].set(i, k, field[2].get(i, k) * koef_e
                     - (j_z.get(i, k) - (h_phi.get(i, k) - h_phi.get(i-1, k)) / dr
                        - (h_phi.get(i, k) + h_phi.get(i-1, k)) / (2. * dr * i))
                     * koef_h);
      }

// #pragma omp for
    for(int i=1; i < (geometry->r_grid_amount-1); i++)
    {
      int k = 0;
      double epsilonx2 = 2 * epsilon.get(i, k);
      double sigma_t = sigma.get(i, k) * time->step;

      double koef_e = (epsilonx2 - sigma_t) / (epsilonx2 + sigma_t);
      double koef_h = 2 * time->step / (epsilonx2 + sigma_t);

      field[2].set(i, k, field[2].get(i, k) * koef_e
                   - (j_z.get(i, k) - (h_phi.get(i, k) - h_phi.get(i-1, k)) / dr
                      - (h_phi.get(i, k) + h_phi.get(i-1, k)) / (2. * dr * i)
                     ) * koef_h);
    }
  }
}

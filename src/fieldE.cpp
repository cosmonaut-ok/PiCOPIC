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
void FieldE::set_pml() // double geometry->pml_length[1], double geometry->pml_length[3], double geometry->pml_length[2])
                       // double geometry->pml_sigma[0], double geometry->pml_sigma[1])
{
  // set_pml(pml_length[1], pml_length[3], pml_length[2],
  //         pml_sigma[0], pml_sigma[1]);

  // defining lenght of sigma calculation ion left wall
  double lenght_sigma_left = geometry->z_cell_size
    * (floor(geometry->z_grid_amount * geometry->pml_length[1]));

  // defining lenght of sigma calculation on right wall
  double lenght_sigma_right = geometry->z_cell_size
    * (floor(geometry->z_grid_amount * geometry->pml_length[3]));

  // defining lenght of sigma calculation on z-wall
  double lenght_sigma_extern = geometry->r_cell_size
    * (floor(geometry->r_grid_amount * geometry->pml_length[2]));

  // if pml is only on z walll
  if ((geometry->pml_length[1] == 0) && (geometry->pml_length[3] == 0)
      && (geometry->pml_length[2] != 0))
  {
    for(int i = 0; i < geometry->r_grid_amount; i++)
      for(int k = 0; k < geometry->z_grid_amount; k++)
        if (geometry->r_size - geometry->r_cell_size * i <= lenght_sigma_extern)
          sigma.set(i, k, geometry->pml_sigma[0] +
                    (geometry->pml_sigma[1] - geometry->pml_sigma[0])
                    / pow(lenght_sigma_extern, 2)
                    * pow((geometry->r_cell_size * (i + 1) - geometry->r_size
                           + lenght_sigma_extern), 2));
  }
  else
  {
    // sigma on r=0 and r=r walls
    for(int k = 0; k < geometry->z_grid_amount; k++)
      for(int i = 0; i < geometry->r_grid_amount; i++)
      {
        if (geometry->z_cell_size * k < lenght_sigma_left)
          sigma.set(i, k, geometry->pml_sigma[0]
                    + (geometry->pml_sigma[1] - geometry->pml_sigma[0])
                    / pow(lenght_sigma_left, 2)
                    * pow(lenght_sigma_left - geometry->z_cell_size * k, 2));
        if (geometry->z_size - geometry->z_cell_size * k < lenght_sigma_right)
          sigma.set(i, k, geometry->pml_sigma[0]
                    + (geometry->pml_sigma[1] - geometry->pml_sigma[0])
                    / pow(lenght_sigma_right, 2)
                    * pow(geometry->z_cell_size * (k + 1)
                          - geometry->z_size + lenght_sigma_right, 2));
      }
  }

  if (geometry->pml_length[2] != 0)
  {
    // sigma assigning
    // sigma on z=z wall
    for(int i = 0; i < geometry->r_grid_amount; i++)
      for(int k = 0; k < geometry->z_grid_amount; k++)
        if (geometry->r_size - geometry->r_cell_size * i <= lenght_sigma_extern)
          if (
            (geometry->r_size - geometry->r_cell_size * (i + 1)
             < lenght_sigma_extern * geometry->z_cell_size * k
             / lenght_sigma_left)
            && (geometry->r_size - geometry->r_cell_size * i
                <= lenght_sigma_extern / lenght_sigma_right
                * (geometry->z_size - geometry->z_cell_size * k))
            )
            sigma.set(i, k, geometry->pml_sigma[0]
                      + (geometry->pml_sigma[1] - geometry->pml_sigma[0])
                      / pow(lenght_sigma_extern, 2)
                      * pow(geometry->r_cell_size * (i + 1) - geometry->r_size
                            + lenght_sigma_extern, 2));
  }
}

void FieldE::calc_field_cylindrical()
{
  Grid3D<double> j = current->current;
  Grid3D<double> h = field_h->field_at_et;

  double dr = geometry->r_cell_size;
  double dz = geometry->z_cell_size;

  // Er first[i] value
  for(int k = 1; k < (geometry->z_grid_amount - 1); k++)
  {
    unsigned int i = 0;
    double epsilonx2 = 2 * epsilon.get(i, k) * EPSILON0;
    double sigma_t = sigma.get(i, k) * time->step;

    double koef_e = (epsilonx2 - sigma_t) / (epsilonx2 + sigma_t);

    double koef_h =  2 * time->step / (epsilonx2 + sigma_t);

    field[0].set(i, k, field[0].get(i, k) * koef_e
                 - (j[0].get(i, k) + (h[1].get(i, k) - h[1].get(i, k-1)) / dz)
                 * koef_h);
  }

  // Ez=on axis
  for(int k = 0; k < (geometry->z_grid_amount - 1); k++)
  {
    int i = 0;
    double epsilonx2 = 2. * epsilon.get(i, k) * EPSILON0;
    double sigma_t = sigma.get(i, k) * time->step;

    double koef_e = (epsilonx2 - sigma_t) / (epsilonx2 + sigma_t);
    double koef_h = 2 * time->step / (epsilonx2 + sigma_t);

    field[2].set(i, k, field[2].get(i, k) * koef_e
                 - (j[2].get(i, k) - 4. / dr * h[1].get(i, k))
                 * koef_h);
  }

  for(int i = 1; i < (geometry->r_grid_amount - 1); i++)
    for(int k = 1; k < (geometry->z_grid_amount - 1); k++)
    {
      double epsilonx2 = 2 * epsilon.get(i, k) * EPSILON0;
      double sigma_t = sigma.get(i, k) * time->step;

      double koef_e = (epsilonx2 - sigma_t) / (epsilonx2 + sigma_t);
      double koef_h = 2 * time->step / (epsilonx2 + sigma_t);

      field[0].set(i, k, field[0].get(i, k) * koef_e
                   - (j[0].get(i, k) + (h[1].get(i, k) - h[1].get(i, k-1)) / dz) * koef_h);

      field[1].set(i, k, field[1].get(i, k) * koef_e
                   - (j[1].get(i, k) - (h[0].get(i, k) - h[0].get(i, k-1))
                      / dz + (h[2].get(i, k) - h[2].get(i-1, k)) / dr) * koef_h);

      // if (k == 5 && i == 5)
      //   MSG("FFFFF " << field[1].get(i, k));

      field[2].set(i, k, field[2].get(i, k) * koef_e
                   - (j[2].get(i, k) - (h[1].get(i, k) - h[1].get(i-1, k)) / dr
                      - (h[1].get(i, k) + h[1].get(i-1, k)) / (2.0 * dr * i)) * koef_h);
    }

  for(int i=1; i < (geometry->r_grid_amount - 1); i++)
  {
    int k = 0;
    double epsilonx2 = 2 * epsilon.get(i, k) * EPSILON0;
    double sigma_t = sigma.get(i, k) * time->step;

    double koef_e = (epsilonx2 - sigma_t) / (epsilonx2 + sigma_t);
    double koef_h = 2 * time->step / (epsilonx2 + sigma_t);

    field[2].set(i, k, field[2].get(i, k) * koef_e -
                 (j[2].get(i, k) - (h[1].get(i, k) - h[1].get(i-1, k)) / dr
                  - (h[1].get(i, k) + h[1].get(i-1, k)) / (2. * dr * i)
                   ) * koef_h);
  }
}

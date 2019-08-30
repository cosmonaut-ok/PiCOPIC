#include "charge.hpp"

Charge::Charge(Geometry *geom, vector<SpecieP *> species): geometry(geom)
{
  // ! FIXME: that +1 is really temporary solution
  // ! to avoid Segfaults while weighting
  density = Grid<double>(geometry->r_grid_amount+1, geometry->z_grid_amount+1);
  species_p = species;
}

void Charge::weight_cylindrical ()
{
  double dr = geometry->r_cell_size;
  double dz = geometry->z_cell_size;

  for (auto sp = species_p.begin(); sp != species_p.end(); ++sp)
  {
    vector<vector <double> * > psp = (**sp).particles;

    // MSG(psp.size());
    for (auto pp = psp.begin(); pp != psp.end(); ++pp)
    {
      vector<double> *p = *pp; // single particle

      int r_i = CELL_NUMBER(P_POS_R((*p)), dr);
      int z_k = CELL_NUMBER(P_POS_Z((*p)), dz);

      // TODO: it could be negative for unknown reason
      if (r_i < 0) r_i = 0;
      if (z_k < 0) z_k = 0;

      int r_i_shift = r_i - geometry->bottom_r_grid_number;
      int z_k_shift = z_k - geometry->left_z_grid_number;

      // TODO: it could be negative for unknown reason
      if (r_i_shift < 0) r_i_shift = 0;
      if (z_k_shift < 0) z_k_shift = 0;

      double r_int, r_ext, r_cell_med; // radius of internal and external circles and cell medium point
      double v_cell_int, v_cell_ext; // volume of internal and external cell, related to macroparticle
      double dz_left, dz_right; // temp var.: width of k and k + 1 cell
      double rho_p, rho_p_weighted; // charge density for particle, weighted charge density for particle

      // MSG(geometry->bottom_r_grid_number << " " << geometry->left_z_grid_number);
      // in first cell other alg. of ro_v calc
      if (P_POS_R((*p)) > dr)
      {
        r_int =  P_POS_R((*p)) - 0.5 * dr; // r1
        r_cell_med = (r_i + 0.5) * dr; // r2
        r_ext = P_POS_R((*p)) + 0.5 * dr; // r3
        // charge density in cell-sized macroparticle
        rho_p = P_CHARGE((*p)) / (2. * PI * dz * dr * P_POS_R((*p)));
        // volumes of internal and external cell, regarding macroparticle
        v_cell_int = CELL_VOLUME(r_i, dr, dz); // v1
        v_cell_ext = CELL_VOLUME(r_i + 1, dr, dz); // v2

        dz_left = (z_k + 0.5) * dz - (P_POS_Z((*p)) - 0.5 * dz); // dz1
        dz_right = (P_POS_Z((*p)) + 0.5 * dz) - (z_k + 0.5) * dz; // dz2

        // weighting in rho[i][k] cell
        rho_p_weighted = rho_p * CYL_RNG_VOL(dz_left, r_int, r_cell_med) / v_cell_int;
        density.inc(r_i_shift, z_k_shift, rho_p_weighted);

        // weighting in rho[i+1][k] cell
        rho_p_weighted = rho_p * CYL_RNG_VOL(dz_left, r_cell_med, r_ext) / v_cell_ext;
        density.inc(r_i_shift+1, z_k_shift, rho_p_weighted);

        // weighting in rho[i][k+1] cell
        rho_p_weighted = rho_p * CYL_RNG_VOL(dz_right, r_int, r_cell_med) / v_cell_int;
        density.inc(r_i_shift, z_k_shift+1, rho_p_weighted);

        // weighting in rho[i+1][k+1] cell
        rho_p_weighted = rho_p * CYL_RNG_VOL(dz_right, r_cell_med, r_ext) / v_cell_ext;
        density.inc(r_i_shift+1, z_k_shift+1, rho_p_weighted);
      }
      else if (P_POS_R((*p)) <= dr / 2.)
      {
        //       r_i = 0;
        //       r1 =  0.;
        //       r2 = (r_i + 0.5) * dr;
        //       r3 = pos[i][0] + 0.5 * dr;
        //       dz1 = (z_k + 0.5) * dz - (pos[i][2] - 0.5 * dz);
        //       dz2 = (pos[i][2] + 0.5 * dz) - (z_k + 0.5) * dz;
        //       ro_v = charge_array[i] / (PI * dz * (2. * pos[i][0] * pos[i][0] + dr * dr / 2.));
        //       v_1 = CYL_VOL(dz, dr);
        //       v_2 = CELL_VOLUME(r_i + 1, dr, dz);
        //       ////////////////////////// /

        //       // weighting in ro[i][k] cell
        //       value = ro_v * PI * dz1 * (dr * dr / 2.-pos[i][0] * dr + pos[i][0] * pos[i][0]) / v_1;
        //       ro1->inc_rho(r_i, z_k, value);

        //       // weighting in ro[i + 1][k] cell
        //       value = ro_v * CYL_RNG_VOL(dz1, r2, r3) / v_2;
        //       ro1->inc_rho(r_i + 1,z_k, value);

        //       // weighting in ro[i][k + 1] cell
        //       value = ro_v*PI * dz2 * (dr * dr / 2.-pos[i][0] * dr + pos[i][0] * pos[i][0]) / v_1;
        //       ro1->inc_rho(r_i, z_k + 1, value);

        //       // weighting in ro[i + 1][k + 1] cell
        //       value = ro_v * CYL_RNG_VOL(dz2, r2, r3) / v_2;
        //       ro1->inc_rho(r_i + 1, z_k + 1, value);

      }
      else
      {
        //       ////////////////////////// /
        //       r1 = pos[i][0] - 0.5 * dr;
        //       r2 = (r_i + 0.5) * dr;
        //       r3 = pos[i][0] + 0.5 * dr;
        //       dz1 = (z_k + 0.5) * dz - (pos[i][2] - 0.5 * dz);
        //       dz2 = (pos[i][2] + 0.5 * dz) - (z_k + 0.5) * dz;
        //       ro_v = charge_array[i] / (2. * PI * dz * dr * pos[i][0]);
        //       v_1 = CYL_VOL(dz, dr);
        //       v_2 = CELL_VOLUME(r_i + 1, dr, dz);
        //       ////////////////////////// /

        //       // weighting in ro[i][k] cell
        //       value = ro_v * CYL_RNG_VOL(dz1, r1, r2) / v_1;
        //       ro1->inc_rho(r_i, z_k, value);

        //       // weighting in ro[i + 1][k] cell
        //       value = ro_v * CYL_RNG_VOL(dz1, r2, r3) / v_2;
        //       ro1->inc_rho(r_i + 1,z_k, value);

        //       // weighting in ro[i][k + 1] cell
        //       value = ro_v * CYL_RNG_VOL(dz2, r1, r2) / v_1;
        //       ro1->inc_rho(r_i, z_k + 1, value);

        //       // weighting in ro[i + 1][k + 1] cell
        //       value = ro_v * CYL_RNG_VOL(dz2, r2, r3) / v_2;
        //       ro1->inc_rho(r_i + 1, z_k + 1, value);
        //     }
      }
    }
  }
}

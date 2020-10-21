/*
 * This file is part of the PiCoPiC distribution (https://github.com/cosmonaut-ok/PiCoPiC).
 * Copyright (c) 2020 Alexander Vynnyk.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "densityCounted.hpp"

void DensityCounted::calc_density_cylindrical(string specie)
{
  for (auto ps = species_p.begin(); ps != species_p.end(); ++ps)
    if (specie.compare((**ps).name) == 0)
      for (auto i = (**ps).particles.begin(); i != (**ps).particles.end(); ++i)
        weight_cylindrical<double>(geometry, &density,
                                   P_POS_R((**i)),
                                   P_POS_Z((**i)),
                                   P_WEIGHT((**i))); // amount of particle in macroparticle

#if defined DENSITY_POSTPROC_BILINEAR
  Grid<double> density_src = density;
  lib::bilinear_interpolation<Grid<double>>(density_src, density);
#endif
}

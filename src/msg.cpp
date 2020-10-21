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

#include "msg.hpp"

using namespace std;

namespace msg
{
  int printed_step = 0;

  void print_probe_data(int probe_type, char* component, int step_number, int dump_interval, int* shape, double current_time)
  {
    int print_header_step = 20;
    if (printed_step % print_header_step == 0)
    {

      cout << endl
           << left << setw(8) << "Step"
           << left << setw(13) << "Saved Frame"
           << left << setw(20) << "Dumping Probe Name"
           << left << setw(21) << "Shape"
           << left << setw(18) << "Model Time (sec)"
           << left << setw(21) << "Simulation Duration"
           << endl;
    }

    char type_comp[100];
    char probe_shape[100];

    switch (probe_type)
    {
    case 0:
      sprintf(type_comp, "frame/%s", component);
      sprintf(probe_shape, "[%i,%i,%i,%i]", shape[0], shape[1], shape[2], shape[3]);
      break;
    case 1:
      sprintf(type_comp, "column/%s", component);
      sprintf(probe_shape, "[%i]", shape[1]);
      break;
    case 2:
      sprintf(type_comp, "row/%s", component);
      sprintf(probe_shape, "[%i]", shape[0]);
      break;
    case 3:
      sprintf(type_comp, "dot/%s", component);
      sprintf(probe_shape, "[%i,%i]", shape[0], shape[1]);
      break;
    case 4:
      sprintf(type_comp, "mpframe/%s", component);
      sprintf(probe_shape, "[%i,%i,%i,%i]", shape[0], shape[1], shape[2], shape[3]);
      break;
    }

    cout << left << setw(8) << step_number * dump_interval
         << left << setw(13) << step_number
         << left << setw(20) << type_comp
         << left << setw(21) << probe_shape
         << left << setw(18) << current_time
         << left << setw(18) << lib::get_simulation_duration()
         << endl;
    ++printed_step;
  }
}

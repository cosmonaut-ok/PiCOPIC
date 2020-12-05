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

using namespace constant;

namespace msg
{
  void print_header (int print_header_step, TimeSim *time)
  {
    //// print header
#ifndef ENABLE_DEBUG
    // int print_header_step = 30;
    if (time->print_header_counter % print_header_step == 0)
    {
      MSG("+------------+-------------+-------------------------------------------------+-------+------------------+---------------------+");
      MSG(std::left << std::setw(13) << "|    Step"
          << std::left << std::setw(14) << "| Saved Frame"
          << std::left << std::setw(50) << "|                Dumping Probe Name"
          << std::left << std::setw(8) << "| Shape"
          << std::left << std::setw(19) << "| Model Time (sec)"
          << std::left << std::setw(22) << "| Simulation Duration"
          << std::left << "|"
        );
      MSG("+------------+-------------+-------------------------------------------------+-------+------------------+---------------------+");
    }

    ++time->print_header_counter;

#endif // ENABLE_DEBUG
  }

  void print_values (std::string probe_name, std::string shape_name,
                     unsigned int probe_schedule, TimeSim *time)
  {
#ifndef ENABLE_DEBUG
    int current_time_step = ceil(time->current / time->step);
    std::string dump_step = std::to_string((int)current_time_step / probe_schedule);
    std::ostringstream time_conv;
    time_conv << time->current;
    std::string cur_time_str = time_conv.str();
    std::string sim_duration = algo::common::get_simulation_duration();

    MSG(std::left << std::setw(13) << "| " + std::to_string(current_time_step)
        << std::left << std::setw(14) << "| " + dump_step
        << std::left << std::setw(50) << "| " + probe_name
        << std::left << std::setw(8) << "| " + shape_name
        << std::left << std::setw(19) << "| " + cur_time_str
        << std::left << std::setw(22) << "| " + sim_duration
        << std::left << "|"
      );
#endif
  }

  void print_final ()
  {
#ifndef DEBUG
    MSG("+------------+-------------+-------------------------------------------------+-------+------------------+---------------------+");
#endif
    LOG_S(INFO) << "SIMULATION COMPLETE";
  }
}

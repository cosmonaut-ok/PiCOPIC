/*
 * This file is part of the PiCoPiC distribution
 * (https://github.com/cosmonaut-ok/PiCoPiC).
 * Copyright (c) 2020,2025 Alexander Vynnyk.
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

#ifndef _MSG_HPP_
#define _MSG_HPP_

#include "spdlog/spdlog.h"

#include <iomanip> // std::setw
#include <iostream>

#include "algo/common.hpp"
#include "defines.hpp"
#include "timeSim.hpp"

#define MSG(message) std::cout << message << std::endl;

#ifndef ENABLE_FIXME
#define ENABLE_FIXME false
#endif

#define MSG_FIXME(message)                                                                                             \
  if ( ENABLE_FIXME )                                                                                                  \
    std::cerr << "FIXME: " << message << "." << std::endl;

namespace msg
{
void print_header (int print_header_step, TimeSim* time);
void print_values (std::string probe_name, std::string shape_name, unsigned int probe_schedule, TimeSim* time);
void print_final ();
} // namespace msg
#endif // end of _MSG_HPP_

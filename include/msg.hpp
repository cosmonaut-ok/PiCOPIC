/* 
 * This file is part of the PiCOPIC distribution (https://github.com/cosmonaut-ok/PiCOPIC).
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

#ifndef _MSG_HPP_
#define _MSG_HPP_

#include <iostream>
#include <iomanip> // std::setw

#include "lib.hpp"

#define MSG(message) std::cout << message << std::endl;

#ifndef FIXME
#define FIXME false
#endif

#define MSG_FIXME(message) if (FIXME) std::cerr << "FIXME: " << message << "." << std::endl;

#endif // end of _MSG_HPP_

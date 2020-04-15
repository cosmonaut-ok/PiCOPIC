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

// #ifndef _MSG_HPP_
// #define _MSG_HPP_

#include <iostream>
#include <iomanip> // std::setw

#include "lib.hpp"

#ifndef DEBUG
#define DEBUG false
#endif

#ifndef FIXME
#define FIXME false
#endif

#define MSG(message) std::cout << message << std::endl;

#define LOG_DBG(message) if (DEBUG) std::cerr << "DEBUG: " << message << "." << std::endl;

#define MSG_FIXME(message) if (FIXME) std::cerr << "FIXME: " << message << "." << std::endl;

#define LOG_INFO(message) std::cerr << "INFO: " << message << "." << std::endl;

#define LOG_WARN(message) std::cerr << "WARNING: " << message << "." << std::endl;

#define LOG_ERR(message) std::cerr << "ERROR: " << message << "." << std::endl;

#define LOG_CRIT(message, exitcode) { std::cerr << "ERROR: " << message << ". Can not continue. Exiting." << std::endl; exit(exitcode); }

// #endif // end of _MSG_HPP_

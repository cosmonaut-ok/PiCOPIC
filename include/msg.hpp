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

#include "loguru.hpp"

#ifdef LOGURU_WITH_STREAMS
#undef LOGURU_WITH_STREAMS
#define LOGURU_WITH_STREAMS 1
#endif

#ifdef LOGURU_SCOPE_TIME_PRECISION
#undef LOGURU_SCOPE_TIME_PRECISION
#define LOGURU_SCOPE_TIME_PRECISION 3
#endif

#ifdef LOGURU_REPLACE_GLOG
#undef LOGURU_REPLACE_GLOG
#define LOGURU_REPLACE_GLOG 1
#endif

#ifdef LOGURU_USE_FMTLIB
#undef LOGURU_USE_FMTLIB
#define LOGURU_USE_FMTLIB 1
#endif

#ifndef _MSG_HPP_
#define _MSG_HPP_

#include <iostream>
#include <iomanip> // std::setw

#include "algo/common.hpp"
#include "defines.hpp"

#define MSG(message) std::cout << message << std::endl;

#ifndef ENABLE_FIXME
#define ENABLE_FIXME false
#endif

#define MSG_FIXME(message) if (ENABLE_FIXME) std::cerr << "FIXME: " << message << "." << std::endl;

#endif // end of _MSG_HPP_

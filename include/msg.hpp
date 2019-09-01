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

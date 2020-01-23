#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)

AC_INIT(PicoPIC, esyscmd(git log --pretty=%h -1 | tr -d \\n), cosmonaut.ok@zoho.com)

AC_CONFIG_HEADERS([include/defines.hpp])

AC_OUTPUT

m4_pattern_allow([^AC_])

dnl AC_CONFIG_AUX_DIR([acdir])
dnl AC_CONFIG_MACRO_DIR([acdir])

# AC_CANONICAL_SYSTEM

AC_SUBST(EXPERIMENTAL_OPTION)
AC_SUBST(LEGACY_OPTION)
AC_SUBST(IEEE)
AC_SUBST(DEBUG)
AC_SUBST(SINGLETHREAD)
AC_SUBST(HDF5)
AC_SUBST(CFLAGS_DEFAULT)
AC_SUBST(CFLAGS_ADDITIONAL)
AC_SUBST(LDFLAGS_ADDITIONAL)
AC_SUBST(DOXYGEN_FORMATS)
AC_SUBST(HDF5_OPTION)
AC_SUBST(PUSHER)
AC_SUBST(OMP_DYNAMIC)
AC_SUBST(DEBUG_OPTION)
AC_SUBST(PROFILER)
AC_SUBST(HDF5_URL)
AC_SUBST(HDF5_VERSION)

AC_SUBST(MKDIR)
MKDIR=mkdir

AC_SUBST(RELEASE)
RELEASE=${PACKAGE_VERSION}

HDF5_URL="https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.5/src/hdf5-1.10.5.tar.bz2"
HDF5_VERSION="hdf5-1.10.5"

dnl AC_ARG_ENABLE([experimental], [AC_HELP_STRING([--enable-experimental], [Enable experimental features (only for debug/development)])], [ENABLE_EXPERIMENTAL="$enableval"], [ENABLE_EXPERIMENTAL=no])
dnl AC_ARG_ENABLE([legacy], [AC_HELP_STRING([--enable-legacy], [Use legacy instead of current features])], [ENABLE_LEGACY="$enableval"], [ENABLE_LEGACY=no])
AC_ARG_ENABLE([debug], [AC_HELP_STRING([--enable-debug], [Enable debug (optimized for GDB)])], [DEBUG="$enableval"], [DEBUG=no])
AC_ARG_ENABLE([fixme], [AC_HELP_STRING([--enable-fixme], [Enable FIXME output ])], [FIXME="$enableval"], [FIXME=${DEBUG}])
AC_ARG_ENABLE([profiler], [AC_HELP_STRING([--enable-profiler], [Enable profiler support for performance debug])], [PROFILER="$enableval"], [PROFILER=no])
AC_ARG_ENABLE([ieee], [AC_HELP_STRING([--enable-ieee], [Keep IEEE compliant rounding and disable hardware acceleration. Decrease calculation speed])], [IEEE="$enableval"], [IEEE=${DEBUG}])
AC_ARG_ENABLE([singlethread], [AC_HELP_STRING([--enable-singlethread], [Build application to run in single-thread mode (used in debug by default)])], [SINGLETHREAD="$enableval"], [SINGLETHREAD=${DEBUG}])
AC_ARG_ENABLE([omp_dynamic], [AC_HELP_STRING([--enable-omp-dynamic], [Enable openmp dynamic thread creation])], [OMP_DYNAMIC="$enableval"], [OMP_DYNAMIC=no])

AC_ARG_ENABLE([hdf5], [AC_HELP_STRING([--disable-hdf5], [Build application with hdf5 library support])], [HDF5="$enableval"], [HDF5=yes])

dnl AC_ARG_WITH([doc-formats], [AC_HELP_STRING([--with-doc-formats], [Set list of documentation formats to be generated])], [DOXYGEN_FORMATS="$withval"], [DOXYGEN_FORMATS="latex html rtf"])
AC_ARG_WITH([pusher], [AC_HELP_STRING([--with-pusher], [Set particles pusher (supported: boris (default), vay, hc (n/a), rrll (n/a))])], [PUSHER="$withval"], [PUSHER="boris-adaptive"])

AC_ARG_WITH([plasma-spatial-distribution], [AC_HELP_STRING([--with-plasma-spatial-distribution], [Set plasma particles spatial distribution type (supported: flat (default), random, regular, centered (not implemented) )])], [PLASMA_SPATIAL="$withval"], [PLASMA_SPATIAL="flat"])

AC_ARG_WITH([plasma-velocity-distribution], [AC_HELP_STRING([--with-plasma-velocity-distribution], [Set plasma particles velocity distribution (supported: thermal (default), rectangular, eigen (singular velocities))])], [PLASMA_VELOCITY="$withval"], [PLASMA_VELOCITY="thermal"])


AC_PROG_CXXCPP
AC_LANG(C++)

# set default and customized options
CFLAGS_DEFAULT="-m64 -mcmodel=medium -std=c++17"

# singlethread
if test x$SINGLETHREAD == xyes; then
     CFLAGS_ADDITIONAL+=" -Wno-unknown-pragmas"
  AC_DEFINE([SINGLETHREAD], [true], [Debug option])
else
  AC_OPENMP
  CFLAGS_ADDITIONAL+=" -fopenmp"
fi

if test x$OMP_DYNAMIC == xyes; then
   AC_DEFINE([OPENMP_DYNAMIC_THREADS], [true], [Set OpenMP threads to be created dynamically])
fi

# speedup
if test x$IEEE != xyes; then
  CFLAGS_ADDITIONAL="$CFLAGS_ADDITIONAL -ffast-math -O3 -march=native -funroll-loops -Wno-write-strings -freciprocal-math -fcx-fortran-rules -msse -mavx"
  AC_DEFINE([SPEEDUP], [true], [Directive to increase calculation speed])
fi

dnl debug
if test x$DEBUG == xyes; then
  CFLAGS_ADDITIONAL+=" -O0 -Wall -Wextra -ggdb3 -fvar-tracking -ggnu-pubnames -pedantic -time -ftree-vectorizer-verbose=7"
  DEBUG_OPTION=true
  AC_DEFINE([DEBUG], [true], [Use debug option])
else
  DEBUG_OPTION=false
fi

dnl debug
if test x$FIXME == xyes; then
  AC_DEFINE([FIXME], [true], [Use FIXME output])
fi

# profiler
if test x$PROFILER = xyes; then
  CFLAGS_ADDITIONAL="$CFLAGS_ADDITIONAL -pg -no-pie"
  LDFLAGS_ADDITIONAL="$LDFLAGS_ADDITIONAL -pg -no-pie"
  AC_DEFINE([PROFILER], [true], [Performance debug option])
fi

if test x$HDF5 = xyes; then
   AC_HAVE_LIBRARY(hdf5, , AC_MSG_WARN([No such library hdf5. You can build it locally after configuration with `make hdf5-install` before build and it includes to build paths automatically]))
   HDF5_OPTION=true
   AC_DEFINE([USE_HDF5], [true], [Use HDF5])
else
   HDF5_OPTION=false
fi

# define a pusher
if test x$PUSHER == xboris-adaptive; then
   AC_DEFINE_UNQUOTED([PUSHER_BORIS_ADAPTIVE], [true], [Performance debug option])
   AC_DEFINE([REL_LIMIT], [5e7], [Relativistit calculations limit (for boris-adaptive pusher)])
elif test x$PUSHER == xboris; then
   AC_DEFINE_UNQUOTED([PUSHER_BORIS_CLASSIC], [true], [Performance debug option])
elif test x$PUSHER == xboris-relativistic; then
   AC_DEFINE_UNQUOTED([PUSHER_BORIS_RELATIVISTIC], [true], [Performance debug option])
elif test x$PUSHER == xvay; then
   AC_DEFINE_UNQUOTED([PUSHER_VAY], [true], [Performance debug option])
else
   AC_MSG_ERROR([Pusher is not supported.])
fi

# define plasma spatial distribution type
if test x$PLASMA_SPATIAL == xrandom; then
   AC_DEFINE_UNQUOTED([PLASMA_SPATIAL_RANDOM], [true], [Set plasma spatial distribution type to random])
elif test x$PLASMA_SPATIAL == xregular; then
   AC_DEFINE_UNQUOTED([PLASMA_SPATIAL_REGULAR], [true], [Set plasma spatial distribution type to regular])
elif test x$PLASMA_SPATIAL == xcentered; then
   AC_DEFINE_UNQUOTED([PLASMA_SPATIAL_CENTERED], [true], [Set plasma spatial distribution type to cell-centered])
elif test x$PLASMA_SPATIAL == xflat; then
   AC_DEFINE_UNQUOTED([PLASMA_SPATIAL_FLAT], [true], [Set plasma spatial distribution type to "flat"])
else
   AC_MSG_ERROR([Plasma spatial distribution type is not supported.])
fi

# define plasma spatial distribution type
if test x$PLASMA_VELOCITY == xthermal; then
   AC_DEFINE_UNQUOTED([PLASMA_VELOCITY_THERMAL], [true], [Set plasma spatial distribution type to random])
elif x$PLASMA_VELOCITY == xrectangular; then
   AC_DEFINE_UNQUOTED([PLASMA_VELOCITY_RECTANGULAR], [true], [Set plasma spatial distribution type to regular ("flat" distribution)])
elif x$PLASMA_VELOCITY == xeigen; then
   AC_DEFINE_UNQUOTED([PLASMA_VELOCITY_EIGEN], [true], [Set plasma spatial distribution type to cell-centered])
else
   AC_MSG_ERROR([Plasma velocity distribution type is not supported.])
fi

AC_OUTPUT(Makefile)

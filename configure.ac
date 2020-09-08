#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)

AC_INIT(PiCOPIC, esyscmd(git log --pretty=%h -1 | tr -d \\n), alexander.vynnyk@protonmail.com)

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
AC_SUBST(CFLAGS_DEFAULT)
AC_SUBST(CFLAGS_ADDITIONAL)
AC_SUBST(LDFLAGS_ADDITIONAL)
AC_SUBST(DOXYGEN_FORMATS)
AC_SUBST(HDF5_OPTION)
AC_SUBST(PUSHER)
AC_SUBST(OMP_DYNAMIC)
AC_SUBST(DEBUG_OPTION)
AC_SUBST(PROFILER)
AC_SUBST(HDF5)
AC_SUBST(HDF5_URL)
AC_SUBST(HDF5_VERSION)

AC_SUBST(MKDIR)
MKDIR=mkdir

AC_SUBST(RELEASE)
RELEASE=${PACKAGE_VERSION}

HDF5_URL="https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.0/src/hdf5-1.12.0.tar.bz2"
HDF5_VERSION="hdf5-1.12.0"

dnl AC_ARG_ENABLE([experimental], [AC_HELP_STRING([--enable-experimental], [Enable experimental features (only for debug/development)])], [ENABLE_EXPERIMENTAL="$enableval"], [ENABLE_EXPERIMENTAL=no])
dnl AC_ARG_ENABLE([legacy], [AC_HELP_STRING([--enable-legacy], [Use legacy instead of current features])], [ENABLE_LEGACY="$enableval"], [ENABLE_LEGACY=no])
AC_ARG_ENABLE([debug], [AC_HELP_STRING([--enable-debug], [Enable debug (optimized for GDB)])], [DEBUG="$enableval"], [DEBUG=no])
AC_ARG_ENABLE([fixme], [AC_HELP_STRING([--enable-fixme], [Enable FIXME output ])], [FIXME="$enableval"], [FIXME=${DEBUG}])
AC_ARG_ENABLE([profiler], [AC_HELP_STRING([--enable-profiler], [Enable profiler support for performance debug])], [PROFILER="$enableval"], [PROFILER=no])
AC_ARG_ENABLE([ieee], [AC_HELP_STRING([--enable-ieee], [Keep IEEE compliant rounding and disable hardware acceleration. Decrease calculation speed])], [IEEE="$enableval"], [IEEE=${DEBUG}])
AC_ARG_ENABLE([singlethread], [AC_HELP_STRING([--enable-singlethread], [Build application to run in single-thread mode (used in debug by default)])], [SINGLETHREAD="$enableval"], [SINGLETHREAD=${DEBUG}])
AC_ARG_ENABLE([omp_dynamic], [AC_HELP_STRING([--enable-omp-dynamic], [Enable openmp dynamic thread creation])], [OMP_DYNAMIC="$enableval"], [OMP_DYNAMIC=no])

AC_ARG_ENABLE([hdf5], [AC_HELP_STRING([--disable-hdf5], [Build application with hdf5 library support])], [HDF5="$enableval"], [HDF5=yes])

AC_ARG_WITH([doc-formats], [AC_HELP_STRING([--with-doc-formats], [Set list of documentation formats to be generated])], [DOXYGEN_FORMATS="$withval"], [DOXYGEN_FORMATS="latex html rtf"])
AC_ARG_WITH([pusher], [AC_HELP_STRING([--with-pusher], [Set particles pusher (supported: boris (default), vay, higuera-cary)])], [PUSHER="$withval"], [PUSHER="boris-adaptive"])
AC_ARG_WITH([temperature-calc], [AC_HELP_STRING([--with-temperature-calc], [Set temperature calculation algorithm. Supported: counting (default), weighting])], [TEMP_CALC="$withval"], [TEMP_CALC="counting"])
AC_ARG_WITH([density-calc], [AC_HELP_STRING([--with-density-calc], [Set density calculation algorithm. Supported: counting (default), weighting])], [DENSITY_CALC="$withval"], [DENSITY_CALC="counting"])

AC_ARG_WITH([plasma-spatial-distribution], [AC_HELP_STRING([--with-plasma-spatial-distribution], [Set plasma particles spatial distribution type (supported: flat (default), random, regular, centered)])], [PLASMA_SPATIAL="$withval"], [PLASMA_SPATIAL="flat"])

AC_ARG_WITH([plasma-velocity-distribution], [AC_HELP_STRING([--with-plasma-velocity-distribution], [Set plasma particles velocity distribution (supported: thermal (default), rectangular, eigen (singular velocities))])], [PLASMA_VELOCITY="$withval"], [PLASMA_VELOCITY="thermal"])

AC_ARG_WITH([current-solver], [AC_HELP_STRING([--with-current-solver], [Set charge concervation scheme of current solver (supported: villasenor-buneman (default) and zigzag)])], [CCS="$withval"], [CCS="villasenor-buneman"])

AC_ARG_ENABLE([coulomb-collisions], [AC_HELP_STRING([--enable-coulomb-collisions], [Enable coulomb collisions. (WARNING! This is an experimental unfinished feature. This switch is only for development purposes.)])], [COLLISIONS="$enableval"], [COLLISIONS=no])

AC_ARG_WITH([coulomb-collisions-scheme], [AC_HELP_STRING([--with-coulomb-collisions-scheme], [Set coulomb collisions scheme (supported: ta77s [10.1002/ctpp.201700121] (default), sk98m [10.1063/1.4742167] and sk98s [sk with symmetric weighting])])], [COULOMB_COLLISIONS_SCHEME="$withval"], [COULOMB_COLLISIONS_SCHEME="ta77s"])

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
  CFLAGS_ADDITIONAL+=" -fopenmp -fopenmp-simd"
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

AC_DEFINE([REL_LIMIT], [5e7], [Relativistit calculations limit (for boris-adaptive pusher)])

# define a pusher
if test x$PUSHER == xboris-adaptive; then
   AC_DEFINE_UNQUOTED([PUSHER_BORIS_ADAPTIVE], [true], [Performance debug option])
elif test x$PUSHER == xboris; then
   AC_DEFINE_UNQUOTED([PUSHER_BORIS_CLASSIC], [true], [Performance debug option])
elif test x$PUSHER == xboris-relativistic; then
   AC_DEFINE_UNQUOTED([PUSHER_BORIS_RELATIVISTIC], [true], [Performance debug option])
elif test x$PUSHER == xvay; then
   AC_DEFINE_UNQUOTED([PUSHER_VAY], [true], [Performance debug option])
elif test x$PUSHER == xhiguera-cary; then
   AC_DEFINE_UNQUOTED([PUSHER_HIGUERA_CARY], [true], [Performance debug option])
else
   AC_MSG_ERROR([Pusher is not supported.])
fi

# define a temperature claculation algorithm
if test x$TEMP_CALC == xweighting; then
   AC_DEFINE_UNQUOTED([TEMP_CALC_WEIGHTING], [true], [Performance debug option])
elif test x$TEMP_CALC == xcounting; then
   AC_DEFINE_UNQUOTED([TEMP_CALC_COUNTING], [true], [Performance debug option])
else
   AC_MSG_ERROR([Temperature calculation algorithm is not supported.])
fi

# define a temperature claculation algorithm
if test x$DENSITY_CALC == xweighting; then
   AC_DEFINE_UNQUOTED([DENSITY_CALC_WEIGHTING], [true], [Performance debug option])
elif test x$DENSITY_CALC == xcounting; then
   AC_DEFINE_UNQUOTED([DENSITY_CALC_COUNTING], [true], [Performance debug option])
else
   AC_MSG_ERROR([Density calculation algorithm is not supported.])
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
elif test x$PLASMA_VELOCITY == xrectangular; then
   AC_DEFINE_UNQUOTED([PLASMA_VELOCITY_RECTANGULAR], [true], [Set plasma spatial distribution type to regular ("flat" distribution)])
elif test x$PLASMA_VELOCITY == xeigen; then
   AC_DEFINE_UNQUOTED([PLASMA_VELOCITY_EIGEN], [true], [Set plasma spatial distribution type to cell-centered])
else
   AC_MSG_ERROR([Plasma velocity distribution type is not supported.])
fi

# define charge concervation scheme
if test x$CCS == xvillasenor-buneman; then
   AC_DEFINE_UNQUOTED([CCS_VILLASENOR_BUNEMAN], [true], [Set plasma charge conservation scheme to Villasenor-Buneman])
elif test x$CCS == xzigzag; then
   AC_DEFINE_UNQUOTED([CCS_ZIGZAG], [true], [Set plasma charge conservation scheme to ZigZag])
else
   AC_MSG_ERROR([Plasma charge conservation scheme is not supported.])
fi

if test x$COLLISIONS = xyes; then
   AC_DEFINE_UNQUOTED([COLLISIONS], [true], [Enable coulomb collisions])
   AC_MSG_WARN([Coulomb collisions is an unfinished experimental feature. Currently it is ONLY for development purposes. Ensure, you know what you do by enabling it.])
fi

# define charge concervation scheme
if test x$COULOMB_COLLISIONS_SCHEME == xta77s; then
   AC_DEFINE_UNQUOTED([COULOMB_COLLISIONS_TA77S], [true], [Set coulomb collision scheme to 10.1002/ctpp.201700121])
elif test x$COULOMB_COLLISIONS_SCHEME == xsk98m; then
   AC_DEFINE_UNQUOTED([COULOMB_COLLISIONS_SK98M], [true], [Set coulomb collision scheme to 10.1063/1.4742167])
elif test x$COULOMB_COLLISIONS_SCHEME == xsk98m; then
   AC_DEFINE_UNQUOTED([COULOMB_COLLISIONS_SK98S], [true], [Set coulomb collision scheme to 10.1063/1.4742167 with weighting from 10.1002/ctpp.201700121])
else
   AC_MSG_ERROR([Coulomb collisions scheme is not supported.])
fi

AC_DEFINE_UNQUOTED([CXXFLAGS], ["$CFLAGS_DEFAULT $CFLAGS_ADDITIONAL"], [define used C(XX)FLAGS])

# define loguru serve variables
AC_DEFINE_UNQUOTED([LOGURU_WITH_STREAMS], [1], [Enable loguru stream behavior])
AC_DEFINE_UNQUOTED([LOGURU_SCOPE_TIME_PRECISION], [3], [Resolution of loguru scope timers. 3=ms, 6=us, 9=ns])
AC_DEFINE_UNQUOTED([LOGURU_REPLACE_GLOG], [1], [Replace GLOG by Loguru])
AC_DEFINE_UNQUOTED([LOGURU_USE_FMTLIB], [0], [Use fmtlib])

AC_OUTPUT(Makefile)

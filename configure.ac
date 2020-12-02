#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)

AC_INIT(PiCoPiC, esyscmd(git log --pretty=%h -1 | tr -d \\n), alexander.vynnyk@protonmail.com)

AC_CONFIG_HEADERS([include/defines.hpp])

AC_OUTPUT

m4_pattern_allow([^AC_])

dnl AC_CONFIG_AUX_DIR([acdir])
dnl AC_CONFIG_MACRO_DIR([acdir])

# AC_CANONICAL_SYSTEM

AC_SUBST(ENABLE_EXPERIMENTAL)
AC_SUBST(ENABLE_DEVELOPMENT)
AC_SUBST(ENABLE_DEBUG)
AC_SUBST(ENABLE_FIXME)
AC_SUBST(ENABLE_PROFILER)
AC_SUBST(ENABLE_IEEE)
AC_SUBST(ENABLE_SINGLETHREAD)
AC_SUBST(ENABLE_OMP_DYNAMIC)
AC_SUBST(ENABLE_MPI)
AC_SUBST(ENABLE_HDF5)
AC_SUBST(HDF5_URL)
AC_SUBST(HDF5_VERSION)
AC_SUBST(WITH_DOXYGEN_FORMATS)
AC_SUBST(WITH_PUSHER)
AC_SUBST(WITH_CCS)

AC_SUBST(CFLAGS_DEFAULT)
AC_SUBST(CFLAGS_ADDITIONAL)
AC_SUBST(LDFLAGS_ADDITIONAL)

AC_SUBST(RELEASE)
RELEASE=${PACKAGE_VERSION}

HDF5_URL="https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.0/src/hdf5-1.12.0.tar.bz2"
HDF5_VERSION="hdf5-1.12.0"

AC_ARG_ENABLE([experimental],
              [AC_HELP_STRING([--enable-experimental],
              [Unlock experimental features support (features may work incorrectly)])],
              [EXPERIMENTAL_OPTION="$enableval"], [EXPERIMENTAL_OPTION=no])

AC_ARG_ENABLE([development],
              [AC_HELP_STRING([--enable-development],
              [Unlock development and/or unfinished features (for development purposes only)])],
              [DEVELOPMENT_OPTION="$enableval"], [DEVELOPMENT_OPTION=no])

AC_ARG_ENABLE([debug], [AC_HELP_STRING([--enable-debug],
              [Enable debug build (optimized for GDB)])],
              [DEBUG_OPTION="$enableval"], [DEBUG_OPTION=no])

AC_ARG_ENABLE([fixme], [AC_HELP_STRING([--enable-fixme],
              [Enable FIXME output ])],
              [FIXME_OPTION="$enableval"], [FIXME_OPTION=${DEBUG_OPTION}])

AC_ARG_ENABLE([profiler], [AC_HELP_STRING([--enable-profiler],
              [Enable profiler support for performance debug])],
              [PROFILER_OPTION="$enableval"], [PROFILER_OPTION=no])

AC_ARG_ENABLE([ieee], [AC_HELP_STRING([--enable-ieee],
              [Keep IEEE compliant rounding and disable hardware acceleration. Decrease calculation speed])],
              [IEEE_OPTION="$enableval"], [IEEE_OPTION=${DEBUG_OPTION}])

AC_ARG_ENABLE([singlethread], [AC_HELP_STRING([--enable-singlethread],
              [Build application to run in single-thread mode (used in debug build by default)])],
              [SINGLETHREAD_OPTION="$enableval"], [SINGLETHREAD_OPTION=${DEBUG_OPTION}])

AC_ARG_ENABLE([omp_dynamic], [AC_HELP_STRING([--enable-omp-dynamic],
              [Enable OpenMP dynamic thread creation])],
              [OMP_DYNAMIC_OPTION="$enableval"], [OMP_DYNAMIC_OPTION=no])

AC_ARG_ENABLE([mpi], [AC_HELP_STRING([--enable-mpi],
              [Enable openMPI support])],
              [MPI_OPTION="$enableval"], [MPI_OPTION=no])

dnl AC_ARG_ENABLE([hdf5], [AC_HELP_STRING([--disable-hdf5],
dnl               [Build application with hdf5 library support])],
dnl               [HDF5_OPTION="$enableval"], [HDF5_OPTION=yes])

AC_ARG_WITH([doc-formats], [AC_HELP_STRING([--with-doc-formats],
            [Set list of documentation formats to be generated])],
            [WITH_DOXYGEN_FORMATS="$withval"], [WITH_DOXYGEN_FORMATS="latex html rtf"])

AC_ARG_WITH([pusher], [AC_HELP_STRING([--with-pusher],
            [Set particles pusher (supported: boris (default), boris-adaptive, boris-relativistic, vay, hc [higuera-cary] )])],
            [WITH_PUSHER="$withval"], [WITH_PUSHER="boris-adaptive"])

AC_ARG_WITH([temperature-calc], [AC_HELP_STRING([--with-temperature-calc],
            [Set temperature calculation algorithm. Supported: counting (default), weighting])],
            [WITH_TEMP_CALC="$withval"], [WITH_TEMP_CALC="counting"])

AC_ARG_WITH([density-calc], [AC_HELP_STRING([--with-density-calc],
            [Set density calculation algorithm. Supported: counting (default), weighting])],
            [WITH_DENSITY_CALC="$withval"], [WITH_DENSITY_CALC="counting"])

AC_ARG_WITH([plasma-spatial-distribution], [AC_HELP_STRING([--with-plasma-spatial-distribution],
            [Set plasma particles spatial distribution type (supported: flat (default), random, regular, centered)])],
            [WITH_PLASMA_SPATIAL="$withval"], [WITH_PLASMA_SPATIAL="flat"])

AC_ARG_WITH([plasma-velocity-distribution], [AC_HELP_STRING([--with-plasma-velocity-distribution],
            [Set plasma particles velocity distribution (supported: thermal (default), rectangular, eigen (singular velocities))])],
            [WITH_PLASMA_VELOCITY="$withval"], [WITH_PLASMA_VELOCITY="thermal"])

AC_ARG_WITH([current-solver], [AC_HELP_STRING([--with-current-solver],
            [Set charge conservation scheme of current solver (supported: vb [villasenor-buneman] (default) and zigzag)])],
            [WITH_CCS="$withval"], [WITH_CCS="vb"])

AC_ARG_ENABLE([coulomb-collisions], [AC_HELP_STRING([--enable-coulomb-collisions],
              [Enable coulomb collisions. (WARNING! This is an experimental unfinished feature. This switch is only for development purposes.)])],
              [COULOMB_COLLISIONS_OPTION="$enableval"], [COULOMB_COLLISIONS_OPTION=no])

AC_ARG_WITH([coulomb-collisions-scheme], [AC_HELP_STRING([--with-coulomb-collisions-scheme],
            [Set coulomb collisions scheme (supported: ta77s [10.1002/ctpp.201700121], sk98 [10.1143/JPSJ.67.4084, 10.1016/j.jcp.2008.03.043] and p12 (default) [10.1063/1.4742167])])],
            [WITH_COULOMB_COLLISIONS_SCHEME="$withval"], [WITH_COULOMB_COLLISIONS_SCHEME="p12"])

AC_ARG_ENABLE([pml], [AC_HELP_STRING([--enable-pml],
              [Enable PML (10.1006/jcph.1994.1159)])],
              [PML_OPTION="$enableval"], [PML_OPTION=yes])

AC_PROG_CXXCPP
AC_LANG(C++)

# set default and customized options
CFLAGS_DEFAULT="-m64 -mcmodel=medium -std=c++17"

# experimental
if test x$EXPERIMENTAL_OPTION == xyes; then
  ENABLE_EXPERIMENTAL=true
  AC_DEFINE_UNQUOTED([ENABLE_EXPERIMENTAL], [true], [Enable experimental options])
  AC_MSG_WARN([You are going to enable experimental features, that was not completely tested. It can cause unexpected behavior!])
else
  ENABLE_EXPERIMENTAL=false
fi

# development
if test x$DEVELOPMENT_OPTION == xyes; then
  ENABLE_DEVELOPMENT=true
  AC_DEFINE_UNQUOTED([ENABLE_DEVELOPMENT], [true], [Enable development options])
  AC_MSG_WARN([You are going to unlock development and unfinished features, that are only for development purposes. It can cause unexpected behavior!])
else
  ENABLE_DEVELOPMENT=false
fi

# debug
if test x$DEBUG_OPTION == xyes; then
  CFLAGS_ADDITIONAL+=" -O0 -Wall -Wextra -ggdb3 -fvar-tracking -ggnu-pubnames -pedantic -time -ftree-vectorizer-verbose=7"
  ENABLE_DEBUG=true
  AC_DEFINE_UNQUOTED([ENABLE_DEBUG], [true], [Use debug option])
else
  ENABLE_DEBUG=false
fi

# fixme
if test x$FIXME_OPTION == xyes; then
  AC_DEFINE_UNQUOTED([ENABLE_FIXME], [true], [Enable FIXME messages])
fi

# profiler
if test x$PROFILER_OPTION = xyes; then
  CFLAGS_ADDITIONAL="$CFLAGS_ADDITIONAL -pg -no-pie"
  LDFLAGS_ADDITIONAL="$LDFLAGS_ADDITIONAL -pg -no-pie"
  AC_DEFINE_UNQUOTED([ENABLE_PROFILER], [true], [Performance debug option])
fi

# singlethread
if test x$SINGLETHREAD_OPTION == xyes; then
  CFLAGS_ADDITIONAL+=" -Wno-unknown-pragmas"
  AC_DEFINE_UNQUOTED([ENABLE_SINGLETHREAD], [true], [Run app in singlethread mode])
else
  AC_OPENMP
  CFLAGS_ADDITIONAL+=" -fopenmp -fopenmp-simd"
fi

# set OpenMP threads to change dynamically
if test x$OMP_DYNAMIC_OPTION == xyes; then
  AC_DEFINE_UNQUOTED([ENABLE_OMP_DYNAMIC], [true], [Set OpenMP threads to be created dynamically])
fi

if test x$MPI_OPTION == xyes; then
  if test x$DEVELOPMENT_OPTION == xyes; then
    AC_MSG_WARN([MPI support is still unfinished feature for development purposes only.])
    ENABLE_MPI=true
    AC_DEFINE_UNQUOTED([ENABLE_MPI], [true], [Enable openMPI support])
  else
    AC_MSG_ERROR([MPI support is still unfinished feature for development purposes only. To ensure, that you know, what you do, please, add --enable-development to your configure arguments])
  fi
fi

# enable IEEE complience
if test x$IEEE_OPTION == xyes; then
  AC_DEFINE_UNQUOTED([ENABLE_IEEE], [true], [Directive to enable acurate IEEE calculations, instead of speed-optimized])
else
  CFLAGS_ADDITIONAL="$CFLAGS_ADDITIONAL -ffast-math -O3 -march=native -funroll-loops -Wno-write-strings -freciprocal-math -fcx-fortran-rules -msse -mavx"
fi

# enable HDF5 data format support
HDF5_OPTION=yes # there is no alternative for now. Let's hardcode it
if test x$HDF5_OPTION = xyes; then
  AC_HAVE_LIBRARY(hdf5, , AC_MSG_WARN([No such library hdf5. You can build it locally after configuration with `make hdf5-install` before build and it includes to build paths automatically]))
  ENABLE_HDF5=true
  AC_DEFINE_UNQUOTED([ENABLE_HDF5], [true], [Use HDF5 data format])
else
  ENABLE_HDF5=false
fi

if test x$PML_OPTION = xyes; then
  AC_DEFINE_UNQUOTED([ENABLE_PML], [true], [Use PML (10.1006/jcph.1994.1159)])
fi

AC_DEFINE_UNQUOTED([REL_LIMIT], [5e7], [Relativistit calculations limit (for boris-adaptive pusher)])

# define a pusher
if test x$WITH_PUSHER == xboris-adaptive; then
  AC_DEFINE_UNQUOTED([SWITCH_PUSHER_BORIS_ADAPTIVE], [true], [Set particles pusher to Boris adaptive])
elif test x$WITH_PUSHER == xboris; then
  AC_DEFINE_UNQUOTED([SWITCH_PUSHER_BORIS], [true], [Set particles pusher to Boris classic])
elif test x$WITH_PUSHER == xboris-relativistic; then
  AC_DEFINE_UNQUOTED([SWITCH_PUSHER_BORIS_RELATIVISTIC], [true], [Set particles pusher to Boris relativistic])
elif test x$WITH_PUSHER == xvay; then
  AC_DEFINE_UNQUOTED([SWITCH_PUSHER_VAY], [true], [Set particles pusher to Vay])
elif test x$WITH_PUSHER == xhc; then
  AC_DEFINE_UNQUOTED([SWITCH_PUSHER_HC], [true], [Set particles pusher to Higuera-Cary])
else
  AC_MSG_ERROR([Pusher $WITH_PUSHER is not supported.])
fi

# define Maxwell solver
AC_DEFINE_UNQUOTED([SWITCH_MAXWELL_SOLVER_YEE], [true], [Yee FDTD Maxwellian solver])

# define a temperature claculation algorithm
if test x$WITH_TEMP_CALC == xweighting; then
  AC_DEFINE_UNQUOTED([SWITCH_TEMP_CALC_WEIGHTING], [true], [Set temperature calculation scheme to weighting])
elif test x$WITH_TEMP_CALC == xcounting; then
  AC_DEFINE_UNQUOTED([SWITCH_TEMP_CALC_COUNTING], [true], [Set temperature calculation scheme to counting])
else
  AC_MSG_ERROR([Temperature calculation scheme $WITH_TEMP_CALC is not supported.])
fi

# define a density claculation algorithm
if test x$WITH_DENSITY_CALC == xweighting; then
  AC_DEFINE_UNQUOTED([SWITCH_DENSITY_CALC_WEIGHTING], [true], [Set density calculation scheme to weighting])
elif test x$WITH_DENSITY_CALC == xcounting; then
  AC_DEFINE_UNQUOTED([SWITCH_DENSITY_CALC_COUNTING], [true], [Set density calculation scheme to counting])
else
  AC_MSG_ERROR([Density calculation scheme $WITH_DENSITY_CALC is not supported.])
fi

# define plasma spatial distribution type
if test x$WITH_PLASMA_SPATIAL == xrandom; then
  AC_DEFINE_UNQUOTED([SWITCH_PLASMA_SPATIAL_RANDOM], [true], [Set plasma spatial distribution type to random])
elif test x$WITH_PLASMA_SPATIAL == xregular; then
  AC_DEFINE_UNQUOTED([SWITCH_PLASMA_SPATIAL_REGULAR], [true], [Set plasma spatial distribution type to regular])
elif test x$WITH_PLASMA_SPATIAL == xcentered; then
  AC_DEFINE_UNQUOTED([SWITCH_PLASMA_SPATIAL_CENTERED], [true], [Set plasma spatial distribution type to centered])
elif test x$WITH_PLASMA_SPATIAL == xflat; then
  AC_DEFINE_UNQUOTED([SWITCH_PLASMA_SPATIAL_FLAT], [true], [Set plasma spatial distribution type to flat])
else
  AC_MSG_ERROR([Plasma spatial distribution type $WITH_PLASMA_SPATIAL is not supported.])
fi

# define plasma spatial distribution type
if test x$WITH_PLASMA_VELOCITY == xthermal; then
  AC_DEFINE_UNQUOTED([SWITCH_PLASMA_VELOCITY_THERMAL], [true], [Set plasma velocity distribution type to thermal])
elif test x$WITH_PLASMA_VELOCITY == xrectangular; then
  AC_DEFINE_UNQUOTED([SWITCH_PLASMA_VELOCITY_RECTANGULAR], [true], [Set plasma velocity distribution type to rectangular])
elif test x$WITH_PLASMA_VELOCITY == xeigen; then
  AC_DEFINE_UNQUOTED([SWITCH_PLASMA_VELOCITY_EIGEN], [true], [Set plasma velocity distribution type to eigen])
else
  AC_MSG_ERROR([Plasma velocity distribution type $WITH_PLASMA_VELOCITY is not supported.])
fi

# define charge concervation scheme
if test x$WITH_CCS == xvb; then
  AC_DEFINE_UNQUOTED([SWITCH_CCS_VB], [true], [Set plasma charge conservation current deposition scheme to Villasenor-Buneman])
elif test x$WITH_CCS == xzigzag; then
  AC_DEFINE_UNQUOTED([SWITCH_CCS_ZIGZAG], [true], [Set plasma charge conservation current deposition scheme to ZigZag])
else
  AC_MSG_ERROR([Plasma charge conservation current deposition scheme $WITH_CCS is not supported.])
fi

if test x$COULOMB_COLLISIONS_OPTION = xyes; then
  if test x$EXPERIMENTAL_OPTION == xyes; then
    AC_DEFINE_UNQUOTED([ENABLE_COULOMB_COLLISIONS], [true], [Enable coulomb collisions])
    AC_MSG_WARN([Coulomb collisions is still experimental and can be unstable.])
  else
    AC_MSG_ERROR([Coulomb collisions is an experimental option. To ensure, that you know, what you do, please, add --enable-experimental to your configure arguments])
  fi

# define coulomb collisions scheme
  if test x$WITH_COULOMB_COLLISIONS_SCHEME == xta77s; then
    AC_DEFINE_UNQUOTED([SWITCH_COULOMB_COLLISIONS_TA77S], [true], [Set coulomb collision scheme to TA77S])
  elif test x$WITH_COULOMB_COLLISIONS_SCHEME == xsk98; then
    AC_DEFINE_UNQUOTED([SWITCH_COULOMB_COLLISIONS_SK98], [true], [Set coulomb collision scheme to SK98])
  elif test x$WITH_COULOMB_COLLISIONS_SCHEME == xp12; then
    AC_DEFINE_UNQUOTED([SWITCH_COULOMB_COLLISIONS_P12], [true], [Set coulomb collision scheme to P12])
  else
    AC_MSG_ERROR([Coulomb collisions scheme $WITH_COULOMB_COLLISIONS_SCHEME is not supported.])
  fi
fi

AC_DEFINE_UNQUOTED([CXXFLAGS], ["$CFLAGS_DEFAULT $CFLAGS_ADDITIONAL"], [define used C(XX)FLAGS])

AC_OUTPUT(Makefile)

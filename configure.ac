#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([PiCoPiC],[esyscmd(git log --pretty=%h -1 | tr -d \\n)],[alexander.vynnyk@protonmail.com],[PiCoPiC])

AC_CANONICAL_TARGET
AC_CANONICAL_HOST

AC_PREREQ([2.72])

AM_INIT_AUTOMAKE([subdir-objects])

AM_PROG_CC_C_O
AC_PROG_CXX
AC_LANG(C++)

AC_CONFIG_HEADERS([include/defines.hpp])

m4_pattern_allow([^AC_])

dnl AC_CANONICAL_TARGET
dnl AC_CANONICAL_HOST
AC_CANONICAL_BUILD
LT_INIT
AC_CONFIG_MACRO_DIRS([m4])

AC_SUBST([CFLAGS])
AC_SUBST([CPPFLAGS])
AC_SUBST([CXXFLAGS])
AC_SUBST([LDFLAGS])
AC_SUBST([AM_CFLAGS])
AC_SUBST([AM_CPPFLAGS])
AC_SUBST([AM_LDFLAGS])

dnl check if CFLAGS CPPFLAGS, CXXFLAGS and LDFLAGS are set
AM_CONDITIONAL(IS_CFLAGS_EMPTY, test "x$CFLAGS" = "x")
AM_CONDITIONAL(IS_CPPFLAGS_EMPTY, test "x$CPPFLAGS" = "x")
AM_CONDITIONAL(IS_CXXFLAGS_EMPTY, test "x$CXXFLAGS" = "x")
AM_CONDITIONAL(IS_LDFLAGS_EMPTY, test "x$LDFLAGS" = "x")

CFLAGS=""
CXXFLAGS=""

dnl
dnl add custom CFLAGS/CPPFLAGS
dnl
dnl === TEST_AND_ADD_CFLAGS(var, flag)
dnl ===   Checks whether $CC supports 'flag' and adds it to 'var'
dnl ===   on success.
AC_DEFUN([TEST_AND_ADD_CFLAGS],
         [SAVED_CFLAGS="$CFLAGS"
          CFLAGS="-Werror $2"
          AC_MSG_CHECKING([whether $CC supports $2])
          dnl Note AC_LANG_PROGRAM([]) uses an old-style main definition.
          AC_COMPILE_IFELSE([AC_LANG_SOURCE([int main(void) { return 0; }])],
                            [AC_MSG_RESULT([yes])]
                            dnl Simply append the variable avoiding a
                            dnl compatibility ifdef for AS_VAR_APPEND as this
                            dnl variable shouldn't grow all that large.
                            [$1="${$1} $2"],
                            [AC_MSG_RESULT([no])])
          CFLAGS="$SAVED_CFLAGS"])

dnl === TEST_AND_ADD_LDFLAGS(var, flag)
dnl ===   Checks whether $CC supports 'flag' and adds it to 'var'
dnl ===   on success.
AC_DEFUN([TEST_AND_ADD_LDFLAGS],
         [SAVED_LDFLAGS="$LDFLAGS"
          LDFLAGS="$2"
          AC_MSG_CHECKING([whether linker supports $2])
          dnl Note AC_LANG_PROGRAM([]) uses an old-style main definition.
          AC_LINK_IFELSE([AC_LANG_SOURCE([int main(void) { return 0; }])],
                         [AC_MSG_RESULT([yes])]
                         dnl Simply append the variable avoiding a
                         dnl compatibility ifdef for AS_VAR_APPEND as this
                         dnl variable shouldn't grow all that large.
                         [$1="${$1} $2"],
                         [AC_MSG_RESULT([no])])
          LDFLAGS="$SAVED_LDFLAGS"])

AM_CONDITIONAL(IS_CPU_X64, test "x$host_cpu" = "xx86_64")
AM_COND_IF([IS_CPU_X64], [ TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-m64]) ])

TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-Ilib/picojson])
TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-Ilib/spdlog/include])
TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-Ilib/HighFive/include])
TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-mcmodel=medium])
TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-D_USE_MATH_DEFINES=1])
TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-DPICOJSON_USE_INT64=true])
TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-DFMT_HEADER_ONLY=1]) dnl for spdlog

dnl create static binary
TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-Wl,-Bstatic])

AM_CFLAGS=$AM_CPPFLAGS
TEST_AND_ADD_CFLAGS([AM_CFLAGS], [-std=gnu20])

TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-std=gnu++20])
AM_CXXFLAGS=$AM_CPPFLAGS

TEST_AND_ADD_LDFLAGS([AM_LDFLAGS], [-ldl])
TEST_AND_ADD_LDFLAGS([AM_LDFLAGS], [-all-static])
TEST_AND_ADD_LDFLAGS([AM_LDFLAGS], [-fopenmp])

ENABLE_DEBUG=no
dnl
dnl 'enable-' style switchers
dnl
dnl enable debug with GDB
AC_ARG_ENABLE([debug],
        AS_HELP_STRING([--enable-debug],
                       [Build with debugging support with GDB @<:@default=<depends on build profile>@:>@]),
              [AM_CONDITIONAL(enable_debug, test x$enableval = xyes)],
              [AM_CONDITIONAL(enable_debug, false)])

AM_COND_IF([enable_debug],
    [
        ENABLE_DEBUG=yes
        AC_DEFINE_UNQUOTED([DEBUG], [1], [Compile with debugger support])
        TEST_AND_ADD_CFLAGS([AM_CFLAGS], [-g3])
        TEST_AND_ADD_CFLAGS([AM_CFLAGS], [-O0])
        TEST_AND_ADD_CFLAGS([AM_CFLAGS], [-ggdb])
        TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-g3])
        TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-O0])
        TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-ggdb])
    ],
    [
        ENABLE_DEBUG=no
        TEST_AND_ADD_CFLAGS([AM_CFLAGS], [-g])
        dnl TEST_AND_ADD_CFLAGS([AM_CFLAGS], [-O3])
        TEST_AND_ADD_CFLAGS([AM_CFLAGS], [-DNDEBUG])
        TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-g])
        dnl TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-O3])
        TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-DNDEBUG])
    ]
    )
AC_SUBST(ENABLE_DEBUG)

dnl performance profiler switcher (for debugging purposes, use with valgrind etc)
AC_ARG_ENABLE([profiler],
        AS_HELP_STRING([--enable-profiler],
                       [Build with profiler support @<:@default=<depends on build profile>@:>@]),
              [AM_CONDITIONAL(enable_profiler, test x$enableval = xyes)],
              [AM_CONDITIONAL(enable_profiler, false)])

AM_COND_IF([enable_profiler],
    [
        TEST_AND_ADD_CFLAGS([AM_CFLAGS], [-pg])
        TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-pg])
        TEST_AND_ADD_LDFLAGS([AM_LDFLAGS], [-pg])
    ])

dnl experimental features switcher
AC_ARG_ENABLE([experimental],
        AS_HELP_STRING([--enable-experimental],
                       [Enable experimental features. @<:@default=no@:>@]),
              [AM_CONDITIONAL(enable_experimental, test x$enableval = xyes)],
              [AM_CONDITIONAL(enable_experimental, false)])

ENABLE_EXPERIMENTAL=no
AM_COND_IF([enable_experimental], [
                                  AC_DEFINE_UNQUOTED([ENABLE_EXPERIMENTAL], [1], [Enable experimental options])
                                  ENABLE_EXPERIMENTAL=yes
                                  AC_MSG_WARN([You are going to enable experimental features, that was not completely tested. It can cause unexpected behavior!])
                                  ])
AC_SUBST(ENABLE_EXPERIMENTAL)

dnl enable FIXME- output
AC_ARG_ENABLE([fixme],
        AS_HELP_STRING([--enable-fixme],
                       [Enable fixme output @<:@default=no@:>@]),
              [AM_CONDITIONAL(enable_fixme, test x$enableval = xyes)],
              [AM_CONDITIONAL(enable_fixme, false)])

AM_COND_IF([enable_fixme], [AC_DEFINE_UNQUOTED([ENABLE_FIXME], [1], [Enable fixme output])])

dnl Keep IEEE compliant rounding and disable hardware acceleration
AC_ARG_ENABLE([ieee],
        AS_HELP_STRING([--enable-ieee],
                       [Keep IEEE compliant rounding and disable hardware acceleration. Decrease calculation speed. @<:@default=no@:>@]),
              [AM_CONDITIONAL(enable_ieee, test x$enableval = xyes)],
              [AM_CONDITIONAL(enable_ieee, false)])

AM_COND_IF([enable_ieee], [
                          AC_DEFINE_UNQUOTED([ENABLE_IEEE], [1], [Keep IEEE compliant rounding and disable hardware acceleration. Decrease calculation speed])
                          ],
                          [
                          TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-O3])
                          TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-ffast-math])
                          TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-march=native])
                          TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-funroll-loops])
                          TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-Wno-write-strings])
                          TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-freciprocal-math])
                          TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-fcx-fortran-rules])
                          TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-msse])
                          TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-mavx])
                          ])

dnl Build application to run in single-thread mode (used in debug build by default)
AC_ARG_ENABLE([singlethread],
        AS_HELP_STRING([--enable-singlethread],
                       [Build application to run in single-thread mode (used in debug build by default). @<:@default=no@:>@]),
              [AM_CONDITIONAL(enable_singlethread, test x$enableval = xyes)],
              [AM_CONDITIONAL(enable_singlethread, false)])

AM_COND_IF([enable_singlethread], [
                                    AC_DEFINE_UNQUOTED([ENABLE_SINGLETHREAD], [1], [Build application to run in single-thread mode])
                                    TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-Wno-unknown-pragmas])
                                  ],
                                  [
                                    AC_OPENMP
                                    TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-fopenmp])
                                    TEST_AND_ADD_CFLAGS([AM_CPPFLAGS], [-fopenmp-simd])
                                  ])

dnl Enable OpenMP dynamic thread creation
AC_ARG_ENABLE([omp-dynamic],
        AS_HELP_STRING([--enable-omp-dynamic],
                       [Enable OpenMP dynamic thread creation. @<:@default=no@:>@]),
              [AM_CONDITIONAL(enable_omp_dynamic, test x$enableval = xyes)],
              [AM_CONDITIONAL(enable_omp_dynamic, false)])

AM_COND_IF([enable_omp_dynamic], [AC_DEFINE_UNQUOTED([ENABLE_OMP_DYNAMIC], [1], [Enable OpenMP dynamic thread creation])])

dnl Enable OpenMPI support
AC_ARG_ENABLE([mpi],
        AS_HELP_STRING([--enable-mpi],
                       [Build application to run in single-thread mode (used in debug build by default). @<:@default=${ENABLE_EXPERIMENTAL}@:>@]),
              [AM_CONDITIONAL(enable_mpi, test x$enableval = xyes && test x$ENABLE_EXPERIMENTAL = xyes)],
              [AM_CONDITIONAL(enable_mpi, false)])

ENABLE_MPI=no
HDF5_ADDITIONAL_OPTIONS=
AM_COND_IF([enable_mpi], [
                         AC_DEFINE_UNQUOTED([ENABLE_MPI], [1], [Enable openMPI support])
                         AC_MSG_WARN([MPI support is still unfinished feature for development purposes only.])
                         ENABLE_MPI=yes
                         HDF5_ADDITIONAL_OPTIONS = --enable-parallel --enable-shared
                         CC = mpicc
                         CXX = mpic++
                         ])
AC_SUBST(ENABLE_MPI)
AC_SUBST(HDF5_ADDITIONAL_OPTIONS)

dnl Build application with hdf5 library/data format support
AC_ARG_ENABLE([hdf5],
        AS_HELP_STRING([--enable-hdf5],
                       [Build application with hdf5 library and data format support. @<:@default=no@:>@]),
              [AM_CONDITIONAL(enable_hdf5, test x$enableval = xyes)],
              [AM_CONDITIONAL(enable_hdf5, true)],
              )

AM_COND_IF([enable_hdf5], [
                          AC_DEFINE_UNQUOTED([ENABLE_HDF5], [1], [Enable openHDF5 support])
                          AC_SEARCH_LIBS([H5open], [hdf5],
                          [TEST_AND_ADD_LDFLAGS([AM_LDFLAGS], [-lhdf5])],
                          [AC_MSG_WARN([No such library hdf5. You can build it locally after configuration with `make hdf5-install` before build and it includes to build paths automatically])]
                          )])


dnl Apply external fields
AC_ARG_ENABLE([external-fields],
        AS_HELP_STRING([--enable-external-fields],
                       [Apply external fields. @<:@default=no@:>@]),
              [AM_CONDITIONAL(enable_external_fields, test x$enableval = xyes && test x$ENABLE_EXPERIMENTAL = xyes)],
              [AM_CONDITIONAL(enable_external_fields, false)])

AM_COND_IF([enable_external_fields], [AC_DEFINE_UNQUOTED([ENABLE_EXTERNAL_FIELDS], [1], [Apply external fields])])

dnl Enable coulomb collisions. (WARNING! This is an experimental unfinished feature. This switch is only for development purposes.)
AC_ARG_ENABLE([coulomb-collisions],
        AS_HELP_STRING([--enable-coulomb-collisions],
                       [Enable coulomb collisions. (WARNING! This is an experimental unfinished feature. This switch is only for development purposes). @<:@default=no@:>@]),
              [AM_CONDITIONAL(enable_coulomb_collisions, test x$enableval = xyes && test x$ENABLE_EXPERIMENTAL = xyes)],
              [AM_CONDITIONAL(enable_coulomb_collisions, false)])

AM_COND_IF([enable_coulomb_collisions], [AC_DEFINE_UNQUOTED([ENABLE_COULOMB_COLLISIONS], [1], [Enable coulomb collisions. (WARNING! This is an experimental unfinished feature. This switch is only for development purposes)])])

dnl Enable PML (10.1006/jcph.1994.1159)
AC_ARG_ENABLE([pml],
        AS_HELP_STRING([--enable-pml],
                       [Enable PML (10.1006/jcph.1994.1159). @<:@default=no@:>@]),
              [AM_CONDITIONAL(enable_pml, test x$enableval = xyes)],
              [AM_CONDITIONAL(enable_pml, false)])

AM_COND_IF([enable_pml], [AC_DEFINE_UNQUOTED([ENABLE_PML], [1], [Use PML (10.1006/jcph.1994.1159)])])

dnl
dnl 'with-' style argument parsers
dnl

dnl Set list of documentation formats to be generated
AC_ARG_WITH([doc-formats], [AS_HELP_STRING([--with-doc-formats],[Set list of documentation formats to be generated])],
            [WITH_DOXYGEN_FORMATS="$withval"], [WITH_DOXYGEN_FORMATS="latex html rtf"])
AC_SUBST(WITH_DOXYGEN_FORMATS)

dnl Set motion equation solver (aka "pusher")
AC_ARG_WITH([pusher],
            [AS_HELP_STRING([--with-pusher], [Set particles pusher (supported: boris (non-relativistic), boris-adaptive [default], boris-relativistic, vay, hc [higuera-cary] ).])],
            [case "${withval}" in
              boris) AC_DEFINE_UNQUOTED([SWITCH_PUSHER_BORIS], [1], [Set particles pusher to Boris classic]) ;;
              boris-adaptive) AC_DEFINE_UNQUOTED([SWITCH_PUSHER_BORIS_ADAPTIVE], [1], [Set particles pusher to Boris adaptive]) ;;
              boris-relativistic) AC_DEFINE_UNQUOTED([SWITCH_PUSHER_BORIS_RELATIVISTIC], [1], [Set particles pusher to Boris relativistic]) ;;
              vay) AC_DEFINE_UNQUOTED([SWITCH_PUSHER_VAY], [1], [Set particles pusher to Vay]) ;;
              hc) AC_DEFINE_UNQUOTED([SWITCH_PUSHER_HC], [1], [Set particles pusher to Higuera-Cary]) ;;
              *) AC_MSG_ERROR([bad value ${withval} for --with-pusher. Should be 'boris', 'boris-relativistic', 'boris-adaptive', 'vay' or 'hc']) ;;
            esac],
            dnl default behavior is boris-adaptive
            [ AC_DEFINE_UNQUOTED([SWITCH_PUSHER_BORIS_ADAPTIVE], [1], [Set particles pusher to Boris adaptive]) ])

dnl Set plasma particles spatial distribution type
AC_ARG_WITH([plasma-spatial-distribution],
            [AS_HELP_STRING([--with-plasma-spatial-distribution],[Set plasma particles spatial distribution type (supported: flat (default), random, regular, centered).])],
            [case "${withval}" in
              flat) AC_DEFINE_UNQUOTED([SWITCH_PLASMA_SPATIAL_FLAT], [1], [Set plasma spatial distribution type to flat]) ;;
              random) AC_DEFINE_UNQUOTED([SWITCH_PLASMA_SPATIAL_RANDOM], [1], [Set plasma spatial distribution type to random]) ;;
              regular) AC_DEFINE_UNQUOTED([SWITCH_PLASMA_SPATIAL_REGULAR], [1], [Set plasma spatial distribution type to regular]) ;;
              centered) AC_DEFINE_UNQUOTED([SWITCH_PLASMA_SPATIAL_CENTERED], [1], [Set plasma spatial distribution type to centered]) ;;
              *) AC_MSG_ERROR([bad value ${withval} for --with-plasma-spatial-distribution. Should be 'flat', 'random', 'regular', or 'centered']) ;;
            esac],
            dnl default behavior is flat
            [ AC_DEFINE_UNQUOTED([SWITCH_PLASMA_SPATIAL_FLAT], [1], [Set plasma spatial distribution type to flat (by default)]) ])

dnl Set plasma particles velocity distribution type
AC_ARG_WITH([plasma-velocity-distribution],
            [AS_HELP_STRING([--with-plasma-velocity-distribution],[Set plasma particles velocity distribution type (supported: thermal (default), rectangular, eigen).])],
            [case "${withval}" in
              thermal) AC_DEFINE_UNQUOTED([SWITCH_PLASMA_VELOCITY_THERMAL], [1], [Set plasma velocity distribution type to thermal]) ;;
              rectangular) AC_DEFINE_UNQUOTED([SWITCH_PLASMA_VELOCITY_RECTANGULAR], [1], [Set plasma velocity distribution type to rectangular]) ;;
              eigen) AC_DEFINE_UNQUOTED([SWITCH_PLASMA_VELOCITY_EIGEN], [1], [Set plasma velocity distribution type to eigen]) ;;
              *) AC_MSG_ERROR([bad value ${withval} for --with-plasma-velocity-distribution. Should be 'thermal', 'rectangular', or 'eigen']) ;;
            esac],
            dnl default behavior is flat
            [ AC_DEFINE_UNQUOTED([SWITCH_PLASMA_VELOCITY_THERMAL], [1], [Set plasma velocity distribution type to thermal]) ])

dnl Set charge conservation scheme of current solver (supported: vb [villasenor-buneman] (default) and zigzag)
AC_ARG_WITH([current-solver],
            [AS_HELP_STRING([--with-current-solver],
                        [Set charge conservation scheme of current solver (supported: vb [villasenor-buneman] (default) and zigzag).])],
            [case "${withval}" in
              vb) AC_DEFINE_UNQUOTED([SWITCH_CCS_VB], [1], [Set plasma charge conservation current deposition scheme to Villasenor-Buneman]) ;;
              zigzag) AC_DEFINE_UNQUOTED([SWITCH_CCS_ZIGZAG], [1], [Set plasma charge conservation current deposition scheme to ZigZag]) ;;
              *) AC_MSG_ERROR([bad value ${withval} for --with-current-solver. Should be 'vb', or 'zigzag']) ;;
            esac],
            dnl default behavior is flat
            [ AC_DEFINE_UNQUOTED([SWITCH_CCS_ZIGZAG], [1], [Set plasma charge conservation current deposition scheme to ZigZag]) ])

dnl Set coulomb collisions scheme (supported: ta77s [10.1002/ctpp.201700121], sk98 [10.1143/JPSJ.67.4084, 10.1016/j.jcp.2008.03.043] and p12 (default) [10.1063/1.4742167])
AC_ARG_WITH([coulomb-collisions-scheme],
            [AS_HELP_STRING([--with-coulomb-collisions-scheme],
            [Set coulomb collisions scheme (supported: ta77s [10.1002/ctpp.201700121], sk98 [10.1143/JPSJ.67.4084, 10.1016/j.jcp.2008.03.043] and p12 (default) [10.1063/1.4742167]).])],
            [case "${withval}" in
              ta77s) AC_DEFINE_UNQUOTED([SWITCH_COULOMB_COLLISIONS_TA77S], [1], [Set coulomb collision scheme to TA77S]) ;;
              sk98) AC_DEFINE_UNQUOTED([SWITCH_COULOMB_COLLISIONS_SK98], [1], [Set coulomb collision scheme to SK98]) ;;
              p12) AC_DEFINE_UNQUOTED([SWITCH_COULOMB_COLLISIONS_P12], [1], [Set coulomb collision scheme to P12]) ;;
              *) AC_MSG_ERROR([bad value ${withval} for --with-coulomb-collisions-scheme. Should be 'ta77s', 'sk98', or 'p12']) ;;
            esac],
            dnl default behavior is flat
            [ AC_DEFINE_UNQUOTED([SWITCH_COULOMB_COLLISIONS_P12], [1], [Set coulomb collision scheme to P12 (default)]) ])

AC_DEFINE_UNQUOTED([PACKAGE_DEPS], ["gcc@14+;picojson@1.3.1-dev;python@3+;hdf5@${HDF5_VERSION};HighFive@v2.9.0;spdlog@v1.15.0+;googletest@v1.17+"],
                                   [define used C(XX)FLAGS])

AC_DEFINE_UNQUOTED([REL_LIMIT], [5e7], [Relativistic calculations limit (for boris-adaptive pusher)])

AC_DEFINE_UNQUOTED([CXXFLAGS], ["$CXXFLAGS"], [CPP/CXXFLAGS current application built with])

dnl define Maxwell solver (TODO: implement solvers other than Yee)
AC_DEFINE_UNQUOTED([SWITCH_MAXWELL_SOLVER_YEE], [true], [Yee FDTD Maxwellian solver])

HDF5_VERSION=1.14.6
AC_SUBST(HDF5_VERSION)

dnl
dnl configure doxygen
dnl
AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN"; then
   AC_MSG_WARN([Doxygen not found - continue without Doxygen support])
fi
AC_CHECK_PROGS([DOT], [dot])
if test -z "$DOT"; then
   AC_MSG_ERROR([Doxygen needs dot, please install dot first])
fi
AC_CHECK_PROGS([PDFLATEX], [pdflatex])
if test -z "$PDFLATEX"; then
   AC_MSG_ERROR([Doxygen needs pdflatex program, it is part of TeX http://www.tug.org/texlive/acquire-netinstall.html])
fi
AC_CHECK_PROGS([CONVERT], [convert])
if test -z "$CONVERT"; then
   AC_MSG_ERROR([Doxygen needs convert tool from imagemagick package, please install imagemagick first])
fi
AM_CONDITIONAL([HAVE_DOXYGEN], [test -n "$DOXYGEN"])
AM_COND_IF([HAVE_DOXYGEN], [AC_CONFIG_FILES([doc/Doxyfile])])
AM_COND_IF([HAVE_DOXYGEN], [AC_CONFIG_FILES([doc/Makefile])])

AC_CONFIG_FILES([Makefile test/Makefile])
AC_OUTPUT

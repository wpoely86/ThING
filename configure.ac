#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.66])
AC_INIT([ThING], [1.0], [sebastianwouters@gmail.com])
AM_INIT_AUTOMAKE([foreign -Wall])
AC_CONFIG_SRCDIR([lib/Gauss.cpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

m4_define([ThING_version_major], [1])
m4_define([ThING_version_minor], [0])
m4_define([ThING_version],[ThING_version_major.ThING_version_minor])

# Checks for programs.
AC_PROG_CXX
#AC_PROG_CC
AC_PROG_CXX_C_O
AC_LANG_PUSH([C++])

LT_INIT
LT_LANG([C++])

# Checks for libraries.
AC_SEARCH_LIBS([dgemm_],[blas],[],AC_MSG_ERROR([Unable to find libblas]))
AC_SEARCH_LIBS([dsyev_],[lapack],[],AC_MSG_ERROR([Unable to find liblapack]))
AC_SEARCH_LIBS([__gmpz_fac_ui],[gmp],[],AC_MSG_ERROR([Unable to find libgmp]))
AC_SEARCH_LIBS([__gmpz_init],[gmpxx],[],AC_MSG_ERROR([Unable to find libgmpxx]))
AC_SEARCH_LIBS([pow],[m],[],AC_MSG_ERROR([Unable to find libm]))

# Checks for header files.
AC_HEADER_ASSERT
AC_CHECK_HEADER([gmpxx.h],[],AC_MSG_ERROR([Could not find gmpxx.h]))

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL

# Checks for library functions.
AC_CHECK_FUNC([pow],[],AC_MSG_ERROR([Could not find pow() function ]))
AC_CHECK_FUNC([sqrt],[],AC_MSG_ERROR([Could not find sqrt() function]))

AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test -z "$DOXYGEN";
then AC_MSG_WARN([Doxygen not found - continuing without Doxygen support])
fi

AC_MSG_CHECKING([whether to build with debug information])
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug],
        [enable debug data generation (def=no)])],
    [debugit="$enableval"],
    [debugit=no])
AC_MSG_RESULT([$debugit])

if test x"$debugit" = x"yes"; then
    AC_DEFINE([DEBUG],[],[Debug Mode])
    AM_CXXFLAGS="$AM_CXXFLAGS -g -Wall -Wno-uninitialized"
else
    AC_DEFINE([NDEBUG],[],[No-debug Mode])
    AM_CXXFLAGS="$AM_CXXFLAGS -O2"
fi

AC_SUBST([AM_CXXFLAGS])


AC_CONFIG_FILES([Makefile lib/Makefile doc/Makefile tools/Makefile lib/libThING.pc basissets/Makefile examples/Makefile])


AC_OUTPUT

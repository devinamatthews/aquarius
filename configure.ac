AC_INIT(aquarius, 1.0, dmatthews@utexas.edu, aquarius, http://www.github.com/devinamatthews/aquarius)
AM_INIT_AUTOMAKE([foreign silent-rules subdir-objects])
AM_SILENT_RULES([yes])
AM_MAINTAINER_MODE([disable])

AC_CONFIG_MACRO_DIR([m4])

AM_CONFIG_HEADER(config.h)
AH_TOP([
#ifndef _CONFIG_H_
#define _CONFIG_H_
])
AH_BOTTOM([
#ifndef DEBUG
#define NDEBUG
#endif

#endif
])

AX_PROG_CC_MPI([true], [], [AC_MSG_FAILURE([No MPI C compiler found.])])
AX_PROG_CXX_MPI([true], [], [AC_MSG_FAILURE([No MPI C++ compiler found.])])
AX_CXX_COMPILE_STDCXX_11([noext])
AC_LANG([C++])
AC_C_RESTRICT
AC_OPENMP

CPPFLAGS="$CPPFLAGS $OPENMP_CXXFLAGS"
LDFLAGS="$LDFLAGS $OPENMP_CXXFLAGS"

AX_BLAS([], [AC_MSG_FAILURE([Could not find BLAS library.])])
AX_LAPACK([], [AC_MSG_FAILURE([Could not find LAPACK library.])])

topdir="\"`(cd $srcdir && pwd)`\""
AC_DEFINE_UNQUOTED([TOPDIR], [$topdir], [The top source directory.])

AC_DEFINE([LOWER_NO_UNDERSCORE], [1], [FORTRAN name-mangling scheme.])
AC_DEFINE([LOWER_UNDERSCORE], [2], [FORTRAN name-mangling scheme.])
AC_DEFINE([UPPER_NO_UNDERSCORE], [3], [FORTRAN name-mangling scheme.])
AC_DEFINE([UPPER_UNDERSCORE], [4], [FORTRAN name-mangling scheme.])

AC_F77_FUNC(sgemm)
AS_IF([test $sgemm = sgemm], [AC_DEFINE([F77_NAME], [LOWER_NO_UNDERSCORE], [Type of FORTRAN name-mangling.])])
AS_IF([test $sgemm = sgemm_], [AC_DEFINE([F77_NAME], [LOWER_UNDERSCORE], [Type of FORTRAN name-mangling.])])
AS_IF([test $sgemm = SGEMM], [AC_DEFINE([F77_NAME], [UPPER_NO_UNDERSCORE], [Type of FORTRAN name-mangling.])])
AS_IF([test $sgemm = SGEMM_], [AC_DEFINE([F77_NAME], [UPPER_UNDERSCORE], [Type of FORTRAN name-mangling.])])

#
# Configure LAWrap
#
AQ_WITH_PACKAGE([LAWrap], [https://github.com/devinamatthews/lawrap],
                [blas.h lapack.h], [.])

#
# Configure MArray
#
AQ_WITH_PACKAGE([MArray], [https://github.com/devinamatthews/marray],
                [marray.hpp], [.])

#
# Configure MPIWrap
#
AQ_WITH_PACKAGE([MPIWrap], [https://github.com/devinamatthews/mpiwrap],
                [mpiwrap.hpp], [])

#
# Configure CTF
#
AQ_WITH_PACKAGE([CTF], [https://github.com/solomonik/ctf],
                [ctf.hpp], [],
                [CTF_linked], [-lctf])
AM_CONDITIONAL([CTF_IS_LOCAL], [test x"$ctf_INCLUDES" = x" -Isrc/external/ctf/include"])

#
# Check for Elemental
#
AQ_WITH_PACKAGE([Elemental], [],
                [El.hpp], [],
                [ElInitialize], [-lEl -lpmrrr -lmetis -lkiss_fft])

#
# Check for libint2
#
AQ_WITH_PACKAGE([libint2], [],
                [libint2.h], [],
                [libint2_static_init], [-lint2])
                
mkdir -p bin

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

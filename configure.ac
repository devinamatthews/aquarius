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

AX_BLAS([], [AC_MSG_FAILURE([Could not find BLAS library.])])
AX_LAPACK([], [AC_MSG_FAILURE([Could not find LAPACK library.])])

topdir=\"`(cd $srcdir && pwd)`\"
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
# Configure CTF
#
AC_ARG_WITH([ctf],
            [AS_HELP_STRING([--with-ctf=<dir>],
                            [Specify the location of CTF. Otherwise, it will be downloaded.])],
            [], [with_ctf=./ctf])
            
ctfdir=$with_ctf
if test x"$ctfdir" = "x./ctf" && ! test -d $ctfdir; then
    if (command -v git >/dev/null 2>&1); then
        echo -e "Downloading CTF from external Git repository...\n"
        if ! git clone https://github.com/solomonik/ctf $ctfdir; then
            echo -e "Could not download CTF repository.\n"
            exit 1
        fi
        echo
    else
        echo 'An existing installation of CTF is not found and cannot be downloaded.'
        echo 'Install Git (www.git-scm.com) or download CTF and specify the location.'
        echo
        exit 1
    fi
fi
AC_SUBST([ctfdir], [$ctfdir])

#
# Check for Elemental
#
AC_ARG_WITH([elemental],
            [AS_HELP_STRING([--with-elemental=<dir>],
                            [Enable the use of the distributed-memory Elemental library.])],
            [], [with_elemental=no])
AC_ARG_WITH([elemental-libs],
            [AS_HELP_STRING([--with-elemental-libs=<lib>],
                            [Specify the libraries to link to use Elemental.])],
            [], [with_elemental_libs="-lEl -lpmrrr"])

elemental_libs=$with_elemental_libs
if test x"$with_elemental" != xno; then
    LDFLAGS_SAVE=$LDFLAGS
    CPPFLAGS_SAVE=$CPPFLAGS
    LDFLAGS=
    CPPFLAGS=$CXXFLAGS
        
    if test x"$with_elemental" != xyes; then
        LDFLAGS+=" -L$with_elemental/lib -L$with_elemental/lib64"
        CPPFLAGS+=" -I$with_elemental/include"
    fi
    
    AC_CHECK_HEADER([El.hpp], [], [AC_MSG_FAILURE([Elemental requested, but could not find El.hpp.])])
    AC_CHECK_LIB([El], [ElInitialize], [true], [AC_MSG_FAILURE([Elemental requested, but could not find -lEl.])], [$elemental_libs])
    
    LDFLAGS=$LDFLAGS_SAVE
    CPPFLAGS=$CPPFLAGS_SAVE
        
    have_elemental=yes
    
    if test x"$with_elemental" != xyes; then
        AC_SUBST([elemental_LIBS], ["-L$with_elemental/lib -L$with_elemental/lib64 $elemental_libs"])
        AC_SUBST([elemental_INCLUDES], ["-I$with_elemental/include"])
    else
        AC_SUBST([elemental_LIBS], ["$elemental_libs"])
        AC_SUBST([elemental_INCLUDES], [])
    fi
fi

AS_IF([test x"$have_elemental" = xyes], [AC_DEFINE([HAVE_ELEMENTAL], [1], [Define if Elemental is to be used.])])
AM_CONDITIONAL([HAVE_ELEMENTAL], [test x"$have_elemental" = xyes])

#
# Check for libint2
#
AC_ARG_WITH([libint2],
            [AS_HELP_STRING([--with-libint2=<dir>],
                            [Enable the use of the libint2 atomic integral library.])],
            [], [with_libint2=no])

if test x"$with_libint2" != xno; then
    LDFLAGS_SAVE=$LDFLAGS
    CPPFLAGS_SAVE=$CPPFLAGS
    LDFLAGS=
    CPPFLAGS=$CXXFLAGS
        
    if test x"$with_libint2" != xyes; then
        LDFLAGS+=" -L$with_libint2/lib -L$with_libint2/lib64"
        CPPFLAGS+=" -I$with_libint2/include -I$with_libint2/include/libint2"
    fi
    
    AC_CHECK_HEADER([libint2.h], [], [AC_MSG_FAILURE([libint2 requested, but could not find libint2.h.])])
    AC_CHECK_LIB([int2], [libint2_static_init], [true], [AC_MSG_FAILURE([libint2 requested, but could not find -lint2.])])
    
    LDFLAGS=$LDFLAGS_SAVE
    CPPFLAGS=$CPPFLAGS_SAVE
        
    have_libint2=yes
    
    if test x"$with_libint2" != xyes; then
        AC_SUBST([libint2_LIBS], ["-L$with_libint2/lib -L$with_libint2/lib64 -lint2"])
        AC_SUBST([libint2_INCLUDES], ["-I$with_libint2/include -I$with_libint2/include/libint2"])
    else
        AC_SUBST([libint2_LIBS], ["-lint2"])
        AC_SUBST([libint2_INCLUDES], [])
    fi
fi

AS_IF([test x"$have_libint2" = xyes], [AC_DEFINE([HAVE_LIBINT2], [1], [Define if libint2 is to be used.])])
AM_CONDITIONAL([HAVE_LIBINT2], [test x"$have_libint2" = xyes])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

#ifndef _AQUARIUS_UTIL_LAWRAP_HPP_
#define _AQUARIUS_UTIL_LAWRAP_HPP_

#if F77_NAME == UPPER_NO_UNDERSCORE
#define FC_FUNC(name,NAME) NAME
#elif F77_NAME == UPPER_UNDERSCORE
#define FC_FUNC(name,NAME) NAME##_
#elif F77_NAME == LOWER_NO_UNDERSCORE
#define FC_FUNC(name,NAME) name
#elif F77_NAME == LOWER_UNDERSCORE
#define FC_FUNC(name,NAME) name##_
#else
#error "Unsupported Fortran naming convention"
#endif

#include "fortran.h"
#include "blas.h"
#include "lapack.h"

namespace aquarius
{
using namespace LAWrap;
using LAWrap::copy;
}

#endif

#ifndef _AQUARIUS_UTIL_GLOBAL_HPP_
#define _AQUARIUS_UTIL_GLOBAL_HPP_

#if !(defined(__GXX_EXPERIMENTAL_CXX0X__) || _MSC_VER >= 1600 || __cplusplus >= 201103l)
#error "A C++11-capable compiler is required."
#endif

#include <omp.h>

#include "config.h"

#include "lawrap.hpp"

#define MARRAY_BLAS_PROTOTYPED
#include "marray.hpp"

namespace aquarius
{
using namespace MArray;
}

#include "stl_ext.hpp"
#include "math_ext.hpp"
#include "distributed.hpp"

#endif

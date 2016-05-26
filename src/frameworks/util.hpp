#ifndef _AQUARIUS_FRAMEWORKS_UTIL_HPP_
#define _AQUARIUS_FRAMEWORKS_UTIL_HPP_

#if !(defined(__GXX_EXPERIMENTAL_CXX0X__) || _MSC_VER >= 1600 || __cplusplus >= 201103l)
#error "A C++11-capable compiler is required."
#endif

#include <omp.h>

#include "config.h"

#include "util/stl_ext.hpp"
#include "util/math_ext.hpp"
#include "util/distributed.hpp"

#endif

#ifndef _AQUARIUS_STL_EXT_HPP_
#define _AQUARIUS_STL_EXT_HPP_

#include <cassert>
#include <cctype>
#include <cmath>
#include <complex>
#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <array>
#include <deque>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <iomanip>
#include <istream>
#include <iterator>
#include <limits>
#include <list>
#include <locale>
#include <map>
#include <memory>
#include <ostream>
#include <random>
#include <set>
#include <sstream>
#include <stack>
#include <tuple>
#include <utility>
#include <regex>
#include <stdexcept>

#include "algorithm.hpp"
#include "any.hpp"
#include "complex.hpp"
#include "cosort.hpp"
#include "global_ptr.hpp"
#include "iostream.hpp"
#include "ptr_list.hpp"
#include "ptr_vector.hpp"
#include "string.hpp"
#include "type_traits.hpp"
#include "vector.hpp"
#include "zip.hpp"

#define INSTANTIATE_SPECIALIZATIONS(name) \
template class name<double>;

#define INSTANTIATE_SPECIALIZATIONS_2(name,extra1) \
template class name<double,extra1>;

#define INSTANTIATE_SPECIALIZATIONS_3(name,extra1,extra2) \
template class name<double,extra1,extra2>;

#define CONCAT(...) __VA_ARGS__

#define STRINGIZE(...) #__VA_ARGS__

#ifdef DEBUG

#define DPRINTF(...) \
do \
{ \
    printf("%s(%d): ", __FILE__, __LINE__); \
    printf(__VA_ARGS__); \
} while (0)

#define DPRINTFC(...) \
do \
{ \
    printf(__VA_ARGS__); \
} while (0)

#else

#define DPRINTF(...)

#define DPRINTFC(...)

#endif

namespace aquarius
{
    using std::string;
    using std::tuple;
    using std::map;
    using std::list;
    using std::set;
    using std::array;
    using std::vector;
    using std::deque;
    using std::pair;
    using std::make_pair;
    using std::make_tuple;
    using std::get;
    using std::tuple_element;
    using std::numeric_limits;

    using std::min;
    using std::max;
    using std::abs;
    using std::copy;
    using std::copy_n;
    using std::fill;
    using std::fill_n;
    using std::sort;
    using std::swap;
    using std::remove;
    using std::rotate;

    using std::unique_ptr;
    using std::shared_ptr;
    using std::make_shared;

    using std::runtime_error;
    using std::logic_error;

    using std::move;
    using std::forward;

    using std::iterator;
    using std::iterator_traits;
    using std::reverse_iterator;
    using std::forward_iterator_tag;
    using std::bidirectional_iterator_tag;
    using std::random_access_iterator_tag;

    using std::type_info;
    using std::true_type;
    using std::false_type;
    using std::unary_function;
    using std::binary_function;

    using std::ostream;
    using std::istream;
    using std::stringstream;
    using std::ostringstream;
    using std::istringstream;
    using std::fstream;
    using std::ofstream;
    using std::ifstream;
    using std::cout;
    using std::cin;
    using std::cerr;
    using std::endl;
    using std::setprecision;
    using std::scientific;
    //using std::defaultfloat; //no GCC support until GCC5 :(
    //using std::hexfloat;
    using std::fixed;
    using std::skipws;
    using std::showpos;
    using std::showbase;
    using std::noshowbase;
    using std::noshowpos;
    using std::showpoint;
    using std::noshowpoint;
    using std::dec;
    using std::oct;
    using std::hex;
    using std::setfill;
    using std::setw;
    using std::streamsize;
    using std::streambuf;

    using std::toupper;
    using std::tolower;

    using std::exception;
    using std::logic_error;
    using std::runtime_error;
    using std::out_of_range;

    using std::regex_search;
    using std::regex_match;
    using std::smatch;
    using std::cmatch;
    using std::regex;

    using std::uniform_int_distribution;
    using std::uniform_real_distribution;
    using std::default_random_engine;
    using std::mt19937;

    using namespace stl_ext;
    using stl_ext::sort;
}

#endif

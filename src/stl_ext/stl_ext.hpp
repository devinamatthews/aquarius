/* Copyright (c) 2013, Devin Matthews
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *      * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL DEVIN MATTHEWS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE. */

#ifndef _AQUARIUS_STL_EXT_HPP_
#define _AQUARIUS_STL_EXT_HPP_

#include <vector>
#include <string>
#include <sstream>
#include <ostream>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <cctype>

namespace std
{

template<typename T> std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    os << "[";
    if (!v.empty()) os << v[0];
    for (int i = 1;i < v.size();i++) os << ", " << v[i];
    os << "]";
    return os;
}

#ifndef CPLUSPLUS11

template<typename I1, typename I2, typename Pred>
I2 copy_if(I1 begin, I1 end, I2 result, Pred pred)
{
    return remove_copy_if(begin, end, result, not1(pred));
}

#endif

template<typename T> std::string str(const T& t)
{
    std::ostringstream oss;
    oss << t;
    return oss.str();
}

template<typename T> std::vector<T> operator+(const std::vector<T>& v1, const std::vector<T>& v2)
{
    using std::copy;
    std::vector<T> r(v1.size() + v2.size());
    copy(v1.begin(), v1.end(), r.begin());
    copy(v2.begin(), v2.end(), r.begin() + v1.size());
    return r;
}

template<typename T> std::vector<T> operator+(const std::vector<T>& v, const T& t)
{
    using std::copy;
    std::vector<T> r(v.size()+1);
    copy(v.begin(), v.end(), r.begin());
    r[v.size()] = t;
    return r;
}

template<typename T> std::vector<T>& operator+=(std::vector<T>& v1, const std::vector<T>& v2)
{
    v1.insert(v1.end(), v2.begin(), v2.end());
    return v1;
}

template<typename T> std::vector<T>& operator+=(std::vector<T>& v, const T& t)
{
    v.insert(v.end(), t);
    return v;
}

template<typename T> std::vector<T> slice(const std::vector<T>& v, int e1, int e2)
{
    return std::vector<T>(v.begin()+e1, v.begin()+e2);
}

template<typename T> std::vector<T> slice(const std::vector<T>& v, int e1)
{
    return std::vector<T>(v.begin()+e1, v.end());
}

template<typename T, class Predicate> std::vector<T>& filter(std::vector<T>& v, Predicate pred)
{
    typename std::vector<T>::iterator i1, i2;

    i1 = v.begin();
    for (i2 = v.begin();i2 != v.end();++i2)
    {
        if (pred(*i2))
        {
            if (i1 != i2) *i1 = *i2;
            ++i1;
        }
    }

    v.resize(i1-v.begin());
    return v;
}

template<typename T, typename U, class Functor> std::vector<U> map(std::vector<T>& v, Functor f)
{
    std::vector<U> v2();

    typename std::vector<T>::const_iterator i;

    for (i = v.begin();i != v.end();++i)
    {
        v2.push_back(f(*i));
    }

    return v2;
}

template<typename T, class Predicate> std::vector<T> filter_copy(const std::vector<T>& v, Predicate pred)
{
    typename std::vector<T> v2(v.size());
    typename std::vector<T>::iterator i1;
    typename std::vector<T>::const_iterator i2;

    i1 = v2.begin();
    for (i2 = v.begin();i2 != v.end();++i2)
    {
        if (pred(*i2))
        {
            *i1 = *i2;
            ++i1;
        }
    }

    v2.resize(i1-v2.begin());
    return v2;
}

template<typename T> std::vector<T>& uniq(std::vector<T>& v)
{
    typename std::vector<T>::iterator i1;

    std::sort(v.begin(), v.end());
    i1 = std::unique(v.begin(), v.end());
    v.resize(i1-v.begin());

    return v;
}

template<typename T> std::vector<T> uniq_copy(const std::vector<T>& v)
{
    typename std::vector<T> v2(v.size());
    typename std::vector<T>::iterator i1;

    std::sort(v2.begin(), v2.end());
    i1 = std::unique(v2.begin(), v2.end());
    v2.resize(i1-v2.begin());

    return v2;
}

template<typename T> std::vector<T> intersection(const std::vector<T>& v1, const std::vector<T>& v2)
{
    typename std::vector<T> v(v1.size()+v2.size());
    typename std::vector<T>::iterator end;

    end = std::set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), v.begin());
    v.resize(end-v.begin());

    return v;
}

template<typename T> std::vector<T>& exclude(std::vector<T>& v1, const std::vector<T>& v2)
{
    typename std::vector<T> v3(v2);
    typename std::vector<T>::iterator i1, i2, i3;

    std::sort(v1.begin(), v1.end());
    std::sort(v3.begin(), v3.end());

    i1 = i2 = v1.begin();
    i3 = v3.begin();
    while (i1 != v1.end())
    {
        if (i3 == v3.end() || *i1 < *i3)
        {
            *i2 = *i1;
            ++i1;
            ++i2;
        }
        else if (*i3 < *i1)
        {
            ++i3;
        }
        else
        {
            ++i1;
        }
    }
    v1.resize(i2-v1.begin());

    return v1;
}

template<typename T> std::vector<T> exclude_copy(const std::vector<T>& v1, const std::vector<T>& v2)
{
    typename std::vector<T> v3(v1);
    return exclude(v3, v2);
}

template<typename T, typename U> std::vector<T>& mask(std::vector<T>& v, const std::vector<U>& mask)
{
    typename std::vector<T>::iterator i1;
    typename std::vector<T>::const_iterator i2;
    typename std::vector<U>::const_iterator i3;

    i1 = v.begin();
    i2 = v.begin();
    i3 = mask.begin();
    for (;i2 != v.end();++i2,++i3)
    {
        if (*i3)
        {
            using std::swap;
            swap(*i1, *i2);
            ++i1;
        }
    }

    v.resize(i1-v.begin());
    return v;
}

template<typename T, typename U> std::vector<T> mask_copy(const std::vector<T>& v, const std::vector<U>& mask)
{
    typename std::vector<T> v2(v.size());
    typename std::vector<T>::iterator i1;
    typename std::vector<T>::const_iterator i2;
    typename std::vector<U>::const_iterator i3;

    i1 = v2.begin();
    i2 = v.begin();
    i3 = mask.begin();
    for (;i2 != v.end();++i2,++i3)
    {
        if (*i3)
        {
            *i1 = *i2;
            ++i1;
        }
    }

    v2.resize(i1-v2.begin());
    return v2;
}

template<typename T, typename U> bool compareFirst(std::pair<T,U> p1, std::pair<T,U> p2)
{
    return p1.first < p2.first;
}

template<typename T, typename U> bool compareSecond(std::pair<T,U> p1, std::pair<T,U> p2)
{
    return p1.second < p2.second;
}

template<typename T>
std::vector<T>& translate(std::vector<T>& s, const std::vector<T>& from, const std::vector<T>& to)
{
    typename std::vector<T> fromSorted;
    typename std::vector<T> toSorted;
    typename std::vector< std::pair<T,T> > pairs;

    if (from.size() != to.size()) throw std::logic_error("from and to must be the same size");

    for (typename std::vector<T>::const_iterator f = from.begin(), t = to.begin();f != from.end();++f, ++t)
    {
        pairs.push_back(std::make_pair(*f, *t));
    }
    std::sort(pairs.begin(), pairs.end(), compareFirst<T,T>);
    for (typename std::vector< std::pair<T,T> >::const_iterator i = pairs.begin();i != pairs.end();++i)
    {
        fromSorted.push_back(i->first);
        toSorted.push_back(i->second);
    }

    for (typename std::vector<T>::iterator l = s.begin();l != s.end();++l)
    {
        typename std::vector<T>::const_iterator lb = std::lower_bound(fromSorted.begin(), fromSorted.end(), *l);

        if (lb != fromSorted.end() && *lb == *l)
        {
            //cout << "from " << *l << " to " << to[lb - fromSorted.begin()] << endl;
            *l = toSorted[lb - fromSorted.begin()];
        }
    }

    return s;
}

template<typename T>
std::vector<T> translate_copy(const std::vector<T>& s, const std::vector<T>& from, const std::vector<T>& to)
{
    typename std::vector<T> s_(s);
    translate(s_, from, to);
    return s_;
}

inline std::string& translate(std::string& s, const std::string& from, const std::string& to)
{
    unsigned char trans[256];

    if (from.size() != to.size()) throw std::logic_error("from and to must be the same size");

    if (s.size() < 256)
    {
        for (int i = 0;i < s.size();i++) trans[(unsigned char)s[i]] = s[i];
    }
    else
    {
        for (int i = 0;i < 256;i++) trans[i] = i;
    }

    for (int i = 0;i < from.size();i++) trans[(unsigned char)from[i]] = to[i];
    for (int i = 0;i < s.size();i++) s[i] = trans[(unsigned char)s[i]];

    return s;
}

inline std::string translate_copy(const std::string& s, const std::string& from, const std::string& to)
{
    std::string s_(s);
    translate(s_, from, to);
    return s_;
}

inline std::string toupper(const std::string& s)
{
    std::string S(s);
    for (std::string::iterator C = S.begin();C != S.end();++C) *C = toupper(*C);
    return S;
}

inline std::string tolower(const std::string& S)
{
    std::string s(S);
    for (std::string::iterator c = s.begin();c != s.end();++c) *c = tolower(*c);
    return s;
}

inline float conj(float v) { return v; }
inline double conj(double v) { return v; }

}

#if defined(CPLUSPLUS11)

#include <type_traits>

#elif defined(BOOST)

#include <boost/type_traits.hpp>

namespace std
{
    using namespace boost;
}

#define enable_if enable_if_c

#else

#error "type_traits not available"

#endif

#define ENABLE_IF_CONST(const_type,return_type) \
template <typename _IsConst = const_type > \
typename std::enable_if<std::is_const<_IsConst>::value, return_type >::type

#define ENABLE_IF_NON_CONST(const_type,return_type) \
template <typename _IsConst = const_type > \
typename std::enable_if<!std::is_const<_IsConst>::value, return_type >::type

#define ENABLE_IF_SAME(old_type,new_type,return_type) \
template <typename new_type > \
typename std::enable_if<std::is_same<const old_type, const new_type >::value, return_type >::type

#endif

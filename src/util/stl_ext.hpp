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
#include <iostream>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <cctype>
#include <iterator>
#include <complex>
#include <cstdarg>
#include <cstdio>

#include "memory/memory.h"

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || _MSC_VER >= 1600 || __cplusplus >= 201103l

#include <type_traits>
#include <memory>

#define NON_COPYABLE(name) \
public: \
    name() = delete; \
private:

#define NON_ASSIGNABLE(name) \
public: \
    name& operator=(const name&) = delete; \
private:

#else

#define NON_COPYABLE(name) \
private: \
    name(const name&);

#define NON_ASSIGNABLE(name) \
private: \
    void operator=(const name&);

namespace std
{
    template <bool cond, class return_type = void> struct enable_if {};
    template <class return_type> struct enable_if<true, return_type> { typedef return_type type; };
    template <> struct enable_if<true> { typedef void type; };

    template <class T, class U> struct is_same      { static const bool value = false; };
    template <class T>          struct is_same<T,T> { static const bool value = true; };

    template <class T> struct remove_cv                   {typedef T type; };
    template <class T> struct remove_cv<const T>          {typedef T type; };
    template <class T> struct remove_cv<volatile T>       {typedef T type; };
    template <class T> struct remove_cv<const volatile T> {typedef T type; };

    template <class T> struct make_unsigned__                   { typedef T type; };
    template <>        struct make_unsigned__<signed char>      { typedef unsigned char type; };
    template <>        struct make_unsigned__<signed short>     { typedef unsigned short type; };
    template <>        struct make_unsigned__<signed int>       { typedef unsigned int type; };
    template <>        struct make_unsigned__<signed long>      { typedef unsigned long type; };
    template <>        struct make_unsigned__<signed long long> { typedef unsigned long long type; };

    template <class T> struct make_unsigned                   { typedef typename make_unsigned__<T>::type type; };
    template <class T> struct make_unsigned<const T>          { typedef const typename make_unsigned__<T>::type type; };
    template <class T> struct make_unsigned<volatile T>       { typedef volatile typename make_unsigned__<T>::type type; };
    template <class T> struct make_unsigned<const volatile T> { typedef const volatile typename make_unsigned__<T>::type type; };

    template <class T> struct is_const          { static const bool value = false; };
    template <class T> struct is_const<const T> { static const bool value = true; };

    template <class T, class U = void> struct is_integral { static const bool value = false; };
    template <class T> struct is_integral<T,
        typename enable_if<is_same<unsigned char,typename remove_cv<typename make_unsigned<T>::type>::type>::value>::type>
        { static const bool value = true; };
    template <class T> struct is_integral<T,
        typename enable_if<is_same<unsigned short,typename remove_cv<typename make_unsigned<T>::type>::type>::value>::type>
        { static const bool value = true; };
    template <class T> struct is_integral<T,
        typename enable_if<is_same<unsigned int,typename remove_cv<typename make_unsigned<T>::type>::type>::value>::type>
        { static const bool value = true; };
    template <class T> struct is_integral<T,
        typename enable_if<is_same<unsigned long,typename remove_cv<typename make_unsigned<T>::type>::type>::value>::type>
        { static const bool value = true; };
    template <class T> struct is_integral<T,
        typename enable_if<is_same<unsigned long long,typename remove_cv<typename make_unsigned<T>::type>::type>::value>::type>
        { static const bool value = true; };
    template <class T> struct is_integral<T,
        typename enable_if<is_same<bool,typename remove_cv<T>::type>::value>::type>
        { static const bool value = true; };
    template <class T> struct is_integral<T,
        typename enable_if<is_same<wchar_t,typename remove_cv<T>::type>::value>::type>
        { static const bool value = true; };

    template <class T, class U = void> struct is_floating_point { static const bool value = false; };
    template <class T> struct is_floating_point<T,
        typename enable_if<is_same<float,typename remove_cv<T>::type>::value>::type>
        { static const bool value = true; };
    template <class T> struct is_floating_point<T,
        typename enable_if<is_same<double,typename remove_cv<T>::type>::value>::type>
        { static const bool value = true; };
    template <class T> struct is_floating_point<T,
        typename enable_if<is_same<long double,typename remove_cv<T>::type>::value>::type>
        { static const bool value = true; };

    template <class T, class U = void> struct is_arithmetic { static const bool value = false; };
    template <class T> struct is_arithmetic<T,
        typename enable_if<is_integral<T>::value>::type>
        { static const bool value = true; };
    template <class T> struct is_arithmetic<T,
        typename enable_if<is_floating_point<T>::value>::type>
        { static const bool value = true; };

    template <class T, class U = void> struct is_pod { static const bool value = false; };
    template <class T> struct is_pod<T,
        typename enable_if<is_arithmetic<T>::value>::type>
        { static const bool value = true; };
    template <class T> struct is_pod<complex<T>,
        typename enable_if<is_arithmetic<T>::value>::type>
        { static const bool value = true; };

    template <class T>
    class shared_ptr
    {
        friend void swap(shared_ptr& a, shared_ptr& b)
        {
            a.swap(b);
        }

        friend ostream& operator<<(ostream& os, const shared_ptr& sp)
        {
            os << sp.get();
            return os;
        }

        friend bool operator==(const shared_ptr& a, const shared_ptr& b)
        {
            return a.get() == b.get();
        }

        friend bool operator==(const T* a, const shared_ptr& b)
        {
            return a == b.get();
        }

        friend bool operator==(const shared_ptr& a, const T* b)
        {
            return a.get() == b;
        }

        friend bool operator<(const shared_ptr& a, const shared_ptr& b)
        {
            return a.get() < b.get();
        }

        friend bool operator<(const T* a, const shared_ptr& b)
        {
            return a < b.get();
        }

        friend bool operator<(const shared_ptr& a, const T* b)
        {
            return a.get() < b;
        }

        protected:
            struct ref_ptr
            {
                T* ptr;
                long count;
                ref_ptr(T* ptr)
                : ptr(ptr), count(0) {}
            };
            ref_ptr* ptr;

            void assign_ptr(ref_ptr* ptr)
            {
                this->ptr = ptr;
                if (ptr != NULL) ptr->count++;
            }

            void release_ptr()
            {
                if (ptr == NULL) return;
                ptr->count--;
                if (ptr->count == 0)
                {
                    delete ptr->ptr;
                    delete ptr;
                }
                ptr = NULL;
            }

        public:
            shared_ptr()
            {
                assign_ptr(NULL);
            }

            shared_ptr(T* p)
            {
                if (p == NULL)
                {
                    assign_ptr(NULL);
                }
                else
                {
                    assign_ptr(new ref_ptr(p));
                }
            }

            shared_ptr(const shared_ptr& other)
            {
                assign_ptr(other.ptr);
            }

            ~shared_ptr()
            {
                release_ptr();
            }

            shared_ptr& operator=(const shared_ptr& other)
            {
                if (other.ptr == this->ptr) return *this;
                release_ptr();
                assign_ptr(other.ptr);
                return *this;
            }

            void swap(shared_ptr& other)
            {
                std::swap(ptr, other.ptr);
            }

            void reset()
            {
                release_ptr();
            }

            void reset(T* p)
            {
                release_ptr();
                if (p == NULL)
                {
                    assign_ptr(NULL);
                }
                else
                {
                    assign_ptr(new ref_ptr(p));
                }
            }

            T* get()
            {
                if (ptr == NULL) return NULL;
                return ptr->ptr;
            }

            const T* get() const
            {
                if (ptr == NULL) return NULL;
                return ptr->ptr;
            }

            T& operator*()
            {
                return *get();
            }

            const T& operator*() const
            {
                return *get();
            }

            T* operator->()
            {
                return get();
            }

            const T* operator->() const
            {
                return get();
            }

            long use_count() const
            {
                if (ptr == NULL) return 0;
                return ptr->count;
            }

            bool unique() const
            {
                return use_count() == 1;
            }

            operator bool() const
            {
                return ptr != NULL;
            }
    };

    template <class T>
    class myvector : public vector<T>
    {
        public:
            using typename vector<T>::size_type;
            using typename vector<T>::value_type;
            using typename vector<T>::pointer;
            using typename vector<T>::const_pointer;

            explicit myvector() : vector<T>() {}

            explicit myvector(size_type n, const value_type& val = value_type())
            : vector<T>(n, val) {}

            template <class InputIterator>
            myvector (InputIterator first, InputIterator last)
            : vector<T>(first, last) {}

            myvector (const myvector& x)
            : vector<T>(x) {}

            pointer data()
            {
                if (this->size() == 0)
                {
                    return NULL;
                }
                else
                {
                    return &(*this)[0];
                }
            }

            const_pointer data() const
            {
                if (this->size() == 0)
                {
                    return NULL;
                }
                else
                {
                    return &(*this)[0];
                }
            }
    };
    #define vector myvector

    #define valarray myvalarray
    template <class T>
    class valarray
    {
        friend void swap(valarray& v1, valarray& v2)
        {
            v1.swap(v2);
        }

        friend T* begin(valarray& v)
        {
            return v.ptr;
        }

        friend const T* begin(const valarray& v)
        {
            return v.ptr;
        }

        friend T* end(valarray& v)
        {
            return v.ptr+v.n;
        }

        friend const T* end(const valarray& v)
        {
            return v.ptr+v.n;
        }

        friend valarray operator+(const valarray& v1, const valarray& v2)
        {
            #ifdef _GLIBCXX_DEBUG
            assert(v1.n == v2.n);
            #endif

            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]+v2[i];
            }

            return out;
        }

        friend valarray operator-(const valarray& v1, const valarray& v2)
        {
            #ifdef _GLIBCXX_DEBUG
            assert(v1.n == v2.n);
            #endif

            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]-v2[i];
            }

            return out;
        }

        friend valarray operator*(const valarray& v1, const valarray& v2)
        {
            #ifdef _GLIBCXX_DEBUG
            assert(v1.n == v2.n);
            #endif

            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]*v2[i];
            }

            return out;
        }

        friend valarray operator/(const valarray& v1, const valarray& v2)
        {
            #ifdef _GLIBCXX_DEBUG
            assert(v1.n == v2.n);
            #endif

            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]/v2[i];
            }

            return out;
        }

        friend valarray operator%(const valarray& v1, const valarray& v2)
        {
            #ifdef _GLIBCXX_DEBUG
            assert(v1.n == v2.n);
            #endif

            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]%v2[i];
            }

            return out;
        }

        friend valarray operator^(const valarray& v1, const valarray& v2)
        {
            #ifdef _GLIBCXX_DEBUG
            assert(v1.n == v2.n);
            #endif

            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]^v2[i];
            }

            return out;
        }

        friend valarray operator&(const valarray& v1, const valarray& v2)
        {
            #ifdef _GLIBCXX_DEBUG
            assert(v1.n == v2.n);
            #endif

            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]&v2[i];
            }

            return out;
        }

        friend valarray operator|(const valarray& v1, const valarray& v2)
        {
            #ifdef _GLIBCXX_DEBUG
            assert(v1.n == v2.n);
            #endif

            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]|v2[i];
            }

            return out;
        }

        friend valarray operator<<(const valarray& v1, const valarray& v2)
        {
            #ifdef _GLIBCXX_DEBUG
            assert(v1.n == v2.n);
            #endif

            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]<<v2[i];
            }

            return out;
        }

        friend valarray operator>>(const valarray& v1, const valarray& v2)
        {
            #ifdef _GLIBCXX_DEBUG
            assert(v1.n == v2.n);
            #endif

            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]>>v2[i];
            }

            return out;
        }

        friend valarray operator+(const T& v1, const valarray& v2)
        {
            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1+v2[i];
            }

            return out;
        }

        friend valarray operator-(const T& v1, const valarray& v2)
        {
            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1-v2[i];
            }

            return out;
        }

        friend valarray operator*(const T& v1, const valarray& v2)
        {
            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1*v2[i];
            }

            return out;
        }

        friend valarray operator/(const T& v1, const valarray& v2)
        {
            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1/v2[i];
            }

            return out;
        }

        friend valarray operator%(const T& v1, const valarray& v2)
        {
            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1%v2[i];
            }

            return out;
        }

        friend valarray operator^(const T& v1, const valarray& v2)
        {
            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1^v2[i];
            }

            return out;
        }

        friend valarray operator&(const T& v1, const valarray& v2)
        {
            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1&v2[i];
            }

            return out;
        }

        friend valarray operator|(const T& v1, const valarray& v2)
        {
            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1|v2[i];
            }

            return out;
        }

        friend valarray operator<<(const T& v1, const valarray& v2)
        {
            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1<<v2[i];
            }

            return out;
        }

        friend valarray operator>>(const T& v1, const valarray& v2)
        {
            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1>>v2[i];
            }

            return out;
        }

        friend valarray operator+(const valarray& v1, const T& v2)
        {
            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]+v2;
            }

            return out;
        }

        friend valarray operator-(const valarray& v1, const T& v2)
        {
            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]-v2;
            }

            return out;
        }

        friend valarray operator*(const valarray& v1, const T& v2)
        {
            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]*v2;
            }

            return out;
        }

        friend valarray operator/(const valarray& v1, const T& v2)
        {
            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]/v2;
            }

            return out;
        }

        friend valarray operator%(const valarray& v1, const T& v2)
        {
            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]%v2;
            }

            return out;
        }

        friend valarray operator^(const valarray& v1, const T& v2)
        {
            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]^v2;
            }

            return out;
        }

        friend valarray operator&(const valarray& v1, const T& v2)
        {
            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]&v2;
            }

            return out;
        }

        friend valarray operator|(const valarray& v1, const T& v2)
        {
            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]|v2;
            }

            return out;
        }

        friend valarray operator<<(const valarray& v1, const T& v2)
        {
            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]<<v2;
            }

            return out;
        }

        friend valarray operator>>(const valarray& v1, const T& v2)
        {
            valarray out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]>>v2;
            }

            return out;
        }

        friend valarray<bool> operator&&(const valarray& v1, const valarray& v2)
        {
            #ifdef _GLIBCXX_DEBUG
            assert(v1.n == v2.n);
            #endif

            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]&&v2[i];
            }

            return out;
        }

        friend valarray<bool> operator||(const valarray& v1, const valarray& v2)
        {
            #ifdef _GLIBCXX_DEBUG
            assert(v1.n == v2.n);
            #endif

            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]||v2[i];
            }

            return out;
        }

        friend valarray<bool> operator==(const valarray& v1, const valarray& v2)
        {
            #ifdef _GLIBCXX_DEBUG
            assert(v1.n == v2.n);
            #endif

            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]==v2[i];
            }

            return out;
        }

        friend valarray<bool> operator!=(const valarray& v1, const valarray& v2)
        {
            #ifdef _GLIBCXX_DEBUG
            assert(v1.n == v2.n);
            #endif

            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]!=v2[i];
            }

            return out;
        }

        friend valarray<bool> operator<(const valarray& v1, const valarray& v2)
        {
            #ifdef _GLIBCXX_DEBUG
            assert(v1.n == v2.n);
            #endif

            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]<v2[i];
            }

            return out;
        }

        friend valarray<bool> operator>(const valarray& v1, const valarray& v2)
        {
            #ifdef _GLIBCXX_DEBUG
            assert(v1.n == v2.n);
            #endif

            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]>v2[i];
            }

            return out;
        }

        friend valarray<bool> operator<=(const valarray& v1, const valarray& v2)
        {
            #ifdef _GLIBCXX_DEBUG
            assert(v1.n == v2.n);
            #endif

            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]<=v2[i];
            }

            return out;
        }

        friend valarray<bool> operator>=(const valarray& v1, const valarray& v2)
        {
            #ifdef _GLIBCXX_DEBUG
            assert(v1.n == v2.n);
            #endif

            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]>=v2[i];
            }

            return out;
        }

        friend valarray<bool> operator&&(const T& v1, const valarray& v2)
        {
            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1&&v2[i];
            }

            return out;
        }

        friend valarray<bool> operator||(const T& v1, const valarray& v2)
        {
            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1||v2[i];
            }

            return out;
        }

        friend valarray<bool> operator==(const T& v1, const valarray& v2)
        {
            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1==v2[i];
            }

            return out;
        }

        friend valarray<bool> operator!=(const T& v1, const valarray& v2)
        {
            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1!=v2[i];
            }

            return out;
        }

        friend valarray<bool> operator<(const T& v1, const valarray& v2)
        {
            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1<v2[i];
            }

            return out;
        }

        friend valarray<bool> operator>(const T& v1, const valarray& v2)
        {
            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1>v2[i];
            }

            return out;
        }

        friend valarray<bool> operator<=(const T& v1, const valarray& v2)
        {
            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1<=v2[i];
            }

            return out;
        }

        friend valarray<bool> operator>=(const T& v1, const valarray& v2)
        {
            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1>=v2[i];
            }

            return out;
        }

        friend valarray<bool> operator&&(const valarray& v1, const T& v2)
        {
            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]&&v2;
            }

            return out;
        }

        friend valarray<bool> operator||(const valarray& v1, const T& v2)
        {
            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]||v2;
            }

            return out;
        }

        friend valarray<bool> operator==(const valarray& v1, const T& v2)
        {
            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]==v2;
            }

            return out;
        }

        friend valarray<bool> operator!=(const valarray& v1, const T& v2)
        {
            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]!=v2;
            }

            return out;
        }

        friend valarray<bool> operator<(const valarray& v1, const T& v2)
        {
            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]<v2;
            }

            return out;
        }

        friend valarray<bool> operator>(const valarray& v1, const T& v2)
        {
            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]>v2;
            }

            return out;
        }

        friend valarray<bool> operator<=(const valarray& v1, const T& v2)
        {
            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]<=v2;
            }

            return out;
        }

        friend valarray<bool> operator>=(const valarray& v1, const T& v2)
        {
            valarray<bool> out(v1.n, false);

            for (size_t i = 0;i < v1.n;i++)
            {
                out[i] = v1[i]>=v2;
            }

            return out;
        }

        protected:
            T* ptr;
            size_t n;

            valarray(size_t n, bool nozero)
            : ptr(NULL), n(n)
            {
                if (n > 0) ptr = new T[n];
            }

        public:
            valarray()
            : ptr(NULL), n(0) {}

            explicit valarray (size_t n)
            : ptr(NULL), n(n)
            {
                if (n > 0)
                {
                    ptr = new T[n];
                    fill(ptr, ptr+n, T());
                }
            }

            valarray (const T& val, size_t n)
            : ptr(NULL), n(n)
            {
                if (n > 0)
                {
                    ptr = new T[n];
                    fill(ptr, ptr+n, val);
                }
            }

            valarray (const T* p, size_t n)
            : ptr(NULL), n(n)
            {
                if (n > 0)
                {
                    ptr = new T[n];
                    copy(p, p+n, ptr);
                }
            }

            valarray (const valarray& x)
            : ptr(NULL), n(x.n)
            {
                if (n > 0)
                {
                    ptr = new T[n];
                    copy(x.ptr, x.ptr+n, ptr);
                }
            }

            ~valarray()
            {
                if (ptr != NULL) delete[] ptr;
            }

            valarray apply(T func(T)) const
            {
                valarray out(n, false);

                for (size_t i = 0;i < n;i++)
                {
                    out[i] = func(ptr[i]);
                }

                return out;
            }

            valarray apply(T func(const T&)) const
            {
                valarray out(n, false);

                for (size_t i = 0;i < n;i++)
                {
                    out[i] = func(ptr[i]);
                }

                return out;
            }

            valarray cshift(int s) const
            {
                valarray out(n, false);

                if (s >= n || s <= -n) s = s%n;

                if (s < 0)
                {
                    copy(ptr    , ptr+n+s, out.ptr-s);
                    copy(ptr+n+s, ptr+n  , out.ptr  );
                }
                else
                {
                    copy(ptr+s, ptr+n, out.ptr    );
                    copy(ptr  , ptr+s, out.ptr+n-s);
                }

                return out;
            }

            T max() const
            {
                if (n == 0) return T();

                T mx = ptr[0];

                for (size_t i = 1;i < n;i++)
                {
                    mx = max(mx, ptr[i]);
                }

                return mx;
            }

            T min() const
            {
                if (n == 0) return T();

                T mn = ptr[0];

                for (size_t i = 1;i < n;i++)
                {
                    mn = min(mn, ptr[i]);
                }

                return mn;
            }

            valarray& operator=(const valarray& other)
            {
                if (n != other.n)
                {
                    if (ptr != NULL) delete[] ptr;
                    if (other.n > 0) ptr = new T[other.n];
                }

                n = other.n;
                copy(other.ptr, other.ptr+n, ptr);

                return *this;
            }

            T& operator[](size_t i)
            {
                #ifdef _GLIBCXX_DEBUG
                assert(i < n);
                #endif
                return ptr[i];
            }

            T operator[](size_t i) const
            {
                #ifdef _GLIBCXX_DEBUG
                assert(i < n);
                #endif
                return ptr[i];
            }

            void resize(size_t sz, T c = T())
            {
                if (n > 0) delete[] ptr;
                n = sz;
                if (n > 0) ptr = new T[n];
                fill(ptr, ptr+n, c);
            }

            valarray shift(int s) const
            {
                valarray out(n, false);

                if (s >= n || s <= -n) s = s%n;

                if (s < 0)
                {
                    copy(ptr, ptr+n+s, out.ptr-s);
                    fill(out.ptr, out.ptr-s, T());
                }
                else
                {
                    copy(ptr+s, ptr+n, out.ptr);
                    fill(out.ptr+n-s, out.ptr+n, T());
                }

                return out;
            }

            size_t size() const
            {
                return n;
            }

            T sum() const
            {
                T s = T();

                for (size_t i = 0;i < n;i++)
                {
                    s += ptr[i];
                }

                return s;
            }

            void swap(valarray& other)
            {
                swap(n, other.n);
                swap(ptr, other.ptr);
            }

            valarray operator+() const
            {
                return valarray(*this);
            }

            valarray operator-() const
            {
                valarray out(n, false);

                for (size_t i = 0;i < n;i++)
                {
                    out[i] = -ptr[i];
                }

                return out;
            }

            valarray operator~() const
            {
                valarray out(n, false);

                for (size_t i = 0;i < n;i++)
                {
                    out[i] = ~ptr[i];
                }

                return out;
            }

            valarray<bool> operator!() const
            {
                valarray<bool> out(n, false);

                for (size_t i = 0;i < n;i++)
                {
                    out[i] = !ptr[i];
                }

                return out;
            }

            valarray& operator*= (const valarray& rhs)
            {
                #ifdef _GLIBCXX_DEBUG
                assert(n == rhs.n);
                #endif

                for (size_t i = 0;i < n;i++)
                {
                    ptr[i] *= rhs[i];
                }

                return *this;
            }

            valarray& operator/= (const valarray& rhs)
            {
                #ifdef _GLIBCXX_DEBUG
                assert(n == rhs.n);
                #endif

                for (size_t i = 0;i < n;i++)
                {
                    ptr[i] /= rhs[i];
                }

                return *this;
            }

            valarray& operator%= (const valarray& rhs)
            {
                #ifdef _GLIBCXX_DEBUG
                assert(n == rhs.n);
                #endif

                for (size_t i = 0;i < n;i++)
                {
                    ptr[i] %= rhs[i];
                }

                return *this;
            }

            valarray& operator+= (const valarray& rhs)
            {
                #ifdef _GLIBCXX_DEBUG
                assert(n == rhs.n);
                #endif

                for (size_t i = 0;i < n;i++)
                {
                    ptr[i] += rhs[i];
                }

                return *this;
            }

            valarray& operator-= (const valarray& rhs)
            {
                #ifdef _GLIBCXX_DEBUG
                assert(n == rhs.n);
                #endif

                for (size_t i = 0;i < n;i++)
                {
                    ptr[i] -= rhs[i];
                }

                return *this;
            }

            valarray& operator^= (const valarray& rhs)
            {
                #ifdef _GLIBCXX_DEBUG
                assert(n == rhs.n);
                #endif

                for (size_t i = 0;i < n;i++)
                {
                    ptr[i] ^= rhs[i];
                }

                return *this;
            }

            valarray& operator&= (const valarray& rhs)
            {
                #ifdef _GLIBCXX_DEBUG
                assert(n == rhs.n);
                #endif

                for (size_t i = 0;i < n;i++)
                {
                    ptr[i] &= rhs[i];
                }

                return *this;
            }

            valarray& operator|= (const valarray& rhs)
            {
                #ifdef _GLIBCXX_DEBUG
                assert(n == rhs.n);
                #endif

                for (size_t i = 0;i < n;i++)
                {
                    ptr[i] |= rhs[i];
                }

                return *this;
            }

            valarray& operator<<= (const valarray& rhs)
            {
                #ifdef _GLIBCXX_DEBUG
                assert(n == rhs.n);
                #endif

                for (size_t i = 0;i < n;i++)
                {
                    ptr[i] <<= rhs[i];
                }

                return *this;
            }

            valarray& operator>>= (const valarray& rhs)
            {
                #ifdef _GLIBCXX_DEBUG
                assert(n == rhs.n);
                #endif

                for (size_t i = 0;i < n;i++)
                {
                    ptr[i] >>= rhs[i];
                }

                return *this;
            }

            valarray& operator*= (const T& val)
            {
                for (size_t i = 0;i < n;i++)
                {
                    ptr[i] *= val;
                }

                return *this;
            }

            valarray& operator/= (const T& val)
            {
                for (size_t i = 0;i < n;i++)
                {
                    ptr[i] /= val;
                }

                return *this;
            }

            valarray& operator%= (const T& val)
            {
                for (size_t i = 0;i < n;i++)
                {
                    ptr[i] %= val;
                }

                return *this;
            }

            valarray& operator+= (const T& val)
            {
                for (size_t i = 0;i < n;i++)
                {
                    ptr[i] += val;
                }

                return *this;
            }

            valarray& operator-= (const T& val)
            {
                for (size_t i = 0;i < n;i++)
                {
                    ptr[i] -= val;
                }

                return *this;
            }

            valarray& operator^= (const T& val)
            {
                for (size_t i = 0;i < n;i++)
                {
                    ptr[i] ^= val;
                }

                return *this;
            }

            valarray& operator&= (const T& val)
            {
                for (size_t i = 0;i < n;i++)
                {
                    ptr[i] &= val;
                }

                return *this;
            }

            valarray& operator|= (const T& val)
            {
                for (size_t i = 0;i < n;i++)
                {
                    ptr[i] |= val;
                }

                return *this;
            }

            valarray& operator<<= (const T& val)
            {
                for (size_t i = 0;i < n;i++)
                {
                    ptr[i] <<= val;
                }

                return *this;
            }

            valarray& operator>>= (const T& val)
            {
                for (size_t i = 0;i < n;i++)
                {
                    ptr[i] >>= val;
                }

                return *this;
            }
    };

    template<typename I1, typename I2, typename Pred>
    I2 copy_if(I1 begin, I1 end, I2 result, Pred pred)
    {
        return remove_copy_if(begin, end, result, not1(pred));
    }
}

#endif

#define ENABLE_IF_CONST(const_type,return_type) \
template <typename _IsConst = const_type > \
typename std::enable_if<std::is_const<_IsConst>::value, return_type >::type

#define ENABLE_IF_NON_CONST(const_type,return_type) \
template <typename _IsConst = const_type > \
typename std::enable_if<!std::is_const<_IsConst>::value, return_type >::type

//#define ENABLE_IF_SAME(old_type,new_type,return_type) \
//template <typename new_type > \
//typename std::enable_if<std::is_base_of<const old_type, const new_type >::value, return_type >::type

#define ENABLE_IF_SAME(old_type,new_type,return_type) \
template <typename new_type > return_type

namespace std
{

template <class T, class U>
struct if_exists
{
    typedef U type;
};

template <class T>
class vararray
{
    friend void swap(vararray& v1, vararray& v2)
    {
        v1.swap(v2);
    }

    friend bool operator==(const vararray& lhs, const vararray& rhs)
    {
        return equal(lhs.begin(), lhs.end(), rhs.begin());
    }

    friend  bool operator!=(const vararray& lhs, const vararray& rhs)
    {
        return !equal(lhs.begin(), lhs.end(), rhs.begin());
    }

    friend  bool operator<(const vararray& lhs, const vararray& rhs)
    {
        return lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
    }

    friend  bool operator<=(const vararray& lhs, const vararray& rhs)
    {
        return !lexicographical_compare(rhs.begin(), rhs.end(), lhs.begin(), lhs.end());
    }

    friend  bool operator>(const vararray& lhs, const vararray& rhs)
    {
        return lexicographical_compare(rhs.begin(), rhs.end(), lhs.begin(), lhs.end());
    }

    friend  bool operator>=(const vararray& lhs, const vararray& rhs)
    {
        return !lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
    }

    public:
        typedef T value_type;
        typedef T& reference;
        typedef const T& const_reference;
        typedef T* pointer;
        typedef const T* const_pointer;
        typedef T* iterator;
        typedef const T* const_iterator;
        typedef std::reverse_iterator<iterator> reverse_iterator;
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
        typedef ptrdiff_t difference_type;
        typedef size_t size_type;

    protected:
        pointer ptr;
        size_type n;

    public:
        vararray()
        : ptr(NULL), n(0) {}

        explicit vararray(size_type n)
        : ptr(NULL), n(n)
        {
            if (n > 0)
            {
                ptr = SAFE_MALLOC(T,n);
                if (!is_pod<T>::value) fill(ptr, ptr+n, value_type());
            }
        }

        vararray(size_type n, const value_type& val)
        : ptr(NULL), n(n)
        {
            if (n > 0)
            {
                ptr = SAFE_MALLOC(T,n);
                fill(ptr, ptr+n, val);
            }
        }

        template <class RAIterator>
        vararray(RAIterator first, RAIterator last)
        : ptr(NULL), n(0)
        {
            #ifdef _GLIBCXX_DEBUG
            assert(first < last);
            #endif
            n = last-first;
            if (n > 0) ptr = SAFE_MALLOC(T,n);
            for (size_t i = 0;first != last;i++,++first)
            {
                ptr[i] = *first;
            }
        }

        vararray(const vararray& x)
        : ptr(NULL), n(x.n)
        {
            if (n > 0)
            {
                ptr = SAFE_MALLOC(T,n);
                copy(x.ptr, x.ptr+n, ptr);
            }
        }

        ~vararray()
        {
            if (ptr != NULL) FREE(ptr);
        }

        vararray& operator=(const vararray& other)
        {
            if (ptr != NULL) FREE(ptr);
            n = other.n;
            if (n > 0)
            {
                ptr = SAFE_MALLOC(T,n);
                copy(other.ptr, other.ptr+n, ptr);
            }
            return *this;
        }

        iterator begin() { return ptr; }

        const_iterator begin() const { return ptr; }

        iterator end()  { return ptr+n; }

        const_iterator end() const  { return ptr+n; }

        reverse_iterator rbegin()  { return reverse_iterator(ptr+n); }

        const_reverse_iterator rbegin() const  { return const_reverse_iterator(ptr+n); }

        reverse_iterator rend() { return reverse_iterator(ptr); }

        const_reverse_iterator rend() const { return const_reverse_iterator(ptr); }

        size_type size() const { return n; }

        size_type max_size() const { return n; }

        void resize(size_type n)
        {
            ptr = SAFE_REALLOC(T,ptr,n);
            if (n > this->n && !is_pod<T>::value) fill(ptr+this->n, ptr+n, value_type());
            this->n = n;
        }

        void resize(size_type n, value_type val)
        {
            ptr = SAFE_REALLOC(T,ptr,n);
            if (n > this->n) fill(ptr+this->n, ptr+n, val);
            this->n = n;
        }

        bool empty() const { return n > 0; }

        reference operator[](size_type i) { return at(i); }

        const_reference operator[](size_type i) const { return at(i); }

        reference at(size_type i)
        {
            #ifdef _GLIBCXX_DEBUG
            assert(i < n);
            #endif
            return ptr[i];
        }

        const_reference at(size_type i) const
        {
            #ifdef _GLIBCXX_DEBUG
            assert(i < n);
            #endif
            return ptr[i];
        }

        reference front() { return at(0); }

        const_reference front() const { return at(0); }

        reference back()
        {
            #ifdef _GLIBCXX_DEBUG
            assert(n > 0);
            #endif
            return at(n-1);
        }

        const_reference back() const
        {
            #ifdef _GLIBCXX_DEBUG
            assert(n > 0);
            #endif
            return at(n-1);
        }

        pointer data() { return ptr; }

        const_pointer data() const { return ptr; }

        template <class RAIterator>
        void assign(RAIterator first, RAIterator last)
        {
            #ifdef _GLIBCXX_DEBUG
            assert(first < last);
            #endif
            this->n = last-first;
            ptr = SAFE_REALLOC(T,ptr,n);
            for (size_type i = 0;first != last;i++,++first)
            {
                ptr[i] = *first;
            }
        }

        void assign(size_type n, const value_type& val)
        {
            this->n = n;
            ptr = SAFE_REALLOC(T,ptr,n);
            fill(ptr, ptr+n, val);
        }

        void swap(vararray& other)
        {
            std::swap(ptr, other.ptr);
            std::swap(n, other.n);
        }

        void clear()
        {
            if (ptr != NULL)
            {
                FREE(ptr);
                ptr = NULL;
            }
            n = 0;
        }
};

template<typename T>
class global_ptr : public shared_ptr<T*>
{
    public:
        global_ptr(const global_ptr& other) : shared_ptr<T*>(other) {}

        global_ptr() : shared_ptr<T*>(new T*(NULL)) {}

        global_ptr(T* ptr) : shared_ptr<T*>(new T*(ptr)) {}

        ~global_ptr()
        {
            if (this->unique() && *shared_ptr<T*>::get())
            {
                delete *shared_ptr<T*>::get();
            }
        }

        global_ptr& operator=(const global_ptr& other)
        {
            if (this->unique() && *shared_ptr<T*>::get())
            {
                delete *shared_ptr<T*>::get();
            }
            shared_ptr<T*>::operator=(other);
            return *this;
        }

        void reset()
        {
            if (*shared_ptr<T*>::get()) delete *shared_ptr<T*>::get();
            *shared_ptr<T*>::get() = NULL;
        }

        void reset(T* p)
        {
            if (p == *shared_ptr<T*>::get()) return;
            if (*shared_ptr<T*>::get()) delete *shared_ptr<T*>::get();
            *shared_ptr<T*>::get() = p;
        }

        T* get()
        {
            return *shared_ptr<T*>::get();
        }

        const T* get() const
        {
            return *shared_ptr<T*>::get();
        }

        T& operator*()
        {
            return **shared_ptr<T*>::get();
        }

        const T& operator*() const
        {
            return **shared_ptr<T*>::get();
        }

        T* operator->()
        {
            return *shared_ptr<T*>::get();
        }

        const T* operator->() const
        {
            return *shared_ptr<T*>::get();
        }

        operator bool() const
        {
            return *shared_ptr<T*>::get() != NULL;
        }
};

inline std::string strprintf(const char* fmt, ...)
{
    va_list list;

    va_start(list, fmt);
    char fake[1];
    int n = vsnprintf(fake, 1, fmt, list);
    va_end(list);

    std::vector<char> s(n+1);
    va_start(list, fmt);
    vsnprintf(s.data(), n+1, fmt, list);
    va_end(list);

    return std::string(s.begin(), s.end());
}

template<typename T> std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    os << "[";
    if (!v.empty()) os << v[0];
    for (int i = 1;i < v.size();i++) os << ", " << v[i];
    os << "]";
    return os;
}

template<typename T> std::vector<T> vec(const T& a)
{
    std::vector<T> v;
    v.push_back(a);
    return v;
}

template<typename T> std::vector<T> vec(const T& a, const T& b)
{
    std::vector<T> v;
    v.push_back(a);
    v.push_back(b);
    return v;
}

template<typename T> std::vector<T> vec(const T& a, const T& b, const T& c)
{
    std::vector<T> v;
    v.push_back(a);
    v.push_back(b);
    v.push_back(c);
    return v;
}

template<typename T> std::vector<T> vec(const T& a, const T& b, const T& c, const T& d)
{
    std::vector<T> v;
    v.push_back(a);
    v.push_back(b);
    v.push_back(c);
    v.push_back(d);
    return v;
}

template<typename T> std::vector<T> vec(const T& a, const T& b, const T& c, const T& d, const T& e)
{
    std::vector<T> v;
    v.push_back(a);
    v.push_back(b);
    v.push_back(c);
    v.push_back(d);
    v.push_back(e);
    return v;
}

template<typename T> std::vector<T> vec(const T& a, const T& b, const T& c, const T& d, const T& e,
                                        const T& f)
{
    std::vector<T> v;
    v.push_back(a);
    v.push_back(b);
    v.push_back(c);
    v.push_back(d);
    v.push_back(e);
    v.push_back(f);
    return v;
}

template<typename T> std::vector<T> vec(const T& a, const T& b, const T& c, const T& d, const T& e,
                                        const T& f, const T& g)
{
    std::vector<T> v;
    v.push_back(a);
    v.push_back(b);
    v.push_back(c);
    v.push_back(d);
    v.push_back(e);
    v.push_back(f);
    v.push_back(g);
    return v;
}

template<typename T> std::vector<T> vec(const T& a, const T& b, const T& c, const T& d, const T& e,
                                        const T& f, const T& g, const T& h)
{
    std::vector<T> v;
    v.push_back(a);
    v.push_back(b);
    v.push_back(c);
    v.push_back(d);
    v.push_back(e);
    v.push_back(f);
    v.push_back(g);
    v.push_back(h);
    return v;
}

template<typename T> std::vector<T> vec(const T& a, const T& b, const T& c, const T& d, const T& e,
                                        const T& f, const T& g, const T& h, const T& i)
{
    std::vector<T> v;
    v.push_back(a);
    v.push_back(b);
    v.push_back(c);
    v.push_back(d);
    v.push_back(e);
    v.push_back(f);
    v.push_back(g);
    v.push_back(h);
    v.push_back(i);
    return v;
}

template<typename T> std::vector<T> vec(const T& a, const T& b, const T& c, const T& d, const T& e,
                                        const T& f, const T& g, const T& h, const T& i, const T& j)
{
    std::vector<T> v;
    v.push_back(a);
    v.push_back(b);
    v.push_back(c);
    v.push_back(d);
    v.push_back(e);
    v.push_back(f);
    v.push_back(g);
    v.push_back(h);
    v.push_back(i);
    v.push_back(j);
    return v;
}
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

template<class Pred1, class Pred2>
class binary_or
{
    protected:
        Pred1 p1;
        Pred2 p2;

    public:
        binary_or(Pred1 p1, Pred2 p2) : p1(p1), p2(p2) {}

        template <typename T>
        bool operator()(const T& t)
        {
            return p1(t)||p2(t);
        }
};

template<class Pred1, class Pred2>
class binary_and
{
    protected:
        Pred1 p1;
        Pred2 p2;

    public:
        binary_and(Pred1 p1, Pred2 p2) : p1(p1), p2(p2) {}

        template <typename T>
        bool operator()(const T& t)
        {
            return p1(t)&&p2(t);
        }
};

template<class Pred1, class Pred2>
binary_or<Pred1,Pred2> or1(Pred1 p1, Pred2 p2)
{
    return binary_or<Pred1,Pred2>(p1,p2);
}

template<class Pred1, class Pred2>
binary_and<Pred1,Pred2> and1(Pred1 p1, Pred2 p2)
{
    return binary_and<Pred1,Pred2>(p1,p2);
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

template<typename T, typename U, class Functor> std::vector<U> apply(std::vector<T>& v, Functor f)
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
    typename std::vector<T> v2(v);
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

template<typename T> typename enable_if<!is_arithmetic<T>::value,T>::type sum(const std::vector<T>& v)
{
    T s;
    for (int i = 0;i < v.size();i++) s += v[i];
    return s;
}

template<typename T> typename enable_if<is_arithmetic<T>::value,T>::type sum(const std::vector<T>& v)
{
    T s = (T)0;
    for (int i = 0;i < v.size();i++) s += v[i];
    return s;
}

template<typename T, typename U> bool contains(const T& v, const U& e)
{
    return find(v.begin(), v.end(), e) != v.end();
}

template<typename T> T& uniq(T& v)
{
    typename T::iterator i1;

    std::sort(v.begin(), v.end());
    i1 = std::unique(v.begin(), v.end());
    v.resize(i1-v.begin());

    return v;
}

template<typename T> T uniq_copy(const T& v)
{
    T v2(v);
    typename T::iterator i1;

    std::sort(v2.begin(), v2.end());
    i1 = std::unique(v2.begin(), v2.end());
    v2.resize(i1-v2.begin());

    return v2;
}

template<typename T> T intersection(const T& v1, const T& v2)
{
    T v;
    T v3(v1);
    T v4(v2);
    typename T::iterator end;

    v.resize(v1.size()+v2.size());

    sort(v3.begin(),v3.end());
    sort(v4.begin(),v4.end());

    end = std::set_intersection(v3.begin(), v3.end(), v4.begin(), v4.end(), v.begin());
    v.resize(end-v.begin());

    return v;
}

/*
 * Exclude elements of v2 from v1
 */
template<typename T> T& exclude(T& v1, const T& v2)
{
    T v3(v2);
    typename T::iterator i1, i2, i3;

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

/*
 * Return elements from v1 that are not also in v2
 */
template<typename T> T exclude_copy(const T& v1, const T& v2)
{
    T v3(v1);
    return exclude(v3, v2);
}

template<typename T, typename U> T& mask(T& v, const U& mask)
{
    typename T::iterator i1;
    typename T::const_iterator i2;
    typename U::const_iterator i3;

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

template<typename T, typename U> T mask_copy(const T& v, const U& mask)
{
    T v2(v);
    typename T::iterator i1;
    typename T::const_iterator i2;
    typename U::const_iterator i3;

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

inline float real(float v) { return v; }
inline double real(double v) { return v; }

inline float imag(float v) { return 0.0; }
inline double imag(double v) { return 0.0; }

template <typename T>
struct real_type
{
    typedef T type;
};

template <typename T>
struct real_type<complex<T> >
{
    typedef T type;
};

template <typename T>
struct complex_type
{
    typedef complex<T> type;
};

template <typename T>
struct complex_type<complex<T> >
{
    typedef complex<T> type;
};

template <class T, class U>
struct doublet
{
    T first;
    U second;

    doublet(T first, U second) : first(first), second(second) {}

    friend void swap(doublet<T,U> first, doublet<T,U> second)
    {
        swap(first.first, second.first);
        swap(first.second, second.second);
    }

    template <typename T_, typename U_>
    doublet(doublet<T_,U_> other)
    : first(other.first), second(other.second) {}

    template <typename T_, typename U_>
    doublet<T,U>& operator=(const doublet<T_,U_>& other)
    {
        first = other.first;
        second = other.second;
        return *this;
    }

    doublet<T,U>& operator=(const doublet<T,U>& other)
    {
        first = other.first;
        second = other.second;
        return *this;
    }

    bool operator==(const doublet<T,U>& other) const
    {
        return first == other.first;
    }

    bool operator!=(const doublet<T,U>& other) const
    {
        return first != other.first;
    }

    bool operator<(const doublet<T,U>& other) const
    {
        return first < other.first;
    }

    bool operator>(const doublet<T,U>& other) const
    {
        return first > other.first;
    }

    bool operator<=(const doublet<T,U>& other) const
    {
        return first <= other.first;
    }

    bool operator>=(const doublet<T,U>& other) const
    {
        return first >= other.first;
    }
};

template <class T, class U>
class coiterator : public iterator<random_access_iterator_tag,
                                   doublet<typename iterator_traits<T>::value_type,
                                           typename iterator_traits<U>::value_type>,
                                   ptrdiff_t,
                                   doublet<typename iterator_traits<T>::pointer,
                                           typename iterator_traits<U>::pointer>,
                                   doublet<typename iterator_traits<T>::reference,
                                           typename iterator_traits<U>::reference> >
{
    T it_T;
    U it_U;

    public:
        coiterator(const T& it_T, const U& it_U) : it_T(it_T), it_U(it_U) {}

        bool operator==(const coiterator<T,U>& other) const
        {
            return it_T == other.it_T;
        }

        bool operator!=(const coiterator<T,U>& other) const
        {
            return it_T != other.it_T;
        }

        bool operator<(const coiterator<T,U>& other) const
        {
            return it_T < other.it_T;
        }

        bool operator>(const coiterator<T,U>& other) const
        {
            return it_T > other.it_T;
        }

        bool operator<=(const coiterator<T,U>& other) const
        {
            return it_T <= other.it_T;
        }

        bool operator>=(const coiterator<T,U>& other) const
        {
            return it_T >= other.it_T;
        }

        typename coiterator<T,U>::reference operator*()
        {
            return typename coiterator<T,U>::reference(*it_T,*it_U);
        }

        typename coiterator<T,U>::reference operator[](ptrdiff_t n)
        {
            return typename coiterator<T,U>::reference(it_T[n],it_U[n]);
        }

        coiterator<T,U>& operator++()
        {
            ++it_T;
            ++it_U;
            return *this;
        }

        coiterator<T,U>& operator--()
        {
            --it_T;
            --it_U;
            return *this;
        }

        coiterator<T,U> operator++(int x)
        {
            return coiterator<T,U>(it_T++, it_U++);
        }

        coiterator<T,U> operator--(int x)
        {
            return coiterator<T,U>(it_T--, it_U--);
        }

        coiterator<T,U>& operator+=(ptrdiff_t n)
        {
            it_T += n;
            it_U += n;
            return *this;
        }

        coiterator<T,U>& operator-=(ptrdiff_t n)
        {
            it_T -= n;
            it_U -= n;
            return *this;
        }

        coiterator<T,U> operator+(ptrdiff_t n) const
        {
            return coiterator<T,U>(it_T+n, it_U+n);
        }

        friend coiterator<T,U> operator+(ptrdiff_t n, const coiterator<T,U>& other)
        {
            return coiterator<T,U>(other.it_T+n, other.it_U+n);
        }

        coiterator<T,U> operator-(ptrdiff_t n) const
        {
            return coiterator<T,U>(it_T-n, it_U-n);
        }

        ptrdiff_t operator-(const coiterator<T,U>& other) const
        {
            return it_T-other.it_T;
        }
};

template <class key_iterator, class val_iterator, class Comparator>
class cocomparator
{
    typedef typename coiterator<key_iterator,val_iterator>::value_type val;

    Comparator comp;

    public:
        cocomparator(Comparator comp) : comp(comp) {}

        bool operator()(const val& r1, const val& r2) const
        {
            return comp(r1.first, r2.first);
        }
};

template <class key_iterator, class val_iterator>
void cosort(key_iterator keys_begin, key_iterator keys_end,
            val_iterator vals_begin, val_iterator vals_end)
{
    coiterator<key_iterator,val_iterator> begin = coiterator<key_iterator,val_iterator>(keys_begin, vals_begin);
    coiterator<key_iterator,val_iterator> end = coiterator<key_iterator,val_iterator>(keys_end, vals_end);
    sort(begin, end);
}

template <class key_iterator, class val_iterator, class Comparator>
void cosort(key_iterator keys_begin, key_iterator keys_end,
            val_iterator vals_begin, val_iterator vals_end,
            Comparator comp)
{
    coiterator<key_iterator,val_iterator> begin = coiterator<key_iterator,val_iterator>(keys_begin, vals_begin);
    coiterator<key_iterator,val_iterator> end = coiterator<key_iterator,val_iterator>(keys_end, vals_end);
    sort(begin, end, cocomparator<key_iterator,val_iterator,Comparator>(comp));
}

}

inline std::complex<float> operator*(std::complex<float> f, double d)
{
    return f*(float)d;
}

inline std::complex<float> operator*(double d, std::complex<float> f)
{
    return f*(float)d;
}

inline std::complex<float> operator/(std::complex<float> f, double d)
{
    return f/(float)d;
}

template <class F, class I>
typename std::enable_if<std::is_integral<I>::value,std::complex<F> >::type
operator*(std::complex<F> f, I i)
{
    return f*(F)i;
}

template <class F, class I>
typename std::enable_if<std::is_integral<I>::value,std::complex<F> >::type
operator*(I i, std::complex<F> f)
{
    return f*(F)i;
}

template <class F, class I>
typename std::enable_if<std::is_integral<I>::value,std::complex<F> >::type
operator/(std::complex<F> f, I i)
{
    return f/(F)i;
}

#endif

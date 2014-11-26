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
#include <tuple>

#include "memory/memory.h"

#if !(defined(__GXX_EXPERIMENTAL_CXX0X__) || _MSC_VER >= 1600 || __cplusplus >= 201103l)
#error "A C++11-capable compiler is required."
#endif

#include <type_traits>
#include <memory>

#define NON_COPYABLE(name) \
public: \
    name(const name&) = delete; \
private:

#define NON_ASSIGNABLE(name) \
public: \
    name& operator=(const name&) = delete; \
private:

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

namespace detail
{

    template <size_t I, typename... Args>
    struct min_size_helper
    {
        size_t operator()(const tuple<Args...>& v)
        {
            return min(get<I-1>(v).size(), min_size_helper<I-1, Args...>()(v));
        }
    };

    template <typename... Args>
    struct min_size_helper<1, Args...>
    {
        size_t operator()(const tuple<Args...>& v)
        {
            return get<0>(v).size();
        }
    };

    template <typename... Args>
    struct min_size_helper<0, Args...>
    {
        size_t operator()(const tuple<Args...>& v)
        {
            return 0;
        }
    };

    template <typename... Args>
    size_t min_size(const tuple<Args...>& v)
    {
        return min_size_helper<sizeof...(Args), Args...>()(v);
    }

    template <size_t I, typename... Args>
    struct cbegin_helper
    {
        void operator()(tuple<typename decay<Args>::type::const_iterator...>& i,
                        const tuple<Args...>& v)
        {
            get<I-1>(i) = get<I-1>(v).begin();
            cbegin_helper<I-1, Args...>()(i, v);
        }
    };

    template <typename... Args>
    struct cbegin_helper<1, Args...>
    {
        void operator()(tuple<typename decay<Args>::type::const_iterator...>& i,
                        const tuple<Args...>& v)
        {
            get<0>(i) = get<0>(v).begin();
        }
    };

    template <typename... Args>
    struct cbegin_helper<0, Args...>
    {
        void operator()(tuple<typename decay<Args>::type::const_iterator...>& i,
                        const tuple<Args...>& v) {}
    };

    template <typename... Args>
    tuple<typename decay<Args>::type::const_iterator...> cbegin(const tuple<Args...>& v)
    {
        tuple<typename decay<Args>::type::const_iterator...> i;
        cbegin_helper<sizeof...(Args), Args...>()(i, v);
        return i;
    }

    template <size_t I, typename... Args>
    struct increment_helper
    {
        void operator()(tuple<typename vector<Args>::const_iterator...>& i)
        {
            ++get<I-1>(i);
            increment_helper<I-1, Args...>()(i);
        }
    };

    template <typename... Args>
    struct increment_helper<1, Args...>
    {
        void operator()(tuple<typename vector<Args>::const_iterator...>& i)
        {
            ++get<0>(i);
        }
    };

    template <typename... Args>
    struct increment_helper<0, Args...>
    {
        void operator()(tuple<typename vector<Args>::const_iterator...>& i) {}
    };

    template <typename... Args>
    void increment(tuple<typename vector<Args>::const_iterator...>& i)
    {
        increment_helper<sizeof...(Args), Args...>()(i);
    }

    template <int... S> struct integer_sequence {};
    template <int N, int... S> struct static_range : static_range<N-1, N-1, S...> {};
    template <int... S> struct static_range<0, S...> : integer_sequence<S...> {};

    template <typename... Args>
    struct call_helper
    {
        template <typename Func, int... S>
        call_helper(Func func, const tuple<Args...>& args, integer_sequence<S...> seq)
        {
            func(get<S>(args)...);
        }
    };

    template <size_t I, typename... Args>
    struct not_end_helper
    {
        bool operator()(const tuple<typename decay<Args>::type::const_iterator...>& i,
                    const tuple<Args...>& v)
        {
            return get<I-1>(i) != get<I-1>(v).end() &&
                   not_end_helper<I-1, Args...>()(i, v);
        }
    };

    template <typename... Args>
    struct not_end_helper<1, Args...>
    {
        bool operator()(const tuple<typename decay<Args>::type::const_iterator...>& i,
                    const tuple<Args...>& v)
        {
            return get<0>(i) != get<0>(v).end();
        }
    };

    template <typename... Args>
    struct not_end_helper<0, Args...>
    {
        bool operator()(const tuple<typename decay<Args>::type::const_iterator...>& i,
                        const tuple<Args...>& v)
        {
            return false;
        }
    };

    template <typename... Args>
    bool not_end(const tuple<typename decay<Args>::type::const_iterator...>& i,
                 const tuple<Args...>& v)
    {
        return not_end_helper<sizeof...(Args), Args...>()(i, v);
    }

    template <size_t I, typename... Args>
    struct print_tuple_helper : print_tuple_helper<I-1, Args...>
    {
        print_tuple_helper(ostream& os, const tuple<Args...>& t)
        : print_tuple_helper<I-1, Args...>(os, t)
        {
            os << get<I-1>(t);
            if (I < sizeof...(Args)) os << ", ";
        }
    };

    template <typename... Args>
    struct print_tuple_helper<1, Args...>
    {
        print_tuple_helper(ostream& os, const tuple<Args...>& t)
        {
            os << get<0>(t);
            if (1 < sizeof...(Args)) os << ", ";
        }
    };

    template <typename... Args>
    struct print_tuple_helper<0, Args...>
    {
        print_tuple_helper(ostream& os, const tuple<Args...>& t) {}
    };


    template <size_t I, typename... Args>
    struct reserve_helper
    {
        void operator()(tuple<vector<Args>...>& t, size_t n)
        {
            get<I-1>(t).reserve(n);
            reserve_helper<I-1, Args...>()(t, n);
        }
    };

    template <typename... Args>
    struct reserve_helper<1, Args...>
    {
        void operator()(tuple<vector<Args>...>& t, size_t n)
        {
            get<0>(t).reserve(n);
        }
    };

    template <typename... Args>
    struct reserve_helper<0, Args...>
    {
        void operator()(tuple<vector<Args>...>& t, size_t n) {}
    };

    template <typename... Args>
    void reserve(tuple<vector<Args>...>& t, size_t n)
    {
        reserve_helper<sizeof...(Args), Args...>()(t, n);
    }

    template <size_t I, typename... Args>
    struct emplace_back_helper
    {
        void operator()(tuple<vector<Args>...>& t, const tuple<Args...>& v)
        {
            get<I-1>(t).emplace_back(get<I-1>(v));
            emplace_back_helper<I-1, Args...>()(t, v);
        }

        void operator()(tuple<vector<Args>...>& t, tuple<Args...>&& v)
        {
            get<I-1>(t).emplace_back(move(get<I-1>(v)));
            emplace_back_helper<I-1, Args...>()(t, v);
        }
    };

    template <typename... Args>
    struct emplace_back_helper<1, Args...>
    {
        void operator()(tuple<vector<Args>...>& t, const tuple<Args...>& v)
        {
            get<0>(t).emplace_back(get<0>(v));
        }

        void operator()(tuple<vector<Args>...>& t, tuple<Args...>&& v)
        {
            get<0>(t).emplace_back(move(get<0>(v)));
        }
    };

    template <typename... Args>
    struct emplace_back_helper<0, Args...>
    {
        void operator()(tuple<vector<Args>...>& t, const tuple<Args...>& v) {}

        void operator()(tuple<vector<Args>...>& t, tuple<Args...>&& v) {}
    };

    template <typename... Args>
    void emplace_back(tuple<vector<Args>...>& t, const tuple<Args...>& v)
    {
        emplace_back_helper<sizeof...(Args), Args...>()(t, v);
    }

    template <typename... Args>
    void emplace_back(tuple<vector<Args>...>& t, tuple<Args...>&& v)
    {
        emplace_back_helper<sizeof...(Args), Args...>()(t, move(v));
    }

}

template <class T, class U>
struct if_exists
{
    typedef U type;
};

template<typename T> class global_ptr
{
    private:
        shared_ptr<T*> ptr;

    public:
        global_ptr() : ptr(new T*(NULL)) {}

        global_ptr(T* ptr) : ptr(new T*(ptr)) {}

        ~global_ptr()
        {
            if (unique() && get()) delete get();
        }

        global_ptr& operator=(const global_ptr& other)
        {
            reset(const_cast<T*>(other.get()));
            return *this;
        }

        void swap(global_ptr& other)
        {
            ptr.swap(other.ptr);
        }

        friend void swap(global_ptr& p1, global_ptr& p2)
        {
            p1.swap(p2);
        }

        long use_count() const
        {
            return ptr.use_count();
        }

        bool unique() const
        {
            return ptr.unique();
        }

        void reset()
        {
            if (get()) delete get();
            *ptr == nullptr;
        }

        void reset(T* p)
        {
            if (p == get()) return;
            if (get()) delete get();
            *ptr = p;
        }

        T* get()
        {
            return *ptr;
        }

        const T* get() const
        {
            return *ptr;
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

        operator bool() const
        {
            return get();
        }
};

inline string strprintf(const char* fmt, ...)
{
    va_list list;

    va_start(list, fmt);
    char fake[1];
    int n = vsnprintf(fake, 1, fmt, list);
    va_end(list);

    vector<char> s(n);
    va_start(list, fmt);
    vsnprintf(s.data(), n, fmt, list);
    va_end(list);

    return string(s.begin(), s.end());
}

template<typename T> ostream& operator<<(ostream& os, const vector<T>& v)
{
    os << "[";
    if (!v.empty()) os << v[0];
    for (int i = 1;i < v.size();i++) os << ", " << v[i];
    os << "]";
    return os;
}

template<typename T> string str(const T& t)
{
    ostringstream oss;
    oss << t;
    return oss.str();
}

template <typename T>
typename T::value_type max(const T& t)
{
    typedef typename T::value_type V;

    if (t.empty()) return V();

    typename T::const_iterator i = t.begin();
    V v = *i;
    for (;i != t.end();++i) if (v < *i) v = *i;

    return v;
}

template <typename T>
typename T::value_type min(const T& t)
{
    typedef typename T::value_type V;

    if (t.empty()) return V();

    typename T::const_iterator i = t.begin();
    V v = *i;
    for (;i != t.end();++i) if (*i < v) v = *i;

    return v;
}

template <typename Func, typename... Args>
void call(Func func, const tuple<Args...>& args)
{
    detail::call_helper<Args...>(func, args, detail::static_range<sizeof...(Args)>());
}

template <typename... Args>
ostream& operator<<(ostream& os, const tuple<Args...>& t)
{
    os << '{';
    detail::print_tuple_helper<sizeof...(Args), Args...>(os, t);
    os << '}';
    return os;
}

template <typename... Args>
vector<tuple<Args...>> zip(const tuple<const vector<Args>&...>& v)
{
    vector<tuple<Args...>> t;
    t.reserve(detail::min_size(v));

    auto i = detail::cbegin(v);
    for (;detail::not_end(i,v);detail::increment<Args...>(i))
    {
        call([&t](typename vector<Args>::const_iterator... args) {t.emplace_back(*args...); }, i);
    }

    return t;
}

template <typename... Args>
vector<tuple<Args...>> zip(const tuple<vector<Args>...>& v_)
{
    return zip(tuple<const vector<Args>&...>(v_));
}

template <typename... Args>
vector<tuple<Args...>> zip(const vector<Args>&... v_)
{
    return zip(tuple<const vector<Args>&...>(v_...));
}

template <typename... Args>
vector<tuple<Args...>> zip(const tuple<vector<Args>&&...>& v)
{
    vector<tuple<Args...>> t;
    t.reserve(detail::min_size(v));

    auto i = detail::cbegin(v);
    for (;detail::not_end(i,v);detail::increment<Args...>(i))
    {
        call([&t](typename vector<Args>::const_iterator... args) {t.emplace_back(move(*args)...); }, i);
    }

    return t;
}

template <typename... Args>
vector<tuple<Args...>> zip(tuple<vector<Args>...>&& v)
{
    return zip(tuple<vector<Args>&&...>(move(v)));
}

template <typename... Args>
vector<tuple<Args...>> zip(vector<Args>&&... v)
{
    return zip(forward_as_tuple(move(v)...));
}

template <typename... Args>
tuple<vector<Args>...> unzip(const vector<tuple<Args...>>& v)
{
    tuple<vector<Args>...> t;
    detail::reserve(t, v.size());

    for (auto i = v.begin();i != v.end();++i)
    {
        detail::emplace_back(t, *i);
    }

    return t;
}

template <typename... Args>
tuple<vector<Args>...> unzip(vector<tuple<Args...>>&& v)
{
    tuple<vector<Args>...> t;
    detail::reserve(t, v.size());

    for (auto i = v.begin();i != v.end();++i)
    {
        detail::emplace_back(t, move(*i));
    }

    return t;
}

template <typename Container, typename Functor, typename=void>
struct __erase
{
    static void erase(Container& v, const Functor& f)
    {
        v.erase(remove_if(v.begin(), v.end(), f), v.end());
    }
};

template <typename Container, typename T>
struct __erase<Container, T, typename enable_if<is_same<typename Container::value_type,T>::value>::type>
{
    static void erase(Container& v, const T& e)
    {
        v.erase(remove(v.begin(), v.end(), e), v.end());
    }
};

template <typename Container, typename T_or_Functor>
void erase(Container& v, const T_or_Functor& x)
{
    __erase<Container, T_or_Functor>::erase(v, x);
}

template <typename T> vector<T> range(T to)
{
    vector<T> r;
    r.reserve((typename vector<T>::size_type)to);

    for (T i = T();i < to;++i)
    {
        r.push_back(i);
    }

    return r;
}

template <typename T> vector<T> range(T from, T to)
{
    vector<T> r;
    r.reserve((typename vector<T>::size_type)(to-from));

    for (T i = from;i < to;++i)
    {
        r.push_back(i);
    }

    return r;
}

template<typename T> vector<T> operator+(const vector<T>& v1, const vector<T>& v2)
{
    vector<T> r(v1.size() + v2.size());
    copy(v1.begin(), v1.end(), r.begin());
    copy(v2.begin(), v2.end(), r.begin() + v1.size());
    return r;
}

template<typename T> vector<T> operator+(const vector<T>& v, const T& t)
{
    vector<T> r(v.size()+1);
    copy(v.begin(), v.end(), r.begin());
    r[v.size()] = t;
    return r;
}

template<typename T> vector<T>& operator+=(vector<T>& v1, const vector<T>& v2)
{
    v1.insert(v1.end(), v2.begin(), v2.end());
    return v1;
}

template<typename T> vector<T>& operator+=(vector<T>& v, const T& t)
{
    v.insert(v.end(), t);
    return v;
}

template<typename T> vector<T> slice(const vector<T>& v, int e1, int e2)
{
    return vector<T>(v.begin()+e1, v.begin()+e2);
}

template<typename T> vector<T> slice(const vector<T>& v, int e1)
{
    return vector<T>(v.begin()+e1, v.end());
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

template<typename T, class Predicate> T& filter(T& v, Predicate pred)
{
    auto i1 = v.begin();
    for (auto i2 = v.begin();i2 != v.end();++i2)
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

template<typename T, class Functor> auto apply(T& v, Functor f)
    -> vector<decltype(f(v.back()))>
{
    typedef decltype(f(v.back())) U;
    vector<U> v2();
    for (auto& i : v)
    {
        v2.emplace_back(f(i));
    }
    return v2;
}

template<typename T, class Predicate> typename decay<T>::type filter_copy(const T& v, Predicate pred)
{
    typename decay<T>::type v2;
    for (auto& i : v)
    {
        if (pred(i)) v2.emplace_back(i);
    }
    return v2;
}

template<typename T, class Predicate> typename decay<T>::type filter_copy(T&& v, Predicate pred)
{
    typename decay<T>::type v2;
    for (auto& i : v)
    {
        if (pred(i)) v2.emplace_back(move(i));
    }
    return v2;
}

template<typename T> typename T::value_type sum(const T& v)
{
    typedef typename T::value_type U;
    U s = U();
    for (auto& i : v) s += i;
    return s;
}

template<typename T, typename U> bool contains(const T& v, const U& e)
{
    return find(v.begin(), v.end(), e) != v.end();
}

template<typename T> T& sort(T& v)
{
    sort(v.begin(), v.end());
    return v;
}

template<typename T, typename Compare> T& sort(T& v, Compare comp)
{
    sort(v.begin(), v.end(), comp);
    return v;
}

template<typename T> typename decay<T>::type sorted(T&& v)
{
    typename decay<T>::type v2(forward<T>(v));
    sort(v2);
    return v2;
}

template<typename T, typename Compare> typename decay<T>::type sorted(T&& v, Compare comp)
{
    typename decay<T>::type v2(forward<T>(v));
    sort(v2, comp);
    return v2;
}

template<typename T> T& uniq(T& v)
{
    sort(v.begin(), v.end());
    auto i1 = unique(v.begin(), v.end());
    v.resize(i1-v.begin());

    return v;
}

template<typename T> typename decay<T>::type uniq_copy(T&& v)
{
    typename decay<T>::type v2(forward<T>(v));
    uniq(v2);
    return v2;
}

template<typename T> T intersection(const T& v1, const T& v2)
{
    T v;
    T v3(v1);
    T v4(v2);

    v.resize(v1.size()+v2.size());

    sort(v3.begin(),v3.end());
    sort(v4.begin(),v4.end());

    auto end = set_intersection(v3.begin(), v3.end(), v4.begin(), v4.end(), v.begin());
    v.resize(end-v.begin());

    return v;
}

/*
 * Exclude elements of v2 from v1
 */
template<typename T> T& exclude(T& v1, const T& v2)
{
    T v3(v2);

    sort(v1.begin(), v1.end());
    sort(v3.begin(), v3.end());

    auto i1 = v1.begin();
    auto i2 = v1.begin();
    auto i3 = v3.begin();
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
template<typename T> typename decay<T>::type exclude_copy(T&& v1, const T& v2)
{
    typename decay<T>::type v3(forward<T>(v1));
    exclude(v3, v2);
    return v3;
}

template<typename T, typename U> T& mask(T& v, const U& mask)
{
    auto i1 = v.begin();
    auto i2 = v.begin();
    auto i3 = mask.begin();
    for (;i2 != v.end();++i2,++i3)
    {
        if (*i3)
        {
            swap(*i1, *i2);
            ++i1;
        }
    }

    v.resize(i1-v.begin());
    return v;
}

template<typename T, typename U> typename decay<T>::type mask_copy(const T& v, const U& mask)
{
    typename decay<T>::type v2;

    auto i3 = mask.begin();
    for (auto& i : v)
    {
        if (*i3++)
        {
            v2.emplace_back(i);
        }
    }

    return v2;
}

template<typename T, typename U> typename decay<T>::type mask_copy(T&& v, const U& mask)
{
    typename decay<T>::type v2;

    auto i3 = mask.begin();
    for (auto& i : v)
    {
        if (*i3++)
        {
            v2.emplace_back(move(i));
        }
    }

    return v2;
}

template<typename T>
T& translate(T& s, const T& from, const T& to)
{
    assert(from.size() == to.size());

    //vector<T> fromSorted(from);
    //vector<T> toSorted(to);
    //cosort(from, to);
    // not working?

    T fromSorted, toSorted;
    tie(fromSorted, toSorted) = unzip(sorted(zip(from, to)));

    for (auto& l : s)
    {
        auto lb = lower_bound(fromSorted.begin(), fromSorted.end(), l);

        if (lb != fromSorted.end() && *lb == l)
        {
            l = toSorted[lb - fromSorted.begin()];
        }
    }

    return s;
}

template<typename T>
typename decay<T>::type translate_copy(T&& s, const T& from, const T& to)
{
    typename decay<T>::type s_(forward<T>(s));
    translate(s_, from, to);
    return s_;
}

inline string& translate(string& s, const string& from, const string& to)
{
    assert(from.size() == to.size());

    unsigned char trans[256];
    if (s.size() < 256)
    {
        for (int i = 0;i < s.size();i++) trans[(unsigned char)s[i]] = s[i];
    }
    else
    {
        for (int i = 0;i < 256;i++) trans[i] = i;
    }

    for (int i = 0;i < from.size();i++) trans[(unsigned char)from[i]] = to[i];
    for (int i = 0;i <    s.size();i++) s[i] = trans[(unsigned char)s[i]];

    return s;
}

inline string translate_copy(string&& s, const string& from, const string& to)
{
    string s_(move(s));
    translate(s_, from, to);
    return s_;
}

inline string translate_copy(const string& s, const string& from, const string& to)
{
    string s_(s);
    translate(s_, from, to);
    return s_;
}

inline string toupper(string&& s)
{
    string S(move(s));
    for (auto& C : S) C = ::toupper(C);
    return S;
}

inline string tolower(string&& S)
{
    string s(move(S));
    for (auto& c : s) c = ::tolower(c);
    return s;
}

inline string toupper(const string& s)
{
    string S(s);
    for (auto& C : S) C = ::toupper(C);
    return S;
}

inline string tolower(const string& S)
{
    string s(S);
    for (auto& c : s) c = ::tolower(c);
    return s;
}

inline       float conj(      float v) { return v; }
inline      double conj(     double v) { return v; }
inline long double conj(long double v) { return v; }

inline       float real(      float v) { return v; }
inline      double real(     double v) { return v; }
inline long double real(long double v) { return v; }

inline       float imag(      float v) { return 0.0f; }
inline      double imag(     double v) { return 0.0; }
inline long double imag(long double v) { return 0.0l; }

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

template <class T, int ndim>
class tensor;

template <class T, int ndim, int dim>
class tensor_ref
{
    friend class tensor_ref<T,ndim,dim-1>;
    friend class tensor<T,ndim>;

    private:
        tensor_ref& operator=(const tensor_ref& other);

    protected:
        tensor<T, ndim>& array;
        size_t idx;

        tensor_ref(const tensor_ref& other)
        : array(other.array), idx(other.idx) {}

        tensor_ref(tensor<T, ndim>& array, size_t idx, int i);

    public:
        tensor_ref<T,ndim,dim+1> operator[](int i);

        const tensor_ref<T,ndim,dim+1> operator[](int i) const;

        //operator T*();

        //operator const T*() const;
};

template <class T, int ndim>
class tensor_ref<T, ndim, ndim>
{
    friend class tensor_ref<T,ndim,ndim-1>;
    friend class tensor<T,ndim>;

    private:
        tensor_ref& operator=(const tensor_ref& other);

    protected:
        tensor<T, ndim>& array;
        size_t idx;

        tensor_ref(const tensor_ref& other)
        : array(other.array), idx(other.idx) {}

        tensor_ref(tensor<T, ndim>& array, size_t idx, int i);

    public:
        T& operator[](int i);

        const T& operator[](int i) const;

        //operator T*();

        //operator const T*() const;
};

template <class T, int ndim>
class tensor
{
    template<class T_, int ndim_, int dim_> friend class tensor_ref;

    friend void swap(tensor& a, tensor& b)
    {
        swap(a.data_, b.data_);
        swap(a.len, b.len);
        swap(a.stride, b.stride);
    }

    protected:
        T* data_;
        vector<int> len;
        vector<size_t> stride;

    public:
        enum Layout {COLUMN_MAJOR, ROW_MAJOR};

        tensor(const tensor& other)
        : len(other.len), stride(other.stride)
        {
            size_t num = 1;
            for (int i = 0;i < ndim;i++) num *= len[i];
            data_ = new T[num];
            copy(other.data_, other.data_+num, data_);
        }

        explicit tensor(const vector<int>& len = vector<int>(ndim,0),
                        const T& val = T(), Layout layout = COLUMN_MAJOR)
        : len(len), stride(ndim)
        {
            assert(len.size() == ndim);

            if (layout == ROW_MAJOR)
            {
                stride[ndim-1] = 1;
                for (int i = ndim-2;i >= 0;i--)
                {
                    stride[i] = stride[i+1]*len[i+1];
                }
                data_ = new T[stride[0]*len[0]];
                fill(data_, data_+stride[0]*len[0], val);
            }
            else
            {
                stride[0] = 1;
                for (int i = 1;i < ndim;i++)
                {
                    stride[i] = stride[i-1]*len[i-1];
                }
                data_ = new T[stride[ndim-1]*len[ndim-1]];
                fill(data_, data_+stride[ndim-1]*len[ndim-1], val);
            }
        }

        ~tensor()
        {
            delete[] data_;
        }

        void resize(const vector<int>& len, const T& val = T(), Layout layout = COLUMN_MAJOR)
        {
            delete[] data_;
            assert(len.size() == ndim);
            this->len = len;

            if (layout == ROW_MAJOR)
            {
                stride[ndim-1] = 1;
                for (int i = ndim-2;i >= 0;i--)
                {
                    stride[i] = stride[i+1]*len[i+1];
                }
                data_ = new T[stride[0]*len[0]];
                fill(data_, data_+stride[0]*len[0], val);
            }
            else
            {
                stride[0] = 1;
                for (int i = 1;i < ndim;i++)
                {
                    stride[i] = stride[i-1]*len[i-1];
                }
                data_ = new T[stride[ndim-1]*len[ndim-1]];
                fill(data_, data_+stride[ndim-1]*len[ndim-1], val);
            }
        }

        void clear()
        {
            resize(vector<int>(ndim, 0));
        }

        tensor& operator=(tensor other)
        {
            swap(*this, other);
            return *this;
        }

        tensor_ref<T,ndim,2> operator[](int i)
        {
            assert(i >= 0 && i < len[0]);
            return tensor_ref<T,ndim,2>(*this, (size_t)0, i);
        }

        const tensor_ref<T,ndim,2> operator[](int i) const
        {
            assert(i >= 0 && i < len[0]);
            return tensor_ref<T,ndim,2>(const_cast<tensor&>(*this), (size_t)0, i);
        }

        double* data() { return data_; }

        const double* data() const { return data_; }

        operator T*(){ return data_; }

        operator const T*() const { return data_; }

        const vector<int>& length() const { return len; }

        const int size() const
        {
            int sz = len[0];
            for (int i = 1;i < ndim;i++)
            {
                sz *= len[i];
            }
            return sz;
        }
};

template <class T>
class tensor<T,0>
{
    protected:
        T data_;

    public:
        tensor(const T& val = T()) : data_(val) {}

        double& data() { return data_; }

        const double& data() const { return data_; }

        operator T&() { return data_; }

        operator const T&() const { return data_; }
};

template <class T>
class tensor<T,1>
{
    friend void swap(tensor& a, tensor& b)
    {
        swap(a.data_, b.data_);
        swap(a.len, b.len);
    }

    protected:
        T* data_;
        int len;

    public:
        tensor(const tensor& other)
        : len(other.len)
        {
            data_ = new T[len];
            copy(other.data_, other.data_+len, data_);
        }

        explicit tensor(int n = 0, const T& val = T())
        : len(n)
        {
            data_ = new T[len];
            fill(data_, data_+len, val);
        }

        ~tensor()
        {
            delete[] data_;
        }

        void resize(int n, const T& val = T())
        {
            len = n;
            delete[] data_;
            data_ = new T[len];
            fill(data_, data_+len, val);
        }

        void clear()
        {
            resize(0);
        }

        tensor& operator=(tensor other)
        {
            swap(*this, other);
            return *this;
        }

        T& operator[](int i)
        {
            assert(i >= 0 && i < len);
            return data_[i];
        }

        const T& operator[](int i) const
        {
            assert(i >= 0 && i < len);
            return data_[i];
        }

        double* data() { return data_; }

        const double* data() const { return data_; }

        operator T*() { return data_; }

        operator const T*() const { return data_; }

        int length() const { return len; }

        int size() const { return len; }
};

template <class T, int ndim, int dim>
tensor_ref<T,ndim,dim>::tensor_ref(tensor<T, ndim>& array, size_t idx, int i)
: array(array), idx(idx+i*array.stride[dim-2]) {}

template <class T, int ndim, int dim>
tensor_ref<T,ndim,dim+1> tensor_ref<T,ndim,dim>::operator[](int i)
{
    assert(i >= 0 && i < array.len[dim-1]);
    return tensor_ref<T,ndim,dim+1>(array, idx, i);
}

template <class T, int ndim, int dim>
const tensor_ref<T,ndim,dim+1> tensor_ref<T,ndim,dim>::operator[](int i) const
{
    assert(i >= 0 && i < array.len[dim-1]);
    return tensor_ref<T,ndim,dim+1>(array, idx, i);
}

//template <class T, int ndim, int dim>
//tensor_ref<T,ndim,dim>::operator T*()
//{
//    return array.data_+idx;
//}

//template <class T, int ndim, int dim>
//tensor_ref<T,ndim,dim>::operator const T*() const
//{
//    return array.data_+idx;
//}

template <class T, int ndim>
tensor_ref<T,ndim,ndim>::tensor_ref(tensor<T, ndim>& array, size_t idx, int i)
: array(array), idx(idx+i*array.stride[ndim-2]) {}

template <class T, int ndim>
T& tensor_ref<T,ndim,ndim>::operator[](int i)
{
    assert(i >= 0 && i < array.len[ndim-1]);
    return array.data_[idx+i*array.stride[ndim-1]];
}

template <class T, int ndim>
const T& tensor_ref<T,ndim,ndim>::operator[](int i) const
{
    assert(i >= 0 && i < array.len[ndim-1]);
    return array.data_[idx+i*array.stride[ndim-1]];
}

//template <class T, int ndim>
//tensor_ref<T,ndim,ndim>::operator T*()
//{
//    return array.data_+idx;
//}

//template <class T, int ndim>
//tensor_ref<T,ndim,ndim>::operator const T*() const
//{
//    return array.data_+idx;
//}

template <class T>
class matrix : public tensor<T,2>
{
    public:
        using typename tensor<T,2>::Layout;
        using tensor<T,2>::ROW_MAJOR;
        using tensor<T,2>::COLUMN_MAJOR;

        matrix() {}

        matrix(int m, int n, const T& val = T(), Layout layout = ROW_MAJOR)
        : tensor<T,2>({m,n}, val, layout) {}

        void resize(int m, int n, const T& val = T(), Layout layout = ROW_MAJOR)
        {
            tensor<T,2>::resize({m,n}, val, layout);
        }
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
    coiterator<key_iterator,val_iterator> begin(keys_begin, vals_begin);
    coiterator<key_iterator,val_iterator> end  (keys_end  , vals_end  );
    sort(begin, end);
}

template <class key_iterator, class val_iterator, class Comparator>
void cosort(key_iterator keys_begin, key_iterator keys_end,
            val_iterator vals_begin, val_iterator vals_end,
            Comparator comp)
{
    coiterator<key_iterator,val_iterator> begin(keys_begin, vals_begin);
    coiterator<key_iterator,val_iterator> end  (keys_end  , vals_end  );
    sort(begin, end, cocomparator<key_iterator,val_iterator,Comparator>(comp));
}

template <class Keys, class Values>
void cosort(Keys keys, Values values)
{
    cosort(keys.begin(), keys.end(), values.begin(), values.end());
}

template <class Keys, class Values, class Comparator>
void cosort(Keys keys, Values values, Comparator comp)
{
    cosort(keys.begin(), keys.end(), values.begin(), values.end(), comp);
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

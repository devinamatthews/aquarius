#ifndef _AQUARIUS_STL_EXT_HPP_
#define _AQUARIUS_STL_EXT_HPP_

#include <array>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
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
#include <cassert>
#include <type_traits>
#include <memory>
#include <map>
#include <list>
#include <set>
#include <iomanip>
#include <cstdint>
#include <cstddef>

#if !(defined(__GXX_EXPERIMENTAL_CXX0X__) || _MSC_VER >= 1600 || __cplusplus >= 201103l)
#error "A C++11-capable compiler is required."
#endif

#define ENABLE_IF_CONST(const_type,return_type) \
template <typename _IsConst = const_type > \
typename enable_if<is_const<_IsConst>::value, return_type >::type

#define ENABLE_IF_NON_CONST(const_type,return_type) \
template <typename _IsConst = const_type > \
typename enable_if<!is_const<_IsConst>::value, return_type >::type

//#define ENABLE_IF_SAME(old_type,new_type,return_type) \
//template <typename new_type > \
//typename enable_if<is_base_of<const old_type, const new_type >::value, return_type >::type

#define ENABLE_IF_SAME(old_type,new_type,return_type) \
template <typename new_type > return_type

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
    using std::pair;
    using std::make_pair;
    using std::make_tuple;
    using std::get;

    using std::complex;
    using std::numeric_limits;

    using std::min;
    using std::max;
    using std::abs;
    using std::copy;
    using std::fill;
    using std::sort;
    using std::swap;
    using std::remove;

    using std::unique_ptr;
    using std::shared_ptr;
    using std::make_shared;

    using std::move;
    using std::forward;
    using std::decay;
    using std::declval;
    using std::enable_if;
    using std::is_same;
    using std::is_const;
    using std::is_integral;
    using std::is_floating_point;
    using std::remove_const;
    using std::conditional;
    using std::initializer_list;

    using std::iterator;
    using std::iterator_traits;
    using std::reverse_iterator;
    using std::forward_iterator_tag;
    using std::bidirectional_iterator_tag;
    using std::random_access_iterator_tag;

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
    using std::fixed;
    using std::skipws;
    using std::showpos;
    using std::setw;
    using std::streamsize;
    using std::streambuf;

    using std::toupper;
    using std::tolower;

    using std::exception;
    using std::logic_error;
    using std::runtime_error;
    using std::out_of_range;

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

            global_ptr(const global_ptr&) = default;

            global_ptr(global_ptr&&) = default;

            global_ptr(T* ptr) : ptr(new T*(ptr)) {}

            ~global_ptr()
            {
                if (unique() && get()) delete get();
            }

            global_ptr& operator=(const global_ptr& other) = default;

            global_ptr& operator=(global_ptr&& other) = default;

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

    namespace detail
    {

        template <typename vector_>
        class ptr_vector_
        {
            private:
                template <typename ptr_>
                class iterator_
                {
                    friend class ptr_vector_;

                    public:
                        typedef typename iterator_traits<ptr_>::iterator_category iterator_category;
                        typedef typename iterator_traits<ptr_>::value_type value_type;
                        typedef typename iterator_traits<ptr_>::difference_type difference_type;
                        typedef typename iterator_traits<ptr_>::pointer pointer;
                        typedef typename iterator_traits<ptr_>::reference reference;

                        typedef typename conditional<is_same<ptr_,value_type*>::value,
                                                     typename vector_::iterator,
                                                     typename vector_::const_iterator>::type ptr_iterator;

                        iterator_() = default;

                        iterator_(const iterator_& other)
                        : it_(other.it_) {}

                        template <typename ptr__>
                        iterator_(typename enable_if<is_same<ptr__,value_type*>::value,
                                                     const iterator_<ptr__>&>::type other)
                        : it_(other.it_) {}

                        iterator_& operator=(const iterator_& other)
                        {
                            it_ = other.it_;
                            return *this;
                        }

                        template <typename ptr__>
                        iterator_& operator=(typename enable_if<is_same<ptr__,value_type*>::value,
                                             const iterator_<ptr__>&>::type other)
                        {
                            it_ = other.it_;
                            return *this;
                        }

                        bool operator==(const iterator_<value_type*>& x)
                        {
                            return it_ == x.it_;
                        }

                        bool operator!=(const iterator_<value_type*>& x)
                        {
                            return it_ != x.it_;
                        }

                        bool operator==(const iterator_<const value_type*>& x)
                        {
                            return it_ == x.it_;
                        }

                        bool operator!=(const iterator_<const value_type*>& x)
                        {
                            return it_ != x.it_;
                        }

                        reference operator*() const
                        {
                            return **it_;
                        }

                        pointer operator->() const
                        {
                            return it_->get();
                        }

                        iterator_& operator++()
                        {
                            ++it_;
                            return *this;
                        }

                        iterator_& operator--()
                        {
                            --it_;
                            return *this;
                        }

                        iterator_ operator++(int x)
                        {
                            return iterator_(it_++);
                        }

                        iterator_ operator--(int x)
                        {
                            return iterator_(it_--);
                        }

                        iterator_ operator+(difference_type n) const
                        {
                            return iterator_(it_+n);
                        }

                        iterator_ operator-(difference_type n) const
                        {
                            return iterator_(it_-n);
                        }

                        friend iterator_ operator+(difference_type n, const iterator_& x)
                        {
                            return x+n;
                        }

                        difference_type operator-(const iterator_<value_type*>& other) const
                        {
                            return it_-other.it_;
                        }

                        difference_type operator-(const iterator_<const value_type*>& other) const
                        {
                            return it_-other.it_;
                        }

                        bool operator<(const iterator_<value_type*>& x) const
                        {
                            return it_ < x.it_;
                        }

                        bool operator>(const iterator_<value_type*>& x) const
                        {
                            return it_ > x.it_;
                        }

                        bool operator<=(const iterator_<value_type*>& x) const
                        {
                            return it_ <= x.it_;
                        }

                        bool operator>=(const iterator_<value_type*>& x) const
                        {
                            return it_ >= x.it_;
                        }

                        bool operator<(const iterator_<const value_type*>& x) const
                        {
                            return it_ < x.it_;
                        }

                        bool operator>(const iterator_<const value_type*>& x) const
                        {
                            return it_ > x.it_;
                        }

                        bool operator<=(const iterator_<const value_type*>& x) const
                        {
                            return it_ <= x.it_;
                        }

                        bool operator>=(const iterator_<const value_type*>& x) const
                        {
                            return it_ >= x.it_;
                        }

                        iterator_& operator+=(difference_type n)
                        {
                            it_ += n;
                            return *this;
                        }

                        iterator_& operator-=(difference_type n)
                        {
                            it_ -= n;
                            return *this;
                        }

                        reference operator[](difference_type n) const
                        {
                            return *it_[n];
                        }

                        friend void swap(iterator_& a, iterator_& b)
                        {
                            using aquarius::swap;
                            swap(a.it_, b.it_);
                        }

                        ptr_iterator base()
                        {
                            return it_;
                        }

                    protected:
                        iterator_(ptr_iterator it) : it_(it) {}

                        ptr_iterator it_;
                };

            public:
                typedef typename vector_::value_type ptr_type;
                typedef typename decay<decltype(*declval<ptr_type>())>::type value_type;
                typedef value_type& reference;
                typedef const value_type& const_reference;
                typedef value_type* pointer;
                typedef const value_type* const_pointer;
                typedef iterator_<pointer> iterator;
                typedef iterator_<const_pointer> const_iterator;
                typedef std::reverse_iterator<iterator> reverse_iterator;
                typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
                typedef typename iterator_traits<iterator>::difference_type difference_type;
                typedef typename vector_::size_type size_type;

                typedef typename vector_::iterator ptr_iterator;
                typedef typename vector_::const_iterator const_ptr_iterator;
                typedef typename vector_::reverse_iterator reverse_ptr_iterator;
                typedef typename vector_::const_reverse_iterator const_reverse_ptr_iterator;

                ptr_vector_() {}

                explicit ptr_vector_(size_type n) : impl_(n) {}

                ptr_vector_(size_type n, const value_type& val)
                {
                    assign(n, val);
                }

                ptr_vector_(size_type n, value_type&& val)
                {
                    assign(n, forward<value_type>(val));
                }

                ptr_vector_(const ptr_vector_&) = delete;

                ptr_vector_(ptr_vector_&&) = default;

                ptr_vector_(initializer_list<value_type> il)
                {
                    assign(il);
                }

                ptr_vector_(initializer_list<pointer> il)
                {
                    assign(il);
                }

                ptr_vector_& operator=(const ptr_vector_&) = delete;

                ptr_vector_& operator=(ptr_vector_&&) = default;

                iterator begin()
                {
                    return iterator(impl_.begin());
                }

                iterator end()
                {
                    return iterator(impl_.end());
                }

                const_iterator begin() const
                {
                    return const_iterator(impl_.begin());
                }

                const_iterator end() const
                {
                    return const_iterator(impl_.end());
                }

                reverse_iterator rbegin()
                {
                    return reverse_iterator(impl_.rbegin());
                }

                reverse_iterator rend()
                {
                    return reverse_iterator(impl_.rend());
                }

                const_reverse_iterator rbegin() const
                {
                    return const_reverse_iterator(impl_.rbegin());
                }

                const_reverse_iterator rend() const
                {
                    return const_reverse_iterator(impl_.rend());
                }

                const_iterator cbegin() const
                {
                    return const_iterator(impl_.begin());
                }

                const_iterator cend() const
                {
                    return const_iterator(impl_.end());
                }

                const_reverse_iterator crbegin() const
                {
                    return const_reverse_iterator(impl_.rbegin());
                }

                const_reverse_iterator crend() const
                {
                    return const_reverse_iterator(impl_.rend());
                }

                ptr_iterator pbegin()
                {
                    return impl_.begin();
                }

                ptr_iterator pend()
                {
                    return impl_.end();
                }

                const_ptr_iterator pbegin() const
                {
                    return impl_.begin();
                }

                const_ptr_iterator pend() const
                {
                    return impl_.end();
                }

                reverse_ptr_iterator rpbegin()
                {
                    return impl_.rbegin();
                }

                reverse_ptr_iterator rpend()
                {
                    return impl_.rend();
                }

                const_reverse_ptr_iterator rpbegin() const
                {
                    return impl_.rbegin();
                }

                const_reverse_ptr_iterator rpend() const
                {
                    return impl_.rend();
                }

                const_ptr_iterator cpbegin() const
                {
                    return impl_.begin();
                }

                const_ptr_iterator cpend() const
                {
                    return impl_.end();
                }

                const_reverse_ptr_iterator crpbegin() const
                {
                    return impl_.rbegin();
                }

                const_reverse_ptr_iterator crpend() const
                {
                    return impl_.rend();
                }

                size_type size() const
                {
                    return impl_.size();
                }

                size_type max_size() const
                {
                    return impl_.max_size();
                }

                void resize(size_type n)
                {
                    impl_.resize(n);
                }

                void resize(size_type n, const value_type& x)
                {
                    if (n <= impl_.size())
                    {
                        impl_.resize(n);
                    }
                    else
                    {
                        impl_.reserve(n);

                        for (size_type i = impl_.size();i < n;i++)
                        {
                            impl_.emplace_back(new value_type(x));
                        }
                    }
                }

                void resize(size_type n, value_type&& x)
                {
                    if (n <= impl_.size())
                    {
                        impl_.resize(n);
                    }
                    else
                    {
                        impl_.reserve(n);

                        impl_.emplace_back(new value_type(move(x)));

                        for (size_type i = impl_.size();i < n;i++)
                        {
                            impl_.emplace_back(new value_type(x));
                        }
                    }
                }

                size_type capacity() const
                {
                    return impl_.capacity();
                }

                bool empty() const
                {
                    return impl_.empty();
                }

                void reserve(size_type n)
                {
                    impl_.reserve(n);
                }

                void shrink_to_fit()
                {
                    impl_.shrink_to_fit();
                }

                reference operator[](size_type n)
                {
                    return *impl_[n];
                }

                const_reference operator[](size_type n) const
                {
                    return *impl_[n];
                }

                reference at(size_type n)
                {
                    return *impl_.at(n);
                }

                const reference at(size_type n) const
                {
                    return *impl_.at(n);
                }

                ptr_type& ptr(size_type n)
                {
                    return impl_.at(n);
                }

                const ptr_type& ptr(size_type n) const
                {
                    return impl_.at(n);
                }

                reference front()
                {
                    return *impl_.front();
                }

                reference back()
                {
                    return *impl_.back();
                }

                const_reference front() const
                {
                    return *impl_.front();
                }

                const_reference back() const
                {
                    return *impl_.back();
                }

                ptr_type& pfront()
                {
                    return impl_.front();
                }

                ptr_type& pback()
                {
                    return impl_.back();
                }

                const ptr_type& pfront() const
                {
                    return impl_.front();
                }

                const ptr_type& pback() const
                {
                    return impl_.back();
                }

                void assign(size_type n, const value_type& val)
                {
                    impl_.clear();
                    resize(n, val);
                }

                void assign(size_type n, value_type&& val)
                {
                    impl_.clear();
                    resize(n, move(val));
                }

                void assign(const ptr_vector_& x)
                {
                    *this = x;
                }

                void assign(ptr_vector_&& x)
                {
                    *this = move(x);
                }

                void assign(initializer_list<value_type> il)
                {
                    impl_.clear();
                    impl_.reserve(il.size());

                    for (auto& val : il)
                    {
                        impl_.emplace_back(new value_type(val));
                    }
                }

                void assign(initializer_list<pointer> il)
                {
                    impl_.clear();
                    impl_.reserve(il.size());

                    for (auto& ptr : il)
                    {
                        impl_.emplace_back(ptr);
                    }
                }

                void push_back(const value_type& x)
                {
                    impl_.emplace_back(new value_type(x));
                }

                void push_back(value_type&& x)
                {
                    impl_.emplace_back(new value_type(move(x)));
                }

                void push_back(const ptr_type& x)
                {
                    impl_.push_back(x);
                }

                void push_back(ptr_type&& x)
                {
                    impl_.push_back(move(x));
                }

                template <typename Pointer>
                typename enable_if< is_same<pointer, Pointer>::value &&
                                   !is_same<pointer,ptr_type>::value>::type
                push_back(Pointer x)
                {
                    impl_.emplace_back(x);
                }

                void pop_back()
                {
                    impl_.pop_back();
                }

                iterator insert(const_iterator position, const value_type& val)
                {
                    impl_.emplace(position.it_, new value_type(val));
                    return *this;
                }

                iterator insert(const_iterator position, value_type&& val)
                {
                    impl_.emplace(position.it_, new value_type(move(val)));
                    return *this;
                }

                iterator insert(const_iterator position, const ptr_type& val)
                {
                    impl_.insert(position.it_, val);
                    return *this;
                }

                iterator insert(const_iterator position, ptr_type&& val)
                {
                    impl_.insert(position.it_, move(val));
                    return *this;
                }

                iterator insert(const_iterator position, pointer val)
                {
                    impl_.emplace(position.it_, val);
                    return *this;
                }

                iterator insert(const_iterator position, size_type n, const value_type& val)
                {
                    impl_.reserve(impl_.size()+n);
                    auto middle = impl_.end();
                    impl_.resize(impl_.size()+n, val);
                    rotate(impl_.begin(), middle, impl_.end());
                    return *this;
                }

                iterator insert(const_iterator position, size_type n, value_type&& val)
                {
                    impl_.reserve(impl_.size()+n);
                    auto middle = impl_.end();
                    impl_.resize(impl_.size()+n, move(val));
                    rotate(impl_.begin(), middle, impl_.end());
                    return *this;
                }

                iterator insert(const_iterator position, initializer_list<value_type> il)
                {
                    impl_.reserve(impl_.size()+il.size());
                    auto middle = impl_.end();

                    for (auto& val : il)
                    {
                        impl_.emplace_back(new value_type(val));
                    }

                    rotate(impl_.begin(), middle, impl_.end());
                    return *this;
                }

                iterator insert(const_iterator position, initializer_list<pointer> il)
                {
                    impl_.reserve(impl_.size()+il.size());
                    auto middle = impl_.end();

                    for (auto& ptr : il)
                    {
                        impl_.emplace_back(ptr);
                    }

                    rotate(impl_.begin(), middle, impl_.end());
                    return *this;
                }

                iterator erase(const_iterator position)
                {
                    return iterator(impl_.erase(position.it_));
                }

                iterator erase(const_iterator first, const_iterator last)
                {
                    return iterator(impl_.erase(first.it_, last.it_));
                }

                ptr_iterator perase(const_ptr_iterator position)
                {
                    /*
                     * Work around stupid bug in libstdc++ and libc++
                     */
                    ptr_iterator pos = impl_.begin()+(position-impl_.cbegin());
                    return impl_.erase(pos);
                }

                ptr_iterator perase(const_ptr_iterator first, const_ptr_iterator last)
                {
                    /*
                     * Work around stupid bug in libstdc++ and libc++
                     */
                    ptr_iterator f = impl_.begin()+(first-impl_.cbegin());
                    ptr_iterator l = impl_.begin()+(last-impl_.cbegin());
                    return impl_.erase(f, l);
                }

                void swap(ptr_vector_& x)
                {
                    impl_.swap(x.impl_);
                }

                void clear()
                {
                    impl_.clear();
                }

                template <typename... Args>
                void emplace(const_iterator position, Args&&... args)
                {
                    impl_.emplace(position.it_, new value_type(forward<Args>(args)...));
                }

                template <typename... Args>
                void emplace_back(Args&&... args)
                {
                    impl_.emplace_back(new value_type(forward<Args>(args)...));
                }

                friend bool operator==(const ptr_vector_& lhs, const ptr_vector_& rhs)
                {
                    if (lhs.size() != rhs.size()) return false;
                    return equal(lhs.begin(), lhs.end(), rhs.begin());
                }

                friend bool operator!=(const ptr_vector_& lhs, const ptr_vector_& rhs)
                {
                    return !(lhs == rhs);
                }

                friend bool operator<(const ptr_vector_& lhs, const ptr_vector_& rhs)
                {
                    return lexicographical_compare(lhs.begin(), lhs.end(),
                                                   rhs.begin(), rhs.end());
                }

                friend bool operator>(const ptr_vector_& lhs, const ptr_vector_& rhs)
                {
                    return rhs < lhs;
                }

                friend bool operator<=(const ptr_vector_& lhs, const ptr_vector_& rhs)
                {
                    return !(rhs < lhs);
                }

                friend bool operator>=(const ptr_vector_& lhs, const ptr_vector_& rhs)
                {
                    return !(lhs < rhs);
                }

                friend void swap(ptr_vector_&& a, ptr_vector_&& b)
                {
                    a.swap(b);
                }

                friend void swap(ptr_vector_& a, ptr_vector_& b)
                {
                    a.swap(b);
                }

            private:
                vector_ impl_;
        };

        template <typename list_>
        class ptr_list_
        {
            private:
                template <typename ptr_, typename=void>
                class iterator_
                {
                    friend class ptr_list_;

                    public:
                        typedef typename iterator_traits<ptr_>::iterator_category iterator_category;
                        typedef typename iterator_traits<ptr_>::value_type value_type;
                        typedef typename iterator_traits<ptr_>::difference_type difference_type;
                        typedef typename iterator_traits<ptr_>::pointer pointer;
                        typedef typename iterator_traits<ptr_>::reference reference;

                        typedef typename conditional<is_same<ptr_,value_type*>::value,
                                                     typename list_::iterator,
                                                     typename list_::const_iterator>::type ptr_iterator;

                        iterator_() = default;

                        iterator_(const iterator_& other)
                        : it_(other.it_) {}

                        template <typename ptr__>
                        iterator_(const iterator_<ptr__, typename enable_if<is_same<ptr__,value_type*>::value>::type>& other)
                        : it_(other.it_) {}

                        iterator_& operator=(const iterator_& other)
                        {
                            it_ = other.it_;
                            return *this;
                        }

                        template <typename ptr__>
                        iterator_& operator=(const iterator_<ptr__, typename enable_if<is_same<ptr__,value_type*>::value>::type>& other)
                        {
                            it_ = other.it_;
                            return *this;
                        }

                        bool operator==(const iterator_<value_type*>& x)
                        {
                            return it_ == x.it_;
                        }

                        bool operator!=(const iterator_<value_type*>& x)
                        {
                            return it_ != x.it_;
                        }

                        bool operator==(const iterator_<const value_type*>& x)
                        {
                            return it_ == x.it_;
                        }

                        bool operator!=(const iterator_<const value_type*>& x)
                        {
                            return it_ != x.it_;
                        }

                        reference operator*() const
                        {
                            return **it_;
                        }

                        pointer operator->() const
                        {
                            return it_->get();
                        }

                        iterator_& operator++()
                        {
                            ++it_;
                            return *this;
                        }

                        iterator_& operator--()
                        {
                            --it_;
                            return *this;
                        }

                        iterator_ operator++(int x)
                        {
                            return iterator_(it_++);
                        }

                        iterator_ operator--(int x)
                        {
                            return iterator_(it_--);
                        }

                        friend void swap(iterator_& a, iterator_& b)
                        {
                            using aquarius::swap;
                            swap(a.it_, b.it_);
                        }

                        ptr_iterator base()
                        {
                            return it_;
                        }

                    protected:
                        iterator_(ptr_iterator it) : it_(it) {}

                        ptr_iterator it_;
                };

                /*
                 * Workaround for lack of support for const_iterator in
                 * erase(), splice(), etc.
                 */
                typename list_::iterator ci2i(list_& l, typename list_::const_iterator ci)
                {
                    auto i = l.begin();
                    while (&*i != &*ci) ++i;
                    return i;
                }

            public:
                typedef typename list_::value_type ptr_type;
                typedef typename decay<decltype(*declval<ptr_type>())>::type value_type;
                typedef value_type& reference;
                typedef const value_type& const_reference;
                typedef value_type* pointer;
                typedef const value_type* const_pointer;
                typedef iterator_<pointer> iterator;
                typedef iterator_<const_pointer> const_iterator;
                typedef std::reverse_iterator<iterator> reverse_iterator;
                typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
                typedef typename iterator_traits<iterator>::difference_type difference_type;
                typedef typename list_::size_type size_type;

                typedef typename list_::iterator ptr_iterator;
                typedef typename list_::const_iterator const_ptr_iterator;
                typedef typename list_::reverse_iterator reverse_ptr_iterator;
                typedef typename list_::const_reverse_iterator const_reverse_ptr_iterator;

                ptr_list_() {}

                explicit ptr_list_(size_type n) : impl_(n) {}

                ptr_list_(size_type n, const value_type& val)
                {
                    assign(n, val);
                }

                ptr_list_(size_type n, value_type&& val)
                {
                    assign(n, forward<value_type>(val));
                }

                ptr_list_(const ptr_list_&) = default;

                ptr_list_(ptr_list_&&) = default;

                ptr_list_(initializer_list<value_type> il)
                {
                    assign(il);
                }

                ptr_list_(initializer_list<pointer> il)
                {
                    assign(il);
                }

                ptr_list_& operator=(const ptr_list_&) = default;

                ptr_list_& operator=(ptr_list_&&) = default;

                iterator begin()
                {
                    return iterator(impl_.begin());
                }

                iterator end()
                {
                    return iterator(impl_.end());
                }

                const_iterator begin() const
                {
                    return const_iterator(impl_.begin());
                }

                const_iterator end() const
                {
                    return const_iterator(impl_.end());
                }

                reverse_iterator rbegin()
                {
                    return reverse_iterator(impl_.rbegin());
                }

                reverse_iterator rend()
                {
                    return reverse_iterator(impl_.rend());
                }

                const_reverse_iterator rbegin() const
                {
                    return const_reverse_iterator(impl_.rbegin());
                }

                const_reverse_iterator rend() const
                {
                    return const_reverse_iterator(impl_.rend());
                }

                const_iterator cbegin() const
                {
                    return const_iterator(impl_.begin());
                }

                const_iterator cend() const
                {
                    return const_iterator(impl_.end());
                }

                const_reverse_iterator crbegin() const
                {
                    return const_reverse_iterator(impl_.rbegin());
                }

                const_reverse_iterator crend() const
                {
                    return const_reverse_iterator(impl_.rend());
                }

                ptr_iterator pbegin()
                {
                    return impl_.begin();
                }

                ptr_iterator pend()
                {
                    return impl_.end();
                }

                const_ptr_iterator pbegin() const
                {
                    return impl_.begin();
                }

                const_ptr_iterator pend() const
                {
                    return impl_.end();
                }

                reverse_ptr_iterator rpbegin()
                {
                    return impl_.rbegin();
                }

                reverse_ptr_iterator rpend()
                {
                    return impl_.rend();
                }

                const_reverse_ptr_iterator rpbegin() const
                {
                    return impl_.rbegin();
                }

                const_reverse_ptr_iterator rpend() const
                {
                    return impl_.rend();
                }

                const_ptr_iterator cpbegin() const
                {
                    return impl_.begin();
                }

                const_ptr_iterator cpend() const
                {
                    return impl_.end();
                }

                const_reverse_ptr_iterator crpbegin() const
                {
                    return impl_.rbegin();
                }

                const_reverse_ptr_iterator crpend() const
                {
                    return impl_.rend();
                }

                size_type size() const
                {
                    return impl_.size();
                }

                size_type max_size() const
                {
                    return impl_.max_size();
                }

                void resize(size_type n)
                {
                    impl_.resize(n);
                }

                void resize(size_type n, const value_type& x)
                {
                    if (n <= impl_.size())
                    {
                        impl_.resize(n);
                    }
                    else
                    {
                        impl_.reserve(n);

                        for (size_type i = impl_.size();i < n;i++)
                        {
                            impl_.emplace_back(new value_type(x));
                        }
                    }
                }

                void resize(size_type n, value_type&& x)
                {
                    if (n <= impl_.size())
                    {
                        impl_.resize(n);
                    }
                    else
                    {
                        impl_.reserve(n);

                        impl_.emplace_back(new value_type(move(x)));

                        for (size_type i = impl_.size();i < n;i++)
                        {
                            impl_.emplace_back(new value_type(x));
                        }
                    }
                }

                bool empty() const
                {
                    return impl_.empty();
                }

                reference front()
                {
                    return *impl_.front();
                }

                reference back()
                {
                    return *impl_.back();
                }

                const_reference front() const
                {
                    return *impl_.front();
                }

                const_reference back() const
                {
                    return *impl_.back();
                }

                ptr_type& pfront()
                {
                    return impl_.front();
                }

                ptr_type& pback()
                {
                    return impl_.back();
                }

                const ptr_type& pfront() const
                {
                    return impl_.front();
                }

                const ptr_type& pback() const
                {
                    return impl_.back();
                }

                void assign(size_type n, const value_type& val)
                {
                    impl_.clear();
                    resize(n, val);
                }

                void assign(size_type n, value_type&& val)
                {
                    impl_.clear();
                    resize(n, move(val));
                }

                void assign(const ptr_list_& x)
                {
                    *this = x;
                }

                void assign(ptr_list_&& x)
                {
                    *this = move(x);
                }

                void assign(initializer_list<value_type> il)
                {
                    impl_.clear();

                    for (auto& val : il)
                    {
                        impl_.emplace_back(new value_type(val));
                    }
                }

                void assign(initializer_list<pointer> il)
                {
                    impl_.clear();

                    for (auto& ptr : il)
                    {
                        impl_.emplace_back(ptr);
                    }
                }

                void push_front(const value_type& x)
                {
                    impl_.emplace_front(new value_type(x));
                }

                void push_front(value_type&& x)
                {
                    impl_.emplace_front(new value_type(move(x)));
                }

                void push_front(const ptr_type& x)
                {
                    impl_.push_front(x);
                }

                void push_front(ptr_type&& x)
                {
                    impl_.push_front(move(x));
                }

                template <typename Pointer>
                typename enable_if< is_same<pointer, Pointer>::value &&
                                   !is_same<pointer,ptr_type>::value>::type
                push_front(Pointer x)
                {
                    impl_.emplace_front(x);
                }

                void pop_front()
                {
                    impl_.pop_front();
                }

                void push_back(const value_type& x)
                {
                    impl_.emplace_back(new value_type(x));
                }

                void push_back(value_type&& x)
                {
                    impl_.emplace_back(new value_type(move(x)));
                }

                void push_back(const ptr_type& x)
                {
                    impl_.push_back(x);
                }

                void push_back(ptr_type&& x)
                {
                    impl_.push_back(move(x));
                }

                template <typename Pointer>
                typename enable_if< is_same<pointer, Pointer>::value &&
                                   !is_same<pointer,ptr_type>::value>::type
                push_back(Pointer x)
                {
                    impl_.emplace_back(x);
                }

                void pop_back()
                {
                    impl_.pop_back();
                }

                iterator insert(const_iterator position, const value_type& val)
                {
                    impl_.emplace(ci2i(impl_, position.it_), new value_type(val));
                    return *this;
                }

                iterator insert(const_iterator position, value_type&& val)
                {
                    impl_.emplace(ci2i(impl_, position.it_), new value_type(move(val)));
                    return *this;
                }

                iterator insert(const_iterator position, const ptr_type& val)
                {
                    impl_.insert(ci2i(impl_, position.it_), val);
                    return *this;
                }

                iterator insert(const_iterator position, ptr_type&& val)
                {
                    impl_.insert(ci2i(impl_, position.it_), move(val));
                    return *this;
                }

                iterator insert(const_iterator position, pointer val)
                {
                    impl_.emplace(ci2i(impl_, position.it_), val);
                    return *this;
                }

                iterator insert(const_iterator position, size_type n, const value_type& val)
                {
                    while (n-- > 0) impl_.emplace(ci2i(impl_, position.it_), val);
                    return *this;
                }

                iterator insert(const_iterator position, size_type n, value_type&& val)
                {
                    if (n > 0)
                    {
                        while (n-- > 1) impl_.emplace(ci2i(impl_, position.it_), val);
                        impl_.emplace(ci2i(impl_, position.it_), move(val));
                    }
                    return *this;
                }

                iterator insert(const_iterator position, initializer_list<value_type> il)
                {
                    for (auto& val : il)
                    {
                        impl_.emplace(ci2i(impl_, position.it_), new value_type(val));
                    }
                    return *this;
                }

                iterator insert(const_iterator position, initializer_list<pointer> il)
                {
                    for (auto& ptr : il)
                    {
                        impl_.emplace(ci2i(impl_, position.it_), ptr);
                    }
                    return *this;
                }

                iterator erase(const_iterator position)
                {
                    return iterator(impl_.erase(ci2i(impl_, position.it_)));
                }

                iterator erase(const_iterator first, const_iterator last)
                {
                    return iterator(impl_.erase(ci2i(impl_, first.it_), ci2i(impl_, last.it_)));
                }

                ptr_iterator perase(const_ptr_iterator position)
                {
                    return impl_.erase(ci2i(impl_, position));
                }

                ptr_iterator perase(const_ptr_iterator first, const_ptr_iterator last)
                {
                    return impl_.erase(ci2i(impl_, first), ci2i(impl_, last));
                }

                void swap(ptr_list_& x)
                {
                    impl_.swap(x.impl_);
                }

                void clear()
                {
                    impl_.clear();
                }

                template <typename... Args>
                void emplace(const_iterator position, Args&&... args)
                {
                    impl_.emplace(ci2i(impl_, position.it_), new value_type(forward<Args>(args)...));
                }

                template <typename... Args>
                void emplace_front(Args&&... args)
                {
                    impl_.emplace_front(new value_type(forward<Args>(args)...));
                }

                template <typename... Args>
                void emplace_back(Args&&... args)
                {
                    impl_.emplace_back(new value_type(forward<Args>(args)...));
                }

                friend bool operator==(const ptr_list_& lhs, const ptr_list_& rhs)
                {
                    if (lhs.size() != rhs.size()) return false;
                    return equal(lhs.begin(), lhs.end(), rhs.begin());
                }

                friend bool operator!=(const ptr_list_& lhs, const ptr_list_& rhs)
                {
                    return !(lhs == rhs);
                }

                friend bool operator<(const ptr_list_& lhs, const ptr_list_& rhs)
                {
                    return lexicographical_compare(lhs.begin(), lhs.end(),
                                                   rhs.begin(), rhs.end());
                }

                friend bool operator>(const ptr_list_& lhs, const ptr_list_& rhs)
                {
                    return rhs < lhs;
                }

                friend bool operator<=(const ptr_list_& lhs, const ptr_list_& rhs)
                {
                    return !(rhs < lhs);
                }

                friend bool operator>=(const ptr_list_& lhs, const ptr_list_& rhs)
                {
                    return !(lhs < rhs);
                }

                friend void swap(ptr_list_&& a, ptr_list_&& b)
                {
                    a.swap(b);
                }

                friend void swap(ptr_list_& a, ptr_list_& b)
                {
                    a.swap(b);
                }

                void splice(const_iterator position, ptr_list_& x)
                {
                    impl_.splice(ci2i(impl_, position.it_), x.impl_);
                }

                void splice (const_iterator position, ptr_list_&& x)
                {
                    impl_.splice(ci2i(impl_, position.it_), move(x.impl_));
                }

                void splice (const_iterator position, ptr_list_& x, const_iterator i)
                {
                    impl_.splice(ci2i(impl_, position.it_), x.impl_, ci2i(x, i.it_));
                }

                void splice (const_iterator position, ptr_list_&& x, const_iterator i)
                {
                    impl_.splice(ci2i(impl_, position.it_), move(x.impl_), ci2i(x, i.it_));
                }

                void splice (const_iterator position, ptr_list_& x,
                             const_iterator first, const_iterator last)
                {
                    impl_.splice(ci2i(impl_, position.it_), x.impl_, ci2i(x, first.it_), ci2i(x, last.it_));
                }

                void splice (const_iterator position, ptr_list_&& x,
                             const_iterator first, const_iterator last)
                {
                    impl_.splice(ci2i(impl_, position.it_), move(x.impl_), ci2i(x, first.it_), ci2i(x, last.it_));
                }

                void remove(const value_type& val)
                {
                    impl_.remove_if([&val](const ptr_type& ptr) { return *ptr == val; });
                }

                template <typename Predicate>
                void remove_if(Predicate pred)
                {
                    impl_.remove_if([&pred](const ptr_type& ptr) { return pred(*ptr); });
                }

                void unique()
                {
                    impl_.unique([](const ptr_type& a, const ptr_type& b) { return *a == *b; });
                }

                template <class BinaryPredicate>
                void unique(BinaryPredicate binary_pred)
                {
                    impl_.unique([&binary_pred](const ptr_type& a, const ptr_type& b)
                                     { return binary_pred(*a, *b); });
                }

                void merge(ptr_list_& x)
                {
                    impl_.merge(x.impl_, [](const ptr_type& a, const ptr_type& b)
                                            { return *a < *b; });
                }

                void merge(ptr_list_&& x)
                {
                    impl_.merge(move(x.impl_), [](const ptr_type& a, const ptr_type& b)
                                                    { return *a < *b; });
                }

                template <class Compare>
                void merge(ptr_list_& x, Compare comp)
                {
                    impl_.merge(x.impl_, [&comp](const ptr_type& a, const ptr_type& b)
                                                { return comp(*a, *b); });
                }

                template <class Compare>
                void merge(ptr_list_&& x, Compare comp)
                {
                    impl_.merge(move(x.impl_), [&comp](const ptr_type& a, const ptr_type& b)
                                                    { return comp(*a, *b); });
                }

                void sort()
                {
                    impl_.sort([](const ptr_type& a, const ptr_type& b)
                                   { return *a < *b; });
                }

                template <class Compare>
                void sort(Compare comp)
                {
                    impl_.sort([&comp](const ptr_type& a, const ptr_type& b)
                                    { return comp(*a, *b); });
                }

                void reverse() noexcept
                {
                    impl_.reverse();
                }

            private:
                list_ impl_;
        };

    }

    template <typename T>
    using ptr_vector = detail::ptr_vector_<vector<T*>>;

    template <typename T>
    using unique_vector = detail::ptr_vector_<vector<unique_ptr<T>>>;

    template <typename T>
    using shared_vector = detail::ptr_vector_<vector<shared_ptr<T>>>;

    template <typename T>
    using global_vector = detail::ptr_vector_<vector<global_ptr<T>>>;

    template <typename T>
    using ptr_list = detail::ptr_list_<list<T*>>;

    template <typename T>
    using unique_list = detail::ptr_list_<list<unique_ptr<T>>>;

    template <typename T>
    using shared_list = detail::ptr_list_<list<shared_ptr<T>>>;

    template <typename T>
    using global_list = detail::ptr_list_<list<global_ptr<T>>>;

    namespace detail
    {

    template <size_t N, typename T, typename... Ts>
    struct vec_helper
    {
        vec_helper(vector<typename decay<T>::type>& v, T&& t, Ts&&... ts)
        {
            v.push_back(forward<T>(t));
            vec_helper<N-1, Ts...>(v, forward<Ts>(ts)...);
        }
    };

    template <typename T>
    struct vec_helper<1, T>
    {
        vec_helper(vector<typename decay<T>::type>& v, T&& t)
        {
            v.push_back(forward<T>(t));
        }
    };

    }

    template <typename T, typename... Ts>
    vector<typename decay<T>::type> vec(T&& t, Ts&&... ts)
    {
        vector<typename decay<T>::type> v;
        v.reserve(1+sizeof...(Ts));
        detail::vec_helper<1+sizeof...(Ts), T, Ts...>(v, forward<T>(t), forward<Ts>(ts)...);
        return v;
    }

    inline string strprintf(const char* fmt, ...)
    {
        va_list list;

        va_start(list, fmt);
        char fake[1];
        int n = vsnprintf(fake, 1, fmt, list);
        va_end(list);

        vector<char> s(n+1);
        va_start(list, fmt);
        vsnprintf(s.data(), n, fmt, list);
        va_end(list);

        return string(s.begin(), s.end());
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

    template <typename T1, typename T2>
    void swap(T1&& a, T2&& b)
    {
        swap(a, b);
    }

    namespace detail
    {

        template <class T, class U>
        struct doublet
        {
            T first;
            U second;

            doublet(const T& first, const U& second) : first(first), second(second) {}

            friend void swap(doublet& first, doublet& second)
            {
                swap(first.first, second.first);
                swap(first.second, second.second);
            }

            friend void swap(doublet&& first, doublet&& second)
            {
                swap(first.first, second.first);
                swap(first.second, second.second);
            }

            doublet(const doublet<T&,U&>& other)
            : first(other.first), second(other.second) {}

            doublet(doublet<T&,U&>&& other)
            : first(move(other.first)), second(move(other.second)) {}

            doublet(const doublet<T,U>& other)
            : first(other.first), second(other.second) {}

            doublet(doublet<T,U>&& other)
            : first(move(other.first)), second(move(other.second)) {}

            doublet& operator=(const doublet<T&,U&>& other)
            {
                first = other.first;
                second = other.second;
                return *this;
            }

            doublet& operator=(doublet<T&,U&>&& other)
            {
                first = move(other.first);
                second = move(other.second);
                return *this;
            }

            doublet& operator=(const doublet<T,U>& other)
            {
                first = other.first;
                second = other.second;
                return *this;
            }

            doublet& operator=(doublet<T,U>&& other)
            {
                first = move(other.first);
                second = move(other.second);
                return *this;
            }

            bool operator==(const doublet<T&,U&>& other) const
            {
                return first == other.first;
            }

            bool operator!=(const doublet<T&,U&>& other) const
            {
                return first != other.first;
            }

            bool operator==(const doublet& other) const
            {
                return first == other.first;
            }

            bool operator!=(const doublet& other) const
            {
                return first != other.first;
            }

            bool operator<(const doublet<T&,U&>& other) const
            {
                return first < other.first;
            }

            bool operator>(const doublet<T&,U&>& other) const
            {
                return first > other.first;
            }

            bool operator<=(const doublet<T&,U&>& other) const
            {
                return first <= other.first;
            }

            bool operator>=(const doublet<T&,U&>& other) const
            {
                return first >= other.first;
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
        struct doublet<T&,U&>
        {
            T& first;
            U& second;

            doublet(T& first, U& second) : first(first), second(second) {}

            friend void swap(doublet& first, doublet& second)
            {
                swap(first.first, second.first);
                swap(first.second, second.second);
            }

            friend void swap(doublet&& first, doublet&& second)
            {
                swap(first.first, second.first);
                swap(first.second, second.second);
            }

            doublet(doublet<T&,U&>& other)
            : first(other.first), second(other.second) {}

            doublet(doublet<T&,U&>&& other)
            : first(other.first), second(other.second) {}

            doublet(doublet<T,U>& other)
            : first(other.first), second(other.second) {}

            doublet& operator=(const doublet<T&,U&>& other)
            {
                first = other.first;
                second = other.second;
                return *this;
            }

            doublet& operator=(doublet<T&,U&>&& other)
            {
                first = move(other.first);
                second = move(other.second);
                return *this;
            }

            doublet& operator=(const doublet<T,U>& other)
            {
                first = other.first;
                second = other.second;
                return *this;
            }

            doublet& operator=(doublet<T,U>&& other)
            {
                first = move(other.first);
                second = move(other.second);
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

            bool operator==(const doublet& other) const
            {
                return first == other.first;
            }

            bool operator!=(const doublet& other) const
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

            bool operator<(const doublet& other) const
            {
                return first < other.first;
            }

            bool operator>(const doublet& other) const
            {
                return first > other.first;
            }

            bool operator<=(const doublet& other) const
            {
                return first <= other.first;
            }

            bool operator>=(const doublet& other) const
            {
                return first >= other.first;
            }
        };

        template <class T, class U>
        class coiterator : public iterator<random_access_iterator_tag,
        //tuple<typename iterator_traits<T>::value_type,
        //      typename iterator_traits<U>::value_type>,
        //ptrdiff_t,
        //tuple<typename iterator_traits<T>::pointer,
        //      typename iterator_traits<U>::pointer>,
        //tuple<typename iterator_traits<T>::reference,
        //      typename iterator_traits<U>::reference>>
        doublet<typename iterator_traits<T>::value_type,
                typename iterator_traits<U>::value_type>,
        ptrdiff_t,
        doublet<typename iterator_traits<T>::pointer,
                typename iterator_traits<U>::pointer>,
        doublet<typename iterator_traits<T>::reference,
                typename iterator_traits<U>::reference>>
        {
            T it_T;
            U it_U;

            public:
                coiterator(const T& it_T, const U& it_U) : it_T(it_T), it_U(it_U) {}

                bool operator==(const coiterator& other) const
                {
                    return it_T == other.it_T;
                }

                bool operator!=(const coiterator& other) const
                {
                    return it_T != other.it_T;
                }

                bool operator<(const coiterator& other) const
                {
                    return it_T < other.it_T;
                }

                bool operator>(const coiterator& other) const
                {
                    return it_T > other.it_T;
                }

                bool operator<=(const coiterator& other) const
                {
                    return it_T <= other.it_T;
                }

                bool operator>=(const coiterator& other) const
                {
                    return it_T >= other.it_T;
                }

                typename coiterator::reference operator*()
                {
                    //return typename coiterator::reference(move(*it_T),move(*it_U));
                    return typename coiterator::reference(*it_T,*it_U);
                }

                typename coiterator::reference operator[](ptrdiff_t n)
                {
                    //return typename coiterator::reference(move(it_T[n]),move(it_U[n]));
                    return typename coiterator::reference(it_T[n],it_U[n]);
                }

                coiterator& operator++()
                {
                    ++it_T;
                    ++it_U;
                    return *this;
                }

                coiterator& operator--()
                {
                    --it_T;
                    --it_U;
                    return *this;
                }

                coiterator operator++(int x)
                {
                    return coiterator(it_T++, it_U++);
                }

                coiterator operator--(int x)
                {
                    return coiterator(it_T--, it_U--);
                }

                coiterator& operator+=(ptrdiff_t n)
                {
                    it_T += n;
                    it_U += n;
                    return *this;
                }

                coiterator& operator-=(ptrdiff_t n)
                {
                    it_T -= n;
                    it_U -= n;
                    return *this;
                }

                coiterator operator+(ptrdiff_t n) const
                {
                    return coiterator(it_T+n, it_U+n);
                }

                friend coiterator operator+(ptrdiff_t n, const coiterator& other)
                {
                    return coiterator(other.it_T+n, other.it_U+n);
                }

                coiterator operator-(ptrdiff_t n) const
                {
                    return coiterator(it_T-n, it_U-n);
                }

                ptrdiff_t operator-(const coiterator& other) const
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

    }

    template <class key_iterator, class val_iterator>
    void cosort(key_iterator keys_begin, key_iterator keys_end,
                val_iterator vals_begin, val_iterator vals_end)
    {
        detail::coiterator<key_iterator,val_iterator> begin(keys_begin, vals_begin);
        detail::coiterator<key_iterator,val_iterator> end  (keys_end  , vals_end  );
        sort(begin, end);
    }

    template <class key_iterator, class val_iterator, class Comparator>
    void cosort(key_iterator keys_begin, key_iterator keys_end,
                val_iterator vals_begin, val_iterator vals_end,
                Comparator comp)
    {
        detail::coiterator<key_iterator,val_iterator> begin(keys_begin, vals_begin);
        detail::coiterator<key_iterator,val_iterator> end  (keys_end  , vals_end  );
        sort(begin, end, detail::cocomparator<key_iterator,val_iterator,Comparator>(comp));
    }

    template <class Keys, class Values>
    void cosort(Keys& keys, Values& values)
    {
        cosort(keys.begin(), keys.end(), values.begin(), values.end());
    }

    template <class Keys, class Values, class Comparator>
    void cosort(Keys& keys, Values& values, Comparator comp)
    {
        cosort(keys.begin(), keys.end(), values.begin(), values.end(), comp);
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
        vector<U> v2;
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

        T fromSorted2(from);
        T toSorted2(to);
        cosort(fromSorted2, toSorted2);

        T fromSorted, toSorted;
        tie(fromSorted, toSorted) = unzip(sorted(zip(from, to)));

        assert(fromSorted2 == fromSorted);
        assert(toSorted2 == toSorted);

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

}

namespace std
{

    inline std::complex<float> operator*(const std::complex<float>& f, double d)
    {
        return f*(float)d;
    }

    inline std::complex<float> operator*(double d, const std::complex<float>& f)
    {
        return f*(float)d;
    }

    inline std::complex<float> operator/(const std::complex<float>& f, double d)
    {
        return f/(float)d;
    }

    inline std::complex<float> operator/(double d, const std::complex<float>& f)
    {
        return std::complex<float>((float)d)/f;
    }

    template <class F, class I>
    typename std::enable_if<std::is_integral<I>::value,std::complex<F> >::type
    operator*(const std::complex<F>& f, I i)
    {
        return f*(F)i;
    }

    template <class F, class I>
    typename std::enable_if<std::is_integral<I>::value,std::complex<F> >::type
    operator*(I i, const std::complex<F>& f)
    {
        return f*(F)i;
    }

    template <class F, class I>
    typename std::enable_if<std::is_integral<I>::value,std::complex<F> >::type
    operator/(const std::complex<F>& f, I i)
    {
        return f/(F)i;
    }

    template <class F, class I>
    typename std::enable_if<std::is_integral<I>::value,std::complex<F> >::type
    operator/(I i, const std::complex<F>& f)
    {
        return std::complex<F>((F)i)/f;
    }

    template<typename T> std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
    {
        os << "[";
        if (!v.empty()) os << v[0];
        for (int i = 1;i < v.size();i++) os << ", " << v[i];
        os << "]";
        return os;
    }

    template<typename T> std::vector<T> operator+(const std::vector<T>& v1, const std::vector<T>& v2)
    {
        std::vector<T> r(v1.size() + v2.size());
        std::copy(v1.begin(), v1.end(), r.begin());
        std::copy(v2.begin(), v2.end(), r.begin() + v1.size());
        return r;
    }

    template<typename T> std::vector<T> operator+(const std::vector<T>& v, const T& t)
    {
        std::vector<T> r(v.size()+1);
        std::copy(v.begin(), v.end(), r.begin());
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

}

#endif

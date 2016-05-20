#ifndef _MPIWRAP_WINDOW_HPP_
#define _MPIWRAP_WINDOW_HPP_

#if MPIWRAP_HAVE_MPI_WIN

#include "mpiwrap_common.hpp"
#include "mpiwrap_info.hpp"

namespace MPIWrap
{

namespace internal
{
    template <typename Derived> class Comm;
}

class Window
{
    friend class Intercomm;
    friend class Intracomm;
    template <typename Derived> friend class internal::Comm;

    protected:
        MPI_Win win;
        MPI_Aint size_;
        void* data_;

        explicit Window(const MPI_Win& win)
        : win(win)
        {
            MPI_Int flag;
            MPI_Win_get_attr(win, MPI_WIN_BASE, &data_, &flag);
            MPI_Win_get_attr(win, MPI_WIN_SIZE, &size_, &flag);
        }

    public:
#if MPIWRAP_CXX11

        Window(Window&& other)
        : win(other.win), size_(other.size_), data_(other.data_)
        {
            other.win = MPI_WIN_NULL;
            other.size_ = 0;
            other.data_ = NULL;
        }

#endif

        ~Window()
        {
            if (win != MPI_WIN_NULL)
            {
                MPIWRAP_CALL(MPI_Win_free(&win));
            }
        }

        operator MPI_Win&() { return win; }

        operator const MPI_Win&() const { return win; }

        operator MPI_Win*() { return &win; }

        operator const MPI_Win*() const { return &win; }

        MPI_Aint size() const { return size_; }

        void* data() { return data_; }

        const void* data() const { return data_; }

        void fence(MPI_Int assertion = 0) const
        {
            MPIWRAP_CALL(MPI_Win_fence(assertion, win));
        }

        void start(const Group& group, MPI_Int assertion = 0) const
        {
            MPIWRAP_CALL(MPI_Win_start(group, assertion, win));
        }

        void post(const Group& group, MPI_Int assertion = 0) const
        {
            MPIWRAP_CALL(MPI_Win_post(group, assertion, win));
        }

        void complete() const
        {
            MPIWRAP_CALL(MPI_Win_complete(win));
        }

        void wait() const
        {
            MPIWRAP_CALL(MPI_Win_wait(win));
        }

        bool test() const
        {
            MPI_Int flag;
            MPIWRAP_CALL(MPI_Win_test(win, &flag));
            return flag;
        }

        void lock(MPI_Int lock_type, MPI_Int rank, MPI_Int assertion = 0) const
        {
            MPIWRAP_CALL(MPI_Win_lock(lock_type, rank, assertion, win));
        }

        void unlock(MPI_Int rank) const
        {
            MPIWRAP_CALL(MPI_Win_unlock(rank, win));
        }

        template <typename T>
        void Put(const T* buffer, MPI_Int count, MPI_Int rank, MPI_Aint offset) const
        {
            Put(buffer, count, rank, offset, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Put(const std::vector<T>& buffer, MPI_Int rank, MPI_Aint offset) const
        {
            Put(buffer, rank, offset, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Put(const T* buffer, MPI_Int count, MPI_Int rank, MPI_Aint offset, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Put(nc(buffer), count, type, rank, offset, count, type, win));
        }

        template <typename T>
        void Put(const std::vector<T>& buffer, MPI_Int rank, MPI_Aint offset, const Datatype& type) const
        {
            Put(&buffer[0], buffer.size(), rank, offset, type);
        }

        template <typename T>
        void Get(T* buffer, MPI_Int count, MPI_Int rank, MPI_Aint offset) const
        {
            Get(buffer, count, rank, offset, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Get(std::vector<T>& buffer, MPI_Int rank, MPI_Aint offset) const
        {
            Get(buffer, rank, offset, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Get(T* buffer, MPI_Int count, MPI_Int rank, MPI_Aint offset, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Get(buffer, count, type, rank, offset, count, type, win));
        }

        template <typename T>
        void Get(std::vector<T>& buffer, MPI_Int rank, MPI_Aint offset, const Datatype& type) const
        {
            Get(&buffer[0], buffer.size(), rank, offset, type);
        }

        template <typename T>
        void Accumulate(const T* buffer, MPI_Int count, MPI_Int rank, MPI_Aint offset, const MPI_Op& op) const
        {
            Accumulate(buffer, count, rank, offset, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Accumulate(const std::vector<T>& buffer, MPI_Int rank, MPI_Aint offset, const MPI_Op& op) const
        {
            Accumulate(buffer, rank, offset, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Accumulate(const T* buffer, MPI_Int count, MPI_Int rank, MPI_Aint offset, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Accumulate(nc(buffer), count, type, rank, offset, count, type, op, win));
        }

        template <typename T>
        void Accumulate(const std::vector<T>& buffer, MPI_Int rank, MPI_Aint offset, const MPI_Op& op, const Datatype& type) const
        {
            Accumulate(&buffer[0], buffer.size(), rank, offset, op, type);
        }

#if MPIWRAP_VERSION_AT_LEAST(3,0)

        void lock(MPI_Int assertion = 0) const
        {
            MPIWRAP_CALL(MPI_Win_lock_all(assertion, win));
        }

        void unlock() const
        {
            MPIWRAP_CALL(MPI_Win_unlock_all(win));
        }

        void flush(MPI_Int rank) const
        {
            MPIWRAP_CALL(MPI_Win_flush(rank, win));
        }

        void flush() const
        {
            MPIWRAP_CALL(MPI_Win_flush_all(win));
        }

        void flush_local(MPI_Int rank) const
        {
            MPIWRAP_CALL(MPI_Win_flush_local(rank, win));
        }

        void flush_local() const
        {
            MPIWRAP_CALL(MPI_Win_flush_local_all(win));
        }

        void sync() const
        {
            MPIWRAP_CALL(MPI_Win_sync(win));
        }

        template <typename T>
        void Get_accumulate(const T* put, T* get, MPI_Int count, MPI_Int rank, MPI_Aint offset, const MPI_Op& op) const
        {
            Get_accumulate(put, get, count, rank, offset, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Get_accumulate(const std::vector<T>& put, std::vector<T>& get, MPI_Int rank, MPI_Aint offset, const MPI_Op& op) const
        {
            Get_accumulate(put, get, rank, offset, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Get_accumulate(const T* put, T* get, MPI_Int count, MPI_Int rank, MPI_Aint offset, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Get_accumulate(put, count, type, get, count, type, rank, offset, count, type, op, win));
        }

        template <typename T>
        void Get_accumulate(const std::vector<T>& put, std::vector<T>& get, MPI_Int rank, MPI_Aint offset, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_ASSERT(put.size() == get.size(),
                           "put and get buffers must have the same size.");
            Get_accumulate(&put[0], &get[0], put.size(), rank, offset, op, type);
        }

        template <typename T>
        T Fetch_and_op(const T& put, MPI_Int rank, MPI_Aint offset, const MPI_Op& op) const
        {
            T get;
            Datatype type = MPI_TYPE_<T>::value();
            MPIWRAP_CALL(MPI_Fetch_and_op(&put, &get, type, rank, offset, op, win));
            return get;
        }

        template <typename T>
        T Compare_and_swap(const T& put, const T& compare, MPI_Int rank, MPI_Aint offset) const
        {
            T get;
            Datatype type = MPI_TYPE_<T>::value();
            MPIWRAP_CALL(MPI_Compare_and_swap(&put, &compare, &get, type, rank, offset, win));
            return get;
        }

        /*
         * TODO: Rput, etc.
         */

#endif
};

}

#endif

#endif

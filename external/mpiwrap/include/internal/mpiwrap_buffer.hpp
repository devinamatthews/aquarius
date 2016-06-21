#ifndef _MPIWRAP_BUFFER_HPP_
#define _MPIWRAP_BUFFER_HPP_

#include "mpiwrap_common.hpp"

namespace MPIWrap
{

class Buffer
{
    protected:
        void* data;
        MPI_Int size;
        bool isAlloced;

    public:
        explicit Buffer(MPI_Int n)
        : data(::operator new(n)), size(n), isAlloced(true)
        {
            MPIWRAP_CALL(MPI_Buffer_attach(data, size));
        }

        Buffer(void* data, MPI_Int n)
        : data(data), size(n), isAlloced(false)
        {
            MPIWRAP_CALL(MPI_Buffer_attach(data, size));
        }

        template <typename Container>
        explicit Buffer(Container& c)
        : data(c.data()), size(c.size()), isAlloced(false)
        {
            MPIWRAP_CALL(MPI_Buffer_attach(data, size));
        }

        ~Buffer()
        {
            MPIWRAP_CALL(MPI_Buffer_detach(&data, &size));
            if (isAlloced) ::operator delete(data);
        }
};

}

#endif

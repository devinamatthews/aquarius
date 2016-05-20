#ifndef _MPIWRAP_STATUS_HPP_
#define _MPIWRAP_STATUS_HPP_

#include "mpiwrap_common.hpp"
#include "mpiwrap_datatype.hpp"

namespace MPIWrap
{

class Status
{
    friend class Intracomm;

    protected:
        MPI_Status status;

        explicit Status(const MPI_Status& status) : status(status) {}

    public:
        Status() {}

        operator MPI_Status&() { return status; }

        operator const MPI_Status&() const { return status; }

        operator MPI_Status*() { return &status; }

        operator const MPI_Status*() const { return &status; }

        MPI_Int source() const
        {
            return status.MPI_SOURCE;
        }

        MPI_Int tag() const
        {
            return status.MPI_TAG;
        }

        template <typename T>
        MPI_Int count() const
        {
            return count(MPI_TYPE_<T>::value());
        }

        MPI_Int count(const MPI_Datatype& type) const
        {
            MPI_Int n;
            MPIWRAP_CALL(MPI_Get_count(nc(&status), type, &n));
            return n;
        }

#if MPIWRAP_HAVE_MPI_COUNT

        template <typename T>
        MPI_Count elements() const
        {
            return elements(MPI_TYPE_<T>::value());
        }

        MPI_Count elements(const MPI_Datatype& type) const
        {
        	MPI_Count n;
            MPIWRAP_CALL(MPI_Get_elements_x(nc(&status), type, &n));
            return n;
        }

#else

        template <typename T>
        MPI_Int elements() const
        {
            return elements(MPI_TYPE_<T>::value());
        }

        MPI_Int elements(const MPI_Datatype& type) const
        {
            MPI_Int n;
            MPIWRAP_CALL(MPI_Get_elements(nc(&status), type, &n));
            return n;
        }

#endif

        bool cancelled() const
        {
            MPI_Int flag;
            MPIWRAP_CALL(MPI_Test_cancelled(nc(&status), &flag));
            return flag;
        }

        // TODO: pack and unpack (incl. external)

        // TODO: envelope, contents
};

}

#endif

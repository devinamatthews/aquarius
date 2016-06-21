#ifndef _MPIWRAP_MESSAGE_HPP_
#define _MPIWRAP_MESSAGE_HPP_

#if MPIWRAP_HAVE_MPI_MESSAGE

#include "mpiwrap_common.hpp"
#include "mpiwrap_datatype.hpp"

namespace MPIWrap
{

class Intercomm;
class Intracomm;
namespace internal
{
    template <typename Derived> class Comm;
}

class Message
{
    friend class internal::Comm<Intercomm>;
    friend class internal::Comm<Intracomm>;

    protected:
        MPI_Message msg;

        Message(const MPI_Message m = MPI_MESSAGE_NULL) : msg(m) {}

    public:
#if MPIWRAP_CXX11

        Message(Message&& other) : msg(other.msg)
        {
            other.msg = MPI_MESSAGE_NULL;
        }

#endif

        operator MPI_Message&() { return msg; }

        operator const MPI_Message&() const { return msg; }

        operator MPI_Message*() { return &msg; }

        operator const MPI_Message*() const { return &msg; }

        operator bool() const { return msg != MPI_MESSAGE_NULL; }

        /*
         * MPI_Mrecv
         */

        template <typename T>
        void Recv(T* buf, MPI_Int count) const
        {
            Recv(buf, count, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Recv(std::vector<T>& buf) const
        {
            Recv(buf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Recv(T* buf, MPI_Int count, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Mrecv(buf, count, type, msg, MPI_STATUS_IGNORE));
        }

        template <typename T>
        void Recv(std::vector<T>& buf, const Datatype& type) const
        {
            Recv(&buf.front(), buf.size(), type);
        }

        template <typename T>
        void Recv(T* buf, MPI_Int count, Status& status) const
        {
            Recv(buf, count, status, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Recv(std::vector<T>& buf, Status& status) const
        {
            Recv(buf, status, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Recv(T* buf, MPI_Int count, Status& status, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Mrecv(buf, count, type, msg, status));
        }

        template <typename T>
        void Recv(std::vector<T>& buf, Status& status, const Datatype& type) const
        {
            Recv(&buf.front(), buf.size(), status, type);
        }

        /*
         * MPI_Imrecv
         */

        template <typename T>
        Request Irecv(T* buf, MPI_Int count) const
        {
            return Irecv(buf, count, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Irecv(std::vector<T>& buf) const
        {
            return Irecv(buf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Irecv(T* buf, MPI_Int count, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Imrecv(buf, count, type, msg, &req));
            return Request(req);
        }

        template <typename T>
        Request Irecv(std::vector<T>& buf, const Datatype& type) const
        {
            return Irecv(&buf.front(), buf.size(), type);
        }
};

}

#endif

#endif

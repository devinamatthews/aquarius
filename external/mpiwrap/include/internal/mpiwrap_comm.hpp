#ifndef _MPIWRAP_COMM_HPP_
#define _MPIWRAP_COMM_HPP_

#include "mpiwrap_common.hpp"
#include "mpiwrap_group.hpp"
#include "mpiwrap_message.hpp"
#include "mpiwrap_window.hpp"
#include "mpiwrap_info.hpp"

#include <limits>

namespace MPIWrap
{

class Intercomm;
class Intracomm;

namespace internal
{

template <typename Derived>
class Comm
{
    friend class Intercomm;
    friend class Intracomm;

    protected:
        MPI_Comm comm;
        const MPI_Int npeers;

        static int getRank(const MPI_Comm& comm)
        {
            MPI_Int rank;
            MPIWRAP_CALL(MPI_Comm_rank(comm, &rank));
            return rank;
        }

        static int getSize(const MPI_Comm& comm)
        {
            MPI_Int size;
            MPIWRAP_CALL(MPI_Comm_size(comm, &size));
            return size;
        }

        std::vector<MPI_Int> displacements(const MPI_Int* counts) const
        {
            std::vector<MPI_Int> displs(size);
            for (int i = 1;i < size;i++)
            {
                if (displs[i-1] > std::numeric_limits<MPI_Int>::max()-counts[i-1])
                {
                    printf("displacement %d too large\n", i);
                    abort();
                }
                displs[i] = displs[i-1]+counts[i-1];
            }
            return displs;
        }

        std::vector<MPI_Int> displacements(const std::vector<MPI_Int>& counts) const
        {
            return displacements(&counts.front());
        }

        MPI_Int sum(const std::vector<MPI_Int>& counts) const
        {
            return std::accumulate(counts.begin(), counts.end(), 0);
        }

        explicit Comm(const MPI_Comm& comm, MPI_Int npeers)
        : comm(comm), npeers(npeers), rank(getRank(comm)), size(getSize(comm)) {}

    public:
        const MPI_Int rank;
        const MPI_Int size;

#if MPIWRAP_CXX11

        Comm(Comm&& other)
        : comm(other.comm), npeers(other.npeers), rank(other.rank), size(other.size)
        {
            other.comm = MPI_COMM_NULL;
        }

#endif

        ~Comm()
        {
            if (comm != MPI_COMM_NULL &&
                comm != MPI_COMM_WORLD &&
                comm != MPI_COMM_SELF)
            {
                MPIWRAP_CALL(MPI_Comm_free(&comm));
            }
        }

        operator MPI_Comm&() { return comm; }

        operator const MPI_Comm&() const { return comm; }

        operator MPI_Comm*() { return &comm; }

        operator const MPI_Comm*() const { return &comm; }

        Group group() const
        {
            MPI_Group g;
            MPIWRAP_CALL(MPI_Comm_group(comm, &g));
            return Group(g);
        }

        MPI_Int compare(const Comm& other) const
        {
            MPI_Int result;
            MPIWRAP_CALL(MPI_Comm_compare(comm, other, &result));
            return result;
        }

        Derived duplicate() const
        {
            MPI_Comm c;
            MPIWRAP_CALL(MPI_Comm_dup(comm, &c));
            return Derived(c);
        }

        Derived subset(const Group& group) const
        {
            MPI_Comm c;
            MPIWRAP_CALL(MPI_Comm_create(comm, group, &c));
            return Derived(c);
        }

        Derived split(MPI_Int color, MPI_Int key)
        {
            MPI_Comm c;
            MPIWRAP_CALL(MPI_Comm_split(comm, color, key, &c));
            return Derived(c);
        }

        bool isIntercommunicator() const
        {
            MPI_Int flag;
            MPIWRAP_CALL(MPI_Comm_test_inter(comm, &flag));
            return flag;
        }

        bool isIntracommunicator() const
        {
            return !isIntercommunicator();
        }

#if MPIWRAP_HAVE_MPI_WIN

        Window window(void* data, MPI_Aint size)
        {
            return window(data, size, Info::null());
        }

        Window window(void* data, MPI_Aint size, const Info& info)
        {
            MPI_Win w;
            MPIWRAP_CALL(MPI_Win_create(data, size, 1, info, comm, &w));
            return Window(w);
        }

        template <typename Container>
        Window window(Container& c)
        {
            return window(&c[0], c.size());
        }

        template <typename Container>
        Window window(Container& c, const Info& info)
        {
            return window(&c[0], c.size(), info);
        }

        Window window(MPI_Aint size)
        {
            return window(size, Info::null());
        }

        Window window(MPI_Aint size, const Info& info)
        {
            MPI_Win w;
            void* foo;
            MPIWRAP_CALL(MPI_Win_allocate(size, 1, info, comm, &foo, &w));
            return Window(w);
        }

#endif

        //TODO: attributes

#if MPIWRAP_HAVE_MPI_MESSAGE

        /*
         * MPI_Mprobe
         */

        Message Mprobe(MPI_Int source, MPI_Int tag, Status& status) const
        {
            Message m;
            MPIWRAP_CALL(MPI_Mprobe(source, tag, comm, m, status));
            return m;
        }

        Message Mprobe(MPI_Int source, MPI_Int tag) const
        {
            Message m;
            MPIWRAP_CALL(MPI_Mprobe(source, tag, comm, m, MPI_STATUS_IGNORE));
            return m;
        }

#endif

        /*
         * MPI_Probe
         */

        void Probe(MPI_Int source, MPI_Int tag, Status& status) const
        {
            MPIWRAP_CALL(MPI_Probe(source, tag, comm, status));
        }

        void Probe(MPI_Int source, MPI_Int tag) const
        {
            MPIWRAP_CALL(MPI_Probe(source, tag, comm, MPI_STATUS_IGNORE));
        }

        /*
         * MPI_Recv
         */

        template <typename T>
        void Recv(T* buf, MPI_Int count, MPI_Int source, MPI_Int tag) const
        {
            Recv(buf, count, source, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Recv(std::vector<T>& buf, MPI_Int source, MPI_Int tag) const
        {
            Recv(buf, source, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Recv(T* buf, MPI_Int count, MPI_Int source, MPI_Int tag, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Recv(buf, count, type, source, tag, comm, MPI_STATUS_IGNORE));
        }

        template <typename T>
        void Recv(std::vector<T>& buf, MPI_Int source, MPI_Int tag, const Datatype& type) const
        {
            Recv(&buf.front(), buf.size(), source, tag, type);
        }

        template <typename T>
        void Recv(T* buf, MPI_Int count, MPI_Int source, MPI_Int tag, Status& status) const
        {
            Recv(buf, count, source, tag, status, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Recv(std::vector<T>& buf, MPI_Int source, MPI_Int tag, Status& status) const
        {
            Recv(buf, source, tag, status, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Recv(T* buf, MPI_Int count, MPI_Int source, MPI_Int tag, Status& status, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Recv(buf, count, type, source, tag, comm, status));
        }

        template <typename T>
        void Recv(std::vector<T>& buf, MPI_Int source, MPI_Int tag, Status& status, const Datatype& type) const
        {
            Recv(&buf.front(), buf.size(), source, tag, status, type);
        }

        /*
         * MPI_Recv_init
         */

        template <typename T>
        Request Recv_init(T* buf, MPI_Int count, MPI_Int source, MPI_Int tag) const
        {
            return Recv_init(buf, count, source, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Recv_init(std::vector<T>& buf, MPI_Int source, MPI_Int tag) const
        {
            return Recv_init(buf, source, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Recv_init(T* buf, MPI_Int count, MPI_Int source, MPI_Int tag, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Recv_init(buf, count, type, source, tag, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Recv_init(std::vector<T>& buf, MPI_Int source, MPI_Int tag, const Datatype& type) const
        {
            return Recv_init(&buf.front(), buf.size(), source, tag, type);
        }

        /*
         * MPI_Bsend
         */

        template <typename T>
        void Bsend(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag) const
        {
            Bsend(buf, count, dest, tag);
        }

        template <typename T>
        void Bsend(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag) const
        {
            Bsend(buf, dest, tag);
        }

        template <typename T>
        void Bsend(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Bsend(buf, count, type, dest, tag, comm));
        }

        template <typename T>
        void Bsend(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            Bsend(&buf.front(), buf.size(), dest, tag, type);
        }

        /*
         * MPI_Bsend_init
         */

        template <typename T>
        Request Bsend_init(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag) const
        {
            return Bsend_init(buf, count, dest, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Bsend_init(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag) const
        {
            return Bsend_init(buf, dest, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Bsend_init(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Bsend_init(buf, count, type, dest, tag, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Bsend_init(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            return Bsend_init(&buf.front(), buf.size(), dest, tag, type);
        }

        /*
         * MPI_Rsend
         */

        template <typename T>
        void Rsend(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag) const
        {
            Rsend(buf, count, dest, tag);
        }

        template <typename T>
        void Rsend(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag) const
        {
            Rsend(buf, dest, tag);
        }

        template <typename T>
        void Rsend(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Rsend(buf, count, type, dest, tag, comm));
        }

        template <typename T>
        void Rsend(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            Rsend(&buf.front(), buf.size(), dest, tag, type);
        }

        /*
         * MPI_Rsend_init
         */

        template <typename T>
        Request Rsend_init(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag) const
        {
            return Rsend_init(buf, count, dest, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Rsend_init(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag) const
        {
            return Rsend_init(buf, dest, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Rsend_init(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Rsend_init(buf, count, type, dest, tag, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Rsend_init(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            return Rsend_init(&buf.front(), buf.size(), dest, tag, type);
        }

        /*
         * MPI_Send
         */

        template <typename T>
        void Send(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag) const
        {
            Send(buf, count, dest, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Send(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag) const
        {
            Send(buf, dest, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Send(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Send(buf, count, type, dest, tag, comm));
        }

        template <typename T>
        void Send(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            Send(&buf.front(), buf.size(), dest, tag, type);
        }

        /*
         * MPI_Send_init
         */

        template <typename T>
        Request Send_init(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag) const
        {
            return Send_init(buf, count, dest, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Send_init(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag) const
        {
            return Send_init(buf, dest, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Send_init(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Send_init(buf, count, type, dest, tag, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Send_init(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            return Send_init(&buf.front(), buf.size(), dest, tag, type);
        }

        /*
         * MPI_Sendrecv
         */

        template <typename T>
        void Sendrecv(const T* sendbuf, MPI_Int sendcount, MPI_Int dest,
                            T* recvbuf, MPI_Int recvcount, MPI_Int source, MPI_Int tag) const
        {
            Sendrecv(sendbuf, sendcount, dest,
                     recvbuf, recvcount, source, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Sendrecv(const std::vector<T>& sendbuf, MPI_Int dest,
                            std::vector<T>& recvbuf, MPI_Int source, MPI_Int tag) const
        {
            Sendrecv(sendbuf, dest,
                     recvbuf, source, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Sendrecv(const T* sendbuf, MPI_Int sendcount, MPI_Int dest,
                            T* recvbuf, MPI_Int recvcount, MPI_Int source, MPI_Int tag,
                      const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Sendrecv(sendbuf, sendcount, type,   dest, tag,
                                      recvbuf, recvcount, type, source, tag, comm, MPI_STATUS_IGNORE));
        }

        template <typename T>
        void Sendrecv(const std::vector<T>& sendbuf, MPI_Int dest,
                            std::vector<T>& recvbuf, MPI_Int source, MPI_Int tag,
                      const Datatype& type) const
        {
            Sendrecv(&sendbuf.front(), sendbuf.size(), dest,
                     &recvbuf.front(), recvbuf.size(), source, tag, type);
        }

        template <typename T>
        void Sendrecv(const T* sendbuf, MPI_Int sendcount, MPI_Int dest,
                            T* recvbuf, MPI_Int recvcount, MPI_Int source, MPI_Int tag, Status& status) const
        {
            Sendrecv(sendbuf, sendcount, dest,
                     recvbuf, recvcount, source, tag, status, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Sendrecv(const std::vector<T>& sendbuf, MPI_Int dest,
                            std::vector<T>& recvbuf, MPI_Int source, MPI_Int tag, Status& status) const
        {
            Sendrecv(sendbuf, dest,
                     recvbuf, source, tag, status, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Sendrecv(const T* sendbuf, MPI_Int sendcount, MPI_Int dest,
                            T* recvbuf, MPI_Int recvcount, MPI_Int source, MPI_Int tag, Status& status,
                      const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Sendrecv(sendbuf, sendcount, type,   dest, tag,
                                      recvbuf, recvcount, type, source, tag, comm, status));
        }

        template <typename T>
        void Sendrecv(const std::vector<T>& sendbuf, MPI_Int dest,
                            std::vector<T>& recvbuf, MPI_Int source, MPI_Int tag, Status& status,
                      const Datatype& type) const
        {
            Sendrecv(&sendbuf.front(), sendbuf.size(), dest,
                     &recvbuf.front(), recvbuf.size(), source, tag, status, type);
        }

        /*
         * MPI_Sendrecv_replace
         */

        template <typename T>
        void Sendrecv(T* buf, MPI_Int count, MPI_Int dest, MPI_Int source, MPI_Int tag) const
        {
            Sendrecv(buf, count, dest, source, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Sendrecv(std::vector<T>& buf, MPI_Int dest, MPI_Int source, MPI_Int tag) const
        {
            Sendrecv(buf, dest, source, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Sendrecv(T* buf, MPI_Int count, MPI_Int dest, MPI_Int source, MPI_Int tag, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Sendrecv_replace(buf, count, type, dest, tag, source, tag, comm, MPI_STATUS_IGNORE));
        }

        template <typename T>
        void Sendrecv(std::vector<T>& buf, MPI_Int dest, MPI_Int source, MPI_Int tag, const Datatype& type) const
        {
            Sendrecv(&buf.front(), buf.size(), dest, source, tag, type);
        }

        template <typename T>
        void Sendrecv(T* buf, MPI_Int count, MPI_Int dest, MPI_Int source, MPI_Int tag, Status& status) const
        {
            Sendrecv(buf, count, dest, source, tag, status, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Sendrecv(std::vector<T>& buf, MPI_Int dest, MPI_Int source, MPI_Int tag, Status& status) const
        {
            Sendrecv(buf, dest, source, tag, status, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Sendrecv(T* buf, MPI_Int count, MPI_Int dest, MPI_Int source, MPI_Int tag, Status& status,
                      const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Sendrecv_replace(buf, count, type, dest, tag, source, tag, comm, status));
        }

        template <typename T>
        void Sendrecv(std::vector<T>& buf, MPI_Int dest, MPI_Int source, MPI_Int tag, Status& status,
                      const Datatype& type) const
        {
            Sendrecv(&buf.front(), buf.size(), dest, source, tag, status, type);
        }

        /*
         * MPI_Ssend
         */

        template <typename T>
        void Ssend(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag) const
        {
            Ssend(buf, count, dest, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Ssend(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag) const
        {
            Ssend(buf, dest, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Ssend(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Ssend(buf, count, type, dest, tag, comm));
        }

        template <typename T>
        void Ssend(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            Ssend(&buf.front(). buf.size(), dest, tag, type);
        }

        /*
         * MPI_Ssend_init
         */

        template <typename T>
        Request Ssend_init(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag) const
        {
            return Ssend_init(buf, count, dest, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ssend_init(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag) const
        {
            return Ssend_init(buf, dest, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ssend_init(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Ssend_init(buf, count, type, dest, tag, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Ssend_init(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            return Ssend_init(&buf.front(), buf.size(), dest, tag, type);
        }


#if MPIWRAP_HAVE_MPI_MESSAGE

        /*
         * MPI_Improbe
         */

        Message Improbe(MPI_Int source, MPI_Int tag, Status& status) const
        {
            Message m;
            MPI_Int flag;
            MPIWRAP_CALL(MPI_Improbe(source, tag, comm, &flag, m, status));
            if (flag) return m;
            else return Message();
        }

        Message Improbe(MPI_Int source, MPI_Int tag) const
        {
            Message m;
            MPI_Int flag;
            MPIWRAP_CALL(MPI_Improbe(source, tag, comm, &flag, m, MPI_STATUS_IGNORE));
            if (flag) return m;
            else return Message();
        }

#endif

        /*
         * MPI_Iprobe
         */

        bool Iprobe(MPI_Int source, MPI_Int tag, Status& status) const
        {
            MPI_Int flag;
            MPIWRAP_CALL(MPI_Iprobe(source, tag, comm, &flag, status));
            return flag;
        }

        bool Iprobe(MPI_Int source, MPI_Int tag) const
        {
            MPI_Int flag;
            MPIWRAP_CALL(MPI_Iprobe(source, tag, comm, &flag, MPI_STATUS_IGNORE));
            return flag;
        }

        /*
         * MPI_Irecv
         */

        template <typename T>
        Request Irecv(T* buf, MPI_Int count, MPI_Int source, MPI_Int tag) const
        {
            return Irecv(buf, count, source, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Irecv(std::vector<T>& buf, MPI_Int source, MPI_Int tag) const
        {
            return Irecv(buf, source, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Irecv(T* buf, MPI_Int count, MPI_Int source, MPI_Int tag, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Irecv(buf, count, type, source, tag, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Irecv(std::vector<T>& buf, MPI_Int source, MPI_Int tag, const Datatype& type) const
        {
            return Irecv(&buf.front(), buf.size(), source, tag, type);
        }

        /*
         * MPI_Ibsend
         */

        template <typename T>
        Request Ibsend(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag) const
        {
            return Ibsend(buf, count, dest, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ibsend(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag) const
        {
            return Ibsend(buf, dest, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ibsend(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Ibsend(buf, count, type, dest, tag, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Ibsend(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            return Ibsend(&buf.front(), buf.size(), dest, tag, type);
        }

        /*
         * MPI_Irsend
         */

        template <typename T>
        Request Irsend(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag) const
        {
            return Irsend(buf, count, dest, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Irsend(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag) const
        {
            return Irsend(buf, dest, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Irsend(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Irsend(buf, count, type, dest, tag, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Irsend(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            return Irsend(&buf.front(), buf.size(), dest, tag, type);
        }

        /*
         * MPI_Isend
         */

        template <typename T>
        Request Isend(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag) const
        {
            return Isend(buf, count, dest, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Isend(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag) const
        {
            return Isend(buf, dest, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Isend(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Isend(buf, count, type, dest, tag, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Isend(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            return Isend(&buf.front(), buf.size(), dest, tag, type);
        }

        /*
         * MPI_Issend
         */

        template <typename T>
        Request Issend(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag) const
        {
            return Issend(buf, count, dest, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Issend(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag) const
        {
            return Issend(buf, dest, tag, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Issend(const T* buf, MPI_Int count, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Issend(buf, count, type, dest, tag, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Issend(const std::vector<T>& buf, MPI_Int dest, MPI_Int tag, const Datatype& type) const
        {
            return Issend(&buf.front(), buf.size(), dest, tag, type);
        }

        /*
         * MPI_Allgather
         */

        template <typename T>
        void Allgather(const T* sendbuf, T* recvbuf, MPI_Int count) const
        {
            Allgather(sendbuf, recvbuf, count, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Allgather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf) const
        {
            Allgather(sendbuf, recvbuf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Allgather(const T* sendbuf, T* recvbuf, MPI_Int count, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Allgather(nc(sendbuf), count, type, recvbuf, count, type, comm));
        }

        template <typename T>
        void Allgather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendbuf.size() == recvbuf.size(),
                           "Send and receive buffers must be the same size.");
            Allgather(&sendbuf.front(), &recvbuf.front(), sendbuf.size(), type);
        }

        /*
         * MPI_Allgatherv
         */

        template <typename T>
        void Allgather(const T* sendbuf, MPI_Int sendcount, T* recvbuf, const MPI_Int* recvcounts) const
        {
            Allgather(sendbuf, sendcount, recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Allgather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf,
                       const std::vector<MPI_Int>& recvcounts) const
        {
            Allgather(sendbuf, recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Allgather(const T* sendbuf, MPI_Int sendcount, T* recvbuf, const MPI_Int* recvcounts,
                       const Datatype& type) const
        {
            std::vector<MPI_Int> displs = displacements(recvcounts);
            Allgather(sendbuf, sendcount, recvbuf, recvcounts, displs, type);
        }

        template <typename T>
        void Allgather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts,
                       const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvbuf.size() == sum(recvcounts),
                           "The receive buffer size must equal the sum of the receive counts.");
            Allgather(sendbuf, recvbuf, recvcounts, displacements(recvcounts), type);
        }

        template <typename T>
        void Allgather(const T* sendbuf, MPI_Int sendcount, T* recvbuf, const MPI_Int* recvcounts, const MPI_Int* recvdispls) const
        {
            Allgather(sendbuf, sendcount, recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Allgather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf,
                       const std::vector<MPI_Int>& recvcounts, const std::vector<MPI_Int>& recvdispls) const
        {
            Allgather(sendbuf, recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Allgather(const T* sendbuf, MPI_Int sendcount, T* recvbuf, const MPI_Int* recvcounts,
                       const MPI_Int* recvdispls, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Allgatherv(nc(sendbuf), sendcount, type, recvbuf, nc(recvcounts), nc(recvdispls), type, comm));
        }

        template <typename T>
        void Allgather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts,
                       const std::vector<MPI_Int>& recvdispls, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvcounts.size() == npeers,
                           "There must be exactly one receive count for each process.");
            MPIWRAP_ASSERT(recvcounts.size() == npeers,
                           "There must be exactly one receive displacement for each process.");
            MPIWRAP_ASSERT(recvbuf.size() >= recvdispls[npeers-1]+recvcounts[npeers-1],
                           "The receive buffer size must be at least as large the sum of last receive count and the last receive displacement.");
            Allgather(&sendbuf.front(), sendbuf.size(), &recvbuf.front(), &recvcounts.front(), &recvdispls.front(), type);
        }

        /*
         * MPI_Allreduce
         */

        template <typename T>
        void Allreduce(const T* sendbuf, T* recvbuf, MPI_Int count, const MPI_Op& op) const
        {
            Allreduce(sendbuf, recvbuf, count, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Allreduce(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI_Op& op) const
        {
            Allreduce(sendbuf, recvbuf, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Allreduce(const T* sendbuf, T* recvbuf, MPI_Int count, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Allreduce(nc(sendbuf), recvbuf, count, type, op, comm));
        }

        template <typename T>
        void Allreduce(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI_Op& op,
                       const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendbuf.size() == recvbuf.size(),
                           "Send and receive buffers must be the same size.");
            Allreduce(&sendbuf.front(), &recvbuf.front(), sendbuf.size(), op, type);
        }

        /*
         * MPI_Alltoall
         */

        template <typename T>
        void Alltoall(const T* sendbuf, T* recvbuf, MPI_Int count) const
        {
            Alltoall(sendbuf, recvbuf, count, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Alltoall(const std::vector<T>& sendbuf, std::vector<T>& recvbuf) const
        {
            Alltoall(sendbuf, recvbuf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Alltoall(const T* sendbuf, T* recvbuf, MPI_Int count, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Alltoall(nc(sendbuf), count, type, recvbuf, count, type, comm));
        }

        template <typename T>
        void Alltoall(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendbuf.size() == recvbuf.size(),
                           "Send and receive buffers must be the same size.");
            MPIWRAP_ASSERT(recvbuf.size()%size == 0,
                           "Receive buffer size must be a multiple of the communicator size.");
            Alltoall(&sendbuf.front(), &recvbuf.front(), recvbuf.size()/size, type);
        }

        /*
         * MPI_Alltoallv
         */

        template <typename T>
        void Alltoall(const T* sendbuf, const MPI_Int* sendcounts, T* recvbuf, const MPI_Int* recvcounts) const
        {
            Alltoall(sendbuf, sendcounts, recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Alltoall(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts,
                            std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts) const
        {
            Alltoall(sendbuf, sendcounts, recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Alltoall(const T* sendbuf, const MPI_Int* sendcounts,
                            T* recvbuf, const MPI_Int* recvcounts, const Datatype& type) const
        {
            std::vector<MPI_Int> senddispls = displacements(sendcounts);
            std::vector<MPI_Int> recvdispls = displacements(recvcounts);
            Alltoall(sendbuf, sendcounts, senddispls,
                     recvbuf, recvcounts, recvdispls, type);
        }

        template <typename T>
        void Alltoall(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts,
                            std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts,
                      const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendbuf.size() == sum(sendcounts),
                           "The send buffer size must equal the sum of the send counts.");
            MPIWRAP_ASSERT(recvbuf.size() == sum(recvcounts),
                           "The receive buffer size must equal the sum of the receive counts.");
            Alltoall(sendbuf, sendcounts, displacements(sendcounts),
                     recvbuf, recvcounts, displacements(recvcounts), type);
        }

        template <typename T>
        void Alltoall(const T* sendbuf, const MPI_Int* sendcounts, const MPI_Int* senddispls,
                            T* recvbuf, const MPI_Int* recvcounts, const MPI_Int* recvdispls) const
        {
            Alltoall(sendbuf, sendcounts, senddispls,
                     recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Alltoall(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts, const std::vector<MPI_Int>& senddispls,
                            std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts, const std::vector<MPI_Int>& recvdispls) const
        {
            Alltoall(sendbuf, sendcounts, senddispls,
                     recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Alltoall(const T* sendbuf, const MPI_Int* sendcounts, const MPI_Int* senddispls,
                            T* recvbuf, const MPI_Int* recvcounts, const MPI_Int* recvdispls, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Alltoallv(nc(sendbuf), nc(sendcounts), nc(senddispls), type,
                                          recvbuf , nc(recvcounts), nc(recvdispls), type, comm));
        }

        template <typename T>
        void Alltoall(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts, const std::vector<MPI_Int>& senddispls,
                            std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts, const std::vector<MPI_Int>& recvdispls,
                      const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendcounts.size() == npeers,
                           "There must be exactly one send count for each process.");
            MPIWRAP_ASSERT(recvcounts.size() == npeers,
                           "There must be exactly one receive count for each process.");
            MPIWRAP_ASSERT(senddispls.size() == npeers,
                           "There must be exactly one send displacement for each process.");
            MPIWRAP_ASSERT(recvdispls.size() == npeers,
                           "There must be exactly one receive displacement for each process.");
            MPIWRAP_ASSERT(sendbuf.size() >= senddispls[npeers-1]+sendcounts[npeers-1],
                           "The send buffer size must be at least as large the sum of last send count and the last send displacement.");
            MPIWRAP_ASSERT(recvbuf.size() >= recvdispls[npeers-1]+recvcounts[npeers-1],
                           "The receive buffer size must be at least as large the sum of last receive count and the last receive displacement.");
            Alltoall(&sendbuf.front(), &sendcounts.front(), &senddispls.front(),
                     &recvbuf.front(), &recvcounts.front(), &recvdispls.front(), type);
        }

        /*
         * MPI_Alltoallw
         */

        //TODO

        /*
         * MPI_Barrier
         */

        void Barrier() const
        {
            MPIWRAP_CALL(MPI_Barrier(comm));
        }

        /*
         * MPI_Reduce_scatter
         */

        template <typename T>
        void Reduce_scatter(const T* sendbuf, T* recvbuf, const MPI_Int* recvcounts, const MPI_Op& op) const
        {
            Reduce_scatter(sendbuf, recvbuf, recvcounts, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Reduce_scatter(const std::vector<T>& sendbuf, std::vector<T>& recvbuf,
                            const std::vector<MPI_Int>& recvcounts, const MPI_Op& op) const
        {
            Reduce_scatter(sendbuf, recvbuf, recvcounts, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Reduce_scatter(const T* sendbuf, T* recvbuf, const MPI_Int* recvcounts, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Reduce_scatter(nc(sendbuf), recvbuf, nc(recvcounts), type, op, comm));
        }

        template <typename T>
        void Reduce_scatter(const std::vector<T>& sendbuf, std::vector<T>& recvbuf,
                            const std::vector<MPI_Int>& recvcounts, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvcounts.size() == npeers,
                           "There must be exactly one receive count for each process.");
            MPIWRAP_ASSERT(sendbuf.size() == sum(recvcounts),
                           "The send buffer size must equal the sum of the receive counts.");
            if (!isIntercommunicator())
            MPIWRAP_ASSERT(recvbuf.size() == recvcounts[rank],
                           "The receive buffer size must equal the receive count for this process.");
            Reduce_scatter(&sendbuf.front(), &recvbuf.front(), &recvcounts.front(), op, type);
        }

        /*
         * MPI_Reduce_scatter_block
         */

        template <typename T>
        void Reduce_scatter(const T* sendbuf, T* recvbuf, MPI_Int count, const MPI_Op& op) const
        {
            Reduce_scatter(sendbuf, recvbuf, count, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Reduce_scatter(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI_Op& op) const
        {
            Reduce_scatter(sendbuf, recvbuf, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Reduce_scatter(const T* sendbuf, T* recvbuf, MPI_Int count, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Reduce_scatter_block(nc(sendbuf), recvbuf, count, type, op, comm));
        }

        template <typename T>
        void Reduce_scatter(const std::vector<T>& sendbuf, std::vector<T>& recvbuf,
                            const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendbuf.size() == recvbuf.size()*npeers,
                           "The send buffer must equal the receive buffer size times the number of processes.");
            Reduce_scatter(&sendbuf.front(), &recvbuf.front(), recvbuf.size(), op, type);
        }

#if MPIWRAP_HAVE_MPI_ICOLLECTIVES

        /*
         * MPI_Iallgather
         */

        template <typename T>
        Request Iallgather(const T* sendbuf, T* recvbuf, MPI_Int count) const
        {
            return Iallgather(sendbuf, recvbuf, count, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iallgather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf) const
        {
            return Iallgather(sendbuf, recvbuf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iallgather(const T* sendbuf, T* recvbuf, MPI_Int count, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Iallgather(nc(sendbuf), count, type, recvbuf, count, type, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Iallgather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendbuf.size() == recvbuf.size(),
                           "Send and receive buffers must be the same size.");
            return Iallgather(&sendbuf.front(), &recvbuf.front(), sendbuf.size(), type);
        }

        /*
         * MPI_Iallgatherv
         */

        template <typename T>
        Request Iallgather(const T* sendbuf, MPI_Int sendcount, T* recvbuf, const MPI_Int* recvcounts) const
        {
            return Iallgather(sendbuf, sendcount, recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iallgather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf,
                       const std::vector<MPI_Int>& recvcounts) const
        {
            return Iallgather(sendbuf, recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iallgather(const T* sendbuf, MPI_Int sendcount, T* recvbuf, const MPI_Int* recvcounts,
                       const Datatype& type) const
        {
            std::vector<MPI_Int> displs = displacements(recvcounts);
            return Iallgather(sendbuf, sendcount, recvbuf, recvcounts, displs, type);
        }

        template <typename T>
        Request Iallgather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts,
                       const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvbuf.size() == sum(recvcounts),
                           "The receive buffer size must equal the sum of the receive counts.");
            return Iallgather(sendbuf, recvbuf, recvcounts, displacements(recvcounts), type);
        }

        template <typename T>
        Request Iallgather(const T* sendbuf, MPI_Int sendcount, T* recvbuf, const MPI_Int* recvcounts, const MPI_Int* recvdispls) const
        {
            return Iallgather(sendbuf, sendcount, recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iallgather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf,
                       const std::vector<MPI_Int>& recvcounts, const std::vector<MPI_Int>& recvdispls) const
        {
            return Iallgather(sendbuf, recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iallgather(const T* sendbuf, MPI_Int sendcount, T* recvbuf, const MPI_Int* recvcounts,
                       const MPI_Int* recvdispls, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Iallgatherv(nc(sendbuf), sendcount, type, recvbuf, nc(recvcounts), nc(recvdispls), type, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Iallgather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts,
                       const std::vector<MPI_Int>& recvdispls, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvcounts.size() == npeers,
                           "There must be exactly one receive count for each process.");
            MPIWRAP_ASSERT(recvcounts.size() == npeers,
                           "There must be exactly one receive displacement for each process.");
            MPIWRAP_ASSERT(recvbuf.size() >= recvdispls[npeers-1]+recvcounts[npeers-1],
                           "The receive buffer size must be at least as large the sum of last receive count and the last receive displacement.");
            return Iallgather(&sendbuf.front(), sendbuf.size(), &recvbuf.front(), &recvcounts.front(), &recvdispls.front(), type);
        }

        /*
         * MPI_Iallreduce
         */

        template <typename T>
        Request Iallreduce(const T* sendbuf, T* recvbuf, MPI_Int count, const MPI_Op& op) const
        {
            return Iallreduce(sendbuf, recvbuf, count, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iallreduce(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI_Op& op) const
        {
            return Iallreduce(sendbuf, recvbuf, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iallreduce(const T* sendbuf, T* recvbuf, MPI_Int count, const MPI_Op& op, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Iallreduce(nc(sendbuf), recvbuf, count, type, op, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Iallreduce(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI_Op& op,
                       const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendbuf.size() == recvbuf.size(),
                           "Send and receive buffers must be the same size.");
            return Iallreduce(&sendbuf.front(), &recvbuf.front(), sendbuf.size(), op, type);
        }

        /*
         * MPI_Ialltoall
         */

        template <typename T>
        Request Ialltoall(const T* sendbuf, T* recvbuf, MPI_Int count) const
        {
            return Ialltoall(sendbuf, recvbuf, count, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ialltoall(const std::vector<T>& sendbuf, std::vector<T>& recvbuf) const
        {
            return Ialltoall(sendbuf, recvbuf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ialltoall(const T* sendbuf, T* recvbuf, MPI_Int count, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Ialltoall(nc(sendbuf), count, type, recvbuf, count, type, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Ialltoall(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendbuf.size() == recvbuf.size(),
                           "Send and receive buffers must be the same size.");
            MPIWRAP_ASSERT(recvbuf.size()%size == 0,
                           "Receive buffer size must be a multiple of the communicator size.");
            return Ialltoall(&sendbuf.front(), &recvbuf.front(), recvbuf.size()/size, type);
        }

        /*
         * MPI_Ialltoallv
         */

        template <typename T>
        Request Ialltoall(const T* sendbuf, const MPI_Int* sendcounts, T* recvbuf, const MPI_Int* recvcounts) const
        {
            return Ialltoall(sendbuf, sendcounts, recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ialltoall(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts,
                                std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts) const
        {
            return Ialltoall(sendbuf, sendcounts, recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ialltoall(const T* sendbuf, const MPI_Int* sendcounts,
                                T* recvbuf, const MPI_Int* recvcounts, const Datatype& type) const
        {
            std::vector<MPI_Int> senddispls = displacements(sendcounts);
            std::vector<MPI_Int> recvdispls = displacements(recvcounts);
            return Ialltoall(sendbuf, sendcounts, senddispls,
                             recvbuf, recvcounts, recvdispls, type);
        }

        template <typename T>
        Request Ialltoall(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts,
                                std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts,
                          const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendbuf.size() == sum(sendcounts),
                           "The send buffer size must equal the sum of the send counts.");
            MPIWRAP_ASSERT(recvbuf.size() == sum(recvcounts),
                           "The receive buffer size must equal the sum of the receive counts.");
            return Ialltoall(sendbuf, sendcounts, displacements(sendcounts),
                             recvbuf, recvcounts, displacements(recvcounts), type);
        }

        template <typename T>
        Request Ialltoall(const T* sendbuf, const MPI_Int* sendcounts, const MPI_Int* senddispls,
                                T* recvbuf, const MPI_Int* recvcounts, const MPI_Int* recvdispls) const
        {
            return Ialltoall(sendbuf, sendcounts, senddispls,
                             recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ialltoall(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts, const std::vector<MPI_Int>& senddispls,
                                std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts, const std::vector<MPI_Int>& recvdispls) const
        {
            return Ialltoall(sendbuf, sendcounts, senddispls,
                             recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ialltoall(const T* sendbuf, const MPI_Int* sendcounts, const MPI_Int* senddispls,
                                T* recvbuf, const MPI_Int* recvcounts, const MPI_Int* recvdispls, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Ialltoallv(nc(sendbuf), nc(sendcounts), nc(senddispls), type,
                                           recvbuf , nc(recvcounts), nc(recvdispls), type, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Ialltoall(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts, const std::vector<MPI_Int>& senddispls,
                                std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts, const std::vector<MPI_Int>& recvdispls,
                          const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendcounts.size() == npeers,
                           "There must be exactly one send count for each process.");
            MPIWRAP_ASSERT(recvcounts.size() == npeers,
                           "There must be exactly one receive count for each process.");
            MPIWRAP_ASSERT(senddispls.size() == npeers,
                           "There must be exactly one send displacement for each process.");
            MPIWRAP_ASSERT(recvdispls.size() == npeers,
                           "There must be exactly one receive displacement for each process.");
            MPIWRAP_ASSERT(sendbuf.size() >= senddispls[npeers-1]+sendcounts[npeers-1],
                           "The send buffer size must be at least as large the sum of last send count and the last send displacement.");
            MPIWRAP_ASSERT(recvbuf.size() >= recvdispls[npeers-1]+recvcounts[npeers-1],
                           "The receive buffer size must be at least as large the sum of last receive count and the last receive displacement.");
            return Ialltoall(&sendbuf.front(), &sendcounts.front(), &senddispls.front(),
                             &recvbuf.front(), &recvcounts.front(), &recvdispls.front(), type);
        }

        /*
         * MPI_Ialltoallw
         */

        //TODO

        /*
         * MPI_Ibarrier
         */

        Request Ibarrier() const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Ibarrier(comm, &req));
            return Request(req);
        }

        /*
         * MPI_Ireduce_scatter
         */

        template <typename T>
        Request Ireduce_scatter(const T* sendbuf, T* recvbuf, const MPI_Int* recvcounts, const MPI_Op& op) const
        {
            return Ireduce_scatter(sendbuf, recvbuf, recvcounts, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ireduce_scatter(const std::vector<T>& sendbuf, std::vector<T>& recvbuf,
                                const std::vector<MPI_Int>& recvcounts, const MPI_Op& op) const
        {
            return Ireduce_scatter(sendbuf, recvbuf, recvcounts, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ireduce_scatter(const T* sendbuf, T* recvbuf, const MPI_Int* recvcounts, const MPI_Op& op, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Ireduce_scatter(nc(sendbuf), recvbuf, nc(recvcounts), type, op, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Ireduce_scatter(const std::vector<T>& sendbuf, std::vector<T>& recvbuf,
                                const std::vector<MPI_Int>& recvcounts, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvcounts.size() == npeers,
                           "There must be exactly one receive count for each process.");
            MPIWRAP_ASSERT(sendbuf.size() == sum(recvcounts),
                           "The send buffer size must equal the sum of the receive counts.");
            if (!isIntercommunicator())
            MPIWRAP_ASSERT(recvbuf.size() == recvcounts[rank],
                           "The receive buffer size must equal the receive count for this process.");
            return Ireduce_scatter(&sendbuf.front(), &recvbuf.front(), &recvcounts.front(), op, type);
        }

        /*
         * MPI_Ireduce_scatter_block
         */

        template <typename T>
        Request Ireduce_scatter(const T* sendbuf, T* recvbuf, MPI_Int count, const MPI_Op& op) const
        {
            return Ireduce_scatter(sendbuf, recvbuf, count, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ireduce_scatter(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI_Op& op) const
        {
            return Ireduce_scatter(sendbuf, recvbuf, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ireduce_scatter(const T* sendbuf, T* recvbuf, MPI_Int count, const MPI_Op& op, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Ireduce_scatter_block(nc(sendbuf), recvbuf, count, type, op, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Ireduce_scatter(const std::vector<T>& sendbuf, std::vector<T>& recvbuf,
                            const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendbuf.size() == recvbuf.size()*npeers,
                           "The send buffer must equal the receive buffer size times the number of processes.");
            return Ireduce_scatter(&sendbuf.front(), &recvbuf.front(), recvbuf.size(), op, type);
        }

#endif
};

}
}

#endif

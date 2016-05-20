#ifndef _MPIWRAP_INTERCOMM_HPP_
#define _MPIWRAP_INTERCOMM_HPP_

#include "mpiwrap_common.hpp"
#include "mpiwrap_datatype.hpp"
#include "mpiwrap_status.hpp"
#include "mpiwrap_request.hpp"
#include "mpiwrap_comm.hpp"
#include "mpiwrap_group.hpp"
#include "mpiwrap_intracomm.hpp"

namespace MPIWrap
{

class Intercomm : protected internal::Comm<Intercomm>
{
    friend class Intracomm;
    template <typename Derived> friend class internal::Comm;

    protected:
        static int getRemoteSize(const MPI_Comm& comm)
        {
            MPI_Int size;
            MPIWRAP_CALL(MPI_Comm_remote_size(comm, &size));
            return size;
        }

        explicit Intercomm(const MPI_Comm& comm)
        : internal::Comm<Intercomm>(comm, getRemoteSize(comm)), remoteSize(getRemoteSize(comm))
        {
            MPIWRAP_ASSERT(isIntercommunicator(), "comm must be an inter-communicator.");
        }

    public:
        const MPI_Int remoteSize;

#if MPIWRAP_CXX11

        Intercomm(Intercomm&& other)
        : internal::Comm<Intercomm>(std::move(other)), remoteSize(other.remoteSize) {}

#endif

        using internal::Comm<Intercomm>::rank;
        using internal::Comm<Intercomm>::size;
        using internal::Comm<Intercomm>::group;
        using internal::Comm<Intercomm>::compare;
        using internal::Comm<Intercomm>::duplicate;
        using internal::Comm<Intercomm>::subset;
        using internal::Comm<Intercomm>::split;
        using internal::Comm<Intercomm>::isIntercommunicator;
        using internal::Comm<Intercomm>::isIntracommunicator;
        using internal::Comm<Intercomm>::Recv_init;
        using internal::Comm<Intercomm>::Bsend_init;
        using internal::Comm<Intercomm>::Rsend_init;
        using internal::Comm<Intercomm>::Send_init;
        using internal::Comm<Intercomm>::Sendrecv;
        using internal::Comm<Intercomm>::Ssend_init;
        using internal::Comm<Intercomm>::Probe;
        using internal::Comm<Intercomm>::Recv;
        using internal::Comm<Intercomm>::Bsend;
        using internal::Comm<Intercomm>::Rsend;
        using internal::Comm<Intercomm>::Send;
        using internal::Comm<Intercomm>::Ssend;
        using internal::Comm<Intercomm>::Iprobe;
        using internal::Comm<Intercomm>::Irecv;
        using internal::Comm<Intercomm>::Ibsend;
        using internal::Comm<Intercomm>::Irsend;
        using internal::Comm<Intercomm>::Isend;
        using internal::Comm<Intercomm>::Issend;
        using internal::Comm<Intercomm>::Allgather;
        using internal::Comm<Intercomm>::Allreduce;
        using internal::Comm<Intercomm>::Alltoall;
        using internal::Comm<Intercomm>::Barrier;
        using internal::Comm<Intercomm>::Reduce_scatter;
        using internal::Comm<Intercomm>::operator MPI_Comm&;
        using internal::Comm<Intercomm>::operator const MPI_Comm&;
        using internal::Comm<Intercomm>::operator MPI_Comm*;
        using internal::Comm<Intercomm>::operator const MPI_Comm*;

#if MPIWRAP_HAVE_MPI_MESSAGE

        using internal::Comm<Intercomm>::Mprobe;
        using internal::Comm<Intercomm>::Improbe;

#endif

#if MPIWRAP_HAVE_MPI_ICOLLECTIVES

        using internal::Comm<Intercomm>::Iallgather;
        using internal::Comm<Intercomm>::Iallreduce;
        using internal::Comm<Intercomm>::Ialltoall;
        using internal::Comm<Intercomm>::Ibarrier;
        using internal::Comm<Intercomm>::Ireduce_scatter;

#endif

        Group remoteGroup() const
        {
            MPI_Group g;
            MPIWRAP_CALL(MPI_Comm_remote_group(comm, &g));
            return Group(g);
        }

        Intracomm merge(bool high = false) const
        {
            MPI_Comm c;
            MPIWRAP_CALL(MPI_Intercomm_merge(comm, high, &c));
            return Intracomm(c);
        }

        Intracomm localComm() const
        {
            return merge().subset(group());
        }

        /*
         * MPI_Bcast group A non-root
         */

        void Bcast() const
        {
            MPIWRAP_CALL(MPI_Bcast(NULL, 0, MPI_BYTE, MPI_PROC_NULL, comm));
        }

        /*
         * MPI_Bcast group A root
         */

        template <typename T>
        void Bcast(const T* buffer, MPI_Int count) const
        {
            Bcast(buffer, count, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Bcast(const std::vector<T>& buffer) const
        {
            Bcast(buffer, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Bcast(const T* buffer, MPI_Int count, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Bcast(nc(buffer), count, type, MPI_ROOT, comm));
        }

        template <typename T>
        void Bcast(const std::vector<T>& buffer, const Datatype& type) const
        {
            Bcast(&buffer.front(), buffer.size(), type);
        }

        /*
         * MPI_Bcast group B
         */

        template <typename T>
        void Bcast(T* buffer, MPI_Int count, MPI_Int root) const
        {
            Bcast(buffer, count, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Bcast(std::vector<T>& buffer, MPI_Int root) const
        {
            Bcast(buffer, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Bcast(T* buffer, MPI_Int count, MPI_Int root, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Bcast(buffer, count, type, root, comm));
        }

        template <typename T>
        void Bcast(std::vector<T>& buffer, MPI_Int root, const Datatype& type) const
        {
            Bcast(&buffer.front(), buffer.size(), root, type);
        }

        /*
         * MPI_Gather non-root group A
         */

        void Gather() const
        {
            MPIWRAP_CALL(MPI_Gather(NULL, 0, MPI_BYTE, NULL, 0, MPI_BYTE, MPI_PROC_NULL, comm));
        }

        /*
         * MPI_Gather root group A
         */

        template <typename T>
        void Gather(T* recvbuf, MPI_Int count) const
        {
            Gather(recvbuf, count, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Gather(std::vector<T>& recvbuf) const
        {
            Gather(recvbuf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Gather(T* recvbuf, MPI_Int count, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Gather(NULL, 0, type, recvbuf, count, type, MPI_ROOT, comm));
        }

        template <typename T>
        void Gather(std::vector<T>& recvbuf, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvbuf.size()%npeers == 0,
                           "Receive buffer size must be a multiple of the remote communicator size.");
            Gather(&recvbuf.front(), recvbuf.size()/npeers, type);
        }

        /*
         * MPI_Gather group B
         */

        template <typename T>
        void Gather(const T* sendbuf, MPI_Int count, MPI_Int root) const
        {
            Gather(sendbuf, count, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Gather(const std::vector<T>& sendbuf, MPI_Int root) const
        {
            Gather(sendbuf, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Gather(const T* sendbuf, MPI_Int count, MPI_Int root, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Gather(sendbuf, count, type, NULL, 0, type, root, comm));
        }

        template <typename T>
        void Gather(const std::vector<T>& sendbuf, MPI_Int root, const Datatype& type) const
        {
            Gather(&sendbuf.front(), sendbuf.size(), root, type);
        }

        /*
         * MPI_Gatherv non-root group A
         */

        void Gatherv() const
        {
            MPIWRAP_CALL(MPI_Gatherv(NULL, 0, MPI_BYTE, NULL, NULL, NULL, MPI_BYTE, MPI_PROC_NULL, comm));
        }

        /*
         * MPI_Gather root group A
         */

        template <typename T>
        void Gather(T* recvbuf, const MPI_Int* recvcounts) const
        {
            Gather(recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Gather(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts) const
        {
            Gather(recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Gather(T* recvbuf, const MPI_Int* recvcounts, const Datatype& type) const
        {
            std::vector<MPI_Int> recvdispls = displacements(recvcounts);
            Gather(recvbuf, recvcounts, &recvdispls.front(), type);
        }

        template <typename T>
        void Gather(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvbuf.size() == sum(recvcounts),
                           "The receive buffer size must equal the sum of the receive counts.");
            Gather(recvbuf, recvcounts, displacements(recvcounts), type);
        }

        template <typename T>
        void Gather(T* recvbuf, const MPI_Int* recvcounts, const MPI_Int* recvdispls) const
        {
            Gather(recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Gather(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts, const std::vector<MPI_Int>& recvdispls) const
        {
            Gather(recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Gather(T* recvbuf, const MPI_Int* recvcounts, const MPI_Int* recvdispls, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Gatherv(NULL, 0, type, recvbuf, recvcounts, recvdispls, type, MPI_ROOT, comm));
        }

        template <typename T>
        void Gather(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts,
                    const std::vector<MPI_Int>& recvdispls, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvcounts.size() == npeers,
                           "There must be exactly one receive count for each remote process.");
            MPIWRAP_ASSERT(recvdispls.size() == npeers,
                           "There must be exactly one receive displacement for each remote process.");
            MPIWRAP_ASSERT(recvbuf.size() >= recvdispls[npeers-1]+recvcounts[npeers-1],
                           "The receive buffer size must be at least as large the sum of last receive count and the last receive displacement.");
            Gather(&recvbuf.front(), &recvcounts.front(), &recvdispls.front(), type);
        }

        /*
         * MPI_Gather group B
         */

        template <typename T>
        void Gatherv(const T* sendbuf, MPI_Int count, MPI_Int root) const
        {
            Gather(sendbuf, count, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Gatherv(const std::vector<T>& sendbuf, MPI_Int root) const
        {
            Gather(sendbuf, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Gatherv(const T* sendbuf, MPI_Int count, MPI_Int root, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Gatherv(sendbuf, count, type, NULL, NULL, NULL, type, root, comm));
        }

        template <typename T>
        void Gatherv(const std::vector<T>& sendbuf, MPI_Int root, const Datatype& type) const
        {
            Gather(&sendbuf.front(), sendbuf.size(), root, type);
        }

        /*
         * MPI_Reduce group A non-root
         */

        void Reduce() const
        {
            MPIWRAP_CALL(MPI_Reduce(NULL, NULL, 0, MPI_BYTE, MPI_OP_NULL, MPI_PROC_NULL, comm));
        }

        /*
         * MPI_Reduce group A root
         */

        template <typename T>
        void Reduce(T* recvbuf, MPI_Int count, const MPI_Op& op) const
        {
            Reduce(recvbuf, count, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Reduce(std::vector<T>& recvbuf, const MPI_Op& op) const
        {
            Reduce(recvbuf, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Reduce(T* recvbuf, MPI_Int count, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Reduce(NULL, recvbuf, count, type, op, MPI_ROOT, comm));
        }

        template <typename T>
        void Reduce(std::vector<T>& recvbuf, const MPI_Op& op, const Datatype& type) const
        {
            Reduce(&recvbuf.front(), recvbuf.size(), op, type);
        }

        /*
         * MPI_Reduce group B
         */

        template <typename T>
        void Reduce(const T* sendbuf, MPI_Int count, const MPI_Op& op, MPI_Int root) const
        {
            Reduce(sendbuf, count, op, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Reduce(const std::vector<T>& sendbuf, const MPI_Op& op, MPI_Int root) const
        {
            Reduce(sendbuf, op, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Reduce(const T* sendbuf, MPI_Int count, const MPI_Op& op, MPI_Int root, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Reduce(sendbuf, NULL, count, type, op, root, comm));
        }

        template <typename T>
        void Reduce(const std::vector<T>& sendbuf, const MPI_Op& op, MPI_Int root, const Datatype& type) const
        {
            Reduce(&sendbuf.front(), sendbuf.size(), op, root, type);
        }

        /*
         * MPI_Scatter group A non-root
         */

        void Scatter() const
        {
            MPIWRAP_CALL(MPI_Scatter(NULL, 0, MPI_BYTE, NULL, 0, MPI_BYTE, MPI_PROC_NULL, comm));
        }

        /*
         * MPI_Scatter group A root
         */

        template <typename T>
        void Scatter(const T* sendbuf, MPI_Int recvcount) const
        {
            Scatter(sendbuf, recvcount, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Scatter(const std::vector<T>& sendbuf) const
        {
            Scatter(sendbuf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Scatter(const T* sendbuf, MPI_Int recvcount, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Scatter(sendbuf, recvcount, type, NULL, 0, type, MPI_ROOT, comm));
        }

        template <typename T>
        void Scatter(const std::vector<T>& sendbuf, const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendbuf.size()%npeers != 0,
                           "Send buffer size must a multiple of the remote group size.");
            Scatter(&sendbuf.front(), sendbuf.size()/npeers, type);
        }

        /*
         * MPI_Scatter group B
         */

        template <typename T>
        void Scatter(T* recvbuf, MPI_Int recvcount, MPI_Int root) const
        {
            Scatter(recvbuf, recvcount, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Scatter(std::vector<T>& recvbuf, MPI_Int root) const
        {
            Scatter(recvbuf, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Scatter(T* recvbuf, MPI_Int recvcount, MPI_Int root, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Scatter(NULL, 0, type, recvbuf, recvcount, type, root, comm));
        }

        template <typename T>
        void Scatter(std::vector<T>& recvbuf, MPI_Int root, const Datatype& type) const
        {
            Scatter(&recvbuf.front(), recvbuf.size(), root, type);
        }

        /*
         * MPI_Scatterv group A non-root
         */

        void Scatterv() const
        {
            MPIWRAP_CALL(MPI_Scatterv(NULL, NULL, NULL, MPI_BYTE, NULL, 0, MPI_BYTE, MPI_PROC_NULL, comm));
        }

        /*
         * MPI_Scatterv group A root
         */

        template <typename T>
        void Scatter(const T* sendbuf, const MPI_Int* sendcounts) const
        {
            Scatter(sendbuf, sendcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Scatter(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts) const
        {
            Scatter(sendbuf, sendcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Scatter(const T* sendbuf, const MPI_Int* sendcounts,
                      const Datatype& type) const
        {
            std::vector<MPI_Int> senddispls = displacements(sendcounts);
            Scatter(sendbuf, sendcounts, &senddispls.front(), type);
        }

        template <typename T>
        void Scatter(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts,
                     const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendbuf.size() == sum(sendcounts),
                           "Send buffer size must equal the sum of the send counts.");
            Scatter(sendbuf, sendcounts, displacements(sendcounts), type);
        }

        template <typename T>
        void Scatter(const T* sendbuf, const MPI_Int* sendcounts, const MPI_Int* senddispls) const
        {
            Scatter(sendbuf, sendcounts, senddispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Scatter(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts,
                     const std::vector<MPI_Int>& senddispls) const
        {
            Scatter(sendbuf, sendcounts, senddispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Scatter(const T* sendbuf, const MPI_Int* sendcounts, const MPI_Int* senddispls,
                      const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Scatterv(sendbuf, sendcounts, senddispls, type, NULL, 0, type, MPI_ROOT, comm));
        }

        template <typename T>
        void Scatter(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts,
                     const std::vector<MPI_Int>& senddispls, const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendcounts.size() == npeers,
                           "There must be exactly one send count for each remote process.");
            MPIWRAP_ASSERT(senddispls.size() == npeers,
                           "There must be exactly one send displacement for each remote process.");
            MPIWRAP_ASSERT(sendbuf.size() >= senddispls[npeers-1]+sendcounts[npeers],
                           "Send buffer size must be at least the sum of the last send count and the last send displacement.");
            Scatter(&sendbuf.front(), &sendcounts.front(), &senddispls.front(), type);
        }

        /*
         * MPI_Scatterv group B
         */

        template <typename T>
        void Scatterv(T* recvbuf, MPI_Int recvcount, MPI_Int root) const
        {
            Scatterv(recvbuf, recvcount, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Scatterv(std::vector<T>& recvbuf, MPI_Int root) const
        {
            Scatterv(recvbuf, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Scatterv(T* recvbuf, MPI_Int recvcount, MPI_Int root, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Scatterv(NULL, NULL, NULL, type, recvbuf, recvcount, type, root, comm));
        }

        template <typename T>
        void Scatterv(std::vector<T>& recvbuf, MPI_Int root, const Datatype& type) const
        {
            Scatterv(&recvbuf.front(), recvbuf.size(), root, type);
        }

#if MPIWRAP_HAVE_MPI_ICOLLECTIVES

        /*
         * MPI_Ibcast group A non-root
         */

        Request Ibcast() const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Ibcast(NULL, 0, MPI_BYTE, MPI_PROC_NULL, comm, &req));
            return Request(req);
        }

        /*
         * MPI_Ibcast group A root
         */

        template <typename T>
        Request Ibcast(const T* buffer, MPI_Int count) const
        {
            return Ibcast(buffer, count, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ibcast(const std::vector<T>& buffer) const
        {
            return Ibcast(buffer, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ibcast(const T* buffer, MPI_Int count, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Ibcast(nc(buffer), count, type, MPI_ROOT, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Ibcast(const std::vector<T>& buffer, const Datatype& type) const
        {
            return Ibcast(&buffer.front(), buffer.size(), type);
        }

        /*
         * MPI_Ibcast group B
         */

        template <typename T>
        Request Ibcast(T* buffer, MPI_Int count, MPI_Int root) const
        {
            return Ibcast(buffer, count, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ibcast(std::vector<T>& buffer, MPI_Int root) const
        {
            return Ibcast(buffer, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ibcast(T* buffer, MPI_Int count, MPI_Int root, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Ibcast(buffer, count, type, root, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Ibcast(std::vector<T>& buffer, MPI_Int root, const Datatype& type) const
        {
            return Ibcast(&buffer.front(), buffer.size(), root, type);
        }

        /*
         * MPI_Igather non-root group A
         */

        Request Igather() const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Igather(NULL, 0, MPI_BYTE, NULL, 0, MPI_BYTE, MPI_PROC_NULL, comm, &req));
            return Request(req);
        }

        /*
         * MPI_Igather root group A
         */

        template <typename T>
        Request Igather(T* recvbuf, MPI_Int count) const
        {
            return Igather(recvbuf, count, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Igather(std::vector<T>& recvbuf) const
        {
            return Igather(recvbuf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Igather(T* recvbuf, MPI_Int count, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Igather(NULL, 0, type, recvbuf, count, type, MPI_ROOT, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Igather(std::vector<T>& recvbuf, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_ASSERT(recvbuf.size()%npeers == 0,
                           "Receive buffer size must be a multiple of the remote communicator size.");
            return Igather(&recvbuf.front(), recvbuf.size()/npeers, type);
        }

        /*
         * MPI_Igather group B
         */

        template <typename T>
        Request Igather(const T* sendbuf, MPI_Int count, MPI_Int root) const
        {
            return Igather(sendbuf, count, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Igather(const std::vector<T>& sendbuf, MPI_Int root) const
        {
            return Igather(sendbuf, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Igather(const T* sendbuf, MPI_Int count, MPI_Int root, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Igather(sendbuf, count, type, NULL, 0, type, root, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Igather(const std::vector<T>& sendbuf, MPI_Int root, const Datatype& type) const
        {
            return Igather(&sendbuf.front(), sendbuf.size(), root, type);
        }

        /*
         * MPI_Igatherv non-root group A
         */

        Request Igatherv() const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Igatherv(NULL, 0, MPI_BYTE, NULL, NULL, NULL, MPI_BYTE, MPI_PROC_NULL, comm, &req));
            return Request(req);
        }

        /*
         * MPI_Igatherv root group A
         */

        template <typename T>
        Request Igather(T* recvbuf, const MPI_Int* recvcounts) const
        {
            return Igather(recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Igather(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts) const
        {
            return Igather(recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Igather(T* recvbuf, const MPI_Int* recvcounts, const Datatype& type) const
        {
            std::vector<MPI_Int> recvdispls = displacements(recvcounts);
            return Igather(recvbuf, recvcounts, &recvdispls.front(), type);
        }

        template <typename T>
        Request Igather(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_ASSERT(recvbuf.size() == sum(recvcounts),
                           "The receive buffer size must equal the sum of the receive counts.");
            return Igather(recvbuf, recvcounts, displacements(recvcounts), type);
        }

        template <typename T>
        Request Igather(T* recvbuf, const MPI_Int* recvcounts, const MPI_Int* recvdispls) const
        {
            return Igather(recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Igather(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts, const std::vector<MPI_Int>& recvdispls) const
        {
            return Igather(recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Igather(T* recvbuf, const MPI_Int* recvcounts, const MPI_Int* recvdispls, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Igatherv(NULL, 0, type, recvbuf, recvcounts, recvdispls, type, MPI_ROOT, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Igather(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts,
                    const std::vector<MPI_Int>& recvdispls, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvcounts.size() == npeers,
                           "There must be exactly one receive count for each remote process.");
            MPIWRAP_ASSERT(recvdispls.size() == npeers,
                           "There must be exactly one receive displacement for each remote process.");
            MPIWRAP_ASSERT(recvbuf.size() >= recvdispls[npeers-1]+recvcounts[npeers-1],
                           "The receive buffer size must be at least as large the sum of last receive count and the last receive displacement.");
            return Igather(&recvbuf.front(), &recvcounts.front(), &recvdispls.front(), type);
        }

        /*
         * MPI_Igatherv group B
         */

        template <typename T>
        Request Igatherv(const T* sendbuf, MPI_Int count, MPI_Int root) const
        {
            return Igather(sendbuf, count, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Igatherv(const std::vector<T>& sendbuf, MPI_Int root) const
        {
            return Igather(sendbuf, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Igatherv(const T* sendbuf, MPI_Int count, MPI_Int root, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Igatherv(sendbuf, count, type, NULL, NULL, NULL, type, root, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Igatherv(const std::vector<T>& sendbuf, MPI_Int root, const Datatype& type) const
        {
            return Igather(&sendbuf.front(), sendbuf.size(), root, type);
        }

        /*
         * MPI_Ireduce group A non-root
         */

        Request Ireduce() const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Ireduce(NULL, NULL, 0, MPI_BYTE, MPI_OP_NULL, MPI_PROC_NULL, comm, &req));
            return Request(req);
        }

        /*
         * MPI_Ireduce group A root
         */

        template <typename T>
        Request Ireduce(T* recvbuf, MPI_Int count, const MPI_Op& op) const
        {
            return Ireduce(recvbuf, count, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ireduce(std::vector<T>& recvbuf, const MPI_Op& op) const
        {
            return Ireduce(recvbuf, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ireduce(T* recvbuf, MPI_Int count, const MPI_Op& op, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Ireduce(NULL, recvbuf, count, type, op, MPI_ROOT, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Ireduce(std::vector<T>& recvbuf, const MPI_Op& op, const Datatype& type) const
        {
            return Ireduce(&recvbuf.front(), recvbuf.size(), op, type);
        }

        /*
         * MPI_Ireduce group B
         */

        template <typename T>
        Request Ireduce(const T* sendbuf, MPI_Int count, const MPI_Op& op, MPI_Int root) const
        {
            return Ireduce(sendbuf, count, op, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ireduce(const std::vector<T>& sendbuf, const MPI_Op& op, MPI_Int root) const
        {
            return Ireduce(sendbuf, op, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ireduce(const T* sendbuf, MPI_Int count, const MPI_Op& op, MPI_Int root, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Ireduce(sendbuf, NULL, count, type, op, root, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Ireduce(const std::vector<T>& sendbuf, const MPI_Op& op, MPI_Int root, const Datatype& type) const
        {
            return Ireduce(&sendbuf.front(), sendbuf.size(), op, root, type);
        }

        /*
         * MPI_Iscatter group A non-root
         */

        Request Iscatter() const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Iscatter(NULL, 0, MPI_BYTE, NULL, 0, MPI_BYTE, MPI_PROC_NULL, comm, &req));
            return Request(req);
        }

        /*
         * MPI_Iscatter group A root
         */

        template <typename T>
        Request Iscatter(const T* sendbuf, MPI_Int recvcount) const
        {
            return Iscatter(sendbuf, recvcount, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iscatter(const std::vector<T>& sendbuf) const
        {
            return Iscatter(sendbuf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iscatter(const T* sendbuf, MPI_Int recvcount, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Iscatter(sendbuf, recvcount, type, NULL, 0, type, MPI_ROOT, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Iscatter(const std::vector<T>& sendbuf, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_ASSERT(sendbuf.size()%npeers != 0,
                           "Send buffer size must a multiple of the remote group size.");
            return Iscatter(&sendbuf.front(), sendbuf.size()/npeers, type);
        }

        /*
         * MPI_Iscatter group B
         */

        template <typename T>
        Request Iscatter(T* recvbuf, MPI_Int recvcount, MPI_Int root) const
        {
            return Iscatter(recvbuf, recvcount, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iscatter(std::vector<T>& recvbuf, MPI_Int root) const
        {
            return Iscatter(recvbuf, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iscatter(T* recvbuf, MPI_Int recvcount, MPI_Int root, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Iscatter(NULL, 0, type, recvbuf, recvcount, type, root, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Iscatter(std::vector<T>& recvbuf, MPI_Int root, const Datatype& type) const
        {
            return Iscatter(&recvbuf.front(), recvbuf.size(), root, type);
        }

        /*
         * MPI_Iscatterv group A non-root
         */

        Request Iscatterv() const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Iscatterv(NULL, NULL, NULL, MPI_BYTE, NULL, 0, MPI_BYTE, MPI_PROC_NULL, comm, &req));
            return Request(req);
        }

        /*
         * MPI_Scatterv group A root
         */

        template <typename T>
        Request Iscatter(const T* sendbuf, const MPI_Int* sendcounts) const
        {
            return Iscatter(sendbuf, sendcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iscatter(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts) const
        {
            return Iscatter(sendbuf, sendcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iscatter(const T* sendbuf, const MPI_Int* sendcounts,
                      const Datatype& type) const
        {
            std::vector<MPI_Int> senddispls = displacements(sendcounts);
            return Iscatter(sendbuf, sendcounts, &senddispls.front(), type);
        }

        template <typename T>
        Request Iscatter(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts,
                     const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_ASSERT(sendbuf.size() == sum(sendcounts),
                           "Send buffer size must equal the sum of the send counts.");
            return Iscatter(sendbuf, sendcounts, displacements(sendcounts), type);
        }

        template <typename T>
        Request Iscatter(const T* sendbuf, const MPI_Int* sendcounts, const MPI_Int* senddispls) const
        {
            return Iscatter(sendbuf, sendcounts, senddispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iscatter(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts,
                     const std::vector<MPI_Int>& senddispls) const
        {
            return Iscatter(sendbuf, sendcounts, senddispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iscatter(const T* sendbuf, const MPI_Int* sendcounts, const MPI_Int* senddispls,
                      const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Iscatterv(sendbuf, sendcounts, senddispls, type, NULL, 0, type, MPI_ROOT, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Iscatter(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts,
                     const std::vector<MPI_Int>& senddispls, const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendcounts.size() == npeers,
                           "There must be exactly one send count for each remote process.");
            MPIWRAP_ASSERT(senddispls.size() == npeers,
                           "There must be exactly one send displacement for each remote process.");
            MPIWRAP_ASSERT(sendbuf.size() >= senddispls[npeers-1]+sendcounts[npeers],
                           "Send buffer size must be at least the sum of the last send count and the last send displacement.");
            return Iscatter(&sendbuf.front(), &sendcounts.front(), &senddispls.front(), type);
        }

        /*
         * MPI_Iscatterv group B
         */

        template <typename T>
        Request Iscatterv(T* recvbuf, MPI_Int recvcount, MPI_Int root) const
        {
            return Iscatterv(recvbuf, recvcount, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iscatterv(std::vector<T>& recvbuf, MPI_Int root) const
        {
            return Iscatterv(recvbuf, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iscatterv(T* recvbuf, MPI_Int recvcount, MPI_Int root, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Iscatterv(NULL, NULL, NULL, type, recvbuf, recvcount, type, root, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Iscatterv(std::vector<T>& recvbuf, MPI_Int root, const Datatype& type) const
        {
            return Iscatterv(&recvbuf.front(), recvbuf.size(), root, type);
        }

#endif
};

inline Intercomm Intracomm::intercomm(MPI_Int leader) const
{
    MPI_Comm c;
    MPIWRAP_CALL(MPI_Intercomm_create(comm, leader, MPI_COMM_NULL, MPI_PROC_NULL, MPI_ANY_TAG, &c));
    return Intercomm(c);
}

inline Intercomm Intracomm::intercomm(MPI_Int leader, const Intracomm& peerComm, MPI_Int remoteLeader, MPI_Int tag) const
{
    MPI_Comm c;
    MPIWRAP_CALL(MPI_Intercomm_create(comm, leader, peerComm, remoteLeader, tag, &c));
    return Intercomm(c);
}

}

#endif

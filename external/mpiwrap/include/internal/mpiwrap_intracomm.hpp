#ifndef _MPIWRAP_INTRACOMM_HPP_
#define _MPIWRAP_INTRACOMM_HPP_

#include "mpiwrap_common.hpp"
#include "mpiwrap_datatype.hpp"
#include "mpiwrap_status.hpp"
#include "mpiwrap_request.hpp"
#include "mpiwrap_comm.hpp"

namespace MPIWrap
{

class Intracomm : protected internal::Comm<Intracomm>
{
    friend class Intercomm;
    template <typename Derived> friend class internal::Comm;

    protected:
        explicit Intracomm(const MPI_Comm& comm)
        : internal::Comm<Intracomm>(comm, getSize(comm))
        {
            MPIWRAP_ASSERT(!isIntercommunicator(), "comm must be an intra-communicator.");
        }

    public:
#if MPIWRAP_CXX11

        Intracomm(Intracomm&& other) : internal::Comm<Intracomm>(std::move(other)) {}

#endif

        using internal::Comm<Intracomm>::rank;
        using internal::Comm<Intracomm>::size;
        using internal::Comm<Intracomm>::group;
        using internal::Comm<Intracomm>::compare;
        using internal::Comm<Intracomm>::duplicate;
        using internal::Comm<Intracomm>::subset;
        using internal::Comm<Intracomm>::split;
        using internal::Comm<Intracomm>::isIntercommunicator;
        using internal::Comm<Intracomm>::isIntracommunicator;
        using internal::Comm<Intracomm>::Recv_init;
        using internal::Comm<Intracomm>::Bsend_init;
        using internal::Comm<Intracomm>::Rsend_init;
        using internal::Comm<Intracomm>::Send_init;
        using internal::Comm<Intracomm>::Sendrecv;
        using internal::Comm<Intracomm>::Ssend_init;
        using internal::Comm<Intracomm>::Probe;
        using internal::Comm<Intracomm>::Recv;
        using internal::Comm<Intracomm>::Bsend;
        using internal::Comm<Intracomm>::Rsend;
        using internal::Comm<Intracomm>::Send;
        using internal::Comm<Intracomm>::Ssend;
        using internal::Comm<Intracomm>::Iprobe;
        using internal::Comm<Intracomm>::Irecv;
        using internal::Comm<Intracomm>::Ibsend;
        using internal::Comm<Intracomm>::Irsend;
        using internal::Comm<Intracomm>::Isend;
        using internal::Comm<Intracomm>::Issend;
        using internal::Comm<Intracomm>::Allgather;
        using internal::Comm<Intracomm>::Allreduce;
        using internal::Comm<Intracomm>::Alltoall;
        using internal::Comm<Intracomm>::Barrier;
        using internal::Comm<Intracomm>::Reduce_scatter;
        using internal::Comm<Intracomm>::operator MPI_Comm&;
        using internal::Comm<Intracomm>::operator const MPI_Comm&;
        using internal::Comm<Intracomm>::operator MPI_Comm*;
        using internal::Comm<Intracomm>::operator const MPI_Comm*;

#if MPIWRAP_HAVE_MPI_MESSAGE

        using internal::Comm<Intracomm>::Mprobe;
        using internal::Comm<Intracomm>::Improbe;

#endif

#if MPIWRAP_HAVE_MPI_ICOLLECTIVES

        using internal::Comm<Intracomm>::Iallgather;
        using internal::Comm<Intracomm>::Iallreduce;
        using internal::Comm<Intracomm>::Ialltoall;
        using internal::Comm<Intracomm>::Ibarrier;
        using internal::Comm<Intracomm>::Ireduce_scatter;

#endif

        /*
         * OpenMPI 1.7.2 doesn't seem to support this.
         *
        Intracomm subset(const Group& group, MPI_Int tag) const
        {
            MPI_Comm c;
            MPIWRAP_CALL(MPI_Comm_create_group(comm, group, tag, &c));
            return Intracomm(c);
        }
        */

        Intercomm intercomm(MPI_Int leader) const;

        Intercomm intercomm(MPI_Int leader, const Intracomm& peerComm, MPI_Int remoteLeader, MPI_Int tag) const;

        static Intracomm world()
        {
            return Intracomm(MPI_COMM_WORLD);
        }

        static Intracomm self()
        {
            return Intracomm(MPI_COMM_SELF);
        }

        /*
         * MPI_Allgather in-place
         */

        template <typename T>
        void Allgather(T* recvbuf, MPI_Int count) const
        {
            Allgather(recvbuf, count, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Allgather(std::vector<T>& recvbuf) const
        {
            Allgather(recvbuf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Allgather(T* recvbuf, MPI_Int count, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Allgather(MPI_IN_PLACE, 0, type, recvbuf, count, type, comm));
        }

        template <typename T>
        void Allgather(std::vector<T>& recvbuf, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvbuf.size()%size == 0,
                           "Receive buffer size must be a multiple of the communicator size.");
            Allgather(&recvbuf.front(), recvbuf.size()/size, type);
        }

        /*
         * MPI_Allgatherv in-place
         */

        template <typename T>
        void Allgather(T* recvbuf, const MPI_Int* recvcounts) const
        {
            Allgather(recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Allgather(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts) const
        {
            Allgather(recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Allgather(T* recvbuf, const MPI_Int* recvcounts, const Datatype& type) const
        {
            std::vector<MPI_Int> displs = displacements(recvcounts);
            Allgather(recvbuf, recvcounts, displs, type);
        }

        template <typename T>
        void Allgather(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvbuf.size() == sum(recvcounts),
                           "The receive buffer size must equal the sum of the receive counts.");
            Allgather(recvbuf, recvcounts, displacements(recvcounts), type);
        }

        template <typename T>
        void Allgather(T* recvbuf, const MPI_Int* recvcounts, const MPI_Int* recvdispls) const
        {
            Allgather(recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Allgather(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts, const std::vector<MPI_Int>& recvdispls) const
        {
            Allgather(recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Allgather(T* recvbuf, const MPI_Int* recvcounts, const MPI_Int* recvdispls, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Allgatherv(MPI_IN_PLACE, 0, type, recvbuf, recvcounts, recvdispls, type, comm));
        }

        template <typename T>
        void Allgather(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts, const std::vector<MPI_Int>& recvdispls, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvcounts.size() == size,
                           "There must be exactly one receive count for each process.");
            MPIWRAP_ASSERT(recvdispls.size() == size,
                           "There must be exactly one receive displacement for each process.");
            MPIWRAP_ASSERT(recvbuf.size() >= recvdispls[size-1]+recvcounts[size-1],
                           "The receive buffer size must be at least as large the sum of last receive count and the last receive displacement.");
            Allgather(&recvbuf.front(), &recvcounts.front(), &recvdispls.front(), type);
        }

        /*
         * MPI_Allreduce in-place
         */

        template <typename T>
        void Allreduce(T* buf, MPI_Int count, const MPI_Op& op) const
        {
            Allreduce(buf, count, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Allreduce(std::vector<T>& buf, const MPI_Op& op) const
        {
            Allreduce(buf, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Allreduce(T* buf, MPI_Int count, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Allreduce(MPI_IN_PLACE, buf, count, type, op, comm));
        }

        template <typename T>
        void Allreduce(std::vector<T>& buf, const MPI_Op& op, const Datatype& type) const
        {
            Allreduce(&buf.front(), buf.size(), op, type);
        }

        /*
         * MPI_Alltoall in-place
         */

        template <typename T>
        void Alltoall(T* recvbuf, MPI_Int count) const
        {
            Alltoall(recvbuf, count, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Alltoall(std::vector<T>& recvbuf) const
        {
            Alltoall(recvbuf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Alltoall(T* recvbuf, MPI_Int count, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Alltoall(MPI_IN_PLACE, 0, type, recvbuf, count, type, comm));
        }

        template <typename T>
        void Alltoall(std::vector<T>& recvbuf, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvbuf.size()%size == 0,
                           "Receive buffer size must be a multiple of the communicator size.");
            Alltoall(&recvbuf.front(), recvbuf.size()/size, type);
        }

        /*
         * MPI_Alltoallv in-place
         */

        template <typename T>
        void Alltoall(T* recvbuf, const MPI_Int* recvcounts) const
        {
            Alltoall(recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Alltoall(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts) const
        {
            Alltoall(recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Alltoall(T* recvbuf, const MPI_Int* recvcounts, const Datatype& type) const
        {
            std::vector<MPI_Int> recvdispls = displacements(recvcounts);
            Alltoall(recvbuf, recvcounts, recvdispls, type);
        }

        template <typename T>
        void Alltoall(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts,
                      const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvbuf.size() == sum(recvcounts),
                           "The receive buffer size must equal the sum of the receive counts.");
            Alltoall(recvbuf, recvcounts, displacements(recvcounts), type);
        }

        template <typename T>
        void Alltoall(T* recvbuf, const MPI_Int* recvcounts, const MPI_Int* recvdispls) const
        {
            Alltoall(recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Alltoall(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts, const std::vector<MPI_Int>& recvdispls) const
        {
            Alltoall(recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Alltoall(T* recvbuf, const MPI_Int* recvcounts, const MPI_Int* recvdispls, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Alltoallv(MPI_IN_PLACE, NULL, NULL, type,
                                       recvbuf, recvcounts, recvdispls, type, comm));
        }

        template <typename T>
        void Alltoall(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts, const std::vector<MPI_Int>& recvdispls,
                      const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvcounts.size() == npeers,
                           "There must be exactly one receive count for each process.");
            MPIWRAP_ASSERT(recvdispls.size() == npeers,
                           "There must be exactly one receive displacement for each process.");
            MPIWRAP_ASSERT(recvbuf.size() >= recvdispls[npeers-1]+recvcounts[npeers-1],
                           "The receive buffer size must be at least as large the sum of last receive count and the last receive displacement.");
            Alltoall(&recvbuf.front(), &recvcounts.front(), &recvdispls.front(), type);
        }

        /*
         * MPI_Alltoallw in-place
         */

        //TODO

        /*
         * MPI_Bcast non-root
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
         * MPI_Bcast root
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
            MPIWRAP_CALL(MPI_Bcast(const_cast<T*>(buffer), count, type, rank, comm));
        }

        template <typename T>
        void Bcast(const std::vector<T>& buffer, const Datatype& type) const
        {
            Bcast(&buffer.front(), buffer.size(), type);
        }

        /*
         * MPI_Gather non-root
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
         * MPI_Gather root
         */

        template <typename T>
        void Gather(const T* sendbuf, T* recvbuf, MPI_Int count) const
        {
            Gather(sendbuf, recvbuf, count, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Gather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf) const
        {
            Gather(sendbuf, recvbuf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Gather(const T* sendbuf, T* recvbuf, MPI_Int count, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Gather(sendbuf, count, type, recvbuf, count, type, rank, comm));
        }

        template <typename T>
        void Gather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const Datatype& type) const
        {
            MPIWRAP_ASSERT(size*sendbuf.size() == recvbuf.size(),
                           "Receive buffers size must be equal communicator size times send buffer size.");
            Gather(&sendbuf.front(), &recvbuf.front(), sendbuf.size(), type);
        }

        /*
         * MPI_Gather root in-place
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
            MPIWRAP_CALL(MPI_Gather(MPI_IN_PLACE, 0, type, recvbuf, count, type, rank, comm));
        }

        template <typename T>
        void Gather(std::vector<T>& recvbuf, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvbuf.size()%size == 0,
                           "Receive buffer size must be a multiple of the communicator size.");
            Gather(&recvbuf.front(), recvbuf.size()/size, type);
        }

        /*
         * MPI_Gatherv non-root
         */

        template <typename T>
        void Gatherv(const T* sendbuf, MPI_Int sendcount, MPI_Int root) const
        {
            Gatherv(sendbuf, sendcount, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Gatherv(const std::vector<T>& sendbuf, MPI_Int root) const
        {
            Gatherv(sendbuf, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Gatherv(const T* sendbuf, MPI_Int sendcount, MPI_Int root, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Gatherv(sendbuf, sendcount, type, NULL, NULL, NULL, type, root, comm));
        }

        template <typename T>
        void Gatherv(const std::vector<T>& sendbuf, MPI_Int root, const Datatype& type) const
        {
            Gatherv(&sendbuf.front(), sendbuf.size(), root, type);
        }

        /*
         * MPI_Gatherv root
         */

        template <typename T>
        void Gather(const T* sendbuf, MPI_Int sendcount, T* recvbuf, const MPI_Int* recvcounts) const
        {
            Gather(sendbuf, sendcount, recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Gather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf,
                    const std::vector<MPI_Int>& recvcounts) const
        {
            Gather(sendbuf, recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Gather(const T* sendbuf, MPI_Int sendcount, T* recvbuf, const MPI_Int* recvcounts,
                     const Datatype& type) const
        {
            std::vector<MPI_Int> recvdispls = displacements(recvcounts);
            Gather(sendbuf, recvbuf, recvcounts, &recvdispls.front(), type);
        }

        template <typename T>
        void Gather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts,
                    const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvbuf.size() == sum(recvcounts),
                           "The receive buffer size must equal the sum of the receive counts.");
            Gather(sendbuf, recvbuf, recvcounts, displacements(recvcounts), type);
        }

        template <typename T>
        void Gather(const T* sendbuf, MPI_Int sendcount, T* recvbuf, const MPI_Int* recvcounts, const MPI_Int* recvdispls) const
        {
            Gather(sendbuf, sendcount, recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Gather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf,
                    const std::vector<MPI_Int>& recvcounts, const std::vector<MPI_Int>& recvdispls) const
        {
            Gather(sendbuf, recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Gather(const T* sendbuf, MPI_Int sendcount, T* recvbuf, const MPI_Int* recvcounts, const MPI_Int* recvdispls,
                     const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Gatherv(sendbuf, sendcount, type, recvbuf, recvcounts, recvdispls, type, rank, comm));
        }

        template <typename T>
        void Gather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts,
                    const std::vector<MPI_Int>& recvdispls, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvcounts.size() == size,
                           "There must be exactly one receive count for each process.");
            MPIWRAP_ASSERT(recvdispls.size() == size,
                           "There must be exactly one receive displacement for each process.");
            MPIWRAP_ASSERT(recvbuf.size() >= recvdispls[size-1]+recvcounts[size-1],
                           "The receive buffer size must be at least as large the sum of last receive count and the last receive displacement.");
            MPIWRAP_ASSERT(sendbuf.size() == recvcounts[rank],
                           "The send buffer size must equal the receive count for this process.");
            Gather(&sendbuf.front(), sendbuf.size(), &recvbuf.front(), &recvcounts.front(), &recvdispls.front(), type);
        }

        /*
         * MPI_Gatherv root in-place
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
            Gather(recvbuf, recvcounts, recvdispls, type);
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
            MPIWRAP_CALL(MPI_Gatherv(MPI_IN_PLACE, 0, type, recvbuf, recvcounts, recvdispls, type, rank, comm));
        }

        template <typename T>
        void Gather(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts, const std::vector<MPI_Int>& recvdispls, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvcounts.size() == size,
                           "There must be exactly one receive count for each process.");
            MPIWRAP_ASSERT(recvdispls.size() == size,
                           "There must be exactly one receive displacement for each process.");
            MPIWRAP_ASSERT(recvbuf.size() >= recvdispls[size-1]+recvcounts[size-1],
                           "The receive buffer size must be at least as large the sum of last receive count and the last receive displacement.");
            Gather(&recvbuf.front(), &recvcounts.front(), &recvdispls.front(), type);
        }

        /*
         * MPI_Reduce non-root
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
            MPIWRAP_CALL(MPI_Reduce(nc(sendbuf), NULL, count, type, op, root, comm));
        }

        template <typename T>
        void Reduce(const std::vector<T>& sendbuf, const MPI_Op& op, MPI_Int root, const Datatype& type) const
        {
            Reduce(&sendbuf.front(), sendbuf.size(), op, root, type);
        }

        /*
         * MPI_Reduce root
         */

        template <typename T>
        void Reduce(const T* sendbuf, T* recvbuf, MPI_Int count, const MPI_Op& op) const
        {
            Reduce(sendbuf, recvbuf, count, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Reduce(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI_Op& op) const
        {
            Reduce(sendbuf, recvbuf, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Reduce(const T* sendbuf, T* recvbuf, MPI_Int count, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Reduce(nc(sendbuf), recvbuf, count, type, op, rank, comm));
        }

        template <typename T>
        void Reduce(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendbuf.size() == recvbuf.size(),
                           "Send and receive buffers must be the same size.");
            Reduce(&sendbuf.front(), &recvbuf.front(), sendbuf.size(), op, type);
        }

        /*
         * MPI_Reduce root in-place
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
            MPIWRAP_CALL(MPI_Reduce(MPI_IN_PLACE, recvbuf, count, type, op, rank, comm));
        }

        template <typename T>
        void Reduce(std::vector<T>& recvbuf, const MPI_Op& op, const Datatype& type) const
        {
            Reduce(&recvbuf.front(), recvbuf.size(), op, type);
        }

        /*
         * MPI_Reduce_scatter in-place
         */

        template <typename T>
        void Reduce_scatter(T* recvbuf, MPI_Int* recvcounts, const MPI_Op& op) const
        {
            Reduce_scatter(recvbuf, recvcounts, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Reduce_scatter(std::vector<T>& recvbuf, std::vector<MPI_Int>& recvcounts, const MPI_Op& op) const
        {
            Reduce_scatter(recvbuf, recvcounts, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Reduce_scatter(T* recvbuf, MPI_Int* recvcounts, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Reduce_scatter(MPI_IN_PLACE, recvbuf, recvcounts, type, op, comm));
        }

        template <typename T>
        void Reduce_scatter(std::vector<T>& recvbuf, std::vector<MPI_Int>& recvcounts, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvcounts.size() == size,
                           "There must be exactly one receive count for each process.");
            MPIWRAP_ASSERT(recvbuf.size() == sum(recvcounts),
                           "The send/receive buffer size must equal the sum of the receive counts.");
            Reduce_scatter(&recvbuf.front(), &recvcounts.front(), op, type);
        }

        /*
         * MPI_Reduce_scatter_block in-place
         */

        template <typename T>
        void Reduce_scatter(T* recvbuf, MPI_Int count, const MPI_Op& op) const
        {
            Reduce_scatter(recvbuf, count, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Reduce_scatter(std::vector<T>& recvbuf, const MPI_Op& op) const
        {
            Reduce_scatter(recvbuf, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Reduce_scatter(T* recvbuf, MPI_Int count, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Reduce_scatter(MPI_IN_PLACE, recvbuf, count, type, op, comm));
        }

        template <typename T>
        void Reduce_scatter(std::vector<T>& recvbuf, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvbuf.size()%npeers == 0, "Receive buffer size must be a multiple of the number of processes.");
            Reduce_scatter(&recvbuf.front(), recvbuf.size()/npeers, op, type);
        }

        /*
         * MPI_Scatter non-root
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
         * MPI_Scatter root
         */

        template <typename T>
        void Scatter(const T* sendbuf, T* recvbuf, MPI_Int recvcount) const
        {
            Scatter(sendbuf, recvbuf, recvcount, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Scatter(const std::vector<T>& sendbuf, std::vector<T>& recvbuf) const
        {
            Scatter(sendbuf, recvbuf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Scatter(const T* sendbuf, T* recvbuf, MPI_Int recvcount, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Scatter(sendbuf, recvcount, type, recvbuf, recvcount, type, rank, comm));
        }

        template <typename T>
        void Scatter(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const Datatype& type) const
        {
            MPIWRAP_ASSERT(size*recvbuf.size() == sendbuf.size(),
                           "Send buffer size must be equal communicator size times send receive size.");
            Scatter(&sendbuf.front(), &recvbuf.front(), recvbuf.size(), type);
        }

        /*
         * MPI_Scatter root in-place
         */

        template <typename T>
        void Scatter(const T* sendbuf, MPI_Int sendcount) const
        {
            Scatter(sendbuf, sendcount, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Scatter(const std::vector<T>& sendbuf) const
        {
            Scatter(sendbuf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Scatter(const T* sendbuf, MPI_Int count, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Scatter(sendbuf, count, type, MPI_IN_PLACE, 0, type, rank, comm));
        }

        template <typename T>
        void Scatter(const std::vector<T>& sendbuf, const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendbuf.size()%size == 0,
                           "Send buffer size must be a multiple of the communicator size.");
            Scatter(&sendbuf.front(), sendbuf.size()/size, type);
        }

        /*
         * MPI_Scatterv non-root
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

        /*
         * MPI_Scatterv root
         */

        template <typename T>
        void Scatter(const T* sendbuf, const MPI_Int* sendcounts, T* recvbuf, MPI_Int recvcount) const
        {
            Scatter(sendbuf, sendcounts, recvbuf, recvcount, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Scatter(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts,
                     std::vector<T>& recvbuf) const
        {
            Scatter(sendbuf, sendcounts, recvbuf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Scatter(const T* sendbuf, const MPI_Int* sendcounts, T* recvbuf, MPI_Int recvcount,
                      const Datatype& type) const
        {
            std::vector<MPI_Int> senddispls = displacements(sendcounts);
            Scatter(sendbuf, sendcounts, &senddispls.front(), recvbuf, recvcount, type);
        }

        template <typename T>
        void Scatter(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts,
                     std::vector<T>& recvbuf, const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendbuf.size() == sum(sendcounts),
                           "Send buffer size must equal the sum of the send counts.");
            Scatter(sendbuf, sendcounts, displacements(sendcounts), recvbuf, type);
        }

        template <typename T>
        void Scatter(const T* sendbuf, const MPI_Int* sendcounts, const MPI_Int* senddispls, T* recvbuf, MPI_Int recvcount) const
        {
            Scatter(sendbuf, sendcounts, senddispls, recvbuf, recvcount, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Scatter(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts,
                     const std::vector<MPI_Int>& senddispls, std::vector<T>& recvbuf) const
        {
            Scatter(sendbuf, sendcounts, senddispls, recvbuf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Scatter(const T* sendbuf, const MPI_Int* sendcounts, const MPI_Int* senddispls, T* recvbuf, MPI_Int recvcount,
                      const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Scatterv(sendbuf, sendcounts, senddispls, type, recvbuf, recvcount, type, rank, comm));
        }

        template <typename T>
        void Scatter(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts,
                     const std::vector<MPI_Int>& senddispls, std::vector<T>& recvbuf, const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendcounts.size() == size,
                           "There must be exactly one send count for each process.");
            MPIWRAP_ASSERT(senddispls.size() == size,
                           "There must be exactly one send displacement for each process.");
            MPIWRAP_ASSERT(recvbuf.size() == sendcounts[rank],
                           "Receive buffer size must equal send count for this process.");
            MPIWRAP_ASSERT(sendbuf.size() >= senddispls[size-1]+sendcounts[size-1],
                           "Send buffer size must be at least the sum of the last send count and the last send displacement.");
            Scatter(&sendbuf.front(), &sendcounts.front(), &senddispls.front(), &recvbuf.front(), recvbuf.size(), type);
        }

        /*
         * MPI_Scatterv root in-place
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
        void Scatter(const T* sendbuf, const MPI_Int* sendcounts, const Datatype& type) const
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
        void Scatter(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts, const std::vector<MPI_Int>& senddispls) const
        {
            Scatter(sendbuf, sendcounts, senddispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Scatter(const T* sendbuf, const MPI_Int* sendcounts, const MPI_Int* senddispls, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Scatterv(sendbuf, sendcounts, senddispls, type, MPI_IN_PLACE, 0, type, rank, comm));
        }

        template <typename T>
        void Scatter(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts, const std::vector<MPI_Int>& senddispls,
                     const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendcounts.size() == size,
                           "There must be exactly one send count for each process.");
            MPIWRAP_ASSERT(senddispls.size() == size,
                           "There must be exactly one send displacement for each process.");
            MPIWRAP_ASSERT(sendbuf.size() >= senddispls[size-1]+sendcounts[size-1],
                           "Send buffer size must be at least the sum of the last send count and the last send displacement.");
            Scatter(&sendbuf.front(), &sendcounts.front(), &senddispls.front(), type);
        }

        /*
         * MPI_Exscan
         */

        template <typename T>
        void Exscan(const T* sendbuf, T* recvbuf, MPI_Int count, const MPI_Op& op) const
        {
            Exscan(sendbuf, recvbuf, count, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Exscan(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI_Op& op) const
        {
            Exscan(sendbuf, recvbuf, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Exscan(const T* sendbuf, T* recvbuf, MPI_Int count, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Exscan(sendbuf, recvbuf, count, type, op, comm));
        }

        template <typename T>
        void Exscan(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendbuf.size() == recvbuf.size(),
                           "Send and receive buffers must be the same size.");
            Exscan(&sendbuf.front(), &recvbuf.begin(), sendbuf.size(), op, type);
        }

        /*
         * MPI_Exscan in-place
         */

        template <typename T>
        void Exscan(T* recvbuf, MPI_Int count, const MPI_Op& op) const
        {
            Exscan(recvbuf, count, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Exscan(std::vector<T>& recvbuf, const MPI_Op& op) const
        {
            Exscan(recvbuf, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Exscan(T* recvbuf, MPI_Int count, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Exscan(MPI_IN_PLACE, recvbuf, count, type, op, comm));
        }

        template <typename T>
        void Exscan(std::vector<T>& recvbuf, const MPI_Op& op, const Datatype& type) const
        {
            Exscan(&recvbuf.front(), recvbuf.size(), op, type);
        }

        /*
         * MPI_Scan
         */

        template <typename T>
        void Scan(const T* sendbuf, T* recvbuf, MPI_Int count, const MPI_Op& op) const
        {
            Scan(sendbuf, recvbuf, count, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Scan(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI_Op& op) const
        {
            Scan(sendbuf, recvbuf, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Scan(const T* sendbuf, T* recvbuf, MPI_Int count, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Scan(sendbuf, recvbuf, count, type, op, comm));
        }

        template <typename T>
        void Scan(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendbuf.size() == recvbuf.size(),
                           "Send and receive buffers must be the same size.");
            Scan(&sendbuf.front(), &recvbuf.front(), sendbuf.size(), op, type);
        }

        /*
         * MPI_Scan in-place
         */

        template <typename T>
        void Scan(T* recvbuf, MPI_Int count, const MPI_Op& op) const
        {
            Scan(recvbuf, count, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Scan(std::vector<T>& recvbuf, const MPI_Op& op) const
        {
            Scan(recvbuf, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        void Scan(T* recvbuf, MPI_Int count, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_CALL(MPI_Scan(MPI_IN_PLACE, recvbuf, count, type, op, comm));
        }

        template <typename T>
        void Scan(std::vector<T>& recvbuf, const MPI_Op& op, const Datatype& type) const
        {
            Scan(&recvbuf.front(), recvbuf.size(), op, type);
        }

#if MPIWRAP_HAVE_MPI_ICOLLECTIVES

        /*
         * MPI_Iallgather in-place
         */

        template <typename T>
        Request Iallgather(T* recvbuf, MPI_Int count) const
        {
            return Iallgather(recvbuf, count, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iallgather(std::vector<T>& recvbuf) const
        {
            return Iallgather(recvbuf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iallgather(T* recvbuf, MPI_Int count, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Iallgather(MPI_IN_PLACE, 0, type, recvbuf, count, type, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Iallgather(std::vector<T>& recvbuf, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvbuf.size()%size == 0,
                           "Receive buffer size must be a multiple of the communicator size.");
            return Iallgather(&recvbuf.front(), recvbuf.size()/size, type);
        }

        /*
         * MPI_Iallgatherv in-place
         */

        template <typename T>
        Request Iallgather(T* recvbuf, const MPI_Int* recvcounts) const
        {
            return Iallgather(recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iallgather(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts) const
        {
            return Iallgather(recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iallgather(T* recvbuf, const MPI_Int* recvcounts, const Datatype& type) const
        {
            std::vector<MPI_Int> displs = displacements(recvcounts);
            return Iallgather(recvbuf, recvcounts, displs, type);
        }

        template <typename T>
        Request Iallgather(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvbuf.size() == sum(recvcounts),
                           "The receive buffer size must equal the sum of the receive counts.");
            return Iallgather(recvbuf, recvcounts, displacements(recvcounts), type);
        }

        template <typename T>
        Request Iallgather(T* recvbuf, const MPI_Int* recvcounts, const MPI_Int* recvdispls) const
        {
            return Iallgather(recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iallgather(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts, const std::vector<MPI_Int>& recvdispls) const
        {
            return Iallgather(recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iallgather(T* recvbuf, const MPI_Int* recvcounts, const MPI_Int* recvdispls, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Iallgatherv(MPI_IN_PLACE, 0, type, recvbuf, recvcounts, recvdispls, type, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Iallgather(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts, const std::vector<MPI_Int>& recvdispls, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvcounts.size() == size,
                           "There must be exactly one receive count for each process.");
            MPIWRAP_ASSERT(recvdispls.size() == size,
                           "There must be exactly one receive displacement for each process.");
            MPIWRAP_ASSERT(recvbuf.size() >= recvdispls[size-1]+recvcounts[size-1],
                           "The receive buffer size must be at least as large the sum of last receive count and the last receive displacement.");
            return Iallgather(&recvbuf.front(), &recvcounts.front(), &recvdispls.front(), type);
        }

        /*
         * MPI_Iallreduce in-place
         */

        template <typename T>
        Request Iallreduce(T* buf, MPI_Int count, const MPI_Op& op) const
        {
            return Iallreduce(buf, count, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iallreduce(std::vector<T>& buf, const MPI_Op& op) const
        {
            return Iallreduce(buf, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iallreduce(T* buf, MPI_Int count, const MPI_Op& op, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Iallreduce(MPI_IN_PLACE, buf, count, type, op, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Iallreduce(std::vector<T>& buf, const MPI_Op& op, const Datatype& type) const
        {
            return Iallreduce(&buf.front(), buf.size(), op, type);
        }

        /*
         * MPI_Ialltoall in-place
         */

        template <typename T>
        Request Ialltoall(T* recvbuf, MPI_Int count) const
        {
            return Ialltoall(recvbuf, count, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ialltoall(std::vector<T>& recvbuf) const
        {
            return Ialltoall(recvbuf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ialltoall(T* recvbuf, MPI_Int count, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Ialltoall(MPI_IN_PLACE, 0, type, recvbuf, count, type, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Ialltoall(std::vector<T>& recvbuf, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvbuf.size()%size == 0,
                           "Receive buffer size must be a multiple of the communicator size.");
            return Ialltoall(&recvbuf.front(), recvbuf.size()/size, type);
        }

        /*
         * MPI_Ialltoallv in-place
         */

        template <typename T>
        Request Ialltoall(T* recvbuf, const MPI_Int* recvcounts) const
        {
            return Ialltoall(recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ialltoall(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts) const
        {
            return Ialltoall(recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ialltoall(T* recvbuf, const MPI_Int* recvcounts, const Datatype& type) const
        {
            std::vector<MPI_Int> recvdispls = displacements(recvcounts);
            return Ialltoall(recvbuf, recvcounts, recvdispls, type);
        }

        template <typename T>
        Request Ialltoall(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts,
                       const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvbuf.size() == sum(recvcounts),
                           "The receive buffer size must equal the sum of the receive counts.");
            return Ialltoall(recvbuf, recvcounts, displacements(recvcounts), type);
        }

        template <typename T>
        Request Ialltoall(T* recvbuf, const MPI_Int* recvcounts, const MPI_Int* recvdispls) const
        {
            return Ialltoall(recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ialltoall(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts, const std::vector<MPI_Int>& recvdispls) const
        {
            return Ialltoall(recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ialltoall(T* recvbuf, const MPI_Int* recvcounts, const MPI_Int* recvdispls, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Ialltoallv(MPI_IN_PLACE, NULL, NULL, type,
                                        recvbuf, recvcounts, recvdispls, type, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Ialltoall(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts, const std::vector<MPI_Int>& recvdispls,
                       const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvcounts.size() == npeers,
                           "There must be exactly one receive count for each process.");
            MPIWRAP_ASSERT(recvdispls.size() == npeers,
                           "There must be exactly one receive displacement for each process.");
            MPIWRAP_ASSERT(recvbuf.size() >= recvdispls[npeers-1]+recvcounts[npeers-1],
                           "The receive buffer size must be at least as large the sum of last receive count and the last receive displacement.");
            return Ialltoall(&recvbuf.front(), &recvcounts.front(), &recvdispls.front(), type);
        }

        /*
         * MPI_Ialltoallw in-place
         */

        //TODO

        /*
         * MPI_Bcast non-root
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
         * MPI_Ibcast root
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
            MPIWRAP_CALL(MPI_Ibcast(const_cast<T*>(buffer), count, type, rank, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Ibcast(const std::vector<T>& buffer, const Datatype& type) const
        {
            return Ibcast(&buffer.front(), buffer.size(), type);
        }

        /*
         * MPI_Igather non-root
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
         * MPI_Igather root
         */

        template <typename T>
        Request Igather(const T* sendbuf, T* recvbuf, MPI_Int count) const
        {
            return Igather(sendbuf, recvbuf, count, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Igather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf) const
        {
            return Igather(sendbuf, recvbuf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Igather(const T* sendbuf, T* recvbuf, MPI_Int count, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Igather(sendbuf, count, type, recvbuf, count, type, rank, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Igather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const Datatype& type) const
        {
            MPIWRAP_ASSERT(size*sendbuf.size() == recvbuf.size(),
                           "Receive buffers size must be equal communicator size times send buffer size.");
            return Igather(&sendbuf.front(), &recvbuf.front(), sendbuf.size(), type);
        }

        /*
         * MPI_Igather root in-place
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
            MPIWRAP_CALL(MPI_Igather(MPI_IN_PLACE, 0, type, recvbuf, count, type, rank, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Igather(std::vector<T>& recvbuf, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvbuf.size()%size == 0,
                           "Receive buffer size must be a multiple of the communicator size.");
            return Igather(&recvbuf.front(), recvbuf.size()/size, type);
        }

        /*
         * MPI_Igatherv non-root
         */

        template <typename T>
        Request Igatherv(const T* sendbuf, MPI_Int sendcount, MPI_Int root) const
        {
            return Igatherv(sendbuf, sendcount, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Igatherv(const std::vector<T>& sendbuf, MPI_Int root) const
        {
            return Igatherv(sendbuf, root, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Igatherv(const T* sendbuf, MPI_Int sendcount, MPI_Int root, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Igatherv(sendbuf, sendcount, type, NULL, NULL, NULL, type, root, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Igatherv(const std::vector<T>& sendbuf, MPI_Int root, const Datatype& type) const
        {
            return Igatherv(&sendbuf.front(), sendbuf.size(), root, type);
        }

        /*
         * MPI_Igatherv root
         */

        template <typename T>
        Request Igather(const T* sendbuf, MPI_Int sendcount, T* recvbuf, const MPI_Int* recvcounts) const
        {
            return Igather(sendbuf, sendcount, recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Igather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf,
                     const std::vector<MPI_Int>& recvcounts) const
        {
            return Igather(sendbuf, recvbuf, recvcounts, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Igather(const T* sendbuf, MPI_Int sendcount, T* recvbuf, const MPI_Int* recvcounts,
                     const Datatype& type) const
        {
            std::vector<MPI_Int> recvdispls = displacements(recvcounts);
            return Igather(sendbuf, recvbuf, recvcounts, &recvdispls.front(), type);
        }

        template <typename T>
        Request Igather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts,
                     const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvbuf.size() == sum(recvcounts),
                           "The receive buffer size must equal the sum of the receive counts.");
            return Igather(sendbuf, recvbuf, recvcounts, displacements(recvcounts), type);
        }

        template <typename T>
        Request Igather(const T* sendbuf, MPI_Int sendcount, T* recvbuf, const MPI_Int* recvcounts, const MPI_Int* recvdispls) const
        {
            return Igather(sendbuf, sendcount, recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Igather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf,
                     const std::vector<MPI_Int>& recvcounts, const std::vector<MPI_Int>& recvdispls) const
        {
            return Igather(sendbuf, recvbuf, recvcounts, recvdispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Igather(const T* sendbuf, MPI_Int sendcount, T* recvbuf, const MPI_Int* recvcounts, const MPI_Int* recvdispls,
                     const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Igatherv(sendbuf, sendcount, type, recvbuf, recvcounts, recvdispls, type, rank, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Igather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts,
                     const std::vector<MPI_Int>& recvdispls, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvcounts.size() == size,
                           "There must be exactly one receive count for each process.");
            MPIWRAP_ASSERT(recvdispls.size() == size,
                           "There must be exactly one receive displacement for each process.");
            MPIWRAP_ASSERT(recvbuf.size() >= recvdispls[size-1]+recvcounts[size-1],
                           "The receive buffer size must be at least as large the sum of last receive count and the last receive displacement.");
            MPIWRAP_ASSERT(sendbuf.size() == recvcounts[rank],
                           "The send buffer size must equal the receive count for this process.");
            return Igather(&sendbuf.front(), sendbuf.size(), &recvbuf.front(), &recvcounts.front(), &recvdispls.front(), type);
        }

        /*
         * MPI_Igatherv root in-place
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
            return Igather(recvbuf, recvcounts, recvdispls, type);
        }

        template <typename T>
        Request Igather(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts, const Datatype& type) const
        {
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
            MPIWRAP_CALL(MPI_Igatherv(MPI_IN_PLACE, 0, type, recvbuf, recvcounts, recvdispls, type, rank, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Igather(std::vector<T>& recvbuf, const std::vector<MPI_Int>& recvcounts, const std::vector<MPI_Int>& recvdispls, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvcounts.size() == size,
                           "There must be exactly one receive count for each process.");
            MPIWRAP_ASSERT(recvdispls.size() == size,
                           "There must be exactly one receive displacement for each process.");
            MPIWRAP_ASSERT(recvbuf.size() >= recvdispls[size-1]+recvcounts[size-1],
                           "The receive buffer size must be at least as large the sum of last receive count and the last receive displacement.");
            return Igather(&recvbuf.front(), &recvcounts.front(), &recvdispls.front(), type);
        }

        /*
         * MPI_Ireduce non-root
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
         * MPI_Ireduce root
         */

        template <typename T>
        Request Ireduce(const T* sendbuf, T* recvbuf, MPI_Int count, const MPI_Op& op) const
        {
            return Ireduce(sendbuf, recvbuf, count, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ireduce(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI_Op& op) const
        {
            return Ireduce(sendbuf, recvbuf, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ireduce(const T* sendbuf, T* recvbuf, MPI_Int count, const MPI_Op& op, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Ireduce(sendbuf, recvbuf, count, type, op, rank, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Ireduce(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendbuf.size() == recvbuf.size(),
                           "Send and receive buffers must be the same size.");
            return Ireduce(&sendbuf.front(), &recvbuf.front(), sendbuf.size(), op, type);
        }

        /*
         * MPI_Ireduce root in-place
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
            MPIWRAP_CALL(MPI_Ireduce(MPI_IN_PLACE, recvbuf, count, type, op, rank, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Ireduce(std::vector<T>& recvbuf, const MPI_Op& op, const Datatype& type) const
        {
            return Ireduce(&recvbuf.front(), recvbuf.size(), op, type);
        }

        /*
         * MPI_Ireduce_scatter in-place
         */

        template <typename T>
        Request Ireduce_scatter(T* recvbuf, MPI_Int* recvcounts, const MPI_Op& op) const
        {
            return Ireduce_scatter(recvbuf, recvcounts, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ireduce_scatter(std::vector<T>& recvbuf, std::vector<MPI_Int>& recvcounts, const MPI_Op& op) const
        {
            return Ireduce_scatter(recvbuf, recvcounts, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ireduce_scatter(T* recvbuf, MPI_Int* recvcounts, const MPI_Op& op, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Ireduce_scatter(MPI_IN_PLACE, recvbuf, recvcounts, type, op, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Ireduce_scatter(std::vector<T>& recvbuf, std::vector<MPI_Int>& recvcounts, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvcounts.size() == size,
                           "There must be exactly one receive count for each process.");
            MPIWRAP_ASSERT(recvbuf.size() == sum(recvcounts),
                           "The send/receive buffer size must equal the sum of the receive counts.");
            return Ireduce_scatter(&recvbuf.front(), &recvcounts.front(), op, type);
        }

        /*
         * MPI_Ireduce_scatter_block in-place
         */

        template <typename T>
        Request Ireduce_scatter(T* recvbuf, MPI_Int count, const MPI_Op& op) const
        {
            return Ireduce_scatter(recvbuf, count, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ireduce_scatter(std::vector<T>& recvbuf, const MPI_Op& op) const
        {
            return Ireduce_scatter(recvbuf, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Ireduce_scatter(T* recvbuf, MPI_Int count, const MPI_Op& op, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Ireduce_scatter(MPI_IN_PLACE, recvbuf, count, type, op, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Ireduce_scatter(std::vector<T>& recvbuf, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_ASSERT(recvbuf.size()%npeers == 0, "Receive buffer size must be a multiple of the number of processes.");
            return Ireduce_scatter(&recvbuf.front(), recvbuf.size()/npeers, op, type);
        }

        /*
         * MPI_Scatter non-root
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
         * MPI_Iscatter root
         */

        template <typename T>
        Request Iscatter(const T* sendbuf, T* recvbuf, MPI_Int recvcount) const
        {
            return Iscatter(sendbuf, recvbuf, recvcount, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iscatter(const std::vector<T>& sendbuf, std::vector<T>& recvbuf) const
        {
            return Iscatter(sendbuf, recvbuf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iscatter(const T* sendbuf, T* recvbuf, MPI_Int recvcount, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Iscatter(sendbuf, recvcount, type, recvbuf, recvcount, type, rank, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Iscatter(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const Datatype& type) const
        {
            MPIWRAP_ASSERT(size*recvbuf.size() == sendbuf.size(),
                           "Send buffer size must be equal communicator size times send receive size.");
            return Iscatter(&sendbuf.front(), &recvbuf.front(), recvbuf.size(), type);
        }

        /*
         * MPI_Iscatter root in-place
         */

        template <typename T>
        Request Iscatter(const T* sendbuf, MPI_Int sendcount) const
        {
            return Iscatter(sendbuf, sendcount, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iscatter(const std::vector<T>& sendbuf) const
        {
            return Iscatter(sendbuf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iscatter(const T* sendbuf, MPI_Int count, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Iscatter(sendbuf, count, type, MPI_IN_PLACE, 0, type, rank, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Iscatter(const std::vector<T>& sendbuf, const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendbuf.size()%size == 0,
                           "Send buffer size must be a multiple of the communicator size.");
            return Iscatter(&sendbuf.front(), sendbuf.size()/size, type);
        }

        /*
         * MPI_Iscatterv non-root
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

        /*
         * MPI_Iscatterv root
         */

        template <typename T>
        Request Iscatter(const T* sendbuf, const MPI_Int* sendcounts, T* recvbuf, MPI_Int recvcount) const
        {
            return Iscatter(sendbuf, sendcounts, recvbuf, recvcount, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iscatter(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts,
                      std::vector<T>& recvbuf) const
        {
            return Iscatter(sendbuf, sendcounts, recvbuf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iscatter(const T* sendbuf, const MPI_Int* sendcounts, T* recvbuf, MPI_Int recvcount,
                      const Datatype& type) const
        {
            std::vector<MPI_Int> senddispls = displacements(sendcounts);
            return Iscatter(sendbuf, sendcounts, &senddispls.front(), recvbuf, recvcount, type);
        }

        template <typename T>
        Request Iscatter(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts,
                      std::vector<T>& recvbuf, const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendbuf.size() == sum(sendcounts),
                           "Send buffer size must equal the sum of the send counts.");
            return Iscatter(sendbuf, sendcounts, displacements(sendcounts), recvbuf, type);
        }

        template <typename T>
        Request Iscatter(const T* sendbuf, const MPI_Int* sendcounts, const MPI_Int* senddispls, T* recvbuf, MPI_Int recvcount) const
        {
            return Iscatter(sendbuf, sendcounts, senddispls, recvbuf, recvcount, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iscatter(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts,
                      const std::vector<MPI_Int>& senddispls, std::vector<T>& recvbuf) const
        {
            return Iscatter(sendbuf, sendcounts, senddispls, recvbuf, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iscatter(const T* sendbuf, const MPI_Int* sendcounts, const MPI_Int* senddispls, T* recvbuf, MPI_Int recvcount,
                      const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Iscatterv(sendbuf, sendcounts, senddispls, type, recvbuf, recvcount, type, rank, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Iscatter(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts,
                      const std::vector<MPI_Int>& senddispls, std::vector<T>& recvbuf, const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendcounts.size() == size,
                           "There must be exactly one send count for each process.");
            MPIWRAP_ASSERT(senddispls.size() == size,
                           "There must be exactly one send displacement for each process.");
            MPIWRAP_ASSERT(recvbuf.size() == sendcounts[rank],
                           "Receive buffer size must equal send count for this process.");
            MPIWRAP_ASSERT(sendbuf.size() >= senddispls[size-1]+sendcounts[size-1],
                           "Send buffer size must be at least the sum of the last send count and the last send displacement.");
            return Iscatter(&sendbuf.front(), &sendcounts.front(), &senddispls.front(), &recvbuf.front(), recvbuf.size(), type);
        }

        /*
         * MPI_Scatterv root in-place
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
        Request Iscatter(const T* sendbuf, const MPI_Int* sendcounts, const Datatype& type) const
        {
            std::vector<MPI_Int> senddispls = displacements(sendcounts);
            return Iscatter(sendbuf, sendcounts, &senddispls.front(), type);
        }

        template <typename T>
        Request Iscatter(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts,
                     const Datatype& type) const
        {
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
        Request Iscatter(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts, const std::vector<MPI_Int>& senddispls) const
        {
            return Iscatter(sendbuf, sendcounts, senddispls, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iscatter(const T* sendbuf, const MPI_Int* sendcounts, const MPI_Int* senddispls, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Iscatterv(sendbuf, sendcounts, senddispls, type, MPI_IN_PLACE, 0, type, rank, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Iscatter(const std::vector<T>& sendbuf, const std::vector<MPI_Int>& sendcounts, const std::vector<MPI_Int>& senddispls,
                      const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendcounts.size() == size,
                           "There must be exactly one send count for each process.");
            MPIWRAP_ASSERT(senddispls.size() == size,
                           "There must be exactly one send displacement for each process.");
            MPIWRAP_ASSERT(sendbuf.size() >= senddispls[size-1]+sendcounts[size-1],
                           "Send buffer size must be at least the sum of the last send count and the last send displacement.");
            return Iscatter(&sendbuf.front(), &sendcounts.front(), &senddispls.front(), type);
        }

        /*
         * MPI_IExscan
         */

        template <typename T>
        Request IExscan(const T* sendbuf, T* recvbuf, MPI_Int count, const MPI_Op& op) const
        {
            return IExscan(sendbuf, recvbuf, count, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request IExscan(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI_Op& op) const
        {
            return IExscan(sendbuf, recvbuf, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request IExscan(const T* sendbuf, T* recvbuf, MPI_Int count, const MPI_Op& op, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Iexscan(sendbuf, recvbuf, count, type, op, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request IExscan(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendbuf.size() == recvbuf.size(),
                           "Send and receive buffers must be the same size.");
            return IExscan(&sendbuf.front(), &recvbuf.begin(), sendbuf.size(), op, type);
        }

        /*
         * MPI_IExscan in-place
         */

        template <typename T>
        Request IExscan(T* recvbuf, MPI_Int count, const MPI_Op& op) const
        {
            return IExscan(recvbuf, count, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request IExscan(std::vector<T>& recvbuf, const MPI_Op& op) const
        {
            return IExscan(recvbuf, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request IExscan(T* recvbuf, MPI_Int count, const MPI_Op& op, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Iexscan(MPI_IN_PLACE, recvbuf, count, type, op, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request IExscan(std::vector<T>& recvbuf, const MPI_Op& op, const Datatype& type) const
        {
            return IExscan(&recvbuf.front(), recvbuf.size(), op, type);
        }

        /*
         * MPI_Iscan
         */

        template <typename T>
        Request Iscan(const T* sendbuf, T* recvbuf, MPI_Int count, const MPI_Op& op) const
        {
            return Iscan(sendbuf, recvbuf, count, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iscan(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI_Op& op) const
        {
            return Iscan(sendbuf, recvbuf, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iscan(const T* sendbuf, T* recvbuf, MPI_Int count, const MPI_Op& op, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Iscan(sendbuf, recvbuf, count, type, op, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Iscan(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI_Op& op, const Datatype& type) const
        {
            MPIWRAP_ASSERT(sendbuf.size() == recvbuf.size(),
                           "Send and receive buffers must be the same size.");
            return Iscan(&sendbuf.front(), &recvbuf.front(), sendbuf.size(), op, type);
        }

        /*
         * MPI_Iscan in-place
         */

        template <typename T>
        Request Iscan(T* recvbuf, MPI_Int count, const MPI_Op& op) const
        {
            return Iscan(recvbuf, count, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iscan(std::vector<T>& recvbuf, const MPI_Op& op) const
        {
            return Iscan(recvbuf, op, MPI_TYPE_<T>::value());
        }

        template <typename T>
        Request Iscan(T* recvbuf, MPI_Int count, const MPI_Op& op, const Datatype& type) const
        {
            MPI_Request req;
            MPIWRAP_CALL(MPI_Iscan(MPI_IN_PLACE, recvbuf, count, type, op, comm, &req));
            return Request(req);
        }

        template <typename T>
        Request Iscan(std::vector<T>& recvbuf, const MPI_Op& op, const Datatype& type) const
        {
            return Iscan(&recvbuf.front(), recvbuf.size(), op, type);
        }

#endif
};

}

#endif

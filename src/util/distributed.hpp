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

#ifndef _AQUARIUS_UTIL_DISTRIBUTED_HPP_
#define _AQUARIUS_UTIL_DISTRIBUTED_HPP_

#include <iostream>
#include <complex>
#include <fstream>

#include "mpi.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#ifdef ELEMENTAL
#include "El.hpp"
#endif

#include "ctf.hpp"
#include "util/stl_ext.hpp"

namespace aquarius
{

namespace tensor
{

template <typename T> class CTFTensor;
template <typename T> class SymmetryBlockedTensor;
template <typename T> class SpinorbitalTensor;

}

template <typename T>
struct MPI_TYPE_ {};

template <>
struct MPI_TYPE_<float>
{
    static const MPI::Datatype& value() { return MPI::FLOAT; }
};

template <>
struct MPI_TYPE_<double>
{
    static const MPI::Datatype& value() { return MPI::DOUBLE; }
};

template <>
struct MPI_TYPE_< std::complex<float> >
{
    static const MPI::Datatype& value() { return MPI::COMPLEX; }
};

template <>
struct MPI_TYPE_< std::complex<double> >
{
    static const MPI::Datatype& value() { return MPI::DOUBLE_COMPLEX; }
};

template <>
struct MPI_TYPE_<short>
{
    static const MPI::Datatype& value() { return MPI::SHORT; }
};

template <>
struct MPI_TYPE_<unsigned short>
{
    static const MPI::Datatype& value() { return MPI::UNSIGNED_SHORT; }
};

template <>
struct MPI_TYPE_<int>
{
    static const MPI::Datatype& value() { return MPI::INT; }
};

template <>
struct MPI_TYPE_<unsigned int>
{
    static const MPI::Datatype& value() { return MPI::UNSIGNED; }
};

template <>
struct MPI_TYPE_<long>
{
    static const MPI::Datatype& value() { return MPI::LONG; }
};

template <>
struct MPI_TYPE_<unsigned long>
{
    static const MPI::Datatype& value() { return MPI::UNSIGNED_LONG; }
};

template <>
struct MPI_TYPE_<long long>
{
    static const MPI::Datatype& value() { return MPI::LONG_LONG; }
};

template <>
struct MPI_TYPE_<unsigned long long>
{
    static const MPI::Datatype& value() { return MPI::UNSIGNED_LONG_LONG; }
};

class Arena
{
    protected:
        MPI::Intracomm comm;
        //std::global_ptr<tCTF_World<float> > ctfs;
        std::global_ptr<tCTF_World<double> > ctfd;
        //std::global_ptr<tCTF_World<std::complex<float> > > ctfc;
        //std::global_ptr<tCTF_World<std::complex<double> > > ctfz;

    public:
        const int rank;
        const int nproc;

        Arena(MPI::Intracomm& comm = MPI::COMM_WORLD)
        : comm(comm), rank(comm.Get_rank()), nproc(comm.Get_size()) {}

        const MPI::Intracomm& getCommunicator() const { return comm; }

        template <typename T>
        tCTF_World<T>& ctf();

        template <typename T>
        const tCTF_World<T>& ctf() const
        {
            return const_cast<const tCTF_World<T>&>(const_cast<Arena&>(*this).ctf<T>());
        }

        template <typename T>
        void Allgather(const T* sendbuf, T* recvbuf, int count) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Allgather(sendbuf, count, type, recvbuf, count, type);
        }

        template <typename T>
        void Allgather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Allgather(sendbuf.data(), sendbuf.size(), type, recvbuf.data(), sendbuf.size(), type);
        }

        template <typename T>
        void Allgather(const T* sendbuf, T* recvbuf, int count, const MPI::Datatype& type) const
        {
            comm.Allgather(sendbuf, count, type, recvbuf, count, type);
        }

        template <typename T>
        void Allgather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI::Datatype& type) const
        {
            comm.Allgather(sendbuf.data(), sendbuf.size(), type, recvbuf.data(), sendbuf.size(), type);
        }

        template <typename T>
        void Allgather(T* recvbuf, int count) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Allgather(MPI::IN_PLACE, 0, type, recvbuf, count, type);
        }

        template <typename T>
        void Allgather(std::vector<T>& recvbuf) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Allgather(MPI::IN_PLACE, 0, type, recvbuf.data(), recvbuf.size()/nproc, type);
        }

        template <typename T>
        void Allgather(T* recvbuf, int count, const MPI::Datatype& type) const
        {
            comm.Allgather(MPI::IN_PLACE, 0, type, recvbuf, count, type);
        }

        template <typename T>
        void Allgather(std::vector<T>& recvbuf, const MPI::Datatype& type) const
        {
            comm.Allgather(MPI::IN_PLACE, 0, type, recvbuf.data(), recvbuf.size()/nproc, type);
        }

        template <typename T>
        void Allgatherv(const T* sendbuf, int sendcount, T* recvbuf, const int* recvcounts) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+recvcounts[i-1];
            comm.Allgatherv(sendbuf, sendcount, type, recvbuf, recvcounts, displs.data(), type);
        }

        template <typename T>
        void Allgatherv(const std::vector<T>& sendbuf, std::vector<T>& recvbuf,
                        const std::vector<int>& recvcounts) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+recvcounts[i-1];
            comm.Allgatherv(sendbuf.data(), sendbuf.size(), type, recvbuf.data(), recvcounts.data(), displs.data(), type);
        }

        template <typename T>
        void Allgatherv(const T* sendbuf, int sendcount, T* recvbuf, const int* recvcounts,
                        const MPI::Datatype& type) const
        {
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+recvcounts[i-1];
            comm.Allgatherv(sendbuf, sendcount, type, recvbuf, recvcounts, displs.data(), type);
        }

        template <typename T>
        void Allgatherv(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const std::vector<int>& recvcounts,
                        const MPI::Datatype& type) const
        {
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+recvcounts[i-1];
            comm.Allgatherv(sendbuf.data(), sendbuf.size(), type, recvbuf.data(), recvcounts.data(), displs.data(), type);
        }

        template <typename T>
        void Allgatherv(T* recvbuf, const int* recvcounts) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+recvcounts[i-1];
            comm.Allgatherv(MPI::IN_PLACE, 0, type, recvbuf, recvcounts, displs.data(), type);
        }

        template <typename T>
        void Allgatherv(std::vector<T>& recvbuf, const std::vector<int>& recvcounts) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+recvcounts[i-1];
            comm.Allgatherv(MPI::IN_PLACE, 0, type, recvbuf.data(), recvcounts.data(), displs.data(), type);
        }

        template <typename T>
        void Allgatherv(T* recvbuf, const int* recvcounts, const MPI::Datatype& type) const
        {
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+recvcounts[i-1];
            comm.Allgatherv(MPI::IN_PLACE, 0, type, recvbuf, recvcounts, displs.data(), type);
        }

        template <typename T>
        void Allgatherv(std::vector<T>& recvbuf, const std::vector<int>& recvcounts, const MPI::Datatype& type) const
        {
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+recvcounts[i-1];
            comm.Allgatherv(MPI::IN_PLACE, 0, type, recvbuf.data(), recvcounts.data(), displs.data(), type);
        }

        template <typename T>
        void Allreduce(const T* sendbuf, T* recvbuf, int count, const MPI::Op& op) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Allreduce(sendbuf, recvbuf, count, type, op);
        }

        template <typename T>
        void Allreduce(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI::Op& op) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Allreduce(sendbuf.data(), recvbuf.data(), sendbuf.size(), type, op);
        }

        template <typename T>
        void Allreduce(const T* sendbuf, T* recvbuf, int count, const MPI::Op& op, const MPI::Datatype& type) const
        {
            comm.Allreduce(sendbuf, recvbuf, count, type, op);
        }

        template <typename T>
        void Allreduce(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI::Op& op,
                       const MPI::Datatype& type) const
        {
            comm.Allreduce(sendbuf.data(), recvbuf.data(), sendbuf.size(), type, op);
        }

        template <typename T>
        void Allreduce(T* buf, int count, const MPI::Op& op) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Allreduce(MPI::IN_PLACE, buf, count, type, op);
        }

        template <typename T>
        void Allreduce(std::vector<T>& buf, const MPI::Op& op) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Allreduce(MPI::IN_PLACE, buf.data(), buf.size(), type, op);
        }

        template <typename T>
        void Allreduce(T* buf, int count, const MPI::Op& op, const MPI::Datatype& type) const
        {
            comm.Allreduce(MPI::IN_PLACE, buf, count, type, op);
        }

        template <typename T>
        void Allreduce(std::vector<T>& buf, const MPI::Op& op, const MPI::Datatype& type) const
        {
            comm.Allreduce(MPI::IN_PLACE, buf.data(), buf.size(), type, op);
        }

        template <typename T>
        void Alltoall(const T* sendbuf, T* recvbuf, int count) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Alltoall(sendbuf, count, type, recvbuf, count, type);
        }

        template <typename T>
        void Alltoall(const std::vector<T>& sendbuf, std::vector<T>& recvbuf) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Alltoall(sendbuf.data(), sendbuf.size()/nproc, type, recvbuf.data(), recvbuf.size()/nproc, type);
        }

        template <typename T>
        void Alltoall(const T* sendbuf, T* recvbuf, int count, const MPI::Datatype& type) const
        {
            comm.Alltoall(sendbuf, count, type, recvbuf, count, type);
        }

        template <typename T>
        void Alltoall(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI::Datatype& type) const
        {
            comm.Alltoall(sendbuf.data(), sendbuf.size()/nproc, type, recvbuf.data(), recvbuf.size()/nproc, type);
        }

        template <typename T>
        void Alltoallv(const T* sendbuf, const int* sendcounts, T* recvbuf, const int* recvcounts) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            std::vector<int> sdispls(nproc);
            sdispls[0] = 0;
            for (int i = 1;i < nproc;i++) sdispls[i] = sdispls[i-1]+sendcounts[i-1];
            std::vector<int> rdispls(nproc);
            rdispls[0] = 0;
            for (int i = 1;i < nproc;i++) rdispls[i] = rdispls[i-1]+recvcounts[i-1];
            comm.Alltoallv(sendbuf, sendcounts, sdispls.data(), type,
                           recvbuf, recvcounts, rdispls.data(), type);
        }

        template <typename T>
        void Alltoallv(const std::vector<T>& sendbuf, const std::vector<int>& sendcounts,
                             std::vector<T>& recvbuf, const std::vector<int>& recvcounts) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            std::vector<int> sdispls(nproc);
            sdispls[0] = 0;
            for (int i = 1;i < nproc;i++) sdispls[i] = sdispls[i-1]+sendcounts[i-1];
            std::vector<int> rdispls(nproc);
            rdispls[0] = 0;
            for (int i = 1;i < nproc;i++) rdispls[i] = rdispls[i-1]+recvcounts[i-1];
            comm.Alltoallv(sendbuf.data(), sendcounts.data(), sdispls.data(), type,
                           recvbuf.data(), recvcounts.data(), rdispls.data(), type);
        }

        template <typename T>
        void Alltoallv(const T* sendbuf, const int* sendcounts,
                             T* recvbuf, const int* recvcounts, const MPI::Datatype& type) const
        {
            std::vector<int> sdispls(nproc);
            sdispls[0] = 0;
            for (int i = 1;i < nproc;i++) sdispls[i] = sdispls[i-1]+sendcounts[i-1];
            std::vector<int> rdispls(nproc);
            rdispls[0] = 0;
            for (int i = 1;i < nproc;i++) rdispls[i] = rdispls[i-1]+recvcounts[i-1];
            comm.Alltoallv(sendbuf, sendcounts, sdispls.data(), type,
                           recvbuf, recvcounts, rdispls.data(), type);
        }

        template <typename T>
        void Alltoallv(const std::vector<T>& sendbuf, const std::vector<int>& sendcounts,
                             std::vector<T>& recvbuf, const std::vector<int>& recvcounts,
                       const MPI::Datatype& type) const
        {
            std::vector<int> sdispls(nproc);
            sdispls[0] = 0;
            for (int i = 1;i < nproc;i++) sdispls[i] = sdispls[i-1]+sendcounts[i-1];
            std::vector<int> rdispls(nproc);
            rdispls[0] = 0;
            for (int i = 1;i < nproc;i++) rdispls[i] = rdispls[i-1]+recvcounts[i-1];
            comm.Alltoallv(sendbuf.data(), sendcounts.data(), sdispls.data(), type,
                           recvbuf.data(), recvcounts.data(), rdispls.data(), type);
        }

        void Barrier() const
        {
            comm.Barrier();
        }

        template <typename T>
        void Bcast(T* buffer, int count, int root) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Bcast(buffer, count, type, root);
        }

        template <typename T>
        void Bcast(std::vector<T>& buffer, int root) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Bcast(buffer.data(), buffer.size(), type, root);
        }

        template <typename T>
        void Bcast(T* buffer, int count, int root, const MPI::Datatype& type) const
        {
            comm.Bcast(buffer, count, type, root);
        }

        template <typename T>
        void Bcast(std::vector<T>& buffer, int root, const MPI::Datatype& type) const
        {
            comm.Bcast(buffer.data(), buffer.size(), type, root);
        }

        template <typename T>
        void Exscan(const T* sendbuf, T* recvbuf, int count, const MPI::Op& op) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Exscan(sendbuf, recvbuf, count, type, op);
        }

        template <typename T>
        void Exscan(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI::Op& op) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Exscan(sendbuf.data(), recvbuf.data(), sendbuf.size(), type, op);
        }

        template <typename T>
        void Gather(const T* sendbuf, int count, int root) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Gather(sendbuf, count, type, NULL, 0, type, root);
        }

        template <typename T>
        void Gather(const T* sendbuf, T* recvbuf, int count) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Gather(sendbuf, count, type, recvbuf, count, type, rank);
        }

        template <typename T>
        void Gather(T* recvbuf, int count) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Gather(MPI::IN_PLACE, 0, type, recvbuf, count, type, rank);
        }

        template <typename T>
        void Gather(const std::vector<T>& sendbuf, int root) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Gather(sendbuf.data(), sendbuf.size(), type, NULL, 0, type, root);
        }

        template <typename T>
        void Gather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Gather(sendbuf.data(), sendbuf.size(), type, recvbuf.data(), sendbuf.size(), type, rank);
        }

        template <typename T>
        void Gather(std::vector<T>& recvbuf) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Gather(MPI::IN_PLACE, 0, type, recvbuf.data(), recvbuf.size()/nproc, type, rank);
        }

        template <typename T>
        void Gather(const T* sendbuf, int count, int root, const MPI::Datatype& type) const
        {
            comm.Gather(sendbuf, count, type, NULL, 0, type, root);
        }

        template <typename T>
        void Gather(const T* sendbuf, T* recvbuf, int count, const MPI::Datatype& type) const
        {
            comm.Gather(sendbuf, count, type, recvbuf, count, type, rank);
        }

        template <typename T>
        void Gather(T* recvbuf, int count, const MPI::Datatype& type) const
        {
            comm.Gather(MPI::IN_PLACE, 0, type, recvbuf, count, type, rank);
        }

        template <typename T>
        void Gather(const std::vector<T>& sendbuf, int root, const MPI::Datatype& type) const
        {
            comm.Gather(sendbuf.data(), sendbuf.size(), type, NULL, 0, type, root);
        }

        template <typename T>
        void Gather(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI::Datatype& type) const
        {
            comm.Gather(sendbuf.data(), sendbuf.size(), type, recvbuf.data(), sendbuf.size(), type, rank);
        }

        template <typename T>
        void Gather(std::vector<T>& recvbuf, const MPI::Datatype& type) const
        {
            comm.Gather(MPI::IN_PLACE, 0, type, recvbuf.data(), recvbuf.size()/nproc, type, rank);
        }

        template <typename T>
        void Gatherv(const T* sendbuf, int sendcount, int root) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Gatherv(sendbuf, sendcount, type, NULL, NULL, NULL, type, root);
        }

        template <typename T>
        void Gatherv(const T* sendbuf, int sendcount, T* recvbuf, const int* recvcounts) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+recvcounts[i-1];
            comm.Gatherv(sendbuf, sendcount, type, recvbuf, recvcounts, displs.data(), type, rank);
        }

        template <typename T>
        void Gatherv(T* recvbuf, const int* recvcounts) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+recvcounts[i-1];
            comm.Gatherv(MPI::IN_PLACE, 0, type, recvbuf, recvcounts, displs.data(), type, rank);
        }

        template <typename T>
        void Gatherv(const std::vector<T>& sendbuf, int root) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Gatherv(sendbuf.data(), sendbuf.size(), type, NULL, NULL, NULL, type, root);
        }

        template <typename T>
        void Gatherv(const std::vector<T>& sendbuf, std::vector<T>& recvbuf,
                     const std::vector<int>& recvcounts) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+recvcounts[i-1];
            comm.Gatherv(sendbuf.data(), sendbuf.size(), type, recvbuf.data(), recvcounts.data(), displs.data(), type, rank);
        }

        template <typename T>
        void Gatherv(std::vector<T>& recvbuf, const std::vector<int>& recvcounts) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+recvcounts[i-1];
            comm.Gatherv(MPI::IN_PLACE, 0, type, recvbuf.data(), recvcounts.data(), displs.data(), type, rank);
        }

        template <typename T>
        void Gatherv(const T* sendbuf, int sendcount, int root, const MPI::Datatype& type) const
        {
            comm.Gatherv(sendbuf, sendcount, type, NULL, NULL, NULL, type, root);
        }

        template <typename T>
        void Gatherv(const T* sendbuf, int sendcount, T* recvbuf, const int* recvcounts,
                     const MPI::Datatype& type) const
        {
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+recvcounts[i-1];
            comm.Gatherv(sendbuf, sendcount, type, recvbuf, recvcounts, displs.data(), type, rank);
        }

        template <typename T>
        void Gatherv(T* recvbuf, const int* recvcounts, const MPI::Datatype& type) const
        {
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+recvcounts[i-1];
            comm.Gatherv(MPI::IN_PLACE, 0, type, recvbuf, recvcounts, displs.data(), type, rank);
        }

        template <typename T>
        void Gatherv(const std::vector<T>& sendbuf, int root, const MPI::Datatype& type) const
        {
            comm.Gatherv(sendbuf.data(), sendbuf.size(), type, NULL, NULL, NULL, type, root);
        }

        template <typename T>
        void Gatherv(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const std::vector<int>& recvcounts,
                     const MPI::Datatype& type) const
        {
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+recvcounts[i-1];
            comm.Gatherv(sendbuf.data(), sendbuf.size(), type, recvbuf.data(), recvcounts.data(), displs.data(), type, rank);
        }

        template <typename T>
        void Gatherv(std::vector<T>& recvbuf, const std::vector<int>& recvcounts, const MPI::Datatype& type) const
        {
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+recvcounts[i-1];
            comm.Gatherv(MPI::IN_PLACE, 0, type, recvbuf.data(), recvcounts.data(), displs.data(), type, rank);
        }

        bool Iprobe(int source, int tag, MPI::Status& status) const
        {
            return comm.Iprobe(source, tag, status);
        }

        bool Iprobe(int source, int tag) const
        {
            return comm.Iprobe(source, tag);
        }

        template <typename T>
        MPI::Request Irecv(T* buf, int count, int source, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            return comm.Irecv(buf, count, type, source, tag);
        }

        template <typename T>
        MPI::Request Irecv(std::vector<T>& buf, int source, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            return comm.Irecv(buf.data(), buf.size(), type, source, tag);
        }

        template <typename T>
        MPI::Request Irecv(T* buf, int count, int source, int tag, const MPI::Datatype& type) const
        {
            return comm.Irecv(buf, count, type, source, tag);
        }

        template <typename T>
        MPI::Request Irecv(std::vector<T>& buf, int source, int tag, const MPI::Datatype& type) const
        {
            return comm.Irecv(buf.data(), buf.size(), type, source, tag);
        }

        template <typename T>
        MPI::Request Irsend(const T* buf, int count, int dest, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            return comm.Irsend(buf, count, type, dest, tag);
        }

        template <typename T>
        MPI::Request Irsend(const std::vector<T>& buf, int dest, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            return comm.Irsend(buf.data(), buf.size(), type, dest, tag);
        }

        template <typename T>
        MPI::Request Irsend(const T* buf, int count, int dest, int tag, const MPI::Datatype& type) const
        {
            return comm.Irsend(buf, count, type, dest, tag);
        }

        template <typename T>
        MPI::Request Irsend(const std::vector<T>& buf, int dest, int tag, const MPI::Datatype& type) const
        {
            return comm.Irsend(buf.data(), buf.size(), type, dest, tag);
        }

        template <typename T>
        MPI::Request Isend(const T* buf, int count, int dest, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            return comm.Isend(buf, count, type, dest, tag);
        }

        template <typename T>
        MPI::Request Isend(const std::vector<T>& buf, int dest, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            return comm.Isend(buf.data(), buf.size(), type, dest, tag);
        }

        template <typename T>
        MPI::Request Isend(const T* buf, int count, int dest, int tag, const MPI::Datatype& type) const
        {
            return comm.Isend(buf, count, type, dest, tag);
        }

        template <typename T>
        MPI::Request Isend(const std::vector<T>& buf, int dest, int tag, const MPI::Datatype& type) const
        {
            return comm.Isend(buf.data(), buf.size(), type, dest, tag);
        }

        template <typename T>
        MPI::Request Issend(const T* buf, int count, int dest, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            return comm.Issend(buf, count, type, dest, tag);
        }

        template <typename T>
        MPI::Request Issend(const std::vector<T>& buf, int dest, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            return comm.Issend(buf.data(), buf.size(), type, dest, tag);
        }

        template <typename T>
        MPI::Request Issend(const T* buf, int count, int dest, int tag, const MPI::Datatype& type) const
        {
            return comm.Issend(buf, count, type, dest, tag);
        }

        template <typename T>
        MPI::Request Issend(const std::vector<T>& buf, int dest, int tag, const MPI::Datatype& type) const
        {
            return comm.Issend(buf.data(), buf.size(), type, dest, tag);
        }

        void Probe(int source, int tag, MPI::Status& status) const
        {
            comm.Probe(source, tag, status);
        }

        void Probe(int source, int tag) const
        {
            comm.Probe(source, tag);
        }

        template <typename T>
        void Recv(T* buf, int count, int source, int tag, MPI::Status& status) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Recv(buf, count, type, source, tag, status);
        }

        template <typename T>
        void Recv(T* buf, int count, int source, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Recv(buf, count, type, source, tag);
        }

        template <typename T>
        void Recv(std::vector<T>& buf, int source, int tag, MPI::Status& status) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Recv(buf.data(), buf.size(), type, source, tag, status);
        }

        template <typename T>
        void Recv(std::vector<T>& buf, int source, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Recv(buf.data(), buf.size(), type, source, tag);
        }

        template <typename T>
        void Recv(T* buf, int count, int source, int tag, MPI::Status& status, const MPI::Datatype& type) const
        {
            comm.Recv(buf, count, type, source, tag, status);
        }

        template <typename T>
        void Recv(T* buf, int count, int source, int tag, const MPI::Datatype& type) const
        {
            comm.Recv(buf, count, type, source, tag);
        }

        template <typename T>
        void Recv(std::vector<T>& buf, int source, int tag, MPI::Status& status, const MPI::Datatype& type) const
        {
            comm.Recv(buf.data(), buf.size(), type, source, tag, status);
        }

        template <typename T>
        void Recv(std::vector<T>& buf, int source, int tag, const MPI::Datatype& type) const
        {
            comm.Recv(buf.data(), buf.size(), type, source, tag);
        }

        template <typename T>
        MPI::Prequest Recv_init(T* buf, int count, int source, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            return comm.Recv_init(buf, count, type, source, tag);
        }

        template <typename T>
        MPI::Prequest Recv_init(std::vector<T>& buf, int source, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            return comm.Recv_init(buf.data(), buf.size(), type, source, tag);
        }

        template <typename T>
        MPI::Prequest Recv_init(T* buf, int count, int source, int tag, const MPI::Datatype& type) const
        {
            return comm.Recv_init(buf, count, type, source, tag);
        }

        template <typename T>
        MPI::Prequest Recv_init(std::vector<T>& buf, int source, int tag, const MPI::Datatype& type) const
        {
            return comm.Recv_init(buf.data(), buf.size(), type, source, tag);
        }

        template <typename T>
        void Reduce(const T* sendbuf, int count, const MPI::Op& op, int root) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Reduce(sendbuf, NULL, count, type, op, root);
        }

        template <typename T>
        void Reduce(const T* sendbuf, T* recvbuf, int count, const MPI::Op& op) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Reduce(sendbuf, recvbuf, count, type, op, rank);
        }

        template <typename T>
        void Reduce(T* recvbuf, int count, const MPI::Op& op) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Reduce(MPI::IN_PLACE, recvbuf, count, type, op, rank);
        }

        template <typename T>
        void Reduce(const std::vector<T>& sendbuf, const MPI::Op& op, int root) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Reduce(sendbuf.data(), NULL, sendbuf.size(), type, op, root);
        }

        template <typename T>
        void Reduce(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI::Op& op) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Reduce(sendbuf.data(), recvbuf.data(), sendbuf.size(), type, op, rank);
        }

        template <typename T>
        void Reduce(std::vector<T>& recvbuf, const MPI::Op& op) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Reduce(MPI::IN_PLACE, recvbuf.data(), recvbuf.size(), type, op, rank);
        }

        template <typename T>
        void Reduce_scatter(const T* sendbuf, T* recvbuf, int* recvcounts, const MPI::Op& op) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Reduce_scatter(sendbuf, recvbuf, recvcounts, type, op);
        }

        template <typename T>
        void Reduce_scatter(const std::vector<T>& sendbuf, std::vector<T>& recvbuf,
                            std::vector<int>& recvcounts, const MPI::Op& op) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Reduce_scatter(sendbuf.data(), recvbuf.data(), recvcounts.data(), type, op);
        }

        template <typename T>
        void Reduce_scatter(T* recvbuf, int* recvcounts, const MPI::Op& op) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Reduce_scatter(MPI::IN_PLACE, recvbuf, recvcounts, type, op);
        }

        template <typename T>
        void Reduce_scatter(std::vector<T>& recvbuf, std::vector<int>& recvcounts, const MPI::Op& op) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Reduce_scatter(MPI::IN_PLACE, recvbuf.data(), recvcounts.data(), type, op);
        }

        template <typename T>
        void Rsend(const T* buf, int count, int dest, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Rsend(buf, count, type, dest, tag);
        }

        template <typename T>
        void Rsend(const std::vector<T>& buf, int dest, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Rsend(buf.data(), buf.size(), type, dest, tag);
        }

        template <typename T>
        void Rsend(const T* buf, int count, int dest, int tag, const MPI::Datatype& type) const
        {
            comm.Rsend(buf, count, type, dest, tag);
        }

        template <typename T>
        void Rsend(const std::vector<T>& buf, int dest, int tag, const MPI::Datatype& type) const
        {
            comm.Rsend(buf.data(), buf.size(), type, dest, tag);
        }

        template <typename T>
        MPI::Prequest Rsend_init(const T* buf, int count, int dest, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            return comm.Rsend_init(buf, count, type, dest, tag);
        }

        template <typename T>
        MPI::Prequest Rsend_init(const std::vector<T>& buf, int dest, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            return comm.Rsend_init(buf.data(), buf.size(), type, dest, tag);
        }

        template <typename T>
        MPI::Prequest Rsend_init(const T* buf, int count, int dest, int tag, const MPI::Datatype& type) const
        {
            return comm.Rsend_init(buf, count, type, dest, tag);
        }

        template <typename T>
        MPI::Prequest Rsend_init(const std::vector<T>& buf, int dest, int tag, const MPI::Datatype& type) const
        {
            return comm.Rsend_init(buf.data(), buf.size(), type, dest, tag);
        }

        template <typename T>
        void Scan(const T* sendbuf, T* recvbuf, int count, const MPI::Op& op) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Scan(sendbuf, recvbuf, count, type, op);
        }

        template <typename T>
        void Scan(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI::Op& op) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Scan(sendbuf.data(), recvbuf.data(), sendbuf.size(), type, op);
        }

        template <typename T>
        void Scatter(T* recvbuf, int recvcount, int root) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Scatter(NULL, 0, type, recvbuf, recvcount, type, root);
        }

        template <typename T>
        void Scatter(const T* sendbuf, int sendcount, T* recvbuf, int recvcount) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Scatter(sendbuf, sendcount, type, recvbuf, recvcount, type, rank);
        }

        template <typename T>
        void Scatter(const T* sendbuf, int sendcount) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Scatter(sendbuf, sendcount, type, MPI::IN_PLACE, 0, type, rank);
        }

        template <typename T>
        void Scatter(std::vector<T>& recvbuf, int root) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Scatter(NULL, 0, type, recvbuf.data(), recvbuf.size(), type, root);
        }

        template <typename T>
        void Scatter(const std::vector<T>& sendbuf, std::vector<T>& recvbuf) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Scatter(sendbuf.data(), recvbuf.size(), type, recvbuf.data(), recvbuf.size(), type, rank);
        }

        template <typename T>
        void Scatter(const std::vector<T>& sendbuf) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Scatter(sendbuf.data(), sendbuf.size()/nproc, type, MPI::IN_PLACE, 0, type, rank);
        }

        template <typename T>
        void Scatter(T* recvbuf, int recvcount, int root, const MPI::Datatype& type) const
        {
            comm.Scatter(NULL, 0, type, recvbuf, recvcount, type, root);
        }

        template <typename T>
        void Scatter(const T* sendbuf, int sendcount, T* recvbuf, int recvcount, const MPI::Datatype& type) const
        {
            comm.Scatter(sendbuf, sendcount, type, recvbuf, recvcount, type, rank);
        }

        template <typename T>
        void Scatter(const T* sendbuf, int sendcount, const MPI::Datatype& type) const
        {
            comm.Scatter(sendbuf, sendcount, type, MPI::IN_PLACE, 0, type, rank);
        }

        template <typename T>
        void Scatter(std::vector<T>& recvbuf, int root, const MPI::Datatype& type) const
        {
            comm.Scatter(NULL, 0, type, recvbuf.data(), recvbuf.size(), type, root);
        }

        template <typename T>
        void Scatter(const std::vector<T>& sendbuf, std::vector<T>& recvbuf, const MPI::Datatype& type) const
        {
            comm.Scatter(sendbuf.data(), recvbuf.size(), type, recvbuf.data(), recvbuf.size(), type, rank);
        }

        template <typename T>
        void Scatter(const std::vector<T>& sendbuf, const MPI::Datatype& type) const
        {
            comm.Scatter(sendbuf.data(), sendbuf.size()/nproc, type, MPI::IN_PLACE, 0, type, rank);
        }

        template <typename T>
        void Scatterv(const T* sendbuf, const int* sendcounts, T* recvbuf, int recvcount, int root) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+sendcounts[i-1];
            comm.Scatterv(NULL, sendcounts, displs.data(), type, recvbuf, recvcount, type, root);
        }

        template <typename T>
        void Scatterv(const T* sendbuf, const int* sendcounts, T* recvbuf, int recvcount) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+sendcounts[i-1];
            comm.Scatterv(sendbuf, sendcounts, displs.data(), type, recvbuf, recvcount, type, rank);
        }

        template <typename T>
        void Scatterv(const T* sendbuf, const int* sendcounts) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+sendcounts[i-1];
            comm.Scatterv(sendbuf, sendcounts, displs.data(), type, MPI::IN_PLACE, 0, type, rank);
        }

        template <typename T>
        void Scatterv(const std::vector<T>& sendbuf, const std::vector<int>& sendcounts,
                      std::vector<T>& recvbuf, int root) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+sendcounts[i-1];
            comm.Scatterv(NULL, sendcounts.data(), displs.data(), type,
                          recvbuf.data(), recvbuf.size(), type, root);
        }

        template <typename T>
        void Scatterv(const std::vector<T>& sendbuf, const std::vector<int>& sendcounts,
                      std::vector<T>& recvbuf) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+sendcounts[i-1];
            comm.Scatterv(sendbuf.data(), sendcounts.data(), displs.data(), type,
                          recvbuf.data(), recvbuf.size(), type, rank);
        }

        template <typename T>
        void Scatterv(const std::vector<T>& sendbuf, const std::vector<int>& sendcounts) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+sendcounts[i-1];
            comm.Scatterv(sendbuf.data(), sendcounts.data(), displs.data(), type,
                          MPI::IN_PLACE, 0, type, rank);
        }

        template <typename T>
        void Scatterv(const T* sendbuf, const int* sendcounts, T* recvbuf, int recvcount, int root,
                      const MPI::Datatype& type) const
        {
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+sendcounts[i-1];
            comm.Scatterv(NULL, sendcounts, displs.data(), type, recvbuf, recvcount, type, root);
        }

        template <typename T>
        void Scatterv(const T* sendbuf, const int* sendcounts, T* recvbuf, int recvcount,
                      const MPI::Datatype& type) const
        {
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+sendcounts[i-1];
            comm.Scatterv(sendbuf, sendcounts, displs.data(), type, recvbuf, recvcount, type, rank);
        }

        template <typename T>
        void Scatterv(const T* sendbuf, const int* sendcounts, const MPI::Datatype& type) const
        {
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+sendcounts[i-1];
            comm.Scatterv(sendbuf, sendcounts, displs.data(), type, MPI::IN_PLACE, 0, type, rank);
        }

        template <typename T>
        void Scatterv(const std::vector<T>& sendbuf, const std::vector<int>& sendcounts,
                      std::vector<T>& recvbuf, int root, const MPI::Datatype& type) const
        {
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+sendcounts[i-1];
            comm.Scatterv(NULL, sendcounts.data(), displs.data(), type,
                          recvbuf.data(), recvbuf.size(), type, root);
        }

        template <typename T>
        void Scatterv(const std::vector<T>& sendbuf, const std::vector<int>& sendcounts,
                      std::vector<T>& recvbuf, const MPI::Datatype& type) const
        {
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+sendcounts[i-1];
            comm.Scatterv(sendbuf.data(), sendcounts.data(), displs.data(), type,
                          recvbuf.data(), recvbuf.size(), type, rank);
        }

        template <typename T>
        void Scatterv(const std::vector<T>& sendbuf, const std::vector<int>& sendcounts,
                      const MPI::Datatype& type) const
        {
            std::vector<int> displs(nproc);
            displs[0] = 0;
            for (int i = 1;i < nproc;i++) displs[i] = displs[i-1]+sendcounts[i-1];
            comm.Scatterv(sendbuf.data(), sendcounts.data(), displs.data(), type,
                          MPI::IN_PLACE, 0, type, rank);
        }

        template <typename T>
        void Send(const T* buf, int count, int dest, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Send(buf, count, type, dest, tag);
        }

        template <typename T>
        void Send(const std::vector<T>& buf, int dest, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Send(buf.data(), buf.size(), type, dest, tag);
        }

        template <typename T>
        void Send(const T* buf, int count, int dest, int tag, const MPI::Datatype& type) const
        {
            comm.Send(buf, count, type, dest, tag);
        }

        template <typename T>
        void Send(const std::vector<T>& buf, int dest, int tag, const MPI::Datatype& type) const
        {
            comm.Send(buf.data(), buf.size(), type, dest, tag);
        }

        template <typename T>
        MPI::Prequest Send_init(const T* buf, int count, int dest, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            return comm.Send_init(buf, count, type, dest, tag);
        }

        template <typename T>
        MPI::Prequest Send_init(const std::vector<T>& buf, int dest, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            return comm.Send_init(buf.data(), buf.size(), type, dest, tag);
        }

        template <typename T>
        MPI::Prequest Send_init(const T* buf, int count, int dest, int tag, const MPI::Datatype& type) const
        {
            return comm.Send_init(buf, count, type, dest, tag);
        }

        template <typename T>
        MPI::Prequest Send_init(const std::vector<T>& buf, int dest, int tag, const MPI::Datatype& type) const
        {
            return comm.Send_init(buf.data(), buf.size(), type, dest, tag);
        }

        template <typename T>
        void Sendrecv(const T* sendbuf, int sendcount, int dest,
                            T* recvbuf, int recvcount, int source, int tag, MPI::Status& status) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Sendrecv(sendbuf, sendcount, type, dest, tag,
                          recvbuf, recvcount, type, source, tag, status);
        }

        template <typename T>
        void Sendrecv(const std::vector<T>& sendbuf, int dest,
                            std::vector<T>& recvbuf, int source, int tag, MPI::Status& status) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Sendrecv(sendbuf.data(), sendbuf.size(), type, dest, tag,
                          recvbuf.data(), recvbuf.size(), type, source, tag, status);
        }

        template <typename T>
        void Sendrecv(const T* sendbuf, int sendcount, int dest,
                            T* recvbuf, int recvcount, int source, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Sendrecv(sendbuf, sendcount, type, dest, tag,
                          recvbuf, recvcount, type, source, tag);
        }

        template <typename T>
        void Sendrecv(const std::vector<T>& sendbuf, int dest,
                            std::vector<T>& recvbuf, int source, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Sendrecv(sendbuf.data(), sendbuf.size(), type, dest, tag,
                          recvbuf.data(), recvbuf.size(), type, source, tag);
        }

        template <typename T>
        void Sendrecv(const T* sendbuf, int sendcount, int dest,
                            T* recvbuf, int recvcount, int source, int tag, MPI::Status& status,
                      const MPI::Datatype& type) const
        {
            comm.Sendrecv(sendbuf, sendcount, type, dest, tag,
                          recvbuf, recvcount, type, source, tag, status);
        }

        template <typename T>
        void Sendrecv(const std::vector<T>& sendbuf, int dest,
                            std::vector<T>& recvbuf, int source, int tag, MPI::Status& status,
                      const MPI::Datatype& type) const
        {
            comm.Sendrecv(sendbuf.data(), sendbuf.size(), type, dest, tag,
                          recvbuf.data(), recvbuf.size(), type, source, tag, status);
        }

        template <typename T>
        void Sendrecv(const T* sendbuf, int sendcount, int dest,
                            T* recvbuf, int recvcount, int source, int tag,
                      const MPI::Datatype& type) const
        {
            comm.Sendrecv(sendbuf, sendcount, type, dest, tag,
                          recvbuf, recvcount, type, source, tag);
        }

        template <typename T>
        void Sendrecv(const std::vector<T>& sendbuf, int dest,
                            std::vector<T>& recvbuf, int source, int tag,
                      const MPI::Datatype& type) const
        {
            comm.Sendrecv(sendbuf.data(), sendbuf.size(), type, dest, tag,
                          recvbuf.data(), recvbuf.size(), type, source, tag);
        }

        template <typename T>
        void Sendrecv(T* buf, int count, int dest, int source, int tag, MPI::Status& status) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Sendrecv_replace(buf, count, type, dest, tag, source, tag, status);
        }

        template <typename T>
        void Sendrecv(std::vector<T>& buf, int dest, int source, int tag, MPI::Status& status) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Sendrecv_replace(buf.data(), buf.size(), type, dest, tag, source, tag, status);
        }

        template <typename T>
        void Sendrecv(T* buf, int count, int dest, int source, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Sendrecv_replace(buf, count, type, dest, tag, source, tag);
        }

        template <typename T>
        void Sendrecv(std::vector<T>& buf, int dest, int source, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Sendrecv_replace(buf.data(), buf.size(), type, dest, tag, source, tag);
        }

        template <typename T>
        void Sendrecv(T* buf, int count, int dest, int source, int tag, MPI::Status& status,
                      const MPI::Datatype& type) const
        {
            comm.Sendrecv_replace(buf, count, type, dest, tag, source, tag, status);
        }

        template <typename T>
        void Sendrecv(std::vector<T>& buf, int dest, int source, int tag, MPI::Status& status,
                      const MPI::Datatype& type) const
        {
            comm.Sendrecv_replace(buf.data(), buf.size(), type, dest, tag, source, tag, status);
        }

        template <typename T>
        void Sendrecv(T* buf, int count, int dest, int source, int tag, const MPI::Datatype& type) const
        {
            comm.Sendrecv_replace(buf, count, type, dest, tag, source, tag);
        }

        template <typename T>
        void Sendrecv(std::vector<T>& buf, int dest, int source, int tag, const MPI::Datatype& type) const
        {
            comm.Sendrecv_replace(buf.data(), buf.size(), type, dest, tag, source, tag);
        }

        template <typename T>
        void Ssend(const T* buf, int count, int dest, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Ssend(buf, count, type, dest, tag);
        }

        template <typename T>
        void Ssend(const std::vector<T>& buf, int dest, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            comm.Ssend(buf.data(), buf.size(), type, dest, tag);
        }

        template <typename T>
        void Ssend(const T* buf, int count, int dest, int tag, const MPI::Datatype& type) const
        {
            comm.Ssend(buf, count, type, dest, tag);
        }

        template <typename T>
        void Ssend(const std::vector<T>& buf, int dest, int tag, const MPI::Datatype& type) const
        {
            comm.Ssend(buf.data(), buf.size(), type, dest, tag);
        }

        template <typename T>
        MPI::Prequest Ssend_init(const T* buf, int count, int dest, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            return comm.Ssend_init(buf, count, type, dest, tag);
        }

        template <typename T>
        MPI::Prequest Ssend_init(const std::vector<T>& buf, int dest, int tag) const
        {
            const MPI::Datatype& type = MPI_TYPE_<T>::value();
            return comm.Ssend_init(buf.data(), buf.size(), type, dest, tag);
        }

        template <typename T>
        MPI::Prequest Ssend_init(const T* buf, int count, int dest, int tag, const MPI::Datatype& type) const
        {
            return comm.Ssend_init(buf, count, type, dest, tag);
        }

        template <typename T>
        MPI::Prequest Ssend_init(const std::vector<T>& buf, int dest, int tag, const MPI::Datatype& type) const
        {
            return comm.Ssend_init(buf.data(), buf.size(), type, dest, tag);
        }
};

/*
template <>
inline tCTF_World<float>& Arena::ctf<float>()
{
    if (!ctfs) ctfs = new tCTF_World<float>(comm);
    return *ctfs;
}
*/

template <>
inline tCTF_World<double>& Arena::ctf<double>()
{
    if (!ctfd) ctfd.reset(new tCTF_World<double>(comm));
    return *ctfd;
}

/*
template <>
inline tCTF_World<std::complex<float> >& Arena::ctf<std::complex<float> >()
{
    if (!ctfc) ctfc = new tCTF_World<std::complex<float> >(comm);
    return *ctfc;
}

template <>
inline tCTF_World<std::complex<double> >& Arena::ctf<std::complex<double> >()
{
    if (!ctfz) ctfz = new tCTF_World<std::complex<double> >(comm);
    return *ctfz;
}
*/

class Distributed
{
    public:
        Arena arena;
        const int rank;
        const int nproc;

        Distributed(const Arena& arena) : arena(arena), rank(arena.rank), nproc(arena.nproc) {}

        static std::ostream& debug()
        {
            static std::fstream fs(std::strprintf("debug.%d", MPI::COMM_WORLD.Get_rank()).c_str(),
                                   std::fstream::out|std::fstream::trunc);
            return fs;
        }
};

}

#endif

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

#include <complex>

#include "mpi.h"
#ifdef USE_ELEMENTAL
#include "elemental.hpp"
#endif

#include "ctf.hpp"

namespace aquarius
{

template <typename T>
struct MPI_TYPE_ {};

template <>
struct MPI_TYPE_<float>
{
    static MPI::Datatype value() { return MPI::FLOAT; }
};

template <>
struct MPI_TYPE_<double>
{
    static MPI::Datatype value() { return MPI::DOUBLE; }
};

template <>
struct MPI_TYPE_< std::complex<float> >
{
    static MPI::Datatype value() { return MPI::COMPLEX; }
};

template <>
struct MPI_TYPE_< std::complex<double> >
{
    static MPI::Datatype value() { return MPI::DOUBLE_COMPLEX; }
};

template <>
struct MPI_TYPE_<short>
{
    static MPI::Datatype value() { return MPI::SHORT; }
};

template <>
struct MPI_TYPE_<unsigned short>
{
    static MPI::Datatype value() { return MPI::UNSIGNED_SHORT; }
};

template <>
struct MPI_TYPE_<int>
{
    static MPI::Datatype value() { return MPI::INT; }
};

template <>
struct MPI_TYPE_<unsigned int>
{
    static MPI::Datatype value() { return MPI::UNSIGNED; }
};

template <>
struct MPI_TYPE_<long>
{
    static MPI::Datatype value() { return MPI::LONG; }
};

template <>
struct MPI_TYPE_<unsigned long>
{
    static MPI::Datatype value() { return MPI::UNSIGNED_LONG; }
};

template <>
struct MPI_TYPE_<long long>
{
    static MPI::Datatype value() { return MPI::LONG_LONG; }
};

template <>
struct MPI_TYPE_<unsigned long long>
{
    static MPI::Datatype value() { return MPI::UNSIGNED_LONG_LONG; }
};

template <typename T>
class Arena
{
    public:
        tCTF_World<T> ctf;
        MPI::Intracomm comm;
        const int rank;
        const int nproc;

        Arena(MPI::Intracomm& comm = MPI::COMM_WORLD)
        : ctf(comm), comm(comm), rank(comm.Get_rank()), nproc(comm.Get_size()) {}

        Arena(int argc, char** argv)
        : ctf(MPI::COMM_WORLD, argc, argv), comm(MPI::COMM_WORLD), rank(MPI::COMM_WORLD.Get_rank()), nproc(MPI::COMM_WORLD.Get_size()) {}
};

template <typename T>
class Distributed
{
    public:
        Arena<T>& arena;
        MPI::Intracomm comm;
        const int rank;
        const int nproc;
        MPI::Datatype type;

        Distributed(Arena<T>& arena)
        : arena(arena), comm(arena.comm), rank(arena.rank), nproc(arena.nproc), type(MPI_TYPE_<T>::value()) {}
};

}

#endif

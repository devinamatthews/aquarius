#ifndef _AQUARIUS_UTIL_DISTRIBUTED_HPP_
#define _AQUARIUS_UTIL_DISTRIBUTED_HPP_

#define MPIWRAP_LONG_LONG
#include "mpiwrap.hpp"

#ifdef HAVE_ELEMENTAL
#include "El.hpp"
#endif

#include "ctf.hpp"
#include "stl_ext.hpp"

namespace aquarius
{

using namespace MPIWrap;

class Arena
{
    protected:
        //global_ptr<tCTF_World<float>> ctfs;
        global_ptr<tCTF_World<double>> ctfd;
        //global_ptr<tCTF_World<complex<float>>> ctfc;
        //global_ptr<tCTF_World<complex<double>>> ctfz;
        shared_ptr<Intracomm> comm_;

    public:
        const int rank;
        const int size;

        Arena(Intracomm&& comm)
        : comm_(new Intracomm(move(comm))), rank(comm.rank), size(comm.size) {}

        Arena() : Arena(Intracomm::world()) {}

        Intracomm& comm() { return *comm_; }

        const Intracomm& comm() const { return *comm_; }

        template <typename T>
        tCTF_World<T>& ctf();

        template <typename T>
        const tCTF_World<T>& ctf() const
        {
            return const_cast<Arena&>(*this).ctf<T>();
        }
};

Arena& arena();

class Distributed
{
    protected:
        Arena arena;

    public:
        Distributed(const Arena& arena) : arena(arena) {}

        const Arena& getArena() const { return arena; }
};

}

#endif

#ifndef _AQUARIUS_UTIL_DISTRIBUTED_HPP_
#define _AQUARIUS_UTIL_DISTRIBUTED_HPP_

#define MPIWRAP_LONG_LONG
#include "external/mpiwrap/mpiwrap.hpp"

#ifdef HAVE_ELEMENTAL
#include "El.hpp"
#endif

#include "ctf.hpp"
#include "util/stl_ext.hpp"

namespace aquarius
{

using namespace MPIWrap;

class Arena
{
    protected:
        global_ptr<CTF::World> ctf_;
        shared_ptr<Intracomm> comm_;

    public:
        const int rank;
        const int size;

        Arena(Intracomm&& comm)
        : comm_(new Intracomm(move(comm))), rank(comm.rank), size(comm.size) {}

        Arena() : Arena(Intracomm::world()) {}

        Intracomm& comm() { return *comm_; }

        const Intracomm& comm() const { return *comm_; }

        CTF::World& ctf()
        {
            if (!ctf_) ctf_.reset(new CTF::World(*comm_));
            return *ctf_;
        }

        const CTF::World& ctf() const
        {
            return const_cast<const CTF::World&>(const_cast<Arena&>(*this).ctf());
        }
};

class Distributed
{
    public:
        Arena arena;

        Distributed(const Arena& arena) : arena(arena) {}

        static ostream& debug()
        {
            static Arena world;
            static fstream fs(strprintf("debug.%d", world.rank).c_str(),
                                   fstream::out|fstream::trunc);
            return fs;
        }
};

}

#endif

#ifndef _AQUARIUS_UTIL_DISTRIBUTED_HPP_
#define _AQUARIUS_UTIL_DISTRIBUTED_HPP_

#include "external/mpiwrap/mpiwrap.hpp"

#ifdef ELEMENTAL
#include "El.hpp"
#endif

#include "ctf.hpp"
#include "util/stl_ext.hpp"

namespace aquarius
{

using namespace MPIWrap;

namespace tensor
{

template <typename T> class CTFTensor;
template <typename T> class SymmetryBlockedTensor;
template <typename T> class SpinorbitalTensor;

}

class Arena : public Intracomm
{
    protected:
        //global_ptr<tCTF_World<float> > ctfs;
        global_ptr<tCTF_World<double> > ctfd;
        //global_ptr<tCTF_World<complex<float> > > ctfc;
        //global_ptr<tCTF_World<complex<double> > > ctfz;

    public:
        Arena(const MPI_Comm& comm = MPI_COMM_WORLD) : Intracomm(comm) {}

        template <typename T>
        tCTF_World<T>& ctf();

        template <typename T>
        const tCTF_World<T>& ctf() const
        {
            return const_cast<const tCTF_World<T>&>(const_cast<Arena&>(*this).ctf<T>());
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
inline tCTF_World<complex<float> >& Arena::ctf<complex<float> >()
{
    if (!ctfc) ctfc = new tCTF_World<complex<float> >(comm);
    return *ctfc;
}

template <>
inline tCTF_World<complex<double> >& Arena::ctf<complex<double> >()
{
    if (!ctfz) ctfz = new tCTF_World<complex<double> >(comm);
    return *ctfz;
}
*/

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

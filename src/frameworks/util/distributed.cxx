#include "distributed.hpp"

namespace aquarius
{

Arena& arena()
{
    static Arena world;
    return world;
}

template <>
tCTF_World<double>& Arena::ctf<double>()
{
    if (!ctfd) ctfd.set(new tCTF_World<double>(*comm_));
    return *ctfd;
}

}

#endif

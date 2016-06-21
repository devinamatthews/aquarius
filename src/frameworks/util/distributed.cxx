#include "distributed.hpp"

namespace aquarius
{

Arena& arena()
{
    static Arena world;
    return world;
}

CTF::World& Arena::ctf()
{
    if (!ctf_) ctf_.set(new CTF::World(*comm_));
    return *ctf_;
}

}

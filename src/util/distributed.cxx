#include "distributed.hpp"

namespace aquarius
{

const Arena& world()
{
    static Arena world_;
    return world_;
}

}

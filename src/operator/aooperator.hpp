#ifndef _AQUARIUS_OPERATOR_AOOPERATOR_HPP_
#define _AQUARIUS_OPERATOR_AOOPERATOR_HPP_

#include "util/global.hpp"

namespace aquarius
{
namespace op
{

class AOOperator : public Distributed
{
    public:
        const int norb;

        AOOperator(const Arena& arena, int norb) : Distributed(arena), norb(norb) {}
};

}
}

#endif

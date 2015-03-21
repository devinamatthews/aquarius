#ifndef _AQUARIUS_OPERATOR_MOOPERATOR_HPP_
#define _AQUARIUS_OPERATOR_MOOPERATOR_HPP_

#include "util/global.hpp"

#include "task/task.hpp"

#include "space.hpp"

namespace aquarius
{
namespace op
{

class MOOperator : public Distributed
{
    public:
        const Space occ;
        const Space vrt;

        MOOperator(const Arena& arena, const Space& occ, const Space& vrt)
        : Distributed(arena), occ(occ), vrt(vrt) {}
};

}
}

#endif

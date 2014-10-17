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

#ifndef _AQUARIUS_OPERATOR_SPACE_HPP_
#define _AQUARIUS_OPERATOR_SPACE_HPP_

#include "tensor/symblocked_tensor.hpp"
#include "util/stl_ext.hpp"
#include "util/distributed.hpp"
#include "task/task.hpp"

namespace aquarius
{
namespace op
{

struct Space
{
    const symmetry::PointGroup& group;
    const std::vector<int> nalpha;
    const std::vector<int> nbeta;

    Space(const symmetry::PointGroup& group)
    : group(group), nalpha(group.getNumIrreps(), 0), nbeta(group.getNumIrreps(), 0) {}

    Space(const symmetry::PointGroup& group, const std::vector<int>& nalpha, const std::vector<int>& nbeta)
    : group(group), nalpha(nalpha), nbeta(nbeta)
    {
        assert(nalpha.size() == group.getNumIrreps());
        assert(nbeta.size() == group.getNumIrreps());
    }

    Space& operator=(const Space& other)
    {
        assert(group == other.group);
        const_cast<std::vector<int>&>(nalpha).assign(other.nalpha.begin(), other.nalpha.end());
        const_cast<std::vector<int>&>(nbeta).assign(other.nbeta.begin(), other.nbeta.end());
        return *this;
    }

    bool operator==(const Space& other) const
    {
        return nalpha == other.nalpha && nbeta == other.nbeta;
    }
};

template <class T>
struct MOSpace : public Space, public Distributed
{
    const std::vector<int> nao;
    const tensor::SymmetryBlockedTensor<T> Calpha;
    const tensor::SymmetryBlockedTensor<T> Cbeta;

    MOSpace(const tensor::SymmetryBlockedTensor<T>& Calpha, const tensor::SymmetryBlockedTensor<T>& Cbeta)
    : Space(Calpha.getGroup(), Calpha.getLengths()[1], Cbeta.getLengths()[1]),
      Distributed(Calpha.arena),
      nao(Calpha.getLengths()[0]),
      Calpha(Calpha),
      Cbeta(Cbeta) {}

    MOSpace(tensor::SymmetryBlockedTensor<T>* Calpha, tensor::SymmetryBlockedTensor<T>* Cbeta)
    : Space(Calpha->getGroup(), Calpha->getLengths()[1], Cbeta->getLengths()[1]),
      Distributed(Calpha->arena),
      nao(Calpha->getLengths()[0]),
      Calpha(Calpha),
      Cbeta(Cbeta) {}
};

}
}

#endif

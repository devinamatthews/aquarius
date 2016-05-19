#ifndef _AQUARIUS_OPERATOR_SPACE_HPP_
#define _AQUARIUS_OPERATOR_SPACE_HPP_

#include "util/global.hpp"

#include "tensor/symblocked_tensor.hpp"
#include "task/task.hpp"

namespace aquarius
{
namespace op
{

struct Space
{
    const symmetry::PointGroup& group;
    const vector<int> nalpha;
    const vector<int> nbeta;

    Space(const symmetry::PointGroup& group)
    : group(group), nalpha(group.getNumIrreps(), 0), nbeta(group.getNumIrreps(), 0) {}

    Space(const symmetry::PointGroup& group, const vector<int>& nalpha, const vector<int>& nbeta)
    : group(group), nalpha(nalpha), nbeta(nbeta)
    {
        assert(nalpha.size() == group.getNumIrreps());
        assert(nbeta.size() == group.getNumIrreps());
    }

    Space& operator=(const Space& other)
    {
        assert(group == other.group);
        const_cast<vector<int>&>(nalpha).assign(other.nalpha.begin(), other.nalpha.end());
        const_cast<vector<int>&>(nbeta).assign(other.nbeta.begin(), other.nbeta.end());
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
    const vector<int> nao;
    const tensor::SymmetryBlockedTensor<T> Calpha;
    const tensor::SymmetryBlockedTensor<T> Cbeta;

    MOSpace(const tensor::SymmetryBlockedTensor<T>& Calpha, const tensor::SymmetryBlockedTensor<T>& Cbeta)
    : Space(Calpha.getGroup(), Calpha.getLengths()[1], Cbeta.getLengths()[1]),
      Distributed(Calpha.arena),
      nao(Calpha.getLengths()[0]),
      Calpha(Calpha),
      Cbeta(Cbeta) {}

    MOSpace(tensor::SymmetryBlockedTensor<T>&& Calpha, tensor::SymmetryBlockedTensor<T>&& Cbeta)
    : Space(Calpha.getGroup(), Calpha.getLengths()[1], Cbeta.getLengths()[1]),
      Distributed(Calpha.arena),
      nao(Calpha.getLengths()[0]),
      Calpha(move(Calpha)),
      Cbeta(move(Cbeta)) {}
};

}
}

#endif

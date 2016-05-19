#ifndef _AQUARIUS_OPERATOR_MULTIPOLE_HPP_
#define _AQUARIUS_OPERATOR_MULTIPOLE_HPP_

#include "../scf/uhf.hpp"
#include "util/global.hpp"

#include "1eoperator.hpp"

namespace aquarius
{
namespace op
{

template <typename T>
class Multipole : public MOOperator<T>, public tensor::CompositeTensor<Multipole<T>,OneElectronOperator<T>,T>
{
    protected:
        int Lmin, Lmax;

    public:
        Multipole(const string& name, const scf::UHF<T>& uhf, int Lmin, int Lmax=-1);

        const OneElectronOperator<T>& operator()(int L, int xyz) const;

        const OneElectronOperator<T>& operator()(int x, int y, int z) const;
};

}
}

#endif

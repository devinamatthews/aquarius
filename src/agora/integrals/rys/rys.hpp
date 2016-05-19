#include "util/global.hpp"

namespace aquarius
{
namespace integrals
{

class Rys
{
    public:
        /**
         * generate the roots and weights of the Rys quadrature
         *
         * see Golub, G. H.; Welsch, J. H. Math. Comput. 23, 221-230 (1969)
         *     K. Ishida, J. Chem. Phys. 95, 5198-205 (1991)
         */
        void operator()(double T, int n, double* rt, double* wt);
};

}
}

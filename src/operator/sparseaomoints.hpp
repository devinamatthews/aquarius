#ifndef _AQUARIUS_OPERATOR_SPARSEAOMOINTS_HPP_
#define _AQUARIUS_OPERATOR_SPARSEAOMOINTS_HPP_

#include "util/global.hpp"

#include "scf/aouhf.hpp"
#include "integrals/2eints.hpp"

#include "moints.hpp"

namespace aquarius
{
namespace op
{

template <typename T>
class SparseAOMOIntegrals : public MOIntegrals<T>
{
    public:
        SparseAOMOIntegrals(const string& name, input::Config& config);

    private:
        enum Side {NONE, PQ, RS};
        enum Index {A, B};

        struct pqrs_integrals : Distributed
        {
            const symmetry::PointGroup& group;
            vector<int> np, nq, nr, ns;
            deque<T> ints;
            deque<idx4_t> idxs;

            pqrs_integrals(const Arena& arena, const symmetry::PointGroup& group) : Distributed(arena), group(group) {}

            /*
             * Read integrals in and break (pq|rs)=(rs|pq) symmetry
             */
            pqrs_integrals(const vector<int>& norb, const integrals::ERI& aoints);

            void free();

            void sortInts(bool rles, size_t& nrs, vector<size_t>& rscount);

            /*
             * Redistribute integrals such that each node has all pq for each rs pair
             */
            void collect(bool rles);

            /*
             * Transform (ab|rs) -> (cb|rs) (index = A) or (ab|rs) -> (ac|rs) (index = B)
             */
            pqrs_integrals transform(Index index, const vector<int>& nc, const vector<vector<T> >& C, bool pleq);

            void transcribe(tensor::SymmetryBlockedTensor<T>& tensor, bool assymij, bool assymkl, Side swap);

            void transpose();
        };

    protected:
        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif

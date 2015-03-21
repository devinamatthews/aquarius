#ifndef _AQUARIUS_OPERATOR_AOMOINTS_HPP_
#define _AQUARIUS_OPERATOR_AOMOINTS_HPP_

#include "util/global.hpp"

#include "scf/aouhf.hpp"
#include "integrals/2eints.hpp"

#include "moints.hpp"

namespace aquarius
{
namespace op
{

template <typename T>
class AOMOIntegrals : public MOIntegrals<T>
{
    public:
        AOMOIntegrals(const string& name, input::Config& config);

    private:
        enum Side {NONE, PQ, RS};
        enum Index {A, B};

        struct abrs_integrals;

        struct pqrs_integrals : Distributed
        {
            const symmetry::PointGroup& group;
            vector<int> np, nq, nr, ns;
            vector<T> ints;
            vector<idx4_t> idxs;

            /*
             * Read integrals in and break (pq|rs)=(rs|pq) symmetry
             */
            pqrs_integrals(const vector<int>& norb, const integrals::ERI& aoints);

            pqrs_integrals(abrs_integrals& abrs);

            void free();

            void sortInts(bool rles, size_t& nrs, vector<size_t>& rscount);

            /*
             * Redistribute integrals such that each node has all pq for each rs pair
             */
            void collect(bool rles);

            /*
             * Transform (ab|rs) -> (cb|rs) (index = A) or (ab|rs) -> (ac|rs) (index = B)
             *
             * C is ldc*nc if trans = 'N' and ldc*[na|nb] if trans = 'T'
             */
            pqrs_integrals transform(Index index, const vector<int>& nc, const vector<vector<T> >& C, bool pleq);

            void transcribe(tensor::SymmetryBlockedTensor<T>& tensor, bool assymij, bool assymkl, Side swap);

            void transpose();
        };

        struct abrs_integrals : Distributed
        {
            const symmetry::PointGroup& group;
            vector<int> na, nb, nr, ns;
            vector<T> ints;
            vector<idx2_t> rs;

            abrs_integrals(const Arena& arena, const symmetry::PointGroup& group)
            : Distributed(arena), group(group) {}

            /*
             * Expand (p_i q_j|r_k s_l) into (pq|r_k s_l), where pleq = true indicates
             * that (pq|rs) = (qp|rs) and only p_i <= q_j is stored in the input
             */
            abrs_integrals(pqrs_integrals& pqrs, const bool pleq);

            /*
             * Transform (ab|rs) -> (cb|rs) (index = A) or (ab|rs) -> (ac|rs) (index = B)
             *
             * C is ldc*nc if trans = 'N' and ldc*[na|nb] if trans = 'T'
             */
            abrs_integrals transform(Index index, const vector<int>& nc, const vector<vector<T> >& C);

            void transcribe(tensor::SymmetryBlockedTensor<T>& tensor, bool assymij, bool assymkl, Side swap);

            void free();

            size_t getNumAB(idx2_t rs);

            size_t getNumAB(idx2_t rs, vector<size_t>& offab);
        };

    protected:
        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif

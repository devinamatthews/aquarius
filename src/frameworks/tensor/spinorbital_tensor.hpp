#ifndef _AQUARIUS_TENSOR_SPINORBITAL_TENSOR_HPP_
#define _AQUARIUS_TENSOR_SPINORBITAL_TENSOR_HPP_

#include "util/global.hpp"

#include "autocc/autocc.hpp"
#include "operator/space.hpp"
#include "task/task.hpp"

#include "symblocked_tensor.hpp"
#include "composite_tensor.hpp"

namespace aquarius
{
namespace tensor
{

template<class T>
class SpinorbitalTensor : public IndexableCompositeTensor<SpinorbitalTensor<T>,SymmetryBlockedTensor<T>,T>,
                          public Distributed
{
    INHERIT_FROM_INDEXABLE_COMPOSITE_TENSOR(SpinorbitalTensor<T>,SymmetryBlockedTensor<T>,T)

    public:
        SpinorbitalTensor(const string& name, const SpinorbitalTensor<T>& t, const T val);

        SpinorbitalTensor(const SpinorbitalTensor<T>& other);

        SpinorbitalTensor(const string& name, const SpinorbitalTensor<T>& other);

        SpinorbitalTensor(const string& name, const Arena& arena,
                          const symmetry::PointGroup& group,
                          const vector<op::Space>& spaces,
                          const vector<int>& nout,
                          const vector<int>& nin, int spin=0);

        SpinorbitalTensor(const string& name, const Arena& arena,
                          const symmetry::PointGroup& group,
                          const symmetry::Representation& rep,
                          const vector<op::Space>& spaces,
                          const vector<int>& nout,
                          const vector<int>& nin, int spin=0);

        ~SpinorbitalTensor();

        const vector<int>& getNumOut() const { return nout; }

        const vector<int>& getNumIn() const { return nin; }

        int getSpin() const { return spin; }

        const symmetry::PointGroup& getGroup() const { return group; }

        SymmetryBlockedTensor<T>& operator()(const vector<int>& alpha_out,
                                             const vector<int>& alpha_in);

        const SymmetryBlockedTensor<T>& operator()(const vector<int>& alpha_out,
                                                   const vector<int>& alpha_in) const;

        void mult(const T alpha, bool conja, const SpinorbitalTensor<T>& A_, const string& idx_A,
                                 bool conjb, const SpinorbitalTensor<T>& B_, const string& idx_B,
                  const T beta_,                                             const string& idx_C);

        void sum(const T alpha, bool conja, const SpinorbitalTensor<T>& A_, const string& idx_A,
                 const T beta_,                                             const string& idx_B);

        void scale(const T alpha, const string& idx_A);

        void weight(const vector<const vector<vector<T>>*>& da,
                    const vector<const vector<vector<T>>*>& db,
                    double shift = 0);

        T dot(bool conja, const SpinorbitalTensor<T>& A, const string& idx_A,
              bool conjb,                                const string& idx_B) const;

        real_type_t<T> norm(int p) const;

    protected:
        struct SpinCase
        {
            SymmetryBlockedTensor<T> *tensor;
            vector<int> alpha_out, alpha_in;

            void construct(SpinorbitalTensor<T>& t,
                           const symmetry::Representation& rep,
                           const vector<int>& alpha_out,
                           const vector<int>& alpha_in);
        };

        const symmetry::PointGroup& group;
        vector<op::Space> spaces;
        vector<int> nout, nin;
        int spin;
        vector<SpinCase> cases;
        static map<const tCTF_World<T>*,map<const symmetry::PointGroup*,pair<int,SpinorbitalTensor<T>*>>> scalars;

        void register_scalar();

        void unregister_scalar();

        SpinorbitalTensor<T>& scalar() const;
};

}
}

#endif

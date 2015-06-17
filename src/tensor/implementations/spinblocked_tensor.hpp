#ifndef _AQUARIUS_TENSOR_SPINBLOCKED_TENSOR_HPP_
#define _AQUARIUS_TENSOR_SPINBLOCKED_TENSOR_HPP_

#include "util/global.hpp"

#include "autocc/autocc.hpp"

#include "tensor/tensor.hpp"

namespace aquarius
{
namespace tensor
{

template <capability_type C> class SpinBlockedTensor;

namespace detail
{
    template <typename SubFactory, capability_type C, typename dummy=void>
    struct ConstructSpinBlock;

    template <typename SubFactory, capability_type C>
    struct ConstructSpinBlock<SubFactory, C, enable_if_t<!ARE_DISTINCT(PGSYMMETRIC,C)>>
    {
        ConstructSpinBlock(SpinBlockedTensor<C>& parent,
                  const INITIALIZER_TYPE(C&~(SPINORBITAL_|BOUNDED|IPSYMMETRIC|PGSYMMETRIC))& subinit,
                  SubFactory f, const vector<int>& alpha_creation,
                  const vector<int>& alpha_annihilation)
        {
            matrix<int> sublen(parent.ndim, parent.group.getNumIrreps());
            vector<int> subsym(parent.ndim, AS);
            int nclass = parent.nalpha.size();
            int nirrep = parent.group.getNumIrreps();
            vector<int> b(parent.ndim);

            int dim = 0;
            for (int i = 0;i < nclass;i++)
            {
                for (int a = 0;a < alpha_creation[i];a++, dim++)
                {
                    b[dim] = 0;
                    for (int j = 0;j < nirrep;j++)
                    {
                        sublen[dim][j] = parent.nalpha_per_irrep[i][j];
                    }
                }
                if (alpha_creation[i] > 0) subsym[dim-1] = NS;

                for (int a = 0;a < parent.ncreation[i]-alpha_creation[i];a++, dim++)
                {
                    b[dim] = 1;
                    for (int j = 0;j < nirrep;j++)
                    {
                        sublen[dim][j] = parent.nbeta_per_irrep[i][j];
                    }
                }
                if (parent.ncreation[i]-alpha_creation[i] > 0) subsym[dim-1] = NS;
            }

            for (int i = 0;i < nclass;i++)
            {
                for (int a = 0;a < alpha_annihilation[i];a++, dim++)
                {
                    b[dim] = 0;
                    for (int j = 0;j < nirrep;j++)
                    {
                        sublen[dim][j] = parent.nalpha_per_irrep[i][j];
                    }
                }
                if (alpha_annihilation[i] > 0) subsym[dim-1] = NS;

                for (int a = 0;a < parent.nannihilation[i]-alpha_annihilation[i];a++, dim++)
                {
                    b[dim] = 1;
                    for (int j = 0;j < nirrep;j++)
                    {
                        sublen[dim][j] = parent.nbeta_per_irrep[i][j];
                    }
                }
                if (parent.nannihilation[i]-alpha_annihilation[i] > 0) subsym[dim-1] = NS;
            }

            Tensor<C&~SPINORBITAL_> *t = new Tensor<C&~SPINORBITAL_>(f(subinit <<
                TensorInitializer<PGSYMMETRIC|BOUNDED|IPSYMMETRIC>(parent.group, sublen, subsym, parent.rep)));
            parent.put(b, t);

            parent.cases.push_back(SpinCase(*t, alpha_creation, alpha_annihilation));
        }
    };

    template <typename SubFactory, capability_type C>
    struct ConstructSpinBlock<SubFactory, C, enable_if_t<ARE_DISTINCT(PGSYMMETRIC,C)>>
    {
        ConstructSpinBlock(SpinBlockedTensor<C>& parent,
                  const INITIALIZER_TYPE(C&~(BOUNDED|IPSYMMETRIC))& subinit,
                  SubFactory f, const vector<int>& alpha_creation,
                  const vector<int>& alpha_annihilation)
        {
            vector<int> sublen(parent.ndim);
            vector<int> subsym(parent.ndim, AS);
            int nclass = parent.nalpha.size();
            vector<int> b(parent.ndim);


            int dim = 0;
            for (int i = 0;i < nclass;i++)
            {
                for (int a = 0;a < alpha_creation[i];a++, dim++)
                {
                    b[dim] = 0;
                    sublen[dim] = parent.nalpha[i];
                }
                if (alpha_creation[i] > 0) subsym[dim-1] = NS;

                for (int a = 0;a < parent.ncreation[i]-alpha_creation[i];a++, dim++)
                {
                    b[dim] = 1;
                    sublen[dim] = parent.nbeta[i];
                }
                if (parent.ncreation[i]-alpha_creation[i] > 0) subsym[dim-1] = NS;
            }

            for (int i = 0;i < nclass;i++)
            {
                for (int a = 0;a < alpha_annihilation[i];a++, dim++)
                {
                    b[dim] = 0;
                    sublen[dim] = parent.nalpha[i];
                }
                if (alpha_annihilation[i] > 0) subsym[dim-1] = NS;

                for (int a = 0;a < parent.nannihilation[i]-alpha_annihilation[i];a++, dim++)
                {
                    b[dim] = 1;
                    sublen[dim] = parent.nbeta[i];
                }
                if (parent.nannihilation[i]-alpha_annihilation[i] > 0) subsym[dim-1] = NS;
            }

            Tensor<C> *t = new Tensor<C>(f(subinit <<
                TensorInitializer<BOUNDED|IPSYMMETRIC>(sublen, subsym)));
            parent.put(b, t);

            parent.cases.push_back(SpinCase(*t, alpha_creation, alpha_annihilation));
        }
    };

    inline int conv_idx(const vector<int>& cidx_A, string& iidx_A)
    {
        iidx_A.resize(cidx_A.size());

        int n = 0;
        for (int i = 0;i < cidx_A.size();i++)
        {
            int j;
            for (j = 0;j < i;j++)
            {
                if (cidx_A[i] == cidx_A[j])
                {
                    iidx_A[i] = iidx_A[j];
                    break;
                }
            }
            if (j == i)
            {
                iidx_A[i] = (char)('A'+n++);
            }
        }

        return n;
    }

    inline int conv_idx(const vector<int>& cidx_A, string& iidx_A,
                        const vector<int>& cidx_B, string& iidx_B)
    {
        iidx_B.resize(cidx_B.size());

        int n = conv_idx(cidx_A, iidx_A);

        for (int i = 0;i < cidx_B.size();i++)
        {
            int j;
            for (j = 0;j < cidx_A.size();j++)
            {
                if (cidx_B[i] == cidx_A[j])
                {
                    iidx_B[i] = iidx_A[j];
                    break;
                }
            }
            if (j == cidx_A.size())
            {
                for (j = 0;j < i;j++)
                {
                    if (cidx_B[i] == cidx_B[j])
                    {
                        iidx_B[i] = iidx_B[j];
                        break;
                    }
                }
                if (j == i)
                {
                    iidx_B[i] = (char)('A'+n++);
                }
            }
        }

        return n;
    }

    inline int conv_idx(const vector<int>& cidx_A, string& iidx_A,
                        const vector<int>& cidx_B, string& iidx_B,
                        const vector<int>& cidx_C, string& iidx_C)
    {
        iidx_C.resize(cidx_C.size());

        int n = conv_idx(cidx_A, iidx_A,
                         cidx_B, iidx_B);

        for (int i = 0;i < cidx_C.size();i++)
        {
            int j;
            for (j = 0;j < cidx_B.size();j++)
            {
                if (cidx_C[i] == cidx_B[j])
                {
                    iidx_C[i] = iidx_B[j];
                    break;
                }
            }
            if (j == cidx_B.size())
            {
                for (j = 0;j < cidx_A.size();j++)
                {
                    if (cidx_C[i] == cidx_A[j])
                    {
                        iidx_C[i] = iidx_A[j];
                        break;
                    }
                }
                if (j == cidx_A.size())
                {
                    for (j = 0;j < i;j++)
                    {
                        if (cidx_C[i] == cidx_C[j])
                        {
                            iidx_C[i] = iidx_C[j];
                            break;
                        }
                    }
                    if (j == i)
                    {
                        iidx_C[i] = (char)('A'+n++);
                    }
                }
            }
        }

        return n;
    }

    template <capability_type C>
    vector<vector<int>>& getBlockLengths(const INITIALIZER_TYPE(C)& init)
    {
        auto& idxinit = init.template as<INDEXABLE>();
        auto& soinit = init.template as<SPINORBITAL_>();

        vector<vector<int>> len(idxinit.ndim, vector<int>(2));

        int p = 0;
        for (int i = 0;i < soinit.ncreation.size();i++)
        {
            for (int j = 0;j < soinit.ncreation[i];j++, p++)
            {
                len[i][0] = soinit.nalpha[i];
                len[i][1] = soinit.nbeta[i];
            }
        }

        for (int i = 0;i < soinit.nannihilation.size();i++)
        {
            for (int j = 0;j < soinit.nannihilation[i];j++, p++)
            {
                len[i][0] = soinit.nalpha[i];
                len[i][1] = soinit.nbeta[i];
            }
        }
        assert(p == len.size());

        return len;
    }
}

template <capability_type C>
class SpinBlockedTensor : public BlockedTensor<C&~SPINORBITAL_,TensorImplementation<C>>
{
    protected:
        struct SpinCase
        {
            Tensor<C&~SPINORBITAL_>& tensor;
            vector<int> alpha_creation,
                        alpha_annihilation;

            SpinCase(Tensor<C&~SPINORBITAL_>& tensor,
                     const vector<int>& alpha_creation,
                     const vector<int>& alpha_annihilation)
            : tensor(tensor), alpha_creation(alpha_creation),
              alpha_annihilation(alpha_annihilation) {}
        };

        vector<SpinCase> cases;

        template <typename SubFactory>
        void initialize(const INITIALIZER_TYPE(C)& init, SubFactory f)
        {
            INITIALIZER_TYPE(C&~(SPINORBITAL_|PGSYMMETRIC|BOUNDED|IPSYMMETRIC)) subinit(init);

            int nspaces = this->nalpha.size();
            int ncreationtot = aquarius::sum(this->ncreation);
            int nannihilationtot = aquarius::sum(this->nannihilation);
            vector<int> whichout(ncreationtot), whichin(nannihilationtot);
            vector<int> alpha_creation(nspaces), alpha_annihilation(nspaces);

            assert(abs(this->spin) <= ncreationtot+nannihilationtot);
            assert(abs(this->spin) >= abs(ncreationtot-nannihilationtot));
            assert(abs(this->spin)%2 == abs(ncreationtot-nannihilationtot)%2);

            for (int alphaout = 0;alphaout <= ncreationtot;alphaout++)
            {
                int alphain = alphaout + (ncreationtot-nannihilationtot-this->spin)/2;
                if (alphain < 0 || alphain > nannihilationtot) continue;

                fill(whichout.begin(), whichout.end(), 0);

                for (bool doneout = false;!doneout;)
                {
                    fill(alpha_creation.begin(), alpha_creation.end(), 0);

                    for (int i = 0;i < alphaout;i++)
                    {
                        alpha_creation[whichout[i]]++;
                    }

                    fill(whichin.begin(), whichin.end(), 0);

                    for (bool donein = false;!donein;)
                    {
                        fill(alpha_annihilation.begin(), alpha_annihilation.end(), 0);

                        for (int i = 0;i < alphain;i++)
                        {
                            alpha_annihilation[whichin[i]]++;
                        }

                        bool ok = true;
                        for (int i = 0;i < this->ncreation.size();i++)
                        {
                            if (alpha_creation[i] > this->ncreation[i] ||
                                alpha_annihilation[i] > this->nannihilation[i])
                            {
                                ok = false;
                                break;
                            }
                        }

                        if (ok)
                        {
                            detail::ConstructSpinBlock<SubFactory, C&~SPINORBITAL_>(*this, subinit, f, alpha_creation, alpha_annihilation);
                        }

                        for (int i = 0;i < alphain;i++)
                        {
                            whichin[i]++;
                            if (i < alphain-1)
                            {
                                if (whichin[i] <= whichin[i+1]) break;
                                whichin[i] = 0;
                            }
                            else
                            {
                                if (whichin[i] < nspaces) break;
                                if (i == alphain-1) donein = true;
                            }
                        }

                        if (alphain == 0) donein = true;
                    }

                    for (int i = 0;i < alphaout;i++)
                    {
                        whichout[i]++;
                        if (i < alphaout-1)
                        {
                            if (whichout[i] <= whichout[i+1]) break;
                            whichout[i] = 0;
                        }
                        else
                        {
                            if (whichout[i] < nspaces) break;
                            if (i == alphaout-1) doneout = true;
                        }
                    }

                    if (alphaout == 0) doneout = true;
                }
            }
        }

    public:
        template <typename SubFactory>
        class Factory
        {
            protected:
                SubFactory f;

            public:
                Factory(SubFactory f) : f(f) {}

                Tensor<C> operator()(const INITIALIZER_TYPE(C)& init)
                {
                    return new SpinBlockedTensor(init, f);
                }
        };

        SpinBlockedTensor(const INITIALIZER_TYPE(C)& init)
        : BlockedTensor<C&~SPINORBITAL_,TensorImplementation<C>>(init, detail::getBlockLengths(init))
        {
            initialize(init, typename Tensor<C&~SPINORBITAL_>::Factory());
        }

        template <typename SubFactory>
        SpinBlockedTensor(const INITIALIZER_TYPE(C)& init, SubFactory f)
        : BlockedTensor<C&~SPINORBITAL_,TensorImplementation<C>>(init, detail::getBlockLengths(init))
        {
            initialize(init, f);
        }

        bool exists(const vector<int>& alpha_creation,
                    const vector<int>& alpha_annihilation) const
        {
            for (int i = 0;i < this->cases.size();i++)
            {
                if (this->cases[i].alpha_creation      == alpha_creation &&
                    this->cases[i].alpha_annihilation  == alpha_annihilation)
                {
                    return true;
                }
            }

            return false;
        }

        Tensor<C&~SPINORBITAL_>& operator()(const vector<int>& alpha_creation,
                                            const vector<int>& alpha_annihilation)
        {
            return const_cast<Tensor<C&~SPINORBITAL_>&>(const_cast<const SpinBlockedTensor&>(*this)(alpha_creation, alpha_annihilation));
        }

        const Tensor<C&~SPINORBITAL_>& operator()(const vector<int>& alpha_creation,
                                                  const vector<int>& alpha_annihilation) const
        {
            for (int i = 0;i < this->cases.size();i++)
            {
                if (this->cases[i].alpha_creation      == alpha_creation &&
                    this->cases[i].alpha_annihilation  == alpha_annihilation)
                {
                    return this->cases[i].tensor;
                }
            }

            return *static_cast<Tensor<C&~SPINORBITAL_>*>(NULL);
        }

        void mult(const Scalar& alpha, bool conja, const TensorImplementation<>& A_, const string& idx_A,
                                       bool conjb, const TensorImplementation<>& B_, const string& idx_B,
                  const Scalar& beta_,                                               const string& idx_C)
        {
            using namespace aquarius::autocc;

            auto& A = static_cast<const SpinBlockedTensor&>(A_);
            auto& B = static_cast<const SpinBlockedTensor&>(B_);

            vector<Scalar> beta(this->cases.size(), beta_);

            for (int sc = 0;sc < this->cases.size();sc++)
            {
                SpinCase& scC = this->cases[sc];

                int nouttot_C = aquarius::sum(this->ncreation);

                string ext;
                vector<Line> lines_C_out(aquarius::sum(this->ncreation));
                vector<Line> lines_C_in(aquarius::sum(this->nannihilation));
                for (int i = 0, s = 0;s < this->ncreation.size();s++)
                {
                    for (int a = 0;a < scC.alpha_creation[s];a++,i++)
                    {
                        if (!contains(ext, idx_C[i])) ext += idx_C[i];
                        lines_C_out[i] = Line(idx_C[i], s, Line::VIRTUAL, Line::ALPHA);
                    }
                    for (int b = 0;b < this->ncreation[s]-scC.alpha_creation[s];b++,i++)
                    {
                        if (!contains(ext, idx_C[i])) ext += idx_C[i];
                        lines_C_out[i] = Line(idx_C[i], s, Line::VIRTUAL, Line::BETA);
                    }
                }
                for (int i = 0,s = 0;s < this->ncreation.size();s++)
                {
                    for (int a = 0;a < scC.alpha_annihilation[s];a++,i++)
                    {
                        if (!contains(ext, idx_C[i+nouttot_C])) ext += idx_C[i+nouttot_C];
                        lines_C_in[i] = Line(idx_C[i+nouttot_C], s, Line::VIRTUAL, Line::ALPHA);
                    }
                    for (int b = 0;b < this->nannihilation[s]-scC.alpha_annihilation[s];b++,i++)
                    {
                        if (!contains(ext, idx_C[i+nouttot_C])) ext += idx_C[i+nouttot_C];
                        lines_C_in[i] = Line(idx_C[i+nouttot_C], s, Line::VIRTUAL, Line::BETA);
                    }
                }

                int nouttot_A = aquarius::sum(A.ncreation);

                string sum;
                vector<Line> lines_A_out(sum(A.ncreation));
                vector<Line> lines_A_in(sum(A.nannihilation));
                vector<Line> lines_AandC_out, lines_AandC_in;
                for (int i = 0, s = 0;s < A.ncreation.size();s++)
                {
                    for (int a = 0;a < A.ncreation[s];a++,i++)
                    {
                        if (contains(ext, idx_A[i]))
                        {
                            int j;for (j = 0;idx_C[j] != idx_A[i];j++);
                            if (j < nouttot_C)
                            {
                                lines_A_out[i] = lines_C_out[j];
                                lines_AandC_out += lines_C_out[j];
                            }
                            else
                            {
                                lines_A_out[i] = lines_C_in[j-nouttot_C];
                                lines_AandC_out += lines_C_in[j-nouttot_C];
                            }
                        }
                        else
                        {
                            if (!contains(sum, idx_A[i])) sum += idx_A[i];
                            lines_A_out[i] = Line(idx_A[i], s, Line::VIRTUAL, Line::BETA);
                        }
                    }
                }
                for (int i = 0, s = 0;s < A.ncreation.size();s++)
                {
                    for (int a = 0;a < A.nannihilation[s];a++,i++)
                    {
                        if (contains(ext, idx_A[i+nouttot_A]))
                        {
                            int j; for (j = 0;idx_C[j] != idx_A[i+nouttot_A];j++);
                            if (j < nouttot_C)
                            {
                                lines_A_in[i] = lines_C_out[j];
                                lines_AandC_in += lines_C_out[j];
                            }
                            else
                            {
                                lines_A_in[i] = lines_C_in[j-nouttot_C];
                                lines_AandC_in += lines_C_in[j-nouttot_C];
                            }
                        }
                        else
                        {
                            if (!contains(sum, idx_A[i+nouttot_A])) sum += idx_A[i+nouttot_A];
                            lines_A_in[i] = Line(idx_A[i+nouttot_A], s, Line::VIRTUAL, Line::BETA);
                        }
                    }
                }

                int nouttot_B = sum(B.ncreation);

                vector<Line> lines_B_out(sum(B.ncreation));
                vector<Line> lines_B_in(sum(B.nannihilation));
                vector<Line> lines_BandC_out, lines_BandC_in;
                for (int i = 0, s = 0;s < B.ncreation.size();s++)
                {
                    for (int a = 0;a < B.ncreation[s];a++,i++)
                    {
                        if (contains(ext, idx_B[i]))
                        {
                            int j; for (j = 0;idx_C[j] != idx_B[i];j++);
                            if (j < nouttot_C)
                            {
                                lines_B_out[i] = lines_C_out[j];
                                lines_BandC_out += lines_C_out[j];
                            }
                            else
                            {
                                lines_B_out[i] = lines_C_in[j-nouttot_C];
                                lines_BandC_out += lines_C_in[j-nouttot_C];
                            }
                        }
                        else
                        {
                            if (!contains(sum, idx_B[i])) sum += idx_B[i];
                            lines_B_out[i] = Line(idx_B[i], s, Line::VIRTUAL, Line::BETA);
                        }
                    }
                }
                for (int i = 0, s = 0;s < B.ncreation.size();s++)
                {
                    for (int a = 0;a < B.nannihilation[s];a++,i++)
                    {
                        if (contains(ext, idx_B[i+nouttot_B]))
                        {
                            int j; for (j = 0;idx_C[j] != idx_B[i+nouttot_B];j++);
                            if (j < nouttot_C)
                            {
                                lines_B_in[i] = lines_C_out[j];
                                lines_BandC_in += lines_C_out[j];
                            }
                            else
                            {
                                lines_B_in[i] = lines_C_in[j-nouttot_C];
                                lines_BandC_in += lines_C_in[j-nouttot_C];
                            }
                        }
                        else
                        {
                            if (!contains(sum, idx_B[i+nouttot_B])) sum += idx_B[i+nouttot_B];
                            lines_B_in[i] = Line(idx_B[i+nouttot_B], s, Line::VIRTUAL, Line::BETA);
                        }
                    }
                }

                vector<Line> lines_CnotAB_out = lines_C_out;
                uniq(lines_AandC_out);
                uniq(lines_BandC_out);
                uniq(lines_CnotAB_out);
                exclude(lines_CnotAB_out, lines_AandC_out);
                exclude(lines_CnotAB_out, lines_BandC_out);

                vector<Line> lines_CnotAB_in = lines_C_in;
                uniq(lines_AandC_in);
                uniq(lines_BandC_in);
                uniq(lines_CnotAB_in);
                exclude(lines_CnotAB_in, lines_AandC_in);
                exclude(lines_CnotAB_in, lines_BandC_in);

                Diagram d = Diagram(Diagram::SPINORBITAL,
                                    vec(Term(Diagram::SPINORBITAL)*
                                        Fragment("A", lines_A_out, lines_A_in)*
                                        Fragment("B", lines_B_out, lines_B_in)));

                for (int s = 0;s < max(max(A.ncreation.size(),B.ncreation.size()),this->ncreation.size());s++)
                {
                    vector<vector<Line> > assym(3);

                    for (vector<Line>::iterator i = lines_AandC_out.begin();i != lines_AandC_out.end();++i)
                        if (i->getType() == s) assym[0].push_back(*i);
                    for (vector<Line>::iterator i = lines_BandC_out.begin();i != lines_BandC_out.end();++i)
                        if (i->getType() == s) assym[1].push_back(*i);
                    for (vector<Line>::iterator i = lines_CnotAB_out.begin();i != lines_CnotAB_out.end();++i)
                        if (i->getType() == s) assym[2].push_back(*i);

                    if (assym[2].empty()) assym.erase(assym.begin()+2);
                    if (assym[1].empty()) assym.erase(assym.begin()+1);
                    if (assym[0].empty()) assym.erase(assym.begin()+0);
                    if (!assym.empty()) d.antisymmetrize(assym);
                }

                for (int s = 0;s < max(max(A.ncreation.size(),B.ncreation.size()),this->ncreation.size());s++)
                {
                    vector<vector<Line> > assym(3);

                    for (vector<Line>::iterator i = lines_AandC_in.begin();i != lines_AandC_in.end();++i)
                        if (i->getType() == s) assym[0].push_back(*i);
                    for (vector<Line>::iterator i = lines_BandC_in.begin();i != lines_BandC_in.end();++i)
                        if (i->getType() == s) assym[1].push_back(*i);
                    for (vector<Line>::iterator i = lines_CnotAB_in.begin();i != lines_CnotAB_in.end();++i)
                        if (i->getType() == s) assym[2].push_back(*i);

                    if (assym[2].empty()) assym.erase(assym.begin()+2);
                    if (assym[1].empty()) assym.erase(assym.begin()+1);
                    if (assym[0].empty()) assym.erase(assym.begin()+0);
                    if (!assym.empty()) d.antisymmetrize(assym);
                }

                d.convert(Diagram::UHF);

                /*
                 * Remove terms which are antisymmetrizations of same-spin groups
                 */
                for (int s = 0;s < max(max(A.ncreation.size(),B.ncreation.size()),this->ncreation.size());s++)
                {
                    for (int spin = 0;spin < 2;spin++)
                    {
                        vector<Term> terms = d.getTerms();
                        for (vector<Term>::iterator t1 = terms.begin();t1 != terms.end();++t1)
                        {
                            for (vector<Term>::iterator t2 = t1+1;t2 != terms.end();++t2)
                            {
                                if (Term(*t1).fixorder(filter_copy(t1->indices(), and1(isSpin(spin),isType(s)))) ==
                                    Term(*t2).fixorder(filter_copy(t2->indices(), and1(isSpin(spin),isType(s)))))
                                {
                                    d -= *t1;
                                    break;
                                }
                            }
                        }
                    }
                }

                d *= Term(Diagram::UHF)*Fragment("C", lines_C_out, lines_C_in);
                d.fixorder(true);

                for (vector<Term>::const_iterator t = d.getTerms().begin();t != d.getTerms().end();++t)
                {
                    double diagFactor = t->getFactor();

                    vector<Fragment>::const_iterator fA, fB, fC;
                    for (vector<Fragment>::const_iterator f = t->getFragments().begin();f != t->getFragments().end();++f)
                    {
                        if (f->getOp() == "A") fA = f;
                        if (f->getOp() == "B") fB = f;
                        if (f->getOp() == "C") fC = f;
                    }

                    vector<Line> out_A = fA->getIndicesOut();
                    vector<Line>  in_A = fA->getIndicesIn();
                    vector<Line> out_B = fB->getIndicesOut();
                    vector<Line>  in_B = fB->getIndicesIn();
                    vector<Line> out_C = fC->getIndicesOut();
                    vector<Line>  in_C = fC->getIndicesIn();

                    vector<int> alpha_creation_A(A.ncreation.size(), 0);
                    vector<int> alpha_annihilation_A(A.ncreation.size(), 0);
                    vector<int> alpha_creation_B(B.ncreation.size(), 0);
                    vector<int> alpha_annihilation_B(B.ncreation.size(), 0);

                    for (vector<Line>::iterator i = out_A.begin();i != out_A.end();++i)
                        if (i->isAlpha()) alpha_creation_A[i->getType()]++;
                    for (vector<Line>::iterator i =  in_A.begin();i !=  in_A.end();++i)
                        if (i->isAlpha()) alpha_annihilation_A[i->getType()]++;
                    for (vector<Line>::iterator i = out_B.begin();i != out_B.end();++i)
                        if (i->isAlpha()) alpha_creation_B[i->getType()]++;
                    for (vector<Line>::iterator i =  in_B.begin();i !=  in_B.end();++i)
                        if (i->isAlpha()) alpha_annihilation_B[i->getType()]++;

                    vector<int> idx_A_(A.ndim);
                    {
                        int i = 0;
                        for (int j = 0;j < out_A.size();j++,i++)
                            idx_A_[i] = out_A[j].asInt();
                        for (int j = 0;j < in_A.size();j++,i++)
                            idx_A_[i] = in_A[j].asInt();
                    }

                    vector<int> idx_B_(B.ndim);
                    {
                        int i = 0;
                        for (int j = 0;j < out_B.size();j++,i++)
                            idx_B_[i] = out_B[j].asInt();
                        for (int j = 0;j < in_B.size();j++,i++)
                            idx_B_[i] = in_B[j].asInt();
                    }

                    vector<int> idx_C_(this->ndim);
                    {
                        int i = 0;
                        for (int j = 0;j < out_C.size();j++,i++)
                            idx_C_[i] = out_C[j].asInt();
                        for (int j = 0;j < in_C.size();j++,i++)
                            idx_C_[i] = in_C[j].asInt();
                    }

                    string idx_A__, idx_B__, idx_C__;
                    detail::conv_idx(idx_A_, idx_A__,
                                     idx_B_, idx_B__,
                                     idx_C_, idx_C__);

                    auto& tensor_A = A(alpha_creation_A, alpha_annihilation_A);
                    auto& tensor_B = B(alpha_creation_B, alpha_annihilation_B);

                    /*
                    cout << alpha << " " << beta[sc] << " " << *t << endl;
                    cout <<    tensor_A.getSymmetry() << " " <<   alpha_creation_A << " " <<   alpha_annihilation_A << endl;
                    cout <<    tensor_B.getSymmetry() << " " <<   alpha_creation_B << " " <<   alpha_annihilation_B << endl;
                    cout << scC.tensor->getSymmetry() << " " << scC.alpha_creation << " " << scC.alpha_annihilation << endl;
                    cout << idx_A__ << " " << idx_B__ << " " << idx_C__ << endl;

                    T before;
                    if (this->ndim == 0)
                    {
                        int64_t n;
                        before = (*scC.tensor)(0).getRawData(n)[0];
                    }
                    */

                    scC.tensor.mult(alpha*diagFactor, conja, tensor_A, idx_A__,
                                                      conjb, tensor_B, idx_B__,
                                            beta[sc],                  idx_C__);

                    /*
                    if (this->ndim == 0)
                    {
                        int64_t n;
                        cout << "? " << (*scC.tensor)(0).getRawData(n)[0]-before << endl;
                    }
                    */

                    beta[sc] = 1.0;
                }
            }
        }

        void sum(const Scalar& alpha, bool conja, const TensorImplementation<>& A_, const string& idx_A,
                 const Scalar& beta_,                                               const string& idx_B)
        {
            using namespace aquarius::autocc;

            auto& A = static_cast<const SpinBlockedTensor&>(A_);

            vector<Scalar> beta(this->cases.size(), beta_);

            for (int sc = 0;sc < this->cases.size();sc++)
            {
                SpinCase& scB = this->cases[sc];

                int nouttot_B = aquarius::sum(this->ncreation);

                string ext;
                vector<Line> lines_B_out(aquarius::sum(this->ncreation));
                vector<Line> lines_B_in(aquarius::sum(this->nannihilation));
                for (int i = 0, s = 0;s < this->ncreation.size();s++)
                {
                    for (int a = 0;a < scB.alpha_creation[s];a++,i++)
                    {
                        if (!contains(ext, idx_B[i])) ext += idx_B[i];
                        lines_B_out[i] = Line(idx_B[i], s, Line::VIRTUAL, Line::ALPHA);
                    }
                    for (int b = 0;b < this->ncreation[s]-scB.alpha_creation[s];b++,i++)
                    {
                        if (!contains(ext, idx_B[i])) ext += idx_B[i];
                        lines_B_out[i] = Line(idx_B[i], s, Line::VIRTUAL, Line::BETA);
                    }
                }
                for (int i = 0,s = 0;s < this->ncreation.size();s++)
                {
                    for (int a = 0;a < scB.alpha_annihilation[s];a++,i++)
                    {
                        if (!contains(ext, idx_B[i+nouttot_B])) ext += idx_B[i+nouttot_B];
                        lines_B_in[i] = Line(idx_B[i+nouttot_B], s, Line::VIRTUAL, Line::ALPHA);
                    }
                    for (int b = 0;b < this->nannihilation[s]-scB.alpha_annihilation[s];b++,i++)
                    {
                        if (!contains(ext, idx_B[i+nouttot_B])) ext += idx_B[i+nouttot_B];
                        lines_B_in[i] = Line(idx_B[i+nouttot_B], s, Line::VIRTUAL, Line::BETA);
                    }
                }

                int nouttot_A = aquarius::sum(A.ncreation);

                string sum;
                vector<Line> lines_A_out(sum(A.ncreation));
                vector<Line> lines_A_in(sum(A.nannihilation));
                vector<Line> lines_AandB_out, lines_AandB_in;
                for (int i = 0, s = 0;s < A.ncreation.size();s++)
                {
                    for (int a = 0;a < A.ncreation[s];a++,i++)
                    {
                        if (contains(ext, idx_A[i]))
                        {
                            int j; for (j = 0;idx_B[j] != idx_A[i];j++);
                            if (j < nouttot_B)
                            {
                                lines_A_out[i] = lines_B_out[j];
                                lines_AandB_out += lines_B_out[j];
                            }
                            else
                            {
                                lines_A_out[i] = lines_B_in[j-nouttot_B];
                                lines_AandB_out += lines_B_in[j-nouttot_B];
                            }
                        }
                        else
                        {
                            if (!contains(sum, idx_A[i])) sum += idx_A[i];
                            lines_A_out[i] = Line(idx_A[i], s, Line::VIRTUAL, Line::BETA);
                        }
                    }
                }
                for (int i = 0, s = 0;s < A.ncreation.size();s++)
                {
                    for (int a = 0;a < A.nannihilation[s];a++,i++)
                    {
                        if (contains(ext, idx_A[i+nouttot_A]))
                        {
                            int j; for (j = 0;idx_B[j] != idx_A[i+nouttot_A];j++);
                            if (j < nouttot_B)
                            {
                                lines_A_in[i] = lines_B_out[j];
                                lines_AandB_in += lines_B_out[j];
                            }
                            else
                            {
                                lines_A_in[i] = lines_B_in[j-nouttot_B];
                                lines_AandB_in += lines_B_in[j-nouttot_B];
                            }
                        }
                        else
                        {
                            if (!contains(sum, idx_A[i+nouttot_A])) sum += idx_A[i+nouttot_A];
                            lines_A_in[i] = Line(idx_A[i+nouttot_A], s, Line::VIRTUAL, Line::BETA);
                        }
                    }
                }

                vector<Line> lines_BnotA_out = lines_B_out;
                uniq(lines_AandB_out);
                uniq(lines_BnotA_out);
                exclude(lines_BnotA_out, lines_AandB_out);

                vector<Line> lines_BnotA_in = lines_B_in;
                uniq(lines_AandB_in);
                uniq(lines_BnotA_in);
                exclude(lines_BnotA_in, lines_AandB_in);

                Diagram d = Diagram(Diagram::SPINORBITAL,
                                    vec(Term(Diagram::SPINORBITAL)*
                                        Fragment("A", lines_A_out, lines_A_in)));

                for (int s = 0;s < max(A.ncreation.size(),this->ncreation.size());s++)
                {
                    vector<vector<Line> > assym(2);

                    for (vector<Line>::iterator i = lines_AandB_out.begin();i != lines_AandB_out.end();++i)
                        if (i->getType() == s) assym[0].push_back(*i);
                    for (vector<Line>::iterator i = lines_BnotA_out.begin();i != lines_BnotA_out.end();++i)
                        if (i->getType() == s) assym[1].push_back(*i);

                    if (!assym[0].empty() && !assym[1].empty()) d.antisymmetrize(assym);
                }

                for (int s = 0;s < max(A.ncreation.size(),this->ncreation.size());s++)
                {
                    vector<vector<Line> > assym(3);

                    for (vector<Line>::iterator i = lines_AandB_in.begin();i != lines_AandB_in.end();++i)
                        if (i->getType() == s) assym[0].push_back(*i);
                    for (vector<Line>::iterator i = lines_BnotA_in.begin();i != lines_BnotA_in.end();++i)
                        if (i->getType() == s) assym[1].push_back(*i);

                    if (!assym[0].empty() && !assym[1].empty()) d.antisymmetrize(assym);
                }

                d.convert(Diagram::UHF);

                /*
                 * Remove terms which are antisymmetrizations of same-spin groups
                 */
                for (int s = 0;s < max(A.ncreation.size(),this->ncreation.size());s++)
                {
                    for (int spin = 0;spin < 2;spin++)
                    {
                        vector<Term> terms = d.getTerms();
                        for (vector<Term>::iterator t1 = terms.begin();t1 != terms.end();++t1)
                        {
                            for (vector<Term>::iterator t2 = t1+1;t2 != terms.end();++t2)
                            {
                                if (Term(*t1).fixorder(filter_copy(t1->indices(), and1(isSpin(spin),isType(s)))) ==
                                    Term(*t2).fixorder(filter_copy(t2->indices(), and1(isSpin(spin),isType(s)))))
                                {
                                    d -= *t1;
                                    break;
                                }
                            }
                        }
                    }
                }

                d *= Term(Diagram::UHF)*Fragment("B", lines_B_out, lines_B_in);
                d.fixorder(true);

                for (vector<Term>::const_iterator t = d.getTerms().begin();t != d.getTerms().end();++t)
                {
                    double diagFactor = t->getFactor();

                    vector<Fragment>::const_iterator fA, fB;
                    for (vector<Fragment>::const_iterator f = t->getFragments().begin();f != t->getFragments().end();++f)
                    {
                        if (f->getOp() == "A") fA = f;
                        if (f->getOp() == "B") fB = f;
                    }

                    vector<Line> out_A = fA->getIndicesOut();
                    vector<Line>  in_A = fA->getIndicesIn();
                    vector<Line> out_B = fB->getIndicesOut();
                    vector<Line>  in_B = fB->getIndicesIn();

                    vector<int> alpha_creation_A(A.ncreation.size(), 0);
                    vector<int> alpha_annihilation_A(A.ncreation.size(), 0);

                    for (vector<Line>::iterator i = out_A.begin();i != out_A.end();++i)
                        if (i->isAlpha()) alpha_creation_A[i->getType()]++;
                    for (vector<Line>::iterator i =  in_A.begin();i !=  in_A.end();++i)
                        if (i->isAlpha()) alpha_annihilation_A[i->getType()]++;

                    vector<int> idx_A_(A.ndim);
                    {
                        int i = 0;
                        for (int j = 0;j < out_A.size();j++,i++)
                            idx_A_[i] = out_A[j].asInt();
                        for (int j = 0;j < in_A.size();j++,i++)
                            idx_A_[i] = in_A[j].asInt();
                    }

                    vector<int> idx_B_(this->ndim);
                    {
                        int i = 0;
                        for (int j = 0;j < out_B.size();j++,i++)
                            idx_B_[i] = out_B[j].asInt();
                        for (int j = 0;j < in_B.size();j++,i++)
                            idx_B_[i] = in_B[j].asInt();
                    }

                    string idx_A__, idx_B__;
                    detail::conv_idx(idx_A_, idx_A__,
                                     idx_B_, idx_B__);

                    //cout << alpha << " " << beta[sc] << " " << *t << endl;

                    auto& tensor_A = A(alpha_creation_A, alpha_annihilation_A);
                    scB.tensor.sum(alpha*diagFactor, conja, tensor_A, idx_A__,
                                           beta[sc],                  idx_B__);

                    beta[sc] = 1.0;
                }
            }
        }

        void scale(const Scalar& alpha, const string& idx_A)
        {
            for (int i = 0;i < idx_A.size();i++)
            {
                for (int j = i+1;j < idx_A.size();j++)
                {
                    assert(idx_A[i] == idx_A[j]);
                }
            }

            for (typename vector<SpinCase>::const_iterator scA = this->cases.begin();scA != this->cases.end();++scA)
            {
                scA->tensor.scale(alpha, idx_A);
            }
        }

        Scalar dot(bool conja, const TensorImplementation<>& A_, const string& idx_A,
                   bool conjb,                                   const string& idx_B) const
        {
            using namespace aquarius::autocc;

            auto& A = static_cast<const SpinBlockedTensor&>(A_);

            Scalar accum(this->F.type);

            int nouttot_A = aquarius::sum(A.ncreation);

            string sum;
            vector<Line> lines_A_out(sum(A.ncreation));
            vector<Line> lines_A_in(sum(A.nannihilation));
            for (int i = 0, s = 0;s < A.ncreation.size();s++)
            {
                for (int a = 0;a < A.ncreation[s];a++,i++)
                {
                    if (!contains(sum, idx_A[i])) sum += idx_A[i];
                    lines_A_out[i] = Line(idx_A[i], s, Line::VIRTUAL, Line::BETA);
                }
            }
            for (int i = 0, s = 0;s < A.ncreation.size();s++)
            {
                for (int a = 0;a < A.nannihilation[s];a++,i++)
                {
                    if (!contains(sum, idx_A[i+nouttot_A])) sum += idx_A[i+nouttot_A];
                    lines_A_in[i] = Line(idx_A[i+nouttot_A], s, Line::VIRTUAL, Line::BETA);
                }
            }

            int nouttot_B = sum(this->ncreation);

            vector<Line> lines_B_out(sum(this->ncreation));
            vector<Line> lines_B_in(sum(this->nannihilation));
            for (int i = 0, s = 0;s < this->ncreation.size();s++)
            {
                for (int a = 0;a < this->ncreation[s];a++,i++)
                {
                    if (!contains(sum, idx_B[i])) sum += idx_B[i];
                    lines_B_out[i] = Line(idx_B[i], s, Line::VIRTUAL, Line::BETA);
                }
            }
            for (int i = 0, s = 0;s < this->ncreation.size();s++)
            {
                for (int a = 0;a < this->nannihilation[s];a++,i++)
                {
                    if (!contains(sum, idx_B[i+nouttot_B])) sum += idx_B[i+nouttot_B];
                    lines_B_in[i] = Line(idx_B[i+nouttot_B], s, Line::VIRTUAL, Line::BETA);
                }
            }

            Diagram d = Diagram(Diagram::SPINORBITAL,
                                vec(Term(Diagram::SPINORBITAL)*
                                    Fragment("A", lines_A_out, lines_A_in)*
                                    Fragment("B", lines_B_out, lines_B_in)));

            d.convert(Diagram::UHF);

            /*
             * Remove terms which are antisymmetrizations of same-spin groups
             */
            for (int s = 0;s < max(A.ncreation.size(),this->ncreation.size());s++)
            {
                for (int spin = 0;spin < 2;spin++)
                {
                    vector<Term> terms = d.getTerms();
                    for (vector<Term>::iterator t1 = terms.begin();t1 != terms.end();++t1)
                    {
                        for (vector<Term>::iterator t2 = t1+1;t2 != terms.end();++t2)
                        {
                            if (Term(*t1).fixorder(filter_copy(t1->indices(), and1(isSpin(spin),isType(s)))) ==
                                Term(*t2).fixorder(filter_copy(t2->indices(), and1(isSpin(spin),isType(s)))))
                            {
                                d -= *t1;
                                break;
                            }
                        }
                    }
                }
            }

            for (vector<Term>::const_iterator t = d.getTerms().begin();t != d.getTerms().end();++t)
            {
                double diagFactor = t->getFactor();

                vector<Fragment>::const_iterator fA, fB;
                for (vector<Fragment>::const_iterator f = t->getFragments().begin();f != t->getFragments().end();++f)
                {
                    if (f->getOp() == "A") fA = f;
                    if (f->getOp() == "B") fB = f;
                }

                vector<Line> out_A = fA->getIndicesOut();
                vector<Line>  in_A = fA->getIndicesIn();
                vector<Line> out_B = fB->getIndicesOut();
                vector<Line>  in_B = fB->getIndicesIn();

                vector<int> alpha_creation_A(    A.ncreation.size(), 0);
                vector<int> alpha_annihilation_A (    A.ncreation.size(), 0);
                vector<int> alpha_creation_B(this->ncreation.size(), 0);
                vector<int> alpha_annihilation_B (this->ncreation.size(), 0);

                for (vector<Line>::iterator i = out_A.begin();i != out_A.end();++i)
                    if (i->isAlpha()) alpha_creation_A[i->getType()]++;
                for (vector<Line>::iterator i =  in_A.begin();i !=  in_A.end();++i)
                    if (i->isAlpha()) alpha_annihilation_A[i->getType()]++;
                for (vector<Line>::iterator i = out_B.begin();i != out_B.end();++i)
                    if (i->isAlpha()) alpha_creation_B[i->getType()]++;
                for (vector<Line>::iterator i =  in_B.begin();i !=  in_B.end();++i)
                    if (i->isAlpha()) alpha_annihilation_B[i->getType()]++;

                vector<int> idx_A_(A.ndim);
                {
                    int i = 0;
                    for (int j = 0;j < out_A.size();j++,i++)
                        idx_A_[i] = out_A[j].asInt();
                    for (int j = 0;j < in_A.size();j++,i++)
                        idx_A_[i] = in_A[j].asInt();
                }

                vector<int> idx_B_(this->ndim);
                {
                    int i = 0;
                    for (int j = 0;j < out_B.size();j++,i++)
                        idx_B_[i] = out_B[j].asInt();
                    for (int j = 0;j < in_B.size();j++,i++)
                        idx_B_[i] = in_B[j].asInt();
                }

                string idx_A__, idx_B__;
                detail::conv_idx(idx_A_, idx_A__,
                                 idx_B_, idx_B__);

                auto& tensor_A =       A(alpha_creation_A, alpha_annihilation_A);
                auto& tensor_B = (*this)(alpha_creation_B, alpha_annihilation_B);

                /*
                cout << alpha << " " << beta[sc] << " " << *t << endl;
                cout <<    tensor_A.getSymmetry() << " " <<   alpha_creation_A << " " <<   alpha_annihilation_A << endl;
                cout <<    tensor_B.getSymmetry() << " " <<   alpha_creation_B << " " <<   alpha_annihilation_B << endl;
                cout << scC.tensor->getSymmetry() << " " << scC.alpha_creation << " " << scC.alpha_annihilation << endl;
                cout << idx_A__ << " " << idx_B__ << " " << idx_C__ << endl;

                T before;
                if (this->ndim == 0)
                {
                    int64_t n;
                    before = (*scC.tensor)(0).getRawData(n)[0];
                }
                */

                accum += tensor_B.dot(conja, tensor_A, idx_A__,
                                      conjb,           idx_B__);

                /*
                if (this->ndim == 0)
                {
                    int64_t n;
                    cout << "? " << (*scC.tensor)(0).getRawData(n)[0]-before << endl;
                }
                */
            }

            return accum;
        }

        Scalar norm(int p) const
        {
            Scalar nrm(this->F.type);

            for (typename vector<SpinCase>::const_iterator sc = this->cases.begin();sc != this->cases.end();++sc)
            {
                double factor = 1;
                for (int s = 0;s < this->ncreation.size();s++)
                {
                    factor *= binom(    this->ncreation[s], sc->alpha_creation[s]);
                    factor *= binom(this->nannihilation[s],  sc->alpha_annihilation[s]);
                }

                Scalar subnrm(sc->tensor.norm(p));

                if (p == 2)
                {
                    nrm += factor*subnrm*subnrm;
                }
                else if (p == 0)
                {
                    nrm = max(nrm,subnrm);
                }
                else if (p == 1)
                {
                    nrm += factor*subnrm;
                }
                else
                {
                    assert(0);
                }
            }

            if (p == 2) nrm = sqrt(nrm);

            return nrm;
        }

        void slice(const Scalar& alpha, bool conja, const TensorImplementation<>& A, const vector<int>& start_A,
                   const Scalar&  beta,                                              const vector<int>& start_B,
                                                                                     const vector<int>& length)
        {
            //TODO
            assert(0);
        }

        void div(const Scalar& alpha, bool conja, const TensorImplementation<>& A_,
                                      bool conjb, const TensorImplementation<>& B_,
                 const Scalar& beta)
        {
            auto& A = static_cast<const SpinBlockedTensor&>(A_);
            auto& B = static_cast<const SpinBlockedTensor&>(B_);

            assert(this->nalpha_per_irrep == A.nalpha_per_irrep &&
                   this->nbeta_per_irrep  == A.nbeta_per_irrep);
            assert(this->nalpha_per_irrep == B.nalpha_per_irrep &&
                   this->nbeta_per_irrep  == B.nbeta_per_irrep);
            assert(this->ncreation     == A.ncreation &&
                   this->nannihilation == A.nannihilation);
            assert(this->ncreation     == B.ncreation &&
                   this->nannihilation == B.nannihilation);

            for (int i = 0;i < this->cases.size();i++)
            {
                this->cases[i].tensor.div(alpha, conja, A.cases[i].tensor,
                                                 conjb, B.cases[i].tensor,
                                          beta);
            }
        }

        void invert(const Scalar& alpha, bool conja, const TensorImplementation<>& A_,
                    const Scalar& beta)
        {
            auto& A = static_cast<const SpinBlockedTensor&>(A_);

            assert(this->nalpha_per_irrep == A.nalpha_per_irrep &&
                   this->nbeta_per_irrep  == A.nbeta_per_irrep);
            assert(this->ncreation     == A.ncreation &&
                   this->nannihilation == A.nannihilation);

            for (int i = 0;i < this->cases.size();i++)
            {
                this->cases[i].tensor.invert(alpha, conja, A.cases[i].tensor,
                                             beta);
            }
        }
};

}
}

#endif

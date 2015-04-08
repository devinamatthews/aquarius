#include "spinorbital_tensor.hpp"

using namespace aquarius::op;
using namespace aquarius::autocc;
using namespace aquarius::task;
using namespace aquarius::symmetry;

namespace aquarius
{
namespace tensor
{

static int conv_idx(const vector<int>& cidx_A, string& iidx_A)
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

static int conv_idx(const vector<int>& cidx_A, string& iidx_A,
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

static int conv_idx(const vector<int>& cidx_A, string& iidx_A,
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

template<class T>
map<const tCTF_World<T>*,map<const PointGroup*,pair<int,SpinorbitalTensor<T>*>>> SpinorbitalTensor<T>::scalars;

template<class T>
SpinorbitalTensor<T>::SpinorbitalTensor(const string& name, const SpinorbitalTensor<T>& t, const T val)
: IndexableCompositeTensor<SpinorbitalTensor<T>,SymmetryBlockedTensor<T>,T >(name, 0, 0),
  Distributed(t.arena), group(t.group), spin(0)
{
    cases.push_back(SpinCase());
    cases.back().construct(*this, group.totallySymmetricIrrep(), vector<int>(), vector<int>());
    *cases.back().tensor = val;
    register_scalar();
}

template<class T>
SpinorbitalTensor<T>::SpinorbitalTensor(const SpinorbitalTensor<T>& other)
: IndexableCompositeTensor<SpinorbitalTensor<T>,SymmetryBlockedTensor<T>,T >(other),
  Distributed(other.arena), group(other.group), spaces(other.spaces),
  nout(other.nout), nin(other.nin), spin(other.spin), cases(other.cases)
{
    assert(tensors.size() == cases.size());
    for (int i = 0;i < tensors.size();i++)
    {
        cases[i].tensor = tensors[i].tensor;
    }
    register_scalar();
}

template<class T>
SpinorbitalTensor<T>::SpinorbitalTensor(const string& name, const SpinorbitalTensor<T>& other)
: IndexableCompositeTensor<SpinorbitalTensor<T>,SymmetryBlockedTensor<T>,T >(name, other),
  Distributed(other.arena), group(other.group), spaces(other.spaces),
  nout(other.nout), nin(other.nin), spin(other.spin), cases(other.cases)
{
    assert(tensors.size() == cases.size());
    for (int i = 0;i < tensors.size();i++)
    {
        cases[i].tensor = tensors[i].tensor;
    }
    register_scalar();
}

template<class T>
SpinorbitalTensor<T>::SpinorbitalTensor(const string& name, const Arena& arena,
                                        const PointGroup& group,
                                        const vector<Space>& spaces,
                                        const vector<int>& nout,
                                        const vector<int>& nin, int spin)
: IndexableCompositeTensor<SpinorbitalTensor<T>,SymmetryBlockedTensor<T>,T>(name, aquarius::sum(nout)+aquarius::sum(nin), 0),
  Distributed(arena), group(group), spaces(spaces), nout(nout), nin(nin), spin(spin)
{
    int nspaces = spaces.size();
    int nouttot = aquarius::sum(nout);
    int nintot = aquarius::sum(nin);
    vector<int> whichout(nouttot), whichin(nintot);
    vector<int> alpha_out(nspaces), alpha_in(nspaces);

    for (int i = 0;i < nspaces;i++) assert(group == spaces[i].group);

    assert(abs(spin) <= nouttot+nintot);
    assert(abs(spin) >= abs(nouttot-nintot));
    assert(abs(spin)%2 == abs(nouttot-nintot)%2);

    for (int alphaout = 0;alphaout <= nouttot;alphaout++)
    {
        int alphain = alphaout + (nouttot-nintot-spin)/2;
        if (alphain < 0 || alphain > nintot) continue;

        fill(whichout.begin(), whichout.end(), 0);

        for (bool doneout = false;!doneout;)
        {
            fill(alpha_out.begin(), alpha_out.end(), 0);

            for (int i = 0;i < alphaout;i++)
            {
                alpha_out[whichout[i]]++;
            }

            fill(whichin.begin(), whichin.end(), 0);

            for (bool donein = false;!donein;)
            {
                fill(alpha_in.begin(), alpha_in.end(), 0);

                for (int i = 0;i < alphain;i++)
                {
                    alpha_in[whichin[i]]++;
                }

                bool ok = true;
                for (int i = 0;i < spaces.size();i++)
                {
                    if (alpha_out[i] > nout[i] ||
                        alpha_in[i] > nin[i])
                    {
                        ok = false;
                        break;
                    }
                }

                if (ok)
                {
                    cases.push_back(SpinCase());
                    cases.back().construct(*this, group.totallySymmetricIrrep(), alpha_out, alpha_in);
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

    register_scalar();
}

template<class T>
SpinorbitalTensor<T>::SpinorbitalTensor(const string& name, const Arena& arena,
                                        const PointGroup& group,
                                        const Representation& rep,
                                        const vector<Space>& spaces,
                                        const vector<int>& nout,
                                        const vector<int>& nin, int spin)
: IndexableCompositeTensor<SpinorbitalTensor<T>,SymmetryBlockedTensor<T>,T>(name, aquarius::sum(nout)+aquarius::sum(nin), 0),
  Distributed(arena), group(group), spaces(spaces), nout(nout), nin(nin), spin(spin)
{
    int nspaces = spaces.size();
    int nouttot = aquarius::sum(nout);
    int nintot = aquarius::sum(nin);
    vector<int> whichout(nouttot), whichin(nintot);
    vector<int> alpha_out(nspaces), alpha_in(nspaces);

    for (int i = 0;i < nspaces;i++) assert(group == spaces[i].group);

    assert(abs(spin) <= nouttot+nintot);
    assert(abs(spin) >= abs(nouttot-nintot));
    assert(abs(spin)%2 == abs(nouttot-nintot)%2);

    for (int alphaout = 0;alphaout <= nouttot;alphaout++)
    {
        int alphain = alphaout + (nouttot-nintot-spin)/2;
        if (alphain < 0 || alphain > nintot) continue;

        fill(whichout.begin(), whichout.end(), 0);

        for (bool doneout = false;!doneout;)
        {
            fill(alpha_out.begin(), alpha_out.end(), 0);

            for (int i = 0;i < alphaout;i++)
            {
                alpha_out[whichout[i]]++;
            }

            fill(whichin.begin(), whichin.end(), 0);

            for (bool donein = false;!donein;)
            {
                fill(alpha_in.begin(), alpha_in.end(), 0);

                for (int i = 0;i < alphain;i++)
                {
                    alpha_in[whichin[i]]++;
                }

                bool ok = true;
                for (int i = 0;i < spaces.size();i++)
                {
                    if (alpha_out[i] > nout[i] ||
                        alpha_in[i] > nin[i])
                    {
                        ok = false;
                        break;
                    }
                }

                if (ok)
                {
                    cases.push_back(SpinCase());
                    cases.back().construct(*this, rep, alpha_out, alpha_in);
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

    register_scalar();
}

template<class T>
SpinorbitalTensor<T>::~SpinorbitalTensor()
{
    unregister_scalar();
}

template<class T>
void SpinorbitalTensor<T>::SpinCase::construct(SpinorbitalTensor<T>& t,
                                               const Representation& rep,
                                               const vector<int>& alpha_out,
                                               const vector<int>& alpha_in)
{
    vector<vector<int>> len(t.ndim);
    vector<int> sym(t.ndim, AS);

    this->alpha_out = alpha_out;
    this->alpha_in = alpha_in;

    int i = 0;

    for (int s = 0;s < t.spaces.size();s++)
    {
        for (int a = 0;a < alpha_out[s];a++)
        {
            len[i++] = t.spaces[s].nalpha;
        }
        if (i > 0) sym[i-1] = NS;

        for (int b = 0;b < t.nout[s]-alpha_out[s];b++)
        {
            len[i++] = t.spaces[s].nbeta;
        }
        if (i > 0) sym[i-1] = NS;
    }

    for (int s = 0;s < t.spaces.size();s++)
    {
        for (int a = 0;a < alpha_in[s];a++)
        {
            len[i++] = t.spaces[s].nalpha;
        }
        if (i > 0) sym[i-1] = NS;

        for (int b = 0;b < t.nin[s]-alpha_in[s];b++)
        {
            len[i++] = t.spaces[s].nbeta;
        }
        if (i > 0) sym[i-1] = NS;
    }

    tensor = new SymmetryBlockedTensor<T>(t.name, t.arena, t.group, rep, t.ndim, len, sym, true);
    t.addTensor(tensor);
}

template<class T>
SymmetryBlockedTensor<T>& SpinorbitalTensor<T>::operator()(const vector<int>& alpha_out,
                                                           const vector<int>& alpha_in)
{
    return const_cast<SymmetryBlockedTensor<T>&>(const_cast<const SpinorbitalTensor<T>&>(*this)(alpha_out, alpha_in));
}

template<class T>
const SymmetryBlockedTensor<T>& SpinorbitalTensor<T>::operator()(const vector<int>& alpha_out,
                                                                 const vector<int>& alpha_in) const
{
    for (typename vector<SpinCase>::const_iterator sc = cases.begin();sc != cases.end();++sc)
    {
        if (sc->alpha_out == alpha_out &&
            sc->alpha_in  == alpha_in) return *(sc->tensor);
    }

    throw logic_error("spin case not found");
}

template<class T>
void SpinorbitalTensor<T>::mult(const T alpha, bool conja, const SpinorbitalTensor<T>& A, const string& idx_A,
                                               bool conjb, const SpinorbitalTensor<T>& B, const string& idx_B,
                                const T beta_,                                            const string& idx_C)
{
    assert(group == A.group);
    assert(group == B.group);
    assert(idx_A.size() == A.ndim);
    assert(idx_B.size() == B.ndim);
    assert(idx_C.size() == this->ndim);
    assert(spaces == A.spaces || this->ndim == 0 || A.ndim == 0);
    assert(spaces == B.spaces || this->ndim == 0 || B.ndim == 0);

    vector<T> beta(cases.size(), beta_);

    for (int sc = 0;sc < cases.size();sc++)
    {
        SpinCase& scC = cases[sc];

        int nouttot_C = aquarius::sum(nout);

        string ext;
        vector<Line> lines_C_out(aquarius::sum(nout));
        vector<Line> lines_C_in(aquarius::sum(nin));
        for (int i = 0, s = 0;s < spaces.size();s++)
        {
            for (int a = 0;a <         scC.alpha_out[s];a++,i++)
            {
                if (!contains(ext, idx_C[i])) ext += idx_C[i];
                lines_C_out[i] = Line(idx_C[i], s, Line::VIRTUAL, Line::ALPHA);
            }
            for (int b = 0;b < nout[s]-scC.alpha_out[s];b++,i++)
            {
                if (!contains(ext, idx_C[i])) ext += idx_C[i];
                lines_C_out[i] = Line(idx_C[i], s, Line::VIRTUAL, Line::BETA);
            }
        }
        for (int i = 0,s = 0;s < spaces.size();s++)
        {
            for (int a = 0;a <        scC.alpha_in[s];a++,i++)
            {
                if (!contains(ext, idx_C[i+nouttot_C])) ext += idx_C[i+nouttot_C];
                lines_C_in[i] = Line(idx_C[i+nouttot_C], s, Line::VIRTUAL, Line::ALPHA);
            }
            for (int b = 0;b < nin[s]-scC.alpha_in[s];b++,i++)
            {
                if (!contains(ext, idx_C[i+nouttot_C])) ext += idx_C[i+nouttot_C];
                lines_C_in[i] = Line(idx_C[i+nouttot_C], s, Line::VIRTUAL, Line::BETA);
            }
        }

        int nouttot_A = aquarius::sum(A.nout);

        string sum;
        vector<Line> lines_A_out(aquarius::sum(A.nout));
        vector<Line> lines_A_in(aquarius::sum(A.nin));
        vector<Line> lines_AandC_out, lines_AandC_in;
        for (int i = 0, s = 0;s < A.spaces.size();s++)
        {
            for (int a = 0;a < A.nout[s];a++,i++)
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
        for (int i = 0, s = 0;s < A.spaces.size();s++)
        {
            for (int a = 0;a < A.nin[s];a++,i++)
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

        int nouttot_B = aquarius::sum(B.nout);

        vector<Line> lines_B_out(aquarius::sum(B.nout));
        vector<Line> lines_B_in(aquarius::sum(B.nin));
        vector<Line> lines_BandC_out, lines_BandC_in;
        for (int i = 0, s = 0;s < B.spaces.size();s++)
        {
            for (int a = 0;a < B.nout[s];a++,i++)
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
        for (int i = 0, s = 0;s < B.spaces.size();s++)
        {
            for (int a = 0;a < B.nin[s];a++,i++)
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
                            {Term(Diagram::SPINORBITAL)*
                             Fragment("A", lines_A_out, lines_A_in)*
                             Fragment("B", lines_B_out, lines_B_in)});

        for (int s = 0;s < max(max(A.spaces.size(),B.spaces.size()),spaces.size());s++)
        {
            vector<vector<Line>> assym(3);

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

        for (int s = 0;s < max(max(A.spaces.size(),B.spaces.size()),spaces.size());s++)
        {
            vector<vector<Line>> assym(3);

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
        for (int s = 0;s < max(max(A.spaces.size(),B.spaces.size()),spaces.size());s++)
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

            vector<int> alpha_out_A(A.spaces.size(), 0);
            vector<int> alpha_in_A(A.spaces.size(), 0);
            vector<int> alpha_out_B(B.spaces.size(), 0);
            vector<int> alpha_in_B(B.spaces.size(), 0);

            for (vector<Line>::iterator i = out_A.begin();i != out_A.end();++i)
                if (i->isAlpha()) alpha_out_A[i->getType()]++;
            for (vector<Line>::iterator i =  in_A.begin();i !=  in_A.end();++i)
                if (i->isAlpha()) alpha_in_A[i->getType()]++;
            for (vector<Line>::iterator i = out_B.begin();i != out_B.end();++i)
                if (i->isAlpha()) alpha_out_B[i->getType()]++;
            for (vector<Line>::iterator i =  in_B.begin();i !=  in_B.end();++i)
                if (i->isAlpha()) alpha_in_B[i->getType()]++;

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
            conv_idx(idx_A_, idx_A__,
                     idx_B_, idx_B__,
                     idx_C_, idx_C__);

            const SymmetryBlockedTensor<T>& tensor_A = A(alpha_out_A, alpha_in_A);
            const SymmetryBlockedTensor<T>& tensor_B = B(alpha_out_B, alpha_in_B);

            if (0)
            {
                cout << alpha << " " << beta[sc] << " " << *t << endl;
                cout <<    tensor_A.getSymmetry() << " " <<   alpha_out_A << " " <<   alpha_in_A << endl;
                cout <<    tensor_B.getSymmetry() << " " <<   alpha_out_B << " " <<   alpha_in_B << endl;
                cout << scC.tensor->getSymmetry() << " " << scC.alpha_out << " " << scC.alpha_in << endl;
                cout << idx_A__ << " " << idx_B__ << " " << idx_C__ << endl;
            }

            scC.tensor->mult(alpha*diagFactor, conja, tensor_A, idx_A__,
                                               conjb, tensor_B, idx_B__,
                                     beta[sc],                  idx_C__);

            beta[sc] = 1.0;
        }
    }
}

template<class T>
void SpinorbitalTensor<T>::sum(const T alpha, bool conja, const SpinorbitalTensor<T>& A, const string& idx_A,
                               const T beta_,                                            const string& idx_B)
{
    assert(group == A.group);
    assert(idx_A.size() == A.ndim);
    assert(idx_B.size() == this->ndim);
    assert(spaces == A.spaces || this->ndim == 0 || A.ndim == 0);

    vector<T> beta(cases.size(), beta_);

    for (int sc = 0;sc < cases.size();sc++)
    {
        SpinCase& scB = cases[sc];

        int nouttot_B = aquarius::sum(this->nout);

        string ext;
        vector<Line> lines_B_out(aquarius::sum(nout));
        vector<Line> lines_B_in(aquarius::sum(nin));
        for (int i = 0, s = 0;s < spaces.size();s++)
        {
            for (int a = 0;a <         scB.alpha_out[s];a++,i++)
            {
                if (!contains(ext, idx_B[i])) ext += idx_B[i];
                lines_B_out[i] = Line(idx_B[i], s, Line::VIRTUAL, Line::ALPHA);
            }
            for (int b = 0;b < nout[s]-scB.alpha_out[s];b++,i++)
            {
                if (!contains(ext, idx_B[i])) ext += idx_B[i];
                lines_B_out[i] = Line(idx_B[i], s, Line::VIRTUAL, Line::BETA);
            }
        }
        for (int i = 0,s = 0;s < spaces.size();s++)
        {
            for (int a = 0;a <        scB.alpha_in[s];a++,i++)
            {
                if (!contains(ext, idx_B[i+nouttot_B])) ext += idx_B[i+nouttot_B];
                lines_B_in[i] = Line(idx_B[i+nouttot_B], s, Line::VIRTUAL, Line::ALPHA);
            }
            for (int b = 0;b < nin[s]-scB.alpha_in[s];b++,i++)
            {
                if (!contains(ext, idx_B[i+nouttot_B])) ext += idx_B[i+nouttot_B];
                lines_B_in[i] = Line(idx_B[i+nouttot_B], s, Line::VIRTUAL, Line::BETA);
            }
        }

        int nouttot_A = aquarius::sum(A.nout);

        string sum;
        vector<Line> lines_A_out(aquarius::sum(A.nout));
        vector<Line> lines_A_in(aquarius::sum(A.nin));
        vector<Line> lines_AandB_out, lines_AandB_in;
        for (int i = 0, s = 0;s < A.spaces.size();s++)
        {
            for (int a = 0;a < A.nout[s];a++,i++)
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
        for (int i = 0, s = 0;s < A.spaces.size();s++)
        {
            for (int a = 0;a < A.nin[s];a++,i++)
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
                            {Term(Diagram::SPINORBITAL)*
                             Fragment("A", lines_A_out, lines_A_in)});

        for (int s = 0;s < max(A.spaces.size(),spaces.size());s++)
        {
            vector<vector<Line>> assym(2);

            for (vector<Line>::iterator i = lines_AandB_out.begin();i != lines_AandB_out.end();++i)
                if (i->getType() == s) assym[0].push_back(*i);
            for (vector<Line>::iterator i = lines_BnotA_out.begin();i != lines_BnotA_out.end();++i)
                if (i->getType() == s) assym[1].push_back(*i);

            if (!assym[0].empty() && !assym[1].empty()) d.antisymmetrize(assym);
        }

        for (int s = 0;s < max(A.spaces.size(),spaces.size());s++)
        {
            vector<vector<Line>> assym(2);

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
        for (int s = 0;s < max(A.spaces.size(),spaces.size());s++)
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

            vector<int> alpha_out_A(A.spaces.size(), 0);
            vector<int> alpha_in_A(A.spaces.size(), 0);

            for (vector<Line>::iterator i = out_A.begin();i != out_A.end();++i)
                if (i->isAlpha()) alpha_out_A[i->getType()]++;
            for (vector<Line>::iterator i =  in_A.begin();i !=  in_A.end();++i)
                if (i->isAlpha()) alpha_in_A[i->getType()]++;

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
            conv_idx(idx_A_, idx_A__,
                     idx_B_, idx_B__);

            //cout << alpha << " " << beta[sc] << " " << *t << endl;

            const SymmetryBlockedTensor<T>& tensor_A = A(alpha_out_A, alpha_in_A);

            scB.tensor->sum(alpha*diagFactor, conja, tensor_A, idx_A__,
                                    beta[sc],                  idx_B__);

            beta[sc] = 1.0;
        }
    }
}

template<class T>
void SpinorbitalTensor<T>::scale(const T alpha, const string& idx_A)
{
    for (int i = 0;i < idx_A.size();i++)
    {
        for (int j = i+1;j < idx_A.size();j++)
        {
            //assert(idx_A[i] == idx_A[j]);
        }
    }

    for (typename vector<SpinCase>::const_iterator scA = cases.begin();scA != cases.end();++scA)
    {
        scA->tensor->scale(alpha);
    }
}

template<class T>
void SpinorbitalTensor<T>::weight(const vector<const vector<vector<T>>*>& da,
                                  const vector<const vector<vector<T>>*>& db,
                                  double shift)
{
    vector<const vector<vector<T>>*> d(this->ndim);

    for (typename vector<SpinCase>::iterator sc = cases.begin();sc != cases.end();++sc)
    {
        int i = 0;
        for (int s = 0;s < spaces.size();s++)
        {
            for (int a = 0;a <         sc->alpha_out[s];a++,i++) d[i] = da[s];
            for (int b = 0;b < nout[s]-sc->alpha_out[s];b++,i++) d[i] = db[s];
        }
        for (int s = 0;s < spaces.size();s++)
        {
            for (int a = 0;a <        sc->alpha_in[s];a++,i++) d[i] = da[s];
            for (int b = 0;b < nin[s]-sc->alpha_in[s];b++,i++) d[i] = db[s];
        }

        sc->tensor->weight(d, shift);
    }
}

template<class T>
T SpinorbitalTensor<T>::dot(bool conja, const SpinorbitalTensor<T>& A, const string& idx_A,
                            bool conjb,                                const string& idx_B) const
{
    SpinorbitalTensor<T>& sodt = scalar();
    sodt.mult(1, conja,     A, idx_A,
                 conjb, *this, idx_B,
              0,                  "");
    vector<T> vals;
    sodt(0)(0).getAllData(vals);
    assert(vals.size() == 1);
    return vals[0];
}

template<class T>
real_type_t<T> SpinorbitalTensor<T>::norm(int p) const
{
    real_type_t<T> nrm = 0;

    for (typename vector<SpinCase>::const_iterator sc = cases.begin();sc != cases.end();++sc)
    {
        double factor = 1;
        for (int s = 0;s < spaces.size();s++)
        {
            factor *= binom(nout[s], sc->alpha_out[s]);
            factor *= binom( nin[s],  sc->alpha_in[s]);
        }

        real_type_t<T> subnrm = sc->tensor->norm(p);

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
    }

    if (p == 2) nrm = sqrt(nrm);

    return nrm;
}


template <typename T>
void SpinorbitalTensor<T>::register_scalar()
{
    if (scalars.find(&arena.ctf<T>()) == scalars.end() ||
        scalars[&arena.ctf<T>()].find(&group) == scalars[&arena.ctf<T>()].end())
    {
        /*
         * If we are the first, make a new entry and put a new
         * scalar in it. The entry in scalars must be made FIRST,
         * since the new scalar will call this constructor too.
         */
        scalars[&arena.ctf<T>()][&group].first = -1;
        scalars[&arena.ctf<T>()][&group].second = new SpinorbitalTensor<T>("scalar", *this, (T)0);
    }

    scalars[&arena.ctf<T>()][&group].first++;
}

template <typename T>
void SpinorbitalTensor<T>::unregister_scalar()
{
    /*
     * The last tensor (besides the scalar in scalars)
     * will delete the entry, so if it does not exist
     * then we must be that scalar and nothing needs to be done.
     */
    if (scalars.find(&arena.ctf<T>()) == scalars.end() ||
        scalars[&arena.ctf<T>()].find(&group) == scalars[&arena.ctf<T>()].end()) return;

    if (--scalars[&arena.ctf<T>()][&group].first == 0)
    {
        /*
         * The entry must be deleted FIRST, so that the scalar
         * knows to do nothing.
         */
        SpinorbitalTensor<T>* scalar = scalars[&arena.ctf<T>()][&group].second;
        scalars[&arena.ctf<T>()].erase(&group);
        if (scalars[&arena.ctf<T>()].empty()) scalars.erase(&arena.ctf<T>());
        delete scalar;
    }
}

template <typename T>
SpinorbitalTensor<T>& SpinorbitalTensor<T>::scalar() const
{
    return *scalars[&arena.ctf<T>()][&group].second;
}

INSTANTIATE_SPECIALIZATIONS(SpinorbitalTensor);

}
}

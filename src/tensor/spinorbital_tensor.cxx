/* Copyin (c) 2013, Devin Matthews
 * All ins reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *      * Redistributions of source code must retain the above copyin
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyin
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

#include "spinorbital_tensor.hpp"

using namespace std;
using namespace aquarius;
using namespace aquarius::op;
using namespace aquarius::tensor;
using namespace aquarius::autocc;
using namespace aquarius::task;

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
SpinorbitalTensor<T>::SpinorbitalTensor(const SpinorbitalTensor<T>& t, const T val)
: IndexableCompositeTensor<SpinorbitalTensor<T>,DistTensor<T>,T >(0, 0),
  Resource(t.arena), spin(0)
{
    cases.push_back(SpinCase());
    cases.back().construct(*this, vector<int>(), vector<int>());
    *cases.back().tensor = val;
}

template<class T>
SpinorbitalTensor<T>::SpinorbitalTensor(const SpinorbitalTensor<T>& other)
: IndexableCompositeTensor<SpinorbitalTensor<T>,DistTensor<T>,T >(other),
  Resource(other.arena), spaces(other.spaces),
  nout(other.nout), nin(other.nin), spin(other.spin), cases(other.cases)
{
    assert(tensors.size() == cases.size());
    for (int i = 0;i < tensors.size();i++)
    {
        cases[i].tensor = tensors[i].tensor;
    }
}

template<class T>
SpinorbitalTensor<T>::SpinorbitalTensor(const Arena& arena,
                                        const vector<Space>& spaces,
                                        const vector<int>& nout,
                                        const vector<int>& nin, int spin)
: IndexableCompositeTensor<SpinorbitalTensor<T>,DistTensor<T>,T>(std::sum(nout)+std::sum(nin), 0),
  Resource(arena), spaces(spaces), nout(nout), nin(nin), spin(spin)
{
    int nspaces = spaces.size();
    int nouttot = std::sum(nout);
    int nintot = std::sum(nin);
    vector<int> whichout(nouttot), whichin(nintot);
    vector<int> alpha_out(nspaces), alpha_in(nspaces);

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
                    cases.back().construct(*this, alpha_out, alpha_in);
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

template<class T>
void SpinorbitalTensor<T>::SpinCase::construct(SpinorbitalTensor<T>& t,
                                               const vector<int>& alpha_out,
                                               const vector<int>& alpha_in)
{
    vector<int> len(t.ndim);
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

    tensor = new DistTensor<T>(t.arena, t.ndim, len, sym, true);
    t.addTensor(tensor);
}

/*
template<class T>
void SpinorbitalTensor<T>::matchTypes(const int nin_A, const int nout_A, const vector<Line>& log_A, const string& idx_A,
                                      const int nin_B, const int nout_B, const vector<Line>& log_B, const string& idx_B,
                                      vector<Line>& out_A, vector<Line>& in_A,
                                      vector<Line>& pout_A, vector<Line>& hout_A,
                                      vector<Line>& pin_A, vector<Line>& hin_A,
                                      vector<Line>& sum_A)
{
    for (int i = 0;i < nout_A;i++)
    {
        int j;
        for (j = 0;j < nout_B+nin_B;j++)
        {
            if (idx_A[i] == idx_B[j])
            {
                bool aisp = log_A[i].isParticle();
                bool bisp = log_B[j].isParticle();
                if (aisp != bisp)
                    throw logic_error("types do not match");
                out_A[i] = Line(idx_A[i], log_B[j].getType());
                if (log_B[j].isParticle())
                {
                    pin_A += Line(idx_A[i], log_B[j].getType());
                }
                else
                {
                    hin_A += Line(idx_A[i], log_B[j].getType());
                }
                break;
            }
        }
        if (j == nout_B+nin_B)
        {
            out_A[i] = Line(idx_A[i], log_A[i].getType());
            sum_A += Line(idx_A[i], log_A[i].getType());
        }
    }

    for (int i = nout_A;i < nout_A+nin_A;i++)
    {
        int j;
        for (j = 0;j < nout_B+nin_B;j++)
        {
            if (idx_A[i] == idx_B[j])
            {
                bool aisp = log_A[i].isParticle();
                bool bisp = log_B[j].isParticle();
                if (aisp != bisp)
                    throw logic_error("types do not match");
                in_A[i-nout_A] = Line(idx_A[i], log_B[j].getType());
                if (log_B[j].isParticle())
                {
                    pout_A += Line(idx_A[i], log_B[j].getType());
                }
                else
                {
                    hout_A += Line(idx_A[i], log_B[j].getType());
                }
                break;
            }
        }
        if (j == nout_B+nin_B)
        {
            in_A[i-nout_A] = Line(idx_A[i], log_A[i].getType());
            sum_A += Line(idx_A[i], log_A[i].getType());
        }
    }
}
*/

template<class T>
DistTensor<T>& SpinorbitalTensor<T>::operator()(const vector<int>& alpha_out,
                                                const vector<int>& alpha_in)
{
    return const_cast<DistTensor<T>&>(const_cast<const SpinorbitalTensor<T>&>(*this)(alpha_out, alpha_in));
}

template<class T>
const DistTensor<T>& SpinorbitalTensor<T>::operator()(const vector<int>& alpha_out,
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
    assert(idx_A.size() == A.ndim);
    assert(idx_B.size() == B.ndim);
    assert(idx_C.size() == this->ndim);
    assert(spaces == A.spaces || this->ndim == 0 || A.ndim == 0);
    assert(spaces == B.spaces || this->ndim == 0 || B.ndim == 0);

    vector<T> beta(cases.size(), beta_);

    for (int sc = 0;sc < cases.size();sc++)
    {
        SpinCase& scC = cases[sc];

        int nouttot_C = std::sum(nout);

        string ext;
        vector<Line> lines_C_out(std::sum(nout));
        vector<Line> lines_C_in(std::sum(nin));
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

        int nouttot_A = std::sum(A.nout);

        string sum;
        vector<Line> lines_A_out(std::sum(A.nout));
        vector<Line> lines_A_in(std::sum(A.nin));
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

        int nouttot_B = std::sum(B.nout);

        vector<Line> lines_B_out(std::sum(B.nout));
        vector<Line> lines_B_in(std::sum(B.nin));
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
                            vec(Term(Diagram::SPINORBITAL)*
                                Fragment("A", lines_A_out, lines_A_in)*
                                Fragment("B", lines_B_out, lines_B_in)));

        for (int s = 0;s < max(max(A.spaces.size(),B.spaces.size()),spaces.size());s++)
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

        for (int s = 0;s < max(max(A.spaces.size(),B.spaces.size()),spaces.size());s++)
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

            const DistTensor<T>& tensor_A = A(alpha_out_A, alpha_in_A);
            const DistTensor<T>& tensor_B = B(alpha_out_B, alpha_in_B);

            //cout << alpha << " " << beta[sc] << " " << *t << endl;
            //cout <<    tensor_A.getSymmetry() << " " << alpha_out_A << " " << alpha_in_A << endl;
            //cout <<    tensor_B.getSymmetry() << " " << alpha_out_B << " " << alpha_in_B << endl;
            //cout << scC.tensor->getSymmetry() << endl;

            scC.tensor->mult(alpha*diagFactor, conja, tensor_A, idx_A__,
                                               conjb, tensor_B, idx_B__,
                                     beta[sc],                  idx_C__);

            beta[sc] = 1.0;
        }
    }

//    /*
//    cout << "contracting: " << alpha << " * " << A.logical << "[";
//    for (int i = 0;i < A.ndim;i++) cout << idx_A[i] << ' ';
//    cout << "] " << B.logical << "[";
//    for (int i = 0;i < B.ndim;i++) cout << idx_B[i] << ' ';
//    cout << "] -> " << beta_ << " " << logical << "[";
//    for (int i = 0;i <   ndim;i++) cout << idx_C[i] << ' ';
//    cout << "]\n" << endl;
//    */
//
//    vector<int> idx_A_(idx_A.size());
//    vector<int> idx_B_(idx_B.size());
//    vector<int> idx_C_(idx_C.size());
//
//    int nouttot_A = sum(A.nout);
//    int nintot_A  = sum(A.nin);
//    int nouttot_B = sum(A.nout);
//    int nintot_B  = sum(A.nin);
//    int nouttot_C = sum(nout);
//    int nintot_C  = sum(nin);
//
//    vector<T> beta(cases.size(), beta_);
//
//    vector<Line> sAin(nintot_A);
//    vector<Line> sAout(nouttot_A);
//    vector<Line> sBin(nintot_B);
//    vector<Line> sBout(nouttot_B);
//    vector<Line> sCin(nintot_C);
//    vector<Line> sCout(nouttot_C);
//
//    for (typename vector<SpinCase>::const_iterator scC = cases.begin();scC != cases.end();++scC)
//    {
//        vector<Line> sum;
//        vector<Line> sum2;
//        vector<Line> apo, bpo;
//        vector<Line> aho, bho;
//        vector<Line> api, bpi;
//        vector<Line> ahi, bhi;
//
//        matchTypes(A.nA+A.nM, A.nE+A.nI,    A.logical, idx_A,
//                     nA+  nM,   nE+  nI, scC->logical, idx_C,
//                   sAout, sAin, apo, aho, api, ahi, sum);
//
//        matchTypes(B.nA+B.nM, B.nE+B.nI,    B.logical, idx_B,
//                     nA+  nM,   nE+  nI, scC->logical, idx_C,
//                   sBout, sBin, bpo, bho, bpi, bhi, sum);
//
//        uniq(sum);
//        for (int i = 0;i < sum.size();i++)
//        {
//            sum2 += sum[i].toAlpha();
//            sum2 += sum[i].toBeta();
//        }
//
//        {
//            int i;
//            for (i = 0;i < A.nA+A.nM;i++)
//            {
//                int j;
//                for (j = 0;j < B.nA+B.nM;j++)
//                {
//                    if (idx_A[i] == idx_B[j])
//                    {
//                        bool aisp = A.logical[i].isParticle();
//                        bool bisp = B.logical[j].isParticle();
//                        if (aisp != bisp)
//                                throw logic_error("types do not match");
//                        sAout[i] = sBout[j];
//                    }
//                }
//                for (;j < B.ndim;j++)
//                {
//                    if (idx_A[i] == idx_B[j])
//                    {
//                        bool aisp = A.logical[i].isParticle();
//                        bool bisp = B.logical[j].isParticle();
//                        if (aisp != bisp)
//                                throw logic_error("types do not match");
//                        sAout[i] = sBin[j-B.nA-B.nM];
//                    }
//                }
//            }
//            for (;i < A.ndim;i++)
//            {
//                int j;
//                for (j = 0;j < B.nA+B.nM;j++)
//                {
//                    if (idx_A[i] == idx_B[j])
//                    {
//                        bool aisp = A.logical[i].isParticle();
//                        bool bisp = B.logical[j].isParticle();
//                        if (aisp != bisp)
//                                throw logic_error("types do not match");
//                        sAin[i-A.nA-A.nM] = sBout[j];
//                    }
//                }
//                for (;j < B.ndim;j++)
//                {
//                    if (idx_A[i] == idx_B[j])
//                    {
//                        bool aisp = A.logical[i].isParticle();
//                        bool bisp = B.logical[j].isParticle();
//                        if (aisp != bisp)
//                                throw logic_error("types do not match");
//                        sAin[i-A.nA-A.nM] = sBin[j-B.nA-B.nM];
//                    }
//                }
//            }
//        }
//
//        for (int i = 0;i < nA+nM;i++)
//        {
//            sCout[i] = Line(idx_C[i], scC->logical[i].getType());
//        }
//
//        for (int i = 0;i < nE+nI;i++)
//        {
//            sCin[i] = Line(idx_C[nA+nM+i], scC->logical[nA+nM+i].getType());
//        }
//
//        Fragment uhfsA("T", sAout, sAin);
//        Fragment uhfsB("U", sBout, sBin);
//        Fragment uhfsC("V", sCout, sCin);
//
//        vector<Line> isect;
//
//        uniq(apo);
//        uniq(bpo);
//        isect = intersection(apo, bpo);
//        exclude(apo, isect);
//        exclude(bpo, isect);
//        uniq(aho);
//        uniq(bho);
//        isect = intersection(aho, bho);
//        exclude(aho, isect);
//        exclude(bho, isect);
//        uniq(api);
//        uniq(bpi);
//        isect = intersection(api, bpi);
//        exclude(api, isect);
//        exclude(bpi, isect);
//        uniq(ahi);
//        uniq(bhi);
//        isect = intersection(ahi, bhi);
//        exclude(ahi, isect);
//        exclude(bhi, isect);
//
//        //cout << apo << "," << bpo << "|" << aho << "," << bho << "|" <<
//        //             api << "," << bpi << "|" << ahi << "," << bhi << endl;
//
//        vector< vector<Line> > assym(2);
//
//        Diagram d(Diagram::UHF);
//        d += Term(Diagram::UHF)*uhfsA*uhfsB;
//        assym[0] = apo;
//        assym[1] = bpo;
//        if (!apo.empty() && !bpo.empty()) d.antisymmetrize(assym);
//        assym[0] = aho;
//        assym[1] = bho;
//        if (!aho.empty() && !bho.empty()) d.antisymmetrize(assym);
//        assym[0] = api;
//        assym[1] = bpi;
//        if (!api.empty() && !bpi.empty()) d.antisymmetrize(assym);
//        assym[0] = ahi;
//        assym[1] = bhi;
//        if (!ahi.empty() && !bhi.empty()) d.antisymmetrize(assym);
//        d.sum(sum);
//        d.fixorder(sum2);
//
//        /*
//         * Remove terms which are antisymmetrizations of same-spin groups
//         */
//        vector<Term> terms = d.getTerms();
//        for (vector<Term>::iterator t1 = terms.begin();t1 != terms.end();++t1)
//        {
//            for (vector<Term>::iterator t2 = t1+1;t2 != terms.end();++t2)
//            {
//                if (Term(*t1).fixorder(true) == Term(*t2).fixorder(true))
//                {
//                    d -= *t1;
//                    break;
//                }
//            }
//        }
//
//        terms = d.getTerms();
//        for (vector<Term>::iterator t = terms.begin();t != terms.end();++t)
//        {
//            *t *= uhfsC;
//
//            //cout << *t << endl;
//
//            double diagFactor = t->getFactor();
//
//            vector<Fragment>::iterator fA, fB, fC;
//            for (vector<Fragment>::iterator f = t->getFragments().begin();f != t->getFragments().end();++f)
//            {
//                if (f->getOp() == "T") fA = f;
//                if (f->getOp() == "U") fB = f;
//                if (f->getOp() == "V") fC = f;
//            }
//
//            vector<Line> out_A = fA->getIndicesOut();
//            vector<Line>  in_A = fA->getIndicesIn();
//            vector<Line> out_B = fB->getIndicesOut();
//            vector<Line>  in_B = fB->getIndicesIn();
//            vector<Line> out_C = fC->getIndicesOut();
//            vector<Line>  in_C = fC->getIndicesIn();
//
//            int nA_A = count_if(out_A.begin(), out_A.end(), isType<PARTICLE+ALPHA>());
//            int nM_A = count_if(out_A.begin(), out_A.end(), isType<    HOLE+ALPHA>());
//            int nE_A = count_if( in_A.begin(),  in_A.end(), isType<PARTICLE+ALPHA>());
//            int nI_A = count_if( in_A.begin(),  in_A.end(), isType<    HOLE+ALPHA>());
//
//            int nA_B = count_if(out_B.begin(), out_B.end(), isType<PARTICLE+ALPHA>());
//            int nM_B = count_if(out_B.begin(), out_B.end(), isType<    HOLE+ALPHA>());
//            int nE_B = count_if( in_B.begin(),  in_B.end(), isType<PARTICLE+ALPHA>());
//            int nI_B = count_if( in_B.begin(),  in_B.end(), isType<    HOLE+ALPHA>());
//
//            typename vector<SpinCase>::const_iterator scA = A.cases.end();
//            for (typename vector<SpinCase>::const_iterator sc = A.cases.begin();sc != A.cases.end();++sc)
//            {
//                if (nA_A == sc->nA && nM_A == sc->nM &&
//                    nE_A == sc->nE && nI_A == sc->nI) scA = sc;
//            }
//            if (scA == A.cases.end()) throw logic_error("no matching spin case for tensor A");
//
//            typename vector<SpinCase>::const_iterator scB = B.cases.end();
//            for (typename vector<SpinCase>::const_iterator sc = B.cases.begin();sc != B.cases.end();++sc)
//            {
//                if (nA_B == sc->nA && nM_B == sc->nM &&
//                    nE_B == sc->nE && nI_B == sc->nI) scB = sc;
//            }
//            if (scB == B.cases.end()) throw logic_error("no matching spin case for tensor B");
//
//            vector<Line> lA(A.ndim);
//            for (int i = 0;i < A.nA+A.nM;i++)
//            {
//                idx_A_[scA->log_to_phys[i]] = out_A[i].asInt();
//                lA[scA->log_to_phys[i]] = out_A[i];
//            }
//
//            for (int i = A.nA+A.nM;i < A.ndim;i++)
//            {
//                idx_A_[scA->log_to_phys[i]] = in_A[i-A.nA-A.nM].asInt();
//                lA[scA->log_to_phys[i]] = in_A[i-A.nA-A.nM];
//            }
//
//            vector<Line> lB(B.ndim);
//            for (int i = 0;i < B.nA+B.nM;i++)
//            {
//                idx_B_[scB->log_to_phys[i]] = out_B[i].asInt();
//                lB[scB->log_to_phys[i]] = out_B[i];
//            }
//
//            for (int i = B.nA+B.nM;i < B.ndim;i++)
//            {
//                idx_B_[scB->log_to_phys[i]] = in_B[i-B.nA-B.nM].asInt();
//                lB[scB->log_to_phys[i]] = in_B[i-B.nA-B.nM];
//            }
//
//            vector<Line> lC(ndim);
//            for (int i = 0;i < nA+nM;i++)
//            {
//                idx_C_[scC->log_to_phys[i]] = out_C[i].asInt();
//                lC[scC->log_to_phys[i]] = out_C[i];
//            }
//
//            for (int i = nA+nM;i < ndim;i++)
//            {
//                idx_C_[scC->log_to_phys[i]] = in_C[i-nA-nM].asInt();
//                lC[scC->log_to_phys[i]] = in_C[i-nA-nM];
//            }
//
//            /*
//            //if (A.logical == "aijk")
//            //{
//                cout << scA->log_to_phys << ' ' << scB->log_to_phys <<
//                        ' ' << scC->log_to_phys << endl;
//            cout << "A(" << lA << ") B(" << lB << ") C(" << lC << ") " <<
//                    alpha*scA->permFactor*scB->permFactor*scC->permFactor*diagFactor <<
//                    ' ' << beta[(int)(scC-cases.begin())] << endl << endl;
//            //}
//            */
//
//            string idx_A__, idx_iteratorB__, idx_C__;
//            conv_idx(idx_A_, idx_A__,
//                     idx_B_, idx_B__,
//                     idx_C_, idx_C__);
//
//            scC->tensor->mult(alpha*(T)(scA->permFactor*scB->permFactor*scC->permFactor*diagFactor),
//            //scC->tensor.mult(alpha*diagFactor,
//                              conja, *scA->tensor, idx_A__, conjb, *scB->tensor, idx_B__,
//                              beta[(int)(scC-cases.begin())], idx_C__);
//
//            beta[(int)(scC-cases.begin())] = 1.0;
//        }
//    }
}

template<class T>
void SpinorbitalTensor<T>::sum(const T alpha, bool conja, const SpinorbitalTensor<T>& A, const string& idx_A,
                               const T beta_,                                            const string& idx_B)
{
    assert(idx_A.size() == A.ndim);
    assert(idx_B.size() == this->ndim);
    assert(spaces == A.spaces || this->ndim == 0 || A.ndim == 0);

    vector<T> beta(cases.size(), beta_);

    for (int sc = 0;sc < cases.size();sc++)
    {
        SpinCase& scB = cases[sc];

        int nouttot_B = std::sum(this->nout);

        string ext;
        vector<Line> lines_B_out(std::sum(nout));
        vector<Line> lines_B_in(std::sum(nin));
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

        int nouttot_A = std::sum(A.nout);

        string sum;
        vector<Line> lines_A_out(std::sum(A.nout));
        vector<Line> lines_A_in(std::sum(A.nin));
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
                            vec(Term(Diagram::SPINORBITAL)*
                                Fragment("A", lines_A_out, lines_A_in)));

        for (int s = 0;s < max(A.spaces.size(),spaces.size());s++)
        {
            vector<vector<Line> > assym(2);

            for (vector<Line>::iterator i = lines_AandB_out.begin();i != lines_AandB_out.end();++i)
                if (i->getType() == s) assym[0].push_back(*i);
            for (vector<Line>::iterator i = lines_BnotA_out.begin();i != lines_BnotA_out.end();++i)
                if (i->getType() == s) assym[1].push_back(*i);

            if (!assym[0].empty() && !assym[1].empty()) d.antisymmetrize(assym);
        }

        for (int s = 0;s < max(A.spaces.size(),spaces.size());s++)
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

            const DistTensor<T>& tensor_A = A(alpha_out_A, alpha_in_A);

            scB.tensor->sum(alpha*diagFactor, conja, tensor_A, idx_A__,
                                    beta[sc],                  idx_B__);

            beta[sc] = 1.0;
        }
    }

//    vector<int> idx_A_(idx_A.size());
//    vector<int> idx_B_(idx_B.size());
//
//    //cout << "summing [";
//    //for (int i = 0;i < A.ndim;i++) cout << idx_A[i] << ' ';
//    //cout << "] [";
//    //for (int i = 0;i <   ndim;i++) cout << idx_B[i] << ' ';
//    //cout << "]" << endl;
//
//    vector<T> beta(cases.size(), beta_);
//
//    vector<Line> sAin(A.nA+A.nM);
//    vector<Line> sAout(A.nE+A.nI);
//
//    vector<Line> sBin(nA+nM);
//    vector<Line> sBout(nE+nI);
//
//    for (typename vector<SpinCase>::const_iterator scB = cases.begin();scB != cases.end();++scB)
//    {
//        vector<Line> sum;
//        vector<Line> sum2;
//        vector<Line> apo, aho, api, ahi;
//
//        matchTypes(A.nA+A.nM, A.nE+A.nI,    A.logical, idx_A,
//                     nA+  nM,   nE+  nI, scB->logical, idx_B,
//                   sAout, sAin, apo, aho, api, ahi, sum);
//
//        for (int i = 0;i < sum.size();i++)
//        {
//            sum2 += sum[i].toAlpha();
//            sum2 += sum[i].toBeta();
//        }
//
//        for (int i = 0;i < nA+nM;i++)
//        {
//            sBout[i] = Line(idx_B[i], scB->logical[i].getType());
//        }
//
//        for (int i = 0;i < nE+nI;i++)
//        {
//            sBin[i] = Line(idx_B[nA+nM+i], scB->logical[nA+nM+i].getType());
//        }
//
//        Fragment uhfsA("T", sAout, sAin);
//        Fragment uhfsB("U", sBout, sBin);
//
//        Diagram d(Diagram::UHF);
//        d += Term(Diagram::UHF)*uhfsA;
//        d.sum(sum);
//        d.fixorder(sum2);
//
//        vector<Term> terms = d.getTerms();
//        for (vector<Term>::iterator t = terms.begin();t != terms.end();++t)
//        {
//            *t *= uhfsB;
//
//            double diagFactor = t->getFactor();
//
//            //cout << *t << endl;
//
//            vector<Fragment>::iterator fA, fB;
//            for (vector<Fragment>::iterator f = t->getFragments().begin();f != t->getFragments().end();++f)
//            {
//                if (f->getOp() == "T") fA = f;
//                if (f->getOp() == "U") fB = f;
//            }
//
//            vector<Line> out_A = fA->getIndicesOut();
//            vector<Line>  in_A = fA->getIndicesIn();
//            vector<Line> out_B = fB->getIndicesOut();
//            vector<Line>  in_B = fB->getIndicesIn();
//
//            int nA_A = count_if(out_A.begin(), out_A.end(), isType<PARTICLE+ALPHA>());
//            int nM_A = count_if(out_A.begin(), out_A.end(), isType<    HOLE+ALPHA>());
//            int nE_A = count_if( in_A.begin(),  in_A.end(), isType<PARTICLE+ALPHA>());
//            int nI_A = count_if( in_A.begin(),  in_A.end(), isType<    HOLE+ALPHA>());
//
//            typename vector<SpinCase>::const_iterator scA = A.cases.end();
//            for (typename vector<SpinCase>::const_iterator sc = A.cases.begin();sc != A.cases.end();++sc)
//            {
//                //printf("%d %d %d %d - %d %d %d %d\n", nA_A, nM_A, nE_A, nI_A,
//                //       sc->nA, sc->nM, sc->nE, sc->nI);
//                if (nA_A == sc->nA && nM_A == sc->nM &&
//                    nE_A == sc->nE && nI_A == sc->nI) scA = sc;
//            }
//            if (scA == A.cases.end()) throw logic_error("no matching spin case for tensor A");
//
//            vector<Line> lA(A.ndim);
//            for (int i = 0;i < A.nA+A.nM;i++)
//            {
//                idx_A_[scA->log_to_phys[i]] = out_A[i].asInt();
//                lA[scA->log_to_phys[i]] = out_A[i];
//            }
//
//            for (int i = A.nA+A.nM;i < A.ndim;i++)
//            {
//                idx_A_[scA->log_to_phys[i]] = in_A[i-A.nA-A.nM].asInt();
//                lA[scA->log_to_phys[i]] = in_A[i-A.nA-A.nM];
//            }
//
//            vector<Line> lB(ndim);
//            for (int i = 0;i < nA+nM;i++)
//            {
//                idx_B_[scB->log_to_phys[i]] = out_B[i].asInt();
//                lB[scB->log_to_phys[i]] = out_B[i];
//            }
//
//            for (int i = nA+nM;i < ndim;i++)
//            {
//                idx_B_[scB->log_to_phys[i]] = in_B[i-nA-nM].asInt();
//                lB[scB->log_to_phys[i]] = in_B[i-nA-nM];
//            }
//
//            //cout << scA->log_to_phys << ' ' << scB->log_to_phys << endl;
//
//            //if (logical == "aijk")
//            //{
//            //cout << Fragment("A", lA, vector<Line>()) << ' ' <<
//            //             Fragment("B", lB, vector<Line>()) << ' ' <<
//            //             alpha << '*' << scA->permFactor << '*' <<
//            //             scB->permFactor << '*' << diagFactor <<
//            //             //alpha << '*' << diagFactor <<
//            //             ' ' << beta[(int)(scB-cases.begin())] << endl;
//            //}
//
//            string idx_A__, idx_B__;
//            conv_idx(idx_A_, idx_A__,
//                     idx_B_, idx_B__);
//
//            scB->tensor->sum(alpha*(T)(scA->permFactor*scB->permFactor*diagFactor),
//            //scB->tensor.sum(alpha*diagFactor,
//                              conja, *scA->tensor, idx_A__,
//                              beta[(int)(scB-cases.begin())], idx_B__);
//
//            beta[(int)(scB-cases.begin())] = 1.0;
//        }
//    }
}

template<class T>
void SpinorbitalTensor<T>::scale(const T alpha, const string& idx_A)
{
    for (int i = 0;i < idx_A.size();i++)
    {
        for (int j = i+1;j < idx_A.size();j++)
        {
            assert(idx_A[i] == idx_A[j]);
        }
    }

    for (typename vector<SpinCase>::const_iterator scA = cases.begin();scA != cases.end();++scA)
    {
        scA->tensor->scale(alpha);
    }
}

template<class T>
void SpinorbitalTensor<T>::weight(const vector<const vector<T>*>& da,
                                  const vector<const vector<T>*>& db)
{
    vector<const vector<T>*> d(this->ndim);

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

        sc->tensor->weight(d);
    }
}

template<class T>
T SpinorbitalTensor<T>::dot(bool conja, const SpinorbitalTensor<T>& A, const string& idx_A,
                            bool conjb,                                const string& idx_B) const
{
    SpinorbitalTensor<T> sodt(A, (T)0);
    sodt.mult(1, conja,     A, idx_A,
                 conjb, *this, idx_B,
              0,                  "");
    vector<T> vals;
    sodt(0).getAllData(vals);
    assert(vals.size() == 1);
    return vals[0];
}

INSTANTIATE_SPECIALIZATIONS(SpinorbitalTensor);

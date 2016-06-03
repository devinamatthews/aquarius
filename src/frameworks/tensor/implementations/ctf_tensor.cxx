#include "ctf_tensor.hpp"

namespace CTF_int { int cdealloc(void * ptr); }

namespace aquarius
{
namespace tensor
{

static const CTF_int::algstrct& getCTFRing(Field F)
{
    static CTF::Ring<      float, true>   sr;
    static CTF::Ring<     double, true>   dr;
    static CTF::Ring<   scomplex,false>  scr;
    static CTF::Ring<   dcomplex,false>  dcr;

    switch (F.type)
    {
        case Field::SINGLE:    return   sr;
        case Field::DOUBLE:    return   dr;
        case Field::SCOMPLEX:  return  scr;
        case Field::DCOMPLEX:  return  dcr;
    }

    return dr;
}

CTFTensor::CTFTensor(const INITIALIZER_TYPE(BOUNDED|IPSYMMETRIC|INDEXABLE|DISTRIBUTED|DIVISIBLE)& init)
: TensorImplementation<BOUNDED|IPSYMMETRIC|INDEXABLE|DISTRIBUTED|DIVISIBLE>(init)
{
    auto& ring = getCTFRing(F);
    ctf.reset(new CTF_int::tensor(&ring, ndim, len.data(), sym.data(), &arena.ctf()));

    strides.resize(ndim, 1);
    for (int i = 1;i < ndim;i++)
    {
        strides[i] = strides[i-1]*len[i-1];
    }

    /*
    auto it = scalars.find(make_pair(F.type, &arena.ctf()));

    if (it == scalars.end())
    {
        void* scalar;

        switch (F.type)
        {
            case Field::SINGLE:
                scalar = new tCTF_Tensor<float>(0, NULL, NULL, arena.ctf<float>(), "scalar");
                break;
            case Field::DOUBLE:
                scalar = new tCTF_Tensor<double>(0, NULL, NULL, arena.ctf<double>(), "scalar");
                break;
            case Field::SCOMPLEX:
                scalar = new tCTF_Tensor<scomplex>(0, NULL, NULL, arena.ctf<scomplex>(), "scalar");
                break;
            case Field::DCOMPLEX:
                scalar = new tCTF_Tensor<dcomplex>(0, NULL, NULL, arena.ctf<dcomplex>(), "scalar");
                break;
       }

       it = scalars.insert(make_pair(p, make_pair(0, scalar))).first;
    }

    it->second.first++;
    */
}

CTFTensor::~CTFTensor()
{
    /*
    pair<int,const void*> p(F.type, ctf_world);
    map<pair<int,const void*>, pair<int,void*> >::iterator it = scalars.find(p);
    assert(it != scalars.end());

    if (--it->second.first == 0)
    {
        switch (F.type)
        {
            case Field::SINGLE:
                delete static_cast<tCTF_Tensor<float>*>(it->second.second);
                break;
            case Field::DOUBLE:
                delete static_cast<tCTF_Tensor<double>*>(it->second.second);
                break;
            case Field::SCOMPLEX:
                delete static_cast<tCTF_Tensor<scomplex>*>(it->second.second);
                break;
            case Field::DCOMPLEX:
                delete static_cast<tCTF_Tensor<dcomplex>*>(it->second.second);
                break;
       }

       scalars.erase(it);
    }
    */
}

void CTFTensor::mult(const Scalar& alpha, bool conja, const TensorImplementation<>& A_, const string& idxA,
                                          bool conjb, const TensorImplementation<>& B_, const string& idxB,
                     const Scalar& beta,                                                const string& idxC)
{
    const CTFTensor& A = static_cast<const CTFTensor&>(A_);
    const CTFTensor& B = static_cast<const CTFTensor&>(B_);
    CTF_int::contraction(const_cast<CTF_int::tensor*>(A.ctf.get()), idxA.data(),
                         const_cast<CTF_int::tensor*>(B.ctf.get()), idxB.data(),
                         (const char*)alpha.data(), ctf.get(), idxC.data(), (const char*)beta.data()).execute();
}

void CTFTensor::sum(const Scalar& alpha, bool conja, const TensorImplementation<>& A_, const string& idxA,
                    const Scalar& beta,                                                const string& idxB)
{
    const CTFTensor& A = static_cast<const CTFTensor&>(A_);
    CTF_int::summation(const_cast<CTF_int::tensor*>(A.ctf.get()), idxA.data(),
                       (const char*)alpha.data(), ctf.get(), idxB.data(), (const char*)beta.data()).execute();
}

void CTFTensor::scale(const Scalar& alpha, const string& idxA)
{
    CTF_int::scaling(ctf.get(), idxA.data(), (const char*)alpha.data()).execute();
}

Scalar CTFTensor::dot(bool conja, const TensorImplementation<>& A_, const string& idxA,
                      bool conjb,                                   const string& idxB) const
{
    const CTFTensor& A = static_cast<const CTFTensor&>(A_);

    auto& ring = getCTFRing(F);
    CTF_int::tensor C(&ring, 0, NULL, NULL, &arena.ctf());
    Scalar alpha((Field::field)F.type, 1);
    Scalar beta((Field::field)F.type, 0);
    CTF_int::contraction(const_cast<CTF_int::tensor*>(A.ctf.get()), idxA.data(),
                         const_cast<CTF_int::tensor*>(  ctf.get()), idxB.data(),
                         (const char*)alpha.data(), &C, NULL, (const char*)beta.data()).execute();

    switch (F.type)
    {
        case Field::SINGLE:
        {
            CTF::Pair<float> p;
            C.read(1, (char*)&p);
            return Scalar(p.d);
        }
        case Field::DOUBLE:
        {
            CTF::Pair<double> p;
            C.read(1, (char*)&p);
            return Scalar(p.d);
        }
        case Field::SCOMPLEX:
        {
            CTF::Pair<scomplex> p;
            C.read(1, (char*)&p);
            return Scalar(p.d);
        }
        case Field::DCOMPLEX:
        {
            CTF::Pair<dcomplex> p;
            C.read(1, (char*)&p);
            return Scalar(p.d);
        }
    }

    return alpha;
}

template <typename T> void div_(key_type n, T alpha, bool conja, const CTF::Pair<T>* A,
                                                     bool conjb, const CTF::Pair<T>* B,
                                            T  beta,                   CTF::Pair<T>* C)
{
    if (conja) if (conjb)
        for (key_type i = 0;i < n;i++) C[i].d = beta*C[i].d + alpha*conj(A[i].d)/conj(B[i].d);
    else
        for (key_type i = 0;i < n;i++) C[i].d = beta*C[i].d + alpha*conj(A[i].d)/     B[i].d;
    else if (conjb)
        for (key_type i = 0;i < n;i++) C[i].d = beta*C[i].d + alpha*     A[i].d /conj(B[i].d);
    else
        for (key_type i = 0;i < n;i++) C[i].d = beta*C[i].d + alpha*     A[i].d /     B[i].d;
}

void CTFTensor::div(const Scalar& alpha, bool conja, const TensorImplementation<>& A_,
                                         bool conjb, const TensorImplementation<>& B_,
                    const Scalar& beta)
{
    const CTFTensor& A = static_cast<const CTFTensor&>(A_);
    const CTFTensor& B = static_cast<const CTFTensor&>(B_);

    const_cast<CTFTensor&>(A).ctf->align(ctf.get());
    const_cast<CTFTensor&>(B).ctf->align(ctf.get());

    int64_t nA, nB, nC;
    void *pA, *pB, *pC;
    A.ctf->read_local(&nA, (char**)&pA);
    B.ctf->read_local(&nB, (char**)&pB);
      ctf->read_local(&nC, (char**)&pC);
    assert(nA == nB && nA == nC);

    switch (F.type)
    {
        case Field::SINGLE:
            div_<float>(nA, (float)alpha, conja, (CTF::Pair<float>*)pA,
                                          conjb, (CTF::Pair<float>*)pB,
                            (float) beta,        (CTF::Pair<float>*)pC);
            break;
        case Field::DOUBLE:
            div_<double>(nA, (double)alpha, conja, (CTF::Pair<double>*)pA,
                                            conjb, (CTF::Pair<double>*)pB,
                             (double) beta,        (CTF::Pair<double>*)pC);
            break;
        case Field::SCOMPLEX:
            div_<scomplex>(nA, (scomplex)alpha, conja, (CTF::Pair<scomplex>*)pA,
                                                conjb, (CTF::Pair<scomplex>*)pB,
                               (scomplex) beta,        (CTF::Pair<scomplex>*)pC);
            break;
        case Field::DCOMPLEX:
            div_<dcomplex>(nA, (dcomplex)alpha, conja, (CTF::Pair<dcomplex>*)pA,
                                                conjb, (CTF::Pair<dcomplex>*)pB,
                               (dcomplex) beta,        (CTF::Pair<dcomplex>*)pC);
            break;
    }

    ctf->write(nC, NULL, NULL, (char*)pC);
    CTF_int::cdealloc(pA);
    CTF_int::cdealloc(pB);
    CTF_int::cdealloc(pC);
}

template <typename T> void invert_(key_type n, T alpha, bool conja, const CTF::Pair<T>* A,
                                               T  beta,                   CTF::Pair<T>* B)
{
    if (conja)
        for (key_type i = 0;i < n;i++) B[i].d = beta*B[i].d + alpha/conj(A[i].d);
    else
        for (key_type i = 0;i < n;i++) B[i].d = beta*B[i].d + alpha/     A[i].d;
}

void CTFTensor::invert(const Scalar& alpha, bool conja, const TensorImplementation<>& A_,
                       const Scalar& beta)
{
    const CTFTensor& A = static_cast<const CTFTensor&>(A_);

    const_cast<CTFTensor&>(A).ctf->align(ctf.get());

    int64_t nA, nB;
    void *pA, *pB;
    A.ctf->read_local(&nA, (char**)&pA);
      ctf->read_local(&nB, (char**)&pB);
    assert(nA == nB);

    switch (F.type)
    {
        case Field::SINGLE:
            invert_<float>(nA, (float)alpha, conja, (CTF::Pair<float>*)pA,
                               (float) beta,        (CTF::Pair<float>*)pB);
            break;
        case Field::DOUBLE:
            invert_<double>(nA, (double)alpha, conja, (CTF::Pair<double>*)pA,
                                (double) beta,        (CTF::Pair<double>*)pB);
            break;
        case Field::SCOMPLEX:
            invert_<scomplex>(nA, (scomplex)alpha, conja, (CTF::Pair<scomplex>*)pA,
                                  (scomplex) beta,        (CTF::Pair<scomplex>*)pB);
            break;
        case Field::DCOMPLEX:
            invert_<dcomplex>(nA, (dcomplex)alpha, conja, (CTF::Pair<dcomplex>*)pA,
                                  (dcomplex) beta,        (CTF::Pair<dcomplex>*)pB);
            break;
    }

    ctf->write(nB, NULL, NULL, (char*)pB);
    CTF_int::cdealloc(pA);
    CTF_int::cdealloc(pB);
}

const vector<key_type>& CTFTensor::getKeyStrides() const
{
    return strides;
}

void CTFTensor::getAllKeys(KeyVector& keys) const
{
    keys.clear();

    //TODO: precompute size

    key_type key = 0;
    vector<int> idx(ndim);
    for (int i = 0;i < ndim;i++)
    {
        if (i < ndim-1 && (sym[i] == AS || sym[i] == SH))
        {
            idx[i+1] = idx[i]+1;
        }

        key += idx[i]*strides[i];
    }

    for (bool done = false;!done;)
    {
        keys.push_back(key);

        for (int i = 0;i < ndim;i++)
        {
            if (idx[i] == len[i]-1 ||
                (i < ndim-1 && (((sym[i] == AS || sym[i] == SH) && idx[i] == idx[i+1]-1) ||
                                 (sym[i] == SY                  && idx[i] == idx[i+1]))))
            {
                int bottom = 0;
                if (i > 0 && (sym[i-1] == AS || sym[i-1] == SH))
                {
                    bottom = idx[i-1]+1;
                }
                else if (i > 0 && sym[i-1] == SY)
                {
                    bottom = idx[i-1];
                }

                key -= (idx[i]-bottom)*strides[i];
                idx[i] = bottom;

                if (i == ndim-1) done = true;
            }
            else
            {
                idx[i]++;
                key += strides[i];
                break;
            }
        }
        if (ndim == 0) done = true;
    }
}

void CTFTensor::getAllData(KeyValueVector& kv) const
{
    //TODO: precompute size and avoid copy

    int64_t n;
    void* d;
    ctf->allread(&n, (char**)&d, false);
    kv.resize(n);

    switch (F.type)
    {
        case Field::SINGLE:    copy_n((      float*)d, n, kv.data<      float>()); break;
        case Field::DOUBLE:    copy_n((     double*)d, n, kv.data<     double>()); break;
        case Field::SCOMPLEX:  copy_n((   scomplex*)d, n, kv.data<   scomplex>()); break;
        case Field::DCOMPLEX:  copy_n((   dcomplex*)d, n, kv.data<   dcomplex>()); break;
    }

    CTF_int::cdealloc(d);

    getAllKeys(kv.keys());
    assert(kv.keys().size() == n);
}

void CTFTensor::getLocalKeys(KeyVector& keys) const
{
    int64_t n;
    void* p;
    ctf->read_local(&n, (char**)&p);
    keys.resize(n);

    switch (F.type)
    {
        case Field::SINGLE:
            for (key_type i = 0;i < n;i++) keys[i] = ((CTF::Pair<      float>*)p)[i].k;
            break;
        case Field::DOUBLE:
            for (key_type i = 0;i < n;i++) keys[i] = ((CTF::Pair<     double>*)p)[i].k;
            break;
        case Field::SCOMPLEX:
            for (key_type i = 0;i < n;i++) keys[i] = ((CTF::Pair<   scomplex>*)p)[i].k;
            break;
        case Field::DCOMPLEX:
            for (key_type i = 0;i < n;i++) keys[i] = ((CTF::Pair<   dcomplex>*)p)[i].k;
            break;
    }

    CTF_int::cdealloc(p);
}

void CTFTensor::getLocalData(KeyValueVector& kv) const
{
    getLocalKeys(kv.keys());
    kv.resize(kv.keys().size());
    getRemoteData(kv.keys().size(), kv.keys().data(), kv.data<void>());
}

void CTFTensor::getRemoteData(key_type n, key_type* keys, void* values) const
{
    switch (F.type)
    {
        case Field::SINGLE:
        {
            vector<CTF::Pair<float>> p(n);
            for (key_type i = 0;i < n;i++) p[i].k = keys[i];
            ctf->read(n, (char*)p.data());
            for (key_type i = 0;i < n;i++) { keys[i] = p[i].k; ((float*)values)[i] = p[i].d; };
        }
        break;
        case Field::DOUBLE:
        {
            vector<CTF::Pair<double>> p(n);
            for (key_type i = 0;i < n;i++) p[i].k = keys[i];
            ctf->read(n, (char*)p.data());
            for (key_type i = 0;i < n;i++) { keys[i] = p[i].k; ((double*)values)[i] = p[i].d; };
        }
        break;
        case Field::SCOMPLEX:
        {
            vector<CTF::Pair<scomplex>> p(n);
            for (key_type i = 0;i < n;i++) p[i].k = keys[i];
            ctf->read(n, (char*)p.data());
            for (key_type i = 0;i < n;i++) { keys[i] = p[i].k; ((scomplex*)values)[i] = p[i].d; };
        }
        break;
        case Field::DCOMPLEX:
        {
            vector<CTF::Pair<dcomplex>> p(n);
            for (key_type i = 0;i < n;i++) p[i].k = keys[i];
            ctf->read(n, (char*)p.data());
            for (key_type i = 0;i < n;i++) { keys[i] = p[i].k; ((dcomplex*)values)[i] = p[i].d; };
        }
        break;
    }
}

void CTFTensor::getRemoteData() const
{
    getRemoteData(0, NULL, NULL);
}

void CTFTensor::setRemoteData(key_type n, const key_type* keys, const void* values)
{
    addRemoteData(n, 1.0, keys, values, 0.0);
}

void CTFTensor::setRemoteData()
{
    setRemoteData(0, NULL, NULL);
}

void CTFTensor::addRemoteData(key_type n, const Scalar& alpha, const key_type* keys,
                              const void* values, const Scalar& beta)
{
    switch (F.type)
    {
        case Field::SINGLE:
        {
            vector<CTF::Pair<float>> p(n);
            for (key_type i = 0;i < n;i++) { p[i].k = keys[i]; p[i].d = ((float*)values)[i]; };
            ctf->write(n, (const char*)alpha.data(), (const char*)beta.data(), (char*)p.data());
        }
        break;
        case Field::DOUBLE:
        {
            vector<CTF::Pair<double>> p(n);
            for (key_type i = 0;i < n;i++) { p[i].k = keys[i]; p[i].d = ((double*)values)[i]; };
            ctf->write(n, (const char*)alpha.data(), (const char*)beta.data(), (char*)p.data());
        }
        break;
        case Field::SCOMPLEX:
        {
            vector<CTF::Pair<scomplex>> p(n);
            for (key_type i = 0;i < n;i++) { p[i].k = keys[i]; p[i].d = ((scomplex*)values)[i]; };
            ctf->write(n, (const char*)alpha.data(), (const char*)beta.data(), (char*)p.data());
        }
        break;
        case Field::DCOMPLEX:
        {
            vector<CTF::Pair<dcomplex>> p(n);
            for (key_type i = 0;i < n;i++) { p[i].k = keys[i]; p[i].d = ((dcomplex*)values)[i]; };
            ctf->write(n, (const char*)alpha.data(), (const char*)beta.data(), (char*)p.data());
        }
        break;
    }
}

void CTFTensor::addRemoteData(const Scalar& alpha, const Scalar& beta)
{
    addRemoteData(0, alpha, NULL, NULL, beta);
}

Scalar CTFTensor::norm(int p) const
{
    Scalar nrm((Field::field)F.type);

    switch (p)
    {
        case 00:
            switch (F.type)
            {
                case Field::SINGLE:
                {
                    CTF::Monoid<float,true> mmax(nrm.to<float>(), CTF_int::default_max<float,true>, MPI_MAX);
                    ctf->reduce_sumabs((char*)nrm.data(), &mmax);
                }
                break;
                case Field::DOUBLE:
                {
                    CTF::Monoid<double,true> mmax(nrm.to<double>(), CTF_int::default_max<double,true>, MPI_MAX);
                    ctf->reduce_sumabs((char*)nrm.data(), &mmax);
                }
                break;
                case Field::SCOMPLEX:
                {
                    CTF::Monoid<scomplex,true> mmax(nrm.to<scomplex>(), CTF_int::default_max<scomplex,true>, MPI_MAX);
                    ctf->reduce_sumabs((char*)nrm.data(), &mmax);
                }
                break;
                case Field::DCOMPLEX:
                {
                    CTF::Monoid<dcomplex,true> mmax(nrm.to<dcomplex>(), CTF_int::default_max<dcomplex,true>, MPI_MAX);
                    ctf->reduce_sumabs((char*)nrm.data(), &mmax);
                }
                break;
            }
            break;
        case 1:
            ctf->reduce_sumabs((char*)nrm.data());
            break;
        case 2:
            ctf->reduce_sumsq((char*)nrm.data());
            nrm = sqrt(nrm);
            break;
        default:
            assert(0);
            break;
    }

    return nrm;
}

void CTFTensor::slice(const Scalar& alpha, bool conja, const vector<int>& start_A, const TensorImplementation<>& A_,
                      const Scalar&  beta,             const vector<int>& start_B, const vector<int>& len)
{
    auto& A = static_cast<const CTFTensor&>(A_);

    vector<int> end_A(start_A);
    vector<int> end_B(start_B);

    for (int i = 0;i < ndim;i++)
    {
        end_A[i] += len[i];
        end_B[i] += len[i];
    }

    ctf->slice(             start_B.data(), end_B.data(), (const char*)beta.data(),
               A.ctf.get(), start_A.data(), end_B.data(), (const char*)alpha.data());
}

}
}

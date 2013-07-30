/* Copyright (c) 2013, Devin Matthews
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following
 * conditions are met:
 *      * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *      * Redistributions in binary form must reproduce the above copyright
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

#include <iostream>

#include "stl_ext/stl_ext.hpp"
#include "util/util.h"

using namespace std;
//using namespace aquarius::tensor;
//using namespace aquarius::slide;

template <class T, class U>
struct if_exists
{
    typedef U type;
};

template <class Derived, class T> class Tensor;
template <class Derived, class T> class ScaledTensor;
template <class Derived, class T> class TensorMult;

template <class Derived, typename T>
class Tensor
{
    public:
        typedef T dtype;

        virtual ~Tensor() {}

        Derived& getDerived() { return static_cast<Derived&>(*this); }

        const Derived& getDerived() const { return static_cast<const Derived&>(*this); }

        //template <typename cvDerived> typename if_exists<typename cvDerived::dtype, TensorMult<Derived,T> >::type
        ENABLE_IF_SAME(Derived,cvDerived,CONCAT(TensorMult<Derived,T>))
        operator*(const cvDerived& other) const
        {
            return TensorMult<Derived,T>(ScaledTensor<const Derived,T>(getDerived(), (T)1),
                                         ScaledTensor<const Derived,T>(other.getDerived(), (T)1));
        }

        virtual T dot(const Derived& A) const = 0;
};

template <class Derived, typename T>
class ScaledTensor
{
    public:
        Derived& tensor_;
        T factor_;

        template <typename cvDerived>
        ScaledTensor(const ScaledTensor<cvDerived,T>& other)
        : tensor_(other.tensor_), factor_(other.factor_) {}

        ScaledTensor(Derived& tensor, const T factor)
        : tensor_(tensor), factor_(factor) {}
};

template <class Derived, typename T>
class TensorMult
{
    private:
        const TensorMult& operator=(const TensorMult<Derived,T>& other);

    public:
        ScaledTensor<const Derived,T> A_;
        ScaledTensor<const Derived,T> B_;
        T factor_;

        template <class Derived1, class Derived2>
        TensorMult(const ScaledTensor<Derived1,T>& A, const ScaledTensor<Derived2,T>& B)
        : A_(A), B_(B), factor_(B.factor_) {}
};

template <class Derived, typename T>
T scalar(const TensorMult<Derived,T>& tm)
{
    return tm.factor_*tm.B_.tensor_.dot(tm.A_.tensor_);
}

template <typename T>
class DistTensor : public Tensor<DistTensor<T>,T>
{
    public:
        virtual ~DistTensor() {}

        T dot(const DistTensor<T>& A) const { return (T)0; }
};

int main(int argc, char **argv)
{
    //MPI::Init(argc, argv);
    //SLIDE::init();
    #ifdef USE_ELEMENTAL
    elem::Initialize(argc, argv);
    #endif

    DistTensor<double> Delta;

    TensorMult<DistTensor<double>,double> tm = Delta*Delta;
    cout << &tm << endl;
    cout << tm.factor_ << endl;
    cout << &tm.A_ << endl;
    cout << &tm.A_.tensor_ << endl;
    cout << &tm.B_ << endl;
    cout << &tm.B_.tensor_ << endl;
    double S2 = scalar(tm);

    #ifdef USE_ELEMENTAL
    elem::Finalize();
    #endif
    //SLIDE::finish();
    //MPI::Finalize();
}

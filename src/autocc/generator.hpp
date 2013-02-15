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

#ifndef _AQUARIUS_AUTOCC_GENERATOR_HPP_
#define _AQUARIUS_AUTOCC_GENERATOR_HPP_

#include <vector>
#include <utility>

#include "line.hpp"

namespace aquarius
{
namespace autocc
{

class Operator;
class BasicOperator;
class ScaleOperator;
class ExponentialOperator;
class OperatorProduct;
class OperatorSum;

class Generator;
class ManifoldGenerator;
class BasicManifoldGenerator;
class ProductManifoldGenerator;
class SumManifoldGenerator;

class Generator
{
    protected:
        int _line;
        Generator() : _line(0) {}
};

class ManifoldGenerator : public Generator
{
    protected:
        Manifold leftMin, leftMax;
        Manifold rightMin, rightMax;

    public:
        ManifoldGenerator(const Manifold& leftMin, const Manifold& leftMax,
                          const Manifold& rightMin, const Manifold& rightMax)
        : leftMin(leftMin), leftMax(leftMax), rightMin(rightMin), rightMax(rightMax) {}

        /**
         * Return the next set of left and right Manifolds such that
         * left_i > left_i-1 || (left_i == left_i-1 && right_i > right_i-1)
         */
        virtual bool next(Manifold& left, Manifold& right) = 0;

        virtual ~ManifoldGenerator() {}
};

class BasicManifoldGenerator : public ManifoldGenerator
{
    friend class BasicOperator;

    private:
        const BasicOperator& op;
        std::vector< std::pair<Manifold,Manifold> > choices;
        std::vector< std::pair<Manifold,Manifold> >::iterator it;
        BasicManifoldGenerator(BasicManifoldGenerator& other);
        BasicManifoldGenerator& operator=(BasicManifoldGenerator& other);

    protected:
        BasicManifoldGenerator(const Manifold& leftMin, const Manifold& leftMax,
                               const Manifold& rightMin, const Manifold& rightMax,
                               const BasicOperator& op)
        : ManifoldGenerator(leftMin, leftMax, rightMin, rightMax), op(op) {}

    public:
        virtual bool next(Manifold& left, Manifold& right);
};

class ProductManifoldGenerator : public ManifoldGenerator
{
    friend class OperatorProduct;

    private:
        const OperatorProduct& op;
        ManifoldGenerator *lGen, *rGen;
        Manifold ll, lr, rl, rr, c;
        bool _ldone, _lhdone, _rdone, _rhdone, cdone, chdone;
        Manifold _l, _r, __l;
        ProductManifoldGenerator(ProductManifoldGenerator& other);
        ProductManifoldGenerator& operator=(ProductManifoldGenerator& other);

    protected:
        ProductManifoldGenerator(const Manifold& leftMin, const Manifold& leftMax,
                                 const Manifold& rightMin, const Manifold& rightMax,
                                 const OperatorProduct& op)
        : ManifoldGenerator(leftMin, leftMax, rightMin, rightMax), op(op) {}

    public:
        virtual bool next(Manifold& left, Manifold& right);
};

class SumManifoldGenerator : public ManifoldGenerator
{
    friend class OperatorSum;

    private:
        const OperatorSum& op;
        bool lHasMore, rHasMore;
        ManifoldGenerator *lGenerator, *rGenerator;
        std::pair<Manifold,Manifold> lNext, rNext;
        SumManifoldGenerator(SumManifoldGenerator& other);
        SumManifoldGenerator& operator=(SumManifoldGenerator& other);

    protected:
        SumManifoldGenerator(const Manifold& leftMin, const Manifold& leftMax,
                             const Manifold& rightMin, const Manifold& rightMax,
                             const OperatorSum& op)
        : ManifoldGenerator(leftMin, leftMax, rightMin, rightMax), op(op) {}

    public:
        virtual bool next(Manifold& left, Manifold& right);
};

}
}

/*
 * Generator macros derived from http://www.codeproject.com/Articles/29524/Generators-in-C
 */
#define generator_next(...) next(__VA_ARGS__) { switch(_line) { case 0:;

#define generator_stop  } _line = __LINE__; case __LINE__:; } return false

#define generator_yield(...) \
        do { \
            _line=__LINE__; \
            __VA_ARGS__; return true; case __LINE__:; \
        } while (0)

#endif

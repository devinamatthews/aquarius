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

#ifndef _AQUARIUS_AUTOCC_OPERATOR_HPP_
#define _AQUARIUS_AUTOCC_OPERATOR_HPP_

#include <vector>
#include <string>
#include <climits>

#include "diagram.hpp"
#include "fraction.hpp"
#include "generator.hpp"

namespace aquarius
{
namespace autocc
{

class Term;
class Line;
class Manifold;

class Operator;
class BasicOperator;
class ExponentialOperator;
class OperatorProduct;
class OperatorSum;

ExponentialOperator& exp(Operator& op, const int order);

enum {CLOSED = 0x1, OPEN = 0x2, CONNECTED = 0x4, DISCONNECTED = 0x8};

class Operator
{
    friend class BasicOperator;
    friend class ComplexOperator;
    friend class ExponentialOperator;
    friend class OperatorProduct;
    friend class OperatorSum;

    protected:
        int ref_count;

        Operator();

        void increment();

        void decrement();

    public:
        virtual ~Operator() {}

        OperatorProduct& operator*(Operator& other);

        OperatorSum& operator+(Operator& other);

        virtual Diagram resolve(Manifold& left, Manifold& right) const = 0;

        virtual ManifoldGenerator* matching(const Manifold& leftMin, const Manifold& leftMax,
                                            const Manifold& rightMin, const Manifold& rightMax) const = 0;

        static void canonicalize(Term& t);
};

BasicOperator& operator<<(BasicOperator& op, std::string term);

class BasicOperator : public Operator
{
    friend class BasicManifoldGenerator;

    friend BasicOperator& operator<<(BasicOperator& op, std::string term)
    {
        Term t(Diagram::SPINORBITAL, term);
        Operator::canonicalize(t);
        op.terms.push_back(t);
        return op;
    }

    protected:
        std::vector<Term> terms;

        void canonicalize();

    public:
        BasicOperator();

        BasicOperator(const Diagram& diagram);

        BasicOperator(const std::string& term);

        BasicOperator(const Term& term);

        BasicOperator(const std::vector<Term>& terms);

        virtual Diagram resolve(Manifold& left, Manifold& right) const;

        virtual ManifoldGenerator* matching(const Manifold& leftMin, const Manifold& leftMax,
                                            const Manifold& rightMin, const Manifold& rightMax) const;
};

class ComplexOperator : public Operator
{
    protected:
        Operator* operand;

    public:
        ComplexOperator(Operator& other);

        ComplexOperator& operator=(const ComplexOperator& other);

        virtual ~ComplexOperator();

        virtual Diagram resolve(Manifold& left, Manifold& right) const;

        virtual ManifoldGenerator* matching(const Manifold& leftMin, const Manifold& leftMax,
                                            const Manifold& rightMin, const Manifold& rightMax) const;
};

class ExponentialOperator : public Operator
{
    friend ExponentialOperator& exp(Operator& op, const int order)
    {
        return *(new ExponentialOperator(op, order));
    }

    friend class ExponentialManifoldGenerator;

    protected:
        Operator* expansion;

        ExponentialOperator(const ExponentialOperator& other);

        ExponentialOperator(Operator& op, int order);

        ExponentialOperator& operator=(const ExponentialOperator& other);

    public:
        virtual ~ExponentialOperator();

        virtual Diagram resolve(Manifold& left, Manifold& right) const;

        virtual ManifoldGenerator* matching(const Manifold& leftMin, const Manifold& leftMax,
                                            const Manifold& rightMin, const Manifold& rightMax) const;
};

class OperatorProduct : public Operator
{
    friend class Operator;
    friend class ProductManifoldGenerator;

    protected:
        Operator* l;
        Operator* r;
        int flags;

        OperatorProduct(const OperatorProduct& other);

        OperatorProduct(Operator& l, Operator& r, int flags = 0xF);

        OperatorProduct& operator=(const OperatorProduct& other);

        static bool isConnected(const Term& t);

        static bool isOpen(const Term& t);

    public:
        virtual ~OperatorProduct();

        OperatorProduct& operator()(int flags);

        void doProduct(Diagram& res, Term lTerm, Term rTerm, Manifold& c) const;

        virtual Diagram resolve(Manifold& left, Manifold& right) const;

        virtual ManifoldGenerator* matching(const Manifold& leftMin, const Manifold& leftMax,
                                            const Manifold& rightMin, const Manifold& rightMax) const;

        class ProductGenerator : public Generator
        {
            public:
                enum Side { LEFT, RIGHT };

            private:
                Manifold& c;
                const Term& term;
                const Side side;

                std::vector< std::vector<Line> > pindices;
                std::vector< std::vector<Line> > hindices;
                std::vector<int> which_bin_p;
                std::vector<int> which_bin_h;
                std::vector<int> how_many_p;
                std::vector<int> how_many_h;
                bool pdone, hdone, valid;
                int pi, hi, idx;

            public:
                ProductGenerator(Manifold& c, const Term& term, const Side side);

                bool next(Line which_p[], Line which_h[]);
        };
};

class OperatorSum : public Operator
{
    friend class Operator;
    friend class SumManifoldGenerator;

    protected:
        Operator* l;
        Operator* r;

        OperatorSum(const OperatorSum& other);

        OperatorSum(Operator& l, Operator& r);

        OperatorSum& operator=(const OperatorSum& other);

    public:
        virtual ~OperatorSum();

        virtual Diagram resolve(Manifold& left, Manifold& right) const;

        virtual ManifoldGenerator* matching(const Manifold& leftMin, const Manifold& leftMax,
                                            const Manifold& rightMin, const Manifold& rightMax) const;
};

}
}

#endif

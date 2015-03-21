#ifndef _AQUARIUS_AUTOCC_OPERATOR_HPP_
#define _AQUARIUS_AUTOCC_OPERATOR_HPP_

#include "util/global.hpp"

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

BasicOperator& operator<<(BasicOperator& op, string term);

class BasicOperator : public Operator
{
    friend class BasicManifoldGenerator;

    friend BasicOperator& operator<<(BasicOperator& op, string term)
    {
        Term t(Diagram::SPINORBITAL, term);
        Operator::canonicalize(t);
        op.terms.push_back(t);
        return op;
    }

    protected:
        vector<Term> terms;

        void canonicalize();

    public:
        BasicOperator();

        BasicOperator(const Diagram& diagram);

        BasicOperator(const string& term);

        BasicOperator(const Term& term);

        BasicOperator(const vector<Term>& terms);

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

                vector< vector<Line> > pindices;
                vector< vector<Line> > hindices;
                vector<int> which_bin_p;
                vector<int> which_bin_h;
                vector<int> how_many_p;
                vector<int> how_many_h;
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

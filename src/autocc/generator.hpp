#ifndef _AQUARIUS_AUTOCC_GENERATOR_HPP_
#define _AQUARIUS_AUTOCC_GENERATOR_HPP_

#include "util/global.hpp"

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
        vector<pair<Manifold,Manifold>> choices;
        vector<pair<Manifold,Manifold>>::iterator it;
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
        Manifold _zero, _one, _max, _next;
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
        pair<Manifold,Manifold> lNext, rNext;
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

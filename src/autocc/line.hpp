#ifndef _AQUARIUS_AUTOCC_LINE_HPP_
#define _AQUARIUS_AUTOCC_LINE_HPP_

#include "util/global.hpp"

namespace aquarius
{
namespace autocc
{

class Manifold;
class Line;

ostream& operator<<(ostream& out, const vector<Line>& v);
ostream& operator<<(ostream& out, const Line& l);

class Manifold
{
    public:
        map<int,int> np, nh;

        Manifold& operator+=(const Manifold& other);

        Manifold& operator-=(const Manifold& other);

        Manifold operator+(const Manifold& other) const;

        Manifold operator-(const Manifold& other) const;

        bool operator==(const Manifold& other) const;

        bool operator<(const Manifold& other) const;
};

class Line
{
    protected:
        uint16_t index;
        static const int TYPE_MAX = 26;
        static const int INDEX_MAX = 512;
        static const unsigned int SPIN_MASK  = 0x0200;
        static const unsigned int OCC_MASK   = 0x8000;
        static const unsigned int TYPE_MASK  = 0x7C00;
        static const unsigned int INDEX_MASK = 0x01FF;
        static const unsigned int _BETA      = 0x0200;
        static const unsigned int _OCCUPIED  = 0x8000;
        static const int SPIN_SHIFT =  9;
        static const int OCC_SHIFT  = 15;
        static const int TYPE_SHIFT = 10;

        Line(uint16_t index) : index(index) {}

    public:
        enum { ALPHA = 0, BETA = 1 };
        enum { OCCUPIED = 1, VIRTUAL = 0 };

        Line() : index(0) {}

        Line(const string& label);

        //Line(char label);

        Line(int index, int type, int occupancy, int spin);

        Line(int index, char type, int occupancy, int spin);

        static vector<Line> parse(const string& s);

        Line& operator=(const Line& line)
        {
            index = line.index;
            return *this;
        }

        bool operator<(const Line& other) const { return index < other.index; }

        bool operator==(const Line& other) const { return index == other.index; }

        bool operator!=(const Line& other) const { return index != other.index; }

        int getType() const { return (index&TYPE_MASK)>>TYPE_SHIFT; }

        int getOccupancy() const { return (index&OCC_MASK)>>OCC_SHIFT; }

        int getSpin() const { return (index&SPIN_MASK)>>SPIN_SHIFT; }

        int getIndex() const { return index&INDEX_MASK; }

        int asInt() const { return index; }

        bool isAlpha() const { return !(index&_BETA); }

        bool isBeta() const { return index&_BETA; }

        bool isVirtual() const { return !(index&_OCCUPIED); }

        bool isOccupied() const { return index&_OCCUPIED; }

        Line asIndex(int index) const
        {
            return Line(*this).toIndex(index);
        }

        Line& toIndex(int index)
        {
            assert(index >= 0 && index < INDEX_MAX);
            this->index = (this->index&(SPIN_MASK|OCC_MASK|TYPE_MASK))|index;
            return *this;
        }

        Line asType(int type) const
        {
            return Line(*this).toType(type);
        }

        Line& toType(int type)
        {
            assert(type >= 0 && type < TYPE_MAX);
            index = (index&(SPIN_MASK|OCC_MASK|INDEX_MASK))|(type<<TYPE_SHIFT);
            return *this;
        }

        Line asOccupancy(int occupancy) const
        {
            return Line(*this).toOccupancy(occupancy);
        }

        Line& toOccupancy(int occupancy)
        {
            assert(occupancy == OCCUPIED || occupancy == VIRTUAL);
            index = (index&(SPIN_MASK|TYPE_MAX|INDEX_MASK))|(occupancy<<OCC_SHIFT);
            return *this;
        }

        Line asSpin(int spin) const
        {
            return Line(*this).toSpin(spin);
        }

        Line& toSpin(int spin)
        {
            assert(spin == ALPHA || spin == BETA);
            index = (index&(OCC_MASK|TYPE_MASK|INDEX_MASK))|(spin<<SPIN_SHIFT);
            return *this;
        }

        Line asAlpha() const
        {
            return Line(index&(~_BETA));
        }

        Line& toAlpha()
        {
            index &= ~_BETA;
            return *this;
        }

        Line asBeta() const
        {
            return Line(index|_BETA);
        }

        Line& toBeta()
        {
            index |= _BETA;
            return *this;
        }

        Line asVirtual() const
        {
            return Line(index&(~_OCCUPIED));
        }

        Line& toVirtual()
        {
            index &= ~_OCCUPIED;
            return *this;
        }

        Line asOccupied() const
        {
            return Line(index|_OCCUPIED);
        }

        Line& toOccupied()
        {
            index |= _OCCUPIED;
            return *this;
        }
};

struct isAlpha : unary_function<Line,bool>
{
    bool operator()(const Line& line) const { return line.isAlpha(); }
};

struct isBeta : unary_function<Line,bool>
{
    bool operator()(const Line& line) const { return line.isBeta(); }
};

struct isOccupied : unary_function<Line,bool>
{
    bool operator()(const Line& line) const { return line.isOccupied(); }
};

struct isVirtual : unary_function<Line,bool>
{
    bool operator()(const Line& line) const { return line.isVirtual(); }
};

class isSpin : public unary_function<Line,bool>
{
    protected:
        int spin;

    public:
        isSpin(int spin) : spin(spin) {}

        bool operator()(const Line& line) const { return line.getSpin() == spin; }
};

class isOccupancy : public unary_function<Line,bool>
{
    protected:
        int occupancy;

    public:
        isOccupancy(int occupancy) : occupancy(occupancy) {}

        bool operator()(const Line& line) const { return line.getOccupancy() == occupancy; }
};

class isType : public unary_function<Line,bool>
{
    protected:
        int type;

    public:
        isType(int type) : type(type) {}

        bool operator()(const Line& line) const { return line.getType() == type; }
};

}
}

#endif

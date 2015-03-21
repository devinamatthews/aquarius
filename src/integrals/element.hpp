#ifndef _AQUARIUS_INTEGRALS_ELEMENT_HPP_
#define _AQUARIUS_INTEGRALS_ELEMENT_HPP_

#include "util/global.hpp"

#include "symmetry/symmetry.hpp"

namespace aquarius
{
namespace input
{
class Molecule;
}

namespace integrals
{

class Element
{
    protected:
        double charge;
        int Z;
        int A;
        int S;
        string symbol;
        string name;
        double mass;

        static const map<string,vector<Element>>& elements();

        Element(double charge, int Z, int A, int S, const string& symbol,
                const string& name, double mass)
        : charge(charge), Z(Z), A(A), S(S), symbol(symbol), name(name), mass(mass) {}

    public:
        static const Element& getElement(string symbol, int A);

        static const Element& getElement(string symbol);

        void setCharge(double charge) { this->charge = charge; }

        double getCharge() const { return charge; }

        int getAtomicNumber() const { return Z; }

        int getNucleonNumber() const { return A; }

        void setMass(double mass) { this->mass = mass; }

        double getMass() const { return mass; }

        int getSpin() const { return S; }

        const string& getName() const { return name; }

        const string& getSymbol() const { return symbol; }

        bool operator==(const Element& other) const
        {
            return Z == other.Z && A == other.A;
        }
};

}
}

#endif

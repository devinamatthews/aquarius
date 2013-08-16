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

#include "util/math.hpp"

#include "molecule.hpp"
#include "basis.hpp"

#include <map>
#include <string>
#include <stdexcept>
#include <cmath>

using namespace std;
using namespace aquarius;
using namespace aquarius::integrals;
using namespace aquarius::symmetry;
using namespace aquarius::input;
using namespace aquarius::task;

struct AtomSpec
{
    string symbol;
    string basisSet;
    string truncation;
    AtomSpec() {}
    AtomSpec(const string& symbol, const string& basisSet, const string& truncation)
    : symbol(symbol), basisSet(basisSet), truncation(truncation) {}
};

struct AtomZmatSpec : AtomSpec
{
    int distanceFrom, angleFrom, dihedralFrom;
    double distance, angle, dihedral;
    AtomZmatSpec() : distanceFrom(-1), angleFrom(-1), dihedralFrom(-1) {}
};

struct AtomCartSpec : AtomSpec
{
    vec3 pos;
    AtomCartSpec() {}
    AtomCartSpec(const string& symbol, const string& basisSet, const string& truncation, const vec3& pos)
    : AtomSpec(symbol, basisSet, truncation), pos(pos) {}
};

template<>
class Config::Extractor<AtomZmatSpec>
{
    private:
    static std::string nextSpec(node_t*& c)
    {
        for (;c && (c->data == "basis_set" || c->data == "truncation");c = c->next);

        if (c)
        {
            std::string s = c->data;
            c = c->next;
            return s;
        }

        return "";
    }

    public:
    inline static AtomZmatSpec extract(node_t* node, int which = 0)
    {
        AtomZmatSpec s;
        node_t *c;
        std::string str;

        if ((c = resolve(node,"basis_set"))) s.basisSet = c->data;
        if ((c = resolve(node,"truncation"))) s.truncation = c->data;

        c = node->children;

        str = nextSpec(c);
        if (str == "") throw BadValueError(path(node));
        s.symbol = str;

        str = nextSpec(c);
        if (str == "") return s;
        s.distanceFrom = Parser<int>::parse(str)-1;

        str = nextSpec(c);
        if (str == "") throw BadValueError(path(node));
        s.distance = Parser<double>::parse(str);

        str = nextSpec(c);
        if (str == "") return s;
        s.angleFrom = Parser<int>::parse(str)-1;

        str = nextSpec(c);
        if (str == "") throw BadValueError(path(node));
        s.angle = Parser<double>::parse(str);

        str = nextSpec(c);
        if (str == "") return s;
        s.dihedralFrom = Parser<int>::parse(str)-1;

        str = nextSpec(c);
        if (str == "") throw BadValueError(path(node));
        s.dihedral = Parser<double>::parse(str);

        return s;
    }
};

template<>
class Config::Extractor<AtomCartSpec>
{
    private:
    static std::string nextSpec(node_t*& c)
    {
        for (;c && (c->data == "basis_set" || c->data == "truncation");c = c->next);

        if (c)
        {
            std::string s = c->data;
            c = c->next;
            return s;
        }

        return "";
    }

    public:
    inline static AtomCartSpec extract(node_t* node, int which = 0)
    {
        AtomCartSpec s;
        node_t *c;
        std::string str;

        if ((c = resolve(node,"basis_set"))) s.basisSet = c->data;
        if ((c = resolve(node,"truncation"))) s.truncation = c->data;

        c = node->children;

        str = nextSpec(c);
        if (str == "") throw BadValueError(path(node));
        s.symbol = str;

        str = nextSpec(c);
        if (str == "") throw BadValueError(path(node));
        s.pos[0] = Parser<double>::parse(str);

        str = nextSpec(c);
        if (str == "") throw BadValueError(path(node));
        s.pos[1] = Parser<double>::parse(str);

        str = nextSpec(c);
        if (str == "") throw BadValueError(path(node));
        s.pos[2] = Parser<double>::parse(str);

        return s;
    }
};

Molecule::Molecule(const Arena& arena, const Config& config)
: Resource(arena)
{
    string name;
    bool hasDefaultBasis;
    bool contaminants = config.get<bool>("molecule.basis.contaminants");
    bool spherical = config.get<bool>("molecule.basis.spherical");
    bool angstrom = (config.get<string>("molecule.units") == "angstrom");
    bool zmat = (config.get<string>("molecule.coords") == "internal");
    BasisSet defaultBasis;

    try
    {
        name = config.get<string>("molecule.basis.basis_set");
        defaultBasis = BasisSet(TOPDIR "/basis/" + name);
        hasDefaultBasis = true;
    }
    catch (EntryNotFoundError& e)
    {
        hasDefaultBasis = false;
    }

    nelec = -config.get<int>("molecule.charge");
    multiplicity = config.get<int>("molecule.multiplicity");

    vector<AtomCartSpec> cartpos;

    if (zmat)
    {
        vector< pair<string,AtomZmatSpec> > atomspecs = config.find<AtomZmatSpec>("molecule.atom");
        for (vector< pair<string,AtomZmatSpec> >::iterator it = atomspecs.begin();it != atomspecs.end();++it)
        {
            vec3 pos, posb, posc, posd;
            AtomZmatSpec a = it->second;
            double bohr = (angstrom ? config.get<double>("constants.angstrom2bohr") : 1);
            double rad = 0.01745329251994329576923690768489;

            if (a.dihedralFrom != -1)
            {
                posb = cartpos[a.distanceFrom].pos;
                posc = cartpos[a.angleFrom].pos;
                posd = cartpos[a.dihedralFrom].pos;

                vec3 cb = unit(posb-posc);
                vec3 dcxcb = unit((posc-posd)^(posb-posc));
                vec3 dcperp = unit((posc-posd)/cb);

                pos = posb + bohr*a.distance*(sin(rad*a.angle)*(sin(rad*a.dihedral)*dcxcb -
                                                                cos(rad*a.dihedral)*dcperp) -
                                                                cos(rad*a.angle)*cb);
            }
            else if (a.angleFrom != -1)
            {
                posb = cartpos[a.distanceFrom].pos;
                pos[2] = posb[2] - bohr*a.distance*cos(rad*a.angle);
                pos[1] = bohr*a.distance*sin(rad*a.angle);
            }
            else if (a.distanceFrom != -1)
            {
                pos[2] = bohr*a.distance;
            }

            cartpos.push_back(AtomCartSpec(a.symbol, a.basisSet, a.truncation, pos));
        }
    }
    else
    {
        vector< pair<string,AtomCartSpec> > atomspecs = config.find<AtomCartSpec>("molecule.atom");
        for (vector< pair<string,AtomCartSpec> >::iterator it = atomspecs.begin();it != atomspecs.end();++it)
        {
            cartpos.push_back(it->second);
        }
    }

    vec3 avg;
    double totmass = 0;

    for (vector<AtomCartSpec>::iterator it = cartpos.begin();it != cartpos.end();++it)
    {
        Element e = Element::getElement(it->symbol.c_str());
        nelec += e.getAtomicNumber();
        avg += it->pos*e.getMass();
        totmass += e.getMass();
    }

    avg /= totmass;

    for (vector<AtomCartSpec>::iterator it = cartpos.begin();it != cartpos.end();++it)
    {
        it->pos -= avg;
    }

    //TODO: getSymmetry();

    PRINT("Molecular Geometry:\n\n");
    for (vector<AtomCartSpec>::iterator it = cartpos.begin();it != cartpos.end();++it)
    {
        PRINT("%3s % 20.15f % 20.15f % 20.15f\n", it->symbol.c_str(), it->pos[0], it->pos[1], it->pos[2]);
    }
    PRINT("\n");

    norb = 0;
    for (vector<AtomCartSpec>::iterator it = cartpos.begin();it != cartpos.end();++it)
    {
        Atom a(Center(PointGroup::C1(), it->pos, Element::getElement(it->symbol)));
        if (it->basisSet != "")
        {
            BasisSet(TOPDIR "/basis/" + it->basisSet).apply(a, spherical, contaminants);
        }
        else if (hasDefaultBasis)
        {
            defaultBasis.apply(a, spherical, contaminants);
        }

        for (vector<Shell>::iterator s = a.getShellsBegin();s != a.getShellsEnd();++s)
        {
            norb += s->getNFunc()*s->getNContr()*s->getDegeneracy();
        }

        atoms.push_back(a);
    }

    if ((nelec+multiplicity)%2 != 1 || multiplicity > nelec+1)
        throw logic_error("incompatible number of electrons and spin multiplicity");

    nucrep = 0;
    for (vector<Atom>::const_iterator a = atoms.begin();a != atoms.end();++a)
    {
        for (vector<Atom>::const_iterator b = atoms.begin();b != a;++b)
        {
            nucrep += a->getCenter().getElement().getCharge() *
                      b->getCenter().getElement().getCharge() /
                      dist(a->getCenter().getCenter(0), b->getCenter().getCenter(0));
        }
    }
}

int Molecule::getNumElectrons() const
{
    return nelec;
}

int Molecule::getNumAlphaElectrons() const
{
    return (nelec+multiplicity)/2;
}

int Molecule::getNumBetaElectrons() const
{
    return (nelec-multiplicity+1)/2;
}

int Molecule::getMultiplicity() const
{
    return multiplicity;
}

int Molecule::getNumOrbitals() const
{
    return norb;
}

double Molecule::getNuclearRepulsion() const
{
    return nucrep;
}

void Molecule::addAtom(const Atom& atom)
{
    atoms.push_back(atom);
}

Molecule::shell_iterator Molecule::getShellsBegin()
{
    if (!atoms.empty())
    {
        return shell_iterator(atoms.begin(), atoms.end(), atoms.front().getShellsBegin());
    }
    else
    {
        return shell_iterator(atoms.begin(), atoms.end(), vector<Shell>::iterator());
    }
}

Molecule::shell_iterator Molecule::getShellsEnd()
{
    if (!atoms.empty())
    {
        return shell_iterator(atoms.end(), atoms.end(), atoms.back().getShellsEnd());
    }
    else
    {
        return shell_iterator(atoms.begin(), atoms.end(), vector<Shell>::iterator());
    }
}

Molecule::const_shell_iterator Molecule::getShellsBegin() const
{
    if (!atoms.empty())
    {
        return const_shell_iterator(atoms.begin(), atoms.end(), atoms.front().getShellsBegin());
    }
    else
    {
        return const_shell_iterator(atoms.begin(), atoms.end(), vector<Shell>::const_iterator());
    }
}

Molecule::const_shell_iterator Molecule::getShellsEnd() const
{
    if (!atoms.empty())
    {
        return const_shell_iterator(atoms.end(), atoms.end(), atoms.back().getShellsEnd());
    }
    else
    {
        return const_shell_iterator(atoms.begin(), atoms.end(), vector<Shell>::const_iterator());
    }
}

vector<Atom>::iterator Molecule::getAtomsBegin()
{
    return atoms.begin();
}

vector<Atom>::iterator Molecule::getAtomsEnd()
{
    return atoms.end();
}

vector<Atom>::const_iterator Molecule::getAtomsBegin() const
{
    return atoms.begin();
}

vector<Atom>::const_iterator Molecule::getAtomsEnd() const
{
    return atoms.end();
}

Molecule::shell_iterator::shell_iterator(const vector<Atom>::iterator& atom_it,
               const vector<Atom>::iterator& atom_it_end,
               const vector<Shell>::iterator& shell_it)
: atom_it(atom_it), atom_it_end(atom_it_end), shell_it(shell_it)
{
}

Molecule::shell_iterator::shell_iterator()
: atom_it(), atom_it_end(), shell_it()
{
}

Molecule::shell_iterator::shell_iterator(const shell_iterator& other)
{
    this->atom_it = other.atom_it;
    this->atom_it_end = other.atom_it_end;
    this->shell_it = other.shell_it;
}

Molecule::shell_iterator& Molecule::shell_iterator::operator=(const shell_iterator& other)
{
    this->atom_it = other.atom_it;
    this->atom_it_end = other.atom_it_end;
    this->shell_it = other.shell_it;

    return *this;
}

Molecule::shell_iterator& Molecule::shell_iterator::operator++()
{
    if (atom_it != atom_it_end)
    {
        ++shell_it;

        if (shell_it == atom_it->getShellsEnd())
        {
            ++atom_it;
            if (atom_it != atom_it_end) shell_it = atom_it->getShellsBegin();
        }
    }

    return *this;
}

Molecule::shell_iterator Molecule::shell_iterator::operator++(int x)
{
    shell_iterator save = *this;

    ++(*this);

    return save;
}

Molecule::shell_iterator Molecule::shell_iterator::operator+(int x)
{
    shell_iterator r(*this);
    for (int i = 0;i < x;i++) ++r;
    return r;
}

Shell& Molecule::shell_iterator::operator*()
{
    return *shell_it;
}

Shell* Molecule::shell_iterator::operator->()
{
    return &(*shell_it);
}

bool Molecule::shell_iterator::operator<(const shell_iterator& other) const
{
    return atom_it < other.atom_it || (atom_it == other.atom_it && shell_it < other.shell_it);
}

bool Molecule::shell_iterator::operator==(const shell_iterator& other) const
{
    return atom_it == other.atom_it && shell_it == other.shell_it;
}

bool Molecule::shell_iterator::operator!=(const shell_iterator& other) const
{
    return atom_it != other.atom_it || shell_it != other.shell_it;
}

Molecule::const_shell_iterator::const_shell_iterator(const vector<Atom>::const_iterator& atom_it,
               const vector<Atom>::const_iterator& atom_it_end,
               const vector<Shell>::const_iterator& shell_it)
: atom_it(atom_it), atom_it_end(atom_it_end), shell_it(shell_it)
{
}

Molecule::const_shell_iterator::const_shell_iterator()
: atom_it(), atom_it_end(), shell_it()
{
}

Molecule::const_shell_iterator::const_shell_iterator(const const_shell_iterator& other)
{
    this->atom_it = other.atom_it;
    this->atom_it_end = other.atom_it_end;
    this->shell_it = other.shell_it;
}

Molecule::const_shell_iterator& Molecule::const_shell_iterator::operator=(const const_shell_iterator& other)
{
    this->atom_it = other.atom_it;
    this->atom_it_end = other.atom_it_end;
    this->shell_it = other.shell_it;

    return *this;
}

Molecule::const_shell_iterator& Molecule::const_shell_iterator::operator++()
{
    if (atom_it != atom_it_end)
    {
        ++shell_it;

        if (shell_it == atom_it->getShellsEnd())
        {
            ++atom_it;
            if (atom_it != atom_it_end) shell_it = atom_it->getShellsBegin();
        }
    }

    return *this;
}

Molecule::const_shell_iterator Molecule::const_shell_iterator::operator++(int x)
{
    const_shell_iterator save = *this;

    ++(*this);

    return save;
}

Molecule::const_shell_iterator Molecule::const_shell_iterator::operator+(int x)
{
    const_shell_iterator r(*this);
    for (int i = 0;i < x;i++) ++r;
    return r;
}

const Shell& Molecule::const_shell_iterator::operator*()
{
    return *shell_it;
}

const Shell* Molecule::const_shell_iterator::operator->()
{
    return &(*shell_it);
}

bool Molecule::const_shell_iterator::operator<(const const_shell_iterator& other) const
{
    return atom_it < other.atom_it || (atom_it == other.atom_it && shell_it < other.shell_it);
}

bool Molecule::const_shell_iterator::operator==(const const_shell_iterator& other) const
{
    return atom_it == other.atom_it && shell_it == other.shell_it;
}

bool Molecule::const_shell_iterator::operator!=(const const_shell_iterator& other) const
{
    return atom_it != other.atom_it || shell_it != other.shell_it;
}

Atom::Atom(const Center& center)
: center(center)
{
    shells = vector<Shell>();
}

void Atom::addShell(const Shell& shell)
{
    shells.push_back(shell);
}

Center& Atom::getCenter()
{
    return center;
}

const Center& Atom::getCenter() const
{
    return center;
}

vector<Shell>::iterator Atom::getShellsBegin()
{
    return shells.begin();
}

vector<Shell>::iterator Atom::getShellsEnd()
{
    return shells.end();
}

vector<Shell>::const_iterator Atom::getShellsBegin() const
{
    return shells.begin();
}

vector<Shell>::const_iterator Atom::getShellsEnd() const
{
    return shells.end();
}

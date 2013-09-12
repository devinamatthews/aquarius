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

#include "util/math_ext.h"

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

namespace aquarius
{
namespace input
{

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
    friend class Config::Extractor<AtomCartSpec>;

    protected:
        static std::string nextSpec(Config::node_t*& c)
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
        static AtomZmatSpec extract(node_t* node, int which = 0)
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
    public:
        static AtomCartSpec extract(node_t* node, int which = 0)
        {
            std::string (*nextSpec)(Config::node_t*& c) = Config::Extractor<AtomZmatSpec>::nextSpec;
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

}
}


MoleculeTask::MoleculeTask(const std::string& name, const input::Config& config)
: Task("molecule", name), config(config)
{
    addProduct(Product("molecule", "molecule"));
}

void MoleculeTask::run(task::TaskDAG& dag, const Arena& arena)
{
    put("molecule", new Molecule(arena, config));
}

REGISTER_TASK(MoleculeTask,"molecule");

Molecule::Molecule(const Arena& arena, const Config& config)
: Resource(arena)
{
    string name;
    bool hasDefaultBasis;
    bool contaminants = config.get<bool>("basis.contaminants");
    bool spherical = config.get<bool>("basis.spherical");
    bool angstrom = (config.get<string>("units") == "angstrom");
    bool zmat = (config.get<string>("coords") == "internal");
    BasisSet defaultBasis;

    try
    {
        name = config.get<string>("basis.basis_set");
        defaultBasis = BasisSet(TOPDIR "/basis/" + name);
        hasDefaultBasis = true;
    }
    catch (EntryNotFoundError& e)
    {
        hasDefaultBasis = false;
    }

    nelec = -config.get<int>("charge");
    multiplicity = config.get<int>("multiplicity");

    vector<AtomCartSpec> cartpos;
    double bohr = (angstrom ? config.get<double>("angstrom2bohr") : 1);

    if (zmat)
    {
        vector< pair<string,AtomZmatSpec> > atomspecs = config.find<AtomZmatSpec>("atom");
        for (vector< pair<string,AtomZmatSpec> >::iterator it = atomspecs.begin();it != atomspecs.end();++it)
        {
            vec3 pos, posb, posc, posd;
            AtomZmatSpec a = it->second;
            double rad = 0.01745329251994329576923690768489;

            if (a.dihedralFrom != -1)
            {
                posb = cartpos[a.distanceFrom].pos;
                posc = cartpos[a.angleFrom].pos;
                posd = cartpos[a.dihedralFrom].pos;

                vec3 cb = unit(posb-posc);
                vec3 dcxcb = unit((posc-posd)^(posb-posc));
                vec3 dcperp = unit((posc-posd)/cb);

                pos = posb + a.distance*(sin(rad*a.angle)*(sin(rad*a.dihedral)*dcxcb -
                                                           cos(rad*a.dihedral)*dcperp) -
                                                           cos(rad*a.angle)*cb);
            }
            else if (a.angleFrom != -1)
            {
                posb = cartpos[a.distanceFrom].pos;
                pos[2] = posb[2] + a.distance*cos(rad*a.angle);
                pos[1] = a.distance*sin(rad*a.angle);
            }
            else if (a.distanceFrom != -1)
            {
                pos[2] = a.distance;
            }

            cartpos.push_back(AtomCartSpec(a.symbol, a.basisSet, a.truncation, pos));
        }
    }
    else
    {
        vector< pair<string,AtomCartSpec> > atomspecs = config.find<AtomCartSpec>("atom");
        for (vector< pair<string,AtomCartSpec> >::iterator it = atomspecs.begin();it != atomspecs.end();++it)
        {
            cartpos.push_back(it->second);
        }
    }

    vec3 com;
    double totmass = 0;

    for (vector<AtomCartSpec>::iterator it = cartpos.begin();it != cartpos.end();++it)
    {
        Element e = Element::getElement(it->symbol.c_str());
        nelec += e.getAtomicNumber();
				it->pos *= bohr;
        com += it->pos*e.getMass();
        totmass += e.getMass();
    }

    com /= totmass;

    for (vector<AtomCartSpec>::iterator it = cartpos.begin();it != cartpos.end();++it)
    {
        it->pos -= com;
    }

    //TODO: getSymmetry();

    if (arena.rank == 0)
    {
        printf("\nMolecular Geometry:\n\n");
        for (vector<AtomCartSpec>::iterator it = cartpos.begin();it != cartpos.end();++it)
        {
						if (it->symbol != "X")
            printf("%3s % 20.15f % 20.15f % 20.15f\n", it->symbol.c_str(), it->pos[0], it->pos[1], it->pos[2]);
        }
    }

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

    if (arena.rank == 0)
    {
        printf("\nThere are %d atomic orbitals\n", norb);
        printf("There are %d alpha and %d beta electrons\n\n", getNumAlphaElectrons(), getNumBetaElectrons());
    }

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

Molecule::shell_iterator Molecule::getShellsBegin()
{
		vector<Atom>::iterator non_zero = atoms.begin();

		while (non_zero != atoms.end() && non_zero->getShellsBegin() == non_zero->getShellsEnd()) ++non_zero;

    if (non_zero != atoms.end())
    {
        return shell_iterator(non_zero, atoms.end(), non_zero->getShellsBegin());
    }
    else
    {
        return shell_iterator(atoms.begin(), atoms.end(), vector<Shell>::iterator());
    }
}

Molecule::shell_iterator Molecule::getShellsEnd()
{
		vector<Atom>::iterator non_zero = atoms.begin();

		while (non_zero != atoms.end() && non_zero->getShellsBegin() == non_zero->getShellsEnd()) ++non_zero;

    if (non_zero != atoms.end())
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
    return const_shell_iterator(const_cast<Molecule&>(*this).getShellsBegin());
}

Molecule::const_shell_iterator Molecule::getShellsEnd() const
{
    return const_shell_iterator(const_cast<Molecule&>(*this).getShellsEnd());
}

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
    nelec = -config.get<int>("charge");
    multiplicity = config.get<int>("multiplicity");

    vector<AtomCartSpec> cartpos;
    initGeometry(config, cartpos);
    initSymmetry(config, cartpos);
    initBasis(config, cartpos);

    if (arena.rank == 0)
    {
        printf("\nMolecular Geometry:\n\n");
        for (vector<Atom>::iterator it = atoms.begin();it != atoms.end();++it)
        {
            if (it->getCenter().getElement().getSymbol() != "X")
            {
                for (vector<vec3>::const_iterator pos = it->getCenter().getCenters().begin();
                     pos != it->getCenter().getCenters().end();++pos)
                    printf("%3s % 20.15f % 20.15f % 20.15f\n", it->getCenter().getElement().getSymbol().c_str(),
                            (*pos)[0], (*pos)[1], (*pos)[2]);
            }
        }
        cout << endl;
    }

    if ((nelec+multiplicity)%2 != 1 || multiplicity > nelec+1)
        throw logic_error("incompatible number of electrons and spin multiplicity");

    if (arena.rank == 0)
    {
        printf("Rotation constants (MHz): %15.6f %15.6f %15.6f\n", 29979.246*rota[0], 29979.246*rota[1], 29979.246*rota[2]);
        cout << "The molecular point group is " << group->getName() << endl;
        cout << "There are " << norb << " atomic orbitals by irrep\n";
        cout << "There are " << getNumAlphaElectrons() << " alpha and "
                             << getNumBetaElectrons() << " beta electrons\n\n";
    }
}

void Molecule::initGeometry(const Config& config, vector<AtomCartSpec>& cartpos)
{
    bool angstrom = (config.get<string>("units") == "angstrom");
    bool zmat = (config.get<string>("coords") == "internal");
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

    nucrep = 0;
    for (vector<AtomCartSpec>::iterator a = cartpos.begin();a != cartpos.end();++a)
    {
        for (vector<AtomCartSpec>::iterator b = cartpos.begin();b != a;++b)
        {
            Element ea = Element::getElement(a->symbol.c_str());
            Element eb = Element::getElement(b->symbol.c_str());
            nucrep += ea.getCharge()*eb.getCharge()/norm(a->pos-b->pos);
        }
    }
}

bool Molecule::isSymmetric(const vector<AtomCartSpec>& cartpos, const mat3x3& op)
{
    for (vector<AtomCartSpec>::const_iterator i1 = cartpos.begin();i1 != cartpos.end();++i1)
    {
        vec3 newpos = i1->pos*op;

        bool found = false;
        for (vector<AtomCartSpec>::const_iterator i2 = cartpos.begin();i2 != cartpos.end();++i2)
        {
            if (norm(newpos-i2->pos) < 1e-8) found = true;
        }
        if (!found) return false;
    }

    return true;
}

void Molecule::initSymmetry(const Config& config, vector<AtomCartSpec>& cartpos)
{
    mat3x3 I;

    for (vector<AtomCartSpec>::iterator it = cartpos.begin();it != cartpos.end();++it)
    {
        vec3& r = it->pos;
        Element e = Element::getElement(it->symbol.c_str());
        double m = e.getMass();
        I += m*(r*r-(r|r));
    }

    vec3 A;
    mat3x3 R;
    I.diagonalize(A, R);

    for (vector<AtomCartSpec>::iterator it = cartpos.begin();it != cartpos.end();++it)
    {
        it->pos = it->pos*R;
    }

    I = 0;
    for (vector<AtomCartSpec>::iterator it = cartpos.begin();it != cartpos.end();++it)
    {
        vec3& r = it->pos;
        Element e = Element::getElement(it->symbol.c_str());
        double m = e.getMass();
        I += m*(r*r-(r|r));
    }

    for (int i = 0;i < 3;i++)
    {
        for (int j = 0;j < 3;j++)
        {
            if (i != j) assert(I[i][j] < 1e-12);
        }
    }

    vec3 x(1,0,0);
    vec3 y(0,1,0);
    vec3 z(0,0,1);
    mat3x3 O = Identity();

    if (config.get<bool>("symmetry"))
    {
        if (A[1] < 1e-8)
        {
            /*
             * Atom: K (treat as D6h)
             */
            //group = &PointGroup::D6h();
            group = &PointGroup::D2h();
        }
        else if (A[0] < 1e-8)
        {
            /*
             * Linear molecules: CXv, DXh (treat as C6h and D6h)
             */
            if (isSymmetric(cartpos, Inversion()))
            {
                //group = &PointGroup::D6h();
                group = &PointGroup::D2h();
            }
            else
            {
                //group = &PointGroup::C6v();
                group = &PointGroup::C2v();
            }
            /*
             * Orient so that the molecule is along z
             */
            for (vector<AtomCartSpec>::iterator it = cartpos.begin();it != cartpos.end();++it)
            {
                if (it->pos.norm() > 1e-8)
                {
                    O = Rotation(it->pos, z);
                    break;
                }
            }
        }
        else if (2*abs(A[0]-A[2])/(A[0]+A[2]) < 1e-8)
        {
            /*
             * Spherical rotors: Td, Oh, Ih
             */
            assert(0);
        }
        else if (2*abs(A[0]-A[1])/(A[0]+A[1]) < 1e-8 ||
                 2*abs(A[1]-A[2])/(A[1]+A[2]) < 1e-8)
        {
            /*
             * Symmetric rotors: Cn, Cnv, Cnh, Dn, Dnh (all n>2), S2n, Dnd
             */
            assert(0);
        }
        else
        {
            /*
             * Asymmetric rotors: C1, Cs, Ci, C2, C2v, C2h, D2, D2h
             */
            bool c2x = isSymmetric(cartpos, C<2>(x));
            bool c2y = isSymmetric(cartpos, C<2>(y));
            bool c2z = isSymmetric(cartpos, C<2>(z));
            bool sx = isSymmetric(cartpos, Reflection(x));
            bool sy = isSymmetric(cartpos, Reflection(y));
            bool sz = isSymmetric(cartpos, Reflection(z));
            bool inv = isSymmetric(cartpos, Inversion());

            if (c2x && c2y && c2z)
            {
                /*
                 * D2, D2h
                 */
                if (inv)
                {
                    group = &PointGroup::D2h();
                }
                else
                {
                    group = &PointGroup::D2();
                }
                /*
                 * No reorientation needed
                 */
            }
            else if (c2x || c2y || c2z)
            {
                /*
                 * C2, C2v, C2h
                 */
                if (inv)
                {
                    group = &PointGroup::C2h();
                }
                else if (sx || sy || sz)
                {
                    group = &PointGroup::C2v();
                }
                else
                {
                    group = &PointGroup::C2();
                }
                /*
                 * Put C2 along z and make X > Y
                 */
                if (c2x)
                {
                    // x,y,z -> y,z,x
                    O = C<3>(vec3(1,1,1));
                }
                if (c2y)
                {
                    // x,y,z -> x,z,y
                    O = C<4>(x);
                }
            }
            else if (sx || sy || sz)
            {
                /*
                 * Cs
                 */
                group = &PointGroup::Cs();
                /*
                 * Put reflection plane orthogonal to z and make X > Y
                 */
                if (sx)
                {
                    // x,y,z -> y,z,x
                    O = C<3>(vec3(1,1,1));
                }
                if (sy)
                {
                    // x,y,z -> x,z,y
                    O = C<4>(x);
                }
            }
            else
            {
                /*
                 * C1, Ci
                 */
                if (inv)
                {
                    group = &PointGroup::Ci();
                }
                else
                {
                    group = &PointGroup::C1();
                }
                /*
                 * No reorientation needed
                 */
            }
        }
    }
    else
    {
        group = &PointGroup::C1();
    }

    for (vector<AtomCartSpec>::iterator it = cartpos.begin();it != cartpos.end();++it)
    {
        it->pos = O*it->pos;
    }

    I = 0;
    for (vector<AtomCartSpec>::iterator it = cartpos.begin();it != cartpos.end();++it)
    {
        vec3& r = it->pos;
        Element e = Element::getElement(it->symbol.c_str());
        double m = e.getMass();
        I += m*(r*r-(r|r));
    }

    for (int i = 0;i < 3;i++)
    {
        for (int j = 0;j < 3;j++)
        {
            if (i != j) assert(I[i][j] < 1e-12);
        }
    }

    rota[0] = 60.199687/I[0][0];
    rota[1] = 60.199687/I[1][1];
    rota[2] = 60.199687/I[2][2];

    for (vector<AtomCartSpec>::iterator i1 = cartpos.begin();;++i1)
    {
        if (i1 == cartpos.end()) break;

        vector<vec3> otherpos;
        for (int op = 0;op < group->getOrder();op++)
        {
            vec3 afterop = group->getOp(op)*i1->pos;
            if (norm(i1->pos-afterop) > 1e-8) otherpos.push_back(afterop);
        }

        for (vector<AtomCartSpec>::iterator i2 = i1;;)
        {
            if (i2 == cartpos.end()) break;

            bool match = false;
            for (vector<vec3>::iterator pos = otherpos.begin();pos != otherpos.end();++pos)
            {
                if (norm(*pos-i2->pos) < 1e-8) match = true;
            }

            if (match)
            {
                i2 = cartpos.erase(i2);
            }
            else
            {
                ++i2;
            }
        }
    }
}

void Molecule::initBasis(const Config& config, const vector<AtomCartSpec>& cartpos)
{
    bool contaminants = config.get<bool>("basis.contaminants");
    bool spherical = config.get<bool>("basis.spherical");

    BasisSet defaultBasis;
    bool hasDefaultBasis;
    try
    {
        string name = config.get<string>("basis.basis_set");
        defaultBasis = BasisSet(TOPDIR "/basis/" + name);
        hasDefaultBasis = true;
    }
    catch (EntryNotFoundError& e)
    {
        hasDefaultBasis = false;
    }

    norb.resize(group->getNumIrreps(), 0);

    for (vector<AtomCartSpec>::const_iterator it = cartpos.begin();it != cartpos.end();++it)
    {
        Atom a(Center(*group, it->pos, Element::getElement(it->symbol)));
        if (it->basisSet != "")
        {
            BasisSet(TOPDIR "/basis/" + it->basisSet).apply(a, spherical, contaminants);
        }
        else if (hasDefaultBasis)
        {
            defaultBasis.apply(a, spherical, contaminants);
        }

        atoms.push_back(a);

        for (vector<Shell>::iterator s = a.getShellsBegin();s != a.getShellsEnd();++s)
        {
            for (int i = 0;i < group->getNumIrreps();i++)
            {
                norb[i] += s->getNFuncInIrrep(i)*s->getNContr();
            }
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

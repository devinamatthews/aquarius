#include "molecule.hpp"

#include "util/math_ext.hpp"
#include "basis.hpp"

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
    double charge_from_input;
    AtomSpec() : charge_from_input(0.0) {}
    AtomSpec(const string& symbol, const string& basisSet,
             const string& truncation, double charge_from_input)
    : symbol(symbol), basisSet(basisSet), truncation(truncation),
      charge_from_input(charge_from_input) {}
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
    AtomCartSpec(const string& symbol, const string& basisSet, const string& truncation, const double& charge_from_input, const vec3& pos)
    : AtomSpec(symbol, basisSet, truncation, charge_from_input), pos(pos) {}
};

template<>
class Config::Extractor<AtomZmatSpec>
{
    protected:
        static bool nextSpec(shared_list<Node>::iterator& i, shared_list<Node>::iterator end)
        {
            for (;i != end && (i->data == "basis_set" || i->data == "truncation");++i);
            return (i != end);
        }

    public:
        static AtomZmatSpec extract(Node& node, int which = 0)
        {
            AtomZmatSpec s;

            Node* c;
            if ((c = resolve(node,"basis_set"))) s.basisSet = c->data;
            if ((c = resolve(node,"truncation"))) s.truncation = c->data;

            auto i = node.children.begin();
            auto e = node.children.end();

            if (!nextSpec(i, e)) throw BadValueError(node.fullName());
            s.symbol = i->data; ++i;

            if (!nextSpec(i, e)) return s;
            s.distanceFrom = Parser<int>::parse(i->data)-1; ++i;

            if (!nextSpec(i, e)) throw BadValueError(node.fullName());
            s.distance = Parser<double>::parse(i->data); ++i;

            if (!nextSpec(i, e)) return s;
            s.angleFrom = Parser<int>::parse(i->data)-1; ++i;

            if (!nextSpec(i, e)) throw BadValueError(node.fullName());
            s.angle = Parser<double>::parse(i->data); ++i;

            if (!nextSpec(i, e)) return s;
            s.dihedralFrom = Parser<int>::parse(i->data)-1; ++i;

            if (!nextSpec(i, e)) throw BadValueError(node.fullName());
            s.dihedral = Parser<double>::parse(i->data);

            return s;
        }
};

template<>
class Config::Extractor<AtomCartSpec>
{
    protected:
        static bool nextSpec(shared_list<Node>::iterator& i, shared_list<Node>::iterator end)
        {
            for (;i != end && (i->data == "basis_set" || i->data == "truncation");++i);
            return (i != end);
        }

    public:
        static AtomCartSpec extract(Node& node, int which = 0)
        {
            AtomCartSpec s;

            Node* c;
            if ((c = resolve(node,"basis_set"))) s.basisSet = c->data;
            if ((c = resolve(node,"truncation"))) s.truncation = c->data;

            auto i = node.children.begin();
            auto e = node.children.end();

            if (!nextSpec(i, e)) throw BadValueError(node.fullName());
            s.symbol = i->data; ++i;

            if (!nextSpec(i, e)) throw BadValueError(node.fullName());
            s.pos[0] = Parser<double>::parse(i->data); ++i;

            if (!nextSpec(i, e)) throw BadValueError(node.fullName());
            s.pos[1] = Parser<double>::parse(i->data); ++i;

            if (!nextSpec(i, e)) throw BadValueError(node.fullName());
            s.pos[2] = Parser<double>::parse(i->data); ++i;

            if (nextSpec(i, e))
                s.charge_from_input = Parser<double>::parse(i->data);

            return s;
        }
};

MoleculeTask::MoleculeTask(const string& name, input::Config& config)
: Task(name, config), config(config)
{
    addProduct(Product("molecule", "molecule"));
}

bool MoleculeTask::run(task::TaskDAG& dag, const Arena& arena)
{
    put("molecule", new Molecule(config, arena));
    return true;
}

Molecule::Molecule(Config& config, const Arena& arena)
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
            for (vector<vec3>::const_iterator pos = it->getCenter().getCenters().begin();
                 pos != it->getCenter().getCenters().end();++pos)
                printf("%3s % 20.15f % 20.15f % 20.15f\n", it->getCenter().getElement().getSymbol().c_str(),
                        (*pos)[0], (*pos)[1], (*pos)[2]);
        }
        cout << endl;

        if ((nelec+multiplicity)%2 != 1 || multiplicity > nelec+1)
            throw logic_error("incompatible number of electrons and spin multiplicity");

        printf("Rotation constants (MHz): %15.6f %15.6f %15.6f\n", 29979.246*rota[0], 29979.246*rota[1], 29979.246*rota[2]);
        cout << "The molecular point group is " << group->getName() << endl;
        cout << "There are " << norb << " atomic orbitals by irrep\n";
        cout << "There are " << getNumAlphaElectrons() << " alpha and "
                             << getNumBetaElectrons() << " beta electrons\n\n";
    }
}

void Molecule::initGeometry(Config& config, vector<AtomCartSpec>& cartpos)
{
    bool angstrom = (config.get<string>("units") == "angstrom");
    bool zmat = (config.get<string>("coords") == "internal");
    double bohr = (angstrom ? config.get<double>("angstrom2bohr") : 1);

    if (zmat)
    {
        for (auto& p : config.find<AtomZmatSpec>("atom"))
        {
            vec3 pos, posb, posc, posd;
            AtomZmatSpec a = p.second;
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
                posc = cartpos[a.angleFrom].pos;

                vec3 cb = unit(posb-posc);
                pos[2] = posb[2]-cb[2]*a.distance*cos(rad*a.angle);
                pos[1] = a.distance*sin(rad*a.angle);
            }
            else if (a.distanceFrom != -1)
            {
                pos[2] = a.distance;
            }

            cartpos.push_back(AtomCartSpec(a.symbol, a.basisSet, a.truncation, a.charge_from_input, pos));
        }
    }
    else
    {
        for (auto& p : config.find<AtomCartSpec>("atom"))
        {
            cartpos.push_back(p.second);
        }
    }

    vec3 com;
    double totmass = 0;

    for (auto it = cartpos.begin();it != cartpos.end();)
    {
        if (it->symbol == "X")
        {
            it = cartpos.erase(it);
        }
        else
        {
            Element e = Element::getElement(it->symbol);
            nelec += e.getAtomicNumber();
            it->pos *= bohr;
            com += it->pos*e.getMass();
            totmass += e.getMass();
            ++it;
        }
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
            if (a->charge_from_input != 0.0)
            {
                ea.setCharge(a->charge_from_input);
            }
            if (b->charge_from_input != 0.0)
            {
                eb.setCharge(b->charge_from_input);
            }
            nucrep += ea.getCharge()*eb.getCharge()/norm(a->pos-b->pos);
        }
    }
}

bool Molecule::isSymmetric(const vector<AtomCartSpec>& cartpos, const mat3x3& op, double tol)
{
    for (auto& a : cartpos)
    {
        vec3 newpos = op*a.pos;

        bool found = false;
        for (auto& b : cartpos)
        {
            if (a.symbol == b.symbol &&
                a.charge_from_input == b.charge_from_input &&
                norm(newpos-b.pos) < tol) found = true;
        }
        if (!found) return false;
    }

    return true;
}

void Molecule::rotate(vector<AtomCartSpec>& cartpos, const mat3x3& R)
{
    for (auto& a : cartpos)
    {
        a.pos = R*a.pos;
    }
}

mat3x3 Molecule::inertia(const vector<AtomCartSpec>& cartpos)
{
    mat3x3 I;

    for (auto& a : cartpos)
    {
        auto& r = a.pos;
        double m = Element::getElement(a.symbol).getMass();
        I += m*(r*r-(r|r));
    }

    return I;
}

void Molecule::initSymmetry(Config& config, vector<AtomCartSpec>& cartpos)
{
    double geomtol = config.get<double>("tolerance.geometry");
    double momtol = config.get<double>("tolerance.moments");

    vec3 A;
    mat3x3 R;
    mat3x3 I = inertia(cartpos);

    I.diagonalize(A, R);
    rotate(cartpos, R);
    I = inertia(cartpos);

    for (int i = 0;i < 3;i++)
    {
        for (int j = 0;j < 3;j++)
        {
            if (i != j) assert(I[i][j] < 1e-10);
        }
    }

    vec3 x(1,0,0);
    vec3 y(0,1,0);
    vec3 z(0,0,1);
    vec3 xyz(1,1,1);

    string subgrp = config.get<string>("subgroup");

    if (A[1] < momtol)
    {
        /*
         * Atom: K (treat as Ih)
         */
        //group = &PointGroup::Ih();
        group = &PointGroup::D2h();

        //if (subgrp == "full" || subgrp == "Ih") {}
        if (subgrp == "full" || subgrp == "D2h") {}
        else if (subgrp == "D2h") group = &PointGroup::D2h();
        else if (subgrp == "C2v") group = &PointGroup::C2v();
        else if (subgrp == "C2h") group = &PointGroup::C2h();
        else if (subgrp == "D2") group = &PointGroup::D2();
        else if (subgrp == "C2") group = &PointGroup::C2();
        else if (subgrp == "Ci") group = &PointGroup::Ci();
        else if (subgrp == "Cs") group = &PointGroup::Cs();
        else if (subgrp == "C1") group = &PointGroup::C1();
        else throw runtime_error(subgrp + " is not a valid subgroup of D2h");
    }
    else if (A[0] < momtol)
    {
        /*
         * Put molecule along z (initially it will be along x)
         */
        rotate(cartpos, C<4>(y));

        /*
         * Linear molecules: CXv, DXh (treat as C6v and D6h)
         */
        if (isSymmetric(cartpos, Inversion(), geomtol))
        {
            //group = &PointGroup::D6h();
            group = &PointGroup::D2h();

            //if (subgrp == "full" || subgrp == "D6h") {}
            if (subgrp == "full" || subgrp == "D2h") {}
            else if (subgrp == "C2v") group = &PointGroup::C2v();
            else if (subgrp == "C2h") group = &PointGroup::C2h();
            else if (subgrp == "D2") group = &PointGroup::D2();
            else if (subgrp == "C2") group = &PointGroup::C2();
            else if (subgrp == "Ci") group = &PointGroup::Ci();
            else if (subgrp == "Cs") group = &PointGroup::Cs();
            else if (subgrp == "C1") group = &PointGroup::C1();
            else throw runtime_error(subgrp + " is not a valid subgroup of D2h");
        }
        else
        {
            //group = &PointGroup::C6v();
            group = &PointGroup::C2v();

            //if (subgrp == "full" || subgrp == "C6v") {}
            if (subgrp == "full" || subgrp == "C2v") {}
            else if (subgrp == "C2") group = &PointGroup::C2();
            else if (subgrp == "Cs")
            {
                group = &PointGroup::Cs();
                rotate(cartpos, C<4>(x));
            }
            else if (subgrp == "C1") group = &PointGroup::C1();
            else throw runtime_error(subgrp + " is not a valid subgroup of C2v");
        }
    }
    else if (2*aquarius::abs(A[0]-A[2])/(A[0]+A[2]) < momtol)
    {
        /*
         * Spherical rotors: Td, Oh, Ih
         */

        /*
         * Search for a pair of C4 axes (Oh)
         */
        bool c4 = false;
        vec3 axis1, axis2;
        for (auto& a : cartpos)
        {
            if (norm(a.pos) < geomtol) continue;

            axis1 = unit(a.pos);
            if (isSymmetric(cartpos, C<4>(axis1), geomtol))
            {
                c4 = true;
                /*
                 * This C4 axis lies along an atom, search for a non-colinear
                 * atom which has another C4 axis.
                 */
                for (auto& b : cartpos)
                {
                    if (a.symbol != b.symbol ||
                        a.charge_from_input != b.charge_from_input ||
                        norm(b.pos) < geomtol ||
                        norm(unit(a.pos) + unit(b.pos)) < geomtol ||
                        norm(unit(a.pos) - unit(b.pos)) < geomtol) continue;

                    axis2 = unit(b.pos);
                    if (isSymmetric(cartpos, C<4>(axis2), geomtol)) break;
                }
                break;
            }

            for (auto& b : cartpos)
            {
                if (a.symbol != b.symbol ||
                    a.charge_from_input != b.charge_from_input ||
                    norm(b.pos) < geomtol ||
                    norm(unit(a.pos) + unit(b.pos)) < geomtol ||
                    norm(unit(a.pos) - unit(b.pos)) < geomtol) continue;

                axis1 = unit(a.pos + b.pos);
                if (isSymmetric(cartpos, C<4>(axis1), geomtol))
                {
                    c4 = true;
                    /*
                     * This C4 axis lies between two atoms. Another C4 axis is
                     * either 90 deg. in the direction of one of these atoms
                     * or 90 deg. but 45 deg. askew (i.e. between two of the
                     * C4-degenerate atoms).
                     */
                    axis2 = unit(b.pos%axis1);
                    if (isSymmetric(cartpos, C<4>(axis2), geomtol)) break;

                    axis2 = C<8>(axis1)*unit(b.pos%axis1);
                    if (isSymmetric(cartpos, C<4>(axis2), geomtol)) break;

                    assert(0);
                }
            }

            if (c4) break;
        }

        /*
         * Otherwise, search for a pair of S4 axes (Td)
         */
        bool s4 = false;
        if (!c4)
        {
            for (auto& a : cartpos)
            {
                if (norm(a.pos) < geomtol) continue;

                axis1 = unit(a.pos);
                if (isSymmetric(cartpos, S<4>(axis1), geomtol))
                {
                    s4 = true;
                    /*
                     * This S4 axis lies along an atom, search for a non-colinear
                     * atom which has another S4 axis.
                     */
                    for (auto& b : cartpos)
                    {
                        if (a.symbol != b.symbol ||
                            a.charge_from_input != b.charge_from_input ||
                            norm(b.pos) < geomtol ||
                            norm(unit(a.pos) + unit(b.pos)) < geomtol ||
                            norm(unit(a.pos) - unit(b.pos)) < geomtol) continue;

                        axis2 = unit(b.pos);
                        if (isSymmetric(cartpos, S<4>(axis2), geomtol)) break;
                    }
                    break;
                }

                for (auto& b : cartpos)
                {
                    if (a.symbol != b.symbol ||
                        a.charge_from_input != b.charge_from_input ||
                        norm(b.pos) < geomtol ||
                        norm(unit(a.pos) + unit(b.pos)) < geomtol ||
                        norm(unit(a.pos) - unit(b.pos)) < geomtol) continue;

                    axis1 = unit(a.pos + b.pos);
                    if (isSymmetric(cartpos, S<4>(axis1), geomtol))
                    {
                        s4 = true;
                        /*
                         * This S4 axis lies between two atoms. Another S4 axis is
                         * either 90 deg. in the direction of one of these atoms
                         * or 90 deg. but 45 deg. askew.
                         */
                        axis2 = unit(b.pos%axis1);
                        if (isSymmetric(cartpos, S<4>(axis2), geomtol)) break;

                        axis2 = C<8>(axis1)*unit(b.pos%axis1);
                        if (isSymmetric(cartpos, S<4>(axis2), geomtol)) break;

                        assert(0);
                    }
                }

                if (s4) break;
            }
        }

        /*
         * Otherwise it's Ih
         */

        if (c4)
        {
            rotate(cartpos, Rotation(axis1,z));
            axis2 = Rotation(axis1,z)*axis2;
            rotate(cartpos, Rotation(axis2,y));

            group = &PointGroup::Oh();

            if (subgrp == "full" || subgrp == "Oh") {}
            else if (subgrp == "Td") group = &PointGroup::Td();
            else if (subgrp == "S6")
            {
                group = &PointGroup::S6();
                rotate(cartpos, Rotation(vec3(1,1,1), z));
            }
            else if (subgrp == "S4") group = &PointGroup::S4();
            else if (subgrp == "D4h") group = &PointGroup::D4h();
            else if (subgrp == "D4") group = &PointGroup::D4();
            else if (subgrp == "C4h") group = &PointGroup::C4h();
            else if (subgrp == "C4v") group = &PointGroup::C4v();
            else if (subgrp == "C4") group = &PointGroup::C4();
            else if (subgrp == "D3d")
            {
                group = &PointGroup::D3d();
                rotate(cartpos, C<24>(-z)*Rotation(vec3(1,1,1), z));
            }
            else if (subgrp == "D3")
            {
                group = &PointGroup::D3();
                rotate(cartpos, C<24>(-z)*Rotation(vec3(1,1,1), z));
            }
            else if (subgrp == "C3v")
            {
                group = &PointGroup::C3v();
                rotate(cartpos, C<8>(-z)*Rotation(vec3(1,1,1), z));
            }
            else if (subgrp == "C3")
            {
                group = &PointGroup::C3();
                rotate(cartpos, Rotation(vec3(1,1,1), z));
            }
            else if (subgrp == "D2h") group = &PointGroup::D2h();
            else if (subgrp == "D2d") group = &PointGroup::D2d();
            else if (subgrp == "D2") group = &PointGroup::D2();
            else if (subgrp == "C2v") group = &PointGroup::C2v();
            else if (subgrp == "C2h") group = &PointGroup::C2h();
            else if (subgrp == "C2") group = &PointGroup::C2();
            else if (subgrp == "Ci") group = &PointGroup::Ci();
            else if (subgrp == "Cs") group = &PointGroup::Cs();
            else if (subgrp == "C1") group = &PointGroup::C1();
            else throw runtime_error(subgrp + " is not a valid subgroup of Oh");
        }
        else if (s4)
        {
            group = &PointGroup::Td();

            rotate(cartpos, Rotation(axis1,z));
            axis2 = Rotation(axis1,z)*axis2;
            rotate(cartpos, Rotation(axis2,y));

            if (subgrp == "full" || subgrp == "Td") {}
            else if (subgrp == "S4") group = &PointGroup::S4();
            else if (subgrp == "C3v")
            {
                group = &PointGroup::C3v();
                rotate(cartpos, C<8>(-z)*Rotation(vec3(1,1,1), z));
            }
            else if (subgrp == "C3")
            {
                group = &PointGroup::C3();
                rotate(cartpos, Rotation(vec3(1,1,1), z));
            }
            else if (subgrp == "D2d") group = &PointGroup::D2d();
            else if (subgrp == "D2") group = &PointGroup::D2();
            else if (subgrp == "C2v")
            {
                group = &PointGroup::C2v();
                rotate(cartpos, C<8>(z));
            }
            else if (subgrp == "C2") group = &PointGroup::C2();
            else if (subgrp == "Cs")
            {
                group = &PointGroup::Cs();
                rotate(cartpos, C<4>(x)*C<8>(z));
            }
            else if (subgrp == "C1") group = &PointGroup::C1();
            else throw runtime_error(subgrp + " is not a valid subgroup of Td");
        }
        else
        {
            /*
             * Look for a C2 axis first
             */
            bool c2 = false;
            vec3 axis1, axis2;
            for (auto& a : cartpos)
            {
                if (norm(a.pos) < geomtol) continue;

                for (auto& b : cartpos)
                {
                    if (a.symbol != b.symbol ||
                        a.charge_from_input != b.charge_from_input ||
                        norm(b.pos) < geomtol ||
                        norm(unit(a.pos) + unit(b.pos)) < geomtol ||
                        norm(unit(a.pos) - unit(b.pos)) < geomtol) continue;

                    axis1 = unit(a.pos+b.pos);
                    if (isSymmetric(cartpos, C<2>(axis1), geomtol))
                    {
                        c2 = true;
                        break;
                    }
                }
                if (c2) break;
            }
            assert(c2);

            /*
             * Then look for an orthogonal C2 axis
             */
            c2 = false;
            for (auto& a : cartpos)
            {
                if (norm(a.pos) < geomtol) continue;

                for (auto& b : cartpos)
                {
                    if (a.symbol != b.symbol ||
                        a.charge_from_input != b.charge_from_input ||
                        norm(b.pos) < geomtol ||
                        norm(unit(a.pos) + unit(b.pos)) < geomtol ||
                        norm(unit(a.pos) - unit(b.pos)) < geomtol) continue;

                    axis2 = unit(a.pos+b.pos);
                    if (axis1*axis2 < geomtol &&
                        isSymmetric(cartpos, C<2>(axis2), geomtol))
                    {
                        c2 = true;
                        break;
                    }
                }
                if (c2) break;
            }
            assert(c2);

            /*
             * Put the C2 axis along z and the C5 axes in the xy, xz, and yz planes
             */
            rotate(cartpos, Rotation(axis1,z));
            axis2 = Rotation(axis1,z)*(axis2%axis1);
            rotate(cartpos, Rotation(axis2,y));

            /*
             * Put the closest C5 axis to +z in the xz plane
             */
            axis1 = Rotation(x,31.71747441146100532424)*z;
            if (isSymmetric(cartpos, C<5>(axis1), geomtol))
                rotate(cartpos, C<4>(z));

            group = &PointGroup::Ih();

            if (subgrp == "full" || subgrp == "Ih") {}
            else if (subgrp == "S10")
            {
                group = &PointGroup::S10();
                rotate(cartpos, Rotation(y,31.71747441146100532424));
            }
            else if (subgrp == "D5d")
            {
                group = &PointGroup::D5d();
                rotate(cartpos, C<4>(z)*Rotation(y,31.71747441146100532424));
            }
            else if (subgrp == "D5")
            {
                group = &PointGroup::D5();
                rotate(cartpos, C<4>(z)*Rotation(y,31.71747441146100532424));
            }
            else if (subgrp == "C5v")
            {
                group = &PointGroup::C5v();
                rotate(cartpos, Rotation(y,31.71747441146100532424));
            }
            else if (subgrp == "C5")
            {
                group = &PointGroup::C5();
                rotate(cartpos, Rotation(y,31.71747441146100532424));
            }
            else if (subgrp == "S6")
            {
                group = &PointGroup::S6();
                rotate(cartpos, Rotation(vec3(1,1,1), z));
            }
            else if (subgrp == "D3d")
            {
                group = &PointGroup::D3d();
                rotate(cartpos, Rotation(x,20.9051574478892990329289550));
            }
            else if (subgrp == "D3")
            {
                group = &PointGroup::D3();
                rotate(cartpos, Rotation(x,20.9051574478892990329289550));
            }
            else if (subgrp == "C3v")
            {
                group = &PointGroup::C3v();
                rotate(cartpos, C<4>(z)*Rotation(x,20.9051574478892990329289550));
            }
            else if (subgrp == "C3")
            {
                group = &PointGroup::C3();
                rotate(cartpos, Rotation(x,20.9051574478892990329289550));
            }
            else if (subgrp == "D2h") group = &PointGroup::D2h();
            else if (subgrp == "D2") group = &PointGroup::D2();
            else if (subgrp == "C2h") group = &PointGroup::C2h();
            else if (subgrp == "C2v") group = &PointGroup::C2v();
            else if (subgrp == "C2") group = &PointGroup::C2();
            else if (subgrp == "Cs") group = &PointGroup::Cs();
            else if (subgrp == "Ci") group = &PointGroup::Ci();
            else if (subgrp == "C1") group = &PointGroup::C1();
            else throw runtime_error(subgrp + " is not a valid subgroup of Ih");
        }
    }
    else if (2*aquarius::abs(A[0]-A[1])/(A[0]+A[1]) < momtol ||
             2*aquarius::abs(A[1]-A[2])/(A[1]+A[2]) < momtol)
    {
        /*
         * Symmetric rotors: Cn, Cnv, Cnh, Dn, Dnh (all n>2), S2n, Dnd (both n>1)
         *
         * Place top axis along z
         */
        if (2*aquarius::abs(A[2]-A[1])/(A[2]+A[1]) < momtol)
            rotate(cartpos, Rotation(x, z));

        /*
         * Find the highest-order rotation axis (up to C6, C7 and higher must be treated
         * by a subgroup)
         */
        int order;
        for (order = 6;order > 1;order--)
        {
            if (isSymmetric(cartpos, Rotation(z, 360.0/order), geomtol)) break;
        }

        /*
         * Look for differentiators between Cnv, Cnh, Dn, etc.
         */
        bool inv = isSymmetric(cartpos, Inversion(), geomtol);
        bool sigmah = isSymmetric(cartpos, Reflection(z), geomtol);
        bool s2n = isSymmetric(cartpos, Rotation(z, 360.0/(2*order))*Reflection(z), geomtol);
        bool c2prime = false, sigmav = false, sigmad = false;
        vec3 axisv, axisd;

        /*
         * Find an atom not on the top axis, then look for C2',
         * sigma_v and/or sigma_d elements
         *
         * The axis for such an element (or the axis in the plane
         * and orthogonal to the top axis for sigma) must bisect
         * the line between two off-axis atoms
         */
        for (int atom1 = 0;atom1 < cartpos.size();atom1++)
        {
            if (norm(cartpos[atom1].pos%z) > geomtol)
            {
                for (int atom2 = 0;atom2 < cartpos.size();atom2++)
                {
                    axisv = unit(cartpos[atom1].pos+cartpos[atom2].pos);
                    if (norm(axisv%z) < geomtol) continue;
                    if ((c2prime = isSymmetric(cartpos, C<2>(axisv), geomtol))) break;
                }

                if (c2prime)
                {
                    sigmav = isSymmetric(cartpos, Reflection(axisv^z), geomtol);
                    axisd = Rotation(z, 360.0/(4*order))*axisv;
                    sigmad = isSymmetric(cartpos, Reflection(axisd^z), geomtol);
                }
                else
                {
                    for (int atom2 = 0;atom2 < cartpos.size();atom2++)
                    {
                        axisv = unit(cartpos[atom1].pos+cartpos[atom2].pos);
                        if (norm(axisv%z) < geomtol) continue;
                        if ((sigmav = isSymmetric(cartpos, Reflection(axisv^z), geomtol))) break;
                    }
                }
                break;
            }
        }

        /*
         * If there is no Cn, n<=6, then there must be Cn, n>6 and odd
         */
        if (order == 1)
        {
            /*
             * Dnd or S2n, n>6, treat as Ci
             */
            if (inv)
            {
                group = &PointGroup::Ci();

                if (subgrp == "full" || subgrp == "Ci") {}
                else if (subgrp == "C1") group = &PointGroup::C1();
                else throw runtime_error(subgrp + " is not a valid subgroup of Ci");
            }
            else if (c2prime)
            {
                /*
                 * Put the sigma_v axis on z and the top axis on x or y
                 * to preserve ordering of the rotatonal constants
                 */
                rotate(cartpos, Rotation(axisv, y));

                /*
                 * Dnh, n>6, treat as C2v
                 */
                if (sigmah)
                {
                    group = &PointGroup::C2v();

                    if (subgrp == "full" || subgrp == "C2v") {}
                    else if (subgrp == "C2") group = &PointGroup::C2();
                    else if (subgrp == "Cs")
                    {
                        group = &PointGroup::Cs();
                        rotate(cartpos, C<4>(x));
                    }
                    else if (subgrp == "C1") group = &PointGroup::C1();
                    else throw runtime_error(subgrp + " is not a valid subgroup of C2v");
                }
                /*
                 * Dn, n>6, treat as C2
                 */
                else
                {
                    group = &PointGroup::C2();

                    if (subgrp == "full" || subgrp == "C2") {}
                    else if (subgrp == "C1") group = &PointGroup::C1();
                    else throw runtime_error(subgrp + " is not a valid subgroup of C2");
                }
            }
            /*
             * Cnh, n>6, treat as Cs
             */
            else if (sigmah)
            {
                group = &PointGroup::Cs();

                if (subgrp == "full" || subgrp == "Cs") {}
                else if (subgrp == "C1") group = &PointGroup::C1();
                else throw runtime_error(subgrp + " is not a valid subgroup of Cs");
            }
            /*
             * Cnv, n>6, treat as Cs
             */
            else if (sigmav)
            {
                /*
                 * Put the sigma-axis on z and the top axis on x or y
                 * to preserve ordering of the rotatonal constants
                 */
                rotate(cartpos, C<4>(x)*Rotation(axisv, y));

                group = &PointGroup::Cs();

                if (subgrp == "full" || subgrp == "Cs") {}
                else if (subgrp == "C1") group = &PointGroup::C1();
                else throw runtime_error(subgrp + " is not a valid subgroup of Cs");
            }
            /*
             * Cn, n>6, no usable symmetry
             */
            else
            {
                group = &PointGroup::C1();

                if (subgrp == "full" || subgrp == "C1") {}
                else throw runtime_error(subgrp + " is not a valid subgroup of Cs");
            }
        }
        else
        {
            /*
             * Put the top axis on z and C2' or sigma_v axis
             * on y
             */
            if (c2prime || sigmav) rotate(cartpos, Rotation(axisv, y));

            if (s2n)
            {
                /*
                 * Dnd, treat D6d as D6
                 */
                if (c2prime)
                {
                    switch (order)
                    {
                        case 2: group = &PointGroup::D2d(); break;
                        case 3: group = &PointGroup::D3d(); break;
                        case 4: group = &PointGroup::D4d(); break;
                        case 5: group = &PointGroup::D5d(); break;
                        case 6: group = &PointGroup::D6(); break;
                    }

                    if (subgrp == "full" || subgrp == group->getName()) {}
                    else if (group == &PointGroup::D3d() && subgrp == "S6") group = &PointGroup::S6();
                    else if (group == &PointGroup::D2d() && subgrp == "S4") group = &PointGroup::S4();
                    else if (group == &PointGroup::D3d() && subgrp == "C3v")
                    {
                        /*
                         * Put the sigma_d axis on y
                         */
                        rotate(cartpos, C<12>(z));

                        group = &PointGroup::C2v();
                    }
                    else if (group == &PointGroup::D2d() && subgrp == "C2v")
                    {
                        /*
                         * Put the sigma_d axis on y
                         */
                        rotate(cartpos, C<8>(z));

                        group = &PointGroup::C2v();
                    }
                    else if ((group == &PointGroup::D3d() ||
                              group == &PointGroup::D6()) && subgrp == "D3") group = &PointGroup::D3();
                    else if ((group == &PointGroup::D2d() ||
                              group == &PointGroup::D4d() ||
                              group == &PointGroup::D6()) && subgrp == "D2") group = &PointGroup::D2();
                    else if ( group == &PointGroup::D5d() && subgrp == "C5") group = &PointGroup::C5();
                    else if ((group == &PointGroup::D3d() ||
                              group == &PointGroup::D6()) && subgrp == "C3") group = &PointGroup::C3();
                    else if (subgrp == "C2")
                    {
                        /*
                         * Put the sigma_v axis on z
                         */
                        rotate(cartpos, C<4>(x));

                        group = &PointGroup::C2();
                    }
                    else if ((group == &PointGroup::D2d() ||
                              group == &PointGroup::D3d()) && subgrp == "Cs")
                    {
                        /*
                         * Put the sigma_d axis on z
                         */
                        rotate(cartpos, C<4>(x)*Rotation(z, 360.0/(4*order)));

                        group = &PointGroup::Cs();
                    }
                    else if (group == &PointGroup::D3d() && subgrp == "Ci") group = &PointGroup::Ci();
                    else if (                               subgrp == "C1") group = &PointGroup::C1();
                    else throw runtime_error(subgrp + " is not a valid subgroup of " + group->getName());
                }
                /*
                 * S2n, treat S8 as C4, S10 as C5, and S12 as C6
                 */
                else
                {
                    switch (order)
                    {
                        case 2: group = &PointGroup::S4(); break;
                        case 3: group = &PointGroup::S6(); break;
                        case 4: group = &PointGroup::C4(); break;
                        case 5: group = &PointGroup::C5(); break;
                        case 6: group = &PointGroup::C6(); break;
                    }

                    if (subgrp == "full" || subgrp == group->getName()) {}
                    else if ((group == &PointGroup::S6() ||
                              group == &PointGroup::C6()) && subgrp == "C3") group = &PointGroup::C3();
                    else if ((group != &PointGroup::C5()  &&
                              group != &PointGroup::S6()) && subgrp == "C2") group = &PointGroup::C2();
                    else if ( group == &PointGroup::S6()  && subgrp == "Ci") group = &PointGroup::Ci();
                    else if (                                subgrp == "C1") group = &PointGroup::C1();
                    else throw runtime_error(subgrp + " is not a valid subgroup of " + group->getName());
                }
            }
            else if (c2prime)
            {
                /*
                 * Dnh
                 */
                if (sigmav)
                {
                    switch (order)
                    {
                        case 2: group = &PointGroup::D2h(); break;
                        case 3: group = &PointGroup::D3h(); break;
                        case 4: group = &PointGroup::D4h(); break;
                        case 5: group = &PointGroup::D5h(); break;
                        case 6: group = &PointGroup::D6h(); break;
                    }

                    if (subgrp == "full" || subgrp == group->getName()) {}
                    else if ( group == &PointGroup::D6h()  && subgrp == "C6h") group = &PointGroup::C6h();
                    else if ( group == &PointGroup::D5h()  && subgrp == "C5h") group = &PointGroup::C5h();
                    else if ( group == &PointGroup::D4h()  && subgrp == "C4h") group = &PointGroup::C4h();
                    else if ((group == &PointGroup::D6h() ||
                              group == &PointGroup::D3h()) && subgrp == "C3h") group = &PointGroup::C3h();
                    else if ((group == &PointGroup::D6h() ||
                              group == &PointGroup::D4h() ||
                              group == &PointGroup::D2h()) && subgrp == "C2h") group = &PointGroup::C2h();
                    else if ( group == &PointGroup::D6h()  && subgrp == "D3h") group = &PointGroup::D3h();
                    else if ((group == &PointGroup::D6h() ||
                              group == &PointGroup::D4h()) && subgrp == "D2h") group = &PointGroup::D2h();
                    else if ( group == &PointGroup::D6h()  && subgrp == "D3d") group = &PointGroup::D3d();
                    else if ( group == &PointGroup::D4h()  && subgrp == "D2d") group = &PointGroup::D2d();
                    else if ( group == &PointGroup::D6h()  && subgrp ==  "S6") group = &PointGroup::S6();
                    else if ( group == &PointGroup::D4h()  && subgrp ==  "S4") group = &PointGroup::S4();
                    else if ( group == &PointGroup::D6h()  && subgrp == "C6v") group = &PointGroup::C6v();
                    else if ( group == &PointGroup::D5h()  && subgrp == "C5v") group = &PointGroup::C5v();
                    else if ( group == &PointGroup::D4h()  && subgrp == "C4v") group = &PointGroup::C4v();
                    else if ((group == &PointGroup::D6h() ||
                              group == &PointGroup::D3h()) && subgrp == "C3v") group = &PointGroup::C3v();
                    else if (                                 subgrp == "C2v")
                    {
                        /*
                         * Put the sigma_v axis on z
                         */
                        rotate(cartpos, C<4>(x));

                        group = &PointGroup::C2v();
                    }
                    else if ( group == &PointGroup::D6h()  && subgrp ==  "D6") group = &PointGroup::D6();
                    else if ( group == &PointGroup::D5h()  && subgrp ==  "D5") group = &PointGroup::D5();
                    else if ( group == &PointGroup::D4h()  && subgrp ==  "D4") group = &PointGroup::D4();
                    else if ((group == &PointGroup::D6h() ||
                              group == &PointGroup::D3h()) && subgrp ==  "D3") group = &PointGroup::D3();
                    else if (                                 subgrp ==  "D2")
                    {
                        /*
                         * Put the C2' axis on z
                         */
                        rotate(cartpos, C<4>(x));

                        group = &PointGroup::D2();
                    }
                    else if ( group == &PointGroup::D6h()  && subgrp ==  "C6") group = &PointGroup::C6();
                    else if ( group == &PointGroup::D5h()  && subgrp ==  "C5") group = &PointGroup::C5();
                    else if ( group == &PointGroup::D4h()  && subgrp ==  "C4") group = &PointGroup::C4();
                    else if ((group == &PointGroup::D6h() ||
                              group == &PointGroup::D3h()) && subgrp ==  "C3") group = &PointGroup::C3();
                    else if (                                 subgrp ==  "C2")
                    {
                        /*
                         * Put the C2' axis on z
                         */
                        rotate(cartpos, C<4>(x));

                        group = &PointGroup::C2();
                    }
                    else if (                                 subgrp ==  "Cs") group = &PointGroup::Cs();
                    else if ((group == &PointGroup::D6h() ||
                              group == &PointGroup::D4h() ||
                              group == &PointGroup::D2h()) && subgrp ==  "Ci") group = &PointGroup::Ci();
                    else if (                                 subgrp ==  "C1") group = &PointGroup::C1();
                    else throw runtime_error(subgrp + " is not a valid subgroup of " + group->getName());
                }
                /*
                 * Dn
                 */
                else
                {
                    switch (order)
                    {
                        case 2: group = &PointGroup::D2(); break;
                        case 3: group = &PointGroup::D3(); break;
                        case 4: group = &PointGroup::D4(); break;
                        case 5: group = &PointGroup::D5(); break;
                        case 6: group = &PointGroup::D6(); break;
                    }

                    if (subgrp == "full" || subgrp == group->getName()) {}
                    else if ( group == &PointGroup::D6()  && subgrp ==  "D3") group = &PointGroup::D3();
                    else if ((group == &PointGroup::D6() ||
                              group == &PointGroup::D4()) && subgrp ==  "D2")
                    {
                        /*
                         * Put the C2' axis on z
                         */
                        rotate(cartpos, C<4>(x));

                        group = &PointGroup::D2();
                    }
                    else if ( group == &PointGroup::D6()  && subgrp ==  "C6") group = &PointGroup::C6();
                    else if ( group == &PointGroup::D5()  && subgrp ==  "C5") group = &PointGroup::C5();
                    else if ( group == &PointGroup::D4()  && subgrp ==  "C4") group = &PointGroup::C4();
                    else if ((group == &PointGroup::D6() ||
                              group == &PointGroup::D3()) && subgrp ==  "C3") group = &PointGroup::C3();
                    else if (                                subgrp ==  "C2")
                    {
                        /*
                         * Put the C2' axis on z
                         */
                        rotate(cartpos, C<4>(x));

                        group = &PointGroup::C2();
                    }
                    else if (                                subgrp ==  "C1") group = &PointGroup::C1();
                    else throw runtime_error(subgrp + " is not a valid subgroup of " + group->getName());
                }
            }
            /*
             * Cnh
             */
            else if (sigmah)
            {
                switch (order)
                {
                    case 2: group = &PointGroup::C2h(); break;
                    case 3: group = &PointGroup::C3h(); break;
                    case 4: group = &PointGroup::C4h(); break;
                    case 5: group = &PointGroup::C5h(); break;
                    case 6: group = &PointGroup::C6h(); break;
                }

                if (subgrp == "full" || subgrp == group->getName()) {}
                else if ( group == &PointGroup::C6h()  && subgrp == "C3h") group = &PointGroup::C3h();
                else if ((group == &PointGroup::C6h() ||
                          group == &PointGroup::C4h()) && subgrp == "C2h") group = &PointGroup::C2h();
                else if ( group == &PointGroup::C6h()  && subgrp ==  "S6") group = &PointGroup::S6();
                else if ( group == &PointGroup::C4h()  && subgrp ==  "S4") group = &PointGroup::S4();
                else if ( group == &PointGroup::C6h()  && subgrp ==  "C6") group = &PointGroup::C6();
                else if ( group == &PointGroup::C5h()  && subgrp ==  "C5") group = &PointGroup::C5();
                else if ( group == &PointGroup::C4h()  && subgrp ==  "C4") group = &PointGroup::C4();
                else if ((group == &PointGroup::C6h() ||
                          group == &PointGroup::C3h()) && subgrp ==  "C3") group = &PointGroup::C3();
                else if ((group == &PointGroup::C6h() ||
                          group == &PointGroup::C4h() ||
                          group == &PointGroup::C2h()) && subgrp ==  "C2") group = &PointGroup::C2();
                else if (                                 subgrp ==  "Cs") group = &PointGroup::Cs();
                else if ((group == &PointGroup::C6h() ||
                          group == &PointGroup::C4h() ||
                          group == &PointGroup::C2h()) && subgrp ==  "Ci") group = &PointGroup::Ci();
                else if (                                 subgrp ==  "C1") group = &PointGroup::C1();
                else throw runtime_error(subgrp + " is not a valid subgroup of " + group->getName());
            }
            /*
             * Cnv
             */
            else if (sigmav)
            {
                switch (order)
                {
                    case 2: group = &PointGroup::C2v(); break;
                    case 3: group = &PointGroup::C3v(); break;
                    case 4: group = &PointGroup::C4v(); break;
                    case 5: group = &PointGroup::C5v(); break;
                    case 6: group = &PointGroup::C6v(); break;
                }

                if (subgrp == "full" || subgrp == group->getName()) {}
                else if ( group == &PointGroup::C6v()  && subgrp == "C3v") group = &PointGroup::C3v();
                else if ((group == &PointGroup::C6v() ||
                          group == &PointGroup::C4v()) && subgrp == "C2v") group = &PointGroup::C2v();
                else if ( group == &PointGroup::C6v()  && subgrp ==  "C6") group = &PointGroup::C6();
                else if ( group == &PointGroup::C5v()  && subgrp ==  "C5") group = &PointGroup::C5();
                else if ( group == &PointGroup::C4v()  && subgrp ==  "C4") group = &PointGroup::C4();
                else if ((group == &PointGroup::C6v() ||
                          group == &PointGroup::C3v()) && subgrp ==  "C3") group = &PointGroup::C3();
                else if ((group == &PointGroup::C6v() ||
                          group == &PointGroup::C4v() ||
                          group == &PointGroup::C2v()) && subgrp ==  "C2") group = &PointGroup::C2();
                else if (                                 subgrp ==  "C1") group = &PointGroup::C1();
                else throw runtime_error(subgrp + " is not a valid subgroup of " + group->getName());
            }
            /*
             * Cn
             */
            else
            {
                switch (order)
                {
                    case 2: group = &PointGroup::C2(); break;
                    case 3: group = &PointGroup::C3(); break;
                    case 4: group = &PointGroup::C4(); break;
                    case 5: group = &PointGroup::C5(); break;
                    case 6: group = &PointGroup::C6(); break;
                }

                if (subgrp == "full" || subgrp == group->getName()) {}
                else if ( group == &PointGroup::C6()  && subgrp ==  "C6") group = &PointGroup::C6();
                else if ( group == &PointGroup::C5()  && subgrp ==  "C5") group = &PointGroup::C5();
                else if ( group == &PointGroup::C4()  && subgrp ==  "C4") group = &PointGroup::C4();
                else if ((group == &PointGroup::C6() ||
                          group == &PointGroup::C3()) && subgrp ==  "C3") group = &PointGroup::C3();
                else if ((group == &PointGroup::C6() ||
                          group == &PointGroup::C4() ||
                          group == &PointGroup::C2()) && subgrp ==  "C2") group = &PointGroup::C2();
                else if (                                subgrp ==  "C1") group = &PointGroup::C1();
                else throw runtime_error(subgrp + " is not a valid subgroup of " + group->getName());
            }
        }
    }
    else
    {
        /*
         * Asymmetric rotors: C1, Cs, Ci, C2, C2v, C2h, D2, D2h
         */
        bool c2x = isSymmetric(cartpos, C<2>(x), geomtol);
        bool c2y = isSymmetric(cartpos, C<2>(y), geomtol);
        bool c2z = isSymmetric(cartpos, C<2>(z), geomtol);
        bool sx = isSymmetric(cartpos, Reflection(x), geomtol);
        bool sy = isSymmetric(cartpos, Reflection(y), geomtol);
        bool sz = isSymmetric(cartpos, Reflection(z), geomtol);
        bool inv = isSymmetric(cartpos, Inversion(), geomtol);

        if (c2x && c2y && c2z)
        {
            /*
             * D2, D2h
             */
            if (inv)
            {
                group = &PointGroup::D2h();

                if (subgrp == "full" || subgrp == "D2h") {}
                else if (subgrp == "C2v") group = &PointGroup::C2v();
                else if (subgrp == "C2h") group = &PointGroup::C2h();
                else if (subgrp == "D2") group = &PointGroup::D2();
                else if (subgrp == "C2") group = &PointGroup::C2();
                else if (subgrp == "Ci") group = &PointGroup::Ci();
                else if (subgrp == "Cs") group = &PointGroup::Cs();
                else if (subgrp == "C1") group = &PointGroup::C1();
                else throw runtime_error(subgrp + " is not a valid subgroup of D2h");
            }
            else
            {
                group = &PointGroup::D2();

                if (subgrp == "full" || subgrp == "D2") {}
                else if (subgrp == "C2") group = &PointGroup::C2();
                else if (subgrp == "C1") group = &PointGroup::C1();
                else throw runtime_error(subgrp + " is not a valid subgroup of D2");
            }
        }
        else if (c2x || c2y || c2z)
        {
            /*
             * Put C2 along z and make X > Y
             */
            if (c2x)
            {
                // x,y,z -> y,z,x
                rotate(cartpos, C<3>(xyz));
            }
            if (c2y)
            {
                // x,y,z -> x,z,y
                rotate(cartpos, C<4>(x));
            }

            /*
             * C2, C2v, C2h
             */
            if (inv)
            {
                group = &PointGroup::C2h();

                if (subgrp == "full" || subgrp == "C2h") {}
                else if (subgrp == "C2") group = &PointGroup::C2();
                else if (subgrp == "Ci") group = &PointGroup::Ci();
                else if (subgrp == "Cs") group = &PointGroup::Cs();
                else if (subgrp == "C1") group = &PointGroup::C1();
                else throw runtime_error(subgrp + " is not a valid subgroup of C2h");
            }
            else if (sx || sy || sz)
            {
                group = &PointGroup::C2v();

                if (subgrp == "full" || subgrp == "C2v") {}
                else if (subgrp == "C2") group = &PointGroup::C2();
                else if (subgrp == "Cs")
                {
                    group = &PointGroup::Cs();
                    /*
                     * Put sigma along z and make X > Y
                     */
                    rotate(cartpos, C<4>(x));
                }
                else if (subgrp == "C1") group = &PointGroup::C1();
                else throw runtime_error(subgrp + " is not a valid subgroup of C2v");
            }
            else
            {
                group = &PointGroup::C2();

                if (subgrp == "full" || subgrp == "C2") {}
                else if (subgrp == "C1") group = &PointGroup::C1();
                else throw runtime_error(subgrp + " is not a valid subgroup of C2");
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
                rotate(cartpos, C<3>(xyz));
            }
            if (sy)
            {
                // x,y,z -> x,z,y
                rotate(cartpos, C<4>(x));
            }

            if (subgrp == "full" || subgrp == "Cs") {}
            else if (subgrp == "C1") group = &PointGroup::C1();
            else throw runtime_error(subgrp + " is not a valid subgroup of Cs");
        }
        else
        {
            /*
             * C1, Ci
             */
            if (inv)
            {
                group = &PointGroup::Ci();

                if (subgrp == "full" || subgrp == "Ci") {}
                else if (subgrp == "C1") group = &PointGroup::C1();
                else throw runtime_error(subgrp + " is not a valid subgroup of Ci");
            }
            else
            {
                group = &PointGroup::C1();

                if (subgrp == "full" || subgrp == "C1") {}
                else throw runtime_error(subgrp + " is not a valid subgroup of C1");
            }
        }
    }

    for (int i = 0;i < group->getOrder();i++)
    {
        if (!isSymmetric(cartpos, group->getOp(i), geomtol))
        {
            for (auto& a : cartpos) cout << a.symbol << " " << a.pos << endl;
            cout << group->getOpName(i) << endl;
            cout << group->getOp(i) << endl;
            throw runtime_error("molecule is not symmetric");
        }
    }

    I = inertia(cartpos);

    for (int i = 0;i < 3;i++)
    {
        for (int j = 0;j < 3;j++)
        {
            if (i != j) assert(I[i][j] < 1e-10);
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
            if (norm(i1->pos-afterop) > geomtol) otherpos.push_back(afterop);
        }

        for (vector<AtomCartSpec>::iterator i2 = i1;;)
        {
            if (i2 == cartpos.end()) break;

            bool match = false;
            for (vector<vec3>::iterator pos = otherpos.begin();pos != otherpos.end();++pos)
            {
                if (norm(*pos-i2->pos) < geomtol) match = true;
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

void Molecule::initBasis(Config& config, const vector<AtomCartSpec>& cartpos)
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
        Element myelem = Element::getElement(it->symbol);
        if (it->charge_from_input != 0.0)
        {
            myelem.setCharge(it->charge_from_input);
        }
        Atom a(Center(*group, it->pos, myelem));
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

class CheckSymmetryTask : public task::Task
{
    protected:
        string group;

    public:
        CheckSymmetryTask(const string& name, input::Config& config)
        : Task(name, config)
        {
            group = config.get<string>("group");

            vector<Requirement> reqs;
            reqs.push_back(Requirement("molecule", "molecule"));
            addProduct(Product("bool", "match", reqs));
        }

        bool run(TaskDAG& dag, const Arena& arena)
        {
            const Molecule& mol = get<Molecule>("molecule");

            bool match = mol.getGroup().getName() == group;

            if (match)
            {
                log(arena) << "passed" << endl;
            }
            else
            {
                error(arena) << "failed: " << mol.getGroup().getName() << " vs " << group << endl;
            }

            put("match", new bool(match));

            return true;
        }
};

}
}

static const char* spec = R"(

angstrom2bohr?
    double 1.88972612456506198632428439,
units?
    enum { angstrom, bohr },
coords?
    enum { internal, cartesian },
multiplicity?
    int 1,
charge?
    int 0,
subgroup?
    enum
    {
        C1, full,
        Cs, Ci, C2, C2v, C2h, D2, D2h,
        C3, C4, C5, C6, C3v, C4v, C5v, C6v,
        C3h, C4h, C5h, C6h, D3, D4, D5, D6,
        D3h, D4h, D5h, D6h, D2d, D3d, D4d, D5d,
        S4, S6, S8, S10, Td, Oh, Ih
    },
tolerance?
{
    geometry?
        double 1e-8,
    moments?
        double 1e-8,
},
atom+
{
    basis_set? string,
    truncation? string,
    # ZMAT specification, e.g. H 2 R 1 A
    # or xyz position
    *+
},
basis?
{
    contaminants?
        bool false,
    spherical?
        bool true,
    basis_set? string,
    truncation*
    {
        # elements affected, e.g. H-Ne
        *,
        # same format as molecule.atom.truncation
        *
    }
}

)";

REGISTER_TASK(aquarius::input::MoleculeTask, "molecule", spec);

REGISTER_TASK(aquarius::input::CheckSymmetryTask, "checksym", "group string");

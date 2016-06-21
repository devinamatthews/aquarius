#ifndef _AQUARIUS_FRAMEWORKS_MOLECULE_MOLECULE_HPP_
#define _AQUARIUS_FRAMEWORKS_MOLECULE_MOLECULE_HPP_

#include "frameworks/util.hpp"
#include "frameworks/task.hpp"
#include "frameworks/symmetry.hpp"

#include "atom.hpp"
#include "basis.hpp"
#include "center.hpp"
#include "element.hpp"
#include "shell.hpp"

namespace aquarius
{
namespace molecule
{

class Molecule
{
    protected:
        vector<Atom> atoms;
        int multiplicity;
        int nelec;
        vector<int> norb;
        double nucrep;
        const symmetry::PointGroup* group;
        double rota[3];

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

        friend class task::Config::Extractor<AtomZmatSpec>;
        friend class task::Config::Extractor<AtomCartSpec>;

        template <typename shell_type, typename atom_iterator_type, typename shell_iterator_type>
        class shell_iterator_ : public iterator<forward_iterator_tag, shell_type>
        {
            friend class Molecule;

            private:
                atom_iterator_type atom_it;
                atom_iterator_type atom_it_end;
                shell_iterator_type shell_it;

            protected:
                shell_iterator_(const atom_iterator_type& atom_it,
                                const atom_iterator_type& atom_it_end,
                                const shell_iterator_type& shell_it)
                : atom_it(atom_it), atom_it_end(atom_it_end), shell_it(shell_it) {}

            public:
                shell_iterator_() {}

                template <typename other_shell_type, typename other_atom_iterator_type,
                          typename other_shell_iterator_type>
                shell_iterator_(const shell_iterator_<other_shell_type,
                                                      other_atom_iterator_type,
                                                      other_shell_iterator_type>& other)
                : atom_it(other.atom_it), atom_it_end(other.atom_it_end), shell_it(other.shell_it) {}

                template <typename other_shell_type, typename other_atom_iterator_type,
                          typename other_shell_iterator_type>
                shell_iterator_ operator=(const shell_iterator_<other_shell_type,
                                                                other_atom_iterator_type,
                                                                other_shell_iterator_type>& other)
                {
                    atom_it = other.atom_it;
                    atom_it_end = other.atom_it_end;
                    shell_it = other.shell_it;
                    return *this;
                }

                shell_iterator_& operator++()
                {
                    if (atom_it != atom_it_end)
                    {
                        ++shell_it;

                        while (shell_it == atom_it->getShells().end())
                        {
                            ++atom_it;
                            if (atom_it != atom_it_end) shell_it = atom_it->getShells().begin();
														else break;
                        }
                    }

                    return *this;
                }

                shell_iterator_ operator++(int x)
                {
                    shell_iterator save = *this;
                    ++(*this);
                    return save;
                }

                shell_iterator_ operator+(int x)
                {
                    shell_iterator r(*this);
                    for (int i = 0;i < x;i++) ++r;
                    return r;
                }

                shell_type& operator*()
                {
                    return *shell_it;
                }

                shell_type* operator->()
                {
                    return &(*shell_it);
                }

                bool operator<(const shell_iterator_& other) const
                {
                    return atom_it < other.atom_it || (atom_it == other.atom_it && atom_it != atom_it_end && shell_it < other.shell_it);
                }

                bool operator==(const shell_iterator_& other) const
                {
                    return atom_it == other.atom_it && (atom_it == atom_it_end || shell_it == other.shell_it);
                }

                bool operator!=(const shell_iterator_& other) const
                {
                    return atom_it != other.atom_it || (atom_it != atom_it_end && shell_it != other.shell_it);
                }
        };

        static bool isSymmetric(const vector<AtomCartSpec>& cartpos, const mat3x3& op);

        vector<AtomCartSpec> initGeometry(const task::Config& config);

        void initSymmetry(const task::Config& config, vector<AtomCartSpec>& cartpos);

        void initBasis(const task::Config& config, const vector<AtomCartSpec>& cartpos);

    public:
        Molecule(task::Config& config, const Arena& arena);

        int getNumElectrons() const { return nelec; }

        int getNumAlphaElectrons() const { return (nelec+multiplicity)/2; }

        int getNumBetaElectrons() const { return (nelec-multiplicity+1)/2; }

        int getMultiplicity() const { return multiplicity; }

        const vector<int>& getNumOrbitals() const { return norb; }

        double getNuclearRepulsion() const { return nucrep; }

        const symmetry::PointGroup& getGroup() const { return *group; }

        typedef shell_iterator_<Shell,
                                vector<Atom>::iterator,
                                vector<Shell>::iterator > shell_iterator;
        typedef shell_iterator_<const Shell,
                                vector<Atom>::const_iterator,
                                vector<Shell>::const_iterator > const_shell_iterator;

        shell_iterator getShellsBegin();

        shell_iterator getShellsEnd();

        const_shell_iterator getShellsBegin() const;

        const_shell_iterator getShellsEnd() const;

        vector<Shell> getShells() const { return vector<Shell>(getShellsBegin(), getShellsEnd()); }

        vector<Atom>& getAtoms() { return atoms; }

        const vector<Atom>& getAtoms() const { return atoms; }
};

}
}

#endif

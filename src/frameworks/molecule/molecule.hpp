#ifndef _AQUARIUS_INPUT_MOLECULE_HPP_
#define _AQUARIUS_INPUT_MOLECULE_HPP_

#include "util/global.hpp"

#include "integrals/shell.hpp"
#include "symmetry/symmetry.hpp"
#include "task/task.hpp"

#include "config.hpp"

namespace aquarius
{
namespace input
{

class Atom;
class AtomCartSpec;

class Molecule
{
    friend class MoleculeTask;

    protected:
        vector<Atom> atoms;
        int multiplicity;
        int nelec;
        vector<int> norb;
        double nucrep;
        const symmetry::PointGroup *group;
        double rota[3];

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

                        while (shell_it == atom_it->getShellsEnd())
                        {
                            ++atom_it;
                            if (atom_it != atom_it_end) shell_it = atom_it->getShellsBegin();
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

        void initGeometry(input::Config& config, vector<AtomCartSpec>& cartpos);

        void initSymmetry(input::Config& config, vector<AtomCartSpec>& cartpos);

        void initBasis(input::Config& config, const vector<AtomCartSpec>& cartpos);

    public:
        Molecule(Config& config, const Arena& arena);

        void print(task::Printer& p) const {}

        int getNumElectrons() const { return nelec; }

        int getNumAlphaElectrons() const { return (nelec+multiplicity)/2; }

        int getNumBetaElectrons() const { return (nelec-multiplicity+1)/2; }

        int getMultiplicity() const { return multiplicity; }

        const vector<int>& getNumOrbitals() const { return norb; }

        double getNuclearRepulsion() const { return nucrep; }

        const symmetry::PointGroup& getGroup() const { return *group; }

        typedef shell_iterator_<integrals::Shell,
                                vector<Atom>::iterator,
                                vector<integrals::Shell>::iterator > shell_iterator;
        typedef shell_iterator_<const integrals::Shell,
                                vector<Atom>::const_iterator,
                                vector<integrals::Shell>::const_iterator > const_shell_iterator;

        shell_iterator getShellsBegin();

        shell_iterator getShellsEnd();

        const_shell_iterator getShellsBegin() const;

        const_shell_iterator getShellsEnd() const;

        vector<Atom>& getAtoms() { return atoms; }

        const vector<Atom>& getAtoms() const { return atoms; }

        vector<Atom>::iterator getAtomsBegin() { return atoms.begin(); }

        vector<Atom>::iterator getAtomsEnd() { return atoms.end(); }

        vector<Atom>::const_iterator getAtomsBegin() const { return atoms.begin(); }

        vector<Atom>::const_iterator getAtomsEnd() const { return atoms.end(); }
};

class Atom
{
    private:
        integrals::Center center;
        vector<integrals::Shell> shells;

    public:
        Atom(const integrals::Center& center) : center(center) {}

        void addShell(const integrals::Shell& shell) { shells.push_back(shell); }

        integrals::Center& getCenter() { return center; }

        const integrals::Center& getCenter() const { return center; }

        vector<integrals::Shell>::iterator getShellsBegin() { return shells.begin(); }

        vector<integrals::Shell>::iterator getShellsEnd() { return shells.end(); }

        vector<integrals::Shell>::const_iterator getShellsBegin() const { return shells.begin(); }

        vector<integrals::Shell>::const_iterator getShellsEnd() const { return shells.end(); }
};

class MoleculeTask : public task::Task
{
    protected:
        Config config;

    public:
        MoleculeTask(const string& name, input::Config& config);

        bool run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif

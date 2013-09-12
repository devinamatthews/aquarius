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

#ifndef _AQUARIUS_INPUT_MOLECULE_HPP_
#define _AQUARIUS_INPUT_MOLECULE_HPP_

#include <vector>
#include <iterator>

#include "integrals/shell.hpp"
#include "symmetry/symmetry.hpp"
#include "util/util.h"
#include "task/task.hpp"

#include "config.hpp"

namespace aquarius
{
namespace input
{

class Atom;

class Molecule : public task::Resource
{
    friend class MoleculeTask;

    protected:
        std::vector<Atom> atoms;
        int multiplicity;
        int nelec;
        int norb;
        double nucrep;

        template <typename shell_type, typename atom_iterator_type, typename shell_iterator_type>
        class shell_iterator_ : public std::iterator<std::forward_iterator_tag, shell_type>
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

    public:
        Molecule(const Arena& arena, const Config& config);

        void print(task::Printer& p) const {}

        int getNumElectrons() const { return nelec; }

        int getNumAlphaElectrons() const { return (nelec+multiplicity)/2; }

        int getNumBetaElectrons() const { return (nelec-multiplicity+1)/2; }

        int getMultiplicity() const { return multiplicity; }

        int getNumOrbitals() const { return norb; }

        double getNuclearRepulsion() const { return nucrep; }

        const symmetry::PointGroup& getPointGroup() const { return symmetry::PointGroup::C1(); }

        typedef shell_iterator_<integrals::Shell,
                                std::vector<Atom>::iterator,
                                std::vector<integrals::Shell>::iterator > shell_iterator;
        typedef shell_iterator_<const integrals::Shell,
                                std::vector<Atom>::const_iterator,
                                std::vector<integrals::Shell>::const_iterator > const_shell_iterator;

        shell_iterator getShellsBegin();

        shell_iterator getShellsEnd();

        const_shell_iterator getShellsBegin() const;

        const_shell_iterator getShellsEnd() const;

        std::vector<Atom>::iterator getAtomsBegin() { return atoms.begin(); }

        std::vector<Atom>::iterator getAtomsEnd() { return atoms.end(); }

        std::vector<Atom>::const_iterator getAtomsBegin() const { return atoms.begin(); }

        std::vector<Atom>::const_iterator getAtomsEnd() const { return atoms.end(); }
};

class Atom
{
    private:
        integrals::Center center;
        std::vector<integrals::Shell> shells;

    public:
        Atom(const integrals::Center& center) : center(center) {}

        void addShell(const integrals::Shell& shell) { shells.push_back(shell); }

        integrals::Center& getCenter() { return center; }

        const integrals::Center& getCenter() const { return center; }

        std::vector<integrals::Shell>::iterator getShellsBegin() { return shells.begin(); }

        std::vector<integrals::Shell>::iterator getShellsEnd() { return shells.end(); }

        std::vector<integrals::Shell>::const_iterator getShellsBegin() const { return shells.begin(); }

        std::vector<integrals::Shell>::const_iterator getShellsEnd() const { return shells.end(); }
};

class MoleculeTask : public task::Task
{
    protected:
        Config config;

    public:
        MoleculeTask(const std::string& name, const input::Config& config);

        void run(task::TaskDAG& dag, const Arena& arena);
};

}
}

#endif

#ifndef _AQUARIUS_FRAMEWORKS_MOLECULE_ATOM_HPP_
#define _AQUARIUS_FRAMEWORKS_MOLECULE_ATOM_HPP_

#include "center.hpp"
#include "shell.hpp"

namespace aquarius
{
namespace molecule
{

class Atom
{
    private:
        Center center;
        vector<Shell> shells;

    public:
        explicit Atom(const Center& center) : center(center) {}

        Center& getCenter() { return center; }

        const Center& getCenter() const { return center; }

        vector<Shell> getShells() { return shells; }

        const vector<Shell>& getShells() const { return shells; }
};

}
}

#endif

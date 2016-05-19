#ifndef _AQUARIUS_INTEGRALS_CENTER_HPP_
#define _AQUARIUS_INTEGRALS_CENTER_HPP_

#include "util/global.hpp"

#include "symmetry/symmetry.hpp"

#include "element.hpp"

namespace aquarius
{
namespace integrals
{

class Center
{
    protected:
        const symmetry::PointGroup* group;
        vector<int> stabilizer;
        vector<int> centermap;
        vector<vec3> centers;
        Element element;

    public:
        Center(const symmetry::PointGroup& group, const vec3& pos, const Element& element);

        const symmetry::PointGroup& getPointGroup() const { return *group; }

        const vector<int>& getStabilizer() const { return stabilizer; }

        const Element& getElement() const { return element; }

        const vector<vec3>& getCenters() const { return centers; }

        const vec3& getCenter(int degen) const { return centers[degen]; }

        int getCenterAfterOp(int op) const { return centermap[op]; }
};

}
}

#endif

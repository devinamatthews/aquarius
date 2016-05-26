#ifndef _AQUARIUS_FRAMEWORKS_MOLECULE_CENTER_HPP_
#define _AQUARIUS_FRAMEWORKS_MOLECULE_CENTER_HPP_

#include "frameworks/util.hpp"
#include "frameworks/symmetry.hpp"

#include "element.hpp"

namespace aquarius
{
namespace molecule
{

class Center
{
    protected:
        const symmetry::PointGroup& group;
        vector<int> stabilizer;
        vector<int> centermap;
        vector<vec3> centers;
        Element element;

    public:
        Center(const symmetry::PointGroup& group, const vec3& pos, const Element& element);

        const symmetry::PointGroup& getPointGroup() const { return group; }

        const vector<int>& getStabilizer() const { return stabilizer; }

        const Element& getElement() const { return element; }

        const vector<vec3>& getCenters() const { return centers; }

        const vec3& getCenter(int degen) const { return centers[degen]; }

        const vec3& getCenterAfterOp(int op) const { return centers[getIdxAfterOp(op)]; }

        int getIdxAfterOp(int op) const { return centermap[op]; }
};

}
}

#endif

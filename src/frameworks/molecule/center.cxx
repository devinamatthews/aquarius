#include "center.hpp"

using namespace aquarius::input;
using namespace aquarius::symmetry;

namespace aquarius
{
namespace integrals
{

Center::Center(const PointGroup& group, const vec3& pos, const Element& element)
: group(&group), centermap(group.getOrder(), -1), centers(1, pos), element(element)
{
    for (int g = 0;g < group.getOrder();g++)
    {
        vec3 image = group.getOp(g)*pos;

        for (int c = 0;c < centers.size();c++)
        {
            if (norm(centers[c]-image) < 1e-6) centermap[g] = c;
        }

        if (centermap[g] == -1)
        {
            centermap[g] = centers.size();
            centers.push_back(image);
        }

        if (centermap[g] == 0) stabilizer.push_back(g);
    }

    assert(group.getOrder() == centers.size()*stabilizer.size());
}

}
}

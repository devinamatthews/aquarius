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

#ifndef _AQUARIUS_INTEGRALS_CENTER_HPP_
#define _AQUARIUS_INTEGRALS_CENTER_HPP_

#include <cstddef>
#include <string>
#include <vector>
#include <stdexcept>
#include <cstring>

#include "memory/memory.h"
#include "util/math_ext.h"
#include "util/blas.h"
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
        std::vector<int> stabilizer;
        std::vector<int> centermap;
        std::vector<vec3> centers;
        Element element;

    public:
        Center(const symmetry::PointGroup& group, const vec3& pos, const Element& element);

        const symmetry::PointGroup& getPointGroup() const { return *group; }

        const std::vector<int>& getStabilizer() const { return stabilizer; }

        const Element& getElement() const { return element; }

        const std::vector<vec3>& getCenters() const { return centers; }

        const vec3& getCenter(int degen) const { return centers[degen]; }

        int getCenterAfterOp(int op) const { return centermap[op]; }
};

}
}

#endif

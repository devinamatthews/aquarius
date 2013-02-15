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

#include "slide.hpp"

#include "internal.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>

using namespace std;

namespace aquarius
{
namespace slide
{

SLIDE* SLIDE::slide = NULL;

SLIDE& SLIDE::init()
{
    if (slide == NULL)
    {
        slide = new SLIDE();
        if (SLIDE_init() != 0) throw SLIDEError("Could not initialize libslide");
    }

    return *slide;
}

void SLIDE::finish()
{
    if (slide != NULL)
    {
        delete slide;
        slide = NULL;
        if (SLIDE_finish() != 0) throw SLIDEError("Could not finish up libslide");
    }
}

void SLIDE::setCartesianOrdering(const int L, const int* ordering)
{
    if (SLIDE_set_cartesian_ordering(L, ordering) != 0) throw SLIDEError("Error setting cartesian ordering");
}

const int* SLIDE::getCartesianOrdering(const int L)
{
    return SLIDE_get_cartesian_ordering(L);
}

void SLIDE::setSphericalOrdering(const int L, const int* ordering)
{
    if (SLIDE_set_spherical_ordering(L, ordering) != 0) throw SLIDEError("Error setting spherical ordering");
}

const int* SLIDE::getSphericalOrdering(const int L)
{
    return SLIDE_get_spherical_ordering(L);
}

void SLIDE::setGroup(const int group)
{
    if (SLIDE_set_group(group) != 0) throw SLIDEError("Error setting group");
}

int SLIDE::getGroup()
{
    return SLIDE_get_group();
}

int SLIDE::getGroupOrder()
{
    return SLIDE_group_order();
}

void SLIDE::setOrdering(const int ordering)
{
    SLIDE_set_ordering(ordering);
}

int SLIDE::getOrdering()
{
    return SLIDE_get_ordering();
}

void SLIDE::setGeometryTolerance(const double tolerance)
{
    SLIDE_set_geometry_tolerance(tolerance);
}

double SLIDE::getGeometryTolerance()
{
    return SLIDE_get_geometry_tolerance();
}

const char* SLIDE::getGroupLabel(const int group)
{
    return SLIDE_group_label(group);
}

const char* SLIDE::getLLabel(const int L)
{
    return SLIDE_L_label(L);
}

const char* SLIDE::getOpLabel(const int op)
{
    return SLIDE_op_label(op);
}

const char* SLIDE::getIrrepLabel(const int irrep)
{
    return SLIDE_irrep_label(irrep);
}

const char* SLIDE::getFuncLabel(const int L, const bool spherical, const int func)
{
    return SLIDE_func_label(L, spherical, func);
}

Context::Context()
{
    context = SLIDE_new_context();

    if (context == NULL)
    {
        throw SLIDEError("Error creating context");
    }
}

Context::~Context()
{
    SLIDE_free_context(context);
}

Context::Context(const Context& other)
{
    context = SLIDE_new_context();

    if (context == NULL)
    {
        throw SLIDEError("Error creating context");
    }

    copy(other.context, other.context+1, context);

    if (context->aosize > 0)
    {
        context->aobuf1 = MALLOC(double, context->aosize);
        context->aobuf2 = MALLOC(double, context->aosize);
        if (context->aobuf1 == NULL || context->aobuf2 == NULL) throw SLIDEError("Could not allocate memory");
    }

    if (context->sosize > 0)
    {
        context->integrals = MALLOC(double, context->sosize);
        if (context->integrals == NULL) throw SLIDEError("Could not allocate memory");
        copy(other.context->integrals, other.context->integrals+other.context->num_integrals, context->integrals);
    }
}

Context& Context::operator=(Context other)
{
    swap(*this, other);
    return *this;
}

size_t Context::process2eInts(const size_t nprocess, double* integrals, idx4_t* indices, const double cutoff)
{
    return SLIDE_process_2e_ints(context, nprocess, integrals, indices, cutoff);
}

size_t Context::process1eInts(const size_t nprocess, double* integrals, idx2_t* indices, const double cutoff)
{
    return SLIDE_process_1e_ints(context, nprocess, integrals, indices, cutoff);
}

void Context::setNumProcessed(size_t num_processed)
{
    context->num_processed = num_processed;
}

size_t Context::getNumProcessed() const
{
    return context->num_processed;
}

size_t Context::getNumIntegrals() const
{
    return context->num_integrals;
}

double* Context::getIntegrals()
{
    return context->integrals;
}

const double* Context::getIntegrals() const
{
    return context->integrals;
}

void Context::calcERI(const double alpha, const double beta,
             const Shell& a, const Shell& b, const Shell& c, const Shell& d)
{
    if (SLIDE_calc_eri(context, alpha, beta, a.shell, b.shell, c.shell, d.shell) < 0)
    {
        throw SLIDEError("Error calculating ERIs");
    }
}

void Context::calcOVI(const double alpha, const double beta,
             const Shell& a, const Shell& b)
{
    if (SLIDE_calc_ovi(context, alpha, beta, a.shell, b.shell) < 0)
    {
        throw SLIDEError("Error calculating OVIs");
    }
}

void Context::calcKEI(const double alpha, const double beta,
             const Shell& a, const Shell& b)
{
    if (SLIDE_calc_kei(context, alpha, beta, a.shell, b.shell) < 0)
    {
        throw SLIDEError("Error calculating KEIs");
    }
}

void Context::calcNAI(const double alpha, const double beta,
             const Shell& a, const Shell& b, const Center centers[], const int ncenters)
{
    center_t** center_ts = new center_t*[ncenters];
    int i;

    for (i = 0;i < ncenters;i++) center_ts[i] = centers[i].center;

    if (SLIDE_calc_nai(context, alpha, beta, a.shell, b.shell, const_cast<const center_t**>(center_ts), ncenters) < 0)
    {
        delete[] center_ts;
        throw SLIDEError("Error calculating NAIs");
    }

    delete[] center_ts;
}

Element::~Element()
{
    SLIDE_free_element(element);
}

Element::Element(const Element& other)
{
    element = SLIDE_copy_element(other.element);
    if (element == NULL) throw SLIDEError("Error creating element");
}

Element& Element::operator=(Element other)
{
    swap(*this, other);
    return *this;
}

Element Element::getElement(const string& symbol)
{
    const element_t* e1 = SLIDE_get_element(symbol.c_str());
    if (e1 == NULL) throw SLIDEError("No such element " + symbol);
    element_t* e2 = SLIDE_copy_element(e1);
    if (e2 == NULL) throw SLIDEError("Error creating element");
    return Element(e2);
}

void Element::setCharge(const double charge)
{
    element->charge = charge;
}

double Element::getCharge() const
{
    return element->charge;
}

int Element::getAtomicNumber() const
{
    return element->Z;
}

int Element::getNucleonNumber() const
{
    return element->nucleon[element->isotope];
}

double Element::getMass() const
{
    return element->mass[element->isotope];
}

int Element::getSpin() const
{
    return element->spin[element->isotope];
}

const char* Element::getName() const
{
    return element->name;
}

const char* Element::getSymbol() const
{
    return element->symbol;
}

bool Element::operator==(const Element& other) const
{
    return (string(element->symbol) == string(other.element->symbol) && element->charge == other.element->charge);
}

Center::Center(const double pos[3], const Element& element)
: element(element)
{
    center = SLIDE_new_center(pos, element.element);
    if (center == NULL) throw SLIDEError("Error creating center");
}

Center::~Center()
{
    SLIDE_free_center(center);
}

Center::Center(const Center& other)
: element(other.element)
{
    center = SLIDE_new_center(other.center->centers[0], element.element);
    if (center == NULL) throw SLIDEError("Error creating center");
}

Center& Center::operator=(Center other)
{
    swap(*this, other);
    return *this;
}

double Center::nuclearRepulsion(const vector<Center>& centers)
{
    vector<Center>::const_iterator it;
    center_t** center_ts = new center_t*[centers.size()];
    int i;
    double nucrep;

    for (it = centers.begin();it != centers.end();++i) center_ts[i] = it->center;

    nucrep = SLIDE_nuclear_repulsion(const_cast<const center_t**>(center_ts), centers.size());

    delete[] center_ts;

    return nucrep;
}

const int* Center::getStabilizer() const
{
    return SLIDE_stabilizer(center);
}

int Center::getDegeneracy() const
{
    return center->degeneracy;
}

const Element& Center::getElement() const
{
    return element;
}

const double* Center::getCenter(const int degen) const
{
    return center->centers[degen];
}

const double* Center::getCenterAfterOp(const int op) const
{
    return center->centers[center->centermap[op]];
}

Shell::Shell(const Center& pos, const int L, const int nprim, const int ncontr, const bool spherical, const bool contaminants,
      const double* exponents, const double* coefficients, const int idx[8])
: center(pos)
{
    shell = SLIDE_new_shell(pos.center, L, nprim, ncontr, spherical, contaminants,
                            exponents, coefficients, idx);
    if (shell == NULL) throw SLIDEError("Error creating shell");
}

Shell::~Shell()
{
    SLIDE_free_shell(shell);
}

Shell::Shell(const Shell& other)
: center(other.center)
{
    shell = SLIDE_copy_shell(other.shell, other.shell->spherical, other.shell->contaminants);
    if (shell == NULL) throw SLIDEError("Error creating shell");
    shell->pos = center.center;
}

Shell::Shell(const Shell& other, const bool spherical, const bool contaminants)
: center(other.center)
{
    shell = SLIDE_copy_shell(other.shell, spherical, contaminants);
    if (shell == NULL) throw SLIDEError("Error creating shell");
    shell->pos = center.center;
}

Shell& Shell::operator=(Shell other)
{
    swap(*this, other);
    return *this;
}

void Shell::funcCenters(const int func, const int irrep, int* proj) const
{
    SLIDE_shell_func_centers(shell, func, irrep, proj);
}

void Shell::aoToSo(const int primitive_ordering, double* aoso, const int ld) const
{
    SLIDE_shell_ao_to_so(shell, primitive_ordering, aoso, ld);
}

int Shell::getL() const
{
    return shell->L;
}

int Shell::getNPrim() const
{
    return shell->nprim;
}

int Shell::getNContr() const
{
    return shell->ncontr;
}

int Shell::getNFunc() const
{
    return shell->nfunc;
}

int Shell::getNFuncInIrrep(const int irrep) const
{
    return shell->nfunc_irrep[irrep];
}

const int* Shell::getNFuncInEachIrrep() const
{
    return shell->nfunc_irrep;
}

const int* Shell::getIdx() const
{
    return shell->idx;
}

int* Shell::getIdx()
{
    return shell->idx;
}

bool Shell::getSpherical() const
{
    return shell->spherical;
}

bool Shell::getContaminants() const
{
    return shell->contaminants;
}

int Shell::getIrrepOfFunc(const int func, const int degen) const
{
    return shell->irreps[func][degen];
}

const int* Shell::getIrrepsOfFunc(const int func) const
{
    return shell->irreps[func];
}

const double* Shell::getExponents() const
{
    return shell->exponents;
}

const double* Shell::getCoefficients() const
{
    return shell->coefficients;
}

const Center& Shell::getCenter() const
{
    return center;
}

}
}

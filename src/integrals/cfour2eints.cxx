/* Copyright (c) 2014, Devin Matthews
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

#include "cfour2eints.hpp"
#include "internal.h"

#include <fstream>

using namespace std;
using namespace aquarius;
using namespace aquarius::integrals;
using namespace aquarius::input;
using namespace aquarius::symmetry;
using namespace aquarius::task;
using namespace aquarius::tensor;

CFOURTwoElectronIntegralsTask::CFOURTwoElectronIntegralsTask(const string& name, const Config& config)
: Task("cfour2eints", name)
{
    vector<Requirement> reqs;
    reqs.push_back(Requirement("molecule", "molecule"));
    addProduct(Product("eri", "I", reqs));
}

void CFOURTwoElectronIntegralsTask::run(TaskDAG& dag, const Arena& arena)
{
    const Molecule& molecule = get<Molecule>("molecule");

    ERI* eri = new ERI(arena, molecule.getGroup());

    const vector<int>& N = molecule.getNumOrbitals();
    int n = molecule.getGroup().getNumIrreps();

    ifstream ifs("IIII");

    int intsize = 8;
    int32_t junk;

    ifs.read((char*)&junk, 4);
    assert(junk == 224);

    ifs.read((char*)&junk, 4);
    if (junk != 0) intsize = 4;

    int batchsize = 600;
    vector<double> ints(batchsize);
    vector<idx4_t> idxs(batchsize);

    ifs.seekg(0);
    while (ifs)
    {
        int64_t recsize = 0;
        ifs.read((char*)&recsize, intsize);
        if (recsize != 8)
        {
            ifs.seekg(recsize+intsize, ifstream::cur);
            continue;
        }

        char label[9] = {};
        ifs.read(label, 8);
        ifs.seekg(intsize, ifstream::cur);
        SymmetryBlockedTensor<double> *tensor = NULL;
        if (strcmp(label, "TWOELSUP") != 0) continue;

        while (true)
        {
            int64_t numints;
            ifs.read((char*)&recsize, intsize);
            assert(recsize == 600*16+8);
            ifs.read((char*)ints.data(), 600*8);
            ifs.read((char*)idxs.data(), 600*8);
            ifs.read((char*)&numints, 8);
            ifs.seekg(intsize, ifstream::cur);

            if (numints < 0) break;

            eri->ints.insert(eri->ints.end(), ints.begin(), ints.begin()+numints);
            eri->idxs.insert(eri->idxs.end(), idxs.begin(), idxs.begin()+numints);
        }
    }

    for (size_t i = 0;i < eri->idxs.size();i++)
    {
        eri->idxs[i].i--;
        eri->idxs[i].j--;
        eri->idxs[i].k--;
        eri->idxs[i].l--;
    }

    put("I", eri);
}

REGISTER_TASK(CFOURTwoElectronIntegralsTask,"cfour2eints");

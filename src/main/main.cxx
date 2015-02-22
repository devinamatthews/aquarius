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

#include <cstdlib>
#include "tensor/symblocked_tensor.hpp"

#include "time/time.hpp"
#include "task/task.hpp"

#include <exception>
#include "mpi.h"
#include "omp.h"

using namespace std;
using namespace aquarius;
using namespace aquarius::time;
using namespace aquarius::task;

int main(int argc, char **argv)
{
    #ifdef ELEMENTAL
    El::Initialize(argc, argv);
    #else
    MPI::Init(argc, argv);
    #endif

    if (getenv("OMP_NUM_THREADS") == NULL)
    {
        omp_set_num_threads(1);
    }

    if (MPI::COMM_WORLD.Get_rank() == 0)
    {
        printf("================================================================================\n");
        printf("                                                                                \n");
        printf("                                                          ii                    \n");
        printf("            aaaaa        qqq    u   u    aaaa    rr rrr        u   u   sssss    \n");
        printf("        aaaaaaaaaaa    qq  qq  uu  uu   aa  aa    rr  rr  ii  uu  uu  ss  ss    \n");
        printf("      aaaa       aaa  qq   qq  uu  uu  aa  aaa   rr       ii  uu  uu  sss       \n");
        printf("    aaa          aaa   qqqqq    uuuu    aaaa aa  r         ii  uuuu     sss     \n");
        printf("   aaa           aaa      qq  q                                           ss    \n");
        printf("   aaa           aaa       qqq                                       sss  ss    \n");
        printf("    aaa        aaaaa                                                   sss      \n");
        printf("     aaaaaaaaaaaa aaa                                                           \n");
        printf("       aaaaaaa     aaa         Advanced Quantum Molecular Iterative             \n");
        printf("                                    Equation Solver                             \n");
        printf("                                                                                \n");
        printf("================================================================================\n");
        printf("\n");
        printf("Running on %d process%s with %d thread%s%s\n\n",
               MPI::COMM_WORLD.Get_size(),
               (MPI::COMM_WORLD.Get_size() > 1 ? "es" : ""),
               omp_get_max_threads(),
               (omp_get_max_threads() > 1 ? "s" : ""),
               (MPI::COMM_WORLD.Get_size() > 1 ? " each" : ""));
    }

    int status = 0;

    do
    {
        Arena world;

        if (argc < 2)
        {
            Logger::error(world) << "No input file specified." << endl << endl;
        }
        else
        {
            try
            {
                TaskDAG dag(argv[1]);
                dag.execute(world);
            }
            catch (const runtime_error& e)
            {
                Logger::error(world) << e.what() << endl << endl;
                status = 1;
            }
        }

        Timer::printTimers(world);
    }
    while (false);

    #ifdef ELEMENTAL
    El::Finalize();
    #else
    MPI::Finalize();
    #endif

    return status;
}

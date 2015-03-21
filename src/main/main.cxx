#include "tensor/symblocked_tensor.hpp"
#include "time/time.hpp"
#include "task/task.hpp"
#include "util/distributed.hpp"

using namespace aquarius;
using namespace aquarius::time;
using namespace aquarius::task;

int main(int argc, char **argv)
{
    #ifdef ELEMENTAL
    El::Initialize(argc, argv);
    #else
    MPI_Init(&argc, &argv);
    #endif

    if (getenv("OMP_NUM_THREADS") == NULL)
    {
        omp_set_num_threads(1);
    }

    int status = 0;

    {
        Arena world;

        if (world.rank == 0)
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
            printf("       aaaaaaa     aaa                Welcome to the New Age                    \n");
            printf("                                                                                \n");
            printf("                                                                                \n");
            printf("================================================================================\n");
            printf("\n");
            printf("Running on %d process%s with %d thread%s%s\n\n",
                   world.size,
                   (world.size > 1 ? "es" : ""),
                   omp_get_max_threads(),
                   (omp_get_max_threads() > 1 ? "s" : ""),
                   (world.size > 1 ? " each" : ""));
        }

        if (argc < 2)
        {
            Logger::error(world) << "No input file specified." << endl << endl;
        }
        else
        {
            //try
            //{
                TaskDAG dag(argv[1]);
                dag.execute(world);
            //}
            //catch (const runtime_error& e)
            //{
            //    Logger::error(world) << e.what() << endl << endl;
            //    status = 1;
            //}
        }

        Timer::printTimers(world);
    }

    #ifdef ELEMENTAL
    El::Finalize();
    #else
    MPI_Finalize();
    #endif

    return status;
}

#include "util/global.hpp"

#include "tensor/symblocked_tensor.hpp"
#include "time/time.hpp"
#include "task/task.hpp"

#ifdef HAVE_LIBINT2
#include "libint2.h"
#endif

using namespace aquarius;
using namespace aquarius::time;
using namespace aquarius::task;

int main(int argc, char **argv)
{
    #ifdef HAVE_ELEMENTAL
    El::Initialize(argc, argv);
    #else
    MPI_Init(&argc, &argv);
    #endif

    #ifdef HAVE_LIBINT2
    libint2_static_init();
    printf("sdflkjsdf\n");
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
            srand(::time(NULL));
            switch (rand()%5)
            {
                case 0:
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
                break;
                case 1:
                    printf("================================================================================\n");
                    printf("                                                                                \n");
                    printf("                 \\ \\ /\\ / /                                                     \n");
                    printf("               ---- /__\\ ----                                                   \n");
                    printf("             ----  /|_O|\\ ----                              ||                  \n");
                    printf("              --- /______\\ ---                                                  \n");
                    printf("               / __________ \\    /==\\   || ||  /==\\   \\/=\\  ||  || ||  /==\\     \n");
                    printf("                /_|____|___\\    ||  ||  || || ||  ||  ||    ||  || || |\\        \n");
                    printf("               /____|____|__\\    \\==||  \\|=|/  \\==/\\\\ ||    ||  \\|=|/  \\=\\      \n");
                    printf("              /__|____|____|_\\      //                                   \\|     \n");
                    printf("             /|____|____|____|\\     \\=/                                \\==/     \n");
                    printf("            /___|____|____|____\\                                                \n");
                    printf("           /___|____|____|____|_\\                                               \n");
                    printf("          /_|____|____|____|____|\\                                              \n");
                    printf("         /____|____|____|____|____\\          Novus Ordo Seclorum                \n");
                    printf("        /_|____|____|____|____|____\\                                            \n");
                    printf("       /|___|____|____|____|____|___\\                                           \n");
                    printf("      /__|____|____|____|____|____|__\\                                          \n");
                    printf("     /__|____|____|____|____|____|____\\                                         \n");
                    printf("    /|____|____|____|____|____|____|___\\                                        \n");
                    printf("                                                                                \n");
                    printf("================================================================================\n");
                break;
                case 2:
                    printf("================================================================================\n");
                    printf("                                                                                \n");
                    printf("                                     __                                         \n");
                    printf("                                    /  \\                                        \n");
                    printf("                 ____   ____ ____   \\__/    ______    _____ ____ ____  ____     \n");
                    printf("        /\\      / __ \\   | | | |     /\\      |  _ \\    | |   | | | |  / ___\\|   \n");
                    printf("       /  \\    / /  \\ \\  | | | |    /  \\     | |/ /    | |   | | | |  | |__     \n");
                    printf("      / /\\ \\   | |  | |  | | | |   / /\\ \\    |   \\     | |   | | | |  \\__  \\    \n");
                    printf("     / ____ \\  \\ \\__/ /  | \\_/ /  / /  \\ \\   | |\\ \\    | |   \\ \\_/ /  ___| |    \n");
                    printf("   _/_/_  _\\_\\_ \\__  /    \\___/  /__\\  /__\\ _|_| \\_\\_ _|_|_   \\___/  |\\____/    \n");
                    printf("                   \\_\\                                                          \n");
                    printf("                                                                                \n");
                    printf("                                  KREE!                                         \n");
                    printf("                                                                                \n");
                    printf("================================================================================\n");
                break;
                case 3:
                    printf("================================================================================\n");
                    printf("                                                                                \n");
                    printf("                                                                                \n");
                    printf("      ____        ____      __  __   ____        _____     __   __  __   _____  \n");
                    printf("     |__  \\      /___ \\    / / / /  |___ \\      /____ \\   / /  / / / /  / ___/  \n");
                    printf("     __ | |     __  / /   / / / /   __ | |     _____/ /  / /  / / / /  | |__    \n");
                    printf("    / / | |    / / / /   / / / /   / / | |    / __  _/  / /  / / / /    \\__ \\   \n");
                    printf("   / /  | |   / / / /   / / / /   / /  | |   / / | |   / /  / / / /        \\ \\  \n");
                    printf("  / /___| |  / /_/ /_  / /_/ /   / /___| |  / /  / /  / /  / /_/ /  _______/ |  \n");
                    printf(" /________|  \\______/  \\____/   /________| /_/  /_/  /_/   \\____/  /________/   \n");
                    printf("                                                                                \n");
                    printf("                                                                                \n");
                    printf("                      Where No Computation Has Gone Before                      \n");
                    printf("                                                                                \n");
                    printf("                                                                                \n");
                    printf("================================================================================\n");
                break;
                case 4:
                    printf("================================================================================\n");
                    printf("                                                                                \n");
                    printf("                ============                                                    \n");
                    printf("              //            \\\\                                                  \n");
                    printf("             //              \\\\                                                 \n");
                    printf("     /\\     //      ====      \\\\  ||  ||     /\\      ====    ||   ||  ||  /===  \n");
                    printf("    //\\\\    ||     //  \\\\     ||  ||  ||    //\\\\    ||  \\\\   ||   ||  ||  \\\\    \n");
                    printf("   //==\\\\   ||     ||  ||     ||  ||  ||   //==\\\\   ||==/    ||   ||  ||    \\\\  \n");
                    printf("  //    \\\\  ||     \\\\__//     ||   \\==/   //    \\\\  ||  \\\\   ||    \\==/   ===/  \n");
                    printf("            \\\\        \\\\      //                                                \n");
                    printf("             \\\\              //                                                 \n");
                    printf("              \\\\            //           I WANT TO COMPUTE                      \n");
                    printf("                ============                                                    \n");
                    printf("                                                                                \n");
                    printf("================================================================================\n");
                break;
            }
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

    #ifdef HAVE_LIBINT2
    libint2_static_cleanup();
    #endif

    #ifdef HAVE_ELEMENTAL
    El::Finalize();
    #else
    MPI_Finalize();
    #endif

    return status;
}

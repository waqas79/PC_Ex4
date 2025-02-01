#include <mpi.h>
#include <iostream>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank == 0) {
        std::cout << "Hello from Master " << world_rank << std::endl;
    } else {
        std::cout << "Hello from Worker " << world_rank << std::endl;
    }

    MPI_Finalize();
    return 0;
}

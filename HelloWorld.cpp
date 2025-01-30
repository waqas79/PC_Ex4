#include <mpi.h>
#include <iostream>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);  // Initialize MPI

    int world_rank;
    int world_size;

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  // Get process rank
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);  // Get number of processes

    std::cout << "Hello from process " << world_rank << " out of " << world_size << " processes." << std::endl;

    MPI_Finalize();  // Finalize MPI
    return 0;
}

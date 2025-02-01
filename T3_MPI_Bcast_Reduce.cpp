#include <mpi.h>
#include <iostream>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int number;
    if (rank == 0) {
        // Root process initializes the number
        number = 10;
        std::cout << "Process 0 broadcasting number: " << number << std::endl;
    }

    // Broadcast number to all processes
    MPI_Bcast(&number, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Each process increments the number by its rank
    number += rank;

    // Output the incremented number from each process
    std::cout << "Process " << rank << " incremented number: " << number << std::endl;

    int sum;
    // Reduce all numbers to the sum at the root process
    MPI_Reduce(&number, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // Only the root process prints the sum
        std::cout << "Sum of all incremented numbers: " << sum << std::endl;
    }

    MPI_Finalize();
    return 0;
}


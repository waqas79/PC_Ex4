
#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <vector>
#include <cmath>

#define LOWER 1          // Lower bound of search range
#define UPPER 10000000   // Upper bound of search range

// Function to check if a number is prime
bool is_prime(int n) {
    if (n < 2) return false;
    if (n == 2 || n == 3) return true;
    if (n % 2 == 0) return false;

    for (int i = 3; i <= sqrt(n); i += 2) {
        if (n % i == 0) return false;
    }
    return true;
}

int main(int argc, char* argv[]) {
    int rank, size;
    int local_start, local_end;
    double start_time, end_time;
    std::vector<int> local_primes, global_primes;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Divide the range among processes
    int range_size = (UPPER - LOWER + 1) / size;
    local_start = LOWER + rank * range_size;
    local_end = (rank == size - 1) ? UPPER : local_start + range_size - 1;

    // Start timer
    MPI_Barrier(MPI_COMM_WORLD);
    start_time = MPI_Wtime();

    // Parallel prime search using OpenMP
    #pragma omp parallel for schedule(dynamic)
    for (int i = local_start; i <= local_end; i++) {
        if (is_prime(i)) {
            #pragma omp critical
            local_primes.push_back(i);
        }
    }

    // Gather prime numbers from all processes at rank 0
    if (rank == 0) {
        global_primes.insert(global_primes.end(), local_primes.begin(), local_primes.end());
        for (int i = 1; i < size; i++) {
            int recv_size;
            MPI_Recv(&recv_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::vector<int> recv_primes(recv_size);
            MPI_Recv(recv_primes.data(), recv_size, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            global_primes.insert(global_primes.end(), recv_primes.begin(), recv_primes.end());
        }
    } else {
        int send_size = local_primes.size();
        MPI_Send(&send_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(local_primes.data(), send_size, MPI_INT, 0, 1, MPI_COMM_WORLD);
    }

    // Stop timer
    end_time = MPI_Wtime();

    // Print results
    if (rank == 0) {
        std::cout << "Total Primes Found: " << global_primes.size() << std::endl;
        std::cout << "Execution Time: " << (end_time - start_time) << " seconds" << std::endl;
    }

    MPI_Finalize();
    return 0;
}

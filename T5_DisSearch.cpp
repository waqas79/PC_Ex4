#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstring>

#define N 99999999  // Array size (large for testing)
#define SEARCH_VALUE 9999  // Target value to search for

// Function to initialize the array with random values
void init_array(std::vector<int> &array) {
    for (size_t i = 0; i < array.size(); i++) {
        array[i] = rand() % 9999;  
    }
}

// Parallel search with OpenMP inside MPI process
bool parallel_search(const std::vector<int> &array, int start, int end, int value) {
    bool found = false;
    #pragma omp parallel for shared(found)
    for (int i = start; i < end; i++) {
        if (array[i] == value) {
            found = true;
        }
    }
    return found;
}

// Sequential search (no MPI, no OpenMP)
bool sequential_search(const std::vector<int> &array, int value) {
    for (size_t i = 0; i < array.size(); i++) {
        if (array[i] == value) {
            return true;
        }
    }
    return false;
}

int main(int argc, char *argv[]) {
    int rank = 0, size = 1;
    bool is_sequential = false;
    bool local_found = false, global_found = false;
    std::vector<int> array(N);
    std::vector<int> local_array;
    double start_time, end_time;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Check if "s" argument is passed for sequential execution
    if (argc > 1 && strcmp(argv[1], "s") == 0) {
        is_sequential = true;
    }

    if (is_sequential && rank == 0) {
        // Run Sequential Search (No MPI, No OpenMP)
        std::cout << "Running Sequential Search..." << std::endl;
        init_array(array);
        start_time = omp_get_wtime();
        global_found = sequential_search(array, SEARCH_VALUE);
        end_time = omp_get_wtime();
    } else {
        // Run Parallel Search (MPI + OpenMP)
        if (rank == 0) {
            init_array(array);
        }

        int local_size = N / size;
        local_array.resize(local_size);

        MPI_Scatter(array.data(), local_size, MPI_INT, local_array.data(), local_size, MPI_INT, 0, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);
        start_time = MPI_Wtime();
        local_found = parallel_search(local_array, 0, local_size, SEARCH_VALUE);
        end_time = MPI_Wtime();

        // Reduce results to check if any process found the value
        MPI_Reduce(&local_found, &global_found, 1, MPI_CXX_BOOL, MPI_LOR, 0, MPI_COMM_WORLD);
    }

    // Print execution time in root process
    if (rank == 0) {
        std::cout << "Search Completed. Execution Time: " << (end_time - start_time) << " seconds" << std::endl;
        if (global_found) {
            std::cout << "Value " << SEARCH_VALUE << " found in the array!" << std::endl;
        } else {
            std::cout << "Value " << SEARCH_VALUE << " not found in the array." << std::endl;
        }
    }

    MPI_Finalize();
    return 0;
}


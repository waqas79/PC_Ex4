#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 1200  // Matrix size

// Allocate a 1D contiguous memory block for a matrix
double* allocate_matrix(int size) {
    double* matrix = (double*) malloc(size * size * sizeof(double));
    if (matrix == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    return matrix;
}

// Free allocated matrix
void free_matrix(double* matrix) {
    free(matrix);
}

// Parallel matrix multiplication using OpenMP
void mat_mult_parallel(double* A, double* B, double* C, int size, int row_start, int row_end) {
    #pragma omp parallel for collapse(2)
    for (int i = row_start; i < row_end; i++) {
        for (int j = 0; j < size; j++) {
            double sum = 0.0;
            for (int k = 0; k < size; k++) {
                sum += A[i * size + k] * B[k * size + j];
            }
            C[i * size + j] = sum;
        }
    }
}

// Sequential matrix multiplication
void mat_mult_sequential(double* A, double* B, double* C, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            double sum = 0.0;
            for (int k = 0; k < size; k++) {
                sum += A[i * size + k] * B[k * size + j];
            }
            C[i * size + j] = sum;
        }
    }
}

int main(int argc, char* argv[]) {
    int rank = 0, size = 1;
    double *A, *B, *C;
    double start_time, end_time;
    int is_sequential = 0;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Check if "s" argument is passed for sequential execution
    if (argc > 1 && strcmp(argv[1], "s") == 0) {
        is_sequential = 1;
    }

    A = allocate_matrix(N);
    B = allocate_matrix(N);
    C = allocate_matrix(N);

    // Initialize matrices A and B
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i * N + j] = (i == j) ? 1.0 : 0.0;  // Identity matrices
            B[i * N + j] = (i == j) ? 1.0 : 0.0;
            C[i * N + j] = 0.0;
        }
    }

    if (is_sequential && rank == 0) {
        printf("Running Sequential Matrix Multiplication...\n");
        start_time = omp_get_wtime();
        mat_mult_sequential(A, B, C, N);
        end_time = omp_get_wtime();
    } else {
        // Parallel Execution (MPI + OpenMP)
        int rows_per_proc = N / size;
        int row_start = rank * rows_per_proc;
        int row_end = (rank + 1) * rows_per_proc;

        MPI_Barrier(MPI_COMM_WORLD);  // Synchronize before timing
        start_time = MPI_Wtime();
        mat_mult_parallel(A, B, C, N, row_start, row_end);
        end_time = MPI_Wtime();

        // Gather results using correct memory layout
        MPI_Gather(C + row_start * N, rows_per_proc * N, MPI_DOUBLE, C, rows_per_proc * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Print execution time
    if (rank == 0) {
        printf("Matrix multiplication completed. C[100][100] = %f\n", C[100 * N + 100]);
        printf("Execution Time: %f seconds\n", end_time - start_time);
    }

    free_matrix(A);
    free_matrix(B);
    free_matrix(C);

    MPI_Finalize();
    return 0;
}

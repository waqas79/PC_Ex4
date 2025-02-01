#include <mpi.h>
#include <iostream>
#include <iomanip> // Required for std::setprecision
#include <unistd.h> // For usleep function

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int ping_pong_count = 0;
    int partner_rank = (world_rank + 1) % 2;

    if (world_rank == 0) {
        // Process 0 sends two initial pings
        double current_time = MPI_Wtime();
        MPI_Send(&ping_pong_count, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD);
        std::cout << "Process " << world_rank << " sent initial ping " << ping_pong_count
                  << " to " << partner_rank << " at " << std::fixed << std::setprecision(4) << current_time << std::endl;

        ping_pong_count++;
        current_time = MPI_Wtime();
        MPI_Send(&ping_pong_count, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD);
        std::cout << "Process " << world_rank << " sent second ping " << ping_pong_count
                  << " to " << partner_rank << " at " << std::fixed << std::setprecision(4) << current_time << std::endl;

        ping_pong_count++;
    }

    while (ping_pong_count < 10) {
        double current_time;
        if (world_rank == ping_pong_count % 2) {
            usleep(1000); // Artificial delay for clearer time differentiation
            current_time = MPI_Wtime();
            MPI_Send(&ping_pong_count, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD);
            std::cout << "Process " << world_rank << " sent ping_pong_count " << ping_pong_count
                      << " to " << partner_rank << " at " << std::fixed << std::setprecision(4) << current_time << std::endl;
        } else {
            usleep(1000); // Artificial delay
            current_time = MPI_Wtime();
            MPI_Recv(&ping_pong_count, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::cout << "Process " << world_rank << " received ping_pong_count " << ping_pong_count
                      << " from " << partner_rank << " at " << std::fixed << std::setprecision(4) << current_time << std::endl;
        }
        ping_pong_count++;
    }

    MPI_Finalize();
    return 0;
}


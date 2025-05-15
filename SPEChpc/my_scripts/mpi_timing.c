#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include "mpi_timing.h"

#define MAX_MPI_FUNCTIONS 20

static MPITimingInfo timing_data[MAX_MPI_FUNCTIONS];
static int timing_entries = 0;
static int timing_initialized = 0;

static void update_timing(const char *function_name, double elapsed_time) {
    int i;
    
    // Find existing entry or create new one
    for (i = 0; i < timing_entries; i++) {
        if (strcmp(timing_data[i].function_name, function_name) == 0) {
            timing_data[i].total_time += elapsed_time;
            timing_data[i].call_count++;
            if (elapsed_time < timing_data[i].min_time)
                timing_data[i].min_time = elapsed_time;
            if (elapsed_time > timing_data[i].max_time)
                timing_data[i].max_time = elapsed_time;
            return;
        }
    }
    
    // Add new entry
    if (timing_entries < MAX_MPI_FUNCTIONS) {
        strncpy(timing_data[timing_entries].function_name, function_name, 63);
        timing_data[timing_entries].function_name[63] = '\0';
        timing_data[timing_entries].total_time = elapsed_time;
        timing_data[timing_entries].call_count = 1;
        timing_data[timing_entries].min_time = elapsed_time;
        timing_data[timing_entries].max_time = elapsed_time;
        timing_entries++;
    }
}

void init_mpi_timing(void) {
    if (!timing_initialized) {
        memset(timing_data, 0, sizeof(timing_data));
        timing_entries = 0;
        timing_initialized = 1;
    }
}

double get_total_mpi_time(void) {
    double total = 0.0;
    int i;
    
    for (i = 0; i < timing_entries; i++) {
        total += timing_data[i].total_time;
    }
    return total;
}

void finalize_mpi_timing(int rank, double total_execution_time) {
    int num_ranks;
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    
    // Collect data from all ranks
    MPITimingInfo all_timing[MAX_MPI_FUNCTIONS * num_ranks];
    int all_entries[num_ranks];
    
    // Gather number of entries from each rank
    MPI_Gather(&timing_entries, 1, MPI_INT, all_entries, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Prepare send data
    MPITimingInfo send_buffer[MAX_MPI_FUNCTIONS];
    memcpy(send_buffer, timing_data, timing_entries * sizeof(MPITimingInfo));
    
    // Gather timing data from all ranks
    if (rank == 0) {
        // Rank 0 collects data from all ranks
        int displacements[num_ranks];
        int recv_counts[num_ranks];
        displacements[0] = 0;
        recv_counts[0] = all_entries[0] * sizeof(MPITimingInfo);
        
        for (int i = 1; i < num_ranks; i++) {
            displacements[i] = displacements[i-1] + recv_counts[i-1];
            recv_counts[i] = all_entries[i] * sizeof(MPITimingInfo);
        }
        
        MPI_Gatherv(send_buffer, timing_entries * sizeof(MPITimingInfo), MPI_BYTE,
                    all_timing, recv_counts, displacements, MPI_BYTE, 0, MPI_COMM_WORLD);
    } else {
        MPI_Gatherv(send_buffer, timing_entries * sizeof(MPITimingInfo), MPI_BYTE,
                    NULL, NULL, NULL, MPI_BYTE, 0, MPI_COMM_WORLD);
    }
    
    // Rank 0 processes and reports aggregated data
    if (rank == 0) {
        // Aggregate data across all ranks
        MPITimingInfo aggregated[MAX_MPI_FUNCTIONS];
        int aggregated_entries = 0;
        
        // Process data from all ranks
        int offset = 0;
        for (int r = 0; r < num_ranks; r++) {
            for (int i = 0; i < all_entries[r]; i++) {
                MPITimingInfo *current = &all_timing[offset + i];
                
                // Find or create entry in aggregated data
                int found = -1;
                for (int j = 0; j < aggregated_entries; j++) {
                    if (strcmp(aggregated[j].function_name, current->function_name) == 0) {
                        found = j;
                        break;
                    }
                }
                
                if (found >= 0) {
                    // Update existing entry
                    aggregated[found].total_time += current->total_time;
                    aggregated[found].call_count += current->call_count;
                    if (current->min_time < aggregated[found].min_time)
                        aggregated[found].min_time = current->min_time;
                    if (current->max_time > aggregated[found].max_time)
                        aggregated[found].max_time = current->max_time;
                } else {
                    // Add new entry
                    aggregated[aggregated_entries] = *current;
                    aggregated_entries++;
                }
            }
            offset += all_entries[r];
        }
        
        // Calculate total MPI time across all ranks
        double total_mpi_time = 0.0;
        for (int i = 0; i < aggregated_entries; i++) {
            total_mpi_time += aggregated[i].total_time;
        }
        
        // Calculate average MPI time per rank
        double avg_mpi_time_per_rank = total_mpi_time / num_ranks;
        double percentage = (avg_mpi_time_per_rank / total_execution_time) * 100.0;
        
        // Print aggregated report
        printf("\n========================================\n");
        printf("AGGREGATED MPI TIMING REPORT (ALL RANKS)\n");
        printf("========================================\n");
        printf("Number of MPI Ranks: %d\n", num_ranks);
        printf("Total Execution Time: %.6f seconds\n", total_execution_time);
        printf("Total MPI Time (all ranks): %.6f seconds\n", total_mpi_time);
        printf("Average MPI Time per rank: %.6f seconds (%.2f%% of execution time)\n\n", 
               avg_mpi_time_per_rank, percentage);
        
        printf("MPI Call Summary (Aggregated Across All Ranks):\n");
        printf("%-20s | %12s | %14s | %14s | %14s | %14s | %8s\n", 
               "Function", "Total Calls", "Total Time (s)", "Avg Time (s)", 
               "Min Time (s)", "Max Time (s)", "% of MPI");
        printf("--------------------------------------------------------------------------------------------------------\n");
        
        // Sort by total time for better readability
        for (int i = 0; i < aggregated_entries - 1; i++) {
            for (int j = i + 1; j < aggregated_entries; j++) {
                if (aggregated[j].total_time > aggregated[i].total_time) {
                    MPITimingInfo temp = aggregated[i];
                    aggregated[i] = aggregated[j];
                    aggregated[j] = temp;
                }
            }
        }
        
        // Print sorted results
        for (int i = 0; i < aggregated_entries; i++) {
            double avg_time = aggregated[i].total_time / aggregated[i].call_count;
            double percentage_of_mpi = (aggregated[i].total_time / total_mpi_time) * 100.0;
            
            printf("%-20s | %12d | %14.6f | %14.9f | %14.9f | %14.9f | %7.4f%%\n",
                   aggregated[i].function_name,
                   aggregated[i].call_count,
                   aggregated[i].total_time,
                   avg_time,
                   aggregated[i].min_time,
                   aggregated[i].max_time,
                   percentage_of_mpi);
        }
        
        printf("\nNOTE: Total time represents the sum across all ranks.\n");
        printf("      Average time is per individual call.\n");
        printf("      Min/Max times are across all calls on all ranks.\n");
    }
}

// Wrapper implementations remain the same
int timed_MPI_Irecv(void *buf, int count, MPI_Datatype datatype, 
                    int source, int tag, MPI_Comm comm, MPI_Request *request) {
    double start_time = MPI_Wtime();
    int result = MPI_Irecv(buf, count, datatype, source, tag, comm, request);
    double end_time = MPI_Wtime();
    
    update_timing("MPI_Irecv", end_time - start_time);
    return result;
}

int timed_MPI_Isend(const void *buf, int count, MPI_Datatype datatype, 
                    int dest, int tag, MPI_Comm comm, MPI_Request *request) {
    double start_time = MPI_Wtime();
    int result = MPI_Isend(buf, count, datatype, dest, tag, comm, request);
    double end_time = MPI_Wtime();
    
    update_timing("MPI_Isend", end_time - start_time);
    return result;
}

int timed_MPI_Wait(MPI_Request *request, MPI_Status *status) {
    double start_time = MPI_Wtime();
    int result = MPI_Wait(request, status);
    double end_time = MPI_Wtime();
    
    update_timing("MPI_Wait", end_time - start_time);
    return result;
}

int timed_MPI_Reduce(const void *sendbuf, void *recvbuf, int count, 
                     MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
    double start_time = MPI_Wtime();
    int result = MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
    double end_time = MPI_Wtime();
    
    update_timing("MPI_Reduce", end_time - start_time);
    return result;
}

int timed_MPI_Barrier(MPI_Comm comm) {
    double start_time = MPI_Wtime();
    int result = MPI_Barrier(comm);
    double end_time = MPI_Wtime();
    
    update_timing("MPI_Barrier", end_time - start_time);
    return result;
}

int timed_MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
    double start_time = MPI_Wtime();
    int result = MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
    double end_time = MPI_Wtime();

    update_timing("MPI_Allreduce", end_time - start_time);
    return result;
}

int timed_MPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                        void *recvbuf, int recvcount, MPI_Datatype recvtype,
                        MPI_Comm comm) {
    double start_time = MPI_Wtime();
    int result = MPI_Allgather(sendbuf, sendcount, sendtype,
                              recvbuf, recvcount, recvtype, comm);
    double end_time = MPI_Wtime();

    update_timing("MPI_Allgather", end_time - start_time);
    return result;
}

int timed_MPI_Recv(void *buf, int count, MPI_Datatype datatype,
                   int source, int tag, MPI_Comm comm, MPI_Status *status) {
    double start_time = MPI_Wtime();
    int result = MPI_Recv(buf, count, datatype, source, tag, comm, status);
    double end_time = MPI_Wtime();

    update_timing("MPI_Recv", end_time - start_time);
    return result;
}

int timed_MPI_Send(const void *buf, int count, MPI_Datatype datatype,
                   int dest, int tag, MPI_Comm comm) {
    double start_time = MPI_Wtime();
    int result = MPI_Send(buf, count, datatype, dest, tag, comm);
    double end_time = MPI_Wtime();

    update_timing("MPI_Send", end_time - start_time);
    return result;
}

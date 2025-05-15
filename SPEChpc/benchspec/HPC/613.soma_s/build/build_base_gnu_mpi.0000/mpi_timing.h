#ifndef MPI_TIMING_H
#define MPI_TIMING_H

#include <mpi.h>

typedef struct {
    char function_name[64];
    double total_time;
    int call_count;
    double min_time;
    double max_time;
} MPITimingInfo;

// Initialize the timing system
void init_mpi_timing(void);

// Finalize and report timing
void finalize_mpi_timing(int rank, double total_execution_time);

// Wrapper function declarations
int timed_MPI_Irecv(void *buf, int count, MPI_Datatype datatype, 
                    int source, int tag, MPI_Comm comm, MPI_Request *request);

int timed_MPI_Isend(const void *buf, int count, MPI_Datatype datatype, 
                    int dest, int tag, MPI_Comm comm, MPI_Request *request);

int timed_MPI_Wait(MPI_Request *request, MPI_Status *status);

int timed_MPI_Reduce(const void *sendbuf, void *recvbuf, int count, 
                     MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);

int timed_MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

int timed_MPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                        void *recvbuf, int recvcount, MPI_Datatype recvtype,
                        MPI_Comm comm);

int timed_MPI_Recv(void *buf, int count, MPI_Datatype datatype,
                   int source, int tag, MPI_Comm comm, MPI_Status *status);

int timed_MPI_Send(const void *buf, int count, MPI_Datatype datatype,
                   int dest, int tag, MPI_Comm comm);

int timed_MPI_Barrier(MPI_Comm comm);

// Get total MPI time for this process
double get_total_mpi_time(void);

#endif

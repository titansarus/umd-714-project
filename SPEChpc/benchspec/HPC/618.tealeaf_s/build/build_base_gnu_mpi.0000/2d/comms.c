#include "comms.h"
#include "settings.h"

#ifdef SPEC_OPENACC
#include <openacc.h>
#endif
#if defined(SPEC_OPENMP) || defined(SPEC_OPENMP_TARGET)
#include "omp.h"
#endif

#ifndef NO_MPI

// Initialise MPI
void initialise_comms(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
}

// Initialise the rank information
void initialise_ranks(Settings* settings)
{
  MPI_Comm_rank(MPI_COMM_WORLD, &settings->rank);
  MPI_Comm_size(MPI_COMM_WORLD, &settings->num_ranks);
#if defined(SPEC_CUDA) || defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
  int num_devices;
  int gpuId;
  MPI_Comm shmcomm;
  int local_rank;
#endif
#ifdef SPEC_OPENACC
  acc_device_t my_device_type;
#endif

// Set the device number to the mod of the local rank and number of devices.
#if defined(SPEC_CUDA) || defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                    MPI_INFO_NULL, &shmcomm);
    MPI_Comm_rank(shmcomm, &local_rank);
#endif

#ifdef SPEC_CUDA
  cudaGetDeviceCount ( &num_devices );
  gpuId = local_rank % num_devices;
  cudaSetDevice( gpuId );
#endif
#ifdef SPEC_OPENACC
  my_device_type = acc_get_device_type();
  num_devices = acc_get_num_devices(my_device_type);
  gpuId = local_rank % num_devices;
  acc_set_device_num(gpuId, my_device_type);
#endif
#ifdef SPEC_OPENMP_TARGET
  num_devices = omp_get_num_devices();
  if (num_devices > 0) {
    gpuId = local_rank % num_devices;
    omp_set_default_device(gpuId);
  }
#endif


  if(settings->rank == MASTER)
  {
    printf("Successfully initialised %d MPI ranks.\n", settings->num_ranks);
  }
}

// Teardown MPI
void finalise_comms()
{
  MPI_Finalize();
}

// Sends a message out and receives a message in
void send_recv_message(Settings* settings, double* send_buffer, 
    double* recv_buffer, int buffer_len, int neighbour, int send_tag, 
    int recv_tag, MPI_Request* send_request, MPI_Request* recv_request)
{
#ifndef SPEC
  START_PROFILING(settings->kernel_profile);
#endif

#if defined(SPEC_OPENACC) && defined(SPEC_ACCEL_AWARE_MPI)
#pragma acc host_data use_device(send_buffer,recv_buffer)
{
#endif
  MPI_Isend(send_buffer, buffer_len, MPI_DOUBLE, 
      neighbour, send_tag, MPI_COMM_WORLD, send_request);
  MPI_Irecv(recv_buffer, buffer_len, MPI_DOUBLE,
      neighbour, recv_tag, MPI_COMM_WORLD, recv_request);

#if defined(SPEC_OPENACC) && defined(SPEC_ACCEL_AWARE_MPI)
}
#endif
#ifndef SPEC
  STOP_PROFILING(settings->kernel_profile, __func__);
#endif
}

// Waits for all requests to complete
void wait_for_requests(
    Settings* settings, int num_requests, MPI_Request* requests)
{
#ifndef SPEC
  START_PROFILING(settings->kernel_profile);
#endif
  MPI_Waitall(num_requests, requests, MPI_STATUSES_IGNORE);
#ifndef SPEC
  STOP_PROFILING(settings->kernel_profile, __func__);
#endif
}

// Reduce over all ranks to get sum
void sum_over_ranks(Settings* settings, double* a)
{
#ifndef SPEC
  START_PROFILING(settings->kernel_profile);
#endif
  double temp = *a;
  MPI_Allreduce(&temp, a, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#ifndef SPEC
  STOP_PROFILING(settings->kernel_profile, __func__);
#endif
}

// Reduce across all ranks to get minimum value
void min_over_ranks(Settings* settings, double* a)
{
#ifndef SPEC
  START_PROFILING(settings->kernel_profile);
#endif
  double temp = *a;
  MPI_Allreduce(&temp, a, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#ifndef SPEC
  STOP_PROFILING(settings->kernel_profile, __func__);
#endif
}

// Synchronise all ranks
void barrier()
{
  MPI_Barrier(MPI_COMM_WORLD);
}

// End the application
void abort_comms()
{
  MPI_Abort(MPI_COMM_WORLD, 1);
}

#else

void initialise_comms(int argc, char** argv) { }
void initialise_ranks(Settings* settings) 
{ 
  settings->rank = MASTER;
  settings->num_ranks = 1;
}
void finalise_comms() { }
void sum_over_ranks(Settings* settings, double* a) { }
void min_over_ranks(Settings* settings, double* a) { }
void barrier() { }
void abort_comms() 
{ 
  exit(1);
}

#endif

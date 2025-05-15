#include "device.h"

#include <stdio.h>
#include <mpi.h>

#ifdef SPEC_OPENMP_TARGET
#include <omp.h>
#endif
#ifdef SPEC_OPENACC
#include <openacc.h>
#endif

int soma_device_num;

int soma_get_device()
{
#if defined SPEC_OPENMP_TARGET && _OPENMP >= 201811 && !defined(SPEC_USE_HOST_THREADS)
    return omp_get_default_device();
#else
    return soma_device_num;
#endif
}

int set_openmp_devices(const struct Phase*const p) {
    int ret=0;
#ifdef SPEC_OPENMP_TARGET
// Determine the local rank ID
// The device type, number of devices of that type on the node.
// Set the device number to the mod of the local rank and number of devices. 
    MPI_Comm shmcomm;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                    MPI_INFO_NULL, &shmcomm);
    int local_rank;
    MPI_Comm_rank(shmcomm, &local_rank);
    int num_devices = omp_get_num_devices();
    if (num_devices > 0) {
      soma_device_num = local_rank % num_devices;
    }
    else
      soma_device_num = omp_get_initial_device();
    // then following will quietly fail if soma_device_num == omp_get_initial_device()
    // which is why we will use soma_get_device() instead of omp_get_default_device()
    omp_set_default_device(soma_device_num);
    printf("Using OpenMP target device %d\n", soma_device_num);
#endif
    return ret;
    }

int set_openacc_devices(const struct Phase*const p)
    {
    int ret=0;
#ifdef SPEC_OPENACC
// Determine the local rank ID
// The device type, number of devices of that type on the node.
// Set the device number to the mod of the local rank and number of devices. 
    MPI_Comm shmcomm;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
                    MPI_INFO_NULL, &shmcomm);
    int local_rank;
    MPI_Comm_rank(shmcomm, &local_rank);
    acc_device_t devtype = acc_get_device_type();
    int numdev = acc_get_num_devices(devtype);
    soma_device_num = local_rank % numdev;
    acc_set_device_num(soma_device_num,devtype); 
    acc_init(devtype);    
    printf("Using OpenACC target device %d\n", soma_device_num);
#endif//SPEC_OPENACC
    return ret;
    }

#include <stdlib.h>
#include "../shared.h"

// Packs left data into buffer.
void pack_left(
        const int x,
        const int y,
        const int depth,
        const int halo_depth,
        double* field,
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
        double* buffer,
        bool is_offload)
#else    //openmp and mpi
        double* buffer)
#endif 
{
#ifdef SPEC_OPENACC
    const int y_inner = y - 2*halo_depth;

#pragma acc kernels loop independent collapse(2) async(0) if(is_offload) \
    present(field[:x*y]) copyout(buffer[:depth*y_inner])
#endif
#ifdef SPEC_OPENMP_TARGET
    const int y_inner = y - 2*halo_depth;
#pragma omp target teams distribute parallel for \
  map(tofrom: field[:x*y]), map(from: buffer[:depth*y_inner]) collapse(2)
#else
#ifdef SPEC_OPENMP
#pragma omp parallel for
#endif
#endif
    for(int jj = halo_depth; jj < y-halo_depth; ++jj)
    {
        for(int kk = halo_depth; kk < halo_depth+depth; ++kk)
        {
            int buf_index = (kk-halo_depth) + (jj-halo_depth)*depth;
            buffer[buf_index] = field[jj*x+kk];
        }
    }
}

// Packs right data into buffer.
void pack_right(
        const int x,
        const int y,
        const int depth,
        const int halo_depth,
        double* field,
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
        double* buffer,
        bool is_offload)
#else    //openmp and mpi
        double* buffer)
#endif 
{
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
    const int y_inner = y - 2*halo_depth;
#endif

#ifdef SPEC_OPENACC
#pragma acc kernels loop independent collapse(2) async(0) if(is_offload) \
    present(field[:x*y]) copyout(buffer[:depth*y_inner])
#endif
#ifdef SPEC_OPENMP_TARGET
#pragma omp target teams distribute parallel for \
  map(tofrom: field[:x*y]), map(from: buffer[:depth*y_inner]) collapse(2)
#else
#ifdef SPEC_OPENMP
#pragma omp parallel for
#endif
#endif
    for(int jj = halo_depth; jj < y-halo_depth; ++jj)
    {
        for(int kk = x-halo_depth-depth; kk < x-halo_depth; ++kk)
        {
            int buf_index = (kk-(x-halo_depth-depth)) + (jj-halo_depth)*depth;
            buffer[buf_index] = field[jj*x+kk];
        }
    }
}

// Packs top data into buffer.
void pack_top(
        const int x,
        const int y,
        const int depth,
        const int halo_depth,
        double* field,
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
        double* buffer,
        bool is_offload)
#else    //openmp and mpi
        double* buffer)
#endif
{
    const int x_inner = x-2*halo_depth;

#ifdef SPEC_OPENACC
#pragma acc kernels loop independent collapse(2) async(0) if(is_offload) \
    present(field[:x*y]) copyout(buffer[:depth*x_inner])
#endif
#ifdef SPEC_OPENMP_TARGET
#pragma omp target teams distribute parallel for\
  map(tofrom: field[:x*y]), map(from: buffer[:depth*x_inner]) collapse(2)
#else
#ifdef SPEC_OPENMP
#pragma omp parallel for
#endif
#endif
    for(int jj = y-halo_depth-depth; jj < y-halo_depth; ++jj)
    {
        for(int kk = halo_depth; kk < x-halo_depth; ++kk)
        {
            int buf_index = (kk-halo_depth) + (jj-(y-halo_depth-depth))*x_inner;
            buffer[buf_index] = field[jj*x+kk];
        }
    }
}

// Packs bottom data into buffer.
void pack_bottom(
        const int x,
        const int y,
        const int depth,
        const int halo_depth,
        double* field,
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
        double* buffer,
        bool is_offload)
#else    //openmp and mpi
        double* buffer)
#endif 
{
    const int x_inner = x-2*halo_depth;

#ifdef SPEC_OPENACC
#pragma acc kernels loop independent collapse(2) async(0) if(is_offload) \
    present(field[:x*y]) copyout(buffer[:depth*x_inner])
#endif
#ifdef SPEC_OPENMP_TARGET
#pragma omp target teams distribute parallel for \
  map(tofrom: field[:x*y]), map(from: buffer[:depth*x_inner]) collapse(2)
#else
#ifdef SPEC_OPENMP
#pragma omp parallel for 
#endif
#endif
    for(int jj = halo_depth; jj < halo_depth+depth; ++jj)
    {
        for(int kk = halo_depth; kk < x-halo_depth; ++kk)
        {
            int buf_index = (kk-halo_depth) + (jj-halo_depth)*x_inner;
            buffer[buf_index] = field[jj*x+kk];
        }
    }
}

// Unpacks left data from buffer.
void unpack_left(
        const int x,
        const int y,
        const int depth,
        const int halo_depth,
        double* field,
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
        double* buffer,
        bool is_offload)
#else    //openmp and mpi
        double* buffer)
#endif
{
#ifdef SPEC_OPENACC
    const int y_inner = y - 2*halo_depth;

#pragma acc kernels loop independent collapse(2) async(0) if(is_offload) \
    present(field[:x*y]) copyin(buffer[:depth*y_inner])
#endif
#ifdef SPEC_OPENMP_TARGET
    const int y_inner = y - 2*halo_depth;
#pragma omp target teams distribute parallel for \
  map(tofrom: field[:x*y]), map(to: buffer[:depth*y_inner]) collapse(2)
#else
#ifdef SPEC_OPENMP
#pragma omp parallel for 
#endif
#endif
    for(int jj = halo_depth; jj < y-halo_depth; ++jj)
    {
        for(int kk = halo_depth-depth; kk < halo_depth; ++kk)
        {
            int buf_index = (kk-(halo_depth-depth)) + (jj-halo_depth)*depth;
            field[jj*x+kk] = buffer[buf_index];
        }
    }
}

// Unpacks right data from buffer.
void unpack_right(
        const int x,
        const int y,
        const int depth,
        const int halo_depth,
        double* field,
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
        double* buffer,
        bool is_offload)
#else    //openmp and mpi
        double* buffer)
#endif
{
#ifdef SPEC_OPENACC
    const int y_inner = y - 2*halo_depth;

#pragma acc kernels loop independent collapse(2) async(0) if(is_offload) \
    present(field[:x*y]) copyin(buffer[:depth*y_inner])
#endif
#ifdef SPEC_OPENMP_TARGET
    const int y_inner = y - 2*halo_depth;
#pragma omp target teams distribute parallel for \
  map(tofrom: field[:x*y]), map(to: buffer[:depth*y_inner]) collapse(2)
#else
#ifdef SPEC_OPENMP
#pragma omp parallel for
#endif
#endif
    for(int jj = halo_depth; jj < y-halo_depth; ++jj)
    {
        for(int kk = x-halo_depth; kk < x-halo_depth+depth; ++kk)
        {
            int buf_index = (kk-(x-halo_depth)) + (jj-halo_depth)*depth;
            field[jj*x+kk] = buffer[buf_index];
        }
    }
}

// Unpacks top data from buffer.
void unpack_top(
        const int x,
        const int y,
        const int depth,
        const int halo_depth,
        double* field,
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
        double* buffer,
        bool is_offload)
#else    //openmp and mpi
        double* buffer)
#endif 
{
    const int x_inner = x-2*halo_depth;

#ifdef SPEC_OPENACC
#pragma acc kernels loop independent collapse(2) async(0) if(is_offload) \
    present(field[:x*y]) copyin(buffer[:depth*x_inner])
#endif
#ifdef SPEC_OPENMP_TARGET
#pragma omp target teams distribute parallel for \
  map(tofrom: field[:x*y]), map(to: buffer[:depth*x_inner]) collapse(2)
#else
#ifdef SPEC_OPENMP
#pragma omp parallel for
#endif
#endif
    for(int jj = y-halo_depth; jj < y-halo_depth+depth; ++jj)
    {
        for(int kk = halo_depth; kk < x-halo_depth; ++kk)
        {
            int buf_index = (kk-halo_depth) + (jj-(y-halo_depth))*x_inner;
            field[jj*x+kk] = buffer[buf_index];
        }
    }
}

// Unpacks bottom data from buffer.
void unpack_bottom(
        const int x,
        const int y,
        const int depth,
        const int halo_depth,
        double* field,
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
        double* buffer,
        bool is_offload)
#else    //openmp and mpi
        double* buffer)
#endif
{
    const int x_inner = x-2*halo_depth;

#ifdef SPEC_OPENACC
#pragma acc kernels loop independent collapse(2) async(0) if(is_offload) \
    present(field[:x*y]) copyin(buffer[:depth*x_inner])
#endif
#ifdef SPEC_OPENMP_TARGET
#pragma omp target teams distribute parallel for \
  map(tofrom: field[:x*y]), map(to: buffer[:depth*x_inner]) collapse(2)
#else
#ifdef SPEC_OPENMP
#pragma omp parallel for
#endif
#endif
    for(int jj = halo_depth-depth; jj < halo_depth; ++jj)
    {
        for(int kk = halo_depth; kk < x-halo_depth; ++kk)
        {
            int buf_index = (kk-halo_depth) + (jj-(halo_depth-depth))*x_inner;
            field[jj*x+kk] = buffer[buf_index];
        }
    }
}

#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
typedef void (*pack_kernel_f)(int,int,int,int,double*,double*,bool);
#else    //openmp and mpi
typedef void (*pack_kernel_f)(int,int,int,int,double*,double*);
#endif 

// Either packs or unpacks data from/to buffers.
void pack_or_unpack(
        const int x,
        const int y,
        const int depth,
        const int halo_depth,
        const int face,
        bool pack,
        double *field,
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
        double* buffer,
        bool is_offload)
#else    //openmp and mpi
        double* buffer)
#endif
{
    pack_kernel_f kernel = NULL;

    switch(face)
    {
        case CHUNK_LEFT:
            kernel = pack ? pack_left : unpack_left;
            break;
        case CHUNK_RIGHT:
            kernel = pack ? pack_right : unpack_right;
            break;
        case CHUNK_TOP:
            kernel = pack ? pack_top : unpack_top;
            break;
        case CHUNK_BOTTOM:
            kernel = pack ? pack_bottom : unpack_bottom;
            break;
        default:
            die(__LINE__, __FILE__, "Incorrect face provided: %d.\n", face);
    }

#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
    kernel(x, y, depth, halo_depth, field, buffer, is_offload);
#else    //openmp and mpi
    kernel(x, y, depth, halo_depth, field, buffer);
#endif
}

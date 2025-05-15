#include <stdlib.h>
#include "../shared.h"

/*
 * 		LOCAL HALOS KERNEL
 */	

void update_left(
        const int x, const int y, const int halo_depth, const int depth, 
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
        double* buffer, bool is_offload);
#else    //openmp and mpi
        double* buffer);
#endif
void update_right(
        const int x, const int y, const int halo_depth, const int depth, 
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
        double* buffer, bool is_offload);
#else    //openmp and mpi
        double* buffer);
#endif
void update_top(
        const int x, const int y, const int halo_depth, const int depth, 
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
        double* buffer, bool is_offload); 
#else    //openmp and mpi
        double* buffer); 
#endif

void update_bottom(
        const int x, const int y, const int halo_depth, const int depth, 
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
        double* buffer, bool is_offload);
#else    //openmp and mpi
        double* buffer);
#endif
void update_face(
        const int x, const int y, const int halo_depth, 
        const int* chunk_neighbours, const int depth, 
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
        double* buffer, bool is_offload);
#else    //openmp and mpi
        double* buffer);
#endif

typedef void (*update_kernel)(int,double*);

// The kernel for updating halos locally
void local_halos(
        const int x,
        const int y,
        const int depth,
        const int halo_depth,
        const int* chunk_neighbours,
        const bool* fields_to_exchange,
        double* density,
        double* energy0,
        double* energy,
        double* u,
        double* p,
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
        double* sd,
        bool is_offload)
#else    //openmp and mpi
        double* sd)
#endif 

{
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
#define LAUNCH_UPDATE(index, buffer)\
    if(fields_to_exchange[index])\
    {\
        update_face(x, y, halo_depth, chunk_neighbours, depth, buffer, is_offload);\
    }
#else    //openmp and mpi
#define LAUNCH_UPDATE(index, buffer)\
    if(fields_to_exchange[index])\
    {\
        update_face(x, y, halo_depth, chunk_neighbours, depth, buffer);\
    }
#endif


    LAUNCH_UPDATE(FIELD_DENSITY, density);
    LAUNCH_UPDATE(FIELD_P, p);
    LAUNCH_UPDATE(FIELD_ENERGY0, energy0);
    LAUNCH_UPDATE(FIELD_ENERGY1, energy);
    LAUNCH_UPDATE(FIELD_U, u);
    LAUNCH_UPDATE(FIELD_SD, sd);
#undef LAUNCH_UPDATE
}

// Updates faces in turn.
void update_face(
        const int x,
        const int y, 
        const int halo_depth,
        const int* chunk_neighbours,
        const int depth,
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
        double* buffer, 
        bool is_offload)
#else    //openmp and mpi
        double* buffer)
#endif
{
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
#define UPDATE_FACE(face, updateKernel) \
    if(chunk_neighbours[face] == EXTERNAL_FACE)\
    {\
        updateKernel(x, y, halo_depth, depth, buffer, is_offload);\
    }
#else    //openmp and mpi
#define UPDATE_FACE(face, updateKernel) \
    if(chunk_neighbours[face] == EXTERNAL_FACE)\
    {\
        updateKernel(x, y, halo_depth, depth,buffer);\
    }
#endif


    UPDATE_FACE(CHUNK_LEFT, update_left);
    UPDATE_FACE(CHUNK_RIGHT, update_right);
    UPDATE_FACE(CHUNK_TOP, update_top);
    UPDATE_FACE(CHUNK_BOTTOM, update_bottom);
}

// Update left halo.
void update_left(
        const int x,
        const int y,
        const int halo_depth,
        const int depth, 
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
        double* buffer,
        bool is_offload)
#else    //openmp and mpi
        double* buffer)
#endif
{
#ifdef SPEC_OPENACC
#pragma acc kernels loop independent collapse(2) async(0) \
    if(is_offload) present(buffer[:x*y])
#endif
#ifdef SPEC_OPENMP_TARGET
#pragma omp target teams distribute parallel for collapse(2) map(tofrom: buffer[:x*y])
#else
#ifdef SPEC_OPENMP
#pragma omp parallel for 
#endif
#endif
    for(int jj = halo_depth; jj < y-halo_depth; ++jj)
    {
        for(int kk = 0; kk < depth; ++kk)
        {
            int base = jj*x;
            buffer[base+(halo_depth-kk-1)] = buffer[base+(halo_depth+kk)];			
        }
    }
}

// Update right halo.
void update_right(
        const int x,
        const int y,
        const int halo_depth,
        const int depth,
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
        double* buffer,
        bool is_offload)
#else    //openmp and mpi
        double* buffer)
#endif
{
#ifdef SPEC_OPENACC
#pragma acc kernels loop independent collapse(2) async(0) \
    if(is_offload) present(buffer[:x*y])
#endif
#ifdef SPEC_OPENMP_TARGET
#pragma omp target teams distribute parallel for collapse(2) map(tofrom: buffer[:x*y])
#else
#ifdef SPEC_OPENMP
#pragma omp parallel for 
#endif
#endif
    for(int jj = halo_depth; jj < y-halo_depth; ++jj)
    {
        for(int kk = 0; kk < depth; ++kk)
        {
            int base = jj*x;
            buffer[base+(x-halo_depth+kk)] 
                = buffer[base+(x-halo_depth-1-kk)];
        }
    }
}

// Update top halo.
void update_top(
        const int x,
        const int y,
        const int halo_depth,
        const int depth, 
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
        double* buffer,
        bool is_offload)
#else    //openmp and mpi
        double* buffer)
#endif
{
#ifdef SPEC_OPENACC
#pragma acc kernels loop independent collapse(2) async(0)\
  if(is_offload) present(buffer[:x*y])
#endif
#ifdef SPEC_OPENMP_TARGET
#pragma omp target teams distribute parallel for map(tofrom: buffer[:x*y]) collapse(2)
#else
#ifdef SPEC_OPENMP
#pragma omp parallel for collapse(2)
#endif
#endif
    for(int jj = 0; jj < depth; ++jj)
    {
        for(int kk = halo_depth; kk < x-halo_depth; ++kk)
        {
            int base = kk;
            buffer[base+(y-halo_depth+jj)*x] 
                = buffer[base+(y-halo_depth-1-jj)*x];
        }
    }
}

// Updates bottom halo.
void update_bottom(
        const int x,
        const int y,
        const int halo_depth,
        const int depth, 
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
        double* buffer,
        bool is_offload)
#else    //openmp and mpi
        double* buffer)
#endif
{
//printf("junjie debug %d\n", is_offload);
#ifdef SPEC_OPENACC
 #pragma acc kernels loop independent collapse(2) async(0)\
  if(is_offload) present(buffer[:x*y])
#endif
#ifdef SPEC_OPENMP_TARGET
 #pragma omp target teams distribute parallel for map(tofrom: buffer[:x*y]) collapse(2)
#else
 #ifdef SPEC_OPENMP
   #pragma omp parallel for collapse(2)
 #endif
#endif
    for(int jj = 0; jj < depth; ++jj)
    {
        for(int kk = halo_depth; kk < x-halo_depth; ++kk)
        {
            int base = kk;
            buffer[base+(halo_depth-jj-1)*x] 
                = buffer[base+(halo_depth+jj)*x];
        }
    }
}

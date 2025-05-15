#include <stdlib.h>
#include "../shared.h"

/*
 *		SHARED SOLVER METHODS
 */

// Copies the current u into u0
void copy_u(
        const int x,
        const int y,
        const int halo_depth,
        double* u0,
        double* u)
{
#ifdef SPEC_OPENACC
#pragma acc kernels loop independent async(0) collapse(2) \
    present(u0[:x*y], u[:x*y])
#endif
#ifdef SPEC_OPENMP_TARGET
#pragma omp target teams distribute parallel for collapse(2) map(tofrom: u0[:x*y], u[:x*y])
#else
#ifdef SPEC_OPENMP
#pragma omp parallel for
#endif
#endif
    for(int jj = halo_depth; jj < y-halo_depth; ++jj)
    {
        for(int kk = halo_depth; kk < x-halo_depth; ++kk)
        {
            const int index = kk + jj*x;
            u0[index] = u[index];	
        }
    }
}

// Calculates the current value of r
void calculate_residual(
        const int x,
        const int y,
        const int halo_depth,
        double* u,
        double* u0,
        double* r,
        double* kx,
        double* ky)
{
#ifdef SPEC_OPENACC
#pragma acc kernels loop independent async(0) collapse(2) \
    present(kx[:x*y], ky[:x*y], u[:x*y], u0[:x*y], r[:x*y])
#endif
#ifdef SPEC_OPENMP_TARGET
#pragma omp target teams distribute parallel for  collapse(2) map(tofrom: kx[:x*y], ky[:x*y], u[:x*y], u0[:x*y], r[:x*y])
#else
#ifdef SPEC_OPENMP
#pragma omp parallel for
#endif
#endif
    for(int jj = halo_depth; jj < y-halo_depth; ++jj)
    {
        for(int kk = halo_depth; kk < x-halo_depth; ++kk)
        {
            const int index = kk + jj*x;
            const double smvp = SMVP(u);
            r[index] = u0[index] - smvp;
        }
    }
}

// Calculates the 2 norm of a given buffer
void calculate_2norm(
        const int x,
        const int y,
        const int halo_depth,
        double* buffer,
        double* norm)
{
    double norm_temp = 0.0;

#ifdef SPEC_OPENACC
#pragma acc kernels loop independent async(0) collapse(2) \
    present(buffer[:x*y])
#endif
#ifdef SPEC_OPENMP_TARGET
#pragma omp target teams distribute parallel for  reduction(+:norm_temp) collapse(2) map(tofrom: buffer[:x*y])
#else
#ifdef SPEC_OPENMP
#pragma omp parallel for reduction(+:norm_temp)
#endif
#endif
    for(int jj = halo_depth; jj < y-halo_depth; ++jj)
    {
        for(int kk = halo_depth; kk < x-halo_depth; ++kk)
        {
            const int index = kk + jj*x;
            norm_temp += buffer[index]*buffer[index];			
        }
    }

#ifdef SPEC_OPENACC
#pragma acc wait
#endif
    *norm += norm_temp;
}

// Finalises the solution
void finalise(
        const int x,
        const int y,
        const int halo_depth,
        double* energy,
        double* density,
        double* u)
{
#ifdef SPEC_OPENACC
#pragma acc kernels loop independent async(0) collapse(2) \
    present(energy[:x*y], u[:x*y], density[:x*y])
#endif
#ifdef SPEC_OPENMP_TARGET
#pragma omp target teams distribute parallel for collapse(2) map(tofrom: energy[:x*y], u[:x*y], density[:x*y])
#else
#ifdef SPEC_OPENMP
#pragma omp parallel for
#endif
#endif
    for(int jj = halo_depth; jj < y-halo_depth; ++jj)
    {
        for(int kk = halo_depth; kk < x-halo_depth; ++kk)
        {
            const int index = kk + jj*x;
            energy[index] = u[index]/density[index];
        }
    }
}

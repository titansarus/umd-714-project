#include "../shared.h"

/*
 *		PPCG SOLVER KERNEL
 */

// Initialises the PPCG solver
void ppcg_init(
        const int x,
        const int y,
        const int halo_depth,
        double theta,
        double* r,
        double* sd)
{
#ifdef SPEC_OPENACC
#pragma acc kernels loop independent collapse(2) \
    present(sd[:x*y], r[:x*y])
#endif
#ifdef SPEC_OPENMP_TARGET
#pragma omp target teams distribute parallel for collapse(2) map(tofrom: sd[:x*y], r[:x*y])
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
            sd[index] = r[index] / theta;
        }
    }
}

// The PPCG inner iteration
void ppcg_inner_iteration(
        const int x,
        const int y,
        const int halo_depth,
        double alpha,
        double beta,
        double* u,
        double* r,
        double* kx,
        double* ky,
        double* sd)
{
#ifdef SPEC_OPENACC
#pragma acc kernels loop independent collapse(2) \
    present(r[:x*y], u[:x*y], kx[:x*y], ky[:x*y], sd[:x*y])
#endif
#ifdef SPEC_OPENMP_TARGET
#pragma omp target teams distribute parallel for collapse(2) map(tofrom: r[:x*y], u[:x*y], kx[:x*y], ky[:x*y], sd[:x*y])
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
            const double smvp = SMVP(sd);
            r[index] -= smvp;
            u[index] += sd[index];
        }
    }

#ifdef SPEC_OPENACC
#pragma acc kernels loop independent collapse(2) \
    present(sd[:x*y], r[:x*y])
#endif
#ifdef SPEC_OPENMP_TARGET
#pragma omp target teams distribute parallel for collapse(2) map(tofrom: sd[:x*y], r[:x*y])
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
            sd[index] = alpha*sd[index] + beta*r[index];
        }
    }
}


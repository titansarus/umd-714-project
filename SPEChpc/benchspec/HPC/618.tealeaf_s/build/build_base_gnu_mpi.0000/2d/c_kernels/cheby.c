#include "../shared.h"

/*
 *		CHEBYSHEV SOLVER KERNEL
 */

// Calculates the new value for u.
void cheby_calc_u(
        const int x,
        const int y,
        const int halo_depth,
        double* u,
        double* p)
{
#ifdef SPEC_OPENACC
    int xy = x*y;

#pragma acc kernels loop independent collapse(2) \
    present(u[:xy], p[:xy])
#endif
#ifdef SPEC_OPENMP_TARGET
    int xy = x*y;
#pragma omp target teams distribute parallel for collapse(2) map(tofrom: u[:xy], p[:xy])
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
            u[index] += p[index];
        }
    }
}

// Initialises the Chebyshev solver
void cheby_init(
        const int x,
        const int y,
        const int halo_depth,
        const double theta,
        double* u,
        double* u0,
        double* p,
        double* r,
        double* w,
        double* kx,
        double* ky)
{
#ifdef SPEC_OPENACC
    int xy = x*y;

#pragma acc kernels loop independent collapse(2) \
    present(w[:xy], r[:xy], kx[:xy], ky[:xy], p[:xy], u0[:xy], u[:xy])
#endif
#ifdef SPEC_OPENMP_TARGET
    int xy = x*y;
#pragma omp target teams distribute parallel for collapse(2) map(tofrom: w[:xy], r[:xy], kx[:xy], ky[:xy], p[:xy], u0[:xy], u[:xy])
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
            w[index] = smvp;
            r[index] = u0[index]-w[index];
            p[index] = r[index] / theta;
        }
    }

    cheby_calc_u(x, y, halo_depth, u, p);
}

// The main chebyshev iteration
void cheby_iterate(
        const int x,
        const int y,
        const int halo_depth,
        double alpha,
        double beta,
        double* u,
        double* u0,
        double* p,
        double* r,
        double* w,
        double* kx,
        double* ky)
{
#ifdef SPEC_OPENACC
    int xy = x*y;

#pragma acc kernels loop independent collapse(2) \
    present(w[:xy], r[:xy], kx[:xy], ky[:xy], p[:xy], u0[:xy], u[:xy])
#endif
#ifdef SPEC_OPENMP_TARGET
    int xy = x*y;
#pragma omp target teams distribute parallel for collapse(2) map(tofrom: w[:xy], r[:xy], kx[:xy], ky[:xy], p[:xy], u0[:xy], u[:xy])
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
            w[index] = smvp;
            r[index] = u0[index]-w[index];
            p[index] = alpha*p[index] + beta*r[index];
        }
    }

    cheby_calc_u(x, y, halo_depth, u, p);
}


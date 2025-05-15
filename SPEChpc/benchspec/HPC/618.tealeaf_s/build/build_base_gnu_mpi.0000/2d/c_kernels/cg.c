#include <stdlib.h>
#include "../shared.h"

/*
 *		CONJUGATE GRADIENT SOLVER KERNEL
 */

// Initialises the CG solver
void cg_init(
        const int x,
        const int y,
        const int halo_depth,
        const int coefficient,
        double rx,
        double ry,
        double* rro,
        double* density,
        double* energy,
        double* u,
        double* p,
        double* r,
        double* w,
        double* kx,
        double* ky)
{
    if(coefficient != CONDUCTIVITY && coefficient != RECIP_CONDUCTIVITY)
    {
        die(__LINE__, __FILE__, "Coefficient %d is not valid.\n", coefficient);
    }
#ifdef SPEC_OPENACC
    int xy = x*y;

#pragma acc kernels loop independent collapse(2) async(0) \
    present(p[:xy], r[:xy], u[:xy], energy[:xy], density[:xy])
#endif
#ifdef SPEC_OPENMP_TARGET
    int xy = x*y;
 #ifdef SPEC_USE_INNER_SIMD
 #pragma omp target teams distribute parallel for map(tofrom: p[:xy], r[:xy], u[:xy], energy[:xy], density[:xy])
 #else
 #pragma omp target teams distribute parallel for simd collapse(2) map(tofrom: p[:xy], r[:xy], u[:xy], energy[:xy], density[:xy])
 #endif
#else
 #ifdef SPEC_OPENMP
 #pragma omp parallel for
 #endif
#endif
    for(int jj = 0; jj < y; ++jj)
    {
#ifdef SPEC_USE_INNER_SIMD
	#pragma omp simd
#endif
        for(int kk = 0; kk < x; ++kk)
        {
            const int index = kk + jj*x;
            p[index] = 0.0;
            r[index] = 0.0;
            u[index] = energy[index]*density[index];
        }
    }

#ifdef SPEC_OPENACC
#pragma acc kernels loop independent collapse(2) async(0) \
    present(w[:xy], density[:xy])
#endif
#ifdef SPEC_OPENMP_TARGET
#ifdef SPEC_USE_INNER_SIMD
#pragma omp target teams distribute parallel for map(tofrom: w[:xy], density[:xy])
#else
#pragma omp target teams distribute parallel for simd collapse(2) map(tofrom: w[:xy], density[:xy])
#endif
#else
#ifdef SPEC_OPENMP
#pragma omp parallel for
#endif
#endif
    for(int jj = 1; jj < y-1; ++jj)
    {
#ifdef SPEC_USE_INNER_SIMD
	#pragma omp simd
#endif
        for(int kk = 1; kk < x-1; ++kk)
        {
            const int index = kk + jj*x;
            w[index] = (coefficient == CONDUCTIVITY) 
                ? density[index] : 1.0/density[index];
        }
    }

#ifdef SPEC_OPENACC
#pragma acc kernels loop independent collapse(2) async(0) \
    present(kx[:xy], ky[:xy], w[:xy])
#endif
#ifdef SPEC_OPENMP_TARGET
#ifdef SPEC_USE_INNER_SIMD
#pragma omp target teams distribute parallel for map(kx[:xy], ky[:xy], w[:xy])
#else
#pragma omp target teams distribute parallel for simd collapse(2) map(kx[:xy], ky[:xy], w[:xy])
#endif
#else
#ifdef SPEC_OPENMP
#pragma omp parallel for
#endif
#endif
    for(int jj = halo_depth; jj < y-1; ++jj)
    {
#ifdef SPEC_USE_INNER_SIMD
	#pragma omp simd
#endif
        for(int kk = halo_depth; kk < x-1; ++kk)
        {
            const int index = kk + jj*x;
            kx[index] = rx*(w[index-1]+w[index]) /
                (2.0*w[index-1]*w[index]);
            ky[index] = ry*(w[index-x]+w[index]) /
                (2.0*w[index-x]*w[index]);
        }
    }

    double rro_temp = 0.0;
#ifdef SPEC_OPENACC
#pragma acc kernels loop independent collapse(2) async(0) \
    present(w[:xy], u[:xy], r[:xy], p[:xy], kx[:xy], ky[:xy])
#endif
#ifdef SPEC_OPENMP_TARGET
#ifdef SPEC_USE_INNER_SIMD
	#pragma omp target teams distribute parallel for map(tofrom: w[:xy], u[:xy], r[:xy], p[:xy], kx[:xy], ky[:xy], rro[:1]) reduction(+:rro_temp)
#else
	#pragma omp target teams distribute parallel for simd collapse(2) map(tofrom: w[:xy], u[:xy], r[:xy], p[:xy], kx[:xy], ky[:xy], rro[:1]) reduction(+:rro_temp)
#endif
#else
#ifdef SPEC_OPENMP
#pragma omp parallel for reduction(+:rro_temp)
#endif
#endif
    for(int jj = halo_depth; jj < y-halo_depth; ++jj)
    {
#ifdef SPEC_USE_INNER_SIMD
	#pragma omp simd
#endif
        for(int kk = halo_depth; kk < x-halo_depth; ++kk)
        {
            const int index = kk + jj*x;
            const double smvp = SMVP(u);
            w[index] = smvp;
            r[index] = u[index]-w[index];
            p[index] = r[index];
            rro_temp += r[index]*p[index];
        }
    }

#ifdef SPEC_OPENACC
#pragma acc wait
#endif

    // Sum locally
    *rro += rro_temp;
}

// Calculates w
void cg_calc_w(
        const int x,
        const int y,
        const int halo_depth,
        double* pw,
        double* p,
        double* w,
        double* kx,
        double* ky)
{
    double pw_temp = 0.0;

#ifdef SPEC_OPENACC
    int xy = x*y;

#pragma acc kernels loop independent collapse(2) async(0) \
    present(w[:xy], p[:xy], kx[:xy], ky[:xy])
#endif
#ifdef SPEC_OPENMP_TARGET
    int xy = x*y;
#ifdef SPEC_USE_INNER_SIMD
	#pragma omp target teams distribute parallel for map(tofrom:pw_temp, w[:xy], p[:xy], kx[:xy], ky[:xy]) reduction(+:pw_temp)
#else
	#pragma omp target teams distribute parallel for simd map(tofrom:pw_temp, w[:xy], p[:xy], kx[:xy], ky[:xy]) reduction(+:pw_temp) collapse(2)
#endif
#else
#ifdef SPEC_OPENMP
#pragma omp parallel for reduction(+:pw_temp)
#endif
#endif
    for(int jj = halo_depth; jj < y-halo_depth; ++jj)
    {
#ifdef SPEC_USE_INNER_SIMD
	#pragma omp simd
#endif
        for(int kk = halo_depth; kk < x-halo_depth; ++kk)
        {
            const int index = kk + jj*x;
            const double smvp = SMVP(p);
            w[index] = smvp;
            pw_temp += w[index]*p[index];
        }
    }

#ifdef SPEC_OPENACC
#pragma acc wait
#endif
    *pw += pw_temp;
}

// Calculates u and r
void cg_calc_ur(
        const int x,
        const int y,
        const int halo_depth,
        const double alpha,
        double* rrn,
        double* u,
        double* p,
        double* r,
        double* w)
{
    double rrn_temp = 0.0;

#ifdef SPEC_OPENACC
    int xy = x*y;

#pragma acc kernels loop independent async(0) \
    present(u[:xy], p[:xy], w[:xy], r[:xy])
#endif
#ifdef SPEC_OPENMP_TARGET
    int xy = x*y;
#ifdef SPEC_USE_INNER_SIMD
	#pragma omp target teams distribute parallel for map(tofrom:rrn_temp) reduction(+:rrn_temp) map(tofrom: u[:xy], p[:xy], w[:xy], r[:xy])
#else
	#pragma omp target teams distribute parallel for simd map(tofrom:rrn_temp) reduction(+:rrn_temp) collapse(2) map(tofrom: u[:xy], p[:xy], w[:xy], r[:xy])
#endif
#else
#ifdef SPEC_OPENMP
#pragma omp parallel for reduction(+:rrn_temp)
#endif
#endif
    for(int jj = halo_depth; jj < y-halo_depth; ++jj)
    {
#ifdef SPEC_OPENACC
#pragma acc loop independent
#endif
#ifdef SPEC_USE_INNER_SIMD
	#pragma omp simd
#endif
        for(int kk = halo_depth; kk < x-halo_depth; ++kk)
        {
            const int index = kk + jj*x;

            u[index] += alpha*p[index];
            r[index] -= alpha*w[index];
            rrn_temp += r[index]*r[index];
        }
    }

#ifdef SPEC_OPENACC
#pragma acc wait
#endif
    *rrn += rrn_temp;
}

// Calculates p
void cg_calc_p(
        const int x,
        const int y,
        const int halo_depth,
        const double beta,
        double* p,
        double* r)
{
#ifdef SPEC_OPENACC
    int xy = x*y;

#pragma acc kernels loop independent async(0) \
    present(p[:xy], r[:xy])
#endif
#ifdef SPEC_OPENMP_TARGET
    int xy = x*y;
#ifdef SPEC_USE_INNER_SIMD
	#pragma omp target teams distribute parallel for map(tofrom: p[:xy], r[:xy])
#else
	#pragma omp target teams distribute parallel for simd collapse(2) map(tofrom: p[:xy], r[:xy])
#endif
#else
#ifdef SPEC_OPENMP
#pragma omp parallel for
#endif
#endif
    for(int jj = halo_depth; jj < y-halo_depth; ++jj)
    {
#ifdef SPEC_OPENACC
#pragma acc loop independent
#endif
#ifdef SPEC_USE_INNER_SIMD
	#pragma omp simd
#endif
        for(int kk = halo_depth; kk < x-halo_depth; ++kk)
        {
            const int index = kk + jj*x;

            p[index] = beta*p[index] + r[index];
        }
    }
}

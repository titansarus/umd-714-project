#pragma once

#include <cmath>
#include <vector>
#include "utils.hpp"

#include "kernels.hpp"

namespace sphexa
{
namespace sph
{
template <typename T, class Dataset>
void computeDensityImpl(const std::vector<int> &l, Dataset &d)
{
    const size_t n = l.size();
    const size_t ngmax = d.ngmax;
    const int *clist = l.data();

    // I assume here that indexes to compute in l are sequentially increasing, e.g 0, 1, 2, 3
    const size_t neighborsOffset = l.front() * ngmax;
    const int *neighbors = d.neighbors.data() + neighborsOffset;

    const size_t nOffset = l.front();
    const int *neighborsCount = d.neighborsCount.data() + nOffset;

    const T *h = d.h.data();
    const T *m = d.m.data();
    const T *x = d.x.data();
    const T *y = d.y.data();
    const T *z = d.z.data();

    const T *whLt = d.wharmonicLookupTable;
    const T *whDerLt = d.wharmonicDerivativeLookupTable;
    const size_t whSize = d.wharmonicLookupTableSize;
    T *ro = d.ro.data() + nOffset;

    const BBox<T> bbox = d.bbox;

    const T K = d.K;
    const T sincIndex = d.sincIndex;

#ifdef SPEC_OPENMP_TARGET
    const int np = d.x.size();
    const int64_t allNeighbors = n * ngmax;
    // Apparently Cray with -O2 has a bug when calling target regions in a loop. (and computeDensityImpl can be called in a loop).
    // A workaround is to call some method or allocate memory to either prevent buggy optimization or other side effect.
    // with -O1 there is no problem
    // Tested with Cray 8.7.3 with NVIDIA Tesla P100 on PizDaint
    std::vector<T> imHereBecauseOfCrayCompilerO2Bug(4, 10);

    #ifdef SPEC_USE_LT_IN_KERNELS
    #pragma omp target teams distribute parallel for map(to                                                                                \
                             : clist [0:n], neighbors [0:allNeighbors], neighborsCount [0:n], m [0:np], h [0:np], x [0:np], y [0:np],      \
                             whLt [0:whSize], whDerLt [0:whSize],                                                                          \
                             z [0:np]) map(from                                                                                            \
                                             : ro [0:n])
    #else
    #pragma omp target teams distribute parallel for map(to                                                                                \
                             : clist [0:n], neighbors [0:allNeighbors], neighborsCount [0:n], m [0:np], h [0:np], x [0:np], y [0:np],      \
                             z [0:np]) map(from                                                                                            \
                                             : ro [0:n])
    #endif
#elif defined(SPEC_OPENACC)
    const int np = d.x.size();
    const int64_t allNeighbors = n * ngmax;
    #ifdef SPEC_USE_LT_IN_KERNELS
    #pragma acc parallel loop copyin(clist [0:n], neighbors [0:allNeighbors], neighborsCount [0:n], m [0:np], h [0:np], x [0:np], y [0:np], \
                            whLt [0:whSize], whDerLt [0:whSize],                                                                            \
                            z[0:np], bbox) copyout(ro [0:n]) default(present)
    #else
    #pragma acc parallel loop copyin(clist [0:n], neighbors [0:allNeighbors], neighborsCount [0:n], m [0:np], h [0:np], x [0:np], y [0:np], \
                            z[0:np], bbox) copyout(ro [0:n]) default(present)
    #endif
#else
#ifdef SPEC_OPENMP
#pragma omp parallel for
#endif
#endif
    for (size_t pi = 0; pi < n; pi++)
    {
        const int i = clist[pi];
        const int nn = neighborsCount[pi];

        T roloc = 0.0;

        // int converstion to avoid a bug that prevents vectorization with some compilers
        for (int pj = 0; pj < nn; pj++)
        {
            const int j = neighbors[pi * ngmax + pj];

            // later can be stores into an array per particle
            T dist = distancePBC(bbox, h[i], x[i], y[i], z[i], x[j], y[j], z[j]); // store the distance from each neighbor

            // calculate the v as ratio between the distance and the smoothing length
            T vloc = dist / h[i];

#ifndef NDEBUG
            if (vloc > 2.0 + 1e-6 || vloc < 0.0)
                printf("ERROR:Density(%d,%d) vloc %f -- x %f %f %f -- %f %f %f -- dist %f -- hi %f\n", i, j, vloc, x[i], y[i], z[i], x[j],
                       y[j], z[j], dist, h[i]);
#endif

	    #if defined(SPEC_USE_LT_IN_KERNELS)
            const T w = K * math_namespace::pow(wharmonic(vloc, whSize, whLt, whDerLt), (int)sincIndex);
	    #else
	    const T w = K * math_namespace::pow(wharmonic(vloc), (int)sincIndex);
	    #endif
	    
            const T value = w / (h[i] * h[i] * h[i]);
            roloc += value * m[j];
        }

        ro[pi] = roloc + m[i] * K / (h[i] * h[i] * h[i])
;
#ifndef NDEBUG
        if (std::isnan(ro[i])) printf("ERROR::Density(%d) density %f, position: (%f %f %f), h: %f\n", i, ro[i], x[i], y[i], z[i], h[i]);
#endif
    }

}

template <typename T, class Dataset>
void computeDensity(const std::vector<int> &l, Dataset &d)
{
#if defined(SPEC_OPENMP_TARGET) || defined(SPEC_OPENACC)
    const T *h = d.h.data();
    const T *m = d.m.data();
    const T *x = d.x.data();
    const T *y = d.y.data();
    const T *z = d.z.data();

    const T *whLt = d.wharmonicLookupTable;
    const T *whDerLt = d.wharmonicDerivativeLookupTable;
    const size_t whSize = d.wharmonicLookupTableSize;
    const int np = d.x.size();

# ifdef SPEC_OPENMP_TARGET
    #ifdef SPEC_USE_LT_IN_KERNELS
    #pragma omp target data map(to: m [0:np], h [0:np], x [0:np], y [0:np],    \
                                    whLt [0:whSize], whDerLt [0:whSize],       \
                                    z [0:np])
    #else
    #pragma omp target data map(to: m [0:np], h [0:np], x [0:np], y [0:np],    \
                                    z [0:np])
    #endif
# else
    #ifdef SPEC_USE_LT_IN_KERNELS
    #pragma acc data copyin(m[0:np], h[0:np], x[0:np], y[0:np],    \
                            whLt[0:whSize], whDerLt[0:whSize],     \
                            z[0:np])
    #else
    #pragma acc data copyin(m[0:np], h[0:np], x[0:np], y[0:np],    \
                            z[0:np])
    #endif
# endif
#endif // SPEC_OPENMP_TARGET && SPEC_OPENACC
    {
        for (const auto &clist : utils::partition(l, d.noOfGpuLoopSplits))
        {
            computeDensityImpl<T>(clist, d);
        }
    }

    // for (size_t i = 0; i < d.ro.size(); ++i)
    // {
    //     printf(" %lu: %.15f", i, d.ro[i]);
    //     if (i % 10 == 0) printf("\n");
    // }
}

template <typename T, class Dataset>
void initFluidDensityAtRest(const std::vector<int> &clist, Dataset &d)
{
    const T *ro = d.ro.data();
    T *ro_0 = d.ro_0.data();

#if defined(SPEC_OPENMP) || defined(SPEC_OPENMP_TARGET)
#pragma omp parallel for
#endif
    for (size_t pi = 0; pi < clist.size(); ++pi)
    {
        const int i = clist[pi];
        ro_0[i] = ro[i];
    }
}

} // namespace sph
} // namespace sphexa

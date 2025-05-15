#include "../shared.h"

// Store original energy state
void store_energy(
        int x,
        int y,
        double* energy0,
        double* energy)
{
#if defined(SPEC_OPENMP) || defined(SPEC_OPENMP_TARGET)
#pragma omp parallel for
#endif
    for(int ii = 0; ii < x*y; ++ii)
    {
        energy[ii] = energy0[ii];
    }
}


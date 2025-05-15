#include "allocator.h"

Allocator* global_allocator;

void free_allocator(Allocator* a)
{
    free(a->all_Monomer.buf);
    free(a->all_Monomer_msd.buf);
    free(a->all_RNG_STATE.buf);
    free(a->all_MERSENNE_TWISTER_STATE.buf);
    free(a->all_MTTSTATE.buf);
    free(a->all_uint_t.buf);
#if defined(SPEC_OPENMP_TARGET) && !defined(SPEC_USE_HOST_THREADS)
    FREE_TYPE_BUF_DEVICE(Monomer)
    FREE_TYPE_BUF_DEVICE(Monomer, _msd)
    FREE_TYPE_BUF_DEVICE(RNG_STATE)
    FREE_TYPE_BUF_DEVICE(MERSENNE_TWISTER_STATE)
    FREE_TYPE_BUF_DEVICE(MTTSTATE)
    FREE_TYPE_BUF_DEVICE(uint_t)
#endif
}

void init_allocator(Allocator* a)
{
    INIT_TYPE_BUF(Monomer, p);
    INIT_TYPE_BUF(Monomer, p, _msd);
    INIT_TYPE_BUF(RNG_STATE, p);
    INIT_TYPE_BUF(MERSENNE_TWISTER_STATE, p);
    INIT_TYPE_BUF(MTTSTATE, p);
    INIT_TYPE_BUF(uint_t, p);
}

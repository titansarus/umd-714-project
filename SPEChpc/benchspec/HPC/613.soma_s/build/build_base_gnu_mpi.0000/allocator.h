#pragma once

#include "phase.h"
#include "device.h"

#ifdef SPEC_OPENMP_TARGET
#include <omp.h>
#endif

typedef unsigned int uint_t;

typedef struct Allocator {
#define DECL_TYPE_BUF(TYPE, ...) \
    struct { \
        size_t capacity, size; \
        TYPE* buf, *device_buf; \
    } all_ ## TYPE ## __VA_ARGS__ ;
#define FUNCT_TYPE_ALLOC(TYPE, ...) \
    static inline size_t alloc_ ## TYPE ## __VA_ARGS__ (size_t n) \
    { \
        global_allocator->all_ ## TYPE ## __VA_ARGS__ .size += n; \
        assert (global_allocator->all_ ## TYPE ## __VA_ARGS__ .size <= global_allocator->all_ ## TYPE ## __VA_ARGS__ .capacity); \
        return global_allocator->all_ ## TYPE ## __VA_ARGS__ .size - n; \
    } \
    static inline TYPE * alloc_ ## TYPE ## __VA_ARGS__ ## _ptr (size_t n) \
    { \
        return global_allocator->all_ ## TYPE ## __VA_ARGS__ .buf + alloc_ ## TYPE ## __VA_ARGS__(n); \
    }
#define INIT_TYPE_BUF(TYPE, p, ...) \
    global_allocator->all_ ## TYPE ## __VA_ARGS__ .capacity = 0; \
    global_allocator->all_ ## TYPE ## __VA_ARGS__ .size = 0; \
    global_allocator->all_ ## TYPE ## __VA_ARGS__ .buf = 0; \
    global_allocator->all_ ## TYPE ## __VA_ARGS__ .device_buf = 0;
#define SET_TYPE_BUF_HOST(TYPE, n, ...) \
    global_allocator->all_ ## TYPE ## __VA_ARGS__ .capacity = n; \
    global_allocator->all_ ## TYPE ## __VA_ARGS__ .size = 0; \
    if(global_allocator->all_ ## TYPE ## __VA_ARGS__ .buf) \
        free(global_allocator->all_ ## TYPE ## __VA_ARGS__ .buf); \
    global_allocator->all_ ## TYPE ## __VA_ARGS__ .buf = \
        global_allocator->all_ ## TYPE ## __VA_ARGS__ .device_buf = \
        n ? malloc(sizeof(TYPE) * global_allocator->all_ ## TYPE ## __VA_ARGS__ .capacity) : 0;

#if defined(SPEC_OPENMP_TARGET) && !defined(SPEC_USE_HOST_THREADS)
#define FREE_TYPE_BUF_DEVICE(TYPE, ...) \
    if(global_allocator->all_ ## TYPE ## __VA_ARGS__ .device_buf != global_allocator->all_ ## TYPE ## __VA_ARGS__ .buf ) \
    { \
        /* omp_target_associate_ptr is broken regarding initial_device, using memcpy instead */ \
        omp_target_disassociate_ptr(global_allocator->all_ ## TYPE ## __VA_ARGS__ .buf, soma_get_device()); \
        omp_target_free(global_allocator->all_ ## TYPE ## __VA_ARGS__ .device_buf, soma_get_device()); \
    }
#define SET_TYPE_BUF(TYPE, n, ...) \
    FREE_TYPE_BUF_DEVICE(TYPE, __VA_ARGS__) \
    SET_TYPE_BUF_HOST(TYPE, n, __VA_ARGS__) \
    if ( n > 0 ) \
    { \
        global_allocator->all_ ## TYPE ## __VA_ARGS__ .device_buf = \
            omp_target_alloc(sizeof(TYPE) * global_allocator->all_ ## TYPE ## __VA_ARGS__ .capacity \
                , soma_get_device()); \
        omp_target_associate_ptr(global_allocator->all_ ## TYPE ## __VA_ARGS__ .buf \
            , global_allocator->all_ ## TYPE ## __VA_ARGS__ .device_buf \
            , global_allocator->all_ ## TYPE ## __VA_ARGS__ .capacity * sizeof(TYPE) \
            , 0, soma_get_device() \
            ); \
    }
#else
#define SET_TYPE_BUF(TYPE, n, ...) SET_TYPE_BUF_HOST(TYPE, n, __VA_ARGS__)
#endif

    DECL_TYPE_BUF(Monomer)
    DECL_TYPE_BUF(Monomer, _msd)
    DECL_TYPE_BUF(RNG_STATE)
    DECL_TYPE_BUF(MERSENNE_TWISTER_STATE)
    DECL_TYPE_BUF(MTTSTATE)
    DECL_TYPE_BUF(uint_t)
} Allocator;

extern Allocator* global_allocator;

void init_allocator(Allocator* a);
void free_allocator(Allocator* a);

FUNCT_TYPE_ALLOC(Monomer)
FUNCT_TYPE_ALLOC(Monomer, _msd)
FUNCT_TYPE_ALLOC(RNG_STATE)
FUNCT_TYPE_ALLOC(MERSENNE_TWISTER_STATE)
FUNCT_TYPE_ALLOC(MTTSTATE)
FUNCT_TYPE_ALLOC(uint_t)

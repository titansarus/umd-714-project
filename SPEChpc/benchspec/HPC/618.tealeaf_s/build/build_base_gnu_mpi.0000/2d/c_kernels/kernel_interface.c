#include "../kernel_interface.h"
#include "c_kernels.h"
#include <stdlib.h>

// Initialisation kernels
void run_set_chunk_data(Chunk* chunk, Settings* settings)
{
    set_chunk_data(
            settings, chunk->x, chunk->y, chunk->left, chunk->bottom, 
            chunk->cell_x, chunk->cell_y, chunk->vertex_x, chunk->vertex_y,
            chunk->volume, chunk->x_area, chunk->y_area);
}

void run_set_chunk_state(Chunk* chunk, Settings* settings, State* states)
{
    set_chunk_state(chunk->x, chunk->y, chunk->vertex_x, chunk->vertex_y, 
            chunk->cell_x, chunk->cell_y, chunk->density, chunk->energy0, 
            chunk->u, settings->num_states, states);
}

void run_kernel_initialise(Chunk* chunk, Settings* settings)
{
    kernel_initialise(settings, chunk->x, chunk->y, &(chunk->density0), 
            &(chunk->density), &(chunk->energy0), &(chunk->energy), 
            &(chunk->u), &(chunk->u0), &(chunk->p), &(chunk->r), 
            &(chunk->mi), &(chunk->w), &(chunk->kx), &(chunk->ky), 
            &(chunk->sd), &(chunk->volume), 
            &(chunk->x_area), &(chunk->y_area), &(chunk->cell_x), 
            &(chunk->cell_y), &(chunk->cell_dx), &(chunk->cell_dy),
            &(chunk->vertex_dx), &(chunk->vertex_dy), &(chunk->vertex_x), 
            &(chunk->vertex_y), &(chunk->cg_alphas), &(chunk->cg_betas), 
            &(chunk->cheby_alphas), &(chunk->cheby_betas));
}

void run_kernel_finalise(
        Chunk* chunk, Settings* settings)
{
    kernel_finalise(
            chunk->density0, chunk->density, chunk->energy0, chunk->energy,
            chunk->u, chunk->u0, chunk->p, chunk->r, chunk->mi, chunk->w,
            chunk->kx, chunk->ky, chunk->sd, chunk->volume, chunk->x_area,
            chunk->y_area, chunk->cell_x, chunk->cell_y, chunk->cell_dx,
            chunk->cell_dy, chunk->vertex_dx, chunk->vertex_dy, chunk->vertex_x,
            chunk->vertex_y, chunk->cg_alphas, chunk->cg_betas,
            chunk->cheby_alphas, chunk->cheby_betas);
}

// Solver-wide kernels
void run_local_halos(
        Chunk* chunk, Settings* settings, int depth)
{
#ifndef SPEC
    START_PROFILING(settings->kernel_profile);
#endif
    local_halos(chunk->x, chunk->y, depth, settings->halo_depth, 
            chunk->neighbours, settings->fields_to_exchange, chunk->density,
            chunk->energy0, chunk->energy, chunk->u, chunk->p, 
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
            chunk->sd, settings->is_offload);
#else    //openmp and mpi
            chunk->sd);
#endif
#ifndef SPEC
    STOP_PROFILING(settings->kernel_profile, __func__);
#endif
}

void run_pack_or_unpack(
        Chunk* chunk, Settings* settings, int depth,
        int face, bool pack, double* field, double* buffer)
{
#ifndef SPEC
    START_PROFILING(settings->kernel_profile);
#endif
    pack_or_unpack(
            chunk->x, chunk->y, depth, settings->halo_depth, 
            face, pack, field, 
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
            buffer, settings->is_offload);
#else    //openmp and mpi
            buffer);
#endif
#ifndef SPEC
    STOP_PROFILING(settings->kernel_profile, __func__);
#endif
}

void run_store_energy(Chunk* chunk, Settings* settings)
{
#ifndef SPEC
    START_PROFILING(settings->kernel_profile);
#endif
    store_energy(chunk->x, chunk->y, chunk->energy0, chunk->energy);
#ifndef SPEC
    STOP_PROFILING(settings->kernel_profile, __func__);
#endif
}

void run_field_summary(
        Chunk* chunk, Settings* settings, 
        double* vol, double* mass, double* ie, double* temp)
{
#ifndef SPEC
    START_PROFILING(settings->kernel_profile);
#endif
    field_summary(chunk->x, chunk->y,
            settings->halo_depth, chunk->volume, chunk->density,
            chunk->energy0, chunk->u, vol, mass, ie, temp);
#ifndef SPEC
    STOP_PROFILING(settings->kernel_profile, __func__);
#endif
}

// CG solver kernels
void run_cg_init(
        Chunk* chunk, Settings* settings, 
        double rx, double ry, double* rro)
{
#ifndef SPEC
    START_PROFILING(settings->kernel_profile);
#endif
    cg_init(chunk->x, chunk->y, 
            settings->halo_depth, settings->coefficient, rx, ry, 
            rro, chunk->density, chunk->energy, chunk->u, 
            chunk->p, chunk->r, chunk->w, 
            chunk->kx, chunk->ky);
#ifndef SPEC
    STOP_PROFILING(settings->kernel_profile, __func__);
#endif
}

void run_cg_calc_w(Chunk* chunk, Settings* settings, double* pw)
{
#ifndef SPEC
    START_PROFILING(settings->kernel_profile);
#endif
    cg_calc_w(chunk->x, chunk->y, 
            settings->halo_depth, pw, chunk->p, 
            chunk->w, chunk->kx,
            chunk->ky);
#ifndef SPEC
    STOP_PROFILING(settings->kernel_profile, __func__);
#endif
}

void run_cg_calc_ur(
        Chunk* chunk, Settings* settings, double alpha, double* rrn)
{
#ifndef SPEC
    START_PROFILING(settings->kernel_profile);
#endif
    cg_calc_ur(chunk->x, chunk->y, 
            settings->halo_depth, alpha, rrn, chunk->u, 
            chunk->p, chunk->r, chunk->w);
#ifndef SPEC
    STOP_PROFILING(settings->kernel_profile, __func__);
#endif
}

void run_cg_calc_p(Chunk* chunk, Settings* settings, double beta)
{
#ifndef SPEC
    START_PROFILING(settings->kernel_profile);
#endif
    cg_calc_p(chunk->x, chunk->y, 
            settings->halo_depth, beta, chunk->p, 
            chunk->r);
#ifndef SPEC
    STOP_PROFILING(settings->kernel_profile, __func__);
#endif
}


// Chebyshev solver kernels
void run_cheby_init(Chunk* chunk, Settings* settings)
{
#ifndef SPEC
    START_PROFILING(settings->kernel_profile);
#endif
    cheby_init(
            chunk->x, chunk->y, settings->halo_depth, 
            chunk->theta, chunk->u, chunk->u0, 
            chunk->p, chunk->r, chunk->w, 
            chunk->kx, chunk->ky);
#ifndef SPEC
    STOP_PROFILING(settings->kernel_profile, __func__);
#endif
}

void run_cheby_iterate(
        Chunk* chunk, Settings* settings, double alpha, double beta)
{
#ifndef SPEC
    START_PROFILING(settings->kernel_profile);
#endif
    cheby_iterate(
            chunk->x, chunk->y, settings->halo_depth, alpha, beta, 
            chunk->u, chunk->u0, chunk->p, chunk->r, chunk->w, 
            chunk->kx, chunk->ky); 
#ifndef SPEC
    STOP_PROFILING(settings->kernel_profile, __func__);
#endif
}


// Jacobi solver kernels
void run_jacobi_init(
        Chunk* chunk, Settings* settings, double rx, double ry)
{
#ifndef SPEC
    START_PROFILING(settings->kernel_profile);
#endif
    jacobi_init(chunk->x, chunk->y, settings->halo_depth, 
            settings->coefficient, rx, ry, chunk->density, chunk->energy, 
            chunk->u0, chunk->u, chunk->kx, chunk->ky);
#ifndef SPEC
    STOP_PROFILING(settings->kernel_profile, __func__);
#endif
}

void run_jacobi_iterate(
        Chunk* chunk, Settings* settings, double* error)
{
#ifndef SPEC
    START_PROFILING(settings->kernel_profile);
#endif
    jacobi_iterate(
            chunk->x, chunk->y, settings->halo_depth, error, chunk->kx, 
            chunk->ky, chunk->u0, chunk->u, chunk->r);
#ifndef SPEC
    STOP_PROFILING(settings->kernel_profile, __func__);
#endif
}


// PPCG solver kernels
void run_ppcg_init(Chunk* chunk, Settings* settings)
{
#ifndef SPEC
    START_PROFILING(settings->kernel_profile);
#endif
    ppcg_init(chunk->x, chunk->y, settings->halo_depth, chunk->theta, 
            chunk->r, chunk->sd);
#ifndef SPEC
    STOP_PROFILING(settings->kernel_profile, __func__);
#endif
}

void run_ppcg_inner_iteration(
        Chunk* chunk, Settings* settings, double alpha, double beta)
{
#ifndef SPEC
    START_PROFILING(settings->kernel_profile);
#endif
    ppcg_inner_iteration(
            chunk->x, chunk->y, settings->halo_depth, alpha, beta, chunk->u, 
            chunk->r, chunk->kx, chunk->ky, chunk->sd);
#ifndef SPEC
    STOP_PROFILING(settings->kernel_profile, __func__);
#endif
}

// Shared solver kernels
void run_copy_u(Chunk* chunk, Settings* settings)
{
#ifndef SPEC
    START_PROFILING(settings->kernel_profile);
#endif
    copy_u(
            chunk->x, chunk->y, settings->halo_depth, chunk->u0, chunk->u);
#ifndef SPEC
    STOP_PROFILING(settings->kernel_profile, __func__);
#endif
}

void run_calculate_residual(Chunk* chunk, Settings* settings)
{
#ifndef SPEC
    START_PROFILING(settings->kernel_profile);
#endif
    calculate_residual(chunk->x, chunk->y, settings->halo_depth, chunk->u, 
            chunk->u0, chunk->r, chunk->kx, chunk->ky);
#ifndef SPEC
    STOP_PROFILING(settings->kernel_profile, __func__);
#endif
}

void run_calculate_2norm(
        Chunk* chunk, Settings* settings, double* buffer, double* norm)
{
#ifndef SPEC
    START_PROFILING(settings->kernel_profile);
#endif
    calculate_2norm(
            chunk->x, chunk->y, settings->halo_depth, buffer, norm);
#ifndef SPEC
    STOP_PROFILING(settings->kernel_profile, __func__);
#endif
}

void run_finalise(Chunk* chunk, Settings* settings)
{
#ifndef SPEC
    START_PROFILING(settings->kernel_profile);
#endif
    finalise(
            chunk->x, chunk->y, settings->halo_depth, chunk->energy, 
            chunk->density, chunk->u);
#ifndef SPEC
    STOP_PROFILING(settings->kernel_profile, __func__);
#endif
}


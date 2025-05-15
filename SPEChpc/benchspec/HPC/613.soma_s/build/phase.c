/* Copyright (C) 2016-2017 Ludwig Schneider

 This file is part of SOMA.

 SOMA is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SOMA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with SOMA.  If not, see <http://www.gnu.org/licenses/>.
*/

//! \file phase.c
//! \brief Implementation of phase.h


#include "phase.h"
#include <stdbool.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include "init.h"
#include "independent_sets.h"
#include "mesh.h"
#include "allocator.h"
#include "device.h"

#ifdef SPEC_OPENACC
#include <openacc.h>
#endif
#ifdef SPEC_OPENMP_TARGET
#include <omp.h>
#endif

#include "mpi_timing.h"

#define MPI_Irecv timed_MPI_Irecv
#define MPI_Isend timed_MPI_Isend
#define MPI_Wait timed_MPI_Wait
#define MPI_Reduce timed_MPI_Reduce
#define MPI_Barrier timed_MPI_Barrier
#define MPI_Allreduce timed_MPI_Allreduce
#define MPI_Allgather timed_MPI_Allgather
#define MPI_Recv timed_MPI_Recv
#define MPI_Send timed_MPI_Send

int init_phase(struct Phase * const p)
    {
    print_version(p->info_MPI.current_core);
    p->start_time = p->time;
    p->start_clock = time(NULL);
    p->n_accepts = 0;
    p->n_moves =0;
    p->n_tries_cm = 0;
    p->n_acc_cm = 0;
    p->end_mono = NULL;
    p->tps_elapsed_time = 1; //Bla default, bigger 0
    p->tps_elapsed_steps = 1; //Bla default, bigger 0

    uint64_t n_polymer_offset;
    MPI_Scan( &(p->n_polymers), &n_polymer_offset, 1,MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_MPI_Comm);
    n_polymer_offset -= p->n_polymers;

    switch(p->args.pseudo_random_number_generator_arg)
    {
    case pseudo_random_number_generator__NULL     : break;
    case pseudo_random_number_generator_arg_PCG32 : break;
    case pseudo_random_number_generator_arg_MT    : 
	SET_TYPE_BUF(MERSENNE_TWISTER_STATE, p->n_polymers);
	break;
    case pseudo_random_number_generator_arg_TT800 :
	SET_TYPE_BUF(MTTSTATE, p->n_polymers);
	break;
    }
    for(uint64_t i=0; i < p->n_polymers;i++)
	{
	p->polymers[i].set_states = NULL;
	p->polymers[i].set_permutation = NULL;

	allocate_rng_state(&(p->polymers[i].poly_state), p->args.pseudo_random_number_generator_arg);
	seed_rng_state(&(p->polymers[i].poly_state), p->args.rng_seed_arg,
		       i+n_polymer_offset, p->args.pseudo_random_number_generator_arg);
	}

    // Max safe move distance
    p->max_safe_jump = p->Lx/p->nx < p->Ly / p->ny ? p->Lx/p->nx : p->Ly / p->ny;
    p->max_safe_jump = p->max_safe_jump < p->Lz / p->nz ? p->max_safe_jump : p->Lz/p->nz;
    p->max_safe_jump *= 0.95;

    // Reference Harmonic Spring Cste
    const soma_scalar_t harmonic_spring_Cste =
	1.0 / sqrt(3.0 * (p->reference_Nbeads - 1.0));
    //Reference energy scale for harmonic springs.
    p->harmonic_normb =
	1.0 / (2.0 * harmonic_spring_Cste * harmonic_spring_Cste);


    p->n_cells = p->nx * p->ny * p->nz;
    uint64_t n_polymers_global_sum;
    MPI_Allreduce(&(p->n_polymers), &n_polymers_global_sum, 1,
		  MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_MPI_Comm);
    assert(p->n_polymers_global == n_polymers_global_sum);
    //Allocate Fields
    p->fields_unified =     (uint16_t *) malloc(p->n_cells*p->n_types*sizeof(uint16_t));
	if (p->fields_unified == NULL) {
	    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	    return -1;
	}

    p->old_fields_unified =     (uint16_t *) malloc(p->n_cells*p->n_types*sizeof(uint16_t));
	if (p->old_fields_unified == NULL) {
	    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	    return -1;
	}
    p->fields_32 = (uint32_t*)malloc(p->n_types * p->n_cells * sizeof(uint32_t));
    if (p->fields_32 == NULL) {
	    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	    return -1;
	}
    p->omega_field_unified = malloc(p->n_cells * p->n_types * sizeof(soma_scalar_t));
    if (p->omega_field_unified == NULL) {
	    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	    return -1;
	}

    p->tempfield =
	(soma_scalar_t *) malloc(p->nx * p->ny * p->nz * sizeof(soma_scalar_t));
    if (p->tempfield == NULL) {
	fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	return -1;
    }

    p->num_bead_type =
	(uint64_t *) malloc(p->n_types * sizeof(uint64_t));
    if (p->num_bead_type == NULL) {
	fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	return -1;
    }
    p->num_bead_type_local =
	(uint64_t *) malloc(p->n_types * sizeof(uint64_t));
    if (p->num_bead_type_local == NULL) {
	fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	return -1;
    }
    p->field_scaling_type = (soma_scalar_t *) malloc(p->n_types * sizeof(soma_scalar_t));
    if (p->field_scaling_type == NULL) {
	fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	return -1;
    }

    // Set all values to zero
    p->num_all_beads = 0;
    p->num_all_beads_local = 0;
    for (unsigned int i = 0; i < p->n_types; i++)
	p->num_bead_type_local[i] = 0;
    // Determine number of  different bead types
    for (uint64_t j = 0; j < p->n_polymers; j++) {	/*Loop over polymers */
      const unsigned int N = p->poly_arch[ p->poly_type_offset[p->polymers[j].type] ];
	for (unsigned int k = 0; k < N; k++) {	/*Loop over monomers */
	    const unsigned int type = get_particle_type(
		p->poly_arch[ p->poly_type_offset[p->polymers[j].type]+1+k]);
	    p->num_bead_type_local[type] += 1;
	    p->num_all_beads_local += 1;
	}
    }
    // Share p->num_all_beads
    MPI_Allreduce(&(p->num_all_beads_local), &(p->num_all_beads), 1,
		  MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_MPI_Comm);

    // Share p->num_bead_type
    for (unsigned int i = 0; i < p->n_types; i++) {
	MPI_Allreduce(&(p->num_bead_type_local[i]), &(p->num_bead_type[i]),
		      1, MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_MPI_Comm);
    }
    // Check if uint16_t density field is enough
    soma_scalar_t check_short = p->num_all_beads/p->n_cells;

    if (check_short > ((USHRT_MAX/100)*95)){
      fprintf(stderr, "ERROR: Cell-density above 95 Percent of USHRT_MAX\n");
	return -1;
    }
    // setting the geometry field, or if it's already initialized measure the free space
    uint64_t ncells = p->n_cells;
    if (p->area51 != NULL) {
        // substract the number of non free cells for the correct density scaling
        for (uint64_t i = 0; i < p->n_cells; i++)
           if ( p->area51[i] > 0 ) ncells--;
    }
    // Loop to calculate scaling parameter
    for (unsigned int i = 0; i < p->n_types; i++)
	p->field_scaling_type[i] =
	    (ncells / ((soma_scalar_t) p->num_all_beads));
    // Info for Ulrich: programm does take excluded volume into account now!
    p->n_accepts = 0;
    p->n_moves = 0;

    p->R = (soma_scalar_t*) malloc( p->n_types * sizeof(soma_scalar_t));
    for (unsigned int i = 0; i < p->n_types; i++){
      //! \todo kBT required.
      p->R[i] = sqrt( p->A[i] * 2);
    }

    // initialize inverse simulation cell parameters
    p->iLx = 1.0/p->Lx;
    p->iLy = 1.0/p->Ly;
    p->iLz = 1.0/p->Lz;

    p->sets = NULL; // Default init of the sets
    p->max_set_members = 0;
    if( p->args.iteration_alg_arg == iteration_alg_arg_SET)
      generate_independet_sets(p);


    init_autotuner(&(p->mc_autotuner));
    init_autotuner(&(p->cm_mc_autotuner));

    copyin_phase(p, false);
    // call update_fields routine
    if(p->bead_data_read)
	{
	update_density_fields(p);
	memcpy(p->old_fields_unified, p->fields_unified, p->n_cells*p->n_types*sizeof(uint16_t));
	}

    return 0;
}

int copyin_phase(struct Phase*const p, bool copyMSD)
    {
#ifdef SPEC_OPENACC
#ifdef SPEC_OPENACC 
#pragma acc enter data copyin(p[0:1])
#pragma acc enter data copyin(p->xn[0:p->n_types*p->n_types])
#pragma acc enter data copyin(p->polymers[0:p->n_polymers_storage])
#pragma acc enter data copyin(p->fields_unified[0:p->n_types*p->n_cells])
#pragma acc enter data copyin(p->old_fields_unified[0:p->n_types*p->n_cells])
#pragma acc enter data copyin(p->fields_32[0:p->n_types*p->n_cells])
#endif 
    if (p->area51 != NULL){
#ifdef SPEC_OPENACC 
#pragma acc enter data copyin(p->area51[0:p->n_cells])
#endif 
	}
#ifdef SPEC_OPENACC 
#pragma acc enter data copyin(p->omega_field_unified[0:p->n_cells*p->n_types])
#endif 
    if (p->external_field_unified != NULL){
#ifdef SPEC_OPENACC 
#pragma acc enter data copyin(p->external_field_unified[0:p->n_cells*p->n_types])
#endif 
	}
#ifdef SPEC_OPENACC 
#pragma acc enter data copyin(p->tempfield[0:p->n_cells])
#pragma acc enter data copyin(p->num_bead_type[0:p->n_types])
#pragma acc enter data copyin(p->num_bead_type_local[0:p->n_types])
#pragma acc enter data copyin(p->A[0:p->n_types])
#pragma acc enter data copyin(p->R[0:p->n_types])
#pragma acc enter data copyin(p->field_scaling_type[0:p->n_types])
#pragma acc enter data copyin(p->poly_type_offset[0:p->n_poly_type])
#pragma acc enter data copyin(p->poly_arch[0:p->poly_arch_length])
#endif 

    if(p->cm_a != NULL)
	{
#ifdef SPEC_OPENACC 
#pragma acc enter data copyin(p->cm_a[0:p->n_poly_type])
#endif 
	}
    if( p->sets != NULL)
	{
#ifdef SPEC_OPENACC 
#pragma acc enter data copyin(p->sets[0:p->n_poly_type])
#endif 
	for(unsigned int i=0; i < p->n_poly_type; i++)
	    {
#ifdef SPEC_OPENACC 
#pragma acc enter data copyin(p->sets[i].set_length[0:p->sets[i].n_sets])
#pragma acc enter data copyin(p->sets[i].sets[0:p->sets[i].n_sets*p->sets[i].max_member])
#endif 
	    }
	}

    // copy bulk polymer data
#pragma acc enter data copyin(p->allocator->all_Monomer.buf[0:p->allocator->all_Monomer.size])
    p->allocator->all_Monomer.device_buf = acc_deviceptr(p->allocator->all_Monomer.buf);
#if _OPENACC >= 201811
    if(copyMSD)
    {
#pragma acc enter data copyin(p->allocator->all_Monomer_msd.buf[0:p->allocator->all_Monomer_msd.size])
    }
    else
    {
#pragma acc enter data create(p->allocator->all_Monomer_msd.buf[0:p->allocator->all_Monomer_msd.size])
    }
    p->allocator->all_Monomer_msd.device_buf = acc_deviceptr(p->allocator->all_Monomer_msd.buf);
#endif
    if (p->allocator->all_MERSENNE_TWISTER_STATE.buf) {
#pragma acc enter data copyin(p->allocator->all_MERSENNE_TWISTER_STATE.buf[0:p->allocator->all_MERSENNE_TWISTER_STATE.size])
    p->allocator->all_MERSENNE_TWISTER_STATE.device_buf = acc_deviceptr(p->allocator->all_MERSENNE_TWISTER_STATE.buf);
    }
    if (p->allocator->all_MTTSTATE.buf) {
#pragma acc enter data copyin(p->allocator->all_MTTSTATE.buf[0:p->allocator->all_MTTSTATE.size])
    p->allocator->all_MTTSTATE.device_buf = acc_deviceptr(p->allocator->all_MTTSTATE.buf);
    }
    if (p->allocator->all_RNG_STATE.buf) {
#pragma acc enter data copyin(p->allocator->all_RNG_STATE.buf[0:p->allocator->all_RNG_STATE.size])
    p->allocator->all_RNG_STATE.device_buf = acc_deviceptr(p->allocator->all_RNG_STATE.buf);
    }
    if (p->allocator->all_uint_t.buf) {
#pragma acc enter data copyin(p->allocator->all_uint_t.buf[0:p->allocator->all_uint_t.size])
    p->allocator->all_uint_t.device_buf = acc_deviceptr(p->allocator->all_uint_t.buf);
    }
#pragma acc enter data copyin(p->allocator[0:1])
#if 0
    for(uint64_t i=0; i < p->n_polymers; i++)
	{
	Polymer* const poly = &(p->polymers[i]);
	copyin_polymer(p, poly);
	}
#endif

    return 0;
#else  // OPENMP TARGET
#ifdef SPEC_OPENMP_TARGET
#ifdef SPEC_OPENMP_TARGET 
#pragma omp target enter data map(to:p[0:1])
#pragma omp target enter data map(to:p->xn[0:p->n_types*p->n_types])
#pragma omp target enter data map(to:p->polymers[0:p->n_polymers_storage])
#pragma omp target enter data map(to:p->fields_unified[0:p->n_types*p->n_cells])
#pragma omp target enter data map(to:p->old_fields_unified[0:p->n_types*p->n_cells])
//FIXME Issue if enabled
#pragma omp target enter data map(to:p->fields_32[0:p->n_types*p->n_cells])
#endif 
    if (p->area51 != NULL){
#ifdef SPEC_OPENMP_TARGET 
#pragma omp target enter data map(to:p->area51[0:p->n_cells])
#endif 
	}
#ifdef SPEC_OPENMP_TARGET 
#pragma omp target enter data map(to:p->omega_field_unified[0:p->n_cells*p->n_types])
#endif 
    if (p->external_field_unified != NULL){
#ifdef SPEC_OPENMP_TARGET 
#pragma omp target enter data map(to:p->external_field_unified[0:p->n_cells*p->n_types])
#endif 
	}
#ifdef SPEC_OPENMP_TARGET 
#pragma omp target enter data map(to:p->tempfield[0:p->n_cells])
#pragma omp target enter data map(to:p->num_bead_type[0:p->n_types])
#pragma omp target enter data map(to:p->num_bead_type_local[0:p->n_types])
#pragma omp target enter data map(to:p->A[0:p->n_types])
#pragma omp target enter data map(to:p->R[0:p->n_types])
#pragma omp target enter data map(to:p->field_scaling_type[0:p->n_types])
#pragma omp target enter data map(to:p->poly_type_offset[0:p->n_poly_type])
#pragma omp target enter data map(to:p->poly_arch[0:p->poly_arch_length])
#endif 

    if(p->cm_a != NULL)
	{
#ifdef SPEC_OPENMP_TARGET 
#pragma omp target enter data map(to:p->cm_a[0:p->n_poly_type])
#endif 
	}
    if( p->sets != NULL)
	{
#ifdef SPEC_OPENMP_TARGET 
#pragma omp target enter data map(to:p->sets[0:p->n_poly_type])
#endif 
	for(unsigned int i=0; i < p->n_poly_type; i++)
	    {
#ifdef SPEC_OPENMP_TARGET 
#pragma omp target enter data map(to:p->sets[i].set_length[0:p->sets[i].n_sets])
#pragma omp target enter data map(to:p->sets[i].sets[0:p->sets[i].n_sets*p->sets[i].max_member])
#endif 
	    }
	}
    // copy bulk polymer data
#pragma omp target enter data map(to:p->allocator[0:1])
#if _OPENMP >= 201811
#pragma omp target update to(p->allocator->all_Monomer.buf[0:p->allocator->all_Monomer.size])
#else
// omp_target_associate_ptr is broken regarding initial_device before OpenMP 5.0, do memcpy
    omp_target_memcpy(p->allocator->all_Monomer.device_buf, p->allocator->all_Monomer.buf
	    , p->allocator->all_Monomer.size*sizeof(Monomer), 0, 0
	    , soma_get_device(), omp_get_initial_device()
	);
#endif
    if(copyMSD)
    {
#if _OPENMP >= 201811
#pragma omp target update to(p->allocator->all_Monomer_msd.buf[0:p->allocator->all_Monomer_msd.size])
#else
	omp_target_memcpy(p->allocator->all_Monomer_msd.device_buf, p->allocator->all_Monomer_msd.buf
		, p->allocator->all_Monomer_msd.size*sizeof(Monomer), 0, 0
		, soma_get_device(), omp_get_initial_device()
	    );
#endif
    }
#if 0 //! FIXME: would need omp_target_memcpy and device ptr translation/indices but is not used anyway
    if (p->allocator->all_MERSENNE_TWISTER_STATE.buf) {
#pragma omp target update to(p->allocator->all_MERSENNE_TWISTER_STATE.buf[0:p->allocator->all_MERSENNE_TWISTER_STATE.size])
    }
    if (p->allocator->all_MTTSTATE.buf) {
#pragma omp target update to(p->allocator->all_MTTSTATE.buf[0:p->allocator->all_MERSENNE_TWISTER_STATE.size])
    }
    if (p->allocator->all_RNG_STATE.buf) {
#pragma omp target update to(p->allocator->all_RNG_STATE.buf[0:p->allocator->all_MERSENNE_TWISTER_STATE.size])
    }
    if (p->allocator->all_uint_t.buf) {
#pragma omp target update to(p->allocator->all_uint_t.buf[0:p->allocator->all_MERSENNE_TWISTER_STATE.size])
    }
#endif
#if 0
    for(uint64_t i=0; i < p->n_polymers; i++)
	{
	Polymer*const poly = &(p->polymers[i]);
	copyin_polymer(p, poly);
	}
#endif
    return 0;
#else
    return p->n_polymers*0+1;
#endif // SPEC_OPENMP_TARGET
#endif//SPEC_OPENACC
    }

int copyout_phase(struct Phase*const p)
    {
#ifdef SPEC_OPENACC
#ifdef SPEC_OPENACC 
#pragma acc exit data delete(p[0:1])
#pragma acc exit data delete(p->xn[0:p->n_types*p->n_types])
#pragma acc exit data delete(p->polymers[0:p->n_polymers_storage])
#pragma acc exit data delete(p->fields_unified[0:p->n_types*p->n_cells])
#pragma acc exit data delete(p->old_fields_unified[0:p->n_types*p->n_cells])
#pragma acc exit data delete(p->fields_32[0:p->n_types*p->n_cells])
#endif 
    if (p->area51 != NULL){
#ifdef SPEC_OPENACC 
#pragma acc exit data delete(p->area51[0:p->n_cells])
#endif 
	}
#ifdef SPEC_OPENACC 
#pragma acc exit data delete(p->omega_field_unified[0:p->n_cells*p->n_types])
#endif 

    if (p->external_field_unified != NULL){
#ifdef SPEC_OPENACC 
#pragma acc exit data delete(p->external_field_unified[0:p->n_cells*p->n_types])
#endif 
	}
#ifdef SPEC_OPENACC 
#pragma acc exit data delete(p->tempfield[0:p->n_cells])
#pragma acc exit data delete(p->num_bead_type[0:p->n_types])
#pragma acc exit data delete(p->num_bead_type_local[0:p->n_types])
#pragma acc exit data delete(p->A[0:p->n_types])
#pragma acc exit data delete(p->R[0:p->n_types])
#pragma acc exit data delete(p->field_scaling_type[0:p->n_types])
#pragma acc exit data delete(p->poly_type_offset[0:p->n_poly_type])
#pragma acc exit data delete(p->poly_arch[0:p->poly_arch_length])
#endif 

    if(p->cm_a != NULL)
	{
#ifdef SPEC_OPENACC 
#pragma acc exit data delete(p->cm_a[0:p->n_poly_type])
#endif 
	}
    if( p->sets != NULL)
	{
#ifdef SPEC_OPENACC 
#pragma acc exit data delete(p->sets[0:p->n_poly_type])
#endif 
	for(unsigned int i=0; i < p->n_poly_type; i++)
	    {
#ifdef SPEC_OPENACC 
#pragma acc exit data delete(p->sets[i].set_length[0:p->sets[i].n_sets])
#pragma acc exit data delete(p->sets[i].sets[0:p->sets[i].n_sets*p->sets[i].max_member])
#endif 
	    }
	}

    // copy bulk polymer data
#pragma acc exit data copyout(p->allocator->all_Monomer.buf[0:p->allocator->all_Monomer.size])
    if (p->allocator->all_MERSENNE_TWISTER_STATE.buf) {
#pragma acc exit data copyout(p->allocator->all_MERSENNE_TWISTER_STATE.buf[0:p->allocator->all_MERSENNE_TWISTER_STATE.size])
    }
    if (p->allocator->all_MTTSTATE.buf) {
#pragma acc exit data copyout(p->allocator->all_MTTSTATE.buf[0:p->allocator->all_MERSENNE_TWISTER_STATE.size])
    }
    if (p->allocator->all_RNG_STATE.buf) {
#pragma acc exit data copyout(p->allocator->all_RNG_STATE.buf[0:p->allocator->all_MERSENNE_TWISTER_STATE.size])
    }
    if (p->allocator->all_uint_t.buf) {
#pragma acc exit data copyout(p->allocator->all_uint_t.buf[0:p->allocator->all_MERSENNE_TWISTER_STATE.size])
    }
#if 0
    for(uint64_t i=0; i < p->n_polymers; i++)
	{
	Polymer*const poly = &(p->polymers[i]);
	copyout_polymer(p, poly);
	}
#endif
    return 0;
#else
#ifdef SPEC_OPENMP_TARGET
#ifdef SPEC_OPENMP_TARGET 
#pragma omp target exit data map(delete:p[0:1])
#pragma omp target exit data map(delete:p->xn[0:p->n_types*p->n_types])
#pragma omp target exit data map(delete:p->polymers[0:p->n_polymers_storage])
#pragma omp target exit data map(delete:p->fields_unified[0:p->n_types*p->n_cells])
#pragma omp target exit data map(delete:p->old_fields_unified[0:p->n_types*p->n_cells])
#pragma omp target exit data map(delete:p->fields_32[0:p->n_types*p->n_cells])
#endif 
    if (p->area51 != NULL){
#ifdef SPEC_OPENMP_TARGET 
#pragma omp target exit data map(delete:p->area51[0:p->n_cells])
#endif 
	}
#ifdef SPEC_OPENMP_TARGET 
#pragma omp target exit data map(delete:p->omega_field_unified[0:p->n_cells*p->n_types])
#endif 

    if (p->external_field_unified != NULL){
#ifdef SPEC_OPENMP_TARGET 
#pragma omp target exit data map(delete:p->external_field_unified[0:p->n_cells*p->n_types])
#endif 
	}
#ifdef SPEC_OPENMP_TARGET 
#pragma omp target exit data map(delete:p->tempfield[0:p->n_cells])
#pragma omp target exit data map(delete:p->num_bead_type[0:p->n_types])
#pragma omp target exit data map(delete:p->num_bead_type_local[0:p->n_types])
#pragma omp target exit data map(delete:p->A[0:p->n_types])
#pragma omp target exit data map(delete:p->R[0:p->n_types])
#pragma omp target exit data map(delete:p->field_scaling_type[0:p->n_types])
#pragma omp target exit data map(delete:p->poly_type_offset[0:p->n_poly_type])
#pragma omp target exit data map(delete:p->poly_arch[0:p->poly_arch_length])
#endif 

    if(p->cm_a != NULL)
	{
#ifdef SPEC_OPENMP_TARGET 
#pragma omp target exit data map(delete:p->cm_a[0:p->n_poly_type])
#endif 
	}
    if( p->sets != NULL)
	{
#ifdef SPEC_OPENMP_TARGET 
#pragma omp target exit data map(delete:p->sets[0:p->n_poly_type])
#endif 
	for(unsigned int i=0; i < p->n_poly_type; i++)
	    {
#ifdef SPEC_OPENMP_TARGET 
#pragma omp target exit data map(delete:p->sets[i].set_length[0:p->sets[i].n_sets])
#pragma omp target exit data map(delete:p->sets[i].sets[0:p->sets[i].n_sets*p->sets[i].max_member])
#endif 
	    }
	}

    // copy bulk polymer data
#if _OPENMP >= 201811
#pragma omp target update from(p->allocator->all_Monomer.buf[0:p->allocator->all_Monomer.size])
#else
// omp_target_associate_ptr is broken regarding initial_device before OpenMP 5.0, do memcpy
    omp_target_memcpy(p->allocator->all_Monomer.buf, p->allocator->all_Monomer.device_buf
	    , p->allocator->all_Monomer.size*sizeof(Monomer), 0, 0
	    , omp_get_initial_device(), soma_get_device()
	);
#endif
#if 0 //! FIXME: would need omp_target_memcpy and device ptr translation/indices but is not used anyway
    if (p->allocator->all_MERSENNE_TWISTER_STATE.buf) {
#pragma omp target exit data map(from:p->allocator->all_MERSENNE_TWISTER_STATE.buf[0:p->allocator->all_MERSENNE_TWISTER_STATE.size])
    }
    if (p->allocator->all_MTTSTATE.buf) {
#pragma omp target exit data map(from:p->allocator->all_MTTSTATE.buf[0:p->allocator->all_MERSENNE_TWISTER_STATE.size])
    }
    if (p->allocator->all_RNG_STATE.buf) {
#pragma omp target exit data map(from:p->allocator->all_RNG_STATE.buf[0:p->allocator->all_MERSENNE_TWISTER_STATE.size])
    }
    if (p->allocator->all_uint_t.buf) {
#pragma omp target exit data map(from:p->allocator->all_uint_t.buf[0:p->allocator->all_MERSENNE_TWISTER_STATE.size])
    }
#endif
#if 0
    for(uint64_t i=0; i < p->n_polymers; i++)
	{
	Polymer*const poly = &(p->polymers[i]);
	copyout_polymer(p, poly);
	}
#endif
    return 0;
#else
    return p->n_polymers*0 +1;
#endif//SPEC_OPENMP_TARGET
#endif//SPEC_OPENACC
    }

int free_phase(struct Phase * const p)
    {
    copyout_phase(p);

    /* de-allocate fields */
    free(p->omega_field_unified);
    free(p->tempfield);
    free(p->fields_unified);
    free(p->old_fields_unified);
    free(p->fields_32);
    free(p->num_bead_type);
    free(p->num_bead_type_local);
    free(p->field_scaling_type);
    free(p->A);
    free(p->R);
    free(p->end_mono);
    free(p->cm_a);

    cmdline_parser_free(&(p->args));

    if(p->sets != NULL)
      {
	for(unsigned int i=0; i < p->n_poly_type; i++)
	  {
	    free(p->sets[i].set_length);
	    free(p->sets[i].sets);
	  }
	free(p->sets);
      }

    /* free polymers */
    for (uint64_t i = 0; i < p->n_polymers; i++)
	{
	Polymer*const poly = &(p->polymers[i]);
	free_polymer(p, poly);
	}

    free(p->polymers);

    free(p->poly_type_offset);
    free(p->poly_arch);

    /* de allocate XN interaction matrix */
    free(p->xn);

    if (p->area51 != NULL) {
	free(p->area51);
    }

    if (p->external_field_unified != NULL){
	free(p->external_field_unified);
    }

    return 0;
}

int update_self_phase(const Phase * const p)
    {
    static unsigned int last_time_call = 0;
    if (last_time_call == 0 || p->time > last_time_call)
	last_time_call = p->time;
    else			//Quick exit, because the property has already been calculated for the time step.
	return 1;

    // Not pointer members are expected to not change on device

#ifdef SPEC_OPENACC 
#pragma acc update self(p->xn[0:p->n_types*p->n_types])
#endif 
#ifdef SPEC_OPENMP_TARGET
#pragma omp target update from(p->xn[0:p->n_types*p->n_types])
#endif 

#ifdef SPEC_OPENACC
#pragma acc update self(p->allocator->all_Monomer.buf[0:p->allocator->all_Monomer.size])
#elif defined SPEC_OPENMP_TARGET
#if _OPENMP >= 201811
#pragma omp target update from(p->allocator->all_Monomer.buf[0:p->allocator->all_Monomer.size])
#else
// omp_target_associate_ptr is broken regarding initial_device before OpenMP 5.0, do memcpy
    omp_target_memcpy(p->allocator->all_Monomer.buf, p->allocator->all_Monomer.device_buf
	    , p->allocator->all_Monomer.size*sizeof(Monomer), 0, 0
	    , omp_get_initial_device(), soma_get_device()
	);
#endif
#endif
#if 0 //! FIXME: would need omp_target_memcpy but is not used anyway
    if (p->allocator->all_MERSENNE_TWISTER_STATE.buf) {
#ifdef SPEC_OPENACC
#pragma acc update self(p->allocator->all_MERSENNE_TWISTER_STATE.buf[0:p->allocator->all_MERSENNE_TWISTER_STATE.size])
#elif defined SPEC_OPENMP_TARGET
#pragma omp target update from(p->allocator->all_MERSENNE_TWISTER_STATE.buf[0:p->allocator->all_MERSENNE_TWISTER_STATE.size])
#endif
    }
    if (p->allocator->all_MTTSTATE.buf) {
#ifdef SPEC_OPENACC
#pragma acc update self(p->allocator->all_MTTSTATE.buf[0:p->allocator->all_MERSENNE_TWISTER_STATE.size])
#elif defined SPEC_OPENMP_TARGET
#pragma omp target update from(p->allocator->all_MTTSTATE.buf[0:p->allocator->all_MERSENNE_TWISTER_STATE.size])
#endif
    }
    if (p->allocator->all_RNG_STATE.buf) {
#ifdef SPEC_OPENACC
#pragma acc update self(p->allocator->all_RNG_STATE.buf[0:p->allocator->all_MERSENNE_TWISTER_STATE.size])
#elif defined SPEC_OPENMP_TARGET
#pragma omp target update from(p->allocator->all_RNG_STATE.buf[0:p->allocator->all_MERSENNE_TWISTER_STATE.size])
#endif
    }
    if (p->allocator->all_uint_t.buf) {
#ifdef SPEC_OPENACC
#pragma acc update self(p->allocator->all_uint_t.buf[0:p->allocator->all_MERSENNE_TWISTER_STATE.size])
#elif defined SPEC_OPENMP_TARGET
#pragma omp target update from(p->allocator->all_uint_t.buf[0:p->allocator->all_MERSENNE_TWISTER_STATE.size])
#endif
    }
#endif
#if 0
    for(uint64_t i=0; i< p->n_polymers; i++)
	update_self_polymer(p, p->polymers+i);
#endif

#ifdef SPEC_OPENACC 
#pragma acc update self(p->fields_unified[0:p->n_cells*p->n_types])
#pragma acc update self(p->old_fields_unified[0:p->n_types*p->n_cells])
#pragma acc update self(p->fields_32[0:p->n_types*p->n_cells])
#endif 
#ifdef SPEC_OPENMP_TARGET
#pragma omp target update from(p->fields_unified[0:p->n_cells*p->n_types])
#pragma omp target update from(p->old_fields_unified[0:p->n_types*p->n_cells])
#pragma omp target update from(p->fields_32[0:p->n_types*p->n_cells])
#endif 


    if (p->area51 != NULL)
	{
#ifdef SPEC_OPENACC 
#pragma acc update self(p->area51[0:p->n_cells])
#endif 
#ifdef SPEC_OPENMP_TARGET
#pragma omp target update from(p->area51[0:p->n_cells])
#endif 
	}
#ifdef SPEC_OPENACC 
#pragma acc update self(p->omega_field_unified[0:p->n_cells*p->n_types])
#endif 
#ifdef SPEC_OPENMP_TARGET
#pragma omp target update from(p->omega_field_unified[0:p->n_cells*p->n_types])
#endif 
    if (p->external_field_unified != NULL)
	{
#ifdef SPEC_OPENACC 
#pragma acc update self(p->external_field_unified[0:p->n_cells*p->n_types])
#endif 
#ifdef SPEC_OPENMP_TARGET
#pragma omp target update from(p->external_field_unified[0:p->n_cells*p->n_types])
#endif 
	}
#ifdef SPEC_OPENACC 
#pragma acc update self(p->tempfield[0:p->n_cells])
#pragma acc update self(p->num_bead_type[0:p->n_types])
#pragma acc update self(p->num_bead_type_local[0:p->n_types])
#pragma acc update self(p->A[0:p->n_types])
#pragma acc update self(p->R[0:p->n_types])
#pragma acc update self(p->field_scaling_type[0:p->n_types])
#pragma acc update self(p->poly_type_offset[0:p->n_poly_type])
#pragma acc update self(p->poly_arch[0:p->poly_arch_length])
#endif 
#ifdef SPEC_OPENMP_TARGET
#pragma omp target update from(p->tempfield[0:p->n_cells])
#pragma omp target update from(p->num_bead_type[0:p->n_types])
#pragma omp target update from(p->num_bead_type_local[0:p->n_types])
#pragma omp target update from(p->A[0:p->n_types])
#pragma omp target update from(p->R[0:p->n_types])
#pragma omp target update from(p->field_scaling_type[0:p->n_types])
#pragma omp target update from(p->poly_type_offset[0:p->n_poly_type])
#pragma omp target update from(p->poly_arch[0:p->poly_arch_length])
#endif 
    if(p->cm_a != NULL)
	{
#ifdef SPEC_OPENACC 
#pragma acc update self(p->cm_a[0:p->n_poly_type])
#endif 
#ifdef SPEC_OPENMP_TARGET
#pragma omp target update from(p->cm_a[0:p->n_poly_type])
#endif 
	}

    //SETS are not updated to host

    return p->n_polymers*0+1;
    }

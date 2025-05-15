/* Copyright (C) 2016-2017 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016-2017 Marcel Langenberg
   Copyright (C) 2016 Fabien Leonforte
   Copyright (C) 2016 Juan Orozco
   Copyright (C) 2016 Yongzhi Ren
   Copyright (C) 2016 N. Harshavardhan Reddy

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

//! \file ana.c
//! \brief Implementation of ana.h

#include "ana.h"
#include <time.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include "phase.h"
#include "polymer.h"
#include "mesh.h"
#include "io.h"
#include "allocator.h"

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

#if ( defined SPEC_OPENACC && ( ! defined SPEC_NO_VAR_ARRAY_REDUCE ) ) \
	|| ( defined SPEC_OPENMP_TARGET && ( ! defined SPEC_NO_VAR_ARRAY_REDUCE ) )
#define MONOMER_BUF device_buf
#else
#define MONOMER_BUF buf
#endif

void calc_Re(const struct Phase * p, soma_scalar_t *const result)
    {
    uint64_t *const counter = (uint64_t*)malloc( p->n_poly_type*sizeof(uint64_t));
    if(counter == NULL)
	{
	fprintf(stderr,"MALLOC ERROR: %s:%d\n",__FILE__,__LINE__);
	return;
	}
    memset(counter,0, p->n_poly_type*sizeof(uint64_t));
    const int n_result = 4*p->n_poly_type;
    memset(result, 0 , n_result * sizeof(soma_scalar_t));

#if defined SPEC_OPENACC && ( ! defined SPEC_NO_VAR_ARRAY_REDUCE )
#pragma acc parallel loop independent present(p) \
    reduction(+:counter[0:p->n_poly_type],result[0:n_result])
#endif
#if defined SPEC_OPENMP_TARGET && ( ! defined SPEC_NO_VAR_ARRAY_REDUCE )
#pragma omp target teams distribute parallel for \
    map(tofrom:counter[0:p->n_poly_type],result[0:n_result]) \
    reduction(+:counter[0:p->n_poly_type],result[0:n_result])
#endif
#if defined SPEC_OPENMP && ( ! defined SPEC_NO_VAR_ARRAY_REDUCE )
#pragma omp parallel for \
    reduction(+:counter[0:p->n_poly_type],result[0:n_result])
#endif
    for (uint64_t npoly = 0; npoly < p->n_polymers; npoly++) {
	const unsigned int type = p->polymers[npoly].type;

	const unsigned int start = p->end_mono[ type*2 + 0];
	const unsigned int end = p->end_mono[ type*2 + 1];
	const soma_scalar_t dx = p->allocator->all_Monomer.MONOMER_BUF[p->polymers[npoly].beads+start].x - p->allocator->all_Monomer.MONOMER_BUF[p->polymers[npoly].beads+end].x;
	const soma_scalar_t dy = p->allocator->all_Monomer.MONOMER_BUF[p->polymers[npoly].beads+start].y - p->allocator->all_Monomer.MONOMER_BUF[p->polymers[npoly].beads+end].y;
	const soma_scalar_t dz = p->allocator->all_Monomer.MONOMER_BUF[p->polymers[npoly].beads+start].z - p->allocator->all_Monomer.MONOMER_BUF[p->polymers[npoly].beads+end].z;

	result[type*4 + 0 ] += dx*dx + dy*dy + dz*dz;
	result[type*4 + 1 ] += dx*dx ;
	result[type*4 + 2 ] += dy*dy ;
	result[type*4 + 3 ] += dz*dz ;
	counter[type] += 1;
	}

    MPI_Allreduce(MPI_IN_PLACE, result, 4*p->n_poly_type, MPI_SOMA_SCALAR, MPI_SUM,
		  p->info_MPI.SOMA_MPI_Comm);
    MPI_Allreduce(MPI_IN_PLACE, counter, p->n_poly_type, MPI_UINT64_T, MPI_SUM,
		  p->info_MPI.SOMA_MPI_Comm);

    for(unsigned int type=0 ; type < p->n_poly_type; type++)
	for(unsigned int i=0; i < 4; i++)
	    if( counter[type] > 0)
		result[type*4+i] /= (soma_scalar_t) counter[type];

    free(counter);
    }

void calc_dvar(const struct Phase * p, soma_scalar_t *dvar)
    {
    soma_scalar_t var = 0.0;
    for (uint64_t index=0; index < p->n_cells*p->n_types; index++) {
    var += (p->fields_unified[index] - p->old_fields_unified[index])*(p->fields_unified[index] - p->old_fields_unified[index]);
    }
    memcpy(p->old_fields_unified, p->fields_unified, p->n_cells*p->n_types*sizeof(uint16_t));
    var = var/(p->n_cells*p->n_types);
    *dvar = var;
    }

void calc_Rg(const struct Phase *p, soma_scalar_t *const result)
    {
    uint64_t *const counter = (uint64_t*)malloc( p->n_poly_type*sizeof(uint64_t));
    if(counter == NULL)
	{
	fprintf(stderr,"MALLOC ERROR: %s:%d\n",__FILE__,__LINE__);
	return;
	}
    memset(counter,0, p->n_poly_type*sizeof(uint64_t));
    const int n_result = 4*p->n_poly_type;
    memset(result, 0 , n_result * sizeof(soma_scalar_t));

#if defined SPEC_OPENACC && ( ! defined SPEC_NO_VAR_ARRAY_REDUCE )
#pragma acc parallel loop independent present(p) \
    reduction(+:counter[0:p->n_poly_type],result[0:n_result])
#endif
#if defined SPEC_OPENMP_TARGET && ( ! defined SPEC_NO_VAR_ARRAY_REDUCE )
#pragma omp target teams distribute parallel for \
    map(tofrom:counter[0:p->n_poly_type],result[0:n_result]) \
    reduction(+:counter[0:p->n_poly_type],result[0:n_result])
#endif
#if defined SPEC_OPENMP && ( ! defined SPEC_NO_VAR_ARRAY_REDUCE )
#pragma omp parallel for \
    reduction(+:counter[0:p->n_poly_type],result[0:n_result])
#endif
    for (uint64_t ipoly = 0; ipoly < p->n_polymers; ipoly++) {
	const unsigned int type =p->polymers[ipoly].type;
	const unsigned int N = p->poly_arch[ p->poly_type_offset[ type ] ];
	soma_scalar_t xcm = 0.;
	soma_scalar_t ycm = 0.;
	soma_scalar_t zcm = 0.;
	soma_scalar_t x2  = 0.;
	soma_scalar_t y2  = 0.;
	soma_scalar_t z2  = 0.;
	for (unsigned int ibead = 0; ibead < N; ibead++) {
	    const soma_scalar_t x1   = p->allocator->all_Monomer.MONOMER_BUF[p->polymers[ipoly].beads+ibead].x;
	    const soma_scalar_t y1   = p->allocator->all_Monomer.MONOMER_BUF[p->polymers[ipoly].beads+ibead].y;
	    const soma_scalar_t z1   = p->allocator->all_Monomer.MONOMER_BUF[p->polymers[ipoly].beads+ibead].z;
	    xcm += x1;
	    ycm += y1;
	    zcm += z1;
	    x2  += x1*x1;
	    y2  += y1*y1;
	    z2  += z1*z1;
	    }
	xcm /= (soma_scalar_t)(N);
	ycm /= (soma_scalar_t)(N);
	zcm /= (soma_scalar_t)(N);

	result[type*4 + 0] += (x2/(soma_scalar_t)(N) - xcm*xcm) +
	    (y2/(soma_scalar_t)(N) - ycm*ycm) +
	    (z2/(soma_scalar_t)(N) - zcm*zcm);
	result[type*4 + 1] += (x2/(soma_scalar_t)(N) - xcm*xcm);
	result[type*4 + 2] += (y2/(soma_scalar_t)(N) - ycm*ycm);
	result[type*4 + 3] += (z2/(soma_scalar_t)(N) - zcm*zcm);
	counter[type] += 1;
	}

    MPI_Allreduce(MPI_IN_PLACE, result, 4*p->n_poly_type, MPI_SOMA_SCALAR, MPI_SUM,
		  p->info_MPI.SOMA_MPI_Comm);
    MPI_Allreduce(MPI_IN_PLACE, counter, p->n_poly_type, MPI_UINT64_T, MPI_SUM,
		  p->info_MPI.SOMA_MPI_Comm);

    for(unsigned int type=0 ; type < p->n_poly_type; type++)
	for(unsigned int i=0; i < 4; i++)
	    if( counter[type] > 0)
		result[type*4 +i] /= (soma_scalar_t) counter[type];
    free(counter);
    }

void calc_anisotropy(const struct Phase * p, soma_scalar_t *const result)
    {
    uint64_t *const counter = (uint64_t*)malloc( p->n_poly_type*sizeof(uint64_t));
    if(counter == NULL)
	{
	fprintf(stderr,"MALLOC ERROR: %s:%d\n",__FILE__,__LINE__);
	return;
	}
    memset(counter,0, p->n_poly_type*sizeof(uint64_t));
    const int n_result = 6*p->n_poly_type;
    memset(result, 0 , n_result * sizeof(soma_scalar_t));

    // loop over local chains
#if defined SPEC_OPENACC && ( ! defined SPEC_NO_VAR_ARRAY_REDUCE )
#pragma acc parallel loop independent present(p) \
    reduction(+:counter[0:p->n_poly_type],result[0:n_result])
#endif
#if defined SPEC_OPENMP_TARGET && ( ! defined SPEC_NO_VAR_ARRAY_REDUCE )
#pragma omp target teams distribute parallel for \
    map(tofrom:counter[0:p->n_poly_type],result[0:n_result]) \
    reduction(+:counter[0:p->n_poly_type],result[0:n_result])
#endif
#if defined SPEC_OPENMP && ( ! defined SPEC_NO_VAR_ARRAY_REDUCE )
#pragma omp parallel for \
    reduction(+:counter[0:p->n_poly_type],result[0:n_result])
#endif
    for (uint64_t ipoly = 0; ipoly < p->n_polymers; ipoly++){
	const Polymer*const poly = &(p->polymers[ipoly]);
	const unsigned int type = poly->type;
	const unsigned int N = p->poly_arch[ p->poly_type_offset[ type ] ];

	// loop over beads in this chain
#ifdef SPEC_OPENACC
#pragma acc loop seq
#endif
	for (unsigned int ibead = 0; ibead < N; ibead++) {

	    const soma_scalar_t x1 = p->allocator->all_Monomer.MONOMER_BUF[poly->beads+ibead].x;
	    const soma_scalar_t y1 = p->allocator->all_Monomer.MONOMER_BUF[poly->beads+ibead].y;
	    const soma_scalar_t z1 = p->allocator->all_Monomer.MONOMER_BUF[poly->beads+ibead].z;

	    // loop over bonds of this bead
	    const int start = get_bondlist_offset(p->poly_arch[p->poly_type_offset[type] + ibead + 1]);
	    if(start > 0){
		int i = start;
		int end;

		do{
		    const uint32_t bn = p->poly_arch[i++];
		    const int info = bn;
		    end = get_end(info);
		    const int offset = get_offset(info);
		    const int neighbour_id = ibead + offset;
		    const unsigned int jbead = neighbour_id;

		    const soma_scalar_t x2 = p->allocator->all_Monomer.MONOMER_BUF[p->polymers[ipoly].beads+jbead].x;
		    const soma_scalar_t y2 = p->allocator->all_Monomer.MONOMER_BUF[p->polymers[ipoly].beads+jbead].y;
		    const soma_scalar_t z2 = p->allocator->all_Monomer.MONOMER_BUF[p->polymers[ipoly].beads+jbead].z;

		    const soma_scalar_t bx = x2 - x1;
		    const soma_scalar_t by = y2 - y1;
		    const soma_scalar_t bz = z2 - z1;

		    result[type*6 + 0] += bx * bx;
		    result[type*6 + 1] += by * by;
		    result[type*6 + 2] += bz * bz;
		    result[type*6 + 3] += bx * by;
		    result[type*6 + 4] += by * bz;
		    result[type*6 + 5] += bz * bx;
		    counter[type] += 1;
		    }while( end == 0);
		}
	    }
	}

    MPI_Allreduce(MPI_IN_PLACE, result, 6*p->n_poly_type, MPI_SOMA_SCALAR, MPI_SUM,
		  p->info_MPI.SOMA_MPI_Comm);
    MPI_Allreduce(MPI_IN_PLACE, counter, p->n_poly_type, MPI_UINT64_T, MPI_SUM,
		  p->info_MPI.SOMA_MPI_Comm);

    for(unsigned int type=0 ; type < p->n_poly_type; type++)
	for(unsigned int i=0; i < 6; i++)
	    if( counter[type] > 0)
		result[type*6 +i] /= (soma_scalar_t) counter[type];
    free(counter);
    }

void calc_MSD(const struct Phase * p, soma_scalar_t *const result)
    {
    uint64_t *const counter = (uint64_t*)malloc( 2*p->n_poly_type*sizeof(uint64_t));
    if(counter == NULL)
	{
	fprintf(stderr,"MALLOC ERROR: %s:%d\n",__FILE__,__LINE__);
	return;
	}
    const int n_counter = 2*p->n_poly_type;
    memset(counter,0, n_counter*sizeof(uint64_t));
    const int n_result = 8*p->n_poly_type;
    memset(result, 0 , n_result * sizeof(soma_scalar_t));

    // Add up local displacement
#if defined SPEC_OPENACC && ( ! defined SPEC_NO_VAR_ARRAY_REDUCE )
#pragma acc parallel loop independent present(p) \
    reduction(+:counter[0:n_counter],result[0:n_result])
#endif
#if defined SPEC_OPENMP_TARGET && ( ! defined SPEC_NO_VAR_ARRAY_REDUCE )
#pragma omp target teams distribute parallel for \
    map(tofrom:counter[0:n_counter],result[0:n_result]) \
    reduction(+:counter[0:n_counter],result[0:n_result])
#endif
#if defined SPEC_OPENMP && ( ! defined SPEC_NO_VAR_ARRAY_REDUCE )
#pragma omp parallel for \
    reduction(+:counter[0:n_counter],result[0:n_result])
#endif
    for (uint64_t j = 0; j < p->n_polymers; j++) {	/*Loop over polymers */

        soma_scalar_t tx_c = 0.;
        soma_scalar_t ty_c = 0.;
        soma_scalar_t tz_c = 0.;
        soma_scalar_t mx_c = 0.;
        soma_scalar_t my_c = 0.;
        soma_scalar_t mz_c = 0.;
	const unsigned int type = p->polymers[j].type;
	const unsigned int N = p->poly_arch[ p->poly_type_offset[ type ] ];
	for (unsigned int k = 0; k < N; k++) {	/*Loop over monomers */
	    const soma_scalar_t tx = p->allocator->all_Monomer_msd.MONOMER_BUF[p->polymers[j].msd_beads+k].x - p->allocator->all_Monomer.MONOMER_BUF[p->polymers[j].beads+k].x;
	    const soma_scalar_t ty = p->allocator->all_Monomer_msd.MONOMER_BUF[p->polymers[j].msd_beads+k].y - p->allocator->all_Monomer.MONOMER_BUF[p->polymers[j].beads+k].y;
	    const soma_scalar_t tz = p->allocator->all_Monomer_msd.MONOMER_BUF[p->polymers[j].msd_beads+k].z - p->allocator->all_Monomer.MONOMER_BUF[p->polymers[j].beads+k].z;

            // Add up values for chain diffusion
            tx_c += p->allocator->all_Monomer.MONOMER_BUF[p->polymers[j].beads+k].x;
            ty_c += p->allocator->all_Monomer.MONOMER_BUF[p->polymers[j].beads+k].y;
            tz_c += p->allocator->all_Monomer.MONOMER_BUF[p->polymers[j].beads+k].z;
            mx_c += p->allocator->all_Monomer_msd.MONOMER_BUF[p->polymers[j].msd_beads+k].x;
            my_c += p->allocator->all_Monomer_msd.MONOMER_BUF[p->polymers[j].msd_beads+k].y;
            mz_c += p->allocator->all_Monomer_msd.MONOMER_BUF[p->polymers[j].msd_beads+k].z;

	    counter[type*2 +0] += 1;
	    result[type*8 + 0] += tx * tx;
	    result[type*8 + 1] += ty * ty;
	    result[type*8 + 2] += tz * tz;
	    result[type*8 + 3] += (tx * tx + ty * ty + tz * tz);
	    }
	counter[type*2 + 1] += 1;
        result[type*8 + 4] += (tx_c-mx_c)*(tx_c-mx_c)/(N*N);
        result[type*8 + 5] += (ty_c-my_c)*(ty_c-my_c)/(N*N);
        result[type*8 + 6] += (tz_c-mz_c)*(tz_c-mz_c)/(N*N);
        result[type*8 + 7] += ((tx_c-mx_c)*(tx_c-mx_c) +(ty_c-my_c)*(ty_c-my_c) + (tz_c-mz_c)*(tz_c-mz_c) )/(N*N);
	}

    MPI_Allreduce(MPI_IN_PLACE, result, 8*p->n_poly_type, MPI_SOMA_SCALAR, MPI_SUM,
		  p->info_MPI.SOMA_MPI_Comm);
    MPI_Allreduce(MPI_IN_PLACE, counter, 2*p->n_poly_type, MPI_UINT64_T, MPI_SUM,
		  p->info_MPI.SOMA_MPI_Comm);

    //Looping over twice the number of poly types. But loop over half the elements
    // 8/2 = 4, because the norm for first and second half differ.
    for(unsigned int type=0 ; type < 2*p->n_poly_type; type++)
	for(unsigned int i=0; i < 4; i++)
	    if( counter[type] > 0)
		result[type*4 +i] /= (soma_scalar_t) counter[type];
    free(counter);
    }

void calc_non_bonded_energy(const struct Phase*const p, soma_scalar_t*const non_bonded_energy)
    {
    memset(non_bonded_energy,0,p->n_types*sizeof(soma_scalar_t));

#if defined SPEC_OPENACC && ( ! defined SPEC_NO_VAR_ARRAY_REDUCE )
#pragma acc parallel loop independent present(p) \
    reduction(+:non_bonded_energy[0:p->n_types])
#endif
#if defined SPEC_OPENMP_TARGET && ( ! defined SPEC_NO_VAR_ARRAY_REDUCE )
#pragma omp target teams distribute parallel for collapse(2) \
    map(tofrom:non_bonded_energy[0:p->n_types]) \
    reduction(+:non_bonded_energy[0:p->n_types])
#endif
#if defined SPEC_OPENMP && ( ! defined SPEC_NO_VAR_ARRAY_REDUCE )
#pragma omp parallel for collapse(2) \
    reduction(+:non_bonded_energy[0:p->n_types])
#endif
    for(unsigned int type=0; type < p->n_types; type++)
	{
	for(uint64_t cell=0; cell < p->n_cells; cell++)
	    {
	    non_bonded_energy[type] +=
		p->omega_field_unified[cell + type*p->n_cells]
		* p->fields_unified[cell + type*p->n_cells ];
	    }
	}
    }

void calc_bonded_energy(const struct Phase*const p, soma_scalar_t*const bonded_energy)
    {
    const int number_soma_bond_types = NUMBER_SOMA_BOND_TYPES;
    memset(bonded_energy,0,number_soma_bond_types*sizeof(soma_scalar_t));

#if defined SPEC_OPENACC && ( ! defined SPEC_NO_VAR_ARRAY_REDUCE )
#pragma acc parallel loop independent present(p) \
    reduction(+:bonded_energy[0:number_soma_bond_types])
#endif
#if defined SPEC_OPENMP_TARGET && ( ! defined SPEC_NO_VAR_ARRAY_REDUCE )
#pragma omp target teams distribute parallel for \
    map(tofrom:bonded_energy[0:number_soma_bond_types]) \
    reduction(+:bonded_energy[0:number_soma_bond_types])
#endif
#if defined SPEC_OPENMP && ( ! defined SPEC_NO_VAR_ARRAY_REDUCE )
#pragma omp parallel for \
    reduction(+:bonded_energy[0:number_soma_bond_types])
#endif
    for(unsigned int poly=0; poly < p->n_polymers; poly++)
	{
	const unsigned int type =p->polymers[poly].type;
	const unsigned int N = p->poly_arch[ p->poly_type_offset[ type ] ];

	for(unsigned int mono=0; mono < N ; mono++)
	    {

	    const int start = get_bondlist_offset(p->poly_arch[p->poly_type_offset[type] + mono + 1]);

	    if(start > 0)
		{
		int i = start;
		unsigned int end;
		do{
		    const uint32_t info = p->poly_arch[i++];
		    end = get_end(info);
		    const unsigned int bond_type = get_bond_type(info);
		    const int offset = get_offset(info);
		    if( offset > 0) //Select each bond only once, i<j
			{
			const int mono_j = mono + offset;

			const soma_scalar_t dx = p->allocator->all_Monomer.MONOMER_BUF[p->polymers[poly].beads+mono].x - p->allocator->all_Monomer.MONOMER_BUF[p->polymers[poly].beads+mono_j].x;
			const soma_scalar_t dy = p->allocator->all_Monomer.MONOMER_BUF[p->polymers[poly].beads+mono].y - p->allocator->all_Monomer.MONOMER_BUF[p->polymers[poly].beads+mono_j].y;
			const soma_scalar_t dz = p->allocator->all_Monomer.MONOMER_BUF[p->polymers[poly].beads+mono].z - p->allocator->all_Monomer.MONOMER_BUF[p->polymers[poly].beads+mono_j].z;
			const soma_scalar_t r2 = dx*dx + dy*dy + dz*dz;

			soma_scalar_t energy = 0;
			soma_scalar_t scale = 1.;
			switch (bond_type) {
			case HARMONICVARIABLESCALE:
			    scale = p->harmonic_normb_variable_scale;
			    /* intentionally falls through */
			case HARMONIC:
			    energy = p->harmonic_normb * r2 *scale;
			    break;

#if ! ( defined SPEC_OPENACC || defined SPEC_OPENMP_TARGET )
			case STIFF:
			    fprintf(stderr,
				    "ERROR: %s:%d stiff bond not yet implemented.\n",
				    __FILE__, __LINE__);
			    break;
			default:
			    fprintf(stderr, "ERROR: %s:%d unknow bond type appeared %d\n",
				    __FILE__, __LINE__,bond_type);
			    break;
#endif
			    }
#ifndef SPEC_OPENACC
			assert( bond_type < NUMBER_SOMA_BOND_TYPES);
#endif
			bonded_energy[ bond_type ] += energy;
			}
		    }while( end == 0);
		}
	    }
	}
    if( p->info_MPI.current_core == 0)
	MPI_Reduce(MPI_IN_PLACE,bonded_energy,NUMBER_SOMA_BOND_TYPES,MPI_SOMA_SCALAR,MPI_SUM,0,p->info_MPI.SOMA_MPI_Comm);
    else
	MPI_Reduce(bonded_energy, NULL ,NUMBER_SOMA_BOND_TYPES,MPI_SOMA_SCALAR,MPI_SUM,0,p->info_MPI.SOMA_MPI_Comm);
    }

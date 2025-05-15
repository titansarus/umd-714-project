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

//! \file generate_positions.c
//! \brief Implementation of generate_positions.h

#include "generate_positions.h"
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "polymer.h"
#include "monomer.h"
#include "bond.h"
#include "phase.h"
#include "mc.h"
#include "mesh.h"
#include "allocator.h"
#include "device.h"

#ifdef SPEC_OPENACC
#include <openacc.h>
#endif

#define SOMA_POLY_CHAIN

static inline bool is_already_set(const uint64_t*const already_set, unsigned int i)
{
    bool ret = (already_set[i/64]>>(i%64))&1;
    return ret;
}

static inline void set_already_set(uint64_t* already_set, unsigned int i)
{
    already_set[i/64] |= (((uint64_t)1)<<(i%64));
}

//! Helper to get the next not set index in a molecule
//!
//! \private
//! \param already_set Array indicating unset particels
//! \param N Number of particles in molecule
//! \return next index.
int get_next_index(const uint64_t*const already_set,const unsigned int N)
    {
    unsigned int i;
    for(i=0; i < N; i++)
	if(! is_already_set(already_set, i) )
	    break;
    if( i==N )
	return -1;
    return i;
    }

//! Helper to set the neighbour of a particle.
//! \private
//!
//! \param jbead Particle index to set.
//! \param neigh Bonded Neighbour
//! \param bond_type Type of bond to neighbour
//! \param already_set Array of already set particles
//! \param poly Polymer of the particle
//! \param p System
//! \return Errorcode
int set_neighbour(const unsigned int jbead,const Monomer*const neigh,
		  const unsigned int bond_type,uint64_t* already_set,
		  Polymer*const poly,const struct Phase*const p)
    {
#ifndef SOMA_POLY_CHAIN
    if(is_already_set(already_set, jbead))
	return 0;
#endif
    Monomer dx;
    Monomer new_m;
    int move_allowed;

    do{
	soma_scalar_t scale = 1.;
	dx.x=dx.y=dx.z=0;
    switch(bond_type)
	{
	case HARMONICVARIABLESCALE:;
	    scale = p->harmonic_normb_variable_scale;
	    /* intentionally falls through */
	case HARMONIC: ;
	    soma_normal_vector(&(poly->poly_state),p->args.pseudo_random_number_generator_arg
			       , &(dx.x), &(dx.y), &(dx.z));
	    dx.x /= sqrt(2*p->harmonic_normb*scale);
	    dx.y /= sqrt(2*p->harmonic_normb*scale);
	    dx.z /= sqrt(2*p->harmonic_normb*scale);
	    new_m.x = neigh->x+dx.x; new_m.y = neigh->y+dx.y; new_m.z = neigh->z+dx.z;
	    break;
case STIFF:
	default:
#if ( ! defined SPEC_OPENACC ) && ( ! defined SPEC_OPENMP_TARGET )
	fprintf(stderr, "ERROR: %s:%d unknow bond type appeared %d\n",
		__FILE__, __LINE__,bond_type);
#endif
	new_m.x=new_m.y=new_m.z=0; //Shut up compiler warning
	}
    move_allowed = ! possible_move_area51(p,neigh->x,neigh->y,neigh->z,dx.x,dx.y,dx.z,true);
	}while( move_allowed );

    p->allocator->all_Monomer.device_buf[poly->beads+jbead].x = new_m.x;
    p->allocator->all_Monomer.device_buf[poly->beads+jbead].y = new_m.y;
    p->allocator->all_Monomer.device_buf[poly->beads+jbead].z = new_m.z;
    set_already_set(already_set, jbead);

    //recursively add all connected neighbors
    const int start = get_bondlist_offset(
	p->poly_arch[p->poly_type_offset[poly->type] + jbead + 1]);

#ifndef SOMA_POLY_CHAIN
    if(start > 0){
	int i = start;
	unsigned int end;
	do{
	    const uint32_t info = p->poly_arch[i++];
	    end = get_end(info);
	    const unsigned int bond_type = get_bond_type(info);
	    const int offset = get_offset(info);

	    const int neighbour_id = jbead + offset;
	    set_neighbour(neighbour_id,&(p->allocator->all_Monomer.device_buf[poly->beads+jbead]),bond_type,already_set,poly,p);

	    }while( end == 0);
	}
#endif
    return 0;
    }

int generate_new_beads(struct Phase*const p)
    {

    if( fabs(p->harmonic_normb_variable_scale) < 1e-5 )
	{
	fprintf(stderr,"WARNING: p->harmonic_normb_variable_scale < 1e-5, this may result in unreasonable generated position or even causes divisions by 0.\n");
	}

    uint64_t* already_set_all = 0;
    size_t* already_set_poly_offset = malloc(p->n_polymers * sizeof(size_t));
    size_t maxN = 0, sum = 0;
    for( uint64_t i= 0; i < p->n_polymers; i++)
    {
	already_set_poly_offset[i] = sum;
	unsigned int N = p->poly_arch[ p->poly_type_offset[p->polymers[i].type] ];
	N = (N/64 + ((N%64) > 0));
	if(N > maxN) maxN = N;
	sum += N;
    }
#if (defined SPEC_OPENMP ) || (defined SPEC_OPENACC ) || (defined SPEC_OPENMP_TARGET )
    already_set_all = malloc( sum*sizeof(uint64_t));
    if(already_set_all == NULL)
	{
	fprintf(stderr,"ERROR: %s:%d Malloc problem.\n",__FILE__,__LINE__);
	return -1;
	}
    memset(already_set_all,0,sum*sizeof(uint64_t));
#else
    already_set_all = malloc( maxN*sizeof(uint64_t) );
#endif

#ifdef SPEC_OPENACC
#pragma acc enter data copyin(already_set_all[0:sum],already_set_poly_offset[0:p->n_polymers])
#pragma acc parallel loop independent present(p,already_set_poly_offset,already_set_all)
#endif
#ifdef SPEC_OPENMP_TARGET
#pragma omp target teams distribute parallel for \
    map(to:already_set_all[0:sum],already_set_poly_offset[0:p->n_polymers])
#endif
#ifdef SPEC_OPENMP
#pragma omp parallel for
#endif
    for( uint64_t i= 0; i < p->n_polymers; i++)
	{
	Polymer *const poly = &(p->polymers[i]);
	const unsigned int N = p->poly_arch[ p->poly_type_offset[poly->type] ];
#if (defined SPEC_OPENMP ) || (defined SPEC_OPENACC ) || (defined SPEC_OPENMP_TARGET )
	uint64_t *already_set = &already_set_all[already_set_poly_offset[i]];
#else
	uint64_t *already_set = already_set_all;
	memset(already_set,0,maxN*sizeof(uint64_t));
#endif

	int free_index;
	while( (free_index=get_next_index(already_set, N)) >= 0 )
	{
	    //Set a free monomer
	    soma_scalar_t x,y,z;
	    do{
		x = soma_rng_soma_scalar(&(poly->poly_state),p->args.pseudo_random_number_generator_arg)*p->Lx;
		y = soma_rng_soma_scalar(&(poly->poly_state),p->args.pseudo_random_number_generator_arg)*p->Ly;
		z = soma_rng_soma_scalar(&(poly->poly_state),p->args.pseudo_random_number_generator_arg)*p->Lz;
		}while( p->area51 != NULL &&
			p->area51[coord_to_index(p, x, y, z)] == 1);

	    p->allocator->all_Monomer.device_buf[poly->beads+free_index].x = x;
	    p->allocator->all_Monomer.device_buf[poly->beads+free_index].y = y;
	    p->allocator->all_Monomer.device_buf[poly->beads+free_index].z = z;
	    set_already_set(already_set, free_index);

#ifdef SOMA_POLY_CHAIN
	    // only supporting linear chain in favor of offloading: set first free neighor
	    while(free_index != -1)
	    {
		const int start = get_bondlist_offset(
		    p->poly_arch[p->poly_type_offset[poly->type] + free_index + 1]);
		if(start > 0)
		{
		    int i = start;
		    unsigned int end;
		    int index = -1;

		    do
		    {
			const int info = p->poly_arch[i++];
			end = get_end(info);
			const unsigned int bond_type = get_bond_type(info);
			const int offset = get_offset(info);

			const int neighbour_id = free_index + offset;
			const unsigned int jbead = neighbour_id;

			if(!is_already_set(already_set, jbead))
			{
			    index = neighbour_id;
			    set_neighbour(jbead,&(p->allocator->all_Monomer.device_buf[poly->beads+free_index]),bond_type,already_set,poly,p);
			    /* set_neighbour(jbead,&tmp,bond_type,already_set,poly,p); */
			    end = 1;
			}
		    } while( end == 0);
		    free_index = index;
		}
	    }
#else
	    const int start = get_bondlist_offset(
		p->poly_arch[p->poly_type_offset[poly->type] + free_index + 1]);
	    if(start > 0){
		int i = start;
		unsigned int end;
		do{
		    const uint32_t info = p->poly_arch[i++];
		    end = get_end(info);
		    const unsigned int bond_type = get_bond_type(info);
		    const int offset = get_offset(info);

		    const int neighbour_id = free_index + offset;
		    const unsigned int jbead = neighbour_id;

		    set_neighbour(jbead,&(p->allocator->all_Monomer.device_buf[poly->beads+free_index]),bond_type,already_set,poly,p);

		    }while( end == 0);
		}
#endif
	}

	}

#ifdef SPEC_OPENACC
#ifndef SPEC_NO_VAR_ARRAY_REDUCE
    acc_memcpy_device(p->allocator->all_Monomer_msd.device_buf, p->allocator->all_Monomer.device_buf
	    , p->allocator->all_Monomer_msd.size*sizeof(Monomer)
	);
#else
    acc_memcpy_from_device(p->allocator->all_Monomer_msd.buf, p->allocator->all_Monomer.device_buf
	    , p->allocator->all_Monomer_msd.size*sizeof(Monomer)
	);
#endif
#elif defined SPEC_OPENMP_TARGET
#ifndef SPEC_NO_VAR_ARRAY_REDUCE
    omp_target_memcpy(p->allocator->all_Monomer_msd.device_buf, p->allocator->all_Monomer.device_buf
	    , p->allocator->all_Monomer_msd.size*sizeof(Monomer), 0, 0
	    , soma_get_device(), soma_get_device()
	);
#else
    omp_target_memcpy(p->allocator->all_Monomer_msd.buf, p->allocator->all_Monomer.device_buf
	    , p->allocator->all_Monomer_msd.size*sizeof(Monomer), 0, 0
	    , omp_get_initial_device(), soma_get_device()
	);
#endif
#else
    memcpy( p->allocator->all_Monomer_msd.buf, p->allocator->all_Monomer.buf
	, p->allocator->all_Monomer.size*sizeof(Monomer));
#endif

    free(already_set_all);
    free(already_set_poly_offset);
    update_density_fields(p);
    memcpy(p->old_fields_unified, p->fields_unified, p->n_cells*p->n_types*sizeof(uint16_t));
    return 0;
    }

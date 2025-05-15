/* Copyright (C) 2016-2017 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016-2017 Marcel Langenberg

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

/* */

//! \file io.c
//! \brief Implementation of io.h

#include "io.h"
#include <stdio.h>
#include <assert.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#ifdef SPEC_OPENMP
#include <omp.h>
#endif//SPEC_OPENMP
#include "mesh.h"
#include "ana.h"
#include "allocator.h"

int init_cripple(struct Phase*p,const unsigned int n)
    {
    p->reference_Nbeads = 64;
    const unsigned int NmonoA=p->reference_Nbeads;
    const uint64_t nN = p->reference_Nbeads *n;
    const unsigned int particlesPerCell = 35;
    if( nN*p->info_MPI.Ncores < 6*6*6*particlesPerCell )
    	{
    	fprintf(stderr,"Too few polymers %d to init a propersystem. %s:%d\n!",n,__FILE__,__LINE__);
    	return -1;
    	}
    p->nx = pow( nN / (soma_scalar_t) particlesPerCell, 1/3.);
    if( p->nx < 6)
	p->nx = 6;
    p->Lx = p->nx/6.;
    p->ny = p->nx;
    p->Ly = p->ny/6.;
    p->nz = p->nx;
    p->Lz = p->nz/6.;

    p->n_cells = p->nx * p->ny * p->nz;
    p->time = 0;

    p->n_polymers = n/p->info_MPI.Ncores;

    if( (unsigned int) p->info_MPI.current_core < n % p->info_MPI.Ncores )
      p->n_polymers += 1;

    p->n_polymers_storage = p->n_polymers;
#ifndef NDEBUG
    /* MPI_Allreduce( &(p->n_polymers), &(p->n_polymers_global), 1 , MPI_UNSIGNED, MPI_SUM,p->info_MPI.SOMA_MPI_Comm); */
    /* assert( n == p->n_polymers_global ); */
#endif//NDEBUG

    p->n_polymers_global = n;
    p->n_types = 1;

    const soma_scalar_t dt = 0.17;
    /* Assuming same diffusion constant for all monomers for this input file, but value can vary! */
    p->A = (soma_scalar_t*) malloc(p->n_types * sizeof(soma_scalar_t));
    if(p->A == NULL){
	fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
	return -4;
	}
    p->A[0] = dt / p->reference_Nbeads;

    /* allocate XN interaction matrix */
    p->xn = (soma_scalar_t *) malloc(p->n_types * p->n_types * sizeof(soma_scalar_t));
    if(p->xn == NULL){
	fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
	return -4;
	}

    /* set parameters */
    p->xn[0] = 30;

        //Create the polymer architecture list.
    //The old formate supports only 1 architecture -- a linear chain.
    p->n_poly_type= 1;
    p->poly_type_offset = (int*) malloc( p->n_poly_type*sizeof(int));
    if(p->poly_type_offset == NULL){
	fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
	return -4;
	}
    p->poly_type_offset[0] = 0;
    //Length = 1(N) + N(mono) + 4 (+1,-1,+1-1)
    p->poly_arch_length = 1 + NmonoA + 4;

    p->poly_arch = (uint32_t*) malloc( p->poly_arch_length * sizeof(uint32_t));
    if(p->poly_arch == NULL){
	fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
	return -4;
	}

    //Set up the architecture
    p->poly_arch[0] = NmonoA;

    p->poly_arch[ 0 + 1 ] = get_info_bl(NmonoA +4,0 );
    for(unsigned int i=1; i < NmonoA-1; i++)
	{
	p->poly_arch[ i + 1] = get_info_bl( NmonoA+3, 0);
	}
    p->poly_arch[NmonoA] = get_info_bl( NmonoA +1, 0);

    //Last monomer
    int end,offset,bond_type,info;
    bond_type = HARMONIC;
    end = 1;
    offset = -1;
    info = get_info(offset, bond_type, end);
    p->poly_arch[ NmonoA + 1 ] = info;
    //First monomer
    end = 1;
    offset = 1;
    info = get_info(offset, bond_type, end);
    p->poly_arch[ NmonoA + 2 ] = info;
    //Middle monomers
    end = 0;
    offset = -1;
    info = get_info(offset, bond_type, end);
    p->poly_arch[ NmonoA + 3 ] = info;

    end = 1;
    offset = 1;
    info = get_info(offset, bond_type, end);
    assert( NmonoA + 4 == p->poly_arch_length-1);
    p->poly_arch[ NmonoA + 4 ] = info;

    /* allocate space for polymers */
    p->polymers = (Polymer *) malloc(p->n_polymers_storage * sizeof(Polymer));
    if(p->polymers == NULL){
	fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
	return -4;
	}
    {
	size_t nMonomers = 0;
	for(uint64_t i =0; i < p->n_polymers; i++)
	{
	    p->polymers[i].type=0;
	    nMonomers += p->poly_arch[ p->poly_type_offset[ p->polymers[i].type ] ];
	}
	SET_TYPE_BUF(Monomer, nMonomers);
	SET_TYPE_BUF(Monomer, nMonomers, _msd);
    }
    for(uint64_t i =0; i < p->n_polymers; i++)
	{
	const unsigned int N = p->poly_arch[ p->poly_type_offset[ p->polymers[i].type ] ];
	/* p->polymers[i].beads = (Monomer*) malloc( N*sizeof(Monomer) ); */
	p->polymers[i].beads = alloc_Monomer( N );
	/* p->polymers[i].msd_beads = (Monomer*) malloc( N*sizeof(Monomer) ); */
	p->polymers[i].msd_beads = alloc_Monomer_msd( N );
	assert(p->polymers[i].beads == p->polymers[i].msd_beads);
	}
    p->bead_data_read = false;
    p->area51 = NULL;
    p->external_field_unified = NULL;
    p->hamiltonian = SCMF0;
    p->harmonic_normb_variable_scale = 1;
    p->cm_a = NULL; // Deactivate CM movement with the old file format.

    return 0;
    }

int screen_output(struct Phase*const p,const unsigned int Nsteps)
    {
    static time_t last_print = 0;
    static unsigned int last_time = 0;
    static double last_sec = 0;
    if(last_time == 0) last_time = p->start_time;
    if(last_sec == 0) last_sec = MPI_Wtime();

    const time_t now = time(NULL); if(last_print == 0) last_print = now;

    const double second = difftime(now,p->start_clock);
    const unsigned int steps_done = p->time - p->start_time;
    const unsigned int steps_to_go = Nsteps-steps_done;

    const double now_sec = MPI_Wtime();
    const double sec =  now_sec - last_sec;
    p->tps_elapsed_time += sec;
    p->tps_elapsed_steps += p->time-last_time;

    if( p->args.screen_output_interval_arg > 0 &&  now - last_print >= p->args.screen_output_interval_arg)
	{
	const double tps = p->tps_elapsed_steps / p->tps_elapsed_time;
	p->tps_elapsed_time = 1./tps;
	p->tps_elapsed_steps = 1;

	struct tm future = * localtime( &now );
	future.tm_sec += steps_to_go/tps;
	const time_t end = mktime( &future ) ;

	if(p->info_MPI.current_core == 0)
	    {
	    fprintf(stdout,"Running for %.0f [s] | TPS %g | steps-to-go: %u | ETA: %s",
		    second,tps,steps_to_go,ctime(&end));
	    fflush(stdout);
	    }
	last_print = now;
	last_time = p->time;
	last_sec = now_sec;
	}
    return 0;
    }

#define MAX_FILENAME_LENGTH 1024
int write_state(struct Phase*const p)
    {
    char filename[MAX_FILENAME_LENGTH];
    if( sprintf(filename,"stateR%dtime%dNpoly%d.txt",p->args.rng_seed_arg,p->time,p->args.npoly_arg) > MAX_FILENAME_LENGTH)
	{
	fprintf(stderr, "ERROR invalid filename length %s:%d\n", __FILE__,__LINE__);
	return -1;
	}

    FILE*f = NULL;
    if( p->info_MPI.current_core == 0)
	{
	f = fopen(filename, "w");
	if( f == NULL)
	    {
	    fprintf(stderr, "ERROR cannot open file to write state file.\n");
	    return -1;
	    }
	}

    update_omega_fields(p);
#if ! ( ( defined SPEC_OPENACC && ( ! defined SPEC_NO_VAR_ARRAY_REDUCE ) ) \
	|| ( defined SPEC_OPENMP_TARGET && ( ! defined SPEC_NO_VAR_ARRAY_REDUCE ) ) )
    update_self_phase(p);
#endif

    //Rg
    soma_scalar_t*const Rg=(soma_scalar_t*)malloc(4*p->n_poly_type*sizeof(soma_scalar_t));
    if(Rg == NULL){fprintf(stderr,"ERROR: Malloc %s:%d \n",__FILE__,__LINE__);return -2;}
    calc_Rg(p,Rg);
    if( f != NULL)
	fprintf(f,"Rg: %f %f %f %f\n",Rg[0],Rg[1],Rg[2],Rg[3]);
    free(Rg);


    //Bond anisotropy
    soma_scalar_t*const a=(soma_scalar_t*)malloc(6*p->n_poly_type*sizeof(soma_scalar_t));
    if(a == NULL){fprintf(stderr,"ERROR: Malloc %s:%d \n",__FILE__,__LINE__);return -2;}
    calc_anisotropy(p, a);
    if (p->info_MPI.current_core == 0)
	fprintf(f,"ba: %f %f %f %f %f %f \n",a[0],a[1],a[2],a[3],a[4],a[5]);
    free(a);

    //MSD
    soma_scalar_t*const msd = (soma_scalar_t*)malloc(8*p->n_poly_type*sizeof(soma_scalar_t));
    if(msd == NULL){fprintf(stderr,"ERROR: Malloc %s:%d \n",__FILE__,__LINE__);return -2;}
    calc_MSD(p, msd);
    if (p->info_MPI.current_core == 0)
	fprintf(f,"msd: %f %f %f %f %f %f %f %f \n",msd[0],msd[1],msd[2],msd[3],msd[4],msd[5],msd[6],msd[7]);
    free(msd);

    // bonded energy calculation
    soma_scalar_t*const b_energy = (soma_scalar_t*)malloc(NUMBER_SOMA_BOND_TYPES*sizeof(soma_scalar_t));
    if(b_energy == NULL){fprintf(stderr,"ERROR: Malloc %s:%d \n",__FILE__,__LINE__);return -2;}
    calc_bonded_energy(p, b_energy);
    if(p->info_MPI.current_core == 0)
	fprintf(f,"BE: %f\n",b_energy[0]);
    free(b_energy);


    if( f!= NULL)
	fclose(f);

    return 0;
    }

int compare_list(soma_scalar_t const * const old,soma_scalar_t const * const new, const unsigned int N,soma_scalar_t epsilon)
    {
    for(unsigned int i=0; i < N; i++)
	{
	const soma_scalar_t diff = fabs(old[i]-new[i]);
	//const soma_scalar_t av = (old[i]+new[i])/2.;
	if( fabs(diff)  > epsilon )
	    return i+1;
	}
    return 0;
    }

int compare_state(struct Phase*const p,soma_scalar_t epsilon)
    {
    char filename[MAX_FILENAME_LENGTH] = {0};
    if( sprintf(filename,"stateR%dtime%dNpoly%d.txt",p->args.rng_seed_arg,p->time,p->args.npoly_arg) > MAX_FILENAME_LENGTH)
	{
	fprintf(stderr, "ERROR invalid filename length %s:%d\n", __FILE__,__LINE__);
	return -1;
	}

    FILE*f = NULL;
    int ret = 0;
    if( p->info_MPI.current_core == 0)
	{
	f = fopen(filename, "r");
	if( f == NULL)
	    {
	    fprintf(stderr, "ERROR cannot open file to read state file.\n");
	    ret |= -42;
	    }
	}

    update_omega_fields(p);
    update_self_phase(p);

    //Rg
    soma_scalar_t*const Rg=(soma_scalar_t*)malloc(4*p->n_poly_type*sizeof(soma_scalar_t));
    if(Rg == NULL){fprintf(stderr,"ERROR: Malloc %s:%d \n",__FILE__,__LINE__);return -2;}
    soma_scalar_t*const Rg_old=(soma_scalar_t*)malloc(4*p->n_poly_type*sizeof(soma_scalar_t));
    if(Rg_old == NULL){fprintf(stderr,"ERROR: Malloc %s:%d \n",__FILE__,__LINE__);return -2;}
    calc_Rg(p,Rg);

    if( f != NULL && ret == 0)
	{
	if( fscanf(f,"Rg: %lf %lf %lf %lf\n",Rg_old +0 ,Rg_old +1 ,Rg_old + 2 ,Rg_old + 3) != 4 )
	    {
	    fprintf(stderr,"Improperly configured state file %s:%d %s\n",__FILE__,__LINE__,filename);
	    ret |= -1;
	    }
	int item = compare_list(Rg_old,Rg,4,epsilon);
	if( item != 0 )
	    {
	    fprintf(stderr, "Failed to match Rg state @ item %d %f %f\n",item,Rg[item-1],Rg_old[item-1]);
	    ret |= 1;
	    }
	}
    free(Rg);
    free(Rg_old);


    //Bond anisotropy
    soma_scalar_t*const a=(soma_scalar_t*)malloc(6*p->n_poly_type*sizeof(soma_scalar_t));
    if(a == NULL){fprintf(stderr,"ERROR: Malloc %s:%d \n",__FILE__,__LINE__);return -2;}
    soma_scalar_t*const a_old=(soma_scalar_t*)malloc(6*p->n_poly_type*sizeof(soma_scalar_t));
    if(a_old == NULL){fprintf(stderr,"ERROR: Malloc %s:%d \n",__FILE__,__LINE__);return -2;}

    calc_anisotropy(p, a);
    if( f!= NULL && ret == 0)
	{
	if( fscanf(f,"ba: %lf %lf %lf %lf %lf %lf \n",a_old+0,a_old+1,a_old+2,a_old+3,a_old+4,a_old+5) != 6)
	    {
	    fprintf(stderr,"Improperly configured state file %s:%d %s\n",__FILE__,__LINE__,filename);
	    ret |= -1;
	    }
	int item = compare_list(a_old,a,6,epsilon);
	if( item != 0 )
	    {
	    fprintf(stderr, "Failed to match ba state @ item %d\n",item);
	    ret |= 1;
	    }
	}
    free(a);
    free(a_old);

    //MSD
    soma_scalar_t*const msd = (soma_scalar_t*)malloc(8*p->n_poly_type*sizeof(soma_scalar_t));
    if(msd == NULL){fprintf(stderr,"ERROR: Malloc %s:%d \n",__FILE__,__LINE__);return -2;}
    soma_scalar_t*const msd_old = (soma_scalar_t*)malloc(8*p->n_poly_type*sizeof(soma_scalar_t));
    if(msd_old == NULL){fprintf(stderr,"ERROR: Malloc %s:%d \n",__FILE__,__LINE__);return -2;}

    calc_MSD(p, msd);
    if( f!= NULL && ret == 0)
	{
	if( fscanf(f,"msd: %lf %lf %lf %lf %lf %lf %lf %lf \n",msd_old+0,msd_old+1,msd_old+2,msd_old+3,msd_old+4,msd_old+5,msd_old+6,msd_old+7) != 8)
	    {
	    fprintf(stderr,"Improperly configured state file %s:%d %s\n",__FILE__,__LINE__,filename);
	    ret |= -1;
	    }
	int item = compare_list(msd_old,msd,8,epsilon);
	if( item != 0 )
	    {
	    fprintf(stderr, "Failed to match msd state @ item %d %f %f\n",item,msd_old[item-1],msd[item-1]);
	    ret |= 1;
	    }
	}
    free(msd);
    free(msd_old);

    // bonded energy calculation
    soma_scalar_t*const b_energy = (soma_scalar_t*)malloc(NUMBER_SOMA_BOND_TYPES*sizeof(soma_scalar_t));
    if(b_energy == NULL){fprintf(stderr,"ERROR: Malloc %s:%d \n",__FILE__,__LINE__);return -2;}
    soma_scalar_t*const b_energy_old = (soma_scalar_t*)malloc(NUMBER_SOMA_BOND_TYPES*sizeof(soma_scalar_t));
    if(b_energy_old == NULL){fprintf(stderr,"ERROR: Malloc %s:%d \n",__FILE__,__LINE__);return -2;}

    calc_bonded_energy(p, b_energy);
    if( f!= NULL && ret == 0)
	{
	if( fscanf(f,"BE: %lf \n",b_energy_old+0) != 1)
	    {
	    fprintf(stderr,"Improperly configured state file %s:%d %s\n",__FILE__,__LINE__,filename);
	    ret |= -1;
	    }
	int item = compare_list(b_energy_old,b_energy,1,epsilon);
	if( item != 0 )
	    {
	    fprintf(stderr, "Failed to match BE state @ item %d %f %f\n",item,b_energy[0],b_energy_old[0]);
	    ret |= 1;
	    }
	}
    free(b_energy);
    free(b_energy_old);

    if( f!= NULL)
	{
	if( ret == 0)
	    printf("State comparison successfull.\n");
	fclose(f);
	}


    return ret;
    }

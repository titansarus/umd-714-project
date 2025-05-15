/* Copyright (C) 2016 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016 Marcel Langenberg
   Copyright (C) 2016 Fabien Leonforte
   Copyright (C) 2016 Juan Orozco
   Copyright (C) 2016 Yongzhi Ren

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

/* SOMA  */
/* TO do Typedef for curand */
/* Clean up the hdf5 output*/

//! \file soma.c
//! \brief Implementation of the main executable SOMA.

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "io.h"
#include "mc.h"
#include "mesh.h"
#include "mpiroutines.h"
#include "init.h"
#include "test.h"
#include "ana.h"
#include "signal.h"
#include "rng.h"
#include "generate_positions.h"
#include "allocator.h"
#include "device.h"

#ifdef SPEC
#include "specmpitime.h"
#endif

#include "mpi_timing.h"

double total_start_time, total_end_time, total_execution_time;

//! Main Function of the Executable SOMA
//! \private
//!
//! \param argc Argument counter
//! \param argv Argument vector
//! \return Errorcode
int main(int argc, char *argv[])
{

    if (MPI_Init(NULL, NULL) != MPI_SUCCESS) {
      fprintf(stderr, "MPI_ERROR (start)\n");
      return -1;
    }
SPEC_TIME_START
SPEC_TIME_START_INIT
    
    // Initialize MPI timing system
    init_mpi_timing();

    // Start total execution timer
    total_start_time = MPI_Wtime();

    Allocator allocator;
    global_allocator = &allocator;
    init_allocator(&allocator);

    Phase phase;
    Phase *const p = &phase;
    p->allocator = &allocator;
    //Set the communicator for this process to world.
    p->info_MPI.SOMA_MPI_Comm = MPI_COMM_WORLD;
    /* initialize MPI */
    init_MPI(p);

    const int args_success = cmdline_parser(argc,argv,&(p->args));
    if( args_success < 0 )
	fprintf(stderr,"Process %d failed to read the cmdline. Exiting.\n",p->info_MPI.current_core);
    const int mpi_args = check_status_on_mpi(p,args_success);

    if( mpi_args != 0 )
	{
	finalize_MPI();
	return 0;
	}
    post_process_args( &(p->args), &(p->info_MPI));
    if(p->args.move_type_arg == move_type_arg_TRIAL)
	{MPI_ERROR_CHECK(1,"ERROR: Trial move type is currently not working.");}

    if(p->args.pseudo_random_number_generator_arg == pseudo_random_number_generator_arg_TT800
	    || p->args.pseudo_random_number_generator_arg == pseudo_random_number_generator_arg_MT)
	{MPI_ERROR_CHECK(1,"ERROR: TT800 PRNG is currently not working.");}

#ifdef SPEC_OPENACC
    const int open_acc = set_openacc_devices(p);
    if( check_status_on_mpi(p,open_acc) != 0){
	if(p->info_MPI.current_core == 0)
	    fprintf(stderr,"ERROR: cannot set openacc devices.\n");
	finalize_MPI();
	return 1;
	}
#endif
#ifdef SPEC_OPENMP_TARGET
    const int open_mp = set_openmp_devices(p);
    if( check_status_on_mpi(p,open_mp) != 0){
	if(p->info_MPI.current_core == 0)
	    fprintf(stderr,"ERROR: cannot set openmp devices.\n");
	finalize_MPI();
	return 1;
	}
#endif

    const unsigned int N_steps = p->args.timesteps_arg;
    /* read in configuration with routine from io */
    const int read = init_cripple(p,p->args.npoly_arg);
    MPI_ERROR_CHECK(read, "Cannot init cripple conf.");

    /* initialize phase */
    const int init = init_phase(p);
    MPI_ERROR_CHECK(init, "Cannot init values.");

    if( !p->bead_data_read )
	{
	const int new_beads = generate_new_beads(p);
	MPI_ERROR_CHECK(new_beads, "Cannot genrate new bead data.");
	}

    if( ! p->args.skip_tests_flag)
	{
	const int test_p = test_particle_types(p);
	MPI_ERROR_CHECK(test_p, "Partile type test failed.");

	const int test51 = test_area51_violation(p);
	MPI_ERROR_CHECK(test51, "Area51 test failed.");

	const int test51_exact = test_area51_exact(p);
	if( ! p->args.nonexact_area51_flag )
	    MPI_ERROR_CHECK(test51_exact, "Area51 exact test failed.");

	const int indepent_sets = test_independet_sets(p);
	MPI_ERROR_CHECK(indepent_sets, "Indepent Set test failed.");
	}

    //Reset the RNG to initial starting conditions.
    reseed(p, p->args.rng_seed_arg);
SPEC_TIME_STOP_INIT
SPEC_TIME_START_COMP
    for (unsigned int i = 0; i < N_steps; i++) {
	const int mc_error = monte_carlo_propagation(p, 1);
        if( mc_error != 0)
	    fprintf(stderr,"ERROR %d in monte_carlo_propagation on rank %d.\n"
		    ,mc_error,p->info_MPI.current_core);

	screen_output(p,N_steps);
    }
SPEC_TIME_STOP_COMP

    if( ! p->args.skip_tests_flag)
	{
	const int test51 = test_area51_violation(p);
	MPI_ERROR_CHECK(test51, "Area51 test failed.");
	const int test51_exact = test_area51_exact(p);
	if(! p->args.nonexact_area51_flag )
	    MPI_ERROR_CHECK(test51_exact, "Area51 exact test failed.");
	}

    int ret =0;
    if( p->args.gen_state_file_flag )
	    {
	    if( write_state(p) != 0)
	        {
	        fprintf(stderr,"ERROR: Writing state file.\n");
	        ret = -1;
	        }
	    else
	        {
	        ret = 0;
	        printf("State file generated\n");
	        }
	    }
    else
	    {
	    if(compare_state(p,1e-5) != 0)
	       {
	       fprintf(stderr, "Unable to verify state. FAILURE.\n");
	       ret = -11;
	       }
	    }


    /* deallocate all memory */
    free_phase(p);

    printf("Rank: %d \t polymers %ld\n",p->info_MPI.current_core, p->n_polymers);
    
    // Stop total execution timer
    total_end_time = MPI_Wtime();
    total_execution_time = total_end_time - total_start_time;

    // Report MPI timing statistics
    finalize_mpi_timing(p->info_MPI.current_core, total_execution_time);

SPEC_TIME_STOP
SPEC_TIME_FINAL(true,0.0)
    finalize_MPI();
    if(p->info_MPI.current_core == 0 && !( ret == 0 || ret == 1) )
	printf("Execution finished with errors.\n");

    free_allocator(&allocator);
    return ret;
}

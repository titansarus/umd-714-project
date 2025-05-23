/*---------------------------------------------------------------------------*/
/*!
 * \file   sweep.c
 * \author Wayne Joubert
 * \date   Wed May 22 11:22:14 EDT 2013
 * \brief  Main driver for sweep miniapp.
 * \note   Copyright (C) 2013 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arguments.h"
#include "env.h"
#include "definitions.h"
#include "dimensions.h"
#include "pointer.h"
#include "quantities.h"
#include "array_operations.h"
#include "sweeper.h"

#include "runner.h"

#ifdef SPEC
#include "specmpitime.h"
#endif

/*===========================================================================*/
/*---Main---*/

int main( int argc, char** argv )
{
  /*---Declarations---*/
  Env env = Env_null();

  /*---Initialize for execution---*/

  Env_initialize( &env, argc, argv );
SPEC_TIME_START
SPEC_TIME_START_INIT
  
  Arguments args = Arguments_null();
  Runner runner = Runner_null();

  Arguments_create( &args, argc, argv );
  Runner_create( &runner );

  Env_set_values( &env, &args );

  /*---Perform run---*/
SPEC_TIME_STOP_INIT
SPEC_TIME_START_COMP
  if( Env_is_proc_active( &env ) )
  {
    Runner_run_case( &runner, &args, &env );
  }
  /* if( Env_is_proc_master( &env ) ) */
  /* { */
  /*   printf( "Normsq result: %.8e  diff: %.3e  %s  time: %.3f  GF/s: %.3f\n", */
  /*           (double)runner.normsq, (double)runner.normsqdiff, */
  /*           runner.normsqdiff==P_zero() ? "PASS" : "FAIL", */
  /*           (double)runner.time, runner.floprate ); */
  /*   /\*---If invoked with no arguments as part of tester, then ouptut */
  /*        pass/fail count banner to be parsed by testing script---*\/ */
  /*   if( argc == 1 ) */
  /*   { */
  /*       const int ntest = 1; */
  /*       const int ntest_passed = runner.normsqdiff==P_zero() ? 1 : 0; */
  /*       printf( "TESTS %i    PASSED %i    FAILED %i\n", */
  /*           ntest, ntest_passed, ntest-ntest_passed ); */
  /*   } */
  /* } */

  /*---Deallocations---*/
SPEC_TIME_STOP_COMP
  Runner_destroy( &runner );
  Arguments_destroy( &args );

  /*---Finalize execution---*/
SPEC_TIME_STOP
SPEC_TIME_FINAL(true,0.0)
  Env_finalize( &env );
  //  printf("blah\n");
 
  exit(0);
} /*---main---*/

/*---------------------------------------------------------------------------*/

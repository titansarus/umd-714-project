/*---------------------------------------------------------------------------*/
/*!
 * \file   array_operations.c
 * \author Wayne Joubert
 * \date   Thu Jan 16 15:39:53 EST 2014
 * \brief  Functions to operate on multi-dim arrays.
 * \note   Copyright (C) 2014 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include "env.h"
#include "definitions.h"
#include "dimensions.h"
#include "array_accessors.h"
#include "quantities.h"
#include "array_operations.h"
#include <stdio.h>

#ifdef __cplusplus
extern "C"
{
#endif

/*===========================================================================*/
/*---Initialize state vector to required input value---*/

void initialize_state( P* const __restrict__   v,
                       const Dimensions        dims,
                       const int               nu,
                       const Quantities* const quan )
{
  int ix = 0;
  int iy = 0;
  int iz = 0;
  int ie = 0;
  int im = 0;
  int iu = 0;

  for( iz=0; iz<dims.ncell_z; ++iz )
  for( iy=0; iy<dims.ncell_y; ++iy )
  for( ix=0; ix<dims.ncell_x; ++ix )
  for( ie=0; ie<dims.ne; ++ie )
  for( im=0; im<dims.nm; ++im )
  for( iu=0; iu<nu; ++iu )
  {
    *ref_state( v, dims, nu, ix, iy, iz, ie, im, iu )
                = Quantities_init_state( quan, ix, iy, iz, ie, im, iu, dims );
  }
}

/*===========================================================================*/
/*---Initialize state vector to zero---*/

void initialize_state_zero( P* const __restrict__ v,
                            const Dimensions      dims,
                            const int             nu )
{
  size_t i = 0;
  size_t n = Dimensions_size_state( dims, nu );

  for( i=0; i<n; ++i )
  {
    v[i] = P_zero();
  }
}

/*===========================================================================*/
/*---Compute vector norm info for state vector---*/

void get_state_norms( const P* const __restrict__ vi,
                      const P* const __restrict__ vo,
		      double* rank_totals,
                      const Dimensions            dims,
                      const int                   nu,
                      P* const __restrict__       normsqp,
                      P* const __restrict__       normsqdiffp,
                      Env* const                  env )
{
  Assert( normsqp     != NULL ? "Null pointer encountered" : 0 );
  Assert( normsqdiffp != NULL ? "Null pointer encountered" : 0 );

  int ix = 0;
  int iy = 0;
  int iz = 0;
  int ie = 0;
  int im = 0;
  int iu = 0;

  P normsq     = P_zero();
  P normsqdiff = P_zero();

  P local_total = 0.0;;  
  P* local_totals = (P*)malloc(dims.ne * sizeof(P));
  int i;
  for (i=0; i<dims.ne; i++) {
    local_totals[i] = 0.0;
  }
  
/* #ifdef USE_MPI */
/*   int world_size; */
/*   MPI_Comm_size(MPI_COMM_WORLD, &world_size); */
/*   P* local_totals = (P*)malloc(world_size * sizeof(P)); */
/*   int i; */
/*   for (i=0; i<world_size; i++){ */
/*     local_totals[i] = 0.0; */
/*   } */
/* #endif */
  
  for( iz=0; iz<dims.ncell_z; ++iz )
  for( iy=0; iy<dims.ncell_y; ++iy )
  for( ix=0; ix<dims.ncell_x; ++ix )
  for( ie=0; ie<dims.ne; ++ie )
    {
  for( im=0; im<dims.nm; ++im )
  for( iu=0; iu<nu; ++iu )
  {
    const P val_vi = *const_ref_state( vi, dims, nu, ix, iy, iz, ie, im, iu );
    const P val_vo = *const_ref_state( vo, dims, nu, ix, iy, iz, ie, im, iu );
    const P diff   = val_vi - val_vo;
    normsq        += val_vo * val_vo;
    normsqdiff    += diff   * diff;
#ifdef USE_MPI
  local_total += val_vo;
#endif
  }
  local_totals[ie] += local_total;
    }
  
  Assert( normsq     >= P_zero() );
  Assert( normsqdiff >= P_zero() );
  normsq     = Env_sum_P( env, normsq );
  normsqdiff = Env_sum_P( env, normsqdiff );

  /*--- sum totals ---*/
#ifdef USE_MPI
  const int mpi_code = MPI_Reduce(local_totals, rank_totals, dims.ne, MPI_DOUBLE, MPI_SUM, 0, Env_mpi_active_comm_( env ));
  MPI_Barrier(Env_mpi_active_comm_(env));
  Assert( mpi_code == MPI_SUCCESS);
#endif

  //free(local_totals);
  *normsqp     = normsq;
  *normsqdiffp = normsqdiff;
}

/*===========================================================================*/
/*---Copy vector---*/

void copy_vector(       P* const __restrict__ vo,
                  const P* const __restrict__ vi,
                  const size_t                n )
{
  Assert( n+1 >= 1 );
  size_t i = 0;

  for( i=0; i<n; ++i )
  {
    vo[i] = vi[i];
  }
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/

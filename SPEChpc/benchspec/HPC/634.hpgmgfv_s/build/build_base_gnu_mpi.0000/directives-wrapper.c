#if defined(SPEC_OPENMP_TARGET) || defined(SPEC_OPENACC)

# include "directives/operators-directives.fv4.c"

#else

# ifdef USE_MPI
#  include <mpi.h>
# endif
# include <stdlib.h>
# include <stdio.h>
# include "level.h"
# include "directives/common.h"

void directives_stub_error(void);

#define FUNCTION_SHOULD_NOT_BE_USED \
  fprintf(stderr, "\nError! The directives stub '%s' should not be used. Re-compile using SPEC_OPENMP_TARGET or SPEC_OPENACC\n", __func__); \
  directives_stub_error();

void directives_stub_error(void)
{
#ifdef USE_MPI
  MPI_Abort(MPI_COMM_WORLD, -1);
#endif
  exit(-1);
}


/* directives/misc.h */

void directives_zero_vector(level_type *d_level, int id)
{
  FUNCTION_SHOULD_NOT_BE_USED
}

void directives_scale_vector(level_type *d_level, int id_c, double scale_a, int id_a)
{
  FUNCTION_SHOULD_NOT_BE_USED
}

void directives_shift_vector(level_type *d_level, int id_c, double shift_a, int id_a)
{
  FUNCTION_SHOULD_NOT_BE_USED
}

void directives_mul_vectors(level_type *d_level, int id_c, double scale, int id_a, int id_b)
{
  FUNCTION_SHOULD_NOT_BE_USED
}

void directives_add_vectors(level_type *d_level, int id_c, double scale_a, int id_a, double scale_b, int id_b)
{
  FUNCTION_SHOULD_NOT_BE_USED
}

double directives_sum(level_type *d_level, int id)
{
  FUNCTION_SHOULD_NOT_BE_USED
  return -0.0;
}

double directives_max_abs(level_type *d_level, int id)
{
  FUNCTION_SHOULD_NOT_BE_USED
  return -0.0;
}

void directives_color_vector(level_type *d_level, int id_a, int colors_in_each_dim, int icolor, int jcolor, int kcolor)
{
  FUNCTION_SHOULD_NOT_BE_USED
}


/* directives/blockCopy.h */

void directives_copy_block(level_type *level, int id, communicator_type *exchange_ghosts, int block_type)
{
  FUNCTION_SHOULD_NOT_BE_USED
}

void directives_increment_block(level_type *level, int id, double prescale, communicator_type *exchange_ghosts, int block_type)
{
  FUNCTION_SHOULD_NOT_BE_USED
}


/* directives/boundary_fv.h */

void directives_apply_BCs_v1(level_type *level, int x_id, int shape)
{
  FUNCTION_SHOULD_NOT_BE_USED
}

void directives_apply_BCs_v2(level_type *level, int x_id, int shape)
{
  FUNCTION_SHOULD_NOT_BE_USED
}

void directives_apply_BCs_v4(level_type *level, int x_id, int shape)
{
  FUNCTION_SHOULD_NOT_BE_USED
}

void directives_extrapolate_betas(level_type *level, int shape)
{
  FUNCTION_SHOULD_NOT_BE_USED
}


/* directives/stencils/smooth.reg.fv4.h */

void directives_smooth_unoptimized(level_type *level, int x_id, int rhs_id, double a, double b, int s, double *c, double *d)
{
  FUNCTION_SHOULD_NOT_BE_USED
}

void directives_smooth_optimized(level_type *level, int x_id, int rhs_id, double a, double b, int s, double *c, double *d)
{
  FUNCTION_SHOULD_NOT_BE_USED
}

/* directives/stencils/residual.reg.fv4.h */

void directives_residual(level_type *level, int res_id, int x_id, int rhs_id, double a, double b, int REBUILD)
{
  FUNCTION_SHOULD_NOT_BE_USED
}

void directives_rebuild(level_type *level, int x_id, int Aii_id, int sumAbsAij_id, double a, double b)
{
  FUNCTION_SHOULD_NOT_BE_USED
}


/* directives/restriction.h */

void directives_restriction(level_type *level_c, int id_c, level_type *level_f, int id_f, communicator_type *restriction, int restrictionType, int block_type)
{
  FUNCTION_SHOULD_NOT_BE_USED
}


/* directives/interpolation_v2.h */

void directives_interpolation_v2(level_type *level_f, int id_f, double prescale_f, level_type *level_c, int id_c, communicator_type *interpolation, int block_type)
{
  FUNCTION_SHOULD_NOT_BE_USED
}


/* directives/interpolation_v4.h */

void directives_interpolation_v4(level_type *level_f, int id_f, double prescale_f, level_type *level_c, int id_c, communicator_type *interpolation, int block_type)
{
  FUNCTION_SHOULD_NOT_BE_USED
}

#endif

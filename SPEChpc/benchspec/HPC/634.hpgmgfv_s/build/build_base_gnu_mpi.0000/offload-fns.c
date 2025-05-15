#include "offload-fns.h"

#include "cuda/common.h"
#include "directives/common.h"

#if defined(SPEC_CUDA) || defined(SPEC_CUDA_API)
#include <cuda_runtime_api.h>
#endif

/* Use CUDA functions when we are not
   using OpenMP target offload or OpenACC. */
# if defined(SPEC_CUDA) && \
  (!defined(SPEC_OPENMP_TARGET) && !defined(SPEC_OPENACC))
int use_cuda = 1;
# else
int use_cuda = 0;
# endif

/* Fallback to specific CUDA functions if we know
   the OpenMP target offload or OpenACC version of
   the function is incorrect for this compiler */
# if defined(SPEC_CUDA) && \
  (defined(SPEC_OPENMP_TARGET) || defined(SPEC_OPENACC)) && \
  defined(USE_CUDA_FALLBACK)
int use_cuda_fallback = 1;
#else
int use_cuda_fallback = 0;
#endif

#define NUM_DIRECTIVE_FUNCTIONS 20
static int bad_directive_function[NUM_DIRECTIVE_FUNCTIONS];
static const char* device_functions[NUM_DIRECTIVE_FUNCTIONS] = {
  "device_increment_block",
  "device_copy_block",
  "device_extrapolate_betas",
  "device_residual",
  "device_rebuild",
  "device_color_vector",
  "device_smooth",
  "device_restriction",
  "device_apply_BCs_v1",
  "device_apply_BCs_v2",
  "device_apply_BCs_v4",
  "device_interpolation_v4",
  "device_interpolation_v2",
  "device_zero_vector",
  "device_scale_vector",
  "device_shift_vector",
  "device_mul_vectors",
  "device_add_vectors",
  "device_sum",
  "device_max_abs"
};

void set_cuda_execution_array(void)
{
  /* Assume all functions are correct for unknown compiler versions */
  for (int i=0; i<NUM_DIRECTIVE_FUNCTIONS; ++i) {
    bad_directive_function[i] = 0;
  }

#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
# if !defined(__ibmxl__) && defined(__clang__) && (__clang_major__ < 11)
  /* LLVM/Clang compiler < 11 and Cray CCE >= 9.0 compiler

     Functions affected by LLVM/Clang bug 44390
     Includes the open-source LLVM/Clang compiler and Cray CCE >= 9.0 compiler
     directives_residual (for residual - 3 and rebuild - 4)
     directives_smooth_optimized and directives_smooth_unoptimized - 6
     directives_interpolation_v4 - 11
     directives_interpolation_v2 - 12 */
#  ifdef OPENMP_TARGET_NOOPT
  /* This code path assumes files are compiled at -O0 optimization */
  static const int bad_openmp_optimization[NUM_DIRECTIVE_FUNCTIONS] =
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
#  else
  static const int bad_openmp_optimization[NUM_DIRECTIVE_FUNCTIONS] =
    {0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0};
#  endif

  for (int i=0; i<NUM_DIRECTIVE_FUNCTIONS; ++i) {
    bad_directive_function[i] = bad_openmp_optimization[i];
  }

# elif defined(__ibmxl__)
  /* IBM compiler
     Exact result match against CUDA when setting
     -Xllvm2ptx -nvvm-compile-options=-fma=0

     Correct but numerically different values when using
     high optimization - this is just here for reference
     directives_extrapolate_betas - 2
     directives_residual (for residual - 3 and rebuild - 4)
     directives_smooth_optimized and directives_smooth_unoptimized - 6
     directives_apply_BCs_v4 - 10
     directives_interpolation_v4 - 11 */

# elif defined(__PGIC__)
  /* PGI compiler
     Functions giving different values to the CUDA version:
     directives_apply_BCs_v2 - 9
     directives_apply_BCs_v4 - 10

     This is a bug in the switch statement where the case
     variable is off-by one, so value 1 is matching case 2,
     value 3 matching case 4, etc. It is fixed in PGI-20.4
  */
#  if defined(__PGLLVM__) && (__PGIC__ <= 20) && (__PGIC_MINOR__ < 4)
  static const int bad_switch_statement[NUM_DIRECTIVE_FUNCTIONS] =
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
#  else
  static const int bad_switch_statement[NUM_DIRECTIVE_FUNCTIONS] =
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
#  endif

  for (int i=0; i<NUM_DIRECTIVE_FUNCTIONS; ++i) {
    bad_directive_function[i] = bad_switch_statement[i];
  }
# endif

# ifdef VERBOSE
  int rank = 0;
#  ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#  endif
  for (int i=0; i<NUM_DIRECTIVE_FUNCTIONS; ++i) {
    if (rank == 0 && bad_directive_function[i]) {
      fprintf(stderr, "[Warning]: compiler issue in directives version of %s",
	     device_functions[i]);
      if (use_cuda_fallback) {
	fprintf(stderr, "... falling back to CUDA\n");
      } else {
	fprintf(stderr, "... ignoring issue and continuing\n");
      }
    }
  }
# endif
#endif
}

int use_cuda_function(int fid)
{
#ifdef DO_TRACE
  const int do_trace = 1;
#else
  const int do_trace = 0;
#endif
  int rank;
  if (fid >= NUM_DIRECTIVE_FUNCTIONS) {
    printf("Array will be accessed out of bounds\n");
    exit(-1);
  }
  if (do_trace) {
    rank = 0;
#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    printf("[Rank %d]: Executing %s\n", rank, device_functions[fid]);
    fflush(stdout);
  }
  return (use_cuda ||
	  (use_cuda_fallback && bad_directive_function[fid]));
}

void check_consistent_build(void) {
#if (defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)) && !defined(SPEC_CUDA)
  if (use_cuda_fallback == 1) {
    printf("Error: you are using a pure OpenACC/OpenMP-target build but asking to use CUDA functions\n");
    exit(-1);
  }
#endif

#if defined(SPEC_CUDA) && !defined(SPEC_OPENACC) && !defined(SPEC_OPENMP_TARGET)
  if (use_cuda == 0) {
    printf("Error: you are using a pure CUDA build but asking to use OpenACC/OpenMP-target functions\n");
    exit(-1);
  }
#endif
}

//---------------------------------------------------------------------------------------------------------------------------------------------------
void *um_malloc(size_t size, int access_policy)
{
#if defined(CUDA_UM_ALLOC)
  void *ptr;
  switch (access_policy) {
  case UM_ACCESS_GPU:
    CUDA_API_ERROR( cudaMallocManaged(&ptr, size, cudaMemAttachGlobal) )
    break;
  case UM_ACCESS_BOTH:
#ifdef CUDA_UM_ZERO_COPY
    // assumes that the direct access to sysmem is supported on this OS/GPU
    CUDA_API_ERROR( cudaMallocHost(&ptr, size) )
#else
    // default is the managed allocation with global attach
    CUDA_API_ERROR( cudaMallocManaged(&ptr, size, cudaMemAttachGlobal) )
#endif
    break;
  case UM_ACCESS_CPU:
    return malloc(size);
    break;
  }
  return ptr;
#else
  // note that currently regular heap allocations are not accessible by GPU
  return malloc(size);
#endif
}

// allocate data that is not migratable
void *um_malloc_pinned(size_t size, int access_policy)
{
#if defined(SPEC_CUDA) || defined(SPEC_CUDA_API)
  void *ptr;
  switch (access_policy) {
  case UM_ACCESS_GPU:
    CUDA_API_ERROR( cudaMalloc(&ptr, size) )
    break;
  case UM_ACCESS_BOTH:
    CUDA_API_ERROR( cudaMallocHost(&ptr, size) )
    break;
  case UM_ACCESS_CPU:
    return malloc(size);
    break;
  }
  return ptr;
#else
  return malloc(size);
#endif
}

void *um_realloc(void *ptr, size_t size, int access_policy)
{
#if defined(CUDA_UM_ALLOC)
  void *new_ptr;
  switch (access_policy) {
  case UM_ACCESS_GPU:
  case UM_ACCESS_BOTH:
    new_ptr = um_malloc(size, access_policy);
    // realloc always happen from size/2 to size in HPGMG
    CUDA_API_ERROR( cudaMemcpy(new_ptr, ptr, (size/2), cudaMemcpyDefault) )
    um_free(ptr, access_policy);
    break;
  case UM_ACCESS_CPU:
    new_ptr = realloc(ptr, size);
    break;
  }
  return new_ptr;
#else
  // note that currently regular heap allocations are not accessible by GPU
  return realloc(ptr, size);
#endif
}

void um_free(void *ptr, int access_policy)
{
#if defined(CUDA_UM_ALLOC)
  switch(access_policy) {
  case UM_ACCESS_GPU:
    CUDA_API_ERROR( cudaFree(ptr) )
    break;
  case UM_ACCESS_BOTH:
#ifdef CUDA_UM_ZERO_COPY
    CUDA_API_ERROR( cudaFreeHost(ptr) )
#else
    CUDA_API_ERROR( cudaFree(ptr) )
#endif
    break;
  case UM_ACCESS_CPU:
    free(ptr);
    break;
  }
#else
  free(ptr);
#endif
}

void um_free_pinned(void *ptr, int access_policy)
{
#if defined (SPEC_CUDA) || defined(SPEC_CUDA_API)
  switch(access_policy) {
  case UM_ACCESS_GPU:
    CUDA_API_ERROR( cudaFree(ptr) )
    break;
  case UM_ACCESS_BOTH:
    CUDA_API_ERROR( cudaFreeHost(ptr) )
    break;
  case UM_ACCESS_CPU:
    free(ptr);
    break;
  }
#else
  free(ptr);
#endif
}

//---------------------------------------------------------------------------------------------------------------------------------------------------

void do_sync(void) {
  /* This function is called by functions outside of this file */
#if defined(SPEC_CUDA) || defined(SPEC_CUDA_API)
  cudaDeviceSynchronize();
#endif
}

void do_sync_if_mixed(void) {
  /* We do the synchronization if we are mixing CUDA compute
     kernels with either OpenMP target offload or OpenACC */
#if defined(SPEC_CUDA)
# if defined(SPEC_OPENMP_TARGET) || defined(SPEC_OPENACC)
  do_sync();
# endif
#endif
}

int device_runtime_init(void)
{
#if defined(SPEC_OPENACC) || defined(SPEC_OPENMP_TARGET)
  return directives_runtime_init();
#else
  return 0;
#endif
}

void device_increment_block(level_type *level, int id, double prescale, communicator_type *exchange_ghosts, int block_type)
{
  /* This function is only used in multi-rank runs.
     Compile with -DHOST_LEVEL_SIZE_THRESHOLD=1000 to test this function
     when using 2 MPI ranks for ./hpgmg-fv 5 8 1 */
  if (use_cuda_function(0)) {
    cuda_increment_block(*level, id, prescale, *exchange_ghosts, block_type);
  } else {
    do_sync_if_mixed();
    directives_increment_block(level, id, prescale, exchange_ghosts, block_type);
    do_sync_if_mixed();
  }
}

void device_copy_block(level_type *level, int id, communicator_type *exchange_ghosts, int block_type)
{
  if (use_cuda_function(1)) {
    cuda_copy_block(*level, id, *exchange_ghosts, block_type);
  } else {
    do_sync_if_mixed();
    directives_copy_block(level, id, exchange_ghosts, block_type);
    do_sync_if_mixed();
  }
}

void device_extrapolate_betas(level_type *level, int shape)
{
  if (use_cuda_function(2)) {
    cuda_extrapolate_betas(*level, shape);
  } else {
    do_sync_if_mixed();
    directives_extrapolate_betas(level, shape);
    do_sync_if_mixed();
  }
}

void device_residual(level_type *level, int res_id, int x_id, int rhs_id, double a, double b) {
  if (use_cuda_function(3)) {
    cuda_residual(*level, res_id, x_id, rhs_id, a, b);
  } else {
    do_sync_if_mixed();
    int rebuild = 0; // Same function for residual and rebuild computations
    directives_residual(level, res_id, x_id, rhs_id, a, b, rebuild);
    do_sync_if_mixed();
  }
}

void device_rebuild(level_type *level, int x_id, int Aii_id, int sumAbsAij_id, double a, double b) {
  if (use_cuda_function(4)) {
    cuda_rebuild(*level, x_id, Aii_id, sumAbsAij_id, a, b);
  } else {
    do_sync_if_mixed();
    int rebuild = 1; // Same function for residual and rebuild computations
    directives_residual(level, sumAbsAij_id, x_id, Aii_id, a, b, rebuild);
    do_sync_if_mixed();
  }
}

void device_color_vector(level_type *d_level, int id_a, int colors_in_each_dim, int icolor, int jcolor, int kcolor) {
  if (use_cuda_function(5)) {
    cuda_color_vector(*d_level, id_a, colors_in_each_dim, icolor, jcolor, kcolor);
  } else {
    do_sync_if_mixed();
    directives_color_vector(d_level, id_a, colors_in_each_dim, icolor, jcolor, kcolor);
    do_sync_if_mixed();
  }
}

void device_smooth(level_type *level, int x_id, int rhs_id, double a, double b, int s, double *c, double *d) {
  if (use_cuda_function(6)) {
    cuda_smooth(*level, x_id, rhs_id, a, b, s, c, d);
  } else {
    do_sync_if_mixed();
#ifdef USE_SMOOTH_UNOPTIMIZED
    directives_smooth_unoptimized(level, x_id, rhs_id, a, b, s, c, d);
#else
    directives_smooth_optimized(level, x_id, rhs_id, a, b, s, c, d);
#endif
    do_sync_if_mixed();
  }
}

void device_restriction(level_type *level_c, int id_c, level_type *level_f, int id_f, communicator_type *restriction, int restrictionType, int block_type) {
  if (use_cuda_function(7)) {
    cuda_restriction(*level_c, id_c, *level_f, id_f, *restriction, restrictionType, block_type);
  } else {
    do_sync_if_mixed();
    directives_restriction(level_c, id_c, level_f, id_f, restriction, restrictionType, block_type);
    do_sync_if_mixed();
  }
}

void device_apply_BCs_v1(level_type *level, int x_id, int shape) {
  if (use_cuda_function(8)) {
    cuda_apply_BCs_v1(*level, x_id, shape);
  } else {
    do_sync_if_mixed();
    // CD: This does not seem to be used
    printf("Not implemented!\n");
    exit(-1);
  }
}

void device_apply_BCs_v2(level_type *level, int x_id, int shape) {
  if (use_cuda_function(9)) {
    cuda_apply_BCs_v2(*level, x_id, shape);
  } else {
    do_sync_if_mixed();
    directives_apply_BCs_v2(level, x_id, shape);
    // directives_apply_BCs_v2_unoptimized(level, x_id, shape);
    do_sync_if_mixed();
  }
}

void device_apply_BCs_v4(level_type *level, int x_id, int shape) {
  if (use_cuda_function(10)) {
    cuda_apply_BCs_v4(*level, x_id, shape);
  } else {
    do_sync_if_mixed();
    directives_apply_BCs_v4(level, x_id, shape);
    do_sync_if_mixed();
  }
}

void device_interpolation_v4(level_type *level_f, int id_f, double prescale_f, level_type *level_c, int id_c, communicator_type *interpolation, int block_type) {
  if (use_cuda_function(11)) {
    cuda_interpolation_v4(*level_f, id_f, prescale_f, *level_c, id_c, *interpolation, block_type);
  } else {
    do_sync_if_mixed();
    directives_interpolation_v4(level_f, id_f, prescale_f, level_c, id_c, interpolation, block_type);
    do_sync_if_mixed();
  }
}

void device_interpolation_v2(level_type *level_f, int id_f, double prescale_f, level_type *level_c, int id_c, communicator_type *interpolation, int block_type) {
  if (use_cuda_function(12)) {
    cuda_interpolation_v2(*level_f, id_f, prescale_f, *level_c, id_c, *interpolation, block_type);
  } else {
    do_sync_if_mixed();
    directives_interpolation_v2(level_f, id_f, prescale_f, level_c, id_c, interpolation, block_type);
    do_sync_if_mixed();
  }
}

void device_zero_vector(level_type *d_level, int id)
{
  if (use_cuda_function(13)) {
    cuda_zero_vector(*d_level, id);
  } else {
    do_sync_if_mixed();
    directives_zero_vector(d_level, id);
    do_sync_if_mixed();
  }
}

void device_scale_vector(level_type *d_level, int id_c, double scale_a, int id_a) {
  if (use_cuda_function(14)) {
    cuda_scale_vector(*d_level, id_c, scale_a, id_a);
  } else {
    do_sync_if_mixed();
    directives_scale_vector(d_level, id_c, scale_a, id_a);
    do_sync_if_mixed();
  }
}

void device_shift_vector(level_type *d_level, int id_c, double shift_a, int id_a) {
  if (use_cuda_function(15)) {
    cuda_shift_vector(*d_level, id_c, shift_a, id_a);
  } else {
    do_sync_if_mixed();
    directives_shift_vector(d_level, id_c, shift_a, id_a);
    do_sync_if_mixed();
  }
}

void device_mul_vectors(level_type *d_level, int id_c, double scale, int id_a, int id_b) {
  if (use_cuda_function(16)) {
    cuda_mul_vectors(*d_level, id_c, scale, id_a, id_b);
  } else {
    do_sync_if_mixed();
    directives_mul_vectors(d_level, id_c, scale, id_a, id_b);
    do_sync_if_mixed();
  }
}

void device_add_vectors(level_type *d_level, int id_c, double scale_a, int id_a, double scale_b, int id_b)
{
  if (use_cuda_function(17)) {
    cuda_add_vectors(*d_level, id_c, scale_a, id_a, scale_b, id_b);
  } else {
    do_sync_if_mixed();
    directives_add_vectors(d_level, id_c, scale_a, id_a, scale_b, id_b);
    do_sync_if_mixed();
  }
}

double device_sum(level_type *d_level, int id) {
  double dsum;
  if (use_cuda_function(18)) {
    dsum = cuda_sum(*d_level, id);
  } else {
    do_sync_if_mixed();
    // CD: This function does not seem to be used.
    dsum = directives_sum(d_level, id);
    do_sync_if_mixed();
  }
  return dsum;
}

double device_max_abs(level_type *d_level, int id) {
  double dmax_abs;
  if (use_cuda_function(19)) {
    dmax_abs = cuda_max_abs(*d_level, id);
  } else {
    do_sync_if_mixed();
    dmax_abs = directives_max_abs(d_level, id);
    do_sync_if_mixed();
  }
  return dmax_abs;
}

//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
// Lu = a*alpha[]*u[] - b*divergence beta[]*gradient u[]
//------------------------------------------------------------------------------------------------------------------------------
#ifndef DEFINES_H
#define DEFINES_H
//------------------------------------------------------------------------------------------------------------------------------
#define  VECTOR_TEMP         0 // 
#define  VECTOR_UTRUE        1 // exact solution used to generate f
#define  VECTOR_F_MINUS_AV   2 // cell centered residual (f-Av)
//------------------------------------------------------------------------------------------------------------------------------
#define  VECTOR_F            3 // original right-hand side (Au=f), cell centered
#define  VECTOR_U            4 // numerical solution
#define  VECTOR_ALPHA        5 // cell centered coefficient
#define  VECTOR_BETA_I       6 // face centered coefficient (n.b. element 0 is the left face of the ghost zone element)
#define  VECTOR_BETA_J       7 // face centered coefficient (n.b. element 0 is the back face of the ghost zone element)
#define  VECTOR_BETA_K       8 // face centered coefficient (n.b. element 0 is the bottom face of the ghost zone element)
//------------------------------------------------------------------------------------------------------------------
#define  VECTOR_DINV         9 // cell centered relaxation parameter (e.g. inverse of the diagonal)
#define  VECTOR_L1INV       10 // cell centered relaxation parameter (e.g. inverse of the L1 norm of each row)
#define  VECTOR_VALID       11 // cell centered array noting which cells are actually present
//------------------------------------------------------------------------------------------------------------------
#define VECTORS_RESERVED    12 // total number of vectors and the starting location for any auxillary bottom solver vectors
//------------------------------------------------------------------------------------------------------------------------------

#ifdef SPEC_CUDA
# define USE_REG
# define USE_TEX
#endif

#ifndef HOST_LEVEL_SIZE_THRESHOLD
# define HOST_LEVEL_SIZE_THRESHOLD 10000
#endif

/* Here to map level data: map(level[:LEVEL_ITEMS]) */
#ifdef SPEC_OPENMP_TARGET
# ifdef SPEC_MANAGED_MEMORY
/* Map the pointed-to level data structure */
#  define LEVEL_ITEMS 1
# else
/* The level data structure is already on the device - attach pointers only */
#  define LEVEL_ITEMS 0
# endif
#endif

/* Several possibilities to use managed memory (SPEC_MANAGED_MEMORY):
   1. Compile with the PGI compiler and the option -ta=tesla:managed on NVIDIA GPUs.
   2. Rely on a hardware capability to access system allocated data on the GPUs,
      e.g. IBM AC922 system or an x86 system with a HMM-enabled kernel.
   3. Use the CUDA API on NVIDIA GPUs to allocate data with cudaMallocManaged. */
#if defined(SPEC_CUDA) || defined(SPEC_CUDA_API)
# ifndef SPEC_MANAGED_MEMORY
#  error "Must specify SPEC_MANAGED_MEMORY when using CUDA or the CUDA API"
# endif
/* The CUDA_UM_ALLOC definition is required.
   The CUDA_UM_ZERO_COPY definition is not required - the user can optionally
   compile with -DCUDA_UM_ZERO_COPY to potentially improve performance */
# ifndef CUDA_UM_ALLOC
#  define CUDA_UM_ALLOC
# endif
#endif

#endif
